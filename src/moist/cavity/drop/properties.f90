!> DROP geometric property routines
submodule(moist_cavity_drop) moist_cavity_drop_properties
   !$ use omp_lib, only: omp_get_thread_num
   use mctc_io_constants, only: pi
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_tile, only: surface_tile_width, surface_tile_candidates, &
      & surface_jet_tile
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   use moist_math_sorter_quicksort, only: qsort
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                             Grid adjacency list setup                             *!
   !* ================================================================================= *!

   !> Build grid point neighbour list for density computation
   !>
   !> @param[inout] self  Cavity instance (projected surface rows `xyz` must be filled)
   !> @param[out]   error Error object
   module subroutine setup_grid_adj_list(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: cutoff, rgrid
      type(prettyprinter) :: pp

      cutoff = self%param%adj_list_grid_cutoff

      call self%grid_adj_list%init(cutoff=cutoff)
      call self%grid_adj_list%update(self%xyz(:, 1:self%ngrid), self%ctx%get_num_threads())

      if (self%ctx%verbosity(4)) then
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%blank()
         call pp%push("Grid neighbour list:")
         call pp%kv("Cutoff distance", cutoff, "bohr")
         call pp%kv("Number of points", self%ngrid)
         call pp%kv("Maximum neighbours", maxval(self%grid_adj_list%nnl))
         call pp%kv("Minimum neighbours", minval(self%grid_adj_list%nnl))
         ! Real point count: its square overflows a default integer past 46340 points
         rgrid = real(max(1, self%ngrid), wp)
         call pp%kv("Average neighbours", &
                    real(sum(self%grid_adj_list%nnl), wp)/rgrid, fmt="(f14.4)")
         call pp%kv("Average neighbours (%)", &
                    real(sum(self%grid_adj_list%nnl), wp)/rgrid/rgrid*100.0_wp, "%", &
                    fmt="(f14.4)")
         call pp%kv("Total pairs", size(self%grid_adj_list%nlat))
         call pp%kv("Total pairs (%)", &
                    real(size(self%grid_adj_list%nlat), wp)/rgrid/rgrid*100.0_wp, "%", fmt="(f14.4)")
         call pp%pop()
      end if

   end subroutine setup_grid_adj_list

   !* ================================================================================= *!
   !*                  Grid point density (reference cavity and cavity)                 *!
   !* ================================================================================= *!

   !> Compute grid point densities for soft and hard cavities
   !>
   !> Uses Wendland kernel to compute densities at each grid point by summing
   !> contributions from neighboring points
   !>
   !> Computes both soft cavity (actual cavity) density (rho_grid) and hard cavity
   !> reference density (rho_grid_anchor)
   !>
   !> The soft density uses the prebuilt grid adjacency list (grid_adj_list)
   !>
   !> The hard density is evaluated against the full (unfiltered) Lebedev grid of its owner sphere
   !>
   !> @param[inout] self Cavity instance with grid points and grid_adj_list
   module subroutine compute_grid_point_density(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid, jgrid, jj, isph
      real(wp) :: d, d2, h, mean_wleb, four_h2
      real(wp) :: xi, yi, zi, xj, yj, zj, dx, dy, dz
      real(wp) :: kval
      type(smoothing_kernel_wendland_type) :: kernel
      integer :: nthreads

      nthreads = self%ctx%get_num_threads()
      ! Allocate output arrays
      if (allocated(self%rho_grid)) deallocate (self%rho_grid)
      allocate (self%rho_grid(self%ngrid), source=0.0_wp)
      if (allocated(self%rho_scal0)) deallocate (self%rho_scal0)
      allocate (self%rho_scal0(self%ngrid), source=0.0_wp)
      if (allocated(self%rho_grid_anchor)) deallocate (self%rho_grid_anchor)
      allocate (self%rho_grid_anchor(self%ngrid), source=0.0_wp)

      ! Compute mean weight for given Lebedev level
      mean_wleb = 4.0_wp*pi/real(self%param%num_leb, wp)

      ! Set kernel parameters from cavity parameters
      h = self%param%rho_grid_h
      ! `rho_grid_h` is a user-tunable cavity parameter, so a non-positive value
      ! reaches the kernel as an ordinary bad input rather than a bug
      call kernel%init(order=2, dimension=2, h=h, error=error)
      if (allocated(error)) return
      four_h2 = self%param%adj_list_grid_cutoff*self%param%rho_grid_h

      !> Soft cavity density (adjacency-list-accelerated)

      ! Each pair (i,j) is stored in both directions in the CSR list, so we
      ! simply accumulate all neighbour contributions into each point
      ! Self-contribution (diagonal) is added separately

      !$omp parallel do num_threads(nthreads) default(shared) &
      !$omp& private(igrid, jj, jgrid, xi, yi, zi, dx, dy, dz, d2, d, kval) &
      !$omp& schedule(dynamic)
      do igrid = 1, self%ngrid
         xi = self%xyz(1, igrid)
         yi = self%xyz(2, igrid)
         zi = self%xyz(3, igrid)

         ! Self-contribution (d=0 => kernel at origin)
         kval = kernel%f0(0.0_wp)
         self%rho_grid(igrid) = kval*self%f(igrid)*self%wleb(igrid)

         ! Neighbour contributions from adjacency list
         do jj = self%grid_adj_list%inl(igrid) + 1, &
            self%grid_adj_list%inl(igrid) + self%grid_adj_list%nnl(igrid)
            jgrid = self%grid_adj_list%nlat(jj)

            dx = xi - self%xyz(1, jgrid)
            dy = yi - self%xyz(2, jgrid)
            dz = zi - self%xyz(3, jgrid)
            d2 = dx*dx + dy*dy + dz*dz

            if (d2 <= four_h2) then
               d = sqrt(d2)
               kval = kernel%f0(d)
               self%rho_grid(igrid) = self%rho_grid(igrid) &
                                      + kval*self%f(jgrid)*self%wleb(jgrid)
            end if
         end do
      end do
      !$omp end parallel do

      ! Apply the constant 1/mean_wleb factor once
      self%rho_grid(1:self%ngrid) = self%rho_grid(1:self%ngrid)/mean_wleb

      !> Hard cavity reference density (full Lebedev grid per owner sphere)

      ! Evaluate the anchor density at each active point against the complete
      ! (unfiltered) Lebedev grid of its owner sphere
      !$omp parallel do num_threads(nthreads) default(shared) &
      !$omp& private(igrid, isph, jj, xi, yi, zi, xj, yj, zj, &
      !$omp&   dx, dy, dz, d2, d, kval) &
      !$omp& schedule(dynamic)
      do igrid = 1, self%ngrid
         isph = self%owner(igrid)

         ! Anchor coordinates of active point i
         xi = self%anchor%xyz(1, self%anchor_of_row(igrid))
         yi = self%anchor%xyz(2, self%anchor_of_row(igrid))
         zi = self%anchor%xyz(3, self%anchor_of_row(igrid))

         do jj = 1, self%param%num_leb
            ! Full-grid Lebedev point j on owner sphere
            xj = self%mol%xyz(1, isph) + self%radii(isph)*self%ang_grid(1, jj)
            yj = self%mol%xyz(2, isph) + self%radii(isph)*self%ang_grid(2, jj)
            zj = self%mol%xyz(3, isph) + self%radii(isph)*self%ang_grid(3, jj)

            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            d2 = dx*dx + dy*dy + dz*dz

            if (d2 <= four_h2) then
               d = sqrt(d2)
               kval = kernel%f0(d)
               self%rho_grid_anchor(igrid) = self%rho_grid_anchor(igrid) &
                                             + kval*self%ang_weight(jj)
            end if
         end do
      end do
      !$omp end parallel do

      ! Apply the constant 1/mean_wleb factor once
      self%rho_grid_anchor(1:self%ngrid) = self%rho_grid_anchor(1:self%ngrid)/mean_wleb

      ! Check for very small soft densities and replace with hard sphere density
      block
         logical, allocatable :: mask(:)
         integer :: n_small_density
         character(len=256) :: report

         allocate (mask(self%ngrid))
         mask = self%rho_grid(1:self%ngrid) < 1.0e-10_wp
         n_small_density = count(mask)

         where (mask)
            self%rho_grid = self%rho_grid_anchor
         end where

         if (n_small_density > 0) then
            write (report, "(i0,a)") n_small_density, &
               " grid points have very small soft densities. "// &
               "Using hard sphere density for these points. Results may be unreliable."
            call self%ctx%warning(trim(report))
         end if

      end block

      self%rho_scal0 = self%rho_grid_anchor/self%rho_grid

   end subroutine compute_grid_point_density

   !* ================================================================================= *!
   !*                               Cavity area and volume                              *!
   !* ================================================================================= *!

   !> Compute surface area from projected grid
   !>
   !> Calculates area elements a = w_leb * f * r^2 and accumulates per-sphere
   !> areas (asph) and the total area, also computing the Gaussian widths xi
   !>
   !> @param[inout] self Cavity instance
   module subroutine compute_area_volume(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: i, isph
      real(wp) :: area_i, volume_i
      real(wp), allocatable :: asph_local(:), vsph_local(:)
      integer :: nthreads
      !> Volume origin moving with the molecule
      real(wp) :: origin(3)

      nthreads = self%ctx%get_num_threads()
      origin = self%get_volume_origin()
      if (allocated(self%a)) deallocate (self%a)
      allocate (self%a(self%ngrid), source=0.0_wp)
      if (allocated(self%asph)) deallocate (self%asph)
      allocate (self%asph(self%nsph), source=0.0_wp)
      if (.not. allocated(self%total_area)) allocate (self%total_area)
      if (allocated(self%v)) deallocate (self%v)
      allocate (self%v(self%ngrid), source=0.0_wp)
      if (allocated(self%vsph)) deallocate (self%vsph)
      allocate (self%vsph(self%nsph), source=0.0_wp)
      if (.not. allocated(self%total_volume)) allocate (self%total_volume)

      !$omp parallel num_threads(nthreads) default(shared) private(i, isph, area_i, volume_i, asph_local, vsph_local)
      allocate (asph_local(self%nsph), source=0.0_wp)
      allocate (vsph_local(self%nsph), source=0.0_wp)

      !$omp do schedule(static)
      do i = 1, self%ngrid
         isph = self%owner(i)

         ! Compute area element for projected surface
         area_i = self%wleb(i)*self%radii(isph)**2*self%f(i)
         self%a(i) = area_i

         ! Compute volume element using divergence theorem
         volume_i = 1.0_wp/3.0_wp*area_i &
                    *dot_product(self%xyz(:, i) - origin, self%normal0(:, i))
         self%v(i) = volume_i

         ! Accumulate per sphere in thread-local buffers
         asph_local(isph) = asph_local(isph) + area_i
         vsph_local(isph) = vsph_local(isph) + volume_i

      end do
      !$omp end do

      !$omp critical (compute_area_volume_reduce)
      self%asph = self%asph + asph_local
      self%vsph = self%vsph + vsph_local
      !$omp end critical (compute_area_volume_reduce)

      deallocate (asph_local, vsph_local)
      !$omp end parallel

      ! Total area
      self%total_area = sum(self%a)

      ! Total volume
      self%total_volume = sum(self%v)

   end subroutine compute_area_volume

   !* ================================================================================= *!
   !*                              Atom-partitioned area                                *!
   !* ================================================================================= *!

   !> Partition the surface area among the atoms by the level set's atomic shares
   !>
   !> The owner of a grid point is the sphere of its anchor, while the projection
   !> can carry the point onto surface shaped by other atoms; each area element
   !> `a_i` therefore goes to the atoms by their shares of the level set at the
   !> projected point (see [[moist_cavity_drop_lsf_type:atom_partition]])
   !>
   !> - the shares of one point sum to one, so the partitioned areas sum to the
   !>   total area, which is checked
   !> - a level set without atomic shares skips this with a warning
   !> - printed per atom by [[analyze_cavity]]
   !>
   !> @param[inout] self  Cavity instance after `compute_area_volume`
   !> @param[out]   error Error handling
   module subroutine compute_atom_area(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      !> Relative tolerance of the sum of the partitioned areas, above the
      !> round-off of a threaded sum over many points
      real(wp), parameter :: closure_tol = 1.0e-8_wp

      type(drop_worker_slots_type) :: slots
      !> Loop-wide abort latch and the LSF evaluation failure handed to it
      type(drop_abort_latch_type) :: abort
      type(error_type), allocatable :: lsf_error
      !> Thread-local per-atom areas
      real(wp), allocatable :: asph_local(:)
      integer :: i, thread_slot

      if (.not. self%lsf_model%has_atom_partition()) then
         call self%ctx%warning("DROP atom areas skipped: the level set has no atomic partition")
         return
      end if

      allocate (self%asph_partition(self%nsph), source=0.0_wp)

      ! Only the value-level state of the LSF is read
      call slots%init(self%ctx, self%lsf_model, 0)
      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) &
      !$omp& private(i, thread_slot, asph_local, lsf_error)
      allocate (asph_local(self%nsph), source=0.0_wp)
      thread_slot = 1
      !$ thread_slot = omp_get_thread_num() + 1

      !$omp do schedule(dynamic)
      do i = 1, self%ngrid
         if (abort%requested) cycle
         call slots%lsf(thread_slot)%lsf%prepare(self%xyz(:, i), lsf_error)
         if (allocated(lsf_error)) then
            call abort%latch_error(lsf_error, i)
            cycle
         end if
         call slots%lsf(thread_slot)%lsf%atom_partition(self%a(i), asph_local)
      end do
      !$omp end do

      !$omp critical (compute_atom_area_reduce)
      self%asph_partition = self%asph_partition + asph_local
      !$omp end critical (compute_atom_area_reduce)

      deallocate (asph_local)
      !$omp end parallel

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call fatal_error(error, &
                             "Atom area partition aborted. (unreachable in normal execution)")
         end if
         return
      end if

      ! A point without atomic shares would drop its area silently, so this is an
      ! error and not a warning
      if (abs(sum(self%asph_partition) - self%total_area) > closure_tol*abs(self%total_area)) then
         call fatal_error(error, "DROP atom areas do not sum to the total area")
         return
      end if

   end subroutine compute_atom_area

   !* ================================================================================= *!
   !*                                     Curvature                                     *!
   !* ================================================================================= *!

   !> Compute principal, mean, and Gaussian curvatures at all grid points
   !>
   !> Computes the principal curvatures (k1, k2) at each grid point by projecting
   !> the LSF Hessian onto the tangent plane of the implicit surface and solving
   !> the resulting 2x2 eigenvalue problem in closed form
   !>
   !> Mean and Gaussian curvatures are then derived as
   !>   K_M = (k1 + k2) / 2
   !>   K_G = k1 * k2
   !>
   !> The tangent plane is spanned by an orthonormal frame {t1, t2} constructed
   !> from the unit surface normal n = g/||g||; The 2x2 shape operator is given as
   !>   S_ab = (1/||g||) * t_a^T * H * t_b
   !> whose eigenvalues are the principal curvatures
   !>
   !> @param[inout] self Cavity with projected grid points
   module subroutine compute_curvature(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      type(drop_worker_slots_type) :: slots
      integer :: igrid, thread_slot
      !> Tile index, tile count, tile width, first grid point, size and column
      integer :: itile, ntile, tile, igrid0, np, p
      !> Shared candidate list of the tile
      integer, allocatable :: cand(:)
      integer :: nc
      !> Tile points and level-set jet, packed Hessian (order of [[moist_math_packed]])
      real(wp), allocatable :: tile_xyz(:, :), tile_lsf0(:), tile_lsf1(:, :), tile_lsf2(:, :)
      !> Packed Hessian of the current point
      real(wp) :: lsf2_rr_loc(6)
      real(wp) :: g_vec(3)
      real(wp) :: g_norm, g_norm_sq, inv_g_norm
      real(wp) :: n_vec(3), t1(3), t2(3), Ht1(3), Ht2(3)
      real(wp) :: S11, S12, S22, half_trace, half_diff, disc
      !> Loop-wide abort latch
      type(drop_abort_latch_type) :: abort

      ! Initialize curvature arrays
      if (allocated(self%k1)) deallocate (self%k1)
      allocate (self%k1(self%ngrid), source=0.0_wp)
      if (allocated(self%k2)) deallocate (self%k2)
      allocate (self%k2(self%ngrid), source=0.0_wp)
      if (allocated(self%KM)) deallocate (self%KM)
      allocate (self%KM(self%ngrid), source=0.0_wp)
      if (allocated(self%KG)) deallocate (self%KG)
      allocate (self%KG(self%ngrid), source=0.0_wp)

      ! Set up thread-local LSF evaluators
      ! The curvature loop below reads the LSF Hessian (lsf2_rr), so the cached
      ! callback LSF must be told to compute up to second order
      call slots%init(self%ctx, self%lsf_model, 2)
      tile = surface_tile_width
      ntile = (self%ngrid + tile - 1)/tile

      call abort%reset()

      ! Tiles of `surface_tile_width` points: one batched order-2 jet per tile
      ! over one shared candidate list; num_threads pins the team to the size the
      ! slots were built for, a larger live team would index slots%lsf out of bounds
      !$omp parallel num_threads(slots%nthreads) default(shared) &
      !$omp& private(itile, igrid0, np, p, igrid, thread_slot, cand, nc, &
      !$omp&   tile_xyz, tile_lsf0, tile_lsf1, tile_lsf2, lsf2_rr_loc, &
      !$omp&   g_vec, g_norm, g_norm_sq, inv_g_norm, &
      !$omp&   n_vec, t1, t2, Ht1, Ht2, &
      !$omp&   S11, S12, S22, half_trace, half_diff, disc)
      thread_slot = 1
      !$ thread_slot = omp_get_thread_num() + 1

      allocate (cand(max(1, slots%lsf(thread_slot)%lsf%ncenters)))
      allocate (tile_xyz(3, tile), tile_lsf0(tile), tile_lsf1(3, tile), tile_lsf2(6, tile))

      !$omp do schedule(dynamic, 1)
      do itile = 1, ntile
         if (abort%requested) cycle

         igrid0 = (itile - 1)*tile + 1
         np = min(tile, self%ngrid - igrid0 + 1)
         tile_xyz(:, 1:np) = self%xyz(:, igrid0:igrid0 + np - 1)

         ! Gradient (g) and Hessian (H) of the whole tile; a failure is parked
         ! for the post-region promotion and the tile is dropped
         call surface_tile_candidates(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand, nc)
         if (.not. surface_jet_tile(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand(1:nc), &
                                    igrid0, "compute_curvature", 2, tile_lsf0, tile_lsf1, &
                                    lsf2_rr=tile_lsf2, abort=abort)) cycle

         do p = 1, np
            igrid = igrid0 + p - 1

            g_vec = tile_lsf1(:, p)
            lsf2_rr_loc = tile_lsf2(:, p)

            ! Gradient magnitude
            g_norm_sq = g_vec(1)*g_vec(1) + g_vec(2)*g_vec(2) + g_vec(3)*g_vec(3)
            g_norm = sqrt(g_norm_sq)
            inv_g_norm = 1.0_wp/g_norm

            ! Unit surface normal
            n_vec(1) = g_vec(1)*inv_g_norm
            n_vec(2) = g_vec(2)*inv_g_norm
            n_vec(3) = g_vec(3)*inv_g_norm

            ! Build orthonormal tangent frame {t1, t2}
            ! Choose reference axis least aligned with n to avoid cancellation
            if (abs(n_vec(1)) < 0.9_wp) then
               ! t1 = normalize(e_x x n) = normalize(0, n_z, -n_y)
               t1(1) = 0.0_wp
               t1(2) = n_vec(3)
               t1(3) = -n_vec(2)
            else
               ! t1 = normalize(e_y x n) = normalize(-n_z, 0, n_x)
               t1(1) = -n_vec(3)
               t1(2) = 0.0_wp
               t1(3) = n_vec(1)
            end if
            t1 = t1*(1.0_wp/sqrt(t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3)))

            ! t2 = n x t1 (unit length by construction)
            t2(1) = n_vec(2)*t1(3) - n_vec(3)*t1(2)
            t2(2) = n_vec(3)*t1(1) - n_vec(1)*t1(3)
            t2(3) = n_vec(1)*t1(2) - n_vec(2)*t1(1)

            ! Shape operator matrix elements: S_ab = (1/||g||) * t_a^T H t_b
            Ht1(1) = lsf2_rr_loc(1)*t1(1) + lsf2_rr_loc(2)*t1(2) + lsf2_rr_loc(4)*t1(3)
            Ht1(2) = lsf2_rr_loc(2)*t1(1) + lsf2_rr_loc(3)*t1(2) + lsf2_rr_loc(5)*t1(3)
            Ht1(3) = lsf2_rr_loc(4)*t1(1) + lsf2_rr_loc(5)*t1(2) + lsf2_rr_loc(6)*t1(3)

            Ht2(1) = lsf2_rr_loc(1)*t2(1) + lsf2_rr_loc(2)*t2(2) + lsf2_rr_loc(4)*t2(3)
            Ht2(2) = lsf2_rr_loc(2)*t2(1) + lsf2_rr_loc(3)*t2(2) + lsf2_rr_loc(5)*t2(3)
            Ht2(3) = lsf2_rr_loc(4)*t2(1) + lsf2_rr_loc(5)*t2(2) + lsf2_rr_loc(6)*t2(3)

            S11 = (t1(1)*Ht1(1) + t1(2)*Ht1(2) + t1(3)*Ht1(3))*inv_g_norm
            S12 = (t2(1)*Ht1(1) + t2(2)*Ht1(2) + t2(3)*Ht1(3))*inv_g_norm
            S22 = (t2(1)*Ht2(1) + t2(2)*Ht2(2) + t2(3)*Ht2(3))*inv_g_norm

            ! Closed-form eigenvalues of the 2x2 shape operator
            half_trace = 0.5_wp*(S11 + S22)
            half_diff = 0.5_wp*(S11 - S22)
            disc = sqrt(half_diff*half_diff + S12*S12)

            ! Principal curvatures (k1 >= k2 by construction)
            self%k1(igrid) = half_trace + disc
            self%k2(igrid) = half_trace - disc

            ! Mean and Gaussian curvature from principal curvatures
            self%KM(igrid) = half_trace
            self%KG(igrid) = self%k1(igrid)*self%k2(igrid)
         end do
      end do
      !$omp end do

      deallocate (cand, tile_xyz, tile_lsf0, tile_lsf1, tile_lsf2)
      !$omp end parallel

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call fatal_error(error, &
                             "Curvature computation aborted. (unreachable in normal execution)")
         end if
         return
      end if

   end subroutine compute_curvature

   !* ================================================================================= *!
   !*                              Gaussian surface charges                             *!
   !* ================================================================================= *!

   !> Recompute gaussian surface charge widths with projected wleb, after CP
   !> Jacobian scaling
   !>
   !> - the "true" xi used for CPCM energy and gradients
   module subroutine compute_gaussians(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: i

      if (allocated(self%xi0)) deallocate (self%xi0)
      allocate (self%xi0(self%ngrid), source=0.0_wp)

      do i = 1, self%ngrid
         self%xi0(i) = self%iswig%xi0(self%owner(i), self%wleb(i))
      end do

   end subroutine compute_gaussians

   !* ================================================================================= *!
   !*                                  Diagnostics                                      *!
   !* ================================================================================= *!

   !> Report curvature-based diagnostics of the computed cavity grid
   !>
   !> With curvatures present
   !>  - k1 and k2 spread (min/mean/max)
   !>  - surface area by curvature regime (convex, saddle, pit)
   !>  - effective probe radius from points with k2 < -1/r_A_max
   !>  - watershed patches of k2, with a per-patch table at verbosity 3
   !>
   !> With `do_fine` and `do_atomarea`, also lists the owner-summed (naive) and
   !> partitioned area of every atom
   !>
   !> @param[inout] self   Cavity instance (after full update)
   !> @param[out]   error  Error handling structure
   module subroutine analyze_cavity(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      type(prettyprinter) :: pp
      integer :: i, ngrid

      ! Effective probe sphere analysis
      real(wp) :: mean_radius, r_A_max, r_A_min, k2_bound
      real(wp) :: sum_area, w_sum
      real(wp) :: r_probe_harmonic, r_probe_mean, r_probe_min, r_probe_max, r_probe_stdev
      real(wp) :: k1_lo, k1_hi, k1_mean, k2_lo, k2_hi, k2_mean
      real(wp) :: area_convex, area_saddle, area_pit, area_total
      real(wp) :: pct_convex, pct_saddle, pct_pit
      real(wp), allocatable :: r_probe_vals(:), r_probe_weights(:)
      integer  :: n_probe, n_k2_any_neg
      integer  :: n_convex, n_saddle, n_pit

      ! Watershed-based patch probe-radius analysis
      integer, parameter :: PATCH_RIDGE = -1
      integer, allocatable :: patch_id(:), eligible_idx(:)
      real(wp), allocatable :: patch_area(:), patch_k2_min(:), patch_r(:), eligible_k2(:)
      real(wp), allocatable :: patch_k2_wsum(:), patch_k2_mean(:)
      integer :: npatch, n_eligible, n_ridge, k_idx, cur_i, jj_p, jneigh, ip, ip_largest
      integer :: first_p, j_pid
      logical :: multi_patches
      real(wp) :: rp_amean, rp_harmonic, rp_min, rp_max, patch_area_total
      type(error_type), allocatable :: sort_error
      type(prettylistprinter) :: plp
      !> Neighbour list of the watershed, local to the diagnostic
      type(adjacency_list_type) :: adj
      integer :: iat

      ngrid = self%ngrid

      !> Effective probe sphere radius routine
      if (allocated(self%k1) .and. allocated(self%k2) .and. &
          allocated(self%a) .and. self%nsph > 0) then

         mean_radius = sum(self%radii)/real(self%nsph, wp)
         r_A_max = maxval(self%radii)
         r_A_min = minval(self%radii)

         ! Global probe-radius upper bound
         k2_bound = -1.0_wp/r_A_max

         ! Distribution summary (min/mean/max) for k1 and k2; quick reality check on the sign convention
         k1_lo = minval(self%k1(1:ngrid))
         k1_hi = maxval(self%k1(1:ngrid))
         k1_mean = sum(self%k1(1:ngrid))/real(ngrid, wp)
         k2_lo = minval(self%k2(1:ngrid))
         k2_hi = maxval(self%k2(1:ngrid))
         k2_mean = sum(self%k2(1:ngrid))/real(ngrid, wp)

         n_probe = 0
         n_k2_any_neg = 0
         n_convex = 0; n_saddle = 0; n_pit = 0
         area_convex = 0.0_wp; area_saddle = 0.0_wp; area_pit = 0.0_wp

         allocate (r_probe_vals(ngrid), r_probe_weights(ngrid))

         do i = 1, ngrid

            ! Four-quadrant classification of (k1, k2) signs:
            !  - convex  (k1>0, k2>0): atom-like (convex in both directions)
            !  - saddle  (k1>0, k2<0): probe rolling in a 2-atom crease
            !  - pit     (k1<0, k2<0): spherical pocket at a 3+ atom junction
            if (self%k1(i) > 0.0_wp .and. self%k2(i) > 0.0_wp) then
               n_convex = n_convex + 1
               area_convex = area_convex + self%a(i)
            else if (self%k1(i) > 0.0_wp .and. self%k2(i) < 0.0_wp) then
               n_saddle = n_saddle + 1
               area_saddle = area_saddle + self%a(i)
            else if (self%k1(i) < 0.0_wp .and. self%k2(i) < 0.0_wp) then
               n_pit = n_pit + 1
               area_pit = area_pit + self%a(i)
            else
               call fatal_error(error, "Invalid curvature regime")
               return
            end if

            if (self%k2(i) < 0.0_wp) n_k2_any_neg = n_k2_any_neg + 1

            !> Global probe filter: concavity sharper than 1/r_A_max
            ! (= probe radius smaller than the largest atom in the system)
            if (self%k2(i) < k2_bound) then
               n_probe = n_probe + 1
               r_probe_vals(n_probe) = -1.0_wp/self%k2(i)
               r_probe_weights(n_probe) = self%a(i)
            end if
         end do

         !> Area-weighted probe-radius statistics
         if (n_probe > 0) then
            w_sum = sum(r_probe_weights(1:n_probe))
            sum_area = w_sum
         else
            w_sum = 0.0_wp
            sum_area = 0.0_wp
         end if

         if (n_probe > 0 .and. w_sum > 0.0_wp) then
            r_probe_mean = sum(r_probe_weights(1:n_probe) &
                               *r_probe_vals(1:n_probe))/w_sum
            r_probe_harmonic = w_sum &
                               /sum(r_probe_weights(1:n_probe) &
                                    /r_probe_vals(1:n_probe))
            r_probe_min = minval(r_probe_vals(1:n_probe))
            r_probe_max = maxval(r_probe_vals(1:n_probe))
            r_probe_stdev = sqrt(sum(r_probe_weights(1:n_probe) &
                                     *(r_probe_vals(1:n_probe) - r_probe_mean)**2)/w_sum)
         end if

         area_total = area_convex + area_saddle + area_pit
         if (area_total > 0.0_wp) then
            pct_convex = 100.0_wp*area_convex/area_total
            pct_saddle = 100.0_wp*area_saddle/area_total
            pct_pit = 100.0_wp*area_pit/area_total
         else
            pct_convex = 0.0_wp; pct_saddle = 0.0_wp; pct_pit = 0.0_wp
         end if

         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%blank()
         call pp%push("Effective probe sphere diagnostics")
         call pp%push("Curvature spread")
         call pp%kvvv("k1 [min/mean/max]", k1_lo, k1_mean, k1_hi, unit="Bohr^-1")
         call pp%kvvv("k2 [min/mean/max]", k2_lo, k2_mean, k2_hi, unit="Bohr^-1")
         call pp%pop()

         call pp%push("Surface area by curvature regime")
         call pp%kv2("Convex  (k1>0, k2>0)", area_convex, "Bohr^2", pct_convex, "%")
         call pp%kv2("Saddle  (k1>0, k2<0)", area_saddle, "Bohr^2", pct_saddle, "%")
         call pp%kv2("Pit     (k1<0, k2<0)", area_pit, "Bohr^2", pct_pit, "%")
         call pp%pop()

         call pp%push("Probe-region filter  k2 < -1/r_A_max")
         call pp%kv2("r_A (min, max)", r_A_min, "Bohr", r_A_max, "Bohr")
         if (n_k2_any_neg == 0) then
            call pp%kv("Status", &
                       "cavity is fully convex -- no probe-rolling regions detected")
         else
            if (n_probe > 0 .and. w_sum > 0.0_wp) then
               call pp%kv2("Probe-region", sum_area, "Bohr^2", &
                           100.0_wp*sum_area/area_total, "%")
               call pp%kv("r_probe (min)", &
                          r_probe_min, unit="Bohr")
               call pp%kv("r_probe (mean)", &
                          r_probe_mean, unit="Bohr")
               call pp%kv("r_probe (-1/<k2>)", &
                          r_probe_harmonic, unit="Bohr")
               call pp%kv("r_probe (stdev)", &
                          r_probe_stdev, unit="Bohr")
            end if
         end if
         call pp%pop()

         !> Watershed segmentation of the non-convex region into basins
         n_eligible = count(self%k2(1:ngrid) < 0.0_wp)
         if (n_eligible > 0) then

            ! The diagnostic owns its neighbour list; the cavity keeps none unless
            ! the grid density asked for one
            call adj%init(cutoff=self%param%adj_list_grid_cutoff)
            call adj%update(self%xyz(:, 1:ngrid), self%ctx%get_num_threads())

            ! Collect eligible (concave) points and their k2 values
            allocate (eligible_idx(n_eligible))
            allocate (eligible_k2(n_eligible))
            ip = 0
            do cur_i = 1, ngrid
               if (self%k2(cur_i) < 0.0_wp) then
                  ip = ip + 1
                  eligible_idx(ip) = cur_i
                  eligible_k2(ip) = self%k2(cur_i)
               end if
            end do

            ! Sort eligible points by ascending k2 (deepest concavity first)
            if (n_eligible > 1) call qsort(eligible_k2, eligible_idx, sort_error)

            allocate (patch_id(ngrid), source=0)
            allocate (patch_area(ngrid), source=0.0_wp)
            allocate (patch_k2_min(ngrid), source=0.0_wp)
            allocate (patch_r(ngrid), source=0.0_wp)
            allocate (patch_k2_wsum(ngrid), source=0.0_wp)
            allocate (patch_k2_mean(ngrid), source=0.0_wp)

            npatch = 0
            n_ridge = 0
            do k_idx = 1, n_eligible
               cur_i = eligible_idx(k_idx)

               ! Inspect neighbour labels: collect first labelled patch and detect whether any
               ! second distinct patch shows up
               first_p = 0
               multi_patches = .false.
               do jj_p = adj%inl(cur_i) + 1, adj%inl(cur_i) + adj%nnl(cur_i)
                  jneigh = adj%nlat(jj_p)
                  j_pid = patch_id(jneigh)
                  if (j_pid > 0) then
                     if (first_p == 0) then
                        first_p = j_pid
                     else if (j_pid /= first_p) then
                        multi_patches = .true.
                        exit
                     end if
                  end if
               end do

               if (first_p == 0) then
                  ! Local k2 minimum within the eligible set -> new basin
                  npatch = npatch + 1
                  patch_id(cur_i) = npatch
                  patch_k2_min(npatch) = self%k2(cur_i)
                  patch_area(npatch) = self%a(cur_i)
                  patch_k2_wsum(npatch) = self%k2(cur_i)*self%a(cur_i)
               else if (.not. multi_patches) then
                  ! Single neighbouring basin -> extend it
                  patch_id(cur_i) = first_p
                  patch_area(first_p) = patch_area(first_p) + self%a(cur_i)
                  patch_k2_wsum(first_p) = patch_k2_wsum(first_p) &
                                           + self%k2(cur_i)*self%a(cur_i)
               else
                  ! Multiple neighbouring basins meet here -> ridge point
                  patch_id(cur_i) = PATCH_RIDGE
                  n_ridge = n_ridge + 1
               end if
            end do

            ! Per-patch probe radius set by the deepest point in the basin, and an
            ! area-weighted mean k2 across each basin; Lebedev weights vary over
            ! the surface, so a plain arithmetic average would over-emphasise
            ! low-weight points near sphere intersections
            do ip = 1, npatch
               patch_r(ip) = -1.0_wp/patch_k2_min(ip)
               if (patch_area(ip) > 0.0_wp) then
                  patch_k2_mean(ip) = patch_k2_wsum(ip)/patch_area(ip)
               end if
            end do

            patch_area_total = sum(patch_area(1:npatch))

            call pp%push("Watershed patch analysis (basins of k2)")
            call pp%kv("Ridge points (excluded)", n_ridge)
            call pp%kv("Patch count", npatch)
            if (npatch > 0 .and. patch_area_total > 0.0_wp) then
               ip_largest = maxloc(patch_area(1:npatch), dim=1)
               rp_min = minval(patch_r(1:npatch))
               rp_max = maxval(patch_r(1:npatch))
               rp_amean = sum(patch_r(1:npatch)*patch_area(1:npatch)) &
                          /patch_area_total
               ! Area-weighted harmonic mean: average k2 in curvature space,
               ! then invert
               rp_harmonic = -1.0_wp/(sum(patch_k2_min(1:npatch) &
                                          *patch_area(1:npatch))/patch_area_total)

               if (area_total > 0.0_wp) then
                  call pp%kv2("Patch area total", patch_area_total, "Bohr^2", &
                              100.0_wp*patch_area_total/area_total, "%")
               else
                  call pp%kv("Patch area total", patch_area_total, unit="Bohr^2")
               end if
               ! r_probe (mean): area-weighted harmonic mean of the per-patch
               ! min-curvature probe radius r_p = -1/k2_min(p)
               call pp%kv("r_probe (mean)", rp_harmonic, unit="Bohr")
               call pp%kv("r_probe (sharpest)", rp_min, unit="Bohr")
               call pp%kv("r_probe (broadest)", rp_max, unit="Bohr")
               call pp%kv("Largest patch", patch_area(ip_largest), unit="Bohr^2")
               call pp%kv("Effective probe radius", rp_amean, unit="Bohr")
            end if
            call pp%pop()

            if (self%ctx%verbosity(3)) then
               if (npatch > 0) then
                  call pp%blank()
                  call pp%push("Per-patch curvatures")
                  plp = new_prettylistprinter( &
                        [8, 16, 18, 18, 16], &
                        [character(len=18) :: "Patch", "Area/Bohr^2", &
                         "k2_mean/Bohr^-1", "k2_min/Bohr^-1", "r_probe/Bohr"], &
                        unit=self%ctx%unit, offset=2, column_gap=2)
                  call plp%print_header()
                  call plp%separator()
                  do ip = 1, npatch
                     call plp%begin_row()
                     call plp%add(ip)
                     call plp%add(patch_area(ip))
                     call plp%add(patch_k2_mean(ip))
                     call plp%add(patch_k2_min(ip))
                     call plp%add(patch_r(ip))
                     call plp%end_row()
                  end do
                  call pp%pop()
               end if
            end if
         end if

      end if

      !> Owner-summed against level-set-partitioned area per atom
      if (self%param%do_atomarea .and. allocated(self%asph_partition)) then
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%blank()
         call pp%push("Atom-resolved surface areas")
         plp = new_prettylistprinter( &
               [6, 6, 16, 16], &
               [character(len=16) :: "Atom", "Sym", "Naive/Bohr^2", "Partition/Bohr^2"], &
               unit=self%ctx%unit, offset=2, column_gap=2)
         call plp%print_header()
         call plp%separator()
         do iat = 1, self%nsph
            call plp%begin_row()
            call plp%add(iat)
            call plp%add(trim(self%mol%sym(self%mol%id(iat))))
            call plp%add(self%asph(iat), fmt="f16.6")
            call plp%add(self%asph_partition(iat), fmt="f16.6")
            call plp%end_row()
         end do
         call pp%pop()
      end if

   end subroutine analyze_cavity

end submodule moist_cavity_drop_properties
