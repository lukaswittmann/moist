!> Solute-solvent potential tables of the 3D models
!>
!> - Ng-split arrays per solvent atom v on the volume grid, from `potential_sites_type`
!> - u_sr(:, v) = sum_a u_pair(|r - r_a|) + q_v (phi(r) - phi_lr(r))
!> - Short-range shape: (ngrid, nv)
!> - ur_lr(:, v) = q_v phi_lr(r)
!> - phi_lr = sum_a q_a erf(alpha |r - r_a|)/|r - r_a|
!> - uk_lr(:, v) = q_v 4 pi exp(-k^2/(4 alpha^2))/k^2 sum_a q_a exp(-i k.(r_a - r0))
!> - Analytic reciprocal shape: (npts_k, nv)
!> - u_pair: 12-6 potential with interaction mixing, once per atom pair
!> - q_a, q_v: solute and solvent tail charges
!> - Full solute field phi: supplied grid field, otherwise point-charge field of tail charges
!> - Phase reference r0: first grid point; zero k = 0 mode
!> - No split or no solute tail charges: zero tails, full field in u_sr
!> - Category rules from 1D tables; caller validation of alpha and mixing rule
!> - Pair potential or point-charge field: grid points at least 1e-6 bohr from solute atoms
!> - L = sum w_sr u_sr + sum w_lr ur_lr + Re sum conj(w_k) uk_lr
!> - Complex adjoint pairing: Re(conj(w_k) d uk_lr)
!> - Field adjoint: w_phi = sum_v q_v w_sr(:, v)
!> - w_q: solute tail-charge adjoint
!> - Geometry derivatives: dL/dr_a to solute gradient, dL/dr_g to grid adjoint for site-built pieces
!> - Grid point 1 contribution through r0; fixed k-points
!> - Further adjoint propagation: `moz_potential_type%compute_adjoint`
module moist_model_moz_potential_kernel_table_3d
   use mctc_env, only: wp, error_type, fatal_error
   !$ use omp_lib, only: omp_get_thread_num
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_lj_base, only: lj_12_6_mix, lj_12_6_evaluate
   use moist_model_moz_potential_kernel_category, only: resolve_categories
   use moist_model_moz_potential_kernel_ng, only: ng_real_tail, ng_fourier_envelope
   implicit none(type, external)
   private

   public :: evaluate_table_3d, evaluate_table_3d_adjoint

   !> Grid points per block of the real-space loops
   integer, parameter :: block_size = 256
   !> Smallest grid point to atom distance at which singular pieces are evaluated, bohr
   real(wp), parameter :: min_distance = 1.0e-6_wp

   !> Resolved inputs of one evaluation
   type :: table_3d_inputs_type
      !> Whether the pair category enters
      logical :: do_pair = .false.
      !> Whether the electrostatic category enters
      logical :: do_elec = .false.
      !> Whether the Ng tails are on
      logical :: tails = .false.
      !> Whether the solute tail charges build a point-charge field (no grid field)
      logical :: has_point = .false.
      !> Whether a singular piece (pair or point field) is evaluated on the grid
      logical :: check_separation = .false.
      !> Number of solute atoms
      integer :: nu = 0
      !> Number of solvent atoms
      integer :: nv = 0
      !> Charges of the point-charge field, e (nu); zero-sized without one
      real(wp), allocatable :: q_point(:)
      !> Solute tail charges, e (nu); zero-sized without tails
      real(wp), allocatable :: q_tail(:)
      !> Solvent charges, e (nv)
      real(wp), allocatable :: q_v(:)
      !> Mixed Lennard-Jones sigma per atom pair, bohr (nu, nv)
      real(wp), allocatable :: sigma_uv(:, :)
      !> Mixed Lennard-Jones epsilon per atom pair, Hartree (nu, nv)
      real(wp), allocatable :: eps_uv(:, :)
   end type table_3d_inputs_type

   !> Scratch of one thread, sized once per evaluation to a full block
   type :: block_scratch_type
      !> Distances to every solute atom (block_size, nu)
      real(wp), allocatable :: dist(:, :)
      !> Separations to one atom (3, block_size)
      real(wp), allocatable :: diff(:, :)
      !> Pair potential of one atom pair and its radial derivative (block_size)
      real(wp), allocatable :: u(:), du(:)
      !> Accumulated pair sum, full field and tail field (block_size)
      real(wp), allocatable :: acc(:), phi(:), phi_lr(:)
      !> Ng tail and its radial factor (block_size)
      real(wp), allocatable :: f(:), g(:)
      !> Radial adjoint factor and adjoint of phi_lr (block_size)
      real(wp), allocatable :: s(:), w_tail(:)
   end type block_scratch_type

contains

   !> Assemble the 3D Ng-split tables
   !>
   !> @param[in]  grid      volume grid with reciprocal points
   !> @param[in]  coord_u   solute coordinates, bohr (3, nu)
   !> @param[in]  site_u    resolved solute sites
   !> @param[in]  site_v    resolved solvent sites
   !> @param[in]  mixing    lennard-jones mixing rule of the interaction
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[in]  alpha     ng split parameter, 1/bohr; validated by the caller
   !> @param[out] u_sr      short-range potential, Hartree (ngrid, nv)
   !> @param[out] ur_lr     long-range real-space potential, Hartree (ngrid, nv)
   !> @param[out] uk_lr     long-range reciprocal potential, Hartree bohr^3 (npts_k, nv)
   !> @param[out] error     category or input failure
   !> @param[in]  phi_u     optional solute grid field, Hartree/e (ngrid); required when the solute sites have one
   subroutine evaluate_table_3d(grid, coord_u, site_u, site_v, mixing, ng_split, alpha, u_sr, ur_lr, uk_lr, &
         & error, phi_u)
      !> Volume grid with reciprocal points
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved solute sites
      type(potential_sites_type), intent(in) :: site_u
      !> Resolved solvent sites
      type(potential_sites_type), intent(in) :: site_v
      !> Lennard-Jones mixing rule of the interaction
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Short-range potential (ngrid, nv)
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      !> Long-range real-space potential (ngrid, nv)
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      !> Long-range reciprocal potential (npts_k, nv)
      complex(wp), allocatable, intent(out) :: uk_lr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional solute grid field (ngrid)
      real(wp), intent(in), optional :: phi_u(:)
      !> Resolved inputs
      type(table_3d_inputs_type) :: inp
      !> Solute grid field, Hartree/e (ngrid); unallocated when the solute has none
      real(wp), allocatable :: phi_field(:)
      !> Scratch per thread
      type(block_scratch_type), allocatable :: scratch(:)
      !> Lowest solute atom too close to a grid point, per thread
      integer, allocatable :: bad_t(:)
      integer :: ib, nblock, g0, g1, nt, it

      call prepare_inputs(grid, coord_u, site_u, site_v, mixing, ng_split, inp, error)
      if (allocated(error)) return
      if (inp%do_elec .and. site_u%has_field) then
         if (.not. present(phi_u)) then
            call fatal_error(error, "3D tables need the solute grid field of a side with a field")
            return
         end if
         if (size(phi_u) /= grid%ngrid) then
            call fatal_error(error, "3D solute grid field does not match the grid")
            return
         end if
         phi_field = phi_u
      end if

      allocate (u_sr(grid%ngrid, inp%nv), ur_lr(grid%ngrid, inp%nv))
      allocate (uk_lr(grid%npts_k, inp%nv), source=(0.0_wp, 0.0_wp))
      nblock = (grid%ngrid + block_size - 1)/block_size
      nt = grid%team_size()
      call new_scratch(scratch, nt, inp%nu)
      allocate (bad_t(nt), source=huge(1))
      !$omp parallel do default(shared) schedule(static) num_threads(nt) private(ib, g0, g1, it)
      do ib = 1, nblock
         g0 = (ib - 1)*block_size + 1
         g1 = min(ib*block_size, grid%ngrid)
         it = 1
         !$ it = omp_get_thread_num() + 1
         call forward_block(grid%xyz(:, g0:g1), coord_u, inp, alpha, phi_field, g0, g1, &
            & u_sr(g0:g1, :), ur_lr(g0:g1, :), scratch(it), bad_t(it))
      end do
      !$omp end parallel do
      call refuse_coincident(minval(bad_t), error)
      if (allocated(error)) then
         deallocate (u_sr, ur_lr, uk_lr)
         return
      end if
      if (inp%tails) call reciprocal_tail(grid, coord_u, inp, alpha, uk_lr)
   end subroutine evaluate_table_3d

   !> Form adjoints of `evaluate_table_3d`
   !>
   !> - Solute field adjoint: w_phi = sum_v q_v w_sr(:, v)
   !> - w_q: tail-charge adjoint, including point-field contribution without a supplied grid field
   !> - No electrostatics: w_phi and w_q unallocated
   !> - Gradient phase with `gradient` and `grid_adjoint`: nuclear and grid-point derivatives of site-built pieces
   !> - L = sum w_sr u_sr + sum w_lr ur_lr + Re sum conj(w_k) uk_lr
   !>
   !> @param[in]     grid          volume grid with reciprocal points
   !> @param[in]     coord_u       solute coordinates, bohr (3, nu)
   !> @param[in]     site_u        resolved solute sites
   !> @param[in]     site_v        resolved solvent sites
   !> @param[in]     mixing        lennard-jones mixing rule of the interaction
   !> @param[in]     ng_split      whether the Coulomb tail is split off
   !> @param[in]     alpha         ng split parameter, 1/bohr; validated by the caller
   !> @param[in]     w_sr          adjoint of u_sr (ngrid, nv)
   !> @param[in]     w_lr          adjoint of ur_lr (ngrid, nv)
   !> @param[in]     w_k           adjoint of uk_lr (npts_k, nv)
   !> @param[out]    w_phi         adjoint of the solute field (ngrid); unallocated without electrostatics
   !> @param[out]    w_q           tail-charge adjoint (nu), zero-sized without charges; unallocated without electrostatics
   !> @param[out]    error         category or input failure
   !> @param[in,out] gradient      optional solute nuclear gradient to accumulate, Hartree/bohr (3, nu)
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator of the model; passed with `gradient`
   subroutine evaluate_table_3d_adjoint(grid, coord_u, site_u, site_v, mixing, ng_split, alpha, &
         & w_sr, w_lr, w_k, w_phi, w_q, error, gradient, grid_adjoint)
      !> Volume grid with reciprocal points
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved solute sites
      type(potential_sites_type), intent(in) :: site_u
      !> Resolved solvent sites
      type(potential_sites_type), intent(in) :: site_v
      !> Lennard-Jones mixing rule of the interaction
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Adjoint of u_sr (ngrid, nv)
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr (ngrid, nv)
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of uk_lr (npts_k, nv)
      complex(wp), intent(in) :: w_k(:, :)
      !> Adjoint of the solute field (ngrid)
      real(wp), allocatable, intent(out) :: w_phi(:)
      !> Adjoint of the solute tail charges (nu)
      real(wp), allocatable, intent(out) :: w_q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional solute nuclear gradient to accumulate (3, nu)
      real(wp), intent(inout), optional :: gradient(:, :)
      !> Optional volume adjoint accumulator of the model
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      !> Resolved inputs
      type(table_3d_inputs_type) :: inp
      !> Adjoint of the full field (ngrid)
      real(wp), allocatable :: w_field(:)
      !> Grid-point adjoint of the pieces built from the sites (3, ngrid)
      real(wp), allocatable :: w_xyz(:, :)
      !> Per-thread solute gradient (3, nu, nt) and charge adjoints (nu, nt)
      real(wp), allocatable :: grad_t(:, :, :), wq_point_t(:, :), wq_tail_t(:, :)
      !> Solute gradient of the reciprocal tail (3, nu)
      real(wp), allocatable :: grad_k(:, :)
      !> Summed solute gradient and charge adjoints
      real(wp), allocatable :: grad(:, :), wq_point(:), wq_tail(:)
      !> Scratch per thread
      type(block_scratch_type), allocatable :: scratch(:)
      !> Lowest solute atom too close to a grid point, per thread
      integer, allocatable :: bad_t(:)
      logical :: geom
      integer :: ib, nblock, g0, g1, nt, it

      if (present(gradient) .neqv. present(grid_adjoint)) then
         call fatal_error(error, "3D table adjoint takes gradient and grid_adjoint together or neither")
         return
      end if
      geom = present(gradient)
      call prepare_inputs(grid, coord_u, site_u, site_v, mixing, ng_split, inp, error)
      if (allocated(error)) return
      if (any(shape(w_sr) /= [grid%ngrid, inp%nv]) .or. any(shape(w_lr) /= [grid%ngrid, inp%nv])) then
         call fatal_error(error, "3D table adjoints w_sr and w_lr must be (ngrid, nv)")
         return
      end if
      if (any(shape(w_k) /= [grid%npts_k, inp%nv])) then
         call fatal_error(error, "3D table adjoint w_k must be (npts_k, nv)")
         return
      end if
      if (geom) then
         if (any(shape(gradient) /= [3, inp%nu])) then
            call fatal_error(error, "3D table gradient must be (3, nu)")
            return
         end if
         if (.not. grid_adjoint%is_initialized()) then
            call fatal_error(error, "3D table adjoint needs an initialized grid adjoint accumulator")
            return
         end if
         if (grid_adjoint%size() /= grid%ngrid) then
            call fatal_error(error, "3D table grid adjoint accumulator does not match the grid")
            return
         end if
      end if

      nt = grid%team_size()
      nblock = (grid%ngrid + block_size - 1)/block_size
      allocate (w_field(grid%ngrid), w_xyz(3, grid%ngrid), source=0.0_wp)
      allocate (grad_t(3, inp%nu, nt), wq_point_t(inp%nu, nt), wq_tail_t(inp%nu, nt), source=0.0_wp)
      call new_scratch(scratch, nt, inp%nu)
      allocate (bad_t(nt), source=huge(1))
      !$omp parallel do default(shared) schedule(static) num_threads(nt) private(ib, g0, g1, it)
      do ib = 1, nblock
         g0 = (ib - 1)*block_size + 1
         g1 = min(ib*block_size, grid%ngrid)
         it = 1
         !$ it = omp_get_thread_num() + 1
         call adjoint_block(grid%xyz(:, g0:g1), coord_u, inp, alpha, geom, w_sr(g0:g1, :), w_lr(g0:g1, :), &
            & w_field(g0:g1), w_xyz(:, g0:g1), grad_t(:, :, it), wq_point_t(:, it), wq_tail_t(:, it), &
            & scratch(it), bad_t(it))
      end do
      !$omp end parallel do
      call refuse_coincident(minval(bad_t), error)
      if (allocated(error)) return
      grad = sum(grad_t, dim=3)
      wq_point = sum(wq_point_t, dim=2)
      wq_tail = sum(wq_tail_t, dim=2)

      if (inp%tails) then
         grad_t = 0.0_wp
         wq_tail_t = 0.0_wp
         call reciprocal_tail_adjoint(grid, coord_u, inp, alpha, geom, w_k, nt, grad_t, wq_tail_t)
         wq_tail = wq_tail + sum(wq_tail_t, dim=2)
         if (geom) then
            grad_k = sum(grad_t, dim=3)
            grad = grad + grad_k
            ! Reciprocal phase reference r0 at grid point 1
            w_xyz(:, 1) = w_xyz(:, 1) - sum(grad_k, dim=2)
         end if
      end if

      if (geom) then
         call grid_adjoint%add_weights(error, w_xyz=w_xyz)
         if (allocated(error)) return
         gradient = gradient + grad
      end if

      if (.not. inp%do_elec) return
      call move_alloc(w_field, w_phi)
      if (allocated(site_u%q)) then
         w_q = wq_tail
         if (inp%has_point) w_q = w_q + wq_point
      else
         allocate (w_q(0))
      end if
   end subroutine evaluate_table_3d_adjoint

   !> Resolve categories, charges and mixed parameters, check the inputs
   !>
   !> - Supplied solute grid field: tail charges for Ng tails only
   !> - No supplied solute field: tail charges also for point-charge field
   !> - No solvent-side field
   !>
   !> @param[in]  grid      volume grid
   !> @param[in]  coord_u   solute coordinates, bohr (3, nu)
   !> @param[in]  site_u    resolved solute sites
   !> @param[in]  site_v    resolved solvent sites
   !> @param[in]  mixing    lennard-jones mixing rule of the interaction
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[out] inp       resolved inputs
   !> @param[out] error     category or input failure
   subroutine prepare_inputs(grid, coord_u, site_u, site_v, mixing, ng_split, inp, error)
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved solute sites
      type(potential_sites_type), intent(in) :: site_u
      !> Resolved solvent sites
      type(potential_sites_type), intent(in) :: site_v
      !> Lennard-Jones mixing rule of the interaction
      integer, intent(in) :: mixing
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Resolved inputs
      type(table_3d_inputs_type), intent(out) :: inp
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: a, v

      inp%nu = site_u%natom
      inp%nv = site_v%natom
      if (inp%nu < 1 .or. inp%nv < 1) then
         call fatal_error(error, "3D tables need resolved sites on both sides")
         return
      end if
      if (site_v%has_field) then
         call fatal_error(error, "3D solvent sites must not carry a field")
         return
      end if
      if (any(shape(coord_u) /= [3, inp%nu])) then
         call fatal_error(error, "3D table solute coordinates must be (3, natom) of the solute sites")
         return
      end if
      if (grid%ngrid < 1 .or. .not. allocated(grid%xyz)) then
         call fatal_error(error, "3D tables need a grid with geometry")
         return
      end if
      call resolve_categories(allocated(site_u%pair), allocated(site_v%pair), &
         & site_u%has_field .or. allocated(site_u%q), allocated(site_v%q), inp%do_pair, inp%do_elec, error)
      if (allocated(error)) return

      if (inp%do_elec) then
         inp%tails = ng_split .and. allocated(site_u%q)
         if (inp%tails .and. grid%npts_k < 1) then
            call fatal_error(error, "3D Ng tails need a grid with reciprocal points")
            return
         end if
         inp%has_point = allocated(site_u%q) .and. .not. site_u%has_field
         if (inp%has_point) inp%q_point = site_u%q
         if (inp%tails) inp%q_tail = site_u%q
         inp%q_v = site_v%q
      end if
      if (.not. allocated(inp%q_point)) allocate (inp%q_point(0))
      if (.not. allocated(inp%q_tail)) allocate (inp%q_tail(0))
      if (.not. allocated(inp%q_v)) allocate (inp%q_v(0))
      inp%check_separation = inp%do_pair .or. inp%has_point

      if (inp%do_pair) then
         allocate (inp%sigma_uv(inp%nu, inp%nv), inp%eps_uv(inp%nu, inp%nv))
         do v = 1, inp%nv
            do a = 1, inp%nu
               call lj_12_6_mix(mixing, site_u%pair%sigma(a), site_u%pair%epsilon(a), &
                  & site_v%pair%sigma(v), site_v%pair%epsilon(v), inp%sigma_uv(a, v), inp%eps_uv(a, v))
            end do
         end do
      end if
   end subroutine prepare_inputs

   !> Allocate the scratch of every thread once, sized to a full block
   !>
   !> @param[out] scratch  scratch per thread (nt)
   !> @param[in]  nt       team size
   !> @param[in]  nu       number of solute atoms
   subroutine new_scratch(scratch, nt, nu)
      !> Scratch per thread
      type(block_scratch_type), allocatable, intent(out) :: scratch(:)
      !> Team size
      integer, intent(in) :: nt
      !> Number of solute atoms
      integer, intent(in) :: nu
      integer :: it
      allocate (scratch(nt))
      do it = 1, nt
         allocate (scratch(it)%dist(block_size, nu), scratch(it)%diff(3, block_size))
         allocate (scratch(it)%u(block_size), scratch(it)%du(block_size), scratch(it)%acc(block_size), &
            & scratch(it)%phi(block_size), scratch(it)%phi_lr(block_size), scratch(it)%f(block_size), &
            & scratch(it)%g(block_size), scratch(it)%s(block_size), scratch(it)%w_tail(block_size))
      end do
   end subroutine new_scratch

   !> Refuse a grid point within `min_distance` of a solute atom
   !>
   !> @param[in]  bad    lowest offending solute atom; huge when none
   !> @param[out] error  coinciding grid point and atom
   subroutine refuse_coincident(bad, error)
      !> Lowest offending solute atom
      integer, intent(in) :: bad
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=64) :: label
      if (bad == huge(bad)) return
      write (label, "(a, i0)") "solute atom ", bad
      call fatal_error(error, "A grid point lies within 1e-6 bohr of "//trim(label)// &
         & "; the pair and Coulomb potentials are singular there")
   end subroutine refuse_coincident

   !> Assemble real-space tables for one block of grid points
   !>
   !> - Singular piece within `min_distance`: lowest offending atom in `bad`, block skipped
   !> - Refusal by caller
   !>
   !> @param[in]     xyz        grid points of the block, bohr (3, n)
   !> @param[in]     coord_u    solute coordinates, bohr (3, nu)
   !> @param[in]     inp        resolved inputs
   !> @param[in]     alpha      ng split parameter, 1/bohr
   !> @param[in]     phi_field  summed field of the field-supplying terms (ngrid), or unallocated
   !> @param[in]     g0         first grid point of the block
   !> @param[in]     g1         last grid point of the block
   !> @param[out]    u_sr       short-range potential of the block (n, nv)
   !> @param[out]    ur_lr      long-range real-space potential of the block (n, nv)
   !> @param[in,out] sc         scratch of this thread
   !> @param[in,out] bad        lowest solute atom too close to a grid point of this thread
   subroutine forward_block(xyz, coord_u, inp, alpha, phi_field, g0, g1, u_sr, ur_lr, sc, bad)
      !> Grid points of the block, bohr (3, n)
      real(wp), intent(in) :: xyz(:, :)
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved inputs
      type(table_3d_inputs_type), intent(in) :: inp
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Summed field of the field-supplying terms (ngrid), or unallocated
      real(wp), allocatable, intent(in) :: phi_field(:)
      !> First grid point of the block
      integer, intent(in) :: g0
      !> Last grid point of the block
      integer, intent(in) :: g1
      !> Short-range potential of the block (n, nv)
      real(wp), intent(out) :: u_sr(:, :)
      !> Long-range real-space potential of the block (n, nv)
      real(wp), intent(out) :: ur_lr(:, :)
      !> Scratch of this thread
      type(block_scratch_type), intent(inout) :: sc
      !> Lowest solute atom too close to a grid point
      integer, intent(inout) :: bad
      integer :: n, a, v, i
      n = size(xyz, 2)
      associate (dist => sc%dist(1:n, :), u => sc%u(1:n), acc => sc%acc(1:n), phi => sc%phi(1:n), &
            & phi_lr => sc%phi_lr(1:n), f => sc%f(1:n))
         do a = 1, inp%nu
            do i = 1, n
               dist(i, a) = norm2(xyz(:, i) - coord_u(:, a))
            end do
         end do
         if (inp%check_separation) then
            do a = 1, inp%nu
               if (minval(dist(:, a)) < min_distance) then
                  bad = min(bad, a)
                  return
               end if
            end do
         end if
         if (inp%do_pair) then
            do v = 1, inp%nv
               acc = 0.0_wp
               do a = 1, inp%nu
                  call lj_12_6_evaluate(inp%sigma_uv(a, v), inp%eps_uv(a, v), dist(:, a), u)
                  acc = acc + u
               end do
               u_sr(:, v) = acc
            end do
         else
            u_sr = 0.0_wp
         end if
         ur_lr = 0.0_wp
         if (.not. inp%do_elec) return
         phi = 0.0_wp
         if (inp%has_point) then
            do a = 1, inp%nu
               phi = phi + inp%q_point(a)/dist(:, a)
            end do
         end if
         if (allocated(phi_field)) phi = phi + phi_field(g0:g1)
         phi_lr = 0.0_wp
         if (inp%tails) then
            do a = 1, inp%nu
               call ng_real_tail(dist(:, a), alpha, f)
               phi_lr = phi_lr + inp%q_tail(a)*f
            end do
         end if
         do v = 1, inp%nv
            ur_lr(:, v) = inp%q_v(v)*phi_lr
            u_sr(:, v) = (u_sr(:, v) + inp%q_v(v)*phi) - ur_lr(:, v)
         end do
      end associate
   end subroutine forward_block

   !> Form adjoints of `forward_block`
   !>
   !> - Radial piece p(|r_g - r_a|): s = w p'(d)/d
   !> - Geometry with `geom`: add s (r_g - r_a) to grid, subtract from atom
   !> - Charge adjoints and w_phi in every case
   !>
   !> @param[in]     xyz       grid points of the block, bohr (3, n)
   !> @param[in]     coord_u   solute coordinates, bohr (3, nu)
   !> @param[in]     inp       resolved inputs
   !> @param[in]     alpha     ng split parameter, 1/bohr
   !> @param[in]     geom      whether the geometric adjoint is accumulated
   !> @param[in]     w_sr      adjoint of u_sr of the block (n, nv)
   !> @param[in]     w_lr      adjoint of ur_lr of the block (n, nv)
   !> @param[out]    w_phi     adjoint of the full field of the block (n)
   !> @param[in,out] w_xyz     grid-point adjoint of the block (3, n)
   !> @param[in,out] grad      solute gradient of this thread (3, nu)
   !> @param[in,out] wq_point  adjoint of the point-field charges of this thread (nu)
   !> @param[in,out] wq_tail   adjoint of the tail charges of this thread (nu)
   !> @param[in,out] sc        scratch of this thread
   !> @param[in,out] bad       lowest solute atom too close to a grid point of this thread
   subroutine adjoint_block(xyz, coord_u, inp, alpha, geom, w_sr, w_lr, w_phi, w_xyz, grad, wq_point, wq_tail, &
         & sc, bad)
      !> Grid points of the block, bohr (3, n)
      real(wp), intent(in) :: xyz(:, :)
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved inputs
      type(table_3d_inputs_type), intent(in) :: inp
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Whether the geometric adjoint is accumulated
      logical, intent(in) :: geom
      !> Adjoint of u_sr of the block (n, nv)
      real(wp), intent(in) :: w_sr(:, :)
      !> Adjoint of ur_lr of the block (n, nv)
      real(wp), intent(in) :: w_lr(:, :)
      !> Adjoint of the full field of the block (n)
      real(wp), intent(out) :: w_phi(:)
      !> Grid-point adjoint of the block (3, n)
      real(wp), intent(inout) :: w_xyz(:, :)
      !> Solute gradient of this thread (3, nu)
      real(wp), intent(inout) :: grad(:, :)
      !> Adjoint of the point-field charges of this thread (nu)
      real(wp), intent(inout) :: wq_point(:)
      !> Adjoint of the tail charges of this thread (nu)
      real(wp), intent(inout) :: wq_tail(:)
      !> Scratch of this thread
      type(block_scratch_type), intent(inout) :: sc
      !> Lowest solute atom too close to a grid point
      integer, intent(inout) :: bad
      real(wp) :: inv
      integer :: n, a, v, i
      n = size(xyz, 2)
      associate (diff => sc%diff(:, 1:n), dist => sc%dist(1:n, 1), u => sc%u(1:n), du => sc%du(1:n), &
            & w_tail => sc%w_tail(1:n), s => sc%s(1:n), f => sc%f(1:n), g => sc%g(1:n))
         w_phi = 0.0_wp
         w_tail = 0.0_wp
         if (inp%do_elec) then
            do v = 1, inp%nv
               w_phi = w_phi + inp%q_v(v)*w_sr(:, v)
            end do
            if (inp%tails) then
               do v = 1, inp%nv
                  w_tail = w_tail + inp%q_v(v)*(w_lr(:, v) - w_sr(:, v))
               end do
            end if
         end if
         do a = 1, inp%nu
            do i = 1, n
               diff(:, i) = xyz(:, i) - coord_u(:, a)
               dist(i) = norm2(diff(:, i))
            end do
            if (inp%check_separation) then
               if (minval(dist) < min_distance) then
                  bad = min(bad, a)
                  return
               end if
            end if
            s = 0.0_wp
            if (geom .and. inp%do_pair) then
               do v = 1, inp%nv
                  call lj_12_6_evaluate(inp%sigma_uv(a, v), inp%eps_uv(a, v), dist, u, du)
                  s = s + w_sr(:, v)*du/dist
               end do
            end if
            if (inp%has_point) then
               do i = 1, n
                  inv = 1.0_wp/dist(i)
                  wq_point(a) = wq_point(a) + w_phi(i)*inv
                  s(i) = s(i) - w_phi(i)*inp%q_point(a)*inv*inv*inv
               end do
            end if
            if (inp%tails) then
               if (geom) then
                  call ng_real_tail(dist, alpha, f, g)
                  do i = 1, n
                     wq_tail(a) = wq_tail(a) + w_tail(i)*f(i)
                     s(i) = s(i) + w_tail(i)*inp%q_tail(a)*g(i)
                  end do
               else
                  call ng_real_tail(dist, alpha, f)
                  do i = 1, n
                     wq_tail(a) = wq_tail(a) + w_tail(i)*f(i)
                  end do
               end if
            end if
            if (.not. geom) cycle
            do i = 1, n
               w_xyz(:, i) = w_xyz(:, i) + s(i)*diff(:, i)
               grad(:, a) = grad(:, a) - s(i)*diff(:, i)
            end do
         end do
      end associate
   end subroutine adjoint_block

   !> Form the analytic reciprocal tail, phased over solute centres
   !>
   !> @param[in]     grid     volume grid with reciprocal points
   !> @param[in]     coord_u  solute coordinates, bohr (3, nu)
   !> @param[in]     inp      resolved inputs with tails on
   !> @param[in]     alpha    ng split parameter, 1/bohr
   !> @param[in,out] uk_lr    long-range reciprocal potential (npts_k, nv)
   subroutine reciprocal_tail(grid, coord_u, inp, alpha, uk_lr)
      !> Volume grid with reciprocal points
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved inputs with tails on
      type(table_3d_inputs_type), intent(in) :: inp
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Long-range reciprocal potential (npts_k, nv)
      complex(wp), intent(inout) :: uk_lr(:, :)
      real(wp) :: r0(3), kvec(3), env, kdotr, re_sum, im_sum, k_phase
      integer :: j, a, v, nt
      r0 = grid%point(1)
      nt = grid%team_size()
      !$omp parallel do default(shared) schedule(static) num_threads(nt) &
      !$omp private(j, kvec, env, kdotr, re_sum, im_sum, k_phase, a, v)
      do j = 1, grid%npts_k
         kvec = grid%kpoint(j)
         env = ng_fourier_envelope(dot_product(kvec, kvec), alpha)
         re_sum = 0.0_wp
         im_sum = 0.0_wp
         do a = 1, inp%nu
            kdotr = dot_product(kvec, coord_u(:, a) - r0)
            re_sum = re_sum + inp%q_tail(a)*cos(kdotr)
            im_sum = im_sum - inp%q_tail(a)*sin(kdotr)
         end do
         do v = 1, inp%nv
            k_phase = inp%q_v(v)*env
            uk_lr(j, v) = cmplx(k_phase*re_sum, k_phase*im_sum, kind=wp)
         end do
      end do
      !$omp end parallel do
   end subroutine reciprocal_tail

   !> Form adjoints of `reciprocal_tail`
   !>
   !> - C_j = E(k_j) sum_v q_v conj(w_k(j, v))
   !> - e_aj = exp(-i k_j.(r_a - r0))
   !> - With `geom`: dL/dr_a = q_a sum_j Im(C_j e_aj) k_j
   !> - dL/dq_a = sum_j Re(C_j e_aj)
   !>
   !> @param[in]     grid       volume grid with reciprocal points
   !> @param[in]     coord_u    solute coordinates, bohr (3, nu)
   !> @param[in]     inp        resolved inputs with tails on
   !> @param[in]     alpha      ng split parameter, 1/bohr
   !> @param[in]     geom       whether the solute gradient is accumulated
   !> @param[in]     w_k        adjoint of uk_lr (npts_k, nv)
   !> @param[in]     nt         team size, the last extent of the per-thread buffers
   !> @param[in,out] grad_t     per-thread solute gradient (3, nu, nt)
   !> @param[in,out] wq_tail_t  per-thread tail charge adjoint (nu, nt)
   subroutine reciprocal_tail_adjoint(grid, coord_u, inp, alpha, geom, w_k, nt, grad_t, wq_tail_t)
      !> Volume grid with reciprocal points
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Solute coordinates, bohr (3, nu)
      real(wp), intent(in) :: coord_u(:, :)
      !> Resolved inputs with tails on
      type(table_3d_inputs_type), intent(in) :: inp
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Whether the solute gradient is accumulated
      logical, intent(in) :: geom
      !> Adjoint of uk_lr (npts_k, nv)
      complex(wp), intent(in) :: w_k(:, :)
      !> Team size
      integer, intent(in) :: nt
      !> Per-thread solute gradient (3, nu, nt)
      real(wp), intent(inout) :: grad_t(:, :, :)
      !> Per-thread tail charge adjoint (nu, nt)
      real(wp), intent(inout) :: wq_tail_t(:, :)
      real(wp) :: r0(3), kvec(3), env, kdotr
      complex(wp) :: c, ce
      integer :: j, a, v, it
      r0 = grid%point(1)
      !$omp parallel do default(shared) schedule(static) num_threads(nt) &
      !$omp private(j, kvec, env, kdotr, c, ce, a, v, it)
      do j = 1, grid%npts_k
         it = 1
         !$ it = omp_get_thread_num() + 1
         kvec = grid%kpoint(j)
         env = ng_fourier_envelope(dot_product(kvec, kvec), alpha)
         if (env == 0.0_wp) cycle
         c = (0.0_wp, 0.0_wp)
         do v = 1, inp%nv
            c = c + inp%q_v(v)*conjg(w_k(j, v))
         end do
         c = env*c
         do a = 1, inp%nu
            kdotr = dot_product(kvec, coord_u(:, a) - r0)
            ce = c*cmplx(cos(kdotr), -sin(kdotr), kind=wp)
            wq_tail_t(a, it) = wq_tail_t(a, it) + real(ce, wp)
            if (geom) grad_t(:, a, it) = grad_t(:, a, it) + inp%q_tail(a)*aimag(ce)*kvec
         end do
      end do
      !$omp end parallel do
   end subroutine reciprocal_tail_adjoint

end module moist_model_moz_potential_kernel_table_3d
