!> DROP projection driver: two chunked passes over per-thread workers
!>
!> Chunks of `projection.batch_size` anchors, one OpenMP loop per pass,
!> `schedule(dynamic)`; every batch holds the same anchors at any thread count:
!>
!> 1. local pass: local Newton and next step per anchor
!>    ([[projection_worker_type]]); the chunk's [[local_pass_type]] kept
!> 2. lattice seeder built once from the geometry and the largest search
!>    ball of every seeded anchor ([[search_geometry]]), its lattice along the
!>    grid model's lattice axes (molecular frame or point group)
!> 3. search pass, heaviest chunks first ([[search_order]]): each worker seeds
!>    from the built lattice, solves and finalizes into its row buffer
!>
!> Then row counts prefix-summed into `offset`, buffers scattered into the CSR
!> row arrays
!>
!> - rows below the weight cut never appended, no post-filter
!> - anchor without an accepted point lost, no row; warned when it had seeds,
!>   only reported without any seed up to the reach (no surface there)
!> - LSF failure, singular tangent block: abort through
!>   [[drop_abort_latch_type]] with the first failing chunk's error in pass
!>   order (anchor order, heaviest first in the search pass) at any thread
!>   count; chunks behind it skipped, those ahead still run
!> - warning when the reach gate removes more than [[reach_warn_share]] of the
!>   kept rows' area (`a/f_reach` without the gate); surface beyond the closed
!>   edge is never searched, so more is likely lost
submodule(moist_cavity_drop) moist_cavity_drop_projection
   !$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_threads, only: drop_abort_latch_type
   use moist_utils_prettyprint, only: prettyprinter
   use moist_cavity_drop_projection_worker, only: projection_worker_type, local_pass_type, search_geometry, &
                                                  search_order
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seed_build_type
   implicit none(type, external)

   !> Area share removed by the reach gate above which the build warns
   real(wp), parameter :: reach_warn_share = 1.0e-3_wp

contains

   !* ================================================================================= *!
   !*                                     Projection                                    *!
   !* ================================================================================= *!

   !> Project every anchor onto the level set surface into the CSR row arrays
   !>
   !> @param[inout] self  Cavity after the pre-filter (anchor phase)
   !> @param[out]   error Build-fatal projection failure
   module subroutine project_all_points(self, error)
      class(cavity_type_drop), intent(inout), target :: self
      type(error_type), allocatable, intent(out) :: error

      type(prettyprinter) :: pp
      type(drop_abort_latch_type) :: abort
      type(projection_worker_type), allocatable :: workers(:)
      !> Local pass of every chunk, kept for its search pass
      type(local_pass_type), allocatable :: local(:)
      !> Seeder of the level, built between the passes and read by every worker
      type(projection_seeder_type) :: seeder
      type(seed_build_type) :: geometry
      type(error_type), allocatable :: local_error
      integer, allocatable :: nrows_of(:), chunk_order(:)
      logical, allocatable :: lost(:), unseeded(:)
      integer :: nmax_anchor, nthreads, thread_slot, chunk, nchunk, idx, ichunk, i, i0, i1, ipass
      integer :: ithread, irow, ia, dest, k, nout, number_base, n_lost, n_unseeded, n_seeded, n_escalated, n_captured
      integer :: n_extended, n_full, n_seeds, n_skipped, n_window_skipped
      real(wp) :: pdist_search_max
      !> Total area, the area the reach gate removes from the rows, one row's
      !> area and reach gate
      real(wp) :: area_tot, area_cut, area_row, f_reach
      logical :: has_seeded
      character(len=*), parameter :: pass_name(2) = [character(len=16) :: "Local projection", "Branch search"]
      character(len=160) :: report

      nmax_anchor = self%anchor%ngrid
      nthreads = self%ctx%get_num_threads()
      chunk = max(1, self%param%batch_size)
      nchunk = (nmax_anchor + chunk - 1)/chunk

      allocate (workers(nthreads))
      allocate (local(nchunk))
      allocate (chunk_order(nchunk))
      do ichunk = 1, nchunk
         chunk_order(ichunk) = ichunk
      end do
      allocate (nrows_of(nmax_anchor), source=0)
      allocate (lost(nmax_anchor), source=.false.)
      allocate (unseeded(nmax_anchor), source=.false.)
      call abort%reset()

      call seeder%init(self%param%seed, error)
      if (allocated(error)) return

      pp = self%ctx%printer(level=4)
      call pp%blank()
      call pp%push("Projection:")

      has_seeded = .false.
      do ipass = 1, 2
         if (ipass == 2) then
            call search_geometry(local, self%mol%xyz, self%radii, self%param%branch%reach_off, nthreads, geometry)
            ! The seed lattice turns with the sphere grids (lab axes without a grid model)
            if (allocated(self%grid_model)) call self%grid_model%lattice_axes(geometry%axes)
            has_seeded = size(geometry%radius_max) > 0
            if (has_seeded) then
               call search_order(local, chunk_order)
               call self%ctx%timer%start("Seed lattice")
               call seeder%build(self%lsf_model, geometry, error)
               call self%ctx%timer%stop("Seed lattice")
               if (allocated(error)) return
               call pp%step("Seed lattice", self%ctx%timer%last())
               if (self%ctx%verbosity(4)) call seeder%report(pp)
            end if
         end if

         call self%ctx%timer%start(trim(pass_name(ipass)))

         !$omp parallel num_threads(nthreads) default(shared) private(thread_slot, idx, ichunk, i, i0, i1, &
         !$omp& local_error)
         thread_slot = 1
         !$ thread_slot = omp_get_thread_num() + 1

         ! First-touch init inside the region so the worker scratch is core-local
         if (.not. allocated(workers(thread_slot)%lsf)) then
            call workers(thread_slot)%init(self%param, self%lsf_model, self%mol%xyz, local_error)
            if (allocated(local_error)) call abort%latch_error(local_error, 0)
         end if
         ! Shared setup failure on any thread before the first chunk starts
         !$omp barrier

         !$omp do schedule(dynamic, 1)
         do idx = 1, nchunk
            ichunk = chunk_order(idx)
            ! Failures latch by their place in the pass's order
            if (abort%skips(idx)) cycle
            i0 = (ichunk - 1)*chunk + 1
            i1 = min(nmax_anchor, ichunk*chunk)
            if (ipass == 1) then
               call workers(thread_slot)%solve_local(self%ctx, [(i, i=i0, i1)], self%anchor%xyz(:, i0:i1), &
                                                     self%anchor%owner(i0:i1), self%anchor%wleb0(i0:i1), &
                                                     self%anchor%f0(i0:i1), local(ichunk), local_error)
            else
               call workers(thread_slot)%search(self%ctx, local(ichunk), seeder, [(i, i=i0, i1)], &
                                                self%anchor%xyz(:, i0:i1), self%anchor%owner(i0:i1), &
                                                self%anchor%wleb0(i0:i1), self%anchor%f0(i0:i1), &
                                                nrows_of(i0:i1), lost(i0:i1), local_error, &
                                                unseeded=unseeded(i0:i1))
            end if
            if (allocated(local_error)) then
               call abort%latch_error(local_error, idx)
               cycle
            end if

         end do
         !$omp end do
         !$omp end parallel

         call self%ctx%timer%stop(trim(pass_name(ipass)))
         if (.not. abort%is_set()) call pp%step(trim(pass_name(ipass)), self%ctx%timer%last())

         if (abort%is_set()) exit
      end do
      call pp%pop()

      ! Worker clones may borrow read-only data of the level set model: release
      ! them before the model can change
      do ithread = 1, size(workers)
         if (allocated(workers(ithread)%lsf)) deallocate (workers(ithread)%lsf)
      end do

      if (abort%is_set()) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call fatal_error(error, "Projection aborted without an error. (unreachable in normal execution)")
         end if
         return
      end if

      if (self%ctx%do_profile) call self%ctx%timer%start("Scatter")

      !> Row layout: prefix sum of the per-anchor row counts
      if (allocated(self%offset)) deallocate (self%offset)
      allocate (self%offset(nmax_anchor + 1))
      self%offset(1) = 0
      do ia = 1, nmax_anchor
         self%offset(ia + 1) = self%offset(ia) + nrows_of(ia)
      end do
      nout = self%offset(nmax_anchor + 1)

      ! Row arrays of the previous surface were freed before the anchor setup
      allocate (self%xyz(3, nout), self%normal0(3, nout), self%wleb(nout), self%lambda0(nout), &
                self%f(nout), self%pdist(nout), self%wbranch(nout), self%phi0(nout), &
                self%owner(nout), self%comp(nout), self%branch(nout), self%anchor_of_row(nout), self%cpjac_scal0(nout), &
                self%w_f0(nout), self%numbering(nout))
      self%ngrid = nout

      !> Scatter each worker's rows to their anchor's CSR slots; rows of one
      !> anchor are contiguous in the buffer that produced them
      number_base = max(1, self%grid%nnode)
      !$omp parallel do num_threads(nthreads) default(shared) private(ithread, irow, ia, k, dest) &
      !$omp& schedule(static)
      do ithread = 1, nthreads
         k = 0
         do irow = 1, workers(ithread)%rows%n
            ia = workers(ithread)%rows%anchor(irow)
            if (irow == 1) then
               k = 0
            else if (workers(ithread)%rows%anchor(irow - 1) /= ia) then
               k = 0
            end if
            k = k + 1
            dest = self%offset(ia) + k
            self%xyz(:, dest) = workers(ithread)%rows%xyz(:, irow)
            self%normal0(:, dest) = workers(ithread)%rows%normal(:, irow)
            self%wleb(dest) = workers(ithread)%rows%wleb(irow)
            self%lambda0(dest) = workers(ithread)%rows%lambda(irow)
            self%pdist(dest) = workers(ithread)%rows%pdist(irow)
            self%wbranch(dest) = workers(ithread)%rows%wbranch(irow)
            self%phi0(dest) = workers(ithread)%rows%phi(irow)
            self%branch(dest) = workers(ithread)%rows%branch(irow)
            self%cpjac_scal0(dest) = workers(ithread)%rows%cpjac(irow)
            self%w_f0(dest) = workers(ithread)%rows%w_f(irow)
            self%anchor_of_row(dest) = ia
            self%owner(dest) = self%anchor%owner(ia)
            self%comp(dest) = self%anchor%comp(ia)
            ! The row's switching value carries its share of the anchor
            self%f(dest) = self%anchor%f0(ia)*workers(ithread)%rows%wbranch(irow)
            self%numbering(dest) = self%anchor%id(ia) + number_base*(workers(ithread)%rows%branch(irow) - 1)
         end do
      end do
      !$omp end parallel do

      n_seeded = 0
      n_escalated = 0
      n_captured = 0
      n_extended = 0
      n_seeds = 0
      n_skipped = 0
      n_full = 0
      n_window_skipped = 0
      pdist_search_max = 0.0_wp
      do ithread = 1, nthreads
         n_seeded = n_seeded + workers(ithread)%n_seeded
         n_escalated = n_escalated + workers(ithread)%n_escalated
         n_captured = n_captured + workers(ithread)%n_captured
         n_extended = n_extended + workers(ithread)%n_extended
         n_seeds = n_seeds + workers(ithread)%n_seeds
         n_skipped = n_skipped + workers(ithread)%n_skipped
         n_full = n_full + workers(ithread)%n_full
         n_window_skipped = n_window_skipped + workers(ithread)%n_window_skipped
         pdist_search_max = max(pdist_search_max, workers(ithread)%pdist_search_max)
      end do

      self%n_searched = n_seeded
      self%n_extended = n_extended
      self%n_full_reach = n_full
      self%pdist_search_max = pdist_search_max

      ! Lost anchors are reported, never silently dropped; only seeds that all
      ! failed are a defect, an anchor without surface in reach is not
      n_unseeded = count(unseeded)
      n_lost = count(lost) - n_unseeded
      if (n_lost > 0) then
         write (report, "(i0,a,i0,a)") n_lost, " of ", nmax_anchor, &
            " anchors had seeds but no converged surface point and were dropped"
         call self%ctx%warning("DROP projection: "//trim(report))
      end if
      if (n_lost + n_unseeded > 0 .and. self%ctx%verbosity(6)) then
         do ia = 1, nmax_anchor
            if (.not. lost(ia)) cycle
            if (unseeded(ia)) then
               write (report, "(a,i0,a,3es12.4)") "Anchor without surface in reach ", self%anchor%id(ia), &
                  ": ", self%anchor%xyz(:, ia)
            else
               write (report, "(a,i0,a,3es12.4)") "Lost anchor ", self%anchor%id(ia), ": ", &
                  self%anchor%xyz(:, ia)
            end if
            call self%ctx%message(trim(report), level=6)
         end do
      end if

      ! Area the reach gate takes from the rows it keeps; surface beyond the
      ! closed edge is never searched, so more is likely lost
      area_tot = 0.0_wp
      area_cut = 0.0_wp
      do i = 1, nout
         area_row = self%wleb(i)*self%radii(self%owner(i))**2*self%f(i)
         area_tot = area_tot + area_row
         f_reach = self%f_reach%f0(self%pdist(i))
         if (f_reach > 0.0_wp .and. f_reach < 1.0_wp) area_cut = area_cut + area_row*(1.0_wp/f_reach - 1.0_wp)
      end do
      if (area_cut > reach_warn_share*area_tot) then
         write (report, "(a,f0.2,a,f0.2,a,f0.2,a)") "the reach gate removes ", 100.0_wp*area_cut/area_tot, &
            " % of the area between reach_on = ", self%param%branch%reach_on, " and reach_off = ", &
            self%param%branch%reach_off, " bohr"
         call self%ctx%warning("DROP projection: "//trim(report)//"; consider a longer reach")
      end if

      if (self%ctx%verbosity(4)) then
         call pp%blank()
         call pp%push("Projection results:")
         call pp%kv("Total points", nout)
         call pp%kv("Anchors searched", n_seeded)
         call pp%kv("Anchors escalated", n_escalated)
         call pp%kv("Extended searches", n_extended)
         call pp%kv("Searched to the reach", n_full)
         if (self%param%proj_level == 2) call pp%kv("Close windows not seeded", n_window_skipped)
         call pp%kv("Largest search radius", pdist_search_max, "Bohr")
         if (area_tot > 0.0_wp) call pp%kv("Area removed by the reach gate", 100.0_wp*area_cut/area_tot, "%")
         call pp%kv("Seeds solved", n_seeds)
         call pp%kv("Seeds solved before, skipped", n_skipped)
         call pp%kv("Duplicate seeds stopped", n_captured)
         call pp%kv("Anchors with branches", count(nrows_of > 1))
         call pp%kv("Branched points", sum(nrows_of, mask=nrows_of > 1))
         call pp%kv("Anchors without surface in reach", n_unseeded)
         call pp%kv("Anchors lost", n_lost)
         if (seeder%n_unresolved > 0) call pp%kv("Seeder cells unresolved", seeder%n_unresolved)
         call pp%pop()
      end if

      if (self%ctx%do_profile) call self%ctx%timer%stop("Scatter")

   end subroutine project_all_points

end submodule moist_cavity_drop_projection
