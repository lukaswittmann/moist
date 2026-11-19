!> DROP projection driver: chunked loop over per-thread workers
!>
!> One OpenMP loop, `schedule(dynamic)` over chunks of `projection.batch_size`
!> anchors; each thread's [[projection_worker_type]] runs local Newton,
!> gate of the local point, seeder, batched Newton and finalize on its chunk and appends rows
!> to its own buffer. Afterwards the per-anchor row counts are prefix-summed
!> into `offset` and the buffers are scattered into the CSR row arrays
!>
!> - rows below the weight cut are never appended, so no post-filter follows
!> - an anchor without an accepted point is lost and keeps no row; with seeds
!>   that all failed it is counted in a warning, without any seed up to the
!>   reach (no surface there) only in the report
!> - LSF failures, an exhausted octree budget and a singular tangent block
!>   abort the build through [[drop_abort_latch_type]], lowest anchor first
!> - a warning when the reach gate removes more than [[reach_warn_share]] of
!>   the area of the rows it keeps: their area without the gate is
!>   `a/f_reach`; surface beyond the closed edge is never searched, so more of
!>   it is likely lost
submodule(moist_cavity_drop) moist_cavity_drop_projection
   !$ use omp_lib, only: omp_get_thread_num, omp_get_wtime
   use, intrinsic :: iso_fortran_env, only: int64
   use moist_cavity_drop_threads, only: drop_abort_latch_type
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_cavity_drop_projection_worker, only: projection_worker_type
   implicit none(type, external)

   !> Area share removed by the reach gate above which the build warns
   real(wp), parameter :: reach_warn_share = 1.0e-3_wp

contains

   !> Wall-clock seconds for the projector's progress report
   !>
   !> `omp_get_wtime` when the build has OpenMP, `system_clock` otherwise, so a
   !> serial build neither references the OpenMP runtime nor loses its timings
   function drop_wall_time() result(seconds)
      !> Elapsed wall time in seconds, measured from an unspecified origin
      real(wp) :: seconds

      !> Clock reading and tick rate of the serial fallback
      integer(int64) :: ticks, rate

      call system_clock(ticks, rate)
      seconds = 0.0_wp
      if (rate > 0_int64) seconds = real(ticks, wp)/real(rate, wp)
      !$ seconds = omp_get_wtime()
   end function drop_wall_time

   !* ================================================================================= *!
   !*                                     Projection                                    *!
   !* ================================================================================= *!

   !> Project every anchor onto the level set surface into the CSR row arrays
   !>
   !> @param[inout] self  Cavity after the pre-filter (anchor phase)
   !> @param[out]   error Build-fatal projection failure
   module subroutine project_all_points(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      type(prettyprinter) :: pp
      type(prettylistprinter) :: plp
      type(drop_abort_latch_type) :: abort
      type(projection_worker_type), allocatable :: workers(:)
      type(error_type), allocatable :: local_error
      integer, allocatable :: nrows_of(:)
      logical, allocatable :: lost(:), unseeded(:)
      integer :: nmax_anchor, nthreads, thread_slot, chunk, nchunk, ichunk, i, i0, i1
      integer :: ithread, irow, ia, dest, k, nout, number_base, n_lost, n_unseeded, n_seeded, n_escalated, n_captured
      integer :: n_extended, n_full, n_restarted
      real(wp) :: pdist_search_max
      !> Total area, the area the reach gate removes from the rows, one row's
      !> area and reach gate
      real(wp) :: area_tot, area_cut, area_row, f_reach
      integer :: chunks_done, next_progress_pct, est_done
      real(wp) :: wall_start, wall_end, now_time, progress_last_time, progress_elapsed, progress_eta
      logical :: trigger_time, trigger_pct
      character(len=160) :: report

      nmax_anchor = self%anchor%ngrid
      nthreads = self%ctx%get_num_threads()
      chunk = max(1, self%param%batch_size)
      nchunk = (nmax_anchor + chunk - 1)/chunk

      allocate (workers(nthreads))
      allocate (nrows_of(nmax_anchor), source=0)
      allocate (lost(nmax_anchor), source=.false.)
      allocate (unseeded(nmax_anchor), source=.false.)
      call abort%reset()

      if (self%ctx%verbosity(4)) then
         plp = new_prettylistprinter([12, 8, 12, 12], &
                                     [character(len=8) :: "Point", "% done", "Elapsed (s)", "Left (s)"], &
                                     unit=self%ctx%unit)
         call plp%blank()
         call plp%header("DROP Projector")
         call plp%print_header()
         call plp%separator()
      end if

      wall_start = drop_wall_time()

      !$omp parallel num_threads(nthreads) default(shared) private(thread_slot, ichunk, i, i0, i1, &
      !$omp& local_error, chunks_done, next_progress_pct, est_done, now_time, progress_last_time, &
      !$omp& progress_elapsed, progress_eta, trigger_time, trigger_pct)
      thread_slot = 1
      !$ thread_slot = omp_get_thread_num() + 1

      ! First-touch init inside the region so the worker scratch is core-local
      call workers(thread_slot)%init(self%param, self%lsf_model, self%mol%xyz, local_error, &
                                     debug=self%ctx%debug, unit=self%ctx%unit)
      if (allocated(local_error)) call abort%latch_error(local_error, 0)

      chunks_done = 0
      next_progress_pct = 10
      progress_last_time = wall_start

      ! Publish a setup failure on any thread before the first chunk starts
      !$omp barrier

      !$omp do schedule(dynamic)
      do ichunk = 1, nchunk
         if (abort%requested) cycle
         i0 = (ichunk - 1)*chunk + 1
         i1 = min(nmax_anchor, ichunk*chunk)
         call workers(thread_slot)%process_chunk(self%ctx, [(i, i=i0, i1)], self%anchor%xyz(:, i0:i1), &
                                                 self%anchor%owner(i0:i1), self%anchor%wleb0(i0:i1), &
                                                 self%anchor%f0(i0:i1), nrows_of(i0:i1), lost(i0:i1), &
                                                 local_error, unseeded=unseeded(i0:i1))
         if (allocated(local_error)) then
            call abort%latch_error(local_error, i0)
            cycle
         end if

         ! Progress is reported by thread 1, extrapolating its own count
         chunks_done = chunks_done + 1
         if (thread_slot == 1 .and. self%ctx%verbosity(4)) then
            est_done = min(nmax_anchor, chunks_done*chunk*nthreads)
            now_time = drop_wall_time()
            trigger_time = (now_time - progress_last_time) >= 5.0_wp
            trigger_pct = nmax_anchor > 0 .and. (100*est_done >= next_progress_pct*nmax_anchor)
            if (trigger_time .or. trigger_pct) then
               progress_elapsed = now_time - wall_start
               progress_eta = 0.0_wp
               if (progress_elapsed > 0.0_wp .and. est_done > 0) then
                  progress_eta = real(max(0, nmax_anchor - est_done), wp)*progress_elapsed/real(est_done, wp)
               end if
               !$omp critical (projection_warning_print)
               call plp%begin_row()
               call plp%add(est_done)
               call plp%add(100.0_wp*real(est_done, wp)/real(max(1, nmax_anchor), wp), fmt="f6.2")
               call plp%add(progress_elapsed, fmt="f8.2")
               call plp%add(progress_eta, fmt="f8.2")
               call plp%end_row()
               !$omp end critical (projection_warning_print)
               progress_last_time = now_time
               do while (next_progress_pct <= 100 .and. 100*est_done >= next_progress_pct*nmax_anchor)
                  next_progress_pct = next_progress_pct + 10
               end do
            end if
         end if
      end do
      !$omp end do
      !$omp end parallel

      wall_end = drop_wall_time()

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call fatal_error(error, "Projection aborted without an error. (unreachable in normal execution)")
         end if
         return
      end if

      if (self%ctx%verbosity(4)) then
         call plp%separator()
         call plp%begin_row()
         call plp%add(nmax_anchor)
         call plp%add("")
         call plp%add(wall_end - wall_start, fmt="f8.2")
         call plp%add(0.0_wp, fmt="f8.2")
         call plp%end_row()
      end if

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
                self%owner(nout), self%branch(nout), self%anchor_of_row(nout), self%cpjac_scal0(nout), &
                self%w_f0(nout), self%numbering(nout))
      self%ngrid = nout

      !> Scatter each worker's rows to their anchor's CSR slots; rows of one
      !> anchor are contiguous in the buffer that produced them
      number_base = max(1, self%nsph*self%param%num_leb)
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
      n_restarted = 0
      n_full = 0
      pdist_search_max = 0.0_wp
      do ithread = 1, nthreads
         n_seeded = n_seeded + workers(ithread)%n_seeded
         n_escalated = n_escalated + workers(ithread)%n_escalated
         n_captured = n_captured + workers(ithread)%n_captured
         n_extended = n_extended + workers(ithread)%n_extended
         n_restarted = n_restarted + workers(ithread)%n_restarted
         n_full = n_full + workers(ithread)%n_full
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
      if (n_lost + n_unseeded > 0 .and. self%ctx%verbosity(5)) then
         do ia = 1, nmax_anchor
            if (.not. lost(ia)) cycle
            if (unseeded(ia)) then
               write (report, "(a,i0,a,3es12.4)") "Anchor without surface in reach ", self%anchor%id(ia), &
                  ": ", self%anchor%xyz(:, ia)
            else
               write (report, "(a,i0,a,3es12.4)") "Lost anchor ", self%anchor%id(ia), ": ", &
                  self%anchor%xyz(:, ia)
            end if
            call self%ctx%debug_message(trim(report))
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
         pp = new_prettyprinter(unit=self%ctx%unit)
         call pp%blank()
         call pp%push("Projection results:")
         call pp%kv("Total points", nout)
         call pp%kv("Anchors searched", n_seeded)
         call pp%kv("Anchors escalated", n_escalated)
         call pp%kv("Extended searches", n_extended)
         call pp%kv("Searched to the reach", n_full)
         call pp%kv("Largest search radius", pdist_search_max, "Bohr")
         if (area_tot > 0.0_wp) call pp%kv("Area removed by the reach gate", 100.0_wp*area_cut/area_tot, "%")
         call pp%kv("Duplicate seeds stopped", n_captured)
         call pp%kv("Anchors with branches", count(nrows_of > 1))
         call pp%kv("Branched points", sum(nrows_of, mask=nrows_of > 1))
         call pp%kv("Anchors solved again without cheap steps", n_restarted)
         call pp%kv("Anchors without surface in reach", n_unseeded)
         call pp%kv("Anchors lost", n_lost)
         call pp%pop()
      end if

   end subroutine project_all_points

end submodule moist_cavity_drop_projection
