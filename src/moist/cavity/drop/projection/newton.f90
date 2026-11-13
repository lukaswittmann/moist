!> Batched two-stage KKT Newton for the DROP closest-point projection
!>
!> Solves, for every pending column of a [[projection_batch_type]],
!>
!>     F(x, lambda) = [ w (x - a) - lambda grad S(x) ]  = 0,   w = phi_alpha
!>                    [         -S(x)                ]
!>
!>     J = [ w I - lambda Hess S   -grad S ]
!>         [      -grad S^T            0   ]
!>
!> - stage 1 drops `lambda Hess S`, a closed-form 4x4 needing only an order-1
!>   jet; it is exact at the anchor, where `lambda = 0`
!> - stage 2 is the full system, an unrolled pivoted 4x4 per point on an
!>   order-2 jet
!> - a point leaves stage 1 once `||F||` drops below `switch_fnorm`, inside its
!>   basin, or after `max_cheap_steps` cheap steps; far from the surface the
!>   full block `w I - lambda Hess S` can be indefinite and walk into a farther
!>   basin, the cheap direction cannot
!> - the direction is scaled so its spatial part is at most `max_step` long: at
!>   a weak-gradient anchor the first-order distance `S / ||grad S||` is
!>   unreliable, and an unscaled step lands outside the molecule where the
!>   screened level set is flat
!> - globalisation: backtracking on the SQP exact-penalty merit
!>   `phi = w/2 ||x - a||^2 + mu |S|`, `mu = max(mu, |lambda+| + merit_margin)`
!>   from `max(1, |lambda|)`; a full-system direction that is no merit descent
!>   (indefinite `w I - lambda Hess S` near saddles) is replaced by the cheap
!>   direction for that iteration, which always descends for `mu > |lambda+|`
!> - Maratos guard: once `||F|| <= switch_fnorm` the acceptance test is the
!>   `||F||` decrease, so the full steps near the solution are not rejected
!> - every iteration gathers the active set stage-1 first, stage-2 second, so a
!>   line-search round is one order-1 and one order-2 `evaluate_batch` call
!> - an accepted trial's jet is the next iterate's jet; nothing is
!>   re-evaluated except one order-2 refresh when a point changes stage, and
!>   one when it converges in stage 1, since the final jet is always order 2
!>
!> Stopping, per point:
!>   - converged: `||F|| < tol` or `max|dz| < max(0.1 tol, 8 eps max|z|)`
!>     (max-norm over x and lambda: the old nlesolver `tolx` test, floored at
!>     the round-off of the iterate so a tolerance below it cannot stall a point)
!>   - failed: `max_iter` steps, singular system, vanishing gradient, or an
!>     iterate beyond the level set's screening reach
!>   - captured, with `capture_duplicates`: an unscaled step whose prediction
!>     `x + dx` lies within `capture_distance` of an earlier accepted column of
!>     the same run; the column is a duplicate of that foot and stops unevaluated
!>
!> Capture:
!>   - captors are accepted columns of the same run, whose status no solve
!>     changes, so a decision depends neither on other runs nor on other
!>     pending columns, and the first column of a foot in run order survives
!>   - a merely converged column captures nothing: several starts can meet at
!>     a saddle that the polish escape sends to different feet
!>
!> Unreached columns (no atom within the screening reach, outputs unspecified):
!>   - a start point: `point_failed_unreached`, routine for far shell seeds; an
!>     unreached anchor itself is build-fatal
!>   - a line-search trial: rejected and backtracked; at the step floor the
!>     column is `point_failed_unreached` instead of the unconditional accept
!>
module moist_cavity_drop_projection_newton
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_math_packed, only: packed_index
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, &
                                                point_captured, &
                                                point_failed_iterations, point_failed_singular, &
                                                point_failed_gradient, &
                                                point_failed_unreached
   implicit none(type, external)
   private

   public :: newton_config_type, newton_scratch_type, newton_solve_batch
   public :: newton_solve_4x4

   !> Armijo coefficient of the backtracking on `||F||`
   real(wp), parameter :: armijo_c = 1.0e-4_wp
   !> Squared gradient norm below which the constraint is degenerate
   real(wp), parameter :: grad2_floor = 1.0e-28_wp
   !> Distance (Bohr) within which a Newton prediction lands on an accepted foot
   !>
   !> - far above `branch_sep_cut` (default `1e-8`), within which the finalize
   !>   merges two feet, so a duplicate start stops rounds before it converges
   !> - well below the separation of distinct feet on near-degenerate arcs,
   !>   which can be a few `1e-3`
   real(wp), parameter :: capture_distance = 1.0e-4_wp

   !> Newton settings
   ! TODO: good idea to introduce a new type here with drop parameters again?
   type :: newton_config_type
      !> Anchor objective weight `w`
      real(wp) :: phi_alpha = 0.5_wp
      !> Residual tolerance; `0.1 tol` is the step tolerance
      real(wp) :: tol = 1.0e-10_wp
      !> Step budget per point
      integer :: max_iter = 150
      !> `||F||` below which a point switches to the full system
      real(wp) :: switch_fnorm = 1.0e-3_wp
      !> Cheap steps after which a point switches to the full system
      integer :: max_cheap_steps = 50
      !> Smallest line-search step, accepted unconditionally
      real(wp) :: ls_floor = 0.1_wp
      !> Longest spatial Newton step before the line search (Bohr)
      real(wp) :: max_step = 1.0_wp
      !> Margin of the merit penalty over the multiplier, `mu >= |lambda+| + margin`
      real(wp) :: merit_margin = 1.0_wp
      !> Stop a column whose prediction lands on an earlier accepted column of its run
      logical :: capture_duplicates = .false.
   end type newton_config_type

   !> Per-thread Newton scratch, column- and list-indexed, grown geometrically
   type :: newton_scratch_type
      !> Allocated columns
      integer :: capacity = 0
      !> Stage per column (1 cheap, 2 full)
      integer, allocatable :: stage(:)
      !> Cheap steps taken, steps taken
      integer, allocatable :: ncheap(:), iter(:)
      !> Steps taken in stage 1 and stage 2
      integer, allocatable :: n1(:), n2(:)
      !> Whether the column's stored jet includes the Hessian
      logical, allocatable :: has_hess(:)
      !> Column was pending when the solve started
      logical, allocatable :: processed(:)
      !> Current residual norm
      real(wp), allocatable :: fnorm(:)
      !> Newton direction and line-search step
      real(wp), allocatable :: dz(:, :), alpha(:)
      !> Merit penalty, merit at the iterate, merit directional derivative
      real(wp), allocatable :: mu(:), merit(:), dmerit(:)
      !> Polish: column escaped a saddle and needs the re-polish (grown by the polish)
      logical, allocatable :: escaped(:)
      !> Active list, line-search list, partition workspace
      integer, allocatable :: list(:), lsl(:), ptmp(:)
      !> Gathered trial points, multipliers and jets (list-indexed), packed Hessian `p2(6, :)`
      real(wp), allocatable :: px(:, :), plam(:), p0(:), p1(:, :), p2(:, :)
      !> Gathered point unreached, `lsf%unreached(k)` read right after its batch (list-indexed)
      logical, allocatable :: punreached(:)
      !> Previous accepted column of the same run, 0 for none (column-indexed)
      integer, allocatable :: prev_accepted(:)
   end type newton_scratch_type

contains

   !* ================================================================================= *!
   !*                                 Batched KKT Newton                                *!
   !* ================================================================================= *!

   !> Solve the KKT system for every pending column in `first:last`
   !>
   !> - columns start from their stored `x` and `lambda`; a column whose start
   !>   is exactly its run's anchor also records the anchor jet in the run
   !> - on return each processed column is `point_converged` with `S`, `gradS`,
   !>   `hessS`, `phi`, `rho` and `normal` of its final iterate, `point_captured`
   !>   (with `capture_duplicates`), or one of the `point_failed_*` codes
   !>
   !> @param[inout] batch        Points to solve
   !> @param[inout] lsf          Thread-local level set
   !> @param[in]    cfg          Newton settings
   !> @param[inout] scratch      Per-thread scratch
   !> @param[out]   error        LSF evaluation failure or unreached anchor (build-fatal)
   !> @param[in]    first        First column (default 1)
   !> @param[in]    last         Last column (default `batch%n`)
   !> @param[out]   stage1_steps Cheap steps per column (diagnostic)
   !> @param[out]   stage2_steps Full steps per column (diagnostic)
   subroutine newton_solve_batch(batch, lsf, cfg, scratch, error, first, last, &
                                 stage1_steps, stage2_steps)
      type(projection_batch_type), intent(inout) :: batch
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      type(newton_config_type), intent(in) :: cfg
      type(newton_scratch_type), intent(inout) :: scratch
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: first
      integer, intent(in), optional :: last
      integer, intent(inout), optional :: stage1_steps(:)
      integer, intent(inout), optional :: stage2_steps(:)

      integer :: j0, j1, j, k, nact, nkeep, nls, m1, m, irun
      real(wp) :: w, a(3), fnew, dzmax, dzfloor, step(4), merit_new, mu_full
      logical :: ok, accept
      character(len=16) :: buffer

      j0 = 1
      j1 = batch%n
      if (present(first)) j0 = first
      if (present(last)) j1 = last
      if (j1 < j0) return
      call reserve_scratch(scratch, j1)
      w = cfg%phi_alpha

      !> Start: one order-1 evaluation of every pending column
      nact = 0
      do j = j0, j1
         scratch%processed(j) = batch%status(j) == point_pending
         if (.not. scratch%processed(j)) cycle
         scratch%stage(j) = 1
         scratch%ncheap(j) = 0
         scratch%iter(j) = 0
         scratch%n1(j) = 0
         scratch%n2(j) = 0
         scratch%has_hess(j) = .false.
         scratch%mu(j) = max(1.0_wp, abs(batch%lambda(j)))
         nact = nact + 1
         scratch%list(nact) = j
      end do
      if (nact == 0) return
      if (cfg%capture_duplicates) then
         call link_accepted(batch, batch%run_start(batch%run_of(j0)), j1, scratch%prev_accepted)
      end if

      do k = 1, nact
         scratch%px(:, k) = batch%x(:, scratch%list(k))
      end do
      call lsf%evaluate_batch(scratch%px(:, 1:nact), 1, scratch%p0(1:nact), &
                              lsf1_r=scratch%p1(:, 1:nact), error=error)
      if (allocated(error)) return
      nkeep = 0
      do k = 1, nact
         j = scratch%list(k)
         irun = batch%run_of(j)
         if (lsf%unreached(k)) then
            ! The anchor jet is read later by every gate; without it the run is lost
            if (all(batch%x(:, j) == batch%run_xyz(:, irun))) then
               write (buffer, "(i0)") batch%run_anchor(irun)
               call fatal_error(error, "Projection Newton: anchor "//trim(buffer)// &
                                " lies beyond the level set reach, no atom screens it")
               return
            end if
            ! Far shell seeds are routine; the column is dropped unread
            batch%status(j) = point_failed_unreached
            cycle
         end if
         nkeep = nkeep + 1
         scratch%list(nkeep) = j
         batch%S(j) = scratch%p0(k)
         batch%gradS(:, j) = scratch%p1(:, k)
         scratch%fnorm(j) = residual_norm(w, batch%x(:, j), batch%run_xyz(:, irun), &
                                          batch%lambda(j), batch%S(j), batch%gradS(:, j))
         if (all(batch%x(:, j) == batch%run_xyz(:, irun))) then
            call batch%set_anchor_jet(irun, batch%S(j), batch%gradS(:, j))
         end if
      end do
      nact = nkeep

      do
         !> Convergence, budget and stage decisions; stage-1 points first
         nkeep = 0
         do k = 1, nact
            j = scratch%list(k)
            if (batch%status(j) /= point_pending) cycle
            if (scratch%fnorm(j) < cfg%tol) then
               batch%status(j) = point_converged
               cycle
            end if
            if (scratch%iter(j) >= cfg%max_iter) then
               batch%status(j) = point_failed_iterations
               cycle
            end if
            if (scratch%stage(j) == 1) then
               if (scratch%fnorm(j) < cfg%switch_fnorm .or. &
                   scratch%ncheap(j) >= cfg%max_cheap_steps) scratch%stage(j) = 2
            end if
            nkeep = nkeep + 1
            scratch%list(nkeep) = j
         end do
         nact = nkeep
         if (nact == 0) exit
         call partition_by_stage(scratch%list, nact, scratch%stage, m1, scratch%ptmp)

         !> Order-2 refresh of stage-2 points that carry only a gradient
         call refresh_hessians(batch, lsf, scratch, scratch%list(m1 + 1:nact), error)
         if (allocated(error)) return

         !> Newton directions, merit penalty and descent check
         nls = 0
         do k = 1, nact
            j = scratch%list(k)
            ! Failed by the refresh
            if (batch%status(j) /= point_pending) cycle
            a = batch%run_xyz(:, batch%run_of(j))
            if (scratch%stage(j) == 2) then
               call full_step(w, batch%x(:, j), a, batch%lambda(j), batch%S(j), batch%gradS(:, j), &
                              batch%hessS(:, j), step, ok)
               if (.not. ok) then
                  batch%status(j) = point_failed_singular
                  cycle
               end if
               ! Penalty the full step asks for, stored below only with the step taken:
               ! the multiplier of a rejected step must not inflate it for good
               mu_full = max(scratch%mu(j), abs(batch%lambda(j) + step(4)) + cfg%merit_margin)
               ! No merit descent (indefinite tangent block): the cheap direction instead
               if (merit_slope(w, batch%x(:, j), a, batch%S(j), batch%gradS(:, j), mu_full, &
                               step(1:3)) >= 0.0_wp) then
                  call cheap_step(w, batch%x(:, j), a, batch%lambda(j), batch%S(j), batch%gradS(:, j), &
                                  step, ok)
                  if (.not. ok) then
                     batch%status(j) = point_failed_gradient
                     cycle
                  end if
               end if
            else
               call cheap_step(w, batch%x(:, j), a, batch%lambda(j), batch%S(j), batch%gradS(:, j), &
                               step, ok)
               if (.not. ok) then
                  batch%status(j) = point_failed_gradient
                  cycle
               end if
            end if
            scratch%mu(j) = max(scratch%mu(j), abs(batch%lambda(j) + step(4)) + cfg%merit_margin)
            ! A scaled step is no prediction of the foot
            if (cfg%capture_duplicates .and. norm2(step(1:3)) < cfg%max_step) then
               if (lands_on_accepted(batch, scratch%prev_accepted, j, batch%x(:, j) + step(1:3))) then
                  batch%status(j) = point_captured
                  cycle
               end if
            end if
            if (norm2(step(1:3)) > cfg%max_step) step = step*(cfg%max_step/norm2(step(1:3)))
            scratch%dz(:, j) = step
            scratch%alpha(j) = 1.0_wp
            scratch%merit(j) = 0.5_wp*w*sum((batch%x(:, j) - a)**2) + scratch%mu(j)*abs(batch%S(j))
            scratch%dmerit(j) = merit_slope(w, batch%x(:, j), a, batch%S(j), batch%gradS(:, j), &
                                            scratch%mu(j), step(1:3))
            nls = nls + 1
            scratch%lsl(nls) = j
         end do

         !> Batched backtracking on the merit (||F|| near the solution): one
         !> evaluation per stage per round
         do while (nls > 0)
            call partition_by_stage(scratch%lsl, nls, scratch%stage, m1, scratch%ptmp)
            do k = 1, nls
               j = scratch%lsl(k)
               scratch%px(:, k) = batch%x(:, j) + scratch%alpha(j)*scratch%dz(1:3, j)
               scratch%plam(k) = batch%lambda(j) + scratch%alpha(j)*scratch%dz(4, j)
            end do
            m = nls
            ! Two batches on one LSF: each slice's mask read before the next call
            call lsf%evaluate_batch(scratch%px(:, 1:m1), 1, scratch%p0(1:m1), &
                                    lsf1_r=scratch%p1(:, 1:m1), error=error)
            if (allocated(error)) return
            do k = 1, m1
               scratch%punreached(k) = lsf%unreached(k)
            end do
            call lsf%evaluate_batch(scratch%px(:, m1 + 1:m), 2, scratch%p0(m1 + 1:m), &
                                    lsf1_r=scratch%p1(:, m1 + 1:m), &
                                    lsf2_rr_packed=scratch%p2(:, m1 + 1:m), error=error)
            if (allocated(error)) return
            do k = m1 + 1, m
               scratch%punreached(k) = lsf%unreached(k - m1)
            end do

            nkeep = 0
            do k = 1, m
               j = scratch%lsl(k)
               a = batch%run_xyz(:, batch%run_of(j))
               ! Unreached trial: rejected unread; at the floor the column fails
               ! instead of committing an unspecified jet
               if (scratch%punreached(k)) then
                  if (scratch%alpha(j) <= cfg%ls_floor) then
                     batch%status(j) = point_failed_unreached
                  else
                     scratch%alpha(j) = max(0.5_wp*scratch%alpha(j), cfg%ls_floor)
                     nkeep = nkeep + 1
                     scratch%lsl(nkeep) = j
                  end if
                  cycle
               end if
               fnew = residual_norm(w, scratch%px(:, k), a, scratch%plam(k), scratch%p0(k), &
                                    scratch%p1(:, k))
               if (scratch%fnorm(j) <= cfg%switch_fnorm) then
                  ! Maratos guard: near the solution the residual decides
                  accept = fnew <= (1.0_wp - armijo_c*scratch%alpha(j))*scratch%fnorm(j)
               else
                  merit_new = 0.5_wp*w*sum((scratch%px(:, k) - a)**2) + scratch%mu(j)*abs(scratch%p0(k))
                  accept = merit_new <= scratch%merit(j) + armijo_c*scratch%alpha(j)*scratch%dmerit(j)
               end if
               accept = accept .or. scratch%alpha(j) <= cfg%ls_floor
               if (.not. accept) then
                  scratch%alpha(j) = max(0.5_wp*scratch%alpha(j), cfg%ls_floor)
                  nkeep = nkeep + 1
                  scratch%lsl(nkeep) = j
                  cycle
               end if

               !> Commit the trial and its jet
               dzmax = max(maxval(abs(scratch%px(:, k) - batch%x(:, j))), &
                           abs(scratch%plam(k) - batch%lambda(j)))
               dzfloor = max(0.1_wp*cfg%tol, 8.0_wp*epsilon(1.0_wp)* &
                             max(maxval(abs(batch%x(:, j))), abs(batch%lambda(j))))
               batch%x(:, j) = scratch%px(:, k)
               batch%lambda(j) = scratch%plam(k)
               batch%S(j) = scratch%p0(k)
               batch%gradS(:, j) = scratch%p1(:, k)
               if (scratch%stage(j) == 2) then
                  batch%hessS(:, j) = scratch%p2(:, k)
                  scratch%has_hess(j) = .true.
                  scratch%n2(j) = scratch%n2(j) + 1
               else
                  scratch%has_hess(j) = .false.
                  scratch%n1(j) = scratch%n1(j) + 1
                  scratch%ncheap(j) = scratch%ncheap(j) + 1
               end if
               scratch%fnorm(j) = fnew
               scratch%iter(j) = scratch%iter(j) + 1

               ! A residual below tolerance converges at the top of the next round
               if (fnew >= cfg%tol .and. dzmax <= dzfloor) batch%status(j) = point_converged
            end do
            nls = nkeep
         end do
      end do

      !> Final jet is always order 2
      nact = 0
      do j = j0, j1
         if (batch%status(j) /= point_converged) cycle
         if (scratch%has_hess(j)) cycle
         nact = nact + 1
         scratch%list(nact) = j
      end do
      call refresh_hessians(batch, lsf, scratch, scratch%list(1:nact), error)
      if (allocated(error)) return

      !> Derived per-point quantities of converged columns
      do j = j0, j1
         if (.not. scratch%processed(j)) cycle
         if (present(stage1_steps)) stage1_steps(j) = scratch%n1(j)
         if (present(stage2_steps)) stage2_steps(j) = scratch%n2(j)
         if (batch%status(j) /= point_converged) cycle
         a = batch%run_xyz(:, batch%run_of(j))
         batch%rho(j) = norm2(batch%x(:, j) - a)
         batch%phi(j) = 0.5_wp*w*sum((batch%x(:, j) - a)**2)
         batch%normal(:, j) = batch%gradS(:, j)/norm2(batch%gradS(:, j))
      end do
   end subroutine newton_solve_batch

   !> Re-evaluate the listed columns at order 2 where only a gradient is stored
   !>
   !> - an unreached column is `point_failed_unreached`, its jet left unstored
   !>
   !> @param[inout] batch     Points
   !> @param[inout] lsf       Thread-local level set
   !> @param[inout] scratch   Newton scratch
   !> @param[in]    cols      Candidate columns
   !> @param[out]   error     LSF evaluation failure
   subroutine refresh_hessians(batch, lsf, scratch, cols, error)
      type(projection_batch_type), intent(inout) :: batch
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      type(newton_scratch_type), intent(inout) :: scratch
      integer, intent(in) :: cols(:)
      type(error_type), allocatable, intent(out) :: error

      integer :: k, j, m

      m = 0
      do k = 1, size(cols)
         j = cols(k)
         if (scratch%has_hess(j)) cycle
         m = m + 1
         scratch%lsl(m) = j
         scratch%px(:, m) = batch%x(:, j)
      end do
      if (m == 0) return
      call lsf%evaluate_batch(scratch%px(:, 1:m), 2, scratch%p0(1:m), &
                              lsf1_r=scratch%p1(:, 1:m), lsf2_rr_packed=scratch%p2(:, 1:m), &
                              error=error)
      if (allocated(error)) return
      do k = 1, m
         j = scratch%lsl(k)
         ! A stored iterate was reached when committed; defensive only
         if (lsf%unreached(k)) then
            batch%status(j) = point_failed_unreached
            cycle
         end if
         batch%S(j) = scratch%p0(k)
         batch%gradS(:, j) = scratch%p1(:, k)
         batch%hessS(:, j) = scratch%p2(:, k)
         scratch%has_hess(j) = .true.
      end do
   end subroutine refresh_hessians

   !> Stable partition of `list(1:n)` into stage-1 columns, then stage-2 columns
   !>
   !> @param[inout] list  Column list
   !> @param[in]    n     Entries in use
   !> @param[in]    stage Stage per column
   !> @param[out]   m1    Number of stage-1 entries
   !> @param[inout] tmp   Workspace of at least `n` entries
   pure subroutine partition_by_stage(list, n, stage, m1, tmp)
      integer, intent(inout) :: list(:)
      integer, intent(in) :: n
      integer, intent(in) :: stage(:)
      integer, intent(out) :: m1
      integer, intent(inout) :: tmp(:)

      integer :: k, m2

      m1 = 0
      m2 = 0
      do k = 1, n
         if (stage(list(k)) == 1) then
            m1 = m1 + 1
            list(m1) = list(k)
         else
            m2 = m2 + 1
            tmp(m2) = list(k)
         end if
      end do
      list(m1 + 1:m1 + m2) = tmp(1:m2)
   end subroutine partition_by_stage

   !* ================================================================================= *!
   !*                                 Duplicate capture                                 *!
   !* ================================================================================= *!

   !> Link every column in `first:last` to the previous accepted column of its run
   !>
   !> - runs are contiguous column ranges, so one sweep in column order suffices
   !>
   !> @param[in]  batch         Points
   !> @param[in]  first         First column, the start of a run
   !> @param[in]  last          Last column
   !> @param[out] prev_accepted Previous accepted column of the same run, 0 for none
   pure subroutine link_accepted(batch, first, last, prev_accepted)
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: first
      integer, intent(in) :: last
      integer, intent(inout) :: prev_accepted(:)

      integer :: j, irun, latest

      irun = 0
      latest = 0
      do j = first, last
         if (batch%run_of(j) /= irun) then
            irun = batch%run_of(j)
            latest = 0
         end if
         prev_accepted(j) = latest
         if (batch%status(j) == point_accepted) latest = j
      end do
   end subroutine link_accepted

   !> Whether a prediction lies within `capture_distance` of an accepted column of the run before `j`
   !>
   !> @param[in] batch         Points
   !> @param[in] prev_accepted Previous accepted column of the same run, 0 for none
   !> @param[in] j             Column
   !> @param[in] x_pred        Predicted foot of column `j`
   pure function lands_on_accepted(batch, prev_accepted, j, x_pred) result(landed)
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: prev_accepted(:)
      integer, intent(in) :: j
      real(wp), intent(in) :: x_pred(3)
      logical :: landed

      integer :: c

      landed = .false.
      c = prev_accepted(j)
      do while (c > 0)
         if (norm2(x_pred - batch%x(:, c)) < capture_distance) then
            landed = .true.
            return
         end if
         c = prev_accepted(c)
      end do
   end function lands_on_accepted

   !* ================================================================================= *!
   !*                              Per-point linear algebra                             *!
   !* ================================================================================= *!

   !> KKT residual norm `||F||`
   !>
   !> @param[in] w      Objective weight
   !> @param[in] x      Point
   !> @param[in] a      Anchor
   !> @param[in] lambda Multiplier
   !> @param[in] S      Level set value
   !> @param[in] g      Level set gradient
   pure function residual_norm(w, x, a, lambda, S, g) result(fnorm)
      real(wp), intent(in) :: w, x(3), a(3), lambda, S, g(3)
      real(wp) :: fnorm

      fnorm = sqrt(sum((w*(x - a) - lambda*g)**2) + S*S)
   end function residual_norm

   !> Directional derivative of the merit `w/2 ||x - a||^2 + mu |S|` along `dx`
   !>
   !> - on the surface `|S|` has the one-sided slope `|grad S . dx|`
   !>
   !> @param[in] w  Objective weight
   !> @param[in] x  Point
   !> @param[in] a  Anchor
   !> @param[in] S  Level set value
   !> @param[in] g  Level set gradient
   !> @param[in] mu Penalty
   !> @param[in] dx Spatial step
   pure function merit_slope(w, x, a, S, g, mu, dx) result(slope)
      real(wp), intent(in) :: w, x(3), a(3), S, g(3), mu, dx(3)
      real(wp) :: slope

      slope = w*dot_product(x - a, dx)
      if (S > 0.0_wp) then
         slope = slope + mu*dot_product(g, dx)
      else if (S < 0.0_wp) then
         slope = slope - mu*dot_product(g, dx)
      else
         slope = slope + mu*abs(dot_product(g, dx))
      end if
   end function merit_slope

   !> Stage-1 step: `J` without `lambda Hess S`, in closed form
   !>
   !>     dlambda = (g.r - w S) / g.g,   dx = (g dlambda - r) / w,   r = w (x - a) - lambda g
   !>
   !> @param[in]  w      Objective weight
   !> @param[in]  x      Point
   !> @param[in]  a      Anchor
   !> @param[in]  lambda Multiplier
   !> @param[in]  S      Level set value
   !> @param[in]  g      Level set gradient
   !> @param[out] dz     Step (dx, dlambda)
   !> @param[out] ok     `.false.` on a vanishing gradient
   pure subroutine cheap_step(w, x, a, lambda, S, g, dz, ok)
      real(wp), intent(in) :: w, x(3), a(3), lambda, S, g(3)
      real(wp), intent(out) :: dz(4)
      logical, intent(out) :: ok

      real(wp) :: r(3), g2, dlam

      g2 = sum(g*g)
      ok = g2 > grad2_floor
      if (.not. ok) then
         dz = 0.0_wp
         return
      end if
      r = w*(x - a) - lambda*g
      dlam = (sum(g*r) - w*S)/g2
      dz(1:3) = (g*dlam - r)/w
      dz(4) = dlam
   end subroutine cheap_step

   !> Stage-2 step: the full KKT system
   !>
   !> @param[in]  w      Objective weight
   !> @param[in]  x      Point
   !> @param[in]  a      Anchor
   !> @param[in]  lambda Multiplier
   !> @param[in]  S      Level set value
   !> @param[in]  g      Level set gradient
   !> @param[in]  H      Level set Hessian, packed (6)
   !> @param[out] dz     Step (dx, dlambda)
   !> @param[out] ok     `.false.` on a singular system
   pure subroutine full_step(w, x, a, lambda, S, g, H, dz, ok)
      real(wp), intent(in) :: w, x(3), a(3), lambda, S, g(3), H(:)
      real(wp), intent(out) :: dz(4)
      logical, intent(out) :: ok

      real(wp) :: jac(4, 4), rhs(4)
      integer :: i, k

      do k = 1, 3
         do i = 1, 3
            jac(i, k) = -lambda*H(packed_index(i, k))
         end do
      end do
      do i = 1, 3
         jac(i, i) = jac(i, i) + w
      end do
      jac(1:3, 4) = -g
      jac(4, 1:3) = -g
      jac(4, 4) = 0.0_wp
      rhs(1:3) = -(w*(x - a) - lambda*g)
      rhs(4) = S
      call newton_solve_4x4(jac, rhs, dz, ok)
   end subroutine full_step

   !> Solve a 4x4 system by Gaussian elimination with partial pivoting
   !>
   !> - singular when a pivot falls below `16 eps` times the largest entry
   !>
   !> @param[in]  a  Matrix
   !> @param[in]  b  Right-hand side
   !> @param[out] x  Solution
   !> @param[out] ok `.false.` on a singular matrix
   pure subroutine newton_solve_4x4(a, b, x, ok)
      real(wp), intent(in) :: a(4, 4), b(4)
      real(wp), intent(out) :: x(4)
      logical, intent(out) :: ok

      real(wp) :: m(4, 4), r(4), tmp(4), t, f, scale
      integer :: k, i, p

      m = a
      r = b
      x = 0.0_wp
      scale = maxval(abs(m))
      ok = scale > 0.0_wp
      if (.not. ok) return

      do k = 1, 3
         p = k
         do i = k + 1, 4
            if (abs(m(i, k)) > abs(m(p, k))) p = i
         end do
         if (abs(m(p, k)) <= 16.0_wp*epsilon(1.0_wp)*scale) then
            ok = .false.
            return
         end if
         if (p /= k) then
            tmp = m(k, :)
            m(k, :) = m(p, :)
            m(p, :) = tmp
            t = r(k)
            r(k) = r(p)
            r(p) = t
         end if
         do i = k + 1, 4
            f = m(i, k)/m(k, k)
            m(i, k:4) = m(i, k:4) - f*m(k, k:4)
            r(i) = r(i) - f*r(k)
         end do
      end do
      if (abs(m(4, 4)) <= 16.0_wp*epsilon(1.0_wp)*scale) then
         ok = .false.
         return
      end if

      x(4) = r(4)/m(4, 4)
      x(3) = (r(3) - m(3, 4)*x(4))/m(3, 3)
      x(2) = (r(2) - m(2, 3)*x(3) - m(2, 4)*x(4))/m(2, 2)
      x(1) = (r(1) - m(1, 2)*x(2) - m(1, 3)*x(3) - m(1, 4)*x(4))/m(1, 1)
   end subroutine newton_solve_4x4

   !* ================================================================================= *!
   !*                                      Scratch                                      *!
   !* ================================================================================= *!

   !> Grow the Newton scratch to `n` columns
   !>
   !> @param[inout] sc Scratch
   !> @param[in]    n  Columns to hold
   subroutine reserve_scratch(sc, n)
      type(newton_scratch_type), intent(inout) :: sc
      integer, intent(in) :: n

      integer :: cap

      if (n <= sc%capacity) return
      cap = max(16, 2*sc%capacity, n)
      if (allocated(sc%stage)) then
         deallocate (sc%stage, sc%ncheap, sc%iter, sc%n1, sc%n2, sc%has_hess, &
                     sc%processed, sc%fnorm, sc%dz, sc%alpha, sc%mu, sc%merit, sc%dmerit, &
                     sc%list, sc%lsl, sc%ptmp, sc%px, sc%plam, sc%p0, sc%p1, sc%p2, &
                     sc%punreached, sc%prev_accepted)
      end if
      allocate (sc%stage(cap), sc%ncheap(cap), sc%iter(cap), sc%n1(cap), sc%n2(cap))
      allocate (sc%has_hess(cap), sc%processed(cap), sc%fnorm(cap))
      allocate (sc%dz(4, cap), sc%alpha(cap), sc%mu(cap), sc%merit(cap), sc%dmerit(cap))
      allocate (sc%list(cap), sc%lsl(cap), sc%ptmp(cap))
      allocate (sc%px(3, cap), sc%plam(cap), sc%p0(cap), sc%p1(3, cap), sc%p2(6, cap))
      allocate (sc%punreached(cap), sc%prev_accepted(cap))
      sc%capacity = cap
   end subroutine reserve_scratch

end module moist_cavity_drop_projection_newton
