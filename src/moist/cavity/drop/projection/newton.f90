!> Batched KKT Newton for the DROP closest-point projection
!>
!> Solves for every pending column of a [[projection_batch_type]], `w = phi_alpha`:
!>
!>     F(x, lambda) = [ w (x - a) - lambda grad S(x) ] = 0
!>                    [            -S(x)             ]
!>
!> - full 4x4 KKT step on order-2 jets; closed-form flat step at `lambda = 0`
!>   or on singular tangent blocks and, above `residual_test_fnorm`, on non-descent
!>   ones; spatial step capped at `max_step`
!> - backtracking on merit `w/2 ||x - a||^2 + mu |S|`, switching to `||F||`
!>   below `residual_test_fnorm`
!> - outcomes per column: `point_converged` (order-2 jet), `point_captured`
!>   (duplicate prediction near a run captor), or `point_failed` (`max_iter`,
!>   vanishing gradient, or unreached)
module moist_cavity_drop_projection_newton
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric
   use moist_math_packed, only: packed_index, packed_contract2
   use moist_utils_format, only: to_string
   use moist_cavity_drop_projection_parameters, only: drop_solve_parameters_type, drop_polish_parameters_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, &
                                                point_captured, point_failed
   implicit none(type, external)
   private

   public :: newton_scratch_type, newton_solve_batch, captor_list_type
   public :: kkt_residuals, reduced_hessian

   !> Armijo coefficient of the backtracking on `||F||`
   real(wp), parameter :: armijo_c = 1.0e-4_wp
   !> Smallest line-search step, accepted unconditionally
   real(wp), parameter :: ls_floor = 0.1_wp
   !> `||F||` at and below which the line search tests the residual instead of the merit
   real(wp), parameter :: residual_test_fnorm = 1.0e-3_wp
   !> Squared gradient norm below which the constraint is degenerate
   real(wp), parameter :: grad2_floor = 1.0e-28_wp
   !> Distance (Bohr) within which a Newton prediction lands on an accepted point
   !>
   !> - far above `branch%sep_cut` (at least `1e-6`), within which the finalize
   !>   merges two points, so a duplicate start stops rounds before it converges
   !> - well below the separation of distinct points on near-degenerate arcs,
   !>   which can be a few `1e-3`
   real(wp), parameter :: capture_distance = 1.0e-4_wp
   !> Factor on the polish's `saddle_eig_tol` a captor's smaller reduced-Hessian
   !> eigenvalue must clear: a strict minimum, never one the third-order test decides
   real(wp), parameter :: captor_eig_margin = 10.0_wp

   !> Floor of the reduced-Hessian eigenvalue threshold
   real(wp), parameter :: tol_eig_base = 1.0e-10_wp
   !> Diagonal and off-diagonal slots of a packed 3x3 Hessian
   integer, parameter :: diag_slots(3) = [1, 3, 6], offdiag_slots(3) = [2, 4, 5]

   !> Captors of the runs of a batch: one list per run, linked through the columns
   !>
   !> - a list keeps one captor per point within `capture_distance`, so it stays
   !>   O(distinct points)
   type :: captor_list_type
      !> First captor column of each run, 0 for none (run-indexed)
      integer, allocatable :: head(:)
      !> Next captor column of the same run, 0 for none (column-indexed)
      integer, allocatable :: next(:)
   contains
      !> Size the lists and fill them with the accepted columns of a column range
      procedure :: init => captors_init
      !> Register a column as a captor of its run
      procedure :: add => captors_add
      !> Whether a point lands on a captor of a run
      procedure :: lands => captors_lands
   end type captor_list_type

   !> Per-thread Newton scratch, column- and list-indexed, grown geometrically
   type :: newton_scratch_type
      !> Allocated columns
      integer :: capacity = 0
      !> Steps taken
      integer, allocatable :: iter(:)
      !> Whether the column's stored jet includes the gradient and Hessian
      logical, allocatable :: has_grad(:), has_hess(:)
      !> Column was pending when the solve started
      logical, allocatable :: processed(:)
      !> Current residual norm
      real(wp), allocatable :: fnorm(:)
      !> Newton direction
      real(wp), allocatable :: dz(:, :)
      !> Merit penalty, merit at the iterate, merit directional derivative
      real(wp), allocatable :: mu(:), merit(:), dmerit(:)
      !> Active list, line-search list, refresh and trial list
      integer, allocatable :: list(:), lsl(:), ptmp(:)
      !> Gathered trial points, multipliers and jets (list-indexed), packed Hessian `p2(6, :)`
      real(wp), allocatable :: px(:, :), plam(:), p0(:), p1(:, :), p2(:, :)
      !> Captors of the runs, with `capture_duplicates`; the polish reads them too
      type(captor_list_type) :: captors
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
   !>   `hessS`, `phi`, `pdist` and `normal` of its final iterate, `point_captured`
   !>   (with `capture_duplicates`), or `point_failed`
   !>
   !> @param[inout] batch     Points to solve
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    cfg       Newton settings
   !> @param[in]    pcfg      Polish settings; its KKT and `saddle_eig_tol` tolerances
   !>                         set the captor test (read with `capture_duplicates` only)
   !> @param[in]    phi_alpha Anchor objective weight `w`
   !> @param[inout] scratch   Per-thread scratch
   !> @param[out]   error     LSF evaluation failure or unreached anchor (build-fatal)
   !> @param[in]    first     First column (default 1)
   !> @param[in]    last      Last column (default `batch%n`)
   !> @param[inout] steps     Steps taken per column (diagnostic)
   subroutine newton_solve_batch(batch, lsf, cfg, pcfg, phi_alpha, scratch, error, first, last, steps)
      type(projection_batch_type), intent(inout) :: batch
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      type(drop_solve_parameters_type), intent(in) :: cfg
      type(drop_polish_parameters_type), intent(in) :: pcfg
      real(wp), intent(in) :: phi_alpha
      type(newton_scratch_type), intent(inout) :: scratch
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: first
      integer, intent(in), optional :: last
      integer, intent(inout), optional :: steps(:)

      integer :: j0, j1, j, k, nact, nkeep, nls, m, m0, ls_round, irun
      real(wp) :: w, a(3), step(4), mu_full, g_len, alpha, captor_tol
      logical :: ok, predicted

      j0 = 1
      j1 = batch%n
      if (present(first)) j0 = first
      if (present(last)) j1 = last
      if (j1 < j0) return
      call reserve_scratch(scratch, j1)
      w = phi_alpha
      captor_tol = captor_eig_margin*pcfg%saddle_eig_tol

      !> Start: one order-1 evaluation of every pending column
      nact = 0
      do j = j0, j1
         scratch%processed(j) = batch%status(j) == point_pending
         if (.not. scratch%processed(j)) cycle
         scratch%iter(j) = 0
         scratch%has_grad(j) = .true.
         scratch%has_hess(j) = .false.
         scratch%mu(j) = max(1.0_wp, abs(batch%lambda(j)))
         nact = nact + 1
         scratch%list(nact) = j
      end do
      if (nact == 0) return
      if (cfg%capture_duplicates) then
         call scratch%captors%init(batch, batch%run_start(batch%run_of(j0)), j1)
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
               call fatal_error(error, "Projection Newton: anchor "//to_string(batch%run_id(irun))// &
                                " lies beyond the level set reach, no atom screens it")
               return
            end if
            ! Far shell seeds are routine; the column is dropped unread
            batch%status(j) = point_failed
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
         !> Convergence and budget decisions
         nkeep = 0
         do k = 1, nact
            j = scratch%list(k)
            if (batch%status(j) /= point_pending) cycle
            if (scratch%fnorm(j) < cfg%tol) then
               call converge(j)
               cycle
            end if
            if (scratch%iter(j) >= cfg%max_iter) then
               batch%status(j) = point_failed
               cycle
            end if
            nkeep = nkeep + 1
            scratch%list(nkeep) = j
         end do
         nact = nkeep
         if (nact == 0) exit

         !> Order-2 refresh of points that carry only a gradient and need the
         !> Hessian: at `lambda = 0` the system holds none
         m = 0
         do k = 1, nact
            j = scratch%list(k)
            if (scratch%has_hess(j) .or. batch%lambda(j) == 0.0_wp) cycle
            m = m + 1
            scratch%ptmp(m) = j
         end do
         call refresh_hessians(batch, lsf, scratch, scratch%ptmp(1:m), error)
         if (allocated(error)) return

         !> Newton directions, merit penalty and descent check
         nls = 0
         do k = 1, nact
            j = scratch%list(k)
            ! Failed by the refresh
            if (batch%status(j) /= point_pending) cycle
            a = batch%run_xyz(:, batch%run_of(j))
            ok = .false.
            if (scratch%has_hess(j)) then
               call full_step(w, batch%x(:, j), a, batch%lambda(j), batch%S(j), batch%gradS(:, j), &
                              batch%hessS(:, j), step, ok)
               if (ok) then
                  ! Merit descent only away from the solution: near it the residual line
                  ! search decides, and a full step toward a saddle is no merit descent
                  if (.not. (scratch%fnorm(j) <= residual_test_fnorm)) then
                     ! Penalty the full step asks for, stored below only with the step taken:
                     ! the multiplier of a rejected step must not inflate it for good
                     mu_full = max(scratch%mu(j), abs(batch%lambda(j) + step(4)) + cfg%merit_margin)
                     ok = merit_slope(w, batch%x(:, j), a, batch%S(j), batch%gradS(:, j), mu_full, &
                                      step(1:3)) < 0.0_wp
                  end if
               else if (.not. singular_tangent_block(w, batch%lambda(j), batch%gradS(:, j), batch%hessS(:, j))) then
                  ! Singular with a regular tangent block: the constraint gradient vanishes against `lambda Hess S`
                  batch%status(j) = point_failed
                  cycle
               end if
            end if
            ! `lambda = 0`: the full system itself, a prediction; a singular tangent
            ! block or no merit descent away from the solution: the fallback
            ! without the constraint curvature, no prediction
            predicted = ok .or. .not. scratch%has_hess(j)
            if (.not. ok) then
               call flat_step(w, batch%x(:, j), a, batch%lambda(j), batch%S(j), batch%gradS(:, j), &
                              step, ok)
               if (.not. ok) then
                  batch%status(j) = point_failed
                  cycle
               end if
            end if
            scratch%mu(j) = max(scratch%mu(j), abs(batch%lambda(j) + step(4)) + cfg%merit_margin)
            ! A scaled step is no prediction of the point either
            if (cfg%capture_duplicates .and. predicted .and. norm2(step(1:3)) < cfg%max_step) then
               if (scratch%captors%lands(batch, batch%run_of(j), batch%x(:, j) + step(1:3))) then
                  batch%status(j) = point_captured
                  cycle
               end if
            end if
            if (norm2(step(1:3)) > cfg%max_step) step = step*(cfg%max_step/norm2(step(1:3)))
            scratch%dz(:, j) = step
            scratch%merit(j) = 0.5_wp*w*sum((batch%x(:, j) - a)**2) + scratch%mu(j)*abs(batch%S(j))
            scratch%dmerit(j) = merit_slope(w, batch%x(:, j), a, batch%S(j), batch%gradS(:, j), &
                                            scratch%mu(j), step(1:3))
            nls = nls + 1
            scratch%lsl(nls) = j
         end do

         !> Batched backtracking on the merit (||F|| near the solution), one step
         !> length per round for every column still searching: round 1 (alpha = 1)
         !> and the floor round (unconditional accept) evaluate order 2 directly;
         !> intermediate rounds evaluate order 0 for the merit test and order 1
         !> only where the Maratos guard reads ||F||
         alpha = 1.0_wp
         ls_round = 0
         do while (nls > 0)
            ls_round = ls_round + 1
            nkeep = 0
            if (ls_round == 1 .or. alpha <= ls_floor) then
               scratch%ptmp(1:nls) = scratch%lsl(1:nls)
               call try_steps(2, 1, nls)
            else
               ! Merit columns first, then residual columns, column order kept in both;
               ! the split [[try_steps]] makes, total: a NaN `||F||` is a merit column
               m0 = 0
               do k = 1, nls
                  j = scratch%lsl(k)
                  if (.not. (scratch%fnorm(j) <= residual_test_fnorm)) then
                     m0 = m0 + 1
                     scratch%ptmp(m0) = j
                  end if
               end do
               m = m0
               do k = 1, nls
                  j = scratch%lsl(k)
                  if (scratch%fnorm(j) <= residual_test_fnorm) then
                     m = m + 1
                     scratch%ptmp(m) = j
                  end if
               end do
               call try_steps(0, 1, m0)
               if (allocated(error)) return
               call try_steps(1, m0 + 1, nls)
            end if
            if (allocated(error)) return
            nls = nkeep
            alpha = max(0.5_wp*alpha, ls_floor)
         end do

         !> Refresh pending columns accepted at order 0 so fnorm and the order-2
         !> jet are ready before the next convergence test
         m = 0
         do k = 1, nact
            j = scratch%list(k)
            if (batch%status(j) /= point_pending) cycle
            if (scratch%has_grad(j)) cycle
            m = m + 1
            scratch%ptmp(m) = j
         end do
         if (m > 0) then
            call refresh_hessians(batch, lsf, scratch, scratch%ptmp(1:m), error)
            if (allocated(error)) return
            do k = 1, m
               j = scratch%ptmp(k)
               if (batch%status(j) /= point_pending) cycle
               a = batch%run_xyz(:, batch%run_of(j))
               scratch%fnorm(j) = residual_norm(w, batch%x(:, j), a, batch%lambda(j), &
                                                batch%S(j), batch%gradS(:, j))
            end do
         end if
      end do

      !> Final jet is always order 2
      nact = 0
      do j = j0, j1
         if (.not. scratch%processed(j)) cycle
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
         if (present(steps)) steps(j) = scratch%iter(j)
         if (batch%status(j) /= point_converged) cycle
         a = batch%run_xyz(:, batch%run_of(j))
         batch%pdist(j) = norm2(batch%x(:, j) - a)
         batch%phi(j) = 0.5_wp*w*sum((batch%x(:, j) - a)**2)
         ! Zero gradient: anchor on a critical point of S, no normal; f_crit closes the point
         g_len = norm2(batch%gradS(:, j))
         if (g_len > 0.0_wp) then
            batch%normal(:, j) = batch%gradS(:, j)/g_len
         else
            batch%normal(:, j) = 0.0_wp
         end if
      end do

   contains

      !> Mark a column converged; with `capture_duplicates`, an order-2 strict
      !> minimum becomes a captor of its run
      !>
      !> @param[in] j Column
      subroutine converge(j)
         integer, intent(in) :: j

         batch%status(j) = point_converged
         if (.not. (cfg%capture_duplicates .and. scratch%has_hess(j))) return
         if (is_strict_minimum(w, batch%x(:, j), batch%run_xyz(:, batch%run_of(j)), batch%lambda(j), &
                               batch%S(j), batch%gradS(:, j), batch%hessS(:, j), pcfg, captor_tol)) then
            call scratch%captors%add(batch, j)
         end if
      end subroutine converge

      !> One line-search round of the columns `ptmp(k0:k1)` at step `alpha`
      !>
      !> - order 2: merit or residual test as the column's `||F||` asks, and the
      !>   unconditional accept at the floor; order 0: merit test; order 1: residual test
      !> - an accepted trial is committed with its jet, a rejected one is kept in `lsl`
      !>
      !> @param[in] order Jet order of the trials (0, 1, 2)
      !> @param[in] k0    First trial in `ptmp`
      !> @param[in] k1    Last trial in `ptmp`
      subroutine try_steps(order, k0, k1)
         integer, intent(in) :: order
         integer, intent(in) :: k0
         integer, intent(in) :: k1

         integer :: i, j, n
         real(wp) :: a(3), fnew, merit_new
         logical :: floor, accept

         n = k1 - k0 + 1
         if (n < 1) return
         do i = 1, n
            j = scratch%ptmp(k0 + i - 1)
            scratch%px(:, i) = batch%x(:, j) + alpha*scratch%dz(1:3, j)
            scratch%plam(i) = batch%lambda(j) + alpha*scratch%dz(4, j)
         end do
         select case (order)
         case (0)
            call lsf%evaluate_batch(scratch%px(:, 1:n), 0, scratch%p0(1:n), error=error)
         case (1)
            call lsf%evaluate_batch(scratch%px(:, 1:n), 1, scratch%p0(1:n), &
                                    lsf1_r=scratch%p1(:, 1:n), error=error)
         case default
            call lsf%evaluate_batch(scratch%px(:, 1:n), 2, scratch%p0(1:n), &
                                    lsf1_r=scratch%p1(:, 1:n), &
                                    lsf2_rr_packed=scratch%p2(:, 1:n), error=error)
         end select
         if (allocated(error)) return

         floor = alpha <= ls_floor
         fnew = 0.0_wp
         do i = 1, n
            j = scratch%ptmp(k0 + i - 1)
            ! Unreached trial: rejected unread; at the floor the column fails
            ! instead of committing an unspecified jet
            if (lsf%unreached(i)) then
               if (floor) then
                  batch%status(j) = point_failed
               else
                  nkeep = nkeep + 1
                  scratch%lsl(nkeep) = j
               end if
               cycle
            end if
            a = batch%run_xyz(:, batch%run_of(j))
            if (order >= 1) then
               fnew = residual_norm(w, scratch%px(:, i), a, scratch%plam(i), scratch%p0(i), &
                                    scratch%p1(:, i))
            end if
            if (scratch%fnorm(j) <= residual_test_fnorm) then
               ! Maratos guard: near the solution the residual decides
               accept = fnew <= (1.0_wp - armijo_c*alpha)*scratch%fnorm(j)
            else
               merit_new = 0.5_wp*w*sum((scratch%px(:, i) - a)**2) + scratch%mu(j)*abs(scratch%p0(i))
               accept = merit_new <= scratch%merit(j) + armijo_c*alpha*scratch%dmerit(j)
            end if
            if (.not. (accept .or. floor)) then
               nkeep = nkeep + 1
               scratch%lsl(nkeep) = j
               cycle
            end if
            call commit(j, i, order, fnew)
         end do
      end subroutine try_steps

      !> Commit trial `i` of column `j` with its jet of order `order`
      !>
      !> - a step below the round-off floor converges the column; at order 1 and
      !>   2 a residual below tolerance converges at the top of the next round
      !>
      !> @param[in] j     Column
      !> @param[in] i     Trial slot in `px`, `plam` and the jets
      !> @param[in] order Jet order of the trial
      !> @param[in] fnew  Residual norm of the trial, read at order 1 and 2
      subroutine commit(j, i, order, fnew)
         integer, intent(in) :: j
         integer, intent(in) :: i
         integer, intent(in) :: order
         real(wp), intent(in) :: fnew

         real(wp) :: dzmax, dzfloor

         dzmax = max(maxval(abs(scratch%px(:, i) - batch%x(:, j))), &
                     abs(scratch%plam(i) - batch%lambda(j)))
         dzfloor = max(0.1_wp*cfg%tol, 8.0_wp*epsilon(1.0_wp)* &
                       max(maxval(abs(batch%x(:, j))), abs(batch%lambda(j))))
         batch%x(:, j) = scratch%px(:, i)
         batch%lambda(j) = scratch%plam(i)
         batch%S(j) = scratch%p0(i)
         scratch%has_grad(j) = order >= 1
         scratch%has_hess(j) = order >= 2
         if (order >= 1) then
            batch%gradS(:, j) = scratch%p1(:, i)
            scratch%fnorm(j) = fnew
         end if
         if (order >= 2) batch%hessS(:, j) = scratch%p2(:, i)
         scratch%iter(j) = scratch%iter(j) + 1

         if (dzmax > dzfloor) return
         if (order == 0) then
            batch%status(j) = point_converged
         else if (fnew >= cfg%tol) then
            call converge(j)
         end if
      end subroutine commit

   end subroutine newton_solve_batch

   !> Re-evaluate the listed columns at order 2 where only a gradient is stored
   !>
   !> - an unreached column is `point_failed`, its jet left unstored
   !> - `cols` may be a section of `scratch%ptmp` or `scratch%list`, which are
   !>   therefore never written here
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
            batch%status(j) = point_failed
            cycle
         end if
         batch%S(j) = scratch%p0(k)
         batch%gradS(:, j) = scratch%p1(:, k)
         batch%hessS(:, j) = scratch%p2(:, k)
         scratch%has_grad(j) = .true.
         scratch%has_hess(j) = .true.
      end do
   end subroutine refresh_hessians

   !* ================================================================================= *!
   !*                                 Duplicate capture                                 *!
   !* ================================================================================= *!

   !> Size the lists for `batch` and fill them with the `point_accepted` columns
   !> of `first:last`
   !>
   !> - runs are contiguous column ranges, so one sweep in column order suffices
   !>
   !> @param[inout] self  Captor lists
   !> @param[in]    batch Points
   !> @param[in]    first First column, the start of a run
   !> @param[in]    last  Last column
   subroutine captors_init(self, batch, first, last)
      class(captor_list_type), intent(inout) :: self
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: first
      integer, intent(in) :: last

      integer :: j, irun

      call grow_links(self%head, batch%nrun)
      call grow_links(self%next, last)
      irun = 0
      do j = first, last
         if (batch%run_of(j) /= irun) then
            irun = batch%run_of(j)
            self%head(irun) = 0
         end if
         self%next(j) = 0
         if (batch%status(j) == point_accepted) call self%add(batch, j)
      end do
   end subroutine captors_init

   !> Register column `j` as a captor of its run
   !>
   !> - skips a column within `capture_distance` of a captor already listed
   !>
   !> @param[inout] self  Captor lists
   !> @param[in]    batch Points
   !> @param[in]    j     Accepted column, or a strict minimum the polish will accept
   pure subroutine captors_add(self, batch, j)
      class(captor_list_type), intent(inout) :: self
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: j

      integer :: irun

      irun = batch%run_of(j)
      if (self%lands(batch, irun, batch%x(:, j))) return
      self%next(j) = self%head(irun)
      self%head(irun) = j
   end subroutine captors_add

   !> Whether `x_pred` lies within `capture_distance` of a captor of run `irun`
   !>
   !> @param[in] self   Captor lists
   !> @param[in] batch  Points
   !> @param[in] irun   Run index
   !> @param[in] x_pred Predicted point
   pure function captors_lands(self, batch, irun, x_pred) result(landed)
      class(captor_list_type), intent(in) :: self
      type(projection_batch_type), intent(in) :: batch
      integer, intent(in) :: irun
      real(wp), intent(in) :: x_pred(3)
      logical :: landed

      integer :: c

      landed = .false.
      c = self%head(irun)
      do while (c > 0)
         if (norm2(x_pred - batch%x(:, c)) < capture_distance) then
            landed = .true.
            return
         end if
         c = self%next(c)
      end do
   end function captors_lands

   !> Grow a link array geometrically to at least `n` entries; contents not kept
   !>
   !> @param[inout] links Link array
   !> @param[in]    n     Entries to hold
   subroutine grow_links(links, n)
      integer, allocatable, intent(inout) :: links(:)
      integer, intent(in) :: n

      integer :: cap

      if (allocated(links)) then
         if (size(links) >= n) return
         cap = max(n, 2*size(links))
         deallocate (links)
      else
         cap = max(16, n)
      end if
      allocate (links(cap))
   end subroutine grow_links

   !> Whether a converged order-2 column is a strict constrained local minimum
   !>
   !> - uses the KKT tolerances of `pcfg` and the reduced Hessian of
   !>   [[check_constrained_optimality]], through the same [[kkt_residuals]] and
   !>   [[reduced_hessian]]; with `rel_eig_tol` at least `pcfg%saddle_eig_tol`,
   !>   every column passing here is accepted by the check
   !> - saddles (maxima included), flat and degenerate points return `.false.`
   !>
   !> @param[in] w           Objective weight
   !> @param[in] x           Converged point
   !> @param[in] a           Anchor
   !> @param[in] lambda      Multiplier
   !> @param[in] S           Level set value
   !> @param[in] g           Level set gradient
   !> @param[in] H           Level set Hessian, packed (6)
   !> @param[in] pcfg        Polish settings (KKT tolerances)
   !> @param[in] rel_eig_tol Reduced-Hessian eigenvalue threshold relative to `||H_L||_F`
   pure function is_strict_minimum(w, x, a, lambda, S, g, H, pcfg, rel_eig_tol) result(strict_min)
      real(wp), intent(in) :: w, x(3), a(3), lambda, S, g(3), H(6), rel_eig_tol
      type(drop_polish_parameters_type), intent(in) :: pcfg
      logical :: strict_min

      real(wp) :: res_feas, res_stat, t1(3), t2(3), mu_min, mu_max, v_min(2), v_max(2), eig_tol

      strict_min = .false.
      call kkt_residuals(w, x, a, lambda, S, g, res_feas, res_stat)
      if (res_feas > pcfg%kkt_tol_feas) return
      if (norm2(g) < pcfg%kkt_tol_grad) return
      if (res_stat > pcfg%kkt_tol_stat) return
      call reduced_hessian(w, lambda, g, H, rel_eig_tol, t1, t2, mu_min, mu_max, v_min, v_max, eig_tol)
      strict_min = mu_min > eig_tol
   end function is_strict_minimum

   !* ================================================================================= *!
   !*                          Optimality of a converged point                          *!
   !* ================================================================================= *!

   !> KKT residuals `|S|` and `||w (x - a) - lambda grad S||`
   !>
   !> - one routine for the captor test and the polish check, so both read the
   !>   same values
   !>
   !> @param[in]  w        Objective weight
   !> @param[in]  x        Point
   !> @param[in]  a        Anchor
   !> @param[in]  lambda   Multiplier
   !> @param[in]  S        Level set value
   !> @param[in]  g        Level set gradient
   !> @param[out] res_feas Feasibility residual
   !> @param[out] res_stat Stationarity residual
   pure subroutine kkt_residuals(w, x, a, lambda, S, g, res_feas, res_stat)
      real(wp), intent(in) :: w, x(3), a(3), lambda, S, g(3)
      real(wp), intent(out) :: res_feas, res_stat

      res_feas = abs(S)
      res_stat = norm2(w*(x - a) - lambda*g)
   end subroutine kkt_residuals

   !> Reduced Hessian `T^T (w I - lambda Hess S) T` on the tangent plane of `g`,
   !> its eigen pairs and the saddle threshold
   !>
   !>     eig_tol = max(tol_eig_base, rel_tol ||H_L||_F)
   !>
   !> - one routine for the captor test, the polish check and the escape walk,
   !>   so all three agree on what a saddle and a minimum are
   !> - Frobenius norm of the packed Hessian, off-diagonal slots twice
   !> - `g` must not vanish
   !>
   !> @param[in]  w       Objective weight
   !> @param[in]  lambda  Multiplier
   !> @param[in]  g       Level set gradient
   !> @param[in]  H       Level set Hessian, packed (6)
   !> @param[in]  rel_tol Threshold relative to `||H_L||_F`
   !> @param[out] t1      First tangent
   !> @param[out] t2      Second tangent
   !> @param[out] mu_min  Smaller eigenvalue
   !> @param[out] mu_max  Larger eigenvalue
   !> @param[out] v_min   Eigenvector of `mu_min` in the tangent frame
   !> @param[out] v_max   Eigenvector of `mu_max` in the tangent frame
   !> @param[out] eig_tol Saddle threshold: below `-eig_tol` a saddle, above `eig_tol` a minimum
   pure subroutine reduced_hessian(w, lambda, g, H, rel_tol, t1, t2, mu_min, mu_max, v_min, v_max, eig_tol)
      real(wp), intent(in) :: w, lambda, g(3), H(6), rel_tol
      real(wp), intent(out) :: t1(3), t2(3), mu_min, mu_max, v_min(2), v_max(2), eig_tol

      real(wp) :: H_L(6), Ht1(3), Ht2(3)

      call setup_tangent_frame(g/norm2(g), t1, t2)
      H_L = -lambda*H
      H_L(diag_slots) = w - lambda*H(diag_slots)
      call packed_contract2(H_L, t1, Ht1, 3)
      call packed_contract2(H_L, t2, Ht2, 3)
      call eig_2x2_symmetric(dot_product(t1, Ht1), dot_product(t1, Ht2), dot_product(t2, Ht2), &
                             mu_min, mu_max, v_min, v_max)
      eig_tol = max(tol_eig_base, rel_tol*sqrt(sum(H_L(diag_slots)**2) + 2.0_wp*sum(H_L(offdiag_slots)**2)))
   end subroutine reduced_hessian

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

   !> Step of `J` without the constraint curvature `lambda Hess S`, in closed form
   !>
   !> - the full step at `lambda = 0`
   !> - a merit descent for `mu > |lambda + dlambda|` at any `lambda`
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
   pure subroutine flat_step(w, x, a, lambda, S, g, dz, ok)
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
   end subroutine flat_step

   !> Whether the tangent block `T^T (w I - lambda Hess S) T` of `J` is singular
   !>
   !> - `det J = -||g||^2 det(tangent block)`: a singular full system has a
   !>   singular tangent block or a gradient that vanishes against the block
   !> - singular: an eigenvalue within `16 eps ||H_L||_F` of zero, floored at
   !>   `tol_eig_base`
   !> - a gradient below the flat step's floor has no tangent plane: `.false.`
   !>
   !> @param[in] w      Objective weight
   !> @param[in] lambda Multiplier
   !> @param[in] g      Level set gradient
   !> @param[in] H      Level set Hessian, packed (6)
   pure function singular_tangent_block(w, lambda, g, H) result(singular)
      real(wp), intent(in) :: w, lambda, g(3), H(6)
      logical :: singular

      real(wp) :: t1(3), t2(3), mu_min, mu_max, v_min(2), v_max(2), eig_tol

      singular = .false.
      if (sum(g*g) <= grad2_floor) return
      call reduced_hessian(w, lambda, g, H, 16.0_wp*epsilon(1.0_wp), t1, t2, mu_min, mu_max, v_min, v_max, eig_tol)
      singular = min(abs(mu_min), abs(mu_max)) <= eig_tol
   end function singular_tangent_block

   !> Step of the full KKT system
   !>
   !> - solved unscaled on purpose: [[newton_solve_4x4]] fails once `||g||` falls
   !>   below about `4 sqrt(eps) ||w I - lambda Hess S||`, a multiplier blown up at
   !>   a flat gradient; such a start reaches no accepted point and ends here
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
      if (allocated(sc%iter)) then
         deallocate (sc%iter, sc%has_grad, sc%has_hess, &
                     sc%processed, sc%fnorm, sc%dz, sc%mu, sc%merit, sc%dmerit, &
                     sc%list, sc%lsl, sc%ptmp, sc%px, sc%plam, sc%p0, sc%p1, sc%p2)
      end if
      allocate (sc%iter(cap))
      allocate (sc%has_grad(cap), sc%has_hess(cap), sc%processed(cap), sc%fnorm(cap))
      allocate (sc%dz(4, cap), sc%mu(cap), sc%merit(cap), sc%dmerit(cap))
      allocate (sc%list(cap), sc%lsl(cap), sc%ptmp(cap))
      allocate (sc%px(3, cap), sc%plam(cap), sc%p0(cap), sc%p1(3, cap), sc%p2(6, cap))
      sc%capacity = cap
   end subroutine reserve_scratch

end module moist_cavity_drop_projection_newton
