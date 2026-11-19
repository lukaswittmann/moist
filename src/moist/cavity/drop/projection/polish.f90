!> Optimality verification and saddle escape for converged DROP projections
!>
!> Per-point routines ported from the single-anchor projector:
!>   - [[check_constrained_optimality]] -- KKT residuals plus the reduced
!>     Hessian on the tangent plane, third order when that is flat
!>   - [[riemannian_newton_escape]] -- walk off a saddle along the surface
!>   - [[retract_to_surface]] -- Newton along the normal back onto `S = 0`
!>
!> [[polish_batch]] applies them to every converged column of a batch, reusing
!> the jet Newton left there, and re-polishes escaped saddles with one more
!> batched Newton pass
!>
!> Run context read only, shared by the worker threads: traces gated by
!> verbosity 5 or `debug`, written straight to `ctx%unit` inside the
!> context's own `moist_context_io` critical; no timers, warnings or debug
!> messages
!>
!> Every level set evaluation goes through the LSF's `evaluate_batch`, the third
!> derivative of the degenerate test included; a point beyond the level set
!> reach is a recoverable failure of the routine that met it, never build-fatal
module moist_cavity_drop_projection_polish
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric
   use moist_math_packed, only: packed_contract2, packed_full_contract3
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, &
                                                point_failed_not_kkt, point_failed_saddle
   use moist_cavity_drop_projection_newton, only: newton_scratch_type, newton_solve_batch
   use moist_cavity_drop_projection_parameters, only: drop_solve_parameters_type, drop_polish_parameters_type
   implicit none(type, external)
   private

   public :: polish_batch
   public :: check_constrained_optimality, riemannian_newton_escape, retract_to_surface

   !> Diagonal and off-diagonal slots of a packed 3x3 Hessian
   integer, parameter :: diag_slots(3) = [1, 3, 6], offdiag_slots(3) = [2, 4, 5]
   !> Floor of the reduced-Hessian eigenvalue threshold
   real(wp), parameter :: tol_eig_base = 1.0e-10_wp
   !> Third directional derivative of the Lagrangian above which a flat point is an inflection
   real(wp), parameter :: tol_D3 = 1.0e-8_wp
   !> Largest objective drop an escape step must achieve; less is ok for a step with a smaller
   !> curvature-based gain at the threshold
   real(wp), parameter :: escape_drop_max = 1.0e-10_wp

contains

   !* ================================================================================= *!
   !*                                    Batch polish                                   *!
   !* ================================================================================= *!

   !> Verify every converged column in `first:last`, escaping saddles
   !>
   !> - minimum: `point_accepted`; no KKT point: `point_failed_not_kkt`
   !> - saddle: escaped, re-solved by one batched Newton pass and re-checked; a
   !>   point that is still a saddle afterwards is `point_failed_saddle`
   !> - a classification the third-order test cannot settle keeps the point,
   !>   as the single-anchor projector did
   !>
   !> @param[in]    ctx       Run context, read only; traces at verbosity 5
   !> @param[inout] batch     Points with Newton's final jets
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    cfg       Polish settings
   !> @param[in]    ncfg      Newton settings for the re-polish
   !> @param[in]    phi_alpha Anchor objective weight
   !> @param[inout] nscratch  Newton scratch
   !> @param[out]   error     LSF evaluation failure (build-fatal)
   !> @param[in]    first     First column (default 1)
   !> @param[in]    last      Last column (default `batch%n`)
   subroutine polish_batch(ctx, batch, lsf, cfg, ncfg, phi_alpha, nscratch, error, first, last)
      type(moist_context_type), intent(in) :: ctx
      type(projection_batch_type), intent(inout) :: batch
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      type(drop_polish_parameters_type), intent(in) :: cfg
      type(drop_solve_parameters_type), intent(in) :: ncfg
      real(wp), intent(in) :: phi_alpha
      type(newton_scratch_type), intent(inout) :: nscratch
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: first
      integer, intent(in), optional :: last

      type(error_type), allocatable :: class_error, escape_error
      real(wp) :: a(3), r_out(3), lambda_out
      integer :: j0, j1, j, n_escaped, index, ncap
      logical :: saddle, kkt

      j0 = 1
      j1 = batch%n
      if (present(first)) j0 = first
      if (present(last)) j1 = last
      if (j1 < j0) return
      if (.not. allocated(nscratch%escaped)) allocate (nscratch%escaped(max(16, j1)))
      if (size(nscratch%escaped) < j1) then
         ncap = max(j1, 2*size(nscratch%escaped))
         deallocate (nscratch%escaped)
         allocate (nscratch%escaped(ncap))
      end if
      nscratch%escaped(j0:j1) = .false.

      n_escaped = 0
      do j = j0, j1
         if (batch%status(j) /= point_converged) cycle
         a = batch%run_xyz(:, batch%run_of(j))
         index = batch%run_anchor(batch%run_of(j))
         call check_constrained_optimality(ctx, cfg, phi_alpha, lsf, batch%x(:, j), a, batch%lambda(j), class_error, &
                                           is_saddle=saddle, is_kkt=kkt, index=index, &
                                           jet_S=batch%S(j), jet_grad=batch%gradS(:, j), &
                                           jet_hess=batch%hessS(:, j), lsf_error=error)
         if (allocated(error)) return
         if (allocated(class_error)) then
            if (ctx%verbosity(5)) then
               !$omp critical (moist_context_io)
               write (ctx%unit, "(1x,i12,1x,a,a)") index, "[optimality] ", trim(class_error%message)
               !$omp end critical (moist_context_io)
            end if
            deallocate (class_error)
         end if
         if (.not. kkt) then
            batch%status(j) = point_failed_not_kkt
            cycle
         end if
         if (.not. saddle) then
            batch%status(j) = point_accepted
            cycle
         end if

         if (ctx%verbosity(5)) then
            !$omp critical (moist_context_io)
            write (ctx%unit, "(1x,i12,1x,a)") index, "[optimality] Saddle point detected, attempting escape ..."
            !$omp end critical (moist_context_io)
         end if
         call riemannian_newton_escape(ctx, cfg, phi_alpha, lsf, batch%x(:, j), a, batch%lambda(j), r_out, &
                                       lambda_out, escape_error, lsf_error=error)
         if (allocated(error)) return
         ! A walk that ran out of iterations still hands back its best point;
         ! Newton and the re-check below decide whether it is usable
         if (allocated(escape_error)) deallocate (escape_error)
         batch%x(:, j) = r_out
         batch%lambda(j) = lambda_out
         batch%status(j) = point_pending
         nscratch%escaped(j) = .true.
         n_escaped = n_escaped + 1
      end do
      if (n_escaped == 0) return

      call newton_solve_batch(batch, lsf, ncfg, phi_alpha, nscratch, error, first=j0, last=j1)
      if (allocated(error)) return

      do j = j0, j1
         if (.not. nscratch%escaped(j)) cycle
         if (batch%status(j) /= point_converged) cycle
         a = batch%run_xyz(:, batch%run_of(j))
         index = batch%run_anchor(batch%run_of(j))
         call check_constrained_optimality(ctx, cfg, phi_alpha, lsf, batch%x(:, j), a, batch%lambda(j), class_error, &
                                           is_saddle=saddle, is_kkt=kkt, index=index, &
                                           jet_S=batch%S(j), jet_grad=batch%gradS(:, j), &
                                           jet_hess=batch%hessS(:, j), lsf_error=error)
         if (allocated(error)) return
         if (allocated(class_error)) deallocate (class_error)
         if (.not. kkt) then
            batch%status(j) = point_failed_not_kkt
         else if (saddle) then
            if (ctx%verbosity(5)) then
               !$omp critical (moist_context_io)
               write (ctx%unit, "(1x,i12,1x,a)") index, "[optimality] Failed to escape saddle point!"
               !$omp end critical (moist_context_io)
            end if
            batch%status(j) = point_failed_saddle
         else
            batch%status(j) = point_accepted
         end if
      end do
   end subroutine polish_batch

   !* ================================================================================= *!
   !*                            Check constrained optimality                           *!
   !* ================================================================================= *!

   !> Whether a converged KKT point is a constrained local minimum or a saddle
   !>
   !> - reduced Hessian eigenvalues on the tangent plane, not the bordered KKT
   !>   matrix
   !> - a flat second-order test falls through to `D^3 L[d,d,d]` along the flat
   !>   direction; nonzero means an inflection, treated as an escapable saddle
   !> - flat to third order raises `error`, classification unreliable
   !> - `r_star` beyond the level set reach: `is_kkt = .false.` and `error`
   !>
   !> @param[in]    ctx          Run context, read only
   !> @param[in]    cfg          Polish settings
   !> @param[in]    phi_alpha    Anchor objective weight
   !> @param[inout] lsf          Thread-local level set
   !> @param[in]    r_star       Candidate point
   !> @param[in]    anchor       Anchor
   !> @param[in]    lambda       Multiplier at `r_star`
   !> @param[out]   error        Unreliable classification, or saddle under strict mode
   !> @param[out]   is_saddle    Saddle or inflection
   !> @param[out]   is_kkt       Feasible and stationary
   !> @param[out]   is_converged Strict local minimum
   !> @param[in]    index        Anchor index for traces
   !> @param[in]    jet_S        Level set value at `r_star`, skips the evaluation
   !> @param[in]    jet_grad     Level set gradient at `r_star`
   !> @param[in]    jet_hess     Level set Hessian at `r_star`, packed (6)
   !> @param[out]   lsf_error    LSF evaluation failure; folded into `error` when absent
   subroutine check_constrained_optimality(ctx, cfg, phi_alpha, lsf, r_star, anchor, lambda, error, &
                                           is_saddle, is_kkt, is_converged, index, &
                                           jet_S, jet_grad, jet_hess, lsf_error)
      type(moist_context_type), intent(in) :: ctx
      type(drop_polish_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_star(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lambda
      type(error_type), allocatable, intent(out) :: error
      logical, intent(out), optional :: is_saddle
      logical, intent(out), optional :: is_kkt
      logical, intent(out), optional :: is_converged
      integer, intent(in), optional :: index
      real(wp), intent(in), optional :: jet_S
      real(wp), intent(in), optional :: jet_grad(3)
      real(wp), intent(in), optional :: jet_hess(:)
      type(error_type), allocatable, intent(out), optional :: lsf_error

      real(wp), parameter :: tol_feas = 1.0e-6_wp
      real(wp), parameter :: tol_stat = 1.0e-6_wp
      real(wp), parameter :: tol_grad = 1.0e-10_wp

      type(error_type), allocatable :: eval_error
      real(wp) :: S_val, grad_S(3), hess_S(6), grad_phi(3)
      real(wp) :: normal_vec(3), norm_grad_S, res_feas, res_stat, grad_L(3)
      real(wp) :: H_lagrangian(6), R_reduced(2, 2), T_mat(3, 2), t1(3), t2(3), Ht1(3), Ht2(3)
      real(wp) :: mu_min, mu_max, v_min(2), v_max(2), eig_tol
      real(wp) :: d_degen(3), D3_ddd
      integer :: id
      logical :: unreached

      if (present(is_saddle)) is_saddle = .false.
      if (present(is_kkt)) is_kkt = .true.
      if (present(is_converged)) is_converged = .false.
      id = 0
      if (present(index)) id = index

      if (present(jet_S) .and. present(jet_grad) .and. present(jet_hess)) then
         S_val = jet_S
         grad_S = jet_grad
         hess_S = jet_hess
      else
         call point_jet(lsf, r_star, 2, S_val, unreached, eval_error, grad_S, hess_S)
         if (allocated(eval_error)) then
            call hand_over(eval_error, error, lsf_error)
            if (present(is_kkt)) is_kkt = .false.
            return
         end if
         if (unreached) then
            call fatal_error(error, "KKT check: candidate point beyond the level set reach")
            if (present(is_kkt)) is_kkt = .false.
            return
         end if
      end if

      res_feas = abs(S_val)
      norm_grad_S = norm2(grad_S)
      if (norm_grad_S < tol_grad) then
         if (ctx%verbosity(5)) then
            !$omp critical (moist_context_io)
            write (ctx%unit, "(1x,i12,1x,a,es12.4)") id, &
               "[KKT check] Degenerate LSF gradient ||grad_S|| = ", norm_grad_S
            !$omp end critical (moist_context_io)
         end if
         ! Cannot classify, but an infeasible point is no KKT point
         if (present(is_kkt)) is_kkt = res_feas <= tol_feas
         return
      end if

      grad_phi = phi_alpha*(r_star - anchor)
      grad_L = grad_phi - lambda*grad_S
      res_stat = norm2(grad_L)

      if (ctx%debug) then
         ! One block, so no other thread's line lands inside it
         !$omp critical (moist_context_io)
         write (ctx%unit, "(1x,a)") "KKT verification:"
         write (ctx%unit, "(3x,a,es12.4,a,es12.4,a)") &
            "Feasibility   |S(r*)| = ", res_feas, "  (tol: ", tol_feas, ")"
         write (ctx%unit, "(3x,a,es12.4,a,es12.4,a)") &
            "Stationarity  ||dL||  = ", res_stat, "  (tol: ", tol_stat, ")"
         !$omp end critical (moist_context_io)
      end if

      if (res_feas > tol_feas .or. res_stat > tol_stat) then
         if (ctx%verbosity(5) .or. ctx%debug) then
            !$omp critical (moist_context_io)
            write (ctx%unit, "(1x,i12,1x,a,es12.4,a,es12.4)") id, &
               "[KKT check] Not at KKT point: |S| = ", res_feas, ", ||dL|| = ", res_stat
            !$omp end critical (moist_context_io)
         end if
         if (present(is_kkt)) is_kkt = .false.
         return
      end if

      !> Reduced Hessian R = T^T (Hess phi - lambda Hess S) T on the tangent plane,
      !> packed `H_L` with `Hess phi = phi_alpha I` on the diagonal slots
      normal_vec = grad_S/norm_grad_S
      call setup_tangent_frame(normal_vec, t1, t2)
      T_mat(:, 1) = t1
      T_mat(:, 2) = t2
      H_lagrangian = -lambda*hess_S
      H_lagrangian(diag_slots) = phi_alpha - lambda*hess_S(diag_slots)
      call packed_contract2(H_lagrangian, t1, Ht1, 3)
      call packed_contract2(H_lagrangian, t2, Ht2, 3)
      R_reduced(1, 1) = dot_product(t1, Ht1)
      R_reduced(1, 2) = dot_product(t1, Ht2)
      R_reduced(2, 2) = dot_product(t2, Ht2)
      R_reduced(2, 1) = R_reduced(1, 2)
      call eig_2x2_symmetric(R_reduced(1, 1), R_reduced(1, 2), R_reduced(2, 2), &
                             mu_min, mu_max, v_min, v_max)

      eig_tol = saddle_eig_threshold(cfg%saddle_eig_tol, H_lagrangian)

      if (ctx%debug) then
         !$omp critical (moist_context_io)
         write (ctx%unit, "(3x,a,3es14.6)") "Reduced Hessian eigenvalues: ", mu_min, mu_max
         !$omp end critical (moist_context_io)
      end if

      if (mu_min < -eig_tol) then
         if (present(is_saddle)) is_saddle = .true.
         if (cfg%strict_minimum_required .and. .not. present(is_saddle)) then
            call fatal_error(error, "Saddle point detected but strict minimum required")
         end if
         return
      end if

      if (mu_min > eig_tol) then
         if (present(is_converged)) is_converged = .true.
         return
      end if

      !> Flat second order: D^3 L[d,d,d] = -lambda S_ijk d_i d_j d_k along d = T v_min
      d_degen = matmul(T_mat, v_min)
      call lagrangian_third_derivative(lsf, r_star, lambda, d_degen, D3_ddd, unreached, eval_error)
      if (allocated(eval_error)) then
         call hand_over(eval_error, error, lsf_error)
         if (present(is_kkt)) is_kkt = .false.
         return
      end if
      if (unreached) then
         call fatal_error(error, "KKT check: candidate point beyond the level set reach")
         if (present(is_kkt)) is_kkt = .false.
         return
      end if

      if (ctx%debug) then
         !$omp critical (moist_context_io)
         write (ctx%unit, "(3x,a,2es16.8,a,es16.8)") "Reduced eigs = ", mu_min, mu_max, &
            ", D3 L[d,d,d] = ", D3_ddd
         !$omp end critical (moist_context_io)
      end if

      if (abs(D3_ddd) > tol_D3) then
         if (present(is_saddle)) is_saddle = .true.
         if (cfg%strict_minimum_required .and. .not. present(is_saddle)) then
            call fatal_error(error, "Inflection point detected but strict minimum required")
         end if
      else
         call fatal_error(error, "Degenerate to third order; classification unreliable")
      end if
   end subroutine check_constrained_optimality

   !* ================================================================================= *!
   !*                    Riemannian Newton escape from saddle points                    *!
   !* ================================================================================= *!

   !> Walk off a saddle along the surface, then descend to a minimum
   !>
   !> - at a saddle, or at a flat point whose third derivative makes it an inflection
   !>   (the classification of [[check_constrained_optimality]]): fixed step along
   !>   the negative-curvature or flat tangent direction; the first-order change
   !>   there is roundoff at a KKT point, so both signs are tried and the lower
   !>   objective kept
   !> - a flat point with zero third derivative ends the walk, as the check keeps it
   !> - elsewhere: Riemannian Newton (or steepest descent when nearly singular)
   !> - stops in the basin of a minimum, at a tangential gradient of `1e-6`;
   !>   [[polish_batch]] converges the point with the batched Newton like any other
   !>   column
   !> - every trial is retracted onto `S = 0` before the acceptance test; an escape
   !>   step must lower `phi`, a Newton or descent step must pass the Armijo test
   !> - `r_out`/`lambda_out` always hold the last accepted point, also on error
   !> - an iterate beyond the level set reach ends the walk with `error`
   !>
   !> @param[in]    ctx         Run context, read only
   !> @param[in]    cfg         Polish settings
   !> @param[in]    phi_alpha   Anchor objective weight
   !> @param[inout] lsf         Thread-local level set
   !> @param[in]    r_init      Start point
   !> @param[in]    anchor      Anchor
   !> @param[in]    lambda_init Multiplier at `r_init`
   !> @param[out]   r_out       Final point
   !> @param[out]   lambda_out  Multiplier at `r_out`
   !> @param[out]   error       Degenerate gradient, unreached iterate or iteration budget exhausted
   !> @param[out]   lsf_error   LSF evaluation failure; folded into `error` when absent
   subroutine riemannian_newton_escape(ctx, cfg, phi_alpha, lsf, r_init, anchor, lambda_init, r_out, lambda_out, &
                                       error, lsf_error)
      type(moist_context_type), intent(in) :: ctx
      type(drop_polish_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_init(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lambda_init
      real(wp), intent(out) :: r_out(3)
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable, intent(out), optional :: lsf_error

      !> Hand-over gradient: below about sqrt(eps) the objective-based line search
      !> cannot resolve progress, and a point handed over below Newton's residual
      !> tolerance would be kept unrefined, off by `||F|| / mu` on a soft point
      real(wp), parameter :: tol_grad_tan = 1.0e-6_wp

      type(error_type), allocatable :: eval_error, retract_error
      real(wp) :: r_curr(3), r_trial(3), r_retracted(3), lambda_curr, lambda_trial
      real(wp) :: S_val, grad_S(3), hess_S(6), grad_phi(3)
      real(wp) :: phi_curr, phi_trial, normal_vec(3), norm_grad_S
      real(wp) :: T_mat(3, 2), t1(3), t2(3), g_tan(2), H_lagrangian(6), R_reduced(2, 2), Ht1(3), Ht2(3)
      real(wp) :: delta_r(3), alpha, directional_deriv, norm_g_tan
      real(wp) :: mu_min, mu_max, v_min(2), v_max(2), s_tan(2), eig_tol, D3_ddd
      real(wp) :: r_cand(3), lambda_cand, phi_cand, drop_min
      integer :: iter, ls_iter, isign, nsign
      logical :: converged, ls_success, at_saddle, escaping, have_trial, unreached

      r_curr = r_init
      lambda_curr = lambda_init
      r_out = r_curr
      lambda_out = lambda_curr
      converged = .false.
      phi_curr = 0.5_wp*phi_alpha*sum((r_curr - anchor)**2)

      do iter = 0, cfg%escape_max_iter
         call point_jet(lsf, r_curr, 2, S_val, unreached, eval_error, grad_S, hess_S)
         if (allocated(eval_error)) then
            call hand_over(eval_error, error, lsf_error)
            return
         end if
         if (unreached) then
            call fatal_error(error, "Riemannian escape: iterate beyond the level set reach")
            return
         end if
         grad_phi = phi_alpha*(r_curr - anchor)

         norm_grad_S = norm2(grad_S)
         if (norm_grad_S < 1.0e-10_wp) then
            call fatal_error(error, "Riemannian escape: degenerate LSF gradient")
            return
         end if
         normal_vec = grad_S/norm_grad_S
         call setup_tangent_frame(normal_vec, t1, t2)
         T_mat(:, 1) = t1
         T_mat(:, 2) = t2

         g_tan(1) = dot_product(t1, grad_phi)
         g_tan(2) = dot_product(t2, grad_phi)
         norm_g_tan = norm2(g_tan)

         H_lagrangian = -lambda_curr*hess_S
         H_lagrangian(diag_slots) = phi_alpha - lambda_curr*hess_S(diag_slots)
         call packed_contract2(H_lagrangian, t1, Ht1, 3)
         call packed_contract2(H_lagrangian, t2, Ht2, 3)
         R_reduced(1, 1) = dot_product(t1, Ht1)
         R_reduced(1, 2) = dot_product(t1, Ht2)
         R_reduced(2, 2) = dot_product(t2, Ht2)
         R_reduced(2, 1) = R_reduced(1, 2)
         call eig_2x2_symmetric(R_reduced(1, 1), R_reduced(1, 2), R_reduced(2, 2), &
                                mu_min, mu_max, v_min, v_max)
         ! Same threshold as the optimality check, so both agree on what a saddle is
         eig_tol = saddle_eig_threshold(cfg%saddle_eig_tol, H_lagrangian)
         at_saddle = mu_min < -eig_tol

         escaping = at_saddle
         if (norm_g_tan < tol_grad_tan .and. .not. at_saddle) then
            if (mu_min > eig_tol) then
               converged = .true.
               exit
            end if
            ! Flat second order: the third order decides, as in the check
            call lagrangian_third_derivative(lsf, r_curr, lambda_curr, matmul(T_mat, v_min), D3_ddd, &
                                             unreached, eval_error)
            if (allocated(eval_error)) then
               call hand_over(eval_error, error, lsf_error)
               return
            end if
            if (unreached) then
               call fatal_error(error, "Riemannian escape: iterate beyond the level set reach")
               return
            end if
            if (abs(D3_ddd) <= tol_D3) then
               converged = .true.
               exit
            end if
            escaping = .true.
         end if

         if (escaping) then
            ! Unsigned: `grad phi . delta_r` is roundoff at a KKT point
            delta_r = T_mat(:, 1)*v_min(1) + T_mat(:, 2)*v_min(2)
            directional_deriv = 0.0_wp
            alpha = cfg%escape_step
         else
            if (mu_min > eig_tol) then
               ! s = -V diag(1/mu) V^T g
               s_tan(1) = -(v_min(1)*g_tan(1) + v_min(2)*g_tan(2))/mu_min*v_min(1) &
                          - (v_max(1)*g_tan(1) + v_max(2)*g_tan(2))/mu_max*v_max(1)
               s_tan(2) = -(v_min(1)*g_tan(1) + v_min(2)*g_tan(2))/mu_min*v_min(2) &
                          - (v_max(1)*g_tan(1) + v_max(2)*g_tan(2))/mu_max*v_max(2)
            else
               s_tan = -g_tan/max(norm_g_tan, 1.0e-10_wp)
            end if
            delta_r = T_mat(:, 1)*s_tan(1) + T_mat(:, 2)*s_tan(2)
            directional_deriv = dot_product(grad_phi, delta_r)
            alpha = cfg%escape_alpha_init
         end if

         ! An escape step tries both signs and keeps the lower retracted objective
         nsign = 1
         if (escaping) nsign = 2
         ls_success = .false.
         do ls_iter = 1, cfg%escape_max_ls_iter
            have_trial = .false.
            do isign = 1, nsign
               r_trial = r_curr + merge(alpha, -alpha, isign == 1)*delta_r
               call retract_to_surface(cfg, phi_alpha, lsf, r_trial, anchor, r_cand, lambda_cand, retract_error, &
                                       lsf_error=eval_error)
               if (allocated(eval_error)) then
                  call hand_over(eval_error, error, lsf_error)
                  return
               end if
               if (allocated(retract_error)) then
                  deallocate (retract_error)
                  cycle
               end if
               phi_cand = 0.5_wp*phi_alpha*sum((r_cand - anchor)**2)
               if (have_trial) then
                  if (phi_cand >= phi_trial) cycle
               end if
               r_retracted = r_cand
               lambda_trial = lambda_cand
               phi_trial = phi_cand
               have_trial = .true.
            end do
            if (.not. have_trial) then
               alpha = alpha*cfg%escape_alpha_reduce
               cycle
            end if

            if (escaping) then
               ! No first-order decrease to demand off a stationary point: the objective must drop
               drop_min = max(16.0_wp*epsilon(phi_curr)*abs(phi_curr), &
                              min(escape_drop_max, 0.25_wp*eig_tol*alpha**2))
               ls_success = phi_trial < phi_curr - drop_min
            else
               ! Sufficient decrease; a bare decrease accepts the overshoot of a long Newton step
               ls_success = phi_trial <= phi_curr + cfg%escape_armijo_c*alpha*directional_deriv
            end if
            if (ls_success) exit
            alpha = alpha*cfg%escape_alpha_reduce
         end do

         if (.not. ls_success) exit

         r_curr = r_retracted
         lambda_curr = lambda_trial
         phi_curr = phi_trial
         r_out = r_curr
         lambda_out = lambda_curr

         if (ctx%debug) then
            !$omp critical (moist_context_io)
            write (ctx%unit, "(i4,1x,es12.4,1x,f12.6,1x,es12.4,1x,es10.2,1x,es12.4,1x,l1)") &
               iter, phi_curr, norm2(r_curr - anchor), norm_g_tan, alpha, mu_min, escaping
            !$omp end critical (moist_context_io)
         end if
      end do

      if (.not. converged .and. iter >= cfg%escape_max_iter) then
         call fatal_error(error, "Riemannian Newton did not converge within max iterations")
      end if
   end subroutine riemannian_newton_escape

   !* ================================================================================= *!
   !*                                 Surface retraction                                *!
   !* ================================================================================= *!

   !> Retract a point onto `S = 0` by Newton steps along the gradient
   !>
   !> - the multiplier is the least-squares stationarity estimate
   !>   `grad phi . grad S / ||grad S||^2` at the retracted point
   !>
   !> @param[in]    cfg        Polish settings
   !> @param[in]    phi_alpha  Anchor objective weight
   !> @param[inout] lsf        Thread-local level set
   !> @param[in]    r_init     Start point
   !> @param[in]    anchor     Anchor
   !> @param[out]   r_out      Point with `|S| < 1e-12`
   !> @param[out]   lambda_out Multiplier at `r_out`
   !> @param[out]   error      Degenerate gradient, unreached trial or iteration budget exhausted
   !> @param[out]   lsf_error  LSF evaluation failure; folded into `error` when absent
   subroutine retract_to_surface(cfg, phi_alpha, lsf, r_init, anchor, r_out, lambda_out, error, lsf_error)
      type(drop_polish_parameters_type), intent(in) :: cfg
      real(wp), intent(in) :: phi_alpha
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_init(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(out) :: r_out(3)
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable, intent(out), optional :: lsf_error

      real(wp), parameter :: tol_S = 1.0e-12_wp

      type(error_type), allocatable :: eval_error
      real(wp) :: r_curr(3), S_val, grad_S(3), norm_grad_S2
      integer :: iter
      logical :: unreached

      r_curr = r_init
      r_out = r_init
      lambda_out = 0.0_wp

      do iter = 1, cfg%retract_max_iter
         call point_jet(lsf, r_curr, 1, S_val, unreached, eval_error, grad_S)
         if (allocated(eval_error)) then
            call hand_over(eval_error, error, lsf_error)
            return
         end if
         ! Recoverable: the escape walk backs off its step
         if (unreached) then
            call fatal_error(error, "Retraction failed: trial point beyond the level set reach")
            return
         end if
         norm_grad_S2 = dot_product(grad_S, grad_S)

         if (abs(S_val) < tol_S) then
            if (norm_grad_S2 < 1.0e-14_wp) then
               call fatal_error(error, "Retraction failed: degenerate LSF gradient at converged point")
               return
            end if
            r_out = r_curr
            lambda_out = dot_product(phi_alpha*(r_curr - anchor), grad_S)/norm_grad_S2
            return
         end if

         if (norm_grad_S2 < 1.0e-14_wp) then
            call fatal_error(error, "Retraction failed: degenerate LSF gradient")
            return
         end if
         r_curr = r_curr - (S_val/norm_grad_S2)*grad_S
      end do

      call fatal_error(error, "Retraction to surface did not converge")
   end subroutine retract_to_surface

   !> Third directional derivative of the Lagrangian, `D^3 L[d,d,d] = -lambda S_ijk d_i d_j d_k`
   !>
   !> - `Hess phi` is constant, so only the level set contributes
   !>
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    r         Surface point
   !> @param[in]    lambda    Multiplier at `r`
   !> @param[in]    d         Tangent direction
   !> @param[out]   D3        Third directional derivative
   !> @param[out]   unreached Point beyond the level set reach; `D3` is then zero and meaningless
   !> @param[out]   error     LSF evaluation failure
   subroutine lagrangian_third_derivative(lsf, r, lambda, d, D3, unreached, error)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r(3)
      real(wp), intent(in) :: lambda
      real(wp), intent(in) :: d(3)
      real(wp), intent(out) :: D3
      logical, intent(out) :: unreached
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: S3, grad3(3), hess3(6), third3(10)

      D3 = 0.0_wp
      call point_jet(lsf, r, 3, S3, unreached, error, grad3, hess3, third3)
      if (allocated(error) .or. unreached) return
      D3 = -lambda*packed_full_contract3(third3, d, 3)
   end subroutine lagrangian_third_derivative

   !> Reduced-Hessian eigenvalue threshold of the saddle classification
   !>
   !>     max(tol_eig_base, rel_tol ||H_L||_F)
   !>
   !> - shared by [[check_constrained_optimality]] and [[riemannian_newton_escape]]
   !> - Frobenius norm of the packed Hessian, off-diagonal slots twice
   !>
   !> @param[in] rel_tol      Threshold relative to `||H_L||_F`
   !> @param[in] H_lagrangian Lagrangian Hessian, packed (6)
   pure function saddle_eig_threshold(rel_tol, H_lagrangian) result(eig_tol)
      real(wp), intent(in) :: rel_tol
      real(wp), intent(in) :: H_lagrangian(6)
      real(wp) :: eig_tol

      real(wp) :: HL_norm

      HL_norm = sqrt(sum(H_lagrangian(diag_slots)**2) + 2.0_wp*sum(H_lagrangian(offdiag_slots)**2))
      eig_tol = max(tol_eig_base, rel_tol*HL_norm)
   end function saddle_eig_threshold

   !* ================================================================================= *!
   !*                            Single-point level set jets                            *!
   !* ================================================================================= *!

   !> Level set jet at one point through the batched entry
   !>
   !> - `unreached`: no atom within the screening reach; the outputs are then
   !>   unspecified and must not be read
   !>
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    point     Evaluation point
   !> @param[in]    order     Highest derivative (0..3)
   !> @param[out]   S         Value
   !> @param[out]   unreached Point beyond the level set reach
   !> @param[out]   error     LSF evaluation failure
   !> @param[out]   grad      Gradient, required for order >= 1
   !> @param[out]   hess      Hessian, packed (6), required for order >= 2
   !> @param[out]   third     Third derivatives, packed (10), required for order 3
   subroutine point_jet(lsf, point, order, S, unreached, error, grad, hess, third)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: order
      real(wp), intent(out) :: S
      logical, intent(out) :: unreached
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(out), optional :: grad(3)
      real(wp), intent(out), optional :: hess(:)
      real(wp), intent(out), optional :: third(10)

      !> One-column batch; Hessian packed (6, 1), third derivatives packed (10, 1)
      real(wp) :: p(3, 1), v(1), g(3, 1), h(6, 1), t3(10, 1)

      unreached = .false.
      p(:, 1) = point
      select case (order)
      case (0)
         call lsf%evaluate_batch(p, 0, v, error=error)
      case (1)
         call lsf%evaluate_batch(p, 1, v, lsf1_r=g, error=error)
      case (2)
         call lsf%evaluate_batch(p, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=error)
      case default
         call lsf%evaluate_batch(p, 3, v, lsf1_r=g, lsf2_rr_packed=h, lsf3_rrr=t3, error=error)
      end select
      if (allocated(error)) return
      unreached = lsf%unreached(1)
      if (unreached) return
      S = v(1)
      if (order >= 1 .and. present(grad)) grad = g(:, 1)
      if (order >= 2 .and. present(hess)) hess = h(:, 1)
      if (order >= 3 .and. present(third)) third = t3(:, 1)
   end subroutine point_jet

   !> Route an LSF failure to `lsf_error` when the caller asked for it
   !>
   !> @param[inout] eval_error LSF failure, consumed
   !> @param[out]   error      Classification channel
   !> @param[out]   lsf_error  Dedicated LSF channel
   subroutine hand_over(eval_error, error, lsf_error)
      type(error_type), allocatable, intent(inout) :: eval_error
      type(error_type), allocatable, intent(inout) :: error
      type(error_type), allocatable, intent(inout), optional :: lsf_error

      if (present(lsf_error)) then
         call move_alloc(eval_error, lsf_error)
      else
         call move_alloc(eval_error, error)
      end if
   end subroutine hand_over

end module moist_cavity_drop_projection_polish
