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
!> Every level set evaluation goes through the LSF's `evaluate_batch`, the third
!> derivative of the degenerate test included; a point beyond the level set
!> reach is a recoverable failure of the routine that met it, never build-fatal
module moist_cavity_drop_projection_polish
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric
   use moist_math_packed, only: packed_contract2, unpack_sym_tensor3
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, &
                                                point_failed_not_kkt, point_failed_saddle
   use moist_cavity_drop_projection_newton, only: newton_config_type, newton_scratch_type, &
                                                 newton_solve_batch
   implicit none(type, external)
   private

   public :: polish_config_type, polish_batch
   public :: check_constrained_optimality, riemannian_newton_escape, retract_to_surface

   !> Diagonal and off-diagonal slots of a packed 3x3 Hessian
   integer, parameter :: diag_slots(3) = [1, 3, 6], offdiag_slots(3) = [2, 4, 5]

   ! TODO: we should also think about making this a general solver or at least some type of Riemannian Newton module

   !> Polish settings
   ! TODO: good idea to introduce a new type here with drop parameters again?
   type :: polish_config_type
      !> Anchor objective weight
      real(wp) :: phi_alpha = 0.5_wp
      !> Riemannian iterations
      integer :: riemann_max_iter = 150
      !> Line-search trials per Riemannian iteration
      integer :: riemann_max_ls_iter = 20
      !> Initial line-search step of a non-saddle Riemannian update
      real(wp) :: riemann_alpha_init = 2.0_wp
      !> Reduction of a rejected line-search step
      real(wp) :: riemann_alpha_reduce = 0.75_wp
      !> Armijo coefficient of the Riemannian line search
      real(wp) :: riemann_armijo_c = 1.0e-4_wp
      !> Fixed step along negative curvature when leaving a saddle
      real(wp) :: riemann_escape_step = 0.5_wp
      !> Reduced-Hessian eigenvalue below which a point counts as a saddle
      real(wp) :: riemann_saddle_eig_tol = 1.0e-8_wp
      !> Retraction iterations
      integer :: retract_max_iter = 50
      !> Report a saddle as an error when the caller does not ask for `is_saddle`
      logical :: strict_minimum_required = .true.
      !> Stability tables and Riemannian traces
      logical :: debug = .false.
      !> Run-context verbosity; per-point traces print at 5
      integer :: verbosity = 1
      !> Output unit of traces
      integer :: unit = output_unit
   end type polish_config_type

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
   !> @param[inout] batch     Points with Newton's final jets
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    cfg       Polish settings
   !> @param[in]    ncfg      Newton settings for the re-polish
   !> @param[inout] nscratch  Newton scratch
   !> @param[out]   error     LSF evaluation failure (build-fatal)
   !> @param[in]    first     First column (default 1)
   !> @param[in]    last      Last column (default `batch%n`)
   subroutine polish_batch(batch, lsf, cfg, ncfg, nscratch, error, first, last)
      type(projection_batch_type), intent(inout) :: batch
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      type(polish_config_type), intent(in) :: cfg
      type(newton_config_type), intent(in) :: ncfg
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
         call check_constrained_optimality(cfg, lsf, batch%x(:, j), a, batch%lambda(j), class_error, &
                                           is_saddle=saddle, is_kkt=kkt, index=index, &
                                           jet_S=batch%S(j), jet_grad=batch%gradS(:, j), &
                                           jet_hess=batch%hessS(:, j), lsf_error=error)
         if (allocated(error)) return
         if (allocated(class_error)) then
            if (cfg%verbosity >= 5) write (cfg%unit, "(1x,i12,1x,a,a)") index, "[optimality] ", &
               trim(class_error%message)
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

         if (cfg%verbosity >= 5) write (cfg%unit, "(1x,i12,1x,a)") index, &
            "[optimality] Saddle point detected, attempting escape ..."
         call riemannian_newton_escape(cfg, lsf, batch%x(:, j), a, batch%lambda(j), r_out, lambda_out, &
                                       escape_error, lsf_error=error)
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

      call newton_solve_batch(batch, lsf, ncfg, nscratch, error, first=j0, last=j1)
      if (allocated(error)) return

      do j = j0, j1
         if (.not. nscratch%escaped(j)) cycle
         if (batch%status(j) /= point_converged) cycle
         a = batch%run_xyz(:, batch%run_of(j))
         index = batch%run_anchor(batch%run_of(j))
         call check_constrained_optimality(cfg, lsf, batch%x(:, j), a, batch%lambda(j), class_error, &
                                           is_saddle=saddle, is_kkt=kkt, index=index, &
                                           jet_S=batch%S(j), jet_grad=batch%gradS(:, j), &
                                           jet_hess=batch%hessS(:, j), lsf_error=error)
         if (allocated(error)) return
         if (allocated(class_error)) deallocate (class_error)
         if (.not. kkt) then
            batch%status(j) = point_failed_not_kkt
         else if (saddle) then
            if (cfg%verbosity >= 5) write (cfg%unit, "(1x,i12,1x,a)") index, &
               "[optimality] Failed to escape saddle point!"
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
   !> @param[in]    cfg          Polish settings
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
   subroutine check_constrained_optimality(cfg, lsf, r_star, anchor, lambda, error, &
                                           is_saddle, is_kkt, is_converged, index, &
                                           jet_S, jet_grad, jet_hess, lsf_error)
      type(polish_config_type), intent(in) :: cfg
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

      real(wp), parameter :: tol_D3 = 1.0e-8_wp
      real(wp), parameter :: tol_feas = 1.0e-6_wp
      real(wp), parameter :: tol_stat = 1.0e-6_wp
      real(wp), parameter :: tol_grad = 1.0e-10_wp
      real(wp), parameter :: tol_eig_base = 1.0e-10_wp

      type(error_type), allocatable :: eval_error
      real(wp) :: S_val, grad_S(3), hess_S(6), grad_phi(3)
      real(wp) :: normal_vec(3), norm_grad_S, res_feas, res_stat, grad_L(3)
      real(wp) :: H_lagrangian(6), R_reduced(2, 2), T_mat(3, 2), t1(3), t2(3), Ht1(3), Ht2(3)
      real(wp) :: mu_min, mu_max, v_min(2), v_max(2), eig_tol, HL_norm
      real(wp) :: lsf3_S(3, 3, 3), d_degen(3), D3_ddd, S3, grad3(3), hess3(6)
      integer :: i, j, k, id
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
         if (cfg%verbosity >= 5) then
            write (cfg%unit, "(1x,i12,1x,a,es12.4)") id, &
               "[KKT check] Degenerate LSF gradient ||grad_S|| = ", norm_grad_S
         end if
         ! Cannot classify, but an infeasible point is no KKT point
         if (present(is_kkt)) is_kkt = res_feas <= tol_feas
         return
      end if

      grad_phi = cfg%phi_alpha*(r_star - anchor)
      grad_L = grad_phi - lambda*grad_S
      res_stat = norm2(grad_L)

      if (cfg%debug) then
         write (cfg%unit, "(1x,a)") "KKT verification:"
         write (cfg%unit, "(3x,a,es12.4,a,es12.4,a)") &
            "Feasibility   |S(r*)| = ", res_feas, "  (tol: ", tol_feas, ")"
         write (cfg%unit, "(3x,a,es12.4,a,es12.4,a)") &
            "Stationarity  ||dL||  = ", res_stat, "  (tol: ", tol_stat, ")"
      end if

      if (res_feas > tol_feas .or. res_stat > tol_stat) then
         if (cfg%verbosity >= 5 .or. cfg%debug) then
            write (cfg%unit, "(1x,i12,1x,a,es12.4,a,es12.4)") id, &
               "[KKT check] Not at KKT point: |S| = ", res_feas, ", ||dL|| = ", res_stat
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
      H_lagrangian(diag_slots) = cfg%phi_alpha - lambda*hess_S(diag_slots)
      call packed_contract2(H_lagrangian, t1, Ht1, 3)
      call packed_contract2(H_lagrangian, t2, Ht2, 3)
      R_reduced(1, 1) = dot_product(t1, Ht1)
      R_reduced(1, 2) = dot_product(t1, Ht2)
      R_reduced(2, 2) = dot_product(t2, Ht2)
      R_reduced(2, 1) = R_reduced(1, 2)
      call eig_2x2_symmetric(R_reduced(1, 1), R_reduced(1, 2), R_reduced(2, 2), &
                             mu_min, mu_max, v_min, v_max)

      ! Frobenius norm, off-diagonal slots twice
      HL_norm = sqrt(sum(H_lagrangian(diag_slots)**2) + 2.0_wp*sum(H_lagrangian(offdiag_slots)**2))
      eig_tol = max(tol_eig_base, 1.0e-8_wp*HL_norm)

      if (cfg%debug) then
         write (cfg%unit, "(3x,a,3es14.6)") "Reduced Hessian eigenvalues: ", mu_min, mu_max
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
      call point_jet(lsf, r_star, 3, S3, unreached, eval_error, grad3, hess3, lsf3_S)
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
      D3_ddd = 0.0_wp
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               D3_ddd = D3_ddd + lsf3_S(i, j, k)*d_degen(i)*d_degen(j)*d_degen(k)
            end do
         end do
      end do
      D3_ddd = -lambda*D3_ddd

      if (cfg%debug) then
         write (cfg%unit, "(3x,a,2es16.8,a,es16.8)") "Reduced eigs = ", mu_min, mu_max, &
            ", D3 L[d,d,d] = ", D3_ddd
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
   !> - at a saddle: fixed step along the negative-curvature tangent direction,
   !>   signed to decrease `phi`
   !> - elsewhere: Riemannian Newton (or steepest descent when nearly singular)
   !> - every trial is retracted onto `S = 0` before the acceptance test
   !> - `r_out`/`lambda_out` always hold the last accepted point, also on error
   !> - an iterate beyond the level set reach ends the walk with `error`
   !>
   !> @param[in]    cfg         Polish settings
   !> @param[inout] lsf         Thread-local level set
   !> @param[in]    r_init      Start point
   !> @param[in]    anchor      Anchor
   !> @param[in]    lambda_init Multiplier at `r_init`
   !> @param[out]   r_out       Final point
   !> @param[out]   lambda_out  Multiplier at `r_out`
   !> @param[out]   error       Degenerate gradient, unreached iterate or iteration budget exhausted
   !> @param[out]   lsf_error   LSF evaluation failure; folded into `error` when absent
   subroutine riemannian_newton_escape(cfg, lsf, r_init, anchor, lambda_init, r_out, lambda_out, &
                                       error, lsf_error)
      type(polish_config_type), intent(in) :: cfg
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: r_init(3)
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lambda_init
      real(wp), intent(out) :: r_out(3)
      real(wp), intent(out) :: lambda_out
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable, intent(out), optional :: lsf_error

      real(wp), parameter :: tol_grad_tan = 1.0e-8_wp

      type(error_type), allocatable :: eval_error, retract_error
      real(wp) :: r_curr(3), r_trial(3), r_retracted(3), lambda_curr, lambda_trial
      real(wp) :: S_val, grad_S(3), hess_S(6), grad_phi(3)
      real(wp) :: phi_curr, phi_trial, normal_vec(3), norm_grad_S
      real(wp) :: T_mat(3, 2), t1(3), t2(3), g_tan(2), H_lagrangian(6), R_reduced(2, 2), Ht1(3), Ht2(3)
      real(wp) :: delta_r(3), alpha, directional_deriv, norm_g_tan
      real(wp) :: mu_min, mu_max, v_min(2), v_max(2), s_tan(2)
      integer :: iter, ls_iter
      logical :: converged, ls_success, at_saddle, unreached

      r_curr = r_init
      lambda_curr = lambda_init
      r_out = r_curr
      lambda_out = lambda_curr
      converged = .false.
      phi_curr = 0.5_wp*cfg%phi_alpha*sum((r_curr - anchor)**2)

      do iter = 0, cfg%riemann_max_iter
         call point_jet(lsf, r_curr, 2, S_val, unreached, eval_error, grad_S, hess_S)
         if (allocated(eval_error)) then
            call hand_over(eval_error, error, lsf_error)
            return
         end if
         if (unreached) then
            call fatal_error(error, "Riemannian escape: iterate beyond the level set reach")
            return
         end if
         grad_phi = cfg%phi_alpha*(r_curr - anchor)

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
         H_lagrangian(diag_slots) = cfg%phi_alpha - lambda_curr*hess_S(diag_slots)
         call packed_contract2(H_lagrangian, t1, Ht1, 3)
         call packed_contract2(H_lagrangian, t2, Ht2, 3)
         R_reduced(1, 1) = dot_product(t1, Ht1)
         R_reduced(1, 2) = dot_product(t1, Ht2)
         R_reduced(2, 2) = dot_product(t2, Ht2)
         R_reduced(2, 1) = R_reduced(1, 2)
         call eig_2x2_symmetric(R_reduced(1, 1), R_reduced(1, 2), R_reduced(2, 2), &
                                mu_min, mu_max, v_min, v_max)
         at_saddle = mu_min < -cfg%riemann_saddle_eig_tol

         if (norm_g_tan < tol_grad_tan .and. .not. at_saddle) then
            converged = .true.
            exit
         end if

         if (at_saddle) then
            s_tan = v_min
            delta_r = T_mat(:, 1)*s_tan(1) + T_mat(:, 2)*s_tan(2)
            directional_deriv = dot_product(grad_phi, delta_r)
            if (directional_deriv > 0.0_wp) then
               s_tan = -s_tan
               delta_r = -delta_r
               directional_deriv = -directional_deriv
            end if
            alpha = cfg%riemann_escape_step
         else
            if (mu_min > cfg%riemann_saddle_eig_tol) then
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
            alpha = cfg%riemann_alpha_init
         end if

         ls_success = .false.
         do ls_iter = 1, cfg%riemann_max_ls_iter
            r_trial = r_curr + alpha*delta_r
            call retract_to_surface(cfg, lsf, r_trial, anchor, r_retracted, lambda_trial, retract_error, &
                                    lsf_error=eval_error)
            if (allocated(eval_error)) then
               call hand_over(eval_error, error, lsf_error)
               return
            end if
            if (allocated(retract_error)) then
               deallocate (retract_error)
               alpha = alpha*cfg%riemann_alpha_reduce
               cycle
            end if

            phi_trial = 0.5_wp*cfg%phi_alpha*sum((r_retracted - anchor)**2)
            if (phi_trial < phi_curr - 1.0e-10_wp) then
               ls_success = .true.
               exit
            end if
            if (directional_deriv < 0.0_wp) then
               if (phi_trial <= phi_curr + cfg%riemann_armijo_c*alpha*directional_deriv) then
                  ls_success = .true.
                  exit
               end if
            end if
            alpha = alpha*cfg%riemann_alpha_reduce
         end do

         if (.not. ls_success) exit

         r_curr = r_retracted
         lambda_curr = lambda_trial
         phi_curr = phi_trial
         r_out = r_curr
         lambda_out = lambda_curr

         if (cfg%debug) then
            write (cfg%unit, "(i4,1x,es12.4,1x,f12.6,1x,es12.4,1x,es10.2,1x,es12.4,1x,l1)") &
               iter, phi_curr, norm2(r_curr - anchor), norm_g_tan, alpha, mu_min, at_saddle
         end if
      end do

      if (.not. converged .and. iter >= cfg%riemann_max_iter) then
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
   !> @param[inout] lsf        Thread-local level set
   !> @param[in]    r_init     Start point
   !> @param[in]    anchor     Anchor
   !> @param[out]   r_out      Point with `|S| < 1e-12`
   !> @param[out]   lambda_out Multiplier at `r_out`
   !> @param[out]   error      Degenerate gradient, unreached trial or iteration budget exhausted
   !> @param[out]   lsf_error  LSF evaluation failure; folded into `error` when absent
   subroutine retract_to_surface(cfg, lsf, r_init, anchor, r_out, lambda_out, error, lsf_error)
      type(polish_config_type), intent(in) :: cfg
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
            lambda_out = dot_product(cfg%phi_alpha*(r_curr - anchor), grad_S)/norm_grad_S2
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
   !> @param[out]   third     Third derivative tensor, required for order 3
   subroutine point_jet(lsf, point, order, S, unreached, error, grad, hess, third)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      real(wp), intent(in) :: point(3)
      integer, intent(in) :: order
      real(wp), intent(out) :: S
      logical, intent(out) :: unreached
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(out), optional :: grad(3)
      real(wp), intent(out), optional :: hess(:)
      real(wp), intent(out), optional :: third(3, 3, 3)

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
      if (order >= 3 .and. present(third)) call unpack_sym_tensor3(t3(:, 1), third, 3)
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
