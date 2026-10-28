!> Newton outer step: forcing terms, inner solve and backtracking line search
submodule(moist_model_moz_solver_jfnk) moist_model_moz_solver_jfnk_step
   use moist_model_moz_solver_core_base, only: validate_solver_arrays
   use moist_model_moz_solver_jfnk_kernels, only: blocked_dot, copy_vector, field_square_sum, trial_field, unpack_vector
   implicit none(type, external)

contains

   module subroutine jfnk_outer_step(self, residual_fn, gamma, gamma_new, rmse, feasible, error)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_fn
      real(wp), intent(in) :: gamma(:, :)
      real(wp), intent(out) :: gamma_new(:, :)
      real(wp), intent(out) :: rmse
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error

      integer :: npts, npair, nflat, i, j
      real(wp) :: alpha, armijo_lhs, armijo_rhs, slope
      real(wp) :: phi0, phi_a, abs_slope, denom, alpha_new
      !> Number of trial steps taken by the backtracking line search
      integer :: ls_iter
      !> .true. once a trial step satisfies the Armijo sufficient-decrease test
      logical :: accepted, any_feasible
      !> Merit of the best extended step and feasibility of an extension trial
      real(wp) :: phi_best
      logical :: expand_feasible
      logical :: gmres_feasible
      type(error_type), allocatable :: gmres_error

      feasible = .false.
      call validate_solver_arrays(gamma, gamma_new, error)
      if (allocated(error)) return
      self%gmres_iters_last = 0
      self%jvp_ready = .false.
      call self%ensure_workspace(gamma)
      npts = size(gamma, 1)
      npair = size(gamma, 2)
      nflat = npts*npair

      ! Initialize Newton state on first call (sentinel: eta_current < 0)
      if (self%eta_current < 0.0_wp) then
         self%eta_current = self%params%eta0
         self%precond_valid_state = .false.
         self%rmse_prev_state = huge(1.0_wp)
         self%newton_iter_count = 0
         self%gmres_iters_total = 0
         self%nrecycled = 0
      end if
      self%newton_iter_count = self%newton_iter_count + 1

      ! Evaluate residual F(gamma)
      call residual_fn%residual(gamma, self%workspace%ws_f_curr, feasible, error)
      self%nresidual_total = self%nresidual_total + 1
      if (allocated(error)) feasible = .false.
      if (.not. feasible) return

      rmse = sqrt(field_square_sum(self%workspace%ws_f_curr, self%team_size())/real(nflat, wp))
      if (.not. ieee_is_finite(rmse)) then
         feasible = .false.
         call fatal_error(error, "Solver residual norm must be finite")
         return
      end if
      if (rmse == 0.0_wp) then
         gamma_new = gamma
         self%prev_step_maxabs = 0.0_wp
         self%prev_alpha = 0.0_wp
         return
      end if

      ! Linearize at gamma for exact Jacobian-vector products
      self%jvp_ready = .false.
      if (self%params%use_exact_jvp) then
         call self%ctx%timer%start("Linearization")
         if (residual_fn%has_jvp()) call residual_fn%linearize(gamma, self%jvp_ready, gmres_error)
         call self%ctx%timer%stop("Linearization")
         if (allocated(gmres_error)) then
            feasible = .false.
            call move_alloc(gmres_error, error)
            return
         end if
         if (.not. self%jvp_ready) then
            feasible = .false.
            call fatal_error(error, self%name()//" needs exact Jacobian-vector products at every iterate")
            return
         end if
      end if

      ! Coarse-space deflation of this step's linear solve, with the exact products
      ! validate requires; a failed or singular coarse matrix solves undeflated
      call self%deflation%deactivate()
      if (self%params%deflation .and. rmse < self%params%deflation_rmse) then
         call self%ctx%timer%start("Deflation setup")
         call residual_fn%deflation_setup(self%deflation, gmres_error)
         call self%ctx%timer%stop("Deflation setup")
         if (allocated(gmres_error)) then
            call self%deflation%deactivate()
            deallocate (gmres_error)
         end if
         if (.not. self%deflation%active) self%nundeflated_stage = self%nundeflated_stage + 1
      end if

      ! The operator decides from the previous step's residual whether it is used
      self%precond_valid_state = .false.
      if (associated(self%precond)) self%precond_valid_state = self%precond%ready(self%rmse_prev_state)
      call self%note_status()

      ! Set up linear solve
      !$omp parallel do collapse(2) schedule(static) &
      !$omp & num_threads(merge(self%team_size(), 1, nflat >= parallel_vector_size))
      do j = 1, npair
         do i = 1, npts
            self%workspace%ws_rhs_flat(i + (j - 1)*npts) = self%workspace%ws_f_curr(i, j)
         end do
      end do
      !$omp end parallel do

      call self%ctx%timer%start("Krylov solve")
      call self%linear_solve(residual_fn, gamma, self%workspace%ws_f_curr, &
                             self%workspace%ws_rhs_flat, self%workspace%ws_delta_flat, &
                             self%eta_current, self%precond_valid_state, &
                             gmres_feasible, gmres_error, &
                             gmres_iters=self%gmres_iters_last)
      call self%ctx%timer%stop("Krylov solve")
      self%gmres_iters_total = self%gmres_iters_total + self%gmres_iters_last
      self%nkrylov_total = self%nkrylov_total + self%gmres_iters_last
      if (allocated(gmres_error)) gmres_feasible = .false.
      if (.not. gmres_feasible) then
         feasible = .false.
         if (allocated(gmres_error)) call move_alloc(gmres_error, error)
         return
      end if

      call unpack_vector(self%workspace%ws_delta_flat, self%workspace%ws_delta, self%team_size())

      ! Compute true directional derivative for Armijo condition via one extra JVP
      ! An unusable product falls back to slope = -1, and a
      ! direction that is not a descent direction is replaced by the residual,
      ! searched with slope -||F||^2 before the Picard fallback below
      call self%apply_jacobian(residual_fn, gamma, self%workspace%ws_f_curr, self%workspace%ws_delta_flat, &
                               self%workspace%ws_w, gmres_feasible, gmres_error)
      if (allocated(gmres_error)) then
         gmres_feasible = .false.
         deallocate (gmres_error)
      end if
      slope = -1.0_wp
      if (gmres_feasible) slope = -blocked_dot(self%workspace%ws_rhs_flat, self%workspace%ws_w, self%team_size())
      if (.not. ieee_is_finite(slope)) slope = -1.0_wp
      if (slope >= 0.0_wp) then
         call copy_vector(self%workspace%ws_rhs_flat, self%workspace%ws_delta_flat, self%team_size())
         call unpack_vector(self%workspace%ws_delta_flat, self%workspace%ws_delta, self%team_size())
         slope = -blocked_dot(self%workspace%ws_rhs_flat, self%workspace%ws_rhs_flat, self%team_size())
      end if

      ! Backtracking line search with safeguarded quadratic interpolation of
      ! the merit function phi(a) = 0.5*||F(gamma - a*delta)||^2
      !
      ! - a rejected step's next trial minimises the quadratic model through
      !   phi(0), phi'(0) = slope and the latest phi(a), clamped to
      !   [0.1, tau]*a (Dennis-Schnabel safeguards)
      ! - the interpolation can shrink faster than fixed tau-halving
      ! - line_search_max_iter caps rejected backtracking trials, independently
      !   of the alpha floor; accepted full steps may be extended below
      phi0 = 0.5_wp*field_square_sum(self%workspace%ws_f_curr, self%team_size())
      abs_slope = abs(slope)
      armijo_rhs = self%params%line_search_c*abs_slope
      alpha = 1.0_wp
      accepted = .false.
      ls_iter = 0
      ! slope < 0 here unless the residual underflows to zero
      do while (slope < 0.0_wp)
         ls_iter = ls_iter + 1
         call trial_field(gamma, self%workspace%ws_delta, alpha, self%workspace%ws_gamma_trial, self%team_size())
         call residual_fn%residual(self%workspace%ws_gamma_trial, self%workspace%ws_f_trial, feasible, gmres_error)
         self%nresidual_total = self%nresidual_total + 1
         if (allocated(gmres_error)) feasible = .false.
         if (feasible) then
            phi_a = 0.5_wp*field_square_sum(self%workspace%ws_f_trial, self%team_size())
            feasible = ieee_is_finite(phi_a)
         end if
         if (feasible) then
            armijo_lhs = phi0 - phi_a
            if (armijo_lhs >= alpha*armijo_rhs) then
               ! Sufficient decrease achieved: keep this Newton step regardless of
               ! how small alpha is, so a valid short step is not thrown away in
               ! favour of the weaker Picard fallback below
               accepted = .true.
               exit
            end if
            ! Minimiser of the quadratic model through phi(0), phi'(0), phi(a)
            ! Armijo failure implies denom > 0, but guard against roundoff
            denom = 2.0_wp*(phi_a - phi0 + alpha*abs_slope)
            if (denom > 0.0_wp) then
               alpha_new = abs_slope*alpha*alpha/denom
               alpha_new = min(max(alpha_new, 0.1_wp*alpha), self%params%line_search_tau*alpha)
            else
               alpha_new = self%params%line_search_tau*alpha
            end if
            alpha = alpha_new
         else
            ! Infeasible trial: no merit value to interpolate; shrink by tau
            alpha = self%params%line_search_tau*alpha
         end if
         if (alpha < self%params%line_search_min_alpha) exit
         if (ls_iter >= self%params%line_search_max_iter) exit
      end do

      ! Extend a full step that still leaves much of the merit: double alpha up
      ! to alpha_max while each longer step is feasible, meets the Armijo test,
      ! and lowers the merit below the best step so far; the first worse trial
      ! ends the extension and the best step is kept
      if (accepted .and. ls_iter == 1 .and. self%params%alpha_max > 1.0_wp) then
         if (phi_a > self%params%line_search_expand_ratio*phi0) then
            phi_best = phi_a
            do while (alpha < self%params%alpha_max)
               alpha_new = min(2.0_wp*alpha, self%params%alpha_max)
               call trial_field(gamma, self%workspace%ws_delta, alpha_new, self%workspace%ws_gamma_trial, &
                  & self%team_size())
               call residual_fn%residual(self%workspace%ws_gamma_trial, self%workspace%ws_f_trial, expand_feasible, &
                  & gmres_error)
               self%nresidual_total = self%nresidual_total + 1
               if (allocated(gmres_error)) then
                  expand_feasible = .false.
                  deallocate (gmres_error)
               end if
               if (.not. expand_feasible) exit
               phi_a = 0.5_wp*field_square_sum(self%workspace%ws_f_trial, self%team_size())
               if (.not. ieee_is_finite(phi_a)) exit
               if (phi_a >= phi_best .or. phi0 - phi_a < alpha_new*armijo_rhs) exit
               alpha = alpha_new
               phi_best = phi_a
            end do
         end if
      end if

      ! If the GMRES line search never met the Armijo test (floor reached, trial
      ! budget exhausted, or every trial infeasible), fall back to a
      ! spectrally-filtered Picard step
      if (.not. accepted) then
         self%nfallback_stage = self%nfallback_stage + 1
         call copy_vector(self%workspace%ws_rhs_flat, self%workspace%ws_delta_flat, self%team_size())
         ! The lent fallback filter, not the active preconditioner
         if (associated(self%smoother)) call self%smoother%apply(self%workspace%ws_delta_flat)
         call unpack_vector(self%workspace%ws_delta_flat, self%workspace%ws_delta, self%team_size())
         ! Backtrack until a trial is both feasible and decreasing; take no step if none is
         ! Only an evaluated, decreasing trial is accepted
         ! `feasible` is restored afterwards because the trial evaluation
         ! overwrites this routine's intent(out) flag
         alpha = self%params%fallback_alpha
         any_feasible = .false.
         do
            call trial_field(gamma, self%workspace%ws_delta, alpha, self%workspace%ws_gamma_trial, self%team_size())
            call residual_fn%residual(self%workspace%ws_gamma_trial, self%workspace%ws_f_trial, feasible, gmres_error)
            self%nresidual_total = self%nresidual_total + 1
            if (allocated(gmres_error)) feasible = .false.
            if (feasible) then
               phi_a = 0.5_wp*field_square_sum(self%workspace%ws_f_trial, self%team_size())
               feasible = ieee_is_finite(phi_a)
            end if
            if (feasible) then
               any_feasible = .true.
               if (phi_a < phi0) then
                  accepted = .true.
                  exit
               end if
            end if
            alpha = alpha*self%params%line_search_tau
            if (alpha < self%params%line_search_min_alpha) exit
         end do
         if (.not. accepted) alpha = 0.0_wp
         ! Every trial non-finite: report infeasible (with the last error, if
         ! any) instead of pretending the unchanged iterate is fine
         feasible = any_feasible
         if (.not. any_feasible) then
            if (allocated(gmres_error)) call move_alloc(gmres_error, error)
            if (.not. allocated(error)) then
               call fatal_error(error, "Newton and Picard trial evaluations were infeasible")
            end if
            return
         end if
         ! A non-decreasing search returns a zero step
         ! The next update retries with a tightened forcing term
         ! At eta_min the first zero step is still returned, so a driver can accept
         ! an input that meets its tolerance
         ! A second one in a row would repeat the same solve and is reported as
         ! infeasible without an error, after a warning that names the stall
         ! The accepted-step diagnostics are kept as on any failure
         if (.not. accepted) self%nzero_stage = self%nzero_stage + 1
         if (.not. accepted .and. self%eta_current <= self%params%eta_min) then
            if (self%stalled_at_min) then
               feasible = .false.
               call self%ctx%warning(self%name()//": stalled at eta_min after two zero steps")
               return
            end if
            self%stalled_at_min = .true.
         end if
      end if
      if (accepted) self%stalled_at_min = .false.

      self%prev_step_maxabs = alpha*maxval(abs(self%workspace%ws_delta_flat))
      self%prev_alpha = alpha
      call trial_field(gamma, self%workspace%ws_delta, alpha, gamma_new, self%team_size())

      ! Forcing term for the next update: Eisenstat-Walker Choice 2 on this and
      ! the previous input residual, so it lags one update
      ! After a zero step the retry starts from the same iterate; tighten the
      ! tolerance so it gets a more accurate Newton direction
      if (alpha <= 0.0_wp) then
         self%eta_current = max(self%params%eta_min, 0.1_wp*self%eta_current)
      else if (self%rmse_prev_state < huge(1.0_wp)*0.5_wp) then
         block
            real(wp) :: eta_trial
            eta_trial = self%params%forcing_beta*(rmse/self%rmse_prev_state)**self%params%forcing_power
            self%eta_current = max(self%params%eta_min, min(self%params%eta_max, eta_trial))
         end block
      end if
      self%rmse_prev_state = rmse
   end subroutine jfnk_outer_step

end submodule moist_model_moz_solver_jfnk_step
