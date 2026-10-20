!> Exact or finite-difference Jacobian-vector products
submodule(moist_model_moz_solver_jfnk) moist_model_moz_solver_jfnk_jacobian
   use moist_model_moz_solver_jfnk_kernels, only: blocked_dot, field_square_sum
   implicit none(type, external)

contains

   module subroutine jfnk_apply_jacobian(self, residual_cb, gamma_curr, f_curr, direction, jv, &
                                         feasible, error)
      class(moist_model_moz_solver_jfnk_type), intent(inout) :: self
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      real(wp), intent(in) :: gamma_curr(:, :)
      real(wp), intent(in) :: f_curr(:, :)
      real(wp), intent(in) :: direction(:)
      real(wp), intent(out) :: jv(:)
      logical, intent(out) :: feasible
      type(error_type), allocatable, intent(out) :: error

      if (self%jvp_ready) then
         call residual_cb%jvp(direction, jv, feasible, error)
         self%njvp_total = self%njvp_total + 1
      else
         call self%workspace%ensure_fd(size(gamma_curr, 1), size(gamma_curr, 2))
         call jfnk_fd_jvp(residual_cb, gamma_curr, f_curr, direction, jv, feasible, error, &
                          self%params%fd_eps, self%workspace%ws_gamma_pert, self%workspace%ws_f_pert, &
                          self%nresidual_total, self%team_size())
      end if
      if (allocated(error)) feasible = .false.
   end subroutine jfnk_apply_jacobian

   !> Finite-difference Jacobian-vector product: jv ~= (F(x + sigma*d) - F(x))/sigma
   !>
   !> Perturbation magnitude follows the JFNK (NITSOL / Brown-Saad) scaling
   !> sigma = sqrt(eps_mach)*(1 + ||x||)/||d||, but with ||x|| taken as the
   !> root-mean-square norm sqrt(sum(x**2)/n) rather than the raw L2 norm
   !>
   !> - RMS scaling keeps ||sigma*d|| independent of the field size at fixed RMS state
   !> - Per-entry increments and accuracy still depend on the direction and residual noise
   !> - Exact products avoid this probe-size choice
   !>
   !> @param[in]  residual_cb Residual evaluator
   !> @param[in]  gamma_curr  current state x
   !> @param[in]  f_curr      residual F(x) at the current state
   !> @param[in]  direction   Krylov direction d (flat, un-normalized)
   !> @param[out] jv          approximation to J(x)*d
   !> @param[out] feasible    Whether the probe succeeded
   !> @param[out] error       Probe failure information
   !> @param[in]  fd_eps      manual step override; if > 0 the probe is a fixed
   !>                         step fd_eps along the unit direction
   !> @param[inout] gamma_pert Perturbed field scratch (npts, npair)
   !> @param[inout] f_pert     Perturbed residual scratch (npts, npair)
   !> @param[inout] nresidual residual-evaluation counter, incremented per probe
   !> @param[in] nt OpenMP team size
   subroutine jfnk_fd_jvp(residual_cb, gamma_curr, f_curr, direction, jv, &
                           feasible, error, fd_eps, gamma_pert, f_pert, nresidual, nt)
      !> Model residual evaluator
      class(moist_model_moz_solver_evaluator_type), intent(in) :: residual_cb
      !> Current point (npts, npair)
      real(wp), intent(in) :: gamma_curr(:, :)
      !> Residual at the current point (npts, npair)
      real(wp), intent(in) :: f_curr(:, :)
      !> Flat direction
      real(wp), intent(in) :: direction(:)
      !> Flat approximate product
      real(wp), intent(out) :: jv(:)
      !> Whether the probe succeeded
      logical, intent(out) :: feasible
      !> Probe failure information
      type(error_type), allocatable, intent(out) :: error
      !> Positive manual step or zero for automatic scaling
      real(wp), intent(in) :: fd_eps
      !> Perturbed field scratch
      real(wp), intent(inout) :: gamma_pert(:, :)
      !> Perturbed residual scratch
      real(wp), intent(inout) :: f_pert(:, :)
      !> Residual-evaluation counter
      integer, intent(inout) :: nresidual
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt

      integer :: npts, npair, i, j
      !> Perturbation scale applied to `direction`, L2 norm of `direction`
      real(wp) :: sigma, norm_d
      !> Root-mean-square norm of the current state vector
      real(wp) :: rms_gamma

      feasible = .true.
      npts = size(gamma_curr, 1)
      npair = size(gamma_curr, 2)

      norm_d = sqrt(blocked_dot(direction, direction, nt))
      if (.not. ieee_is_finite(norm_d)) then
         feasible = .false.
         call fatal_error(error, "Jacobian direction norm must be finite")
         return
      end if
      if (norm_d == 0.0_wp) then
         !$omp parallel do schedule(static) num_threads(merge(nt, 1, size(jv) >= parallel_vector_size))
         do i = 1, size(jv)
            jv(i) = 0.0_wp
         end do
         !$omp end parallel do
         return
      end if

      !> sigma scales the (un-normalized) direction, so
      !> jv = (F(x+sigma*d)-F(x))/sigma recovers J*direction
      !> A positive fd_eps overrides the auto scaling with a fixed step
      !> along the unit direction (sigma = fd_eps/||d||), so that
      !> jv = ||d||*(F(x + fd_eps*d/||d||) - F(x))/fd_eps as before
      if (fd_eps > 0.0_wp) then
         sigma = fd_eps/norm_d
      else
         !> Grid-size-independent RMS norm of the current state
         rms_gamma = sqrt(field_square_sum(gamma_curr, nt)/real(size(gamma_curr), wp))
         sigma = sqrt(epsilon(1.0_wp))*(1.0_wp + rms_gamma)/norm_d
      end if

      !$omp parallel do collapse(2) schedule(static) &
      !$omp& num_threads(merge(nt, 1, size(direction) >= parallel_vector_size))
      do j = 1, npair
         do i = 1, npts
            gamma_pert(i, j) = gamma_curr(i, j) + sigma*direction(i + (j - 1)*npts)
         end do
      end do
      !$omp end parallel do
      call residual_cb%residual(gamma_pert, f_pert, feasible, error)
      nresidual = nresidual + 1
      if (allocated(error)) feasible = .false.
      if (.not. feasible) return

      !$omp parallel do collapse(2) schedule(static) num_threads(merge(nt, 1, size(jv) >= parallel_vector_size))
      do j = 1, npair
         do i = 1, npts
            jv(i + (j - 1)*npts) = (f_pert(i, j) - f_curr(i, j))/sigma
         end do
      end do
      !$omp end parallel do
   end subroutine jfnk_fd_jvp

end submodule moist_model_moz_solver_jfnk_jacobian
