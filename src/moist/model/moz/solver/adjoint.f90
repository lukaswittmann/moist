!> Solve a frozen transpose system with the shared JFNK/GMRES implementation
module moist_model_moz_solver_adjoint
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_model_moz_solver_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_base, only: moist_model_moz_solver_type
   use moist_model_moz_solver_factory, only: new_solver
   use moist_model_moz_solver_parameters, only: moist_model_moz_solver_parameters_type
   use moist_model_moz_solver_preconditioner, only: moist_model_moz_solver_preconditioner_type
   implicit none(type, external)
   private
   public :: solve_adjoint

   !> Linear residual J**T*x - rhs borrowing the frozen primal evaluator
   type, extends(moist_model_moz_solver_evaluator_type) :: transpose_evaluator_type
      !> Primal evaluator; must stay linearized throughout the solve
      class(moist_model_moz_solver_evaluator_type), pointer :: primal => null()
      !> Right-hand side with the primal field shape
      real(wp), allocatable :: rhs(:, :)
      !> OpenMP team size; zero takes omp_get_max_threads
      integer :: nthreads = 0
   contains
      procedure :: residual => transpose_residual
      procedure :: has_jvp => transpose_has_jvp
      procedure :: linearize => transpose_linearize
      procedure :: jvp => transpose_jvp
   end type transpose_evaluator_type

contains

   !> Solve J**T*x = rhs without assembling J or changing the primal linearization
   !>
   !> Uses the shared GMRES through a linear JFNK evaluator. The tolerance is an
   !> absolute RMS of the unpreconditioned transpose residual. A supplied left
   !> preconditioner must act on the transpose system and remain fixed per inner solve.
   !> A guess, e.g. the solution of a nearby system, replaces the zero start.
   subroutine solve_adjoint(evaluator, rhs, solution, error, tolerance, max_iter, preconditioner, nthreads, guess, &
      & workspace_owner, restart, max_krylov)
      !> Frozen primal operator, borrowed for this call
      class(moist_model_moz_solver_evaluator_type), target, intent(in) :: evaluator
      !> Euclidean RHS (npts, ncol)
      real(wp), intent(in) :: rhs(:, :)
      !> Solution; allocated only on success
      real(wp), allocatable, intent(out) :: solution(:, :)
      !> Invalid operator, controls or failure to converge
      type(error_type), allocatable, intent(out) :: error
      !> Positive absolute RMS threshold; default 1e-10
      real(wp), optional, intent(in) :: tolerance
      !> Maximum linear refinement steps; default 32
      integer, optional, intent(in) :: max_iter
      !> Optional borrowed left preconditioner for J**T
      class(moist_model_moz_solver_preconditioner_type), target, optional, intent(inout) :: preconditioner
      !> OpenMP team size, e.g. ctx%get_num_threads(); absent or zero takes omp_get_max_threads
      integer, intent(in), optional :: nthreads
      !> Optional starting point with the shape of rhs
      real(wp), intent(in), optional :: guess(:, :)
      !> Completed primal solver; its scratch is returned on success or failure
      class(moist_model_moz_solver_type), optional, intent(inout) :: workspace_owner
      !> Krylov basis restart length; default 50
      integer, optional, intent(in) :: restart
      !> Inner iterations per linear refinement; default 50
      integer, optional, intent(in) :: max_krylov
      type(transpose_evaluator_type) :: linear
      type(moist_model_moz_solver_parameters_type) :: parameters
      class(moist_model_moz_solver_type), allocatable :: solver
      real(wp), allocatable :: x(:, :), next(:, :), residual(:, :)
      real(wp) :: tol, rms, forcing
      integer :: limit
      logical :: feasible

      tol = 1.0e-10_wp
      limit = 32
      if (present(tolerance)) tol = tolerance
      if (present(max_iter)) limit = max_iter
      if (present(restart)) parameters%jfnk_restart = restart
      if (present(max_krylov)) parameters%jfnk_max_krylov = max_krylov
      if (parameters%jfnk_restart < 1 .or. parameters%jfnk_max_krylov < 1) then
         call fatal_error(error, "Adjoint needs positive Krylov limits")
         return
      end if
      if (.not. evaluator%has_adjoint()) then
         call fatal_error(error, "Adjoint solve requires an exact transpose operator")
         return
      end if
      if (size(rhs) == 0 .or. .not. all(ieee_is_finite(rhs)) .or. .not. ieee_is_finite(tol)) then
         call fatal_error(error, "Adjoint RHS must be nonempty and finite, with finite tolerance")
         return
      end if
      if (tol <= 0.0_wp .or. limit < 1) then
         call fatal_error(error, "Adjoint tolerance and refinement limit must be positive")
         return
      end if
      if (present(guess)) then
         if (any(shape(guess) /= shape(rhs)) .or. .not. all(ieee_is_finite(guess))) then
            call fatal_error(error, "Adjoint guess must be finite with the shape of the RHS")
            return
         end if
      end if
      linear%primal => evaluator
      if (present(nthreads)) then
         if (nthreads < 0) then
            call fatal_error(error, "Adjoint thread count must be nonnegative")
            return
         end if
         linear%nthreads = nthreads
      end if
      linear%rhs = rhs
      allocate (x, mold=rhs)
      allocate (next, residual, mold=rhs)
      x = 0.0_wp
      rms = sqrt(sum(rhs**2)/real(size(rhs), wp))
      if (.not. ieee_is_finite(rms)) then
         call fatal_error(error, "Adjoint RHS norm is non-finite")
         return
      end if
      if (present(guess)) then
         x = guess
         call linear%residual(x, residual, feasible, error)
         if (allocated(error)) return
         if (.not. feasible) then
            call fatal_error(error, "Adjoint residual of the guess is not usable")
            return
         end if
         rms = sqrt(sum(residual**2)/real(size(residual), wp))
      end if
      if (rms <= tol) then
         call move_alloc(x, solution)
         return
      end if
      forcing = max(1.0e-12_wp, min(0.1_wp, 0.1_wp*tol/rms))
      parameters%jfnk_eta0 = forcing
      parameters%jfnk_eta_min = forcing
      parameters%jfnk_eta_max = forcing
      parameters%jfnk_alpha_max = 1.0_wp
      parameters%jfnk_cgs2 = .true.
      parameters%nthreads = linear%nthreads
      call new_solver(solver, "JFNK", error, parameters)
      if (allocated(error)) return
      if (present(preconditioner)) call solver%set_preconditioner(preconditioner)
      if (present(workspace_owner)) call workspace_owner%exchange_workspace(solver%workspace)
      call solver%begin_stage()
      call adjoint_refine(solver, linear, x, next, residual, tol, limit, error)
      ! Every callback/convergence failure returns through this ownership boundary.
      if (present(workspace_owner)) call workspace_owner%exchange_workspace(solver%workspace)
      if (allocated(error)) return
      call move_alloc(x, solution)
   end subroutine solve_adjoint

   !> Linear refinements and final true residual check, with one cleanup boundary
   !>
   !> @param[inout] solver Transpose solver with its own settings and borrowed scratch
   !> @param[in] linear Frozen transpose evaluator
   !> @param[inout] x Current solution
   !> @param[out] next Proposed solution
   !> @param[out] residual Scratch for the true residual
   !> @param[in] tol Absolute RMS residual goal
   !> @param[in] limit Maximum refinements
   !> @param[out] error Operator failure or unconverged transpose system
   subroutine adjoint_refine(solver, linear, x, next, residual, tol, limit, error)
      class(moist_model_moz_solver_type), intent(inout) :: solver
      type(transpose_evaluator_type), intent(in) :: linear
      real(wp), intent(inout) :: x(:, :)
      real(wp), intent(out) :: next(:, :), residual(:, :)
      real(wp), intent(in) :: tol
      integer, intent(in) :: limit
      type(error_type), allocatable, intent(out) :: error
      integer :: iter
      real(wp) :: rms
      logical :: feasible

      do iter = 1, limit
         call solver%outer_step(linear, x, next, rms, feasible, error)
         if (allocated(error)) return
         if (.not. feasible) exit
         if (rms <= tol) exit
         x = next
      end do
      call linear%residual(x, residual, feasible, error)
      if (allocated(error)) return
      rms = sqrt(sum(residual**2)/real(size(residual), wp))
      if (.not. feasible .or. rms > tol) then
         call fatal_error(error, "Transpose system did not reach the requested adjoint tolerance")
         return
      end if
   end subroutine adjoint_refine

   !> Evaluate the frozen linear transpose residual
   subroutine transpose_residual(self, gamma, residual, feasible, error)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      !> Current linear iterate
      real(wp), intent(in) :: gamma(:, :)
      !> J**T*gamma - rhs
      real(wp), intent(out) :: residual(:, :)
      !> Product usability
      logical, intent(out) :: feasible
      !> Transpose operator failure
      type(error_type), allocatable, intent(out) :: error
      feasible = .false.
      if (any(shape(gamma) /= shape(self%rhs)) .or. any(shape(residual) /= shape(gamma))) then
         call fatal_error(error, "Adjoint field shapes do not match the RHS")
         return
      end if
      call flat_transpose_residual(self, size(gamma), gamma, self%rhs, residual, feasible, error)
   end subroutine transpose_residual

   !> Transpose residual on flat storage, so no field-sized temporary is needed
   !>
   !> @param[in] self Linear evaluator
   !> @param[in] n Number of field values
   !> @param[in] gamma Current linear iterate, flat (n)
   !> @param[in] rhs Right-hand side, flat (n)
   !> @param[out] residual J**T*gamma - rhs, flat (n)
   !> @param[out] feasible Product usability
   !> @param[out] error Transpose operator failure
   subroutine flat_transpose_residual(self, n, gamma, rhs, residual, feasible, error)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      !> Number of field values
      integer, intent(in) :: n
      !> Current linear iterate and right-hand side, sequence-associated with (npts, ncol)
      real(wp), intent(in) :: gamma(n), rhs(n)
      !> J**T*gamma - rhs
      real(wp), intent(out) :: residual(n)
      !> Product usability
      logical, intent(out) :: feasible
      !> Transpose operator failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i, nt

      call self%primal%adjoint(gamma, residual, feasible, error)
      if (allocated(error) .or. .not. feasible) return
      nt = merge(resolve_num_threads(self%nthreads), 1, n >= 65536)
      !$omp parallel do schedule(static) num_threads(nt)
      do i = 1, n
         residual(i) = residual(i) - rhs(i)
      end do
      !$omp end parallel do
   end subroutine flat_transpose_residual

   !> The frozen transpose product is an exact linear JVP
   logical function transpose_has_jvp(self) result(available)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      available = .true.
   end function transpose_has_jvp

   !> Keep the primal linearization fixed while solving its transpose system
   subroutine transpose_linearize(self, gamma, feasible, error)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      !> Linear iterate; does not change the operator
      real(wp), intent(in) :: gamma(:, :)
      !> Operator is already frozen
      logical, intent(out) :: feasible
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      feasible = .true.
   end subroutine transpose_linearize

   !> Forward an exact product through the frozen primal transpose
   subroutine transpose_jvp(self, direction, jv, feasible, error)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      !> Flat direction and output product
      real(wp), intent(in) :: direction(:)
      real(wp), intent(out) :: jv(:)
      !> Product usability
      logical, intent(out) :: feasible
      !> Transpose operator failure
      type(error_type), allocatable, intent(out) :: error
      call self%primal%adjoint(direction, jv, feasible, error)
   end subroutine transpose_jvp

end module moist_model_moz_solver_adjoint
