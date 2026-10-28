!> Solve a frozen transpose system with the shared JFNK/GMRES implementation
module moist_model_moz_solver_jfnk_adjoint
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type, resolve_num_threads
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_model_moz_solver_core_evaluator, only: moist_model_moz_solver_evaluator_type
   use moist_model_moz_solver_core_adjoint, only: moist_model_moz_solver_jfnk_adjoint_parameters_type
   use moist_model_moz_solver_core_base, only: moist_model_moz_solver_type
   use moist_model_moz_solver_jfnk, only: moist_model_moz_solver_jfnk_type, &
      & moist_model_moz_solver_jfnk_parameters_type, new_jfnk_solver
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   implicit none(type, external)
   private
   public :: solve_adjoint

   !> Linear residual J**T*x - rhs borrowing the frozen primal evaluator
   type, extends(moist_model_moz_solver_evaluator_type) :: transpose_evaluator_type
      !> Primal evaluator; must stay linearized throughout the solve
      class(moist_model_moz_solver_evaluator_type), pointer :: primal => null()
      !> Right-hand side with the primal field shape
      real(wp), allocatable :: rhs(:, :)
      !> OpenMP team size, ctx%get_num_threads() of the solve
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
   !> - Uses the shared GMRES through a linear JFNK evaluator
   !> - Tolerance is an absolute RMS of the unpreconditioned transpose residual
   !> - A supplied left preconditioner acts on the transpose system and stays fixed per inner solve
   !> - A guess, e.g. the solution of a nearby system, replaces the zero start
   !>
   !> @param[in]    ctx              run context; verbosity and OpenMP team size of the inner solver
   !> @param[in]    evaluator        frozen primal operator, borrowed for this call
   !> @param[in]    rhs              Euclidean right-hand side (npts, ncol)
   !> @param[out]   solution         solution (npts, ncol); allocated only on success
   !> @param[out]   error            invalid operator, controls or failure to converge
   !> @param[in]    tolerance        positive absolute RMS threshold; default 1e-10
   !> @param[in]    max_iter         maximum linear refinement steps; default 32
   !> @param[inout] preconditioner   borrowed left preconditioner for J**T
   !> @param[in]    guess            starting point (npts, ncol)
   !> @param[inout] workspace_owner  completed primal solver; its scratch is returned on success or failure,
   !>                                and its transpose controls apply when parameters is absent
   !> @param[in]    parameters       transpose Krylov restart and inner iterations; otherwise the owner's,
   !>                                else 50 and 50
   subroutine solve_adjoint(ctx, evaluator, rhs, solution, error, tolerance, max_iter, preconditioner, guess, &
      & workspace_owner, parameters)
      !> Run context; verbosity and OpenMP team size of the inner solver
      type(moist_context_type), intent(in), target :: ctx
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
      !> Optional starting point with the shape of rhs
      real(wp), intent(in), optional :: guess(:, :)
      !> Completed primal solver; its scratch is returned on success or failure,
      !> and its transpose controls apply when parameters is absent
      class(moist_model_moz_solver_type), optional, intent(inout) :: workspace_owner
      !> Transpose Krylov controls; otherwise the owner's, else defaults
      type(moist_model_moz_solver_jfnk_adjoint_parameters_type), optional, intent(in) :: parameters
      !> Frozen transpose evaluator
      type(transpose_evaluator_type) :: linear
      !> Effective transpose Krylov controls
      type(moist_model_moz_solver_jfnk_adjoint_parameters_type) :: controls
      !> Controls of the inner JFNK solver
      type(moist_model_moz_solver_jfnk_parameters_type) :: inner
      !> Inner JFNK solver
      type(moist_model_moz_solver_jfnk_type) :: solver
      !> Current solution, proposed solution and true residual, each (npts, ncol)
      real(wp), allocatable :: x(:, :), next(:, :), residual(:, :)
      !> Absolute RMS goal, current RMS residual and inner forcing term
      real(wp) :: tol, rms, forcing
      !> Maximum refinements
      integer :: limit
      !> Whether the last evaluation was usable
      logical :: feasible

      tol = 1.0e-10_wp
      limit = 32
      if (present(tolerance)) tol = tolerance
      if (present(max_iter)) limit = max_iter
      if (present(parameters)) then
         controls = parameters
      else if (present(workspace_owner)) then
         controls = workspace_owner%adjoint
      end if
      call controls%validate(error)
      if (allocated(error)) return
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
      linear%nthreads = ctx%get_num_threads()
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
      ! Intel does not like structure constructor with nested parameter defaults
      inner%restart = controls%restart
      inner%max_krylov = controls%max_krylov
      inner%eta0 = forcing
      inner%eta_min = forcing
      inner%eta_max = forcing
      inner%alpha_max = 1.0_wp
      inner%cgs2 = .true.
      call inner%validate(error)
      if (allocated(error)) return
      solver = new_jfnk_solver(ctx, inner)
      if (present(preconditioner)) call solver%set_preconditioner(preconditioner)
      if (present(workspace_owner)) call workspace_owner%exchange_workspace(solver%workspace)
      call solver%begin_stage()
      call adjoint_refine(solver, linear, x, next, residual, tol, limit, error)
      ! Every callback/convergence failure returns through this ownership boundary
      if (present(workspace_owner)) call workspace_owner%exchange_workspace(solver%workspace)
      if (allocated(error)) return
      call move_alloc(x, solution)
   end subroutine solve_adjoint

   !> Linear refinements and final true residual check, with one cleanup boundary
   !>
   !> @param[inout] solver    transpose solver with its own settings and borrowed scratch
   !> @param[in]    linear    frozen transpose evaluator
   !> @param[inout] x         current solution (npts, ncol)
   !> @param[out]   next      proposed solution (npts, ncol)
   !> @param[out]   residual  scratch for the true residual (npts, ncol)
   !> @param[in]    tol       absolute RMS residual goal
   !> @param[in]    limit     maximum refinements
   !> @param[out]   error     operator failure or unconverged transpose system
   subroutine adjoint_refine(solver, linear, x, next, residual, tol, limit, error)
      !> Transpose solver with its own settings and borrowed scratch
      class(moist_model_moz_solver_type), intent(inout) :: solver
      !> Frozen transpose evaluator
      type(transpose_evaluator_type), intent(in) :: linear
      !> Current solution
      real(wp), intent(inout) :: x(:, :)
      !> Proposed solution and scratch for the true residual
      real(wp), intent(out) :: next(:, :), residual(:, :)
      !> Absolute RMS residual goal
      real(wp), intent(in) :: tol
      !> Maximum refinements
      integer, intent(in) :: limit
      !> Operator failure or unconverged transpose system
      type(error_type), allocatable, intent(out) :: error
      !> Refinement counter
      integer :: iter
      !> RMS residual
      real(wp) :: rms
      !> Whether the last step or evaluation was usable
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
   !>
   !> @param[in]  self      linear evaluator
   !> @param[in]  gamma     current linear iterate (npts, ncol)
   !> @param[out] residual  J**T*gamma - rhs (npts, ncol)
   !> @param[out] feasible  product usability
   !> @param[out] error     transpose operator failure or shape mismatch
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
   !> @param[in]  self      linear evaluator
   !> @param[in]  n         number of field values
   !> @param[in]  gamma     current linear iterate, flat (n)
   !> @param[in]  rhs       right-hand side, flat (n)
   !> @param[out] residual  J**T*gamma - rhs, flat (n)
   !> @param[out] feasible  product usability
   !> @param[out] error     transpose operator failure
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
      !> Value index and OpenMP team size
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

   !> Whether the frozen transpose product is an exact linear JVP; always true
   !>
   !> @param[in] self  linear evaluator
   logical function transpose_has_jvp(self) result(available)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      available = .true.
   end function transpose_has_jvp

   !> Keep the primal linearization fixed while solving its transpose system
   !>
   !> @param[in]  self      linear evaluator
   !> @param[in]  gamma     linear iterate; does not change the operator
   !> @param[out] feasible  always true, operator is already frozen
   !> @param[out] error     unused, left unallocated
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
   !>
   !> @param[in]  self       linear evaluator
   !> @param[in]  direction  flat direction
   !> @param[out] jv         flat product J**T*direction
   !> @param[out] feasible   product usability
   !> @param[out] error      transpose operator failure
   subroutine transpose_jvp(self, direction, jv, feasible, error)
      !> Linear evaluator
      class(transpose_evaluator_type), intent(in) :: self
      !> Flat direction
      real(wp), intent(in) :: direction(:)
      !> Flat product J**T*direction
      real(wp), intent(out) :: jv(:)
      !> Product usability
      logical, intent(out) :: feasible
      !> Transpose operator failure
      type(error_type), allocatable, intent(out) :: error
      call self%primal%adjoint(direction, jv, feasible, error)
   end subroutine transpose_jvp

end module moist_model_moz_solver_jfnk_adjoint
