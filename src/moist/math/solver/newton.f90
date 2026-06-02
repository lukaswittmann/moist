!> Newton-Raphson nonlinear equation solver interface for MOIST
!>
!> self module provides a clean interface to the nlesolver_module,
!> wrapping it with MOIST conventions and error handling.
module moist_math_solver_newton
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_type, only: solver_base_type
   use nlesolver_module
   implicit none
   private

   public :: moist_math_solver_newton_type
   public :: new_newton_solver  ! Factory function (unified interface)
   public :: newton_linear_solver_dense, newton_linear_solver_lsqr
   public :: newton_linear_solver_lusol, newton_linear_solver_lsmr
   public :: newton_linesearch_simple, newton_linesearch_backtrack
   public :: newton_linesearch_exact, newton_linesearch_fixedpoint
   public :: newton_bounds_ignore, newton_bounds_scalar
   public :: newton_bounds_vector, newton_bounds_wall

   !> Linear solver options
   integer, parameter :: newton_linear_solver_dense = NLESOLVER_SPARSITY_DENSE
   integer, parameter :: newton_linear_solver_lsqr = NLESOLVER_SPARSITY_LSQR
   integer, parameter :: newton_linear_solver_lusol = NLESOLVER_SPARSITY_LUSOL
   integer, parameter :: newton_linear_solver_lsmr = NLESOLVER_SPARSITY_LSMR

   !> Line search options
   integer, parameter :: newton_linesearch_simple = NLESOLVER_LINESEARCH_SIMPLE
   integer, parameter :: newton_linesearch_backtrack = NLESOLVER_LINESEARCH_BACKTRACKING
   integer, parameter :: newton_linesearch_exact = NLESOLVER_LINESEARCH_EXACT
   integer, parameter :: newton_linesearch_fixedpoint = NLESOLVER_LINESEARCH_FIXEDPOINT

   !> Bounds handling options
   integer, parameter :: newton_bounds_ignore = NLESOLVER_IGNORE_BOUNDS
   integer, parameter :: newton_bounds_scalar = NLESOLVER_SCALAR_BOUNDS
   integer, parameter :: newton_bounds_vector = NLESOLVER_VECTOR_BOUNDS
   integer, parameter :: newton_bounds_wall = NLESOLVER_WALL_BOUNDS

   !> User function interfaces (clean, without class argument)
   abstract interface
      !> Compute the function vector
      subroutine func_interface(x, f)
         import :: wp
         real(wp), dimension(:), intent(in) :: x  !> variables
         real(wp), dimension(:), intent(out) :: f !> function values
      end subroutine func_interface

      !> Compute the Jacobian matrix (dense)
      subroutine grad_interface(x, jac)
         import :: wp
         real(wp), dimension(:), intent(in) :: x      !> variables
         real(wp), dimension(:, :), intent(out) :: jac !> Jacobian matrix
      end subroutine grad_interface

      !> Compute the Jacobian matrix (sparse)
      subroutine grad_sparse_interface(x, jac_sparse)
         import :: wp
         real(wp), dimension(:), intent(in) :: x           !> variables
         real(wp), dimension(:), intent(out) :: jac_sparse !> sparse Jacobian
      end subroutine grad_sparse_interface

      !> Debug callback for iteration monitoring
      subroutine debug_callback_interface(iter, x, f, jac, jac_sparse)
         import :: wp
         integer, intent(in) :: iter                              !> iteration number
         real(wp), dimension(:), intent(in) :: x                  !> current variables
         real(wp), dimension(:), intent(in) :: f                  !> current residuals
         real(wp), dimension(:, :), intent(in), optional :: jac    !> Jacobian (dense)
         real(wp), dimension(:), intent(in), optional :: jac_sparse !> Jacobian (sparse)
      end subroutine debug_callback_interface

      !> User input check callback for early termination
      function user_input_check_interface(x) result(stop_solver)
         import :: wp
         real(wp), dimension(:), intent(in) :: x  !> current variables
         logical :: stop_solver                    !> true to stop solver
      end function user_input_check_interface
   end interface

   !> Context-aware user function interfaces (with context parameter for thread safety)
   abstract interface
      !> Compute the function vector (with context)
      subroutine func_context_interface(x, f, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x     !> variables
         real(wp), dimension(:), intent(out) :: f    !> function values
         class(*), intent(in) :: context             !> user context data
      end subroutine func_context_interface

      !> Compute the Jacobian matrix (dense, with context)
      subroutine grad_context_interface(x, jac, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x      !> variables
         real(wp), dimension(:, :), intent(out) :: jac !> Jacobian matrix
         class(*), intent(in) :: context              !> user context data
      end subroutine grad_context_interface

      !> Compute the Jacobian matrix (sparse, with context)
      subroutine grad_sparse_context_interface(x, jac_sparse, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x           !> variables
         real(wp), dimension(:), intent(out) :: jac_sparse !> sparse Jacobian
         class(*), intent(in) :: context                   !> user context data
      end subroutine grad_sparse_context_interface

      !> Debug callback for iteration monitoring (with context)
      subroutine debug_callback_context_interface(iter, x, f, context, jac, jac_sparse)
         import :: wp
         integer, intent(in) :: iter                              !> iteration number
         real(wp), dimension(:), intent(in) :: x                  !> current variables
         real(wp), dimension(:), intent(in) :: f                  !> current residuals
         class(*), intent(in) :: context                          !> user context data
         real(wp), dimension(:, :), intent(in), optional :: jac    !> Jacobian (dense)
         real(wp), dimension(:), intent(in), optional :: jac_sparse !> Jacobian (sparse)
      end subroutine debug_callback_context_interface

      !> User input check callback for early termination (with context)
      function user_input_check_context_interface(x, context) result(stop_solver)
         import :: wp
         real(wp), dimension(:), intent(in) :: x  !> current variables
         class(*), intent(in) :: context          !> user context data
         logical :: stop_solver                    !> true to stop solver
      end function user_input_check_context_interface
   end interface

   !> Newton-Raphson solver wrapper for MOIST
   type, extends(solver_base_type) :: moist_math_solver_newton_type
      private
      !> Underlying nlesolver instance
      class(nlesolver_type), allocatable :: solver

   contains
      !> Solve the nonlinear system
      procedure :: solve => newton_solve

      !> Get status information
      procedure :: status => newton_status

      !> Destroy the solver
      procedure :: destroy => newton_destroy
   end type moist_math_solver_newton_type

   !> Internal per-instance callback payload for Newton wrapper callbacks.
   type, extends(nlesolver_type) :: moist_newton_bridge_type
      private
      !> Legacy callbacks
      procedure(func_interface), pointer, nopass :: user_func => null()
      procedure(grad_interface), pointer, nopass :: user_grad => null()
      procedure(grad_sparse_interface), pointer, nopass :: user_grad_sparse => null()
      procedure(debug_callback_interface), pointer, nopass :: user_debug => null()
      procedure(user_input_check_interface), pointer, nopass :: user_check => null()

      !> Context-aware callbacks
      procedure(func_context_interface), pointer, nopass :: user_func_ctx => null()
      procedure(grad_context_interface), pointer, nopass :: user_grad_ctx => null()
      procedure(grad_sparse_context_interface), pointer, nopass :: user_grad_sparse_ctx => null()
      procedure(debug_callback_context_interface), pointer, nopass :: user_debug_ctx => null()
      procedure(user_input_check_context_interface), pointer, nopass :: user_check_ctx => null()

      !> Per-instance context and temporary callback state.
      class(*), allocatable :: user_context
      logical :: debug_mode = .false.
      integer :: debug_unit = -1
      real(wp), dimension(:), allocatable :: current_x
   end type moist_newton_bridge_type

contains

   !> Factory function to create and initialize a Newton solver
   !>
   !> This is a standalone constructor that allocates and initializes a Newton solver.
   !> Use this with polymorphic allocation:
   !>
   !>   class(solver_base_type), allocatable :: solver
   !>   call new_newton_solver(solver, n, m, func, grad, error, ...)
   !>
   !> Factory function to create and initialize a Newton solver (unified interface)
   !>
   !> This is a standalone constructor that allocates and initializes a Newton solver.
   !> Supports both legacy (no context) and context-aware (thread-safe) interfaces.
   !> Use this with polymorphic allocation:
   !>
   !>   ! Context-aware (thread-safe):
   !>   call new_newton_solver(solver, n, m, error, &
   !>      func_ctx=my_residual, grad_ctx=my_jacobian, context=my_context, ...)
   !>
   !>   ! Legacy:
   !>   call new_newton_solver(solver, n, m, error, &
   !>      func=my_residual, grad=my_jacobian, ...)
   !>
   subroutine new_newton_solver(solver, n, m, error, &
                                func, grad, grad_sparse, debug_callback, user_input_check, &
                                func_ctx, grad_ctx, grad_sparse_ctx, debug_callback_ctx, user_input_check_ctx, context, &
                                max_iter, tol, tolx, linear_solver, linesearch, use_broyden, broyden_update_n, &
                                alpha, alpha_min, alpha_max, bounds_mode, xlow, xupp, verbose, iunit, &
                                irow, icol, debug_unit)
      class(solver_base_type), allocatable, intent(out) :: solver
      integer, intent(in) :: n
      integer, intent(in) :: m
      type(error_type), allocatable, intent(out) :: error
      ! Legacy interface callbacks (all optional)
      procedure(func_interface), optional :: func
      procedure(grad_interface), optional :: grad
      procedure(grad_sparse_interface), optional :: grad_sparse
      procedure(debug_callback_interface), optional :: debug_callback
      procedure(user_input_check_interface), optional :: user_input_check
      ! Context-aware interface callbacks (all optional)
      procedure(func_context_interface), optional :: func_ctx
      procedure(grad_context_interface), optional :: grad_ctx
      procedure(grad_sparse_context_interface), optional :: grad_sparse_ctx
      procedure(debug_callback_context_interface), optional :: debug_callback_ctx
      procedure(user_input_check_context_interface), optional :: user_input_check_ctx
      class(*), intent(in), optional :: context
      ! Solver options
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol
      real(wp), intent(in), optional :: tolx
      integer, intent(in), optional :: linear_solver
      integer, intent(in), optional :: linesearch
      logical, intent(in), optional :: use_broyden
      integer, intent(in), optional :: broyden_update_n
      real(wp), intent(in), optional :: alpha
      real(wp), intent(in), optional :: alpha_min
      real(wp), intent(in), optional :: alpha_max
      integer, intent(in), optional :: bounds_mode
      real(wp), dimension(n), intent(in), optional :: xlow
      real(wp), dimension(n), intent(in), optional :: xupp
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: iunit
      integer, dimension(:), intent(in), optional :: irow
      integer, dimension(:), intent(in), optional :: icol
      integer, intent(in), optional :: debug_unit

      type(moist_math_solver_newton_type), allocatable, target :: tmp
      !> Local variables
      integer :: max_iter_use, linear_solver_use, linesearch_use
      integer :: broyden_update_n_use, bounds_mode_use
      real(wp) :: tol_use, tolx_use, alpha_use, alpha_min_use, alpha_max_use
      logical :: use_broyden_use, verbose_use, use_context_interface
      integer :: istat
      character(len=:), allocatable :: message

      ! Allocate concrete type
      allocate (tmp)
      allocate (moist_newton_bridge_type :: tmp%solver)

      ! Set defaults
      max_iter_use = 100
      if (present(max_iter)) max_iter_use = max_iter

      tol_use = 1.0e-8_wp
      if (present(tol)) tol_use = tol

      tolx_use = 1.0e-10_wp
      if (present(tolx)) tolx_use = tolx

      linear_solver_use = newton_linear_solver_dense
      if (present(linear_solver)) linear_solver_use = linear_solver

      linesearch_use = newton_linesearch_backtrack
      if (present(linesearch)) linesearch_use = linesearch

      use_broyden_use = .false.
      if (present(use_broyden)) use_broyden_use = use_broyden

      broyden_update_n_use = 4
      if (present(broyden_update_n)) broyden_update_n_use = broyden_update_n

      alpha_use = 1.0_wp
      if (present(alpha)) alpha_use = alpha

      alpha_min_use = 0.1_wp
      if (present(alpha_min)) alpha_min_use = alpha_min

      alpha_max_use = 1.0_wp
      if (present(alpha_max)) alpha_max_use = alpha_max

      bounds_mode_use = newton_bounds_ignore
      if (present(bounds_mode)) bounds_mode_use = bounds_mode

      verbose_use = .false.
      if (present(verbose)) verbose_use = verbose

      use_context_interface = present(func_ctx) .or. present(grad_ctx) .or. present(grad_sparse_ctx)

      ! Determine which interface to use (context-aware takes precedence)
      if (use_context_interface) then
         ! Context-aware interface
         if (.not. present(context)) then
            call fatal_error(error, "Context-aware callbacks require a context to be provided")
            return
         end if
         if (.not. present(func_ctx)) then
            call fatal_error(error, "func_ctx must be provided when using context-aware interface")
            return
         end if
      else
         ! Legacy interface (no context)
         if (.not. present(func)) then
            call fatal_error(error, "func must be provided when not using context-aware interface")
            return
         end if
      end if

      select type (bridge => tmp%solver)
      type is (moist_newton_bridge_type)
         if (use_context_interface) then
            bridge%user_func_ctx => func_ctx
            if (present(grad_ctx)) bridge%user_grad_ctx => grad_ctx
            if (present(grad_sparse_ctx)) bridge%user_grad_sparse_ctx => grad_sparse_ctx
            if (present(debug_callback_ctx)) then
               bridge%user_debug_ctx => debug_callback_ctx
               bridge%debug_mode = .true.
            end if
            if (present(user_input_check_ctx)) bridge%user_check_ctx => user_input_check_ctx
            allocate (bridge%user_context, source=context)
         else
            bridge%user_func => func
            if (present(grad)) bridge%user_grad => grad
            if (present(grad_sparse)) bridge%user_grad_sparse => grad_sparse
            if (present(debug_callback)) then
               bridge%user_debug => debug_callback
               bridge%debug_mode = .true.
            end if
            if (present(user_input_check)) bridge%user_check => user_input_check
         end if
         if (present(debug_unit)) bridge%debug_unit = debug_unit
      class default
         call fatal_error(error, "Newton internal error: failed to allocate callback bridge")
         return
      end select

      ! Initialize the underlying nlesolver

      call tmp%solver%initialize( &
         n=n, m=m, max_iter=max_iter_use, tol=tol_use, &
         func=wrapper_func_module, grad=wrapper_grad_module, grad_sparse=wrapper_grad_sparse_module, &
         step_mode=linesearch_use, alpha=alpha_use, &
         alpha_min=alpha_min_use, alpha_max=alpha_max_use, &
         tolx=tolx_use, use_broyden=use_broyden_use, &
         broyden_update_n=broyden_update_n_use, &
         verbose=verbose_use, iunit=iunit, &
         sparsity_mode=linear_solver_use, &
         irow=irow, icol=icol, &
         bounds_mode=bounds_mode_use, xlow=xlow, xupp=xupp, &
         export_iteration=wrapper_export_module, &
         user_input_check=wrapper_check_module)

      ! Check initialization status
      call tmp%solver%status(istat, message)
      if (istat /= 0) then
         call fatal_error(error, "Newton solver initialization failed: "//message)
         return
      end if

      ! Move to polymorphic output
      call move_alloc(tmp, solver)

   end subroutine new_newton_solver

   !> Solve the nonlinear system using Newton-Raphson iteration.
   !>
   !> Attempts to find a solution x* such that f(x*) ~= 0 using Newton's method
   !> with optional line search and Broyden updates. The iteration continues until
   !> convergence criteria are met or maximum iterations reached.
   !>
   !> Convergence is achieved when:
   !>  - ||f(x)|| < tol (function norm below tolerance), or
   !>  - ||x_new - x_old|| < tolx (change in x below tolerance)
   !>
   !> @param[inout] self   Solver instance (must be initialized)
   !> @param[inout] x      Initial guess on input [n], solution on output [n]
   !> @param[out]   error  Error object, allocated if solver fails or doesn't converge
   subroutine newton_solve(self, x, error)
      class(moist_math_solver_newton_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      !> Local variables
      integer :: istat
      character(len=:), allocatable :: message

      if (.not. allocated(self%solver)) then
         call fatal_error(error, "Newton solver is not initialized")
         return
      end if

      ! Call the underlying solver
      call self%solver%solve(x)

      ! Check status and convert to error_type
      call self%solver%status(istat, message)

      select case (istat)
      case (1)
         ! Success: Required accuracy achieved
         return
      case (2)
         ! Success: Solution cannot be improved (likely converged)
         return
      case (3)
         ! Warning: Maximum iterations reached
         call fatal_error(error, "Newton solver [code 3]: Maximum iterations reached")
      case (4)
         ! Stopped by user
         write (*, *) "Newton solver stopped by user."
         ! call fatal_error(error, "Newton solver [code 4]: Stopped by user")
      case (5)
         ! Too many uphill steps
         call fatal_error(error, "Newton solver [code 5]: Too many steps in uphill direction")
      case default
         if (istat < 0) then
            ! Error condition - decode the error code
            call fatal_error(error, decode_error_message(istat, message))
         else
            ! Unknown status
            call fatal_error(error, "Newton solver [code "//trim(int_to_str(istat))//"]: Unknown status")
         end if
      end select

   end subroutine newton_solve

   !> Get status information from the solver.
   !>
   !> Retrieves the current status code and descriptive message from the last
   !> solver operation. Useful for diagnosing convergence issues or checking
   !> solver state.
   !>
   !> Status codes:
   !>  - 0: Successfully initialized
   !>  - 1: Required accuracy achieved
   !>  - 2: Solution cannot be improved
   !>  - 3: Maximum iterations reached
   !>  - 4: Stopped by user
   !>  - 5: Too many uphill steps
   !>  - Negative: Error condition (see message)
   !>
   !> @param[inout] self     Solver instance
   !> @param[out]   istat    Status code (optional)
   !> @param[out]   message  Human-readable status message (optional)
   subroutine newton_status(self, istat, message)
      class(moist_math_solver_newton_type), intent(inout) :: self
      integer, intent(out), optional :: istat
      character(len=:), allocatable, intent(out), optional :: message

      if (allocated(self%solver)) then
         call self%solver%status(istat, message)
      else
         if (present(istat)) istat = -999
         if (present(message)) message = "Newton solver is not initialized"
      end if
   end subroutine newton_status

   !> Destroy the solver and release resources.
   !>
   !> Cleans up internal state, deallocates memory, and nullifies function pointers.
   !> Should be called when the solver is no longer needed.
   !>
   !> @param[inout] self  Solver instance to destroy
   subroutine newton_destroy(self)
      class(moist_math_solver_newton_type), intent(inout), target :: self

      if (allocated(self%solver)) then
         call self%solver%destroy()
         select type (bridge => self%solver)
         type is (moist_newton_bridge_type)
            bridge%user_func => null()
            bridge%user_grad => null()
            bridge%user_grad_sparse => null()
            bridge%user_debug => null()
            bridge%user_check => null()

            bridge%user_func_ctx => null()
            bridge%user_grad_ctx => null()
            bridge%user_grad_sparse_ctx => null()
            bridge%user_debug_ctx => null()
            bridge%user_check_ctx => null()

            if (allocated(bridge%user_context)) deallocate (bridge%user_context)
            if (allocated(bridge%current_x)) deallocate (bridge%current_x)
            bridge%debug_mode = .false.
            bridge%debug_unit = -1
         end select
         deallocate (self%solver)
      end if
   end subroutine newton_destroy

   !> Internal wrapper for user function callback.
   !>
   !> Adapts the clean user interface (without class argument) to the
   !> nlesolver_module interface and dispatches through the instance-owned bridge.
   !>
   !> @param[inout] me  nlesolver instance (unused, required by interface)
   !> @param[in]    x   Variables [n]
   !> @param[out]   f   Function values (residuals) [m]
   subroutine wrapper_func_module(me, x, f)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f

      f = 0.0_wp
      select type (bridge => me)
      type is (moist_newton_bridge_type)
         if (associated(bridge%user_func_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_func_ctx(x, f, bridge%user_context)
         else if (associated(bridge%user_func)) then
            call bridge%user_func(x, f)
         end if
      end select
   end subroutine wrapper_func_module

   !> Internal wrapper for user Jacobian callback (dense).
   !>
   !> Adapts the clean user interface to nlesolver_module interface for dense
   !> Jacobian computation. Supports both legacy and context-aware interfaces.
   !>
   !> @param[inout] me   nlesolver instance (unused, required by interface)
   !> @param[in]    x    Variables [n]
   !> @param[out]   jac  Dense Jacobian matrix df/dx [m,n]
   subroutine wrapper_grad_module(me, x, jac)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: jac

      jac = 0.0_wp
      select type (bridge => me)
      type is (moist_newton_bridge_type)
         if (associated(bridge%user_grad_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_grad_ctx(x, jac, bridge%user_context)
         else if (associated(bridge%user_grad)) then
            call bridge%user_grad(x, jac)
         end if
      end select
   end subroutine wrapper_grad_module

   !> Decode error code into a detailed error message.
   !>
   !> Provides human-readable descriptions for all Newton solver error codes.
   !>
   !> @param[in]  istat    Error status code
   !> @param[in]  message  Original message from nlesolver
   !> @return     Detailed error message with code and description
   function decode_error_message(istat, message) result(msg)
      integer, intent(in) :: istat
      character(len=*), intent(in) :: message
      character(len=:), allocatable :: msg
      character(len=:), allocatable :: code_desc

      select case (istat)
      case (-1)
         code_desc = "Invalid alpha (step length must be in (0,1])"
      case (-2)
         code_desc = "Invalid alpha_min (must be in (0,1])"
      case (-3)
         code_desc = "Invalid alpha_max (must be in (0,1])"
      case (-4)
         code_desc = "Alpha_min must be < alpha_max"
      case (-5)
         code_desc = "Invalid step_mode (linesearch method)"
      case (-6)
         code_desc = "Error solving linear system (singular or ill-conditioned Jacobian)"
      case (-7)
         code_desc = "More than 5 consecutive uphill steps"
      case (-8)
         code_desc = "Divide by zero in Broyden update"
      case (-9)
         code_desc = "Out of memory"
      case (-10)
         code_desc = "Function routine not associated"
      case (-11)
         code_desc = "Gradient routine not associated"
      case (-12)
         code_desc = "Backtracking linesearch parameter c must be in (0,1)"
      case (-13)
         code_desc = "Backtracking linesearch parameter tau must be in (0,1)"
      case (-14)
         code_desc = "Must specify grad_sparse, irow, and icol for sparse mode"
      case (-15)
         code_desc = "Sparse pattern arrays irow and icol must be same length"
      case (-16)
         code_desc = "Lower bounds xlow > upper bounds xupp"
      case (-17)
         code_desc = "Error adjusting search direction for bounds"
      case (-18)
         code_desc = "Invalid norm_mode"
      case (-999)
         code_desc = "Solver not initialized"
      case (-1004)
         code_desc = "LSQR/LSMR: System appears ill-conditioned"
      case (-1005)
         code_desc = "LSQR/LSMR: Iteration limit reached"
      case (-1006)
         code_desc = "Custom sparse solver not provided"
      case default
         code_desc = "Unknown error"
      end select

      msg = "Newton solver [code "//trim(int_to_str(istat))//"]: "// &
            trim(code_desc) ! // " | Details: " // trim(message)

   end function decode_error_message

   !> Convert integer to string (helper function).
   !>
   !> @param[in]  i  Integer to convert
   !> @return     String representation
   function int_to_str(i) result(s)
      integer, intent(in) :: i
      character(len=20) :: s
      write (s, '(I0)') i
   end function int_to_str

   !> Internal wrapper for user Jacobian callback (sparse).
   !>
   !> Adapts the clean user interface to nlesolver_module interface for sparse
   !> Jacobian computation. Elements
   !> correspond to (irow, icol) sparsity pattern provided during initialization.
   !> Supports both legacy and context-aware interfaces.
   !>
   !> @param[inout] me          nlesolver instance (unused, required by interface)
   !> @param[in]    x           Variables [n]
   !> @param[out]   jac_sparse  Sparse Jacobian nonzero elements [n_nonzeros]
   subroutine wrapper_grad_sparse_module(me, x, jac_sparse)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: jac_sparse

      jac_sparse = 0.0_wp
      select type (bridge => me)
      type is (moist_newton_bridge_type)
         if (associated(bridge%user_grad_sparse_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_grad_sparse_ctx(x, jac_sparse, bridge%user_context)
         else if (associated(bridge%user_grad_sparse)) then
            call bridge%user_grad_sparse(x, jac_sparse)
         end if
      end select
   end subroutine wrapper_grad_sparse_module

   !> Internal wrapper for debug/export callback.
   !>
   !> Called after each Newton iteration to export debug information.
   !> Adapts nlesolver export interface to clean user callback.
   !>
   !> @param[in] x     Current variables [n]
   !> @param[in] f     Current residuals [m]
   !> @param[in] iter  Iteration number
   subroutine wrapper_export_module(me, x, f, iter)
      use nlesolver_module, only: nlesolver_type
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(in) :: f
      integer, intent(in) :: iter

      real(wp), dimension(:, :), allocatable :: jac_dense
      integer :: n, m

      select type (bridge => me)
      type is (moist_newton_bridge_type)
         ! Store current x for user_input_check callback
         if (allocated(bridge%current_x)) then
            if (size(bridge%current_x) /= size(x)) deallocate (bridge%current_x)
         end if
         if (.not. allocated(bridge%current_x)) allocate (bridge%current_x(size(x)))
         bridge%current_x = x

         if (.not. bridge%debug_mode) return

         n = size(x)
         m = size(f)

         ! Recompute Jacobian for debug output if callback is set
         if (associated(bridge%user_debug_ctx) .and. allocated(bridge%user_context)) then
            ! Context-aware debug callback
            if (associated(bridge%user_grad_ctx)) then
               ! Dense mode
               allocate (jac_dense(m, n))
               call bridge%user_grad_ctx(x, jac_dense, bridge%user_context)
               call bridge%user_debug_ctx(iter, x, f, bridge%user_context, jac=jac_dense)
            else if (associated(bridge%user_grad_sparse_ctx)) then
               ! Sparse mode - skip Jacobian in export for now
               call bridge%user_debug_ctx(iter, x, f, bridge%user_context)
            else
               ! No Jacobian available, just pass x and f
               call bridge%user_debug_ctx(iter, x, f, bridge%user_context)
            end if
         else if (associated(bridge%user_debug)) then
            ! Legacy debug callback
            if (associated(bridge%user_grad)) then
               ! Dense mode
               allocate (jac_dense(m, n))
               call bridge%user_grad(x, jac_dense)
               call bridge%user_debug(iter, x, f, jac=jac_dense)
            else if (associated(bridge%user_grad_sparse)) then
               ! Sparse mode - Note: cannot determine sparse size from nlesolver_type (private)
               ! User must pass consistent size; we skip Jacobian in export for now
               call bridge%user_debug(iter, x, f)
            else
               ! No Jacobian available, just pass x and f
               call bridge%user_debug(iter, x, f)
            end if
         end if
      end select

   end subroutine wrapper_export_module

   !> Internal wrapper for user_input_check callback.
   !>
   !> Adapts the clean user interface to nlesolver_module interface.
   !> Supports both legacy and context-aware interfaces.
   subroutine wrapper_check_module(me, user_stop)
      class(nlesolver_type), intent(inout) :: me
      logical, intent(out) :: user_stop

      user_stop = .false.
      select type (bridge => me)
      type is (moist_newton_bridge_type)
         if (allocated(bridge%current_x)) then
            if (associated(bridge%user_check_ctx) .and. allocated(bridge%user_context)) then
               user_stop = bridge%user_check_ctx(bridge%current_x, bridge%user_context)
            else if (associated(bridge%user_check)) then
               user_stop = bridge%user_check(bridge%current_x)
            end if
         end if
      end select
   end subroutine wrapper_check_module

end module moist_math_solver_newton
