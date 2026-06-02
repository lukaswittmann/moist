!> MOIST wrapper for SLSQP constrained optimization solver
!>
!> Provides a clean interface to the SLSQP (Sequential Least Squares Programming)
!> solver for constrained nonlinear optimization problems:
!>
!>   minimize f(x)
!>   subject to:  g_i(x) = 0  (equality constraints, i=1..meq)
!>                g_i(x) >= 0 (inequality constraints, i=meq+1..m)
!>                xl <= x <= xu (bounds)
module moist_math_solver_slsqp
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_type, only: solver_base_type
   use slsqp_module, only: slsqp_solver
   use iso_fortran_env, only: output_unit

   implicit none
   private

   public :: moist_math_solver_slsqp_type
   public :: new_slsqp_solver  ! Factory function (unified interface)

   ! User function interfaces (legacy, no context)
   abstract interface
      !> Compute objective function
      subroutine objective_interface(x, f)
         import :: wp
         real(wp), dimension(:), intent(in) :: x   !> variables
         real(wp), intent(out) :: f                !> objective value
      end subroutine objective_interface

      !> Compute gradient of objective
      subroutine objective_grad_interface(x, df)
         import :: wp
         real(wp), dimension(:), intent(in) :: x   !> variables
         real(wp), dimension(:), intent(out) :: df !> gradient
      end subroutine objective_grad_interface

      !> Compute constraints
      subroutine constraints_interface(x, c)
         import :: wp
         real(wp), dimension(:), intent(in) :: x   !> variables
         real(wp), dimension(:), intent(out) :: c  !> constraint values
      end subroutine constraints_interface

      !> Compute constraint Jacobian
      subroutine constraints_grad_interface(x, dc)
         import :: wp
         real(wp), dimension(:), intent(in) :: x      !> variables
         real(wp), dimension(:, :), intent(out) :: dc  !> constraint Jacobian (m x n)
      end subroutine constraints_grad_interface

      !> Iteration callback for debugging
      subroutine iteration_callback_interface(iter, x, f, c)
         import :: wp
         integer, intent(in) :: iter                  !> iteration number
         real(wp), dimension(:), intent(in) :: x      !> current variables
         real(wp), intent(in) :: f                    !> objective value
         real(wp), dimension(:), intent(in) :: c      !> constraint values
      end subroutine iteration_callback_interface
   end interface

   ! Context-aware user function interfaces (new interface for thread safety)
   abstract interface
      !> Compute objective function (with context)
      subroutine objective_context_interface(x, f, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x   !> variables
         real(wp), intent(out) :: f                !> objective value
         class(*), intent(in) :: context           !> user context data
      end subroutine objective_context_interface

      !> Compute gradient of objective (with context)
      subroutine objective_grad_context_interface(x, df, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x   !> variables
         real(wp), dimension(:), intent(out) :: df !> gradient
         class(*), intent(in) :: context           !> user context data
      end subroutine objective_grad_context_interface

      !> Compute constraints (with context)
      subroutine constraints_context_interface(x, c, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x   !> variables
         real(wp), dimension(:), intent(out) :: c  !> constraint values
         class(*), intent(in) :: context           !> user context data
      end subroutine constraints_context_interface

      !> Compute constraint Jacobian (with context)
      subroutine constraints_grad_context_interface(x, dc, context)
         import :: wp
         real(wp), dimension(:), intent(in) :: x      !> variables
         real(wp), dimension(:, :), intent(out) :: dc  !> constraint Jacobian (m x n)
         class(*), intent(in) :: context              !> user context data
      end subroutine constraints_grad_context_interface

      !> Iteration callback for debugging (with context)
      subroutine iteration_callback_context_interface(iter, x, f, c, context)
         import :: wp
         integer, intent(in) :: iter                  !> iteration number
         real(wp), dimension(:), intent(in) :: x      !> current variables
         real(wp), intent(in) :: f                    !> objective value
         real(wp), dimension(:), intent(in) :: c      !> constraint values
         class(*), intent(in) :: context              !> user context data
      end subroutine iteration_callback_context_interface
   end interface

   !> SLSQP solver wrapper for MOIST
   type, extends(solver_base_type) :: moist_math_solver_slsqp_type
      private
      !> Underlying SLSQP solver instance
      class(slsqp_solver), allocatable :: solver

   contains
      !> Solve the optimization problem
      procedure :: solve => slsqp_solve

      !> Clean up resources
      procedure :: destroy => slsqp_destroy
   end type moist_math_solver_slsqp_type

   !> Internal per-instance callback payload for the SLSQP wrappers.
   type, extends(slsqp_solver) :: moist_slsqp_bridge_type
      private
      !> User-provided function pointers (legacy, no context)
      procedure(objective_interface), pointer, nopass :: user_obj => null()
      procedure(objective_grad_interface), pointer, nopass :: user_obj_grad => null()
      procedure(constraints_interface), pointer, nopass :: user_con => null()
      procedure(constraints_grad_interface), pointer, nopass :: user_con_grad => null()
      procedure(iteration_callback_interface), pointer, nopass :: user_iter_callback => null()

      !> Context-aware function pointers
      procedure(objective_context_interface), pointer, nopass :: user_obj_ctx => null()
      procedure(objective_grad_context_interface), pointer, nopass :: user_obj_grad_ctx => null()
      procedure(constraints_context_interface), pointer, nopass :: user_con_ctx => null()
      procedure(constraints_grad_context_interface), pointer, nopass :: user_con_grad_ctx => null()
      procedure(iteration_callback_context_interface), pointer, nopass :: user_iter_callback_ctx => null()

      !> User context copied into the solver instance.
      class(*), allocatable :: user_context
   end type moist_slsqp_bridge_type

contains

   !> Factory function to create and initialize an SLSQP solver (unified interface)
   !>
   !> This is a standalone constructor that allocates and initializes an SLSQP solver.
   !> Supports both legacy (no context) and context-aware (thread-safe) interfaces.
   !> Use this with polymorphic allocation:
   !>
   !>   ! Context-aware (thread-safe):
   !>   call new_slsqp_solver(solver, n, m, meq, error, &
   !>      obj_ctx=my_obj, obj_grad_ctx=my_grad, context=my_context, ...)
   !>
   !>   ! Legacy:
   !>   call new_slsqp_solver(solver, n, m, meq, error, &
   !>      obj=my_obj, obj_grad=my_grad, ...)
   !>
   subroutine new_slsqp_solver(solver, n, m, meq, error, &
                               obj, obj_grad, con, con_grad, iter_callback, &
                               obj_ctx, obj_grad_ctx, con_ctx, con_grad_ctx, iter_callback_ctx, context, &
                               max_iter, tol, xl, xu, verbose, iunit, toldx, toldf)
      class(solver_base_type), allocatable, intent(out) :: solver
      integer, intent(in) :: n
      integer, intent(in) :: m
      integer, intent(in) :: meq
      type(error_type), allocatable, intent(out) :: error
      ! Legacy interface callbacks (all optional)
      procedure(objective_interface), optional :: obj
      procedure(objective_grad_interface), optional :: obj_grad
      procedure(constraints_interface), optional :: con
      procedure(constraints_grad_interface), optional :: con_grad
      procedure(iteration_callback_interface), optional :: iter_callback
      ! Context-aware interface callbacks (all optional)
      procedure(objective_context_interface), optional :: obj_ctx
      procedure(objective_grad_context_interface), optional :: obj_grad_ctx
      procedure(constraints_context_interface), optional :: con_ctx
      procedure(constraints_grad_context_interface), optional :: con_grad_ctx
      procedure(iteration_callback_context_interface), optional :: iter_callback_ctx
      class(*), intent(in), optional :: context
      ! Solver options
      integer, intent(in), optional :: max_iter
      real(wp), intent(in), optional :: tol
      real(wp), dimension(n), intent(in), optional :: xl
      real(wp), dimension(n), intent(in), optional :: xu
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: iunit
      real(wp), intent(in), optional :: toldx
      real(wp), intent(in), optional :: toldf

      type(moist_math_solver_slsqp_type), allocatable, target :: tmp
      integer :: max_iter_use, iunit_use
      real(wp) :: tol_use, toldx_use, toldf_use
      logical :: verbose_use, status_ok, use_context_interface

      ! Allocate concrete type
      allocate (tmp)
      allocate (moist_slsqp_bridge_type :: tmp%solver)

      ! Set defaults
      max_iter_use = 100
      if (present(max_iter)) max_iter_use = max_iter

      tol_use = 1.0e-8_wp
      if (present(tol)) tol_use = tol

      ! Set stagnation detection tolerances (disabled by default in SLSQP)
      toldx_use = -1.0_wp  ! Disabled
      if (present(toldx)) toldx_use = toldx

      toldf_use = -1.0_wp  ! Disabled
      if (present(toldf)) toldf_use = toldf

      verbose_use = .false.
      if (present(verbose)) verbose_use = verbose

      iunit_use = output_unit
      if (present(iunit)) iunit_use = iunit
      if (.not. verbose_use) iunit_use = 0  ! 0 = no printing

      use_context_interface = present(obj_ctx) .or. present(obj_grad_ctx) .or. present(con_ctx)

      select type (bridge => tmp%solver)
      type is (moist_slsqp_bridge_type)
         ! Validate interface inputs first. Callback pointers are assigned only
         ! after initialize(), because initialize() calls destroy() internally.
         if (use_context_interface) then
            ! Context-aware interface
            if (.not. present(context)) then
               call fatal_error(error, "Context-aware callbacks require a context to be provided")
               return
            end if
            if (.not. present(obj_ctx)) then
               call fatal_error(error, "obj_ctx must be provided when using context-aware interface")
               return
            end if
            if (.not. present(obj_grad_ctx)) then
               call fatal_error(error, &
                                "obj_grad_ctx must be provided when using context-aware interface")
               return
            end if
         else
            ! Legacy interface (no context)
            if (.not. present(obj)) then
               call fatal_error(error, "obj must be provided when not using context-aware interface")
               return
            end if
            if (.not. present(obj_grad)) then
               call fatal_error(error, &
                                "obj_grad must be provided when not using context-aware interface")
               return
            end if
         end if
      class default
         call fatal_error(error, "SLSQP internal error: failed to allocate callback bridge")
         return
      end select

      ! Initialize SLSQP solver
      if (present(iter_callback) .or. present(iter_callback_ctx)) then
         call tmp%solver%initialize(n=n, m=m, meq=meq, max_iter=max_iter_use, &
                                    acc=tol_use, toldx=toldx_use, toldf=toldf_use, &
                                    f=wrapper_objective, g=wrapper_gradient, linesearch_mode=1, &
                                    xl=xl, xu=xu, status_ok=status_ok, iprint=iunit_use, report=wrapper_iteration)
      else
         call tmp%solver%initialize(n=n, m=m, meq=meq, max_iter=max_iter_use, &
                                    acc=tol_use, toldx=toldx_use, toldf=toldf_use, &
                                    f=wrapper_objective, g=wrapper_gradient, linesearch_mode=1, &
                                    xl=xl, xu=xu, status_ok=status_ok, iprint=iunit_use)
      end if

      if (.not. status_ok) then
         call fatal_error(error, "SLSQP initialization failed")
         return
      end if

      select type (bridge => tmp%solver)
      type is (moist_slsqp_bridge_type)
         if (use_context_interface) then
            bridge%user_obj_ctx => obj_ctx
            bridge%user_obj_grad_ctx => obj_grad_ctx
            if (present(con_ctx)) bridge%user_con_ctx => con_ctx
            if (present(con_grad_ctx)) bridge%user_con_grad_ctx => con_grad_ctx
            if (present(iter_callback_ctx)) bridge%user_iter_callback_ctx => iter_callback_ctx
            allocate (bridge%user_context, source=context)
         else
            bridge%user_obj => obj
            bridge%user_obj_grad => obj_grad
            if (present(con)) bridge%user_con => con
            if (present(con_grad)) bridge%user_con_grad => con_grad
            if (present(iter_callback)) bridge%user_iter_callback => iter_callback
         end if
      class default
         call fatal_error(error, "SLSQP internal error: callback bridge unavailable after initialize")
         return
      end select

      ! Move to polymorphic output
      call move_alloc(tmp, solver)

   end subroutine new_slsqp_solver

   !> Solve the optimization problem
   subroutine slsqp_solve(self, x, error)
      use slsqp_module, only: mode_to_status_message
      class(moist_math_solver_slsqp_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      integer :: istat
      character(len=:), allocatable :: status_msg

      if (.not. allocated(self%solver)) then
         call fatal_error(error, "SLSQP solver is not initialized")
         return
      end if

      ! Call SLSQP solver
      call self%solver%optimize(x, istat)

      ! Check status
      if (istat /= 0) then
         status_msg = mode_to_status_message(istat)
         call fatal_error(error, "SLSQP solver failed with status "//trim(int_to_str(istat))//": "//status_msg)
      end if

   end subroutine slsqp_solve

   !> Clean up solver resources
   subroutine slsqp_destroy(self)
      class(moist_math_solver_slsqp_type), intent(inout), target :: self

      if (allocated(self%solver)) then
         call self%solver%destroy()
         select type (bridge => self%solver)
         type is (moist_slsqp_bridge_type)
            bridge%user_obj => null()
            bridge%user_obj_grad => null()
            bridge%user_con => null()
            bridge%user_con_grad => null()
            bridge%user_iter_callback => null()

            bridge%user_obj_ctx => null()
            bridge%user_obj_grad_ctx => null()
            bridge%user_con_ctx => null()
            bridge%user_con_grad_ctx => null()
            bridge%user_iter_callback_ctx => null()

            if (allocated(bridge%user_context)) deallocate (bridge%user_context)
         end select
         deallocate (self%solver)
      end if

   end subroutine slsqp_destroy

   !> Wrapper for objective function (SLSQP interface)
   !> Supports both legacy and context-aware interfaces
   subroutine wrapper_objective(me, x, f, c)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      real(wp), dimension(:), intent(out) :: c

      f = 0.0_wp
      c = 0.0_wp

      select type (bridge => me)
      type is (moist_slsqp_bridge_type)
         ! Compute objective
         if (associated(bridge%user_obj_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_obj_ctx(x, f, bridge%user_context)
         else if (associated(bridge%user_obj)) then
            call bridge%user_obj(x, f)
         end if

         ! Compute constraints
         if (associated(bridge%user_con_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_con_ctx(x, c, bridge%user_context)
         else if (associated(bridge%user_con)) then
            call bridge%user_con(x, c)
         end if
      end select

   end subroutine wrapper_objective

   !> Wrapper for gradient computation (SLSQP interface)
   !> Supports both legacy and context-aware interfaces
   subroutine wrapper_gradient(me, x, df, dc)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df
      real(wp), dimension(:, :), intent(out) :: dc

      df = 0.0_wp
      dc = 0.0_wp

      select type (bridge => me)
      type is (moist_slsqp_bridge_type)
         ! Compute objective gradient
         if (associated(bridge%user_obj_grad_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_obj_grad_ctx(x, df, bridge%user_context)
         else if (associated(bridge%user_obj_grad)) then
            call bridge%user_obj_grad(x, df)
         end if

         ! Compute constraint Jacobian
         if (associated(bridge%user_con_grad_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_con_grad_ctx(x, dc, bridge%user_context)
         else if (associated(bridge%user_con_grad)) then
            call bridge%user_con_grad(x, dc)
         end if
      end select

   end subroutine wrapper_gradient

   !> Wrapper for iteration callback (SLSQP interface)
   !> Supports both legacy and context-aware interfaces
   subroutine wrapper_iteration(me, iter, x, f, c)
      class(slsqp_solver), intent(inout) :: me
      integer, intent(in) :: iter
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(in) :: f
      real(wp), dimension(:), intent(in) :: c

      select type (bridge => me)
      type is (moist_slsqp_bridge_type)
         ! Call user iteration callback if provided
         if (associated(bridge%user_iter_callback_ctx) .and. allocated(bridge%user_context)) then
            call bridge%user_iter_callback_ctx(iter, x, f, c, bridge%user_context)
         else if (associated(bridge%user_iter_callback)) then
            call bridge%user_iter_callback(iter, x, f, c)
         end if
      end select

   end subroutine wrapper_iteration

   !> Convert integer to string
   function int_to_str(i) result(s)
      integer, intent(in) :: i
      character(len=20) :: s
      write (s, '(I0)') i
   end function int_to_str

end module moist_math_solver_slsqp
