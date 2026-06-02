!> MOIST wrapper for L-BFGS-B bound-constrained optimization solver
!>
!> Provides a clean interface to the L-BFGS-B (Limited-memory Broyden-Fletcher-Goldfarb-Shanno
!> with Bounds) solver for large-scale bound-constrained optimization:
!>
!>   minimize f(x)
!>   subject to: l <= x <= u
!>
!> L-BFGS-B is particularly efficient for large problems where storing the full Hessian
!> is impractical. It uses a limited-memory BFGS update with reverse communication.
module moist_math_solver_lbfgsb
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_type, only: solver_base_type
   use lbfgsb_module, only: setulb, lbfgsp_wp
   use iso_fortran_env, only: output_unit

   implicit none
   private

   public :: moist_math_solver_lbfgsb_type
   public :: new_lbfgsb_solver  ! Factory function (unified interface)

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

      !> Iteration callback for debugging
      subroutine iteration_callback_interface(iter, x, f)
         import :: wp
         integer, intent(in) :: iter                  !> iteration number
         real(wp), dimension(:), intent(in) :: x      !> current variables
         real(wp), intent(in) :: f                    !> current objective value
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

      !> Iteration callback for debugging (with context)
      subroutine iteration_callback_context_interface(iter, x, f, context)
         import :: wp
         integer, intent(in) :: iter                  !> iteration number
         real(wp), dimension(:), intent(in) :: x      !> current variables
         real(wp), intent(in) :: f                    !> current objective value
         class(*), intent(in) :: context              !> user context data
      end subroutine iteration_callback_context_interface
   end interface

   !> L-BFGS-B solver wrapper for MOIST
   type, extends(solver_base_type) :: moist_math_solver_lbfgsb_type
      private
      !> Problem dimension
      integer :: n

      !> Number of limited-memory corrections (3 <= m <= 20 recommended)
      integer :: m

      !> Bounds type for each variable (n)
      !>  0 = unbounded
      !>  1 = lower bound only
      !>  2 = both bounds
      !>  3 = upper bound only
      integer, allocatable :: nbd(:)

      !> Lower bounds (n)
      real(wp), allocatable :: l(:)

      !> Upper bounds (n)
      real(wp), allocatable :: u(:)

      !> Current objective value
      real(wp) :: f

      !> Current gradient (n)
      real(wp), allocatable :: g(:)

      !> Tolerance parameters
      real(wp) :: factr  !> Relative function tolerance (scaled by machine epsilon)
      real(wp) :: pgtol  !> Projected gradient tolerance

      !> Maximum iterations
      integer :: max_iter

      !> Working arrays for L-BFGS-B
      real(wp), allocatable :: wa(:)   !> Real working array
      integer, allocatable :: iwa(:)   !> Integer working array

      !> L-BFGS-B state arrays
      character(len=60) :: task        !> Communication string
      character(len=60) :: csave       !> Character state
      logical :: lsave(4)              !> Logical state
      integer :: isave(44)             !> Integer state
      real(wp) :: dsave(29)            !> Real state

      !> Verbosity control
      integer :: iprint                !> Print control flag

      !> User-provided function pointers (legacy, no context)
      procedure(objective_interface), pointer, nopass :: user_obj => null()
      procedure(objective_grad_interface), pointer, nopass :: user_obj_grad => null()
      procedure(iteration_callback_interface), pointer, nopass :: user_iter_callback => null()

      !> Context-aware function pointers (new interface for thread safety)
      procedure(objective_context_interface), pointer, nopass :: user_obj_ctx => null()
      procedure(objective_grad_context_interface), pointer, nopass :: user_obj_grad_ctx => null()
      procedure(iteration_callback_context_interface), pointer, nopass :: &
         user_iter_callback_ctx => null()

      !> User context (stored as unlimited polymorphic for maximum flexibility)
      class(*), allocatable :: user_context

   contains
      !> Solve the optimization problem
      procedure :: solve => lbfgsb_solve

      !> Clean up resources
      procedure :: destroy => lbfgsb_destroy
   end type moist_math_solver_lbfgsb_type

contains

   !> Factory function to create and initialize an L-BFGS-B solver (unified interface)
   !>
   !> This is a standalone constructor that allocates and initializes an L-BFGS-B solver.
   !> Supports both legacy (no context) and context-aware (thread-safe) interfaces.
   !> Use this with polymorphic allocation:
   !>
   !>   ! Context-aware (thread-safe):
   !>   call new_lbfgsb_solver(solver, n, error, &
   !>      obj_ctx=my_obj, obj_grad_ctx=my_grad, context=my_context, ...)
   !>
   !>   ! Legacy:
   !>   call new_lbfgsb_solver(solver, n, error, &
   !>      obj=my_obj, obj_grad=my_grad, ...)
   !>
   !> @param[out]   solver             Polymorphic solver instance
   !> @param[in]    n                  Number of variables
   !> @param[out]   error              Error object, allocated on failure
   !> @param[in]    obj                Legacy objective function (optional)
   !> @param[in]    obj_grad           Legacy gradient function (optional)
   !> @param[in]    iter_callback      Legacy iteration callback (optional)
   !> @param[in]    obj_ctx            Context-aware objective (optional)
   !> @param[in]    obj_grad_ctx       Context-aware gradient (optional)
   !> @param[in]    iter_callback_ctx  Context-aware callback (optional)
   !> @param[in]    context            User context data (optional)
   !> @param[in]    m                  Memory parameter (default: 10)
   !> @param[in]    l                  Lower bounds [n] (optional)
   !> @param[in]    u                  Upper bounds [n] (optional)
   !> @param[in]    nbd                Bound type codes [n] (optional)
   !> @param[in]    factr              Function tolerance (default: 1.0e7)
   !> @param[in]    pgtol              Gradient tolerance (default: 1.0e-5)
   !> @param[in]    max_iter           Maximum iterations (default: 15000)
   !> @param[in]    verbose            Enable output (default: false)
   !> @param[in]    iunit              Output unit (default: output_unit)
   subroutine new_lbfgsb_solver(solver, n, error, &
                                obj, obj_grad, iter_callback, &
                                obj_ctx, obj_grad_ctx, iter_callback_ctx, context, &
                                m, l, u, nbd, factr, pgtol, max_iter, verbose, iunit)
      class(solver_base_type), allocatable, intent(out) :: solver
      integer, intent(in) :: n
      type(error_type), allocatable, intent(out) :: error
      ! Legacy interface callbacks (all optional)
      procedure(objective_interface), optional :: obj
      procedure(objective_grad_interface), optional :: obj_grad
      procedure(iteration_callback_interface), optional :: iter_callback
      ! Context-aware interface callbacks (all optional)
      procedure(objective_context_interface), optional :: obj_ctx
      procedure(objective_grad_context_interface), optional :: obj_grad_ctx
      procedure(iteration_callback_context_interface), optional :: iter_callback_ctx
      class(*), intent(in), optional :: context
      ! Solver options
      integer, intent(in), optional :: m
      real(wp), dimension(n), intent(in), optional :: l
      real(wp), dimension(n), intent(in), optional :: u
      integer, dimension(n), intent(in), optional :: nbd
      real(wp), intent(in), optional :: factr
      real(wp), intent(in), optional :: pgtol
      integer, intent(in), optional :: max_iter
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: iunit

      type(moist_math_solver_lbfgsb_type), allocatable, target :: tmp
      integer :: m_use, max_iter_use, iunit_use, iprint_use
      real(wp) :: factr_use, pgtol_use
      logical :: verbose_use
      integer :: wa_size, i

      ! Allocate concrete type
      allocate (tmp)

      ! Store dimension
      tmp%n = n

      ! Set defaults for parameters
      m_use = 10
      if (present(m)) m_use = m
      if (m_use < 3 .or. m_use > 20) then
         call fatal_error(error, "L-BFGS-B memory parameter m must be in range [3,20], got "// &
                          int_to_str(m_use))
         return
      end if
      tmp%m = m_use

      factr_use = 1.0e7_wp
      if (present(factr)) factr_use = factr
      tmp%factr = factr_use

      pgtol_use = 1.0e-5_wp
      if (present(pgtol)) pgtol_use = pgtol
      tmp%pgtol = pgtol_use

      max_iter_use = 15000
      if (present(max_iter)) max_iter_use = max_iter
      tmp%max_iter = max_iter_use

      verbose_use = .false.
      if (present(verbose)) verbose_use = verbose

      iunit_use = output_unit
      if (present(iunit)) iunit_use = iunit

      ! Set iprint based on verbosity
      if (verbose_use) then
         iprint_use = 1  ! Print summary every iteration
      else
         iprint_use = -1 ! No output
      end if
      tmp%iprint = iprint_use

      ! Determine which interface to use (context-aware takes precedence)
      if (present(obj_ctx) .or. present(obj_grad_ctx)) then
         ! Context-aware interface
         if (.not. present(obj_ctx)) then
            call fatal_error(error, "Context-aware interface requires obj_ctx")
            return
         end if
         if (.not. present(obj_grad_ctx)) then
            call fatal_error(error, "Context-aware interface requires obj_grad_ctx")
            return
         end if
         if (.not. present(context)) then
            call fatal_error(error, "Context-aware interface requires context")
            return
         end if

         tmp%user_obj_ctx => obj_ctx
         tmp%user_obj_grad_ctx => obj_grad_ctx
         if (present(iter_callback_ctx)) tmp%user_iter_callback_ctx => iter_callback_ctx
         allocate (tmp%user_context, source=context)

      else if (present(obj) .or. present(obj_grad)) then
         ! Legacy interface
         if (.not. present(obj)) then
            call fatal_error(error, "Legacy interface requires obj")
            return
         end if
         if (.not. present(obj_grad)) then
            call fatal_error(error, "Legacy interface requires obj_grad")
            return
         end if

         tmp%user_obj => obj
         tmp%user_obj_grad => obj_grad
         if (present(iter_callback)) tmp%user_iter_callback => iter_callback

      else
         call fatal_error(error, "Must provide either legacy (obj, obj_grad) or "// &
                          "context-aware (obj_ctx, obj_grad_ctx, context) callbacks")
         return
      end if

      ! Allocate and initialize bounds arrays
      allocate (tmp%nbd(n))
      allocate (tmp%l(n))
      allocate (tmp%u(n))

      if (present(nbd)) then
         tmp%nbd = nbd
      else
         ! Default: no bounds
         tmp%nbd = 0
      end if

      if (present(l)) then
         tmp%l = l
      else
         tmp%l = 0.0_wp  ! Unused if nbd=0 or 3
      end if

      if (present(u)) then
         tmp%u = u
      else
         tmp%u = 0.0_wp  ! Unused if nbd=0 or 1
      end if

      ! Validate bounds
      do i = 1, n
         if (tmp%nbd(i) < 0 .or. tmp%nbd(i) > 3) then
            call fatal_error(error, "Invalid nbd("//int_to_str(i)//") = "// &
                             int_to_str(tmp%nbd(i))//" (must be 0, 1, 2, or 3)")
            return
         end if
         if (tmp%nbd(i) == 2) then
            if (tmp%l(i) > tmp%u(i)) then
               call fatal_error(error, &
                                "Lower bound exceeds upper bound for variable "//int_to_str(i))
               return
            end if
         end if
      end do

      ! Allocate gradient
      allocate (tmp%g(n))

      ! Allocate working arrays
      ! wa size: (2*m*n + 5*n + 11*m^2 + 8*m)
      wa_size = 2*tmp%m*n + 5*n + 11*tmp%m**2 + 8*tmp%m
      allocate (tmp%wa(wa_size))
      allocate (tmp%iwa(3*n))

      ! Initialize task to START
      tmp%task = 'START'

      ! Move to polymorphic output
      call move_alloc(tmp, solver)

   end subroutine new_lbfgsb_solver

   !> Solve the optimization problem using L-BFGS-B
   !>
   !> Uses reverse communication to iteratively compute objective and gradient
   !> until convergence or failure.
   !>
   !> @param[inout] self   Solver instance (must be initialized)
   !> @param[inout] x      Initial guess on input [n], solution on output [n]
   !> @param[out]   error  Error object, allocated if solver fails
   subroutine lbfgsb_solve(self, x, error)
      class(moist_math_solver_lbfgsb_type), intent(inout), target :: self
      real(wp), dimension(:), intent(inout) :: x
      type(error_type), allocatable, intent(out) :: error

      integer :: iter
      logical :: converged

      ! Validate input
      if (size(x) /= self%n) then
         call fatal_error(error, "Input x size mismatch: expected "//int_to_str(self%n)// &
                          ", got "//int_to_str(size(x)))
         return
      end if

      ! Initialize iteration counter
      iter = 0
      converged = .false.

      ! Main reverse-communication loop
      do while (.not. converged)
         ! Call L-BFGS-B setulb
         call setulb(self%n, self%m, x, self%l, self%u, self%nbd, self%f, self%g, &
                     self%factr, self%pgtol, self%wa, self%iwa, self%task, &
                     self%iprint, self%csave, self%lsave, self%isave, self%dsave)

         ! Check task status
         if (self%task(1:2) == 'FG') then
            ! Compute function and gradient
            call compute_objective_and_gradient(self, x)

         else if (self%task(1:5) == 'NEW_X') then
            ! New iteration completed
            iter = self%isave(30)  ! Current iteration from isave(30)

            ! Call iteration callback if provided
            if (associated(self%user_iter_callback)) then
               call self%user_iter_callback(iter, x, self%f)
            else if (associated(self%user_iter_callback_ctx)) then
               call self%user_iter_callback_ctx(iter, x, self%f, self%user_context)
            end if

            ! Check iteration limit
            if (iter >= self%max_iter) then
               call fatal_error(error, "L-BFGS-B reached maximum iterations ("// &
                                int_to_str(self%max_iter)//")")
               return
            end if

         else if (self%task(1:4) == 'CONV') then
            ! Convergence achieved
            converged = .true.

         else if (self%task(1:5) == 'ERROR' .or. self%task(1:5) == 'ABNOR') then
            ! Error or abnormal termination
            call fatal_error(error, "L-BFGS-B error: "//trim(self%task))
            return

         else
            ! Unknown task
            call fatal_error(error, "L-BFGS-B unknown task: "//trim(self%task))
            return
         end if

      end do

   end subroutine lbfgsb_solve

   !> Clean up solver resources
   !>
   !> Deallocates all working arrays and nullifies function pointers.
   !>
   !> @param[inout] self  Solver instance to destroy
   subroutine lbfgsb_destroy(self)
      class(moist_math_solver_lbfgsb_type), intent(inout), target :: self

      ! Deallocate arrays
      if (allocated(self%nbd)) deallocate (self%nbd)
      if (allocated(self%l)) deallocate (self%l)
      if (allocated(self%u)) deallocate (self%u)
      if (allocated(self%g)) deallocate (self%g)
      if (allocated(self%wa)) deallocate (self%wa)
      if (allocated(self%iwa)) deallocate (self%iwa)

      ! Nullify function pointers (legacy)
      if (associated(self%user_obj)) nullify (self%user_obj)
      if (associated(self%user_obj_grad)) nullify (self%user_obj_grad)
      if (associated(self%user_iter_callback)) nullify (self%user_iter_callback)

      ! Nullify function pointers (context-aware)
      if (associated(self%user_obj_ctx)) nullify (self%user_obj_ctx)
      if (associated(self%user_obj_grad_ctx)) nullify (self%user_obj_grad_ctx)
      if (associated(self%user_iter_callback_ctx)) nullify (self%user_iter_callback_ctx)

      ! Deallocate context
      if (allocated(self%user_context)) deallocate (self%user_context)

   end subroutine lbfgsb_destroy

   !> Internal helper: compute objective and gradient
   !>
   !> Dispatches to user-provided callbacks (legacy or context-aware).
   !>
   !> @param[inout] self  Solver instance
   !> @param[in]    x     Current variables [n]
   subroutine compute_objective_and_gradient(self, x)
      type(moist_math_solver_lbfgsb_type), intent(inout) :: self
      real(wp), dimension(:), intent(in) :: x

      ! Use context-aware interface if available
      if (associated(self%user_obj_ctx) .and. associated(self%user_obj_grad_ctx)) then
         call self%user_obj_ctx(x, self%f, self%user_context)
         call self%user_obj_grad_ctx(x, self%g, self%user_context)
      else if (associated(self%user_obj) .and. associated(self%user_obj_grad)) then
         ! Use legacy interface
         call self%user_obj(x, self%f)
         call self%user_obj_grad(x, self%g)
      else
         ! Should never happen if initialize() validated correctly
         error stop "L-BFGS-B: No valid callback functions"
      end if

   end subroutine compute_objective_and_gradient

   !> Convert integer to string (helper function)
   !>
   !> @param[in]  i  Integer to convert
   !> @return     String representation
   function int_to_str(i) result(s)
      integer, intent(in) :: i
      character(len=20) :: s
      write (s, '(I0)') i
   end function int_to_str

end module moist_math_solver_lbfgsb
