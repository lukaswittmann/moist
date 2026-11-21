!> Comprehensive tests for mathematical solvers with cross-solver comparison
!>
!> This module tests the optimization solvers (Newton, L-BFGS-B, fmin) on the same
!> problems using different formulations:
!>  - Newton: Solves nonlinear systems f(x)=0, handles constraints via Lagrangian
!>  - L-BFGS-B: Minimizes objectives with bound constraints only
!>
!> Test problems:
!>  1. Rosenbrock: Unconstrained optimization (equivalent formulations)
!>  2. Constrained circle: Minimize distance with constraint (Lagrangian/penalty)
module test_math_solvers
   use mctc_env, only: wp
   use mctc_io_utils, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_solver_type, only: solver_base_type

   use moist_math_solver_newton, only: new_newton_solver, newton_linesearch_backtrack, newton_linesearch_exact, &
                                       newton_linear_solver_custom_step
   use moist_math_solver_lbfgsb, only: new_lbfgsb_solver
   use mctc_env_error, only: moist_error_type => error_type
   ! Raw vendored solver APIs (upstream test ports folded into this suite)
   use lbfgsb_module, only: setulb
   use moist_math_solver_fmin, only: fmin
   use nlesolver_module, only: nlesolver_type, &
                               NLESOLVER_SCALAR_BOUNDS, &
                               NLESOLVER_SPARSITY_LSQR, NLESOLVER_SPARSITY_LUSOL, NLESOLVER_SPARSITY_LSMR
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use test_helpers, only: fd4_scalar, fd4_offsets
   implicit none(type, external)
   private

   public :: collect_math_solvers

   !> Tolerance for Newton solver
   real(wp), parameter :: newton_thr = 1.0e-10_wp
   !> Tolerance for L-BFGS-B solver
   real(wp), parameter :: lbfgsb_thr = 1.0e-6_wp

   !* ================================================================================= *!
   !*     Raw-kernel test constants (vendored APIs ported from upstream test suites)    *!
   !* ================================================================================= *!

   !> L-BFGS-B raw-kernel problem size and convergence controls
   integer, parameter :: lbfgsb_n = 25, lbfgsb_m = 5
   real(wp), parameter :: lbfgsb_factr = 1.0e7_wp, lbfgsb_pgtol = 1.0e-5_wp
   !> fmin raw-kernel tolerances (upstream accepts 10*tol on the minimizer)
   real(wp), parameter :: fmin_tol = 1.0e-8_wp
   real(wp), parameter :: fmin_thr = 10.0_wp*fmin_tol
   !> Newton/nlesolver raw-kernel problem definition
   integer, parameter :: newton_n = 2, newton_m = 2, newton_max_iter = 100
   real(wp), parameter :: newton_tol = 1.0e-8_wp
   !> Root of the test system (positive branch)
   real(wp), parameter :: newton_x_ref(2) = [0.5477225575051661_wp, -0.2_wp]
   !> Acceptance tolerance on the converged Newton solution
   real(wp), parameter :: newton_x_thr = 1.0e-4_wp
   !> COO sparsity pattern of the Jacobian (structural zero at (2,1) omitted)
   integer, parameter :: newton_irow(3) = [1, 1, 2]
   integer, parameter :: newton_icol(3) = [1, 2, 2]

contains

   !> Collect all solver comparison tests
   subroutine collect_math_solvers(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("rosenbrock_newton", test_rosenbrock_newton), &
         & new_unittest("newton_custom_step", test_newton_custom_step), &
         & new_unittest("newton_custom_step_linesearch", test_newton_custom_step_linesearch), &
         & new_unittest("newton_custom_step_requires_callback", test_newton_custom_step_requires_callback), &
         & new_unittest("circle_newton_lagrangian", test_circle_newton_lagrangian), &
         & new_unittest("rosenbrock_lbfgsb", test_rosenbrock_lbfgsb), &
         & new_unittest("circle_lbfgsb_penalty", test_circle_lbfgsb_penalty), &
         & new_unittest("lbfgsb_kernel_driver1_default", test_lbfgsb_driver1), &
         & new_unittest("lbfgsb_kernel_driver2_eval_limit", test_lbfgsb_driver2), &
         & new_unittest("lbfgsb_kernel_driver3_large_budget", test_lbfgsb_driver3), &
         & new_unittest("fmin_kernel_sin", test_fmin_sin), &
         & new_unittest("fmin_kernel_parabola", test_fmin_parabola), &
         & new_unittest("newton_kernel_dense_sweep", test_newton_dense), &
         & new_unittest("newton_kernel_sparse_sweep", test_newton_sparse), &
         & new_unittest("lbfgsb_invalid_memory", test_lbfgsb_invalid_memory, should_fail=.true.), &
         & new_unittest("newton_maxiter_fail", test_newton_maxiter_fail, should_fail=.true.), &
         & new_unittest("lbfgsb_kernel_bound_active", test_lbfgsb_bound_active) &
         & ]
   end subroutine collect_math_solvers

   !> Solve a linear residual with a caller-provided Newton step.
   subroutine test_newton_custom_step(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp) :: x(1)
      integer :: context

      context = 0
      x = 0.0_wp
      call new_newton_solver(solver, n=1, m=1, error=solver_error, &
                             func_ctx=custom_step_residual, step_ctx=custom_step_direction, &
                             context=context, tol=1.0e-12_wp)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if
      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if
      call check(error, x(1), 2.0_wp, thr=1.0e-12_wp)
   end subroutine test_newton_custom_step

   !> Custom-step mode without a step callback is refused before initialization
   subroutine test_newton_custom_step_requires_callback(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      integer :: context

      context = 0
      call new_newton_solver(solver, n=1, m=1, error=solver_error, &
                             func_ctx=custom_step_residual, context=context, &
                             linear_solver=newton_linear_solver_custom_step, tol=1.0e-12_wp)
      call check(error, allocated(solver_error), message="Custom step without step_ctx accepted")
      if (allocated(error)) return
      call check(error, index(solver_error%message, "step_ctx") > 0, &
                 message="Unexpected error: "//solver_error%message)
      if (allocated(error)) return
      call check(error, .not. allocated(solver), message="Solver allocated despite the error")
   end subroutine test_newton_custom_step_requires_callback

   !> Residual for the custom step test.
   subroutine custom_step_residual(x, f, context)
      real(wp), intent(in) :: x(:)
      real(wp), intent(out) :: f(:)
      class(*), intent(in) :: context
      f = x - 2.0_wp
   end subroutine custom_step_residual

   !> Exact Newton direction for the custom step test.
   subroutine custom_step_direction(x, f, p, context)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: f(:)
      real(wp), intent(out) :: p(:)
      class(*), intent(in) :: context
      p = -f
   end subroutine custom_step_direction

   !> Nonlinear custom steps with backtracking and exact line searches
   subroutine test_newton_custom_step_linesearch(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      integer, parameter :: modes(2) = [newton_linesearch_backtrack, newton_linesearch_exact]
      real(wp) :: x(1), f(1)
      integer :: context, mode

      context = 0
      do mode = 1, size(modes)
         ! Full Newton step overshoots the root, exercising step-length selection
         x = 0.5_wp
         call new_newton_solver(solver, n=1, m=1, error=solver_error, &
                                func_ctx=custom_nonlinear_residual, step_ctx=custom_nonlinear_direction, &
                                context=context, linesearch=modes(mode), tol=1.0e-12_wp)
         if (allocated(solver_error)) then
            call check(error, .false., message=solver_error%message)
            return
         end if
         call solver%solve(x, solver_error)
         if (allocated(solver_error)) then
            call check(error, .false., message=solver_error%message)
            return
         end if
         call check(error, x(1), 2.0_wp, thr=1.0e-10_wp, &
                    message="Custom Newton root, line search="//to_string(modes(mode)))
         if (allocated(error)) return
         call custom_nonlinear_residual(x, f, context)
         call check(error, abs(f(1)) < 1.0e-10_wp, &
                    message="Custom Newton residual, line search="//to_string(modes(mode)))
         if (allocated(error)) return
      end do
   end subroutine test_newton_custom_step_linesearch

   !> Quadratic residual with a positive root at x = 2
   subroutine custom_nonlinear_residual(x, f, context)
      real(wp), intent(in) :: x(:)
      real(wp), intent(out) :: f(:)
      class(*), intent(in) :: context

      f = x**2 - 4.0_wp
   end subroutine custom_nonlinear_residual

   !> Exact Newton direction for the quadratic residual
   subroutine custom_nonlinear_direction(x, f, p, context)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: f(:)
      real(wp), intent(out) :: p(:)
      class(*), intent(in) :: context

      p = -f/(2.0_wp*x)
   end subroutine custom_nonlinear_direction

   !* ================================================================================= *!
   !*                Test Problem 1: Rosenbrock Function (Unconstrained)                *!
   !* ================================================================================= *!
   !> Physical problem: Minimize Rosenbrock function
   !>   f(x,y) = (1-x)^2 + 100(y-x^2 )^2
   !> Solution: x = 1, y = 1
   !>
   !> Newton formulation: Solve gradient = 0
   !>   f_1 = df/ dx = -2(1-x) - 400x(y-x^2 ) = 0
   !>   f_2 = df/ dy = 200(y-x^2 ) = 0
   !>
   !> L-BFGS-B formulation: Minimize f with bounds (equivalent)

   !> Test Rosenbrock with Newton solver
   subroutine test_rosenbrock_newton(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x

      ! Initial guess (same for all solvers)
      x = [0.5_wp, 0.5_wp]

      ! Initialize Newton solver with gradient = 0 formulation
      call new_newton_solver(solver, &
                             n=2, m=2, &
                             func=rosenbrock_gradient_residual, &
                             grad=rosenbrock_hessian, &
                             error=solver_error, &
                             tol=1.0e-12_wp)

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Solve
      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Verify solution
      call check(error, x(1), 1.0_wp, thr=newton_thr)
      if (allocated(error)) return
      call check(error, x(2), 1.0_wp, thr=newton_thr)

   end subroutine test_rosenbrock_newton

   !> Test Rosenbrock with L-BFGS-B solver
   subroutine test_rosenbrock_lbfgsb(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x
      real(wp), dimension(2) :: xl, xu

      ! Initial guess (same as Newton)
      x = [0.5_wp, 0.5_wp]

      ! Unbounded problem (set very large bounds)
      xl = -1.0e10_wp
      xu = 1.0e10_wp

      ! Initialize L-BFGS-B solver with direct minimization
      call new_lbfgsb_solver( &
         solver=solver, &
         n=2, &
         error=solver_error, &
         obj=rosenbrock_objective, &
         obj_grad=rosenbrock_objective_gradient, &
         l=xl, u=xu, &
         factr=1.0e4_wp, &
         pgtol=1.0e-10_wp)

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Solve
      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Verify solution (use L-BFGS-B tolerance)
      call check(error, x(1), 1.0_wp, thr=lbfgsb_thr)
      if (allocated(error)) return
      call check(error, x(2), 1.0_wp, thr=lbfgsb_thr)

   end subroutine test_rosenbrock_lbfgsb

   !* ================================================================================= *!
   !*                   Test Problem 2: Constrained Circle Projection                   *!
   !* ================================================================================= *!
   !> Physical problem: Find point on unit circle closest to target (2, 3)
   !>   minimize: (x-2)^2 + (y-3)^2
   !>   subject to: x^2 + y^2 = 1
   !> Solution: x = 2/ sqrt 13 ~= 0.5547, y = 3/ sqrt 13 ~= 0.8321
   !>
   !> Newton formulation: Augmented Lagrangian (3 variables: x, y, lambda )
   !>   f_1 = 2(x-2) + 2 lambda x = 0
   !>   f_2 = 2(y-3) + 2 lambda y = 0
   !>   f_3 = x^2 + y^2 - 1 = 0
   !>
   !> L-BFGS-B formulation: Penalty method
   !>   minimize: (x-2)^2 + (y-3)^2 + 10^6 (x^2 + y^2 - 1)^2

   !> Test circle projection with Newton (Lagrangian formulation)
   subroutine test_circle_newton_lagrangian(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(3) :: x  ! [x, y, lambda]
      real(wp) :: x_expected, y_expected

      ! Expected solution
      x_expected = 2.0_wp/sqrt(13.0_wp)  ! ~ 0.5547
      y_expected = 3.0_wp/sqrt(13.0_wp)  ! ~ 0.8321

      ! Initial guess: point on circle + Lagrange multiplier estimate
      ! Start from normalized direction toward target
      x = [0.5547_wp, 0.8321_wp, -2.5_wp]  ! [x, y, lambda_guess]

      ! Initialize Newton solver with Lagrangian formulation
      call new_newton_solver(solver, &
                             n=3, m=3, &
                             func=circle_lagrangian_residual, &
                             grad=circle_lagrangian_jacobian, &
                             error=solver_error, &
                             tol=1.0e-12_wp)

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Solve
      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Verify solution (only check x, y; lambda is auxiliary)
      call check(error, x(1), x_expected, thr=newton_thr)
      if (allocated(error)) return
      call check(error, x(2), y_expected, thr=newton_thr)
      if (allocated(error)) return

      ! Verify constraint satisfaction
      call check(error, x(1)**2 + x(2)**2, 1.0_wp, thr=newton_thr)

   end subroutine test_circle_newton_lagrangian

   !> Test circle projection with L-BFGS-B (penalty formulation)
   subroutine test_circle_lbfgsb_penalty(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x  ! [x, y]
      real(wp), dimension(2) :: xl, xu
      real(wp) :: x_expected, y_expected

      ! Expected solution
      x_expected = 2.0_wp/sqrt(13.0_wp)
      y_expected = 3.0_wp/sqrt(13.0_wp)

      ! Set reasonable bounds
      xl = [-2.0_wp, -2.0_wp]
      xu = [2.0_wp, 2.0_wp]

      ! Initial guess
      x = [0.5_wp, 0.75_wp]

      ! Initialize L-BFGS-B solver with penalty objective
      call new_lbfgsb_solver( &
         solver=solver, &
         n=2, &
         error=solver_error, &
         obj=circle_penalty_objective, &
         obj_grad=circle_penalty_objective_gradient, &
         l=xl, u=xu, &
         factr=1.0e4_wp, &
         pgtol=1.0e-12_wp)

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Solve
      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      ! Verify solution
      call check(error, x(1), x_expected, thr=lbfgsb_thr)
      if (allocated(error)) return
      call check(error, x(2), y_expected, thr=lbfgsb_thr)
      if (allocated(error)) return

      ! Verify constraint satisfaction
      call check(error, x(1)**2 + x(2)**2, 1.0_wp, thr=lbfgsb_thr)

   end subroutine test_circle_lbfgsb_penalty

   !* ================================================================================= *!
   !*                            Rosenbrock Problem Functions                           *!
   !* ================================================================================= *!

   !> Rosenbrock gradient residual for Newton: df = 0
   subroutine rosenbrock_gradient_residual(x, f)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f

      ! f_1 = df/ dx = -2(1-x) - 400x(y-x^2 )
      f(1) = -2.0_wp*(1.0_wp - x(1)) - 400.0_wp*x(1)*(x(2) - x(1)**2)

      ! f_2 = df/ dy = 200(y-x^2 )
      f(2) = 200.0_wp*(x(2) - x(1)**2)
   end subroutine rosenbrock_gradient_residual

   !> Rosenbrock Hessian for Newton
   subroutine rosenbrock_hessian(x, jac)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: jac

      !  d^2 f/ dx^2 = 2 - 400(y - 3x^2 )
      jac(1, 1) = 2.0_wp - 400.0_wp*(x(2) - 3.0_wp*x(1)**2)

      !  d^2 f/ dx dy = -400x
      jac(1, 2) = -400.0_wp*x(1)

      !  d^2 f/ dy dx = -400x
      jac(2, 1) = -400.0_wp*x(1)

      !  d^2 f/ dy^2 = 200
      jac(2, 2) = 200.0_wp
   end subroutine rosenbrock_hessian

   !> Rosenbrock objective
   subroutine rosenbrock_objective(x, f)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f

      ! f = (1-x)^2 + 100(y-x^2 )^2
      f = (1.0_wp - x(1))**2 + 100.0_wp*(x(2) - x(1)**2)**2
   end subroutine rosenbrock_objective

   !> Rosenbrock gradient
   subroutine rosenbrock_objective_gradient(x, df)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df

      !  df/ dx = -2(1-x) - 400x(y-x^2 )
      df(1) = -2.0_wp*(1.0_wp - x(1)) - 400.0_wp*x(1)*(x(2) - x(1)**2)

      !  df/ dy = 200(y-x^2 )
      df(2) = 200.0_wp*(x(2) - x(1)**2)
   end subroutine rosenbrock_objective_gradient

   !* ================================================================================= *!
   !*                              Circle Problem Functions                             *!
   !* ================================================================================= *!

   !> Circle Lagrangian residual for Newton: [ dL, constraint] = 0
   !> Variables: [x, y, lambda ]
   subroutine circle_lagrangian_residual(x, f)
      real(wp), dimension(:), intent(in) :: x  ! [x, y, lambda]
      real(wp), dimension(:), intent(out) :: f

      ! f_1 = dL/ dx = 2(x-2) + 2 lambda x = 0
      f(1) = 2.0_wp*(x(1) - 2.0_wp) + 2.0_wp*x(3)*x(1)

      ! f_2 = dL/ dy = 2(y-3) + 2 lambda y = 0
      f(2) = 2.0_wp*(x(2) - 3.0_wp) + 2.0_wp*x(3)*x(2)

      ! f_3 = constraint = x^2 + y^2 - 1 = 0
      f(3) = x(1)**2 + x(2)**2 - 1.0_wp
   end subroutine circle_lagrangian_residual

   !> Circle Lagrangian Jacobian for Newton
   subroutine circle_lagrangian_jacobian(x, jac)
      real(wp), dimension(:), intent(in) :: x  ! [x, y, lambda]
      real(wp), dimension(:, :), intent(out) :: jac

      ! Row 1: derivatives of f_1
      jac(1, 1) = 2.0_wp + 2.0_wp*x(3)  ! df1/dx
      jac(1, 2) = 0.0_wp                   ! df1/dy
      jac(1, 3) = 2.0_wp*x(1)            ! df1/dlambda

      ! Row 2: derivatives of f_2
      jac(2, 1) = 0.0_wp                   ! df2/dx
      jac(2, 2) = 2.0_wp + 2.0_wp*x(3)  ! df2/dy
      jac(2, 3) = 2.0_wp*x(2)            ! df2/dlambda

      ! Row 3: derivatives of constraint
      jac(3, 1) = 2.0_wp*x(1)            ! df3/dx
      jac(3, 2) = 2.0_wp*x(2)            ! df3/dy
      jac(3, 3) = 0.0_wp                   ! df3/dlambda
   end subroutine circle_lagrangian_jacobian

   !> Circle penalty objective for L-BFGS-B: distance + constraint penalty
   subroutine circle_penalty_objective(x, f)
      real(wp), dimension(:), intent(in) :: x  ! [x, y]
      real(wp), intent(out) :: f
      real(wp) :: constraint
      real(wp), parameter :: penalty = 1.0e7_wp

      ! Constraint violation: c = x^2 + y^2 - 1
      constraint = x(1)**2 + x(2)**2 - 1.0_wp

      ! f = (x-2)^2 + (y-3)^2 + penalty . c^2
      f = (x(1) - 2.0_wp)**2 + (x(2) - 3.0_wp)**2 + penalty*constraint**2
   end subroutine circle_penalty_objective

   !> Circle penalty objective gradient for L-BFGS-B
   subroutine circle_penalty_objective_gradient(x, df)
      real(wp), dimension(:), intent(in) :: x  ! [x, y]
      real(wp), dimension(:), intent(out) :: df
      real(wp) :: constraint
      real(wp), parameter :: penalty = 1.0e7_wp

      ! Constraint violation: c = x^2 + y^2 - 1
      constraint = x(1)**2 + x(2)**2 - 1.0_wp

      !  df/ dx = 2(x-2) + penalty . 2c . 2x
      df(1) = 2.0_wp*(x(1) - 2.0_wp) + penalty*2.0_wp*constraint*2.0_wp*x(1)

      !  df/ dy = 2(y-3) + penalty . 2c . 2y
      df(2) = 2.0_wp*(x(2) - 3.0_wp) + penalty*2.0_wp*constraint*2.0_wp*x(2)
   end subroutine circle_penalty_objective_gradient

   !* ================================================================================= *!
   !*       Raw-kernel ports: jacobwilliams/lbfgsb (setulb reverse communication)       *!
   !* ================================================================================= *!
   ! All three drivers minimize the same bound-constrained 25-variable extended
   ! Rosenbrock problem (optimum f = 0); they differ only in termination policy
   ! File/console output is suppressed via iprint = -1

   !> driver1: run to the default convergence test (f -> 0)
   subroutine test_lbfgsb_driver1(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: f
      character(len=60) :: final_task

      call run_lbfgsb_kernel(0, f, final_task)
      call check(error, final_task(1:4) == "CONV", &
                 "driver1 did not converge: "//trim(final_task))
      if (allocated(error)) return
      call check(error, f < 1.0e-4_wp, "driver1 final objective not near zero")
   end subroutine test_lbfgsb_driver1

   !> driver2: stop after at most 99 function/gradient evaluations; the custom
   !> STOP (or early convergence) must fire and the objective must decrease
   subroutine test_lbfgsb_driver2(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: f, f0, g0(lbfgsb_n), x0(lbfgsb_n)
      character(len=60) :: final_task

      x0 = 3.0_wp
      call lbfgsb_rosenbrock_eval(x0, f0, g0)

      call run_lbfgsb_kernel(99, f, final_task)
      call check(error, final_task(1:4) == "STOP" .or. final_task(1:4) == "CONV", &
                 "driver2 did not stop as expected: "//trim(final_task))
      if (allocated(error)) return
      call check(error, f < f0, "driver2 made no progress")
   end subroutine test_lbfgsb_driver2

   !> driver3: with a large evaluation budget the solver should reach f -> 0
   subroutine test_lbfgsb_driver3(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: f
      character(len=60) :: final_task

      call run_lbfgsb_kernel(900, f, final_task)
      call check(error, final_task(1:4) == "CONV" .or. final_task(1:4) == "STOP", &
                 "driver3 terminated abnormally: "//trim(final_task))
      if (allocated(error)) return
      call check(error, f < 1.0e-4_wp, "driver3 final objective not near zero")
   end subroutine test_lbfgsb_driver3

   !> Drive the setulb reverse-communication loop on the sample problem. If
   !> max_nfg > 0 a custom STOP is issued once isave(34) reaches max_nfg or the
   !> projected gradient becomes negligible; otherwise only built-in tests stop
   subroutine run_lbfgsb_kernel(max_nfg, f, final_task)
      integer, intent(in) :: max_nfg
      real(wp), intent(out) :: f
      character(len=60), intent(out) :: final_task

      integer, parameter :: nws = 2*lbfgsb_m*lbfgsb_n + 5*lbfgsb_n + 11*lbfgsb_m*lbfgsb_m + 8*lbfgsb_m
      real(wp) :: x(lbfgsb_n), l(lbfgsb_n), u(lbfgsb_n), g(lbfgsb_n), wa(nws), dsave(29)
      integer :: nbd(lbfgsb_n), iwa(3*lbfgsb_n), isave(44), i
      character(len=60) :: task, csave
      logical :: lsave(4)

      ! Odd variables: bounds [1, 100]; even variables: bounds [-100, 100]
      do i = 1, lbfgsb_n, 2
         nbd(i) = 2
         l(i) = 1.0_wp
         u(i) = 100.0_wp
      end do
      do i = 2, lbfgsb_n, 2
         nbd(i) = 2
         l(i) = -100.0_wp
         u(i) = 100.0_wp
      end do

      x = 3.0_wp
      f = 0.0_wp
      g = 0.0_wp

      task = "START"
      do while (task(1:2) == "FG" .or. task(1:5) == "NEW_X" .or. task(1:5) == "START")
         call setulb(lbfgsb_n, lbfgsb_m, x, l, u, nbd, f, g, lbfgsb_factr, lbfgsb_pgtol, &
                     wa, iwa, task, -1, csave, lsave, isave, dsave)
         if (task(1:2) == "FG") then
            call lbfgsb_rosenbrock_eval(x, f, g)
         else if (task(1:5) == "NEW_X" .and. max_nfg > 0) then
            if (isave(34) >= max_nfg) then
               task = "STOP: total f and g evaluation limit reached"
            end if
            if (dsave(13) <= 1.0e-10_wp*(1.0_wp + abs(f))) then
               task = "STOP: projected gradient sufficiently small"
            end if
         end if
      end do

      final_task = task
   end subroutine run_lbfgsb_kernel

   !> Extended Rosenbrock objective and gradient (the L-BFGS-B sample problem)
   subroutine lbfgsb_rosenbrock_eval(x, f, g)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      real(wp), dimension(:), intent(out) :: g

      real(wp) :: t1, t2
      integer :: i

      f = 0.25_wp*(x(1) - 1.0_wp)**2
      do i = 2, lbfgsb_n
         f = f + (x(i) - x(i - 1)**2)**2
      end do
      f = 4.0_wp*f

      t1 = x(2) - x(1)**2
      g(1) = 2.0_wp*(x(1) - 1.0_wp) - 16.0_wp*x(1)*t1
      do i = 2, lbfgsb_n - 1
         t2 = t1
         t1 = x(i + 1) - x(i)**2
         g(i) = 8.0_wp*t2 - 16.0_wp*x(i)*t1
      end do
      g(lbfgsb_n) = 8.0_wp*t1
   end subroutine lbfgsb_rosenbrock_eval

   !* ================================================================================= *!
   !*        Raw-kernel ports: jacobwilliams/fmin (1D derivative-free minimizer)        *!
   !* ================================================================================= *!

   !> fmin test.f90: minimize sin(x) on [-4, 0]; the minimum is at x = -pi/2
   subroutine test_fmin_sin(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: pi = acos(-1.0_wp)
      real(wp) :: xmin

      xmin = fmin(fmin_sin_func, -4.0_wp, 0.0_wp, fmin_tol)
      call check(error, xmin, -pi/2.0_wp, thr=fmin_thr)
   end subroutine test_fmin_sin

   !> Minimize (x - 2)^2 on [0, 5]; the minimum is at x = 2
   subroutine test_fmin_parabola(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: xmin

      xmin = fmin(fmin_parabola_func, 0.0_wp, 5.0_wp, fmin_tol)
      call check(error, xmin, 2.0_wp, thr=fmin_thr)
   end subroutine test_fmin_parabola

   !> f(x) = sin(x)
   function fmin_sin_func(x) result(f)
      real(wp), intent(in) :: x
      real(wp) :: f
      f = sin(x)
   end function fmin_sin_func

   !> f(x) = (x - 2)^2
   function fmin_parabola_func(x) result(f)
      real(wp), intent(in) :: x
      real(wp) :: f
      f = (x - 2.0_wp)**2
   end function fmin_parabola_func

   !* ================================================================================= *!
   !*         Raw-kernel ports: jacobwilliams/nlesolver-fortran (nlesolver_type)        *!
   !* ================================================================================= *!
   ! Both solve f1 = x1^2 + x2 - 0.1 = 0, f2 = x2 + 0.2 = 0; root (sqrt(0.3),
   ! -0.2). Dense sweeps line-search/Broyden; sparse drives LSQR/LUSOL/LSMR

   !> nlesolver_test_1.f90: dense solve across all four line-search step modes
   !> with the Broyden update both off and on, using scalar variable bounds
   subroutine test_newton_dense(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: step_mode, b
      logical :: use_broyden

      do step_mode = 1, 4
         do b = 0, 1
            use_broyden = (b == 1)
            call run_newton_dense(step_mode, use_broyden, &
                                  "dense step_mode="//to_string(step_mode)//" broyden="//merge("T", "F", use_broyden), error)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_newton_dense

   !> sparse_test.f90: same line-search/Broyden sweep, run through each of the
   !> sparse linear-solver backends (LSQR, LUSOL, LSMR)
   subroutine test_newton_sparse(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: modes(3), k, step_mode, b
      logical :: use_broyden
      character(len=16) :: names(3)

      modes = [NLESOLVER_SPARSITY_LSQR, NLESOLVER_SPARSITY_LUSOL, NLESOLVER_SPARSITY_LSMR]
      names = ["LSQR ", "LUSOL", "LSMR "]

      do k = 1, 3
         do step_mode = 1, 4
            do b = 0, 1
               use_broyden = (b == 1)
               call run_newton_sparse(modes(k), step_mode, use_broyden, &
                                      "sparse "//trim(names(k))//" step_mode="//to_string(step_mode)// &
                                      " broyden="//merge("T", "F", use_broyden), error)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_newton_sparse

   !> Run one dense configuration and verify convergence to the root
   subroutine run_newton_dense(step_mode, use_broyden, label, error)
      integer, intent(in) :: step_mode
      logical, intent(in) :: use_broyden
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      type(nlesolver_type) :: solver
      real(wp) :: x(newton_n)
      integer :: istat
      character(len=:), allocatable :: msg

      call solver%initialize(n=newton_n, m=newton_m, max_iter=newton_max_iter, tol=newton_tol, &
                             func=newton_func, grad=newton_grad, step_mode=step_mode, &
                             use_broyden=use_broyden, n_intervals=2, fmin_tol=1.0e-2_wp, &
                             verbose=.false., bounds_mode=NLESOLVER_SCALAR_BOUNDS, &
                             xlow=[0.0_wp, -5.0_wp], xupp=[1.0_wp, 0.0_wp])
      call solver%status(istat, msg)
      call check(error, istat == 0, label//": initialize failed: "//msg)
      if (allocated(error)) return

      x = [1.0_wp, 2.0_wp]
      call solver%solve(x)
      call solver%status(istat, msg)
      call check(error, istat == 1, label//": solve did not converge: "//msg)
      if (allocated(error)) return
      call check_newton_root(x, label, error)

   end subroutine run_newton_dense

   !> Run one sparse configuration and verify convergence to the root
   subroutine run_newton_sparse(sparsity_mode, step_mode, use_broyden, label, error)
      integer, intent(in) :: sparsity_mode
      integer, intent(in) :: step_mode
      logical, intent(in) :: use_broyden
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      type(nlesolver_type) :: solver
      real(wp) :: x(newton_n)
      integer :: istat
      character(len=:), allocatable :: msg

      call solver%initialize(n=newton_n, m=newton_m, max_iter=newton_max_iter, tol=newton_tol, &
                             atol=newton_tol, btol=newton_tol, func=newton_func, &
                             grad_sparse=newton_grad_sparse, step_mode=step_mode, &
                             use_broyden=use_broyden, n_intervals=2, fmin_tol=1.0e-2_wp, &
                             verbose=.false., sparsity_mode=sparsity_mode, &
                             irow=newton_irow, icol=newton_icol, damp=0.0_wp)
      call solver%status(istat, msg)
      call check(error, istat == 0, label//": initialize failed: "//msg)
      if (allocated(error)) return

      x = [1.0_wp, 2.0_wp]
      call solver%solve(x)
      call solver%status(istat, msg)
      call check(error, istat == 1, label//": solve did not converge: "//msg)
      if (allocated(error)) return
      call check_newton_root(x, label, error)

   end subroutine run_newton_sparse

   !> Assert that x is the expected root and the residual norm is small
   subroutine check_newton_root(x, label, error)
      real(wp), dimension(newton_n), intent(in) :: x
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: f(newton_m)

      f(1) = x(1)**2 + x(2) - 0.1_wp
      f(2) = x(2) + 0.2_wp

      call check(error, norm2(f) < 1.0e-6_wp, label//": residual not converged")
      if (allocated(error)) return
      call check(error, x(1), newton_x_ref(1), thr=newton_x_thr, message=label//": x(1) mismatch")
      if (allocated(error)) return
      call check(error, x(2), newton_x_ref(2), thr=newton_x_thr, message=label//": x(2) mismatch")
   end subroutine check_newton_root

   !> Residual vector of the test system (nlesolver_module func interface)
   subroutine newton_func(me, x, f)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f

      f(1) = x(1)**2 + x(2) - 0.1_wp
      f(2) = x(2) + 0.2_wp
   end subroutine newton_func

   !> Dense Jacobian of the test system
   subroutine newton_grad(me, x, g)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: g

      g(1, 1) = 2.0_wp*x(1)
      g(2, 1) = 0.0_wp
      g(1, 2) = 1.0_wp
      g(2, 2) = 1.0_wp
   end subroutine newton_grad

   !> Sparse Jacobian packed into the COO pattern (newton_irow, newton_icol)
   subroutine newton_grad_sparse(me, x, g)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g

      real(wp) :: g_dense(newton_m, newton_n)

      call newton_grad(me, x, g_dense)
      g(1) = g_dense(1, 1)
      g(2) = g_dense(1, 2)
      g(3) = g_dense(2, 2)
   end subroutine newton_grad_sparse

   !* ================================================================================= *!
   !*       Failure-path tests for the solver wrappers (error-reporting coverage)       *!
   !* ================================================================================= *!

   !> L-BFGS-B construction must reject an out-of-range memory parameter m
   !> The valid range is [3, 20] (lbfgsb.f90 lines 225-229); m=1 is below the
   !> minimum, so new_lbfgsb_solver must allocate solver_error and must NOT
   !> allocate the solver. Pure construction-time validation: fully deterministic
   subroutine test_lbfgsb_invalid_memory(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: xl, xu
      integer, parameter :: bad_m = 1   ! below the documented minimum of 3

      xl = -1.0e10_wp
      xu = 1.0e10_wp

      call new_lbfgsb_solver( &
         solver=solver, &
         n=2, &
         error=solver_error, &
         obj=rosenbrock_objective, &
         obj_grad=rosenbrock_objective_gradient, &
         l=xl, u=xu, &
         m=bad_m)

      ! Expected failure (registered should_fail=.true.): the out-of-range memory
      ! parameter must be rejected at construction. Propagating solver_error makes
      ! testdrive report EXPECTED FAIL; an UNEXPECTED PASS flags a regression
      if (allocated(solver_error)) call test_failed(error, trim(solver_error%message))
   end subroutine test_lbfgsb_invalid_memory

   !> Newton must surface non-convergence when iterations are exhausted
   subroutine test_newton_maxiter_fail(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x

      ! Far start so a single step cannot converge or stagnate
      x = [-1.2_wp, 1.0_wp]

      call new_newton_solver(solver, &
                             n=2, m=2, &
                             func=rosenbrock_gradient_residual, &
                             grad=rosenbrock_hessian, &
                             error=solver_error, &
                             max_iter=1, &
                             tol=1.0e-12_wp, &
                             tolx=1.0e-14_wp)

      ! Construction is expected to succeed; the one-iteration solve from a far
      ! start must then report non-convergence (registered should_fail=.true.)
      if (allocated(solver)) then
         call solver%solve(x, solver_error)
      end if
      if (allocated(solver_error)) call test_failed(error, trim(solver_error%message))
   end subroutine test_newton_maxiter_fail

   !* ================================================================================= *!
   !*             L-BFGS-B raw-kernel: active-bound and error-state coverage            *!
   !* ================================================================================= *!

   !> Bound-active coverage
   !> Minimization of the convex quadratic sum_i (x_i - 3)^2 with an upper bound
   !> u_i = 1 (nbd = 3, upper-only) -> the unconstrained optimum x_i = 3 lies OUTSIDE
   !> the box, so the constrained optimum is clipped onto  the active upper bound x_i = 1
   !> with f = n*(1-3)^2 = 12
   subroutine test_lbfgsb_bound_active(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 3, m = 5
      integer, parameter :: nws = 2*m*n + 5*n + 11*m*m + 8*m
      real(wp) :: x(n), l(n), u(n), g(n), f, wa(nws), dsave(29)
      integer :: nbd(n), iwa(3*n), isave(44), i
      character(len=60) :: task, csave
      logical :: lsave(4)

      ! Upper bound only: nbd(i) = 3, u(i) = 1; l(i) is unused for nbd = 3
      do i = 1, n
         nbd(i) = 3
         u(i) = 1.0_wp
         l(i) = 0.0_wp
      end do

      x = -2.0_wp     ! feasible start, below the active bound
      f = 0.0_wp
      g = 0.0_wp

      task = "START"
      do while (task(1:2) == "FG" .or. task(1:5) == "NEW_X" .or. task(1:5) == "START")
         call setulb(n, m, x, l, u, nbd, f, g, lbfgsb_factr, lbfgsb_pgtol, &
                     wa, iwa, task, -1, csave, lsave, isave, dsave)
         if (task(1:2) == "FG") then
            f = sum((x - 3.0_wp)**2)
            g = 2.0_wp*(x - 3.0_wp)
         end if
      end do

      call check(error, task(1:4) == "CONV", &
                 "bound-active problem did not converge: "//trim(task))
      if (allocated(error)) return
      ! Every variable must pin to its active upper bound u = 1
      call check(error, maxval(abs(x - 1.0_wp)) < 1.0e-5_wp, &
                 "bound-active solution not on the active upper bound")
      if (allocated(error)) return
      call check(error, abs(f - 12.0_wp) < 1.0e-8_wp, &
                 "bound-active objective not n*(1-3)^2 = 12")
      if (allocated(error)) return
   end subroutine test_lbfgsb_bound_active

end module test_math_solvers
