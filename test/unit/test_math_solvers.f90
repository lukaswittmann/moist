!> Comprehensive tests for mathematical solvers with cross-solver comparison
!>
!> This module tests multiple optimization solvers (Newton, SLSQP, L-BFGS-B) on the same
!> problems using different formulations:
!>  - Newton: Solves nonlinear systems f(x)=0, handles constraints via Lagrangian
!>  - SLSQP: Minimizes objectives, handles constraints explicitly
!>  - L-BFGS-B: Minimizes objectives with bound constraints only
!>
!> Test problems:
!>  1. Rosenbrock: Unconstrained optimization (equivalent formulations)
!>  2. Constrained circle: Minimize distance with constraint (Lagrangian/explicit/penalty)
module test_math_solvers
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_type, only: solver_base_type
   use moist_math_solver_newton
   use moist_math_solver_slsqp
   use moist_math_solver_lbfgsb
   use moist_math_solver_slsqp_multi_tangent
   use moist_math_solver_slsqp_deflation, only: new_slsqp_deflation_solver, &
                                                moist_math_solver_slsqp_deflation_type
   use moist_math_solver_newton_deflation, only: new_newton_deflation_solver, &
                                                 moist_math_solver_newton_deflation_type
   use mctc_env_error, only: moist_error_type => error_type
   ! Raw vendored solver APIs (upstream test ports folded into this suite).
   use slsqp_module, only: slsqp_solver
   use lbfgsb_module, only: setulb
   use moist_math_solver_fmin, only: fmin
   use nlesolver_module, only: nlesolver_type, &
      NLESOLVER_SCALAR_BOUNDS, &
      NLESOLVER_SPARSITY_LSQR, NLESOLVER_SPARSITY_LUSOL, NLESOLVER_SPARSITY_LSMR
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private

   public :: collect_math_solvers

   !> Tolerance for Newton solver
   real(wp), parameter :: newton_thr = 1.0e-10_wp
   !> Tolerance for SLSQP solver
   real(wp), parameter :: slsqp_thr = 1.0e-7_wp
   !> Tolerance for L-BFGS-B solver
   real(wp), parameter :: lbfgsb_thr = 1.0e-6_wp

   !===========================================================================
   ! Raw-kernel test constants (vendored APIs ported from upstream test suites)
   !===========================================================================
   !> SLSQP raw-kernel acceptance tolerance.
   real(wp), parameter :: slsqp_kernel_thr = 1.0e-4_wp
   !> Reference solution of the constrained Rosenbrock problem (upstream slsqp).
   real(wp), parameter :: slsqp_rosen_x(2) = &
      [0.78641515097183889_wp, 0.61769831659541152_wp]
   !> L-BFGS-B raw-kernel problem size and convergence controls.
   integer, parameter :: lbfgsb_n = 25, lbfgsb_m = 5
   real(wp), parameter :: lbfgsb_factr = 1.0e7_wp, lbfgsb_pgtol = 1.0e-5_wp
   !> fmin raw-kernel tolerances (upstream accepts 10*tol on the minimizer).
   real(wp), parameter :: fmin_tol = 1.0e-8_wp
   real(wp), parameter :: fmin_thr = 10.0_wp*fmin_tol
   !> Newton/nlesolver raw-kernel problem definition.
   integer, parameter :: newton_n = 2, newton_m = 2, newton_max_iter = 100
   real(wp), parameter :: newton_tol = 1.0e-8_wp
   !> Root of the test system (positive branch).
   real(wp), parameter :: newton_x_ref(2) = [0.5477225575051661_wp, -0.2_wp]
   !> Acceptance tolerance on the converged Newton solution.
   real(wp), parameter :: newton_x_thr = 1.0e-4_wp
   !> COO sparsity pattern of the Jacobian (structural zero at (2,1) omitted).
   integer, parameter :: newton_irow(3) = [1, 1, 2]
   integer, parameter :: newton_icol(3) = [1, 2, 2]

contains

   !> Collect all solver comparison tests
   subroutine collect_math_solvers(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("rosenbrock-newton", test_rosenbrock_newton), &
         & new_unittest("circle-newton-lagrangian", test_circle_newton_lagrangian), &
         & new_unittest("rosenbrock-slsqp", test_rosenbrock_slsqp), &
         & new_unittest("circle-slsqp-explicit", test_circle_slsqp_explicit), &
         & new_unittest("rosenbrock-slsqp-multistart", test_rosenbrock_slsqp_multistart), &
         & new_unittest("circle-slsqp-multistart", test_circle_slsqp_multistart), &
         & new_unittest("rosenbrock-lbfgsb", test_rosenbrock_lbfgsb), &
         & new_unittest("circle-lbfgsb-penalty", test_circle_lbfgsb_penalty), &
         & new_unittest("slsqp-deflation-two-circle-union", test_slsqp_deflation_two_circle_union), &
         & new_unittest("newton-deflation-cubic", test_newton_deflation_cubic), &
         & new_unittest("slsqp-kernel-rosenbrock-inequality", test_slsqp_rosenbrock), &
         & new_unittest("slsqp-kernel-quadratic-eq-ineq", test_slsqp_quadratic), &
         & new_unittest("slsqp-kernel-fd-gradients", test_slsqp_fd_gradients), &
         & new_unittest("slsqp-kernel-hock-schittkowski-71", test_slsqp_hs71), &
         & new_unittest("slsqp-kernel-stopping-nan-bounds", test_slsqp_stopping), &
         & new_unittest("lbfgsb-kernel-driver1-default", test_lbfgsb_driver1), &
         & new_unittest("lbfgsb-kernel-driver2-eval-limit", test_lbfgsb_driver2), &
         & new_unittest("lbfgsb-kernel-driver3-large-budget", test_lbfgsb_driver3), &
         & new_unittest("fmin-kernel-sin", test_fmin_sin), &
         & new_unittest("fmin-kernel-parabola", test_fmin_parabola), &
         & new_unittest("newton-kernel-dense-sweep", test_newton_dense), &
         & new_unittest("newton-kernel-sparse-sweep", test_newton_sparse) &
         & ]
   end subroutine collect_math_solvers

   !==============================================================================
   ! Test Problem 1: Rosenbrock Function (Unconstrained)
   !==============================================================================
   !> Physical problem: Minimize Rosenbrock function
   !>   f(x,y) = (1-x)^2 + 100(y-x^2 )^2
   !> Solution: x = 1, y = 1
   !>
   !> Newton formulation: Solve gradient = 0
   !>   f_1 = df/ dx = -2(1-x) - 400x(y-x^2 ) = 0
   !>   f_2 = df/ dy = 200(y-x^2 ) = 0
   !>
   !> SLSQP formulation: Minimize f directly (equivalent)
   !>
   !> L-BFGS-B formulation: Minimize f with bounds (equivalent)

   !> Test Rosenbrock with Newton solver
   subroutine test_rosenbrock_newton(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x
      integer :: istat

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

      call solver%destroy()

   end subroutine test_rosenbrock_newton

   !> Test Rosenbrock with SLSQP solver
   subroutine test_rosenbrock_slsqp(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x
      real(wp), dimension(2) :: xl, xu
      integer :: istat

      ! Initial guess (same as Newton)
      x = [0.5_wp, 0.5_wp]

      ! Unbounded problem (set very large bounds)
      xl = -1.0e10_wp
      xu = 1.0e10_wp

      ! Initialize SLSQP solver with direct minimization (no constraints)
      call new_slsqp_solver( &
         solver=solver, &
         n=2, m=0, meq=0, &
         error=solver_error, &
         obj=rosenbrock_objective, &
         obj_grad=rosenbrock_objective_gradient, &
         xl=xl, xu=xu, &
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

      ! Verify solution (use looser tolerance for SLSQP)
      call check(error, x(1), 1.0_wp, thr=1.0e-5_wp)
      if (allocated(error)) return
      call check(error, x(2), 1.0_wp, thr=1.0e-5_wp)

      call solver%destroy()
      deallocate(solver)

   end subroutine test_rosenbrock_slsqp

   !> Test Rosenbrock with L-BFGS-B solver
   subroutine test_rosenbrock_lbfgsb(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x
      real(wp), dimension(2) :: xl, xu

      ! Initial guess (same as Newton/SLSQP)
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

      call solver%destroy()
      deallocate(solver)

   end subroutine test_rosenbrock_lbfgsb

   !> Test Rosenbrock with SLSQP multi-start (unconstrained)
   !> Note: multistart uses 3D Lebedev seeds, so we wrap the 2D Rosenbrock
   !> objective/gradient into 3D by ignoring z.
   subroutine test_rosenbrock_slsqp_multistart(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(3) :: x  ! [x, y, z] (z unused)
      real(wp), dimension(3) :: xl, xu

      ! Initial guess (ignored by multistart solver)
      x = [0.0_wp, 0.0_wp, 0.0_wp]

      ! Unbounded problem (set very large bounds)
      xl = -1.0e10_wp
      xu = 1.0e10_wp

      call new_slsqp_multi_tangent_solver( &
         anchor=[0.5_wp, 0.5_wp, 0.0_wp], &
         owner_pos=[0.5_wp, 0.5_wp, -1.0_wp], &
         solver=solver, &
         n=3, m=0, meq=0, &
         obj_ctx=rosenbrock_objective_ctx3, &
         obj_grad_ctx=rosenbrock_objective_gradient_ctx3, &
         con_ctx=empty_constraint_ctx3, &
         con_grad_ctx=empty_constraint_gradient_ctx3, &
         context=1, &
         max_iter=200, &
         tol=1.0e-14_wp, &
         toldx=1.0e-14_wp, &
         toldf=1.0e-14_wp, &
         radii=[0.0001_wp, 0.0_wp], &
         n_points=[6, 6], &
         xl=xl, &
         xu=xu, &
         debug=.false., &
         error=solver_error &
      )

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      call check(error, x(1), 1.0_wp, thr=slsqp_thr)
      if (allocated(error)) return
      call check(error, x(2), 1.0_wp, thr=slsqp_thr)

      call solver%destroy()
      deallocate(solver)
   end subroutine test_rosenbrock_slsqp_multistart


   !==============================================================================
   ! Test Problem 2: Constrained Circle Projection
   !==============================================================================
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
   !> SLSQP formulation: Explicit constraint
   !>   minimize: (x-2)^2 + (y-3)^2
   !>   subject to: x^2 + y^2 - 1 = 0 (equality)
   !>
   !> L-BFGS-B formulation: Penalty method
   !>   minimize: (x-2)^2 + (y-3)^2 + 10^6 (x^2 + y^2 - 1)^2

   !> Test circle projection with Newton (Lagrangian formulation)
   subroutine test_circle_newton_lagrangian(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(3) :: x  ! [x, y, λ]
      real(wp) :: x_expected, y_expected
      integer :: istat

      ! Expected solution
      x_expected = 2.0_wp / sqrt(13.0_wp)  ! ≈ 0.5547
      y_expected = 3.0_wp / sqrt(13.0_wp)  ! ≈ 0.8321

      ! Initial guess: point on circle + Lagrange multiplier estimate
      ! Start from normalized direction toward target
      x = [0.5547_wp, 0.8321_wp, -2.5_wp]  ! [x, y, λ_guess]

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

      call solver%destroy()

   end subroutine test_circle_newton_lagrangian

   !> Test circle projection with SLSQP (explicit constraint)
   subroutine test_circle_slsqp_explicit(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x  ! [x, y]
      real(wp), dimension(2) :: xl, xu
      real(wp) :: x_expected, y_expected
      integer :: istat

      ! Expected solution
      x_expected = 2.0_wp / sqrt(13.0_wp)
      y_expected = 3.0_wp / sqrt(13.0_wp)

      ! Set reasonable bounds for constrained optimization
      xl = [-2.0_wp, -2.0_wp]
      xu = [2.0_wp, 2.0_wp]

      ! Initial guess (same point on circle as Newton)
      x = [0.5_wp, 0.75_wp]

      ! Initialize SLSQP solver with explicit constraint
      call new_slsqp_solver( &
         solver=solver, &
         n=2, m=1, meq=1, &
         error=solver_error, &
         obj=circle_objective, &
         obj_grad=circle_objective_gradient, &
         con=circle_constraint, &
         con_grad=circle_constraint_gradient, &
         xl=xl, xu=xu, &
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
      call check(error, x(1), x_expected, thr=slsqp_thr)
      if (allocated(error)) return
      call check(error, x(2), y_expected, thr=slsqp_thr)
      if (allocated(error)) return

      ! Verify constraint satisfaction
      call check(error, x(1)**2 + x(2)**2, 1.0_wp, thr=slsqp_thr)

      call solver%destroy()
      deallocate(solver)

   end subroutine test_circle_slsqp_explicit

   !> Test circle projection with SLSQP multi-start (explicit constraint)
   subroutine test_circle_slsqp_multistart(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(3) :: x  ! [x, y, z]
      real(wp), dimension(3) :: xl, xu
      real(wp) :: x_expected, y_expected

      ! Expected solution
      x_expected = 2.0_wp / sqrt(13.0_wp)
      y_expected = 3.0_wp / sqrt(13.0_wp)

      ! Set reasonable bounds
      xl = [-2.0_wp, -2.0_wp, -2.0_wp]
      xu = [2.0_wp, 2.0_wp, 2.0_wp]

      ! Initial guess (ignored by multistart solver)
      x = [0.0_wp, 0.0_wp, 0.0_wp]

      call new_slsqp_multi_tangent_solver( &
         anchor=[2.0_wp, 3.0_wp, 0.0_wp], &
         owner_pos=[2.0_wp, 3.0_wp, -1.0_wp], &
         solver=solver, &
         n=3, m=1, meq=1, &
         obj_ctx=circle_objective_ctx3, &
         obj_grad_ctx=circle_objective_gradient_ctx3, &
         con_ctx=circle_constraint_ctx3, &
         con_grad_ctx=circle_constraint_gradient_ctx3, &
         context=1, &
         max_iter=200, &
         tol=1.0e-12_wp, &
         toldx=1.0e-12_wp, &
         toldf=1.0e-12_wp, &
         xl=xl, &
         xu=xu, &
         debug=.false., &
         error=solver_error &
      )

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      call check(error, x(1), x_expected, thr=slsqp_thr)
      if (allocated(error)) return
      call check(error, x(2), y_expected, thr=slsqp_thr)
      if (allocated(error)) return
      call check(error, x(1)**2 + x(2)**2, 1.0_wp, thr=slsqp_thr)

      call solver%destroy()
      deallocate(solver)
   end subroutine test_circle_slsqp_multistart

   !> Test circle projection with L-BFGS-B (penalty formulation)
   subroutine test_circle_lbfgsb_penalty(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x  ! [x, y]
      real(wp), dimension(2) :: xl, xu
      real(wp) :: x_expected, y_expected

      ! Expected solution
      x_expected = 2.0_wp / sqrt(13.0_wp)
      y_expected = 3.0_wp / sqrt(13.0_wp)

      ! Set reasonable bounds
      xl = [-2.0_wp, -2.0_wp]
      xu = [2.0_wp, 2.0_wp]

      ! Initial guess (same point as SLSQP)
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

      call solver%destroy()
      deallocate(solver)

   end subroutine test_circle_lbfgsb_penalty

   !==============================================================================
   ! Rosenbrock Problem Functions
   !==============================================================================

   !> Rosenbrock gradient residual for Newton: df = 0
   subroutine rosenbrock_gradient_residual(x, f)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f

      ! f_1 = df/ dx = -2(1-x) - 400x(y-x^2 )
      f(1) = -2.0_wp * (1.0_wp - x(1)) - 400.0_wp * x(1) * (x(2) - x(1)**2)

      ! f_2 = df/ dy = 200(y-x^2 )
      f(2) = 200.0_wp * (x(2) - x(1)**2)
   end subroutine rosenbrock_gradient_residual

   !> Rosenbrock Hessian for Newton
   subroutine rosenbrock_hessian(x, jac)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:,:), intent(out) :: jac

      !  d^2 f/ dx^2 = 2 - 400(y - 3x^2 )
      jac(1, 1) = 2.0_wp - 400.0_wp * (x(2) - 3.0_wp * x(1)**2)

      !  d^2 f/ dx dy = -400x
      jac(1, 2) = -400.0_wp * x(1)

      !  d^2 f/ dy dx = -400x
      jac(2, 1) = -400.0_wp * x(1)

      !  d^2 f/ dy^2 = 200
      jac(2, 2) = 200.0_wp
   end subroutine rosenbrock_hessian

   !> Rosenbrock objective for SLSQP
   subroutine rosenbrock_objective(x, f)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f

      ! f = (1-x)^2 + 100(y-x^2 )^2
      f = (1.0_wp - x(1))**2 + 100.0_wp * (x(2) - x(1)**2)**2
   end subroutine rosenbrock_objective

   !> Rosenbrock objective with context for 3D variables (ignore z).
   !> Needed because multistart seeds are 3D.
   subroutine rosenbrock_objective_ctx3(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      class(*), intent(in) :: context

      call rosenbrock_objective(x(1:2), f)
   end subroutine rosenbrock_objective_ctx3

   !> Rosenbrock gradient for SLSQP
   subroutine rosenbrock_objective_gradient(x, df)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df

      !  df/ dx = -2(1-x) - 400x(y-x^2 )
      df(1) = -2.0_wp * (1.0_wp - x(1)) - 400.0_wp * x(1) * (x(2) - x(1)**2)

      !  df/ dy = 200(y-x^2 )
      df(2) = 200.0_wp * (x(2) - x(1)**2)
   end subroutine rosenbrock_objective_gradient

   !> Rosenbrock objective gradient with context for 3D variables (ignore z).
   !> Needed because multistart seeds are 3D.
   subroutine rosenbrock_objective_gradient_ctx3(x, df, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df
      class(*), intent(in) :: context

      df = 0.0_wp
      call rosenbrock_objective_gradient(x(1:2), df(1:2))
   end subroutine rosenbrock_objective_gradient_ctx3

   !> Empty constraint with context for 3D variables (m=0).
   !> Required by multistart's SLSQP interface even when m=0.
   subroutine empty_constraint_ctx3(x, c, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: c
      class(*), intent(in) :: context

      if (size(c) > 0) c = 0.0_wp
   end subroutine empty_constraint_ctx3

   !> Empty constraint gradient with context for 3D variables (m=0).
   !> Required by multistart's SLSQP interface even when m=0.
   subroutine empty_constraint_gradient_ctx3(x, dc, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:,:), intent(out) :: dc
      class(*), intent(in) :: context

      if (size(dc) > 0) dc = 0.0_wp
   end subroutine empty_constraint_gradient_ctx3

   !==============================================================================
   ! Circle Problem Functions
   !==============================================================================

   !> Circle Lagrangian residual for Newton: [ dL, constraint] = 0
   !> Variables: [x, y, lambda ]
   subroutine circle_lagrangian_residual(x, f)
      real(wp), dimension(:), intent(in) :: x  ! [x, y, λ]
      real(wp), dimension(:), intent(out) :: f

      ! f_1 = dL/ dx = 2(x-2) + 2 lambda x = 0
      f(1) = 2.0_wp * (x(1) - 2.0_wp) + 2.0_wp * x(3) * x(1)

      ! f_2 = dL/ dy = 2(y-3) + 2 lambda y = 0
      f(2) = 2.0_wp * (x(2) - 3.0_wp) + 2.0_wp * x(3) * x(2)

      ! f_3 = constraint = x^2 + y^2 - 1 = 0
      f(3) = x(1)**2 + x(2)**2 - 1.0_wp
   end subroutine circle_lagrangian_residual

   !> Circle Lagrangian Jacobian for Newton
   subroutine circle_lagrangian_jacobian(x, jac)
      real(wp), dimension(:), intent(in) :: x  ! [x, y, λ]
      real(wp), dimension(:,:), intent(out) :: jac

      ! Row 1: derivatives of f_1
      jac(1, 1) = 2.0_wp + 2.0_wp * x(3)  ! ∂f₁/∂x
      jac(1, 2) = 0.0_wp                   ! ∂f₁/∂y
      jac(1, 3) = 2.0_wp * x(1)            ! ∂f₁/∂λ

      ! Row 2: derivatives of f_2
      jac(2, 1) = 0.0_wp                   ! ∂f₂/∂x
      jac(2, 2) = 2.0_wp + 2.0_wp * x(3)  ! ∂f₂/∂y
      jac(2, 3) = 2.0_wp * x(2)            ! ∂f₂/∂λ

      ! Row 3: derivatives of constraint
      jac(3, 1) = 2.0_wp * x(1)            ! ∂f₃/∂x
      jac(3, 2) = 2.0_wp * x(2)            ! ∂f₃/∂y
      jac(3, 3) = 0.0_wp                   ! ∂f₃/∂λ
   end subroutine circle_lagrangian_jacobian

   !> Circle objective for SLSQP: distance squared
   subroutine circle_objective(x, f)
      real(wp), dimension(:), intent(in) :: x  ! [x, y]
      real(wp), intent(out) :: f

      ! f = (x-2)^2 + (y-3)^2
      f = (x(1) - 2.0_wp)**2 + (x(2) - 3.0_wp)**2
   end subroutine circle_objective

   !> Circle objective with context (for multi-start SLSQP)
   subroutine circle_objective_ctx(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      class(*), intent(in) :: context

      call circle_objective(x, f)
   end subroutine circle_objective_ctx

   !> Circle objective gradient for SLSQP
   subroutine circle_objective_gradient(x, df)
      real(wp), dimension(:), intent(in) :: x  ! [x, y]
      real(wp), dimension(:), intent(out) :: df

      !  df/ dx = 2(x-2)
      df(1) = 2.0_wp * (x(1) - 2.0_wp)

      !  df/ dy = 2(y-3)
      df(2) = 2.0_wp * (x(2) - 3.0_wp)
   end subroutine circle_objective_gradient

   !> Circle objective gradient with context (for multi-start SLSQP)
   subroutine circle_objective_gradient_ctx(x, df, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df
      class(*), intent(in) :: context

      call circle_objective_gradient(x, df)
   end subroutine circle_objective_gradient_ctx

   !> Circle constraint for SLSQP: x^2 + y^2 = 1
   subroutine circle_constraint(x, c)
      real(wp), dimension(:), intent(in) :: x  ! [x, y]
      real(wp), dimension(:), intent(out) :: c

      ! c = x^2 + y^2 - 1 = 0 (equality)
      c(1) = x(1)**2 + x(2)**2 - 1.0_wp
   end subroutine circle_constraint

   !> Circle constraint with context (for multi-start SLSQP)
   subroutine circle_constraint_ctx(x, c, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: c
      class(*), intent(in) :: context

      call circle_constraint(x, c)
   end subroutine circle_constraint_ctx

   !> Circle constraint gradient for SLSQP
   subroutine circle_constraint_gradient(x, dc)
      real(wp), dimension(:), intent(in) :: x   ! [x, y]
      real(wp), dimension(:,:), intent(out) :: dc  ! [m, n]

      !  dc/ dx = 2x
      dc(1, 1) = 2.0_wp * x(1)

      !  dc/ dy = 2y
      dc(1, 2) = 2.0_wp * x(2)
   end subroutine circle_constraint_gradient

   !> Circle constraint gradient with context (for multi-start SLSQP)
   subroutine circle_constraint_gradient_ctx(x, dc, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:,:), intent(out) :: dc
      class(*), intent(in) :: context

      call circle_constraint_gradient(x, dc)
   end subroutine circle_constraint_gradient_ctx

   !> Circle objective with context for 3D variables (ignore z).
   subroutine circle_objective_ctx3(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      class(*), intent(in) :: context

      call circle_objective(x(1:2), f)
   end subroutine circle_objective_ctx3

   !> Circle objective gradient with context for 3D variables (ignore z).
   subroutine circle_objective_gradient_ctx3(x, df, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df
      class(*), intent(in) :: context

      df = 0.0_wp
      call circle_objective_gradient(x(1:2), df(1:2))
   end subroutine circle_objective_gradient_ctx3

   !> Circle constraint with context for 3D variables (ignore z).
   subroutine circle_constraint_ctx3(x, c, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: c
      class(*), intent(in) :: context

      call circle_constraint(x(1:2), c)
   end subroutine circle_constraint_ctx3

   !> Circle constraint gradient with context for 3D variables (ignore z).
   subroutine circle_constraint_gradient_ctx3(x, dc, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:,:), intent(out) :: dc
      class(*), intent(in) :: context

      dc = 0.0_wp
      call circle_constraint_gradient(x(1:2), dc(1:1, 1:2))
   end subroutine circle_constraint_gradient_ctx3

   !> Circle penalty objective for L-BFGS-B: distance + constraint penalty
   subroutine circle_penalty_objective(x, f)
      real(wp), dimension(:), intent(in) :: x  ! [x, y]
      real(wp), intent(out) :: f
      real(wp) :: constraint
      real(wp), parameter :: penalty = 1.0e7_wp

      ! Constraint violation: c = x^2 + y^2 - 1
      constraint = x(1)**2 + x(2)**2 - 1.0_wp

      ! f = (x-2)^2 + (y-3)^2 + penalty . c^2
      f = (x(1) - 2.0_wp)**2 + (x(2) - 3.0_wp)**2 + penalty * constraint**2
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
      df(1) = 2.0_wp * (x(1) - 2.0_wp) + penalty * 2.0_wp * constraint * 2.0_wp * x(1)

      !  df/ dy = 2(y-3) + penalty . 2c . 2y
      df(2) = 2.0_wp * (x(2) - 3.0_wp) + penalty * 2.0_wp * constraint * 2.0_wp * x(2)
   end subroutine circle_penalty_objective_gradient

   !==============================================================================
   ! Test Problem: SLSQP-deflation two-circle union
   !==============================================================================
   !> Projection of an anchor onto the union of two unit circles centered at
   !> (-2, 0) and (+2, 0), represented by the smooth product-form constraint
   !>   h(x, y) = (||r - c_L||^2 - 1) * (||r - c_R||^2 - 1) = 0.
   !> The set of KKT points for min 0.5||r - a||^2 s.t. h = 0 has FOUR members
   !> (nearest + farthest point on each circle). With anchor a = (0, 0.5),
   !> symmetry about the y-axis pairs them.
   !>
   !> Deflation is expected to enumerate all four KKT points from a single seed
   !> at the anchor. Each point must lie on one of the two unit circles.
   subroutine test_slsqp_deflation_two_circle_union(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(2) :: x
      real(wp), dimension(2) :: xl, xu
      real(wp), allocatable :: roots(:, :)
      integer :: n_roots, i, n_on_left, n_on_right
      integer :: dummy_context

      xl = [-3.0_wp, -3.0_wp]
      xu = [3.0_wp, 3.0_wp]

      ! Seed at the anchor.
      x = [0.0_wp, 0.5_wp]

      dummy_context = 0

      call new_slsqp_deflation_solver( &
         solver=solver, &
         n=2, m=1, meq=1, &
         obj_ctx=two_circle_objective, &
         obj_grad_ctx=two_circle_objective_grad, &
         con_ctx=two_circle_constraint, &
         con_grad_ctx=two_circle_constraint_grad, &
         context=dummy_context, &
         xl=xl, xu=xu, &
         max_iter=200, &
         tol=1.0e-12_wp, &
         toldx=1.0e-12_wp, &
         toldf=1.0e-12_wp, &
         max_roots=6, &
         dedup_tol=1.0e-4_wp, &
         error=solver_error)

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         call solver%destroy()
         deallocate (solver)
         return
      end if

      select type (sd => solver)
      type is (moist_math_solver_slsqp_deflation_type)
         call sd%get_raw_candidates(roots, n_roots)
      class default
         call check(error, .false., message="unexpected solver type")
         call solver%destroy()
         deallocate (solver)
         return
      end select

      ! Deflation should find at least 2 distinct branches (nearest point on
      ! each circle). Finding all 4 KKT points is also acceptable.
      call check(error, n_roots >= 2, .true., &
                 message="deflation must enumerate at least 2 distinct roots")
      if (allocated(error)) then
         call solver%destroy()
         deallocate (solver)
         return
      end if

      ! Every enumerated root must lie on one of the circles.
      n_on_left = 0
      n_on_right = 0
      do i = 1, n_roots
         call check(error, on_two_circle_union(roots(:, i)), .true., &
                    message="root not on circle union")
         if (allocated(error)) then
            call solver%destroy()
            deallocate (solver)
            return
         end if
         if (on_left_circle(roots(:, i))) n_on_left = n_on_left + 1
         if (on_right_circle(roots(:, i))) n_on_right = n_on_right + 1
      end do

      ! Deflation must discover both branches (at least one point on each circle).
      call check(error, n_on_left >= 1, .true., &
                 message="deflation missed the left circle branch")
      if (allocated(error)) then
         call solver%destroy()
         deallocate (solver)
         return
      end if
      call check(error, n_on_right >= 1, .true., &
                 message="deflation missed the right circle branch")

      call solver%destroy()
      deallocate (solver)
   end subroutine test_slsqp_deflation_two_circle_union

   !==============================================================================
   ! Test Problem: Newton-deflation on a scalar cubic
   !==============================================================================
   !> Enumerate all three real roots of F(x) = (x-1)(x+1)(x-3) = x^3 - 3x^2 - x + 3
   !> from a single seed x0 = 0. Newton from x0=0 converges to one root;
   !> iterated deflation must recover the other two.
   subroutine test_newton_deflation_cubic(error)
      type(error_type), allocatable, intent(out) :: error
      class(solver_base_type), allocatable :: solver
      type(moist_error_type), allocatable :: solver_error
      real(wp), dimension(1) :: x
      real(wp), allocatable :: roots(:, :)
      integer :: n_roots, i
      integer :: dummy_context
      logical :: found_m1, found_p1, found_p3

      x(1) = 0.0_wp
      dummy_context = 0

      call new_newton_deflation_solver( &
         solver=solver, &
         n=1, m=1, &
         func_ctx=cubic_residual_ctx, &
         grad_ctx=cubic_jacobian_ctx, &
         context=dummy_context, &
         max_iter=200, &
         tol=1.0e-12_wp, &
         tolx=1.0e-12_wp, &
         max_roots=6, &
         dedup_tol=1.0e-4_wp, &
         error=solver_error)

      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         return
      end if

      call solver%solve(x, solver_error)
      if (allocated(solver_error)) then
         call check(error, .false., message=solver_error%message)
         call solver%destroy()
         deallocate (solver)
         return
      end if

      select type (nd => solver)
      type is (moist_math_solver_newton_deflation_type)
         call nd%get_raw_candidates(roots, n_roots)
      class default
         call check(error, .false., message="unexpected solver type")
         call solver%destroy()
         deallocate (solver)
         return
      end select

      ! Must enumerate all three distinct real roots.
      call check(error, n_roots, 3)
      if (allocated(error)) then
         call solver%destroy()
         deallocate (solver)
         return
      end if

      found_m1 = .false.; found_p1 = .false.; found_p3 = .false.
      do i = 1, n_roots
         if (abs(roots(1, i) + 1.0_wp) < 1.0e-4_wp) found_m1 = .true.
         if (abs(roots(1, i) - 1.0_wp) < 1.0e-4_wp) found_p1 = .true.
         if (abs(roots(1, i) - 3.0_wp) < 1.0e-4_wp) found_p3 = .true.
      end do

      call check(error, found_m1 .and. found_p1 .and. found_p3, .true., &
                 message="Newton-deflation must enumerate all three roots (-1, +1, +3)")

      call solver%destroy()
      deallocate (solver)
   end subroutine test_newton_deflation_cubic

   !------------------------------------------------------------------
   ! SLSQP-deflation helpers (two-circle-union problem)
   !------------------------------------------------------------------

   !> Objective: 0.5 * ||x - anchor||^2 with anchor = (0, 0.5).
   subroutine two_circle_objective(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      class(*), intent(in) :: context

      associate (dummy => context)
      end associate
      f = 0.5_wp * ((x(1) - 0.0_wp)**2 + (x(2) - 0.5_wp)**2)
   end subroutine two_circle_objective

   subroutine two_circle_objective_grad(x, df, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: df
      class(*), intent(in) :: context

      associate (dummy => context)
      end associate
      df(1) = x(1) - 0.0_wp
      df(2) = x(2) - 0.5_wp
   end subroutine two_circle_objective_grad

   !> Constraint h(x) = (||x - c_L||^2 - 1) * (||x - c_R||^2 - 1) = 0
   subroutine two_circle_constraint(x, c, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: c
      class(*), intent(in) :: context

      real(wp) :: u, v

      associate (dummy => context)
      end associate
      u = (x(1) + 2.0_wp)**2 + x(2)**2 - 1.0_wp
      v = (x(1) - 2.0_wp)**2 + x(2)**2 - 1.0_wp
      c(1) = u * v
   end subroutine two_circle_constraint

   !> Constraint Jacobian: grad_h = v * grad_u + u * grad_v
   subroutine two_circle_constraint_grad(x, dc, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: dc
      class(*), intent(in) :: context

      real(wp) :: u, v

      associate (dummy => context)
      end associate
      u = (x(1) + 2.0_wp)**2 + x(2)**2 - 1.0_wp
      v = (x(1) - 2.0_wp)**2 + x(2)**2 - 1.0_wp
      ! grad_u = (2(x+2), 2y);  grad_v = (2(x-2), 2y)
      dc(1, 1) = v * 2.0_wp * (x(1) + 2.0_wp) + u * 2.0_wp * (x(1) - 2.0_wp)
      dc(1, 2) = v * 2.0_wp * x(2) + u * 2.0_wp * x(2)
   end subroutine two_circle_constraint_grad

   !> Test helper: is r on one of the two unit circles (within tol)?
   pure function on_two_circle_union(r) result(hit)
      real(wp), dimension(2), intent(in) :: r
      logical :: hit

      hit = on_left_circle(r) .or. on_right_circle(r)
   end function on_two_circle_union

   pure function on_left_circle(r) result(hit)
      real(wp), dimension(2), intent(in) :: r
      logical :: hit
      real(wp), parameter :: tol = 1.0e-4_wp

      hit = abs(sqrt((r(1) + 2.0_wp)**2 + r(2)**2) - 1.0_wp) < tol
   end function on_left_circle

   pure function on_right_circle(r) result(hit)
      real(wp), dimension(2), intent(in) :: r
      logical :: hit
      real(wp), parameter :: tol = 1.0e-4_wp

      hit = abs(sqrt((r(1) - 2.0_wp)**2 + r(2)**2) - 1.0_wp) < tol
   end function on_right_circle

   !------------------------------------------------------------------
   ! Newton-deflation helpers (scalar cubic problem)
   !------------------------------------------------------------------

   !> Residual F(x) = (x-1)(x+1)(x-3) = x^3 - 3x^2 - x + 3
   subroutine cubic_residual_ctx(x, f, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f
      class(*), intent(in) :: context

      associate (dummy => context)
      end associate
      f(1) = x(1)**3 - 3.0_wp*x(1)**2 - x(1) + 3.0_wp
   end subroutine cubic_residual_ctx

   !> Jacobian F'(x) = 3x^2 - 6x - 1
   subroutine cubic_jacobian_ctx(x, jac, context)
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: jac
      class(*), intent(in) :: context

      associate (dummy => context)
      end associate
      jac(1, 1) = 3.0_wp*x(1)**2 - 6.0_wp*x(1) - 1.0_wp
   end subroutine cubic_jacobian_ctx

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/slsqp (slsqp_module%slsqp_solver)
   !===========================================================================
   ! These exercise the raw SLSQP class directly; the wrapper layer above goes
   ! through new_slsqp_solver/solver%solve. Reference solutions are upstream.

   !> slsqp_test.f90: minimize Rosenbrock subject to x1^2 + x2^2 <= 1.
   subroutine test_slsqp_rosenbrock(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 2, m = 1, meq = 0, max_iter = 100
      real(wp), parameter :: acc = 1.0e-8_wp
      real(wp), parameter :: xl(n) = [-1.0_wp, -1.0_wp]
      real(wp), parameter :: xu(n) = [1.0_wp, 1.0_wp]
      type(slsqp_solver) :: solver
      real(wp) :: x(n)
      integer :: istat, iterations
      logical :: status_ok

      x = [0.1_wp, 0.1_wp]
      call solver%initialize(n, m, meq, max_iter, acc, slsqp_rosenbrock_func, &
                             slsqp_rosenbrock_grad, xl, xu, status_ok=status_ok, &
                             linesearch_mode=1, iprint=0)
      call check(error, status_ok, "slsqp initialize failed")
      if (allocated(error)) return

      call solver%optimize(x, istat, iterations)
      call check(error, istat, 0)
      if (allocated(error)) return
      call check(error, x(1), slsqp_rosen_x(1), thr=slsqp_kernel_thr)
      if (allocated(error)) return
      call check(error, x(2), slsqp_rosen_x(2), thr=slsqp_kernel_thr)

      call solver%destroy()
   end subroutine test_slsqp_rosenbrock

   !> slsqp_test_2.f90: minimize x1^2 + x2^2 + x3 subject to x1*x2 - x3 = 0
   !> (equality) and x3 - 1 >= 0 (inequality). Optimum x = [1,1,1], f = 3.
   subroutine test_slsqp_quadratic(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 3, m = 2, meq = 1, max_iter = 100
      real(wp), parameter :: acc = 1.0e-7_wp
      real(wp), parameter :: xl(n) = [-10.0_wp, -10.0_wp, -10.0_wp]
      real(wp), parameter :: xu(n) = [10.0_wp, 10.0_wp, 10.0_wp]
      type(slsqp_solver) :: solver
      real(wp) :: x(n), f
      integer :: istat, iterations
      logical :: status_ok

      x = [1.0_wp, 2.0_wp, 3.0_wp]
      call solver%initialize(n, m, meq, max_iter, acc, slsqp_quadratic_func, &
                             slsqp_quadratic_grad, xl, xu, status_ok=status_ok, &
                             linesearch_mode=1, iprint=0, alphamin=0.1_wp, alphamax=0.5_wp)
      call check(error, status_ok, "slsqp initialize failed")
      if (allocated(error)) return

      call solver%optimize(x, istat, iterations)
      call check(error, istat, 0)
      if (allocated(error)) return

      f = x(1)**2 + x(2)**2 + x(3)
      call check(error, f, 3.0_wp, thr=1.0e-5_wp)
      if (allocated(error)) return
      call check(error, x(1)*x(2) - x(3), 0.0_wp, thr=1.0e-5_wp)
      if (allocated(error)) return
      call check(error, x(3) >= 1.0_wp - slsqp_kernel_thr, "x3 >= 1 constraint violated")

      call solver%destroy()
   end subroutine test_slsqp_quadratic

   !> slsqp_test_3.f90: Rosenbrock solved with finite-difference gradients in
   !> all three modes (1=backward, 2=forward, 3=central).
   subroutine test_slsqp_fd_gradients(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 2, m = 1, meq = 0, max_iter = 100
      real(wp), parameter :: acc = 1.0e-8_wp
      real(wp), parameter :: gradient_delta = 1.0e-5_wp
      real(wp), parameter :: xl(n) = [-1.0_wp, -1.0_wp]
      real(wp), parameter :: xu(n) = [1.0_wp, 1.0_wp]
      type(slsqp_solver) :: solver
      real(wp) :: x(n)
      integer :: istat, iterations, gradient_mode
      logical :: status_ok
      character(len=64) :: msg

      do gradient_mode = 1, 3
         x = [0.1_wp, 0.1_wp]
         ! grad is passed to satisfy the argument but is never called because
         ! gradient_mode /= 0 selects finite differences.
         call solver%initialize(n, m, meq, max_iter, acc, slsqp_rosenbrock_func, &
                                slsqp_rosenbrock_grad, xl, xu, status_ok=status_ok, &
                                linesearch_mode=1, iprint=0, gradient_mode=gradient_mode, &
                                gradient_delta=gradient_delta)
         write (msg, '(a,i0)') "slsqp initialize failed, gradient_mode=", gradient_mode
         call check(error, status_ok, trim(msg))
         if (allocated(error)) return

         call solver%optimize(x, istat, iterations)
         write (msg, '(a,i0)') "slsqp did not converge, gradient_mode=", gradient_mode
         call check(error, istat == 0, trim(msg))
         if (allocated(error)) return

         write (msg, '(a,i0)') "x(1) mismatch, gradient_mode=", gradient_mode
         call check(error, x(1), slsqp_rosen_x(1), thr=1.0e-3_wp, message=trim(msg))
         if (allocated(error)) return
         write (msg, '(a,i0)') "x(2) mismatch, gradient_mode=", gradient_mode
         call check(error, x(2), slsqp_rosen_x(2), thr=1.0e-3_wp, message=trim(msg))
         if (allocated(error)) return

         call solver%destroy()
      end do
   end subroutine test_slsqp_fd_gradients

   !> slsqp_test_71.f90: Hock-Schittkowski problem 71, solved with both NNLS
   !> modes (1=nnls, 2=bvls). Optimum x = (1, 4.743, 3.821, 1.379, 0), f = 17.014.
   subroutine test_slsqp_hs71(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 5, m = 2, meq = 2, max_iter = 100
      real(wp), parameter :: acc = 1.0e-8_wp
      real(wp), parameter :: xl(n) = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      real(wp), parameter :: xu(n) = [5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 1.0e10_wp]
      real(wp), parameter :: x_ref(4) = [1.0_wp, 4.74299963_wp, 3.82114998_wp, 1.37940829_wp]
      real(wp), parameter :: f_ref = 17.0140173_wp
      type(slsqp_solver) :: solver
      real(wp) :: x(n), f
      integer :: istat, iterations, nnls_mode, i
      logical :: status_ok
      character(len=64) :: msg

      do nnls_mode = 1, 2
         x = [1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp, -24.0_wp]
         call solver%initialize(n, m, meq, max_iter, acc, slsqp_hs71_func, slsqp_hs71_grad, &
                                xl, xu, status_ok=status_ok, linesearch_mode=1, iprint=0, &
                                nnls_mode=nnls_mode)
         write (msg, '(a,i0)') "slsqp initialize failed, nnls_mode=", nnls_mode
         call check(error, status_ok, trim(msg))
         if (allocated(error)) return

         call solver%optimize(x, istat, iterations)
         write (msg, '(a,i0)') "slsqp did not converge, nnls_mode=", nnls_mode
         call check(error, istat == 0, trim(msg))
         if (allocated(error)) return

         do i = 1, 4
            write (msg, '(a,i0,a,i0)') "x(", i, ") mismatch, nnls_mode=", nnls_mode
            call check(error, x(i), x_ref(i), thr=1.0e-3_wp, message=trim(msg))
            if (allocated(error)) return
         end do

         f = x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3)
         write (msg, '(a,i0)') "objective mismatch, nnls_mode=", nnls_mode
         call check(error, f, f_ref, thr=1.0e-4_wp, message=trim(msg))
         if (allocated(error)) return

         call solver%destroy()
      end do
   end subroutine test_slsqp_hs71

   !> slsqp_test_stopping_criterion.f90: Rosenbrock with the constraint treated
   !> as equality (meq=1), missing bounds passed as NaN, all stop tols zero.
   ! NOTE: This test crashes if compiled with -ffpe-trap=invalid
   subroutine test_slsqp_stopping(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 2, m = 1, meq = 1, max_iter = 100
      real(wp), parameter :: acc = 1.0e-8_wp
      real(wp), parameter :: tolf = 0.0_wp, toldf = 0.0_wp, toldx = 0.0_wp
      type(slsqp_solver) :: solver
      real(wp) :: x(n), xl(n), xu(n), nan
      integer :: istat, iterations
      logical :: status_ok

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      xl = [-1.0_wp, nan]
      xu = [nan, 1.0_wp]

      x = [0.1_wp, 0.1_wp]
      call solver%initialize(n, m, meq, max_iter, acc, slsqp_rosenbrock_func, &
                             slsqp_rosenbrock_grad, xl, xu, status_ok=status_ok, &
                             linesearch_mode=1, iprint=0, tolf=tolf, toldf=toldf, toldx=toldx)
      call check(error, status_ok, "slsqp initialize failed")
      if (allocated(error)) return

      call solver%optimize(x, istat, iterations)
      call check(error, istat, 0)
      if (allocated(error)) return
      call check(error, x(1), slsqp_rosen_x(1), thr=slsqp_kernel_thr)
      if (allocated(error)) return
      call check(error, x(2), slsqp_rosen_x(2), thr=slsqp_kernel_thr)

      call solver%destroy()
   end subroutine test_slsqp_stopping

   !> Rosenbrock objective f = 100*(x2 - x1^2)^2 + (1 - x1)^2 with the
   !> constraint c = 1 - x1^2 - x2^2 (slsqp_module func interface).
   subroutine slsqp_rosenbrock_func(me, x, f, c)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      real(wp), dimension(:), intent(out) :: c

      f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2
      c(1) = 1.0_wp - x(1)**2 - x(2)**2
   end subroutine slsqp_rosenbrock_func

   !> Analytic gradients for slsqp_rosenbrock_func.
   subroutine slsqp_rosenbrock_grad(me, x, g, a)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g
      real(wp), dimension(:, :), intent(out) :: a

      g(1) = -400.0_wp*(x(2) - x(1)**2)*x(1) - 2.0_wp*(1.0_wp - x(1))
      g(2) = 200.0_wp*(x(2) - x(1)**2)
      a(1, 1) = -2.0_wp*x(1)
      a(1, 2) = -2.0_wp*x(2)
   end subroutine slsqp_rosenbrock_grad

   !> Objective f = x1^2 + x2^2 + x3 with c1 = x1*x2 - x3 (eq), c2 = x3 - 1 (ineq).
   subroutine slsqp_quadratic_func(me, x, f, c)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      real(wp), dimension(:), intent(out) :: c

      f = x(1)**2 + x(2)**2 + x(3)
      c(1) = x(1)*x(2) - x(3)
      c(2) = x(3) - 1.0_wp
   end subroutine slsqp_quadratic_func

   !> Analytic gradients for slsqp_quadratic_func.
   subroutine slsqp_quadratic_grad(me, x, g, a)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g
      real(wp), dimension(:, :), intent(out) :: a

      g(1) = 2.0_wp*x(1)
      g(2) = 2.0_wp*x(2)
      g(3) = 1.0_wp
      a(1, 1) = x(2)
      a(1, 2) = x(1)
      a(1, 3) = -1.0_wp
      a(2, 1) = 0.0_wp
      a(2, 2) = 0.0_wp
      a(2, 3) = 1.0_wp
   end subroutine slsqp_quadratic_grad

   !> Hock-Schittkowski problem 71 objective and constraints.
   subroutine slsqp_hs71_func(me, x, f, c)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      real(wp), dimension(:), intent(out) :: c

      f = x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3)
      c(1) = x(1)*x(2)*x(3)*x(4) - x(5) - 25.0_wp
      c(2) = x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2 - 40.0_wp
   end subroutine slsqp_hs71_func

   !> Analytic gradients for slsqp_hs71_func.
   subroutine slsqp_hs71_grad(me, x, g, a)
      class(slsqp_solver), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g
      real(wp), dimension(:, :), intent(out) :: a

      g(1) = x(4)*(2.0_wp*x(1) + x(2) + x(3))
      g(2) = x(1)*x(4)
      g(3) = x(1)*x(4) + 1.0_wp
      g(4) = x(1)*(x(1) + x(2) + x(3))
      g(5) = 0.0_wp

      a = 0.0_wp
      a(1, 1) = x(2)*x(3)*x(4)
      a(1, 2) = x(1)*x(3)*x(4)
      a(1, 3) = x(1)*x(2)*x(4)
      a(1, 4) = x(1)*x(2)*x(3)
      a(1, 5) = -1.0_wp
      a(2, 1) = 2.0_wp*x(1)
      a(2, 2) = 2.0_wp*x(2)
      a(2, 3) = 2.0_wp*x(3)
      a(2, 4) = 2.0_wp*x(4)
   end subroutine slsqp_hs71_grad

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/lbfgsb (setulb reverse communication)
   !===========================================================================
   ! All three drivers minimize the same bound-constrained 25-variable extended
   ! Rosenbrock problem (optimum f = 0); they differ only in termination policy.
   ! File/console output is suppressed via iprint = -1.

   !> driver1: run to the default convergence test (f -> 0).
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
   !> STOP (or early convergence) must fire and the objective must decrease.
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

   !> driver3: with a large evaluation budget the solver should reach f -> 0.
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
   !> projected gradient becomes negligible; otherwise only built-in tests stop.
   subroutine run_lbfgsb_kernel(max_nfg, f, final_task)
      integer, intent(in) :: max_nfg
      real(wp), intent(out) :: f
      character(len=60), intent(out) :: final_task

      integer, parameter :: nws = 2*lbfgsb_m*lbfgsb_n + 5*lbfgsb_n + 11*lbfgsb_m*lbfgsb_m + 8*lbfgsb_m
      real(wp) :: x(lbfgsb_n), l(lbfgsb_n), u(lbfgsb_n), g(lbfgsb_n), wa(nws), dsave(29)
      integer :: nbd(lbfgsb_n), iwa(3*lbfgsb_n), isave(44), i
      character(len=60) :: task, csave
      logical :: lsave(4)

      ! Odd variables: bounds [1, 100]; even variables: bounds [-100, 100].
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

      task = 'START'
      do while (task(1:2) == 'FG' .or. task(1:5) == 'NEW_X' .or. task(1:5) == 'START')
         call setulb(lbfgsb_n, lbfgsb_m, x, l, u, nbd, f, g, lbfgsb_factr, lbfgsb_pgtol, &
                     wa, iwa, task, -1, csave, lsave, isave, dsave)
         if (task(1:2) == 'FG') then
            call lbfgsb_rosenbrock_eval(x, f, g)
         else if (task(1:5) == 'NEW_X' .and. max_nfg > 0) then
            if (isave(34) >= max_nfg) &
               task = 'STOP: total f and g evaluation limit reached'
            if (dsave(13) <= 1.0e-10_wp*(1.0_wp + abs(f))) &
               task = 'STOP: projected gradient sufficiently small'
         end if
      end do

      final_task = task
   end subroutine run_lbfgsb_kernel

   !> Extended Rosenbrock objective and gradient (the L-BFGS-B sample problem).
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

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/fmin (1D derivative-free minimizer)
   !===========================================================================

   !> fmin test.f90: minimize sin(x) on [-4, 0]; the minimum is at x = -pi/2.
   subroutine test_fmin_sin(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: pi = acos(-1.0_wp)
      real(wp) :: xmin

      xmin = fmin(fmin_sin_func, -4.0_wp, 0.0_wp, fmin_tol)
      call check(error, xmin, -pi/2.0_wp, thr=fmin_thr)
   end subroutine test_fmin_sin

   !> Minimize (x - 2)^2 on [0, 5]; the minimum is at x = 2.
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

   !===========================================================================
   ! Raw-kernel ports: jacobwilliams/nlesolver-fortran (nlesolver_type)
   !===========================================================================
   ! Both solve f1 = x1^2 + x2 - 0.1 = 0, f2 = x2 + 0.2 = 0; root (sqrt(0.3),
   ! -0.2). Dense sweeps line-search/Broyden; sparse drives LSQR/LUSOL/LSMR.

   !> nlesolver_test_1.f90: dense solve across all four line-search step modes
   !> with the Broyden update both off and on, using scalar variable bounds.
   subroutine test_newton_dense(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: step_mode, b
      logical :: use_broyden
      character(len=64) :: label

      do step_mode = 1, 4
         do b = 0, 1
            use_broyden = (b == 1)
            write (label, '(a,i0,a,l1)') "dense step_mode=", step_mode, " broyden=", use_broyden
            call run_newton_dense(step_mode, use_broyden, trim(label), error)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_newton_dense

   !> sparse_test.f90: same line-search/Broyden sweep, run through each of the
   !> sparse linear-solver backends (LSQR, LUSOL, LSMR).
   subroutine test_newton_sparse(error)
      type(error_type), allocatable, intent(out) :: error

      integer :: modes(3), k, step_mode, b
      logical :: use_broyden
      character(len=16) :: names(3)
      character(len=80) :: label

      modes = [NLESOLVER_SPARSITY_LSQR, NLESOLVER_SPARSITY_LUSOL, NLESOLVER_SPARSITY_LSMR]
      names = ["LSQR ", "LUSOL", "LSMR "]

      do k = 1, 3
         do step_mode = 1, 4
            do b = 0, 1
               use_broyden = (b == 1)
               write (label, '(a,a,a,i0,a,l1)') "sparse ", trim(names(k)), &
                  " step_mode=", step_mode, " broyden=", use_broyden
               call run_newton_sparse(modes(k), step_mode, use_broyden, trim(label), error)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_newton_sparse

   !> Run one dense configuration and verify convergence to the root.
   subroutine run_newton_dense(step_mode, use_broyden, label, error)
      integer, intent(in) :: step_mode
      logical, intent(in) :: use_broyden
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      type(nlesolver_type) :: solver
      real(wp) :: x(newton_n)
      integer :: istat
      character(len=:), allocatable :: message

      call solver%initialize(n=newton_n, m=newton_m, max_iter=newton_max_iter, tol=newton_tol, &
                             func=newton_func, grad=newton_grad, step_mode=step_mode, &
                             use_broyden=use_broyden, n_intervals=2, fmin_tol=1.0e-2_wp, &
                             verbose=.false., bounds_mode=NLESOLVER_SCALAR_BOUNDS, &
                             xlow=[0.0_wp, -5.0_wp], xupp=[1.0_wp, 0.0_wp])
      call solver%status(istat, message)
      call check(error, istat == 0, label//": initialize failed: "//message)
      if (allocated(error)) return

      x = [1.0_wp, 2.0_wp]
      call solver%solve(x)
      call check_newton_root(x, label, error)

      call solver%destroy()
   end subroutine run_newton_dense

   !> Run one sparse configuration and verify convergence to the root.
   subroutine run_newton_sparse(sparsity_mode, step_mode, use_broyden, label, error)
      integer, intent(in) :: sparsity_mode
      integer, intent(in) :: step_mode
      logical, intent(in) :: use_broyden
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error

      type(nlesolver_type) :: solver
      real(wp) :: x(newton_n)
      integer :: istat
      character(len=:), allocatable :: message

      call solver%initialize(n=newton_n, m=newton_m, max_iter=newton_max_iter, tol=newton_tol, &
                             atol=newton_tol, btol=newton_tol, func=newton_func, &
                             grad_sparse=newton_grad_sparse, step_mode=step_mode, &
                             use_broyden=use_broyden, n_intervals=2, fmin_tol=1.0e-2_wp, &
                             verbose=.false., sparsity_mode=sparsity_mode, &
                             irow=newton_irow, icol=newton_icol, damp=0.0_wp)
      call solver%status(istat, message)
      call check(error, istat == 0, label//": initialize failed: "//message)
      if (allocated(error)) return

      x = [1.0_wp, 2.0_wp]
      call solver%solve(x)
      call check_newton_root(x, label, error)

      call solver%destroy()
   end subroutine run_newton_sparse

   !> Assert that x is the expected root and the residual norm is small.
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

   !> Residual vector of the test system (nlesolver_module func interface).
   subroutine newton_func(me, x, f)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: f

      f(1) = x(1)**2 + x(2) - 0.1_wp
      f(2) = x(2) + 0.2_wp
   end subroutine newton_func

   !> Dense Jacobian of the test system.
   subroutine newton_grad(me, x, g)
      class(nlesolver_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:, :), intent(out) :: g

      g(1, 1) = 2.0_wp*x(1)
      g(2, 1) = 0.0_wp
      g(1, 2) = 1.0_wp
      g(2, 2) = 1.0_wp
   end subroutine newton_grad

   !> Sparse Jacobian packed into the COO pattern (newton_irow, newton_icol).
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

end module test_math_solvers
