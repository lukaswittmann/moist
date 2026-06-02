module moist_math_solver
   ! Re-export solver types and functions from submodules
   use moist_math_solver_newton, only: &
      & moist_math_solver_newton_type, &
      & new_newton_solver, &
      & newton_linear_solver_dense, &
      & newton_linear_solver_lsqr, &
      & newton_linear_solver_lusol, &
      & newton_linear_solver_lsmr, &
      & newton_linesearch_simple, &
      & newton_linesearch_backtrack, &
      & newton_linesearch_exact, &
      & newton_linesearch_fixedpoint, &
      & newton_bounds_ignore, &
      & newton_bounds_scalar, &
      & newton_bounds_vector, &
      & newton_bounds_wall

   use moist_math_solver_slsqp, only: &
      & moist_math_solver_slsqp_type, &
      & new_slsqp_solver

   use moist_math_solver_lbfgsb, only: &
      & moist_math_solver_lbfgsb_type, &
      & new_lbfgsb_solver

   use moist_math_lapack_solver, only: &
      & solver_type, &
      & lapack_algorithm, &
      & lapack_solver

   implicit none
   private

   ! Re-export Newton solver
   public :: moist_math_solver_newton_type
   public :: new_newton_solver
   public :: newton_linear_solver_dense, newton_linear_solver_lsqr
   public :: newton_linear_solver_lusol, newton_linear_solver_lsmr
   public :: newton_linesearch_simple, newton_linesearch_backtrack
   public :: newton_linesearch_exact, newton_linesearch_fixedpoint
   public :: newton_bounds_ignore, newton_bounds_scalar
   public :: newton_bounds_vector, newton_bounds_wall

   ! Re-export SLSQP solver
   public :: moist_math_solver_slsqp_type
   public :: new_slsqp_solver

   ! Re-export L-BFGS-B solver
   public :: moist_math_solver_lbfgsb_type
   public :: new_lbfgsb_solver

   ! Re-export LAPACK solver
   public :: solver_type
   public :: lapack_algorithm
   public :: lapack_solver

contains

end module moist_math_solver
