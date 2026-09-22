!> Shared abstract interface for nonlinear solvers
module moist_math_solver_type
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type
   implicit none(type, external)
   private
   public :: solver_base_type

   !> Abstract base type for nonlinear solvers
   type, abstract :: solver_base_type
   contains
      !> Solve the problem starting from initial guess
      procedure(solve_solver), deferred :: solve

      !> Clean up solver resources
      procedure(destroy_solver), deferred :: destroy
   end type solver_base_type

   !> Deferred solver operations
   abstract interface

      !> Solve the system or optimization problem
      !>
      !> @param[in,out] self Solver instance
      !> @param[in,out] x Initial guess, replaced by the solution
      !> @param[out] error Solver failure
      subroutine solve_solver(self, x, error)
         import :: solver_base_type, wp, error_type
         implicit none(type, external)
         !> Solver instance
         class(solver_base_type), intent(inout), target :: self
         !> Initial guess, replaced by the solution
         real(wp), intent(inout) :: x(:)
         !> Solver failure
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_solver

      !> Destroy solver and free resources
      !>
      !> @param[in,out] self Solver instance
      subroutine destroy_solver(self)
         import :: solver_base_type
         implicit none(type, external)
         !> Solver instance
         class(solver_base_type), intent(inout), target :: self
      end subroutine destroy_solver
   end interface

end module moist_math_solver_type
