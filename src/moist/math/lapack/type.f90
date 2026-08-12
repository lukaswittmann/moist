module moist_math_lapack_type
   use mctc_env, only: sp, dp, error_type
   implicit none
   private

   public :: eigen_solver_type, context_solver

   !> Abstract base class for general solvers
   type, abstract :: eigen_solver_type
   contains
      generic :: solve => solve_sp, solve_dp
      procedure(solve_sp), deferred :: solve_sp
      procedure(solve_dp), deferred :: solve_dp
   end type eigen_solver_type

   abstract interface
      subroutine solve_sp(self, hmat, smat, eval, error)
         import :: eigen_solver_type, error_type, sp
         class(eigen_solver_type), intent(inout) :: self
         real(sp), contiguous, intent(inout) :: hmat(:, :)
         real(sp), contiguous, intent(in) :: smat(:, :)
         real(sp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_sp
      subroutine solve_dp(self, hmat, smat, eval, error)
         import :: eigen_solver_type, error_type, dp
         class(eigen_solver_type), intent(inout) :: self
         real(dp), contiguous, intent(inout) :: hmat(:, :)
         real(dp), contiguous, intent(in) :: smat(:, :)
         real(dp), contiguous, intent(inout) :: eval(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_dp
   end interface

   !> Abstract base class for creating solver instances
   type, abstract :: context_solver
   contains
      !> Create new instance of solver
      procedure(new), deferred :: new
      !> Delete an solver instance
      procedure(delete), deferred :: delete
   end type context_solver

   abstract interface
      !> Create new solver
      subroutine new(self, solver, ndim)
         import :: context_solver, eigen_solver_type
         !> Instance of the solver factory
         class(context_solver), intent(inout) :: self
         !> New solver
         class(eigen_solver_type), allocatable, intent(out) :: solver
         !> Dimension of the eigenvalue problem
         integer, intent(in) :: ndim
      end subroutine new

      !> Delete solver instance
      subroutine delete(self, solver)
         import :: context_solver, eigen_solver_type
         !> Instance of the solver factory
         class(context_solver), intent(inout) :: self
         !> Solver instance
         class(eigen_solver_type), allocatable, intent(inout) :: solver
      end subroutine delete
   end interface

end module moist_math_lapack_type
