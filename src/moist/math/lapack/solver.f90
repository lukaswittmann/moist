!> Provides a wrapper for the eigenvalue solvers provided by LAPACK

!> LAPACK based eigenvalue solvers
module moist_math_lapack_solver
   use moist_math_lapack_type, only: solver_type, context_solver
   use moist_math_lapack_sygvd, only: sygvd_solver, new_sygvd
   use moist_math_lapack_sygvr, only: sygvr_solver, new_sygvr
   use mctc_env, only: sp, dp, error_type
   implicit none
   private

   public :: solver_type, lapack_algorithm

   !> Possible solvers provided by LAPACK
   type :: enum_lapack
      !> Divide-and-conquer solver
      integer :: gvd = 1
      !> Relatively robust solver
      integer :: gvr = 2
   end type enum_lapack

   !> Actual enumerator of possible solvers
   type(enum_lapack), parameter :: lapack_algorithm = enum_lapack()

   !> Generator for LAPACK based electronic solvers
   type, public, extends(context_solver) :: lapack_solver
      !> Selected electronic solver algorithm
      integer :: algorithm = lapack_algorithm%gvd
   contains
      !> Create new instance of electronic solver
      procedure :: new
      !> Delete an electronic solver instance
      procedure :: delete
   end type lapack_solver

contains

!> Create new electronic solver
   subroutine new(self, solver, ndim)
      !> Instance of the solver factory
      class(lapack_solver), intent(inout) :: self
      !> New electronic solver
      class(solver_type), allocatable, intent(out) :: solver
      !> Dimension of the eigenvalue problem
      integer, intent(in) :: ndim

      select case (self%algorithm)
      case (lapack_algorithm%gvd)
         block
            type(sygvd_solver), allocatable :: tmp
            allocate (tmp)
            call new_sygvd(tmp, ndim)
            call move_alloc(tmp, solver)
         end block
      case (lapack_algorithm%gvr)
         block
            type(sygvr_solver), allocatable :: tmp
            allocate (tmp)
            call new_sygvr(tmp, ndim)
            call move_alloc(tmp, solver)
         end block
      end select
   end subroutine new

!> Delete electronic solver instance
   subroutine delete(self, solver)
      !> Instance of the solver factory
      class(lapack_solver), intent(inout) :: self
      !> Electronic solver instance
      class(solver_type), allocatable, intent(inout) :: solver

      if (allocated(solver)) deallocate (solver)
   end subroutine delete

end module moist_math_lapack_solver
