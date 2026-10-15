!> Abstract left preconditioner for Newton-Krylov linear solves
!>
!> Model drivers own the operators and lend them to the solver. An operator
!> used inside GMRES must be linear, nonsingular and fixed throughout the
!> inner solve, including restarts. Scratch and counters may change.
!> Fallback filters share this interface but act outside the inner solve.
!>
!> Concrete 1D operators live in model/moz/1d/preconditioner; TODO: port the 3D operators.
!> A coarse-grid Jacobian implementation needs transfers between the grids,
!> a fixed linear coarse solve and treatment of the remaining fine-grid modes.
module moist_model_moz_solver_preconditioner
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   implicit none(type, external)
   private
   public :: moist_model_moz_solver_preconditioner_type, check_preconditioner_name

   !> Left preconditioner lent to a Krylov solver for the duration of a solve
   type, abstract :: moist_model_moz_solver_preconditioner_type
   contains
      procedure(ready_interface), deferred :: ready
      procedure(apply_interface), deferred :: apply
      procedure(describe_interface), deferred :: describe
   end type moist_model_moz_solver_preconditioner_type

   abstract interface
      !> Whether the operator is built and should be used for the next linear solve
      !> Called once per Newton update; building and refreshing belong to the model.
      !>
      !> @param[in] self  Preconditioner
      !> @param[in] rmse  Previous input residual RMS; huge(wp) at the start of a stage
      !> @return True when `apply` should precondition the linear solve
      logical function ready_interface(self, rmse) result(usable)
         import :: moist_model_moz_solver_preconditioner_type, wp
         implicit none(type, external)
         !> Preconditioner
         class(moist_model_moz_solver_preconditioner_type), intent(in) :: self
         !> Previous input residual RMS; huge(wp) at the start of a stage
         real(wp), intent(in) :: rmse
      end function ready_interface

      !> Replace vflat by M^{-1} vflat in place, producing finite values
      !>
      !> Flattening uses Fortran column order: i + (pair - 1)*npts.
      !> The model selects the field representation; this interface adds no
      !> potential splitting or grid conversion. Do not retain a pointer to vflat.
      !>
      !> @param[inout] self  Preconditioner and its scratch
      !> @param[inout] vflat Flat Krylov vector
      subroutine apply_interface(self, vflat)
         import :: moist_model_moz_solver_preconditioner_type, wp
         implicit none(type, external)
         !> Preconditioner and its scratch
         class(moist_model_moz_solver_preconditioner_type), intent(inout) :: self
         !> Flat Krylov vector
         real(wp), contiguous, target, intent(inout) :: vflat(:)
      end subroutine apply_interface

      !> Status for the solver notes: the operator `apply` uses at this residual,
      !> or why none is used; the solver reports every change of it
      !>
      !> @param[in] self Preconditioner
      !> @param[in] rmse residual RMS of the previous Newton step
      function describe_interface(self, rmse) result(text)
         import :: moist_model_moz_solver_preconditioner_type, wp
         implicit none(type, external)
         !> Preconditioner
         class(moist_model_moz_solver_preconditioner_type), intent(in) :: self
         !> Residual RMS of the previous Newton step
         real(wp), intent(in) :: rmse
         !> Status text, e.g. "jacobi" or "kspace, off until RMSR < 1.0E-03"
         character(len=:), allocatable :: text
      end function describe_interface
   end interface

contains

   !> Reject a preconditioner name the model does not implement
   !>
   !> @param[in]  name    requested name, case-insensitive
   !> @param[in]  allowed lowercase names the model implements
   !> @param[in]  label   model label for the message, e.g. "1D VV"
   !> @param[out] error   unknown or unsupported name
   subroutine check_preconditioner_name(name, allowed, label, error)
      !> Requested name
      character(len=*), intent(in) :: name
      !> Names the model implements
      character(len=*), intent(in) :: allowed(:)
      !> Model label
      character(len=*), intent(in) :: label
      !> Unknown or unsupported name
      type(error_type), allocatable, intent(out) :: error
      !> Comma-separated list of the allowed names
      character(len=:), allocatable :: choices
      !> Allowed-name counter
      integer :: i

      if (any(allowed == to_lower(trim(adjustl(name))))) return
      if (size(allowed) == 0) then
         call fatal_error(error, "No preconditioners are available for "//label)
         return
      end if
      choices = trim(allowed(1))
      do i = 2, size(allowed)
         choices = choices//", "//trim(allowed(i))
      end do
      call fatal_error(error, "Preconditioner '"//trim(adjustl(name))//"' is not available for "// &
         & label//"; choose one of: "//choices)
   end subroutine check_preconditioner_name

end module moist_model_moz_solver_preconditioner
