!> 1D MOZ preconditioners: radial spectral filter and block k-space inverse
!>
!> The model owns both operators; `select_1d_preconditioner` points at the one
!> lent to the solver, after the name was checked with `check_preconditioner_name`
!> against the model's list
module moist_model_moz_1d_preconditioner
   use mctc_io_utils, only: to_lower
   use moist_model_moz_solver_preconditioner, only: moist_model_moz_solver_preconditioner_type
   use moist_model_moz_1d_preconditioner_spectral, only: moist_model_moz_1d_preconditioner_spectral_type
   use moist_model_moz_1d_preconditioner_kspace, only: moist_model_moz_1d_preconditioner_kspace_type
   implicit none(type, external)
   private
   public :: moist_model_moz_1d_preconditioner_spectral_type, moist_model_moz_1d_preconditioner_kspace_type
   public :: select_1d_preconditioner, preconditioners_1d_uv, preconditioners_1d_vv

   !> Names the 1D UV model implements
   character(len=*), parameter :: preconditioners_1d_uv(2) = [character(len=8) :: "none", "spectral"]
   !> Names the 1D VV model implements
   character(len=*), parameter :: preconditioners_1d_vv(3) = [character(len=8) :: "none", "spectral", "kspace"]

contains

   !> Point at the operator selected by name; "none" leaves the pointer null
   !>
   !> @param[in]    name     Configured name, already checked against the model's list
   !> @param[inout] spectral Model-owned spectral filter, already set up
   !> @param[out]   precond  Operator to lend to the solver, null for "none"
   !> @param[inout] kspace   Model-owned block inverse (VV only), bound to `spectral` here
   subroutine select_1d_preconditioner(name, spectral, precond, kspace)
      !> Configured name
      character(len=*), intent(in) :: name
      !> Model-owned spectral filter
      type(moist_model_moz_1d_preconditioner_spectral_type), target, intent(inout) :: spectral
      !> Operator to lend, null for "none"
      class(moist_model_moz_solver_preconditioner_type), pointer, intent(out) :: precond
      !> Model-owned block inverse
      type(moist_model_moz_1d_preconditioner_kspace_type), target, optional, intent(inout) :: kspace

      precond => null()
      select case (to_lower(trim(adjustl(name))))
      case ("spectral")
         precond => spectral
      case ("kspace")
         if (.not. present(kspace)) return
         call kspace%bind(spectral)
         precond => kspace
      case default
         ! "none" and unchecked names lend no operator
      end select
   end subroutine select_1d_preconditioner

end module moist_model_moz_1d_preconditioner
