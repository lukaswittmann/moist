!> Controls of the transpose (adjoint) Krylov solves of MOZ response calculations
!>
!> - Below the solver base type, so every solver can carry them
!> - Read from the jfnk.adjoint group of the solver parameters
module moist_model_moz_solver_core_adjoint
   use mctc_env, only: error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_jfnk_adjoint_parameters_type

   !> Transpose Krylov controls, independent of the nonlinear solve
   !>
   !> JSON/TOML keys: restart, max_krylov
   type, extends(moist_model_parameters_type) :: moist_model_moz_solver_jfnk_adjoint_parameters_type
      !> Transpose GMRES restart length
      integer :: restart = 50
      !> Inner GMRES iterations per transpose refinement
      integer :: max_krylov = 50
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => adjoint_init_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => adjoint_register_entries
      !> Check the Krylov limits
      procedure :: validate => adjoint_validate
   end type moist_model_moz_solver_jfnk_adjoint_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine adjoint_init_defaults(self)
      !> Parameter values
      class(moist_model_moz_solver_jfnk_adjoint_parameters_type), intent(inout) :: self
      type(moist_model_moz_solver_jfnk_adjoint_parameters_type) :: defaults

      self%restart = defaults%restart
      self%max_krylov = defaults%max_krylov
   end subroutine adjoint_init_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine adjoint_register_entries(self)
      !> Parameter values
      class(moist_model_moz_solver_jfnk_adjoint_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("restart", self%restart)
      call self%register_int_scalar("max_krylov", self%max_krylov)
   end subroutine adjoint_register_entries

   !> Check the Krylov limits
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error Invalid limits
   subroutine adjoint_validate(self, error)
      !> Parameter values
      class(moist_model_moz_solver_jfnk_adjoint_parameters_type), intent(inout) :: self
      !> Invalid limits
      type(error_type), allocatable, intent(out) :: error

      if (self%restart < 1 .or. self%max_krylov < 1) then
         call fatal_error(error, "Adjoint needs positive Krylov limits")
      end if
   end subroutine adjoint_validate

end module moist_model_moz_solver_core_adjoint
