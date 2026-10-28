!> Settings of the 3D MOZ (UV) solve
!>
!> JSON/TOML keys: closure, correction, finite_size, susceptibility, convergence.*, solver.* with
!> the solver controls as a nested group plus solver.jfnk.precond,
!> and solvent.* with the settings of the 1D VV solve of the model's solvent copy
!> Verbosity and thread count come from the run context
module moist_model_moz_3d_parameters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_solver, only: moist_model_moz_solver_parameters_type, check_preconditioner_name
   use moist_model_moz_solvent_vv_parameters, only: solvent_vv_parameters_type
   use moist_model_moz_thermodynamics, only: check_correction_name, check_finite_size_name
   implicit none(type, external)
   private

   public :: model_moz_3d_parameters_type, preconditioners_3d_uv, susceptibility_modes, cold_guess_modes

   !> Preconditioners the 3D UV solve implements
   character(len=*), parameter :: preconditioners_3d_uv(1) = [character(len=8) :: "none"]
   !> Susceptibility modes: OZ reconstruction, direct transform, or linear interpolation
   character(len=*), parameter :: susceptibility_modes(3) = [character(len=8) :: "oz", "direct", "linear"]
   !> Cold starts: map c_s = 0 through OZ, or start with t_s = 0
   character(len=*), parameter :: cold_guess_modes(2) = [character(len=8) :: "cs0", "ts0"]

   !> Settings of the 3D UV solve and of the solvent's VV solve
   type, extends(moist_model_parameters_type) :: model_moz_3d_parameters_type
      !> Closure: HNC, KH or PSE2-4
      character(len=8) :: closure = "KH"
      !> Pressure correction added to the solvation free energy: none, PC or PC+
      character(len=8) :: correction = "none"
      !> Independent leading electrostatic finite-size term: none or huenenberger
      character(len=16) :: finite_size = "none"
      !> Solvent response on the grid's k-points
      !>
      !> - oz reconstructs it from cubic short-range c(k)
      !> - direct and linear sample the exported h-based response
      character(len=8) :: susceptibility = "oz"
      !> RMS residual below which the UV solve is converged
      real(wp) :: tolerance = 1.0e-7_wp
      !> Outer iterations of the UV solve
      integer :: maxiter = 128
      !> Cold solves ramp the closure up from KH: KH, PSE2, PSE3, PSE4, HNC
      logical :: closure_ramp = .true.
      !> RMS residual below which the ramp steps to the next closure
      real(wp) :: closure_ramp_rmse = 1.0e-5_wp
      !> Use the previous converged UV fields as an initial guess at fixed geometry
      logical :: reuse_initial_guess = .true.
      !> Cold initial field: cs0 maps c_s = 0; ts0 sets t_s = 0 directly
      character(len=8) :: cold_guess = "cs0"
      !> Preconditioner of the Newton-Krylov solves: none
      character(len=8) :: precond = "none"
      !> Solver controls of the UV solve
      type(moist_model_moz_solver_parameters_type) :: solver
      !> Settings of the VV solve of the model's solvent copy
      type(solvent_vv_parameters_type) :: solvent
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_3d_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => register_3d_entries
      !> Check names and ranges, including the solver and solvent settings
      procedure :: validate => validate_3d_parameters
   end type model_moz_3d_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_3d_defaults(self)
      !> Parameter values
      class(model_moz_3d_parameters_type), intent(inout) :: self
      type(model_moz_3d_parameters_type) :: defaults

      self%closure = defaults%closure
      self%correction = defaults%correction
      self%finite_size = defaults%finite_size
      self%susceptibility = defaults%susceptibility
      self%tolerance = defaults%tolerance
      self%maxiter = defaults%maxiter
      self%closure_ramp = defaults%closure_ramp
      self%closure_ramp_rmse = defaults%closure_ramp_rmse
      self%reuse_initial_guess = defaults%reuse_initial_guess
      self%cold_guess = defaults%cold_guess
      self%precond = defaults%precond
      call self%solver%init_defaults()
      call self%solvent%init_defaults()
   end subroutine init_3d_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_3d_entries(self)
      !> Parameter values
      class(model_moz_3d_parameters_type), intent(inout), target :: self

      call self%register_string("closure", self%closure, level=2)
      call self%register_string("correction", self%correction, level=2)
      call self%register_string("finite_size", self%finite_size, level=2)
      call self%register_string("susceptibility", self%susceptibility)
      call self%register_real_scalar("convergence.tolerance", self%tolerance, level=2)
      call self%register_int_scalar("convergence.maxiter", self%maxiter, level=2)
      call self%register_logical("convergence.closure_ramp", self%closure_ramp)
      call self%register_real_scalar("convergence.closure_ramp_rmse", self%closure_ramp_rmse)
      call self%register_logical("convergence.reuse_initial_guess", self%reuse_initial_guess)
      call self%register_string("convergence.cold_guess", self%cold_guess)
      call self%register_group("solver", self%solver)
      ! Preconditioner settings print with the JFNK controls
      call self%register_string("solver.jfnk.precond", self%precond, &
         & level=self%solver%group_level("jfnk", 3))
      call self%register_group("solvent", self%solvent)
   end subroutine register_3d_entries

   !> Check names and ranges, including the solver and solvent settings
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error First invalid setting
   subroutine validate_3d_parameters(self, error)
      !> Parameter values
      class(model_moz_3d_parameters_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure

      call new_closure(closure, self%closure, error)
      if (allocated(error)) return
      call check_correction_name(self%correction, error)
      if (allocated(error)) return
      call check_finite_size_name(self%finite_size, error)
      if (allocated(error)) return
      if (.not. any(to_lower(trim(adjustl(self%cold_guess))) == cold_guess_modes)) then
         call fatal_error(error, "Unknown 3D MOZ cold guess: "//trim(self%cold_guess)//" (choose cs0 or ts0)")
         return
      end if
      if (.not. any(to_lower(trim(adjustl(self%susceptibility))) == susceptibility_modes)) then
         call fatal_error(error, "Unknown 3D MOZ susceptibility mode: "//trim(self%susceptibility)// &
            & " (choose oz, direct or linear)")
         return
      end if
      if (to_lower(trim(adjustl(self%finite_size))) /= "none" .and. &
         & to_lower(trim(adjustl(self%susceptibility))) /= "oz") then
         call fatal_error(error, "huenenberger correction requires the OZ solvent susceptibility")
         return
      end if
      call check_preconditioner_name(self%precond, preconditioners_3d_uv, "3D UV", error)
      if (allocated(error)) return
      if (self%solver%runs("jfnk") .and. self%solver%jfnk%deflation) then
         call fatal_error(error, "jfnk.deflation is only available for the 1D VV solve")
         return
      end if
      if (.not. positive(self%tolerance) .or. self%maxiter < 1 .or. .not. positive(self%closure_ramp_rmse)) then
         call fatal_error(error, "3D MOZ convergence needs a positive tolerance, iteration limit and ramp residual")
      else
         call self%solver%validate(error)
         if (.not. allocated(error)) call self%solvent%validate(error)
      end if
   end subroutine validate_3d_parameters

   !> Finite positive value
   !>
   !> @param[in] value Scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

end module moist_model_moz_3d_parameters
