!> Settings of the 1D MOZ (UV) solve
!>
!> JSON/TOML keys: closure, convergence.*, damping.*, solver.* with the solver
!> controls as a nested group plus solver.jfnk.precond*, and
!> solvent.* with the settings of the 1D VV solve of the model's solvent copy
!> Verbosity and thread count come from the run context
module moist_model_moz_1d_parameters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_solver, only: moist_model_moz_solver_parameters_type, check_preconditioner_name
   use moist_model_moz_1d_preconditioner, only: preconditioners_1d_uv
   use moist_model_moz_solvent_vv_parameters, only: solvent_vv_parameters_type
   implicit none(type, external)
   private

   public :: model_moz_1d_parameters_type

   !> Settings of the 1D UV solve and of the solvent's VV solve
   type, extends(moist_model_parameters_type) :: model_moz_1d_parameters_type
      !> Closure: HNC, KH or PSE2-4
      character(len=8) :: closure = "HNC"
      !> RMS residual below which the UV solve is converged
      real(wp) :: tolerance = 1.0e-7_wp
      !> Outer iterations of the UV solve
      integer :: maxiter = 128
      !> Cold solves ramp the closure up from KH: KH, PSE2, PSE3, then PSE4 or HNC
      logical :: closure_ramp = .true.
      !> RMS residual below which the ramp steps to the next closure
      real(wp) :: closure_ramp_rmse = 1.0e-5_wp
      !> Reuse converged UV fields across host updates at fixed geometry
      logical :: reuse_initial_guess = .true.
      !> Small-k damping of updates the solver allows to post-process
      logical :: damp_small_k = .true.
      !> Wave number and k = 0 weight of the damping, 1/bohr
      real(wp) :: damp_small_k_k0 = 1.0_wp, damp_small_k_min_weight = 0.15_wp
      !> Preconditioner of the Newton-Krylov solves: none or spectral
      character(len=8) :: precond = "spectral"
      !> Wave number, 1/bohr, and k = 0 weight of the spectral filter
      real(wp) :: precond_k0 = 1.0_wp, precond_w0 = 0.1_wp
      !> Solver controls of the UV solve
      type(moist_model_moz_solver_parameters_type) :: solver
      !> Settings of the VV solve of the model's solvent copy
      type(solvent_vv_parameters_type) :: solvent
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_1d_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => register_1d_entries
      !> Check names and ranges, including the solver and solvent settings
      procedure :: validate => validate_1d_parameters
   end type model_moz_1d_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_1d_defaults(self)
      !> Parameter values
      class(model_moz_1d_parameters_type), intent(inout) :: self
      type(model_moz_1d_parameters_type) :: defaults

      self%closure = defaults%closure
      self%tolerance = defaults%tolerance
      self%maxiter = defaults%maxiter
      self%closure_ramp = defaults%closure_ramp
      self%closure_ramp_rmse = defaults%closure_ramp_rmse
      self%reuse_initial_guess = defaults%reuse_initial_guess
      self%damp_small_k = defaults%damp_small_k
      self%damp_small_k_k0 = defaults%damp_small_k_k0
      self%damp_small_k_min_weight = defaults%damp_small_k_min_weight
      self%precond = defaults%precond
      self%precond_k0 = defaults%precond_k0
      self%precond_w0 = defaults%precond_w0
      call self%solver%init_defaults()
      call self%solvent%init_defaults()
   end subroutine init_1d_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_1d_entries(self)
      !> Parameter values
      class(model_moz_1d_parameters_type), intent(inout), target :: self

      !> Print levels: damping acts on Picard and DIIS updates; preconditioner settings print with JFNK
      integer :: damping_level, precond_level, filter_level

      damping_level = self%solver%group_level("picard", 3)
      precond_level = self%solver%group_level("jfnk", 2)
      filter_level = self%solver%group_level("jfnk", 3)
      call self%register_string("closure", self%closure, level=2)
      call self%register_real_scalar("convergence.tolerance", self%tolerance, level=2)
      call self%register_int_scalar("convergence.maxiter", self%maxiter, level=2)
      call self%register_logical("convergence.closure_ramp", self%closure_ramp)
      call self%register_real_scalar("convergence.closure_ramp_rmse", self%closure_ramp_rmse)
      call self%register_logical("convergence.reuse_initial_guess", self%reuse_initial_guess)
      call self%register_logical("damping.small_k", self%damp_small_k, level=damping_level)
      call self%register_real_scalar("damping.small_k_k0", self%damp_small_k_k0, level=damping_level)
      call self%register_real_scalar("damping.small_k_min_weight", self%damp_small_k_min_weight, level=damping_level)
      call self%register_group("solver", self%solver)
      call self%register_string("solver.jfnk.precond", self%precond, level=precond_level)
      call self%register_real_scalar("solver.jfnk.precond_k0", self%precond_k0, level=filter_level)
      call self%register_real_scalar("solver.jfnk.precond_w0", self%precond_w0, level=filter_level)
      call self%register_group("solvent", self%solvent)
   end subroutine register_1d_entries

   !> Check names and ranges, including the solver and solvent settings
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error First invalid setting
   subroutine validate_1d_parameters(self, error)
      !> Parameter values
      class(model_moz_1d_parameters_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure

      call new_closure(closure, self%closure, error)
      if (allocated(error)) return
      call check_preconditioner_name(self%precond, preconditioners_1d_uv, "1D UV", error)
      if (allocated(error)) return
      if (self%solver%runs("jfnk") .and. self%solver%jfnk%deflation) then
         call fatal_error(error, "jfnk.deflation is only available for the 1D VV solve")
         return
      end if
      if (.not. positive(self%tolerance) .or. self%maxiter < 1 .or. .not. positive(self%closure_ramp_rmse)) then
         call fatal_error(error, "1D MOZ convergence needs a positive tolerance, iteration limit and ramp residual")
      else if (.not. positive(self%damp_small_k_k0) .or. .not. weight(self%damp_small_k_min_weight)) then
         call fatal_error(error, "1D MOZ damping needs a positive k0 and a k = 0 weight in [0, 1]")
      else if (.not. positive(self%precond_k0) .or. .not. weight(self%precond_w0)) then
         call fatal_error(error, "1D MOZ spectral filter needs a positive k0 and a k = 0 weight in [0, 1]")
      else
         call self%solver%validate(error)
         if (.not. allocated(error)) call self%solvent%validate(error)
      end if
   end subroutine validate_1d_parameters

   !> Finite positive value
   !>
   !> @param[in] value Scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

   !> Finite weight in [0, 1]
   !>
   !> @param[in] value Weight to check
   pure logical function weight(value)
      !> Weight to check
      real(wp), intent(in) :: value
      weight = ieee_is_finite(value) .and. value >= 0.0_wp .and. value <= 1.0_wp
   end function weight

end module moist_model_moz_1d_parameters
