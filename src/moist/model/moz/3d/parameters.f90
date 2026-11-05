!> Settings of the 3D MOZ (UV) solve
!>
!> JSON/TOML keys: closure, correction, finite_size, susceptibility, convergence.*, solver.* with
!> the solver controls as a nested group plus solver.jfnk.precond and its precond_* settings,
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
   character(len=*), parameter :: preconditioners_3d_uv(2) = [character(len=16) :: "none", "cartesian-aux"]
   !> Susceptibility modes: OZ reconstruction, direct transform, or linear interpolation
   character(len=*), parameter :: susceptibility_modes(3) = [character(len=8) :: "oz", "direct", "linear"]
   !> Cold starts: map c_s = 0 through OZ, or start with t_s = 0
   character(len=*), parameter :: cold_guess_modes(2) = [character(len=8) :: "cs0", "ts0"]

   !> Settings of the 3D UV solve and of the solvent's VV solve
   type, extends(moist_model_parameters_type) :: model_moz_3d_parameters_type
      !> Closure: HNC, KH or PSE2-4
      character(len=8) :: closure = "HNC"
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
      real(wp) :: tolerance = 1.0e-9_wp
      !> Outer iterations of the UV solve
      integer :: maxiter = 128
      !> Cold solves ramp the closure up from KH: KH, PSE2, PSE3, then PSE4 or HNC
      logical :: closure_ramp = .true.
      !> RMS residual below which the ramp steps to the next closure
      real(wp) :: closure_ramp_rmse = 1.0e-5_wp
      !> Use the previous converged UV fields as an initial guess at fixed geometry
      logical :: reuse_initial_guess = .true.
      !> Cold initial field: cs0 maps c_s = 0; ts0 sets t_s = 0 directly
      character(len=8) :: cold_guess = "cs0"
      !> Preconditioner of the Newton-Krylov solves: none or cartesian-aux
      character(len=16) :: precond = "none"
      !> cartesian-aux: requested cells per axis of the auxiliary lattice, even, at least 4, capped by the grid
      integer :: precond_aux_n = 64
      !> cartesian-aux: precond_aux_n/2 adds a coarse lattice of half the capped fine cells, precond_aux_n none
      integer :: precond_coarse_n = 32
      !> cartesian-aux: most Krylov steps of each lattice solve, at least 1
      integer :: precond_inner_max_iter = 2
      !> cartesian-aux: residual reduction that ends a lattice solve, in (0, 1)
      real(wp) :: precond_inner_rtol = 0.1_wp
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
      self%precond_aux_n = defaults%precond_aux_n
      self%precond_coarse_n = defaults%precond_coarse_n
      self%precond_inner_max_iter = defaults%precond_inner_max_iter
      self%precond_inner_rtol = defaults%precond_inner_rtol
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
      call self%register_int_scalar("solver.jfnk.precond_aux_n", self%precond_aux_n, &
         & level=self%solver%group_level("jfnk", 3))
      call self%register_int_scalar("solver.jfnk.precond_coarse_n", self%precond_coarse_n, &
         & level=self%solver%group_level("jfnk", 3))
      call self%register_int_scalar("solver.jfnk.precond_inner_max_iter", self%precond_inner_max_iter, &
         & level=self%solver%group_level("jfnk", 3))
      call self%register_real_scalar("solver.jfnk.precond_inner_rtol", self%precond_inner_rtol, &
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
      if (to_lower(trim(adjustl(self%precond))) == "cartesian-aux") then
         call check_cartesian_aux(self, error)
         if (allocated(error)) return
      end if
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

   !> Settings the cartesian-aux preconditioner needs
   !>
   !> - Lattice sizes and inner controls in range
   !> - A solver that applies preconditioners (jfnk or hybrid), with exact Jacobian
   !>   products for the closure derivative and without recycling, which assumes a
   !>   fixed preconditioner
   !>
   !> @param[in]  self   parameter values
   !> @param[out] error  first invalid setting
   subroutine check_cartesian_aux(self, error)
      !> Parameter values
      class(model_moz_3d_parameters_type), intent(in) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error

      if (self%precond_aux_n < 4 .or. mod(self%precond_aux_n, 2) /= 0) then
         call fatal_error(error, "cartesian-aux: jfnk.precond_aux_n must be even and at least 4")
      else if (self%precond_coarse_n /= self%precond_aux_n .and. 2*self%precond_coarse_n /= self%precond_aux_n) then
         call fatal_error(error, "cartesian-aux: jfnk.precond_coarse_n must be precond_aux_n/2, "// &
            & "or precond_aux_n without a coarse lattice")
      else if (self%precond_inner_max_iter < 1) then
         call fatal_error(error, "cartesian-aux: jfnk.precond_inner_max_iter must be at least 1")
      else if (.not. positive(self%precond_inner_rtol) .or. self%precond_inner_rtol >= 1.0_wp) then
         call fatal_error(error, "cartesian-aux: jfnk.precond_inner_rtol must lie in (0, 1)")
      else if (.not. self%solver%runs("jfnk")) then
         call fatal_error(error, "cartesian-aux needs the jfnk or hybrid solver")
      else if (.not. self%solver%jfnk%use_exact_jvp) then
         call fatal_error(error, "cartesian-aux needs jfnk.use_exact_jvp = true")
      else if (self%solver%jfnk%recycle > 0) then
         call fatal_error(error, "cartesian-aux changes with each application; set jfnk.recycle = 0")
      end if
   end subroutine check_cartesian_aux

   !> Finite positive value
   !>
   !> @param[in] value Scalar to check
   pure logical function positive(value)
      !> Scalar to check
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

end module moist_model_moz_3d_parameters
