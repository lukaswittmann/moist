!> Settings of the 1D VV solve of the bulk solvent
!>
!> - JSON/TOML keys: theory, drism.*, closure, grid.*, convergence.*, damping.*, solver.*
!> - solver.* holds the solver controls as a nested group plus
!>   solver.jfnk.precond* and solver.jfnk.deflation_*
!> - Verbosity and thread count come from the run context
module moist_model_moz_solvent_vv_parameters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type, print_hidden
   use moist_model_moz_closure, only: moist_model_moz_closure_type, new_closure
   use moist_model_moz_solver, only: moist_model_moz_solver_parameters_type, &
      & moist_model_moz_solver_picard_parameters_type, moist_model_moz_solver_diis_parameters_type, &
      & moist_model_moz_solver_jfnk_parameters_type, moist_model_moz_solver_jfnk_adjoint_parameters_type, &
      & moist_model_moz_solver_hybrid_parameters_type, check_preconditioner_name
   use moist_model_moz_1d_preconditioner, only: preconditioners_1d_vv
   use moist_model_moz_solvent_theory, only: canonical_bulk_theory
   implicit none(type, external)
   private

   public :: solvent_vv_parameters_type

   !> Settings of the 1D VV solve
   type, extends(moist_model_parameters_type) :: solvent_vv_parameters_type
      !> Bulk theory: ordinary RISM by default; XRISM is an alias
      character(len=8) :: theory = "RISM"
      !> DRISM Gaussian width, bohr
      real(wp) :: drism_smear = 0.6682_wp
      !> Closure: HNC, KH or PSE2-4
      character(len=8) :: closure = "KH"
      !> Points of the uniform radial pair
      integer :: npts = 8192
      !> Radial spacing, bohr
      real(wp) :: dr = 0.015_wp
      !> RMS residual below which a stage is converged
      real(wp) :: tolerance = 1.0e-11_wp
      !> Outer iterations per stage
      integer :: maxiter = 128
      !> Coupling stages lambda = 1/n, ..., 1 of the short-range potential
      integer :: nlambda = 1
      !> Small-k damping of updates the solver allows to post-process
      logical :: damp_small_k = .true.
      !> Wave number and k = 0 weight of the damping, 1/bohr
      real(wp) :: damp_small_k_k0 = 1.0_wp, damp_small_k_min_weight = 0.15_wp
      !> Preconditioner of the Newton-Krylov solves: none, spectral or kspace
      character(len=8) :: precond = "spectral"
      !> Wave number, 1/bohr, and k = 0 weight of the spectral filter
      real(wp) :: precond_k0 = 1.0_wp, precond_w0 = 0.1_wp
      !> Gaussian width of the deflation profiles, bohr
      real(wp) :: deflation_sigma = 2.35_wp
      !> Profiles per site pair of the deflation space
      integer :: deflation_nfun = 2
      !> Solver controls; the VV solve deflates its Newton-Krylov steps by default
      ! Intel requires explicit constructors for the nested parameter defaults
      type(moist_model_moz_solver_parameters_type) :: solver = &
         & moist_model_moz_solver_parameters_type( &
         & picard=moist_model_moz_solver_picard_parameters_type(), &
         & diis=moist_model_moz_solver_diis_parameters_type(), &
         & jfnk=moist_model_moz_solver_jfnk_parameters_type( &
         & adjoint=moist_model_moz_solver_jfnk_adjoint_parameters_type(), deflation=.true.), &
         & hybrid=moist_model_moz_solver_hybrid_parameters_type())
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_vv_defaults
      !> Declare fields for JSON/TOML input, output and printing
      procedure :: register_entries => register_vv_entries
      !> Check names and ranges
      procedure :: validate => validate_vv_parameters
      !> Whether two settings define the same bulk equations and radial grid
      procedure :: same_state => same_vv_state
   end type solvent_vv_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_vv_defaults(self)
      !> Parameter values
      class(solvent_vv_parameters_type), intent(inout) :: self
      type(solvent_vv_parameters_type) :: defaults

      self%theory = defaults%theory
      self%drism_smear = defaults%drism_smear
      self%closure = defaults%closure
      self%npts = defaults%npts
      self%dr = defaults%dr
      self%tolerance = defaults%tolerance
      self%maxiter = defaults%maxiter
      self%nlambda = defaults%nlambda
      self%damp_small_k = defaults%damp_small_k
      self%damp_small_k_k0 = defaults%damp_small_k_k0
      self%damp_small_k_min_weight = defaults%damp_small_k_min_weight
      self%precond = defaults%precond
      self%precond_k0 = defaults%precond_k0
      self%precond_w0 = defaults%precond_w0
      self%deflation_sigma = defaults%deflation_sigma
      self%deflation_nfun = defaults%deflation_nfun
      call self%solver%init_defaults()
      self%solver%jfnk%deflation = defaults%solver%jfnk%deflation
   end subroutine init_vv_defaults

   !> Declare parameter fields for JSON/TOML input, output and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_vv_entries(self)
      !> Parameter values
      class(solvent_vv_parameters_type), intent(inout), target :: self

      !> Print levels: DRISM settings print only for DRISM; damping acts on Picard and DIIS
      !> updates; preconditioner and deflation settings print with JFNK
      integer :: drism_level, damping_level, precond_level, filter_level

      drism_level = merge(2, print_hidden, canonical_bulk_theory(self%theory) == "DRISM")
      damping_level = self%solver%group_level("picard", 3)
      precond_level = self%solver%group_level("jfnk", 2)
      filter_level = self%solver%group_level("jfnk", 3)
      call self%register_string("theory", self%theory, level=2)
      call self%register_real_scalar("drism.smear", self%drism_smear, level=max(3, drism_level))
      call self%register_string("closure", self%closure, level=2)
      call self%register_int_scalar("grid.npts", self%npts, level=2)
      call self%register_real_scalar("grid.dr", self%dr, level=2)
      call self%register_real_scalar("convergence.tolerance", self%tolerance, level=2)
      call self%register_int_scalar("convergence.maxiter", self%maxiter, level=2)
      call self%register_int_scalar("convergence.nlambda", self%nlambda)
      call self%register_logical("damping.small_k", self%damp_small_k, level=damping_level)
      call self%register_real_scalar("damping.small_k_k0", self%damp_small_k_k0, level=damping_level)
      call self%register_real_scalar("damping.small_k_min_weight", self%damp_small_k_min_weight, level=damping_level)
      call self%register_group("solver", self%solver)
      call self%register_string("solver.jfnk.precond", self%precond, level=precond_level)
      call self%register_real_scalar("solver.jfnk.precond_k0", self%precond_k0, level=filter_level)
      call self%register_real_scalar("solver.jfnk.precond_w0", self%precond_w0, level=filter_level)
      call self%register_real_scalar("solver.jfnk.deflation_sigma", self%deflation_sigma, level=filter_level)
      call self%register_int_scalar("solver.jfnk.deflation_nfun", self%deflation_nfun, level=filter_level)
   end subroutine register_vv_entries

   !> Check names and ranges, including the solver controls
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error First invalid setting
   subroutine validate_vv_parameters(self, error)
      !> Parameter values
      class(solvent_vv_parameters_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      class(moist_model_moz_closure_type), allocatable :: closure

      if (canonical_bulk_theory(self%theory) == "") then
         call fatal_error(error, "Unknown bulk theory: "//trim(self%theory)//" (choose RISM or DRISM)")
         return
      end if
      if (canonical_bulk_theory(self%theory) == "DRISM") then
         if (.not. ieee_is_finite(self%drism_smear) .or. self%drism_smear < 0.0_wp) then
            call fatal_error(error, "DRISM smear must be finite and >= 0 (bohr)")
            return
         end if
      end if
      call new_closure(closure, self%closure, error)
      if (allocated(error)) return
      call check_preconditioner_name(self%precond, preconditioners_1d_vv, "1D VV", error)
      if (allocated(error)) return
      if (self%npts < 2 .or. .not. positive(self%dr)) then
         call fatal_error(error, "1D VV grid needs at least two points and a finite positive spacing")
      else if (.not. positive(self%tolerance) .or. self%maxiter < 1 .or. self%nlambda < 1) then
         call fatal_error(error, "1D VV convergence needs a positive tolerance, iteration limit and stage count")
      else if (.not. positive(self%damp_small_k_k0) .or. .not. weight(self%damp_small_k_min_weight)) then
         call fatal_error(error, "1D VV damping needs a positive k0 and a k = 0 weight in [0, 1]")
      else if (.not. positive(self%precond_k0) .or. .not. weight(self%precond_w0)) then
         call fatal_error(error, "1D VV spectral filter needs a positive k0 and a k = 0 weight in [0, 1]")
      else if (.not. positive(self%deflation_sigma) .or. self%deflation_nfun < 1) then
         call fatal_error(error, "1D VV deflation needs a positive width and at least one profile")
      else
         call self%solver%validate(error)
      end if
   end subroutine validate_vv_parameters

   !> Compare the physical bulk equations and the discretization
   !>
   !> - Solver controls do not change the state
   !> - Achieved tolerance is checked by the solvent
   !>
   !> @param[in] self Accepted bulk settings
   !> @param[in] other Requested bulk settings
   pure logical function same_vv_state(self, other) result(same)
      !> Accepted bulk settings
      class(solvent_vv_parameters_type), intent(in) :: self
      !> Requested bulk settings
      type(solvent_vv_parameters_type), intent(in) :: other
      same = canonical_bulk_theory(self%theory) == canonical_bulk_theory(other%theory) .and. &
         & to_lower(trim(self%closure)) == to_lower(trim(other%closure)) .and. &
         & self%npts == other%npts .and. self%dr == other%dr
      if (canonical_bulk_theory(self%theory) == "DRISM") then
         same = same .and. self%drism_smear == other%drism_smear
      end if
   end function same_vv_state

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

end module moist_model_moz_solvent_vv_parameters
