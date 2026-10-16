!> Application settings layered over the native UV and bulk-VV parameters
module moist_input
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   implicit none(type, external)
   private
   public :: rism_input_type, new_rism_input, rism_grid_parameters_type, rism_potential_parameters_type

   !> State point and optional total molecular charge
   type, extends(moist_model_parameters_type) :: rism_system_parameters_type
      character(len=64) :: solvent = ""
      real(wp) :: temperature = 298.15_wp, pressure = 101325.0_wp
      real(wp), allocatable :: charge
   contains
      procedure :: init_defaults => system_defaults
      procedure :: register_entries => system_entries
      procedure :: validate => system_validate
   end type rism_system_parameters_type

   !> UV grid configuration; only keys used by the selected grid are registered
   type, extends(moist_model_parameters_type) :: rism_grid_parameters_type
      integer, private :: dimension = 3
      character(len=16) :: kind = "cartesian"
      real(wp) :: dr = 0.6_wp
      integer :: npts = 8192, nx = 128, ny = 128, nz = 128
      real(wp) :: margin = 10.0_wp
      integer :: nrad = 50, angular_degree = 17, becke_k = 3
      real(wp) :: rmax = 10.0_wp, kbuffer = 2.0_wp, nufft_tolerance = 1.0e-10_wp
      !> Dimensionless real exponent of the bounded HandyMod radial mapping
      real(wp) :: handymod_exponent = 2.0_wp
   contains
      procedure :: init_defaults => grid_defaults
      procedure :: register_entries => grid_entries
      procedure :: validate => grid_validate
   end type rism_grid_parameters_type

   !> Choices for assembling native, host-independent potential terms
   type, extends(moist_model_parameters_type) :: rism_potential_parameters_type
      logical, private :: solvent_side = .false.
      character(len=16) :: lj = "uff", lj_fallback = "", charges = "fixed", environment = "gas"
      character(len=24) :: mixing = "lorentz-berthelot"
      logical :: ng_split = .true.
      real(wp) :: alpha = 1.0_wp
      !> Optional explicit per-atom values: electrons, bohr and Hartree
      real(wp), allocatable :: q(:), sigma(:), epsilon(:)
   contains
      procedure :: init_defaults => potential_defaults
      procedure :: register_entries => potential_entries
      procedure :: validate => potential_validate
   end type rism_potential_parameters_type

   !> One owned input object; UV implementation selected by the CLI command
   type, extends(moist_model_parameters_type) :: rism_input_type
      class(moist_model_parameters_type), allocatable :: model_parameters
      type(rism_system_parameters_type) :: system
      type(rism_grid_parameters_type) :: grid
      type(rism_potential_parameters_type) :: potential, solvent_potential
   contains
      procedure :: init_defaults => input_defaults
      procedure :: register_entries => input_entries
      procedure :: validate => input_validate
   end type rism_input_type

contains

   !> Select the native UV settings and the matching application defaults
   !>
   !> @param[out] self Input object
   !> @param[in] mode CLI model command: rism1d or rism3d
   !> @param[out] error Unsupported command
   subroutine new_rism_input(self, mode, error)
      type(rism_input_type), intent(out) :: self
      character(len=*), intent(in) :: mode
      type(error_type), allocatable, intent(out) :: error

      select case (trim(mode))
      case ("rism1d")
         allocate(model_moz_1d_parameters_type :: self%model_parameters)
         self%grid%dimension = 1
      case ("rism3d")
         allocate(model_moz_3d_parameters_type :: self%model_parameters)
      case default
         call fatal_error(error, "Unknown RISM model command: "//trim(mode))
         return
      end select
      self%solvent_potential%solvent_side = .true.
      call self%init_defaults()
   end subroutine new_rism_input

   !> Reset values while preserving the selected UV implementation
   !> @param[inout] self Input object
   subroutine input_defaults(self)
      class(rism_input_type), intent(inout) :: self
      if (allocated(self%model_parameters)) call self%model_parameters%init_defaults()
      call self%system%init_defaults()
      call self%grid%init_defaults()
      call self%potential%init_defaults()
      call self%solvent_potential%init_defaults()
   end subroutine input_defaults

   !> Register native UV keys at the root and application groups alongside them
   !> @param[inout] self Input object
   subroutine input_entries(self)
      class(rism_input_type), intent(inout), target :: self
      if (allocated(self%model_parameters)) call self%register_group("", self%model_parameters)
      call self%register_group("system", self%system)
      call self%register_group("grid", self%grid)
      call self%register_group("potential", self%potential)
      call self%register_group("solvent.potential", self%solvent_potential)
   end subroutine input_entries

   !> Validate independent settings; structure-dependent checks belong to setup
   !> @param[inout] self Input object
   !> @param[out] error First invalid setting
   subroutine input_validate(self, error)
      class(rism_input_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%model_parameters)) then
         call fatal_error(error, "Select a RISM model before reading its input")
         return
      end if
      call self%model_parameters%validate(error)
      if (.not. allocated(error)) call self%system%validate(error)
      if (.not. allocated(error)) call self%grid%validate(error)
      if (.not. allocated(error)) call self%potential%validate(error)
      if (.not. allocated(error)) call self%solvent_potential%validate(error)
      if (allocated(error)) return
      if ((self%potential%lj == "oplsaa" .or. self%potential%lj_fallback == "oplsaa" &
         & .or. self%solvent_potential%lj == "oplsaa" .or. self%solvent_potential%lj_fallback == "oplsaa") &
         & .and. self%solvent_potential%mixing /= "geometric") then
         call fatal_error(error, "OPLS-AA needs solvent.potential.mixing = 'geometric'")
      end if
   end subroutine input_validate

   !> Restore the default state point, with no explicit charge override
   !> @param[inout] self State point
   subroutine system_defaults(self)
      class(rism_system_parameters_type), intent(inout) :: self
      self%solvent = ""
      self%temperature = 298.15_wp
      self%pressure = 101325.0_wp
      if (allocated(self%charge)) deallocate(self%charge)
   end subroutine system_defaults

   !> Register state-point fields
   !> @param[inout] self State point
   subroutine system_entries(self)
      class(rism_system_parameters_type), intent(inout), target :: self
      call self%register_string("solvent", self%solvent)
      call self%register_real_scalar("temperature", self%temperature)
      call self%register_real_scalar("pressure", self%pressure)
      call self%register_alloc_real_scalar("charge", self%charge)
   end subroutine system_entries

   !> Check state-point ranges; solvent identity is resolved after CLI overrides
   !> @param[inout] self State point
   !> @param[out] error Invalid temperature, pressure or charge
   subroutine system_validate(self, error)
      class(rism_system_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      if (.not. positive(self%temperature) .or. .not. nonnegative(self%pressure)) then
         call fatal_error(error, "System temperature must be positive and pressure nonnegative (K, Pa)")
      else if (allocated(self%charge)) then
         if (.not. ieee_is_finite(self%charge)) call fatal_error(error, "System charge must be finite")
      end if
   end subroutine system_validate

   !> Restore UV grid defaults for the selected dimension
   !> @param[inout] self Grid settings
   subroutine grid_defaults(self)
      class(rism_grid_parameters_type), intent(inout) :: self
      type(rism_grid_parameters_type) :: defaults
      self%kind = "cartesian"
      self%dr = defaults%dr
      if (self%dimension == 1) then
         self%kind = "radial"
         self%dr = 0.015_wp
      end if
      self%npts = defaults%npts
      self%nx = defaults%nx
      self%ny = defaults%ny
      self%nz = defaults%nz
      self%margin = defaults%margin
      self%nrad = defaults%nrad
      self%angular_degree = defaults%angular_degree
      self%becke_k = defaults%becke_k
      self%rmax = defaults%rmax
      self%handymod_exponent = defaults%handymod_exponent
      self%kbuffer = defaults%kbuffer
      self%nufft_tolerance = defaults%nufft_tolerance
   end subroutine grid_defaults

   !> Register the selected grid's keys; unused grid keys fail strict input
   !> @param[inout] self Grid settings
   subroutine grid_entries(self)
      class(rism_grid_parameters_type), intent(inout), target :: self
      call self%register_string("kind", self%kind)
      call self%register_real_scalar("dr", self%dr)
      select case (to_lower(trim(adjustl(self%kind))))
      case ("radial")
         call self%register_int_scalar("npts", self%npts)
      case ("cartesian")
         call self%register_int_scalar("nx", self%nx)
         call self%register_int_scalar("ny", self%ny)
         call self%register_int_scalar("nz", self%nz)
         call self%register_real_scalar("margin", self%margin)
      case ("molecular")
         call self%register_int_scalar("nrad", self%nrad)
         call self%register_int_scalar("angular_degree", self%angular_degree)
         call self%register_real_scalar("rmax", self%rmax)
         call self%register_real_scalar("handymod_exponent", self%handymod_exponent)
         call self%register_int_scalar("becke_k", self%becke_k)
         call self%register_real_scalar("kbuffer", self%kbuffer)
         call self%register_real_scalar("nufft_tolerance", self%nufft_tolerance)
      case default
         ! Unknown grid kinds are rejected by grid_validate.
         continue
      end select
   end subroutine grid_entries

   !> Check grid names and dimensions without realizing geometry
   !> @param[inout] self Grid settings
   !> @param[out] error Unsupported grid or invalid size
   subroutine grid_validate(self, error)
      class(rism_grid_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      self%kind = to_lower(trim(adjustl(self%kind)))
      if (.not. positive(self%dr)) then
         call fatal_error(error, "UV grid.dr must be positive (bohr)")
         return
      end if
      if (self%dimension == 1 .and. self%kind /= "radial") then
         call fatal_error(error, "rism1d needs grid.kind = 'radial'")
         return
      end if
      select case (self%kind)
      case ("radial")
         if (self%dimension /= 1 .or. self%npts < 2) then
            call fatal_error(error, "Radial UV grids need rism1d and at least two points")
         end if
      case ("cartesian")
         if (min(self%nx, self%ny, self%nz) < 2 .or. .not. nonnegative(self%margin)) then
            call fatal_error(error, "Cartesian UV grids need at least two points per axis and a nonnegative margin")
         end if
      case ("molecular")
         if (self%nrad < 1 .or. self%angular_degree < 0 .or. self%becke_k < 1 .or. &
            & .not. positive(self%rmax) .or. .not. positive(self%handymod_exponent) &
            & .or. .not. nonnegative(self%kbuffer) &
            & .or. .not. positive(self%nufft_tolerance)) then
            call fatal_error(error, "Invalid molecular UV grid settings")
         end if
      case default
         call fatal_error(error, "Unknown UV grid.kind: "//trim(self%kind))
      end select
   end subroutine grid_validate

   !> Restore potential choices; explicit solute charges have no guessed default
   !> @param[inout] self Potential settings
   subroutine potential_defaults(self)
      class(rism_potential_parameters_type), intent(inout) :: self
      self%lj = "uff"
      self%charges = "fixed"
      if (self%solvent_side) then
         self%lj = "cspce"
         self%charges = "cspce"
      end if
      self%lj_fallback = ""
      self%environment = "gas"
      self%mixing = "lorentz-berthelot"
      self%ng_split = .true.
      self%alpha = 1.0_wp
      if (allocated(self%q)) deallocate(self%q)
      if (allocated(self%sigma)) deallocate(self%sigma)
      if (allocated(self%epsilon)) deallocate(self%epsilon)
   end subroutine potential_defaults

   !> Register potential selection and optional per-atom arrays
   !> @param[inout] self Potential settings
   subroutine potential_entries(self)
      class(rism_potential_parameters_type), intent(inout), target :: self
      call self%register_string("lj", self%lj)
      call self%register_string("lj_fallback", self%lj_fallback)
      call self%register_string("charges", self%charges)
      call self%register_alloc_real_vector("q", self%q)
      call self%register_alloc_real_vector("sigma", self%sigma)
      call self%register_alloc_real_vector("epsilon", self%epsilon)
      call self%register_logical("ng_split", self%ng_split)
      call self%register_real_scalar("alpha", self%alpha)
      if (self%solvent_side) then
         call self%register_string("mixing", self%mixing)
         select case (to_lower(trim(adjustl(self%charges))))
         case ("hirshfeld", "resp", "mbis", "chelpg")
            call self%register_string("environment", self%environment)
         case default
            ! Other charge models do not use an environment setting.
            continue
         end select
      end if
   end subroutine potential_entries

   !> Check choices and array values before native potential construction
   !> @param[inout] self Potential settings
   !> @param[out] error Unknown choice, missing data or an unused array
   subroutine potential_validate(self, error)
      class(rism_potential_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      self%lj = to_lower(trim(adjustl(self%lj)))
      self%lj_fallback = to_lower(trim(adjustl(self%lj_fallback)))
      self%charges = to_lower(trim(adjustl(self%charges)))
      self%environment = to_lower(trim(adjustl(self%environment)))
      self%mixing = to_lower(trim(adjustl(self%mixing)))
      if (.not. valid_lj(self%lj) .or. &
         & (len_trim(self%lj_fallback) > 0 .and. .not. valid_lj(self%lj_fallback)) &
         & .or. self%lj_fallback == "custom") then
         call fatal_error(error, "LJ choices are uff, dreiding, tm, gaff, oplsaa, spce, cspce or custom (primary only)")
         return
      end if
      if (self%lj == "custom") then
         if (.not. allocated(self%sigma) .or. .not. allocated(self%epsilon)) then
            call fatal_error(error, "Custom LJ needs per-atom sigma (bohr) and epsilon (Hartree)")
            return
         end if
      else if (allocated(self%sigma) .or. allocated(self%epsilon)) then
         call fatal_error(error, "Explicit sigma/epsilon require lj = 'custom'")
         return
      end if
      select case (self%charges)
      case ("fixed")
         if (.not. allocated(self%q)) then
            call fatal_error(error, "Fixed charges need potential.q: one charge per atom in electrons")
            return
         end if
      case ("spce", "cspce", "none")
         continue
      case ("hirshfeld", "resp", "mbis", "chelpg")
         if (.not. self%solvent_side) then
            call fatal_error(error, "Table charge models belong to solvent.potential; supply explicit solute q")
            return
         end if
         select case (self%environment)
         case ("gas", "solvent", "conductor")
            continue
         case default
            call fatal_error(error, "Charge environment must be gas, solvent or conductor")
            return
         end select
      case default
         call fatal_error(error, "Unknown potential.charges: "//trim(self%charges))
         return
      end select
      if (self%charges /= "fixed" .and. allocated(self%q)) then
         call fatal_error(error, "Explicit q requires charges = 'fixed'")
      else if (.not. positive(self%alpha)) then
         call fatal_error(error, "Potential alpha must be positive (1/bohr)")
      else if (self%mixing /= "lorentz-berthelot" .and. self%mixing /= "geometric") then
         call fatal_error(error, "Solvent mixing must be lorentz-berthelot or geometric")
      end if
   end subroutine potential_validate

   !> Whether an LJ factory is available
   !> @param[in] name Normalized choice
   pure logical function valid_lj(name)
      character(len=*), intent(in) :: name
      select case (name)
      case ("uff", "dreiding", "tm", "gaff", "oplsaa", "spce", "cspce", "custom")
         valid_lj = .true.
      case default
         valid_lj = .false.
      end select
   end function valid_lj

   !> Finite positive scalar
   !> @param[in] value Scalar to check
   pure logical function positive(value)
      real(wp), intent(in) :: value
      positive = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive

   !> Finite nonnegative scalar
   !> @param[in] value Scalar to check
   pure logical function nonnegative(value)
      real(wp), intent(in) :: value
      nonnegative = ieee_is_finite(value) .and. value >= 0.0_wp
   end function nonnegative
end module moist_input
