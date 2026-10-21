!> Application settings layered over the native UV and bulk-VV parameters
module moist_app_rism_input
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_math_grid_3d_factory, only: rism_grid_parameters_type
   use moist_model_moz_potential_factory, only: rism_potential_parameters_type
   implicit none(type, external)
   private
   public :: rism_input_type, new_rism_input

   !> State point and optional total molecular charge
   type, extends(moist_model_parameters_type) :: rism_system_parameters_type
      !> Solvent name or id, resolved against the solvent table; empty selects the default
      character(len=64) :: solvent = ""
      !> Temperature in K and pressure in Pa; unallocated selects the solvent-system default
      real(wp), allocatable :: temperature, pressure
      !> Total molecular charge in e; unallocated means no override
      real(wp), allocatable :: charge
   contains
      procedure :: init_defaults => system_defaults
      procedure :: register_entries => system_entries
      procedure :: validate => system_validate
   end type rism_system_parameters_type

   !> One owned input object; UV implementation selected by the CLI command
   type, extends(moist_model_parameters_type) :: rism_input_type
      !> Native UV parameters, allocated by new_rism_input (1D or 3D)
      class(moist_model_parameters_type), allocatable :: model_parameters
      !> State point and charge
      type(rism_system_parameters_type) :: system
      !> Radial or spatial grid settings
      type(rism_grid_parameters_type) :: grid
      !> Solute-side and solvent-side potential settings
      type(rism_potential_parameters_type) :: potential, solvent_potential
   contains
      procedure :: init_defaults => input_defaults
      procedure :: register_entries => input_entries
      procedure :: validate => input_validate
   end type rism_input_type

contains

   !> Select the native UV settings and the matching application defaults
   !>
   !> @param[out] self   input object
   !> @param[in] mode    CLI model command: rism1d or rism3d
   !> @param[out] error  unsupported command
   subroutine new_rism_input(self, mode, error)
      !> Input object
      type(rism_input_type), intent(out) :: self
      !> CLI model command: rism1d or rism3d
      character(len=*), intent(in) :: mode
      !> Unsupported command
      type(error_type), allocatable, intent(out) :: error

      select case (trim(mode))
      case ("rism1d")
         allocate(model_moz_1d_parameters_type :: self%model_parameters)
         call self%grid%set_dimension(1)
      case ("rism3d")
         allocate(model_moz_3d_parameters_type :: self%model_parameters)
      case default
         call fatal_error(error, "Unknown RISM model command: "//trim(mode))
         return
      end select
      call self%solvent_potential%set_solvent_side(.true.)
      call self%init_defaults()
   end subroutine new_rism_input

   !> Reset values while preserving the selected UV implementation
   !>
   !> @param[inout] self  input object
   subroutine input_defaults(self)
      !> Input object
      class(rism_input_type), intent(inout) :: self
      if (allocated(self%model_parameters)) call self%model_parameters%init_defaults()
      call self%system%init_defaults()
      call self%grid%init_defaults()
      call self%potential%init_defaults()
      call self%solvent_potential%init_defaults()
   end subroutine input_defaults

   !> Register native UV keys at the root and application groups alongside them
   !>
   !> - Order of tables and files: system, grid, solute potential, UV settings, solvent
   !>
   !> @param[inout] self  input object
   subroutine input_entries(self)
      !> Input object
      class(rism_input_type), intent(inout), target :: self
      call self%register_group("system", self%system)
      call self%register_group("grid", self%grid)
      call self%register_group("potential", self%potential)
      if (allocated(self%model_parameters)) call self%register_group("", self%model_parameters)
      call self%register_group("solvent.potential", self%solvent_potential)
   end subroutine input_entries

   !> Validate independent settings; structure-dependent checks belong to setup
   !>
   !> @param[inout] self  input object
   !> @param[out] error   first invalid setting
   subroutine input_validate(self, error)
      !> Input object
      class(rism_input_type), intent(inout) :: self
      !> First invalid setting
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

   !> Restore the default state point, with no explicit overrides
   !>
   !> @param[inout] self  state point
   subroutine system_defaults(self)
      !> State point
      class(rism_system_parameters_type), intent(inout) :: self
      self%solvent = ""
      if (allocated(self%temperature)) deallocate(self%temperature)
      if (allocated(self%pressure)) deallocate(self%pressure)
      if (allocated(self%charge)) deallocate(self%charge)
   end subroutine system_defaults

   !> Register state-point fields
   !>
   !> @param[inout] self  state point
   subroutine system_entries(self)
      !> State point
      class(rism_system_parameters_type), intent(inout), target :: self
      call self%register_string("solvent", self%solvent, level=2)
      call self%register_alloc_real_scalar("temperature", self%temperature, level=2)
      call self%register_alloc_real_scalar("pressure", self%pressure, level=2)
      call self%register_alloc_real_scalar("charge", self%charge, level=2)
   end subroutine system_entries

   !> Check the charge; the solvent system checks the state point after CLI overrides
   !>
   !> @param[inout] self  state point
   !> @param[out] error   invalid charge
   subroutine system_validate(self, error)
      !> State point
      class(rism_system_parameters_type), intent(inout) :: self
      !> Invalid charge
      type(error_type), allocatable, intent(out) :: error
      if (allocated(self%charge)) then
         if (.not. ieee_is_finite(self%charge)) call fatal_error(error, "System charge must be finite")
      end if
   end subroutine system_validate

end module moist_app_rism_input
