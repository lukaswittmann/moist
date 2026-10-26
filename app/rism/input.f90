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
   public :: rism_input_type, new_rism_input, rism_checkpoint_parameters_type, rism_checkpoint_files_type

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

   !> Files and save interval of the checkpoints of one solve
   type, extends(moist_model_parameters_type) :: rism_checkpoint_files_type
      !> Saved file; unallocated selects the automatic name
      character(len=:), allocatable :: file
      !> File read as initial guess; unallocated reads `file`, and for an automatic
      !> VV name also the other closures' automatic files
      character(len=:), allocatable :: guess
      !> Iterations between intermediate saves; 0 saves only the final field
      integer :: interval = 5
   contains
      procedure :: init_defaults => files_defaults
      procedure :: register_entries => files_entries
      procedure :: validate => files_validate
   end type rism_checkpoint_files_type

   !> Automatic checkpoints of the bulk VV and the UV solves
   type, extends(moist_model_parameters_type) :: rism_checkpoint_parameters_type
      !> Read initial guesses and save states; false reads and writes nothing
      logical :: enabled = .true.
      !> Bulk VV and UV checkpoint files
      type(rism_checkpoint_files_type) :: vv, uv
   contains
      procedure :: init_defaults => checkpoint_defaults
      procedure :: register_entries => checkpoint_entries
      procedure :: validate => checkpoint_validate
   end type rism_checkpoint_parameters_type

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
      !> Checkpoint files of the VV and UV solves
      type(rism_checkpoint_parameters_type) :: checkpoint
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
      call self%checkpoint%init_defaults()
   end subroutine input_defaults

   !> Register native UV keys at the root and application groups alongside them
   !>
   !> - Order of tables and files: system, grid, solute potential, UV settings, solvent,
   !>   checkpoint
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
      call self%register_group("checkpoint", self%checkpoint)
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
      if (.not. allocated(error)) call self%checkpoint%validate(error)
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

   !> Restore automatic checkpoints with the default interval
   !>
   !> @param[inout] self  checkpoint settings
   subroutine checkpoint_defaults(self)
      !> Checkpoint settings
      class(rism_checkpoint_parameters_type), intent(inout) :: self
      type(rism_checkpoint_parameters_type) :: defaults
      self%enabled = defaults%enabled
      call self%vv%init_defaults()
      call self%uv%init_defaults()
   end subroutine checkpoint_defaults

   !> Register the switch and the VV and UV file groups
   !>
   !> @param[inout] self  checkpoint settings
   subroutine checkpoint_entries(self)
      !> Checkpoint settings
      class(rism_checkpoint_parameters_type), intent(inout), target :: self
      call self%register_logical("enabled", self%enabled, level=2)
      call self%register_group("vv", self%vv)
      call self%register_group("uv", self%uv)
   end subroutine checkpoint_entries

   !> Check both file groups
   !>
   !> @param[inout] self  checkpoint settings
   !> @param[out] error   first invalid setting
   subroutine checkpoint_validate(self, error)
      !> Checkpoint settings
      class(rism_checkpoint_parameters_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      call self%vv%validate(error)
      if (.not. allocated(error)) call self%uv%validate(error)
   end subroutine checkpoint_validate

   !> Restore automatic names and the default interval
   !>
   !> @param[inout] self  checkpoint files
   subroutine files_defaults(self)
      !> Checkpoint files
      class(rism_checkpoint_files_type), intent(inout) :: self
      type(rism_checkpoint_files_type) :: defaults
      if (allocated(self%file)) deallocate (self%file)
      if (allocated(self%guess)) deallocate (self%guess)
      self%interval = defaults%interval
   end subroutine files_defaults

   !> Register the saved file, the guess file and the interval
   !>
   !> @param[inout] self  checkpoint files
   subroutine files_entries(self)
      !> Checkpoint files
      class(rism_checkpoint_files_type), intent(inout), target :: self
      call self%register_alloc_string("file", self%file, level=2)
      call self%register_alloc_string("guess", self%guess, level=2)
      call self%register_int_scalar("interval", self%interval)
   end subroutine files_entries

   !> Require nonblank file names and a nonnegative interval
   !>
   !> @param[inout] self  checkpoint files
   !> @param[out] error   invalid name or interval
   subroutine files_validate(self, error)
      !> Checkpoint files
      class(rism_checkpoint_files_type), intent(inout) :: self
      !> Invalid name or interval
      type(error_type), allocatable, intent(out) :: error
      if (self%interval < 0) then
         call fatal_error(error, "Checkpoint interval must be >= 0 (0 saves only the final field)")
         return
      end if
      if (allocated(self%file)) then
         if (len_trim(self%file) == 0) then
            call fatal_error(error, "Checkpoint file name must not be blank")
            return
         end if
      end if
      if (allocated(self%guess)) then
         if (len_trim(self%guess) == 0) call fatal_error(error, "Checkpoint guess file name must not be blank")
      end if
   end subroutine files_validate

end module moist_app_rism_input
