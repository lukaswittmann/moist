!> Settings of a whole RISM input file and the shared steps of building a UV model from them
!>
!> - One file: UV keys at the root with the nested `[solvent]` VV settings, `[grid]`,
!>   `[potential]`, `[solvent.potential]` and `[checkpoint]`
!> - The solvent and its state point come from `[solvent]` name, temperature and pressure,
!>   or from the host as a solvent system, which takes precedence; the total charge comes
!>   from the structure at update
!> - 1D and 3D extensions add the UV settings and overload `new_moz_1d_model` / `new_moz_3d_model`
!> - Restore the defaults with `init_defaults` or `read_file` before use
!>
!> Usage:
!>
!>    type(model_moz_3d_input_type) :: input
!>    type(model_moz_3d_type) :: model
!>    call input%read_file("rism.toml", error, strict=.true.)
!>    if (.not. allocated(error)) call new_moz_3d_model(model, ctx, input, error)
!>    if (.not. allocated(error)) call model%update(mol, error)
module moist_model_moz_input
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_context, only: moist_context_type
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_3d_factory, only: rism_grid_parameters_type
   use moist_model_moz_potential_factory, only: rism_potential_parameters_type, configure_potential
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent, solvent_vv_parameters_type
   use moist_model_moz_solvent_theory, only: canonical_bulk_theory
   use moist_model_moz_solvent_water, only: new_spce_water_geometry, is_water
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type, new_file_checkpoint
   use moist_model_moz_closure_ramp, only: closure_ladder
   use moist_data_solvents, only: get_solvent_id
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   implicit none(type, external)
   private
   public :: model_moz_input_type, moz_solvent_state_parameters_type
   public :: moz_checkpoint_parameters_type, moz_checkpoint_files_type
   public :: moz_input_defaults, moz_input_entries, moz_input_validate
   public :: new_input_solvent, new_vv_checkpoint, new_uv_checkpoint, vv_checkpoint_name

   !> Closures whose automatic VV files can warm-start another closure, in order of preference
   character(len=4), parameter :: fallback_closures(size(closure_ladder)) = &
      & closure_ladder(size(closure_ladder):1:-1)

   !> Error of settings used before their defaults were restored
   character(len=*), parameter :: undefaulted = &
      & "Restore the RISM input defaults with init_defaults or read_file before use"

   !> Solvent name and state point, read from the `[solvent]` table beside the VV settings
   type, extends(moist_model_parameters_type) :: moz_solvent_state_parameters_type
      !> Solvent name or alias of the solvent table; unallocated leaves the solvent to the host
      character(len=:), allocatable :: name
      !> Temperature in K and pressure in Pa; unallocated selects the solvent's defaults
      real(wp), allocatable :: temperature, pressure
   contains
      procedure :: init_defaults => state_defaults
      procedure :: register_entries => state_entries
      procedure :: validate => state_validate
   end type moz_solvent_state_parameters_type

   !> Files, directory and save interval of the checkpoints of one solve
   type, extends(moist_model_parameters_type) :: moz_checkpoint_files_type
      !> Saved file; unallocated selects the automatic name
      character(len=:), allocatable :: file
      !> File read as initial guess; unallocated reads `file`, and for an automatic
      !> VV name also the other closures' automatic files
      character(len=:), allocatable :: guess
      !> Directory of the automatic names; unallocated selects the working directory
      character(len=:), allocatable :: directory
      !> Iterations between intermediate saves; 0 saves only the final field
      integer :: interval = 0
   contains
      procedure :: init_defaults => files_defaults
      procedure :: register_entries => files_entries
      procedure :: validate => files_validate
   end type moz_checkpoint_files_type

   !> Automatic checkpoints of the bulk VV and the UV solves
   type, extends(moist_model_parameters_type) :: moz_checkpoint_parameters_type
      !> Read initial guesses and save states; false reads and writes nothing
      logical :: enabled = .true.
      !> Bulk VV and UV checkpoint files
      type(moz_checkpoint_files_type) :: vv, uv
   contains
      procedure :: init_defaults => checkpoint_defaults
      procedure :: register_entries => checkpoint_entries
      procedure :: validate => checkpoint_validate
   end type moz_checkpoint_parameters_type

   !> Settings shared by the 1D and 3D RISM input files
   !>
   !> - Extensions bind moz_input_defaults, moz_input_entries and moz_input_validate
   !>   into their own init_defaults, register_entries and validate
   type, abstract, extends(moist_model_parameters_type) :: model_moz_input_type
      !> Solvent name and state point
      type(moz_solvent_state_parameters_type) :: solvent_state
      !> Radial or spatial grid settings
      type(rism_grid_parameters_type) :: grid
      !> Solute-side and solvent-side potential settings
      type(rism_potential_parameters_type) :: potential, solvent_potential
      !> Checkpoint files of the VV and UV solves
      type(moz_checkpoint_parameters_type) :: checkpoint
      !> Whether the defaults were restored; settings are unusable before
      logical, private :: defaulted = .false.
   contains
      !> Solvent system of the `[solvent]` name and state point
      procedure :: new_system => input_new_system
      !> Construct the UV model of the input with a host solvent system, not yet updated
      procedure(new_model_ifc), deferred :: new_model_system
      !> Construct the UV model of the input with the `[solvent]` solvent, not yet updated
      procedure :: new_model_file => input_new_model_file
      generic :: new_model => new_model_system, new_model_file
   end type model_moz_input_type

   abstract interface
      !> Construct the UV model of an input, not yet updated
      !>
      !> @param[in]  self    input settings
      !> @param[in]  ctx     borrowed context, outlives the model
      !> @param[in]  system  solvent at its state point
      !> @param[out] model   constructed model; update it with the solute next
      !> @param[out] error   invalid settings, solvent or construction failure
      subroutine new_model_ifc(self, ctx, system, model, error)
         import :: model_moz_input_type, moist_context_type, solvent_system_type, solvation_model_type, error_type
         implicit none(type, external)
         !> Input settings
         class(model_moz_input_type), intent(in) :: self
         !> Borrowed context, outlives the model
         type(moist_context_type), intent(in), target :: ctx
         !> Solvent at its state point
         type(solvent_system_type), intent(in) :: system
         !> Constructed model
         class(solvation_model_type), allocatable, intent(out) :: model
         !> Invalid settings, solvent or construction failure
         type(error_type), allocatable, intent(out) :: error
      end subroutine new_model_ifc
   end interface

contains

   !> Restore the shared defaults; the solvent potential keeps its solvent side
   !>
   !> - Extensions select the grid dimension before calling this
   !>
   !> @param[inout] self  input settings
   subroutine moz_input_defaults(self)
      !> Input settings
      class(model_moz_input_type), intent(inout) :: self
      call self%solvent_state%init_defaults()
      call self%grid%init_defaults()
      call self%potential%init_defaults()
      call self%solvent_potential%set_solvent_side(.true.)
      call self%solvent_potential%init_defaults()
      call self%checkpoint%init_defaults()
      self%defaulted = .true.
   end subroutine moz_input_defaults

   !> Register the shared groups; extensions register their UV keys at the root first
   !>
   !> @param[inout] self  input settings
   subroutine moz_input_entries(self)
      !> Input settings
      class(model_moz_input_type), intent(inout), target :: self
      call self%register_group("solvent", self%solvent_state)
      call self%register_group("grid", self%grid)
      call self%register_group("potential", self%potential)
      call self%register_group("solvent.potential", self%solvent_potential)
      call self%register_group("checkpoint", self%checkpoint)
   end subroutine moz_input_entries

   !> Validate the shared settings; structure-dependent checks belong to the model update
   !>
   !> @param[inout] self   input settings
   !> @param[out]   error  missing defaults or first invalid setting
   subroutine moz_input_validate(self, error)
      !> Input settings
      class(model_moz_input_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%defaulted) then
         call fatal_error(error, undefaulted)
         return
      end if
      call self%solvent_state%validate(error)
      if (.not. allocated(error)) call self%grid%validate(error)
      if (.not. allocated(error)) call self%potential%validate(error)
      if (.not. allocated(error)) call self%solvent_potential%validate(error)
      if (.not. allocated(error)) call self%checkpoint%validate(error)
      if (allocated(error)) return
      if ((any(self%potential%lj == "oplsaa") .or. any(self%solvent_potential%lj == "oplsaa")) &
         & .and. self%solvent_potential%mixing /= "geometric") then
         call fatal_error(error, "OPLS-AA needs solvent.potential.mixing = 'geometric'")
      end if
   end subroutine moz_input_validate

   !> Solvent system of the `[solvent]` name and state point
   !>
   !> - Unset temperature and pressure select the solvent's defaults
   !>
   !> @param[in]  self    input settings with restored defaults
   !> @param[out] system  solvent at its state point
   !> @param[out] error   missing defaults or name, invalid state point or unknown solvent
   subroutine input_new_system(self, system, error)
      !> Input settings with restored defaults
      class(model_moz_input_type), intent(in) :: self
      !> Solvent at its state point
      type(solvent_system_type), intent(out) :: system
      !> Missing or invalid solvent
      type(error_type), allocatable, intent(out) :: error
      type(moz_solvent_state_parameters_type) :: state
      integer :: solvent_id

      if (.not. self%defaulted) then
         call fatal_error(error, undefaulted)
         return
      end if
      state = self%solvent_state
      call state%validate(error)
      if (allocated(error)) return
      if (.not. allocated(state%name)) then
         call fatal_error(error, "No solvent: set solvent.name or pass a solvent system")
         return
      end if
      call get_solvent_id(state%name, solvent_id, error)
      if (.not. allocated(error)) call new_solvent_system(system, solvent_id, temperature=state%temperature, &
         & pressure_si=state%pressure, error=error)
   end subroutine input_new_system

   !> Construct the UV model of the input with the `[solvent]` solvent, not yet updated
   !>
   !> @param[in]  self   input settings
   !> @param[in]  ctx    borrowed context, outlives the model
   !> @param[out] model  constructed model; update it with the solute next
   !> @param[out] error  missing solvent, invalid settings or construction failure
   subroutine input_new_model_file(self, ctx, model, error)
      !> Input settings
      class(model_moz_input_type), intent(in) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Constructed model
      class(solvation_model_type), allocatable, intent(out) :: model
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(solvent_system_type) :: system

      call self%new_system(system, error)
      if (.not. allocated(error)) call self%new_model_system(ctx, system, model, error)
   end subroutine input_new_model_file

   !> Build the updated bulk solvent of validated settings, with its VV checkpoint attached
   !>
   !> - A first LJ table of SPC/E or cSPC/E selects the SPC/E geometry for water
   !>
   !> @param[in]  input    validated input settings
   !> @param[in]  ctx      borrowed context
   !> @param[in]  system   solvent at its state point
   !> @param[in]  params   settings of the VV solve
   !> @param[out] solvent  updated solvent
   !> @param[out] error    solvent construction or potential failure
   subroutine new_input_solvent(input, ctx, system, params, solvent, error)
      !> Validated input settings
      class(model_moz_input_type), intent(in) :: input
      !> Borrowed context
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent at its state point
      type(solvent_system_type), intent(in) :: system
      !> Settings of the VV solve
      type(solvent_vv_parameters_type), intent(in) :: params
      !> Updated solvent
      type(solvent_vv_type), intent(out) :: solvent
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(solvent_system_type) :: state
      type(moz_file_checkpoint_type) :: store

      state = system
      if (input%solvent_potential%lj(1) == "spce" .or. input%solvent_potential%lj(1) == "cspce") then
         if (is_water(state%solv_mol, state%solvent%id)) call new_spce_water_geometry(state%solv_mol)
      end if
      call new_vv_solvent(solvent, ctx, state, error)
      if (.not. allocated(error)) call configure_potential(solvent%potential, input%solvent_potential, error)
      if (.not. allocated(error)) call solvent%update(error)
      if (allocated(error)) return
      if (input%checkpoint%enabled) then
         call new_vv_checkpoint(input%checkpoint%vv, state%solvent%id, params, store)
         call solvent%set_checkpoint(store)
      end if
   end subroutine new_input_solvent

   !> VV file store
   !>
   !> - Saved file: `file`, else the automatic name (vv_checkpoint_name) in `directory`
   !> - Guess files: `guess`, else `file`; with neither, the automatic files of the
   !>   requested closure and then of the other closures, which can only warm-start
   !>
   !> @param[in]  files       VV checkpoint settings
   !> @param[in]  solvent_id  solvent id; 0 for a custom solvent
   !> @param[in]  params      settings of the VV solve
   !> @param[out] store       file store
   subroutine new_vv_checkpoint(files, solvent_id, params, store)
      !> VV checkpoint settings
      type(moz_checkpoint_files_type), intent(in) :: files
      !> Solvent id
      integer, intent(in) :: solvent_id
      !> Settings of the VV solve
      type(solvent_vv_parameters_type), intent(in) :: params
      !> File store
      type(moz_file_checkpoint_type), intent(out) :: store
      character(len=:), allocatable :: file, guesses(:)
      character(len=:), allocatable :: theory
      integer :: i, n, width

      theory = canonical_bulk_theory(params%theory)
      if (allocated(files%file)) then
         file = trim(files%file)
      else
         file = automatic_file(files, vv_checkpoint_name(solvent_id, theory, params%closure))
      end if
      if (allocated(files%guess)) then
         guesses = [character(len=len_trim(files%guess)) :: trim(files%guess)]
      else if (allocated(files%file)) then
         guesses = [character(len=len(file)) :: file]
      else
         width = len(automatic_file(files, vv_checkpoint_name(solvent_id, theory, "PSE4"))) + len_trim(params%closure)
         allocate (character(len=width) :: guesses(1 + size(fallback_closures)))
         guesses(1) = file
         n = 1
         do i = 1, size(fallback_closures)
            if (to_lower(trim(fallback_closures(i))) == to_lower(trim(adjustl(params%closure)))) cycle
            n = n + 1
            guesses(n) = automatic_file(files, vv_checkpoint_name(solvent_id, theory, fallback_closures(i)))
         end do
         guesses = guesses(:n)
      end if
      call new_file_checkpoint(store, file, guesses, files%interval)
   end subroutine new_vv_checkpoint

   !> UV file store: `file` or the automatic name in `directory`, guessing from `guess` or the saved file
   !>
   !> @param[in]  files      UV checkpoint settings
   !> @param[in]  automatic  automatic file name of the model
   !> @param[out] store      file store
   subroutine new_uv_checkpoint(files, automatic, store)
      !> UV checkpoint settings
      type(moz_checkpoint_files_type), intent(in) :: files
      !> Automatic file name
      character(len=*), intent(in) :: automatic
      !> File store
      type(moz_file_checkpoint_type), intent(out) :: store
      character(len=:), allocatable :: file

      if (allocated(files%file)) then
         file = trim(files%file)
      else
         file = automatic_file(files, automatic)
      end if
      if (allocated(files%guess)) then
         call new_file_checkpoint(store, file, [character(len=len_trim(files%guess)) :: trim(files%guess)], &
            & files%interval)
      else
         call new_file_checkpoint(store, file, interval=files%interval)
      end if
   end subroutine new_uv_checkpoint

   !> Automatic VV checkpoint name, vv_rism_solvent_<id>_<theory>_<closure>.mbin
   !>
   !> @param[in] solvent_id  solvent id; 0 for a custom solvent
   !> @param[in] theory      bulk theory name
   !> @param[in] closure     closure name
   pure function vv_checkpoint_name(solvent_id, theory, closure) result(file)
      !> Solvent id
      integer, intent(in) :: solvent_id
      !> Bulk theory name
      character(len=*), intent(in) :: theory
      !> Closure name
      character(len=*), intent(in) :: closure
      !> File name
      character(len=:), allocatable :: file
      character(len=16) :: id

      write (id, "(i0)") solvent_id
      file = "vv_rism_solvent_"//trim(id)//"_"//to_lower(trim(adjustl(theory)))//"_"// &
         & to_lower(trim(adjustl(closure)))//".mbin"
   end function vv_checkpoint_name

   !> Automatic file name in the directory of the settings
   !>
   !> @param[in] files  checkpoint settings
   !> @param[in] name   automatic file name
   pure function automatic_file(files, name) result(file)
      !> Checkpoint settings
      type(moz_checkpoint_files_type), intent(in) :: files
      !> Automatic file name
      character(len=*), intent(in) :: name
      !> File path
      character(len=:), allocatable :: file
      character(len=:), allocatable :: dir

      file = name
      if (.not. allocated(files%directory)) return
      dir = trim(files%directory)
      if (len(dir) == 0) return
      if (dir(len(dir):) == "/") then
         file = dir//name
      else
         file = dir//"/"//name
      end if
   end function automatic_file

   !> Leave the solvent to the host, at the solvent's default state point
   !>
   !> @param[inout] self  solvent state
   subroutine state_defaults(self)
      !> Solvent state
      class(moz_solvent_state_parameters_type), intent(inout) :: self
      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%temperature)) deallocate (self%temperature)
      if (allocated(self%pressure)) deallocate (self%pressure)
   end subroutine state_defaults

   !> Register the solvent name, temperature and pressure
   !>
   !> @param[inout] self  solvent state
   subroutine state_entries(self)
      !> Solvent state
      class(moz_solvent_state_parameters_type), intent(inout), target :: self
      call self%register_alloc_string("name", self%name, level=2)
      call self%register_alloc_real_scalar("temperature", self%temperature, level=2)
      call self%register_alloc_real_scalar("pressure", self%pressure, level=2)
   end subroutine state_entries

   !> Require a nonblank name and a finite positive temperature and pressure
   !>
   !> @param[inout] self   solvent state
   !> @param[out]   error  invalid name or state point
   subroutine state_validate(self, error)
      !> Solvent state
      class(moz_solvent_state_parameters_type), intent(inout) :: self
      !> Invalid name or state point
      type(error_type), allocatable, intent(out) :: error
      if (is_blank(self%name)) then
         call fatal_error(error, "Solvent name must not be blank")
      else if (.not. positive_or_unset(self%temperature)) then
         call fatal_error(error, "Solvent temperature must be positive (K)")
      else if (.not. positive_or_unset(self%pressure)) then
         call fatal_error(error, "Solvent pressure must be positive (Pa)")
      end if
   end subroutine state_validate

   !> Whether an optional value is unset or finite and positive
   !>
   !> @param[in] value  optional value
   pure function positive_or_unset(value) result(valid)
      !> Optional value
      real(wp), allocatable, intent(in) :: value
      !> Unset, or finite and positive
      logical :: valid
      valid = .true.
      if (allocated(value)) valid = ieee_is_finite(value) .and. value > 0.0_wp
   end function positive_or_unset

   !> Restore automatic checkpoints
   !>
   !> @param[inout] self  checkpoint settings
   subroutine checkpoint_defaults(self)
      !> Checkpoint settings
      class(moz_checkpoint_parameters_type), intent(inout) :: self
      type(moz_checkpoint_parameters_type) :: defaults
      self%enabled = defaults%enabled
      call self%vv%init_defaults()
      call self%uv%init_defaults()
   end subroutine checkpoint_defaults

   !> Register the switch and the VV and UV file groups
   !>
   !> @param[inout] self  checkpoint settings
   subroutine checkpoint_entries(self)
      !> Checkpoint settings
      class(moz_checkpoint_parameters_type), intent(inout), target :: self
      call self%register_logical("enabled", self%enabled, level=2)
      call self%register_group("vv", self%vv)
      call self%register_group("uv", self%uv)
   end subroutine checkpoint_entries

   !> Check both file groups
   !>
   !> @param[inout] self   checkpoint settings
   !> @param[out]   error  first invalid setting
   subroutine checkpoint_validate(self, error)
      !> Checkpoint settings
      class(moz_checkpoint_parameters_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      call self%vv%validate(error)
      if (.not. allocated(error)) call self%uv%validate(error)
   end subroutine checkpoint_validate

   !> Restore automatic names in the working directory and the default interval
   !>
   !> @param[inout] self  checkpoint files
   subroutine files_defaults(self)
      !> Checkpoint files
      class(moz_checkpoint_files_type), intent(inout) :: self
      type(moz_checkpoint_files_type) :: defaults
      if (allocated(self%file)) deallocate (self%file)
      if (allocated(self%guess)) deallocate (self%guess)
      if (allocated(self%directory)) deallocate (self%directory)
      self%interval = defaults%interval
   end subroutine files_defaults

   !> Register the saved file, the guess file, the directory and the interval
   !>
   !> @param[inout] self  checkpoint files
   subroutine files_entries(self)
      !> Checkpoint files
      class(moz_checkpoint_files_type), intent(inout), target :: self
      call self%register_alloc_string("file", self%file, level=2)
      call self%register_alloc_string("guess", self%guess, level=2)
      call self%register_alloc_string("directory", self%directory, level=2)
      call self%register_int_scalar("interval", self%interval)
   end subroutine files_entries

   !> Require nonblank names and a nonnegative interval
   !>
   !> @param[inout] self   checkpoint files
   !> @param[out]   error  invalid name or interval
   subroutine files_validate(self, error)
      !> Checkpoint files
      class(moz_checkpoint_files_type), intent(inout) :: self
      !> Invalid name or interval
      type(error_type), allocatable, intent(out) :: error
      if (self%interval < 0) then
         call fatal_error(error, "Checkpoint interval must be >= 0 (0 saves only the final field)")
         return
      end if
      if (is_blank(self%file)) then
         call fatal_error(error, "Checkpoint file name must not be blank")
      else if (is_blank(self%guess)) then
         call fatal_error(error, "Checkpoint guess file name must not be blank")
      else if (is_blank(self%directory)) then
         call fatal_error(error, "Checkpoint directory must not be blank")
      end if
   end subroutine files_validate

   !> Whether an optional name is set but blank
   !>
   !> @param[in] name  optional name
   pure function is_blank(name) result(blank)
      !> Optional name
      character(len=:), allocatable, intent(in) :: name
      !> Set and blank
      logical :: blank
      blank = .false.
      if (allocated(name)) blank = len_trim(name) == 0
   end function is_blank

end module moist_model_moz_input
