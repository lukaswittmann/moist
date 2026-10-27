!> RISM input file of a 1D UV model and its constructor overload
!>
!> - UV keys at the root, nested `[solvent]` VV settings; shared groups from moist_model_moz_input
!> - Radial grid defaults; EC charges refused, their fit needs a volume grid
!> - `new_moz_1d_model(model, ctx, input, error)`: solvent from `[solvent]`, radial grids, solute
!>   potential and checkpoint stores from the input; the model still needs `update(mol)`
!> - `new_moz_1d_model(model, ctx, input, system, error)`: the same with a host solvent system
module moist_model_moz_1d_input
   use mctc_env, only: error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_uniform_radial_pair
   use moist_model_moz_potential_factory, only: configure_potential
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type
   use moist_model_moz_input, only: model_moz_input_type, moz_input_defaults, moz_input_entries, &
      & moz_input_validate, new_input_solvent, new_uv_checkpoint
   use moist_model_moz_1d_parameters, only: model_moz_1d_parameters_type
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_solvent_system, only: solvent_system_type
   implicit none(type, external)
   private
   public :: model_moz_1d_input_type, new_moz_1d_model

   !> Automatic name of the UV checkpoint file
   character(len=*), parameter :: uv_checkpoint_file = "uv_rism_1d.mbin"

   !> Settings of a 1D RISM input file
   type, extends(model_moz_input_type) :: model_moz_1d_input_type
      !> Settings of the UV and VV solves
      type(model_moz_1d_parameters_type) :: uv
   contains
      procedure :: init_defaults => input_1d_defaults
      procedure :: register_entries => input_1d_entries
      procedure :: validate => input_1d_validate
      procedure :: new_model_system => input_1d_new_model
   end type model_moz_1d_input_type

   !> Construct from an input file's settings, with the solvent of the file or of the host
   interface new_moz_1d_model
      module procedure new_moz_1d_model_file
      module procedure new_moz_1d_model_input
   end interface new_moz_1d_model

contains

   !> Construct a 1D UV model from input settings with the `[solvent]` solvent; update it with the solute next
   !>
   !> @param[out] self   model
   !> @param[in]  ctx    borrowed context, outlives the model
   !> @param[in]  input  input settings with restored defaults and a solvent name
   !> @param[out] error  missing solvent, invalid settings or construction failure
   subroutine new_moz_1d_model_file(self, ctx, input, error)
      !> Model
      type(model_moz_1d_type), intent(out) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Input settings with restored defaults and a solvent name
      type(model_moz_1d_input_type), intent(in) :: input
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(solvent_system_type) :: system

      call input%new_system(system, error)
      if (.not. allocated(error)) call new_moz_1d_model(self, ctx, input, system, error)
   end subroutine new_moz_1d_model_file

   !> Construct a 1D UV model from input settings and a host solvent system; update it with the solute next
   !>
   !> - The settings are validated on a copy; the `[solvent]` name and state point are unused
   !>
   !> @param[out] self    model
   !> @param[in]  ctx     borrowed context, outlives the model
   !> @param[in]  input   input settings with restored defaults
   !> @param[in]  system  solvent at its state point
   !> @param[out] error   invalid settings, solvent, grid or construction failure
   subroutine new_moz_1d_model_input(self, ctx, input, system, error)
      !> Model
      type(model_moz_1d_type), intent(out) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Input settings with restored defaults
      type(model_moz_1d_input_type), intent(in) :: input
      !> Solvent at its state point
      type(solvent_system_type), intent(in) :: system
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(model_moz_1d_input_type) :: settings
      type(solvent_vv_type) :: solvent
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moz_file_checkpoint_type) :: store

      settings = input
      call settings%validate(error)
      if (.not. allocated(error)) call new_input_solvent(settings, ctx, system, settings%uv%solvent, solvent, error)
      if (.not. allocated(error)) call new_uniform_radial_pair(rgrid, kgrid, settings%grid%npts, settings%grid%dr, error)
      if (allocated(error)) return
      call new_moz_1d_model(self, ctx, rgrid, kgrid, solvent, error, settings%uv)
      if (.not. allocated(error)) call configure_potential(self%potential, settings%potential, error, &
         & mixing=solvent%potential%mixing())
      if (allocated(error)) return
      if (settings%checkpoint%enabled) then
         call new_uv_checkpoint(settings%checkpoint%uv, uv_checkpoint_file, store)
         call self%set_checkpoint(store)
      end if
   end subroutine new_moz_1d_model_input

   !> Restore the 1D defaults, radial grid included
   !>
   !> @param[inout] self  input settings
   subroutine input_1d_defaults(self)
      !> Input settings
      class(model_moz_1d_input_type), intent(inout) :: self
      call self%uv%init_defaults()
      call self%grid%set_dimension(1)
      call moz_input_defaults(self)
   end subroutine input_1d_defaults

   !> Register the UV keys at the root, then the shared groups
   !>
   !> @param[inout] self  input settings
   subroutine input_1d_entries(self)
      !> Input settings
      class(model_moz_1d_input_type), intent(inout), target :: self
      call self%register_group("", self%uv)
      call moz_input_entries(self)
   end subroutine input_1d_entries

   !> Validate the shared and the UV settings; EC charges need a volume grid
   !>
   !> @param[inout] self   input settings
   !> @param[out]   error  first invalid setting
   subroutine input_1d_validate(self, error)
      !> Input settings
      class(model_moz_1d_input_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      call moz_input_validate(self, error)
      if (.not. allocated(error)) call self%uv%validate(error)
      if (allocated(error)) return
      if (any(self%potential%charges == "ec")) then
         call fatal_error(error, "EC charges need a volume grid; use the 3D model")
      end if
   end subroutine input_1d_validate

   !> Construct the 1D UV model of the input with a host solvent system, not yet updated
   !>
   !> @param[in]  self    input settings
   !> @param[in]  ctx     borrowed context, outlives the model
   !> @param[in]  system  solvent at its state point
   !> @param[out] model   constructed model
   !> @param[out] error   invalid settings, solvent or construction failure
   subroutine input_1d_new_model(self, ctx, system, model, error)
      !> Input settings
      class(model_moz_1d_input_type), intent(in) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent at its state point
      type(solvent_system_type), intent(in) :: system
      !> Constructed model
      class(solvation_model_type), allocatable, intent(out) :: model
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(model_moz_1d_type), allocatable :: uv

      allocate (uv)
      call new_moz_1d_model(uv, ctx, self, system, error)
      if (.not. allocated(error)) call move_alloc(uv, model)
   end subroutine input_1d_new_model

end module moist_model_moz_1d_input
