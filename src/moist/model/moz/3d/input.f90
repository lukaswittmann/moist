!> RISM input file of a 3D UV model and its constructor overload
!>
!> - UV keys at the root, nested `[solvent]` VV settings; shared groups from moist_model_moz_input
!> - `new_moz_3d_model(model, ctx, input, error)`: solvent from `[solvent]`, grid, solute potential
!>   and checkpoint stores from the input; the model still needs `update(mol)`
!> - `new_moz_3d_model(model, ctx, input, system, error)`: the same with a host solvent system
module moist_model_moz_3d_input
   use mctc_env, only: error_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_factory, only: new_spatial_grid
   use moist_model_moz_potential_factory, only: configure_potential
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type
   use moist_model_moz_input, only: model_moz_input_type, moz_input_defaults, moz_input_entries, &
      & moz_input_validate, new_input_solvent, new_uv_checkpoint
   use moist_model_moz_3d_parameters, only: model_moz_3d_parameters_type
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_solvent_system, only: solvent_system_type
   implicit none(type, external)
   private
   public :: model_moz_3d_input_type, new_moz_3d_model

   !> Automatic name of the UV checkpoint file
   character(len=*), parameter :: uv_checkpoint_file = "uv_rism_3d.mbin"

   !> Settings of a 3D RISM input file
   type, extends(model_moz_input_type) :: model_moz_3d_input_type
      !> Settings of the UV and VV solves
      type(model_moz_3d_parameters_type) :: uv
   contains
      procedure :: init_defaults => input_3d_defaults
      procedure :: register_entries => input_3d_entries
      procedure :: validate => input_3d_validate
      procedure :: new_model_system => input_3d_new_model
   end type model_moz_3d_input_type

   !> Construct from an input file's settings, with the solvent of the file or of the host
   interface new_moz_3d_model
      module procedure new_moz_3d_model_file
      module procedure new_moz_3d_model_input
   end interface new_moz_3d_model

contains

   !> Construct a 3D UV model from input settings with the `[solvent]` solvent; update it with the solute next
   !>
   !> @param[out] self   model
   !> @param[in]  ctx    borrowed context, outlives the model
   !> @param[in]  input  input settings with restored defaults and a solvent name
   !> @param[out] error  missing solvent, invalid settings or construction failure
   subroutine new_moz_3d_model_file(self, ctx, input, error)
      !> Model
      type(model_moz_3d_type), intent(out) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Input settings with restored defaults and a solvent name
      type(model_moz_3d_input_type), intent(in) :: input
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(solvent_system_type) :: system

      call input%new_system(system, error)
      if (.not. allocated(error)) call new_moz_3d_model(self, ctx, input, system, error)
   end subroutine new_moz_3d_model_file

   !> Construct a 3D UV model from input settings and a host solvent system; update it with the solute next
   !>
   !> - The settings are validated on a copy; the `[solvent]` name and state point are unused
   !>
   !> @param[out] self    model
   !> @param[in]  ctx     borrowed context, outlives the model
   !> @param[in]  input   input settings with restored defaults
   !> @param[in]  system  solvent at its state point
   !> @param[out] error   invalid settings, solvent, grid or construction failure
   subroutine new_moz_3d_model_input(self, ctx, input, system, error)
      !> Model
      type(model_moz_3d_type), intent(out) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Input settings with restored defaults
      type(model_moz_3d_input_type), intent(in) :: input
      !> Solvent at its state point
      type(solvent_system_type), intent(in) :: system
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(model_moz_3d_input_type) :: settings
      type(solvent_vv_type) :: solvent
      class(moist_math_grid_3d_type), allocatable :: grid
      type(moz_file_checkpoint_type) :: store

      settings = input
      call settings%validate(error)
      if (.not. allocated(error)) call new_input_solvent(settings, ctx, system, settings%uv%solvent, solvent, error)
      if (.not. allocated(error)) call new_spatial_grid(settings%grid, grid, error)
      if (allocated(error)) return
      call new_moz_3d_model(self, ctx, grid, solvent, error, settings%uv)
      if (.not. allocated(error)) call configure_potential(self%potential, settings%potential, error, &
         & mixing=solvent%potential%mixing())
      if (allocated(error)) return
      if (settings%checkpoint%enabled) then
         call new_uv_checkpoint(settings%checkpoint%uv, uv_checkpoint_file, store)
         call self%set_checkpoint(store)
      end if
   end subroutine new_moz_3d_model_input

   !> Restore the 3D defaults
   !>
   !> @param[inout] self  input settings
   subroutine input_3d_defaults(self)
      !> Input settings
      class(model_moz_3d_input_type), intent(inout) :: self
      call self%uv%init_defaults()
      call moz_input_defaults(self)
   end subroutine input_3d_defaults

   !> Register the UV keys at the root, then the shared groups
   !>
   !> @param[inout] self  input settings
   subroutine input_3d_entries(self)
      !> Input settings
      class(model_moz_3d_input_type), intent(inout), target :: self
      call self%register_group("", self%uv)
      call moz_input_entries(self)
   end subroutine input_3d_entries

   !> Validate the shared and the UV settings
   !>
   !> @param[inout] self   input settings
   !> @param[out]   error  first invalid setting
   subroutine input_3d_validate(self, error)
      !> Input settings
      class(model_moz_3d_input_type), intent(inout) :: self
      !> First invalid setting
      type(error_type), allocatable, intent(out) :: error
      call moz_input_validate(self, error)
      if (.not. allocated(error)) call self%uv%validate(error)
   end subroutine input_3d_validate

   !> Construct the 3D UV model of the input with a host solvent system, not yet updated
   !>
   !> @param[in]  self    input settings
   !> @param[in]  ctx     borrowed context, outlives the model
   !> @param[in]  system  solvent at its state point
   !> @param[out] model   constructed model
   !> @param[out] error   invalid settings, solvent or construction failure
   subroutine input_3d_new_model(self, ctx, system, model, error)
      !> Input settings
      class(model_moz_3d_input_type), intent(in) :: self
      !> Borrowed context, outlives the model
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent at its state point
      type(solvent_system_type), intent(in) :: system
      !> Constructed model
      class(solvation_model_type), allocatable, intent(out) :: model
      !> Construction failure
      type(error_type), allocatable, intent(out) :: error
      type(model_moz_3d_type), allocatable :: uv

      allocate (uv)
      call new_moz_3d_model(uv, ctx, self, system, error)
      if (.not. allocated(error)) call move_alloc(uv, model)
   end subroutine input_3d_new_model

end module moist_model_moz_3d_input
