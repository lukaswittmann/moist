!> Construction of a grid model from the settings tables
!>
!> One call builds the model a cavity takes as its optional `grid_model`
!> argument; the cavity never knows which concrete model, or combination of
!> models, it received
module moist_cavity_grid_model_factory
   use mctc_env, only: error_type
   use moist_cavity_grid_model_type, only: grid_model_type
   use moist_cavity_grid_model_symmetry_parameters, only: moist_cavity_grid_model_symmetry_parameters_type
   use moist_cavity_grid_model_symmetry_type, only: grid_model_symmetry_type, new_grid_model_symmetry
   implicit none(type, external)
   private

   public :: new_grid_model

contains

   !> Build the grid model of the given settings
   !>
   !> - `symmetry`: a [[grid_model_symmetry_type]] (sphere orientations by molecular
   !>   frame and point group)
   !> - no settings: `model` stays unallocated, which a cavity takes as no
   !>   model at all
   !>
   !> @param[out] model    Grid model, nothing computed
   !> @param[out] error    Invalid setting
   !> @param[in]  symmetry Settings of the symmetry model
   subroutine new_grid_model(model, error, symmetry)
      !> Grid model
      class(grid_model_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Symmetry settings
      type(moist_cavity_grid_model_symmetry_parameters_type), intent(in), optional :: symmetry

      !> Symmetry model under construction
      type(grid_model_symmetry_type), allocatable :: sym

      if (.not. present(symmetry)) return
      allocate (sym)
      call new_grid_model_symmetry(sym, error, param=symmetry)
      if (allocated(error)) return
      call move_alloc(sym, model)

   end subroutine new_grid_model

end module moist_cavity_grid_model_factory
