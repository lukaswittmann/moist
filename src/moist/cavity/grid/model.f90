!> Grid models: providers of the components of every atomic sphere
!>
!> - [[moist_cavity_grid_model_type]]: abstract model and its components
!> - [[moist_cavity_grid_model_factory]]: one constructor from the settings tables
!> - [[moist_cavity_grid_model_symmetry]]: point-group symmetry model
module moist_cavity_grid_model
   use moist_cavity_grid_model_type, only: grid_model_type, grid_components_type, grid_size_default
   use moist_cavity_grid_model_factory, only: new_grid_model
   use moist_cavity_grid_model_symmetry, only: grid_model_symmetry_type, new_grid_model_symmetry, &
      & moist_cavity_grid_model_symmetry_parameters_type, check_symmetry_parameters
   implicit none(type, external)
   private

   public :: grid_model_type, grid_components_type, grid_size_default
   public :: new_grid_model
   public :: grid_model_symmetry_type, new_grid_model_symmetry
   public :: moist_cavity_grid_model_symmetry_parameters_type, check_symmetry_parameters

end module moist_cavity_grid_model
