!> Angular grids of the atomic spheres: components, models and Jacobians
!>
!> - [[moist_cavity_grid_model]]: providers of the components of every sphere
!> - [[moist_cavity_grid_type]]: flat container of the components a cavity loops over
!> - [[moist_cavity_grid_frame]]: rotation and weight Jacobians for forward-mode terms
module moist_cavity_grid
   use moist_cavity_grid_model, only: grid_model_type, grid_components_type, grid_size_default, &
      & new_grid_model, grid_model_symmetry_type, new_grid_model_symmetry, &
      & moist_cavity_grid_model_symmetry_parameters_type, check_symmetry_parameters
   use moist_cavity_grid_type, only: cavity_grid_type
   use moist_cavity_grid_frame, only: grid_frame_type, new_grid_frame
   implicit none(type, external)
   private

   public :: grid_model_type, grid_components_type, grid_size_default
   public :: cavity_grid_type, grid_frame_type, new_grid_frame, new_grid_model
   public :: grid_model_symmetry_type, new_grid_model_symmetry
   public :: moist_cavity_grid_model_symmetry_parameters_type, check_symmetry_parameters

end module moist_cavity_grid
