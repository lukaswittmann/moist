!> Point-group symmetry model: settings, construction and the model type
!>
!> - [[moist_cavity_grid_model_symmetry_parameters]]: the top-level `symmetry` table
!> - [[moist_cavity_grid_model_symmetry_type]]: detection, orientations, report and fields
module moist_cavity_grid_model_symmetry
   use moist_cavity_grid_model_symmetry_parameters, only: moist_cavity_grid_model_symmetry_parameters_type, &
      & check_symmetry_parameters
   use moist_cavity_grid_model_symmetry_type, only: grid_model_symmetry_type, new_grid_model_symmetry
   implicit none(type, external)
   private

   public :: moist_cavity_grid_model_symmetry_parameters_type, check_symmetry_parameters
   public :: grid_model_symmetry_type, new_grid_model_symmetry

end module moist_cavity_grid_model_symmetry
