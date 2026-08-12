!> Re-export of all solvation models
module moist_model
   use moist_type, only: solvation_model_type
   use moist_model_general, only: solvation_model_general, new_model_general
   implicit none
   private

   public :: solvation_model_type
   public :: solvation_model_general, new_model_general

end module moist_model
