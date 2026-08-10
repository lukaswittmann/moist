!> Re-export of all solvation models
module moist_model
   use moist_type, only: solvation_model
   use moist_model_general, only: general_solvation_model, new_general_model
   implicit none
   private

   public :: solvation_model
   public :: general_solvation_model, new_general_model

end module moist_model
