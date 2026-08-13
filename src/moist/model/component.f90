module moist_model_components

   use moist_model_component_gostshyp, only: solvation_model_component_gostshyp, &
      & new_component_gostshyp
   use moist_model_component_pcm, only: solvation_model_component_cpcm, &
      & new_component_cpcm, solvation_model_component_cosmo, new_component_cosmo, &
      & solver_type, potential_source
   use moist_model_component_pv, only: solvation_model_component_pv, new_component_pv

   implicit none
   private

   public :: solvation_model_component_cpcm, new_component_cpcm
   public :: solvation_model_component_cosmo, new_component_cosmo
   public :: solver_type, potential_source
   public :: solvation_model_component_gostshyp, new_component_gostshyp
   public :: solvation_model_component_pv, new_component_pv

end module moist_model_components
