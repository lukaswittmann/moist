!> Continuum component constructors and their common typed interface
module moist_model_continuum_component

   use moist_model_continuum_component_gostshyp, only: model_continuum_component_gostshyp, &
      & new_component_gostshyp, moist_gostshyp_parameters_type
   use moist_model_continuum_component_pcm, only: model_continuum_component_cpcm, &
      & new_component_cpcm, model_continuum_component_cosmo, new_component_cosmo, &
      & solver_type, moist_pcm_parameters_type
   use moist_model_continuum_component_pv, only: model_continuum_component_pv, new_component_pv

   use moist_model_continuum_component_type, only: model_continuum_component_type
   implicit none(type, external)
   private
   public :: model_continuum_component_type

   public :: model_continuum_component_cpcm, new_component_cpcm
   public :: model_continuum_component_cosmo, new_component_cosmo
   public :: solver_type, moist_pcm_parameters_type
   public :: model_continuum_component_gostshyp, new_component_gostshyp, moist_gostshyp_parameters_type
   public :: model_continuum_component_pv, new_component_pv

end module moist_model_continuum_component
