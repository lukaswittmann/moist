!> Continuum models and components evaluated on cavity surfaces
module moist_model_continuum
   use moist_model_continuum_type, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component, only: model_continuum_component_type, &
      & model_continuum_component_cpcm, new_component_cpcm, &
      & model_continuum_component_cosmo, new_component_cosmo, &
      & model_continuum_component_pv, new_component_pv, &
      & model_continuum_component_gostshyp, new_component_gostshyp, &
      & solver_type, moist_pcm_parameters_type, moist_gostshyp_parameters_type
   implicit none(type, external)
   private
   public :: model_continuum_type, new_continuum_model, model_continuum_component_type
   public :: model_continuum_component_cpcm, new_component_cpcm
   public :: model_continuum_component_cosmo, new_component_cosmo
   public :: model_continuum_component_pv, new_component_pv
   public :: model_continuum_component_gostshyp, new_component_gostshyp
   public :: solver_type, moist_pcm_parameters_type, moist_gostshyp_parameters_type
end module moist_model_continuum
