!> Public continuum and MOZ model families
module moist_model
   use moist_model_type, only: solvation_model_type
   use moist_model_continuum, only: model_continuum_type, new_continuum_model
   use moist_model_moz, only: model_moz_1d_type, new_moz_1d_model, &
      & model_moz_3d_type, new_moz_3d_model
   implicit none(type, external)
   private
   public :: solvation_model_type
   public :: model_continuum_type, new_continuum_model
   public :: model_moz_1d_type, new_moz_1d_model, model_moz_3d_type, new_moz_3d_model
end module moist_model
