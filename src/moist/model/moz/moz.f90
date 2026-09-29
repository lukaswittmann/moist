!> Molecular Ornstein-Zernike model families
module moist_model_moz
   use moist_model_moz_1d, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_3d, only: model_moz_3d_type, new_moz_3d_model
   implicit none(type, external)
   private
   public :: model_moz_1d_type, new_moz_1d_model
   public :: model_moz_3d_type, new_moz_3d_model
end module moist_model_moz
