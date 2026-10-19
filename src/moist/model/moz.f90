!> Molecular Ornstein-Zernike model families and the bulk solvent
!>
!> - Potential terms from `moist_model_moz_potential`
module moist_model_moz
   use moist_model_moz_1d_type, only: model_moz_1d_type, new_moz_1d_model
   use moist_model_moz_3d_type, only: model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_solvent_vv, only: solvent_vv_type, new_vv_solvent
   implicit none(type, external)
   private
   public :: model_moz_1d_type, new_moz_1d_model
   public :: model_moz_3d_type, new_moz_3d_model
   public :: solvent_vv_type, new_vv_solvent
end module moist_model_moz
