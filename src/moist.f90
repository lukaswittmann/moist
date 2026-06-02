
module moist
   use mctc_io, only: structure_type, new
   use mctc_env, only: error_type, fatal_error, wp
   use moist_version, only: get_moist_version
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_radii, only: radius_type, new_radii
   use moist_build_info

   use moist_model, only: solvation_model
   ! use moist_model_gems, only : gems_model, new_gems_model
   implicit none
   public

end module moist
