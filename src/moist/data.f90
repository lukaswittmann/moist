module moist_data
   use moist_data_radii_legacy, only: get_radius
   use moist_data_en, only: get_electronegativity
   use moist_data_hardness, only: get_hardness
   use moist_data_solvents, only: solvation_system_type, new_solvation_system
   use moist_data_solvents, only: get_solvent_id
   implicit none(type, external)
   private

   public :: get_radius
   public :: get_electronegativity
   public :: get_hardness
   public :: solvation_system_type, new_solvation_system, get_solvent_id
end module moist_data
