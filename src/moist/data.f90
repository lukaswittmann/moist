module moist_data
   use moist_data_radii_legacy, only: get_radius
   use moist_data_en, only: get_electronegativity
   use moist_data_hardness, only: get_hardness
   use moist_data_solvents, only: solvent_system_type, new_solvent_system
   use moist_data_solvents, only: get_solvent_id
   use moist_data_solvents, only: solvent_descriptors_type, get_solvent_descriptors, &
      & get_solvent_descriptors_by_id
   implicit none(type, external)
   private

   public :: get_radius
   public :: get_electronegativity
   public :: get_hardness
   public :: solvent_system_type, new_solvent_system, get_solvent_id
   public :: solvent_descriptors_type, get_solvent_descriptors, get_solvent_descriptors_by_id
end module moist_data
