!> Solvents at a thermodynamic state and their bulk descriptors
module moist_solvent
   use moist_solvent_descriptors, only: solvent_descriptors_type, get_solvent_descriptors, &
      & get_solvent_descriptors_by_id
   use moist_solvent_system, only: solvent_system_type, new_solvent_system
   implicit none(type, external)
   private

   public :: solvent_descriptors_type, get_solvent_descriptors, get_solvent_descriptors_by_id
   public :: solvent_system_type, new_solvent_system
end module moist_solvent
