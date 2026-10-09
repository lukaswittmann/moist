!> Category rules shared by the 1D and 3D potential tables
!>
!> - Categories: short-range pair and electrostatics
!> - Category on one side only: error
!> - Category absent on both sides: no contribution
!> - No category on either side: error
module moist_model_moz_potential_kernel_category
   use mctc_env, only: error_type, fatal_error
   implicit none(type, external)
   private

   public :: resolve_categories

contains

   !> Decide which categories enter the tables from their presence on each side
   !>
   !> @param[in]  pair_u   whether the first side has Lennard-Jones data
   !> @param[in]  pair_v   whether the second side has Lennard-Jones data
   !> @param[in]  elec_u   whether the first side has electrostatics (a field or tail charges)
   !> @param[in]  elec_v   whether the second side has electrostatics
   !> @param[out] do_pair  whether the pair category enters
   !> @param[out] do_elec  whether the electrostatic category enters
   !> @param[out] error    one-sided category or no category
   subroutine resolve_categories(pair_u, pair_v, elec_u, elec_v, do_pair, do_elec, error)
      !> Whether the first side has Lennard-Jones data
      logical, intent(in) :: pair_u
      !> Whether the second side has Lennard-Jones data
      logical, intent(in) :: pair_v
      !> Whether the first side has electrostatics
      logical, intent(in) :: elec_u
      !> Whether the second side has electrostatics
      logical, intent(in) :: elec_v
      !> Whether the pair category enters
      logical, intent(out) :: do_pair
      !> Whether the electrostatic category enters
      logical, intent(out) :: do_elec
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      do_pair = .false.
      do_elec = .false.
      if (.not. (pair_u .or. pair_v .or. elec_u .or. elec_v)) then
         call fatal_error(error, "No potential terms on either side of the interaction")
         return
      end if
      if (pair_u .neqv. pair_v) then
         call fatal_error(error, "Pair potential on one side of the interaction only; "// &
            & "both sides need one or neither")
         return
      end if
      if (elec_u .neqv. elec_v) then
         call fatal_error(error, "Electrostatic potential on one side of the interaction only; "// &
            & "both sides need one or neither")
         return
      end if
      do_pair = pair_u
      do_elec = elec_u
   end subroutine resolve_categories

end module moist_model_moz_potential_kernel_category
