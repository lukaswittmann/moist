!> Water recognition and rigid geometry shared by SPC/E and cSPC/E
module moist_model_moz_solvent_water
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use moist_data_solvents, only: get_solvent_id
   implicit none(type, external)
   private
   public :: new_spce_water_geometry, is_water

contains

   !> Construct neutral O,H,H water with 1 Angstrom bonds and a 109.47 degree angle
   !>
   !> @param[out] mol Rigid SPC/E or cSPC/E water structure, in bohr
   subroutine new_spce_water_geometry(mol)
      !> Water structure
      type(structure_type), intent(out) :: mol
      real(wp), parameter :: half_angle = 109.47_wp*acos(-1.0_wp)/360.0_wp
      real(wp), parameter :: x = aatoau*sin(half_angle), z = aatoau*cos(half_angle)

      call new(mol, [8, 1, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, x, 0.0_wp, z, -x, 0.0_wp, z], [3, 3]))
   end subroutine new_spce_water_geometry

   !> Whether the structure is one neutral, closed-shell water molecule
   !>
   !> - Non-water table solvent id: false
   !>
   !> @param[in] mol         structure with initialized atom identities
   !> @param[in] solvent_id  optional solvent id of the table; 0 for a custom solvent
   function is_water(mol, solvent_id) result(water)
      !> Structure
      class(structure_type), intent(in) :: mol
      !> Optional solvent id of the table
      integer, intent(in), optional :: solvent_id
      !> Whether the structure is neutral H2O
      logical :: water
      type(error_type), allocatable :: error
      integer :: water_id
      water = .false.
      if (present(solvent_id)) then
         if (solvent_id /= 0) then
            call get_solvent_id("water", water_id, error)
            if (allocated(error)) return
            if (solvent_id /= water_id) return
         end if
      end if
      if (mol%nat /= 3) return
      if (count(mol%num(mol%id) == 8) /= 1 .or. count(mol%num(mol%id) == 1) /= 2) return
      if (.not. ieee_is_finite(mol%charge) .or. abs(mol%charge) > 1.0e-8_wp .or. mol%uhf /= 0) return
      water = .true.
   end function is_water

end module moist_model_moz_solvent_water
