!> Legacy SPC/E-style water Lennard-Jones parameters
!>
!> Berendsen et al., J. Phys. Chem. 91, 6269 (1987), doi:10.1021/j100308a038
module moist_model_moz_potential_lj_spce_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kjtoau
   implicit none(type, external)
   private

   public :: lookup_spce_lj

   !> Oxygen sigma, nm
   real(wp), parameter :: spce_sigma_ow_nm = 0.316572_wp
   !> Hydrogen sigma, nm
   real(wp), parameter :: spce_sigma_hw_nm = 0.116572_wp
   !> Oxygen epsilon, kJ/mol
   real(wp), parameter :: spce_epsilon_ow_kj = 0.6497752516070757_wp
   !> Hydrogen epsilon, kJ/mol
   real(wp), parameter :: spce_epsilon_hw_kj = 0.06497752516070757_wp

contains

   !> LJ parameters of a water site by element
   !>
   !> @param[in] number Atomic number, 8 for oxygen or 1 for hydrogen
   !> @param[out] sigma Sigma, bohr; zero when not found
   !> @param[out] epsilon Epsilon, Hartree; zero when not found
   !> @param[out] found Whether the element is a water site
   pure subroutine lookup_spce_lj(number, sigma, epsilon, found)
      !> Atomic number
      integer, intent(in) :: number
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Epsilon, Hartree
      real(wp), intent(out) :: epsilon
      !> Whether the element is a water site
      logical, intent(out) :: found

      sigma = 0.0_wp
      epsilon = 0.0_wp
      found = .true.
      select case (number)
      case (8)
         sigma = spce_sigma_ow_nm*10.0_wp*aatoau
         epsilon = spce_epsilon_ow_kj*kjtoau
      case (1)
         sigma = spce_sigma_hw_nm*10.0_wp*aatoau
         epsilon = spce_epsilon_hw_kj*kjtoau
      case default
         found = .false.
      end select
   end subroutine lookup_spce_lj

end module moist_model_moz_potential_lj_spce_parameters
