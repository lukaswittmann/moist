!> Transition-metal element Lennard-Jones parameters
module moist_model_moz_potential_lj_tm_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kjtoau
   implicit none(type, external)
   private

   public :: lookup_tm_lj, n_tm

   !> Highest supported atomic number (Hg); the tables are indexed by atomic number
   integer, parameter :: n_tm = 80

   !> Sentinel of an unsupported element
   real(wp), parameter :: miss = -1.0_wp

   !> Sigma per atomic number, nm; `miss` where unsupported
   real(wp), parameter :: tm_sigma_nm(n_tm) = [ &
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 1-10
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 11-20
      & 0.6660_wp, 0.5738_wp, 0.5534_wp, 0.5468_wp, 0.5808_wp, & ! Sc Ti V Cr Mn
      & 0.5938_wp, 0.5610_wp, 0.5440_wp, 0.5336_wp, 0.5694_wp, & ! Fe Co Ni Cu Zn
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 31-40
      & miss, miss, miss, 0.5878_wp, 0.5430_wp, miss, miss, miss, miss, miss, & ! 41-50, Ru Rh
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 51-60
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 61-70
      & miss, miss, miss, miss, miss, miss, miss, 0.5340_wp, miss, 0.5602_wp & ! 71-80, Pt Hg
      & ]

   !> Epsilon per atomic number, kJ/mol; `miss` where unsupported
   real(wp), parameter :: tm_epsilon_kj(n_tm) = [ &
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 1-10
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 11-20
      & 0.46024_wp, 3.83672804_wp, 7.966336_wp, 6.351312_wp, 4.556376_wp, & ! Sc Ti V Cr Mn
      & 2.920432_wp, 5.004064_wp, 11.0876_wp, 8.987232_wp, 4.58148_wp, & ! Fe Co Ni Cu Zn
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 31-40
      & miss, miss, miss, 1.748912_wp, 16.510064_wp, miss, miss, miss, miss, miss, & ! 41-50, Ru Rh
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 51-60
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 61-70
      & miss, miss, miss, miss, miss, miss, miss, 21.225432_wp, miss, 8.1797201_wp & ! 71-80, Pt Hg
      & ]

contains

   !> Transition-metal Lennard-Jones parameters of an element
   !>
   !> @param[in] number Atomic number
   !> @param[out] sigma Sigma, bohr; zero when unsupported
   !> @param[out] epsilon Well depth, Hartree; zero when unsupported
   !> @param[out] found Whether the element has a row
   pure subroutine lookup_tm_lj(number, sigma, epsilon, found)
      !> Atomic number
      integer, intent(in) :: number
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Well depth, Hartree
      real(wp), intent(out) :: epsilon
      !> Whether the element has a row
      logical, intent(out) :: found

      sigma = 0.0_wp
      epsilon = 0.0_wp
      found = .false.
      if (number < 1 .or. number > n_tm) return
      if (tm_sigma_nm(number) == miss) return
      sigma = tm_sigma_nm(number)*10.0_wp*aatoau
      epsilon = tm_epsilon_kj(number)*kjtoau
      found = .true.
   end subroutine lookup_tm_lj

end module moist_model_moz_potential_lj_tm_parameters
