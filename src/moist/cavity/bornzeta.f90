!> Born-fitted Gaussian width scale shared by the iSwiG and DROP surfaces
!>
!> - Surface charge i has the exponent xi_i = zeta / (R_i sqrt(w_i)), with w_i
!>   the Lebedev weight (summing to 4 pi) and R_i the sphere radius
!> - zeta is fitted so that C-PCM gives the exact Born energy of a central
!>   charge in a sphere; the energy scales as 1/R, so one radius fixes it
!> - Values of Table II of J. Chem. Phys. 122, 194110 (2005) in full
!>   precision, as tabulated in PySCF (pyscf/solvent/pcm.py, XI); 6 points
!>   refitted with test/dev/test_cavity_drop_born_fit.f90, where the
!>   published value leaves a Born error of 4.4e-6
module moist_cavity_bornzeta
   use mctc_env, only: wp

   implicit none(type, external)
   private

   public :: born_zeta_nleb, get_born_zeta

   !> Lebedev sizes with a fitted zeta: every rule without negative weights
   integer, parameter :: born_zeta_nleb(29) = [ &
                         6, 14, 26, 38, 50, 86, 110, 146, &
                         170, 194, 302, 350, 434, 590, 770, 974, &
                         1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, &
                         3890, 4334, 4802, 5294, 5810]

   !> Fitted zeta matching `born_zeta_nleb`
   real(wp), parameter :: born_zeta(29) = [ &
                          4.84570840675_wp, 4.86458714334_wp, 4.85478226219_wp, 4.90105812685_wp, &
                          4.89250673295_wp, 4.89741372580_wp, 4.90101060987_wp, 4.89825187392_wp, &
                          4.90685517725_wp, 4.90337644248_wp, 4.90498088169_wp, 4.86879474832_wp, &
                          4.90567349080_wp, 4.90624071359_wp, 4.90656435779_wp, 4.90685167998_wp, &
                          4.90704098216_wp, 4.90721023869_wp, 4.90733270691_wp, 4.90744499142_wp, &
                          4.90753082825_wp, 4.90760972766_wp, 4.90767282394_wp, 4.90773141371_wp, &
                          4.90777965981_wp, 4.90782469526_wp, 4.90749125553_wp, 4.90762073452_wp, &
                          4.90792902522_wp]

contains

   !> Look up the fitted zeta of a Lebedev size
   !>
   !> @param[in]  num_leb Lebedev points per sphere
   !> @param[out] zeta    Fitted zeta; zero when `found` is false
   !> @param[out] found   Whether `num_leb` has a fitted value
   pure subroutine get_born_zeta(num_leb, zeta, found)
      !> Lebedev points per sphere
      integer, intent(in) :: num_leb
      !> Fitted Gaussian width scale
      real(wp), intent(out) :: zeta
      !> Whether a fitted value exists
      logical, intent(out) :: found

      integer :: igrid

      zeta = 0.0_wp
      found = .false.
      do igrid = 1, size(born_zeta_nleb)
         if (born_zeta_nleb(igrid) == num_leb) then
            zeta = born_zeta(igrid)
            found = .true.
            return
         end if
      end do
   end subroutine get_born_zeta

end module moist_cavity_bornzeta
