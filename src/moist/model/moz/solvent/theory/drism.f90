!> Perkyns-Pettitt (1992) dielectric correction for one rigid polar solvent
!>
!> - Port of the Perkyns-Pettitt zeta as implemented in RISMiCal (drism.f / mathlib.f:diprot)
!> - Original centre-of-absolute-charge and dipole projection retained
!> - Complete molecular frame removes orientation ambiguity
!> - Small-argument series give the analytical zeta/k^2 limit
!>
!> zeta_ab = A exp(-smear^2 k^2/2) F_a F_b,
!> F_a = j0(k x_a) j0(k y_a) j1(k z_a),
!> A = ((epsilon - 1)/y - 3)/rho, y = 4 pi beta rho mu^2/9
module moist_model_moz_solvent_theory_drism
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_moz_solvent_theory_type, only: bulk_theory_type, check_bulk_charge_neutrality
   implicit none(type, external)
   private
   public :: drism_theory_type, new_drism_theory

   !> Dielectric correction at a fixed solvent state point; atomic units
   type, extends(bulk_theory_type) :: drism_theory_type
      !> Site coordinates in the dipole-aligned molecular frame, bohr
      real(wp), allocatable, private :: coord(:, :)
      !> Correction amplitude (bohr^3) and Gaussian width (bohr)
      real(wp), private :: amplitude = 0.0_wp, smear = 0.0_wp
   contains
      procedure :: name => drism_name
      procedure :: correction_at => drism_correction_at
   end type drism_theory_type
contains
   !> Prepare the dielectric correction once from the solvent's site charges
   !>
   !> @param[out] self Fixed DRISM theory
   !> @param[in] xyz Site coordinates, bohr (3, ns)
   !> @param[in] q Site charges, e (ns)
   !> @param[in] rho Molecular number density, 1/bohr^3
   !> @param[in] beta Inverse thermal energy, 1/Hartree
   !> @param[in] dielectric Requested relative dielectric constant, at least one
   !> @param[in] smear Nonnegative Gaussian width, bohr
   !> @param[out] error Invalid settings or a vanishing molecular dipole
   subroutine new_drism_theory(self, xyz, q, rho, beta, dielectric, smear, error)
      !> Fixed DRISM theory
      type(drism_theory_type), intent(out) :: self
      !> Site coordinates (bohr), charges (e), density, beta, DRISM dielectric and smear
      real(wp), intent(in) :: xyz(:, :), q(:), rho, beta, dielectric, smear
      !> Invalid settings or vanishing molecular dipole
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: centre(3), dipole(3), ez(3), ex(3), ey(3), transverse(3)
      real(wp) :: qabs, mu, y, radius2, best_radius2
      real(wp), allocatable :: centred(:, :)
      integer :: i

      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= size(q) .or. size(q) < 1) then
         call fatal_error(error, "DRISM requires site coordinates (3, ns) and charges (ns)")
         return
      end if
      if (.not. all(ieee_is_finite(xyz)) .or. .not. all(ieee_is_finite(q)) .or. &
         & .not. all(ieee_is_finite([rho, beta, dielectric, smear]))) then
         call fatal_error(error, "DRISM requires finite solvent data and dielectric settings")
         return
      end if
      if (rho <= 0.0_wp .or. beta <= 0.0_wp .or. dielectric < 1.0_wp .or. smear < 0.0_wp) then
         call fatal_error(error, "DRISM needs positive density and beta, dielectric >= 1 and smear >= 0")
         return
      end if
      call check_bulk_charge_neutrality(q, error)
      if (allocated(error)) return
      qabs = sum(abs(q))
      if (qabs <= 0.0_wp) then
         call fatal_error(error, "DRISM requires a polar solvent with a nonzero site dipole")
         return
      end if
      centre = matmul(xyz, abs(q))/qabs
      centred = xyz - spread(centre, 2, size(q))
      dipole = matmul(centred, q)
      mu = norm2(dipole)
      ! RISMiCal's dipole floor, e bohr; excludes nonpolar models
      if (mu < 1.0e-4_wp) then
         call fatal_error(error, "DRISM requires a polar solvent with a nonzero site dipole")
         return
      end if
      ez = dipole/mu
      ex = 0.0_wp
      best_radius2 = 0.0_wp
      ! Complete the frame with the longest transverse site displacement
      ! A relative tie tolerance keeps symmetric sites stable under rotation;
      ! choosing the other water hydrogen only reverses the transverse axes
      do i = 1, size(q)
         transverse = centred(:, i) - dot_product(centred(:, i), ez)*ez
         radius2 = sum(transverse**2)
         if (radius2 > best_radius2*(1.0_wp + 64.0_wp*epsilon(1.0_wp))) then
            best_radius2 = radius2
            ex = transverse
         end if
      end do
      allocate (self%coord(3, size(q)))
      if (best_radius2 <= epsilon(1.0_wp)**2*max(1.0_wp, maxval(abs(centred))**2)) then
         ! A collinear molecule has no transverse factors
         self%coord(1:2, :) = 0.0_wp
      else
         ex = ex/sqrt(best_radius2)
         ey = [ez(2)*ex(3) - ez(3)*ex(2), ez(3)*ex(1) - ez(1)*ex(3), &
            & ez(1)*ex(2) - ez(2)*ex(1)]
         do i = 1, size(q)
            self%coord(1, i) = dot_product(centred(:, i), ex)
            self%coord(2, i) = dot_product(centred(:, i), ey)
         end do
      end if
      do i = 1, size(q)
         self%coord(3, i) = dot_product(centred(:, i), ez)
      end do
      y = (4.0_wp*acos(-1.0_wp)/9.0_wp)*beta*rho*mu**2
      self%amplitude = ((dielectric - 1.0_wp)/y - 3.0_wp)/rho
      self%smear = smear
      if (.not. ieee_is_finite(self%amplitude)) then
         call fatal_error(error, "DRISM dielectric amplitude is not finite")
      end if
   end subroutine new_drism_theory

   !> Canonical name
   !>
   !> @param[in] self Fixed theory
   pure function drism_name(self) result(name)
      !> Fixed theory
      class(drism_theory_type), intent(in) :: self
      !> Canonical name, "DRISM"
      character(len=8) :: name
      name = "DRISM"
   end function drism_name

   !> Evaluate the correction through F/k; zeta/k^2 is finite at zero
   !>
   !> @param[in] self Fixed theory
   !> @param[in] k Nonnegative wave number, 1/bohr
   !> @param[out] zeta Dielectric correction, bohr^3 (ns, ns)
   !> @param[out] zeta_over_k2 Regularized correction, bohr^5 (ns, ns)
   pure subroutine drism_correction_at(self, k, zeta, zeta_over_k2)
      !> Fixed theory
      class(drism_theory_type), intent(in) :: self
      !> Nonnegative wave number, 1/bohr
      real(wp), intent(in) :: k
      !> Dielectric correction (bohr^3) and regularized correction (bohr^5), both (ns, ns)
      real(wp), intent(out) :: zeta(:, :), zeta_over_k2(:, :)
      real(wp) :: f(size(zeta, 1)), amplitude
      integer :: i, j

      amplitude = self%amplitude*exp(-0.5_wp*(self%smear*k)**2)
      do i = 1, size(f)
         f(i) = j0(k*self%coord(1, i))*j0(k*self%coord(2, i))*self%coord(3, i)* &
            & j1_over_x(k*self%coord(3, i))
      end do
      do j = 1, size(f)
         do i = 1, size(f)
            zeta_over_k2(i, j) = amplitude*f(i)*f(j)
         end do
      end do
      zeta = k*k*zeta_over_k2
   end subroutine drism_correction_at

   !> Spherical j0, with a cancellation-free small-argument series
   !>
   !> @param[in] x Dimensionless argument
   pure elemental function j0(x) result(value)
      !> Dimensionless argument
      real(wp), intent(in) :: x
      !> Value of j0(x)
      real(wp) :: value, x2
      x2 = x*x
      if (x2 < 0.01_wp) then
         value = 1.0_wp + x2*(-1.0_wp/6.0_wp + x2*(1.0_wp/120.0_wp + &
            & x2*(-1.0_wp/5040.0_wp + x2/362880.0_wp)))
      else
         value = sin(x)/x
      end if
   end function j0

   !> Spherical j1(x)/x, including the exact limit 1/3
   !>
   !> @param[in] x Dimensionless argument
   pure elemental function j1_over_x(x) result(value)
      !> Dimensionless argument
      real(wp), intent(in) :: x
      !> Value of j1(x)/x
      real(wp) :: value, x2
      x2 = x*x
      if (x2 < 0.01_wp) then
         value = 1.0_wp/3.0_wp + x2*(-1.0_wp/30.0_wp + x2*(1.0_wp/840.0_wp + &
            & x2*(-1.0_wp/45360.0_wp + x2/3991680.0_wp)))
      else
         value = (sin(x)/x - cos(x))/x2
      end if
   end function j1_over_x
end module moist_model_moz_solvent_theory_drism
