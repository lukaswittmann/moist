!> Ng split of the Coulomb interaction
!>
!> - Unit-charge real tail: f(r) = erf(alpha r)/r
!> - Boys form: f(r) = (2 alpha/sqrt(pi)) F0(alpha^2 r^2)
!> - Origin limit: f(0) = 2 alpha/sqrt(pi)
!> - Reciprocal envelope: E(k) = 4 pi exp(-k^2/(4 alpha^2))/k^2
!> - Dropped k = 0 mode: E(0) = 0
!> - Radial gradient factor: g(r) = f'(r)/r = -(4 alpha^3/sqrt(pi)) F1(alpha^2 r^2)
!> - Separation-vector multiplier in gradients: g(r)
!> - F0, F1 from `moist_math_boys: boys01`; accurate near r = 0 without cancellation
!> - `require_split_alpha` validation at model setter or solvent tables
!> - Reference: K.-C. Ng, J. Chem. Phys. 61, 2680 (1974); doi:10.1063/1.1682399
module moist_model_moz_potential_kernel_ng
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_boys, only: boys01
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none(type, external)
   private

   public :: ng_real_tail, ng_fourier_envelope, require_split_alpha

   !> 2/sqrt(pi)
   real(wp), parameter :: two_over_sqrt_pi = 2.0_wp/sqrt(acos(-1.0_wp))
   !> 4 pi
   real(wp), parameter :: four_pi = 4.0_wp*acos(-1.0_wp)
   !> Squared wavenumber treated as the k = 0 mode, 1/bohr^2
   real(wp), parameter :: ng_k2_zero = 1.0e-20_wp

contains

   !> Evaluate the real-space Ng tail of unit charges and optional radial factor
   !>
   !> @param[in]  r      distance, bohr
   !> @param[in]  alpha  split parameter, 1/bohr
   !> @param[out] f      erf(alpha r)/r, 1/bohr
   !> @param[out] g      optional f'(r)/r, 1/bohr^3
   elemental subroutine ng_real_tail(r, alpha, f, g)
      !> Distance, bohr
      real(wp), intent(in) :: r
      !> Split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> erf(alpha r)/r
      real(wp), intent(out) :: f
      !> f'(r)/r
      real(wp), intent(out), optional :: g
      real(wp) :: x, f0, f1, ex
      x = alpha*r
      call boys01(x*x, f0, f1, ex)
      f = two_over_sqrt_pi*alpha*f0
      if (present(g)) g = -2.0_wp*two_over_sqrt_pi*alpha**3*f1
   end subroutine ng_real_tail

   !> Reciprocal Ng envelope of unit charges, zero for the k = 0 mode
   !>
   !> @param[in] k2     squared wavenumber, 1/bohr^2
   !> @param[in] alpha  split parameter, 1/bohr
   elemental function ng_fourier_envelope(k2, alpha) result(env)
      !> Squared wavenumber, 1/bohr^2
      real(wp), intent(in) :: k2
      !> Split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> 4 pi exp(-k^2/(4 alpha^2))/k^2, bohr^2
      real(wp) :: env
      if (k2 < ng_k2_zero) then
         env = 0.0_wp
      else
         env = four_pi*exp(-0.25_wp*k2/(alpha*alpha))/k2
      end if
   end function ng_fourier_envelope

   !> Refuse a non-positive or non-finite split parameter when the split is on
   !>
   !> @param[in]  ng_split  whether the Coulomb tail is split off
   !> @param[in]  alpha     ng split parameter, 1/bohr
   !> @param[out] error     invalid split parameter
   subroutine require_split_alpha(ng_split, alpha, error)
      !> Whether the Coulomb tail is split off
      logical, intent(in) :: ng_split
      !> Ng split parameter, 1/bohr
      real(wp), intent(in) :: alpha
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. ng_split) return
      if (.not. ieee_is_finite(alpha)) then
         call fatal_error(error, "Ng split parameter alpha must be finite")
      else if (alpha <= 0.0_wp) then
         call fatal_error(error, "Ng split parameter alpha must be positive")
      end if
   end subroutine require_split_alpha

end module moist_model_moz_potential_kernel_ng
