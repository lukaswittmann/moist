!> Chebyshev-2 radial quadrature for atom-centered integration grids.
!>
!> Maps the Chebyshev-Gauss nodes `x_i = cos(i*pi/(n+1))` on [-1, 1] to
!> radii on [0, inf) via `r = p * (1 + x) / (1 - x)`, where `p` is an
!> atom-dependent scale (typically a fraction of the covalent radius,
!> in bohr). This is Becke's 1988 second-kind Chebyshev quadrature.
!>
!> Weights returned by `chebyshev2_radii` include the `r^2 dr` Jacobian
!> so that, for a spherically symmetric f,
!>    integral_0^inf f(r) r^2 dr  approximately  sum_i w_i f(r_i)
!>
!> Units convention: `p` is supplied in bohr; returned radii and weights
!> are in bohr (and bohr^3 respectively). The caller is responsible for
!> ensuring the scale matches the rest of the atomic units workflow.
module moist_math_grid_radial
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   implicit none
   private

   public :: chebyshev2_radii

contains

   !> Generate nr Chebyshev-2 radial nodes and weights with `r^2 dr`
   !> Jacobian folded into the weights.
   !>
   !> @param[in]  nr      Number of radial points (nr >= 1).
   !> @param[in]  p       Radial scale (bohr); typically a fraction of
   !>                     the atomic covalent radius.
   !> @param[out] radii   Radii in bohr, shape (nr).
   !> @param[out] weights Weights in bohr^3, shape (nr); include r^2 dr.
   pure subroutine chebyshev2_radii(nr, p, radii, weights)
      !> Number of radial points
      integer, intent(in)  :: nr
      !> Radial scale in bohr
      real(wp), intent(in)  :: p
      !> Radial nodes in bohr
      real(wp), intent(out) :: radii(:)
      !> Radial weights in bohr^3 (carrying the r^2 dr Jacobian)
      real(wp), intent(out) :: weights(:)

      integer  :: ir
      real(wp) :: x_i, one_minus_x

      do ir = 1, nr
         x_i = cos(real(ir, wp)*pi/real(nr + 1, wp))
         one_minus_x = 1.0_wp - x_i
         radii(ir) = (1.0_wp + x_i)/one_minus_x*p
         ! Chebyshev-2 weight (2*pi/(nr+1)) times Jacobian dr/dx = 2*p/(1-x)^2
         ! combined with r^2 and simplified:
         !   w_i = (2*pi/(nr+1)) * p^3 * (1+x)^2.5 / (1-x)^3.5
         weights(ir) = (2.0_wp*pi/real(nr + 1, wp)) &
                      & *p**3*(1.0_wp + x_i)**2.5_wp/one_minus_x**3.5_wp
      end do
   end subroutine chebyshev2_radii

end module moist_math_grid_radial
