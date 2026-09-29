!> Chebyshev-2 radial quadrature for atom-centered integration grids
module moist_math_quadrature_chebyshev
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   implicit none(type, external)
   private

   public :: chebyshev2_radii

contains

   !> Generate nr Chebyshev-2 radial nodes and weights
   !>
   !> The Jacobian `r^2 dr` is folded into the weights
   !>
   !> @param[in]  nr      Number of radial points (nr >= 1)
   !> @param[in]  p       Radial scale (bohr); typically a fraction of
   !>                     the atomic covalent radius
   !> @param[out] radii   Radii in bohr, shape (nr)
   !> @param[out] weights Weights in bohr^3, shape (nr); include r^2 dr
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

end module moist_math_quadrature_chebyshev
