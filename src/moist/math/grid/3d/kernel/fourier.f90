!> Geometry kernels of quadrature-weighted Fourier sums at a fixed k lattice
module moist_math_grid_3d_kernel_fourier
   use mctc_env, only: wp
   implicit none(type, external)
   private
   public :: fourier_derivative_mode, fourier_measure_adjoint, fourier_position_adjoint

contains

   !> Multiplying a mode by i*k differentiates exp(+i*k*(r-kref)).
   pure elemental function fourier_derivative_mode(k, value) result(derivative)
      !> Wavevector component, bohr^-1
      real(wp), intent(in) :: k
      !> Reciprocal coefficient
      complex(wp), intent(in) :: value
      !> Spatial derivative coefficient
      complex(wp) :: derivative
      derivative = cmplx(0.0_wp, k, wp)*value
   end function fourier_derivative_mode

   !> Forward measure derivative: f_j*Re(sum_k seed_k*exp(+i*k*x_j)).
   pure elemental function fourier_measure_adjoint(field, backward_seed) result(weight)
      !> Primal real field, held fixed
      real(wp), intent(in) :: field
      !> Unnormalized backward sum of the complex output seed
      complex(wp), intent(in) :: backward_seed
      !> Sensitivity to the quadrature measure
      real(wp) :: weight
      weight = field*real(backward_seed, wp)
   end function fourier_measure_adjoint

   !> Forward multiplier is w_j*f_j; backward multiplier is norm*seed_j.
   !> The reference receives the negative sum of these position sensitivities.
   pure elemental function fourier_position_adjoint(multiplier, backward_derivative) result(weight)
      !> Real multiplier for this point
      real(wp), intent(in) :: multiplier
      !> Unnormalized backward sum of i*k times the reciprocal seed/primal
      complex(wp), intent(in) :: backward_derivative
      !> Sensitivity to the point coordinate
      real(wp) :: weight
      weight = multiplier*real(backward_derivative, wp)
   end function fourier_position_adjoint

end module moist_math_grid_3d_kernel_fourier
