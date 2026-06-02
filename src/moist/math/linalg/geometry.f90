!> Geometric helpers and numerically stable scalar utilities.
!>
!> Tangent-frame construction on a sphere and the softplus (smooth maximum)
!> function used by the smooth cavity construction.
module moist_math_linalg_geometry
   use mctc_env_accuracy, only: wp
   implicit none
   private

   public :: cross_product
   public :: setup_tangent_frame
   public :: logaddexp

contains

   !> Right-handed cross product of two 3-vectors, c = a x b.
   !>
   !> Direct 6-multiply/3-subtract implementation. At this fixed size a hand
   !> written kernel beats any BLAS call (no general-n analogue exists anyway)
   !> and stays inlinable by the compiler.
   !> @param[in] a  First vector [3]
   !> @param[in] b  Second vector [3]
   !> @returns   c  Cross product a x b [3]
   pure function cross_product(a, b) result(c)
      !> First vector
      real(wp), intent(in) :: a(3)
      !> Second vector
      real(wp), intent(in) :: b(3)
      !> Cross product a x b
      real(wp) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross_product

   !> Construct orthonormal tangent frame on a sphere
   !>
   !> Given a normal vector n (e.g., radial direction on sphere), constructs
   !> two orthonormal tangent vectors t1, t2 such that {t1, t2, n} form a
   !> right-handed orthonormal basis. Uses Gram-Schmidt orthogonalization
   !> starting from the coordinate axis least aligned with n.
   !>
   !> @param[in]  normal Normal vector (will be normalized internally)
   !> @param[out] t1     First tangent vector (orthogonal to normal)
   !> @param[out] t2     Second tangent vector (orthogonal to normal and t1)
   pure subroutine setup_tangent_frame(normal, t1, t2)
      !> Normal vector (input, will be normalized)
      real(wp), intent(in) :: normal(3)
      !> First tangent vector (output)
      real(wp), intent(out) :: t1(3)
      !> Second tangent vector (output)
      real(wp), intent(out) :: t2(3)

      real(wp) :: n_normalized(3), n_norm
      real(wp) :: v_tmp(3), v_norm, proj
      integer :: min_axis

      ! Normalize the input normal vector
      n_norm = sqrt(dot_product(normal, normal))
      n_normalized = normal/n_norm

      ! Construct t1 via Gram-Schmidt: pick coordinate axis least aligned with normal
      min_axis = minloc(abs(n_normalized), dim=1)
      t1 = 0.0_wp
      t1(min_axis) = 1.0_wp

      ! Project out normal component and normalize
      proj = dot_product(t1, n_normalized)
      v_tmp = t1 - proj*n_normalized
      v_norm = sqrt(dot_product(v_tmp, v_tmp))
      t1 = v_tmp/v_norm

      ! Construct t2 as cross product: t2 = n x t1
      t2 = cross_product(n_normalized, t1)

   end subroutine setup_tangent_frame

   !> Numerically stable computation of the softplus function
   !> Computes the smooth maximum (softplus) function
   !>
   !> f(x) = \log(1 + e^x)
   !>
   !> Uses numerically stable formulation to avoid overflow for large x
   !> and underflow for large negative x:
   !> - For x > 0: compute as x + log(1 + exp(-x))
   !> - For x <= 0: compute directly as log(1 + exp(x))
   !>
   !> @param[in] x   Input value
   !> @returns   val Softplus value
   pure elemental function logaddexp(x) result(val)
      !> Input value
      real(wp), intent(in) :: x
      !> Result: log(1 + exp(x))
      real(wp) :: val

      ! Numerically stable formulation
      if (x > 0.0_wp) then
         val = x + log(1.0_wp + exp(-x))
      else
         val = log(1.0_wp + exp(x))
      end if

   end function logaddexp

end module moist_math_linalg_geometry
