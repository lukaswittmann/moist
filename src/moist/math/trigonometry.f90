!> Trigonometric and rotational utilities for 3D geometry
!>
!> Provides rotation matrix construction for aligning coordinate frames
!> to arbitrary directions, used e.g. for orienting Lebedev grids along
!> surface normals.
module moist_math_trigonometry
   use mctc_env_accuracy, only: wp
   use moist_math_linalg, only: cross_product
   implicit none
   private

   public :: rotation_z_to_n

contains

   !> Build rotation matrix that maps the z-axis [0,0,1] to a given unit vector.
   !>
   !> Uses Rodrigues' formula. Falls back to trivial matrices when the
   !> target is (anti-)parallel to z, where the cross-product axis is
   !> undefined.
   !>
   !> @param[in]  n_hat  Target unit vector (must be normalised)
   !> @param[out] R      3x3 rotation matrix satisfying R * [0,0,1] = n_hat
   pure subroutine rotation_z_to_n(n_hat, R)
      !> Target unit vector
      real(wp), intent(in) :: n_hat(3)
      !> Resulting 3x3 rotation matrix
      real(wp), intent(out) :: R(3, 3)

      real(wp) :: ez(3), v(3), c, s, t
      real(wp), parameter :: tol = 1.0e-14_wp

      ez = [0.0_wp, 0.0_wp, 1.0_wp]
      c = dot_product(ez, n_hat)

      if (c > 1.0_wp - tol) then
         ! n_hat ~ +z, identity
         R = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
         return
      end if

      if (c < -1.0_wp + tol) then
         ! n_hat ~ -z, rotate 180 degrees around x
         R = reshape([1, 0, 0, 0, -1, 0, 0, 0, -1], [3, 3])
         return
      end if

      ! v = ez x n_hat  (rotation axis)
      v = cross_product(ez, n_hat)
      s = norm2(v)
      t = (1.0_wp - c)/(s*s)

      ! Rodrigues: R = I + [v]_x + [v]_x^2 * (1-c)/s^2
      R(1, 1) = 1.0_wp + t*(-v(3)**2 - v(2)**2)
      R(2, 1) = v(3) + t*(v(1)*v(2))
      R(3, 1) = -v(2) + t*(v(1)*v(3))
      R(1, 2) = -v(3) + t*(v(1)*v(2))
      R(2, 2) = 1.0_wp + t*(-v(3)**2 - v(1)**2)
      R(3, 2) = v(1) + t*(v(2)*v(3))
      R(1, 3) = v(2) + t*(v(1)*v(3))
      R(2, 3) = -v(1) + t*(v(2)*v(3))
      R(3, 3) = 1.0_wp + t*(-v(2)**2 - v(1)**2)
   end subroutine rotation_z_to_n

end module moist_math_trigonometry
