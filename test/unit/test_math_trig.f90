!> Test suite for trigonometric utilities in moist_math_trigonometry
module test_math_trig
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_trigonometry, only: rotation_z_to_n
   implicit none (type, external)
   private

   public :: collect_math_trig

   real(wp), parameter :: thr = 100.0_wp * epsilon(1.0_wp)

contains

   !> Collect all trigonometry tests
   subroutine collect_math_trig(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("rot_z_to_z", test_rot_z_to_z), &
         new_unittest("rot_z_to_neg_z", test_rot_z_to_neg_z), &
         new_unittest("rot_z_to_x", test_rot_z_to_x), &
         new_unittest("rot_z_to_y", test_rot_z_to_y), &
         new_unittest("rot_z_to_arbitrary", test_rot_z_to_arbitrary), &
         new_unittest("rot_orthogonality", test_rot_orthogonality), &
         new_unittest("rot_determinant", test_rot_determinant) &
      ]
   end subroutine collect_math_trig

   !> Identity case: rotating z to z should give the identity matrix
   subroutine test_rot_z_to_z(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), ez(3)
      integer :: i, j

      ez = [0.0_wp, 0.0_wp, 1.0_wp]
      call rotation_z_to_n(ez, R)

      do i = 1, 3
         do j = 1, 3
            if (i == j) then
               call check(error, abs(R(i,j) - 1.0_wp) < thr, "Diagonal should be 1")
            else
               call check(error, abs(R(i,j)) < thr, "Off-diagonal should be 0")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_rot_z_to_z

   !> Anti-parallel case: rotating z to -z
   subroutine test_rot_z_to_neg_z(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), neg_z(3), result(3)

      neg_z = [0.0_wp, 0.0_wp, -1.0_wp]
      call rotation_z_to_n(neg_z, R)

      ! R * [0,0,1] should give [0,0,-1]
      result = matmul(R, [0.0_wp, 0.0_wp, 1.0_wp])
      call check(error, abs(result(1)) < thr, "x component should be 0")
      if (allocated(error)) return
      call check(error, abs(result(2)) < thr, "y component should be 0")
      if (allocated(error)) return
      call check(error, abs(result(3) + 1.0_wp) < thr, "z component should be -1")
   end subroutine test_rot_z_to_neg_z

   !> Rotate z-axis to x-axis
   subroutine test_rot_z_to_x(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), ex(3), result(3)

      ex = [1.0_wp, 0.0_wp, 0.0_wp]
      call rotation_z_to_n(ex, R)

      result = matmul(R, [0.0_wp, 0.0_wp, 1.0_wp])
      call check(error, abs(result(1) - 1.0_wp) < thr, "x component should be 1")
      if (allocated(error)) return
      call check(error, abs(result(2)) < thr, "y component should be 0")
      if (allocated(error)) return
      call check(error, abs(result(3)) < thr, "z component should be 0")
   end subroutine test_rot_z_to_x

   !> Rotate z-axis to y-axis
   subroutine test_rot_z_to_y(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), ey(3), result(3)

      ey = [0.0_wp, 1.0_wp, 0.0_wp]
      call rotation_z_to_n(ey, R)

      result = matmul(R, [0.0_wp, 0.0_wp, 1.0_wp])
      call check(error, abs(result(1)) < thr, "x component should be 0")
      if (allocated(error)) return
      call check(error, abs(result(2) - 1.0_wp) < thr, "y component should be 1")
      if (allocated(error)) return
      call check(error, abs(result(3)) < thr, "z component should be 0")
   end subroutine test_rot_z_to_y

   !> Rotate z-axis to an arbitrary normalised direction
   subroutine test_rot_z_to_arbitrary(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), n(3), result(3), nrm

      n = [1.0_wp, 2.0_wp, 3.0_wp]
      nrm = norm2(n)
      n = n / nrm

      call rotation_z_to_n(n, R)

      result = matmul(R, [0.0_wp, 0.0_wp, 1.0_wp])
      call check(error, abs(result(1) - n(1)) < thr, "x component mismatch")
      if (allocated(error)) return
      call check(error, abs(result(2) - n(2)) < thr, "y component mismatch")
      if (allocated(error)) return
      call check(error, abs(result(3) - n(3)) < thr, "z component mismatch")
   end subroutine test_rot_z_to_arbitrary

   !> Verify R is orthogonal: R^T R = I
   subroutine test_rot_orthogonality(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), RtR(3, 3), n(3)
      integer :: i, j

      n = [0.3_wp, -0.7_wp, 0.5_wp]
      n = n / norm2(n)
      call rotation_z_to_n(n, R)

      RtR = matmul(transpose(R), R)
      do i = 1, 3
         do j = 1, 3
            if (i == j) then
               call check(error, abs(RtR(i,j) - 1.0_wp) < thr, "R^T R diagonal should be 1")
            else
               call check(error, abs(RtR(i,j)) < thr, "R^T R off-diagonal should be 0")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_rot_orthogonality

   !> Verify det(R) = +1 (proper rotation, not reflection)
   subroutine test_rot_determinant(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: R(3, 3), n(3), det

      n = [-0.6_wp, 0.1_wp, -0.8_wp]
      n = n / norm2(n)
      call rotation_z_to_n(n, R)

      det = R(1,1)*(R(2,2)*R(3,3) - R(2,3)*R(3,2)) &
          - R(1,2)*(R(2,1)*R(3,3) - R(2,3)*R(3,1)) &
          + R(1,3)*(R(2,1)*R(3,2) - R(2,2)*R(3,1))

      call check(error, abs(det - 1.0_wp) < thr, "Determinant should be +1")
   end subroutine test_rot_determinant

end module test_math_trig
