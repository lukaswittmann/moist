!> Test suite for the scalar rendering of moist_utils_format
module test_utils_format
   use, intrinsic :: iso_fortran_env, only: int32, int64, real64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_utils_format

contains

   !> Collect all scalar rendering tests
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_utils_format(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("int_as_i0", test_int_as_i0), &
                  new_unittest("int64_as_i0", test_int64_as_i0), &
                  new_unittest("real_as_g0", test_real_as_g0), &
                  new_unittest("real_with_format", test_real_with_format), &
                  new_unittest("failed_write", test_failed_write) &
                  ]

   end subroutine collect_utils_format

   !> Default integers render as `(i0)` does, digit boundaries and extremes included
   !>
   !> @param[out] error Test failure
   subroutine test_int_as_i0(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer(int32), parameter :: vals(14) = [0_int32, 7_int32, -7_int32, 9_int32, 10_int32, &
         & -10_int32, 99_int32, 100_int32, -100_int32, 12345_int32, -99999_int32, 1000000000_int32, &
         & huge(0_int32), -huge(0_int32) - 1_int32]
      character(len=32) :: ref
      integer :: k

      do k = 1, size(vals)
         write (ref, "(i0)") vals(k)
         call check(error, to_string(vals(k)), trim(ref), "int32 differs from (i0)")
         if (allocated(error)) return
         call check(error, len(to_string(vals(k))), len_trim(ref), "int32 rendering has blanks")
         if (allocated(error)) return
      end do
   end subroutine test_int_as_i0

   !> 64-bit integers render as `(i0)` does, both extremes included
   !>
   !> @param[out] error Test failure
   subroutine test_int64_as_i0(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer(int64), parameter :: vals(7) = [0_int64, -1_int64, 10_int64**18, -10_int64**18, &
         & 9223372036854775806_int64, huge(0_int64), -huge(0_int64) - 1_int64]
      character(len=32) :: ref
      integer :: k

      do k = 1, size(vals)
         write (ref, "(i0)") vals(k)
         call check(error, to_string(vals(k)), trim(ref), "int64 differs from (i0)")
         if (allocated(error)) return
         call check(error, len(to_string(vals(k))), len_trim(ref), "int64 rendering has blanks")
         if (allocated(error)) return
      end do
   end subroutine test_int64_as_i0

   !> Reals without a format render as `(g0)` does, blanks stripped
   !>
   !> @param[out] error Test failure
   subroutine test_real_as_g0(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(real64), parameter :: vals(7) = [0.0_real64, 1.5_real64, -2.5e-300_real64, 1.0e100_real64, &
         & tiny(0.0_real64), huge(0.0_real64), -1.0_real64/3.0_real64]
      character(len=128) :: ref
      integer :: k

      do k = 1, size(vals)
         write (ref, "(g0)") vals(k)
         ref = adjustl(ref)
         call check(error, to_string(vals(k)), trim(ref), "real differs from (g0)")
         if (allocated(error)) return
         call check(error, len(to_string(vals(k))), len_trim(ref), "real rendering has blanks")
         if (allocated(error)) return
      end do
   end subroutine test_real_as_g0

   !> Reals with a format: the edit as written, leading blanks stripped
   !>
   !> @param[out] error Test failure
   subroutine test_real_with_format(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check(error, to_string(1.0e-3_real64, "(es9.2)"), "1.00E-03", "es9.2")
      if (allocated(error)) return
      call check(error, to_string(-3.0_real64, "(f6.2)"), "-3.00", "f6.2 of a negative value")
      if (allocated(error)) return
      call check(error, len(to_string(-3.0_real64, "(f6.2)")), 5, "f6.2 keeps no leading blank")
      if (allocated(error)) return
      call check(error, to_string(12.5_real64, "(f0.1)"), "12.5", "f0.1")
      if (allocated(error)) return
      call check(error, to_string(2.0_real64, "('x = ', f3.1)"), "x = 2.0", "literal text in the format")
   end subroutine test_real_with_format

   !> A failed internal write renders as `"*"`
   !>
   !> @param[out] error Test failure
   subroutine test_failed_write(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check(error, to_string(1.0_real64, "(i5)"), "*", "integer edit for a real")
      if (allocated(error)) return
      call check(error, to_string(1.0_real64, "(f6.2"), "*", "unbalanced parentheses")
   end subroutine test_failed_write

end module test_utils_format
