!> Concurrent scalar rendering with moist_utils_format
!>
!> Run as a selected test, outside test-drive's OpenMP team
module test_utils_format_threaded
   use, intrinsic :: iso_fortran_env, only: int32, real64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_utils_format_threaded

contains

   !> Collect concurrent rendering tests
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_utils_format_threaded(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("concurrent_rendering", test_concurrent_rendering)]
   end subroutine collect_utils_format_threaded

   !> Threads rendering numbers of varying length get the serial results
   !>
   !> - a deferred-length result races here on gfortran (GCC PR 113797); the
   !>   race hits only a fraction of calls, so a pass does not prove its absence
   !> - integers take the I/O-free path; reals go through serialized internal writes
   !>
   subroutine test_concurrent_rendering(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nval = 8, nint_calls = 400000, nreal_calls = 20000
      integer(int32), parameter :: ints(nval) = [0, -7, 42, 1000, -123456, 99999999, &
         & huge(0_int32), -huge(0_int32) - 1]
      real(real64), parameter :: reals(nval) = [0.0_real64, -1.5_real64, 1.0e-3_real64, &
         & 2.5e10_real64, -7.0e-300_real64, 3.0_real64, 1.0e100_real64, -42.0_real64]
      character(len=32) :: int_ref(nval), real_ref(nval)
      integer :: i, k, failures

      do k = 1, nval
         int_ref(k) = to_string(ints(k))
         real_ref(k) = to_string(reals(k), "(es12.4)")
      end do

      failures = 0
      !$omp parallel do num_threads(4) schedule(static, 1) private(k) reduction(+:failures)
      do i = 1, nint_calls
         k = modulo(i, nval) + 1
         if (.not. same_text(to_string(ints(k)), int_ref(k))) failures = failures + 1
      end do
      !$omp end parallel do
      call check(error, failures, 0, "concurrent integer renderings differ from the serial ones")
      if (allocated(error)) return

      !$omp parallel do num_threads(4) schedule(static, 1) private(k) reduction(+:failures)
      do i = 1, nreal_calls
         k = modulo(i, nval) + 1
         if (.not. same_text(to_string(reals(k), "(es12.4)"), real_ref(k))) failures = failures + 1
      end do
      !$omp end parallel do
      call check(error, failures, 0, "concurrent real renderings differ from the serial ones")
   end subroutine test_concurrent_rendering

   !> Whether a rendering equals its reference, length included
   !>
   !> @param[in] text Rendering
   !> @param[in] ref  Reference, blank-padded
   pure function same_text(text, ref) result(same)
      !> Rendering
      character(len=*), intent(in) :: text
      !> Reference, blank-padded
      character(len=*), intent(in) :: ref
      !> Whether they agree
      logical :: same

      same = len(text) == len_trim(ref) .and. text == ref
   end function same_text

end module test_utils_format_threaded
