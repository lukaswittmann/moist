!> Test suite for releasing the FFT backend's worker threads
!>
!> - `moist_fft_release_threads` succeeds, a repeated call too, and the
!>   backend keeps transforming afterwards with an unchanged result
!> - One test only: the release must not run beside another transform, and
!>   test-drive runs the tests of a suite concurrently
module test_math_fft_threads
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_fft, only: moist_fft_dst4, moist_fft_release_threads
   implicit none(type, external)
   private

   public :: collect_math_fft_threads

contains

   !> Collect all math_fft_threads tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_fft_threads(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("release_keeps_backend_working", test_release_keeps_backend)]

   end subroutine collect_math_fft_threads

   !> Release twice between two identical transforms
   !>
   !> @param[out] error  Test failure
   subroutine test_release_keeps_backend(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Transform length and batch size
      integer(c_int), parameter :: npts = 16, nbatch = 3
      !> Input columns, transforms before and after the release
      real(c_double) :: input(npts*nbatch), before(npts*nbatch), after(npts*nbatch)
      !> Value index
      integer :: i

      do i = 1, size(input)
         input(i) = sin(0.37_wp*real(i, wp))
      end do
      call check(error, moist_fft_dst4(npts, nbatch, input, before, 1.0_c_double) == 0, &
         & more="transform before the release failed")
      if (allocated(error)) return
      call check(error, moist_fft_release_threads() == 0, more="release failed")
      if (allocated(error)) return
      call check(error, moist_fft_release_threads() == 0, more="repeated release failed")
      if (allocated(error)) return
      call check(error, moist_fft_dst4(npts, nbatch, input, after, 1.0_c_double) == 0, &
         & more="transform after the release failed")
      if (allocated(error)) return
      call check(error, all(after == before), more="release changed the transform")

   end subroutine test_release_keeps_backend

end module test_math_fft_threads
