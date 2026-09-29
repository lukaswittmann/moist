!> Suite proving the Cartesian batched FFT takes the threaded ducc0 path
!>
!> - Ordinary suites run inside test-drive's `!$omp parallel` region, where
!>   `omp_in_parallel()` is always true, so `cartesian_trafo_fft_r2k`/`fft_k2r`
!>   force nthreads = 1 and only the single-field fallback runs
!> - Meson/CMake test registration invokes each test of this suite directly
!>   (tester `<suite> <test>`), which test-drive's `run_selected` dispatches
!>   with no enclosing `!$omp parallel` region, so the threaded ducc0 path is
!>   actually exercised
!> - Companion of the `math_grid_3d` suite, whose batched test cannot reach
!>   that path
module test_math_grid_3d_threaded
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_in_parallel
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, skip_test
   use moist_math_fft, only: moist_fft_r2c_3d, moist_fft_c2r_3d
   use moist_math_grid_3d_base, only: moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_grid_3d
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_threaded

contains

   !> Collect all math_grid_3d_threaded tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_3d_threaded(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("batched_fft_uses_threaded_backend", test_threaded_batch), &
                  new_unittest("batched_ifft_uses_threaded_backend", test_threaded_round_trip) &
                  ]
   end subroutine collect_math_grid_3d_threaded

   !> Check that the threaded ducc0 path is reachable, then pin the thread count
   !>
   !> - Skips when OpenMP is disabled or only one thread is available
   !> - Fails when called inside an enclosing OpenMP team, i.e. when the test
   !>   was run by the suite-level parallel loop instead of as a selected test
   !> - Otherwise sets up to four threads; the caller restores `max_threads`
   !>
   !> @param[out] error        test failure, or set by `skip_test`
   !> @param[out] max_threads  thread count to restore afterwards
   subroutine enter_threaded(error, max_threads)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error
      !> Thread count to restore afterwards
      integer, intent(out) :: max_threads

      logical :: has_omp, in_parallel

      has_omp = .false.
      max_threads = 1
      in_parallel = .false.
!$    has_omp = .true.
!$    max_threads = omp_get_max_threads()
!$    in_parallel = omp_in_parallel()

      if (.not. has_omp .or. max_threads < 2) then
         call skip_test(error, "OpenMP disabled or only one thread available")
         return
      end if

      call check(error, .not. in_parallel, &
         & "threaded FFT check must run as a selected test, not inside the suite-level parallel run")
      if (allocated(error)) return

!$    call omp_set_num_threads(min(4, max_threads))
   end subroutine enter_threaded

   !> Outside any enclosing OpenMP team, batched forward FFT matches the reference
   !>
   !> - `cartesian_trafo_fft_r2k` selects nthreads = omp_get_max_threads() and
   !>   routes through the batched ducc0 backend (see its `nthreads == 1` branch)
   !> - Preconditions forcing that path are checked by `enter_threaded` before
   !>   the transform runs
   !> - Batched result is then compared with the single-field backend call
   !>
   !> @param[out] error  test failure, or set by `skip_test`
   subroutine test_threaded_batch(error)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f_r(:, :)
      complex(wp), allocatable :: f_k(:, :), f_k_ref(:, :)
      integer :: max_threads, i, status

      call enter_threaded(error, max_threads)
      if (allocated(error)) return

      call new_cartesian_grid_3d(grid, 15, 14, 13, 0.4_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call grid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      allocate (f_r(grid%ngrid, 1), f_k(grid%npts_k, 1), f_k_ref(grid%npts_k, 1))
      do i = 1, grid%ngrid
         f_r(i, 1) = sin(0.013_wp*i) + cos(0.007_wp*i)
      end do

      status = moist_fft_r2c_3d(int(grid%nz, c_int), int(grid%ny, c_int), int(grid%nx, c_int), &
         & f_r(:, 1), f_k_ref(:, 1), grid%dv)
      call check(error, status == 0, "reference forward FFT backend call failed")
      if (allocated(error)) then
         call trafo%destroy(); call grid%destroy(); return
      end if

      call trafo%fft_r2k(f_r, f_k, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
      else
         call check(error, maxval(abs(f_k(:, 1) - f_k_ref(:, 1))) < 1.0e-11_wp, &
            & "threaded batched forward transform deviates from the single-field backend call")
      end if

!$    call omp_set_num_threads(max_threads)

      deallocate (f_r, f_k, f_k_ref)
      call trafo%destroy()
      call grid%destroy()
   end subroutine test_threaded_batch

   !> Outside any enclosing OpenMP team, a multi-site batch transforms both ways
   !>
   !> - Three sites on an odd/even 15 x 14 x 13 grid, so the threaded
   !>   `moist_fft_r2c_3d_batch` and `moist_fft_c2r_3d_batch` paths both run
   !>   with sites sharing one worker pool
   !> - Forward block matches per-site single-field forward calls
   !> - Backward block matches per-site single-field backward calls, including
   !>   the 1/Vbox scaling; C2R destroys its input, so every call gets its own
   !>   copy of the spectrum
   !> - Backward of forward reproduces the input
   !>
   !> @param[out] error  test failure, or set by `skip_test`
   subroutine test_threaded_round_trip(error)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      integer, parameter :: ns = 3
      real(wp), allocatable :: f_r(:, :), g_r(:, :), g_r_ref(:, :)
      complex(wp), allocatable :: f_k(:, :), f_k_ref(:, :), scratch(:)
      integer :: max_threads, i, a, status

      call enter_threaded(error, max_threads)
      if (allocated(error)) return

      call new_cartesian_grid_3d(grid, 15, 14, 13, 0.4_wp, error=merr)
      if (.not. allocated(merr)) call grid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
!$       call omp_set_num_threads(max_threads)
         return
      end if
      allocate (f_r(grid%ngrid, ns), g_r(grid%ngrid, ns), g_r_ref(grid%ngrid, ns))
      allocate (f_k(grid%npts_k, ns), f_k_ref(grid%npts_k, ns), scratch(grid%npts_k))
      do a = 1, ns
         do i = 1, grid%ngrid
            f_r(i, a) = sin(0.013_wp*i*a) + cos(0.007_wp*i) + 0.1_wp*a
         end do
      end do

      ! Per-site single-field references in both directions
      do a = 1, ns
         status = moist_fft_r2c_3d(int(grid%nz, c_int), int(grid%ny, c_int), int(grid%nx, c_int), &
            & f_r(:, a), f_k_ref(:, a), grid%dv)
         if (status == 0) then
            scratch = f_k_ref(:, a)
            status = moist_fft_c2r_3d(int(grid%nz, c_int), int(grid%ny, c_int), int(grid%nx, c_int), &
               & scratch, g_r_ref(:, a), 1.0_wp/grid%vbox)
         end if
         call check(error, status == 0, "reference single-field FFT backend call failed")
         if (allocated(error)) exit
      end do

      if (.not. allocated(error)) then
         call trafo%fft_r2k(f_r, f_k, merr)
         if (allocated(merr)) call test_failed(error, merr%message)
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(f_k - f_k_ref)) < 1.0e-11_wp, &
            & "threaded multi-site forward transform deviates from the single-field backend calls")
      end if
      if (.not. allocated(error)) then
         call trafo%fft_k2r(f_k, g_r, merr)
         if (allocated(merr)) call test_failed(error, merr%message)
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(g_r - g_r_ref)) < 1.0e-12_wp, &
            & "threaded multi-site backward transform deviates from the single-field backend calls")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(g_r - f_r)) < 1.0e-12_wp, &
            & "threaded multi-site backward transform does not invert the forward one")
      end if

!$    call omp_set_num_threads(max_threads)

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_threaded_round_trip

end module test_math_grid_3d_threaded
