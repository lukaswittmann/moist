!> Suite proving the Cartesian batched FFT takes the threaded ducc0 path
!>
!> - Ordinary suites run inside test-drive's `!$omp parallel` region, which
!>   forces nthreads = 1 in the Cartesian FFT
!> - Each test here is invoked directly (`<suite> <test>`), outside any team
!> - Companion of `math_grid_3d` (batched FFT) and `math_grid_nufft` (NUFFT)
!> - NUFFT: type 1/2 and type 3, both directions, one column and a batch,
!>   against the analytic Gaussian FT and a serial run
module test_math_grid_3d_threaded
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_in_parallel, &
!$    & omp_get_max_active_levels
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, skip_test
   use moist_math_fft, only: moist_fft_r2c_3d, moist_fft_c2r_3d
   use moist_math_grid_3d_base, only: moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_grid_3d
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_grid_uniform, molecular_grid_set_kgrid, &
      & moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo, default_nufft_tol
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_threaded

   !> Radial shells per atom of the threaded NUFFT grid
   !>
   !> Grid sized so FINUFFT takes its multithreaded bin sort
   !>
   !> - FINUFFT sorts multithreaded only when `10*M` exceeds the fine grid `N`
   !> - Water at (nrad, nang, rmax, dr) = (40, 194, 6, 0.7): M = 18932 points,
   !>   type-1 fine grid 60 x 48 x 48 = 138240, so 10*M = 189320 > N
   !> - The `math_grid_nufft` probe grid (50, 110, 16, 0.7) sorts
   !>   single-threaded
   integer, parameter :: nufft_nrad = 40
   !> Raw Lebedev point count per shell of the threaded NUFFT grid
   integer, parameter :: nufft_nang = 194
   !> Outer radial clamp of the threaded NUFFT grid (bohr)
   real(wp), parameter :: nufft_rmax = 6.0_wp
   !> Reciprocal resolution of the threaded NUFFT grid (bohr); k_max = pi/dr
   real(wp), parameter :: nufft_dr = 0.7_wp
   !> Batch widths per route: one column and a batch
   !>
   !> Different FINUFFT paths: a single column spreads and interpolates with
   !> the full team; a batch runs one column per thread, each spreading alone
   integer, parameter :: nufft_widths(2) = [1, 3]
   !> Acceptance bound for the forward transform against the analytic FT
   !>
   !> Quadrature-limited, relative to each column's peak `(pi/a)^(3/2)`;
   !> measured 1.4e-5 (one column) and 6.9e-5 (three columns) on both routes
   real(wp), parameter :: nufft_analytic_tol = 5.0e-4_wp

contains

   !> Collect all math_grid_3d_threaded tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_3d_threaded(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("batched_fft_uses_threaded_backend", test_threaded_batch), &
                  new_unittest("batched_ifft_uses_threaded_backend", test_threaded_round_trip), &
                  new_unittest("molecular_nufft_type12_threaded", test_nufft_type12_threaded), &
                  new_unittest("molecular_nufft_type3_threaded", test_nufft_type3_threaded) &
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
   !> - Batched result is then compared with the single-field backend call, or set by `skip_test`
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
      if (.not. allocated(merr)) call grid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
!$       call omp_set_num_threads(max_threads)
         return
      end if
      allocate (f_r(grid%ngrid, 1), f_k(grid%npts_k, 1), f_k_ref(grid%npts_k, 1))
      do i = 1, grid%ngrid
         f_r(i, 1) = sin(0.013_wp*i) + cos(0.007_wp*i)
      end do

      status = moist_fft_r2c_3d(int(grid%nz, c_int), int(grid%ny, c_int), int(grid%nx, c_int), &
         & f_r(:, 1), f_k_ref(:, 1), grid%dv)
      call check(error, status == 0, "reference forward FFT backend call failed")
      if (allocated(error)) then
!$       call omp_set_num_threads(max_threads)
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
   !> - Backward of forward reproduces the input, or set by `skip_test`
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

   !> Threaded molecular NUFFT on the type-1/2 route, or set by `skip_test`
   subroutine test_nufft_type12_threaded(error)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error

      call run_threaded_nufft(error, .true.)
   end subroutine test_nufft_type12_threaded

   !> Threaded molecular NUFFT on the type-3 route, or set by `skip_test`
   subroutine test_nufft_type3_threaded(error)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error

      call run_threaded_nufft(error, .false.)
   end subroutine test_nufft_type3_threaded

   !> Outside any enclosing OpenMP team, the molecular NUFFT runs multithreaded
   !>
   !> - Preconditions checked by `enter_threaded`, plus a nonzero OpenMP
   !>   max-active-levels: FINUFFT opens its own parallel regions, which an
   !>   `OMP_MAX_ACTIVE_LEVELS=0` environment would serialise
   !> - The trafo plans with `omp_get_max_threads()` workers, so
   !>   `omp_set_num_threads` sizes FINUFFT's team
   !> - Thread count and grid are restored on every exit path
   !>
   !> @param[out] error       test failure, or set by `skip_test`
   !> @param[in]  use_type12  transform route, `.false.` for type 3
   subroutine run_threaded_nufft(error, use_type12)
      !> Test failure, or set by `skip_test`
      type(error_type), allocatable, intent(out) :: error
      !> Transform route, `.false.` for type 3
      logical, intent(in) :: use_type12

      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(mctc_error), allocatable :: merr
      integer :: max_threads, max_levels, k

      call enter_threaded(error, max_threads)
      if (allocated(error)) return

      max_levels = 1
!$    max_levels = omp_get_max_active_levels()
      if (max_levels < 1) then
         call skip_test(error, "OpenMP max active levels is 0, FINUFFT would run serially")
      else
         call new(mol, [8, 1, 1], reshape([ &
                                          0.0_wp, 0.0_wp, 0.0_wp, &
                                          1.43_wp, 0.0_wp, 1.11_wp, &
                                          -1.43_wp, 0.0_wp, 1.11_wp], [3, 3]))
         call new_molecular_grid_uniform(mg, mol, nufft_nrad, nufft_nang, merr, rmax=nufft_rmax)
         if (.not. allocated(merr)) call molecular_grid_set_kgrid(mg, nufft_dr, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
         else
            do k = 1, size(nufft_widths)
               call check_threaded_nufft(error, mg, mol%xyz(:, 1), use_type12, nufft_widths(k))
               if (allocated(error)) exit
            end do
         end if
         call mg%destroy()
      end if

!$    call omp_set_num_threads(max_threads)
   end subroutine run_threaded_nufft

   !> Threaded forward and backward NUFFT of an `nv`-column Gaussian block
   !>
   !> - Column `iv` holds `exp(-a_iv |r - c_iv|^2)`, `a_iv = 1 + 0.3 (iv - 1)`,
   !>   centre stepped 0.2 bohr per column along (1, -1, 1) from `origin`;
   !>   every column differs in width and phase
   !> - Forward against the analytic FT (see `analytic_ft_error`)
   !> - Forward and backward against a serial run of the same transform,
   !>   relative to the serial maximum, at the requested NUFFT tolerance
   !> - Serial run under `omp_set_num_threads(1)`, the host-side throttle the
   !>   trafo must honour when it plans
   !> - Measured threaded vs serial: 7e-16 for one column (summation order in
   !>   the spreader), bitwise for three (one column per thread either way)
   !>
   !> @param[out]    error       test failure
   !> @param[in,out] mg          molecular grid with its reciprocal grid set
   !> @param[in]     origin      centre of the first Gaussian (bohr)
   !> @param[in]     use_type12  transform route, `.false.` for type 3
   !> @param[in]     nv          number of columns
   subroutine check_threaded_nufft(error, mg, origin, use_type12, nv)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular grid with its reciprocal grid set
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mg
      !> Centre of the first Gaussian
      real(wp), intent(in) :: origin(3)
      !> Transform route
      logical, intent(in) :: use_type12
      !> Number of columns
      integer, intent(in) :: nv

      real(wp), allocatable :: f(:, :), g_thr(:, :), g_ser(:, :)
      complex(wp), allocatable :: fk_thr(:, :), fk_ser(:, :)
      real(wp) :: cen(3, nv), alpha(nv), scale, dev
      integer :: iv, j, nthreads

      do iv = 1, nv
         alpha(iv) = 1.0_wp + 0.3_wp*real(iv - 1, wp)
         cen(:, iv) = origin + 0.2_wp*real(iv - 1, wp)*[1.0_wp, -1.0_wp, 1.0_wp]
      end do
      allocate (f(mg%ngrid, nv), g_thr(mg%ngrid, nv), g_ser(mg%ngrid, nv))
      allocate (fk_thr(mg%npts_k, nv), fk_ser(mg%npts_k, nv))
      do iv = 1, nv
         do j = 1, mg%ngrid
            f(j, iv) = exp(-alpha(iv)*sum((mg%xyz(:, j) - cen(:, iv))**2))
         end do
      end do

      call nufft_both_ways(error, mg, use_type12, f, fk_thr, g_thr)
      if (allocated(error)) return

      nthreads = 1
!$    nthreads = omp_get_max_threads()
!$    call omp_set_num_threads(1)
      call nufft_both_ways(error, mg, use_type12, f, fk_ser, g_ser)
!$    call omp_set_num_threads(nthreads)
      if (allocated(error)) return

      dev = 0.0_wp
      do iv = 1, nv
         dev = max(dev, analytic_ft_error(mg, fk_thr(:, iv), cen(:, iv), alpha(iv)))
      end do
      call check(error, dev <= nufft_analytic_tol, &
         & "threaded molecular forward NUFFT disagrees with the analytic Gaussian FT")
      if (allocated(error)) return

      scale = maxval(abs(fk_ser))
      dev = maxval(abs(fk_thr - fk_ser))/scale
      call check(error, ieee_is_finite(dev) .and. scale > 0.0_wp .and. dev <= default_nufft_tol, &
         & "threaded molecular forward NUFFT disagrees with the serial run")
      if (allocated(error)) return

      scale = maxval(abs(g_ser))
      dev = maxval(abs(g_thr - g_ser))/scale
      call check(error, ieee_is_finite(dev) .and. scale > 0.0_wp .and. dev <= default_nufft_tol, &
         & "threaded molecular backward NUFFT disagrees with the serial run")
   end subroutine check_threaded_nufft

   !> Forward transform of a block and backward transform of the result
   !>
   !> Fresh plans per call, so each run sorts and plans under the OpenMP
   !> settings in force at the call
   !>
   !> @param[out]    error       test failure
   !> @param[in,out] mg          molecular grid with its reciprocal grid set
   !> @param[in]     use_type12  transform route, `.false.` for type 3
   !> @param[in]     f           real-space block, ngrid x nv
   !> @param[out]    f_k         forward transform, npts_k x nv
   !> @param[out]    g           backward transform of `f_k`, ngrid x nv
   subroutine nufft_both_ways(error, mg, use_type12, f, f_k, g)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular grid with its reciprocal grid set
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mg
      !> Transform route
      logical, intent(in) :: use_type12
      !> Real-space block
      real(wp), intent(in) :: f(:, :)
      !> Forward transform
      complex(wp), intent(out) :: f_k(:, :)
      !> Backward transform of `f_k`
      real(wp), intent(out) :: g(:, :)

      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f_in(:, :), g_out(:, :)
      complex(wp), allocatable :: fk_out(:, :), fk_in(:, :)

      f_k = (0.0_wp, 0.0_wp)
      g = 0.0_wp
      call new_molecular_grid_trafo(trafo, mg, use_type12=use_type12)
      call trafo%prepare(size(f, 2), merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         call trafo%destroy()
         return
      end if
      if (trafo%is_type12 .neqv. use_type12) then
         call test_failed(error, "molecular NUFFT transform-route switch did not take effect")
         call trafo%destroy()
         return
      end if

      f_in = f
      allocate (fk_out(size(f_k, 1), size(f_k, 2)), g_out(size(g, 1), size(g, 2)))
      call trafo%fft_r2k(f_in, fk_out, merr)
      if (.not. allocated(merr)) then
         ! Backward destroys its input, so it gets a copy
         fk_in = fk_out
         call trafo%fft_k2r(fk_in, g_out, merr)
      end if
      if (allocated(merr)) then
         call test_failed(error, merr%message)
      else
         f_k = fk_out
         g = g_out
      end if
      call trafo%destroy()
   end subroutine nufft_both_ways

   !> Largest `|F(k) - F_exact(k)| / F_exact(0)` inside the ball `|k| <= pi/dr`
   !>
   !> - `F_exact(k) = (pi/a)^(3/2) exp(-k^2/(4a)) exp(-i k.(c - kref))`, the
   !>   convention of `test_forward_analytic` in the `math_grid_nufft` suite
   !> - Only the inscribed ball is resolved by the grid spacing
   !> - Non-finite deviation returns `huge(1.0_wp)`
   !>
   !> @param[in] mg     molecular grid supplying k-points and `kref`
   !> @param[in] fk     numerical transform of one column, length npts_k
   !> @param[in] cen    Gaussian centre (bohr)
   !> @param[in] alpha  Gaussian exponent (bohr^-2)
   function analytic_ft_error(mg, fk, cen, alpha) result(emax)
      !> Molecular grid supplying k-points and `kref`
      type(moist_math_grid_3d_molecular_type), intent(in) :: mg
      !> Numerical transform of one column
      complex(wp), intent(in) :: fk(:)
      !> Gaussian centre
      real(wp), intent(in) :: cen(3)
      !> Gaussian exponent
      real(wp), intent(in) :: alpha
      !> Largest deviation relative to the peak
      real(wp) :: emax

      real(wp) :: kvec(3), f0, phase, dev
      complex(wp) :: fex
      integer :: j

      f0 = (pi/alpha)**1.5_wp
      emax = 0.0_wp
      do j = 1, mg%npts_k
         kvec = mg%kpoint(j)
         if (norm2(kvec) > pi/mg%dr) cycle
         phase = dot_product(kvec, cen - mg%kref)
         fex = f0*exp(-sum(kvec**2)/(4.0_wp*alpha))*cmplx(cos(phase), -sin(phase), wp)
         dev = abs(fk(j) - fex)/f0
         if (.not. ieee_is_finite(dev)) then
            emax = huge(1.0_wp)
            return
         end if
         emax = max(emax, dev)
      end do
   end function analytic_ft_error

end module test_math_grid_3d_threaded
