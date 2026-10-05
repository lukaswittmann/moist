!> Suite proving the Cartesian batched FFT takes the threaded ducc0 path
!>
!> - Ordinary suites run inside test-drive's `!$omp parallel` region, which
!>   forces nthreads = 1 in the Cartesian FFT
!> - Each test here is invoked directly (`<suite> <test>`), outside any team
!> - Companion of `math_grid_3d` (batched FFT) and `math_grid_nufft` (NUFFT)
!> - NUFFT: type 1/2 and type 3, both directions, one column and a batch,
!>   against the analytic Gaussian FT, a direct inverse sum and a serial run
!> - Radial transforms: one trafo per thread, cloned from a shared template
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
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_chebyshev2_type, &
      & new_chebyshev2_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_becke_type, &
      & new_becke_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, &
      & moist_math_grid_radial_recipe_type, new_radial_grid, new_uniform_radial_pair
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_recipe_override_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_grid, molecular_grid_set_kgrid, &
      & moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo, default_nufft_tol
   use test_helpers, only: get_uniform_recipe
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
   !> Acceptance bound for the backward transform against a direct inverse sum
   real(wp), parameter :: inverse_sum_tol = 100.0_wp*default_nufft_tol

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
                  new_unittest("molecular_nufft_type3_threaded", test_nufft_type3_threaded), &
                  new_unittest("radial_trafo_one_per_thread", test_radial_trafo_per_thread) &
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
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
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
         call get_uniform_recipe(recipe, overrides, nufft_nrad, nufft_nang, merr, rmax=nufft_rmax)
         if (.not. allocated(merr)) call new_molecular_grid(mg, merr, recipe=recipe, &
                                                            overrides=overrides, reciprocal=.false.)
         if (.not. allocated(merr)) call mg%update(mol, merr)
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
   !>   center stepped 0.2 bohr per column along (1, -1, 1) from `origin`;
   !>   every column differs in width and phase
   !> - Forward against the analytic FT (see `analytic_ft_error`)
   !> - Backward against a direct type-2 sum over all k-points of the threaded
   !>   spectrum, at three nodes per column: the ones nearest the Gaussian
   !>   center and nearest +-0.3 bohr from it along (1, -1, 1)
   !> - Forward and backward against a serial run of the same transform,
   !>   relative to the serial maximum, at the requested NUFFT tolerance
   !> - Serial run under `omp_set_num_threads(1)`, the host-side throttle the
   !>   trafo must honour when it plans
   !> - Measured threaded vs serial: 7e-16 for one column (summation order in
   !>   the spreader), bitwise for three (one column per thread either way)
   !>
   !> @param[out]    error       test failure
   !> @param[in,out] mg          molecular grid with its reciprocal grid set
   !> @param[in]     origin      center of the first Gaussian (bohr)
   !> @param[in]     use_type12  transform route, `.false.` for type 3
   !> @param[in]     nv          number of columns
   subroutine check_threaded_nufft(error, mg, origin, use_type12, nv)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular grid with its reciprocal grid set
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mg
      !> Center of the first Gaussian
      real(wp), intent(in) :: origin(3)
      !> Transform route
      logical, intent(in) :: use_type12
      !> Number of columns
      integer, intent(in) :: nv

      real(wp), allocatable :: f(:, :), g_thr(:, :), g_ser(:, :)
      complex(wp), allocatable :: fk_thr(:, :), fk_ser(:, :)
      real(wp) :: cen(3, nv), alpha(nv), scale, dev, phase, inverse_ref, inverse_scale
      real(wp) :: point(3), kvec(3), probe_xyz(3)
      integer :: probe, ik, ip
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

      ! Direct inverse sum at the Gaussian center and +-0.3 bohr along (1, -1, 1)
      inverse_scale = mg%dkx*mg%dky*mg%dkz/(2.0_wp*pi)**3
      scale = maxval(abs(g_thr))
      do iv = 1, nv
         do probe = 1, 3
            probe_xyz = cen(:, iv) + 0.3_wp*real(probe - 2, wp)*[1.0_wp, -1.0_wp, 1.0_wp]
            ip = minloc(sum((mg%xyz - spread(probe_xyz, 2, mg%ngrid))**2, dim=1), dim=1)
            point = mg%xyz(:, ip) - mg%kref
            inverse_ref = 0.0_wp
            do ik = 1, mg%npts_k
               kvec = mg%kpoint(ik)
               phase = dot_product(kvec, point)
               inverse_ref = inverse_ref + real(fk_thr(ik, iv)*cmplx(cos(phase), sin(phase), wp), wp)
            end do
            dev = abs(g_thr(ip, iv) - inverse_scale*inverse_ref)/scale
            call check(error, ieee_is_finite(dev) .and. scale > 0.0_wp .and. dev < inverse_sum_tol, &
               & "threaded molecular backward NUFFT disagrees with a direct inverse sum")
            if (allocated(error)) return
         end do
      end do

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
   !> @param[in] cen    Gaussian center (bohr)
   !> @param[in] alpha  Gaussian exponent (bohr^-2)
   function analytic_ft_error(mg, fk, cen, alpha) result(emax)
      !> Molecular grid supplying k-points and `kref`
      type(moist_math_grid_3d_molecular_type), intent(in) :: mg
      !> Numerical transform of one column
      complex(wp), intent(in) :: fk(:)
      !> Gaussian center
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
         if (norm2(kvec) > pi/mg%get_dr()) cycle
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

   !> One Chebyshev-II + Becke radial grid with a fixed scale
   !>
   !> @param[out] grid  Radial grid, tagged for the quadrature trafo
   !> @param[in]  n     Number of nodes
   !> @param[in]  p     Becke scale
   !> @param[out] merr  Construction error
   subroutine chebyshev_becke_grid(grid, n, p, merr)
      !> Radial grid
      type(moist_math_grid_radial_type), intent(out) :: grid
      !> Number of nodes
      integer, intent(in) :: n
      !> Becke scale
      real(wp), intent(in) :: p
      !> Construction error
      type(mctc_error), allocatable, intent(out) :: merr

      type(moist_math_grid_radial_rule_chebyshev2_type) :: rule
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_recipe_type) :: recipe

      call new_chebyshev2_rule(rule)
      call new_becke_mapping(becke, merr, scale=p)
      if (allocated(merr)) return
      allocate (recipe%rule, source=rule)
      allocate (recipe%mapping, source=becke)
      recipe%npts = n
      call new_radial_grid(grid, recipe, 1, merr)
   end subroutine chebyshev_becke_grid

   !> Worker of the radial per-thread test: clone the template once, transform every nteam-th column
   !>
   !> @param[in]     template  Shared template trafo, read only
   !> @param[in]     it        Worker index, 1..nteam
   !> @param[in]     nteam     Number of workers
   !> @param[in]     f         r-space fields, shape (nr, nb)
   !> @param[in,out] g         k-space results, shape (nk, nb); this worker's columns written
   !> @param[in,out] h         r-space results of the forward adjoint, shape (nr, nb)
   !> @param[out]    stat      0 on success
   subroutine radial_trafo_worker(template, it, nteam, f, g, h, stat)
      !> Shared template trafo
      class(moist_math_grid_radial_trafo_type), intent(in) :: template
      !> Worker index
      integer, intent(in) :: it
      !> Number of workers
      integer, intent(in) :: nteam
      !> r-space fields
      real(wp), intent(in) :: f(:, :)
      !> k-space results
      real(wp), intent(inout) :: g(:, :)
      !> r-space results of the forward adjoint
      real(wp), intent(inout) :: h(:, :)
      !> 0 on success
      integer, intent(out) :: stat

      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      integer :: j

      stat = 0
      allocate (trafo, source=template)
      do j = it, size(f, 2), nteam
         call trafo%fbt_r2k(f(:, j), g(:, j), merr)
         if (.not. allocated(merr)) call trafo%fbt_r2k_adj(g(:, j), h(:, j), merr)
         if (allocated(merr)) then
            stat = 1
            return
         end if
      end do
   end subroutine radial_trafo_worker

   !> One radial trafo per OpenMP thread, cloned from a shared template, matches the serial result
   !>
   !> - DST-IV on a uniform pair and the quadrature trafo on a
   !>   Chebyshev-II + Becke pair
   !> - The same trafo code runs serially and per thread, so the results
   !>   must be equal, not merely close
   subroutine test_radial_trafo_per_thread(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nteam = 4, nb = 12
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: template
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :), g(:, :), h(:, :), gref(:, :), href(:, :)
      integer :: icase, it, i, j, max_threads, stat(nteam)

      call enter_threaded(error, max_threads)
      if (allocated(error)) return

      do icase = 1, 2
         if (icase == 1) then
            call new_uniform_radial_pair(rgrid, kgrid, 96, 0.1_wp, merr)
         else
            call chebyshev_becke_grid(rgrid, 64, 1.2_wp, merr)
            if (.not. allocated(merr)) call chebyshev_becke_grid(kgrid, 80, 1.5_wp, merr)
         end if
         if (.not. allocated(merr)) call new_radial_trafo(template, rgrid, kgrid, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            exit
         end if

         allocate (f(template%nr, nb), g(template%nk, nb), h(template%nr, nb))
         allocate (gref(template%nk, nb), href(template%nr, nb))
         do j = 1, nb
            do i = 1, template%nr
               f(i, j) = sin(0.37_wp*real(i, wp) + real(j, wp))
            end do
         end do
         do j = 1, nb
            call template%fbt_r2k(f(:, j), gref(:, j), merr)
            if (.not. allocated(merr)) call template%fbt_r2k_adj(gref(:, j), href(:, j), merr)
            if (allocated(merr)) exit
         end do
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            exit
         end if

         !$omp parallel do schedule(static) num_threads(nteam) default(none) &
         !$omp shared(template, f, g, h, stat)
         do it = 1, nteam
            call radial_trafo_worker(template, it, nteam, f, g, h, stat(it))
         end do
         !$omp end parallel do

         call check(error, all(stat == 0), "Per-thread radial transform reported an error")
         if (allocated(error)) exit
         call check(error, all(g == gref) .and. all(h == href), &
            & "Per-thread radial trafo clones do not reproduce the serial result")
         if (allocated(error)) exit
         deallocate (template, f, g, h, gref, href)
      end do

!$    call omp_set_num_threads(max_threads)
   end subroutine test_radial_trafo_per_thread

end module test_math_grid_3d_threaded
