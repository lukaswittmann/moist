!> Test suite for the 3D volume grids and their Fourier transform engines
!>
!> - Volume adjoint: storage and channel accumulation
!> - Cartesian grid: measure, integrate, destroy
!> - Cartesian transform: FFT round trip, batched call, inside an OpenMP team
!> - Analytic Gaussian: forward, inverse and spectral gradient vs the exact pair
!> - Threaded ducc0 backend: see the `math_grid_3d_threaded` suite
!> - Cross-validation: molecular NUFFT vs Cartesian FFT of an LJ + Yukawa field
!> - Geometry updates: motion, stretches, rebuild, pruning, copies, widths
!> - Settings: rebuild, reciprocal reference, recovery after a rejected update
!> - Error paths (`should_fail`): invalid settings, calls out of order,
!>   non-finite input, mis-sized blocks
module test_math_grid_3d
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use mctc_env, only: wp, fatal_error
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: center_at_origin, get_uniform_recipe, get_handymod_recipe, &
      & get_qc_handymod_recipe
   use moist_math_fft, only: moist_fft_r2c_3d, moist_fft_r2c_3d_batch
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_grid_3d
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_recipe_override_type, default_element_recipes
   use moist_math_grid_3d_partition, only: becke_partition_weights
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_grid, molecular_grid_set_kgrid, moist_math_grid_3d_molecular_trafo_type, &
      & new_molecular_grid_trafo
   implicit none(type, external)
   private

   public :: collect_math_grid_3d

   !> Exponent of the analytic test Gaussian exp(-a |r - c|^2), 1/bohr^2
   real(wp), parameter :: gauss_a = 1.0_wp
   !> Off-center, off-axis Gaussian center c, bohr
   real(wp), parameter :: gauss_c(3) = [0.3_wp, -0.5_wp, 0.7_wp]

   !> Requested NUFFT accuracy for the copied-domain DC identities
   real(wp), parameter :: NUFFT_COPY_TOL = 1.0e-13_wp
   !> DC residual bound, relative to the l1 norm of the transform input
   real(wp), parameter :: NUFFT_COPY_THR = 1.0e-12_wp

contains

   !> Collect all math_grid_3d tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_3d(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("volume_adjoint_storage", test_volume_adjoint), &
                  new_unittest("measure_is_cell_volume", test_measure), &
                  new_unittest("integrate_constant_is_box_volume", test_integrate_const), &
                  new_unittest("destroy_idempotent", test_destroy), &
                  new_unittest("fft_round_trip_identity", test_fft_round_trip), &
                  new_unittest("fft_batched_matches_single", test_fft_batched), &
                  new_unittest("fft_gaussian_analytic", test_fft_gaussian), &
                  new_unittest("ifft_gaussian_analytic", test_ifft_gaussian), &
                  new_unittest("fft_gaussian_spectral_gradient", test_fft_gaussian_gradient), &
                  new_unittest("ft_lj_coulomb_nufft_vs_fft", test_ft_lj_coulomb), &
                  new_unittest("cartesian_motion_fft", test_cartesian), &
                  new_unittest("molecular_motion_period", test_molecular_motion), &
                  new_unittest("molecular_nufft_copy", test_molecular_trafo), &
                  new_unittest("molecular_gaussian_widths", test_molecular_gaussian), &
                  new_unittest("molecular_pruning_deformation", test_molecular_pruning), &
                  new_unittest("direct_grid_updates", test_direct_updates), &
                  new_unittest("direct_grid_period", test_direct_period), &
                  new_unittest("grid_settings", test_grid_settings), &
                  new_unittest("molecular_update_recovers_after_error", test_update_recovers), &
                  new_unittest("molecular_rebuild_no_kgrid", test_rebuild_no_kgrid), &
                  new_unittest("molecular_destroy_resets_reference", test_destroy_kref), &
                  new_unittest("volume_adjoint_bad_uninitialized", test_volume_adjoint_uninit_fails, &
                     & should_fail=.true.), &
                  new_unittest("volume_adjoint_bad_xyz_shape", test_volume_adjoint_xyz_fails, &
                     & should_fail=.true.), &
                  new_unittest("volume_adjoint_bad_weight_size", test_volume_adjoint_w_fails, &
                     & should_fail=.true.), &
                  new_unittest("volume_adjoint_bad_xi_size", test_volume_adjoint_xi_fails, &
                     & should_fail=.true.), &
                  new_unittest("cartesian_bad_size", test_cartesian_size_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_point_count", test_cartesian_count_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_spacing", test_cartesian_dr_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_spacing_nan", test_cartesian_dr_nan_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_xi0_factor", test_cartesian_xi0_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_origin", test_cartesian_origin_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_no_atoms", test_cartesian_no_atoms_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_coordinates", test_cartesian_coords_fails, should_fail=.true.), &
                  new_unittest("cartesian_bad_rebuild_before_update", test_cartesian_rebuild_fails, &
                     & should_fail=.true.), &
                  new_unittest("cartesian_bad_trafo_before_update", test_cartesian_trafo_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_bad_ssf_a", test_molecular_ssf_a_fails, should_fail=.true.), &
                  new_unittest("molecular_bad_pruning_threshold", test_molecular_pruning_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_bad_rebuild_before_update", test_molecular_rebuild_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_bad_coordinates", test_molecular_coords_fails, should_fail=.true.), &
                  new_unittest("molecular_bad_reciprocal_count", test_molecular_kcount_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_bad_point_count", test_point_count_overflow_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_bad_nonpositive_handymod_m", test_handymod_nonpositive_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_bad_huge_handymod_m", test_handymod_huge_fails, should_fail=.true.), &
                  new_unittest("molecular_trafo_bad_unprepared_forward", test_trafo_unprepared_r2k_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_trafo_bad_unprepared_backward", test_trafo_unprepared_k2r_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_trafo_bad_batch_width", test_trafo_width_fails, should_fail=.true.), &
                  new_unittest("molecular_trafo_bad_block_size", test_trafo_block_fails, should_fail=.true.) &
                  ]
   end subroutine collect_math_grid_3d

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !> - Callers must therefore return cleanly (never `test_failed`) when their
   !>   own setup fails, or a broken fixture would pass as the expected error
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(mctc_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> Volume adjoint allocates, accumulates and zeroes every channel
   !>
   !> - init sizes all channels and zeroes them; a second init resizes them
   !> - add_weights sums repeated contributions channel by channel and leaves
   !>   omitted channels untouched
   !> - zero clears every channel
   !> - A channel whose shape disagrees with the others is not initialized
   subroutine test_volume_adjoint(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp) :: w_xyz(3, 4), w_w(4), w_xi(4)
      integer :: i

      call check(error, .not. acc%is_initialized(), "default volume adjoint claims to be initialized")
      if (allocated(error)) return
      call acc%init(4)
      call check(error, acc%is_initialized() .and. acc%size() == 4)
      if (allocated(error)) return
      call check(error, all(acc%w_xyz == 0.0_wp) .and. all(acc%w_w == 0.0_wp) .and. all(acc%w_xi == 0.0_wp))
      if (allocated(error)) return

      w_xyz = reshape([(0.5_wp*real(i, wp), i=1, 12)], [3, 4])
      w_w = [1.0_wp, -2.0_wp, 3.0_wp, -4.0_wp]
      w_xi = [0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp]
      call acc%add_weights(err, w_xyz=w_xyz, w_w=w_w, w_xi=w_xi)
      call require_success(error, err)
      if (allocated(error)) return
      call acc%add_weights(err, w_w=w_w)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, all(acc%w_xyz == w_xyz) .and. all(acc%w_w == 2.0_wp*w_w) &
         & .and. all(acc%w_xi == w_xi), "add_weights did not accumulate channel by channel")
      if (allocated(error)) return

      call acc%add_weights(err, w_xyz=w_xyz, w_xi=w_xi)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, all(acc%w_xyz == 2.0_wp*w_xyz) .and. all(acc%w_w == 2.0_wp*w_w) &
         & .and. all(acc%w_xi == 2.0_wp*w_xi), "position and width adjoints must sum repeated contributions")
      if (allocated(error)) return
      call acc%zero()
      call check(error, all(acc%w_xyz == 0.0_wp) .and. all(acc%w_w == 0.0_wp) .and. all(acc%w_xi == 0.0_wp))
      if (allocated(error)) return
      call acc%init(2)
      call check(error, acc%is_initialized() .and. acc%size() == 2, "reinit must resize all channels")
      if (allocated(error)) return
      call check(error, all(acc%w_xyz == 0.0_wp) .and. all(acc%w_w == 0.0_wp) .and. all(acc%w_xi == 0.0_wp))
      if (allocated(error)) return
      acc%w_xyz = w_xyz
      call check(error, .not. acc%is_initialized(), "inconsistent position channel shape must be rejected")
      if (allocated(error)) return
      call acc%init(2)
      acc%w_w = w_w
      call check(error, .not. acc%is_initialized(), "inconsistent integration channel size must be rejected")
      if (allocated(error)) return
      call acc%init(2)
      acc%w_xi = w_xi
      call check(error, .not. acc%is_initialized(), "inconsistent width channel size must be rejected")
   end subroutine test_volume_adjoint

   !> Constant integrand f(r) = 1 (a volume probe)
   !>
   !> @param[in] r  point in bohr
   pure function one(r) result(val)
      !> Point in bohr
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val
      val = 1.0_wp + 0.0_wp*r(1)
   end function one

   !> Forward transform through the abstract engine, failing the test on error
   !>
   !> @param[out]    error  test error
   !> @param[in,out] trafo  transform engine, prepared where the backend needs it
   !> @param[in,out] f_r    real-space block
   !> @param[out]    f_k    reciprocal-space block
   subroutine forward(error, trafo, f_r, f_k)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Transform engine (prepared where the backend needs it)
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Real-space block
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      !> Reciprocal-space block
      complex(wp), intent(out), contiguous, target :: f_k(:, :)

      type(mctc_error), allocatable :: merr

      call trafo%fft_r2k(f_r, f_k, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine forward

   !> Backward transform through the abstract engine, failing the test on error
   !>
   !> @param[out]    error  test error
   !> @param[in,out] trafo  transform engine
   !> @param[in,out] f_k    reciprocal-space block, destroyed
   !> @param[out]    f_r    real-space block
   subroutine backward(error, trafo, f_k, f_r)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Reciprocal-space block (destroyed)
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      !> Real-space block
      real(wp), intent(out), contiguous, target :: f_r(:, :)

      type(mctc_error), allocatable :: merr

      call trafo%fft_k2r(f_k, f_r, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine backward

   !> Uniform cell volume dV = dr**3 from measure(i) at any index; point(i) returns the stored node
   subroutine test_measure(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_cartesian_grid_3d(grid, 4, 5, 6, 0.3_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call check(error, abs(grid%measure(1) - grid%dv) < epsilon(1.0_wp))
      if (allocated(error)) return
      call check(error, abs(grid%measure(grid%ngrid) - 0.3_wp**3) < 1.0e-14_wp)
      if (allocated(error)) return
      call check(error, all(grid%point(1) == grid%xyz(:, 1)) &
         & .and. all(grid%point(grid%ngrid) == grid%xyz(:, grid%ngrid)), "point must return stored coordinates")
      call grid%destroy()
   end subroutine test_measure

   !> Integral of 1 dV over the box equals the box volume, dispatched polymorphically
   !>
   !> - Field integral of 1 + x^2 over the 8^3 box of spacing 0.25 bohr, nodes
   !>   at the cell centers +-0.125 .. +-0.875 bohr: box volume 8 times
   !>   (1 + mean x^2) = 8*(1 + 0.328125) = 10.625
   subroutine test_integrate_const(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_type), pointer :: g
      real(wp) :: res
      type(mctc_error), allocatable :: merr

      call new_cartesian_grid_3d(grid, 8, 8, 8, 0.25_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      g => grid
      call g%integrate(one, res)
      call check(error, abs(res - grid%vbox) < 1.0e-10_wp*grid%vbox)
      if (allocated(error)) return
      call check(error, g%ngrid == grid%ngrid)
      if (allocated(error)) return
      call g%integrate_field(1.0_wp + grid%xyz(1, :)**2, res)
      call check(error, abs(res - 10.625_wp) < 1.0e-12_wp, &
         & "Cartesian field integral must use the cell volume")
      call grid%destroy()
   end subroutine test_integrate_const

   !> Repeated destroy() leaves the grid empty and releases every geometry channel
   subroutine test_destroy(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_cartesian_grid_3d(grid, 4, 4, 4, 0.5_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call grid%destroy()
      call grid%destroy()
      call check(error, grid%ngrid == 0)
      if (allocated(error)) return
      call check(error, .not. allocated(grid%xyz) .and. .not. allocated(grid%w) &
         & .and. .not. allocated(grid%owner) .and. .not. allocated(grid%xi0), "destroy must release every geometry channel")
   end subroutine test_destroy

   !> Round trip fft_k2r(fft_r2k(f)) reproduces f to machine precision
   !>
   !> - Batched interface operates on (ntot, nv) blocks, here nv = 1
   subroutine test_fft_round_trip(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: f_r(:, :), f_back(:, :)
      complex(wp), allocatable :: f_k(:, :)
      integer :: i
      type(mctc_error), allocatable :: merr

      call new_cartesian_grid_3d(grid, 6, 6, 6, 0.4_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call grid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      allocate (f_r(grid%ngrid, 1), f_back(grid%ngrid, 1), f_k(grid%npts_k, 1))
      do i = 1, grid%ngrid
         f_r(i, 1) = sin(0.1_wp*real(i, wp)) + 0.5_wp
      end do

      call forward(error, trafo, f_r, f_k)
      if (allocated(error)) return
      call backward(error, trafo, f_k, f_back)   ! note: overwrites f_k
      if (allocated(error)) return
      call check(error, maxval(abs(f_back - f_r)) < 1.0e-10_wp)

      deallocate (f_r, f_back, f_k)
      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fft_round_trip

   !> Batched multi-site transform matches the single-field backend
   !>
   !> - Odd and even non-cubic grids
   !> - An existing outer OpenMP team must not spawn a second FFT worker team
   !> - Runs inside test-drive's own `!$omp parallel do` (see CONTRIBUTING.md:47),
   !>   so `cartesian_trafo_fft_r2k`/`fft_k2r` always see `omp_in_parallel() == .true.`
   !>   and force nthreads = 1 regardless of `omp_set_num_threads`
   !> - Looping over a requested thread count gains nothing here, every value
   !>   would hit the same single-field fallback path
   !> - Threaded ducc0 backend (nthreads > 1, auto-selected outside any enclosing
   !>   team) is proven in the `math_grid_3d_threaded` suite, invoked through
   !>   `run_selected` so it never enters that region
   !> - Nested-team check below passes nthreads=4 itself instead of relying on
   !>   auto-detection, so the restriction does not affect it and it stays here
   subroutine test_fft_batched(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: cr(:, :), tr(:, :)
      complex(wp), allocatable :: ck(:, :), ck_ref(:, :)
      integer, allocatable :: stat_batch(:)
      integer :: nx, ns, i, a, status

      do nx = 15, 16
         do ns = 1, 3, 2
            call new_cartesian_grid_3d(grid, nx, 14, 13, 0.4_wp, error=merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message); return
            end if
            call grid%new_trafo(trafo, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message); return
            end if
            allocate (cr(grid%ngrid, ns), tr(grid%ngrid, ns))
            allocate (ck(grid%npts_k, ns), ck_ref(grid%npts_k, ns))
            allocate (stat_batch(ns))
            do a = 1, ns
               do i = 1, grid%ngrid
                  cr(i, a) = sin(0.013_wp*i*a) + cos(0.007_wp*i)
               end do
               status = moist_fft_r2c_3d(grid%nz, grid%ny, grid%nx, cr(:, a), ck_ref(:, a), grid%dv)
               call check(error, status == 0, "reference forward FFT backend call failed")
               if (allocated(error)) return
            end do

            call forward(error, trafo, cr, ck)
            if (allocated(error)) return
            call check(error, maxval(abs(ck - ck_ref)) < 1.0e-11_wp, &
               & "batched forward transform deviates from the single-field backend call")
            if (allocated(error)) return

            ! An existing outer team must not create another FFT worker team
            !$omp parallel do schedule(static)
            do a = 1, ns
               stat_batch(a) = moist_fft_r2c_3d_batch(grid%nz, grid%ny, grid%nx, 1, &
                  & cr(:, a), ck(:, a), grid%dv, 4)
            end do
            !$omp end parallel do
            call check(error, all(stat_batch == 0), &
               & "batched FFT backend call failed inside an outer team")
            if (allocated(error)) return
            call check(error, maxval(abs(ck - ck_ref)) < 1.0e-11_wp, &
               & "batched transform inside an outer team deviates from the reference")
            if (allocated(error)) return

            call backward(error, trafo, ck, tr)
            if (allocated(error)) return
            call check(error, maxval(abs(tr - cr)) < 1.0e-12_wp, &
               & "batched backward transform does not invert the forward one")
            if (allocated(error)) return

            deallocate (cr, tr, ck, ck_ref, stat_batch)
            call trafo%destroy()
            call grid%destroy()
         end do
      end do
   end subroutine test_fft_batched

   !> Off-center test Gaussian exp(-a |r - c|^2)
   !>
   !> @param[in] r  point, bohr
   pure function gaussian(r) result(val)
      !> Point, bohr
      real(wp), intent(in) :: r(3)
      !> Field value
      real(wp) :: val

      val = exp(-gauss_a*sum((r - gauss_c)**2))
   end function gaussian

   !> Continuous Fourier transform of the test Gaussian, phase-referenced to r1
   !>
   !> - F(k) = (pi/a)^(3/2) exp(-|k|^2/(4a)) exp(-i k.c) in the convention
   !>   F(k) = integral f(r) exp(-i k.r) dV
   !> - Cartesian FFT phases are referenced to grid point 1, so the grid
   !>   transform equals F(k) exp(i k.r1)
   !>
   !> @param[in] k   wave vector, 1/bohr
   !> @param[in] r1  phase reference, bohr
   pure function gaussian_ft(k, r1) result(val)
      !> Wave vector, 1/bohr
      real(wp), intent(in) :: k(3)
      !> Phase reference, bohr
      real(wp), intent(in) :: r1(3)
      !> Transform value
      complex(wp) :: val

      real(wp) :: ph

      ph = -dot_product(k, gauss_c - r1)
      val = (pi/gauss_a)**1.5_wp*exp(-dot_product(k, k)/(4.0_wp*gauss_a)) &
         & *cmplx(cos(ph), sin(ph), wp)
   end function gaussian_ft

   !> Mixed odd/even Cartesian grid and engine for the analytic Gaussian tests
   !>
   !> - 45 x 48 x 47 points at dr = 0.25 bohr
   !> - Odd R2C-halved x axis, even y axis with a Nyquist plane
   !> - Negative frequencies on y and z
   !> - Aliasing error ~ exp(-(pi/dr)^2/(4a)) ~ 1e-17 of F(0)
   !> - Truncation error ~ exp(-a (L/2 - |c_i|)^2), largest on the z faces at
   !>   ~ 2e-12
   !>
   !> @param[out] grid   Cartesian grid, box centered on zero
   !> @param[out] trafo  transform engine bound to grid
   !> @param[out] error  test failure
   subroutine gaussian_grid(grid, trafo, error)
      !> Cartesian grid, box centered on zero
      type(moist_math_grid_3d_cartesian_type), intent(out), target :: grid
      !> Transform engine bound to grid
      class(moist_math_grid_3d_trafo_type), allocatable, intent(out) :: trafo
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: merr

      call new_cartesian_grid_3d(grid, 45, 48, 47, 0.25_wp, error=merr)
      if (.not. allocated(merr)) call grid%new_trafo(trafo, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine gaussian_grid

   !> Forward FFT of an off-center Gaussian equals its analytic transform
   !>
   !> - Compared as complex values on every k-point of the half spectrum
   !> - Pins the kpoint() layout and signs on all three axes, the phase
   !>   reference (grid point 1 = origin + dr/2) and the dV normalization
   subroutine test_fft_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: f_r(:, :)
      complex(wp), allocatable :: f_k(:, :), ref(:)
      real(wp) :: r1(3)
      integer :: i, j

      call gaussian_grid(grid, trafo, error)
      if (allocated(error)) return
      r1 = grid%origin + 0.5_wp*grid%dr
      allocate (f_r(grid%ngrid, 1), f_k(grid%npts_k, 1), ref(grid%npts_k))
      do i = 1, grid%ngrid
         f_r(i, 1) = gaussian(grid%xyz(:, i))
      end do
      do j = 1, grid%npts_k
         ref(j) = gaussian_ft(grid%kpoint(j), r1)
      end do

      call forward(error, trafo, f_r, f_k)
      if (allocated(error)) return
      call check(error, all(abs(f_k(:, 1) - ref) <= 1.0e-10_wp*(pi/gauss_a)**1.5_wp), &
         & "forward FFT of the Gaussian deviates from its analytic transform")

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fft_gaussian

   !> Inverse FFT of the analytic Gaussian spectrum reproduces the Gaussian
   !>
   !> - Only check on the absolute 1/Vbox normalization; a round trip cannot
   !>   see a common scale error of the two directions
   !> - Analytic half spectrum is Hermitian by construction, as C2R requires
   subroutine test_ifft_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: f_r(:, :), ref(:)
      complex(wp), allocatable :: f_k(:, :)
      real(wp) :: r1(3)
      integer :: i, j

      call gaussian_grid(grid, trafo, error)
      if (allocated(error)) return
      r1 = grid%origin + 0.5_wp*grid%dr
      allocate (f_r(grid%ngrid, 1), f_k(grid%npts_k, 1), ref(grid%ngrid))
      do j = 1, grid%npts_k
         f_k(j, 1) = gaussian_ft(grid%kpoint(j), r1)
      end do
      do i = 1, grid%ngrid
         ref(i) = gaussian(grid%xyz(:, i))
      end do

      call backward(error, trafo, f_k, f_r)
      if (allocated(error)) return
      call check(error, all(abs(f_r(:, 1) - ref) <= 1.0e-10_wp), &
         & "inverse FFT of the analytic spectrum deviates from the Gaussian")

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_ifft_gaussian

   !> Spectral gradient of the Gaussian matches its analytic gradient
   !>
   !> - d/dx_d f = k2r(i k_d r2k(f)) against -2a (x_d - c_d) f, per axis
   !> - A sign or ordering error in any axis's frequencies flips or scrambles
   !>   that gradient component
   !> - Even-axis Nyquist modes break Hermitian symmetry under i k_d, but the
   !>   Gaussian spectrum there is ~ 1e-17 of F(0)
   subroutine test_fft_gaussian_gradient(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: f_r(:, :), g_r(:, :), ref(:)
      complex(wp), allocatable :: f_k(:, :), g_k(:, :)
      real(wp) :: k(3)
      integer :: i, j, d
      character(len=*), parameter :: axis = "xyz"

      call gaussian_grid(grid, trafo, error)
      if (allocated(error)) return
      allocate (f_r(grid%ngrid, 1), g_r(grid%ngrid, 1), ref(grid%ngrid))
      allocate (f_k(grid%npts_k, 1), g_k(grid%npts_k, 1))
      do i = 1, grid%ngrid
         f_r(i, 1) = gaussian(grid%xyz(:, i))
      end do
      call forward(error, trafo, f_r, f_k)
      if (allocated(error)) return

      do d = 1, 3
         do j = 1, grid%npts_k
            k = grid%kpoint(j)
            g_k(j, 1) = cmplx(0.0_wp, k(d), wp)*f_k(j, 1)
         end do
         call backward(error, trafo, g_k, g_r)
         if (allocated(error)) return
         do i = 1, grid%ngrid
            ref(i) = -2.0_wp*gauss_a*(grid%xyz(d, i) - gauss_c(d))*gaussian(grid%xyz(:, i))
         end do
         call check(error, all(abs(g_r(:, 1) - ref) <= 1.0e-9_wp), &
            & "spectral gradient along "//axis(d:d)//" deviates from the analytic gradient")
         if (allocated(error)) return
      end do

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fft_gaussian_gradient

   !> Small H2 molecular grid with a reciprocal grid and an unprepared engine
   !>
   !> @param[out] mgrid   molecular grid, must outlive mtrafo
   !> @param[out] mtrafo  unprepared engine bound to mgrid
   !> @param[out] err     library error from the construction
   subroutine molecular_trafo_fixture(mgrid, mtrafo, err)
      !> Molecular grid, must outlive mtrafo
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mgrid
      !> Unprepared engine bound to mgrid
      type(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: mtrafo
      !> Library error from the construction
      type(mctc_error), allocatable, intent(out) :: err

      type(structure_type) :: mol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call get_uniform_recipe(recipe, overrides, 12, 26, err, rmax=4.0_wp)
      if (allocated(err)) return
      call new_molecular_grid(mgrid, err, recipe=recipe, overrides=overrides, reciprocal=.false.)
      if (allocated(err)) return
      call mgrid%update(mol, err)
      if (allocated(err)) return
      call molecular_grid_set_kgrid(mgrid, 0.8_wp, err)
      if (allocated(err)) return
      call new_molecular_grid_trafo(mtrafo, mgrid)
   end subroutine molecular_trafo_fixture

   !> Forward transform on an unprepared molecular engine fails
   subroutine test_trafo_unprepared_r2k_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call molecular_trafo_fixture(mgrid, mtrafo, err)
      if (allocated(err)) return
      allocate (f(mgrid%ngrid, 1), fk(mgrid%npts_k, 1))
      f = 1.0_wp
      call mtrafo%fft_r2k(f, fk, err)
      call expect_error(error, err, "call prepare(nv)")
      call mtrafo%destroy()
   end subroutine test_trafo_unprepared_r2k_fails

   !> Backward transform on an unprepared molecular engine fails
   subroutine test_trafo_unprepared_k2r_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call molecular_trafo_fixture(mgrid, mtrafo, err)
      if (allocated(err)) return
      allocate (f(mgrid%ngrid, 1), fk(mgrid%npts_k, 1))
      fk = (0.0_wp, 0.0_wp)
      call mtrafo%fft_k2r(fk, f, err)
      call expect_error(error, err, "call prepare(nv)")
      call mtrafo%destroy()
   end subroutine test_trafo_unprepared_k2r_fails

   !> Molecular engine prepared for one column fails on a two-column block
   subroutine test_trafo_width_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call molecular_trafo_fixture(mgrid, mtrafo, err)
      if (.not. allocated(err)) call mtrafo%prepare(1, err)
      if (allocated(err)) return
      allocate (f(mgrid%ngrid, 2), fk(mgrid%npts_k, 2))
      f = 1.0_wp
      call mtrafo%fft_r2k(f, fk, err)
      call expect_error(error, err, "batch width")
      call mtrafo%destroy()
   end subroutine test_trafo_width_fails

   !> Molecular engine fails on a real-space block longer than the grid
   subroutine test_trafo_block_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call molecular_trafo_fixture(mgrid, mtrafo, err)
      if (.not. allocated(err)) call mtrafo%prepare(1, err)
      if (allocated(err)) return
      allocate (f(mgrid%ngrid + 1, 1), fk(mgrid%npts_k, 1))
      f = 1.0_wp
      call mtrafo%fft_r2k(f, fk, err)
      call expect_error(error, err, "field block size")
      call mtrafo%destroy()
   end subroutine test_trafo_block_fails

   !> Lennard-Jones + screened-Coulomb (Yukawa) field at `r`, summed over sources
   !>
   !> - Sources at `coord(:,a)` with per-atom well depth `eps`, size `sig`, charge `q`
   !> - Source distance softened to `dsoft`, so the field is bounded at the grid
   !>   points nearest a core
   !> - Coulomb part screened by `kappa`, so the field is compact
   !> - Together they keep the transform well resolved on both the uniform
   !>   Cartesian and the atom-centered grid, which makes a cross-grid comparison
   !>   meaningful
   !>
   !> @param[in] r      field point in bohr
   !> @param[in] coord  source coordinates, shape (3, nat)
   !> @param[in] eps    per-source LJ well depth
   !> @param[in] sig    per-source LJ size
   !> @param[in] q      per-source charge
   !> @param[in] dsoft  distance softening in bohr
   !> @param[in] kappa  Coulomb screening in 1/bohr
   pure function lj_screened_coulomb(r, coord, eps, sig, q, dsoft, kappa) result(val)
      !> Field point (bohr)
      real(wp), intent(in) :: r(3)
      !> Source coordinates, shape (3, nat)
      real(wp), intent(in) :: coord(:, :)
      !> Per-source LJ well depth / size and charge
      real(wp), intent(in) :: eps(:), sig(:), q(:)
      !> Distance softening and Coulomb screening (bohr, 1/bohr)
      real(wp), intent(in) :: dsoft, kappa
      !> Field value at r
      real(wp) :: val

      integer :: a
      real(wp) :: d, sr6

      val = 0.0_wp
      do a = 1, size(eps)
         d = sqrt((r(1) - coord(1, a))**2 + (r(2) - coord(2, a))**2 &
                  + (r(3) - coord(3, a))**2)
         d = max(d, dsoft)
         sr6 = (sig(a)/d)**6
         val = val + 4.0_wp*eps(a)*(sr6*sr6 - sr6) + q(a)*exp(-kappa*d)/d
      end do
   end function lj_screened_coulomb

   !> Finiteness of a complex value, true iff both its parts are finite
   !>
   !> - Guards the `max(maxerr, abs(...))` accumulators below
   !> - A NaN sample would otherwise vanish silently, an IEEE comparison against
   !>   NaN is always false so gfortran's `max` can keep the running finite
   !>   accumulator
   !>
   !> @param[in] z  value to test
   pure function complex_is_finite(z) result(ok)
      !> Value to test
      complex(wp), intent(in) :: z
      !> True if both the real and imaginary parts are finite
      logical :: ok

      ok = ieee_is_finite(real(z, wp)) .and. ieee_is_finite(aimag(z))
   end function complex_is_finite

   !> Cross-validate the molecular-grid NUFFT against the Cartesian-grid FFT
   !>
   !> Common Lennard-Jones + screened-Coulomb field sourced at the atoms
   !> - Molecular NUFFT forward transform reproduces the brute-force direct DFT
   !>   `sum_j w_j f(r_j) exp(-i k.(r_j-kref))` at sampled k-points (tight,
   !>   FINUFFT correctness check independent of grid accuracy)
   !> - DC mode (k=0) agrees between the two grids, both equal the grid
   !>   quadrature of the field, `integral f dV`
   !> - At the lowest k along an axis the Cartesian FFT magnitude matches the
   !>   molecular grid's transform magnitude (same Fourier coefficient from two
   !>   independent grids and transforms, to grid accuracy)
   subroutine test_ft_lj_coulomb(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(moist_math_grid_3d_cartesian_type), target :: cgrid
      class(moist_math_grid_3d_trafo_type), allocatable :: ctrafo
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      !> Softened LJ + screened Coulomb model parameters (see field helper)
      real(wp), parameter :: dsoft = 3.5_wp, kappa = 0.5_wp
      integer, parameter :: nsamp = 8
      real(wp), allocatable :: eps(:), sig(:), q(:), coord(:, :)
      real(wp), allocatable :: f_mol(:, :), f_cart(:, :)
      complex(wp), allocatable :: fk_mol(:, :), fk_cart(:, :)
      real(wp) :: kvec(3), ph, fscale, maxerr, dc_mol, dc_cart, dcref
      complex(wp) :: acc
      integer :: nat, a, j, s, jk, k0_mol, mx

      ! --- Solute structure + synthetic LJ/charge parameters ---
      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      nat = mol%nat
      allocate (coord(3, nat), eps(nat), sig(nat), q(nat))
      coord = mol%xyz
      do a = 1, nat
         eps(a) = 0.01_wp
         sig(a) = 5.0_wp
         q(a) = real(1 - 2*mod(a, 2), wp)*0.5_wp   ! alternating +/- 0.5 charges
      end do

      ! --- Molecular grid + NUFFT of the field ---
      call get_uniform_recipe(recipe, overrides, 40, 110, merr, rmax=8.0_wp)
      if (.not. allocated(merr)) call new_molecular_grid(mgrid, merr, recipe=recipe, &
         & overrides=overrides, reciprocal=.false.)
      if (.not. allocated(merr)) call mgrid%update(mol, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call molecular_grid_set_kgrid(mgrid, 0.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call new_molecular_grid_trafo(mtrafo, mgrid)
      call mtrafo%prepare(1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if

      allocate (f_mol(mgrid%ngrid, 1), fk_mol(mgrid%npts_k, 1))
      do j = 1, mgrid%ngrid
         f_mol(j, 1) = lj_screened_coulomb(mgrid%xyz(:, j), coord, eps, sig, q, dsoft, kappa)
      end do
      call forward(error, mtrafo, f_mol, fk_mol)
      if (allocated(error)) return
      fscale = maxval(abs(fk_mol(:, 1)))
      call check(error, fscale > 0.0_wp, "molecular NUFFT transform magnitude is zero")
      if (allocated(error)) return

      ! --- Cartesian grid (centered, enclosing the field) + FFT of same field ---
      call new_cartesian_grid_3d(cgrid, 64, 64, 64, 0.3_wp, error=merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call cgrid%new_trafo(ctrafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      allocate (f_cart(cgrid%ngrid, 1), fk_cart(cgrid%npts_k, 1))
      do j = 1, cgrid%ngrid
         f_cart(j, 1) = lj_screened_coulomb(cgrid%point(j), coord, eps, sig, q, dsoft, kappa)
      end do
      call forward(error, ctrafo, f_cart, fk_cart)
      if (allocated(error)) return

      ! --- Check 1: molecular NUFFT == brute-force direct DFT (tight) ---
      maxerr = 0.0_wp
      do s = 1, nsamp
         jk = 1 + (s - 1)*(mgrid%npts_k/nsamp)
         kvec = mgrid%kpoint(jk)
         acc = (0.0_wp, 0.0_wp)
         do j = 1, mgrid%ngrid
            ph = -dot_product(kvec, mgrid%xyz(:, j) - mgrid%kref)
            acc = acc + mgrid%w(j)*f_mol(j, 1)*cmplx(cos(ph), sin(ph), wp)
         end do
         if (.not. complex_is_finite(fk_mol(jk, 1)) .or. .not. complex_is_finite(acc)) then
            call test_failed(error, "molecular NUFFT check produced a non-finite value")
            return
         end if
         maxerr = max(maxerr, abs(fk_mol(jk, 1) - acc))
      end do
      call check(error, maxerr <= 1.0e-6_wp*fscale, &
         & "molecular NUFFT deviates from the direct DFT of the LJ+Coulomb field")

      ! --- Check 2: DC mode (k=0) agrees between the two grids ---
      if (.not. allocated(error)) then
         k0_mol = 0
         do j = 1, mgrid%npts_k
            if (maxval(abs(mgrid%kpoint(j))) < 1.0e-12_wp) then
               k0_mol = j; exit
            end if
         end do
         call check(error, k0_mol > 0, "molecular k-grid has no DC mode")
         if (.not. allocated(error)) then
            dc_mol = real(fk_mol(k0_mol, 1), wp)
            dc_cart = real(fk_cart(1, 1), wp)          ! Cartesian kpoint(1) = (0,0,0)
            call check(error, ieee_is_finite(dc_mol) .and. ieee_is_finite(dc_cart), &
               & "DC-mode comparison produced a non-finite value")
            if (.not. allocated(error)) then
               dcref = max(abs(dc_cart), 1.0e-30_wp)
               call check(error, abs(dc_mol - dc_cart) <= 2.0e-2_wp*dcref, &
                  & "molecular and Cartesian grids disagree on the field integral (DC)")
            end if
         end if
      end if

      ! --- Check 3: low-k Fourier magnitudes agree across the two grids ---
      ! Cartesian FFT value at q=(mx*dkx,0,0) (R2C layout, kx fastest) vs the
      ! molecular grid's transform at the same physical q
      ! Magnitudes are phase-reference independent; with check 1 this ties
      ! the NUFFT and FFT spectra of the same field together at low k
      if (.not. allocated(error)) then
         call check(error, abs(fk_cart(1, 1)) > 0.0_wp, "Cartesian DC Fourier coefficient is zero")
      end if
      if (.not. allocated(error)) then
         maxerr = 0.0_wp
         do mx = 1, 3
            jk = mx + 1
            kvec = cgrid%kpoint(jk)               ! = (mx*dkx, 0, 0)
            acc = (0.0_wp, 0.0_wp)
            do j = 1, mgrid%ngrid
               ph = -dot_product(kvec, mgrid%xyz(:, j))
               acc = acc + mgrid%w(j)*f_mol(j, 1)*cmplx(cos(ph), sin(ph), wp)
            end do
            if (.not. complex_is_finite(fk_cart(jk, 1)) .or. .not. complex_is_finite(acc)) then
               call test_failed(error, "low-k cross-grid check produced a non-finite value")
               return
            end if
            maxerr = max(maxerr, abs(abs(fk_cart(jk, 1)) - abs(acc)))
         end do
         call check(error, maxerr <= 3.0e-2_wp*abs(fk_cart(1, 1)), &
            & "Cartesian FFT and molecular transform disagree at low k")
      end if

      call mtrafo%destroy()
      call mgrid%destroy()
      call ctrafo%destroy()
      call cgrid%destroy()
   end subroutine test_ft_lj_coulomb

   !> Molecular grids retain their recipe across geometry changes
   !>
   !> - Every recipe: the field integral is the weighted sum over all stored nodes
   !> - Every recipe: translation moves points rigidly and keeps the weights
   !> - After a small stretch, update matches a freshly constructed grid
   !> - destroy() empties all arrays and counts
   subroutine test_direct_updates(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid, reference
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :), weights(:)
      real(wp) :: shift(3)
      integer :: recipe
      real(wp) :: integral, expected, scale

      shift = [0.2_wp, -0.3_wp, 0.1_wp]
      do recipe = 1, 5
         call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
         call construct_direct_grid(grid, mol, recipe, err)
         call require_success(error, err)
         if (allocated(error)) return
         xyz = grid%xyz
         weights = grid%w
         expected = sum(weights*(xyz(1, :) - 0.5_wp*xyz(2, :) + 0.3_wp*xyz(3, :) + 1.0_wp))
         ! The linear field changes sign, so scale by the sum of absolute contributions
         scale = sum(abs(weights*(xyz(1, :) - 0.5_wp*xyz(2, :) + 0.3_wp*xyz(3, :) + 1.0_wp)))
         call grid%integrate_field(xyz(1, :) - 0.5_wp*xyz(2, :) + 0.3_wp*xyz(3, :) + 1.0_wp, integral)
         call check(error, abs(integral - expected) < 1.0e-12_wp*max(1.0_wp, scale), &
            & "molecular field integral must use every stored quadrature weight")
         if (allocated(error)) return
         call check(error, grid%natom, 2)
         if (allocated(error)) return
         call check(error, all(grid%owner >= 1 .and. grid%owner <= 2))
         if (allocated(error)) return
         mol%xyz = mol%xyz + spread(shift, 2, mol%nat)
         call grid%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, grid%ngrid, size(weights))
         if (allocated(error)) return
         call check(error, maxval(abs(grid%xyz - xyz - spread(shift, 2, grid%ngrid))) < 1.0e-10_wp)
         if (allocated(error)) return
         call check(error, maxval(abs(grid%w - weights)) < 1.0e-10_wp*max(1.0_wp, maxval(abs(weights))))
         if (allocated(error)) return
         mol%xyz(1, 2) = mol%xyz(1, 2) + 0.2_wp
         call grid%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         call construct_direct_grid(reference, mol, recipe, err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, grid%ngrid, reference%ngrid)
         if (allocated(error)) return
         call check(error, maxval(abs(grid%xyz - reference%xyz)) < epsilon(1.0_wp))
         if (allocated(error)) return
         call check(error, maxval(abs(grid%w - reference%w)) < epsilon(1.0_wp))
         if (allocated(error)) return
         call grid%destroy()
         call check(error, grid%ngrid == 0 .and. grid%natom == 0)
         if (allocated(error)) return
         call check(error,.not. allocated(grid%xyz) .and. .not. allocated(grid%w) &
            & .and. .not. allocated(grid%owner) .and. .not. allocated(grid%xi0))
         if (allocated(error)) return
      end do
   end subroutine test_direct_updates

   !> Molecular updates retain an explicitly set reciprocal reference convention
   !>
   !> - Grid built without a reciprocal grid, then `molecular_grid_set_kgrid`
   !>   with the first-point reference and with an explicit reference
   !> - update keeps kref, mode counts and spacings, translation moves kref with
   !>   the molecule, rebuild keeps the translated kref
   subroutine test_direct_period(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol
      real(wp) :: reference(3), dk(3), shift(3)
      integer :: modes(3), policy

      shift = [0.2_wp, -0.3_wp, 0.1_wp]
      do policy = 1, 2
         call new(mol, [1, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
            & 0.0_wp, 1.0_wp, 0.0_wp, 10.0_wp, 0.0_wp, 0.0_wp], [3, 3]))
         call get_qc_handymod_recipe(recipe, err, nrad=8, degree=17, rmin=0.0_wp, rmax=4.0_wp, m=2.0_wp)
         if (.not. allocated(err)) call new_molecular_grid(grid, err, recipe=recipe, reciprocal=.false.)
         if (.not. allocated(err)) call grid%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         if (policy == 1) then
            call molecular_grid_set_kgrid(grid, 1.0_wp, err)
         else
            reference = sum(mol%xyz, dim=2)/real(mol%nat, wp) + [0.5_wp, 0.2_wp, 0.0_wp]
            call molecular_grid_set_kgrid(grid, 1.0_wp, err, reference=reference)
         end if
         call require_success(error, err)
         if (allocated(error)) return
         reference = grid%kref
         modes = [grid%nkx, grid%nky, grid%nkz]
         dk = [grid%dkx, grid%dky, grid%dkz]
         call grid%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(grid%kref - reference)) < 1.0e-14_wp)
         if (allocated(error)) return
         mol%xyz = mol%xyz + spread(shift, 2, mol%nat)
         call grid%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(grid%kref - reference - shift)) < 1.0e-14_wp)
         if (allocated(error)) return
         call check(error, all([grid%nkx, grid%nky, grid%nkz] == modes) &
            & .and. maxval(abs([grid%dkx, grid%dky, grid%dkz] - dk)) < epsilon(1.0_wp))
         if (allocated(error)) return
         call grid%rebuild(err)
         call require_success(error, err)
         if (allocated(error)) return
         call check(error, maxval(abs(grid%kref - reference - shift)) < 1.0e-14_wp)
         if (allocated(error)) return
      end do
   end subroutine test_direct_period

   !> Construct and realize molecular grids with non-default partitioning and angular settings
   !>
   !> 1. Per-element defaults, Becke k = 2, pruning threshold 0.05
   !> 2. Chebyshev-II x Becke, radial cutoffs [0.1, 5], SSF
   !> 3. Chebyshev-II x HandyMod, Becke k = 1, pruning threshold 0.1
   !> 4. Midpoint x HandyMod, constant degree 5, SSF
   !> 5. Midpoint x HandyMod, arc bands with floor and cap, Becke k = 1
   !>
   !> @param[out] grid    molecular grid, realized at `mol`
   !> @param[in]  mol     molecular structure
   !> @param[in]  recipe  recipe selection, 1 to 5
   !> @param[out] error   library error, set for an unknown recipe
   subroutine construct_direct_grid(grid, mol, recipe, error)
      !> Molecular grid
      type(moist_math_grid_3d_molecular_type), intent(out) :: grid
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Recipe selection
      integer, intent(in) :: recipe
      !> Library error
      type(mctc_error), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_recipe_type) :: atomic
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      select case (recipe)
      case (1)
         call default_element_recipes(atomic, overrides, error)
         if (.not. allocated(error)) call new_molecular_grid(grid, error, recipe=atomic, &
            & overrides=overrides, becke_k=2, pruning_threshold=0.05_wp, reciprocal=.false.)
      case (2)
         call get_uniform_recipe(atomic, overrides, 12, 26, error, rmin=0.1_wp, rmax=5.0_wp)
         if (.not. allocated(error)) call new_molecular_grid(grid, error, recipe=atomic, &
            & overrides=overrides, ssf_a=0.64_wp, reciprocal=.false.)
      case (3)
         call get_handymod_recipe(atomic, 12, 26, 0.0_wp, 5.0_wp, 2.0_wp, error)
         if (.not. allocated(error)) call new_molecular_grid(grid, error, recipe=atomic, &
            & becke_k=1, pruning_threshold=0.1_wp, reciprocal=.false.)
      case (4)
         call get_qc_handymod_recipe(atomic, error, nrad=12, degree=5, rmin=0.0_wp, rmax=5.0_wp, &
            & m=2.0_wp)
         if (.not. allocated(error)) call new_molecular_grid(grid, error, recipe=atomic, &
            & ssf_a=0.64_wp, reciprocal=.false.)
      case (5)
         call get_qc_handymod_recipe(atomic, error, nrad=12, rmin=0.0_wp, rmax=5.0_wp, m=2.0_wp, &
            & arc_r=[1.0_wp], arc_a=[0.5_wp, 1.0_wp], nang_min=110, nang_max=302)
         if (.not. allocated(error)) call new_molecular_grid(grid, error, recipe=atomic, &
            & becke_k=1, reciprocal=.false.)
      case default
         call fatal_error(error, "unknown test quadrature recipe")
      end select
      if (allocated(error)) return
      call grid%update(mol, error)
   end subroutine construct_direct_grid

   !> Forward a library error into the test framework
   !>
   !> @param[out] error  test error
   !> @param[in]  err    library error, ignored if unallocated
   subroutine require_success(error, err)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Library error
      type(mctc_error), allocatable, intent(in) :: err
      if (allocated(err)) call test_failed(error, err%message)
   end subroutine require_success

   !> Cartesian quadrature, FFT round trip and centroid translation
   !>
   !> - Gaussian integral matches pi**1.5, FFT round trip restores the field
   !> - Translation moves points rigidly, all owners are 0
   !> - Volume gradient reports its phase-2 stub through the error channel
   subroutine test_cartesian(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type), target :: domain
      type(structure_type) :: mol
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(volume_adjoint_type) :: acc
      real(wp), allocatable :: xyz(:, :), values(:, :), restored(:, :)
      complex(wp), allocatable :: spectrum(:, :)
      real(wp) :: shift(3), integral, gradient(3, 1)
      real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)

      domain%nx = 24
      domain%ny = 24
      domain%nz = 24
      domain%dr = 0.5_wp
      call domain%validate(err)
      call require_success(error, err)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, allocated(domain%xi0) .and. .not. domain%has_geometry_dependent_xi0())
      if (allocated(error)) return
      xyz = domain%xyz
      values = reshape(exp(-sum(xyz**2, dim=1)), [domain%ngrid, 1])
      integral = sum(domain%w*values(:, 1))
      call check(error, integral, pi**1.5_wp, thr=1.0e-10_wp)
      if (allocated(error)) return
      allocate (spectrum(domain%npts_k, 1), restored(domain%ngrid, 1))
      call domain%new_trafo(trafo, err)
      call require_success(error, err)
      if (allocated(error)) return
      call trafo%fft_r2k(values, spectrum, err)
      if (.not. allocated(err)) call trafo%fft_k2r(spectrum, restored, err)
      call trafo%destroy()
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(restored - values)) < 1.0e-12_wp)
      if (allocated(error)) return
      shift = [0.2_wp, -0.4_wp, 1.1_wp]
      mol%xyz(:, 1) = shift
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(domain%xyz - xyz - spread(shift, 2, domain%ngrid))) < 2.0e-15_wp)
      if (allocated(error)) return
      call check(error, all(domain%owner == 0))
      if (allocated(error)) return
      call acc%init(domain%ngrid)
      gradient = 0.0_wp
      call domain%get_volume_gradient(acc, gradient, err)
      call check(error, allocated(err), "phase-2 volume gradient must report its stub")
      if (allocated(error)) return
      call check(error, index(err%message, "cartesian") > 0)
   end subroutine test_cartesian

   !> Translation preserves quadrature and modes, stretching needs explicit rebuild
   !>
   !> - Translation keeps weights, owners and reciprocal modes
   !> - Small stretch keeps modes and moves surviving nodes with their owner
   !> - Large stretch fails naming rebuild and leaves the grid unchanged, explicit
   !>   rebuild enlarges the modes
   subroutine test_molecular_motion(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: domain
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :), weights(:)
      integer, allocatable :: owner(:)
      real(wp) :: shift(3), dk(3), center(3)
      integer :: modes(3), i
      !> Outer HandyMod radius (bohr)
      real(wp), parameter :: rmax = 4.0_wp

      call get_qc_handymod_recipe(recipe, err, nrad=12, rmax=rmax)
      if (.not. allocated(err)) call new_molecular_grid(domain, err, recipe=recipe, dr=1.0_wp, &
         & kbuffer=1.0_wp)
      call require_success(error, err)
      if (allocated(error)) return
      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      xyz = domain%xyz
      weights = domain%w
      owner = domain%owner
      modes = [domain%nkx, domain%nky, domain%nkz]
      dk = [domain%dkx, domain%dky, domain%dkz]
      call check(error, .not. allocated(domain%xi0) .and. .not. domain%has_geometry_dependent_xi0())
      if (allocated(error)) return
      shift = [0.125_wp, -0.25_wp, 0.5_wp]
      mol%xyz = mol%xyz + spread(shift, 2, mol%nat)
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, domain%ngrid, size(weights))
      if (allocated(error)) return
      call check(error, maxval(abs(domain%xyz - xyz - spread(shift, 2, domain%ngrid))) < 2.0e-15_wp &
         & .and. maxval(abs(domain%w - weights)) < 1.0e-13_wp .and. all(domain%owner == owner))
      if (allocated(error)) return
      ! A small stretch keeps modes and moves surviving nodes with their owner
      mol%xyz(1, 2) = mol%xyz(1, 2) + 0.02_wp
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, all([domain%nkx, domain%nky, domain%nkz] == modes) &
         & .and. maxval(abs([domain%dkx, domain%dky, domain%dkz] - dk)) < epsilon(1.0_wp))
      if (allocated(error)) return
      center = sum(mol%xyz, dim=2)/real(mol%nat, wp)
      call check(error, maxval(abs(domain%kref - center)) < epsilon(1.0_wp))
      if (allocated(error)) return
      call check(error, all(domain%owner >= 1 .and. domain%owner <= mol%nat))
      if (allocated(error)) return
      do i = 1, domain%ngrid
         if (sqrt(sum((domain%xyz(:, i) - mol%xyz(:, domain%owner(i)))**2)) >= rmax) then
            call test_failed(error, "molecular point does not follow its owning atom")
            return
         end if
      end do
      xyz = domain%xyz
      mol%xyz(1, 2) = mol%xyz(1, 2) + 8.0_wp
      call domain%update(mol, err)
      call check(error, allocated(err), "large stretch must not silently resize the reciprocal grid")
      if (allocated(error)) return
      call check(error, index(err%message, "molecular") > 0 .and. index(err%message, "rebuild") > 0)
      if (allocated(error)) return
      call check(error, maxval(abs(domain%xyz - xyz)) < epsilon(1.0_wp))
      if (allocated(error)) return
      call domain%rebuild(err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, domain%nkx > modes(1))
      if (allocated(error)) return
      call domain%update(mol, err)
      call require_success(error, err)
   end subroutine test_molecular_motion

   !> Copied domain owns its geometry and transforms through the fast NUFFT path
   !>
   !> - Copy keeps its kref after the source moves
   !> - Type 1/2 NUFFT, DC mode equals the weighted field sum
   !> - Inverse of a DC-only spectrum restores a constant
   !> - Both within `NUFFT_COPY_THR` of the input l1 norm
   subroutine test_molecular_trafo(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type), target :: domain, copied
      type(structure_type) :: mol
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: values(:, :), restored(:, :)
      complex(wp), allocatable :: spectrum(:, :)
      real(wp) :: volume
      integer :: zero_mode
      real(wp), parameter :: two_pi = 8.0_wp*atan(1.0_wp)
      type(moist_math_grid_atomic_recipe_type) :: recipe

      call get_qc_handymod_recipe(recipe, err, nrad=12, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(domain, err, recipe=recipe, dr=1.0_wp, &
         & nufft_tol=NUFFT_COPY_TOL)
      call require_success(error, err)
      if (allocated(error)) return
      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      copied = domain
      mol%xyz(:, 1) = [1.0_wp, 2.0_wp, 3.0_wp]
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(copied%kref)) < epsilon(1.0_wp))
      if (allocated(error)) return
      call copied%new_trafo(trafo, err)
      if (.not. allocated(err)) call trafo%prepare(1, err)
      call require_success(error, err)
      if (allocated(error)) return
      select type (trafo)
      type is (moist_math_grid_3d_molecular_trafo_type)
         call check(error, trafo%is_type12, "production domain must use NUFFT type 1/2")
      end select
      if (allocated(error)) then
         call trafo%destroy()
         return
      end if
      allocate (values(copied%ngrid, 1), restored(copied%ngrid, 1), spectrum(copied%npts_k, 1))
      values(:, 1) = exp(-sum(copied%xyz**2, dim=1))
      associate (g => copied)
         zero_mode = g%nkx/2 + 1 + g%nkx*(g%nky/2 + g%nky*(g%nkz/2))
         volume = two_pi**3/(g%dkx*g%dky*g%dkz)
      end associate
      call trafo%fft_r2k(values, spectrum, err)
      call require_success(error, err)
      if (.not. allocated(error)) then
         call check(error, abs(spectrum(zero_mode, 1) - sum(copied%w*values(:, 1))) &
            & < NUFFT_COPY_THR*sum(abs(copied%w*values(:, 1))))
      end if
      if (.not. allocated(error)) then
         spectrum = cmplx(0.0_wp, 0.0_wp, wp)
         spectrum(zero_mode, 1) = cmplx(volume, 0.0_wp, wp)
         call trafo%fft_k2r(spectrum, restored, err)
         call require_success(error, err)
         if (.not. allocated(error)) call check(error, maxval(abs(restored - 1.0_wp)) < NUFFT_COPY_THR)
      end if
      call trafo%destroy()
   end subroutine test_molecular_trafo

   !> Weight-derived widths advertise their dependence and satisfy the scale law
   !>
   !> - xi0 allocated and geometry dependent, sized as the grid
   !> - Scale law xi0**3*w = xi0_factor**3
   subroutine test_molecular_gaussian(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: domain
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol

      call get_qc_handymod_recipe(recipe, err, nrad=8, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(domain, err, recipe=recipe, gaussian=.true., &
         & xi0_factor=1.2_wp, ssf_a=0.64_wp)
      call require_success(error, err)
      if (allocated(error)) return
      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, allocated(domain%xi0) .and. domain%has_geometry_dependent_xi0())
      if (allocated(error)) return
      call check(error, maxval(abs(domain%xi0**3*domain%w - domain%get_xi0_factor()**3)) < 1.0e-13_wp)
      if (allocated(error)) return
      call check(error, size(domain%xi0), domain%ngrid)
   end subroutine test_molecular_gaussian

   !> Pruned molecular weights match a reconstruction from a fixed isolated-atom quadrature
   !>
   !> - Retained points and updated Becke weights agree with the reference at two geometries
   !> - Deformation changes the weights and crosses the pruning threshold, reciprocal modes stay fixed
   subroutine test_molecular_pruning(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: domain, raw
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol, atom
      real(wp), allocatable :: before(:), expected(:), points(:, :), partition(:)
      real(wp) :: point(3), weight
      integer :: i, j, site, step, previous_count, modes(3)

      call get_qc_handymod_recipe(recipe, err, nrad=8, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(raw, err, recipe=recipe, dr=1.0_wp, becke_k=1)
      call new(atom, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      if (.not. allocated(err)) call raw%update(atom, err)
      call require_success(error, err)
      if (allocated(error)) return
      call new_molecular_grid(domain, err, recipe=recipe, dr=1.0_wp, becke_k=1, pruning_threshold=0.1_wp)
      call require_success(error, err)
      if (allocated(error)) return
      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
      allocate (before(2*raw%ngrid), expected(2*raw%ngrid), points(3, raw%ngrid), partition(raw%ngrid))
      previous_count = 0
      do step = 1, 2
         call domain%update(mol, err)
         call require_success(error, err)
         if (allocated(error)) return
         j = 0
         expected = 0.0_wp
         do site = 1, 2
            do i = 1, raw%ngrid
               points(:, i) = raw%xyz(:, i) + mol%xyz(:, site)
            end do
            call becke_partition_weights(site, points, mol%xyz, [1, 1], partition, &
               & stiffness=domain%get_becke_k())
            do i = 1, raw%ngrid
               point = points(:, i)
               weight = raw%w(i)*partition(i)
               expected((site - 1)*raw%ngrid + i) = weight
               if (partition(i) < domain%get_pruning_threshold() .or. abs(weight) < 1.0e-14_wp) cycle
               j = j + 1
               call check(error, j <= domain%ngrid, "domain omitted a retained quadrature point")
               if (allocated(error)) return
               call check(error, domain%owner(j) == site &
                  & .and. maxval(abs(domain%xyz(:, j) - point)) < 1.0e-13_wp &
                  & .and. abs(domain%w(j) - weight) < 1.0e-13_wp, &
                  & "retained local point or updated Becke weight differs from reference")
               if (allocated(error)) return
            end do
         end do
         call check(error, j, domain%ngrid, "domain retained a point below the threshold")
         if (allocated(error)) return
         if (step == 1) then
            before = expected
            previous_count = domain%ngrid
            modes = [domain%nkx, domain%nky, domain%nkz]
            mol%xyz(1, 2) = mol%xyz(1, 2) + 0.4_wp
         else
            call check(error, maxval(abs(expected - before)) > 1.0e-4_wp, "deformation must change Becke weights")
            if (allocated(error)) return
            call check(error, previous_count /= domain%ngrid, "deformation must cross the pruning threshold")
            if (allocated(error)) return
            call check(error, all([domain%nkx, domain%nky, domain%nkz] == modes))
         end if
      end do
   end subroutine test_molecular_pruning

   !> Update and rebuild consume the single set of construction fields
   !>
   !> - Cartesian: rebuild picks up new nx, dr and xi0_factor, centered on the
   !>   molecule, settings survive destroy and update
   !> - Molecular: quadrature settings are fixed at construction, so rebuild
   !>   keeps them and a new grid picks up nrad, gaussian and xi0_factor;
   !>   `molecular_grid_set_kgrid` replaces the reciprocal settings, rebuild
   !>   re-sizes with them and keeps nufft_tol
   subroutine test_grid_settings(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_cartesian_type) :: cart
      type(moist_math_grid_3d_molecular_type) :: molecular, refined
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol
      integer :: count, coarse(3), modes(3)

      call new(mol, [1], reshape([0.2_wp, -0.4_wp, 0.6_wp], [3, 1]))
      cart%nx = 4
      cart%ny = 6
      cart%nz = 8
      cart%dr = 0.5_wp
      call cart%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      cart%nx = 10
      cart%dr = 0.3_wp
      cart%xi0_factor = 1.2_wp
      call cart%rebuild(err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, cart%ngrid, 10*6*8)
      if (allocated(error)) return
      call check(error, maxval(abs(cart%w - 0.3_wp**3)) < 1.0e-14_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(cart%xi0 - 4.0_wp)) < 1.0e-14_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(sum(cart%xyz, dim=2)/real(cart%ngrid, wp) - mol%xyz(:, 1))) < 1.0e-13_wp)
      if (allocated(error)) return
      call cart%destroy()
      call cart%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, cart%ngrid, 10*6*8)
      if (allocated(error)) return
      call check(error, maxval(abs(cart%xi0 - 4.0_wp)) < 1.0e-14_wp)
      if (allocated(error)) return
      call get_qc_handymod_recipe(recipe, err, nrad=8, degree=5, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(molecular, err, recipe=recipe, dr=1.0_wp, &
         & nufft_tol=1.0e-8_wp)
      if (.not. allocated(err)) call molecular%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      count = molecular%ngrid
      coarse = [molecular%nkx, molecular%nky, molecular%nkz]
      call molecular%rebuild(err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, molecular%ngrid, count, "rebuild must keep the constructed quadrature")
      if (allocated(error)) return
      call check(error, .not. allocated(molecular%xi0), "rebuild must keep point potentials")
      if (allocated(error)) return
      call get_qc_handymod_recipe(recipe, err, nrad=16, degree=5, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(refined, err, recipe=recipe, dr=1.0_wp, &
         & nufft_tol=1.0e-8_wp, gaussian=.true., xi0_factor=1.3_wp)
      if (.not. allocated(err)) call refined%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, refined%ngrid, 2*count)
      if (allocated(error)) return
      call check(error, maxval(abs(refined%xi0**3*refined%w - 1.3_wp**3)) < 1.0e-13_wp)
      if (allocated(error)) return
      call molecular_grid_set_kgrid(molecular, 0.5_wp, err, buffer=1.5_wp, reference=mol%xyz(:, 1))
      call require_success(error, err)
      if (allocated(error)) return
      modes = [molecular%nkx, molecular%nky, molecular%nkz]
      call check(error, all(modes > coarse), "a finer dr must enlarge the reciprocal grid")
      if (allocated(error)) return
      call molecular%rebuild(err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, all([molecular%nkx, molecular%nky, molecular%nkz] == modes), &
         & "rebuild must re-size with the stored reciprocal settings")
      if (allocated(error)) return
      call check(error, molecular%get_dr(), 0.5_wp, thr=epsilon(1.0_wp))
      if (allocated(error)) return
      call check(error, molecular%get_kbuffer(), 1.5_wp, thr=epsilon(1.0_wp))
      if (allocated(error)) return
      call check(error, molecular%get_nufft_tol(), 1.0e-8_wp, thr=epsilon(1.0_wp))
   end subroutine test_grid_settings

   !> Molecular grid stays usable after a rejected update
   !>
   !> - An oxygen override whose radial cutoff removes every node fails the
   !>   first update that meets oxygen, inside the assembly
   !> - The rejected update leaves the realized points unchanged
   !> - The same grid then updates cleanly on a geometry without oxygen
   subroutine test_update_recovers(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: domain
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type) :: overrides(1)
      type(structure_type) :: mol, water
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: shift(3)

      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call new(water, [8, 1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.43_wp, -0.98_wp, &
         & 0.0_wp, -1.43_wp, -0.98_wp], [3, 3]))
      call get_handymod_recipe(recipe, 8, 26, 0.0_wp, 4.0_wp, 2.0_wp, err)
      if (.not. allocated(err)) call get_handymod_recipe(overrides(1)%recipe, 8, 26, 0.0_wp, 4.0_wp, &
         & 2.0_wp, err)
      call require_success(error, err)
      if (allocated(error)) return
      ! Every HandyMod node lies below rmax = 4, so this cutoff empties the oxygen grid
      overrides(1)%elements = [8]
      overrides(1)%recipe%radial%rcut_lower = 5.0_wp
      call new_molecular_grid(domain, err, recipe=recipe, overrides=overrides)
      if (.not. allocated(err)) call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      xyz = domain%xyz
      call domain%update(water, err)
      call check(error, allocated(err), "an oxygen grid without radial nodes must fail the update")
      if (allocated(error)) return
      call check(error, index(err%message, "cutoffs removed every node") > 0, err%message)
      if (allocated(error)) return
      call check(error, domain%ngrid == size(xyz, 2) .and. all(domain%xyz == xyz), &
         & "a rejected update must leave the realized points unchanged")
      if (allocated(error)) return
      shift = [0.1_wp, -0.2_wp, 0.3_wp]
      mol%xyz(:, 1) = shift
      call domain%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, domain%ngrid == size(xyz, 2), "recovered grid lost points")
      if (allocated(error)) return
      call check(error, maxval(abs(domain%xyz - xyz - spread(shift, 2, domain%ngrid))) < 2.0e-15_wp, &
         & "recovered grid does not follow its atom")
   end subroutine test_update_recovers

   !> Single-atom structure at the origin
   !>
   !> @param[out] mol  structure
   subroutine single_atom(mol)
      !> Structure
      type(structure_type), intent(out) :: mol

      call new(mol, [1], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
   end subroutine single_atom

   !> add_weights on an uninitialized volume adjoint fails
   subroutine test_volume_adjoint_uninit_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%add_weights(err, w_w=[1.0_wp])
      call expect_error(error, err, "not initialized")
   end subroutine test_volume_adjoint_uninit_fails

   !> Position weights of the wrong shape fail
   subroutine test_volume_adjoint_xyz_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: err
      real(wp) :: w_xyz(2, 4)

      call acc%init(4)
      w_xyz = 1.0_wp
      call acc%add_weights(err, w_xyz=w_xyz)
      call expect_error(error, err, "xyz weight shape mismatch")
   end subroutine test_volume_adjoint_xyz_fails

   !> Integration-weight weights of the wrong length fail
   subroutine test_volume_adjoint_w_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%init(4)
      call acc%add_weights(err, w_w=[1.0_wp, 2.0_wp, 3.0_wp])
      call expect_error(error, err, "integration weight size mismatch")
   end subroutine test_volume_adjoint_w_fails

   !> Gaussian-width weights of the wrong length fail
   subroutine test_volume_adjoint_xi_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: err

      call acc%init(4)
      call acc%add_weights(err, w_xi=[1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp])
      call expect_error(error, err, "xi weight size mismatch")
   end subroutine test_volume_adjoint_xi_fails

   !> Cartesian grid with a zero extent fails
   subroutine test_cartesian_size_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err

      call new_cartesian_grid_3d(grid, 0, 4, 4, 0.5_wp, error=err)
      call expect_error(error, err, "must be positive")
   end subroutine test_cartesian_size_fails

   !> Cartesian point count beyond the default integer range fails
   subroutine test_cartesian_count_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err

      grid%nx = huge(0)
      call grid%validate(err)
      call expect_error(error, err, "exceeds the integer range")
   end subroutine test_cartesian_count_fails

   !> Negative Cartesian spacing fails
   subroutine test_cartesian_dr_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err

      call new_cartesian_grid_3d(grid, 4, 4, 4, -0.5_wp, error=err)
      call expect_error(error, err, "dr must be finite and positive")
   end subroutine test_cartesian_dr_fails

   !> NaN Cartesian spacing fails
   subroutine test_cartesian_dr_nan_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err

      grid%dr = ieee_value(1.0_wp, ieee_quiet_nan)
      call grid%validate(err)
      call expect_error(error, err, "dr and xi0_factor must be finite")
   end subroutine test_cartesian_dr_nan_fails

   !> Non-positive Cartesian Gaussian-width factor fails
   subroutine test_cartesian_xi0_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err

      grid%xi0_factor = 0.0_wp
      call grid%validate(err)
      call expect_error(error, err, "xi0_factor must be finite and positive")
   end subroutine test_cartesian_xi0_fails

   !> Non-finite Cartesian origin fails
   subroutine test_cartesian_origin_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err
      real(wp) :: origin(3)

      origin = [0.0_wp, ieee_value(1.0_wp, ieee_quiet_nan), 0.0_wp]
      call new_cartesian_grid_3d(grid, 4, 4, 4, 0.5_wp, origin, err)
      call expect_error(error, err, "origin must be finite")
   end subroutine test_cartesian_origin_fails

   !> Cartesian update on a structure without atoms fails
   subroutine test_cartesian_no_atoms_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err
      type(structure_type) :: mol

      call grid%update(mol, err)
      call expect_error(error, err, "at least one solute atom")
   end subroutine test_cartesian_no_atoms_fails

   !> Cartesian update on non-finite coordinates fails
   subroutine test_cartesian_coords_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err
      type(structure_type) :: mol

      call single_atom(mol)
      mol%xyz(2, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call grid%update(mol, err)
      call expect_error(error, err, "coordinates must be finite")
   end subroutine test_cartesian_coords_fails

   !> Cartesian rebuild before any update fails
   subroutine test_cartesian_rebuild_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err

      call grid%rebuild(err)
      call expect_error(error, err, "before rebuild")
   end subroutine test_cartesian_rebuild_fails

   !> Cartesian transform engine before any update fails
   subroutine test_cartesian_trafo_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: err

      call grid%new_trafo(trafo, err)
      call expect_error(error, err, "before new_trafo")
   end subroutine test_cartesian_trafo_fails

   !> Zero Becke-SSF parameter fails molecular construction
   subroutine test_molecular_ssf_a_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: err

      call new_molecular_grid(grid, err, ssf_a=0.0_wp)
      call expect_error(error, err, "ssf_a must be in (0, 1]")
   end subroutine test_molecular_ssf_a_fails

   !> Pruning threshold of one fails molecular construction
   subroutine test_molecular_pruning_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: err

      call new_molecular_grid(grid, err, pruning_threshold=1.0_wp)
      call expect_error(error, err, "pruning_threshold must be in [0, 1)")
   end subroutine test_molecular_pruning_fails

   !> Molecular rebuild of a constructed grid before any update fails
   subroutine test_molecular_rebuild_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(mctc_error), allocatable :: err

      call new_molecular_grid(grid, err)
      if (allocated(err)) return
      call grid%rebuild(err)
      call expect_error(error, err, "update before rebuild")
   end subroutine test_molecular_rebuild_fails

   !> Molecular update on non-finite coordinates fails
   subroutine test_molecular_coords_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(mctc_error), allocatable :: err
      type(structure_type) :: mol

      call single_atom(mol)
      mol%xyz(3, 1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call get_qc_handymod_recipe(recipe, err, nrad=4, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(grid, err, recipe=recipe)
      if (allocated(err)) return
      call grid%update(mol, err)
      call expect_error(error, err, "coordinates must be finite")
   end subroutine test_molecular_coords_fails

   !> Reciprocal mode count overflow fails before the integer conversion
   subroutine test_molecular_kcount_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(mctc_error), allocatable :: err
      type(structure_type) :: mol

      call single_atom(mol)
      call get_qc_handymod_recipe(recipe, err, nrad=4)
      if (.not. allocated(err)) call new_molecular_grid(grid, err, recipe=recipe, dr=1.0e-12_wp)
      if (allocated(err)) return
      call grid%update(mol, err)
      call expect_error(error, err, "reciprocal grid is too large")
   end subroutine test_molecular_kcount_fails

   !> Per-atom point count beyond the default integer range fails
   !>
   !> - Fails before the raw-buffer allocation is attempted
   subroutine test_point_count_overflow_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(mctc_error), allocatable :: err
      type(structure_type) :: mol

      call single_atom(mol)
      call get_uniform_recipe(recipe, overrides, 1000000, 5810, err, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(grid, err, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (allocated(err)) return
      call grid%update(mol, err)
      call expect_error(error, err, "point count exceeds")
   end subroutine test_point_count_overflow_fails

   !> Non-positive HandyMod parameter fails the recipe
   subroutine test_handymod_nonpositive_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(mctc_error), allocatable :: err

      call get_handymod_recipe(recipe, 8, 26, 0.0_wp, 4.0_wp, 0.0_wp, err)
      call expect_error(error, err, "HandyMod mapping: m must be > 0")
   end subroutine test_handymod_nonpositive_fails

   !> HandyMod parameter too large for the radial interval fails the recipe
   !>
   !> - 2**m - 1 must stay below rmax - rmin; m = 1000 keeps 2**m finite
   subroutine test_handymod_huge_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(mctc_error), allocatable :: err

      call get_handymod_recipe(recipe, 8, 26, 0.0_wp, 4.0_wp, 1000.0_wp, err)
      call expect_error(error, err, "HandyMod mapping: rmax - rmin must be greater than 2^m - 1")
   end subroutine test_handymod_huge_fails

   !> Grid constructed without a reciprocal grid keeps none across update and rebuild
   subroutine test_rebuild_no_kgrid(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(structure_type) :: mol

      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call get_uniform_recipe(recipe, overrides, 8, 26, err, rmin=0.0_wp, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(grid, err, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(err)) call grid%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, .not. grid%has_kgrid .and. grid%npts_k == 0, &
         & "a grid constructed with reciprocal = .false. must not create a k-grid on its own")
      if (allocated(error)) return
      mol%xyz(1, 2) = mol%xyz(1, 2) + 0.2_wp
      call grid%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, .not. grid%has_kgrid, "update of a k-gridless grid must not create one")
      if (allocated(error)) return
      call grid%rebuild(err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, .not. grid%has_kgrid .and. grid%npts_k == 0, &
         & "rebuild of a grid without a k-grid must only rebuild geometry")
   end subroutine test_rebuild_no_kgrid

   !> destroy() resets the explicit reference offset, so a later update matches a fresh grid
   subroutine test_destroy_kref(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol
      real(wp) :: centroid(3), offset(3)

      call get_qc_handymod_recipe(recipe, err, nrad=8, rmax=4.0_wp)
      if (.not. allocated(err)) call new_molecular_grid(grid, err, recipe=recipe, dr=1.0_wp)
      call require_success(error, err)
      if (allocated(error)) return
      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call grid%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      centroid = sum(mol%xyz, dim=2)/real(mol%nat, wp)
      offset = [0.5_wp, 0.2_wp, 0.0_wp]
      call molecular_grid_set_kgrid(grid, grid%get_dr(), err, reference=centroid + offset)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(grid%kref - centroid - offset)) < 1.0e-13_wp)
      if (allocated(error)) return
      call grid%destroy()
      call grid%update(mol, err)
      call require_success(error, err)
      if (allocated(error)) return
      call check(error, maxval(abs(grid%kref - centroid)) < 1.0e-13_wp, &
         & "destroy must reset the explicit reference offset before the next update")
   end subroutine test_destroy_kref

end module test_math_grid_3d
