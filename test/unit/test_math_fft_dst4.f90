!> Test suite for the DST-IV behind the 1D-RISM Fourier-Bessel transform
!>
!> - Pins the ducc0 type-4 DST (`ortho` off) to FFTW's `RODFT11` convention,
!>   which the radial grid prefactors assume
!> - Backend: DST-IV against the defining sum with an exactly reduced phase
!> - Identities: batched call matches single column, DST-IV is its own inverse
!>   up to 2n
!> - Radial grid: Fourier-Bessel round trip is exact, adjoint identity, and a
!>   Gaussian with a closed-form transform pins the absolute prefactors
!> - Error paths (`should_fail`): uniform radial grid rejects zero points and
!>   zero, NaN or infinite spacing
module test_math_fft_dst4
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, &
      & ieee_positive_inf
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_fft_dst4, only: dst4_plan_type, dst4_work_type
   use moist_math_grid_1d_base, only: moist_math_grid_1d_trafo_type
   use moist_math_grid_1d_uniform, only: moist_math_grid_1d_uniform_type, &
      & new_uniform_radial_grid
   implicit none(type, external)
   private

   public :: collect_math_fft_dst4

   !> Transform lengths exercised by the kernel tests: powers of two, an odd length
   !>
   !> A prime is also included -- a backend that only handled smooth sizes, or
   !> that fell back to a different algorithm for awkward ones, would show up
   !> here
   integer, parameter :: test_lengths(7) = [8, 9, 16, 31, 64, 67, 69]

contains

   !> Collect all math_fft_dst4 tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_fft_dst4(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dst4_matches_reference", test_dst4_reference_match), &
                  new_unittest("dst4_batched_matches_single", test_dst4_batched), &
                  new_unittest("dst4_is_its_own_inverse", test_dst4_involution), &
                  new_unittest("bad_uniform_npts_zero", test_bad_uniform_npts_zero, should_fail=.true.), &
                  new_unittest("bad_uniform_dr_zero", test_bad_uniform_dr_zero, should_fail=.true.), &
                  new_unittest("bad_uniform_dr_nan", test_bad_uniform_dr_nan, should_fail=.true.), &
                  new_unittest("bad_uniform_dr_inf", test_bad_uniform_dr_inf, should_fail=.true.), &
                  new_unittest("uniform_fbt_round_trip", test_fbt_round_trip), &
                  new_unittest("uniform_fbt_batched_matches_single", test_fbt_batched), &
                  new_unittest("uniform_fbt_batched_empty_batch", test_fbt_batched_empty), &
                  new_unittest("uniform_fbt_adjoint", test_fbt_adjoint), &
                  new_unittest("uniform_fbt_matches_analytic_gaussian", test_fbt_analytic) &
                  ]
   end subroutine collect_math_fft_dst4

   ! --------------------------------------------------------------------------
   ! Helpers
   ! --------------------------------------------------------------------------

   !> Deterministic, well-conditioned test signal (no RNG state to seed)
   !>
   !> @param[in] i  1-based sample index
   pure function signal(i) result(x)
      !> Sample index
      integer, intent(in) :: i
      !> Sample value
      real(wp) :: x

      x = sin(1.7_wp*real(i, wp)) + 0.3_wp*cos(0.31_wp*real(i*i, wp)) &
         & + 0.1_wp*real(modulo(i, 7) - 3, wp)
   end function signal

   !> Unnormalised DST-IV evaluated directly from its definition
   !>
   !> `Y_k = 2 sum_j X_j sin(pi (j+1/2)(k+1/2)/n)`; the phase is
   !> `pi*(2j+1)(2k+1)/(4n)`, so reducing the integer product modulo `8n`
   !> reduces the argument modulo `2*pi` exactly; without that, arguments of
   !> order `pi*n` would lose several digits in `sin` and the reference would
   !> be less accurate than the transform it checks
   !>
   !> @param[in]  x  Input samples (size n)
   !> @param[out] y  Transformed samples (size n)
   subroutine dst4_reference(x, y)
      !> Input samples
      real(wp), intent(in) :: x(:)
      !> Transformed samples
      real(wp), intent(out) :: y(:)

      integer :: n, j, k, m
      real(wp) :: pi, acc

      n = size(x)
      pi = acos(-1.0_wp)
      do k = 1, n
         acc = 0.0_wp
         do j = 1, n
            m = modulo((2*(j - 1) + 1)*(2*(k - 1) + 1), 8*n)
            acc = acc + x(j)*sin(pi*real(m, wp)/(4.0_wp*real(n, wp)))
         end do
         y(k) = 2.0_wp*acc
      end do
   end subroutine dst4_reference

   !> Largest elementwise deviation, relative to the scale of `ref`
   !>
   !> A non-finite entry in either array returns `huge(1.0_wp)` rather than
   !> feeding NaN into `maxval`/`max`: gfortran's pairwise reduction can
   !> silently drop a NaN comparison and keep the running (finite) extremum,
   !> which would let a corrupted value pass every `check` built on this
   !> function
   !>
   !> @param[in] got  Values under test
   !> @param[in] ref  Reference values
   pure function relerr(got, ref) result(err)
      !> Values under test
      real(wp), intent(in) :: got(:)
      !> Reference values
      real(wp), intent(in) :: ref(:)
      !> Relative max-norm deviation
      real(wp) :: err

      real(wp) :: scal

      if (.not. all(ieee_is_finite(got)) .or. .not. all(ieee_is_finite(ref))) then
         err = huge(1.0_wp)
         return
      end if
      scal = maxval(abs(ref))
      if (scal <= 0.0_wp) scal = 1.0_wp
      err = maxval(abs(got - ref))/scal
   end function relerr

   !> Run one single-column DST-IV through a freshly built plan
   !>
   !> @param[out] error  Set if the DST-IV backend failed
   !> @param[in]  x      Input samples (size n)
   !> @param[out] y      Transformed samples (size n)
   subroutine dst4_once(error, x, y)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Input samples
      real(wp), intent(in) :: x(:)
      !> Transformed samples
      real(wp), intent(out) :: y(:)

      type(dst4_plan_type) :: plan
      type(dst4_work_type) :: work
      real(wp), allocatable :: buf(:)
      type(mctc_error), allocatable :: merr

      call plan%init(size(x), 1)
      call plan%new_work(work)
      buf = x
      call plan%execute(work, buf, y, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call work%destroy()
      call plan%destroy()
   end subroutine dst4_once

   !> Build a uniform radial grid and its trafo
   !>
   !> @param[out] error  Set on failure
   !> @param[out] grid   Initialised grid (must be a target)
   !> @param[out] trafo  Trafo bound to `grid`
   !> @param[in]  npts   Node count
   !> @param[in]  dr     Node spacing (bohr)
   subroutine make_grid(error, grid, trafo, npts, dr)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Initialised grid
      type(moist_math_grid_1d_uniform_type), intent(out), target :: grid
      !> Trafo bound to `grid`
      class(moist_math_grid_1d_trafo_type), allocatable, intent(out) :: trafo
      !> Node count
      integer, intent(in) :: npts
      !> Node spacing
      real(wp), intent(in) :: dr

      type(mctc_error), allocatable :: merr

      call new_uniform_radial_grid(grid, npts, dr, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call grid%new_trafo(trafo, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
   end subroutine make_grid

   ! --------------------------------------------------------------------------
   ! Layer 1: the backend against the defining sum
   ! --------------------------------------------------------------------------

   !> The backend's DST-IV must be the unnormalised RODFT11 convention
   !>
   !> The radial prefactors assume this, for even, odd and prime lengths alike
   !>
   !> @param[out] error  Test failure
   subroutine test_dst4_reference_match(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: it, n, i
      real(wp), allocatable :: x(:), y(:), ref(:)

      do it = 1, size(test_lengths)
         n = test_lengths(it)
         allocate (x(n), y(n), ref(n))
         do i = 1, n
            x(i) = signal(i)
         end do
         call dst4_once(error, x, y)
         if (allocated(error)) return
         call dst4_reference(x, ref)
         call check(error, relerr(y, ref), 0.0_wp, thr=1.0e-14_wp)
         if (allocated(error)) return
         deallocate (x, y, ref)
      end do
   end subroutine test_dst4_reference_match

   ! --------------------------------------------------------------------------
   ! Layer 2: structural identities
   ! --------------------------------------------------------------------------

   !> The batched call the radial trafo uses must reproduce the single-column one
   !>
   !> A mismatch would silently change results whenever a caller switches to
   !> the batched entry points
   !>
   !> @param[out] error  Test failure
   subroutine test_dst4_batched(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 48, nb = 5
      type(dst4_plan_type) :: plan
      type(dst4_work_type) :: work
      type(mctc_error), allocatable :: merr
      integer :: ic, i
      real(wp) :: xall(n, nb), xin(n, nb), yall(n, nb), col(n), yref(n)

      do ic = 1, nb
         do i = 1, n
            xall(i, ic) = signal(i + 13*ic)
         end do
      end do

      call plan%init(n, nb)
      call plan%new_work(work)
      xin = xall
      call plan%execute(work, xin, yall, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call work%destroy()
      call plan%destroy()

      do ic = 1, nb
         col = xall(:, ic)
         call dst4_reference(col, yref)
         call check(error, relerr(yall(:, ic), yref), 0.0_wp, thr=1.0e-14_wp)
         if (allocated(error)) return
      end do
   end subroutine test_dst4_batched

   !> DST-IV is its own inverse up to the factor 2n
   !>
   !> The radial grid relies on exactly this when it uses one transform for
   !> both directions
   !>
   !> @param[out] error  Test failure
   subroutine test_dst4_involution(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 64
      integer :: i
      real(wp) :: x(n), y(n), z(n), scaled(n)

      do i = 1, n
         x(i) = signal(i)
      end do

      call dst4_once(error, x, y)
      if (allocated(error)) return
      call dst4_once(error, y, z)
      if (allocated(error)) return
      scaled = 2.0_wp*real(n, wp)*x
      call check(error, relerr(z, scaled), 0.0_wp, thr=1.0e-14_wp)
   end subroutine test_dst4_involution

   ! --------------------------------------------------------------------------
   ! Uniform radial grid constructor: invalid parameters
   ! --------------------------------------------------------------------------

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
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

   !> A grid without nodes is refused
   !>
   !> @param[out] error  Test failure, set on the expected error
   subroutine test_bad_uniform_npts_zero(error)
      !> Test failure, set on the expected error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_uniform_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_grid(grid, 0, 0.1_wp, merr)
      call expect_error(error, merr, "at least one point")
   end subroutine test_bad_uniform_npts_zero

   !> A zero spacing is refused
   !>
   !> @param[out] error  Test failure, set on the expected error
   subroutine test_bad_uniform_dr_zero(error)
      !> Test failure, set on the expected error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_uniform_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_grid(grid, 16, 0.0_wp, merr)
      call expect_error(error, merr, "positive spacing")
   end subroutine test_bad_uniform_dr_zero

   !> A NaN spacing is refused rather than building NaN nodes
   !>
   !> `dr <= 0` is false for NaN, so only an explicit finiteness guard catches it
   !>
   !> @param[out] error  Test failure, set on the expected error
   subroutine test_bad_uniform_dr_nan(error)
      !> Test failure, set on the expected error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_uniform_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_grid(grid, 16, ieee_value(1.0_wp, ieee_quiet_nan), merr)
      call expect_error(error, merr, "finite spacing")
   end subroutine test_bad_uniform_dr_nan

   !> An infinite spacing is refused rather than building Inf nodes
   !>
   !> @param[out] error  Test failure, set on the expected error
   subroutine test_bad_uniform_dr_inf(error)
      !> Test failure, set on the expected error
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_1d_uniform_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_grid(grid, 16, ieee_value(1.0_wp, ieee_positive_inf), merr)
      call expect_error(error, merr, "finite spacing")
   end subroutine test_bad_uniform_dr_inf

   ! --------------------------------------------------------------------------
   ! Layer 3: the radial grid built on the DST-IV
   ! --------------------------------------------------------------------------

   !> `fbt_k2r . fbt_r2k` is the identity, not merely an approximation
   !>
   !> The weight diagonals telescope through the involution to
   !> (dk*dr*2n)/(2*pi), which is exactly 1 for dk = pi/(n*dr); this is a
   !> machine-precision check on the whole radial pipeline, and it fails if
   !> the backend's normalisation drifts by any factor at all
   !>
   !> @param[out] error  Test failure
   subroutine test_fbt_round_trip(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts = 256
      real(wp), parameter :: dr = 0.05_wp
      type(moist_math_grid_1d_uniform_type), target :: grid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      integer :: i
      real(wp) :: f(npts), g(npts), back(npts)

      call make_grid(error, grid, trafo, npts, dr)
      if (allocated(error)) return

      do i = 1, npts
         f(i) = exp(-0.7_wp*grid%r(i)**2)*(1.0_wp + 0.4_wp*grid%r(i))
      end do

      call trafo%fbt_r2k(f, g, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call trafo%fbt_k2r(g, back, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call check(error, relerr(back, f), 0.0_wp, thr=1.0e-13_wp)
      if (allocated(error)) return

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fbt_round_trip

   !> The batched `fbt_*_all` entry points must agree with the single-column ones
   !>
   !> The solvers use both
   !>
   !> @param[out] error  Test failure
   subroutine test_fbt_batched(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts = 128, nb = 4
      real(wp), parameter :: dr = 0.08_wp
      type(moist_math_grid_1d_uniform_type), target :: grid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      integer :: i, ic
      real(wp) :: f(npts, nb), gall(npts, nb), one(npts)

      call make_grid(error, grid, trafo, npts, dr)
      if (allocated(error)) return

      do ic = 1, nb
         do i = 1, npts
            f(i, ic) = exp(-(0.4_wp + 0.2_wp*ic)*grid%r(i)**2)
         end do
      end do

      call trafo%fbt_r2k_all(f, gall, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      do ic = 1, nb
         call trafo%fbt_r2k(f(:, ic), one, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message); return
         end if
         call check(error, relerr(gall(:, ic), one), 0.0_wp, thr=1.0e-14_wp)
         if (allocated(error)) return
      end do

      call trafo%fbt_k2r_all(f, gall, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      do ic = 1, nb
         call trafo%fbt_k2r(f(:, ic), one, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message); return
         end if
         call check(error, relerr(gall(:, ic), one), 0.0_wp, thr=1.0e-14_wp)
         if (allocated(error)) return
      end do

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fbt_batched

   !> A zero-width batch on a freshly built trafo must not crash
   !>
   !> The cached batched buffers start at width 0, so the very first
   !> `fbt_*_all` call with an empty batch used to match that cached width
   !> and skip allocating `work_all`/`tmp_all`, which were then handed
   !> unallocated to the non-allocatable `execute` dummies
   !>
   !> @param[out] error  Test failure
   subroutine test_fbt_batched_empty(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts = 32
      real(wp), parameter :: dr = 0.1_wp
      type(moist_math_grid_1d_uniform_type), target :: grid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp) :: f(npts, 0), g(npts, 0)

      call make_grid(error, grid, trafo, npts, dr)
      if (allocated(error)) return

      call trafo%fbt_r2k_all(f, g, merr)
      call check(error, .not. allocated(merr), "empty-batch forward transform reported an error")
      if (allocated(error)) return

      call trafo%fbt_k2r_all(f, g, merr)
      call check(error, .not. allocated(merr), "empty-batch backward transform reported an error")
      if (allocated(error)) return

      call trafo%fbt_r2k_adj_all(f, g, merr)
      call check(error, .not. allocated(merr), "empty-batch forward-adjoint transform reported an error")
      if (allocated(error)) return

      call trafo%fbt_k2r_adj_all(f, g, merr)
      call check(error, .not. allocated(merr), "empty-batch backward-adjoint transform reported an error")
      if (allocated(error)) return

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fbt_batched_empty

   !> `<FBT_rk a, b> = <a, FBT_rk^T b>`
   !>
   !> The RISM Newton/Krylov solvers apply the adjoint transforms and would
   !> drift silently if the two stopped being exact transposes of one another
   !>
   !> The residual is scaled by `||FBT a|| ||b||`, the natural size of the
   !> bilinear form; scaling by `|<FBT a, b>|` instead would measure the
   !> cancellation in the dot product (~2e-14) rather than how well the
   !> transpose holds, which is ~3e-16
   !>
   !> @param[out] error  Test failure
   subroutine test_fbt_adjoint(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts = 96
      real(wp), parameter :: dr = 0.1_wp
      type(moist_math_grid_1d_uniform_type), target :: grid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      integer :: i
      real(wp) :: a(npts), b(npts), fa(npts), atb(npts)
      real(wp) :: lhs, rhs, scal

      call make_grid(error, grid, trafo, npts, dr)
      if (allocated(error)) return

      do i = 1, npts
         a(i) = signal(i)*exp(-0.2_wp*grid%r(i))
         b(i) = signal(i + 41)*exp(-0.15_wp*grid%k(i))
      end do

      call trafo%fbt_r2k(a, fa, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call trafo%fbt_r2k_adj(b, atb, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      lhs = dot_product(fa, b)
      rhs = dot_product(a, atb)
      scal = max(norm2(fa)*norm2(b), 1.0e-30_wp)
      call check(error, abs(lhs - rhs)/scal, 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return

      call trafo%fbt_k2r(a, fa, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call trafo%fbt_k2r_adj(b, atb, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      lhs = dot_product(fa, b)
      rhs = dot_product(a, atb)
      scal = max(norm2(fa)*norm2(b), 1.0e-30_wp)
      call check(error, abs(lhs - rhs)/scal, 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fbt_adjoint

   !> Both transforms of a Gaussian match its closed-form 3D Fourier transform
   !>
   !> - `f(r) = exp(-a r^2)` has `F(k) = (pi/a)^(3/2) exp(-k^2/(4a))`;
   !>   `fbt_r2k(f)` is checked against `F`, `fbt_k2r(F)` against `f`
   !> - Pins the absolute prefactors: a factor `c` in `fbt_r2k` and `1/c` in
   !>   `fbt_k2r` passes the round trip and the adjoint test but fails here
   !> - Batched entry points on two widths, `a` and `1.3 a`
   !> - Deviation relative to the reference maximum; measured at most 4.2e-16
   !>   forward and 4.4e-16 backward (npts = 256, dr = 0.05, a = 0.7), the
   !>   Gaussian being resolved to round-off on both grids
   !>
   !> @param[out] error  Test failure
   subroutine test_fbt_analytic(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: npts = 256
      real(wp), parameter :: dr = 0.05_wp
      !> Gaussian exponents of the two columns (bohr^-2)
      real(wp), parameter :: expo(2) = [0.7_wp, 0.91_wp]
      !> Acceptance bound, relative to the reference maximum
      real(wp), parameter :: thr = 3.0e-15_wp
      real(wp), parameter :: pi = acos(-1.0_wp)
      type(moist_math_grid_1d_uniform_type), target :: grid
      class(moist_math_grid_1d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      integer :: ic
      real(wp) :: f(npts, 2), fk(npts, 2), got(npts, 2)

      call make_grid(error, grid, trafo, npts, dr)
      if (allocated(error)) return

      do ic = 1, 2
         f(:, ic) = exp(-expo(ic)*grid%r**2)
         fk(:, ic) = (pi/expo(ic))**1.5_wp*exp(-grid%k**2/(4.0_wp*expo(ic)))
      end do

      call trafo%fbt_r2k(f(:, 1), got(:, 1), merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call check(error, relerr(got(:, 1), fk(:, 1)), 0.0_wp, thr=thr)
      if (allocated(error)) return

      call trafo%fbt_k2r(fk(:, 1), got(:, 1), merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      call check(error, relerr(got(:, 1), f(:, 1)), 0.0_wp, thr=thr)
      if (allocated(error)) return

      call trafo%fbt_r2k_all(f, got, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      do ic = 1, 2
         call check(error, relerr(got(:, ic), fk(:, ic)), 0.0_wp, thr=thr)
         if (allocated(error)) return
      end do

      call trafo%fbt_k2r_all(fk, got, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message); return
      end if
      do ic = 1, 2
         call check(error, relerr(got(:, ic), f(:, ic)), 0.0_wp, thr=thr)
         if (allocated(error)) return
      end do

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fbt_analytic

end module test_math_fft_dst4
