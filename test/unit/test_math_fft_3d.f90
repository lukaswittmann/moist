!> Test suite for the ducc0-backed 3D FFT entry points in moist_math_fft
!>
!> Backend-level coverage, no grid object in between; a Fortran field `x(nx, ny, nz)`
!> is the C row-major `(n0, n1, n2) = (nz, ny, nx)`, so axis 0 is z (third Fortran
!> index) and axis 2 is x (first Fortran index)
!>
module test_math_fft_3d
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_fft, only: moist_fft_c2c_3d, moist_fft_c2c_3d_pass, &
      & moist_fft_c2c_3d_pass_inplace, moist_fft_r2c_3d, moist_fft_c2r_3d, &
      & moist_fft_r2c_3d_batch, moist_fft_c2r_3d_batch, moist_fft_dst4
   implicit none(type, external)
   private

   public :: collect_math_fft_3d

   !> Working precision, the C interface's double
   integer, parameter :: wp = c_double
   !> Complex kind of the C interface
   integer, parameter :: cwp = c_double_complex
   !> Circle constant
   real(wp), parameter :: pi = acos(-1.0_wp)

   !> Exponent of the analytic test Gaussian exp(-a |r - c|^2), 1/bohr^2
   real(wp), parameter :: gauss_a = 1.0_wp
   !> Off-center, off-axis Gaussian center c, bohr
   real(wp), parameter :: gauss_c(3) = [0.3_wp, -0.5_wp, 0.7_wp]

   !> Status returned by the backend for a ducc0 exception other than bad_alloc
   integer(c_int), parameter :: status_error = 1_c_int

contains

   !> Collect all math_fft_3d tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_fft_3d(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("c2c_plane_wave_spike", test_plane_wave_spike), &
                  new_unittest("c2c_matches_direct_dft", test_direct_dft), &
                  new_unittest("c2c_gaussian_analytic", test_gaussian), &
                  new_unittest("c2c_round_trip", test_round_trip), &
                  new_unittest("pass_round_trip", test_pass_round_trip), &
                  new_unittest("passes_compose_to_full", test_passes_compose), &
                  new_unittest("pass_inplace_matches_pass", test_inplace_matches_pass), &
                  new_unittest("pass_strided_subblock", test_strided_subblock), &
                  new_unittest("c2c_linearity", test_linearity), &
                  new_unittest("c2c_parseval", test_parseval), &
                  new_unittest("r2c_matches_c2c_half", test_r2c_half), &
                  new_unittest("c2r_matches_c2c_backward", test_c2r_backward), &
                  new_unittest("batch_matches_single", test_batch), &
                  new_unittest("zero_extent_is_noop", test_zero_extent), &
                  new_unittest("negative_extent_is_rejected", test_negative_extent), &
                  new_unittest("rejected_pass_leaves_data_untouched", test_rejected_pass), &
                  new_unittest("bad_pass_axis_too_large", test_pass_axis_large_fails, should_fail=.true.), &
                  new_unittest("bad_pass_axis_negative", test_pass_axis_negative_fails, should_fail=.true.), &
                  new_unittest("bad_pass_inplace_axis", test_pass_inplace_axis_fails, should_fail=.true.) &
                  ]
   end subroutine collect_math_fft_3d

   !* ======================================================================================== *!
   !*                                         Helpers                                          *!
   !* ======================================================================================== *!

   !> Deterministic complex test sample, no RNG state to seed
   !>
   !> @param[in] l  sample index
   pure function sample(l) result(z)
      !> Sample index
      integer, intent(in) :: l
      !> Sample value
      complex(wp) :: z

      real(wp) :: t

      t = real(l, wp)
      z = cmplx(sin(1.7_wp*t) + 0.1_wp*real(modulo(l, 7) - 3, wp), &
         & cos(0.31_wp*t*t) - 0.2_wp*sin(0.9_wp*t), wp)
   end function sample

   !> Fill a field with deterministic complex samples
   !>
   !> @param[out] x     field to fill
   !> @param[in]  seed  offset of the sample index
   subroutine fill(x, seed)
      !> Field to fill
      complex(wp), intent(out) :: x(:, :, :)
      !> Offset of the sample index
      integer, intent(in) :: seed

      integer :: i, j, k, l

      l = seed
      do k = 1, size(x, 3)
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               l = l + 1
               x(i, j, k) = sample(l)
            end do
         end do
      end do
   end subroutine fill

   !> Backend direction flag
   !>
   !> @param[in] fwd  forward transform, exp(-i k x)
   pure function dir(fwd) result(flag)
      !> Forward transform
      logical, intent(in) :: fwd
      !> Nonzero for forward
      integer(c_int) :: flag

      flag = merge(1_c_int, 0_c_int, fwd)
   end function dir

   !> Roots of unity exp(-+2 pi i m/n), m = 0..n-1 with n = size(tw)
   !>
   !> @param[in]  fwd  negative exponent (forward) if true
   !> @param[out] tw   roots of unity, zero-based
   pure subroutine twiddles(fwd, tw)
      !> Negative exponent if true
      logical, intent(in) :: fwd
      !> Roots of unity
      complex(wp), intent(out) :: tw(0:)

      integer :: m, n
      real(wp) :: sgn, ang

      n = size(tw)
      sgn = merge(-1.0_wp, 1.0_wp, fwd)
      do m = 0, n - 1
         ang = 2.0_wp*pi*real(m, wp)/real(n, wp)
         tw(m) = cmplx(cos(ang), sgn*sin(ang), wp)
      end do
   end subroutine twiddles

   !> Direct O(N^2) 3D DFT from its definition
   !>
   !> - y(k) = fct sum_j x(j) exp(-+2 pi i sum_d k_d j_d/n_d), k and j zero-based
   !> - Phase written as m/N with N = nx*ny*nz and m reduced modulo N in integer
   !>   arithmetic, so the reference is accurate to round-off
   !>
   !> @param[in]  x    input field (nx, ny, nz)
   !> @param[out] y    transformed field (nx, ny, nz)
   !> @param[in]  fwd  forward (negative exponent) if true
   !> @param[in]  fct  factor applied to every output value
   subroutine dft3(x, y, fwd, fct)
      !> Input field
      complex(wp), intent(in) :: x(:, :, :)
      !> Transformed field
      complex(wp), intent(out) :: y(:, :, :)
      !> Forward if true
      logical, intent(in) :: fwd
      !> Output factor
      real(wp), intent(in) :: fct

      complex(wp), allocatable :: tw(:)
      complex(wp) :: acc
      integer :: nx, ny, nz, ntot, kx, ky, kz, jx, jy, jz, m

      nx = size(x, 1); ny = size(x, 2); nz = size(x, 3)
      ntot = nx*ny*nz
      allocate (tw(0:ntot - 1))
      call twiddles(fwd, tw)
      do kz = 0, nz - 1
         do ky = 0, ny - 1
            do kx = 0, nx - 1
               acc = (0.0_wp, 0.0_wp)
               do jz = 0, nz - 1
                  do jy = 0, ny - 1
                     do jx = 0, nx - 1
                        m = modulo(modulo(kx*jx, nx)*ny*nz + modulo(ky*jy, ny)*nx*nz &
                           & + modulo(kz*jz, nz)*nx*ny, ntot)
                        acc = acc + x(jx + 1, jy + 1, jz + 1)*tw(m)
                     end do
                  end do
               end do
               y(kx + 1, ky + 1, kz + 1) = fct*acc
            end do
         end do
      end do
   end subroutine dft3

   !> Direct 1D DFT along one backend axis of a 3D block
   !>
   !> @param[in]  x     input block (mx, my, mz)
   !> @param[out] y     transformed block (mx, my, mz)
   !> @param[in]  axis  backend axis, 0 (z, dim 3), 1 (y, dim 2) or 2 (x, dim 1)
   !> @param[in]  fwd   forward (negative exponent) if true
   !> @param[in]  fct   factor applied to every output value
   subroutine dft_axis(x, y, axis, fwd, fct)
      !> Input block
      complex(wp), intent(in) :: x(:, :, :)
      !> Transformed block
      complex(wp), intent(out) :: y(:, :, :)
      !> Backend axis
      integer, intent(in) :: axis
      !> Forward if true
      logical, intent(in) :: fwd
      !> Output factor
      real(wp), intent(in) :: fct

      complex(wp), allocatable :: tw(:)
      complex(wp) :: acc
      integer :: dim, n, i, j, k, jj, idx(3), src(3)

      dim = 3 - axis
      n = size(x, dim)
      allocate (tw(0:n - 1))
      call twiddles(fwd, tw)
      do k = 1, size(x, 3)
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               idx = [i, j, k]
               src = idx
               acc = (0.0_wp, 0.0_wp)
               do jj = 1, n
                  src(dim) = jj
                  acc = acc + x(src(1), src(2), src(3))*tw(modulo((idx(dim) - 1)*(jj - 1), n))
               end do
               y(i, j, k) = fct*acc
            end do
         end do
      end do
   end subroutine dft_axis

   !> Largest elementwise deviation, relative to the largest reference magnitude
   !>
   !> - Non-finite entries return huge(1.0_wp), so a NaN cannot slip through a
   !>   reduction that drops it
   !>
   !> @param[in] got  values under test
   !> @param[in] ref  reference values
   pure function relerr(got, ref) result(err)
      !> Values under test
      complex(wp), intent(in) :: got(:)
      !> Reference values
      complex(wp), intent(in) :: ref(:)
      !> Relative max-norm deviation
      real(wp) :: err

      real(wp) :: scal

      if (.not. (all(ieee_is_finite(real(got))) .and. all(ieee_is_finite(aimag(got))) &
         & .and. all(ieee_is_finite(real(ref))) .and. all(ieee_is_finite(aimag(ref))))) then
         err = huge(1.0_wp)
         return
      end if
      scal = maxval(abs(ref))
      if (scal <= 0.0_wp) scal = 1.0_wp
      err = maxval(abs(got - ref))/scal
   end function relerr

   !> Real-valued counterpart of `relerr`
   !>
   !> @param[in] got  values under test
   !> @param[in] ref  reference values
   pure function relerr_real(got, ref) result(err)
      !> Values under test
      real(wp), intent(in) :: got(:)
      !> Reference values
      real(wp), intent(in) :: ref(:)
      !> Relative max-norm deviation
      real(wp) :: err

      err = relerr(cmplx(got, 0.0_wp, wp), cmplx(ref, 0.0_wp, wp))
   end function relerr_real

   !> Full 3D transform through `moist_fft_c2c_3d`
   !>
   !> @param[out] error  test failure on a nonzero status
   !> @param[in]  x      input field (nx, ny, nz)
   !> @param[out] y      transformed field (nx, ny, nz)
   !> @param[in]  fwd    forward if true
   !> @param[in]  fct    factor applied to every output value
   subroutine c2c(error, x, y, fwd, fct)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Input field
      complex(cwp), intent(in) :: x(:, :, :)
      !> Transformed field
      complex(cwp), intent(out) :: y(:, :, :)
      !> Forward if true
      logical, intent(in) :: fwd
      !> Output factor
      real(wp), intent(in) :: fct

      integer(c_int) :: status

      status = moist_fft_c2c_3d(int(size(x, 3), c_int), int(size(x, 2), c_int), &
         & int(size(x, 1), c_int), x, y, dir(fwd), fct)
      call check(error, status == 0, "moist_fft_c2c_3d returned a nonzero status")
   end subroutine c2c

   !> One out-of-place single-axis pass over a whole contiguous field
   !>
   !> @param[out]    error  test failure on a nonzero status
   !> @param[in]     x      input field (nx, ny, nz)
   !> @param[in,out] y      transformed field (nx, ny, nz)
   !> @param[in]     axis   backend axis, 0 (z), 1 (y) or 2 (x)
   !> @param[in]     fwd    forward if true
   !> @param[in]     fct    factor applied to every transformed value
   subroutine pass(error, x, y, axis, fwd, fct)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Input field
      complex(cwp), intent(in) :: x(:, :, :)
      !> Transformed field
      complex(cwp), intent(inout) :: y(:, :, :)
      !> Backend axis
      integer, intent(in) :: axis
      !> Forward if true
      logical, intent(in) :: fwd
      !> Output factor
      real(wp), intent(in) :: fct

      integer(c_int) :: status

      status = moist_fft_c2c_3d_pass(int(size(x, 3), c_int), int(size(x, 2), c_int), &
         & int(size(x, 1), c_int), int(size(x, 1)*size(x, 2), c_int), int(size(x, 1), c_int), &
         & x, y, int(axis, c_int), dir(fwd), fct)
      call check(error, status == 0, "moist_fft_c2c_3d_pass returned a nonzero status")
   end subroutine pass

   !> One in-place single-axis pass over a whole contiguous field
   !>
   !> @param[out]    error  test failure on a nonzero status
   !> @param[in,out] x      field (nx, ny, nz), transformed in place
   !> @param[in]     axis   backend axis, 0 (z), 1 (y) or 2 (x)
   !> @param[in]     fwd    forward if true
   !> @param[in]     fct    factor applied to every transformed value
   subroutine pass_inplace(error, x, axis, fwd, fct)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Field, transformed in place
      complex(cwp), intent(inout) :: x(:, :, :)
      !> Backend axis
      integer, intent(in) :: axis
      !> Forward if true
      logical, intent(in) :: fwd
      !> Output factor
      real(wp), intent(in) :: fct

      integer(c_int) :: status

      status = moist_fft_c2c_3d_pass_inplace(int(size(x, 3), c_int), int(size(x, 2), c_int), &
         & int(size(x, 1), c_int), int(size(x, 1)*size(x, 2), c_int), int(size(x, 1), c_int), &
         & x, int(axis, c_int), dir(fwd), fct)
      call check(error, status == 0, "moist_fft_c2c_3d_pass_inplace returned a nonzero status")
   end subroutine pass_inplace

   !* ======================================================================================== *!
   !*                                  Discrete analytic                                       *!
   !* ======================================================================================== *!

   !> Single plane wave maps to a single spike of height fct*N
   !>
   !> - Grid 5 x 6 x 7 (all extents distinct, so a transposed layout misplaces
   !>   the spike), wave vectors with negative components, zero and the ky = 3
   !>   Nyquist mode of the even axis
   !> - Forward of exp(+2 pi i k.n/N) and backward of exp(-2 pi i k.n/N) both put
   !>   fct*N at zero-based index k mod n on each axis; a flipped sign lands on
   !>   -k mod n instead
   !> - Every other bin is zero to round-off
   !>
   !> @param[out] error  test failure
   subroutine test_plane_wave_spike(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 5, ny = 6, nz = 7, ncase = 5
      integer, parameter :: kvec(3, ncase) = reshape([0, 0, 0, 2, -1, 3, -2, 3, -3, 1, 0, -1, &
         & -1, -2, 2], [3, ncase])
      ! Measured 1.4e-16 relative to fct*N
      real(wp), parameter :: tol = 1.0e-15_wp
      real(wp), parameter :: fct = 0.625_wp
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :), ref(:, :, :)
      complex(wp), allocatable :: tw(:)
      integer :: ntot, ic, idir, jx, jy, jz, m
      logical :: fwd

      ntot = nx*ny*nz
      allocate (x(nx, ny, nz), y(nx, ny, nz), ref(nx, ny, nz), tw(0:ntot - 1))
      do idir = 1, 2
         fwd = idir == 1
         ! Input carries the opposite sign of the transform kernel
         call twiddles(.not. fwd, tw)
         do ic = 1, ncase
            do jz = 0, nz - 1
               do jy = 0, ny - 1
                  do jx = 0, nx - 1
                     m = modulo(modulo(kvec(1, ic)*jx, nx)*ny*nz + modulo(kvec(2, ic)*jy, ny)*nx*nz &
                        & + modulo(kvec(3, ic)*jz, nz)*nx*ny, ntot)
                     x(jx + 1, jy + 1, jz + 1) = tw(m)
                  end do
               end do
            end do
            ref = (0.0_wp, 0.0_wp)
            ref(modulo(kvec(1, ic), nx) + 1, modulo(kvec(2, ic), ny) + 1, modulo(kvec(3, ic), nz) + 1) = &
               & cmplx(fct*real(ntot, wp), 0.0_wp, wp)

            call c2c(error, x, y, fwd, fct)
            if (allocated(error)) return
            call check(error, relerr([y], [ref]) < tol, &
               & "plane wave does not map to a single spike of height fct*N at its wave vector")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_plane_wave_spike

   !> c2c_3d matches a direct O(N^2) DFT in both directions
   !>
   !> - 5 x 6 x 7 grid: odd and even, non-power-of-two, all extents distinct
   !> - fct = 1 and fct = 0.37 in each direction
   !>
   !> @param[out] error  test failure
   subroutine test_direct_dft(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 5, ny = 6, nz = 7
      ! Measured 6.3e-16
      real(wp), parameter :: tol = 5.0e-15_wp
      real(wp), parameter :: fcts(2) = [1.0_wp, 0.37_wp]
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :), ref(:, :, :)
      integer :: idir, ifct
      logical :: fwd

      allocate (x(nx, ny, nz), y(nx, ny, nz), ref(nx, ny, nz))
      call fill(x, 0)
      do idir = 1, 2
         fwd = idir == 1
         do ifct = 1, size(fcts)
            call dft3(x, ref, fwd, fcts(ifct))
            call c2c(error, x, y, fwd, fcts(ifct))
            if (allocated(error)) return
            call check(error, relerr([y], [ref]) < tol, &
               & "moist_fft_c2c_3d deviates from the direct DFT")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_direct_dft

   !* ======================================================================================== *!
   !*                                 Continuous analytic                                      *!
   !* ======================================================================================== *!

   !> Signed frequency index of the zero-based bin m on an n-point axis
   !>
   !> @param[in] m  zero-based bin, any integer
   !> @param[in] n  axis length
   pure function signed_bin(m, n) result(ms)
      !> Zero-based bin
      integer, intent(in) :: m
      !> Axis length
      integer, intent(in) :: n
      !> Signed index in [-n/2, (n-1)/2]
      integer :: ms

      ms = modulo(m, n)
      if (ms >= (n + 1)/2) ms = ms - n
   end function signed_bin

   !> Continuous Fourier transform of the test Gaussian, phase-referenced to r1
   !>
   !> - F(k) = (pi/a)^(3/2) exp(-|k|^2/(4a)) exp(-i k.c) in the convention
   !>   F(k) = integral f(r) exp(-i k.r) dV
   !> - The DFT measures phases from grid point 1, so it equals F(k) exp(i k.r1)
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

   !> Modulated off-center Gaussian matches its closed-form transform both ways
   !>
   !> - Field f(r) = exp(-a |r - c|^2) exp(i q.(r - r1)), q a grid wave vector with
   !>   bins (3, -2, 1), so the input is genuinely complex and the spectrum is the
   !>   Gaussian's shifted by q: y(m) = F(k_{m - q}) exp(i k_{m - q}.r1)
   !> - 46 x 45 x 48 points at dr = 0.25 bohr, box centered on zero; points
   !>   r_i = (i - 1 - n/2) dr, so r1 is the first point on each axis
   !> - Forward with fct = dV against the analytic spectrum; backward of the
   !>   analytic spectrum with fct = 1/Vbox against f
   !> - Aliasing error ~ exp(-(pi/dr)^2/(4a)) ~ 1e-17 of F(0); truncation error
   !>   ~ exp(-a (L/2 - |c_i|)^2), largest on the y faces (5.6 - 0.5 bohr), ~ 5e-12
   !>
   !> @param[out] error  test failure
   subroutine test_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n(3) = [46, 45, 48]
      ! n/2, spelled out to keep the truncation explicit
      integer, parameter :: nhalf(3) = [23, 22, 24]
      integer, parameter :: mq(3) = [3, -2, 1]
      real(wp), parameter :: dr = 0.25_wp
      ! Measured 2.8e-13 (forward, relative to F(0)) and 1.1e-12 (backward, absolute)
      real(wp), parameter :: tol_fwd = 2.0e-12_wp, tol_bwd = 1.0e-11_wp
      complex(cwp), allocatable :: f(:, :, :), y(:, :, :), spec(:, :, :), back(:, :, :)
      real(wp) :: r1(3), r(3), k(3), ang, dv, vbox, err
      integer :: i(3), ix, iy, iz, d

      r1 = real(-nhalf, wp)*dr
      dv = dr**3
      vbox = real(product(n), wp)*dv
      allocate (f(n(1), n(2), n(3)), y(n(1), n(2), n(3)))
      allocate (spec(n(1), n(2), n(3)), back(n(1), n(2), n(3)))
      do iz = 1, n(3)
         do iy = 1, n(2)
            do ix = 1, n(1)
               i = [ix, iy, iz]
               r = real(i - 1 - nhalf, wp)*dr
               ang = 0.0_wp
               do d = 1, 3
                  ang = ang + 2.0_wp*pi*real(modulo(mq(d)*(i(d) - 1), n(d)), wp)/real(n(d), wp)
               end do
               f(ix, iy, iz) = exp(-gauss_a*sum((r - gauss_c)**2))*cmplx(cos(ang), sin(ang), wp)
               do d = 1, 3
                  k(d) = 2.0_wp*pi*real(signed_bin(i(d) - 1 - mq(d), n(d)), wp)/(real(n(d), wp)*dr)
               end do
               spec(ix, iy, iz) = gaussian_ft(k, r1)
            end do
         end do
      end do

      call c2c(error, f, y, .true., dv)
      if (allocated(error)) return
      err = maxval(abs(y - spec))/(pi/gauss_a)**1.5_wp
      call check(error, err < tol_fwd, &
         & "forward c2c of the modulated Gaussian deviates from its analytic transform")
      if (allocated(error)) return

      call c2c(error, spec, back, .false., 1.0_wp/vbox)
      if (allocated(error)) return
      err = maxval(abs(back - f))
      call check(error, err < tol_bwd, &
         & "backward c2c of the analytic spectrum deviates from the modulated Gaussian")
   end subroutine test_gaussian

   !* ======================================================================================== *!
   !*                                      Round trips                                         *!
   !* ======================================================================================== *!

   !> backward(forward(x)) reproduces x for the normalized factor pairs
   !>
   !> - 12 x 10 x 9 grid; (1, 1/N) and the unitary (1/sqrt(N), 1/sqrt(N))
   !>
   !> @param[out] error  test failure
   subroutine test_round_trip(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 12, ny = 10, nz = 9
      ! Measured 7.2e-16
      real(wp), parameter :: tol = 5.0e-15_wp
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
      real(wp) :: fct_f(2), fct_b(2), ntot
      integer :: ip

      ntot = real(nx*ny*nz, wp)
      fct_f = [1.0_wp, 1.0_wp/sqrt(ntot)]
      fct_b = [1.0_wp/ntot, 1.0_wp/sqrt(ntot)]
      allocate (x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz))
      call fill(x, 17)
      do ip = 1, 2
         call c2c(error, x, y, .true., fct_f(ip))
         if (allocated(error)) return
         call c2c(error, y, z, .false., fct_b(ip))
         if (allocated(error)) return
         call check(error, relerr([z], [x]) < tol, &
            & "backward(forward(x)) does not reproduce x")
         if (allocated(error)) return
      end do
   end subroutine test_round_trip

   !> Six single-axis passes with fct = 1/n_axis on the way back reproduce x
   !>
   !> - Forward with pass (x to y on axis 2) then pass_inplace on axes 1 and 0
   !> - Backward with pass_inplace on axes 0, 1, 2, each scaled by its own 1/n_axis
   !> - Only the product of the factors is visible here; the per-pass factor is
   !>   pinned by `pass_inplace_matches_pass` and `pass_strided_subblock`
   !>
   !> @param[out] error  test failure
   subroutine test_pass_round_trip(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 12, ny = 10, nz = 9
      ! Measured 6.9e-16
      real(wp), parameter :: tol = 5.0e-15_wp
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :)
      integer :: axis, nax(0:2)

      nax = [nz, ny, nx]
      allocate (x(nx, ny, nz), y(nx, ny, nz))
      call fill(x, 29)
      y = (0.0_wp, 0.0_wp)
      call pass(error, x, y, 2, .true., 1.0_wp)
      if (allocated(error)) return
      do axis = 1, 0, -1
         call pass_inplace(error, y, axis, .true., 1.0_wp)
         if (allocated(error)) return
      end do
      do axis = 0, 2
         call pass_inplace(error, y, axis, .false., 1.0_wp/real(nax(axis), wp))
         if (allocated(error)) return
      end do
      call check(error, relerr([y], [x]) < tol, &
         & "backward passes of the forward passes do not reproduce x")
   end subroutine test_pass_round_trip

   !* ======================================================================================== *!
   !*                                      Consistency                                         *!
   !* ======================================================================================== *!

   !> Three single-axis passes compose to the full 3D transform
   !>
   !> - Chain A mirrors ducc0's own axis order for an out-of-place call: pass on
   !>   axis 2 (x to y), then pass_inplace on axes 1 and 0, fct on the last pass
   !> - Chain B runs the opposite order with three out-of-place passes, fct on the
   !>   first, so it also checks that the passes commute
   !> - Both directions, 7 x 6 x 5 grid
   !>
   !> @param[out] error  test failure
   subroutine test_passes_compose(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 7, ny = 6, nz = 5
      real(wp), parameter :: fct = 0.73_wp
      ! Measured 2.3e-16 (chain A) and 1.6e-16 (chain B)
      real(wp), parameter :: tol = 2.0e-15_wp
      complex(cwp), allocatable :: x(:, :, :), full(:, :, :), a(:, :, :), t1(:, :, :), t2(:, :, :)
      integer :: idir
      logical :: fwd

      allocate (x(nx, ny, nz), full(nx, ny, nz), a(nx, ny, nz), t1(nx, ny, nz), t2(nx, ny, nz))
      call fill(x, 41)
      do idir = 1, 2
         fwd = idir == 1
         call c2c(error, x, full, fwd, fct)
         if (allocated(error)) return

         a = (0.0_wp, 0.0_wp)
         call pass(error, x, a, 2, fwd, 1.0_wp)
         if (.not. allocated(error)) call pass_inplace(error, a, 1, fwd, 1.0_wp)
         if (.not. allocated(error)) call pass_inplace(error, a, 0, fwd, fct)
         if (allocated(error)) return
         call check(error, relerr([a], [full]) < tol, &
            & "pass axis 2 then pass_inplace axes 1, 0 differ from moist_fft_c2c_3d")
         if (allocated(error)) return

         t1 = (0.0_wp, 0.0_wp); t2 = (0.0_wp, 0.0_wp); a = (0.0_wp, 0.0_wp)
         call pass(error, x, t1, 0, fwd, fct)
         if (.not. allocated(error)) call pass(error, t1, t2, 1, fwd, 1.0_wp)
         if (.not. allocated(error)) call pass(error, t2, a, 2, fwd, 1.0_wp)
         if (allocated(error)) return
         call check(error, relerr([a], [full]) < tol, &
            & "out-of-place passes on axes 0, 1, 2 differ from moist_fft_c2c_3d")
         if (allocated(error)) return
      end do
   end subroutine test_passes_compose

   !> pass_inplace gives the same result as the out-of-place pass
   !>
   !> - Every axis, both directions, fct = 1.3; each pass also matches a direct
   !>   1D DFT along its axis
   !>
   !> @param[out] error  test failure
   subroutine test_inplace_matches_pass(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 7, ny = 6, nz = 5
      real(wp), parameter :: fct = 1.3_wp
      ! Measured 0 (in place vs out of place; bound at a few ulp) and 4.0e-16 (vs DFT)
      real(wp), parameter :: tol_same = 1.0e-15_wp, tol_ref = 3.0e-15_wp
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :), ref(:, :, :)
      integer :: axis, idir
      logical :: fwd

      allocate (x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz), ref(nx, ny, nz))
      call fill(x, 53)
      do idir = 1, 2
         fwd = idir == 1
         do axis = 0, 2
            y = (0.0_wp, 0.0_wp)
            call pass(error, x, y, axis, fwd, fct)
            if (allocated(error)) return
            z = x
            call pass_inplace(error, z, axis, fwd, fct)
            if (allocated(error)) return
            call dft_axis(x, ref, axis, fwd, fct)
            call check(error, relerr([z], [y]) < tol_same, &
               & "pass_inplace differs from the out-of-place pass")
            if (allocated(error)) return
            call check(error, relerr([y], [ref]) < tol_ref, &
               & "single-axis pass deviates from the direct 1D DFT along its axis")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_inplace_matches_pass

   !> Strided sub-block passes touch exactly the addressed elements
   !>
   !> - Parent field 9 x 8 x 7, sub-block 5 x 3 x 4 starting at (3, 4, 2), so
   !>   s1 = 9 and s0 = 72; the block start is passed by sequence association
   !> - Out of place: the out parent is prefilled with a sentinel; its block equals
   !>   the direct 1D DFT of the in block, every other element keeps the sentinel
   !>   bit for bit, and the in parent is unchanged bit for bit
   !> - In place: the block is transformed, everything else is unchanged bit for bit
   !> - Every axis, both directions
   !>
   !> @param[out] error  test failure
   subroutine test_strided_subblock(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: px = 9, py = 8, pz = 7
      integer, parameter :: mx = 5, my = 3, mz = 4
      integer, parameter :: i0 = 3, j0 = 4, k0 = 2
      real(wp), parameter :: fct = 0.9_wp
      complex(cwp), parameter :: sentinel = (-7.25_wp, 3.5_wp)
      ! Measured 3.8e-16
      real(wp), parameter :: tol = 3.0e-15_wp
      complex(cwp), allocatable :: a(:, :, :), a0(:, :, :), b(:, :, :), ref(:, :, :)
      logical :: inblock(px, py, pz)
      integer(c_int) :: status
      integer :: axis, idir
      logical :: fwd

      allocate (a(px, py, pz), a0(px, py, pz), b(px, py, pz), ref(mx, my, mz))
      call fill(a0, 67)
      inblock = .false.
      inblock(i0:i0 + mx - 1, j0:j0 + my - 1, k0:k0 + mz - 1) = .true.

      do idir = 1, 2
         fwd = idir == 1
         do axis = 0, 2
            call dft_axis(a0(i0:i0 + mx - 1, j0:j0 + my - 1, k0:k0 + mz - 1), ref, axis, fwd, fct)

            a = a0
            b = sentinel
            status = moist_fft_c2c_3d_pass(int(mz, c_int), int(my, c_int), int(mx, c_int), &
               & int(px*py, c_int), int(px, c_int), a(i0, j0, k0), b(i0, j0, k0), &
               & int(axis, c_int), dir(fwd), fct)
            call check(error, status == 0, "strided moist_fft_c2c_3d_pass returned a nonzero status")
            if (allocated(error)) return
            call check(error, relerr([b(i0:i0 + mx - 1, j0:j0 + my - 1, k0:k0 + mz - 1)], [ref]) < tol, &
               & "strided out-of-place pass deviates from the direct 1D DFT of the block")
            if (allocated(error)) return
            call check(error, all(pack(b, .not. inblock) == sentinel), &
               & "strided out-of-place pass wrote outside its sub-block")
            if (allocated(error)) return
            call check(error, all(a == a0), "strided out-of-place pass modified its input")
            if (allocated(error)) return

            a = a0
            status = moist_fft_c2c_3d_pass_inplace(int(mz, c_int), int(my, c_int), int(mx, c_int), &
               & int(px*py, c_int), int(px, c_int), a(i0, j0, k0), int(axis, c_int), dir(fwd), fct)
            call check(error, status == 0, "strided moist_fft_c2c_3d_pass_inplace returned a nonzero status")
            if (allocated(error)) return
            call check(error, relerr([a(i0:i0 + mx - 1, j0:j0 + my - 1, k0:k0 + mz - 1)], [ref]) < tol, &
               & "strided in-place pass deviates from the direct 1D DFT of the block")
            if (allocated(error)) return
            call check(error, all(pack(a, .not. inblock) == pack(a0, .not. inblock)), &
               & "strided in-place pass wrote outside its sub-block")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_strided_subblock

   !> c2c_3d is linear: F(alpha x + beta y) = alpha F(x) + beta F(y)
   !>
   !> - Complex alpha, real beta, both directions, 8 x 5 x 6 grid
   !>
   !> @param[out] error  test failure
   subroutine test_linearity(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 8, ny = 5, nz = 6
      complex(wp), parameter :: alpha = (0.3_wp, -1.2_wp)
      real(wp), parameter :: beta = 2.5_wp
      ! Measured 1.8e-16
      real(wp), parameter :: tol = 1.0e-15_wp
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :), fx(:, :, :), fy(:, :, :), fxy(:, :, :)
      integer :: idir
      logical :: fwd

      allocate (x(nx, ny, nz), y(nx, ny, nz), fx(nx, ny, nz), fy(nx, ny, nz), fxy(nx, ny, nz))
      call fill(x, 71)
      call fill(y, 911)
      do idir = 1, 2
         fwd = idir == 1
         call c2c(error, x, fx, fwd, 1.0_wp)
         if (.not. allocated(error)) call c2c(error, y, fy, fwd, 1.0_wp)
         if (.not. allocated(error)) call c2c(error, alpha*x + beta*y, fxy, fwd, 1.0_wp)
         if (allocated(error)) return
         call check(error, relerr([fxy], [alpha*fx + beta*fy]) < tol, &
            & "moist_fft_c2c_3d is not linear")
         if (allocated(error)) return
      end do
   end subroutine test_linearity

   !> Parseval: sum |X|^2 = N sum |x|^2 at fct = 1, equal norms at fct = 1/sqrt(N)
   !>
   !> - Both directions, 9 x 7 x 4 grid
   !>
   !> @param[out] error  test failure
   subroutine test_parseval(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 9, ny = 7, nz = 4
      ! Measured 4.3e-16 (fct = 1) and 2.2e-16 (fct = 1/sqrt(N)), relative
      real(wp), parameter :: tol = 3.0e-15_wp
      complex(cwp), allocatable :: x(:, :, :), y(:, :, :)
      real(wp) :: ntot, ex, ey
      integer :: idir
      logical :: fwd

      ntot = real(nx*ny*nz, wp)
      allocate (x(nx, ny, nz), y(nx, ny, nz))
      call fill(x, 83)
      ex = sum(abs(x)**2)
      do idir = 1, 2
         fwd = idir == 1
         call c2c(error, x, y, fwd, 1.0_wp)
         if (allocated(error)) return
         ey = sum(abs(y)**2)
         call check(error, abs(ey - ntot*ex) < tol*ntot*ex, &
            & "unnormalized transform violates sum |X|^2 = N sum |x|^2")
         if (allocated(error)) return
         call c2c(error, x, y, fwd, 1.0_wp/sqrt(ntot))
         if (allocated(error)) return
         ey = sum(abs(y)**2)
         call check(error, abs(ey - ex) < tol*ex, &
            & "unitary-scaled transform does not preserve the norm")
         if (allocated(error)) return
      end do
   end subroutine test_parseval

   !* ======================================================================================== *!
   !*                                    Real transforms                                       *!
   !* ======================================================================================== *!

   !> r2c output equals the non-redundant half of c2c on real input
   !>
   !> - Half spectrum along x: (nx/2 + 1, ny, nz); odd and even nx
   !> - fct = 0.9 on both
   !>
   !> @param[out] error  test failure
   subroutine test_r2c_half(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ny = 6, nz = 5
      integer, parameter :: nxs(2) = [7, 8]
      real(wp), parameter :: fct = 0.9_wp
      ! Measured 1.0e-16
      real(wp), parameter :: tol = 1.0e-15_wp
      real(c_double), allocatable :: xr(:, :, :)
      complex(cwp), allocatable :: xc(:, :, :), full(:, :, :), half(:, :, :)
      integer(c_int) :: status
      integer :: in, nx, nh

      do in = 1, size(nxs)
         nx = nxs(in)
         nh = nx/2 + 1
         allocate (xc(nx, ny, nz), full(nx, ny, nz), half(nh, ny, nz))
         call fill(xc, 97*in)
         xr = real(xc, wp)
         xc = cmplx(xr, 0.0_wp, wp)
         call c2c(error, xc, full, .true., fct)
         if (allocated(error)) return
         status = moist_fft_r2c_3d(int(nz, c_int), int(ny, c_int), int(nx, c_int), xr, half, fct)
         call check(error, status == 0, "moist_fft_r2c_3d returned a nonzero status")
         if (allocated(error)) return
         call check(error, relerr([half], [full(1:nh, :, :)]) < tol, &
            & "r2c differs from the non-redundant half of c2c")
         if (allocated(error)) return
         deallocate (xc, full, half)
      end do
   end subroutine test_r2c_half

   !> c2r equals the real c2c backward of a Hermitian spectrum and inverts r2c
   !>
   !> - Hermitian spectrum: forward c2c of a real field; c2r of its half with
   !>   fct = 0.4 against the real part of the full c2c backward with fct = 0.4
   !> - c2r(r2c(x)) with fct = 1/N reproduces x
   !> - c2r destroys its input, so each call gets its own copy
   !> - Odd and even nx
   !>
   !> @param[out] error  test failure
   subroutine test_c2r_backward(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ny = 5, nz = 6
      integer, parameter :: nxs(2) = [9, 8]
      real(wp), parameter :: fct = 0.4_wp
      ! Measured 4.6e-16 (vs c2c) and 6.8e-16 (round trip)
      real(wp), parameter :: tol = 5.0e-15_wp
      real(c_double), allocatable :: xr(:, :, :), yr(:, :, :)
      complex(cwp), allocatable :: xc(:, :, :), full(:, :, :), back(:, :, :), half(:, :, :)
      integer(c_int) :: status
      integer :: in, nx, nh

      do in = 1, size(nxs)
         nx = nxs(in)
         nh = nx/2 + 1
         allocate (xc(nx, ny, nz), full(nx, ny, nz), back(nx, ny, nz), half(nh, ny, nz))
         allocate (yr(nx, ny, nz))
         call fill(xc, 131*in)
         xr = real(xc, wp)
         xc = cmplx(xr, 0.0_wp, wp)
         call c2c(error, xc, full, .true., 1.0_wp)
         if (allocated(error)) return
         call c2c(error, full, back, .false., fct)
         if (allocated(error)) return

         half = full(1:nh, :, :)
         status = moist_fft_c2r_3d(int(nz, c_int), int(ny, c_int), int(nx, c_int), half, yr, fct)
         call check(error, status == 0, "moist_fft_c2r_3d returned a nonzero status")
         if (allocated(error)) return
         call check(error, relerr_real([yr], [real(back, wp)]) < tol, &
            & "c2r differs from the real c2c backward of the Hermitian spectrum")
         if (allocated(error)) return

         status = moist_fft_r2c_3d(int(nz, c_int), int(ny, c_int), int(nx, c_int), xr, half, 1.0_wp)
         if (status == 0) status = moist_fft_c2r_3d(int(nz, c_int), int(ny, c_int), int(nx, c_int), &
            & half, yr, 1.0_wp/real(nx*ny*nz, wp))
         call check(error, status == 0, "r2c/c2r round trip returned a nonzero status")
         if (allocated(error)) return
         call check(error, relerr_real([yr], [xr]) < tol, "c2r(r2c(x))/N does not reproduce x")
         if (allocated(error)) return
         deallocate (xc, full, back, half, yr)
      end do
   end subroutine test_c2r_backward

   !> Batched r2c/c2r reproduce per-site single-field calls
   !>
   !> - Layout (nx, ny, nz, nb), site slowest; half spectrum (nx/2 + 1, ny, nz, nb)
   !> - nb = 3 on a 7 x 6 x 5 grid, nthreads = 1 and 4; inside test-drive's
   !>   parallel region nthreads = 4 falls back to the serial transform, the
   !>   threaded path is covered by `math_grid_3d_threaded`
   !> - c2r destroys its input, so each call gets its own copy
   !>
   !> @param[out] error  test failure
   subroutine test_batch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 7, ny = 6, nz = 5, nb = 3
      ! nx/2 + 1
      integer, parameter :: nh = 4
      integer, parameter :: nthreads(2) = [1, 4]
      real(wp), parameter :: fct_f = 0.8_wp, fct_b = 1.7_wp
      ! Measured 0 (forward and backward; bound at a few ulp)
      real(wp), parameter :: tol = 1.0e-15_wp
      real(c_double), allocatable :: xr(:, :, :, :), yr(:, :, :, :), yr_ref(:, :, :, :)
      complex(cwp), allocatable :: xc(:, :, :), h(:, :, :, :), h_ref(:, :, :, :), scratch(:, :, :, :)
      integer(c_int) :: status
      integer :: s, it

      allocate (xr(nx, ny, nz, nb), yr(nx, ny, nz, nb), yr_ref(nx, ny, nz, nb))
      allocate (xc(nx, ny, nz), h(nh, ny, nz, nb), h_ref(nh, ny, nz, nb), scratch(nh, ny, nz, nb))
      do s = 1, nb
         call fill(xc, 1000*s)
         xr(:, :, :, s) = real(xc, wp)
      end do

      status = 0
      do s = 1, nb
         if (status == 0) status = moist_fft_r2c_3d(int(nz, c_int), int(ny, c_int), int(nx, c_int), &
            & xr(:, :, :, s), h_ref(:, :, :, s), fct_f)
         scratch(:, :, :, s) = h_ref(:, :, :, s)
         if (status == 0) status = moist_fft_c2r_3d(int(nz, c_int), int(ny, c_int), int(nx, c_int), &
            & scratch(:, :, :, s), yr_ref(:, :, :, s), fct_b)
      end do
      call check(error, status == 0, "per-site reference transform returned a nonzero status")
      if (allocated(error)) return

      do it = 1, size(nthreads)
         status = moist_fft_r2c_3d_batch(int(nz, c_int), int(ny, c_int), int(nx, c_int), int(nb, c_int), &
            & xr, h, fct_f, int(nthreads(it), c_int))
         call check(error, status == 0, "moist_fft_r2c_3d_batch returned a nonzero status")
         if (allocated(error)) return
         call check(error, relerr([h], [h_ref]) < tol, &
            & "batched r2c differs from per-site single-field calls")
         if (allocated(error)) return

         scratch = h_ref
         status = moist_fft_c2r_3d_batch(int(nz, c_int), int(ny, c_int), int(nx, c_int), int(nb, c_int), &
            & scratch, yr, fct_b, int(nthreads(it), c_int))
         call check(error, status == 0, "moist_fft_c2r_3d_batch returned a nonzero status")
         if (allocated(error)) return
         call check(error, relerr_real([yr], [yr_ref]) < tol, &
            & "batched c2r differs from per-site single-field calls")
         if (allocated(error)) return
      end do
   end subroutine test_batch

   !* ======================================================================================== *!
   !*                              Degenerate and error paths                                  *!
   !* ======================================================================================== *!

   !> Zero extents return success and touch nothing
   !>
   !> - ducc0 returns early on an empty view; one zero extent on each entry
   !>   point, output buffers prefilled with a sentinel stay bit-identical
   !>
   !> @param[out] error  test failure
   subroutine test_zero_extent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      complex(cwp), parameter :: sentinel = (-7.25_wp, 3.5_wp)
      real(c_double), parameter :: rsentinel = -7.25_wp
      complex(cwp) :: x(8), y(8)
      real(c_double) :: xr(8), yr(8)
      integer(c_int) :: status(7)

      x = (1.0_wp, 2.0_wp)
      xr = 1.0_wp
      y = sentinel
      yr = rsentinel
      status(1) = moist_fft_c2c_3d(0_c_int, 2_c_int, 2_c_int, x, y, 1_c_int, 1.0_wp)
      status(2) = moist_fft_c2c_3d_pass(2_c_int, 0_c_int, 2_c_int, 4_c_int, 2_c_int, x, y, &
         & 0_c_int, 1_c_int, 1.0_wp)
      status(3) = moist_fft_c2c_3d_pass_inplace(2_c_int, 2_c_int, 0_c_int, 0_c_int, 0_c_int, y, &
         & 1_c_int, 0_c_int, 1.0_wp)
      status(4) = moist_fft_r2c_3d(2_c_int, 0_c_int, 2_c_int, xr, y, 1.0_wp)
      status(5) = moist_fft_c2r_3d(0_c_int, 2_c_int, 2_c_int, x, yr, 1.0_wp)
      status(6) = moist_fft_r2c_3d_batch(2_c_int, 2_c_int, 2_c_int, 0_c_int, xr, y, 1.0_wp, 1_c_int)
      status(7) = moist_fft_c2r_3d_batch(2_c_int, 2_c_int, 2_c_int, 0_c_int, x, yr, 1.0_wp, 1_c_int)
      call check(error, all(status == 0), "zero extent did not return a zero status")
      if (allocated(error)) return
      call check(error, all(y == sentinel) .and. all(yr == rsentinel), &
         & "zero-extent transform wrote to its output")
      if (allocated(error)) return
      call check(error, all(x == (1.0_wp, 2.0_wp)) .and. all(xr == 1.0_wp), &
         & "zero-extent transform modified its input")
   end subroutine test_zero_extent

   !> Negative extents return status 1 and touch nothing
   !>
   !> - Checked in the shim before any view is built; a negative extent would
   !>   otherwise wrap to a huge size_t and address far outside the buffers
   !> - One negative extent on each entry point, buffers prefilled with a
   !>   sentinel stay bit-identical
   !>
   !> @param[out] error  test failure
   subroutine test_negative_extent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      complex(cwp), parameter :: sentinel = (-7.25_wp, 3.5_wp)
      real(c_double), parameter :: rsentinel = -7.25_wp
      complex(cwp) :: x(8), y(8)
      real(c_double) :: xr(8), yr(8)
      integer(c_int) :: status(9)

      x = (1.0_wp, 2.0_wp)
      xr = 1.0_wp
      y = sentinel
      yr = rsentinel
      status(1) = moist_fft_c2c_3d(-1_c_int, 2_c_int, 2_c_int, x, y, 1_c_int, 1.0_wp)
      status(2) = moist_fft_c2c_3d_pass(2_c_int, -2_c_int, 2_c_int, 4_c_int, 2_c_int, x, y, &
         & 0_c_int, 1_c_int, 1.0_wp)
      status(3) = moist_fft_c2c_3d_pass_inplace(2_c_int, 2_c_int, -2_c_int, 4_c_int, 2_c_int, y, &
         & 1_c_int, 0_c_int, 1.0_wp)
      status(4) = moist_fft_r2c_3d(2_c_int, -1_c_int, 2_c_int, xr, y, 1.0_wp)
      status(5) = moist_fft_c2r_3d(2_c_int, 2_c_int, -2_c_int, x, yr, 1.0_wp)
      status(6) = moist_fft_r2c_3d_batch(2_c_int, 2_c_int, 2_c_int, -1_c_int, xr, y, 1.0_wp, 1_c_int)
      status(7) = moist_fft_c2r_3d_batch(-2_c_int, 2_c_int, 2_c_int, 1_c_int, x, yr, 1.0_wp, 2_c_int)
      status(8) = moist_fft_dst4(-4_c_int, 1_c_int, xr, yr, 1.0_wp)
      status(9) = moist_fft_dst4(4_c_int, -1_c_int, xr, yr, 1.0_wp)
      call check(error, all(status == status_error), "negative extent did not return status 1")
      if (allocated(error)) return
      call check(error, all(y == sentinel) .and. all(yr == rsentinel), &
         & "negative-extent transform wrote to its output")
      if (allocated(error)) return
      call check(error, all(x == (1.0_wp, 2.0_wp)) .and. all(xr == 1.0_wp), &
         & "negative-extent transform modified its input")
   end subroutine test_negative_extent

   !> Rejected pass reports status 1 and leaves both buffers untouched
   !>
   !> - Axis 3 on pass, axis -1 on pass_inplace; ducc0 validates axes before
   !>   touching data
   !>
   !> @param[out] error  test failure
   subroutine test_rejected_pass(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nx = 4, ny = 3, nz = 2
      complex(cwp), parameter :: sentinel = (-7.25_wp, 3.5_wp)
      complex(cwp), allocatable :: x(:, :, :), x0(:, :, :), y(:, :, :)
      integer(c_int) :: status

      allocate (x(nx, ny, nz), x0(nx, ny, nz), y(nx, ny, nz))
      call fill(x0, 5)
      x = x0
      y = sentinel
      status = moist_fft_c2c_3d_pass(int(nz, c_int), int(ny, c_int), int(nx, c_int), &
         & int(nx*ny, c_int), int(nx, c_int), x, y, 3_c_int, 1_c_int, 1.0_wp)
      call check(error, status == status_error, "pass with axis 3 did not return status 1")
      if (allocated(error)) return
      call check(error, all(x == x0) .and. all(y == sentinel), &
         & "rejected pass modified its buffers")
      if (allocated(error)) return

      status = moist_fft_c2c_3d_pass_inplace(int(nz, c_int), int(ny, c_int), int(nx, c_int), &
         & int(nx*ny, c_int), int(nx, c_int), x, -1_c_int, 1_c_int, 1.0_wp)
      call check(error, status == status_error, "pass_inplace with axis -1 did not return status 1")
      if (allocated(error)) return
      call check(error, all(x == x0), "rejected pass_inplace modified its data")
   end subroutine test_rejected_pass

   !> Out-of-place pass with axis 3 fails with status 1
   !>
   !> - `should_fail`: fails the test only on that status, so a clean return is
   !>   reported as an unexpected pass
   !>
   !> @param[out] error  test failure, set only on the expected status
   subroutine test_pass_axis_large_fails(error)
      !> Test failure, set only on the expected status
      type(error_type), allocatable, intent(out) :: error

      complex(cwp) :: x(2, 2, 2), y(2, 2, 2)
      integer(c_int) :: status

      x = (1.0_wp, 0.0_wp)
      y = (0.0_wp, 0.0_wp)
      status = moist_fft_c2c_3d_pass(2_c_int, 2_c_int, 2_c_int, 4_c_int, 2_c_int, x, y, &
         & 3_c_int, 1_c_int, 1.0_wp)
      if (status == status_error) call test_failed(error, "pass rejected axis 3 with status 1")
   end subroutine test_pass_axis_large_fails

   !> Out-of-place pass with axis -1 fails with status 1
   !>
   !> - Negative axis wraps to a huge size_t in the shim, caught by ducc0's axis check
   !>
   !> @param[out] error  test failure, set only on the expected status
   subroutine test_pass_axis_negative_fails(error)
      !> Test failure, set only on the expected status
      type(error_type), allocatable, intent(out) :: error

      complex(cwp) :: x(2, 2, 2), y(2, 2, 2)
      integer(c_int) :: status

      x = (1.0_wp, 0.0_wp)
      y = (0.0_wp, 0.0_wp)
      status = moist_fft_c2c_3d_pass(2_c_int, 2_c_int, 2_c_int, 4_c_int, 2_c_int, x, y, &
         & -1_c_int, 0_c_int, 1.0_wp)
      if (status == status_error) call test_failed(error, "pass rejected axis -1 with status 1")
   end subroutine test_pass_axis_negative_fails

   !> In-place pass with axis 3 fails with status 1
   !>
   !> @param[out] error  test failure, set only on the expected status
   subroutine test_pass_inplace_axis_fails(error)
      !> Test failure, set only on the expected status
      type(error_type), allocatable, intent(out) :: error

      complex(cwp) :: x(2, 2, 2)
      integer(c_int) :: status

      x = (1.0_wp, 0.0_wp)
      status = moist_fft_c2c_3d_pass_inplace(2_c_int, 2_c_int, 2_c_int, 4_c_int, 2_c_int, x, &
         & 3_c_int, 1_c_int, 1.0_wp)
      if (status == status_error) call test_failed(error, "pass_inplace rejected axis 3 with status 1")
   end subroutine test_pass_inplace_axis_fails

end module test_math_fft_3d
