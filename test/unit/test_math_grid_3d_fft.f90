!> Cartesian 3D grid FFT accuracy against analytic Gaussian transforms
!>
!> - Forward FFT of off-center Gaussians vs their analytic transforms, three columns
!> - Inverse FFT of the analytic spectra vs the Gaussians, three columns
!> - Spectral gradient vs the analytic gradient, per axis
!> - k-space convolution of two off-center Gaussians vs the analytic convolution
module test_math_grid_3d_fft
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: get_cartesian_gaussian_grid
   use moist_math_grid_3d_base, only: moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_fft

   !> Exponent of the analytic test Gaussian exp(-a |r - c|^2), 1/bohr^2
   real(wp), parameter :: gauss_a = 1.0_wp
   !> Off-center, off-axis Gaussian center c, bohr
   real(wp), parameter :: gauss_c(3) = [0.3_wp, -0.5_wp, 0.7_wp]
   !> Distinct centers of the batched columns, bohr; a skipped or swapped column fails
   real(wp), parameter :: batch_c(3, 3) = reshape([0.3_wp, -0.5_wp, 0.7_wp, &
                                                   -0.6_wp, 0.4_wp, -0.2_wp, &
                                                   0.1_wp, 0.6_wp, -0.5_wp], [3, 3])

contains

   !> Collect all math_grid_3d_fft tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_3d_fft(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fft_gaussian_analytic", test_fft_gaussian), &
                  new_unittest("ifft_gaussian_analytic", test_ifft_gaussian), &
                  new_unittest("fft_gaussian_spectral_gradient", test_fft_gaussian_gradient), &
                  new_unittest("fft_convolution_analytic", test_fft_convolution) &
                  ]
   end subroutine collect_math_grid_3d_fft

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

   !> Gaussian exp(-a |r - c|^2) of free exponent and center
   !>
   !> @param[in] r  point, bohr
   !> @param[in] c  center, bohr
   !> @param[in] a  exponent, 1/bohr^2
   pure function gaussian_at(r, c, a) result(val)
      !> Point, bohr
      real(wp), intent(in) :: r(3)
      !> Center, bohr
      real(wp), intent(in) :: c(3)
      !> Exponent, 1/bohr^2
      real(wp), intent(in) :: a
      !> Field value
      real(wp) :: val

      val = exp(-a*sum((r - c)**2))
   end function gaussian_at

   !> Continuous Fourier transform of exp(-a |r - c|^2), phase-referenced to r1
   !>
   !> - F(k) = (pi/a)^(3/2) exp(-|k|^2/(4a)) exp(-i k.c) in the convention
   !>   F(k) = integral f(r) exp(-i k.r) dV
   !> - Cartesian FFT phases are referenced to grid point 1, so the grid
   !>   transform equals F(k) exp(i k.r1)
   !>
   !> @param[in] k   wave vector, 1/bohr
   !> @param[in] r1  phase reference, bohr
   !> @param[in] c   Gaussian center, bohr
   pure function gaussian_ft(k, r1, c) result(val)
      !> Wave vector, 1/bohr
      real(wp), intent(in) :: k(3)
      !> Phase reference, bohr
      real(wp), intent(in) :: r1(3)
      !> Gaussian center, bohr
      real(wp), intent(in) :: c(3)
      !> Transform value
      complex(wp) :: val

      real(wp) :: ph

      ph = -dot_product(k, c - r1)
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
   !> - One worker selects the per-column backend; `math_grid_3d_threaded`
   !>   covers the batched one
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

      call get_cartesian_gaussian_grid(grid, 45, 48, 47, 0.25_wp, error=merr)
      grid%nthreads = 1
      if (.not. allocated(merr)) call grid%new_trafo(trafo, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine gaussian_grid

   !> Forward FFT of off-center Gaussians equals their analytic transforms
   !>
   !> - Compared as complex values on every k-point of the half spectrum
   !> - Pins the kpoint() layout and signs on all three axes, the phase
   !>   reference (grid point 1 = origin + dr/2) and the dV normalization
   !> - Three columns with distinct centers in one serial call
   subroutine test_fft_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: f_r(:, :)
      complex(wp), allocatable :: f_k(:, :), ref(:)
      real(wp) :: r1(3)
      integer :: i, j, iv
      character(len=1) :: col

      call gaussian_grid(grid, trafo, error)
      if (allocated(error)) return
      r1 = grid%get_origin() + 0.5_wp*grid%dr
      allocate (f_r(grid%ngrid, 3), f_k(grid%npts_k, 3), ref(grid%npts_k))
      do iv = 1, 3
         do i = 1, grid%ngrid
            f_r(i, iv) = gaussian_at(grid%xyz(:, i), batch_c(:, iv), gauss_a)
         end do
      end do

      call forward(error, trafo, f_r, f_k)
      if (allocated(error)) return
      do iv = 1, 3
         do j = 1, grid%npts_k
            ref(j) = gaussian_ft(grid%kpoint(j), r1, batch_c(:, iv))
         end do
         write (col, "(i1)") iv
         call check(error, all(abs(f_k(:, iv) - ref) <= 1.0e-10_wp*(pi/gauss_a)**1.5_wp), &
            & "forward FFT of Gaussian "//col//" deviates from its analytic transform")
         if (allocated(error)) return
      end do

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fft_gaussian

   !> Inverse FFT of the analytic Gaussian spectra reproduces the Gaussians
   !>
   !> - Only check on the absolute 1/Vbox normalization; a round trip cannot
   !>   see a common scale error of the two directions
   !> - Analytic half spectrum is Hermitian by construction, as C2R requires
   !> - Three columns with distinct centers in one serial call
   subroutine test_ifft_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: f_r(:, :), ref(:)
      complex(wp), allocatable :: f_k(:, :)
      real(wp) :: r1(3)
      integer :: i, j, iv
      character(len=1) :: col

      call gaussian_grid(grid, trafo, error)
      if (allocated(error)) return
      r1 = grid%get_origin() + 0.5_wp*grid%dr
      allocate (f_r(grid%ngrid, 3), f_k(grid%npts_k, 3), ref(grid%ngrid))
      do iv = 1, 3
         do j = 1, grid%npts_k
            f_k(j, iv) = gaussian_ft(grid%kpoint(j), r1, batch_c(:, iv))
         end do
      end do

      call backward(error, trafo, f_k, f_r)
      if (allocated(error)) return
      do iv = 1, 3
         do i = 1, grid%ngrid
            ref(i) = gaussian_at(grid%xyz(:, i), batch_c(:, iv), gauss_a)
         end do
         write (col, "(i1)") iv
         call check(error, all(abs(f_r(:, iv) - ref) <= 1.0e-10_wp), &
            & "inverse FFT of spectrum "//col//" deviates from its Gaussian")
         if (allocated(error)) return
      end do

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

   !> k-space convolution of two off-center Gaussians matches the analytic convolution
   !>
   !> - g_a = exp(-a |r - c_a|^2), g_b = exp(-b |r - c_b|^2), both forward transformed
   !> - Product times exp(-i k.r1): each forward carries exp(i k.r1), the backward expects one
   !> - (g_a * g_b)(r) = (pi/(a+b))^(3/2) exp(-a b/(a+b) |r - (c_a + c_b)|^2), point by point
   !> - Cartesian counterpart of the molecular `convolution_matches_analytic` NUFFT test
   !> - Pins the dV and 1/Vbox normalizations together with the phase reference: a missing
   !>   r1 correction shifts the result by r1 (half the box), a wrong measure scales it
   !> - Expected error ~ 1e-10 absolute, from the periodic image of the convolution
   !>   across the x faces, exponent ab/(a+b) = 0.75:
   !>   exp(-0.75 (L_x - 5.5 - 0.2)^2) ~ 1e-10
   !> - Band-limit and aliasing terms below 1e-17: product spectrum exp(-|k|^2/3) at
   !>   k = pi/dr, cross term F_b(k - 2 pi/dr) F_a(k) <= exp(-39)
   !> - Single-field truncation exp(-a (L/2 - |c_i|)^2) <= 1e-12
   !> - Tolerance 3e-10 absolute against a peak of 0.70; measured maximum 3e-11 to 1e-10
   subroutine test_fft_convolution(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Exponent of the second Gaussian, 1/bohr^2
      real(wp), parameter :: conv_b = 3.0_wp
      !> Center of the second Gaussian, bohr; c_a + c_b = (-0.2, -0.3, 0.3)
      real(wp), parameter :: conv_c(3) = [-0.5_wp, 0.2_wp, -0.4_wp]
      !> Absolute acceptance bound, see the error budget above
      real(wp), parameter :: conv_tol = 3.0e-10_wp
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      real(wp), allocatable :: fa_r(:, :), fb_r(:, :), g_r(:, :)
      complex(wp), allocatable :: fa_k(:, :), fb_k(:, :)
      real(wp) :: r1(3), cab(3), k(3), ab, peak, ph, ref
      integer :: i, j

      call gaussian_grid(grid, trafo, error)
      if (allocated(error)) return
      r1 = grid%get_origin() + 0.5_wp*grid%dr
      cab = gauss_c + conv_c
      ab = gauss_a*conv_b/(gauss_a + conv_b)
      peak = (pi/(gauss_a + conv_b))**1.5_wp

      allocate (fa_r(grid%ngrid, 1), fb_r(grid%ngrid, 1), g_r(grid%ngrid, 1))
      allocate (fa_k(grid%npts_k, 1), fb_k(grid%npts_k, 1))
      do i = 1, grid%ngrid
         fa_r(i, 1) = gaussian(grid%xyz(:, i))
         fb_r(i, 1) = gaussian_at(grid%xyz(:, i), conv_c, conv_b)
      end do
      call forward(error, trafo, fa_r, fa_k)
      if (allocated(error)) return
      call forward(error, trafo, fb_r, fb_k)
      if (allocated(error)) return

      ! Drop the second exp(i k.r1) the product carries
      do j = 1, grid%npts_k
         k = grid%kpoint(j)
         ph = -dot_product(k, r1)
         fa_k(j, 1) = fa_k(j, 1)*fb_k(j, 1)*cmplx(cos(ph), sin(ph), wp)
      end do
      call backward(error, trafo, fa_k, g_r)
      if (allocated(error)) return

      do i = 1, grid%ngrid
         ref = peak*gaussian_at(grid%xyz(:, i), cab, ab)
         call check(error, ieee_is_finite(ref), "analytic convolution reference is not finite")
         if (allocated(error)) return
         call check(error, g_r(i, 1), ref, thr=conv_tol, &
            & more="FFT convolution of two Gaussians deviates from the analytic convolution")
         if (allocated(error)) return
      end do

      call trafo%destroy()
      call grid%destroy()
   end subroutine test_fft_convolution

end module test_math_grid_3d_fft
