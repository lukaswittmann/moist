!> Test suite verifying the vendored FINUFFT (non-uniform FFT) integration
!>
!> Exercises the bundled FINUFFT C++ library through its Fortran wrapper and
!> the finufft_mod options module, confirming the full build/link chain works;
!> every output mode is checked against the direct sum
module test_finufft
   use mctc_env, only: wp
   use iso_fortran_env, only: int64
   use iso_c_binding, only: c_int, c_int64_t, c_double, c_double_complex, c_ptr, c_null_ptr
   use finufft_mod, only: finufft_opts
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_finufft

   !> Number of nonuniform source points used by the transforms
   integer(int64), parameter :: npts = 1000_int64
   !> Number of output Fourier modes
   integer(int64), parameter :: nmodes = 100_int64
   !> Requested FINUFFT tolerance
   real(wp), parameter :: tol = 1.0e-9_wp
   !> Acceptance threshold for the all-mode error, relative to the largest direct-sum mode
   !>
   !> Measured 2.2e-10 with default options and 7.3e-10 with upsampfac = 1.25
   real(wp), parameter :: thr = 5.0e-9_wp
   !> pi = acos(-1.0_wp)
   real(wp), parameter :: pi = 3.14159265358979323846_wp

   !> FINUFFT entry points with the options argument as a C pointer, so a null
   !> pointer (default options) is passed without an unassociated Fortran pointer
   interface
      subroutine finufft1d1_default(nj, xj, cj, iflag, eps, ms, fk, opts, ier) &
         & bind(c, name="finufft1d1_")
         import :: c_int, c_int64_t, c_double, c_double_complex, c_ptr
         integer(c_int64_t), intent(in) :: nj
         real(c_double), intent(in) :: xj(*)
         complex(c_double_complex), intent(in) :: cj(*)
         integer(c_int), intent(in) :: iflag
         real(c_double), intent(in) :: eps
         integer(c_int64_t), intent(in) :: ms
         complex(c_double_complex), intent(inout) :: fk(*)
         type(c_ptr), value :: opts
         integer(c_int), intent(out) :: ier
      end subroutine finufft1d1_default

      subroutine finufft_makeplan_default(ttype, dim, n_modes, iflag, ntrans, eps, plan, opts, ier) &
         & bind(c, name="finufft_makeplan_")
         import :: c_int, c_int64_t, c_double, c_ptr
         integer(c_int), intent(in) :: ttype
         integer(c_int), intent(in) :: dim
         integer(c_int64_t), intent(in) :: n_modes(*)
         integer(c_int), intent(in) :: iflag
         integer(c_int), intent(in) :: ntrans
         real(c_double), intent(in) :: eps
         integer(c_int64_t), intent(out) :: plan
         type(c_ptr), value :: opts
         integer(c_int), intent(out) :: ier
      end subroutine finufft_makeplan_default
   end interface

contains

   !> Collect all FINUFFT integration tests
   !>
   !> @param[out] testsuite  Collection of tests
   subroutine collect_finufft(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("finufft1d1_default_opts", test_1d1_default), &
         new_unittest("finufft1d1_custom_opts", test_1d1_custom), &
         new_unittest("finufft_makeplan_eps_too_small", test_makeplan_eps_too_small) &
      ]
   end subroutine collect_finufft

   !> 1D type-1 transform with default options, checked against a direct DFT
   !>
   !> @param[out] error  Test failure
   subroutine test_1d1_default(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: xj(:)
      complex(wp), allocatable :: cj(:), fk(:)
      integer :: iflag, ier

      call make_problem(xj, cj, fk, iflag)
      ! Null options pointer selects FINUFFT's default options
      call finufft1d1_default(npts, xj, cj, iflag, tol, nmodes, fk, c_null_ptr, ier)

      call check(error, ier == 0, "finufft1d1 (default opts) returned nonzero status")
      if (allocated(error)) return
      call check_mode(error, xj, cj, fk, iflag)
   end subroutine test_1d1_default

   !> Same transform, but driving the finufft_opts derived type explicitly
   !>
   !> @param[out] error  Test failure
   subroutine test_1d1_custom(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(finufft_opts) :: opts
      external :: finufft1d1, finufft_default_opts

      real(wp), allocatable :: xj(:)
      complex(wp), allocatable :: cj(:), fk(:)
      integer :: iflag, ier

      call make_problem(xj, cj, fk, iflag)

      call finufft_default_opts(opts)
      opts%debug = 0
      opts%upsampfac = 1.25_wp

      call finufft1d1(npts, xj, cj, iflag, tol, nmodes, fk, opts, ier)

      call check(error, ier == 0, "finufft1d1 (custom opts) returned nonzero status")
      if (allocated(error)) return
      call check_mode(error, xj, cj, fk, iflag)
   end subroutine test_1d1_custom

   !> Guru makeplan with an unreasonably small tolerance must fail cleanly
   !>
   !> Must report `FINUFFT_ERR_EPS_TOO_SMALL` instead of aborting the
   !> process; regression test for a GCC/libc++ unwinder clash: linking
   !> libgcc_eh into libmoist for its heap-trampoline symbols also pulled in
   !> GCC's own _Unwind_Resume, which aborted every C++ exception
   !> FINUFFT/ducc0 threw across a frame with cleanup handlers (see
   !> config/meson.build and CMakeLists.txt, the heapt_w fix)
   !>
   !> @param[out] error  Test failure
   subroutine test_makeplan_eps_too_small(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      external :: finufft_destroy

      integer :: ttype, dim, ntrans, iflag, ier
      integer(int64) :: n_modes(3)
      !> Opaque pointer-to-plan, passed by reference like FFTW3's legacy
      !> Fortran interface (see finufft/fortran/finufftfort.cpp)
      integer(int64) :: plan
      real(wp) :: bad_tol

      ttype = 1
      dim = 3
      ntrans = 1
      iflag = 1
      n_modes = [4_int64, 4_int64, 4_int64]
      bad_tol = 1.0e-20_wp

      call finufft_makeplan_default(ttype, dim, n_modes, iflag, ntrans, bad_tol, plan, c_null_ptr, ier)

      call check(error, ier /= 0, &
         & "finufft_makeplan with eps=1e-20 must return a nonzero status, not abort")
      if (ier == 0) call finufft_destroy(plan, ier)
   end subroutine test_makeplan_eps_too_small

   !> Build a reproducible quasi-random nonuniform problem in [-pi, pi)
   !>
   !> @param[out] xj     Nonuniform source coordinates
   !> @param[out] cj     Complex source strengths
   !> @param[out] fk     Output mode coefficients
   !> @param[out] iflag  Sign of the imaginary unit in the transform
   subroutine make_problem(xj, cj, fk, iflag)
      !> Nonuniform source coordinates
      real(wp), allocatable, intent(out) :: xj(:)
      !> Complex source strengths
      complex(wp), allocatable, intent(out) :: cj(:)
      !> Output mode coefficients
      complex(wp), allocatable, intent(out) :: fk(:)
      !> Sign of the imaginary unit in the transform
      integer, intent(out) :: iflag
      integer(int64) :: j

      allocate(xj(npts), cj(npts), fk(nmodes))
      do j = 1, npts
         xj(j) = pi * cos(pi * real(j, wp) / real(npts, wp))
         cj(j) = cmplx(sin(100.0_wp * j / npts), &
            & cos(1.0_wp + 50.0_wp * j / npts), wp)
      end do
      iflag = 1
   end subroutine make_problem

   !> Compare every FINUFFT output mode against its direct summation
   !>
   !> Largest deviation over all `nmodes` modes, relative to the largest
   !> direct-sum mode, so a spurious output entry cannot inflate the scale
   !>
   !> @param[out] error  Test failure
   !> @param[in]  xj     Nonuniform source coordinates
   !> @param[in]  cj     Complex source strengths
   !> @param[in]  fk     Output mode coefficients
   !> @param[in]  iflag  Sign of the imaginary unit in the transform
   subroutine check_mode(error, xj, cj, fk, iflag)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Nonuniform source coordinates
      real(wp), intent(in) :: xj(:)
      !> Complex source strengths
      complex(wp), intent(in) :: cj(:)
      !> Output mode coefficients
      complex(wp), intent(in) :: fk(:)
      !> Sign of the imaginary unit in the transform
      integer, intent(in) :: iflag
      complex(wp) :: fkref
      real(wp) :: fmax, errmax
      integer(int64) :: j, k, kindex

      fmax = 0.0_wp
      errmax = 0.0_wp
      do k = -nmodes / 2, nmodes / 2 - 1
         ! direct evaluation of the mode at frequency k
         fkref = (0.0_wp, 0.0_wp)
         do j = 1, npts
            fkref = fkref + cj(j) * cmplx(cos(k * xj(j)), &
               & sin(iflag * k * xj(j)), wp)
         end do
         ! FINUFFT stores modes as -N/2 .. N/2-1, so frequency k is at this index
         kindex = k + nmodes / 2 + 1
         errmax = max(errmax, abs(fk(kindex) - fkref))
         fmax = max(fmax, abs(fkref))
      end do
      call check(error, fmax > 0.0_wp .and. errmax / fmax < thr, &
         & "FINUFFT 1D type-1 mode error too large")
   end subroutine check_mode

end module test_finufft
