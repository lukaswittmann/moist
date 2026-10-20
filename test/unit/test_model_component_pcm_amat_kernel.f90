!> Kernel-level unit tests for the auto-generated Gaussian PCM pair kernel
module test_model_component_pcm_amat_kernel
   use mctc_env_accuracy, only: wp
   use moist_math_boys, only: boys_x_far, boys_x_taylor, boys012
   use mctc_env_error, only: moist_error => error_type
   use test_pcm_conductor_probe, only: probe_conductor_matrix, probe_conductor_bilinear_weights, &
      & probe_conductor_fixed_seed_tangent
   use test_helpers, only: fd4_scalar, fd4_offsets
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_amat_kernel

   !> Sample points (xi_i, xi_j, r2)
   integer, parameter :: npts = 10
   real(wp), parameter :: pts(3, npts) = reshape([ &
                                                 3.0_wp, 4.0_wp, 1.7361111111111111e-4_wp, &
                                                 3.0_wp, 4.0_wp, 1.0e-2_wp, &
                                                 3.0_wp, 4.0_wp, 8.6805555555555556e-2_wp, &
                                                 3.0_wp, 4.0_wp, 0.2_wp, &
                                                 3.0_wp, 4.0_wp, 1.0_wp, &
                                                 3.0_wp, 4.0_wp, 13.0_wp, &
                                                 3.0_wp, 4.0_wp, 14.0625_wp, &
                                                 3.0_wp, 4.0_wp, 40.0_wp, &
                                                 1.5_wp, 12.0_wp, 0.3_wp, &
                                                 9.0_wp, 11.0_wp, 2.5_wp], [3, npts])

   !> Extreme sample points used for the value and symmetry checks only
   integer, parameter :: nextreme = 4
   real(wp), parameter :: extreme_pts(3, nextreme) = reshape([ &
                                                             3.0_wp, 4.0_wp, 1.0e-11_wp, &
                                                             3.0_wp, 4.0_wp, 1.0e-6_wp, &
                                                             0.5_wp, 0.5_wp, 1.0e-8_wp, &
                                                             3.0_wp, 4.0_wp, 1.0e6_wp], [3, nextreme])

   !> Diagonal sample points (xi, f)
   integer, parameter :: ndiag = 4
   real(wp), parameter :: diag_pts(2, ndiag) = reshape([ &
                                                       3.0_wp, 1.0_wp, &
                                                       7.5_wp, 0.35_wp, &
                                                       0.8_wp, 0.02_wp, &
                                                       12.0_wp, 0.99_wp], [2, ndiag])

   !> Finite-difference steps
   real(wp), parameter :: xi_step = 1.0e-3_wp
   real(wp), parameter :: r2_rel_step = 3.0e-3_wp

   !> 4-point central FD tolerances, matching the cfc_kernel convention
   !> Analytic value checked as `actual` against the finite FD value, thr = atol + rtol*|FD|
   real(wp), parameter :: grad_atol = 1.0e-10_wp
   real(wp), parameter :: grad_rtol = 1.0e-9_wp
   real(wp), parameter :: hess_atol = 1.0e-10_wp
   real(wp), parameter :: hess_rtol = 1.0e-9_wp

   !> Tolerance for identities that must hold to (near) machine precision
   real(wp), parameter :: exact_thr = 1.0e-14_wp

   !> Relative tolerance for identities
   real(wp), parameter :: ulp_thr = 8.0_wp*epsilon(1.0_wp)

contains

   !> Collect the generated-kernel test suite
   subroutine collect_model_component_pcm_amat_kernel(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("near_value_vs_erf", test_near_value_vs_erf), &
                  new_unittest("boys_vs_reference", test_boys_vs_reference), &
                  new_unittest("near_grad_fd", test_near_grad_fd), &
                  new_unittest("near_hess_fd", test_near_hess_fd), &
                  new_unittest("near_exchange_symmetry", test_near_exchange_symmetry), &
                  new_unittest("far_is_coulomb", test_far_is_coulomb), &
                  new_unittest("far_matches_near_at_seam", test_far_matches_near_at_seam), &
                  new_unittest("diag_derivatives_fd", test_diag_derivatives_fd) &
                  ]

   end subroutine collect_model_component_pcm_amat_kernel

   !> Reference kernel value evaluated straight from its closed form
   !>
   !> @param[in]  xi_i  Gaussian width of grid point i
   !> @param[in]  xi_j  Gaussian width of grid point j
   !> @param[in]  r2    Squared separation
   !> @return     a     Reference kernel value erf(p*r)/r
   pure function reference_value(xi_i, xi_j, r2) result(a)
      !> Gaussian widths of the two  grid points
      real(wp), intent(in) :: xi_i, xi_j
      !> Squared separation
      real(wp), intent(in) :: r2
      !> Reference kernel value
      real(wp) :: a

      !> Pair width and separation
      real(wp) :: p, r

      p = xi_i*xi_j/sqrt(xi_i*xi_i + xi_j*xi_j)
      r = sqrt(r2)
      a = erf(p*r)/r

   end function reference_value

   !> Boys argument x = (p*r)**2 of a sample point
   !>
   !> @param[in]  xi_i  Gaussian width of grid point i
   !> @param[in]  xi_j  Gaussian width of grid point j
   !> @param[in]  r2    Squared separation
   !> @return     x     Boys argument
   pure function boys_argument(xi_i, xi_j, r2) result(x)
      !> Gaussian widths of the two  grid points
      real(wp), intent(in) :: xi_i, xi_j
      !> Squared separation
      real(wp), intent(in) :: r2
      !> Boys argument
      real(wp) :: x

      x = r2*(xi_i*xi_i*xi_j*xi_j)/(xi_i*xi_i + xi_j*xi_j)

   end function boys_argument

   !> The full kernel reproduces erf(p*r)/r on every branch
   subroutine test_near_value_vs_erf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample-point index
      integer :: i
      !> Generated and reference kernel values
      real(wp) :: a, a_ref

      do i = 1, npts
         call probe_pair_value(pts(1, i), pts(2, i), pts(3, i), a, error)
         if (allocated(error)) return
         a_ref = reference_value(pts(1, i), pts(2, i), pts(3, i))
         call check(error, a, a_ref, thr=exact_thr*max(1.0_wp, abs(a_ref)), &
                    more="generated kernel value deviates from erf(p*r)/r")
         if (allocated(error)) return
      end do

      do i = 1, nextreme
         call probe_pair_value(extreme_pts(1, i), extreme_pts(2, i), extreme_pts(3, i), a, error)
         if (allocated(error)) return
         a_ref = reference_value(extreme_pts(1, i), extreme_pts(2, i), &
                                 extreme_pts(3, i))
         call check(error, a, a_ref, thr=exact_thr*max(1.0_wp, abs(a_ref)), &
                    more="generated kernel value deviates from erf(p*r)/r")
         if (allocated(error)) return
      end do

   end subroutine test_near_value_vs_erf

   !> The Boys evaluator reproduces F_n(x) = int_0^1 t**(2n) exp(-x t**2) dt
   subroutine test_boys_vs_reference(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample-point index
      integer :: i
      !> Boys argument and generated Boys values
      real(wp) :: x, f0, f1, f2, ex
      !> Reference Boys values
      real(wp) :: f0_ref, f1_ref, f2_ref
      !> Arguments where the downward recursion is unconditionally stable
      real(wp), parameter :: stable_x(5) = &
                             [0.5_wp, 1.0_wp, 4.0_wp, 20.0_wp, 60.0_wp]

      ! Exact values at the origin
      call boys012(0.0_wp, f0, f1, f2, ex)
      call check(error, f0, 1.0_wp, thr=exact_thr, more="F0(0) is wrong")
      if (allocated(error)) return
      call check(error, f1, 1.0_wp/3.0_wp, thr=exact_thr, more="F1(0) is wrong")
      if (allocated(error)) return
      call check(error, f2, 0.2_wp, thr=exact_thr, more="F2(0) is wrong")
      if (allocated(error)) return
      call check(error, ex, 1.0_wp, thr=exact_thr, more="exp(-0) is wrong")
      if (allocated(error)) return

      do i = 1, size(stable_x)
         x = stable_x(i)
         call boys012(x, f0, f1, f2, ex)
         f0_ref = 0.88622692545275801_wp*erf(sqrt(x))/sqrt(x)
         f1_ref = 0.5_wp*(f0_ref - exp(-x))/x
         f2_ref = 0.5_wp*(3.0_wp*f1_ref - exp(-x))/x
         call check(error, f0, f0_ref, thr=exact_thr*max(1.0_wp, abs(f0_ref)), &
                    more="F0 deviates from its closed form")
         if (allocated(error)) return
         call check(error, f1, f1_ref, thr=exact_thr*max(1.0_wp, abs(f1_ref)), &
                    more="F1 deviates from the stable recursion")
         if (allocated(error)) return
         call check(error, f2, f2_ref, thr=exact_thr*max(1.0_wp, abs(f2_ref)), &
                    more="F2 deviates from the stable recursion")
         if (allocated(error)) return
      end do

      ! The Taylor branch must join the recursion branch continuously at the
      ! crossover; a short series would show up here first
      x = boys_x_taylor
      call boys012(x*(1.0_wp - epsilon(1.0_wp)), f0, f1, f2, ex)
      call boys012(x*(1.0_wp + 8.0_wp*epsilon(1.0_wp)), f0_ref, f1_ref, &
                   f2_ref, ex)
      call check(error, f0, f0_ref, thr=1.0e-15_wp, &
                 more="Boys F0 is discontinuous at the Taylor crossover")
      if (allocated(error)) return
      call check(error, f1, f1_ref, thr=1.0e-15_wp, &
                 more="Boys F1 is discontinuous at the Taylor crossover")
      if (allocated(error)) return
      call check(error, f2, f2_ref, thr=1.0e-15_wp, &
                 more="Boys F2 is discontinuous at the Taylor crossover")

   end subroutine test_boys_vs_reference

   !> First derivatives of the full kernel against finite differences
   subroutine test_near_grad_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample-point and stencil indices
      integer :: i, k
      !> Sample point
      real(wp) :: xi_i, xi_j, r2, step
      !> Analytic value and first derivatives
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      !> Stencil buffer
      real(wp) :: vals(4)
      !> Finite-difference derivative
      real(wp) :: fd

      do i = 1, npts
         xi_i = pts(1, i)
         xi_j = pts(2, i)
         r2 = pts(3, i)
         call probe_pair_grad(xi_i, xi_j, r2, a, a_xi_i, a_xi_j, a_r2, error)
         if (allocated(error)) return

         do k = 1, 4
            vals(k) = reference_value(xi_i + fd4_offsets(k)*xi_step, xi_j, r2)
         end do
         call fd4_scalar(vals(1), vals(2), vals(3), vals(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_xi_i, fd, thr=grad_atol + grad_rtol*abs(fd), &
                    more="dA/dxi_i vs FD")
         if (allocated(error)) return

         do k = 1, 4
            vals(k) = reference_value(xi_i, xi_j + fd4_offsets(k)*xi_step, r2)
         end do
         call fd4_scalar(vals(1), vals(2), vals(3), vals(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_xi_j, fd, thr=grad_atol + grad_rtol*abs(fd), &
                    more="dA/dxi_j vs FD")
         if (allocated(error)) return

         step = r2_rel_step*r2
         do k = 1, 4
            vals(k) = reference_value(xi_i, xi_j, r2 + fd4_offsets(k)*step)
         end do
         call fd4_scalar(vals(1), vals(2), vals(3), vals(4), step, fd, error)
         if (allocated(error)) return
         call check(error, a_r2, fd, thr=grad_atol + grad_rtol*abs(fd), &
                    more="dA/dr2 vs FD")
         if (allocated(error)) return
      end do

   end subroutine test_near_grad_fd

   !> Second derivatives against finite differences of the analytic gradient
   subroutine test_near_hess_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample-point and stencil indices
      integer :: i, k
      !> Sample point
      real(wp) :: xi_i, xi_j, r2, step
      !> Analytic value, first and second derivatives
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      real(wp) :: a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2
      !> Stencil buffers for the three first-derivative channels
      real(wp) :: gi(4), gj(4), gr(4)
      !> Finite-difference derivative
      real(wp) :: fd
      !> Unused stencil outputs
      real(wp) :: dummy

      do i = 1, npts
         xi_i = pts(1, i)
         xi_j = pts(2, i)
         r2 = pts(3, i)
         call probe_pair_hess(xi_i, xi_j, r2, a, a_xi_i, a_xi_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2, error)
         if (allocated(error)) return

         do k = 1, 4
            call probe_pair_grad(xi_i + fd4_offsets(k)*xi_step, xi_j, r2, dummy, gi(k), gj(k), gr(k), error)
            if (allocated(error)) return
         end do
         call fd4_scalar(gi(1), gi(2), gi(3), gi(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_ii, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A/dxi_i**2 vs FD")
         if (allocated(error)) return
         call fd4_scalar(gj(1), gj(2), gj(3), gj(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_ij, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A/dxi_i dxi_j vs FD")
         if (allocated(error)) return
         call fd4_scalar(gr(1), gr(2), gr(3), gr(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_ir2, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A/dxi_i dr2 vs FD")
         if (allocated(error)) return

         do k = 1, 4
            call probe_pair_grad(xi_i, xi_j + fd4_offsets(k)*xi_step, r2, dummy, gi(k), gj(k), gr(k), error)
            if (allocated(error)) return
         end do
         call fd4_scalar(gj(1), gj(2), gj(3), gj(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_jj, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A/dxi_j**2 vs FD")
         if (allocated(error)) return
         call fd4_scalar(gr(1), gr(2), gr(3), gr(4), xi_step, fd, error)
         if (allocated(error)) return
         call check(error, a_jr2, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A/dxi_j dr2 vs FD")
         if (allocated(error)) return

         step = r2_rel_step*r2
         do k = 1, 4
            call probe_pair_grad(xi_i, xi_j, r2 + fd4_offsets(k)*step, dummy, gi(k), gj(k), gr(k), error)
            if (allocated(error)) return
         end do
         call fd4_scalar(gr(1), gr(2), gr(3), gr(4), step, fd, error)
         if (allocated(error)) return
         call check(error, a_r2r2, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A/dr2**2 vs FD")
         if (allocated(error)) return
      end do

   end subroutine test_near_hess_fd

   !> Swapping the two  grid points mirrors the width channels and leaves the rest
   subroutine test_near_exchange_symmetry(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample-point index
      integer :: i
      !> Kernel data in the original and swapped orderings
      real(wp) :: a, a_i, a_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2
      real(wp) :: b, b_i, b_j, b_r2, b_ii, b_ij, b_jj, b_ir2, b_jr2, b_r2r2

      do i = 1, npts
         call probe_pair_hess(pts(1, i), pts(2, i), pts(3, i), a, a_i, a_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, &
            & a_r2r2, error)
         if (allocated(error)) return
         call probe_pair_hess(pts(2, i), pts(1, i), pts(3, i), b, b_i, b_j, b_r2, b_ii, b_ij, b_jj, b_ir2, b_jr2, &
            & b_r2r2, error)
         if (allocated(error)) return

         call check(error, a, b, thr=exact_thr*max(1.0_wp, abs(a)), &
                    more="kernel value is not symmetric under i <-> j")
         if (allocated(error)) return
         call check(error, a_i, b_j, thr=exact_thr*max(1.0_wp, abs(a_i)), &
                    more="dA/dxi_i is not the mirror of dA/dxi_j")
         if (allocated(error)) return
         call check(error, a_ii, b_jj, thr=exact_thr*max(1.0_wp, abs(a_ii)), &
                    more="d2A/dxi_i**2 is not the mirror of d2A/dxi_j**2")
         if (allocated(error)) return
         call check(error, a_ij, b_ij, thr=exact_thr*max(1.0_wp, abs(a_ij)), &
                    more="the mixed width Hessian is not symmetric")
         if (allocated(error)) return
         call check(error, a_ir2, b_jr2, thr=exact_thr*max(1.0_wp, abs(a_ir2)), &
                    more="the mixed width/separation channels do not mirror")
         if (allocated(error)) return
         call check(error, a_r2r2, b_r2r2, thr=exact_thr*max(1.0_wp, abs(a_r2r2)), &
                    more="d2A/dr2**2 is not symmetric under i <-> j")
         if (allocated(error)) return
      end do

   end subroutine test_near_exchange_symmetry

   !> The saturated branch is exactly the bare Coulomb kernel
   subroutine test_far_is_coulomb(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample index
      integer :: i
      !> Squared separation and Boys argument
      real(wp) :: r2
      !> Saturated-branch and full-kernel data
      real(wp) :: a_far, a_r2_far, a_r2r2_far
      real(wp) :: a, a_xi_i, a_xi_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2
      !> Widths whose pair width is 1 by construction (p = 1/sqrt(2)*xi)
      real(wp), parameter :: xi = 3.0_wp
      !> Separations placing x = (p*r)**2 comfortably past the cutoff
      real(wp), parameter :: factors(4) = [1.0_wp, 2.0_wp, 10.0_wp, 100.0_wp]

      do i = 1, size(factors)
         ! p**2 = xi**2/2, so x = p**2*r2; pick r2 from the requested x
         r2 = factors(i)*boys_x_far*2.0_wp/(xi*xi)

         call coulomb_reference(r2, a_far, a_r2_far, a_r2r2_far)
         call probe_pair_hess(xi, xi, r2, a, a_xi_i, a_xi_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2, error)
         if (allocated(error)) return

         ! Compare the generated full kernel to the analytic Coulomb limit
         call check(error, a, a_far, thr=ulp_thr*abs(a_far), &
                    more="full kernel does not reduce to 1/r past the cutoff")
         if (allocated(error)) return
         call check(error, a_r2, a_r2_far, thr=ulp_thr*abs(a_r2_far), &
                    more="full dA/dr2 does not reduce to the Coulomb form")
         if (allocated(error)) return
         call check(error, a_r2r2, a_r2r2_far, thr=ulp_thr*abs(a_r2r2_far), &
                    more="full d2A/dr2**2 does not reduce to the Coulomb form")
         if (allocated(error)) return

         call check(error, a_xi_i, 0.0_wp, thr=0.0_wp, &
                    more="dA/dxi_i is not identically zero past the cutoff")
         if (allocated(error)) return
         call check(error, a_xi_j, 0.0_wp, thr=0.0_wp, &
                    more="dA/dxi_j is not identically zero past the cutoff")
         if (allocated(error)) return
         call check(error, abs(a_ii) + abs(a_ij) + abs(a_jj), 0.0_wp, thr=0.0_wp, &
                    more="the width Hessian is not identically zero past the cutoff")
         if (allocated(error)) return
         call check(error, abs(a_ir2) + abs(a_jr2), 0.0_wp, thr=0.0_wp, &
                    more="the mixed channels are not identically zero past the cutoff")
         if (allocated(error)) return
      end do

   end subroutine test_far_is_coulomb

   !> The two branches agree to machine precision immediately below the seam
   subroutine test_far_matches_near_at_seam(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Squared separation just inside the near branch
      real(wp) :: r2
      !> Saturated-branch and full-kernel data
      real(wp) :: a_far, a_r2_far
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      !> Widths whose pair width squared is xi**2/2
      real(wp), parameter :: xi = 3.0_wp

      ! Step just below the cutoff so the erf/exp branch is taken
      r2 = (1.0_wp - 1.0e-12_wp)*boys_x_far*2.0_wp/(xi*xi)
      call coulomb_reference(r2, a_far, a_r2_far)
      call probe_pair_grad(xi, xi, r2, a, a_xi_i, a_xi_j, a_r2, error)
      if (allocated(error)) return

      call check(error, a, a_far, thr=ulp_thr*abs(a_far), &
                 more="value is discontinuous across the saturation cutoff")
      if (allocated(error)) return
      call check(error, a_r2, a_r2_far, thr=ulp_thr*abs(a_r2_far), &
                 more="dA/dr2 is discontinuous across the saturation cutoff")
      if (allocated(error)) return
      call check(error, abs(a_xi_i) + abs(a_xi_j), 0.0_wp, thr=1.0e-30_wp, &
                 more="width channels are not negligible at the cutoff")

   end subroutine test_far_matches_near_at_seam

   !> Diagonal self-interaction derivatives against finite differences
   subroutine test_diag_derivatives_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sample-point and stencil indices
      integer :: i, k
      !> Sample point and step sizes
      real(wp) :: xi, f, step_xi, step_f
      !> Analytic value and derivatives
      real(wp) :: a, a_xi, a_f, a_xi_xi, a_xi_f, a_f_f
      !> Stencil buffers
      real(wp) :: vals(4), gxi(4), gf(4)
      !> Finite-difference derivative
      real(wp) :: fd
      !> Unused stencil output
      real(wp) :: dummy

      do i = 1, ndiag
         xi = diag_pts(1, i)
         f = diag_pts(2, i)
         step_xi = 1.0e-4_wp*xi
         step_f = 1.0e-4_wp*f
         call probe_diag_hess(xi, f, a, a_xi, a_f, a_xi_xi, a_xi_f, a_f_f, error)
         if (allocated(error)) return

         call check(error, a, sqrt(2.0_wp/(4.0_wp*atan(1.0_wp)))*xi/f, &
                    thr=exact_thr*max(1.0_wp, abs(a)), &
                    more="diagonal value deviates from sqrt(2/pi)*xi/f")
         if (allocated(error)) return

         do k = 1, 4
            call probe_diag_grad(xi + fd4_offsets(k)*step_xi, f, vals(k), gxi(k), gf(k), error)
            if (allocated(error)) return
         end do
         call fd4_scalar(vals(1), vals(2), vals(3), vals(4), step_xi, fd, error)
         if (allocated(error)) return
         call check(error, a_xi, fd, thr=grad_atol + grad_rtol*abs(fd), &
                    more="dA_ii/dxi vs FD")
         if (allocated(error)) return
         call fd4_scalar(gxi(1), gxi(2), gxi(3), gxi(4), step_xi, fd, error)
         if (allocated(error)) return
         call check(error, a_xi_xi, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A_ii/dxi**2 vs FD")
         if (allocated(error)) return
         call fd4_scalar(gf(1), gf(2), gf(3), gf(4), step_xi, fd, error)
         if (allocated(error)) return
         call check(error, a_xi_f, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A_ii/dxi df vs FD")
         if (allocated(error)) return

         do k = 1, 4
            call probe_diag_grad(xi, f + fd4_offsets(k)*step_f, vals(k), gxi(k), gf(k), error)
            if (allocated(error)) return
         end do
         call fd4_scalar(vals(1), vals(2), vals(3), vals(4), step_f, fd, error)
         if (allocated(error)) return
         call check(error, a_f, fd, thr=grad_atol + grad_rtol*abs(fd), &
                    more="dA_ii/df vs FD")
         if (allocated(error)) return
         call fd4_scalar(gf(1), gf(2), gf(3), gf(4), step_f, fd, error)
         if (allocated(error)) return
         call check(error, a_f_f, fd, thr=hess_atol + hess_rtol*abs(fd), &
                    more="d2A_ii/df**2 vs FD")
         if (allocated(error)) return

         call probe_diag_value(xi, f, dummy, error)
         if (allocated(error)) return
         call check(error, dummy, a, thr=ulp_thr*abs(a), &
                    more="diagonal value routine disagrees with the Hessian routine")
         if (allocated(error)) return
      end do

   end subroutine test_diag_derivatives_fd

   !> Read a two-point entry through the generated operator
   !>
   !> @param[in] xi_i First Gaussian width, bohr^-1
   !> @param[in] xi_j Second Gaussian width, bohr^-1
   !> @param[in] r2 Squared separation, bohr^2
   !> @param[out] a Matrix entry
   !> @param[out] error Test failure
   subroutine probe_pair_value(xi_i, xi_j, r2, a, error)
      !> First Gaussian width, bohr^-1
      real(wp), intent(in) :: xi_i
      !> Second Gaussian width, bohr^-1
      real(wp), intent(in) :: xi_j
      !> Squared separation, bohr^2
      real(wp), intent(in) :: r2
      !> Matrix entry
      real(wp), intent(out) :: a
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp) :: xyz(3, 2), matrix(2, 2)

      xyz = 0.0_wp
      xyz(1:2, 2) = sqrt(0.5_wp*r2)
      call probe_conductor_matrix([xi_i, xi_j], [1.0_wp, 1.0_wp], xyz, matrix, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      a = matrix(1, 2)

   end subroutine probe_pair_value

   !> Read generated pair derivatives in width and squared-distance coordinates
   !>
   !> @param[in] xi_i First Gaussian width, bohr^-1
   !> @param[in] xi_j Second Gaussian width, bohr^-1
   !> @param[in] r2 Squared separation, bohr^2
   !> @param[out] a Matrix entry
   !> @param[out] a_i First width derivative
   !> @param[out] a_j Second width derivative
   !> @param[out] a_r2 Squared-separation derivative
   !> @param[out] error Test failure
   subroutine probe_pair_grad(xi_i, xi_j, r2, a, a_i, a_j, a_r2, error)
      !> First Gaussian width, bohr^-1
      real(wp), intent(in) :: xi_i
      !> Second Gaussian width, bohr^-1
      real(wp), intent(in) :: xi_j
      !> Squared separation, bohr^2
      real(wp), intent(in) :: r2
      !> Matrix entry
      real(wp), intent(out) :: a
      !> First width derivative
      real(wp), intent(out) :: a_i
      !> Second width derivative
      real(wp), intent(out) :: a_j
      !> Squared-separation derivative
      real(wp), intent(out) :: a_r2
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp) :: xyz(3, 2)
      real(wp) :: w_xi(2), w_f(2), w_r(3, 2)

      xyz = 0.0_wp
      xyz(1:2, 2) = sqrt(0.5_wp*r2)
      call probe_pair_value(xi_i, xi_j, r2, a, error)
      if (allocated(error)) return
      call probe_conductor_bilinear_weights([xi_i, xi_j], [1.0_wp, 1.0_wp], xyz, &
         & [1.0_wp, 0.0_wp], [0.0_wp, 1.0_wp], w_xi, w_f, w_r, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      a_i = w_xi(1)
      a_j = w_xi(2)
      a_r2 = w_r(1, 2)/(2.0_wp*xyz(1, 2))

   end subroutine probe_pair_grad

   !> Read generated pair Hessians without subtracting diagonal terms
   !>
   !> With f = 1, diagonal terms are linear in xi and independent of position.
   !> Their Hessians vanish, so seed = [1, 1] isolates the pair Hessian directly.
   !> For A(r2), the Cartesian cross derivative A_xy = 4*x*y*A_r2r2 avoids
   !> cancellation between radial terms near coincident points.
   !>
   !> @param[in] xi_i First Gaussian width, bohr^-1
   !> @param[in] xi_j Second Gaussian width, bohr^-1
   !> @param[in] r2 Squared separation, bohr^2
   !> @param[out] a Matrix entry
   !> @param[out] a_i First width derivative
   !> @param[out] a_j Second width derivative
   !> @param[out] a_r2 Squared-separation derivative
   !> @param[out] a_ii First width Hessian entry
   !> @param[out] a_ij Mixed width Hessian entry
   !> @param[out] a_jj Second width Hessian entry
   !> @param[out] a_ir2 First width/separation Hessian entry
   !> @param[out] a_jr2 Second width/separation Hessian entry
   !> @param[out] a_r2r2 Squared-separation Hessian entry
   !> @param[out] error Test failure
   subroutine probe_pair_hess(xi_i, xi_j, r2, a, a_i, a_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2, error)
      !> First Gaussian width, bohr^-1
      real(wp), intent(in) :: xi_i
      !> Second Gaussian width, bohr^-1
      real(wp), intent(in) :: xi_j
      !> Squared separation, bohr^2
      real(wp), intent(in) :: r2
      !> Matrix entry
      real(wp), intent(out) :: a
      !> First width derivative
      real(wp), intent(out) :: a_i
      !> Second width derivative
      real(wp), intent(out) :: a_j
      !> Squared-separation derivative
      real(wp), intent(out) :: a_r2
      !> First width Hessian entry
      real(wp), intent(out) :: a_ii
      !> Mixed width Hessian entry
      real(wp), intent(out) :: a_ij
      !> Second width Hessian entry
      real(wp), intent(out) :: a_jj
      !> First width/separation Hessian entry
      real(wp), intent(out) :: a_ir2
      !> Second width/separation Hessian entry
      real(wp), intent(out) :: a_jr2
      !> Squared-separation Hessian entry
      real(wp), intent(out) :: a_r2r2
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp) :: xyz(3, 2)
      real(wp) :: dr(3, 2, 3), dxi(2, 3), df(2, 3)
      real(wp) :: h_r(3, 2, 3), h_xi(2, 3), h_f(2, 3)

      xyz = 0.0_wp
      xyz(1:2, 2) = sqrt(0.5_wp*r2)
      call probe_pair_grad(xi_i, xi_j, r2, a, a_i, a_j, a_r2, error)
      if (allocated(error)) return
      dr = 0.0_wp
      dxi = 0.0_wp
      df = 0.0_wp
      dxi(1, 1) = 1.0_wp
      dxi(2, 2) = 1.0_wp
      dr(1, 2, 3) = 1.0_wp
      call probe_conductor_fixed_seed_tangent([xi_i, xi_j], [1.0_wp, 1.0_wp], xyz, &
         & [1.0_wp, 1.0_wp], dr, dxi, df, h_r, h_xi, h_f, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      a_ii = h_xi(1, 1)
      a_ij = h_xi(1, 2)
      a_jj = h_xi(2, 2)
      a_ir2 = h_xi(1, 3)/(2.0_wp*xyz(1, 2))
      a_jr2 = h_xi(2, 3)/(2.0_wp*xyz(1, 2))
      a_r2r2 = h_r(2, 2, 3)/(4.0_wp*xyz(1, 2)*xyz(2, 2))

   end subroutine probe_pair_hess

   !> Read a diagonal entry through the generated operator
   !>
   !> @param[in] xi Gaussian width, bohr^-1
   !> @param[in] f Switching factor
   !> @param[out] a Matrix entry
   !> @param[out] error Test failure
   subroutine probe_diag_value(xi, f, a, error)
      !> Gaussian width, bohr^-1
      real(wp), intent(in) :: xi
      !> Switching factor
      real(wp), intent(in) :: f
      !> Matrix entry
      real(wp), intent(out) :: a
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp) :: xyz(3, 1), matrix(1, 1)

      xyz = 0.0_wp
      call probe_conductor_matrix([xi], [f], xyz, matrix, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      a = matrix(1, 1)

   end subroutine probe_diag_value

   !> Read generated diagonal width and switching derivatives
   !>
   !> @param[in] xi Gaussian width, bohr^-1
   !> @param[in] f Switching factor
   !> @param[out] a Matrix entry
   !> @param[out] a_xi Width derivative
   !> @param[out] a_f Switching derivative
   !> @param[out] error Test failure
   subroutine probe_diag_grad(xi, f, a, a_xi, a_f, error)
      !> Gaussian width, bohr^-1
      real(wp), intent(in) :: xi
      !> Switching factor
      real(wp), intent(in) :: f
      !> Matrix entry
      real(wp), intent(out) :: a
      !> Width derivative
      real(wp), intent(out) :: a_xi
      !> Switching derivative
      real(wp), intent(out) :: a_f
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp) :: xyz(3, 1)
      real(wp) :: w_xi(1), w_f(1), w_r(3, 1)

      xyz = 0.0_wp
      call probe_diag_value(xi, f, a, error)
      if (allocated(error)) return
      call probe_conductor_bilinear_weights([xi], [f], xyz, [1.0_wp], [1.0_wp], w_xi, w_f, w_r, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      a_xi = w_xi(1)
      a_f = w_f(1)

   end subroutine probe_diag_grad

   !> Read generated diagonal Hessians at a fixed unit seed
   !>
   !> @param[in] xi Gaussian width, bohr^-1
   !> @param[in] f Switching factor
   !> @param[out] a Matrix entry
   !> @param[out] a_xi Width derivative
   !> @param[out] a_f Switching derivative
   !> @param[out] a_xi_xi Width Hessian entry
   !> @param[out] a_xi_f Mixed width/switching Hessian entry
   !> @param[out] a_f_f Switching Hessian entry
   !> @param[out] error Test failure
   subroutine probe_diag_hess(xi, f, a, a_xi, a_f, a_xi_xi, a_xi_f, a_f_f, error)
      !> Gaussian width, bohr^-1
      real(wp), intent(in) :: xi
      !> Switching factor
      real(wp), intent(in) :: f
      !> Matrix entry
      real(wp), intent(out) :: a
      !> Width derivative
      real(wp), intent(out) :: a_xi
      !> Switching derivative
      real(wp), intent(out) :: a_f
      !> Width Hessian entry
      real(wp), intent(out) :: a_xi_xi
      !> Mixed width/switching Hessian entry
      real(wp), intent(out) :: a_xi_f
      !> Switching Hessian entry
      real(wp), intent(out) :: a_f_f
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp) :: xyz(3, 1)
      real(wp) :: dr(3, 1, 2), dxi(1, 2), df(1, 2)
      real(wp) :: h_r(3, 1, 2), h_xi(1, 2), h_f(1, 2)

      xyz = 0.0_wp
      call probe_diag_grad(xi, f, a, a_xi, a_f, error)
      if (allocated(error)) return
      dr = 0.0_wp
      dxi(1, :) = [1.0_wp, 0.0_wp]
      df(1, :) = [0.0_wp, 1.0_wp]
      call probe_conductor_fixed_seed_tangent([xi], [f], xyz, [1.0_wp], dr, dxi, df, h_r, h_xi, h_f, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      a_xi_xi = 2.0_wp*h_xi(1, 1)
      a_xi_f = 2.0_wp*h_xi(1, 2)
      a_f_f = 2.0_wp*h_f(1, 2)

   end subroutine probe_diag_hess

   !> Analytic Coulomb limit and squared-separation derivatives
   !>
   !> @param[in] r2 Positive squared separation
   !> @param[out] a Coulomb value
   !> @param[out] a_r First squared-separation derivative
   !> @param[out] a_rr Second squared-separation derivative
   subroutine coulomb_reference(r2, a, a_r, a_rr)
      !> Positive squared separation
      real(wp), intent(in) :: r2
      !> Coulomb value
      real(wp), intent(out) :: a
      !> First squared-separation derivative
      real(wp), intent(out) :: a_r
      !> Second squared-separation derivative
      real(wp), intent(out), optional :: a_rr

      a = 1.0_wp/sqrt(r2)
      a_r = -0.5_wp*a/r2
      if (present(a_rr)) a_rr = 0.75_wp*a/r2**2

   end subroutine coulomb_reference

end module test_model_component_pcm_amat_kernel
