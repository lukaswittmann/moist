!> Numerical checks of HNC, KH and PSE2-4 kernels
! SPDX-License-Identifier: MPL-2.0
!>
!> Checks closure sweeps and their pointwise derivatives:
!>   * dc/dt and d2c/dt2 match a 4-point (Richardson) central
!>     difference of c and of dc respectively,
!>   * values agree with independent exponential and polynomial references,
!>   * Phi(h, c) stays finite over the sampled domain
!>
!> Sample points stay away from the branch boundaries so each
!> finite-difference stencil stays within one branch
module test_moz_closure_kernel
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_model_moz_closure_kernel, only: &
      & closure_hnc_c, &
      & closure_hnc_dc, &
      & closure_hnc_d2c, &
      & closure_hnc_phi, &
      & closure_kh_c, &
      & closure_kh_dc, &
      & closure_kh_d2c, &
      & closure_kh_phi, &
      & closure_pse2_c, &
      & closure_pse2_dc, &
      & closure_pse2_d2c, &
      & closure_pse2_phi, &
      & closure_pse3_c, &
      & closure_pse3_dc, &
      & closure_pse3_d2c, &
      & closure_pse3_phi, &
      & closure_pse4_c, &
      & closure_pse4_dc, &
      & closure_pse4_d2c, &
      & closure_pse4_phi, &
      & closure_dc_sweep, &
      & closure_d2c_sweep, &
      & closure_c_dc_sweep
   implicit none(type, external)
   private

   public :: collect_moz_closure_kernel

   real(wp), parameter :: beta = 1.3_wp
   real(wp), parameter :: lambda = 0.7_wp
   !> Near-optimal step for the 4-point formula: truncation ~h^4 and
   !> roundoff ~eps/h both land well below fd_tol here
   real(wp), parameter :: fd_h = 1.0e-3_wp
   real(wp), parameter :: fd_tol = 1.0e-10_wp
   !> Kernel and independent reference formulas can differ through rounding
   real(wp), parameter :: id_tol = 1.0e-13_wp

   integer, parameter :: npt = 33
   real(wp), parameter :: pu(33) = [ &
         -1.600000_wp, -1.600000_wp, -1.600000_wp, -1.600000_wp, &
         -1.600000_wp, -0.900000_wp, -0.900000_wp, -0.900000_wp, &
         -0.900000_wp, -0.900000_wp, -0.350000_wp, -0.350000_wp, &
         -0.350000_wp, -0.350000_wp, -0.350000_wp, -0.350000_wp, &
         0.250000_wp, 0.250000_wp, 0.250000_wp, 0.250000_wp, &
         0.250000_wp, 0.250000_wp, 0.800000_wp, 0.800000_wp, &
         0.800000_wp, 0.800000_wp, 0.800000_wp, 0.800000_wp, &
         1.500000_wp, 1.500000_wp, 1.500000_wp, 1.500000_wp, &
         1.500000_wp]
   real(wp), parameter :: pt(33) = [ &
         -2.000000_wp, -1.150000_wp, -0.400000_wp, 0.450000_wp, &
         1.100000_wp, -2.000000_wp, -1.150000_wp, -0.400000_wp, &
         0.450000_wp, 1.100000_wp, -2.000000_wp, -1.150000_wp, &
         -0.400000_wp, 0.450000_wp, 1.100000_wp, 2.200000_wp, &
         -2.000000_wp, -1.150000_wp, -0.400000_wp, 0.450000_wp, &
         1.100000_wp, 2.200000_wp, -2.000000_wp, -1.150000_wp, &
         -0.400000_wp, 0.450000_wp, 1.100000_wp, 2.200000_wp, &
         -1.150000_wp, -0.400000_wp, 0.450000_wp, 1.100000_wp, &
         2.200000_wp]
   integer, parameter :: nph = 20
   real(wp), parameter :: ph(20) = [ &
         -0.850000_wp, -0.850000_wp, -0.850000_wp, -0.850000_wp, &
         -0.300000_wp, -0.300000_wp, -0.300000_wp, -0.300000_wp, &
         0.400000_wp, 0.400000_wp, 0.400000_wp, 0.400000_wp, &
         1.600000_wp, 1.600000_wp, 1.600000_wp, 1.600000_wp, &
         3.000000_wp, 3.000000_wp, 3.000000_wp, 3.000000_wp]
   real(wp), parameter :: pc(20) = [ &
         -2.500000_wp, -0.600000_wp, 0.200000_wp, 1.400000_wp, &
         -2.500000_wp, -0.600000_wp, 0.200000_wp, 1.400000_wp, &
         -2.500000_wp, -0.600000_wp, 0.200000_wp, 1.400000_wp, &
         -2.500000_wp, -0.600000_wp, 0.200000_wp, 1.400000_wp, &
         -2.500000_wp, -0.600000_wp, 0.200000_wp, 1.400000_wp]

contains

   !> Collect the closure-kernel tests
   !>
   !> @param[out] testsuite Collected tests
   subroutine collect_moz_closure_kernel(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
               new_unittest("hnc_derivatives", test_hnc_derivatives), &
               new_unittest("hnc_consistency", test_hnc_consistency), &
               new_unittest("kh_derivatives", test_kh_derivatives), &
               new_unittest("kh_consistency", test_kh_consistency), &
               new_unittest("pse2_derivatives", test_pse2_derivatives), &
               new_unittest("pse2_consistency", test_pse2_consistency), &
               new_unittest("pse3_derivatives", test_pse3_derivatives), &
               new_unittest("pse3_consistency", test_pse3_consistency), &
               new_unittest("pse4_derivatives", test_pse4_derivatives), &
               new_unittest("pse4_consistency", test_pse4_consistency), &
               new_unittest("sweep_references", test_sweeps), &
               new_unittest("sweep_sections", test_sweep_sections), &
               new_unittest("sweep_invalid_kind", test_sweep_invalid_kind), &
               new_unittest("fast_exp_accuracy", test_fast_exp)]
   end subroutine collect_moz_closure_kernel

   !> HNC: dc/dt and d2c/dt2 against a 4-point central
   !> difference of c and of dc, at physical input values
   !>
   !> @param[out] error Test failure
   subroutine test_hnc_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp) :: hs, fd, ana, fp2, fp1, fm1, fm2, values(1)

      do i = 1, npt
         hs = fd_h*max(1.0_wp, abs(pt(i)))
         call closure_hnc_c([pu(i)], [pt(i) + 2.0_wp*hs], beta, lambda, values)
         fp2 = values(1)
         call closure_hnc_c([pu(i)], [pt(i) + hs], beta, lambda, values)
         fp1 = values(1)
         call closure_hnc_c([pu(i)], [pt(i) - hs], beta, lambda, values)
         fm1 = values(1)
         call closure_hnc_c([pu(i)], [pt(i) - 2.0_wp*hs], beta, lambda, values)
         fm2 = values(1)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_hnc_dc(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "HNC dc/dt vs FD")
         if (allocated(error)) return
         fp2 = closure_hnc_dc(pu(i), pt(i) + 2.0_wp*hs, beta, lambda)
         fp1 = closure_hnc_dc(pu(i), pt(i) + hs, beta, lambda)
         fm1 = closure_hnc_dc(pu(i), pt(i) - hs, beta, lambda)
         fm2 = closure_hnc_dc(pu(i), pt(i) - 2.0_wp*hs, beta, lambda)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_hnc_d2c(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "HNC d2c/dt2 vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_hnc_derivatives

   !> HNC: sweeps agree with the closure formula and Phi is finite
   !>
   !> @param[out] error Test failure
   subroutine test_hnc_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: values(npt), reference
      integer :: i

      call closure_hnc_c(pu, pt, beta, lambda, values)
      do i = 1, npt
         reference = reference_c(0, pu(i), pt(i), beta, lambda)
         call check(error, dev(values(i), reference) < id_tol, "HNC closure reference")
         if (allocated(error)) return
      end do
      do i = 1, nph
         reference = closure_hnc_phi(ph(i), pc(i))
         call check(error, abs(reference) < huge(1.0_wp), "HNC Phi finite")
         if (allocated(error)) return
      end do
   end subroutine test_hnc_consistency

   !> KH: dc/dt and d2c/dt2 against a 4-point central
   !> difference of c and of dc, at physical input values
   !>
   !> @param[out] error Test failure
   subroutine test_kh_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp) :: hs, fd, ana, fp2, fp1, fm1, fm2, values(1)

      do i = 1, npt
         hs = fd_h*max(1.0_wp, abs(pt(i)))
         call closure_kh_c([pu(i)], [pt(i) + 2.0_wp*hs], beta, lambda, values)
         fp2 = values(1)
         call closure_kh_c([pu(i)], [pt(i) + hs], beta, lambda, values)
         fp1 = values(1)
         call closure_kh_c([pu(i)], [pt(i) - hs], beta, lambda, values)
         fm1 = values(1)
         call closure_kh_c([pu(i)], [pt(i) - 2.0_wp*hs], beta, lambda, values)
         fm2 = values(1)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_kh_dc(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "KH dc/dt vs FD")
         if (allocated(error)) return
         fp2 = closure_kh_dc(pu(i), pt(i) + 2.0_wp*hs, beta, lambda)
         fp1 = closure_kh_dc(pu(i), pt(i) + hs, beta, lambda)
         fm1 = closure_kh_dc(pu(i), pt(i) - hs, beta, lambda)
         fm2 = closure_kh_dc(pu(i), pt(i) - 2.0_wp*hs, beta, lambda)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_kh_d2c(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "KH d2c/dt2 vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_kh_derivatives

   !> KH: sweeps agree with the closure formula and Phi is finite
   !>
   !> @param[out] error Test failure
   subroutine test_kh_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: values(npt), reference
      integer :: i

      call closure_kh_c(pu, pt, beta, lambda, values)
      do i = 1, npt
         reference = reference_c(1, pu(i), pt(i), beta, lambda)
         call check(error, dev(values(i), reference) < id_tol, "KH closure reference")
         if (allocated(error)) return
      end do
      do i = 1, nph
         reference = closure_kh_phi(ph(i), pc(i))
         call check(error, abs(reference) < huge(1.0_wp), "KH Phi finite")
         if (allocated(error)) return
      end do
   end subroutine test_kh_consistency

   !> PSE2: dc/dt and d2c/dt2 against a 4-point central
   !> difference of c and of dc, at physical input values
   !>
   !> @param[out] error Test failure
   subroutine test_pse2_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp) :: hs, fd, ana, fp2, fp1, fm1, fm2, values(1)

      do i = 1, npt
         hs = fd_h*max(1.0_wp, abs(pt(i)))
         call closure_pse2_c([pu(i)], [pt(i) + 2.0_wp*hs], beta, lambda, values)
         fp2 = values(1)
         call closure_pse2_c([pu(i)], [pt(i) + hs], beta, lambda, values)
         fp1 = values(1)
         call closure_pse2_c([pu(i)], [pt(i) - hs], beta, lambda, values)
         fm1 = values(1)
         call closure_pse2_c([pu(i)], [pt(i) - 2.0_wp*hs], beta, lambda, values)
         fm2 = values(1)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_pse2_dc(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "PSE2 dc/dt vs FD")
         if (allocated(error)) return
         fp2 = closure_pse2_dc(pu(i), pt(i) + 2.0_wp*hs, beta, lambda)
         fp1 = closure_pse2_dc(pu(i), pt(i) + hs, beta, lambda)
         fm1 = closure_pse2_dc(pu(i), pt(i) - hs, beta, lambda)
         fm2 = closure_pse2_dc(pu(i), pt(i) - 2.0_wp*hs, beta, lambda)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_pse2_d2c(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "PSE2 d2c/dt2 vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_pse2_derivatives

   !> PSE2: sweeps agree with the closure formula and Phi is finite
   !>
   !> @param[out] error Test failure
   subroutine test_pse2_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: values(npt), reference
      integer :: i

      call closure_pse2_c(pu, pt, beta, lambda, values)
      do i = 1, npt
         reference = reference_c(2, pu(i), pt(i), beta, lambda)
         call check(error, dev(values(i), reference) < id_tol, "PSE2 closure reference")
         if (allocated(error)) return
      end do
      do i = 1, nph
         reference = closure_pse2_phi(ph(i), pc(i))
         call check(error, abs(reference) < huge(1.0_wp), "PSE2 Phi finite")
         if (allocated(error)) return
      end do
   end subroutine test_pse2_consistency

   !> PSE3: dc/dt and d2c/dt2 against a 4-point central
   !> difference of c and of dc, at physical input values
   !>
   !> @param[out] error Test failure
   subroutine test_pse3_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp) :: hs, fd, ana, fp2, fp1, fm1, fm2, values(1)

      do i = 1, npt
         hs = fd_h*max(1.0_wp, abs(pt(i)))
         call closure_pse3_c([pu(i)], [pt(i) + 2.0_wp*hs], beta, lambda, values)
         fp2 = values(1)
         call closure_pse3_c([pu(i)], [pt(i) + hs], beta, lambda, values)
         fp1 = values(1)
         call closure_pse3_c([pu(i)], [pt(i) - hs], beta, lambda, values)
         fm1 = values(1)
         call closure_pse3_c([pu(i)], [pt(i) - 2.0_wp*hs], beta, lambda, values)
         fm2 = values(1)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_pse3_dc(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "PSE3 dc/dt vs FD")
         if (allocated(error)) return
         fp2 = closure_pse3_dc(pu(i), pt(i) + 2.0_wp*hs, beta, lambda)
         fp1 = closure_pse3_dc(pu(i), pt(i) + hs, beta, lambda)
         fm1 = closure_pse3_dc(pu(i), pt(i) - hs, beta, lambda)
         fm2 = closure_pse3_dc(pu(i), pt(i) - 2.0_wp*hs, beta, lambda)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_pse3_d2c(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "PSE3 d2c/dt2 vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_pse3_derivatives

   !> PSE3: sweeps agree with the closure formula and Phi is finite
   !>
   !> @param[out] error Test failure
   subroutine test_pse3_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: values(npt), reference
      integer :: i

      call closure_pse3_c(pu, pt, beta, lambda, values)
      do i = 1, npt
         reference = reference_c(3, pu(i), pt(i), beta, lambda)
         call check(error, dev(values(i), reference) < id_tol, "PSE3 closure reference")
         if (allocated(error)) return
      end do
      do i = 1, nph
         reference = closure_pse3_phi(ph(i), pc(i))
         call check(error, abs(reference) < huge(1.0_wp), "PSE3 Phi finite")
         if (allocated(error)) return
      end do
   end subroutine test_pse3_consistency

   !> PSE4: dc/dt and d2c/dt2 against a 4-point central
   !> difference of c and of dc, at physical input values
   !>
   !> @param[out] error Test failure
   subroutine test_pse4_derivatives(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      real(wp) :: hs, fd, ana, fp2, fp1, fm1, fm2, values(1)

      do i = 1, npt
         hs = fd_h*max(1.0_wp, abs(pt(i)))
         call closure_pse4_c([pu(i)], [pt(i) + 2.0_wp*hs], beta, lambda, values)
         fp2 = values(1)
         call closure_pse4_c([pu(i)], [pt(i) + hs], beta, lambda, values)
         fp1 = values(1)
         call closure_pse4_c([pu(i)], [pt(i) - hs], beta, lambda, values)
         fm1 = values(1)
         call closure_pse4_c([pu(i)], [pt(i) - 2.0_wp*hs], beta, lambda, values)
         fm2 = values(1)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_pse4_dc(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "PSE4 dc/dt vs FD")
         if (allocated(error)) return
         fp2 = closure_pse4_dc(pu(i), pt(i) + 2.0_wp*hs, beta, lambda)
         fp1 = closure_pse4_dc(pu(i), pt(i) + hs, beta, lambda)
         fm1 = closure_pse4_dc(pu(i), pt(i) - hs, beta, lambda)
         fm2 = closure_pse4_dc(pu(i), pt(i) - 2.0_wp*hs, beta, lambda)
         fd = (-fp2 + 8.0_wp*fp1 - 8.0_wp*fm1 + fm2)/(12.0_wp*hs)
         ana = closure_pse4_d2c(pu(i), pt(i), beta, lambda)
         call check(error, dev(fd, ana) < fd_tol, "PSE4 d2c/dt2 vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_pse4_derivatives

   !> PSE4: sweeps agree with the closure formula and Phi is finite
   !>
   !> @param[out] error Test failure
   subroutine test_pse4_consistency(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: values(npt), reference
      integer :: i

      call closure_pse4_c(pu, pt, beta, lambda, values)
      do i = 1, npt
         reference = reference_c(4, pu(i), pt(i), beta, lambda)
         call check(error, dev(values(i), reference) < id_tol, "PSE4 closure reference")
         if (allocated(error)) return
      end do
      do i = 1, nph
         reference = closure_pse4_phi(ph(i), pc(i))
         call check(error, abs(reference) < huge(1.0_wp), "PSE4 Phi finite")
         if (allocated(error)) return
      end do
   end subroutine test_pse4_consistency

   !> fast_exp, seen through HNC d2c = exp(t) at u = 0, stays within
   !> 1 ulp of the libm exponential up to ln(huge), and stays finite in
   !> the HNC kernels between 709.44, where n = 1024, and ln(huge)
   !>
   !> @param[out] error Test failure
   subroutine test_fast_exp(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 200001
      real(wp) :: x, e, x_max
      real(wp) :: t(2), u(2), c(2)
      integer :: i

      x_max = log(huge(1.0_wp))
      do i = 1, n
         x = -708.0_wp + (x_max + 708.0_wp)*real(i - 1, wp)/real(n - 1, wp)
         x = min(x, x_max)
         e = closure_hnc_d2c(0.0_wp, x, beta, lambda)
         call check(error, abs(e - exp(x)) <= spacing(exp(x)), "fast_exp within 1 ulp")
         if (allocated(error)) return
      end do

      e = closure_hnc_dc(0.0_wp, 709.5_wp, 1.0_wp, 1.0_wp)
      call check(error, abs(e) <= huge(1.0_wp), "HNC dc finite at t = 709.5")
      if (allocated(error)) return
      call check(error, abs(e - (exp(709.5_wp) - 1.0_wp)) <= spacing(exp(709.5_wp)), &
         & "HNC dc within 1 ulp at t = 709.5")
      if (allocated(error)) return

      u = 0.0_wp
      t = [709.5_wp, x_max]
      call closure_hnc_c(u, t, beta, lambda, c)
      call check(error, all(abs(c) <= huge(1.0_wp)), "HNC c finite up to ln(huge)")
      if (allocated(error)) return
      call check(error, all(abs(c - (exp(t) - 1.0_wp - t)) <= spacing(exp(t))), &
         & "HNC c within 1 ulp up to ln(huge)")
      if (allocated(error)) return
   end subroutine test_fast_exp

   !> Mixed absolute/relative deviation: `b` spans many orders of
   !> magnitude across the sampled domain, and may pass through
   !> zero, so neither a pure relative nor a pure absolute
   !> measure is usable on its own
   !>
   !> @param[in] a Computed value
   !> @param[in] b Reference value
   !> @return Scaled deviation
   pure function dev(a, b) result(d)
      !> Computed and reference values
      real(wp), intent(in) :: a, b
      !> Scaled deviation
      real(wp) :: d

      d = abs(a - b)/(1.0_wp + abs(b))
   end function dev


   !> Independent formulas cover both branches, the join and a strongly repulsive potential
   !>
   !> @param[out] error Test failure
   subroutine test_sweeps(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: u(9) = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 0.3_wp, -0.2_wp, 900.0_wp]
      real(wp), parameter :: t(9) = [-1.0_wp, -1.0e-10_wp, 0.0_wp, 1.0e-10_wp, &
         & 1.0_wp, 4.0_wp, 0.4_wp, -0.3_wp, 0.2_wp]
      real(wp) :: values(9)
      integer :: order, i

      do order = 0, 4
         call evaluate_c(order, u, t, beta, lambda, values)
         do i = 1, size(t)
            call check(error, dev(values(i), reference_c(order, u(i), t(i), beta, lambda)) < id_tol, &
               & "Closure sweep agrees with its mathematical reference")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_sweeps

   !> An unknown closure kind gives NaN instead of silently evaluating PSE4
   !>
   !> @param[out] error Test failure
   subroutine test_sweep_invalid_kind(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: u(3) = [0.0_wp, 0.3_wp, -0.2_wp], t(3) = [-0.5_wp, 0.0_wp, 0.7_wp]
      real(wp) :: c(3), dc(3), d2c(3)
      integer, parameter :: bad_kinds(3) = [-1, 5, 7]
      integer :: i

      do i = 1, size(bad_kinds)
         call closure_dc_sweep(u, t, beta, lambda, bad_kinds(i), dc)
         call check(error, all(ieee_is_nan(dc)), "dc sweep of an unknown kind is NaN")
         if (allocated(error)) return
         call closure_d2c_sweep(u, t, beta, lambda, bad_kinds(i), d2c)
         call check(error, all(ieee_is_nan(d2c)), "d2c sweep of an unknown kind is NaN")
         if (allocated(error)) return
         call closure_c_dc_sweep(u, t, beta, lambda, bad_kinds(i), c, dc)
         call check(error, all(ieee_is_nan(c)) .and. all(ieee_is_nan(dc)), "fused sweep of an unknown kind is NaN")
         if (allocated(error)) return
      end do
      call closure_dc_sweep(u, t, beta, lambda, 4, dc)
      call check(error, .not. any(ieee_is_nan(dc)), "PSE4 sweep stays finite")
   end subroutine test_sweep_invalid_kind

   !> Small tails, reversed sections and empty arrays use the same block interface
   !>
   !> @param[out] error Test failure
   subroutine test_sweep_sections(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u(17), t(17), values(17), strided(34)
      integer :: order, i

      do i = 1, size(t)
         u(i) = 0.1_wp*real(i - 9, wp)
         t(i) = 0.2_wp*real(i - 8, wp)
      end do
      do order = 0, 4
         call evaluate_c(order, u, t, beta, lambda, values)
         strided = 42.0_wp
         call evaluate_c(order, u(17:1:-1), t(17:1:-1), beta, lambda, strided(1::2))
         call check(error, all(strided(1::2) == values(17:1:-1)), "Reversed sweep sections")
         if (allocated(error)) return
         call check(error, all(strided(2::2) == 42.0_wp), "Sweep preserves excluded entries")
         if (allocated(error)) return
         call evaluate_c(order, u(:0), t(:0), beta, lambda, strided(:0))
         call check(error, all(strided(1::2) == values(17:1:-1)) .and. &
            & all(strided(2::2) == 42.0_wp), "Empty sweep preserves storage")
         if (allocated(error)) return
      end do
   end subroutine test_sweep_sections

   !> Dispatch test inputs to the five production sweeps
   subroutine evaluate_c(order, u, t, beta, lambda, c)
      integer, intent(in) :: order
      real(wp), intent(in) :: u(:), t(:), beta, lambda
      real(wp), intent(out) :: c(:)
      select case (order)
      case (0)
         call closure_hnc_c(u, t, beta, lambda, c)
      case (1)
         call closure_kh_c(u, t, beta, lambda, c)
      case (2)
         call closure_pse2_c(u, t, beta, lambda, c)
      case (3)
         call closure_pse3_c(u, t, beta, lambda, c)
      case (4)
         call closure_pse4_c(u, t, beta, lambda, c)
      case default
         error stop "Unknown test closure order"
      end select
   end subroutine evaluate_c

   !> Reference formula uses intrinsic exp and a factorial series, independent of kernel helpers
   pure function reference_c(order, u, t, beta, lambda) result(c)
      integer, intent(in) :: order
      real(wp), intent(in) :: u, t, beta, lambda
      real(wp) :: c, x, series, term
      integer :: k

      x = -beta*lambda*u + t
      if (order == 0 .or. x <= 0.0_wp) then
         c = exp(x) - 1.0_wp - t
      else
         series = 0.0_wp
         term = 1.0_wp
         do k = 1, order
            term = term*x/real(k, wp)
            series = series + term
         end do
         c = series - t
      end if
   end function reference_c

end module test_moz_closure_kernel
