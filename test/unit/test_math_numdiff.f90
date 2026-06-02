!> Faithful ports of the upstream tests
!> (tests/test1.f90, test2.f90, dsm_test.f90), adapted to the moist test-drive
!> harness. The upstream programs only print results; this port keeps their
!> exact problems but adds assertions against the known analytic answer.
!>   * test1: 6-function, 10-variable problem; finite-difference Jacobian for
!>            every available method is compared to the analytic Jacobian.
!>   * test2: scalar f(x) = x + sin(x); derivative compared to 1 + cos(x).
!>   * dsm_test: MINPACK DSM/FDJS sparse-Jacobian partitioning on the
!>            Coleman-Garbow-More neutron-kinetics sparsity pattern.
module test_math_numdiff
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use numerical_differentiation_module, only: numdiff_type, get_finite_diff_formula
   use dsm_module, only: dsm, fdjs
   implicit none (type, external)
   private

   public :: collect_math_numdiff

   !> Tolerance on finite-difference vs analytic derivatives.
   real(wp), parameter :: thr = 5.0e-5_wp

contains

   !> Collect all NumDiff tests
   subroutine collect_math_numdiff(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("jacobian-all-methods", test_methods_dense), &
         new_unittest("scalar-derivative", test_scalar_accuracy), &
         new_unittest("dsm-fdjs-partitioning", test_dsm_partitioning) &
         ]
   end subroutine collect_math_numdiff

   !> test1: compute the dense Jacobian of the 6x10 problem with every available
   !> finite-difference method and compare each to the analytic Jacobian.
   subroutine test_methods_dense(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 10, m = 6
      real(wp), parameter :: x(n) = 1.0_wp
      real(wp), parameter :: xlow(n) = -10.0_wp, xhigh(n) = 10.0_wp
      real(wp), parameter :: dpert(n) = 1.0e-5_wp
      type(numdiff_type) :: prob
      real(wp), allocatable :: jac(:, :)
      real(wp) :: jan(m, n)
      character(len=:), allocatable :: formula
      character(len=64) :: msg
      integer :: ids(48), k, i, nfound

      call analytic_jac(x, jan)

      ids(1:44) = [(i, i=1, 44)]
      ids(45:48) = [500, 600, 700, 800]
      nfound = 0

      do k = 1, size(ids)
         i = ids(k)
         call get_finite_diff_formula(i, formula)
         if (formula == '') cycle

         call prob%destroy()
         call prob%initialize(n, m, xlow, xhigh, 1, dpert, problem_func=prob6_func, &
                              sparsity_mode=1, jacobian_method=i, &
                              partition_sparsity_pattern=.false., cache_size=0)
         if (prob%failed()) cycle

         call prob%compute_jacobian_dense(x, jac)
         if (prob%failed()) cycle

         write (msg, '(a,i0)') "jacobian mismatch, method ", i
         call check(error, maxval(abs(jac - jan)) < thr, trim(msg))
         if (allocated(error)) return
         nfound = nfound + 1
      end do

      call check(error, nfound >= 1, "no valid finite-difference methods were exercised")
   end subroutine test_methods_dense

   !> test2: derivative of x + sin(x) at x = 1 via forward and central
   !> differences, compared to the analytic value 1 + cos(1).
   subroutine test_scalar_accuracy(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 1, m = 1
      real(wp), parameter :: x(n) = 1.0_wp
      real(wp), parameter :: xlow(n) = -1.0e5_wp, xhigh(n) = 1.0e5_wp
      real(wp), parameter :: dpert(n) = 1.0e-5_wp
      real(wp) :: deriv_true
      type(numdiff_type) :: prob
      real(wp), allocatable :: jac(:, :)
      integer :: methods(2), k, i
      character(len=64) :: msg

      deriv_true = 1.0_wp + cos(1.0_wp)
      methods = [1, 3]   ! forward, central

      do k = 1, size(methods)
         i = methods(k)
         call prob%destroy()
         call prob%initialize(n, m, xlow, xhigh, 1, dpert, problem_func=scalar_func, &
                              sparsity_mode=1, jacobian_method=i, &
                              partition_sparsity_pattern=.false., cache_size=0)
         write (msg, '(a,i0)') "scalar init failed, method ", i
         call check(error, .not. prob%failed(), trim(msg))
         if (allocated(error)) return

         call prob%compute_jacobian_dense(x, jac)
         write (msg, '(a,i0)') "scalar derivative mismatch, method ", i
         call check(error, jac(1, 1), deriv_true, thr=thr, message=trim(msg))
         if (allocated(error)) return
      end do

      call prob%destroy()
   end subroutine test_scalar_accuracy

   !> dsm_test: partition a sparse Jacobian with DSM and reconstruct it with
   !> FDJS on the neutron-kinetics pattern, checking the relative error for a
   !> range of problem sizes (alternating column- and row-oriented patterns).
   subroutine test_dsm_partitioning(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: indcol(6000), indrow(6000), ipntr(1201), jpntr(1201), ngrp(1200)
      real(wp) :: d(1200), fjac(6000), fjacd(1200), fvec(1200), x(1200), xd(1200)
      integer :: i, info, ip, j, jp, l, m, maxgrp, mingrp, n, nnz, numgrp
      logical :: col
      real(wp) :: errij, errmax, fjact, fjactr, s
      character(len=64) :: msg

      col = .true.

      do n = 300, 1200, 300
         m = n
         l = n/3

         ! Build the neutron-kinetics sparsity pattern.
         nnz = 0
         do j = 1, n
            nnz = nnz + 1
            indrow(nnz) = j
            indcol(nnz) = j
            if (mod(j, l) /= 0) then
               nnz = nnz + 1
               indrow(nnz) = j + 1
               indcol(nnz) = j
            end if
            if (j <= 2*l) then
               nnz = nnz + 1
               indrow(nnz) = j + l
               indcol(nnz) = j
               if (mod(j, l) /= 1) then
                  nnz = nnz + 1
                  indrow(nnz) = j - 1
                  indcol(nnz) = j
               end if
            end if
            nnz = nnz + 1
            if (j > l) then
               indrow(nnz) = j - l
            else
               indrow(nnz) = j + 2*l
            end if
            indcol(nnz) = j
         end do

         call dsm(m, n, nnz, indrow, indcol, ngrp, maxgrp, mingrp, info, ipntr, jpntr)
         write (msg, '(a,i0)') "dsm input error, n=", n
         call check(error, info > 0, trim(msg))
         if (allocated(error)) return

         do j = 1, n
            x(j) = real(j, wp)/real(n, wp)
         end do
         call dsm_fcn(n, x, indcol, ipntr, fvec)

         ! Approximate the Jacobian one column-group at a time.
         do numgrp = 1, maxgrp
            do j = 1, n
               d(j) = 0.0_wp
               if (ngrp(j) == numgrp) d(j) = 1.0e-6_wp
               xd(j) = x(j) + d(j)
            end do
            call dsm_fcn(n, xd, indcol, ipntr, fjacd)
            do i = 1, m
               fjacd(i) = fjacd(i) - fvec(i)
            end do
            if (col) then
               call fdjs(m, n, col, indrow, jpntr, ngrp, numgrp, d, fjacd, fjac)
            else
               call fdjs(m, n, col, indcol, ipntr, ngrp, numgrp, d, fjacd, fjac)
            end if
         end do

         ! Compare against the analytic Jacobian of dsm_fcn.
         errmax = 0.0_wp
         if (col) then
            do j = 1, n
               do jp = jpntr(j), jpntr(j + 1) - 1
                  i = indrow(jp)
                  s = 0.0_wp
                  do ip = ipntr(i), ipntr(i + 1) - 1
                     s = s + x(indcol(ip))
                  end do
                  s = s + x(i)
                  fjact = 1.0_wp + 2.0_wp*s
                  if (i == j) fjact = 2.0_wp*fjact
                  errij = fjac(jp) - fjact
                  if (fjact /= 0.0_wp) errij = errij/fjact
                  errmax = max(errmax, abs(errij))
               end do
            end do
         else
            do i = 1, m
               s = 0.0_wp
               do ip = ipntr(i), ipntr(i + 1) - 1
                  s = s + x(indcol(ip))
               end do
               s = s + x(i)
               fjactr = 1.0_wp + 2.0_wp*s
               do ip = ipntr(i), ipntr(i + 1) - 1
                  j = indcol(ip)
                  fjact = fjactr
                  if (i == j) fjact = 2.0_wp*fjact
                  errij = fjac(ip) - fjact
                  if (fjact /= 0.0_wp) errij = errij/fjact
                  errmax = max(errmax, abs(errij))
               end do
            end do
         end if

         write (msg, '(a,i0)') "dsm/fdjs reconstruction error too large, n=", n
         call check(error, errmax < 1.0e-4_wp, trim(msg))
         if (allocated(error)) return

         col = .not. col
      end do
   end subroutine test_dsm_partitioning

   !===========================================================================
   ! Problem definitions
   !===========================================================================

   !> The 6-function, 10-variable problem from upstream test1.
   subroutine prob6_func(me, x, f, funcs_to_compute)
      !> NumDiff instance (unused)
      class(numdiff_type), intent(inout) :: me
      !> Variables
      real(wp), dimension(:), intent(in) :: x
      !> Function values
      real(wp), dimension(:), intent(out) :: f
      !> Indices of the functions to evaluate
      integer, dimension(:), intent(in) :: funcs_to_compute

      if (any(funcs_to_compute == 1)) f(1) = x(1)*x(2) - x(3)**3
      if (any(funcs_to_compute == 2)) f(2) = x(3) - 1.0_wp
      if (any(funcs_to_compute == 3)) f(3) = x(4)*x(5)
      if (any(funcs_to_compute == 4)) f(4) = 2.0_wp*x(6) + sin(x(7))
      if (any(funcs_to_compute == 5)) f(5) = cos(x(8)) + sqrt(abs(x(9)))
      if (any(funcs_to_compute == 6)) f(6) = 1.0_wp/(1.0_wp + exp(x(10)))
   end subroutine prob6_func

   !> Analytic Jacobian of prob6_func.
   subroutine analytic_jac(x, jan)
      !> Variables
      real(wp), dimension(:), intent(in) :: x
      !> Analytic Jacobian (6,10)
      real(wp), dimension(:, :), intent(out) :: jan

      jan = 0.0_wp
      jan(1, 1) = x(2)
      jan(1, 2) = x(1)
      jan(1, 3) = -3.0_wp*x(3)**2
      jan(2, 3) = 1.0_wp
      jan(3, 4) = x(5)
      jan(3, 5) = x(4)
      jan(4, 6) = 2.0_wp
      jan(4, 7) = cos(x(7))
      jan(5, 8) = -sin(x(8))
      jan(5, 9) = sign(0.5_wp/sqrt(abs(x(9))), x(9))
      jan(6, 10) = -exp(x(10))/(1.0_wp + exp(x(10)))**2
   end subroutine analytic_jac

   !> Scalar problem f(x) = x + sin(x) from upstream test2.
   subroutine scalar_func(me, x, f, funcs_to_compute)
      !> NumDiff instance (unused)
      class(numdiff_type), intent(inout) :: me
      !> Variables
      real(wp), dimension(:), intent(in) :: x
      !> Function values
      real(wp), dimension(:), intent(out) :: f
      !> Indices of the functions to evaluate
      integer, dimension(:), intent(in) :: funcs_to_compute

      if (any(funcs_to_compute == 1)) f(1) = x(1) + sin(x(1))
   end subroutine scalar_func

   !> Quadratic test function for the DSM/FDJS check (from upstream dsm_test).
   subroutine dsm_fcn(n, x, indcol, ipntr, fvec)
      !> Number of variables/functions
      integer, intent(in) :: n
      !> Variables
      real(wp), intent(in) :: x(n)
      !> Column indices of the row-oriented sparsity pattern
      integer, intent(in) :: indcol(*)
      !> Row pointers into indcol
      integer, intent(in) :: ipntr(n + 1)
      !> Function values
      real(wp), intent(out) :: fvec(n)

      integer :: i, ip
      real(wp) :: s

      do i = 1, n
         s = 0.0_wp
         do ip = ipntr(i), ipntr(i + 1) - 1
            s = s + x(indcol(ip))
         end do
         s = s + x(i)
         fvec(i) = s*(1.0_wp + s) + 1.0_wp
      end do
   end subroutine dsm_fcn

end module test_math_numdiff
