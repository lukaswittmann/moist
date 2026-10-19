!> Test suite for the Wendland smoothing kernels in moist_math_smoothing_kernels
!>
!> `init` is the only fallible operation: it selects a normalization and binds
!> the evaluation procedures for one (order, dimension) pair, and there is no
!> kernel to bind for an unsupported combination or a non-positive smoothing
!> length. Those used to terminate the process; they now report, which also
!> makes the guarantee testable that a rejected `init` leaves the kernel
!> detached rather than half-configured
module test_math_smoothing_kernels
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_smoothing_kernels, only: smoothing_kernel_wendland_type
   use test_helpers, only: fd4_scalar, fd4_offsets
   implicit none(type, external)
   private

   public :: collect_math_smoothing_kernels

   !> Smoothing length used by every well-formed init below
   real(wp), parameter :: h_ref = 0.5_wp

contains

   !> Collect all smoothing-kernel tests
   subroutine collect_math_smoothing_kernels(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("supported_combinations", test_supported_combinations), &
                  new_unittest("numeric_channels", test_numeric_channels), &
                  new_unittest("normalization_quadrature", test_normalization_quadrature), &
                  new_unittest("unsupported_dimension", test_unsupported_dimension), &
                  new_unittest("unsupported_order", test_unsupported_order), &
                  new_unittest("nonpositive_h", test_nonpositive_h), &
                  new_unittest("failed_init_detaches", test_failed_init_detaches) &
                  ]

   end subroutine collect_math_smoothing_kernels

   !> Every documented (order, dimension) pair initializes and evaluates
   subroutine test_supported_combinations(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel
      type(moist_error_type), allocatable :: rejected
      integer :: orders(3), iorder, idim

      orders = [2, 4, 6]

      do iorder = 1, size(orders)
         do idim = 1, 3
            call kernel%init(order=orders(iorder), dimension=idim, h=h_ref, error=rejected)
            call check(error,.not. allocated(rejected), "supported combination accepted")
            if (allocated(error)) return

            call check(error, associated(kernel%compute), "kernel bound its evaluator")
            if (allocated(error)) return
            call check(error, associated(kernel%compute_deriv), "kernel bound its derivative")
            if (allocated(error)) return

            !> A Wendland kernel is positive at the origin and vanishes at its
            !> support radius of 2h, which is enough to tell a bound evaluator
            !> from a stale one
            call check(error, kernel%f0(0.0_wp) > 0.0_wp, "kernel is positive at r = 0")
            if (allocated(error)) return
            call check(error, kernel%f0(2.0_wp*h_ref), 0.0_wp, "kernel vanishes at r = 2h")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_supported_combinations

   !> Interior values, radial and smoothing-length derivatives, and vector channels
   !>
   !> Reference polynomials are Wendland's psi_{d,k}(s) at s = q/2 (support 2h),
   !> scaled to W(0) = 1: psi_{1,k} in 1D, psi_{3,k} in 2D/3D, k = order/2
   !> (H. Wendland, Adv. Comput. Math. 4, 389 (1995)); `normalization_quadrature`
   !> checks the alpha table independently
   !>
   !> Derivatives are compared with `fd4_scalar` at step 1e-3 h, relative tolerance
   !> 1e-7: measured worst relative deviation 3.2e-11 (radial) and 1.4e-10 (h), macOS
   !> arm64 gfortran 14.3, one thread; both pre-fix C6 dW/dq formulas fail it
   subroutine test_numeric_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel, shifted
      type(moist_error_type), allocatable :: rejected
      integer :: orders(3), io, dim, ih, iq, is, k
      real(wp) :: h, q, s, r, a, poly, expected, fd_r, fd_h, step, f(4), grad(3), x(3)
      real(wp), parameter :: lengths(2) = [0.5_wp, 1.3_wp]
      real(wp), parameter :: alpha(3, 3) = reshape([ &
                                                   5.0_wp/8.0_wp, 0.75_wp, 55.0_wp/64.0_wp, &
                                                   7.0_wp/(4.0_wp*pi), 9.0_wp/(4.0_wp*pi), 39.0_wp/(14.0_wp*pi), &
                                                   21.0_wp/(16.0_wp*pi), 495.0_wp/(256.0_wp*pi), 1365.0_wp/(512.0_wp*pi)], [3, 3])
      !> Finite-difference step, in units of h
      real(wp), parameter :: fd_step = 1.0e-3_wp
      !> Relative tolerance of the derivative checks, scaled by the FD value
      real(wp), parameter :: fd_rel = 1.0e-7_wp

      orders = [2, 4, 6]
      do io = 1, 3
         do dim = 1, 3
            do iq = 1, 3
               q = 0.2_wp + 0.6_wp*real(iq - 1, wp)
               s = 0.5_wp*q
               do ih = 1, 2
                  h = lengths(ih)
                  r = q*h
                  call kernel%init(orders(io), dim, h, rejected)
                  call check(error,.not. allocated(rejected), "numeric init accepted")
                  if (allocated(error)) return
                  call check(error, kernel%order == orders(io), "stored order")
                  if (allocated(error)) return
                  call check(error, kernel%dimension == dim, "stored dimension")
                  if (allocated(error)) return
                  call check(error, kernel%h, h, "stored smoothing length")
                  if (allocated(error)) return
                  a = alpha(io, dim)/h**dim
                  call check(error, kernel%f0(0.0_wp), a, "normalization", thr=1.0e-12_wp)
                  if (allocated(error)) return
                  select case (orders(io))
                  case (2)
                     if (dim == 1) then
                        poly = (1.0_wp - s)**3*(3.0_wp*s + 1.0_wp)
                     else
                        poly = (1.0_wp - s)**4*(4.0_wp*s + 1.0_wp)
                     end if
                  case (4)
                     if (dim == 1) then
                        poly = (1.0_wp - s)**5*(8.0_wp*s*s + 5.0_wp*s + 1.0_wp)
                     else
                        poly = (1.0_wp - s)**6*(35.0_wp*s*s + 18.0_wp*s + 3.0_wp)/3.0_wp
                     end if
                  case (6)
                     if (dim == 1) then
                        poly = (1.0_wp - s)**7*(21.0_wp*s**3 + 19.0_wp*s*s + 7.0_wp*s + 1.0_wp)
                     else
                        poly = (1.0_wp - s)**8*(32.0_wp*s**3 + 25.0_wp*s*s + 8.0_wp*s + 1.0_wp)
                     end if
                  case default
                     call check(error, .false., "unknown test order")
                     return
                  end select
                  call check(error, kernel%f0(r), a*poly, "interior kernel", thr=1.0e-12_wp)
                  if (allocated(error)) return
                  call check(error, kernel%f0(2.5_wp*h), 0.0_wp, "outside support value")
                  if (allocated(error)) return
                  call check(error, kernel%f1(2.5_wp*h), 0.0_wp, "outside support derivative")
                  if (allocated(error)) return
                  call check(error, kernel%f1(0.0_wp), 0.0_wp, "origin derivative")
                  if (allocated(error)) return
                  call check(error, kernel%f1(r) < 0.0_wp, "interior derivative sign")
                  if (allocated(error)) return

                  step = fd_step*h
                  do is = 1, 4
                     f(is) = kernel%f0(r + fd4_offsets(is)*step)
                  end do
                  call fd4_scalar(f(1), f(2), f(3), f(4), step, fd_r, error)
                  if (allocated(error)) return
                  call check(error, kernel%f1(r), fd_r, "radial finite difference", thr=fd_rel*abs(fd_r))
                  if (allocated(error)) return
                  do is = 1, 4
                     call shifted%init(orders(io), dim, h + fd4_offsets(is)*step, rejected)
                     call check(error,.not. allocated(rejected), "shifted smoothing length accepted")
                     if (allocated(error)) return
                     f(is) = shifted%f0(r)
                  end do
                  call fd4_scalar(f(1), f(2), f(3), f(4), step, fd_h, error)
                  if (allocated(error)) return
                  call check(error, kernel%gradient_h(r), fd_h, "h finite difference", thr=fd_rel*abs(fd_h))
                  if (allocated(error)) return

                  expected = -real(dim, wp)*a/h
                  call check(error, kernel%gradient_h(0.0_wp), expected, "h derivative at origin", thr=1.0e-12_wp)
                  if (allocated(error)) return
                  x = r*[2.0_wp/3.0_wp, -1.0_wp/3.0_wp, 2.0_wp/3.0_wp]
                  call kernel%gradient(r, x, grad)
                  do k = 1, 3
                     call check(error, grad(k), fd_r*x(k)/r, "gradient components", thr=fd_rel*abs(fd_r))
                     if (allocated(error)) return
                  end do
                  call kernel%gradient(0.0_wp, x, grad)
                  call check(error, maxval(abs(grad)), 0.0_wp, "zero radius gradient")
                  if (allocated(error)) return
                  call check(error, kernel%gradient_h(2.5_wp*h), 0.0_wp, "outside support h derivative")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_numeric_channels

   !> Every kernel integrates to one over R^d, independently of the alpha table
   !>
   !> Integrates S_d r^(d-1) W(r) over the support [0, 2h], S_d = 2, 2 pi, 4 pi,
   !> with 7-point Gauss-Legendre (nodes and weights from mpmath to 20 digits, as
   !> in Abramowitz & Stegun Table 25.4), exact for the integrands here (degree
   !> <= 13); measured worst deviation 4.4e-16, macOS arm64 gfortran 14.3, one thread
   subroutine test_normalization_quadrature(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel
      type(moist_error_type), allocatable :: rejected
      integer :: orders(3), io, dim, ih, i
      real(wp) :: h, r, total
      real(wp), parameter :: lengths(2) = [0.5_wp, 1.3_wp]
      !> Surface of the unit sphere in 1, 2 and 3 dimensions
      real(wp), parameter :: sphere(3) = [2.0_wp, 2.0_wp*pi, 4.0_wp*pi]
      !> Gauss-Legendre nodes on [-1, 1]
      real(wp), parameter :: nodes(7) = [ &
                             -0.94910791234275852453_wp, -0.74153118559939443986_wp, &
                             -0.40584515137739716691_wp, 0.0_wp, 0.40584515137739716691_wp, &
                             0.74153118559939443986_wp, 0.94910791234275852453_wp]
      !> Gauss-Legendre weights on [-1, 1]
      real(wp), parameter :: weights(7) = [ &
                             0.12948496616886969327_wp, 0.27970539148927666790_wp, &
                             0.38183005050511894495_wp, 0.41795918367346938776_wp, 0.38183005050511894495_wp, &
                             0.27970539148927666790_wp, 0.12948496616886969327_wp]

      orders = [2, 4, 6]
      do io = 1, 3
         do dim = 1, 3
            do ih = 1, 2
               h = lengths(ih)
               call kernel%init(orders(io), dim, h, rejected)
               call check(error,.not. allocated(rejected), "quadrature init accepted")
               if (allocated(error)) return
               ! Map [-1, 1] onto the support [0, 2h]
               total = 0.0_wp
               do i = 1, 7
                  r = h*(1.0_wp + nodes(i))
                  total = total + h*weights(i)*sphere(dim)*r**(dim - 1)*kernel%f0(r)
               end do
               call check(error, total, 1.0_wp, "kernel integrates to one", thr=1.0e-12_wp)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_normalization_quadrature

   !> Each order rejects a dimension outside 1..3 and says which one
   subroutine test_unsupported_dimension(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel
      type(moist_error_type), allocatable :: rejected
      integer :: orders(3), iorder

      orders = [2, 4, 6]

      do iorder = 1, size(orders)
         call kernel%init(order=orders(iorder), dimension=4, h=h_ref, error=rejected)

         call check(error, allocated(rejected), "dimension 4 rejected")
         if (allocated(error)) return
         call check(error, index(rejected%message, "unsupported dimension 4") > 0, &
                    "message names the offending dimension")
         if (allocated(error)) return
         deallocate (rejected)
      end do

   end subroutine test_unsupported_dimension

   !> An order with no Wendland form is rejected and named
   subroutine test_unsupported_order(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel
      type(moist_error_type), allocatable :: rejected

      call kernel%init(order=3, dimension=2, h=h_ref, error=rejected)

      call check(error, allocated(rejected), "order 3 rejected")
      if (allocated(error)) return
      call check(error, index(rejected%message, "unsupported order 3") > 0, &
                 "message names the offending order")

   end subroutine test_unsupported_order

   !> A non-positive smoothing length is refused rather than divided by
   subroutine test_nonpositive_h(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel
      type(moist_error_type), allocatable :: rejected

      call kernel%init(order=2, dimension=2, h=0.0_wp, error=rejected)
      call check(error, allocated(rejected), "zero smoothing length rejected")
      if (allocated(error)) return
      call check(error, index(rejected%message, "must be positive") > 0, &
                 "message explains the requirement")
      if (allocated(error)) return
      deallocate (rejected)

      call kernel%init(order=2, dimension=2, h=-1.0_wp, error=rejected)
      call check(error, allocated(rejected), "negative smoothing length rejected")

   end subroutine test_nonpositive_h

   !> A rejected re-init must not leave the previous kernel in place
   !>
   !> Without this the caller could ignore the error and keep evaluating a
   !> kernel whose normalization belongs to the previous, unrelated request
   subroutine test_failed_init_detaches(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(smoothing_kernel_wendland_type) :: kernel
      type(moist_error_type), allocatable :: rejected

      call kernel%init(order=2, dimension=2, h=h_ref, error=rejected)
      call check(error,.not. allocated(rejected), "first init accepted")
      if (allocated(error)) return
      call check(error, associated(kernel%compute), "first init bound an evaluator")
      if (allocated(error)) return

      call kernel%init(order=5, dimension=2, h=h_ref, error=rejected)
      call check(error, allocated(rejected), "second init rejected")
      if (allocated(error)) return
      call check(error,.not. associated(kernel%compute), &
                 "rejected init detached the evaluator")
      if (allocated(error)) return
      call check(error,.not. associated(kernel%compute_deriv), &
                 "rejected init detached the derivative")

   end subroutine test_failed_init_detaches

end module test_math_smoothing_kernels
