!> Test suite for the Wendland smoothing kernels in moist_math_smoothing_kernels.
!>
!> `init` is the only fallible operation: it selects a normalization and binds
!> the evaluation procedures for one (order, dimension) pair, and there is no
!> kernel to bind for an unsupported combination or a non-positive smoothing
!> length. Those used to terminate the process; they now report, which also
!> makes the guarantee testable that a rejected `init` leaves the kernel
!> detached rather than half-configured.
module test_math_smoothing_kernels
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_smoothing_kernels, only: wendland_kernel_type
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
                  new_unittest("supported-combinations", test_supported_combinations), &
                  new_unittest("unsupported-dimension", test_unsupported_dimension), &
                  new_unittest("unsupported-order", test_unsupported_order), &
                  new_unittest("nonpositive-h", test_nonpositive_h), &
                  new_unittest("failed-init-detaches", test_failed_init_detaches) &
                  ]

   end subroutine collect_math_smoothing_kernels

   !> Every documented (order, dimension) pair initializes and evaluates.
   subroutine test_supported_combinations(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(wendland_kernel_type) :: kernel
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
            !> from a stale one.
            call check(error, kernel%f0(0.0_wp) > 0.0_wp, "kernel is positive at r = 0")
            if (allocated(error)) return
            call check(error, kernel%f0(2.0_wp*h_ref), 0.0_wp, "kernel vanishes at r = 2h")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_supported_combinations

   !> Each order rejects a dimension outside 1..3 and says which one.
   subroutine test_unsupported_dimension(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(wendland_kernel_type) :: kernel
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

   !> An order with no Wendland form is rejected and named.
   subroutine test_unsupported_order(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(wendland_kernel_type) :: kernel
      type(moist_error_type), allocatable :: rejected

      call kernel%init(order=3, dimension=2, h=h_ref, error=rejected)

      call check(error, allocated(rejected), "order 3 rejected")
      if (allocated(error)) return
      call check(error, index(rejected%message, "unsupported order 3") > 0, &
                 "message names the offending order")

   end subroutine test_unsupported_order

   !> A non-positive smoothing length is refused rather than divided by.
   subroutine test_nonpositive_h(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(wendland_kernel_type) :: kernel
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

   !> A rejected re-init must not leave the previous kernel in place.
   !>
   !> Without this the caller could ignore the error and keep evaluating a
   !> kernel whose normalization belongs to the previous, unrelated request.
   subroutine test_failed_init_detaches(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(wendland_kernel_type) :: kernel
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
