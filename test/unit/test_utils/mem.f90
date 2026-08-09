!> Test suite for the array growth helpers in moist_utils_mem
module test_utils_mem
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_mem, only: grow_array
   implicit none(type, external)
   private

   public :: collect_utils_mem

contains

   !> Collect all array-growth tests
   subroutine collect_utils_mem(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("grow-preserves-and-fills", test_grow_preserves_and_fills), &
                  new_unittest("grow-from-unallocated", test_grow_from_unallocated), &
                  new_unittest("same-size-is-noop", test_same_size_is_noop), &
                  new_unittest("shrink-real-1d-reports", test_shrink_real_1d_reports), &
                  new_unittest("shrink-int-1d-reports", test_shrink_int_1d_reports), &
                  new_unittest("shrink-logical-1d-reports", test_shrink_logical_1d_reports), &
                  new_unittest("shrink-real-2d-reports", test_shrink_real_2d_reports), &
                  new_unittest("dim1-change-reports", test_dim1_change_reports) &
                  ]

   end subroutine collect_utils_mem

   !> Growing keeps the old contents and fills the new tail
   subroutine test_grow_preserves_and_fills(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:)
      real(wp), allocatable :: m(:, :)
      integer, allocatable :: n(:)
      logical, allocatable :: l(:)
      type(moist_error_type), allocatable :: refused

      a = [1.0_wp, 2.0_wp, 3.0_wp]
      call grow_array(a, 5, fill_value=-1.0_wp, error=refused)
      call check(error,.not. allocated(refused), "growing a 1d real array succeeds")
      if (allocated(error)) return
      call check(error, size(a), 5, "1d real grew to the requested size")
      if (allocated(error)) return
      call check(error, a(1), 1.0_wp, "leading element preserved")
      if (allocated(error)) return
      call check(error, a(3), 3.0_wp, "last old element preserved")
      if (allocated(error)) return
      call check(error, a(5), -1.0_wp, "new tail took the fill value")
      if (allocated(error)) return

      allocate (m(3, 2), source=7.0_wp)
      call grow_array(m, 3, 4, fill_value=0.0_wp, error=refused)
      call check(error,.not. allocated(refused), "growing a 2d real array succeeds")
      if (allocated(error)) return
      call check(error, size(m, 2), 4, "2d real grew along the second extent")
      if (allocated(error)) return
      call check(error, size(m, 1), 3, "2d real kept its first extent")
      if (allocated(error)) return
      call check(error, m(2, 2), 7.0_wp, "old column preserved")
      if (allocated(error)) return
      call check(error, m(2, 4), 0.0_wp, "new column took the fill value")
      if (allocated(error)) return

      n = [4, 5]
      call grow_array(n, 3, fill_value=9, error=refused)
      call check(error,.not. allocated(refused), "growing a 1d integer array succeeds")
      if (allocated(error)) return
      call check(error, n(2), 5, "1d integer preserved")
      if (allocated(error)) return
      call check(error, n(3), 9, "1d integer filled")
      if (allocated(error)) return

      l = [.true.]
      call grow_array(l, 2, fill_value=.true., error=refused)
      call check(error,.not. allocated(refused), "growing a 1d logical array succeeds")
      if (allocated(error)) return
      call check(error, l(1), "1d logical preserved")
      if (allocated(error)) return
      call check(error, l(2), "1d logical filled")

   end subroutine test_grow_preserves_and_fills

   !> An unallocated array is a valid starting point, not an error
   subroutine test_grow_from_unallocated(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:)
      real(wp), allocatable :: m(:, :)
      type(moist_error_type), allocatable :: refused

      call grow_array(a, 4, fill_value=2.0_wp, error=refused)
      call check(error,.not. allocated(refused), "growing from unallocated succeeds")
      if (allocated(error)) return
      call check(error, size(a), 4, "grew from unallocated")
      if (allocated(error)) return
      call check(error, a(1), 2.0_wp, "fill applied throughout")
      if (allocated(error)) return

      !> The rank-2 case takes its first extent from the request when there is
      !> no existing array to match, so this must not trip the dim1 guard.
      call grow_array(m, 3, 2, fill_value=1.0_wp, error=refused)
      call check(error,.not. allocated(refused), "2d growth from unallocated succeeds")
      if (allocated(error)) return
      call check(error, size(m, 1), 3, "2d first extent taken from the request")
      if (allocated(error)) return
      call check(error, size(m, 2), 2, "2d second extent taken from the request")

   end subroutine test_grow_from_unallocated

   !> Asking for the size an array already has changes nothing
   subroutine test_same_size_is_noop(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:)
      type(moist_error_type), allocatable :: refused

      a = [1.0_wp, 2.0_wp]
      call grow_array(a, 2, fill_value=99.0_wp, error=refused)
      call check(error,.not. allocated(refused), "a no-op resize is not an error")
      if (allocated(error)) return
      call check(error, size(a), 2, "size unchanged")
      if (allocated(error)) return
      call check(error, a(2), 2.0_wp, "contents untouched, fill not applied")

   end subroutine test_same_size_is_noop

   !> A shrink request is reported and leaves the array intact
   subroutine test_shrink_real_1d_reports(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:)
      type(moist_error_type), allocatable :: refused

      a = [1.0_wp, 2.0_wp, 3.0_wp]
      call grow_array(a, 1, error=refused)

      call check(error, allocated(refused), "shrink reported instead of terminating")
      if (allocated(error)) return
      call check(error, index(refused%message, "Cannot shrink") > 0, &
                 "message names the refused operation")
      if (allocated(error)) return

      !> The caller unwinds on this error, so the array it was handed must still
      !> be the one it had: same size, same contents.
      call check(error, size(a), 3, "array kept its size")
      if (allocated(error)) return
      call check(error, a(3), 3.0_wp, "array kept its contents")

   end subroutine test_shrink_real_1d_reports

   !> The integer specific reports too, and names itself
   subroutine test_shrink_int_1d_reports(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: n(:)
      type(moist_error_type), allocatable :: refused

      n = [1, 2, 3, 4]
      call grow_array(n, 2, error=refused)

      call check(error, allocated(refused), "shrink reported")
      if (allocated(error)) return
      call check(error, index(refused%message, "grow_array_int_1d") > 0, &
                 "message identifies the integer specific")
      if (allocated(error)) return
      call check(error, size(n), 4, "array kept its size")

   end subroutine test_shrink_int_1d_reports

   !> The logical specific reports too
   subroutine test_shrink_logical_1d_reports(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      logical, allocatable :: l(:)
      type(moist_error_type), allocatable :: refused

      l = [.true., .false., .true.]
      call grow_array(l, 1, error=refused)

      call check(error, allocated(refused), "shrink reported")
      if (allocated(error)) return
      call check(error, index(refused%message, "grow_array_logical_1d") > 0, &
                 "message identifies the logical specific")
      if (allocated(error)) return
      call check(error, size(l), 3, "array kept its size")

   end subroutine test_shrink_logical_1d_reports

   !> The rank-2 specific reports a shrink of the second extent
   subroutine test_shrink_real_2d_reports(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: m(:, :)
      type(moist_error_type), allocatable :: refused

      allocate (m(3, 5), source=1.0_wp)
      call grow_array(m, 3, 2, error=refused)

      call check(error, allocated(refused), "shrink reported")
      if (allocated(error)) return
      call check(error, index(refused%message, "Cannot shrink") > 0, &
                 "message names the refused operation")
      if (allocated(error)) return
      call check(error, size(m, 2), 5, "array kept its second extent")

   end subroutine test_shrink_real_2d_reports

   !> Re-shaping the leading extent is refused even when the array would grow
   subroutine test_dim1_change_reports(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: m(:, :)
      type(moist_error_type), allocatable :: refused

      allocate (m(3, 2), source=1.0_wp)

      !> Second extent grows here, so only the first-extent guard can fire.
      call grow_array(m, 4, 6, error=refused)

      call check(error, allocated(refused), "first-extent change reported")
      if (allocated(error)) return
      call check(error, index(refused%message, "first dimension") > 0, &
                 "message names the first dimension")
      if (allocated(error)) return
      call check(error, size(m, 1), 3, "array kept its first extent")
      if (allocated(error)) return
      call check(error, size(m, 2), 2, "array kept its second extent")

   end subroutine test_dim1_change_reports

end module test_utils_mem
