!> Test suite for the array growth helpers in moist_utils_mem
module test_utils_mem
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_mem, only: grow_array, filter_array
   implicit none(type, external)
   private

   public :: collect_utils_mem

contains

   !> Collect all array-growth tests
   subroutine collect_utils_mem(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("grow_preserves_and_fills", test_grow_preserves_and_fills), &
                  new_unittest("grow_from_unallocated", test_grow_from_unallocated), &
                  new_unittest("same_size_is_noop", test_same_size_is_noop), &
                  new_unittest("shrink_real_1d_reports", test_shrink_real_1d_reports), &
                  new_unittest("shrink_int_1d_reports", test_shrink_int_1d_reports), &
                  new_unittest("shrink_logical_1d_reports", test_shrink_logical_1d_reports), &
                  new_unittest("shrink_real_2d_reports", test_shrink_real_2d_reports), &
                  new_unittest("dim1_change_reports", test_dim1_change_reports), &
                  new_unittest("grow_real_1d_contracts", test_grow_real_1d_contracts), &
                  new_unittest("filter_real_1d_contracts", test_filter_real_1d_contracts), &
                  new_unittest("grow_real_2d_contracts", test_grow_real_2d_contracts), &
                  new_unittest("filter_real_2d_contracts", test_filter_real_2d_contracts), &
                  new_unittest("grow_int_1d_contracts", test_grow_int_1d_contracts), &
                  new_unittest("filter_int_1d_contracts", test_filter_int_1d_contracts), &
                  new_unittest("grow_logical_1d_contracts", test_grow_logical_1d_contracts), &
                  new_unittest("filter_logical_1d_contracts", test_filter_logical_1d_contracts) &
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
      !> no existing array to match, so this must not trip the dim1 guard
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
      !> be the one it had: same size, same contents
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
      if (allocated(error)) return
      call check(error, all(n == [1, 2, 3, 4]), "refusal preserves every element")

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
      if (allocated(error)) return
      call check(error, all(l .eqv. [.true., .false., .true.]), "refusal preserves every element")

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
      if (allocated(error)) return
      call check(error, all(abs(m - 1.0_wp) <= 1.0e-12_wp), "refusal preserves every element")

   end subroutine test_shrink_real_2d_reports

   !> Re-shaping the leading extent is refused even when the array would grow
   subroutine test_dim1_change_reports(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: m(:, :)
      type(moist_error_type), allocatable :: refused

      allocate (m(3, 2), source=1.0_wp)

      !> Second extent grows here, so only the first-extent guard can fire
      call grow_array(m, 4, 6, error=refused)

      call check(error, allocated(refused), "first-extent change reported")
      if (allocated(error)) return
      call check(error, index(refused%message, "first dimension") > 0, &
                 "message names the first dimension")
      if (allocated(error)) return
      call check(error, size(m, 1), 3, "array kept its first extent")
      if (allocated(error)) return
      call check(error, size(m, 2), 2, "array kept its second extent")
      if (allocated(error)) return
      call check(error, all(abs(m - 1.0_wp) <= 1.0e-12_wp), "refusal preserves every element")

   end subroutine test_dim1_change_reports

   !> Verify every element, default fill, no-op, and empty growth
   subroutine test_grow_real_1d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:), expected(:)
      type(moist_error_type), allocatable :: refused

      a = [2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp]
      expected = a
      call grow_array(a, 6, fill_value=-9.0_wp, error=refused)
      call check(error,.not. allocated(refused), "explicit growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 6, "explicit growth reaches the requested size")
      if (allocated(error)) return
      call check(error, all(abs(a(1:4) - expected) <= 1.0e-12_wp), "all old values preserved")
      if (allocated(error)) return
      call check(error, all(abs(a(5:6) + 9.0_wp) <= 1.0e-12_wp), "all new values filled")
      if (allocated(error)) return
      expected = a
      call grow_array(a, 6, error=refused)
      call check(error,.not. allocated(refused), "same size succeeds")
      if (allocated(error)) return
      call check(error, size(a), 6, "same size keeps the size")
      if (allocated(error)) return
      call check(error, all(abs(a - expected) <= 1.0e-12_wp), "same size preserves all values")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 3, error=refused)
      call check(error,.not. allocated(refused), "unallocated default growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 3, "unallocated growth reaches the requested size")
      if (allocated(error)) return
      call check(error, all(abs(a - 0.0_wp) <= 1.0e-12_wp), "default fill covers every element")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 0, error=refused)
      call check(error,.not. allocated(refused), "zero request on an unallocated array is accepted")
      if (allocated(error)) return
      if (allocated(a)) deallocate (a)
      allocate (a(0))
      call grow_array(a, 2, error=refused)
      call check(error,.not. allocated(refused), "allocated empty array grows")
      if (allocated(error)) return
      call check(error, size(a), 2, "empty array grows to the requested size")
      if (allocated(error)) return
      call check(error, all(abs(a - 0.0_wp) <= 1.0e-12_wp), "empty growth fills every element")

   end subroutine test_grow_real_1d_contracts

   !> Verify selection, stable order, prefix extent, empty masks, and allocation
   subroutine test_filter_real_1d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:), expected(:)
      logical :: keep(4)

      keep = [.true., .false., .true., .true.]
      call filter_array(a, 3, keep, 2)
      call check(error,.not. allocated(a), "unallocated filter is skipped")
      if (allocated(error)) return
      a = [2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp]
      expected = a([1, 3])
      call filter_array(a, 3, keep, 2)
      call check(error, size(a, 1), 2, "filter uses nvalid extent")
      if (allocated(error)) return
      call check(error, all(abs(a - expected) <= 1.0e-12_wp), "filter preserves selected values and order")
      if (allocated(error)) return
      a = [2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp]
      expected = a
      keep = .true.
      call filter_array(a, 4, keep, 4)
      call check(error, size(a, 1), 4, "all-true mask keeps the extent")
      if (allocated(error)) return
      call check(error, all(abs(a - expected) <= 1.0e-12_wp), "all-true mask keeps every value")
      if (allocated(error)) return
      keep = .false.
      call filter_array(a, 4, keep, 0)
      call check(error, allocated(a), "all-false mask keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 1), 0, "all-false mask gives zero extent")
      if (allocated(error)) return
      call filter_array(a, 0, keep, 0)
      call check(error, allocated(a), "empty filter keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 1), 0, "empty prefix remains empty")

   end subroutine test_filter_real_1d_contracts

   !> Verify every element, default fill, no-op, and empty growth
   subroutine test_grow_real_2d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:, :), expected(:, :)
      type(moist_error_type), allocatable :: refused

      a = reshape([2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp], [2, 2])
      expected = a
      call grow_array(a, 2, 4, fill_value=-9.0_wp, error=refused)
      call check(error,.not. allocated(refused), "explicit growth succeeds")
      if (allocated(error)) return
      call check(error, size(a, 1), 2, "explicit growth keeps the first extent")
      if (allocated(error)) return
      call check(error, size(a, 2), 4, "explicit growth reaches the requested second extent")
      if (allocated(error)) return
      call check(error, all(abs(a(:, 1:2) - expected) <= 1.0e-12_wp), "all old values preserved")
      if (allocated(error)) return
      call check(error, all(abs(a(:, 3:4) + 9.0_wp) <= 1.0e-12_wp), "all new values filled")
      if (allocated(error)) return
      expected = a
      call grow_array(a, 2, 4, error=refused)
      call check(error,.not. allocated(refused), "same size succeeds")
      if (allocated(error)) return
      call check(error, size(a, 1), 2, "same size keeps the first extent")
      if (allocated(error)) return
      call check(error, size(a, 2), 4, "same size keeps the second extent")
      if (allocated(error)) return
      call check(error, all(abs(a - expected) <= 1.0e-12_wp), "same size preserves all values")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 2, 3, error=refused)
      call check(error,.not. allocated(refused), "unallocated default growth succeeds")
      if (allocated(error)) return
      call check(error, size(a, 1), 2, "unallocated growth takes the first extent")
      if (allocated(error)) return
      call check(error, size(a, 2), 3, "unallocated growth takes the second extent")
      if (allocated(error)) return
      call check(error, all(abs(a - 0.0_wp) <= 1.0e-12_wp), "default fill covers every element")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 2, 0, error=refused)
      call check(error,.not. allocated(refused), "zero request on an unallocated array is accepted")
      if (allocated(error)) return
      if (allocated(a)) deallocate (a)
      allocate (a(2, 0))
      call grow_array(a, 2, 2, error=refused)
      call check(error,.not. allocated(refused), "allocated empty array grows")
      if (allocated(error)) return
      call check(error, size(a, 1), 2, "empty array keeps the first extent")
      if (allocated(error)) return
      call check(error, size(a, 2), 2, "empty array grows to the requested second extent")
      if (allocated(error)) return
      call check(error, all(abs(a - 0.0_wp) <= 1.0e-12_wp), "empty growth fills every element")

   end subroutine test_grow_real_2d_contracts

   !> Verify selection, stable order, prefix extent, empty masks, and allocation
   subroutine test_filter_real_2d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:, :), expected(:, :)
      logical :: keep(4)

      keep = [.true., .false., .true., .true.]
      call filter_array(a, 3, keep, 2)
      call check(error,.not. allocated(a), "unallocated filter is skipped")
      if (allocated(error)) return
      a = reshape([2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp, 10.0_wp, 12.0_wp, 14.0_wp, 16.0_wp], [2, 4])
      expected = a(:, [1, 3])
      call filter_array(a, 3, keep, 2)
      call check(error, size(a, 2), 2, "filter uses nvalid extent")
      if (allocated(error)) return
      call check(error, size(a, 1), 2, "filter preserves row extent")
      if (allocated(error)) return
      call check(error, all(abs(a - expected) <= 1.0e-12_wp), "filter preserves selected values and order")
      if (allocated(error)) return
      a = reshape([2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp, 10.0_wp, 12.0_wp, 14.0_wp, 16.0_wp], [2, 4])
      expected = a
      keep = .true.
      call filter_array(a, 4, keep, 4)
      call check(error, size(a, 1), 2, "all-true mask keeps the row extent")
      if (allocated(error)) return
      call check(error, size(a, 2), 4, "all-true mask keeps the column extent")
      if (allocated(error)) return
      call check(error, all(abs(a - expected) <= 1.0e-12_wp), "all-true mask keeps every value")
      if (allocated(error)) return
      keep = .false.
      call filter_array(a, 4, keep, 0)
      call check(error, allocated(a), "all-false mask keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 2), 0, "all-false mask gives zero extent")
      if (allocated(error)) return
      call filter_array(a, 0, keep, 0)
      call check(error, allocated(a), "empty filter keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 2), 0, "empty prefix remains empty")

   end subroutine test_filter_real_2d_contracts

   !> Verify every element, default fill, no-op, and empty growth
   subroutine test_grow_int_1d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: a(:), expected(:)
      type(moist_error_type), allocatable :: refused

      a = [2, 4, 6, 8]
      expected = a
      call grow_array(a, 6, fill_value=-9, error=refused)
      call check(error,.not. allocated(refused), "explicit growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 6, "explicit growth reaches the requested size")
      if (allocated(error)) return
      call check(error, all(a(1:4) == expected), "all old values preserved")
      if (allocated(error)) return
      call check(error, all(a(5:6) == -9), "all new values filled")
      if (allocated(error)) return
      expected = a
      call grow_array(a, 6, error=refused)
      call check(error,.not. allocated(refused), "same size succeeds")
      if (allocated(error)) return
      call check(error, size(a), 6, "same size keeps the size")
      if (allocated(error)) return
      call check(error, all(a == expected), "same size preserves all values")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 3, error=refused)
      call check(error,.not. allocated(refused), "unallocated default growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 3, "unallocated growth reaches the requested size")
      if (allocated(error)) return
      call check(error, all(a == 0), "default fill covers every element")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 0, error=refused)
      call check(error,.not. allocated(refused), "zero request on an unallocated array is accepted")
      if (allocated(error)) return
      if (allocated(a)) deallocate (a)
      allocate (a(0))
      call grow_array(a, 2, error=refused)
      call check(error,.not. allocated(refused), "allocated empty array grows")
      if (allocated(error)) return
      call check(error, size(a), 2, "empty array grows to the requested size")
      if (allocated(error)) return
      call check(error, all(a == 0), "empty growth fills every element")

   end subroutine test_grow_int_1d_contracts

   !> Verify selection, stable order, prefix extent, empty masks, and allocation
   subroutine test_filter_int_1d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: a(:), expected(:)
      logical :: keep(4)

      keep = [.true., .false., .true., .true.]
      call filter_array(a, 3, keep, 2)
      call check(error,.not. allocated(a), "unallocated filter is skipped")
      if (allocated(error)) return
      a = [2, 4, 6, 8]
      expected = a([1, 3])
      call filter_array(a, 3, keep, 2)
      call check(error, size(a, 1), 2, "filter uses nvalid extent")
      if (allocated(error)) return
      call check(error, all(a == expected), "filter preserves selected values and order")
      if (allocated(error)) return
      a = [2, 4, 6, 8]
      expected = a
      keep = .true.
      call filter_array(a, 4, keep, 4)
      call check(error, size(a, 1), 4, "all-true mask keeps the extent")
      if (allocated(error)) return
      call check(error, all(a == expected), "all-true mask keeps every value")
      if (allocated(error)) return
      keep = .false.
      call filter_array(a, 4, keep, 0)
      call check(error, allocated(a), "all-false mask keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 1), 0, "all-false mask gives zero extent")
      if (allocated(error)) return
      call filter_array(a, 0, keep, 0)
      call check(error, allocated(a), "empty filter keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 1), 0, "empty prefix remains empty")

   end subroutine test_filter_int_1d_contracts

   !> Verify every element, default fill, no-op, and empty growth
   subroutine test_grow_logical_1d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      logical, allocatable :: a(:), expected(:)
      type(moist_error_type), allocatable :: refused

      a = [.true., .false., .true., .false.]
      expected = a
      call grow_array(a, 6, fill_value=.true., error=refused)
      call check(error,.not. allocated(refused), "explicit growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 6, "explicit growth reaches the requested size")
      if (allocated(error)) return
      call check(error, all(a(1:4) .eqv. expected), "all old values preserved")
      if (allocated(error)) return
      call check(error, all(a(5:6) .eqv. .true.), "all new values filled")
      if (allocated(error)) return
      expected = a
      call grow_array(a, 6, error=refused)
      call check(error,.not. allocated(refused), "same size succeeds")
      if (allocated(error)) return
      call check(error, size(a), 6, "same size keeps the size")
      if (allocated(error)) return
      call check(error, all(a .eqv. expected), "same size preserves all values")
      if (allocated(error)) return
      call grow_array(a, 8, fill_value=.false., error=refused)
      call check(error,.not. allocated(refused), "explicit false growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 8, "explicit false growth reaches the requested size")
      if (allocated(error)) return
      call check(error,.not. any(a(7:8)), "explicit false fill is honored")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 3, error=refused)
      call check(error,.not. allocated(refused), "unallocated default growth succeeds")
      if (allocated(error)) return
      call check(error, size(a), 3, "unallocated growth reaches the requested size")
      if (allocated(error)) return
      call check(error, all(a .eqv. .false.), "default fill covers every element")
      if (allocated(error)) return
      deallocate (a)
      call grow_array(a, 0, error=refused)
      call check(error,.not. allocated(refused), "zero request on an unallocated array is accepted")
      if (allocated(error)) return
      if (allocated(a)) deallocate (a)
      allocate (a(0))
      call grow_array(a, 2, error=refused)
      call check(error,.not. allocated(refused), "allocated empty array grows")
      if (allocated(error)) return
      call check(error, size(a), 2, "empty array grows to the requested size")
      if (allocated(error)) return
      call check(error, all(a .eqv. .false.), "empty growth fills every element")

   end subroutine test_grow_logical_1d_contracts

   !> Verify selection, stable order, prefix extent, empty masks, and allocation
   subroutine test_filter_logical_1d_contracts(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      logical, allocatable :: a(:), expected(:)
      logical :: keep(4)

      keep = [.true., .false., .true., .true.]
      call filter_array(a, 3, keep, 2)
      call check(error,.not. allocated(a), "unallocated filter is skipped")
      if (allocated(error)) return
      a = [.true., .false., .false., .true.]
      expected = a([1, 3])
      call filter_array(a, 3, keep, 2)
      call check(error, size(a, 1), 2, "filter uses nvalid extent")
      if (allocated(error)) return
      call check(error, all(a .eqv. expected), "filter preserves selected values and order")
      if (allocated(error)) return
      a = [.true., .false., .false., .true.]
      expected = a
      keep = .true.
      call filter_array(a, 4, keep, 4)
      call check(error, size(a, 1), 4, "all-true mask keeps the extent")
      if (allocated(error)) return
      call check(error, all(a .eqv. expected), "all-true mask keeps every value")
      if (allocated(error)) return
      keep = .false.
      call filter_array(a, 4, keep, 0)
      call check(error, allocated(a), "all-false mask keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 1), 0, "all-false mask gives zero extent")
      if (allocated(error)) return
      call filter_array(a, 0, keep, 0)
      call check(error, allocated(a), "empty filter keeps allocation")
      if (allocated(error)) return
      call check(error, size(a, 1), 0, "empty prefix remains empty")

   end subroutine test_filter_logical_1d_contracts

end module test_utils_mem
