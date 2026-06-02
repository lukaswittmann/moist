!> Test suite for sorter utilities in moist_math_sorter
module test_math_sorters
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_sorter, only: qsort
   implicit none (type, external)
   private

   public :: collect_math_sorters

   !> Numerical tolerance for floating-point comparisons in tests.
   real(wp), parameter :: thr = 10.0_wp * epsilon(1.0_wp)

contains

   !> Collect all sorter tests.
   subroutine collect_math_sorters(testsuite)
      !> Collection of tests.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("qsort_values_mixed", test_qsort_values_mixed), &
         new_unittest("qsort_values_descending_long", test_qsort_values_descending_long), &
         new_unittest("qsort_large_random_sorted_order", test_qsort_large_random_sorted_order), &
         new_unittest("qsort_with_indices_tracks_permutation", test_qsort_with_indices_tracks_permutation), &
         new_unittest("qsort_edge_sizes", test_qsort_edge_sizes), &
         new_unittest("qsort_index_size_mismatch_error", test_qsort_index_size_mismatch_error) &
      ]
   end subroutine collect_math_sorters

   !> Sort mixed values including negatives and duplicates.
   subroutine test_qsort_values_mixed(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(9), expected(9)

      a = [3.0_wp, -1.0_wp, 2.0_wp, 2.0_wp, 0.0_wp, -5.0_wp, 4.0_wp, 1.0_wp, -1.0_wp]
      expected = [-5.0_wp, -1.0_wp, -1.0_wp, 0.0_wp, 1.0_wp, 2.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]

      call qsort(a, error=sort_error)

      call check(error, .not. allocated(sort_error), more="qsort returned an unexpected error")
      if (allocated(error)) return
      call check(error, all(a(1:size(a)-1) <= a(2:size(a))), more="Array must be nondecreasing")
      if (allocated(error)) return
      call check(error, maxval(abs(a - expected)) < thr, more="Sorted values do not match expectation")
   end subroutine test_qsort_values_mixed

   !> Sort a longer descending array to exercise quicksort partitioning.
   subroutine test_qsort_values_descending_long(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(64), expected(64)
      integer :: i

      do i = 1, 64
         a(i) = real(65 - i, wp)
         expected(i) = real(i, wp)
      end do

      call qsort(a, error=sort_error)

      call check(error, .not. allocated(sort_error), more="qsort returned an unexpected error")
      if (allocated(error)) return
      call check(error, all(a(1:size(a)-1) <= a(2:size(a))), more="Array must be nondecreasing")
      if (allocated(error)) return
      call check(error, maxval(abs(a - expected)) < thr, more="Descending input was not sorted correctly")
   end subroutine test_qsort_values_descending_long

   !> Sort 10k random values and verify neighboring order.
   subroutine test_qsort_large_random_sorted_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(10000)
      integer, allocatable :: seed(:)
      integer :: seed_size, i

      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      do i = 1, seed_size
         seed(i) = 420 + 69 * i
      end do
      call random_seed(put=seed)
      deallocate(seed)

      call random_number(a)
      a = 6.9_wp * a - 4.2_wp

      call qsort(a, error=sort_error)

      call check(error, .not. allocated(sort_error), more="qsort returned an unexpected error")
      if (allocated(error)) return

      do i = 1, size(a) - 1
         call check(error, a(i) <= a(i+1), more="Large random array must be nondecreasing")
         if (allocated(error)) return
      end do
   end subroutine test_qsort_large_random_sorted_order

   !> Sort values with index tracking and verify value-index consistency.
   subroutine test_qsort_with_indices_tracks_permutation(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(8), original(8)
      integer :: ind(8)
      logical :: seen(8)
      integer :: i

      original = [2.5_wp, -4.0_wp, 1.0_wp, 7.25_wp, 0.0_wp, -1.5_wp, 3.0_wp, 6.5_wp]
      a = original
      ind = [(i, i=1, size(ind))]

      call qsort(a, ind, sort_error)

      call check(error, .not. allocated(sort_error), more="qsort with indices returned an unexpected error")
      if (allocated(error)) return
      call check(error, all(a(1:size(a)-1) <= a(2:size(a))), more="Array must be nondecreasing")
      if (allocated(error)) return
      call check(error, all(ind >= 1 .and. ind <= size(ind)), more="Indices must remain within bounds")
      if (allocated(error)) return

      seen = .false.
      do i = 1, size(ind)
         if (seen(ind(i))) then
            call check(error, .false., "Index array must be a permutation")
            return
         end if
         seen(ind(i)) = .true.
         call check(error, abs(a(i) - original(ind(i))) < thr, &
            "Sorted value and tracked index are inconsistent")
         if (allocated(error)) return
      end do

      call check(error, all(seen), more="Index array must contain each original position exactly once")
   end subroutine test_qsort_with_indices_tracks_permutation

   !> Verify empty, singleton, and two-element edge sizes.
   subroutine test_qsort_edge_sizes(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a0(0), a1(1), a2(2)

      a1 = [42.0_wp]
      a2 = [2.0_wp, -1.0_wp]

      call qsort(a0, error=sort_error)
      call check(error, .not. allocated(sort_error), more="Empty array should not produce an error")
      if (allocated(error)) return

      call qsort(a1, error=sort_error)
      call check(error, .not. allocated(sort_error), more="Singleton array should not produce an error")
      if (allocated(error)) return
      call check(error, abs(a1(1) - 42.0_wp) < thr, more="Singleton array must remain unchanged")
      if (allocated(error)) return

      call qsort(a2, error=sort_error)
      call check(error, .not. allocated(sort_error), more="Two-element sort should not produce an error")
      if (allocated(error)) return
      call check(error, all(a2(1:size(a2)-1) <= a2(2:size(a2))), more="Two-element array was not sorted")
      if (allocated(error)) return
      call check(error, abs(a2(1) + 1.0_wp) < thr .and. abs(a2(2) - 2.0_wp) < thr, &
         more="Two-element array sorted values are incorrect")
   end subroutine test_qsort_edge_sizes

   !> Ensure mismatched index size reports an error.
   subroutine test_qsort_index_size_mismatch_error(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(4)
      integer :: ind(3)

      a = [4.0_wp, 1.0_wp, 3.0_wp, 2.0_wp]
      ind = [1, 2, 3]

      call qsort(a, ind, sort_error)

      call check(error, allocated(sort_error), more="Mismatched index size should produce an error")
   end subroutine test_qsort_index_size_mismatch_error

end module test_math_sorters
