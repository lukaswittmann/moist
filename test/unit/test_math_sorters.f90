!> Test suite for sorter utilities in moist_math_sorter
module test_math_sorters
   use mctc_env, only: wp
   use mctc_io_utils, only: to_string
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_sorter, only: qsort, counting_argsort
   implicit none(type, external)
   private

   public :: collect_math_sorters

   !> Numerical tolerance for floating-point comparisons in tests
   real(wp), parameter :: thr = 10.0_wp*epsilon(1.0_wp)

contains

   !> Collect all sorter tests
   subroutine collect_math_sorters(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("qsort_values_mixed", test_qsort_values_mixed), &
                  new_unittest("qsort_values_descending_long", test_qsort_values_descending_long), &
                  new_unittest("qsort_large_random_sorted_order", test_qsort_large_random_sorted_order), &
                  new_unittest("qsort_with_indices_tracks_permutation", test_qsort_with_indices_tracks_permutation), &
                  new_unittest("qsort_edge_sizes", test_qsort_edge_sizes), &
                  new_unittest("qsort_index_size_mismatch_error", test_qsort_index_size_mismatch_error), &
                  new_unittest("qsort_partition_companions", test_qsort_partition_companions), &
                  new_unittest("qsort_short_index_validation", test_qsort_short_index_validation), &
                  new_unittest("counting_argsort_stable", test_counting_argsort_stable), &
                  new_unittest("counting_argsort_edges", test_counting_argsort_edges) &
                  ]
   end subroutine collect_math_sorters

   !> Sort mixed values including negatives and duplicates
   subroutine test_qsort_values_mixed(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(9), expected(9)

      a = [3.0_wp, -1.0_wp, 2.0_wp, 2.0_wp, 0.0_wp, -5.0_wp, 4.0_wp, 1.0_wp, -1.0_wp]
      expected = [-5.0_wp, -1.0_wp, -1.0_wp, 0.0_wp, 1.0_wp, 2.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]

      call qsort(a, error=sort_error)

      call check(error,.not. allocated(sort_error), more="qsort returned an unexpected error")
      if (allocated(error)) return
      call check(error, all(a(1:size(a) - 1) <= a(2:size(a))), more="Array must be nondecreasing")
      if (allocated(error)) return
      call check(error, maxval(abs(a - expected)) < thr, more="Sorted values do not match expectation")
   end subroutine test_qsort_values_mixed

   !> Sort a longer descending array to exercise quicksort partitioning
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

      call check(error,.not. allocated(sort_error), more="qsort returned an unexpected error")
      if (allocated(error)) return
      call check(error, all(a(1:size(a) - 1) <= a(2:size(a))), more="Array must be nondecreasing")
      if (allocated(error)) return
      call check(error, maxval(abs(a - expected)) < thr, more="Descending input was not sorted correctly")
   end subroutine test_qsort_values_descending_long

   !> Sort 10k random values and verify neighboring order
   subroutine test_qsort_large_random_sorted_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a(10000)
      integer, allocatable :: seed(:)
      integer :: seed_size, i

      call random_seed(size=seed_size)
      allocate (seed(seed_size))
      do i = 1, seed_size
         seed(i) = 420 + 69*i
      end do
      call random_seed(put=seed)
      deallocate (seed)

      call random_number(a)
      a = 6.9_wp*a - 4.2_wp

      call qsort(a, error=sort_error)

      call check(error,.not. allocated(sort_error), more="qsort returned an unexpected error")
      if (allocated(error)) return

      do i = 1, size(a) - 1
         call check(error, a(i) <= a(i + 1), more="Large random array must be nondecreasing")
         if (allocated(error)) return
      end do
   end subroutine test_qsort_large_random_sorted_order

   !> Sort values with index tracking and verify value-index consistency
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

      call check(error,.not. allocated(sort_error), more="qsort with indices returned an unexpected error")
      if (allocated(error)) return
      call check(error, all(a(1:size(a) - 1) <= a(2:size(a))), more="Array must be nondecreasing")
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

   !> Verify empty, singleton, and two-element edge sizes
   subroutine test_qsort_edge_sizes(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a0(0), a1(1), a2(2)

      a1 = [42.0_wp]
      a2 = [2.0_wp, -1.0_wp]

      call qsort(a0, error=sort_error)
      call check(error,.not. allocated(sort_error), more="Empty array should not produce an error")
      if (allocated(error)) return

      call qsort(a1, error=sort_error)
      call check(error,.not. allocated(sort_error), more="Singleton array should not produce an error")
      if (allocated(error)) return
      call check(error, abs(a1(1) - 42.0_wp) < thr, more="Singleton array must remain unchanged")
      if (allocated(error)) return

      call qsort(a2, error=sort_error)
      call check(error,.not. allocated(sort_error), more="Two-element sort should not produce an error")
      if (allocated(error)) return
      call check(error, all(a2(1:size(a2) - 1) <= a2(2:size(a2))), more="Two-element array was not sorted")
      if (allocated(error)) return
      call check(error, abs(a2(1) + 1.0_wp) < thr .and. abs(a2(2) - 2.0_wp) < thr, &
                 more="Two-element array sorted values are incorrect")
   end subroutine test_qsort_edge_sizes

   !> Ensure mismatched index size reports an error
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

   !> Verify both sorter paths at the insertion threshold and with duplicate keys
   subroutine test_qsort_partition_companions(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp), allocatable :: a(:), original(:), plain(:)
      integer, allocatable :: ind(:)
      character(len=:), allocatable :: label
      integer :: n, i, j, pattern, isize
      !> Partitioning runs while hi - lo > insertion_cutoff (24): n = 24, 25 are
      !> insertion-only, n = 26 is the first partition, 64 and 128 recurse
      integer, parameter :: sizes(5) = [24, 25, 26, 64, 128]

      do isize = 1, size(sizes)
         n = sizes(isize)
         allocate (a(n), original(n), plain(n), ind(n))
         do pattern = 1, 3
            do i = 1, n
               select case (pattern)
               case (1)
                  original(i) = real(mod(37*i, 17) - 8, wp)
               case (2)
                  original(i) = real(n - i, wp)
               case default
                  original(i) = -2.0_wp
               end select
               ind(i) = 100 + i
            end do
            label = " (n = "//to_string(n)//", pattern "//to_string(pattern)//")"
            a = original
            plain = original
            call qsort(a, ind, sort_error)
            call check(error,.not. allocated(sort_error), "Indexed partition sort failed"//label)
            if (allocated(error)) return
            call check(error, all(a(:n - 1) <= a(2:)), "Indexed partition order is incorrect"//label)
            if (allocated(error)) return
            call check(error, all(ind >= 101 .and. ind <= 100 + n), "Companion labels out of range"//label)
            if (allocated(error)) return
            do i = 1, n
               call check(error, count(ind == 100 + i) == 1, "Companions must remain a permutation"//label)
               if (allocated(error)) return
               call check(error, abs(a(i) - original(ind(i) - 100)) < thr, "Companion pairing lost"//label)
               if (allocated(error)) return
            end do
            call qsort(plain, error=sort_error)
            call check(error,.not. allocated(sort_error), "Plain partition sort failed"//label)
            if (allocated(error)) return
            call check(error, all(plain(:n - 1) <= plain(2:)), "Plain partition order is incorrect"//label)
            if (allocated(error)) return
            do j = 1, n
               call check(error, count(abs(plain - original(j)) < thr) == count(abs(original - original(j)) < thr), &
                          "Plain sort must preserve the input multiset"//label)
               if (allocated(error)) return
            end do
         end do
         deallocate (a, original, plain, ind)
      end do
   end subroutine test_qsort_partition_companions

   !> Validate optional index sizes even for empty and singleton values
   subroutine test_qsort_short_index_validation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_error_type), allocatable :: sort_error
      real(wp) :: a0(0), a1(1)
      integer :: ind0(0), ind1(1)

      a1 = 3.0_wp
      ind1 = 77
      call qsort(a0, ind0, sort_error)
      call check(error,.not. allocated(sort_error), "Empty paired sort failed")
      if (allocated(error)) return
      call qsort(a1, ind1, sort_error)
      call check(error,.not. allocated(sort_error), "Singleton paired sort failed")
      if (allocated(error)) return
      call check(error, ind1(1) == 77 .and. abs(a1(1) - 3.0_wp) < thr, "Singleton pair changed")
      if (allocated(error)) return
      call qsort(a0, ind1, sort_error)
      call check(error, allocated(sort_error), "Empty values must still validate index size")
      if (allocated(error)) return
      call qsort(a1, ind0, sort_error)
      call check(error, allocated(sort_error), "Singleton values must still validate index size")
   end subroutine test_qsort_short_index_validation

   !> Check stable bucket permutation against an independent explicit answer
   subroutine test_counting_argsort_stable(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: buckets(10), perm(10), expected(10)

      buckets = [5, 0, 3, 5, 0, 1, 3, 1, 5, 0]
      ! Buckets in ascending order, original indices ascending within each:
      ! 0: {2, 5, 10}, 1: {6, 8}, 3: {3, 7}, 5: {1, 4, 9}
      expected = [2, 5, 10, 6, 8, 3, 7, 1, 4, 9]
      perm = -99
      call counting_argsort(buckets, 5, perm)
      call check(error, all(perm == expected), "Counting argsort order, permutation or stability lost")
   end subroutine test_counting_argsort_stable

   !> Check empty, singleton, and all-zero bucket boundary cases
   subroutine test_counting_argsort_edges(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer :: b0(0), p0(0), b1(1), p1(1), b4(4), p4(4)

      call counting_argsort(b0, 0, p0)
      call counting_argsort(b0, 7, p0)
      b1 = 7
      p1 = -99
      call counting_argsort(b1, 7, p1)
      call check(error, p1(1) == 1, "Singleton maximum bucket must map to index one")
      if (allocated(error)) return
      b4 = 0
      p4 = -99
      call counting_argsort(b4, 0, p4)
      call check(error, all(p4 == [1, 2, 3, 4]), "All-zero bucket sort must preserve stable order")
   end subroutine test_counting_argsort_edges

end module test_math_sorters
