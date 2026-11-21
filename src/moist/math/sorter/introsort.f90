! SPDX-License-Identifier: MPL-2.0 AND MIT
! Ported from fortran-lang/stdlib v0.8.1 (commit 0ede3015), src/sorting/stdlib_sorting_sort.fypp;
! MIT, Copyright (c) 2019-2021 stdlib contributors, license text in LICENSE next to this
! file. Introsort of D. R. Musser, included in stdlib with his permission on the condition that
! it is cited: Musser, D.R., "Introspective Sorting and Selection Algorithms",
! Software--Practice and Experience, Vol. 27(8), 983-993 (August 1997)
!
!> Unstable in-place sort, stdlib's `sort` (introsort)
!>
!> - quicksort with a median-of-three pivot, insertion sort on partitions of at
!>   most 16 entries, heap sort once the recursion exceeds `2 floor(log2 n)`
!> - O(n log n) worst case, no scratch, nothing allocated
!> - kinds: default integer and `real(wp)`; with a NaN the order is undefined
!>
!> Changes from stdlib: kinds instantiated by hand, default-integer indices,
!> `contiguous` dummies, depth limit from the bit length of `n` instead of a
!> floating-point logarithm (no `log(0)`), iterative sift-down, `reverse` as a
!> reversal after the ascending sort
module moist_math_sorter_introsort
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: sort

   !> Partitions of at most this many entries are insertion-sorted
   integer, parameter :: insert_size = 16

   !> Unstable sort in place, non-decreasing or with `reverse` non-increasing
   !>
   !> `call sort(array[, reverse])`
   interface sort
      module procedure :: sort_int
      module procedure :: sort_real
   end interface sort

contains

   !> Recursion depth before the heap sort takes over: `2 floor(log2 n)`
   !>
   !> @param[in] n Number of entries
   !> @returns     Depth limit, 0 for `n < 2`
   pure function depth_limit(n) result(limit)
      !> Number of entries
      integer, intent(in) :: n
      !> Depth limit
      integer :: limit

      limit = 0
      if (n >= 2) limit = 2*(bit_size(n) - 1 - leadz(n))
   end function depth_limit

   !* ================================================================================= *!
   !*                                   Default integer                                 *!
   !* ================================================================================= *!

   !> Unstable sort of an integer array
   !>
   !> @param[inout] array   Values, sorted on return
   !> @param[in]    reverse Non-increasing order
   pure subroutine sort_int(array, reverse)
      !> Values, sorted on return
      integer, intent(inout), contiguous :: array(0:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      integer :: lo, hi, tmp

      call introsort_int(array, depth_limit(size(array)))
      if (.not. present(reverse)) return
      if (.not. reverse) return
      lo = 0
      hi = size(array) - 1
      do while (lo < hi)
         tmp = array(lo)
         array(lo) = array(hi)
         array(hi) = tmp
         lo = lo + 1
         hi = hi - 1
      end do
   end subroutine sort_int

   !> Quicksort, insertion sort below `insert_size`, heap sort at the depth limit
   !>
   !> @param[inout] array Values
   !> @param[in]    limit Remaining recursion depth
   pure recursive subroutine introsort_int(array, limit)
      !> Values
      integer, intent(inout), contiguous :: array(0:)
      !> Remaining recursion depth
      integer, intent(in) :: limit

      integer :: pivot

      if (size(array) <= insert_size) then
         call insertion_int(array)
      else if (limit == 0) then
         call heap_sort_int(array)
      else
         call partition_int(array, pivot)
         call introsort_int(array(0:pivot - 1), limit - 1)
         call introsort_int(array(pivot + 1:), limit - 1)
      end if
   end subroutine introsort_int

   !> Lomuto partition about the median of the first, middle and last entry
   !>
   !> @param[inout] array Values, at least two
   !> @param[out]   pivot Final position of the pivot
   pure subroutine partition_int(array, pivot)
      !> Values
      integer, intent(inout), contiguous :: array(0:)
      !> Final position of the pivot
      integer, intent(out) :: pivot

      integer :: u, v, w, x, y, i, j, n

      n = size(array)
      ! Median of three, moved to the end
      u = array(0)
      v = array(n/2 - 1)
      w = array(n - 1)
      if ((u > v) .neqv. (u > w)) then
         x = u
         array(0) = array(n - 1)
         array(n - 1) = u
      else if ((v < u) .neqv. (v < w)) then
         x = v
         array(n/2 - 1) = array(n - 1)
         array(n - 1) = v
      else
         x = w
      end if
      i = -1
      do j = 0, n - 2
         if (array(j) <= x) then
            i = i + 1
            y = array(i)
            array(i) = array(j)
            array(j) = y
         end if
      end do
      y = array(i + 1)
      array(i + 1) = array(n - 1)
      array(n - 1) = y
      pivot = i + 1
   end subroutine partition_int

   !> Insertion sort
   !>
   !> @param[inout] array Values
   pure subroutine insertion_int(array)
      !> Values
      integer, intent(inout), contiguous :: array(0:)

      integer :: i, j, key

      do j = 1, size(array) - 1
         key = array(j)
         i = j - 1
         do while (i >= 0)
            if (array(i) <= key) exit
            array(i + 1) = array(i)
            i = i - 1
         end do
         array(i + 1) = key
      end do
   end subroutine insertion_int

   !> Heap sort on a max-heap
   !>
   !> @param[inout] array Values
   pure subroutine heap_sort_int(array)
      !> Values
      integer, intent(inout), contiguous :: array(0:)

      integer :: i, heap_size, y

      heap_size = size(array)
      do i = (heap_size - 2)/2, 0, -1
         call max_heapify_int(array, i, heap_size)
      end do
      do i = heap_size - 1, 1, -1
         y = array(0)
         array(0) = array(i)
         array(i) = y
         call max_heapify_int(array, 0, i)
      end do
   end subroutine heap_sort_int

   !> Sift `array(root)` down the max-heap `array(0:heap_size-1)`
   !>
   !> @param[inout] array     Heap
   !> @param[in]    root      Entry to sift down
   !> @param[in]    heap_size Entries of the heap
   pure subroutine max_heapify_int(array, root, heap_size)
      !> Heap
      integer, intent(inout), contiguous :: array(0:)
      !> Entry to sift down
      integer, intent(in) :: root
      !> Entries of the heap
      integer, intent(in) :: heap_size

      integer :: i, l, r, largest, y

      i = root
      do
         largest = i
         l = 2*i + 1
         r = l + 1
         if (l < heap_size) then
            if (array(l) > array(largest)) largest = l
         end if
         if (r < heap_size) then
            if (array(r) > array(largest)) largest = r
         end if
         if (largest == i) exit
         y = array(i)
         array(i) = array(largest)
         array(largest) = y
         i = largest
      end do
   end subroutine max_heapify_int

   !* ================================================================================= *!
   !*                                      real(wp)                                     *!
   !* ================================================================================= *!

   !> Unstable sort of a real array
   !>
   !> @param[inout] array   Values, sorted on return
   !> @param[in]    reverse Non-increasing order
   pure subroutine sort_real(array, reverse)
      !> Values, sorted on return
      real(wp), intent(inout), contiguous :: array(0:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      integer :: lo, hi
      real(wp) :: tmp

      call introsort_real(array, depth_limit(size(array)))
      if (.not. present(reverse)) return
      if (.not. reverse) return
      lo = 0
      hi = size(array) - 1
      do while (lo < hi)
         tmp = array(lo)
         array(lo) = array(hi)
         array(hi) = tmp
         lo = lo + 1
         hi = hi - 1
      end do
   end subroutine sort_real

   !> Quicksort, insertion sort below `insert_size`, heap sort at the depth limit
   !>
   !> @param[inout] array Values
   !> @param[in]    limit Remaining recursion depth
   pure recursive subroutine introsort_real(array, limit)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)
      !> Remaining recursion depth
      integer, intent(in) :: limit

      integer :: pivot

      if (size(array) <= insert_size) then
         call insertion_real(array)
      else if (limit == 0) then
         call heap_sort_real(array)
      else
         call partition_real(array, pivot)
         call introsort_real(array(0:pivot - 1), limit - 1)
         call introsort_real(array(pivot + 1:), limit - 1)
      end if
   end subroutine introsort_real

   !> Lomuto partition about the median of the first, middle and last entry
   !>
   !> @param[inout] array Values, at least two
   !> @param[out]   pivot Final position of the pivot
   pure subroutine partition_real(array, pivot)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)
      !> Final position of the pivot
      integer, intent(out) :: pivot

      real(wp) :: u, v, w, x, y
      integer :: i, j, n

      n = size(array)
      ! Median of three, moved to the end
      u = array(0)
      v = array(n/2 - 1)
      w = array(n - 1)
      if ((u > v) .neqv. (u > w)) then
         x = u
         array(0) = array(n - 1)
         array(n - 1) = u
      else if ((v < u) .neqv. (v < w)) then
         x = v
         array(n/2 - 1) = array(n - 1)
         array(n - 1) = v
      else
         x = w
      end if
      i = -1
      do j = 0, n - 2
         if (array(j) <= x) then
            i = i + 1
            y = array(i)
            array(i) = array(j)
            array(j) = y
         end if
      end do
      y = array(i + 1)
      array(i + 1) = array(n - 1)
      array(n - 1) = y
      pivot = i + 1
   end subroutine partition_real

   !> Insertion sort
   !>
   !> @param[inout] array Values
   pure subroutine insertion_real(array)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)

      integer :: i, j
      real(wp) :: key

      do j = 1, size(array) - 1
         key = array(j)
         i = j - 1
         do while (i >= 0)
            if (array(i) <= key) exit
            array(i + 1) = array(i)
            i = i - 1
         end do
         array(i + 1) = key
      end do
   end subroutine insertion_real

   !> Heap sort on a max-heap
   !>
   !> @param[inout] array Values
   pure subroutine heap_sort_real(array)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)

      integer :: i, heap_size
      real(wp) :: y

      heap_size = size(array)
      do i = (heap_size - 2)/2, 0, -1
         call max_heapify_real(array, i, heap_size)
      end do
      do i = heap_size - 1, 1, -1
         y = array(0)
         array(0) = array(i)
         array(i) = y
         call max_heapify_real(array, 0, i)
      end do
   end subroutine heap_sort_real

   !> Sift `array(root)` down the max-heap `array(0:heap_size-1)`
   !>
   !> @param[inout] array     Heap
   !> @param[in]    root      Entry to sift down
   !> @param[in]    heap_size Entries of the heap
   pure subroutine max_heapify_real(array, root, heap_size)
      !> Heap
      real(wp), intent(inout), contiguous :: array(0:)
      !> Entry to sift down
      integer, intent(in) :: root
      !> Entries of the heap
      integer, intent(in) :: heap_size

      integer :: i, l, r, largest
      real(wp) :: y

      i = root
      do
         largest = i
         l = 2*i + 1
         r = l + 1
         if (l < heap_size) then
            if (array(l) > array(largest)) largest = l
         end if
         if (r < heap_size) then
            if (array(r) > array(largest)) largest = r
         end if
         if (largest == i) exit
         y = array(i)
         array(i) = array(largest)
         array(largest) = y
         i = largest
      end do
   end subroutine max_heapify_real

end module moist_math_sorter_introsort
