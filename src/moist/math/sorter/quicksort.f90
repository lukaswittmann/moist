!> Quicksort-based sorting utilities for real arrays.
module moist_math_sorter_quicksort
   use mctc_env, only: wp, error_type, fatal_error
   implicit none
   private

   !> Partition size threshold where insertion sort is used.
   integer, parameter :: insertion_cutoff = 24

   public :: qsort

contains

   !> Sort a real array in ascending order.
   !> @param[inout] a     Values to sort in place.
   !> @param[inout] ind   Optional index array permuted with `a`.
   !> @param[out]   error Error object set when input validation fails.
   subroutine qsort(a, ind, error)
      !> Values to sort in ascending order.
      real(wp), intent(inout), contiguous :: a(:)
      !> Optional index array tracking the permutation applied to `a`.
      integer, intent(inout), contiguous, optional :: ind(:)
      !> Error information for invalid inputs.
      type(error_type), allocatable, intent(out) :: error

      if (present(ind)) then
         if (size(ind) /= size(a)) then
            call fatal_error(error, "qsort: size(ind) must match size(a)")
            return
         end if
      end if

      if (size(a) <= 1) return

      if (present(ind)) then
         call sort_idx(a, ind, 1, size(a))
      else
         call sort_plain(a, 1, size(a))
      end if
   end subroutine qsort

   !> Quicksort core for values only (no coupled index array).
   !> @param[inout] a     Values to sort in place.
   !> @param[in]    left  Left bound (inclusive).
   !> @param[in]    right Right bound (inclusive).
   recursive subroutine sort_plain(a, left, right)
      !> Values to sort in place.
      real(wp), intent(inout) :: a(:)
      !> Left bound of the active partition (inclusive).
      integer, intent(in)    :: left, right

      !> Current left/right bounds for tail-recursive partition processing.
      integer  :: lo, hi, i, j
      !> Pivot value used for partitioning.
      real(wp) :: pivot

      lo = left
      hi = right

      do while (hi - lo > insertion_cutoff)
         pivot = a(median3(a, lo, hi))

         i = lo
         j = hi
         do
            do while (a(i) < pivot)
               i = i + 1
            end do
            do while (pivot < a(j))
               j = j - 1
            end do
            if (i > j) exit
            call swap_r(a(i), a(j))
            i = i + 1
            j = j - 1
         end do

         if (j - lo < hi - i) then
            if (lo < j) call sort_plain(a, lo, j)
            lo = i
         else
            if (i < hi) call sort_plain(a, i, hi)
            hi = j
         end if
      end do

      call insertion_plain(a, lo, hi)
   end subroutine sort_plain

   !> Quicksort core that keeps an index array synchronized with values.
   !> @param[inout] a     Values to sort in place.
   !> @param[inout] ind   Index array permuted with `a`.
   !> @param[in]    left  Left bound (inclusive).
   !> @param[in]    right Right bound (inclusive).
   recursive subroutine sort_idx(a, ind, left, right)
      !> Values to sort in place.
      real(wp), intent(inout) :: a(:)
      !> Index array permuted alongside `a`.
      integer, intent(inout) :: ind(:)
      !> Left and right bounds of the active partition (inclusive).
      integer, intent(in)    :: left, right

      !> Current left/right bounds and scan indices for partitioning.
      integer  :: lo, hi, i, j
      !> Pivot value used for partitioning.
      real(wp) :: pivot

      lo = left
      hi = right

      do while (hi - lo > insertion_cutoff)
         pivot = a(median3(a, lo, hi))

         i = lo
         j = hi
         do
            do while (a(i) < pivot)
               i = i + 1
            end do
            do while (pivot < a(j))
               j = j - 1
            end do
            if (i > j) exit
            call swap_pair(a(i), a(j), ind(i), ind(j))
            i = i + 1
            j = j - 1
         end do

         if (j - lo < hi - i) then
            if (lo < j) call sort_idx(a, ind, lo, j)
            lo = i
         else
            if (i < hi) call sort_idx(a, ind, i, hi)
            hi = j
         end if
      end do

      call insertion_idx(a, ind, lo, hi)
   end subroutine sort_idx

   !> Insertion sort on a bounded segment of a real array.
   !> @param[inout] a     Values to sort in place.
   !> @param[in]    left  Left bound (inclusive).
   !> @param[in]    right Right bound (inclusive).
   subroutine insertion_plain(a, left, right)
      !> Values to sort in place.
      real(wp), intent(inout) :: a(:)
      !> Left and right bounds of the insertion-sort segment.
      integer, intent(in)    :: left, right

      !> Loop index and backward scan index.
      integer  :: i, j
      !> Candidate element inserted into the sorted prefix.
      real(wp) :: key

      do i = left + 1, right
         key = a(i)
         j = i - 1
         do while (j >= left)
            if (.not. (key < a(j))) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = key
      end do
   end subroutine insertion_plain

   !> Insertion sort on values while keeping indices synchronized.
   !> @param[inout] a     Values to sort in place.
   !> @param[inout] ind   Index array permuted with `a`.
   !> @param[in]    left  Left bound (inclusive).
   !> @param[in]    right Right bound (inclusive).
   subroutine insertion_idx(a, ind, left, right)
      !> Values to sort in place.
      real(wp), intent(inout) :: a(:)
      !> Index array permuted alongside `a`.
      integer, intent(inout) :: ind(:)
      !> Left and right bounds of the insertion-sort segment.
      integer, intent(in)    :: left, right

      !> Loop index, backward scan index, and saved index entry.
      integer  :: i, j, ikey
      !> Candidate element inserted into the sorted prefix.
      real(wp) :: key

      do i = left + 1, right
         key = a(i)
         ikey = ind(i)
         j = i - 1
         do while (j >= left)
            if (.not. (key < a(j))) exit
            a(j + 1) = a(j)
            ind(j + 1) = ind(j)
            j = j - 1
         end do
         a(j + 1) = key
         ind(j + 1) = ikey
      end do
   end subroutine insertion_idx

   !> Return the index of the median among left, middle, and right values.
   !> @param[in]  a     Input array.
   !> @param[in]  left  Left bound (inclusive).
   !> @param[in]  right Right bound (inclusive).
   !> @return           Index of the median-of-three pivot candidate.
   pure integer function median3(a, left, right) result(idx)
      !> Input array.
      real(wp), intent(in) :: a(:)
      !> Left and right bounds used for median-of-three selection.
      integer, intent(in) :: left, right

      !> Middle index between `left` and `right`.
      integer :: mid

      mid = left + (right - left)/2

      if (a(left) < a(mid)) then
         if (a(mid) < a(right)) then
            idx = mid
         else if (a(left) < a(right)) then
            idx = right
         else
            idx = left
         end if
      else
         if (a(left) < a(right)) then
            idx = left
         else if (a(mid) < a(right)) then
            idx = right
         else
            idx = mid
         end if
      end if
   end function median3

   !> Swap two real values.
   !> @param[inout] x First value.
   !> @param[inout] y Second value.
   pure subroutine swap_r(x, y)
      !> First value to exchange.
      real(wp), intent(inout) :: x, y
      !> Temporary storage for swapping.
      real(wp) :: tmp

      tmp = x
      x = y
      y = tmp
   end subroutine swap_r

   !> Swap two values and their associated indices.
   !> @param[inout] x  First value.
   !> @param[inout] y  Second value.
   !> @param[inout] ix Index associated with `x`.
   !> @param[inout] iy Index associated with `y`.
   pure subroutine swap_pair(x, y, ix, iy)
      !> Values to exchange.
      real(wp), intent(inout) :: x, y
      !> Index entries to exchange with their corresponding values.
      integer, intent(inout) :: ix, iy
      !> Temporary real storage for swapping.
      real(wp) :: tx
      !> Temporary integer storage for swapping.
      integer  :: ti

      tx = x
      x = y
      y = tx

      ti = ix
      ix = iy
      iy = ti
   end subroutine swap_pair

end module moist_math_sorter_quicksort
