!> Stable counting argsort of small non-negative integer keys
!>
!> - `perm(j)` is the position of the j-th entry in key order; equal keys keep
!>   their input order
!> - O(n + nkey) for keys in `0 .. size(counts) - 1`
!> - the counters are the caller's workspace (a fixed-size local array), so
!>   nothing is allocated per call
!> - stdlib has no counterpart: its `radix_sort` sorts values, not indices
module moist_math_sorter_counting
   implicit none(type, external)
   private

   public :: counting_argsort

contains

   !> Permutation that sorts entries by their integer keys, stably
   !>
   !> @param[in]  keys   Key of each entry, `0 .. size(counts) - 1`
   !> @param[out] perm   Entries in key order, `perm(1:size(keys))`
   !> @param[out] counts Workspace, one counter per key value, undefined on return
   pure subroutine counting_argsort(keys, perm, counts)
      !> Key of each entry
      integer, intent(in), contiguous :: keys(:)
      !> Entries in key order
      integer, intent(out), contiguous :: perm(:)
      !> Workspace, one counter per key value
      integer, intent(out), contiguous :: counts(0:)

      integer :: i, k, top, c, nxt

      if (size(perm) < size(keys)) error stop "counting_argsort: perm is smaller than keys"
      top = size(counts) - 1
      if (minval(keys) < 0 .or. maxval(keys) > top) then
         error stop "counting_argsort: key outside 0 .. size(counts) - 1"
      end if
      counts = 0
      do i = 1, size(keys)
         counts(keys(i)) = counts(keys(i)) + 1
      end do
      ! Exclusive prefix sum: slots before each key
      c = 0
      do k = 0, top
         nxt = c + counts(k)
         counts(k) = c
         c = nxt
      end do
      ! Scatter in input order, so equal keys keep it
      do i = 1, size(keys)
         counts(keys(i)) = counts(keys(i)) + 1
         perm(counts(keys(i))) = i
      end do
   end subroutine counting_argsort

end module moist_math_sorter_counting
