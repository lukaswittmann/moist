!> Counting-sort-based argsort for integer bucket keys.
!>
!> Produces a permutation array `perm` such that `keys(perm(1)) <= keys(perm(2)) <= ...`
!> in a single O(N) pass. Adapted from the int8 radix sort pattern in
!> fortran-lang/stdlib (which degenerates to counting sort for 8-bit keys).
!>
!> Designed for spatial bucket indices (0..n_buckets-1) where n_buckets is
!> small (typically <= 256). The sort is stable: atoms with the same bucket
!> appear in perm in their original order.
module moist_math_sorter_counting_sort
   implicit none
   private

   public :: counting_argsort

contains

   !> Build a permutation that sorts atoms by their integer bucket keys.
   !>
   !> After the call, `perm(j)` is the original index of the j-th atom in
   !> sorted order. Bucket values must lie in `[0, max_bucket]`.
   !>
   !> @param[in]  buckets     Per-element bucket index, range [0, max_bucket]
   !> @param[in]  max_bucket  Largest valid bucket value
   !> @param[out] perm        Output permutation array (same size as buckets)
   pure subroutine counting_argsort(buckets, max_bucket, perm)
      !> Per-element spatial bucket (0-based)
      integer, intent(in) :: buckets(:)
      !> Upper bound of bucket range (inclusive)
      integer, intent(in) :: max_bucket
      !> Output permutation: perm(j) = original index of the j-th sorted element
      integer, intent(out) :: perm(:)

      integer :: i, j, n
      integer, allocatable :: counts(:), offsets(:)

      n = size(buckets)
      allocate (counts(0:max_bucket), offsets(0:max_bucket))

      ! Count occurrences per bucket
      counts = 0
      do i = 1, n
         counts(buckets(i)) = counts(buckets(i)) + 1
      end do

      ! Prefix sum: counts -> scatter offsets
      offsets(0) = 0
      do j = 1, max_bucket
         offsets(j) = offsets(j - 1) + counts(j - 1)
      end do

      ! Scatter: build permutation in sorted order (stable)
      do i = 1, n
         offsets(buckets(i)) = offsets(buckets(i)) + 1
         perm(offsets(buckets(i))) = i
      end do
   end subroutine counting_argsort

end module moist_math_sorter_counting_sort
