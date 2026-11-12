!> Allocation-free in-place heapsort for integer arrays
module moist_math_sorter_heapsort
   implicit none(type, external)
   private

   public :: heapsort

contains

   !> In-place ascending heapsort, allocation-free
   !>
   !> @param[in,out] a  values to sort
   pure subroutine heapsort(a)
      !> Values to sort
      integer, intent(inout) :: a(:)

      integer :: n, i, tmp

      n = size(a)
      do i = n/2, 1, -1
         call sift_down(a, i, n)
      end do
      do i = n, 2, -1
         tmp = a(1)
         a(1) = a(i)
         a(i) = tmp
         call sift_down(a, 1, i - 1)
      end do
   end subroutine heapsort

   !> Restore the max-heap property below `root` within `a(1:last)`
   !>
   !> @param[in,out] a     heap
   !> @param[in]     root  root index
   !> @param[in]     last  last heap index
   pure subroutine sift_down(a, root, last)
      !> Heap
      integer, intent(inout) :: a(:)
      !> Root index
      integer, intent(in) :: root
      !> Last heap index
      integer, intent(in) :: last

      integer :: i, child, tmp

      i = root
      do
         child = 2*i
         if (child > last) exit
         if (child < last) then
            if (a(child + 1) > a(child)) child = child + 1
         end if
         if (a(i) >= a(child)) exit
         tmp = a(i)
         a(i) = a(child)
         a(child) = tmp
         i = child
      end do
   end subroutine sift_down

end module moist_math_sorter_heapsort
