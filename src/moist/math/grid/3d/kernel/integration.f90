!> Compensated scalar quadrature kernels, independent of grids and threading
module moist_math_grid_3d_kernel_integration
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: integration_weight, compensated_sum, compensated_weighted_sum, merge_compensated_sums
   public :: compensated_block_sum

   abstract interface
      !> Integration weight at an absolute field index; context is borrowed
      pure function integration_weight(context, i) result(weight)
         import :: wp
         implicit none(type, external)
         class(*), intent(in) :: context
         integer, intent(in) :: i
         real(wp) :: weight
      end function integration_weight
   end interface

contains

   !> Kahan sum of f(first:last), retaining the rounding correction
   !>
   !> total is the Kahan running sum; correction is carried into the next addition.
   !> Empty ranges return zero for both outputs.
   pure subroutine compensated_sum(f, first, last, total, correction)
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: first, last
      real(wp), intent(out) :: total, correction
      real(wp) :: term, updated
      integer :: i

      total = 0.0_wp
      correction = 0.0_wp
      do i = first, last
         term = f(i) - correction
         updated = total + term
         correction = (updated - total) - term
         total = updated
      end do
   end subroutine compensated_sum

   !> Kahan sum of weight(context,i)*f(i), without a whole-field temporary
   pure subroutine compensated_weighted_sum(f, first, last, weight, context, total, correction)
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: first, last
      procedure(integration_weight) :: weight
      class(*), intent(in) :: context
      real(wp), intent(out) :: total, correction
      real(wp) :: term, updated
      integer :: i

      total = 0.0_wp
      correction = 0.0_wp
      do i = first, last
         term = weight(context, i)*f(i) - correction
         updated = total + term
         correction = (updated - total) - term
         total = updated
      end do
   end subroutine compensated_weighted_sum

   !> Neumaier block sum, retaining both components even with strong cancellation
   !>
   !> The corrected value is total + correction. Weighted callers supply both
   !> optional arguments; unweighted callers omit both. offset shifts the weight
   !> index when f contains a slice of a larger grid field.
   pure subroutine compensated_block_sum(f, first, last, total, correction, weight, context, offset)
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: first, last
      real(wp), intent(out) :: total, correction
      procedure(integration_weight), optional :: weight
      class(*), intent(in), optional :: context
      integer, intent(in), optional :: offset
      integer :: i, origin

      total = 0.0_wp
      correction = 0.0_wp
      origin = 0
      if (present(offset)) origin = offset
      if (present(weight)) then
         do i = first, last
            call add_compensated(weight(context, i + origin)*f(i), total, correction)
         end do
      else
         do i = first, last
            call add_compensated(f(i), total, correction)
         end do
      end if
   end subroutine compensated_block_sum

   !> Merge (total, correction) blocks in index order with Neumaier compensation
   !>
   !> Keeping both components prevents cancellation between large block totals
   !> from discarding the small contributions retained inside each block.
   pure function merge_compensated_sums(partial) result(total)
      real(wp), intent(in) :: partial(:, :)
      real(wp) :: total, correction
      integer :: b, component

      total = 0.0_wp
      correction = 0.0_wp
      do b = 1, size(partial, 2)
         do component = 1, 2
            call add_compensated(partial(component, b), total, correction)
         end do
      end do
      total = total + correction
   end function merge_compensated_sums

   !> Neumaier update: track rounding when either addend dominates
   pure subroutine add_compensated(value, total, correction)
      real(wp), intent(in) :: value
      real(wp), intent(inout) :: total, correction
      real(wp) :: updated

      updated = total + value
      if (abs(total) >= abs(value)) then
         correction = correction + ((total - updated) + value)
      else
         correction = correction + ((value - updated) + total)
      end if
      total = updated
   end subroutine add_compensated

end module moist_math_grid_3d_kernel_integration
