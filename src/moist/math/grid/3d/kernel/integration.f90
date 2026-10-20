!> Compensated scalar quadrature kernels, independent of grids and threading
module moist_math_grid_3d_kernel_integration
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: compensated_sum, compensated_weighted_sum, merge_compensated_sums
   public :: compensated_block_sum

contains

   !> Kahan sum of f(first:last), retaining the rounding correction
   !>
   !> - total is the Kahan running sum
   !> - correction is carried into the next addition
   !> - Empty ranges return zero for both outputs
   !>
   !> @param[in] f Tabulated field values
   !> @param[in] first First field index
   !> @param[in] last Last field index, possibly less than first
   !> @param[out] total Running sum
   !> @param[out] correction Running rounding correction
   pure subroutine compensated_sum(f, first, last, total, correction)
      !> Tabulated field values
      real(wp), intent(in) :: f(:)
      !> First and last field indices
      integer, intent(in) :: first, last
      !> Running sum and rounding correction
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

   !> Kahan sum of weights(i)*f(i), without a whole-field temporary
   !>
   !> @param[in] f Tabulated field values
   !> @param[in] first First field index
   !> @param[in] last Last field index, possibly less than first
   !> @param[in] weights Final volume weights, indexed like f
   !> @param[out] total Running sum
   !> @param[out] correction Running rounding correction
   pure subroutine compensated_weighted_sum(f, first, last, weights, total, correction)
      !> Tabulated field values
      real(wp), intent(in) :: f(:)
      !> First and last field indices
      integer, intent(in) :: first, last
      !> Final volume weights
      real(wp), intent(in) :: weights(:)
      !> Running sum and rounding correction
      real(wp), intent(out) :: total, correction
      real(wp) :: term, updated
      integer :: i

      total = 0.0_wp
      correction = 0.0_wp
      do i = first, last
         term = weights(i)*f(i) - correction
         updated = total + term
         correction = (updated - total) - term
         total = updated
      end do
   end subroutine compensated_weighted_sum

   !> Neumaier block sum, retaining both components even with strong cancellation
   !>
   !> - Corrected value: total + correction
   !> - Optional weights indexed like f, including bounded field slices
   !>
   !> @param[in] f Tabulated field values
   !> @param[in] first First field index
   !> @param[in] last Last field index, possibly less than first
   !> @param[out] total Running sum
   !> @param[out] correction Running rounding correction
   !> @param[in] weights Optional final volume weights, indexed like f
   pure subroutine compensated_block_sum(f, first, last, total, correction, weights)
      !> Tabulated field values
      real(wp), intent(in) :: f(:)
      !> First and last field indices
      integer, intent(in) :: first, last
      !> Running sum and rounding correction
      real(wp), intent(out) :: total, correction
      !> Optional final volume weights
      real(wp), intent(in), optional :: weights(:)
      integer :: i

      total = 0.0_wp
      correction = 0.0_wp
      if (present(weights)) then
         do i = first, last
            call add_compensated(weights(i)*f(i), total, correction)
         end do
      else
         do i = first, last
            call add_compensated(f(i), total, correction)
         end do
      end if
   end subroutine compensated_block_sum

   !> Merge (total, correction) blocks in index order with Neumaier compensation
   !>
   !> - Keeping both components prevents cancellation between large block totals
   !>   from discarding the small contributions retained inside each block
   !>
   !> @param[in] partial (total, correction) per block, (2, nblocks)
   pure function merge_compensated_sums(partial) result(total)
      !> Block totals and corrections
      real(wp), intent(in) :: partial(:, :)
      !> Merged sum, and running correction
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
   !>
   !> @param[in] value Term to add
   !> @param[inout] total Running sum
   !> @param[inout] correction Running rounding correction
   pure subroutine add_compensated(value, total, correction)
      !> Term to add
      real(wp), intent(in) :: value
      !> Running sum and rounding correction
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
