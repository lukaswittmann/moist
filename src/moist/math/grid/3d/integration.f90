!> Fixed-block volume quadrature with compensated, reproducible accumulation
module moist_math_grid_3d_integration
   use mctc_env, only: wp
   use moist_math_grid_3d_kernel_integration, only: compensated_sum, compensated_weighted_sum, &
      & compensated_block_sum, merge_compensated_sums
!$ use omp_lib, only: omp_in_parallel
   implicit none(type, external)
   private

   public :: integrate_values, integrate_weighted_values
   public :: integration_block_size, integration_parallel_threshold

   !> Blocks do not depend on the worker count, including serial large sums
   integer, parameter :: block_size = 8192
   !> Avoid parallel-region and partial-array overhead for small fields
   integer, parameter :: parallel_threshold = 65536
   !> Shared partition for callers which produce fields in bounded scratch
   integer, parameter :: integration_block_size = block_size
   !> Shared field length from which blocks run in parallel
   integer, parameter :: integration_parallel_threshold = parallel_threshold

contains

   !> Unweighted sum of a tabulated field; the caller applies its volume element
   !>
   !> @param[in] f Tabulated field values
   !> @param[in] nthreads Worker budget
   function integrate_values(f, nthreads) result(total)
      !> Tabulated field values
      real(wp), intent(in) :: f(:)
      !> Worker budget
      integer, intent(in) :: nthreads
      !> Unweighted sum
      real(wp) :: total

      total = integrate_blocks(f, nthreads)
   end function integrate_values

   !> Weighted sum from final volume weights, without a whole-field temporary
   !>
   !> @param[in] f Tabulated field values
   !> @param[in] nthreads Worker budget
   !> @param[in] weights Final volume weights, bohr**3, same length as f
   function integrate_weighted_values(f, nthreads, weights) result(total)
      !> Tabulated field values
      real(wp), intent(in) :: f(:)
      !> Worker budget
      integer, intent(in) :: nthreads
      !> Final volume weights, bohr**3
      real(wp), intent(in) :: weights(:)
      !> Weighted integral
      real(wp) :: total

      total = integrate_blocks(f, nthreads, weights)
   end function integrate_weighted_values

   !> Partition, evaluate and merge; weighted callers supply final volume weights
   !>
   !> @param[in] f Tabulated field values
   !> @param[in] nthreads Worker budget
   !> @param[in] weights Optional volume weights, bohr**3, same length as f
   function integrate_blocks(f, nthreads, weights) result(total)
      !> Tabulated field values
      real(wp), intent(in) :: f(:)
      !> Worker budget
      integer, intent(in) :: nthreads
      !> Optional volume weights, bohr**3
      real(wp), optional, intent(in) :: weights(:)
      !> Integral and running correction
      real(wp) :: total, correction
      !> (total, correction) per block, (2, nb)
      real(wp), allocatable :: partial(:, :)
      !> Block count, team size, block index, block bounds
      integer :: nb, nt, b, first, last

      if (size(f) < parallel_threshold) then
         if (present(weights)) then
            call compensated_weighted_sum(f, 1, size(f), weights, total, correction)
         else
            call compensated_sum(f, 1, size(f), total, correction)
         end if
         return
      end if

      nb = 1 + (size(f) - 1)/block_size
      nt = max(1, min(nthreads, nb))
!$    if (omp_in_parallel()) nt = 1
      allocate (partial(2, nb))
      !$omp parallel do schedule(static) num_threads(nt) if(nt > 1) private(first, last)
      do b = 1, nb
         first = 1 + (b - 1)*block_size
         last = first + min(block_size - 1, size(f) - first)
         call compensated_block_sum(f, first, last, partial(1, b), partial(2, b), weights)
      end do
      !$omp end parallel do
      total = merge_compensated_sums(partial)
   end function integrate_blocks

end module moist_math_grid_3d_integration
