!> Fixed-block volume quadrature with compensated, reproducible accumulation
module moist_math_grid_3d_integration
   use mctc_env, only: wp
   use moist_math_grid_3d_kernel_integration, only: integration_weight, &
      & compensated_sum, compensated_weighted_sum, compensated_block_sum, merge_compensated_sums
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
   integer, parameter :: integration_parallel_threshold = parallel_threshold

contains

   !> Unweighted sum of a tabulated field; the caller applies its volume element
   function integrate_values(f, nthreads) result(total)
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: nthreads
      real(wp) :: total

      total = integrate_blocks(f, nthreads)
   end function integrate_values

   !> Weighted sum, retaining the domain's per-point weight dispatch
   !>
   !> The pure callback and context remain borrowed until all workers finish.
   function integrate_weighted_values(f, nthreads, weight, context) result(total)
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: nthreads
      procedure(integration_weight) :: weight
      class(*), intent(in) :: context
      real(wp) :: total

      total = integrate_blocks(f, nthreads, weight, context)
   end function integrate_weighted_values

   !> Partition, evaluate and merge; weighted callers supply both optional arguments
   function integrate_blocks(f, nthreads, weight, context) result(total)
      real(wp), intent(in) :: f(:)
      integer, intent(in) :: nthreads
      procedure(integration_weight), optional :: weight
      class(*), optional, intent(in) :: context
      real(wp) :: total, correction
      real(wp), allocatable :: partial(:, :)
      integer :: nb, nt, b, first, last

      if (size(f) < parallel_threshold) then
         if (present(weight)) then
            call compensated_weighted_sum(f, 1, size(f), weight, context, total, correction)
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
         call compensated_block_sum(f, first, last, partial(1, b), partial(2, b), weight, context)
      end do
      !$omp end parallel do
      total = merge_compensated_sums(partial)
   end function integrate_blocks

end module moist_math_grid_3d_integration
