!> Parallel modified Gram-Schmidt with BLAS and fused streaming kernels.
module moist_math_orthogonalization
   use mctc_env, only: wp
   use moist_math_blas, only: dot, axpy
   !$ use omp_lib, only: omp_get_max_threads
   implicit none(type, external)
   private
   public :: moist_orthogonalize_mgs
   !> Keep small vectors serial.
   integer, parameter :: parallel_vector_size = 65536
   !> Deterministic reduction partition, independent of worker count.
   integer, parameter :: reduction_block_size = 8192
   !> Measured crossover: fuse when each field occupies at least 256 MiB.
   !> Below this size the separate BLAS sweeps benefit more from cache reuse.
   integer, parameter :: fused_vector_size = 32*1024*1024
contains
   !> Modified Gram-Schmidt with fused update/next-dot sweeps on one OpenMP team.
   !> @param[in] v Contiguous basis columns.
   !> @param[inout] w Vector to orthogonalize.
   !> @param[out] h Projection coefficients in basis order.
   !> @param[in] fused Optional kernel override for benchmarking and validation.
   !> @param[in] nthreads OpenMP team size; absent or <= 0 takes omp_get_max_threads.
   subroutine moist_orthogonalize_mgs(v, w, h, fused, nthreads)
      !> Existing basis and working vector.
      real(wp), intent(in), contiguous :: v(:, :)
      real(wp), intent(inout), contiguous :: w(:)
      !> Hessenberg column.
      real(wp), intent(out) :: h(:)
      !> Optional explicit kernel selection; default uses a measured size crossover.
      logical, intent(in), optional :: fused
      !> OpenMP team size, e.g. ctx%get_num_threads(); absent or <= 0 takes omp_get_max_threads.
      integer, intent(in), optional :: nthreads
      !> Effective team size.
      integer :: nt
      !> Whether adjacent projection sweeps share a SIMD pass.
      logical :: combine
      !> Fixed-block partial products; order is independent of team size.
      real(wp), allocatable :: partial(:)
      !> Shared scalar projection.
      real(wp) :: coefficient, product
      !> Basis/block indices and block bounds.
      integer :: j, b, first, last, nb, i
      if (size(v, 2) == 0) return
      if (size(w) < parallel_vector_size) then
         do j = 1, size(v, 2)
            h(j) = dot_product(w, v(:, j))
            w = w - h(j)*v(:, j)
         end do
         return
      end if
      nt = 1
      !$ nt = max(1, omp_get_max_threads())
      if (present(nthreads)) then
         if (nthreads > 0) nt = nthreads
      end if
      combine = size(w) >= fused_vector_size
      if (present(fused)) combine = fused
      if (.not. combine) then
         call blas_mgs(v, w, h, nt)
         return
      end if
      nb = (size(w) - 1)/reduction_block_size + 1
      allocate (partial(nb))
      !$omp parallel num_threads(nt) private(j, b, first, last, i, product)
      do j = 1, size(v, 2)
         !$omp do schedule(static)
         do b = 1, nb
            first = (b - 1)*reduction_block_size + 1
            last = min(size(w), first + reduction_block_size - 1)
            if (j > 1) then
               ! Complete one MGS projection and accumulate the next dot
               ! product with a single vector load/store per element.
               product = 0.0_wp
               !$omp simd reduction(+:product)
               do i = first, last
                  w(i) = w(i) - coefficient*v(i, j - 1)
                  product = product + w(i)*v(i, j)
               end do
               !$omp end simd
               partial(b) = product
            else
               partial(b) = dot(w(first:last), v(first:last, j))
            end if
         end do
         !$omp end do
         !$omp single
         coefficient = sum(partial)
         h(j) = coefficient
         !$omp end single
      end do
      !$omp do schedule(static)
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         last = min(size(w), first + reduction_block_size - 1)
         call axpy(-coefficient, v(first:last, size(v, 2)), w(first:last))
      end do
      !$omp end do
      !$omp end parallel
   end subroutine moist_orthogonalize_mgs

   !> Modified Gram-Schmidt using one OpenMP team for every projection/update.
   !> @param[in] v Contiguous basis columns.
   !> @param[inout] w Vector to orthogonalize.
   !> @param[out] h Projection coefficients in basis order.
   !> @param[in] nt OpenMP team size.
   subroutine blas_mgs(v, w, h, nt)
      !> Existing basis and working vector.
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(inout) :: w(:)
      !> Hessenberg column.
      real(wp), intent(out) :: h(:)
      !> OpenMP team size.
      integer, intent(in) :: nt
      !> Fixed-block partial products; order is independent of team size.
      real(wp), allocatable :: partial(:)
      !> Shared scalar projection.
      real(wp) :: coefficient
      !> Basis/block indices and block bounds.
      integer :: j, b, first, last, nb
      if (size(w) < parallel_vector_size) then
         do j = 1, size(v, 2)
            h(j) = dot_product(w, v(:, j))
            w = w - h(j)*v(:, j)
         end do
         return
      end if
      nb = (size(w) - 1)/reduction_block_size + 1
      allocate (partial(nb))
      !$omp parallel num_threads(nt) private(j, b, first, last)
      do j = 1, size(v, 2)
         !$omp do schedule(static)
         do b = 1, nb
            first = (b - 1)*reduction_block_size + 1
            last = min(size(w), first + reduction_block_size - 1)
            partial(b) = dot(w(first:last), v(first:last, j))
         end do
         !$omp end do
         !$omp single
         coefficient = sum(partial)
         h(j) = coefficient
         !$omp end single
         !$omp do schedule(static)
         do b = 1, nb
            first = (b - 1)*reduction_block_size + 1
            last = min(size(w), first + reduction_block_size - 1)
            call axpy(-coefficient, v(first:last, j), w(first:last))
         end do
         !$omp end do
      end do
      !$omp end parallel
   end subroutine blas_mgs

end module moist_math_orthogonalization
