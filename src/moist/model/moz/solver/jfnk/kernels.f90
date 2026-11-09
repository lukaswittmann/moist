!> Blocked vector kernels of the JFNK solve
module moist_model_moz_solver_jfnk_kernels
   use mctc_env, only: wp
   use moist_math_blas, only: gemv, gemv_rows, dot
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none(type, external)
   private

   public :: parallel_vector_size, reduction_block_size
   public :: blocked_dot, project_basis, update_basis, field_square_sum, &
      & copy_vector, divide_vector, unpack_vector, trial_field
   public :: generate_givens, apply_givens, upper_backsolve

   !> Avoid OpenMP team overhead for small vectors
   integer, parameter :: parallel_vector_size = 65536
   !> Fixed reduction partitions, independent of the OpenMP team size
   integer, parameter :: reduction_block_size = 8192

contains

   !> Dot product with fixed block and combination order across thread counts
   !>
   !> - Small vectors retain the intrinsic path
   !>
   !> @param[in] x   first vector
   !> @param[in] y   second vector, of the same length
   !> @param[in] nt  OpenMP team size
   function blocked_dot(x, y, nt) result(value)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> First vector
      real(wp), intent(in) :: x(:)
      !> Second vector, of the same length
      real(wp), intent(in) :: y(:)
      !> Scalar dot product
      real(wp) :: value
      !> One partial sum per fixed block
      real(wp), allocatable :: partial(:)
      !> Block index, block count, and first/last element
      integer :: b, nb, first, last

      if (size(x) < parallel_vector_size) then
         value = dot_product(x, y)
         return
      end if
      nb = (size(x) - 1)/reduction_block_size + 1
      allocate (partial(nb))
      !$omp parallel do schedule(static) private(first, last) num_threads(nt)
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         last = min(size(x), first + reduction_block_size - 1)
         partial(b) = dot(x(first:last), y(first:last))
      end do
      !$omp end parallel do
      value = sum(partial)
   end function blocked_dot

   !> Parallel BLAS projection without packing strided matrix tiles
   !>
   !> @param[in]  n   rows of the basis
   !> @param[in]  m   columns of the basis
   !> @param[in]  v   basis (n, m)
   !> @param[in]  w   vector to project (n)
   !> @param[out] h   projection coefficients (m), combined in fixed block order
   !> @param[in]  nt  OpenMP team size
   subroutine project_basis(n, m, v, w, h, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Rows of the basis, explicit to preserve the original BLAS leading dimension
      integer, intent(in) :: n
      !> Columns of the basis
      integer, intent(in) :: m
      !> Contiguous basis (n, m)
      real(wp), intent(in) :: v(n, m)
      !> Contiguous vector to project (n)
      real(wp), intent(in) :: w(n)
      !> Projection coefficients
      real(wp), intent(out) :: h(m)
      !> Small block products, not field-sized packing buffers
      real(wp), allocatable :: partial(:, :)
      !> Block bounds and coefficient index
      integer :: nb, b, first, count, j
      if (n < parallel_vector_size) then
         call gemv(v, w, h, trans="t")
         return
      end if
      nb = (n - 1)/reduction_block_size + 1
      allocate (partial(m, nb))
      !$omp parallel do schedule(static) private(first, count) num_threads(nt)
      do b = 1, nb
         first = (b - 1)*reduction_block_size + 1
         count = min(reduction_block_size, n - first + 1)
         call gemv_rows(v, first, count, w(first:first + count - 1), &
            & partial(:, b), 1.0_wp, 0.0_wp, "T")
      end do
      !$omp end parallel do
      do j = 1, m
         h(j) = sum(partial(j, :))
      end do
   end subroutine project_basis

   !> Cache-blocked BLAS basis update, sharing the destination across columns in cache
   !>
   !> @param[in]    n      rows of the basis
   !> @param[in]    m      columns of the basis
   !> @param[in]    v      basis (n, m)
   !> @param[in]    h      coefficients (m)
   !> @param[in]    alpha  scale factor
   !> @param[inout] w      vector (n), updated to w + alpha*V*h
   !> @param[in]    nt     OpenMP team size
   subroutine update_basis(n, m, v, h, alpha, w, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Rows of the basis, explicit to preserve the BLAS leading dimension without packing
      integer, intent(in) :: n
      !> Columns of the basis
      integer, intent(in) :: m
      !> Basis (n, m)
      real(wp), intent(in) :: v(n, m)
      !> Coefficients (m)
      real(wp), intent(in) :: h(m)
      !> Scale factor
      real(wp), intent(in) :: alpha
      !> Output vector
      real(wp), intent(inout) :: w(n)
      !> Block index and bounds
      integer :: b, first, count, j
      if (n < parallel_vector_size) then
         do j = 1, m
            w = w + alpha*h(j)*v(:, j)
         end do
         return
      end if
      !$omp parallel do schedule(static) private(first, count) num_threads(nt)
      do b = 1, (n - 1)/reduction_block_size + 1
         first = (b - 1)*reduction_block_size + 1
         count = min(reduction_block_size, n - first + 1)
         call gemv_rows(v, first, count, h, w(first:first + count - 1), alpha, 1.0_wp, "N")
      end do
      !$omp end parallel do
   end subroutine update_basis

   !> Squared field norm without flattening or allocating a field-sized temporary
   !>
   !> - Sum of squares with fixed block and combination order
   !>
   !> @param[in] x   input field in native Fortran storage order
   !> @param[in] nt  OpenMP team size
   function field_square_sum(x, nt) result(value)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input field
      real(wp), intent(in) :: x(:, :)
      !> Sum of squares
      real(wp) :: value
      !> Partial sums in spatial-block then column order
      real(wp), allocatable :: partial(:, :)
      !> Block/column indices, block count, and element bounds
      integer :: b, j, nb, first, last

      if (size(x) < parallel_vector_size) then
         value = sum(x**2)
         return
      end if
      nb = (size(x, 1) - 1)/reduction_block_size + 1
      allocate (partial(nb, size(x, 2)))
      !$omp parallel do collapse(2) schedule(static) private(first, last) num_threads(nt)
      do j = 1, size(x, 2)
         do b = 1, nb
            first = (b - 1)*reduction_block_size + 1
            last = min(size(x, 1), first + reduction_block_size - 1)
            partial(b, j) = sum(x(first:last, j)**2)
         end do
      end do
      !$omp end parallel do
      value = sum(partial)
   end function field_square_sum

   !> Copy a flat vector with the same static ownership as Krylov updates
   !>
   !> @param[in] source Input vector
   !> @param[out] dest Copied vector
   !> @param[in] nt OpenMP team size
   subroutine copy_vector(source, dest, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input vector
      real(wp), intent(in) :: source(:)
      !> Output vector of the same length
      real(wp), intent(out) :: dest(:)
      !> Element index
      integer :: i
      !$omp parallel do schedule(static) num_threads(merge(nt, 1, size(source) >= parallel_vector_size))
      do i = 1, size(source)
         dest(i) = source(i)
      end do
      !$omp end parallel do
   end subroutine copy_vector

   !> Normalize a vector without changing division to reciprocal multiplication
   !>
   !> @param[in] source Input vector
   !> @param[in] divisor Normalization divisor
   !> @param[out] dest Normalized vector
   !> @param[in] nt OpenMP team size
   subroutine divide_vector(source, divisor, dest, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Input vector
      real(wp), intent(in) :: source(:)
      !> Normalization divisor
      real(wp), intent(in) :: divisor
      !> Output vector of the same length
      real(wp), intent(out) :: dest(:)
      !> Element index
      integer :: i
      !$omp parallel do schedule(static) num_threads(merge(nt, 1, size(source) >= parallel_vector_size))
      do i = 1, size(source)
         dest(i) = source(i)/divisor
      end do
      !$omp end parallel do
   end subroutine divide_vector

   !> Unpack a vector in Fortran storage order without a reshape temporary
   !>
   !> @param[in] source Packed vector
   !> @param[out] dest Unpacked field
   !> @param[in] nt OpenMP team size
   subroutine unpack_vector(source, dest, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Packed input vector
      real(wp), intent(in) :: source(:)
      !> Output field with size equal to the vector length
      real(wp), intent(out) :: dest(:, :)
      !> Spatial and pair indices
      integer :: i, j
      !$omp parallel do collapse(2) schedule(static) num_threads(merge(nt, 1, size(source) >= parallel_vector_size))
      do j = 1, size(dest, 2)
         do i = 1, size(dest, 1)
            dest(i, j) = source(i + (j - 1)*size(dest, 1))
         end do
      end do
      !$omp end parallel do
   end subroutine unpack_vector

   !> Form a trial iterate with independent elementwise updates
   !>
   !> @param[in] gamma Current field
   !> @param[in] delta Newton direction
   !> @param[in] alpha Step length
   !> @param[out] trial Trial field
   !> @param[in] nt OpenMP team size
   subroutine trial_field(gamma, delta, alpha, trial, nt)
      !> OpenMP team size from the solver's thread count
      integer, intent(in) :: nt
      !> Current field
      real(wp), intent(in) :: gamma(:, :)
      !> Newton direction
      real(wp), intent(in) :: delta(:, :)
      !> Step length
      real(wp), intent(in) :: alpha
      !> Trial field of matching shape
      real(wp), intent(out) :: trial(:, :)
      !> Spatial and pair indices
      integer :: i, j
      !$omp parallel do collapse(2) schedule(static) num_threads(merge(nt, 1, size(gamma) >= parallel_vector_size))
      do j = 1, size(gamma, 2)
         do i = 1, size(gamma, 1)
            trial(i, j) = gamma(i, j) - alpha*delta(i, j)
         end do
      end do
      !$omp end parallel do
   end subroutine trial_field

   !> Form a Givens rotation that removes b without an absolute scale cutoff
   !>
   !> @param[in] a First entry
   !> @param[in] b Entry to eliminate
   !> @param[out] c Rotation cosine
   !> @param[out] s Rotation sine
   subroutine generate_givens(a, b, c, s)
      !> First entry
      real(wp), intent(in) :: a
      !> Entry to eliminate
      real(wp), intent(in) :: b
      !> Rotation cosine
      real(wp), intent(out) :: c
      !> Rotation sine
      real(wp), intent(out) :: s
      real(wp) :: r

      if (b == 0.0_wp) then
         c = 1.0_wp
         s = 0.0_wp
      else if (abs(b) > abs(a)) then
         r = -a/b
         s = 1.0_wp/sqrt(1.0_wp + r**2)
         c = s*r
      else
         r = -b/a
         c = 1.0_wp/sqrt(1.0_wp + r**2)
         s = c*r
      end if
   end subroutine generate_givens

   !> Apply a Givens rotation to the pair (a, b)
   !>
   !> - a := c*a - s*b
   !> - b := s*a + c*b, with the old a
   !>
   !> @param[in]    c  rotation cosine
   !> @param[in]    s  rotation sine
   !> @param[inout] a  first component, rotated in place
   !> @param[inout] b  second component, rotated in place
   subroutine apply_givens(c, s, a, b)
      !> Rotation cosine
      real(wp), intent(in) :: c
      !> Rotation sine
      real(wp), intent(in) :: s
      !> First component, rotated in place
      real(wp), intent(inout) :: a
      !> Second component, rotated in place
      real(wp), intent(inout) :: b
      real(wp) :: tmp

      tmp = c*a - s*b
      b = s*a + c*b
      a = tmp
   end subroutine apply_givens

   !> Solve the small triangular Krylov system without dividing by a zero pivot
   !>
   !> - Only the upper triangle of `r` is read
   !> - With `rank_tol`, a pivot at most `rank_tol` times the norm of its column counts as
   !>   zero: that column depends on the earlier ones up to rounding. The test is per
   !>   column, since flexible directions scale each column on its own
   !>
   !> @param[in] r Upper triangular matrix
   !> @param[in] rhs Rotated least-squares right-hand side
   !> @param[out] x Solution; only usable when feasible
   !> @param[out] feasible Nonsingular finite system and finite solution
   !> @param[in] rank_tol Optional relative size below which a pivot counts as zero
   subroutine upper_backsolve(r, rhs, x, feasible, rank_tol)
      !> Upper triangular matrix
      real(wp), intent(in) :: r(:, :)
      !> Rotated right-hand side
      real(wp), intent(in) :: rhs(:)
      !> Solution coefficients
      real(wp), intent(out) :: x(:)
      !> Whether the solution is usable
      logical, intent(out) :: feasible
      !> Relative pivot size that counts as zero
      real(wp), intent(in), optional :: rank_tol

      integer :: n, i, j
      real(wp) :: sum_val

      n = size(rhs)
      feasible = .false.
      if (.not. all(ieee_is_finite(rhs))) return
      do j = 1, n
         if (.not. all(ieee_is_finite(r(1:j, j)))) return
      end do
      do i = n, 1, -1
         if (r(i, i) == 0.0_wp) return
         if (present(rank_tol)) then
            if (abs(r(i, i)) <= rank_tol*norm2(r(1:i, i))) return
         end if
         sum_val = rhs(i)
         do j = i + 1, n
            sum_val = sum_val - r(i, j)*x(j)
         end do
         x(i) = sum_val/r(i, i)
         if (.not. ieee_is_finite(x(i))) return
      end do
      feasible = .true.
   end subroutine upper_backsolve

end module moist_model_moz_solver_jfnk_kernels
