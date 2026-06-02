!> Utilities for packed symmetric matrix storage.
!> Uses column-packed upper triangular format (LAPACK style):
!> For a 3x3 matrix, the packed form stores: (1,1), (2,1), (2,2), (3,1), (3,2), (3,3)
!> This allows use with LAPACK packed routines (DSPTRF, DSPTRS, etc.)
module moist_math_packed_sym
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: npair_from_ns
   public :: packed_index
   public :: packed_index_lut
   public :: unpack_sym_matrix
   public :: pack_sym_matrix
   public :: packed_symmatmul
   public :: packed_symmatmul_lut

   !> Lookup table for packed symmetric matrix indices.
   !> Stores precomputed packed_index(i,j) for all i,j in [1,ns].
   !> Thread-safe: each instance has its own copy.
   type :: packed_index_lut
      !> Number of sites (matrix dimension)
      integer :: ns = 0
      !> Precomputed indices: idx(i,j) gives packed_index(i,j)
      integer, allocatable :: idx(:, :)
   contains
      procedure :: init => packed_index_lut_init
      procedure :: dealloc => packed_index_lut_dealloc
   end type packed_index_lut

contains

   !> Initialize the lookup table for a given matrix dimension.
   !> @param[inout] self  Lookup table to initialize
   !> @param[in]    ns    Matrix dimension (number of sites)
   subroutine packed_index_lut_init(self, ns)
      class(packed_index_lut), intent(inout) :: self
      integer, intent(in) :: ns
      integer :: i, j

      call self%dealloc()
      self%ns = ns
      allocate (self%idx(ns, ns))

      do j = 1, ns
         do i = 1, ns
            self%idx(i, j) = packed_index(i, j)
         end do
      end do
   end subroutine packed_index_lut_init

   !> Deallocate the lookup table.
   !> @param[inout] self  Lookup table to deallocate
   subroutine packed_index_lut_dealloc(self)
      class(packed_index_lut), intent(inout) :: self
      if (allocated(self%idx)) deallocate (self%idx)
      self%ns = 0
   end subroutine packed_index_lut_dealloc

   !> Return number of unique pairs for ns sites (upper triangle including diagonal).
   !> npair = ns * (ns + 1) / 2
   pure function npair_from_ns(ns) result(npair)
      integer, intent(in) :: ns
      integer :: npair
      npair = ns*(ns + 1)/2
   end function npair_from_ns

   !> Convert (i, j) indices to packed column-major upper-triangular index.
   !> Handles i > j by swapping, so packed_index(i,j) == packed_index(j,i).
   !> @param[in] i  Row index (1 to ns)
   !> @param[in] j  Column index (1 to ns)
   !> @return Index into packed array (1 to npair)
   pure function packed_index(i, j) result(idx)
      integer, intent(in) :: i, j
      integer :: idx
      integer :: ii, jj
      if (i <= j) then
         ii = i
         jj = j
      else
         ii = j
         jj = i
      end if
      idx = jj*(jj - 1)/2 + ii
   end function packed_index

   !> Unpack a symmetric matrix from packed to full storage.
   !> @param[in]  packed  Input packed array (npair)
   !> @param[out] full    Output full matrix (ns, ns), symmetric
   !> @param[in]  ns      Matrix dimension
   pure subroutine unpack_sym_matrix(packed, full, ns)
      real(wp), intent(in) :: packed(:)
      real(wp), intent(out) :: full(:, :)
      integer, intent(in) :: ns
      integer :: i, j, idx
      do j = 1, ns
         do i = 1, ns
            idx = packed_index(i, j)
            full(i, j) = packed(idx)
         end do
      end do
   end subroutine unpack_sym_matrix

   !> Pack a symmetric matrix from full to packed storage.
   !> Only the upper triangle of full is read.
   !> @param[in]  full    Input full matrix (ns, ns)
   !> @param[out] packed  Output packed array (npair)
   !> @param[in]  ns      Matrix dimension
   pure subroutine pack_sym_matrix(full, packed, ns)
      real(wp), intent(in) :: full(:, :)
      real(wp), intent(out) :: packed(:)
      integer, intent(in) :: ns
      integer :: i, j, idx
      do j = 1, ns
         do i = 1, j
            idx = packed_index(i, j)
            packed(idx) = full(i, j)
         end do
      end do
   end subroutine pack_sym_matrix

   !> Symmetric matrix multiplication for packed matrices: C = A * B.
   !> Given two symmetric matrices A and B in packed form, computes C = A * B
   !> and stores only the upper triangle of C in packed form.
   !> Note: The product of two symmetric matrices is not generally symmetric,
   !> but for the RISM equation the specific matrix products preserve symmetry.
   !> @param[in]  A_packed  First symmetric matrix in packed form (npair)
   !> @param[in]  B_packed  Second symmetric matrix in packed form (npair)
   !> @param[out] C_packed  Result packed matrix (npair), upper triangle only
   !> @param[in]  ns        Matrix dimension
   pure subroutine packed_symmatmul(A_packed, B_packed, C_packed, ns)
      real(wp), intent(in) :: A_packed(:)
      real(wp), intent(in) :: B_packed(:)
      real(wp), intent(out) :: C_packed(:)
      integer, intent(in) :: ns
      integer :: i, j, k, idx_c
      real(wp) :: a_ik, b_kj, sum_val

      do j = 1, ns
         do i = 1, j
            idx_c = packed_index(i, j)
            sum_val = 0.0_wp
            do k = 1, ns
               a_ik = A_packed(packed_index(i, k))
               b_kj = B_packed(packed_index(k, j))
               sum_val = sum_val + a_ik*b_kj
            end do
            C_packed(idx_c) = sum_val
         end do
      end do
   end subroutine packed_symmatmul

   !> Symmetric matrix multiplication for packed matrices using LUT: C = A * B.
   !> Uses precomputed lookup table for vectorizable index access.
   !> @param[in]  A_packed  First symmetric matrix in packed form (npair)
   !> @param[in]  B_packed  Second symmetric matrix in packed form (npair)
   !> @param[out] C_packed  Result packed matrix (npair), upper triangle only
   !> @param[in]  lut       Precomputed index lookup table
   pure subroutine packed_symmatmul_lut(A_packed, B_packed, C_packed, lut)
      real(wp), intent(in) :: A_packed(:)
      real(wp), intent(in) :: B_packed(:)
      real(wp), intent(out) :: C_packed(:)
      type(packed_index_lut), intent(in) :: lut
      integer :: i, j, k, idx_c, ns
      real(wp) :: sum_val

      ns = lut%ns
      do j = 1, ns
         do i = 1, j
            idx_c = lut%idx(i, j)
            sum_val = 0.0_wp
            do k = 1, ns
               sum_val = sum_val + A_packed(lut%idx(i, k))*B_packed(lut%idx(k, j))
            end do
            C_packed(idx_c) = sum_val
         end do
      end do
   end subroutine packed_symmatmul_lut

end module moist_math_packed_sym
