!> Utilities for packed symmetric matrix and tensor storage
!>
!> - column-packed upper triangular format (LAPACK style); a 3x3 matrix packs
!>   as (1,1), (2,1), (2,2), (3,1), (3,2), (3,3)
!> - usable with the LAPACK packed routines (DSPTRF, DSPTRS, etc.)
!> - ranks 3 and 4: sorted index tuples, last index slowest; each slab of fixed
!>   last index holds the next lower rank in its own packed order
!> - index arithmetic in int64; counts beyond huge(1) error stop
module moist_math_packed_sym
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: npair_from_ns
   public :: packed_index
   public :: packed_index_lut
   public :: unpack_sym_matrix
   public :: pack_sym_matrix
   public :: pack_sym_part
   public :: packed_symmatmul
   public :: packed_symmatmul_lut
   public :: ntriple_from_ns
   public :: nquad_from_ns
   public :: packed_index3
   public :: packed_index4
   public :: multiplicity3
   public :: multiplicity4
   public :: pack_sym_tensor3
   public :: unpack_sym_tensor3
   public :: pack_sym_tensor4
   public :: unpack_sym_tensor4
   public :: packed_contract2
   public :: packed_contract3
   public :: packed_contract4
   public :: packed_full_contract2
   public :: packed_full_contract3
   public :: packed_full_contract4
   public :: packed_idx33, packed_idx333

   !> `packed_index(i, j)` for 3x3: xx, xy, yy, xz, yz, zz
   integer, parameter :: packed_idx33(3, 3) = reshape([1, 2, 4, 2, 3, 5, 4, 5, 6], [3, 3])

   !> `packed_index3(i, j, k)` for 3x3x3, order of [[packed_index3]]
   integer, parameter :: packed_idx333(3, 3, 3) = reshape([ &
                                                         1, 2, 5, 2, 3, 6, 5, 6, 8, &
                                                         2, 3, 6, 3, 4, 7, 6, 7, 9, &
                                                         5, 6, 8, 6, 7, 9, 8, 9, 10], [3, 3, 3])

   !> Lookup table for packed symmetric matrix indices
   !>
   !> - precomputed packed_index(i,j) for all i,j in [1,ns]
   !> - thread-safe: each instance has its own copy
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

   !> Initialize the lookup table for a given matrix dimension
   !>
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

   !> Deallocate the lookup table
   !>
   !> @param[inout] self  Lookup table to deallocate
   subroutine packed_index_lut_dealloc(self)
      class(packed_index_lut), intent(inout) :: self
      if (allocated(self%idx)) deallocate (self%idx)
      self%ns = 0
   end subroutine packed_index_lut_dealloc

   !> Number of unique pairs for ns sites, upper triangle with diagonal
   !>
   !> - error stop when the count exceeds huge(1)
   !>
   !> @param[in] ns  Dimension, nonnegative
   pure function npair_from_ns(ns) result(npair)
      !> Dimension
      integer, intent(in) :: ns
      !> ns (ns + 1) / 2
      integer :: npair
      integer(int64) :: n

      n = int(ns, int64)*(ns + 1)/2
      if (n > huge(npair)) error stop "npair_from_ns: count exceeds huge(1)"
      npair = int(n)
   end function npair_from_ns

   !> Convert (i, j) indices to packed column-major upper-triangular index
   !>
   !> - i > j is handled by swapping, so packed_index(i,j) == packed_index(j,i)
   !>
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
      idx = int(int(jj, int64)*(jj - 1)/2 + ii)
   end function packed_index

   !> Unpack a symmetric matrix from packed to full storage
   !>
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

   !> Pack a symmetric matrix from full to packed storage
   !>
   !> - only the upper triangle of full is read
   !>
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

   !> Pack the symmetric part of a general matrix to packed storage
   !>
   !> - off-diagonal slot (i, j) holds (full(i,j) + full(j,i))/2
   !> - diagonal copied
   !>
   !> @param[in]  full    Input full matrix (ns, ns), not necessarily symmetric
   !> @param[out] packed  Output packed array (npair)
   !> @param[in]  ns      Matrix dimension
   pure subroutine pack_sym_part(full, packed, ns)
      real(wp), intent(in) :: full(:, :)
      real(wp), intent(out) :: packed(:)
      integer, intent(in) :: ns
      integer :: i, j, idx
      do j = 1, ns
         do i = 1, j - 1
            idx = packed_index(i, j)
            packed(idx) = 0.5_wp*(full(i, j) + full(j, i))
         end do
         packed(packed_index(j, j)) = full(j, j)
      end do
   end subroutine pack_sym_part

   !> Symmetric matrix multiplication for packed matrices: C = A * B
   !>
   !> - A and B in packed form; only the upper triangle of C is stored, packed
   !> - note: the product of two symmetric matrices is not generally
   !>   symmetric, but the RISM equation's matrix products preserve symmetry
   !>
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

   !> Symmetric matrix multiplication for packed matrices using LUT: C = A * B
   !>
   !> - precomputed lookup table for vectorizable index access
   !>
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
               sum_val = sum_val + A_packed(lut%idx(k, i))*B_packed(lut%idx(k, j))
            end do
            C_packed(idx_c) = sum_val
         end do
      end do
   end subroutine packed_symmatmul_lut

   !* ================================================================================= *!
   !*                        Rank 3 and 4: counts and index maps                        *!
   !* ================================================================================= *!

   !> Number of sorted index triples for ns sites
   !>
   !> - error stop when the count exceeds huge(1)
   !>
   !> @param[in] ns  Dimension, nonnegative
   pure function ntriple_from_ns(ns) result(ntriple)
      !> Dimension
      integer, intent(in) :: ns
      !> ns (ns + 1) (ns + 2) / 6
      integer :: ntriple
      integer(int64) :: n

      n = int(ns, int64)*(ns + 1)*(ns + 2)/6
      if (n > huge(ntriple)) error stop "ntriple_from_ns: count exceeds huge(1)"
      ntriple = int(n)
   end function ntriple_from_ns

   !> Number of sorted index quadruples for ns sites
   !>
   !> - error stop when the count exceeds huge(1)
   !>
   !> @param[in] ns  Dimension, nonnegative
   pure function nquad_from_ns(ns) result(nquad)
      !> Dimension
      integer, intent(in) :: ns
      !> ns (ns + 1) (ns + 2) (ns + 3) / 24
      integer :: nquad
      integer(int64) :: n

      n = int(ns, int64)*(ns + 1)*(ns + 2)*(ns + 3)/24
      if (n > huge(nquad)) error stop "nquad_from_ns: count exceeds huge(1)"
      nquad = int(n)
   end function nquad_from_ns

   !> Packed index of a symmetric rank-3 entry
   !>
   !> - arguments in any order; sorted to a <= b <= c
   !> - a + b (b - 1) / 2 + (c - 1) c (c + 1) / 6
   !>
   !> @param[in] i  First index, 1 to ns
   !> @param[in] j  Second index, 1 to ns
   !> @param[in] k  Third index, 1 to ns
   pure function packed_index3(i, j, k) result(idx)
      !> First index
      integer, intent(in) :: i
      !> Second index
      integer, intent(in) :: j
      !> Third index
      integer, intent(in) :: k
      !> Index into the packed array, 1 to ntriple
      integer :: idx
      integer :: a, b, c

      call sort3(i, j, k, a, b, c)
      idx = int(a + int(b, int64)*(b - 1)/2 + int(c - 1, int64)*c*(c + 1)/6)
   end function packed_index3

   !> Packed index of a symmetric rank-4 entry
   !>
   !> - arguments in any order; sorted to a <= b <= c <= d
   !> - rank-3 index of (a, b, c) + (d - 1) d (d + 1) (d + 2) / 24
   !>
   !> @param[in] i  First index, 1 to ns
   !> @param[in] j  Second index, 1 to ns
   !> @param[in] k  Third index, 1 to ns
   !> @param[in] l  Fourth index, 1 to ns
   pure function packed_index4(i, j, k, l) result(idx)
      !> First index
      integer, intent(in) :: i
      !> Second index
      integer, intent(in) :: j
      !> Third index
      integer, intent(in) :: k
      !> Fourth index
      integer, intent(in) :: l
      !> Index into the packed array, 1 to nquad
      integer :: idx
      integer :: a, b, c, d

      call sort4(i, j, k, l, a, b, c, d)
      idx = int(a + int(b, int64)*(b - 1)/2 + int(c - 1, int64)*c*(c + 1)/6 &
         & + int(d - 1, int64)*d*(d + 1)*(d + 2)/24)
   end function packed_index4

   !> Number of distinct permutations of a rank-3 index tuple
   !>
   !> - 1, 3 or 6; dense-sum weight of a packed entry
   !>
   !> @param[in] i  First index
   !> @param[in] j  Second index
   !> @param[in] k  Third index
   pure function multiplicity3(i, j, k) result(m)
      !> First index
      integer, intent(in) :: i
      !> Second index
      integer, intent(in) :: j
      !> Third index
      integer, intent(in) :: k
      !> Number of distinct permutations
      integer :: m
      integer :: a, b, c

      call sort3(i, j, k, a, b, c)
      m = sorted_multiplicity3(a, b, c)
   end function multiplicity3

   !> Number of distinct permutations of a rank-4 index tuple
   !>
   !> - 1, 4, 6, 12 or 24; dense-sum weight of a packed entry
   !>
   !> @param[in] i  First index
   !> @param[in] j  Second index
   !> @param[in] k  Third index
   !> @param[in] l  Fourth index
   pure function multiplicity4(i, j, k, l) result(m)
      !> First index
      integer, intent(in) :: i
      !> Second index
      integer, intent(in) :: j
      !> Third index
      integer, intent(in) :: k
      !> Fourth index
      integer, intent(in) :: l
      !> Number of distinct permutations
      integer :: m
      integer :: a, b, c, d

      call sort4(i, j, k, l, a, b, c, d)
      m = sorted_multiplicity4(a, b, c, d)
   end function multiplicity4

   !> Number of distinct permutations of a sorted rank-3 index tuple
   !>
   !> - requires a <= b <= c; no sorting
   !>
   !> @param[in] a  Smallest index
   !> @param[in] b  Middle index
   !> @param[in] c  Largest index
   pure function sorted_multiplicity3(a, b, c) result(m)
      !> Smallest index
      integer, intent(in) :: a
      !> Middle index
      integer, intent(in) :: b
      !> Largest index
      integer, intent(in) :: c
      !> Number of distinct permutations
      integer :: m

      if (a == c) then
         m = 1
      else if (a == b .or. b == c) then
         m = 3
      else
         m = 6
      end if
   end function sorted_multiplicity3

   !> Number of distinct permutations of a sorted rank-4 index tuple
   !>
   !> - requires a <= b <= c <= d; no sorting
   !> - 4! over the product of run-length factorials
   !>
   !> @param[in] a  Smallest index
   !> @param[in] b  Second smallest index
   !> @param[in] c  Second largest index
   !> @param[in] d  Largest index
   pure function sorted_multiplicity4(a, b, c, d) result(m)
      !> Smallest index
      integer, intent(in) :: a
      !> Second smallest index
      integer, intent(in) :: b
      !> Second largest index
      integer, intent(in) :: c
      !> Largest index
      integer, intent(in) :: d
      !> Number of distinct permutations
      integer :: m
      integer :: run, denom

      denom = 1
      run = 1
      if (b == a) then
         run = 2
         denom = 2
      end if
      if (c == b) then
         run = run + 1
         denom = denom*run
      else
         run = 1
      end if
      if (d == c) then
         run = run + 1
         denom = denom*run
      end if
      m = 24/denom
   end function sorted_multiplicity4

   !> Sort three integers ascending
   !>
   !> @param[in]  i  First value
   !> @param[in]  j  Second value
   !> @param[in]  k  Third value
   !> @param[out] a  Smallest
   !> @param[out] b  Middle
   !> @param[out] c  Largest
   pure subroutine sort3(i, j, k, a, b, c)
      !> First value
      integer, intent(in) :: i
      !> Second value
      integer, intent(in) :: j
      !> Third value
      integer, intent(in) :: k
      !> Smallest
      integer, intent(out) :: a
      !> Middle
      integer, intent(out) :: b
      !> Largest
      integer, intent(out) :: c

      a = i
      b = j
      c = k
      call swap_ordered(a, b)
      call swap_ordered(b, c)
      call swap_ordered(a, b)
   end subroutine sort3

   !> Sort four integers ascending
   !>
   !> @param[in]  i  First value
   !> @param[in]  j  Second value
   !> @param[in]  k  Third value
   !> @param[in]  l  Fourth value
   !> @param[out] a  Smallest
   !> @param[out] b  Second smallest
   !> @param[out] c  Second largest
   !> @param[out] d  Largest
   pure subroutine sort4(i, j, k, l, a, b, c, d)
      !> First value
      integer, intent(in) :: i
      !> Second value
      integer, intent(in) :: j
      !> Third value
      integer, intent(in) :: k
      !> Fourth value
      integer, intent(in) :: l
      !> Smallest
      integer, intent(out) :: a
      !> Second smallest
      integer, intent(out) :: b
      !> Second largest
      integer, intent(out) :: c
      !> Largest
      integer, intent(out) :: d

      a = i
      b = j
      c = k
      d = l
      ! Five-comparator network
      call swap_ordered(a, b)
      call swap_ordered(c, d)
      call swap_ordered(a, c)
      call swap_ordered(b, d)
      call swap_ordered(b, c)
   end subroutine sort4

   !> Order two integers so that x <= y
   !>
   !> @param[in,out] x  Becomes the smaller value
   !> @param[in,out] y  Becomes the larger value
   pure subroutine swap_ordered(x, y)
      !> First value
      integer, intent(inout) :: x
      !> Second value
      integer, intent(inout) :: y
      integer :: t

      if (x > y) then
         t = x
         x = y
         y = t
      end if
   end subroutine swap_ordered

   !* ================================================================================= *!
   !*                           Rank 3 and 4: pack and unpack                           *!
   !* ================================================================================= *!

   !> Pack a symmetric rank-3 tensor from full to packed storage
   !>
   !> - reads only entries with i <= j <= k; no averaging
   !> - sorted loop order is packed order; running counter, no index function
   !>
   !> @param[in]  full    Full tensor (ns, ns, ns)
   !> @param[out] packed  Packed tensor (ntriple)
   !> @param[in]  ns      Dimension
   pure subroutine pack_sym_tensor3(full, packed, ns)
      !> Full tensor
      real(wp), intent(in) :: full(:, :, :)
      !> Packed tensor
      real(wp), intent(out) :: packed(:)
      !> Dimension
      integer, intent(in) :: ns
      integer :: i, j, k, p

      p = 0
      do k = 1, ns
         do j = 1, k
            do i = 1, j
               p = p + 1
               packed(p) = full(i, j, k)
            end do
         end do
      end do
   end subroutine pack_sym_tensor3

   !> Unpack a symmetric rank-3 tensor from packed to full storage
   !>
   !> - every permutation of each sorted tuple is filled
   !> - packed-order loop writes all 6 orders; repeated indices rewrite one slot
   !>
   !> @param[in]  packed  Packed tensor (ntriple)
   !> @param[out] full    Full tensor (ns, ns, ns), symmetric
   !> @param[in]  ns      Dimension
   pure subroutine unpack_sym_tensor3(packed, full, ns)
      !> Packed tensor
      real(wp), intent(in) :: packed(:)
      !> Full tensor
      real(wp), intent(out) :: full(:, :, :)
      !> Dimension
      integer, intent(in) :: ns
      integer :: a, b, c, p
      real(wp) :: x

      p = 0
      do c = 1, ns
         do b = 1, c
            do a = 1, b
               p = p + 1
               x = packed(p)
               full(a, b, c) = x
               full(a, c, b) = x
               full(b, a, c) = x
               full(b, c, a) = x
               full(c, a, b) = x
               full(c, b, a) = x
            end do
         end do
      end do
   end subroutine unpack_sym_tensor3

   !> Pack a symmetric rank-4 tensor from full to packed storage
   !>
   !> - reads only entries with i <= j <= k <= l; no averaging
   !> - sorted loop order is packed order; running counter, no index function
   !>
   !> @param[in]  full    Full tensor (ns, ns, ns, ns)
   !> @param[out] packed  Packed tensor (nquad)
   !> @param[in]  ns      Dimension
   pure subroutine pack_sym_tensor4(full, packed, ns)
      !> Full tensor
      real(wp), intent(in) :: full(:, :, :, :)
      !> Packed tensor
      real(wp), intent(out) :: packed(:)
      !> Dimension
      integer, intent(in) :: ns
      integer :: i, j, k, l, p

      p = 0
      do l = 1, ns
         do k = 1, l
            do j = 1, k
               do i = 1, j
                  p = p + 1
                  packed(p) = full(i, j, k, l)
               end do
            end do
         end do
      end do
   end subroutine pack_sym_tensor4

   !> Unpack a symmetric rank-4 tensor from packed to full storage
   !>
   !> - every permutation of each sorted tuple is filled
   !> - packed-order loop writes all 24 orders; repeated indices rewrite one slot
   !>
   !> @param[in]  packed  Packed tensor (nquad)
   !> @param[out] full    Full tensor (ns, ns, ns, ns), symmetric
   !> @param[in]  ns      Dimension
   pure subroutine unpack_sym_tensor4(packed, full, ns)
      !> Packed tensor
      real(wp), intent(in) :: packed(:)
      !> Full tensor
      real(wp), intent(out) :: full(:, :, :, :)
      !> Dimension
      integer, intent(in) :: ns
      integer :: a, b, c, d, p
      real(wp) :: x

      p = 0
      do d = 1, ns
         do c = 1, d
            do b = 1, c
               do a = 1, b
                  p = p + 1
                  x = packed(p)
                  full(a, b, c, d) = x
                  full(a, b, d, c) = x
                  full(a, c, b, d) = x
                  full(a, c, d, b) = x
                  full(a, d, b, c) = x
                  full(a, d, c, b) = x
                  full(b, a, c, d) = x
                  full(b, a, d, c) = x
                  full(b, c, a, d) = x
                  full(b, c, d, a) = x
                  full(b, d, a, c) = x
                  full(b, d, c, a) = x
                  full(c, a, b, d) = x
                  full(c, a, d, b) = x
                  full(c, b, a, d) = x
                  full(c, b, d, a) = x
                  full(c, d, a, b) = x
                  full(c, d, b, a) = x
                  full(d, a, b, c) = x
                  full(d, a, c, b) = x
                  full(d, b, a, c) = x
                  full(d, b, c, a) = x
                  full(d, c, a, b) = x
                  full(d, c, b, a) = x
               end do
            end do
         end do
      end do
   end subroutine unpack_sym_tensor4

   !* ================================================================================= *!
   !*                             Rank 2 to 4: contractions                             *!
   !* ================================================================================= *!

   !> Contract one index of a packed symmetric matrix with a vector, H v
   !>
   !> - one sequential pass over the packed upper triangle (DSPMV style)
   !> - off-diagonal entry updates out(i) and out(j); diagonal read once
   !>
   !> @param[in]  h_packed  Packed matrix (npair)
   !> @param[in]  v         Vector (ns)
   !> @param[out] out       H v (ns)
   !> @param[in]  ns        Dimension
   pure subroutine packed_contract2(h_packed, v, out, ns)
      !> Packed matrix
      real(wp), intent(in) :: h_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Contracted vector
      real(wp), intent(out) :: out(:)
      !> Dimension
      integer, intent(in) :: ns
      integer :: i, j, p
      real(wp) :: acc, hij, vj

      ! out(j) first written in column j; later columns add the k > j terms
      p = 0
      do j = 1, ns
         vj = v(j)
         acc = 0.0_wp
         do i = 1, j - 1
            p = p + 1
            hij = h_packed(p)
            out(i) = out(i) + hij*vj
            acc = acc + hij*v(i)
         end do
         p = p + 1
         out(j) = acc + h_packed(p)*vj
      end do
   end subroutine packed_contract2

   !> Contract one index of a packed symmetric rank-3 tensor with a vector
   !>
   !> - out(i, j) = sum_k T(i, j, k) v(k), symmetric, rank-2 packed
   !> - k split at i and j so each range has a fixed sorted tuple; no sorting
   !>
   !> @param[in]  t_packed    Packed tensor (ntriple)
   !> @param[in]  v           Vector (ns)
   !> @param[out] out_packed  Contracted matrix, packed (npair)
   !> @param[in]  ns          Dimension
   pure subroutine packed_contract3(t_packed, v, out_packed, ns)
      !> Packed tensor
      real(wp), intent(in) :: t_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Contracted matrix, packed
      real(wp), intent(out) :: out_packed(:)
      !> Dimension
      integer, intent(in) :: ns
      integer :: i, j, k, p
      integer(int64) :: tet_j, base
      real(wp) :: acc

      p = 0
      do j = 1, ns
         tet_j = int(j - 1, int64)*j*(j + 1)/6
         do i = 1, j
            acc = 0.0_wp
            ! k <= i: sorted (k, i, j)
            base = int(i, int64)*(i - 1)/2 + tet_j
            do k = 1, i
               acc = acc + t_packed(base + k)*v(k)
            end do
            ! i < k <= j: sorted (i, k, j)
            base = i + tet_j
            do k = i + 1, j
               acc = acc + t_packed(base + int(k, int64)*(k - 1)/2)*v(k)
            end do
            ! k > j: sorted (i, j, k)
            base = i + int(j, int64)*(j - 1)/2
            do k = j + 1, ns
               acc = acc + t_packed(base + int(k - 1, int64)*k*(k + 1)/6)*v(k)
            end do
            p = p + 1
            out_packed(p) = acc
         end do
      end do
   end subroutine packed_contract3

   !> Contract one index of a packed symmetric rank-4 tensor with a vector
   !>
   !> - out(i, j, k) = sum_l F(i, j, k, l) v(l), symmetric, rank-3 packed
   !> - l split at i, j and k so each range has a fixed sorted tuple; no sorting
   !>
   !> @param[in]  f_packed    Packed tensor (nquad)
   !> @param[in]  v           Vector (ns)
   !> @param[out] out_packed  Contracted tensor, packed (ntriple)
   !> @param[in]  ns          Dimension
   pure subroutine packed_contract4(f_packed, v, out_packed, ns)
      !> Packed tensor
      real(wp), intent(in) :: f_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Contracted tensor, packed
      real(wp), intent(out) :: out_packed(:)
      !> Dimension
      integer, intent(in) :: ns
      integer :: i, j, k, l, p
      integer(int64) :: tet_j, tet_k, pent_k, base
      real(wp) :: acc

      p = 0
      do k = 1, ns
         tet_k = int(k - 1, int64)*k*(k + 1)/6
         pent_k = int(k - 1, int64)*k*(k + 1)*(k + 2)/24
         do j = 1, k
            tet_j = int(j - 1, int64)*j*(j + 1)/6
            do i = 1, j
               acc = 0.0_wp
               ! l <= i: sorted (l, i, j, k)
               base = int(i, int64)*(i - 1)/2 + tet_j + pent_k
               do l = 1, i
                  acc = acc + f_packed(base + l)*v(l)
               end do
               ! i < l <= j: sorted (i, l, j, k)
               base = i + tet_j + pent_k
               do l = i + 1, j
                  acc = acc + f_packed(base + int(l, int64)*(l - 1)/2)*v(l)
               end do
               ! j < l <= k: sorted (i, j, l, k)
               base = i + int(j, int64)*(j - 1)/2 + pent_k
               do l = j + 1, k
                  acc = acc + f_packed(base + int(l - 1, int64)*l*(l + 1)/6)*v(l)
               end do
               ! l > k: sorted (i, j, k, l)
               base = i + int(j, int64)*(j - 1)/2 + tet_k
               do l = k + 1, ns
                  acc = acc + f_packed(base + int(l - 1, int64)*l*(l + 1)*(l + 2)/24)*v(l)
               end do
               p = p + 1
               out_packed(p) = acc
            end do
         end do
      end do
   end subroutine packed_contract4

   !> Full contraction of a packed symmetric matrix, H:vv
   !>
   !> - multiplicity-weighted sum over packed entries, i <= j
   !>
   !> @param[in] h_packed  Packed matrix (npair)
   !> @param[in] v         Vector (ns)
   !> @param[in] ns        Dimension
   pure function packed_full_contract2(h_packed, v, ns) result(r)
      !> Packed matrix
      real(wp), intent(in) :: h_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Dimension
      integer, intent(in) :: ns
      !> sum_ij H(i, j) v(i) v(j)
      real(wp) :: r
      integer :: i, j, p
      real(wp) :: vj

      r = 0.0_wp
      p = 0
      do j = 1, ns
         vj = v(j)
         do i = 1, j - 1
            p = p + 1
            r = r + 2.0_wp*h_packed(p)*v(i)*vj
         end do
         p = p + 1
         r = r + h_packed(p)*vj*vj
      end do
   end function packed_full_contract2

   !> Full contraction of a packed symmetric rank-3 tensor, T:vvv
   !>
   !> - multiplicity-weighted sum over packed entries, i <= j <= k
   !>
   !> @param[in] t_packed  Packed tensor (ntriple)
   !> @param[in] v         Vector (ns)
   !> @param[in] ns        Dimension
   pure function packed_full_contract3(t_packed, v, ns) result(r)
      !> Packed tensor
      real(wp), intent(in) :: t_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Dimension
      integer, intent(in) :: ns
      !> sum_ijk T(i, j, k) v(i) v(j) v(k)
      real(wp) :: r
      integer :: i, j, k, p
      real(wp) :: vk, vjk

      r = 0.0_wp
      p = 0
      do k = 1, ns
         vk = v(k)
         do j = 1, k
            vjk = v(j)*vk
            do i = 1, j
               p = p + 1
               r = r + real(sorted_multiplicity3(i, j, k), wp)*t_packed(p)*v(i)*vjk
            end do
         end do
      end do
   end function packed_full_contract3

   !> Full contraction of a packed symmetric rank-4 tensor, F:vvvv
   !>
   !> - multiplicity-weighted sum over packed entries, i <= j <= k <= l
   !>
   !> @param[in] f_packed  Packed tensor (nquad)
   !> @param[in] v         Vector (ns)
   !> @param[in] ns        Dimension
   pure function packed_full_contract4(f_packed, v, ns) result(r)
      !> Packed tensor
      real(wp), intent(in) :: f_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Dimension
      integer, intent(in) :: ns
      !> sum_ijkl F(i, j, k, l) v(i) v(j) v(k) v(l)
      real(wp) :: r
      integer :: i, j, k, l, p
      real(wp) :: vl, vkl, vjkl

      r = 0.0_wp
      p = 0
      do l = 1, ns
         vl = v(l)
         do k = 1, l
            vkl = v(k)*vl
            do j = 1, k
               vjkl = v(j)*vkl
               do i = 1, j
                  p = p + 1
                  r = r + real(sorted_multiplicity4(i, j, k, l), wp)*f_packed(p)*v(i)*vjkl
               end do
            end do
         end do
      end do
   end function packed_full_contract4

end module moist_math_packed_sym
