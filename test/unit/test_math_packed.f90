!> Test suite for packed symmetric matrix utilities in moist_math_packed
module test_math_packed
   use mctc_env, only: wp
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_packed, only: npair_from_ns, packed_index, packed_index_lut, unpack_sym_matrix, &
      & pack_sym_matrix, pack_sym_part, packed_symmatmul, packed_symmatmul_lut, ntriple_from_ns, nquad_from_ns, &
      & packed_index3, packed_index4, multiplicity3, multiplicity4, pack_sym_tensor3, unpack_sym_tensor3, &
      & pack_sym_tensor4, unpack_sym_tensor4, packed_contract2, packed_contract3, packed_contract4, &
      & packed_full_contract2, packed_full_contract3, packed_full_contract4
   use moist_math_lapack_sptrs, only: wrap_dsptrf, wrap_dsptrs
   implicit none(type, external)
   private

   public :: collect_math_packed

   !> Round-off level for sums of up to 256 products of O(1) numbers (dense F:vvvv, ns = 4)
   real(wp), parameter :: thr = 1.0e-13_wp
   integer, parameter :: sizes(4) = [1, 2, 3, 7]
   !> Dimensions for the rank-3 and rank-4 tests
   integer, parameter :: tsizes(4) = [1, 2, 3, 4]

contains

   !> Collect all packed symmetric tests
   subroutine collect_math_packed(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("npair", test_npair), &
                  new_unittest("index_layout_literal", test_index_layout_literal), &
                  new_unittest("index_vs_enumeration", test_index_vs_enumeration), &
                  new_unittest("lut_vs_index", test_lut_vs_index), &
                  new_unittest("lut_reinit_dealloc", test_lut_reinit_dealloc), &
                  new_unittest("pack_unpack_roundtrip", test_pack_unpack_roundtrip), &
                  new_unittest("pack_reads_upper_only", test_pack_reads_upper_only), &
                  new_unittest("pack_part_symmetrises", test_pack_part_symmetrises), &
                  new_unittest("symmatmul_vs_dense", test_symmatmul_vs_dense), &
                  new_unittest("symmatmul_symmetric_case", test_symmatmul_symmetric_case), &
                  new_unittest("symmatmul_lut_vs_dense", test_symmatmul_lut_vs_dense), &
                  new_unittest("layout_vs_lapack_sptrs", test_layout_vs_lapack_sptrs), &
                  new_unittest("ntriple_nquad", test_ntriple_nquad), &
                  new_unittest("index3_bijection", test_index3_bijection), &
                  new_unittest("index4_bijection", test_index4_bijection), &
                  new_unittest("index3_slabs_are_rank2", test_index3_slabs_are_rank2), &
                  new_unittest("index4_slabs_are_rank3", test_index4_slabs_are_rank3), &
                  new_unittest("index3_permutation_invariant", test_index3_permutation_invariant), &
                  new_unittest("index4_permutation_invariant", test_index4_permutation_invariant), &
                  new_unittest("multiplicity_counts_dense", test_multiplicity_counts_dense), &
                  new_unittest("tensor3_roundtrip", test_tensor3_roundtrip), &
                  new_unittest("tensor4_roundtrip", test_tensor4_roundtrip), &
                  new_unittest("tensor_pack_rejects_asymmetric", test_tensor_pack_rejects_asymmetric), &
                  new_unittest("contract2_vs_dense", test_contract2_vs_dense), &
                  new_unittest("contract3_vs_dense", test_contract3_vs_dense), &
                  new_unittest("contract4_vs_dense", test_contract4_vs_dense), &
                  new_unittest("full_contract_vs_dense", test_full_contract_vs_dense), &
                  new_unittest("full_contract_vs_multiplicity", test_full_contract_vs_multiplicity), &
                  new_unittest("contract_sees_perturbation", test_contract_sees_perturbation), &
                  new_unittest("contract2_vs_rowwise", test_contract2_vs_rowwise), &
                  new_unittest("contract2_updates_both", test_contract2_updates_both), &
                  new_unittest("tensor_pack_matches_index", test_tensor_pack_matches_index), &
                  new_unittest("count_boundary", test_count_boundary), &
                  new_unittest("index_boundary", test_index_boundary) &
                  ]
   end subroutine collect_math_packed

   !> Deterministic pseudo-random numbers in [-1, 1) from a linear congruential generator
   subroutine fill_random(x, seed)
      real(wp), intent(out) :: x(:)
      integer, intent(in) :: seed
      integer(int64) :: state
      integer :: i

      state = int(seed, int64)
      do i = 1, size(x)
         state = mod(state*1103515245_int64 + 12345_int64, 2147483648_int64)
         x(i) = 2.0_wp*real(state, wp)/2147483648.0_wp - 1.0_wp
      end do
   end subroutine fill_random

   !> Random packed array of dimension ns
   subroutine random_packed(ns, seed, packed)
      integer, intent(in) :: ns, seed
      real(wp), allocatable, intent(out) :: packed(:)

      allocate (packed(ns*(ns + 1)/2))
      call fill_random(packed, seed)
   end subroutine random_packed

   !> Number of pairs is the triangular number
   subroutine test_npair(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ns

      call check(error, npair_from_ns(0), 0, "npair(0)")
      if (allocated(error)) return
      call check(error, npair_from_ns(1), 1, "npair(1)")
      if (allocated(error)) return
      call check(error, npair_from_ns(3), 6, "npair(3)")
      if (allocated(error)) return
      call check(error, npair_from_ns(7), 28, "npair(7)")
      if (allocated(error)) return
      do ns = 1, 20
         call check(error, npair_from_ns(ns) - npair_from_ns(ns - 1), ns, "npair increment equals ns")
         if (allocated(error)) return
      end do
   end subroutine test_npair

   !> Documented 3x3 order (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), both argument orders
   subroutine test_index_layout_literal(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ri(6) = [1, 2, 2, 3, 3, 3]
      integer, parameter :: rj(6) = [1, 1, 2, 1, 2, 3]
      integer :: k

      do k = 1, 6
         call check(error, packed_index(ri(k), rj(k)), k, "lower-triangle argument order")
         if (allocated(error)) return
         call check(error, packed_index(rj(k), ri(k)), k, "upper-triangle argument order")
         if (allocated(error)) return
      end do
   end subroutine test_index_layout_literal

   !> Index equals the running counter of an explicit column-by-column enumeration
   subroutine test_index_vs_enumeration(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: is, ns, i, j, counter

      do is = 1, size(sizes)
         ns = sizes(is)
         counter = 0
         do j = 1, ns
            do i = 1, j
               counter = counter + 1
               call check(error, packed_index(i, j), counter, "index (i<=j)")
               if (allocated(error)) return
               call check(error, packed_index(j, i), counter, "index symmetric (i>=j)")
               if (allocated(error)) return
            end do
         end do
         call check(error, counter, npair_from_ns(ns), "enumeration count equals npair")
         if (allocated(error)) return
      end do
   end subroutine test_index_vs_enumeration

   !> Lookup table entries equal the index function for every (i,j)
   subroutine test_lut_vs_index(error)
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      integer :: is, ns, i, j

      do is = 1, size(sizes)
         ns = sizes(is)
         call lut%init(ns)
         call check(error, lut%ns, ns, "lut dimension")
         if (allocated(error)) return
         call check(error, all(shape(lut%idx) == [ns, ns]), "lut shape")
         if (allocated(error)) return
         do j = 1, ns
            do i = 1, ns
               call check(error, lut%idx(i, j), packed_index(i, j), "lut entry")
               if (allocated(error)) return
            end do
         end do
      end do
      call lut%dealloc()
   end subroutine test_lut_vs_index

   !> Re-initialization with another size and deallocation reset the table
   subroutine test_lut_reinit_dealloc(error)
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut

      call lut%init(5)
      call lut%init(2)
      call check(error, lut%ns, 2, "re-init dimension")
      if (allocated(error)) return
      call check(error, all(shape(lut%idx) == [2, 2]), "re-init shape")
      if (allocated(error)) return
      call check(error, lut%idx(2, 2), 3, "re-init entry")
      if (allocated(error)) return

      call lut%dealloc()
      call check(error, lut%ns, 0, "dealloc dimension")
      if (allocated(error)) return
      call check(error, .not. allocated(lut%idx), "dealloc table")
      if (allocated(error)) return
      call lut%dealloc()
      call check(error, lut%ns, 0, "dealloc of empty table")
   end subroutine test_lut_reinit_dealloc

   !> Unpack gives an exactly symmetric matrix holding the packed values; pack inverts it
   subroutine test_pack_unpack_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: packed(:), packed2(:), full(:, :)
      integer :: is, ns, i, j, counter

      do is = 1, size(sizes)
         ns = sizes(is)
         call random_packed(ns, 17 + is, packed)
         allocate (full(ns, ns), packed2(size(packed)))
         call unpack_sym_matrix(packed, full, ns)

         counter = 0
         do j = 1, ns
            do i = 1, j
               counter = counter + 1
               call check(error, full(i, j), packed(counter), "upper element equals packed value")
               if (allocated(error)) return
               call check(error, full(j, i), packed(counter), "lower element mirrors packed value")
               if (allocated(error)) return
            end do
         end do

         call pack_sym_matrix(full, packed2, ns)
         call check(error, all(packed2 == packed), "pack(unpack(p)) == p")
         if (allocated(error)) return
         deallocate (full, packed2)
      end do
   end subroutine test_pack_unpack_roundtrip

   !> Only the upper triangle of the full matrix enters pack; the lower one may hold garbage
   subroutine test_pack_reads_upper_only(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: packed(:), packed2(:), full(:, :), dirty(:, :)
      integer :: is, ns, i, j

      do is = 1, size(sizes)
         ns = sizes(is)
         call random_packed(ns, 91 + is, packed)
         allocate (full(ns, ns), dirty(ns, ns), packed2(size(packed)))
         call unpack_sym_matrix(packed, full, ns)
         dirty = full
         do j = 1, ns
            do i = j + 1, ns
               dirty(i, j) = 1.0e30_wp*real(i - j, wp)
            end do
         end do
         call pack_sym_matrix(dirty, packed2, ns)
         call check(error, all(packed2 == packed), "lower triangle ignored by pack")
         if (allocated(error)) return
         deallocate (full, dirty, packed2)
      end do
   end subroutine test_pack_reads_upper_only

   !> pack_sym_part averages the off-diagonal pairs; pack_sym_matrix of the same input does not
   subroutine test_pack_part_symmetrises(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: full(3, 3), packed(6), upper(6), ref(6)

      ! Entries are small integers, so the halves are exact
      full = reshape([1.0_wp, 4.0_wp, 7.0_wp, &
                      2.0_wp, 5.0_wp, 8.0_wp, &
                      3.0_wp, 6.0_wp, 9.0_wp], [3, 3], order=[2, 1])
      ! Slots xx, xy, yy, xz, yz, zz
      ref = [1.0_wp, 3.0_wp, 5.0_wp, 5.0_wp, 7.0_wp, 9.0_wp]

      call pack_sym_part(full, packed, 3)
      call check(error, all(packed == ref), "symmetric part in every slot")
      if (allocated(error)) return

      ! Negative control: the upper-triangle helper keeps full(i,j) only
      call pack_sym_matrix(full, upper, 3)
      call check(error, upper(packed_index(1, 2)) /= packed(packed_index(1, 2)), "upper pack differs at xy")
      if (allocated(error)) return
      call check(error, upper(packed_index(1, 3)) /= packed(packed_index(1, 3)), "upper pack differs at xz")
      if (allocated(error)) return
      call check(error, upper(packed_index(2, 3)) /= packed(packed_index(2, 3)), "upper pack differs at yz")
      if (allocated(error)) return
      call check(error, all(upper([1, 3, 6]) == packed([1, 3, 6])), "diagonals agree")
   end subroutine test_pack_part_symmetrises

   !> C = A*B for independent symmetric A, B: upper triangle equals dense matmul
   subroutine test_symmatmul_vs_dense(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: pa(:), pb(:), pc(:), a(:, :), b(:, :), ref(:, :)
      integer :: is, ns, i, j

      do is = 1, size(sizes)
         ns = sizes(is)
         call random_packed(ns, 5 + is, pa)
         call random_packed(ns, 105 + is, pb)
         allocate (pc(size(pa)), a(ns, ns), b(ns, ns), ref(ns, ns))
         call unpack_sym_matrix(pa, a, ns)
         call unpack_sym_matrix(pb, b, ns)
         ref = matmul(a, b)
         call packed_symmatmul(pa, pb, pc, ns)
         do j = 1, ns
            do i = 1, j
               call check(error, abs(pc(packed_index(i, j)) - ref(i, j)) < thr, "upper triangle of A*B")
               if (allocated(error)) return
            end do
         end do
         deallocate (pc, a, b, ref)
      end do
   end subroutine test_symmatmul_vs_dense

   !> B = A (so A*A is symmetric): packed result unpacks to the full dense square
   subroutine test_symmatmul_symmetric_case(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: pa(:), pc(:), a(:, :), c(:, :), ref(:, :)
      integer :: is, ns

      do is = 1, size(sizes)
         ns = sizes(is)
         call random_packed(ns, 33 + is, pa)
         allocate (pc(size(pa)), a(ns, ns), c(ns, ns), ref(ns, ns))
         call unpack_sym_matrix(pa, a, ns)
         ref = matmul(a, a)
         call packed_symmatmul(pa, pa, pc, ns)
         call unpack_sym_matrix(pc, c, ns)
         call check(error, maxval(abs(c - ref)) < thr, "unpacked A*A equals dense square")
         if (allocated(error)) return
         deallocate (pc, a, c, ref)
      end do
   end subroutine test_symmatmul_symmetric_case

   !> LUT variant agrees with the dense reference and with the plain variant bit for bit
   subroutine test_symmatmul_lut_vs_dense(error)
      type(error_type), allocatable, intent(out) :: error
      type(packed_index_lut) :: lut
      real(wp), allocatable :: pa(:), pb(:), pc(:), pd(:), a(:, :), b(:, :), ref(:, :)
      integer :: is, ns, i, j

      do is = 1, size(sizes)
         ns = sizes(is)
         call random_packed(ns, 7 + is, pa)
         call random_packed(ns, 207 + is, pb)
         allocate (pc(size(pa)), pd(size(pa)), a(ns, ns), b(ns, ns), ref(ns, ns))
         call unpack_sym_matrix(pa, a, ns)
         call unpack_sym_matrix(pb, b, ns)
         ref = matmul(a, b)
         call lut%init(ns)
         call packed_symmatmul_lut(pa, pb, pc, lut)
         do j = 1, ns
            do i = 1, j
               call check(error, abs(pc(packed_index(i, j)) - ref(i, j)) < thr, "LUT upper triangle of A*B")
               if (allocated(error)) return
            end do
         end do
         call packed_symmatmul(pa, pb, pd, ns)
         call check(error, maxval(abs(pc - pd)) < thr, "LUT equals plain variant")
         if (allocated(error)) return
         deallocate (pc, pd, a, b, ref)
      end do
      call lut%dealloc()
   end subroutine test_symmatmul_lut_vs_dense

   !> Packed layout is LAPACK's 'U' layout: DSPTRF/DSPTRS on the packed array solve the unpacked system
   subroutine test_layout_vs_lapack_sptrs(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: pa(:), pf(:), a(:, :), x(:, :), b(:, :)
      integer, allocatable :: ipiv(:)
      integer :: is, ns, info, i, nrhs

      nrhs = 3
      do is = 1, size(sizes)
         ns = sizes(is)
         call random_packed(ns, 61 + is, pa)
         do i = 1, ns
            pa(packed_index(i, i)) = pa(packed_index(i, i)) + real(ns, wp)
         end do
         allocate (a(ns, ns), x(ns, nrhs), b(ns, nrhs), pf(size(pa)), ipiv(ns))
         call unpack_sym_matrix(pa, a, ns)
         call fill_random(b(:, 1), 3)
         call fill_random(b(:, 2), 4)
         call fill_random(b(:, 3), 5)
         x = b
         pf = pa
         call wrap_dsptrf(pf, ipiv, info, "u")
         call check(error, info, 0, "dsptrf info")
         if (allocated(error)) return
         call wrap_dsptrs(pf, x, ipiv, info, "u")
         call check(error, info, 0, "dsptrs info")
         if (allocated(error)) return
         call check(error, maxval(abs(matmul(a, x) - b)) < thr, "unpacked A*x reproduces b")
         if (allocated(error)) return
         deallocate (a, x, b, pf, ipiv)
      end do
   end subroutine test_layout_vs_lapack_sptrs

   !* ================================================================================= *!
   !*                               Rank 3 and 4: helpers                               *!
   !* ================================================================================= *!

   !> Pass if within thr absolutely or 10 thr relative to ref
   pure function within_tol(x, ref) result(ok)
      !> Value under test
      real(wp), intent(in) :: x
      !> Reference
      real(wp), intent(in) :: ref
      !> Within tolerance
      logical :: ok
      ok = abs(x - ref) <= thr .or. abs(x - ref) <= 10.0_wp*thr*abs(ref)
   end function within_tol

   !> Random rank-3 packed tensor of dimension ns
   subroutine random_packed3(ns, seed, packed)
      !> Dimension
      integer, intent(in) :: ns
      !> Generator seed
      integer, intent(in) :: seed
      !> Random packed tensor
      real(wp), allocatable, intent(out) :: packed(:)

      allocate (packed(ntriple_from_ns(ns)))
      call fill_random(packed, seed)
   end subroutine random_packed3

   !> Random rank-4 packed tensor of dimension ns
   subroutine random_packed4(ns, seed, packed)
      !> Dimension
      integer, intent(in) :: ns
      !> Generator seed
      integer, intent(in) :: seed
      !> Random packed tensor
      real(wp), allocatable, intent(out) :: packed(:)

      allocate (packed(nquad_from_ns(ns)))
      call fill_random(packed, seed)
   end subroutine random_packed4

   !* ================================================================================= *!
   !*                        Rank 3 and 4: counts and index maps                        *!
   !* ================================================================================= *!

   !> Tuple counts grow by the next lower rank's count; literal values at ns = 3
   subroutine test_ntriple_nquad(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ns

      call check(error, ntriple_from_ns(0), 0, "ntriple(0)")
      if (allocated(error)) return
      call check(error, nquad_from_ns(0), 0, "nquad(0)")
      if (allocated(error)) return
      call check(error, ntriple_from_ns(3), 10, "ntriple(3)")
      if (allocated(error)) return
      call check(error, nquad_from_ns(3), 15, "nquad(3)")
      if (allocated(error)) return
      do ns = 1, 20
         call check(error, ntriple_from_ns(ns) - ntriple_from_ns(ns - 1), npair_from_ns(ns), &
            & "ntriple increment equals npair")
         if (allocated(error)) return
         call check(error, nquad_from_ns(ns) - nquad_from_ns(ns - 1), ntriple_from_ns(ns), &
            & "nquad increment equals ntriple")
         if (allocated(error)) return
      end do
   end subroutine test_ntriple_nquad

   !> Sorted triples map one to one onto 1..ntriple, in enumeration order (last index slowest)
   subroutine test_index3_bijection(error)
      type(error_type), allocatable, intent(out) :: error
      logical, allocatable :: seen(:)
      integer :: is, ns, i, j, k, idx, counter

      do is = 1, size(tsizes)
         ns = tsizes(is)
         allocate (seen(ntriple_from_ns(ns)))
         seen = .false.
         counter = 0
         do k = 1, ns
            do j = 1, k
               do i = 1, j
                  counter = counter + 1
                  idx = packed_index3(i, j, k)
                  call check(error, idx >= 1 .and. idx <= ntriple_from_ns(ns), "index3 in range")
                  if (allocated(error)) return
                  call check(error, .not. seen(idx), "index3 not hit twice")
                  if (allocated(error)) return
                  seen(idx) = .true.
                  call check(error, idx, counter, "index3 follows enumeration order")
                  if (allocated(error)) return
               end do
            end do
         end do
         call check(error, counter, ntriple_from_ns(ns), "sorted triple count equals ntriple")
         if (allocated(error)) return
         call check(error, all(seen), "index3 onto")
         if (allocated(error)) return
         deallocate (seen)
      end do
   end subroutine test_index3_bijection

   !> Sorted quadruples map one to one onto 1..nquad, in enumeration order (last index slowest)
   subroutine test_index4_bijection(error)
      type(error_type), allocatable, intent(out) :: error
      logical, allocatable :: seen(:)
      integer :: is, ns, i, j, k, l, idx, counter

      do is = 1, size(tsizes)
         ns = tsizes(is)
         allocate (seen(nquad_from_ns(ns)))
         seen = .false.
         counter = 0
         do l = 1, ns
            do k = 1, l
               do j = 1, k
                  do i = 1, j
                     counter = counter + 1
                     idx = packed_index4(i, j, k, l)
                     call check(error, idx >= 1 .and. idx <= nquad_from_ns(ns), "index4 in range")
                     if (allocated(error)) return
                     call check(error, .not. seen(idx), "index4 not hit twice")
                     if (allocated(error)) return
                     seen(idx) = .true.
                     call check(error, idx, counter, "index4 follows enumeration order")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
         call check(error, counter, nquad_from_ns(ns), "sorted quadruple count equals nquad")
         if (allocated(error)) return
         call check(error, all(seen), "index4 onto")
         if (allocated(error)) return
         deallocate (seen)
      end do
   end subroutine test_index4_bijection

   !> Each slab of fixed last index is a contiguous block in rank-2 packed order
   subroutine test_index3_slabs_are_rank2(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ns, i, j, k, offset

      ns = 4
      do k = 1, ns
         offset = packed_index3(1, 1, k) - 1
         call check(error, offset, ntriple_from_ns(k - 1), "slab k starts after all slabs below")
         if (allocated(error)) return
         do j = 1, k
            do i = 1, j
               call check(error, packed_index3(i, j, k) - offset, packed_index(i, j), &
                  & "slab order equals rank-2 order")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_index3_slabs_are_rank2

   !> Each slab of fixed last index is a contiguous block in rank-3 packed order
   subroutine test_index4_slabs_are_rank3(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ns, i, j, k, l, offset

      ns = 4
      do l = 1, ns
         offset = packed_index4(1, 1, 1, l) - 1
         call check(error, offset, nquad_from_ns(l - 1), "slab l starts after all slabs below")
         if (allocated(error)) return
         do k = 1, l
            do j = 1, k
               do i = 1, j
                  call check(error, packed_index4(i, j, k, l) - offset, packed_index3(i, j, k), &
                     & "slab order equals rank-3 order")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_index4_slabs_are_rank3

   !> All six argument orders give the same rank-3 index
   subroutine test_index3_permutation_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: is, ns, i, j, k, ref

      do is = 2, size(tsizes)
         ns = tsizes(is)
         do k = 1, ns
            do j = 1, ns
               do i = 1, ns
                  ref = packed_index3(i, j, k)
                  call check(error, all([packed_index3(i, k, j), packed_index3(j, i, k), &
                     & packed_index3(j, k, i), packed_index3(k, i, j), packed_index3(k, j, i)] == ref), &
                     & "index3 permutation invariant")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_index3_permutation_invariant

   !> All 24 argument orders give the same rank-4 index
   subroutine test_index4_permutation_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: is, ns, i, j, k, l, ref, p1, p2, p3, p4, nperm
      integer :: t(4)

      do is = 2, size(tsizes)
         ns = tsizes(is)
         do l = 1, ns
            do k = 1, ns
               do j = 1, ns
                  do i = 1, ns
                     t = [i, j, k, l]
                     ref = packed_index4(i, j, k, l)
                     nperm = 0
                     do p1 = 1, 4
                        do p2 = 1, 4
                           do p3 = 1, 4
                              do p4 = 1, 4
                                 if (p1 == p2 .or. p1 == p3 .or. p1 == p4 .or. p2 == p3 .or. p2 == p4 &
                                    & .or. p3 == p4) cycle
                                 nperm = nperm + 1
                                 call check(error, packed_index4(t(p1), t(p2), t(p3), t(p4)), ref, &
                                    & "index4 permutation invariant")
                                 if (allocated(error)) return
                              end do
                           end do
                        end do
                     end do
                     call check(error, nperm, 24, "all 24 permutations visited")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_index4_permutation_invariant

   !> Multiplicity of a packed entry equals the number of dense tuples mapping onto it
   subroutine test_multiplicity_counts_dense(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: hits3(:), hits4(:)
      integer :: is, ns, i, j, k, l

      do is = 1, size(tsizes)
         ns = tsizes(is)
         allocate (hits3(ntriple_from_ns(ns)), hits4(nquad_from_ns(ns)))
         hits3 = 0
         hits4 = 0
         do l = 1, ns
            do k = 1, ns
               do j = 1, ns
                  hits3(packed_index3(j, k, l)) = hits3(packed_index3(j, k, l)) + 1
                  do i = 1, ns
                     hits4(packed_index4(i, j, k, l)) = hits4(packed_index4(i, j, k, l)) + 1
                  end do
               end do
            end do
         end do
         do k = 1, ns
            do j = 1, k
               do i = 1, j
                  call check(error, multiplicity3(k, i, j), hits3(packed_index3(i, j, k)), &
                     & "multiplicity3 equals dense hit count")
                  if (allocated(error)) return
                  do l = k, ns
                     call check(error, multiplicity4(l, j, k, i), hits4(packed_index4(i, j, k, l)), &
                        & "multiplicity4 equals dense hit count")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
         deallocate (hits3, hits4)
      end do
   end subroutine test_multiplicity_counts_dense

   !* ================================================================================= *!
   !*                           Rank 3 and 4: pack and unpack                           *!
   !* ================================================================================= *!

   !> Unpack fills every permutation; pack inverts unpack and unpack inverts pack exactly
   subroutine test_tensor3_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: packed(:), packed2(:), full(:, :, :), full2(:, :, :)
      integer :: is, ns, i, j, k

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed3(ns, 23 + is, packed)
         allocate (full(ns, ns, ns), full2(ns, ns, ns), packed2(size(packed)))
         call unpack_sym_tensor3(packed, full, ns)
         do k = 1, ns
            do j = 1, ns
               do i = 1, ns
                  call check(error, all([full(i, k, j), full(j, i, k), full(j, k, i), full(k, i, j), &
                     & full(k, j, i)] == full(i, j, k)), "unpacked rank-3 tensor symmetric")
                  if (allocated(error)) return
               end do
            end do
         end do
         call pack_sym_tensor3(full, packed2, ns)
         call check(error, all(packed2 == packed), "pack(unpack(p)) == p")
         if (allocated(error)) return
         call unpack_sym_tensor3(packed2, full2, ns)
         call check(error, all(full2 == full), "unpack(pack(D)) == D for symmetric D")
         if (allocated(error)) return
         deallocate (full, full2, packed2)
      end do
   end subroutine test_tensor3_roundtrip

   !> Unpack fills every permutation; pack inverts unpack and unpack inverts pack exactly
   subroutine test_tensor4_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: packed(:), packed2(:), full(:, :, :, :), full2(:, :, :, :)
      integer :: is, ns, i, j, k, l

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed4(ns, 29 + is, packed)
         allocate (full(ns, ns, ns, ns), full2(ns, ns, ns, ns), packed2(size(packed)))
         call unpack_sym_tensor4(packed, full, ns)
         do l = 1, ns
            do k = 1, ns
               do j = 1, ns
                  do i = 1, ns
                     ! Generators of S4: adjacent swaps and a cycle
                     call check(error, all([full(j, i, k, l), full(i, k, j, l), full(i, j, l, k), &
                        & full(l, i, j, k)] == full(i, j, k, l)), "unpacked rank-4 tensor symmetric")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
         call pack_sym_tensor4(full, packed2, ns)
         call check(error, all(packed2 == packed), "pack(unpack(p)) == p")
         if (allocated(error)) return
         call unpack_sym_tensor4(packed2, full2, ns)
         call check(error, all(full2 == full), "unpack(pack(D)) == D for symmetric D")
         if (allocated(error)) return
         deallocate (full, full2, packed2)
      end do
   end subroutine test_tensor4_roundtrip

   !> Negative control: an asymmetric dense input does not survive pack then unpack
   subroutine test_tensor_pack_rejects_asymmetric(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: packed(:), full(:, :, :), full2(:, :, :)
      real(wp), allocatable :: packed4(:), full4(:, :, :, :), full42(:, :, :, :)
      integer, parameter :: ns = 3

      call random_packed3(ns, 41, packed)
      allocate (full(ns, ns, ns), full2(ns, ns, ns))
      call unpack_sym_tensor3(packed, full, ns)
      ! Unsorted position, never read by pack
      full(3, 1, 2) = full(3, 1, 2) + 1.0_wp
      call pack_sym_tensor3(full, packed, ns)
      call unpack_sym_tensor3(packed, full2, ns)
      call check(error, maxval(abs(full2 - full)) > thr, "asymmetric rank-3 input detected")
      if (allocated(error)) return
      call check(error, count(full2 /= full), 1, "only the perturbed rank-3 entry differs")
      if (allocated(error)) return

      call random_packed4(ns, 43, packed4)
      allocate (full4(ns, ns, ns, ns), full42(ns, ns, ns, ns))
      call unpack_sym_tensor4(packed4, full4, ns)
      full4(2, 1, 3, 3) = full4(2, 1, 3, 3) + 1.0_wp
      call pack_sym_tensor4(full4, packed4, ns)
      call unpack_sym_tensor4(packed4, full42, ns)
      call check(error, maxval(abs(full42 - full4)) > thr, "asymmetric rank-4 input detected")
      if (allocated(error)) return
      call check(error, count(full42 /= full4), 1, "only the perturbed rank-4 entry differs")
   end subroutine test_tensor_pack_rejects_asymmetric

   !* ================================================================================= *!
   !*                             Rank 2 to 4: contractions                             *!
   !* ================================================================================= *!

   !> H v against the dense matrix-vector product
   subroutine test_contract2_vs_dense(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: ph(:), h(:, :), v(:), out(:), ref(:)
      integer :: is, ns, i

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed(ns, 51 + is, ph)
         allocate (h(ns, ns), v(ns), out(ns), ref(ns))
         call fill_random(v, 151 + is)
         call unpack_sym_matrix(ph, h, ns)
         ref = matmul(h, v)
         call packed_contract2(ph, v, out, ns)
         do i = 1, ns
            call check(error, within_tol(out(i), ref(i)), "H v equals dense product")
            if (allocated(error)) return
         end do
         deallocate (h, v, out, ref)
      end do
   end subroutine test_contract2_vs_dense

   !> T v against the dense contraction over the last index, all (i, j)
   subroutine test_contract3_vs_dense(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: pt(:), t(:, :, :), v(:), pout(:), out(:, :), ref(:, :)
      integer :: is, ns, i, j, k

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed3(ns, 53 + is, pt)
         allocate (t(ns, ns, ns), v(ns), pout(npair_from_ns(ns)), out(ns, ns), ref(ns, ns))
         call fill_random(v, 153 + is)
         call unpack_sym_tensor3(pt, t, ns)
         ref = 0.0_wp
         do k = 1, ns
            do j = 1, ns
               do i = 1, ns
                  ref(i, j) = ref(i, j) + t(i, j, k)*v(k)
               end do
            end do
         end do
         call packed_contract3(pt, v, pout, ns)
         call unpack_sym_matrix(pout, out, ns)
         do j = 1, ns
            do i = 1, ns
               call check(error, within_tol(out(i, j), ref(i, j)), "T v equals dense contraction")
               if (allocated(error)) return
            end do
         end do
         deallocate (t, v, pout, out, ref)
      end do
   end subroutine test_contract3_vs_dense

   !> F v against the dense contraction over the last index, all (i, j, k)
   subroutine test_contract4_vs_dense(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: pf(:), f(:, :, :, :), v(:), pout(:), out(:, :, :), ref(:, :, :)
      integer :: is, ns, i, j, k, l

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed4(ns, 57 + is, pf)
         allocate (f(ns, ns, ns, ns), v(ns), pout(ntriple_from_ns(ns)), out(ns, ns, ns), ref(ns, ns, ns))
         call fill_random(v, 157 + is)
         call unpack_sym_tensor4(pf, f, ns)
         ref = 0.0_wp
         do l = 1, ns
            do k = 1, ns
               do j = 1, ns
                  do i = 1, ns
                     ref(i, j, k) = ref(i, j, k) + f(i, j, k, l)*v(l)
                  end do
               end do
            end do
         end do
         call packed_contract4(pf, v, pout, ns)
         call unpack_sym_tensor3(pout, out, ns)
         do k = 1, ns
            do j = 1, ns
               do i = 1, ns
                  call check(error, within_tol(out(i, j, k), ref(i, j, k)), "F v equals dense contraction")
                  if (allocated(error)) return
               end do
            end do
         end do
         deallocate (f, v, pout, out, ref)
      end do
   end subroutine test_contract4_vs_dense

   !> H:vv, T:vvv, F:vvvv against dense sums and against chained single contractions
   subroutine test_full_contract_vs_dense(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: ph(:), pt(:), pf(:), h(:, :), t(:, :, :), f(:, :, :, :), v(:)
      real(wp), allocatable :: w1(:), w2(:), w3(:)
      real(wp) :: ref2, ref3, ref4
      integer :: is, ns, i, j, k, l

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed(ns, 61 + is, ph)
         call random_packed3(ns, 63 + is, pt)
         call random_packed4(ns, 67 + is, pf)
         allocate (h(ns, ns), t(ns, ns, ns), f(ns, ns, ns, ns), v(ns))
         allocate (w1(ns), w2(npair_from_ns(ns)), w3(ntriple_from_ns(ns)))
         call fill_random(v, 161 + is)
         call unpack_sym_matrix(ph, h, ns)
         call unpack_sym_tensor3(pt, t, ns)
         call unpack_sym_tensor4(pf, f, ns)
         ref2 = 0.0_wp
         ref3 = 0.0_wp
         ref4 = 0.0_wp
         do l = 1, ns
            do k = 1, ns
               ref2 = ref2 + h(k, l)*v(k)*v(l)
               do j = 1, ns
                  ref3 = ref3 + t(j, k, l)*v(j)*v(k)*v(l)
                  do i = 1, ns
                     ref4 = ref4 + f(i, j, k, l)*v(i)*v(j)*v(k)*v(l)
                  end do
               end do
            end do
         end do
         call check(error, within_tol(packed_full_contract2(ph, v, ns), ref2), "H:vv equals dense sum")
         if (allocated(error)) return
         call check(error, within_tol(packed_full_contract3(pt, v, ns), ref3), "T:vvv equals dense sum")
         if (allocated(error)) return
         call check(error, within_tol(packed_full_contract4(pf, v, ns), ref4), "F:vvvv equals dense sum")
         if (allocated(error)) return

         ! Chain of one-index contractions down to a vector
         call packed_contract4(pf, v, w3, ns)
         call packed_contract3(w3, v, w2, ns)
         call packed_contract2(w2, v, w1, ns)
         call check(error, within_tol(dot_product(w1, v), ref4), "chained contractions equal F:vvvv")
         if (allocated(error)) return
         deallocate (h, t, f, v, w1, w2, w3)
      end do
   end subroutine test_full_contract_vs_dense

   !> Full contractions equal the multiplicity-weighted sum over packed entries
   subroutine test_full_contract_vs_multiplicity(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: pt(:), pf(:), v(:)
      real(wp) :: ref3, ref4
      integer :: is, ns, i, j, k, l

      do is = 1, size(tsizes)
         ns = tsizes(is)
         call random_packed3(ns, 71 + is, pt)
         call random_packed4(ns, 73 + is, pf)
         allocate (v(ns))
         call fill_random(v, 171 + is)
         ref3 = 0.0_wp
         ref4 = 0.0_wp
         do l = 1, ns
            do k = 1, l
               do j = 1, k
                  ref3 = ref3 + real(multiplicity3(j, k, l), wp)*pt(packed_index3(j, k, l))*v(j)*v(k)*v(l)
                  do i = 1, j
                     ref4 = ref4 + real(multiplicity4(i, j, k, l), wp)*pf(packed_index4(i, j, k, l)) &
                        & *v(i)*v(j)*v(k)*v(l)
                  end do
               end do
            end do
         end do
         call check(error, within_tol(packed_full_contract3(pt, v, ns), ref3), "T:vvv equals weighted packed sum")
         if (allocated(error)) return
         call check(error, within_tol(packed_full_contract4(pf, v, ns), ref4), "F:vvvv equals weighted packed sum")
         if (allocated(error)) return
         deallocate (v)
      end do
   end subroutine test_full_contract_vs_multiplicity

   !> Negative control: perturbing one packed entry moves the contraction beyond tolerance
   subroutine test_contract_sees_perturbation(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 3
      real(wp), allocatable :: pt(:), pf(:)
      real(wp) :: v(ns), out3(npair_from_ns(ns)), out3b(npair_from_ns(ns))
      real(wp) :: out4(ntriple_from_ns(ns)), out4b(ntriple_from_ns(ns))
      real(wp) :: s3, s4

      call random_packed3(ns, 81, pt)
      call random_packed4(ns, 83, pf)
      call fill_random(v, 181)
      call packed_contract3(pt, v, out3, ns)
      call packed_contract4(pf, v, out4, ns)
      s3 = packed_full_contract3(pt, v, ns)
      s4 = packed_full_contract4(pf, v, ns)

      ! Off-diagonal entries (1,2,3) and (1,2,3,3)
      pt(packed_index3(1, 2, 3)) = pt(packed_index3(1, 2, 3)) + 1.0e-6_wp
      pf(packed_index4(1, 2, 3, 3)) = pf(packed_index4(1, 2, 3, 3)) + 1.0e-6_wp
      call packed_contract3(pt, v, out3b, ns)
      call packed_contract4(pf, v, out4b, ns)
      call check(error, maxval(abs(out3b - out3)) > thr, "rank-3 contraction sees perturbation")
      if (allocated(error)) return
      call check(error, maxval(abs(out4b - out4)) > thr, "rank-4 contraction sees perturbation")
      if (allocated(error)) return
      call check(error, .not. within_tol(packed_full_contract3(pt, v, ns), s3), "T:vvv sees perturbation")
      if (allocated(error)) return
      call check(error, .not. within_tol(packed_full_contract4(pf, v, ns), s4), "F:vvvv sees perturbation")
   end subroutine test_contract_sees_perturbation

   !* ================================================================================= *!
   !*                            Rank 2: one-pass contraction                           *!
   !* ================================================================================= *!

   !> Row-by-row H v through the index function; pre-DSPMV reference
   pure subroutine contract2_rowwise(h_packed, v, out, ns)
      !> Packed matrix
      real(wp), intent(in) :: h_packed(:)
      !> Vector
      real(wp), intent(in) :: v(:)
      !> Contracted vector
      real(wp), intent(out) :: out(:)
      !> Dimension
      integer, intent(in) :: ns
      integer :: i, k
      real(wp) :: acc

      do i = 1, ns
         acc = 0.0_wp
         do k = 1, ns
            acc = acc + h_packed(packed_index(i, k))*v(k)
         end do
         out(i) = acc
      end do
   end subroutine contract2_rowwise

   !> H v against the row-wise reference and dense matmul, 1e-14 of sum_k |H_ik v_k|
   subroutine test_contract2_vs_rowwise(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns_list(7) = [1, 2, 3, 4, 7, 16, 33]
      real(wp), parameter :: rtol = 1.0e-14_wp
      real(wp), allocatable :: ph(:), h(:, :), v(:), out(:), scale(:)
      real(wp), allocatable :: ref_row(:), ref_mm(:)
      integer :: is, ns

      do is = 1, size(ns_list)
         ns = ns_list(is)
         call random_packed(ns, 251 + is, ph)
         allocate (h(ns, ns), v(ns), out(ns), ref_row(ns), ref_mm(ns), scale(ns))
         call fill_random(v, 351 + is)
         call unpack_sym_matrix(ph, h, ns)
         call contract2_rowwise(ph, v, ref_row, ns)
         ref_mm = matmul(h, v)
         scale = matmul(abs(h), abs(v))
         call packed_contract2(ph, v, out, ns)
         call check(error, all(abs(out - ref_row) <= rtol*scale), &
            & "H v equals row-wise reference")
         if (allocated(error)) return
         call check(error, all(abs(out - ref_mm) <= rtol*scale), &
            & "H v equals dense matmul")
         if (allocated(error)) return
         deallocate (h, v, out, ref_row, ref_mm, scale)
      end do
   end subroutine test_contract2_vs_rowwise

   !> Negative control: one off-diagonal packed entry moves both out(i) and out(j)
   subroutine test_contract2_updates_both(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns = 5
      real(wp), parameter :: delta = 1.0e-6_wp
      real(wp), allocatable :: ph(:)
      real(wp) :: v(ns), out(ns), out2(ns), expect(ns)
      integer :: i, j

      call random_packed(ns, 451, ph)
      call fill_random(v, 452)
      call packed_contract2(ph, v, out, ns)
      do j = 2, ns
         do i = 1, j - 1
            ph(packed_index(i, j)) = ph(packed_index(i, j)) + delta
            call packed_contract2(ph, v, out2, ns)
            ph(packed_index(i, j)) = ph(packed_index(i, j)) - delta
            expect = 0.0_wp
            expect(i) = delta*v(j)
            expect(j) = delta*v(i)
            call check(error, all(abs(out2 - out - expect) <= thr), &
               & "off-diagonal entry enters out(i) and out(j) only")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_contract2_updates_both

   !* ================================================================================= *!
   !*                            Rank 3 and 4: pack placement                           *!
   !* ================================================================================= *!

   !> Pseudo-random value from a seed, deterministic, in [-1, 1)
   pure function seeded_value(seed) result(x)
      !> Seed
      integer(int64), intent(in) :: seed
      !> Value
      real(wp) :: x
      integer(int64) :: state
      integer :: m

      state = mod(seed, 2147483648_int64)
      do m = 1, 3
         state = mod(state*1103515245_int64 + 12345_int64, 2147483648_int64)
      end do
      x = 2.0_wp*real(state, wp)/2147483648.0_wp - 1.0_wp
   end function seeded_value

   !> Pack places every dense entry at the index-function slot; unpack restores it
   !>
   !> - symmetric input from elementary symmetric polynomials of the tuple, no sorting
   !> - negative control: an asymmetric input does not survive pack then unpack
   subroutine test_tensor_pack_matches_index(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ns_list(2) = [3, 5]
      real(wp), allocatable :: f3(:, :, :), g3(:, :, :), p3(:)
      real(wp), allocatable :: f4(:, :, :, :), g4(:, :, :, :), p4(:), flat(:)
      integer(int64) :: e1, e2, e3, e4
      integer :: is, ns, i, j, k, l

      do is = 1, size(ns_list)
         ns = ns_list(is)
         allocate (f3(ns, ns, ns), g3(ns, ns, ns), p3(ntriple_from_ns(ns)))
         allocate (f4(ns, ns, ns, ns), g4(ns, ns, ns, ns), p4(nquad_from_ns(ns)))
         do l = 1, ns
            do k = 1, ns
               do j = 1, ns
                  e1 = j + k + l
                  e2 = j*k + k*l + l*j
                  e3 = j*k*l
                  f3(j, k, l) = seeded_value(e1 + 101*e2 + 10007*e3)
                  do i = 1, ns
                     e1 = i + j + k + l
                     e2 = i*j + i*k + i*l + j*k + j*l + k*l
                     e3 = i*j*k + i*j*l + i*k*l + j*k*l
                     e4 = i*j*k*l
                     f4(i, j, k, l) = seeded_value(e1 + 101*e2 + 10007*e3 + 1000003*e4)
                  end do
               end do
            end do
         end do

         call pack_sym_tensor3(f3, p3, ns)
         call pack_sym_tensor4(f4, p4, ns)
         do l = 1, ns
            do k = 1, ns
               do j = 1, ns
                  call check(error, p3(packed_index3(j, k, l)) == f3(j, k, l), &
                     & "rank-3 pack slot equals index function")
                  if (allocated(error)) return
                  do i = 1, ns
                     call check(error, p4(packed_index4(i, j, k, l)) == f4(i, j, k, l), &
                        & "rank-4 pack slot equals index function")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
         call unpack_sym_tensor3(p3, g3, ns)
         call check(error, all(g3 == f3), "rank-3 unpack(pack(D)) == D")
         if (allocated(error)) return
         call unpack_sym_tensor4(p4, g4, ns)
         call check(error, all(g4 == f4), "rank-4 unpack(pack(D)) == D")
         if (allocated(error)) return

         ! Negative control: dense random input, not symmetric
         allocate (flat(ns**4))
         call fill_random(flat, 551 + is)
         f3 = reshape(flat(:ns**3), [ns, ns, ns])
         f4 = reshape(flat, [ns, ns, ns, ns])
         call pack_sym_tensor3(f3, p3, ns)
         call unpack_sym_tensor3(p3, g3, ns)
         call check(error, any(g3 /= f3), "asymmetric rank-3 input not restored")
         if (allocated(error)) return
         call pack_sym_tensor4(f4, p4, ns)
         call unpack_sym_tensor4(p4, g4, ns)
         call check(error, any(g4 /= f4), "asymmetric rank-4 input not restored")
         if (allocated(error)) return
         deallocate (f3, g3, p3, f4, g4, p4, flat)
      end do
   end subroutine test_tensor_pack_matches_index

   !* ================================================================================= *!
   !*                          Large ns: int64 index arithmetic                         *!
   !* ================================================================================= *!

   !> Simplex number C(n + r - 1, r) in int64, one exact division per factor
   pure function simplex64(r, n) result(s)
      !> Rank
      integer, intent(in) :: r
      !> Dimension, nonnegative
      integer, intent(in) :: n
      !> Number of sorted r-tuples over n sites
      integer(int64) :: s
      integer :: m

      s = 1
      do m = 1, r
         s = s*(n + m - 1)/m
      end do
   end function simplex64

   !> Packed index of a sorted tuple, 1 + sum_m C(t(m) + m - 2, m), in int64
   pure function index64(t) result(idx)
      !> Sorted tuple, ascending
      integer, intent(in) :: t(:)
      !> Packed index
      integer(int64) :: idx
      integer :: m

      idx = 1
      do m = 1, size(t)
         idx = idx + simplex64(m, t(m) - 1)
      end do
   end function index64

   !> Largest ns whose rank-r count fits in a default integer
   pure function largest_ns(r) result(ns)
      !> Rank
      integer, intent(in) :: r
      !> Largest dimension
      integer :: ns

      ns = 1
      do while (simplex64(r, ns + 1) <= huge(1))
         ns = ns + 1
      end do
   end function largest_ns

   !> Counts at the largest fitting ns equal the int64 reference; cited ns = 215 value
   subroutine test_count_boundary(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: n2, n3, n4

      n2 = largest_ns(2)
      n3 = largest_ns(3)
      n4 = largest_ns(4)
      call check(error, simplex64(4, n4) <= huge(1) .and. &
         & simplex64(4, n4 + 1) > huge(1), &
         & "reference boundary brackets huge(1)")
      if (allocated(error)) return
      call check(error, int(npair_from_ns(n2), int64), simplex64(2, n2), &
         & "npair at boundary")
      if (allocated(error)) return
      call check(error, int(ntriple_from_ns(n3), int64), simplex64(3, n3), &
         & "ntriple at boundary")
      if (allocated(error)) return
      call check(error, int(nquad_from_ns(n4), int64), simplex64(4, n4), &
         & "nquad at boundary")
      if (allocated(error)) return
      call check(error, nquad_from_ns(215), 91537110, "nquad(215)")
      if (allocated(error)) return
      call check(error, nquad_from_ns(214), int(simplex64(4, 214)), "nquad(214)")
   end subroutine test_count_boundary

   !> Indices at and one past the largest fitting ns equal the int64 reference
   !>
   !> - one past: first entry of the next slab, still within huge(1)
   !> - arguments unsorted to exercise the sort with large values
   subroutine test_index_boundary(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: n2, n3, n4

      n2 = largest_ns(2)
      n3 = largest_ns(3)
      n4 = largest_ns(4)

      call check(error, int(packed_index(n2, n2), int64), index64([n2, n2]), &
         & "packed_index at boundary")
      if (allocated(error)) return
      call check(error, int(packed_index(n2 + 1, 1), int64), index64([1, n2 + 1]), &
         & "packed_index one past boundary")
      if (allocated(error)) return

      call check(error, int(packed_index3(n3, n3, n3), int64), index64([n3, n3, n3]), &
         & "packed_index3 at boundary")
      if (allocated(error)) return
      call check(error, int(packed_index3(n3 + 1, 1, 2), int64), &
         & index64([1, 2, n3 + 1]), &
         & "packed_index3 one past boundary")
      if (allocated(error)) return

      call check(error, int(packed_index4(n4, n4, n4, n4), int64), &
         & index64([n4, n4, n4, n4]), "packed_index4 at boundary")
      if (allocated(error)) return
      call check(error, int(packed_index4(3, n4 + 1, 1, 2), int64), &
         & index64([1, 2, 3, n4 + 1]), "packed_index4 one past boundary")
      if (allocated(error)) return
      call check(error, int(packed_index4(216, 216, 216, 216), int64), &
         & index64([216, 216, 216, 216]), "packed_index4(216, 216, 216, 216)")
      if (allocated(error)) return
      call check(error, int(packed_index4(216, 3, 215, 100), int64), &
         & index64([3, 100, 215, 216]), "packed_index4 mixed tuple at ns = 216")
   end subroutine test_index_boundary

end module test_math_packed
