!> Test suite for packed symmetric matrix utilities in moist_math_packed_sym
module test_math_packed_sym
   use mctc_env, only: wp
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_packed_sym, only: npair_from_ns, packed_index, packed_index_lut, unpack_sym_matrix, &
      & pack_sym_matrix, packed_symmatmul, packed_symmatmul_lut
   use moist_math_lapack_sptrs, only: wrap_dsptrf, wrap_dsptrs
   implicit none(type, external)
   private

   public :: collect_math_packed_sym

   !> Round-off level for sums of at most 7 products of O(1) numbers
   real(wp), parameter :: thr = 1.0e-13_wp
   integer, parameter :: sizes(4) = [1, 2, 3, 7]

contains

   !> Collect all packed symmetric tests
   subroutine collect_math_packed_sym(testsuite)
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
                  new_unittest("symmatmul_vs_dense", test_symmatmul_vs_dense), &
                  new_unittest("symmatmul_symmetric_case", test_symmatmul_symmetric_case), &
                  new_unittest("symmatmul_lut_vs_dense", test_symmatmul_lut_vs_dense), &
                  new_unittest("layout_vs_lapack_sptrs", test_layout_vs_lapack_sptrs) &
                  ]
   end subroutine collect_math_packed_sym

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

end module test_math_packed_sym
