!> Tests of the sorting routines in moist_math_sorter
!>
!> Reference: the definition of each order, no second sorting implementation
!>
!> - stable order: entry `i` at `1 + #{keys ahead of it in the order}
!>   + #{equal keys before it in the input}` ([[stable_reference_int]])
!> - unstable [[sort]]: sorted values only
!>
!> Fixtures:
!>
!> - n = 0 .. 2049; merges from 64 entries on, introsort partitions from 17
!> - random, sorted, reversed, few-valued, constant, sawtooth (many natural
!>   runs), organ-pipe, uneven blocks
!> - heap sort on mixed values: organ-pipe and uneven blocks only; constant,
!>   few-valued and sawtooth reach it with constant slices
!> - backward merges and every collapse rule: uneven blocks only, the one input
!>   with natural runs above the minimum run
!>
!> Error paths: none; a too-small workspace is an `error stop`
module test_math_sorters
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_negative
   use mctc_env, only: wp
   use test_helpers, only: det_rand
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_sorter, only: ord_sort, sort, sort_index, sort_adjoint, counting_argsort
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_math_sorters

   !> Fixture sizes: trivial, insertion-only, around the 64-entry merge threshold, merging
   integer, parameter :: sizes(*) = [0, 1, 2, 3, 16, 17, 31, 63, 64, 65, 127, 200, 1000, 2049]
   !> Number of input patterns
   integer, parameter :: npattern = 8
   !> Input pattern names, for failure messages
   character(len=10), parameter :: pattern_name(npattern) = [character(len=10) :: &
                                                             "random", "sorted", "reversed", "few-valued", &
                                                             "constant", "sawtooth", "organ-pipe", "blocks"]
   !> Block lengths of the uneven-block pattern, cycled
   !>
   !> At n = 2049: backward merges and every collapse rule
   integer, parameter :: block_len(*) = [267, 152, 56, 59, 133, 330, 71, 110, 92, 274, 151, 277]
   !> Sentinel of an unused workspace
   integer, parameter :: unused = -huge(1)

contains

   !> Collect all sorter tests
   subroutine collect_math_sorters(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ord_sort_int_ascending", test_ord_sort_int_ascending), &
                  new_unittest("ord_sort_real_ascending", test_ord_sort_real_ascending), &
                  new_unittest("ord_sort_real_stable", test_ord_sort_real_stable), &
                  new_unittest("ord_sort_reverse", test_ord_sort_reverse), &
                  new_unittest("ord_sort_work", test_ord_sort_work), &
                  new_unittest("sort_int_ascending", test_sort_int_ascending), &
                 new_unittest("sort_int64_ascending", test_sort_int64_ascending), &
                  new_unittest("sort_real_ascending", test_sort_real_ascending), &
                  new_unittest("sort_reverse", test_sort_reverse), &
                  new_unittest("sort_index_real", test_sort_index_real), &
                  new_unittest("sort_index_int", test_sort_index_int), &
                  new_unittest("sort_index_reverse", test_sort_index_reverse), &
                  new_unittest("sort_index_work", test_sort_index_work), &
                  new_unittest("sort_adjoint_companion", test_sort_adjoint_companion), &
                  new_unittest("counting_argsort_stable", test_counting_argsort_stable), &
                  new_unittest("counting_argsort_edges", test_counting_argsort_edges) &
                  ]
   end subroutine collect_math_sorters

   !* ================================================================================= *!
   !*                                     ord_sort                                      *!
   !* ================================================================================= *!

   !> Integer ord_sort: stable ascending order on every fixture
   subroutine test_ord_sort_int_ascending(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: a(:), orig(:)
      integer :: is, ip, n, nmerge
      !> Result of `label`
      character(len=:), allocatable :: tag

      nmerge = 0
      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = int_pattern(n, ip)
            if (n >= 64 .and. .not. is_sorted_int(orig)) nmerge = nmerge + 1
            a = orig
            call ord_sort(a)
            call label(n, ip, tag)
            call check(error, all(a == orig(stable_reference_int(orig, .false.))), &
                       "ord_sort (integer) is not the ascending order"//tag)
            if (allocated(error)) return
         end do
      end do
      call check(error, nmerge > 0, "fixture never reaches the merges of ord_sort")
   end subroutine test_ord_sort_int_ascending

   !> Real ord_sort: stable ascending order on every fixture
   subroutine test_ord_sort_real_ascending(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:), orig(:)
      integer :: is, ip, n, nmerge
      !> Result of `label`
      character(len=:), allocatable :: tag

      nmerge = 0
      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = real_pattern(n, ip)
            if (n >= 64 .and. .not. is_sorted_real(orig)) nmerge = nmerge + 1
            a = orig
            call ord_sort(a)
            call label(n, ip, tag)
            call check(error, all(a == orig(stable_reference_real(orig, .false.))), &
                       "ord_sort (real) is not the ascending order"//tag)
            if (allocated(error)) return
         end do
      end do
      call check(error, nmerge > 0, "fixture never reaches the merges of ord_sort")
   end subroutine test_ord_sort_real_ascending

   !> Equal keys in input order: signed zeros, equal yet distinguishable
   subroutine test_ord_sort_real_stable(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ns(*) = [40, 300, 2049]
      integer :: is

      do is = 1, size(ns)
         call check_stable_zeros(error, signed_zero_fixture(ns(is)), .false., "ord_sort, n = "//to_string(ns(is)))
         if (allocated(error)) return
      end do
   end subroutine test_ord_sort_real_stable

   !> `reverse`: non-increasing, equal keys still in input order
   subroutine test_ord_sort_reverse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: ia(:), iorig(:)
      real(wp), allocatable :: ra(:), rorig(:)
      integer :: is, ip, n
      !> Result of `label`
      character(len=:), allocatable :: tag

      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            iorig = int_pattern(n, ip)
            ia = iorig
            call ord_sort(ia, reverse=.true.)
            call label(n, ip, tag)
            call check(error, all(ia == iorig(stable_reference_int(iorig, .true.))), &
                       "ord_sort (integer, reverse) is not the descending order"//tag)
            if (allocated(error)) return
            rorig = real_pattern(n, ip)
            ra = rorig
            call ord_sort(ra, reverse=.true.)
            call check(error, all(ra == rorig(stable_reference_real(rorig, .true.))), &
                       "ord_sort (real, reverse) is not the descending order"//tag)
            if (allocated(error)) return
         end do
      end do
      call check_stable_zeros(error, signed_zero_fixture(300), .true., "ord_sort, reverse")
   end subroutine test_ord_sort_reverse

   !> Caller workspace of `n/2` entries: same result, workspace written
   subroutine test_ord_sort_work(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ns(*) = [65, 1000, 2049]
      integer, allocatable :: ia(:), ib(:), iwork(:)
      real(wp), allocatable :: ra(:), rb(:), rwork(:)
      integer :: is, n

      do is = 1, size(ns)
         n = ns(is)
         ia = int_pattern(n, 1)
         ib = ia
         allocate (iwork(n/2), source=unused)
         call ord_sort(ia, work=iwork)
         call ord_sort(ib)
         call check(error, all(ia == ib), "ord_sort (integer) with work differs, n = "//to_string(n))
         if (allocated(error)) return
         call check(error, any(iwork /= unused), "ord_sort (integer) left the caller workspace unused, n = " &
                    //to_string(n))
         if (allocated(error)) return

         ra = real_pattern(n, 1)
         rb = ra
         allocate (rwork(n/2), source=real(unused, wp))
         call ord_sort(ra, work=rwork)
         call ord_sort(rb)
         call check(error, all(ra == rb), "ord_sort (real) with work differs, n = "//to_string(n))
         if (allocated(error)) return
         call check(error, any(rwork /= real(unused, wp)), &
                    "ord_sort (real) left the caller workspace unused, n = "//to_string(n))
         if (allocated(error)) return
         deallocate (iwork, rwork)
      end do
   end subroutine test_ord_sort_work

   !* ================================================================================= *!
   !*                                       sort                                        *!
   !* ================================================================================= *!

   !> Integer introsort: ascending values on every fixture
   subroutine test_sort_int_ascending(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: a(:), orig(:)
      integer :: is, ip, n, ndeep
      !> Result of `label`
      character(len=:), allocatable :: tag

      ndeep = 0
      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = int_pattern(n, ip)
            ! Few-valued, constant: depth limit reached, heap sort on constant slices only
            if (n >= 1000 .and. (ip == 4 .or. ip == 5)) ndeep = ndeep + 1
            a = orig
            call sort(a)
            call label(n, ip, tag)
            call check(error, all(a == orig(stable_reference_int(orig, .false.))), &
                       "sort (integer) is not the ascending order"//tag)
            if (allocated(error)) return
         end do
      end do
      call check(error, ndeep > 0, "fixture never drives sort into its heap sort")
   end subroutine test_sort_int_ascending

   !> 64-bit integer introsort: ascending values on every fixture, keys beyond the
   !> default-integer range
   subroutine test_sort_int64_ascending(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer(int64), parameter :: shift = 4294967296_int64
      integer, allocatable :: orig(:)
      integer(int64), allocatable :: a(:)
      integer :: is, ip, n

      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = int_pattern(n, ip)
            a = int(orig, int64)*shift - 7_int64
            call sort(a)
            call check(error, all(a == int(orig(stable_reference_int(orig, .false.)), int64)*shift - 7_int64), &
                       "sort (64-bit integer) is not the ascending order"//label(n, ip))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_sort_int64_ascending

   !> Real introsort: ascending values on every fixture
   subroutine test_sort_real_ascending(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:), orig(:)
      integer :: is, ip, n, npart
      !> Result of `label`
      character(len=:), allocatable :: tag

      npart = 0
      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = real_pattern(n, ip)
            if (n > 16 .and. .not. is_sorted_real(orig)) npart = npart + 1
            a = orig
            call sort(a)
            call label(n, ip, tag)
            call check(error, all(a == orig(stable_reference_real(orig, .false.))), &
                       "sort (real) is not the ascending order"//tag)
            if (allocated(error)) return
         end do
      end do
      call check(error, npart > 0, "fixture never reaches the partitions of sort")
   end subroutine test_sort_real_ascending

   !> `reverse`: non-increasing values
   subroutine test_sort_reverse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: ia(:), iorig(:)
      real(wp), allocatable :: ra(:), rorig(:)
      integer :: is, ip, n
      !> Result of `label`
      character(len=:), allocatable :: tag

      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            iorig = int_pattern(n, ip)
            ia = iorig
            call sort(ia, reverse=.true.)
            call label(n, ip, tag)
            call check(error, all(ia == iorig(stable_reference_int(iorig, .true.))), &
                       "sort (integer, reverse) is not the descending order"//tag)
            if (allocated(error)) return
            rorig = real_pattern(n, ip)
            ra = rorig
            call sort(ra, reverse=.true.)
            call check(error, all(ra == rorig(stable_reference_real(rorig, .true.))), &
                       "sort (real, reverse) is not the descending order"//tag)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_sort_reverse

   !* ================================================================================= *!
   !*                              sort_index, sort_adjoint                             *!
   !* ================================================================================= *!

   !> Real keys: index is the stable permutation and reproduces the sorted keys
   subroutine test_sort_index_real(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: a(:), orig(:)
      integer, allocatable :: idx(:)
      integer :: is, ip, n, nties
      !> Result of `label`
      character(len=:), allocatable :: tag

      nties = 0
      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = real_pattern(n, ip)
            if (n >= 64 .and. has_ties_real(orig) .and. .not. is_sorted_real(orig)) nties = nties + 1
            a = orig
            allocate (idx(n), source=0)
            call sort_index(a, idx)
            call label(n, ip, tag)
            call check_permutation(error, idx, "sort_index (real)"//tag)
            if (allocated(error)) return
            call check(error, all(orig(idx) == a), "sort_index (real): index does not reproduce the keys" &
                       //tag)
            if (allocated(error)) return
            call check(error, all(idx == stable_reference_real(orig, .false.)), &
                       "sort_index (real) is not the stable ascending order"//tag)
            if (allocated(error)) return
            deallocate (idx)
         end do
      end do
      call check(error, nties > 0, "fixture never merges unsorted keys with ties")
   end subroutine test_sort_index_real

   !> Integer keys: index is the stable permutation and reproduces the sorted keys
   subroutine test_sort_index_int(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: a(:), orig(:), idx(:)
      integer :: is, ip, n, nties
      !> Result of `label`
      character(len=:), allocatable :: tag

      nties = 0
      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            orig = int_pattern(n, ip)
            if (n >= 64 .and. has_ties_int(orig) .and. .not. is_sorted_int(orig)) nties = nties + 1
            a = orig
            allocate (idx(n), source=0)
            call sort_index(a, idx)
            call label(n, ip, tag)
            call check_permutation(error, idx, "sort_index (integer)"//tag)
            if (allocated(error)) return
            call check(error, all(orig(idx) == a), "sort_index (integer): index does not reproduce the keys" &
                       //tag)
            if (allocated(error)) return
            call check(error, all(idx == stable_reference_int(orig, .false.)), &
                       "sort_index (integer) is not the stable ascending order"//tag)
            if (allocated(error)) return
            deallocate (idx)
         end do
      end do
      call check(error, nties > 0, "fixture never merges unsorted keys with ties")
   end subroutine test_sort_index_int

   !> `reverse`: stable descending permutation, equal keys in ascending index
   subroutine test_sort_index_reverse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: ra(:), rorig(:)
      integer, allocatable :: ia(:), iorig(:), idx(:)
      integer :: is, ip, n
      !> Result of `label`
      character(len=:), allocatable :: tag

      do is = 1, size(sizes)
         n = sizes(is)
         do ip = 1, npattern
            allocate (idx(n), source=0)
            rorig = real_pattern(n, ip)
            ra = rorig
            call sort_index(ra, idx, reverse=.true.)
            call label(n, ip, tag)
            call check(error, all(rorig(idx) == ra) .and. all(idx == stable_reference_real(rorig, .true.)), &
                       "sort_index (real, reverse) is not the stable descending order"//tag)
            if (allocated(error)) return
            iorig = int_pattern(n, ip)
            ia = iorig
            call sort_index(ia, idx, reverse=.true.)
            call check(error, all(iorig(idx) == ia) .and. all(idx == stable_reference_int(iorig, .true.)), &
                       "sort_index (integer, reverse) is not the stable descending order"//tag)
            if (allocated(error)) return
            deallocate (idx)
         end do
      end do
   end subroutine test_sort_index_reverse

   !> Caller workspaces of `n/2` entries, both or one: same permutation, workspaces written
   subroutine test_sort_index_work(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ns(*) = [65, 1000, 2049]
      real(wp), allocatable :: a(:), orig(:), work(:)
      integer, allocatable :: idx(:), ref(:), iwork(:)
      integer :: is, n, variant
      character(len=:), allocatable :: what

      do is = 1, size(ns)
         n = ns(is)
         orig = real_pattern(n, 4)
         ref = stable_reference_real(orig, .false.)
         allocate (idx(n), work(n/2), iwork(n/2))
         do variant = 1, 3
            a = orig
            idx = 0
            work = real(unused, wp)
            iwork = unused
            select case (variant)
            case (1)
               what = "work and iwork"
               call sort_index(a, idx, work=work, iwork=iwork)
            case (2)
               what = "work"
               call sort_index(a, idx, work=work)
            case default
               what = "iwork"
               call sort_index(a, idx, iwork=iwork)
            end select
            call check(error, all(idx == ref), "sort_index with "//what//" is not the stable order, n = " &
                       //to_string(n))
            if (allocated(error)) return
            if (variant /= 3) then
               call check(error, any(work /= real(unused, wp)), "sort_index left work unused, n = "//to_string(n))
               if (allocated(error)) return
            end if
            if (variant /= 2) then
               call check(error, any(iwork /= unused), "sort_index left iwork unused, n = "//to_string(n))
               if (allocated(error)) return
            end if
         end do
         deallocate (idx, work, iwork)
      end do
   end subroutine test_sort_index_work

   !> Companion array permuted with the keys, stably, for both key kinds
   subroutine test_sort_adjoint_companion(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ns(*) = [17, 200, 1000]
      real(wp), allocatable :: ra(:), rorig(:)
      integer, allocatable :: ia(:), iorig(:), lab(:), lab0(:)
      integer :: is, ip, n, i
      !> Result of `label`
      character(len=:), allocatable :: tag

      do is = 1, size(ns)
         n = ns(is)
         lab0 = [(1000 + 7*i, i=1, n)]
         do ip = 1, npattern
            rorig = real_pattern(n, ip)
            ra = rorig
            lab = lab0
            call sort_adjoint(ra, lab)
            call label(n, ip, tag)
            call check(error, all(lab == lab0(stable_reference_real(rorig, .false.))) .and. &
                       all(ra == rorig(stable_reference_real(rorig, .false.))), &
                       "sort_adjoint (real) lost the companion order"//tag)
            if (allocated(error)) return
            iorig = int_pattern(n, ip)
            ia = iorig
            lab = lab0
            call sort_adjoint(ia, lab)
            call check(error, all(lab == lab0(stable_reference_int(iorig, .false.))) .and. &
                       all(ia == iorig(stable_reference_int(iorig, .false.))), &
                       "sort_adjoint (integer) lost the companion order"//tag)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_sort_adjoint_companion

   !* ================================================================================= *!
   !*                                 counting_argsort                                  *!
   !* ================================================================================= *!

   !> Stable permutation against an explicit answer and the stable order of random keys
   subroutine test_counting_argsort_stable(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 1000, nkey = 64
      integer :: buckets(10), perm10(10), expected(10), counts6(0:5)
      integer :: keys(n), perm(n), counts(0:nkey - 1)
      real(wp) :: r(n)
      integer :: state

      buckets = [5, 0, 3, 5, 0, 1, 3, 1, 5, 0]
      ! Keys ascending, original indices ascending within each:
      ! 0: {2, 5, 10}, 1: {6, 8}, 3: {3, 7}, 5: {1, 4, 9}
      expected = [2, 5, 10, 6, 8, 3, 7, 1, 4, 9]
      perm10 = -99
      call counting_argsort(buckets, perm10, counts6)
      call check(error, all(perm10 == expected), "counting_argsort order, permutation or stability lost")
      if (allocated(error)) return

      state = 64
      call det_rand(r, state, 0.0_wp, real(nkey, wp))
      keys = min(nkey - 1, int(r))
      call check(error, has_ties_int(keys) .and. maxval(keys) == nkey - 1, &
                 "fixture needs ties and the largest key")
      if (allocated(error)) return
      perm = -99
      call counting_argsort(keys, perm, counts)
      call check(error, all(perm == stable_reference_int(keys, .false.)), &
                 "counting_argsort is not the stable ascending order of random keys")
   end subroutine test_counting_argsort_stable

   !> Empty, singleton at the largest key, single key value
   subroutine test_counting_argsort_edges(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer :: b0(0), p0(0), b1(1), p1(1), b4(4), p4(4), c1(0:0), c8(0:7)

      call counting_argsort(b0, p0, c1)
      call counting_argsort(b0, p0, c8)
      b1 = 7
      p1 = -99
      call counting_argsort(b1, p1, c8)
      call check(error, p1(1) == 1, "Singleton at the largest key must map to index one")
      if (allocated(error)) return
      b4 = 0
      p4 = -99
      call counting_argsort(b4, p4, c1)
      call check(error, all(p4 == [1, 2, 3, 4]), "A single key value must keep the input order")
   end subroutine test_counting_argsort_edges

   !* ================================================================================= *!
   !*                                      Helpers                                      *!
   !* ================================================================================= *!

   !> Integer input of a pattern
   !>
   !> @param[in] n       Number of entries
   !> @param[in] pattern Pattern id, 1..npattern
   function int_pattern(n, pattern) result(a)
      !> Number of entries
      integer, intent(in) :: n
      !> Pattern id
      integer, intent(in) :: pattern
      !> Input
      integer :: a(n)

      real(wp) :: r(n)
      integer :: state, i, k, pos, nblock

      select case (pattern)
      case (1)
         ! Random in -n..n, with occasional ties
         state = 1009 + n
         call det_rand(r, state, -1.0_wp, 1.0_wp)
         a = nint(r*real(n, wp))
      case (2)
         a = [(i, i=1, n)]
      case (3)
         a = [(n - i, i=1, n)]
      case (4)
         a = [(mod(7*i, 5), i=1, n)]
      case (5)
         a = 3
      case (6)
         ! Ascending runs of seven, equal keys across runs
         a = [(mod(i, 7), i=1, n)]
      case (7)
         a = [(min(i, n + 1 - i), i=1, n)]
      case default
         ! Ascending blocks of uneven lengths, each restarting at 1, equal keys across blocks
         pos = 0
         k = 0
         do while (pos < n)
            k = mod(k, size(block_len)) + 1
            nblock = min(block_len(k), n - pos)
            a(pos + 1:pos + nblock) = [(i, i=1, nblock)]
            pos = pos + nblock
         end do
      end select
   end function int_pattern

   !> Real input of a pattern, the integer one scaled except for distinct random values
   !>
   !> @param[in] n       Number of entries
   !> @param[in] pattern Pattern id, 1..npattern
   function real_pattern(n, pattern) result(a)
      !> Number of entries
      integer, intent(in) :: n
      !> Pattern id
      integer, intent(in) :: pattern
      !> Input
      real(wp) :: a(n)

      integer :: state

      if (pattern == 1) then
         state = 4099 + n
         call det_rand(a, state, -4.2_wp, 2.7_wp)
      else
         a = 0.25_wp*real(int_pattern(n, pattern), wp) - 1.5_wp
      end if
   end function real_pattern

   !> Values -1, 0, 1 in random order with zeros of random sign
   !>
   !> @param[in] n Number of entries
   function signed_zero_fixture(n) result(a)
      !> Number of entries
      integer, intent(in) :: n
      !> Input
      real(wp) :: a(n)

      real(wp) :: r(n), s(n)
      integer :: state, i

      state = 77 + n
      call det_rand(r, state, 0.0_wp, 4.0_wp)
      call det_rand(s, state)
      do i = 1, n
         select case (int(r(i)))
         case (0)
            a(i) = -1.0_wp
         case (3)
            a(i) = 1.0_wp
         case default
            a(i) = sign(0.0_wp, s(i) - 0.5_wp)
         end select
      end do
   end function signed_zero_fixture

   !> Sort the signed-zero fixture and check the zeros keep their input sign sequence
   !>
   !> - sign via `ieee_is_negative`
   !> - not `sign(1.0_wp, a)`: +1 for -0.0 under ifx defaults (`-assume nominus0`)
   !>
   !> @param[out] error   Test failure
   !> @param[in]  orig    Signed-zero fixture, zeros of both signs
   !> @param[in]  reverse Descending sort
   !> @param[in]  what    Label of the call
   subroutine check_stable_zeros(error, orig, reverse, what)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Signed-zero fixture
      real(wp), intent(in) :: orig(:)
      !> Descending sort
      logical, intent(in) :: reverse
      !> Label of the call
      character(len=*), intent(in) :: what

      real(wp), allocatable :: a(:)
      logical, allocatable :: neg_in(:), neg_out(:)

      neg_in = pack(ieee_is_negative(orig), orig == 0.0_wp)
      call check(error, any(neg_in) .and. .not. all(neg_in), what//": fixture needs zeros of both signs")
      if (allocated(error)) return
      a = orig
      call ord_sort(a, reverse=reverse)
      if (reverse) then
         call check(error, all(a(:size(a) - 1) >= a(2:)), what//": not non-increasing")
      else
         call check(error, all(a(:size(a) - 1) <= a(2:)), what//": not non-decreasing")
      end if
      if (allocated(error)) return
      neg_out = pack(ieee_is_negative(a), a == 0.0_wp)
      call check(error, size(neg_in) == size(neg_out), what//": zeros lost")
      if (allocated(error)) return
      call check(error, all(neg_in .eqv. neg_out), what//": equal keys (signed zeros) left their input order")
   end subroutine check_stable_zeros

   !> Stable order by definition
   !>
   !> - entry `i` at `1 + #{keys ahead} + #{equal keys before i}`
   !>
   !> @param[in] keys       Keys
   !> @param[in] descending Non-increasing order
   pure function stable_reference_int(keys, descending) result(perm)
      !> Keys
      integer, intent(in) :: keys(:)
      !> Non-increasing order
      logical, intent(in) :: descending
      !> Input position of each sorted entry
      integer :: perm(size(keys))

      integer :: i, pos

      do i = 1, size(keys)
         if (descending) then
            pos = 1 + count(keys > keys(i)) + count(keys(:i - 1) == keys(i))
         else
            pos = 1 + count(keys < keys(i)) + count(keys(:i - 1) == keys(i))
         end if
         perm(pos) = i
      end do
   end function stable_reference_int

   !> Stable order by definition
   !>
   !> - entry `i` at `1 + #{keys ahead} + #{equal keys before i}`
   !>
   !> @param[in] keys       Keys
   !> @param[in] descending Non-increasing order
   pure function stable_reference_real(keys, descending) result(perm)
      !> Keys
      real(wp), intent(in) :: keys(:)
      !> Non-increasing order
      logical, intent(in) :: descending
      !> Input position of each sorted entry
      integer :: perm(size(keys))

      integer :: i, pos

      do i = 1, size(keys)
         if (descending) then
            pos = 1 + count(keys > keys(i)) + count(keys(:i - 1) == keys(i))
         else
            pos = 1 + count(keys < keys(i)) + count(keys(:i - 1) == keys(i))
         end if
         perm(pos) = i
      end do
   end function stable_reference_real

   !> Check that `idx` holds each of 1..size(idx) once
   !>
   !> @param[out] error Test failure
   !> @param[in]  idx   Index array
   !> @param[in]  what  Label
   subroutine check_permutation(error, idx, what)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Index array
      integer, intent(in) :: idx(:)
      !> Label
      character(len=*), intent(in) :: what

      logical :: seen(size(idx))
      integer :: i

      seen = .false.
      do i = 1, size(idx)
         if (idx(i) < 1 .or. idx(i) > size(idx)) then
            call check(error, .false., what//": index out of range")
            return
         end if
         seen(idx(i)) = .true.
      end do
      call check(error, all(seen), what//": index is not a permutation")
   end subroutine check_permutation

   !> Whether the values are non-decreasing
   !>
   !> @param[in] a Values
   pure function is_sorted_int(a) result(sorted)
      !> Values
      integer, intent(in) :: a(:)
      !> Non-decreasing
      logical :: sorted

      sorted = .true.
      if (size(a) > 1) sorted = all(a(:size(a) - 1) <= a(2:))
   end function is_sorted_int

   !> Whether the values are non-decreasing
   !>
   !> @param[in] a Values
   pure function is_sorted_real(a) result(sorted)
      !> Values
      real(wp), intent(in) :: a(:)
      !> Non-decreasing
      logical :: sorted

      sorted = .true.
      if (size(a) > 1) sorted = all(a(:size(a) - 1) <= a(2:))
   end function is_sorted_real

   !> Whether any two entries are equal
   !>
   !> @param[in] a Values
   pure function has_ties_int(a) result(ties)
      !> Values
      integer, intent(in) :: a(:)
      !> Two equal entries
      logical :: ties

      integer :: i

      ties = .false.
      do i = 2, size(a)
         if (any(a(:i - 1) == a(i))) then
            ties = .true.
            return
         end if
      end do
   end function has_ties_int

   !> Whether any two entries are equal
   !>
   !> @param[in] a Values
   pure function has_ties_real(a) result(ties)
      !> Values
      real(wp), intent(in) :: a(:)
      !> Two equal entries
      logical :: ties

      integer :: i

      ties = .false.
      do i = 2, size(a)
         if (any(a(:i - 1) == a(i))) then
            ties = .true.
            return
         end if
      end do
   end function has_ties_real

   !> Failure label of a fixture
   !>
   !> @param[in]  n       Number of entries
   !> @param[in]  pattern Pattern id, 1..npattern
   !> @param[out] text    Label
   pure subroutine label(n, pattern, text)
      !> Number of entries
      integer, intent(in) :: n
      !> Pattern id
      integer, intent(in) :: pattern
      !> Label
      character(len=:), allocatable, intent(out) :: text

      text = " (n = "//to_string(n)//", "//trim(pattern_name(pattern))//")"
   end subroutine label

end module test_math_sorters
