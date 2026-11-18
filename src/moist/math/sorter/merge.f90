! SPDX-License-Identifier: MPL-2.0 AND MIT
! Ported from fortran-lang/stdlib v0.8.1 (commit 0ede3015), src/sorting/stdlib_sorting_ord_sort.fypp
! and the run stack of src/sorting/stdlib_sorting.fypp; MIT, Copyright (c) 2019-2021 stdlib
! contributors, license text in LICENSE next to this file. stdlib translates the merge
! sort of Rust's slice.rs (Apache-2.0 OR MIT, Copyright 2012-2015 The Rust Project Developers)
!
!> Stable natural-run merge sort, stdlib's `ord_sort`
!>
!> - Rust `slice.rs` merge sort, a simplified Timsort: natural runs, strictly
!>   descending ones reversed, are extended to a minimum run by insertion sort
!>   and merged under the Timsort stack invariants
!> - O(n) on ascending or descending input, O(n log n) worst case
!> - scratch of `size(array)/2` entries: the caller's `work`, else allocated;
!>   below 64 entries the whole array is insertion-sorted and nothing is allocated
!> - kinds: default integer and `real(wp)`; with a NaN the order is undefined
!>
!> Changes from stdlib: kinds instantiated by hand, default-integer indices,
!> `pure`, `contiguous` dummies, `reverse` as a reversal around the ascending
!> sort (stdlib's `sort_index` scheme), no buffer when no merge runs
module moist_math_sorter_merge
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: ord_sort
   !> Shared with [[moist_math_sorter_argsort]]
   public :: run_type, max_merge_stack, calc_min_run, collapse

   !> Most pending runs, enough for 2**64 elements (Fibonacci growth of the run lengths)
   integer, parameter :: max_merge_stack = int(ceiling(log(2.0_wp**64)/log(1.6180339887_wp)))

   !> Pending run of the merge stack
   type :: run_type
      !> First entry, 0-based
      integer :: base = 0
      !> Number of entries
      integer :: len = 0
   end type run_type

   !> Stable sort in place, non-decreasing or with `reverse` non-increasing
   !>
   !> `call ord_sort(array[, work, reverse])`
   !>
   !> - `work`: scratch of the array's type, at least `size(array)/2`
   !>   entries, undefined on return
   !> - `reverse`: non-increasing order, equal values keep their input order
   interface ord_sort
      module procedure :: ord_sort_int
      module procedure :: ord_sort_real
   end interface ord_sort

contains

   !* ================================================================================= *!
   !*                                  Shared run logic                                 *!
   !* ================================================================================= *!

   !> Minimum run length: 32..63, so that `n/min_run` is at most a power of two;
   !> `n` itself below 64
   !>
   !> @param[in] n Number of entries
   !> @returns     Minimum run length
   pure function calc_min_run(n) result(min_run)
      !> Number of entries
      integer, intent(in) :: n
      !> Minimum run length
      integer :: min_run

      integer :: num, r

      num = n
      r = 0
      do while (num >= 64)
         r = ior(r, iand(num, 1))
         num = ishft(num, -1)
      end do
      min_run = num + r
   end function calc_min_run

   !> Pair of adjacent runs to merge next, -1 when the stack invariants hold
   !>
   !> Invariants, `len(-1)` the newest run:
   !>
   !> 1. `len(-3) > len(-2) + len(-1)`
   !> 2. `len(-2) > len(-1)`
   !>
   !> @param[in] runs Pending runs, oldest first
   !> @returns        Index of the older run of the pair, or -1
   pure function collapse(runs) result(r)
      !> Pending runs, oldest first
      type(run_type), intent(in) :: runs(0:)
      !> Index of the older run of the pair, or -1
      integer :: r

      integer :: n
      logical :: test

      n = size(runs)
      test = .false.
      if (n >= 2) then
         if (runs(n - 1)%base == 0 .or. runs(n - 2)%len <= runs(n - 1)%len) then
            test = .true.
         else if (n >= 3) then
            ! |X| <= |Y| + |Z|: merge for rho1 or rho2
            if (runs(n - 3)%len <= runs(n - 2)%len + runs(n - 1)%len) then
               test = .true.
            else if (n >= 4) then
               ! |W| <= |X| + |Y|: merge for rho1 or rho3
               if (runs(n - 4)%len <= runs(n - 3)%len + runs(n - 2)%len) test = .true.
            end if
         end if
      end if
      if (test) then
         ! |X| < |Z|: merge X and Y (rho1), else Y and Z (rho2..rho4)
         if (n >= 3) then
            if (runs(n - 3)%len < runs(n - 1)%len) then
               r = n - 3
               return
            end if
         end if
         r = n - 2
      else
         r = -1
      end if
   end function collapse

   !* ================================================================================= *!
   !*                                   Default integer                                 *!
   !* ================================================================================= *!

   !> Stable sort of an integer array
   !>
   !> @param[inout] array   Values, sorted on return
   !> @param[out]   work    Scratch, at least `size(array)/2` entries
   !> @param[in]    reverse Non-increasing order
   pure subroutine ord_sort_int(array, work, reverse)
      !> Values, sorted on return
      integer, intent(inout), contiguous :: array(0:)
      !> Scratch, at least `size(array)/2` entries
      integer, intent(out), optional, contiguous :: work(0:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      integer, allocatable :: scratch(:)
      integer :: n, stat
      logical :: rev

      n = size(array)
      if (present(work)) then
         if (size(work) < n/2) error stop "ord_sort: work array is too small"
      end if
      rev = .false.
      if (present(reverse)) rev = reverse

      if (rev) call reverse_int(array)
      if (n <= calc_min_run(n)) then
         if (n >= 2) call insertion_int(array)
      else if (present(work)) then
         call merge_sort_int(array, work)
      else
         allocate (scratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "ord_sort: allocation of the buffer failed"
         call merge_sort_int(array, scratch)
      end if
      if (rev) call reverse_int(array)
   end subroutine ord_sort_int

   !> Natural-run merge sort driven by the run stack
   !>
   !> Runs are found walking backwards, so most merges go forwards
   !>
   !> @param[inout] array Values
   !> @param[inout] buf   Scratch, at least `size(array)/2` entries
   pure subroutine merge_sort_int(array, buf)
      !> Values
      integer, intent(inout), contiguous :: array(0:)
      !> Scratch
      integer, intent(inout), contiguous :: buf(0:)

      integer :: array_size, finish, min_run, r, r_count, start
      type(run_type) :: runs(0:max_merge_stack - 1), left, right

      array_size = size(array)
      min_run = calc_min_run(array_size)
      if (array_size <= min_run) then
         if (array_size >= 2) call insertion_int(array)
         return
      end if

      r_count = 0
      finish = array_size - 1
      do while (finish >= 0)
         ! Next natural run, reversed when strictly descending
         start = finish
         if (start > 0) then
            start = start - 1
            if (array(start + 1) < array(start)) then
               descending: do while (start > 0)
                  if (array(start) >= array(start - 1)) exit descending
                  start = start - 1
               end do descending
               call reverse_int(array(start:finish))
            else
               ascending: do while (start > 0)
                  if (array(start) < array(start - 1)) exit ascending
                  start = start - 1
               end do ascending
            end if
         end if

         ! Short run extended to `min_run` by insertion
         insert: do while (start > 0)
            if (finish - start >= min_run - 1) exit insert
            start = start - 1
            call insert_head_int(array(start:finish))
         end do insert
         if (start == 0 .and. finish == array_size - 1) return

         runs(r_count) = run_type(base=start, len=finish - start + 1)
         finish = start - 1
         r_count = r_count + 1

         merge_loop: do
            r = collapse(runs(0:r_count - 1))
            if (r < 0 .or. r_count <= 1) exit merge_loop
            left = runs(r + 1)
            right = runs(r)
            call merge_int(array(left%base:right%base + right%len - 1), left%len, buf)
            runs(r) = run_type(base=left%base, len=left%len + right%len)
            if (r == r_count - 3) runs(r + 1) = runs(r + 2)
            r_count = r_count - 1
         end do merge_loop
      end do
      if (r_count /= 1) error stop "ord_sort: merge sort completed without one run"
   end subroutine merge_sort_int

   !> Merge the sorted runs `array(0:mid-1)` and `array(mid:)` in place
   !>
   !> The shorter run is copied to `buf`; ties take the left run first
   !>
   !> @param[inout] array Two adjacent sorted runs
   !> @param[in]    mid   Length of the left run, `0 < mid < size(array)`
   !> @param[inout] buf   Scratch, at least the shorter run
   pure subroutine merge_int(array, mid, buf)
      !> Two adjacent sorted runs
      integer, intent(inout), contiguous :: array(0:)
      !> Length of the left run
      integer, intent(in) :: mid
      !> Scratch
      integer, intent(inout), contiguous :: buf(0:)

      integer :: array_len, i, j, k

      array_len = size(array)
      if (mid <= array_len - mid) then
         ! Left run shorter: merge forwards
         buf(0:mid - 1) = array(0:mid - 1)
         i = 0
         j = mid
         merge_lower: do k = 0, array_len - 1
            if (buf(i) <= array(j)) then
               array(k) = buf(i)
               i = i + 1
               if (i >= mid) exit merge_lower
            else
               array(k) = array(j)
               j = j + 1
               if (j >= array_len) then
                  array(k + 1:) = buf(i:mid - 1)
                  exit merge_lower
               end if
            end if
         end do merge_lower
      else
         ! Right run shorter: merge backwards
         buf(0:array_len - mid - 1) = array(mid:array_len - 1)
         i = mid - 1
         j = array_len - mid - 1
         merge_upper: do k = array_len - 1, 0, -1
            if (buf(j) >= array(i)) then
               array(k) = buf(j)
               j = j - 1
               if (j < 0) exit merge_upper
            else
               array(k) = array(i)
               i = i - 1
               if (i < 0) then
                  array(0:k - 1) = buf(0:j)
                  exit merge_upper
               end if
            end if
         end do merge_upper
      end if
   end subroutine merge_int

   !> Stable insertion sort
   !>
   !> @param[inout] array Values
   pure subroutine insertion_int(array)
      !> Values
      integer, intent(inout), contiguous :: array(0:)

      integer :: i, j, key

      do j = 1, size(array) - 1
         key = array(j)
         i = j - 1
         do while (i >= 0)
            if (array(i) <= key) exit
            array(i + 1) = array(i)
            i = i - 1
         end do
         array(i + 1) = key
      end do
   end subroutine insertion_int

   !> Insert `array(0)` into the sorted `array(1:)`, ahead of equal values
   !>
   !> @param[inout] array Values
   pure subroutine insert_head_int(array)
      !> Values
      integer, intent(inout), contiguous :: array(0:)

      integer :: i, tmp

      tmp = array(0)
      find_hole: do i = 1, size(array) - 1
         if (array(i) >= tmp) exit find_hole
         array(i - 1) = array(i)
      end do find_hole
      array(i - 1) = tmp
   end subroutine insert_head_int

   !> Reverse in place
   !>
   !> @param[inout] array Values
   pure subroutine reverse_int(array)
      !> Values
      integer, intent(inout), contiguous :: array(0:)

      integer :: lo, hi, tmp

      lo = 0
      hi = size(array) - 1
      do while (lo < hi)
         tmp = array(lo)
         array(lo) = array(hi)
         array(hi) = tmp
         lo = lo + 1
         hi = hi - 1
      end do
   end subroutine reverse_int

   !* ================================================================================= *!
   !*                                      real(wp)                                     *!
   !* ================================================================================= *!

   !> Stable sort of a real array
   !>
   !> @param[inout] array   Values, sorted on return
   !> @param[out]   work    Scratch, at least `size(array)/2` entries
   !> @param[in]    reverse Non-increasing order
   pure subroutine ord_sort_real(array, work, reverse)
      !> Values, sorted on return
      real(wp), intent(inout), contiguous :: array(0:)
      !> Scratch, at least `size(array)/2` entries
      real(wp), intent(out), optional, contiguous :: work(0:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      real(wp), allocatable :: scratch(:)
      integer :: n, stat
      logical :: rev

      n = size(array)
      if (present(work)) then
         if (size(work) < n/2) error stop "ord_sort: work array is too small"
      end if
      rev = .false.
      if (present(reverse)) rev = reverse

      if (rev) call reverse_real(array)
      if (n <= calc_min_run(n)) then
         if (n >= 2) call insertion_real(array)
      else if (present(work)) then
         call merge_sort_real(array, work)
      else
         allocate (scratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "ord_sort: allocation of the buffer failed"
         call merge_sort_real(array, scratch)
      end if
      if (rev) call reverse_real(array)
   end subroutine ord_sort_real

   !> Natural-run merge sort driven by the run stack
   !>
   !> Runs are found walking backwards, so most merges go forwards
   !>
   !> @param[inout] array Values
   !> @param[inout] buf   Scratch, at least `size(array)/2` entries
   pure subroutine merge_sort_real(array, buf)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)
      !> Scratch
      real(wp), intent(inout), contiguous :: buf(0:)

      integer :: array_size, finish, min_run, r, r_count, start
      type(run_type) :: runs(0:max_merge_stack - 1), left, right

      array_size = size(array)
      min_run = calc_min_run(array_size)
      if (array_size <= min_run) then
         if (array_size >= 2) call insertion_real(array)
         return
      end if

      r_count = 0
      finish = array_size - 1
      do while (finish >= 0)
         ! Next natural run, reversed when strictly descending
         start = finish
         if (start > 0) then
            start = start - 1
            if (array(start + 1) < array(start)) then
               descending: do while (start > 0)
                  if (array(start) >= array(start - 1)) exit descending
                  start = start - 1
               end do descending
               call reverse_real(array(start:finish))
            else
               ascending: do while (start > 0)
                  if (array(start) < array(start - 1)) exit ascending
                  start = start - 1
               end do ascending
            end if
         end if

         ! Short run extended to `min_run` by insertion
         insert: do while (start > 0)
            if (finish - start >= min_run - 1) exit insert
            start = start - 1
            call insert_head_real(array(start:finish))
         end do insert
         if (start == 0 .and. finish == array_size - 1) return

         runs(r_count) = run_type(base=start, len=finish - start + 1)
         finish = start - 1
         r_count = r_count + 1

         merge_loop: do
            r = collapse(runs(0:r_count - 1))
            if (r < 0 .or. r_count <= 1) exit merge_loop
            left = runs(r + 1)
            right = runs(r)
            call merge_real(array(left%base:right%base + right%len - 1), left%len, buf)
            runs(r) = run_type(base=left%base, len=left%len + right%len)
            if (r == r_count - 3) runs(r + 1) = runs(r + 2)
            r_count = r_count - 1
         end do merge_loop
      end do
      if (r_count /= 1) error stop "ord_sort: merge sort completed without one run"
   end subroutine merge_sort_real

   !> Merge the sorted runs `array(0:mid-1)` and `array(mid:)` in place
   !>
   !> The shorter run is copied to `buf`; ties take the left run first
   !>
   !> @param[inout] array Two adjacent sorted runs
   !> @param[in]    mid   Length of the left run, `0 < mid < size(array)`
   !> @param[inout] buf   Scratch, at least the shorter run
   pure subroutine merge_real(array, mid, buf)
      !> Two adjacent sorted runs
      real(wp), intent(inout), contiguous :: array(0:)
      !> Length of the left run
      integer, intent(in) :: mid
      !> Scratch
      real(wp), intent(inout), contiguous :: buf(0:)

      integer :: array_len, i, j, k

      array_len = size(array)
      if (mid <= array_len - mid) then
         ! Left run shorter: merge forwards
         buf(0:mid - 1) = array(0:mid - 1)
         i = 0
         j = mid
         merge_lower: do k = 0, array_len - 1
            if (buf(i) <= array(j)) then
               array(k) = buf(i)
               i = i + 1
               if (i >= mid) exit merge_lower
            else
               array(k) = array(j)
               j = j + 1
               if (j >= array_len) then
                  array(k + 1:) = buf(i:mid - 1)
                  exit merge_lower
               end if
            end if
         end do merge_lower
      else
         ! Right run shorter: merge backwards
         buf(0:array_len - mid - 1) = array(mid:array_len - 1)
         i = mid - 1
         j = array_len - mid - 1
         merge_upper: do k = array_len - 1, 0, -1
            if (buf(j) >= array(i)) then
               array(k) = buf(j)
               j = j - 1
               if (j < 0) exit merge_upper
            else
               array(k) = array(i)
               i = i - 1
               if (i < 0) then
                  array(0:k - 1) = buf(0:j)
                  exit merge_upper
               end if
            end if
         end do merge_upper
      end if
   end subroutine merge_real

   !> Stable insertion sort
   !>
   !> @param[inout] array Values
   pure subroutine insertion_real(array)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)

      integer :: i, j
      real(wp) :: key

      do j = 1, size(array) - 1
         key = array(j)
         i = j - 1
         do while (i >= 0)
            if (array(i) <= key) exit
            array(i + 1) = array(i)
            i = i - 1
         end do
         array(i + 1) = key
      end do
   end subroutine insertion_real

   !> Insert `array(0)` into the sorted `array(1:)`, ahead of equal values
   !>
   !> @param[inout] array Values
   pure subroutine insert_head_real(array)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)

      integer :: i
      real(wp) :: tmp

      tmp = array(0)
      find_hole: do i = 1, size(array) - 1
         if (array(i) >= tmp) exit find_hole
         array(i - 1) = array(i)
      end do find_hole
      array(i - 1) = tmp
   end subroutine insert_head_real

   !> Reverse in place
   !>
   !> @param[inout] array Values
   pure subroutine reverse_real(array)
      !> Values
      real(wp), intent(inout), contiguous :: array(0:)

      integer :: lo, hi
      real(wp) :: tmp

      lo = 0
      hi = size(array) - 1
      do while (lo < hi)
         tmp = array(lo)
         array(lo) = array(hi)
         array(hi) = tmp
         lo = lo + 1
         hi = hi - 1
      end do
   end subroutine reverse_real

end module moist_math_sorter_merge
