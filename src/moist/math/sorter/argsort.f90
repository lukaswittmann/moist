! SPDX-License-Identifier: MPL-2.0 AND MIT
! Ported from fortran-lang/stdlib v0.8.1 (commit 0ede3015), src/sorting/stdlib_sorting_sort_adjoint.fypp
! and `sort_index` of src/sorting/stdlib_sorting.fypp; MIT, Copyright (c) 2019-2021 stdlib
! contributors, license text in LICENSE next to this file. stdlib translates the merge
! sort of Rust's slice.rs (Apache-2.0 OR MIT, Copyright 2012-2015 The Rust Project Developers)
!
!> Stable argsort, stdlib's `sort_index` and `sort_adjoint`
!>
!> - [[ord_sort]]'s natural-run merge sort that moves an integer companion
!>   array along with the keys
!> - `sort_adjoint(array, adjoint_array)`: the caller's companion, e.g. labels
!> - `sort_index(array, index)`: the companion is `1..n`, so `array_in(index)`
!>   is `array_out`; equal keys keep ascending `index`
!> - both sort `array` itself
!> - scratch of `size(array)/2` keys (`work`) and integers (`iwork`), the
!>   caller's or allocated; below 64 entries nothing is allocated
!> - kinds: default-integer and `real(wp)` keys, default-integer companions;
!>   with a NaN the order is undefined
!>
!> Changes from stdlib: kinds instantiated by hand, default-integer indices
!> and companions, `pure`, `contiguous` dummies, companion size checked, no
!> buffers when no merge runs
module moist_math_sorter_argsort
   use mctc_env, only: wp
   use moist_math_sorter_merge, only: run_type, max_merge_stack, calc_min_run, collapse
   implicit none(type, external)
   private

   public :: sort_index, sort_adjoint

   !> Stable argsort: `index` such that `array_in(index)` is `array_out`
   !>
   !> `call sort_index(array, index[, work, iwork, reverse])`
   !>
   !> - `array` is sorted too; `index` needs at least `size(array)` entries
   !> - `work`, `iwork`: scratch, at least `size(array)/2` entries each
   !> - `reverse`: non-increasing order, equal keys keep ascending `index`
   interface sort_index
      module procedure :: sort_index_int
      module procedure :: sort_index_real
   end interface sort_index

   !> Stable sort of `array` that moves `adjoint_array` along
   !>
   !> `call sort_adjoint(array, adjoint_array[, work, iwork, reverse])`
   interface sort_adjoint
      module procedure :: sort_adjoint_int
      module procedure :: sort_adjoint_real
   end interface sort_adjoint

contains

   !* ================================================================================= *!
   !*                                  Default integer                                  *!
   !* ================================================================================= *!

   !> Stable argsort of integer keys
   !>
   !> @param[inout] array   Keys, sorted on return
   !> @param[out]   index   Original position of each sorted key
   !> @param[out]   work    Key scratch, at least `size(array)/2` entries
   !> @param[out]   iwork   Index scratch, at least `size(array)/2` entries
   !> @param[in]    reverse Non-increasing order
   pure subroutine sort_index_int(array, index, work, iwork, reverse)
      !> Keys, sorted on return
      integer, intent(inout), contiguous :: array(:)
      !> Original position of each sorted key
      integer, intent(out), contiguous :: index(:)
      !> Key scratch
      integer, intent(out), optional, contiguous :: work(:)
      !> Index scratch
      integer, intent(out), optional, contiguous :: iwork(:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      integer :: i, n

      n = size(array)
      if (size(index) < n) error stop "sort_index: index is smaller than array"
      do i = 1, n
         index(i) = i
      end do
      call sort_adjoint_int(array, index(1:n), work, iwork, reverse)
   end subroutine sort_index_int

   !> Stable sort of integer keys moving a companion array along
   !>
   !> @param[inout] array         Keys, sorted on return
   !> @param[inout] adjoint_array Companion, permuted with `array`
   !> @param[out]   work          Key scratch, at least `size(array)/2` entries
   !> @param[out]   iwork         Companion scratch, at least `size(array)/2` entries
   !> @param[in]    reverse       Non-increasing order
   pure subroutine sort_adjoint_int(array, adjoint_array, work, iwork, reverse)
      !> Keys, sorted on return
      integer, intent(inout), contiguous :: array(0:)
      !> Companion, permuted with `array`
      integer, intent(inout), contiguous :: adjoint_array(0:)
      !> Key scratch
      integer, intent(out), optional, contiguous :: work(0:)
      !> Companion scratch
      integer, intent(out), optional, contiguous :: iwork(0:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      integer, allocatable :: scratch(:), iscratch(:)
      integer :: n, stat
      logical :: rev

      n = size(array)
      if (size(adjoint_array) < n) error stop "sort_adjoint: adjoint_array is smaller than array"
      if (present(work)) then
         if (size(work) < n/2) error stop "sort_adjoint: work array is too small"
      end if
      if (present(iwork)) then
         if (size(iwork) < n/2) error stop "sort_adjoint: iwork array is too small"
      end if
      rev = .false.
      if (present(reverse)) rev = reverse

      if (rev) call reverse_int(array, adjoint_array)
      if (n <= calc_min_run(n)) then
         if (n >= 2) call insertion_int(array, adjoint_array)
      else if (present(work) .and. present(iwork)) then
         call merge_sort_int(array, adjoint_array, work, iwork)
      else if (present(work)) then
         allocate (iscratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the companion buffer failed"
         call merge_sort_int(array, adjoint_array, work, iscratch)
      else if (present(iwork)) then
         allocate (scratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the key buffer failed"
         call merge_sort_int(array, adjoint_array, scratch, iwork)
      else
         allocate (scratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the key buffer failed"
         allocate (iscratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the companion buffer failed"
         call merge_sort_int(array, adjoint_array, scratch, iscratch)
      end if
      if (rev) call reverse_int(array, adjoint_array)
   end subroutine sort_adjoint_int

   !> Natural-run merge sort driven by the run stack, companion moved along
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   !> @param[inout] buf           Key scratch, at least `size(array)/2` entries
   !> @param[inout] ibuf          Companion scratch, at least `size(array)/2` entries
   pure subroutine merge_sort_int(array, adjoint_array, buf, ibuf)
      !> Keys
      integer, intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)
      !> Key scratch
      integer, intent(inout), contiguous :: buf(0:)
      !> Companion scratch
      integer, intent(inout), contiguous :: ibuf(0:)

      integer :: array_size, finish, min_run, r, r_count, start
      type(run_type) :: runs(0:max_merge_stack - 1), left, right

      array_size = size(array)
      min_run = calc_min_run(array_size)
      if (array_size <= min_run) then
         if (array_size >= 2) call insertion_int(array, adjoint_array)
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
               call reverse_int(array(start:finish), adjoint_array(start:finish))
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
            call insert_head_int(array(start:finish), adjoint_array(start:finish))
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
            call merge_int(array(left%base:right%base + right%len - 1), left%len, buf, &
                           adjoint_array(left%base:right%base + right%len - 1), ibuf)
            runs(r) = run_type(base=left%base, len=left%len + right%len)
            if (r == r_count - 3) runs(r + 1) = runs(r + 2)
            r_count = r_count - 1
         end do merge_loop
      end do
      if (r_count /= 1) error stop "sort_adjoint: merge sort completed without one run"
   end subroutine merge_sort_int

   !> Merge the sorted runs `array(0:mid-1)` and `array(mid:)` in place
   !>
   !> The shorter run is copied to the scratch; ties take the left run first
   !>
   !> @param[inout] array         Two adjacent sorted runs
   !> @param[in]    mid           Length of the left run, `0 < mid < size(array)`
   !> @param[inout] buf           Key scratch, at least the shorter run
   !> @param[inout] adjoint_array Companion of `array`
   !> @param[inout] ibuf          Companion scratch, at least the shorter run
   pure subroutine merge_int(array, mid, buf, adjoint_array, ibuf)
      !> Two adjacent sorted runs
      integer, intent(inout), contiguous :: array(0:)
      !> Length of the left run
      integer, intent(in) :: mid
      !> Key scratch
      integer, intent(inout), contiguous :: buf(0:)
      !> Companion of `array`
      integer, intent(inout), contiguous :: adjoint_array(0:)
      !> Companion scratch
      integer, intent(inout), contiguous :: ibuf(0:)

      integer :: array_len, i, j, k

      array_len = size(array)
      if (mid <= array_len - mid) then
         ! Left run shorter: merge forwards
         buf(0:mid - 1) = array(0:mid - 1)
         ibuf(0:mid - 1) = adjoint_array(0:mid - 1)
         i = 0
         j = mid
         merge_lower: do k = 0, array_len - 1
            if (buf(i) <= array(j)) then
               array(k) = buf(i)
               adjoint_array(k) = ibuf(i)
               i = i + 1
               if (i >= mid) exit merge_lower
            else
               array(k) = array(j)
               adjoint_array(k) = adjoint_array(j)
               j = j + 1
               if (j >= array_len) then
                  array(k + 1:) = buf(i:mid - 1)
                  adjoint_array(k + 1:) = ibuf(i:mid - 1)
                  exit merge_lower
               end if
            end if
         end do merge_lower
      else
         ! Right run shorter: merge backwards
         buf(0:array_len - mid - 1) = array(mid:array_len - 1)
         ibuf(0:array_len - mid - 1) = adjoint_array(mid:array_len - 1)
         i = mid - 1
         j = array_len - mid - 1
         merge_upper: do k = array_len - 1, 0, -1
            if (buf(j) >= array(i)) then
               array(k) = buf(j)
               adjoint_array(k) = ibuf(j)
               j = j - 1
               if (j < 0) exit merge_upper
            else
               array(k) = array(i)
               adjoint_array(k) = adjoint_array(i)
               i = i - 1
               if (i < 0) then
                  array(0:k - 1) = buf(0:j)
                  adjoint_array(0:k - 1) = ibuf(0:j)
                  exit merge_upper
               end if
            end if
         end do merge_upper
      end if
   end subroutine merge_int

   !> Stable insertion sort, companion moved along
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   pure subroutine insertion_int(array, adjoint_array)
      !> Keys
      integer, intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)

      integer :: i, j, key, key_adjoint

      do j = 1, size(array) - 1
         key = array(j)
         key_adjoint = adjoint_array(j)
         i = j - 1
         do while (i >= 0)
            if (array(i) <= key) exit
            array(i + 1) = array(i)
            adjoint_array(i + 1) = adjoint_array(i)
            i = i - 1
         end do
         array(i + 1) = key
         adjoint_array(i + 1) = key_adjoint
      end do
   end subroutine insertion_int

   !> Insert `array(0)` into the sorted `array(1:)`, ahead of equal keys
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   pure subroutine insert_head_int(array, adjoint_array)
      !> Keys
      integer, intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)

      integer :: i, tmp, tmp_adjoint

      tmp = array(0)
      tmp_adjoint = adjoint_array(0)
      find_hole: do i = 1, size(array) - 1
         if (array(i) >= tmp) exit find_hole
         array(i - 1) = array(i)
         adjoint_array(i - 1) = adjoint_array(i)
      end do find_hole
      array(i - 1) = tmp
      adjoint_array(i - 1) = tmp_adjoint
   end subroutine insert_head_int

   !> Reverse keys and companion in place
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   pure subroutine reverse_int(array, adjoint_array)
      !> Keys
      integer, intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)

      integer :: lo, hi, tmp, itmp

      lo = 0
      hi = size(array) - 1
      do while (lo < hi)
         tmp = array(lo)
         array(lo) = array(hi)
         array(hi) = tmp
         itmp = adjoint_array(lo)
         adjoint_array(lo) = adjoint_array(hi)
         adjoint_array(hi) = itmp
         lo = lo + 1
         hi = hi - 1
      end do
   end subroutine reverse_int

   !* ================================================================================= *!
   !*                                      real(wp)                                     *!
   !* ================================================================================= *!

   !> Stable argsort of real keys
   !>
   !> @param[inout] array   Keys, sorted on return
   !> @param[out]   index   Original position of each sorted key
   !> @param[out]   work    Key scratch, at least `size(array)/2` entries
   !> @param[out]   iwork   Index scratch, at least `size(array)/2` entries
   !> @param[in]    reverse Non-increasing order
   pure subroutine sort_index_real(array, index, work, iwork, reverse)
      !> Keys, sorted on return
      real(wp), intent(inout), contiguous :: array(:)
      !> Original position of each sorted key
      integer, intent(out), contiguous :: index(:)
      !> Key scratch
      real(wp), intent(out), optional, contiguous :: work(:)
      !> Index scratch
      integer, intent(out), optional, contiguous :: iwork(:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      integer :: i, n

      n = size(array)
      if (size(index) < n) error stop "sort_index: index is smaller than array"
      do i = 1, n
         index(i) = i
      end do
      call sort_adjoint_real(array, index(1:n), work, iwork, reverse)
   end subroutine sort_index_real

   !> Stable sort of real keys moving a companion array along
   !>
   !> @param[inout] array         Keys, sorted on return
   !> @param[inout] adjoint_array Companion, permuted with `array`
   !> @param[out]   work          Key scratch, at least `size(array)/2` entries
   !> @param[out]   iwork         Companion scratch, at least `size(array)/2` entries
   !> @param[in]    reverse       Non-increasing order
   pure subroutine sort_adjoint_real(array, adjoint_array, work, iwork, reverse)
      !> Keys, sorted on return
      real(wp), intent(inout), contiguous :: array(0:)
      !> Companion, permuted with `array`
      integer, intent(inout), contiguous :: adjoint_array(0:)
      !> Key scratch
      real(wp), intent(out), optional, contiguous :: work(0:)
      !> Companion scratch
      integer, intent(out), optional, contiguous :: iwork(0:)
      !> Non-increasing order
      logical, intent(in), optional :: reverse

      real(wp), allocatable :: scratch(:)
      integer, allocatable :: iscratch(:)
      integer :: n, stat
      logical :: rev

      n = size(array)
      if (size(adjoint_array) < n) error stop "sort_adjoint: adjoint_array is smaller than array"
      if (present(work)) then
         if (size(work) < n/2) error stop "sort_adjoint: work array is too small"
      end if
      if (present(iwork)) then
         if (size(iwork) < n/2) error stop "sort_adjoint: iwork array is too small"
      end if
      rev = .false.
      if (present(reverse)) rev = reverse

      if (rev) call reverse_real(array, adjoint_array)
      if (n <= calc_min_run(n)) then
         if (n >= 2) call insertion_real(array, adjoint_array)
      else if (present(work) .and. present(iwork)) then
         call merge_sort_real(array, adjoint_array, work, iwork)
      else if (present(work)) then
         allocate (iscratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the companion buffer failed"
         call merge_sort_real(array, adjoint_array, work, iscratch)
      else if (present(iwork)) then
         allocate (scratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the key buffer failed"
         call merge_sort_real(array, adjoint_array, scratch, iwork)
      else
         allocate (scratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the key buffer failed"
         allocate (iscratch(0:n/2 - 1), stat=stat)
         if (stat /= 0) error stop "sort_adjoint: allocation of the companion buffer failed"
         call merge_sort_real(array, adjoint_array, scratch, iscratch)
      end if
      if (rev) call reverse_real(array, adjoint_array)
   end subroutine sort_adjoint_real

   !> Natural-run merge sort driven by the run stack, companion moved along
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   !> @param[inout] buf           Key scratch, at least `size(array)/2` entries
   !> @param[inout] ibuf          Companion scratch, at least `size(array)/2` entries
   pure subroutine merge_sort_real(array, adjoint_array, buf, ibuf)
      !> Keys
      real(wp), intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)
      !> Key scratch
      real(wp), intent(inout), contiguous :: buf(0:)
      !> Companion scratch
      integer, intent(inout), contiguous :: ibuf(0:)

      integer :: array_size, finish, min_run, r, r_count, start
      type(run_type) :: runs(0:max_merge_stack - 1), left, right

      array_size = size(array)
      min_run = calc_min_run(array_size)
      if (array_size <= min_run) then
         if (array_size >= 2) call insertion_real(array, adjoint_array)
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
               call reverse_real(array(start:finish), adjoint_array(start:finish))
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
            call insert_head_real(array(start:finish), adjoint_array(start:finish))
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
            call merge_real(array(left%base:right%base + right%len - 1), left%len, buf, &
                            adjoint_array(left%base:right%base + right%len - 1), ibuf)
            runs(r) = run_type(base=left%base, len=left%len + right%len)
            if (r == r_count - 3) runs(r + 1) = runs(r + 2)
            r_count = r_count - 1
         end do merge_loop
      end do
      if (r_count /= 1) error stop "sort_adjoint: merge sort completed without one run"
   end subroutine merge_sort_real

   !> Merge the sorted runs `array(0:mid-1)` and `array(mid:)` in place
   !>
   !> The shorter run is copied to the scratch; ties take the left run first
   !>
   !> @param[inout] array         Two adjacent sorted runs
   !> @param[in]    mid           Length of the left run, `0 < mid < size(array)`
   !> @param[inout] buf           Key scratch, at least the shorter run
   !> @param[inout] adjoint_array Companion of `array`
   !> @param[inout] ibuf          Companion scratch, at least the shorter run
   pure subroutine merge_real(array, mid, buf, adjoint_array, ibuf)
      !> Two adjacent sorted runs
      real(wp), intent(inout), contiguous :: array(0:)
      !> Length of the left run
      integer, intent(in) :: mid
      !> Key scratch
      real(wp), intent(inout), contiguous :: buf(0:)
      !> Companion of `array`
      integer, intent(inout), contiguous :: adjoint_array(0:)
      !> Companion scratch
      integer, intent(inout), contiguous :: ibuf(0:)

      integer :: array_len, i, j, k

      array_len = size(array)
      if (mid <= array_len - mid) then
         ! Left run shorter: merge forwards
         buf(0:mid - 1) = array(0:mid - 1)
         ibuf(0:mid - 1) = adjoint_array(0:mid - 1)
         i = 0
         j = mid
         merge_lower: do k = 0, array_len - 1
            if (buf(i) <= array(j)) then
               array(k) = buf(i)
               adjoint_array(k) = ibuf(i)
               i = i + 1
               if (i >= mid) exit merge_lower
            else
               array(k) = array(j)
               adjoint_array(k) = adjoint_array(j)
               j = j + 1
               if (j >= array_len) then
                  array(k + 1:) = buf(i:mid - 1)
                  adjoint_array(k + 1:) = ibuf(i:mid - 1)
                  exit merge_lower
               end if
            end if
         end do merge_lower
      else
         ! Right run shorter: merge backwards
         buf(0:array_len - mid - 1) = array(mid:array_len - 1)
         ibuf(0:array_len - mid - 1) = adjoint_array(mid:array_len - 1)
         i = mid - 1
         j = array_len - mid - 1
         merge_upper: do k = array_len - 1, 0, -1
            if (buf(j) >= array(i)) then
               array(k) = buf(j)
               adjoint_array(k) = ibuf(j)
               j = j - 1
               if (j < 0) exit merge_upper
            else
               array(k) = array(i)
               adjoint_array(k) = adjoint_array(i)
               i = i - 1
               if (i < 0) then
                  array(0:k - 1) = buf(0:j)
                  adjoint_array(0:k - 1) = ibuf(0:j)
                  exit merge_upper
               end if
            end if
         end do merge_upper
      end if
   end subroutine merge_real

   !> Stable insertion sort, companion moved along
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   pure subroutine insertion_real(array, adjoint_array)
      !> Keys
      real(wp), intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)

      integer :: i, j, key_adjoint
      real(wp) :: key

      do j = 1, size(array) - 1
         key = array(j)
         key_adjoint = adjoint_array(j)
         i = j - 1
         do while (i >= 0)
            if (array(i) <= key) exit
            array(i + 1) = array(i)
            adjoint_array(i + 1) = adjoint_array(i)
            i = i - 1
         end do
         array(i + 1) = key
         adjoint_array(i + 1) = key_adjoint
      end do
   end subroutine insertion_real

   !> Insert `array(0)` into the sorted `array(1:)`, ahead of equal keys
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   pure subroutine insert_head_real(array, adjoint_array)
      !> Keys
      real(wp), intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)

      integer :: i, tmp_adjoint
      real(wp) :: tmp

      tmp = array(0)
      tmp_adjoint = adjoint_array(0)
      find_hole: do i = 1, size(array) - 1
         if (array(i) >= tmp) exit find_hole
         array(i - 1) = array(i)
         adjoint_array(i - 1) = adjoint_array(i)
      end do find_hole
      array(i - 1) = tmp
      adjoint_array(i - 1) = tmp_adjoint
   end subroutine insert_head_real

   !> Reverse keys and companion in place
   !>
   !> @param[inout] array         Keys
   !> @param[inout] adjoint_array Companion
   pure subroutine reverse_real(array, adjoint_array)
      !> Keys
      real(wp), intent(inout), contiguous :: array(0:)
      !> Companion
      integer, intent(inout), contiguous :: adjoint_array(0:)

      integer :: lo, hi, itmp
      real(wp) :: tmp

      lo = 0
      hi = size(array) - 1
      do while (lo < hi)
         tmp = array(lo)
         array(lo) = array(hi)
         array(hi) = tmp
         itmp = adjoint_array(lo)
         adjoint_array(lo) = adjoint_array(hi)
         adjoint_array(hi) = itmp
         lo = lo + 1
         hi = hi - 1
      end do
   end subroutine reverse_real

end module moist_math_sorter_argsort
