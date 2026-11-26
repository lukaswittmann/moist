!> Tests of the DROP OpenMP framework: lowest-index abort latch, per-thread worker slots
module test_cavity_drop_threads
   use mctc_env, only: moist_error_type => error_type, fatal_error
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_cavity_drop_threads, only: drop_abort_latch_type, drop_worker_slots_type
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_cavity_drop_threads

contains

   !> Collect the latch and slot tests
   !>
   !> @param[out] testsuite Thread framework tests
   subroutine collect_cavity_drop_threads(testsuite)
      !> Thread framework tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("latch_lowest_index_wins", test_latch_lowest_index), &
         & new_unittest("latch_skips_behind_failure", test_latch_skips), &
         & new_unittest("latch_reset_clears", test_latch_reset), &
         & new_unittest("latch_error_message_status", test_latch_variants), &
         & new_unittest("latch_status_and_error_replace", test_latch_replace), &
         & new_unittest("latch_parallel_loop", test_latch_parallel), &
         & new_unittest("slots_sized_by_count", test_slots_count) &
         & ]
   end subroutine collect_cavity_drop_threads

   !> The lowest latched index wins, whatever the order the failures arrive in
   subroutine test_latch_lowest_index(error)
      type(error_type), allocatable, intent(out) :: error
      type(drop_abort_latch_type) :: latch
      integer, parameter :: orders(3, 3) = reshape([5, 3, 7, 7, 5, 3, 3, 7, 5], [3, 3])
      integer :: i, k

      do k = 1, size(orders, 2)
         call latch%reset()
         do i = 1, size(orders, 1)
            call latch_index(latch, orders(i, k))
         end do
         call check(error, latch%is_set(), "nothing latched")
         if (allocated(error)) return
         call check(error, latch%igrid, 3, "not the lowest index")
         if (allocated(error)) return
         call check(error, allocated(latch%error), "the error was not kept")
         if (allocated(error)) return
         call check(error, latch%error%message == "failure 3", "not the lowest index's error: "//latch%error%message)
         if (allocated(error)) return
      end do
   end subroutine test_latch_lowest_index

   !> Only indices behind the latched failure are skipped; an index-free failure skips all
   subroutine test_latch_skips(error)
      type(error_type), allocatable, intent(out) :: error
      type(drop_abort_latch_type) :: latch

      call latch%reset()
      call check(error, .not. latch%skips(1) .and. .not. latch%skips(huge(1)), "skips without a failure")
      if (allocated(error)) return
      call latch%latch_message("failure 4", 4)
      call check(error, .not. latch%skips(3) .and. .not. latch%skips(4) .and. latch%skips(5), &
                 "skips not behind index 4")
      if (allocated(error)) return
      ! A failure without an index (setup) displaces the indexed one and skips every index
      call latch%latch_message("setup")
      call check(error, latch%igrid, 0, "an index-free failure did not win")
      if (allocated(error)) return
      call check(error, latch%skips(1), "index 1 not skipped behind a setup failure")
      if (allocated(error)) return
      call latch%latch_message("failure 2", 2)
      call check(error, latch%igrid == 0 .and. latch%error%message == "setup", &
                 "an indexed failure displaced the setup failure")
   end subroutine test_latch_skips

   !> A reset latch holds nothing and skips nothing
   subroutine test_latch_reset(error)
      type(error_type), allocatable, intent(out) :: error
      type(drop_abort_latch_type) :: latch

      call latch%latch_message("failure 2", 2)
      call latch%reset()
      call check(error, .not. latch%is_set() .and. latch%igrid == 0 .and. latch%status == 0 &
                 .and. .not. allocated(latch%error), "reset left a failure")
      if (allocated(error)) return
      call check(error, .not. latch%skips(3), "reset latch skips")
      if (allocated(error)) return
      call latch%latch_status(7, 2)
      call latch%reset()
      call check(error, .not. latch%is_set() .and. latch%status == 0, "reset left a status")
   end subroutine test_latch_reset

   !> Each latch call stores its own representation; a losing error object is freed
   subroutine test_latch_variants(error)
      type(error_type), allocatable, intent(out) :: error
      type(drop_abort_latch_type) :: latch
      type(moist_error_type), allocatable :: err

      call latch%reset()
      call fatal_error(err, "object 3")
      call latch%latch_error(err, 3)
      call check(error, .not. allocated(err), "the latched error object was not taken")
      if (allocated(error)) return
      call check(error, latch%igrid == 3 .and. latch%status == 0 .and. latch%error%message == "object 3", &
                 "error object not latched")
      if (allocated(error)) return
      call fatal_error(err, "object 6")
      call latch%latch_error(err, 6)
      call check(error, .not. allocated(err), "the losing error object was not freed")
      if (allocated(error)) return
      call check(error, latch%error%message == "object 3", "a higher index displaced the error")
      if (allocated(error)) return

      call latch%reset()
      call latch%latch_message("message 2", 2)
      call check(error, latch%igrid == 2 .and. latch%status == 0 .and. latch%error%message == "message 2", &
                 "message not latched")
      if (allocated(error)) return

      call latch%reset()
      call latch%latch_status(-4, 9)
      call check(error, latch%igrid == 9 .and. latch%status == -4 .and. .not. allocated(latch%error), &
                 "status not latched")
   end subroutine test_latch_variants

   !> A winning status drops the latched error, a winning error clears the status
   subroutine test_latch_replace(error)
      type(error_type), allocatable, intent(out) :: error
      type(drop_abort_latch_type) :: latch

      call latch%reset()
      call latch%latch_message("message 5", 5)
      call latch%latch_status(7, 2)
      call check(error, latch%igrid == 2 .and. latch%status == 7 .and. .not. allocated(latch%error), &
                 "a status did not replace the error")
      if (allocated(error)) return
      call latch%latch_message("message 1", 1)
      call check(error, latch%igrid == 1 .and. latch%status == 0 .and. allocated(latch%error), &
                 "an error did not replace the status")
      if (allocated(error)) return
      call check(error, latch%error%message == "message 1", "wrong message after the replacement")
      if (allocated(error)) return
      ! Losers replace nothing
      call latch%latch_status(8, 4)
      call check(error, latch%igrid == 1 .and. latch%status == 0 .and. allocated(latch%error), &
                 "a losing status replaced the error")
   end subroutine test_latch_replace

   !> In a parallel loop, indices taken from the top, the lowest failing index wins
   subroutine test_latch_parallel(error)
      type(error_type), allocatable, intent(out) :: error
      !> Loop length, failures at every multiple of `stride`
      integer, parameter :: n = 400, stride = 37
      type(drop_abort_latch_type) :: latch
      integer :: i

      call latch%reset()
      !$omp parallel do num_threads(4) schedule(dynamic) default(shared) private(i)
      do i = n, 1, -1
         if (latch%skips(i)) cycle
         if (modulo(i, stride) == 0) call latch_index(latch, i)
      end do
      !$omp end parallel do
      call check(error, latch%igrid, stride, "not the lowest failing index")
      if (allocated(error)) return
      call check(error, latch%error%message == "failure 37", "not the lowest index's error")
   end subroutine test_latch_parallel

   !> Slots sized for a given team, at least one
   subroutine test_slots_count(error)
      type(error_type), allocatable, intent(out) :: error
      type(drop_worker_slots_type) :: slots

      call slots%init(3)
      call check(error, slots%nthreads == 3 .and. size(slots%lsf) == 3 .and. .not. allocated(slots%phi), &
                 "three slots without objectives expected")
      if (allocated(error)) return
      call check(error, .not. allocated(slots%lsf(1)%lsf), "a slot was filled at init")
      if (allocated(error)) return
      call slots%init(0)
      call check(error, slots%nthreads == 1 .and. size(slots%lsf) == 1, "no team gives one slot")
      if (allocated(error)) return
      call slots%init(2, objective=.true.)
      call check(error, size(slots%lsf) == 2 .and. size(slots%phi) == 2, "objectives not sized with the slots")
   end subroutine test_slots_count

   !> Latch the failure of loop index `i` with the message "failure i"
   !>
   !> @param[inout] latch Abort latch
   !> @param[in]    i     Failing index
   subroutine latch_index(latch, i)
      !> Abort latch
      type(drop_abort_latch_type), intent(inout) :: latch
      !> Failing index
      integer, intent(in) :: i

      call latch%latch_message("failure "//to_string(i), i)
   end subroutine latch_index

end module test_cavity_drop_threads
