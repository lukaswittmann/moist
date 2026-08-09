!> Test suite for the shared run context in moist_context
module test_utils_context
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_context, only: moist_context_type, new_context
   implicit none(type, external)
   private

   public :: collect_utils_context

contains

   !> Collect all context tests
   subroutine collect_utils_context(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("defaults-and-settings", test_defaults_and_settings), &
                  new_unittest("shared-timer-tree", test_shared_timer_tree), &
                  new_unittest("writes-guard", test_writes_guard), &
                  new_unittest("profile-flag", test_profile_flag), &
                  new_unittest("report-depth", test_report_depth), &
                  new_unittest("threads-default", test_threads_default), &
                  new_unittest("threads-explicit", test_threads_explicit), &
                  new_unittest("set-num-threads", test_set_num_threads), &
                  new_unittest("delete-releases-pin", test_delete_releases_pin), &
                  new_unittest("owned-logfile", test_owned_logfile), &
                  new_unittest("print-settings-runs", test_print_settings_runs), &
                  new_unittest("debug-message-gated", test_debug_message_gated), &
                  new_unittest("delete-is-safe", test_delete_is_safe) &
                  ]

   end subroutine collect_utils_context

   !> new_context applies defaults and honours explicit settings
   subroutine test_defaults_and_settings(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx)
      call check(error, ctx%verbosity == 1, "default verbosity is 1")
      if (allocated(error)) return
      call check(error,.not. ctx%debug, "default debug is false")
      if (allocated(error)) return
      ! a fresh context owns an initialized, empty timer tree
      call check(error, ctx%timer%num_nodes() == 0, "fresh timer has no nodes")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, verbosity=3, debug=.true.)
      call check(error, ctx%verbosity == 3, "explicit verbosity applied")
      if (allocated(error)) return
      call check(error, ctx%debug, "explicit debug applied")
      call ctx%delete()
   end subroutine test_defaults_and_settings

   !> Several holders borrowing one context by pointer accumulate into the same
   !> timer tree - the core reason the context is shared by pointer
   subroutine test_shared_timer_tree(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(moist_context_type), pointer :: holder_a, holder_b

      call new_context(ctx)

      ! two independent holders borrow the very same context
      holder_a => ctx
      holder_b => ctx

      call holder_a%timer%start("A")
      call holder_a%timer%stop()
      call holder_b%timer%start("B")
      call holder_b%timer%stop()

      ! both holders' nodes live in the one shared tree, reachable from ctx
      call check(error, ctx%timer%num_nodes() == 2, "both holders share one tree")
      if (allocated(error)) return
      call check(error, ctx%timer%get("A") >= 0.0_wp, "node A is present")
      if (allocated(error)) return
      call check(error, ctx%timer%get("B") >= 0.0_wp, "node B is present")

      call ctx%delete()
   end subroutine test_shared_timer_tree

   !> writes(level) follows the verbosity>=level (or debug) print guard
   subroutine test_writes_guard(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, verbosity=2)
      call check(error, ctx%writes(1), "level 1 writes at verbosity 2")
      if (allocated(error)) return
      call check(error, ctx%writes(2), "level 2 writes at verbosity 2")
      if (allocated(error)) return
      call check(error,.not. ctx%writes(3), "level 3 is silent at verbosity 2")
      if (allocated(error)) return
      call ctx%delete()

      ! debug unlocks the diagnostic band even at low verbosity
      call new_context(ctx, verbosity=0, debug=.true.)
      call check(error, ctx%writes(3), "debug unlocks level 3")
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_writes_guard

   !> do_profile follows verbosity >= 3 by default and honours an explicit override.
   subroutine test_profile_flag(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      ! below the threshold: detailed profiling off
      call new_context(ctx, verbosity=3)
      call check(error,.not. ctx%do_profile, "verbosity 3 disables profiling")
      if (allocated(error)) return
      call ctx%delete()

      ! at/above the threshold: detailed profiling on
      call new_context(ctx, verbosity=4)
      call check(error, ctx%do_profile, "verbosity 4 enables profiling")
      if (allocated(error)) return
      call ctx%delete()

      ! explicit override wins over the verbosity default
      call new_context(ctx, verbosity=1, do_profile=.true.)
      call check(error, ctx%do_profile, "explicit do_profile overrides low verbosity")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, verbosity=5, do_profile=.false.)
      call check(error,.not. ctx%do_profile, "explicit do_profile overrides high verbosity")
      call ctx%delete()
   end subroutine test_profile_flag

   !> report_depth maps verbosity to the timing-tree print depth: 1->0, 2->1,
   !> 3->2, and >=4 -> unbounded
   subroutine test_report_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, verbosity=1)
      call check(error, ctx%report_depth() == 0, "verbosity 1 -> depth 0")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, verbosity=2)
      call check(error, ctx%report_depth() == 1, "verbosity 2 -> depth 1")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, verbosity=3)
      call check(error, ctx%report_depth() == 2, "verbosity 3 -> depth 2")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, verbosity=4)
      call check(error, ctx%report_depth() > 1000, "verbosity 4 -> unbounded depth")
      call ctx%delete()
   end subroutine test_report_depth

   !> A fresh context reports at least one thread (1 without OpenMP), and records
   !> a start timestamp
   subroutine test_threads_default(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx)
      call check(error, ctx%get_num_threads() >= 1, "default thread count is >= 1")
      if (allocated(error)) return
      call check(error, allocated(ctx%start_time), "start timestamp recorded")
      if (allocated(error)) return
      call check(error, len(ctx%start_time) == 19, "timestamp has YYYY-MM-DD HH:MM:SS length")
      call ctx%delete()
   end subroutine test_threads_default

   !> An explicit thread count is honoured (deterministic with or without OpenMP)
   subroutine test_threads_explicit(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=1)
      call check(error, ctx%get_num_threads() == 1, "explicit nthreads=1 applied")
      call ctx%delete()
   end subroutine test_threads_explicit

   !> set_num_threads retunes the recorded thread count after construction, and
   !> releasing the pin puts the OpenMP runtime back where it was
   !>
   !> The baseline is read through `get_num_threads` with no pin active, which is
   !> the live environment value, so the round trip is asserted exactly rather
   !> than as ">= 1" -- the latter is satisfied by a leaked pin of 1 and would
   !> not detect the pin failing to release
   subroutine test_set_num_threads(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      !> Live environment thread budget, before moist pins anything
      integer :: baseline

      call new_context(ctx)
      baseline = ctx%get_num_threads()
      ! an explicit positive count is honoured deterministically (OMP or not)
      call ctx%set_num_threads(1)
      call check(error, ctx%get_num_threads() == 1, "set_num_threads(1) records 1")
      if (allocated(error)) return
      call check(error, ctx%nthreads_pin == 1, "set_num_threads(1) pins 1")
      if (allocated(error)) return
      ! a non-positive request releases the pin -> back to the environment budget
      call ctx%set_num_threads(0)
      call check(error, ctx%nthreads_pin == 0, "set_num_threads(0) releases the pin")
      if (allocated(error)) return
      call check(error, ctx%get_num_threads() == baseline, &
                 "released pin restores the environment thread budget")
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_set_num_threads

   !> Deleting a context with an active pin restores the environment budget too:
   !> the pin lives in a global OpenMP control, so it must not outlive its owner
   subroutine test_delete_releases_pin(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx, probe
      !> Live environment thread budget, before moist pins anything
      integer :: baseline

      call new_context(probe)
      baseline = probe%get_num_threads()
      call probe%delete()

      call new_context(ctx, nthreads=1)
      call check(error, ctx%get_num_threads() == 1, "constructor pin applied")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(probe)
      call check(error, probe%get_num_threads() == baseline, &
                 "delete restores the environment thread budget")
      call probe%delete()
   end subroutine test_delete_releases_pin

   !> The context opens, writes to, and closes an owned log file
   subroutine test_owned_logfile(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_owned_logfile.tmp"
      integer :: iu, stat
      logical :: is_open
      character(64) :: firstline

      call new_context(ctx, logfile=path)
      call check(error, ctx%io_stat == 0, "log file opened cleanly")
      if (allocated(error)) return
      call check(error, ctx%owns_unit, "context owns the log unit")
      if (allocated(error)) return

      call ctx%message("hello from the context")
      call ctx%delete()

      ! The owned unit is closed after teardown. Ask by file rather than by unit
      ! number: the suite runs its cases concurrently, and a `newunit` in another
      ! case can be handed the number we just freed, which would make a
      ! unit-based inquire report our file as still open.
      inquire (file=path, opened=is_open)
      call check(error,.not. is_open, "owned unit closed on delete")
      if (allocated(error)) return
      call check(error,.not. ctx%owns_unit, "ownership dropped on delete")
      if (allocated(error)) return

      ! and the file holds what we wrote
      open (newunit=iu, file=path, status='old', action='read', iostat=stat)
      call check(error, stat == 0, "log file exists after delete")
      if (allocated(error)) then
         return
      end if
      read (iu, '(a)', iostat=stat) firstline
      call check(error, stat == 0, "log file is non-empty")
      if (.not. allocated(error)) &
         call check(error, trim(firstline) == "hello from the context", "message written to file")
      close (iu, status='delete')
   end subroutine test_owned_logfile

   !> print_settings renders without error and produces output
   subroutine test_print_settings_runs(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_print_settings.tmp"
      integer :: iu, stat
      character(256) :: line
      logical :: saw_verbosity

      call new_context(ctx, verbosity=2, logfile=path)
      call ctx%print_settings()
      call ctx%delete()

      open (newunit=iu, file=path, status='old', action='read', iostat=stat)
      call check(error, stat == 0, "settings file exists")
      if (allocated(error)) return

      saw_verbosity = .false.
      do
         read (iu, '(a)', iostat=stat) line
         if (stat /= 0) exit
         if (index(line, "Verbosity") > 0) saw_verbosity = .true.
      end do
      close (iu, status='delete')
      call check(error, saw_verbosity, "settings block mentions Verbosity")
   end subroutine test_print_settings_runs

   !> debug_message is silent unless debug is enabled
   subroutine test_debug_message_gated(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_debug_message.tmp"
      integer :: iu, stat
      character(64) :: line

      ! debug off: no debug file is opened at all (no orphan file), and the
      ! debug stream falls back to the main output unit.
      call new_context(ctx, debug=.false., debugfile=path)
      call check(error,.not. ctx%owns_debug_unit, "no debug file opened when debug off")
      if (.not. allocated(error)) then
         call check(error, ctx%debug_unit == ctx%unit, "debug stream falls back to main unit")
      end if
      call ctx%debug_message("should not appear")
      call ctx%delete()
      if (allocated(error)) return

      ! debug on: the message lands in the debug file
      call new_context(ctx, debug=.true., debugfile=path)
      call ctx%debug_message("diagnostic line")
      call ctx%delete()
      open (newunit=iu, file=path, status='old', action='read', iostat=stat)
      call check(error, stat == 0, "debug file exists (debug on)")
      if (.not. allocated(error)) then
         read (iu, '(a)', iostat=stat) line
         call check(error, stat == 0, "debug file non-empty when debug on")
         if (.not. allocated(error)) then
            call check(error, trim(line) == "diagnostic line", "debug message written")
         end if
      end if
      close (iu, status='delete')
   end subroutine test_debug_message_gated

   !> delete() is safe to call, including twice, and reports no nodes after
   subroutine test_delete_is_safe(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx)
      call ctx%timer%start("work")
      call ctx%timer%stop()
      call ctx%delete()
      call check(error, ctx%timer%num_nodes() == 0, "delete clears the timer")
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_delete_is_safe

end module test_utils_context
