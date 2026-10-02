!> Test suite for the shared run context in moist_context
module test_utils_context
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
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
                  new_unittest("defaults_and_settings", test_defaults_and_settings), &
                  new_unittest("shared_timer_tree", test_shared_timer_tree), &
                  new_unittest("writes_guard", test_writes_guard), &
                  new_unittest("profile_flag", test_profile_flag), &
                  new_unittest("report_depth", test_report_depth), &
                  new_unittest("threads_default", test_threads_default), &
                  new_unittest("threads_explicit", test_threads_explicit), &
                  new_unittest("set_num_threads", test_set_num_threads), &
                  new_unittest("delete_releases_pin", test_delete_releases_pin), &
                  new_unittest("owned_logfile", test_owned_logfile), &
                  new_unittest("print_settings_runs", test_print_settings_runs), &
                  new_unittest("debug_message_gated", test_debug_message_gated), &
                  new_unittest("delete_is_safe", test_delete_is_safe), &
                  new_unittest("output_contracts", test_output_contracts), &
                  new_unittest("failed_files", test_failed_files), &
                  new_unittest("thread_runtime", test_thread_runtime) &
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
      call check(error,.not. ctx%writes(4), "debug does not unlock profiling band")
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_writes_guard

   !> do_profile follows verbosity >= 4 by default and honours an explicit override
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

      call new_context(ctx, verbosity=0)
      call check(error, ctx%report_depth() == 0, "silent verbosity clamps depth to zero")
      if (allocated(error)) return
      call ctx%delete()

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
      ! unit-based inquire report our file as still open
      inquire (file=path, opened=is_open)
      call check(error,.not. is_open, "owned unit closed on delete")
      if (allocated(error)) return
      call check(error,.not. ctx%owns_unit, "ownership dropped on delete")
      if (allocated(error)) return

      ! and the file holds what we wrote
      open (newunit=iu, file=path, status="old", action="read", iostat=stat)
      call check(error, stat == 0, "log file exists after delete")
      if (allocated(error)) then
         return
      end if
      read (iu, "(a)", iostat=stat) firstline
      call check(error, stat == 0, "log file is non-empty")
      if (.not. allocated(error)) then
         call check(error, trim(firstline) == "hello from the context", "message written to file")
      end if
      close (iu, status="delete")
   end subroutine test_owned_logfile

   !> print_settings renders without error and produces output
   subroutine test_print_settings_runs(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_print_settings.tmp"
      integer :: iu, stat
      character(256) :: line
      logical :: saw_verbosity, saw_threads, saw_log, main_empty, saw_debug, saw_profile, saw_time
      integer :: override
      character(16) :: expected_threads

      call new_context(ctx, verbosity=2, logfile=path)
      open (newunit=override, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      write (expected_threads, "(i0)") ctx%get_num_threads()
      call ctx%print_settings(unit=override)
      flush (ctx%unit, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      rewind (override, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call ctx%delete()

      open (newunit=iu, file=path, status="old", action="read", iostat=stat)
      call check(error, stat == 0, "settings file exists")
      if (allocated(error)) return

      read (iu, "(a)", iostat=stat) line
      main_empty = stat < 0
      saw_verbosity = .false.
      saw_threads = .false.
      saw_log = .false.
      saw_debug = .false.
      saw_profile = .false.
      saw_time = .false.
      do
         read (override, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (index(line, "Verbosity") > 0) saw_verbosity = index(line, "2") > 0
         if (index(line, "OMP threads") > 0) then
            saw_threads = line(scan(trim(line), " ", back=.true.) + 1:len_trim(line)) == trim(expected_threads)
         end if
         if (index(line, "Log file") > 0) saw_log = index(line, path) > 0
         if (index(line, "Debug ") > 0 .and. index(line, "Debug file") == 0) then
            saw_debug = line(len_trim(line):len_trim(line)) == "F"
         end if
         if (index(line, "Detailed profiling") > 0) saw_profile = line(len_trim(line):len_trim(line)) == "F"
         if (index(line, "Start time") > 0) saw_time = index(line, ctx%start_time) > 0
      end do
      close (iu, status="delete", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      close (override, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call check(error, main_empty .and. saw_verbosity .and. saw_threads .and. saw_log .and. &
                 saw_debug .and. saw_profile .and. saw_time, &
                 "settings honors override and renders verbosity, threads and log path")
   end subroutine test_print_settings_runs

   !> debug_message is silent unless debug is enabled
   subroutine test_debug_message_gated(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_debug_message.tmp"
      integer :: iu, stat
      character(64) :: line

      ! debug off: no debug file is opened at all (no orphan file), and the
      ! debug stream falls back to the main output unit
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
      open (newunit=iu, file=path, status="old", action="read", iostat=stat)
      call check(error, stat == 0, "debug file exists (debug on)")
      if (.not. allocated(error)) then
         read (iu, "(a)", iostat=stat) line
         call check(error, stat == 0, "debug file non-empty when debug on")
         if (.not. allocated(error)) then
            call check(error, trim(line) == "diagnostic line", "debug message written")
         end if
      end if
      close (iu, status="delete")
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

   !> Messages honor levels and route debug output to the selected stream
   subroutine test_output_contracts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      integer :: iu, du, stat
      logical :: opened, ok
      character(256) :: line

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call new_context(ctx, verbosity=1, unit=iu, debug=.true.)
      ok = ctx%unit == iu .and. ctx%debug_unit == iu
      call ctx%debug_message("fallback debug")
      call ctx%print_settings()
      call ctx%delete()
      inquire (unit=iu, opened=opened)
      rewind (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      read (iu, "(a)", iostat=stat) line
      ok = ok .and. opened .and. stat == 0 .and. trim(line) == "fallback debug"
      opened = .false.
      do
         read (iu, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (index(line, "Log file") > 0) opened = index(line, "(stdout)") > 0
      end do
      ok = ok .and. opened
      close (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call check(error, ok, "borrowed unit stays open and receives fallback debug")
      if (allocated(error)) return

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call new_context(ctx, verbosity=0, unit=iu)
      call ctx%message("silent default")
      call ctx%debug_message("silent debug")
      rewind (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      read (iu, "(a)", iostat=stat) line
      ok = stat < 0
      call ctx%delete()
      close (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call check(error, ok, "verbosity zero and disabled debug produce no output")
      if (allocated(error)) return

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call new_context(ctx, verbosity=1, unit=iu)
      call ctx%message("default visible")
      call ctx%message("level suppressed", level=2)
      rewind (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      read (iu, "(a)", iostat=stat) line
      ok = stat == 0 .and. trim(line) == "default visible"
      read (iu, "(a)", iostat=stat) line
      ok = ok .and. stat < 0
      call ctx%delete()
      close (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call check(error, ok, "message default and explicit levels are honored")
      if (allocated(error)) return

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call new_context(ctx, unit=iu, debug=.true., debugfile="context_split_debug.tmp")
      ok = ctx%owns_debug_unit .and. ctx%debugfile == "context_split_debug.tmp"
      call ctx%message("main only")
      call ctx%debug_message("debug only")
      call ctx%delete()
      inquire (file="context_split_debug.tmp", opened=opened)
      ok = ok .and. .not. opened .and. .not. ctx%owns_debug_unit
      rewind (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      read (iu, "(a)", iostat=stat) line
      ok = ok .and. stat == 0 .and. trim(line) == "main only"
      read (iu, "(a)", iostat=stat) line
      ok = ok .and. stat < 0
      close (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      open (newunit=du, file="context_split_debug.tmp", status="old", action="read", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      read (du, "(a)", iostat=stat) line
      ok = ok .and. stat == 0 .and. trim(line) == "debug only"
      close (du, status="delete", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call check(error, ok, "debug stream is separate and owned file closes")
   end subroutine test_output_contracts

   !> Failed owned-file opens retain a usable borrowed stream and expose status
   subroutine test_failed_files(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      integer :: iu, stat
      logical :: ok

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call new_context(ctx, unit=iu, logfile="missing_context_directory/log")
      ok = ctx%io_stat /= 0 .and. .not. ctx%owns_unit .and. ctx%unit == iu
      call ctx%delete()
      call new_context(ctx, unit=iu, debug=.true., debugfile="missing_context_directory/debug")
      ok = ok .and. ctx%io_stat /= 0 .and. .not. ctx%owns_debug_unit .and. ctx%debug_unit == iu
      call ctx%delete()
      close (iu, iostat=stat)
      call check(error, stat == 0, "test stream I/O succeeds")
      if (allocated(error)) return
      call check(error, ok, "failed file opens report error and preserve fallback")
   end subroutine test_failed_files

   !> Retuning restores the first host budget and follows later host changes
   subroutine test_thread_runtime(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      integer :: baseline, runtime
      logical :: ok

      call new_context(ctx)
      baseline = ctx%get_num_threads()
      call ctx%set_num_threads(3)
      runtime = 3
!$    runtime = omp_get_max_threads()
      ok = ctx%get_num_threads() == 3 .and. runtime == 3
      call ctx%set_num_threads(4)
      call ctx%set_num_threads(-1)
      runtime = baseline
!$    runtime = omp_get_max_threads()
      ok = ok .and. runtime == baseline .and. ctx%get_num_threads() == baseline
      ok = ok .and. ctx%nthreads_env == 0 .and. ctx%nthreads_pin == 0
!$    call omp_set_num_threads(5)
!$    ok = ok .and. ctx%get_num_threads() == 5
      call ctx%set_num_threads(0)
!$    ok = ok .and. omp_get_max_threads() == 5
!$    call omp_set_num_threads(baseline)
      call ctx%delete()
      call new_context(ctx, nthreads=3)
      ok = ok .and. ctx%get_num_threads() == 3
      call ctx%delete()
      ok = ok .and. ctx%nthreads_pin == 0
!$    ok = ok .and. omp_get_max_threads() == baseline
!$    call omp_set_num_threads(baseline)
      call check(error, ok, "pin applies, retune preserves baseline, release follows host")
   end subroutine test_thread_runtime

end module test_utils_context
