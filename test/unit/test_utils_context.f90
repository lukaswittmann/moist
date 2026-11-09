!> Test suite for the shared run context in moist_context
module test_utils_context
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use mctc_env, only: wp, host_error_type => error_type
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
                  new_unittest("verbosity_guard", test_verbosity_guard), &
                  new_unittest("profile_flag", test_profile_flag), &
                  new_unittest("report_depth", test_report_depth), &
                  new_unittest("host_timer_group", test_host_timer_group), &
                  new_unittest("host_timer_balance", test_host_timer_balance), &
                  new_unittest("host_timer_print_open", test_host_timer_print_open), &
                  new_unittest("host_timer_names", test_host_timer_names), &
                  new_unittest("threads_default", test_threads_default), &
                  new_unittest("threads_explicit", test_threads_explicit), &
                  new_unittest("threads_fixed", test_threads_fixed), &
                  new_unittest("threads_leave_runtime", test_threads_leave_runtime), &
                  new_unittest("owned_logfile", test_owned_logfile), &
                  new_unittest("print_settings_runs", test_print_settings_runs), &
                  new_unittest("debug_message_gated", test_debug_message_gated), &
                  new_unittest("message_flush_gated", test_message_flush_gated), &
                  new_unittest("warning_in_output", test_warning_in_output), &
                  new_unittest("delete_is_safe", test_delete_is_safe), &
                  new_unittest("output_contracts", test_output_contracts), &
                  new_unittest("failed_files", test_failed_files) &
                  ]

   end subroutine collect_utils_context

   !> new_context applies defaults and honours explicit settings
   subroutine test_defaults_and_settings(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=0)
      call check(error, ctx%verbosity() == 2, "default verbosity is 2")
      if (allocated(error)) return
      call check(error,.not. ctx%debug, "default debug is false")
      if (allocated(error)) return
      call check(error,.not. ctx%warnings_to_stderr, "default keeps warnings off stderr")
      if (allocated(error)) return
      ! a fresh context owns an initialized, empty timer tree
      call check(error, ctx%timer%num_nodes() == 0, "fresh timer has no nodes")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=3)
      call check(error, ctx%verbosity() == 3, "explicit verbosity applied")
      if (allocated(error)) return
      call ctx%delete()

      ! debug raises the verbosity to the debug level
      call new_context(ctx, nthreads=0, verbosity=2, debug=.true.)
      call check(error, ctx%debug, "explicit debug applied")
      if (allocated(error)) return
      call check(error, ctx%verbosity() == 6, "debug raises verbosity to 6")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, warnings_to_stderr=.true.)
      call check(error, ctx%warnings_to_stderr, "explicit warnings_to_stderr applied")
      if (allocated(error)) return
      call ctx%delete()

      ! the debug verbosity alone does not enable debug
      call new_context(ctx, nthreads=0, verbosity=6)
      call check(error,.not. ctx%debug, "verbosity 6 leaves debug off")
      call ctx%delete()
   end subroutine test_defaults_and_settings

   !> Several holders borrowing one context by pointer accumulate into the same
   !> timer tree - the core reason the context is shared by pointer
   subroutine test_shared_timer_tree(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type), target :: ctx
      type(moist_context_type), pointer :: holder_a, holder_b

      call new_context(ctx, nthreads=0)

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

   !> verbosity(level) follows the verbosity>=level print guard
   subroutine test_verbosity_guard(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=0, verbosity=2)
      call check(error, ctx%verbosity(1), "level 1 writes at verbosity 2")
      if (allocated(error)) return
      call check(error, ctx%verbosity(2), "level 2 writes at verbosity 2")
      if (allocated(error)) return
      call check(error,.not. ctx%verbosity(3), "level 3 is silent at verbosity 2")
      if (allocated(error)) return
      call ctx%delete()

      ! debug unlocks every level even at zero verbosity
      call new_context(ctx, nthreads=0, verbosity=0, debug=.true.)
      call check(error, ctx%verbosity(6), "debug unlocks level 6")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=5)
      call check(error,.not. ctx%verbosity(6), "level 6 is silent at verbosity 5")
      if (allocated(error)) return
      call ctx%delete()
   end subroutine test_verbosity_guard

   !> do_profile is off by default at any verbosity and set only explicitly
   subroutine test_profile_flag(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=0, verbosity=6)
      call check(error,.not. ctx%do_profile, "verbosity 6 leaves profiling off")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, debug=.true.)
      call check(error,.not. ctx%do_profile, "debug leaves profiling off")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=0, do_profile=.true.)
      call check(error, ctx%do_profile, "explicit do_profile enables profiling")
      call ctx%delete()
   end subroutine test_profile_flag

   !> report_depth maps verbosity to the timing-tree print depth: 1->0, 2->1,
   !> 3->2, and >=4 -> unbounded
   subroutine test_report_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)
      call check(error, ctx%report_depth() == 0, "silent verbosity clamps depth to zero")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=1)
      call check(error, ctx%report_depth() == 0, "verbosity 1 -> depth 0")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=2)
      call check(error, ctx%report_depth() == 1, "verbosity 2 -> depth 1")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=3)
      call check(error, ctx%report_depth() == 2, "verbosity 3 -> depth 2")
      if (allocated(error)) return
      call ctx%delete()

      call new_context(ctx, nthreads=0, verbosity=4)
      call check(error, ctx%report_depth() > 1000, "verbosity 4 -> unbounded depth")
      call ctx%delete()
   end subroutine test_report_depth

   !> Host timers land under `Host` and nest only among themselves: a moist timer
   !> started inside a host frame stays top level, a host timer started inside a
   !> moist frame still lands under `Host`
   subroutine test_host_timer_group(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      type(host_error_type), allocatable :: err

      call new_context(ctx, nthreads=0)
      call ctx%start_timer("Coupling", err)
      if (.not. allocated(err)) call ctx%timer%start("Setup")
      if (.not. allocated(err)) call ctx%start_timer("potential", err)
      if (.not. allocated(err)) call ctx%stop_timer("potential", err)
      if (.not. allocated(err)) call ctx%timer%stop("Setup")
      if (.not. allocated(err)) call ctx%stop_timer("Coupling", err)
      call check(error,.not. allocated(err), "balanced host timers are accepted")
      if (allocated(error)) return

      call check(error, node_depth_of(ctx, "Host") == 0, "Host is a top-level node")
      if (allocated(error)) return
      call check(error, node_depth_of(ctx, "Setup") == 0, "moist timer stays top level in a host frame")
      if (allocated(error)) return
      call check(error, node_depth_of(ctx, "Coupling") == 1, "host timer lands under Host")
      if (allocated(error)) return
      call check(error, node_depth_of(ctx, "potential") == 2, "host timer in a moist frame stays under Host")
      if (allocated(error)) return
      call check(error, ctx%timer%get("Host/Coupling/potential") >= 0.0_wp, "nested host path resolves")
      if (allocated(error)) return
      call check(error,.not. ieee_is_nan(ctx%timer%get("Host")), "closing the last host timer stops Host")
      call ctx%delete()
   end subroutine test_host_timer_group

   !> A stop must name the innermost open host timer; a rejected call changes nothing
   subroutine test_host_timer_balance(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      type(host_error_type), allocatable :: err

      call new_context(ctx, nthreads=0)
      call ctx%stop_timer("A", err)
      call check(error, allocated(err), "stop without an open host timer is refused")
      if (allocated(error)) return

      call ctx%start_timer("", err)
      call check(error, allocated(err), "empty name is refused")
      if (allocated(error)) return
      call ctx%start_timer(repeat("x", 33), err)
      call check(error, allocated(err), "name past 32 characters is refused")
      if (allocated(error)) return
      call check(error, ctx%timer%num_nodes() == 0, "refused starts create no node")
      if (allocated(error)) return

      call ctx%start_timer(repeat("x", 32), err)
      call check(error,.not. allocated(err), "32-character name is accepted")
      if (allocated(error)) return
      call ctx%stop_timer(repeat("x", 32), err)
      call check(error,.not. allocated(err), "32-character name stops")
      if (allocated(error)) return

      call ctx%start_timer("A", err)
      if (.not. allocated(err)) call ctx%start_timer("B", err)
      call check(error,.not. allocated(err), "nested host timers start")
      if (allocated(error)) return
      call ctx%stop_timer("A", err)
      call check(error, allocated(err), "stop of an outer host timer is refused")
      if (allocated(error)) return
      call ctx%stop_timer("B", err)
      if (.not. allocated(err)) call ctx%stop_timer("A", err)
      call check(error,.not. allocated(err), "the refused stop left the stack intact")
      if (allocated(error)) return
      call ctx%stop_timer("A", err)
      call check(error, allocated(err), "a closed host timer cannot stop again")
      if (allocated(error)) return

      call check(error,.not. ieee_is_nan(ctx%timer%get("Host/A/B")), "balanced timers are not poisoned")
      call ctx%delete()
   end subroutine test_host_timer_balance

   !> A report printed while a host timer is open leaves it and `Host` measuring
   subroutine test_host_timer_print_open(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      type(host_error_type), allocatable :: err
      integer :: iu, stat

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "scratch unit for the timer report opens")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0)
      call ctx%start_timer("A", err)
      if (.not. allocated(err)) call ctx%timer%write(iu)
      if (.not. allocated(err)) call ctx%stop_timer("A", err)
      close (iu)
      call check(error,.not. allocated(err), "a host timer stops after a report")
      if (allocated(error)) return
      call check(error,.not. ieee_is_nan(ctx%timer%get("Host/A")), "the report left the host timer measuring")
      if (allocated(error)) return
      call check(error,.not. ieee_is_nan(ctx%timer%get("Host")), "the report left Host measuring")
      call ctx%delete()
   end subroutine test_host_timer_print_open

   !> Every byte of a host timer name counts: 1 to 32 bytes of UTF-8, no control
   !> character, no leading or trailing space, for a start and a stop alike
   subroutine test_host_timer_names(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      type(host_error_type), allocatable :: err
      !> UTF-8 encodings of U+00E9, U+20AC and U+1F600
      character(len=*), parameter :: e_acute = achar(195)//achar(169), &
         & euro = achar(226)//achar(130)//achar(172), &
         & smiley = achar(240)//achar(159)//achar(152)//achar(128)

      call new_context(ctx, nthreads=0)
      call refuse("", "empty")
      call refuse(repeat("x", 32)//"  ", "exceeds 32 bytes")
      call refuse(repeat(e_acute, 17), "exceeds 32 bytes")
      call refuse(" A", "space")
      call refuse("A ", "space")
      call refuse("A"//achar(9)//"B", "control")
      call refuse("A"//achar(0), "control")
      call refuse("A"//achar(194)//achar(133), "control")
      call refuse(achar(233), "UTF-8")
      call refuse(achar(195), "UTF-8")
      call refuse(achar(192)//achar(129), "UTF-8")
      call refuse(achar(224)//achar(128)//achar(128), "UTF-8")
      call refuse(achar(237)//achar(160)//achar(128), "UTF-8")
      call refuse(achar(244)//achar(144)//achar(128)//achar(128), "UTF-8")
      if (allocated(error)) return
      call check(error, ctx%timer%num_nodes() == 0, "refused names create no node")
      if (allocated(error)) return

      call ctx%start_timer(repeat(e_acute, 16), err)
      call check(error,.not. allocated(err), "32 bytes of UTF-8 are accepted")
      if (allocated(error)) return
      call ctx%start_timer("A"//euro//smiley, err)
      call check(error,.not. allocated(err), "three- and four-byte characters are accepted")
      if (allocated(error)) return
      call ctx%stop_timer("A"//euro//smiley//" ", err)
      call check(error, allocated(err), "a stop with a trailing space is refused")
      if (allocated(error)) return
      call ctx%stop_timer("A"//euro//smiley, err)
      if (.not. allocated(err)) call ctx%stop_timer(repeat(e_acute, 16), err)
      call check(error,.not. allocated(err), "the refused stop left the timers open")
      call ctx%delete()

   contains

      !> Check that a start with `name` is refused by the rule named in `rule`
      subroutine refuse(name, rule)
         character(len=*), intent(in) :: name, rule

         if (allocated(error)) return
         call ctx%start_timer(name, err)
         call check(error, allocated(err), "invalid host timer name is refused: "//rule)
         if (allocated(error)) return
         call check(error, index(err%message, rule) > 0, "refusal names the rule '"//rule//"', got: "//err%message)
      end subroutine refuse

   end subroutine test_host_timer_names

   !> Depth of the first timer node with a name, -1 when absent
   !>
   !> @param[in] ctx   context holding the timer
   !> @param[in] name  node name
   integer function node_depth_of(ctx, name) result(depth)
      type(moist_context_type), intent(in) :: ctx
      character(*), intent(in) :: name
      integer :: id

      depth = -1
      do id = 1, ctx%timer%num_nodes()
         if (ctx%timer%node_name(id) == name) then
            depth = ctx%timer%node_depth(id)
            return
         end if
      end do
   end function node_depth_of

   !> A fresh context reports at least one thread (1 without OpenMP), and records
   !> a start timestamp
   subroutine test_threads_default(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=0)
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

   !> The thread count is fixed at construction: an unset count takes the
   !> OpenMP environment of that moment and ignores later host changes
   subroutine test_threads_fixed(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx, later
      !> Environment thread budget before the host changes it
      integer :: baseline
      logical :: ok

      call new_context(ctx, nthreads=0)
      baseline = ctx%get_num_threads()
      ok = .true.
!$    ok = baseline == omp_get_max_threads()
!$    call omp_set_num_threads(baseline + 3)
      call new_context(later, nthreads=0)
      ok = ok .and. ctx%get_num_threads() == baseline
!$    ok = ok .and. later%get_num_threads() == baseline + 3
!$    call omp_set_num_threads(baseline)
      call later%delete()
      call ctx%delete()
      call check(error, ok, "unset count resolved once at construction")
   end subroutine test_threads_fixed

   !> Explicit counts leave the host's OpenMP runtime untouched: each context
   !> reports its own count and the runtime keeps its setting throughout
   subroutine test_threads_leave_runtime(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: first, second
      !> Environment thread budget and the runtime seen while contexts live
      integer :: baseline, runtime
      logical :: ok

      baseline = 1
!$    baseline = omp_get_max_threads()
      call new_context(first, nthreads=baseline + 1)
      call new_context(second, nthreads=baseline + 2)
      runtime = baseline
!$    runtime = omp_get_max_threads()
      ok = runtime == baseline .and. first%get_num_threads() == baseline + 1 .and. &
         & second%get_num_threads() == baseline + 2
      call first%delete()
      call second%delete()
      call check(error, ok, "explicit counts leave the runtime alone")
   end subroutine test_threads_leave_runtime

   !> The context opens, writes to, and closes an owned log file
   subroutine test_owned_logfile(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_owned_logfile.tmp"
      integer :: iu, stat
      logical :: is_open
      character(64) :: firstline

      call new_context(ctx, nthreads=0, logfile=path)
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

   !> print_settings honors a unit override
   subroutine test_print_settings_runs(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_print_settings.tmp"
      integer :: iu, stat, log_stat
      character(256) :: line
      logical :: main_empty
      integer :: override

      call new_context(ctx, nthreads=0, verbosity=2, logfile=path)
      open (newunit=override, status="scratch", action="readwrite", iostat=stat)
      if (stat == 0) call ctx%print_settings(unit=override)
      call ctx%delete()
      call check(error, stat == 0, "scratch override unit opens")
      if (allocated(error)) then
         call remove_file(path)
         return
      end if
      close (override)

      ! All output went to the override, so the owned log file stays empty
      main_empty = .false.
      open (newunit=iu, file=path, status="old", action="read", iostat=log_stat)
      if (log_stat == 0) then
         read (iu, "(a)", iostat=stat) line
         main_empty = stat < 0
         close (iu, status="delete")
      end if

      call check(error, log_stat == 0, "settings file exists")
      if (allocated(error)) return
      call check(error, main_empty, "unit override leaves the owned log file empty")
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
      call new_context(ctx, nthreads=0, debug=.false., debugfile=path)
      call check(error,.not. ctx%owns_debug_unit, "no debug file opened when debug off")
      if (.not. allocated(error)) then
         call check(error, ctx%debug_unit == ctx%unit, "debug stream falls back to main unit")
      end if
      call ctx%debug_message("should not appear")
      call ctx%delete()
      if (allocated(error)) return

      ! debug on: the message lands in the debug file
      call new_context(ctx, nthreads=0, debug=.true., debugfile=path)
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

   !> message_flush writes like message: only at an enabled level
   subroutine test_message_flush_gated(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_message_flush.tmp"
      integer :: iu, stat
      character(64) :: line

      call new_context(ctx, nthreads=0, verbosity=1, logfile=path)
      call ctx%message_flush("level two is silent", level=2)
      call ctx%message_flush("level one is written")
      call ctx%delete()
      open (newunit=iu, file=path, status="old", action="read", iostat=stat)
      call check(error, stat == 0, "log file exists")
      if (allocated(error)) return
      read (iu, "(a)", iostat=stat) line
      call check(error, stat == 0, "log file is non-empty")
      if (.not. allocated(error)) then
         call check(error, trim(line) == "level one is written", "only the enabled message written")
      end if
      if (.not. allocated(error)) then
         read (iu, "(a)", iostat=stat) line
         call check(error, is_iostat_end(stat), "silenced message not written")
      end if
      close (iu, status="delete")
   end subroutine test_message_flush_gated

   !> A warning lands in the output unit even at verbosity 0, in order with the messages
   subroutine test_warning_in_output(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      character(*), parameter :: path = "test_context_warning.tmp"
      integer :: iu, stat
      character(64) :: line

      call new_context(ctx, nthreads=0, verbosity=0, logfile=path)
      call ctx%message("silent at verbosity 0")
      call ctx%warning("always written")
      call ctx%delete()
      open (newunit=iu, file=path, status="old", action="read", iostat=stat)
      call check(error, stat == 0, "log file exists")
      if (allocated(error)) return
      read (iu, "(a)", iostat=stat) line
      call check(error, stat == 0, "log file is non-empty")
      if (.not. allocated(error)) then
         call check(error, trim(line) == "[Warning] always written", "warning written to the log")
      end if
      if (.not. allocated(error)) then
         read (iu, "(a)", iostat=stat) line
         call check(error, is_iostat_end(stat), "only the warning written")
      end if
      close (iu, status="delete")
   end subroutine test_warning_in_output

   !> delete() is safe to call, including twice, and reports no nodes after
   subroutine test_delete_is_safe(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx

      call new_context(ctx, nthreads=0)
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
      logical :: main_routed, debug_routed, fallback_read, unnamed
      logical :: silent, levels_honored, owned, named, closed, released, main_only, debug_only
      character(256) :: line
      character(*), parameter :: split_path = "test_context_split_debug.tmp"

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "scratch unit for the fallback stream opens")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=1, unit=iu, debug=.true.)
      main_routed = ctx%unit == iu
      debug_routed = ctx%debug_unit == iu
      ! A borrowed unit has no owned file name
      unnamed = .not. allocated(ctx%logfile)
      call ctx%debug_message("fallback debug")
      call ctx%delete()
      ! Reading the borrowed unit back after delete proves it stayed open
      fallback_read = .false.
      rewind (iu, iostat=stat)
      if (stat == 0) then
         read (iu, "(a)", iostat=stat) line
         fallback_read = stat == 0 .and. trim(line) == "fallback debug"
      end if
      close (iu)
      call check(error, main_routed, "borrowed unit is the main stream")
      if (allocated(error)) return
      call check(error, debug_routed, "debug stream falls back to the borrowed unit")
      if (allocated(error)) return
      call check(error, fallback_read, "borrowed unit stays open and receives fallback debug")
      if (allocated(error)) return
      call check(error, unnamed, "a borrowed unit names no owned log file")
      if (allocated(error)) return

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "scratch unit for the silent stream opens")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=0, unit=iu)
      call ctx%message("silent default")
      call ctx%debug_message("silent debug")
      call ctx%delete()
      silent = .false.
      rewind (iu, iostat=stat)
      if (stat == 0) then
         read (iu, "(a)", iostat=stat) line
         silent = stat < 0
      end if
      close (iu)
      call check(error, silent, "verbosity zero and disabled debug produce no output")
      if (allocated(error)) return

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "scratch unit for the level stream opens")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, verbosity=1, unit=iu)
      call ctx%message("default visible")
      call ctx%message("level suppressed", level=2)
      call ctx%delete()
      levels_honored = .false.
      rewind (iu, iostat=stat)
      if (stat == 0) then
         read (iu, "(a)", iostat=stat) line
         levels_honored = stat == 0 .and. trim(line) == "default visible"
         read (iu, "(a)", iostat=stat) line
         levels_honored = levels_honored .and. stat < 0
      end if
      close (iu)
      call check(error, levels_honored, "message default and explicit levels are honored")
      if (allocated(error)) return

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "scratch unit for the split stream opens")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, unit=iu, debug=.true., debugfile=split_path)
      owned = ctx%owns_debug_unit
      ! The name is only recorded once the file opened; never read it unallocated
      named = .false.
      if (allocated(ctx%debugfile)) named = ctx%debugfile == split_path
      call ctx%message("main only")
      call ctx%debug_message("debug only")
      call ctx%delete()
      inquire (file=split_path, opened=closed)
      closed = .not. closed
      released = .not. ctx%owns_debug_unit
      main_only = .false.
      rewind (iu, iostat=stat)
      if (stat == 0) then
         read (iu, "(a)", iostat=stat) line
         main_only = stat == 0 .and. trim(line) == "main only"
         read (iu, "(a)", iostat=stat) line
         main_only = main_only .and. stat < 0
      end if
      close (iu)
      debug_only = .false.
      open (newunit=du, file=split_path, status="old", action="read", iostat=stat)
      if (stat == 0) then
         read (du, "(a)", iostat=stat) line
         debug_only = stat == 0 .and. trim(line) == "debug only"
         close (du, status="delete")
      end if
      call check(error, owned, "context owns the split debug file")
      if (allocated(error)) return
      call check(error, named, "split debug file is recorded by name")
      if (allocated(error)) return
      call check(error, closed, "owned debug file closes on delete")
      if (allocated(error)) return
      call check(error, released, "debug ownership dropped on delete")
      if (allocated(error)) return
      call check(error, main_only, "main stream receives only the message")
      if (allocated(error)) return
      call check(error, debug_only, "debug file receives only the debug message")
   end subroutine test_output_contracts

   !> Failed owned-file opens retain a usable borrowed stream and expose status
   subroutine test_failed_files(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_context_type) :: ctx
      integer :: iu, stat
      logical :: log_reported, log_unowned, log_fallback, debug_reported, debug_unowned, debug_fallback

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      call check(error, stat == 0, "scratch unit for the fallback stream opens")
      if (allocated(error)) return
      call new_context(ctx, nthreads=0, unit=iu, logfile="missing_context_directory/log")
      log_reported = ctx%io_stat /= 0
      log_unowned = .not. ctx%owns_unit
      log_fallback = ctx%unit == iu
      call ctx%delete()
      call new_context(ctx, nthreads=0, unit=iu, debug=.true., debugfile="missing_context_directory/debug")
      debug_reported = ctx%io_stat /= 0
      debug_unowned = .not. ctx%owns_debug_unit
      debug_fallback = ctx%debug_unit == iu
      call ctx%delete()
      close (iu)
      call check(error, log_reported, "failed log file open is reported")
      if (allocated(error)) return
      call check(error, log_unowned, "failed log file is not owned")
      if (allocated(error)) return
      call check(error, log_fallback, "failed log file keeps the borrowed unit")
      if (allocated(error)) return
      call check(error, debug_reported, "failed debug file open is reported")
      if (allocated(error)) return
      call check(error, debug_unowned, "failed debug file is not owned")
      if (allocated(error)) return
      call check(error, debug_fallback, "failed debug file keeps the borrowed unit")
   end subroutine test_failed_files

   !> Delete a file left behind by a test that failed before its own cleanup
   !>
   !> @param[in] path File to delete, ignored when absent
   subroutine remove_file(path)
      !> File to delete
      character(*), intent(in) :: path
      integer :: iu, stat

      open (newunit=iu, file=path, status="old", iostat=stat)
      if (stat == 0) close (iu, status="delete")
   end subroutine remove_file

end module test_utils_context
