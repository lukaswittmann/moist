!> Shared run context for a moist calculation
!>
!> A single `moist_context_type` instance is constructed once at the top of a
!> run (the C-API handle wrapper or the CLI driver) and then *borrowed* by
!> every model and cavity through a non-owning pointer set at construction
!>
!> Bundles the run-wide settings that were previously duplicated on each
!> type:
!>
!>   * `verbosity` - the project-wide output level (0=silent, 1=summary,
!>     2=user events, 3=diagnostics),
!>   * `debug`     - a diagnostic flag OR-ed with the top verbosity level,
!>   * `timer`     - the profiling timer; because the context is shared by
!>     pointer, every borrower appends into this one timer tree
!>
module moist_context
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use moist_utils_timer, only: timer_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
!$ use omp_lib, only: omp_get_max_threads
   implicit none(type, external)
   private

   public :: moist_context_type, new_context, resolve_num_threads

   !> Shared run context: owned at the top level, borrowed by every model/cavity
   type :: moist_context_type
      !> Output unit for informational messages (errors still go to error_unit)
      integer :: unit = output_unit
      !> Output unit for debug dumps; defaults to `unit` (a separate `debugfile`
      !> splits the two streams)
      integer :: debug_unit = output_unit
      !> Verbosity level (0=silent, 1=summary, 2=user events, 3=diagnostics)
      integer :: verbosity = 1
      !> Debug flag; OR-ed with the top verbosity level at diagnostic sites
      logical :: debug = .false.
      !> Detailed-profiling flag, default verbosity >= 4
      !>
      !> - when set, fine-grained (per-part) timers are recorded on top of the
      !>   coarse module timings, e.g. the individual gradient sub-steps in
      !>   the DROP hot loop
      logical :: do_profile = .false.
      !> OpenMP thread count, fixed by `new_context`; 0 only before construction
      integer, private :: nthreads = 0
      !> Run start timestamp, formatted `YYYY-MM-DD HH:MM:SS`
      character(:), allocatable :: start_time
      !> Name of the owned main output file, if any (for display)
      character(:), allocatable :: logfile
      !> Name of the owned debug output file, if any (for display)
      character(:), allocatable :: debugfile
      !> Whether the context opened (and must close) `unit`
      logical :: owns_unit = .false.
      !> Whether the context opened (and must close) `debug_unit`
      logical :: owns_debug_unit = .false.
      !> iostat from opening the owned file(s); non-zero signals a failed open
      integer :: io_stat = 0
      !> Central profiling timer; every borrower appends into this one tree
      type(timer_type) :: timer
   contains
      !> Tear down the owned timer and close any owned files
      procedure :: delete => delete_context
      !> Predicate for the "print at level N" guard idiom
      procedure :: writes => context_writes
      !> Write a message to the context output unit, gated by verbosity
      procedure :: message => context_message
      !> Write a debug message to the debug unit when debug is enabled
      procedure :: debug_message => context_debug_message
      !> OpenMP thread count fixed at construction
      procedure :: get_num_threads => context_get_num_threads
      !> Render the run settings via the pretty printer
      procedure :: print_settings => context_print_settings
      !> Maximum timing-tree depth to print for the current verbosity
      procedure :: report_depth => context_report_depth
   end type moist_context_type

contains

   !> Initialize a run context and its timer
   !>
   !> The thread count is fixed here and never changes afterwards: a positive
   !> `nthreads` as given, otherwise the calling thread's `omp_get_max_threads`
   !> at construction, otherwise 1
   !>
   !> The count sizes MOIST's own OpenMP teams through `num_threads` clauses;
   !> the host's OpenMP runtime and the BLAS/LAPACK threading are left alone
   !>
   !> Passing `logfile` makes the context open and own that file as its main
   !> output unit, `debugfile` likewise for the debug stream; ownership follows
   !> the single-owner contract, `delete` closes whatever `new_context` opened
   !>
   !> @param[out] self       context to initialize
   !> @param[in]  nthreads   thread count; absent or <= 0 takes the current OpenMP setting
   !> @param[in]  verbosity  output level (default 1)
   !> @param[in]  debug      diagnostic flag (default .false.)
   !> @param[in]  unit       already-open output unit (default output_unit)
   !> @param[in]  do_profile detailed-profiling flag; defaults to verbosity >= 4
   !> @param[in]  logfile    path of a main output file for the context to own
   !> @param[in]  debugfile  path of a separate debug output file to own
   subroutine new_context(self, nthreads, verbosity, debug, unit, do_profile, logfile, debugfile)
      !> Context to initialize
      type(moist_context_type), intent(out) :: self
      !> OpenMP thread count; absent or <= 0 takes the current OpenMP setting
      integer, intent(in), optional :: nthreads
      !> Output verbosity level
      integer, intent(in), optional :: verbosity
      !> Diagnostic flag
      logical, intent(in), optional :: debug
      !> Output unit for informational messages
      integer, intent(in), optional :: unit
      !> Detailed-profiling flag; defaults to verbosity >= 4
      logical, intent(in), optional :: do_profile
      !> Path of a main output file to open and own
      character(*), intent(in), optional :: logfile
      !> Path of a separate debug output file to open and own
      character(*), intent(in), optional :: debugfile

      !> Values used to format the start timestamp
      character(8) :: date
      character(10) :: time
      !> Scratch unit / iostat for opening owned files
      integer :: iu, stat

      if (present(verbosity)) self%verbosity = verbosity
      if (present(debug)) self%debug = debug
      if (present(unit)) self%unit = unit

      !> Detailed profiling follows the verbosity by default (level >= 4, i.e.
      !> the "show everything" band), but a caller may force it on/off
      self%do_profile = self%verbosity >= 4
      if (present(do_profile)) self%do_profile = do_profile

      !> Thread count, resolved once; the OpenMP runtime itself is never changed
      self%nthreads = resolve_num_threads(nthreads)

      !> Run start timestamp
      call date_and_time(date=date, time=time)
      self%start_time = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "// &
         & time(1:2)//":"//time(3:4)//":"//time(5:6)

      !> Own a main output file when a path is given (takes precedence over unit)
      if (present(logfile)) then
         open (newunit=iu, file=logfile, status="replace", action="write", iostat=stat)
         self%io_stat = stat
         if (stat == 0) then
            self%unit = iu
            self%owns_unit = .true.
            self%logfile = logfile
         end if
      end if

      ! Debug output defaults to the main unit; a separate file splits the stream,
      ! but only when debug is actually enabled -- otherwise no orphan file is
      ! created (nothing would ever be written to it)
      self%debug_unit = self%unit
      if (present(debugfile) .and. self%debug) then
         open (newunit=iu, file=debugfile, status="replace", action="write", iostat=stat)
         if (stat /= 0) self%io_stat = stat
         if (stat == 0) then
            self%debug_unit = iu
            self%owns_debug_unit = .true.
            self%debugfile = debugfile
         end if
      end if

      call self%timer%new(verbose=self%verbosity > 1)

   end subroutine new_context

   !> Release the resources owned by the context (its timer and any owned files)
   subroutine delete_context(self)
      !> Context to tear down
      class(moist_context_type), intent(inout) :: self
      !> Whether the unit is still open at teardown
      logical :: is_open

      if (self%owns_unit) then
         inquire (unit=self%unit, opened=is_open)
         if (is_open) close (self%unit)
         self%owns_unit = .false.
      end if
      if (self%owns_debug_unit) then
         inquire (unit=self%debug_unit, opened=is_open)
         if (is_open) close (self%debug_unit)
         self%owns_debug_unit = .false.
      end if

      call self%timer%delete()

   end subroutine delete_context

   !> Return whether output at the given level should be emitted
   !>
   !> Mirrors the established `if (verbosity >= level .or. debug)` guard: a
   !> message tagged with `level` prints when the run verbosity reaches it, and
   !> the debug flag additionally unlocks the diagnostic band (level <= 3)
   !>
   !> @param[in]  level  verbosity level of the message
   pure function context_writes(self, level) result(do_write)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Verbosity level of the prospective message
      integer, intent(in) :: level
      !> Whether the message should be written
      logical :: do_write

      do_write = self%verbosity >= level .or. (self%debug .and. level <= 3)

   end function context_writes

   !> Maximum timing-tree depth to print at the current verbosity, which
   !> drives how much of the fully recorded tree is shown:
   !>   1 -> top-level phases only (depth 0),
   !>   2 -> + their components (depth 1),
   !>   3 -> + subcomponents (depth 2),
   !>   >= 4 -> everything (unbounded)
   pure function context_report_depth(self) result(d)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Deepest node level to print
      integer :: d

      if (self%verbosity >= 4) then
         d = huge(1)
      else
         d = max(0, self%verbosity - 1)
      end if

   end function context_report_depth

   !> Write a message to the context output unit if the level is enabled
   !>
   !> @param[in]  msg    message text
   !> @param[in]  level  verbosity level required to emit it (default 1)
   subroutine context_message(self, msg, level)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Message text
      character(len=*), intent(in) :: msg
      !> Verbosity level required to emit the message
      integer, intent(in), optional :: level
      !> Effective level
      integer :: lvl

      lvl = 1
      if (present(level)) lvl = level

      if (self%writes(lvl)) then
         !$omp critical(moist_context_io)
         write (self%unit, "(a)") msg
         !$omp end critical(moist_context_io)
      end if

   end subroutine context_message

   !> Write a debug message to the debug unit when debug is enabled
   !>
   !> Intended for diagnostic dumps that may originate inside parallel regions:
   !> the write is serialized with the same named critical as `message`
   !>
   !> @param[in]  msg  message text
   subroutine context_debug_message(self, msg)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Message text
      character(len=*), intent(in) :: msg

      if (self%debug) then
         !$omp critical(moist_context_io)
         write (self%debug_unit, "(a)") msg
         !$omp end critical(moist_context_io)
      end if

   end subroutine context_debug_message

   !> OpenMP thread count for the run
   !>
   !> The count fixed by `new_context`; a context that was never constructed
   !> reads the OpenMP environment, always >= 1
   function context_get_num_threads(self) result(nt)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Effective thread count (>= 1)
      integer :: nt

      nt = resolve_num_threads(self%nthreads)

   end function context_get_num_threads

   !> OpenMP team size for a kernel called with or without a context count
   !>
   !> @param[in] nthreads  explicit count; absent or <= 0 takes omp_get_max_threads
   function resolve_num_threads(nthreads) result(nt)
      !> Explicit thread count
      integer, intent(in), optional :: nthreads
      !> Effective thread count (>= 1)
      integer :: nt

      nt = 1
!$    nt = max(1, omp_get_max_threads())
      if (present(nthreads)) then
         if (nthreads > 0) nt = nthreads
      end if

   end function resolve_num_threads

   !> Render the run-wide settings through the pretty printer
   !>
   !> Runs single-threaded at the top of a run, so it needs no I/O critical
   !> section
   !>
   !> @param[in]  unit  optional output unit override (default `self%unit`)
   subroutine context_print_settings(self, unit)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Optional output unit override
      integer, intent(in), optional :: unit
      !> Pretty printer used to render the block
      type(prettyprinter) :: pp
      !> Effective output unit
      integer :: iu
      !> Effective thread count, materialized into a variable so the polymorphic
      !> `kv` sink sees a concrete integer (a function result is not conveyed)
      integer :: nthreads

      iu = self%unit
      if (present(unit)) iu = unit
      nthreads = self%get_num_threads()

      pp = new_prettyprinter(unit=iu, col_value=28)
      call pp%push("MOIST settings:")
      call pp%kv("Verbosity", self%verbosity)
      call pp%kv("Debug", self%debug)
      call pp%kv("Detailed profiling", self%do_profile)
      call pp%kv("OMP threads", nthreads)
      if (allocated(self%start_time)) call pp%kv("Start time", self%start_time)
      if (allocated(self%logfile)) then
         call pp%kv("Log file", self%logfile)
      else
         call pp%kv("Log file", "(stdout)")
      end if
      if (allocated(self%debugfile)) then
         call pp%kv("Debug file", self%debugfile)
      else
         call pp%kv("Debug file", "(same as log)")
      end if
      call pp%pop()
      call pp%blank()

   end subroutine context_print_settings

end module moist_context
