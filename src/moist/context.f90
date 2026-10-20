!> Shared run context for a moist calculation
!>
!> A single `moist_context_type` instance is constructed once at the top of a
!> run (the C-API handle wrapper or the CLI driver) and then *borrowed* by
!> every model and cavity through a non-owning pointer set at construction
!>
!> Bundles the run-wide settings that were previously duplicated on each
!> type:
!>
!>   * `verbosity` - the project-wide output level, see below,
!>   * `debug`     - the debug-stream flag; it raises `verbosity` to 6,
!>   * `timer`     - the profiling timer; because the context is shared by
!>     pointer, every borrower appends into this one timer tree; the host
!>     times its own work through the context's `start_timer`/`stop_timer`,
!>     under a top-level `Host` node, never through `timer%start`/`timer%stop`,
!>     which belong to moist's own nodes,
!>   * `cite`      - the citation labels of the parts that ran; like the timer,
!>     every borrower records into this one list, printed at the end of a run
!>
!> Verbosity levels; a host calls moist once per SCF iteration and geometry
!> step, so each level also bounds how often it prints:
!>
!>   0  silent          nothing (warnings always print)
!>   1  minimal         header block and results, once per run
!>   2  normal          key parameters once per run, a summary per geometry
!>   3  detailed        all parameters and per-geometry tables
!>   4  solver summary  one block per inner solve (DROP projection, ...)
!>   5  solver trace    per-iteration and per-point lines of inner solves
!>   6  debug           everything; the debug stream itself needs `debug`
!>
module moist_context
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use mctc_env, only: error_type, fatal_error
   use moist_utils_timer, only: timer_type, timer_name_len => name_len, cat_host
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_output_citations, only: citation, print_citations
!$ use omp_lib, only: omp_get_max_threads
   implicit none(type, external)
   private

   public :: moist_context_type, new_context, resolve_num_threads

   !> Top-level timer node that groups every host timer
   character(len=*), parameter :: host_group = "Host"

   !> Shared run context: owned at the top level, borrowed by every model/cavity
   type :: moist_context_type
      !> Output unit for informational messages (errors and warnings go to error_unit)
      integer :: unit = output_unit
      !> Output unit for debug dumps; defaults to `unit` (a separate `debugfile`
      !> splits the two streams)
      integer :: debug_unit = output_unit
      !> Verbosity level (0=silent ... 6=debug, see the module header); read
      !> through `verbosity()`
      integer, private :: level = 2
      !> Debug flag; enables the debug stream and raises the verbosity to 6
      logical :: debug = .false.
      !> Detailed-profiling flag, independent of the verbosity
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
      !> Open host timer nodes, the `Host` group first; a stack apart from the
      !> timer's own, so host and moist nodes never nest into each other
      integer, allocatable, private :: host_frames(:)
      !> Number of open host timer nodes, `Host` included
      integer, private :: nhost = 0
      !> Citation labels recorded by `cite`, in first-use order
      integer, allocatable, private :: cited(:)
      !> Number of recorded citation labels
      integer, private :: ncited = 0
   contains
      !> Tear down the owned timer and close any owned files
      procedure :: delete => delete_context
      !> Verbosity level of the run
      procedure, private :: context_get_verbosity
      !> Whether output at a given level prints
      procedure, private :: context_writes
      !> `verbosity()` returns the level, `verbosity(n)` whether level n prints
      generic :: verbosity => context_get_verbosity, context_writes
      !> Write a message to the context output unit, gated by verbosity
      procedure :: message => context_message
      !> Write a message like `message`, then flush the output unit
      procedure :: message_flush => context_message_flush
      !> Write a debug message to the debug unit when debug is enabled
      procedure :: debug_message => context_debug_message
      !> Write a warning to stderr, copied into the log when it is a file
      procedure :: warning => context_warning
      !> OpenMP thread count fixed at construction
      procedure :: get_num_threads => context_get_num_threads
      !> Render the run settings via the pretty printer
      procedure :: print_settings => context_print_settings
      !> Maximum timing-tree depth to print for the current verbosity
      procedure :: report_depth => context_report_depth
      !> Start a host timer under the `Host` node; a host uses this, never
      !> `timer%start`, which nests into moist's own nodes
      procedure :: start_timer => context_start_timer
      !> Stop the innermost open host timer by name
      procedure :: stop_timer => context_stop_timer
      !> Record that a part with this citation label ran
      procedure :: cite => context_cite
      !> Print the General citations plus those recorded by `cite`
      procedure :: print_citations => context_print_citations
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
   !> @param[in]  verbosity  output level (default 2)
   !> @param[in]  debug      diagnostic flag (default .false.), raises the verbosity to 6
   !> @param[in]  unit       already-open output unit (default output_unit)
   !> @param[in]  do_profile detailed-profiling flag (default .false.)
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
      !> Detailed-profiling flag (default .false.)
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

      if (present(verbosity)) self%level = verbosity
      if (present(debug)) self%debug = debug
      if (self%debug) self%level = max(self%level, 6)
      if (present(unit)) self%unit = unit
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

      ! Debug output defaults to the main unit
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

      ! Cite moist package
      call self%cite(citation%moist)

      ! The timer report follows the verbosity: silent at 0, deeper from 1 on
      call self%timer%new(verbose=self%level > 1, unit=self%unit, &
         & print_depth=merge(self%report_depth(), -1, self%level >= 1))

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
      self%nhost = 0
      if (allocated(self%host_frames)) deallocate (self%host_frames)
      self%ncited = 0
      if (allocated(self%cited)) deallocate (self%cited)

   end subroutine delete_context

   !> Verbosity level of the run, `ctx%verbosity()`
   pure function context_get_verbosity(self) result(level)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Verbosity level
      integer :: level

      level = self%level

   end function context_get_verbosity

   !> Return whether output at the given level should be emitted, `ctx%verbosity(level)`
   !>
   !> A message tagged with `level` prints when the run verbosity reaches it
   !>
   !> @param[in]  level  verbosity level of the message
   pure function context_writes(self, level) result(do_write)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Verbosity level of the prospective message
      integer, intent(in) :: level
      !> Whether the message should be written
      logical :: do_write

      do_write = self%level >= level

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

      if (self%level >= 4) then
         d = huge(1)
      else
         d = max(0, self%level - 1)
      end if

   end function context_report_depth

   !> Start a host timer: a node under the top-level `Host` node
   !>
   !> - host timers nest among themselves on their own stack: a moist timer
   !>   started meanwhile stays where it would be, and a host timer never
   !>   lands under a moist node
   !> - the first open host timer also starts `Host`, the last one closed
   !>   stops it, so `Host` times the spans with any host timer open
   !>
   !> @param[in]  name   timer name under the contract of `check_timer_name`
   !> @param[out] error  invalid name or uninitialized timer; nothing started
   subroutine context_start_timer(self, name, error)
      !> Context instance
      class(moist_context_type), intent(inout) :: self
      !> Timer name
      character(len=*), intent(in) :: name
      !> Error on an invalid name or an uninitialized timer
      type(error_type), allocatable, intent(out) :: error
      !> Parent node and node of the new timer
      integer :: parent, id

      call check_timer_name(name, error)
      if (allocated(error)) return

      if (self%nhost == 0) then
         parent = self%timer%resolve(host_group, 0, cat_host)
         if (parent == 0) then
            call fatal_error(error, "Context timer is not initialized")
            return
         end if
         call push_host_frame(self, parent)
         call self%timer%start(parent)
      else
         parent = self%host_frames(self%nhost)
      end if

      id = self%timer%resolve(name, parent)
      call push_host_frame(self, id)
      call self%timer%start(id)

   end subroutine context_start_timer

   !> Stop the innermost open host timer, which must carry `name`
   !>
   !> Closing the last open host timer also stops `Host`
   !>
   !> @param[in]  name   name of the innermost open host timer
   !> @param[out] error  invalid name, no host timer open or a different
   !>                    innermost one; nothing stopped
   subroutine context_stop_timer(self, name, error)
      !> Context instance
      class(moist_context_type), intent(inout) :: self
      !> Expected timer name
      character(len=*), intent(in) :: name
      !> Error on an invalid name or an unbalanced stop
      type(error_type), allocatable, intent(out) :: error
      !> Innermost open host node
      integer :: id

      call check_timer_name(name, error)
      if (allocated(error)) return
      if (self%nhost < 2) then
         call fatal_error(error, "No host timer is open")
         return
      end if
      id = self%host_frames(self%nhost)
      if (self%timer%node_name(id) /= name) then
         call fatal_error(error, "Host timer '"//name//"' is not the innermost open one, '"// &
            & trim(self%timer%node_name(id))//"'")
         return
      end if

      call self%timer%stop(id)
      self%nhost = self%nhost - 1
      if (self%nhost == 1) then
         call self%timer%stop(self%host_frames(1))
         self%nhost = 0
      end if

   end subroutine context_stop_timer

   !> Check a host timer name against the contract shared by every language
   !>
   !> - 1 to `name_len` bytes of UTF-8: an ASCII character takes one byte,
   !>   `é` two
   !> - no control character (U+0000 to U+001F, U+007F to U+009F)
   !> - no leading or trailing space: every byte is part of the name, so a
   !>   padded Fortran buffer is passed as `trim(buffer)`
   !>
   !> @param[in]  name   candidate name
   !> @param[out] error  the first rule the name breaks
   subroutine check_timer_name(name, error)
      !> Candidate name
      character(len=*), intent(in) :: name
      !> Error naming the broken rule
      type(error_type), allocatable, intent(out) :: error
      !> Byte position, sequence length and continuation index
      integer :: i, nbytes, k
      !> Current byte and decoded code point
      integer :: byte, code
      !> Byte limit, as text for the message
      character(len=8) :: limit

      if (len(name) == 0) then
         call fatal_error(error, "Timer name is empty")
         return
      end if
      if (len(name) > timer_name_len) then
         write (limit, "(i0)") timer_name_len
         call fatal_error(error, "Timer name exceeds "//trim(limit)//" bytes")
         return
      end if
      if (name(1:1) == " " .or. name(len(name):len(name)) == " ") then
         call fatal_error(error, "Timer name starts or ends with a space")
         return
      end if

      i = 1
      do while (i <= len(name))
         byte = ichar(name(i:i))
         select case (byte)
         case (0:127)
            nbytes = 1
            code = byte
         case (194:223)
            nbytes = 2
            code = iand(byte, 31)
         case (224:239)
            nbytes = 3
            code = iand(byte, 15)
         case (240:244)
            nbytes = 4
            code = iand(byte, 7)
         case default
            nbytes = 0
         end select
         if (nbytes == 0 .or. i + nbytes - 1 > len(name)) exit
         do k = i + 1, i + nbytes - 1
            byte = ichar(name(k:k))
            if (iand(byte, 192) /= 128) exit
            code = 64*code + iand(byte, 63)
         end do
         if (k < i + nbytes) exit
         ! overlong forms, UTF-16 surrogates and code points past U+10FFFF
         if ((nbytes == 3 .and. code < 2048) .or. (nbytes == 4 .and. code < 65536) &
            & .or. (code >= 55296 .and. code <= 57343) .or. code > 1114111) exit
         if (code < 32 .or. (code >= 127 .and. code <= 159)) then
            call fatal_error(error, "Timer name contains a control character")
            return
         end if
         i = i + nbytes
      end do
      if (i <= len(name)) call fatal_error(error, "Timer name is not valid UTF-8")

   end subroutine check_timer_name

   !> Push a node onto the host timer stack, growing it on demand
   subroutine push_host_frame(self, id)
      !> Context instance
      class(moist_context_type), intent(inout) :: self
      !> Node to push
      integer, intent(in) :: id
      !> Grown stack
      integer, allocatable :: frames(:)

      if (.not. allocated(self%host_frames)) allocate (self%host_frames(8))
      if (self%nhost == size(self%host_frames)) then
         allocate (frames(2*self%nhost))
         frames(:self%nhost) = self%host_frames(:self%nhost)
         call move_alloc(frames, self%host_frames)
      end if
      self%nhost = self%nhost + 1
      self%host_frames(self%nhost) = id

   end subroutine push_host_frame

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

      if (self%verbosity(lvl)) then
         !$omp critical(moist_context_io)
         write (self%unit, "(a)") msg
         !$omp end critical(moist_context_io)
      end if

   end subroutine context_message

   !> Write a message like `message`, then flush the output unit
   !>
   !> Write and flush share one critical section, so no other message lands
   !> between them
   !>
   !> @param[in]  msg    message text
   !> @param[in]  level  verbosity level required to emit it (default 1)
   subroutine context_message_flush(self, msg, level)
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

      if (self%verbosity(lvl)) then
         !$omp critical(moist_context_io)
         write (self%unit, "(a)") msg
         flush (self%unit)
         !$omp end critical(moist_context_io)
      end if

   end subroutine context_message_flush

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

   !> Write a warning to stderr, copied into the log when it is a file
   !>
   !> Not gated by verbosity; a terminal run shows the warning once on stderr,
   !> an owned logfile or a host-supplied unit receives a second copy in place
   !>
   !> Both writes flush inside the same critical as `message`, so the warning
   !> lands before a later crash
   !>
   !> @param[in]  msg  warning text, without the `[Warning]` prefix
   subroutine context_warning(self, msg)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Warning text
      character(len=*), intent(in) :: msg

      !$omp critical(moist_context_io)
      write (error_unit, '("[Warning]", 1x, a)') msg
      flush (error_unit)
      if (self%unit /= output_unit .and. self%unit /= error_unit) then
         write (self%unit, '("[Warning]", 1x, a)') msg
         flush (self%unit)
      end if
      !$omp end critical(moist_context_io)

   end subroutine context_warning

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
      call pp%kv("Verbosity", self%level)
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

   !> Record a citation label, once; parts call this from `update`, after
   !> `require_context`, with one of the `citation_*` labels of
   !> `moist_output_citations`
   !>
   !> @param[in]  label  registry label of the part that ran
   subroutine context_cite(self, label)
      !> Context instance
      class(moist_context_type), intent(inout) :: self
      !> Registry label
      integer, intent(in) :: label
      !> Grown label list
      integer, allocatable :: cited(:)

      if (self%ncited > 0) then
         if (any(self%cited(:self%ncited) == label)) return
      end if

      if (.not. allocated(self%cited)) allocate (self%cited(8))
      if (self%ncited == size(self%cited)) then
         allocate (cited(2*self%ncited))
         cited(:self%ncited) = self%cited(:self%ncited)
         call move_alloc(cited, self%cited)
      end if
      self%ncited = self%ncited + 1
      self%cited(self%ncited) = label

   end subroutine context_cite

   !> Print the General citations plus those recorded by `cite`; silent at
   !> verbosity 0
   !>
   !> @param[in]  unit  optional output unit override (default `self%unit`)
   subroutine context_print_citations(self, unit)
      !> Context instance
      class(moist_context_type), intent(in) :: self
      !> Optional output unit override
      integer, intent(in), optional :: unit
      !> Effective output unit
      integer :: iu

      if (self%level < 1) return
      iu = self%unit
      if (present(unit)) iu = unit

      if (allocated(self%cited)) then
         call print_citations(iu, self%cited(:self%ncited))
      else
         call print_citations(iu, [integer ::])
      end if

   end subroutine context_print_citations

end module moist_context
