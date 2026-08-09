!> Dynamic, nesting-based profiling timer
!>
!> Timers are addressed by name, never by a pre-allocated integer slot:
!> The registry grows on demand, so callers do not need to know the number of
!> steps up front
!> A timer is identified by the pair (parent, name), where the parent is
!> whichever timer is open when `start` is called
!> Re-entering a name (e.g. inside a loop) accumulates into the same node;
!> arbitrary nesting depth (sub- and subsub-timers) falls out of the runtime
!> nesting with no extra bookkeeping
!>
!> Two call styles share the same node storage:
!>   * name-based `start("X")` / `stop()` (optionally `stop("X")` to assert the
!>     name) uses an internal stack to derive the parent and is convenient for
!>     straight-line code.
!>   * handle-based `start(id)` / `stop(id)`, with `id` obtained once from
!>     `resolve`, skips the hash lookup and the stack; use it inside hot loops.
!>
!> Any timer left running at report time, or closed by a mismatched
!> `stop("name")`, is poisoned and rendered as NaN
!>
!> FIXME: TODO: the timer is NOT thread-safe. Inside OpenMP regions gate start/stop to
!> a single thread (see cavity/drop/gradient.f90).
module moist_utils_timer
   use mctc_env, only: wp, int64 => i8
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none(type, external)
   private

   public :: timer_type
   public :: cat_none, cat_setup, cat_energy, cat_potential, cat_gradient, &
             cat_hessian, cat_solve, cat_properties, cat_io

   !> Timer categories: an orthogonal tag summed across the module tree at
   !> report time. A node inherits its parent's category unless one is given.
   integer, parameter :: cat_none = 0
   integer, parameter :: cat_setup = 1
   integer, parameter :: cat_energy = 2
   integer, parameter :: cat_potential = 3
   integer, parameter :: cat_gradient = 4
   integer, parameter :: cat_hessian = 5
   integer, parameter :: cat_solve = 6
   integer, parameter :: cat_properties = 7
   integer, parameter :: cat_io = 8

   !> Number of named categories
   integer, parameter :: ncat = 8

   !> Human-readable category labels (index 0 = uncategorized, never in pivot)
   character(len=12), parameter :: cat_label(0:ncat) = [character(len=12) :: &
                                                        "n/a", "Setup", "Energy", "Potential", &
                                                        "Gradient", "Hessian", "Solve", &
                                                        "Properties", "I/O"]

   !> Maximum stored length of a timer name
   integer, parameter :: name_len = 32

   !> Initial node capacity and hash-table size (both grow on demand)
   integer, parameter :: init_cap = 32
   integer, parameter :: init_tab = 64

   !> Report table layout. `sect_w` is the display width (columns, not bytes) of
   !> the leftmost "section" column; the numeric tail (time/percent/calls) is a
   !> fixed 32 columns, so the whole rule is `sect_w + 32`.
   integer, parameter :: sect_w = 26
   integer, parameter :: tbl_w = sect_w + 32

   !> Unicode (UTF-8) box-drawing guides for the tree column. NOTE: intentionally
   !> non-ASCII, at explicit request, so this is the one place the module departs
   !> from the ASCII-only house rule. Each glyph is 3 bytes but a single display
   !> column; `disp_len` measures columns (not bytes) so the numeric columns to
   !> their right stay aligned regardless of nesting depth.
   !>   tg_branch = "|--" (a non-last child), tg_last = "`--" (the last child),
   !>   tg_pipe   = "|  " (an ancestor guide), tg_blank = "   " (a spent guide).
   character(len=*), parameter :: tg_branch = "├─ "
   character(len=*), parameter :: tg_last = "└─ "
   character(len=*), parameter :: tg_pipe = "│  "
   character(len=*), parameter :: tg_blank = "  "
   !> Bare vertical guide drawn on a group-separator blank line so the tree stays
   !> visually connected across the gap (aligns with the child connectors).
   character(len=*), parameter :: tg_vert = "│"

   !> FNV-1a (32-bit) constants; arithmetic is masked to 32 bits each step so
   !> the intermediate product stays well inside a signed 64-bit integer.
   integer(int64), parameter :: fnv_offset = 2166136261_int64
   integer(int64), parameter :: fnv_prime = 16777619_int64
   integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)

   !> Hierarchical profiling timer.
   type :: timer_type

      !> Number of registered nodes
      integer, private :: n = 0
      !> Allocated node capacity
      integer, private :: cap = 0
      !> Whether the timer is initialized and measuring
      logical, private :: active = .false.
      !> Whether the report prints the verbose (wall+cpu) footer
      logical, private :: verbose = .false.
      !> Whether the report separates the top-level groups with a blank line
      logical, private :: group_blanks = .true.

      !> Node name
      character(len=name_len), private, allocatable :: names(:)
      !> Parent node index (0 = top level)
      integer, private, allocatable :: parent(:)
      !> Category tag (cat_none = none)
      integer, private, allocatable :: category(:)
      !> Number of start calls
      integer, private, allocatable :: ncalls(:)
      !> Accumulated wall time in clock ticks
      integer(int64), private, allocatable :: t_acc(:)
      !> Start tick while running
      integer(int64), private, allocatable :: t0(:)
      !> Whether the node is currently open
      logical, private, allocatable :: running(:)
      !> Whether the node was left unbalanced (reported as NaN)
      logical, private, allocatable :: poisoned(:)

      !> Open-addressing hash table (values are node indices, 0 = empty)
      integer, private, allocatable :: htab(:)
      !> Nesting stack of open node indices
      integer, private, allocatable :: stack(:)
      !> Stack pointer (number of open frames)
      integer, private :: sp = 0

      !> Wall-clock tick at root start
      integer(int64), private :: t_root0 = 0
      !> Clock ticks per second
      integer(int64), private :: rate = 1
      !> CPU time at root start
      real(wp), private :: cpu_root0 = 0.0_wp

   contains

      procedure :: new => new_timer
      procedure :: delete => delete_timer
      procedure :: reset => reset_timer
      procedure :: set_verbose
      procedure :: set_group_blanks
      procedure :: resolve
      procedure, private :: start_name
      procedure, private :: start_id
      generic :: start => start_name, start_id
      procedure, private :: stop_top
      procedure, private :: stop_name
      procedure, private :: stop_id
      generic :: stop => stop_top, stop_name, stop_id
      procedure, private :: get_total
      procedure, private :: get_name
      generic :: get => get_total, get_name
      procedure :: write => write_report
      procedure :: write_timing
      ! read-only introspection over the registered tree (for programmatic
      ! consumers that enumerate whatever was measured, e.g. benchmarks)
      procedure :: num_nodes
      procedure :: node_name
      procedure :: node_depth
      procedure :: node_time
      !> Currently-open node (top of the nesting stack), 0 if none. Lets a
      !> handle-based hot loop attach its fine timers under whatever node the
      !> caller opened by name, instead of guessing a fixed parent.
      procedure :: current
      !> Current nesting depth (number of open frames). Pair with `unwind` to
      !> keep the shared stack balanced across early error returns.
      procedure :: current_depth
      !> Stop and pop open frames until the stack is back at a saved depth.
      procedure :: unwind

   end type timer_type

contains

!> Initialize the timer, allocating the growable registry and stamping the
!> root wall/CPU clocks.
!> @param[inout] self     Timer instance
!> @param[in]    verbose  Optional: print the verbose report footer
   subroutine new_timer(self, verbose)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Optional verbose flag
      logical, intent(in), optional :: verbose

      call self%delete()

      self%verbose = .false.
      if (present(verbose)) self%verbose = verbose

      self%cap = init_cap
      allocate (character(len=name_len) :: self%names(self%cap))
      allocate (self%parent(self%cap), source=0)
      allocate (self%category(self%cap), source=cat_none)
      allocate (self%ncalls(self%cap), source=0)
      allocate (self%t_acc(self%cap), source=0_int64)
      allocate (self%t0(self%cap), source=0_int64)
      allocate (self%running(self%cap), source=.false.)
      allocate (self%poisoned(self%cap), source=.false.)
      allocate (self%htab(init_tab), source=0)
      allocate (self%stack(init_cap), source=0)

      self%n = 0
      self%sp = 0
      self%active = .true.

      call system_clock(self%t_root0, self%rate)
      call cpu_time(self%cpu_root0)

   end subroutine new_timer

!> Release all timer storage.
!> @param[inout] self  Timer instance
   subroutine delete_timer(self)
      !> Timer instance
      class(timer_type), intent(inout) :: self

      self%n = 0
      self%cap = 0
      self%sp = 0
      self%active = .false.
      self%verbose = .false.
      if (allocated(self%names)) deallocate (self%names)
      if (allocated(self%parent)) deallocate (self%parent)
      if (allocated(self%category)) deallocate (self%category)
      if (allocated(self%ncalls)) deallocate (self%ncalls)
      if (allocated(self%t_acc)) deallocate (self%t_acc)
      if (allocated(self%t0)) deallocate (self%t0)
      if (allocated(self%running)) deallocate (self%running)
      if (allocated(self%poisoned)) deallocate (self%poisoned)
      if (allocated(self%htab)) deallocate (self%htab)
      if (allocated(self%stack)) deallocate (self%stack)

   end subroutine delete_timer

!> Zero all accumulated times while keeping the registered node tree, and
!> restart the root clocks. Used for per-window reporting: `write` then `reset`
!> begins a fresh window without reallocating.
!> @param[inout] self  Timer instance
   subroutine reset_timer(self)
      !> Timer instance
      class(timer_type), intent(inout) :: self

      if (.not. self%active) return

      if (self%n > 0) then
         self%t_acc(1:self%n) = 0_int64
         self%t0(1:self%n) = 0_int64
         self%ncalls(1:self%n) = 0
         self%running(1:self%n) = .false.
         self%poisoned(1:self%n) = .false.
      end if
      self%sp = 0

      call system_clock(self%t_root0)
      call cpu_time(self%cpu_root0)

   end subroutine reset_timer

!> Set the verbose report flag without touching accumulated timings.
!> @param[inout] self  Timer instance
!> @param[in]    flag  New verbose flag
   subroutine set_verbose(self, flag)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> New verbose flag
      logical, intent(in) :: flag

      self%verbose = flag

   end subroutine set_verbose

!> Set the blank-line grouping flag: when on, the report separates the top-level
!> groups with a blank line. Does not touch accumulated timings.
!> @param[inout] self  Timer instance
!> @param[in]    flag  New grouping flag
   subroutine set_group_blanks(self, flag)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> New grouping flag
      logical, intent(in) :: flag

      self%group_blanks = flag

   end subroutine set_group_blanks

!> Return a stable node handle for (name, parent), creating the node if needed.
!> Use the handle with the `start(id)`/`stop(id)` fast path inside hot loops.
!> @param[inout] self      Timer instance
!> @param[in]    name      Timer name
!> @param[in]    parent    Parent node handle (0 = top level)
!> @param[in]    category  Optional category tag (set on first creation)
!> @return       Node handle (>= 1)
   function resolve(self, name, parent, category) result(id)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Timer name
      character(len=*), intent(in) :: name
      !> Parent handle
      integer, intent(in) :: parent
      !> Optional category
      integer, intent(in), optional :: category
      !> Node handle
      integer :: id

      integer(int64) :: h
      integer :: slot

      if (.not. self%active) then
         id = 0
         return
      end if

      h = hash_key(name, parent)
      slot = int(iand(h, int(size(self%htab) - 1, int64)), kind(slot)) + 1

      ! probe for an existing node
      do
         id = self%htab(slot)
         if (id == 0) exit
         if (self%parent(id) == parent .and. self%names(id) == pad_name(name)) then
            if (present(category)) then
               if (self%category(id) == cat_none) self%category(id) = category
            end if
            return
         end if
         slot = slot + 1
         if (slot > size(self%htab)) slot = 1
      end do

      ! create a new node
      call ensure_node_capacity(self)
      self%n = self%n + 1
      id = self%n
      self%names(id) = name
      self%parent(id) = parent
      ! Category precedence: an explicit tag wins; otherwise inherit the parent's
      ! (so the handle path `resolve`+`start_id` categorizes the same way as the
      ! name/stack path, and fine-grained timers are not dropped from the pivot).
      if (present(category)) then
         self%category(id) = category
      else if (parent > 0) then
         self%category(id) = self%category(parent)
      else
         self%category(id) = cat_none
      end if
      self%t_acc(id) = 0_int64
      self%t0(id) = 0_int64
      self%running(id) = .false.
      self%poisoned(id) = .false.
      self%ncalls(id) = 0

      ! insert into the hash table, growing it first if it would be too full
      if (self%n*10 >= size(self%htab)*7) then
         call grow_table(self)
      else
         call table_insert(self, id)
      end if

   end function resolve

!> Start a named timer under the currently open timer (its parent). Pushes the
!> node onto the nesting stack.
!> @param[inout] self      Timer instance
!> @param[in]    name      Timer name
!> @param[in]    category  Optional category tag
   subroutine start_name(self, name, category)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Timer name
      character(len=*), intent(in) :: name
      !> Optional category
      integer, intent(in), optional :: category

      integer :: par, id

      if (.not. self%active) return

      if (self%sp > 0) then
         par = self%stack(self%sp)
      else
         par = 0
      end if

      id = self%resolve(name, par, category)

      ! inherit the parent's category unless one was explicitly given
      if (self%category(id) == cat_none .and. par > 0) self%category(id) = self%category(par)

      call ensure_stack_capacity(self)
      self%sp = self%sp + 1
      self%stack(self%sp) = id

      call start_id(self, id)

   end subroutine start_name

!> Start a timer by handle (hot-path fast form, no lookup, no stack).
!> @param[inout] self  Timer instance
!> @param[in]    id    Node handle from `resolve`
   subroutine start_id(self, id)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Node handle
      integer, intent(in) :: id

      integer(int64) :: c

      if (.not. self%active .or. id < 1 .or. id > self%n) return

      ! A start on an already-running node is an imbalance (a missing stop):
      ! poison it so the interval surfaces as NaN instead of being silently
      ! dropped (the closing stop would only measure from this last start).
      if (self%running(id)) self%poisoned(id) = .true.

      call system_clock(c)
      self%t0(id) = c
      self%running(id) = .true.
      self%ncalls(id) = self%ncalls(id) + 1

   end subroutine start_id

!> Stop the innermost open timer (top of the nesting stack).
!> @param[inout] self  Timer instance
   subroutine stop_top(self)
      !> Timer instance
      class(timer_type), intent(inout) :: self

      integer :: id

      if (.not. self%active) return
      if (self%sp < 1) return

      id = self%stack(self%sp)
      self%sp = self%sp - 1
      call stop_id(self, id)

   end subroutine stop_top

!> Stop a named timer, asserting it matches an open frame. If inner frames were
!> left open (e.g. an early return), they are unwound and poisoned as NaN so the
!> imbalance is visible and the stack recovers.
!> @param[inout] self  Timer instance
!> @param[in]    name  Expected timer name
   subroutine stop_name(self, name)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Expected timer name
      character(len=*), intent(in) :: name

      integer :: k, found, id

      if (.not. self%active) return

      ! locate the matching frame, searching from the top down
      found = 0
      do k = self%sp, 1, -1
         if (self%names(self%stack(k)) == pad_name(name)) then
            found = k
            exit
         end if
      end do
      if (found == 0) return

      ! poison any skipped inner frames
      do k = self%sp, found + 1, -1
         self%poisoned(self%stack(k)) = .true.
         self%running(self%stack(k)) = .false.
      end do

      id = self%stack(found)
      self%sp = found - 1
      call stop_id(self, id)

   end subroutine stop_name

!> Stop a timer by handle (hot-path fast form). A stop without a matching start
!> poisons the node.
!> @param[inout] self  Timer instance
!> @param[in]    id    Node handle
   subroutine stop_id(self, id)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Node handle
      integer, intent(in) :: id

      integer(int64) :: c

      if (.not. self%active .or. id < 1 .or. id > self%n) return

      call system_clock(c)
      if (self%running(id)) then
         self%t_acc(id) = self%t_acc(id) + (c - self%t0(id))
         self%running(id) = .false.
      else
         self%poisoned(id) = .true.
      end if

   end subroutine stop_id

!> Overall elapsed wall time since `new`/`reset`, in seconds.
!> @param[inout] self  Timer instance
!> @return       Elapsed seconds
   function get_total(self) result(sec)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Elapsed seconds
      real(wp) :: sec

      integer(int64) :: c

      if (.not. self%active) then
         sec = 0.0_wp
         return
      end if
      call system_clock(c)
      sec = real(c - self%t_root0, wp)/real(self%rate, wp)

   end function get_total

!> Accumulated wall time of a named timer, in seconds. The name may be a path
!> "Parent/Child/..." to disambiguate a leaf name that occurs under several
!> parents. Returns 0 for an unknown path and NaN for a poisoned/open timer.
!> @param[inout] self  Timer instance
!> @param[in]    path  Timer name or slash-separated path
!> @return       Elapsed seconds (NaN if unbalanced)
   function get_name(self, path) result(sec)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Timer name or path
      character(len=*), intent(in) :: path
      !> Elapsed seconds
      real(wp) :: sec

      integer :: id

      id = find_path(self, path)
      if (id == 0) then
         sec = 0.0_wp
      else
         sec = node_seconds(self, id)
      end if

   end function get_name

!> Number of registered timer nodes. Nodes are indexed 1..num_nodes in
!> registration (pre-order) order: a parent is always registered before its
!> children, so iterating 1..num_nodes yields a valid tree traversal.
!> @param[in] self  Timer instance
!> @return    Node count
   pure function num_nodes(self) result(n)
      !> Timer instance
      class(timer_type), intent(in) :: self
      !> Node count
      integer :: n

      n = self%n

   end function num_nodes

!> Return the currently-open node (top of the nesting stack), or 0 if the stack
!> is empty or the timer is inactive.
!> @param[in] self  Timer instance
!> @return    Node handle of the innermost open timer, else 0
   pure function current(self) result(id)
      !> Timer instance
      class(timer_type), intent(in) :: self
      !> Handle of the innermost open node
      integer :: id

      id = 0
      if (self%active .and. self%sp > 0) id = self%stack(self%sp)

   end function current

!> Current nesting depth (number of open frames on the stack), 0 if inactive.
!>
!> Capture this at the top of a timed routine and pass it to `unwind` on an early
!> (error) return so any frames the routine opened are closed and the run-wide
!> stack is left balanced for whatever runs next.
!> @param[in] self  Timer instance
!> @return    Number of currently-open frames
   pure function current_depth(self) result(d)
      !> Timer instance
      class(timer_type), intent(in) :: self
      !> Open-frame count
      integer :: d

      d = 0
      if (self%active) d = self%sp

   end function current_depth

!> Close open frames until the nesting stack is back at `depth`.
!>
!> Each unwound frame is stopped through the normal `stop_top` path, so it
!> records the partial time measured up to the unwind point. Intended for error
!> paths that `return` between a `start` and its matching `stop`; a no-op when the
!> stack is already at (or below) `depth`.
!> @param[inout] self   Timer instance
!> @param[in]    depth  Target open-frame count to unwind back to
   subroutine unwind(self, depth)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Target depth (as returned by `current_depth`)
      integer, intent(in) :: depth

      if (.not. self%active) return
      do while (self%sp > max(0, depth))
         call stop_top(self)
      end do

   end subroutine unwind

!> Name of a node by index.
!> @param[in] self  Timer instance
!> @param[in] id    Node index (1..num_nodes)
!> @return    Trimmed node name
   pure function node_name(self, id) result(name)
      !> Timer instance
      class(timer_type), intent(in) :: self
      !> Node index
      integer, intent(in) :: id
      !> Node name
      character(len=:), allocatable :: name

      if (id < 1 .or. id > self%n) then
         name = ""
      else
         name = trim(self%names(id))
      end if

   end function node_name

!> Nesting depth of a node (0 = top level), from its parent chain.
!> @param[in] self  Timer instance
!> @param[in] id    Node index (1..num_nodes)
!> @return    Depth
   pure function node_depth(self, id) result(depth)
      !> Timer instance
      class(timer_type), intent(in) :: self
      !> Node index
      integer, intent(in) :: id
      !> Depth
      integer :: depth

      integer :: cur

      depth = 0
      if (id < 1 .or. id > self%n) return
      cur = self%parent(id)
      do while (cur > 0)
         depth = depth + 1
         cur = self%parent(cur)
      end do

   end function node_depth

!> Accumulated seconds of a node by index (NaN if unbalanced).
!> @param[in] self  Timer instance
!> @param[in] id    Node index (1..num_nodes)
!> @return    Seconds
   pure function node_time(self, id) result(sec)
      !> Timer instance
      class(timer_type), intent(in) :: self
      !> Node index
      integer, intent(in) :: id
      !> Seconds
      real(wp) :: sec

      sec = node_seconds(self, id)

   end function node_time

!> Write a single named timer line to `iunit`.
!> @param[inout] self     Timer instance
!> @param[in]    iunit    Output unit
!> @param[in]    path     Timer name or path
!> @param[in]    verbose  Unused; kept for call compatibility
   subroutine write_timing(self, iunit, path, verbose)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Output unit
      integer, intent(in) :: iunit
      !> Timer name or path
      character(len=*), intent(in) :: path
      !> Unused verbose flag
      logical, intent(in), optional :: verbose

      real(wp) :: sec

      sec = self%get(path)
      write (iunit, '(1x,a,1x,"...",f12.4," sec")') trim(path), sec

   end subroutine write_timing

!> Write the hierarchical timing report: a module tree with per-node share and
!> call counts, an "(other)" remainder per parent, a category-pivot summary, and
!> a wall/CPU footer. Any timer still open is poisoned and shown as NaN.
!>
!> The tree is drawn with Unicode box-drawing guides and a fixed, header-labelled
!> column layout; percentages are expressed against the run total ("%tot") so
!> every row shares one stable reference and two runs diff cleanly.
!> @param[inout] self       Timer instance
!> @param[in]    iunit      Output unit
!> @param[in]    inmsg      Optional label for the total row
!> @param[in]    max_depth  Optional deepest node level to print (0 = top-level
!>                          only). Deeper nodes are folded into their ancestor's
!>                          total. Absent = print the whole tree.
   subroutine write_report(self, iunit, inmsg, max_depth)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Output unit
      integer, intent(in) :: iunit
      !> Optional label
      character(len=*), intent(in), optional :: inmsg
      !> Optional maximum print depth
      integer, intent(in), optional :: max_depth

      real(wp) :: total_sec, total_cpu, cpu_now, cat_sec(ncat), pct
      integer :: i, c, depth_cap
      logical :: first_root, prev_expanded
      character(len=32) :: numbuf
      character(len=26) :: msg

      if (.not. self%active) return

      depth_cap = huge(1)
      if (present(max_depth)) depth_cap = max_depth

      ! poison anything still open: reaching a report with an open timer is a bug
      do i = 1, self%n
         if (self%running(i)) self%poisoned(i) = .true.
      end do

      total_sec = self%get()
      call cpu_time(cpu_now)
      total_cpu = cpu_now - self%cpu_root0

      if (present(inmsg)) then
         msg = inmsg//" timings"
      else
         msg = "total time"
      end if

      write (iunit, "(a)")
      write (iunit, "(1x,a)") banner_line()
      write (iunit, "(a)")
      ! column header, then a rule the full table width
      write (iunit, "(1x,a)") header_line()
      write (iunit, "(1x,a)") repeat("-", tbl_w)

      ! module tree, top-level nodes first (each drawn as a bare root; separate
      ! multiple roots with a blank line when grouping is on)
      first_root = .true.
      prev_expanded = .false.
      do i = 1, self%n
         if (self%parent(i) /= 0) cycle
         if (self%ncalls(i) == 0 .and. .not. self%poisoned(i)) cycle
         ! separate top-level groups with a blank line, but only after a group that
         ! actually expanded into printed subcategories, so runs of bare one-line
         ! groups stay compact
         if (self%group_blanks .and. .not. first_root .and. prev_expanded) then
            write (iunit, "(a)")
         end if
         call print_node(self, iunit, i, "", 0, .true., total_sec, depth_cap)
         first_root = .false.
         prev_expanded = node_has_printed_children(self, i, 0, depth_cap)
      end do

      ! category pivot: sum each maximal category region once (a node whose
      ! parent has a different category is the head of that region)
      cat_sec = 0.0_wp
      do i = 1, self%n
         if (self%category(i) == cat_none) cycle
         if (self%poisoned(i) .or. self%ncalls(i) == 0) cycle
         if (self%parent(i) == 0) then
            cat_sec(self%category(i)) = cat_sec(self%category(i)) + node_seconds(self, i)
         else if (self%category(self%parent(i)) /= self%category(i)) then
            ! Head of a category region nested inside a differently-categorized
            ! one. Count it under its own category ...
            cat_sec(self%category(i)) = cat_sec(self%category(i)) + node_seconds(self, i)
            ! ... but the enclosing region's head already counted this child's
            ! time inclusively under the parent's category, so remove the overlap
            ! to keep the totals a partition of wall time (no double counting).
            if (self%category(self%parent(i)) /= cat_none) then
               cat_sec(self%category(self%parent(i))) = &
                  cat_sec(self%category(self%parent(i))) - node_seconds(self, i)
            end if
         end if
      end do

      if (any(cat_sec > 0.0_wp)) then
         write (iunit, "(1x,a)") repeat("-", tbl_w)
         write (iunit, "(1x,a)") "by category"
         do c = 1, ncat
            if (cat_sec(c) <= 0.0_wp) cycle
            pct = 0.0_wp
            if (total_sec > 0.0_wp) pct = 100.0_wp*cat_sec(c)/total_sec
            write (numbuf, "(f13.4,f9.2)") cat_sec(c), pct
            write (iunit, "(1x,a)") pad_disp(trim(cat_label(c)), sect_w)//numbuf(1:22)
         end do
      end if

      write (iunit, "(1x,a)") repeat("=", tbl_w)
      if (self%verbose) then
         write (iunit, '(2x,"Wall-time",4x,f14.4," sec")') total_sec
         if (total_sec > 0.0_wp) then
            write (iunit, '(2x,"CPU-time",5x,f14.4," sec",2x,"(",f0.1,"x)")') &
               total_cpu, total_cpu/total_sec
         else
            write (iunit, '(2x,"CPU-time",5x,f14.4," sec")') total_cpu
         end if
      else
         write (iunit, '(2x,a26,f14.4," sec")') msg, total_sec
      end if
      write (iunit, "(1x,a)") repeat("=", tbl_w)
      write (iunit, "(a)")

   end subroutine write_report

!> Print one node and, recursively, its children, plus an "(other)" remainder.
!>
!> The section column carries the Unicode tree guides: `prefix` is the accumulated
!> ancestor guide string ("|  "/"   " per level), and `is_last` selects this node's
!> own connector ("`--" vs "|--"). A depth-0 node is the tree root and is drawn
!> bare (no connector). Percentages are against the run total (`total_sec`).
!> @param[in] self       Timer instance
!> @param[in] iunit      Output unit
!> @param[in] id         Node index
!> @param[in] prefix     Accumulated ancestor guide string for this node's row
!> @param[in] depth      Nesting depth (0 = tree root, drawn without a connector)
!> @param[in] is_last    Whether this node is the last among its printed siblings
!> @param[in] total_sec  Run total, the shared reference for every "%tot"
!> @param[in] max_depth  Deepest level to print; children below it are folded in
!>
!> Whether the top-level groups are separated by a blank line is read from the
!> timer's own `group_blanks` flag (set via `set_group_blanks`).
   recursive subroutine print_node(self, iunit, id, prefix, depth, is_last, &
      & total_sec, max_depth)
      !> Timer instance
      type(timer_type), intent(in) :: self
      !> Output unit
      integer, intent(in) :: iunit
      !> Node index
      integer, intent(in) :: id
      !> Accumulated ancestor guide string
      character(len=*), intent(in) :: prefix
      !> Nesting depth (0 = root)
      integer, intent(in) :: depth
      !> Whether this is the last printed sibling
      logical, intent(in) :: is_last
      !> Run total for the percentage
      real(wp), intent(in) :: total_sec
      !> Deepest level to print; children below it are folded into this node
      integer, intent(in) :: max_depth

      real(wp) :: sec, pct, child_sum, rem
      integer :: j, k, nkids
      integer, allocatable :: kids(:)
      logical :: show_other, kid_last, prev_expanded
      character(len=:), allocatable :: sect, child_prefix
      character(len=32) :: numbuf

      ! skip nodes that were never measured (created via resolve but unused)
      if (self%ncalls(id) == 0 .and. .not. self%poisoned(id)) return

      sec = node_seconds(self, id)
      pct = 0.0_wp
      if (total_sec > 0.0_wp) pct = 100.0_wp*sec/total_sec

      ! the section cell: root bare, otherwise ancestor guides + this connector
      if (depth == 0) then
         sect = trim(self%names(id))
      else if (is_last) then
         sect = prefix//tg_last//trim(self%names(id))
      else
         sect = prefix//tg_branch//trim(self%names(id))
      end if

      write (numbuf, '(f13.4,f9.2,i9,"x")') sec, pct, self%ncalls(id)
      write (iunit, "(1x,a)") pad_disp(sect, sect_w)//numbuf(1:32)

      ! at the depth cap, fold all descendants into this node (show total only)
      if (depth >= max_depth) return

      ! gather printable children (registration order) and their summed time
      allocate (kids(0))
      child_sum = 0.0_wp
      do j = 1, self%n
         if (self%parent(j) /= id) cycle
         if (self%ncalls(j) == 0 .and. .not. self%poisoned(j)) cycle
         kids = [kids, j]
         if (.not. self%poisoned(j)) child_sum = child_sum + node_seconds(self, j)
      end do
      nkids = size(kids)

      ! is there time not attributed to any child? (drawn as a final "(other)")
      show_other = .false.
      rem = 0.0_wp
      if (nkids > 0 .and. .not. self%poisoned(id)) then
         rem = sec - child_sum
         ! only surface the remainder when it is a meaningful slice of the run
         ! (> 0.1%); tiny leftovers just clutter the tree
         if (total_sec > 0.0_wp) then
            show_other = 100.0_wp*rem/total_sec > 0.1_wp
         else
            show_other = rem > 1.0e-6_wp
         end if
      end if

      ! guide string inherited by this node's children (root's children hang off
      ! an empty prefix so they read as top-level rows)
      if (depth == 0) then
         child_prefix = ""
      else if (is_last) then
         child_prefix = prefix//tg_blank
      else
         child_prefix = prefix//tg_pipe
      end if

      prev_expanded = .false.
      do k = 1, nkids
         ! optionally space out the top-level groups (the root's direct children),
         ! but only after a group that actually expanded into printed subcategories,
         ! so runs of bare one-line groups stay compact; the blank still carries the
         ! vertical guide so the tree stays connected
         if (self%group_blanks .and. depth == 0 .and. k > 1 .and. prev_expanded) then
            write (iunit, "(1x,a)") child_prefix//tg_vert
         end if
         ! the visual last child is "(other)" when present, else the last kid
         kid_last = (k == nkids) .and. .not. show_other
         call print_node(self, iunit, kids(k), child_prefix, depth + 1, kid_last, &
            & total_sec, max_depth)
         prev_expanded = node_has_printed_children(self, kids(k), depth + 1, max_depth)
      end do

      ! remainder not attributed to any child (always the last visual sibling)
      if (show_other) then
         if (self%group_blanks .and. depth == 0 .and. prev_expanded) then
            write (iunit, "(1x,a)") child_prefix//tg_vert
         end if
         pct = 0.0_wp
         if (total_sec > 0.0_wp) pct = 100.0_wp*rem/total_sec
         write (numbuf, "(f13.4,f9.2)") rem, pct
         write (iunit, "(1x,a)") &
            & pad_disp(child_prefix//tg_last//"(other)", sect_w)//numbuf(1:22)
      end if

   end subroutine print_node

!> Accumulated seconds of a node, or NaN if it was left unbalanced.
!> @param[in] self  Timer instance
!> @param[in] id    Node index
!> @return    Seconds (NaN if poisoned/open)
   pure function node_seconds(self, id) result(sec)
      !> Timer instance
      type(timer_type), intent(in) :: self
      !> Node index
      integer, intent(in) :: id
      !> Seconds
      real(wp) :: sec

      if (id < 1 .or. id > self%n) then
         sec = ieee_value(1.0_wp, ieee_quiet_nan)
      else if (self%poisoned(id) .or. self%running(id)) then
         sec = ieee_value(1.0_wp, ieee_quiet_nan)
      else
         sec = real(self%t_acc(id), wp)/real(self%rate, wp)
      end if

   end function node_seconds

!> Whether `print_node` would emit any child rows beneath node `id`: it has at
!> least one printable child (called, or poisoned) that is not folded away by the
!> depth cap. Used to decide group spacing so the separator blank is drawn only
!> after a category that actually expanded into printed subcategories.
!> @param[in] self       Timer instance
!> @param[in] id         Node index
!> @param[in] depth      Nesting depth at which `id` is printed
!> @param[in] max_depth  Deepest level to print (children below it are folded)
!> @return    .true. if any child row would be printed under `id`
   pure function node_has_printed_children(self, id, depth, max_depth) result(has)
      !> Timer instance
      type(timer_type), intent(in) :: self
      !> Node index
      integer, intent(in) :: id
      !> Depth at which the node is printed
      integer, intent(in) :: depth
      !> Deepest level to print
      integer, intent(in) :: max_depth
      !> Whether any child row would be printed
      logical :: has

      integer :: j

      has = .false.
      ! at the depth cap the node's descendants are folded into it (none printed)
      if (depth >= max_depth) return
      if (id < 1 .or. id > self%n) return
      do j = 1, self%n
         if (self%parent(j) /= id) cycle
         if (self%ncalls(j) == 0 .and. .not. self%poisoned(j)) cycle
         has = .true.
         return
      end do

   end function node_has_printed_children

!> Resolve a slash-separated path to a node index, walking from the root.
!> @param[in] self  Timer instance
!> @param[in] path  Name or "Parent/Child/..." path
!> @return    Node index, or 0 if not found
   pure function find_path(self, path) result(id)
      !> Timer instance
      type(timer_type), intent(in) :: self
      !> Path string
      character(len=*), intent(in) :: path
      !> Node index
      integer :: id

      integer :: cur, lo, hi, j
      character(len=name_len) :: seg

      cur = 0
      id = 0
      lo = 1
      do
         if (lo > len_trim(path)) exit
         hi = index(path(lo:), "/")
         if (hi == 0) then
            seg = path(lo:len_trim(path))
            lo = len_trim(path) + 1
         else
            seg = path(lo:lo + hi - 2)
            lo = lo + hi
         end if
         ! find child named seg under cur
         id = 0
         do j = 1, self%n
            if (self%parent(j) == cur .and. self%names(j) == seg) then
               id = j
               exit
            end if
         end do
         if (id == 0) return
         cur = id
      end do

   end function find_path

!> Left-justified, blank-padded copy of a name for fixed-length comparison.
!> @param[in] name  Timer name
!> @return    Name padded/truncated to name_len
   pure function pad_name(name) result(padded)
      !> Timer name
      character(len=*), intent(in) :: name
      !> Padded name
      character(len=name_len) :: padded

      padded = name

   end function pad_name

!> FNV-1a hash of (name, parent), masked to 32 bits.
!> @param[in] name    Timer name
!> @param[in] parent  Parent node index
!> @return    Hash value in [0, 2^32)
   pure function hash_key(name, parent) result(h)
      !> Timer name
      character(len=*), intent(in) :: name
      !> Parent index
      integer, intent(in) :: parent
      !> Hash value
      integer(int64) :: h

      integer :: k

      h = fnv_offset
      do k = 1, len_trim(name)
         h = ieor(h, int(iachar(name(k:k)), int64))
         h = iand(h*fnv_prime, mask32)
      end do
      h = ieor(h, iand(int(parent, int64), mask32))
      h = iand(h*fnv_prime, mask32)

   end function hash_key

!> Insert an existing node into the hash table by linear probing.
!> @param[inout] self  Timer instance
!> @param[in]    id    Node index to insert
   subroutine table_insert(self, id)
      !> Timer instance
      class(timer_type), intent(inout) :: self
      !> Node index
      integer, intent(in) :: id

      integer(int64) :: h
      integer :: slot

      h = hash_key(self%names(id), self%parent(id))
      slot = int(iand(h, int(size(self%htab) - 1, int64)), kind(slot)) + 1
      do
         if (self%htab(slot) == 0) then
            self%htab(slot) = id
            return
         end if
         slot = slot + 1
         if (slot > size(self%htab)) slot = 1
      end do

   end subroutine table_insert

!> Double the hash table and reinsert every registered node.
!> @param[inout] self  Timer instance
   subroutine grow_table(self)
      !> Timer instance
      class(timer_type), intent(inout) :: self

      integer, allocatable :: newtab(:)
      integer :: i

      allocate (newtab(size(self%htab)*2), source=0)
      call move_alloc(newtab, self%htab)
      do i = 1, self%n
         call table_insert(self, i)
      end do

   end subroutine grow_table

!> Ensure room for one more node, doubling the node arrays if full.
!> @param[inout] self  Timer instance
   subroutine ensure_node_capacity(self)
      !> Timer instance
      class(timer_type), intent(inout) :: self

      character(len=name_len), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      integer(int64), allocatable :: lbuf(:)
      logical, allocatable :: bbuf(:)
      integer :: newcap, m

      if (self%n < self%cap) return

      newcap = self%cap*2
      m = self%cap

      allocate (character(len=name_len) :: cbuf(newcap))
      cbuf(1:m) = self%names(1:m)
      call move_alloc(cbuf, self%names)

      allocate (ibuf(newcap), source=0)
      ibuf(1:m) = self%parent(1:m)
      call move_alloc(ibuf, self%parent)

      allocate (ibuf(newcap), source=cat_none)
      ibuf(1:m) = self%category(1:m)
      call move_alloc(ibuf, self%category)

      allocate (ibuf(newcap), source=0)
      ibuf(1:m) = self%ncalls(1:m)
      call move_alloc(ibuf, self%ncalls)

      allocate (lbuf(newcap), source=0_int64)
      lbuf(1:m) = self%t_acc(1:m)
      call move_alloc(lbuf, self%t_acc)

      allocate (lbuf(newcap), source=0_int64)
      lbuf(1:m) = self%t0(1:m)
      call move_alloc(lbuf, self%t0)

      allocate (bbuf(newcap), source=.false.)
      bbuf(1:m) = self%running(1:m)
      call move_alloc(bbuf, self%running)

      allocate (bbuf(newcap), source=.false.)
      bbuf(1:m) = self%poisoned(1:m)
      call move_alloc(bbuf, self%poisoned)

      self%cap = newcap

   end subroutine ensure_node_capacity

!> Ensure room for one more open frame on the nesting stack.
!> @param[inout] self  Timer instance
   subroutine ensure_stack_capacity(self)
      !> Timer instance
      class(timer_type), intent(inout) :: self

      integer, allocatable :: sbuf(:)
      integer :: m

      if (self%sp < size(self%stack)) return

      m = size(self%stack)
      allocate (sbuf(m*2), source=0)
      sbuf(1:m) = self%stack(1:m)
      call move_alloc(sbuf, self%stack)

   end subroutine ensure_stack_capacity

!> Display width (terminal columns) of a possibly-UTF-8 string. Counts every byte
!> except UTF-8 continuation bytes (10xxxxxx); since the only multibyte glyphs the
!> report emits are the single-column box-drawing guides, this equals the number
!> of printed columns and lets fixed-width alignment survive the tree guides.
!> @param[in] s  Byte string (ASCII and/or box-drawing UTF-8)
!> @return    Number of display columns
   pure function disp_len(s) result(n)
      !> Input string
      character(len=*), intent(in) :: s
      !> Column count
      integer :: n

      integer :: k, b

      n = 0
      do k = 1, len(s)
         b = iachar(s(k:k))
         if (iand(b, 192) /= 128) n = n + 1
      end do

   end function disp_len

!> Right-pad `s` with blanks to a display width of `w` columns (no truncation if
!> it is already wider). Uses `disp_len` so a UTF-8 tree prefix pads correctly.
!> @param[in] s  String to pad
!> @param[in] w  Target display width in columns
!> @return    Padded string
   pure function pad_disp(s, w) result(padded)
      !> Input string
      character(len=*), intent(in) :: s
      !> Target width
      integer, intent(in) :: w
      !> Padded result
      character(len=:), allocatable :: padded

      integer :: dl

      dl = disp_len(s)
      if (dl < w) then
         padded = s//repeat(" ", w - dl)
      else
         padded = s
      end if

   end function pad_disp

!> The centersd "T I M I N G S" banner, `tbl_w` columns wide.
!> @return  Banner line
   pure function banner_line() result(s)
      !> Banner line
      character(len=:), allocatable :: s

      character(len=*), parameter :: title = "  T I M I N G S  "
      integer :: left

      left = (tbl_w - len(title))/2
      s = repeat("=", left)//title//repeat("=", tbl_w - left - len(title))

   end function banner_line

!> The column-header row ("section  time [s]  %tot  calls"), aligned to the same
!> fixed columns the data rows use.
!> @return  Header line
   function header_line() result(s)
      !> Header line
      character(len=:), allocatable :: s

      character(len=32) :: tail

      write (tail, "(a13,a9,a10)") "time [s]", "%tot", "calls"
      s = pad_disp("section", sect_w)//tail

   end function header_line

end module moist_utils_timer
