!> Test suite for the timer
module test_utils_timer
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use moist_utils_timer, only: timer_type, cat_setup, cat_gradient
   implicit none(type, external)
   private

   public :: collect_utils_timer

contains

   !> Collect all timer tests
   subroutine collect_utils_timer(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("loop_accumulation", test_loop_accumulation), &
         new_unittest("nesting_and_paths", test_nesting_and_paths), &
         new_unittest("unbalanced_poisons_nan", test_unbalanced_poisons_nan), &
         new_unittest("mismatch_unwind", test_mismatch_unwind), &
         new_unittest("handle_fast_path", test_handle_fast_path), &
         new_unittest("current_open_node", test_current_open_node), &
         new_unittest("introspection", test_introspection), &
         new_unittest("double_start_poisons", test_double_start_poisons), &
         new_unittest("stale_handle_safe", test_stale_handle_safe), &
         new_unittest("unwind_to_depth", test_unwind_to_depth), &
         new_unittest("reset_and_inactive", test_reset_and_inactive) &
         ]

   end subroutine collect_utils_timer

   !> Spin the CPU for a short but measurable interval
   subroutine busy()
      integer :: k
      real(wp) :: x
      x = 0.0_wp
      do k = 1, 50000
         x = x + sqrt(real(k, wp))
      end do
      if (x < 0.0_wp) error stop
   end subroutine busy

   !> Re-entering a name accumulates into one node and counts every call
   subroutine test_loop_accumulation(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t
      real(wp) :: t_once, t_first, t_five
      integer :: i

      call t%new()

      call t%start("once")
      call busy()
      call t%stop()
      t_once = t%get("once")

      call t%start("many")
      call t%start("inner")
      call busy()
      call t%stop()
      t_first = t%get("many/inner")
      do i = 2, 5
         call t%start("inner")
         call busy()
         call t%stop()
      end do
      call t%stop()
      t_five = t%get("many/inner")

      call check(error, t_once > 0.0_wp, "single interval must accrue time")
      if (allocated(error)) return
      call check(error, .not. ieee_is_nan(t_five), "looped timer must not be NaN")
      if (allocated(error)) return
      ! re-entering the name must add to the node, not replace its time
      call check(error, t_five > t_first, "five accumulations must exceed one")
      if (allocated(error)) return
      ! the parent must contain at least its child's time
      call check(error, t%get("many") >= t_five, "parent must contain child time")

      call t%delete()
   end subroutine test_loop_accumulation

   !> Identical leaf names under different parents are distinct nodes
   subroutine test_nesting_and_paths(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t

      call t%new()

      ! same leaf name "solve" under two different parents
      call t%start("A")
      call t%start("solve")
      call busy()
      call t%stop()
      call t%stop()

      call t%start("B")
      call t%start("solve")
      call t%start("deep")
      call busy()
      call t%stop()
      call t%stop()
      call t%stop()

      call check(error, t%get("A/solve") > 0.0_wp, "A/solve must exist")
      if (allocated(error)) return
      call check(error, t%get("B/solve/deep") > 0.0_wp, "3-level path must resolve")
      if (allocated(error)) return
      ! an unknown path returns exactly zero (never created)
      call check(error, t%get("A/nonexistent") == 0.0_wp, "unknown path -> 0")

      call t%delete()
   end subroutine test_nesting_and_paths

   !> A timer started and never stopped reports NaN
   subroutine test_unbalanced_poisons_nan(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t

      call t%new()

      call t%start("closed")
      call busy()
      call t%stop()

      call t%start("leaked")
      ! deliberately no stop

      call check(error, .not. ieee_is_nan(t%get("closed")), "balanced timer is finite")
      if (allocated(error)) return
      call check(error, ieee_is_nan(t%get("leaked")), "unstopped timer must be NaN")

      call t%delete()
   end subroutine test_unbalanced_poisons_nan

   !> stop("outer") while "inner" is still open unwinds and poisons "inner",
   !> while "outer" itself closes cleanly
   subroutine test_mismatch_unwind(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t

      call t%new()

      call t%start("outer")
      call t%start("inner")
      call busy()
      ! close the outer frame by name while inner is still open
      call t%stop("outer")

      call check(error, ieee_is_nan(t%get("outer/inner")), &
                 "skipped inner frame must be poisoned")
      if (allocated(error)) return
      call check(error, .not. ieee_is_nan(t%get("outer")), &
                 "named-stop target must close cleanly")

      call t%delete()
   end subroutine test_mismatch_unwind

   !> The handle (resolve + start/stop id) fast path accumulates into the same
   !> node as the name-based path
   subroutine test_handle_fast_path(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t
      integer :: h_parent, h_child
      integer :: i

      call t%new()

      h_parent = t%resolve("grad", 0, cat_gradient)
      h_child = t%resolve("kernel", h_parent)

      call t%start(h_parent)
      do i = 1, 4
         call t%start(h_child)
         call busy()
         call t%stop(h_child)
      end do
      call t%stop(h_parent)

      call check(error, t%get("grad/kernel") > 0.0_wp, "handle child must accrue time")
      if (allocated(error)) return
      call check(error, .not. ieee_is_nan(t%get("grad")), "handle parent must be finite")

      call t%delete()
   end subroutine test_handle_fast_path

   !> current() returns the innermost open node so a handle-based consumer can
   !> nest its timers under whatever the caller opened by name (0 when nothing is
   !> open). This is what the DROP gradient hot loop uses to attach its fine
   !> timers under the caller's "Gradients" node regardless of nesting depth
   subroutine test_current_open_node(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t
      integer :: outer, inner, child

      call t%new()

      call check(error, t%current() == 0, "no open node initially")
      if (allocated(error)) return

      call t%start("outer")
      outer = t%resolve("outer", 0)
      call check(error, t%current() == outer, "current is the open outer node")
      if (allocated(error)) return

      call t%start("inner")
      inner = t%resolve("inner", outer)
      call check(error, t%current() == inner, "current follows into the inner node")
      if (allocated(error)) return

      ! a handle resolved against current() nests under the open node
      child = t%resolve("child", t%current())
      call check(error, t%node_depth(child) == 2, "child of current sits at depth 2")
      if (allocated(error)) return

      call t%stop()
      call check(error, t%current() == outer, "current pops back to outer")
      if (allocated(error)) return
      call t%stop()
      call check(error, t%current() == 0, "current empty after all stops")

      call t%delete()
   end subroutine test_current_open_node

   !> The read-only introspection API enumerates the registered tree in
   !> pre-order (parent before child), with names, depths and times consistent
   !> with path lookup. This is what programmatic consumers (e.g. the DROP
   !> timing benchmark) use instead of a hardcoded node table
   subroutine test_introspection(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t

      call t%new()

      call t%start("A", category=cat_setup)   ! node 1, depth 0
      call t%start("b1")                       ! node 2, depth 1
      call busy()
      call t%stop()
      call t%start("b2")                       ! node 3, depth 1
      call t%start("c")                        ! node 4, depth 2
      call busy()
      call t%stop()
      call t%stop()
      call t%stop()

      call check(error, t%num_nodes() == 4, "expected four registered nodes")
      if (allocated(error)) return
      ! registration order is pre-order: parent registered before its children
      call check(error, t%node_name(1) == "A" .and. t%node_depth(1) == 0, &
                 "node 1 is top-level A")
      if (allocated(error)) return
      call check(error, t%node_name(4) == "c" .and. t%node_depth(4) == 2, &
                 "node 4 is c at depth 2")
      if (allocated(error)) return
      ! time by index matches time by path
      call check(error, t%node_time(4) == t%get("A/b2/c"), &
                 "node_time(id) must equal get(path)")

      call t%delete()
   end subroutine test_introspection

   !> A second start on an already-running node is an imbalance (a missing stop):
   !> it must poison the node so the dropped interval surfaces as NaN rather than
   !> being silently discarded
   subroutine test_double_start_poisons(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t
      integer :: h

      call t%new()
      h = t%resolve("x", 0)

      call t%start(h)
      call busy()
      call t%start(h)   ! second start without an intervening stop -> imbalance
      call busy()
      call t%stop(h)

      call check(error, ieee_is_nan(t%node_time(h)), "double-started node must be NaN")

      call t%delete()
   end subroutine test_double_start_poisons

   !> Out-of-range / stale handles are a safe no-op on every id-taking entry
   !> point, never an out-of-bounds access
   subroutine test_stale_handle_safe(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t
      integer :: bad

      call t%new()
      call t%start("real")
      call t%stop()

      bad = t%num_nodes() + 100   ! never a registered id

      ! none of these may crash or corrupt the tree
      call t%start(bad)
      call t%stop(bad)
      call t%start(0)
      call t%stop(-3)

      call check(error, t%num_nodes() == 1, "stale ids register no nodes")
      if (allocated(error)) return
      call check(error, ieee_is_nan(t%node_time(bad)), "stale id time is NaN")
      if (allocated(error)) return
      call check(error, len(t%node_name(bad)) == 0, "stale id name is empty")
      if (allocated(error)) return
      call check(error, t%node_depth(bad) == 0, "stale id depth is 0")
      if (allocated(error)) return
      ! the real node is untouched
      call check(error, .not. ieee_is_nan(t%node_time(1)), "real node still finite")

      call t%delete()
   end subroutine test_stale_handle_safe

   !> unwind(depth) closes frames left open by an early (error) return, restoring
   !> the shared stack to a saved depth; the closed frames record partial time and
   !> are NOT poisoned (they ran, the routine just bailed out)
   subroutine test_unwind_to_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t
      integer :: d, outer

      call t%new()

      call t%start("outer")
      outer = t%resolve("outer", 0)
      d = t%current_depth()            ! baseline with only "outer" open

      call t%start("mid")
      call t%start("leaf")             ! two extra frames left open
      call check(error, t%current_depth() == d + 2, "two frames opened above baseline")
      if (allocated(error)) return

      call t%unwind(d)                 ! simulate an early return unwinding to baseline

      call check(error, t%current_depth() == d, "unwind restores the saved depth")
      if (allocated(error)) return
      call check(error, t%current() == outer, "outer is open again after unwind")
      if (allocated(error)) return
      ! unwound frames were stopped cleanly (partial time), not poisoned
      call check(error, .not. ieee_is_nan(t%get("outer/mid")), "unwound frame closed cleanly")
      if (allocated(error)) return

      call t%stop()                    ! close outer -> balanced
      call check(error, t%current_depth() == 0, "stack empty after final stop")

      call t%delete()
   end subroutine test_unwind_to_depth

   !> reset() zeroes accumulated time but keeps the registered tree; an inactive
   !> (deleted) timer is a safe no-op that reports zero
   subroutine test_reset_and_inactive(error)
      type(error_type), allocatable, intent(out) :: error
      type(timer_type) :: t

      call t%new()

      call t%start("step", category=cat_setup)
      call busy()
      call t%stop()
      call check(error, t%get("step") > 0.0_wp, "time accrued before reset")
      if (allocated(error)) return

      call t%reset()
      ! node still resolvable (tree kept) but time zeroed
      call check(error, t%get("step") == 0.0_wp, "reset zeroes accumulated time")
      if (allocated(error)) return

      ! deleted timer: all queries are safe no-ops
      call t%delete()
      call check(error, t%get() == 0.0_wp, "inactive total is zero")
      if (allocated(error)) return
      ! start/stop on an inactive timer must not crash
      call t%start("noop")
      call t%stop()
   end subroutine test_reset_and_inactive

end module test_utils_timer
