!> Certified best-first octree search for all closest-point branches
!>
!> Enumerates every local minimum of the distance to an implicit surface that
!> can lie within a given radius of an anchor, together with a proof that no
!> further one was missed. The proof is purely geometric: minima of the
!> distance live *on* the surface, so a region certified to contain no surface
!> contains no minimum. Covering the search ball with such regions therefore
!> reduces "did we find them all?" to "is the ball covered?".
!>
!> The caller supplies that certificate through a probe callback returning, at
!> a point `x`, both the level set value `S(x)` and a radius `r` such that `S`
!> has no zero in `B(x, r)`. Nothing else about the surface is assumed, so the
!> search is shared by every level set that can bound its own gradient (the
!> SvdW LSF is 1-Lipschitz and returns `r = |S(x)|` exactly).
!>
!> Algorithm -- best-first branch and bound on a cube octree:
!>
!>     root = cube(centre = anchor, half = rho_max)
!>     loop:
!>       pop the box with the smallest lower bound on ||x - anchor||
!>       if that bound exceeds rho_max  -> DONE, the rest of the ball is proven
!>       if |S(centre)| >= circumradius -> discard, certified surface-free
!>       if the box is at seed size     -> keep as a survivor
!>       else                           -> split into 8 and push
!>
!> Popping in order of the distance lower bound means the search sweeps
!> outwards from the anchor, so the boxes that can hold the *closest* branch
!> are examined first and the ones that cannot hold any are never examined at
!> all. When the heap's best box is worse than `rho_max`, every remaining box
!> is too, and the loop stops with the ball fully accounted for.
!>
!> `rho_max` also tightens itself as the search runs. Whenever a probed centre
!> `c` has the opposite sign to the anchor, the segment from the anchor to `c`
!> crosses the surface, so the closest branch is no further out than that
!> crossing. The crossing cannot lie inside `c`'s own surface-free ball either,
!> which puts it a further `r` short of the centre:
!>
!>     rho_min <= ||c - anchor|| - r
!>
!> Both halves of that bound are quantities the probe has already returned, so
!> the ball the search still has to certify shrinks for free.
!>
!> What the certificate covers is the *coverage of the search ball*: every box
!> is either examined or proven to hold no surface, so no branch escapes the
!> surviving leaves. Turning those leaves into seeds is a separate step, and in
!> `octree_seed_cluster` mode it is a heuristic -- one seed per discrete local
!> minimum of an estimated surface distance, which can merge two basins that
!> share a survivor patch. `octree_seed_per_leaf` performs no such reduction and
!> is the mode to fall back on when that matters.
!>
!> This is deliberately not a [[solver_base_type]]: it does not refine a single
!> iterate towards a solution, it partitions a region and returns seeds. The
!> caller runs its own local solver on those seeds afterwards.
module moist_math_solver_octree_branch
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use, intrinsic :: iso_fortran_env, only: int64, output_unit
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none
   private

   public :: moist_math_octree_branch_type
   public :: octree_seed_cluster, octree_seed_per_leaf

   !> Emit one seed per discrete local minimum of the survivor set (default)
   integer, parameter :: octree_seed_cluster = 1
   !> Emit one seed per surviving leaf
   integer, parameter :: octree_seed_per_leaf = 2

   !> Largest octree depth the integer lattice key can address without
   !> overflowing a 64-bit integer (the key packs three coordinates in base
   !> `2**depth`). Far above any depth a sane `seed_size` produces.
   integer, parameter :: max_addressable_depth = 20

   !> Number of radius-tightening events kept for the debug report
   integer, parameter :: max_logged_tightenings = 12

   !> Width of the trace's action column, and of the buffers written into it
   integer, parameter :: trace_action_width = 40

   abstract interface
      !> Probe the level set at a point
      !>
      !> @param[in]  x       Evaluation point (3)
      !> @param[out] lsf0    Level set value at `x`
      !> @param[out] radius  Radius of a ball around `x` free of surface; zero
      !>                     when the level set cannot certify one, which makes
      !>                     the search prune nothing and run into its budget
      !> @param[in]  context Caller context forwarded unchanged
      subroutine octree_probe_context_interface(x, lsf0, radius, context)
         import :: wp
         real(wp), intent(in) :: x(3)
         real(wp), intent(out) :: lsf0
         real(wp), intent(out) :: radius
         class(*), intent(in) :: context
      end subroutine octree_probe_context_interface
   end interface

   !> Best-first octree branch search.
   !>
   !> One instance per thread: [[octree_init]] sizes the scratch once and every
   !> [[octree_run]] reuses it, so the per-anchor hot path allocates nothing.
   type :: moist_math_octree_branch_type

      !*--------------------------- Search configuration -------------------- *!

      !> Edge length at which a surviving box stops splitting and becomes a
      !> seed candidate (Bohr). Sets the resolution at which two minima are
      !> still told apart.
      real(wp) :: seed_size = 0.2_wp
      !> Largest number of boxes examined per anchor before the search gives up
      integer :: max_boxes = 200000
      !> Largest number of surviving leaves retained per anchor
      integer :: max_survivors = 200000
      !> Largest octree depth
      integer :: max_depth = 12
      !> How survivors are turned into seeds (`octree_seed_*`)
      integer :: seed_mode = octree_seed_cluster
      !> Emit a per-anchor trace of what the search did
      logical :: debug = .false.
      !> Output unit for the trace
      integer :: unit = output_unit
      !> Anchor id shown in the trace's first column; set by the caller
      integer :: point_id = 0

      !*------------------------------- Results ------------------------------ *!

      !> Seed points handed back to the caller's local solver (3, n_seeds)
      real(wp), allocatable :: seeds(:, :)
      !> Number of seeds produced by the last run
      integer :: n_seeds = 0
      !> Radius the last run ended up certifying (after self-tightening)
      real(wp) :: rho_max_final = 0.0_wp
      !> Boxes examined by the last run, for cost diagnostics
      integer :: n_boxes_visited = 0

      !*---------------------------- Debug bookkeeping ----------------------- *!

      !> Per-depth tallies of the last run: boxes popped, boxes certified
      !> surface-free, boxes split, and boxes kept as survivors.
      integer :: dbg_popped(0:max_addressable_depth) = 0
      integer :: dbg_free(0:max_addressable_depth) = 0
      integer :: dbg_split(0:max_addressable_depth) = 0
      integer :: dbg_kept(0:max_addressable_depth) = 0
      !> Were tallies collected?
      logical :: dbg_valid = .false.
      !> Radius-tightening events: the box number, the distance to the probed
      !> centre that produced the new upper bound, that distance less the
      !> centre's own surface-free radius, and the resulting cap
      integer :: dbg_tighten_box(max_logged_tightenings) = 0
      real(wp) :: dbg_tighten_raw(max_logged_tightenings) = 0.0_wp
      real(wp) :: dbg_tighten_hit(max_logged_tightenings) = 0.0_wp
      real(wp) :: dbg_tighten_cap(max_logged_tightenings) = 0.0_wp
      integer :: n_tightenings = 0
      !> Run inputs kept for the report
      real(wp) :: dbg_lsf0_anchor = 0.0_wp
      real(wp) :: dbg_rho_start = 0.0_wp
      real(wp) :: dbg_rho2_slack = 0.0_wp
      integer :: dbg_leaf_depth = 0
      !> Survivors before the final cap trimmed them
      integer :: dbg_n_surv_raw = 0
      !> Peak heap occupancy
      integer :: dbg_heap_peak = 0
      !> Table used for the live walk trace
      type(prettylistprinter) :: trace_tbl

      !*----------------------------- Working state -------------------------- *!

      !> Binary min-heap of pending boxes, keyed on the box's lower bound for
      !> ||x - anchor||. A box is stored as its octree depth plus its integer
      !> lattice position at that depth, from which centre and half-width
      !> follow exactly -- no accumulated floating-point drift down the tree.
      real(wp), allocatable :: heap_key(:)
      integer, allocatable :: heap_depth(:)
      integer, allocatable :: heap_lat(:, :)
      integer :: heap_size = 0

      !> Surviving leaves: lattice position, centre, the box's distance lower
      !> bound (used to apply the final cap) and an estimate of how far the
      !> *surface inside that box* is from the anchor (used to rank them).
      !>
      !> Neither of the two obvious rankings works. The lower bound is
      !> quantized, so whole rings of leaves share one value and a smooth patch
      !> becomes a plateau on which dozens of cells look like minima. The centre
      !> distance is worse: survivors form a slab straddling the surface, and a
      !> centre sitting inside the slab is nearer the anchor than the surface it
      !> stands for, which drags the ranking onto the slab's inner face and
      !> spreads the minimum around a ring. The estimate below corrects the
      !> centre distance by the level set value -- the surface lies |S| away
      !> from the centre, beyond it when the centre is on the anchor's side of
      !> the surface and nearer when it is not -- which is exact for a signed
      !> distance field with a locally flat surface and puts the minimum where
      !> the surface actually comes closest.
      integer, allocatable :: surv_lat(:, :)
      real(wp), allocatable :: surv_centre(:, :)
      real(wp), allocatable :: surv_rho(:)
      real(wp), allocatable :: surv_rho_est(:)
      integer :: n_surv = 0

      !> Survivor lattice keys, sorted, with the permutation that produced
      !> them. Used to answer face-neighbour queries by binary search.
      integer(int64), allocatable :: surv_key(:)
      integer(int64), allocatable :: surv_key_sorted(:)
      integer, allocatable :: surv_order(:)

      !> Hard ceiling on the heap size, derived from the survivor budget
      integer :: heap_capacity_max = 0

      !> Anchor and root half-width of the run in progress
      real(wp) :: anchor(3) = 0.0_wp
      real(wp) :: root_half = 0.0_wp
      !> Lattice edge count at the leaf depth, the base of the lattice key
      integer(int64) :: lattice_base = 1_int64

   contains
      !> Size the scratch arrays; call once per thread
      procedure :: init => octree_init
      !> Run the search for one anchor
      procedure :: run => octree_run
      !> Print a summary of the last run
      procedure :: report => octree_report
      !> Release the scratch arrays
      procedure :: destroy => octree_destroy
   end type moist_math_octree_branch_type

contains

!* ================================================================================= *!
!*                                     Lifecycle                                     *!
!* ================================================================================= *!

   !> Configure the search and allocate a starting scratch
   !>
   !> The arrays grow on demand up to the configured budgets, so an anchor
   !> whose ball is nearly empty never pays for the worst case, while the
   !> budgets stay a hard ceiling on a pathological one.
   !>
   !> @param[inout] self          Search instance
   !> @param[in]    seed_size     Edge length at which boxes stop splitting (Bohr)
   !> @param[in]    max_boxes     Box budget per anchor
   !> @param[in]    max_survivors Survivor budget per anchor
   !> @param[in]    max_depth     Depth limit
   !> @param[in]    seed_mode     One of the `octree_seed_*` constants
   !> @param[out]   error         Invalid configuration
   subroutine octree_init(self, seed_size, max_boxes, max_survivors, max_depth, &
                          seed_mode, error)
      class(moist_math_octree_branch_type), intent(inout) :: self
      real(wp), intent(in), optional :: seed_size
      integer, intent(in), optional :: max_boxes
      integer, intent(in), optional :: max_survivors
      integer, intent(in), optional :: max_depth
      integer, intent(in), optional :: seed_mode
      type(error_type), allocatable, intent(out) :: error

      integer :: heap_capacity

      !> Starting size of the growable scratch arrays
      integer, parameter :: initial_capacity = 1024

      if (present(seed_size)) self%seed_size = seed_size
      if (present(max_boxes)) self%max_boxes = max_boxes
      if (present(max_survivors)) self%max_survivors = max_survivors
      if (present(max_depth)) self%max_depth = max_depth
      if (present(seed_mode)) self%seed_mode = seed_mode

      if (self%seed_size <= 0.0_wp) then
         call fatal_error(error, "Octree branch search: seed_size must be positive")
         return
      end if
      if (self%max_boxes < 1 .or. self%max_survivors < 1) then
         call fatal_error(error, "Octree branch search: box and survivor budgets must be positive")
         return
      end if
      if (self%max_depth < 1 .or. self%max_depth > max_addressable_depth) then
         call fatal_error(error, "Octree branch search: max_depth out of range")
         return
      end if
      if (self%seed_mode /= octree_seed_cluster .and. &
          self%seed_mode /= octree_seed_per_leaf) then
         call fatal_error(error, "Octree branch search: unknown seed mode")
         return
      end if

      call self%destroy()

      ! Splitting a box removes one entry and adds eight, so the heap never
      ! needs more than eight slots per surviving leaf plus a small margin.
      self%heap_capacity_max = 8*self%max_survivors + 64

      heap_capacity = min(initial_capacity, self%heap_capacity_max)
      allocate (self%heap_key(heap_capacity))
      allocate (self%heap_depth(heap_capacity))
      allocate (self%heap_lat(3, heap_capacity))

      call grow_survivor_arrays(self, min(initial_capacity, self%max_survivors))
   end subroutine octree_init

   !> Grow the heap arrays to at least `required`, capped by the budget
   !> @param[inout] self     Search instance
   !> @param[in]    required Capacity the heap must reach
   !> @param[out]   ok       `.false.` when the budget forbids the growth
   subroutine grow_heap_arrays(self, required, ok)
      class(moist_math_octree_branch_type), intent(inout) :: self
      integer, intent(in) :: required
      logical, intent(out) :: ok

      real(wp), allocatable :: new_key(:)
      integer, allocatable :: new_depth(:), new_lat(:, :)
      integer :: new_capacity, old_capacity

      old_capacity = size(self%heap_key)
      ok = required <= self%heap_capacity_max
      if (.not. ok) return
      if (required <= old_capacity) return

      new_capacity = min(max(2*old_capacity, required), self%heap_capacity_max)

      allocate (new_key(new_capacity))
      allocate (new_depth(new_capacity))
      allocate (new_lat(3, new_capacity))
      new_key(1:old_capacity) = self%heap_key(1:old_capacity)
      new_depth(1:old_capacity) = self%heap_depth(1:old_capacity)
      new_lat(:, 1:old_capacity) = self%heap_lat(:, 1:old_capacity)

      call move_alloc(new_key, self%heap_key)
      call move_alloc(new_depth, self%heap_depth)
      call move_alloc(new_lat, self%heap_lat)
   end subroutine grow_heap_arrays

   !> Allocate or grow the survivor arrays to `capacity`, preserving contents
   !> @param[inout] self     Search instance
   !> @param[in]    capacity Number of survivors the arrays must hold
   subroutine grow_survivor_arrays(self, capacity)
      class(moist_math_octree_branch_type), intent(inout) :: self
      integer, intent(in) :: capacity

      integer, allocatable :: new_lat(:, :), new_order(:)
      real(wp), allocatable :: new_centre(:, :), new_rho(:), new_rho_est(:)
      integer(int64), allocatable :: new_key(:), new_key_sorted(:)
      integer :: kept

      kept = 0
      if (allocated(self%surv_rho)) kept = min(self%n_surv, size(self%surv_rho))

      allocate (new_lat(3, capacity))
      allocate (new_centre(3, capacity))
      allocate (new_rho(capacity))
      allocate (new_rho_est(capacity))
      allocate (new_key(capacity))
      allocate (new_key_sorted(capacity))
      allocate (new_order(capacity))

      if (kept > 0) then
         new_lat(:, 1:kept) = self%surv_lat(:, 1:kept)
         new_centre(:, 1:kept) = self%surv_centre(:, 1:kept)
         new_rho(1:kept) = self%surv_rho(1:kept)
         new_rho_est(1:kept) = self%surv_rho_est(1:kept)
      end if

      call move_alloc(new_lat, self%surv_lat)
      call move_alloc(new_centre, self%surv_centre)
      call move_alloc(new_rho, self%surv_rho)
      call move_alloc(new_rho_est, self%surv_rho_est)
      call move_alloc(new_key, self%surv_key)
      call move_alloc(new_key_sorted, self%surv_key_sorted)
      call move_alloc(new_order, self%surv_order)

      if (allocated(self%seeds)) deallocate (self%seeds)
      allocate (self%seeds(3, capacity))
   end subroutine grow_survivor_arrays

   !> Release the scratch arrays
   !> @param[inout] self Search instance
   subroutine octree_destroy(self)
      class(moist_math_octree_branch_type), intent(inout) :: self

      if (allocated(self%heap_key)) deallocate (self%heap_key)
      if (allocated(self%heap_depth)) deallocate (self%heap_depth)
      if (allocated(self%heap_lat)) deallocate (self%heap_lat)
      if (allocated(self%surv_lat)) deallocate (self%surv_lat)
      if (allocated(self%surv_centre)) deallocate (self%surv_centre)
      if (allocated(self%surv_rho)) deallocate (self%surv_rho)
      if (allocated(self%surv_rho_est)) deallocate (self%surv_rho_est)
      if (allocated(self%surv_key)) deallocate (self%surv_key)
      if (allocated(self%surv_key_sorted)) deallocate (self%surv_key_sorted)
      if (allocated(self%surv_order)) deallocate (self%surv_order)
      if (allocated(self%seeds)) deallocate (self%seeds)

      self%heap_size = 0
      self%n_surv = 0
      self%n_seeds = 0
   end subroutine octree_destroy

!* ================================================================================= *!
!*                                    Box geometry                                   *!
!* ================================================================================= *!

   !> Half-width of a box at a given depth
   !> @param[in] self  Search instance (carries the root half-width)
   !> @param[in] depth Octree depth (0 = root)
   !> @returns   half  Half edge length at that depth
   pure function box_half(self, depth) result(half)
      class(moist_math_octree_branch_type), intent(in) :: self
      integer, intent(in) :: depth
      real(wp) :: half

      half = self%root_half/real(2**depth, wp)
   end function box_half

   !> Depth at which a root of half-width `rho_max` first reaches `seed_size`
   !>
   !> Only used to tell the caller what to set `max_depth` to; the search itself
   !> walks the depth down rather than computing it in floating point.
   !>
   !> @param[in] self    Search instance
   !> @param[in] rho_max Root half-width (Bohr)
   !> @returns   depth   Smallest depth with edge length at or below `seed_size`
   pure function required_depth(self, rho_max) result(depth)
      class(moist_math_octree_branch_type), intent(in) :: self
      real(wp), intent(in) :: rho_max
      integer :: depth

      depth = 0
      do while (2.0_wp*rho_max/real(2**depth, wp) > self%seed_size)
         depth = depth + 1
         if (depth >= max_addressable_depth) exit
      end do
   end function required_depth

   !> Offset of a box centre from the anchor
   !>
   !> The root is centred on the anchor, so at depth `d` with half-width `h`
   !> the lattice cell `lat` sits at `h*(2*lat + 1 - 2**d)` -- exact in integer
   !> arithmetic, unlike accumulating half-steps while descending the tree.
   !>
   !> @param[in] self  Search instance
   !> @param[in] depth Octree depth
   !> @param[in] lat   Lattice position at that depth
   !> @returns   d     Centre offset from the anchor (3)
   pure function box_offset(self, depth, lat) result(d)
      class(moist_math_octree_branch_type), intent(in) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      real(wp) :: d(3)

      d = box_half(self, depth)*real(2*lat + 1 - 2**depth, wp)
   end function box_offset

   !> Smallest distance from the anchor to any point of a box
   !> @param[in] self  Search instance
   !> @param[in] depth Octree depth
   !> @param[in] lat   Lattice position at that depth
   !> @returns   rho   Lower bound on ||x - anchor|| over the box
   pure function box_rho_lower(self, depth, lat) result(rho)
      class(moist_math_octree_branch_type), intent(in) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      real(wp) :: rho

      real(wp) :: d(3), half, gap(3)

      half = box_half(self, depth)
      d = box_offset(self, depth, lat)
      gap = max(abs(d) - half, 0.0_wp)
      rho = norm2(gap)
   end function box_rho_lower

!* ================================================================================= *!
!*                                   Min-heap                                        *!
!* ================================================================================= *!

   !> Push a box onto the heap, sifting it up into place
   !> @param[inout] self  Search instance
   !> @param[in]    key   Distance lower bound of the box
   !> @param[in]    depth Octree depth
   !> @param[in]    lat   Lattice position
   !> @param[out]   ok    `.false.` when the heap is full
   subroutine heap_push(self, key, depth, lat, ok)
      class(moist_math_octree_branch_type), intent(inout) :: self
      real(wp), intent(in) :: key
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      logical, intent(out) :: ok

      integer :: child, parent

      if (self%heap_size >= size(self%heap_key)) then
         call grow_heap_arrays(self, self%heap_size + 1, ok)
         if (.not. ok) return
      end if
      ok = .true.

      self%heap_size = self%heap_size + 1
      child = self%heap_size
      self%heap_key(child) = key
      self%heap_depth(child) = depth
      self%heap_lat(:, child) = lat

      do while (child > 1)
         parent = child/2
         if (self%heap_key(parent) <= self%heap_key(child)) exit
         call heap_swap(self, parent, child)
         child = parent
      end do
   end subroutine heap_push

   !> Pop the box with the smallest distance lower bound
   !> @param[inout] self  Search instance
   !> @param[out]   key   Distance lower bound of the popped box
   !> @param[out]   depth Octree depth
   !> @param[out]   lat   Lattice position
   subroutine heap_pop(self, key, depth, lat)
      class(moist_math_octree_branch_type), intent(inout) :: self
      real(wp), intent(out) :: key
      integer, intent(out) :: depth
      integer, intent(out) :: lat(3)

      integer :: parent, child

      key = self%heap_key(1)
      depth = self%heap_depth(1)
      lat = self%heap_lat(:, 1)

      self%heap_key(1) = self%heap_key(self%heap_size)
      self%heap_depth(1) = self%heap_depth(self%heap_size)
      self%heap_lat(:, 1) = self%heap_lat(:, self%heap_size)
      self%heap_size = self%heap_size - 1

      parent = 1
      do
         child = 2*parent
         if (child > self%heap_size) exit
         if (child < self%heap_size) then
            if (self%heap_key(child + 1) < self%heap_key(child)) child = child + 1
         end if
         if (self%heap_key(parent) <= self%heap_key(child)) exit
         call heap_swap(self, parent, child)
         parent = child
      end do
   end subroutine heap_pop

   !> Exchange two heap entries
   !> @param[inout] self Search instance
   !> @param[in]    i    First heap position
   !> @param[in]    j    Second heap position
   pure subroutine heap_swap(self, i, j)
      class(moist_math_octree_branch_type), intent(inout) :: self
      integer, intent(in) :: i, j

      real(wp) :: tmp_key
      integer :: tmp_depth, tmp_lat(3)

      tmp_key = self%heap_key(i)
      tmp_depth = self%heap_depth(i)
      tmp_lat = self%heap_lat(:, i)

      self%heap_key(i) = self%heap_key(j)
      self%heap_depth(i) = self%heap_depth(j)
      self%heap_lat(:, i) = self%heap_lat(:, j)

      self%heap_key(j) = tmp_key
      self%heap_depth(j) = tmp_depth
      self%heap_lat(:, j) = tmp_lat
   end subroutine heap_swap

!* ================================================================================= *!
!*                                    Main search                                    *!
!* ================================================================================= *!

   !> Enumerate all branch seeds within the admissible radius of one anchor
   !>
   !> On return `self%seeds(:, 1:self%n_seeds)` holds the seed points and
   !> `self%rho_max_final` the radius actually certified. An empty seed list is
   !> a legitimate answer -- it means the whole ball was proven surface-free.
   !>
   !> @param[inout] self         Search instance (must be initialized)
   !> @param[in]    anchor       Anchor point (3)
   !> @param[in]    lsf0_anchor  Level set value at the anchor, for the sign test
   !> @param[in]    rho_max      Starting admissible radius (Bohr)
   !> @param[in]    rho2_slack   Squared-distance slack defining the admissible
   !>                            set: rho_max^2 = rho_min^2 + rho2_slack. Used to
   !>                            re-derive rho_max whenever the search improves
   !>                            its upper bound on rho_min
   !> @param[in]    probe        Level set probe callback
   !> @param[in]    context      Context forwarded to `probe`
   !> @param[out]   error        Budget exhausted, i.e. no certificate obtained
   subroutine octree_run(self, anchor, lsf0_anchor, rho_max, rho2_slack, &
                         probe, context, error)
      class(moist_math_octree_branch_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lsf0_anchor
      real(wp), intent(in) :: rho_max
      real(wp), intent(in) :: rho2_slack
      procedure(octree_probe_context_interface) :: probe
      class(*), intent(in) :: context
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: rho_cap, rho_hit, key, half, centre(3), lsf0, excl_radius
      real(wp) :: child_key, excl_ratio, rho_centre, rho_cross
      integer :: depth, lat(3), child_lat(3), leaf_depth
      integer :: ix, iy, iz, n_pushed
      logical :: ok, crosses, trace
      !> Trace text for the box currently being examined. Sized to the trace's
      !> action column; an internal write past it is a runtime "End of record".
      character(len=trace_action_width) :: action
      !> Scratch for composing error messages
      character(len=256) :: buffer

      if (.not. allocated(self%heap_key)) then
         call fatal_error(error, "Octree branch search: run before init")
         return
      end if
      if (rho_max <= 0.0_wp) then
         call fatal_error(error, "Octree branch search: rho_max must be positive")
         return
      end if
      ! A negative slack would put a NaN into the radius cap, after which every
      ! comparison against it is false and the search silently stops bounding
      ! anything. Callers derive it from a weight floor, so a bad floor lands here.
      if (rho2_slack < 0.0_wp) then
         call fatal_error(error, "Octree branch search: rho2_slack must not be negative")
         return
      end if

      self%anchor = anchor
      self%root_half = rho_max
      self%heap_size = 0
      self%n_surv = 0
      self%n_seeds = 0
      self%n_boxes_visited = 0

      trace = self%debug
      self%dbg_valid = trace
      self%dbg_popped = 0
      self%dbg_free = 0
      self%dbg_split = 0
      self%dbg_kept = 0
      self%n_tightenings = 0
      self%dbg_heap_peak = 0
      self%dbg_n_surv_raw = 0
      self%dbg_lsf0_anchor = lsf0_anchor
      self%dbg_rho_start = rho_max
      self%dbg_rho2_slack = rho2_slack

      ! Depth at which the edge length first drops to the seed size. Every
      ! survivor stops there, so they all share one lattice -- which is what
      ! lets face adjacency be an integer test further down.
      !
      ! Stopping early because `max_depth` ran out would hand back leaves
      ! coarser than asked for, merging minima that the configured resolution
      ! separates -- and still report a completed search. The other caps are
      ! errors for the same reason.
      leaf_depth = 0
      do while (2.0_wp*box_half(self, leaf_depth) > self%seed_size)
         leaf_depth = leaf_depth + 1
         if (leaf_depth >= self%max_depth) exit
      end do
      if (2.0_wp*box_half(self, leaf_depth) > self%seed_size) then
         write (buffer, "(a,i0,a,f0.4,a,f0.4,a,i0,a)") &
            "Octree branch search: max_depth (", self%max_depth, &
            ") bottoms out at a leaf edge of ", 2.0_wp*box_half(self, leaf_depth), &
            " Bohr, coarser than the requested seed_size of ", self%seed_size, &
            " Bohr. Raise max_depth to at least ", required_depth(self, rho_max), &
            " or loosen seed_size."
         call fatal_error(error, trim(buffer))
         return
      end if
      self%lattice_base = 2_int64**leaf_depth
      self%dbg_leaf_depth = leaf_depth

      excl_ratio = 0.0_wp
      rho_cap = rho_max
      ! An anchor sitting exactly on the surface is its own closest branch.
      rho_hit = huge(1.0_wp)
      if (lsf0_anchor == 0.0_wp) then
         rho_hit = 0.0_wp
         rho_cap = min(rho_cap, sqrt(rho2_slack))
      end if

      if (self%debug) then
         call trace_banner(self, anchor, leaf_depth)
         call trace_open(self, leaf_depth)
      end if

      call heap_push(self, 0.0_wp, 0, [0, 0, 0], ok)
      if (.not. ok) then
         call fatal_error(error, "Octree branch search: heap capacity too small for the root")
         return
      end if

      do while (self%heap_size > 0)
         call heap_pop(self, key, depth, lat)

         ! Best-first: every box still queued is at least this far out, so a
         ! key beyond the cap proves the remainder of the ball carries no
         ! branch weight and the search is finished.
         if (key > rho_cap) then
            if (trace) then
               write (action, "(a)") "beyond rho_max: search complete"
               call trace_row(self, depth, lat, leaf_depth, 0.0_wp, 0.0_wp, key, action)
            end if
            exit
         end if

         self%n_boxes_visited = self%n_boxes_visited + 1
         if (trace) self%dbg_popped(depth) = self%dbg_popped(depth) + 1
         if (self%n_boxes_visited > self%max_boxes) then
            call fatal_error(error, &
                             "Octree branch search: box budget exhausted before the "// &
                             "search ball was certified. Either the level set supplies "// &
                             "no exclusion radius (pruning nothing), or seed_size is too "// &
                             "small for the admissible radius.")
            return
         end if

         half = box_half(self, depth)
         centre = anchor + box_offset(self, depth, lat)

         call probe(centre, lsf0, excl_radius, context)

         if (trace) then
            excl_ratio = 0.0_wp
            if (half > 0.0_wp) excl_ratio = excl_radius/(half*sqrt(3.0_wp))
         end if

         crosses = lsf0_anchor /= 0.0_wp .and. lsf0*lsf0_anchor <= 0.0_wp
         if (crosses) then
            rho_centre = norm2(centre - anchor)
            rho_cross = max(0.0_wp, rho_centre - excl_radius)
            if (rho_cross < rho_hit) then
               rho_hit = rho_cross
               rho_cap = min(rho_cap, sqrt(rho_hit*rho_hit + rho2_slack))
               if (self%n_tightenings < max_logged_tightenings) then
                  self%n_tightenings = self%n_tightenings + 1
                  self%dbg_tighten_box(self%n_tightenings) = self%n_boxes_visited
                  self%dbg_tighten_raw(self%n_tightenings) = rho_centre
                  self%dbg_tighten_hit(self%n_tightenings) = rho_hit
                  self%dbg_tighten_cap(self%n_tightenings) = rho_cap
               end if
            end if
         end if

         ! Certified surface-free: the exclusion ball around the centre
         ! swallows the box, so no zero of S lies inside it.
         if (excl_radius >= half*sqrt(3.0_wp)) then
            if (trace) then
               self%dbg_free(depth) = self%dbg_free(depth) + 1
               write (action, "(a)") "surface-free: subtree pruned"
               call trace_row(self, depth, lat, leaf_depth, lsf0, excl_ratio, key, action)
            end if
            cycle
         end if

         if (depth >= leaf_depth) then
            if (self%n_surv >= size(self%surv_rho)) then
               if (self%n_surv >= self%max_survivors) then
                  call fatal_error(error, &
                                   "Octree branch search: survivor budget exhausted; "// &
                                   "seed_size is too small for the admissible radius.")
                  return
               end if
               call grow_survivor_arrays(self, &
                                         min(2*size(self%surv_rho), self%max_survivors))
            end if
            self%n_surv = self%n_surv + 1
            self%surv_lat(:, self%n_surv) = lat
            self%surv_centre(:, self%n_surv) = centre
            self%surv_rho(self%n_surv) = key
            if (trace) then
               self%dbg_kept(depth) = self%dbg_kept(depth) + 1
               write (action, "(a,i0)") "leaf: survivor ", self%n_surv
               call trace_row(self, depth, lat, leaf_depth, lsf0, excl_ratio, key, action)
            end if
            ! Estimated distance from the anchor to the surface *through* this
            ! box. The offset is the certified radius, not `abs(lsf0)`: the
            ! probe contract gives units only to the radius, whereas the level
            ! set value is a value whose sign alone is meaningful. The two
            ! coincide for a 1-Lipschitz LSF, but an LSF that must divide by a
            ! gradient bound (or works in density units) would rank nonsense.
            self%surv_rho_est(self%n_surv) = norm2(centre - anchor)
            if (lsf0*lsf0_anchor >= 0.0_wp) then
               self%surv_rho_est(self%n_surv) = self%surv_rho_est(self%n_surv) + excl_radius
            else
               self%surv_rho_est(self%n_surv) = self%surv_rho_est(self%n_surv) - excl_radius
            end if
            cycle
         end if

         if (trace) self%dbg_split(depth) = self%dbg_split(depth) + 1
         n_pushed = 0
         do iz = 0, 1
            do iy = 0, 1
               do ix = 0, 1
                  child_lat = 2*lat + [ix, iy, iz]
                  child_key = box_rho_lower(self, depth + 1, child_lat)
                  if (child_key > rho_cap) cycle
                  n_pushed = n_pushed + 1
                  call heap_push(self, child_key, depth + 1, child_lat, ok)
                  if (trace) self%dbg_heap_peak = max(self%dbg_heap_peak, self%heap_size)
                  if (.not. ok) then
                     call fatal_error(error, &
                                      "Octree branch search: heap capacity exhausted")
                     return
                  end if
               end do
            end do
         end do

         if (trace) then
            if (n_pushed == 8) then
               write (action, "(a)") "split into 8"
            else
               write (action, "(a,i0,a)") "split ", n_pushed, &
                  "/8, rest beyond rho_max"
            end if
            call trace_row(self, depth, lat, leaf_depth, lsf0, excl_ratio, key, action)
         end if
      end do

      if (trace) call self%trace_tbl%separator()

      self%rho_max_final = rho_cap
      self%dbg_n_surv_raw = self%n_surv

      ! Survivors were collected against a cap that may have tightened after
      ! they were accepted; drop the ones the final cap rules out.
      call drop_survivors_beyond(self, rho_cap)

      call collect_seeds(self, leaf_depth)

      if (trace) call self%report(anchor)
   end subroutine octree_run

   !> Discard survivors whose distance lower bound exceeds the final cap
   !> @param[inout] self    Search instance
   !> @param[in]    rho_cap Final admissible radius
   pure subroutine drop_survivors_beyond(self, rho_cap)
      class(moist_math_octree_branch_type), intent(inout) :: self
      real(wp), intent(in) :: rho_cap

      integer :: i, n_kept

      n_kept = 0
      do i = 1, self%n_surv
         if (self%surv_rho(i) > rho_cap) cycle
         n_kept = n_kept + 1
         if (n_kept /= i) then
            self%surv_lat(:, n_kept) = self%surv_lat(:, i)
            self%surv_centre(:, n_kept) = self%surv_centre(:, i)
            self%surv_rho(n_kept) = self%surv_rho(i)
            self%surv_rho_est(n_kept) = self%surv_rho_est(i)
         end if
      end do
      self%n_surv = n_kept
   end subroutine drop_survivors_beyond

!* ================================================================================= *!
!*                              Survivors to seed points                             *!
!* ================================================================================= *!

   !> Turn surviving leaves into seed points
   !>
   !> In `per_leaf` mode every survivor is a seed. In `cluster` mode a survivor
   !> is a seed when none of its 26 lattice neighbours present in the survivor
   !> set is closer to the anchor -- the discrete local minima of the distance
   !> field restricted to the survivors. That is one seed per basin rather than
   !> one per box, and no connected patch can be skipped: the smallest survivor
   !> of a connected group is a local minimum of that group by construction.
   !>
   !> The full 26-neighbourhood matters. Survivors are not a surface but a slab
   !> roughly two cells thick (a leaf survives whenever the surface passes
   !> within its circumradius), and across a slab the six face neighbours alone
   !> leave a cell looking like a minimum whenever the descent direction runs
   !> diagonally -- which on a curved patch is most of them. Testing against
   !> face neighbours only produced ~100 seeds for a patch holding one minimum.
   !>
   !> Ties are broken by lattice key so a symmetric pair of neighbours produces
   !> one seed rather than two or none.
   !>
   !> @param[inout] self       Search instance
   !> @param[in]    leaf_depth Depth all survivors share
   subroutine collect_seeds(self, leaf_depth)
      class(moist_math_octree_branch_type), intent(inout) :: self
      integer, intent(in) :: leaf_depth

      integer :: i, dx, dy, dz, neighbour, lattice_edge
      integer :: neighbour_lat(3)
      logical :: is_minimum

      self%n_seeds = 0
      if (self%n_surv <= 0) return

      if (self%seed_mode == octree_seed_per_leaf) then
         self%n_seeds = self%n_surv
         self%seeds(:, 1:self%n_seeds) = self%surv_centre(:, 1:self%n_surv)
         return
      end if

      do i = 1, self%n_surv
         self%surv_key(i) = lattice_key(self, self%surv_lat(:, i))
         self%surv_key_sorted(i) = self%surv_key(i)
      end do
      call heapsort_keys(self%n_surv, self%surv_key_sorted, self%surv_order)

      lattice_edge = 2**leaf_depth

      do i = 1, self%n_surv
         is_minimum = .true.
         neighbour_loop: do dz = -1, 1
            do dy = -1, 1
               do dx = -1, 1
                  if (dx == 0 .and. dy == 0 .and. dz == 0) cycle
                  neighbour_lat = self%surv_lat(:, i) + [dx, dy, dz]
                  if (any(neighbour_lat < 0) .or. any(neighbour_lat >= lattice_edge)) cycle
                  neighbour = find_survivor(self, lattice_key(self, neighbour_lat))
                  if (neighbour <= 0) cycle
                  if (precedes(self, neighbour, i)) then
                     is_minimum = .false.
                     exit neighbour_loop
                  end if
               end do
            end do
         end do neighbour_loop

         if (.not. is_minimum) cycle
         self%n_seeds = self%n_seeds + 1
         self%seeds(:, self%n_seeds) = self%surv_centre(:, i)
      end do
   end subroutine collect_seeds

   !> Whether survivor `a` sorts strictly before survivor `b`
   !>
   !> Estimated surface distance first (see `surv_rho_est`), with the lattice
   !> key as the tie-break, so the order is total and two exactly-tied
   !> neighbours cannot each rule the other out.
   !>
   !> @param[in] self Search instance
   !> @param[in] a    First survivor index
   !> @param[in] b    Second survivor index
   !> @returns        `.true.` when `a` precedes `b`
   pure function precedes(self, a, b) result(before)
      class(moist_math_octree_branch_type), intent(in) :: self
      integer, intent(in) :: a, b
      logical :: before

      if (self%surv_rho_est(a) /= self%surv_rho_est(b)) then
         before = self%surv_rho_est(a) < self%surv_rho_est(b)
      else
         before = self%surv_key(a) < self%surv_key(b)
      end if
   end function precedes

   !> Pack a lattice position into a single sortable integer
   !> @param[in] self Search instance
   !> @param[in] lat  Lattice position at the leaf depth
   !> @returns   key  Packed key
   pure function lattice_key(self, lat) result(key)
      class(moist_math_octree_branch_type), intent(in) :: self
      integer, intent(in) :: lat(3)
      integer(int64) :: key

      key = int(lat(1), int64) &
            + self%lattice_base*(int(lat(2), int64) &
                                 + self%lattice_base*int(lat(3), int64))
   end function lattice_key

   !> Sort `n` keys ascending, carrying an index permutation alongside
   !>
   !> Heapsort: in place, no recursion, and no worst case that a pathological
   !> survivor layout could trigger.
   !>
   !> @param[in]    n     Number of entries to sort
   !> @param[inout] keys  Keys, sorted ascending on exit
   !> @param[out]   order Permutation mapping sorted position to input position
   pure subroutine heapsort_keys(n, keys, order)
      integer, intent(in) :: n
      integer(int64), intent(inout) :: keys(:)
      integer, intent(out) :: order(:)

      integer :: i

      do i = 1, n
         order(i) = i
      end do

      do i = n/2, 1, -1
         call sift_down(i, n, keys, order)
      end do
      do i = n, 2, -1
         call swap_entries(1, i, keys, order)
         call sift_down(1, i - 1, keys, order)
      end do
   end subroutine heapsort_keys

   !> Restore the max-heap property below `start` within `bound`
   !> @param[in]    start Position to sift down from
   !> @param[in]    bound Last position belonging to the heap
   !> @param[inout] keys  Key array
   !> @param[inout] order Permutation moved alongside the keys
   pure subroutine sift_down(start, bound, keys, order)
      integer, intent(in) :: start, bound
      integer(int64), intent(inout) :: keys(:)
      integer, intent(inout) :: order(:)

      integer :: root, child

      root = start
      do
         child = 2*root
         if (child > bound) exit
         if (child < bound) then
            if (keys(child + 1) > keys(child)) child = child + 1
         end if
         if (keys(root) >= keys(child)) exit
         call swap_entries(root, child, keys, order)
         root = child
      end do
   end subroutine sift_down

   !> Exchange two entries of a key array and its permutation
   !> @param[in]    i1    First position
   !> @param[in]    i2    Second position
   !> @param[inout] keys  Key array
   !> @param[inout] order Permutation moved alongside the keys
   pure subroutine swap_entries(i1, i2, keys, order)
      integer, intent(in) :: i1, i2
      integer(int64), intent(inout) :: keys(:)
      integer, intent(inout) :: order(:)

      integer(int64) :: tmp_key
      integer :: tmp_idx

      tmp_key = keys(i1)
      keys(i1) = keys(i2)
      keys(i2) = tmp_key

      tmp_idx = order(i1)
      order(i1) = order(i2)
      order(i2) = tmp_idx
   end subroutine swap_entries

   !> Find the survivor holding a given lattice key
   !> @param[in] self Search instance
   !> @param[in] key  Packed lattice key to look for
   !> @returns   idx  Survivor index, or 0 when that cell did not survive
   pure function find_survivor(self, key) result(idx)
      class(moist_math_octree_branch_type), intent(in) :: self
      integer(int64), intent(in) :: key
      integer :: idx

      integer :: low, high, mid

      idx = 0
      low = 1
      high = self%n_surv
      do while (low <= high)
         mid = (low + high)/2
         if (self%surv_key_sorted(mid) == key) then
            idx = self%surv_order(mid)
            return
         else if (self%surv_key_sorted(mid) < key) then
            low = mid + 1
         else
            high = mid - 1
         end if
      end do
   end function find_survivor

!* ================================================================================= *!
!*                                   Debug report                                    *!
!* ================================================================================= *!

   !> Print a per-anchor trace of the last run
   !>
   !> Three tables, because three questions come up when a search misbehaves.
   !> The depth table answers "where did the boxes go?" -- a level whose boxes
   !> are nearly all split rather than certified surface-free is a level where
   !> the exclusion radius is buying nothing. The tightening table answers "did
   !> the search shrink its own ball, and how early?" -- the first sign change
   !> should arrive within the first handful of boxes, and the gap between its
   !> two distance columns is what the exclusion radius bought. The seed table
   !> answers "what came out, and is it plausible?" -- seeds should sit at a
   !> distance close to the certified radius, one per branch.
   !>
   !> Callers inside an OpenMP region should wrap this in a critical section;
   !> the report is one contiguous block per anchor so it stays readable.
   !>
   !> @param[inout] self   Search instance holding the last run's bookkeeping
   !> @param[in]    anchor Anchor the run was for
   subroutine octree_report(self, anchor)
      class(moist_math_octree_branch_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)

      type(prettylistprinter) :: plp
      integer :: depth, i, total_popped, total_free, total_split, total_kept
      real(wp) :: rho_seed

      write (self%unit, "(a)") ""
      write (self%unit, "(4x,a,f13.5)") "rho_max (certified)    ", self%rho_max_final
      write (self%unit, "(a)") ""

      !> Where the boxes went, level by level
      !>
      !> Collected only on a debug run, so a report asked for after an ordinary
      !> one says so rather than printing a table of zeros.
      if (.not. self%dbg_valid) then
         write (self%unit, "(4x,a)") &
            "Per-depth tallies were not collected: set `debug` before `run`."
         write (self%unit, "(a)") ""
      else
         plp = new_prettylistprinter( &
               widths=[7, 12, 12, 12, 12, 12], &
               headers=[character(len=12) :: "depth", "edge", "popped", &
                        "surf-free", "split", "kept"], &
               unit=self%unit, offset=4, column_gap=1)
         call plp%print_header()
         call plp%separator()

         total_popped = 0
         total_free = 0
         total_split = 0
         total_kept = 0
         do depth = 0, self%dbg_leaf_depth
            total_popped = total_popped + self%dbg_popped(depth)
            total_free = total_free + self%dbg_free(depth)
            total_split = total_split + self%dbg_split(depth)
            total_kept = total_kept + self%dbg_kept(depth)
            if (self%dbg_popped(depth) == 0) cycle
            call plp%begin_row()
            call plp%add(depth)
            call plp%add(2.0_wp*self%dbg_rho_start/real(2**depth, wp), fmt="F12.4")
            call plp%add(self%dbg_popped(depth))
            call plp%add(self%dbg_free(depth))
            call plp%add(self%dbg_split(depth))
            call plp%add(self%dbg_kept(depth))
            call plp%end_row()
         end do
         call plp%separator()
         call plp%begin_row()
         call plp%add("total")
         call plp%skip()
         call plp%add(total_popped)
         call plp%add(total_free)
         call plp%add(total_split)
         call plp%add(total_kept)
         call plp%end_row()
         write (self%unit, "(a)") ""
      end if

      !> How the admissible radius came down
      if (self%n_tightenings > 0) then
         plp = new_prettylistprinter( &
               widths=[10, 14, 14, 14, 14], &
               headers=[character(len=14) :: "event", "at box", "centre at", &
                        "surface within", "new rho_max"], &
               unit=self%unit, offset=4, column_gap=1)
         call plp%print_header()
         call plp%separator()
         do i = 1, self%n_tightenings
            call plp%begin_row()
            call plp%add(i)
            call plp%add(self%dbg_tighten_box(i))
            call plp%add(self%dbg_tighten_raw(i), fmt="F14.5")
            call plp%add(self%dbg_tighten_hit(i), fmt="F14.5")
            call plp%add(self%dbg_tighten_cap(i), fmt="F14.5")
            call plp%end_row()
         end do
         call plp%separator()
         write (self%unit, "(a)") ""
      else
         write (self%unit, "(4x,a)") &
            "No sign change met: the run kept the radius it started with."
         write (self%unit, "(a)") ""
      end if

      !> What came out
      write (self%unit, "(4x,a,i13)") "boxes examined         ", self%n_boxes_visited
      write (self%unit, "(4x,a,i13)") "peak heap size         ", self%dbg_heap_peak
      write (self%unit, "(4x,a,i13)") "survivors (raw)        ", self%dbg_n_surv_raw
      write (self%unit, "(4x,a,i13)") "survivors (in radius)  ", self%n_surv
      write (self%unit, "(4x,a,i13)") "seeds                  ", self%n_seeds
      write (self%unit, "(a)") ""

      if (self%n_seeds > 0) then
         plp = new_prettylistprinter( &
               widths=[8, 13, 13, 13, 13], &
               headers=[character(len=13) :: "seed", "x", "y", "z", "rho"], &
               unit=self%unit, offset=4, column_gap=1)
         call plp%print_header()
         call plp%separator()
         do i = 1, self%n_seeds
            rho_seed = norm2(self%seeds(:, i) - anchor)
            call plp%begin_row()
            call plp%add(i)
            call plp%add(self%seeds(1, i), fmt="F13.5")
            call plp%add(self%seeds(2, i), fmt="F13.5")
            call plp%add(self%seeds(3, i), fmt="F13.5")
            call plp%add(rho_seed, fmt="F13.5")
            call plp%end_row()
         end do
         call plp%separator()
      end if
      write (self%unit, "(a)") ""
   end subroutine octree_report

!* ================================================================================= *!
!*                                    Walk trace                                     *!
!* ================================================================================= *!

   !> Open the per-anchor trace with the run's inputs
   !>
   !> @param[in]    self       Search instance
   !> @param[in]    anchor     Anchor the run is for
   !> @param[in]    leaf_depth Deepest level this run can reach
   subroutine trace_banner(self, anchor, leaf_depth)
      class(moist_math_octree_branch_type), intent(in) :: self
      real(wp), intent(in) :: anchor(3)
      integer, intent(in) :: leaf_depth

      write (self%unit, "(a)") ""
      write (self%unit, "(2x,a)") repeat("=", 74)
      write (self%unit, "(2x,a,i0,a,3(1x,f10.5))") &
         "== Octree branch search, point ", self%point_id, ", anchor", anchor
      write (self%unit, "(2x,a)") repeat("=", 74)

      write (self%unit, "(4x,a,es13.5)") "S(anchor)              ", self%dbg_lsf0_anchor
      write (self%unit, "(4x,a,f13.5)") "rho_max (start)        ", self%dbg_rho_start
      write (self%unit, "(4x,a,f13.5)") "rho^2 slack            ", self%dbg_rho2_slack
      write (self%unit, "(4x,a,f13.5)") "seed size (Bohr)       ", self%seed_size
      write (self%unit, "(4x,a,i13)") "leaf depth             ", leaf_depth
      write (self%unit, "(4x,a,f13.5)") "leaf edge (Bohr)       ", &
         2.0_wp*self%root_half/real(2**leaf_depth, wp)
      write (self%unit, "(a)") ""
   end subroutine trace_banner

   !> Octant this box descended into at a given level of the tree
   !>
   !> The lattice position at depth `depth` is the concatenation of the octant
   !> choices made on the way down: each level appends one bit per axis. Level
   !> `level` therefore reads bit `depth - level` of each component, and the
   !> three bits pack into the familiar 1..8 octant index (+x fastest).
   !>
   !> @param[in] lat   Lattice position at `depth`
   !> @param[in] depth Depth the box sits at
   !> @param[in] level Level to read, 1 .. depth
   !> @returns   oct   Octant index 1..8
   pure function child_octant(lat, depth, level) result(oct)
      integer, intent(in) :: lat(3)
      integer, intent(in) :: depth
      integer, intent(in) :: level
      integer :: oct

      integer :: bit

      bit = depth - level
      oct = 1 + ibits(lat(1), bit, 1) &
            + 2*ibits(lat(2), bit, 1) &
            + 4*ibits(lat(3), bit, 1)
   end function child_octant

   !> Open the walk trace for one anchor
   !>
   !> One column per tree level, two characters wide, so the depth a box sits
   !> at is simply the column its octant digit lands in: reading down the table
   !> shows the search descending and backing out again.
   !>
   !> @param[inout] self       Search instance
   !> @param[in]    leaf_depth Deepest level this run can reach
   subroutine trace_open(self, leaf_depth)
      class(moist_math_octree_branch_type), intent(inout) :: self
      integer, intent(in) :: leaf_depth

      integer :: widths(leaf_depth + 5)
      character(len=6) :: headers(leaf_depth + 5)
      integer :: level

      widths(1) = 8
      headers(1) = "point"
      do level = 1, leaf_depth
         widths(level + 1) = 2
         write (headers(level + 1), "(i0)") level
      end do
      widths(leaf_depth + 2) = 12
      headers(leaf_depth + 2) = "S"
      widths(leaf_depth + 3) = 8
      headers(leaf_depth + 3) = "excl/r"
      widths(leaf_depth + 4) = 9
      headers(leaf_depth + 4) = "rho_lo"
      widths(leaf_depth + 5) = trace_action_width
      headers(leaf_depth + 5) = "action"

      self%trace_tbl = new_prettylistprinter(widths=widths, headers=headers, &
                                             unit=self%unit, offset=4, column_gap=1)
      call self%trace_tbl%print_header()
      call self%trace_tbl%separator()
   end subroutine trace_open

   !> Emit one row of the walk trace
   !>
   !> @param[inout] self       Search instance
   !> @param[in]    depth      Depth of the box just examined
   !> @param[in]    lat        Its lattice position
   !> @param[in]    leaf_depth Deepest level this run can reach
   !> @param[in]    lsf0       Level set value at the box centre
   !> @param[in]    excl_ratio Exclusion radius over the box circumradius; the
   !>                          box is certified surface-free once this reaches 1
   !> @param[in]    rho_lo     Distance lower bound of the box
   !> @param[in]    action     What the search decided to do with it
   subroutine trace_row(self, depth, lat, leaf_depth, lsf0, excl_ratio, rho_lo, action)
      class(moist_math_octree_branch_type), intent(inout) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      integer, intent(in) :: leaf_depth
      real(wp), intent(in) :: lsf0
      real(wp), intent(in) :: excl_ratio
      real(wp), intent(in) :: rho_lo
      character(len=*), intent(in) :: action

      integer :: level

      call self%trace_tbl%begin_row()
      call self%trace_tbl%add(self%point_id)
      do level = 1, leaf_depth
         if (level <= depth) then
            call self%trace_tbl%add(child_octant(lat, depth, level))
         else
            call self%trace_tbl%skip()
         end if
      end do
      call self%trace_tbl%add(lsf0, fmt="ES12.3")
      call self%trace_tbl%add(excl_ratio, fmt="F8.3")
      call self%trace_tbl%add(rho_lo, fmt="F9.4")
      call self%trace_tbl%add(trim(action))
      call self%trace_tbl%end_row()
   end subroutine trace_row

end module moist_math_solver_octree_branch
