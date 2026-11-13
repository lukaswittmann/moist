!> Certified octree surface coverage seeds for the batched DROP projection
!>
!> Searched level-synchronously: every depth's surviving boxes are probed in
!> one batched call, so one `evaluate_batch` per depth replaces the one LSF
!> call per box of a best-first search
!>
!> Certificate unchanged -- a box is discarded only when
!>   - its distance lower bound exceeds the admissible cap `rho_cap`
!>   - the probe's surface-free ball around its center swallows it
!>   - the optional region certificate excludes it
!>
!> Seeds come from sign changes, not from every surviving leaf:
!>   - a surviving leaf has evidence when its center and a face neighbour's
!>     center (one leaf edge away) lie on opposite sides of the surface;
!>     neighbours not probed at the leaf depth are probed in one extra batch
!>   - every sign-change segment between face neighbours gives one seed at its
!>     crossing, emitted once although both of its leaves see it
!>   - a surviving leaf without evidence is refined `refine_levels` times: its
!>     eight children are probed in one batch over all such leaves, and a child
!>     on the other side of the leaf center gives a seed on that segment; a
!>     leaf whose children show no sign change is dropped, as is every leaf
!>     without evidence when `refine_levels` is 0
!>   - a crossing is the linear interpolation of the two probe values
!>   - every seed lies within `reach` of the capped ball, so a later, tighter
!>     ball around the same anchor can reuse the seeds instead of searching again
!>
!> `rho_cap` tightens once per depth from that depth's sign changes,
!> `rho_min <= ||c - anchor|| - r` for a probed center `c` of opposite sign
!> and surface-free radius `r`; a best-first search tightens it after every
!> box, so the level-synchronous search prunes slightly less and may probe
!> more boxes
!>
!> A center beyond the level set reach takes the anchor's value and the
!> LSF's `far_exclusion_radius()`, so it never tightens the cap, and its
!> unreached flag keeps it out of every sign-change segment
!>
!> The probe stays the single abstraction: [[octree_run]] takes any batched
!> probe, the seeder entry plugs in the level set
module moist_cavity_drop_projection_seed_octree
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seed_request_type
   implicit none(type, external)
   private

   public :: octree_seeder_type, octree_probe_batch_interface, octree_exclude_region_interface

   !> Supported depth limit, keeping integer lattice coordinates representable
   integer, parameter :: max_addressable_depth = 20
   !> |S| at which a sphere-traced point counts as on the surface
   real(wp), parameter :: trace_feasible_tol = 1.0e-6_wp
   !> Sphere-trace steps before the ray counts as stalled
   integer, parameter :: trace_max_steps = 48
   !> Gradient norm below which the trace direction is meaningless
   real(wp), parameter :: trace_grad_floor = 1.0e-8_wp
   !> Longest single step of the sphere trace (Bohr)
   real(wp), parameter :: trace_step_limit = 15.0_wp
   !> Refinement levels below the leaf depth for a surviving leaf without evidence
   integer, parameter :: default_refine_levels = 1
   !> Face-neighbour lattice steps, the three positive directions first
   integer, parameter :: face_step(3, 6) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1, &
                                                    -1, 0, 0, 0, -1, 0, 0, 0, -1], [3, 6])
   !> Spatial hash multipliers of the leaf table
   integer(int64), parameter :: hash_mult(3) = [73856093_int64, 19349663_int64, 83492791_int64]

   abstract interface
      !> Probe the level set at a batch of points
      !>
      !> @param[in]    points  Evaluation points (3, n)
      !> @param[out]   lsf0    Level set values (n); a point the probe cannot evaluate
      !>                       reports a value of the anchor's sign, so it never tightens
      !>                       the cap
      !> @param[out]   radius  Surface-free radius around each point; zero certifies nothing
      !> @param[out]   reached Whether the probe evaluated the point; an unreached point
      !>                       never ends a sign-change segment
      !> @param[inout] context Caller context forwarded unchanged
      !> @param[out]   error   Probe failure (build-fatal)
      subroutine octree_probe_batch_interface(points, lsf0, radius, reached, context, error)
         import :: wp, error_type
         implicit none(type, external)
         real(wp), intent(in), contiguous :: points(:, :)
         real(wp), intent(out), contiguous :: lsf0(:)
         real(wp), intent(out), contiguous :: radius(:)
         logical, intent(out), contiguous :: reached(:)
         class(*), intent(inout) :: context
         type(error_type), allocatable, intent(out) :: error
      end subroutine octree_probe_batch_interface

      !> Certify a whole box irrelevant to the caller's root search
      !>
      !> @param[in] center     Box midpoint
      !> @param[in] half_width Cube half-width
      !> @param[in] context    Caller context forwarded unchanged
      function octree_exclude_region_interface(center, half_width, context) result(excluded)
         import :: wp
         implicit none(type, external)
         real(wp), intent(in) :: center(3)
         real(wp), intent(in) :: half_width
         class(*), intent(in) :: context
         logical :: excluded
      end function octree_exclude_region_interface
   end interface

   !> Level set probe context of the seeder entry
   type :: octree_lsf_context_type
      !> Thread-local level set
      class(moist_cavity_drop_lsf_type), pointer :: lsf => null()
      !> Anchor of the search
      real(wp) :: anchor(3) = 0.0_wp
      !> Level set value at the anchor, reported at unreached probe points
      real(wp) :: lsf0_anchor = 0.0_wp
   end type octree_lsf_context_type

   !> Level-synchronous octree seeder; one instance per thread
   type, extends(projection_seeder_type) :: octree_seeder_type

      !* ---------------------------- Search configuration ---------------------------- *!

      !> Edge length of the leaves whose sign changes become seeds (Bohr)
      real(wp) :: seed_size = 0.2_wp
      !> Tree-descent boxes probed per anchor before the search gives up, counted at
      !> every depth; face-neighbour and refinement probes are not counted, the survivor
      !> budget bounds them (six per surviving leaf, eight per refined parent and level)
      integer :: max_boxes = 500000
      !> Surviving leaves per anchor
      integer :: max_survivors = 200000
      !> Depth limit
      integer :: max_depth = 12
      !> Refinement levels below the leaf depth for a surviving leaf without evidence
      integer :: refine_levels = default_refine_levels
      !> Per-anchor report after each search
      logical :: debug = .false.
      !> Report unit
      integer :: unit = output_unit
      !> Anchor id shown in the report
      integer :: point_id = 0

      !* ---------------------------------- Results ----------------------------------- *!

      !> Admissible radius after self-tightening
      real(wp) :: rho_max_final = 0.0_wp
      !> Tree-descent boxes probed by the last search, the count `max_boxes` bounds
      integer :: n_boxes_visited = 0
      !> Boxes discarded by the region certificate
      integer :: n_regions_excluded = 0

      !* ------------------------------- Working state -------------------------------- *!

      !> Current and next frontier: lattice position and distance lower bound
      integer, allocatable :: cur_lat(:, :), nxt_lat(:, :)
      real(wp), allocatable :: cur_key(:), nxt_key(:)
      !> Probe inputs and outputs of the current depth, leaf face neighbours appended
      real(wp), allocatable :: centers(:, :), probe0(:), probe_r(:)
      logical, allocatable :: probe_ok(:)
      !> Leaf passed the cap, the surface-free ball and the region certificate
      logical, allocatable :: survivor(:)
      !> Open-addressing table of the leaf lattice: key and box index (0 = empty)
      integer(int64), allocatable :: table_key(:)
      integer, allocatable :: table_box(:)
      !> Slots of the table in use, a power of two
      integer :: table_size = 0
      !> Sign-change segments: end points (3, :) and their values
      real(wp), allocatable :: seg_a(:, :), seg_b(:, :), seg_va(:), seg_vb(:)
      !> Refinement parents: centers (3, :) and values
      real(wp), allocatable :: ref_c(:, :), ref_v(:)
      !> Refinement children: centers (3, :), values, radii, probe status
      real(wp), allocatable :: kid_c(:, :), kid_v(:), kid_r(:)
      logical, allocatable :: kid_ok(:)
      !> Root half-width of the search in progress
      real(wp) :: root_half = 0.0_wp
      !> Whether `init` ran
      logical :: initialized = .false.

      !* ----------------------------- Debug bookkeeping ------------------------------ *!

      !> Per-depth tallies: probed, surface-free, region-excluded, split, surviving leaves
      integer :: dbg_probed(0:max_addressable_depth) = 0
      integer :: dbg_free(0:max_addressable_depth) = 0
      integer :: dbg_excluded(0:max_addressable_depth) = 0
      integer :: dbg_split(0:max_addressable_depth) = 0
      integer :: dbg_leaves(0:max_addressable_depth) = 0
      !> Cap after each depth
      real(wp) :: dbg_cap(0:max_addressable_depth) = 0.0_wp
      !> Face neighbours probed after the leaf depth
      integer :: dbg_neighbours = 0
      !> Surviving leaves with sign-change evidence
      integer :: dbg_evidence = 0
      !> Surviving leaves without evidence (refined unless `refine_levels` is 0), and the
      !> children probed for them
      integer :: dbg_refined = 0
      integer :: dbg_children = 0
      !> Seeds from face-neighbour segments and from refinement children
      integer :: dbg_face_seeds = 0
      integer :: dbg_child_seeds = 0
      !> Leaf depth and starting radius of the last search
      integer :: dbg_leaf_depth = 0
      real(wp) :: dbg_rho_start = 0.0_wp
   contains
      !> Validate the configuration and size the scratch
      procedure :: init => octree_init
      !> Search one anchor with any batched probe
      procedure :: run => octree_run
      !> Seeder entry: search one anchor with the level set probe
      procedure :: seed => octree_seed
      !> Print the last search's tallies
      procedure :: report => octree_report
      !> Release the scratch
      procedure :: destroy => octree_destroy
   end type octree_seeder_type

contains

   !* ================================================================================= *!
   !*                                     Lifecycle                                     *!
   !* ================================================================================= *!

   !> Validate the configuration and size the scratch
   !>
   !> @param[inout] self          Seeder
   !> @param[out]   error         Invalid configuration
   !> @param[in]    seed_size     Leaf edge length (Bohr)
   !> @param[in]    max_boxes     Tree-descent box budget per anchor
   !> @param[in]    max_survivors Survivor budget per anchor
   !> @param[in]    max_depth     Depth limit
   !> @param[in]    refine_levels Refinement levels of a leaf without evidence
   subroutine octree_init(self, error, seed_size, max_boxes, max_survivors, max_depth, refine_levels)
      class(octree_seeder_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: seed_size
      integer, intent(in), optional :: max_boxes
      integer, intent(in), optional :: max_survivors
      integer, intent(in), optional :: max_depth
      integer, intent(in), optional :: refine_levels

      integer, parameter :: initial_capacity = 1024

      if (present(seed_size)) self%seed_size = seed_size
      if (present(max_boxes)) self%max_boxes = max_boxes
      if (present(max_survivors)) self%max_survivors = max_survivors
      if (present(max_depth)) self%max_depth = max_depth
      if (present(refine_levels)) self%refine_levels = refine_levels

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
      if (self%refine_levels < 0) then
         call fatal_error(error, "Octree branch search: refine_levels must not be negative")
         return
      end if
      call self%destroy()

      call grow_frontier(self%cur_lat, self%cur_key, min(initial_capacity, self%max_boxes), 0)
      call grow_frontier(self%nxt_lat, self%nxt_key, min(initial_capacity, self%max_boxes), 0)
      call reserve_boxes(self, size(self%cur_key), 0)
      self%initialized = .true.
   end subroutine octree_init

   !> Release the scratch
   !>
   !> @param[inout] self Seeder
   subroutine octree_destroy(self)
      class(octree_seeder_type), intent(inout) :: self

      if (allocated(self%cur_lat)) deallocate (self%cur_lat)
      if (allocated(self%nxt_lat)) deallocate (self%nxt_lat)
      if (allocated(self%cur_key)) deallocate (self%cur_key)
      if (allocated(self%nxt_key)) deallocate (self%nxt_key)
      if (allocated(self%centers)) deallocate (self%centers)
      if (allocated(self%probe0)) deallocate (self%probe0)
      if (allocated(self%probe_r)) deallocate (self%probe_r)
      if (allocated(self%probe_ok)) deallocate (self%probe_ok)
      if (allocated(self%survivor)) deallocate (self%survivor)
      if (allocated(self%table_key)) deallocate (self%table_key)
      if (allocated(self%table_box)) deallocate (self%table_box)
      if (allocated(self%seg_a)) deallocate (self%seg_a, self%seg_b, self%seg_va, self%seg_vb)
      if (allocated(self%ref_c)) deallocate (self%ref_c, self%ref_v)
      if (allocated(self%kid_c)) deallocate (self%kid_c, self%kid_v, self%kid_r, self%kid_ok)
      self%table_size = 0
      self%n_points = 0
      self%n_boxes_visited = 0
      self%n_regions_excluded = 0
      self%initialized = .false.
   end subroutine octree_destroy

   !* ================================================================================= *!
   !*                                    Seeder entry                                   *!
   !* ================================================================================= *!

   !> Cover the admissible surface region of one anchor and emit crossing seeds
   !>
   !> - the starting radius is the request's `rho_start` when a verified foot
   !>   gave one, else `sqrt(rho^2 + slack)` for the distance `rho` of the
   !>   surface point a sphere trace along the normal reaches, else the
   !>   request's `rho_fallback`
   !> - a radius from a foot or the trace bounds the whole window and is used
   !>   as it is, however large
   !> - with no slack, an anchor on the surface (exactly, or within the trace
   !>   tolerance, so the starting radius is zero) is its own single seed
   !> - `nearest_first` orders the seeds by their distance to the anchor
   !>
   !> @param[inout] self      Seeder (initialized)
   !> @param[inout] lsf       Thread-local level set
   !> @param[in]    request   Anchor request
   !> @param[out]   error     LSF failure or exhausted budget (build-fatal)
   subroutine octree_seed(self, lsf, request, error)
      class(octree_seeder_type), intent(inout), target :: self
      class(moist_cavity_drop_lsf_type), intent(inout), target :: lsf
      type(seed_request_type), intent(in) :: request
      type(error_type), allocatable, intent(out) :: error

      type(octree_lsf_context_type) :: ctx
      real(wp) :: rho_start, rho_traced
      logical :: traced

      self%n_points = 0
      if (.not. self%initialized) then
         call fatal_error(error, "Octree branch search: run before init")
         return
      end if

      if (request%lsf0_anchor == 0.0_wp .and. request%rho2_slack == 0.0_wp) then
         call reset_tallies(self, 0.0_wp)
         call self%push(request%anchor)
         return
      end if

      ctx%lsf => lsf
      ctx%anchor = request%anchor
      ctx%lsf0_anchor = request%lsf0_anchor

      if (request%rho_start > 0.0_wp) then
         rho_start = request%rho_start
      else
         call trace_rho_bound(ctx, request, rho_traced, traced, error)
         if (allocated(error)) return
         if (traced) then
            rho_start = sqrt(rho_traced*rho_traced + request%rho2_slack)
         else
            rho_start = request%rho_fallback
         end if
      end if
      if (rho_start <= 0.0_wp) then
         call reset_tallies(self, 0.0_wp)
         call self%push(request%anchor)
         return
      end if

      self%point_id = request%index
      call self%run(request%anchor, request%lsf0_anchor, rho_start, request%rho2_slack, &
                    lsf_probe, ctx, error, exclude_region=lsf_exclude_region)
      if (allocated(error)) return
      if (request%nearest_first) call sort_points_by_distance(self, request%anchor)
   end subroutine octree_seed

   !> Order the seeds by distance to the anchor, ties in emission order
   !>
   !> @param[inout] self   Seeder holding the seeds
   !> @param[in]    anchor Anchor
   subroutine sort_points_by_distance(self, anchor)
      type(octree_seeder_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)

      real(wp), allocatable :: dist(:), sorted(:, :)
      integer, allocatable :: order(:), work(:)
      integer :: k, n

      n = self%n_points
      if (n < 2) return
      allocate (dist(n), order(n), work(n), sorted(3, n))
      do k = 1, n
         dist(k) = norm2(self%points(:, k) - anchor)
         order(k) = k
      end do
      call merge_sort_by_key(dist, order, work, 1, n)
      do k = 1, n
         sorted(:, k) = self%points(:, order(k))
      end do
      self%points(:, 1:n) = sorted
   end subroutine sort_points_by_distance

   !> Stable merge sort of an index list by a real key
   !>
   !> @param[in]    key   Sort keys, indexed by the list entries
   !> @param[inout] order Index list, sorted on `lo:hi`
   !> @param[inout] work  Scratch of the list's size
   !> @param[in]    lo    First entry
   !> @param[in]    hi    Last entry
   recursive subroutine merge_sort_by_key(key, order, work, lo, hi)
      real(wp), intent(in) :: key(:)
      integer, intent(inout) :: order(:)
      integer, intent(inout) :: work(:)
      integer, intent(in) :: lo, hi

      integer :: mid, i, j, k

      if (hi <= lo) return
      mid = (lo + hi)/2
      call merge_sort_by_key(key, order, work, lo, mid)
      call merge_sort_by_key(key, order, work, mid + 1, hi)
      i = lo
      j = mid + 1
      do k = lo, hi
         if (j > hi) then
            work(k) = order(i)
            i = i + 1
         else if (i > mid) then
            work(k) = order(j)
            j = j + 1
         else if (key(order(j)) < key(order(i))) then
            work(k) = order(j)
            j = j + 1
         else
            work(k) = order(i)
            i = i + 1
         end if
      end do
      order(lo:hi) = work(lo:hi)
   end subroutine merge_sort_by_key

   !* ================================================================================= *!
   !*                              Level-synchronous search                             *!
   !* ================================================================================= *!

   !> Cover unresolved regions within the admissible radius of one anchor
   !>
   !> - on return `points(:, 1:n_points)` holds the crossing seeds and
   !>   `rho_max_final` the radius actually certified
   !>
   !> @param[inout] self           Seeder (initialized)
   !> @param[in]    anchor         Anchor
   !> @param[in]    lsf0_anchor    Level set value at the anchor, for the sign test
   !> @param[in]    rho_max        Starting admissible radius (Bohr)
   !> @param[in]    rho2_slack     Squared-distance window, `rho_max^2 = rho_min^2 + slack`
   !> @param[in]    probe          Batched probe
   !> @param[inout] context        Probe context
   !> @param[out]   error          Exhausted budget, coarse depth limit, probe failure
   !> @param[in]    exclude_region Optional whole-box certificate
   subroutine octree_run(self, anchor, lsf0_anchor, rho_max, rho2_slack, probe, context, error, &
                         exclude_region)
      class(octree_seeder_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: lsf0_anchor
      real(wp), intent(in) :: rho_max
      real(wp), intent(in) :: rho2_slack
      procedure(octree_probe_batch_interface) :: probe
      class(*), intent(inout) :: context
      type(error_type), allocatable, intent(out) :: error
      procedure(octree_exclude_region_interface), optional :: exclude_region

      real(wp) :: rho_cap, rho_hit, half, child_key, rho_cross
      integer :: depth, leaf_depth, ncur, nnext, nleaves, k, ix, iy, iz, child_lat(3)
      character(len=256) :: buffer

      call reset_tallies(self, rho_max)

      if (.not. self%initialized) then
         call fatal_error(error, "Octree branch search: run before init")
         return
      end if
      if (rho_max <= 0.0_wp) then
         call fatal_error(error, "Octree branch search: rho_max must be positive")
         return
      end if
      ! A negative slack puts a NaN into the cap, after which nothing is bounded
      if (rho2_slack < 0.0_wp) then
         call fatal_error(error, "Octree branch search: rho2_slack must not be negative")
         return
      end if

      self%root_half = rho_max

      ! Depth at which the edge first drops to the seed size; stopping short of
      ! it would merge minima the configured resolution separates
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
      self%dbg_leaf_depth = leaf_depth
      ! A seed lies within one leaf edge of a surviving leaf's center, and that
      ! center within a half-diagonal of the capped ball
      self%reach = (1.0_wp + 0.5_wp*sqrt(3.0_wp))*2.0_wp*box_half(self, leaf_depth)

      rho_cap = rho_max
      rho_hit = huge(1.0_wp)
      ! An anchor on the surface is its own closest branch
      if (lsf0_anchor == 0.0_wp) then
         rho_hit = 0.0_wp
         rho_cap = min(rho_cap, sqrt(rho2_slack))
      end if

      ncur = 1
      self%cur_lat(:, 1) = 0
      self%cur_key(1) = 0.0_wp

      ! Children are queued against the cap of their parent's depth, which only
      ! tightens after the next probe, so every queued box is probed
      do depth = 0, leaf_depth
         if (ncur == 0) exit

         ! Within the box budget: the root is one box, and child generation
         ! refuses a frontier that would exceed it
         self%n_boxes_visited = self%n_boxes_visited + ncur
         self%dbg_probed(depth) = ncur

         half = box_half(self, depth)
         if (size(self%probe0) < ncur) call reserve_boxes(self, size(self%cur_key), 0)
         do k = 1, ncur
            self%centers(:, k) = anchor + box_offset(self, depth, self%cur_lat(:, k))
         end do
         call probe(self%centers(:, 1:ncur), self%probe0(1:ncur), self%probe_r(1:ncur), &
                    self%probe_ok(1:ncur), context, error)
         if (allocated(error)) return

         ! Sign changes of this depth tighten the cap for all of its boxes
         if (lsf0_anchor /= 0.0_wp) then
            do k = 1, ncur
               if (self%probe0(k)*lsf0_anchor > 0.0_wp) cycle
               rho_cross = max(0.0_wp, norm2(self%centers(:, k) - anchor) - self%probe_r(k))
               rho_hit = min(rho_hit, rho_cross)
            end do
            if (rho_hit < huge(1.0_wp)) rho_cap = min(rho_cap, sqrt(rho_hit*rho_hit + rho2_slack))
         end if
         self%dbg_cap(depth) = rho_cap

         nnext = 0
         nleaves = 0
         do k = 1, ncur
            self%survivor(k) = .false.
            if (self%cur_key(k) > rho_cap) cycle

            ! The surface-free ball swallows the box
            if (self%probe_r(k) >= half*sqrt(3.0_wp)) then
               self%dbg_free(depth) = self%dbg_free(depth) + 1
               cycle
            end if

            if (present(exclude_region)) then
               if (exclude_region(self%centers(:, k), half, context)) then
                  self%n_regions_excluded = self%n_regions_excluded + 1
                  self%dbg_excluded(depth) = self%dbg_excluded(depth) + 1
                  cycle
               end if
            end if

            ! The cap tightens no further after the leaf depth, so a leaf that
            ! passes it is searched for sign changes
            if (depth >= leaf_depth) then
               if (nleaves >= self%max_survivors) then
                  call fatal_error(error, &
                                   "Octree branch search: survivor budget exhausted; "// &
                                   "seed_size is too small for the admissible radius.")
                  return
               end if
               nleaves = nleaves + 1
               self%survivor(k) = .true.
               cycle
            end if

            self%dbg_split(depth) = self%dbg_split(depth) + 1
            do iz = 0, 1
               do iy = 0, 1
                  do ix = 0, 1
                     child_lat = 2*self%cur_lat(:, k) + [ix, iy, iz]
                     child_key = box_rho_lower(self, depth + 1, child_lat)
                     if (child_key > rho_cap) cycle
                     if (self%n_boxes_visited + nnext + 1 > self%max_boxes) then
                        call budget_error(error)
                        return
                     end if
                     if (nnext + 1 > size(self%nxt_key)) then
                        call grow_frontier(self%nxt_lat, self%nxt_key, &
                                           min(2*size(self%nxt_key), self%max_boxes), nnext)
                     end if
                     nnext = nnext + 1
                     self%nxt_lat(:, nnext) = child_lat
                     self%nxt_key(nnext) = child_key
                  end do
               end do
            end do
         end do

         if (depth >= leaf_depth) then
            self%dbg_leaves(depth) = nleaves
            if (nleaves > 0) call seed_leaves(self, anchor, depth, ncur, probe, context, error)
            if (allocated(error)) return
            exit
         end if

         call swap_frontiers(self)
         ncur = nnext
      end do

      self%rho_max_final = rho_cap

      if (self%debug) then
         !$omp critical (drop_octree_report)
         call self%report(anchor)
         !$omp end critical (drop_octree_report)
      end if
   end subroutine octree_run

   !> Exhausted box budget
   !>
   !> @param[out] error Budget refusal
   subroutine budget_error(error)
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, &
                       "Octree branch search: box budget exhausted before the "// &
                       "search ball was certified. Either the level set supplies "// &
                       "no exclusion radius (pruning nothing), or seed_size is too "// &
                       "small for the admissible radius.")
   end subroutine budget_error

   !> Exchange the current and next frontier storage
   !>
   !> @param[inout] self Seeder
   subroutine swap_frontiers(self)
      type(octree_seeder_type), intent(inout) :: self

      integer, allocatable :: lat(:, :)
      real(wp), allocatable :: key(:)

      call move_alloc(self%cur_lat, lat)
      call move_alloc(self%cur_key, key)
      call move_alloc(self%nxt_lat, self%cur_lat)
      call move_alloc(self%nxt_key, self%cur_key)
      call move_alloc(lat, self%nxt_lat)
      call move_alloc(key, self%nxt_key)
   end subroutine swap_frontiers

   !> Grow a frontier to `cap` entries, keeping the first `keep`
   !>
   !> @param[inout] lat  Lattice positions (3, :)
   !> @param[inout] key  Distance lower bounds
   !> @param[in]    cap  New capacity
   !> @param[in]    keep Entries to keep
   subroutine grow_frontier(lat, key, cap, keep)
      integer, allocatable, intent(inout) :: lat(:, :)
      real(wp), allocatable, intent(inout) :: key(:)
      integer, intent(in) :: cap
      integer, intent(in) :: keep

      integer, allocatable :: tmp_lat(:, :)
      real(wp), allocatable :: tmp_key(:)

      allocate (tmp_lat(3, cap), tmp_key(cap))
      if (keep > 0) then
         tmp_lat(:, 1:keep) = lat(:, 1:keep)
         tmp_key(1:keep) = key(1:keep)
      end if
      call move_alloc(tmp_lat, lat)
      call move_alloc(tmp_key, key)
   end subroutine grow_frontier

   !> Grow the per-box probe storage to at least `cap` boxes, keeping the first `keep`
   !>
   !> @param[inout] self Seeder
   !> @param[in]    cap  Boxes to hold
   !> @param[in]    keep Leading boxes to keep
   subroutine reserve_boxes(self, cap, keep)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: cap
      integer, intent(in) :: keep

      real(wp), allocatable :: c(:, :), v(:), r(:)
      logical, allocatable :: ok(:), surv(:)
      integer :: n

      if (allocated(self%probe0)) then
         if (size(self%probe0) >= cap) return
         n = max(cap, 2*size(self%probe0))
      else
         n = cap
      end if
      allocate (c(3, n), v(n), r(n), ok(n), surv(n))
      if (keep > 0) then
         c(:, 1:keep) = self%centers(:, 1:keep)
         v(1:keep) = self%probe0(1:keep)
         r(1:keep) = self%probe_r(1:keep)
         ok(1:keep) = self%probe_ok(1:keep)
         surv(1:keep) = self%survivor(1:keep)
      end if
      call move_alloc(c, self%centers)
      call move_alloc(v, self%probe0)
      call move_alloc(r, self%probe_r)
      call move_alloc(ok, self%probe_ok)
      call move_alloc(surv, self%survivor)
   end subroutine reserve_boxes

   !> Clear the per-search results and tallies
   !>
   !> @param[inout] self      Seeder
   !> @param[in]    rho_start Starting radius of the search (Bohr)
   subroutine reset_tallies(self, rho_start)
      type(octree_seeder_type), intent(inout) :: self
      real(wp), intent(in) :: rho_start

      self%n_points = 0
      self%reach = 0.0_wp
      self%rho_max_final = 0.0_wp
      self%n_boxes_visited = 0
      self%n_regions_excluded = 0
      self%dbg_probed = 0
      self%dbg_free = 0
      self%dbg_excluded = 0
      self%dbg_split = 0
      self%dbg_leaves = 0
      self%dbg_cap = 0.0_wp
      self%dbg_neighbours = 0
      self%dbg_evidence = 0
      self%dbg_refined = 0
      self%dbg_children = 0
      self%dbg_face_seeds = 0
      self%dbg_child_seeds = 0
      self%dbg_leaf_depth = 0
      self%dbg_rho_start = rho_start
   end subroutine reset_tallies

   !* ================================================================================= *!
   !*                                Sign-change seeding                                *!
   !* ================================================================================= *!

   !> Seeds from the surviving leaves of the leaf depth
   !>
   !> - face neighbours not probed at this depth are probed in one batch, outside
   !>   the box budget
   !> - a segment between opposite-sign face neighbours is kept once: by the
   !>   leaf on its negative side unless that neighbour is a survivor too
   !> - surviving leaves without evidence go to [[refine_leaves]]
   !>
   !> @param[inout] self    Seeder; `survivor(1:ncur)` flags the surviving leaves
   !> @param[in]    anchor  Anchor
   !> @param[in]    depth   Leaf depth
   !> @param[in]    ncur    Leaves probed at the leaf depth
   !> @param[in]    probe   Batched probe
   !> @param[inout] context Probe context
   !> @param[out]   error   Probe failure
   subroutine seed_leaves(self, anchor, depth, ncur, probe, context, error)
      type(octree_seeder_type), intent(inout) :: self
      real(wp), intent(in) :: anchor(3)
      integer, intent(in) :: depth
      integer, intent(in) :: ncur
      procedure(octree_probe_batch_interface) :: probe
      class(*), intent(inout) :: context
      type(error_type), allocatable, intent(out) :: error

      integer :: k, d, j, ntot, nseg, nref, nleaves, nbr(3)
      logical :: evidence

      nleaves = self%dbg_leaves(depth)
      call table_reset(self, ncur + 6*nleaves)
      do k = 1, ncur
         call table_insert(self, depth, self%cur_lat(:, k), k)
      end do

      ! Face neighbours never probed at this depth
      ntot = ncur
      do k = 1, ncur
         if (.not. self%survivor(k)) cycle
         do d = 1, 6
            nbr = self%cur_lat(:, k) + face_step(:, d)
            if (table_find(self, depth, nbr) > 0) cycle
            ntot = ntot + 1
            if (ntot > size(self%probe0)) call reserve_boxes(self, ntot, ntot - 1)
            self%centers(:, ntot) = anchor + box_offset(self, depth, nbr)
            self%survivor(ntot) = .false.
            call table_insert(self, depth, nbr, ntot)
         end do
      end do
      if (ntot > ncur) then
         self%dbg_neighbours = ntot - ncur
         call probe(self%centers(:, ncur + 1:ntot), self%probe0(ncur + 1:ntot), &
                    self%probe_r(ncur + 1:ntot), self%probe_ok(ncur + 1:ntot), context, error)
         if (allocated(error)) return
      end if

      ! Sign-change segments, and the leaves without any
      nseg = 0
      nref = 0
      do k = 1, ncur
         if (.not. self%survivor(k)) cycle
         evidence = .false.
         do d = 1, 6
            j = table_find(self, depth, self%cur_lat(:, k) + face_step(:, d))
            if (j == 0) error stop "moist DROP octree: face neighbour missing from the leaf table"
            if (.not. (self%probe_ok(k) .and. self%probe_ok(j))) cycle
            if ((self%probe0(k) > 0.0_wp) .eqv. (self%probe0(j) > 0.0_wp)) cycle
            evidence = .true.
            if (d > 3 .and. self%survivor(j)) cycle
            nseg = nseg + 1
            call reserve_segments(self, nseg)
            self%seg_a(:, nseg) = self%centers(:, k)
            self%seg_b(:, nseg) = self%centers(:, j)
            self%seg_va(nseg) = self%probe0(k)
            self%seg_vb(nseg) = self%probe0(j)
         end do
         if (evidence) then
            self%dbg_evidence = self%dbg_evidence + 1
         else if (self%probe_ok(k)) then
            nref = nref + 1
            call reserve_parents(self, nref)
            self%ref_c(:, nref) = self%centers(:, k)
            self%ref_v(nref) = self%probe0(k)
         end if
      end do
      self%dbg_face_seeds = nseg
      self%dbg_refined = nref

      call refine_leaves(self, nref, box_half(self, depth), nseg, probe, context, error)
      if (allocated(error)) return
      self%dbg_child_seeds = nseg - self%dbg_face_seeds

      do k = 1, nseg
         call self%push(segment_crossing(self%seg_a(:, k), self%seg_b(:, k), self%seg_va(k), &
                                         self%seg_vb(k)))
      end do
   end subroutine seed_leaves

   !> Probe the children of leaves without evidence for features thinner than a leaf
   !>
   !> - one batch per level over all parents, outside the box budget; a child on
   !>   the other side of its parent's center adds the segment between the two centers
   !> - the reached children of a parent without such a child become the next
   !>   level's parents
   !>
   !> @param[inout] self    Seeder; parents in `ref_c`, `ref_v`
   !> @param[in]    nref    Parents of the first level
   !> @param[in]    half    Half-width of the parents (Bohr)
   !> @param[inout] nseg    Sign-change segments, appended to
   !> @param[in]    probe   Batched probe
   !> @param[inout] context Probe context
   !> @param[out]   error   Probe failure
   subroutine refine_leaves(self, nref, half, nseg, probe, context, error)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: nref
      real(wp), intent(in) :: half
      integer, intent(inout) :: nseg
      procedure(octree_probe_batch_interface) :: probe
      class(*), intent(inout) :: context
      type(error_type), allocatable, intent(out) :: error

      integer :: level, npar, nkid, nnext, i, c, m
      real(wp) :: h
      logical :: found

      npar = nref
      h = half
      do level = 1, self%refine_levels
         if (npar == 0) exit
         nkid = 8*npar
         self%dbg_children = self%dbg_children + nkid
         call reserve_children(self, nkid)
         do i = 1, npar
            do c = 0, 7
               m = 8*(i - 1) + c + 1
               self%kid_c(1, m) = self%ref_c(1, i) + 0.5_wp*h*real(2*iand(c, 1) - 1, wp)
               self%kid_c(2, m) = self%ref_c(2, i) + 0.5_wp*h*real(2*iand(ishft(c, -1), 1) - 1, wp)
               self%kid_c(3, m) = self%ref_c(3, i) + 0.5_wp*h*real(2*iand(ishft(c, -2), 1) - 1, wp)
            end do
         end do
         call probe(self%kid_c(:, 1:nkid), self%kid_v(1:nkid), self%kid_r(1:nkid), &
                    self%kid_ok(1:nkid), context, error)
         if (allocated(error)) return

         nnext = 0
         do i = 1, npar
            found = .false.
            do m = 8*(i - 1) + 1, 8*i
               if (.not. self%kid_ok(m)) cycle
               if ((self%kid_v(m) > 0.0_wp) .eqv. (self%ref_v(i) > 0.0_wp)) cycle
               found = .true.
               nseg = nseg + 1
               call reserve_segments(self, nseg)
               self%seg_a(:, nseg) = self%ref_c(:, i)
               self%seg_b(:, nseg) = self%kid_c(:, m)
               self%seg_va(nseg) = self%ref_v(i)
               self%seg_vb(nseg) = self%kid_v(m)
            end do
            if (found .or. level == self%refine_levels) cycle
            ! Compacted in place: the write position never passes the read position
            do m = 8*(i - 1) + 1, 8*i
               if (.not. self%kid_ok(m)) cycle
               nnext = nnext + 1
               self%kid_c(:, nnext) = self%kid_c(:, m)
               self%kid_v(nnext) = self%kid_v(m)
            end do
         end do
         if (nnext > 0) then
            call reserve_parents(self, nnext)
            do i = 1, nnext
               self%ref_c(:, i) = self%kid_c(:, i)
               self%ref_v(i) = self%kid_v(i)
            end do
         end if
         npar = nnext
         h = 0.5_wp*h
      end do
   end subroutine refine_leaves

   !> Crossing of a sign-change segment by linear interpolation of its end values
   !>
   !> @param[in] a  First end point
   !> @param[in] b  Second end point
   !> @param[in] va Value at `a`
   !> @param[in] vb Value at `b`, on the other side of the surface
   pure function segment_crossing(a, b, va, vb) result(x)
      real(wp), intent(in) :: a(3)
      real(wp), intent(in) :: b(3)
      real(wp), intent(in) :: va
      real(wp), intent(in) :: vb
      real(wp) :: x(3)

      x = a + (va/(va - vb))*(b - a)
   end function segment_crossing

   !> Grow the segment storage to at least `n`
   !>
   !> @param[inout] self Seeder
   !> @param[in]    n    Segments to hold
   subroutine reserve_segments(self, n)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: n

      real(wp), allocatable :: a(:, :), b(:, :), va(:), vb(:)
      integer :: cap, keep

      keep = 0
      if (allocated(self%seg_va)) then
         if (size(self%seg_va) >= n) return
         keep = size(self%seg_va)
      end if
      cap = max(64, 2*keep, n)
      allocate (a(3, cap), b(3, cap), va(cap), vb(cap))
      if (keep > 0) then
         a(:, 1:keep) = self%seg_a
         b(:, 1:keep) = self%seg_b
         va(1:keep) = self%seg_va
         vb(1:keep) = self%seg_vb
      end if
      call move_alloc(a, self%seg_a)
      call move_alloc(b, self%seg_b)
      call move_alloc(va, self%seg_va)
      call move_alloc(vb, self%seg_vb)
   end subroutine reserve_segments

   !> Grow the refinement parent storage to at least `n`
   !>
   !> @param[inout] self Seeder
   !> @param[in]    n    Parents to hold
   subroutine reserve_parents(self, n)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: n

      real(wp), allocatable :: c(:, :), v(:)
      integer :: cap, keep

      keep = 0
      if (allocated(self%ref_v)) then
         if (size(self%ref_v) >= n) return
         keep = size(self%ref_v)
      end if
      cap = max(64, 2*keep, n)
      allocate (c(3, cap), v(cap))
      if (keep > 0) then
         c(:, 1:keep) = self%ref_c
         v(1:keep) = self%ref_v
      end if
      call move_alloc(c, self%ref_c)
      call move_alloc(v, self%ref_v)
   end subroutine reserve_parents

   !> Grow the refinement child storage to at least `n`; contents are not kept
   !>
   !> @param[inout] self Seeder
   !> @param[in]    n    Children to hold
   subroutine reserve_children(self, n)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: n

      integer :: cap

      cap = n
      if (allocated(self%kid_v)) then
         if (size(self%kid_v) >= n) return
         cap = max(n, 2*size(self%kid_v))
         deallocate (self%kid_c, self%kid_v, self%kid_r, self%kid_ok)
      end if
      allocate (self%kid_c(3, cap), self%kid_v(cap), self%kid_r(cap), self%kid_ok(cap))
   end subroutine reserve_children

   !* ================================================================================= *!
   !*                                 Leaf lattice table                                *!
   !* ================================================================================= *!

   !> Empty the table for up to `n` boxes
   !>
   !> @param[inout] self Seeder
   !> @param[in]    n    Boxes the table will hold
   subroutine table_reset(self, n)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: n

      integer :: size_needed

      size_needed = 64
      do while (size_needed < 2*n)
         size_needed = 2*size_needed
      end do
      if (allocated(self%table_box)) then
         if (size(self%table_box) < size_needed) deallocate (self%table_key, self%table_box)
      end if
      if (.not. allocated(self%table_box)) allocate (self%table_key(size_needed), self%table_box(size_needed))
      self%table_size = size_needed
      self%table_box(1:size_needed) = 0
   end subroutine table_reset

   !> Table key of a lattice position, the face neighbours outside the root included
   !>
   !> @param[in] depth Octree depth
   !> @param[in] lat   Lattice position, each in `-1 .. 2^depth`
   pure function table_key_of(depth, lat) result(key)
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      integer(int64) :: key

      integer(int64) :: m

      m = 2_int64**depth + 2_int64
      key = int(lat(1) + 1, int64) + m*(int(lat(2) + 1, int64) + m*int(lat(3) + 1, int64))
   end function table_key_of

   !> First probe slot of a lattice position
   !>
   !> @param[in] self Seeder
   !> @param[in] lat  Lattice position
   pure function table_slot(self, lat) result(slot)
      type(octree_seeder_type), intent(in) :: self
      integer, intent(in) :: lat(3)
      integer :: slot

      integer(int64) :: h

      h = ieor(ieor(int(lat(1) + 1, int64)*hash_mult(1), int(lat(2) + 1, int64)*hash_mult(2)), &
               int(lat(3) + 1, int64)*hash_mult(3))
      slot = int(iand(h, int(self%table_size - 1, int64))) + 1
   end function table_slot

   !> Record box `box` at a lattice position
   !>
   !> @param[inout] self  Seeder
   !> @param[in]    depth Octree depth
   !> @param[in]    lat   Lattice position
   !> @param[in]    box   Box index
   subroutine table_insert(self, depth, lat, box)
      type(octree_seeder_type), intent(inout) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      integer, intent(in) :: box

      integer(int64) :: key
      integer :: slot

      key = table_key_of(depth, lat)
      slot = table_slot(self, lat)
      do while (self%table_box(slot) /= 0)
         slot = slot + 1
         if (slot > self%table_size) slot = 1
      end do
      self%table_key(slot) = key
      self%table_box(slot) = box
   end subroutine table_insert

   !> Box at a lattice position, 0 when none
   !>
   !> @param[in] self  Seeder
   !> @param[in] depth Octree depth
   !> @param[in] lat   Lattice position
   function table_find(self, depth, lat) result(box)
      type(octree_seeder_type), intent(in) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      integer :: box

      integer(int64) :: key
      integer :: slot

      key = table_key_of(depth, lat)
      slot = table_slot(self, lat)
      box = 0
      do while (self%table_box(slot) /= 0)
         if (self%table_key(slot) == key) then
            box = self%table_box(slot)
            return
         end if
         slot = slot + 1
         if (slot > self%table_size) slot = 1
      end do
   end function table_find

   !* ================================================================================= *!
   !*                                    Box geometry                                   *!
   !* ================================================================================= *!

   !> Half-width of a box at `depth`
   !>
   !> @param[in] self  Seeder (root half-width)
   !> @param[in] depth Octree depth (0 = root)
   pure function box_half(self, depth) result(half)
      type(octree_seeder_type), intent(in) :: self
      integer, intent(in) :: depth
      real(wp) :: half

      half = self%root_half/real(2**depth, wp)
   end function box_half

   !> Depth at which a root of half-width `rho_max` reaches `seed_size`
   !>
   !> @param[in] self    Seeder
   !> @param[in] rho_max Root half-width (Bohr)
   pure function required_depth(self, rho_max) result(depth)
      type(octree_seeder_type), intent(in) :: self
      real(wp), intent(in) :: rho_max
      integer :: depth

      depth = 0
      do while (2.0_wp*rho_max/real(2**depth, wp) > self%seed_size)
         depth = depth + 1
         if (depth >= max_addressable_depth) exit
      end do
   end function required_depth

   !> Box center offset from the anchor, `h (2 lat + 1 - 2^d)`, exact in integers
   !>
   !> @param[in] self  Seeder
   !> @param[in] depth Octree depth
   !> @param[in] lat   Lattice position
   pure function box_offset(self, depth, lat) result(d)
      type(octree_seeder_type), intent(in) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      real(wp) :: d(3)

      d = box_half(self, depth)*real(2*lat + 1 - 2**depth, wp)
   end function box_offset

   !> Smallest anchor distance over a box
   !>
   !> @param[in] self  Seeder
   !> @param[in] depth Octree depth
   !> @param[in] lat   Lattice position
   pure function box_rho_lower(self, depth, lat) result(rho)
      type(octree_seeder_type), intent(in) :: self
      integer, intent(in) :: depth
      integer, intent(in) :: lat(3)
      real(wp) :: rho

      rho = norm2(max(abs(box_offset(self, depth, lat)) - box_half(self, depth), 0.0_wp))
   end function box_rho_lower

   !* ================================================================================= *!
   !*                                  Level set probe                                  *!
   !* ================================================================================= *!

   !> Level set values and surface-free radii at a batch of box centers
   !>
   !> - a center beyond the level set reach has no value: it reports the
   !>   anchor's own value (no cap tightening), the LSF's
   !>   `far_exclusion_radius()`, the surface-free radius shared by every
   !>   unreached point (zero when the LSF certifies none), and `reached` false
   !> - a reached center keeps `exclusion_radius(lsf0)`; the far radius is
   !>   not a valid floor there
   !>
   !> @param[in]    points  Box centers (3, n)
   !> @param[out]   lsf0    Values
   !> @param[out]   radius  Surface-free radii
   !> @param[out]   reached Whether the level set reaches the center
   !> @param[inout] context [[octree_lsf_context_type]]
   !> @param[out]   error   LSF failure
   subroutine lsf_probe(points, lsf0, radius, reached, context, error)
      real(wp), intent(in), contiguous :: points(:, :)
      real(wp), intent(out), contiguous :: lsf0(:)
      real(wp), intent(out), contiguous :: radius(:)
      logical, intent(out), contiguous :: reached(:)
      class(*), intent(inout) :: context
      type(error_type), allocatable, intent(out) :: error

      integer :: k, n
      real(wp) :: far_radius

      n = size(points, 2)
      select type (context)
      type is (octree_lsf_context_type)
         call context%lsf%evaluate_batch(points, 0, lsf0, error=error)
         if (allocated(error)) return
         far_radius = context%lsf%far_exclusion_radius()
         do k = 1, n
            ! Unreached columns hold `+Inf`: replaced before any sign test reads them
            reached(k) = .not. context%lsf%unreached(k)
            if (.not. reached(k)) then
               lsf0(k) = context%lsf0_anchor
               radius(k) = far_radius
            else
               radius(k) = context%lsf%exclusion_radius(lsf0(k))
            end if
         end do
      class default
         error stop "moist DROP octree: level set probe given a foreign context"
      end select
   end subroutine lsf_probe

   !> Whole-box stationarity exclusion from the level set
   !>
   !> @param[in] center     Box midpoint
   !> @param[in] half_width Cube half-width
   !> @param[in] context    [[octree_lsf_context_type]]
   function lsf_exclude_region(center, half_width, context) result(excluded)
      real(wp), intent(in) :: center(3)
      real(wp), intent(in) :: half_width
      class(*), intent(in) :: context
      logical :: excluded

      excluded = .false.
      select type (context)
      type is (octree_lsf_context_type)
         excluded = context%lsf%stationary_box_excluded(context%anchor, center, half_width)
      class default
         error stop "moist DROP octree: region certificate given a foreign context"
      end select
   end function lsf_exclude_region

   !> Upper bound on the closest branch distance by a Newton-step sphere trace
   !>
   !> - marches along the normal by `|S| / ||grad S||`, at most
   !>   `trace_step_limit` per step; a sign flip or a point within
   !>   `trace_feasible_tol` of the surface bounds `rho_min`
   !> - a ray that runs out of steps, meets a vanishing gradient or leaves the
   !>   level set reach bounds nothing
   !>
   !> @param[inout] ctx       Probe context
   !> @param[in]    request   Anchor request
   !> @param[out]   rho_bound Upper bound on `rho_min` when `found`
   !> @param[out]   found     Whether the ray reached the surface
   !> @param[out]   error     LSF failure
   subroutine trace_rho_bound(ctx, request, rho_bound, found, error)
      type(octree_lsf_context_type), intent(inout) :: ctx
      type(seed_request_type), intent(in) :: request
      real(wp), intent(out) :: rho_bound
      logical, intent(out) :: found
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x(3, 1), val(1), grad(3, 1), gnorm, dir(3), step
      integer :: istep

      rho_bound = 0.0_wp
      found = .false.
      x(:, 1) = request%anchor
      call ctx%lsf%evaluate_batch(x, 1, val, lsf1_r=grad, error=error)
      if (allocated(error)) return
      if (ctx%lsf%unreached(1)) return

      do istep = 1, trace_max_steps
         if (abs(val(1)) <= trace_feasible_tol) then
            rho_bound = norm2(x(:, 1) - request%anchor)
            found = .true.
            return
         end if
         gnorm = norm2(grad(:, 1))
         if (gnorm < trace_grad_floor) return
         dir = grad(:, 1)/gnorm
         if (val(1) > 0.0_wp) dir = -dir
         step = min(abs(val(1))/gnorm, trace_step_limit)
         x(:, 1) = x(:, 1) + step*dir
         call ctx%lsf%evaluate_batch(x, 1, val, lsf1_r=grad, error=error)
         if (allocated(error)) return
         if (ctx%lsf%unreached(1)) return
         if (val(1)*request%lsf0_anchor < 0.0_wp) then
            rho_bound = norm2(x(:, 1) - request%anchor)
            found = .true.
            return
         end if
      end do

   end subroutine trace_rho_bound

   !* ================================================================================= *!
   !*                                    Debug report                                   *!
   !* ================================================================================= *!

   !> Print the last search's per-depth tallies and seeds
   !>
   !> - callers inside an OpenMP region should serialize it
   !>
   !> @param[in] self   Seeder
   !> @param[in] anchor Anchor of the last search
   subroutine octree_report(self, anchor)
      class(octree_seeder_type), intent(in) :: self
      real(wp), intent(in) :: anchor(3)

      integer :: depth, k

      write (self%unit, "(a)") ""
      write (self%unit, "(2x,a,i0,a,3(1x,f10.5))") "== Octree seeds, point ", self%point_id, ", anchor", anchor
      write (self%unit, "(4x,a,f13.5)") "rho_max (start)        ", self%dbg_rho_start
      write (self%unit, "(4x,a,f13.5)") "rho_max (certified)    ", self%rho_max_final
      write (self%unit, "(4x,a,i13)") "leaf depth             ", self%dbg_leaf_depth
      write (self%unit, "(4x,7a10)") "depth", "probed", "free", "excluded", "split", "leaves", "cap"
      do depth = 0, self%dbg_leaf_depth
         if (self%dbg_probed(depth) == 0) cycle
         write (self%unit, "(4x,6i10,f10.4)") depth, self%dbg_probed(depth), self%dbg_free(depth), &
            self%dbg_excluded(depth), self%dbg_split(depth), self%dbg_leaves(depth), self%dbg_cap(depth)
      end do
      write (self%unit, "(4x,a,i13)") "face neighbours probed ", self%dbg_neighbours
      write (self%unit, "(4x,a,i13)") "leaves with evidence   ", self%dbg_evidence
      write (self%unit, "(4x,a,i13)") "leaves without evidence", self%dbg_refined
      write (self%unit, "(4x,a,i13)") "children probed        ", self%dbg_children
      write (self%unit, "(4x,a,i13)") "tree boxes probed      ", self%n_boxes_visited
      write (self%unit, "(4x,a,i13)") "face-segment seeds     ", self%dbg_face_seeds
      write (self%unit, "(4x,a,i13)") "child-segment seeds    ", self%dbg_child_seeds
      write (self%unit, "(4x,a,i13)") "seeds                  ", self%n_points
      do k = 1, self%n_points
         write (self%unit, "(6x,i6,4f13.5)") k, self%points(:, k), norm2(self%points(:, k) - anchor)
      end do
   end subroutine octree_report

end module moist_cavity_drop_projection_seed_octree

