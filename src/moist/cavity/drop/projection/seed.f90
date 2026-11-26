!> Lattice seeder of the batched DROP projection
!>
!> - `init` from the level's seed-step settings
!> - `build` once per geometry, between the local and the search pass
!>   ([[seed_build_type]])
!> - `seed` per anchor ([[seed_request_type]]), seeder read only, shared by
!>   the threads; results in the thread's [[seeder_scratch_type]]
!> - seeds `points(:, 1:n_points)`, count known before the caller appends, so
!>   an anchor's run is never split
!> - seeds within `radius + reach` stand for the surface inside `radius`
!> - same seeds at any thread count: fixed level set slices, clones without
!>   history
!> - lattice frame ([[lattice_frame_type]]): cells, edge midpoints and bins
!>   live in `o + A^T (x - o)`, `o` the centroid and `A` the axes; positions
!>   handed to the level set, samples and seeds are lab positions; lab axes
!>   (`A = I`) skip the map, so their seeds are bitwise those of a seeder
!>   without axes
!>
!> Build ([[seeder_build]], `seed/`):
!>
!> 1. region: union of the seeded anchors' largest balls (`radius_max`)
!> 2. lattice on the atom centroid along `axes` (the cavity's grid model,
!>    [[grid_model_type:lattice_axes]]), spacing `h 2^(m - s)`;
!>    `h = seed_size`, coarsest spacing at most [[discovery_size]],
!>    `refine_levels` below `h`
!> 3. `discover.f90`: coarse nodes within a cell diagonal of the region
!> 4. `refine.f90`: cells split while they change sign or lack a surface-free
!>    certificate (corner `exclusion_radius` balls); below `h` only cells
!>    without a sign change, `refine_levels` times
!> 5. `samples.f90`: one sample per sign-change edge, retracted onto `S = 0`,
!>    with normal and second fundamental form; binned
!> - sign change only towards unreached corners: no sample, `n_unresolved`
!> - deepest cell without sign change: no sample
!>
!> Seed ([[seeder_seed]], `seed/query.f90`), no level set evaluation:
!>
!> - samples of the cells meeting the ball; bins cleared by the stationarity
!>   certificate skipped
!> - level's filters ([[step_candidate]]): all four at level 2, cone at
!>   level 3, none at levels 1 and 4
!> - osculating pre-step towards the anchor, at most one spacing
!> - `reach`: covering radius plus longest retraction and pre-step
module moist_cavity_drop_projection_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int8, int64
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_projection_parameters, only: drop_seed_parameters_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type
   use moist_cavity_drop_projection_seed_lattice, only: key_bits, key_bias, key_table_type, table_init, table_find, &
                                                       table_insert, lattice_key, segment_crossing, &
                                                       evaluate_slices, lattice_frame_type
   use moist_utils_mem, only: grow_array
   use moist_math_linalg_geometry, only: identity3, is_rotation
   use moist_math_linalg, only: setup_tangent_frame
   use moist_math_packed, only: unpack_sym_matrix
   use moist_math_sorter, only: sort_index
   use moist_utils_prettyprint, only: prettyprinter
   implicit none(type, external)
   private

   public :: projection_seeder_type, seeder_scratch_type, seed_request_type, seed_build_type

   !> Coarsest lattice spacing without a strong exclusion certificate (Bohr)
   real(wp), parameter :: discovery_size = 1.0_wp
   !> Coarsest lattice spacing with a strong exclusion certificate (Bohr)
   real(wp), parameter :: discovery_size_certified = 4.0_wp
   !> Maximum dense 3D mask entries (bytes) used during coarse volume discovery
   integer(int64), parameter :: max_mask_cells = 67108864_int64
   !> Points per level set call of the build
   integer, parameter :: build_chunk = 256
   !> Capped normal Newton steps of the retraction before the order-2 jet
   integer, parameter :: retract_steps = 2
   !> Grazing filter: a sample farther than its spacing whose offset from the
   !> anchor makes `|d . n| / |d|` smaller than 0.20 is dropped (squared threshold)
   real(wp), parameter :: grazing_cos2 = 0.04_wp
   !> Back-facing filter: a sample whose normal makes `n . n_sph` at most this
   !> with the anchor's outward sphere direction is dropped
   real(wp), parameter :: back_facing_cos = -0.30_wp
   !> Corner pairs of the 12 cell edges, corner `1+i+2j+4k`
   integer, parameter :: edge_corner(2, 12) = reshape([1, 2, 3, 4, 5, 6, 7, 8, &
                                                       1, 3, 2, 4, 5, 7, 6, 8, &
                                                       1, 5, 2, 6, 3, 7, 4, 8], [2, 12])

   !> Per-anchor seeding request: the anchor and its ball
   !>
   !> - what depends on the projection level alone (the seed filters) is the
   !>   seeder's setting, from `drop_seed_parameters_type`
   type :: seed_request_type
      !> Anchor position
      real(wp) :: anchor(3) = 0.0_wp
      !> Owning atom, 0 when unknown
      integer :: owner = 0
      !> Owning atom center
      real(wp) :: owner_center(3) = 0.0_wp
      !> Radius of the ball to search (Bohr), positive; the whole ball is searched,
      !> surface met on the way does not shrink it
      real(wp) :: pdist_start = -1.0_wp
      !> Inner ball radius already searched (Bohr); non-negative requests only the
      !> outer shell `(pdist_inner, pdist_start]`; negative searches the whole ball
      real(wp) :: pdist_inner = -1.0_wp
      !> Seeds ordered by distance to the anchor, nearest first
      logical :: nearest_first = .false.
      !> Whether the anchor already has an accepted local point; the tangent-cone
      !> and back-facing filters apply only then
      logical :: has_accepted = .false.
   end type seed_request_type

   !> What the seeder's build sees: the geometry and the search ahead
   type :: seed_build_type
      !> Atom centers (3, nat), finite
      real(wp), allocatable :: centers(:, :)
      !> Atom radii (nat, Bohr)
      real(wp), allocatable :: radii(:)
      !> Anchors the search pass will seed (3, n), in anchor order, finite
      real(wp), allocatable :: anchors(:, :)
      !> Largest radius each such anchor's search may ask for (n, Bohr), finite:
      !> `radius` for a search of one ball, `reach` for an extended search;
      !> unallocated, `reach` for every anchor
      real(wp), allocatable :: radius_max(:)
      !> Largest search radius (Bohr), positive and finite
      real(wp) :: reach = 0.0_wp
      !> Threads the build may use; the result must not depend on it
      integer :: nthreads = 1
      !> Lattice axes as columns, a proper rotation; the lab axes by default
      real(wp) :: axes(3, 3) = identity3
   end type seed_build_type

   !> Lattice nodes of a build: finest-unit coordinates, values, certificates
   type :: node_set_type
      !> Nodes
      integer :: n = 0
      !> Nodes evaluated so far; `nval+1:n` wait for their values
      integer :: nval = 0
      !> Coordinates in finest lattice units (3, n)
      integer, allocatable :: ijk(:, :)
      !> Level set values
      real(wp), allocatable :: val(:)
      !> Surface-free radius around each node (Bohr)
      real(wp), allocatable :: free(:)
      !> Whether the level set reaches the node
      logical, allocatable :: reached(:)
      !> Key to node
      type(key_table_type) :: table
   end type node_set_type

   !> Working state of one build
   type :: build_type
      !> Level set clones of the build, one per thread, filled on first use
      type(drop_worker_slots_type) :: slots
      !> Lattice origin and axes
      type(lattice_frame_type) :: frame
      !> Finest lattice spacing, `h 2^-refine_levels` (Bohr)
      real(wp) :: unit = 1.0_wp
      !> Coarse spacing `h_0` (Bohr)
      real(wp) :: h0 = 1.0_wp
      !> Level of the sample size `h`
      integer :: m = 0
      !> Deepest level
      integer :: depth = 0
      !> Lattice nodes
      type(node_set_type) :: nodes
      !> Crossings
      integer :: ncross = 0
      !> Crossing positions (3, :)
      real(wp), allocatable :: cross(:, :)
      !> Edge length of each crossing (Bohr)
      real(wp), allocatable :: cross_h(:)
      !> Edge midpoint of each crossing (3, :)
      real(wp), allocatable :: cross_mid(:, :)
      !> Edge axis of each crossing
      integer, allocatable :: cross_axis(:)
   end type build_type

   !> One thread's scratch: the seeds and candidates of the last request
   type :: seeder_scratch_type
      !> Seeds of the last request (3, capacity)
      real(wp), allocatable :: points(:, :)
      !> Seeds produced by the last request
      integer :: n_points = 0
      !> Distance beyond the searched ball at which the last request's seeds still
      !> stand for surface inside it (Bohr); 0 without a sample
      real(wp) :: reach = 0.0_wp
      !> Whether the last request returned only seeds outside `request%pdist_inner`,
      !> so no seed of the inner ball is repeated
      logical :: shell_filtered = .false.
      !> Candidate samples
      integer :: ncand = 0
      !> Sample index of each candidate
      integer, allocatable :: cand(:)
      !> Sort permutation of the candidates
      integer, allocatable :: cand_perm(:)
      !> Starting point of each kept candidate (3, :)
      real(wp), allocatable :: cand_x0(:, :)
      !> Squared distance of each kept start to the anchor (Bohr^2)
      real(wp), allocatable :: cand_dist(:)
   end type seeder_scratch_type

   !> Build-first lattice seeder; built once, read by every thread
   type :: projection_seeder_type

      !* ---------------------------------- Settings ---------------------------------- *!

      !> Sample spacing `h`, the finest regular lattice (Bohr)
      real(wp) :: sample_size = 1.0_wp
      !> Lattice levels below `h` for a cell without sign change
      integer :: refine_levels = 0
      !> Seed filter ([[step_candidate]]): tangent cone of an anchor with an
      !> accepted local point; with all four filters off every sample is kept
      logical :: filter_cone = .false.
      !> Seed filter: patches the anchor sees edge-on
      logical :: filter_grazing = .false.
      !> Seed filter: back-facing samples of an anchor with an accepted local point
      logical :: filter_back_facing = .false.
      !> Seed filter: samples beyond both focal centers of their patch
      logical :: filter_focal = .false.
      !> Tangent-cone scale in units of `h (1 + r k_max)`
      real(wp) :: cone_scale = 1.5_wp
      !> Whether `init` ran
      logical :: initialized = .false.

      !* ------------------------------- Built lattice -------------------------------- *!

      !> Whether `build` ran for the current geometry
      logical :: built = .false.
      !> Lattice origin and axes
      type(lattice_frame_type) :: frame
      !> Samples
      integer :: nsample = 0
      !> Retracted positions (3, nsample)
      real(wp), allocatable :: y(:, :)
      !> Unit normals (3, nsample), zero where the retraction failed
      real(wp), allocatable :: normal(:, :)
      !> Second fundamental form in the tangent frame of the normal: 11, 12, 22 (3, nsample)
      real(wp), allocatable :: form(:, :)
      !> Largest principal curvature magnitude `max(|k1|, |k2|)` per sample (1/Bohr)
      real(wp), allocatable :: kmax(:)
      !> Lattice spacing of each sample's edge (Bohr)
      real(wp), allocatable :: spacing(:)
      !> Midpoint of each sample's edge (3, nsample); the four cells around the
      !> edge decide whether a ball sees the sample
      real(wp), allocatable :: edge_mid(:, :)
      !> Axis of each sample's edge
      integer, allocatable :: edge_axis(:)
      !> Covering radius of the samples (Bohr)
      real(wp) :: cover = 0.0_wp
      !> Largest sample spacing (Bohr)
      real(wp) :: spacing_max = 0.0_wp
      !> Bin origin (Bohr)
      real(wp) :: bin_origin(3) = 0.0_wp
      !> Bin edge (Bohr)
      real(wp) :: bin_size = 1.0_wp
      !> Key table of the bins
      type(key_table_type) :: bins
      !> First sample of each bin (bin order) and its end
      integer, allocatable :: bin_first(:)

      !* ---------------------------------- Tallies ---------------------------------- *!

      !> Cells of the last build that may hold surface but gave no seed, reported
      integer :: n_unresolved = 0
      !> Order-0 node probes of the build
      integer :: n_probes = 0
      !> Order-1 evaluations of the build
      integer :: n_order1 = 0
      !> Order-2 evaluations of the build
      integer :: n_order2 = 0
      !> Bytes held by the built lattice
      integer(int64) :: lattice_bytes = 0_int64
   contains
      !> Validate the settings
      procedure :: init => seeder_init
      !> Build the seed lattice for one geometry
      procedure :: build => seeder_build
      !> Seeds of one anchor
      procedure :: seed => seeder_seed
      !> Print the build tallies
      procedure :: report => seeder_report
   end type projection_seeder_type

   interface

      !> [seed/build.f90] Build the seed lattice for one geometry
      module subroutine seeder_build(self, lsf, geometry, error)
         implicit none(type, external)
         class(projection_seeder_type), intent(inout) :: self
         class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
         type(seed_build_type), intent(in) :: geometry
         type(error_type), allocatable, intent(out) :: error
      end subroutine seeder_build

      !> [seed/build.f90] Drop the built lattice
      module subroutine clear_lattice(self)
         implicit none(type, external)
         type(projection_seeder_type), intent(inout) :: self
      end subroutine clear_lattice

      !> [seed/nodes.f90] Lattice step of level `s` in finest units
      pure module function level_step(work, s) result(step)
         implicit none(type, external)
         type(build_type), intent(in) :: work
         integer, intent(in) :: s
         integer :: step
      end function level_step

      !> [seed/nodes.f90] Lattice-frame position of a lattice node
      pure module function node_position(work, ijk) result(x)
         implicit none(type, external)
         type(build_type), intent(in) :: work
         integer, intent(in) :: ijk(3)
         real(wp) :: x(3)
      end function node_position

      !> [seed/nodes.f90] Values and certificates of the nodes waiting for them
      module subroutine evaluate_nodes(self, lsf, work, error)
         implicit none(type, external)
         type(projection_seeder_type), intent(inout) :: self
         class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
         type(build_type), intent(inout) :: work
         type(error_type), allocatable, intent(out) :: error
      end subroutine evaluate_nodes

      !> [seed/nodes.f90] Add a lattice node unless known; its value waits for [[evaluate_nodes]]
      module subroutine add_node(nodes, ijk)
         implicit none(type, external)
         type(node_set_type), intent(inout) :: nodes
         integer, intent(in) :: ijk(3)
      end subroutine add_node

      !> [seed/nodes.f90] Lattice coordinates out of the key range
      module subroutine range_error(error)
         implicit none(type, external)
         type(error_type), allocatable, intent(out) :: error
      end subroutine range_error

      !> [seed/discover.f90] Probe every coarse node within a cell diagonal of the searched anchor balls
      module subroutine discover_volume(self, lsf, geometry, work, error)
         implicit none(type, external)
         type(projection_seeder_type), intent(inout) :: self
         class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
         type(seed_build_type), intent(in) :: geometry
         type(build_type), intent(inout) :: work
         type(error_type), allocatable, intent(out) :: error
      end subroutine discover_volume

      !> [seed/refine.f90] Split cells level by level and collect the crossings of the unsplit ones
      module subroutine refine_cells(self, lsf, work, error)
         implicit none(type, external)
         type(projection_seeder_type), intent(inout) :: self
         class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
         type(build_type), intent(inout) :: work
         type(error_type), allocatable, intent(out) :: error
      end subroutine refine_cells

      !> [seed/samples.f90] Retract the crossings onto `S = 0` and store normal and second fundamental form
      module subroutine retract_samples(self, lsf, work, error)
         implicit none(type, external)
         type(projection_seeder_type), intent(inout) :: self
         class(moist_cavity_drop_lsf_type), intent(in), target :: lsf
         type(build_type), intent(inout) :: work
         type(error_type), allocatable, intent(out) :: error
      end subroutine retract_samples

      !> [seed/samples.f90] Sort the samples into cubic bins of edge `edge` by their edge midpoint
      module subroutine build_bins(self, edge)
         implicit none(type, external)
         type(projection_seeder_type), intent(inout) :: self
         real(wp), intent(in) :: edge
      end subroutine build_bins

      !> [seed/samples.f90] Bin coordinates of a point
      pure module function bin_coords(self, x) result(j)
         implicit none(type, external)
         type(projection_seeder_type), intent(in) :: self
         real(wp), intent(in) :: x(3)
         integer :: j(3)
      end function bin_coords

      !> [seed/query.f90] Seeds of one anchor from the built lattice, no level set evaluation
      module subroutine seeder_seed(self, scratch, lsf, request, error)
         implicit none(type, external)
         class(projection_seeder_type), intent(in) :: self
         type(seeder_scratch_type), intent(inout) :: scratch
         class(moist_cavity_drop_lsf_type), intent(in) :: lsf
         type(seed_request_type), intent(in) :: request
         type(error_type), allocatable, intent(out) :: error
      end subroutine seeder_seed

   end interface

contains

   !* ================================================================================= *!
   !*                                     Lifecycle                                     *!
   !* ================================================================================= *!

   !> Take the seed-step settings of the level and validate them
   !>
   !> - a previous build is dropped, whether the settings pass or not
   !>
   !> @param[inout] self  Seeder
   !> @param[in]    seed  Seed-step settings: sample spacing positive and finite,
   !>                     refinement not negative, cone scale positive and finite
   !> @param[out]   error Invalid settings
   subroutine seeder_init(self, seed, error)
      class(projection_seeder_type), intent(inout) :: self
      type(drop_seed_parameters_type), intent(in) :: seed
      type(error_type), allocatable, intent(out) :: error

      self%sample_size = seed%seed_size
      self%refine_levels = seed%refine_levels
      self%filter_cone = seed%filter_cone
      self%filter_grazing = seed%filter_grazing
      self%filter_back_facing = seed%filter_back_facing
      self%filter_focal = seed%filter_focal
      self%cone_scale = seed%cone_scale
      self%initialized = .false.
      call clear_lattice(self)
      if (.not. (self%sample_size > 0.0_wp)) then
         call fatal_error(error, "Lattice seeder: sample size must be positive")
         return
      end if
      if (.not. ieee_is_finite(self%sample_size)) then
         call fatal_error(error, "Lattice seeder: sample size must be finite")
         return
      end if
      if (self%refine_levels < 0) then
         call fatal_error(error, "Lattice seeder: refine_levels must not be negative")
         return
      end if
      ! A negative scale acts as its magnitude, a non-finite one switches the cone off
      if (.not. (self%cone_scale > 0.0_wp)) then
         call fatal_error(error, "Lattice seeder: cone scale must be positive")
         return
      end if
      if (.not. ieee_is_finite(self%cone_scale)) then
         call fatal_error(error, "Lattice seeder: cone scale must be finite")
         return
      end if
      self%initialized = .true.
   end subroutine seeder_init

   !> Print the build tallies
   !>
   !> @param[in]    self Seeder
   !> @param[inout] pp   Printer, its gate decides whether anything prints
   subroutine seeder_report(self, pp)
      class(projection_seeder_type), intent(in) :: self
      type(prettyprinter), intent(inout) :: pp

      call pp%push("Lattice seeder:")
      call pp%kv("Sample spacing", self%sample_size, "Bohr")
      call pp%kv("Node probes (order 0)", self%n_probes)
      call pp%kv("Samples", self%nsample)
      call pp%kv("Unresolved cells", self%n_unresolved)
      call pp%kv("Retractions to order 1", self%n_order1)
      call pp%kv("Retractions to order 2", self%n_order2)
      call pp%kv("Covering radius", self%cover, "Bohr")
      call pp%kv("Lattice memory", 1.0e-6_wp*real(self%lattice_bytes, wp), "MB")
      call pp%pop()
   end subroutine seeder_report

end module moist_cavity_drop_projection_seed
