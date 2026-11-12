!> Abstract level set function (LSF) base type for the DROP scheme
!>
!> - DROP's discretization is independent of *which* LSF defines the surface
!> - every concrete LSF extends this base
!>
!> A concrete LSF supplies:
!>   - its own constructor with the per-LSF parameters, e.g.
!>     [[moist_cavity_drop_lsf_svdw_type]]%new
!>   - a `candidate_space` declaration in that constructor (see below)
!>   - implementations of the deferred procedures
!>   - optionally, overrides of the higher-order and direction-contracted derivatives the base
!>     declares with aborting defaults (`f4_*`, `f2_rArB`, `f3_r_rArB`, `f4_rr_rArB`,
!>     `f*_radrad`, `f*_rA_rad`, `normalized_f01_rA`, `tangent_*`, `hvp_*`, `vjp_f1_rA`,
!>     `vjp_f1_rad`)
!>
!> Coverage of that optional family:
!>   - SvdW and CFC: all of it
!>   - isodensity: only the identically zero `vjp_f1_rA`; clears `radius_dependent` instead of
!>     implementing the radius family
!>   - fourth derivatives and radius queries (including `vjp_f1_rad`) abort there rather than
!>     return a zero the caller cannot tell from a real answer
!>
!> Packed spatial groups
!> ---------------------
!>   - symmetric spatial groups of the per-point accessors packed as [[moist_math_packed_sym]]:
!>     `rr` (6) xx, xy, yy, xz, yz, zz; `rrr` (10); `rrrr` (15)
!>   - other indices keep their dense extents
!>   - `vjp_f1_rA` / `vjp_f1_rad` seed `w2(6)` packed symmetric, off-diagonals not doubled
!>
!> Batched families
!> ----------------
!> Point-last, packed symmetric groups, one shared candidate list per tile, per-column `status`
!> (0 ok, 1 unreached, 2 failed); looping default over the single-point path, overridable by
!> a kernel:
!>
!>   - `evaluate_batch`: the jet, orders 0..4
!>   - `vjp_f1_rA_batch`: seed-contracted nuclear derivative of the jet (batched `vjp_f1_rA`),
!>     output `(3, ncand, np)`, exact zeros out of reach
!>
!> Atomic partition
!> ----------------
!> `atom_partition` splits the level set at a point into atomic shares summing to one, from the
!> LSF's own atomic terms; `has_atom_partition` reports whether the override exists
!>
!>   - SvdW and CFC: each one-, two- and three-body term shared equally among its atoms
!>   - internal isodensity: Mulliken split of the density
!>   - callback isodensity: no atomic terms, aborting default
!>
!> Geometric bounds
!> ----------------
!> Two analytic "how far" queries that let a caller skip work without knowing the concrete LSF;
!> opposite safe directions, so neither is expressible through the other:
!>
!>   - [[exclusion_radius]]`(lsf0)`: minimum distance from the point to any *surface*, given
!>     the value there
!>     - consumer: octree surface coverage, "can I discard this region?"
!>     - safe answer **zero** (prune nothing); an over-estimate silently loses branches
!>     - companion [[far_exclusion_radius]]`()`: same promise for every point no atom reaches
!>       (batch `status = 1`), a model constant, zero by default
!>   - [[screening_offset]]`(radius)`: maximum distance from the point of an *atom* still
!>     contributing above `screening_threshold` (SvdW with fade: the band end, threshold unread)
!>     - consumer: screening, "can I skip this atom?"
!>     - safe answer **huge** (screen nothing); an under-estimate silently drops contributions
!>
!> Screening contract
!> ------------------
!>    `screening_offset(radius)`: radial offset outward from the atom surface at which the
!>    atom's contribution equals `screening_threshold`
!>
!>   - exact analytic inverse, never a bound or a search
!>   - single source of the reach `R + offset(R)` for both consumers, rebuilt together in
!>     [[lsf_base_refresh_bounds]] so they cannot drift apart:
!>     - per-point gate ([[lsf_base_screen_candidates]]): squared-distance compare against
!>       `(R + offset(R))**2`, cached per `update` beside the center in `cand_screen`; no `exp`,
!>       no `sqrt` in the reject path
!>     - candidate cell list ([[lsf_base_candidates_in_box]],
!>       [[lsf_base_candidates_in_sphere]]): reach without margin, a superset of the gate's
!>       survivors; the gate stays the exact filter
!>
!> Candidate index space
!> ---------------------
!> Space of the candidate ids reaching `prepare_subset` (from the cell list or a caller's
!> `candidate_indices`), declared per LSF, never inferred:
!>
!>   - [[lsf_candidate_space_sorted]]: the base's spatial sort; near-sequential geometry gathers
!>     in the screen loop, ~15% of a projection update
!>   - [[lsf_candidate_space_user]]: the caller's own ids, for per-atom data indexed by them
!>   - set in the constructor; `update` builds the candidate mirror (`cand_screen`,
!>     `cand_radii`, `cand_to_user`) and the cell list in that space, and aborts on an
!>     undeclared one
!>
!> Thread-local clones
!> -------------------
!>   - the cavity holds one `class(moist_cavity_drop_lsf_type), allocatable` model and clones
!>     it per thread
!>   - one clone per [[lsf_thread_slot]], each allocated independently
module moist_cavity_drop_lsf_base
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type
   use moist_math_cell_list, only: moist_cell_list_type, heapsort
   use moist_math_sorter_counting_sort, only: counting_argsort
   implicit none(type, external)
   private

   public :: moist_cavity_drop_lsf_type
   public :: lsf_thread_slot
   public :: lsf_base_update
   public :: lsf_candidate_space_undeclared
   public :: lsf_candidate_space_user
   public :: lsf_candidate_space_sorted
   public :: lsf_cand_bound_sq
   public :: lsf_batch_require_outputs
   public :: lsf_vjp_batch_shape_message
   public :: lsf_sorted_candidates
   public :: lsf_screen_records

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Row range of the atom center inside a `cand_screen` record
   integer, parameter :: lsf_cand_center(2) = [1, 3]
   !> Row of the squared screening reach inside a `cand_screen` record
   integer, parameter :: lsf_cand_bound_sq = 4
   !> Number of rows in a `cand_screen` record
   integer, parameter :: lsf_cand_screen_rows = 4

   !> `candidate_space` sentinel: the concrete LSF never declared one
   !>
   !> - reaching `update` in this state is a programming error and aborts
   integer, parameter :: lsf_candidate_space_undeclared = 0
   !> `candidate_space`: user-space (caller) atom ids, handed to `prepare_subset`
   !> untouched
   !>
   !> - required by LSFs indexing per-atom data (basis shells, ...) with the
   !>   caller's numbering
   integer, parameter :: lsf_candidate_space_user = 1
   !> `candidate_space`: candidate ids relabelled into the base's spatially-sorted
   !> ordering
   !>
   !> - lets `prepare_subset` index `cand_screen` / `cand_radii` directly, with
   !>   near-sequential access
   integer, parameter :: lsf_candidate_space_sorted = 2

   !> Default candidate count below which the candidate cell list is a single cell
   !>
   !> - initial value of `cell_full_scan_below`, changed through
   !>   [[lsf_base_set_full_scan_below]]
   !> - both modes return the same candidates, so the switch only affects speed
   integer, parameter :: lsf_cell_list_full_scan_below = 200

   !> Abstract LSF base
   type, abstract :: moist_cavity_drop_lsf_type
      !> Short name of the level-set model, e.g. `SvdW`; set by the constructor
      character(len=:), allocatable :: name
      !> Number of atomic centers
      integer :: ncenters = 0
      !> Molecular structure data
      type(structure_type) :: mol
      !> Lower corner of the atomic center enclosure, refreshed with geometry
      real(wp) :: center_lower(3) = 0.0_wp
      !> Upper corner of the atomic center enclosure
      real(wp) :: center_upper(3) = 0.0_wp
      !> Whether the enclosure contains finite, nonempty geometry
      logical :: center_bounds_valid = .false.
      !> Atomic radii (per-atom, user-space)
      real(wp), allocatable :: radii(:)
      !> Screening threshold below which the LSF contribution is treated as zero
      !>
      !> - owned by the cavity, coupled to the projection tolerance
      !> - read by the concrete LSF to size its screening caches on every `update`
      !> - direct users (tests calling the concrete without a cavity) may set it
      !>   before calling `update`
      real(wp) :: screening_threshold = 0.0_wp
      !> Highest spatial derivative order actually present in the per-point cache after
      !> the latest `prepare`/`prepare_subset`
      integer :: prepared_deriv = -1
      !> Point of the latest `prepare`/`prepare_subset` is unreached (escape)
      !>
      !> - set by every concrete `prepare`; SvdW and CFC only, never for the
      !>   isodensity LSFs or SSES
      !> - every accessor aborts on it ([[lsf_base_require_reached]]), except
      !>   `active_count` and `active_atom`
      logical :: prepared_unreached = .false.

      !> Index space the concrete LSF wants its `prepare_subset` candidate ids in
      !> Set by the concrete's constructor to [[lsf_candidate_space_user]] or
      !> [[lsf_candidate_space_sorted]]; `update` aborts while undeclared
      integer :: candidate_space = lsf_candidate_space_undeclared

      !> LSF dependency on atomic radii
      logical :: radius_dependent = .true.

      !* ------------------------- Spatial sort machinery -------------------------- *!

      !> Permutation mapping a user-space atom index i to its spatially-sorted
      !> position, rebuilt by `update` / `set_centers` [ncenters]
      integer, allocatable :: orig_to_sorted(:)
      !> Inverse of `orig_to_sorted`: the user-space atom index stored at
      !> spatially-sorted position j [ncenters]
      integer, allocatable :: sorted_to_orig(:)

      !* ------------------- Candidate-space geometry and screening ---------------- *!
      !> The arrays below are indexed by *candidate id in the declared
      !> `candidate_space`*, so a screen loop consumes `prepare_subset`'s ids with
      !> no per-candidate translation:
      !>
      !> - `sorted` LSF -- the spatially-sorted mirror of the geometry
      !> - `user` LSF -- the user-space arrays themselves

      !> Everything the per-point reject test reads, interleaved per candidate
      !> [4, ncenters]:
      !>   `cand_screen(1:3, c)` -- the atom center
      !>   `cand_screen(4, c)`   -- the squared screening reach `(R + offset(R))**2`,
      !>                            saturated at `huge` rather than overflowed when
      !>                            screening is disabled
      !> One 32-byte record, so rejecting a candidate touches a single cache line
      !> instead of striding through separate centers and bounds arrays
      !>
      !> - index with [[lsf_cand_center]] / [[lsf_cand_bound_sq]], never bare literals
      real(wp), allocatable :: cand_screen(:, :)
      !> Atom radii in candidate space, read only for candidates that pass the gate
      !> and hence kept out of `cand_screen` [ncenters]
      real(wp), allocatable :: cand_radii(:)
      !> User-space atom id of candidate c [ncenters]
      integer, allocatable :: cand_to_user(:)
      !> Candidate ids covering every atom, ordered so the resulting active list
      !> comes out in ascending user-space order, fed to the screen loop by
      !> `prepare` (the unscreened full scan) [ncenters]
      integer, allocatable :: full_scan_cand(:)
      !> Candidate-space cell list over `cand_screen` centers with the gate's exact
      !> reach `R + offset(R)`, no margin, rebuilt with the screening bounds
      type(moist_cell_list_type) :: cand_cells
      !> Candidate count below which `cand_cells` is a single cell
      !>
      !> - set through [[lsf_base_set_full_scan_below]], which the cavity calls
      !>   with `screening.cell_grid_full_scan_below`
      integer :: cell_full_scan_below = lsf_cell_list_full_scan_below
   contains

      !* --------------------------- Geometry and screening --------------------------- *!

      !> Bind molecular geometry and rebuild the screening caches
      !>
      !> - concrete LSFs may override to refresh additional caches
      !> - an override *must* call `lsf_base_update(self, mol, radii)` first
      procedure :: update => lsf_base_update
      !> Move the atom centers without touching radii, parameters or derivative
      !> storage (used by nuclear finite-difference drivers)
      procedure, non_overridable :: set_centers => lsf_base_set_centers
      !> Reject candidates whose contribution is below `screening_threshold`
      procedure, non_overridable :: screen_candidates => lsf_base_screen_candidates
      !> Candidate ids whose screening reach sphere meets an axis-aligned box
      procedure, non_overridable :: candidates_in_box => lsf_base_candidates_in_box
      !> Candidate ids whose screening reach sphere meets a ball
      procedure, non_overridable :: candidates_in_sphere => lsf_base_candidates_in_sphere
      !> Change `screening_threshold` in place, rewriting only the cached reach
      procedure, non_overridable :: set_screening_threshold => lsf_base_set_screening_threshold
      !> Change the candidate count below which the candidate cell list is one cell
      procedure, non_overridable :: set_full_scan_below => lsf_base_set_full_scan_below
      !> Whether the cached reach still equals the live `screening_offset`
      procedure, non_overridable :: reach_current => lsf_base_reach_current
      !> Exact radial offset from the atom surface at which this LSF's
      !> contribution equals `screening_threshold`
      procedure(lsf_screening_offset_iface), deferred :: screening_offset
      !> Geometric bound query 1 of 2: how far a *surface* can be from a point
      !> with a known LSF value (see the "Two geometric-bound queries" note above)
      procedure :: exclusion_radius => lsf_base_exclusion_radius
      !> Surface-free radius around every point no atom reaches, a model constant
      procedure :: far_exclusion_radius => lsf_base_far_exclusion_radius
      !> Whether a box provably contains no stationary distance point
      procedure :: stationary_box_excluded => lsf_base_stationary_box_excluded

      !* -------------------------- Density and build phase --------------------------- *!

      !> Chain-rule factor `dS/drho` of a density-backed LSF, and whether this
      !> LSF has a density at all (see [[lsf_base_density_adjoint_factor]])
      procedure :: density_adjoint_factor => lsf_base_density_adjoint_factor
      !> Reject projected surface points this LSF cannot represent (see
      !> [[lsf_base_check_surface]])
      procedure :: check_surface => lsf_base_check_surface
      !> Declare the host requests of the build phase, at the cavity scope (none by default)
      procedure :: declare_coupling => lsf_base_declare_coupling
      !> Whether the level set requires build-phase answers before `update` (no by default)
      procedure :: requires_build_inputs => lsf_base_requires_build_inputs
      !> Read the build-phase answers before `update` (nothing by default)
      procedure :: read_build => lsf_base_read_build

      !* ------------------------------ Point evaluation ------------------------------ *!

      !> Cache per-point screening / state
      procedure(lsf_prepare_iface), deferred :: prepare
      !> Cache per-point state with a caller-provided candidate list
      procedure(lsf_prepare_subset_iface), deferred :: prepare_subset
      !> Configure the highest spatial derivative order required
      procedure(lsf_set_max_deriv_iface), deferred :: set_max_deriv
      !> Abort when an accessor is asked for an order `prepare` did not compute
      procedure :: require_deriv => lsf_base_require_deriv
      !> Abort when an accessor reads the cache of an unreached prepare
      procedure :: require_reached => lsf_base_require_reached
      !> Whether the prepared point is unreached, no atom in reach
      procedure :: unreached => lsf_base_unreached
      !> Report an unreached prepare through `status`, or as `error` without it
      procedure :: report_unreached => lsf_base_report_unreached
      !> Batched jet at many points, orders 0..4; the looping default is correct,
      !> concrete LSFs override it with a vectorized kernel (SvdW: strip kernel)
      procedure :: evaluate_batch => lsf_base_evaluate_batch
      !> Number of atoms currently active after `prepare`/`prepare_subset`
      procedure(lsf_active_count_iface), deferred :: active_count
      !> User-space atom id of the i-th currently active atom
      procedure(lsf_active_atom_iface), deferred :: active_atom
      !> LSF value only (lowest-cost path; used by marching cubes)
      procedure(lsf_f0_iface), deferred :: f0
      !> Combined value/gradient/Hessian (any subset via optional args)
      procedure(lsf_f012_r_iface), deferred :: f012_r
      !> Third spatial derivative (plus optionally lower-order outputs)
      procedure(lsf_f3_rrr_iface), deferred :: f3_rrr
      !> Mixed third derivative: spatial Hessian w.r.t. nuclear positions
      procedure(lsf_f3_rr_rA_iface), deferred :: f3_rr_rA

      !* ------------------ Higher-order and contracted derivatives ------------------- *!

      !> Pure spatial fourth derivative d^4S / dr^4
      procedure :: f4_rrrr => lsf_base_f4_rrrr
      !> Mixed fourth derivative d^4S / (dr^3 dR_A)
      procedure :: f4_rrr_rA => lsf_base_f4_rrr_rA
      !> Radius derivative d^3S / (dr^2 dR_a) and its lower orders
      procedure :: f3_rr_rad => lsf_base_f3_rr_rad
      !> Pure nuclear Hessian d^2S / (dR_A dR_B)
      procedure :: f2_rArB => lsf_base_f2_rArB
      !> Mixed third derivative d^3S / (dr dR_A dR_B)
      procedure :: f3_r_rArB => lsf_base_f3_r_rArB
      !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_B)
      procedure :: f4_rr_rArB => lsf_base_f4_rr_rArB
      !> Pure radius Hessian d^2S / (dR_a dR_b)
      procedure :: f2_radrad => lsf_base_f2_radrad
      !> Mixed third derivative d^3S / (dr dR_a dR_b)
      procedure :: f3_r_radrad => lsf_base_f3_r_radrad
      !> Mixed fourth derivative d^4S / (dr^2 dR_a dR_b)
      procedure :: f4_rr_radrad => lsf_base_f4_rr_radrad
      !> Mixed nuclear-radius Hessian d^2S / (dR_A dR_b)
      procedure :: f2_rA_rad => lsf_base_f2_rA_rad
      !> Mixed third derivative d^3S / (dr dR_A dR_b)
      procedure :: f3_r_rA_rad => lsf_base_f3_r_rA_rad
      !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_b)
      procedure :: f4_rr_rA_rad => lsf_base_f4_rr_rA_rad
      !> Normalized level set S/||grad S|| and its nuclear gradient
      procedure :: normalized_f01_rA => lsf_base_normalized_f01_rA

      !* ------------------------------ Atomic partition ------------------------------ *!

      !> Whether this level set splits into atomic shares (no by default)
      procedure :: has_atom_partition => lsf_base_has_atom_partition
      !> Add the atomic shares of the level set at the prepared point
      procedure :: atom_partition => lsf_base_atom_partition

      !* ---------------------- Directional and adjoint products ---------------------- *!

      !> Directional nuclear derivative of the value
      procedure :: tangent_f0 => lsf_base_tangent_f0
      !> Directional nuclear derivative of the spatial gradient
      procedure :: tangent_f1_r => lsf_base_tangent_f1_r
      !> Directional nuclear derivative of the spatial Hessian
      procedure :: tangent_f2_rr => lsf_base_tangent_f2_rr
      !> Directional nuclear derivative of the third spatial derivative
      procedure :: tangent_f3_rrr => lsf_base_tangent_f3_rrr
      !> Nuclear Hessian-vector product
      procedure :: hvp_f1_rA => lsf_base_hvp_f1_rA
      !> Directional nuclear derivative of `f2_r_rA`
      procedure :: hvp_f2_r_rA => lsf_base_hvp_f2_r_rA
      !> Directional nuclear derivative of `f3_rr_rA`
      procedure :: hvp_f3_rr_rA => lsf_base_hvp_f3_rr_rA
      !> Radius row of the joint Hessian-vector product
      procedure :: hvp_f1_rad => lsf_base_hvp_f1_rad
      !> Joint directional derivative of `f2_r_rad`
      procedure :: hvp_f2_r_rad => lsf_base_hvp_f2_r_rad
      !> Joint directional derivative of `f3_rr_rad`
      procedure :: hvp_f3_rr_rad => lsf_base_hvp_f3_rr_rad
      !> Jet-contracted nuclear vector-Jacobian product: the nuclear gradient of
      !> the value/gradient/Hessian jet contracted against per-point adjoint
      !> weights over the spatial indices
      procedure :: vjp_f1_rA => lsf_base_vjp_f1_rA
      !> Jet-contracted radius vector-Jacobian product: the radius gradient of
      !> that same jet contracted against the same per-point adjoint weights
      procedure :: vjp_f1_rad => lsf_base_vjp_f1_rad
      !> Batched `vjp_f1_rA` over a tile, packed symmetric seeds, shared list;
      !> the looping default is correct, concrete LSFs may override it
      procedure :: vjp_f1_rA_batch => lsf_base_vjp_f1_rA_batch
   end type moist_cavity_drop_lsf_type

   !> Slot holding one polymorphic per-thread LSF clone
   !>
   !> - `type(lsf_thread_slot) :: arr(:)`: each element allocated and cloned on its own, unlike a
   !>   `class(...), allocatable :: arr(:)` allocated as a whole with one dynamic type
   type :: lsf_thread_slot
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
   end type lsf_thread_slot

   abstract interface

      !> Escape contract shared with `prepare_subset`, mirroring `evaluate_batch`:
      !>   - `status` present: 0 ok, 1 unreached (no atom within reach, or the
      !>     level-set argument `Z` / `PD` not positive); no `error` for it
      !>   - `status` absent: an unreached point is an `error` ("point is
      !>     unreached (no atom within reach)")
      !>   - both ways the cache is left readable but unreached: every accessor
      !>     then aborts, `active_count` / `active_atom` keep working
      !>   - `error` otherwise only for genuine failures (refused configuration,
      !>     host callback, core fill), `status` then 0
      !>   - isodensity and SSES are never unreached (`status = 0`)
      !>
      !> @param[inout] self   LSF instance
      !> @param[in]    point  Evaluation point (3,)
      !> @param[out]   error  Evaluation failure at this point
      !> @param[out]   status 0 ok, 1 unreached (optional)
      subroutine lsf_prepare_iface(self, point, error, status)
         import :: wp, error_type, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         real(wp), intent(in) :: point(3)
         type(error_type), allocatable, intent(out) :: error
         integer, intent(out), optional :: status
      end subroutine lsf_prepare_iface

      !> Escape contract as [[lsf_prepare_iface]]
      !>
      !> @param[inout] self              LSF instance
      !> @param[in]    point             Evaluation point (3,)
      !> @param[in]    candidate_indices Atom ids to consider, always in the
      !>                                 LSF's declared `candidate_space`; unique
      !>                                 ids in 1..ncenters
      !> @param[out]   error             Evaluation failure at this point
      !> @param[out]   status            0 ok, 1 unreached (optional)
      subroutine lsf_prepare_subset_iface(self, point, candidate_indices, error, status)
         import :: wp, error_type, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         real(wp), intent(in) :: point(3)
         integer, intent(in) :: candidate_indices(:)
         type(error_type), allocatable, intent(out) :: error
         integer, intent(out), optional :: status
      end subroutine lsf_prepare_subset_iface

      !> @param[inout] self LSF instance
      !> @param[in]    n    Requested max derivative order
      subroutine lsf_set_max_deriv_iface(self, n)
         import :: moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         integer, intent(in) :: n
      end subroutine lsf_set_max_deriv_iface

      !> @param[in]  self LSF instance
      !> @returns         Number of currently active atoms
      pure function lsf_active_count_iface(self) result(n)
         import :: moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         integer :: n
      end function lsf_active_count_iface

      !> @param[in]  self LSF instance
      !> @param[in]  i    Active-list index (1 <= i <= active_count())
      !> @returns         User-space atom id of the i-th active atom
      pure function lsf_active_atom_iface(self, i) result(idx)
         import :: moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         integer, intent(in) :: i
         integer :: idx
      end function lsf_active_atom_iface

      !> @param[in]  self LSF instance
      !> @param[out] val  LSF value at the current evaluation point
      subroutine lsf_f0_iface(self, val)
         import :: wp, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out) :: val
      end subroutine lsf_f0_iface

      !> Packed Hessian (6): xx, xy, yy, xz, yz, zz
      !>
      !> @param[in]  self     LSF instance
      !> @param[out] lsf0     LSF value (optional)
      !> @param[out] lsf1_r   Gradient w.r.t. spatial coords (optional)
      !> @param[out] lsf2_rr  Hessian w.r.t. spatial coords (optional)
      subroutine lsf_f012_r_iface(self, lsf0, lsf1_r, lsf2_rr)
         import :: wp, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf0
         real(wp), intent(out), optional :: lsf1_r(:)
         real(wp), intent(out), optional :: lsf2_rr(:)
      end subroutine lsf_f012_r_iface

      !> Packed: `lsf2_rr` (6) xx, xy, yy, xz, yz, zz; `lsf3_rrr` (10) in `packed_index3` order
      !>
      !> @param[in]  self     LSF instance
      !> @param[out] lsf0     LSF value (optional)
      !> @param[out] lsf1_r   Spatial gradient (optional)
      !> @param[out] lsf2_rr  Spatial Hessian (optional)
      !> @param[out] lsf3_rrr Third spatial derivative tensor
      subroutine lsf_f3_rrr_iface(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         import :: wp, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf0
         real(wp), intent(out), optional :: lsf1_r(:)
         real(wp), intent(out), optional :: lsf2_rr(:)
         real(wp), intent(out) :: lsf3_rrr(:)
      end subroutine lsf_f3_rrr_iface

      !> Every nuclear index is an active-list index: slot `i` belongs to the
      !> atom `active_atom(i)`. The caller owns all three buffers and sizes their
      !> nuclear extent from `active_count()`. `lsf3_rr_rA` is packed (6, 3, n_active):
      !> rr slot xx, xy, yy, xz, yz, zz; beta; nuclear slot
      !>
      !> @param[in]  self       LSF instance
      !> @param[out] lsf1_rA    LSF gradient w.r.t. nuclear positions (optional)
      !> @param[out] lsf2_r_rA  Mixed second derivative (optional)
      !> @param[out] lsf3_rr_rA Mixed third derivative
      subroutine lsf_f3_rr_rA_iface(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         import :: wp, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf1_rA(:, :)
         real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
         real(wp), intent(out) :: lsf3_rr_rA(:, :, :)
      end subroutine lsf_f3_rr_rA_iface

      !> Exact radial offset, measured outward from the atom surface, at which the
      !> LSF's contribution from an atom of the given radius equals
      !> `screening_threshold`
      !>
      !> Contract:
      !>   - analytic and exact -- the closed-form inverse of the LSF's own decay,
      !>     not a bound and not a numerical search
      !>   - `>= 0`
      !>   - `huge(0.0_wp)` when screening is disabled (`screening_threshold <= 0`)
      !>     or the parameters make the criterion ill-defined
      !>
      !> Both screening consumers are derived from it: the per-point gate caches
      !> `(radius + screening_offset(radius))**2`, and the candidate cell list
      !> takes `radius + screening_offset(radius)` as its per-atom reach
      !>
      !> @param[in] self    LSF instance
      !> @param[in] radius  Atom radius (Bohr)
      !> @returns           Radial offset from the atom surface (Bohr)
      pure function lsf_screening_offset_iface(self, radius) result(offset)
         import :: wp, moist_cavity_drop_lsf_type
         implicit none(type, external)
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(in) :: radius
         real(wp) :: offset
      end function lsf_screening_offset_iface

   end interface

contains

   !* ================================================================================= *!
   !*                                    Lifecycle                                      *!
   !* ================================================================================= *!

   !> Default `update` implementation: copy geometry and radii into the common LSF
   !> state, then rebuild the spatial sort, the candidate-space geometry mirror and
   !> the screening bounds
   !>
   !> Concrete LSFs that need extra work (e.g. resize their own per-atom caches)
   !> override this and call back to it via `call lsf_base_update(self, mol, radii)`
   !> *first*, since the override typically sizes itself from the state this
   !> routine establishes. (The `self%parent_type%update` form is not available:
   !> the base is abstract.)
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    mol   Molecular structure
   !> @param[in]    radii Per-atom radii (size mol%nat)
   subroutine lsf_base_update(self, mol, radii)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)

      self%mol = mol
      self%ncenters = mol%nat
      if (allocated(self%radii)) deallocate (self%radii)
      self%radii = radii

      call lsf_base_rebuild_screening(self)
   end subroutine lsf_base_update

   !> Move the atom centers in place
   !>
   !> Refreshes the spatial sort, the candidate-space geometry mirror and the
   !> screening bounds, but leaves radii, LSF parameters, derivative order and any
   !> concrete per-atom storage untouched
   !>
   !> - lets nuclear finite-difference drivers displace one coordinate without
   !>   re-initialising the whole LSF
   !>
   !> @param[inout] self    LSF instance (must have been `update`d)
   !> @param[in]    centers Atom centers [ndim, ncenters]
   subroutine lsf_base_set_centers(self, centers)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      real(wp), intent(in) :: centers(:, :)

      if (.not. allocated(self%mol%xyz)) then
         error stop "moist DROP LSF: set_centers before update"
      end if
      if (size(centers, 2) /= self%ncenters) then
         error stop "moist DROP LSF: set_centers atom count mismatch"
      end if

      self%mol%xyz = centers
      call lsf_base_rebuild_screening(self)
   end subroutine lsf_base_set_centers

   !> Rebuild the spatial sort, the candidate-space geometry mirror and the
   !> per-atom screening bounds from the current `mol%xyz` / `radii`
   !>
   !> Runs once per geometry change, never on the per-point path
   !>
   !> @param[inout] self LSF instance
   subroutine lsf_base_rebuild_screening(self)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self

      !> Loop indices and atom count
      integer :: i, j, n

      if (self%candidate_space /= lsf_candidate_space_user .and. &
          self%candidate_space /= lsf_candidate_space_sorted) then
         error stop "moist DROP LSF: candidate_space was never declared -- the "// &
            "concrete LSF's constructor must set it to lsf_candidate_space_user "// &
            "or lsf_candidate_space_sorted"
      end if

      n = self%ncenters
      self%center_lower = 0.0_wp
      self%center_upper = 0.0_wp
      self%center_bounds_valid = n > 0 .and. all(ieee_is_finite(self%mol%xyz))
      if (self%center_bounds_valid) then
         self%center_lower = minval(self%mol%xyz, dim=2)
         self%center_upper = maxval(self%mol%xyz, dim=2)
      end if

      call lsf_base_spatial_sort(self%mol%xyz, n, self%sorted_to_orig, self%orig_to_sorted)

      if (allocated(self%cand_screen)) deallocate (self%cand_screen)
      if (allocated(self%cand_radii)) deallocate (self%cand_radii)
      if (allocated(self%cand_to_user)) deallocate (self%cand_to_user)
      if (allocated(self%full_scan_cand)) deallocate (self%full_scan_cand)
      allocate (self%cand_screen(lsf_cand_screen_rows, n))
      allocate (self%cand_radii(n))
      allocate (self%cand_to_user(n))
      allocate (self%full_scan_cand(n))

      select case (self%candidate_space)
      case (lsf_candidate_space_sorted)
         do j = 1, n
            i = self%sorted_to_orig(j)
            self%cand_screen(lsf_cand_center(1):lsf_cand_center(2), j) = self%mol%xyz(:, i)
            self%cand_radii(j) = self%radii(i)
            self%cand_to_user(j) = i
         end do
         ! Visiting the full scan through orig_to_sorted keeps the resulting
         ! active list in ascending user-space order, so an unscreened `prepare`
         ! sums in the same order regardless of how the spatial sort came out
         do i = 1, n
            self%full_scan_cand(i) = self%orig_to_sorted(i)
         end do
      case default
         do i = 1, n
            self%cand_screen(lsf_cand_center(1):lsf_cand_center(2), i) = self%mol%xyz(:, i)
            self%cand_radii(i) = self%radii(i)
            self%cand_to_user(i) = i
            self%full_scan_cand(i) = i
         end do
      end select

      call lsf_base_refresh_bounds(self)
   end subroutine lsf_base_rebuild_screening

   !> Recompute the cached squared screening reach of every candidate and rebuild
   !> the candidate cell list from it
   !>
   !> Split out of [[lsf_base_rebuild_screening]] because it is the only part
   !> that depends on `screening_threshold`: it does not re-sort, so
   !> [[lsf_base_set_screening_threshold]] can re-run it alone
   !>
   !> - the cell list takes the unsquared reach and saturates it at `sqrt(huge)`
   !>   exactly as the gate does, so both see bitwise the same squared reach
   !>
   !> @param[inout] self LSF instance with the candidate mirror already built
   subroutine lsf_base_refresh_bounds(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: self

      !> Screening reach of one atom
      real(wp) :: bound
      !> Screening reach per candidate, handed to the cell list
      real(wp), allocatable :: reach(:)
      integer :: j

      if (.not. allocated(self%cand_screen)) return
      allocate (reach(size(self%cand_radii)))

      do j = 1, size(self%cand_radii)
         bound = self%cand_radii(j) + self%screening_offset(self%cand_radii(j))
         reach(j) = bound
         self%cand_screen(lsf_cand_bound_sq, j) = saturated_square(bound)
      end do

      call self%cand_cells%build(self%cand_screen(lsf_cand_center(1):lsf_cand_center(2), :), &
                                 reach, self%cell_full_scan_below)
   end subroutine lsf_base_refresh_bounds

   !> Squared screening reach as the gate caches it
   !>
   !> - saturates at `huge` from `sqrt(huge)` on, so "screening disabled" never
   !>   turns into an overflow
   !>
   !> @param[in] bound Screening reach `R + offset(R)`
   !> @returns   sq    Cached squared reach
   elemental function saturated_square(bound) result(sq)
      !> Screening reach
      real(wp), intent(in) :: bound
      !> Cached squared reach
      real(wp) :: sq

      if (bound >= sqrt(huge(0.0_wp))) then
         sq = huge(0.0_wp)
      else
         sq = bound*bound
      end if
   end function saturated_square

   !> Whether the cached reach still equals the live `screening_offset`
   !>
   !> - recomputes candidate 1's squared reach from the live parameters exactly
   !>   as [[lsf_base_refresh_bounds]] does and compares bitwise
   !> - one candidate suffices: every cached reach comes from one refresh of one
   !>   formula in `param` and `R_j`, so a stale refresh shows on any of them;
   !>   O(1), cheap enough for every `prepare`
   !> - catches a refresh made under other parameters, e.g. `param` changed,
   !>   `set_screening_threshold` called, `param` restored
   !> - `.true.` without candidates, `.false.` before `update`
   !>
   !> @param[in] self    LSF instance
   !> @returns   current Cached reach built for the live parameters
   pure function lsf_base_reach_current(self) result(current)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Cached reach built for the live parameters
      logical :: current

      current = .false.
      if (.not. allocated(self%cand_screen)) return
      current = .true.
      if (size(self%cand_radii) == 0) return
      current = self%cand_screen(lsf_cand_bound_sq, 1) == &
                saturated_square(self%cand_radii(1) + self%screening_offset(self%cand_radii(1)))
   end function lsf_base_reach_current

   !> Change the screening threshold in place
   !>
   !> Only the cached per-candidate reach depends on the threshold, so this
   !> rewrites one row of `cand_screen`, rebuilds the candidate cell list and
   !> leaves the spatial sort and the candidate mirror untouched
   !>
   !> - cheap enough to bracket a solve with: O(ncenters) work
   !> - tighter and looser thresholds both safe, the cell list taking the new reach
   !> - candidate lists queried before the call are stale
   !>
   !> @param[inout] self      LSF instance
   !> @param[in]    threshold New screening threshold
   subroutine lsf_base_set_screening_threshold(self, threshold)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      !> New screening threshold
      real(wp), intent(in) :: threshold

      if (threshold == self%screening_threshold) return
      self%screening_threshold = threshold
      call lsf_base_refresh_bounds(self)
   end subroutine lsf_base_set_screening_threshold

   !> Change the candidate count below which the candidate cell list is one cell
   !>
   !> - stored, and the candidate cell list rebuilt when `update` already built
   !>   one; otherwise read by the next `update`
   !> - both modes return the same candidates, so only speed changes
   !> - `error stop` on a negative count
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    n    Candidate count threshold, >= 0; 0 never collapses to one cell
   subroutine lsf_base_set_full_scan_below(self, n)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      !> Candidate count threshold
      integer, intent(in) :: n

      if (n < 0) error stop "moist DROP LSF: set_full_scan_below needs n >= 0"
      if (n == self%cell_full_scan_below) return
      self%cell_full_scan_below = n
      call lsf_base_refresh_bounds(self)
   end subroutine lsf_base_set_full_scan_below

   !> Coarse spatial bucket sort of the atom centers
   !>
   !> Assigns each atom to one of 6^3 = 216 coarse 3D buckets and builds the
   !> permutation with a single-pass counting sort (O(N), stable), inspired by
   !> stdlib's int8 radix_sort
   !>
   !> - sorting the geometry this way turns the gathers driven by cell-local
   !>   candidate lists into near-sequential reads
   !>
   !> @param[in]    xyz            Atom centers [ndim, n]
   !> @param[in]    n              Number of atoms
   !> @param[inout] sorted_to_orig User-space id of sorted position j [n]
   !> @param[inout] orig_to_sorted Sorted position of user-space id i [n]
   pure subroutine lsf_base_spatial_sort(xyz, n, sorted_to_orig, orig_to_sorted)
      real(wp), intent(in) :: xyz(:, :)
      integer, intent(in) :: n
      integer, allocatable, intent(inout) :: sorted_to_orig(:)
      integer, allocatable, intent(inout) :: orig_to_sorted(:)

      !> Loop indices
      integer :: i, j
      !> Per-axis minima / extents used to discretize coordinates for bucket assignment
      real(wp) :: xmin(ndim), xmax(ndim), span(ndim)
      !> Discretized per-axis bucket index
      integer :: bx, by, bz
      !> Buckets per axis for spatial grouping (6^3 = 216 buckets, fits in 0..255)
      integer, parameter :: nbx = 6, nby = 6, nbz = 6
      !> Total number of spatial buckets
      integer, parameter :: n_buckets = nbx*nby*nbz
      !> Per-atom spatial bucket (0..n_buckets-1)
      integer, allocatable :: buckets(:)

      if (allocated(sorted_to_orig)) deallocate (sorted_to_orig)
      if (allocated(orig_to_sorted)) deallocate (orig_to_sorted)
      allocate (sorted_to_orig(n))
      allocate (orig_to_sorted(n))

      if (n <= 0) return
      if (n == 1) then
         sorted_to_orig(1) = 1
         orig_to_sorted(1) = 1
         return
      end if

      allocate (buckets(n))
      xmin = minval(xyz(:, 1:n), dim=2)
      xmax = maxval(xyz(:, 1:n), dim=2)
      span = xmax - xmin

      do i = 1, n
         if (span(1) > 0.0_wp) then
            bx = min(nbx - 1, int((xyz(1, i) - xmin(1))/span(1)*nbx))
         else
            bx = 0
         end if
         if (span(2) > 0.0_wp) then
            by = min(nby - 1, int((xyz(2, i) - xmin(2))/span(2)*nby))
         else
            by = 0
         end if
         if (span(3) > 0.0_wp) then
            bz = min(nbz - 1, int((xyz(3, i) - xmin(3))/span(3)*nbz))
         else
            bz = 0
         end if
         buckets(i) = bx + by*nbx + bz*nbx*nby
      end do

      call counting_argsort(buckets, n_buckets - 1, sorted_to_orig)
      do j = 1, n
         orig_to_sorted(sorted_to_orig(j)) = j
      end do
      deallocate (buckets)
   end subroutine lsf_base_spatial_sort

   !> Abort when an accessor is asked for a derivative order the latest
   !> `prepare` did not compute
   !>
   !> @param[in] self   LSF instance
   !> @param[in] order  Derivative order the caller is about to read
   !> @param[in] caller Accessor name, for the diagnostic
   subroutine lsf_base_require_deriv(self, order, caller)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      integer, intent(in) :: order
      character(len=*), intent(in) :: caller

      character(len=16) :: got, want

      call self%require_reached(caller)
      if (self%prepared_deriv < 0) then
         error stop "moist DROP LSF: "//caller//" read the per-point cache before "// &
            "any successful prepare -- the level-set jet is not a result"
      end if
      if (order <= self%prepared_deriv) return

      write (want, "(i0)") order
      write (got, "(i0)") self%prepared_deriv
      error stop "moist DROP LSF: "//caller//" needs derivative order "// &
         trim(want)//" but prepare ran at "//trim(got)// &
         " -- raise set_max_deriv before prepare"
   end subroutine lsf_base_require_deriv

   !> Abort when an accessor reads the cache of an unreached prepare
   !>
   !> - the escape has no value, so an accessor returns nothing rather than a
   !>   zero the caller cannot tell from a result
   !> - [[lsf_base_require_deriv]] calls it first; accessors without an order
   !>   check call it directly
   !>
   !> @param[in] self   LSF instance
   !> @param[in] caller Accessor name, for the diagnostic
   pure subroutine lsf_base_require_reached(self, caller)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Accessor name
      character(len=*), intent(in) :: caller

      if (.not. self%prepared_unreached) return
      error stop "moist DROP LSF: accessor called after an unreached prepare ("//caller//")"
   end subroutine lsf_base_require_reached

   !> Whether the point of the latest `prepare` / `prepare_subset` is unreached
   !>
   !> - unreached: no atom within reach, or `Z` / `PD` not positive; the LSF has
   !>   no value there (escape)
   !> - never for isodensity (an empty set is a defined value there) or SSES
   !> - the same flag `prepare`'s `status = 1` reports
   !>
   !> @param[in] self      LSF instance after `prepare`
   !> @returns   unreached Point has no atom in reach
   pure function lsf_base_unreached(self) result(unreached)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Point has no atom in reach
      logical :: unreached

      unreached = self%prepared_unreached
   end function lsf_base_unreached

   !> Report the latest prepare's escape: `status = 1`, or an `error` without `status`
   !>
   !> - called last by every concrete `prepare` that can be unreached, after
   !>   `prepared_unreached` is set
   !> - the message says "no atom within reach", qualified when atoms are active
   !>   and the lifted jet is not finite instead
   !> - leaves `status` and `error` alone for a reached point
   !>
   !> @param[in]    self   LSF instance after the screen
   !> @param[inout] error  Set for an unreached point without `status`
   !> @param[inout] status Set to 1 for an unreached point (optional)
   subroutine lsf_base_report_unreached(self, error, status)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Set for an unreached point without `status`
      type(error_type), allocatable, intent(inout) :: error
      !> Set to 1 for an unreached point
      integer, intent(inout), optional :: status

      if (.not. self%prepared_unreached) return
      if (present(status)) then
         status = 1
         return
      end if
      if (self%active_count() > 0) then
         call fatal_error(error, "prepare: point is unreached (no atom within reach) or the "// &
                          "value is not finite with atoms in reach ("//self%name//")")
      else
         call fatal_error(error, "prepare: point is unreached (no atom within reach) ("// &
                          self%name//")")
      end if
   end subroutine lsf_base_report_unreached

   !> Level set jet at a batch of points, up to spatial derivative `order`
   !>
   !> Contract shared by the default and every vectorized override:
   !>   - `points(3, n)` contiguous; outputs are filled for columns `1:n` only
   !>   - `order` 0..4 selects which outputs are written: 0 -> `lsf0`, 1 -> `lsf1_r`
   !>     too, 2 -> packed `lsf2_rr_packed(6, n)` too, 3 -> packed `lsf3_rrr(10, n)` too,
   !>     4 -> packed `lsf4_rrrr(15, n)` too; an output above `order` is left
   !>     untouched, a missing or misshaped output for the requested order aborts
   !>   - Hessian packed `lsf2_rr_packed(6, n)`: xx, xy, yy, xz, yz, zz;
   !>     `unpack_sym_matrix` from [[moist_math_packed_sym]] for a dense view
   !>   - packed slots: sorted tuples `a <= b <= c <= d`, last index slowest,
   !>     [[moist_math_packed_sym]]'s `packed_index3` / `packed_index4`
   !>   - `candidate_indices`, when present, is the atom list for *all* points
   !>     (octree children reuse the parent's list, seeds of one anchor share the
   !>     anchor's neighbourhood): unique ids in 1..ncenters (`error stop`
   !>     otherwise), sorted ascending before use; absent means the full screen
   !>     per point
   !>   - optional `status(n)`: 0 ok, 1 unreached (no atom within reach), 2 failed
   !>     or not computed (a host error, or a column at or behind the one an
   !>     `error` stopped at); every entry `1:n` is defined on every return;
   !>     unreached and failed columns leave their outputs unspecified; with
   !>     `status` absent an unreached column is an `error`
   !>   - `error` names the first failing column and leaves the remaining columns
   !>     unspecified, their `status` 2
   !>   - the single-point accessor cache (`f0`, `f012_r`, ...) is unspecified
   !>     after this call; re-`prepare` before reading it
   !>
   !> Looping default:
   !>   - `prepare` / `prepare_subset` (on one ascending copy of the list) plus
   !>     the per-point accessors for every column; `order` above the prepared `set_max_deriv` aborts
   !>     (`require_deriv`), so a caller raises `set_max_deriv(order)` first
   !>   - orders 2..4 call the packed `f012_r` / `f3_rrr` / `f4_rrrr` straight into the columns
   !>   - unreached as the prepare's `status = 1` says (CFC: no atom in reach;
   !>     isodensity: never, an empty set is a defined value there): `status = 1`
   !>     and every output of the column through `order` a quiet NaN, or `error`
   !>     naming the column when `status` is absent
   !>   - `status` starts at 2 for every column, so a `prepare` error leaves the
   !>     failed column and those behind it at 2
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    points            Evaluation points (3, n)
   !> @param[in]    order             Highest spatial derivative to write (0..4)
   !> @param[out]   lsf0              Values (n)
   !> @param[out]   lsf1_r            Gradients (3, n), written for order >= 1
   !> @param[out]   lsf2_rr_packed    Packed Hessians (6, n), written for order >= 2
   !> @param[out]   lsf3_rrr          Packed third derivatives (10, n), written for order >= 3
   !> @param[out]   lsf4_rrrr         Packed fourth derivatives (15, n), written for order 4
   !> @param[out]   status            Per-column status (n), 0 ok, 1 unreached, 2 failed
   !>                                 or not computed
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   !> @param[out]   error             First failing point, by column
   subroutine lsf_base_evaluate_batch(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                      lsf3_rrr, lsf4_rrrr, status, candidate_indices, error)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      real(wp), intent(in), contiguous :: points(:, :)
      integer, intent(in) :: order
      real(wp), intent(out), contiguous :: lsf0(:)
      real(wp), intent(out), optional, contiguous :: lsf1_r(:, :)
      real(wp), intent(out), optional, contiguous :: lsf2_rr_packed(:, :)
      real(wp), intent(out), optional, contiguous :: lsf3_rrr(:, :)
      real(wp), intent(out), optional, contiguous :: lsf4_rrrr(:, :)
      integer, intent(out), optional, contiguous :: status(:)
      integer, intent(in), optional :: candidate_indices(:)
      type(error_type), allocatable, intent(out) :: error

      !> Fill of an unreached column
      real(wp) :: nan
      !> Ascending copy of `candidate_indices(1:nc)`
      integer, allocatable :: sorted(:)
      !> Column, point count, list size and the prepare's escape status
      integer :: i, n, nc, pstat
      character(len=16) :: col

      n = size(points, 2)
      call lsf_batch_require_outputs(self%name, points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                     lsf3_rrr, lsf4_rrrr, status)
      if (present(status)) status(1:n) = 2
      if (present(candidate_indices)) then
         allocate (sorted(size(candidate_indices)))
         call lsf_sorted_candidates(candidate_indices, sorted, nc, self%ncenters)
      end if

      do i = 1, n
         if (present(candidate_indices)) then
            call self%prepare_subset(points(:, i), sorted(1:nc), error, pstat)
         else
            call self%prepare(points(:, i), error, pstat)
         end if
         if (allocated(error)) then
            write (col, "(i0)") i
            error%message = "evaluate_batch failed at point "//trim(col)//": "//error%message
            return
         end if
         if (pstat /= 0) then
            write (col, "(i0)") i
            if (.not. present(status)) then
               call fatal_error(error, "evaluate_batch: point "//trim(col)// &
                                " is unreached, no atom within reach ("//self%name//")")
               return
            end if
            status(i) = 1
            nan = ieee_value(0.0_wp, ieee_quiet_nan)
            lsf0(i) = nan
            if (order >= 1) lsf1_r(:, i) = nan
            if (order >= 2) lsf2_rr_packed(:, i) = nan
            if (order >= 3) lsf3_rrr(:, i) = nan
            if (order >= 4) lsf4_rrrr(:, i) = nan
            cycle
         end if
         call self%require_deriv(order, "evaluate_batch")
         select case (order)
         case (0)
            call self%f012_r(lsf0=lsf0(i))
         case (1)
            call self%f012_r(lsf0=lsf0(i), lsf1_r=lsf1_r(:, i))
         case (2)
            call self%f012_r(lsf0=lsf0(i), lsf1_r=lsf1_r(:, i), lsf2_rr=lsf2_rr_packed(:, i))
         case (3)
            call self%f3_rrr(lsf0=lsf0(i), lsf1_r=lsf1_r(:, i), lsf2_rr=lsf2_rr_packed(:, i), &
                             lsf3_rrr=lsf3_rrr(:, i))
         case default
            call self%f3_rrr(lsf0=lsf0(i), lsf1_r=lsf1_r(:, i), lsf2_rr=lsf2_rr_packed(:, i), &
                             lsf3_rrr=lsf3_rrr(:, i))
            call self%f4_rrrr(lsf4_rrrr(:, i))
         end select
         if (present(status)) status(i) = 0
      end do
   end subroutine lsf_base_evaluate_batch

   !> Abort on an `evaluate_batch` call whose outputs cannot hold `order`
   !>
   !> Shared by the looping default and every override, so the aborts read the
   !> same everywhere
   !>
   !> - `points` with a leading extent other than 3
   !> - `order` outside 0..4
   !> - an output of the requested order absent, or with a leading extent other
   !>   than its width (3, 6, 10, 15), or with fewer than `n` columns
   !> - `status` present with fewer than `n` entries
   !>
   !> @param[in] who       LSF name for the diagnostic
   !> @param[in] points    Evaluation points (3, n); only the extents are read
   !> @param[in] order     Requested derivative order
   !> @param[in] lsf0      Values
   !> @param[in] lsf1_r    Gradients (optional)
   !> @param[in] lsf2_rr_packed Packed Hessians (optional)
   !> @param[in] lsf3_rrr  Packed third derivatives (optional)
   !> @param[in] lsf4_rrrr Packed fourth derivatives (optional)
   !> @param[in] status    Per-column status (optional)
   subroutine lsf_batch_require_outputs(who, points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                        lsf3_rrr, lsf4_rrrr, status)
      !> LSF name for the diagnostic
      character(len=*), intent(in) :: who
      !> Evaluation points; only the extents are read
      real(wp), intent(in) :: points(:, :)
      !> Requested derivative order
      integer, intent(in) :: order
      !> Values; only the extent is read
      real(wp), intent(in) :: lsf0(:)
      !> Gradients; only the extents are read
      real(wp), intent(in), optional :: lsf1_r(:, :)
      !> Packed Hessians; only the extents are read
      real(wp), intent(in), optional :: lsf2_rr_packed(:, :)
      !> Packed third derivatives; only the extents are read
      real(wp), intent(in), optional :: lsf3_rrr(:, :)
      !> Packed fourth derivatives; only the extents are read
      real(wp), intent(in), optional :: lsf4_rrrr(:, :)
      !> Per-column status; only the extent is read
      integer, intent(in), optional :: status(:)

      character(len=*), parameter :: pre = "moist DROP LSF evaluate_batch: "
      !> Number of points
      integer :: n

      if (size(points, 1) /= ndim) error stop pre//"points must be (3, n) ("//who//")"
      n = size(points, 2)
      if (order < 0 .or. order > 4) error stop pre//"order must be 0..4 ("//who//")"
      if (size(lsf0) < n) error stop pre//"lsf0 holds fewer than n columns ("//who//")"
      if (order >= 1) then
         if (.not. present(lsf1_r)) error stop pre//"order >= 1 needs lsf1_r ("//who//")"
         if (size(lsf1_r, 1) /= ndim .or. size(lsf1_r, 2) < n) then
            error stop pre//"lsf1_r must be (3, >= n) ("//who//")"
         end if
      end if
      if (order >= 2) then
         if (.not. present(lsf2_rr_packed)) then
            error stop pre//"order >= 2 needs lsf2_rr_packed ("//who//")"
         end if
         if (size(lsf2_rr_packed, 1) /= 6 .or. size(lsf2_rr_packed, 2) < n) then
            error stop pre//"lsf2_rr_packed must be packed (6, >= n) ("//who//")"
         end if
      end if
      if (order >= 3) then
         if (.not. present(lsf3_rrr)) error stop pre//"order >= 3 needs lsf3_rrr ("//who//")"
         if (size(lsf3_rrr, 1) /= 10 .or. size(lsf3_rrr, 2) < n) then
            error stop pre//"lsf3_rrr must be packed (10, >= n) ("//who//")"
         end if
      end if
      if (order >= 4) then
         if (.not. present(lsf4_rrrr)) error stop pre//"order 4 needs lsf4_rrrr ("//who//")"
         if (size(lsf4_rrrr, 1) /= 15 .or. size(lsf4_rrrr, 2) < n) then
            error stop pre//"lsf4_rrrr must be packed (15, >= n) ("//who//")"
         end if
      end if
      if (present(status)) then
         if (size(status) < n) error stop pre//"status holds fewer than n entries ("//who//")"
      end if
   end subroutine lsf_batch_require_outputs

   !> Copy a shared candidate list into `out(1:n)`, ascending
   !>
   !> The batch contract's canonical order: atoms outside a point's reach add
   !> exact zeros in that order, which is what makes tiling invisible
   !>
   !> - already ascending lists are copied, others heapsorted in place
   !> - contract: unique ids in 1..ncenters; a repeated or out-of-range id is a
   !>   caller bug and aborts (`error stop`), checked after the sort
   !> - allocation-free; `size(out) >= size(ids)`
   !>
   !> @param[in]  ids      Candidate ids in `candidate_space`
   !> @param[out] out      Ascending copy in `out(1:n)`
   !> @param[out] n        Number of ids
   !> @param[in]  ncenters Largest valid id
   pure subroutine lsf_sorted_candidates(ids, out, n, ncenters)
      !> Candidate ids
      integer, intent(in) :: ids(:)
      !> Ascending copy
      integer, intent(out) :: out(:)
      !> Number of ids
      integer, intent(out) :: n
      !> Largest valid id
      integer, intent(in) :: ncenters

      integer :: i
      logical :: ascending

      n = size(ids)
      if (size(out) < n) error stop "moist DROP LSF: candidate buffer smaller than the list"
      ascending = .true.
      do i = 1, n
         out(i) = ids(i)
         if (i > 1) ascending = ascending .and. ids(i - 1) <= ids(i)
      end do
      if (.not. ascending) call heapsort(out(1:n))
      if (n == 0) return
      if (out(1) < 1 .or. out(n) > ncenters) then
         error stop "moist DROP LSF: candidate id outside 1..ncenters"
      end if
      do i = 2, n
         if (out(i) == out(i - 1)) error stop "moist DROP LSF: repeated candidate id"
      end do
   end subroutine lsf_sorted_candidates

   !* ================================================================================= *!
   !*                Erroring defaults of the optional derivative set                   *!
   !* ================================================================================= *!
   !
   ! Returning zeros here would be a silent wrong answer, the caller being unable
   ! to tell a genuine vanishing tensor from a missing implementation:
   !
   ! - each default aborts with its own `error stop`, spelled out in full so the
   !   diagnostic naming the missing derivative sits at the line that raises it
   ! - each result is still zeroed on the line before the abort, with the dummy
   !   arguments folded into that store: neither is ever observable (the next
   ! statement terminates the program) and both exist only so that
   ! `-Wunused-dummy-argument` does not have to be switched off here

   !> Erroring default of the pure spatial fourth derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Fourth spatial derivative, packed (15)
   subroutine lsf_base_f4_rrrr(self, lsf4_rrrr)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rrrr(:)

      lsf4_rrrr = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rrrr derivative"
   end subroutine lsf_base_f4_rrrr

   !> Erroring default of the mixed fourth derivative d^4S / (dr^3 dR_A)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf4_rrr_rA Mixed fourth derivative, packed rrr [10, 3, n_active]
   subroutine lsf_base_f4_rrr_rA(self, lsf4_rrr_rA)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rrr_rA(:, :, :)

      lsf4_rrr_rA = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rrr_rA derivative"
   end subroutine lsf_base_f4_rrr_rA

   !> Erroring default of the radius derivative d^3S / (dr^2 dR_a)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf1_rad    dS/dR_a [>= n_active] (optional)
   !> @param[out] lsf2_r_rad  d^2S/(dr dR_a) [3, >= n_active] (optional)
   !> @param[out] lsf3_rr_rad d^3S/(dr^2 dR_a), packed rr [6, >= n_active]
   subroutine lsf_base_f3_rr_rad(self, lsf1_rad, lsf2_r_rad, lsf3_rr_rad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rad(:)
      real(wp), intent(out), optional :: lsf2_r_rad(:, :)
      real(wp), intent(out) :: lsf3_rr_rad(:, :)

      if (present(lsf1_rad)) lsf1_rad = 0.0_wp
      if (present(lsf2_r_rad)) lsf2_r_rad = 0.0_wp
      lsf3_rr_rad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f3_rr_rad derivative"
   end subroutine lsf_base_f3_rr_rad

   !> Erroring default of the pure nuclear Hessian
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf2_rArB Nuclear Hessian [3, n_active, 3, n_active]
   subroutine lsf_base_f2_rArB(self, lsf2_rArB)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf2_rArB(:, :, :, :)

      lsf2_rArB = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f2_rArB derivative"
   end subroutine lsf_base_f2_rArB

   !> Erroring default of the mixed third derivative d^3S / (dr dR_A dR_B)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf3_r_rArB Mixed third derivative [3, 3, n_active, 3, n_active]
   subroutine lsf_base_f3_r_rArB(self, lsf3_r_rArB)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf3_r_rArB(:, :, :, :, :)

      lsf3_r_rArB = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f3_r_rArB derivative"
   end subroutine lsf_base_f3_r_rArB

   !> Erroring default of the mixed fourth derivative d^4S / (dr^2 dR_A dR_B)
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf4_rr_rArB Mixed fourth derivative, packed rr [6, 3, n_act, 3, n_act]
   subroutine lsf_base_f4_rr_rArB(self, lsf4_rr_rArB)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rr_rArB(:, :, :, :, :)

      lsf4_rr_rArB = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rr_rArB derivative"
   end subroutine lsf_base_f4_rr_rArB

   !* ================================================================================= *!
   !*                     Two-radius and nuclear-radius derivatives                     *!
   !* ================================================================================= *!

   !> Erroring default of the pure radius Hessian
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_radrad d^2S/(dR_a dR_b) [>= n_act, >= n_act]
   subroutine lsf_base_f2_radrad(self, lsf2_radrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf2_radrad(:, :)

      lsf2_radrad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f2_radrad derivative"
   end subroutine lsf_base_f2_radrad

   !> Erroring default of the mixed third derivative d^3S / (dr dR_a dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_radrad d^3S/(dr dR_a dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_base_f3_r_radrad(self, lsf3_r_radrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf3_r_radrad(:, :, :)

      lsf3_r_radrad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f3_r_radrad derivative"
   end subroutine lsf_base_f3_r_radrad

   !> Erroring default of the mixed fourth derivative d^4S / (dr^2 dR_a dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b), packed rr [6, >= n_act, >= n_act]
   subroutine lsf_base_f4_rr_radrad(self, lsf4_rr_radrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rr_radrad(:, :, :)

      lsf4_rr_radrad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rr_radrad derivative"
   end subroutine lsf_base_f4_rr_radrad

   !> Erroring default of the mixed nuclear-radius Hessian
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_rA_rad d^2S/(dR_A dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_base_f2_rA_rad(self, lsf2_rA_rad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf2_rA_rad(:, :, :)

      lsf2_rA_rad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f2_rA_rad derivative"
   end subroutine lsf_base_f2_rA_rad

   !> Erroring default of the mixed third derivative d^3S / (dr dR_A dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_rA_rad d^3S/(dr dR_A dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_base_f3_r_rA_rad(self, lsf3_r_rA_rad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf3_r_rA_rad(:, :, :, :)

      lsf3_r_rA_rad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f3_r_rA_rad derivative"
   end subroutine lsf_base_f3_r_rA_rad

   !> Erroring default of the mixed fourth derivative d^4S / (dr^2 dR_A dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b), packed rr [6, 3, >= n_act, >= n_act]
   subroutine lsf_base_f4_rr_rA_rad(self, lsf4_rr_rA_rad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rr_rA_rad(:, :, :, :)

      lsf4_rr_rA_rad = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rr_rA_rad derivative"
   end subroutine lsf_base_f4_rr_rA_rad

   !> Erroring default of the normalized level set and its nuclear gradient
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] val      S/||grad S||
   !> @param[out] deriv_rA Nuclear gradient of the normalized value (optional)
   subroutine lsf_base_normalized_f01_rA(self, val, deriv_rA)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: val
      real(wp), intent(out), optional :: deriv_rA(:, :)

      val = 0.0_wp*self%ncenters
      if (present(deriv_rA)) deriv_rA = 0.0_wp
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "normalized_f01_rA derivative"
   end subroutine lsf_base_normalized_f01_rA

   !> Erroring default of the contracted value derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted value derivative
   subroutine lsf_base_tangent_f0(self, v, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res

      res = 0.0_wp*size(v, 2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "tangent_f0 derivative"
   end subroutine lsf_base_tangent_f0

   !> Erroring default of the contracted gradient derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted gradient derivative [3]
   subroutine lsf_base_tangent_f1_r(self, v, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:)

      res = 0.0_wp*size(v, 2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "tangent_f1_r derivative"
   end subroutine lsf_base_tangent_f1_r

   !> Erroring default of the contracted Hessian derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted Hessian derivative, packed (6)
   subroutine lsf_base_tangent_f2_rr(self, v, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:)

      res = 0.0_wp*size(v, 2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "tangent_f2_rr derivative"
   end subroutine lsf_base_tangent_f2_rr

   !> Erroring default of the contracted third-derivative derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted third-derivative derivative, packed (10)
   subroutine lsf_base_tangent_f3_rrr(self, v, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:)

      res = 0.0_wp*size(v, 2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "tangent_f3_rrr derivative"
   end subroutine lsf_base_tangent_f3_rrr

   !> Erroring default of the nuclear Hessian-vector product
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted nuclear Hessian [3, n_active]
   !> @param[in]  vrad Radius directions [ncenters] (optional)
   subroutine lsf_base_hvp_f1_rA(self, v, res, vrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:, :)
      real(wp), intent(in), optional :: vrad(:)

      res = 0.0_wp*size(v, 2)*self%ncenters
      if (present(vrad)) res = res*size(vrad)
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f1_rA derivative"
   end subroutine lsf_base_hvp_f1_rA

   !> Erroring default of the contracted mixed third derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted mixed third derivative [3, 3, n_active]
   !> @param[in]  vrad Radius directions [ncenters] (optional)
   subroutine lsf_base_hvp_f2_r_rA(self, v, res, vrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:, :, :)
      real(wp), intent(in), optional :: vrad(:)

      res = 0.0_wp*size(v, 2)*self%ncenters
      if (present(vrad)) res = res*size(vrad)
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f2_r_rA derivative"
   end subroutine lsf_base_hvp_f2_r_rA

   !> Erroring default of the contracted mixed fourth derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted mixed fourth derivative, packed rr [6, 3, n_active]
   !> @param[in]  vrad Radius directions [ncenters] (optional)
   subroutine lsf_base_hvp_f3_rr_rA(self, v, res, vrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:, :, :)
      real(wp), intent(in), optional :: vrad(:)

      res = 0.0_wp*size(v, 2)*self%ncenters
      if (present(vrad)) res = res*size(vrad)
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f3_rr_rA derivative"
   end subroutine lsf_base_hvp_f3_rr_rA

   !* ================================================================================= *!
   !*                      Radius row of the joint Hessian-vector product               *!
   !* ================================================================================= *!

   !> Erroring default of the radius Hessian-vector product
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted radius Hessian row [n_active]
   subroutine lsf_base_hvp_f1_rad(self, v, vrad, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(in) :: vrad(:)
      real(wp), intent(out) :: res(:)

      res = 0.0_wp*size(v, 2)*size(vrad)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f1_rad derivative"
   end subroutine lsf_base_hvp_f1_rad

   !> Erroring default of the contracted mixed spatial-radius third derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed third derivative [3, n_active]
   subroutine lsf_base_hvp_f2_r_rad(self, v, vrad, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(in) :: vrad(:)
      real(wp), intent(out) :: res(:, :)

      res = 0.0_wp*size(v, 2)*size(vrad)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f2_r_rad derivative"
   end subroutine lsf_base_hvp_f2_r_rad

   !> Erroring default of the contracted mixed spatial-radius fourth derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed fourth derivative, packed rr [6, n_active]
   subroutine lsf_base_hvp_f3_rr_rad(self, v, vrad, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(in) :: vrad(:)
      real(wp), intent(out) :: res(:, :)

      res = 0.0_wp*size(v, 2)*size(vrad)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f3_rr_rad derivative"
   end subroutine lsf_base_hvp_f3_rr_rad

   !* ================================================================================= *!
   !*                      Jet-contracted vector-Jacobian products                      *!
   !* ================================================================================= *!

   !> Erroring default of the jet-contracted nuclear vector-Jacobian product
   !>
   !> Reverse-mode mirror of the `tangent_*` family: where those contract the
   !> *nuclear* index against a displacement direction, this one contracts the
   !> *spatial* (jet) indices against per-point adjoint weights and keeps the
   !> nuclear index, so no implementation ever has to materialize the full
   !> `[6, 3, n_active]` packed tensor `f3_rr_rA` returns
   !>
   !> The name follows the family convention: this is the `f1_rA` rung, the
   !> nuclear-gradient row, just as `hvp_f1_rA` is that same rung contracted
   !> with a nuclear direction instead
   !>
   !> The result an implementation must produce is
   !>
   !>     res(beta, i) = w0 * lsf1_rA(beta, i)
   !>                  + sum_a w1(a) * lsf2_r_rA(a, beta, i)
   !>                  + sum_{a<=b} m_ab * w2(ab) * lsf3_rr_rA(ab, beta, i)
   !>
   !> for `i = 1 .. active_count()`, with the same compact active-index
   !> convention as `f3_rr_rA`, slot `i` belonging to atom `active_atom(i)`
   !>
   !> - columns beyond `active_count()` are left untouched, and nothing is written
   !>   when `active_count()` is zero
   !> - `w2` is the packed symmetric seed (6): xx, xy, yy, xz, yz, zz
   !> - off-diagonals are not doubled; the contraction applies `m_ab` = 1 on the
   !>   diagonal and 2 off it
   !> - `pack_sym_part` of [[moist_math_packed_sym]] packs a general seed
   !> - an implementation needs whatever prepared order its own `f3_rr_rA` needs,
   !>   2 for SvdW and 3 for CFC, which asks for the highest *total* order
   !> - must call `require_deriv(<that order>, "vjp_f1_rA")`
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Adjoint weights of the spatial Hessian, packed (6)
   !> @param[out] res  Contracted nuclear gradient [3, >= n_active]
   subroutine lsf_base_vjp_f1_rA(self, w0, w1, w2, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: w0
      real(wp), intent(in) :: w1(3)
      real(wp), intent(in) :: w2(:)
      real(wp), intent(out) :: res(:, :)

      res = 0.0_wp*w0*size(w1)*size(w2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "vjp_f1_rA derivative"
   end subroutine lsf_base_vjp_f1_rA

   !> Erroring default of the jet-contracted radius vector-Jacobian product
   !>
   !> The radius twin of [[lsf_base_vjp_f1_rA]]: the same per-point adjoint
   !> weights contract the same spatial (jet) indices, but the surviving
   !> derivative is taken with respect to the atomic radii
   !>
   !> - a radius is a scalar parameter, so the surviving index is just the atom
   !> - there is no Cartesian slot, and the result is one rank lower than the
   !>   nuclear one
   !>
   !> The result an implementation must produce is
   !>
   !>     res(i) = w0 * lsf1_rad(i)
   !>            + sum_a w1(a) * lsf2_r_rad(a, i)
   !>            + sum_{a<=b} m_ab * w2(ab) * lsf3_rr_rad(ab, i)
   !>
   !> for `i = 1 .. active_count()`, with the same compact active-index
   !> convention as `f3_rr_rad`, slot `i` belonging to atom `active_atom(i)`
   !>
   !> - entries beyond `active_count()` are left untouched, and nothing is written
   !>   when `active_count()` is zero
   !> - `w2` is the packed symmetric seed (6): xx, xy, yy, xz, yz, zz
   !> - off-diagonals are not doubled; the contraction applies `m_ab` = 1 on the
   !>   diagonal and 2 off it
   !> - `pack_sym_part` of [[moist_math_packed_sym]] packs a general seed
   !> - an implementation needs whatever prepared order its own `f3_rr_rad` needs,
   !>   2 for both SvdW and CFC today
   !> - must call `require_deriv(<that order>, "vjp_f1_rad")`
   !>
   !> Unlike `vjp_f1_rA`, which the isodensity LSFs override with an identically
   !> zero row, this default is *not* overridden there:
   !>
   !> - those LSFs do not see the radii at all, clear `radius_dependent`, and
   !>   leave the whole radius family erroring
   !> - a consumer gates on `radius_dependent`
   !> - one that reaches here anyway has asked an LSF for a derivative it
   !>   structurally does not have, and must abort rather than receive a zero it
   !>   cannot distinguish from a real answer
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Adjoint weights of the spatial Hessian, packed (6)
   !> @param[out] res  Contracted radius gradient [>= n_active]
   subroutine lsf_base_vjp_f1_rad(self, w0, w1, w2, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: w0
      real(wp), intent(in) :: w1(3)
      real(wp), intent(in) :: w2(:)
      real(wp), intent(out) :: res(:)

      res = 0.0_wp*w0*size(w1)*size(w2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "vjp_f1_rad derivative"
   end subroutine lsf_base_vjp_f1_rad

   !> Batched jet-contracted nuclear vector-Jacobian product over a tile of points
   !>
   !> Per column `p` the same contraction as [[lsf_base_vjp_f1_rA]] after
   !> `prepare(points(:, p))`, written against the tile's shared candidate list
   !>
   !> Contract shared by the default and every override:
   !>   - `points(3, np)`, `w0(np)`, `w1(3, np)`, `w2(6, np)` contiguous; columns
   !>     `1:np` are read, outputs are filled for columns `1:np` only
   !>   - `w2` packs the *symmetric part* of the `(3, 3)` seed in the order xx, xy, yy,
   !>     xz, yz, zz, off-diagonal entries not doubled; the contraction applies the
   !>     factor 2 on them, so the column equals `vjp_f1_rA` called with the same
   !>     packed seed (`lsf2_rr` is symmetric: only the symmetric part of a seed
   !>     contributes), and `pack_sym_part` from [[moist_math_packed_sym]] packs a
   !>     general seed
   !>   - `candidate_indices` as in `evaluate_batch`: ids in the LSF's `candidate_space`,
   !>     unique and in 1..ncenters (`error stop` otherwise), sorted ascending before
   !>     use, so the order of the list moves no bits; absent means the ids
   !>     1..ncenters, i.e. every atom, and the full screen per point
   !>   - `res(3, ncand, np)` with `ncand = size(candidate_indices)`, or `ncenters`
   !>     when absent; `res(:, k, p)` is the seed-contracted nuclear derivative at
   !>     point `p` with respect to the atom of candidate `candidate_indices(k)` (id
   !>     `k` when absent), in the list's own order; exact zero for a candidate out
   !>     of reach of point `p`, never a skipped write
   !>   - optional `status(np)`: 0 ok, 1 unreached (no atom within reach), 2 failed
   !>     or not computed; every entry `1:np` is defined on every return; an
   !>     unreached column is NaN in `res(:, :, p)`
   !>   - with `status` present an unreached column never stops the tile: the other
   !>     columns are computed
   !>   - with `status` absent an unreached column is an `error` naming the column;
   !>     the columns behind it are unspecified
   !>   - `error` names the first failing column ("vjp_f1_rA_batch failed at point
   !>     <i>: ...", a `prepare` / `prepare_subset` error) and the call returns, with
   !>     or without `status`; the failed and later columns are unspecified, their
   !>     `status` 2, as in `evaluate_batch`
   !>   - a shape that cannot hold the tile is a caller bug: `error stop`, see
   !>     [[lsf_vjp_batch_shape_message]]
   !>   - the single-point accessor cache is unspecified after this call
   !>
   !> Looping default:
   !>   - `prepare` / `prepare_subset` (on one ascending copy of the list), then the
   !>     per-point `vjp_f1_rA` with the packed seed, per column
   !>   - the active slots scatter into the candidate positions through
   !>     `active_atom` and `cand_to_user`; the rest of the column is zero
   !>   - needs the prepared order the per-point `vjp_f1_rA` needs (2 for SvdW, 3 for
   !>     CFC), a caller raises `set_max_deriv` first; the per-point abort stands
   !>   - `status` starts at 2 for every column, so a `prepare` error leaves the
   !>     failed column and those behind it at 2
   !>   - SSES runs this default, one `prepare` per column; both isodensity LSFs
   !>     override it with evaluation-free zero rows
   !>
   !> @param[inout] self              LSF instance, `prepare` mutates its cache
   !> @param[in]    points            Evaluation points (3, np)
   !> @param[in]    w0                Adjoint weights of the value (np)
   !> @param[in]    w1                Adjoint weights of the gradient (3, np)
   !> @param[in]    w2                Packed symmetric part of the Hessian seed (6, np)
   !> @param[out]   res               Contracted nuclear gradient (3, ncand, np)
   !> @param[out]   status            Per-column status (np), 0 ok, 1 unreached, 2 failed
   !>                                 or not computed
   !> @param[out]   error             First failing column
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   subroutine lsf_base_vjp_f1_rA_batch(self, points, w0, w1, w2, res, status, error, &
                                       candidate_indices)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      real(wp), intent(in), contiguous :: points(:, :)
      real(wp), intent(in), contiguous :: w0(:)
      real(wp), intent(in), contiguous :: w1(:, :)
      real(wp), intent(in), contiguous :: w2(:, :)
      real(wp), intent(out), contiguous :: res(:, :, :)
      integer, intent(out), optional, contiguous :: status(:)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: candidate_indices(:)

      !> Fill of an unreached column
      real(wp) :: nan
      !> Per-point result over the active slots
      real(wp), allocatable :: work(:, :)
      !> Ascending copy of `candidate_indices(1:nc)`
      integer, allocatable :: sorted(:)
      !> Candidate position of user-space atom id, 0 when not a candidate [ncenters]
      integer, allocatable :: position(:)
      type(error_type), allocatable :: column_error
      !> Column, point count, list size, slot, candidate position, escape status
      integer :: p, np, ncand, nc, i, k, pstat
      character(len=16) :: col
      character(len=:), allocatable :: msg

      np = size(points, 2)
      ncand = self%ncenters
      if (present(candidate_indices)) ncand = size(candidate_indices)
      if (present(status)) then
         msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand, &
                                           status)
      else
         msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand)
      end if
      if (len(msg) > 0) error stop msg
      if (present(status)) status(1:np) = 2

      ! Candidate position by user-space id; ids validated and sorted as evaluate_batch
      allocate (position(max(self%ncenters, 1)), work(ndim, max(self%ncenters, 1)))
      position = 0
      if (present(candidate_indices)) then
         allocate (sorted(ncand))
         call lsf_sorted_candidates(candidate_indices, sorted, nc, self%ncenters)
         do k = 1, ncand
            position(self%cand_to_user(candidate_indices(k))) = k
         end do
      else
         do k = 1, ncand
            position(self%cand_to_user(k)) = k
         end do
      end if
      nan = ieee_value(0.0_wp, ieee_quiet_nan)

      do p = 1, np
         if (present(candidate_indices)) then
            call self%prepare_subset(points(:, p), sorted(1:nc), column_error, pstat)
         else
            call self%prepare(points(:, p), column_error, pstat)
         end if
         if (allocated(column_error)) then
            write (col, "(i0)") p
            call move_alloc(column_error, error)
            error%message = "vjp_f1_rA_batch failed at point "//trim(col)//": "// &
                            error%message
            return
         end if
         if (pstat /= 0) then
            if (.not. present(status)) then
               write (col, "(i0)") p
               call fatal_error(error, "vjp_f1_rA_batch: point "//trim(col)// &
                                " is unreached, no atom within reach ("//self%name//")")
               return
            end if
            status(p) = 1
            res(:, :, p) = nan
            cycle
         end if

         call self%vjp_f1_rA(w0(p), w1(:, p), w2(:, p), work)
         res(:, :, p) = 0.0_wp
         do i = 1, self%active_count()
            k = position(self%active_atom(i))
            if (k == 0) then
               error stop "moist DROP LSF vjp_f1_rA_batch: active atom outside the list"
            end if
            res(:, k, p) = work(:, i)
         end do
         if (present(status)) status(p) = 0
      end do
   end subroutine lsf_base_vjp_f1_rA_batch

   !> Describe a `vjp_f1_rA_batch` call whose arrays cannot hold the tile
   !>
   !> Shared by the looping default and every override, so the aborts read the
   !> same everywhere; the caller `error stop`s on a nonempty message (`error stop`
   !> is not testable in process, the message is)
   !>
   !> - `points` with a leading extent other than 3
   !> - `w0` / `w1` / `w2` with a leading extent other than 1 / 3 / 6, or fewer
   !>   than `np` columns
   !> - `res` with extents other than `(3, ncand, >= np)`
   !> - `status` present with fewer than `np` entries
   !>
   !> @param[in] who    LSF name for the diagnostic
   !> @param[in] points Evaluation points (3, np); only the extents are read
   !> @param[in] w0     Value weights (np); only the extent is read
   !> @param[in] w1     Gradient weights (3, np); only the extents are read
   !> @param[in] w2     Packed Hessian seeds (6, np); only the extents are read
   !> @param[in] res    Result (3, ncand, np); only the extents are read
   !> @param[in] ncand  Number of candidates, `size(candidate_indices)` or `ncenters`
   !> @param[in] status Per-column status (optional); only the extent is read
   !> @returns   msg   Diagnostic, empty when the shapes fit
   pure function lsf_vjp_batch_shape_message(who, points, w0, w1, w2, res, ncand, &
                                             status) result(msg)
      !> LSF name for the diagnostic
      character(len=*), intent(in) :: who
      !> Evaluation points; only the extents are read
      real(wp), intent(in) :: points(:, :)
      !> Value weights; only the extent is read
      real(wp), intent(in) :: w0(:)
      !> Gradient weights; only the extents are read
      real(wp), intent(in) :: w1(:, :)
      !> Packed Hessian seeds; only the extents are read
      real(wp), intent(in) :: w2(:, :)
      !> Result; only the extents are read
      real(wp), intent(in) :: res(:, :, :)
      !> Number of candidates
      integer, intent(in) :: ncand
      !> Per-column status; only the extent is read
      integer, intent(in), optional :: status(:)
      !> Diagnostic, empty when the shapes fit
      character(len=:), allocatable :: msg

      character(len=*), parameter :: pre = "moist DROP LSF vjp_f1_rA_batch: "
      integer :: np

      msg = ""
      if (size(points, 1) /= ndim) then
         msg = pre//"points must be (3, np) ("//who//")"
         return
      end if
      np = size(points, 2)
      if (size(w0) < np) then
         msg = pre//"w0 holds fewer than np entries ("//who//")"
      else if (size(w1, 1) /= ndim .or. size(w1, 2) < np) then
         msg = pre//"w1 must be (3, >= np) ("//who//")"
      else if (size(w2, 1) /= 6 .or. size(w2, 2) < np) then
         msg = pre//"w2 must be packed (6, >= np) ("//who//")"
      else if (size(res, 1) /= ndim .or. size(res, 2) /= ncand &
               .or. size(res, 3) < np) then
         msg = pre//"res must be (3, ncand, >= np) ("//who//")"
      end if
      if (len(msg) > 0) return
      if (present(status)) then
         if (size(status) < np) then
            msg = pre//"status holds fewer than np entries ("//who//")"
         end if
      end if
   end function lsf_vjp_batch_shape_message

   !* ================================================================================= *!
   !*                            Geometric bound queries                                *!
   !* ================================================================================= *!

   !> Radius of a ball around the evaluation point provably free of surface
   !>
   !> Returns `r` such that `S` has no zero in `B(x, r)`, given `S(x)`
   !>
   !> - the generic construction is `r = |S(x)| / L`, with `L` a bound on the
   !>   gradient norm over the region, since `S` cannot travel from `S(x)` to zero
   !>   in less than that distance
   !> - concrete LSFs override with their own `L` (SvdW has the exact `L = 1`)
   !>
   !> What lets the octree surface coverage (projection levels 3-4) discard boxes
   !> within the search window: minima live on the surface, so a region proven
   !> free of surface holds no minimum
   !>
   !> Zero is the safe answer and the base returns it unconditionally: it
   !> certifies nothing, so a search built on it prunes nothing and runs into
   !> its own box budget rather than silently reporting a false certificate
   !>
   !> @param[in] self  LSF instance
   !> @param[in] lsf0  LSF value at the evaluation point
   !> @returns   r     Surface-free radius (zero when uncertified)
   pure function lsf_base_exclusion_radius(self, lsf0) result(r)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: lsf0
      real(wp) :: r

      r = 0.0_wp
      if (self%ncenters < 0 .and. lsf0 /= 0.0_wp) r = 0.0_wp
   end function lsf_base_exclusion_radius

   !> Surface-free radius around every unreached point
   !>
   !> Valid for each point no atom reaches: the per-point path's empty active
   !> list, a batch column with `status = 1`, i.e. every point whose nearest
   !> atom surface lies at least the reach `R + screening_offset(R)` away
   !>
   !> - a constant of the model and geometry, no point argument; the caller
   !>   applies it only to unreached points
   !> - zero (certify nothing) by default and whenever no such bound exists
   !> - SvdW and CFC with the fade on override it (band end minus the nearest-atom
   !>   bound of any surface point)
   !>
   !> @param[in] self  LSF instance
   !> @returns   r     Surface-free radius around an unreached point (zero when uncertified)
   pure function lsf_base_far_exclusion_radius(self) result(r)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Surface-free radius
      real(wp) :: r

      r = 0.0_wp
      if (self%ncenters < 0) r = 0.0_wp
   end function lsf_base_far_exclusion_radius

   !> Keep boxes when no stationarity bound is available
   !>
   !> @param[in] self       LSF instance
   !> @param[in] anchor     Distance origin (3)
   !> @param[in] center     Box midpoint (3)
   !> @param[in] half_width Nonnegative cube half-width
   pure function lsf_base_stationary_box_excluded(self, anchor, center, half_width) result(excluded)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Distance origin
      real(wp), intent(in) :: anchor(3)
      !> Box midpoint
      real(wp), intent(in) :: center(3)
      !> Cube half-width
      real(wp), intent(in) :: half_width
      !> No exclusion without a concrete LSF bound
      logical :: excluded

      excluded = .false.
   end function lsf_base_stationary_box_excluded

   !> Is this level set a function of a density, and with which constant factor?
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] factor  Chain-rule factor `dS/drho`; zero when unavailable
   !> @returns            Whether this LSF is backed by a density
   function lsf_base_density_adjoint_factor(self, factor) result(available)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: factor
      logical :: available

      factor = 0.0_wp
      available = .false.
   end function lsf_base_density_adjoint_factor

   !> Default: no atomic partition
   !>
   !> @param[in] self LSF instance
   pure function lsf_base_has_atom_partition(self) result(available)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Whether `atom_partition` is implemented
      logical :: available

      available = .false.
   end function lsf_base_has_atom_partition

   !> Erroring default of the atomic partition
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    scale Weight of the point
   !> @param[inout] acc   Per-atom accumulator, user space [ncenters]
   subroutine lsf_base_atom_partition(self, scale, acc)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      real(wp), intent(in) :: scale
      real(wp), intent(inout) :: acc(:)

      acc = acc + 0.0_wp*scale*self%ncenters
      error stop "moist DROP LSF: Chosen level set has no atomic partition"
   end subroutine lsf_base_atom_partition

   !> Default declaration: the level set reads nothing from the host
   !>
   !> @param[in]    self     LSF instance
   !> @param[inout] coupling Coupling being declared, scope set by the cavity
   !> @param[out]   error    Error handling
   subroutine lsf_base_declare_coupling(self, coupling, error)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Coupling being declared
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine lsf_base_declare_coupling

   !> Default: the level set needs no host state beyond the structure
   !>
   !> @param[in] self LSF instance
   pure function lsf_base_requires_build_inputs(self) result(requires)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Whether a build-phase answer is required
      logical :: requires

      requires = .false.
   end function lsf_base_requires_build_inputs

   !> Default build-phase read: nothing to read
   !>
   !> @param[inout] self     LSF instance
   !> @param[in]    coupling Cavity-scoped view of the build phase
   !> @param[in]    mol      Structure of the cavity's latest update
   !> @param[out]   error    Error handling
   subroutine lsf_base_read_build(self, coupling, mol, error)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self
      !> Cavity-scoped view of the build phase
      class(coupling_view_type), intent(in) :: coupling
      !> Structure of the cavity's latest update
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   end subroutine lsf_base_read_build

   !> Reject projected surface points this LSF cannot represent
   !>
   !> Called once per cavity update with the cavity's points, which the cavity
   !> marks as converged surface projections or not; the base accepts every point
   !>
   !> TODO: this has to be combined with the disconnected cav search at some point
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  xyz       Cavity points in Bohr (3, n)
   !> @param[in]  projected Whether each point is a converged surface projection (n)
   !> @param[out] error     Set when a projected point is rejected
   subroutine lsf_base_check_surface(self, xyz, projected, error)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: xyz(:, :)
      logical, intent(in) :: projected(:)
      type(error_type), allocatable, intent(out) :: error

      if (self%ncenters < 0 .and. size(xyz) < 0 .and. any(projected)) return
   end subroutine lsf_base_check_surface

   !* ================================================================================= *!
   !*                                    Screening                                      *!
   !* ================================================================================= *!

   !> Reject the candidates whose contribution is below `screening_threshold`
   !>
   !> The whole reject path is one squared-distance compare against the cached
   !> screening record: no `exp`, no `sqrt`, no allocation. `active_cand` is a
   !> caller-owned buffer (size >= size(candidate_indices)) so the hot path stays
   !> allocation-free
   !>
   !> @param[in]  self              LSF instance (must have been `update`d)
   !> @param[in]  point             Evaluation point [ndim]
   !> @param[in]  candidate_indices Candidate ids in the declared candidate space
   !> @param[out] active_cand       Surviving candidate ids, first `n_active` entries
   !> @param[out] n_active          Number of survivors
   pure subroutine lsf_base_screen_candidates(self, point, candidate_indices, &
                                              active_cand, n_active)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: point(ndim)
      integer, intent(in) :: candidate_indices(:)
      integer, intent(out) :: active_cand(:)
      integer, intent(out) :: n_active

      n_active = 0
      if (.not. allocated(self%cand_screen)) return
      call lsf_screen_records(self%cand_screen, point, candidate_indices, active_cand, n_active)
   end subroutine lsf_base_screen_candidates

   !> The gate of [[lsf_base_screen_candidates]] on the screening records alone
   !>
   !> - for a concrete LSF whose survivor buffer is a component of itself: the
   !>   records and the buffer are disjoint actuals, never `self` plus a part
   !>   of it
   !>
   !> @param[in]  screen            Screening records `cand_screen` [4, ncenters]
   !> @param[in]  point             Evaluation point [ndim]
   !> @param[in]  candidate_indices Candidate ids in the declared candidate space
   !> @param[out] active_cand       Surviving candidate ids, first `n_active` entries
   !> @param[out] n_active          Number of survivors
   pure subroutine lsf_screen_records(screen, point, candidate_indices, active_cand, n_active)
      !> Screening records
      real(wp), intent(in) :: screen(:, :)
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Candidate ids
      integer, intent(in) :: candidate_indices(:)
      !> Surviving candidate ids
      integer, intent(out) :: active_cand(:)
      !> Number of survivors
      integer, intent(out) :: n_active

      integer :: i, c
      real(wp) :: dx, dy, dz

      n_active = 0
      do i = 1, size(candidate_indices)
         c = candidate_indices(i)
         dx = point(1) - screen(1, c)
         dy = point(2) - screen(2, c)
         dz = point(3) - screen(3, c)
         if (dx*dx + dy*dy + dz*dz > screen(lsf_cand_bound_sq, c)) cycle
         n_active = n_active + 1
         active_cand(n_active) = c
      end do
   end subroutine lsf_screen_records

   !> Candidate ids whose screening reach sphere meets the box `[lo, hi]`
   !>
   !> Superset of the per-point gate for every point in the box, so feeding the
   !> result to [[lsf_base_screen_candidates]] reproduces the full-scan active set
   !>
   !> - ascending in candidate space, ready for `prepare_subset`
   !> - empty before `update`
   !> - strict superset: the box test widens the squared reach by a few ulp, so
   !>   FMA contraction in the gate cannot drop an active candidate; the gate is
   !>   the exact filter
   !> - error stop on a non-finite corner
   !>
   !> @param[in]  self  LSF instance
   !> @param[in]  lo    box corner [ndim]
   !> @param[in]  hi    opposite box corner [ndim]
   !> @param[out] cand  candidate ids in `cand(1:n)`, `size(cand) >= ncenters`
   !> @param[out] n     number of candidates
   pure subroutine lsf_base_candidates_in_box(self, lo, hi, cand, n)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Box corner
      real(wp), intent(in) :: lo(ndim)
      !> Opposite box corner
      real(wp), intent(in) :: hi(ndim)
      !> Candidate ids
      integer, intent(out) :: cand(:)
      !> Number of candidates
      integer, intent(out) :: n

      call self%cand_cells%query_box(lo, hi, cand, n)
   end subroutine lsf_base_candidates_in_box

   !> Candidate ids whose screening reach sphere meets the ball `(center, r)`
   !>
   !> Superset of the per-point gate for every point in the ball
   !>
   !> - ascending in candidate space, ready for `prepare_subset`
   !> - empty before `update`
   !> - strict superset as for [[lsf_base_candidates_in_box]]; the gate is the
   !>   exact filter
   !> - error stop on a non-finite center or radius
   !>
   !> @param[in]  self    LSF instance
   !> @param[in]  center  ball center [ndim]
   !> @param[in]  r       ball radius, `>= 0`
   !> @param[out] cand    candidate ids in `cand(1:n)`, `size(cand) >= ncenters`
   !> @param[out] n       number of candidates
   pure subroutine lsf_base_candidates_in_sphere(self, center, r, cand, n)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      !> Ball center
      real(wp), intent(in) :: center(ndim)
      !> Ball radius
      real(wp), intent(in) :: r
      !> Candidate ids
      integer, intent(out) :: cand(:)
      !> Number of candidates
      integer, intent(out) :: n

      call self%cand_cells%query_sphere(center, r, cand, n)
   end subroutine lsf_base_candidates_in_sphere

end module moist_cavity_drop_lsf_base
