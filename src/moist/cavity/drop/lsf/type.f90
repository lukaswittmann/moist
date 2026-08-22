!> Abstract level set function (LSF) base type for the DROP scheme
!>
!> The DROP cavity discretization is independent of *which* LSF defines the cavity surface.
!> This module provides the abstract base any concrete LSF must extend.
!>
!> A concrete LSF supplies:
!>   - its own constructor (per-LSF parameters; see e.g. SvdW's [[moist_cavity_drop_lsf_svdw_type]]%new)
!>   - a `candidate_space` declaration in that constructor (see below)
!>   - implementations of the deferred procedures below
!>   - optionally, overrides for the higher-order and direction-contracted
!>     derivatives the base declares with erroring defaults (`f4_*`, `f2_rArB`,
!>     `f3_r_rArB`, `f4_rr_rArB`, `f*_radrad`, `f*_rA_rad`, `normalized_f01_rA`,
!>     `tangent_*`, `hvp_*`). SvdW and CFC implement all of them; the isodensity
!>     LSFs implement none, so asking one of those for a fourth derivative aborts
!>     rather than returning a zero the caller cannot distinguish from a real
!>     answer
!>
!> Two geometric-bound queries
!> ---------------------------
!> Both of the questions below are of the form "how far can something be", both are
!> answered analytically from the LSF's own mathematics, and both exist so a caller
!> can skip work without knowing which concrete LSF it holds. They are otherwise
!> independent -- different subject, different consumer, different failure mode:
!>
!>   - [[exclusion_radius]]`(lsf0)` -- how far a piece of *surface* must be from the
!>     evaluation point, given the value there. Answers "can I discard this region?"
!>     for the certified branch search. Its safe answer is **zero** (certify
!>     nothing, prune nothing); an over-estimate would silently lose branches.
!>   - [[screening_offset]]`(radius)` -- how far an *atom* can be from the evaluation
!>     point and still contribute above `screening_threshold`. Answers "can I skip
!>     this atom?" for screening. Its safe answer is **huge** (screen nothing); an
!>     under-estimate would silently drop contributions.
!>
!> So the safe direction is opposite for the two, which is why neither can be
!> expressed in terms of the other, and why only the screening one has a documented
!> safety margin (see [[lsf_screening_reach_margin]]).
!>
!> Screening contract
!> ------------------
!> Every LSF answers exactly one screening question:
!>
!>    `screening_offset(radius)` = the radial offset, measured outward from the
!>    atom surface, at which this LSF's contribution from that atom equals
!>    `screening_threshold`.
!>
!> It is an *exact* analytic inverse, never a bound and never a search, and it is
!> the single source of truth for both screening consumers:
!>
!>   - the **per-point gate**: the base caches `(R + offset(R))**2` once per `update`,
!>     interleaved with the atom center in `cand_screen`, and the gate is a
!>     squared-distance compare against that cache -- no `exp` and no `sqrt`
!>     anywhere in the reject path (see [[lsf_base_screen_candidates]]).
!>   - the **cell-grid reach**: `neighbor_cutoff = lsf_screening_reach_margin * offset`,
!>     handed to the cavity so it can size the molecular cell grid without knowing
!>     the concrete LSF. `neighbor_cutoff` is `non_overridable` precisely so the two
!>     consumers cannot drift apart again.
!>
!> Candidate index space
!> ---------------------
!> The cavity's cell grid stores candidate atom ids in user space. Some LSFs want
!> them relabelled into the base's spatially-sorted ordering (which turns the
!> per-candidate geometry gathers in the screen loop into near-sequential reads and
!> is worth roughly 15% of a projection update); others index per-atom data with the
!> caller's own ids and must not be relabelled.
!>
!> That choice is *declared*, not inferred: each concrete LSF sets `candidate_space`
!> to [[lsf_candidate_space_sorted]] or [[lsf_candidate_space_user]] in its
!> constructor, and the base performs (or skips) the relabelling accordingly in
!> [[lsf_base_remap_candidate_grid]]. `update` aborts on an undeclared space, so a
!> silent no-op can no longer masquerade as a deliberate choice.
!>
!> The cavity holds a `class(moist_cavity_drop_lsf_type), allocatable` model from which thread-local clones are made.
!> The [[lsf_thread_slot]] wrapper enables arrays of polymorphic LSF clones, since Fortran has no class(...), allocatable :: arr(:)
module moist_cavity_drop_lsf_base
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_math_sorter_counting_sort, only: counting_argsort
   implicit none ()
   private

   public :: moist_cavity_drop_lsf_type
   public :: lsf_thread_slot
   public :: lsf_base_update
   public :: lsf_candidate_space_undeclared
   public :: lsf_candidate_space_user
   public :: lsf_candidate_space_sorted
   public :: lsf_screening_reach_margin
   public :: lsf_cand_bound_sq

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Row range of the atom center inside a `cand_screen` record
   integer, parameter :: lsf_cand_center(2) = [1, 3]
   !> Row of the squared screening reach inside a `cand_screen` record
   integer, parameter :: lsf_cand_bound_sq = 4
   !> Number of rows in a `cand_screen` record
   integer, parameter :: lsf_cand_screen_rows = 4

   !> `candidate_space` sentinel: the concrete LSF never declared one. Reaching
   !> `update` in this state is a programming error and aborts.
   integer, parameter :: lsf_candidate_space_undeclared = 0
   !> `candidate_space`: candidate ids are user-space (caller) atom ids and are
   !> handed to `prepare_subset` untouched. Required by LSFs that index per-atom
   !> data (basis shells, ...) with the caller's numbering.
   integer, parameter :: lsf_candidate_space_user = 1
   !> `candidate_space`: candidate ids are relabelled into the base's
   !> spatially-sorted ordering, so `prepare_subset` can index `cand_screen` /
   !> `cand_radii` directly with near-sequential access.
   integer, parameter :: lsf_candidate_space_sorted = 2

   !> Safety margin applied to `screening_offset` when sizing the cavity's
   !> per-atom cell-grid reach.
   !>
   !> The gate is exact, so a reach of exactly `R + offset(R)` already admits every
   !> atom that can pass it for points inside the atom bounding box. The margin buys
   !> headroom for the two places where that argument is not airtight: the cell grid
   !> clamps out-of-box query points into the boundary cell, and cell membership is
   !> decided by a sphere/AABB test in a different floating-point expression than the
   !> gate itself. A larger reach only ever adds candidates the gate then rejects, so
   !> the margin can never change a result -- only the amount of work.
   !>
   !> This replaces SvdW's historical trick of bisecting for `0.1 * threshold`
   !> instead of `threshold`, which was the same margin expressed as an
   !> unnamed fudge on the target value: `-3*ln(0.1*thr)/k` exceeds
   !> `-3*ln(thr)/k` by `3*ln(10)/k`, i.e. by ~8% at the production threshold.
   real(wp), parameter :: lsf_screening_reach_margin = 1.10_wp

   !> Abstract LSF base
   type, abstract :: moist_cavity_drop_lsf_type
      !> Number of atomic centers
      integer :: ncenters = 0
      !> Molecular structure data
      type(structure_type) :: mol
      !> Atomic radii (per-atom, user-space)
      real(wp), allocatable :: radii(:)
      !> Screening threshold below which the LSF contribution is treated
      !> as zero. Owned by the cavity (couples to the projection tolerance);
      !> the LSF concrete reads this value to size its internal screening
      !> caches whenever `update` runs. Direct users (tests calling the
      !> concrete without a cavity) may set it before calling `update`
      real(wp) :: screening_threshold = 0.0_wp
      !> Highest spatial derivative order actually present in the per-point cache after
      !> the latest `prepare`/`prepare_subset`
      integer :: prepared_deriv = -1

      !> Index space the concrete LSF wants its `prepare_subset` candidate ids in.
      !> Set by the concrete's constructor to [[lsf_candidate_space_user]] or
      !> [[lsf_candidate_space_sorted]]; `update` aborts while undeclared.
      integer :: candidate_space = lsf_candidate_space_undeclared

      !> LSF dependency on atomic radii
      logical :: radius_dependent = .true.

      !* ------------------------- Spatial sort machinery -------------------------- *!

      !> Permutation mapping a user-space atom index i to its spatially-sorted
      !> position. Rebuilt by `update` / `set_centers`; size = ncenters.
      integer, allocatable :: orig_to_sorted(:)
      !> Inverse of `orig_to_sorted`: the user-space atom index stored at
      !> spatially-sorted position j. Size = ncenters.
      integer, allocatable :: sorted_to_orig(:)

      !* ------------------- Candidate-space geometry and screening ---------------- *!
      !> The arrays below are indexed by *candidate id in the declared
      !> `candidate_space`*, so a screen loop can consume `prepare_subset`'s ids
      !> without any per-candidate translation. For a `sorted` LSF they are the
      !> spatially-sorted mirror of the geometry; for a `user` LSF they are the
      !> user-space arrays themselves.

      !> Everything the per-point reject test reads, interleaved per candidate
      !> [4, ncenters]:
      !>   `cand_screen(1:3, c)` -- the atom center
      !>   `cand_screen(4, c)`   -- the squared screening reach `(R + offset(R))**2`,
      !>                            saturated at `huge` rather than overflowed when
      !>                            screening is disabled
      !> One 32-byte record, so rejecting a candidate touches a single cache line
      !> instead of striding through a separate centers array and a separate bounds
      !> array. Use the [[lsf_cand_center]] / [[lsf_cand_bound_sq]] index parameters
      !> rather than bare literals.
      real(wp), allocatable :: cand_screen(:, :)
      !> Atom radii in candidate space [ncenters]. Read only for candidates that
      !> pass the gate, hence kept out of `cand_screen`
      real(wp), allocatable :: cand_radii(:)
      !> User-space atom id of candidate c [ncenters]
      integer, allocatable :: cand_to_user(:)
      !> Candidate ids covering every atom, in an order that makes the resulting
      !> active list come out in ascending user-space order. Fed to the screen
      !> loop by `prepare` (the unscreened full scan) [ncenters]
      integer, allocatable :: full_scan_cand(:)
   contains
      !> Bind molecular geometry and rebuild the screening caches. Concrete LSFs
      !> may override to refresh additional caches; the override *must* call the
      !> base implementation first via `call lsf_base_update(self, mol, radii)`
      procedure :: update => lsf_base_update
      !> Move the atom centers without touching radii, parameters or derivative
      !> storage (used by nuclear finite-difference drivers)
      procedure, non_overridable :: set_centers => lsf_base_set_centers
      !> Abort when an accessor is asked for an order `prepare` did not compute
      procedure :: require_deriv => lsf_base_require_deriv
      !> Relabel the cavity's cell-grid candidate lists into the declared
      !> candidate index space. Not overridable: the declaration decides.
      procedure, non_overridable :: remap_candidate_grid => lsf_base_remap_candidate_grid
      !> Reject candidates whose contribution is below `screening_threshold`
      procedure, non_overridable :: screen_candidates => lsf_base_screen_candidates
      !> Change `screening_threshold` in place, rewriting only the cached reach
      procedure, non_overridable :: set_screening_threshold => &
         lsf_base_set_screening_threshold
      !> Geometric bound query 1 of 2: how far a *surface* can be from a point
      !> with a known LSF value (see the "Two geometric-bound queries" note above)
      procedure :: exclusion_radius => lsf_base_exclusion_radius
      !> Geometric bound query 2 of 2, consumer side: how far an *atom* can be
      !> from a point and still matter. Derived from `screening_offset`; not
      !> overridable, so gate and reach stay the same mathematics.
      procedure, non_overridable :: neighbor_cutoff => lsf_base_neighbor_cutoff
      !> Cache per-point screening / state
      procedure(lsf_prepare_iface), deferred :: prepare
      !> Cache per-point state with a caller-provided candidate list.
      procedure(lsf_prepare_subset_iface), deferred :: prepare_subset
      !> Configure the highest spatial derivative order required.
      procedure(lsf_set_max_deriv_iface), deferred :: set_max_deriv
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
      !> Exact radial offset from the atom surface at which this LSF's
      !> contribution equals `screening_threshold`
      procedure(lsf_screening_offset_iface), deferred :: screening_offset

      !* ------------------ Higher-order and contracted derivatives ---------------- *!

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
   end type moist_cavity_drop_lsf_type

   !> Wrapper struct for arrays of polymorphic-allocatable LSF clones.
   !>
   !> Fortran does (afaik) not allow `class(...), allocatable :: arr(:)`, so per-thread LSF copies are
   !> stored as `type(lsf_thread_slot) :: arr(:)` with the polymorphic clone living inside the wrapper
   type :: lsf_thread_slot
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
   end type lsf_thread_slot

   abstract interface

      !> @param[inout] self   LSF instance
      !> @param[in]    point  Evaluation point (3,)
      !> @param[out]   error  Evaluation failure at this point
      subroutine lsf_prepare_iface(self, point, error)
         import :: wp, error_type, moist_cavity_drop_lsf_type
         implicit none ()
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         real(wp), intent(in) :: point(3)
         type(error_type), allocatable, intent(out) :: error
      end subroutine lsf_prepare_iface

      !> @param[inout] self              LSF instance
      !> @param[in]    point             Evaluation point (3,)
      !> @param[in]    candidate_indices Atom ids to consider, always in the
      !>                                 LSF's declared `candidate_space`
      !> @param[out]   error             Evaluation failure at this point
      subroutine lsf_prepare_subset_iface(self, point, candidate_indices, error)
         import :: wp, error_type, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         real(wp), intent(in) :: point(3)
         integer, intent(in) :: candidate_indices(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine lsf_prepare_subset_iface

      !> @param[inout] self LSF instance
      !> @param[in]    n    Requested max derivative order
      subroutine lsf_set_max_deriv_iface(self, n)
         import :: moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(inout) :: self
         integer, intent(in) :: n
      end subroutine lsf_set_max_deriv_iface

      !> @param[in]  self LSF instance
      !> @returns         Number of currently active atoms
      pure function lsf_active_count_iface(self) result(n)
         import :: moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         integer :: n
      end function lsf_active_count_iface

      !> @param[in]  self LSF instance
      !> @param[in]  i    Active-list index (1 <= i <= active_count())
      !> @returns         User-space atom id of the i-th active atom
      pure function lsf_active_atom_iface(self, i) result(idx)
         import :: moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         integer, intent(in) :: i
         integer :: idx
      end function lsf_active_atom_iface

      !> @param[in]  self LSF instance
      !> @param[out] val  LSF value at the current evaluation point
      subroutine lsf_f0_iface(self, val)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out) :: val
      end subroutine lsf_f0_iface

      !> @param[in]  self     LSF instance
      !> @param[out] lsf0     LSF value (optional)
      !> @param[out] lsf1_r   Gradient w.r.t. spatial coords (optional)
      !> @param[out] lsf2_rr  Hessian w.r.t. spatial coords (optional)
      subroutine lsf_f012_r_iface(self, lsf0, lsf1_r, lsf2_rr)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf0
         real(wp), intent(out), optional :: lsf1_r(:)
         real(wp), intent(out), optional :: lsf2_rr(:, :)
      end subroutine lsf_f012_r_iface

      !> @param[in]  self     LSF instance
      !> @param[out] lsf0     LSF value (optional)
      !> @param[out] lsf1_r   Spatial gradient (optional)
      !> @param[out] lsf2_rr  Spatial Hessian (optional)
      !> @param[out] lsf3_rrr Third spatial derivative tensor [3, 3, 3]
      subroutine lsf_f3_rrr_iface(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf0
         real(wp), intent(out), optional :: lsf1_r(:)
         real(wp), intent(out), optional :: lsf2_rr(:, :)
         real(wp), intent(out) :: lsf3_rrr(:, :, :)
      end subroutine lsf_f3_rrr_iface

      !> Every nuclear index is an active-list index: slot `i` belongs to the
      !> atom `active_atom(i)`. The caller owns all three buffers and sizes their
      !> nuclear extent from `active_count()`.
      !>
      !> @param[in]  self       LSF instance
      !> @param[out] lsf1_rA    LSF gradient w.r.t. nuclear positions (optional)
      !> @param[out] lsf2_r_rA  Mixed second derivative (optional)
      !> @param[out] lsf3_rr_rA Mixed third derivative
      subroutine lsf_f3_rr_rA_iface(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         import :: wp, moist_cavity_drop_lsf_type
         class(moist_cavity_drop_lsf_type), intent(in) :: self
         real(wp), intent(out), optional :: lsf1_rA(:, :)
         real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
         real(wp), intent(out) :: lsf3_rr_rA(:, :, :, :)
      end subroutine lsf_f3_rr_rA_iface

      !> Exact radial offset, measured outward from the atom surface, at which
      !> this LSF's contribution from an atom of the given radius equals
      !> `screening_threshold`.
      !>
      !> Contract:
      !>   - analytic and exact -- the closed-form inverse of the LSF's own decay,
      !>     not a bound and not a numerical search
      !>   - `>= 0`
      !>   - `huge(0.0_wp)` when screening is disabled (`screening_threshold <= 0`)
      !>     or the parameters make the criterion ill-defined
      !>
      !> Both screening consumers are derived from it: the per-point gate caches
      !> `radius + screening_offset(radius)`, and the cavity cell-grid reach is
      !> `lsf_screening_reach_margin` times this value.
      !>
      !> @param[in] self    LSF instance
      !> @param[in] radius  Atom radius (Bohr)
      !> @returns           Radial offset from the atom surface (Bohr)
      pure function lsf_screening_offset_iface(self, radius) result(offset)
         import :: wp, moist_cavity_drop_lsf_type
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
   !> concrete per-atom storage untouched. Nuclear finite-difference drivers use
   !> this to displace one coordinate without re-initialising the whole LSF.
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
   !> Runs once per geometry change, never on the per-point path.
   !>
   !> @param[inout] self LSF instance
   subroutine lsf_base_rebuild_screening(self)
      class(moist_cavity_drop_lsf_type), intent(inout) :: self

      !> Loop indices and atom count
      integer :: i, j, n
      !> Screening reach of one atom and the largest bound whose square is finite
      real(wp) :: bound, bound_sq_max

      if (self%candidate_space /= lsf_candidate_space_user .and. &
          self%candidate_space /= lsf_candidate_space_sorted) then
         error stop "moist DROP LSF: candidate_space was never declared -- the "// &
            "concrete LSF's constructor must set it to lsf_candidate_space_user "// &
            "or lsf_candidate_space_sorted"
      end if

      n = self%ncenters

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
         ! sums in the same order regardless of how the spatial sort came out.
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

   !> Recompute the cached squared screening reach of every candidate
   !>
   !> Split out of [[lsf_base_rebuild_screening]] because it is the only part
   !> that depends on `screening_threshold`: it touches no allocation and does
   !> not re-sort, so [[lsf_base_set_screening_threshold]] can re-run it alone.
   !>
   !> @param[inout] self LSF instance with the candidate mirror already built
   subroutine lsf_base_refresh_bounds(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: self

      !> Screening reach of one atom and the largest bound whose square is finite
      real(wp) :: bound, bound_sq_max
      integer :: j

      if (.not. allocated(self%cand_screen)) return

      ! Largest bound whose square still fits; beyond it the gate saturates so
      ! that "screening disabled" never turns into an overflow.
      bound_sq_max = sqrt(huge(0.0_wp))
      do j = 1, size(self%cand_radii)
         bound = self%cand_radii(j) + self%screening_offset(self%cand_radii(j))
         if (bound >= bound_sq_max) then
            self%cand_screen(lsf_cand_bound_sq, j) = huge(0.0_wp)
         else
            self%cand_screen(lsf_cand_bound_sq, j) = bound*bound
         end if
      end do
   end subroutine lsf_base_refresh_bounds

   !> Change the screening threshold in place
   !>
   !> Only the cached per-candidate reach depends on the threshold, so this
   !> rewrites one row of `cand_screen` and leaves the spatial sort, the
   !> candidate mirror and every allocation untouched. That makes it cheap
   !> enough to bracket a solve with -- O(ncenters) stores, no allocation.
   !>
   !> The cell grid is *not* rebuilt. Its reach is derived from the threshold
   !> the grid was built with, so this is safe only for a threshold that is
   !> *looser* than that one (a smaller reach keeps the grid a valid superset).
   !> Tightening the threshold below the grid's own would silently miss atoms.
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

   !> Coarse spatial bucket sort of the atom centers
   !>
   !> Assigns each atom to one of 6^3 = 216 coarse 3D buckets and builds the
   !> permutation with a single-pass counting sort (O(N), stable), inspired by
   !> stdlib's int8 radix_sort. Sorting the geometry this way turns the gathers
   !> driven by cell-local candidate lists into near-sequential reads.
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

      if (self%prepared_deriv < 0) return
      if (order <= self%prepared_deriv) return

      write (want, "(i0)") order
      write (got, "(i0)") self%prepared_deriv
      error stop "moist DROP LSF: "//caller//" needs derivative order "// &
         trim(want)//" but prepare ran at "//trim(got)// &
         " -- raise set_max_deriv before prepare"
   end subroutine lsf_base_require_deriv

   !* ================================================================================= *!
   !*                Erroring defaults of the optional derivative set                   *!
   !* ================================================================================= *!
   !
   ! Returning zeros here would be a silent wrong answer -- the caller cannot tell
   ! a genuine vanishing tensor from a missing implementation. Each default
   ! therefore aborts with its own `error stop`, spelled out in full so the
   ! diagnostic naming the missing derivative sits at the line that raises it.
   ! Each result is still zeroed on the line before the abort, with the dummy
   ! arguments folded into that store: neither is ever observable (the next
   ! statement terminates the program) and both exist only so that
   ! `-Wunused-dummy-argument` does not have to be switched off here.

   !> Erroring default of the pure spatial fourth derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Fourth spatial derivative [3, 3, 3, 3]
   subroutine lsf_base_f4_rrrr(self, lsf4_rrrr)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rrrr(:, :, :, :)

      lsf4_rrrr = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rrrr derivative"
   end subroutine lsf_base_f4_rrrr

   !> Erroring default of the mixed fourth derivative d^4S / (dr^3 dR_A)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf4_rrr_rA Mixed fourth derivative [3, 3, 3, 3, n_active]
   subroutine lsf_base_f4_rrr_rA(self, lsf4_rrr_rA)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rrr_rA(:, :, :, :, :)

      lsf4_rrr_rA = 0.0_wp*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "f4_rrr_rA derivative"
   end subroutine lsf_base_f4_rrr_rA

   !> Erroring default of the radius derivative d^3S / (dr^2 dR_a)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf1_rad    dS/dR_a [>= n_active] (optional)
   !> @param[out] lsf2_r_rad  d^2S/(dr dR_a) [3, >= n_active] (optional)
   !> @param[out] lsf3_rr_rad d^3S/(dr^2 dR_a) [3, 3, >= n_active]
   subroutine lsf_base_f3_rr_rad(self, lsf1_rad, lsf2_r_rad, lsf3_rr_rad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out), optional :: lsf1_rad(:)
      real(wp), intent(out), optional :: lsf2_r_rad(:, :)
      real(wp), intent(out) :: lsf3_rr_rad(:, :, :)

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
   !> @param[out] lsf4_rr_rArB Mixed fourth derivative [3, 3, 3, n_act, 3, n_act]
   subroutine lsf_base_f4_rr_rArB(self, lsf4_rr_rArB)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rr_rArB(:, :, :, :, :, :)

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
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_base_f4_rr_radrad(self, lsf4_rr_radrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rr_radrad(:, :, :, :)

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
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b) [3, 3, 3, >= n_act, >= n_act]
   subroutine lsf_base_f4_rr_rA_rad(self, lsf4_rr_rA_rad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(out) :: lsf4_rr_rA_rad(:, :, :, :, :)

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
   !> @param[out] res  Contracted Hessian derivative [3, 3]
   subroutine lsf_base_tangent_f2_rr(self, v, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:, :)

      res = 0.0_wp*size(v, 2)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "tangent_f2_rr derivative"
   end subroutine lsf_base_tangent_f2_rr

   !> Erroring default of the contracted third-derivative derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  Contracted third-derivative derivative [3, 3, 3]
   subroutine lsf_base_tangent_f3_rrr(self, v, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:, :, :)

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
   !> @param[out] res  Contracted mixed fourth derivative [3, 3, 3, n_active]
   !> @param[in]  vrad Radius directions [ncenters] (optional)
   subroutine lsf_base_hvp_f3_rr_rA(self, v, res, vrad)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(out) :: res(:, :, :, :)
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
   !> @param[out] res  Contracted mixed fourth derivative [3, 3, n_active]
   subroutine lsf_base_hvp_f3_rr_rad(self, v, vrad, res)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: v(:, :)
      real(wp), intent(in) :: vrad(:)
      real(wp), intent(out) :: res(:, :, :)

      res = 0.0_wp*size(v, 2)*size(vrad)*self%ncenters
      error stop "moist DROP LSF: Chosen level set does not support "// &
         "hvp_f3_rr_rad derivative"
   end subroutine lsf_base_hvp_f3_rr_rad

   !* ================================================================================= *!
   !*                            Geometric bound queries                                *!
   !* ================================================================================= *!

   !> Radius of a ball around the evaluation point provably free of surface
   !>
   !> Returns `r` such that `S` has no zero in `B(x, r)`, given `S(x)`. The
   !> generic construction is `r = |S(x)| / L` with `L` a bound on the gradient
   !> norm over the region: `S` cannot travel from `S(x)` to zero in less than
   !> that distance. Concrete LSFs override with their own `L` (SvdW has the
   !> exact `L = 1`).
   !>
   !> This is what lets the certified branch search (`proj_level = 9`)
   !> enumerate *all* branches carrying weight: minima live on the surface, so
   !> a region proven free of surface holds no minimum.
   !>
   !> Zero is the safe answer and the base returns it unconditionally: it
   !> certifies nothing, so a search built on it prunes nothing and runs into
   !> its own box budget rather than silently reporting a false certificate.
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

   !> Per-atom cell-grid reach: the screening offset plus an explicit safety margin
   !>
   !> The cavity adds the atom radius and hands the result to the molecular cell
   !> grid, so a candidate list always contains every atom the per-point gate could
   !> still accept. Derived from the same `screening_offset` the gate uses, which is
   !> why this procedure is `non_overridable`.
   !>
   !> @param[in] self    LSF instance
   !> @param[in] radius  Atom radius (Bohr)
   !> @returns           Radial offset from the atom surface (Bohr)
   pure function lsf_base_neighbor_cutoff(self, radius) result(d)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      real(wp), intent(in) :: radius
      real(wp) :: d

      d = self%screening_offset(radius)
      if (d >= huge(0.0_wp)/lsf_screening_reach_margin) then
         ! Screening disabled: keep the reach unbounded instead of overflowing
         d = huge(0.0_wp)
      else
         d = lsf_screening_reach_margin*d
      end if
   end function lsf_base_neighbor_cutoff

   !* ================================================================================= *!
   !*                                    Screening                                      *!
   !* ================================================================================= *!

   !> Relabel the cavity's cell-grid candidate atom ids into the declared space
   !>
   !> The cavity always builds the grid with user-space ids; a `sorted` LSF gets
   !> them rewritten once here, per projector, so that its per-point screen loop
   !> pays nothing per candidate.
   !>
   !> @param[in]    self      LSF instance (must have been `update`d)
   !> @param[inout] cell_nlat Flat cell-grid candidate atom-id list
   subroutine lsf_base_remap_candidate_grid(self, cell_nlat)
      class(moist_cavity_drop_lsf_type), intent(in) :: self
      integer, intent(inout) :: cell_nlat(:)

      integer :: i

      if (self%candidate_space /= lsf_candidate_space_sorted) return
      if (.not. allocated(self%orig_to_sorted)) return
      do i = 1, size(cell_nlat)
         cell_nlat(i) = self%orig_to_sorted(cell_nlat(i))
      end do
   end subroutine lsf_base_remap_candidate_grid

   !> Reject the candidates whose contribution is below `screening_threshold`
   !>
   !> The whole reject path is one squared-distance compare against the cached
   !> screening record: no `exp`, no `sqrt`, no allocation. `active_cand` is a
   !> caller-owned buffer (size >= size(candidate_indices)) so the hot path stays
   !> allocation-free.
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

      integer :: i, c
      real(wp) :: dx, dy, dz

      n_active = 0
      if (.not. allocated(self%cand_screen)) return

      do i = 1, size(candidate_indices)
         c = candidate_indices(i)
         dx = point(1) - self%cand_screen(1, c)
         dy = point(2) - self%cand_screen(2, c)
         dz = point(3) - self%cand_screen(3, c)
         if (dx*dx + dy*dy + dz*dz > self%cand_screen(lsf_cand_bound_sq, c)) cycle
         n_active = n_active + 1
         active_cand(n_active) = c
      end do
   end subroutine lsf_base_screen_candidates

end module moist_cavity_drop_lsf_base
