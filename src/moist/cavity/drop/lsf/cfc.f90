!> COSMO Fine Cavity (CFC) level set function
!>
!> The mathematics lives in the code-generated module
!> [[moist_cavity_drop_lsf_cfc_kernel]]; this module is the orchestration layer
!> between it and the LSF contract of [[moist_cavity_drop_lsf_base]]. It is the
!> CFC twin of [[moist_cavity_drop_lsf_svdw]] and follows the same division of
!> labour.
!>
!> The Diedenhofen-Klamt 2018 pseudo-density is
!>
!>    PD(r) = sum_a       exp{ a1 (s_a - 1) }
!>          + sum_{a<b} c (1 - vec(s_a).vec(s_b))^m exp{ a2 (s_a + s_b - 2) }
!>
!> with `s_a = ||r - r_a|| / R_a` and `vec(s_a) = (r - r_a) / R_a`, and the level
!> set this module returns is
!>
!>    S(r) = -log PD(r)
!>
!> matching the SvdW sign convention (interior negative). The kernel returns
!> derivatives of `S` itself, so nothing here negates anything.
!>
!> Why this one is O(n**2) and SvdW is not
!> ---------------------------------------
!> SvdW factorizes: its whole derivative set follows from power sums of a
!> per-atom screening factor, so an evaluation point costs O(n_active). CFC does
!> not. The `(1 - vec(s_a).vec(s_b))^m` factor couples the two centers of a pair,
!> there is no power-sum rewrite, and the pair sweep is irreducibly a double
!> loop. That is a property of the level set, not of this implementation. What
!> the kernel *does* avoid is O(n**2) *storage*: the direction-contracted
!> families (`tangent_*`, `hvp_*`) contract a nuclear index inside the kernel, so
!> the pairwise nuclear tensors are never formed unless a caller explicitly asks
!> for one of the uncontracted `f*_rArB` accessors -- whose result is itself
!> O(n**2), so there the cost is the answer.
!>
!> Division of labour
!> ------------------
!>   - `prepare` / `prepare_subset` run the base screening gate, cache the
!>     per-atom geometry of the survivors, and accumulate the two families that
!>     need no direction vector: the spatial pseudo-density tensors `pd*`
!>     (aggregate) and the one-nuclear-index tensors `qn*` (per active atom).
!>     Both sweeps are fused into the same atom / pair loops.
!>   - every accessor lifts those accumulators to derivatives of `S` on demand.
!>     The direction-contracted families (`tg*`, `hvp*`) and the two-nucleus
!>     family (`qq*`) depend on arguments `prepare` does not have, so those
!>     accessors run their own sweep.
!>
!> Buffers, never allocations
!> --------------------------
!> No accessor returns an `allocatable, intent(out)` result. Every result is a
!> caller-provided buffer whose nuclear extent the caller sizes from
!> [[active_count]]; the accessor writes its first `active_count()` slots and
!> leaves the rest alone. The per-point buffers are sized once per `update`.
!>
!> Index space of the nuclear outputs
!> ----------------------------------
!> Every nuclear index is an *active-list* index: slot `i` belongs to the atom
!> `active_atom(i)`. Screened-away atoms have no slot at all.
!>
!> Derivative-order contract
!> -------------------------
!> `max_deriv` is the highest *total* derivative order `prepare` provisions. It
!> buys the spatial family up to that order and the one-nuclear-index family up
!> to one spatial order less, so both stop at the same total order and a
!> spatial-only caller never pays for a nuclear ladder it will not read. Each
!> accessor therefore asks [[require_deriv]] for the highest *total* order of the
!> accumulator it reads, not for the order of the tensor it returns: `f3_rr_rA`
!> returns a rank-3 object but reads `qn2_rr`, which is a total order 3, so it
!> requires 3. Anything an accessor builds itself (`tg*`, `hvp*`, `qq*`, `rd*`,
!> `rh*`) does not enter that count. An accessor asked for an order `prepare`
!> did not provision aborts; none of them returns zeros, which a caller could
!> not tell from a real answer.
!>
!> Cavity radii
!> ------------
!> The radii are parameters of the level set, so a model that ties them to the
!> geometry differentiates through them as well. `f3_rr_rad` is that channel's
!> gradient ladder and `hvp_*_rad` its Hessian row, the radius counterparts of
!> `f3_rr_rA` and `hvp_*_rA`. Because a radius is a scalar parameter and not a
!> coordinate, the `_rad` results carry no trailing derivative index -- their
!> rank *is* the spatial order -- and the channel consumes no derivative order,
!> so a `_rad` accessor requires one order less than its `_rA` twin.
!>
!> The direction of the `hvp_*` accessors correspondingly becomes a pair
!> `(v, vrad)`: passing `vrad` to `hvp_*_rA` promotes its contraction from
!> `sum_B v_B . d/dR_B` to `sum_B (v_B . d/dR_B + vr_B d/dRad_B)`, and
!> `hvp_*_rad` is the row that retains a radius instead of a position. Between
!> them the two cover all four blocks of the joint Hessian. A fixed-radius
!> caller omits `vrad` and is unaffected, down to the bit -- the radius terms
!> live in their own kernel routines and its sweeps are untouched.
!>
!> Reverse mode
!> ------------
!> `vjp_f1_rA` is the adjoint of the `f3_rr_rA` ladder: instead of returning the
!> three mixed tensors it contracts their spatial indices against per-point
!> adjoint weights `(w0, w1, w2)` and returns the nuclear gradient row alone,
!> three numbers per atom instead of 3 + 9 + 27. It is a pure accessor like
!> `f3_rr_rA` -- same cached `pd*` and `qn*`, same prepared order, no sweep of
!> its own -- and the saving is entirely inside the generated kernel, where the
!> weights are folded in before the common subexpressions are taken.
!>
!> `vjp_f1_rad` is the same adjoint on the radius channel, the twin of
!> `vjp_f1_rA` exactly as `f3_rr_rad` is the twin of `f3_rr_rA`: one adjoint jet
!> in, and the `f1_rad` / `f2_r_rad` / `f3_rr_rad` ladder collapsed to a single
!> number per atom rather than 1 + 3 + 9. It carries the radius channel's two
!> deviations with it -- no trailing derivative index, and one prepared order
!> less than its `_rA` twin -- and, like `f3_rr_rad`, runs the `rd*` sweep
!> itself.
!>
!> `rd*` is deliberately not cached by `prepare` even though it depends on
!> nothing but the geometry: caching it would charge every fixed-radius caller
!> an extra O(n_active**2) pair sweep at every grid point for a family it never
!> reads.
!>
module moist_cavity_drop_lsf_cfc
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_sorted
   use moist_cavity_drop_lsf_cfc_param, only: moist_cavity_drop_lsf_cfc_param_type
   use moist_cavity_drop_lsf_cfc_kernel, only: &
      cfc_atomic_term_eval, cfc_pair_spatial_eval, &
      cfc_atomic_nuclear_eval, cfc_pair_nuclear_eval, &
      cfc_atomic_hessian_eval, cfc_pair_hessian_eval, &
      cfc_atomic_tangent_eval, cfc_pair_tangent_eval, &
      cfc_atomic_hvp_eval, cfc_pair_hvp_eval, &
      cfc_atomic_radius_eval, cfc_pair_radius_eval, &
      cfc_atomic_radius2_eval, cfc_pair_radius2_eval, &
      cfc_atomic_nucrad_eval, cfc_pair_nucrad_eval, &
      cfc_atomic_radius_hvp_eval, cfc_pair_radius_hvp_eval, &
      cfc_spatial_eval, cfc_nuclear_eval, cfc_hessian_eval, &
      cfc_tangent_eval, cfc_hvp_eval, cfc_vjp_eval, &
      cfc_radius_eval, cfc_radius2_eval, cfc_nucrad_eval, cfc_radius_hvp_eval, &
      cfc_radius_vjp_eval
   implicit none (type, external)
   private

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Smallest spatial-gradient norm [[lsf_normalized_f01_rA]] treats as
   !> non-degenerate. Below it the normalized level set has no defined value and
   !> the routine reports zero, matching `svdw_normalized_eval`
   real(wp), parameter :: normalized_eps = 1.0e-14_wp

   !> COSMO Fine Cavity LSF
   !>
   !> Concrete level set function: takes its four shape parameters through
   !> [[new]], caches the screened per-atom geometry plus the `pd*` / `qn*`
   !> accumulators via [[prepare]], and lifts every derivative tensor on demand
   !> from the generated kernel. Inherits the common atom-LSF state (ncenters,
   !> mol, radii, screening caches) from [[moist_cavity_drop_lsf_type]].
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_cfc_type
      !> CFC shape parameters (a1, a2, c, m)
      type(moist_cavity_drop_lsf_cfc_param_type) :: param

      !> Highest total derivative order `prepare` provisions (0..4)
      integer :: max_deriv = 2

      !* --------------------------- Per-point active list -------------------------- *!

      !> Number of atoms that survived screening at the cached point
      integer :: n_active = 0
      !> Candidate-space ids of the survivors, first `n_active` valid [ncenters]
      integer, allocatable :: active_cand(:)
      !> User-space atom id of active slot i [ncenters]
      integer, allocatable :: act_atom(:)
      !> Displacement `r - R_A` of active slot i [ndim, ncenters]
      real(wp), allocatable :: act_d(:, :)
      !> Radius of active slot i [ncenters]
      real(wp), allocatable :: act_radius(:)

      !* --------------------- Aggregate pseudo-density tensors ---------------------- *!
      !> `d^k PD / dr^k` summed over every active atom and pair. Fixed-size
      !> components, refreshed up to `max_deriv` by every `prepare`; the orders
      !> above it are stale but defined, and no lift branch reads them.

      !> Pseudo-density value
      real(wp) :: pd0 = 0.0_wp
      !> Pseudo-density spatial gradient
      real(wp) :: pd1_r(ndim) = 0.0_wp
      !> Pseudo-density spatial Hessian
      real(wp) :: pd2_rr(ndim, ndim) = 0.0_wp
      !> Pseudo-density third spatial derivative
      real(wp) :: pd3_rrr(ndim, ndim, ndim) = 0.0_wp
      !> Pseudo-density fourth spatial derivative
      real(wp) :: pd4_rrrr(ndim, ndim, ndim, ndim) = 0.0_wp

      !* ------------------ Per-atom one-nuclear-index tensors (qn) ------------------ *!
      !> `d/dR_A of d^k PD / dr^k`, summed over atom A's own term and every pair
      !> containing A. The trailing index is the active slot; the one before it
      !> is the retained nuclear component. Filled to spatial order
      !> `max_deriv - 1`, i.e. to the same *total* order as the `pd*` family.

      !> Order-0 nuclear tensor [ndim, ncenters]
      real(wp), allocatable :: qn0(:, :)
      !> Order-1 nuclear tensor [ndim, ndim, ncenters]
      real(wp), allocatable :: qn1_r(:, :, :)
      !> Order-2 nuclear tensor [ndim, ndim, ndim, ncenters]
      real(wp), allocatable :: qn2_rr(:, :, :, :)
      !> Order-3 nuclear tensor [ndim, ndim, ndim, ndim, ncenters]
      real(wp), allocatable :: qn3_rrr(:, :, :, :, :)
   contains
      !> Constructor: configure shape parameters and declare the candidate space
      procedure, public :: new => lsf_new
      !> Bind molecular geometry and resize the per-atom caches
      procedure, public :: update => lsf_update
      !> Point preparation: screening plus `pd*` / `qn*` accumulation
      procedure, public :: prepare => lsf_prepare
      !> Point preparation with a caller-provided candidate atom list
      procedure, public :: prepare_subset => lsf_prepare_subset
      !> Configure the highest total derivative order `prepare` provisions
      procedure, public :: set_max_deriv => lsf_set_max_deriv
      !> Number of atoms active after the latest prepare/prepare_subset
      procedure, public :: active_count => lsf_active_count
      !> User-space atom index of the i-th active atom
      procedure, public :: active_atom => lsf_active_atom
      !> Value only (lowest-cost path)
      procedure, public :: f0 => lsf_f0
      !> Combined value/gradient/Hessian
      procedure, public :: f012_r => lsf_f012_r
      !> Third spatial derivative
      procedure, public :: f3_rrr => lsf_f3_rrr
      !> Mixed third derivative: spatial Hessian w.r.t. nuclear positions
      procedure, public :: f3_rr_rA => lsf_f3_rr_rA
      !> Pure spatial fourth derivative
      procedure, public :: f4_rrrr => lsf_f4_rrrr
      !> Mixed fourth derivative d^4S / dr^3 dR_A
      procedure, public :: f4_rrr_rA => lsf_f4_rrr_rA
      !> Radius derivative d^3S / (dr^2 dR_a) and its lower orders
      procedure, public :: f3_rr_rad => lsf_f3_rr_rad
      !> Pure radius Hessian d^2S / dR_a dR_b
      procedure, public :: f2_radrad => lsf_f2_radrad
      !> Mixed third derivative d^3S / dr dR_a dR_b
      procedure, public :: f3_r_radrad => lsf_f3_r_radrad
      !> Mixed fourth derivative d^4S / dr^2 dR_a dR_b
      procedure, public :: f4_rr_radrad => lsf_f4_rr_radrad
      !> Mixed nuclear-radius Hessian d^2S / dR_A dR_b
      procedure, public :: f2_rA_rad => lsf_f2_rA_rad
      !> Mixed third derivative d^3S / dr dR_A dR_b
      procedure, public :: f3_r_rA_rad => lsf_f3_r_rA_rad
      !> Mixed fourth derivative d^4S / dr^2 dR_A dR_b
      procedure, public :: f4_rr_rA_rad => lsf_f4_rr_rA_rad
      !> Pure nuclear Hessian d^2S / dR_A dR_B
      procedure, public :: f2_rArB => lsf_f2_rArB
      !> Mixed third derivative d^3S / dr dR_A dR_B
      procedure, public :: f3_r_rArB => lsf_f3_r_rArB
      !> Mixed fourth derivative d^4S / dr^2 dR_A dR_B
      procedure, public :: f4_rr_rArB => lsf_f4_rr_rArB
      !> Normalized level set f0 / ||f1_r|| and its nuclear derivative
      procedure, public :: normalized_f01_rA => lsf_normalized_f01_rA
      !> Directional nuclear derivative of the value
      procedure, public :: tangent_f0 => lsf_tangent_f0
      !> Directional nuclear derivative of the spatial gradient
      procedure, public :: tangent_f1_r => lsf_tangent_f1_r
      !> Directional nuclear derivative of the spatial Hessian
      procedure, public :: tangent_f2_rr => lsf_tangent_f2_rr
      !> Directional nuclear derivative of the third spatial derivative
      procedure, public :: tangent_f3_rrr => lsf_tangent_f3_rrr
      !> Nuclear Hessian-vector product
      procedure, public :: hvp_f1_rA => lsf_hvp_f1_rA
      !> Directional nuclear derivative of `f2_r_rA`
      procedure, public :: hvp_f2_r_rA => lsf_hvp_f2_r_rA
      !> Directional nuclear derivative of `f3_rr_rA`
      procedure, public :: hvp_f3_rr_rA => lsf_hvp_f3_rr_rA
      !> Radius row of the joint Hessian-vector product
      procedure, public :: hvp_f1_rad => lsf_hvp_f1_rad
      !> Joint directional derivative of `f2_r_rad`
      procedure, public :: hvp_f2_r_rad => lsf_hvp_f2_r_rad
      !> Joint directional derivative of `f3_rr_rad`
      procedure, public :: hvp_f3_rr_rad => lsf_hvp_f3_rr_rad
      !> Adjoint jet contracted onto the nuclear gradient
      procedure, public :: vjp_f1_rA => lsf_vjp_f1_rA
      !> Adjoint jet contracted onto the radius gradient
      procedure, public :: vjp_f1_rad => lsf_vjp_f1_rad
      !> Exact radial offset where the CFC contribution equals the threshold
      procedure, public :: screening_offset => lsf_screening_offset
      !> Finalizer
      final :: finalize_lsf_cfc
   end type moist_cavity_drop_lsf_cfc_type

   public :: moist_cavity_drop_lsf_cfc_type

contains

   !* ================================================================================= *!
   !*                              LSF lifecycle methods                                *!
   !* ================================================================================= *!

   !> Configure CFC shape parameters and declare the candidate index space
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    a1    Atomic-term exponent (optional, default -15)
   !> @param[in]    a2    Pair-term exponent (optional, default -9)
   !> @param[in]    c     Pair-term coupling (optional, default 5)
   !> @param[in]    m     Pair-term power (optional, default 4; the kernel bakes
   !>                     m = 4 into its derivatives, so a different value is
   !>                     only honoured by the screening offset)
   subroutine lsf_new(self, a1, a2, c, m)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Atomic-term exponent override (optional)
      real(wp), intent(in), optional :: a1
      !> Pair-term exponent override (optional)
      real(wp), intent(in), optional :: a2
      !> Pair-term coupling override (optional)
      real(wp), intent(in), optional :: c
      !> Pair-term power override (optional; screening only)
      integer, intent(in), optional :: m

      ! The atom / pair sweeps index the base's candidate-space geometry mirror
      ! directly, so candidate ids must arrive spatially sorted.
      self%candidate_space = lsf_candidate_space_sorted

      call self%param%new(a1=a1, a2=a2, c=c, m=m)
   end subroutine lsf_new

   !> Bind molecular geometry and resize the per-atom caches
   !>
   !> The base handles the spatial sort, the candidate-space geometry mirror and
   !> the screening bounds; this override only sizes the per-point buffers, once,
   !> to the molecule's atom count.
   !>
   !> The `qn*` buffers are zeroed here in full. Per point only the provisioned
   !> orders are refreshed, so this is what keeps the unprovisioned ones exact
   !> zeros instead of undefined memory.
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    mol    Molecular structure
   !> @param[in]    radii  Per-atom radii (size mol%nat)
   subroutine lsf_update(self, mol, radii)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)

      !> Atom capacity of the per-point buffers
      integer :: n_alloc

      call lsf_base_update(self, mol, radii)

      n_alloc = mol%nat
      if (allocated(self%active_cand)) deallocate (self%active_cand)
      if (allocated(self%act_atom)) deallocate (self%act_atom)
      if (allocated(self%act_d)) deallocate (self%act_d)
      if (allocated(self%act_radius)) deallocate (self%act_radius)
      if (allocated(self%qn0)) deallocate (self%qn0)
      if (allocated(self%qn1_r)) deallocate (self%qn1_r)
      if (allocated(self%qn2_rr)) deallocate (self%qn2_rr)
      if (allocated(self%qn3_rrr)) deallocate (self%qn3_rrr)

      allocate (self%active_cand(n_alloc))
      allocate (self%act_atom(n_alloc))
      allocate (self%act_d(ndim, n_alloc))
      allocate (self%act_radius(n_alloc))
      allocate (self%qn0(ndim, n_alloc), source=0.0_wp)
      allocate (self%qn1_r(ndim, ndim, n_alloc), source=0.0_wp)
      allocate (self%qn2_rr(ndim, ndim, ndim, n_alloc), source=0.0_wp)
      allocate (self%qn3_rrr(ndim, ndim, ndim, ndim, n_alloc), source=0.0_wp)

      self%n_active = 0
      self%prepared_deriv = -1
   end subroutine lsf_update

   !> Configure the highest total derivative order `prepare` provisions
   !>
   !> Nothing is allocated here: the aggregate tensors are fixed-size components
   !> and the per-atom buffers are sized by [[update]] to the full order set.
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    n    Requested max derivative order (0..4)
   subroutine lsf_set_max_deriv(self, n)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Requested max derivative order
      integer, intent(in) :: n

      self%max_deriv = min(4, max(0, n))
   end subroutine lsf_set_max_deriv

   !> Screen every atom at the evaluation point and refresh the caches
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    point Evaluation point (3,)
   !> @param[out]   error Evaluation failure (never set by CFC)
   subroutine lsf_prepare(self, point, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      call lsf_cfc_screen(self, point, self%full_scan_cand)
   end subroutine lsf_prepare

   !> Screen a caller-provided candidate list and refresh the caches
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point (3,)
   !> @param[in]    candidate_indices Atom ids in the base's sorted candidate space
   !> @param[out]   error             Evaluation failure (never set by CFC)
   subroutine lsf_prepare_subset(self, point, candidate_indices, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      call lsf_cfc_screen(self, point, candidate_indices)
   end subroutine lsf_prepare_subset

   !> Run the base screening gate and rebuild the per-point accumulators
   !>
   !> Pass one is the base's allocation-free reject test; pass two touches only
   !> the survivors and fills the per-atom geometry; pass three is the atom sweep
   !> and pass four the pair sweep, each accumulating the spatial family and --
   !> when a nuclear order was provisioned -- the one-nuclear-index family in the
   !> same iteration.
   !>
   !> The pair sweep is the O(n_active**2) part and is the single most expensive
   !> thing the CFC level set does. It cannot be made linear (see the module
   !> header); what it can be, and is, is a single kernel call per pair that
   !> yields both families already assembled, with no split-to-spatial
   !> recombination left to this module.
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point [ndim]
   !> @param[in]    candidate_indices Candidate ids in the sorted candidate space
   subroutine lsf_cfc_screen(self, point, candidate_indices)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)

      !> Loop counters, candidate id, survivor count
      integer :: i, j, c, n
      !> Provisioned spatial order and the nuclear spatial order derived from it
      integer :: md, nd

      self%n_active = 0
      self%prepared_deriv = self%max_deriv
      if (.not. allocated(self%cand_screen)) return

      call self%screen_candidates(point, candidate_indices, self%active_cand, n)
      self%n_active = n

      md = self%max_deriv
      nd = md - 1

      self%pd0 = 0.0_wp
      if (md >= 1) self%pd1_r = 0.0_wp
      if (md >= 2) self%pd2_rr = 0.0_wp
      if (md >= 3) self%pd3_rrr = 0.0_wp
      if (md >= 4) self%pd4_rrrr = 0.0_wp
      if (n == 0) return

      if (nd >= 0) self%qn0(:, 1:n) = 0.0_wp
      if (nd >= 1) self%qn1_r(:, :, 1:n) = 0.0_wp
      if (nd >= 2) self%qn2_rr(:, :, :, 1:n) = 0.0_wp
      if (nd >= 3) self%qn3_rrr(:, :, :, :, 1:n) = 0.0_wp

      do i = 1, n
         c = self%active_cand(i)
         self%act_atom(i) = self%cand_to_user(c)
         self%act_radius(i) = self%cand_radii(c)
         self%act_d(1, i) = point(1) - self%cand_screen(1, c)
         self%act_d(2, i) = point(2) - self%cand_screen(2, c)
         self%act_d(3, i) = point(3) - self%cand_screen(3, c)
      end do

      do i = 1, n
         call cfc_atomic_term_eval(self%act_d(:, i), self%act_radius(i), &
                                   self%param%a1, md, &
                                   self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                                   self%pd4_rrrr)
         if (nd < 0) cycle
         call cfc_atomic_nuclear_eval(self%act_d(:, i), self%act_radius(i), &
                                      self%param%a1, nd, &
                                      self%qn0(:, i), self%qn1_r(:, :, i), &
                                      self%qn2_rr(:, :, :, i), self%qn3_rrr(:, :, :, :, i))
      end do

      do i = 1, n
         do j = i + 1, n
            call cfc_pair_spatial_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%param%a2, self%param%c, md, &
                                       self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                                       self%pd4_rrrr)
            if (nd < 0) cycle
            call cfc_pair_nuclear_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%param%a2, self%param%c, nd, &
                                       self%qn0(:, i), self%qn0(:, j), &
                                       self%qn1_r(:, :, i), self%qn1_r(:, :, j), &
                                       self%qn2_rr(:, :, :, i), self%qn2_rr(:, :, :, j), &
                                       self%qn3_rrr(:, :, :, :, i), self%qn3_rrr(:, :, :, :, j))
         end do
      end do
   end subroutine lsf_cfc_screen

   !> Number of atoms currently active after the latest prepare/prepare_subset
   !>
   !> @param[in] self LSF instance
   !> @returns   n    Active-atom count
   pure integer function lsf_active_count(self) result(n)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      n = self%n_active
   end function lsf_active_count

   !> User-space atom index of the i-th currently active atom
   !>
   !> @param[in] self LSF instance
   !> @param[in] i    Active-list index (1 <= i <= active_count())
   !> @returns   idx  User-space atom id
   pure integer function lsf_active_atom(self, i) result(idx)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: i
      idx = self%act_atom(i)
   end function lsf_active_atom

   !* ================================================================================= *!
   !*                             Pure spatial derivatives                              *!
   !* ================================================================================= *!

   !> Level-set value only
   !>
   !> @param[in]  self LSF instance
   !> @param[out] val  LSF value at the prepared point
   subroutine lsf_f0(self, val)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out) :: val

      !> Unused higher-order kernel outputs
      real(wp) :: d1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)

      val = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(0, "f0")

      call cfc_spatial_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            self%pd4_rrrr, 0, val, d1, d2, d3, d4)
   end subroutine lsf_f0

   !> Level-set value, spatial gradient and spatial Hessian
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    LSF value (optional)
   !> @param[out] lsf1_r  Spatial gradient (optional)
   !> @param[out] lsf2_rr Spatial Hessian (optional)
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:, :)

      !> Kernel outputs
      real(wp) :: f0, f1_r(ndim), f2_rr(ndim, ndim), d3(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Highest spatial-derivative order requested
      integer :: md

      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%n_active == 0) return

      md = 0
      if (present(lsf1_r)) md = 1
      if (present(lsf2_rr)) md = 2
      call self%require_deriv(md, "f012_r")

      call cfc_spatial_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            self%pd4_rrrr, md, f0, f1_r, f2_rr, d3, d4)

      if (present(lsf0)) lsf0 = f0
      if (present(lsf1_r)) lsf1_r = f1_r
      if (present(lsf2_rr)) lsf2_rr = f2_rr
   end subroutine lsf_f012_r

   !> Third spatial derivative, plus optionally the lower orders
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     LSF value (optional)
   !> @param[out] lsf1_r   Spatial gradient (optional)
   !> @param[out] lsf2_rr  Spatial Hessian (optional)
   !> @param[out] lsf3_rrr Third spatial derivative [3, 3, 3]
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:, :)
      !> Third spatial derivative
      real(wp), intent(out) :: lsf3_rrr(:, :, :)

      !> Kernel outputs
      real(wp) :: f0, f1_r(ndim), f2_rr(ndim, ndim), f3_rrr(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)

      lsf3_rrr = 0.0_wp
      if (present(lsf0)) lsf0 = 0.0_wp
      if (present(lsf1_r)) lsf1_r = 0.0_wp
      if (present(lsf2_rr)) lsf2_rr = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(3, "f3_rrr")

      call cfc_spatial_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            self%pd4_rrrr, 3, f0, f1_r, f2_rr, f3_rrr, d4)

      if (present(lsf0)) lsf0 = f0
      if (present(lsf1_r)) lsf1_r = f1_r
      if (present(lsf2_rr)) lsf2_rr = f2_rr
      lsf3_rrr = f3_rrr
   end subroutine lsf_f3_rrr

   !> Fourth spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Fourth spatial derivative [3, 3, 3, 3]
   subroutine lsf_f4_rrrr(self, lsf4_rrrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Fourth spatial derivative
      real(wp), intent(out) :: lsf4_rrrr(:, :, :, :)

      !> Kernel outputs
      real(wp) :: f0, f1_r(ndim), f2_rr(ndim, ndim), f3_rrr(ndim, ndim, ndim)
      real(wp) :: f4_rrrr(ndim, ndim, ndim, ndim)

      lsf4_rrrr = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(4, "f4_rrrr")

      call cfc_spatial_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            self%pd4_rrrr, 4, f0, f1_r, f2_rr, f3_rrr, f4_rrrr)
      lsf4_rrrr = f4_rrrr
   end subroutine lsf_f4_rrrr

   !* ================================================================================= *!
   !*                          One-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Mixed third derivative d^3S / (dr^2 dR_A) and its lower orders
   !>
   !> All three outputs are active-indexed: slot `i` belongs to `active_atom(i)`.
   !> The order requirement is 3, not 2: the kernel reads `qn2_rr`, whose total
   !> derivative order is three.
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    dS/dR_A [3, >= active_count()] (optional)
   !> @param[out] lsf2_r_rA  d^2S/(dr dR_A) [3, 3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rA d^3S/(dr^2 dR_A) [3, 3, 3, >= active_count()]
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear gradient
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :, :)

      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(3, "f3_rr_rA")

      do ia = 1, self%n_active
         call cfc_nuclear_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                               self%qn0(:, ia), self%qn1_r(:, :, ia), &
                               self%qn2_rr(:, :, :, ia), self%qn3_rrr(:, :, :, :, ia), &
                               2, f1_rA, f2_r_rA, f3_rr_rA, d4)
         if (present(lsf1_rA)) lsf1_rA(:, ia) = f1_rA
         if (present(lsf2_r_rA)) lsf2_r_rA(:, :, ia) = f2_r_rA
         lsf3_rr_rA(:, :, :, ia) = f3_rr_rA
      end do
   end subroutine lsf_f3_rr_rA

   !> Mixed fourth derivative d^4S / (dr^3 dR_A), active-indexed
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf4_rrr_rA d^4S/(dr^3 dR_A) [3, 3, 3, 3, >= active_count()]
   subroutine lsf_f4_rrr_rA(self, lsf4_rrr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rrr_rA(:, :, :, :, :)

      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: f4_rrr_rA(ndim, ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(4, "f4_rrr_rA")

      do ia = 1, self%n_active
         call cfc_nuclear_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                               self%qn0(:, ia), self%qn1_r(:, :, ia), &
                               self%qn2_rr(:, :, :, ia), self%qn3_rrr(:, :, :, :, ia), &
                               3, f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
         lsf4_rrr_rA(:, :, :, :, ia) = f4_rrr_rA
      end do
   end subroutine lsf_f4_rrr_rA

   !* ================================================================================= *!
   !*                               Radius derivatives                                  *!
   !* ================================================================================= *!
   !
   ! The radii are parameters of the level set, so a model that ties them to the
   ! geometry differentiates through them as well. Three things separate this
   ! channel from the `_rA` one, all of them because `R_a` is a scalar parameter
   ! and not a coordinate: no sign flip, no trailing derivative index (the rank
   ! equals the spatial order), and no derivative order consumed.
   !
   ! Unlike `qn*`, the `rd*` family is *not* cached by `prepare`. It could be --
   ! it depends on nothing but the geometry -- but that would charge every
   ! fixed-radius caller an extra O(n_active**2) pair sweep at every grid point
   ! for a family it never reads. It is accumulated on demand instead, the same
   ! way `tg*` is.

   !> Accumulate the per-atom radius-derivative pseudo-density tensors
   !>
   !> The radius counterpart of the `qn*` sweep in [[lsf_cfc_screen]], run on
   !> demand. Slot `i` of each output belongs to `active_atom(i)`.
   !>
   !> `rd3_rrr` is optional because [[cfc_radius_hvp]] needs the family only up
   !> to order two and would otherwise have to materialize an unread
   !> `27 * n_active` buffer to satisfy the kernel signature. Omitting it sinks
   !> the order-3 accumulation into two stack rows instead. The loops are spelled
   !> out twice rather than once behind a pointer: the kernel arguments are
   !> `intent(inout)` accumulators taken one atom at a time, so the only
   !> single-body alternative is an aliasing-opaque pointer with `merge`-computed
   !> slot indices in the pair loop.
   !>
   !> @param[in]  self    LSF instance
   !> @param[in]  level   Highest spatial-derivative order to accumulate (0..3)
   !> @param[out] rd0     Order-0 radius tensor [n_active]
   !> @param[out] rd1_r   Order-1 radius tensor [ndim, n_active]
   !> @param[out] rd2_rr  Order-2 radius tensor [ndim, ndim, n_active]
   !> @param[out] rd3_rrr Order-3 radius tensor [ndim, ndim, ndim, n_active]
   !>                     (optional; discarded into stack scratch when absent)
   subroutine radius_tensors(self, level, rd0, rd1_r, rd2_rr, rd3_rrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: level
      !> Order-0 radius tensor
      real(wp), intent(out) :: rd0(:)
      !> Order-1 radius tensor
      real(wp), intent(out) :: rd1_r(:, :)
      !> Order-2 radius tensor
      real(wp), intent(out) :: rd2_rr(:, :, :)
      !> Order-3 radius tensor
      real(wp), intent(out), optional :: rd3_rrr(:, :, :, :)

      !> Order-3 sink, one row per kernel slot so the pair call has no aliasing
      real(wp) :: sk3_a(ndim, ndim, ndim), sk3_b(ndim, ndim, ndim)
      !> Active-list indices
      integer :: ia, ib

      rd0 = 0.0_wp
      rd1_r = 0.0_wp
      rd2_rr = 0.0_wp

      if (present(rd3_rrr)) then
         rd3_rrr = 0.0_wp

         do ia = 1, self%n_active
            call cfc_atomic_radius_eval(self%act_d(:, ia), self%act_radius(ia), &
                                        self%param%a1, &
                                        level, rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                        rd3_rrr(:, :, :, ia))
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_radius_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                         self%act_radius(ia), self%act_radius(ib), &
                                         self%param%a2, self%param%c, level, &
                                         rd0(ia), rd0(ib), &
                                         rd1_r(:, ia), rd1_r(:, ib), &
                                         rd2_rr(:, :, ia), rd2_rr(:, :, ib), &
                                         rd3_rrr(:, :, :, ia), rd3_rrr(:, :, :, ib))
            end do
         end do
      else
         ! Zeroed once, not per atom: the sinks are `intent(inout)` accumulators
         ! and reading them undefined would trap under the debug FP settings.
         sk3_a = 0.0_wp
         sk3_b = 0.0_wp

         do ia = 1, self%n_active
            call cfc_atomic_radius_eval(self%act_d(:, ia), self%act_radius(ia), &
                                        self%param%a1, &
                                        level, rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                        sk3_a)
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_radius_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                         self%act_radius(ia), self%act_radius(ib), &
                                         self%param%a2, self%param%c, level, &
                                         rd0(ia), rd0(ib), &
                                         rd1_r(:, ia), rd1_r(:, ib), &
                                         rd2_rr(:, :, ia), rd2_rr(:, :, ib), &
                                         sk3_a, sk3_b)
            end do
         end do
      end if
   end subroutine radius_tensors

   !> Radius derivative d^3S / (dr^2 dR_a) and its lower orders
   !>
   !> All three outputs are active-indexed: slot `i` belongs to `active_atom(i)`.
   !> `R_a` is the *radius* of that atom, so the results carry no trailing
   !> derivative index and the rank equals the spatial order.
   !>
   !> The order requirement is 2, not 3: `d/dR_a` consumes no derivative order,
   !> so the deepest input is `rd2_rr`, whose total order is two -- one less than
   !> `lsf_f3_rr_rA` needs for the same spatial order.
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf1_rad    dS/dR_a [>= active_count()] (optional)
   !> @param[out] lsf2_r_rad  d^2S/(dr dR_a) [3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rad d^3S/(dr^2 dR_a) [3, 3, >= active_count()]
   subroutine lsf_f3_rr_rad(self, lsf1_rad, lsf2_r_rad, lsf3_rr_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Radius gradient
      real(wp), intent(out), optional :: lsf1_rad(:)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rad(:, :)
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_rr_rad(:, :, :)

      !> Per-atom radius tensors
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rd3_rrr(:, :, :, :)
      !> Kernel outputs of one atom
      real(wp) :: f1_rad, f2_r_rad(ndim), f3_rr_rad(ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         ! Nothing active means no owned slot, so "writes the first
         ! `active_count()` slots" degenerates to writing none. Zero rather than
         ! return bare: the outputs are `intent(out)`, and a bare return would
         ! hand the caller a buffer it is not allowed to read.
         if (present(lsf1_rad)) lsf1_rad = 0.0_wp
         if (present(lsf2_r_rad)) lsf2_r_rad = 0.0_wp
         lsf3_rr_rad = 0.0_wp
         return
      end if
      call self%require_deriv(2, "f3_rr_rad")

      allocate (rd0(self%n_active))
      allocate (rd1_r(ndim, self%n_active))
      allocate (rd2_rr(ndim, ndim, self%n_active))
      allocate (rd3_rrr(ndim, ndim, ndim, self%n_active))
      call radius_tensors(self, 2, rd0, rd1_r, rd2_rr, rd3_rrr)

      do ia = 1, self%n_active
         call cfc_radius_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                              rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                              rd3_rrr(:, :, :, ia), 2, &
                              f1_rad, f2_r_rad, f3_rr_rad, d4)
         if (present(lsf1_rad)) lsf1_rad(ia) = f1_rad
         if (present(lsf2_r_rad)) lsf2_r_rad(:, ia) = f2_r_rad
         lsf3_rr_rad(:, :, ia) = f3_rr_rad
      end do
   end subroutine lsf_f3_rr_rad

   !> Normalized level set S/||grad S|| and its nuclear gradient
   !>
   !> The normalization itself is level-set agnostic -- it only combines outputs
   !> the kernel already returns -- so it is spelled out here rather than
   !> generated, matching `svdw_normalized_eval` term for term including the
   !> degenerate-gradient guard.
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] val      S/||grad S||
   !> @param[out] deriv_rA d/dR_A of S/||grad S|| [3, >= active_count()] (optional)
   subroutine lsf_normalized_f01_rA(self, val, deriv_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Normalized level-set value
      real(wp), intent(out) :: val
      !> Nuclear gradient of the normalized value
      real(wp), intent(out), optional :: deriv_rA(:, :)

      !> Pure spatial outputs
      real(wp) :: f0, f1_r(ndim), f2_rr(ndim, ndim), d3(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> One-nuclear-index outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), d3n(ndim, ndim, ndim)
      real(wp) :: d4n(ndim, ndim, ndim, ndim)
      !> Gradient norm, its reciprocal powers and the gradient/mixed contraction
      real(wp) :: g, inv_g, inv_g3, gd(ndim)
      !> Active-list index and Cartesian index
      integer :: ia, j

      val = 0.0_wp
      if (self%n_active == 0) return
      ! The nuclear half reads qn1_r (total order 2); the value half only needs
      ! pd1_r, but the two share one entry point so the stricter order rules.
      call self%require_deriv(2, "normalized_f01_rA")

      call cfc_spatial_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            self%pd4_rrrr, 1, f0, f1_r, f2_rr, d3, d4)

      g = sqrt(f1_r(1)*f1_r(1) + f1_r(2)*f1_r(2) + f1_r(3)*f1_r(3))
      if (g < normalized_eps) then
         if (present(deriv_rA)) deriv_rA(:, 1:self%n_active) = 0.0_wp
         return
      end if

      inv_g = 1.0_wp/g
      inv_g3 = inv_g*inv_g*inv_g
      val = f0*inv_g
      if (.not. present(deriv_rA)) return

      do ia = 1, self%n_active
         call cfc_nuclear_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                               self%qn0(:, ia), self%qn1_r(:, :, ia), &
                               self%qn2_rr(:, :, :, ia), self%qn3_rrr(:, :, :, :, ia), &
                               1, f1_rA, f2_r_rA, d3n, d4n)
         gd = 0.0_wp
         do j = 1, ndim
            gd = gd + f1_r(j)*f2_r_rA(j, :)
         end do
         deriv_rA(:, ia) = f1_rA*inv_g - f0*gd*inv_g3
      end do
   end subroutine lsf_normalized_f01_rA

   !* ================================================================================= *!
   !*                          Two-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Fill the uncontracted two-nucleus family at one spatial order
   !>
   !> Shared body of the three `f*_rArB` accessors. Two sweeps:
   !>
   !>   1. the *diagonal* `qq` blocks, one per active atom: atom A's own term
   !>      plus the `aa` / `bb` block of every pair containing A. One pass over
   !>      the atoms and one over the unordered pairs.
   !>   2. the ordered-pair loop, which for `A /= B` gets its `qq` block from the
   !>      single pair term containing both centers (no other term survives two
   !>      derivatives with respect to two different nuclei) and for `A == B`
   !>      reuses the diagonal accumulated in pass one, then lifts.
   !>
   !> The `qq` scratch is a local allocation rather than a persistent buffer for
   !> the same reason SvdW's pair cache is: only these three accessors need it,
   !> none of them sits on the projection hot path, and a persistent order-2
   !> two-nucleus cache would cost 936 bytes per atom per thread.
   !>
   !> @param[in]  self         LSF instance
   !> @param[in]  level        Spatial order to lift (0, 1 or 2)
   !> @param[out] lsf2_rArB    d^2S/(dR_A dR_B) (optional)
   !> @param[out] lsf3_r_rArB  d^3S/(dr dR_A dR_B) (optional)
   !> @param[out] lsf4_rr_rArB d^4S/(dr^2 dR_A dR_B) (optional)
   subroutine hessian_family(self, level, lsf2_rArB, lsf3_r_rArB, lsf4_rr_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Spatial order to lift
      integer, intent(in) :: level
      !> Nuclear Hessian
      real(wp), intent(out), optional :: lsf2_rArB(:, :, :, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: lsf3_r_rArB(:, :, :, :, :)
      !> Mixed fourth derivative
      real(wp), intent(out), optional :: lsf4_rr_rArB(:, :, :, :, :, :)

      !> Diagonal two-nucleus accumulators, one block per active atom
      real(wp), allocatable :: qd0(:, :, :), qd1(:, :, :, :), qd2(:, :, :, :, :)
      !> Cross block of one ordered pair, and the two diagonal blocks it also
      !> produces and this loop discards
      real(wp) :: x0(ndim, ndim), x1(ndim, ndim, ndim), x2(ndim, ndim, ndim, ndim)
      real(wp) :: w0a(ndim, ndim), w1a(ndim, ndim, ndim), w2a(ndim, ndim, ndim, ndim)
      real(wp) :: w0b(ndim, ndim), w1b(ndim, ndim, ndim), w2b(ndim, ndim, ndim, ndim)
      !> Lift outputs of one ordered pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices, active-atom count and tensor components
      integer :: ia, ib, n, s, t, j, k

      n = self%n_active

      allocate (qd0(ndim, ndim, n), source=0.0_wp)
      allocate (qd1(ndim, ndim, ndim, n), source=0.0_wp)
      allocate (qd2(ndim, ndim, ndim, ndim, n), source=0.0_wp)

      do ia = 1, n
         call cfc_atomic_hessian_eval(self%act_d(:, ia), self%act_radius(ia), &
                                      self%param%a1, &
                                      level, qd0(:, :, ia), qd1(:, :, :, ia), &
                                      qd2(:, :, :, :, ia))
      end do
      do ia = 1, n
         do ib = ia + 1, n
            x0 = 0.0_wp
            x1 = 0.0_wp
            x2 = 0.0_wp
            call cfc_pair_hessian_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                       self%act_radius(ia), self%act_radius(ib), &
                                       self%param%a2, self%param%c, level, &
                                       qd0(:, :, ia), x0, qd0(:, :, ib), &
                                       qd1(:, :, :, ia), x1, qd1(:, :, :, ib), &
                                       qd2(:, :, :, :, ia), x2, qd2(:, :, :, :, ib))
         end do
      end do

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               x0 = qd0(:, :, ia)
               x1 = qd1(:, :, :, ia)
               x2 = qd2(:, :, :, :, ia)
            else
               x0 = 0.0_wp
               x1 = 0.0_wp
               x2 = 0.0_wp
               w0a = 0.0_wp
               w1a = 0.0_wp
               w2a = 0.0_wp
               w0b = 0.0_wp
               w1b = 0.0_wp
               w2b = 0.0_wp
               call cfc_pair_hessian_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                          self%act_radius(ia), self%act_radius(ib), &
                                          self%param%a2, self%param%c, level, &
                                          w0a, x0, w0b, w1a, x1, w1b, w2a, x2, w2b)
            end if

            call cfc_hessian_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                  self%qn0(:, ia), self%qn1_r(:, :, ia), &
                                  self%qn2_rr(:, :, :, ia), &
                                  self%qn0(:, ib), self%qn1_r(:, :, ib), &
                                  self%qn2_rr(:, :, :, ib), &
                                  x0, x1, x2, level, b2, b3, b4)

            if (present(lsf2_rArB)) lsf2_rArB(:, ia, :, ib) = b2
            if (present(lsf3_r_rArB)) then
               do t = 1, ndim
                  do s = 1, ndim
                     lsf3_r_rArB(:, s, ia, t, ib) = b3(:, s, t)
                  end do
               end do
            end if
            if (present(lsf4_rr_rArB)) then
               do t = 1, ndim
                  do s = 1, ndim
                     do k = 1, ndim
                        do j = 1, ndim
                           lsf4_rr_rArB(j, k, s, ia, t, ib) = b4(j, k, s, t)
                        end do
                     end do
                  end do
               end do
            end if
         end do
      end do
   end subroutine hessian_family

   !> Pure nuclear Hessian d^2S / (dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf2_rArB d^2S/(dR_A dR_B) [3, >= n_act, 3, >= n_act]
   subroutine lsf_f2_rArB(self, lsf2_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear Hessian
      real(wp), intent(out) :: lsf2_rArB(:, :, :, :)

      if (self%n_active == 0) return
      call self%require_deriv(1, "f2_rArB")

      call hessian_family(self, 0, lsf2_rArB=lsf2_rArB)
   end subroutine lsf_f2_rArB

   !> Mixed third derivative d^3S / (dr dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf3_r_rArB d^3S/(dr dR_A dR_B) [3, 3, >= n_act, 3, >= n_act]
   subroutine lsf_f3_r_rArB(self, lsf3_r_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_rArB(:, :, :, :, :)

      if (self%n_active == 0) return
      call self%require_deriv(2, "f3_r_rArB")

      call hessian_family(self, 1, lsf3_r_rArB=lsf3_r_rArB)
   end subroutine lsf_f3_r_rArB

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf4_rr_rArB d^4S/(dr^2 dR_A dR_B) [3, 3, 3, >= n_act, 3, >= n_act]
   subroutine lsf_f4_rr_rArB(self, lsf4_rr_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rr_rArB(:, :, :, :, :, :)

      if (self%n_active == 0) return
      call self%require_deriv(3, "f4_rr_rArB")

      call hessian_family(self, 2, lsf4_rr_rArB=lsf4_rr_rArB)
   end subroutine lsf_f4_rr_rArB

   !* ================================================================================= *!
   !*                     Two-radius and nuclear-radius derivatives                     *!
   !* ================================================================================= *!
   !
   ! The uncontracted second order of the joint position/radius derivative,
   ! completing the four blocks whose contractions `hvp_*_rA` and `hvp_*_rad`
   ! already supply. Structurally these are [[hessian_family]] with one or both
   ! retained nuclear slots replaced by a retained radius, and cheaper for it: a
   ! radius adds no Cartesian index and consumes no derivative order.
   !
   ! One difference in the loop shape. `hessian_family` visits every *ordered*
   ! pair and re-evaluates the cross block each time, because the `ab` block
   ! transposed is the `ba` one. That holds for the radius-radius family too, but
   ! not for the nuclear-radius one -- its two slots carry different kinds of derivative --
   ! so the kernel emits `ba` alongside `ab` and both families visit each
   ! *unordered* pair once, filling the two ordered slots from a single call.
   !
   ! The six accessors below zero their result before the empty-active return,
   ! matching the `tangent_*` family. Elsewhere they write only the first
   ! `active_count()` slots and leave the rest alone, but with nothing active
   ! there is no owned slot at all, and the results are `intent(out)`: a bare
   ! return would hand the caller a buffer it is not allowed to read.

   !> Fill the uncontracted two-radius family at one spatial order
   !>
   !> Shared body of the three `f*_radrad` accessors; see the section note above.
   !>
   !> @param[in]  self           LSF instance
   !> @param[in]  level          Spatial order to lift (0, 1 or 2)
   !> @param[out] lsf2_radrad    d^2S/(dR_a dR_b) (optional)
   !> @param[out] lsf3_r_radrad  d^3S/(dr dR_a dR_b) (optional)
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b) (optional)
   subroutine radius2_family(self, level, lsf2_radrad, lsf3_r_radrad, lsf4_rr_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Spatial order to lift
      integer, intent(in) :: level
      !> Radius Hessian
      real(wp), intent(out), optional :: lsf2_radrad(:, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: lsf3_r_radrad(:, :, :)
      !> Mixed fourth derivative
      real(wp), intent(out), optional :: lsf4_rr_radrad(:, :, :, :)

      !> Per-atom radius tensors, in both the A and the B slot of the lift
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rd3_rrr(:, :, :, :)
      !> Diagonal two-radius accumulators, one block per active atom
      real(wp), allocatable :: dd0(:), dd1(:, :), dd2(:, :, :)
      !> Cross block of every unordered pair, cached by the single sweep below.
      !> Only the orders `level` will actually lift are allocated.
      real(wp), allocatable :: cx0(:), cx1(:, :), cx2(:, :, :)
      !> Cross block of one pair
      real(wp) :: x0, x1(ndim), x2(ndim, ndim)
      !> Lift outputs of one pair
      real(wp) :: b2, b3(ndim), b4(ndim, ndim)
      !> Active-list indices, active-atom count, pair counter and pair count
      integer :: ia, ib, n, ip, np

      n = self%n_active
      np = n*(n - 1)/2

      allocate (rd0(n), rd1_r(ndim, n), rd2_rr(ndim, ndim, n))
      allocate (rd3_rrr(ndim, ndim, ndim, n))
      call radius_tensors(self, level, rd0, rd1_r, rd2_rr, rd3_rrr)

      allocate (dd0(n), source=0.0_wp)
      allocate (dd1(ndim, n), source=0.0_wp)
      allocate (dd2(ndim, ndim, n), source=0.0_wp)

      allocate (cx0(np))
      allocate (cx1(ndim, merge(np, 0, level >= 1)))
      allocate (cx2(ndim, ndim, merge(np, 0, level >= 2)))

      do ia = 1, n
         call cfc_atomic_radius2_eval(self%act_d(:, ia), self%act_radius(ia), &
                                      self%param%a1, &
                                      level, dd0(ia), dd1(:, ia), dd2(:, :, ia))
      end do
      ! One `cfc_pair_radius2_eval` call yields both the diagonal contributions
      ! and the cross block, so cache the cross here rather than repeating the
      ! whole O(n^2) sweep in the output loop below.
      ip = 0
      do ia = 1, n
         do ib = ia + 1, n
            ip = ip + 1
            x0 = 0.0_wp
            x1 = 0.0_wp
            x2 = 0.0_wp
            call cfc_pair_radius2_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                       self%act_radius(ia), self%act_radius(ib), &
                                       self%param%a2, self%param%c, level, &
                                       dd0(ia), x0, dd0(ib), &
                                       dd1(:, ia), x1, dd1(:, ib), &
                                       dd2(:, :, ia), x2, dd2(:, :, ib))
            cx0(ip) = x0
            if (level >= 1) cx1(:, ip) = x1
            if (level >= 2) cx2(:, :, ip) = x2
         end do
      end do

      ip = 0
      do ia = 1, n
         call cfc_radius2_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                               rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                               rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                               dd0(ia), dd1(:, ia), dd2(:, :, ia), level, b2, b3, b4)
         call radius2_store(ia, ia, b2, b3, b4, &
                            lsf2_radrad, lsf3_r_radrad, lsf4_rr_radrad)

         do ib = ia + 1, n
            ip = ip + 1
            x0 = cx0(ip)
            x1 = 0.0_wp
            x2 = 0.0_wp
            if (level >= 1) x1 = cx1(:, ip)
            if (level >= 2) x2 = cx2(:, :, ip)
            call cfc_radius2_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                  rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                  rd0(ib), rd1_r(:, ib), rd2_rr(:, :, ib), &
                                  x0, x1, x2, level, b2, b3, b4)
            ! Symmetric in its two radii, so one lift fills both ordered slots.
            call radius2_store(ia, ib, b2, b3, b4, &
                               lsf2_radrad, lsf3_r_radrad, lsf4_rr_radrad)
            call radius2_store(ib, ia, b2, b3, b4, &
                               lsf2_radrad, lsf3_r_radrad, lsf4_rr_radrad)
         end do
      end do
   end subroutine radius2_family

   !> Store one lifted two-radius block into whichever outputs are present
   !>
   !> @param[in]    ia             Active-list index of the first radius
   !> @param[in]    ib             Active-list index of the second radius
   !> @param[in]    b2             Order-0 block
   !> @param[in]    b3             Order-1 block
   !> @param[in]    b4             Order-2 block
   !> @param[inout] lsf2_radrad    Radius Hessian (optional)
   !> @param[inout] lsf3_r_radrad  Mixed third derivative (optional)
   !> @param[inout] lsf4_rr_radrad Mixed fourth derivative (optional)
   pure subroutine radius2_store(ia, ib, b2, b3, b4, &
                                 lsf2_radrad, lsf3_r_radrad, lsf4_rr_radrad)
      !> Active-list index of the first radius
      integer, intent(in) :: ia
      !> Active-list index of the second radius
      integer, intent(in) :: ib
      !> Order-0 block
      real(wp), intent(in) :: b2
      !> Order-1 block
      real(wp), intent(in) :: b3(ndim)
      !> Order-2 block
      real(wp), intent(in) :: b4(ndim, ndim)
      !> Radius Hessian
      real(wp), intent(inout), optional :: lsf2_radrad(:, :)
      !> Mixed third derivative
      real(wp), intent(inout), optional :: lsf3_r_radrad(:, :, :)
      !> Mixed fourth derivative
      real(wp), intent(inout), optional :: lsf4_rr_radrad(:, :, :, :)

      if (present(lsf2_radrad)) lsf2_radrad(ia, ib) = b2
      if (present(lsf3_r_radrad)) lsf3_r_radrad(:, ia, ib) = b3
      if (present(lsf4_rr_radrad)) lsf4_rr_radrad(:, :, ia, ib) = b4
   end subroutine radius2_store

   !> Fill the uncontracted nuclear-radius family at one spatial order
   !>
   !> Shared body of the three `f*_rA_rad` accessors. The first atom index
   !> carries the position derivative and the second the radius one, so the
   !> result is *not* symmetric under exchanging them; each unordered pair is
   !> visited once and its `ab` and `ba` kernel blocks fill the two ordered slots.
   !>
   !> @param[in]  self           LSF instance
   !> @param[in]  level          Spatial order to lift (0, 1 or 2)
   !> @param[out] lsf2_rA_rad    d^2S/(dR_A dR_b) (optional)
   !> @param[out] lsf3_r_rA_rad  d^3S/(dr dR_A dR_b) (optional)
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b) (optional)
   subroutine nucrad_family(self, level, lsf2_rA_rad, lsf3_r_rA_rad, lsf4_rr_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Spatial order to lift
      integer, intent(in) :: level
      !> Nuclear-radius Hessian
      real(wp), intent(out), optional :: lsf2_rA_rad(:, :, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: lsf3_r_rA_rad(:, :, :, :)
      !> Mixed fourth derivative
      real(wp), intent(out), optional :: lsf4_rr_rA_rad(:, :, :, :, :)

      !> Per-atom radius tensors, filling the lift's B slot
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rd3_rrr(:, :, :, :)
      !> Diagonal nuclear-radius accumulators, one block per active atom
      real(wp), allocatable :: dd0(:, :), dd1(:, :, :), dd2(:, :, :, :)
      !> Both cross blocks of every unordered pair, cached by the single sweep
      !> below. Only the orders `level` will actually lift are allocated.
      real(wp), allocatable :: cx0(:, :), cx1(:, :, :), cx2(:, :, :, :)
      real(wp), allocatable :: cy0(:, :), cy1(:, :, :), cy2(:, :, :, :)
      !> The two cross blocks of one pair
      real(wp) :: x0(ndim), x1(ndim, ndim), x2(ndim, ndim, ndim)
      real(wp) :: y0(ndim), y1(ndim, ndim), y2(ndim, ndim, ndim)
      !> Lift outputs of one ordered pair
      real(wp) :: b2(ndim), b3(ndim, ndim), b4(ndim, ndim, ndim)
      !> Active-list indices, active-atom count, pair counter and pair count
      integer :: ia, ib, n, ip, np

      n = self%n_active
      np = n*(n - 1)/2

      allocate (rd0(n), rd1_r(ndim, n), rd2_rr(ndim, ndim, n))
      allocate (rd3_rrr(ndim, ndim, ndim, n))
      call radius_tensors(self, level, rd0, rd1_r, rd2_rr, rd3_rrr)

      allocate (dd0(ndim, n), source=0.0_wp)
      allocate (dd1(ndim, ndim, n), source=0.0_wp)
      allocate (dd2(ndim, ndim, ndim, n), source=0.0_wp)

      allocate (cx0(ndim, np), cy0(ndim, np))
      allocate (cx1(ndim, ndim, merge(np, 0, level >= 1)))
      allocate (cy1(ndim, ndim, merge(np, 0, level >= 1)))
      allocate (cx2(ndim, ndim, ndim, merge(np, 0, level >= 2)))
      allocate (cy2(ndim, ndim, ndim, merge(np, 0, level >= 2)))

      do ia = 1, n
         call cfc_atomic_nucrad_eval(self%act_d(:, ia), self%act_radius(ia), &
                                     self%param%a1, &
                                     level, dd0(:, ia), dd1(:, :, ia), dd2(:, :, :, ia))
      end do
      ! `cfc_pair_nucrad_eval` is the most expensive routine in this family, and
      ! one call yields *both* the two diagonal contributions and the two cross
      ! blocks. Cache the cross blocks here so the output loop below can lift
      ! them without a second, identical O(n^2) sweep.
      ip = 0
      do ia = 1, n
         do ib = ia + 1, n
            ip = ip + 1
            x0 = 0.0_wp
            x1 = 0.0_wp
            x2 = 0.0_wp
            y0 = 0.0_wp
            y1 = 0.0_wp
            y2 = 0.0_wp
            call cfc_pair_nucrad_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                      self%act_radius(ia), self%act_radius(ib), &
                                      self%param%a2, self%param%c, level, &
                                      dd0(:, ia), x0, y0, dd0(:, ib), &
                                      dd1(:, :, ia), x1, y1, dd1(:, :, ib), &
                                      dd2(:, :, :, ia), x2, y2, dd2(:, :, :, ib))
            cx0(:, ip) = x0
            cy0(:, ip) = y0
            if (level >= 1) then
               cx1(:, :, ip) = x1
               cy1(:, :, ip) = y1
            end if
            if (level >= 2) then
               cx2(:, :, :, ip) = x2
               cy2(:, :, :, ip) = y2
            end if
         end do
      end do

      ip = 0
      do ia = 1, n
         call cfc_nucrad_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                              self%qn0(:, ia), self%qn1_r(:, :, ia), &
                              self%qn2_rr(:, :, :, ia), &
                              rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                              dd0(:, ia), dd1(:, :, ia), dd2(:, :, :, ia), &
                              level, b2, b3, b4)
         call nucrad_store(ia, ia, b2, b3, b4, &
                           lsf2_rA_rad, lsf3_r_rA_rad, lsf4_rr_rA_rad)

         do ib = ia + 1, n
            ip = ip + 1
            x0 = cx0(:, ip)
            y0 = cy0(:, ip)
            x1 = 0.0_wp
            x2 = 0.0_wp
            y1 = 0.0_wp
            y2 = 0.0_wp
            if (level >= 1) then
               x1 = cx1(:, :, ip)
               y1 = cy1(:, :, ip)
            end if
            if (level >= 2) then
               x2 = cx2(:, :, :, ip)
               y2 = cy2(:, :, :, ip)
            end if
            ! `ab` is position on ia and radius on ib; `ba` is the other way.
            call cfc_nucrad_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                 self%qn0(:, ia), self%qn1_r(:, :, ia), &
                                 self%qn2_rr(:, :, :, ia), &
                                 rd0(ib), rd1_r(:, ib), rd2_rr(:, :, ib), &
                                 x0, x1, x2, level, b2, b3, b4)
            call nucrad_store(ia, ib, b2, b3, b4, &
                              lsf2_rA_rad, lsf3_r_rA_rad, lsf4_rr_rA_rad)
            call cfc_nucrad_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                 self%qn0(:, ib), self%qn1_r(:, :, ib), &
                                 self%qn2_rr(:, :, :, ib), &
                                 rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                 y0, y1, y2, level, b2, b3, b4)
            call nucrad_store(ib, ia, b2, b3, b4, &
                              lsf2_rA_rad, lsf3_r_rA_rad, lsf4_rr_rA_rad)
         end do
      end do
   end subroutine nucrad_family

   !> Store one lifted nuclear-radius block into whichever outputs are present
   !>
   !> @param[in]    ia             Active-list index of the position derivative
   !> @param[in]    ib             Active-list index of the radius derivative
   !> @param[in]    b2             Order-0 block
   !> @param[in]    b3             Order-1 block
   !> @param[in]    b4             Order-2 block
   !> @param[inout] lsf2_rA_rad    Nuclear-radius Hessian (optional)
   !> @param[inout] lsf3_r_rA_rad  Mixed third derivative (optional)
   !> @param[inout] lsf4_rr_rA_rad Mixed fourth derivative (optional)
   pure subroutine nucrad_store(ia, ib, b2, b3, b4, &
                                lsf2_rA_rad, lsf3_r_rA_rad, lsf4_rr_rA_rad)
      !> Active-list index of the position derivative
      integer, intent(in) :: ia
      !> Active-list index of the radius derivative
      integer, intent(in) :: ib
      !> Order-0 block
      real(wp), intent(in) :: b2(ndim)
      !> Order-1 block
      real(wp), intent(in) :: b3(ndim, ndim)
      !> Order-2 block
      real(wp), intent(in) :: b4(ndim, ndim, ndim)
      !> Nuclear-radius Hessian
      real(wp), intent(inout), optional :: lsf2_rA_rad(:, :, :)
      !> Mixed third derivative
      real(wp), intent(inout), optional :: lsf3_r_rA_rad(:, :, :, :)
      !> Mixed fourth derivative
      real(wp), intent(inout), optional :: lsf4_rr_rA_rad(:, :, :, :, :)

      if (present(lsf2_rA_rad)) lsf2_rA_rad(:, ia, ib) = b2
      if (present(lsf3_r_rA_rad)) lsf3_r_rA_rad(:, :, ia, ib) = b3
      if (present(lsf4_rr_rA_rad)) lsf4_rr_rA_rad(:, :, :, ia, ib) = b4
   end subroutine nucrad_store

   !> Pure radius Hessian d^2S / (dR_a dR_b), active-indexed in both atoms
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_radrad d^2S/(dR_a dR_b) [>= n_act, >= n_act]
   subroutine lsf_f2_radrad(self, lsf2_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Radius Hessian
      real(wp), intent(out) :: lsf2_radrad(:, :)

      if (self%n_active == 0) then
         lsf2_radrad = 0.0_wp
         return
      end if
      ! One order below `f2_rArB`: neither retained radius consumes one.
      call self%require_deriv(0, "f2_radrad")

      call radius2_family(self, 0, lsf2_radrad=lsf2_radrad)
   end subroutine lsf_f2_radrad

   !> Mixed third derivative d^3S / (dr dR_a dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_radrad d^3S/(dr dR_a dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_f3_r_radrad(self, lsf3_r_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_radrad(:, :, :)

      if (self%n_active == 0) then
         lsf3_r_radrad = 0.0_wp
         return
      end if
      call self%require_deriv(1, "f3_r_radrad")

      call radius2_family(self, 1, lsf3_r_radrad=lsf3_r_radrad)
   end subroutine lsf_f3_r_radrad

   !> Mixed fourth derivative d^4S / (dr^2 dR_a dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_f4_rr_radrad(self, lsf4_rr_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rr_radrad(:, :, :, :)

      if (self%n_active == 0) then
         lsf4_rr_radrad = 0.0_wp
         return
      end if
      call self%require_deriv(2, "f4_rr_radrad")

      call radius2_family(self, 2, lsf4_rr_radrad=lsf4_rr_radrad)
   end subroutine lsf_f4_rr_radrad

   !> Mixed nuclear-radius Hessian d^2S / (dR_A dR_b)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_rA_rad d^2S/(dR_A dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_f2_rA_rad(self, lsf2_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear-radius Hessian
      real(wp), intent(out) :: lsf2_rA_rad(:, :, :)

      if (self%n_active == 0) then
         lsf2_rA_rad = 0.0_wp
         return
      end if
      ! The retained *position* still costs an order; the radius does not, so
      ! this matches `f2_rArB` rather than sitting one above it.
      call self%require_deriv(1, "f2_rA_rad")

      call nucrad_family(self, 0, lsf2_rA_rad=lsf2_rA_rad)
   end subroutine lsf_f2_rA_rad

   !> Mixed third derivative d^3S / (dr dR_A dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_rA_rad d^3S/(dr dR_A dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_f3_r_rA_rad(self, lsf3_r_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_rA_rad(:, :, :, :)

      if (self%n_active == 0) then
         lsf3_r_rA_rad = 0.0_wp
         return
      end if
      call self%require_deriv(2, "f3_r_rA_rad")

      call nucrad_family(self, 1, lsf3_r_rA_rad=lsf3_r_rA_rad)
   end subroutine lsf_f3_r_rA_rad

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b) [3, 3, 3, >= n_act, >= n_act]
   subroutine lsf_f4_rr_rA_rad(self, lsf4_rr_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative
      real(wp), intent(out) :: lsf4_rr_rA_rad(:, :, :, :, :)

      if (self%n_active == 0) then
         lsf4_rr_rA_rad = 0.0_wp
         return
      end if
      call self%require_deriv(3, "f4_rr_rA_rad")

      call nucrad_family(self, 2, lsf4_rr_rA_rad=lsf4_rr_rA_rad)
   end subroutine lsf_f4_rr_rA_rad

   !* ================================================================================= *!
   !*                     Direction-contracted nuclear derivatives                      *!
   !* ================================================================================= *!

   !> Accumulate the direction-contracted pseudo-density tensors
   !>
   !> The `tg*` family depends on the caller's displacement field, so `prepare`
   !> cannot own it; this is its atom-plus-pair sweep. It is used only by the
   !> `tangent_*` accessors -- the `hvp_*` ones get the same family for free from
   !> `cfc_*_hvp_eval` and must never call this as well, or every contribution
   !> would be counted twice.
   !>
   !> @param[in]  self    LSF instance
   !> @param[in]  v       Nuclear displacement directions [ndim, ncenters]
   !> @param[in]  level   Highest spatial-derivative order to accumulate (0..3)
   !> @param[out] tg0     Order-0 contracted tensor
   !> @param[out] tg1_r   Order-1 contracted tensor
   !> @param[out] tg2_rr  Order-2 contracted tensor
   !> @param[out] tg3_rrr Order-3 contracted tensor
   subroutine tangent_tensors(self, v, level, tg0, tg1_r, tg2_rr, tg3_rrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: level
      !> Order-0 contracted tensor
      real(wp), intent(out) :: tg0
      !> Order-1 contracted tensor
      real(wp), intent(out) :: tg1_r(ndim)
      !> Order-2 contracted tensor
      real(wp), intent(out) :: tg2_rr(ndim, ndim)
      !> Order-3 contracted tensor
      real(wp), intent(out) :: tg3_rrr(ndim, ndim, ndim)

      !> Active-list indices
      integer :: ia, ib

      tg0 = 0.0_wp
      tg1_r = 0.0_wp
      tg2_rr = 0.0_wp
      tg3_rrr = 0.0_wp

      do ia = 1, self%n_active
         call cfc_atomic_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                      self%param%a1, &
                                      v(:, self%act_atom(ia)), level, &
                                      tg0, tg1_r, tg2_rr, tg3_rrr)
      end do
      do ia = 1, self%n_active
         do ib = ia + 1, self%n_active
            call cfc_pair_tangent_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                       self%act_radius(ia), self%act_radius(ib), &
                                       self%param%a2, self%param%c, &
                                       v(:, self%act_atom(ia)), v(:, self%act_atom(ib)), &
                                       level, tg0, tg1_r, tg2_rr, tg3_rrr)
         end do
      end do
   end subroutine tangent_tensors

   !> Directional nuclear derivative of the level-set value
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . dS/dR_B
   subroutine lsf_tangent_f0(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted value derivative
      real(wp), intent(out) :: res

      !> Contracted pseudo-density tensors and unused kernel outputs
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: d1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(0, "tangent_f0")

      call tangent_tensors(self, v, 0, tg0, tg1_r, tg2_rr, tg3_rrr)
      call cfc_tangent_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            tg0, tg1_r, tg2_rr, tg3_rrr, 0, res, d1, d2, d3)
   end subroutine lsf_tangent_f0

   !> Directional nuclear derivative of the spatial gradient
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of dS/dr [3]
   subroutine lsf_tangent_f1_r(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted gradient derivative
      real(wp), intent(out) :: res(:)

      !> Contracted pseudo-density tensors and kernel outputs
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: t0, t1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(1, "tangent_f1_r")

      call tangent_tensors(self, v, 1, tg0, tg1_r, tg2_rr, tg3_rrr)
      call cfc_tangent_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            tg0, tg1_r, tg2_rr, tg3_rrr, 1, t0, t1, d2, d3)
      res = t1
   end subroutine lsf_tangent_f1_r

   !> Directional nuclear derivative of the spatial Hessian
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^2S/dr^2 [3, 3]
   subroutine lsf_tangent_f2_rr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted Hessian derivative
      real(wp), intent(out) :: res(:, :)

      !> Contracted pseudo-density tensors and kernel outputs
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), d3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(2, "tangent_f2_rr")

      call tangent_tensors(self, v, 2, tg0, tg1_r, tg2_rr, tg3_rrr)
      call cfc_tangent_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            tg0, tg1_r, tg2_rr, tg3_rrr, 2, t0, t1, t2, d3)
      res = t2
   end subroutine lsf_tangent_f2_rr

   !> Directional nuclear derivative of the third spatial derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^3S/dr^3 [3, 3, 3]
   subroutine lsf_tangent_f3_rrr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted third-derivative derivative
      real(wp), intent(out) :: res(:, :, :)

      !> Contracted pseudo-density tensors and kernel outputs
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)

      res = 0.0_wp
      if (self%n_active == 0) return
      call self%require_deriv(3, "tangent_f3_rrr")

      call tangent_tensors(self, v, 3, tg0, tg1_r, tg2_rr, tg3_rrr)
      call cfc_tangent_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            tg0, tg1_r, tg2_rr, tg3_rrr, 3, t0, t1, t2, t3)
      res = t3
   end subroutine lsf_tangent_f3_rrr

   !> Fill the per-atom Hessian-vector-product family (and its tangent ladder)
   !>
   !> One atom sweep and one pair sweep of `cfc_*_hvp_eval`, which fills the
   !> aggregate `tg*` accumulators and the per-atom `hvp*` ones together. The
   !> tangent ladder comes out of the same derivative structure, so it is free
   !> here -- and is exactly why [[tangent_tensors]] must not also run.
   !>
   !> Supplying `vrad` promotes the contraction from `sum_B v_B . d/dR_B` to
   !> `sum_B (v_B . d/dR_B + vr_B d/dRad_B)`: a second kernel sweep adds the
   !> radius half to `tg*` and `hvp*` and fills `rh*` in the same pass. `vrad`
   !> is the switch: the three `rh*` buffers must be present whenever it is, and
   !> supplying them without it is harmless (they come back zeroed). Without
   !> `vrad` nothing changes, down to the bit -- the position-only sweep is
   !> untouched.
   !>
   !> The `hv*` family is optional for the mirror-image reason: [[cfc_radius_hvp]]
   !> wants the joint `tg*` and `rh*` but not the position row, and materializing
   !> `39 * n_active` doubles to throw away is the one part of that trade worth
   !> avoiding. Absent, the row accumulates into stack scratch instead. That is
   !> what doubles the sweeps below into `present`-guarded pairs; a single body
   !> would need an aliasing-opaque pointer with `merge`-computed slot indices,
   !> which is a worse thing to put in the `O(n_active**2)` pair loop than a
   !> duplicated call statement. `hv0`, `hv1` and `hv2` are all-or-none.
   !>
   !> @param[in]  self    LSF instance
   !> @param[in]  v       Nuclear displacement directions [ndim, ncenters]
   !> @param[in]  level   Highest spatial-derivative order to accumulate (0..2)
   !> @param[out] tg0     Order-0 contracted tensor
   !> @param[out] tg1_r   Order-1 contracted tensor
   !> @param[out] tg2_rr  Order-2 contracted tensor
   !> @param[out] tg3_rrr Order-3 contracted tensor (unused above level 2)
   !> @param[out] hv0     Order-0 per-atom HVP tensor [ndim, n_active] (optional)
   !> @param[out] hv1     Order-1 per-atom HVP tensor [ndim, ndim, n_active] (optional)
   !> @param[out] hv2     Order-2 per-atom HVP tensor [ndim, ndim, ndim, n_active]
   !>                     (optional)
   !> @param[in]  vrad    Radius directions [ncenters] (optional; see above)
   !> @param[out] rh0     Order-0 per-atom radius-row tensor [n_active] (optional)
   !> @param[out] rh1_r   Order-1 per-atom radius-row tensor [ndim, n_active] (optional)
   !> @param[out] rh2_rr  Order-2 per-atom radius-row tensor [ndim, ndim, n_active]
   !>                     (optional)
   subroutine hvp_tensors(self, v, level, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2, &
                          vrad, rh0, rh1_r, rh2_rr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: level
      !> Order-0 contracted tensor
      real(wp), intent(out) :: tg0
      !> Order-1 contracted tensor
      real(wp), intent(out) :: tg1_r(ndim)
      !> Order-2 contracted tensor
      real(wp), intent(out) :: tg2_rr(ndim, ndim)
      !> Order-3 contracted tensor
      real(wp), intent(out) :: tg3_rrr(ndim, ndim, ndim)
      !> Order-0 per-atom HVP tensor
      real(wp), intent(out), optional :: hv0(:, :)
      !> Order-1 per-atom HVP tensor
      real(wp), intent(out), optional :: hv1(:, :, :)
      !> Order-2 per-atom HVP tensor
      real(wp), intent(out), optional :: hv2(:, :, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)
      !> Order-0 per-atom radius-row tensor
      real(wp), intent(out), optional :: rh0(:)
      !> Order-1 per-atom radius-row tensor
      real(wp), intent(out), optional :: rh1_r(:, :)
      !> Order-2 per-atom radius-row tensor
      real(wp), intent(out), optional :: rh2_rr(:, :, :)

      !> Position-row sinks, one set per kernel slot so the pair calls do not alias
      real(wp) :: sk_a0(ndim), sk_a1(ndim, ndim), sk_a2(ndim, ndim, ndim)
      real(wp) :: sk_b0(ndim), sk_b1(ndim, ndim), sk_b2(ndim, ndim, ndim)
      !> Whether the caller wants the position row at all
      logical :: want_hv
      !> Active-list indices
      integer :: ia, ib

      want_hv = present(hv0)
      if (want_hv .neqv. (present(hv1) .and. present(hv2))) then
         error stop "moist DROP CFC LSF: hvp_tensors takes hv0, hv1 and hv2 "// &
            "together or not at all"
      end if

      tg0 = 0.0_wp
      tg1_r = 0.0_wp
      tg2_rr = 0.0_wp
      tg3_rrr = 0.0_wp
      if (want_hv) then
         hv0 = 0.0_wp
         hv1 = 0.0_wp
         hv2 = 0.0_wp
      else
         ! Zeroed once, not per atom: the kernel slots are `intent(inout)`
         ! accumulators, so leaving the sinks undefined would have the kernel
         ! read uninitialized stack and trap under the debug FP settings.
         sk_a0 = 0.0_wp
         sk_a1 = 0.0_wp
         sk_a2 = 0.0_wp
         sk_b0 = 0.0_wp
         sk_b1 = 0.0_wp
         sk_b2 = 0.0_wp
      end if
      ! Zeroed before the `vrad` gate, not after it: a caller may hand over the
      ! buffers without asking for the radius channel, and an untouched
      ! `intent(out)` array is the one thing it could not inspect safely.
      if (present(rh0)) rh0 = 0.0_wp
      if (present(rh1_r)) rh1_r = 0.0_wp
      if (present(rh2_rr)) rh2_rr = 0.0_wp

      if (want_hv) then
         do ia = 1, self%n_active
            call cfc_atomic_hvp_eval(self%act_d(:, ia), self%act_radius(ia), &
                                     self%param%a1, &
                                     v(:, self%act_atom(ia)), level, &
                                     tg0, tg1_r, tg2_rr, tg3_rrr, &
                                     hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia))
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                      self%act_radius(ia), self%act_radius(ib), &
                                      self%param%a2, self%param%c, &
                                      v(:, self%act_atom(ia)), v(:, self%act_atom(ib)), &
                                      level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                      hv0(:, ia), hv0(:, ib), &
                                      hv1(:, :, ia), hv1(:, :, ib), &
                                      hv2(:, :, :, ia), hv2(:, :, :, ib))
            end do
         end do
      else
         do ia = 1, self%n_active
            call cfc_atomic_hvp_eval(self%act_d(:, ia), self%act_radius(ia), &
                                     self%param%a1, &
                                     v(:, self%act_atom(ia)), level, &
                                     tg0, tg1_r, tg2_rr, tg3_rrr, &
                                     sk_a0, sk_a1, sk_a2)
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                      self%act_radius(ia), self%act_radius(ib), &
                                      self%param%a2, self%param%c, &
                                      v(:, self%act_atom(ia)), v(:, self%act_atom(ib)), &
                                      level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                      sk_a0, sk_b0, &
                                      sk_a1, sk_b1, &
                                      sk_a2, sk_b2)
            end do
         end do
      end if

      if (.not. present(vrad)) return
      ! `vrad` is the switch, but the loops below write `rh*` unconditionally, so
      ! a caller that turns the radius pass on without handing over all three
      ! buffers would write through an absent dummy. Refuse instead.
      if (.not. (present(rh0) .and. present(rh1_r) .and. present(rh2_rr))) then
         error stop "moist DROP CFC LSF: hvp_tensors needs rh0, rh1_r and "// &
            "rh2_rr whenever vrad is supplied"
      end if

      if (want_hv) then
         do ia = 1, self%n_active
            call cfc_atomic_radius_hvp_eval(self%act_d(:, ia), self%act_radius(ia), &
                                            self%param%a1, &
                                            v(:, self%act_atom(ia)), vrad(self%act_atom(ia)), &
                                            level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                            hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia), &
                                            rh0(ia), rh1_r(:, ia), rh2_rr(:, :, ia))
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_radius_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                             self%act_radius(ia), self%act_radius(ib), &
                                             self%param%a2, self%param%c, &
                                             v(:, self%act_atom(ia)), v(:, self%act_atom(ib)), &
                                             vrad(self%act_atom(ia)), vrad(self%act_atom(ib)), &
                                             level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                             hv0(:, ia), hv0(:, ib), &
                                             hv1(:, :, ia), hv1(:, :, ib), &
                                             hv2(:, :, :, ia), hv2(:, :, :, ib), &
                                             rh0(ia), rh0(ib), &
                                             rh1_r(:, ia), rh1_r(:, ib), &
                                             rh2_rr(:, :, ia), rh2_rr(:, :, ib))
            end do
         end do
      else
         do ia = 1, self%n_active
            call cfc_atomic_radius_hvp_eval(self%act_d(:, ia), self%act_radius(ia), &
                                            self%param%a1, &
                                            v(:, self%act_atom(ia)), vrad(self%act_atom(ia)), &
                                            level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                            sk_a0, sk_a1, sk_a2, &
                                            rh0(ia), rh1_r(:, ia), rh2_rr(:, :, ia))
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_radius_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                             self%act_radius(ia), self%act_radius(ib), &
                                             self%param%a2, self%param%c, &
                                             v(:, self%act_atom(ia)), v(:, self%act_atom(ib)), &
                                             vrad(self%act_atom(ia)), vrad(self%act_atom(ib)), &
                                             level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                             sk_a0, sk_b0, &
                                             sk_a1, sk_b1, &
                                             sk_a2, sk_b2, &
                                             rh0(ia), rh0(ib), &
                                             rh1_r(:, ia), rh1_r(:, ib), &
                                             rh2_rr(:, :, ia), rh2_rr(:, :, ib))
            end do
         end do
      end if
   end subroutine hvp_tensors

   !> Nuclear Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters] (optional); supplying it
   !>                  promotes the contraction to the joint position/radius one
   !> @param[out] res  sum_B v_B . d^2S/(dR_A dR_B) [3, >= active_count()]
   subroutine lsf_hvp_f1_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted nuclear Hessian
      real(wp), intent(out) :: res(:, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors and the radius-row by-product
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      real(wp), allocatable :: rh0(:), rh1_r(:, :), rh2_rr(:, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         res = 0.0_wp
         return
      end if
      call self%require_deriv(1, "hvp_f1_rA")

      allocate (hv0(ndim, self%n_active))
      allocate (hv1(ndim, ndim, self%n_active))
      allocate (hv2(ndim, ndim, ndim, self%n_active))
      if (present(vrad)) then
         call alloc_radius_row(self%n_active, rh0, rh1_r, rh2_rr)
         call hvp_tensors(self, v, 0, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2, &
                          vrad, rh0, rh1_r, rh2_rr)
      else
         call hvp_tensors(self, v, 0, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)
      end if

      do ia = 1, self%n_active
         call cfc_hvp_eval(self%pd0, self%pd1_r, self%pd2_rr, tg0, tg1_r, tg2_rr, &
                           self%qn0(:, ia), self%qn1_r(:, :, ia), self%qn2_rr(:, :, :, ia), &
                           hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia), 0, h1, d2, d3)
         res(:, ia) = h1
      end do
   end subroutine lsf_hvp_f1_rA

   !> Directional nuclear derivative of `f2_r_rA`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters] (optional)
   !> @param[out] res  sum_B v_B . d^3S/(dr dR_A dR_B) [3, 3, >= active_count()]
   subroutine lsf_hvp_f2_r_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors and the radius-row by-product
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      real(wp), allocatable :: rh0(:), rh1_r(:, :), rh2_rr(:, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         res = 0.0_wp
         return
      end if
      call self%require_deriv(2, "hvp_f2_r_rA")

      allocate (hv0(ndim, self%n_active))
      allocate (hv1(ndim, ndim, self%n_active))
      allocate (hv2(ndim, ndim, ndim, self%n_active))
      if (present(vrad)) then
         call alloc_radius_row(self%n_active, rh0, rh1_r, rh2_rr)
         call hvp_tensors(self, v, 1, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2, &
                          vrad, rh0, rh1_r, rh2_rr)
      else
         call hvp_tensors(self, v, 1, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)
      end if

      do ia = 1, self%n_active
         call cfc_hvp_eval(self%pd0, self%pd1_r, self%pd2_rr, tg0, tg1_r, tg2_rr, &
                           self%qn0(:, ia), self%qn1_r(:, :, ia), self%qn2_rr(:, :, :, ia), &
                           hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia), 1, h1, h2, d3)
         res(:, :, ia) = h2
      end do
   end subroutine lsf_hvp_f2_r_rA

   !> Directional nuclear derivative of `f3_rr_rA`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters] (optional)
   !> @param[out] res  sum_B v_B . d^4S/(dr^2 dR_A dR_B) [3, 3, 3, >= active_count()]
   subroutine lsf_hvp_f3_rr_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: res(:, :, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors and the radius-row by-product
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      real(wp), allocatable :: rh0(:), rh1_r(:, :), rh2_rr(:, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), h3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         res = 0.0_wp
         return
      end if
      call self%require_deriv(3, "hvp_f3_rr_rA")

      allocate (hv0(ndim, self%n_active))
      allocate (hv1(ndim, ndim, self%n_active))
      allocate (hv2(ndim, ndim, ndim, self%n_active))
      if (present(vrad)) then
         call alloc_radius_row(self%n_active, rh0, rh1_r, rh2_rr)
         call hvp_tensors(self, v, 2, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2, &
                          vrad, rh0, rh1_r, rh2_rr)
      else
         call hvp_tensors(self, v, 2, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)
      end if

      do ia = 1, self%n_active
         call cfc_hvp_eval(self%pd0, self%pd1_r, self%pd2_rr, tg0, tg1_r, tg2_rr, &
                           self%qn0(:, ia), self%qn1_r(:, :, ia), self%qn2_rr(:, :, :, ia), &
                           hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia), 2, h1, h2, h3)
         res(:, :, :, ia) = h3
      end do
   end subroutine lsf_hvp_f3_rr_rA

   !* ================================================================================= *!
   !*                     Radius row of the joint Hessian-vector product                *!
   !* ================================================================================= *!
   !
   ! The `hvp_*_rA` procedures above retain a nuclear index; these three retain a
   ! *radius* index, so between them they cover all four blocks of the joint
   ! position/radius Hessian-vector product. Both rows read the same `tg*`
   ! aggregate, because the joint direction lands in a single family: nothing in
   ! the generated lift distinguishes its two halves.
   !
   ! The radius-radius coupling between *different* atoms is nonzero and is not
   ! recoverable from per-atom data: `rd*` carries no dependence on another
   ! atom's radius, so all of it arrives through the pair term and through the
   ! log, exactly as for `f2_rArB`.

   !> Allocate the per-atom radius-row buffers
   !>
   !> Only ever called on a path that supplies a radius direction, because the
   !> rows exist only for a joint `(v, vrad)` contraction: live output for
   !> [[cfc_radius_hvp]], write-only scratch for the three `hvp_*_rA` entry
   !> points, which need `hvp_tensors` to have somewhere to put the row it fills
   !> alongside the `tg*` they do read. A pure-nuclear contraction never runs the
   !> radius sweep at all, so it must not pay for these -- hence the `present`
   !> guard at those three call sites rather than an unconditional allocation.
   !>
   !> @param[in]  n      Active-atom count
   !> @param[out] rh0    Order-0 radius-row buffer [n]
   !> @param[out] rh1_r  Order-1 radius-row buffer [ndim, n]
   !> @param[out] rh2_rr Order-2 radius-row buffer [ndim, ndim, n]
   subroutine alloc_radius_row(n, rh0, rh1_r, rh2_rr)
      !> Active-atom count
      integer, intent(in) :: n
      !> Order-0 radius-row buffer
      real(wp), allocatable, intent(out) :: rh0(:)
      !> Order-1 radius-row buffer
      real(wp), allocatable, intent(out) :: rh1_r(:, :)
      !> Order-2 radius-row buffer
      real(wp), allocatable, intent(out) :: rh2_rr(:, :, :)

      allocate (rh0(n))
      allocate (rh1_r(ndim, n))
      allocate (rh2_rr(ndim, ndim, n))
   end subroutine alloc_radius_row

   !> Radius row of the joint Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  sum_B (v_B . d/dR_B + vr_B d/dR_b) dS/dR_a
   !>                  [>= active_count()]
   subroutine lsf_hvp_f1_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted radius Hessian row
      real(wp), intent(out) :: res(:)

      call cfc_radius_hvp(self, v, vrad, 0, "hvp_f1_rad", res1=res)
   end subroutine lsf_hvp_f1_rad

   !> Joint directional derivative of `f2_r_rad`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed third derivative [3, >= active_count()]
   subroutine lsf_hvp_f2_r_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :)

      call cfc_radius_hvp(self, v, vrad, 1, "hvp_f2_r_rad", res2=res)
   end subroutine lsf_hvp_f2_r_rad

   !> Joint directional derivative of `f3_rr_rad`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed fourth derivative [3, 3, >= active_count()]
   subroutine lsf_hvp_f3_rr_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: res(:, :, :)

      call cfc_radius_hvp(self, v, vrad, 2, "hvp_f3_rr_rad", res3=res)
   end subroutine lsf_hvp_f3_rr_rad

   !> Shared driver of the three radius-row Hessian-vector products
   !>
   !> The three public entry points differ only in `max_deriv` and in which
   !> output they keep, and unlike the `_rA` ladder there is no reason to spell
   !> the loop out three times: the radius channel has one code path. Exactly one
   !> of `res1`/`res2`/`res3` must be present, matching `max_deriv`.
   !>
   !> Two layer-1 sweeps feed it. `hvp_tensors` with a radius direction gives the
   !> joint `tg*` and the `rh*` family; `radius_tensors` gives `rd*` up to order
   !> two. Both also produce a row this driver has no use for -- `hvp*`, which
   !> belongs to the *position* row, and `rd3_rrr` -- and both are asked to sink
   !> it into stack scratch by omitting the buffer. The work itself is not
   !> avoidable and not worth avoiding: the generated kernel produces the two
   !> rows' radius terms from one shared CSE block, cheaper than computing them
   !> apart even counting the waste (6111 temporaries against 4427 + 2230). Only
   !> the `66 * n_active` doubles of throwaway storage were.
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  v         Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad      Radius directions [ncenters]
   !> @param[in]  max_deriv Highest spatial-derivative order (0..2)
   !> @param[in]  caller    Name used by `require_deriv` on failure
   !> @param[out] res1      Order-0 result [>= active_count()] (optional)
   !> @param[out] res2      Order-1 result [3, >= active_count()] (optional)
   !> @param[out] res3      Order-2 result [3, 3, >= active_count()] (optional)
   subroutine cfc_radius_hvp(self, v, vrad, max_deriv, caller, res1, res2, res3)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Highest spatial-derivative order
      integer, intent(in) :: max_deriv
      !> Caller name for the derivative-order check
      character(len=*), intent(in) :: caller
      !> Order-0 result
      real(wp), intent(out), optional :: res1(:)
      !> Order-1 result
      real(wp), intent(out), optional :: res2(:, :)
      !> Order-2 result
      real(wp), intent(out), optional :: res3(:, :, :)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom radius and radius-row tensors
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rh0(:), rh1_r(:, :), rh2_rr(:, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1, h2(ndim), h3(ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) then
         ! Nothing is active, so no slot is owned and "writes the first
         ! `active_count()` slots" degenerates to writing none. Zero rather than
         ! return bare: the results are `intent(out)`, so a bare return hands the
         ! caller an undefined buffer.
         if (present(res1)) res1 = 0.0_wp
         if (present(res2)) res2 = 0.0_wp
         if (present(res3)) res3 = 0.0_wp
         return
      end if
      ! One order less than the `hvp_*_rA` entry point of the same spatial
      ! order, for two independent reasons: `d/dRad_a` consumes no derivative
      ! order, and the radius row reads no cached `qn*` (whose spatial ladder
      ! `prepare` fills one short of `max_deriv`).
      call self%require_deriv(max_deriv, caller)

      allocate (rd0(self%n_active))
      allocate (rd1_r(ndim, self%n_active))
      allocate (rd2_rr(ndim, ndim, self%n_active))
      call alloc_radius_row(self%n_active, rh0, rh1_r, rh2_rr)

      call hvp_tensors(self, v, max_deriv, tg0, tg1_r, tg2_rr, tg3_rrr, &
                       vrad=vrad, rh0=rh0, rh1_r=rh1_r, rh2_rr=rh2_rr)
      call radius_tensors(self, max_deriv, rd0, rd1_r, rd2_rr)

      do ia = 1, self%n_active
         call cfc_radius_hvp_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                  tg0, tg1_r, tg2_rr, &
                                  rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                  rh0(ia), rh1_r(:, ia), rh2_rr(:, :, ia), &
                                  max_deriv, h1, h2, h3)
         if (present(res1)) res1(ia) = h1
         if (present(res2)) res2(:, ia) = h2
         if (present(res3)) res3(:, :, ia) = h3
      end do
   end subroutine cfc_radius_hvp

   !* ================================================================================= *!
   !*                  Jet-contracted nuclear vector-Jacobian product                   *!
   !* ================================================================================= *!
   !
   ! The reverse-mode mirror of the `hvp_*_rA` block above. Those contract the
   ! *nuclear* index of the nuclear Hessian against a displacement field and keep
   ! the spatial ones; this one contracts the *spatial* (jet) indices of the
   ! `f1_rA` / `f2_r_rA` / `f3_rr_rA` ladder against per-point adjoint weights
   ! and keeps the nuclear one.
   !
   ! It reads exactly what `f3_rr_rA` reads -- the `pd*` aggregate and the `qn*`
   ! per-atom family, both cached by `prepare` -- so it needs no sweep of its own
   ! and requires the same prepared order. What it saves is in the kernel, not
   ! here: `cfc_vjp_eval` folds the weights in before the generated code takes
   ! its common subexpressions, so three numbers per atom come out where the
   ! uncontracted accessor produces 3 + 9 + 27, and the caller's
   ! `[3, 3, 3, n_active]` buffer disappears.

   !> Adjoint jet contracted onto the nuclear gradient, active-indexed
   !>
   !> Returns, for every active atom `i`,
   !>
   !>    res(beta, i) = w0 f1_rA(beta, i)
   !>                 + sum_a w1(a) f2_r_rA(a, beta, i)
   !>                 + sum_a sum_b w2(a, b) f3_rr_rA(a, b, beta, i)
   !>
   !> i.e. the `f1_rA` rung with the jet indices summed away, just as
   !> `hvp_f1_rA` is that same rung contracted with a nuclear direction instead.
   !> `w2` is a general 3x3: all nine entries are contracted, with no symmetry
   !> assumption and no folded factor of two.
   !>
   !> Slot `i` belongs to `active_atom(i)`; columns past `active_count()` are
   !> left untouched and nothing is written at all when no atom is active. The
   !> order requirement is 3, not 2, for the same reason as `f3_rr_rA`: the
   !> kernel reads `qn2_rr`, whose total derivative order is three, and `prepare`
   !> fills that family only one spatial order below `max_deriv`.
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of `f1_rA`
   !> @param[in]  w1   Adjoint weights of `f2_r_rA` [3]
   !> @param[in]  w2   Adjoint weights of `f3_rr_rA` [3, 3]
   !> @param[out] res  Contracted nuclear gradient [3, >= active_count()]
   subroutine lsf_vjp_f1_rA(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Adjoint weight of `f1_rA`
      real(wp), intent(in) :: w0
      !> Adjoint weights of `f2_r_rA`
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of `f3_rr_rA`
      real(wp), intent(in) :: w2(3, 3)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: res(:, :)

      !> Kernel output of one atom
      real(wp) :: vjp1(ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(3, "vjp_f1_rA")

      do ia = 1, self%n_active
         call cfc_vjp_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                           self%qn0(:, ia), self%qn1_r(:, :, ia), &
                           self%qn2_rr(:, :, :, ia), w0, w1, w2, vjp1)
         res(:, ia) = vjp1
      end do
   end subroutine lsf_vjp_f1_rA

   !* ================================================================================= *!
   !*                   Jet-contracted radius vector-Jacobian product                   *!
   !* ================================================================================= *!

   !> Adjoint jet contracted onto the radius gradient, active-indexed
   !>
   !> Returns, for every active atom `i`,
   !>
   !>    res(i) = w0 f1_rad(i)
   !>           + sum_a w1(a) f2_r_rad(a, i)
   !>           + sum_a sum_b w2(a, b) f3_rr_rad(a, b, i)
   !>
   !> i.e. the `f1_rad` rung with the jet indices summed away, the radius twin of
   !> `vjp_f1_rA`. The result is a scalar per atom and not a Cartesian row,
   !> because `R_a` is a radius; the jet itself is unchanged. `w2` is a general
   !> 3x3: all nine entries are contracted, with no symmetry assumption and no
   !> folded factor of two.
   !>
   !> Slot `i` belongs to `active_atom(i)`; entries past `active_count()` are
   !> left untouched and nothing is written at all when no atom is active. The
   !> order requirement is 2, one less than `vjp_f1_rA` needs, for the same
   !> reason `f3_rr_rad` needs one less than `f3_rr_rA`: `d/dR_a` consumes no
   !> derivative order, so the deepest input is `rd2_rr` at total order two.
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of `f1_rad`
   !> @param[in]  w1   Adjoint weights of `f2_r_rad` [3]
   !> @param[in]  w2   Adjoint weights of `f3_rr_rad` [3, 3]
   !> @param[out] res  Contracted radius gradient [>= active_count()]
   subroutine lsf_vjp_f1_rad(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Adjoint weight of `f1_rad`
      real(wp), intent(in) :: w0
      !> Adjoint weights of `f2_r_rad`
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of `f3_rr_rad`
      real(wp), intent(in) :: w2(3, 3)
      !> Contracted radius gradient
      real(wp), intent(out) :: res(:)

      !> Per-atom radius tensors
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      !> Kernel output of one atom
      real(wp) :: vjp1
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(2, "vjp_f1_rad")

      allocate (rd0(self%n_active))
      allocate (rd1_r(ndim, self%n_active))
      allocate (rd2_rr(ndim, ndim, self%n_active))
      call radius_tensors(self, 2, rd0, rd1_r, rd2_rr)

      do ia = 1, self%n_active
         call cfc_radius_vjp_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                  rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                  w0, w1, w2, vjp1)
         res(ia) = vjp1
      end do
   end subroutine lsf_vjp_f1_rad

   !* ================================================================================= *!
   !*                                 Screening offset                                  *!
   !* ================================================================================= *!

   !> Radial offset from the atom surface where the CFC contribution hits the threshold
   !>
   !> The pseudo-density gets two contributions from atom `a`, both decaying in the
   !> reduced distance `s_a = 1 + delta/R_a`, so the offset scales with the radius.
   !>
   !> @param[in] self    LSF instance
   !> @param[in] radius  Atom radius (Bohr)
   !> @returns           Radial offset from the atom surface (Bohr)
   pure function lsf_screening_offset(self, radius) result(offset)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Atom radius (Bohr)
      real(wp), intent(in) :: radius
      !> Radial offset from the atom surface (Bohr)
      real(wp) :: offset

      !> Threshold, its logarithm, the two branch offsets and the pair envelope cap
      real(wp) :: thr, log_thr, delta_atomic, delta_pair, cap
      !> `2**m` of the pair envelope
      real(wp) :: two_pow_m

      thr = self%screening_threshold
      if (thr <= 0.0_wp .or. self%param%a1 >= 0.0_wp .or. self%param%a2 >= 0.0_wp) then
         offset = huge(0.0_wp)
         return
      end if
      ! The pair envelope needs m <= |a2|; without it no closed form is available
      if (real(self%param%m, wp) > abs(self%param%a2)) then
         offset = huge(0.0_wp)
         return
      end if

      log_thr = log(thr)
      delta_atomic = -log_thr/abs(self%param%a1)*radius

      two_pow_m = 2.0_wp**self%param%m
      cap = self%param%c*two_pow_m
      delta_pair = 2.0_wp*radius*(log_thr - log(cap))/self%param%a2

      offset = max(0.0_wp, delta_atomic, delta_pair)
   end function lsf_screening_offset

   !* ================================================================================= *!
   !*                                    Finalizer                                      *!
   !* ================================================================================= *!

   !> Release the allocatable components of a CFC LSF instance
   !>
   !> @param[inout] self LSF instance being finalized
   subroutine finalize_lsf_cfc(self)
      !> LSF instance
      type(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self

      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%active_cand)) deallocate (self%active_cand)
      if (allocated(self%act_atom)) deallocate (self%act_atom)
      if (allocated(self%act_d)) deallocate (self%act_d)
      if (allocated(self%act_radius)) deallocate (self%act_radius)
      if (allocated(self%qn0)) deallocate (self%qn0)
      if (allocated(self%qn1_r)) deallocate (self%qn1_r)
      if (allocated(self%qn2_rr)) deallocate (self%qn2_rr)
      if (allocated(self%qn3_rrr)) deallocate (self%qn3_rrr)
      if (allocated(self%cand_screen)) deallocate (self%cand_screen)
      if (allocated(self%cand_radii)) deallocate (self%cand_radii)
      if (allocated(self%cand_to_user)) deallocate (self%cand_to_user)
      if (allocated(self%full_scan_cand)) deallocate (self%full_scan_cand)
      if (allocated(self%orig_to_sorted)) deallocate (self%orig_to_sorted)
      if (allocated(self%sorted_to_orig)) deallocate (self%sorted_to_orig)

      ! Deallocate structure_type allocatable components
      if (allocated(self%mol%id)) deallocate (self%mol%id)
      if (allocated(self%mol%num)) deallocate (self%mol%num)
      if (allocated(self%mol%sym)) deallocate (self%mol%sym)
      if (allocated(self%mol%xyz)) deallocate (self%mol%xyz)
      if (allocated(self%mol%lattice)) deallocate (self%mol%lattice)
      if (allocated(self%mol%periodic)) deallocate (self%mol%periodic)
      if (allocated(self%mol%bond)) deallocate (self%mol%bond)
      if (allocated(self%mol%comment)) deallocate (self%mol%comment)
      if (allocated(self%mol%sdf)) deallocate (self%mol%sdf)
      if (allocated(self%mol%pdb)) deallocate (self%mol%pdb)
   end subroutine finalize_lsf_cfc

end module moist_cavity_drop_lsf_cfc
