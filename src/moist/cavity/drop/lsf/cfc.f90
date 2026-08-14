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
!> requires 3. Anything an accessor builds itself (`tg*`, `hvp*`, `qq*`) does not
!> enter that count. An accessor asked for an order `prepare` did not provision
!> aborts; none of them returns zeros, which a caller could not tell from a real
!> answer.
!>
module moist_cavity_drop_lsf_cfc
   use mctc_env, only: error_type
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_sorted
   use moist_cavity_drop_lsf_cfc_kernel, only: &
      cfc_atomic_term_eval, cfc_pair_spatial_eval, &
      cfc_atomic_nuclear_eval, cfc_pair_nuclear_eval, &
      cfc_atomic_hessian_eval, cfc_pair_hessian_eval, &
      cfc_atomic_tangent_eval, cfc_pair_tangent_eval, &
      cfc_atomic_hvp_eval, cfc_pair_hvp_eval, &
      cfc_spatial_eval, cfc_nuclear_eval, cfc_hessian_eval, &
      cfc_tangent_eval, cfc_hvp_eval
   implicit none (type, external)
   private

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Diedenhofen-Klamt 2018 defaults for the four CFC shape parameters
   real(wp), parameter :: a1_default = -15.0_wp
   real(wp), parameter :: a2_default = -9.0_wp
   real(wp), parameter :: c_default = 5.0_wp
   integer, parameter :: m_default = 4

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
      !> Atomic-term exponent (Diedenhofen-Klamt 2018: -15)
      real(wp) :: a1 = a1_default
      !> Pair-term exponent (Diedenhofen-Klamt 2018: -9)
      real(wp) :: a2 = a2_default
      !> Pair-term coupling constant (Diedenhofen-Klamt 2018: 5)
      real(wp) :: c = c_default
      !> Pair-term polynomial power (Diedenhofen-Klamt 2018: 4)
      !>
      !> The kernel bakes `m = 4` into its symbolic differentiation, so this
      !> field only feeds [[lsf_screening_offset]]; a value other than 4 is
      !> inconsistent with the derivatives and is not supported
      integer :: m = m_default

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

      if (present(a1)) self%a1 = a1
      if (present(a2)) self%a2 = a2
      if (present(c)) self%c = c
      if (present(m)) self%m = m
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
         call cfc_atomic_term_eval(self%act_d(:, i), self%act_radius(i), self%a1, md, &
                                   self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                                   self%pd4_rrrr)
         if (nd < 0) cycle
         call cfc_atomic_nuclear_eval(self%act_d(:, i), self%act_radius(i), self%a1, nd, &
                                      self%qn0(:, i), self%qn1_r(:, :, i), &
                                      self%qn2_rr(:, :, :, i), self%qn3_rrr(:, :, :, :, i))
      end do

      do i = 1, n
         do j = i + 1, n
            call cfc_pair_spatial_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%a2, self%c, md, &
                                       self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                                       self%pd4_rrrr)
            if (nd < 0) cycle
            call cfc_pair_nuclear_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%a2, self%c, nd, &
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
         call cfc_atomic_hessian_eval(self%act_d(:, ia), self%act_radius(ia), self%a1, &
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
                                       self%a2, self%c, level, &
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
                                          self%a2, self%c, level, &
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
         call cfc_atomic_tangent_eval(self%act_d(:, ia), self%act_radius(ia), self%a1, &
                                      v(:, self%act_atom(ia)), level, &
                                      tg0, tg1_r, tg2_rr, tg3_rrr)
      end do
      do ia = 1, self%n_active
         do ib = ia + 1, self%n_active
            call cfc_pair_tangent_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                       self%act_radius(ia), self%act_radius(ib), &
                                       self%a2, self%c, &
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
   !> @param[in]  self    LSF instance
   !> @param[in]  v       Nuclear displacement directions [ndim, ncenters]
   !> @param[in]  level   Highest spatial-derivative order to accumulate (0..2)
   !> @param[out] tg0     Order-0 contracted tensor
   !> @param[out] tg1_r   Order-1 contracted tensor
   !> @param[out] tg2_rr  Order-2 contracted tensor
   !> @param[out] tg3_rrr Order-3 contracted tensor (unused above level 2)
   !> @param[out] hv0     Order-0 per-atom HVP tensor [ndim, n_active]
   !> @param[out] hv1     Order-1 per-atom HVP tensor [ndim, ndim, n_active]
   !> @param[out] hv2     Order-2 per-atom HVP tensor [ndim, ndim, ndim, n_active]
   subroutine hvp_tensors(self, v, level, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)
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
      real(wp), intent(out) :: hv0(:, :)
      !> Order-1 per-atom HVP tensor
      real(wp), intent(out) :: hv1(:, :, :)
      !> Order-2 per-atom HVP tensor
      real(wp), intent(out) :: hv2(:, :, :, :)

      !> Active-list indices
      integer :: ia, ib

      tg0 = 0.0_wp
      tg1_r = 0.0_wp
      tg2_rr = 0.0_wp
      tg3_rrr = 0.0_wp
      hv0 = 0.0_wp
      hv1 = 0.0_wp
      hv2 = 0.0_wp

      do ia = 1, self%n_active
         call cfc_atomic_hvp_eval(self%act_d(:, ia), self%act_radius(ia), self%a1, &
                                  v(:, self%act_atom(ia)), level, &
                                  tg0, tg1_r, tg2_rr, tg3_rrr, &
                                  hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia))
      end do
      do ia = 1, self%n_active
         do ib = ia + 1, self%n_active
            call cfc_pair_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                   self%act_radius(ia), self%act_radius(ib), &
                                   self%a2, self%c, &
                                   v(:, self%act_atom(ia)), v(:, self%act_atom(ib)), &
                                   level, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                   hv0(:, ia), hv0(:, ib), &
                                   hv1(:, :, ia), hv1(:, :, ib), &
                                   hv2(:, :, :, ia), hv2(:, :, :, ib))
         end do
      end do
   end subroutine hvp_tensors

   !> Nuclear Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^2S/(dR_A dR_B) [3, >= active_count()]
   subroutine lsf_hvp_f1_rA(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted nuclear Hessian
      real(wp), intent(out) :: res(:, :)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(1, "hvp_f1_rA")

      allocate (hv0(ndim, self%n_active))
      allocate (hv1(ndim, ndim, self%n_active))
      allocate (hv2(ndim, ndim, ndim, self%n_active))
      call hvp_tensors(self, v, 0, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)

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
   !> @param[out] res  sum_B v_B . d^3S/(dr dR_A dR_B) [3, 3, >= active_count()]
   subroutine lsf_hvp_f2_r_rA(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :, :)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(2, "hvp_f2_r_rA")

      allocate (hv0(ndim, self%n_active))
      allocate (hv1(ndim, ndim, self%n_active))
      allocate (hv2(ndim, ndim, ndim, self%n_active))
      call hvp_tensors(self, v, 1, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)

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
   !> @param[out] res  sum_B v_B . d^4S/(dr^2 dR_A dR_B) [3, 3, 3, >= active_count()]
   subroutine lsf_hvp_f3_rr_rA(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed fourth derivative
      real(wp), intent(out) :: res(:, :, :, :)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), h3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      if (self%n_active == 0) return
      call self%require_deriv(3, "hvp_f3_rr_rA")

      allocate (hv0(ndim, self%n_active))
      allocate (hv1(ndim, ndim, self%n_active))
      allocate (hv2(ndim, ndim, ndim, self%n_active))
      call hvp_tensors(self, v, 2, tg0, tg1_r, tg2_rr, tg3_rrr, hv0, hv1, hv2)

      do ia = 1, self%n_active
         call cfc_hvp_eval(self%pd0, self%pd1_r, self%pd2_rr, tg0, tg1_r, tg2_rr, &
                           self%qn0(:, ia), self%qn1_r(:, :, ia), self%qn2_rr(:, :, :, ia), &
                           hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia), 2, h1, h2, h3)
         res(:, :, :, ia) = h3
      end do
   end subroutine lsf_hvp_f3_rr_rA

   !* ================================================================================= *!
   !*                                 Screening offset                                  *!
   !* ================================================================================= *!

   !> Radial offset from the atom surface where the CFC contribution hits the threshold
   !>
   !> The pseudo-density gets two contributions from atom `a`, both decaying in the
   !> reduced distance `s_a = 1 + delta/R_a`, so the offset scales with the radius.
   !>
   !> **Atomic term** -- exact:
   !>     `exp(a1 (s_a - 1)) = thr`   with `a1 < 0`
   !>     => `delta_atomic = -log(thr) / |a1| * R`
   !>
   !> **Pair term** -- the pair contribution of `a` with a partner `b` is
   !>     `c (1 - vec(s_a).vec(s_b))^m exp(a2 (s_a + s_b - 2))`.
   !>     `vec(s_a)` has length `s_a`, *not* one, so
   !>     `(1 - vec(s_a).vec(s_b))^m <= (1 + s_a s_b)^m` -- a polynomial in `s_a`
   !>     that grows without bound. Fixing the partner at its own surface
   !>     (`s_b = 1`, the dominant configuration for a point near the cavity) leaves
   !>     `c (1 + s_a)^m exp(a2 (s_a - 1))`, which has no closed-form inverse.
   !>     It does have a clean closed-form envelope: for `m <= |a2|` and `s >= 1`,
   !>     `(1 + s)^m exp(a2 (s-1)/2) <= 2^m` (the log-derivative
   !>     `m/(1+s) + a2/2` is <= 0 there), hence
   !>         `c (1 + s)^m exp(a2 (s-1)) <= c 2^m exp((a2/2)(s-1))`
   !>     and setting the right-hand side to `thr` gives
   !>         `delta_pair = 2 R (log(thr) - log(c 2^m)) / a2`
   !>     -- exactly twice the offset this routine's ancestor used. That ancestor
   !>     bounded `(1 - vec(s_a).vec(s_b))^m` by `2^m`, which is only valid while
   !>     `s_a <= 1`, i.e. never in the screening regime; measured against an
   !>     unscreened reference it under-bounded the dropped pseudo-density mass by
   !>     two to three orders of magnitude. When `m > |a2|` the envelope argument
   !>     fails and the routine disables screening rather than guess.
   !>
   !> The larger of the two branches is the offset. It is the single source of
   !> truth for both the per-point gate and the cell-grid reach; the per-point gate
   !> used to borrow SvdW's radius-independent `exp(-(k/3)(x-R)) >= thr` criterion
   !> with `k = 3`, which is unrelated to CFC's own decay and was the second half of
   !> the screening defect fixed here.
   !>
   !> Caveat, deliberately left standing: the `s_b = 1` choice is an assumption, not
   !> a bound. A partner whose interior contains the evaluation point has `s_b < 1`
   !> and enhances the pair term by up to `exp(|a2|)`. Covering that case rigorously
   !> costs another `2R` of reach for no accuracy where it matters (deep inside the
   !> cavity the pseudo-density is large, so the relative error is negligible).
   !>
   !> If `threshold <= 0` or any exponent is non-negative the routine
   !> returns `huge(0.0_wp)`, i.e. screening is disabled.
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
      if (thr <= 0.0_wp .or. self%a1 >= 0.0_wp .or. self%a2 >= 0.0_wp) then
         offset = huge(0.0_wp)
         return
      end if
      ! The pair envelope needs m <= |a2|; without it no closed form is available
      if (real(self%m, wp) > abs(self%a2)) then
         offset = huge(0.0_wp)
         return
      end if

      log_thr = log(thr)
      delta_atomic = -log_thr/abs(self%a1)*radius

      two_pow_m = 2.0_wp**self%m
      cap = self%c*two_pow_m
      delta_pair = 2.0_wp*radius*(log_thr - log(cap))/self%a2

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
