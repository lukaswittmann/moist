!> COSMO Fine Cavity (CFC) level set function
!>
!> Two evaluation paths over the generated [[moist_cavity_drop_lsf_cfc_kernel]]:
!>
!> - per point: `prepare` / `prepare_subset` run the base screening gate, accumulate `pd*` /
!>   `qn*` per atom and pair, lift the spatial jet once at `max_deriv`; spatial accessors hand
!>   out that dense jet packed
!> - batch: [[lsf_evaluate_batch]] runs stage 1 `cfc_strip_powersums_o{L}` and stage 2
!>   `cfc_strip_lift_o{L}` per strip of `cfc_strip_width(order)` points over one shared atom
!>   list; spatial jet only, per-atom families (`qn*`, tangent, hvp, radius) stay per point
!> - nuclear vjp: [[lsf_vjp_f1_rA_batch]] runs stage 1 at order 2 with geometry cache, the
!>   order-2 lift, `cfc_strip_lift_vjp` and `cfc_strip_vjp` per strip of `cfc_strip_width(2)`
!>   points; no `qn*` formed
!>
!> Fade
!> ----
!> - default on, 27/3 e-folds
!> - atomic term `w_A e_A`, pair term `w_A w_B c (1 - q)**4 u_A u_B`, `w = S_4(1 - t)` over a
!>   band in e-folds of `u_A = exp(a2 sigma_A)`, `sigma_A = ||r - R_A||/R_a - 1`, each atom's
!>   own radius as unit ([[fade_band]])
!> - [[lsf_screening_offset]] returns the band end: gate mask and fade agree,
!>   `screening_threshold` changes nothing
!> - every per-atom family faded alike, band cached at `prepare` (`act_fade_start`,
!>   `act_fade_width`)
!> - certificates from [[certificate_setup]]: [[lsf_exclusion_radius]] `= |tau(S)|`,
!>   [[lsf_far_exclusion_radius]] for points with no atom in reach
!>
!> Batch order
!> -----------
!> - Morton Z-order sort of the Newton active list (2^-4 Bohr cells): full update 8-13%
!>   faster at projection levels 1-3 (mb16_43_01, upu23, nleb 302); SvdW only paid the
!>   ~11 ns/point sort
!> - projection does not sort; sort per LSF if revisited
!>
!> Not handled
!> -----------
!> - point on a nucleus: both kernels divide by `||r - R_A||`, finite value but non-finite
!>   derivatives (cusp of the atom's term)
!>   - batch: cusp kept in the nucleus lane, sub-strip lifted as is, column unreached
!>     (status 1, NaN) from order 1
!>   - per point: `prepare` reports unreached from `max_deriv = 1`
!>   - batched nuclear vjp: column unreached (status 1, NaN), strip-mates re-run without it
module moist_cavity_drop_lsf_cfc
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
                                              ieee_quiet_nan
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_sorted, &
                                         lsf_cand_bound_sq, lsf_batch_require_outputs, &
                                         lsf_sorted_candidates, lsf_screen_records, &
                                         lsf_vjp_batch_shape_message
   use moist_math_packed_sym, only: pack_sym_matrix, pack_sym_tensor3, pack_sym_tensor4, &
                                    unpack_sym_matrix
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
      cfc_radius_vjp_eval, &
      cfc_strip_width, cfc_npack, cfc_strip_geo_len, &
      cfc_strip_powersums_o0, cfc_strip_lift_o0, &
      cfc_strip_powersums_o1, cfc_strip_lift_o1, &
      cfc_strip_powersums_o2, cfc_strip_lift_o2, &
      cfc_strip_powersums_o3, cfc_strip_lift_o3, &
      cfc_strip_powersums_o4, cfc_strip_lift_o4, &
      cfc_strip_lift_vjp, cfc_strip_vjp
   implicit none(type, external)
   private

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Packed components of the rank-2, -3 and -4 spatial tensors
   integer, parameter :: np2 = cfc_npack(2), np3 = cfc_npack(3), np4 = cfc_npack(4)
   !> Widest strip of any order
   integer, parameter :: max_strip = maxval(cfc_strip_width)
   !> Strip width of the batched nuclear vjp, the order-2 width
   integer, parameter :: vjp_strip = cfc_strip_width(2)

   !> Fade band at infinity: the untruncated model, every atom on the plain body
   !>
   !> - `fade_start + fade_width` must stay finite, hence a finite width
   !> - fade on: the band of [[fade_band]], from `param`
   !> - both kernels: no `sigma` exceeds it, so every per-atom routine and the
   !>   strip run their unfaded bodies, bit for bit the unfaded kernel
   real(wp), parameter :: fade_off_start = huge(1.0_wp)
   !> Fade band width while the fade is off
   real(wp), parameter :: fade_off_width = 1.0_wp

   !> Largest `max(-a1, -a2)` the strip kernel's inline exp admits
   real(wp), parameter :: kernel_exp_max = 708.0_wp

   !> Largest box diagonal of one strip-kernel call (Bohr)
   !>
   !> Stage 1 classifies atoms and pairs against the strip's bounding box, not
   !> per point:
   !>
   !> - an atom whose reach meets the box is evaluated at every point of the
   !>   strip, a pair of two such atoms whose reach spheres meet likewise; each
   !>   point keeps only the atoms within its own reach (the gate mask, exact
   !>   zeros otherwise)
   !> - so a strip of box diagonal `D` pays for the pairs of a ball of radius
   !>   ~`r + D/2` around every point instead of `r` (reach `r` ~10-20 Bohr), a
   !>   waste ~`(1 + D/(2r))**6` that is quadratic in the atom count: <= 1.8x at
   !>   2 Bohr, 21-233x measured for strips scattered over a 45 Bohr lattice
   !> - a wider strip is cut into compact sub-strips ([[split_strip]])
   real(wp), parameter :: split_diagonal = 2.0_wp

   !> Smoothstep order of the kernel's fade `w = S_n(1 - t)` (generator constant)
   integer, parameter :: fade_order = 4

   !> Pair power `m` of the generated kernel, which bakes in `(1 - q)**4`
   !>
   !> - `param%m` reaches only the screening offset and the certificates; at any
   !>   other `m` both certificates certify nothing
   integer, parameter :: kernel_pair_power = 4

   !> Setup-time constants of the exclusion certificates
   !>
   !> - built on two bounds of the faded model: `|grad S| <= Lam(S)` at every point
   !>   of level `S` ([[lsf_exclusion_radius]]) and `PD_f <= F_n(min_A x_A)`
   !>   ([[lsf_far_exclusion_radius]]); `n = ncenters`
   !> - refreshed by `update`; `key`, `fade` and `n` record the parameters they
   !>   were built from, and a mismatch with the live `param` / `ncenters` makes
   !>   both radii return zero (safe)
   type :: cfc_certificate_type
      !> Reached certificate available (fade on, guards of [[certificate_setup]])
      logical :: valid = .false.
      !> Far certificate available (fade on, guards of [[certificate_setup]])
      logical :: far_valid = .false.
      !> Fade switch the constants belong to
      logical :: fade = .false.
      !> Parameters the constants belong to: a1, a2, c, m, a, b
      real(wp) :: key(6) = 0.0_wp
      !> Atom count the constants belong to
      integer :: n = -1
      !> Fade onset and width in sigma ([[fade_band]])
      real(wp) :: x_on = 0.0_wp, delta = 0.0_wp
      !> Fade slope maxima `mu_a`, `mu_b` ([[fade_slope_max]])
      real(wp) :: mu_a = 0.0_wp, mu_b = 0.0_wp
      !> Lipschitz bound `L` of the unfaded model (1/Bohr)
      real(wp) :: lip = 0.0_wp
      !> Level-only slope terms `K1`, `Kq`, `K0` of `Lam(S) = L + K1 + Kq e^(S/2) + K0 e^S`
      real(wp) :: k1 = 0.0_wp, kq = 0.0_wp, k0 = 0.0_wp
      !> `L' = L + K1`
      real(wp) :: lp = 0.0_wp
      !> `1/Lam(0)`, the radius per unit of `-S` inside
      real(wp) :: inv_lam0 = 0.0_wp
      !> Level where the closed-form radius peaks
      real(wp) :: s_cap = 0.0_wp
      !> Radius at `s_cap` (Bohr)
      real(wp) :: r_cap = 0.0_wp
      !> Root `xi*` of `F_n(xi) = 1`, reduced distance of the nearest atom of any surface point
      real(wp) :: xi_star = 0.0_wp
      !> `max(R_min (x_on + delta - xi*), r_cap)`, clipped at zero (Bohr)
      real(wp) :: far_radius = 0.0_wp
   end type cfc_certificate_type

   !> Smallest spatial-gradient norm [[lsf_normalized_f01_rA]] treats as
   !> non-degenerate
   !>
   !> - below it the normalized level set has no defined value and the routine
   !>   reports zero, matching `svdw_normalized_eval`
   real(wp), parameter :: normalized_eps = 1.0e-14_wp

   !> COSMO Fine Cavity LSF
   !>
   !> Concrete level set function: takes its four shape parameters through
   !> [[new]], caches the screened per-atom geometry plus the `pd*` / `qn*`
   !> accumulators via [[prepare]], and lifts every derivative tensor on demand
   !> from the generated kernel
   !>
   !> - inherits the common atom-LSF state (ncenters, mol, radii, screening
   !>   caches) from [[moist_cavity_drop_lsf_type]]
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_cfc_type
      !> CFC shape parameters (a1, a2, c, m)
      type(moist_cavity_drop_lsf_cfc_param_type) :: param

      !> Highest total derivative order `prepare` provisions (0..4)
      integer :: max_deriv = 2

      !> Reciprocal Lipschitz bound of the unfaded `S`
      real(wp) :: excl_inv_lip = 0.0_wp
      !> Exclusion-certificate constants, refreshed by `update`
      type(cfc_certificate_type) :: cert

      !> Fade switch the screening reach was built with by `update`
      logical :: bound_fade_enabled = .false.
      !> Fade band `[fade_start, fade_width, a2]` the reach was built with (fade on)
      real(wp) :: bound_fade(3) = 0.0_wp

      !* --------------------------- Per-point active list ---------------------------- *!

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
      !> Fade onset of the prepared point in sigma ([[fade_band]]); every per-atom
      !> family reads it, so all of them see the band the `pd*` / `qn*` were built with
      real(wp) :: act_fade_start = fade_off_start
      !> Fade band width of the prepared point in sigma
      real(wp) :: act_fade_width = fade_off_width

      !* ---------------------- Aggregate pseudo-density tensors ---------------------- *!
      !> `d^k PD / dr^k` summed over every active atom and pair; fixed-size
      !> components, refreshed up to `max_deriv` by every `prepare`

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

      !* --------------------------------- Lifted jet --------------------------------- *!
      !> `cfc_spatial_eval` of the `pd*` at `max_deriv`, once per `prepare`; the
      !> spatial accessors read it, orders above `max_deriv` are exact zeros

      !> Level-set value
      real(wp) :: jet0 = 0.0_wp
      !> Level-set gradient
      real(wp) :: jet1_r(ndim) = 0.0_wp
      !> Level-set Hessian
      real(wp) :: jet2_rr(ndim, ndim) = 0.0_wp
      !> Level-set third derivative
      real(wp) :: jet3_rrr(ndim, ndim, ndim) = 0.0_wp
      !> Level-set fourth derivative
      real(wp) :: jet4_rrrr(ndim, ndim, ndim, ndim) = 0.0_wp

      !* ------------------ Per-atom one-nuclear-index tensors (qn) ------------------- *!
      !> `d/dR_A of d^k PD / dr^k`, summed over atom A's own term and every pair
      !> containing A
      !>
      !> - trailing index is the active slot, the one before it the retained
      !>   nuclear component
      !> - filled to spatial order `max_deriv - 1`, i.e. to the same *total* order
      !>   as the `pd*` family

      !> Order-0 nuclear tensor [ndim, ncenters]
      real(wp), allocatable :: qn0(:, :)
      !> Order-1 nuclear tensor [ndim, ndim, ncenters]
      real(wp), allocatable :: qn1_r(:, :, :)
      !> Order-2 nuclear tensor [ndim, ndim, ndim, ncenters]
      real(wp), allocatable :: qn2_rr(:, :, :, :)
      !> Order-3 nuclear tensor [ndim, ndim, ndim, ndim, ncenters]
      real(wp), allocatable :: qn3_rrr(:, :, :, :, :)

      !* ---------------------------- Batch candidate list ---------------------------- *!

      !> Shared atom list of the current batch sub-strip, ascending candidate ids [ncenters]
      integer, allocatable :: bat_cand(:)
      !> Centers of the shared list [ndim, ncenters]
      real(wp), allocatable :: bat_center(:, :)
      !> Radii of the shared list [ncenters]
      real(wp), allocatable :: bat_radius(:)
      !> Squared screening reach of the shared list [ncenters]
      real(wp), allocatable :: bat_bsq(:)
      !> Stage-1 geometry cache and pair workspace [cfc_strip_geo_len, ncenters]
      !>
      !> - the kernel's `geo(cfc_strip_width(L), cfc_strip_geo_slots(L), nat)`
      !>   holds 160/256/448/384/312 doubles per atom at L = 0..4, so
      !>   `cfc_strip_geo_len = 448` per atom covers every order and list length
      real(wp), allocatable :: strip_geo(:, :)
      !> Stage-1 survivor list workspace of the shared list [ncenters]
      integer, allocatable :: bat_sv(:)
      !> Stage-1 per-atom reach workspace of the shared list [ncenters]
      real(wp), allocatable :: bat_reach(:)

      !* ------------------------- Batched vjp strip scratch -------------------------- *!
      !> Grown on demand by [[lsf_vjp_f1_rA_batch]]; the kernel sees them as
      !> `(vjp_strip, ...)` arrays by sequence association

      !> Rank-3 jet workspace of `cfc_strip_vjp` [10*vjp_strip*ncenters]
      real(wp), allocatable :: vjp_jet3(:)
      !> Accumulated-gradient workspace of `cfc_strip_vjp` [3*vjp_strip*ncenters]
      real(wp), allocatable :: vjp_gacc(:)
      !> Strip result when not written into the tile directly [3*vjp_strip*ncenters]
      real(wp), allocatable :: vjp_res(:)
      !> Position in `candidate_indices` by candidate id [ncenters]; defined only at
      !> the ids of the current call's list
      integer, allocatable :: vjp_pos(:)
      !> Position in `candidate_indices` of shared-list atom `a` [ncenters]
      integer, allocatable :: vjp_apos(:)
   contains
      !> Constructor: configure shape parameters and declare the candidate space
      procedure, public :: new => lsf_new
      !> Bind molecular geometry and resize the per-atom caches
      procedure, public :: update => lsf_update
      !> Point preparation: screening plus `pd*` / `qn*` accumulation
      procedure, public :: prepare => lsf_prepare
      !> Batched spatial jet through the strip kernel
      procedure, public :: evaluate_batch => lsf_evaluate_batch
      !> Batched jet-contracted nuclear gradient through the strip kernel
      procedure, public :: vjp_f1_rA_batch => lsf_vjp_f1_rA_batch
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
      !> CFC splits into atomic shares
      procedure, public :: has_atom_partition => lsf_has_atom_partition
      !> Add the atomic shares of PD at the prepared point
      procedure, public :: atom_partition => lsf_atom_partition
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
      !> Surface-free radius from the value: `|S|/L`, fade on `|tau(S)|`
      procedure, public :: exclusion_radius => lsf_exclusion_radius
      !> Surface-free radius around every unreached point (fade on)
      procedure, public :: far_exclusion_radius => lsf_far_exclusion_radius
      !> Finalizer
      final :: finalize_lsf_cfc
   end type moist_cavity_drop_lsf_cfc_type

   public :: moist_cavity_drop_lsf_cfc_type

contains

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[inout] self Object to initialize
   !> @param[in] param Configuration copied by value
   subroutine lsf_new(self, param)
      !> LSF to initialize
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Configuration; omitted means compiled defaults
      type(moist_cavity_drop_lsf_cfc_param_type), intent(in), optional :: param

      ! The atom / pair sweeps index the base's candidate-space geometry mirror
      ! directly, so candidate ids must arrive spatially sorted
      self%candidate_space = lsf_candidate_space_sorted
      self%name = "CFC"

      call self%param%init_defaults()
      if (present(param)) self%param = param
   end subroutine lsf_new

   !> Bind molecular geometry and resize the per-atom caches
   !>
   !> The base handles the spatial sort, the candidate-space geometry mirror and
   !> the screening bounds; this override only sizes the per-point buffers, once,
   !> to the molecule's atom count
   !>
   !> The `qn*` buffers are zeroed here in full; per point only the provisioned
   !> orders are refreshed, so this is what keeps the unprovisioned ones exact
   !> zeros instead of undefined memory
   !>
   !> Also sizes the batch buffers (shared list, strip geometry cache), records
   !> the fade band the screening reach was built with ([[fade_current]]) and
   !> builds the exclusion certificates ([[certificate_setup]])
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
      if (allocated(self%bat_cand)) deallocate (self%bat_cand)
      if (allocated(self%bat_center)) deallocate (self%bat_center)
      if (allocated(self%bat_radius)) deallocate (self%bat_radius)
      if (allocated(self%bat_bsq)) deallocate (self%bat_bsq)
      if (allocated(self%strip_geo)) deallocate (self%strip_geo)
      if (allocated(self%bat_sv)) deallocate (self%bat_sv)
      if (allocated(self%bat_reach)) deallocate (self%bat_reach)

      allocate (self%active_cand(n_alloc))
      allocate (self%act_atom(n_alloc))
      allocate (self%act_d(ndim, n_alloc))
      allocate (self%act_radius(n_alloc))
      allocate (self%qn0(ndim, n_alloc), source=0.0_wp)
      allocate (self%qn1_r(ndim, ndim, n_alloc), source=0.0_wp)
      allocate (self%qn2_rr(ndim, ndim, ndim, n_alloc), source=0.0_wp)
      allocate (self%qn3_rrr(ndim, ndim, ndim, ndim, n_alloc), source=0.0_wp)
      allocate (self%bat_cand(n_alloc))
      allocate (self%bat_center(ndim, n_alloc))
      allocate (self%bat_radius(n_alloc))
      allocate (self%bat_bsq(n_alloc))
      allocate (self%strip_geo(cfc_strip_geo_len, n_alloc))
      allocate (self%bat_sv(n_alloc))
      allocate (self%bat_reach(n_alloc))

      self%n_active = 0
      self%prepared_deriv = -1
      self%prepared_unreached = .false.
      ! The reach `lsf_base_update` just built belongs to this band
      self%bound_fade_enabled = self%param%fade_enabled
      self%bound_fade = [self%param%fade_start, self%param%fade_width, self%param%a2]

      call cfc_set_exclusion_bound(self)
      call certificate_setup(self)
   end subroutine lsf_update

   !> Configure the highest total derivative order `prepare` provisions
   !>
   !> Nothing is allocated here: the aggregate tensors are fixed-size components
   !> and the per-atom buffers are sized by [[update]] to the full order set
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
   !> - unreached ([[lsf_cfc_screen]]): `status = 1`, or `error` without
   !>   `status`; the cache stays readable-but-unreached, accessors abort
   !> - `error stop` before `update`
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    point  Evaluation point (3,)
   !> @param[out]   error  Invalid fade band, or parameters changed since `update`
   !>                      ([[check_point]]); or an unreached point without `status`
   !> @param[out]   status 0 ok, 1 unreached (optional)
   subroutine lsf_prepare(self, point, error, status)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Invalid fade band, or parameters changed since `update`
      type(error_type), allocatable, intent(out) :: error
      !> 0 ok, 1 unreached
      integer, intent(out), optional :: status

      if (present(status)) status = 0
      if (.not. allocated(self%cand_screen)) error stop "moist DROP LSF: prepare before update"
      call check_point(self, error)
      if (allocated(error)) return
      call lsf_cfc_screen(self, point, self%full_scan_cand)
      call self%report_unreached(error, status)
   end subroutine lsf_prepare

   !> Screen a caller-provided candidate list and refresh the caches
   !>
   !> - unreached and before `update` as in [[lsf_prepare]]
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point (3,)
   !> @param[in]    candidate_indices Atom ids in the base's sorted candidate space
   !> @param[out]   error             Invalid fade band, or parameters changed since `update`;
   !>                                 or an unreached point without `status`
   !> @param[out]   status            0 ok, 1 unreached (optional)
   subroutine lsf_prepare_subset(self, point, candidate_indices, error, status)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)
      !> Invalid fade band, or parameters changed since `update`
      type(error_type), allocatable, intent(out) :: error
      !> 0 ok, 1 unreached
      integer, intent(out), optional :: status

      if (present(status)) status = 0
      if (.not. allocated(self%cand_screen)) error stop "moist DROP LSF: prepare before update"
      call check_point(self, error)
      if (allocated(error)) return
      call lsf_cfc_screen(self, point, candidate_indices)
      call self%report_unreached(error, status)
   end subroutine lsf_prepare_subset

   !> Refuse the per-point path while it cannot evaluate the configured model
   !>
   !> - the parameter validation, as [[check_batch]]: `m = 4` always, and with
   !>   the fade on a finite, positive band with `a2 < 0`, so the kernel's own
   !>   `error stop` is never reached
   !> - fade switched, band changed, or (fade off) shape parameters changed since
   !>   `update` ([[fade_current]]): the gate reach still belongs to the old ones
   !> - on refusal the per-point cache holds no atom (`n_active = 0`)
   !>
   !> @param[inout] self  LSF instance
   !> @param[out]   error Refusal
   subroutine check_point(self, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Refusal
      type(error_type), allocatable, intent(out) :: error

      call self%param%validate(error)
      if (.not. allocated(error)) then
         if (.not. fade_current(self)) then
            call fatal_error(error, "CFC parameters changed after update; call update")
         end if
      end if
      if (allocated(error)) then
         self%n_active = 0
         self%prepared_deriv = -1
         self%prepared_unreached = .false.
      end if
   end subroutine check_point

   !> Whether the fade switch and band match those the screening reach was built with
   !>
   !> - fade switch unchanged since `update`
   !> - fade on: `[fade_start, fade_width, a2]` unchanged; the band is not read
   !>   with the fade off
   !> - both modes: the cached reach equals the live [[lsf_screening_offset]]
   !>   ([[lsf_base_reach_current]]); fade on the band end, fade off the
   !>   threshold reach of `a1, a2, c, m`, so a parameter change after `update`
   !>   or a `set_screening_threshold` refresh under other parameters is caught
   !>
   !> @param[in] self    LSF instance
   !> @returns   current Reach built for the live fade parameters
   pure function fade_current(self) result(current)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Reach built for the live fade parameters
      logical :: current

      current = self%param%fade_enabled .eqv. self%bound_fade_enabled
      if (current .and. self%param%fade_enabled) then
         current = all([self%param%fade_start, self%param%fade_width, self%param%a2] == &
                       self%bound_fade)
      end if
      if (current) current = self%reach_current()
   end function fade_current

   !> Fade band of both kernels, in `sigma_A = ||r - R_A||/R_a - 1`
   !>
   !> - stage 1 of the strip reads it per call, the per-atom routines through
   !>   `act_fade_start`, `act_fade_width` cached at `prepare`
   !>
   !> - fade off: [[fade_off_start]], [[fade_off_width]], the untruncated model
   !> - fade on: `fade_start/|a2|`, `fade_width/|a2|`, e-folds of `u_A` turned
   !>   into reduced distance (each atom's own radius as unit)
   !>
   !> @param[in]  self  LSF instance
   !> @param[out] start Fade onset in sigma
   !> @param[out] width Fade band width in sigma
   pure subroutine fade_band(self, start, width)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Fade onset in sigma
      real(wp), intent(out) :: start
      !> Fade band width in sigma
      real(wp), intent(out) :: width

      if (.not. self%param%fade_enabled) then
         start = fade_off_start
         width = fade_off_width
         return
      end if
      start = self%param%fade_start/abs(self%param%a2)
      width = self%param%fade_width/abs(self%param%a2)
   end subroutine fade_band

   !> Run the base screening gate and rebuild the per-point accumulators
   !>
   !> Pass one is the base's allocation-free reject test; pass two touches only
   !> the survivors and fills the per-atom geometry; pass three is the atom sweep
   !> and pass four the pair sweep, each accumulating the spatial family and --
   !> when a nuclear order was provisioned -- the one-nuclear-index family in the
   !> same iteration
   !>
   !> The pair sweep is the O(n_active**2) part and the single most expensive
   !> thing the CFC level set does
   !>
   !> Pass five lifts the spatial jet at `max_deriv` into the cache (`jet*`)
   !>
   !> Unreached (`prepared_unreached`): no survivor, `PD` not positive (the lift
   !> divides by it), or any lifted output through `max_deriv` non-finite (a
   !> tiny positive `PD` far out with screening off, a point on a nucleus from
   !> order 1), the per-point twin of the batch's [[poison_unreached]]
   !>
   !> - cannot be made linear (see the module header)
   !> - what it can be, and is, is a single kernel call per pair yielding both
   !>   families already assembled, with no split-to-spatial recombination left
   !>   to this module
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
      !> Fade band in sigma
      real(wp) :: fs, fw

      self%n_active = 0
      self%prepared_deriv = self%max_deriv
      self%prepared_unreached = .true.
      if (.not. allocated(self%cand_screen)) error stop "moist DROP LSF: prepare before update"
      call fade_band(self, fs, fw)
      self%act_fade_start = fs
      self%act_fade_width = fw

      call lsf_screen_records(self%cand_screen, point, candidate_indices, self%active_cand, n)
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
                                   self%param%a1, &
                                   self%act_fade_start, self%act_fade_width, md, &
                                   self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                                   self%pd4_rrrr)
         if (nd < 0) cycle
         call cfc_atomic_nuclear_eval(self%act_d(:, i), self%act_radius(i), &
                                      self%param%a1, &
                                      self%act_fade_start, self%act_fade_width, nd, &
                                      self%qn0(:, i), self%qn1_r(:, :, i), &
                                      self%qn2_rr(:, :, :, i), self%qn3_rrr(:, :, :, :, i))
      end do

      do i = 1, n
         do j = i + 1, n
            call cfc_pair_spatial_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%param%a2, self%param%c, &
                                       self%act_fade_start, self%act_fade_width, md, &
                                       self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                                       self%pd4_rrrr)
            if (nd < 0) cycle
            call cfc_pair_nuclear_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%param%a2, self%param%c, &
                                       self%act_fade_start, self%act_fade_width, nd, &
                                       self%qn0(:, i), self%qn0(:, j), &
                                       self%qn1_r(:, :, i), self%qn1_r(:, :, j), &
                                       self%qn2_rr(:, :, :, i), self%qn2_rr(:, :, :, j), &
                                       self%qn3_rrr(:, :, :, :, i), self%qn3_rrr(:, :, :, :, j))
         end do
      end do
      if (.not. (self%pd0 > 0.0_wp)) return
      call cfc_spatial_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, self%pd4_rrrr, md, &
                            self%jet0, self%jet1_r, self%jet2_rr, self%jet3_rrr, self%jet4_rrrr)
      self%prepared_unreached = .not. jet_finite(md, self%jet0, self%jet1_r, self%jet2_rr, &
                                                 self%jet3_rrr, self%jet4_rrrr)
   end subroutine lsf_cfc_screen

   !> Whether every output of the dense jet through `order` is finite
   !>
   !> @param[in] order    Highest spatial order (0..4)
   !> @param[in] f0       Value
   !> @param[in] f1_r     Gradient
   !> @param[in] f2_rr    Hessian
   !> @param[in] f3_rrr   Third derivative
   !> @param[in] f4_rrrr  Fourth derivative
   !> @returns   finite   Every output through `order` finite
   pure function jet_finite(order, f0, f1_r, f2_rr, f3_rrr, f4_rrrr) result(finite)
      !> Highest spatial order
      integer, intent(in) :: order
      !> Value
      real(wp), intent(in) :: f0
      !> Gradient
      real(wp), intent(in) :: f1_r(ndim)
      !> Hessian
      real(wp), intent(in) :: f2_rr(ndim, ndim)
      !> Third derivative
      real(wp), intent(in) :: f3_rrr(ndim, ndim, ndim)
      !> Fourth derivative
      real(wp), intent(in) :: f4_rrrr(ndim, ndim, ndim, ndim)
      !> Every output through `order` finite
      logical :: finite

      finite = ieee_is_finite(f0)
      if (order >= 1) finite = finite .and. all(ieee_is_finite(f1_r))
      if (order >= 2) finite = finite .and. all(ieee_is_finite(f2_rr))
      if (order >= 3) finite = finite .and. all(ieee_is_finite(f3_rrr))
      if (order >= 4) finite = finite .and. all(ieee_is_finite(f4_rrrr))
   end function jet_finite

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

      call self%require_deriv(0, "f0")
      val = self%jet0
   end subroutine lsf_f0

   !> Level-set value, spatial gradient and spatial Hessian
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    LSF value (optional)
   !> @param[out] lsf1_r  Spatial gradient (optional)
   !> @param[out] lsf2_rr Spatial Hessian, packed [6] (optional)
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian, packed
      real(wp), intent(out), optional :: lsf2_rr(:)

      !> Highest spatial-derivative order requested
      integer :: md

      md = 0
      if (present(lsf1_r)) md = 1
      if (present(lsf2_rr)) md = 2
      call self%require_deriv(md, "f012_r")

      if (present(lsf0)) lsf0 = self%jet0
      if (present(lsf1_r)) lsf1_r = self%jet1_r
      if (present(lsf2_rr)) call pack_sym_matrix(self%jet2_rr, lsf2_rr, ndim)
   end subroutine lsf_f012_r

   !> Third spatial derivative, plus optionally the lower orders
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     LSF value (optional)
   !> @param[out] lsf1_r   Spatial gradient (optional)
   !> @param[out] lsf2_rr  Spatial Hessian, packed [6] (optional)
   !> @param[out] lsf3_rrr Third spatial derivative, packed [10]
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Spatial Hessian, packed
      real(wp), intent(out), optional :: lsf2_rr(:)
      !> Third spatial derivative, packed
      real(wp), intent(out) :: lsf3_rrr(:)

      call self%require_deriv(3, "f3_rrr")

      if (present(lsf0)) lsf0 = self%jet0
      if (present(lsf1_r)) lsf1_r = self%jet1_r
      if (present(lsf2_rr)) call pack_sym_matrix(self%jet2_rr, lsf2_rr, ndim)
      call pack_sym_tensor3(self%jet3_rrr, lsf3_rrr, ndim)
   end subroutine lsf_f3_rrr

   !> Fourth spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Fourth spatial derivative, packed [15]
   subroutine lsf_f4_rrrr(self, lsf4_rrrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Fourth spatial derivative, packed
      real(wp), intent(out) :: lsf4_rrrr(:)

      call self%require_deriv(4, "f4_rrrr")

      call pack_sym_tensor4(self%jet4_rrrr, lsf4_rrrr, ndim)
   end subroutine lsf_f4_rrrr

   !* ================================================================================= *!
   !*                                 Batched evaluation                                *!
   !* ================================================================================= *!

   !> Level-set jet at a batch of points through the strip kernel
   !>
   !> Override of [[lsf_base_evaluate_batch]], same contract:
   !>
   !> - points cut into strips of `cfc_strip_width(order)`, each strip into
   !>   compact sub-strips of box diagonal <= [[split_diagonal]]
   !>   ([[split_strip]]); per sub-strip stage 1 over one shared atom list, then
   !>   the lift; `order` alone governs (no `require_deriv`, `max_deriv` untouched)
   !> - `candidate_indices` present: the shared list of every sub-strip, sorted
   !>   ascending in candidate space before use; absent: per sub-strip
   !>   `candidates_in_box` over its bounding box, a superset of every point's
   !>   active set, the kernel's per-point gate mask making the result the full
   !>   screen per point
   !> - unreached column: `status = 1`, or `error` naming the column when
   !>   `status` is absent; every output of the column through `order` holds a
   !>   quiet NaN, so it cannot pass for a surface point; status 1 means
   !>   - `PD <= 0` after stage 1, every atom of the list masked by the gate (no
   !>     atom within reach); its sums are patched to `PD = 1` so the lift
   !>     stays finite, its derivative sums are exact zeros already
   !>   - or any non-finite output through `order` after the lift (the lift
   !>     divides by `PD`), e.g. a point on a nucleus, or a tiny `PD` whose
   !>     reciprocal powers overflow at the higher orders
   !>   - with screening off the inline exp clamp keeps `PD >= exp(-708)` per
   !>     atom, so a far point has a large finite value, its higher orders
   !>     possibly not
   !> - `status = 2` only on a refusal of the whole batch, never per column
   !> - on-nucleus sub-strip (stage 1's `n_zero > 0`): lifted as is, the strip
   !>   kernel keeping the cusp in the nucleus lane; that column keeps its finite
   !>   value at order 0 and comes back unreached from order 1, every other
   !>   column of the sub-strip is unaffected
   !> - the per-point cache is invalidated (`prepared_deriv = -1`)
   !> - refused with `error`, no column written, `status` 2 ([[check_batch]]):
   !>   `m /= 4`, invalid fade band, parameters changed since `update`, exponents
   !>   outside the kernel's range
   !> - the single-point cache is unspecified after the call
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    points            Evaluation points (3, n)
   !> @param[in]    order             Highest spatial derivative to write (0..4)
   !> @param[out]   lsf0              Values (n)
   !> @param[out]   lsf1_r            Gradients (3, n), written for order >= 1
   !> @param[out]   lsf2_rr_packed    Packed Hessians (6, n), the lift's own slots,
   !>                                 written for order >= 2
   !> @param[out]   lsf3_rrr          Packed third derivatives (10, n), written for order >= 3
   !> @param[out]   lsf4_rrrr         Packed fourth derivatives (15, n), written for order 4
   !> @param[out]   status            Per column 0 ok, 1 unreached, 2 refused (optional)
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   !> @param[out]   error             Refused configuration, or first unreached column
   !>                                 when `status` is absent
   subroutine lsf_evaluate_batch(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, lsf3_rrr, &
                                 lsf4_rrrr, status, candidate_indices, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Highest spatial derivative to write
      integer, intent(in) :: order
      !> Values
      real(wp), intent(out), contiguous :: lsf0(:)
      !> Gradients
      real(wp), intent(out), optional, contiguous :: lsf1_r(:, :)
      !> Packed Hessians
      real(wp), intent(out), optional, contiguous :: lsf2_rr_packed(:, :)
      !> Packed third derivatives
      real(wp), intent(out), optional, contiguous :: lsf3_rrr(:, :)
      !> Packed fourth derivatives
      real(wp), intent(out), optional, contiguous :: lsf4_rrrr(:, :)
      !> Per-column status
      integer, intent(out), optional, contiguous :: status(:)
      !> Shared atom list
      integer, intent(in), optional :: candidate_indices(:)
      !> Refused configuration or first unreached column
      type(error_type), allocatable, intent(out) :: error

      !> Lift outputs of one strip in sub-strip order, packed
      real(wp) :: t0(max_strip), t1(ndim, max_strip), t2(np2, max_strip)
      real(wp) :: t3(np3, max_strip), t4(np4, max_strip)
      !> Column reached at least one atom, sub-strip order
      logical :: reached(max_strip)
      !> Points of one sub-strip, gathered
      real(wp) :: pbuf(ndim, max_strip)
      !> Strip column of each sub-strip slot, sub-strip starts and lengths
      integer :: perm(max_strip), seg0(max_strip), segn(max_strip)
      !> Sub-strip bounding box
      real(wp) :: lo(ndim), hi(ndim)
      !> Point count, strip width, strip bounds and size, list size, zero-distance pairs
      integer :: n, width, s0, s1, np, nc, n_zero
      !> Candidate id and displacement of the nucleus check on the error path
      integer :: c
      real(wp) :: dc(ndim)
      !> Sub-strip count, index, first slot and size
      integer :: nseg, k, q0, q1, m
      !> Slot and column indices
      integer :: j, i
      !> Column label and list length for the error message
      character(len=16) :: col, ncol

      n = size(points, 2)
      call lsf_batch_require_outputs("CFC", points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                     lsf3_rrr, lsf4_rrrr, status)
      if (n == 0) return
      if (.not. allocated(self%cand_screen)) then
         error stop "moist DROP LSF CFC: evaluate_batch before update"
      end if
      call check_batch(self, error)
      if (allocated(error)) then
         if (present(status)) status(1:n) = 2
         return
      end if
      ! The batch leaves the per-point cache alone but stale: accessors abort until `prepare`
      self%prepared_deriv = -1

      width = cfc_strip_width(order)
      nc = 0
      if (present(candidate_indices)) then
         call lsf_sorted_candidates(candidate_indices, self%bat_cand, nc, self%ncenters)
         call gather_list(self, nc)
      end if

      do s0 = 1, n, width
         s1 = min(n, s0 + width - 1)
         np = s1 - s0 + 1
         call split_strip(np, points(:, s0:s1), perm, seg0, segn, nseg)

         do k = 1, nseg
            q0 = seg0(k)
            m = segn(k)
            q1 = q0 + m - 1
            do j = 1, m
               pbuf(:, j) = points(:, s0 - 1 + perm(q0 + j - 1))
            end do
            if (.not. present(candidate_indices)) then
               lo = pbuf(:, 1)
               hi = pbuf(:, 1)
               do j = 2, m
                  lo = min(lo, pbuf(:, j))
                  hi = max(hi, pbuf(:, j))
               end do
               ! The cell list itself, a subobject disjoint from `bat_cand`
               call self%cand_cells%query_box(lo, hi, self%bat_cand, nc)
               call gather_list(self, nc)
            end if
            call strip_eval(self, order, m, pbuf, nc, t0(q0:q1), t1(:, q0:q1), t2(:, q0:q1), &
                            t3(:, q0:q1), t4(:, q0:q1), reached(q0:q1), n_zero)
         end do
         call poison_unreached(order, np, t0, t1, t2, t3, t4, reached)

         do j = 1, np
            i = s0 - 1 + perm(j)
            lsf0(i) = t0(j)
            if (order >= 1) lsf1_r(:, i) = t1(:, j)
            if (order >= 2) lsf2_rr_packed(:, i) = t2(:, j)
            if (order >= 3) lsf3_rrr(:, i) = t3(:, j)
            if (order >= 4) lsf4_rrrr(:, i) = t4(:, j)
         end do

         if (present(status)) then
            do j = 1, np
               status(s0 - 1 + perm(j)) = merge(0, 1, reached(j))
            end do
         else if (.not. all(reached(1:np))) then
            i = n + 1
            do j = 1, np
               if (.not. reached(j)) i = min(i, s0 - 1 + perm(j))
            end do
            write (col, "(i0)") i
            write (ncol, "(i0)") nc
            call fatal_error(error, "evaluate_batch: point "//trim(col)// &
                             " is unreached, no atom within the screening reach or PD "// &
                             "underflowed (CFC, "//trim(ncol)//" candidates "// &
                             merge("given   ", "from box", present(candidate_indices))//")")
            ! Name the cusp when the column sits on a nucleus; error path only, O(ncenters)
            if (order >= 1) then
               do c = 1, self%ncenters
                  dc = points(:, i) - self%cand_screen(1:ndim, c)
                  if (dc(1)*dc(1) + dc(2)*dc(2) + dc(3)*dc(3) /= 0.0_wp) cycle
                  error%message = error%message//" (on-nucleus column at order >= 1)"
                  exit
               end do
            end if
            return
         end if
      end do
   end subroutine lsf_evaluate_batch

   !> Refuse a batch the strip kernel cannot evaluate as configured
   !>
   !> - the parameter validation: `m = 4` always, and with the fade on a
   !>   finite, positive band with `a2 < 0`
   !> - fade switched, band changed, or (fade off) shape parameters changed since
   !>   `update` ([[fade_current]])
   !> - exponents outside the kernel's range: `a1 <= 0`, `a2 <= 0`,
   !>   `max(-a1, -a2) <= 708` (inline exp), checked here so the kernel's own
   !>   `error stop` is never reached
   !>
   !> @param[inout] self  LSF instance
   !> @param[out]   error Refusal
   subroutine check_batch(self, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Refusal
      type(error_type), allocatable, intent(out) :: error

      call self%param%validate(error)
      if (allocated(error)) return
      if (.not. fade_current(self)) then
         call fatal_error(error, "CFC parameters changed after update; call update")
         return
      end if
      associate (a1 => self%param%a1, a2 => self%param%a2)
         if (.not. (a1 <= 0.0_wp .and. a2 <= 0.0_wp .and. max(-a1, -a2) <= kernel_exp_max)) then
            call fatal_error(error, "CFC evaluate_batch needs a1 <= 0, a2 <= 0 and "// &
                             "max(-a1, -a2) <= 708 (strip kernel)")
         end if
      end associate
   end subroutine check_batch

   !> Gather centers, radii and gate reach of the shared list `bat_cand(1:nc)`
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    nc   Length of the shared list
   pure subroutine gather_list(self, nc)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Length of the shared list
      integer, intent(in) :: nc

      !> List index and candidate id
      integer :: i, c

      do i = 1, nc
         c = self%bat_cand(i)
         self%bat_center(1, i) = self%cand_screen(1, c)
         self%bat_center(2, i) = self%cand_screen(2, c)
         self%bat_center(3, i) = self%cand_screen(3, c)
         self%bat_radius(i) = self%cand_radii(c)
         self%bat_bsq(i) = self%cand_screen(lsf_cand_bound_sq, c)
      end do
   end subroutine gather_list

   !> Cut a strip into compact sub-strips by recursive median bisection
   !>
   !> - a segment whose box diagonal exceeds [[split_diagonal]] is sorted along
   !>   its longest box axis (stable insertion sort) and cut at the median; a
   !>   segment within the cap, or of one point, is a sub-strip
   !> - sub-strip k holds the strip columns `perm(seg0(k) : seg0(k) + segn(k) - 1)`
   !> - a compact strip is one sub-strip in column order, so the kernel call is
   !>   the unsplit one
   !> - the result depends on the points only, never on the candidate list
   !> - allocation-free: at most `np <= max_strip` pending segments, all disjoint
   !>
   !> @param[in]  np     Points in the strip
   !> @param[in]  points Strip points [ndim, np]
   !> @param[out] perm   Strip column of each sub-strip slot [np]
   !> @param[out] seg0   First slot of each sub-strip [np]
   !> @param[out] segn   Points of each sub-strip [np]
   !> @param[out] nseg   Number of sub-strips
   pure subroutine split_strip(np, points, perm, seg0, segn, nseg)
      !> Points in the strip
      integer, intent(in) :: np
      !> Strip points
      real(wp), intent(in) :: points(ndim, np)
      !> Strip column of each sub-strip slot
      integer, intent(out) :: perm(np)
      !> First slot of each sub-strip
      integer, intent(out) :: seg0(np)
      !> Points of each sub-strip
      integer, intent(out) :: segn(np)
      !> Number of sub-strips
      integer, intent(out) :: nseg

      !> Pending segments as slot ranges, fixed size (`np <= max_strip`, no automatic arrays)
      integer :: stack_lo(max_strip), stack_hi(max_strip)
      !> Stack top, segment bounds, median, sort axis, slot indices, sorted column
      integer :: top, lo, hi, mid, axis, i, j, key
      !> Segment bounding box, its extent, the sort key
      real(wp) :: blo(ndim), bhi(ndim), ext(ndim), x

      do i = 1, np
         perm(i) = i
      end do
      nseg = 0
      top = 1
      stack_lo(1) = 1
      stack_hi(1) = np
      do while (top > 0)
         lo = stack_lo(top)
         hi = stack_hi(top)
         top = top - 1
         blo = points(:, perm(lo))
         bhi = blo
         do i = lo + 1, hi
            blo = min(blo, points(:, perm(i)))
            bhi = max(bhi, points(:, perm(i)))
         end do
         ext = bhi - blo
         if (hi == lo .or. ext(1)*ext(1) + ext(2)*ext(2) + ext(3)*ext(3) <= &
             split_diagonal*split_diagonal) then
            nseg = nseg + 1
            seg0(nseg) = lo
            segn(nseg) = hi - lo + 1
            cycle
         end if

         axis = maxloc(ext, dim=1)
         do i = lo + 1, hi
            key = perm(i)
            x = points(axis, key)
            j = i - 1
            do while (j >= lo)
               if (points(axis, perm(j)) <= x) exit
               perm(j + 1) = perm(j)
               j = j - 1
            end do
            perm(j + 1) = key
         end do

         ! Right half below the left on the stack, so the left half is cut first
         mid = (lo + hi)/2
         top = top + 1
         stack_lo(top) = mid + 1
         stack_hi(top) = hi
         top = top + 1
         stack_lo(top) = lo
         stack_hi(top) = mid
      end do
   end subroutine split_strip

   !> Both stages of the strip kernel on one sub-strip over the gathered list
   !>
   !> - `n_zero > 0` (a point on a nucleus): lifted all the same; stage 1 keeps
   !>   the cusp in the nucleus lane, whose non-finite derivative sums give
   !>   non-finite outputs for [[poison_unreached]] (its value is finite), every
   !>   other lane untouched
   !> - columns with `PD <= 0` marked unreached and patched to `PD = 1` before
   !>   the lift
   !>
   !> Stage-1 sums live in fixed buffers of the widest strip; the explicit-shape
   !> kernel dummies see them as `(np, ...)` by sequence association, so nothing
   !> is allocated per strip
   !>
   !> @param[inout] self    LSF instance with `bat_*` gathered (`strip_geo` written)
   !> @param[in]    order   Highest spatial order (0..4)
   !> @param[in]    np      Points in the sub-strip
   !> @param[in]    points  Sub-strip points [ndim, np]
   !> @param[in]    nc      Length of the shared list
   !> @param[out]   l0      Values [np]
   !> @param[out]   l1      Gradients [ndim, np]
   !> @param[out]   l2      Packed Hessians [6, np]
   !> @param[out]   l3      Packed third derivatives [10, np]
   !> @param[out]   l4      Packed fourth derivatives [15, np]
   !> @param[out]   reached Column has `PD > 0` [np]
   !> @param[out]   n_zero  Point-atom pairs at zero distance
   subroutine strip_eval(self, order, np, points, nc, l0, l1, l2, l3, l4, reached, n_zero)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Highest spatial order
      integer, intent(in) :: order
      !> Points in the sub-strip
      integer, intent(in) :: np
      !> Sub-strip points
      real(wp), intent(in) :: points(ndim, np)
      !> Length of the shared list
      integer, intent(in) :: nc
      !> Values
      real(wp), intent(out) :: l0(np)
      !> Gradients
      real(wp), intent(out) :: l1(ndim, np)
      !> Packed Hessians
      real(wp), intent(out) :: l2(np2, np)
      !> Packed third derivatives
      real(wp), intent(out) :: l3(np3, np)
      !> Packed fourth derivatives
      real(wp), intent(out) :: l4(np4, np)
      !> Column has `PD > 0`
      logical, intent(out) :: reached(np)
      !> Point-atom pairs at zero distance
      integer, intent(out) :: n_zero

      !> Stage-1 sums of the widest strip, seen as `(np, npack)` by the kernel
      real(wp) :: b0(max_strip), b1(max_strip*ndim), b2(max_strip*np2)
      real(wp) :: b3(max_strip*np3), b4(max_strip*np4)
      !> Fade band in sigma
      real(wp) :: fs, fw
      !> Point index
      integer :: p

      call fade_band(self, fs, fw)
      ! `strip_geo` is the kernel's rank-3 explicit-shape `geo` by sequence
      ! association; `bat_sv` and `bat_reach` are its `(nat)` workspaces
      associate (c => self%bat_center, r => self%bat_radius, q => self%bat_bsq, &
                 a1 => self%param%a1, a2 => self%param%a2, cp => self%param%c)
         select case (order)
         case (0)
            call cfc_strip_powersums_o0(np, points, nc, c, r, q, a1, a2, cp, fs, fw, b0, &
                                        n_zero, self%strip_geo, self%bat_sv, &
                                        self%bat_reach)
         case (1)
            call cfc_strip_powersums_o1(np, points, nc, c, r, q, a1, a2, cp, fs, fw, b0, b1, &
                                        n_zero, self%strip_geo, self%bat_sv, &
                                        self%bat_reach)
         case (2)
            call cfc_strip_powersums_o2(np, points, nc, c, r, q, a1, a2, cp, fs, fw, b0, b1, &
                                        b2, n_zero, self%strip_geo, self%bat_sv, &
                                        self%bat_reach)
         case (3)
            call cfc_strip_powersums_o3(np, points, nc, c, r, q, a1, a2, cp, fs, fw, b0, b1, &
                                        b2, b3, n_zero, self%strip_geo, self%bat_sv, &
                                        self%bat_reach)
         case default
            call cfc_strip_powersums_o4(np, points, nc, c, r, q, a1, a2, cp, fs, fw, b0, b1, &
                                        b2, b3, b4, n_zero, self%strip_geo, self%bat_sv, &
                                        self%bat_reach)
         end select
      end associate

      do p = 1, np
         reached(p) = b0(p) > 0.0_wp
         if (.not. reached(p)) b0(p) = 1.0_wp
      end do

      select case (order)
      case (0)
         call cfc_strip_lift_o0(np, b0, l0)
      case (1)
         call cfc_strip_lift_o1(np, b0, b1, l0, l1)
      case (2)
         call cfc_strip_lift_o2(np, b0, b1, b2, l0, l1, l2)
      case (3)
         call cfc_strip_lift_o3(np, b0, b1, b2, b3, l0, l1, l2, l3)
      case default
         call cfc_strip_lift_o4(np, b0, b1, b2, b3, b4, l0, l1, l2, l3, l4)
      end select
   end subroutine strip_eval

   !> Mark non-finite columns unreached and fill every unreached column with NaN
   !>
   !> - a reached column whose outputs through `order` are not all finite
   !>   becomes unreached
   !> - every output through `order` of an unreached column becomes a quiet NaN,
   !>   so a caller ignoring `status` cannot read it as a surface point
   !> - outputs above `order` are left alone
   !>
   !> @param[in]    order   Highest spatial order (0..4)
   !> @param[in]    np      Points in the strip
   !> @param[inout] l0      Values [np]
   !> @param[inout] l1      Gradients [ndim, np]
   !> @param[inout] l2      Packed Hessians [6, np]
   !> @param[inout] l3      Packed third derivatives [10, np]
   !> @param[inout] l4      Packed fourth derivatives [15, np]
   !> @param[inout] reached Column reached at least one atom [np]
   pure subroutine poison_unreached(order, np, l0, l1, l2, l3, l4, reached)
      !> Highest spatial order
      integer, intent(in) :: order
      !> Points in the strip
      integer, intent(in) :: np
      !> Values
      real(wp), intent(inout) :: l0(np)
      !> Gradients
      real(wp), intent(inout) :: l1(ndim, np)
      !> Packed Hessians
      real(wp), intent(inout) :: l2(np2, np)
      !> Packed third derivatives
      real(wp), intent(inout) :: l3(np3, np)
      !> Packed fourth derivatives
      real(wp), intent(inout) :: l4(np4, np)
      !> Column reached at least one atom
      logical, intent(inout) :: reached(np)

      !> Quiet NaN
      real(wp) :: qnan
      !> Point index
      integer :: p
      !> Every output of the column through `order` is finite
      logical :: finite

      qnan = ieee_value(qnan, ieee_quiet_nan)
      do p = 1, np
         if (reached(p)) then
            finite = ieee_is_finite(l0(p))
            if (order >= 1) finite = finite .and. all(ieee_is_finite(l1(:, p)))
            if (order >= 2) finite = finite .and. all(ieee_is_finite(l2(:, p)))
            if (order >= 3) finite = finite .and. all(ieee_is_finite(l3(:, p)))
            if (order >= 4) finite = finite .and. all(ieee_is_finite(l4(:, p)))
            reached(p) = finite
         end if
         if (reached(p)) cycle
         l0(p) = qnan
         if (order >= 1) l1(:, p) = qnan
         if (order >= 2) l2(:, p) = qnan
         if (order >= 3) l3(:, p) = qnan
         if (order >= 4) l4(:, p) = qnan
      end do
   end subroutine poison_unreached

   !* ================================================================================= *!
   !*                     Batched jet-contracted nuclear gradient                       *!
   !* ================================================================================= *!

   !> Batched jet-contracted nuclear gradient through the strip kernel
   !>
   !> Override of [[lsf_base_vjp_f1_rA_batch]], same contract:
   !>
   !> - per strip of `cfc_strip_width(2)` points cut into compact sub-strips of box
   !>   diagonal <= [[split_diagonal]] ([[split_strip]]): stage 1 at order 2 over one
   !>   shared atom list with the geometry cache, the order-2 lift (the reached
   !>   decision), `cfc_strip_lift_vjp`, then `cfc_strip_vjp` on the cache
   !>   ([[strip_vjp_pass]])
   !> - the shared list as `evaluate_batch`: `candidate_indices` sorted ascending in
   !>   candidate space, or per sub-strip `candidates_in_box` over its box
   !> - the ids are validated ([[lsf_sorted_candidates]]) before the position map is built
   !> - list present, already ascending and the strip one sub-strip: the kernel writes
   !>   the tile columns in place; any other case scatters the sub-strip result into
   !>   the caller's positions after zeroing the column
   !> - `max_deriv` is not read on the strip path, `qn*` are never formed; the on-nucleus
   !>   fallback sets 3 for its duration and restores it
   !> - unreached column (`PD <= 0` after stage 1, a non-finite value, gradient or
   !>   Hessian of the order-2 lift, or a non-finite result column, which catches an
   !>   adjoint that overflows where the lift does not): NaN in `res(:, :, p)` and
   !>   `status = 1`, or `error` naming the column when `status` is absent; the other
   !>   columns are computed
   !> - on-nucleus column (stage 1's `n_zero > 0`): the sub-strip is discarded, the
   !>   column is unreached ([[nucleus_unreached]]: the cusp makes the per-point jet
   !>   non-finite from order 1) and the remaining columns are re-run as a sub-strip of
   !>   their own, so they match the tile without the on-nucleus columns bit for bit
   !> - refused with `error` "vjp_f1_rA_batch failed at point 1: ...", no column written,
   !>   `status` 2 ([[check_batch]]): `m /= 4`, invalid fade band, parameters changed
   !>   since `update`, exponents outside the kernel's range
   !> - the per-point cache is invalidated (`prepared_deriv = -1`)
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    points            Evaluation points (3, np)
   !> @param[in]    w0                Adjoint weights of the value (np)
   !> @param[in]    w1                Adjoint weights of the gradient (3, np)
   !> @param[in]    w2                Packed symmetric part of the Hessian seed (6, np)
   !> @param[out]   res               Contracted nuclear gradient (3, ncand, np)
   !> @param[out]   status            Per-column status (np), 0 ok, 1 unreached, 2 refused
   !>                                 (optional)
   !> @param[out]   error             Refused configuration, or first unreached column
   !>                                 when `status` is absent
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   subroutine lsf_vjp_f1_rA_batch(self, points, w0, w1, w2, res, status, error, &
                                  candidate_indices)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Value weights
      real(wp), intent(in), contiguous :: w0(:)
      !> Gradient weights
      real(wp), intent(in), contiguous :: w1(:, :)
      !> Packed symmetric Hessian seeds
      real(wp), intent(in), contiguous :: w2(:, :)
      !> Contracted nuclear gradient
      real(wp), intent(out), contiguous :: res(:, :, :)
      !> Per-column status
      integer, intent(out), optional, contiguous :: status(:)
      !> Refused configuration or first unreached column
      type(error_type), allocatable, intent(out) :: error
      !> Shared atom list
      integer, intent(in), optional :: candidate_indices(:)

      !> Sub-strip points and seeds, gathered in sub-strip order
      real(wp) :: pb(ndim, vjp_strip), b0(vjp_strip), b1(ndim, vjp_strip)
      real(wp) :: b2(np2, vjp_strip)
      !> Sub-strip points and seeds without the on-nucleus columns
      real(wp) :: cp(ndim, vjp_strip), c0(vjp_strip), c1(ndim, vjp_strip)
      real(wp) :: c2(np2, vjp_strip)
      !> Sub-strip bounding box and fill of an unreached column
      real(wp) :: lo(ndim), hi(ndim), nan
      !> Column reached at least one atom: strip slots, sub-strip, reduced sub-strip
      logical :: reached(vjp_strip), rsub(vjp_strip), rc(vjp_strip)
      !> Column holds a point on a nucleus
      logical :: nuc(vjp_strip)
      !> Write the tile columns in place (list ascending), and for this strip
      logical :: direct, in_place
      !> List absent: the box query `bat_cand` is its own tile position
      logical :: box
      !> Strip column of each sub-strip slot, sub-strip starts and lengths
      integer :: perm(vjp_strip), seg0(vjp_strip), segn(vjp_strip)
      !> Sub-strip column of each reduced column
      integer :: lane(vjp_strip)
      !> Point count, list sizes, strip bounds and size
      integer :: n, ncand, nc, s0, s1, np
      !> Sub-strip count, index, first slot and size, reduced size
      integer :: nseg, k, q0, m, nk
      !> Sub-strip column, tile column, list atom, zero-distance pair counts
      integer :: j, i, a, n_zero, n_zero_c
      !> Column label for the error message
      character(len=16) :: col
      character(len=:), allocatable :: msg

      n = size(points, 2)
      ncand = self%ncenters
      if (present(candidate_indices)) ncand = size(candidate_indices)
      if (present(status)) then
         msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand, &
                                           status)
      else
         msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand)
      end if
      if (len(msg) > 0) error stop msg
      if (n == 0) return
      if (.not. allocated(self%cand_screen)) then
         error stop "moist DROP LSF CFC: vjp_f1_rA_batch before update"
      end if
      call check_batch(self, error)
      if (allocated(error)) then
         error%message = "vjp_f1_rA_batch failed at point 1: "//error%message
         if (present(status)) status(1:n) = 2
         return
      end if
      ! The batch overwrites the per-point cache: accessors abort until `prepare`
      self%prepared_deriv = -1
      call grow_vjp_scratch(self)
      nan = ieee_value(0.0_wp, ieee_quiet_nan)

      ! Tile position of each shared-list atom: its index in `candidate_indices`, or
      ! the candidate id itself without a list (the full id range is the tile)
      box = .not. present(candidate_indices)
      direct = .false.
      nc = 0
      if (.not. box) then
         ! Validates the ids, so the position map below may index with them
         call lsf_sorted_candidates(candidate_indices, self%bat_cand, nc, self%ncenters)
         do k = 1, ncand
            self%vjp_pos(candidate_indices(k)) = k
         end do
         call gather_list(self, nc)
         do a = 1, nc
            self%vjp_apos(a) = self%vjp_pos(self%bat_cand(a))
         end do
         direct = all(self%bat_cand(1:nc) == candidate_indices)
      end if

      do s0 = 1, n, vjp_strip
         s1 = min(n, s0 + vjp_strip - 1)
         np = s1 - s0 + 1
         call split_strip(np, points(:, s0:s1), perm, seg0, segn, nseg)
         ! One sub-strip is the strip in column order: the kernel may write the tile
         in_place = direct .and. nseg == 1

         do k = 1, nseg
            q0 = seg0(k)
            m = segn(k)
            do j = 1, m
               i = s0 - 1 + perm(q0 + j - 1)
               pb(:, j) = points(:, i)
               b0(j) = w0(i)
               b1(:, j) = w1(:, i)
               b2(:, j) = w2(:, i)
            end do
            if (.not. present(candidate_indices)) then
               lo = pb(:, 1)
               hi = pb(:, 1)
               do j = 2, m
                  lo = min(lo, pb(:, j))
                  hi = max(hi, pb(:, j))
               end do
               ! The cell list itself, a subobject disjoint from `bat_cand`
               call self%cand_cells%query_box(lo, hi, self%bat_cand, nc)
               call gather_list(self, nc)
            end if

            if (in_place) then
               call strip_vjp_pass(self, m, pb, b0, b1, b2, nc, res(:, :, s0:s1), rsub, &
                                   n_zero)
            else
               call strip_vjp_pass(self, m, pb, b0, b1, b2, nc, self%vjp_res, rsub, &
                                   n_zero)
            end if

            if (n_zero > 0) then
               ! On-nucleus columns apart, the rest re-run as a sub-strip of its own
               call nucleus_columns(m, nc, pb, self%bat_center, nuc)
               nk = 0
               do j = 1, m
                  if (nuc(j)) cycle
                  nk = nk + 1
                  lane(nk) = j
                  cp(:, nk) = pb(:, j)
                  c0(nk) = b0(j)
                  c1(:, nk) = b1(:, j)
                  c2(:, nk) = b2(:, j)
               end do
               if (nk > 0) then
                  call strip_vjp_pass(self, nk, cp, c0, c1, c2, nc, self%vjp_res, rc, &
                                      n_zero_c)
                  if (n_zero_c > 0) then
                     error stop "moist DROP LSF CFC: vjp_f1_rA_batch nucleus left"
                  end if
               end if
               rsub = .false.
               do a = 1, nk
                  j = lane(a)
                  rsub(j) = rc(a)
                  if (.not. rsub(j)) cycle
                  i = s0 - 1 + perm(q0 + j - 1)
                  if (box) then
                     call scatter_column(nc, self%bat_cand, self%vjp_res(3*nc*(a - 1) + 1), &
                                         res(:, :, i))
                  else
                     call scatter_column(nc, self%vjp_apos, self%vjp_res(3*nc*(a - 1) + 1), &
                                         res(:, :, i))
                  end if
               end do
               do j = 1, m
                  if (.not. nuc(j)) cycle
                  i = s0 - 1 + perm(q0 + j - 1)
                  call nucleus_unreached(self, points(:, i), nc)
               end do
            else if (.not. in_place) then
               do j = 1, m
                  if (.not. rsub(j)) cycle
                  i = s0 - 1 + perm(q0 + j - 1)
                  if (box) then
                     call scatter_column(nc, self%bat_cand, self%vjp_res(3*nc*(j - 1) + 1), &
                                         res(:, :, i))
                  else
                     call scatter_column(nc, self%vjp_apos, self%vjp_res(3*nc*(j - 1) + 1), &
                                         res(:, :, i))
                  end if
               end do
            end if
            reached(q0:q0 + m - 1) = rsub(1:m)
         end do

         do j = 1, np
            if (.not. reached(j)) res(:, :, s0 - 1 + perm(j)) = nan
         end do
         if (present(status)) then
            do j = 1, np
               status(s0 - 1 + perm(j)) = merge(0, 1, reached(j))
            end do
         else if (.not. all(reached(1:np))) then
            i = n + 1
            do j = 1, np
               if (.not. reached(j)) i = min(i, s0 - 1 + perm(j))
            end do
            write (col, "(i0)") i
            call fatal_error(error, "vjp_f1_rA_batch: point "//trim(col)// &
                             " is unreached, no atom within reach ("//self%name//")")
            return
         end if
      end do
   end subroutine lsf_vjp_f1_rA_batch

   !> Size the vjp strip scratch for the current atom count
   !>
   !> @param[inout] self LSF instance with `ncenters` set
   pure subroutine grow_vjp_scratch(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self

      !> Atoms the scratch must serve
      integer :: need

      if (allocated(self%vjp_pos)) then
         if (size(self%vjp_pos) < self%ncenters) deallocate (self%vjp_pos, self%vjp_apos)
      end if
      if (.not. allocated(self%vjp_pos)) then
         allocate (self%vjp_pos(max(self%ncenters, 1)), self%vjp_apos(max(self%ncenters, 1)))
      end if
      need = vjp_strip*max(self%ncenters, 1)
      if (allocated(self%vjp_res)) then
         if (size(self%vjp_res) >= 3*need) return
         deallocate (self%vjp_jet3, self%vjp_gacc, self%vjp_res)
      end if
      allocate (self%vjp_jet3(10*need), self%vjp_gacc(3*need), self%vjp_res(3*need))
   end subroutine grow_vjp_scratch

   !> Stages 1 to 3 of the nuclear vjp on one sub-strip over the gathered list
   !>
   !> - stage 1 at order 2 with the geometry cache into `strip_geo`
   !> - `n_zero > 0` (a point on a nucleus): returns after stage 1, `reached` false,
   !>   `rs` undefined
   !> - the order-2 lift decides `reached` as [[lsf_evaluate_batch]]: `PD <= 0`, or a
   !>   non-finite value, gradient or Hessian
   !> - substitute of an unreached column before the adjoint: `PD = 1`, derivative sums
   !>   zero; `lam` is then finite and the lift never divides by a non-positive `PD`
   !>   (`cfc_strip_lift_vjp` requires `PD > 0`), the column's `rs` is the caller's to
   !>   discard
   !> - the adjoint carries one more power of `1/PD` than the lift (`2/PD^3`), so it
   !>   overflows for `PD` in about 1e-154 to 1e-103 where the lift is finite: a column
   !>   whose final `rs` is not finite becomes unreached, whatever the cause
   !>
   !> Stage buffers live in fixed arrays of the strip width; the explicit-shape
   !> kernel dummies see them as `(np, ...)` by sequence association
   !>
   !> @param[inout] self    LSF instance with `bat_*` gathered
   !> @param[in]    np      Points in the sub-strip, 1..vjp_strip
   !> @param[in]    points  Sub-strip points [ndim, np]
   !> @param[in]    w0      Value weights [np]
   !> @param[in]    w1      Gradient weights [ndim, np]
   !> @param[in]    w2      Packed symmetric Hessian seeds [6, np]
   !> @param[in]    nc      Length of the shared list
   !> @param[out]   rs      Nuclear gradient per list atom and point [ndim, nc, np]
   !> @param[out]   reached Column reached at least one atom, lift and `rs` finite [np]
   !> @param[out]   n_zero  Point-atom pairs at zero distance
   subroutine strip_vjp_pass(self, np, points, w0, w1, w2, nc, rs, reached, n_zero)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Points in the sub-strip
      integer, intent(in) :: np
      !> Sub-strip points
      real(wp), intent(in) :: points(ndim, np)
      !> Value weights
      real(wp), intent(in) :: w0(np)
      !> Gradient weights
      real(wp), intent(in) :: w1(ndim, np)
      !> Packed symmetric Hessian seeds
      real(wp), intent(in) :: w2(np2, np)
      !> Length of the shared list
      integer, intent(in) :: nc
      !> Nuclear gradient per list atom and point
      real(wp), intent(out) :: rs(ndim, nc, np)
      !> Column reached at least one atom, lift and result finite
      logical, intent(out) :: reached(np)
      !> Point-atom pairs at zero distance
      integer, intent(out) :: n_zero

      !> Stage-1 sums and their adjoints, seen as `(np, ...)` by the kernel
      real(wp) :: pd0(vjp_strip), pd1(vjp_strip*ndim), pd2(vjp_strip*np2)
      real(wp) :: m0(vjp_strip), m1(vjp_strip*ndim), m2(vjp_strip*np2)
      !> Order-2 lift of the sub-strip
      real(wp) :: l0(vjp_strip), l1(ndim, vjp_strip), l2(np2, vjp_strip)
      !> Fade band in sigma
      real(wp) :: fs, fw
      !> Point index and stage 3's zero-distance count
      integer :: p, n_zero3

      call fade_band(self, fs, fw)
      associate (c => self%bat_center, r => self%bat_radius, q => self%bat_bsq, &
                 a1 => self%param%a1, a2 => self%param%a2, cp => self%param%c)
         call cfc_strip_powersums_o2(np, points, nc, c, r, q, a1, a2, cp, fs, fw, pd0, &
                                     pd1, pd2, n_zero, self%strip_geo, self%bat_sv, &
                                     self%bat_reach)
         if (n_zero > 0) then
            reached = .false.
            return
         end if

         do p = 1, np
            reached(p) = pd0(p) > 0.0_wp
            if (.not. reached(p)) pd0(p) = 1.0_wp
         end do
         call cfc_strip_lift_o2(np, pd0, pd1, pd2, l0, l1, l2)
         do p = 1, np
            if (.not. reached(p)) cycle
            reached(p) = ieee_is_finite(l0(p)) .and. all(ieee_is_finite(l1(:, p))) .and. &
                         all(ieee_is_finite(l2(:, p)))
         end do
         call clear_unreached(np, reached, pd0, pd1, pd2)

         call cfc_strip_lift_vjp(np, pd0, pd1, pd2, w0, w1, w2, m0, m1, m2)
         call cfc_strip_vjp(np, points, nc, c, r, q, a1, a2, cp, fs, fw, m0, m1, m2, &
                            self%strip_geo, self%vjp_jet3, self%vjp_gacc, self%bat_sv, &
                            self%bat_reach, rs, n_zero3)
      end associate
      n_zero = n_zero3
      ! The adjoint can overflow where the lift does not: the result decides
      do p = 1, np
         if (reached(p)) reached(p) = all(ieee_is_finite(rs(:, :, p)))
      end do
   end subroutine strip_vjp_pass

   !> Sums of the unreached columns of a sub-strip set to `PD = 1`, no derivatives
   !>
   !> @param[in]    np      Points in the sub-strip
   !> @param[in]    reached Column reached at least one atom and is finite [np]
   !> @param[inout] ps0     Order-0 sums [np]
   !> @param[inout] ps1     Order-1 sums [np, ndim]
   !> @param[inout] ps2     Order-2 sums [np, 6]
   pure subroutine clear_unreached(np, reached, ps0, ps1, ps2)
      !> Points in the sub-strip
      integer, intent(in) :: np
      !> Column reached at least one atom and is finite
      logical, intent(in) :: reached(np)
      !> Order-0 sums
      real(wp), intent(inout) :: ps0(np)
      !> Order-1 sums
      real(wp), intent(inout) :: ps1(np, ndim)
      !> Order-2 sums
      real(wp), intent(inout) :: ps2(np, np2)

      !> Point index
      integer :: p

      do p = 1, np
         if (reached(p)) cycle
         ps0(p) = 1.0_wp
         ps1(p, :) = 0.0_wp
         ps2(p, :) = 0.0_wp
      end do
   end subroutine clear_unreached

   !> Columns of a sub-strip that sit exactly on an atom of the shared list
   !>
   !> The kernel's own test, `dx*dx + dy*dy + dz*dz == 0`
   !>
   !> @param[in]  np     Points in the sub-strip
   !> @param[in]  nc     Length of the shared list
   !> @param[in]  points Sub-strip points [ndim, np]
   !> @param[in]  center Centers of the shared list [ndim, nc]
   !> @param[out] nuc    Column has a zero-distance atom [np]
   pure subroutine nucleus_columns(np, nc, points, center, nuc)
      !> Points in the sub-strip
      integer, intent(in) :: np
      !> Length of the shared list
      integer, intent(in) :: nc
      !> Sub-strip points
      real(wp), intent(in) :: points(ndim, np)
      !> Centers of the shared list
      real(wp), intent(in) :: center(ndim, *)
      !> Column has a zero-distance atom
      logical, intent(out) :: nuc(np)

      !> Point and atom index
      integer :: p, a
      !> Displacement
      real(wp) :: dx, dy, dz

      nuc = .false.
      do p = 1, np
         do a = 1, nc
            dx = points(1, p) - center(1, a)
            dy = points(2, p) - center(2, a)
            dz = points(3, p) - center(3, a)
            if (dx*dx + dy*dy + dz*dz == 0.0_wp) then
               nuc(p) = .true.
               exit
            end if
         end do
      end do
   end subroutine nucleus_columns

   !> Write one strip column into the caller's candidate positions
   !>
   !> @param[in]  nc   Length of the shared list
   !> @param[in]  apos Candidate position of list atom j [nc]
   !> @param[in]  rs   Strip column [ndim, nc]
   !> @param[out] out  Tile column, exact zero where no list atom lands [ndim, ncand]
   pure subroutine scatter_column(nc, apos, rs, out)
      !> Length of the shared list
      integer, intent(in) :: nc
      !> Candidate position of list atom j
      integer, intent(in) :: apos(nc)
      !> Strip column
      real(wp), intent(in) :: rs(ndim, nc)
      !> Tile column
      real(wp), intent(out) :: out(:, :)

      !> List atom
      integer :: j

      out = 0.0_wp
      do j = 1, nc
         out(:, apos(j)) = rs(:, j)
      end do
   end subroutine scatter_column

   !> On-nucleus column of the batched vjp, checked on the guarded per-point path
   !>
   !> - [[lsf_cfc_screen]] over the sub-strip's list at `max_deriv = 3`: the cusp makes
   !>   the per-point jet non-finite from order 1, so the column is unreached and has
   !>   no value to return
   !> - `error stop` if the point comes back reached, never the case for a point on a
   !>   nucleus
   !> - `max_deriv` restored and the per-point cache invalidated on exit
   !>
   !> @param[inout] self  LSF instance with `bat_cand(1:nc)` set
   !> @param[in]    point Evaluation point [ndim]
   !> @param[in]    nc    Length of the shared list
   subroutine nucleus_unreached(self, point, nc)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Length of the shared list
      integer, intent(in) :: nc

      !> Caller's order
      integer :: md_saved

      md_saved = self%max_deriv
      self%max_deriv = 3
      ! `bat_cand` aliases a component of `self`: legal as the screen never touches it
      call lsf_cfc_screen(self, point, self%bat_cand(1:nc))
      if (.not. self%prepared_unreached) then
         error stop "moist DROP LSF CFC: vjp_f1_rA_batch on-nucleus column reached"
      end if
      self%max_deriv = md_saved
      self%prepared_deriv = -1
   end subroutine nucleus_unreached

   !* ================================================================================= *!
   !*                          One-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Mixed third derivative d^3S / (dr^2 dR_A) and its lower orders
   !>
   !> All three outputs are active-indexed: slot `i` belongs to `active_atom(i)`
   !> The order requirement is 3, not 2: the kernel reads `qn2_rr`, whose total
   !> derivative order is three
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    dS/dR_A [3, >= active_count()] (optional)
   !> @param[out] lsf2_r_rA  d^2S/(dr dR_A) [3, 3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rA d^3S/(dr^2 dR_A) [6, 3, >= active_count()], rr packed
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear gradient
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      !> Mixed third derivative, rr packed
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :)

      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Active-list index and nuclear component
      integer :: ia, b

      call self%require_deriv(3, "f3_rr_rA")

      do ia = 1, self%n_active
         call cfc_nuclear_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                               self%qn0(:, ia), self%qn1_r(:, :, ia), &
                               self%qn2_rr(:, :, :, ia), self%qn3_rrr(:, :, :, :, ia), &
                               2, f1_rA, f2_r_rA, f3_rr_rA, d4)
         if (present(lsf1_rA)) lsf1_rA(:, ia) = f1_rA
         if (present(lsf2_r_rA)) lsf2_r_rA(:, :, ia) = f2_r_rA
         do b = 1, ndim
            call pack_sym_matrix(f3_rr_rA(:, :, b), lsf3_rr_rA(:, b, ia), ndim)
         end do
      end do
   end subroutine lsf_f3_rr_rA

   !> Mixed fourth derivative d^4S / (dr^3 dR_A), active-indexed
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf4_rrr_rA d^4S/(dr^3 dR_A) [10, 3, >= active_count()], rrr packed
   subroutine lsf_f4_rrr_rA(self, lsf4_rrr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative, rrr packed
      real(wp), intent(out) :: lsf4_rrr_rA(:, :, :)

      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: f4_rrr_rA(ndim, ndim, ndim, ndim)
      !> Active-list index and nuclear component
      integer :: ia, b

      call self%require_deriv(4, "f4_rrr_rA")

      do ia = 1, self%n_active
         call cfc_nuclear_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                               self%qn0(:, ia), self%qn1_r(:, :, ia), &
                               self%qn2_rr(:, :, :, ia), self%qn3_rrr(:, :, :, :, ia), &
                               3, f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
         do b = 1, ndim
            call pack_sym_tensor3(f4_rrr_rA(:, :, :, b), lsf4_rrr_rA(:, b, ia), ndim)
         end do
      end do
   end subroutine lsf_f4_rrr_rA

   !* ================================================================================= *!
   !*                               Radius derivatives                                  *!
   !* ================================================================================= *!
   !
   ! The radii are parameters of the level set, so a model tying them to the
   ! geometry differentiates through them as well, and three things separate
   ! this channel from the `_rA` one, all because `R_a` is a scalar parameter
   ! rather than a coordinate:
   !
   ! - no sign flip
   ! - no trailing derivative index, the rank equalling the spatial order
   ! - no derivative order consumed
   !
   ! Unlike `qn*`, the `rd*` family is *not* cached by `prepare`; it could be,
   ! depending on nothing but the geometry, but that would charge every
   ! fixed-radius caller an extra O(n_active**2) pair sweep at every grid point
   ! for a family it never reads, so it is accumulated on demand instead, the
   ! same way `tg*` is

   !> Accumulate the per-atom radius-derivative pseudo-density tensors
   !>
   !> The radius counterpart of the `qn*` sweep in [[lsf_cfc_screen]], run on
   !> demand, with slot `i` of each output belonging to `active_atom(i)`
   !>
   !> `rd3_rrr` is optional because [[cfc_radius_hvp]] needs the family only up to
   !> order two and would otherwise have to materialize an unread `27 * n_active`
   !> buffer to satisfy the kernel signature:
   !>
   !> - omitting it sinks the order-3 accumulation into two stack rows instead
   !> - the loops are spelled out twice rather than once behind a pointer, the
   !>   kernel arguments being `intent(inout)` accumulators taken one atom at a
   !>   time
   !> - the only single-body alternative is an aliasing-opaque pointer with
   !>   `merge`-computed slot indices in the pair loop
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
                                        self%act_fade_start, self%act_fade_width, &
                                        level, rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                        rd3_rrr(:, :, :, ia))
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_radius_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                         self%act_radius(ia), self%act_radius(ib), &
                                         self%param%a2, self%param%c, &
                                         self%act_fade_start, self%act_fade_width, level, &
                                         rd0(ia), rd0(ib), &
                                         rd1_r(:, ia), rd1_r(:, ib), &
                                         rd2_rr(:, :, ia), rd2_rr(:, :, ib), &
                                         rd3_rrr(:, :, :, ia), rd3_rrr(:, :, :, ib))
            end do
         end do
      else
         ! Zeroed once, not per atom: the sinks are `intent(inout)` accumulators
         ! and reading them undefined would trap under the debug FP settings
         sk3_a = 0.0_wp
         sk3_b = 0.0_wp

         do ia = 1, self%n_active
            call cfc_atomic_radius_eval(self%act_d(:, ia), self%act_radius(ia), &
                                        self%param%a1, &
                                        self%act_fade_start, self%act_fade_width, &
                                        level, rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                        sk3_a)
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_radius_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                         self%act_radius(ia), self%act_radius(ib), &
                                         self%param%a2, self%param%c, &
                                         self%act_fade_start, self%act_fade_width, level, &
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
   !> All three outputs are active-indexed: slot `i` belongs to `active_atom(i)`
   !> `R_a` is the *radius* of that atom, so the results carry no trailing
   !> derivative index and the rank equals the spatial order
   !>
   !> The order requirement is 2, not 3: `d/dR_a` consumes no derivative order,
   !> so the deepest input is `rd2_rr`, whose total order is two -- one less than
   !> `lsf_f3_rr_rA` needs for the same spatial order
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf1_rad    dS/dR_a [>= active_count()] (optional)
   !> @param[out] lsf2_r_rad  d^2S/(dr dR_a) [3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rad d^3S/(dr^2 dR_a) [6, >= active_count()], rr packed
   subroutine lsf_f3_rr_rad(self, lsf1_rad, lsf2_r_rad, lsf3_rr_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Radius gradient
      real(wp), intent(out), optional :: lsf1_rad(:)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rad(:, :)
      !> Mixed third derivative, rr packed
      real(wp), intent(out) :: lsf3_rr_rad(:, :)

      !> Per-atom radius tensors
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rd3_rrr(:, :, :, :)
      !> Kernel outputs of one atom
      real(wp) :: f1_rad, f2_r_rad(ndim), f3_rr_rad(ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

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
         call pack_sym_matrix(f3_rr_rad, lsf3_rr_rad(:, ia), ndim)
      end do
   end subroutine lsf_f3_rr_rad

   !> Normalized level set S/||grad S|| and its nuclear gradient
   !>
   !> The normalization itself is level-set agnostic -- it only combines outputs
   !> the kernel already returns -- so it is spelled out here rather than
   !> generated, matching `svdw_normalized_eval` term for term including the
   !> degenerate-gradient guard
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
      real(wp) :: f0, f1_r(ndim)
      !> One-nuclear-index outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), d3n(ndim, ndim, ndim)
      real(wp) :: d4n(ndim, ndim, ndim, ndim)
      !> Gradient norm, its reciprocal powers and the gradient/mixed contraction
      real(wp) :: g, inv_g, inv_g3, gd(ndim)
      !> Active-list index and Cartesian index
      integer :: ia, j

      val = 0.0_wp
      call self%require_deriv(2, "normalized_f01_rA")

      f0 = self%jet0
      f1_r = self%jet1_r

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
   !*                                 Atomic partition                                  *!
   !* ================================================================================= *!

   !> CFC splits into atomic shares
   !>
   !> @param[in] self LSF instance
   pure function lsf_has_atom_partition(self) result(available)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Always available
      logical :: available

      available = .true.
   end function lsf_has_atom_partition

   !> Add the atomic shares of PD at the prepared point
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    scale Weight of the point
   !> @param[inout] acc   Per-atom accumulator, user space [ncenters]
   subroutine lsf_atom_partition(self, scale, acc)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self
      !> Weight of the point
      real(wp), intent(in) :: scale
      !> Per-atom accumulator
      real(wp), intent(inout) :: acc(:)

      !> Atomic shares of PD, active-indexed
      real(wp) :: pd_atom(self%n_active)
      !> Value of one term, and its unused derivative accumulators
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: t4(ndim, ndim, ndim, ndim)
      integer :: i, j

      call self%require_deriv(0, "atom_partition")
      t1 = 0.0_wp
      t2 = 0.0_wp
      t3 = 0.0_wp
      t4 = 0.0_wp

      do i = 1, self%n_active
         t0 = 0.0_wp
         call cfc_atomic_term_eval(self%act_d(:, i), self%act_radius(i), &
                                   self%param%a1, self%act_fade_start, self%act_fade_width, &
                                   0, t0, t1, t2, t3, t4)
         pd_atom(i) = t0
      end do
      do i = 1, self%n_active
         do j = i + 1, self%n_active
            t0 = 0.0_wp
            call cfc_pair_spatial_eval(self%act_d(:, i), self%act_d(:, j), &
                                       self%act_radius(i), self%act_radius(j), &
                                       self%param%a2, self%param%c, &
                                       self%act_fade_start, self%act_fade_width, &
                                       0, t0, t1, t2, t3, t4)
            pd_atom(i) = pd_atom(i) + 0.5_wp*t0
            pd_atom(j) = pd_atom(j) + 0.5_wp*t0
         end do
      end do
      pd_atom = pd_atom*(scale/sum(pd_atom))
      do i = 1, self%n_active
         acc(self%act_atom(i)) = acc(self%act_atom(i)) + pd_atom(i)
      end do
   end subroutine lsf_atom_partition

   !* ================================================================================= *!
   !*                          Two-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Fill the uncontracted two-nucleus family at one spatial order
   !>
   !> Shared body of the three `f*_rArB` accessors, in two sweeps:
   !>
   !>   1. the *diagonal* `qq` blocks, one per active atom: atom A's own term
   !>      plus the `aa` / `bb` block of every pair containing A. One pass over
   !>      the atoms and one over the unordered pairs
   !>   2. the ordered-pair loop, which for `A /= B` gets its `qq` block from the
   !>      single pair term containing both centers (no other term survives two
   !>      derivatives with respect to two different nuclei) and for `A == B`
   !>      reuses the diagonal accumulated in pass one, then lifts
   !>
   !> The `qq` scratch is a local allocation rather than a persistent buffer for
   !> the same reason SvdW's pair cache is: only these three accessors need it,
   !> none of them sits on the projection hot path, and a persistent order-2
   !> two-nucleus cache would cost 936 bytes per atom per thread
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
      !> Mixed fourth derivative, rr packed
      real(wp), intent(out), optional :: lsf4_rr_rArB(:, :, :, :, :)

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
      integer :: ia, ib, n, s, t

      n = self%n_active

      allocate (qd0(ndim, ndim, n), source=0.0_wp)
      allocate (qd1(ndim, ndim, ndim, n), source=0.0_wp)
      allocate (qd2(ndim, ndim, ndim, ndim, n), source=0.0_wp)

      do ia = 1, n
         call cfc_atomic_hessian_eval(self%act_d(:, ia), self%act_radius(ia), &
                                      self%param%a1, &
                                      self%act_fade_start, self%act_fade_width, &
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
                                       self%param%a2, self%param%c, &
                                       self%act_fade_start, self%act_fade_width, level, &
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
                                          self%param%a2, self%param%c, &
                                          self%act_fade_start, self%act_fade_width, level, &
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
                     call pack_sym_matrix(b4(:, :, s, t), &
                                          lsf4_rr_rArB(:, s, ia, t, ib), ndim)
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

      call self%require_deriv(2, "f3_r_rArB")

      call hessian_family(self, 1, lsf3_r_rArB=lsf3_r_rArB)
   end subroutine lsf_f3_r_rArB

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf4_rr_rArB d^4S/(dr^2 dR_A dR_B) [6, 3, >= n_act, 3, >= n_act], rr packed
   subroutine lsf_f4_rr_rArB(self, lsf4_rr_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative, rr packed
      real(wp), intent(out) :: lsf4_rr_rArB(:, :, :, :, :)

      call self%require_deriv(3, "f4_rr_rArB")

      call hessian_family(self, 2, lsf4_rr_rArB=lsf4_rr_rArB)
   end subroutine lsf_f4_rr_rArB

   !* ================================================================================= *!
   !*                     Two-radius and nuclear-radius derivatives                     *!
   !* ================================================================================= *!

   !> Fill the uncontracted two-radius family at one spatial order
   !>
   !> Shared body of the three `f*_radrad` accessors; see the section note above
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
      !> Mixed fourth derivative, rr packed
      real(wp), intent(out), optional :: lsf4_rr_radrad(:, :, :)

      !> Per-atom radius tensors, in both the A and the B slot of the lift
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rd3_rrr(:, :, :, :)
      !> Diagonal two-radius accumulators, one block per active atom
      real(wp), allocatable :: dd0(:), dd1(:, :), dd2(:, :, :)
      !> Cross block of every unordered pair, cached by the single sweep below
      !> Only the orders `level` will actually lift are allocated
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
                                      self%act_fade_start, self%act_fade_width, &
                                      level, dd0(ia), dd1(:, ia), dd2(:, :, ia))
      end do
      ip = 0
      do ia = 1, n
         do ib = ia + 1, n
            ip = ip + 1
            x0 = 0.0_wp
            x1 = 0.0_wp
            x2 = 0.0_wp
            call cfc_pair_radius2_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                       self%act_radius(ia), self%act_radius(ib), &
                                       self%param%a2, self%param%c, &
                                       self%act_fade_start, self%act_fade_width, level, &
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
            ! Symmetric in its two radii, so one lift fills both ordered slots
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
      !> Mixed fourth derivative, rr packed
      real(wp), intent(inout), optional :: lsf4_rr_radrad(:, :, :)

      if (present(lsf2_radrad)) lsf2_radrad(ia, ib) = b2
      if (present(lsf3_r_radrad)) lsf3_r_radrad(:, ia, ib) = b3
      if (present(lsf4_rr_radrad)) then
         call pack_sym_matrix(b4, lsf4_rr_radrad(:, ia, ib), ndim)
      end if
   end subroutine radius2_store

   !> Fill the uncontracted nuclear-radius family at one spatial order
   !>
   !> Shared body of the three `f*_rA_rad` accessors
   !>
   !> - the first atom index carries the position derivative and the second the
   !>   radius one, so the result is *not* symmetric under exchanging them
   !> - each unordered pair is visited once, its `ab` and `ba` kernel blocks
   !>   filling the two ordered slots
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
      !> Mixed fourth derivative, rr packed
      real(wp), intent(out), optional :: lsf4_rr_rA_rad(:, :, :, :)

      !> Per-atom radius tensors, filling the lift's B slot
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rd3_rrr(:, :, :, :)
      !> Diagonal nuclear-radius accumulators, one block per active atom
      real(wp), allocatable :: dd0(:, :), dd1(:, :, :), dd2(:, :, :, :)
      !> Both cross blocks of every unordered pair, cached by the single sweep
      !> below; only the orders `level` will actually lift are allocated
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
                                     self%act_fade_start, self%act_fade_width, &
                                     level, dd0(:, ia), dd1(:, :, ia), dd2(:, :, :, ia))
      end do
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
                                      self%param%a2, self%param%c, &
                                      self%act_fade_start, self%act_fade_width, level, &
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
            ! `ab` is position on ia and radius on ib; `ba` is the other way
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
      !> Mixed fourth derivative, rr packed
      real(wp), intent(inout), optional :: lsf4_rr_rA_rad(:, :, :, :)

      !> Nuclear component
      integer :: c

      if (present(lsf2_rA_rad)) lsf2_rA_rad(:, ia, ib) = b2
      if (present(lsf3_r_rA_rad)) lsf3_r_rA_rad(:, :, ia, ib) = b3
      if (present(lsf4_rr_rA_rad)) then
         do c = 1, ndim
            call pack_sym_matrix(b4(:, :, c), lsf4_rr_rA_rad(:, c, ia, ib), ndim)
         end do
      end if
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

      ! One order below `f2_rArB`: neither retained radius consumes one
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

      call self%require_deriv(1, "f3_r_radrad")

      call radius2_family(self, 1, lsf3_r_radrad=lsf3_r_radrad)
   end subroutine lsf_f3_r_radrad

   !> Mixed fourth derivative d^4S / (dr^2 dR_a dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b) [6, >= n_act, >= n_act], rr packed
   subroutine lsf_f4_rr_radrad(self, lsf4_rr_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative, rr packed
      real(wp), intent(out) :: lsf4_rr_radrad(:, :, :)

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

      ! The retained *position* still costs an order while the radius does not,
      ! so this matches `f2_rArB` rather than sitting one above it
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

      call self%require_deriv(2, "f3_r_rA_rad")

      call nucrad_family(self, 1, lsf3_r_rA_rad=lsf3_r_rA_rad)
   end subroutine lsf_f3_r_rA_rad

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b) [6, 3, >= n_act, >= n_act], rr packed
   subroutine lsf_f4_rr_rA_rad(self, lsf4_rr_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Mixed fourth derivative, rr packed
      real(wp), intent(out) :: lsf4_rr_rA_rad(:, :, :, :)

      call self%require_deriv(3, "f4_rr_rA_rad")

      call nucrad_family(self, 2, lsf4_rr_rA_rad=lsf4_rr_rA_rad)
   end subroutine lsf_f4_rr_rA_rad

   !* ================================================================================= *!
   !*                     Direction-contracted nuclear derivatives                      *!
   !* ================================================================================= *!

   !> Accumulate the direction-contracted pseudo-density tensors
   !>
   !> The `tg*` family depends on the caller's displacement field, so `prepare`
   !> cannot own it and this is its atom-plus-pair sweep
   !>
   !> - used only by the `tangent_*` accessors
   !> - the `hvp_*` ones get the same family for free from `cfc_*_hvp_eval` and
   !>   must never call this as well, or every contribution is counted twice
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
                                      self%act_fade_start, self%act_fade_width, &
                                      v(:, self%act_atom(ia)), level, &
                                      tg0, tg1_r, tg2_rr, tg3_rrr)
      end do
      do ia = 1, self%n_active
         do ib = ia + 1, self%n_active
            call cfc_pair_tangent_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                       self%act_radius(ia), self%act_radius(ib), &
                                       self%param%a2, self%param%c, &
                                       self%act_fade_start, self%act_fade_width, &
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
   !> @param[out] res  sum_B v_B . d/dR_B of d^2S/dr^2, packed [6]
   subroutine lsf_tangent_f2_rr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted Hessian derivative, packed
      real(wp), intent(out) :: res(:)

      !> Contracted pseudo-density tensors and kernel outputs
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), d3(ndim, ndim, ndim)

      call self%require_deriv(2, "tangent_f2_rr")

      call tangent_tensors(self, v, 2, tg0, tg1_r, tg2_rr, tg3_rrr)
      call cfc_tangent_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            tg0, tg1_r, tg2_rr, tg3_rrr, 2, t0, t1, t2, d3)
      call pack_sym_matrix(t2, res, ndim)
   end subroutine lsf_tangent_f2_rr

   !> Directional nuclear derivative of the third spatial derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^3S/dr^3, packed [10]
   subroutine lsf_tangent_f3_rrr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted third-derivative derivative, packed
      real(wp), intent(out) :: res(:)

      !> Contracted pseudo-density tensors and kernel outputs
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)

      call self%require_deriv(3, "tangent_f3_rrr")

      call tangent_tensors(self, v, 3, tg0, tg1_r, tg2_rr, tg3_rrr)
      call cfc_tangent_eval(self%pd0, self%pd1_r, self%pd2_rr, self%pd3_rrr, &
                            tg0, tg1_r, tg2_rr, tg3_rrr, 3, t0, t1, t2, t3)
      call pack_sym_tensor3(t3, res, ndim)
   end subroutine lsf_tangent_f3_rrr

   !> Fill the per-atom Hessian-vector-product family (and its tangent ladder)
   !>
   !> One atom sweep and one pair sweep of `cfc_*_hvp_eval`, which fills the
   !> aggregate `tg*` accumulators and the per-atom `hvp*` ones together
   !>
   !> - the tangent ladder comes out of the same derivative structure, so it is
   !>   free here, and is exactly why [[tangent_tensors]] must not also run
   !>
   !> Supplying `vrad` promotes the contraction from `sum_B v_B . d/dR_B` to
   !> `sum_B (v_B . d/dR_B + vr_B d/dRad_B)`:
   !>
   !> - a second kernel sweep adds the radius half to `tg*` and `hvp*` and fills
   !>   `rh*` in the same pass
   !> - `vrad` is the switch, so the three `rh*` buffers must be present with it,
   !>   while supplying them without it is harmless (they come back zeroed)
   !> - without `vrad` nothing changes, down to the bit, the position-only sweep
   !>   being untouched
   !>
   !> The `hv*` family is optional for the mirror-image reason:
   !>
   !> - [[cfc_radius_hvp]] wants the joint `tg*` and `rh*` but not the position
   !>   row, and materializing `39 * n_active` doubles to throw away is the one
   !>   part of that trade worth avoiding
   !> - absent, the row accumulates into stack scratch instead
   !> - that is what doubles the sweeps below into `present`-guarded pairs, a
   !>   single body needing an aliasing-opaque pointer with `merge`-computed slot
   !>   indices, worse in the `O(n_active**2)` pair loop than a duplicated call
   !> - `hv0`, `hv1` and `hv2` are all-or-none
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
         ! read uninitialized stack and trap under the debug FP settings
         sk_a0 = 0.0_wp
         sk_a1 = 0.0_wp
         sk_a2 = 0.0_wp
         sk_b0 = 0.0_wp
         sk_b1 = 0.0_wp
         sk_b2 = 0.0_wp
      end if
      ! Zeroed before the `vrad` gate, not after it: a caller may hand over the
      ! buffers without asking for the radius channel, and an untouched
      ! `intent(out)` array is the one thing it could not inspect safely
      if (present(rh0)) rh0 = 0.0_wp
      if (present(rh1_r)) rh1_r = 0.0_wp
      if (present(rh2_rr)) rh2_rr = 0.0_wp

      if (want_hv) then
         do ia = 1, self%n_active
            call cfc_atomic_hvp_eval(self%act_d(:, ia), self%act_radius(ia), &
                                     self%param%a1, &
                                     self%act_fade_start, self%act_fade_width, &
                                     v(:, self%act_atom(ia)), level, &
                                     tg0, tg1_r, tg2_rr, tg3_rrr, &
                                     hv0(:, ia), hv1(:, :, ia), hv2(:, :, :, ia))
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                      self%act_radius(ia), self%act_radius(ib), &
                                      self%param%a2, self%param%c, &
                                      self%act_fade_start, self%act_fade_width, &
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
                                     self%act_fade_start, self%act_fade_width, &
                                     v(:, self%act_atom(ia)), level, &
                                     tg0, tg1_r, tg2_rr, tg3_rrr, &
                                     sk_a0, sk_a1, sk_a2)
         end do
         do ia = 1, self%n_active
            do ib = ia + 1, self%n_active
               call cfc_pair_hvp_eval(self%act_d(:, ia), self%act_d(:, ib), &
                                      self%act_radius(ia), self%act_radius(ib), &
                                      self%param%a2, self%param%c, &
                                      self%act_fade_start, self%act_fade_width, &
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
      ! a caller turning the radius pass on without handing over all three
      ! buffers would write through an absent dummy -- refuse instead
      if (.not. (present(rh0) .and. present(rh1_r) .and. present(rh2_rr))) then
         error stop "moist DROP CFC LSF: hvp_tensors needs rh0, rh1_r and "// &
            "rh2_rr whenever vrad is supplied"
      end if

      if (want_hv) then
         do ia = 1, self%n_active
            call cfc_atomic_radius_hvp_eval(self%act_d(:, ia), self%act_radius(ia), &
                                            self%param%a1, &
                                            self%act_fade_start, self%act_fade_width, &
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
                                             self%act_fade_start, self%act_fade_width, &
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
                                            self%act_fade_start, self%act_fade_width, &
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
                                             self%act_fade_start, self%act_fade_width, &
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
   !> @param[out] res  sum_B v_B . d^4S/(dr^2 dR_A dR_B) [6, 3, >= active_count()], rr packed
   subroutine lsf_hvp_f3_rr_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed fourth derivative, rr packed
      real(wp), intent(out) :: res(:, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom HVP tensors and the radius-row by-product
      real(wp), allocatable :: hv0(:, :), hv1(:, :, :), hv2(:, :, :, :)
      real(wp), allocatable :: rh0(:), rh1_r(:, :), rh2_rr(:, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), h3(ndim, ndim, ndim)
      !> Active-list index and nuclear component
      integer :: ia, b

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
         do b = 1, ndim
            call pack_sym_matrix(h3(:, :, b), res(:, b, ia), ndim)
         end do
      end do
   end subroutine lsf_hvp_f3_rr_rA

   !* ================================================================================= *!
   !*                     Radius row of the joint Hessian-vector product                *!
   !* ================================================================================= *!

   !> Allocate the per-atom radius-row buffers
   !>
   !> Only ever called on a path that supplies a radius direction, the rows
   !> existing only for a joint `(v, vrad)` contraction:
   !>
   !> - live output for [[cfc_radius_hvp]]
   !> - write-only scratch for the three `hvp_*_rA` entry points, which need
   !>   `hvp_tensors` to have somewhere to put the row it fills alongside the
   !>   `tg*` they do read
   !> - a pure-nuclear contraction never runs the radius sweep at all and must
   !>   not pay for these, hence the `present` guard at those three call sites
   !>   rather than an unconditional allocation
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
   !> @param[out] res  Contracted mixed fourth derivative [6, >= active_count()], rr packed
   subroutine lsf_hvp_f3_rr_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed fourth derivative, rr packed
      real(wp), intent(out) :: res(:, :)

      call cfc_radius_hvp(self, v, vrad, 2, "hvp_f3_rr_rad", res3=res)
   end subroutine lsf_hvp_f3_rr_rad

   !> Shared driver of the three radius-row Hessian-vector products
   !>
   !> The three public entry points differ only in `max_deriv` and in which output
   !> they keep, and unlike the `_rA` ladder there is no reason to spell the loop
   !> out three times, the radius channel having one code path
   !>
   !> - exactly one of `res1`/`res2`/`res3` must be present, matching `max_deriv`
   !>
   !> Two layer-1 sweeps feed it:
   !>
   !> - `hvp_tensors` with a radius direction gives the joint `tg*` and the `rh*`
   !>   family, `radius_tensors` gives `rd*` up to order two
   !> - both also produce a row this driver has no use for -- `hvp*`, belonging to
   !>   the *position* row, and `rd3_rrr` -- and both are asked to sink it into
   !>   stack scratch by omitting the buffer
   !> - the work itself is not avoidable and not worth avoiding, the generated
   !>   kernel producing the two rows' radius terms from one shared CSE block,
   !>   cheaper than computing them apart even counting the waste (6111
   !>   temporaries against 4427 + 2230)
   !> - only the `66 * n_active` doubles of throwaway storage were worth avoiding
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  v         Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad      Radius directions [ncenters]
   !> @param[in]  max_deriv Highest spatial-derivative order (0..2)
   !> @param[in]  caller    Name used by `require_deriv` on failure
   !> @param[out] res1      Order-0 result [>= active_count()] (optional)
   !> @param[out] res2      Order-1 result [3, >= active_count()] (optional)
   !> @param[out] res3      Order-2 result [6, >= active_count()], packed (optional)
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
      !> Order-2 result, packed
      real(wp), intent(out), optional :: res3(:, :)

      !> Contracted pseudo-density tensors
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      !> Per-atom radius and radius-row tensors
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      real(wp), allocatable :: rh0(:), rh1_r(:, :), rh2_rr(:, :, :)
      !> Kernel outputs of one atom
      real(wp) :: h1, h2(ndim), h3(ndim, ndim)
      !> Active-list index
      integer :: ia

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
         if (present(res3)) call pack_sym_matrix(h3, res3(:, ia), ndim)
      end do
   end subroutine cfc_radius_hvp

   !* ================================================================================= *!
   !*                  Jet-contracted nuclear vector-Jacobian product                   *!
   !* ================================================================================= *!

   !> Adjoint jet contracted onto the nuclear gradient, active-indexed
   !>
   !> Returns, for every active atom `i`,
   !>
   !>    res(beta, i) = w0 f1_rA(beta, i)
   !>                 + sum_a w1(a) f2_r_rA(a, beta, i)
   !>                 + sum_a sum_b w2(a, b) f3_rr_rA(a, b, beta, i)
   !>
   !> i.e. the `f1_rA` rung with the jet indices summed away, just as
   !> `hvp_f1_rA` is that same rung contracted with a nuclear direction instead
   !>
   !> - `w2` is the packed symmetric seed [6], off-diagonals not doubled: each
   !>   off-diagonal slot stands for both (a, b) and (b, a), the same convention
   !>   as `vjp_f1_rA_batch`
   !> - slot `i` belongs to `active_atom(i)`, columns past `active_count()` are
   !>   left untouched, and nothing is written at all when no atom is active
   !> - the order requirement is 3, not 2, for the same reason as `f3_rr_rA`: the
   !>   kernel reads `qn2_rr`, of total derivative order three, and `prepare`
   !>   fills that family only one spatial order below `max_deriv`
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of `f1_rA`
   !> @param[in]  w1   Adjoint weights of `f2_r_rA` [3]
   !> @param[in]  w2   Adjoint weights of `f3_rr_rA`, packed symmetric [6]
   !> @param[out] res  Contracted nuclear gradient [3, >= active_count()]
   subroutine lsf_vjp_f1_rA(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Adjoint weight of `f1_rA`
      real(wp), intent(in) :: w0
      !> Adjoint weights of `f2_r_rA`
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of `f3_rr_rA`, packed symmetric
      real(wp), intent(in) :: w2(:)
      !> Contracted nuclear gradient
      real(wp), intent(out) :: res(:, :)

      !> Kernel output of one atom
      real(wp) :: vjp1(ndim)
      !> Seed unpacked to a symmetric matrix, both off-diagonal slots = w2(slot)
      real(wp) :: wsym(ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(3, "vjp_f1_rA")

      call unpack_sym_matrix(w2, wsym, ndim)
      do ia = 1, self%n_active
         call cfc_vjp_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                           self%qn0(:, ia), self%qn1_r(:, :, ia), &
                           self%qn2_rr(:, :, :, ia), w0, w1, wsym, vjp1)
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
   !> because `R_a` is a radius; the jet itself is unchanged. `w2` is the packed
   !> symmetric seed [6], off-diagonals not doubled: each off-diagonal slot stands
   !> for both (a, b) and (b, a), the same convention as `vjp_f1_rA`
   !>
   !> Slot `i` belongs to `active_atom(i)`, entries past `active_count()` are left
   !> untouched, and nothing is written at all when no atom is active
   !>
   !> - the order requirement is 2, one less than `vjp_f1_rA` needs, for the same
   !>   reason `f3_rr_rad` needs one less than `f3_rr_rA`
   !> - `d/dR_a` consumes no derivative order, so the deepest input is `rd2_rr` at
   !>   total order two
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of `f1_rad`
   !> @param[in]  w1   Adjoint weights of `f2_r_rad` [3]
   !> @param[in]  w2   Adjoint weights of `f3_rr_rad`, packed symmetric [6]
   !> @param[out] res  Contracted radius gradient [>= active_count()]
   subroutine lsf_vjp_f1_rad(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Adjoint weight of `f1_rad`
      real(wp), intent(in) :: w0
      !> Adjoint weights of `f2_r_rad`
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of `f3_rr_rad`, packed symmetric
      real(wp), intent(in) :: w2(:)
      !> Contracted radius gradient
      real(wp), intent(out) :: res(:)

      !> Per-atom radius tensors
      real(wp), allocatable :: rd0(:), rd1_r(:, :), rd2_rr(:, :, :)
      !> Kernel output of one atom
      real(wp) :: vjp1
      !> Seed unpacked to a symmetric matrix, both off-diagonal slots = w2(slot)
      real(wp) :: wsym(ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(2, "vjp_f1_rad")

      allocate (rd0(self%n_active))
      allocate (rd1_r(ndim, self%n_active))
      allocate (rd2_rr(ndim, ndim, self%n_active))
      call radius_tensors(self, 2, rd0, rd1_r, rd2_rr)

      call unpack_sym_matrix(w2, wsym, ndim)
      do ia = 1, self%n_active
         call cfc_radius_vjp_eval(self%pd0, self%pd1_r, self%pd2_rr, &
                                  rd0(ia), rd1_r(:, ia), rd2_rr(:, :, ia), &
                                  w0, w1, wsym, vjp1)
         res(ia) = vjp1
      end do
   end subroutine lsf_vjp_f1_rad

   !* ================================================================================= *!
   !*                                 Screening offset                                  *!
   !* ================================================================================= *!

   !> Radial offset from the atom surface where the CFC contribution hits the threshold
   !>
   !> The pseudo-density gets two contributions from atom `a`, both decaying in the
   !> reduced distance `s_a = 1 + delta/R_a`, so the offset scales with the radius
   !>
   !> - fade on: the band end `R_a (fade_start + fade_width)/|a2|`, where every
   !>   term of the atom is exactly zero; `screening_threshold` is not read
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

      if (self%param%fade_enabled) then
         offset = huge(0.0_wp)
         if (self%param%a2 >= 0.0_wp) return
         offset = radius*(self%param%fade_start + self%param%fade_width)/abs(self%param%a2)
         ! An invalid band (refused by `prepare` and `evaluate_batch`) screens nothing
         if (.not. (offset >= 0.0_wp .and. offset < huge(0.0_wp))) offset = huge(0.0_wp)
         return
      end if

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

      ! Exact at m = kernel_pair_power only, the kernel's pair term is (1 - q)**4
      two_pow_m = 2.0_wp**self%param%m
      cap = self%param%c*two_pow_m
      delta_pair = 2.0_wp*radius*(log_thr - log(cap))/self%param%a2

      offset = max(0.0_wp, delta_atomic, delta_pair)
   end function lsf_screening_offset

   !* ================================================================================= *!
   !*                        Surface-free exclusion certificate                         *!
   !* ================================================================================= *!

   !> Cache the global Lipschitz bound of the CFC level set for the bound geometry
   !>
   !> - zero (certify nothing) for `a1 >= 0`, `a2 >= 0`, `c < 0`, odd `m`,
   !>   `m` /= [[kernel_pair_power]], a non-positive radius, or a bound that
   !>   is not finite
   !>
   !> @param[inout] self LSF instance (reads `param`, `ncenters`, `radii`)
   subroutine cfc_set_exclusion_bound(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self

      !> Bracket of the `log eps` bisection; the crossing sits near `log eps = 0`
      !> for any sane parameterization, and this span reaches it even for absurd
      !> ones
      real(wp), parameter :: eps_log_span = 40.0_wp
      !> Bisection steps, closing the bracket to far below roundoff, with the
      !> whole search running once per `update`
      integer, parameter :: eps_steps = 80

      !> Exponent magnitudes and the tail-domination margin `beta - alpha/2`
      real(wp) :: alpha, beta, gam
      !> The two smallest radii in the structure
      real(wp) :: rmin, rmin2
      !> Pair power as a real, its Young conjugate `m/(m-1)`, and `log B`
      real(wp) :: mreal, conj, log_b
      !> The atomic and pair brackets, and the bound itself
      real(wp) :: lip_atom, lip_pair, lip
      !> Bisection bracket and midpoint, in `log eps`
      real(wp) :: lo, hi, mid
      !> Atom index and bisection counter
      integer :: i, it

      self%excl_inv_lip = 0.0_wp
      if (self%ncenters < 1) return
      if (.not. allocated(self%radii)) return
      if (size(self%radii) < self%ncenters) return

      alpha = -self%param%a1
      beta = -self%param%a2
      mreal = real(self%param%m, wp)

      if (alpha <= 0.0_wp .or. beta <= 0.0_wp) return
      if (self%param%c < 0.0_wp) return
      if (self%param%m < 2 .or. modulo(self%param%m, 2) /= 0) return
      ! The kernel evaluates m = 4 whatever `param%m` says
      if (self%param%m /= kernel_pair_power) return

      rmin = huge(0.0_wp)
      rmin2 = huge(0.0_wp)
      do i = 1, self%ncenters
         if (.not. (self%radii(i) > 0.0_wp)) return
         if (self%radii(i) < rmin) then
            rmin2 = rmin
            rmin = self%radii(i)
         else if (self%radii(i) < rmin2) then
            rmin2 = self%radii(i)
         end if
      end do

      lip_atom = alpha/rmin
      lip = lip_atom

      if (self%ncenters >= 2 .and. self%param%c > 0.0_wp) then
         lip_pair = beta*(1.0_wp/rmin + 1.0_wp/rmin2)
         gam = beta - 0.5_wp*alpha
         if (gam <= 0.0_wp) return

         log_b = log(self%param%c) + (2.0_wp*beta - alpha) &
                 + mreal*(log(mreal/gam) - 1.0_wp) &
                 + log(0.5_wp*real(self%ncenters - 1, wp)) &
                 - mreal*log(rmin)
         conj = mreal/(mreal - 1.0_wp)

         lo = -eps_log_span
         hi = eps_log_span
         do it = 1, eps_steps
            mid = 0.5_wp*(lo + hi)
            if (lip_atom + exp(log_b + mreal*mid) < &
                lip_pair + (mreal - 1.0_wp)*exp(-conj*mid)) then
               lo = mid
            else
               hi = mid
            end if
         end do

         mid = 0.5_wp*(lo + hi)
         lip = max(lip_atom + exp(log_b + mreal*mid), &
                   lip_pair + (mreal - 1.0_wp)*exp(-conj*mid))
      end if

      ! Catches NaN as well as a bound that overflowed on its way here
      if (.not. (lip > 0.0_wp)) return
      if (lip > sqrt(huge(0.0_wp))) return

      self%excl_inv_lip = 1.0_wp/lip
   end subroutine cfc_set_exclusion_bound

   !> Radius of a ball around the evaluation point free of surface
   !>
   !> Fade off: `|S|/L` with the cached Lipschitz bound `L` of the untruncated model
   !>
   !> Fade on: `r = |tau(S)|`, `tau(phi) = int_0^phi ds/Lam(s)`
   !>
   !> - `|grad S| <= Lam(S) = L' + Kq e^(S/2) + K0 e^S` at every point of level
   !>   `S`: today's `L` for the weighted terms, the band slope `|w'|` bounded by
   !>   the level alone (`n` band atoms, Holder exponent 1/2 in the partners),
   !>   so `T = tau(S)` is exactly 1-Lipschitz
   !> - closed form per point, constants from [[certificate_setup]]:
   !>   - `S <= 0`: `r = -S/Lam(0)`
   !>   - `S > 0`: `r = tau_lb(min(S, s_cap))` ([[tau_lower]]); `S = +Inf` gives `r_cap`
   !> - at the defaults equals `|S|/L` up to S ~ 10 to ~1e-4 relative for
   !>   n <~ 100, ~1e-3 at n = 1000
   !> - zero (certify nothing) outside the guards of [[certificate_setup]] and
   !>   while the constants are stale ([[cfc_certificate_type]]); raw `|S|/L`
   !>   over-claims inside a band (1.01-1.63x on constructed cases)
   !>
   !> A NaN value (an unreached batch column) or `-Inf` certifies nothing in both
   !> modes, `+Inf` nothing with the fade off
   !>
   !> @param[in] self  LSF instance
   !> @param[in] lsf0  LSF value at the evaluation point
   !> @returns   r     Surface-free radius (zero when uncertified)
   pure function lsf_exclusion_radius(self, lsf0) result(r)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> LSF value at the evaluation point
      real(wp), intent(in) :: lsf0
      !> Surface-free radius
      real(wp) :: r

      r = 0.0_wp
      if (ieee_is_nan(lsf0)) return
      if (lsf0 < -huge(lsf0)) return
      if (.not. certificate_current(self)) return
      if (.not. self%param%fade_enabled) then
         if (lsf0 > huge(lsf0)) return
         r = abs(lsf0)*self%excl_inv_lip
         return
      end if

      if (.not. self%cert%valid) return
      if (lsf0 <= 0.0_wp) then
         r = -lsf0*self%cert%inv_lam0
      else
         r = tau_lower(self%cert, min(lsf0, self%cert%s_cap))
      end if
      r = max(r, 0.0_wp)
   end function lsf_exclusion_radius

   !> Surface-free radius around every unreached point, fade on
   !>
   !> `max(R_min (x_on + delta - xi*), r_cap)`, `xi*` the root of
   !> `F_n(xi) = n e^(-alpha xi) + n(n-1)/2 c (1 + (1 + xi)^2)^m e^(-2 beta xi) = 1`
   !>
   !> - `PD_f <= F_n(min_A x_A)` (`F_n` bounds the `n` atom and `n(n-1)/2` pair
   !>   terms at the smallest reduced distance), so every surface point has an
   !>   atom with `d_A <= R_A xi*`; a
   !>   point with no atom within the band end has `d_A >= R_A (x_on + delta)`
   !>   for every atom, and `d_A` is 1-Lipschitz
   !> - `r_cap`: `tau(+Inf)` of the reached certificate, valid where `PD_f = 0`
   !> - status 1 means no atom within the band end OR a non-finite lift; the
   !>   certificate covers the first, which every status-1 column is at order 0
   !>   (the octree's probe; `w_A` underflowing to zero needs the last ~1e-60 of
   !>   a band); from order 1 a column on a nucleus is status 1 too and is not
   !>   covered
   !> - per-point path: every point with an empty active list
   !> - zero with the fade off (no reach, no such bound), outside the guards of
   !>   [[certificate_setup]] (`F_n(50) > 1` included), or while the constants
   !>   are stale; `r_cap` alone once `xi*` exceeds the band end
   !>
   !> @param[in] self LSF instance
   !> @returns   r    Surface-free radius around an unreached point
   pure function lsf_far_exclusion_radius(self) result(r)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Surface-free radius
      real(wp) :: r

      r = 0.0_wp
      if (.not. self%param%fade_enabled) return
      if (.not. certificate_current(self)) return
      if (.not. self%cert%far_valid) return
      r = self%cert%far_radius
   end function lsf_far_exclusion_radius

   !> Whether the certificate constants belong to the live parameters and atom count
   !>
   !> @param[in] self    LSF instance
   !> @returns   current Constants match `param` and `ncenters`
   pure function certificate_current(self) result(current)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Constants match `param` and `ncenters`
      logical :: current

      current = self%cert%n == self%ncenters .and. &
                (self%cert%fade .eqv. self%param%fade_enabled) .and. &
                all(self%cert%key == [self%param%a1, self%param%a2, self%param%c, &
                                      real(self%param%m, wp), self%param%fade_start, &
                                      self%param%fade_width])
   end function certificate_current

   !> Setup-time constants of the faded exclusion certificates
   !>
   !> The constants of the faded certificate bounds, once per `update`, after
   !> [[cfc_set_exclusion_bound]] (`L`):
   !>
   !>    alpha = -a1, beta = -a2, gam = beta - alpha/2, x_on, delta = [[fade_band]]
   !>    mu_a = max_t S'(t) e^(-alpha delta t)
   !>    mu_b = max_t S'(t) (2 + x_on + delta t)^m e^(-beta delta t)   ([[fade_slope_max]])
   !>    C_h = max(1, (m/gam)^m e^(gam - m)),  c_in = e^(max(0, beta - alpha))
   !>    pre = n/(delta R_min),  e_p = c mu_b e^(-beta x_on)
   !>    K1 = pre e_p c_in,  Kq = pre e_p C_h sqrt(n - 1),  K0 = pre mu_a e^(-alpha x_on)
   !>    L' = L + K1,  Lam(0) = L' + Kq + K0
   !>    s_cap = argmax tau_lb,  r_cap = tau_lb(s_cap)                  ([[certificate_cap]])
   !>    xi*: bisection of F_n(xi) = 1 on [0, 50]
   !>    far radius = max(R_min (x_on + delta - xi*), r_cap)
   !>
   !> Guards (else that certificate is off):
   !>
   !> - both: fade on, valid band, `alpha (x_on + delta) <= 600` (the bounds'
   !>   exponentials underflow beyond), `alpha > 0`, `beta > 0`, `c >= 0`, even
   !>   `m >= 2` and `m` = [[kernel_pair_power]] (the kernel's pair power),
   !>   `m <= beta` (the pair envelope decreases beyond the surface),
   !>   `beta <= alpha`, every radius positive, `n > 0`
   !> - `beta <= alpha` is a conservative choice, not required by the derivation
   !>   (Lemma 3 covers `beta > alpha` through `c_in`); it mirrors [[fade_band]],
   !>   which counts the band in units of `|a2|`
   !> - far: `F_n(50) <= 1`; otherwise only the far block is skipped
   !>   (`far_valid` off, far radius 0), the reached certificate is still built
   !> - reached: also `gam > 0` and a finite `L` (today's conditions)
   !>
   !> @param[inout] self LSF instance after [[cfc_set_exclusion_bound]]
   subroutine certificate_setup(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(inout) :: self

      !> Upper end of the `xi*` bisection and its steps
      real(wp), parameter :: xi_hi = 50.0_wp
      integer, parameter :: xi_steps = 200
      !> Largest `alpha (x_on + delta)` certified; `e^(-600)` is ~1e-261, still
      !> normal in double precision
      real(wp), parameter :: band_exp_max = 600.0_wp
      !> Exponent magnitudes, coupling, pair power and `beta - alpha/2`
      real(wp) :: alpha, beta, c, mreal, gam
      !> Fade band in sigma
      real(wp) :: x_on, delta
      !> Smallest radius, atom count, partner constants and prefactors
      real(wp) :: rmin, nreal, c_h, c_in, pre, e_p
      !> Bisection bracket
      real(wp) :: lo, hi, mid
      integer :: i

      self%cert = cfc_certificate_type()
      self%cert%key = [self%param%a1, self%param%a2, self%param%c, real(self%param%m, wp), &
                       self%param%fade_start, self%param%fade_width]
      self%cert%fade = self%param%fade_enabled
      self%cert%n = self%ncenters

      if (.not. self%param%fade_enabled) return
      if (.not. fade_band_valid(self)) return
      if (self%ncenters <= 0) return
      if (.not. allocated(self%radii)) return
      if (size(self%radii) < self%ncenters) return

      alpha = -self%param%a1
      beta = -self%param%a2
      c = self%param%c
      mreal = real(self%param%m, wp)
      if (.not. all(ieee_is_finite([alpha, beta, c]))) return
      if (.not. (alpha > 0.0_wp .and. beta > 0.0_wp)) return
      if (c < 0.0_wp) return
      if (self%param%m < 2 .or. modulo(self%param%m, 2) /= 0) return
      if (self%param%m /= kernel_pair_power) return
      if (mreal > beta) return
      ! Conservative, as `fade_band` counts in |a2|; Lemma 3 covers beta > alpha via c_in
      if (beta > alpha) return
      if (.not. all(self%radii(1:self%ncenters) > 0.0_wp)) return
      rmin = minval(self%radii(1:self%ncenters))
      nreal = real(self%ncenters, wp)

      call fade_band(self, x_on, delta)
      if (.not. alpha*(x_on + delta) <= band_exp_max) return
      self%cert%x_on = x_on
      self%cert%delta = delta

      ! Far certificate: F_n decreases from F_n(0) >= n >= 1; skipped for F_n(xi_hi) > 1
      if (far_bound(0.0_wp) <= 1.0_wp) then
         self%cert%xi_star = 0.0_wp
         self%cert%far_valid = .true.
      else if (far_bound(xi_hi) <= 1.0_wp) then
         lo = 0.0_wp
         hi = xi_hi
         do i = 1, xi_steps
            mid = 0.5_wp*(lo + hi)
            if (far_bound(mid) > 1.0_wp) then
               lo = mid
            else
               hi = mid
            end if
         end do
         self%cert%xi_star = hi
         self%cert%far_valid = .true.
      end if
      if (self%cert%far_valid) then
         self%cert%far_radius = max(0.0_wp, rmin*(x_on + delta - self%cert%xi_star))
      end if

      ! Reached certificate: today's L (needs gam > 0)
      gam = beta - 0.5_wp*alpha
      if (.not. gam > 0.0_wp) return
      if (.not. self%excl_inv_lip > 0.0_wp) return

      associate (k => self%cert)
         k%lip = 1.0_wp/self%excl_inv_lip
         k%mu_a = fade_slope_max(alpha*delta, 0.0_wp, 0.0_wp, 0)
         k%mu_b = fade_slope_max(beta*delta, 2.0_wp + x_on, delta, self%param%m)
         c_h = max(1.0_wp, (mreal/gam)**self%param%m*exp(gam - mreal))
         c_in = exp(max(0.0_wp, beta - alpha))
         pre = nreal/(delta*rmin)
         e_p = c*k%mu_b*exp(-beta*x_on)
         k%k1 = pre*e_p*c_in
         k%kq = pre*e_p*c_h*sqrt(nreal - 1.0_wp)
         k%k0 = pre*k%mu_a*exp(-alpha*x_on)
         k%lp = k%lip + k%k1
         k%inv_lam0 = 1.0_wp/(k%lp + k%kq + k%k0)
         if (.not. (ieee_is_finite(k%lp) .and. k%inv_lam0 > 0.0_wp)) return
         call certificate_cap(k)
         if (k%far_valid) k%far_radius = max(k%far_radius, k%r_cap)
         k%valid = .true.
      end associate

   contains

      !> `F_n(xi)`, the bound of `PD_f` by the nearest atom's reduced distance
      !>
      !> @param[in] xi Reduced distance, >= 0
      !> @returns   f  `F_n(xi)`
      pure function far_bound(xi) result(f)
         !> Reduced distance
         real(wp), intent(in) :: xi
         !> Bound
         real(wp) :: f

         f = nreal*exp(-alpha*xi) + 0.5_wp*nreal*(nreal - 1.0_wp)*c* &
             (1.0_wp + (1.0_wp + xi)**2)**self%param%m*exp(-2.0_wp*beta*xi)
      end function far_bound

   end subroutine certificate_setup

   !> Whether the configured fade band is finite, positive and has `a2 < 0`
   !>
   !> - the tests of the parameter validation without an error object, for pure
   !>   callers
   !>
   !> @param[in] self  LSF instance
   !> @returns   valid Band usable
   pure function fade_band_valid(self) result(valid)
      !> LSF instance
      class(moist_cavity_drop_lsf_cfc_type), intent(in) :: self
      !> Band usable
      logical :: valid

      associate (a => self%param%fade_start, b => self%param%fade_width, &
                 a2 => self%param%a2)
         valid = ieee_is_finite(a) .and. ieee_is_finite(b) .and. ieee_is_finite(a2)
         if (valid) valid = a > 0.0_wp .and. b > 0.0_wp .and. a2 < 0.0_wp
      end associate
   end function fade_band_valid

   !> Closed-form lower bound `tau_lb(phi)` of `tau(phi)`, `phi >= 0`
   !>
   !> `(phi - 2 ln((L' + Gam e)/(L' + Gam)))/L'`, `e = e^(phi/2)`, `Gam = Kq + K0 e`
   !> (from `Kq e^(s/2) + K0 e^s <= Gam(phi) e^(s/2)` on `[0, phi]`)
   !>
   !> @param[in] k   Certificate constants
   !> @param[in] phi Level, >= 0
   !> @returns   t   Lower bound of tau (Bohr)
   pure function tau_lower(k, phi) result(t)
      !> Certificate constants
      type(cfc_certificate_type), intent(in) :: k
      !> Level
      real(wp), intent(in) :: phi
      !> Lower bound of tau
      real(wp) :: t

      !> `e^(phi/2)` and `Gam(phi)`
      real(wp) :: e, gam

      e = exp(0.5_wp*phi)
      gam = k%kq + k%k0*e
      t = (phi - 2.0_wp*log((k%lp + gam*e)/(k%lp + gam)))/k%lp
   end function tau_lower

   !> Cap level `s_cap` (argmax of `tau_lb` over `phi >= 0`) and `r_cap = tau_lb(s_cap)`
   !>
   !> - tightness only: `tau` increases and `tau_lb <= tau` at every `phi`, so
   !>   `tau_lb(min(S, s_cap))` is a valid radius for any `s_cap >= 0`
   !> - search on `[0, phi_hi]`, `phi_hi` where `Kq e^(phi/2)` or `K0 e^phi`
   !>   reaches 1e8 (`tau` has stopped growing there), at most 1000: 2048-point
   !>   grid, then golden section on the two cells around the grid maximum
   !>
   !> @param[inout] k Certificate constants with `L'`, `Kq`, `K0` set
   pure subroutine certificate_cap(k)
      !> Certificate constants
      type(cfc_certificate_type), intent(inout) :: k

      !> Grid size and golden-section steps
      integer, parameter :: ngrid = 2048, ngold = 80
      !> Slope term that ends the search, and the largest level searched
      real(wp), parameter :: k_end = 1.0e8_wp, phi_max = 1000.0_wp
      !> Inverse golden ratio
      real(wp), parameter :: ginv = 0.6180339887498949_wp
      !> Search range and grid step
      real(wp) :: phi_hi, h
      !> Golden-section bracket and probes
      real(wp) :: lo, hi, x1, x2, f1, f2, t, best, phi_best
      integer :: i

      phi_hi = phi_max
      if (k%kq > 0.0_wp) phi_hi = min(phi_hi, 2.0_wp*log(k_end/k%kq))
      if (k%k0 > 0.0_wp) phi_hi = min(phi_hi, log(k_end/k%k0))
      phi_hi = max(0.0_wp, phi_hi)

      h = phi_hi/real(ngrid, wp)
      phi_best = 0.0_wp
      best = tau_lower(k, 0.0_wp)
      do i = 1, ngrid
         t = tau_lower(k, h*real(i, wp))
         if (t > best) then
            best = t
            phi_best = h*real(i, wp)
         end if
      end do

      lo = max(0.0_wp, phi_best - h)
      hi = min(phi_hi, phi_best + h)
      x1 = hi - ginv*(hi - lo)
      x2 = lo + ginv*(hi - lo)
      f1 = tau_lower(k, x1)
      f2 = tau_lower(k, x2)
      do i = 1, ngold
         if (f1 < f2) then
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + ginv*(hi - lo)
            f2 = tau_lower(k, x2)
         else
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - ginv*(hi - lo)
            f1 = tau_lower(k, x1)
         end if
      end do
      if (max(f1, f2) > best) then
         best = max(f1, f2)
         phi_best = merge(x1, x2, f1 >= f2)
      end if

      k%s_cap = phi_best
      k%r_cap = max(0.0_wp, best)
   end subroutine certificate_cap

   !> Fade slope maximum `max_t S'(t) (p0 + p1 t)^m e^(-q t)` on `[0, 1]`
   !>
   !> - `mu_a`: `q = alpha delta`, `m = 0`; `mu_b`: `q = beta delta`,
   !>   `p0 = 2 + x_on`, `p1 = delta`, `m` the pair power
   !> - `S = S_n`, `n` = [[fade_order]], `S'(t) = (2n+1)!/(n!)^2 t^n (1-t)^n`;
   !>   every factor is log-concave, so the product is unimodal and golden
   !>   section on the bracket of a 1024-cell grid finds the maximum; the result
   !>   is raised by 1e-10 relative against rounding
   !> - grid endpoints skipped (`S' = 0` there), so no integer power of zero
   !>
   !> @param[in] q  Exponential rate over the band
   !> @param[in] p0 Polynomial base at `t = 0`
   !> @param[in] p1 Polynomial base slope
   !> @param[in] m  Polynomial power (0: no polynomial)
   !> @returns   mu Maximum
   pure function fade_slope_max(q, p0, p1, m) result(mu)
      !> Exponential rate
      real(wp), intent(in) :: q
      !> Polynomial base at t = 0
      real(wp), intent(in) :: p0
      !> Polynomial base slope
      real(wp), intent(in) :: p1
      !> Polynomial power
      integer, intent(in) :: m
      !> Maximum
      real(wp) :: mu

      !> Grid size and golden-section steps
      integer, parameter :: ngrid = 1024, ngold = 80
      !> Inverse golden ratio
      real(wp), parameter :: ginv = 0.6180339887498949_wp
      !> Slope prefactor of `S_n'`
      real(wp) :: cs
      !> Golden-section bracket and probes
      real(wp) :: lo, hi, x1, x2, f1, f2, t, best, t_best
      integer :: i

      cs = 1.0_wp
      do i = fade_order + 1, 2*fade_order + 1
         cs = cs*real(i, wp)
      end do
      do i = 1, fade_order
         cs = cs/real(i, wp)
      end do

      ! Interior points only: S' is exactly 0 at both ends, where `0**n` would
      ! raise divide-by-zero in arm64 `__powidf2`
      best = 0.0_wp
      t_best = 0.5_wp
      do i = 1, ngrid - 1
         t = real(i, wp)/real(ngrid, wp)
         f1 = slope_term(t)
         if (f1 > best) then
            best = f1
            t_best = t
         end if
      end do
      lo = max(0.0_wp, t_best - 1.0_wp/real(ngrid, wp))
      hi = min(1.0_wp, t_best + 1.0_wp/real(ngrid, wp))
      x1 = hi - ginv*(hi - lo)
      x2 = lo + ginv*(hi - lo)
      f1 = slope_term(x1)
      f2 = slope_term(x2)
      do i = 1, ngold
         if (f1 < f2) then
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + ginv*(hi - lo)
            f2 = slope_term(x2)
         else
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - ginv*(hi - lo)
            f1 = slope_term(x1)
         end if
      end do
      mu = (1.0_wp + 1.0e-10_wp)*max(best, f1, f2)

   contains

      !> `S'(t) (p0 + p1 t)^m e^(-q t)` at one band coordinate
      !>
      !> @param[in] tt Band coordinate in [0, 1]
      !> @returns   f  Term value
      pure function slope_term(tt) result(f)
         !> Band coordinate
         real(wp), intent(in) :: tt
         !> Term value
         real(wp) :: f

         f = cs*(tt*(1.0_wp - tt))**fade_order*exp(-q*tt)
         if (m > 0) f = f*(p0 + p1*tt)**m
      end function slope_term

   end function fade_slope_max

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
      if (allocated(self%bat_cand)) deallocate (self%bat_cand)
      if (allocated(self%bat_center)) deallocate (self%bat_center)
      if (allocated(self%bat_radius)) deallocate (self%bat_radius)
      if (allocated(self%bat_bsq)) deallocate (self%bat_bsq)
      if (allocated(self%strip_geo)) deallocate (self%strip_geo)
      if (allocated(self%bat_sv)) deallocate (self%bat_sv)
      if (allocated(self%bat_reach)) deallocate (self%bat_reach)
      if (allocated(self%vjp_jet3)) deallocate (self%vjp_jet3)
      if (allocated(self%vjp_gacc)) deallocate (self%vjp_gacc)
      if (allocated(self%vjp_res)) deallocate (self%vjp_res)
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
