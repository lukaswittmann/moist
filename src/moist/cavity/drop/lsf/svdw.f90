!> Smooth van der Waals (SvdW) level set function
!>
!> Orchestration between the code-generated [[moist_cavity_drop_lsf_svdw_kernel]] and the LSF
!> contract of [[moist_cavity_drop_lsf_base]]
!>
!> Division of labour
!> ------------------
!> - kernel symbolic in the power sums P_m = sum_A u_A**m, m = 1, 2, 3, 3/2, with
!>   u_A = exp(-(k/3)(||r - R_A|| - R_a)); Z, S = -(1/k) log Z and every derivative in
!>   O(n_active)
!> - `prepare` / `prepare_subset`: base screening gate, then cache only per active atom its
!>   user-space id, `r - R_A`, radius and `||r - R_A||`, plus the power-sum tensors (packed for
!>   the lift, dense for the per-atom families)
!> - buffers sized once per `update`, never reallocated per point
!> - accessors rebuild the tensors they need from the kernel; nothing derivative-shaped stored
!>   between points, no SSD derivative cache
!>
!> Buffers, never allocations
!> --------------------------
!> - no `allocatable, intent(out)` results: caller buffers, nuclear extent from
!>   [[active_count]], first `active_count()` slots written, rest untouched
!> - keeps the OpenMP projection path free of per-call heap traffic
!> - symmetric spatial groups packed (rank 2: 6, rank 3: 10, rank 4: 15) as
!>   [[moist_math_packed]]; nuclear and radius indices dense
!> - nuclear indices are active-list slots, slot `i` the atom `active_atom(i)`; screened-away
!>   atoms have no slot, so no zeroing over the whole molecule
!>
!> Derivative-order contract
!> -------------------------
!> - `prepare` provisions only the power sums, up to `max_deriv`
!> - each accessor asks [[require_deriv]] for the power-sum order its kernel routine reads, not
!>   the order of the tensor it returns, and rebuilds that tensor itself
!> - e.g. `f4_rrr_rA` needs `ps3` and builds order-4 per-atom tensors on the spot, so it runs
!>   at `set_max_deriv(3)`
!>
!> Strip kernel
!> ------------
!> - stage 1 `svdw_strip_powersums_o{L}`: power sums of a strip over one shared atom list,
!>   inline `exp`, the gate's squared compare as mask; stage 2 `svdw_strip_lift_o{L}`: packed
!>   jet of S
!> - `prepare`: stage 1 as a strip of one at `max_deriv`, lifted once (`np = 1`); `f0`,
!>   `f012_r`, `f3_rrr`, `f4_rrrr` and the value half of `normalized_f01_rA` read that jet
!> - [[lsf_evaluate_batch]]: both stages per strip of `svdw_strip_width(order)` points; `order`
!>   alone governs, `max_deriv` only the per-point cache
!> - [[lsf_vjp_f1_rA_batch]]: stage 1 at order 2 with geometry cache, order-2 lift,
!>   `svdw_strip_lift_vjp`, `svdw_strip_vjp`, per strip of `svdw_strip_width(2)` points
!> - per-atom families (nuclear, radius, pair, tangent, hvp, vjp): dense unpack of the same
!>   sums and stage 1's geometry cache (`x`, `1/x`, `u`); one exp per atom and point, fade band
!>   of [[fade_band]] in both
!>
!> Fade and reach
!> --------------
!> - fade is the model: default on, `fade_start = 21`, `fade_width = 3` e-folds
!> - `u_A -> w(d_A) u_A`, `w = S_4(1 - t)` over `X_on = 3a/k`, `Delta = 3b/k` ([[fade_band]])
!> - [[lsf_screening_offset]] returns the band end `X_on + Delta`: gate mask `gr2 <= bsq` and
!>   fade agree, masked atoms add exact zeros, `screening_threshold` changes nothing
!> - exclusion certificates from [[certificate_setup]]: `exclusion_radius = |tau(S)|`,
!>   `far_exclusion_radius` for points with no atom in reach (needs `blend_1b > 0`)
!> - fade off: untruncated model, hard screening at the threshold
!>
!> Handled here, not in the kernel
!> -------------------------------
!> - point on a nucleus (kernel divides by `||r - R_A||`):
!>   - stage 1 counts zero-distance pairs (`n_zero`); `prepare` rebuilds the sums in
!>     [[powersums_guarded]], `evaluate_batch` redoes that strip point by point
!>   - [[atom_tensors]] and [[atom_tangent_tensors]] intercept it: at `x == 0` the kernel's own
!>     `u_A` as value, every derivative order zero (`n = (r - R_A)/x` taken as zero)
!> - no atom in reach:
!>   - per point: `n_active == 0` or a non-finite lift through `max_deriv` with `S > 0` or NaN
!>     (`Z` underflow) marks the point unreached (`unreached()`, no `error`); later
!>     accessors abort (`require_reached`), except `active_count` / `active_atom`
!>   - batch: `p1 == 0` after stage 1 or a non-finite lift with `S > 0` or NaN marks the column
!>     unreached (`unreached(column)`); `S` kept where finite, else `+Inf`, derivatives NaN
!>   - a non-finite lift with `S <= 0` is an `error`, never unreached
!> - blending weights used as given, no low-`n_active` guard (a shared batch list leaves the
!>   atom count unknown); the blocks vanish there up to the eps-relative residue of
!>   `q**2 - p3` and `p1**3 - 3 p1 p2 + 2 p3`
module moist_cavity_drop_lsf_svdw
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
                                            ieee_quiet_nan, ieee_positive_inf
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_base_update, lsf_candidate_space_sorted, &
                                         lsf_cand_bound_sq, lsf_batch_require_outputs, &
                                         lsf_sorted_candidates, lsf_screen_records, &
                                         lsf_vjp_batch_shape_message
   use moist_math_packed, only: pack_sym_matrix, pack_sym_tensor3, pack_sym_tensor4, &
                                    unpack_sym_matrix
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_svdw_kernel, only: nkind, &
                                                svdw_atom_eval, svdw_atom_tangent_eval, &
                                                svdw_nuclear_eval, &
                                                svdw_radius_eval, svdw_radius_hvp_eval, &
                                                svdw_atom_radius_tangent_eval, &
                                                svdw_atom_rad_eval, svdw_atom_radrad_eval, &
                                                svdw_atom_tangent_rad_eval, &
                                                svdw_pair_eval, svdw_pair_diag_eval, &
                                                svdw_radpair_eval, svdw_radpair_diag_eval, &
                                                svdw_nucrad_eval, svdw_nucrad_diag_eval, &
                                                svdw_tangent_eval, svdw_hvp_eval, &
                                                svdw_vjp_eval, svdw_radius_vjp_eval, &
                                                svdw_normalized_eval, &
                                                svdw_strip_width, svdw_npack, &
                                                svdw_strip_powersums_o0, svdw_strip_lift_o0, &
                                                svdw_strip_powersums_o1, svdw_strip_lift_o1, &
                                                svdw_strip_powersums_o2, svdw_strip_lift_o2, &
                                                svdw_strip_powersums_o3, svdw_strip_lift_o3, &
                                                svdw_strip_powersums_o4, svdw_strip_lift_o4, &
                                                svdw_strip_lift_vjp, svdw_strip_vjp
   implicit none(type, external)
   private

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Packed components of the rank-2, -3 and -4 spatial tensors
   integer, parameter :: np2 = svdw_npack(2), np3 = svdw_npack(3), np4 = svdw_npack(4)
   !> Widest strip of any order
   integer, parameter :: max_strip = maxval(svdw_strip_width)
   !> Strip width of the nuclear vjp pipeline, the order-2 width
   integer, parameter :: vjp_strip = svdw_strip_width(2)

   !> Packed slot of a dense rank-2 index, the kernel's and `moist_math_packed`'s order
   integer, parameter :: slot2(ndim, ndim) = reshape([1, 2, 4, 2, 3, 5, 4, 5, 6], [ndim, ndim])
   !> Packed slot of a dense rank-3 index
   integer, parameter :: slot3(ndim, ndim, ndim) = reshape([ &
                         1, 2, 5, 2, 3, 6, 5, 6, 8, &
                         2, 3, 6, 3, 4, 7, 6, 7, 9, &
                         5, 6, 8, 6, 7, 9, 8, 9, 10], [ndim, ndim, ndim])

   !> Fade band at infinity: the untruncated model, every atom on the plain body
   !>
   !> - `fade_start + fade_width` must stay finite, hence a finite width
   !> - fade on: the band of [[fade_band]], from `param`
   real(wp), parameter :: fade_off_start = huge(1.0_wp)
   !> Fade band width while the fade is off
   real(wp), parameter :: fade_off_width = 1.0_wp

   !> Smoothstep order of the kernel's fade `w = S_n(1 - t)` (generator constant)
   integer, parameter :: fade_order = 4

   !> Setup-time constants of the faded exclusion certificates
   !>
   !> - built on the band bound `||grad S|| <= 1 + K(S)` of the faded model
   !>   ([[lsf_svdw_exclusion_radius]]) and the far bound on `Z` behind
   !>   [[lsf_svdw_far_exclusion_radius]]; constants in [[certificate_setup]];
   !>   `n = ncenters`
   !> - refreshed by `update`; `key` and `n` record the parameters they were
   !>   built from, and a mismatch with the live `param` / `ncenters` makes both
   !>   radii return zero (safe), so `param` mutated after `update` certifies
   !>   nothing until the next `update`
   type :: svdw_certificate_type
      !> Reached certificate available: fade on, `s >= 0`, `s1 > 0`, `s3 < 24 s1`
      logical :: valid = .false.
      !> Far certificate available: fade on, `s >= 0`, `s1 > 0`, `n > 0`
      logical :: far_valid = .false.
      !> Parameters the constants belong to: k, s1, s2, s3, a, b
      real(wp) :: key(6) = 0.0_wp
      !> Atom count the constants belong to
      integer :: n = -1
      !> Blending sharpness k
      real(wp) :: k = 0.0_wp
      !> Fade onset in e-folds `a`
      real(wp) :: a = 0.0_wp
      !> `eps = exp(-a)`, the largest `u` beyond the onset
      real(wp) :: eps = 0.0_wp
      !> Slope bound `K(eta) = A eta**3 + B eta**(3/2) + G eta`, coefficient A
      real(wp) :: coef_a = 0.0_wp
      !> Coefficient B of `K`
      real(wp) :: coef_b = 0.0_wp
      !> Coefficient G of `K`
      real(wp) :: coef_g = 0.0_wp
      !> `K` at `S <= 0`, `K0 = A eps**3 + B eps**(3/2) + G eps`
      real(wp) :: k0 = 0.0_wp
      !> Level where the closed-form radius peaks (Bohr)
      real(wp) :: s_cap = 0.0_wp
      !> Radius at `s_cap` (Bohr)
      real(wp) :: r_cap = 0.0_wp
      !> `Zmax`, `Z <= Zmax e1**3` on `g >= 0`
      real(wp) :: zmax = 0.0_wp
      !> Nearest-atom surface distance bound of any surface point (Bohr)
      real(wp) :: dstar = 0.0_wp
      !> `X_on + Delta - dstar`, clipped at zero (Bohr)
      real(wp) :: far_radius = 0.0_wp
   end type svdw_certificate_type

   !> Smooth van der Waals LSF
   !>
   !> Concrete level set function: takes its blending parameters directly through
   !> [[new]], caches the screened per-atom geometry and the power sums via
   !> [[prepare]], and rebuilds every derivative tensor on demand from the
   !> generated kernel
   !>
   !> - inherits the common atom-LSF state (ncenters, mol, radii, screening
   !>   caches) from [[moist_cavity_drop_lsf_type]]
   type, extends(moist_cavity_drop_lsf_type) :: moist_cavity_drop_lsf_svdw_type

      !> SvdW blending parameters (k, one-, two- and three-body weights)
      type(moist_cavity_drop_lsf_svdw_param_type) :: param

      !> Highest spatial-derivative order the power sums are accumulated to
      integer :: max_deriv = 2

      !> Exclusion-certificate constants of the faded model, refreshed by `update`
      type(svdw_certificate_type) :: cert

      !* --------------------------- Per-point active list ---------------------------- *!

      !? Belongs in the base LSF type; the isodensity LSFs still need them added

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
      !> Distance `||r - R_A||` of active slot i [ncenters]
      !>
      !> - stage 1's geometry cache, the same `x` its sums saw
      !> - read by the on-nucleus guard and by every per-atom kernel
      real(wp), allocatable :: act_x(:)
      !> Reciprocal distance of active slot i, stage 1's cache [ncenters]
      !>
      !> - `Inf` on a nucleus, where only order 0 of an atom is evaluated
      real(wp), allocatable :: act_inv_x(:)
      !> Unfaded screening factor `u_A` of active slot i, stage 1's cache [ncenters]
      !>
      !> - the inline exp of the sums, so sums and per-atom tensors share it
      real(wp), allocatable :: act_u(:)
      !> Center of active slot i, gathered for stage 1 [ndim, ncenters]
      real(wp), allocatable :: act_center(:, :)
      !> Squared screening reach of active slot i, gathered for stage 1 [ncenters]
      real(wp), allocatable :: act_bsq(:)

      !* ---------------------------- Batch candidate list ---------------------------- *!

      !> Shared atom list of the current batch strip, ascending candidate ids [ncenters]
      integer, allocatable :: bat_cand(:)
      !> Centers of the shared list [ndim, ncenters]
      real(wp), allocatable :: bat_center(:, :)
      !> Radii of the shared list [ncenters]
      real(wp), allocatable :: bat_radius(:)
      !> Squared screening reach of the shared list [ncenters]
      real(wp), allocatable :: bat_bsq(:)

      !* ------------------------- Batched vjp strip scratch -------------------------- *!
      !> Grown on demand by [[lsf_vjp_f1_rA_batch]]; the kernel sees them as
      !> `(np, nc)` and `(3, nc, np)` by sequence association

      !> Distance of stage 1's geometry cache [vjp_strip*ncenters]
      real(wp), allocatable :: vjp_x(:)
      !> Reciprocal distance of stage 1's geometry cache [vjp_strip*ncenters]
      real(wp), allocatable :: vjp_ix(:)
      !> Unfaded screening factor of stage 1's geometry cache [vjp_strip*ncenters]
      real(wp), allocatable :: vjp_u(:)
      !> Strip result when not written into the tile directly [3*vjp_strip*ncenters]
      real(wp), allocatable :: vjp_res(:)
      !> Position in `candidate_indices` by candidate id [ncenters]; defined only at
      !> the ids of the current call's list
      integer, allocatable :: vjp_pos(:)
      !> Position in `candidate_indices` of shared-list atom `j` [ncenters]
      integer, allocatable :: vjp_apos(:)
      !> Tile position by user-space atom id for the on-nucleus fallback [ncenters];
      !> zero between calls, set and cleared around each [[vjp_point]]
      integer, allocatable :: vjp_upos(:)

      !* --------------------------------- Power sums --------------------------------- *!
      !> Orders above `max_deriv` are stale after `prepare`, the accessors guarding
      !> on [[require_deriv]] rather than on the contents
      !>
      !> - `pk*`: stage-1 layout of a strip of one, packed, read by the lift
      !> - `ps*`: dense unpack through order 3, read by the per-atom families

      !> Order-0 sums, strip layout `(1, kind)` (1 = p1, 2 = p2, 3 = p3, 4 = q)
      real(wp) :: pk0(1, nkind) = 0.0_wp
      !> Order-1 sums, strip layout
      real(wp) :: pk1(1, ndim, nkind) = 0.0_wp
      !> Order-2 sums, strip layout, packed
      real(wp) :: pk2(1, np2, nkind) = 0.0_wp
      !> Order-3 sums, strip layout, packed
      real(wp) :: pk3(1, np3, nkind) = 0.0_wp
      !> Order-4 sums, strip layout, packed
      real(wp) :: pk4(1, np4, nkind) = 0.0_wp
      !> Order-0 power sums, last index the kind
      real(wp) :: ps0(nkind) = 0.0_wp
      !> Order-1 power-sum tensors
      real(wp) :: ps1(ndim, nkind) = 0.0_wp
      !> Order-2 power-sum tensors
      real(wp) :: ps2(ndim, ndim, nkind) = 0.0_wp
      !> Order-3 power-sum tensors
      real(wp) :: ps3(ndim, ndim, ndim, nkind) = 0.0_wp

      !* --------------------------------- Lifted jet --------------------------------- *!
      !> The lift of the per-point sums at `max_deriv`, once per `prepare`; the
      !> spatial accessors read it, orders above `max_deriv` are stale

      !> Value, strip layout `(1)`
      real(wp) :: jet0(1) = 0.0_wp
      !> Gradient
      real(wp) :: jet1(ndim, 1) = 0.0_wp
      !> Packed Hessian
      real(wp) :: jet2(np2, 1) = 0.0_wp
      !> Packed third derivative
      real(wp) :: jet3(np3, 1) = 0.0_wp
      !> Packed fourth derivative
      real(wp) :: jet4(np4, 1) = 0.0_wp
   contains
      !> Constructor: configure blending parameters and declare the candidate space
      procedure, public :: new => lsf_new
      !> Bind molecular geometry and resize the per-atom caches
      procedure, public :: update => lsf_update
      !> Point preparation: screening plus power-sum accumulation
      procedure, public :: prepare => lsf_prepare
      !> Point preparation with a caller-provided candidate atom list
      procedure, public :: prepare_subset => lsf_prepare_subset
      !> Configure the highest power-sum order `prepare` accumulates
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
      !> Batched spatial jet through the strip kernel
      procedure, public :: evaluate_batch => lsf_evaluate_batch
      !> Batched jet-contracted nuclear gradient through the strip kernel
      procedure, public :: vjp_f1_rA_batch => lsf_vjp_f1_rA_batch
      !> Mixed fourth derivative d^4S / dr^3 dR_A
      procedure, public :: f4_rrr_rA => lsf_f4_rrr_rA
      !> Radius derivative d^3S / dr^2 dR_a
      procedure, public :: f3_rr_rad => lsf_f3_rr_rad
      !> Pure nuclear Hessian d^2S / dR_A dR_B
      procedure, public :: f2_rArB => lsf_f2_rArB
      !> Mixed third derivative d^3S / dr dR_A dR_B
      procedure, public :: f3_r_rArB => lsf_f3_r_rArB
      !> Mixed fourth derivative d^4S / dr^2 dR_A dR_B
      procedure, public :: f4_rr_rArB => lsf_f4_rr_rArB
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
      !> Normalized level set f0 / ||f1_r|| and its nuclear derivative
      procedure, public :: normalized_f01_rA => lsf_normalized_f01_rA
      !> SvdW splits into atomic shares
      procedure, public :: has_atom_partition => lsf_has_atom_partition
      !> Add the atomic shares of Z at the prepared point
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
      !> Jet-contracted nuclear gradient (reverse mode)
      procedure, public :: vjp_f1_rA => lsf_vjp_f1_rA
      !> Jet-contracted radius gradient (reverse mode)
      procedure, public :: vjp_f1_rad => lsf_vjp_f1_rad
      !> Atom reach: the fade band end, or the threshold offset with the fade off
      procedure, public :: screening_offset => lsf_screening_offset
      !> Surface-free radius from the value: `|S|`, or `|tau(S)|` with the fade on
      procedure, public :: exclusion_radius => lsf_svdw_exclusion_radius
      !> Surface-free radius around every unreached point (fade on)
      procedure, public :: far_exclusion_radius => lsf_svdw_far_exclusion_radius
      !> Exclude boxes using the positive cone of atomic surface normals
      procedure, public :: stationary_box_excluded => lsf_svdw_stationary_box_excluded
      !> Finalizer
      final :: finalize_lsf_svdw
   end type moist_cavity_drop_lsf_svdw_type

   public :: moist_cavity_drop_lsf_svdw_type

contains

   !> Construct from parameter values; omission uses compiled defaults
   !>
   !> @param[inout] self Object to initialize
   !> @param[in] param Configuration copied by value
   subroutine lsf_new(self, param)
      !> LSF to initialize
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Configuration; omitted means compiled defaults
      type(moist_cavity_drop_lsf_svdw_param_type), intent(in), optional :: param

      ! The screen loop indexes the base's candidate-space geometry mirror
      ! directly, so candidate ids must arrive spatially sorted
      self%candidate_space = lsf_candidate_space_sorted
      self%name = "SvdW"

      call self%param%init_defaults()
      if (present(param)) self%param = param
   end subroutine lsf_new

   !> Bind molecular geometry and resize the per-atom caches
   !>
   !> The base handles the spatial sort, the candidate-space geometry mirror and
   !> the screening bounds; this override only sizes the per-point buffers, once,
   !> to the molecule's atom count
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    mol    Molecular structure
   !> @param[in]    radii  Per-atom radii (size mol%nat)
   subroutine lsf_update(self, mol, radii)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
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
      if (allocated(self%act_x)) deallocate (self%act_x)
      if (allocated(self%act_inv_x)) deallocate (self%act_inv_x)
      if (allocated(self%act_u)) deallocate (self%act_u)
      if (allocated(self%act_center)) deallocate (self%act_center)
      if (allocated(self%act_bsq)) deallocate (self%act_bsq)
      if (allocated(self%bat_cand)) deallocate (self%bat_cand)
      if (allocated(self%bat_center)) deallocate (self%bat_center)
      if (allocated(self%bat_radius)) deallocate (self%bat_radius)
      if (allocated(self%bat_bsq)) deallocate (self%bat_bsq)

      allocate (self%active_cand(n_alloc))
      allocate (self%act_atom(n_alloc))
      allocate (self%act_d(ndim, n_alloc))
      allocate (self%act_radius(n_alloc))
      allocate (self%act_x(n_alloc))
      allocate (self%act_inv_x(n_alloc))
      allocate (self%act_u(n_alloc))
      allocate (self%act_center(ndim, n_alloc))
      allocate (self%act_bsq(n_alloc))
      allocate (self%bat_cand(n_alloc))
      allocate (self%bat_center(ndim, n_alloc))
      allocate (self%bat_radius(n_alloc))
      allocate (self%bat_bsq(n_alloc))

      self%n_active = 0
      self%prepared_deriv = -1
      self%prepared_unreached = .false.
      call certificate_setup(self)
   end subroutine lsf_update

   !> Configure the highest power-sum order `prepare` accumulates
   !>
   !> Nothing is allocated here: the power sums are fixed-size components and the
   !> per-atom buffers are sized by [[update]]
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    n    Requested max derivative order (0..4)
   subroutine lsf_set_max_deriv(self, n)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Requested max derivative order
      integer, intent(in) :: n

      self%max_deriv = min(4, max(0, n))
   end subroutine lsf_set_max_deriv

   !> Screen every atom at the evaluation point and refresh the caches
   !>
   !> - unreached ([[lsf_svdw_screen]]): no `error`, `unreached()` reports it; the
   !>   cache stays readable-but-unreached, accessors abort
   !> - a jet not finite with `S <= 0`: `error`, the cache invalidated
   !>   ([[report_singular]])
   !> - `error stop` before `update`
   !>
   !> @param[inout] self   LSF instance
   !> @param[in]    point  Evaluation point (3,)
   !> @param[out]   error  Invalid parameters or fade band changed since `update`
   !>                      ([[check_fade]]), cache cleared; or a jet not finite inside
   subroutine lsf_prepare(self, point, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      !> Jet not finite inside
      logical :: singular

      if (.not. allocated(self%cand_screen)) error stop "moist DROP LSF: prepare before update"
      call self%param%validate(error)
      if (.not. allocated(error) .and. self%param%fade_enabled) call check_fade(self, error)
      if (allocated(error)) then
         self%n_active = 0
         self%prepared_deriv = -1
         self%prepared_unreached = .false.
         return
      end if
      call lsf_svdw_screen(self, point, self%full_scan_cand, singular)
      if (singular) call report_singular(self, error)
   end subroutine lsf_prepare

   !> Screen a caller-provided candidate list and refresh the caches
   !>
   !> - unreached, a jet not finite inside and before `update` as in [[lsf_prepare]]
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point (3,)
   !> @param[in]    candidate_indices Atom ids in the base's sorted candidate space
   !> @param[out]   error             Invalid parameters or fade band changed since
   !>                                 `update` ([[check_fade]]), cache cleared; or a
   !>                                 jet not finite inside
   subroutine lsf_prepare_subset(self, point, candidate_indices, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error

      !> Jet not finite inside
      logical :: singular

      if (.not. allocated(self%cand_screen)) error stop "moist DROP LSF: prepare before update"
      call self%param%validate(error)
      if (.not. allocated(error) .and. self%param%fade_enabled) call check_fade(self, error)
      if (allocated(error)) then
         self%n_active = 0
         self%prepared_deriv = -1
         self%prepared_unreached = .false.
         return
      end if
      call lsf_svdw_screen(self, point, candidate_indices, singular)
      if (singular) call report_singular(self, error)
   end subroutine lsf_prepare_subset

   !> Refuse a prepared point whose jet is not finite with `S <= 0`
   !>
   !> - not reached in practice: a point on a nucleus takes the guarded sums
   !>   ([[powersums_guarded]]), which stay finite
   !> - the per-point cache invalidated (`prepared_deriv = -1`): accessors abort
   !>
   !> @param[inout] self  LSF instance after [[lsf_svdw_screen]]
   !> @param[out]   error Refusal
   subroutine report_singular(self, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Refusal
      type(error_type), allocatable, intent(out) :: error

      self%prepared_deriv = -1
      call fatal_error(error, "prepare: level-set jet not finite with atoms in reach and "// &
                       "S <= 0 ("//self%name//")")
   end subroutine report_singular

   !> Run the base screening gate and rebuild the per-point caches
   !>
   !> - pass one: the base's allocation-free reject test
   !> - pass two: per survivor its id, radius, center, gate reach and
   !>   displacement
   !> - pass three: stage 1 of the strip kernel as a strip of one at
   !>   `max_deriv`, the same statements as [[lsf_evaluate_batch]]; packed sums
   !>   kept for the lift, dense unpack through order 3 for the per-atom families
   !> - a survivor at zero distance (stage 1's `n_zero`) sends the point to the
   !>   guarded loop of [[powersums_guarded]]; taken per point, never in
   !>   production
   !> - pass four: the lift at `max_deriv` into the cached jet (`jet*`), read by
   !>   the spatial accessors
   !> - unreached (`prepared_unreached`): no survivor, or any cached output
   !>   through `max_deriv` non-finite without `S <= 0`, the per-point twin of
   !>   [[poison_unreached]]; `Z` not positive (underflow, or `blend_1b = 0`
   !>   with one or two atoms) fails the value, a small positive `Z` only the
   !>   higher orders (powers of `1/Z`: normal `Z` below ~1e-77 at order 4, near
   !>   the fade band end with the defaults)
   !> - singular (`singular`): a cached output non-finite with `S <= 0`; not
   !>   unreached (a NaN `S`, from `Z <= 0` by round-off, is unreached)
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    point             Evaluation point [ndim]
   !> @param[in]    candidate_indices Candidate ids in the sorted candidate space
   !> @param[out]   singular          Jet not finite with `S <= 0`
   subroutine lsf_svdw_screen(self, point, candidate_indices, singular)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Candidate ids in the sorted candidate space
      integer, intent(in) :: candidate_indices(:)
      !> Jet not finite with `S <= 0`
      logical, intent(out) :: singular

      !> Loop counter, candidate id, survivor count, derivative order and
      !> zero-distance pair count
      integer :: i, c, n, md, n_zero

      singular = .false.
      self%n_active = 0
      self%prepared_deriv = self%max_deriv
      self%prepared_unreached = .true.
      if (.not. allocated(self%cand_screen)) error stop "moist DROP LSF: prepare before update"

      call lsf_screen_records(self%cand_screen, point, candidate_indices, self%active_cand, n)
      self%n_active = n

      md = self%max_deriv
      if (n == 0) then
         call clear_powersums(self, md)
         return
      end if

      do i = 1, n
         c = self%active_cand(i)
         self%act_atom(i) = self%cand_to_user(c)
         self%act_radius(i) = self%cand_radii(c)
         self%act_center(1, i) = self%cand_screen(1, c)
         self%act_center(2, i) = self%cand_screen(2, c)
         self%act_center(3, i) = self%cand_screen(3, c)
         self%act_bsq(i) = self%cand_screen(lsf_cand_bound_sq, c)
         self%act_d(1, i) = point(1) - self%cand_screen(1, c)
         self%act_d(2, i) = point(2) - self%cand_screen(2, c)
         self%act_d(3, i) = point(3) - self%cand_screen(3, c)
      end do

      call point_powersums(self, point, md, n_zero)
      if (n_zero > 0) then
         call powersums_guarded(self, md)
      else
         call unpack_powersums(self, md)
      end if
      call point_lift(self, md)
      if (jet_finite(md, self%jet0(1), self%jet1(:, 1), self%jet2(:, 1), self%jet3(:, 1), &
                     self%jet4(:, 1))) then
         self%prepared_unreached = .false.
      else
         self%prepared_unreached = .not. self%jet0(1) <= 0.0_wp
         singular = .not. self%prepared_unreached
      end if
   end subroutine lsf_svdw_screen

   !> Stage 1 of the strip kernel for the prepared point, a strip of one
   !>
   !> @param[inout] self      LSF instance with the active gather filled
   !> @param[in]    point     Evaluation point [ndim]
   !> @param[in]    md        Highest spatial order to accumulate (0..4)
   !> @param[out]   n_zero    Active atoms at zero distance
   subroutine point_powersums(self, point, md, n_zero)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Highest spatial order
      integer, intent(in) :: md
      !> Active atoms at zero distance
      integer, intent(out) :: n_zero

      !> Active-atom count
      integer :: n
      !> Accumulate the kind-4 sums
      logical :: iq
      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      n = self%n_active
      iq = self%param%blend_2b /= 0.0_wp
      call fade_band(self, x_on, delta)
      select case (md)
      case (0)
         call svdw_strip_powersums_o0(1, point, n, self%act_center, self%act_radius, &
                                      self%act_bsq, self%param%blend_k, iq, x_on, &
                                      delta, self%pk0, n_zero, cache_x=self%act_x, &
                                      cache_inv_x=self%act_inv_x, cache_u=self%act_u)
      case (1)
         call svdw_strip_powersums_o1(1, point, n, self%act_center, self%act_radius, &
                                      self%act_bsq, self%param%blend_k, iq, x_on, &
                                      delta, self%pk0, self%pk1, n_zero, &
                                      cache_x=self%act_x, cache_inv_x=self%act_inv_x, &
                                      cache_u=self%act_u)
      case (2)
         call svdw_strip_powersums_o2(1, point, n, self%act_center, self%act_radius, &
                                      self%act_bsq, self%param%blend_k, iq, x_on, &
                                      delta, self%pk0, self%pk1, self%pk2, n_zero, &
                                      cache_x=self%act_x, cache_inv_x=self%act_inv_x, &
                                      cache_u=self%act_u)
      case (3)
         call svdw_strip_powersums_o3(1, point, n, self%act_center, self%act_radius, &
                                      self%act_bsq, self%param%blend_k, iq, x_on, &
                                      delta, self%pk0, self%pk1, self%pk2, self%pk3, &
                                      n_zero, cache_x=self%act_x, &
                                      cache_inv_x=self%act_inv_x, cache_u=self%act_u)
      case default
         call svdw_strip_powersums_o4(1, point, n, self%act_center, self%act_radius, &
                                      self%act_bsq, self%param%blend_k, iq, x_on, &
                                      delta, self%pk0, self%pk1, self%pk2, self%pk3, &
                                      self%pk4, n_zero, cache_x=self%act_x, &
                                      cache_inv_x=self%act_inv_x, cache_u=self%act_u)
      end select
   end subroutine point_powersums

   !> Dense copy of the packed per-point sums through `min(md, 3)`
   !>
   !> @param[inout] self LSF instance with `pk*` filled
   !> @param[in]    md   Highest spatial order accumulated
   pure subroutine unpack_powersums(self, md)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Highest spatial order accumulated
      integer, intent(in) :: md

      !> Kind and spatial indices
      integer :: k, i, j, l

      do k = 1, nkind
         self%ps0(k) = self%pk0(1, k)
         if (md >= 1) self%ps1(:, k) = self%pk1(1, :, k)
         if (md >= 2) then
            do j = 1, ndim
               do i = 1, ndim
                  self%ps2(i, j, k) = self%pk2(1, slot2(i, j), k)
               end do
            end do
         end if
         if (md >= 3) then
            do l = 1, ndim
               do j = 1, ndim
                  do i = 1, ndim
                     self%ps3(i, j, l, k) = self%pk3(1, slot3(i, j, l), k)
                  end do
               end do
            end do
         end if
      end do
   end subroutine unpack_powersums

   !> Zero both copies of the per-point sums through `md`
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    md   Highest spatial order
   pure subroutine clear_powersums(self, md)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Highest spatial order
      integer, intent(in) :: md

      self%pk0 = 0.0_wp
      self%ps0 = 0.0_wp
      if (md >= 1) self%pk1 = 0.0_wp
      if (md >= 1) self%ps1 = 0.0_wp
      if (md >= 2) self%pk2 = 0.0_wp
      if (md >= 2) self%ps2 = 0.0_wp
      if (md >= 3) self%pk3 = 0.0_wp
      if (md >= 3) self%ps3 = 0.0_wp
      if (md >= 4) self%pk4 = 0.0_wp
   end subroutine clear_powersums

   !> Accumulate the power sums with the on-nucleus guard in the loop
   !>
   !> Only reached when some active atom sits exactly at the evaluation point;
   !> see [[atom_tensors]] for the convention applied there
   !>
   !> - accumulates the dense sums from [[atom_tensors]] (stage 1's cached `u`)
   !>   and packs them into the strip layout the lift reads
   !> - always accumulates the `q` kind, even when the two-body weight is zero
   !>   and stage 1 would leave those slots at zero; every consumer reaches that
   !>   kind only through the two-body weight, so both are correct
   !>
   !> @param[inout] self      LSF instance with the per-atom geometry filled
   !> @param[in]    max_deriv Highest spatial-derivative order to accumulate
   subroutine powersums_guarded(self, max_deriv)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: max_deriv

      !> Per-atom kind tensors of one active atom
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Dense order-4 sums, packed on exit
      real(wp) :: ps4(ndim, ndim, ndim, ndim, nkind)
      !> Active-list index and kind
      integer :: ia, k

      call clear_powersums(self, max_deriv)
      if (max_deriv >= 4) ps4 = 0.0_wp
      do ia = 1, self%n_active
         call atom_tensors(self, ia, max_deriv, at0, at1, at2, at3, at4)
         self%ps0 = self%ps0 + at0
         if (max_deriv >= 1) self%ps1 = self%ps1 + at1
         if (max_deriv >= 2) self%ps2 = self%ps2 + at2
         if (max_deriv >= 3) self%ps3 = self%ps3 + at3
         if (max_deriv >= 4) ps4 = ps4 + at4
      end do

      do k = 1, nkind
         self%pk0(1, k) = self%ps0(k)
         if (max_deriv >= 1) self%pk1(1, :, k) = self%ps1(:, k)
         if (max_deriv >= 2) call pack_sym_matrix(self%ps2(:, :, k), self%pk2(1, :, k), ndim)
         if (max_deriv >= 3) call pack_sym_tensor3(self%ps3(:, :, :, k), self%pk3(1, :, k), ndim)
         if (max_deriv >= 4) call pack_sym_tensor4(ps4(:, :, :, :, k), self%pk4(1, :, k), ndim)
      end do
   end subroutine powersums_guarded

   !> Number of atoms currently active after the latest prepare/prepare_subset
   !>
   !> @param[in] self LSF instance
   !> @returns   n    Active-atom count
   pure integer function lsf_active_count(self) result(n)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      n = self%n_active
   end function lsf_active_count

   !> User-space atom index of the i-th currently active atom
   !>
   !> @param[in] self LSF instance
   !> @param[in] i    Active-list index (1 <= i <= active_count())
   !> @returns   idx  User-space atom id
   pure integer function lsf_active_atom(self, i) result(idx)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: i
      idx = self%act_atom(i)
   end function lsf_active_atom

   !* ================================================================================= *!
   !*                          Kernel-facing private helpers                            *!
   !* ================================================================================= *!

   !> Blending weights the kernel is called with, as configured
   !>
   !> - no low-`n_active` guard (see the module header): a shared batch list
   !>   leaves the atom count of a point unknown, and the per-point path uses the
   !>   same weights so both paths stay one model
   !>
   !> @param[in]  self LSF instance
   !> @param[out] s_1  One-body weight
   !> @param[out] s_2  Two-body weight
   !> @param[out] s_3  Three-body weight
   pure subroutine svdw_weights(self, s_1, s_2, s_3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> One-body weight
      real(wp), intent(out) :: s_1
      !> Two-body weight
      real(wp), intent(out) :: s_2
      !> Three-body weight
      real(wp), intent(out) :: s_3

      s_1 = self%param%blend_1b
      s_2 = self%param%blend_2b
      s_3 = self%param%blend_3b
   end subroutine svdw_weights

   !> Fade band handed to stage 1, in `d_A = ||r - R_A|| - R_a`
   !>
   !> - fade off: [[fade_off_start]], [[fade_off_width]], the untruncated model
   !> - fade on: `X_on = 3 fade_start/blend_k`, `Delta = 3 fade_width/blend_k`,
   !>   e-folds of `u_A` turned into distance
   !>
   !> @param[in]  self  LSF instance
   !> @param[out] x_on  Fade onset X_on (Bohr)
   !> @param[out] delta Fade band width Delta (Bohr)
   pure subroutine fade_band(self, x_on, delta)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Fade onset
      real(wp), intent(out) :: x_on
      !> Fade band width
      real(wp), intent(out) :: delta

      if (.not. self%param%fade_enabled) then
         x_on = fade_off_start
         delta = fade_off_width
         return
      end if
      x_on = 3.0_wp*self%param%fade_start/self%param%blend_k
      delta = 3.0_wp*self%param%fade_width/self%param%blend_k
   end subroutine fade_band

   !> Whether the parameters pass their validation, without an error object
   !>
   !> - the pure twin of the parameter `validate`, for pure callers: finite
   !>   `blend_k > 0`, and with the fade on a finite, positive band
   !>
   !> @param[in] self  LSF instance
   !> @returns   valid Parameters usable
   pure function fade_band_valid(self) result(valid)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Parameters usable
      logical :: valid

      associate (a => self%param%fade_start, b => self%param%fade_width, &
                 k => self%param%blend_k)
         valid = ieee_is_finite(k)
         if (valid) valid = k > 0.0_wp
         if (valid .and. self%param%fade_enabled) then
            valid = ieee_is_finite(a) .and. ieee_is_finite(b)
            if (valid) valid = a > 0.0_wp .and. b > 0.0_wp
         end if
      end associate
   end function fade_band_valid

   !> Refuse a fade band changed since `update`
   !>
   !> - called by `prepare`, `prepare_subset`, `evaluate_batch` and
   !>   `vjp_f1_rA_batch` while the fade is on, after the parameter validation
   !> - `param` changed since `update` ([[certificate_current]]): the gate reach
   !>   and the certificates still belong to the old band
   !> - cached reach not the live band end ([[lsf_base_reach_current]]): a
   !>   refresh (`set_screening_threshold`) under other parameters since, the
   !>   recorded band alone does not see it
   !> - with the fade on [[lsf_screening_offset]] is the band end itself,
   !>   whatever `screening_threshold` says
   !>
   !> @param[inout] self  LSF instance
   !> @param[out]   error Band changed since `update`
   subroutine check_fade(self, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Band changed since `update`
      type(error_type), allocatable, intent(out) :: error

      if (.not. (certificate_current(self) .and. self%reach_current())) then
         call fatal_error(error, "SvdW fade parameters changed after update; call update")
         return
      end if
   end subroutine check_fade

   !> Per-atom kind tensors of active slot `ia`, guarded against `x == 0`
   !>
   !> The kernel divides by the distance from the first derivative onwards, and on
   !> a nucleus that distance is zero with the direction `(r - R_A)/x` undefined:
   !>
   !> - the convention here, the same one the retired SSD fill used, takes it as
   !>   zero
   !> - the value tensor still carries the kernel's own atom function, while every
   !>   derivative order of this atom drops out of the power sums
   !> - geometry from stage 1's cache, fade band from [[fade_band]]
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..4)
   !> @param[out] at0       Order-0 kind tensor
   !> @param[out] at1       Order-1 kind tensor
   !> @param[out] at2       Order-2 kind tensor
   !> @param[out] at3       Order-3 kind tensor
   !> @param[out] at4       Order-4 kind tensor
   pure subroutine atom_tensors(self, ia, max_deriv, at0, at1, at2, at3, at4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 kind tensor
      real(wp), intent(out) :: at0(nkind)
      !> Order-1 kind tensor
      real(wp), intent(out) :: at1(ndim, nkind)
      !> Order-2 kind tensor
      real(wp), intent(out) :: at2(ndim, ndim, nkind)
      !> Order-3 kind tensor
      real(wp), intent(out) :: at3(ndim, ndim, ndim, nkind)
      !> Order-4 kind tensor
      real(wp), intent(out) :: at4(ndim, ndim, ndim, ndim, nkind)

      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      call fade_band(self, x_on, delta)
      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                             max_deriv, at0, at1, at2, at3, at4, self%act_x(ia), &
                             self%act_inv_x(ia), self%act_u(ia), x_on, delta)
         return
      end if

      ! On the nucleus: order 0 is still the kernel's own atom function, the rest is zero
      call svdw_atom_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                          0, at0, at1, at2, at3, at4, self%act_x(ia), self%act_inv_x(ia), &
                          self%act_u(ia), x_on, delta)
      if (max_deriv >= 1) at1 = 0.0_wp
      if (max_deriv >= 2) at2 = 0.0_wp
      if (max_deriv >= 3) at3 = 0.0_wp
      if (max_deriv >= 4) at4 = 0.0_wp
   end subroutine atom_tensors

   !> Direction-contracted per-atom kind tensors of active slot `ia`
   !>
   !> Same `x == 0` convention as [[atom_tensors]]: every order of the contracted
   !> tensor carries a nuclear derivative, so all of them vanish there
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  v_a       Nuclear displacement direction of this atom [ndim]
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..3)
   !> @param[out] aw0       Order-0 contracted kind tensor
   !> @param[out] aw1       Order-1 contracted kind tensor
   !> @param[out] aw2       Order-2 contracted kind tensor
   !> @param[out] aw3       Order-3 contracted kind tensor
   pure subroutine atom_tangent_tensors(self, ia, v_a, max_deriv, aw0, aw1, aw2, aw3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Nuclear displacement direction of this atom
      real(wp), intent(in) :: v_a(ndim)
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 contracted kind tensor
      real(wp), intent(out) :: aw0(nkind)
      !> Order-1 contracted kind tensor
      real(wp), intent(out) :: aw1(ndim, nkind)
      !> Order-2 contracted kind tensor
      real(wp), intent(out) :: aw2(ndim, ndim, nkind)
      !> Order-3 contracted kind tensor
      real(wp), intent(out) :: aw3(ndim, ndim, ndim, nkind)

      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      if (self%act_x(ia) > 0.0_wp) then
         call fade_band(self, x_on, delta)
         call svdw_atom_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                     self%param%blend_k, v_a, max_deriv, aw0, aw1, aw2, aw3, &
                                     self%act_x(ia), self%act_inv_x(ia), self%act_u(ia), &
                                     x_on, delta)
         return
      end if

      aw0 = 0.0_wp
      if (max_deriv >= 1) aw1 = 0.0_wp
      if (max_deriv >= 2) aw2 = 0.0_wp
      if (max_deriv >= 3) aw3 = 0.0_wp
   end subroutine atom_tangent_tensors

   !> Radius-contracted per-atom kind tensors of active slot `ia`
   !>
   !> The radius half of the joint `(v, vr)` contraction: `awr{n} = vr_a atd{n}`,
   !> `(blend_k m_kind / 3) vr_a at{n}` on the plateau. It lands in the same
   !> family as [[atom_tangent_tensors]], so the caller accumulates both into one
   !> `ws*` and every aggregate kernel routine picks the radius term up unchanged
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  vr_a      Radius direction of this atom
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..3)
   !> @param[out] awr0      Order-0 radius-contracted kind tensor
   !> @param[out] awr1      Order-1 radius-contracted kind tensor
   !> @param[out] awr2      Order-2 radius-contracted kind tensor
   !> @param[out] awr3      Order-3 radius-contracted kind tensor
   pure subroutine atom_radius_tangent_tensors(self, ia, vr_a, max_deriv, &
                                               awr0, awr1, awr2, awr3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Radius direction of this atom
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 radius-contracted kind tensor
      real(wp), intent(out) :: awr0(nkind)
      !> Order-1 radius-contracted kind tensor
      real(wp), intent(out) :: awr1(ndim, nkind)
      !> Order-2 radius-contracted kind tensor
      real(wp), intent(out) :: awr2(ndim, ndim, nkind)
      !> Order-3 radius-contracted kind tensor
      real(wp), intent(out) :: awr3(ndim, ndim, ndim, nkind)

      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      call fade_band(self, x_on, delta)
      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_radius_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                            self%param%blend_k, vr_a, max_deriv, &
                                            awr0, awr1, awr2, awr3, self%act_x(ia), &
                                            self%act_inv_x(ia), self%act_u(ia), x_on, delta)
         return
      end if

      ! On the nucleus: order 0 is the kernel's own scaled atom function, the rest is zero
      call svdw_atom_radius_tangent_eval(self%act_d(:, ia), self%act_radius(ia), &
                                         self%param%blend_k, vr_a, 0, &
                                         awr0, awr1, awr2, awr3, self%act_x(ia), &
                                         self%act_inv_x(ia), self%act_u(ia), x_on, delta)
      if (max_deriv >= 1) awr1 = 0.0_wp
      if (max_deriv >= 2) awr2 = 0.0_wp
      if (max_deriv >= 3) awr3 = 0.0_wp
   end subroutine atom_radius_tangent_tensors

   !> Per-atom radius tensors `atd = d/dR_a at` of active slot `ia`
   !>
   !> The spatial derivatives of `-dg/dd`, what every radius aggregate reads in
   !> place of `at`; `(blend_k m_kind / 3) at` on the plateau. Same `x == 0`
   !> convention as [[atom_tensors]]
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..3)
   !> @param[out] atd0      Order-0 radius kind tensor
   !> @param[out] atd1      Order-1 radius kind tensor
   !> @param[out] atd2      Order-2 radius kind tensor
   !> @param[out] atd3      Order-3 radius kind tensor
   pure subroutine atom_rad_tensors(self, ia, max_deriv, atd0, atd1, atd2, atd3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 radius kind tensor
      real(wp), intent(out) :: atd0(nkind)
      !> Order-1 radius kind tensor
      real(wp), intent(out) :: atd1(ndim, nkind)
      !> Order-2 radius kind tensor
      real(wp), intent(out) :: atd2(ndim, ndim, nkind)
      !> Order-3 radius kind tensor
      real(wp), intent(out) :: atd3(ndim, ndim, ndim, nkind)

      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      call fade_band(self, x_on, delta)
      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_rad_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                                 max_deriv, atd0, atd1, atd2, atd3, self%act_x(ia), &
                                 self%act_inv_x(ia), self%act_u(ia), x_on, delta)
         return
      end if

      ! On the nucleus: order 0 from the kernel, the rest is zero
      call svdw_atom_rad_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                              0, atd0, atd1, atd2, atd3, self%act_x(ia), &
                              self%act_inv_x(ia), self%act_u(ia), x_on, delta)
      if (max_deriv >= 1) atd1 = 0.0_wp
      if (max_deriv >= 2) atd2 = 0.0_wp
      if (max_deriv >= 3) atd3 = 0.0_wp
   end subroutine atom_rad_tensors

   !> Per-atom second radius tensors `atdd = d^2/dR_a^2 at` of active slot `ia`
   !>
   !> The spatial derivatives of `d^2g/dd^2`, read by the diagonal radius-radius
   !> block. Same `x == 0` convention as [[atom_tensors]]
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..2)
   !> @param[out] atdd0     Order-0 second radius kind tensor
   !> @param[out] atdd1     Order-1 second radius kind tensor
   !> @param[out] atdd2     Order-2 second radius kind tensor
   pure subroutine atom_radrad_tensors(self, ia, max_deriv, atdd0, atdd1, atdd2)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 second radius kind tensor
      real(wp), intent(out) :: atdd0(nkind)
      !> Order-1 second radius kind tensor
      real(wp), intent(out) :: atdd1(ndim, nkind)
      !> Order-2 second radius kind tensor
      real(wp), intent(out) :: atdd2(ndim, ndim, nkind)

      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      call fade_band(self, x_on, delta)
      if (self%act_x(ia) > 0.0_wp) then
         call svdw_atom_radrad_eval(self%act_d(:, ia), self%act_radius(ia), &
                                    self%param%blend_k, max_deriv, atdd0, atdd1, atdd2, &
                                    self%act_x(ia), self%act_inv_x(ia), self%act_u(ia), &
                                    x_on, delta)
         return
      end if

      ! On the nucleus: order 0 from the kernel, the rest is zero
      call svdw_atom_radrad_eval(self%act_d(:, ia), self%act_radius(ia), self%param%blend_k, &
                                 0, atdd0, atdd1, atdd2, self%act_x(ia), &
                                 self%act_inv_x(ia), self%act_u(ia), x_on, delta)
      if (max_deriv >= 1) atdd1 = 0.0_wp
      if (max_deriv >= 2) atdd2 = 0.0_wp
   end subroutine atom_radrad_tensors

   !> Radius derivative `awd = d/dR_a (aw + awr)` of active slot `ia`'s joint share
   !>
   !> `awd{n} = -v_a . atd{n+1} + vr_a atdd{n}`, what [[svdw_radius_hvp]] reads
   !> for the radius derivative of the joint-direction `ws`
   !>
   !> - on a nucleus the nuclear half vanishes with every derivative order of
   !>   [[atom_rad_tensors]], leaving `awd0 = vr_a atdd0` and zero above
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  ia        Active-list index
   !> @param[in]  v_a       Nuclear displacement direction of this atom [ndim]
   !> @param[in]  vr_a      Radius direction of this atom
   !> @param[in]  max_deriv Highest spatial-derivative order to set (0..2)
   !> @param[out] awd0      Order-0 tensor
   !> @param[out] awd1      Order-1 tensor
   !> @param[out] awd2      Order-2 tensor
   pure subroutine atom_tangent_rad_tensors(self, ia, v_a, vr_a, max_deriv, awd0, awd1, awd2)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Active-list index
      integer, intent(in) :: ia
      !> Nuclear displacement direction of this atom
      real(wp), intent(in) :: v_a(ndim)
      !> Radius direction of this atom
      real(wp), intent(in) :: vr_a
      !> Highest spatial-derivative order to set
      integer, intent(in) :: max_deriv
      !> Order-0 tensor
      real(wp), intent(out) :: awd0(nkind)
      !> Order-1 tensor
      real(wp), intent(out) :: awd1(ndim, nkind)
      !> Order-2 tensor
      real(wp), intent(out) :: awd2(ndim, ndim, nkind)

      !> Fade band in `d_A`
      real(wp) :: x_on, delta
      !> Second radius tensors of the nucleus fallback
      real(wp) :: atdd0(nkind), atdd1(ndim, nkind), atdd2(ndim, ndim, nkind)

      if (self%act_x(ia) > 0.0_wp) then
         call fade_band(self, x_on, delta)
         call svdw_atom_tangent_rad_eval(self%act_d(:, ia), self%act_radius(ia), &
                                         self%param%blend_k, v_a, vr_a, max_deriv, &
                                         awd0, awd1, awd2, self%act_x(ia), &
                                         self%act_inv_x(ia), self%act_u(ia), x_on, delta)
         return
      end if

      call atom_radrad_tensors(self, ia, 0, atdd0, atdd1, atdd2)
      awd0 = vr_a*atdd0
      if (max_deriv >= 1) awd1 = 0.0_wp
      if (max_deriv >= 2) awd2 = 0.0_wp
   end subroutine atom_tangent_rad_tensors

   !> Accumulate the direction-contracted power sums for a nuclear direction field
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  v         Nuclear displacement directions [ndim, ncenters]
   !> @param[in]  max_deriv Highest spatial-derivative order to accumulate (0..3)
   !> @param[out] ws0       Order-0 contracted power sums
   !> @param[out] ws1       Order-1 contracted power sums
   !> @param[out] ws2       Order-2 contracted power sums
   !> @param[out] ws3       Order-3 contracted power sums
   !> @param[in]  vrad      Radius directions [ncenters] (optional)
   pure subroutine tangent_powersums(self, v, max_deriv, ws0, ws1, ws2, ws3, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Highest spatial-derivative order to accumulate
      integer, intent(in) :: max_deriv
      !> Order-0 contracted power sums
      real(wp), intent(out) :: ws0(nkind)
      !> Order-1 contracted power sums
      real(wp), intent(out) :: ws1(ndim, nkind)
      !> Order-2 contracted power sums
      real(wp), intent(out) :: ws2(ndim, ndim, nkind)
      !> Order-3 contracted power sums
      real(wp), intent(out) :: ws3(ndim, ndim, ndim, nkind)
      !> Radius half of the joint direction; absent leaves every operation below
      !> untouched, so a geometry-independent radius model stays bit-for-bit
      real(wp), intent(in), optional :: vrad(:)

      !> Per-atom contracted kind tensors
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      !> Per-atom radius-contracted kind tensors
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Active-list index
      integer :: ia

      ws0 = 0.0_wp
      if (max_deriv >= 1) ws1 = 0.0_wp
      if (max_deriv >= 2) ws2 = 0.0_wp
      if (max_deriv >= 3) ws3 = 0.0_wp

      do ia = 1, self%n_active
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), max_deriv, &
                                   aw0, aw1, aw2, aw3)
         ws0 = ws0 + aw0
         if (max_deriv >= 1) ws1 = ws1 + aw1
         if (max_deriv >= 2) ws2 = ws2 + aw2
         if (max_deriv >= 3) ws3 = ws3 + aw3
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), &
                                             max_deriv, awr0, awr1, awr2, awr3)
            ws0 = ws0 + awr0
            if (max_deriv >= 1) ws1 = ws1 + awr1
            if (max_deriv >= 2) ws2 = ws2 + awr2
            if (max_deriv >= 3) ws3 = ws3 + awr3
         end if
      end do
   end subroutine tangent_powersums

   !> Build the per-atom kind-tensor cache the pair loops consume
   !>
   !> The pair kernels read atom A's tensors once per partner, so the cache turns
   !> an O(n_active**2) exponential count back into O(n_active)
   !>
   !> A local allocation rather than a persistent buffer, because:
   !>
   !> - only the three pair accessors need it
   !> - none of them sits on the projection hot path
   !> - a persistent order-4 cache would cost 3.9 kB per atom per thread
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Highest spatial-derivative order to fill (2..4)
   !> @param[out] ca0       Order-0 kind tensors of every active atom
   !> @param[out] ca1       Order-1 kind tensors of every active atom
   !> @param[out] ca2       Order-2 kind tensors of every active atom
   !> @param[out] ca3       Order-3 kind tensors of every active atom
   !> @param[out] ca4       Order-4 kind tensors of every active atom
   subroutine pair_atom_cache(self, max_deriv, ca0, ca1, ca2, ca3, ca4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Highest spatial-derivative order to fill
      integer, intent(in) :: max_deriv
      !> Order-0 kind tensors
      real(wp), allocatable, intent(out) :: ca0(:, :)
      !> Order-1 kind tensors
      real(wp), allocatable, intent(out) :: ca1(:, :, :)
      !> Order-2 kind tensors
      real(wp), allocatable, intent(out) :: ca2(:, :, :, :)
      !> Order-3 kind tensors
      real(wp), allocatable, intent(out) :: ca3(:, :, :, :, :)
      !> Order-4 kind tensors
      real(wp), allocatable, intent(out) :: ca4(:, :, :, :, :, :)

      !> Active-list index and active-atom count
      integer :: ia, n

      n = self%n_active
      allocate (ca0(nkind, n))
      allocate (ca1(ndim, nkind, n))
      allocate (ca2(ndim, ndim, nkind, n))
      allocate (ca3(ndim, ndim, ndim, nkind, n))
      allocate (ca4(ndim, ndim, ndim, ndim, nkind, n))

      do ia = 1, n
         call atom_tensors(self, ia, max_deriv, ca0(:, ia), ca1(:, :, ia), &
                           ca2(:, :, :, ia), ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia))
      end do
   end subroutine pair_atom_cache

   !> Build the per-atom radius-tensor cache the radius pair loops consume
   !>
   !> The `atd` twin of [[pair_atom_cache]], for [[svdw_radpair_block]] and
   !> [[svdw_nucrad_block]]
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Highest spatial-derivative order to fill (0..3)
   !> @param[out] cd0       Order-0 radius kind tensors of every active atom
   !> @param[out] cd1       Order-1 radius kind tensors of every active atom
   !> @param[out] cd2       Order-2 radius kind tensors of every active atom
   !> @param[out] cd3       Order-3 radius kind tensors of every active atom
   subroutine rad_atom_cache(self, max_deriv, cd0, cd1, cd2, cd3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Highest spatial-derivative order to fill
      integer, intent(in) :: max_deriv
      !> Order-0 radius kind tensors
      real(wp), allocatable, intent(out) :: cd0(:, :)
      !> Order-1 radius kind tensors
      real(wp), allocatable, intent(out) :: cd1(:, :, :)
      !> Order-2 radius kind tensors
      real(wp), allocatable, intent(out) :: cd2(:, :, :, :)
      !> Order-3 radius kind tensors
      real(wp), allocatable, intent(out) :: cd3(:, :, :, :, :)

      !> Active-list index and active-atom count
      integer :: ia, n

      n = self%n_active
      allocate (cd0(nkind, n))
      allocate (cd1(ndim, nkind, n))
      allocate (cd2(ndim, ndim, nkind, n))
      allocate (cd3(ndim, ndim, ndim, nkind, n))

      do ia = 1, n
         call atom_rad_tensors(self, ia, max_deriv, cd0(:, ia), cd1(:, :, ia), &
                               cd2(:, :, :, ia), cd3(:, :, :, :, ia))
      end do
   end subroutine rad_atom_cache

   !* ================================================================================= *!
   !*                             Pure spatial derivatives                              *!
   !* ================================================================================= *!

   !> Lift the cached per-point sums into the cached jet, `np = 1`
   !>
   !> - stage 2 of the strip kernel at `np = 1`, the statements of the batch path
   !> - writes `jet*` through `order`; outputs above `order` are left stale
   !>
   !> @param[inout] self  LSF instance with the sums of `max_deriv >= order`
   !> @param[in]    order Highest spatial order to write (0..4)
   subroutine point_lift(self, order)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Highest spatial order
      integer, intent(in) :: order

      !> Blending weights
      real(wp) :: s_1, s_2, s_3

      call svdw_weights(self, s_1, s_2, s_3)
      associate (k => self%param%blend_k)
         select case (order)
         case (0)
            call svdw_strip_lift_o0(1, k, s_1, s_2, s_3, self%pk0, self%jet0)
         case (1)
            call svdw_strip_lift_o1(1, k, s_1, s_2, s_3, self%pk0, self%pk1, self%jet0, self%jet1)
         case (2)
            call svdw_strip_lift_o2(1, k, s_1, s_2, s_3, self%pk0, self%pk1, self%pk2, &
                                    self%jet0, self%jet1, self%jet2)
         case (3)
            call svdw_strip_lift_o3(1, k, s_1, s_2, s_3, self%pk0, self%pk1, self%pk2, &
                                    self%pk3, self%jet0, self%jet1, self%jet2, self%jet3)
         case default
            call svdw_strip_lift_o4(1, k, s_1, s_2, s_3, self%pk0, self%pk1, self%pk2, &
                                    self%pk3, self%pk4, self%jet0, self%jet1, self%jet2, &
                                    self%jet3, self%jet4)
         end select
      end associate
   end subroutine point_lift

   !> Whether every output of one jet column through `order` is finite
   !>
   !> - the reached rule shared by `prepare` and [[poison_unreached]]
   !>
   !> @param[in] order  Highest spatial order (0..4)
   !> @param[in] l0     Value
   !> @param[in] l1     Gradient [ndim]
   !> @param[in] l2     Packed Hessian [6]
   !> @param[in] l3     Packed third derivative [10]
   !> @param[in] l4     Packed fourth derivative [15]
   !> @returns   finite Every output through `order` finite
   pure function jet_finite(order, l0, l1, l2, l3, l4) result(finite)
      !> Highest spatial order
      integer, intent(in) :: order
      !> Value
      real(wp), intent(in) :: l0
      !> Gradient
      real(wp), intent(in) :: l1(ndim)
      !> Packed Hessian
      real(wp), intent(in) :: l2(np2)
      !> Packed third derivative
      real(wp), intent(in) :: l3(np3)
      !> Packed fourth derivative
      real(wp), intent(in) :: l4(np4)
      !> Every output through `order` finite
      logical :: finite

      finite = ieee_is_finite(l0)
      if (order >= 1) finite = finite .and. all(ieee_is_finite(l1))
      if (order >= 2) finite = finite .and. all(ieee_is_finite(l2))
      if (order >= 3) finite = finite .and. all(ieee_is_finite(l3))
      if (order >= 4) finite = finite .and. all(ieee_is_finite(l4))
   end function jet_finite

   !> Level-set value only
   !>
   !> @param[in]  self LSF instance
   !> @param[out] val  LSF value at the prepared point
   subroutine lsf_f0(self, val)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out) :: val

      call self%require_deriv(0, "f0")
      val = self%jet0(1)
   end subroutine lsf_f0

   !> Level-set value, spatial gradient and spatial Hessian
   !>
   !> @param[in]  self    LSF instance
   !> @param[out] lsf0    LSF value (optional)
   !> @param[out] lsf1_r  Spatial gradient (optional)
   !> @param[out] lsf2_rr Packed spatial Hessian [6] (optional)
   subroutine lsf_f012_r(self, lsf0, lsf1_r, lsf2_rr)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Packed spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:)

      !> Highest spatial-derivative order requested
      integer :: md

      md = 0
      if (present(lsf1_r)) md = 1
      if (present(lsf2_rr)) md = 2
      call self%require_deriv(md, "f012_r")

      if (present(lsf0)) lsf0 = self%jet0(1)
      if (present(lsf1_r)) lsf1_r = self%jet1(:, 1)
      if (present(lsf2_rr)) lsf2_rr = self%jet2(:, 1)
   end subroutine lsf_f012_r

   !> Third spatial derivative, plus optionally the lower orders
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] lsf0     LSF value (optional)
   !> @param[out] lsf1_r   Spatial gradient (optional)
   !> @param[out] lsf2_rr  Packed spatial Hessian [6] (optional)
   !> @param[out] lsf3_rrr Packed third spatial derivative [10]
   subroutine lsf_f3_rrr(self, lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value
      real(wp), intent(out), optional :: lsf0
      !> Spatial gradient
      real(wp), intent(out), optional :: lsf1_r(:)
      !> Packed spatial Hessian
      real(wp), intent(out), optional :: lsf2_rr(:)
      !> Packed third spatial derivative
      real(wp), intent(out) :: lsf3_rrr(:)

      call self%require_deriv(3, "f3_rrr")

      if (present(lsf0)) lsf0 = self%jet0(1)
      if (present(lsf1_r)) lsf1_r = self%jet1(:, 1)
      if (present(lsf2_rr)) lsf2_rr = self%jet2(:, 1)
      lsf3_rrr = self%jet3(:, 1)
   end subroutine lsf_f3_rrr

   !> Fourth spatial derivative
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf4_rrrr Packed fourth spatial derivative [15]
   subroutine lsf_f4_rrrr(self, lsf4_rrrr)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Packed fourth spatial derivative
      real(wp), intent(out) :: lsf4_rrrr(:)

      call self%require_deriv(4, "f4_rrrr")

      lsf4_rrrr = self%jet4(:, 1)
   end subroutine lsf_f4_rrrr

   !* ================================================================================= *!
   !*                                 Batched evaluation                                *!
   !* ================================================================================= *!

   !> Level-set jet at a batch of points through the strip kernel
   !>
   !> Override of [[lsf_base_evaluate_batch]], same contract:
   !>
   !> - points cut into strips of `svdw_strip_width(order)`; per strip stage 1
   !>   over one shared atom list, then the lift; `order` alone governs (no
   !>   `require_deriv`, `max_deriv` untouched)
   !> - `candidate_indices` present: the shared list of every strip, sorted
   !>   ascending in candidate space before use; absent: per strip
   !>   `candidates_in_box` over the strip's bounding box, a superset of every
   !>   point's active set, the kernel's per-point gate mask making the result
   !>   the full screen per point
   !> - unreached column (`unreached(column)`): `lsf0` the finite `S` where only a
   !>   derivative overflowed, `+Inf` otherwise, every derivative through `order` a
   !>   quiet NaN ([[poison_unreached]]); unreached means no atom
   !>   within the band end OR a non-finite lift with `S > 0` or NaN:
   !>   - before the lift: `p1 == 0` after stage 1, no atom within reach; its
   !>     sums are patched to three unit atoms so the lift stays finite
   !>     ([[patch_unreached]])
   !>   - after the lift: a non-finite output through `order` with `S > 0` or NaN, e.g.
   !>     `Z` underflowing to zero while `p1 > 0` (`u_A` below ~1e-103, reachable
   !>     with screening disabled; the inline `exp` clamp at `exp(-708)` gives
   !>     `Z = 0` too), or `Z = s3 e3 = 0` at `blend_1b = 0` with one or two
   !>     atoms in reach
   !>   - [[lsf_svdw_far_exclusion_radius]] covers the first case only, which
   !>     every unreached column is when `blend_1b > 0` (see there)
   !> - a non-finite output with `S <= 0`: `error` naming the first such
   !>   column of the strip
   !> - on-nucleus strip (stage 1's `n_zero > 0`): every point of the strip
   !>   re-evaluated on the guarded per-point path with `max_deriv = order`
   !> - the parameters are validated first, and with the fade on the band checked
   !>   against `update` ([[check_fade]]); a refusal is an `error` with no column
   !>   written
   !> - the per-point cache is invalidated (`prepared_deriv = -1`): accessors
   !>   abort until the next `prepare`
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
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   !> @param[out]   error             Invalid fade band, or first column with a jet not
   !>                                 finite at `S <= 0`
   subroutine lsf_evaluate_batch(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, lsf3_rrr, &
                                 lsf4_rrrr, candidate_indices, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
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
      !> Shared atom list
      integer, intent(in), optional :: candidate_indices(:)
      !> Invalid fade band or first column with a jet not finite at `S <= 0`
      type(error_type), allocatable, intent(out) :: error

      !> Lift outputs of one strip, packed
      real(wp) :: t0(max_strip), t1(ndim, max_strip), t2(np2, max_strip)
      real(wp) :: t3(np3, max_strip), t4(np4, max_strip)
      !> Column has a finite jet, column jet not finite at `S <= 0`
      logical :: reached(max_strip), singular(max_strip)
      !> Strip bounding box
      real(wp) :: lo(ndim), hi(ndim)
      !> Point count, strip width, strip bounds and size, list size, pair count
      integer :: n, width, s0, s1, np, nc, n_zero
      !> Point and column indices
      integer :: p, i
      !> Column label and list length for the error message
      character(len=16) :: col, ncol

      n = size(points, 2)
      call lsf_batch_require_outputs("SvdW", points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                     lsf3_rrr, lsf4_rrrr)
      call self%begin_batch(n)
      if (n == 0) return
      if (.not. allocated(self%cand_screen)) then
         error stop "moist DROP LSF SvdW: evaluate_batch before update"
      end if
      call self%param%validate(error)
      if (.not. allocated(error) .and. self%param%fade_enabled) call check_fade(self, error)
      if (allocated(error)) return
      ! The batch overwrites the per-point cache: accessors abort until `prepare`
      self%prepared_deriv = -1

      width = svdw_strip_width(order)
      nc = 0
      if (present(candidate_indices)) then
         call lsf_sorted_candidates(candidate_indices, self%bat_cand, nc, self%ncenters)
         call gather_list(self, nc)
      end if

      do s0 = 1, n, width
         s1 = min(n, s0 + width - 1)
         np = s1 - s0 + 1
         if (.not. present(candidate_indices)) then
            lo = points(:, s0)
            hi = points(:, s0)
            do i = s0 + 1, s1
               lo = min(lo, points(:, i))
               hi = max(hi, points(:, i))
            end do
            ! The cell list itself, a subobject disjoint from `bat_cand`
            call self%cand_cells%query_box(lo, hi, self%bat_cand, nc)
            call gather_list(self, nc)
         end if

         call strip_eval(self, order, np, points(:, s0:s1), nc, t0, t1, t2, t3, t4, reached, &
                         n_zero)
         if (n_zero > 0) call strip_guarded(self, order, np, points(:, s0:s1), nc, &
                                            t0, t1, t2, t3, t4, reached)
         call poison_unreached(order, np, t0, t1, t2, t3, t4, reached, singular)

         do p = 1, np
            i = s0 + p - 1
            lsf0(i) = t0(p)
            if (order >= 1) lsf1_r(:, i) = t1(:, p)
            if (order >= 2) lsf2_rr_packed(:, i) = t2(:, p)
            if (order >= 3) lsf3_rrr(:, i) = t3(:, p)
            if (order >= 4) lsf4_rrrr(:, i) = t4(:, p)
            self%batch_unreached(i) = .not. (reached(p) .or. singular(p))
         end do

         if (any(singular(1:np))) then
            p = findloc(singular(1:np), .true., dim=1)
            write (col, "(i0)") s0 + p - 1
            write (ncol, "(i0)") nc
            call fatal_error(error, "evaluate_batch: point "//trim(col)// &
                             " has a level-set jet not finite with S <= 0 (SvdW, "// &
                             trim(ncol)//" candidates "// &
                             merge("given   ", "from box", present(candidate_indices))//")")
            return
         end if
      end do
   end subroutine lsf_evaluate_batch

   !> Gather centers, radii and gate reach of the shared list `bat_cand(1:nc)`
   !>
   !> @param[inout] self LSF instance
   !> @param[in]    nc   Length of the shared list
   pure subroutine gather_list(self, nc)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
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

   !> Both stages of the strip kernel on one strip over the gathered list
   !>
   !> - `n_zero > 0` (a point on a nucleus): returns after stage 1, outputs
   !>   undefined; the caller re-evaluates the strip
   !> - unreached columns (`p1 == 0`) patched to three unit atoms before the lift
   !>
   !> Stage-1 sums live in fixed buffers of the widest strip; the explicit-shape
   !> kernel dummies see them as `(np, ...)` by sequence association, so nothing
   !> is allocated per strip
   !>
   !> @param[in]  self    LSF instance with `bat_*` gathered
   !> @param[in]  order   Highest spatial order (0..4)
   !> @param[in]  np      Points in the strip
   !> @param[in]  points  Strip points [ndim, np]
   !> @param[in]  nc      Length of the shared list
   !> @param[out] l0      Values [np]
   !> @param[out] l1      Gradients [ndim, np]
   !> @param[out] l2      Packed Hessians [6, np]
   !> @param[out] l3      Packed third derivatives [10, np]
   !> @param[out] l4      Packed fourth derivatives [15, np]
   !> @param[out] reached Column reached at least one atom [np]
   !> @param[out] n_zero  Point-atom pairs at zero distance
   subroutine strip_eval(self, order, np, points, nc, l0, l1, l2, l3, l4, reached, n_zero)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Highest spatial order
      integer, intent(in) :: order
      !> Points in the strip
      integer, intent(in) :: np
      !> Strip points
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
      !> Column reached at least one atom
      logical, intent(out) :: reached(np)
      !> Point-atom pairs at zero distance
      integer, intent(out) :: n_zero

      !> Stage-1 sums of the widest strip, seen as `(np, ..., nkind)` by the kernel
      real(wp) :: b0(max_strip*nkind), b1(max_strip*ndim*nkind)
      real(wp) :: b2(max_strip*np2*nkind), b3(max_strip*np3*nkind), b4(max_strip*np4*nkind)
      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Accumulate the kind-4 sums
      logical :: iq
      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      iq = self%param%blend_2b /= 0.0_wp
      call fade_band(self, x_on, delta)
      associate (k => self%param%blend_k, c => self%bat_center, r => self%bat_radius, &
                 q => self%bat_bsq)
         select case (order)
         case (0)
            call svdw_strip_powersums_o0(np, points, nc, c, r, q, k, iq, x_on, &
                                         delta, b0, n_zero)
         case (1)
            call svdw_strip_powersums_o1(np, points, nc, c, r, q, k, iq, x_on, &
                                         delta, b0, b1, n_zero)
         case (2)
            call svdw_strip_powersums_o2(np, points, nc, c, r, q, k, iq, x_on, &
                                         delta, b0, b1, b2, n_zero)
         case (3)
            call svdw_strip_powersums_o3(np, points, nc, c, r, q, k, iq, x_on, &
                                         delta, b0, b1, b2, b3, n_zero)
         case default
            call svdw_strip_powersums_o4(np, points, nc, c, r, q, k, iq, x_on, &
                                         delta, b0, b1, b2, b3, b4, n_zero)
         end select
      end associate
      if (n_zero > 0) return

      call patch_unreached(np, b0, reached)

      call svdw_weights(self, s_1, s_2, s_3)
      associate (k => self%param%blend_k)
         select case (order)
         case (0)
            call svdw_strip_lift_o0(np, k, s_1, s_2, s_3, b0, l0)
         case (1)
            call svdw_strip_lift_o1(np, k, s_1, s_2, s_3, b0, b1, l0, l1)
         case (2)
            call svdw_strip_lift_o2(np, k, s_1, s_2, s_3, b0, b1, b2, l0, l1, l2)
         case (3)
            call svdw_strip_lift_o3(np, k, s_1, s_2, s_3, b0, b1, b2, b3, l0, l1, l2, l3)
         case default
            call svdw_strip_lift_o4(np, k, s_1, s_2, s_3, b0, b1, b2, b3, b4, l0, l1, l2, &
                                    l3, l4)
         end select
      end associate
   end subroutine strip_eval

   !> Mark the unreached columns of a strip and patch their sums to three unit atoms
   !>
   !> - unreached: `p1 == 0`, every atom of the shared list masked by the gate
   !> - the derivative sums of such a column are exact zeros already; the value
   !>   sums become `p1 = p2 = p3 = q = 3`, so `Z = 3 (s_1 + s_2) + s_3 > 0`
   !>   for any admissible weights (one unit atom gave `Z = s_1`, zero at
   !>   `blend_1b = 0`) and the lift divides by nothing that vanishes
   !> - the patched outputs are overwritten by [[poison_unreached]]
   !>
   !> @param[in]    np      Points in the strip
   !> @param[inout] ps0     Order-0 sums of the strip [np, nkind]
   !> @param[out]   reached Column reached at least one atom [np]
   pure subroutine patch_unreached(np, ps0, reached)
      !> Points in the strip
      integer, intent(in) :: np
      !> Order-0 sums of the strip
      real(wp), intent(inout) :: ps0(np, nkind)
      !> Column reached at least one atom
      logical, intent(out) :: reached(np)

      !> Point index
      integer :: p

      do p = 1, np
         reached(p) = ps0(p, 1) /= 0.0_wp
         if (reached(p)) cycle
         ps0(p, :) = 3.0_wp
      end do
   end subroutine patch_unreached

   !> Classify the non-finite columns and fill every unreached column with the exterior limit
   !>
   !> - a reached column whose outputs through `order` are not all finite becomes
   !>   singular when `S <= 0`, unreached otherwise (`S > 0`, `+Inf` from `Z`
   !>   underflowed to zero while `p1 > 0`, or NaN from a round-off `Z <= 0`);
   !>   `reached` false for both
   !> - every unreached column: every derivative through `order` a quiet NaN; `S` kept
   !>   where finite (`> 0`, only a derivative overflowed), else `+Inf` (`-log 0`, also
   !>   for a column with no atom in reach); a singular column is left as lifted
   !> - outputs above `order` are left alone
   !>
   !> @param[in]    order    Highest spatial order (0..4)
   !> @param[in]    np       Points in the strip
   !> @param[inout] l0       Values [np]
   !> @param[inout] l1       Gradients [ndim, np]
   !> @param[inout] l2       Packed Hessians [6, np]
   !> @param[inout] l3       Packed third derivatives [10, np]
   !> @param[inout] l4       Packed fourth derivatives [15, np]
   !> @param[inout] reached  Column reached at least one atom on entry, has a finite jet
   !>                        on exit [np]
   !> @param[out]   singular Column jet not finite with `S <= 0` [np]
   pure subroutine poison_unreached(order, np, l0, l1, l2, l3, l4, reached, singular)
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
      !> Column reached at least one atom on entry, has a finite jet on exit
      logical, intent(inout) :: reached(np)
      !> Column jet not finite with `S <= 0`
      logical, intent(out) :: singular(np)

      !> Quiet NaN
      real(wp) :: qnan
      !> Point index
      integer :: p
      !> Column keeps its finite `S`
      logical :: keep

      qnan = ieee_value(qnan, ieee_quiet_nan)
      do p = 1, np
         singular(p) = .false.
         keep = .false.
         if (reached(p)) then
            reached(p) = jet_finite(order, l0(p), l1(:, p), l2(:, p), l3(:, p), l4(:, p))
            singular(p) = .not. reached(p) .and. l0(p) <= 0.0_wp
            keep = ieee_is_finite(l0(p))
         end if
         if (reached(p) .or. singular(p)) cycle
         if (.not. keep) l0(p) = ieee_value(qnan, ieee_positive_inf)
         if (order >= 1) l1(:, p) = qnan
         if (order >= 2) l2(:, p) = qnan
         if (order >= 3) l3(:, p) = qnan
         if (order >= 4) l4(:, p) = qnan
      end do
   end subroutine poison_unreached

   !> Re-evaluate a strip point by point on the guarded per-point path
   !>
   !> - taken when stage 1 found a point on a nucleus
   !> - [[lsf_svdw_screen]] over the strip's shared list at `max_deriv = order`,
   !>   its cached jet copied out; `max_deriv` restored and the per-point cache
   !>   invalidated on exit
   !>
   !> @param[inout] self    LSF instance with `bat_cand(1:nc)` set
   !> @param[in]    order   Highest spatial order (0..4)
   !> @param[in]    np      Points in the strip
   !> @param[in]    points  Strip points [ndim, np]
   !> @param[in]    nc      Length of the shared list
   !> @param[out]   l0      Values [np]
   !> @param[out]   l1      Gradients [ndim, np]
   !> @param[out]   l2      Packed Hessians [6, np]
   !> @param[out]   l3      Packed third derivatives [10, np]
   !> @param[out]   l4      Packed fourth derivatives [15, np]
   !> @param[out]   reached Column reached at least one atom [np]
   subroutine strip_guarded(self, order, np, points, nc, l0, l1, l2, l3, l4, reached)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Highest spatial order
      integer, intent(in) :: order
      !> Points in the strip
      integer, intent(in) :: np
      !> Strip points
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
      !> Column reached at least one atom
      logical, intent(out) :: reached(np)

      !> Caller's per-point order and point index
      integer :: md_saved, p
      !> Jet not finite inside; the copied jet lets [[poison_unreached]] see it
      logical :: singular

      md_saved = self%max_deriv
      self%max_deriv = order
      do p = 1, np
         ! `bat_cand` aliases a component of `self`: legal as the screen never touches it
         call lsf_svdw_screen(self, points(:, p), self%bat_cand(1:nc), singular)
         reached(p) = .not. self%prepared_unreached
         if (.not. reached(p)) cycle
         l0(p) = self%jet0(1)
         l1(:, p) = self%jet1(:, 1)
         l2(:, p) = self%jet2(:, 1)
         l3(:, p) = self%jet3(:, 1)
         l4(:, p) = self%jet4(:, 1)
      end do
      self%max_deriv = md_saved
      self%prepared_deriv = -1
   end subroutine strip_guarded

   !* ================================================================================= *!
   !*                     Batched jet-contracted nuclear gradient                       *!
   !* ================================================================================= *!

   !> Batched jet-contracted nuclear gradient through the strip kernel
   !>
   !> Override of [[lsf_base_vjp_f1_rA_batch]], same contract:
   !>
   !> - per strip of `svdw_strip_width(2)` points: stage 1 at order 2 over one shared
   !>   atom list with the geometry cache, the order-2 lift (the reached decision of
   !>   [[lsf_evaluate_batch]]), `svdw_strip_lift_vjp`, then `svdw_strip_vjp` on the
   !>   cache; no second `exp` ([[strip_vjp_pass]])
   !> - the shared list as `evaluate_batch`: `candidate_indices` sorted ascending in
   !>   candidate space, or per strip `candidates_in_box` over the strip's box
   !> - the ids are validated ([[lsf_sorted_candidates]]) before the position map is built
   !> - list present and already ascending: the kernel writes the tile columns in
   !>   place; any other list, and the box mode, scatter the strip result into the
   !>   caller's positions after zeroing the column
   !> - `max_deriv` is not read on the strip path; the on-nucleus fallback sets 2 for
   !>   its duration and restores it
   !> - unreached column (`unreached(column)`: `p1 == 0`, or with `S > 0` or NaN a
   !>   non-finite order-2 lift or a non-finite result column, which catches an adjoint
   !>   that overflows where the lift does not): NaN in `res(:, :, p)`; the other
   !>   columns are computed
   !> - a non-finite lift or result column with `S <= 0` (`1/d` terms of a point within
   !>   about 1e-60 Bohr of a nucleus): `error` as `evaluate_batch`
   !> - on-nucleus column (stage 1's `n_zero > 0`): the strip is discarded, the
   !>   on-nucleus columns go through the guarded per-point path
   !>   ([[vjp_point]], `max_deriv = 2` for its duration) and the remaining columns
   !>   are re-run as a strip of their own, so they match the tile without the
   !>   on-nucleus columns bit for bit
   !> - the parameters are validated first, and with the fade on the band checked
   !>   against `update` ([[check_fade]]); a refusal is an `error` "vjp_f1_rA_batch
   !>   failed at point 1: ..." with no column written
   !> - the per-point cache is invalidated (`prepared_deriv = -1`)
   !>
   !> @param[inout] self              LSF instance
   !> @param[in]    points            Evaluation points (3, np)
   !> @param[in]    w0                Adjoint weights of the value (np)
   !> @param[in]    w1                Adjoint weights of the gradient (3, np)
   !> @param[in]    w2                Packed symmetric part of the Hessian seed (6, np)
   !> @param[out]   res               Contracted nuclear gradient (3, ncand, np)
   !> @param[out]   error             Invalid fade band
   !> @param[in]    candidate_indices Shared atom list in `candidate_space` (optional)
   subroutine lsf_vjp_f1_rA_batch(self, points, w0, w1, w2, res, error, candidate_indices)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
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
      !> Invalid fade band
      type(error_type), allocatable, intent(out) :: error
      !> Shared atom list
      integer, intent(in), optional :: candidate_indices(:)

      !> Strip seeds and points without the on-nucleus columns
      real(wp) :: cp(ndim, vjp_strip), c0(vjp_strip), c1(ndim, vjp_strip)
      real(wp) :: c2(np2, vjp_strip)
      !> Strip bounding box and fill of an unreached column
      real(wp) :: lo(ndim), hi(ndim), nan
      !> Column reached at least one atom, of the strip and of the reduced strip
      logical :: reached(vjp_strip), reached_c(vjp_strip)
      !> Column holds a point on a nucleus
      logical :: nuc(vjp_strip)
      !> Column not finite with `S <= 0`, of the strip and of the reduced strip
      logical :: singular(vjp_strip), singular_c(vjp_strip)
      !> Write the tile columns in place
      logical :: direct
      !> List absent: the box query `bat_cand` is its own tile position
      logical :: box
      !> Strip column of each reduced column
      integer :: lane(vjp_strip)
      !> Point count, list sizes, strip bounds and size
      integer :: n, ncand, nc, s0, s1, np, nk
      !> Column, candidate and atom indices, zero-distance pair counts
      integer :: p, i, k, j, n_zero, n_zero_c
      !> First singular column of the strip
      integer :: first
      !> Column label for the error message
      character(len=16) :: col
      character(len=:), allocatable :: msg

      n = size(points, 2)
      ncand = self%ncenters
      if (present(candidate_indices)) ncand = size(candidate_indices)
      msg = lsf_vjp_batch_shape_message(self%name, points, w0, w1, w2, res, ncand)
      if (len(msg) > 0) error stop msg
      if (n == 0) return
      if (.not. allocated(self%cand_screen)) then
         error stop "moist DROP LSF SvdW: vjp_f1_rA_batch before update"
      end if
      call self%param%validate(error)
      if (.not. allocated(error) .and. self%param%fade_enabled) call check_fade(self, error)
      if (allocated(error)) then
         error%message = "vjp_f1_rA_batch failed at point 1: "//error%message
         return
      end if
      ! The batch overwrites the per-point cache: accessors abort until `prepare`
      self%prepared_deriv = -1
      call self%begin_batch(n)
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
         do j = 1, nc
            self%vjp_apos(j) = self%vjp_pos(self%bat_cand(j))
         end do
         direct = all(self%bat_cand(1:nc) == candidate_indices)
      end if

      do s0 = 1, n, vjp_strip
         s1 = min(n, s0 + vjp_strip - 1)
         np = s1 - s0 + 1
         if (.not. present(candidate_indices)) then
            lo = points(:, s0)
            hi = points(:, s0)
            do i = s0 + 1, s1
               lo = min(lo, points(:, i))
               hi = max(hi, points(:, i))
            end do
            ! The cell list itself, a subobject disjoint from `bat_cand`
            call self%cand_cells%query_box(lo, hi, self%bat_cand, nc)
            call gather_list(self, nc)
         end if

         if (direct) then
            call strip_vjp_pass(self, np, points(:, s0:s1), w0(s0:s1), w1(:, s0:s1), &
                                w2(:, s0:s1), nc, res(:, :, s0:s1), reached, singular, n_zero)
         else
            call strip_vjp_pass(self, np, points(:, s0:s1), w0(s0:s1), w1(:, s0:s1), &
                                w2(:, s0:s1), nc, self%vjp_res, reached, singular, n_zero)
         end if
         first = n + 1
         if (n_zero == 0) then
            do p = 1, np
               if (singular(p)) first = min(first, s0 + p - 1)
            end do
         end if

         if (n_zero > 0) then
            ! On-nucleus columns apart, the rest re-run as a strip of its own
            call nucleus_columns(np, nc, self%vjp_x, nuc)
            nk = 0
            do p = 1, np
               if (nuc(p)) cycle
               nk = nk + 1
               lane(nk) = p
               cp(:, nk) = points(:, s0 + p - 1)
               c0(nk) = w0(s0 + p - 1)
               c1(:, nk) = w1(:, s0 + p - 1)
               c2(:, nk) = w2(:, s0 + p - 1)
            end do
            if (nk > 0) then
               call strip_vjp_pass(self, nk, cp, c0, c1, c2, nc, self%vjp_res, &
                                   reached_c, singular_c, n_zero_c)
               if (n_zero_c > 0) then
                  error stop "moist DROP LSF SvdW: vjp_f1_rA_batch nucleus left in strip"
               end if
               do k = 1, nk
                  if (singular_c(k)) first = min(first, s0 + lane(k) - 1)
               end do
            end if
         end if
         if (first <= n) then
            write (col, "(i0)") first
            call fatal_error(error, "vjp_f1_rA_batch: point "//trim(col)// &
                             " has a level-set jet not finite with S <= 0 ("// &
                             self%name//")")
            return
         end if

         if (n_zero > 0) then
            reached = .false.
            do k = 1, nk
               p = lane(k)
               reached(p) = reached_c(k)
               i = s0 + p - 1
               if (.not. reached(p)) cycle
               if (box) then
                  call scatter_column(nc, self%bat_cand, &
                                      self%vjp_res(3*nc*(k - 1) + 1), res(:, :, i))
               else
                  call scatter_column(nc, self%vjp_apos, &
                                      self%vjp_res(3*nc*(k - 1) + 1), res(:, :, i))
               end if
            end do
            ! Tile position by user-space id over the strip's list, cleared after
            do j = 1, nc
               if (box) then
                  self%vjp_upos(self%cand_to_user(self%bat_cand(j))) = self%bat_cand(j)
               else
                  self%vjp_upos(self%cand_to_user(self%bat_cand(j))) = self%vjp_apos(j)
               end if
            end do
            do p = 1, np
               if (.not. nuc(p)) cycle
               i = s0 + p - 1
               call vjp_point(self, points(:, i), w0(i), w1(:, i), w2(:, i), nc, &
                              self%vjp_upos, res(:, :, i), reached(p))
            end do
            do j = 1, nc
               self%vjp_upos(self%cand_to_user(self%bat_cand(j))) = 0
            end do
         else if (.not. direct) then
            do p = 1, np
               if (.not. reached(p)) cycle
               if (box) then
                  call scatter_column(nc, self%bat_cand, &
                                      self%vjp_res(3*nc*(p - 1) + 1), res(:, :, s0 + p - 1))
               else
                  call scatter_column(nc, self%vjp_apos, &
                                      self%vjp_res(3*nc*(p - 1) + 1), res(:, :, s0 + p - 1))
               end if
            end do
         end if

         do p = 1, np
            if (reached(p)) cycle
            res(:, :, s0 + p - 1) = nan
            self%batch_unreached(s0 + p - 1) = .true.
         end do
      end do
   end subroutine lsf_vjp_f1_rA_batch

   !> Size the vjp strip scratch for the current atom count
   !>
   !> @param[inout] self LSF instance with `ncenters` set
   pure subroutine grow_vjp_scratch(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self

      !> Elements of one cache for a full strip
      integer :: need

      if (allocated(self%vjp_pos)) then
         if (size(self%vjp_pos) < self%ncenters) deallocate (self%vjp_pos, self%vjp_apos, &
                                                             self%vjp_upos)
      end if
      if (.not. allocated(self%vjp_pos)) then
         allocate (self%vjp_pos(max(self%ncenters, 1)), self%vjp_apos(max(self%ncenters, 1)))
         allocate (self%vjp_upos(max(self%ncenters, 1)), source=0)
      end if
      need = vjp_strip*max(self%ncenters, 1)
      if (allocated(self%vjp_x)) then
         if (size(self%vjp_x) >= need) return
         deallocate (self%vjp_x, self%vjp_ix, self%vjp_u, self%vjp_res)
      end if
      allocate (self%vjp_x(need), self%vjp_ix(need), self%vjp_u(need), &
                self%vjp_res(3*need))
   end subroutine grow_vjp_scratch

   !> Stages 1 to 3 of the nuclear vjp on one strip over the gathered list
   !>
   !> - stage 1 at order 2 with the geometry cache into `vjp_x`, `vjp_ix`, `vjp_u`
   !> - `n_zero > 0` (a point on a nucleus): returns after stage 1, `reached` false,
   !>   `rs` undefined, the cache readable by [[nucleus_columns]]
   !> - the order-2 lift decides `reached` as [[lsf_evaluate_batch]]: `p1 == 0`, or a
   !>   non-finite value, gradient or Hessian (`singular` instead with `S <= 0`, as
   !>   [[poison_unreached]]); the sums of such a column are set to
   !>   three unit atoms before the adjoint lift, so it stays finite, and its `rs`
   !>   column is the caller's to discard
   !> - the adjoint carries one more power of `1/Z` than the lift (`1/Z^3`), so it
   !>   overflows for `Z` in about 1e-154 to 1e-103 where the lift is finite: a column
   !>   whose final `rs` is not finite becomes unreached with `S > 0` or NaN, singular
   !>   with `S <= 0` (powers of `1/d` within about 1e-60 Bohr of a nucleus)
   !>
   !> Stage buffers live in fixed arrays of the strip width; the explicit-shape
   !> kernel dummies see them as `(np, ...)` by sequence association
   !>
   !> @param[inout] self    LSF instance with `bat_*` gathered
   !> @param[in]    np      Points in the strip, 1..vjp_strip
   !> @param[in]    points  Strip points [ndim, np]
   !> @param[in]    w0      Value weights [np]
   !> @param[in]    w1      Gradient weights [ndim, np]
   !> @param[in]    w2      Packed symmetric Hessian seeds [6, np]
   !> @param[in]    nc      Length of the shared list
   !> @param[out]   rs      Nuclear gradient per list atom and point [ndim, nc, np]
   !> @param[out]   reached  Column reached at least one atom, lift and `rs` finite [np]
   !> @param[out]   singular Column lift or `rs` not finite with `S <= 0` [np]
   !> @param[out]   n_zero   Point-atom pairs at zero distance
   subroutine strip_vjp_pass(self, np, points, w0, w1, w2, nc, rs, reached, singular, n_zero)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Points in the strip
      integer, intent(in) :: np
      !> Strip points
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
      !> Column lift or `rs` not finite with `S <= 0`
      logical, intent(out) :: singular(np)
      !> Point-atom pairs at zero distance
      integer, intent(out) :: n_zero

      !> Stage-1 sums and their adjoints, seen as `(np, ..., nkind)` by the kernel
      real(wp) :: b0(vjp_strip*nkind), b1(vjp_strip*ndim*nkind), b2(vjp_strip*np2*nkind)
      real(wp) :: m0(vjp_strip*nkind), m1(vjp_strip*ndim*nkind), m2(vjp_strip*np2*nkind)
      !> Order-2 lift of the strip
      real(wp) :: l0(vjp_strip), l1(ndim, vjp_strip), l2(np2, vjp_strip)
      !> Orders above 2 for the shared finiteness rule, never read
      real(wp), parameter :: z3(np3) = 0.0_wp, z4(np4) = 0.0_wp
      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Accumulate the kind-4 sums
      logical :: iq
      !> Fade band in `d_A`
      real(wp) :: x_on, delta
      !> Point index and stage 3's zero-distance count
      integer :: p, n_zero3

      iq = self%param%blend_2b /= 0.0_wp
      call fade_band(self, x_on, delta)
      call svdw_weights(self, s_1, s_2, s_3)
      associate (k => self%param%blend_k, c => self%bat_center, r => self%bat_radius, &
                 q => self%bat_bsq)
         call svdw_strip_powersums_o2(np, points, nc, c, r, q, k, iq, x_on, delta, b0, &
                                      b1, b2, n_zero, cache_x=self%vjp_x, &
                                      cache_inv_x=self%vjp_ix, cache_u=self%vjp_u)
         singular = .false.
         if (n_zero > 0) then
            reached = .false.
            return
         end if

         call patch_unreached(np, b0, reached)
         call svdw_strip_lift_o2(np, k, s_1, s_2, s_3, b0, b1, b2, l0, l1, l2)
         do p = 1, np
            if (.not. reached(p)) cycle
            reached(p) = jet_finite(2, l0(p), l1(:, p), l2(:, p), z3, z4)
            singular(p) = .not. reached(p) .and. l0(p) <= 0.0_wp
         end do
         call clear_unreached(np, reached, b0, b1, b2)

         call svdw_strip_lift_vjp(np, k, s_1, s_2, s_3, b0, b1, b2, w0, w1, w2, m0, m1, &
                                  m2)
         call svdw_strip_vjp(np, points, nc, c, r, k, iq, x_on, delta, m0, m1, m2, &
                             self%vjp_x, self%vjp_ix, self%vjp_u, rs, n_zero3)
      end associate
      n_zero = n_zero3
      ! The adjoint can overflow where the lift does not: the result decides
      do p = 1, np
         if (.not. reached(p)) cycle
         reached(p) = all(ieee_is_finite(rs(:, :, p)))
         singular(p) = .not. reached(p) .and. l0(p) <= 0.0_wp
      end do
   end subroutine strip_vjp_pass

   !> Sums of the unreached columns of a strip set to three unit atoms
   !>
   !> - value sums `p1 = p2 = p3 = q = 3` as [[patch_unreached]], derivative sums zero
   !> - `Z = 3 (s_1 + s_2) + s_3` then, positive for nonnegative weights not all zero
   !>   as in [[patch_unreached]]; other sets may leave Inf or NaN in discarded lanes
   !>
   !> @param[in]    np      Points in the strip
   !> @param[in]    reached Column reached at least one atom and is finite [np]
   !> @param[inout] ps0     Order-0 sums [np, nkind]
   !> @param[inout] ps1     Order-1 sums [np, ndim, nkind]
   !> @param[inout] ps2     Order-2 sums [np, 6, nkind]
   pure subroutine clear_unreached(np, reached, ps0, ps1, ps2)
      !> Points in the strip
      integer, intent(in) :: np
      !> Column reached at least one atom and is finite
      logical, intent(in) :: reached(np)
      !> Order-0 sums
      real(wp), intent(inout) :: ps0(np, nkind)
      !> Order-1 sums
      real(wp), intent(inout) :: ps1(np, ndim, nkind)
      !> Order-2 sums
      real(wp), intent(inout) :: ps2(np, np2, nkind)

      !> Point index
      integer :: p

      do p = 1, np
         if (reached(p)) cycle
         ps0(p, :) = 3.0_wp
         ps1(p, :, :) = 0.0_wp
         ps2(p, :, :) = 0.0_wp
      end do
   end subroutine clear_unreached

   !> Columns of a strip that sit exactly on an atom of the shared list
   !>
   !> @param[in]  np    Points in the strip
   !> @param[in]  nc    Length of the shared list
   !> @param[in]  cache_x Stage 1's distance cache [np, nc]
   !> @param[out] nuc   Column has a zero-distance atom [np]
   pure subroutine nucleus_columns(np, nc, cache_x, nuc)
      !> Points in the strip
      integer, intent(in) :: np
      !> Length of the shared list
      integer, intent(in) :: nc
      !> Stage 1's distance cache
      real(wp), intent(in) :: cache_x(np, nc)
      !> Column has a zero-distance atom
      logical, intent(out) :: nuc(np)

      !> Point index
      integer :: p

      do p = 1, np
         nuc(p) = any(cache_x(p, :) == 0.0_wp)
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

   !> One column of the batched vjp on the guarded per-point path
   !>
   !> - taken for a point on a nucleus: [[lsf_svdw_screen]] over the strip's list at
   !>   `max_deriv = 2`, then [[lsf_vjp_f1_rA]] with the packed seed, the
   !>   active slots scattered by user-space id; `max_deriv` restored and the
   !>   per-point cache invalidated on exit
   !>
   !> @param[inout] self     LSF instance with `bat_cand(1:nc)` set
   !> @param[in]    point    Evaluation point [ndim]
   !> @param[in]    w0       Value weight
   !> @param[in]    w1       Gradient weights [ndim]
   !> @param[in]    w2       Packed symmetric Hessian seed [6]
   !> @param[in]    nc       Length of the shared list
   !> @param[in]    position Candidate position of a user-space atom id, 0 when none
   !> @param[out]   out      Tile column [ndim, ncand], exact zero off the active atoms
   !> @param[out]   reached  The point reached an atom and its jet is finite
   subroutine vjp_point(self, point, w0, w1, w2, nc, position, out, reached)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Value weight
      real(wp), intent(in) :: w0
      !> Gradient weights
      real(wp), intent(in) :: w1(ndim)
      !> Packed symmetric Hessian seed
      real(wp), intent(in) :: w2(:)
      !> Length of the shared list
      integer, intent(in) :: nc
      !> Candidate position of a user-space atom id
      integer, intent(in) :: position(:)
      !> Tile column
      real(wp), intent(out) :: out(:, :)
      !> The point reached an atom and its jet is finite
      logical, intent(out) :: reached

      !> Per-point result over the active slots
      real(wp), allocatable :: work(:, :)
      !> Caller's per-point order, active slot and candidate position
      integer :: md_saved, i, k
      !> Jet not finite inside
      logical :: singular

      md_saved = self%max_deriv
      self%max_deriv = 2
      ! `bat_cand` aliases a component of `self`: legal as the screen never touches it
      call lsf_svdw_screen(self, point, self%bat_cand(1:nc), singular)
      ! The guarded sums keep an on-nucleus jet finite
      if (singular) error stop "moist DROP LSF SvdW: vjp_f1_rA_batch on-nucleus jet not finite"
      reached = .not. self%prepared_unreached
      if (reached) then
         allocate (work(ndim, max(self%ncenters, 1)))
         call lsf_vjp_f1_rA(self, w0, w1, w2, work)
         out = 0.0_wp
         do i = 1, self%n_active
            k = position(self%act_atom(i))
            if (k == 0) error stop "moist DROP LSF SvdW: vjp_f1_rA_batch active atom " &
               //"outside the list"
            out(:, k) = work(:, i)
         end do
      end if
      self%max_deriv = md_saved
      self%prepared_deriv = -1
   end subroutine vjp_point

   !* ================================================================================= *!
   !*                          One-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Mixed third derivative d^3S / (dr^2 dR_A) and its lower orders
   !>
   !> All three outputs are active-indexed: slot `i` belongs to `active_atom(i)`
   !>
   !> @param[in]  self       LSF instance
   !> @param[out] lsf1_rA    dS/dR_A [3, >= active_count()] (optional)
   !> @param[out] lsf2_r_rA  d^2S/(dr dR_A) [3, 3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rA d^3S/(dr^2 dR_A) [6, 3, >= active_count()], packed in the rr pair
   subroutine lsf_f3_rr_rA(self, lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear gradient
      real(wp), intent(out), optional :: lsf1_rA(:, :)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rA(:, :, :)
      !> Mixed third derivative, packed in the rr pair
      real(wp), intent(out) :: lsf3_rr_rA(:, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim, ndim)
      !> Active-list index and nuclear component
      integer :: ia, s

      call self%require_deriv(2, "f3_rr_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 3, at0, at1, at2, at3, at4)
         call svdw_nuclear_eval(self%param%blend_k, s_1, s_2, s_3, &
                                self%ps0, self%ps1, self%ps2, self%ps3, &
                                at0, at1, at2, at3, at4, 2, &
                                f1_rA, f2_r_rA, f3_rr_rA, d4)
         if (present(lsf1_rA)) lsf1_rA(:, ia) = f1_rA
         if (present(lsf2_r_rA)) lsf2_r_rA(:, :, ia) = f2_r_rA
         do s = 1, ndim
            call pack_sym_matrix(f3_rr_rA(:, :, s), lsf3_rr_rA(:, s, ia), ndim)
         end do
      end do
   end subroutine lsf_f3_rr_rA

   !* ================================================================================= *!
   !*                              Radius derivatives                                   *!
   !* ================================================================================= *!

   !> Radius derivative d^3S / (dr^2 dR_a) and its lower orders
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf1_rad     dS/dR_a [>= active_count()] (optional)
   !> @param[out] lsf2_r_rad   d^2S/(dr dR_a) [3, >= active_count()] (optional)
   !> @param[out] lsf3_rr_rad  d^3S/(dr^2 dR_a) [6, >= active_count()], packed in the rr pair
   subroutine lsf_f3_rr_rad(self, lsf1_rad, lsf2_r_rad, lsf3_rr_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Radius gradient
      real(wp), intent(out), optional :: lsf1_rad(:)
      !> Mixed second derivative
      real(wp), intent(out), optional :: lsf2_r_rad(:, :)
      !> Mixed third derivative, packed in the rr pair
      real(wp), intent(out) :: lsf3_rr_rad(:, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom radius tensors
      real(wp) :: atd0(nkind), atd1(ndim, nkind), atd2(ndim, ndim, nkind)
      real(wp) :: atd3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: f1_rad, f2_r_rad(ndim), f3_rr_rad(ndim, ndim)
      real(wp) :: d4(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(2, "f3_rr_rad")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         ! `d/dR_a` consumes no derivative order, so order 2 is all that is
         ! needed here -- one order less than `f3_rr_rA`. `atom_rad_tensors`
         ! declares `atd3` `intent(out)` and then leaves it unset at
         ! `max_deriv = 2`; `svdw_radius_eval` never reads it, but it is still
         ! passed by reference, so define it rather than hand an `intent(in)`
         ! dummy undefined memory
         call atom_rad_tensors(self, ia, 2, atd0, atd1, atd2, atd3)
         atd3 = 0.0_wp
         call svdw_radius_eval(self%param%blend_k, s_1, s_2, s_3, &
                               self%ps0, self%ps1, self%ps2, self%ps3, &
                               atd0, atd1, atd2, atd3, 2, &
                               f1_rad, f2_r_rad, f3_rr_rad, d4)
         if (present(lsf1_rad)) lsf1_rad(ia) = f1_rad
         if (present(lsf2_r_rad)) lsf2_r_rad(:, ia) = f2_r_rad
         call pack_sym_matrix(f3_rr_rad, lsf3_rr_rad(:, ia), ndim)
      end do
   end subroutine lsf_f3_rr_rad

   !> Mixed fourth derivative d^4S / (dr^3 dR_A), active-indexed
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf4_rrr_rA d^4S/(dr^3 dR_A) [10, 3, >= active_count()], packed in the rrr triple
   subroutine lsf_f4_rrr_rA(self, lsf4_rrr_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative, packed in the rrr triple
      real(wp), intent(out) :: lsf4_rrr_rA(:, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), f3_rr_rA(ndim, ndim, ndim)
      real(wp) :: f4_rrr_rA(ndim, ndim, ndim, ndim)
      !> Active-list index and nuclear component
      integer :: ia, s

      call self%require_deriv(3, "f4_rrr_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 4, at0, at1, at2, at3, at4)
         call svdw_nuclear_eval(self%param%blend_k, s_1, s_2, s_3, &
                                self%ps0, self%ps1, self%ps2, self%ps3, &
                                at0, at1, at2, at3, at4, 3, &
                                f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA)
         do s = 1, ndim
            call pack_sym_tensor3(f4_rrr_rA(:, :, :, s), lsf4_rrr_rA(:, s, ia), ndim)
         end do
      end do
   end subroutine lsf_f4_rrr_rA

   !> Normalized level set S/||grad S|| and its nuclear gradient
   !>
   !> @param[in]  self     LSF instance
   !> @param[out] val      S/||grad S||
   !> @param[out] deriv_rA d/dR_A of S/||grad S|| [3, >= active_count()] (optional)
   subroutine lsf_normalized_f01_rA(self, val, deriv_rA)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Normalized level-set value
      real(wp), intent(out) :: val
      !> Nuclear gradient of the normalized value
      real(wp), intent(out), optional :: deriv_rA(:, :)

      !> Blending weights and pure spatial outputs
      real(wp) :: s_1, s_2, s_3, f0, f1_r(ndim)
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> One-nuclear-index outputs of one atom
      real(wp) :: f1_rA(ndim), f2_r_rA(ndim, ndim), d3n(ndim, ndim, ndim)
      real(wp) :: d4n(ndim, ndim, ndim, ndim)
      !> Normalized outputs of one atom
      real(wp) :: norm0, norm1_rA(ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(1, "normalized_f01_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      f0 = self%jet0(1)
      f1_r = self%jet1(:, 1)

      ! The value half needs no nuclear input; the zero dummies keep the single
      ! shared gradient-norm branch of the kernel in charge of the degenerate case
      f1_rA = 0.0_wp
      f2_r_rA = 0.0_wp
      call svdw_normalized_eval(f0, f1_r, f1_rA, f2_r_rA, norm0, norm1_rA)
      val = norm0
      if (.not. present(deriv_rA)) return

      do ia = 1, self%n_active
         call atom_tensors(self, ia, 2, at0, at1, at2, at3, at4)
         call svdw_nuclear_eval(self%param%blend_k, s_1, s_2, s_3, &
                                self%ps0, self%ps1, self%ps2, self%ps3, &
                                at0, at1, at2, at3, at4, 1, &
                                f1_rA, f2_r_rA, d3n, d4n)
         call svdw_normalized_eval(f0, f1_r, f1_rA, f2_r_rA, norm0, norm1_rA)
         deriv_rA(:, ia) = norm1_rA
      end do
   end subroutine lsf_normalized_f01_rA

   !* ================================================================================= *!
   !*                                 Atomic partition                                  *!
   !* ================================================================================= *!

   !> SvdW splits into atomic shares
   !>
   !> @param[in] self LSF instance
   pure function lsf_has_atom_partition(self) result(available)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Always available
      logical :: available

      available = .true.
   end function lsf_has_atom_partition

   !> Add the atomic shares of Z at the prepared point
   !>
   !> @param[inout] self  LSF instance
   !> @param[in]    scale Weight of the point
   !> @param[inout] acc   Per-atom accumulator, user space [ncenters]
   subroutine lsf_atom_partition(self, scale, acc)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self
      !> Weight of the point
      real(wp), intent(in) :: scale
      !> Per-atom accumulator
      real(wp), intent(inout) :: acc(:)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors, only order 0 read
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Atomic shares of Z, active-indexed
      real(wp) :: z_atom(self%n_active)
      !> Power sums p1, p2 and q, and p1 without the atom
      real(wp) :: p1, p2, q, p1_other
      integer :: ia

      call self%require_deriv(0, "atom_partition")

      call svdw_weights(self, s_1, s_2, s_3)
      p1 = self%ps0(1)
      p2 = self%ps0(2)
      q = self%ps0(4)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 0, at0, at1, at2, at3, at4)
         p1_other = p1 - at0(1)
         z_atom(ia) = s_1*at0(3) + 0.5_wp*s_2*at0(4)*(q - at0(4)) &
                      + s_3/6.0_wp*at0(1)*(p1_other*p1_other - (p2 - at0(2)))
      end do
      z_atom = z_atom*(scale/sum(z_atom))
      do ia = 1, self%n_active
         acc(self%act_atom(ia)) = acc(self%act_atom(ia)) + z_atom(ia)
      end do
   end subroutine lsf_atom_partition

   !* ================================================================================= *!
   !*                          Two-nuclear-index derivatives                            *!
   !* ================================================================================= *!

   !> Pure nuclear Hessian d^2S / (dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self      LSF instance
   !> @param[out] lsf2_rArB d^2S/(dR_A dR_B) [3, >= n_act, 3, >= n_act]
   subroutine lsf_f2_rArB(self, lsf2_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear Hessian
      real(wp), intent(out) :: lsf2_rArB(:, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices and active-atom count
      integer :: ia, ib, n

      call self%require_deriv(0, "f2_rArB")

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, 2, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_pair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                        self%ps0, self%ps1, self%ps2, &
                                        ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                        ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia), &
                                        0, b2, b3, b4)
            else
               call svdw_pair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                   ca3(:, :, :, :, ia), &
                                   ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                   ca3(:, :, :, :, ib), 0, b2, b3, b4)
            end if
            lsf2_rArB(:, ia, :, ib) = b2
         end do
      end do
   end subroutine lsf_f2_rArB

   !> Mixed third derivative d^3S / (dr dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf3_r_rArB d^3S/(dr dR_A dR_B) [3, 3, >= n_act, 3, >= n_act]
   subroutine lsf_f3_r_rArB(self, lsf3_r_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_rArB(:, :, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices, active-atom count and nuclear components
      integer :: ia, ib, n, s, t

      call self%require_deriv(1, "f3_r_rArB")

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, 3, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_pair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                        self%ps0, self%ps1, self%ps2, &
                                        ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                        ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia), &
                                        1, b2, b3, b4)
            else
               call svdw_pair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                   ca3(:, :, :, :, ia), &
                                   ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                   ca3(:, :, :, :, ib), 1, b2, b3, b4)
            end if
            do t = 1, ndim
               do s = 1, ndim
                  lsf3_r_rArB(:, s, ia, t, ib) = b3(:, s, t)
               end do
            end do
         end do
      end do
   end subroutine lsf_f3_r_rArB

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_B), active-indexed in both nuclei
   !>
   !> @param[in]  self         LSF instance
   !> @param[out] lsf4_rr_rArB d^4S/(dr^2 dR_A dR_B) [6, 3, >= n_act, 3, >= n_act], packed in the rr pair
   subroutine lsf_f4_rr_rArB(self, lsf4_rr_rArB)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative, packed in the rr pair
      real(wp), intent(out) :: lsf4_rr_rArB(:, :, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Result blocks of one pair
      real(wp) :: b2(ndim, ndim), b3(ndim, ndim, ndim), b4(ndim, ndim, ndim, ndim)
      !> Active-list indices, active-atom count and nuclear components
      integer :: ia, ib, n, s, t

      call self%require_deriv(2, "f4_rr_rArB")

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, 4, ca0, ca1, ca2, ca3, ca4)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_pair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                        self%ps0, self%ps1, self%ps2, &
                                        ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                        ca3(:, :, :, :, ia), ca4(:, :, :, :, :, ia), &
                                        2, b2, b3, b4)
            else
               call svdw_pair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                   ca3(:, :, :, :, ia), &
                                   ca0(:, ib), ca1(:, :, ib), ca2(:, :, :, ib), &
                                   ca3(:, :, :, :, ib), 2, b2, b3, b4)
            end if
            do t = 1, ndim
               do s = 1, ndim
                  call pack_sym_matrix(b4(:, :, s, t), lsf4_rr_rArB(:, s, ia, t, ib), ndim)
               end do
            end do
         end do
      end do
   end subroutine lsf_f4_rr_rArB

   !* ================================================================================= *!
   !*                     Two-radius and nuclear-radius derivatives                     *!
   !* ================================================================================= *!

   !> Shared pair loop of the two-radius block
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Spatial-derivative order of the requested output
   !> @param[in]  caller    Name reported by `require_deriv` on an under-request
   !> @param[out] res2      d^2S/(dR_a dR_b) (optional)
   !> @param[out] res3      d^3S/(dr dR_a dR_b) (optional)
   !> @param[out] res4      d^4S/(dr^2 dR_a dR_b), packed in the rr pair (optional)
   subroutine svdw_radpair_block(self, max_deriv, caller, res2, res3, res4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Spatial-derivative order of the requested output
      integer, intent(in) :: max_deriv
      !> Name reported on an under-request
      character(len=*), intent(in) :: caller
      !> Radius Hessian
      real(wp), intent(out), optional :: res2(:, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: res3(:, :, :)
      !> Mixed fourth derivative, packed in the rr pair
      real(wp), intent(out), optional :: res4(:, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom radius tensors, orders 0..3
      real(wp), allocatable :: cd0(:, :), cd1(:, :, :), cd2(:, :, :, :)
      real(wp), allocatable :: cd3(:, :, :, :, :)
      !> Second radius tensors of the diagonal atom
      real(wp) :: atdd0(nkind), atdd1(ndim, nkind), atdd2(ndim, ndim, nkind)
      !> Result blocks of one pair; a radius slot adds no Cartesian index
      real(wp) :: b2, b3(ndim), b4(ndim, ndim)
      !> Active-list indices and active-atom count
      integer :: ia, ib, n

      call self%require_deriv(max_deriv, caller)

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call rad_atom_cache(self, max_deriv, cd0, cd1, cd2, cd3)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call atom_radrad_tensors(self, ia, max_deriv, atdd0, atdd1, atdd2)
               call svdw_radpair_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                           self%ps0, self%ps1, self%ps2, &
                                           cd0(:, ia), cd1(:, :, ia), cd2(:, :, :, ia), &
                                           atdd0, atdd1, atdd2, max_deriv, b2, b3, b4)
            else
               call svdw_radpair_eval(self%param%blend_k, s_1, s_2, s_3, &
                                      self%ps0, self%ps1, self%ps2, &
                                      cd0(:, ia), cd1(:, :, ia), cd2(:, :, :, ia), &
                                      cd0(:, ib), cd1(:, :, ib), cd2(:, :, :, ib), &
                                      max_deriv, b2, b3, b4)
            end if
            if (present(res2)) res2(ia, ib) = b2
            if (present(res3)) res3(:, ia, ib) = b3
            if (present(res4)) call pack_sym_matrix(b4, res4(:, ia, ib), ndim)
         end do
      end do
   end subroutine svdw_radpair_block

   !> Shared pair loop of the nuclear-radius block
   !>
   !> The first atom index carries the *position* derivative and the second the
   !> *radius* one, so unlike [[svdw_radpair_block]] the result is not symmetric
   !> under exchanging them; the loop therefore visits every ordered pair
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  max_deriv Spatial-derivative order of the requested output
   !> @param[in]  caller    Name reported by `require_deriv` on an under-request
   !> @param[out] res2      d^2S/(dR_A dR_b) (optional)
   !> @param[out] res3      d^3S/(dr dR_A dR_b) (optional)
   !> @param[out] res4      d^4S/(dr^2 dR_A dR_b), packed in the rr pair (optional)
   subroutine svdw_nucrad_block(self, max_deriv, caller, res2, res3, res4)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Spatial-derivative order of the requested output
      integer, intent(in) :: max_deriv
      !> Name reported on an under-request
      character(len=*), intent(in) :: caller
      !> Nuclear-radius Hessian
      real(wp), intent(out), optional :: res2(:, :, :)
      !> Mixed third derivative
      real(wp), intent(out), optional :: res3(:, :, :, :)
      !> Mixed fourth derivative, packed in the rr pair
      real(wp), intent(out), optional :: res4(:, :, :, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Cached per-atom kind tensors, orders 0..4
      real(wp), allocatable :: ca0(:, :), ca1(:, :, :), ca2(:, :, :, :)
      real(wp), allocatable :: ca3(:, :, :, :, :), ca4(:, :, :, :, :, :)
      !> Cached per-atom radius tensors, orders 0..3
      real(wp), allocatable :: cd0(:, :), cd1(:, :, :), cd2(:, :, :, :)
      real(wp), allocatable :: cd3(:, :, :, :, :)
      !> Result blocks of one pair; only the position slot adds an index
      real(wp) :: b2(ndim), b3(ndim, ndim), b4(ndim, ndim, ndim)
      !> Active-list indices, active-atom count and nuclear component
      integer :: ia, ib, n, s

      call self%require_deriv(max_deriv, caller)

      n = self%n_active
      call svdw_weights(self, s_1, s_2, s_3)
      call pair_atom_cache(self, max_deriv + 1, ca0, ca1, ca2, ca3, ca4)
      ! the diagonal reads atom A's radius tensors one order above the result
      call rad_atom_cache(self, max_deriv + 1, cd0, cd1, cd2, cd3)

      do ib = 1, n
         do ia = 1, n
            if (ia == ib) then
               call svdw_nucrad_diag_eval(self%param%blend_k, s_1, s_2, s_3, &
                                          self%ps0, self%ps1, self%ps2, &
                                          ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                          ca3(:, :, :, :, ia), &
                                          cd0(:, ia), cd1(:, :, ia), cd2(:, :, :, ia), &
                                          cd3(:, :, :, :, ia), max_deriv, b2, b3, b4)
            else
               call svdw_nucrad_eval(self%param%blend_k, s_1, s_2, s_3, &
                                     self%ps0, self%ps1, self%ps2, &
                                     ca0(:, ia), ca1(:, :, ia), ca2(:, :, :, ia), &
                                     ca3(:, :, :, :, ia), &
                                     cd0(:, ib), cd1(:, :, ib), cd2(:, :, :, ib), &
                                     max_deriv, b2, b3, b4)
            end if
            if (present(res2)) res2(:, ia, ib) = b2
            if (present(res3)) res3(:, :, ia, ib) = b3
            if (present(res4)) then
               do s = 1, ndim
                  call pack_sym_matrix(b4(:, :, s), res4(:, s, ia, ib), ndim)
               end do
            end if
         end do
      end do
   end subroutine svdw_nucrad_block

   !> Pure radius Hessian d^2S / (dR_a dR_b), active-indexed in both atoms
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_radrad d^2S/(dR_a dR_b) [>= n_act, >= n_act]
   subroutine lsf_f2_radrad(self, lsf2_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Radius Hessian
      real(wp), intent(out) :: lsf2_radrad(:, :)

      call svdw_radpair_block(self, 0, "f2_radrad", res2=lsf2_radrad)
   end subroutine lsf_f2_radrad

   !> Mixed third derivative d^3S / (dr dR_a dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_radrad d^3S/(dr dR_a dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_f3_r_radrad(self, lsf3_r_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_radrad(:, :, :)

      call svdw_radpair_block(self, 1, "f3_r_radrad", res3=lsf3_r_radrad)
   end subroutine lsf_f3_r_radrad

   !> Mixed fourth derivative d^4S / (dr^2 dR_a dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_radrad d^4S/(dr^2 dR_a dR_b) [6, >= n_act, >= n_act], packed in the rr pair
   subroutine lsf_f4_rr_radrad(self, lsf4_rr_radrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative, packed in the rr pair
      real(wp), intent(out) :: lsf4_rr_radrad(:, :, :)

      call svdw_radpair_block(self, 2, "f4_rr_radrad", res4=lsf4_rr_radrad)
   end subroutine lsf_f4_rr_radrad

   !> Mixed nuclear-radius Hessian d^2S / (dR_A dR_b)
   !>
   !> @param[in]  self        LSF instance
   !> @param[out] lsf2_rA_rad d^2S/(dR_A dR_b) [3, >= n_act, >= n_act]
   subroutine lsf_f2_rA_rad(self, lsf2_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear-radius Hessian
      real(wp), intent(out) :: lsf2_rA_rad(:, :, :)

      call svdw_nucrad_block(self, 0, "f2_rA_rad", res2=lsf2_rA_rad)
   end subroutine lsf_f2_rA_rad

   !> Mixed third derivative d^3S / (dr dR_A dR_b)
   !>
   !> @param[in]  self          LSF instance
   !> @param[out] lsf3_r_rA_rad d^3S/(dr dR_A dR_b) [3, 3, >= n_act, >= n_act]
   subroutine lsf_f3_r_rA_rad(self, lsf3_r_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed third derivative
      real(wp), intent(out) :: lsf3_r_rA_rad(:, :, :, :)

      call svdw_nucrad_block(self, 1, "f3_r_rA_rad", res3=lsf3_r_rA_rad)
   end subroutine lsf_f3_r_rA_rad

   !> Mixed fourth derivative d^4S / (dr^2 dR_A dR_b)
   !>
   !> @param[in]  self           LSF instance
   !> @param[out] lsf4_rr_rA_rad d^4S/(dr^2 dR_A dR_b) [6, 3, >= n_act, >= n_act], packed in the rr pair
   subroutine lsf_f4_rr_rA_rad(self, lsf4_rr_rA_rad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Mixed fourth derivative, packed in the rr pair
      real(wp), intent(out) :: lsf4_rr_rA_rad(:, :, :, :)

      call svdw_nucrad_block(self, 2, "f4_rr_rA_rad", res4=lsf4_rr_rA_rad)
   end subroutine lsf_f4_rr_rA_rad

   !* ================================================================================= *!
   !*                     Direction-contracted nuclear derivatives                      *!
   !* ================================================================================= *!

   !> Directional nuclear derivative of the level-set value
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . dS/dR_B
   subroutine lsf_tangent_f0(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted value derivative
      real(wp), intent(out) :: res

      !> Blending weights, contracted power sums and unused kernel outputs
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: d1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)

      call self%require_deriv(0, "tangent_f0")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 0, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 0, res, d1, d2, d3)
   end subroutine lsf_tangent_f0

   !> Directional nuclear derivative of the spatial gradient
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of dS/dr [3]
   subroutine lsf_tangent_f1_r(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted gradient derivative
      real(wp), intent(out) :: res(:)

      !> Blending weights, contracted power sums and kernel outputs
      real(wp) :: s_1, s_2, s_3, t0
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: t1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)

      call self%require_deriv(1, "tangent_f1_r")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 1, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 1, t0, t1, d2, d3)
      res = t1
   end subroutine lsf_tangent_f1_r

   !> Directional nuclear derivative of the spatial Hessian
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^2S/dr^2, packed [6]
   subroutine lsf_tangent_f2_rr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted Hessian derivative, packed
      real(wp), intent(out) :: res(:)

      !> Blending weights, contracted power sums and kernel outputs
      real(wp) :: s_1, s_2, s_3, t0
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: t1(ndim), t2(ndim, ndim), d3(ndim, ndim, ndim)

      call self%require_deriv(2, "tangent_f2_rr")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 2, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 2, t0, t1, t2, d3)
      call pack_sym_matrix(t2, res, ndim)
   end subroutine lsf_tangent_f2_rr

   !> Directional nuclear derivative of the third spatial derivative
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d/dR_B of d^3S/dr^3, packed [10]
   subroutine lsf_tangent_f3_rrr(self, v, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted third-derivative derivative, packed
      real(wp), intent(out) :: res(:)

      !> Blending weights, contracted power sums and kernel outputs
      real(wp) :: s_1, s_2, s_3, t0
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      real(wp) :: t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)

      call self%require_deriv(3, "tangent_f3_rrr")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 3, ws0, ws1, ws2, ws3)
      call svdw_tangent_eval(self%param%blend_k, s_1, s_2, s_3, &
                             self%ps0, self%ps1, self%ps2, self%ps3, &
                             ws0, ws1, ws2, ws3, 3, t0, t1, t2, t3)
      call pack_sym_tensor3(t3, res, ndim)
   end subroutine lsf_tangent_f3_rrr

   !> Nuclear Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^2S/(dR_A dR_B) [3, >= active_count()]
   !> @param[in]  vrad Radius directions [ncenters], optional; supplying them
   !>                  promotes the contraction to the joint direction
   !>                  `(v_B, vr_B)`, turning `res` from the nuclear-nuclear
   !>                  block into the nuclear row of the joint Hessian
   subroutine lsf_hvp_f1_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted nuclear Hessian
      real(wp), intent(out) :: res(:, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), d2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(0, "hvp_f1_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 0, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 1, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), 1, aw0, aw1, aw2, aw3)
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), 1, &
                                             awr0, awr1, awr2, awr3)
            aw0 = aw0 + awr0
            aw1 = aw1 + awr1
         end if
         call svdw_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                            at0, at1, at2, at3, aw0, aw1, aw2, aw3, 0, h1, d2, d3)
         res(:, ia) = h1
      end do
   end subroutine lsf_hvp_f1_rA

   !> Directional nuclear derivative of `f2_r_rA`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^3S/(dr dR_A dR_B) [3, 3, >= active_count()]
   !> @param[in]  vrad Radius directions [ncenters], optional; supplying them
   !>                  promotes the contraction to the joint direction
   !>                  `(v_B, vr_B)`, exactly as for [[lsf_hvp_f1_rA]]
   subroutine lsf_hvp_f2_r_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), d3(ndim, ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(1, "hvp_f2_r_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 1, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 2, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), 2, aw0, aw1, aw2, aw3)
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), 2, &
                                             awr0, awr1, awr2, awr3)
            aw0 = aw0 + awr0
            aw1 = aw1 + awr1
            aw2 = aw2 + awr2
         end if
         call svdw_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                            at0, at1, at2, at3, aw0, aw1, aw2, aw3, 1, h1, h2, d3)
         res(:, :, ia) = h2
      end do
   end subroutine lsf_hvp_f2_r_rA

   !> Directional nuclear derivative of `f3_rr_rA`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[out] res  sum_B v_B . d^4S/(dr^2 dR_A dR_B) [6, 3, >= active_count()],
   !>                  packed in the rr pair
   !> @param[in]  vrad Radius directions [ncenters], optional; supplying them
   !>                  promotes the contraction to the joint direction
   !>                  `(v_B, vr_B)`, exactly as for [[lsf_hvp_f1_rA]]
   subroutine lsf_hvp_f3_rr_rA(self, v, res, vrad)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Contracted mixed fourth derivative, packed in the rr pair
      real(wp), intent(out) :: res(:, :, :)
      !> Radius directions
      real(wp), intent(in), optional :: vrad(:)

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom kind tensors and their contracted counterparts
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      real(wp) :: aw0(nkind), aw1(ndim, nkind), aw2(ndim, ndim, nkind)
      real(wp) :: aw3(ndim, ndim, ndim, nkind)
      real(wp) :: awr0(nkind), awr1(ndim, nkind), awr2(ndim, ndim, nkind)
      real(wp) :: awr3(ndim, ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1(ndim), h2(ndim, ndim), h3(ndim, ndim, ndim)
      !> Active-list index and nuclear component
      integer :: ia, s

      call self%require_deriv(2, "hvp_f3_rr_rA")

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, 2, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 3, at0, at1, at2, at3, at4)
         call atom_tangent_tensors(self, ia, v(:, self%act_atom(ia)), 3, aw0, aw1, aw2, aw3)
         if (present(vrad)) then
            call atom_radius_tangent_tensors(self, ia, vrad(self%act_atom(ia)), 3, &
                                             awr0, awr1, awr2, awr3)
            aw0 = aw0 + awr0
            aw1 = aw1 + awr1
            aw2 = aw2 + awr2
            aw3 = aw3 + awr3
         end if
         call svdw_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                            at0, at1, at2, at3, aw0, aw1, aw2, aw3, 2, h1, h2, h3)
         do s = 1, ndim
            call pack_sym_matrix(h3(:, :, s), res(:, s, ia), ndim)
         end do
      end do
   end subroutine lsf_hvp_f3_rr_rA

   !* ================================================================================= *!
   !*                    Radius row of the joint Hessian-vector product                 *!
   !* ================================================================================= *!

   !> Radius row of the joint Hessian-vector product, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  sum_B (v_B . d/dR_B + vr_B d/dR_b) dS/dR_a
   !>                  [>= active_count()]
   subroutine lsf_hvp_f1_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted radius Hessian row
      real(wp), intent(out) :: res(:)

      call svdw_radius_hvp(self, v, vrad, 0, "hvp_f1_rad", res1=res)
   end subroutine lsf_hvp_f1_rad

   !> Joint directional derivative of `f2_r_rad`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed third derivative [3, >= active_count()]
   subroutine lsf_hvp_f2_r_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed third derivative
      real(wp), intent(out) :: res(:, :)

      call svdw_radius_hvp(self, v, vrad, 1, "hvp_f2_r_rad", res2=res)
   end subroutine lsf_hvp_f2_r_rad

   !> Joint directional derivative of `f3_rr_rad`, active-indexed
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  v    Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad Radius directions [ncenters]
   !> @param[out] res  Contracted mixed fourth derivative [6, >= active_count()], packed in the rr pair
   subroutine lsf_hvp_f3_rr_rad(self, v, vrad, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Radius directions
      real(wp), intent(in) :: vrad(:)
      !> Contracted mixed fourth derivative, packed in the rr pair
      real(wp), intent(out) :: res(:, :)

      call svdw_radius_hvp(self, v, vrad, 2, "hvp_f3_rr_rad", res3=res)
   end subroutine lsf_hvp_f3_rr_rad

   !> Shared driver of the three radius-row Hessian-vector products
   !>
   !> The three public entry points differ only in `max_deriv` and in which output
   !> they keep, and unlike the `_rA` ladder there is no reason to spell the loop
   !> out three times, the radius channel having one code path
   !>
   !> - exactly one of `res1`/`res2`/`res3` must be present, matching `max_deriv`
   !>
   !> @param[in]  self      LSF instance
   !> @param[in]  v         Nuclear displacement directions [3, ncenters]
   !> @param[in]  vrad      Radius directions [ncenters]
   !> @param[in]  max_deriv Highest spatial-derivative order (0..2)
   !> @param[in]  caller    Name used by `require_deriv` on failure
   !> @param[out] res1      Order-0 result [>= active_count()] (optional)
   !> @param[out] res2      Order-1 result [3, >= active_count()] (optional)
   !> @param[out] res3      Order-2 result [6, >= active_count()], packed (optional)
   subroutine svdw_radius_hvp(self, v, vrad, max_deriv, caller, res1, res2, res3)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
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

      !> Blending weights and contracted power sums
      real(wp) :: s_1, s_2, s_3
      real(wp) :: ws0(nkind), ws1(ndim, nkind), ws2(ndim, ndim, nkind)
      real(wp) :: ws3(ndim, ndim, ndim, nkind)
      !> Per-atom radius tensors and the radius derivative of the joint share
      real(wp) :: atd0(nkind), atd1(ndim, nkind), atd2(ndim, ndim, nkind)
      real(wp) :: atd3(ndim, ndim, ndim, nkind)
      real(wp) :: awd0(nkind), awd1(ndim, nkind), awd2(ndim, ndim, nkind)
      !> Kernel outputs of one atom
      real(wp) :: h1, h2(ndim), h3(ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(max_deriv, caller)

      call svdw_weights(self, s_1, s_2, s_3)
      call tangent_powersums(self, v, max_deriv, ws0, ws1, ws2, ws3, vrad)
      do ia = 1, self%n_active
         ! `d/dR_a` adds no index, so unlike the `_rA` ladder the per-atom
         ! tensors are needed only to `max_deriv`, not `max_deriv + 1`; `awd`
         ! is `d/dR_a` of this atom's whole share `aw + awr` of `ws`
         call atom_rad_tensors(self, ia, max_deriv, atd0, atd1, atd2, atd3)
         call atom_tangent_rad_tensors(self, ia, v(:, self%act_atom(ia)), &
                                       vrad(self%act_atom(ia)), max_deriv, &
                                       awd0, awd1, awd2)
         ! `atd3` exists only because the tensor builder takes it as a
         ! non-optional `intent(out)` dummy; the eval below never sees it, so it
         ! is left undefined on purpose
         call svdw_radius_hvp_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, ws0, ws1, ws2, &
                                   atd0, atd1, atd2, awd0, awd1, awd2, max_deriv, &
                                   h1, h2, h3)
         if (present(res1)) res1(ia) = h1
         if (present(res2)) res2(:, ia) = h2
         if (present(res3)) call pack_sym_matrix(h3, res3(:, ia), ndim)
      end do
   end subroutine svdw_radius_hvp

   !* ================================================================================= *!
   !*                        Jet-contracted nuclear derivative                          *!
   !* ================================================================================= *!

   !> Adjoint-weighted nuclear gradient, active-indexed
   !>
   !> Contracts one evaluation point's adjoint weights against the nuclear
   !> ladder and keeps only the nuclear index:
   !>
   !>    res(s, i) = w0*lsf1_rA(s, i) + sum_a w1(a)*lsf2_r_rA(a, s, i)
   !>                + sum_a sum_b w2(a, b)*lsf3_rr_rA(a, b, s, i)
   !>
   !> The `f1_rA` rung with the jet indices contracted away, exactly as
   !> [[lsf_hvp_f1_rA]] is that same rung contracted with a nuclear direction --
   !> the reverse-mode mirror of the `tangent_*` family, which contracts the
   !> nuclear index and keeps the spatial ones
   !>
   !> `w2` is the packed symmetric seed [6], off-diagonals not doubled; the kernel
   !> sees it unpacked to a symmetric 3x3, so each off-diagonal slot counts twice
   !>
   !> The contraction happens inside the kernel, on the weighted jet rather than
   !> on its 3 + 9 + 27 components, so an adjoint pass never has to materialize
   !> the `(3, 3, 3, >= active_count())` tensor [[lsf_f3_rr_rA]] would hand it
   !>
   !> - the prepared order is the same as for [[lsf_f3_rr_rA]], the same
   !>   `lsf3_rr_rA` rung entering the contraction
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the level-set value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Packed symmetric adjoint weights of the spatial Hessian [6]
   !> @param[out] res  Jet-contracted nuclear gradient [3, >= active_count()]
   subroutine lsf_vjp_f1_rA(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Packed symmetric adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(:)
      !> Jet-contracted nuclear gradient
      real(wp), intent(out) :: res(:, :)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom kind tensors
      real(wp) :: at0(nkind), at1(ndim, nkind), at2(ndim, ndim, nkind)
      real(wp) :: at3(ndim, ndim, ndim, nkind), at4(ndim, ndim, ndim, ndim, nkind)
      !> Kernel output of one atom
      real(wp) :: vjp_f1_rA(ndim)
      !> Symmetric seed unpacked from `w2`
      real(wp) :: wsym(ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(2, "vjp_f1_rA")
      call unpack_sym_matrix(w2, wsym, ndim)

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         call atom_tensors(self, ia, 3, at0, at1, at2, at3, at4)
         call svdw_vjp_eval(self%param%blend_k, s_1, s_2, s_3, &
                            self%ps0, self%ps1, self%ps2, &
                            at0, at1, at2, at3, w0, w1, wsym, vjp_f1_rA)
         res(:, ia) = vjp_f1_rA
      end do
   end subroutine lsf_vjp_f1_rA

   !> Adjoint-weighted radius gradient, active-indexed
   !>
   !> The radius twin of [[lsf_vjp_f1_rA]]: the same evaluation point's adjoint
   !> weights, contracted against the radius ladder instead of the nuclear one,
   !>
   !>    res(i) = w0*lsf1_rad(i) + sum_a w1(a)*lsf2_r_rad(a, i)
   !>             + sum_a sum_b w2(a, b)*lsf3_rr_rad(a, b, i)
   !>
   !> A radius is a scalar, so -- exactly as [[lsf_f3_rr_rad]] carries one rank
   !> less than [[lsf_f3_rr_rA]] -- there is no index left once the jet indices
   !> are contracted away: one number per active atom, against the three of
   !> [[lsf_vjp_f1_rA]] and the 1 + 3 + 9 of [[lsf_f3_rr_rad]]
   !>
   !> `w2` is the packed symmetric seed [6], off-diagonals not doubled; the kernel
   !> sees it unpacked to a symmetric 3x3, so each off-diagonal slot counts twice
   !>
   !> The contraction happens inside the kernel, on the weighted jet rather than
   !> on its components
   !>
   !> - the prepared order is the same as for [[lsf_f3_rr_rad]], since `d/dR_a`
   !>   consumes no derivative order, so the kernel reads the power sums and
   !>   per-atom tensors only up to order 2
   !>
   !> @param[in]  self LSF instance
   !> @param[in]  w0   Adjoint weight of the level-set value
   !> @param[in]  w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]  w2   Packed symmetric adjoint weights of the spatial Hessian [6]
   !> @param[out] res  Jet-contracted radius gradient [>= active_count()]
   subroutine lsf_vjp_f1_rad(self, w0, w1, w2, res)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Packed symmetric adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(:)
      !> Jet-contracted radius gradient
      real(wp), intent(out) :: res(:)

      !> Blending weights
      real(wp) :: s_1, s_2, s_3
      !> Per-atom radius tensors
      real(wp) :: atd0(nkind), atd1(ndim, nkind), atd2(ndim, ndim, nkind)
      real(wp) :: atd3(ndim, ndim, ndim, nkind)
      !> Kernel output of one atom
      real(wp) :: vjp_f1_rad
      !> Symmetric seed unpacked from `w2`
      real(wp) :: wsym(ndim, ndim)
      !> Active-list index
      integer :: ia

      call self%require_deriv(2, "vjp_f1_rad")
      call unpack_sym_matrix(w2, wsym, ndim)

      call svdw_weights(self, s_1, s_2, s_3)
      do ia = 1, self%n_active
         ! Order 2 is enough: the radius derivative bumps no index, so the kernel
         ! never reaches `atd3` and does not take it as an argument -- unlike
         ! `svdw_radius_eval`, which does and therefore needs `atd3` defined
         call atom_rad_tensors(self, ia, 2, atd0, atd1, atd2, atd3)
         call svdw_radius_vjp_eval(self%param%blend_k, s_1, s_2, s_3, &
                                   self%ps0, self%ps1, self%ps2, &
                                   atd0, atd1, atd2, w0, w1, wsym, vjp_f1_rad)
         res(ia) = vjp_f1_rad
      end do
   end subroutine lsf_vjp_f1_rad

   !* ================================================================================= *!
   !*                                Screening offset                                   *!
   !* ================================================================================= *!

   !> Atom reach of the SvdW model beyond the atom surface
   !>
   !> Fade on: the band end `X_on + Delta = 3(fade_start + fade_width)/k`
   !>
   !> - `screening_threshold` is not read; the reach is part of the model
   !> - the same `x_on + delta` stage 1 forms as `band_end`, so gate and fade
   !>   agree on one number: beyond it `w = S_4(1 - t)` has the factor
   !>   `clamp(1 - t)**5 = 0` exactly, and the gate's squared compare masks the
   !>   rest
   !> - the gate (`r2 <= (R + X)**2`) and the fade (`sqrt(r2) - R >= X`) can
   !>   disagree in a rounding sliver at the edge, where the term is a gate zero
   !>   or at most ~1e-75 `u`: numerically negligible
   !> - band invalid (validation, `blend_k <= 0`): `huge`, screen nothing;
   !>   [[check_fade]] refuses the point before any sum
   !>
   !> Fade off: the offset where `u = exp(-(k/3)(x - R))` equals the threshold,
   !>
   !>    `x - R = -3*ln(threshold)/k`
   !>
   !> - radius-independent; `huge` with `threshold <= 0`, `k <= 0` or a
   !>   negative radius (screening disabled)
   !>
   !> The *screening* bound query, whose companion `exclusion_radius` below is the
   !> *surface* bound query; see the module header of
   !> [[moist_cavity_drop_lsf_base]] for why the two cannot be expressed in terms
   !> of each other
   !>
   !> @param[in] self    SvdW LSF instance (reads param, screening_threshold)
   !> @param[in] radius  Atom radius (Bohr, unused: the SvdW offset is radius-independent)
   !> @returns   offset  Radial offset from the atom surface (Bohr)
   pure function lsf_screening_offset(self, radius) result(offset)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Atom radius (Bohr)
      real(wp), intent(in) :: radius
      !> Radial offset from atom surface (Bohr)
      real(wp) :: offset

      !> Blending sharpness and screening threshold
      real(wp) :: k_local, threshold
      !> Fade band in `d_A`
      real(wp) :: x_on, delta

      k_local = self%param%blend_k
      threshold = self%screening_threshold

      if (self%param%fade_enabled) then
         offset = huge(0.0_wp)
         if (.not. fade_band_valid(self)) return
         if (radius < 0.0_wp) return
         call fade_band(self, x_on, delta)
         offset = x_on + delta
         return
      end if

      ! With threshold <= 0, k <= 0 or a nonsensical radius the criterion is
      ! ill-defined; report an unbounded reach so screening is effectively disabled
      if (threshold <= 0.0_wp .or. k_local <= 0.0_wp .or. radius < 0.0_wp) then
         offset = huge(0.0_wp)
         return
      end if

      offset = max(0.0_wp, -3.0_wp*log(threshold)/k_local)
   end function lsf_screening_offset

   !* ================================================================================= *!
   !*                        Surface-free exclusion certificate                         *!
   !* ================================================================================= *!

   !> Surface-free radius around a point, from its value `S`
   !>
   !> Fade off: `r = |S|`, the 1-Lipschitz property of the log-sum-exp soft minimum
   !>
   !> - `grad S` is the blend-weighted mean of unit vectors `(r - r_I)/||r - r_I||`,
   !>   so `||grad S|| <= 1` and `S` cannot reach zero within `|S(x)|`
   !> - needs non-negative blend coefficients; a negative one turns the mean into
   !>   an extrapolation and certifies nothing (zero)
   !>
   !> Fade on: `r = |tau(S)|`, `tau(phi) = int_0^phi ds/(1 + K(s))`
   !>
   !> - `Z` is homogeneous of degree 3 in `g_A = w_A u_A`, so inside a band
   !>   `||grad S|| <= 1 + K(S)`, `K = A eta**3 + B eta**(3/2) + G eta`,
   !>   `eta = exp(k S/3 - a)`; `T = tau(S)` is exactly 1-Lipschitz
   !> - closed form per point, constants from [[certificate_setup]]:
   !>   - `S <= 0`: `r = -S/(1 + K0)`
   !>   - `S > 0`: `phi = min(S, s_cap)`, `gam = A eta(phi)**2 + B sqrt(eta(phi)) + G`,
   !>     `r = phi - (3/k) ln((1 + gam eta(phi))/(1 + gam eps))`
   !> - equals `|S|` to ~1e-4 Bohr up to S ~ 6 Bohr at the defaults, capped at
   !>   `r_cap` (~8-10 Bohr); `S = +Inf` gives `r_cap`
   !> - zero (certify nothing) outside the validity guard (`s >= 0`, `s1 > 0`,
   !>   `s3 < 24 s1`) and while the constants are stale (see
   !>   [[svdw_certificate_type]])
   !>
   !> A NaN value or `-Inf` certifies nothing in both modes; `+Inf` (`PD`/`Z` zero, the
   !> value of an unreached batch column without a finite `S`) nothing with the fade off
   !>
   !> @param[in] self LSF instance
   !> @param[in] lsf0 LSF value at the evaluation point
   !> @returns   r    Radius of a ball around the point free of surface
   pure function lsf_svdw_exclusion_radius(self, lsf0) result(r)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> LSF value at the evaluation point
      real(wp), intent(in) :: lsf0
      !> Surface-free radius
      real(wp) :: r

      !> Capped level, its `eta` and the slope factor `Gamma`
      real(wp) :: phi, eta, gam

      r = 0.0_wp
      if (ieee_is_nan(lsf0)) return
      if (lsf0 < -huge(lsf0)) return
      if (self%param%blend_1b < 0.0_wp) return
      if (self%param%blend_2b < 0.0_wp) return
      if (self%param%blend_3b < 0.0_wp) return
      if (.not. self%param%fade_enabled) then
         if (lsf0 > huge(lsf0)) return
         r = abs(lsf0)
         return
      end if

      if (.not. certificate_current(self)) return
      if (.not. self%cert%valid) return
      associate (c => self%cert)
         if (lsf0 <= 0.0_wp) then
            r = -lsf0/(1.0_wp + c%k0)
         else
            phi = min(lsf0, c%s_cap)
            eta = exp(c%k*phi/3.0_wp - c%a)
            gam = c%coef_a*eta*eta + c%coef_b*sqrt(eta) + c%coef_g
            r = phi - (3.0_wp/c%k)*log((1.0_wp + gam*eta)/(1.0_wp + gam*c%eps))
         end if
      end associate
      r = max(r, 0.0_wp)
   end function lsf_svdw_exclusion_radius

   !> Surface-free radius around every unreached point, fade on
   !>
   !> `X_on + Delta - d*`, `d* = (3/k) ln n + ln(Zmax)/k`, `n = ncenters`
   !>
   !> - every surface point `y` has an atom with `d_A(y) <= d*` (`Z(y) = 1 <=
   !>   Zmax e1**3`, `e1 <= n max u_A`); a point with no atom within the band
   !>   end has every atom surface at least the reach `X_on + Delta` away, and
   !>   the nearest-surface distance is 1-Lipschitz
   !> - unreached means no atom within the band end OR a non-finite lift with
   !>   `S > 0`; the certificate covers only the first, which `blend_1b > 0`
   !>   guarantees at order 0: the one-body term `s1 sum g_A**3` keeps `Z > 0` for
   !>   any atom in reach (`u_A >= exp(-(a + b))`; underflow needs `w_A < 1e-92`,
   !>   the last ~1e-19 of the band)
   !> - valid then for every unreached column of an order-0 batch and every point
   !>   with an empty active list on the per-point path
   !> - remark: the same argument gives `r >= nearest - d*` for any point whose
   !>   nearest-atom surface distance `nearest` is known; not exposed
   !> - zero with the fade off (no reach, no such bound), with a negative blend
   !>   coefficient, at `blend_1b = 0` (an unreached column can have one or two
   !>   atoms in reach, `Z = s3 e3 = 0`), while the constants are stale, or when
   !>   `d*` exceeds the band end
   !>
   !> @param[in] self LSF instance
   !> @returns   r    Surface-free radius around an unreached point
   pure function lsf_svdw_far_exclusion_radius(self) result(r)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Surface-free radius
      real(wp) :: r

      r = 0.0_wp
      if (.not. self%param%fade_enabled) return
      if (.not. certificate_current(self)) return
      if (.not. self%cert%far_valid) return
      r = self%cert%far_radius
   end function lsf_svdw_far_exclusion_radius

   !> Whether the certificate constants belong to the live parameters and atom count
   !>
   !> @param[in] self    LSF instance
   !> @returns   current Constants match `param` and `ncenters`
   pure function certificate_current(self) result(current)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Constants match `param` and `ncenters`
      logical :: current

      current = self%cert%n == self%ncenters .and. &
                all(self%cert%key == [self%param%blend_k, self%param%blend_1b, &
                                      self%param%blend_2b, self%param%blend_3b, &
                                      self%param%fade_start, self%param%fade_width])
   end function certificate_current

   !> Setup-time constants of the faded exclusion certificates
   !>
   !> The constants of the faded certificate bounds, once per `update`:
   !>
   !>    mu1  = max_t S'(t) e^(-b t)
   !>    mu32 = max_t S'(t) e^(-3bt/2) w(t)^(1/2)
   !>    mu3  = max_t S'(t) e^(-3bt) w(t)^2              ([[fade_slope_max]])
   !>    zmin = s3/6 - s3^2/(16 (s1 + s3/3)),  C3 = 1/zmin,  Cq = 1/min(s1, s2/2)
   !>    A = (n/(3b)) 3 s1 mu3,  B = (n/(3b)) 1.5 s2 mu32 sqrt(min(C3, Cq)),
   !>    G = (n/(3b)) 0.5 s3 mu1 C3^(2/3),  K0 = A eps^3 + B eps^(3/2) + G eps
   !>    Zmax = max(s1, s3/6) if s1 - s2/2 + s3/3 >= 0, else s1 + s2/2 + s3/6
   !>    d* = (3/k) ln n + ln(Zmax)/k,  far radius = X_on + Delta - d*
   !>    s_cap = argmax tau_lb,  r_cap = tau_lb(s_cap)    ([[certificate_cap]])
   !>
   !> - `B = 0` when `s2 = 0`, `G = 0` and `C3` unused when `s3 = 0`
   !> - fade off, invalid band, non-finite or negative weights, `s1 = 0`, no
   !>   atoms: both certificates off
   !>
   !> @param[inout] self LSF instance after `lsf_base_update`
   subroutine certificate_setup(self)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self

      !> Blending weights, sharpness and band in e-folds
      real(wp) :: s1, s2, s3, k, a, b
      !> Fade band in `d_A`
      real(wp) :: x_on, delta
      !> Fade slope maxima
      real(wp) :: mu1, mu32, mu3
      !> Power-sum control constants and the slope prefactor `n/(3b)`
      real(wp) :: zmin, c3, cq, pre

      self%cert = svdw_certificate_type()
      s1 = self%param%blend_1b
      s2 = self%param%blend_2b
      s3 = self%param%blend_3b
      k = self%param%blend_k
      a = self%param%fade_start
      b = self%param%fade_width
      self%cert%key = [k, s1, s2, s3, a, b]
      self%cert%n = self%ncenters

      if (.not. self%param%fade_enabled) return
      if (.not. fade_band_valid(self)) return
      if (self%ncenters <= 0) return
      if (.not. all(ieee_is_finite([s1, s2, s3]))) return
      if (min(s1, s2, s3) < 0.0_wp) return

      call fade_band(self, x_on, delta)
      self%cert%k = k
      self%cert%a = a
      self%cert%eps = exp(-a)

      ! Both need s1 > 0: reached, no level-only bound without it; far, with
      ! s1 = 0 an unreached column (Z = s3 e3 = 0) can have atoms in reach
      if (.not. s1 > 0.0_wp) return

      ! Far certificate, Zmax >= s1 > 0
      if (s1 - 0.5_wp*s2 + s3/3.0_wp >= 0.0_wp) then
         self%cert%zmax = max(s1, s3/6.0_wp)
      else
         self%cert%zmax = s1 + 0.5_wp*s2 + s3/6.0_wp
      end if
      self%cert%dstar = (3.0_wp/k)*log(real(self%ncenters, wp)) + log(self%cert%zmax)/k
      self%cert%far_radius = max(0.0_wp, x_on + delta - self%cert%dstar)
      self%cert%far_valid = .true.

      ! Reached certificate: s3 < 24 s1
      if (s3 > 0.0_wp .and. .not. s3 < 24.0_wp*s1) return

      mu1 = fade_slope_max(b, 1)
      mu32 = fade_slope_max(b, 2)
      mu3 = fade_slope_max(b, 3)
      c3 = huge(1.0_wp)
      if (s3 > 0.0_wp) then
         zmin = s3/6.0_wp - s3*s3/(16.0_wp*(s1 + s3/3.0_wp))
         c3 = 1.0_wp/zmin
      end if
      cq = huge(1.0_wp)
      if (s2 > 0.0_wp) cq = 1.0_wp/min(s1, 0.5_wp*s2)
      pre = real(self%ncenters, wp)/(3.0_wp*b)

      associate (c => self%cert)
         c%coef_a = pre*3.0_wp*s1*mu3
         c%coef_b = 0.0_wp
         if (s2 > 0.0_wp) c%coef_b = pre*1.5_wp*s2*mu32*sqrt(min(c3, cq))
         c%coef_g = 0.0_wp
         if (s3 > 0.0_wp) c%coef_g = pre*0.5_wp*s3*mu1*c3**(2.0_wp/3.0_wp)
         c%k0 = c%coef_a*c%eps**3 + c%coef_b*c%eps**1.5_wp + c%coef_g*c%eps
         call certificate_cap(c)
         c%valid = .true.
      end associate
   end subroutine certificate_setup

   !> Closed-form lower bound `tau_lb(phi)` of `tau(phi)`, `phi >= 0`
   !>
   !> @param[in] c   Certificate constants
   !> @param[in] phi Level (Bohr), >= 0
   !> @returns   t   `phi - (3/k) ln((1 + gam eta)/(1 + gam eps))`
   pure function tau_lower(c, phi) result(t)
      !> Certificate constants
      type(svdw_certificate_type), intent(in) :: c
      !> Level
      real(wp), intent(in) :: phi
      !> Lower bound of tau
      real(wp) :: t

      !> `eta(phi)` and `Gamma(phi)`
      real(wp) :: eta, gam

      eta = exp(c%k*phi/3.0_wp - c%a)
      gam = c%coef_a*eta*eta + c%coef_b*sqrt(eta) + c%coef_g
      t = phi - (3.0_wp/c%k)*log((1.0_wp + gam*eta)/(1.0_wp + gam*c%eps))
   end function tau_lower

   !> Cap level `s_cap` (argmax of `tau_lb` over `phi >= 0`) and `r_cap = tau_lb(s_cap)`
   !>
   !> - tightness only: `tau` increases and `tau_lb <= tau` at every `phi`, so
   !>   `tau_lb(min(S, s_cap))` is a valid radius for any `s_cap >= 0`
   !> - search over `phi = (3/k)(a + L)`, `L = ln eta` from `-a` up to where one
   !>   term of `K` reaches 1e8 (`tau` has stopped growing there): 2048-point grid,
   !>   then golden section on the two cells around the grid maximum
   !>
   !> @param[inout] c Certificate constants with A, B, G, k, a, eps set
   pure subroutine certificate_cap(c)
      !> Certificate constants
      type(svdw_certificate_type), intent(inout) :: c

      !> Grid size and golden-section steps
      integer, parameter :: ngrid = 2048, ngold = 80
      !> `K` level that ends the search
      real(wp), parameter :: k_end = 1.0e8_wp
      !> Inverse golden ratio
      real(wp), parameter :: ginv = 0.6180339887498949_wp
      !> Upper `ln eta`, search range in `phi`, grid step
      real(wp) :: l_hi, phi_hi, h
      !> Golden-section bracket and probes
      real(wp) :: lo, hi, x1, x2, f1, f2, t, best, phi_best
      integer :: i

      l_hi = log(k_end) - log(c%coef_a)/3.0_wp
      if (c%coef_b > 0.0_wp) l_hi = min(l_hi, 2.0_wp*(log(k_end) - log(c%coef_b)))
      if (c%coef_g > 0.0_wp) l_hi = min(l_hi, log(k_end) - log(c%coef_g))
      l_hi = min(l_hi, 200.0_wp)
      phi_hi = max(0.0_wp, (3.0_wp/c%k)*(c%a + l_hi))

      h = phi_hi/real(ngrid, wp)
      phi_best = 0.0_wp
      best = tau_lower(c, 0.0_wp)
      do i = 1, ngrid
         t = tau_lower(c, h*real(i, wp))
         if (t > best) then
            best = t
            phi_best = h*real(i, wp)
         end if
      end do

      lo = max(0.0_wp, phi_best - h)
      hi = min(phi_hi, phi_best + h)
      x1 = hi - ginv*(hi - lo)
      x2 = lo + ginv*(hi - lo)
      f1 = tau_lower(c, x1)
      f2 = tau_lower(c, x2)
      do i = 1, ngold
         if (f1 < f2) then
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + ginv*(hi - lo)
            f2 = tau_lower(c, x2)
         else
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - ginv*(hi - lo)
            f1 = tau_lower(c, x1)
         end if
      end do
      if (max(f1, f2) > best) then
         best = max(f1, f2)
         phi_best = merge(x1, x2, f1 >= f2)
      end if

      c%s_cap = phi_best
      c%r_cap = max(0.0_wp, best)
   end subroutine certificate_cap

   !> Fade slope maximum `mu_j = max_t S'(t) e^(-m b t) w(t)^p` on `[0, 1]`
   !>
   !> - `kind_j = 1`: `m = 1`, `p = 0` (`mu1`); `2`: `m = 3/2`, `p = 1/2` (`mu32`);
   !>   `3`: `m = 3`, `p = 2` (`mu3`)
   !> - `S = S_n`, `n` = [[fade_order]], `S'(t) = c t^n (1-t)^n`,
   !>   `w(t) = S_n(1 - t)`; every factor is log-concave, so the product is
   !>   unimodal and golden section on the bracket of a 1024-cell grid finds the
   !>   maximum; the result is raised by 1e-10 relative against rounding
   !> - grid endpoints skipped (`S' = 0` there), so no integer power of zero
   !>
   !> @param[in] b      Fade width in e-folds
   !> @param[in] kind_j Which maximum (1, 2, 3)
   !> @returns   mu     Maximum
   pure function fade_slope_max(b, kind_j) result(mu)
      !> Fade width in e-folds
      real(wp), intent(in) :: b
      !> Which maximum
      integer, intent(in) :: kind_j
      !> Maximum
      real(wp) :: mu

      !> Grid size and golden-section steps
      integer, parameter :: ngrid = 1024, ngold = 80
      !> Inverse golden ratio
      real(wp), parameter :: ginv = 0.6180339887498949_wp
      !> Golden-section bracket and probes
      real(wp) :: lo, hi, x1, x2, f1, f2, t, best, t_best
      integer :: i

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

      !> `S'(t) e^(-m b t) w(t)^p` of this maximum
      !>
      !> @param[in] tt Band coordinate in [0, 1]
      !> @returns   f  Term value
      pure function slope_term(tt) result(f)
         !> Band coordinate
         real(wp), intent(in) :: tt
         !> Term value
         real(wp) :: f

         !> Smoothstep slope and fade weight
         real(wp) :: ds, w

         call smoothstep_slope_weight(tt, ds, w)
         select case (kind_j)
         case (1)
            f = ds*exp(-b*tt)
         case (2)
            f = ds*exp(-1.5_wp*b*tt)*sqrt(w)
         case default
            f = ds*exp(-3.0_wp*b*tt)*w*w
         end select
      end function slope_term

   end function fade_slope_max

   !> Slope `S_n'(t)` and fade weight `w = S_n(1 - t)` of the order-[[fade_order]] smoothstep
   !>
   !> - `S_n'(t) = (2n+1)!/(n!)^2 t^n (1-t)^n`
   !> - `w = (1-t)^(n+1) sum_k C(n+k, k) t^k`, the kernel's mirrored Horner form
   !>
   !> @param[in]  t  Band coordinate in (0, 1); `t = 0` or `1` traps in arm64 `__powidf2`
   !> @param[out] ds Slope `S_n'(t)`
   !> @param[out] w  Fade weight `1 - S_n(t)`
   pure subroutine smoothstep_slope_weight(t, ds, w)
      !> Band coordinate
      real(wp), intent(in) :: t
      !> Slope
      real(wp), intent(out) :: ds
      !> Fade weight
      real(wp), intent(out) :: w

      !> Binomial coefficient, slope prefactor and Horner accumulator
      real(wp) :: binom, c, acc
      integer :: j, m

      c = 1.0_wp
      do j = fade_order + 1, 2*fade_order + 1
         c = c*real(j, wp)
      end do
      do j = 1, fade_order
         c = c/real(j, wp)
      end do
      ds = c*(t*(1.0_wp - t))**fade_order

      acc = 0.0_wp
      do j = fade_order, 0, -1
         binom = 1.0_wp
         do m = 1, j
            binom = binom*real(fade_order + m, wp)/real(m, wp)
         end do
         acc = acc*t + binom
      end do
      w = (1.0_wp - t)**(fade_order + 1)*acc
   end subroutine smoothstep_slope_weight

   !> Exclude stationary points by a common sign of atomic cross products
   !>
   !> For nonnegative blend coefficients, Z has non-negative coefficients in
   !> `g_A = w_A u_A`, each `g_A` non-increasing in `d_A`, so grad S is a
   !> nonnegative combination of atom normals
   !>
   !> - coefficients `Z_A ((k/3) g_A + |w_A'| u_A)/(k Z) >= 0`; needs a monotone
   !>   fade (`w' <= 0`), which every polynomial smoothstep `S_n` gives; fade off
   !>   is `w = 1`
   !> - stationarity needs `(x - anchor) x grad S = 0`; each atom gives
   !>   `(x - anchor) x (anchor - atom)` over a positive distance, and a component
   !>   with one strict sign for every atom throughout the box cannot cancel
   !> - a cached enclosure of all atom centers bounds the whole positive cone in
   !>   constant work per box, covering every possible subset of atoms in reach
   !> - outward-rounded interval arithmetic includes coordinate rounding
   !> - an inconclusive enclosure keeps the box for subdivision or refinement
   !>
   !> @param[in] self       Positive-mixture SvdW model
   !> @param[in] anchor     Distance origin (3)
   !> @param[in] center     Box midpoint (3)
   !> @param[in] half_width Nonnegative cube half-width
   pure function lsf_svdw_stationary_box_excluded(self, anchor, center, half_width) result(excluded)
      !> LSF instance
      class(moist_cavity_drop_lsf_svdw_type), intent(in) :: self
      !> Distance origin
      real(wp), intent(in) :: anchor(3)
      !> Box midpoint
      real(wp), intent(in) :: center(3)
      !> Cube half-width
      real(wp), intent(in) :: half_width
      !> Whole box contains no stationary distance point
      logical :: excluded
      !> Outward enclosures of x-anchor and anchor-atom
      real(wp) :: ulo(3), uhi(3), vlo(3), vhi(3)
      !> Product and cross-product enclosures
      real(wp) :: plo, phi, qlo, qhi, cross_lo, cross_hi
      !> Cyclic component indices
      integer :: k, j, l

      excluded = .false.
      if (self%ncenters <= 0) return
      if (.not. self%center_bounds_valid) return
      if (.not. all(ieee_is_finite([self%param%blend_k, self%param%blend_1b, &
                                    self%param%blend_2b, self%param%blend_3b]))) return
      if (self%param%blend_k <= 0.0_wp) return
      if (min(self%param%blend_1b, self%param%blend_2b, self%param%blend_3b) < 0.0_wp) return
      if (max(self%param%blend_1b, self%param%blend_2b, self%param%blend_3b) <= 0.0_wp) return
      if (.not. all(ieee_is_finite(anchor))) return
      if (.not. all(ieee_is_finite(center))) return
      if (.not. ieee_is_finite(half_width)) return
      if (half_width < 0.0_wp) return

      ulo = nearest(nearest(center - half_width, -1.0_wp) - anchor, -1.0_wp)
      uhi = nearest(nearest(center + half_width, 1.0_wp) - anchor, 1.0_wp)
      if (.not. all(ieee_is_finite(ulo))) return
      if (.not. all(ieee_is_finite(uhi))) return
      if (all(ulo <= 0.0_wp .and. uhi >= 0.0_wp)) return
      vlo = nearest(anchor - self%center_upper, -1.0_wp)
      vhi = nearest(anchor - self%center_lower, 1.0_wp)
      if (.not. all(ieee_is_finite(vlo))) return
      if (.not. all(ieee_is_finite(vhi))) return
      do k = 1, 3
         j = modulo(k, 3) + 1
         l = modulo(k + 1, 3) + 1
         call interval_product(ulo(j), uhi(j), vlo(l), vhi(l), plo, phi)
         call interval_product(ulo(l), uhi(l), vlo(j), vhi(j), qlo, qhi)
         cross_lo = nearest(plo - qhi, -1.0_wp)
         cross_hi = nearest(phi - qlo, 1.0_wp)
         if (.not. ieee_is_finite(cross_lo)) return
         if (.not. ieee_is_finite(cross_hi)) return
         if (cross_lo > 0.0_wp .or. cross_hi < 0.0_wp) then
            excluded = .true.
            return
         end if
      end do
   end function lsf_svdw_stationary_box_excluded

   !> Outward-rounded product of two finite real intervals
   !>
   !> @param[in]  alo First lower endpoint
   !> @param[in]  ahi First upper endpoint
   !> @param[in]  blo Second lower endpoint
   !> @param[in]  bhi Second upper endpoint
   !> @param[out] lo  Product lower endpoint
   !> @param[out] hi  Product upper endpoint
   pure subroutine interval_product(alo, ahi, blo, bhi, lo, hi)
      !> First lower endpoint
      real(wp), intent(in) :: alo
      !> First upper endpoint
      real(wp), intent(in) :: ahi
      !> Second lower endpoint
      real(wp), intent(in) :: blo
      !> Second upper endpoint
      real(wp), intent(in) :: bhi
      !> Product lower endpoint
      real(wp), intent(out) :: lo
      !> Product upper endpoint
      real(wp), intent(out) :: hi
      !> Endpoint products
      real(wp) :: products(4)

      products = [alo*blo, alo*bhi, ahi*blo, ahi*bhi]
      lo = nearest(minval(products), -1.0_wp)
      hi = nearest(maxval(products), 1.0_wp)
   end subroutine interval_product

   !* ================================================================================= *!
   !*                                    Finalizer                                      *!
   !* ================================================================================= *!

   !> Release the allocatable components of an SvdW LSF instance
   !>
   !> @param[inout] self LSF instance being finalized
   subroutine finalize_lsf_svdw(self)
      !> LSF instance
      type(moist_cavity_drop_lsf_svdw_type), intent(inout) :: self

      if (allocated(self%radii)) deallocate (self%radii)
      if (allocated(self%active_cand)) deallocate (self%active_cand)
      if (allocated(self%act_atom)) deallocate (self%act_atom)
      if (allocated(self%act_d)) deallocate (self%act_d)
      if (allocated(self%act_radius)) deallocate (self%act_radius)
      if (allocated(self%act_x)) deallocate (self%act_x)
      if (allocated(self%act_inv_x)) deallocate (self%act_inv_x)
      if (allocated(self%act_u)) deallocate (self%act_u)
      if (allocated(self%act_center)) deallocate (self%act_center)
      if (allocated(self%act_bsq)) deallocate (self%act_bsq)
      if (allocated(self%bat_cand)) deallocate (self%bat_cand)
      if (allocated(self%bat_center)) deallocate (self%bat_center)
      if (allocated(self%bat_radius)) deallocate (self%bat_radius)
      if (allocated(self%bat_bsq)) deallocate (self%bat_bsq)
      if (allocated(self%vjp_x)) deallocate (self%vjp_x)
      if (allocated(self%vjp_ix)) deallocate (self%vjp_ix)
      if (allocated(self%vjp_u)) deallocate (self%vjp_u)
      if (allocated(self%vjp_res)) deallocate (self%vjp_res)

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
   end subroutine finalize_lsf_svdw

end module moist_cavity_drop_lsf_svdw
