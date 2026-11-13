!> DROP level-set unit tests
!>
!> Shared checks run one [[moist_cavity_drop_lsf_type]] runner for SvdW, CFC and SSES
!> SvdW sweeps `blend_k` x (`blend_3b`, `blend_2b`); CFC and SSES use compiled defaults
!>
!> Groups follow the registry; a tag in parentheses limits a check to those LSFs
!>
!> Common:
!>   * `sign_convention`                       interior < 0, exterior > 0
!>   * `screening_separation`                  screened vs unscreened pair (SvdW, CFC)
!>
!> Spatial and nuclear derivatives:
!>   * `f1_r_fd`                               grad vs 4-point central FD of f0
!>   * `f2_rr_fd`                              Hessian vs FD of grad
!>   * `f3_rrr_fd`                             third spatial deriv vs FD of Hessian
!>   * `f1_rA_fd`                              nuclear grad vs FD of f0
!>   * `f2_r_rA_fd`                            mixed deriv vs FD of grad
!>   * `f3_rr_rA_fd`                           mixed third deriv vs FD of Hessian
!>   * `f2_rArB`, `f3_r_rArB`                  pure/mixed nuclear seconds
!>   * `f4_rrrr`, `f4_rrr_rA`, `f4_rr_rArB`    fourth-order derivatives
!>
!> Radius derivatives:
!> SvdW rescales one atom's tensors; CFC differentiates a two-center kernel
!>   * `f1_rad_fd`                             radius grad vs FD of f0
!>   * `f2_r_rad_fd`                           mixed spatial-radius deriv vs FD of grad
!>   * `f3_rr_rad_fd`                          mixed third radius deriv vs FD of Hessian
!>   * `radrad_fd`                             `f{2,3,4}_*_radrad` vs a one-radius FD
!>                                             of the `f*_rad` ladder, element by element
!>   * `ra_rad_fd`                             `f{2,3,4}_*_rA_rad` vs a one-radius FD
!>                                             of the `f*_rA` ladder; pins asymmetric
!>                                             atom-slot order
!>   * `hvp_rad_fd`                            radius row of the joint position/radius
!>                                             HVP vs FD of the `f*_rad` ladder along
!>                                             `(v, vr)`
!>   * `hvp_ra_joint_fd`                       nuclear row of the same HVP, `vrad`
!>                                             supplied, vs FD of the `f*_rA` ladder
!>                                             along `(v, vr)`
!>   * `radius_pairwise`                       exact HVP contraction of `f*_radrad` and
!>                                             `f*_rA_rad` at roundoff
!> FD thresholds: [[fd_abs_thr]], [[radius_fd_thr]], [[high_order_thr]]
!>
!> Normalization:
!>   * `normalized_f1_rA`                      normalized LSF nuclear grad
!>
!> Partition:
!>   * `atom_partition`                        SvdW: shares vs `-dS/dR_a`
!>                                             CFC: shares vs sub-structure
!>                                             pseudo-densities
!>                                             SSES: shares vs `(dF/dR_a)/(dq_a/dR_a)`
!>
!> Consistency:
!>   * `f012_hessfree_value_grad`              Hessian-free value+grad vs the full path
!>   * `contracted_vs_pairwise`                every `tangent_*` / `hvp_*`
!>                                             against the uncontracted tensors
!>
!> Exclusion radius (SvdW, CFC):
!>   * `exclusion_radius_never_overclaims`     surface-free certified balls; fade
!>                                             off and on (default, biting band),
!>                                             molecule and lattice, raw `|S|`
!>                                             (CFC `|S|/L`) over-claim as negative
!>                                             control
!>   * `exclusion_radius_gated_on_*`           no certificate outside the bound's
!>                                             assumptions
!>   * `*_certificate_constants`               faded constants vs the derivation
!>   * `*_far_exclusion_radius`                no surface near unreached points
!>   * `svdw_far_exclusion_radius_blend_1b_zero`
!>                                             far certificate off where an
!>                                             unreached point has atoms in reach
!>
!> Unreached prepare (SvdW, CFC):
!>   * `empty_active`                          no atom in reach: `unreached()` from
!>                                             `prepare` / `prepare_subset`, no
!>                                             `error`; SvdW also `Z` underflow with
!>                                             atoms in reach; both a finite value
!>                                             with non-finite higher orders:
!>                                             unreached at `max_deriv = 4`, reached
!>                                             at 0, as the batch
!>   * an accessor after an unreached prepare aborts (`error stop`), not
!>     testable in process
!>
!> Batched evaluation:
!>   * `svdw_batch_vs_accessors`               strip kernel vs per-point accessors,
!>                                             orders 0..4, molecule and lattice,
!>                                             ragged strips, packed slot for slot
!>   * `svdw_batch_candidate_list`             shared list present (sorted, reversed)
!>                                             vs absent, bit for bit
!>   * `svdw_batch_unreached`                  `unreached(column)`, `+Inf` and NaN
!>                                             derivatives on far columns, near
!>                                             columns bit for bit, `blend_1b = 0`,
!>                                             `Z` underflow with screening off
!>   * `svdw_batch_on_nucleus`                 on-nucleus strip vs the per-point path
!>   * `sses_batch_default_packed`             looping default (SSES has no override),
!>                                             orders 3/4 packed, packed Hessian
!>   * `cfc_batch_vs_accessors`                CFC strip kernel vs per-point accessors,
!>                                             orders 0..4, molecule and lattice,
!>                                             packed Hessian as for SvdW
!>   * `cfc_batch_candidate_list`              shared list present (sorted, reversed)
!>                                             vs absent, bit for bit
!>   * `cfc_batch_unreached`                   `unreached(column)`, `+Inf` and NaN
!>                                             derivatives on far columns, near
!>                                             columns at tolerance
!>   * `cfc_batch_scattered_strips`            strips cut into compact sub-strips
!>                                             equal the compact batches
!>   * `cfc_batch_on_nucleus`                  nucleus value finite at order 0,
!>                                             strip-mates bit for bit; from order 1
!>                                             an `error` naming the cusp column
!>   * an `order` the outputs cannot hold aborts (`error stop`), not testable in
!>     process; so does `order >= 2` without `lsf2_rr_packed`
!>   * `svdw_set_full_scan_below`              cell list mode changes, candidates not
!>
!> Batched nuclear vjp (`vjp_f1_rA_batch`, SvdW and CFC; isodensity's identically
!> zero `vjp_f1_rA` has no fixture in this suite):
!>   * `*_vjp_batch_vs_point`                  eight-point tile (plateau, in-band,
!>                                             far) vs `prepare` + `vjp_f1_rA`
!>                                             with the packed seed;
!>                                             exact zeros out of reach, unreached
!>                                             and NaN far out; controls: a packed
!>                                             off-diagonal perturbation and the
!>                                             upper triangle without the factor 2
!>   * `*_vjp_batch_candidate_subset`          strict list vs full list bit for bit,
!>                                             reversed list, sentinel past the
!>                                             tile; control: a dropped reached
!>                                             candidate
!>   * `svdw_vjp_batch_failed_column`          `error` names the column and the call
!>                                             returns
!>   * `svdw_vjp_batch_strips`                 70 points over three strips of a
!>                                             1000-atom lattice, list absent /
!>                                             ascending / reversed, an unreached
!>                                             column mid-strip, the same point in
!>                                             two tiles; control: a dropped
!>                                             reached candidate
!>   * `svdw_vjp_batch_on_nucleus`             two on-nucleus columns equal the
!>                                             per-point path, strip-mates bit for
!>                                             bit those of the tile without them;
!>                                             controls: a point 1e-3 off, a
!>                                             changed seed
!>   * `svdw_vjp_batch_on_nucleus_strips`      on-nucleus columns in two strips, next
!>                                             to an unreached column, an all-nuclei
!>                                             strip, scatter mode, `max_deriv`
!>                                             restored; control: a point 1e-3 off
!>   * `*_vjp_batch_empty_list`                an empty `candidate_indices`: every
!>                                             column unreached, no `error`; control:
!>                                             the full list reaches
!>   * `svdw_vjp_batch_fade_classes`           plateau, in-band and beyond-band atoms
!>                                             in one strip vs per point; control:
!>                                             the band onset moved
!>   * `cfc_vjp_batch_strips`                  70 points over three strips of a
!>                                             1000-atom lattice cut into sub-strips,
!>                                             list absent / ascending / reversed, an
!>                                             unreached column, the same point in two
!>                                             tiles; control: a dropped reached
!>                                             candidate
!>   * `cfc_vjp_batch_on_nucleus_strips`       nuclei (refused by the cusp) in two
!>                                             strips, next to an unreached column, an
!>                                             all-nuclei strip, scatter mode: an
!>                                             `error` naming a nucleus column, the
!>                                             tile without them vs per point;
!>                                             control: the nuclei moved 1e-3
!>   * `cfc_vjp_batch_on_nucleus_mates`        a nucleus in a compact strip is an
!>                                             `error`, the mates alone reached;
!>                                             control: a changed seed
!>   * `cfc_vjp_batch_fade_classes`            plateau, in-band and beyond-band atoms
!>                                             and pairs of every class in one strip
!>                                             vs per point, two bands and the fade
!>                                             off; control: the band onset moved
!>   * `*_vjp_batch_tiny_pd`                   level-set sum in 1e-154..1e-103, lift
!>                                             finite and adjoint not: the column is
!>                                             unreached and NaN, never reached with
!>                                             Inf; controls reached, equal per point
!>   * `bad_cfc_vjp_batch_check_batch`         the CFC override's own refusal names
!>                                             point 1
!>   * `bad_*_vjp_batch_near_nucleus`          a point 1e-100 off a nucleus: lift
!>                                             finite, adjoint not, `S <= 0`; an
!>                                             `error` naming that column, never
!>                                             unreached
!>   * `vjp_batch_shape_messages`              every misshaped argument named; the
!>                                             `error stop` itself is not testable
!>                                             in process
!>
!> Fade band (CFC):
!>   * `cfc_batch_fade_band`                   fade limit bit for bit, biting band
!>                                             moves values, orders 1..4 vs FD
!>   * `cfc_fade_tiling`                       a point's bits independent of its
!>                                             tile, the list honoured, reach = band
!>                                             end, one in-band atom vs `S_4` formula
!>   * `cfc_fade_prepare_vs_batch`             `prepare` + accessors vs the batch,
!>                                             orders 0..4, bands 27/3 and 5/2,
!>                                             unreached column by column
!>   * `cfc_fade_fd_atom_families`             nuclear, pair, radius, radius-pair,
!>                                             nuclear-radius, vjp and joint hvp
!>                                             families vs FD inside one band piece
!>   * `cfc_fade_contracted`                   tangent, hvp, vjp vs full-tensor
!>                                             contractions inside the band
!>
!> Fade band (SvdW):
!>   * `svdw_fade_tiling_invariance`           a point's bits independent of its tile
!>   * `svdw_fade_limit`                       band beyond every atom = fade off, bits
!>   * `svdw_fade_band_reach`                  reach = band end at every threshold
!>   * `svdw_fade_fd`                          orders 1..3 vs FD inside one band piece
!>   * `svdw_fade_fd_atom_families`            nuclear, pair, radius, radius-pair,
!>                                             nuclear-radius, vjp and joint hvp
!>                                             families vs FD inside one band piece
!>
!> Error paths (`should_fail`):
!>   * `bad_svdw_vjp_batch_failed_column`      failed `vjp_f1_rA_batch` column
!>   * `bad_svdw_fade_width_nonpositive`       fade width 0
!>   * `bad_svdw_param_changed_after_update`   fade band changed after `update`
!>   * `bad_svdw_stale_reach`                  reach refreshed under another band
!>   * `bad_cfc_prepare_fade_width_nonpositive` CFC fade width 0, `prepare`
!>   * `bad_cfc_fade_changed_after_update`     CFC fade band changed after `update`
!>   * `bad_cfc_stale_reach`                   CFC reach refreshed under another band
!>   * `bad_cfc_stale_reach_fade_off`          CFC fade off, `c` changed after `update`
!>   * `bad_cfc_fade_width_nonpositive`        CFC fade width 0
!>   * `bad_cfc_batch_exponent_range`          CFC `evaluate_batch` with `-a1 > 708`
!>   * `bad_svdw_blend_k_fade_off`             SvdW `blend_k = 0` with the fade off
!>
!> SvdW-only:
!>   * `point_on_nucleus_is_finite`            finite field and derivatives on a nucleus
!>   * `body_order_scaling`                    1b/2b/3b weight reduction
!>
!> SSES-only:
!>   * `sphere`                                exact single-atom override
!>   * `independent_reference`                 value vs an mpmath KKT solve
!>   * `separable_distant_atom`                a distant atom changes nothing
!>   * `neck_and_invalid_input`                bridged neck, open gap, eps = 0 refused
!>   * `entrant_switch`                        analytic midplane zero vs DROP switching
!>   * `rigid_motion_and_subset`               rigid-motion invariance, unscreened subset
!>   * `field_scale`                           `scale` multiplies every channel
!>   * `drop_surface`                          DROP on the sphere and under a
!>                                             rescaled field
!>   * `envelope_matches_implicit`             envelope atom derivatives vs implicit jets
!>   * `generated_spatial_response`            generated kernels vs seeded jets
!>   * `cache_invalidation`, `thread_clone`    native cache rebuild and pointer repair
!>
!> Structures: MB16-43/Heavy28/Amino20x4/But14diol/UPU23 via `get_test_structures`
!> CPCM radii via `get_test_radii`; samples via `get_test_points`
!>
!> Nuclear FDs use `set_centers` to preserve caches and `max_deriv`
!> Rebuild spatial sort and screening bounds after atom moves
module test_cavity_drop_lsf
   use moist_cavity_drop_lsf_cfc_param, only: moist_cavity_drop_lsf_cfc_param_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use mctc_env_accuracy, only: wp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, &
                                            ieee_quiet_nan, ieee_negative_inf, ieee_positive_inf
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, &
                           center_at_origin, fd4_scalar, fd6_scalar, det_rand
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, &
                                         lsf_vjp_batch_shape_message
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_kernel, only: strip_width_svdw => svdw_strip_width
   use moist_math_packed, only: packed_index, packed_index3, packed_index4, &
                                    unpack_sym_matrix, pack_sym_matrix, pack_sym_part, &
                                    unpack_sym_tensor3, unpack_sym_tensor4
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_cfc_kernel, only: strip_width_cfc => cfc_strip_width
   use moist_cavity_drop_lsf_sses, only: moist_cavity_drop_lsf_sses_type
   use moist_cavity_drop_lsf_sses_algebra, only: sses_jet_type, jet_derivative
   use moist_cavity_drop_lsf_sses_kernel, only: seeded_evaluate => sses_evaluate
   use moist_cavity_drop_lsf_sses_cache, only: sses_cache_type, sses_cache_build, sses_cached_evaluate
   use moist_cavity_drop_lsf_sses_param, only: moist_cavity_drop_lsf_sses_param_type, temperature_length
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid, &
                                    moist_math_smoothstep_polynomial_type, new_smoothstep_polynomial
   use moist_context, only: moist_context_type, new_context
   use moist_radii_custom, only: radius_type_custom, new_custom_radii_atoms
   use moist_model_continuum_gems_utils, only: BuildSuperStructure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_lsf

   integer, parameter :: ndim = 3

   !> Hessian-free consistency tolerance
   !>
   !> SvdW derivative branches use different expression schedules; roundoff agreement
   real(wp), parameter :: HESSFREE_ABS = 1.0e-14_wp
   real(wp), parameter :: HESSFREE_REL = 1.0e-12_wp

   !> Batch vs per-point tolerance, normwise per tensor
   !>
   !> Same lift statements at different `np`, power sums in candidate vs user order
   !>
   !> - orders 0..3: 1e-11; measured 2026-11-03 at most 2.3e-12 (third order)
   !> - fourth order: 1e-10; measured up to 3.6e-11 at a far point, where a
   !>   1-ulp perturbation of the power sums alone moves the per-point lift by
   !>   up to 2.3e-10 normwise (conditioning probe on the strip kernel)
   real(wp), parameter :: BATCH_REL = 1.0e-11_wp
   real(wp), parameter :: BATCH_REL4 = 1.0e-10_wp
   real(wp), parameter :: BATCH_ABS = 1.0e-13_wp

   !> Batched vs per-point nuclear vjp, normwise relative per column
   !>
   !> The looping default runs the per-point `vjp_f1_rA` on the same prepared point,
   !> so the two differ only by the candidate order of the power sums; the SvdW
   !> strip override (adjoint lift plus nuclear stage) measured at most 4.7e-14
   !> over the 70-point, three-strip tile
   real(wp), parameter :: VJP_BATCH_REL = 1.0e-12_wp
   !> Smallest deviation a deliberately wrong input must produce in the vjp controls
   real(wp), parameter :: VJP_BATCH_MISS = 1.0e-9_wp

   real(wp), parameter :: STEP_SIZE = 1.0e-3_wp

   !> FD thresholds
   real(wp), parameter :: ABS_THR = 2.0e-10_wp
   real(wp), parameter :: REL_THR = 1.0e-9_wp
   real(wp), parameter :: CFC_ABS_THR = 5.0e-9_wp
   real(wp), parameter :: CFC_REL_THR = 5.0e-8_wp
   !> SSES FD thresholds
   !>
   !> The moment solver stops at a residual of `1e-11 (1 + |m|)`; derivatives
   !> linear in the moment carry that error, which the stencil divides by
   !> `STEP_SIZE`. A tighter solver passes the SvdW thresholds
   real(wp), parameter :: SSES_ABS_THR = 5.0e-8_wp
   real(wp), parameter :: SSES_REL_THR = 5.0e-7_wp

   !> SvdW fade in the user-space nuclear and radius FD sweeps: off
   !>
   !> - those sweeps index nuclear outputs by user id and need every atom
   !>   active; with the fade on, atoms beyond the band end have no slot
   !> - the faded per-atom families are FD-checked inside a band by
   !>   `svdw_fade_fd_atom_families`; spatial sweeps keep the fade on
   logical, parameter :: fd_user_space_fade = .false.

   !> Tolerance of the atomic shares, which are O(1)
   real(wp), parameter :: PARTITION_THR = 1.0e-12_wp

   !> FD sampling points per structure
   integer, parameter :: n_points_svdw = 4
   integer, parameter :: n_points_cfc = 1
   integer, parameter :: n_points_sses = 1
   !> SvdW high-order sampling points, the [[get_test_points]] default
   integer, parameter :: n_points_svdw_high_order = 7

   integer, parameter :: n_svdw_blends = 3
   integer, parameter :: n_svdw_gammas = 2
   real(wp), parameter :: svdw_blend_k_values(n_svdw_blends) = [0.1_wp, 1.0_wp, 4.0_wp]
   real(wp), parameter :: svdw_gamma_values(n_svdw_gammas) = [0.0_wp, 0.7_wp]
   real(wp), parameter :: svdw_2b_values(n_svdw_gammas) = [0.5_wp, 0.0_wp]

   ! Legacy SvdW smoothing
   real(wp), parameter :: svdw_legacy_blend_k = 3.0_wp
   real(wp), parameter :: svdw_legacy_blend_2b = 1.0_wp
   real(wp), parameter :: svdw_legacy_blend_3b = 1.0_wp

contains

   subroutine collect_cavity_drop_lsf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  !> Common
                  new_unittest("svdw_sign_convention", test_svdw_sign_convention), &
                  new_unittest("cfc_sign_convention", test_cfc_sign_convention), &
                  new_unittest("sses_sign_convention", test_sses_sign_convention), &
                  new_unittest("svdw_screening_separation", test_svdw_screening_separation), &
                  new_unittest("cfc_screening_separation", test_cfc_screening_separation), &
                  !> Spatial and nuclear derivatives
                  new_unittest("svdw_f1_r_fd", test_svdw_f1_r_fd), &
                  new_unittest("cfc_f1_r_fd", test_cfc_f1_r_fd), &
                  new_unittest("sses_f1_r_fd", test_sses_f1_r_fd), &
                  new_unittest("svdw_f2_rr_fd", test_svdw_f2_rr_fd), &
                  new_unittest("cfc_f2_rr_fd", test_cfc_f2_rr_fd), &
                  new_unittest("sses_f2_rr_fd", test_sses_f2_rr_fd), &
                  new_unittest("svdw_f3_rrr_fd", test_svdw_f3_rrr_fd), &
                  new_unittest("cfc_f3_rrr_fd", test_cfc_f3_rrr_fd), &
                  new_unittest("sses_f3_rrr_fd", test_sses_f3_rrr_fd), &
                  new_unittest("svdw_f1_ra_fd", test_svdw_f1_rA_fd), &
                  new_unittest("cfc_f1_ra_fd", test_cfc_f1_rA_fd), &
                  new_unittest("sses_f1_ra_fd", test_sses_f1_rA_fd), &
                  new_unittest("svdw_f2_r_ra_fd", test_svdw_f2_r_rA_fd), &
                  new_unittest("cfc_f2_r_ra_fd", test_cfc_f2_r_rA_fd), &
                  new_unittest("sses_f2_r_ra_fd", test_sses_f2_r_rA_fd), &
                  new_unittest("svdw_f3_rr_ra_fd", test_svdw_f3_rr_rA_fd), &
                  new_unittest("cfc_f3_rr_ra_fd", test_cfc_f3_rr_rA_fd), &
                  new_unittest("sses_f3_rr_ra_fd", test_sses_f3_rr_rA_fd), &
                  new_unittest("svdw_f2_rarb", test_svdw_f2_rArB), &
                  new_unittest("cfc_f2_rarb", test_cfc_f2_rArB), &
                  new_unittest("sses_f2_rarb", test_sses_f2_rArB), &
                  new_unittest("svdw_f3_r_rarb", test_svdw_f3_r_rArB), &
                  new_unittest("cfc_f3_r_rarb", test_cfc_f3_r_rArB), &
                  new_unittest("sses_f3_r_rarb", test_sses_f3_r_rArB), &
                  new_unittest("svdw_f4_rrrr", test_svdw_f4_rrrr), &
                  new_unittest("cfc_f4_rrrr", test_cfc_f4_rrrr), &
                  new_unittest("sses_f4_rrrr", test_sses_f4_rrrr), &
                  new_unittest("svdw_f4_rrr_ra", test_svdw_f4_rrr_rA), &
                  new_unittest("cfc_f4_rrr_ra", test_cfc_f4_rrr_rA), &
                  new_unittest("sses_f4_rrr_ra", test_sses_f4_rrr_rA), &
                  new_unittest("svdw_f4_rr_rarb", test_svdw_f4_rr_rArB), &
                  new_unittest("cfc_f4_rr_rarb", test_cfc_f4_rr_rArB), &
                  new_unittest("sses_f4_rr_rarb", test_sses_f4_rr_rArB), &
                  !> Radius derivatives
                  new_unittest("svdw_f1_rad_fd", test_svdw_f1_rad_fd), &
                  new_unittest("cfc_f1_rad_fd", test_cfc_f1_rad_fd), &
                  new_unittest("sses_f1_rad_fd", test_sses_f1_rad_fd), &
                  new_unittest("svdw_f2_r_rad_fd", test_svdw_f2_r_rad_fd), &
                  new_unittest("cfc_f2_r_rad_fd", test_cfc_f2_r_rad_fd), &
                  new_unittest("sses_f2_r_rad_fd", test_sses_f2_r_rad_fd), &
                  new_unittest("svdw_f3_rr_rad_fd", test_svdw_f3_rr_rad_fd), &
                  new_unittest("cfc_f3_rr_rad_fd", test_cfc_f3_rr_rad_fd), &
                  new_unittest("sses_f3_rr_rad_fd", test_sses_f3_rr_rad_fd), &
                  new_unittest("svdw_radrad_fd", test_svdw_radrad_fd), &
                  new_unittest("cfc_radrad_fd", test_cfc_radrad_fd), &
                  new_unittest("sses_radrad_fd", test_sses_radrad_fd), &
                  new_unittest("svdw_ra_rad_fd", test_svdw_rA_rad_fd), &
                  new_unittest("cfc_ra_rad_fd", test_cfc_rA_rad_fd), &
                  new_unittest("sses_ra_rad_fd", test_sses_rA_rad_fd), &
                  new_unittest("svdw_hvp_rad_fd", test_svdw_hvp_rad_fd), &
                  new_unittest("cfc_hvp_rad_fd", test_cfc_hvp_rad_fd), &
                  new_unittest("sses_hvp_rad_fd", test_sses_hvp_rad_fd), &
                  new_unittest("svdw_hvp_ra_joint_fd", test_svdw_hvp_rA_joint_fd), &
                  new_unittest("cfc_hvp_ra_joint_fd", test_cfc_hvp_rA_joint_fd), &
                  new_unittest("sses_hvp_ra_joint_fd", test_sses_hvp_rA_joint_fd), &
                  new_unittest("svdw_radius_pairwise", test_svdw_radius_pairwise), &
                  new_unittest("cfc_radius_pairwise", test_cfc_radius_pairwise), &
                  new_unittest("sses_radius_pairwise", test_sses_radius_pairwise), &
                  !> Normalization
                  new_unittest("svdw_normalized_f1_ra", test_svdw_normalized_f1_rA), &
                  new_unittest("cfc_normalized_f1_ra", test_cfc_normalized_f1_rA), &
                  new_unittest("sses_normalized_f1_ra", test_sses_normalized_f1_rA), &
                  !> Partition
                  new_unittest("svdw_atom_partition", test_svdw_atom_partition), &
                  new_unittest("cfc_atom_partition", test_cfc_atom_partition), &
                  new_unittest("sses_atom_partition", test_sses_atom_partition), &
                  !> Consistency
                  new_unittest("svdw_f012_hessfree_value_grad", test_svdw_f012_hessfree), &
                  new_unittest("cfc_f012_hessfree_value_grad", test_cfc_f012_hessfree), &
                  new_unittest("sses_f012_hessfree_value_grad", test_sses_f012_hessfree), &
                  new_unittest("svdw_contracted_vs_pairwise", test_svdw_contracted), &
                  new_unittest("cfc_contracted_vs_pairwise", test_cfc_contracted), &
                  new_unittest("sses_contracted_vs_pairwise", test_sses_contracted), &
                  !> Exclusion radius
                  new_unittest("svdw_exclusion_radius_never_overclaims", &
                               test_svdw_exclusion_radius), &
                  new_unittest("svdw_exclusion_radius_gated_on_blends", &
                               test_svdw_exclusion_radius_gate), &
                  new_unittest("svdw_certificate_constants", test_svdw_certificate_constants), &
                  new_unittest("svdw_far_exclusion_radius", test_svdw_far_exclusion_radius), &
                  new_unittest("svdw_far_exclusion_radius_blend_1b_zero", &
                               test_svdw_far_blend_1b_zero), &
                  new_unittest("cfc_exclusion_radius_never_overclaims", &
                               test_cfc_exclusion_radius), &
                  new_unittest("cfc_exclusion_radius_gated_on_exponents", &
                               test_cfc_exclusion_radius_gate), &
                  new_unittest("cfc_certificate_constants", test_cfc_certificate_constants), &
                  new_unittest("cfc_far_exclusion_radius", test_cfc_far_exclusion_radius), &
                  !> Empty active lists
                  new_unittest("svdw_empty_active", test_svdw_empty_active), &
                  new_unittest("cfc_empty_active", test_cfc_empty_active), &
                  !> Candidate cell list
                  new_unittest("svdw_box_candidates_match_full_scan", &
                               test_svdw_box_candidates), &
                  !> Batched evaluation
                  new_unittest("svdw_batch_vs_accessors", test_svdw_batch_vs_accessors), &
                  new_unittest("svdw_batch_candidate_list", test_svdw_batch_candidate_list), &
                  new_unittest("svdw_batch_unreached", test_svdw_batch_unreached), &
                  new_unittest("svdw_batch_on_nucleus", test_svdw_batch_on_nucleus), &
                  new_unittest("svdw_fade_tiling_invariance", test_svdw_fade_tiling), &
                  new_unittest("svdw_fade_limit", test_svdw_fade_limit), &
                  new_unittest("svdw_fade_band_reach", test_svdw_fade_reach), &
                  new_unittest("bad_svdw_fade_width_nonpositive", test_bad_fade_width, &
                               should_fail=.true.), &
                  new_unittest("bad_svdw_param_changed_after_update", test_bad_param_changed, &
                               should_fail=.true.), &
                  new_unittest("bad_svdw_stale_reach", test_bad_svdw_stale_reach, &
                               should_fail=.true.), &
                  new_unittest("svdw_fade_fd", test_svdw_fade_fd), &
                  new_unittest("svdw_fade_fd_atom_families", test_svdw_fade_fd_atom_families), &
                  new_unittest("svdw_set_full_scan_below", test_svdw_set_full_scan_below), &
                  !> Batched nuclear vjp
                  new_unittest("svdw_vjp_batch_vs_point", test_svdw_vjp_batch_vs_point), &
                  new_unittest("cfc_vjp_batch_vs_point", test_cfc_vjp_batch_vs_point), &
                  new_unittest("svdw_vjp_batch_candidate_subset", &
                               test_svdw_vjp_batch_candidate_subset), &
                  new_unittest("cfc_vjp_batch_candidate_subset", &
                               test_cfc_vjp_batch_candidate_subset), &
                  new_unittest("svdw_vjp_batch_failed_column", &
                               test_svdw_vjp_batch_failed_column), &
                  new_unittest("svdw_vjp_batch_strips", test_svdw_vjp_batch_strips), &
                  new_unittest("svdw_vjp_batch_on_nucleus", &
                               test_svdw_vjp_batch_on_nucleus), &
                  new_unittest("svdw_vjp_batch_on_nucleus_strips", &
                               test_svdw_vjp_batch_on_nucleus_strips), &
                  new_unittest("svdw_vjp_batch_empty_list", &
                               test_svdw_vjp_batch_empty_list), &
                  new_unittest("cfc_vjp_batch_empty_list", &
                               test_cfc_vjp_batch_empty_list), &
                  new_unittest("svdw_vjp_batch_fade_classes", &
                               test_svdw_vjp_batch_fade_classes), &
                  new_unittest("cfc_vjp_batch_strips", test_cfc_vjp_batch_strips), &
                  new_unittest("cfc_vjp_batch_on_nucleus_strips", &
                               test_cfc_vjp_batch_on_nucleus_strips), &
                  new_unittest("cfc_vjp_batch_on_nucleus_mates", &
                               test_cfc_vjp_batch_on_nucleus_mates), &
                  new_unittest("cfc_vjp_batch_fade_classes", &
                               test_cfc_vjp_batch_fade_classes), &
                  new_unittest("svdw_vjp_batch_tiny_pd", test_svdw_vjp_batch_tiny_pd), &
                  new_unittest("cfc_vjp_batch_tiny_pd", test_cfc_vjp_batch_tiny_pd), &
                  new_unittest("vjp_batch_shape_messages", &
                               test_vjp_batch_shape_messages), &
                  new_unittest("bad_svdw_vjp_batch_failed_column", &
                               test_bad_svdw_vjp_batch_failed_column, &
                               should_fail=.true.), &
                  new_unittest("bad_cfc_vjp_batch_check_batch", &
                               test_bad_cfc_vjp_batch_check_batch, should_fail=.true.), &
                  new_unittest("bad_svdw_vjp_batch_near_nucleus", &
                               test_bad_svdw_vjp_batch_near_nucleus, should_fail=.true.), &
                  new_unittest("bad_cfc_vjp_batch_near_nucleus", &
                               test_bad_cfc_vjp_batch_near_nucleus, should_fail=.true.), &
                  new_unittest("sses_batch_default_packed", test_sses_batch_default_packed), &
                  new_unittest("cfc_batch_vs_accessors", test_cfc_batch_vs_accessors), &
                  new_unittest("cfc_batch_candidate_list", test_cfc_batch_candidate_list), &
                  new_unittest("cfc_batch_unreached", test_cfc_batch_unreached), &
                  new_unittest("cfc_batch_scattered_strips", test_cfc_batch_scattered_strips), &
                  new_unittest("cfc_batch_on_nucleus", test_cfc_batch_on_nucleus), &
                  new_unittest("cfc_batch_fade_band", test_cfc_batch_fade_band), &
                  new_unittest("cfc_fade_tiling", test_cfc_fade_tiling), &
                  new_unittest("cfc_fade_prepare_vs_batch", test_cfc_fade_prepare_vs_batch), &
                  new_unittest("cfc_fade_fd_atom_families", test_cfc_fade_fd_atom_families), &
                  new_unittest("cfc_fade_contracted", test_cfc_fade_contracted), &
                  new_unittest("bad_cfc_prepare_fade_width_nonpositive", &
                               test_bad_cfc_prepare_fade_width, should_fail=.true.), &
                  new_unittest("bad_cfc_fade_changed_after_update", test_bad_cfc_fade_changed, &
                               should_fail=.true.), &
                  new_unittest("bad_cfc_stale_reach", test_bad_cfc_stale_reach, &
                               should_fail=.true.), &
                  new_unittest("bad_cfc_stale_reach_fade_off", test_bad_cfc_stale_reach_off, &
                               should_fail=.true.), &
                  new_unittest("bad_cfc_fade_width_nonpositive", test_bad_cfc_fade_width, &
                               should_fail=.true.), &
                  new_unittest("bad_cfc_batch_exponent_range", test_bad_cfc_exponent_range, &
                               should_fail=.true.), &
                  new_unittest("bad_svdw_blend_k_fade_off", test_bad_svdw_blend_k_fade_off, &
                               should_fail=.true.), &
                  !> SvdW-only
                  new_unittest("svdw_point_on_nucleus_is_finite", test_svdw_on_nucleus), &
                  new_unittest("svdw_body_order_scaling", test_svdw_body_order_scaling), &
                  !> SSES-only
                  new_unittest("sses_sphere", test_sses_sphere), &
                  new_unittest("sses_independent_reference", test_sses_reference), &
                  new_unittest("sses_separable_distant_atom", test_sses_separable), &
                  new_unittest("sses_neck_and_invalid_input", test_sses_geometry), &
                  new_unittest("sses_entrant_switch", test_sses_entrant_switch), &
                  new_unittest("sses_rigid_motion_and_subset", test_sses_invariance), &
                  new_unittest("sses_field_scale", test_sses_field_scale), &
                  new_unittest("sses_drop_surface", test_sses_drop), &
                  new_unittest("sses_envelope_matches_implicit", test_sses_envelope), &
                  new_unittest("sses_generated_spatial_response", test_sses_generated_spatial), &
                  new_unittest("sses_cache_invalidation", test_sses_cache_invalidation), &
                  new_unittest("sses_thread_clone", test_sses_thread_clone) &
                  ]
   end subroutine collect_cavity_drop_lsf

   !* ================================================================================= *!
   !*                              Local helpers                                        *!
   !* ================================================================================= *!

   !> Create and bind a polymorphic LSF
   !>
   !> SvdW-only blend overrides; CFC uses fixed shape parameters, the fade switch
   !> and band of SvdW and CFC from the fade optionals
   !> Set screening threshold before `new`/`update` propagates it to SSD
   !> Default threshold 0 disables screening
   !>
   !> @param[out] lsf                 polymorphic LSF allocatable
   !> @param[in]  mol                 molecular structure to bind
   !> @param[in]  radii               per-atom radii (size mol%nat)
   !> @param[in]  max_deriv           highest spatial derivative order to enable
   !> @param[in]  kind                "svdw", "cfc" or "sses"
   !> @param[in]  blend_k             optional svdw blend sharpness override
   !> @param[in]  blend_2b            optional svdw 2-body weight override
   !> @param[in]  blend_3b            optional svdw 3-body weight override
   !> @param[in]  screening_threshold optional screening cutoff (default 0)
   subroutine init_lsf(lsf, mol, radii, max_deriv, kind, blend_k, blend_2b, &
                       blend_3b, screening_threshold, fade_enabled, fade_start, fade_width)
      class(moist_cavity_drop_lsf_type), allocatable, intent(out) :: lsf
      !> Bound molecular structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii (size mol%nat)
      real(wp), intent(in) :: radii(:)
      !> Maximum spatial derivative order
      integer, intent(in) :: max_deriv
      !> Concrete kind selector: "svdw", "cfc" or "sses"
      character(len=*), intent(in) :: kind
      !> Optional SvdW blend_k override
      real(wp), intent(in), optional :: blend_k
      !> Optional SvdW blend_2b override
      real(wp), intent(in), optional :: blend_2b
      !> Optional SvdW blend_3b override
      real(wp), intent(in), optional :: blend_3b
      !> Optional screening threshold (default 0; SvdW with the fade on ignores it)
      real(wp), intent(in), optional :: screening_threshold
      !> Optional fade switch (default: the parameter default, on for SvdW and CFC)
      logical, intent(in), optional :: fade_enabled
      !> Optional SvdW / CFC fade onset in e-folds
      real(wp), intent(in), optional :: fade_start
      !> Optional SvdW / CFC fade width in e-folds
      real(wp), intent(in), optional :: fade_width

      type(moist_cavity_drop_lsf_svdw_type) :: tmp_svdw
      type(moist_cavity_drop_lsf_cfc_type)  :: tmp_cfc
      type(moist_cavity_drop_lsf_sses_type) :: tmp_sses
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc_param
      real(wp) :: thr

      thr = 0.0_wp
      if (present(screening_threshold)) thr = screening_threshold

      select case (kind)
      case ("svdw")
         tmp_svdw%screening_threshold = thr
         call param%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b, &
                        fade_enabled=fade_enabled, fade_start=fade_start, fade_width=fade_width)
         call tmp_svdw%new(param)
         call tmp_svdw%update(mol, radii)
         call tmp_svdw%set_max_deriv(max_deriv)
         allocate (lsf, source=tmp_svdw)
      case ("cfc")
         tmp_cfc%screening_threshold = thr
         call cfc_param%new(fade_enabled=fade_enabled, fade_start=fade_start, &
                            fade_width=fade_width)
         call tmp_cfc%new(cfc_param)
         call tmp_cfc%update(mol, radii)
         call tmp_cfc%set_max_deriv(max_deriv)
         allocate (lsf, source=tmp_cfc)
      case ("sses")
         tmp_sses%screening_threshold = thr
         call tmp_sses%new()
         call tmp_sses%update(mol, radii)
         call tmp_sses%set_max_deriv(max_deriv)
         allocate (lsf, source=tmp_sses)
      case default
         error stop "init_lsf: unknown kind '"//kind//"'"
      end select
   end subroutine init_lsf

   !> Refresh atom positions without resetting LSF caches or `max_deriv`
   !>
   !> Base API `set_centers` needs no concrete dispatch
   !> `radii` retained for call-site clarity
   !>
   !> @param[inout] lsf      polymorphic LSF (must be allocated)
   !> @param[in]    centers  perturbed positions (3, mol%nat)
   !> @param[in]    radii    per-atom radii (size mol%nat), unchanged
   subroutine refresh_ssd(lsf, centers, radii)
      !> Polymorphic LSF whose geometry is refreshed in place
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Perturbed atom positions (3, mol%nat)
      real(wp), intent(in) :: centers(:, :)
      !> Per-atom radii (size mol%nat), unchanged by this call
      real(wp), intent(in) :: radii(:)

      if (size(radii) /= size(centers, 2)) then
         error stop "refresh_ssd: radii/centers size mismatch"
      end if
      call lsf%set_centers(centers)
   end subroutine refresh_ssd

   !> Unpack packed symmetric matrices of trailing indices to full storage
   !>
   !> - the trailing indices of `packed` and `full` are flattened to `m`
   !>
   !> @param[in]  packed  Packed matrices (6, m)
   !> @param[out] full    Full matrices (3, 3, m)
   !> @param[in]  m       Number of flattened trailing entries
   pure subroutine unpack_pair_batch(packed, full, m)
      !> Packed matrices
      real(wp), intent(in) :: packed(6, m)
      !> Full matrices
      real(wp), intent(out) :: full(ndim, ndim, m)
      !> Flattened trailing entries
      integer, intent(in) :: m

      integer :: i

      do i = 1, m
         call unpack_sym_matrix(packed(:, i), full(:, :, i), ndim)
      end do
   end subroutine unpack_pair_batch

   !> Unpack packed symmetric rank-3 tensors of trailing indices to full storage
   !>
   !> - the trailing indices of `packed` and `full` are flattened to `m`
   !>
   !> @param[in]  packed  Packed tensors (10, m)
   !> @param[out] full    Full tensors (3, 3, 3, m)
   !> @param[in]  m       Number of flattened trailing entries
   pure subroutine unpack_triple_batch(packed, full, m)
      !> Packed tensors
      real(wp), intent(in) :: packed(10, m)
      !> Full tensors
      real(wp), intent(out) :: full(ndim, ndim, ndim, m)
      !> Flattened trailing entries
      integer, intent(in) :: m

      integer :: i

      do i = 1, m
         call unpack_sym_tensor3(packed(:, i), full(:, :, :, i), ndim)
      end do
   end subroutine unpack_triple_batch

   !> Prepare a displaced LSF or fail the test
   !>
   !> Base preparation checked beside active-list identity
   !> Displaced preparation failures would reuse prior values in FD quotients
   !> `.true.` signals immediate caller return
   !>
   !> @param[inout] lsf   Polymorphic LSF, prepared in place
   !> @param[in]    point Sampling point
   !> @param[out]   error Test error, allocated iff `prepare` failed
   logical function prepare_failed(lsf, point, error)
      !> Polymorphic LSF to prepare at `point`
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Sampling point
      real(wp), intent(in) :: point(:)
      !> Test error, allocated iff `prepare` failed
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: lsf_err

      call lsf%prepare(point, lsf_err)
      prepare_failed = allocated(lsf_err)
      if (prepare_failed) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
      end if
   end function prepare_failed

   !* ================================================================================= *!
   !*                              Sign + cutoff checks                                 *!
   !* ================================================================================= *!

   !> SvdW dispatch for the sign-convention check
   subroutine test_svdw_sign_convention(error)
      type(error_type), allocatable, intent(out) :: error
      call run_sign_convention(error, "svdw")
   end subroutine test_svdw_sign_convention

   !> CFC dispatch for the sign-convention check
   subroutine test_cfc_sign_convention(error)
      type(error_type), allocatable, intent(out) :: error
      call run_sign_convention(error, "cfc")
   end subroutine test_cfc_sign_convention

   !> SSES dispatch for the sign-convention check
   subroutine test_sses_sign_convention(error)
      type(error_type), allocatable, intent(out) :: error
      call run_sign_convention(error, "sses")
   end subroutine test_sses_sign_convention

   !> Check negative interior and positive exterior LSF on LiH
   subroutine run_sign_convention(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)
      real(wp) :: inside_pt(ndim), outside_pt(ndim)
      real(wp) :: val_inside, val_outside
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 0, kind)

      ! Atom 1 center, inside the cavity
      inside_pt = mol%xyz(:, 1) + [0.05_wp, 0.05_wp, 0.0_wp]
      call lsf%prepare(inside_pt, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      call lsf%f0(val_inside)
      if (val_inside >= 0.0_wp) then
         call test_failed(error, "Expected lsf0 < 0 deep inside the cavity")
         return
      end if

      ! Beyond molecular bounds along +x; SvdW inside its fade plateau (onset
      ! 11.45 Bohr off the surfaces at k = 5.5), CFC inside its reach (4.33 R_A
      ! from the nucleus at the default band), as an unreached point has no value
      ! on the per-point path
      outside_pt = [maxval(mol%xyz(1, :)) + merge(20.0_wp, 8.0_wp, kind == "sses"), &
                    0.0_wp, 0.0_wp]
      call lsf%prepare(outside_pt, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed outside: "//lsf_err%message)
         return
      end if
      call check(error, .not. lsf%unreached(), more="the outside point must be reached")
      if (allocated(error)) return
      call lsf%f0(val_outside)
      if (val_outside <= 0.0_wp) then
         call test_failed(error, "Expected lsf0 > 0 well outside the cavity")
         return
      end if
   end subroutine run_sign_convention

   !> SvdW dispatch for the molecule-pair separation check
   subroutine test_svdw_screening_separation(error)
      type(error_type), allocatable, intent(out) :: error
      call run_screening_separation(error, "svdw")
   end subroutine test_svdw_screening_separation

   !> CFC dispatch for the molecule-pair separation check
   subroutine test_cfc_screening_separation(error)
      type(error_type), allocatable, intent(out) :: error
      call run_screening_separation(error, "cfc")
   end subroutine test_cfc_screening_separation

   !> Check screening as a copied molecule separates from fixed centers
   !>
   !> Sample near fixed centers at `val ~ 0`; vary copy separation
   !> Adjacent copy contributes; copy beyond ~25-30 bohr gets screened
   !> Screened LSF tracks unscreened reference within ~`ntot**2 * thr`
   !>
   !> [[BuildSuperStructure]] appends moving copy after fixed centers
   !> Moving atom IDs > `nc` expose contribute-to-screen transition
   !>
   !> Legacy SvdW weights keep 2b/3b at unity for `ntot**2 * thr` bound
   !>
   !> SvdW: the threshold sweep runs with the fade off (the only mode that reads
   !> the threshold); with the fade on, screened and unscreened agree bit for bit
   subroutine run_screening_separation(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: center_mol, moving_mol, super_mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_ref, lsf_scr
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: thr, gap, shift, val_ref, val_scr
      integer :: icase, ithr, istep, ip, nc, ntot

      real(wp), parameter :: separation_thresholds(6) = [ &
                             1.0e-8_wp, 1.0E-9_wp, 1.0e-10_wp, 1.0e-11_wp, 1.0e-12_wp, 1.0e-13_wp]

      integer, parameter :: separation_n_steps = 100
      real(wp), parameter :: separation_gap_step = 0.25_wp
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols, 10)

      do icase = 1, size(mols)
         ! Fixed centers and translated copy
         center_mol = mols(icase)
         call center_at_origin(center_mol)
         nc = center_mol%nat
         call get_test_points(center_mol, points, 12)

         do ithr = 1, size(separation_thresholds)
            thr = separation_thresholds(ithr)

            do istep = 0, separation_n_steps
               gap = real(istep, wp)*separation_gap_step

               ! Join fixed and moving molecules
               moving_mol = center_mol
               shift = maxval(center_mol%xyz(1, :)) - minval(moving_mol%xyz(1, :)) + gap
               moving_mol%xyz(1, :) = moving_mol%xyz(1, :) + shift
               call BuildSuperStructure(center_mol, moving_mol, super_mol, &
                                        no_displacement=.true.)
               ntot = super_mol%nat

               call get_test_radii(super_mol, radii)
               call init_lsf(lsf_ref, super_mol, radii, 0, kind, &
                             blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                             blend_3b=svdw_legacy_blend_3b, screening_threshold=0.0_wp, &
                             fade_enabled=.false.)
               call init_lsf(lsf_scr, super_mol, radii, 0, kind, &
                             blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                             blend_3b=svdw_legacy_blend_3b, screening_threshold=thr, &
                             fade_enabled=.false.)

               do ip = 1, size(points, 2)
                  call lsf_ref%prepare(points(:, ip), lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call lsf_ref%f0(val_ref)
                  call lsf_scr%prepare(points(:, ip), lsf_err)
                  call lsf_scr%f0(val_scr)

                  ! Three-body error scales as `thr * ntot**2`
                  call check(error, val_scr, val_ref, thr=real(ntot, wp)**2*thr, &
                             more="screened lsf diverged from unscreened during separation (lower than natoms^2*thr)")
                  if (allocated(error)) return

               end do
               deallocate (lsf_ref, lsf_scr)

               ! SvdW fade on: the threshold is not read, screened = unscreened bit for bit
               if (kind /= "svdw" .or. ithr > 1 .or. mod(istep, 10) /= 0) cycle
               call init_lsf(lsf_ref, super_mol, radii, 0, kind, &
                             blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                             blend_3b=svdw_legacy_blend_3b, screening_threshold=0.0_wp)
               call init_lsf(lsf_scr, super_mol, radii, 0, kind, &
                             blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                             blend_3b=svdw_legacy_blend_3b, screening_threshold=thr)
               do ip = 1, size(points, 2)
                  call lsf_ref%prepare(points(:, ip), lsf_err)
                  if (.not. allocated(lsf_err)) call lsf_scr%prepare(points(:, ip), lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call lsf_ref%f0(val_ref)
                  call lsf_scr%f0(val_scr)
                  if (val_scr /= val_ref) then
                     call test_failed(error, "fade on: the screening threshold changes S")
                     return
                  end if
               end do
               deallocate (lsf_ref, lsf_scr)
            end do
         end do
      end do

   end subroutine run_screening_separation

   !> Active moving atom after latest `prepare` (user ID > `nc`)
   logical function any_moving_active(lsf, nc) result(active)
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      integer, intent(in) :: nc
      integer :: i

      active = .false.
      do i = 1, lsf%active_count()
         if (lsf%active_atom(i) > nc) then
            active = .true.
            return
         end if
      end do
   end function any_moving_active

   !> Require unscreened active-list identity
   !>
   !> Accessor slot `i` belongs to `active_atom(i)`; higher slots unwritten
   !> FD buffers use `nat` and user IDs; valid only with screening off
   !> Guard against uninitialised comparisons from screened fixtures
   !>
   !> Check after first `prepare`; FD displacements preserve active list
   !>
   !> @param[out] error Set if the active list is screened or reordered
   !> @param[in]  lsf   LSF instance, after `prepare`
   !> @param[in]  nat   User-space atom count the caller indexes with
   subroutine check_identity_active(error, lsf, nat)
      type(error_type), allocatable, intent(out) :: error
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      integer, intent(in) :: nat

      integer :: i

      call check(error, lsf%active_count(), nat, &
                 more="active list must be unscreened for user-space indexing")
      if (allocated(error)) return
      do i = 1, nat
         call check(error, lsf%active_atom(i), i, &
                    more="unscreened active list must be the identity")
         if (allocated(error)) return
      end do
   end subroutine check_identity_active

   !* ================================================================================= *!
   !*                           Spatial-derivative FD tests                             *!
   !* ================================================================================= *!

   !> SvdW dispatch for the spatial gradient FD check
   subroutine test_svdw_f1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_r_fd(error, "svdw")
   end subroutine test_svdw_f1_r_fd

   !> CFC dispatch for the spatial gradient FD check
   subroutine test_cfc_f1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_r_fd(error, "cfc")
   end subroutine test_cfc_f1_r_fd

   !> SSES dispatch for the spatial gradient FD check
   subroutine test_sses_f1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_r_fd(error, "sses")
   end subroutine test_sses_f1_r_fd

   !> grad vs 4-point central FD of f0
   subroutine run_f1_r_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, axis, i, iblend, igamma, nblend, ngamma
      real(wp) :: analytic(ndim), numeric(ndim), point(ndim), shifted(ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      type(mctc_error), allocatable :: lsf_err

      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 1, kind, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call lsf%f012_r(lsf1_r=analytic)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f0(f_pp)
                     shifted = point; shifted(axis) = point(axis) + STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f0(f_p)
                     shifted = point; shifted(axis) = point(axis) - STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f0(f_m)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f0(f_mm)
                     call fd4_scalar(f_pp, f_p, f_m, f_mm, STEP_SIZE, numeric(axis), error)
                     if (allocated(error)) return
                  end do
                  do i = 1, ndim
                     call check(error, analytic(i), numeric(i), &
                                thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                     if (allocated(error)) return
                  end do
               end do
               deallocate (lsf)
            end do
         end do
      end do
   end subroutine run_f1_r_fd

   !> SvdW dispatch for the spatial Hessian FD check
   subroutine test_svdw_f2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rr_fd(error, "svdw")
   end subroutine test_svdw_f2_rr_fd

   !> CFC dispatch for the spatial Hessian FD check
   subroutine test_cfc_f2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rr_fd(error, "cfc")
   end subroutine test_cfc_f2_rr_fd

   !> SSES dispatch for the spatial Hessian FD check
   subroutine test_sses_f2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rr_fd(error, "sses")
   end subroutine test_sses_f2_rr_fd

   !> Hessian vs 4-point FD of grad
   subroutine run_f2_rr_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, axis, i, j, iblend, igamma, nblend, ngamma
      real(wp) :: analytic(ndim, ndim), numeric(ndim, ndim), point(ndim), shifted(ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim), packed_h(6)
      type(mctc_error), allocatable :: lsf_err

      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 2, kind, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call lsf%f012_r(lsf2_rr=packed_h)
                  call unpack_sym_matrix(packed_h, analytic, ndim)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f012_r(lsf1_r=g_pp)
                     shifted = point; shifted(axis) = point(axis) + STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f012_r(lsf1_r=g_p)
                     shifted = point; shifted(axis) = point(axis) - STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f012_r(lsf1_r=g_m)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*STEP_SIZE
                     call lsf%prepare(shifted, lsf_err); call lsf%f012_r(lsf1_r=g_mm)
                     do i = 1, ndim
                        call fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), STEP_SIZE, numeric(i, axis), error)
                        if (allocated(error)) return
                     end do
                  end do
                  do j = 1, ndim
                     do i = 1, ndim
                        call check(error, analytic(i, j), numeric(i, j), &
                                   thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                        if (allocated(error)) return
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
      end do
   end subroutine run_f2_rr_fd

   !> SvdW dispatch for the spatial third-derivative FD check
   subroutine test_svdw_f3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rrr_fd(error, "svdw")
   end subroutine test_svdw_f3_rrr_fd

   !> CFC dispatch for the spatial third-derivative FD check
   subroutine test_cfc_f3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rrr_fd(error, "cfc")
   end subroutine test_cfc_f3_rrr_fd

   !> SSES dispatch for the spatial third-derivative FD check
   subroutine test_sses_f3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rrr_fd(error, "sses")
   end subroutine test_sses_f3_rrr_fd

   !> Third spatial derivative vs FD of Hessian
   !>
   !> Both branches use `f3_rrr` to share the internal path at -O3
   subroutine run_f3_rrr_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, axis, i, j, k, iblend, igamma, nblend, ngamma
      real(wp) :: analytic(ndim, ndim, ndim), packed_h(6), packed_t(10)
      real(wp) :: numeric(ndim, ndim, ndim), point(ndim), shifted(ndim)
      real(wp) :: hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call lsf%f3_rrr(lsf3_rrr=packed_t)
                  call unpack_sym_tensor3(packed_t, analytic, ndim)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=packed_h, lsf3_rrr=packed_t)
                     call unpack_sym_matrix(packed_h, hess_pp, ndim)
                     shifted = point; shifted(axis) = point(axis) + eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=packed_h, lsf3_rrr=packed_t)
                     call unpack_sym_matrix(packed_h, hess_p, ndim)
                     shifted = point; shifted(axis) = point(axis) - eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=packed_h, lsf3_rrr=packed_t)
                     call unpack_sym_matrix(packed_h, hess_m, ndim)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=packed_h, lsf3_rrr=packed_t)
                     call unpack_sym_matrix(packed_h, hess_mm, ndim)
                     do j = 1, ndim
                        do i = 1, ndim
                           call fd4_scalar(hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), eps, &
                                           numeric(i, j, axis), error)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
                  do k = 1, ndim
                     do j = 1, ndim
                        do i = 1, ndim
                           call check(error, analytic(i, j, k), numeric(i, j, k), &
                                      thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
      end do
   end subroutine run_f3_rrr_fd

   !* ================================================================================= *!
   !*                          Nuclear-derivative FD tests                              *!
   !* ================================================================================= *!

   !> SvdW dispatch for the nuclear gradient FD check
   subroutine test_svdw_f1_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rA_fd(error, "svdw")
   end subroutine test_svdw_f1_rA_fd

   !> CFC dispatch for the nuclear gradient FD check
   subroutine test_cfc_f1_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rA_fd(error, "cfc")
   end subroutine test_cfc_f1_rA_fd

   !> SSES dispatch for the nuclear gradient FD check
   subroutine test_sses_f1_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rA_fd(error, "sses")
   end subroutine test_sses_f1_rA_fd

   !> Nuclear-position gradient vs FD of f0
   subroutine run_f1_rA_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, atom, axis, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :)
      real(wp), allocatable :: dummy_3rd(:, :, :)
      real(wp) :: numeric, f_pp, f_p, f_m, f_mm
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         allocate (dummy_3rd(6, ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate (analytic(ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_ssd(lsf, centers_base, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f3_rr_rA(lsf1_rA=analytic, lsf3_rr_rA=dummy_3rd)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f0(f_pp)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f0(f_p)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f0(f_m)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f0(f_mm)
                        call fd4_scalar(f_pp, f_p, f_m, f_mm, eps, numeric, error)
                        if (allocated(error)) return
                        call check(error, analytic(axis, atom), numeric, &
                                   thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                        if (allocated(error)) return
                     end do
                  end do
               end do
               deallocate (analytic, lsf)
            end do
         end do
         deallocate (centers_base, centers_local, dummy_3rd)
      end do
   end subroutine run_f1_rA_fd

   !> SvdW dispatch for the mixed spatial-nuclear FD check
   subroutine test_svdw_f2_r_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rA_fd(error, "svdw")
   end subroutine test_svdw_f2_r_rA_fd

   !> CFC dispatch for the mixed spatial-nuclear FD check
   subroutine test_cfc_f2_r_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rA_fd(error, "cfc")
   end subroutine test_cfc_f2_r_rA_fd

   !> SSES dispatch for the mixed spatial-nuclear FD check
   subroutine test_sses_f2_r_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rA_fd(error, "sses")
   end subroutine test_sses_f2_r_rA_fd

   !> Mixed spatial-nuclear second derivative vs FD of spatial grad
   subroutine run_f2_r_rA_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, atom, axis, i, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :)
      real(wp), allocatable :: dummy_3rd(:, :, :)
      real(wp) :: numeric, g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         allocate (dummy_3rd(6, ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate (analytic(ndim, ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_ssd(lsf, centers_base, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f3_rr_rA(lsf2_r_rA=analytic, lsf3_rr_rA=dummy_3rd)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf1_r=g_pp)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf1_r=g_p)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf1_r=g_m)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf1_r=g_mm)
                        do i = 1, ndim
                           call fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), eps, numeric, error)
                           if (allocated(error)) return
                           call check(error, analytic(i, axis, atom), numeric, &
                                      thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
               deallocate (analytic, lsf)
            end do
         end do
         deallocate (centers_base, centers_local, dummy_3rd)
      end do
   end subroutine run_f2_r_rA_fd

   !> SvdW dispatch for the mixed third FD check
   subroutine test_svdw_f3_rr_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rA_fd(error, "svdw")
   end subroutine test_svdw_f3_rr_rA_fd

   !> CFC dispatch for the mixed third FD check
   subroutine test_cfc_f3_rr_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rA_fd(error, "cfc")
   end subroutine test_cfc_f3_rr_rA_fd

   !> SSES dispatch for the mixed third FD check
   subroutine test_sses_f3_rr_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rA_fd(error, "sses")
   end subroutine test_sses_f3_rr_rA_fd

   !> Mixed third derivative (Hess w.r.t. spatial coords, grad w.r.t
   !> nuclei) vs FD of spatial Hessian
   subroutine run_f3_rr_rA_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, atom, axis, i, j, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :), packed_a(:, :, :)
      real(wp) :: numeric, hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
      real(wp) :: packed_h(6)
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         allocate (analytic(ndim, ndim, ndim, mol%nat), packed_a(6, ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_ssd(lsf, centers_base, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f3_rr_rA(lsf3_rr_rA=packed_a)
                  call unpack_pair_batch(packed_a, analytic, ndim*mol%nat)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=packed_h)
                        call unpack_sym_matrix(packed_h, hess_pp, ndim)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=packed_h)
                        call unpack_sym_matrix(packed_h, hess_p, ndim)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=packed_h)
                        call unpack_sym_matrix(packed_h, hess_m, ndim)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=packed_h)
                        call unpack_sym_matrix(packed_h, hess_mm, ndim)
                        do j = 1, ndim
                           do i = 1, ndim
                              call fd4_scalar(hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), eps, numeric, &
                                              error)
                              if (allocated(error)) return
                              call check(error, analytic(i, j, axis, atom), numeric, &
                                         thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
         deallocate (centers_base, centers_local, analytic, packed_a)
      end do
   end subroutine run_f3_rr_rA_fd

   !* ================================================================================= *!
   !*                          Radius-derivative FD tests                               *!
   !* ================================================================================= *!
   !
   ! `R_a` denotes radius; `f3_rr_rad` rank equals spatial order
   ! Nuclear FDs use `set_centers`; radius FDs need cached-radii `update`
   ! `refresh_radii` restores `max_deriv` after the update
   !
   ! SvdW radius derivative rescales one atom's tensors
   ! CFC radius derivative differentiates a two-center kernel

   !> SvdW dispatch for the radius gradient FD check
   subroutine test_svdw_f1_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rad_fd(error, "svdw")
   end subroutine test_svdw_f1_rad_fd

   !> CFC dispatch for the radius gradient FD check
   subroutine test_cfc_f1_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rad_fd(error, "cfc")
   end subroutine test_cfc_f1_rad_fd

   !> SSES dispatch for the radius gradient FD check
   subroutine test_sses_f1_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rad_fd(error, "sses")
   end subroutine test_sses_f1_rad_fd

   !> Radius gradient dS/dR_a vs FD of f0
   subroutine run_f1_rad_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      integer  :: icase, ipt, atom, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:)
      real(wp), allocatable :: dummy_2nd(:, :)
      real(wp) :: numeric, f_pp, f_p, f_m, f_mm
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (radii_local(mol%nat))
         allocate (dummy_2nd(6, mol%nat))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate (analytic(mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_radii(lsf, mol, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f3_rr_rad(lsf1_rad=analytic, lsf3_rr_rad=dummy_2nd)
                  do atom = 1, mol%nat
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f0(f_pp)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f0(f_p)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f0(f_m)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f0(f_mm)
                     call fd4_scalar(f_pp, f_p, f_m, f_mm, eps, numeric, error)
                     if (allocated(error)) return
                     call check(error, analytic(atom), numeric, &
                                thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                     if (allocated(error)) return
                  end do
               end do
               deallocate (analytic, lsf)
            end do
         end do
         deallocate (radii_local, dummy_2nd)
      end do
   end subroutine run_f1_rad_fd

   !> SvdW dispatch for the mixed spatial-radius FD check
   subroutine test_svdw_f2_r_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rad_fd(error, "svdw")
   end subroutine test_svdw_f2_r_rad_fd

   !> CFC dispatch for the mixed spatial-radius FD check
   subroutine test_cfc_f2_r_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rad_fd(error, "cfc")
   end subroutine test_cfc_f2_r_rad_fd

   !> SSES dispatch for the mixed spatial-radius FD check
   subroutine test_sses_f2_r_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rad_fd(error, "sses")
   end subroutine test_sses_f2_r_rad_fd

   !> Mixed spatial-radius second derivative vs FD of the spatial gradient
   subroutine run_f2_r_rad_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      integer  :: icase, ipt, atom, i, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :)
      real(wp), allocatable :: dummy_2nd(:, :)
      real(wp) :: numeric, g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (radii_local(mol%nat))
         allocate (dummy_2nd(6, mol%nat))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate (analytic(ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_radii(lsf, mol, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f3_rr_rad(lsf2_r_rad=analytic, lsf3_rr_rad=dummy_2nd)
                  do atom = 1, mol%nat
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf1_r=g_pp)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf1_r=g_p)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf1_r=g_m)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf1_r=g_mm)
                     do i = 1, ndim
                        call fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), eps, numeric, error)
                        if (allocated(error)) return
                        call check(error, analytic(i, atom), numeric, &
                                   thr_abs=fd_abs_thr(kind), thr_rel=fd_rel_thr(kind))
                        if (allocated(error)) return
                     end do
                  end do
               end do
               deallocate (analytic, lsf)
            end do
         end do
         deallocate (radii_local, dummy_2nd)
      end do
   end subroutine run_f2_r_rad_fd

   !> SvdW dispatch for the mixed third radius FD check
   subroutine test_svdw_f3_rr_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rad_fd(error, "svdw")
   end subroutine test_svdw_f3_rr_rad_fd

   !> CFC dispatch for the mixed third radius FD check
   subroutine test_cfc_f3_rr_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rad_fd(error, "cfc")
   end subroutine test_cfc_f3_rr_rad_fd

   !> SSES dispatch for the mixed third radius FD check
   subroutine test_sses_f3_rr_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rad_fd(error, "sses")
   end subroutine test_sses_f3_rr_rad_fd

   !> FD thresholds for deep radius checks
   !>
   !> CFC radius derivatives add quotient terms from `1/R_a`
   !> Nuclear derivatives lack that chain; `f3_rr_rA` clears tighter thresholds
   !>
   !> Deep checks miss tight thresholds by 2-3x at `STEP_SIZE`
   !> Half step clears them, consistent with fourth-order stencil error
   !> `f1_rad` and `f2_r_rad` retain tight thresholds
   !>
   !> @param[in]  kind    LSF kind
   !> @param[out] thr_abs Absolute threshold
   !> @param[out] thr_rel Relative threshold
   subroutine radius_fd_thr(kind, thr_abs, thr_rel)
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Absolute threshold
      real(wp), intent(out) :: thr_abs
      !> Relative threshold
      real(wp), intent(out) :: thr_rel

      select case (kind)
      case ("cfc")
         thr_abs = CFC_ABS_THR
         thr_rel = CFC_REL_THR
      case ("sses")
         thr_abs = SSES_ABS_THR
         thr_rel = SSES_REL_THR
      case default
         thr_abs = ABS_THR
         thr_rel = REL_THR
      end select
   end subroutine radius_fd_thr

   !> Absolute FD threshold of the spatial and nuclear checks
   !>
   !> @param[in] kind LSF kind
   pure function fd_abs_thr(kind) result(thr)
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Threshold
      real(wp) :: thr

      thr = ABS_THR
      if (kind == "sses") thr = SSES_ABS_THR
   end function fd_abs_thr

   !> Relative FD threshold of the spatial and nuclear checks
   !>
   !> @param[in] kind LSF kind
   pure function fd_rel_thr(kind) result(thr)
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Threshold
      real(wp) :: thr

      thr = REL_THR
      if (kind == "sses") thr = SSES_REL_THR
   end function fd_rel_thr

   !> FD thresholds of the shared high-order checks
   !>
   !> @param[in]  kind    LSF kind
   !> @param[out] thr_abs Absolute threshold
   !> @param[out] thr_rel Relative threshold
   subroutine high_order_thr(kind, thr_abs, thr_rel)
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Absolute threshold
      real(wp), intent(out) :: thr_abs
      !> Relative threshold
      real(wp), intent(out) :: thr_rel

      select case (kind)
      case ("cfc")
         thr_abs = CFC_ABS_THR
         thr_rel = CFC_REL_THR
      case ("sses")
         thr_abs = SSES_ABS_THR
         thr_rel = SSES_REL_THR
      case default
         thr_abs = ABS_THR
         thr_rel = REL_THR
      end select
   end subroutine high_order_thr

   !> Sampling points of the high-order checks
   !>
   !> SvdW keeps the default sample of [[get_test_points]]; CFC and SSES use [[fd_points]]
   !>
   !> @param[in] kind LSF kind
   pure function high_order_points(kind) result(npts)
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Number of sampling points
      integer :: npts

      if (kind == "svdw") then
         npts = n_points_svdw_high_order
      else
         npts = fd_points(kind)
      end if
   end function high_order_points

   !> Number of parameter sets of the high-order sweeps
   !>
   !> SvdW sweeps `blend_k` x (`blend_2b`, `blend_3b`); CFC and SSES use one set
   !>
   !> @param[in] kind LSF kind
   pure function sweep_count(kind) result(n)
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Parameter-set count
      integer :: n
      integer :: nblend, ngamma

      call svdw_sweep_sizes(kind, nblend, ngamma)
      n = nblend*ngamma
   end function sweep_count

   !> Create the LSF of one parameter set of the high-order sweep
   !>
   !> @param[out] lsf       Polymorphic LSF
   !> @param[in]  mol       Bound structure
   !> @param[in]  radii     Per-atom radii
   !> @param[in]  max_deriv Highest derivative order
   !> @param[in]  kind      LSF kind
   !> @param[in]  isweep    Parameter-set index in `1:sweep_count(kind)`
   subroutine init_lsf_sweep(lsf, mol, radii, max_deriv, kind, isweep, fade_enabled)
      !> Polymorphic LSF
      class(moist_cavity_drop_lsf_type), allocatable, intent(out) :: lsf
      !> Bound structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Highest derivative order
      integer, intent(in) :: max_deriv
      !> LSF kind
      character(len=*), intent(in) :: kind
      !> Parameter-set index
      integer, intent(in) :: isweep
      !> Optional SvdW fade switch (default on)
      logical, intent(in), optional :: fade_enabled
      integer :: nblend, ngamma, iblend, igamma

      call svdw_sweep_sizes(kind, nblend, ngamma)
      iblend = (isweep - 1)/ngamma + 1
      igamma = mod(isweep - 1, ngamma) + 1
      call init_lsf(lsf, mol, radii, max_deriv, kind, &
                    blend_k=svdw_sweep_blend(kind, iblend), &
                    blend_2b=svdw_sweep_2b(kind, igamma), &
                    blend_3b=svdw_sweep_gamma(kind, igamma), fade_enabled=fade_enabled)
   end subroutine init_lsf_sweep

   !> Mixed third derivative d^3S/(dr^2 dR_a) vs FD of the spatial Hessian
   subroutine run_f3_rr_rad_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      integer  :: icase, ipt, atom, i, j, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :), packed_a(:, :)
      real(wp) :: numeric
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: packed_h(6)
      real(wp) :: eps, thr_abs, thr_rel
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call radius_fd_thr(kind, thr_abs, thr_rel)
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (radii_local(mol%nat))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate (analytic(ndim, ndim, mol%nat), packed_a(6, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_radii(lsf, mol, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f3_rr_rad(lsf3_rr_rad=packed_a)
                  call unpack_pair_batch(packed_a, analytic, mol%nat)
                  do atom = 1, mol%nat
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=packed_h)
                     call unpack_sym_matrix(packed_h, h_pp, ndim)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=packed_h)
                     call unpack_sym_matrix(packed_h, h_p, ndim)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=packed_h)
                     call unpack_sym_matrix(packed_h, h_m, ndim)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=packed_h)
                     call unpack_sym_matrix(packed_h, h_mm, ndim)
                     do i = 1, ndim
                        do j = 1, ndim
                           call fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), h_mm(i, j), eps, numeric, error)
                           if (allocated(error)) return
                           call check(error, analytic(i, j, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
               deallocate (analytic, packed_a, lsf)
            end do
         end do
         deallocate (radii_local)
      end do
   end subroutine run_f3_rr_rad_fd

   !* ================================================================================= *!
   !*                Uncontracted two-radius and nuclear-radius FD tests                *!
   !* ================================================================================= *!
   !
   ! [[run_radius_contracted_vs_pairwise]] checks one direction pair
   ! Shared errors can cancel under that contraction
   ! Perturb one radius and difference every element of the other slot
   !
   ! Radius perturbation of atom `b` yields:
   !   * `f*_rad` FD: `f*_radrad(.., :, b)`
   !   * `f*_rA` FD: `f*_rA_rad(.., :, b)`
   ! Nuclear-radius slots differ: position first, radius second
   ! Elementwise FD catches transposition hidden by contraction

   !> SvdW dispatch for the two-radius block FD check
   subroutine test_svdw_radrad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radrad_fd(error, "svdw")
   end subroutine test_svdw_radrad_fd

   !> CFC dispatch for the two-radius block FD check
   subroutine test_cfc_radrad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radrad_fd(error, "cfc")
   end subroutine test_cfc_radrad_fd

   !> SSES dispatch for the two-radius block FD check
   subroutine test_sses_radrad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radrad_fd(error, "sses")
   end subroutine test_sses_radrad_fd

   !> `f{2,3,4}_*_radrad` vs a one-radius FD of the whole `f*_rad` ladder
   subroutine run_radrad_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      real(wp), allocatable :: rr2(:, :), rr3(:, :, :), rr4(:, :, :, :)
      real(wp), allocatable :: l1(:, :), l2(:, :, :), l3(:, :, :, :)
      real(wp), allocatable :: packed_rr4(:, :, :), packed_l3(:, :, :)
      integer  :: icase, ipt, atom, other, i, j, istep, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim), numeric, eps, thr_abs, thr_rel
      real(wp), parameter :: steps(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call radius_fd_thr(kind, thr_abs, thr_rel)
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (radii_local(mol%nat))
         allocate (rr2(mol%nat, mol%nat), rr3(ndim, mol%nat, mol%nat), &
                   rr4(ndim, ndim, mol%nat, mol%nat))
         allocate (l1(mol%nat, 4), l2(ndim, mol%nat, 4), l3(ndim, ndim, mol%nat, 4))
         allocate (packed_rr4(6, mol%nat, mol%nat), packed_l3(6, mol%nat, 4))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_radii(lsf, mol, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f2_radrad(rr2)
                  call lsf%f3_r_radrad(rr3)
                  call lsf%f4_rr_radrad(packed_rr4)
                  call unpack_pair_batch(packed_rr4, rr4, mol%nat*mol%nat)
                  do atom = 1, mol%nat
                     do istep = 1, 4
                        radii_local = radii
                        radii_local(atom) = radii_local(atom) + steps(istep)*eps
                        call refresh_radii(lsf, mol, radii_local)
                        if (prepare_failed(lsf, point, error)) return
                        call lsf%f3_rr_rad(lsf1_rad=l1(:, istep), &
                                           lsf2_r_rad=l2(:, :, istep), &
                                           lsf3_rr_rad=packed_l3(:, :, istep))
                        call unpack_pair_batch(packed_l3(:, :, istep), l3(:, :, :, istep), mol%nat)
                     end do
                     do other = 1, mol%nat
                        call fd4_scalar(l1(other, 1), l1(other, 2), l1(other, 3), l1(other, 4), eps, numeric, error)
                        if (allocated(error)) return
                        call check(error, rr2(other, atom), numeric, &
                                   thr_abs=thr_abs, thr_rel=thr_rel)
                        if (allocated(error)) return
                        do i = 1, ndim
                           call fd4_scalar(l2(i, other, 1), l2(i, other, 2), l2(i, other, 3), l2(i, other, 4), eps, &
                                           numeric, error)
                           if (allocated(error)) return
                           call check(error, rr3(i, other, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                           do j = 1, ndim
                              call fd4_scalar(l3(i, j, other, 1), l3(i, j, other, 2), l3(i, j, other, 3), &
                                              l3(i, j, other, 4), eps, numeric, error)
                              if (allocated(error)) return
                              call check(error, rr4(i, j, other, atom), numeric, &
                                         thr_abs=thr_abs, thr_rel=thr_rel)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
         deallocate (radii_local, rr2, rr3, rr4, l1, l2, l3, packed_rr4, packed_l3)
      end do
   end subroutine run_radrad_fd

   !> SvdW dispatch for the nuclear-radius block FD check
   subroutine test_svdw_rA_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rA_rad_fd(error, "svdw")
   end subroutine test_svdw_rA_rad_fd

   !> CFC dispatch for the nuclear-radius block FD check
   subroutine test_cfc_rA_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rA_rad_fd(error, "cfc")
   end subroutine test_cfc_rA_rad_fd

   !> SSES dispatch for the nuclear-radius block FD check
   subroutine test_sses_rA_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rA_rad_fd(error, "sses")
   end subroutine test_sses_rA_rad_fd

   !> `f{2,3,4}_*_rA_rad` vs a one-radius FD of the whole `f*_rA` ladder
   subroutine run_rA_rad_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      real(wp), allocatable :: mx2(:, :, :), mx3(:, :, :, :), mx4(:, :, :, :, :)
      real(wp), allocatable :: n1(:, :, :), n2(:, :, :, :), n3(:, :, :, :, :)
      real(wp), allocatable :: packed_mx4(:, :, :, :), packed_n3(:, :, :, :)
      integer  :: icase, ipt, atom, other, s_ax, i, j, istep
      integer  :: iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim), numeric, eps, thr_abs, thr_rel
      real(wp), parameter :: steps(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call radius_fd_thr(kind, thr_abs, thr_rel)
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (radii_local(mol%nat))
         allocate (mx2(ndim, mol%nat, mol%nat), mx3(ndim, ndim, mol%nat, mol%nat), &
                   mx4(ndim, ndim, ndim, mol%nat, mol%nat))
         allocate (n1(ndim, mol%nat, 4), n2(ndim, ndim, mol%nat, 4), &
                   n3(ndim, ndim, ndim, mol%nat, 4))
         allocate (packed_mx4(6, ndim, mol%nat, mol%nat), packed_n3(6, ndim, mol%nat, 4))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_radii(lsf, mol, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%f2_rA_rad(mx2)
                  call lsf%f3_r_rA_rad(mx3)
                  call lsf%f4_rr_rA_rad(packed_mx4)
                  call unpack_pair_batch(packed_mx4, mx4, ndim*mol%nat*mol%nat)
                  do atom = 1, mol%nat
                     do istep = 1, 4
                        radii_local = radii
                        radii_local(atom) = radii_local(atom) + steps(istep)*eps
                        call refresh_radii(lsf, mol, radii_local)
                        if (prepare_failed(lsf, point, error)) return
                        call lsf%f3_rr_rA(lsf1_rA=n1(:, :, istep), &
                                          lsf2_r_rA=n2(:, :, :, istep), &
                                          lsf3_rr_rA=packed_n3(:, :, :, istep))
                        call unpack_pair_batch(packed_n3(:, :, :, istep), n3(:, :, :, :, istep), &
                                               ndim*mol%nat)
                     end do
                     do other = 1, mol%nat
                        do s_ax = 1, ndim
                           call fd4_scalar(n1(s_ax, other, 1), n1(s_ax, other, 2), n1(s_ax, other, 3), &
                                           n1(s_ax, other, 4), eps, numeric, error)
                           if (allocated(error)) return
                           call check(error, mx2(s_ax, other, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                           do i = 1, ndim
                              call fd4_scalar(n2(i, s_ax, other, 1), n2(i, s_ax, other, 2), n2(i, s_ax, other, 3), &
                                              n2(i, s_ax, other, 4), eps, numeric, error)
                              if (allocated(error)) return
                              call check(error, mx3(i, s_ax, other, atom), numeric, &
                                         thr_abs=thr_abs, thr_rel=thr_rel)
                              if (allocated(error)) return
                              do j = 1, ndim
                                 call fd4_scalar(n3(i, j, s_ax, other, 1), n3(i, j, s_ax, other, 2), &
                                                 n3(i, j, s_ax, other, 3), n3(i, j, s_ax, other, 4), eps, numeric, &
                                                 error)
                                 if (allocated(error)) return
                                 call check(error, mx4(i, j, s_ax, other, atom), numeric, &
                                            thr_abs=thr_abs, thr_rel=thr_rel)
                                 if (allocated(error)) return
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
         deallocate (radii_local, mx2, mx3, mx4, n1, n2, n3, packed_mx4, packed_n3)
      end do
   end subroutine run_rA_rad_fd

   !* ================================================================================= *!
   !*                    Joint position/radius Hessian-vector products                  *!
   !* ================================================================================= *!
   !
   ! For `R = R(X)`, HVP direction combines `(v_B, vr_B)`
   ! Step centers by `eps*v` and radii by `eps*vr`
   ! Independent directions expose missing or swapped contraction terms
   !
   ! Two rows are covered:
   !   * `hvp_f*_rad`  radius row, `d/dR_a` retained  -- FD of `f*_rad`
   !   * `hvp_f*_rA`   nuclear row with `vrad` passed -- FD of `f*_rA`
   ! Both rows cover all four joint-Hessian blocks, including cross-atom radii

   !> Deterministic asymmetric joint direction, max norm 1
   !>
   !> FD error scales as `(max|v|*eps)**4`
   !> Unscaled `max|v| ~ 4` at `nat = 16` yields ~1e-9 relative error
   !> Unit scaling meets shared `ABS_THR`/`REL_THR`
   !>
   !> @param[in]  nat  Number of atoms
   !> @param[out] v    Nuclear displacement directions [3, nat]
   !> @param[out] vrad Radius directions [nat]
   subroutine joint_direction(nat, v, vrad)
      !> Number of atoms
      integer, intent(in) :: nat
      !> Nuclear displacement directions [3, nat]
      real(wp), intent(out) :: v(:, :)
      !> Radius directions [nat]
      real(wp), intent(out) :: vrad(:)

      integer :: ia
      real(wp) :: scal

      do ia = 1, nat
         v(1, ia) = 0.31_wp*real(mod(ia, 7) + 1, wp) - 0.7_wp
         v(2, ia) = -0.17_wp*real(mod(ia, 5) + 1, wp) + 0.4_wp
         v(3, ia) = 0.23_wp*real(mod(ia, 3) + 1, wp)
         ! Independent of `v` to expose mixed terms
         vrad(ia) = 0.19_wp*real(mod(ia, 4) + 1, wp) - 0.43_wp
      end do

      scal = max(maxval(abs(v)), maxval(abs(vrad)))
      if (scal > 0.0_wp) then
         v = v/scal
         vrad = vrad/scal
      end if
   end subroutine joint_direction

   !> SvdW dispatch for the uncontracted radius-block cross-check
   subroutine test_svdw_radius_pairwise(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radius_contracted_vs_pairwise(error, "svdw")
   end subroutine test_svdw_radius_pairwise

   !> CFC dispatch for the uncontracted radius-block cross-check
   subroutine test_cfc_radius_pairwise(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radius_contracted_vs_pairwise(error, "cfc")
   end subroutine test_cfc_radius_pairwise

   !> SSES dispatch for the uncontracted radius-block cross-check
   subroutine test_sses_radius_pairwise(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radius_contracted_vs_pairwise(error, "sses")
   end subroutine test_sses_radius_pairwise

   !> Recover both `hvp_*` rows from uncontracted radius blocks
   !>
   !> Exact joint-Hessian contraction for `R = R(X)`; roundoff agreement
   !> Pin `f*_radrad` and `f*_rA_rad` at `REL_THR`
   !>
   !> Both rows retain asymmetric nuclear-radius atom slots:
   !>
   !>    hvp_f{n}_rad(A)  = sum_B [ v_B . f_rA_rad(:, B, A) + vr_B f_radrad(A, B) ]
   !>    hvp_f{n}_rA(A,s) = sum_B[v_B . f_rArB(s,A,:,B) + vr_B f_rA_rad(s,A,B)]
   !>
   !> Position index on B in radius row; on A in nuclear row
   subroutine run_radius_contracted_vs_pairwise(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :), v(:, :), vrad(:)
      real(wp), allocatable :: rr2(:, :), rr3(:, :, :), rr4(:, :, :, :)
      real(wp), allocatable :: mx2(:, :, :), mx3(:, :, :, :), mx4(:, :, :, :, :)
      real(wp), allocatable :: q2(:, :, :, :), q3(:, :, :, :, :), q4(:, :, :, :, :, :)
      real(wp), allocatable :: r1(:), r2(:, :), r3(:, :, :)
      real(wp), allocatable :: h1(:, :), h2(:, :, :), h3(:, :, :, :)
      real(wp), allocatable :: p_rr4(:, :, :), p_mx4(:, :, :, :), p_q4(:, :, :, :, :)
      real(wp), allocatable :: p_r3(:, :), p_h3(:, :, :)
      integer  :: icase, ipt, iA, iB, s_ax, t_ax, j, k, nat
      real(wp) :: point(ndim), acc
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         nat = mol%nat
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))

         allocate (v(ndim, nat), vrad(nat))
         call joint_direction(nat, v, vrad)

         allocate (rr2(nat, nat), rr3(ndim, nat, nat), rr4(ndim, ndim, nat, nat))
         allocate (mx2(ndim, nat, nat), mx3(ndim, ndim, nat, nat))
         allocate (mx4(ndim, ndim, ndim, nat, nat))
         allocate (q2(ndim, nat, ndim, nat), q3(ndim, ndim, nat, ndim, nat))
         allocate (q4(ndim, ndim, ndim, nat, ndim, nat))
         allocate (r1(nat), r2(ndim, nat), r3(ndim, ndim, nat))
         allocate (h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat))
         allocate (p_rr4(6, nat, nat), p_mx4(6, ndim, nat, nat), p_q4(6, ndim, nat, ndim, nat))
         allocate (p_r3(6, nat), p_h3(6, ndim, nat))

         call init_lsf(lsf, mol, radii, 4, kind)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            call lsf%prepare(point, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if

            call lsf%f2_radrad(rr2)
            call lsf%f3_r_radrad(rr3)
            call lsf%f4_rr_radrad(p_rr4)
            call lsf%f2_rA_rad(mx2)
            call lsf%f3_r_rA_rad(mx3)
            call lsf%f4_rr_rA_rad(p_mx4)
            call lsf%f2_rArB(q2)
            call lsf%f3_r_rArB(q3)
            call lsf%f4_rr_rArB(p_q4)

            call lsf%hvp_f1_rad(v, vrad, r1)
            call lsf%hvp_f2_r_rad(v, vrad, r2)
            call lsf%hvp_f3_rr_rad(v, vrad, p_r3)
            call lsf%hvp_f1_rA(v, h1, vrad)
            call lsf%hvp_f2_r_rA(v, h2, vrad)
            call lsf%hvp_f3_rr_rA(v, p_h3, vrad)
            call unpack_pair_batch(p_rr4, rr4, nat*nat)
            call unpack_pair_batch(p_mx4, mx4, ndim*nat*nat)
            call unpack_pair_batch(p_q4, q4, ndim*nat*ndim*nat)
            call unpack_pair_batch(p_r3, r3, nat)
            call unpack_pair_batch(p_h3, h3, ndim*nat)

            !* ---------------------------- the radius row ---------------------------- *!
            do iA = 1, lsf%active_count()
               acc = 0.0_wp
               do iB = 1, lsf%active_count()
                  acc = acc + vrad(lsf%active_atom(iB))*rr2(iA, iB)
                  do t_ax = 1, ndim
                     acc = acc + v(t_ax, lsf%active_atom(iB))*mx2(t_ax, iB, iA)
                  end do
               end do
               call check(error, r1(iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return

               do j = 1, ndim
                  acc = 0.0_wp
                  do iB = 1, lsf%active_count()
                     acc = acc + vrad(lsf%active_atom(iB))*rr3(j, iA, iB)
                     do t_ax = 1, ndim
                        acc = acc + v(t_ax, lsf%active_atom(iB))*mx3(j, t_ax, iB, iA)
                     end do
                  end do
                  call check(error, r2(j, iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return

                  do k = 1, ndim
                     acc = 0.0_wp
                     do iB = 1, lsf%active_count()
                        acc = acc + vrad(lsf%active_atom(iB))*rr4(j, k, iA, iB)
                        do t_ax = 1, ndim
                           acc = acc + v(t_ax, lsf%active_atom(iB)) &
                                 *mx4(j, k, t_ax, iB, iA)
                        end do
                     end do
                     call check(error, r3(j, k, iA), acc, &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do

            !* ------ the position row, its radius half the nuclear-radius block ------ *!
            do iA = 1, lsf%active_count()
               do s_ax = 1, ndim
                  acc = 0.0_wp
                  do iB = 1, lsf%active_count()
                     acc = acc + vrad(lsf%active_atom(iB))*mx2(s_ax, iA, iB)
                     do t_ax = 1, ndim
                        acc = acc + q2(s_ax, iA, t_ax, iB) &
                              *v(t_ax, lsf%active_atom(iB))
                     end do
                  end do
                  call check(error, h1(s_ax, iA), acc, &
                             thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return

                  do j = 1, ndim
                     acc = 0.0_wp
                     do iB = 1, lsf%active_count()
                        acc = acc + vrad(lsf%active_atom(iB))*mx3(j, s_ax, iA, iB)
                        do t_ax = 1, ndim
                           acc = acc + q3(j, s_ax, iA, t_ax, iB) &
                                 *v(t_ax, lsf%active_atom(iB))
                        end do
                     end do
                     call check(error, h2(j, s_ax, iA), acc, &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return

                     do k = 1, ndim
                        acc = 0.0_wp
                        do iB = 1, lsf%active_count()
                           acc = acc + vrad(lsf%active_atom(iB)) &
                                 *mx4(j, k, s_ax, iA, iB)
                           do t_ax = 1, ndim
                              acc = acc + q4(j, k, s_ax, iA, t_ax, iB) &
                                    *v(t_ax, lsf%active_atom(iB))
                           end do
                        end do
                        call check(error, h3(j, k, s_ax, iA), acc, &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do

         deallocate (lsf, v, vrad, rr2, rr3, rr4, mx2, mx3, mx4, q2, q3, q4, &
                     r1, r2, r3, h1, h2, h3, p_rr4, p_mx4, p_q4, p_r3, p_h3)
      end do
   end subroutine run_radius_contracted_vs_pairwise

   !> SvdW dispatch for the radius-row HVP FD check
   subroutine test_svdw_hvp_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rad_fd(error, "svdw")
   end subroutine test_svdw_hvp_rad_fd

   !> CFC dispatch for the radius-row HVP FD check
   subroutine test_cfc_hvp_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rad_fd(error, "cfc")
   end subroutine test_cfc_hvp_rad_fd

   !> SSES dispatch for the radius-row HVP FD check
   subroutine test_sses_hvp_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rad_fd(error, "sses")
   end subroutine test_sses_hvp_rad_fd

   !> Radius row `hvp_f*_rad` vs joint-direction FD of the `f*_rad` ladder
   subroutine run_hvp_rad_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      real(wp), allocatable :: v(:, :), vrad(:)
      real(wp), allocatable :: a1(:), a2(:, :), a3(:, :, :)
      real(wp), allocatable :: l1(:, :), l2(:, :, :), l3(:, :, :, :)
      real(wp), allocatable :: p_a3(:, :), p_l3(:, :, :)
      integer  :: icase, ipt, atom, i, j, istep, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim), numeric, eps, thr_abs, thr_rel
      real(wp), parameter :: steps(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call radius_fd_thr(kind, thr_abs, thr_rel)
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         allocate (radii_local(mol%nat), v(ndim, mol%nat), vrad(mol%nat))
         allocate (a1(mol%nat), a2(ndim, mol%nat), a3(ndim, ndim, mol%nat))
         allocate (l1(mol%nat, 4), l2(ndim, mol%nat, 4), l3(ndim, ndim, mol%nat, 4))
         allocate (p_a3(6, mol%nat), p_l3(6, mol%nat, 4))
         centers_base = mol%xyz
         call joint_direction(mol%nat, v, vrad)
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_joint(lsf, mol, centers_base, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%hvp_f1_rad(v, vrad, a1)
                  call lsf%hvp_f2_r_rad(v, vrad, a2)
                  call lsf%hvp_f3_rr_rad(v, vrad, p_a3)
                  call unpack_pair_batch(p_a3, a3, mol%nat)
                  do istep = 1, 4
                     centers_local = centers_base + steps(istep)*eps*v
                     radii_local = radii + steps(istep)*eps*vrad
                     call refresh_joint(lsf, mol, centers_local, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f3_rr_rad(lsf1_rad=l1(:, istep), &
                                        lsf2_r_rad=l2(:, :, istep), &
                                        lsf3_rr_rad=p_l3(:, :, istep))
                     call unpack_pair_batch(p_l3(:, :, istep), l3(:, :, :, istep), mol%nat)
                  end do
                  do atom = 1, mol%nat
                     call fd4_scalar(l1(atom, 1), l1(atom, 2), l1(atom, 3), l1(atom, 4), eps, numeric, error)
                     if (allocated(error)) return
                     call check(error, a1(atom), numeric, &
                                thr_abs=thr_abs, thr_rel=thr_rel)
                     if (allocated(error)) return
                     do i = 1, ndim
                        call fd4_scalar(l2(i, atom, 1), l2(i, atom, 2), l2(i, atom, 3), l2(i, atom, 4), eps, numeric, &
                                        error)
                        if (allocated(error)) return
                        call check(error, a2(i, atom), numeric, &
                                   thr_abs=thr_abs, thr_rel=thr_rel)
                        if (allocated(error)) return
                        do j = 1, ndim
                           call fd4_scalar(l3(i, j, atom, 1), l3(i, j, atom, 2), l3(i, j, atom, 3), l3(i, j, atom, 4), &
                                           eps, numeric, error)
                           if (allocated(error)) return
                           call check(error, a3(i, j, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
         deallocate (centers_base, centers_local, radii_local, v, vrad, &
                     a1, a2, a3, l1, l2, l3, p_a3, p_l3)
      end do
   end subroutine run_hvp_rad_fd

   !> SvdW dispatch for the joint nuclear-row HVP FD check
   subroutine test_svdw_hvp_rA_joint_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rA_joint_fd(error, "svdw")
   end subroutine test_svdw_hvp_rA_joint_fd

   !> CFC dispatch for the joint nuclear-row HVP FD check
   subroutine test_cfc_hvp_rA_joint_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rA_joint_fd(error, "cfc")
   end subroutine test_cfc_hvp_rA_joint_fd

   !> SSES dispatch for the joint nuclear-row HVP FD check
   subroutine test_sses_hvp_rA_joint_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rA_joint_fd(error, "sses")
   end subroutine test_sses_hvp_rA_joint_fd

   !> Nuclear row `hvp_f*_rA(v, res, vrad)` vs joint-direction FD of `f*_rA`
   subroutine run_hvp_rA_joint_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), radii_local(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      real(wp), allocatable :: v(:, :), vrad(:)
      real(wp), allocatable :: h1(:, :), h2(:, :, :), h3(:, :, :, :)
      real(wp), allocatable :: n1(:, :, :), n2(:, :, :, :), n3(:, :, :, :, :)
      real(wp), allocatable :: p_h3(:, :, :), p_n3(:, :, :, :)
      integer  :: icase, ipt, atom, s_ax, i, j, istep, iblend, igamma, nblend, ngamma
      real(wp) :: point(ndim), numeric, eps, thr_abs, thr_rel
      real(wp), parameter :: steps(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE
      call radius_fd_thr(kind, thr_abs, thr_rel)
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         allocate (radii_local(mol%nat), v(ndim, mol%nat), vrad(mol%nat))
         allocate (h1(ndim, mol%nat), h2(ndim, ndim, mol%nat), &
                   h3(ndim, ndim, ndim, mol%nat))
         allocate (n1(ndim, mol%nat, 4), n2(ndim, ndim, mol%nat, 4), &
                   n3(ndim, ndim, ndim, mol%nat, 4))
         allocate (p_h3(6, ndim, mol%nat), p_n3(6, ndim, mol%nat, 4))
         centers_base = mol%xyz
         call joint_direction(mol%nat, v, vrad)
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                             fade_enabled=fd_user_space_fade, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_joint(lsf, mol, centers_base, radii)
                  call lsf%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call check_identity_active(error, lsf, mol%nat)
                  if (allocated(error)) return
                  call lsf%hvp_f1_rA(v, h1, vrad)
                  call lsf%hvp_f2_r_rA(v, h2, vrad)
                  call lsf%hvp_f3_rr_rA(v, p_h3, vrad)
                  call unpack_pair_batch(p_h3, h3, ndim*mol%nat)
                  do istep = 1, 4
                     centers_local = centers_base + steps(istep)*eps*v
                     radii_local = radii + steps(istep)*eps*vrad
                     call refresh_joint(lsf, mol, centers_local, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f3_rr_rA(lsf1_rA=n1(:, :, istep), &
                                       lsf2_r_rA=n2(:, :, :, istep), &
                                       lsf3_rr_rA=p_n3(:, :, :, istep))
                     call unpack_pair_batch(p_n3(:, :, :, istep), n3(:, :, :, :, istep), &
                                            ndim*mol%nat)
                  end do
                  do atom = 1, mol%nat
                     do s_ax = 1, ndim
                        call fd4_scalar(n1(s_ax, atom, 1), n1(s_ax, atom, 2), n1(s_ax, atom, 3), n1(s_ax, atom, 4), &
                                        eps, numeric, error)
                        if (allocated(error)) return
                        call check(error, h1(s_ax, atom), numeric, &
                                   thr_abs=thr_abs, thr_rel=thr_rel)
                        if (allocated(error)) return
                        do i = 1, ndim
                           call fd4_scalar(n2(i, s_ax, atom, 1), n2(i, s_ax, atom, 2), n2(i, s_ax, atom, 3), &
                                           n2(i, s_ax, atom, 4), eps, numeric, error)
                           if (allocated(error)) return
                           call check(error, h2(i, s_ax, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                           do j = 1, ndim
                              call fd4_scalar(n3(i, j, s_ax, atom, 1), n3(i, j, s_ax, atom, 2), &
                                              n3(i, j, s_ax, atom, 3), n3(i, j, s_ax, atom, 4), eps, numeric, error)
                              if (allocated(error)) return
                              call check(error, h3(i, j, s_ax, atom), numeric, &
                                         thr_abs=thr_abs, thr_rel=thr_rel)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
         deallocate (centers_base, centers_local, radii_local, v, vrad, &
                     h1, h2, h3, n1, n2, n3, p_h3, p_n3)
      end do
   end subroutine run_hvp_rA_joint_fd

   !* ================================================================================= *!
   !*                           Prepare at an unreached point                           *!
   !* ================================================================================= *!
   !
   ! A point beyond every reach has no value: `prepare` reports it through
   ! `unreached()`, never through an `error` or zeros
   !
   ! Every accessor after an unreached prepare aborts (`error stop`), which the
   ! harness cannot catch in process; only the flag is tested

   !> SvdW dispatch for the unreached-prepare check, plus `Z` underflow in reach
   subroutine test_svdw_empty_active(error)
      type(error_type), allocatable, intent(out) :: error
      call run_empty_active(error, "svdw")
      if (allocated(error)) return
      call check_z_underflow_unreached(error)
      if (allocated(error)) return
      call check_svdw_band_end_unreached(error)
   end subroutine test_svdw_empty_active

   !> CFC dispatch for the unreached-prepare check, plus a tiny `PD` in reach
   subroutine test_cfc_empty_active(error)
      type(error_type), allocatable, intent(out) :: error
      call run_empty_active(error, "cfc")
      if (allocated(error)) return
      call check_cfc_tiny_pd_unreached(error)
   end subroutine test_cfc_empty_active

   !> Prepare at a point with no atom in reach
   !>
   !> - `prepare` and `prepare_subset` (empty list): `unreached()`, no `error`,
   !>   `active_count() = 0`
   !> - negative control: a point near atom 1 afterwards is reached and its
   !>   value finite, so the flag is per prepare
   subroutine run_empty_active(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)
      real(wp) :: point(ndim), val
      integer :: none(0)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)

      ! A point beyond every reach (fade band end, or the threshold for CFC)
      call init_lsf(lsf, mol, radii, 3, kind, screening_threshold=1.0e-10_wp)
      point = [5.0e2_wp, 5.0e2_wp, 5.0e2_wp]

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "prepare must not raise at an unreached point: "//lsf_err%message)
         return
      end if
      call check(error, lsf%unreached(), more="prepare: unreached() must report the point")
      if (allocated(error)) return
      call check(error, lsf%active_count(), 0, more="unreached point must leave no active atom")
      if (allocated(error)) return

      call lsf%prepare_subset(point, none, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "prepare_subset must not raise at an unreached point: "// &
                          lsf_err%message)
         return
      end if
      call check(error, lsf%unreached(), more="prepare_subset: empty list must be unreached")
      if (allocated(error)) return
      call check(error, lsf%active_count(), 0, more="active_count must stay readable")
      if (allocated(error)) return

      call lsf%prepare(mol%xyz(:, 1) + [0.3_wp, 0.2_wp, 0.1_wp], lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "prepare near atom 1 failed: "//lsf_err%message)
         return
      end if
      call check(error, .not. lsf%unreached(), more="a reached prepare must clear the flag")
      if (allocated(error)) return
      call lsf%f0(val)
      call check(error, ieee_is_finite(val), more="reached value must be finite")
   end subroutine run_empty_active

   !> SvdW unreached with atoms in reach: `Z` underflows to zero
   !>
   !> - fade off, screening off: a point ~200 Bohr out keeps every atom active
   !>   (`p1 > 0`) but `Z` underflows; `unreached()` with `active_count() > 0`,
   !>   as the batch reports for the same point, `S = +Inf` there
   subroutine check_z_underflow_unreached(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: point(ndim), b0(1)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.false.)
      call init_svdw_param(lsf, mol, radii, 0.0_wp, param)
      point = lsf%center_upper + 200.0_wp

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "prepare must not raise at an unreached point: "//lsf_err%message)
         return
      end if
      call check(error, lsf%active_count() > 0, &
                 more="fixture must keep atoms in reach; the rest says nothing otherwise")
      if (allocated(error)) return
      call check(error, lsf%unreached(), more="Z underflow must be unreached")
      if (allocated(error)) return

      call lsf%evaluate_batch(reshape(point, [ndim, 1]), 0, b0, error=lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, lsf%unreached(1), more="prepare and batch must agree on unreached")
      if (allocated(error)) return
      call check(error, b0(1) > huge(b0(1)), more="an unreached column must hold S = +Inf")
   end subroutine check_z_underflow_unreached

   !> SvdW unreached by its higher orders: a finite value near the band end
   !>
   !> - one carbon, `R = 3`, default fade; the point at `1 - t = 1e-6` from the
   !>   band end keeps the atom active, `Z` is small but positive, and the
   !>   third and fourth derivatives overflow (powers of `1/Z`)
   !> - `max_deriv = 4`: unreached (`S > 0`) with `active_count() = 1`, no `error`
   !> - `max_deriv = 0`: reached, finite value (negative control: the rule is per
   !>   order, not per point)
   !> - the batch agrees at orders 4 (unreached) and 0 (reached); the order-4 column
   !>   keeps the finite order-0 `S` (no `+Inf`), so `exclusion_radius` of it stays
   !>   below `r_cap` with the fade on
   subroutine check_svdw_band_end_unreached(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: point(ndim), x_on, delta, val
      real(wp) :: b0(1), b1(ndim, 1), b2(6, 1), b3(10, 1), b4(15, 1)
      !> Order-4 batch value of the unreached column
      real(wp) :: s4
      type(mctc_error), allocatable :: lsf_err

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]))
      call param%new(fade_enabled=.true.)
      call init_svdw_param(lsf, mol, [3.0_wp], 0.0_wp, param)
      x_on = 3.0_wp*lsf%param%fade_start/lsf%param%blend_k
      delta = 3.0_wp*lsf%param%fade_width/lsf%param%blend_k
      point = [3.0_wp + x_on + delta*(1.0_wp - 1.0e-6_wp), 0.0_wp, 0.0_wp]

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "band end: prepare must not raise: "//lsf_err%message)
         return
      end if
      call check(error, lsf%active_count(), 1, more="band end: the atom must stay in reach")
      if (allocated(error)) return
      call check(error, lsf%unreached(), more="band end: non-finite orders 3/4 must be unreached")
      if (allocated(error)) return
      call lsf%evaluate_batch(reshape(point, [ndim, 1]), 4, b0, b1, b2, b3, b4, error=lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "band end: evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, lsf%unreached(1), more="band end: the order-4 batch must agree (unreached)")
      if (allocated(error)) return
      s4 = b0(1)
      call check(error, ieee_is_finite(s4) .and. s4 > 0.0_wp, &
                 more="band end: only derivatives overflow, the batch must keep the finite S > 0")
      if (allocated(error)) return

      call lsf%set_max_deriv(0)
      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "band end: prepare at order 0 failed: "//lsf_err%message)
         return
      end if
      call check(error, .not. lsf%unreached(), more="band end: the value alone is reached")
      if (allocated(error)) return
      call lsf%f0(val)
      call check(error, ieee_is_finite(val), more="band end: finite value at order 0")
      if (allocated(error)) return
      call lsf%evaluate_batch(reshape(point, [ndim, 1]), 0, b0, error=lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "band end: order-0 evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, .not. lsf%unreached(1), &
                 more="band end: the order-0 batch must agree (reached)")
      if (allocated(error)) return
      ! A `+Inf` stand-in would certify `r_cap` around a point at finite S
      call check(error, s4, b0(1), rel=.true., thr=1.0e-12_wp, &
                 more="band end: the order-4 unreached value must be the order-0 S")
      if (allocated(error)) return
      call check(error, lsf%exclusion_radius(s4), lsf%exclusion_radius(b0(1)), rel=.true., &
                 thr=1.0e-12_wp, more="band end: exclusion radius of the order-4 value")
   end subroutine check_svdw_band_end_unreached

   !> CFC unreached by its higher orders: a tiny `PD` with both atoms in reach
   !>
   !> - two carbons 100 Bohr apart, `R = 3`, fade off, screening off; at
   !>   `x = 51` `PD` is positive and the value finite, the third and fourth
   !>   derivatives overflow
   !> - `max_deriv = 4`: unreached (`S > 0`) with `active_count() = 2`, as the
   !>   order-4 batch; `max_deriv = 0`: reached and a finite value, as the
   !>   order-0 batch; the order-4 column keeps that finite `S`, not `+Inf`
   subroutine check_cfc_tiny_pd_unreached(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: point(ndim), val
      real(wp) :: b0(1), b1(ndim, 1), b2(6, 1), b3(10, 1), b4(15, 1)
      !> Order-4 batch value of the unreached column
      real(wp) :: s4
      type(mctc_error), allocatable :: lsf_err

      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, 100.0_wp, 0.0_wp, 0.0_wp], [ndim, 2]))
      call param%new(fade_enabled=.false.)
      call init_cfc_param(lsf, mol, [3.0_wp, 3.0_wp], 0.0_wp, param)
      point = [51.0_wp, 0.0_wp, 0.0_wp]

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "tiny PD: prepare must not raise: "//lsf_err%message)
         return
      end if
      call check(error, lsf%active_count(), 2, more="tiny PD: both atoms must stay in reach")
      if (allocated(error)) return
      call check(error, lsf%unreached(), more="tiny PD: non-finite orders 3/4 must be unreached")
      if (allocated(error)) return
      call lsf%evaluate_batch(reshape(point, [ndim, 1]), 4, b0, b1, b2, b3, b4, error=lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "tiny PD: evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, lsf%unreached(1), more="tiny PD: the order-4 batch must agree (unreached)")
      if (allocated(error)) return
      s4 = b0(1)
      call check(error, ieee_is_finite(s4) .and. s4 > 0.0_wp, &
                 more="tiny PD: only derivatives overflow, the batch must keep the finite S > 0")
      if (allocated(error)) return

      call lsf%set_max_deriv(0)
      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "tiny PD: prepare at order 0 failed: "//lsf_err%message)
         return
      end if
      call check(error, .not. lsf%unreached(), more="tiny PD: the value alone is reached")
      if (allocated(error)) return
      call lsf%f0(val)
      call check(error, ieee_is_finite(val), more="tiny PD: finite value at order 0")
      if (allocated(error)) return
      call lsf%evaluate_batch(reshape(point, [ndim, 1]), 0, b0, error=lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "tiny PD: order-0 evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, .not. lsf%unreached(1), &
                 more="tiny PD: the order-0 batch must agree (reached)")
      if (allocated(error)) return
      ! A `+Inf` stand-in would certify `r_cap` around a point at finite S
      call check(error, s4, b0(1), rel=.true., thr=1.0e-12_wp, &
                 more="tiny PD: the order-4 unreached value must be the order-0 S")
      if (allocated(error)) return
      call check(error, lsf%exclusion_radius(s4), lsf%exclusion_radius(b0(1)), rel=.true., &
                 thr=1.0e-12_wp, more="tiny PD: exclusion radius of the order-4 value")
   end subroutine check_cfc_tiny_pd_unreached

   !> Rebind LSF to perturbed positions and radii
   !>
   !> Radius changes require `update` to refresh per-atom caches
   !>
   !> @param[inout] lsf     Polymorphic LSF, rebound in place
   !> @param[in]    mol     Molecular structure supplying everything but xyz
   !> @param[in]    centers Perturbed atom positions (3, mol%nat)
   !> @param[in]    radii   Perturbed per-atom radii (size mol%nat)
   subroutine refresh_joint(lsf, mol, centers, radii)
      !> Polymorphic LSF whose geometry and radii are refreshed in place
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Molecular structure to rebind to
      type(structure_type), intent(in) :: mol
      !> Perturbed atom positions
      real(wp), intent(in) :: centers(:, :)
      !> Perturbed per-atom radii
      real(wp), intent(in) :: radii(:)

      type(structure_type) :: mol_local

      if (size(radii) /= mol%nat .or. size(centers, 2) /= mol%nat) then
         error stop "refresh_joint: radii/centers/structure size mismatch"
      end if
      mol_local = mol
      mol_local%xyz = centers
      call lsf%update(mol_local, radii)
      call lsf%set_max_deriv(3)
   end subroutine refresh_joint

   !> Rebind LSF to perturbed radii
   !>
   !> Radius FDs need full `update` for per-atom caches
   !> Restore `max_deriv` afterward to keep `f012_r` valid
   !>
   !> @param[inout] lsf   Polymorphic LSF, rebound in place
   !> @param[in]    mol   Molecular structure (unchanged)
   !> @param[in]    radii Perturbed per-atom radii (size mol%nat)
   subroutine refresh_radii(lsf, mol, radii)
      !> Polymorphic LSF whose radii are refreshed in place
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Molecular structure to rebind to
      type(structure_type), intent(in) :: mol
      !> Perturbed per-atom radii
      real(wp), intent(in) :: radii(:)

      if (size(radii) /= mol%nat) then
         error stop "refresh_radii: radii/structure size mismatch"
      end if
      call lsf%update(mol, radii)
      call lsf%set_max_deriv(3)
   end subroutine refresh_radii

   !* ================================================================================= *!
   !*                     Candidate cell list versus the full scan                      *!
   !* ================================================================================= *!

   !> Point-box candidates, screened, reproduce the full-scan active set
   !>
   !> - 1000-atom jittered lattice (binned cell list) and a small molecule (single cell)
   !> - fade off at two thresholds, the second set in place by
   !>   `set_screening_threshold`, then the default fade (reach = band end)
   !> - points in the padded enclosure plus far points with no atom and no candidate
   !> - one extended box whose candidates serve every point inside it
   !> - ball candidates equal a brute force over the gate's own reach
   !> - active sets compared as sorted user ids: `prepare` and `prepare_subset`
   !>   visit candidates in different orders
   !> - `unreached()` of both prepares equal; the far points unreached
   subroutine test_svdw_box_candidates(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nside = 10
      integer, parameter :: npoint = 150
      integer, parameter :: nextra = 10
      real(wp), parameter :: spacing = 3.0_wp
      real(wp), parameter :: thresholds(3) = [1.0e-4_wp, 1.0e-10_wp, 1.0e-10_wp]

      type(structure_type) :: mols(2)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      type(mctc_error), allocatable :: lsf_err
      real(wp), allocatable :: radii(:), xyz(:, :)
      integer, allocatable :: num(:), cand(:), full(:), boxed(:), ref(:)
      real(wp) :: point(ndim), u(ndim), lo(ndim), hi(ndim), blo(ndim), bhi(ndim)
      real(wp) :: r, reach, d2
      integer :: state
      integer :: imol, ithr, ip, i, j, k, n, ncand, nfull, nbox, nonempty, pruned, nref
      !> Unreached flag of the full-scan and the box-list prepare
      logical :: un_full, un_box

      allocate (num(nside**3), xyz(ndim, nside**3))
      state = 41
      n = 0
      do k = 1, nside
         do j = 1, nside
            do i = 1, nside
               n = n + 1
               call det_rand(u, state)
               xyz(:, n) = spacing*(real([i, j, k], wp) + 0.4_wp*(u - 0.5_wp))
               num(n) = merge(6, 1, mod(n, 3) == 0)
            end do
         end do
      end do
      call new(mols(1), num, xyz)
      call get_structure(mols(2), "MB16-43", "01")

      do imol = 1, size(mols)
         call get_test_radii(mols(imol), radii)
         call init_lsf(lsf, mols(imol), radii, 0, "svdw", screening_threshold=thresholds(1), &
                       fade_enabled=.false.)
         do ithr = 1, size(thresholds)
            ! In-place rebuild of the candidate cell list (fade off), then the fade
            ! on, whose reach is the band end
            if (ithr == 2) call lsf%set_screening_threshold(thresholds(ithr))
            if (ithr == 3) call init_lsf(lsf, mols(imol), radii, 0, "svdw", &
                                         screening_threshold=thresholds(ithr))
            call check(error, lsf%cand_cells%full_scan, imol == 2, &
                       more="cell list mode must follow the 200-candidate switch")
            if (allocated(error)) return

            allocate (cand(lsf%ncenters), full(lsf%ncenters), boxed(lsf%ncenters), &
                      ref(lsf%ncenters))
            lo = lsf%center_lower - 4.0_wp
            hi = lsf%center_upper + 4.0_wp
            nonempty = 0
            pruned = 0
            do ip = 1, npoint
               call det_rand(u, state)
               point = lo + (hi - lo)*u
               if (ip > npoint - 10) point = hi + 40.0_wp*u + 20.0_wp

               call lsf%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "prepare failed: "//lsf_err%message)
                  return
               end if
               un_full = lsf%unreached()
               nfull = lsf%active_count()
               full(1:nfull) = sorted_ids(lsf, nfull)

               call lsf%candidates_in_box(point, point, cand, ncand)
               call lsf%prepare_subset(point, cand(1:ncand), lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "prepare_subset failed: "//lsf_err%message)
                  return
               end if
               un_box = lsf%unreached()
               nbox = lsf%active_count()
               boxed(1:nbox) = sorted_ids(lsf, nbox)

               call check(error, ncand >= nbox, more="box list must be a superset")
               if (allocated(error)) return
               call check(error, nbox, nfull, more="active count, box list vs full scan")
               if (allocated(error)) return
               call check(error, un_box .eqv. un_full, more="unreached, box list vs full scan")
               if (allocated(error)) return
               if (any(boxed(1:nbox) /= full(1:nfull))) then
                  call test_failed(error, "active atoms differ, box list vs full scan")
                  return
               end if
               if (ip > npoint - 10 .and. (nfull /= 0 .or. ncand /= 0 .or. .not. un_full)) then
                  call test_failed(error, "far point must have no active atom, no candidate "// &
                                   "and be unreached")
                  return
               end if
               if (nfull > 0 .and. nfull < lsf%ncenters) nonempty = nonempty + 1
               if (ncand < lsf%ncenters) pruned = pruned + 1
            end do
            call check(error, nonempty > 0, &
                       more="fixture must screen some but not all atoms somewhere")
            if (allocated(error)) return
            if (imol == 1) then
               call check(error, lsf%cand_cells%ncells > 1 .and. pruned > 0, &
                          more="lattice box queries must prune candidates")
               if (allocated(error)) return
            end if

            ! Extended box around the middle of the atoms, candidates reused inside it
            blo = 0.5_wp*(lsf%center_lower + lsf%center_upper) &
                  - 0.1_wp*(lsf%center_upper - lsf%center_lower) - 0.5_wp
            bhi = 0.5_wp*(lsf%center_lower + lsf%center_upper) &
                  + 0.1_wp*(lsf%center_upper - lsf%center_lower) + 0.5_wp
            call lsf%candidates_in_box(blo, bhi, cand, ncand)
            if (imol == 1) then
               call check(error, ncand > 0 .and. ncand < lsf%ncenters, &
                          more="extended box must prune the lattice")
               if (allocated(error)) return
            end if
            do ip = 1, nextra
               call det_rand(u, state)
               point = blo + (bhi - blo)*u

               call lsf%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "prepare failed: "//lsf_err%message)
                  return
               end if
               un_full = lsf%unreached()
               nfull = lsf%active_count()
               full(1:nfull) = sorted_ids(lsf, nfull)

               call lsf%prepare_subset(point, cand(1:ncand), lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "prepare_subset failed: "//lsf_err%message)
                  return
               end if
               un_box = lsf%unreached()
               nbox = lsf%active_count()
               boxed(1:nbox) = sorted_ids(lsf, nbox)

               call check(error, nbox, nfull, more="active count, extended box vs full scan")
               if (allocated(error)) return
               call check(error, un_box .eqv. un_full, more="unreached, extended box vs full scan")
               if (allocated(error)) return
               if (any(boxed(1:nbox) /= full(1:nfull))) then
                  call test_failed(error, "active atoms differ, extended box vs full scan")
                  return
               end if
            end do

            ! Ball candidates against a brute force over the reach handed to the list
            do ip = 1, nextra
               call det_rand(u, state)
               point = lo + (hi - lo)*u
               r = 0.5_wp*real(ip - 1, wp)
               call lsf%candidates_in_sphere(point, r, cand, ncand)
               nref = 0
               do i = 1, size(lsf%cand_radii)
                  reach = lsf%cand_radii(i) + lsf%screening_offset(lsf%cand_radii(i))
                  d2 = sum((lsf%cand_screen(1:ndim, i) - point)**2)
                  if (d2 <= (r + reach)*(r + reach)) then
                     nref = nref + 1
                     ref(nref) = i
                  end if
               end do
               call check(error, ncand, nref, more="ball candidates vs brute force: count")
               if (allocated(error)) return
               if (any(cand(1:ncand) /= ref(1:nref))) then
                  call test_failed(error, "ball candidates differ from brute force")
                  return
               end if
            end do
            deallocate (cand, full, boxed, ref)
         end do
      end do
   end subroutine test_svdw_box_candidates

   !> Active user ids after the latest prepare, ascending
   function sorted_ids(lsf, n) result(ids)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Active count
      integer, intent(in) :: n
      !> Ascending active ids
      integer :: ids(n)

      integer :: i, j, tmp

      do i = 1, n
         ids(i) = lsf%active_atom(i)
      end do
      do i = 2, n
         tmp = ids(i)
         j = i - 1
         do while (j >= 1)
            if (ids(j) <= tmp) exit
            ids(j + 1) = ids(j)
            j = j - 1
         end do
         ids(j + 1) = tmp
      end do
   end function sorted_ids

   !* ================================================================================= *!
   !*                       Hessian-free value+gradient path                            *!
   !* ================================================================================= *!

   !> SvdW dispatch for the Hessian-free value+gradient consistency check
   subroutine test_svdw_f012_hessfree(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f012_hessfree(error, "svdw")
   end subroutine test_svdw_f012_hessfree

   !> CFC dispatch for the Hessian-free value+gradient consistency check
   subroutine test_cfc_f012_hessfree(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f012_hessfree(error, "cfc")
   end subroutine test_cfc_f012_hessfree

   !> SSES dispatch for the Hessian-free value+gradient consistency check
   subroutine test_sses_f012_hessfree(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f012_hessfree(error, "sses")
   end subroutine test_sses_f012_hessfree

   !> `f012_r` skips the Hessian pass when `lsf2_rr` is not requested
   subroutine run_f012_hessfree(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, i, iblend, igamma, nblend, ngamma
      real(wp) :: v0, v1, g0(ndim), g1(ndim), h(6)
      type(mctc_error), allocatable :: lsf_err

      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, fd_points(kind))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 2, kind, &
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  call lsf%prepare(points(:, ipt), lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  ! Hessian skipped (no lsf2_rr) vs Hessian computed
                  !
                  ! SvdW branches use different expression associations
                  ! Value and gradient agree within a few ULP (~1e-14 relative)
                  call lsf%f012_r(lsf0=v0, lsf1_r=g0)
                  call lsf%f012_r(lsf0=v1, lsf1_r=g1, lsf2_rr=h)
                  call check(error, v0, v1, thr_abs=HESSFREE_ABS, thr_rel=HESSFREE_REL, &
                             more="value changed when the Hessian is skipped")
                  if (allocated(error)) return
                  do i = 1, ndim
                     call check(error, g0(i), g1(i), thr_abs=HESSFREE_ABS, &
                                thr_rel=HESSFREE_REL, &
                                more="gradient changed when the Hessian is skipped")
                     if (allocated(error)) return
                  end do
               end do
               deallocate (lsf)
            end do
         end do
      end do
   end subroutine run_f012_hessfree

   !* ================================================================================= *!
   !*                       SvdW-only sweep helpers                                     *!
   !* ================================================================================= *!

   !> Sampling points per structure for one concrete
   !>
   !> CFC count reflects pair-kernel cost; see [[n_points_cfc]]
   pure function fd_points(kind) result(npts)
      !> Concrete kind selector
      character(len=*), intent(in) :: kind
      !> Number of sampling points to request from `get_test_points`
      integer :: npts
      if (kind == "cfc") then
         npts = n_points_cfc
      else if (kind == "sses") then
         npts = n_points_sses
      else
         npts = n_points_svdw
      end if
   end function fd_points

   pure subroutine svdw_sweep_sizes(kind, nblend, ngamma)
      !> Concrete kind selector
      character(len=*), intent(in) :: kind
      !> Number of blend_k values to sweep
      integer, intent(out) :: nblend
      !> Number of gamma values to sweep
      integer, intent(out) :: ngamma
      if (kind == "svdw") then
         nblend = n_svdw_blends
         ngamma = n_svdw_gammas
      else
         nblend = 1
         ngamma = 1
      end if
   end subroutine svdw_sweep_sizes

   !> SvdW blend_k value, or `huge()` sentinel for CFC
   pure function svdw_sweep_blend(kind, i) result(val)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: i
      real(wp) :: val
      if (kind == "svdw") then
         val = svdw_blend_k_values(i)
      else
         ! CFC placeholder for optional argument wrapper
         val = 0.0_wp
      end if
   end function svdw_sweep_blend

   !> i-th SvdW 2-body weight, or 0 for CFC
   pure function svdw_sweep_2b(kind, i) result(val)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: i
      real(wp) :: val
      if (kind == "svdw") then
         val = svdw_2b_values(i)
      else
         val = 0.0_wp
      end if
   end function svdw_sweep_2b

   !> i-th SvdW gamma (3-body) value, or 0 for CFC
   pure function svdw_sweep_gamma(kind, i) result(val)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: i
      real(wp) :: val
      if (kind == "svdw") then
         val = svdw_gamma_values(i)
      else
         val = 0.0_wp
      end if
   end function svdw_sweep_gamma

   !* ================================================================================= *!
   !*                       High-order derivative tests, all LSFs                       *!
   !* ================================================================================= *!
   !
   ! One runner per quantity, dispatched for SvdW, CFC and SSES
   ! SvdW sweeps `blend_k` x (`blend_2b`, `blend_3b`); CFC and SSES use one set
   ! Thresholds from [[high_order_thr]], points from [[high_order_points]]
   ! CFC exponents (-15, -9) increase order-4 FD truncation error
   ! Nuclear outputs are active-indexed; threshold 0 and identity guards align indices

   !> SvdW dispatch of [[run_f2_rArB]]
   subroutine test_svdw_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rArB(error, "svdw")
   end subroutine test_svdw_f2_rArB

   !> CFC dispatch of [[run_f2_rArB]]
   subroutine test_cfc_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rArB(error, "cfc")
   end subroutine test_cfc_f2_rArB

   !> SSES dispatch of [[run_f2_rArB]]
   subroutine test_sses_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rArB(error, "sses")
   end subroutine test_sses_f2_rArB

   !> Nuclear-position second derivative `f2_rArB` vs FD of `f1_rA`
   subroutine run_f2_rArB(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp) :: hi_abs, hi_rel
      integer :: isweep
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, atomA, axisA, atomB, axisB
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :)
      real(wp), allocatable :: rA_fwd(:, :), rA_fwd2(:, :)
      real(wp), allocatable :: rA_bwd(:, :), rA_bwd2(:, :)
      real(wp), allocatable :: dummy_rr_rA(:, :, :)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(rA_fwd)) deallocate (rA_fwd, rA_fwd2, rA_bwd, rA_bwd2)
         allocate (rA_fwd(ndim, mol%nat), rA_fwd2(ndim, mol%nat))
         allocate (rA_bwd(ndim, mol%nat), rA_bwd2(ndim, mol%nat))
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         if (allocated(analytic)) deallocate (analytic)
         allocate (dummy_rr_rA(6, ndim, mol%nat))
         allocate (analytic(ndim, mol%nat, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 3, kind, isweep, &
                                fade_enabled=fd_user_space_fade)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f2_rArB(analytic)
               ! Check active shape and identity before user-ID reads
               call check(error, size(analytic, 2), prim%active_count())
               if (allocated(error)) return
               call check(error, prim%active_count(), mol%nat)
               if (allocated(error)) return
               do atomA = 1, mol%nat
                  call check(error, prim%active_atom(atomA), atomA, &
                             more="unscreened active list must be the identity")
                  if (allocated(error)) return
               end do
               do atomB = 1, mol%nat
                  do axisB = 1, ndim
                     centers_local = centers_base
                     centers_local(axisB, atomB) = centers_local(axisB, atomB) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf1_rA=rA_fwd, lsf3_rr_rA=dummy_rr_rA)
                     centers_local = centers_base
                     centers_local(axisB, atomB) = centers_local(axisB, atomB) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf1_rA=rA_fwd2, lsf3_rr_rA=dummy_rr_rA)
                     centers_local = centers_base
                     centers_local(axisB, atomB) = centers_local(axisB, atomB) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf1_rA=rA_bwd, lsf3_rr_rA=dummy_rr_rA)
                     centers_local = centers_base
                     centers_local(axisB, atomB) = centers_local(axisB, atomB) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf1_rA=rA_bwd2, lsf3_rr_rA=dummy_rr_rA)
                     do atomA = 1, mol%nat
                        do axisA = 1, ndim
                           call fd4_scalar(rA_fwd2(axisA, atomA), rA_fwd(axisA, atomA), rA_bwd(axisA, atomA), &
                                           rA_bwd2(axisA, atomA), STEP_SIZE, numeric, error)
                           if (allocated(error)) return
                           call check(error, analytic(axisA, atomA, axisB, atomB), numeric, &
                                      thr_abs=hi_abs, thr_rel=hi_rel)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine run_f2_rArB

   !> SvdW dispatch of [[run_f3_r_rArB]]
   subroutine test_svdw_f3_r_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_r_rArB(error, "svdw")
   end subroutine test_svdw_f3_r_rArB

   !> CFC dispatch of [[run_f3_r_rArB]]
   subroutine test_cfc_f3_r_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_r_rArB(error, "cfc")
   end subroutine test_cfc_f3_r_rArB

   !> SSES dispatch of [[run_f3_r_rArB]]
   subroutine test_sses_f3_r_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_r_rArB(error, "sses")
   end subroutine test_sses_f3_r_rArB

   !> Mixed third derivative `f3_r_rArB` vs FD of `f2_r_rA`
   subroutine run_f3_r_rArB(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp) :: hi_abs, hi_rel
      integer :: isweep
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iA, iB, axisA, axisB, jdir
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :)
      real(wp), allocatable :: r_rA_fwd(:, :, :), r_rA_fwd2(:, :, :)
      real(wp), allocatable :: r_rA_bwd(:, :, :), r_rA_bwd2(:, :, :)
      real(wp), allocatable :: dummy_rr_rA(:, :, :)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(r_rA_fwd)) deallocate (r_rA_fwd, r_rA_fwd2, r_rA_bwd, r_rA_bwd2)
         allocate (r_rA_fwd(ndim, ndim, mol%nat), r_rA_fwd2(ndim, ndim, mol%nat))
         allocate (r_rA_bwd(ndim, ndim, mol%nat), r_rA_bwd2(ndim, ndim, mol%nat))
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         if (allocated(analytic)) deallocate (analytic)
         allocate (dummy_rr_rA(6, ndim, mol%nat))
         allocate (analytic(ndim, ndim, mol%nat, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 3, kind, isweep, &
                                fade_enabled=fd_user_space_fade)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f3_r_rArB(analytic)
               do iB = 1, mol%nat
                  do axisB = 1, ndim
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf2_r_rA=r_rA_fwd, lsf3_rr_rA=dummy_rr_rA)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf2_r_rA=r_rA_fwd2, lsf3_rr_rA=dummy_rr_rA)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf2_r_rA=r_rA_bwd, lsf3_rr_rA=dummy_rr_rA)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf2_r_rA=r_rA_bwd2, lsf3_rr_rA=dummy_rr_rA)
                     do iA = 1, mol%nat
                        do axisA = 1, ndim
                           do jdir = 1, ndim
                              call fd4_scalar(r_rA_fwd2(jdir, axisA, iA), r_rA_fwd(jdir, axisA, iA), &
                                              r_rA_bwd(jdir, axisA, iA), r_rA_bwd2(jdir, axisA, iA), STEP_SIZE, numeric, &
                                              error)
                              if (allocated(error)) return
                              call check(error, analytic(jdir, axisA, iA, axisB, iB), numeric, &
                                         thr_abs=hi_abs, thr_rel=hi_rel)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine run_f3_r_rArB

   !> SvdW dispatch of [[run_f4_rrrr]]
   subroutine test_svdw_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rrrr(error, "svdw")
   end subroutine test_svdw_f4_rrrr

   !> CFC dispatch of [[run_f4_rrrr]]
   subroutine test_cfc_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rrrr(error, "cfc")
   end subroutine test_cfc_f4_rrrr

   !> SSES dispatch of [[run_f4_rrrr]]
   subroutine test_sses_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rrrr(error, "sses")
   end subroutine test_sses_f4_rrrr

   !> Pure spatial fourth derivative `f4_rrrr` vs FD of `f3_rrr`
   subroutine run_f4_rrrr(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp) :: hi_abs, hi_rel
      integer :: isweep
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, axis, i, j, kk
      real(wp) :: point(ndim), work_point(ndim)
      real(wp) :: t3_fwd(ndim, ndim, ndim), t3_fwd2(ndim, ndim, ndim)
      real(wp) :: t3_bwd(ndim, ndim, ndim), t3_bwd2(ndim, ndim, ndim)
      real(wp) :: analytic(ndim, ndim, ndim, ndim), packed_t(10), packed_q(15)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))
         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 4, kind, isweep)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f4_rrrr(packed_q)
               call unpack_sym_tensor4(packed_q, analytic, ndim)
               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=packed_t)
                  call unpack_sym_tensor3(packed_t, t3_fwd, ndim)
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=packed_t)
                  call unpack_sym_tensor3(packed_t, t3_fwd2, ndim)
                  work_point = point
                  work_point(axis) = point(axis) - STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=packed_t)
                  call unpack_sym_tensor3(packed_t, t3_bwd, ndim)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=packed_t)
                  call unpack_sym_tensor3(packed_t, t3_bwd2, ndim)
                  do i = 1, ndim
                     do j = 1, ndim
                        do kk = 1, ndim
                           call fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), t3_bwd(i, j, kk), t3_bwd2(i, j, kk), &
                                           STEP_SIZE, numeric, error)
                           if (allocated(error)) return
                           call check(error, analytic(i, j, kk, axis), numeric, &
                                      thr_abs=hi_abs, thr_rel=hi_rel)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine run_f4_rrrr

   !> SvdW dispatch of [[run_f4_rrr_rA]]
   subroutine test_svdw_f4_rrr_rA(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rrr_rA(error, "svdw")
   end subroutine test_svdw_f4_rrr_rA

   !> CFC dispatch of [[run_f4_rrr_rA]]
   subroutine test_cfc_f4_rrr_rA(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rrr_rA(error, "cfc")
   end subroutine test_cfc_f4_rrr_rA

   !> SSES dispatch of [[run_f4_rrr_rA]]
   subroutine test_sses_f4_rrr_rA(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rrr_rA(error, "sses")
   end subroutine test_sses_f4_rrr_rA

   !> Mixed fourth derivative `f4_rrr_rA` (3 r-axes + 1 R-axis) vs FD of `f3_rrr`
   subroutine run_f4_rrr_rA(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp) :: hi_abs, hi_rel
      integer :: isweep
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, atom, axis, i, j, kk
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :), packed_a(:, :, :)
      real(wp) :: packed_t(10)
      real(wp) :: t3_fwd(ndim, ndim, ndim), t3_fwd2(ndim, ndim, ndim)
      real(wp) :: t3_bwd(ndim, ndim, ndim), t3_bwd2(ndim, ndim, ndim)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(analytic)) deallocate (analytic, packed_a)
         allocate (analytic(ndim, ndim, ndim, ndim, mol%nat), packed_a(10, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 4, kind, isweep, &
                                fade_enabled=fd_user_space_fade)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f4_rrr_rA(packed_a)
               call unpack_triple_batch(packed_a, analytic, ndim*mol%nat)
               do atom = 1, mol%nat
                  do axis = 1, ndim
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=packed_t)
                     call unpack_sym_tensor3(packed_t, t3_fwd, ndim)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=packed_t)
                     call unpack_sym_tensor3(packed_t, t3_fwd2, ndim)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=packed_t)
                     call unpack_sym_tensor3(packed_t, t3_bwd, ndim)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=packed_t)
                     call unpack_sym_tensor3(packed_t, t3_bwd2, ndim)
                     do i = 1, ndim
                        do j = 1, ndim
                           do kk = 1, ndim
                              call fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), t3_bwd(i, j, kk), t3_bwd2(i, j, kk), &
                                              STEP_SIZE, numeric, error)
                              if (allocated(error)) return
                              call check(error, analytic(i, j, kk, axis, atom), numeric, &
                                         thr_abs=hi_abs, thr_rel=hi_rel)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine run_f4_rrr_rA

   !> SvdW dispatch of [[run_f4_rr_rArB]]
   subroutine test_svdw_f4_rr_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rr_rArB(error, "svdw")
   end subroutine test_svdw_f4_rr_rArB

   !> CFC dispatch of [[run_f4_rr_rArB]]
   subroutine test_cfc_f4_rr_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rr_rArB(error, "cfc")
   end subroutine test_cfc_f4_rr_rArB

   !> SSES dispatch of [[run_f4_rr_rArB]]
   subroutine test_sses_f4_rr_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f4_rr_rArB(error, "sses")
   end subroutine test_sses_f4_rr_rArB

   !> Mixed fourth derivative `f4_rr_rArB` (2 r-axes + 2 R-axes) vs FD of `f3_rr_rA`
   subroutine run_f4_rr_rArB(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp) :: hi_abs, hi_rel
      integer :: isweep
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iA, iB, axisA, axisB, j, kk
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :, :)
      real(wp), allocatable :: rr_rA_fwd(:, :, :, :), rr_rA_fwd2(:, :, :, :)
      real(wp), allocatable :: rr_rA_bwd(:, :, :, :), rr_rA_bwd2(:, :, :, :)
      real(wp), allocatable :: packed_a(:, :, :, :, :), packed_rr_rA(:, :, :)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(analytic)) deallocate (analytic, packed_a, packed_rr_rA)
         if (allocated(rr_rA_fwd)) deallocate (rr_rA_fwd, rr_rA_fwd2, rr_rA_bwd, rr_rA_bwd2)
         allocate (analytic(ndim, ndim, ndim, mol%nat, ndim, mol%nat))
         allocate (packed_a(6, ndim, mol%nat, ndim, mol%nat), packed_rr_rA(6, ndim, mol%nat))
         allocate (rr_rA_fwd(ndim, ndim, ndim, mol%nat), rr_rA_fwd2(ndim, ndim, ndim, mol%nat))
         allocate (rr_rA_bwd(ndim, ndim, ndim, mol%nat), rr_rA_bwd2(ndim, ndim, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 4, kind, isweep, &
                                fade_enabled=fd_user_space_fade)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f4_rr_rArB(packed_a)
               call unpack_pair_batch(packed_a, analytic, ndim*mol%nat*ndim*mol%nat)
               do iB = 1, mol%nat
                  do axisB = 1, ndim
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=packed_rr_rA)
                     call unpack_pair_batch(packed_rr_rA, rr_rA_fwd, ndim*mol%nat)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=packed_rr_rA)
                     call unpack_pair_batch(packed_rr_rA, rr_rA_fwd2, ndim*mol%nat)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=packed_rr_rA)
                     call unpack_pair_batch(packed_rr_rA, rr_rA_bwd, ndim*mol%nat)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=packed_rr_rA)
                     call unpack_pair_batch(packed_rr_rA, rr_rA_bwd2, ndim*mol%nat)
                     do iA = 1, mol%nat
                        do axisA = 1, ndim
                           do j = 1, ndim
                              do kk = 1, ndim
                                 call fd4_scalar(rr_rA_fwd2(j, kk, axisA, iA), rr_rA_fwd(j, kk, axisA, iA), &
                                                 rr_rA_bwd(j, kk, axisA, iA), rr_rA_bwd2(j, kk, axisA, iA), STEP_SIZE, &
                                                 numeric, error)
                                 if (allocated(error)) return
                                 call check(error, analytic(j, kk, axisA, iA, axisB, iB), &
                                            numeric, thr_abs=hi_abs, thr_rel=hi_rel)
                                 if (allocated(error)) return
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine run_f4_rr_rArB

   !> SvdW dispatch of [[run_normalized_f1_rA]]
   subroutine test_svdw_normalized_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      call run_normalized_f1_rA(error, "svdw")
   end subroutine test_svdw_normalized_f1_rA

   !> CFC dispatch of [[run_normalized_f1_rA]]
   subroutine test_cfc_normalized_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      call run_normalized_f1_rA(error, "cfc")
   end subroutine test_cfc_normalized_f1_rA

   !> SSES dispatch of [[run_normalized_f1_rA]]
   subroutine test_sses_normalized_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      call run_normalized_f1_rA(error, "sses")
   end subroutine test_sses_normalized_f1_rA

   !> Normalized-LSF nuclear gradient vs FD of the normalized value
   subroutine run_normalized_f1_rA(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp) :: hi_abs, hi_rel
      integer :: isweep
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      real(wp), allocatable :: deriv_rA(:, :)
      integer  :: icase, ipt, atom, axis
      real(wp) :: point(ndim)
      real(wp) :: normalized_val
      real(wp) :: f_forward, f_backward, f_forward2, f_backward2, numeric
      type(mctc_error), allocatable :: lsf_err

      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(deriv_rA)) deallocate (deriv_rA)
         allocate (deriv_rA(ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 2, kind, isweep, &
                                fade_enabled=fd_user_space_fade)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%normalized_f01_rA(normalized_val, deriv_rA=deriv_rA)
               do atom = 1, mol%nat
                  do axis = 1, ndim
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%normalized_f01_rA(f_forward)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%normalized_f01_rA(f_forward2)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%normalized_f01_rA(f_backward)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%normalized_f01_rA(f_backward2)
                     call fd4_scalar(f_forward2, f_forward, f_backward, f_backward2, STEP_SIZE, numeric, error)
                     if (allocated(error)) return
                     call check(error, deriv_rA(axis, atom), numeric, &
                                thr_abs=hi_abs, thr_rel=hi_rel)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine run_normalized_f1_rA

   !> SvdW dispatch of [[run_contracted]]
   subroutine test_svdw_contracted(error)
      type(error_type), allocatable, intent(out) :: error
      call run_contracted(error, "svdw")
   end subroutine test_svdw_contracted

   !> CFC dispatch of [[run_contracted]]
   subroutine test_cfc_contracted(error)
      type(error_type), allocatable, intent(out) :: error
      call run_contracted(error, "cfc")
   end subroutine test_cfc_contracted

   !> SSES dispatch of [[run_contracted]]
   subroutine test_sses_contracted(error)
      type(error_type), allocatable, intent(out) :: error
      call run_contracted(error, "sses")
   end subroutine test_sses_contracted

   !> `tangent_*`, `hvp_*` and `vjp_*` against explicit full-tensor contractions
   !>
   !> Independent sibling checks at roundoff, without finite differences
   !>
   !> - `fade`: CFC band `a = 5`, `b = 2` e-folds of `u_A`; points beyond every
   !>   band end (unreached) are skipped, and some reached point must have an atom
   !>   inside its band (negative control: else only the plateau bodies run)
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  LSF kind
   !> @param[in]  fade  CFC biting band 5/2 (optional; default: the default band)
   subroutine run_contracted(error, kind, fade)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      logical, intent(in), optional :: fade
      real(wp) :: hi_abs, hi_rel
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: prim
      real(wp), allocatable :: radii(:), points(:, :), v(:, :)
      real(wp), allocatable :: f1_rA(:, :), f2_r_rA(:, :, :), f3_rr_rA(:, :, :, :)
      real(wp), allocatable :: f4_rrr_rA(:, :, :, :, :)
      real(wp), allocatable :: f2_rArB(:, :, :, :), f3_r_rArB(:, :, :, :, :)
      real(wp), allocatable :: f4_rr_rArB(:, :, :, :, :, :)
      real(wp), allocatable :: h1(:, :), h2(:, :, :), h3(:, :, :, :)
      real(wp), allocatable :: vjp(:, :), vjp_rad(:)
      real(wp), allocatable :: p_f3(:, :, :), p_f4a(:, :, :), p_f4q(:, :, :, :, :)
      real(wp), allocatable :: p_h3(:, :, :), p_rad3(:, :)
      real(wp), allocatable :: rad1(:), rad2(:, :), rad3(:, :, :)
      real(wp) :: w0_adj, w1_adj(ndim), w2_adj(ndim, ndim), w2_sym(ndim, ndim), w2_packed(6)
      real(wp) :: p_t2(6), p_t3(10), acc_sym
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: ref0, ref1(ndim), ref2(ndim, ndim), ref3(ndim, ndim, ndim)
      real(wp) :: acc, point(ndim)
      integer  :: icase, ipt, iA, iB, s_ax, t_ax, j, k, l, nat, atomA, n_inband
      logical :: faded
      !> Band distances in sigma, fade on
      real(wp), allocatable :: d(:)
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc_param
      type(mctc_error), allocatable :: lsf_err

      faded = .false.
      if (present(fade)) faded = fade
      n_inband = 0
      call high_order_thr(kind, hi_abs, hi_rel)
      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         nat = mol%nat
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, high_order_points(kind))

         ! Deterministic asymmetric direction field
         if (allocated(v)) deallocate (v)
         allocate (v(ndim, nat))
         do iA = 1, nat
            v(1, iA) = 0.31_wp*real(iA, wp) - 0.7_wp
            v(2, iA) = -0.17_wp*real(iA, wp) + 0.4_wp
            v(3, iA) = 0.23_wp*real(mod(iA, 3) + 1, wp)
         end do

         if (allocated(f1_rA)) deallocate (f1_rA, f2_r_rA, f3_rr_rA, f4_rrr_rA, &
                                           f2_rArB, f3_r_rArB, f4_rr_rArB, h1, h2, h3, &
                                           vjp, vjp_rad, rad1, rad2, rad3, &
                                           p_f3, p_f4a, p_f4q, p_h3, p_rad3)
         allocate (f1_rA(ndim, nat), f2_r_rA(ndim, ndim, nat))
         allocate (f3_rr_rA(ndim, ndim, ndim, nat))
         allocate (f4_rrr_rA(ndim, ndim, ndim, ndim, nat))
         allocate (f2_rArB(ndim, nat, ndim, nat))
         allocate (f3_r_rArB(ndim, ndim, nat, ndim, nat))
         allocate (f4_rr_rArB(ndim, ndim, ndim, nat, ndim, nat))
         allocate (h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat))
         allocate (vjp(ndim, nat), vjp_rad(nat))
         allocate (rad1(nat), rad2(ndim, nat), rad3(ndim, ndim, nat))
         allocate (p_f3(6, ndim, nat), p_f4a(10, ndim, nat), p_f4q(6, ndim, nat, ndim, nat))
         allocate (p_h3(6, ndim, nat), p_rad3(6, nat))

         if (faded) then
            call init_lsf(prim, mol, radii, 4, kind, fade_enabled=.true., fade_start=5.0_wp, &
                          fade_width=2.0_wp)
         else
            call init_lsf(prim, mol, radii, 4, kind, blend_k=svdw_legacy_blend_k, &
                          blend_2b=svdw_legacy_blend_2b, blend_3b=svdw_legacy_blend_3b)
         end if

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            call prim%prepare(point, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if
            if (prim%unreached()) then
               if (faded) cycle
               call test_failed(error, "LSF prepare: point unreached")
               return
            end if
            if (faded) then
               d = band_distances(mol, radii, point)/radii(1:nat)
               if (any(d > 5.0_wp/abs(cfc_param%a2) .and. d < 7.0_wp/abs(cfc_param%a2))) then
                  n_inband = n_inband + 1
               end if
            end if

            call prim%f3_rr_rA(f1_rA, f2_r_rA, p_f3)
            call prim%f4_rrr_rA(p_f4a)
            call prim%f2_rArB(f2_rArB)
            call prim%f3_r_rArB(f3_r_rArB)
            call prim%f4_rr_rArB(p_f4q)
            call unpack_pair_batch(p_f3, f3_rr_rA, ndim*nat)
            call unpack_triple_batch(p_f4a, f4_rrr_rA, ndim*nat)
            call unpack_pair_batch(p_f4q, f4_rr_rArB, ndim*nat*ndim*nat)

            !* ------------------------------ tangent_* ------------------------------- *!
            ref0 = 0.0_wp
            ref1 = 0.0_wp
            ref2 = 0.0_wp
            ref3 = 0.0_wp
            do iA = 1, prim%active_count()
               atomA = prim%active_atom(iA)
               do s_ax = 1, ndim
                  ref0 = ref0 + v(s_ax, atomA)*f1_rA(s_ax, iA)
                  ref1 = ref1 + v(s_ax, atomA)*f2_r_rA(:, s_ax, iA)
                  ref2 = ref2 + v(s_ax, atomA)*f3_rr_rA(:, :, s_ax, iA)
                  ref3 = ref3 + v(s_ax, atomA)*f4_rrr_rA(:, :, :, s_ax, iA)
               end do
            end do

            call prim%tangent_f0(v, t0)
            call check(error, t0, ref0, thr_abs=ABS_THR, thr_rel=REL_THR)
            if (allocated(error)) return
            call prim%tangent_f1_r(v, t1)
            do j = 1, ndim
               call check(error, t1(j), ref1(j), thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
            call prim%tangent_f2_rr(v, p_t2)
            call unpack_sym_matrix(p_t2, t2, ndim)
            do k = 1, ndim
               do j = 1, ndim
                  call check(error, t2(j, k), ref2(j, k), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
            call prim%tangent_f3_rrr(v, p_t3)
            call unpack_sym_tensor3(p_t3, t3, ndim)
            do l = 1, ndim
               do k = 1, ndim
                  do j = 1, ndim
                     call check(error, t3(j, k, l), ref3(j, k, l), &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do

            !* -------------------------------- hvp_* --------------------------------- *!
            call prim%hvp_f1_rA(v, h1)
            call prim%hvp_f2_r_rA(v, h2)
            call prim%hvp_f3_rr_rA(v, p_h3)
            call unpack_pair_batch(p_h3, h3, ndim*nat)
            do iA = 1, prim%active_count()
               do s_ax = 1, ndim
                  acc = 0.0_wp
                  do iB = 1, prim%active_count()
                     do t_ax = 1, ndim
                        acc = acc + f2_rArB(s_ax, iA, t_ax, iB)*v(t_ax, prim%active_atom(iB))
                     end do
                  end do
                  call check(error, h1(s_ax, iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return

                  do j = 1, ndim
                     acc = 0.0_wp
                     do iB = 1, prim%active_count()
                        do t_ax = 1, ndim
                           acc = acc + f3_r_rArB(j, s_ax, iA, t_ax, iB) &
                                 *v(t_ax, prim%active_atom(iB))
                        end do
                     end do
                     call check(error, h2(j, s_ax, iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                     do k = 1, ndim
                        acc = 0.0_wp
                        do iB = 1, prim%active_count()
                           do t_ax = 1, ndim
                              acc = acc + f4_rr_rArB(j, k, s_ax, iA, t_ax, iB) &
                                    *v(t_ax, prim%active_atom(iB))
                           end do
                        end do
                        call check(error, h3(j, k, s_ax, iA), acc, &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do

            !* ------------------------------ vjp_f1_rA ------------------------------- *!
            ! Reverse `tangent_*`: contract jet indices, retain nuclear index
            ! Asymmetric `w2` exposes folded Hessian slots
            w0_adj = 0.37_wp
            w1_adj = [0.19_wp, -0.53_wp, 0.71_wp]
            do k = 1, ndim
               do j = 1, ndim
                  w2_adj(j, k) = 0.11_wp*real(j, wp) - 0.29_wp*real(k, wp) &
                                 + 0.07_wp*real(j*k, wp)
               end do
            end do

            ! The packed seed carries the symmetric part of `w2_adj`; the dense
            ! contraction with the full `w2_adj` (symmetric Hessian) must agree
            call check(error, maxval(abs(w2_adj - transpose(w2_adj))) > 0.1_wp, &
                       more="vjp seed is not asymmetric")
            if (allocated(error)) return
            call pack_sym_part(w2_adj, w2_packed, ndim)
            w2_sym = 0.5_wp*(w2_adj + transpose(w2_adj))

            call prim%vjp_f1_rA(w0_adj, w1_adj, w2_packed, vjp)
            do iA = 1, prim%active_count()
               do s_ax = 1, ndim
                  acc = w0_adj*f1_rA(s_ax, iA) &
                        + dot_product(w1_adj, f2_r_rA(:, s_ax, iA)) &
                        + sum(w2_adj*f3_rr_rA(:, :, s_ax, iA))
                  call check(error, vjp(s_ax, iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
                  acc_sym = w0_adj*f1_rA(s_ax, iA) &
                            + dot_product(w1_adj, f2_r_rA(:, s_ax, iA)) &
                            + sum(w2_sym*f3_rr_rA(:, :, s_ax, iA))
                  call check(error, vjp(s_ax, iA), acc_sym, thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do

            !* ------------------------------ vjp_f1_rad ------------------------------ *!
            ! Radius-ladder contraction; one scalar per active atom
            call prim%f3_rr_rad(rad1, rad2, p_rad3)
            call unpack_pair_batch(p_rad3, rad3, nat)
            call prim%vjp_f1_rad(w0_adj, w1_adj, w2_packed, vjp_rad)
            do iA = 1, prim%active_count()
               acc = w0_adj*rad1(iA) &
                     + dot_product(w1_adj, rad2(:, iA)) &
                     + sum(w2_adj*rad3(:, :, iA))
               call check(error, vjp_rad(iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
               acc_sym = w0_adj*rad1(iA) &
                         + dot_product(w1_adj, rad2(:, iA)) &
                         + sum(w2_sym*rad3(:, :, iA))
               call check(error, vjp_rad(iA), acc_sym, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do

            !* ------------------ nuclear Hessian exchange symmetry ------------------- *!
            ! Swap both nuclear slots in d2S/(dR_A dR_B)
            ! Compare distinct `qq` pair calls with A/B reversed
            do iB = 1, prim%active_count()
               do iA = 1, prim%active_count()
                  do t_ax = 1, ndim
                     do s_ax = 1, ndim
                        call check(error, f2_rArB(s_ax, iA, t_ax, iB), &
                                   f2_rArB(t_ax, iB, s_ax, iA), &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (faded) call check(error, n_inband > 0, more="no reached point with an atom in its band")
   end subroutine run_contracted

   !* ================================================================================= *!
   !*                          A point sitting on a nucleus                             *!
   !* ================================================================================= *!

   !> Require finite accessors at a nucleus
   !>
   !> `||r - R_A||` is nondifferentiable at nucleus; kernel divides by it
   !> `atom_tensors` keeps value and zeros derivative contributions
   !> Exercise guard through every accessor
   !>
   !> Check continuous `f0` against a nearby point
   subroutine test_svdw_on_nucleus(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      integer, parameter :: nat = 3
      !> Displacement used for the continuity probe (Bohr)
      real(wp), parameter :: nudge = 1.0e-7_wp
      integer :: atomic_numbers(nat)
      real(wp) :: centers(ndim, nat), radii(nat), point(ndim)
      real(wp) :: v(ndim, nat)
      real(wp) :: f0_on, f0_near, lsf0, lsf1_r(ndim), lsf2_rr(6)
      real(wp) :: f3_rrr(10), f4_rrrr(15)
      real(wp) :: f1_rA(ndim, nat), f2_r_rA(ndim, ndim, nat)
      real(wp) :: f3_rr_rA(6, ndim, nat)
      real(wp) :: f4_rrr_rA(10, ndim, nat)
      real(wp) :: f2_rArB(ndim, nat, ndim, nat)
      real(wp) :: f3_r_rArB(ndim, ndim, nat, ndim, nat)
      real(wp) :: f4_rr_rArB(6, ndim, nat, ndim, nat)
      real(wp) :: norm0, norm1_rA(ndim, nat)
      real(wp) :: t0, t1(ndim), t2(6), t3(10)
      real(wp) :: h1(ndim, nat), h2(ndim, ndim, nat), h3(6, ndim, nat)
      integer :: iat
      type(mctc_error), allocatable :: lsf_err

      atomic_numbers = 1
      centers = reshape([ &
                        -1.40_wp, 0.20_wp, -0.10_wp, &
                        1.25_wp, -0.30_wp, 0.35_wp, &
                        0.15_wp, 1.10_wp, -0.45_wp], [ndim, nat])
      radii = [0.80_wp, 1.00_wp, 0.75_wp]
      do iat = 1, nat
         v(:, iat) = [0.3_wp*iat, -0.2_wp*iat, 0.11_wp]
      end do

      call new(mol, atomic_numbers, centers)
      prim%screening_threshold = 0.0_wp
      call prim%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=svdw_legacy_blend_k, &
         blend_2b=svdw_legacy_blend_2b, blend_3b=svdw_legacy_blend_3b))
      call prim%update(mol, radii)
      call prim%set_max_deriv(4)

      ! Exactly on nucleus 2
      point = centers(:, 2)
      call prim%prepare(point, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      call check(error, prim%active_count(), nat, "on-nucleus point lost an atom")
      if (allocated(error)) return

      call prim%f0(f0_on)
      call check(error, finite0(f0_on), "f0 is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f012_r(lsf0, lsf1_r, lsf2_rr)
      call check(error, finite0(lsf0) .and. finite1(lsf1_r) .and. finiten(lsf2_rr), &
                 "f012_r is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f3_rrr(lsf3_rrr=f3_rrr)
      call check(error, finiten(f3_rrr), "f3_rrr is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f4_rrrr(f4_rrrr)
      call check(error, finiten(reshape(f4_rrrr, [size(f4_rrrr)])), &
                 "f4_rrrr is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f3_rr_rA(f1_rA, f2_r_rA, f3_rr_rA)
      call check(error, finiten(reshape(f1_rA, [size(f1_rA)])) &
                 .and. finiten(reshape(f2_r_rA, [size(f2_r_rA)])) &
                 .and. finiten(reshape(f3_rr_rA, [size(f3_rr_rA)])), &
                 "f3_rr_rA is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f4_rrr_rA(f4_rrr_rA)
      call check(error, finiten(reshape(f4_rrr_rA, [size(f4_rrr_rA)])), &
                 "f4_rrr_rA is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f2_rArB(f2_rArB)
      call check(error, finiten(reshape(f2_rArB, [size(f2_rArB)])), &
                 "f2_rArB is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f3_r_rArB(f3_r_rArB)
      call check(error, finiten(reshape(f3_r_rArB, [size(f3_r_rArB)])), &
                 "f3_r_rArB is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f4_rr_rArB(f4_rr_rArB)
      call check(error, finiten(reshape(f4_rr_rArB, [size(f4_rr_rArB)])), &
                 "f4_rr_rArB is not finite on a nucleus")
      if (allocated(error)) return

      call prim%normalized_f01_rA(norm0, deriv_rA=norm1_rA)
      call check(error, finite0(norm0) .and. finiten(reshape(norm1_rA, [size(norm1_rA)])), &
                 "normalized_f01_rA is not finite on a nucleus")
      if (allocated(error)) return

      call prim%tangent_f0(v, t0)
      call prim%tangent_f1_r(v, t1)
      call prim%tangent_f2_rr(v, t2)
      call prim%tangent_f3_rrr(v, t3)
      call check(error, finite0(t0) .and. finite1(t1) .and. finiten(t2) .and. finiten(t3), &
                 "tangent_* is not finite on a nucleus")
      if (allocated(error)) return

      call prim%hvp_f1_rA(v, h1)
      call prim%hvp_f2_r_rA(v, h2)
      call prim%hvp_f3_rr_rA(v, h3)
      call check(error, finiten(reshape(h1, [size(h1)])) &
                 .and. finiten(reshape(h2, [size(h2)])) &
                 .and. finiten(reshape(h3, [size(h3)])), &
                 "hvp_* is not finite on a nucleus")
      if (allocated(error)) return

      ! Value continuous across nucleus; slope discontinuous
      point = centers(:, 2) + [nudge, 0.0_wp, 0.0_wp]
      call prim%prepare(point, lsf_err)
      call prim%f0(f0_near)
      call check(error, f0_on, f0_near, thr_abs=1.0e-6_wp, thr_rel=1.0e-6_wp, &
                 more="f0 jumps at a nucleus")
   end subroutine test_svdw_on_nucleus

   !> Finiteness of a scalar
   !>
   !> @param[in] x Value to test
   !> @returns     `.true.` unless `x` is NaN or infinite
   pure logical function finite0(x) result(ok)
      !> Value to test
      real(wp), intent(in) :: x
      ok = (x == x) .and. abs(x) <= huge(1.0_wp)
   end function finite0

   !> Finiteness of a rank-1 array
   !>
   !> @param[in] x Array to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finiten(x) result(ok)
      !> Array to test
      real(wp), intent(in) :: x(:)
      integer :: i
      ok = .true.
      do i = 1, size(x)
         if (.not. finite0(x(i))) ok = .false.
      end do
   end function finiten

   !> Finiteness of a vector
   !>
   !> @param[in] x Vector to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finite1(x) result(ok)
      !> Vector to test
      real(wp), intent(in) :: x(:)
      ok = finiten(x)
   end function finite1

   !> Signed distance from a point to a sphere surface
   !>
   !> Reference for body-order reduction
   !> Local expression; LSF has no signed-distance module
   !>
   !> @param[in] point  Evaluation coordinates [3]
   !> @param[in] center Sphere center [3]
   !> @param[in] radius Sphere radius
   !> @returns   d      Signed distance to the sphere surface
   pure function ssd_ref(point, center, radius) result(d)
      !> Evaluation coordinates
      real(wp), intent(in) :: point(ndim)
      !> Sphere center
      real(wp), intent(in) :: center(ndim)
      !> Sphere radius
      real(wp), intent(in) :: radius
      !> Signed distance
      real(wp) :: d

      d = norm2(point - center) - radius
   end function ssd_ref

   !> Check SvdW body-order weight reductions
   !>
   !> `blend_2b=1`: mean of two atom SSDs
   !> `blend_3b=1`: mean of three atom SSDs
   subroutine test_svdw_body_order_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      integer, parameter :: n2 = 2, n3 = 3
      integer :: atomic_numbers(3)
      real(wp) :: centers2(ndim, n2), radii2(n2), point2(ndim), d2(n2)
      real(wp) :: centers3(ndim, n3), radii3(n3), point3(ndim), d3(n3)
      real(wp) :: lsf0, expected
      type(mctc_error), allocatable :: lsf_err

      atomic_numbers = 1

      centers2 = reshape([ &
                         -1.20_wp, 0.10_wp, 0.00_wp, &
                         1.10_wp, -0.20_wp, 0.30_wp], [ndim, n2])
      radii2 = [0.85_wp, 1.05_wp]
      point2 = [0.35_wp, 0.45_wp, -0.25_wp]

      call new(mol, atomic_numbers(:n2), centers2)
      prim%screening_threshold = 0.0_wp
      call prim%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.4_wp, blend_1b=0.0_wp, &
         blend_2b=1.0_wp, blend_3b=0.0_wp))
      call prim%update(mol, radii2)
      call prim%set_max_deriv(0)
      call prim%set_centers(centers2)
      call prim%prepare(point2, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      call prim%f0(lsf0)
      d2(1) = ssd_ref(point2, centers2(:, 1), radii2(1))
      d2(2) = ssd_ref(point2, centers2(:, 2), radii2(2))
      expected = 0.5_wp*sum(d2)
      call check(error, lsf0, expected, thr_abs=ABS_THR, thr_rel=REL_THR)
      if (allocated(error)) return

      centers3 = reshape([ &
                         -1.40_wp, 0.20_wp, -0.10_wp, &
                         1.25_wp, -0.30_wp, 0.35_wp, &
                         0.15_wp, 1.10_wp, -0.45_wp], [ndim, n3])
      radii3 = [0.80_wp, 1.00_wp, 0.75_wp]
      point3 = [0.10_wp, 0.25_wp, 0.40_wp]

      call new(mol, atomic_numbers(:n3), centers3)
      prim%screening_threshold = 0.0_wp
      call prim%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.7_wp, blend_1b=0.0_wp, &
         blend_2b=0.0_wp, blend_3b=1.0_wp))
      call prim%update(mol, radii3)
      call prim%set_max_deriv(0)
      call prim%set_centers(centers3)
      call prim%prepare(point3, lsf_err)
      call prim%f0(lsf0)
      d3(1) = ssd_ref(point3, centers3(:, 1), radii3(1))
      d3(2) = ssd_ref(point3, centers3(:, 2), radii3(2))
      d3(3) = ssd_ref(point3, centers3(:, 3), radii3(3))
      expected = sum(d3)/3.0_wp
      call check(error, lsf0, expected, thr_abs=ABS_THR, thr_rel=REL_THR)
   end subroutine test_svdw_body_order_scaling

   !* ================================================================================= *!
   !*                                 Atomic partition                                  *!
   !* ================================================================================= *!

   !> SvdW shares against the radius gradient
   !>
   !> Z is homogeneous of degree three in the screening factors and
   !> `du_A/dR_a = (k/3) u_A`, so the equal-split share `Z_A/Z` is `-dS/dR_a`;
   !> checked over the blend sweep, where the shares also sum to one
   subroutine test_svdw_atom_partition(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :), shares(:), f1_rad(:), f3_rr_rad(:, :)
      integer :: icase, ipt, iblend, igamma, nat

      call get_test_structures(mols)
      do icase = 1, size(mols)
         nat = mols(icase)%nat
         call get_test_radii(mols(icase), radii)
         call get_test_points(mols(icase), points, n_points_svdw)
         allocate (shares(nat), f1_rad(nat), f3_rr_rad(6, nat))
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               call init_lsf(lsf, mols(icase), radii, 2, "svdw", blend_k=svdw_blend_k_values(iblend), &
                             fade_enabled=fd_user_space_fade, &
                             blend_2b=svdw_2b_values(igamma), blend_3b=svdw_gamma_values(igamma))
               do ipt = 1, size(points, 2)
                  if (prepare_failed(lsf, points(:, ipt), error)) return
                  call check_identity_active(error, lsf, nat)
                  if (allocated(error)) return
                  shares = 0.0_wp
                  call lsf%atom_partition(1.0_wp, shares)
                  call lsf%f3_rr_rad(lsf1_rad=f1_rad, lsf3_rr_rad=f3_rr_rad)
                  call check(error, maxval(abs(shares + f1_rad)), 0.0_wp, thr=PARTITION_THR, &
                             more="SvdW share differs from -dS/dR_a")
                  if (allocated(error)) return
                  call check(error, sum(shares), 1.0_wp, thr=PARTITION_THR, more="SvdW shares do not sum to one")
                  if (allocated(error)) return
               end do
            end do
         end do
         deallocate (shares, f1_rad, f3_rr_rad)
      end do
   end subroutine test_svdw_atom_partition

   !> CFC shares against the pseudo-densities of sub-structures
   !>
   !> The atomic term of A is the pseudo-density of A alone and a pair term that of
   !> the pair less both atomic terms, all read back as `exp(-S)` from the level set
   !> itself; the share of A is its atomic term plus half of each of its pair terms
   subroutine test_cfc_atom_partition(error)
      type(error_type), allocatable, intent(out) :: error

      !> Leading atoms of each test structure forming the fixture
      integer, parameter :: nsub = 4
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: sub
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: shares(nsub), expect(nsub), pd_atom(nsub), pd_pair, pd
      integer :: icase, ipt, ia, ib

      call get_test_structures(mols)
      do icase = 1, size(mols)
         if (mols(icase)%nat < nsub) cycle
         call new(sub, mols(icase)%num(mols(icase)%id(1:nsub)), mols(icase)%xyz(:, 1:nsub))
         call get_test_radii(sub, radii)
         call get_test_points(sub, points, n_points_svdw)
         call init_lsf(lsf, sub, radii, 0, "cfc")
         do ipt = 1, size(points, 2)
            call cfc_pseudo_density(error, sub, radii, [1, 2, 3, 4], points(:, ipt), pd)
            if (allocated(error)) return
            do ia = 1, nsub
               call cfc_pseudo_density(error, sub, radii, [ia], points(:, ipt), pd_atom(ia))
               if (allocated(error)) return
            end do
            expect = pd_atom
            do ia = 1, nsub
               do ib = ia + 1, nsub
                  call cfc_pseudo_density(error, sub, radii, [ia, ib], points(:, ipt), pd_pair)
                  if (allocated(error)) return
                  pd_pair = pd_pair - pd_atom(ia) - pd_atom(ib)
                  expect(ia) = expect(ia) + 0.5_wp*pd_pair
                  expect(ib) = expect(ib) + 0.5_wp*pd_pair
               end do
            end do

            if (prepare_failed(lsf, points(:, ipt), error)) return
            call check_identity_active(error, lsf, nsub)
            if (allocated(error)) return
            shares = 0.0_wp
            call lsf%atom_partition(1.0_wp, shares)
            call check(error, maxval(abs(shares - expect/pd)), 0.0_wp, thr=PARTITION_THR, &
                       more="CFC share differs from its sub-structure terms")
            if (allocated(error)) return
            call check(error, sum(shares), 1.0_wp, thr=PARTITION_THR, more="CFC shares do not sum to one")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_cfc_atom_partition

   !> SSES shares against the radius gradient
   !>
   !> The shares are the optimized simplex weights. q_a depends on r_a only through
   !> its own offset, so the envelope theorem gives `dF/dR_a = w_a dq_a/dR_a`
   !> with `dq_a/dR_a = -2 R_a - 2 p R_a/sqrt(R_a**2 + eps**2)`
   subroutine test_sses_atom_partition(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      real(wp), allocatable :: radii(:), points(:, :), shares(:), slope(:)
      real(wp), allocatable :: f1_rad(:), f3_rr_rad(:, :)
      real(wp) :: p, eps
      integer :: icase, ipt, nat

      call get_test_structures(mols)
      do icase = 1, size(mols)
         nat = mols(icase)%nat
         if (nat < 2) cycle
         call get_test_radii(mols(icase), radii)
         call get_test_points(mols(icase), points, n_points_sses)
         allocate (shares(nat), slope(nat), f1_rad(nat), f3_rr_rad(6, nat))
         call lsf%new()
         call lsf%update(mols(icase), radii)
         call lsf%set_max_deriv(2)
         p = lsf%param%probe_radius
         eps = lsf%param%smoothing_length
         slope = lsf%param%scale*(-2.0_wp*radii - 2.0_wp*p*radii/sqrt(radii**2 + eps**2))
         do ipt = 1, size(points, 2)
            if (prepare_failed(lsf, points(:, ipt), error)) return
            call check_identity_active(error, lsf, nat)
            if (allocated(error)) return
            shares = 0.0_wp
            call lsf%atom_partition(1.0_wp, shares)
            call lsf%f3_rr_rad(lsf1_rad=f1_rad, lsf3_rr_rad=f3_rr_rad)
            call check(error, maxval(abs(shares - f1_rad/slope)), 0.0_wp, thr=PARTITION_THR, &
                       more="SSES share differs from (dF/dR_a)/(dq_a/dR_a)")
            if (allocated(error)) return
            call check(error, sum(shares), 1.0_wp, thr=PARTITION_THR, more="SSES shares do not sum to one")
            if (allocated(error)) return
            call check(error, minval(shares) >= 0.0_wp, more="negative SSES share")
            if (allocated(error)) return
         end do
         deallocate (shares, slope, f1_rad, f3_rr_rad)
      end do
   end subroutine test_sses_atom_partition

   !> Unscreened CFC pseudo-density `exp(-S)` of a subset of atoms
   !>
   !> @param[out] error Set when the preparation fails
   !> @param[in]  mol   Structure the subset is taken from
   !> @param[in]  radii Per-atom radii of `mol`
   !> @param[in]  atoms Atoms of the subset
   !> @param[in]  point Evaluation point
   !> @param[out] pd    Pseudo-density of the subset at `point`
   subroutine cfc_pseudo_density(error, mol, radii, atoms, point, pd)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: radii(:)
      integer, intent(in) :: atoms(:)
      real(wp), intent(in) :: point(ndim)
      real(wp), intent(out) :: pd

      type(structure_type) :: part
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: s

      pd = 0.0_wp
      call new(part, mol%num(mol%id(atoms)), mol%xyz(:, atoms))
      call init_lsf(lsf, part, radii(atoms), 0, "cfc")
      if (prepare_failed(lsf, point, error)) return
      call lsf%f0(s)
      pd = exp(-s)
   end subroutine cfc_pseudo_density

   !* ================================================================================= *!
   !*                        Surface-free exclusion certificate                         *!
   !* ================================================================================= *!

   !> Probe certified surface-free ball boundaries
   !>
   !> `exclusion_radius(S(x))` promises no zero inside `B(x, r)`
   !> Require unchanged sign across each claimed boundary
   !>
   !> - 14 directions at `0.999 r`, to allow tight `r = |S|` claims
   !> - optional `n_random` uniform interior samples per ball (LCG, seeded)
   !> - a probe no atom reaches (empty active list, the per-point path's zero)
   !>   counts as exterior, `S = +Inf`
   !> - optional `n_near_uncertified`: points with `0 < |S| <= 1` that got no
   !>   radius, so the caller can rule out a vacuous pass near the surface
   !> - optional `n_near_loose`: points with `0 < |S| <= 1` and
   !>   `r < (1 - 1e-3) |S|`, so the caller can rule out a conservative radius
   !>
   !> Caller checks `n_certified > 0` to avoid vacuous success
   !>
   !> @param[inout] lsf                LSF bound to the structure the points came from
   !> @param[in]    points             Sampling points (3, npoint)
   !> @param[out]   n_certified        Points that claimed a positive radius
   !> @param[out]   error              Test failure
   !> @param[in]    n_random           Interior samples per ball (optional, default 0)
   !> @param[out]   n_near_uncertified Near-surface points without a radius (optional)
   !> @param[out]   n_near_loose       Near-surface points below `(1 - 1e-3) |S|` (optional)
   subroutine probe_exclusion_balls(lsf, points, n_certified, error, n_random, &
                                    n_near_uncertified, n_near_loose)
      !> LSF bound to the structure the points came from
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Sampling points
      real(wp), intent(in) :: points(:, :)
      !> Points that claimed a positive radius
      integer, intent(out) :: n_certified
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Interior samples per ball
      integer, intent(in), optional :: n_random
      !> Near-surface points without a radius
      integer, intent(out), optional :: n_near_uncertified
      !> Near-surface points with a loose radius
      integer, intent(out), optional :: n_near_loose

      !> Tightness a near-surface radius must reach, relative to `|S|`
      real(wp), parameter :: near_tight = 1.0_wp - 1.0e-3_wp
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: center(ndim), lsf_center, r
      integer :: ipoint, nrand, nnear, nloose
      logical :: flipped
      integer :: state

      nrand = 0
      if (present(n_random)) nrand = n_random
      n_certified = 0
      nnear = 0
      nloose = 0
      state = 2024

      do ipoint = 1, size(points, 2)
         center = points(:, ipoint)

         call lsf%prepare(center, lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "LSF prepare failed: "//lsf_err%message)
            return
         end if
         if (lsf%unreached()) cycle
         call lsf%f0(lsf_center)

         r = lsf%exclusion_radius(lsf_center)
         if (abs(lsf_center) <= 1.0_wp .and. lsf_center /= 0.0_wp .and. &
             r < near_tight*abs(lsf_center)) nloose = nloose + 1
         if (r <= 0.0_wp) then
            if (abs(lsf_center) <= 1.0_wp .and. lsf_center /= 0.0_wp) nnear = nnear + 1
            cycle
         end if
         n_certified = n_certified + 1

         call ball_sign_change(lsf, center, lsf_center, r, nrand, state, flipped, error)
         if (allocated(error)) return
         if (flipped) then
            call test_failed(error, &
                             "exclusion radius over-claims: S changes sign inside "// &
                             "the ball it certifies as surface-free")
            return
         end if
      end do
      if (present(n_near_uncertified)) n_near_uncertified = nnear
      if (present(n_near_loose)) n_near_loose = nloose
   end subroutine probe_exclusion_balls

   !> Whether `S` changes sign inside `B(center, r)` at the probe set
   !>
   !> - 14 directions at fractions 0.25, 0.5, 0.75, 0.999 of `r`, then
   !>   `n_random` uniform interior samples (rejection from the cube)
   !> - a probe with an empty active list counts as exterior
   !>
   !> @param[inout] lsf      LSF instance
   !> @param[in]    center   Ball center
   !> @param[in]    s_center `S` at the center, nonzero
   !> @param[in]    r        Ball radius
   !> @param[in]    n_random Interior samples
   !> @param[inout] state    LCG state
   !> @param[out]   flipped  A probe has the other sign (or zero)
   !> @param[out]   error    Test failure from `prepare`
   subroutine ball_sign_change(lsf, center, s_center, r, n_random, state, flipped, error)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Ball center
      real(wp), intent(in) :: center(ndim)
      !> `S` at the center
      real(wp), intent(in) :: s_center
      !> Ball radius
      real(wp), intent(in) :: r
      !> Interior samples
      integer, intent(in) :: n_random
      !> LCG state
      integer, intent(inout) :: state
      !> A probe has the other sign
      logical, intent(out) :: flipped
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Fractions of the radius probed along each direction
      real(wp), parameter :: fracs(4) = [0.25_wp, 0.5_wp, 0.75_wp, 0.999_wp]
      ! Deterministic probe directions for reproducible failures
      integer, parameter :: ndir = 14
      real(wp), parameter :: dirs(ndim, ndir) = reshape([ &
                             1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, &
                             0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
                             0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, &
                             1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, &
                             1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, &
                             1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp, &
                             1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp, 1.0_wp], &
                             [ndim, ndir])
      real(wp) :: probe(ndim), u(ndim)
      integer :: idir, ifrac, isample

      flipped = .false.
      do idir = 1, ndir
         do ifrac = 1, size(fracs)
            probe = center + fracs(ifrac)*r*dirs(:, idir)/norm2(dirs(:, idir))
            if (probe_flips(lsf, probe, s_center, error)) flipped = .true.
            if (allocated(error) .or. flipped) return
         end do
      end do
      isample = 0
      do while (isample < n_random)
         call det_rand(u, state, -1.0_wp, 1.0_wp)
         if (sum(u*u) > 1.0_wp) cycle
         isample = isample + 1
         probe = center + 0.999_wp*r*u
         if (probe_flips(lsf, probe, s_center, error)) flipped = .true.
         if (allocated(error) .or. flipped) return
      end do
   end subroutine ball_sign_change

   !> Whether `S` at `probe` lacks the sign of `s_center`; unreached is exterior
   !>
   !> @param[inout] lsf      LSF instance
   !> @param[in]    probe    Probe point
   !> @param[in]    s_center Reference value
   !> @param[out]   error    Test failure from `prepare`
   !> @returns      flips    Sign differs or `S = 0`
   function probe_flips(lsf, probe, s_center, error) result(flips)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Probe point
      real(wp), intent(in) :: probe(ndim)
      !> Reference value
      real(wp), intent(in) :: s_center
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Sign differs or zero
      logical :: flips

      type(mctc_error), allocatable :: lsf_err
      real(wp) :: s

      flips = .false.
      call lsf%prepare(probe, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      if (lsf%unreached()) then
         flips = s_center < 0.0_wp
         return
      end if
      call lsf%f0(s)
      flips = s*s_center <= 0.0_wp
   end function probe_flips

   !> Check SvdW exclusion radius across blend sharpness, then with the fade on
   !>
   !> - `get_test_points` over three `blend_k`, fade off: tight `r = |S|`
   !>   (1-Lipschitz); default fade on: the same loop at `|tau(S)|`
   !> - fade on at the defaults (21/3) and at a biting band (4/2), on MB16-43 01
   !>   and the 1000-atom lattice: points inside, near the surface, inside a band
   !>   and far (up to `S` past the cap), 14 directions x 4 fractions plus 64
   !>   random interior samples per ball; every near-surface point (`|S| <= 1`)
   !>   must be certified, at the defaults with `r >= (1 - 1e-3) |S|` (the
   !>   certificate equals `|S|` to 1e-4 there, so a 2x-conservative radius fails)
   !> - negative control: one atom (R = 3), band 18/7,
   !>   point at `t = 0.97` of the band on the x axis: `S = 20.5`, true distance 13.52; the ball `|S|` reaches
   !>   inside the sphere (sign change found), the certificate's ball does not
   subroutine test_svdw_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol, mols(2)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: x(ndim), s, r, dist, band(2, 2)
      integer :: iblend, ifade, imol, ifd, n_certified, n_near, n_loose
      logical :: flipped
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call get_test_points(mol, points)

      do ifade = 1, 2
         do iblend = 1, n_svdw_blends
            call init_lsf(lsf, mol, radii, 0, "svdw", blend_k=svdw_blend_k_values(iblend), &
                          fade_enabled=ifade == 2)

            call probe_exclusion_balls(lsf, points, n_certified, error)
            if (allocated(error)) return

            if (n_certified == 0) then
               call test_failed(error, "SvdW certified no ball at all: the probe "// &
                                "loop never ran")
               return
            end if
         end do
      end do

      ! Fade on: default and biting band, molecule and lattice
      band = reshape([21.0_wp, 3.0_wp, 4.0_wp, 2.0_wp], [2, 2])
      mols(1) = mol
      call lattice_structure(mols(2), 10, 3.0_wp, 17)
      do imol = 1, size(mols)
         call get_test_radii(mols(imol), radii)
         do ifd = 1, 2
            call init_lsf(lsf, mols(imol), radii, 0, "svdw", fade_enabled=.true., &
                          fade_start=band(1, ifd), fade_width=band(2, ifd))
            call certificate_points(lsf, mols(imol), radii, band(1, ifd), band(2, ifd), &
                                    points, error)
            if (allocated(error)) return
            call probe_exclusion_balls(lsf, points, n_certified, error, n_random=64, &
                                       n_near_uncertified=n_near, n_near_loose=n_loose)
            if (allocated(error)) return
            call check(error, n_certified > size(points, 2)/2, &
                       more="fade on: most category points must be certified")
            if (allocated(error)) return
            call check(error, n_near, 0, more="fade on: a near-surface point got no radius")
            if (allocated(error)) return
            if (ifd == 1) then
               call check(error, n_loose, 0, &
                          more="fade on, defaults: a near-surface radius below (1 - 1e-3) |S|")
               if (allocated(error)) return
            end if
         end do
      end do

      ! Negative control: raw |S| over-claims inside the band of a lone atom
      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]))
      radii = [3.0_wp]
      call init_lsf(lsf, mol, radii, 0, "svdw", fade_enabled=.true., fade_start=18.0_wp, &
                    fade_width=7.0_wp)
      x = [radii(1) + 3.0_wp*(18.0_wp + 0.97_wp*7.0_wp)/5.5_wp, 0.0_wp, 0.0_wp]
      dist = norm2(x) - radii(1)
      call lsf%prepare(x, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      call lsf%f0(s)
      call check(error, s/dist > 1.5_wp, more="control: |S| must exceed the true distance")
      if (allocated(error)) return
      state = 67
      call ball_sign_change(lsf, x, s, s, 64, state, flipped, error)
      if (allocated(error)) return
      call check(error, flipped, more="control: the |S| ball must contain surface")
      if (allocated(error)) return
      r = lsf%exclusion_radius(s)
      call check(error, r > 0.0_wp .and. r < dist, &
                 more="control: certificate must be positive and below the true distance")
      if (allocated(error)) return
      call ball_sign_change(lsf, x, s, r, 2000, state, flipped, error)
      if (allocated(error)) return
      call check(error, .not. flipped, more="control: the certificate ball holds surface")
   end subroutine test_svdw_exclusion_radius

   !> Category points for the faded certificate: inside, near, band, far
   !>
   !> - candidates `R_A + R_a + s` along 14 directions per atom (every 7th atom
   !>   on the lattice), `s` from `-0.6 R_a` to past the band end
   !> - kept by `S` at the candidate, up to 8 per category: inside (`S < 0`), near
   !>   (`0 < S <= 1`), mid (`1 < S <= X_on`), band (`S > X_on`, reached)
   !> - every category must be present
   !>
   !> @param[inout] lsf    Faded SvdW instance bound to `mol`
   !> @param[in]    mol    Structure
   !> @param[in]    radii  Per-atom radii
   !> @param[in]    a      Fade onset in e-folds
   !> @param[in]    b      Fade width in e-folds
   !> @param[out]   points Selected points [ndim, n]
   !> @param[out]   error  Test failure, a category is empty
   subroutine certificate_points(lsf, mol, radii, a, b, points, error)
      !> Faded SvdW instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Fade onset in e-folds
      real(wp), intent(in) :: a
      !> Fade width in e-folds
      real(wp), intent(in) :: b
      !> Selected points
      real(wp), allocatable, intent(out) :: points(:, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: per_cat = 8, nshell = 12
      real(wp) :: dirs(ndim, 14), p(ndim), s, x_on, band_end, shells(nshell)
      real(wp) :: buf(ndim, 4*per_cat)
      integer :: count_cat(4), icat, ia, j, k, n, stride
      type(mctc_error), allocatable :: lsf_err

      x_on = 3.0_wp*a/5.5_wp
      band_end = 3.0_wp*(a + b)/5.5_wp
      shells = [-0.6_wp, -0.2_wp, 0.05_wp, 0.3_wp, 0.8_wp, 0.5_wp*x_on, x_on, &
                0.5_wp*(x_on + band_end), band_end - 0.3_wp, band_end + 0.5_wp, &
                band_end + 2.0_wp, band_end + 5.0_wp]
      dirs = 0.0_wp
      do j = 1, ndim
         dirs(j, 2*j - 1) = 1.0_wp
         dirs(j, 2*j) = -1.0_wp
      end do
      k = 6
      do ia = -1, 1, 2
         do j = -1, 1, 2
            do n = -1, 1, 2
               k = k + 1
               dirs(:, k) = real([ia, j, n], wp)/sqrt(3.0_wp)
            end do
         end do
      end do

      stride = max(1, mol%nat/100)
      count_cat = 0
      n = 0
      do ia = 1, mol%nat, stride
         do k = 1, nshell
            do j = 1, size(dirs, 2)
               p = mol%xyz(:, ia) + (radii(ia) + shells(k)*merge(radii(ia), 1.0_wp, &
                                                                 shells(k) < 0.0_wp))*dirs(:, j)
               call lsf%prepare(p, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "prepare failed: "//lsf_err%message)
                  return
               end if
               if (lsf%unreached()) cycle
               call lsf%f0(s)
               if (s < 0.0_wp) then
                  icat = 1
               else if (s <= 1.0_wp) then
                  icat = 2
               else if (s <= x_on) then
                  icat = 3
               else
                  icat = 4
               end if
               if (count_cat(icat) >= per_cat) cycle
               count_cat(icat) = count_cat(icat) + 1
               n = n + 1
               buf(:, n) = p
            end do
         end do
      end do
      if (any(count_cat == 0)) then
         call test_failed(error, "certificate points: a category (inside, near, mid, band) "// &
                          "is empty")
         return
      end if
      points = buf(:, 1:n)
   end subroutine certificate_points

   !> SvdW exclusion radius gates
   !>
   !> - fade off: non-negative blends certify `r = |S|` exactly, a negative
   !>   coefficient nothing
   !> - fade on: `S = -1` certifies `1/(1 + K0)` (below 1, within 1e-6);
   !>   nothing at `s3 >= 24 s1`, at `blend_1b = 0`, at a negative coefficient,
   !>   for a NaN or `-Inf` value (both modes), and once `param` changed after
   !>   `update` (stale constants)
   !> - the far certificate needs `s >= 0` and `s1 > 0`: zero at `blend_1b = 0`
   !>   ([[test_svdw_far_blend_1b_zero]])
   subroutine test_svdw_exclusion_radius_gate(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)
      real(wp) :: r

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)

      call init_lsf(lsf, mol, radii, 0, "svdw", blend_2b=1.0_wp, fade_enabled=.false.)
      call check(error, lsf%exclusion_radius(-1.0_wp), 1.0_wp, thr=0.0_wp, &
                 message="fade off: non-negative blends should certify r = |S|")
      if (allocated(error)) return
      call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_negative_inf)), 0.0_wp, &
                 thr=0.0_wp, message="fade off: -Inf must certify nothing")
      if (allocated(error)) return
      call init_lsf(lsf, mol, radii, 0, "svdw", blend_2b=-1.0_wp, fade_enabled=.false.)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="fade off: a negative blend coefficient must certify nothing")
      if (allocated(error)) return

      call init_lsf(lsf, mol, radii, 0, "svdw", blend_2b=1.0_wp)
      r = lsf%exclusion_radius(-1.0_wp)
      call check(error, r < 1.0_wp .and. r > 1.0_wp - 1.0e-6_wp, &
                 more="fade on: S = -1 must certify 1/(1 + K0)")
      if (allocated(error)) return
      call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_quiet_nan)), 0.0_wp, &
                 thr=0.0_wp, message="fade on: NaN must certify nothing")
      if (allocated(error)) return
      call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_negative_inf)), 0.0_wp, &
                 thr=0.0_wp, message="fade on: -Inf must certify nothing")
      if (allocated(error)) return
      select type (lsf)
      type is (moist_cavity_drop_lsf_svdw_type)
         lsf%param%blend_k = 6.0_wp
         call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                    message="fade on: stale constants must certify nothing")
         if (allocated(error)) return
         call check(error, lsf%far_exclusion_radius(), 0.0_wp, thr=0.0_wp, &
                    message="fade on: stale constants must certify nothing (far)")
         if (allocated(error)) return
      end select

      call init_lsf(lsf, mol, radii, 0, "svdw", blend_2b=-1.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="fade on: a negative blend coefficient must certify nothing")
      if (allocated(error)) return
      call check(error, lsf%far_exclusion_radius(), 0.0_wp, thr=0.0_wp, &
                 message="fade on: a negative blend coefficient must certify nothing (far)")
      if (allocated(error)) return

      ! s3 >= 24 s1: no level-only bound
      call init_lsf_1b(lsf, mol, radii, 0.125_wp, 3.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="fade on: s3 = 24 s1 must certify nothing")
      if (allocated(error)) return
      call check(error, lsf%exclusion_radius(2.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="fade on: s3 = 24 s1 must certify nothing (S > 0)")
      if (allocated(error)) return
      call init_lsf_1b(lsf, mol, radii, 0.13_wp, 3.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp) > 0.0_wp, &
                 more="fade on: s3 just below 24 s1 must certify")
      if (allocated(error)) return

      ! blend_1b = 0: neither certificate
      call init_lsf_1b(lsf, mol, radii, 0.0_wp, 3.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="fade on: blend_1b = 0 must certify nothing")
      if (allocated(error)) return
      call check(error, lsf%far_exclusion_radius(), 0.0_wp, thr=0.0_wp, &
                 message="fade on: blend_1b = 0 must certify nothing (far)")
   end subroutine test_svdw_exclusion_radius_gate

   !> Faded SvdW with explicit one- and three-body weights (two-body off)
   !>
   !> @param[out] lsf    SvdW instance after `update`
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  s1     One-body weight
   !> @param[in]  s3     Three-body weight
   subroutine init_lsf_1b(lsf, mol, radii, s1, s3)
      !> SvdW instance
      class(moist_cavity_drop_lsf_type), allocatable, intent(out) :: lsf
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> One-body weight
      real(wp), intent(in) :: s1
      !> Three-body weight
      real(wp), intent(in) :: s3

      type(moist_cavity_drop_lsf_svdw_type) :: tmp
      type(moist_cavity_drop_lsf_svdw_param_type) :: param

      call param%new(blend_1b=s1, blend_3b=s3, fade_enabled=.true.)
      call tmp%new(param)
      call tmp%update(mol, radii)
      allocate (lsf, source=tmp)
   end subroutine init_lsf_1b

   !> Certificate constants against independently derived tables
   !>
   !> - k = 5.5, band 18/7 (the tabulated band), S_4, s = (1, 0, 3) and (1, 1, 1), n = 16 and 1000:
   !>   `d*`, far radius, `S_cap`, `r_cap` and `r(S)` at S = 1, 5, 9 within 1e-3 of
   !>   the independent Python reference `validate.py` (grid maxima)
   !> - geometry irrelevant beyond `n`: a cubic grid of carbons
   subroutine test_svdw_certificate_constants(error)
      type(error_type), allocatable, intent(out) :: error

      !> Reference rows: s set, n, d*, far, S_cap, r_cap, r(1), r(5), r(9)
      real(wp), parameter :: ref(7, 4) = reshape([ &
                             1.5123_wp, 12.1240_wp, 10.4425_wp, 9.8452_wp, 1.0000_wp, 5.0000_wp, 8.9451_wp, &
                             3.7679_wp, 9.8685_wp, 9.6905_wp, 7.9359_wp, 1.0000_wp, 4.9976_wp, 7.8921_wp, &
                             1.5123_wp, 12.1240_wp, 10.3975_wp, 9.8953_wp, 1.0000_wp, 5.0000_wp, 8.9671_wp, &
                             3.7679_wp, 9.8685_wp, 9.3757_wp, 8.1599_wp, 1.0000_wp, 4.9988_wp, 8.1343_wp], [7, 4])
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: xyz(ndim, 1000), got(7)
      integer :: num(1000), i, irow, n

      do i = 1, 1000
         xyz(:, i) = 5.0_wp*real([mod(i - 1, 10), mod((i - 1)/10, 10), (i - 1)/100], wp)
      end do
      num = 6
      do irow = 1, 4
         n = merge(16, 1000, mod(irow, 2) == 1)
         if (irow <= 2) then
            call param%new(blend_1b=1.0_wp, blend_2b=0.0_wp, blend_3b=3.0_wp, fade_enabled=.true., &
                           fade_start=18.0_wp, fade_width=7.0_wp)
         else
            call param%new(blend_1b=1.0_wp, blend_2b=1.0_wp, blend_3b=1.0_wp, fade_enabled=.true., &
                           fade_start=18.0_wp, fade_width=7.0_wp)
         end if
         call new(mol, num(1:n), xyz(:, 1:n))
         call lsf%new(param)
         call lsf%update(mol, [(3.0_wp, i=1, n)])
         got = [lsf%cert%dstar, lsf%far_exclusion_radius(), lsf%cert%s_cap, lsf%cert%r_cap, &
                lsf%exclusion_radius(1.0_wp), lsf%exclusion_radius(5.0_wp), &
                lsf%exclusion_radius(9.0_wp)]
         call check(error, maxval(abs(got - ref(:, irow))), 0.0_wp, thr=1.0e-3_wp, &
                    more="certificate constants differ from the derivation tables")
         if (allocated(error)) return
         call check(error, lsf%exclusion_radius(huge(1.0_wp)), lsf%cert%r_cap, thr=0.0_wp, &
                    more="S = huge must give r_cap")
         if (allocated(error)) return
      end do
   end subroutine test_svdw_certificate_constants

   !> Far certificate: no surface within `far_exclusion_radius()` of an unreached point
   !>
   !> - MB16-43 01 and the 1000-atom lattice, default fade: random points in the
   !>   box padded by band end + 8 Bohr, kept when `evaluate_batch` reports
   !>   the column unreached (no atom within the reach)
   !> - `far_exclusion_radius() > 0` at the defaults (11.6 Bohr at n = 16, 9.3 at
   !>   n = 1000)
   !> - each ball probed by [[ball_sign_change]]: 14 directions x 4 fractions and
   !>   64 random interior samples, a probe no atom reaches counting as exterior
   !> - the claim covers unreached points only: a point just inside the reach
   !>   (reached) may sit closer to the surface than the far radius, and is not
   !>   covered (documented, not asserted)
   subroutine test_svdw_far_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncand = 400, nkeep = 24
      type(structure_type) :: mols(2)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)
      real(wp) :: cand(ndim, ncand), b0(ncand), u(ndim), lo(ndim), hi(ndim), pad, r
      integer :: st(ncand), imol, i, nfar
      logical :: flipped
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mols(1), "MB16-43", "01")
      call lattice_structure(mols(2), 10, 3.0_wp, 23)
      state = 3
      do imol = 1, size(mols)
         call get_test_radii(mols(imol), radii)
         call init_lsf(lsf, mols(imol), radii, 0, "svdw", fade_enabled=.true.)
         r = lsf%far_exclusion_radius()
         call check(error, r > 0.0_wp, more="far certificate must be positive at the defaults")
         if (allocated(error)) return
         pad = lsf%screening_offset(radii(1)) + maxval(radii) + 8.0_wp
         lo = lsf%center_lower - pad
         hi = lsf%center_upper + pad
         do i = 1, ncand
            call det_rand(u, state)
            cand(:, i) = lo + (hi - lo)*u
         end do
         call lsf%evaluate_batch(cand, 0, b0, error=lsf_err)
         call batch_status(lsf, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         nfar = 0
         do i = 1, ncand
            if (st(i) /= 1) cycle
            nfar = nfar + 1
            if (nfar > nkeep) exit
            call ball_sign_change(lsf, cand(:, i), 1.0_wp, r, 64, state, flipped, error)
            if (allocated(error)) return
            if (flipped) then
               call test_failed(error, "far exclusion radius over-claims: surface inside "// &
                                "the ball around an unreached point")
               return
            end if
         end do
         call check(error, nfar > 0, more="no unreached point sampled")
         if (allocated(error)) return
      end do
   end subroutine test_svdw_far_exclusion_radius

   !> Far certificate off at `blend_1b = 0`, where an unreached point can have atoms in reach
   !>
   !> - A (R = 10) at the origin, B and C (R = 3) at (11, +-1, 0), weights
   !>   (0, 0, 3), default band (end 13.09 Bohr): at x = (-5.5, 0, 0) only A is in
   !>   reach, `Z = s3 e3 = 0`, `S = +Inf`, so unreached
   !> - y = (6, 0, 0), 11.5 Bohr from x, reaches all three atoms and is inside
   !>   (`S < 0`), so surface lies between x and y; the formula without the
   !>   `s1 > 0` guard gave 12.6 Bohr and claimed y
   !> - `far_exclusion_radius()` must be 0
   subroutine test_svdw_far_blend_1b_zero(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: xyz(ndim, 3), radii(3), pts(ndim, 2), b0(2), band_end, d_unguarded
      integer :: st(2)
      type(mctc_error), allocatable :: lsf_err

      xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 11.0_wp, 1.0_wp, 0.0_wp, &
                     11.0_wp, -1.0_wp, 0.0_wp], [ndim, 3])
      radii = [10.0_wp, 3.0_wp, 3.0_wp]
      call new(mol, [6, 6, 6], xyz)
      call param%new(blend_1b=0.0_wp, blend_2b=0.0_wp, blend_3b=3.0_wp, fade_enabled=.true.)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)

      pts(:, 1) = [-5.5_wp, 0.0_wp, 0.0_wp]
      pts(:, 2) = [6.0_wp, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(pts, 0, b0, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, st(1), 1, more="x must be unreached (Z = 0 with A in reach)")
      if (allocated(error)) return
      call check(error, st(2) == 0 .and. b0(2) < 0.0_wp, &
                 more="y must be reached and inside")
      if (allocated(error)) return

      ! Control: the unguarded radius (Zmax = s3/6, n = 3) would have covered y
      band_end = 3.0_wp*(param%fade_start + param%fade_width)/param%blend_k
      d_unguarded = band_end - (3.0_wp*log(3.0_wp) + log(0.5_wp))/param%blend_k
      call check(error, norm2(pts(:, 2) - pts(:, 1)) < d_unguarded, &
                 more="control: y must lie inside the unguarded far radius")
      if (allocated(error)) return

      call check(error, lsf%far_exclusion_radius(), 0.0_wp, thr=0.0_wp, &
                 message="blend_1b = 0: the far certificate must be off")
   end subroutine test_svdw_far_blend_1b_zero

   !> Check CFC exclusion radius, fade off and on
   !>
   !> - `get_test_structures` (10), fade off (`r = |S|/L`, cached Lipschitz
   !>   bound) and fade on at the defaults (`|tau(S)|`): no sign change in any
   !>   certified ball
   !> - fade on at the defaults (27/3) and at a biting band (5/2), on MB16-43 01
   !>   and a 512-atom lattice: category points ([[cfc_certificate_points]]:
   !>   inside, near, mid, band, far), 14 directions x 4 fractions plus 64
   !>   random interior samples per ball; every near-surface point (`|S| <= 1`)
   !>   certified, and at the defaults every point with `|S| <= 10` at
   !>   `>= (1 - 1e-3) |S|/L` (`|S|/L` from a fade-off twin, same `L`; the
   !>   certificate equals it to ~1e-4 there, so a 2x-conservative radius fails)
   !> - negative control: one carbon (R = 2), default
   !>   band, point on the x axis at `t = 0.97`: raw `|S|/L` exceeds the true
   !>   distance `R sigma` (sign change found in its ball), the certificate is
   !>   positive, below it, and its ball holds no surface
   subroutine test_cfc_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol, cat_mols(2)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf, lsf_off
      type(moist_cavity_drop_lsf_cfc_param_type) :: defaults
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: band(2, 2), x(ndim), s, r, r_raw, dist, sigma
      integer :: icase, ifade, imol, ifd, ip, n_certified, n_total, n_near, n_loose
      logical :: flipped
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols, 10)
      do ifade = 1, 2
         n_total = 0
         do icase = 1, size(mols)
            call get_test_radii(mols(icase), radii)
            call get_test_points(mols(icase), points)
            call init_lsf(lsf, mols(icase), radii, 0, "cfc", fade_enabled=ifade == 2)
            call probe_exclusion_balls(lsf, points, n_certified, error)
            if (allocated(error)) return
            n_total = n_total + n_certified
         end do
         call check(error, n_total > 0, &
                    more="CFC certified no ball at all: the probe loop never ran")
         if (allocated(error)) return
      end do

      ! Fade on: default and biting band, molecule and lattice
      band = reshape([defaults%fade_start, defaults%fade_width, 5.0_wp, 2.0_wp], [2, 2])
      call get_structure(cat_mols(1), "MB16-43", "01")
      call lattice_structure(cat_mols(2), 8, 3.0_wp, 17)
      do imol = 1, size(cat_mols)
         call get_test_radii(cat_mols(imol), radii)
         do ifd = 1, 2
            call init_lsf(lsf, cat_mols(imol), radii, 0, "cfc", fade_enabled=.true., &
                          fade_start=band(1, ifd), fade_width=band(2, ifd))
            call cfc_certificate_points(lsf, cat_mols(imol), radii, &
                                        band(1, ifd)/abs(defaults%a2), &
                                        band(2, ifd)/abs(defaults%a2), points, error)
            if (allocated(error)) return
            call probe_exclusion_balls(lsf, points, n_certified, error, n_random=64, &
                                       n_near_uncertified=n_near)
            if (allocated(error)) return
            call check(error, n_certified > size(points, 2)/2, &
                       more="fade on: most category points must be certified")
            if (allocated(error)) return
            call check(error, n_near, 0, more="fade on: a near-surface point got no radius")
            if (allocated(error)) return
            if (ifd /= 1) cycle

            ! Defaults: the faded certificate stays at today's |S|/L up to S = 10
            call init_lsf(lsf_off, cat_mols(imol), radii, 0, "cfc", fade_enabled=.false.)
            n_loose = 0
            do ip = 1, size(points, 2)
               call lsf%prepare(points(:, ip), lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               if (lsf%unreached()) cycle
               call lsf%f0(s)
               if (abs(s) > 10.0_wp) cycle
               if (lsf%exclusion_radius(s) < (1.0_wp - 1.0e-3_wp)*lsf_off%exclusion_radius(s)) &
                  n_loose = n_loose + 1
            end do
            call check(error, n_loose, 0, &
                       more="fade on, defaults: a radius below (1 - 1e-3) |S|/L at |S| <= 10")
            if (allocated(error)) return
         end do
      end do

      ! Negative control: raw |S|/L over-claims inside the band of a lone atom
      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]))
      radii = [2.0_wp]
      call init_lsf(lsf, mol, radii, 0, "cfc", fade_enabled=.true.)
      call init_lsf(lsf_off, mol, radii, 0, "cfc", fade_enabled=.false.)
      sigma = (defaults%fade_start + 0.97_wp*defaults%fade_width)/abs(defaults%a2)
      x = [radii(1)*(1.0_wp + sigma), 0.0_wp, 0.0_wp]
      dist = radii(1)*sigma
      call lsf%prepare(x, lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "LSF prepare failed: "//lsf_err%message)
         return
      end if
      call lsf%f0(s)
      r_raw = lsf_off%exclusion_radius(s)
      call check(error, r_raw/dist > 1.2_wp, &
                 more="control: |S|/L must exceed the true distance")
      if (allocated(error)) return
      state = 67
      call ball_sign_change(lsf, x, s, r_raw, 64, state, flipped, error)
      if (allocated(error)) return
      call check(error, flipped, more="control: the |S|/L ball must contain surface")
      if (allocated(error)) return
      r = lsf%exclusion_radius(s)
      call check(error, r > 0.0_wp .and. r < dist, &
                 more="control: certificate must be positive and below the true distance")
      if (allocated(error)) return
      call ball_sign_change(lsf, x, s, r, 2000, state, flipped, error)
      if (allocated(error)) return
      call check(error, .not. flipped, more="control: the certificate ball holds surface")
   end subroutine test_cfc_exclusion_radius

   !> Category points for the faded CFC certificate: inside, near, mid, band, far
   !>
   !> - candidates `R_A + R_a (1 + s)` along 14 directions per atom (about 30
   !>   atoms), `s` in sigma from `-0.6` to past the band end, values from one
   !>   order-0 `evaluate_batch`
   !> - kept by `S` and the band, up to 8 per category: inside (`S < 0`), near
   !>   (`0 < S <= 1`), mid (`S > 1`, no atom in its band), band (some atom in
   !>   its band, `S <= s_cap`), far (reached, `S > s_cap`: the radius is `r_cap`)
   !> - every category must be present
   !>
   !> @param[inout] lsf    Faded CFC instance bound to `mol`
   !> @param[in]    mol    Structure
   !> @param[in]    radii  Per-atom radii
   !> @param[in]    x_on   Fade onset in sigma
   !> @param[in]    delta  Fade band width in sigma
   !> @param[out]   points Selected points [ndim, n]
   !> @param[out]   error  Test failure, a category is empty
   subroutine cfc_certificate_points(lsf, mol, radii, x_on, delta, points, error)
      !> Faded CFC instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Fade onset in sigma
      real(wp), intent(in) :: x_on
      !> Fade band width in sigma
      real(wp), intent(in) :: delta
      !> Selected points
      real(wp), allocatable, intent(out) :: points(:, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: per_cat = 8, nshell = 11, ncat = 5
      real(wp) :: dirs(ndim, 14), shells(nshell), sig(mol%nat), s_cap
      real(wp), allocatable :: cand(:, :), b0(:)
      real(wp) :: buf(ndim, ncat*per_cat)
      integer, allocatable :: st(:)
      integer :: count_cat(ncat), icat, ia, j, k, n, nc, stride
      type(mctc_error), allocatable :: lsf_err

      select type (lsf)
      type is (moist_cavity_drop_lsf_cfc_type)
         s_cap = lsf%cert%s_cap
      class default
         error stop "cfc_certificate_points: not a CFC instance"
      end select
      shells = [-0.6_wp, -0.2_wp, 0.02_wp, 0.1_wp, 0.3_wp, 0.5_wp*x_on, &
                x_on + 0.2_wp*delta, x_on + 0.5_wp*delta, x_on + 0.9_wp*delta, &
                x_on + delta + 0.1_wp, x_on + delta + 0.5_wp]
      dirs = 0.0_wp
      do j = 1, ndim
         dirs(j, 2*j - 1) = 1.0_wp
         dirs(j, 2*j) = -1.0_wp
      end do
      k = 6
      do ia = -1, 1, 2
         do j = -1, 1, 2
            do n = -1, 1, 2
               k = k + 1
               dirs(:, k) = real([ia, j, n], wp)/sqrt(3.0_wp)
            end do
         end do
      end do

      stride = max(1, mol%nat/30)
      allocate (cand(ndim, (mol%nat/stride + 1)*nshell*size(dirs, 2)))
      nc = 0
      do ia = 1, mol%nat, stride
         do k = 1, nshell
            do j = 1, size(dirs, 2)
               nc = nc + 1
               cand(:, nc) = mol%xyz(:, ia) + radii(ia)*(1.0_wp + shells(k))*dirs(:, j)
            end do
         end do
      end do
      allocate (b0(nc), st(nc))
      call lsf%evaluate_batch(cand(:, 1:nc), 0, b0, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
         return
      end if

      count_cat = 0
      n = 0
      do k = 1, nc
         if (st(k) /= 0) cycle
         if (b0(k) < 0.0_wp) then
            icat = 1
         else if (b0(k) <= 1.0_wp) then
            icat = 2
         else if (b0(k) > s_cap) then
            icat = 5
         else
            do ia = 1, mol%nat
               sig(ia) = norm2(cand(:, k) - mol%xyz(:, ia))/radii(ia) - 1.0_wp
            end do
            icat = merge(4, 3, any(sig > x_on .and. sig < x_on + delta))
         end if
         if (count_cat(icat) >= per_cat) cycle
         count_cat(icat) = count_cat(icat) + 1
         n = n + 1
         buf(:, n) = cand(:, k)
      end do
      if (any(count_cat == 0)) then
         call test_failed(error, "certificate points: a category (inside, near, mid, band, "// &
                          "far) is empty")
         return
      end if
      points = buf(:, 1:n)
   end subroutine cfc_certificate_points

   !> CFC exclusion radius gates
   !>
   !> - fade off: LiH CPCM radii 2.46 and 3.87 Bohr, `L` = 8.4/Bohr, so 0.119 Bohr
   !>   per unit of `|S|` (catches radius, exponent or unit errors), linear in
   !>   `|S|`; NaN, `-Inf` and `+Inf` certify nothing, `far_exclusion_radius` is 0
   !> - fade on, defaults: `S = -1` certifies `1/Lam(0)`, within 1e-6 below the
   !>   fade-off `1/L`, linear in `|S|`; NaN and `-Inf` nothing, `+Inf` gives
   !>   `r_cap > 0`; far radius positive; both zero once `param` changed after
   !>   `update` (stale constants)
   !> - parameter guards, each with the fade on: `2|a2| <= |a1|` (pair term
   !>   outlives the atomic one: no reached certificate, the far one stays),
   !>   odd `m` and `c < 0` (`PD` may go negative), `a1 > 0` (growing atomic
   !>   term), `|a2| > |a1|` (conservative band guard; the fade-off bound of the
   !>   same parameters certifies, the control), `m > |a2|` (pair envelope not
   !>   decreasing), `m = 2` (not the kernel's pair power 4) and band `a = 400`
   !>   (`alpha (x_on + delta)` = 672 > 600, exponentials underflow): both radii
   !>   zero; controls: `m = 4` and `a = 27` at otherwise the same parameters
   !>   certify both
   subroutine test_cfc_exclusion_radius_gate(error)
      type(error_type), allocatable, intent(out) :: error

      !> Guard cases: a1, a2, c, m, fade_start, expected far radius positive
      integer, parameter :: ncase = 8
      real(wp), parameter :: guard(3, ncase) = reshape([ &
                             -15.0_wp, -7.0_wp, 5.0_wp, &
                             -15.0_wp, -9.0_wp, 5.0_wp, &
                             -15.0_wp, -9.0_wp, -5.0_wp, &
                             1.0_wp, -9.0_wp, 5.0_wp, &
                             -8.0_wp, -9.0_wp, 5.0_wp, &
                             -15.0_wp, -3.5_wp, 5.0_wp, &
                             -15.0_wp, -9.0_wp, 5.0_wp, &
                             -15.0_wp, -9.0_wp, 5.0_wp], [3, ncase])
      integer, parameter :: guard_m(ncase) = [4, 3, 4, 4, 4, 4, 2, 4]
      real(wp), parameter :: guard_a(ncase) = [27.0_wp, 27.0_wp, 27.0_wp, 27.0_wp, 27.0_wp, &
                                               27.0_wp, 27.0_wp, 400.0_wp]
      logical, parameter :: guard_far(ncase) = [.true., .false., .false., .false., .false., &
                                                .false., .false., .false.]
      character(len=*), parameter :: guard_name(ncase) = [character(len=14) :: &
                                     "2|a2| < |a1|", "odd m", "c < 0", "a1 > 0", &
                                     "|a2| > |a1|", "m > |a2|", "m = 2", "a = 400"]
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: cfc, off
      real(wp) :: r, r_off
      integer :: i

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)

      call off%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp, a2=-9.0_wp, &
                                                             c=5.0_wp, m=4, fade_enabled=.false.))
      call off%update(mol, radii)
      r_off = off%exclusion_radius(-1.0_wp)
      call check(error, r_off > 0.05_wp .and. r_off < 0.2_wp, &
                 more="fade off: default parameters must certify ~0.1 Bohr per unit of |S|")
      if (allocated(error)) return
      call check(error, off%exclusion_radius(-3.0_wp), 3.0_wp*r_off, thr=0.0_wp, &
                 message="fade off: the certified radius must scale with |S|")
      if (allocated(error)) return
      call check_nonfinite(off, .false.)
      if (allocated(error)) return
      call check(error, off%far_exclusion_radius(), 0.0_wp, thr=0.0_wp, &
                 message="fade off: no far certificate")
      if (allocated(error)) return

      call cfc%new()
      call cfc%update(mol, radii)
      r = cfc%exclusion_radius(-1.0_wp)
      call check(error, r <= r_off .and. r > (1.0_wp - 1.0e-6_wp)*r_off, &
                 more="fade on: S = -1 must certify 1/Lam(0), just below 1/L")
      if (allocated(error)) return
      call check(error, cfc%exclusion_radius(-3.0_wp), 3.0_wp*r, thr=0.0_wp, &
                 message="fade on: S <= 0 radius must scale with |S|")
      if (allocated(error)) return
      call check_nonfinite(cfc, .true.)
      if (allocated(error)) return
      call check(error, cfc%far_exclusion_radius() > 0.0_wp, &
                 more="fade on: the far certificate must be positive at the defaults")
      if (allocated(error)) return
      cfc%param%fade_start = 20.0_wp
      call check(error, cfc%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="fade on: stale constants must certify nothing")
      if (allocated(error)) return
      call check(error, cfc%far_exclusion_radius(), 0.0_wp, thr=0.0_wp, &
                 message="fade on: stale constants must certify nothing (far)")
      if (allocated(error)) return

      do i = 1, ncase
         call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=guard(1, i), &
                                                                a2=guard(2, i), c=guard(3, i), m=guard_m(i), &
                                                                fade_start=guard_a(i)))
         call cfc%update(mol, radii)
         call check(error, cfc%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                    message="fade on, "//trim(guard_name(i))//": must certify nothing")
         if (allocated(error)) return
         call check(error, cfc%far_exclusion_radius() > 0.0_wp .eqv. guard_far(i), &
                    more="fade on, "//trim(guard_name(i))//": far certificate on/off wrong")
         if (allocated(error)) return
      end do

      ! Control: |a2| > |a1| is a fade-only guard, the fade-off bound certifies
      call off%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-8.0_wp, a2=-9.0_wp, &
                                                             c=5.0_wp, m=4, fade_enabled=.false.))
      call off%update(mol, radii)
      call check(error, off%exclusion_radius(-1.0_wp) > 0.0_wp, &
                 more="control: fade off at |a2| > |a1| must certify")
      if (allocated(error)) return

      ! Controls of the m = 2 and a = 400 rows
      call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp, a2=-9.0_wp, &
                                                             c=5.0_wp, m=4, fade_start=27.0_wp))
      call cfc%update(mol, radii)
      call check(error, cfc%exclusion_radius(-1.0_wp) > 0.0_wp .and. &
                 cfc%far_exclusion_radius() > 0.0_wp, &
                 more="control: m = 4, a = 27 must certify both radii")

   contains

      !> NaN and `-Inf` certify nothing; `+Inf` gives `r_cap` (fade on) or nothing
      !>
      !> @param[in] lsf   CFC instance
      !> @param[in] faded Fade on
      subroutine check_nonfinite(lsf, faded)
         !> CFC instance
         type(moist_cavity_drop_lsf_cfc_type), intent(in) :: lsf
         !> Fade on
         logical, intent(in) :: faded

         call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_quiet_nan)), 0.0_wp, &
                    thr=0.0_wp, message="NaN must certify nothing")
         if (allocated(error)) return
         call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_negative_inf)), &
                    0.0_wp, thr=0.0_wp, message="-Inf must certify nothing")
         if (allocated(error)) return
         if (faded) then
            call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_positive_inf)), &
                       lsf%cert%r_cap, thr=0.0_wp, message="fade on: +Inf must give r_cap")
            if (allocated(error)) return
            call check(error, lsf%cert%r_cap > 0.0_wp, more="fade on: r_cap must be positive")
         else
            call check(error, lsf%exclusion_radius(ieee_value(1.0_wp, ieee_positive_inf)), &
                       0.0_wp, thr=0.0_wp, message="fade off: +Inf must certify nothing")
         end if
      end subroutine check_nonfinite

   end subroutine test_cfc_exclusion_radius_gate

   !> Certificate constants against an independent reference implementation
   !>
   !> - defaults (`a1, a2, c, m` = -15, -9, 5, 4; band 27/3), S_4, cubic grid of
   !>   carbons 5 Bohr apart, R = 2, n = 16 and 1000: `xi*`, far radius, `r_cap`,
   !>   `r(S)` at S = -1, 2, 5, 10 within 1e-9 relative, `s_cap` within 1e-6 (flat
   !>   maximum) of `validate.py`'s `Model` (numpy; slope maxima by a bounded
   !>   optimizer raised 1e-10, as here)
   !> - `S = huge` gives `r_cap`
   !> - negative control: at n = 1000, S = 10 the faded radius sits > 1e-4
   !>   relative below the fade-off `|S|/L`, so the comparison resolves the fade terms
   subroutine test_cfc_certificate_constants(error)
      type(error_type), allocatable, intent(out) :: error

      !> Reference rows: xi*, far, r_cap, r(-1), r(2), r(5), r(10), s_cap
      real(wp), parameter :: ref(8, 2) = reshape([ &
                             0.646920171670627_wp, 5.37282632332541_wp, 2.20769306790104_wp, &
                             0.068891194350184_wp, 0.137782377840153_wp, 0.344455840474173_wp, &
                             0.688909790273216_wp, 45.7229717007_wp, &
                             1.21265586107333_wp, 4.24135494452001_wp, 0.815828471196347_wp, &
                             0.0394144776650353_wp, 0.0788271422162369_wp, &
                             0.197050475260116_wp, 0.393786143953053_wp, 42.1461910903_wp], [8, 2])
      integer, parameter :: ns(2) = [16, 1000]
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: lsf, off
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: xyz(ndim, 1000), got(8)
      integer :: num(1000), i, irow, n

      do i = 1, 1000
         xyz(:, i) = 5.0_wp*real([mod(i - 1, 10), mod((i - 1)/10, 10), (i - 1)/100], wp)
      end do
      num = 6
      do irow = 1, 2
         n = ns(irow)
         call new(mol, num(1:n), xyz(:, 1:n))
         call lsf%new()
         call lsf%update(mol, [(2.0_wp, i=1, n)])
         got = [lsf%cert%xi_star, lsf%far_exclusion_radius(), lsf%cert%r_cap, &
                lsf%exclusion_radius(-1.0_wp), lsf%exclusion_radius(2.0_wp), &
                lsf%exclusion_radius(5.0_wp), lsf%exclusion_radius(10.0_wp), lsf%cert%s_cap]
         call check(error, maxval(abs(got(1:7)/ref(1:7, irow) - 1.0_wp)) < 1.0e-9_wp, &
                    more="certificate constants differ from validate.py (n = "// &
                    trim(merge("16  ", "1000", irow == 1))//")")
         if (allocated(error)) return
         call check(error, abs(got(8)/ref(8, irow) - 1.0_wp) < 1.0e-6_wp, &
                    more="s_cap differs from validate.py")
         if (allocated(error)) return
         call check(error, lsf%exclusion_radius(huge(1.0_wp)), lsf%cert%r_cap, thr=0.0_wp, &
                    more="S = huge must give r_cap")
         if (allocated(error)) return
      end do

      ! Negative control: n = 1000, the fade terms move r(10) off |S|/L
      call param%new(fade_enabled=.false.)
      call off%new(param)
      call off%update(mol, [(2.0_wp, i=1, n)])
      call check(error, abs(lsf%exclusion_radius(10.0_wp)/off%exclusion_radius(10.0_wp) - 1.0_wp) &
                 > 1.0e-4_wp, more="control: faded and fade-off radii must differ at n = 1000")
   end subroutine test_cfc_certificate_constants

   !> Far certificate: no surface within `far_exclusion_radius()` of an unreached point
   !>
   !> - MB16-43 01 and a 512-atom lattice, default fade: random points in the box
   !>   padded by the reach + 8 Bohr, kept when `evaluate_batch` reports
   !>   the column unreached (no atom within the reach); each ball probed by
   !>   [[ball_sign_change]] (14 directions x 4 fractions, 64 random samples, a
   !>   probe no atom reaches counting as exterior)
   !> - control: carbon dimer (R = 2, 3 Bohr apart), x on the axis 0.01 Bohr past
   !>   the reach of the near atom: unreached; marching from x toward the dimer
   !>   finds the surface at `d_ray`, and `0 < far_exclusion_radius() < d_ray`
   !>   (the certificate is positive and below the true distance)
   !> - `xi* > 0` and strictly increasing over n = 2, 16, 1000 (cubic carbon
   !>   grid 5 Bohr apart, R = 2): `F_n` grows with `n`
   !> - loose on real geometry by construction: an unreached point sits more
   !>   than `R (1 + x_end)` from every centre, the surface within `R (1 + xi*)`
   !>   of some centre, so the ball probes cannot see an over-claim of up to
   !>   `R xi*` (the review's `xi*`-dropping mutation passes them);
   !>   `cfc_certificate_constants` pins the far radius itself
   subroutine test_cfc_far_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncand = 400, nkeep = 24
      !> Atom counts of the `xi*` property
      integer, parameter :: ns(3) = [2, 16, 1000]
      type(structure_type) :: mols(2), dimer
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)
      real(wp) :: cand(ndim, ncand), b0(ncand), u(ndim), lo(ndim), hi(ndim), pad, r
      real(wp) :: x(ndim, 1), xb(1), probe(ndim), s, d_ray
      integer :: st(ncand), stx(1), imol, i, nfar
      logical :: flipped
      integer :: state
      type(mctc_error), allocatable :: lsf_err
      !> Concrete instance and lattice of the `xi*` property
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      type(structure_type) :: grid
      real(wp) :: xyz(ndim, 1000), xi(3)
      integer :: num(1000), irow

      do i = 1, 1000
         xyz(:, i) = 5.0_wp*real([mod(i - 1, 10), mod((i - 1)/10, 10), (i - 1)/100], wp)
      end do
      num = 6
      do irow = 1, size(ns)
         call new(grid, num(1:ns(irow)), xyz(:, 1:ns(irow)))
         call cfc%new()
         call cfc%update(grid, [(2.0_wp, i=1, ns(irow))])
         xi(irow) = cfc%cert%xi_star
      end do
      call check(error, xi(1) > 0.0_wp .and. xi(2) > xi(1) .and. xi(3) > xi(2), &
                 more="xi* must be positive and strictly increasing over n = 2, 16, 1000")
      if (allocated(error)) return

      call get_structure(mols(1), "MB16-43", "01")
      call lattice_structure(mols(2), 8, 3.0_wp, 23)
      state = 3
      do imol = 1, size(mols)
         call get_test_radii(mols(imol), radii)
         call init_lsf(lsf, mols(imol), radii, 0, "cfc", fade_enabled=.true.)
         r = lsf%far_exclusion_radius()
         call check(error, r > 0.0_wp, more="far certificate must be positive at the defaults")
         if (allocated(error)) return
         pad = lsf%screening_offset(maxval(radii)) + maxval(radii) + 8.0_wp
         lo = lsf%center_lower - pad
         hi = lsf%center_upper + pad
         do i = 1, ncand
            call det_rand(u, state)
            cand(:, i) = lo + (hi - lo)*u
         end do
         call lsf%evaluate_batch(cand, 0, b0, error=lsf_err)
         call batch_status(lsf, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         nfar = 0
         do i = 1, ncand
            if (st(i) /= 1) cycle
            nfar = nfar + 1
            if (nfar > nkeep) exit
            call ball_sign_change(lsf, cand(:, i), 1.0_wp, r, 64, state, flipped, error)
            if (allocated(error)) return
            if (flipped) then
               call test_failed(error, "far exclusion radius over-claims: surface inside "// &
                                "the ball around an unreached point")
               return
            end if
         end do
         call check(error, nfar > 0, more="no unreached point sampled")
         if (allocated(error)) return
      end do

      ! Control: positive and below the distance to the surface along a ray
      call new(dimer, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.0_wp, 0.0_wp], &
                                      [ndim, 2]))
      call init_lsf(lsf, dimer, [2.0_wp, 2.0_wp], 0, "cfc", fade_enabled=.true.)
      x(:, 1) = [-(2.0_wp + lsf%screening_offset(2.0_wp)) - 0.01_wp, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(x, 0, xb, error=lsf_err)
      call batch_status(lsf, stx)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, stx(1), 1, more="control: x must be unreached")
      if (allocated(error)) return
      d_ray = 0.0_wp
      do i = 1, 20000
         d_ray = 1.0e-3_wp*real(i, wp)
         probe = x(:, 1) + [d_ray, 0.0_wp, 0.0_wp]
         call lsf%prepare(probe, lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "LSF prepare failed: "//lsf_err%message)
            return
         end if
         if (lsf%unreached()) cycle
         call lsf%f0(s)
         if (s <= 0.0_wp) exit
      end do
      call check(error, s <= 0.0_wp, more="control: the ray march must reach the surface")
      if (allocated(error)) return
      r = lsf%far_exclusion_radius()
      call check(error, r > 0.0_wp .and. r < d_ray, &
                 more="control: far radius must be positive and below the distance")
   end subroutine test_cfc_far_exclusion_radius

   !* ================================================================================= *!
   !*                                   SSES-only                                       *!
   !* ================================================================================= *!

   !> Construct a small test geometry
   !>
   !> @param[out] lsf Level set
   !> @param[out] mol Molecular geometry
   !> @param[out] radii Bare atom radii
   !> @param[in] n Atom count
   !> @param[in] order Highest spatial derivative order
   subroutine sses_setup(lsf, mol, radii, n, order)
      !> Level set
      type(moist_cavity_drop_lsf_sses_type), intent(out) :: lsf
      !> Molecular geometry
      type(structure_type), intent(out) :: mol
      !> Bare atom radii
      real(wp), allocatable, intent(out) :: radii(:)
      !> Atom count
      integer, intent(in) :: n
      !> Highest spatial derivative order
      integer, intent(in) :: order
      real(wp) :: xyz(3, 3)
      integer :: i

      xyz(:, 1) = [-1.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [2.3_wp, 0.1_wp, 0.0_wp]
      xyz(:, 3) = [1.0_wp, 1.1_wp, 0.3_wp]
      radii = [1.1_wp, 0.8_wp, 0.9_wp]
      radii = radii(1:n)
      call new(mol, [(1, i=1, size(radii))], xyz(:, 1:size(radii)))
      call lsf%new(moist_cavity_drop_lsf_sses_param_type( &
                   probe_radius=0.5_wp, smoothing_length=0.35_wp, scale=1.0_wp))
      call lsf%set_max_deriv(order)
      call lsf%update(mol, radii)
   end subroutine sses_setup

   !> Prepare a point and translate evaluation failures into test failures
   !>
   !> @param[in,out] lsf Level set
   !> @param[in] point Evaluation point
   !> @param[out] error Unit test failure
   subroutine sses_prepare(lsf, point, error)
      !> Level set
      type(moist_cavity_drop_lsf_sses_type), intent(inout) :: lsf
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: failure

      call lsf%prepare(point, failure)
      if (allocated(failure)) call test_failed(error, failure%message)
   end subroutine sses_prepare

   !> Implicit Taylor jet with atom seeds as jet variables, the quadratic-cost reference
   !>
   !> @param[in] lsf Prepared level set
   !> @param[in] mol Molecular geometry
   !> @param[in] radii Bare atom radii
   !> @param[in] point Evaluation point
   !> @param[in] order Highest total derivative order
   !> @param[in] a First atom, seeded in variable 4
   !> @param[in] axis_a Axis of atom a, zero for its radius
   !> @param[in] b Second atom, seeded in variable 5
   !> @param[in] axis_b Axis of atom b, zero for its radius
   !> @param[in] v Nuclear direction, seeded in variable 5
   !> @param[in] vrad Radius direction, seeded in variable 5
   !> @returns jet Seeded Taylor jet
   function sses_implicit_jet(lsf, mol, radii, point, order, a, axis_a, b, axis_b, v, vrad) result(jet)
      !> Prepared level set
      type(moist_cavity_drop_lsf_sses_type), intent(in) :: lsf
      !> Molecular geometry
      type(structure_type), intent(in) :: mol
      !> Bare atom radii
      real(wp), intent(in) :: radii(:)
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> Highest total derivative order
      integer, intent(in) :: order
      !> First atom, seeded in variable 4
      integer, intent(in), optional :: a
      !> Axis of atom a, zero for its radius
      integer, intent(in), optional :: axis_a
      !> Second atom, seeded in variable 5
      integer, intent(in), optional :: b
      !> Axis of atom b, zero for its radius
      integer, intent(in), optional :: axis_b
      !> Nuclear direction, seeded in variable 5
      real(wp), intent(in), optional :: v(:, :)
      !> Radius direction, seeded in variable 5
      real(wp), intent(in), optional :: vrad(:)
      type(sses_jet_type) :: jet
      real(wp) :: ps(3, 5), cs(3, size(radii), 5), rs(size(radii), 5)
      integer :: k, spatial
      logical :: valid

      ps = 0.0_wp
      cs = 0.0_wp
      rs = 0.0_wp
      spatial = order
      if (present(a)) spatial = spatial - 1
      if (present(b) .or. present(v) .or. present(vrad)) spatial = spatial - 1
      if (spatial > 0) then
         do k = 1, 3
            ps(k, k) = 1.0_wp
         end do
      end if
      if (present(a)) then
         if (axis_a == 0) then
            rs(a, 4) = 1.0_wp
         else
            cs(axis_a, a, 4) = 1.0_wp
         end if
      end if
      if (present(b)) then
         if (axis_b == 0) then
            rs(b, 5) = 1.0_wp
         else
            cs(axis_b, b, 5) = 1.0_wp
         end if
      end if
      if (present(v)) cs(:, :, 5) = v
      if (present(vrad)) rs(:, 5) = vrad
      call seeded_evaluate(point, mol%xyz, radii, lsf%param%probe_radius, lsf%param%smoothing_length, &
                           lsf%param%temperature(), ps, cs, rs, order, jet, valid)
      if (.not. valid) jet%c = huge(1.0_wp)
   end function sses_implicit_jet

   !> Exact sphere value, spatial jet and radius derivatives
   subroutine test_sses_sphere(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: f, g(3), h(3, 3), t(3, 3, 3), fourth(3, 3, 3, 3)
      real(wp) :: hp(6), tp(10), fp(15), hrp(6, 1)
      real(wp) :: dr(1), gr(3, 1), hr(3, 3, 1), rr(1, 1)
      real(wp) :: point(3), expected(3, 3)
      integer :: i

      call sses_setup(lsf, mol, radii, 1, 4)
      point = mol%xyz(:, 1) + [0.2_wp, 0.3_wp, 0.4_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f3_rrr(f, g, hp, tp)
      call unpack_sym_matrix(hp, h, 3)
      call unpack_sym_tensor3(tp, t, 3)
      call check(error, abs(f - (0.29_wp - radii(1)*radii(1))) < 1.0e-14_wp)
      if (allocated(error)) return
      call check(error, maxval(abs(g - [0.4_wp, 0.6_wp, 0.8_wp])) < 1.0e-14_wp)
      if (allocated(error)) return
      expected = 0.0_wp
      do i = 1, 3
         expected(i, i) = 2.0_wp
      end do
      call check(error, maxval(abs(h - expected)) < 1.0e-14_wp)
      if (allocated(error)) return
      call lsf%f4_rrrr(fp)
      call unpack_sym_tensor4(fp, fourth, 3)
      call check(error, maxval(abs(t)) + maxval(abs(fourth)) < 1.0e-14_wp)
      if (allocated(error)) return
      call lsf%f3_rr_rad(dr, gr, hrp)
      call unpack_pair_batch(hrp, hr, 1)
      call lsf%f2_radrad(rr)
      call check(error, abs(dr(1) + 2.0_wp*radii(1)) + abs(rr(1, 1) + 2.0_wp) &
                 + maxval(abs(gr)) + maxval(abs(hr)) < 1.0e-14_wp)
   end subroutine test_sses_sphere

   !> Compiled value against an independent primal simplex optimization
   subroutine test_sses_reference(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: value

      call sses_setup(lsf, mol, radii, 3, 0)
      call sses_prepare(lsf, [1.1_wp, 0.7_wp, 0.2_wp], error)
      if (allocated(error)) return
      call lsf%f0(value)
      ! mpmath KKT solve of the primal simplex problem with tau = 2 eps r0
      call check(error, abs(value - (-2.045379173485545_wp)) < 1.0e-12_wp)
   end subroutine test_sses_reference

   !> A distant atom of any radius leaves the field and its derivatives unchanged
   subroutine test_sses_separable(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: alone, joint
      type(structure_type) :: water, pair
      real(wp), parameter :: thr = 1.0e-13_wp
      real(wp) :: xyz(3, 4), radii(4), point(3), diff
      real(wp) :: f, g(3), h(6), fj, gj(3), hj(6)
      real(wp) :: n1(3, 3), n2(3, 3, 3), n3(6, 3, 3), j1(3, 4), j2(3, 3, 4), j3(6, 3, 4)
      real(wp) :: r1(3), r2(3, 3), r3(6, 3), k1(4), k2(3, 4), k3(6, 4)
      integer :: sample

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.22_wp]
      xyz(:, 2) = [0.0_wp, 1.43_wp, -0.89_wp]
      xyz(:, 3) = [0.0_wp, -1.43_wp, -0.89_wp]
      xyz(:, 4) = [190.0_wp, 0.0_wp, 0.0_wp]
      ! The distant radius differs strongly from the water radii
      radii = [2.9_wp, 2.3_wp, 2.3_wp, 4.0_wp]
      call new(water, [8, 1, 1], xyz(:, 1:3))
      call new(pair, [8, 1, 1, 54], xyz)
      call alone%new(moist_cavity_drop_lsf_sses_param_type())
      call joint%new(moist_cavity_drop_lsf_sses_param_type())
      call alone%set_max_deriv(3)
      call joint%set_max_deriv(3)
      call alone%update(water, radii(1:3))
      call joint%update(pair, radii)
      do sample = 1, 4
         point = [0.3_wp*real(sample, wp), 1.5_wp + 0.7_wp*real(sample, wp), -0.4_wp*real(sample, wp)]
         call sses_prepare(alone, point, error)
         if (allocated(error)) return
         call sses_prepare(joint, point, error)
         if (allocated(error)) return
         call alone%f012_r(f, g, h)
         call joint%f012_r(fj, gj, hj)
         call alone%f3_rr_rA(n1, n2, n3)
         call joint%f3_rr_rA(j1, j2, j3)
         call alone%f3_rr_rad(r1, r2, r3)
         call joint%f3_rr_rad(k1, k2, k3)
         diff = max(abs(f - fj), maxval(abs(g - gj)), maxval(abs(h - hj)))
         diff = max(diff, maxval(abs(n1 - j1(:, 1:3))), maxval(abs(n2 - j2(:, :, 1:3))), &
                    maxval(abs(n3 - j3(:, :, 1:3))))
         diff = max(diff, maxval(abs(r1 - k1(1:3))), maxval(abs(r2 - k2(:, 1:3))), &
                    maxval(abs(r3 - k3(:, 1:3))))
         call check(error, diff < thr*(1.0_wp + abs(f)), "distant atom changes the water field")
         if (allocated(error)) return
         diff = max(maxval(abs(j1(:, 4))), maxval(abs(j2(:, :, 4))), maxval(abs(j3(:, :, 4))), &
                    abs(k1(4)), maxval(abs(k2(:, 4))), maxval(abs(k3(:, 4))))
         call check(error, diff < thr, "distant atom has nonzero derivatives")
         if (allocated(error)) return
      end do
   end subroutine test_sses_separable

   !> SES-like bridging, an open distant gap and invalid smoothing rejection
   subroutine test_sses_geometry(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: failure
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      real(wp) :: xyz(3, 2), f

      xyz = 0.0_wp
      xyz(1, :) = [-1.0_wp, 1.0_wp]
      call new(mol, [1, 1], xyz)
      call lsf%new(moist_cavity_drop_lsf_sses_param_type(probe_radius=0.5_wp, smoothing_length=0.1_wp))
      call lsf%set_max_deriv(0)
      call lsf%update(mol, [1.0_wp, 1.0_wp])
      call sses_prepare(lsf, [0.0_wp, 0.3_wp, 0.0_wp], error)
      if (allocated(error)) return
      call lsf%f0(f)
      call check(error, f < 0.0_wp)
      if (allocated(error)) return
      xyz(1, :) = [-4.0_wp, 4.0_wp]
      call lsf%set_centers(xyz)
      call sses_prepare(lsf, [0.0_wp, 0.0_wp, 0.0_wp], error)
      if (allocated(error)) return
      call lsf%f0(f)
      call check(error, f > 0.0_wp)
      if (allocated(error)) return
      lsf%param%smoothing_length = 0.0_wp
      call lsf%prepare([0.0_wp, 0.0_wp, 0.0_wp], failure)
      call check(error, allocated(failure) .and. lsf%prepared_deriv == -1)
   end subroutine test_sses_geometry

   !> Preserve a narrow regular entrant point through DROP critical switching
   subroutine test_sses_entrant_switch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_math_smoothstep_sigmoid_type) :: critical
      type(structure_type) :: mol
      real(wp) :: centers(3, 2), point(3), a, y, f, g(3), h(6), norm, tau
      y = 0.005_wp
      ! Midplane zero of two unit atoms: the entropy of equal weights adds tau log 2
      tau = 2.0_wp*0.1_wp*temperature_length
      a = sqrt(1.0_wp + 2.0_wp*sqrt(1.01_wp) - 2.0_wp*sqrt(y*y + 0.01_wp) &
               + tau*log(2.0_wp) - y*y)
      centers = 0.0_wp
      centers(1, :) = [-a, a]
      call new(mol, [1, 1], centers)
      call lsf%new(moist_cavity_drop_lsf_sses_param_type(probe_radius=1.0_wp, &
                                                         smoothing_length=0.1_wp, scale=1.0_wp))
      call lsf%update(mol, [1.0_wp, 1.0_wp])
      point = [0.0_wp, y, 0.0_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f012_r(f, g, h)
      call new_smoothstep_sigmoid(critical, param%w_0ls_from, param%w_0ls_to, &
                                 p_hi=param%w_0ls_p, a_hi=param%w_0ls_a, p_lo=param%w_0ls_p, a_lo=param%w_0ls_a)
      norm = sqrt(sum(g*g))
      call check(error, abs(f) < 1.0e-12_wp .and. critical%f0(norm) == 0.0_wp)
      if (allocated(error)) return
      lsf%param%scale = 100.0_wp
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f012_r(f, g, h)
      call check(error, abs(f) < 1.0e-10_wp .and. critical%f0(sqrt(sum(g*g))) == 1.0_wp)
   end subroutine test_sses_entrant_switch

   !> Rigid-motion invariance; without screening a subset prepare keeps every atom
   subroutine test_sses_invariance(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: failure
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: point(3), shift(3), rotation(3, 3), f, g(3), h(3, 3), f2, g2(3), h2(3, 3), hp(6)

      call sses_setup(lsf, mol, radii, 3, 2)
      point = [1.1_wp, 0.7_wp, 0.2_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f012_r(f, g, hp)
      call unpack_sym_matrix(hp, h, 3)
      call lsf%prepare_subset(point, [2], failure)
      call check(error,.not. allocated(failure))
      if (allocated(error)) return
      call lsf%f0(f2)
      call check(error, abs(f - f2) < 1.0e-14_wp .and. lsf%active_count() == 3)
      if (allocated(error)) return
      rotation = reshape([0.0_wp, 0.6_wp, 0.8_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.8_wp, -0.6_wp], [3, 3])
      shift = [2.0_wp, -3.0_wp, 1.0_wp]
      call lsf%set_centers(matmul(rotation, mol%xyz) + spread(shift, 2, 3))
      call sses_prepare(lsf, matmul(rotation, point) + shift, error)
      if (allocated(error)) return
      call lsf%f012_r(f2, g2, hp)
      call unpack_sym_matrix(hp, h2, 3)
      call check(error, abs(f - f2) + maxval(abs(matmul(rotation, g) - g2)) &
                 + maxval(abs(matmul(rotation, matmul(h, transpose(rotation))) - h2)) < 1.0e-11_wp)
   end subroutine test_sses_invariance

   !> Constant scaling preserves normalized geometry and all derivative channels
   subroutine test_sses_field_scale(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: original, scaled
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: f, g(3), h(6), third(10), fourth(15)
      real(wp) :: f2, g2(3), h2(6), t2(10), q2(15)
      real(wp) :: nuclear(3, 3), mixed(3, 3, 3), high(6, 3, 3), nuclear2(3, 3)
      real(wp) :: mixed2(3, 3, 3), high2(6, 3, 3), radial(3), radial2(3)
      real(wp) :: gr(3, 3), gr2(3, 3), hr(6, 3), hr2(6, 3)
      real(wp) :: normalized(3, 3), normalized2(3, 3)
      call sses_setup(original, mol, radii, 3, 4)
      scaled = original
      scaled%param%scale = 100.0_wp
      call sses_prepare(original, [1.1_wp, 0.7_wp, 0.2_wp], error)
      if (allocated(error)) return
      call sses_prepare(scaled, [1.1_wp, 0.7_wp, 0.2_wp], error)
      if (allocated(error)) return
      call original%f3_rrr(f, g, h, third)
      call scaled%f3_rrr(f2, g2, h2, t2)
      call original%f4_rrrr(fourth)
      call scaled%f4_rrrr(q2)
      call check(error, abs(f2/100.0_wp - f) + maxval(abs(g2/100.0_wp - g)) &
                 + maxval(abs(h2/100.0_wp - h)) + maxval(abs(t2/100.0_wp - third)) &
                 + maxval(abs(q2/100.0_wp - fourth)) < 1.0e-12_wp)
      if (allocated(error)) return
      call original%f3_rr_rA(nuclear, mixed, high)
      call scaled%f3_rr_rA(nuclear2, mixed2, high2)
      call original%f3_rr_rad(radial, gr, hr)
      call scaled%f3_rr_rad(radial2, gr2, hr2)
      call check(error, maxval(abs(nuclear2/100.0_wp - nuclear)) &
                 + maxval(abs(mixed2/100.0_wp - mixed)) + maxval(abs(high2/100.0_wp - high)) &
                 + maxval(abs(radial2/100.0_wp - radial)) + maxval(abs(gr2/100.0_wp - gr)) &
                 + maxval(abs(hr2/100.0_wp - hr)) < 1.0e-12_wp)
      if (allocated(error)) return
      call original%normalized_f01_rA(f, normalized)
      call scaled%normalized_f01_rA(f2, normalized2)
      call check(error, abs(f - f2) + maxval(abs(normalized - normalized2)) < 1.0e-12_wp)
   end subroutine test_sses_field_scale

   !> Envelope atom derivatives agree with atom seeds differentiated through the moment equation
   subroutine test_sses_envelope(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      type(sses_jet_type) :: jet
      real(wp), allocatable :: radii(:)
      real(wp), parameter :: thr = 1.0e-10_wp
      real(wp) :: point(3), v(3, 3), vr(3), diff
      real(wp) :: n1(3, 3), n2(3, 3, 3), n3(3, 3, 3, 3), n4(3, 3, 3, 3, 3)
      real(wp) :: r1(3), r2(3, 3), r3(3, 3, 3)
      real(wp) :: nn(3, 3, 3, 3), nn1(3, 3, 3, 3, 3), nn2(3, 3, 3, 3, 3, 3)
      real(wp) :: rn(3, 3, 3), rn1(3, 3, 3, 3), rn2(3, 3, 3, 3, 3)
      real(wp) :: rr(3, 3), rr1(3, 3, 3), rr2(3, 3, 3, 3)
      real(wp) :: h1(3, 3), h2(3, 3, 3), h3(3, 3, 3, 3), hr1(3), hr2(3, 3), hr3(3, 3, 3)
      real(wp) :: t0, t1(3), t2(3, 3), t3(3, 3, 3)
      real(wp) :: pn3(6, 3, 3), pr3(6, 3), pn4(10, 3, 3), pnn2(6, 3, 3, 3, 3), prn2(6, 3, 3, 3)
      real(wp) :: prr2(6, 3, 3), ph3(6, 3, 3), phr3(6, 3), pt2(6), pt3(10)
      integer :: i, j, k, a, b, s, t

      call sses_setup(lsf, mol, radii, 3, 3)
      point = [1.1_wp, 0.7_wp, 0.2_wp]
      v = reshape([0.2_wp, -0.1_wp, 0.3_wp, -0.2_wp, 0.4_wp, 0.1_wp, 0.3_wp, 0.2_wp, -0.4_wp], [3, 3])
      vr = [0.1_wp, -0.2_wp, 0.3_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f3_rr_rA(n1, n2, pn3)
      call lsf%f3_rr_rad(r1, r2, pr3)
      call lsf%f4_rrr_rA(pn4)
      call lsf%f2_rArB(nn)
      call lsf%f3_r_rArB(nn1)
      call lsf%f4_rr_rArB(pnn2)
      call lsf%f2_rA_rad(rn)
      call lsf%f3_r_rA_rad(rn1)
      call lsf%f4_rr_rA_rad(prn2)
      call lsf%f2_radrad(rr)
      call lsf%f3_r_radrad(rr1)
      call lsf%f4_rr_radrad(prr2)
      call lsf%hvp_f1_rA(v, h1, vrad=vr)
      call lsf%hvp_f2_r_rA(v, h2, vrad=vr)
      call lsf%hvp_f3_rr_rA(v, ph3, vrad=vr)
      call lsf%hvp_f1_rad(v, vr, hr1)
      call lsf%hvp_f2_r_rad(v, vr, hr2)
      call lsf%hvp_f3_rr_rad(v, vr, phr3)
      call lsf%tangent_f0(v, t0)
      call lsf%tangent_f1_r(v, t1)
      call lsf%tangent_f2_rr(v, pt2)
      call lsf%tangent_f3_rrr(v, pt3)
      call unpack_pair_batch(pn3, n3, 3*3)
      call unpack_pair_batch(pr3, r3, 3)
      call unpack_triple_batch(pn4, n4, 3*3)
      call unpack_pair_batch(pnn2, nn2, 3*3*3*3)
      call unpack_pair_batch(prn2, rn2, 3*3*3)
      call unpack_pair_batch(prr2, rr2, 3*3)
      call unpack_pair_batch(ph3, h3, 3*3)
      call unpack_pair_batch(phr3, hr3, 3)
      call unpack_sym_matrix(pt2, t2, 3)
      call unpack_sym_tensor3(pt3, t3, 3)
      diff = 0.0_wp
      do a = 1, 3
         do s = 0, 3
            jet = sses_implicit_jet(lsf, mol, radii, point, 4, a=a, axis_a=s)
            if (s == 0) then
               call compare(r1(a), [4])
               do i = 1, 3
                  call compare(r2(i, a), [i, 4])
                  do j = 1, 3
                     call compare(r3(i, j, a), [i, j, 4])
                  end do
               end do
            else
               call compare(n1(s, a), [4])
               do i = 1, 3
                  call compare(n2(i, s, a), [i, 4])
                  do j = 1, 3
                     call compare(n3(i, j, s, a), [i, j, 4])
                     do k = 1, 3
                        call compare(n4(i, j, k, s, a), [i, j, k, 4])
                     end do
                  end do
               end do
            end if
            jet = sses_implicit_jet(lsf, mol, radii, point, 4, a=a, axis_a=s, v=v, vrad=vr)
            if (s == 0) then
               call compare(hr1(a), [4, 5])
               do i = 1, 3
                  call compare(hr2(i, a), [i, 4, 5])
                  do j = 1, 3
                     call compare(hr3(i, j, a), [i, j, 4, 5])
                  end do
               end do
            else
               call compare(h1(s, a), [4, 5])
               do i = 1, 3
                  call compare(h2(i, s, a), [i, 4, 5])
                  do j = 1, 3
                     call compare(h3(i, j, s, a), [i, j, 4, 5])
                  end do
               end do
            end if
            do b = 1, 3
               do t = 0, 3
                  jet = sses_implicit_jet(lsf, mol, radii, point, 4, a=a, axis_a=s, b=b, axis_b=t)
                  if (s > 0 .and. t > 0) then
                     call compare(nn(s, a, t, b), [4, 5])
                     do i = 1, 3
                        call compare(nn1(i, s, a, t, b), [i, 4, 5])
                        do j = 1, 3
                           call compare(nn2(i, j, s, a, t, b), [i, j, 4, 5])
                        end do
                     end do
                  else if (s > 0 .and. t == 0) then
                     call compare(rn(s, a, b), [4, 5])
                     do i = 1, 3
                        call compare(rn1(i, s, a, b), [i, 4, 5])
                        do j = 1, 3
                           call compare(rn2(i, j, s, a, b), [i, j, 4, 5])
                        end do
                     end do
                  else if (s == 0 .and. t == 0) then
                     call compare(rr(a, b), [4, 5])
                     do i = 1, 3
                        call compare(rr1(i, a, b), [i, 4, 5])
                        do j = 1, 3
                           call compare(rr2(i, j, a, b), [i, j, 4, 5])
                        end do
                     end do
                  end if
               end do
            end do
         end do
      end do
      jet = sses_implicit_jet(lsf, mol, radii, point, 4, v=v)
      call compare(t0, [5])
      do i = 1, 3
         call compare(t1(i), [i, 5])
         do j = 1, 3
            call compare(t2(i, j), [i, j, 5])
            do k = 1, 3
               call compare(t3(i, j, k), [i, j, k, 5])
            end do
         end do
      end do
      call check(error, diff < thr, "envelope atom derivatives match implicit atom seeds")
      if (allocated(error)) return
      ! Closed-form native blocks (order >= 2 prepared) against the jet fallback (order 1)
      call lsf%normalized_f01_rA(t0, n1)
      call lsf%set_max_deriv(1)
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%normalized_f01_rA(t0, r2)
      call check(error, maxval(abs(n1 - r2)) < thr, "native and jet normalized center derivatives agree")
   contains
      !> Accumulate the relative deviation from the reference jet
      subroutine compare(actual, indices)
         !> Accessor value
         real(wp), intent(in) :: actual
         !> Jet variables of the reference derivative
         integer, intent(in) :: indices(:)
         real(wp) :: expected
         expected = jet_derivative(jet, indices)
         diff = max(diff, abs(actual - expected)/(1.0_wp + abs(expected)))
      end subroutine compare
   end subroutine test_sses_envelope

   !> Compare generated spatial kernels with implicit Taylor jets
   subroutine test_sses_generated_spatial(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(sses_cache_type) :: cache
      type(sses_jet_type) :: reference, actual
      real(wp) :: centers(3, 3), radii(3), point(3), ps(3, 5), cs(3, 3, 5), rs(3, 5)
      real(wp) :: probe, epsilon, tau
      integer :: geometry, sample, order, k, last
      logical :: valid

      centers(:, 1) = [-1.0_wp, 0.0_wp, 0.0_wp]
      centers(:, 2) = [2.3_wp, 0.1_wp, 0.0_wp]
      centers(:, 3) = [1.0_wp, 1.1_wp, 0.3_wp]
      ps = 0.0_wp
      cs = 0.0_wp
      rs = 0.0_wp
      do k = 1, 3
         ps(k, k) = 1.0_wp
      end do
      do geometry = 1, 3
         radii = 1.1_wp
         if (geometry == 2) radii = [1.1_wp, 0.8_wp, 0.9_wp]
         if (geometry == 3) radii = [1.1_wp, 1.1_wp, 1.1000000000001_wp]
         do sample = 1, 8
            probe = 0.1_wp + real(mod(sample, 3), wp)*0.6_wp
            epsilon = 0.1_wp + real(mod(sample, 4), wp)*0.3_wp
            point = [-1.3_wp + real(sample, wp)*0.57_wp, 0.7_wp, 0.2_wp]
            if (sample == 8) point = [30.0_wp, -20.0_wp, 10.0_wp]
            ! Mean radius as reference radius, the sampling of the retired global temperature
            tau = 2.0_wp*epsilon*sum(radii)/3.0_wp
            call sses_cache_build(cache, centers, radii, probe, epsilon, tau, valid)
            call check(error, valid, "spatial cache valid")
            if (allocated(error)) return
            call seeded_evaluate(point, centers, radii, probe, epsilon, tau, ps, cs, rs, 4, reference, valid)
            call check(error, valid, "spatial seeded reference valid")
            if (allocated(error)) return
            do order = 0, 4
               call sses_cached_evaluate(cache, point, order, actual, valid)
               call check(error, valid, "generated spatial evaluation valid")
               if (allocated(error)) return
               last = (order + 1)*(order + 2)*(order + 3)*(order + 4)*(order + 5)/120 - 1
               call check(error, maxval(real(abs(actual%c(0:last) - reference%c(0:last)) &
                                             /(1.0_wp + abs(reference%c(0:last))), wp)) < 2.0e-12_wp, &
                          "generated spatial kernels agree with seeded jets")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_sses_generated_spatial

   !> Rebuild caches after parameter, radius and geometry changes
   subroutine test_sses_cache_invalidation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf, fresh
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: point(3), actual, expected
      integer :: m

      call sses_setup(lsf, mol, radii, 3, 2)
      point = [1.1_wp, 0.7_wp, 0.2_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      do m = 1, 3
         select case (m)
         case (1)
            lsf%param%smoothing_length = 0.5_wp
         case (2)
            radii(2) = 1.3_wp
            call lsf%update(mol, radii)
         case default
            mol%xyz(1, 1) = mol%xyz(1, 1) + 0.3_wp
            call lsf%update(mol, radii)
         end select
         call fresh%new(lsf%param)
         call fresh%update(mol, radii)
         call sses_prepare(lsf, point, error)
         if (allocated(error)) return
         call sses_prepare(fresh, point, error)
         if (allocated(error)) return
         call lsf%f0(actual)
         call fresh%f0(expected)
         call check(error, actual, expected, thr=1.0e-13_wp)
         if (allocated(error)) return
      end do
   end subroutine test_sses_cache_invalidation

   !> Native cache pointers are repaired after intrinsic thread-clone assignment
   subroutine test_sses_thread_clone(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: original, clone, fresh
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      real(wp) :: expected, actual, point(3)
      integer :: i
      call sses_setup(original, mol, radii, 3, 2)
      call sses_prepare(original, [1.1_wp, 0.7_wp, 0.2_wp], error)
      if (allocated(error)) return
      clone = original
      call sses_setup(fresh, mol, radii, 3, 2)
      do i = 1, 12
         point = [real(i, wp)*0.2_wp, 0.4_wp, -0.1_wp]
         call sses_prepare(clone, point, error)
         if (allocated(error)) return
         call clone%f0(actual)
         call sses_prepare(original, -point, error)
         if (allocated(error)) return
         call sses_prepare(fresh, point, error)
         if (allocated(error)) return
         call fresh%f0(expected)
         call check(error, actual, expected, thr=2.0e-10_wp)
         if (allocated(error)) return
      end do
   end subroutine test_sses_thread_clone

   !> Build a DROP cavity on an SSES field
   !>
   !> @param[out] cavity DROP cavity
   !> @param[in] mol Molecular geometry
   !> @param[in] radii Bare atom radii
   !> @param[in] scale Field scale
   !> @param[in] ctx Context, outlives the cavity
   !> @param[out] error Unit test failure
   subroutine sses_drop_cavity(cavity, mol, radii, scale, ctx, error)
      !> DROP cavity
      type(cavity_type_drop), intent(out) :: cavity
      !> Molecular geometry
      type(structure_type), intent(in) :: mol
      !> Bare atom radii
      real(wp), intent(in) :: radii(:)
      !> Field scale
      real(wp), intent(in) :: scale
      !> Context, outlives the cavity
      type(moist_context_type), intent(in), target :: ctx
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: failure
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(moist_cavity_drop_lsf_sses_param_type) :: lsf_param
      type(moist_cavity_drop_parameters_type) :: param
      type(radius_type_custom) :: radius_model

      call new_custom_radii_atoms(radii, radius_model, failure)
      if (allocated(failure)) then
         call test_failed(error, failure%message)
         return
      end if
      param%num_leb = 26
      ! Single projection with Newton refinement
      param%proj_level = 1
      lsf_param%scale = scale
      call lsf%new(lsf_param)
      call new_cavity_drop(cavity, radius_model=radius_model, lsf_model=lsf, &
                           param=param, ctx=ctx, error=failure)
      if (.not. allocated(failure)) call cavity%update(mol, error=failure)
      if (allocated(failure)) call test_failed(error, failure%message)
   end subroutine sses_drop_cavity

   !> DROP projects onto the exact SSES sphere, and a rescaled field moves no point
   !>
   !> An absolute projection tolerance on a strongly scaled field is unreachable
   !> in double precision and silently drops points; the scale must not matter.
   subroutine test_sses_drop(error)
      !> Unit test failure
      type(error_type), allocatable, intent(out) :: error
      type(cavity_type_drop) :: sphere, default, scaled
      type(moist_context_type), target :: ctx
      type(structure_type) :: atom, water
      real(wp) :: xyz(3, 3)
      integer :: i

      call new_context(ctx, nthreads=0, verbosity=0)
      xyz = 0.0_wp
      call new(atom, [1], xyz(:, 1:1))
      call sses_drop_cavity(sphere, atom, [2.0_wp], 1.0_wp, ctx, error)
      if (allocated(error)) return
      call check(error, sphere%ngrid > 0, "no projected points on the sphere")
      if (allocated(error)) return
      do i = 1, sphere%ngrid
         call check(error, abs(sum(sphere%xyz(:, i)*sphere%xyz(:, i)) - 4.0_wp) < 1.0e-8_wp, &
                    "projected point off the bare sphere")
         if (allocated(error)) return
      end do
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.22_wp]
      xyz(:, 2) = [0.0_wp, 1.43_wp, -0.89_wp]
      xyz(:, 3) = [0.0_wp, -1.43_wp, -0.89_wp]
      call new(water, [8, 1, 1], xyz)
      call sses_drop_cavity(default, water, [2.9_wp, 2.3_wp, 2.3_wp], 1.0_wp, ctx, error)
      if (allocated(error)) return
      call sses_drop_cavity(scaled, water, [2.9_wp, 2.3_wp, 2.3_wp], 4.0_wp, ctx, error)
      if (allocated(error)) return
      call check(error, default%ngrid > 0 .and. default%ngrid == scaled%ngrid, &
                 "field scale changes the projected point count")
      if (allocated(error)) return
      call check(error, maxval(abs(default%xyz(:, 1:default%ngrid) - scaled%xyz(:, 1:scaled%ngrid))) < 1.0e-6_wp, &
                 "field scale moves projected points")
   end subroutine test_sses_drop

   !* ================================================================================= *!
   !*                          Batched evaluation (strip kernel)                        *!
   !* ================================================================================= *!

   !> SvdW `evaluate_batch` vs the per-point accessors, orders 0..4
   !>
   !> - MB16-43 and a 1000-atom lattice, production screening, default and
   !>   legacy weights (`include_q` off and on)
   !> - random points inside the atom box and up to 3 Bohr outside it
   !> - batch sizes 1, width, width + 1 and 37 per order (ragged last strips)
   !> - packed `lsf2_rr` / `lsf3_rrr` / `lsf4_rrrr` against the packed accessors,
   !>   slot for slot
   subroutine test_svdw_batch_vs_accessors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mols(2)
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      integer :: imol, iset, order, isize, n
      integer :: state
      integer :: sizes(4)

      call get_structure(mols(1), "MB16-43", "01")
      call lattice_structure(mols(2), 10, 3.0_wp, 41)
      state = 20261103
      do imol = 1, size(mols)
         call get_test_radii(mols(imol), radii)
         do iset = 1, 2
            if (iset == 1) then
               call init_lsf(lsf, mols(imol), radii, 4, "svdw", screening_threshold=1.0e-11_wp)
            else
               call init_lsf(lsf, mols(imol), radii, 4, "svdw", blend_k=svdw_legacy_blend_k, &
                             blend_2b=svdw_legacy_blend_2b, blend_3b=svdw_legacy_blend_3b, &
                             screening_threshold=1.0e-11_wp)
            end if
            do order = 0, 4
               sizes = [1, strip_width_svdw(order), strip_width_svdw(order) + 1, 37]
               do isize = 1, size(sizes)
                  n = sizes(isize)
                  call check_batch_vs_accessors(error, lsf, order, n, 3.0_wp, state)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_svdw_batch_vs_accessors

   !> One random batch against `prepare` plus the per-point accessors
   !>
   !> @param[out]   error  Test failure, naming the order, column and tensor
   !> @param[inout] lsf    SvdW instance after `update`
   !> @param[in]    order  Batch order
   !> @param[in]    n      Batch size
   !> @param[in]    pad    Box padding beyond the atom enclosure (Bohr)
   !> @param[inout] state  LCG state
   subroutine check_batch_vs_accessors(error, lsf, order, n, pad, state)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> SvdW instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Batch order
      integer, intent(in) :: order
      !> Batch size
      integer, intent(in) :: n
      !> Box padding
      real(wp), intent(in) :: pad
      !> LCG state
      integer, intent(inout) :: state

      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      integer :: st(n), i
      type(mctc_error), allocatable :: lsf_err

      lo = lsf%center_lower - pad
      hi = lsf%center_upper + pad
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
      end do
      call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
         return
      end if
      if (any(st /= 0)) then
         call test_failed(error, "batch reports an unreached column inside the reach")
         return
      end if
      call compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4)
   end subroutine check_batch_vs_accessors

   !> Batch outputs against `prepare` plus the per-point accessors, column by column
   !>
   !> Normwise per tensor: `max|batch - point| <= max(BATCH_ABS, BATCH_REL max|point|)`
   !>
   !> @param[out]   error   Test failure
   !> @param[inout] lsf     LSF instance
   !> @param[in]    order   Batch order
   !> @param[in]    points  Batch points [ndim, n]
   !> @param[in]    b0      Batch values
   !> @param[in]    b1      Batch gradients
   !> @param[in]    b2      Batch packed Hessians
   !> @param[in]    b3      Batch packed third derivatives
   !> @param[in]    b4      Batch packed fourth derivatives
   !> @param[in]    skip    Columns not to compare (optional)
   subroutine compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4, skip)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Batch order
      integer, intent(in) :: order
      !> Batch points
      real(wp), intent(in) :: points(:, :)
      !> Batch values
      real(wp), intent(in) :: b0(:)
      !> Batch gradients
      real(wp), intent(in) :: b1(:, :)
      !> Batch packed Hessians
      real(wp), intent(in) :: b2(:, :)
      !> Batch packed third derivatives
      real(wp), intent(in) :: b3(:, :)
      !> Batch packed fourth derivatives
      real(wp), intent(in) :: b4(:, :)
      !> Columns not to compare
      logical, intent(in), optional :: skip(:)

      real(wp) :: f0, f1(ndim), f2(6), f3(10), f4(15)
      integer :: i
      character(len=96) :: where
      type(mctc_error), allocatable :: lsf_err

      call lsf%set_max_deriv(order)
      do i = 1, size(points, 2)
         if (present(skip)) then
            if (skip(i)) cycle
         end if
         call lsf%prepare(points(:, i), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "prepare failed: "//lsf_err%message)
            return
         end if
         select case (order)
         case (0)
            call lsf%f012_r(lsf0=f0)
         case (1)
            call lsf%f012_r(lsf0=f0, lsf1_r=f1)
         case (2)
            call lsf%f012_r(lsf0=f0, lsf1_r=f1, lsf2_rr=f2)
         case default
            call lsf%f3_rrr(lsf0=f0, lsf1_r=f1, lsf2_rr=f2, lsf3_rrr=f3)
            if (order == 4) call lsf%f4_rrrr(f4)
         end select
         write (where, "(a,i0,a,i0,a,i0)") "order ", order, ", column ", i, " of ", &
            size(points, 2)
         call batch_close(error, [b0(i)], [f0], "value, "//trim(where))
         if (allocated(error)) return
         if (order >= 1) call batch_close(error, b1(:, i), f1, "gradient, "//trim(where))
         if (allocated(error)) return
         if (order >= 2) then
            call batch_close(error, b2(:, i), f2, "Hessian, "//trim(where))
            if (allocated(error)) return
         end if
         if (order >= 3) then
            call batch_close(error, b3(:, i), f3, "third, "//trim(where))
            if (allocated(error)) return
         end if
         if (order >= 4) then
            call batch_close(error, b4(:, i), f4, "fourth, "//trim(where), rel=BATCH_REL4)
            if (allocated(error)) return
         end if
      end do
   end subroutine compare_with_accessors

   !> Normwise closeness of a batch tensor to its per-point reference
   !>
   !> @param[out] error  Test failure with the deviation
   !> @param[in]  got    Batch tensor, flattened
   !> @param[in]  ref    Per-point tensor, flattened
   !> @param[in]  what   Tensor and column, for the message
   !> @param[in]  rel    Relative tolerance, default [[BATCH_REL]] (optional)
   subroutine batch_close(error, got, ref, what, rel)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Batch tensor
      real(wp), intent(in) :: got(:)
      !> Per-point tensor
      real(wp), intent(in) :: ref(:)
      !> Tensor and column
      character(len=*), intent(in) :: what
      !> Relative tolerance
      real(wp), intent(in), optional :: rel

      real(wp) :: dev, scale, tol
      character(len=80) :: msg

      tol = BATCH_REL
      if (present(rel)) tol = rel
      dev = maxval(abs(got - ref))
      scale = maxval(abs(ref))
      if (dev <= max(BATCH_ABS, tol*scale)) return
      write (msg, "(a,es10.3,a,es10.3)") " deviates by ", dev, " at norm ", scale
      call test_failed(error, "batch vs per point, "//what//trim(msg))
   end subroutine batch_close

   !> Shared candidate list present vs absent, bit for bit
   !>
   !> - present: `candidates_in_box` over the whole batch's box, once sorted and
   !>   once reversed (the override sorts it back)
   !> - absent: the override's own per-strip box list
   !> - every list visits the reached atoms in ascending candidate order and the
   !>   extra atoms add exact zeros, so the bits agree
   subroutine test_svdw_batch_candidate_list(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 53
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim)
      real(wp) :: a0(n), a1(ndim, n), a2(6, n), a3(10, n), a4(15, n)
      real(wp) :: c0(n), c1(ndim, n), c2(6, n), c3(10, n), c4(15, n)
      integer, allocatable :: cand(:), rev(:)
      integer :: i, ncand, order, ilist
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call lattice_structure(mol, 10, 3.0_wp, 7)
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 4, "svdw", screening_threshold=1.0e-11_wp)
      allocate (cand(lsf%ncenters), rev(lsf%ncenters))

      ! Points in a corner cube, so the whole-batch box prunes the lattice
      state = 99
      lo = lsf%center_lower + 1.0_wp
      hi = lsf%center_lower + 5.0_wp
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
      end do
      lo = minval(points, dim=2)
      hi = maxval(points, dim=2)
      call lsf%candidates_in_box(lo, hi, cand, ncand)
      call check(error, ncand > 0 .and. ncand < lsf%ncenters, &
                 more="the batch box must prune the lattice")
      if (allocated(error)) return
      rev(1:ncand) = cand(ncand:1:-1)

      do order = 0, 4
         call lsf%evaluate_batch(points, order, a0, a1, a2, a3, a4, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         do ilist = 1, 2
            if (ilist == 1) then
               call lsf%evaluate_batch(points, order, c0, c1, c2, c3, c4, &
                                       candidate_indices=cand(1:ncand), error=lsf_err)
            else
               call lsf%evaluate_batch(points, order, c0, c1, c2, c3, c4, &
                                       candidate_indices=rev(1:ncand), error=lsf_err)
            end if
            if (allocated(lsf_err)) then
               call test_failed(error, "evaluate_batch with a list failed: "//lsf_err%message)
               return
            end if
            if (any(c0 /= a0)) then
               call test_failed(error, "candidate list changes the values")
               return
            end if
            if (order >= 1 .and. any(c1 /= a1)) then
               call test_failed(error, "candidate list changes the gradients")
               return
            end if
            if (order >= 2 .and. any(c2 /= a2)) then
               call test_failed(error, "candidate list changes the Hessians")
               return
            end if
            if (order >= 3 .and. any(c3 /= a3)) then
               call test_failed(error, "candidate list changes the third derivatives")
               return
            end if
            if (order >= 4 .and. any(c4 /= a4)) then
               call test_failed(error, "candidate list changes the fourth derivatives")
               return
            end if
         end do
      end do
   end subroutine test_svdw_batch_candidate_list

   !> Unreached columns: `unreached(column)`, `S = +Inf` and NaN derivatives, no `error`,
   !> neighbours unaffected
   !>
   !> - one strip mixing far points (beyond every reach, column 1 among them)
   !>   with near ones; near columns against the per-point accessors and, at
   !>   [[BATCH_REL]], against the same near points without the far ones (the far
   !>   columns widen the strip box, so plateau atoms run the faded body: same
   !>   model, other loop code)
   !> - a far column holds `+Inf` and every derivative through `order` NaN
   !> - order 0 with only `lsf0`, and `n = 0`, accepted
   !> - `blend_1b = 0`: far columns still unreached with the exterior fill, near
   !>   columns finite
   !> - fade off, screening off: a point ~200 Bohr away has `p1 > 0` but `Z`
   !>   underflows to zero, reported unreached
   !> - negative control: the batch without the far points raises nothing
   subroutine test_svdw_batch_unreached(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 11
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      type(moist_cavity_drop_lsf_svdw_type) :: svdw
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: points(ndim, n), u(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      real(wp) :: c0(n), c1(ndim, n), c2(6, n), c3(10, n), c4(15, n)
      integer :: st(n), i, j, order, iset
      integer, allocatable :: near(:)
      logical :: far(n)
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 4, "svdw", screening_threshold=1.0e-11_wp)

      state = 5
      far = .false.
      far([1, 2, 3, 7, 11]) = .true.
      do i = 1, n
         call det_rand(u, state)
         if (far(i)) then
            points(:, i) = lsf%center_upper + 200.0_wp + 10.0_wp*u
         else
            points(:, i) = lsf%center_lower + (lsf%center_upper - lsf%center_lower)*u
         end if
      end do
      near = pack([(i, i=1, n)], .not. far)

      do order = 0, 4
         st = -1
         call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
         call batch_status(lsf, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "unreached columns must not raise: "//lsf_err%message)
            return
         end if
         if (any(st /= merge(1, 0, far))) then
            call test_failed(error, "unreached(column) must hold exactly on the far columns")
            return
         end if
         call check(error, exterior_columns(order, far, b0, b1, b2, b3, b4), &
                    more="unreached columns must hold +Inf and NaN derivatives")
         if (allocated(error)) return
         call compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4, skip=far)
         if (allocated(error)) return

         call lsf%evaluate_batch(points(:, near), order, c0, c1, c2, c3, c4, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "reached columns alone must not raise: "//lsf_err%message)
            return
         end if
         call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "unreached columns must not raise: "//lsf_err%message)
            return
         end if
         do j = 1, size(near)
            call check_column_close(error, order, c0(j:j), c1(:, j:j), c2(:, j:j), &
                                    c3(:, j:j), c4(:, j:j), near(j), b0, b1, b2, b3, b4, &
                                    "near column next to far ones")
            if (allocated(error)) return
         end do
      end do

      ! Order 0 with only the value, and an empty batch
      st = -1
      call lsf%evaluate_batch(points, 0, b0, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "order 0 with lsf0 alone failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == merge(1, 0, far)) .and. all((b0 > huge(b0)) .eqv. far), &
                 more="order 0 with lsf0 alone: unreached flag or +Inf fill wrong")
      if (allocated(error)) return
      call lsf%evaluate_batch(points(:, 1:0), 4, b0(1:0), b1(:, 1:0), b2(:, 1:0), &
                              b3(:, 1:0), b4(:, 1:0), error=lsf_err)
      call batch_status(lsf, st(1:0))
      if (allocated(lsf_err)) then
         call test_failed(error, "empty batch failed: "//lsf_err%message)
         return
      end if

      ! blend_1b = 0 at the production threshold, then screening off
      do iset = 1, 2
         call param%new(blend_1b=merge(0.0_wp, 1.0_wp, iset == 1), fade_enabled=iset == 1)
         call init_svdw_param(svdw, mol, radii, merge(1.0e-11_wp, 0.0_wp, iset == 1), param)
         do order = 0, 4
            st = -1
            call svdw%evaluate_batch(points, order, b0, b1, b2, b3, b4, &
                                     error=lsf_err)
            call batch_status(svdw, st)
            if (allocated(lsf_err)) then
               call test_failed(error, "unreached columns must not raise: "//lsf_err%message)
               return
            end if
            if (any(st /= merge(1, 0, far))) then
               call test_failed(error, merge("blend_1b = 0  ", "screening off ", iset == 1)// &
                                ": unreached(column) must hold exactly on the far columns")
               return
            end if
            call check(error, exterior_columns(order, far, b0, b1, b2, b3, b4), &
                       more="unreached columns must hold +Inf and NaN derivatives")
            if (allocated(error)) return
            call check(error, finiten(b0(near)), more="near columns must stay finite")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_svdw_batch_unreached

   !> The flagged columns hold the exterior fill: `+Inf` value, NaN derivatives through `order`
   !>
   !> @param[in] order Batch order
   !> @param[in] flag  Columns to test
   !> @param[in] b0    Values
   !> @param[in] b1    Gradients
   !> @param[in] b2    Packed Hessians
   !> @param[in] b3    Packed third derivatives
   !> @param[in] b4    Packed fourth derivatives
   !> @returns   ok    All flagged columns filled
   pure logical function exterior_columns(order, flag, b0, b1, b2, b3, b4) result(ok)
      !> Batch order
      integer, intent(in) :: order
      !> Columns to test
      logical, intent(in) :: flag(:)
      !> Values
      real(wp), intent(in) :: b0(:)
      !> Gradients
      real(wp), intent(in) :: b1(:, :)
      !> Packed Hessians
      real(wp), intent(in) :: b2(:, :)
      !> Packed third derivatives
      real(wp), intent(in) :: b3(:, :)
      !> Packed fourth derivatives
      real(wp), intent(in) :: b4(:, :)

      integer :: i

      ok = .true.
      do i = 1, size(flag)
         if (.not. flag(i)) cycle
         ok = ok .and. b0(i) > huge(b0(i))
         if (order >= 1) ok = ok .and. all(ieee_is_nan(b1(:, i)))
         if (order >= 2) ok = ok .and. all(ieee_is_nan(b2(:, i)))
         if (order >= 3) ok = ok .and. all(ieee_is_nan(b3(:, i)))
         if (order >= 4) ok = ok .and. all(ieee_is_nan(b4(:, i)))
      end do
   end function exterior_columns

   !> A point exactly on a nucleus inside a batch matches the per-point path
   !>
   !> - the strip falls back to the guarded per-point path; every column of that
   !>   strip, the nucleus one included, against the accessors
   !> - the nucleus column finite at every order: the unguarded strip would give
   !>   NaN from `0 * (1/0)`
   subroutine test_svdw_batch_on_nucleus(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 19
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), u(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      integer :: st(n), i, order, iset
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      state = 17
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = mol%xyz(:, 1) + 3.0_wp*(u - 0.5_wp)
      end do
      points(:, 5) = mol%xyz(:, 3)
      points(:, 18) = mol%xyz(:, 7)

      do iset = 1, 2
         if (iset == 1) then
            call init_lsf(lsf, mol, radii, 4, "svdw", screening_threshold=1.0e-11_wp)
         else
            call init_lsf(lsf, mol, radii, 4, "svdw", blend_k=svdw_legacy_blend_k, &
                          blend_2b=svdw_legacy_blend_2b, blend_3b=svdw_legacy_blend_3b, &
                          screening_threshold=1.0e-11_wp)
         end if
         do order = 0, 4
            call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, &
                                    error=lsf_err)
            call batch_status(lsf, st)
            if (allocated(lsf_err)) then
               call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
               return
            end if
            if (any(st /= 0)) then
               call test_failed(error, "on-nucleus batch reports an unreached column")
               return
            end if
            if (.not. (finite0(b0(5)) .and. finite0(b0(18)))) then
               call test_failed(error, "on-nucleus value is not finite")
               return
            end if
            if (order >= 1) then
               if (.not. (finite1(b1(:, 5)) .and. finite1(b1(:, 18)))) then
                  call test_failed(error, "on-nucleus gradient is not finite")
                  return
               end if
            end if
            if (order >= 3) then
               if (.not. (finiten(b3(:, 5)) .and. finiten(b3(:, 18)))) then
                  call test_failed(error, "on-nucleus third derivative is not finite")
                  return
               end if
            end if
            if (order >= 4) then
               if (.not. (finiten(b4(:, 5)) .and. finiten(b4(:, 18)))) then
                  call test_failed(error, "on-nucleus fourth derivative is not finite")
                  return
               end if
            end if
            call compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_svdw_batch_on_nucleus

   !> Looping default `evaluate_batch` at orders 3 and 4: packed accessors
   !>
   !> - SSES: no `evaluate_batch` override, so the base default serves it (CFC
   !>   and SvdW run the strip kernel)
   !> - the default calls `f3_rrr` / `f4_rrrr` and copies the packed slots; against
   !>   the accessors at the batch tolerances, not bit for bit: the
   !>   SSES moment solver warm-starts from the previous point, so a re-prepare
   !>   moves the last bits (measured <= 3e-15)
   !> - no column unreached
   subroutine test_sses_batch_default_packed(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 7
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), u(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      integer :: st(n), i
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 4, "sses", screening_threshold=1.0e-11_wp)
      state = 3
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lsf%center_lower + (lsf%center_upper - lsf%center_lower)*u
      end do
      st = -1
      call lsf%evaluate_batch(points, 4, b0, b1, b2, b3, b4, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0), more="default must leave every column reached")
      if (allocated(error)) return
      call compare_with_accessors(error, lsf, 4, points, b0, b1, b2, b3, b4)
   end subroutine test_sses_batch_default_packed

   !> CFC unreached columns through the strip kernel
   !>
   !> - MB16-43 01, threshold 1e-11 (CFC's reach), far points ~200 Bohr out mixed
   !>   with near ones, column 1 far; orders 0..4
   !> - `unreached(column)` exactly on the far columns, no `error`; `+Inf` value and
   !>   every derivative through `order` NaN there
   !> - near columns against the per-point accessors, and at [[BATCH_REL]] against
   !>   the same near points without the far ones: the far columns put the near
   !>   ones into other sub-strips of other sizes, which moves bits at the FMA
   !>   level (other loop code under FP contraction), not values
   !> - negative control: the far-free batch raises nothing and has no unreached column
   subroutine test_cfc_batch_unreached(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 9
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), u(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      real(wp) :: c0(n), c1(ndim, n), c2(6, n), c3(10, n), c4(15, n)
      integer :: st(n), sn(n), i, j, order
      integer, allocatable :: near(:)
      logical :: far(n)
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, "cfc", screening_threshold=1.0e-11_wp)

      state = 15
      far = .false.
      far([1, 4, 8]) = .true.
      do i = 1, n
         call det_rand(u, state)
         if (far(i)) then
            points(:, i) = lsf%center_upper + 200.0_wp + 10.0_wp*u
         else
            points(:, i) = lsf%center_lower + (lsf%center_upper - lsf%center_lower)*u
         end if
      end do
      near = pack([(i, i=1, n)], .not. far)

      do order = 0, 4
         st = -1
         call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
         call batch_status(lsf, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "unreached columns must not raise: "//lsf_err%message)
            return
         end if
         if (any(st /= merge(1, 0, far))) then
            call test_failed(error, "unreached(column) must hold exactly on the far columns")
            return
         end if
         call check(error, exterior_columns(order, far, b0, b1, b2, b3, b4), &
                    more="unreached columns must hold +Inf and NaN derivatives")
         if (allocated(error)) return
         call compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4, skip=far)
         if (allocated(error)) return

         sn = -1
         call lsf%evaluate_batch(points(:, near), order, c0, c1, c2, c3, c4, &
                                 error=lsf_err)
         call batch_status(lsf, sn(1:size(near)))
         if (allocated(lsf_err)) then
            call test_failed(error, "reached columns alone must not raise: "//lsf_err%message)
            return
         end if
         call check(error, all(sn(1:size(near)) == 0), more="near columns must be reached")
         if (allocated(error)) return
         do j = 1, size(near)
            call check_column_close(error, order, c0(j:j), c1(:, j:j), c2(:, j:j), &
                                    c3(:, j:j), c4(:, j:j), near(j), b0, b1, b2, b3, b4, &
                                    "near column next to far ones")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_cfc_batch_unreached

   !> CFC `evaluate_batch` (strip kernel) vs the per-point accessors, orders 0..4
   !>
   !> - MB16-43 01 at the production threshold 1e-11; 1000-atom lattice at 6 Bohr
   !>   spacing and threshold 1e-6: the per-point reference prepares the nuclear
   !>   families pair by pair, O(n_active**2) per point, and at 3 Bohr / 1e-11
   !>   every point reaches ~1000 atoms (seconds per point at order 4)
   !> - random points inside the atom box and up to 3 Bohr outside it
   !> - batch sizes 1, width, width + 1 and 37 per order (ragged last strips)
   !> - tolerances of [[compare_with_accessors]]; measured <= 3e-14 at orders 0..2
   !>   and the per-point 1-ulp floor at orders 3..4
   subroutine test_cfc_batch_vs_accessors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: thresholds(2) = [1.0e-11_wp, 1.0e-6_wp]
      type(structure_type) :: mols(2)
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      integer :: imol, order, isize, n
      integer :: state
      integer :: sizes(4)

      call get_structure(mols(1), "MB16-43", "01")
      call lattice_structure(mols(2), 10, 6.0_wp, 43)
      state = 20261104
      do imol = 1, size(mols)
         call get_test_radii(mols(imol), radii)
         call init_lsf(lsf, mols(imol), radii, 4, "cfc", screening_threshold=thresholds(imol))
         do order = 0, 4
            sizes = [1, strip_width_cfc(order), strip_width_cfc(order) + 1, 37]
            do isize = 1, size(sizes)
               n = sizes(isize)
               call check_batch_vs_accessors(error, lsf, order, n, 3.0_wp, state)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_cfc_batch_vs_accessors

   !> CFC shared candidate list present vs absent, bit for bit
   !>
   !> - 1000-atom lattice at 6 Bohr spacing, threshold 1e-11, 53 points in a
   !>   corner cube (the batch box prunes the lattice; strips wider than the
   !>   split cap, so sub-strips are cut)
   !> - present: `candidates_in_box` over the whole batch's box, once sorted and
   !>   once reversed (the override sorts it back); absent: per sub-strip box lists
   !> - the sub-strip cut depends on the points only, and the extra atoms of a
   !>   longer list add exact zeros in ascending order, so the bits agree
   subroutine test_cfc_batch_candidate_list(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 53
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim)
      real(wp) :: a0(n), a1(ndim, n), a2(6, n), a3(10, n), a4(15, n)
      real(wp) :: c0(n), c1(ndim, n), c2(6, n), c3(10, n), c4(15, n)
      integer, allocatable :: cand(:), rev(:)
      integer :: i, ncand, order, ilist
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call lattice_structure(mol, 10, 6.0_wp, 7)
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 4, "cfc", screening_threshold=1.0e-11_wp)
      allocate (cand(lsf%ncenters), rev(lsf%ncenters))

      state = 99
      lo = lsf%center_lower + 1.0_wp
      hi = lsf%center_lower + 5.0_wp
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
      end do
      lo = minval(points, dim=2)
      hi = maxval(points, dim=2)
      call lsf%candidates_in_box(lo, hi, cand, ncand)
      call check(error, ncand > 0 .and. ncand < lsf%ncenters, &
                 more="the batch box must prune the lattice")
      if (allocated(error)) return
      rev(1:ncand) = cand(ncand:1:-1)

      do order = 0, 4
         call lsf%evaluate_batch(points, order, a0, a1, a2, a3, a4, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         do ilist = 1, 2
            if (ilist == 1) then
               call lsf%evaluate_batch(points, order, c0, c1, c2, c3, c4, &
                                       candidate_indices=cand(1:ncand), error=lsf_err)
            else
               call lsf%evaluate_batch(points, order, c0, c1, c2, c3, c4, &
                                       candidate_indices=rev(1:ncand), error=lsf_err)
            end if
            if (allocated(lsf_err)) then
               call test_failed(error, "evaluate_batch with a list failed: "//lsf_err%message)
               return
            end if
            do i = 1, n
               if (.not. same_column(order, i, a0, a1, a2, a3, a4, i, c0, c1, c2, c3, c4)) then
                  call test_failed(error, "candidate list changes the batch bits")
                  return
               end if
            end do
         end do
      end do
   end subroutine test_cfc_batch_candidate_list

   !> CFC batch over points scattered across 20 Bohr equals the compact batches
   !>
   !> - MB16-43 01, threshold 1e-11; 8 clusters of 8 points (0.5 Bohr cubes)
   !>   with centers spread over a 20 Bohr cube around the molecule, columns
   !>   interleaved so that every strip holds all clusters: its box diagonal is
   !>   far above the split cap and the strip is cut into sub-strips
   !> - reference: each cluster alone (one compact sub-strip); every order 0..4
   !>   at the batch tolerances, unreached column by column; reached columns also
   !>   against the per-point accessors
   !> - checks the permutation and scatter of the split; on 16 atoms the pair
   !>   waste of an unsplit strip is negligible, its cost is measured by the bench
   subroutine test_cfc_batch_scattered_strips(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: ncl = 8, per = 8, n = ncl*per
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), cluster(ndim, per), centers(ndim, ncl), u(ndim), mid(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      real(wp) :: c0(per), c1(ndim, per), c2(6, per), c3(10, per), c4(15, per)
      integer :: st(n), sc(per), i, k, j, order
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 4, "cfc", screening_threshold=1.0e-11_wp)

      state = 67
      mid = 0.5_wp*(lsf%center_lower + lsf%center_upper)
      do k = 1, ncl
         call det_rand(u, state)
         centers(:, k) = mid + 20.0_wp*(u - 0.5_wp)
      end do
      do i = 1, n
         call det_rand(u, state)
         k = mod(i - 1, ncl) + 1
         points(:, i) = centers(:, k) + 0.5_wp*(u - 0.5_wp)
      end do

      do order = 0, 4
         st = -1
         call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
         call batch_status(lsf, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "scattered batch failed: "//lsf_err%message)
            return
         end if
         call check(error, count(st == 0) > n/2, more="most scattered columns must be reached")
         if (allocated(error)) return
         do k = 1, ncl
            do j = 1, per
               cluster(:, j) = points(:, k + (j - 1)*ncl)
            end do
            sc = -1
            call lsf%evaluate_batch(cluster, order, c0, c1, c2, c3, c4, &
                                    error=lsf_err)
            call batch_status(lsf, sc)
            if (allocated(lsf_err)) then
               call test_failed(error, "compact batch failed: "//lsf_err%message)
               return
            end if
            do j = 1, per
               i = k + (j - 1)*ncl
               if (sc(j) /= st(i)) then
                  call test_failed(error, "scattered and compact batches disagree on unreached")
                  return
               end if
               if (st(i) /= 0) cycle
               call check_column_close(error, order, c0(j:j), c1(:, j:j), c2(:, j:j), &
                                       c3(:, j:j), c4(:, j:j), i, b0, b1, b2, b3, b4, &
                                       "scattered vs compact")
               if (allocated(error)) return
            end do
         end do
         call compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4, &
                                     skip=st /= 0)
         if (allocated(error)) return
      end do
   end subroutine test_cfc_batch_scattered_strips

   !> A point exactly on a nucleus inside a CFC batch
   !>
   !> - both CFC paths give a finite value there but non-finite derivatives (the
   !>   atom's cusp; both kernels divide by `||r - R_A||`): order 0 reached and
   !>   finite; from order 1 the jet is not finite with `S <= 0`, an `error`
   !>   naming column 5 and its cusp, never an unreached column
   !> - order 0: fade off against the accessors, fade on (`a = 5`, `b = 2`)
   !>   against the same batch without the nucleus columns, at the batch tolerances
   !> - lane isolation at order 0, both instances: in one compact strip the other
   !>   columns are bit for bit those of the same strip with the nucleus point
   !>   nudged off it (a per-point fallback for the strip, the design before,
   !>   fails this); negative control: the nudged strip raises nothing at orders 1..4
   subroutine test_cfc_batch_on_nucleus(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 19
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: off, faded
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: points(ndim, n), u(ndim), compact(ndim, n), moved(ndim, n), c(ndim)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      real(wp) :: c0(n), c1(ndim, n), c2(6, n), c3(10, n), c4(15, n)
      integer :: st(n), sc(n), i, j, order
      integer, allocatable :: rest(:), full(:)
      !> Instance of the lane-isolation pass, 1 fade off, 2 fade on
      integer :: inst
      logical :: nucleus(n)
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.false.)
      call init_cfc_param(off, mol, radii, 1.0e-11_wp, param)
      call param%new(fade_enabled=.true., fade_start=5.0_wp, fade_width=2.0_wp)
      call init_cfc_param(faded, mol, radii, 1.0e-11_wp, param)
      state = 19
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = mol%xyz(:, 1) + 1.5_wp*(u - 0.5_wp)
      end do
      points(:, 5) = mol%xyz(:, 1)
      points(:, 12) = mol%xyz(:, 2)
      nucleus = .false.
      nucleus([5, 12]) = .true.
      rest = pack([(i, i=1, n)], .not. nucleus)

      ! Order 0 has no cusp: the nucleus value is a value
      st = -1
      call off%evaluate_batch(points, 0, b0, error=lsf_err)
      call batch_status(off, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "on-nucleus batch at order 0 failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0), more="fade off: no column may be unreached at order 0")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(b0)), &
                 more="every value, the nucleus ones included, must be finite at order 0")
      if (allocated(error)) return
      call compare_with_accessors(error, off, 0, points, b0, b1, b2, b3, b4)
      if (allocated(error)) return

      st = -1
      sc = -1
      call faded%evaluate_batch(points, 0, b0, error=lsf_err)
      call batch_status(faded, st)
      if (.not. allocated(lsf_err)) then
         call faded%evaluate_batch(points(:, rest), 0, c0, error=lsf_err)
         call batch_status(faded, sc(1:size(rest)))
      end if
      if (allocated(lsf_err)) then
         call test_failed(error, "faded on-nucleus batch failed: "//lsf_err%message)
         return
      end if
      if (any(st /= 0) .or. any(sc(1:size(rest)) /= 0)) then
         call test_failed(error, "fade on: no column may be unreached at order 0")
         return
      end if
      call check(error, all(ieee_is_finite(b0)), &
                 more="faded: every value, the nucleus ones included, finite at order 0")
      if (allocated(error)) return
      do j = 1, size(rest)
         call check_column_close(error, 0, c0(j:j), c1(:, j:j), c2(:, j:j), &
                                 c3(:, j:j), c4(:, j:j), rest(j), b0, b1, b2, b3, b4, &
                                 "faded column next to a nucleus")
         if (allocated(error)) return
      end do

      ! From order 1 the first nucleus column is the error, its cusp named
      do order = 1, 4
         do inst = 1, 2
            if (inst == 1) then
               call off%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
            else
               call faded%evaluate_batch(points, order, b0, b1, b2, b3, b4, error=lsf_err)
            end if
            if (.not. allocated(lsf_err)) then
               call test_failed(error, "an on-nucleus column from order 1 must raise "// &
                                merge("(fade off)", "(fade on) ", inst == 1))
               return
            end if
            call check(error, index(lsf_err%message, &
                                    "point 5 has a level-set jet not finite with S <= 0") > 0 &
                       .and. index(lsf_err%message, "(on-nucleus column at order >= 1)") > 0, &
                       more="unexpected on-nucleus message: "//lsf_err%message)
            if (allocated(error)) return
            deallocate (lsf_err)
         end do
      end do

      ! Lane isolation: one compact strip (box diagonal < 2 Bohr, one sub-strip)
      ! with a nucleus at column 5, against the same strip with that point
      ! nudged 1e-3 Bohr off it inside the box; same explicit list, so the
      ! strip-mates are the same lanes of the same kernel call, bit for bit
      c = mol%xyz(:, 1)
      do i = 1, n
         call det_rand(u, state)
         compact(:, i) = c + (u - 0.5_wp)
      end do
      compact(:, 1) = c + 0.45_wp
      compact(:, 2) = c - 0.45_wp
      compact(:, 5) = c
      moved = compact
      moved(:, 5) = c + [1.0e-3_wp, 0.0_wp, 0.0_wp]
      full = [(i, i=1, off%ncenters)]
      do inst = 1, 2
         if (inst == 1) then
            call off%evaluate_batch(compact, 0, b0, candidate_indices=full, error=lsf_err)
            call batch_status(off, st)
            if (.not. allocated(lsf_err)) then
               call off%evaluate_batch(moved, 0, c0, candidate_indices=full, error=lsf_err)
               call batch_status(off, sc)
            end if
         else
            call faded%evaluate_batch(compact, 0, b0, candidate_indices=full, error=lsf_err)
            call batch_status(faded, st)
            if (.not. allocated(lsf_err)) then
               call faded%evaluate_batch(moved, 0, c0, candidate_indices=full, error=lsf_err)
               call batch_status(faded, sc)
            end if
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "lane-isolation batch failed: "//lsf_err%message)
            return
         end if
         if (any(st /= 0) .or. any(sc /= 0)) then
            call test_failed(error, "lane isolation: no column may be unreached at order 0")
            return
         end if
         do j = 1, n
            if (j == 5) cycle
            if (.not. same_column(0, j, b0, b1, b2, b3, b4, j, c0, c1, c2, c3, c4)) then
               call test_failed(error, "a strip-mate of the nucleus column changed its bits "// &
                                merge("(fade off)", "(fade on) ", inst == 1))
               return
            end if
         end do
         do order = 1, 4
            if (inst == 1) then
               call off%evaluate_batch(moved, order, c0, c1, c2, c3, c4, &
                                       candidate_indices=full, error=lsf_err)
               call batch_status(off, sc)
            else
               call faded%evaluate_batch(moved, order, c0, c1, c2, c3, c4, &
                                         candidate_indices=full, error=lsf_err)
               call batch_status(faded, sc)
            end if
            if (allocated(lsf_err)) then
               call test_failed(error, "control: the nudged strip must not raise: "// &
                                lsf_err%message)
               return
            end if
            call check(error, all(sc == 0), more="control: the nudged strip must be reached")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_cfc_batch_on_nucleus

   !> CFC instance from an explicit parameter set, `max_deriv = 4`
   !>
   !> @param[out] lsf       CFC instance after `update`
   !> @param[in]  mol       Structure
   !> @param[in]  radii     Per-atom radii
   !> @param[in]  threshold Screening threshold
   !> @param[in]  param     Shape and fade parameters
   subroutine init_cfc_param(lsf, mol, radii, threshold, param)
      !> CFC instance
      type(moist_cavity_drop_lsf_cfc_type), intent(out) :: lsf
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Screening threshold
      real(wp), intent(in) :: threshold
      !> Shape and fade parameters
      type(moist_cavity_drop_lsf_cfc_param_type), intent(in) :: param

      lsf%screening_threshold = threshold
      call lsf%new(param)
      call lsf%update(mol, radii)
      call lsf%set_max_deriv(4)
   end subroutine init_cfc_param

   !> Sample points off both CFC band edges with one atom inside its band
   !>
   !> - the CFC band is per atom in `sigma_B = ||p - R_B||/R_b - 1`: onset
   !>   `start`, width `width`, i.e. `R_b start` and `R_b width` Bohr in `d_B`
   !> - candidates `R_A + R_a + s` along 14 directions per atom, `s` in 0.4..2.5 Bohr
   !> - kept when every atom B is more than `margin` Bohr (in `d_B`) off both of
   !>   its edges and some atom lies inside its band by more than `margin`
   !>
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  start  Fade onset in sigma
   !> @param[in]  width  Fade band width in sigma
   !> @param[in]  margin Distance kept from both edges (Bohr)
   !> @param[out] pts    Accepted points [ndim, npts]
   subroutine cfc_band_points(mol, radii, start, width, margin, pts)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Fade onset in sigma
      real(wp), intent(in) :: start
      !> Fade band width in sigma
      real(wp), intent(in) :: width
      !> Distance kept from both edges
      real(wp), intent(in) :: margin
      !> Accepted points
      real(wp), allocatable, intent(out) :: pts(:, :)

      real(wp), parameter :: shells(5) = [0.4_wp, 0.8_wp, 1.2_wp, 1.8_wp, 2.5_wp]
      real(wp) :: dirs(ndim, 14), buf(ndim, 14*5*mol%nat), p(ndim), d(mol%nat)
      real(wp) :: x_on(mol%nat), x_end(mol%nat)
      integer :: i, j, k, n, ia

      dirs = 0.0_wp
      do j = 1, ndim
         dirs(j, 2*j - 1) = 1.0_wp
         dirs(j, 2*j) = -1.0_wp
      end do
      k = 6
      do i = -1, 1, 2
         do j = -1, 1, 2
            do n = -1, 1, 2
               k = k + 1
               dirs(:, k) = real([i, j, n], wp)/sqrt(3.0_wp)
            end do
         end do
      end do
      x_on = radii(1:mol%nat)*start
      x_end = radii(1:mol%nat)*(start + width)
      n = 0
      do ia = 1, mol%nat
         do k = 1, size(shells)
            do j = 1, size(dirs, 2)
               p = mol%xyz(:, ia) + (radii(ia) + shells(k))*dirs(:, j)
               d = band_distances(mol, radii, p)
               if (any(abs(d - x_on) <= margin .or. abs(d - x_end) <= margin)) cycle
               if (.not. any(d > x_on + margin .and. d < x_end - margin)) cycle
               n = n + 1
               buf(:, n) = p
            end do
         end do
      end do
      pts = buf(:, 1:n)
   end subroutine cfc_band_points

   !> CFC fade through `evaluate_batch`: fade limit, negative control, FD inside a band
   !>
   !> - fade limit: band `a = 1e5` e-folds (onset ~1.1e4 radii), asserted beyond
   !>   every point and atom, against fade off with screening off (threshold 0,
   !>   so neither instance masks an atom): bit for bit at orders 0..4
   !> - negative control: a biting band (`a = 5`, `b = 2`) moves the value of a
   !>   reached column by more than 1e-6 relative
   !> - FD: the biting band, points from [[cfc_band_points]] (every atom >= 0.05
   !>   Bohr off both edges, one inside its band; stencil `+-2 STEP_SIZE`); order
   !>   L at the point against the 4-point central FD of order L - 1
   !>   on the stencil, L = 1..4, both through `evaluate_batch`; CFC FD
   !>   thresholds
   !> - negative control of the FD: the unfaded gradient misses the faded FD by
   !>   more than 100 times the threshold
   subroutine test_cfc_batch_fade_band(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 37, npts_max = 4
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), pts(:, :)
      type(moist_cavity_drop_lsf_cfc_type) :: off, far, near
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim), smax, h, miss
      real(wp) :: a0(n), a1(ndim, n), a2(6, n), a3(10, n), a4(15, n)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      real(wp) :: point(ndim), sten(ndim, 4*ndim), steps(4)
      real(wp) :: p0(1), p1(ndim, 1), p2(6, 1), p3(10, 1), p4(15, 1)
      real(wp) :: s0(4*ndim), s1(ndim, 4*ndim), s2(6, 4*ndim), s3(10, 4*ndim)
      real(wp) :: s4(15, 4*ndim), fd, ana, f(4)
      integer :: sa(n), sb(n), ss(4*ndim), i, k, order, ipt, axis, istep, ic, a, b, c
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.false.)
      call init_cfc_param(off, mol, radii, 0.0_wp, param)
      call param%new(fade_enabled=.true., fade_start=1.0e5_wp, fade_width=7.0_wp)
      call init_cfc_param(far, mol, radii, 1.0e-11_wp, param)
      call param%new(fade_start=5.0_wp, fade_width=2.0_wp)
      call init_cfc_param(near, mol, radii, 1.0e-11_wp, param)

      ! Fade limit
      state = 37
      lo = off%center_lower - 3.0_wp
      hi = off%center_upper + 3.0_wp
      smax = 0.0_wp
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
         smax = max(smax, maxval(band_distances(mol, radii, points(:, i))/radii(1:mol%nat)))
      end do
      call check(error, 1.0e5_wp/abs(far%param%a2) > smax, &
                 more="the far band must start beyond every atom")
      if (allocated(error)) return
      do order = 0, 4
         call off%evaluate_batch(points, order, a0, a1, a2, a3, a4, error=lsf_err)
         call batch_status(off, sa)
         if (.not. allocated(lsf_err)) then
            call far%evaluate_batch(points, order, b0, b1, b2, b3, b4, &
                                    error=lsf_err)
            call batch_status(far, sb)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         if (any(sa /= 0) .or. any(sb /= 0)) then
            call test_failed(error, "fade limit: a column near the molecule is unreached")
            return
         end if
         do i = 1, n
            if (.not. same_column(order, i, a0, a1, a2, a3, a4, i, b0, b1, b2, b3, b4)) then
               call test_failed(error, "band beyond every atom changes the batch bits")
               return
            end if
         end do
      end do

      ! Negative control of the limit
      call near%evaluate_batch(points, 0, b0, error=lsf_err)
      call batch_status(near, sb)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch with the near band failed: "//lsf_err%message)
         return
      end if
      call off%evaluate_batch(points, 0, a0, error=lsf_err)
      call batch_status(off, sa)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch fade off failed: "//lsf_err%message)
         return
      end if
      call check(error, any(sb == 0 .and. abs(b0 - a0) > 1.0e-6_wp*max(1.0_wp, abs(a0))), &
                 more="a band inside the molecule must change the value")
      if (allocated(error)) return

      ! FD inside the band
      call cfc_band_points(mol, radii, 5.0_wp/abs(near%param%a2), 2.0_wp/abs(near%param%a2), &
                           0.05_wp, pts)
      call check(error, size(pts, 2) >= npts_max, more="too few points with an atom in the band")
      if (allocated(error)) return
      h = STEP_SIZE
      steps = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]
      miss = 0.0_wp
      do ipt = 1, npts_max
         point = pts(:, 1 + (ipt - 1)*(size(pts, 2)/npts_max))
         do axis = 1, ndim
            do istep = 1, 4
               ic = istep + 4*(axis - 1)
               sten(:, ic) = point
               sten(axis, ic) = point(axis) + steps(istep)*h
            end do
         end do
         do order = 1, 4
            call near%evaluate_batch(reshape(point, [ndim, 1]), order, p0, p1, p2, p3, p4, &
                                     error=lsf_err)
            if (.not. allocated(lsf_err)) then
               call near%evaluate_batch(sten, order - 1, s0, s1, s2, s3, s4, &
                                        error=lsf_err)
               call batch_status(near, ss)
            end if
            if (allocated(lsf_err)) then
               call test_failed(error, "evaluate_batch in the band failed: "//lsf_err%message)
               return
            end if
            if (any(ss /= 0)) then
               call test_failed(error, "a stencil point inside the band is unreached")
               return
            end if
            ! component k of order L - 1, differenced along `axis`, is the order-L slot
            do axis = 1, ndim
               do k = 1, cfc_npack_test(order - 1)
                  do istep = 1, 4
                     ic = istep + 4*(axis - 1)
                     select case (order)
                     case (1)
                        f(istep) = s0(ic)
                     case (2)
                        f(istep) = s1(k, ic)
                     case (3)
                        f(istep) = s2(k, ic)
                     case default
                        f(istep) = s3(k, ic)
                     end select
                  end do
                  call fd4_scalar(f(1), f(2), f(3), f(4), h, fd, error)
                  if (allocated(error)) return
                  select case (order)
                  case (1)
                     ana = p1(axis, 1)
                  case (2)
                     ana = p2(packed_index(k, axis), 1)
                  case (3)
                     call packed_pair(k, a, b)
                     ana = p3(packed_index3(a, b, axis), 1)
                  case default
                     call packed_triple(k, a, b, c)
                     ana = p4(packed_index4(a, b, c, axis), 1)
                  end select
                  call check(error, ana, fd, thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
                  if (allocated(error)) then
                     write (error%message, "(a,i0,a,i0,a,i0,a,i0,a,es12.4,a,es12.4)") &
                        "fade FD order ", order, " point ", ipt, " axis ", axis, &
                        " slot ", k, ": analytic ", ana, " FD ", fd
                     return
                  end if
                  if (order == 1) then
                     call off%evaluate_batch(reshape(point, [ndim, 1]), 1, a0(1:1), &
                                             a1(:, 1:1), error=lsf_err)
                     if (allocated(lsf_err)) then
                        call test_failed(error, "unfaded batch failed: "//lsf_err%message)
                        return
                     end if
                     miss = max(miss, abs(a1(axis, 1) - fd))
                  end if
               end do
            end do
         end do
      end do
      call check(error, miss > 100.0_wp*CFC_ABS_THR, &
                 more="the band must change the gradient beyond the FD threshold")
   end subroutine test_cfc_batch_fade_band

   !> Packed component count of a symmetric rank-L spatial tensor, L = 0..3
   !>
   !> @param[in] rank Tensor rank
   !> @returns   n    Packed components
   pure integer function cfc_npack_test(rank) result(n)
      !> Tensor rank
      integer, intent(in) :: rank

      select case (rank)
      case (0)
         n = 1
      case (1)
         n = 3
      case (2)
         n = 6
      case default
         n = 10
      end select
   end function cfc_npack_test

   !> Sorted index pair of packed rank-2 slot `k` (LAPACK upper order)
   !>
   !> @param[in]  k Packed slot 1..6
   !> @param[out] a First index
   !> @param[out] b Second index, `a <= b`
   pure subroutine packed_pair(k, a, b)
      !> Packed slot
      integer, intent(in) :: k
      !> First index
      integer, intent(out) :: a
      !> Second index
      integer, intent(out) :: b

      integer, parameter :: ia(6) = [1, 1, 2, 1, 2, 3], ib(6) = [1, 2, 2, 3, 3, 3]

      a = ia(k)
      b = ib(k)
   end subroutine packed_pair

   !> Sorted index triple of packed rank-3 slot `k` (sorted tuples, last slowest)
   !>
   !> @param[in]  k Packed slot 1..10
   !> @param[out] a First index
   !> @param[out] b Second index
   !> @param[out] c Third index, `a <= b <= c`
   pure subroutine packed_triple(k, a, b, c)
      !> Packed slot
      integer, intent(in) :: k
      !> First index
      integer, intent(out) :: a
      !> Second index
      integer, intent(out) :: b
      !> Third index
      integer, intent(out) :: c

      integer :: i, j, l

      do l = 1, ndim
         do j = 1, l
            do i = 1, j
               if (packed_index3(i, j, l) /= k) cycle
               a = i
               b = j
               c = l
               return
            end do
         end do
      end do
      error stop "packed_triple: slot outside 1..10"
   end subroutine packed_triple

   !> CFC per-point path with the fade on against `evaluate_batch`, orders 0..4
   !>
   !> - MB16-43 01, production threshold; bands `a = 27`, `b = 3` (the shipped
   !>   default) and the biting `a = 5`, `b = 2` e-folds of `u_A`
   !> - 37 random points in the atom box padded by 3 Bohr plus one beyond every
   !>   band end; reached columns at [[BATCH_REL]] / [[BATCH_REL4]] through
   !>   [[compare_with_accessors]] (`prepare` plus the spatial accessors; the
   !>   batch has no per-atom families, `cfc_fade_fd_atom_families` and
   !>   `cfc_fade_contracted` check those)
   !> - unreached column by column: `prepare` reports every column the batch
   !>   leaves unreached through `unreached()`, the far point among them
   !> - negative controls: at both bands some point has an atom inside its band
   !>   (else only the plateau bodies run); at the biting band the unfaded value
   !>   misses the faded batch by more than 1e-6 relative
   subroutine test_cfc_fade_prepare_vs_batch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 38
      !> `[fade_start, fade_width]` per pass, e-folds of `u_A`
      real(wp), parameter :: bands(2, 2) = reshape([27.0_wp, 3.0_wp, 5.0_wp, 2.0_wp], [2, 2])
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), d(:), x_on(:), x_end(:)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf, off
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim), f0, miss
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      integer :: st(n), i, iband, order, n_inband
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.false.)
      call init_cfc_param(off, mol, radii, 1.0e-11_wp, param)
      state = 20261105
      lo = off%center_lower - 3.0_wp
      hi = off%center_upper + 3.0_wp
      do i = 1, n - 1
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
      end do
      points(:, n) = off%center_upper + 50.0_wp

      do iband = 1, size(bands, 2)
         call param%new(fade_enabled=.true., fade_start=bands(1, iband), &
                        fade_width=bands(2, iband))
         call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
         x_on = radii(1:mol%nat)*bands(1, iband)/abs(lsf%param%a2)
         x_end = radii(1:mol%nat)*sum(bands(:, iband))/abs(lsf%param%a2)
         n_inband = 0
         do i = 1, n
            d = band_distances(mol, radii, points(:, i))
            if (any(d > x_on .and. d < x_end)) n_inband = n_inband + 1
         end do
         call check(error, n_inband > 0, more="no point with an atom inside its band")
         if (allocated(error)) return
         call check(error, all(band_distances(mol, radii, points(:, n)) > x_end), &
                    more="the far point must lie beyond every band end")
         if (allocated(error)) return

         do order = 0, 4
            st = -1
            call lsf%evaluate_batch(points, order, b0, b1, b2, b3, b4, &
                                    error=lsf_err)
            call batch_status(lsf, st)
            if (allocated(lsf_err)) then
               call test_failed(error, "faded evaluate_batch failed: "//lsf_err%message)
               return
            end if
            call check(error, count(st == 0) > n/2, more="most columns must be reached")
            if (allocated(error)) return
            call check(error, st(n), 1, more="the far point must be unreached")
            if (allocated(error)) return
            call compare_with_accessors(error, lsf, order, points, b0, b1, b2, b3, b4, &
                                        skip=st /= 0)
            if (allocated(error)) return
            do i = 1, n
               if (st(i) == 0) cycle
               call lsf%prepare(points(:, i), lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "prepare at an unreached column failed: "//lsf_err%message)
                  return
               end if
               call check(error, lsf%unreached(), more="prepare must report the column unreached")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Negative control at the biting band: the order-4 values of the last pass
      call off%set_max_deriv(0)
      miss = 0.0_wp
      do i = 1, n
         if (st(i) /= 0) cycle
         call off%prepare(points(:, i), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "unfaded prepare failed: "//lsf_err%message)
            return
         end if
         call off%f0(f0)
         miss = max(miss, abs(f0 - b0(i))/max(1.0_wp, abs(b0(i))))
      end do
      call check(error, miss > 1.0e-6_wp, more="the biting band must change the value")
   end subroutine test_cfc_fade_prepare_vs_batch

   !> CFC dispatch of [[run_contracted]] with the fade on
   subroutine test_cfc_fade_contracted(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_contracted(error, "cfc", fade=.true.)
   end subroutine test_cfc_fade_contracted

   !> CFC fade width 0 refused by `prepare` (parameter validation)
   subroutine test_bad_cfc_prepare_fade_width(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_width=0.0_wp)
      call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%prepare(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], lsf_err)
      call expect_error(error, lsf_err, "must be positive")
   end subroutine test_bad_cfc_prepare_fade_width

   !> CFC fade band changed after `update` refused by `evaluate_batch` (stale reach)
   subroutine test_bad_cfc_fade_changed(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: b0(1)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_start=5.0_wp, fade_width=2.0_wp)
      call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
      lsf%param%fade_start = 26.0_wp
      call lsf%evaluate_batch(reshape(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]), &
                              0, b0, error=lsf_err)
      call expect_error(error, lsf_err, "CFC parameters changed after update")
   end subroutine test_bad_cfc_fade_changed

   !> CFC fade width 0 refused by `evaluate_batch` (parameter validation)
   subroutine test_bad_cfc_fade_width(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: b0(1)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_width=0.0_wp)
      call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%evaluate_batch(reshape(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]), &
                              0, b0, error=lsf_err)
      call expect_error(error, lsf_err, "must be positive")
   end subroutine test_bad_cfc_fade_width

   !> CFC atomic exponent past the strip kernel's inline exp range: `evaluate_batch` refuses
   !>
   !> - `a1 = -800`, so `max(-a1, -a2) > 708`; the per-point path is not limited by it
   subroutine test_bad_cfc_exponent_range(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: b0(1)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(a1=-800.0_wp)
      call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%evaluate_batch(reshape(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]), &
                              0, b0, error=lsf_err)
      call expect_error(error, lsf_err, "max(-a1, -a2) <= 708")
   end subroutine test_bad_cfc_exponent_range

   !> SvdW `blend_k = 0` with the fade off: `evaluate_batch` refuses
   !>
   !> - the fade-off model reads no band, but `blend_k` is validated whatever the fade
   subroutine test_bad_svdw_blend_k_fade_off(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: b0(1)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(blend_k=0.0_wp, fade_enabled=.false.)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%evaluate_batch(reshape(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]), &
                              0, b0, error=lsf_err)
      call expect_error(error, lsf_err, "SvdW blend_k must be finite and positive")
   end subroutine test_bad_svdw_blend_k_fade_off

   !> CFC reach refreshed under another band, band restored: `evaluate_batch` refuses
   !>
   !> - `fade_start` 26 -> 5, `set_screening_threshold` rebuilds the reach from
   !>   it, `fade_start` back to 26: the recorded band matches `update`'s, the
   !>   cached reach does not
   subroutine test_bad_cfc_stale_reach(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: b0(1)
      type(mctc_error), allocatable :: lsf_err

      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp], [ndim, 2]))
      call param%new(fade_enabled=.true., fade_start=26.0_wp, fade_width=7.0_wp)
      call init_cfc_param(lsf, mol, [3.0_wp, 3.0_wp], 1.0e-11_wp, param)
      lsf%param%fade_start = 5.0_wp
      call lsf%set_screening_threshold(1.0e-10_wp)
      lsf%param%fade_start = 26.0_wp
      call lsf%evaluate_batch(reshape([7.6_wp, 0.0_wp, 0.0_wp], [ndim, 1]), 0, b0, error=lsf_err)
      call expect_error(error, lsf_err, "CFC parameters changed after update")
   end subroutine test_bad_cfc_stale_reach

   !> CFC fade off, `c` changed after `update`: `prepare` refuses
   !>
   !> - fade-off reach: the threshold reach of `a1, a2, c, m`
   !> - `c` 5 -> 50 moves the pair offset, the cached reach goes stale
   subroutine test_bad_cfc_stale_reach_off(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      type(mctc_error), allocatable :: lsf_err

      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 0.0_wp, 0.0_wp], [ndim, 2]))
      call param%new(fade_enabled=.false.)
      call init_cfc_param(lsf, mol, [3.0_wp, 3.0_wp], 1.0e-11_wp, param)
      lsf%param%c = 50.0_wp
      call lsf%prepare([7.6_wp, 0.0_wp, 0.0_wp], lsf_err)
      call expect_error(error, lsf_err, "CFC parameters changed after update")
   end subroutine test_bad_cfc_stale_reach_off

   !> Tiling invariance, list, reach and band formula of the CFC fade
   !>
   !> - MB16-43 01, biting band `a = 5`, `b = 2` e-folds, production threshold;
   !>   P from [[cfc_band_points]] (an atom inside its band, every atom >= 0.05
   !>   Bohr off both edges)
   !> - bit for bit: P at column 1 of a strip jittered by <= 1e-3 Bohr against P
   !>   at columns 7 and `width` of other such strips, shared list ascending and
   !>   reversed; at [[BATCH_REL]]: P alone and a strip spread over 1 Bohr
   !> - the list is honoured: the list without P's in-band atom moves the value
   !> - reach = band end: `screening_offset(R) = R (a + b)/|a2|` for two radii
   !> - one carbon in its band, order 0, against `S = -a1 sigma - log w`,
   !>   `w = 1 - S_4(t)` of `moist_math_smoothstep` over the band in sigma (the
   !>   module docs' atomic term); negative control: the unfaded `-a1 sigma`
   !>   misses it
   subroutine test_cfc_fade_tiling(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: wmax = maxval(strip_width_cfc)
      type(structure_type) :: mol, one
      real(wp), allocatable :: radii(:), pts(:, :), d(:)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf, single
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      type(moist_math_smoothstep_polynomial_type) :: fade
      real(wp) :: strip(ndim, wmax), p(ndim), start, width, r_test(2), expected, sigma, w
      real(wp) :: a0(wmax), a1(ndim, wmax), a2(6, wmax), a3(10, wmax), a4(15, wmax)
      real(wp) :: b0(wmax), b1(ndim, wmax), b2(6, wmax), b3(10, wmax), b4(15, wmax)
      integer, allocatable :: full(:), rev(:), cut(:)
      integer :: order, width_o, icase, col, j, inband, st(1)
      integer :: state
      character(len=64) :: where
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_start=5.0_wp, fade_width=2.0_wp)
      call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
      start = 5.0_wp/abs(lsf%param%a2)
      width = 2.0_wp/abs(lsf%param%a2)
      call cfc_band_points(mol, radii, start, width, 0.05_wp, pts)
      call check(error, size(pts, 2) > 0, more="no point with an atom inside its band")
      if (allocated(error)) return
      p = pts(:, 1)
      d = band_distances(mol, radii, p)
      inband = findloc(d > radii(1:mol%nat)*start .and. d < radii(1:mol%nat)*(start + width), &
                       .true., dim=1)

      full = [(j, j=1, lsf%ncenters)]
      rev = full(lsf%ncenters:1:-1)
      cut = pack(full, full /= lsf%orig_to_sorted(inband))
      state = 13

      do order = 0, 4
         width_o = strip_width_cfc(order)
         call jittered_strip(p, 1.0e-3_wp, width_o, 1, state, strip)
         call lsf%evaluate_batch(strip(:, 1:width_o), order, a0, a1, a2, a3, a4, &
                                 candidate_indices=full, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         ! icase: 1, 2 bit for bit at columns 7 / width; 3 alone; 4 spread strip
         do icase = 1, 4
            select case (icase)
            case (1)
               col = min(7, width_o)
               call jittered_strip(p, 1.0e-3_wp, width_o, col, state, strip)
               call lsf%evaluate_batch(strip(:, 1:width_o), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=rev, error=lsf_err)
            case (2)
               col = width_o
               call jittered_strip(p, 1.0e-3_wp, width_o, col, state, strip)
               call lsf%evaluate_batch(strip(:, 1:width_o), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=full, error=lsf_err)
            case (3)
               col = 1
               call lsf%evaluate_batch(reshape(p, [ndim, 1]), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=full, error=lsf_err)
            case default
               col = min(7, width_o)
               call jittered_strip(p, 1.0_wp, width_o, col, state, strip)
               call lsf%evaluate_batch(strip(:, 1:width_o), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=full, error=lsf_err)
            end select
            if (allocated(lsf_err)) then
               call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
               return
            end if
            write (where, "(a,i0,a,i0,a,i0)") "order ", order, ", case ", icase, ", column ", col
            if (icase <= 2) then
               if (.not. same_column(order, 1, a0, a1, a2, a3, a4, col, b0, b1, b2, b3, b4)) then
                  call test_failed(error, "tile changes the bits of P, "//trim(where))
                  return
               end if
            else
               call check_column_close(error, order, a0, a1, a2, a3, a4, col, b0, b1, b2, b3, &
                                       b4, trim(where))
               if (allocated(error)) return
            end if
         end do

         call jittered_strip(p, 1.0e-3_wp, width_o, 1, state, strip)
         call lsf%evaluate_batch(strip(:, 1:width_o), order, b0, b1, b2, b3, b4, &
                                 candidate_indices=cut, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch on the cut list failed: "//lsf_err%message)
            return
         end if
         if (b0(1) == a0(1)) then
            call test_failed(error, "dropping an atom inside its band must change the value")
            return
         end if
      end do

      ! Reach = band end, per atom in units of its radius
      r_test = [2.3_wp, 3.4_wp]
      do j = 1, size(r_test)
         expected = r_test(j)*(5.0_wp + 2.0_wp)/abs(lsf%param%a2)
         call check(error, abs(lsf%screening_offset(r_test(j)) - expected) <= &
                    4.0_wp*epsilon(1.0_wp)*expected, &
                    more="fade on: screening_offset must be the band end R (a + b)/|a2|")
         if (allocated(error)) return
      end do

      ! One atom inside its band against the closed form
      call new(one, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]))
      call init_cfc_param(single, one, [3.0_wp], 1.0e-11_wp, param)
      sigma = start + 0.4_wp*width
      call new_smoothstep_polynomial(fade, start + width, start, order=4)
      call fade%eval(sigma, w)
      expected = -single%param%a1*sigma - log(w)
      call single%evaluate_batch(reshape([3.0_wp*(1.0_wp + sigma), 0.0_wp, 0.0_wp], [ndim, 1]), &
                                 0, b0, error=lsf_err)
      call batch_status(single, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "single-atom evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, st(1), 0, more="the single in-band atom must be reached")
      if (allocated(error)) return
      call check(error, abs(b0(1) - expected) <= 1.0e-12_wp*abs(expected), &
                 more="one in-band atom: S = -a1 sigma - log(1 - S_4(t))")
      if (allocated(error)) return
      call check(error, abs(-single%param%a1*sigma - expected) > 1.0e-6_wp*abs(expected), &
                 more="negative control: the band must bite at this sigma")
   end subroutine test_cfc_fade_tiling

   !* ================================================================================= *!
   !*                                  SvdW fade band                                   *!
   !* ================================================================================= *!

   !> SvdW instance from an explicit parameter set, `max_deriv = 4`
   !>
   !> @param[out] lsf       SvdW instance after `update`
   !> @param[in]  mol       Structure
   !> @param[in]  radii     Per-atom radii
   !> @param[in]  threshold Screening threshold, 0 disables screening
   !> @param[in]  param     Shape and fade parameters
   subroutine init_svdw_param(lsf, mol, radii, threshold, param)
      !> SvdW instance
      type(moist_cavity_drop_lsf_svdw_type), intent(out) :: lsf
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Screening threshold
      real(wp), intent(in) :: threshold
      !> Shape and fade parameters
      type(moist_cavity_drop_lsf_svdw_param_type), intent(in) :: param

      lsf%screening_threshold = threshold
      call lsf%new(param)
      call lsf%update(mol, radii)
      call lsf%set_max_deriv(4)
   end subroutine init_svdw_param

   !> Sample points off both band edges with one atom inside the band
   !>
   !> - candidates `R_A + R_a + s` along 14 directions per atom, `s` in 0.4..2.5 Bohr
   !> - kept when `|d_B - X_on| > margin` and `|d_B - X_on - Delta| > margin` for
   !>   every atom B and some atom lies inside `(X_on + margin, X_on + Delta - margin)`
   !>
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  x_on   Fade onset in `d` (Bohr)
   !> @param[in]  delta  Fade band width in `d` (Bohr)
   !> @param[in]  margin Distance kept from both edges (Bohr)
   !> @param[out] pts    Accepted points [ndim, npts]
   subroutine fade_band_points(mol, radii, x_on, delta, margin, pts)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Fade onset
      real(wp), intent(in) :: x_on
      !> Fade band width
      real(wp), intent(in) :: delta
      !> Distance kept from both edges
      real(wp), intent(in) :: margin
      !> Accepted points
      real(wp), allocatable, intent(out) :: pts(:, :)

      real(wp), parameter :: shells(5) = [0.4_wp, 0.8_wp, 1.2_wp, 1.8_wp, 2.5_wp]
      real(wp) :: dirs(ndim, 14), buf(ndim, 14*5*mol%nat), p(ndim), d(mol%nat)
      integer :: i, j, k, n, ia

      dirs = 0.0_wp
      do j = 1, ndim
         dirs(j, 2*j - 1) = 1.0_wp
         dirs(j, 2*j) = -1.0_wp
      end do
      k = 6
      do i = -1, 1, 2
         do j = -1, 1, 2
            do n = -1, 1, 2
               k = k + 1
               dirs(:, k) = real([i, j, n], wp)/sqrt(3.0_wp)
            end do
         end do
      end do
      n = 0
      do ia = 1, mol%nat
         do k = 1, size(shells)
            do j = 1, size(dirs, 2)
               p = mol%xyz(:, ia) + (radii(ia) + shells(k))*dirs(:, j)
               do i = 1, mol%nat
                  d(i) = norm2(p - mol%xyz(:, i)) - radii(i)
               end do
               if (any(abs(d - x_on) <= margin .or. abs(d - x_on - delta) <= margin)) cycle
               if (.not. any(d > x_on + margin .and. d < x_on + delta - margin)) cycle
               n = n + 1
               buf(:, n) = p
            end do
         end do
      end do
      pts = buf(:, 1:n)
   end subroutine fade_band_points

   !> Band-relative distance `d_B = ||p - R_B|| - R_b` of every atom
   !>
   !> @param[in] mol   Structure
   !> @param[in] radii Per-atom radii
   !> @param[in] p     Point
   !> @returns   d     Distances [mol%nat]
   pure function band_distances(mol, radii, p) result(d)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Point
      real(wp), intent(in) :: p(ndim)
      !> Distances
      real(wp) :: d(mol%nat)

      integer :: i

      do i = 1, mol%nat
         d(i) = norm2(p - mol%xyz(:, i)) - radii(i)
      end do
   end function band_distances

   !> Bitwise equality of column `i` of one batch and column `j` of another
   !>
   !> @param[in] order Batch order
   !> @param[in] i     Column of the first set
   !> @param[in] a0    First values
   !> @param[in] a1    First gradients
   !> @param[in] a2    First packed Hessians
   !> @param[in] a3    First packed third derivatives
   !> @param[in] a4    First packed fourth derivatives
   !> @param[in] j     Column of the second set
   !> @param[in] b0    Second values
   !> @param[in] b1    Second gradients
   !> @param[in] b2    Second packed Hessians
   !> @param[in] b3    Second packed third derivatives
   !> @param[in] b4    Second packed fourth derivatives
   !> @returns   same  Every output through `order` equal
   pure logical function same_column(order, i, a0, a1, a2, a3, a4, j, b0, b1, b2, b3, b4) &
      result(same)
      !> Batch order
      integer, intent(in) :: order
      !> Column of the first set
      integer, intent(in) :: i
      !> First values
      real(wp), intent(in) :: a0(:)
      !> First gradients
      real(wp), intent(in) :: a1(:, :)
      !> First packed Hessians
      real(wp), intent(in) :: a2(:, :)
      !> First packed third derivatives
      real(wp), intent(in) :: a3(:, :)
      !> First packed fourth derivatives
      real(wp), intent(in) :: a4(:, :)
      !> Column of the second set
      integer, intent(in) :: j
      !> Second values
      real(wp), intent(in) :: b0(:)
      !> Second gradients
      real(wp), intent(in) :: b1(:, :)
      !> Second packed Hessians
      real(wp), intent(in) :: b2(:, :)
      !> Second packed third derivatives
      real(wp), intent(in) :: b3(:, :)
      !> Second packed fourth derivatives
      real(wp), intent(in) :: b4(:, :)

      same = a0(i) == b0(j)
      if (order >= 1) same = same .and. all(a1(:, i) == b1(:, j))
      if (order >= 2) same = same .and. all(a2(:, i) == b2(:, j))
      if (order >= 3) same = same .and. all(a3(:, i) == b3(:, j))
      if (order >= 4) same = same .and. all(a4(:, i) == b4(:, j))
   end function same_column

   !> Tiling invariance of a point's batch outputs under a biting fade
   !>
   !> - MB16-43 01, band `a = 4`, `b = 2` e-folds, production threshold; point P
   !>   with an atom inside the band, every atom >= 0.05 Bohr off both edges
   !> - reference: P at column 1 of a full strip of co-tenants jittered by
   !>   <= 1e-3 Bohr, so the strip box classifies every atom (plateau, band,
   !>   beyond) as P alone does
   !> - bit for bit: P at columns 7 and `width` of a strip with other co-tenants,
   !>   shared list ascending and reversed (the override sorts it back)
   !> - at [[BATCH_REL]] only: P alone (`np = 1`, other loop code under FP
   !>   contraction) and co-tenants spread over 1 Bohr (atoms move
   !>   between the plain and the faded body)
   !> - negative control: dropping one atom inside P's band moves the value
   subroutine test_svdw_fade_tiling(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: wmax = maxval(strip_width_svdw)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), pts(:, :), d(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: strip(ndim, wmax), p(ndim), x_on, delta
      real(wp) :: a0(wmax), a1(ndim, wmax), a2(6, wmax), a3(10, wmax), a4(15, wmax)
      real(wp) :: b0(wmax), b1(ndim, wmax), b2(6, wmax), b3(10, wmax), b4(15, wmax)
      integer, allocatable :: full(:), rev(:), cut(:)
      integer :: order, width, icase, col, j, inband
      integer :: state
      character(len=64) :: where
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_start=4.0_wp, fade_width=2.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      x_on = 3.0_wp*4.0_wp/lsf%param%blend_k
      delta = 3.0_wp*2.0_wp/lsf%param%blend_k
      call fade_band_points(mol, radii, x_on, delta, 0.05_wp, pts)
      call check(error, size(pts, 2) > 0, more="no point with an atom inside the band")
      if (allocated(error)) return
      p = pts(:, 1)
      d = band_distances(mol, radii, p)
      inband = findloc(d > x_on .and. d < x_on + delta, .true., dim=1)

      full = [(j, j=1, lsf%ncenters)]
      rev = full(lsf%ncenters:1:-1)
      cut = pack(full, full /= lsf%orig_to_sorted(inband))
      state = 11

      do order = 0, 4
         width = strip_width_svdw(order)
         call jittered_strip(p, 1.0e-3_wp, width, 1, state, strip)
         call lsf%evaluate_batch(strip(:, 1:width), order, a0, a1, a2, a3, a4, &
                                 candidate_indices=full, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         ! icase: 1, 2 bit for bit at columns 7 / width; 3 alone; 4 spread strip
         do icase = 1, 4
            select case (icase)
            case (1)
               col = 7
               call jittered_strip(p, 1.0e-3_wp, width, col, state, strip)
               call lsf%evaluate_batch(strip(:, 1:width), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=rev, error=lsf_err)
            case (2)
               col = width
               call jittered_strip(p, 1.0e-3_wp, width, col, state, strip)
               call lsf%evaluate_batch(strip(:, 1:width), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=full, error=lsf_err)
            case (3)
               col = 1
               call lsf%evaluate_batch(reshape(p, [ndim, 1]), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=full, error=lsf_err)
            case default
               col = 7
               call jittered_strip(p, 1.0_wp, width, col, state, strip)
               call lsf%evaluate_batch(strip(:, 1:width), order, b0, b1, b2, b3, b4, &
                                       candidate_indices=full, error=lsf_err)
            end select
            if (allocated(lsf_err)) then
               call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
               return
            end if
            write (where, "(a,i0,a,i0,a,i0)") "order ", order, ", case ", icase, ", column ", col
            if (icase <= 2) then
               if (.not. same_column(order, 1, a0, a1, a2, a3, a4, col, b0, b1, b2, b3, b4)) then
                  call test_failed(error, "tile changes the bits of P, "//trim(where))
                  return
               end if
            else
               call check_column_close(error, order, a0, a1, a2, a3, a4, col, b0, b1, b2, b3, &
                                       b4, trim(where))
               if (allocated(error)) return
            end if
         end do

         call jittered_strip(p, 1.0e-3_wp, width, 1, state, strip)
         call lsf%evaluate_batch(strip(:, 1:width), order, b0, b1, b2, b3, b4, &
                                 candidate_indices=cut, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch on the cut list failed: "//lsf_err%message)
            return
         end if
         if (b0(1) == a0(1)) then
            call test_failed(error, "dropping an atom inside the band must change the value")
            return
         end if
      end do
   end subroutine test_svdw_fade_tiling

   !> Strip of `width` points jittered around `p`, `p` itself at column `col`
   !>
   !> @param[in]    p      Centre point
   !> @param[in]    spread Edge of the jitter cube (Bohr)
   !> @param[in]    width  Strip size
   !> @param[in]    col    Column holding `p` exactly
   !> @param[inout] state  LCG state
   !> @param[out]   strip  Strip points [ndim, >= width]
   subroutine jittered_strip(p, spread, width, col, state, strip)
      !> Centre point
      real(wp), intent(in) :: p(ndim)
      !> Edge of the jitter cube
      real(wp), intent(in) :: spread
      !> Strip size
      integer, intent(in) :: width
      !> Column holding `p`
      integer, intent(in) :: col
      !> LCG state
      integer, intent(inout) :: state
      !> Strip points
      real(wp), intent(out) :: strip(:, :)

      real(wp) :: u(ndim)
      integer :: j

      do j = 1, width
         call det_rand(u, state)
         strip(:, j) = p + spread*(u - 0.5_wp)
      end do
      strip(:, col) = p
   end subroutine jittered_strip

   !> One batch column against a reference column at [[BATCH_REL]] / [[BATCH_REL4]]
   !>
   !> @param[out] error Test failure
   !> @param[in]  order Batch order
   !> @param[in]  a0    Reference value [1]
   !> @param[in]  a1    Reference gradient [ndim, 1]
   !> @param[in]  a2    Reference packed Hessian [6, 1]
   !> @param[in]  a3    Reference packed third derivative [10, 1]
   !> @param[in]  a4    Reference packed fourth derivative [15, 1]
   !> @param[in]  j     Column of the batch
   !> @param[in]  b0    Batch values
   !> @param[in]  b1    Batch gradients
   !> @param[in]  b2    Batch packed Hessians
   !> @param[in]  b3    Batch packed third derivatives
   !> @param[in]  b4    Batch packed fourth derivatives
   !> @param[in]  what  Label for the message
   subroutine check_column_close(error, order, a0, a1, a2, a3, a4, j, b0, b1, b2, b3, b4, what)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Batch order
      integer, intent(in) :: order
      !> Reference value
      real(wp), intent(in) :: a0(:)
      !> Reference gradient
      real(wp), intent(in) :: a1(:, :)
      !> Reference packed Hessian
      real(wp), intent(in) :: a2(:, :)
      !> Reference packed third derivative
      real(wp), intent(in) :: a3(:, :)
      !> Reference packed fourth derivative
      real(wp), intent(in) :: a4(:, :)
      !> Column of the batch
      integer, intent(in) :: j
      !> Batch values
      real(wp), intent(in) :: b0(:)
      !> Batch gradients
      real(wp), intent(in) :: b1(:, :)
      !> Batch packed Hessians
      real(wp), intent(in) :: b2(:, :)
      !> Batch packed third derivatives
      real(wp), intent(in) :: b3(:, :)
      !> Batch packed fourth derivatives
      real(wp), intent(in) :: b4(:, :)
      !> Label for the message
      character(len=*), intent(in) :: what

      call batch_close(error, [b0(j)], [a0(1)], "value, "//what)
      if (allocated(error)) return
      if (order >= 1) call batch_close(error, b1(:, j), a1(:, 1), "gradient, "//what)
      if (allocated(error)) return
      if (order >= 2) call batch_close(error, b2(:, j), a2(:, 1), "Hessian, "//what)
      if (allocated(error)) return
      if (order >= 3) call batch_close(error, b3(:, j), a3(:, 1), "third, "//what)
      if (allocated(error)) return
      if (order >= 4) call batch_close(error, b4(:, j), a4(:, 1), "fourth, "//what, &
                                       rel=BATCH_REL4)
   end subroutine check_column_close

   !> A band beyond every atom reproduces the untruncated model bit for bit
   !>
   !> - reference: fade off, screening off (threshold 0)
   !> - band `a = 1000` e-folds (X_on ~ 545 Bohr), asserted beyond every point
   !>   and atom; batch at orders 0..4 and the per-point accessors, against fade off
   !> - negative control: a band inside the molecule (`a = 4`, `b = 2`) moves the
   !>   value of a reached column by more than 1e-6 relative
   subroutine test_svdw_fade_limit(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 37
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: off, far, near
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim), dmax
      real(wp) :: a0(n), a1(ndim, n), a2(6, n), a3(10, n), a4(15, n)
      real(wp) :: b0(n), b1(ndim, n), b2(6, n), b3(10, n), b4(15, n)
      real(wp) :: f0, f1(ndim), f2(6), f3(10), f4(15)
      real(wp) :: g0, g1(ndim), g2(6), g3(10), g4(15)
      integer :: sa(n), sb(n), i, order
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.false.)
      call init_svdw_param(off, mol, radii, 0.0_wp, param)
      call param%new(fade_enabled=.true., fade_start=1000.0_wp, fade_width=7.0_wp)
      call init_svdw_param(far, mol, radii, 0.0_wp, param)
      call param%new(fade_start=4.0_wp, fade_width=2.0_wp)
      call init_svdw_param(near, mol, radii, 0.0_wp, param)

      state = 31
      lo = off%center_lower - 3.0_wp
      hi = off%center_upper + 3.0_wp
      dmax = 0.0_wp
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
         dmax = max(dmax, maxval(band_distances(mol, radii, points(:, i))))
      end do
      call check(error, 3.0_wp*1000.0_wp/far%param%blend_k > dmax, &
                 more="the far band must start beyond every atom")
      if (allocated(error)) return

      do order = 0, 4
         call off%evaluate_batch(points, order, a0, a1, a2, a3, a4, error=lsf_err)
         call batch_status(off, sa)
         if (.not. allocated(lsf_err)) then
            call far%evaluate_batch(points, order, b0, b1, b2, b3, b4, &
                                    error=lsf_err)
            call batch_status(far, sb)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "evaluate_batch failed: "//lsf_err%message)
            return
         end if
         if (any(sa /= 0) .or. any(sb /= 0)) then
            call test_failed(error, "fade limit: a column near the molecule is unreached")
            return
         end if
         do i = 1, n
            if (.not. same_column(order, i, a0, a1, a2, a3, a4, i, b0, b1, b2, b3, b4)) then
               call test_failed(error, "band beyond every atom changes the batch bits")
               return
            end if
         end do
      end do

      do i = 1, n
         call off%prepare(points(:, i), lsf_err)
         if (.not. allocated(lsf_err)) call far%prepare(points(:, i), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "prepare failed: "//lsf_err%message)
            return
         end if
         call off%f3_rrr(lsf0=f0, lsf1_r=f1, lsf2_rr=f2, lsf3_rrr=f3)
         call off%f4_rrrr(f4)
         call far%f3_rrr(lsf0=g0, lsf1_r=g1, lsf2_rr=g2, lsf3_rrr=g3)
         call far%f4_rrrr(g4)
         if (f0 /= g0 .or. any(f1 /= g1) .or. any(f2 /= g2) .or. any(f3 /= g3) .or. &
             any(f4 /= g4)) then
            call test_failed(error, "band beyond every atom changes the per-point bits")
            return
         end if
      end do

      call near%evaluate_batch(points, 0, b0, error=lsf_err)
      call batch_status(near, sb)
      if (allocated(lsf_err)) then
         call test_failed(error, "evaluate_batch with the near band failed: "//lsf_err%message)
         return
      end if
      call check(error, any(sb == 0 .and. abs(b0 - a0) > 1.0e-6_wp*max(1.0_wp, abs(a0))), &
                 more="a band inside the molecule must change the value")
   end subroutine test_svdw_fade_limit

   !> Fade on: the atom reach is the band end at every screening threshold
   !>
   !> - thresholds 1e-11, 1e-9, 0 with the default band: `screening_offset`
   !>   equals `3(a + b)/k` (to rounding of the two divisions), batch and
   !>   `prepare` accepted, every value bit for bit the same
   !> - an atom beyond the band end is screened (exact zero of the fade there)
   !>   while one inside the band is kept
   !> - negative control: fade off at 1e-9 gives the threshold offset
   !>   `-3 ln(1e-9)/k`, shorter than the band end
   subroutine test_svdw_fade_reach(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: thresholds(3) = [1.0e-11_wp, 1.0e-9_wp, 0.0_wp]
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: points(ndim, 2), b0(2), ref(2), band_end, f0, f0_ref
      integer :: icase
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      points(:, 1) = mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp]
      points(:, 2) = mol%xyz(:, 2) + [0.0_wp, 2.0_wp, 0.0_wp]
      do icase = 1, size(thresholds)
         call param%new(fade_enabled=.true.)
         call init_svdw_param(lsf, mol, radii, thresholds(icase), param)
         band_end = 3.0_wp*(lsf%param%fade_start + lsf%param%fade_width)/lsf%param%blend_k
         call check(error, lsf%screening_offset(radii(1)), band_end, thr=8.0_wp*epsilon(1.0_wp), &
                    more="fade on: the reach must be the band end")
         if (allocated(error)) return
         call lsf%evaluate_batch(points, 0, b0, error=lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "batch refused with the fade on: "//lsf_err%message)
            return
         end if
         call lsf%prepare(points(:, 1), lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "prepare refused with the fade on: "//lsf_err%message)
            return
         end if
         call lsf%f0(f0)
         if (icase == 1) then
            ref = b0
            f0_ref = f0
         else if (any(b0 /= ref) .or. f0 /= f0_ref) then
            call test_failed(error, "the screening threshold changes the faded model")
            return
         end if
      end do

      ! One atom just inside, then just beyond the band end of atom 1
      call lsf%prepare(mol%xyz(:, 1) + [radii(1) + band_end - 1.0e-3_wp, 0.0_wp, 0.0_wp], &
                       lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "prepare failed: "//lsf_err%message)
         return
      end if
      call check(error, any([(lsf%active_atom(icase), icase=1, lsf%active_count())] == 1), &
                 more="an atom inside its band must be kept")
      if (allocated(error)) return
      call lsf%prepare(mol%xyz(:, 1) - [radii(1) + band_end + 1.0e-3_wp, 0.0_wp, 0.0_wp], &
                       lsf_err)
      if (allocated(lsf_err)) then
         call test_failed(error, "prepare failed: "//lsf_err%message)
         return
      end if
      if (lsf%active_count() > 0) then
         call check(error, .not. any([(lsf%active_atom(icase), icase=1, lsf%active_count())] == 1), &
                    more="an atom beyond its band end must be screened")
         if (allocated(error)) return
      end if

      call param%new(fade_enabled=.false.)
      call init_svdw_param(lsf, mol, radii, 1.0e-9_wp, param)
      call check(error, lsf%screening_offset(radii(1)), -3.0_wp*log(1.0e-9_wp)/lsf%param%blend_k, &
                 thr=8.0_wp*epsilon(1.0_wp), more="fade off: the threshold offset")
      if (allocated(error)) return
      call check(error, lsf%screening_offset(radii(1)) < band_end, &
                 more="the control offset must differ from the band end")
   end subroutine test_svdw_fade_reach

   !> Fail the test only when `err` carries the expected message
   !>
   !> For `should_fail` tests: a clean return or another error reports the test as failed
   !>
   !> @param[out] error    Test failure, set only on the expected error
   !> @param[in]  err      Library error, possibly unallocated
   !> @param[in]  expected Substring the message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(mctc_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> Per-column flag of the latest batch on `lsf`, read right after the call
   !>
   !> - 0 reached, 1 unreached ([[moist_cavity_drop_lsf_type:unreached]])
   !> - -1 past the latest batch width (a refused call writes no mask)
   !>
   !> @param[in]  lsf LSF right after `evaluate_batch` / `vjp_f1_rA_batch`
   !> @param[out] st  Per-column flag
   subroutine batch_status(lsf, st)
      !> LSF right after the batch call
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Per-column flag
      integer, intent(out) :: st(:)

      !> Column
      integer :: i

      do i = 1, size(st)
         if (i > lsf%batch_width) then
            st(i) = -1
         else
            st(i) = merge(1, 0, lsf%unreached(i))
         end if
      end do
   end subroutine batch_status

   !> Non-positive fade width refused by `prepare` (parameter validation)
   subroutine test_bad_fade_width(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_width=0.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%prepare(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], lsf_err)
      call expect_error(error, lsf_err, "must be positive")
   end subroutine test_bad_fade_width

   !> Fade band changed after `update` refused by `prepare` (stale reach and certificates)
   subroutine test_bad_param_changed(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_start=4.0_wp, fade_width=2.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      lsf%param%fade_start = 21.0_wp
      lsf%param%fade_width = 3.0_wp
      call lsf%prepare(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], lsf_err)
      call expect_error(error, lsf_err, "SvdW fade parameters changed after update")
   end subroutine test_bad_param_changed

   !> SvdW reach refreshed under another band, band restored: `evaluate_batch` refuses
   !>
   !> - `fade_start` 4 -> 2, `set_screening_threshold` rebuilds the reach from
   !>   it, `fade_start` back to 4: the recorded parameters match `update`'s,
   !>   the cached reach does not
   subroutine test_bad_svdw_stale_reach(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: b0(1)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_start=4.0_wp, fade_width=2.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      lsf%param%fade_start = 2.0_wp
      call lsf%set_screening_threshold(1.0e-10_wp)
      lsf%param%fade_start = 4.0_wp
      call lsf%evaluate_batch(reshape(mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp], [ndim, 1]), &
                              0, b0, error=lsf_err)
      call expect_error(error, lsf_err, "SvdW fade parameters changed after update")
   end subroutine test_bad_svdw_stale_reach

   !> Spatial derivatives of orders 1..3 against FD with every stencil inside one band piece
   !>
   !> - band `a = 4`, `b = 2` e-folds, production threshold, MB16-43 01; points
   !>   from [[fade_band_points]] with every atom >= 0.05 Bohr off both edges
   !>   (stencil half-width `2 STEP_SIZE`), one atom inside the band
   !> - per-point path (`prepare` + accessors), 4-point central FD, suite thresholds
   !> - negative control: the unfaded gradient misses the faded FD by more than
   !>   100 times the threshold, so the check sees the fade terms
   subroutine test_svdw_fade_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: npts_max = 6
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), pts(:, :)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf, off
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: point(ndim), shifted(ndim), x_on, delta, h, miss
      real(wp) :: f(-2:2), g(ndim, -2:2), hs(ndim, ndim, -2:2), packed_h(6), packed_t(10)
      real(wp) :: g1(ndim), h2(ndim, ndim), d3(ndim, ndim, ndim), goff(ndim)
      real(wp) :: n1(ndim), n2(ndim, ndim), n3(ndim, ndim, ndim)
      integer :: ipt, axis, step, i, j
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.false.)
      call init_svdw_param(off, mol, radii, 1.0e-11_wp, param)
      call param%new(fade_enabled=.true., fade_start=4.0_wp, fade_width=2.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%set_max_deriv(3)
      x_on = 3.0_wp*4.0_wp/lsf%param%blend_k
      delta = 3.0_wp*2.0_wp/lsf%param%blend_k
      call fade_band_points(mol, radii, x_on, delta, 0.05_wp, pts)
      call check(error, size(pts, 2) >= npts_max, more="too few points with an atom in the band")
      if (allocated(error)) return

      h = STEP_SIZE
      miss = 0.0_wp
      do ipt = 1, npts_max
         point = pts(:, 1 + (ipt - 1)*(size(pts, 2)/npts_max))
         call lsf%prepare(point, lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "prepare failed: "//lsf_err%message)
            return
         end if
         call lsf%f3_rrr(lsf1_r=g1, lsf2_rr=packed_h, lsf3_rrr=packed_t)
         call unpack_sym_matrix(packed_h, h2, ndim)
         call unpack_sym_tensor3(packed_t, d3, ndim)
         do axis = 1, ndim
            do step = -2, 2
               if (step == 0) cycle
               shifted = point
               shifted(axis) = point(axis) + real(step, wp)*h
               call lsf%prepare(shifted, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "displaced prepare failed: "//lsf_err%message)
                  return
               end if
               call lsf%f3_rrr(lsf0=f(step), lsf1_r=g(:, step), lsf2_rr=packed_h, &
                               lsf3_rrr=packed_t)
               call unpack_sym_matrix(packed_h, hs(:, :, step), ndim)
            end do
            call fd4_scalar(f(2), f(1), f(-1), f(-2), h, n1(axis), error)
            if (allocated(error)) return
            do i = 1, ndim
               call fd4_scalar(g(i, 2), g(i, 1), g(i, -1), g(i, -2), h, n2(i, axis), error)
               if (allocated(error)) return
               do j = 1, ndim
                  call fd4_scalar(hs(i, j, 2), hs(i, j, 1), hs(i, j, -1), hs(i, j, -2), h, &
                                  n3(i, j, axis), error)
                  if (allocated(error)) return
               end do
            end do
         end do
         do i = 1, ndim
            call check(error, g1(i), n1(i), thr_abs=fd_abs_thr("svdw"), &
                       thr_rel=fd_rel_thr("svdw"))
            if (allocated(error)) return
            do j = 1, ndim
               call check(error, h2(i, j), n2(i, j), thr_abs=fd_abs_thr("svdw"), &
                          thr_rel=fd_rel_thr("svdw"))
               if (allocated(error)) return
               do step = 1, ndim
                  call check(error, d3(i, j, step), n3(i, j, step), thr_abs=fd_abs_thr("svdw"), &
                             thr_rel=fd_rel_thr("svdw"))
                  if (allocated(error)) return
               end do
            end do
         end do

         call off%prepare(point, lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "unfaded prepare failed: "//lsf_err%message)
            return
         end if
         call off%f012_r(lsf1_r=goff)
         miss = max(miss, maxval(abs(goff - n1)))
      end do
      call check(error, miss > 100.0_wp*fd_abs_thr("svdw"), &
                 more="the band must change the gradient beyond the FD threshold")
   end subroutine test_svdw_fade_fd

   !> SvdW dispatch of [[run_fade_fd_atom_families]]
   subroutine test_svdw_fade_fd_atom_families(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_fade_fd_atom_families(error, "svdw")
   end subroutine test_svdw_fade_fd_atom_families

   !> CFC dispatch of [[run_fade_fd_atom_families]]
   subroutine test_cfc_fade_fd_atom_families(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_fade_fd_atom_families(error, "cfc")
   end subroutine test_cfc_fade_fd_atom_families

   !> Per-atom families with the fade on against FD inside one band piece
   !>
   !> - MB16-43 01, production threshold; every atom >= 0.05 Bohr off both band
   !>   edges, one inside its band; every stencil `+-3 STEP_SIZE` in one piece

   !>   - SvdW: band `a = 4`, `b = 2` e-folds in `d_A` ([[fade_band_points]])
   !>   - CFC: band `a = 5`, `b = 2` e-folds of `u_A`, per atom in sigma
   !>     ([[cfc_band_points]]); the radius rows check the faded radius rule
   !>     `dsigma/dR = -||d||/R**2`, the radius also scaling the band
   !> - one center or radius displaced per stencil, 6-point central FD of
   !>   `S`, `dS/dr`, `d2S/dr2`, `f1_rA`, `f1_rad` and the weighted jet
   !>   (4-point truncation of the CFC radius rows near a band end: 7e-9 at
   !>   `STEP_SIZE`, over `CFC_ABS_THR`)
   !>   `w0 S + w1 . dS/dr + w2 : d2S/dr2`:
   !>   - centers: `f1_rA`, `f2_r_rA`, `f3_rr_rA`, `f2_rArB`, `vjp_f1_rA`
   !>   - radii: `f1_rad`, `f2_r_rad`, `f3_rr_rad`, `f2_radrad`, `f2_rA_rad`,
   !>     `vjp_f1_rad`
   !>   - joint direction of [[joint_direction]]: `hvp_f1_rA` (with `vrad`),
   !>     `hvp_f1_rad`
   !> - results mapped to user space through `active_atom`; atoms beyond the band
   !>   end (no slot, every one >= 0.05 Bohr off it) carry exact zeros on both sides
   !> - thresholds: SvdW `ABS_THR` / `REL_THR`, CFC `CFC_ABS_THR` / `CFC_REL_THR`
   !> - negative control: the unfaded `f1_rA` and `f1_rad` miss the faded FD by
   !>   more than 100 times the threshold
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  "svdw" or "cfc"
   subroutine run_fade_fd_atom_families(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> LSF kind
      character(len=*), intent(in) :: kind

      integer, parameter :: npts_max = 3
      real(wp), parameter :: steps(6) = [3.0_wp, 2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp, -3.0_wp]
      real(wp), parameter :: w0 = 0.73_wp, w1(ndim) = [-0.41_wp, 0.88_wp, 0.17_wp]
      real(wp), parameter :: w2(ndim, ndim) = reshape([0.29_wp, 0.62_wp, -0.36_wp, &
                                                       -0.54_wp, 0.08_wp, 0.45_wp, &
                                                       0.13_wp, -0.77_wp, 0.91_wp], [ndim, ndim])
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), pts(:, :), centers(:, :), rad_local(:), cen_local(:, :)
      real(wp), allocatable :: v(:, :), vrad(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf, off
      type(moist_cavity_drop_lsf_svdw_param_type) :: svdw_param
      type(moist_cavity_drop_lsf_cfc_param_type) :: cfc_param
      ! analytic, user space
      real(wp), allocatable :: a1(:, :), a2(:, :, :), a3(:, :, :, :), ab(:, :, :, :)
      real(wp), allocatable :: r1(:), r2(:, :), r3(:, :, :), rr(:, :), rx(:, :, :)
      real(wp), allocatable :: vja(:, :), vjr(:), hva(:, :), hvr(:), o1(:, :), or1(:)
      ! stencil samples, user space
      real(wp) :: s0(6), sg(ndim, 6), sh(ndim, ndim, 6), sj(6)
      real(wp), allocatable :: sa(:, :, :), sr(:, :)
      real(wp) :: point(ndim), numeric, miss, thr_a, thr_r
      integer :: ipt, atom, other, dir, istep, i, j, nat
      logical :: joint

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      nat = mol%nat
      call init_lsf(off, mol, radii, 4, kind, screening_threshold=1.0e-11_wp, &
                    fade_enabled=.false.)
      if (kind == "svdw") then
         call init_lsf(lsf, mol, radii, 4, kind, screening_threshold=1.0e-11_wp, &
                       fade_enabled=.true., fade_start=4.0_wp, fade_width=2.0_wp)
         call fade_band_points(mol, radii, 3.0_wp*4.0_wp/svdw_param%blend_k, &
                               3.0_wp*2.0_wp/svdw_param%blend_k, 0.05_wp, pts)
         thr_a = ABS_THR
         thr_r = REL_THR
      else
         call init_lsf(lsf, mol, radii, 4, kind, screening_threshold=1.0e-11_wp, &
                       fade_enabled=.true., fade_start=5.0_wp, fade_width=2.0_wp)
         call cfc_band_points(mol, radii, 5.0_wp/abs(cfc_param%a2), &
                              2.0_wp/abs(cfc_param%a2), 0.05_wp, pts)
         thr_a = CFC_ABS_THR
         thr_r = CFC_REL_THR
      end if
      call check(error, size(pts, 2) >= npts_max, more="too few points with an atom in the band")
      if (allocated(error)) return

      allocate (centers(ndim, nat), cen_local(ndim, nat), rad_local(nat), v(ndim, nat), vrad(nat))
      allocate (a1(ndim, nat), a2(ndim, ndim, nat), a3(ndim, ndim, ndim, nat), &
                ab(ndim, nat, ndim, nat), r1(nat), r2(ndim, nat), r3(ndim, ndim, nat), &
                rr(nat, nat), rx(ndim, nat, nat), vja(ndim, nat), vjr(nat), hva(ndim, nat), &
                hvr(nat), o1(ndim, nat), or1(nat), sa(ndim, nat, 6), sr(nat, 6))
      centers = mol%xyz
      call joint_direction(nat, v, vrad)

      miss = 0.0_wp
      do ipt = 1, npts_max
         point = pts(:, 1 + (ipt - 1)*(size(pts, 2)/npts_max))
         ! negative control: the unfaded nuclear and radius gradients at this point
         call refresh_joint(off, mol, centers, radii)
         if (prepare_failed(off, point, error)) return
         call user_f1(off, nat, o1, or1, error)
         if (allocated(error)) return

         call refresh_joint(lsf, mol, centers, radii)
         call fade_family_analytic(lsf, nat, point, v, vrad, w0, w1, w2, a1, a2, a3, ab, &
                                   r1, r2, r3, rr, rx, vja, vjr, hva, hvr, error)
         if (allocated(error)) return

         ! dir 1..3: center axis, 4: radius, 5: joint direction (once, atom 1)
         do atom = 1, nat
            do dir = 1, 5
               joint = dir == 5
               if (joint .and. atom > 1) cycle
               do istep = 1, 6
                  cen_local = centers
                  rad_local = radii
                  if (dir <= ndim) then
                     cen_local(dir, atom) = cen_local(dir, atom) + steps(istep)*STEP_SIZE
                  else if (dir == 4) then
                     rad_local(atom) = rad_local(atom) + steps(istep)*STEP_SIZE
                  else
                     cen_local = cen_local + steps(istep)*STEP_SIZE*v
                     rad_local = rad_local + steps(istep)*STEP_SIZE*vrad
                  end if
                  call refresh_joint(lsf, mol, cen_local, rad_local)
                  call fade_family_sample(lsf, nat, point, w0, w1, w2, s0(istep), &
                                          sg(:, istep), sh(:, :, istep), sj(istep), &
                                          sa(:, :, istep), sr(:, istep), error)
                  if (allocated(error)) return
               end do
               if (joint) then
                  do other = 1, nat
                     call fd6_scalar(sr(other, 1), sr(other, 2), sr(other, 3), &
                                     sr(other, 4), sr(other, 5), sr(other, 6), &
                                     STEP_SIZE, numeric, error)
                     if (allocated(error)) return
                     call check(error, hvr(other), numeric, thr_abs=thr_a, thr_rel=thr_r, &
                                more="hvp_f1_rad")
                     if (allocated(error)) return
                     do i = 1, ndim
                        call fd6_scalar(sa(i, other, 1), sa(i, other, 2), sa(i, other, 3), &
                                        sa(i, other, 4), sa(i, other, 5), sa(i, other, 6), &
                                        STEP_SIZE, numeric, error)
                        if (allocated(error)) return
                        call check(error, hva(i, other), numeric, thr_abs=thr_a, &
                                   thr_rel=thr_r, more="hvp_f1_rA")
                        if (allocated(error)) return
                     end do
                  end do
                  cycle
               end if

               ! value, gradient, Hessian and weighted jet of this displacement
               call fd6_scalar(s0(1), s0(2), s0(3), &
                               s0(4), s0(5), s0(6), &
                               STEP_SIZE, numeric, error)
               if (allocated(error)) return
               if (dir <= ndim) then
                  call check(error, a1(dir, atom), numeric, thr_abs=thr_a, thr_rel=thr_r, &
                             more="f1_rA")
                  miss = max(miss, abs(o1(dir, atom) - numeric))
               else
                  call check(error, r1(atom), numeric, thr_abs=thr_a, thr_rel=thr_r, more="f1_rad")
                  miss = max(miss, abs(or1(atom) - numeric))
               end if
               if (allocated(error)) return
               call fd6_scalar(sj(1), sj(2), sj(3), &
                               sj(4), sj(5), sj(6), &
                               STEP_SIZE, numeric, error)
               if (allocated(error)) return
               if (dir <= ndim) then
                  call check(error, vja(dir, atom), numeric, thr_abs=thr_a, thr_rel=thr_r, &
                             more="vjp_f1_rA")
               else
                  call check(error, vjr(atom), numeric, thr_abs=thr_a, thr_rel=thr_r, &
                             more="vjp_f1_rad")
               end if
               if (allocated(error)) return
               do i = 1, ndim
                  call fd6_scalar(sg(i, 1), sg(i, 2), sg(i, 3), &
                                  sg(i, 4), sg(i, 5), sg(i, 6), &
                                  STEP_SIZE, numeric, error)
                  if (allocated(error)) return
                  if (dir <= ndim) then
                     call check(error, a2(i, dir, atom), numeric, thr_abs=thr_a, &
                                thr_rel=thr_r, more="f2_r_rA")
                  else
                     call check(error, r2(i, atom), numeric, thr_abs=thr_a, thr_rel=thr_r, &
                                more="f2_r_rad")
                  end if
                  if (allocated(error)) return
                  do j = 1, ndim
                     call fd6_scalar(sh(i, j, 1), sh(i, j, 2), sh(i, j, 3), &
                                     sh(i, j, 4), sh(i, j, 5), sh(i, j, 6), &
                                     STEP_SIZE, numeric, error)
                     if (allocated(error)) return
                     if (dir <= ndim) then
                        call check(error, a3(i, j, dir, atom), numeric, thr_abs=thr_a, &
                                   thr_rel=thr_r, more="f3_rr_rA")
                     else
                        call check(error, r3(i, j, atom), numeric, thr_abs=thr_a, &
                                   thr_rel=thr_r, more="f3_rr_rad")
                     end if
                     if (allocated(error)) return
                  end do
               end do

               ! second atom index: f2_rArB / f2_rA_rad from f1_rA, f2_radrad from f1_rad
               do other = 1, nat
                  do i = 1, ndim
                     call fd6_scalar(sa(i, other, 1), sa(i, other, 2), sa(i, other, 3), &
                                     sa(i, other, 4), sa(i, other, 5), sa(i, other, 6), &
                                     STEP_SIZE, numeric, error)
                     if (allocated(error)) return
                     if (dir <= ndim) then
                        call check(error, ab(i, other, dir, atom), numeric, thr_abs=thr_a, &
                                   thr_rel=thr_r, more="f2_rArB")
                     else
                        call check(error, rx(i, other, atom), numeric, thr_abs=thr_a, &
                                   thr_rel=thr_r, more="f2_rA_rad")
                     end if
                     if (allocated(error)) return
                  end do
                  if (dir == 4) then
                     call fd6_scalar(sr(other, 1), sr(other, 2), sr(other, 3), &
                                     sr(other, 4), sr(other, 5), sr(other, 6), &
                                     STEP_SIZE, numeric, error)
                     if (allocated(error)) return
                     call check(error, rr(other, atom), numeric, thr_abs=thr_a, &
                                thr_rel=thr_r, more="f2_radrad")
                     if (allocated(error)) return
                  end if
               end do
            end do
         end do
      end do
      call check(error, miss > 100.0_wp*thr_a, &
                 more="the band must change f1_rA / f1_rad beyond the FD threshold")
   end subroutine run_fade_fd_atom_families

   !> Faded analytic per-atom families at the bound geometry, user space
   !>
   !> @param[inout] lsf   Faded SvdW or CFC instance, prepared here
   !> @param[in]    nat   Atom count; atoms beyond the band end get zeros
   !> @param[in]    point Evaluation point
   !> @param[in]    v     Joint direction, nuclear half [ndim, nat]
   !> @param[in]    vrad  Joint direction, radius half [nat]
   !> @param[in]    w0    Jet weight of S
   !> @param[in]    w1    Jet weights of dS/dr
   !> @param[in]    w2    Jet weights of d2S/dr2
   !> @param[out]   a1    f1_rA
   !> @param[out]   a2    f2_r_rA
   !> @param[out]   a3    f3_rr_rA
   !> @param[out]   ab    f2_rArB
   !> @param[out]   r1    f1_rad
   !> @param[out]   r2    f2_r_rad
   !> @param[out]   r3    f3_rr_rad
   !> @param[out]   rr    f2_radrad
   !> @param[out]   rx    f2_rA_rad
   !> @param[out]   vja   vjp_f1_rA
   !> @param[out]   vjr   vjp_f1_rad
   !> @param[out]   hva   hvp_f1_rA with vrad
   !> @param[out]   hvr   hvp_f1_rad
   !> @param[out]   error Test failure
   subroutine fade_family_analytic(lsf, nat, point, v, vrad, w0, w1, w2, a1, a2, a3, ab, &
                                   r1, r2, r3, rr, rx, vja, vjr, hva, hvr, error)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      integer, intent(in) :: nat
      real(wp), intent(in) :: point(ndim), v(:, :), vrad(:), w0, w1(ndim), w2(ndim, ndim)
      real(wp), intent(out) :: a1(:, :), a2(:, :, :), a3(:, :, :, :), ab(:, :, :, :)
      real(wp), intent(out) :: r1(:), r2(:, :), r3(:, :, :), rr(:, :), rx(:, :, :)
      real(wp), intent(out) :: vja(:, :), vjr(:), hva(:, :), hvr(:)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: t1(ndim, nat), t2(ndim, ndim, nat), t3(ndim, ndim, ndim, nat)
      real(wp) :: tb(ndim, nat, ndim, nat), q1(nat), q2(ndim, nat), q3(ndim, ndim, nat)
      real(wp) :: qq(nat, nat), qx(ndim, nat, nat)
      real(wp) :: p3(6, ndim, nat), pq3(6, nat), w2_packed(6)
      integer :: map(nat), i, j, nact

      if (prepare_failed(lsf, point, error)) return
      call active_map_user(lsf, nat, map, nact)
      call lsf%f3_rr_rA(lsf1_rA=t1, lsf2_r_rA=t2, lsf3_rr_rA=p3)
      call unpack_pair_batch(p3, t3, ndim*nat)
      call lsf%f2_rArB(tb)
      call lsf%f3_rr_rad(lsf1_rad=q1, lsf2_r_rad=q2, lsf3_rr_rad=pq3)
      call unpack_pair_batch(pq3, q3, nat)
      call lsf%f2_radrad(qq)
      call lsf%f2_rA_rad(qx)
      a1 = 0.0_wp
      a2 = 0.0_wp
      a3 = 0.0_wp
      ab = 0.0_wp
      r1 = 0.0_wp
      r2 = 0.0_wp
      r3 = 0.0_wp
      rr = 0.0_wp
      rx = 0.0_wp
      vja = 0.0_wp
      vjr = 0.0_wp
      hva = 0.0_wp
      hvr = 0.0_wp
      do i = 1, nact
         a1(:, map(i)) = t1(:, i)
         a2(:, :, map(i)) = t2(:, :, i)
         a3(:, :, :, map(i)) = t3(:, :, :, i)
         r1(map(i)) = q1(i)
         r2(:, map(i)) = q2(:, i)
         r3(:, :, map(i)) = q3(:, :, i)
         do j = 1, nact
            ab(:, map(i), :, map(j)) = tb(:, i, :, j)
            rr(map(i), map(j)) = qq(i, j)
            rx(:, map(i), map(j)) = qx(:, i, j)
         end do
      end do
      ! The weights are not symmetric; the packed seed carries their symmetric part
      call pack_sym_part(w2, w2_packed, ndim)
      call lsf%vjp_f1_rA(w0, w1, w2_packed, t1)
      call lsf%vjp_f1_rad(w0, w1, w2_packed, q1)
      do i = 1, nact
         vja(:, map(i)) = t1(:, i)
         vjr(map(i)) = q1(i)
      end do
      ! the directions are indexed by user atom, as the accessors read them
      call lsf%hvp_f1_rA(v, t1, vrad)
      call lsf%hvp_f1_rad(v, vrad, q1)
      do i = 1, nact
         hva(:, map(i)) = t1(:, i)
         hvr(map(i)) = q1(i)
      end do
   end subroutine fade_family_analytic

   !> One FD stencil sample of the faded families, user space
   !>
   !> @param[inout] lsf   Faded SvdW or CFC instance at the displaced geometry
   !> @param[in]    nat   Atom count; atoms beyond the band end get zeros
   !> @param[in]    point Evaluation point
   !> @param[in]    w0    Jet weight of S
   !> @param[in]    w1    Jet weights of dS/dr
   !> @param[in]    w2    Jet weights of d2S/dr2
   !> @param[out]   s0    S
   !> @param[out]   sg    dS/dr
   !> @param[out]   sh    d2S/dr2
   !> @param[out]   sj    Weighted jet
   !> @param[out]   sa    f1_rA [ndim, nat]
   !> @param[out]   sr    f1_rad [nat]
   !> @param[out]   error Test failure
   subroutine fade_family_sample(lsf, nat, point, w0, w1, w2, s0, sg, sh, sj, sa, sr, error)
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      integer, intent(in) :: nat
      real(wp), intent(in) :: point(ndim), w0, w1(ndim), w2(ndim, ndim)
      real(wp), intent(out) :: s0, sg(ndim), sh(ndim, ndim), sj, sa(:, :), sr(:)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: t1(ndim, nat), q1(nat), sh_packed(6)

      if (prepare_failed(lsf, point, error)) return
      call lsf%f012_r(lsf0=s0, lsf1_r=sg, lsf2_rr=sh_packed)
      call unpack_sym_matrix(sh_packed, sh, ndim)
      sj = w0*s0 + sum(w1*sg) + sum(w2*sh)
      call user_f1(lsf, nat, t1, q1, error)
      if (allocated(error)) return
      sa = t1
      sr = q1
   end subroutine fade_family_sample

   !> `f1_rA` and `f1_rad` of a prepared SvdW or CFC instance, user space
   !>
   !> @param[in]  lsf   Prepared SvdW or CFC instance
   !> @param[in]  nat   Atom count
   !> @param[out] f1    f1_rA [ndim, nat], zero for atoms beyond the band end
   !> @param[out] r1    f1_rad [nat], zero for atoms beyond the band end
   !> @param[out] error Test failure (none set; kept for the call sites)
   subroutine user_f1(lsf, nat, f1, r1, error)
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      integer, intent(in) :: nat
      real(wp), intent(out) :: f1(:, :), r1(:)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: t1(ndim, nat), q1(nat), q3(6, nat), t3(6, ndim, nat)
      integer :: map(nat), i, nact

      call active_map_user(lsf, nat, map, nact)
      call lsf%f3_rr_rA(lsf1_rA=t1, lsf3_rr_rA=t3)
      call lsf%f3_rr_rad(lsf1_rad=q1, lsf3_rr_rad=q3)
      f1 = 0.0_wp
      r1 = 0.0_wp
      do i = 1, nact
         f1(:, map(i)) = t1(:, i)
         r1(map(i)) = q1(i)
      end do
   end subroutine user_f1

   !> Active slot -> user atom map
   !>
   !> - atoms beyond the fade band end have no slot; their user-space entries
   !>   are exact zeros, which the callers fill
   !>
   !> @param[in]  lsf   Prepared LSF
   !> @param[in]  nat   Atom count
   !> @param[out] map   User atom of active slot i, first `nact` valid
   !> @param[out] nact  Active count
   subroutine active_map_user(lsf, nat, map, nact)
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      integer, intent(in) :: nat
      integer, intent(out) :: map(:)
      integer, intent(out) :: nact

      integer :: i

      map = 0
      nact = min(lsf%active_count(), nat)
      do i = 1, nact
         map(i) = lsf%active_atom(i)
      end do
   end subroutine active_map_user

   !> `set_full_scan_below`: the cell list mode changes, the candidates do not
   !>
   !> - 1000-atom lattice, `candidates_in_box` at the default (binned, 1000 >= 200)
   !>   and after `set_full_scan_below(2000)` (one cell): same ascending ids
   !> - negative control: the box query prunes the lattice in both modes, so the
   !>   comparison is not between two full lists
   subroutine test_svdw_set_full_scan_below(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      integer, allocatable :: binned(:), single(:)
      real(wp) :: lo(ndim), hi(ndim)
      integer :: nb, ns

      call lattice_structure(mol, 10, 3.0_wp, 5)
      call get_test_radii(mol, radii)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      allocate (binned(lsf%ncenters), single(lsf%ncenters))
      lo = lsf%center_lower + 1.0_wp
      hi = lsf%center_lower + 4.0_wp
      call check(error, lsf%cell_full_scan_below, 200)
      if (allocated(error)) return
      call check(error, .not. lsf%cand_cells%full_scan, more="default lattice list must be binned")
      if (allocated(error)) return
      call lsf%candidates_in_box(lo, hi, binned, nb)
      call lsf%set_full_scan_below(2000)
      call check(error, lsf%cand_cells%full_scan, more="setter must collapse the list to one cell")
      if (allocated(error)) return
      call lsf%candidates_in_box(lo, hi, single, ns)
      call check(error, nb > 0 .and. nb < lsf%ncenters, more="box query must prune the lattice")
      if (allocated(error)) return
      call check(error, ns == nb, more="candidate count depends on the cell list mode")
      if (allocated(error)) return
      call check(error, all(single(1:ns) == binned(1:nb)), &
                 more="candidates depend on the cell list mode")
   end subroutine test_svdw_set_full_scan_below

   !* ================================================================================= *!
   !*                         Batched nuclear vjp (`vjp_f1_rA_batch`)                   *!
   !* ================================================================================= *!

   !> Fixed seeds of column `p`: value, gradient, a non-symmetric dense Hessian seed
   !> and the packed symmetric part of it
   !>
   !> - closed-form entries, no RNG; `dense(a, b) /= dense(b, a)`
   !> - `w2` is the upper triangle of `(dense + dense^T)/2`, off-diagonals not doubled
   !>
   !> @param[in]  p     Column
   !> @param[out] w0    Value weight
   !> @param[out] w1    Gradient weights
   !> @param[out] dense Non-symmetric Hessian seed
   !> @param[out] w2    Packed symmetric part (6)
   pure subroutine vjp_seeds(p, w0, w1, dense, w2)
      !> Column
      integer, intent(in) :: p
      !> Value weight
      real(wp), intent(out) :: w0
      !> Gradient weights
      real(wp), intent(out) :: w1(ndim)
      !> Non-symmetric Hessian seed
      real(wp), intent(out) :: dense(ndim, ndim)
      !> Packed symmetric part
      real(wp), intent(out) :: w2(6)

      integer :: a, b

      w0 = 0.35_wp + 0.1_wp*real(p, wp)
      do a = 1, ndim
         w1(a) = cos(1.3_wp*real(a, wp) + 0.7_wp*real(p, wp))
         do b = 1, ndim
            dense(a, b) = sin(0.9_wp*real(a, wp) - 1.7_wp*real(b, wp) + 0.5_wp*real(p, wp)) &
                          + 0.2_wp*real(a*b, wp)
         end do
      end do
      call pack_sym_matrix(0.5_wp*(dense + transpose(dense)), w2, ndim)
   end subroutine vjp_seeds

   !> Per-point `vjp_f1_rA` of the prepared point, scattered to candidate ids
   !>
   !> - candidate id of an active atom: `orig_to_sorted` (SvdW and CFC declare the
   !>   sorted space), independent of the base's `cand_to_user` inverse
   !>
   !> @param[inout] lsf    LSF, `prepare`d at the point
   !> @param[in]    w0     Value weight
   !> @param[in]    w1     Gradient weights
   !> @param[in]    seed   Packed symmetric Hessian seed (6), off-diagonals not doubled
   !> @param[out]   out    Result (3, ncenters), zero for inactive candidates
   !> @param[out]   active Candidate has a slot in the active list (ncenters)
   subroutine point_vjp_dense(lsf, w0, w1, seed, out, active)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Value weight
      real(wp), intent(in) :: w0
      !> Gradient weights
      real(wp), intent(in) :: w1(ndim)
      !> Packed symmetric Hessian seed
      real(wp), intent(in) :: seed(6)
      !> Result by candidate id
      real(wp), intent(out) :: out(:, :)
      !> Candidate is active
      logical, intent(out) :: active(:)

      real(wp), allocatable :: work(:, :)
      integer :: i, k

      allocate (work(ndim, lsf%ncenters))
      call lsf%vjp_f1_rA(w0, w1, seed, work)
      out = 0.0_wp
      active = .false.
      do i = 1, lsf%active_count()
         k = lsf%orig_to_sorted(lsf%active_atom(i))
         out(:, k) = work(:, i)
         active(k) = .true.
      end do
   end subroutine point_vjp_dense

   !> Batched vs per-point nuclear vjp, SvdW (fade band 4 / 3 e-folds)
   subroutine test_svdw_vjp_batch_vs_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_vs_point(error, "svdw", 12.0_wp, 4.0_wp)
   end subroutine test_svdw_vjp_batch_vs_point

   !> Batched vs per-point nuclear vjp, CFC (fade band 5 / 2 e-folds)
   subroutine test_cfc_vjp_batch_vs_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_vs_point(error, "cfc")
   end subroutine test_cfc_vjp_batch_vs_point

   !> `vjp_f1_rA_batch` against the per-point `vjp_f1_rA` on a tile of eight points
   !>
   !> - MB16-43 01, production threshold, fade on (CFC default band, SvdW 12 / 4
   !>   e-folds); the tile holds a plateau point (the fade changes nothing), a
   !>   point with an atom inside the band (the fade
   !>   moves the value), five more reached ones and a far unreached one
   !> - fixed non-symmetric seeds; the batch gets the packed symmetric part
   !> - reference: `prepare` + `vjp_f1_rA` with the packed seed, scattered
   !>   to candidate ids; list absent and the full list reversed, normwise
   !>   [[VJP_BATCH_REL]] per column; reached columns hold exact zeros at inactive
   !>   candidates; the far column is unreached and NaN
   !> - the per-point `vjp_f1_rA` with `pack_sym_part` of the non-symmetric seed
   !>   agrees too: only the symmetric part contributes
   !> - negative controls: one packed off-diagonal entry moved by 1e-6 misses the
   !>   unperturbed reference by at least [[VJP_BATCH_MISS]] and matches the same
   !>   perturbation of the per-point seed; the per-point seed with halved
   !>   off-diagonals (the upper triangle without the factor 2) misses
   subroutine run_vjp_batch_vs_point(error, kind, fade_start, fade_width)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> "svdw" or "cfc"
      character(len=*), intent(in) :: kind
      !> Fade onset in e-folds, the parameter default when absent
      real(wp), intent(in), optional :: fade_start
      !> Fade width in e-folds, the parameter default when absent
      real(wp), intent(in), optional :: fade_width

      integer, parameter :: n = 8, npool = 40
      !> Heap-proof fill of the result and extra columns past the tile
      real(wp), parameter :: sentinel = 12345.0_wp
      integer, parameter :: nextra = 2
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf, off
      real(wp) :: pool(ndim, npool), points(ndim, n), u(ndim), lo(ndim), hi(ndim)
      real(wp) :: o0(npool), f0(npool), diff
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), w2p(6, n), dense(ndim, ndim, n)
      real(wp) :: dense_packed(6), w2_half(6)
      real(wp), allocatable :: ref(:, :, :), refu(:, :, :)
      real(wp), allocatable :: refp(:, :, :), res(:, :, :), resl(:, :, :), resp(:, :, :)
      logical, allocatable :: act(:, :)
      integer, allocatable :: full(:)
      integer :: sp(npool), so(npool), st(n), stl(n), nsel, i, p, k, ip, ib, ncen
      real(wp) :: dev, devl, devp_ref, devp, devu, scale
      integer :: n_inactive
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, kind, screening_threshold=1.0e-11_wp, &
                    fade_start=fade_start, fade_width=fade_width)
      call init_lsf(off, mol, radii, 0, kind, screening_threshold=1.0e-11_wp, &
                    fade_enabled=.false.)
      ncen = lsf%ncenters

      ! Plateau and in-band points by what the fade does to the value
      state = 20261105
      lo = lsf%center_lower - 3.0_wp
      hi = lsf%center_upper + 3.0_wp
      do i = 1, npool
         call det_rand(u, state)
         pool(:, i) = lo + (hi - lo)*u
      end do
      call lsf%evaluate_batch(pool, 0, f0, error=lsf_err)
      call batch_status(lsf, sp)
      if (allocated(lsf_err)) then
         call test_failed(error, "pool evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call off%evaluate_batch(pool, 0, o0, error=lsf_err)
      call batch_status(off, so)
      if (allocated(lsf_err)) then
         call test_failed(error, "unfaded pool evaluate_batch failed: "//lsf_err%message)
         return
      end if
      ip = 0
      ib = 0
      do i = 1, npool
         if (sp(i) /= 0 .or. so(i) /= 0) cycle
         diff = abs(f0(i) - o0(i))/max(1.0_wp, abs(o0(i)))
         if (diff < 1.0e-10_wp .and. ip == 0) ip = i
         if (diff > 1.0e-8_wp .and. ib == 0) ib = i
      end do
      call check(error, ip > 0, more="no plateau point in the pool")
      if (allocated(error)) return
      call check(error, ib > 0, more="no in-band point in the pool")
      if (allocated(error)) return
      points(:, 1) = pool(:, ip)
      points(:, 2) = pool(:, ib)
      nsel = 2
      do i = 1, npool
         if (nsel == n - 1) exit
         if (sp(i) /= 0 .or. i == ip .or. i == ib) cycle
         nsel = nsel + 1
         points(:, nsel) = pool(:, i)
      end do
      call check(error, nsel == n - 1, more="the pool holds too few reached points")
      if (allocated(error)) return
      points(:, n) = lsf%center_upper + 200.0_wp

      allocate (ref(ndim, ncen, n), refu(ndim, ncen, n), &
                refp(ndim, ncen, n), act(ncen, n), &
                res(ndim, ncen, n + nextra), resl(ndim, ncen, n + nextra), &
                resp(ndim, ncen, n + nextra))
      res = sentinel
      resl = sentinel
      resp = sentinel
      ref = 0.0_wp
      refu = 0.0_wp
      refp = 0.0_wp
      act = .false.
      do p = 1, n
         call vjp_seeds(p, w0(p), w1(:, p), dense(:, :, p), w2(:, p))
      end do
      w2p = w2
      w2p(2, :) = w2p(2, :) + 1.0e-6_wp

      do p = 1, n - 1
         call lsf%prepare(points(:, p), lsf_err)
         if (allocated(lsf_err) .or. lsf%unreached()) then
            call test_failed(error, "per-point prepare failed on a reached column")
            return
         end if
         call point_vjp_dense(lsf, w0(p), w1(:, p), w2(:, p), ref(:, :, p), act(:, p))
         ! Off-diagonals halved: the upper triangle without the factor 2 (wrong)
         w2_half = w2(:, p)
         w2_half([2, 4, 5]) = 0.5_wp*w2_half([2, 4, 5])
         call point_vjp_dense(lsf, w0(p), w1(:, p), w2_half, refu(:, :, p), act(:, p))
         call point_vjp_dense(lsf, w0(p), w1(:, p), w2p(:, p), refp(:, :, p), act(:, p))
      end do
      n_inactive = count(.not. act(:, 1:n - 1))
      call check(error, n_inactive > 0, more="every candidate active: no zeros to check")
      if (allocated(error)) return

      ! List absent
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st(1:n - 1) == 0) .and. st(n) == 1, &
                 more="unreached(column) must hold on the far column only")
      if (allocated(error)) return
      call check(error, all(res(:, :, n + 1:) == sentinel), &
                 more="columns past the tile must keep the sentinel, list absent")
      if (allocated(error)) return
      if (.not. all(ieee_is_nan(res(:, :, n)))) then
         call test_failed(error, "the unreached column must be NaN")
         return
      end if
      ! Full list, reversed
      allocate (full(ncen))
      full = [(ncen + 1 - i, i=1, ncen)]
      stl = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resl, error=lsf_err, &
                               candidate_indices=full)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "vjp_f1_rA_batch with a list failed: "//lsf_err%message)
         return
      end if
      call check(error, all(stl == st), more="the full list must not move the unreached flags")
      if (allocated(error)) return
      call check(error, all(resl(:, :, n + 1:) == sentinel), &
                 more="columns past the tile must keep the sentinel, full list")
      if (allocated(error)) return

      dev = 0.0_wp
      devl = 0.0_wp
      do p = 1, n - 1
         scale = max(norm2(ref(:, :, p)), tiny(1.0_wp))
         dev = max(dev, norm2(res(:, :, p) - ref(:, :, p))/scale)
         devl = max(devl, norm2(resl(:, :, p) - ref(:, full, p))/scale)
         do k = 1, ncen
            if (act(k, p)) cycle
            if (any(res(:, k, p) /= 0.0_wp) .or. &
                any(resl(:, ncen + 1 - k, p) /= 0.0_wp)) then
               call test_failed(error, "a candidate out of reach must hold an exact zero")
               return
            end if
         end do
      end do
      call check(error, dev <= VJP_BATCH_REL, more="batch vs per-point vjp, list absent")
      if (allocated(error)) return
      call check(error, devl <= VJP_BATCH_REL, more="batch vs per-point vjp, full list")
      if (allocated(error)) return

      ! Negative controls
      devu = 0.0_wp
      do p = 1, n - 1
         devu = max(devu, maxval(abs(res(:, :, p) - refu(:, :, p))))
      end do
      call check(error, devu >= VJP_BATCH_MISS, &
                 more="the upper triangle without the factor 2 must miss the batch")
      if (allocated(error)) return

      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2p, resp, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "perturbed vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(resp(:, :, n + 1:) == sentinel), &
                 more="columns past the tile must keep the sentinel, perturbed seed")
      if (allocated(error)) return
      devp = 0.0_wp
      devp_ref = 0.0_wp
      do p = 1, n - 1
         devp = max(devp, maxval(abs(resp(:, :, p) - ref(:, :, p))))
         devp_ref = max(devp_ref, norm2(resp(:, :, p) - refp(:, :, p)) &
                        /max(norm2(refp(:, :, p)), tiny(1.0_wp)))
      end do
      call check(error, devp >= VJP_BATCH_MISS, &
                 more="a 1e-6 packed perturbation must miss the unperturbed reference")
      if (allocated(error)) return
      call check(error, devp_ref <= VJP_BATCH_REL, &
                 more="the perturbed batch must match the perturbed per-point seed")
   end subroutine run_vjp_batch_vs_point

   !> Candidate subset of `vjp_f1_rA_batch`, SvdW (fade band 4 / 3 e-folds)
   subroutine test_svdw_vjp_batch_candidate_subset(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_candidate_subset(error, "svdw", 4.0_wp, 3.0_wp)
   end subroutine test_svdw_vjp_batch_candidate_subset

   !> Candidate subset of `vjp_f1_rA_batch`, CFC (fade band 5 / 2 e-folds)
   subroutine test_cfc_vjp_batch_candidate_subset(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_candidate_subset(error, "cfc", 5.0_wp, 2.0_wp)
   end subroutine test_cfc_vjp_batch_candidate_subset

   !> A strict candidate list vs the full list, column by column
   !>
   !> - 1000-atom lattice, six points in a corner cube, all reached
   !> - the list keeps every candidate any point reaches plus every seventh one
   !>   (out of reach of all points), so the jet is the full one and columns at
   !>   listed candidates equal the full-list columns at [[VJP_BATCH_REL]] (the
   !>   full screen and the listed one sum the power sums in different orders)
   !> - the list reversed: the columns come out in the list's own order, bit for
   !>   bit those of the sorted list (the sort inside makes the order invisible)
   !> - `res` carries two extra tile columns holding a sentinel: untouched, nothing
   !>   is written past column `np`
   !> - negative control: dropping a reached candidate from the list moves the other
   !>   columns by at least [[VJP_BATCH_MISS]], so the equality above is not vacuous
   subroutine run_vjp_batch_candidate_subset(error, kind, fade_start, fade_width)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> "svdw" or "cfc"
      character(len=*), intent(in) :: kind
      !> Fade onset in e-folds
      real(wp), intent(in) :: fade_start
      !> Fade width in e-folds
      real(wp), intent(in) :: fade_width

      integer, parameter :: n = 6
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), dense(ndim, ndim, n)
      real(wp), allocatable :: full(:, :, :), sub(:, :, :), rev(:, :, :), drop(:, :, :)
      logical, allocatable :: used(:)
      integer, allocatable :: list(:), rlist(:), dlist(:)
      integer :: st(n), ncen, nsub, nused, i, j, p, idrop
      real(wp) :: miss, big, dev
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call lattice_structure(mol, 10, 3.0_wp, 7)
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, kind, screening_threshold=1.0e-11_wp, &
                    fade_enabled=.true., fade_start=fade_start, fade_width=fade_width)
      ncen = lsf%ncenters
      state = 99
      lo = lsf%center_lower + 1.0_wp
      hi = lsf%center_lower + 5.0_wp
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lo + (hi - lo)*u
         call vjp_seeds(i, w0(i), w1(:, i), dense(:, :, i), w2(:, i))
      end do

      allocate (full(ndim, ncen, n))
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, full, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "full-list vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0), more="every corner-cube point must be reached")
      if (allocated(error)) return

      allocate (used(ncen))
      do j = 1, ncen
         used(j) = any(full(:, j, :) /= 0.0_wp)
      end do
      nused = count(used)
      list = pack([(j, j=1, ncen)], used .or. mod([(j, j=1, ncen)], 7) == 0)
      nsub = size(list)
      call check(error, nused > 0 .and. nsub > nused .and. nsub < ncen, &
                 more="the list must be a strict superset of the reached candidates")
      if (allocated(error)) return

      allocate (sub(ndim, nsub, n + 2), rev(ndim, nsub, n + 2))
      sub = sentinel
      rev = sentinel
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, sub, error=lsf_err, &
                               candidate_indices=list)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "subset vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0), more="the subset must keep every column reached")
      if (allocated(error)) return
      rlist = list(nsub:1:-1)
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, rev, error=lsf_err, &
                               candidate_indices=rlist)
      if (allocated(lsf_err)) then
         call test_failed(error, "reversed vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      dev = 0.0_wp
      do p = 1, n
         dev = max(dev, norm2(sub(:, :, p) - full(:, list, p)) &
                   /max(norm2(full(:, list, p)), tiny(1.0_wp)))
         do j = 1, nsub
            if (.not. used(list(j)) .and. any(sub(:, j, p) /= 0.0_wp)) then
               call test_failed(error, "a candidate out of reach must hold an exact zero")
               return
            end if
            if (any(rev(:, j, p) /= sub(:, nsub + 1 - j, p))) then
               call test_failed(error, "reversed list changes the bits of a column")
               return
            end if
         end do
      end do
      call check(error, dev <= VJP_BATCH_REL, more="subset vs full-list columns")
      if (allocated(error)) return
      if (any(sub(:, :, n + 1:n + 2) /= sentinel) .or. &
          any(rev(:, :, n + 1:n + 2) /= sentinel)) then
         call test_failed(error, "vjp_f1_rA_batch wrote past the last tile column")
         return
      end if

      ! Negative control: drop the reached candidate of the largest column
      big = -1.0_wp
      idrop = 0
      do j = 1, ncen
         if (.not. used(j)) cycle
         if (maxval(abs(full(:, j, :))) > big) then
            big = maxval(abs(full(:, j, :)))
            idrop = j
         end if
      end do
      dlist = pack(list, list /= idrop)
      allocate (drop(ndim, size(dlist), n))
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, drop, error=lsf_err, &
                               candidate_indices=dlist)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "dropped-candidate batch failed: "//lsf_err%message)
         return
      end if
      miss = 0.0_wp
      do p = 1, n
         do j = 1, size(dlist)
            if (st(p) /= 0) cycle
            miss = max(miss, maxval(abs(drop(:, j, p) - full(:, dlist(j), p))))
         end do
      end do
      call check(error, miss >= VJP_BATCH_MISS, &
                 more="dropping a reached candidate must move the other columns")
   end subroutine run_vjp_batch_candidate_subset

   !> A failed `prepare` column: `error` names it and the call returns (SvdW)
   !>
   !> - fade width 0 is refused by `prepare` ("must be positive"), at every column
   !> - as `evaluate_batch`: `error` names column 1 and the call returns; the
   !>   unreached mask is unspecified after an error and not checked
   !> - negative control: with a valid fade the same call raises nothing and
   !>   leaves every column reached
   subroutine test_svdw_vjp_batch_failed_column(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param, valid
      real(wp) :: points(ndim, 2), w0(2), w1(ndim, 2), w2(6, 2), dense(ndim, ndim)
      real(wp), allocatable :: res(:, :, :)
      integer :: st(2), p
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      points(:, 1) = mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp]
      points(:, 2) = mol%xyz(:, 2) + [0.0_wp, 2.0_wp, 0.0_wp]
      do p = 1, 2
         call vjp_seeds(p, w0(p), w1(:, p), dense, w2(:, p))
      end do

      call param%new(fade_enabled=.true., fade_width=0.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%set_max_deriv(2)
      allocate (res(ndim, lsf%ncenters, 2))
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
      call check(error, allocated(lsf_err), more="a failed column must set error")
      if (allocated(error)) return
      call check(error, &
                 index(lsf_err%message, "vjp_f1_rA_batch failed at point 1: ") > 0, &
                 more="error must name column 1: "//lsf_err%message)
      if (allocated(error)) return

      call valid%new(fade_enabled=.true.)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, valid)
      call lsf%set_max_deriv(2)
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
      call batch_status(lsf, st)
      call check(error, .not. allocated(lsf_err), more="a valid fade must not raise")
      if (allocated(error)) return
      call check(error, all(st == 0), more="a valid fade must leave every column reached")
   end subroutine test_svdw_vjp_batch_failed_column

   !> Shape descriptions of `vjp_f1_rA_batch`: each misfit named, a fit silent
   !>
   !> The `error stop` on a nonempty message is not testable in process; the
   !> message function is the single place deciding it
   subroutine test_vjp_batch_shape_messages(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: p3(3, 4), p2(2, 4), w0(4), w0s(3), w1(3, 4), w1b(2, 4), w1s(3, 3)
      real(wp) :: w2(6, 4), w2b(5, 4), w2s(6, 3), r(3, 6, 4), rc(3, 5, 4), rn(3, 6, 3)
      real(wp) :: rd(2, 6, 4), ra(3, 6, 6)

      p3 = 0.0_wp
      p2 = 0.0_wp
      w0 = 0.0_wp
      w0s = 0.0_wp
      w1 = 0.0_wp
      w1b = 0.0_wp
      w1s = 0.0_wp
      w2 = 0.0_wp
      w2b = 0.0_wp
      w2s = 0.0_wp
      r = 0.0_wp
      rc = 0.0_wp
      rn = 0.0_wp
      rd = 0.0_wp
      ra = 0.0_wp

      call check(error, len(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2, r, 6)), &
                 0, more="a fitting call must give no message")
      if (allocated(error)) return
      call check(error, len(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2, ra, 6)), 0, &
                 more="extra tile columns in res must fit")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p2, w0, w1, w2, r, 6), &
                              "points must be (3, np)") > 0, more="points leading extent")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0s, w1, w2, r, 6), &
                              "w0 holds fewer") > 0, more="w0 length")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1b, w2, r, 6), &
                              "w1 must be (3, >= np)") > 0, more="w1 leading extent")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1s, w2, r, 6), &
                              "w1 must be (3, >= np)") > 0, more="w1 columns")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2b, r, 6), &
                              "w2 must be packed (6, >= np)") > 0, &
                 more="w2 packed length")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2s, r, 6), &
                              "w2 must be packed (6, >= np)") > 0, more="w2 columns")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2, rc, 6), &
                              "res must be (3, ncand, >= np)") > 0, more="res candidates")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2, rn, 6), &
                              "res must be (3, ncand, >= np)") > 0, more="res columns")
      if (allocated(error)) return
      call check(error, index(lsf_vjp_batch_shape_message("X", p3, w0, w1, w2, rd, 6), &
                              "res must be (3, ncand, >= np)") > 0, &
                 more="res leading extent")
   end subroutine test_vjp_batch_shape_messages

   !> Failed column raises `error` in `vjp_f1_rA_batch` (SvdW)
   !>
   !> Fade width 0 is refused by `prepare` at column 1
   subroutine test_bad_svdw_vjp_batch_failed_column(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), res(:, :, :)
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: points(ndim, 1), w0(1), w1(ndim, 1), w2(6, 1), dense(ndim, ndim)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_width=0.0_wp)
      call init_svdw_param(lsf, mol, radii, 1.0e-11_wp, param)
      call lsf%set_max_deriv(2)
      points(:, 1) = mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp]
      call vjp_seeds(1, w0(1), w1(:, 1), dense, w2(:, 1))
      allocate (res(ndim, lsf%ncenters, 1))
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
      if (.not. allocated(lsf_err)) return
      if (index(lsf_err%message, "failed at point 1") > 0) then
         call test_failed(error, lsf_err%message)
      end if
   end subroutine test_bad_svdw_vjp_batch_failed_column

   !> Failed column raises `error` in `vjp_f1_rA_batch` (CFC)
   !>
   !> Fade width 0 is refused by the override's own `check_batch` at column 1
   subroutine test_bad_cfc_vjp_batch_check_batch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), res(:, :, :)
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp) :: points(ndim, 1), w0(1), w1(ndim, 1), w2(6, 1), dense(ndim, ndim)
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call param%new(fade_enabled=.true., fade_width=0.0_wp)
      call init_cfc_param(lsf, mol, radii, 1.0e-11_wp, param)
      points(:, 1) = mol%xyz(:, 1) + [2.0_wp, 0.0_wp, 0.0_wp]
      call vjp_seeds(1, w0(1), w1(:, 1), dense, w2(:, 1))
      allocate (res(ndim, lsf%ncenters, 1))
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
      if (.not. allocated(lsf_err)) return
      if (index(lsf_err%message, "failed at point 1") > 0) then
         call test_failed(error, lsf_err%message)
      end if
   end subroutine test_bad_cfc_vjp_batch_check_batch

   !> Non-finite adjoint inside, a sub-normal distance off a nucleus (SvdW)
   subroutine test_bad_svdw_vjp_batch_near_nucleus(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_bad_vjp_near_nucleus(error, "svdw")
   end subroutine test_bad_svdw_vjp_batch_near_nucleus

   !> Non-finite adjoint inside, a sub-normal distance off a nucleus (CFC)
   subroutine test_bad_cfc_vjp_batch_near_nucleus(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_bad_vjp_near_nucleus(error, "cfc")
   end subroutine test_bad_cfc_vjp_batch_near_nucleus

   !> A non-finite vjp column with `S <= 0` is an `error`, never unreached
   !>
   !> - two carbons, `R = 3`, the first at the origin so a tiny offset is
   !>   representable; a reached control 2 Bohr off it, then a point 1e-100 off it:
   !>   no zero distance (the on-nucleus path is not taken), the lift finite, the
   !>   adjoint's powers of `1/d` overflow (measured onset 1e-60..1e-70 for both)
   !>   while `S <= 0`
   !> - expected: `error` naming point 2 as not finite with `S <= 0`; a clean return
   !>   (the column marked unreached) fails the test
   !>
   !> @param[out] error Test failure, set only on the expected error
   !> @param[in]  kind  "svdw" or "cfc"
   subroutine run_bad_vjp_near_nucleus(error, kind)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> "svdw" or "cfc"
      character(len=*), intent(in) :: kind

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), res(:, :, :)
      real(wp) :: points(ndim, 2), w0(2), w1(ndim, 2), w2(6, 2), dense(ndim, ndim)
      integer :: p
      type(mctc_error), allocatable :: lsf_err

      call new(mol, [6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, 2.5_wp, 0.0_wp, 0.0_wp], &
                                    [ndim, 2]))
      radii = [3.0_wp, 3.0_wp]
      call init_lsf(lsf, mol, radii, 2, kind)
      points(:, 1) = [-2.0_wp, 0.0_wp, 0.0_wp]
      do p = 1, 2
         call vjp_seeds(p, w0(p), w1(:, p), dense, w2(:, p))
      end do
      allocate (res(ndim, lsf%ncenters, 2))
      points(:, 2) = [1.0e-100_wp, 0.0_wp, 0.0_wp]
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, lsf_err)
      call expect_error(error, lsf_err, "point 2 has a level-set jet not finite with S <= 0")
   end subroutine run_bad_vjp_near_nucleus

   !> Per-point `vjp_f1_rA` of every reached column of a tile, by candidate id
   !>
   !> - `prepare` at the point, then [[point_vjp_dense]] with the packed seed `w2`
   !> - an unreached column, or one `prepare` refuses (a CFC nucleus), is left at
   !>   zero with `reached` false
   !>
   !> @param[inout] lsf     LSF with `max_deriv >= 2`
   !> @param[in]    points  Tile points (3, n)
   !> @param[in]    w0      Value weights (n)
   !> @param[in]    w1      Gradient weights (3, n)
   !> @param[in]    w2      Packed symmetric Hessian seeds (6, n)
   !> @param[out]   ref     Per-point result (3, ncenters, n)
   !> @param[out]   active  Candidate has an active slot (ncenters, n)
   !> @param[out]   reached Column reached an atom (n)
   subroutine point_vjp_tile(lsf, points, w0, w1, w2, ref, active, reached)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Tile points
      real(wp), intent(in) :: points(:, :)
      !> Value weights
      real(wp), intent(in) :: w0(:)
      !> Gradient weights
      real(wp), intent(in) :: w1(:, :)
      !> Packed symmetric Hessian seeds
      real(wp), intent(in) :: w2(:, :)
      !> Per-point result
      real(wp), intent(out) :: ref(:, :, :)
      !> Candidate has an active slot
      logical, intent(out) :: active(:, :)
      !> Column reached an atom
      logical, intent(out) :: reached(:)

      type(mctc_error), allocatable :: lsf_err
      integer :: p

      ref = 0.0_wp
      active = .false.
      do p = 1, size(points, 2)
         call lsf%prepare(points(:, p), lsf_err)
         reached(p) = .not. allocated(lsf_err) .and. .not. lsf%unreached()
         if (.not. reached(p)) cycle
         call point_vjp_dense(lsf, w0(p), w1(:, p), w2(:, p), ref(:, :, p), active(:, p))
      end do
   end subroutine point_vjp_tile

   !> Normwise relative deviation of a tile column from a reference column
   !>
   !> @param[in] a Column under test
   !> @param[in] b Reference column
   !> @returns   dev `norm2(a - b)/norm2(b)`, floored at `tiny`
   pure function column_dev(a, b) result(dev)
      !> Column under test
      real(wp), intent(in) :: a(:, :)
      !> Reference column
      real(wp), intent(in) :: b(:, :)
      !> Normwise relative deviation
      real(wp) :: dev

      dev = norm2(a - b)/max(norm2(b), tiny(1.0_wp))
   end function column_dev

   !> Atoms of the plateau, in-band and beyond-band classes against a point box
   !>
   !> The classes of stage 1's strip against the box `[lo, hi]`, from the
   !> signed sphere distance `d = |r - R_A| - R_a` at the nearest and farthest
   !> box point: plateau `d_max <= x_on`, beyond `d_min >= x_on + delta`, else in band
   !>
   !> @param[in]  mol     Structure
   !> @param[in]  radii   Per-atom radii
   !> @param[in]  lo      Lower corner of the box
   !> @param[in]  hi      Upper corner of the box
   !> @param[in]  x_on    Fade onset (Bohr)
   !> @param[in]  delta   Fade width (Bohr)
   !> @param[out] counts  Atoms in the plateau, band and beyond classes (3)
   pure subroutine fade_class_counts(mol, radii, lo, hi, x_on, delta, counts)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Lower corner of the box
      real(wp), intent(in) :: lo(ndim)
      !> Upper corner of the box
      real(wp), intent(in) :: hi(ndim)
      !> Fade onset
      real(wp), intent(in) :: x_on
      !> Fade width
      real(wp), intent(in) :: delta
      !> Atoms per class
      integer, intent(out) :: counts(3)

      real(wp) :: e(ndim), f(ndim), dnear, dfar
      integer :: a

      counts = 0
      do a = 1, mol%nat
         e = max(lo - mol%xyz(:, a), 0.0_wp, mol%xyz(:, a) - hi)
         f = max(abs(mol%xyz(:, a) - lo), abs(mol%xyz(:, a) - hi))
         dnear = norm2(e) - radii(a)
         dfar = norm2(f) - radii(a)
         if (dfar <= x_on) then
            counts(1) = counts(1) + 1
         else if (dnear >= x_on + delta) then
            counts(3) = counts(3) + 1
         else
            counts(2) = counts(2) + 1
         end if
      end do
   end subroutine fade_class_counts

   !> Batched vs per-point nuclear vjp across three strips (SvdW)
   !>
   !> - 1000-atom lattice, fade band 4 / 3 e-folds, 70 points drifting along x so
   !>   the three strips (32 + 32 + 6) see different atoms; point 45, in the
   !>   middle of the second strip, is far out: unreached and NaN, the
   !>   columns around it computed
   !> - reference: `prepare` + `vjp_f1_rA` per point; the list absent (per-strip
   !>   box lists), an ascending list (written in place by the kernel) and the same
   !>   list reversed (scattered), normwise [[VJP_BATCH_REL]]; exact zeros at
   !>   candidates without an active slot
   !> - tiling: the tile reversed puts every point in another strip with other
   !>   mates; the columns agree at [[VJP_BATCH_REL]], not bitwise
   !> - negative control: dropping the largest reached candidate from the list moves
   !>   the other columns by at least [[VJP_BATCH_MISS]]
   subroutine test_svdw_vjp_batch_strips(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 70, iun = 45
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), rpoints(ndim, n), u(ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), rw0(n), rw1(ndim, n), rw2(6, n)
      real(wp) :: dense(ndim, ndim)
      real(wp), allocatable :: ref(:, :, :), res(:, :, :), resr(:, :, :), resl(:, :, :)
      real(wp), allocatable :: resv(:, :, :), resd(:, :, :)
      logical, allocatable :: act(:, :), used(:)
      logical :: reached(n)
      integer, allocatable :: list(:), rlist(:), dlist(:)
      integer :: st(n), str(n), stl(n), ncen, nsub, i, j, p, idrop
      real(wp) :: dev, devl, devt, miss, big
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call check(error, n > 2*strip_width_svdw(2) .and. n <= 3*strip_width_svdw(2), &
                 more="the tile must cross exactly three strips")
      if (allocated(error)) return
      call lattice_structure(mol, 10, 3.0_wp, 31)
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, "svdw", screening_threshold=1.0e-11_wp, &
                    fade_enabled=.true., fade_start=4.0_wp, fade_width=3.0_wp)
      ncen = lsf%ncenters

      state = 4711
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lsf%center_lower + 4.0_wp*u &
                        + [1.0_wp + 0.12_wp*real(i, wp), 1.0_wp, 1.0_wp]
         call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
      end do
      points(:, iun) = lsf%center_upper + 200.0_wp
      do i = 1, n
         rpoints(:, i) = points(:, n + 1 - i)
         rw0(i) = w0(n + 1 - i)
         rw1(:, i) = w1(:, n + 1 - i)
         rw2(:, i) = w2(:, n + 1 - i)
      end do

      allocate (ref(ndim, ncen, n), act(ncen, n), res(ndim, ncen, n), &
                resr(ndim, ncen, n), resl(ndim, ncen, n), used(ncen))
      call point_vjp_tile(lsf, points, w0, w1, w2, ref, act, reached)
      call check(error, count(.not. reached) == 1 .and. .not. reached(iun), &
                 more="only the far column may be unreached per point")
      if (allocated(error)) return

      ! List absent: per-strip box lists
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == merge(1, 0, [(i == iun, i=1, n)])), &
                 more="unreached on the far column only")
      if (allocated(error)) return
      if (.not. all(ieee_is_nan(res(:, :, iun)))) then
         call test_failed(error, "the unreached column must be NaN")
         return
      end if
      dev = 0.0_wp
      do p = 1, n
         if (p == iun) cycle
         dev = max(dev, column_dev(res(:, :, p), ref(:, :, p)))
         do j = 1, ncen
            if (.not. act(j, p) .and. any(res(:, j, p) /= 0.0_wp)) then
               call test_failed(error, "a candidate out of reach must hold an exact zero")
               return
            end if
         end do
      end do
      call check(error, dev <= VJP_BATCH_REL, &
                 more="batch vs per-point vjp over three strips")
      if (allocated(error)) return

      ! Ascending strict list (in place) and its reverse (scattered)
      do j = 1, ncen
         used(j) = any(act(j, :))
      end do
      list = pack([(j, j=1, ncen)], used .or. mod([(j, j=1, ncen)], 7) == 0)
      nsub = size(list)
      call check(error, count(used) > 0 .and. nsub > count(used) .and. nsub < ncen, &
                 more="the list must be a strict superset of the reached candidates")
      if (allocated(error)) return
      allocate (resv(ndim, nsub, n + 2), resd(ndim, nsub, n + 2))
      resv = sentinel
      resd = sentinel
      stl = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resv, error=lsf_err, &
                               candidate_indices=list)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "list vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(stl == st), more="the list must not move the unreached flags")
      if (allocated(error)) return
      rlist = list(nsub:1:-1)
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resd, error=lsf_err, &
                               candidate_indices=rlist)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "reversed vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      devl = 0.0_wp
      do p = 1, n
         if (p == iun) cycle
         devl = max(devl, column_dev(resv(:, :, p), ref(:, list, p)))
         do j = 1, nsub
            if (any(resd(:, j, p) /= resv(:, nsub + 1 - j, p))) then
               call test_failed(error, "reversed list changes the bits of a column")
               return
            end if
         end do
      end do
      call check(error, devl <= VJP_BATCH_REL, more="in-place list vs per-point vjp")
      if (allocated(error)) return
      if (any(resv(:, :, n + 1:n + 2) /= sentinel) .or. &
          any(resd(:, :, n + 1:n + 2) /= sentinel)) then
         call test_failed(error, "vjp_f1_rA_batch wrote past the last tile column")
         return
      end if

      ! Tiling: the reversed tile, every point in another strip
      str = -1
      call lsf%vjp_f1_rA_batch(rpoints, rw0, rw1, rw2, resr, error=lsf_err)
      call batch_status(lsf, str)
      if (allocated(lsf_err)) then
         call test_failed(error, "reversed-tile vjp_f1_rA_batch failed: "// &
                          lsf_err%message)
         return
      end if
      call check(error, all(str(n:1:-1) == st), &
                 more="the tile order must not move the unreached flags")
      if (allocated(error)) return
      devt = 0.0_wp
      do p = 1, n
         if (p == iun) cycle
         devt = max(devt, column_dev(res(:, :, p), resr(:, :, n + 1 - p)))
      end do
      call check(error, devt <= VJP_BATCH_REL, &
                 more="the same point in two different tiles")
      if (allocated(error)) return

      ! Negative control: drop the reached candidate of the largest column
      big = -1.0_wp
      idrop = 0
      do j = 1, nsub
         if (.not. used(list(j))) cycle
         if (maxval(abs(resv(:, j, 1:n))) > big) then
            big = maxval(abs(resv(:, j, 1:n)))
            idrop = list(j)
         end if
      end do
      dlist = pack(list, list /= idrop)
      deallocate (resd)
      allocate (resd(ndim, size(dlist), n))
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resd, error=lsf_err, &
                               candidate_indices=dlist)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "dropped-candidate batch failed: "//lsf_err%message)
         return
      end if
      miss = 0.0_wp
      do p = 1, n
         if (stl(p) /= 0) cycle
         do j = 1, size(dlist)
            miss = max(miss, maxval(abs(resd(:, j, p) - ref(:, dlist(j), p))))
         end do
      end do
      call check(error, miss >= VJP_BATCH_MISS, &
                 more="dropping a reached candidate must move the other columns")
   end subroutine test_svdw_vjp_batch_strips

   !> Points on a nucleus inside a vjp tile (SvdW)
   !>
   !> - MB16-43 01, twelve points near atom 1, two of them exactly on atoms 3 and
   !>   7; the list absent, then the full ascending list (the in-place mode)
   !> - on-nucleus columns: [[VJP_BATCH_REL]] against the per-point `prepare` +
   !>   `vjp_f1_rA`, finite, reached
   !> - the other columns equal, bit for bit, those of the tile without the two
   !>   on-nucleus columns (the strip is re-run without them)
   !> - negative controls: a column of the same tile moved off the nucleus by 1e-3
   !>   misses the on-nucleus reference by at least [[VJP_BATCH_MISS]]; a mate with
   !>   a changed seed misses its bit-for-bit twin
   subroutine test_svdw_vjp_batch_on_nucleus(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 12, nuc(2) = [5, 9]
      !> Heap-proof fill of the results, kept by the columns past the tile
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), pn(ndim, n - 2), pm(ndim, n), u(ndim)
      real(wp) :: dense(ndim, ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), v0(n - 2), v1(ndim, n - 2), v2(6, n - 2)
      real(wp) :: x0(n)
      real(wp), allocatable :: ref(:, :, :), res(:, :, :), resn(:, :, :), resm(:, :, :)
      real(wp), allocatable :: refm(:, :, :)
      integer, allocatable :: full(:)
      logical, allocatable :: act(:, :), actm(:, :)
      logical :: reached(n), reachedm(n)
      integer :: st(n), stn(n - 2), mode, i, j, k, ncen
      real(wp) :: dev, miss, mate_miss
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, "svdw", screening_threshold=1.0e-11_wp)
      ncen = lsf%ncenters
      state = 17
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = mol%xyz(:, 1) + 3.0_wp*(u - 0.5_wp)
         call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
      end do
      points(:, nuc(1)) = mol%xyz(:, 3)
      points(:, nuc(2)) = mol%xyz(:, 7)
      k = 0
      do i = 1, n
         if (i == nuc(1) .or. i == nuc(2)) cycle
         k = k + 1
         pn(:, k) = points(:, i)
         v0(k) = w0(i)
         v1(:, k) = w1(:, i)
         v2(:, k) = w2(:, i)
      end do
      pm = points
      pm(:, nuc(1)) = mol%xyz(:, 3) + [1.0e-3_wp, 0.0_wp, 0.0_wp]
      x0 = w0
      x0(2) = x0(2) + 1.0e-3_wp

      allocate (ref(ndim, ncen, n), refm(ndim, ncen, n), res(ndim, ncen, n + 1), &
                resn(ndim, ncen, n - 2), resm(ndim, ncen, n), act(ncen, n), actm(ncen, n))
      call point_vjp_tile(lsf, points, w0, w1, w2, ref, act, reached)
      call point_vjp_tile(lsf, pm, w0, w1, w2, refm, actm, reachedm)
      call check(error, all(reached) .and. all(reachedm), &
                 more="every column must be reached")
      if (allocated(error)) return
      full = [(j, j=1, ncen)]

      do mode = 1, 2
         st = -1
         stn = -1
         res = sentinel
         resn = sentinel
         resm = sentinel
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
            call batch_status(lsf, st)
         else
            call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err, &
                                     candidate_indices=full)
            call batch_status(lsf, st)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
            return
         end if
         call check(error, all(st == 0), more="an on-nucleus column must stay reached")
         if (allocated(error)) return
         call check(error, all(res(:, :, n + 1) == sentinel), &
                    more="the column past the tile must keep the sentinel")
         if (allocated(error)) return
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(pn, v0, v1, v2, resn, error=lsf_err)
            call batch_status(lsf, stn)
         else
            call lsf%vjp_f1_rA_batch(pn, v0, v1, v2, resn, error=lsf_err, &
                                     candidate_indices=full)
            call batch_status(lsf, stn)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "tile without the nucleus failed: "//lsf_err%message)
            return
         end if
         if (.not. all(ieee_is_finite(res(:, :, 1:n)))) then
            call test_failed(error, "an on-nucleus tile must stay finite")
            return
         end if

         dev = 0.0_wp
         do i = 1, n
            dev = max(dev, column_dev(res(:, :, i), ref(:, :, i)))
            do k = 1, ncen
               if (.not. act(k, i) .and. any(res(:, k, i) /= 0.0_wp)) then
                  call test_failed(error, &
                                   "a candidate out of reach must hold an exact zero")
                  return
               end if
            end do
         end do
         call check(error, dev <= VJP_BATCH_REL, &
                    more="tile with on-nucleus columns vs per-point vjp")
         if (allocated(error)) return

         k = 0
         do i = 1, n
            if (i == nuc(1) .or. i == nuc(2)) cycle
            k = k + 1
            if (any(res(:, :, i) /= resn(:, :, k))) then
               call test_failed(error, &
                                "an on-nucleus column moved the bits of a strip-mate")
               return
            end if
         end do

         ! Negative controls
         miss = maxval(abs(res(:, :, nuc(1)) - refm(:, :, nuc(1))))
         call check(error, miss >= VJP_BATCH_MISS, &
                    more="a point 1e-3 off the nucleus must miss the on-nucleus column")
         if (allocated(error)) return
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(points, x0, w1, w2, resm, error=lsf_err)
            call batch_status(lsf, st)
         else
            call lsf%vjp_f1_rA_batch(points, x0, w1, w2, resm, error=lsf_err, &
                                     candidate_indices=full)
            call batch_status(lsf, st)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "changed-seed tile failed: "//lsf_err%message)
            return
         end if
         mate_miss = maxval(abs(resm(:, :, 2) - res(:, :, 2)))
         call check(error, mate_miss >= VJP_BATCH_MISS, &
                    more="a changed seed must move its own column")
         if (allocated(error)) return
      end do
   end subroutine test_svdw_vjp_batch_on_nucleus

   !> Maximum derivative order the SvdW LSF is prepared for, -1 for another LSF
   !>
   !> @param[in] lsf LSF instance
   !> @returns   md  `max_deriv` of an SvdW LSF
   function svdw_max_deriv(lsf) result(md)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Prepared derivative order
      integer :: md

      md = -1
      select type (lsf)
      type is (moist_cavity_drop_lsf_svdw_type)
         md = lsf%max_deriv
      end select
   end function svdw_max_deriv

   !> On-nucleus columns across strips of one vjp tile (SvdW)
   !>
   !> - MB16-43 01 and a 1000-atom lattice (inactive candidates on the on-nucleus
   !>   path), 70 points in three strips (32 + 32 + 6) near atom 1, with
   !>   `max_deriv = 3` so the on-nucleus fallback has an order to restore
   !> - case 1: on-nucleus columns 5 and 40, in two different strips
   !> - case 2: on-nucleus column 5 next to an unreached column 6 in one strip
   !> - case 3: every column of the second strip on an atom (no strip-mate left to
   !>   re-run)
   !> - each case with the list absent, the full ascending list (in place) and the
   !>   full list reversed (scatter); the sentinel-filled result keeps it past the
   !>   tile, holds exact zeros at inactive candidates, and `max_deriv` is the same
   !>   after the call
   !> - reference: the per-point `prepare` + `vjp_f1_rA`, normwise
   !>   [[VJP_BATCH_REL]]; reached, or unreached and NaN for the far column
   !> - negative control: the first on-nucleus point moved off by 1e-3 gives another
   !>   reference, which the on-nucleus column misses by at least [[VJP_BATCH_MISS]]
   subroutine test_svdw_vjp_batch_on_nucleus_strips(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 70, ncase = 3, nextra = 2
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: base(ndim, n), points(ndim, n), pmove(ndim, 1), u(ndim)
      real(wp) :: dense(ndim, ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n)
      real(wp), allocatable :: ref(:, :, :), refm(:, :, :), res(:, :, :)
      logical, allocatable :: act(:, :), actm(:, :)
      logical :: reached(n), reachedm(1), far(n)
      integer, allocatable :: idx(:)
      integer :: atom_of(n), st(n), want(n), ifix, icase, mode, i, j, k, c, ncen, md0, md1
      real(wp) :: dev, miss
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      do ifix = 1, 2
         if (allocated(ref)) deallocate (ref, refm, res, act, actm)
         if (ifix == 1) then
            call get_structure(mol, "MB16-43", "01")
         else
            call lattice_structure(mol, 10, 3.0_wp, 31)
         end if
         call get_test_radii(mol, radii)
         call init_lsf(lsf, mol, radii, 3, "svdw", screening_threshold=1.0e-11_wp)
         ncen = lsf%ncenters
         state = 29
         do i = 1, n
            call det_rand(u, state)
            base(:, i) = mol%xyz(:, 1) + 3.0_wp*(u - 0.5_wp)
            call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
         end do
         allocate (ref(ndim, ncen, n), refm(ndim, ncen, 1), &
                   res(ndim, ncen, n + nextra), act(ncen, n), actm(ncen, 1))

         do icase = 1, ncase
            atom_of = 0
            far = .false.
            select case (icase)
            case (1)
               atom_of(5) = 3
               atom_of(40) = 7
            case (2)
               atom_of(5) = 3
               far(6) = .true.
            case default
               atom_of(33:64) = [(1 + mod(j - 1, mol%nat), j=1, 32)]
            end select
            points = base
            do i = 1, n
               if (atom_of(i) > 0) points(:, i) = mol%xyz(:, atom_of(i))
               if (far(i)) points(:, i) = lsf%center_upper + 200.0_wp
            end do
            want = merge(1, 0, far)
            c = findloc(atom_of > 0, .true., dim=1)

            call point_vjp_tile(lsf, points, w0, w1, w2, ref, act, reached)
            call check(error, all(reached .eqv. .not. far), &
                       more="the per-point reference must reach exactly the near columns")
            if (allocated(error)) return
            ! The lattice must leave inactive candidates on the on-nucleus path
            call check(error, ifix == 1 .or. .not. all(act(:, c)), &
                       more="an on-nucleus lattice column must have inactive candidates")
            if (allocated(error)) return
            ! Negative control reference: the first on-nucleus point displaced
            pmove(:, 1) = points(:, c) + [1.0e-3_wp, 0.0_wp, 0.0_wp]
            call point_vjp_tile(lsf, pmove, w0(c:c), w1(:, c:c), w2(:, c:c), refm, actm, &
                                reachedm)
            call check(error, reachedm(1), more="the displaced point must be reached")
            if (allocated(error)) return

            do mode = 1, 3
               if (mode == 3) then
                  idx = [(ncen + 1 - j, j=1, ncen)]
               else
                  idx = [(j, j=1, ncen)]
               end if
               res = sentinel
               st = -1
               md0 = svdw_max_deriv(lsf)
               if (mode == 1) then
                  call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, &
                                           error=lsf_err)
                  call batch_status(lsf, st)
               else
                  call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, &
                                           error=lsf_err, candidate_indices=idx)
                  call batch_status(lsf, st)
               end if
               md1 = svdw_max_deriv(lsf)
               if (allocated(lsf_err)) then
                  call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
                  return
               end if
               call check(error, md0 == 3 .and. md1 == md0, &
                          more="max_deriv must be restored after the on-nucleus fallback")
               if (allocated(error)) return
               call check(error, all(st == want), &
                          more="unreached exactly on the unreached reference columns")
               if (allocated(error)) return
               call check(error, all(res(:, :, n + 1:) == sentinel), &
                          more="columns past the tile must keep the sentinel")
               if (allocated(error)) return

               dev = 0.0_wp
               do i = 1, n
                  if (far(i)) then
                     if (.not. all(ieee_is_nan(res(:, :, i)))) then
                        call test_failed(error, "the unreached column must be NaN")
                        return
                     end if
                     cycle
                  end if
                  dev = max(dev, column_dev(res(:, :, i), ref(:, idx, i)))
                  do k = 1, ncen
                     if (act(idx(k), i)) cycle
                     if (any(res(:, k, i) /= 0.0_wp)) then
                        call test_failed(error, &
                                         "a candidate out of reach must hold a zero")
                        return
                     end if
                  end do
               end do
               call check(error, dev <= VJP_BATCH_REL, &
                          more="on-nucleus tile vs per-point vjp, normwise")
               if (allocated(error)) return

               ! Negative control: another point must not match this column
               miss = maxval(abs(res(:, :, c) - refm(:, idx, 1)))
               call check(error, miss >= VJP_BATCH_MISS, &
                          more="a point 1e-3 off the nucleus must miss the column")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_svdw_vjp_batch_on_nucleus_strips

   !> Empty `candidate_indices` of `vjp_f1_rA_batch`, SvdW
   subroutine test_svdw_vjp_batch_empty_list(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_empty_list(error, "svdw")
   end subroutine test_svdw_vjp_batch_empty_list

   !> Empty `candidate_indices` of `vjp_f1_rA_batch`, CFC
   subroutine test_cfc_vjp_batch_empty_list(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_empty_list(error, "cfc")
   end subroutine test_cfc_vjp_batch_empty_list

   !> An empty shared list is legal: no atom in reach, every column unreached
   !>
   !> - as `evaluate_batch`, whose list is "unique ids in 1..ncenters" (none is
   !>   vacuously fine): `res` is `(3, 0, np)`, every column unreached and no
   !>   `error`
   !> - negative control: the same points with the full list are reached
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  "svdw" or "cfc"
   subroutine run_vjp_batch_empty_list(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> "svdw" or "cfc"
      character(len=*), intent(in) :: kind

      integer, parameter :: n = 3
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), res(:, :, :), resf(:, :, :)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), w0(n), w1(ndim, n), w2(6, n), dense(ndim, ndim)
      integer, allocatable :: none(:), full(:)
      integer :: st(n), p, j
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, kind, screening_threshold=1.0e-11_wp)
      do p = 1, n
         points(:, p) = mol%xyz(:, p) + [2.0_wp, 0.0_wp, 0.0_wp]
         call vjp_seeds(p, w0(p), w1(:, p), dense, w2(:, p))
      end do
      allocate (none(0), res(ndim, 0, n))
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err, &
                               candidate_indices=none)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "an empty list must not raise: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 1), more="an empty list reaches no atom: every column unreached")
      if (allocated(error)) return

      full = [(j, j=1, lsf%ncenters)]
      allocate (resf(ndim, lsf%ncenters, n))
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resf, error=lsf_err, &
                               candidate_indices=full)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "full-list control failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0), more="the full list must reach these points")
   end subroutine run_vjp_batch_empty_list

   !> Plateau, in-band and beyond-band atoms in one strip (SvdW)
   !>
   !> - 1000-atom lattice, 32 points in a 1 Bohr cube at a corner, the fade band
   !>   set in Bohr (onset 10, width 15) and converted to e-folds through `blend_k`,
   !>   so the strip box meets all three classes ([[fade_class_counts]], asserted)
   !> - batch vs per-point `vjp_f1_rA` at [[VJP_BATCH_REL]], list absent, and again
   !>   with the onset moved to 4 Bohr, where plateau atoms of the first band
   !>   turn in-band
   !> - the fade off (hard screening) against per point as well
   !> - negative control: the second band moves the batch result of the first one
   !>   by at least [[VJP_BATCH_MISS]], so the classes bite
   subroutine test_svdw_vjp_batch_fade_classes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 32
      real(wp), parameter :: x_on(2) = [10.0_wp, 4.0_wp], delta = 15.0_wp
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp), allocatable :: radii(:)
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim), dense(ndim, ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), k
      real(wp), allocatable :: ref(:, :, :), res(:, :, :, :)
      logical, allocatable :: act(:, :)
      logical :: reached(n)
      integer :: st(n), counts(3), ib, i, p, ncen
      real(wp) :: dev, miss
      integer :: state
      type(mctc_error), allocatable :: lsf_err
      class(moist_cavity_drop_lsf_type), allocatable :: one

      call lattice_structure(mol, 10, 3.0_wp, 41)
      call get_test_radii(mol, radii)
      call param%new()
      k = param%blend_k
      ncen = mol%nat
      allocate (res(ndim, ncen, n, 2), ref(ndim, ncen, n), act(ncen, n))

      do ib = 1, 2
         call init_lsf(one, mol, radii, 3, "svdw", screening_threshold=1.0e-11_wp, &
                       fade_enabled=.true., fade_start=x_on(ib)*k/3.0_wp, &
                       fade_width=delta*k/3.0_wp)
         if (ib == 1) then
            state = 67
            do i = 1, n
               call det_rand(u, state)
               points(:, i) = one%center_lower + 1.0_wp + u
               call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
            end do
            lo = one%center_lower + 1.0_wp
            hi = lo + 1.0_wp
         end if
         call fade_class_counts(mol, radii, lo, hi, x_on(ib), delta, counts)
         call check(error, all(counts > 0), &
                    more="the strip must meet plateau, in-band and beyond atoms")
         if (allocated(error)) return

         call point_vjp_tile(one, points, w0, w1, w2, ref, act, reached)
         call check(error, all(reached), more="every corner-cube point must be reached")
         if (allocated(error)) return
         st = -1
         call one%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, :, ib), &
                                  error=lsf_err)
         call batch_status(one, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
            return
         end if
         call check(error, all(st == 0), more="every column must be reached")
         if (allocated(error)) return
         dev = 0.0_wp
         do p = 1, n
            dev = max(dev, column_dev(res(:, :, p, ib), ref(:, :, p)))
         end do
         call check(error, dev <= VJP_BATCH_REL, &
                    more="batch vs per-point vjp, fade classes")
         if (allocated(error)) return
      end do

      miss = maxval(abs(res(:, :, :, 1) - res(:, :, :, 2)))
      call check(error, miss >= VJP_BATCH_MISS, &
                 more="moving the band onset must move the batch result")
      if (allocated(error)) return

      ! Fade off: the untruncated model with hard screening at the threshold
      call init_lsf(one, mol, radii, 3, "svdw", screening_threshold=1.0e-11_wp, &
                    fade_enabled=.false.)
      call point_vjp_tile(one, points, w0, w1, w2, ref, act, reached)
      st = -1
      call one%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, :, 1), &
                               error=lsf_err)
      call batch_status(one, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "fade-off vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0) .and. all(reached), &
                 more="every fade-off column must be reached")
      if (allocated(error)) return
      dev = 0.0_wp
      do p = 1, n
         dev = max(dev, column_dev(res(:, :, p, 1), ref(:, :, p)))
      end do
      call check(error, dev <= VJP_BATCH_REL, more="batch vs per-point vjp, fade off")
   end subroutine test_svdw_vjp_batch_fade_classes

   !> Batched vs per-point nuclear vjp across three strips (CFC)
   !>
   !> - 1000-atom lattice, fade band 9 / 3 e-folds, 70 points drifting along x and
   !>   spread over 4 Bohr, so the three strips (32 + 32 + 6) exceed the 2 Bohr
   !>   sub-strip cap (asserted) and run as several sub-strips each; point 45, in the
   !>   middle of the second strip, is far out: unreached and NaN, the
   !>   columns around it computed
   !> - reference: `prepare` + `vjp_f1_rA` per point; the list absent (per-sub-strip
   !>   box lists), an ascending list and the same list reversed (both scattered
   !>   through the sub-strip permutation), normwise [[VJP_BATCH_REL]]; exact zeros at
   !>   candidates without an active slot; `res` carries sentinel columns past the
   !>   tile
   !> - tiling: the tile reversed puts every point in another strip with other
   !>   mates; the columns agree at [[VJP_BATCH_REL]], not bitwise
   !> - negative control: dropping the largest reached candidate from the list moves
   !>   the other columns by at least [[VJP_BATCH_MISS]]
   subroutine test_cfc_vjp_batch_strips(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 70, iun = 45
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), rpoints(ndim, n), u(ndim), lo(ndim), hi(ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), rw0(n), rw1(ndim, n), rw2(6, n)
      real(wp) :: dense(ndim, ndim)
      real(wp), allocatable :: ref(:, :, :), res(:, :, :), resr(:, :, :)
      real(wp), allocatable :: resv(:, :, :), resd(:, :, :)
      logical, allocatable :: act(:, :), used(:)
      logical :: reached(n)
      integer, allocatable :: list(:), rlist(:), dlist(:)
      integer :: st(n), str(n), stl(n), ncen, nsub, i, j, p, idrop
      real(wp) :: dev, devl, devt, miss, big
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call check(error, n > 2*strip_width_cfc(2) .and. n <= 3*strip_width_cfc(2), &
                 more="the tile must cross exactly three strips")
      if (allocated(error)) return
      call lattice_structure(mol, 10, 3.0_wp, 31)
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, "cfc", screening_threshold=1.0e-11_wp, &
                    fade_enabled=.true., fade_start=9.0_wp, fade_width=3.0_wp)
      ncen = lsf%ncenters

      state = 4711
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = lsf%center_lower + 4.0_wp*u &
                        + [1.0_wp + 0.12_wp*real(i, wp), 1.0_wp, 1.0_wp]
         call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
      end do
      points(:, iun) = lsf%center_upper + 200.0_wp
      lo = minval(points(:, 1:strip_width_cfc(2)), dim=2)
      hi = maxval(points(:, 1:strip_width_cfc(2)), dim=2)
      call check(error, norm2(hi - lo) > 2.0_wp, &
                 more="the first strip must exceed the sub-strip diagonal")
      if (allocated(error)) return
      do i = 1, n
         rpoints(:, i) = points(:, n + 1 - i)
         rw0(i) = w0(n + 1 - i)
         rw1(:, i) = w1(:, n + 1 - i)
         rw2(:, i) = w2(:, n + 1 - i)
      end do

      allocate (ref(ndim, ncen, n), act(ncen, n), res(ndim, ncen, n + 2), &
                resr(ndim, ncen, n), used(ncen))
      call point_vjp_tile(lsf, points, w0, w1, w2, ref, act, reached)
      call check(error, count(.not. reached) == 1 .and. .not. reached(iun), &
                 more="only the far column may be unreached per point")
      if (allocated(error)) return

      ! List absent: per-sub-strip box lists
      res = sentinel
      st = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, 1:n), &
                               error=lsf_err)
      call batch_status(lsf, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == merge(1, 0, [(i == iun, i=1, n)])), &
                 more="unreached on the far column only")
      if (allocated(error)) return
      if (.not. all(ieee_is_nan(res(:, :, iun)))) then
         call test_failed(error, "the unreached column must be NaN")
         return
      end if
      if (any(res(:, :, n + 1:n + 2) /= sentinel)) then
         call test_failed(error, "vjp_f1_rA_batch wrote past the last tile column")
         return
      end if
      dev = 0.0_wp
      do p = 1, n
         if (p == iun) cycle
         dev = max(dev, column_dev(res(:, :, p), ref(:, :, p)))
         do j = 1, ncen
            if (.not. act(j, p) .and. any(res(:, j, p) /= 0.0_wp)) then
               call test_failed(error, "a candidate out of reach must hold an exact zero")
               return
            end if
         end do
      end do
      call check(error, dev <= VJP_BATCH_REL, &
                 more="batch vs per-point vjp over three strips")
      if (allocated(error)) return

      ! Ascending strict list and its reverse
      do j = 1, ncen
         used(j) = any(act(j, :))
      end do
      list = pack([(j, j=1, ncen)], used .or. mod([(j, j=1, ncen)], 7) == 0)
      nsub = size(list)
      call check(error, count(used) > 0 .and. nsub > count(used) .and. nsub < ncen, &
                 more="the list must be a strict superset of the reached candidates")
      if (allocated(error)) return
      allocate (resv(ndim, nsub, n + 2), resd(ndim, nsub, n + 2))
      resv = sentinel
      resd = sentinel
      stl = -1
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resv, error=lsf_err, &
                               candidate_indices=list)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "list vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(stl == st), more="the list must not move the unreached flags")
      if (allocated(error)) return
      rlist = list(nsub:1:-1)
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resd, error=lsf_err, &
                               candidate_indices=rlist)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "reversed vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      devl = 0.0_wp
      do p = 1, n
         if (p == iun) cycle
         devl = max(devl, column_dev(resv(:, :, p), ref(:, list, p)))
         do j = 1, nsub
            if (any(resd(:, j, p) /= resv(:, nsub + 1 - j, p))) then
               call test_failed(error, "reversed list changes the bits of a column")
               return
            end if
         end do
      end do
      call check(error, devl <= VJP_BATCH_REL, more="list vs per-point vjp")
      if (allocated(error)) return
      if (any(resv(:, :, n + 1:n + 2) /= sentinel) .or. &
          any(resd(:, :, n + 1:n + 2) /= sentinel)) then
         call test_failed(error, "vjp_f1_rA_batch wrote past the last tile column")
         return
      end if

      ! Tiling: the reversed tile, every point in another strip
      str = -1
      call lsf%vjp_f1_rA_batch(rpoints, rw0, rw1, rw2, resr, error=lsf_err)
      call batch_status(lsf, str)
      if (allocated(lsf_err)) then
         call test_failed(error, "reversed-tile vjp_f1_rA_batch failed: "// &
                          lsf_err%message)
         return
      end if
      call check(error, all(str(n:1:-1) == st), &
                 more="the tile order must not move the unreached flags")
      if (allocated(error)) return
      devt = 0.0_wp
      do p = 1, n
         if (p == iun) cycle
         devt = max(devt, column_dev(res(:, :, p), resr(:, :, n + 1 - p)))
      end do
      call check(error, devt <= VJP_BATCH_REL, &
                 more="the same point in two different tiles")
      if (allocated(error)) return

      ! Negative control: drop the reached candidate of the largest column
      big = -1.0_wp
      idrop = 0
      do j = 1, nsub
         if (.not. used(list(j))) cycle
         if (maxval(abs(resv(:, j, 1:n))) > big) then
            big = maxval(abs(resv(:, j, 1:n)))
            idrop = list(j)
         end if
      end do
      dlist = pack(list, list /= idrop)
      deallocate (resd)
      allocate (resd(ndim, size(dlist), n))
      call lsf%vjp_f1_rA_batch(points, w0, w1, w2, resd, error=lsf_err, &
                               candidate_indices=dlist)
      call batch_status(lsf, stl)
      if (allocated(lsf_err)) then
         call test_failed(error, "dropped-candidate batch failed: "//lsf_err%message)
         return
      end if
      miss = 0.0_wp
      do p = 1, n
         if (stl(p) /= 0) cycle
         do j = 1, size(dlist)
            miss = max(miss, maxval(abs(resd(:, j, p) - ref(:, dlist(j), p))))
         end do
      end do
      call check(error, miss >= VJP_BATCH_MISS, &
                 more="dropping a reached candidate must move the other columns")
   end subroutine test_cfc_vjp_batch_strips

   !> Maximum derivative order the CFC LSF is prepared for, -1 for another LSF
   !>
   !> @param[in] lsf LSF instance
   !> @returns   md  `max_deriv` of a CFC LSF
   function cfc_max_deriv(lsf) result(md)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Prepared derivative order
      integer :: md

      md = -1
      select type (lsf)
      type is (moist_cavity_drop_lsf_cfc_type)
         md = lsf%max_deriv
      end select
   end function cfc_max_deriv

   !> Per-point cache state of the CFC LSF, -2 for another LSF
   !>
   !> @param[in] lsf LSF instance
   !> @returns   pd  `prepared_deriv` of a CFC LSF, -1 when the cache is invalid
   function cfc_prepared_deriv(lsf) result(pd)
      !> LSF instance
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Order the per-point cache holds
      integer :: pd

      pd = -2
      select type (lsf)
      type is (moist_cavity_drop_lsf_cfc_type)
         pd = lsf%prepared_deriv
      end select
   end function cfc_prepared_deriv

   !> On-nucleus columns across strips of one vjp tile (CFC)
   !>
   !> - MB16-43 01 and a 1000-atom lattice (inactive candidates on the displaced
   !>   path), 70 points in three strips (32 + 32 + 6) near atom 1, fade band 9 / 3
   !>   e-folds, `max_deriv = 4`
   !> - the cusp of a nucleus makes the jet non-finite from order 1: the per-point
   !>   reference refuses an on-nucleus point with an `error`, and the batch raises
   !>   an `error` naming an on-nucleus column, never an unreached one
   !> - case 1: on-nucleus columns 5 and 40, in two different strips
   !> - case 2: on-nucleus column 5 next to an unreached column 6 in one strip
   !> - case 3: every column of the second strip on an atom
   !> - each case with the list absent, the full ascending list and the full list
   !>   reversed (scatter); `max_deriv` unchanged and the per-point cache invalid
   !>   after the refused call
   !> - the tile without its on-nucleus columns against the per-point `prepare` +
   !>   `vjp_f1_rA`, normwise [[VJP_BATCH_REL]]: sentinel past the tile, exact zeros
   !>   at inactive candidates, the far column unreached and NaN
   !> - negative control: the tile with every on-nucleus point moved 1e-3 off the
   !>   atom raises nothing, is reached off the far column, and equals the displaced
   !>   per-point columns
   subroutine test_cfc_vjp_batch_on_nucleus_strips(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 70, ncase = 3, nextra = 2
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: base(ndim, n), points(ndim, n), pmove(ndim, n), u(ndim)
      real(wp) :: dense(ndim, ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n)
      real(wp), allocatable :: ref(:, :, :), refm(:, :, :), res(:, :, :), resm(:, :, :)
      logical, allocatable :: act(:, :), actm(:, :)
      logical :: reached(n), reachedm(n), far(n), named
      integer, allocatable :: idx(:), keep(:)
      integer :: atom_of(n), st(n), stm(n), ifix, icase, mode, i, j, k, c, ncen, m
      integer :: md0, md1
      real(wp) :: dev, devm
      character(len=16) :: col
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      do ifix = 1, 2
         if (allocated(ref)) deallocate (ref, refm, res, resm, act, actm)
         if (ifix == 1) then
            call get_structure(mol, "MB16-43", "01")
         else
            call lattice_structure(mol, 10, 3.0_wp, 31)
         end if
         call get_test_radii(mol, radii)
         call init_lsf(lsf, mol, radii, 4, "cfc", screening_threshold=1.0e-11_wp, &
                       fade_enabled=.true., fade_start=9.0_wp, fade_width=3.0_wp)
         ncen = lsf%ncenters
         state = 29
         do i = 1, n
            call det_rand(u, state)
            base(:, i) = mol%xyz(:, 1) + 3.0_wp*(u - 0.5_wp)
            call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
         end do
         allocate (ref(ndim, ncen, n), refm(ndim, ncen, n), res(ndim, ncen, n + nextra), &
                   resm(ndim, ncen, n + nextra), act(ncen, n), actm(ncen, n))

         do icase = 1, ncase
            atom_of = 0
            far = .false.
            select case (icase)
            case (1)
               atom_of(5) = 3
               atom_of(40) = 7
            case (2)
               atom_of(5) = 3
               far(6) = .true.
            case default
               atom_of(33:64) = [(1 + mod(j - 1, mol%nat), j=1, 32)]
            end select
            points = base
            do i = 1, n
               if (atom_of(i) > 0) points(:, i) = mol%xyz(:, atom_of(i))
               if (far(i)) points(:, i) = lsf%center_upper + 200.0_wp
            end do
            c = findloc(atom_of > 0, .true., dim=1)
            keep = pack([(i, i=1, n)], atom_of == 0)
            m = size(keep)

            call point_vjp_tile(lsf, points, w0, w1, w2, ref, act, reached)
            call check(error, all(reached .eqv. (atom_of == 0 .and. .not. far)), &
                       more="the per-point reference must refuse the nucleus columns "// &
                       "and leave the far ones unreached")
            if (allocated(error)) return
            ! Negative control reference: every on-nucleus point displaced
            pmove = points
            do i = 1, n
               if (atom_of(i) > 0) then
                  pmove(:, i) = points(:, i) + [1.0e-3_wp, 0.0_wp, 0.0_wp]
               end if
            end do
            call point_vjp_tile(lsf, pmove, w0, w1, w2, refm, actm, reachedm)
            call check(error, all(reachedm .eqv. .not. far), &
                       more="the displaced points must be reached")
            if (allocated(error)) return
            ! The lattice must leave inactive candidates on the displaced path
            call check(error, ifix == 1 .or. .not. all(actm(:, c)), &
                       more="a displaced lattice column must have inactive candidates")
            if (allocated(error)) return

            do mode = 1, 3
               if (mode == 3) then
                  idx = [(ncen + 1 - j, j=1, ncen)]
               else
                  idx = [(j, j=1, ncen)]
               end if
               md0 = cfc_max_deriv(lsf)
               if (mode == 1) then
                  call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
               else
                  call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err, &
                                           candidate_indices=idx)
               end if
               md1 = cfc_max_deriv(lsf)
               if (.not. allocated(lsf_err)) then
                  call test_failed(error, "an on-nucleus column must raise")
                  return
               end if
               named = .false.
               do i = 1, n
                  if (atom_of(i) == 0) cycle
                  write (col, "(i0)") i
                  named = named .or. index(lsf_err%message, "vjp_f1_rA_batch: point "// &
                                           trim(col)//" lies on a nucleus") > 0
               end do
               call check(error, named, more="error must name an on-nucleus column: "// &
                          lsf_err%message)
               if (allocated(error)) return
               deallocate (lsf_err)
               call check(error, md0 == 4 .and. md1 == md0, &
                          more="max_deriv must be unchanged by the refused call")
               if (allocated(error)) return
               call check(error, cfc_prepared_deriv(lsf) == -1, &
                          more="the per-point cache must be invalid after the batch")
               if (allocated(error)) return

               ! The tile without its on-nucleus columns
               res = sentinel
               st = -1
               if (mode == 1) then
                  call lsf%vjp_f1_rA_batch(points(:, keep), w0(keep), w1(:, keep), &
                                           w2(:, keep), res, error=lsf_err)
               else
                  call lsf%vjp_f1_rA_batch(points(:, keep), w0(keep), w1(:, keep), &
                                           w2(:, keep), res, error=lsf_err, &
                                           candidate_indices=idx)
               end if
               call batch_status(lsf, st(1:m))
               if (allocated(lsf_err)) then
                  call test_failed(error, "the tile without the nuclei failed: "// &
                                   lsf_err%message)
                  return
               end if
               call check(error, all(st(1:m) == merge(1, 0, far(keep))), &
                          more="unreached exactly on the far columns")
               if (allocated(error)) return
               call check(error, all(res(:, :, m + 1:) == sentinel), &
                          more="columns past the tile must keep the sentinel")
               if (allocated(error)) return

               dev = 0.0_wp
               do k = 1, m
                  i = keep(k)
                  if (.not. reached(i)) then
                     if (.not. all(ieee_is_nan(res(:, :, k)))) then
                        call test_failed(error, "an unreached column must be NaN")
                        return
                     end if
                     cycle
                  end if
                  dev = max(dev, column_dev(res(:, :, k), ref(:, idx, i)))
                  do j = 1, ncen
                     if (act(idx(j), i)) cycle
                     if (any(res(:, j, k) /= 0.0_wp)) then
                        call test_failed(error, &
                                         "a candidate out of reach must hold a zero")
                        return
                     end if
                  end do
               end do
               call check(error, dev <= VJP_BATCH_REL, &
                          more="tile without the nuclei vs per-point vjp, normwise")
               if (allocated(error)) return

               ! Negative control: the same tile with the nuclei 1e-3 off raises nothing,
               ! the displaced columns equal the displaced reference
               resm = sentinel
               stm = -1
               if (mode == 1) then
                  call lsf%vjp_f1_rA_batch(pmove, w0, w1, w2, resm, error=lsf_err)
               else
                  call lsf%vjp_f1_rA_batch(pmove, w0, w1, w2, resm, error=lsf_err, &
                                           candidate_indices=idx)
               end if
               call batch_status(lsf, stm)
               if (allocated(lsf_err)) then
                  call test_failed(error, "displaced tile failed: "//lsf_err%message)
                  return
               end if
               call check(error, all(stm(:) == merge(1, 0, far)), &
                          more="a point 1e-3 off the nucleus must be reached")
               if (allocated(error)) return
               devm = 0.0_wp
               do i = 1, n
                  if (far(i)) cycle
                  devm = max(devm, column_dev(resm(:, :, i), refm(:, idx, i)))
               end do
               call check(error, devm <= VJP_BATCH_REL, &
                          more="displaced tile vs per-point vjp, normwise")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_cfc_vjp_batch_on_nucleus_strips

   !> A nucleus inside a compact vjp strip is refused, its mates alone are not (CFC)
   !>
   !> - MB16-43 01, twelve points within 0.4 Bohr of atom 3, one of them exactly on
   !>   it, so the strip stays one sub-strip (box diagonal <= 2 Bohr)
   !> - the tile with the nucleus column: `error` naming column 5, list absent and
   !>   full list
   !> - the other eleven columns alone: reached, nothing written past the tile
   !> - negative control: a mate with a changed seed misses its unchanged twin
   !>
   !> @param[out] error Test failure
   subroutine test_cfc_vjp_batch_on_nucleus_mates(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 12, nuc = 5
      real(wp), parameter :: sentinel = 12345.0_wp
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), pn(ndim, n - 1), u(ndim), dense(ndim, ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), v0(n - 1), v1(ndim, n - 1)
      real(wp) :: v2(6, n - 1), x0(n - 1)
      real(wp), allocatable :: res(:, :, :), resn(:, :, :), resm(:, :, :)
      integer, allocatable :: full(:)
      integer :: stn(n - 1), mode, i, k, ncen
      real(wp) :: mate_miss
      integer :: state
      type(mctc_error), allocatable :: lsf_err

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, "cfc", screening_threshold=1.0e-11_wp)
      ncen = lsf%ncenters
      state = 17
      do i = 1, n
         call det_rand(u, state)
         points(:, i) = mol%xyz(:, 3) + 0.8_wp*(u - 0.5_wp)
         call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
      end do
      points(:, nuc) = mol%xyz(:, 3)
      k = 0
      do i = 1, n
         if (i == nuc) cycle
         k = k + 1
         pn(:, k) = points(:, i)
         v0(k) = w0(i)
         v1(:, k) = w1(:, i)
         v2(:, k) = w2(:, i)
      end do
      x0 = v0
      x0(2) = x0(2) + 1.0e-3_wp
      full = [(i, i=1, ncen)]
      allocate (res(ndim, ncen, n), resn(ndim, ncen, n), resm(ndim, ncen, n - 1))

      do mode = 1, 2
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
         else
            call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err, &
                                     candidate_indices=full)
         end if
         if (.not. allocated(lsf_err)) then
            call test_failed(error, "an on-nucleus column must raise")
            return
         end if
         call check(error, index(lsf_err%message, &
                                 "vjp_f1_rA_batch: point 5 lies on a nucleus") > 0, &
                    more="error must name the on-nucleus column: "//lsf_err%message)
         if (allocated(error)) return
         deallocate (lsf_err)

         resn = sentinel
         stn = -1
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(pn, v0, v1, v2, resn, error=lsf_err)
         else
            call lsf%vjp_f1_rA_batch(pn, v0, v1, v2, resn, error=lsf_err, &
                                     candidate_indices=full)
         end if
         call batch_status(lsf, stn)
         if (allocated(lsf_err)) then
            call test_failed(error, "tile without the nucleus failed: "//lsf_err%message)
            return
         end if
         call check(error, all(stn == 0), more="the mates alone must be reached")
         if (allocated(error)) return
         call check(error, all(resn(:, :, n) == sentinel), &
                    more="the column past the tile must keep the sentinel")
         if (allocated(error)) return

         ! Negative control: a changed seed moves its own column
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(pn, x0, v1, v2, resm, error=lsf_err)
         else
            call lsf%vjp_f1_rA_batch(pn, x0, v1, v2, resm, error=lsf_err, &
                                     candidate_indices=full)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "changed-seed tile failed: "//lsf_err%message)
            return
         end if
         mate_miss = maxval(abs(resm(:, :, 2) - resn(:, :, 2)))
         call check(error, mate_miss >= VJP_BATCH_MISS, &
                    more="a changed seed must move its own column")
         if (allocated(error)) return
      end do
   end subroutine test_cfc_vjp_batch_on_nucleus_mates

   !> Atoms and atom pairs of the plateau, in-band and beyond classes (CFC)
   !>
   !> The classes of stage 1's strip against the box `[lo, hi]` in the reduced
   !> distance `sigma = (|r - R_A| - R_a)/R_a` at the nearest and farthest box point:
   !> plateau `sigma_far <= x_on`, beyond `sigma_near >= x_on + delta`, else in band;
   !> pairs are counted among the atoms that are not beyond
   !>
   !> @param[in]  mol     Structure
   !> @param[in]  radii   Per-atom radii
   !> @param[in]  lo      Lower corner of the box
   !> @param[in]  hi      Upper corner of the box
   !> @param[in]  x_on    Fade onset in sigma
   !> @param[in]  delta   Fade width in sigma
   !> @param[out] atoms   Atoms in the plateau, band and beyond classes (3)
   !> @param[out] pairs   Pairs plateau-plateau, plateau-band, band-band (3)
   pure subroutine cfc_class_counts(mol, radii, lo, hi, x_on, delta, atoms, pairs)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Lower corner of the box
      real(wp), intent(in) :: lo(ndim)
      !> Upper corner of the box
      real(wp), intent(in) :: hi(ndim)
      !> Fade onset in sigma
      real(wp), intent(in) :: x_on
      !> Fade width in sigma
      real(wp), intent(in) :: delta
      !> Atoms per class
      integer, intent(out) :: atoms(3)
      !> Pairs per class combination
      integer, intent(out) :: pairs(3)

      real(wp) :: e(ndim), f(ndim), snear, sfar
      integer :: cls(mol%nat), a

      atoms = 0
      do a = 1, mol%nat
         e = max(lo - mol%xyz(:, a), 0.0_wp, mol%xyz(:, a) - hi)
         f = max(abs(mol%xyz(:, a) - lo), abs(mol%xyz(:, a) - hi))
         snear = (norm2(e) - radii(a))/radii(a)
         sfar = (norm2(f) - radii(a))/radii(a)
         if (sfar <= x_on) then
            cls(a) = 1
         else if (snear >= x_on + delta) then
            cls(a) = 3
         else
            cls(a) = 2
         end if
         atoms(cls(a)) = atoms(cls(a)) + 1
      end do
      pairs(1) = atoms(1)*(atoms(1) - 1)/2
      pairs(2) = atoms(1)*atoms(2)
      pairs(3) = atoms(2)*(atoms(2) - 1)/2
   end subroutine cfc_class_counts

   !> Plateau, in-band and beyond-band atoms and pairs in one strip (CFC)
   !>
   !> - 1000-atom lattice, 32 points in a 1 Bohr cube at a corner, the band set in
   !>   e-folds 12 / 3 and 8 / 3; the strip box meets plateau, in-band and beyond
   !>   atoms and the plateau-plateau, plateau-band and band-band pairs
   !>   ([[cfc_class_counts]], asserted by counting)
   !> - batch vs per-point `vjp_f1_rA` at [[VJP_BATCH_REL]], list absent, for both
   !>   bands, and with the fade off (hard screening) as well
   !> - negative control: the second band moves the batch result of the first one
   !>   by at least [[VJP_BATCH_MISS]], so the classes bite
   subroutine test_cfc_vjp_batch_fade_classes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 32
      real(wp), parameter :: efold_on(2) = [12.0_wp, 8.0_wp], efold_w = 3.0_wp
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_param_type) :: param
      real(wp), allocatable :: radii(:)
      real(wp) :: points(ndim, n), u(ndim), lo(ndim), hi(ndim), dense(ndim, ndim)
      real(wp) :: w0(n), w1(ndim, n), w2(6, n), a2
      real(wp), allocatable :: ref(:, :, :), res(:, :, :, :)
      logical, allocatable :: act(:, :)
      logical :: reached(n)
      integer :: st(n), atoms(3), pairs(3), ib, i, p, ncen
      real(wp) :: dev, miss
      integer :: state
      type(mctc_error), allocatable :: lsf_err
      class(moist_cavity_drop_lsf_type), allocatable :: one

      call lattice_structure(mol, 10, 3.0_wp, 41)
      call get_test_radii(mol, radii)
      call param%new()
      a2 = abs(param%a2)
      ncen = mol%nat
      allocate (res(ndim, ncen, n, 2), ref(ndim, ncen, n), act(ncen, n))

      do ib = 1, 2
         call init_lsf(one, mol, radii, 3, "cfc", screening_threshold=1.0e-11_wp, &
                       fade_enabled=.true., fade_start=efold_on(ib), fade_width=efold_w)
         if (ib == 1) then
            state = 67
            do i = 1, n
               call det_rand(u, state)
               points(:, i) = one%center_lower + 1.0_wp + u
               call vjp_seeds(i, w0(i), w1(:, i), dense, w2(:, i))
            end do
            lo = one%center_lower + 1.0_wp
            hi = lo + 1.0_wp
         end if
         call cfc_class_counts(mol, radii, lo, hi, efold_on(ib)/a2, efold_w/a2, atoms, &
                               pairs)
         call check(error, all(atoms > 0) .and. all(pairs > 0), &
                    more="the strip must meet every atom and pair class")
         if (allocated(error)) return

         call point_vjp_tile(one, points, w0, w1, w2, ref, act, reached)
         call check(error, all(reached), more="every corner-cube point must be reached")
         if (allocated(error)) return
         st = -1
         call one%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, :, ib), &
                                  error=lsf_err)
         call batch_status(one, st)
         if (allocated(lsf_err)) then
            call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
            return
         end if
         call check(error, all(st == 0), more="every column must be reached")
         if (allocated(error)) return
         dev = 0.0_wp
         do p = 1, n
            dev = max(dev, column_dev(res(:, :, p, ib), ref(:, :, p)))
         end do
         call check(error, dev <= VJP_BATCH_REL, &
                    more="batch vs per-point vjp, fade classes")
         if (allocated(error)) return
      end do

      miss = maxval(abs(res(:, :, :, 1) - res(:, :, :, 2)))
      call check(error, miss >= VJP_BATCH_MISS, &
                 more="moving the band onset must move the batch result")
      if (allocated(error)) return

      ! Fade off: the untruncated model with hard screening at the threshold
      call init_lsf(one, mol, radii, 3, "cfc", screening_threshold=1.0e-11_wp, &
                    fade_enabled=.false.)
      call point_vjp_tile(one, points, w0, w1, w2, ref, act, reached)
      st = -1
      call one%vjp_f1_rA_batch(points, w0, w1, w2, res(:, :, :, 1), &
                               error=lsf_err)
      call batch_status(one, st)
      if (allocated(lsf_err)) then
         call test_failed(error, "fade-off vjp_f1_rA_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(st == 0) .and. all(reached), &
                 more="every fade-off column must be reached")
      if (allocated(error)) return
      dev = 0.0_wp
      do p = 1, n
         dev = max(dev, column_dev(res(:, :, p, 1), ref(:, :, p)))
      end do
      call check(error, dev <= VJP_BATCH_REL, more="batch vs per-point vjp, fade off")
   end subroutine test_cfc_vjp_batch_fade_classes

   !> Point on an outward ray whose nearest atom sits a given number of e-folds out
   !>
   !> - e-folds of atom A at `x`: `rate (|x - R_A|/R_a - 1)` for `scaled`, else
   !>   `rate (|x - R_A| - R_a)`; bisection on the ray `centroid + t direction`
   !>
   !> @param[in]  mol       Structure
   !> @param[in]  radii     Per-atom radii
   !> @param[in]  rate      e-folds per unit of the distance measure
   !> @param[in]  scaled    Distance in units of the atom's own radius
   !> @param[in]  direction Ray direction, any length
   !> @param[in]  efold     e-folds of the nearest atom at the point
   !> @param[out] point     Point on the ray
   subroutine point_at_efolds(mol, radii, rate, scaled, direction, efold, point)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> e-folds per unit of the distance measure
      real(wp), intent(in) :: rate
      !> Distance in units of the atom's own radius
      logical, intent(in) :: scaled
      !> Ray direction
      real(wp), intent(in) :: direction(ndim)
      !> Target e-folds of the nearest atom
      real(wp), intent(in) :: efold
      !> Point on the ray
      real(wp), intent(out) :: point(ndim)

      real(wp) :: centroid(ndim), dir(ndim), x(ndim), t_lo, t_hi, t, d, e
      integer :: it, a

      centroid = sum(mol%xyz, dim=2)/real(mol%nat, wp)
      dir = direction/norm2(direction)
      t_lo = 0.0_wp
      t_hi = 1.0e4_wp
      do it = 1, 200
         t = 0.5_wp*(t_lo + t_hi)
         x = centroid + t*dir
         e = huge(1.0_wp)
         do a = 1, mol%nat
            d = norm2(x - mol%xyz(:, a))
            if (scaled) then
               e = min(e, rate*(d/radii(a) - 1.0_wp))
            else
               e = min(e, rate*(d - radii(a)))
            end if
         end do
         if (e < efold) then
            t_lo = t
         else
            t_hi = t
         end if
      end do
      point = centroid + t_hi*dir
   end subroutine point_at_efolds

   !> Tiny level-set sums where the adjoint overflows and the lift does not (SvdW)
   subroutine test_svdw_vjp_batch_tiny_pd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_tiny_pd(error, "svdw")
   end subroutine test_svdw_vjp_batch_tiny_pd

   !> Tiny level-set sums where the adjoint overflows and the lift does not (CFC)
   subroutine test_cfc_vjp_batch_tiny_pd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call run_vjp_batch_tiny_pd(error, "cfc")
   end subroutine test_cfc_vjp_batch_tiny_pd

   !> A reached column whose adjoint overflows is unreached, never reached with Inf
   !>
   !> - the lift divides by the level-set sum `PD` (SvdW: `Z`) through `1/PD^2`, the
   !>   adjoint through `1/PD^3`: for `PD` between about 1e-154 and 1e-103 the order-2
   !>   lift is finite and the adjoint is not
   !> - fade band moved far out (240 / 3 e-folds), so a point 100s of e-folds from
   !>   every atom is still inside the band; MB16-43 01, one outward ray per column,
   !>   the distance found from the e-fold measure of the nearest atom
   !>   ([[point_at_efolds]]); the tiny column sits at 165 e-folds (CFC, `PD` ~
   !>   1e-120) and 100 e-folds (SvdW, `Z ~ u^3` ~ 1e-130), measured window 150-210
   !>   and 80-110
   !> - tile of three columns: a control at `PD` ~ 1e-30, the tiny column, a second
   !>   control; the precondition is asserted: `evaluate_batch` at order 2 reaches
   !>   the tiny column with a finite jet
   !> - tiny column: unreached (`S > 0`, no `error`) and NaN, list absent and reversed;
   !>   alone in a tile with the ascending list (the kernel writes in place) the same;
   !>   it is unreached for the per-point path too
   !> - controls: reached, equal to the per-point `vjp_f1_rA` at [[VJP_BATCH_REL]],
   !>   the column past the tile keeps its sentinel
   !>
   !> @param[out] error Test failure
   !> @param[in]  kind  "svdw" or "cfc"
   subroutine run_vjp_batch_tiny_pd(error, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> "svdw" or "cfc"
      character(len=*), intent(in) :: kind

      integer, parameter :: n = 3
      real(wp), parameter :: sentinel = 12345.0_wp
      real(wp), parameter :: dirs(ndim, n) = reshape([1.0_wp, 0.3_wp, -0.5_wp, &
                                                      -0.4_wp, 1.0_wp, 0.2_wp, &
                                                      0.2_wp, -0.7_wp, 1.0_wp], [ndim, n])
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_param_type) :: sparam
      type(moist_cavity_drop_lsf_cfc_param_type) :: cparam
      real(wp), allocatable :: radii(:), res(:, :, :), ref(:, :, :), one(:, :, :)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: points(ndim, n), w0(n), w1(ndim, n), w2(6, n), dense(ndim, ndim)
      real(wp) :: f0(n), f1(ndim, n), f2(6, n)
      real(wp) :: rate, fade_start, fade_width, e_tiny, e_ctl, dev
      logical, allocatable :: act(:, :)
      logical :: scaled, reached(n)
      integer, allocatable :: full(:), rev(:)
      integer :: st(n), st1(1), so(n), mode, p, ncen
      type(mctc_error), allocatable :: lsf_err

      select case (kind)
      case ("svdw")
         call sparam%new()
         rate = sparam%blend_k/3.0_wp
         scaled = .false.
         fade_start = 240.0_wp
         fade_width = 3.0_wp
         e_tiny = 100.0_wp
         e_ctl = 23.0_wp
      case default
         call cparam%new()
         rate = abs(cparam%a2)
         scaled = .true.
         fade_start = 240.0_wp
         fade_width = 3.0_wp
         e_tiny = 165.0_wp
         e_ctl = 41.0_wp
      end select
      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 3, kind, screening_threshold=1.0e-11_wp, &
                    fade_enabled=.true., fade_start=fade_start, fade_width=fade_width)
      ncen = lsf%ncenters
      call point_at_efolds(mol, radii, rate, scaled, dirs(:, 1), e_ctl, points(:, 1))
      call point_at_efolds(mol, radii, rate, scaled, dirs(:, 2), e_tiny, points(:, 2))
      call point_at_efolds(mol, radii, rate, scaled, dirs(:, 3), 1.1_wp*e_ctl, &
                           points(:, 3))
      do p = 1, n
         call vjp_seeds(p, w0(p), w1(:, p), dense, w2(:, p))
      end do

      ! Precondition: the order-2 lift reaches the tiny column with a finite jet
      call lsf%evaluate_batch(points, 2, f0, f1, f2, error=lsf_err)
      call batch_status(lsf, so)
      if (allocated(lsf_err)) then
         call test_failed(error, "order-2 evaluate_batch failed: "//lsf_err%message)
         return
      end if
      call check(error, all(so == 0), &
                 more="the tiny column must be reached by the order-2 lift")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(f0)) .and. all(ieee_is_finite(f1)) .and. &
                 all(ieee_is_finite(f2)), more="the order-2 jet must be finite")
      if (allocated(error)) return

      allocate (ref(ndim, ncen, n), act(ncen, n), res(ndim, ncen, n + 1), &
                one(ndim, ncen, 1))
      call point_vjp_tile(lsf, points, w0, w1, w2, ref, act, reached)
      call check(error, reached(1) .and. reached(3) .and. .not. reached(2), &
                 more="the per-point path must reach the controls, not the tiny column")
      if (allocated(error)) return

      full = [(p, p=1, ncen)]
      rev = [(ncen + 1 - p, p=1, ncen)]
      do mode = 1, 2
         res = sentinel
         st = -1
         if (mode == 1) then
            call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, error=lsf_err)
            call batch_status(lsf, st)
         else
            call lsf%vjp_f1_rA_batch(points, w0, w1, w2, res, &
                                     error=lsf_err, candidate_indices=rev)
            call batch_status(lsf, st)
         end if
         if (allocated(lsf_err)) then
            call test_failed(error, "vjp_f1_rA_batch failed: "//lsf_err%message)
            return
         end if
         call check(error, all(st == [0, 1, 0]), &
                    more="tiny column must be unreached, controls reached")
         if (allocated(error)) return
         call check(error, all(ieee_is_nan(res(:, :, 2))), &
                    more="the tiny column must be NaN, not Inf or a number")
         if (allocated(error)) return
         call check(error, all(res(:, :, n + 1) == sentinel), &
                    more="the column past the tile must keep the sentinel")
         if (allocated(error)) return
         dev = 0.0_wp
         do p = 1, n, 2
            if (mode == 1) then
               dev = max(dev, column_dev(res(:, :, p), ref(:, :, p)))
            else
               dev = max(dev, column_dev(res(:, :, p), ref(:, rev, p)))
            end if
         end do
         call check(error, dev <= VJP_BATCH_REL, &
                    more="control columns vs per-point vjp, normwise")
         if (allocated(error)) return
      end do

      ! Alone with the ascending list: the kernel writes the tile column in place
      one = sentinel
      st1 = -1
      call lsf%vjp_f1_rA_batch(points(:, 2:2), w0(2:2), w1(:, 2:2), w2(:, 2:2), one, &
                               error=lsf_err, candidate_indices=full)
      call batch_status(lsf, st1)
      if (allocated(lsf_err)) then
         call test_failed(error, "single-column vjp_f1_rA_batch failed: "// &
                          lsf_err%message)
         return
      end if
      call check(error, st1(1) == 1 .and. all(ieee_is_nan(one)), &
                 more="the tiny column alone, in place: unreached and NaN")
   end subroutine run_vjp_batch_tiny_pd

   !> Jittered cubic lattice of H and C atoms, as in the box-candidate test
   !>
   !> @param[out] mol     Structure, `nside**3` atoms
   !> @param[in]  nside   Atoms per edge
   !> @param[in]  spacing Lattice constant (Bohr)
   !> @param[in]  seed    LCG seed of the jitter
   subroutine lattice_structure(mol, nside, spacing, seed)
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Atoms per edge
      integer, intent(in) :: nside
      !> Lattice constant
      real(wp), intent(in) :: spacing
      !> LCG seed
      integer, intent(in) :: seed

      integer, allocatable :: num(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: u(ndim)
      integer :: state
      integer :: i, j, k, n

      allocate (num(nside**3), xyz(ndim, nside**3))
      state = seed
      n = 0
      do k = 1, nside
         do j = 1, nside
            do i = 1, nside
               n = n + 1
               call det_rand(u, state)
               xyz(:, n) = spacing*(real([i, j, k], wp) + 0.4_wp*(u - 0.5_wp))
               num(n) = merge(6, 1, mod(n, 3) == 0)
            end do
         end do
      end do
      call new(mol, num, xyz)
   end subroutine lattice_structure

end module test_cavity_drop_lsf
