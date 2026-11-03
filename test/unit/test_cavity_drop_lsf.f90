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
!>   * `exclusion_radius_never_overclaims`     surface-free certified balls
!>   * `exclusion_radius_gated_on_*`           no certificate outside the bound's
!>                                             assumptions
!>
!> Empty active lists (SvdW, CFC):
!>   * `empty_active`                          full result buffers with all atoms screened
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
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, &
                           center_at_origin, fd4_scalar
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_sses, only: moist_cavity_drop_lsf_sses_type
   use moist_cavity_drop_lsf_sses_algebra, only: sses_jet_type, jet_derivative
   use moist_cavity_drop_lsf_sses_kernel, only: seeded_evaluate => sses_evaluate
   use moist_cavity_drop_lsf_sses_cache, only: sses_cache_type, sses_cache_build, sses_cached_evaluate
   use moist_cavity_drop_lsf_sses_param, only: moist_cavity_drop_lsf_sses_param_type, temperature_length
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_sigmoid_bump_type, new_swif_sigmoid_bump
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

   !> Result-buffer poison for empty active lists
   real(wp), parameter :: EMPTY_POISON = -1.0e30_wp

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
                  new_unittest("cfc_exclusion_radius_never_overclaims", &
                               test_cfc_exclusion_radius), &
                  new_unittest("cfc_exclusion_radius_gated_on_exponents", &
                               test_cfc_exclusion_radius_gate), &
                  !> Empty active lists
                  new_unittest("svdw_empty_active", test_svdw_empty_active), &
                  new_unittest("cfc_empty_active", test_cfc_empty_active), &
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
   !> SvdW-only blend overrides; CFC uses fixed shape parameters
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
                       blend_3b, screening_threshold)
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
      !> Optional screening threshold (default 0, no screening)
      real(wp), intent(in), optional :: screening_threshold

      type(moist_cavity_drop_lsf_svdw_type) :: tmp_svdw
      type(moist_cavity_drop_lsf_cfc_type)  :: tmp_cfc
      type(moist_cavity_drop_lsf_sses_type) :: tmp_sses
      type(moist_cavity_drop_lsf_svdw_param_type) :: param
      real(wp) :: thr

      thr = 0.0_wp
      if (present(screening_threshold)) thr = screening_threshold

      select case (kind)
      case ("svdw")
         tmp_svdw%screening_threshold = thr
         call param%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
         call tmp_svdw%new(param)
         call tmp_svdw%update(mol, radii)
         call tmp_svdw%set_max_deriv(max_deriv)
         allocate (lsf, source=tmp_svdw)
      case ("cfc")
         tmp_cfc%screening_threshold = thr
         call tmp_cfc%new()
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

      ! Beyond molecular bounds along +x
      outside_pt = [maxval(mol%xyz(1, :)) + 20.0_wp, 0.0_wp, 0.0_wp]
      call lsf%prepare(outside_pt, lsf_err)
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
                             blend_3b=svdw_legacy_blend_3b, screening_threshold=0.0_wp)
               call init_lsf(lsf_scr, super_mol, radii, 0, kind, &
                             blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                             blend_3b=svdw_legacy_blend_3b, screening_threshold=thr)

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
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
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
                  call lsf%f012_r(lsf2_rr=analytic)
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
      real(wp) :: analytic(ndim, ndim, ndim), dummy_third(ndim, ndim, ndim)
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
                  call lsf%f3_rrr(lsf3_rrr=analytic)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=hess_pp, lsf3_rrr=dummy_third)
                     shifted = point; shifted(axis) = point(axis) + eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=hess_p, lsf3_rrr=dummy_third)
                     shifted = point; shifted(axis) = point(axis) - eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=hess_m, lsf3_rrr=dummy_third)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*eps
                     call lsf%prepare(shifted, lsf_err)
                     call lsf%f3_rrr(lsf2_rr=hess_mm, lsf3_rrr=dummy_third)
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
      real(wp), allocatable :: dummy_3rd(:, :, :, :)
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
         allocate (dummy_3rd(ndim, ndim, ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
      real(wp), allocatable :: dummy_3rd(:, :, :, :)
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
         allocate (dummy_3rd(ndim, ndim, ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
      real(wp), allocatable :: analytic(:, :, :, :)
      real(wp) :: numeric, hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
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
         allocate (analytic(ndim, ndim, ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
                  call lsf%f3_rr_rA(lsf3_rr_rA=analytic)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=hess_pp)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=hess_p)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=hess_m)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point, lsf_err); call lsf%f012_r(lsf2_rr=hess_mm)
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
         deallocate (centers_base, centers_local, analytic)
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
      real(wp), allocatable :: dummy_2nd(:, :, :)
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
         allocate (dummy_2nd(ndim, ndim, mol%nat))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
      real(wp), allocatable :: dummy_2nd(:, :, :)
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
         allocate (dummy_2nd(ndim, ndim, mol%nat))
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
   subroutine init_lsf_sweep(lsf, mol, radii, max_deriv, kind, isweep)
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
      integer :: nblend, ngamma, iblend, igamma

      call svdw_sweep_sizes(kind, nblend, ngamma)
      iblend = (isweep - 1)/ngamma + 1
      igamma = mod(isweep - 1, ngamma) + 1
      call init_lsf(lsf, mol, radii, max_deriv, kind, &
                    blend_k=svdw_sweep_blend(kind, iblend), &
                    blend_2b=svdw_sweep_2b(kind, igamma), &
                    blend_3b=svdw_sweep_gamma(kind, igamma))
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
      real(wp), allocatable :: analytic(:, :, :)
      real(wp) :: numeric
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
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
                             blend_k=svdw_sweep_blend(kind, iblend), &
                             blend_2b=svdw_sweep_2b(kind, igamma), &
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate (analytic(ndim, ndim, mol%nat))
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
                  call lsf%f3_rr_rad(lsf3_rr_rad=analytic)
                  do atom = 1, mol%nat
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=h_pp)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) + eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=h_p)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=h_m)
                     radii_local = radii
                     radii_local(atom) = radii_local(atom) - 2.0_wp*eps
                     call refresh_radii(lsf, mol, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f012_r(lsf2_rr=h_mm)
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
               deallocate (analytic, lsf)
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
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
                  call lsf%f4_rr_radrad(rr4)
                  do atom = 1, mol%nat
                     do istep = 1, 4
                        radii_local = radii
                        radii_local(atom) = radii_local(atom) + steps(istep)*eps
                        call refresh_radii(lsf, mol, radii_local)
                        if (prepare_failed(lsf, point, error)) return
                        call lsf%f3_rr_rad(lsf1_rad=l1(:, istep), &
                                           lsf2_r_rad=l2(:, :, istep), &
                                           lsf3_rr_rad=l3(:, :, :, istep))
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
         deallocate (radii_local, rr2, rr3, rr4, l1, l2, l3)
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
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
                  call lsf%f4_rr_rA_rad(mx4)
                  do atom = 1, mol%nat
                     do istep = 1, 4
                        radii_local = radii
                        radii_local(atom) = radii_local(atom) + steps(istep)*eps
                        call refresh_radii(lsf, mol, radii_local)
                        if (prepare_failed(lsf, point, error)) return
                        call lsf%f3_rr_rA(lsf1_rA=n1(:, :, istep), &
                                          lsf2_r_rA=n2(:, :, :, istep), &
                                          lsf3_rr_rA=n3(:, :, :, :, istep))
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
         deallocate (radii_local, mx2, mx3, mx4, n1, n2, n3)
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
            call lsf%f4_rr_radrad(rr4)
            call lsf%f2_rA_rad(mx2)
            call lsf%f3_r_rA_rad(mx3)
            call lsf%f4_rr_rA_rad(mx4)
            call lsf%f2_rArB(q2)
            call lsf%f3_r_rArB(q3)
            call lsf%f4_rr_rArB(q4)

            call lsf%hvp_f1_rad(v, vrad, r1)
            call lsf%hvp_f2_r_rad(v, vrad, r2)
            call lsf%hvp_f3_rr_rad(v, vrad, r3)
            call lsf%hvp_f1_rA(v, h1, vrad)
            call lsf%hvp_f2_r_rA(v, h2, vrad)
            call lsf%hvp_f3_rr_rA(v, h3, vrad)

            !* --------------------------------- the radius row -------------------------------- *!
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

            !* ----- the position row, whose radius half is the nuclear-radius block again ----- *!
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
                     r1, r2, r3, h1, h2, h3)
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
         centers_base = mol%xyz
         call joint_direction(mol%nat, v, vrad)
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
                  call lsf%hvp_f3_rr_rad(v, vrad, a3)
                  do istep = 1, 4
                     centers_local = centers_base + steps(istep)*eps*v
                     radii_local = radii + steps(istep)*eps*vrad
                     call refresh_joint(lsf, mol, centers_local, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f3_rr_rad(lsf1_rad=l1(:, istep), &
                                        lsf2_r_rad=l2(:, :, istep), &
                                        lsf3_rr_rad=l3(:, :, :, istep))
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
                     a1, a2, a3, l1, l2, l3)
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
         centers_base = mol%xyz
         call joint_direction(mol%nat, v, vrad)
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
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
                  call lsf%hvp_f3_rr_rA(v, h3, vrad)
                  do istep = 1, 4
                     centers_local = centers_base + steps(istep)*eps*v
                     radii_local = radii + steps(istep)*eps*vrad
                     call refresh_joint(lsf, mol, centers_local, radii_local)
                     if (prepare_failed(lsf, point, error)) return
                     call lsf%f3_rr_rA(lsf1_rA=n1(:, :, istep), &
                                       lsf2_r_rA=n2(:, :, :, istep), &
                                       lsf3_rr_rA=n3(:, :, :, :, istep))
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
                     h1, h2, h3, n1, n2, n3)
      end do
   end subroutine run_hvp_rA_joint_fd

   !* ================================================================================= *!
   !*                        Accessors at an empty active list                          *!
   !* ================================================================================= *!
   !
   ! Far-field screening may reject every atom
   ! `intent(out)` results still require defined buffers at active_count() = 0
   ! Inspect full buffers; other checks inspect only active prefixes
   !
   ! Surviving buffer poison exposes undefined accessor results

   !> SvdW dispatch for the empty-active-list check
   subroutine test_svdw_empty_active(error)
      type(error_type), allocatable, intent(out) :: error
      call run_empty_active(error, "svdw")
   end subroutine test_svdw_empty_active

   !> CFC dispatch for the empty-active-list check
   subroutine test_cfc_empty_active(error)
      type(error_type), allocatable, intent(out) :: error
      call run_empty_active(error, "cfc")
   end subroutine test_cfc_empty_active

   !> Check radius-accessor results with no active atoms
   subroutine run_empty_active(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), v(:, :), vrad(:)
      real(wp), allocatable :: g1(:), g2(:, :), g3(:, :, :)
      real(wp), allocatable :: rr2(:, :), rr3(:, :, :), rr4(:, :, :, :)
      real(wp), allocatable :: mx2(:, :, :), mx3(:, :, :, :), mx4(:, :, :, :, :)
      real(wp), allocatable :: h1(:, :), h2(:, :, :), h3(:, :, :, :)
      real(wp) :: point(ndim)
      integer  :: nat

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)
      nat = mol%nat

      ! Nonzero threshold permits an empty active list
      call init_lsf(lsf, mol, radii, 3, kind, screening_threshold=1.0e-10_wp)
      point = [5.0e2_wp, 5.0e2_wp, 5.0e2_wp]
      if (prepare_failed(lsf, point, error)) return

      call check(error, lsf%active_count(), 0, &
                 message="fixture must leave no active atom; the rest of this "// &
                 "test says nothing otherwise")
      if (allocated(error)) return

      allocate (v(ndim, nat), vrad(nat))
      call joint_direction(nat, v, vrad)

      allocate (g1(nat), g2(ndim, nat), g3(ndim, ndim, nat))
      allocate (rr2(nat, nat), rr3(ndim, nat, nat), rr4(ndim, ndim, nat, nat))
      allocate (mx2(ndim, nat, nat), mx3(ndim, ndim, nat, nat))
      allocate (mx4(ndim, ndim, ndim, nat, nat))
      allocate (h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat))

      ! Uncontracted radius ladder
      g1 = EMPTY_POISON; g2 = EMPTY_POISON; g3 = EMPTY_POISON
      call lsf%f3_rr_rad(g1, g2, g3)
      call check_all_zero(error, reshape(g1, [size(g1)]), "f1_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(g2, [size(g2)]), "f2_r_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(g3, [size(g3)]), "f3_rr_rad")
      if (allocated(error)) return

      ! Uncontracted two-radius and nuclear-radius blocks
      rr2 = EMPTY_POISON; rr3 = EMPTY_POISON; rr4 = EMPTY_POISON
      mx2 = EMPTY_POISON; mx3 = EMPTY_POISON; mx4 = EMPTY_POISON
      call lsf%f2_radrad(rr2)
      call lsf%f3_r_radrad(rr3)
      call lsf%f4_rr_radrad(rr4)
      call lsf%f2_rA_rad(mx2)
      call lsf%f3_r_rA_rad(mx3)
      call lsf%f4_rr_rA_rad(mx4)
      call check_all_zero(error, reshape(rr2, [size(rr2)]), "f2_radrad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(rr3, [size(rr3)]), "f3_r_radrad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(rr4, [size(rr4)]), "f4_rr_radrad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(mx2, [size(mx2)]), "f2_rA_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(mx3, [size(mx3)]), "f3_r_rA_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(mx4, [size(mx4)]), "f4_rr_rA_rad")
      if (allocated(error)) return

      ! Radius row of the joint HVP
      g1 = EMPTY_POISON; g2 = EMPTY_POISON; g3 = EMPTY_POISON
      call lsf%hvp_f1_rad(v, vrad, g1)
      call lsf%hvp_f2_r_rad(v, vrad, g2)
      call lsf%hvp_f3_rr_rad(v, vrad, g3)
      call check_all_zero(error, reshape(g1, [size(g1)]), "hvp_f1_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(g2, [size(g2)]), "hvp_f2_r_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(g3, [size(g3)]), "hvp_f3_rr_rad")
      if (allocated(error)) return

      ! CFC nuclear-row branches with and without radius direction
      h1 = EMPTY_POISON; h2 = EMPTY_POISON; h3 = EMPTY_POISON
      call lsf%hvp_f1_rA(v, h1)
      call lsf%hvp_f2_r_rA(v, h2)
      call lsf%hvp_f3_rr_rA(v, h3)
      call check_all_zero(error, reshape(h1, [size(h1)]), "hvp_f1_rA")
      if (allocated(error)) return
      call check_all_zero(error, reshape(h2, [size(h2)]), "hvp_f2_r_rA")
      if (allocated(error)) return
      call check_all_zero(error, reshape(h3, [size(h3)]), "hvp_f3_rr_rA")
      if (allocated(error)) return

      h1 = EMPTY_POISON; h2 = EMPTY_POISON; h3 = EMPTY_POISON
      call lsf%hvp_f1_rA(v, h1, vrad)
      call lsf%hvp_f2_r_rA(v, h2, vrad)
      call lsf%hvp_f3_rr_rA(v, h3, vrad)
      call check_all_zero(error, reshape(h1, [size(h1)]), "hvp_f1_rA (joint)")
      if (allocated(error)) return
      call check_all_zero(error, reshape(h2, [size(h2)]), "hvp_f2_r_rA (joint)")
      if (allocated(error)) return
      call check_all_zero(error, reshape(h3, [size(h3)]), "hvp_f3_rr_rA (joint)")
   end subroutine run_empty_active

   !> Fail unless every element of `a` is exactly zero
   !>
   !> @param[out] error Set when a poisoned (or any nonzero) element survives
   !> @param[in]  a     Flattened result buffer
   !> @param[in]  name  Accessor name, for the failure message
   subroutine check_all_zero(error, a, name)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Flattened result buffer
      real(wp), intent(in) :: a(:)
      !> Accessor name
      character(len=*), intent(in) :: name

      call check(error, maxval(abs(a)), 0.0_wp, thr=0.0_wp, &
                 message=name//" left its result undefined at an empty active list")
   end subroutine check_all_zero

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
      real(wp) :: v0, v1, g0(ndim), g1(ndim), h(ndim, ndim)
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
      real(wp), allocatable :: dummy_rr_rA(:, :, :, :)
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
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol%nat))
         allocate (analytic(ndim, mol%nat, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 3, kind, isweep)
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
      real(wp), allocatable :: dummy_rr_rA(:, :, :, :)
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
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol%nat))
         allocate (analytic(ndim, ndim, mol%nat, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 3, kind, isweep)
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
      real(wp) :: analytic(ndim, ndim, ndim, ndim)
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
               call prim%f4_rrrr(analytic)
               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=t3_fwd)
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=t3_fwd2)
                  work_point = point
                  work_point(axis) = point(axis) - STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=t3_bwd)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*STEP_SIZE
                  call prim%prepare(work_point, lsf_err)
                  call prim%f3_rrr(lsf3_rrr=t3_bwd2)
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
      real(wp), allocatable :: analytic(:, :, :, :, :)
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
         if (allocated(analytic)) deallocate (analytic)
         allocate (analytic(ndim, ndim, ndim, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 4, kind, isweep)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f4_rrr_rA(analytic)
               do atom = 1, mol%nat
                  do axis = 1, ndim
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=t3_fwd)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=t3_fwd2)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=t3_bwd)
                     centers_local = centers_base
                     centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rrr(lsf3_rrr=t3_bwd2)
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
         if (allocated(analytic)) deallocate (analytic)
         if (allocated(rr_rA_fwd)) deallocate (rr_rA_fwd, rr_rA_fwd2, rr_rA_bwd, rr_rA_bwd2)
         allocate (analytic(ndim, ndim, ndim, mol%nat, ndim, mol%nat))
         allocate (rr_rA_fwd(ndim, ndim, ndim, mol%nat), rr_rA_fwd2(ndim, ndim, ndim, mol%nat))
         allocate (rr_rA_bwd(ndim, ndim, ndim, mol%nat), rr_rA_bwd2(ndim, ndim, ndim, mol%nat))

         do isweep = 1, sweep_count(kind)
            call init_lsf_sweep(prim, mol, radii, 4, kind, isweep)
            do ipt = 1, size(points, 2)
               point = points(:, ipt)
               call prim%set_centers(centers_base)
               call prim%prepare(point, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call prim%f4_rr_rArB(analytic)
               do iB = 1, mol%nat
                  do axisB = 1, ndim
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) + STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_fwd)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) + 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_fwd2)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) - STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_bwd)
                     centers_local = centers_base
                     centers_local(axisB, iB) = centers_local(axisB, iB) - 2.0_wp*STEP_SIZE
                     call prim%set_centers(centers_local)
                     call prim%prepare(point, lsf_err)
                     call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_bwd2)
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
            call init_lsf_sweep(prim, mol, radii, 2, kind, isweep)
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
   subroutine run_contracted(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
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
      real(wp), allocatable :: rad1(:), rad2(:, :), rad3(:, :, :)
      real(wp) :: w0_adj, w1_adj(ndim), w2_adj(ndim, ndim)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: ref0, ref1(ndim), ref2(ndim, ndim), ref3(ndim, ndim, ndim)
      real(wp) :: acc, point(ndim)
      integer  :: icase, ipt, iA, iB, s_ax, t_ax, j, k, l, nat, atomA
      type(mctc_error), allocatable :: lsf_err

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
                                           vjp, vjp_rad, rad1, rad2, rad3)
         allocate (f1_rA(ndim, nat), f2_r_rA(ndim, ndim, nat))
         allocate (f3_rr_rA(ndim, ndim, ndim, nat))
         allocate (f4_rrr_rA(ndim, ndim, ndim, ndim, nat))
         allocate (f2_rArB(ndim, nat, ndim, nat))
         allocate (f3_r_rArB(ndim, ndim, nat, ndim, nat))
         allocate (f4_rr_rArB(ndim, ndim, ndim, nat, ndim, nat))
         allocate (h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat))
         allocate (vjp(ndim, nat), vjp_rad(nat))
         allocate (rad1(nat), rad2(ndim, nat), rad3(ndim, ndim, nat))

         call init_lsf(prim, mol, radii, 4, kind, blend_k=svdw_legacy_blend_k, &
                       blend_2b=svdw_legacy_blend_2b, blend_3b=svdw_legacy_blend_3b)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            call prim%prepare(point, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if

            call prim%f3_rr_rA(f1_rA, f2_r_rA, f3_rr_rA)
            call prim%f4_rrr_rA(f4_rrr_rA)
            call prim%f2_rArB(f2_rArB)
            call prim%f3_r_rArB(f3_r_rArB)
            call prim%f4_rr_rArB(f4_rr_rArB)

            !* ----------------------------------- tangent_* ----------------------------------- *!
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
            call prim%tangent_f2_rr(v, t2)
            do k = 1, ndim
               do j = 1, ndim
                  call check(error, t2(j, k), ref2(j, k), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
            call prim%tangent_f3_rrr(v, t3)
            do l = 1, ndim
               do k = 1, ndim
                  do j = 1, ndim
                     call check(error, t3(j, k, l), ref3(j, k, l), &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do

            !* ------------------------------------- hvp_* ------------------------------------- *!
            call prim%hvp_f1_rA(v, h1)
            call prim%hvp_f2_r_rA(v, h2)
            call prim%hvp_f3_rr_rA(v, h3)
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

            !* ----------------------------------- vjp_f1_rA ----------------------------------- *!
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

            call prim%vjp_f1_rA(w0_adj, w1_adj, w2_adj, vjp)
            do iA = 1, prim%active_count()
               do s_ax = 1, ndim
                  acc = w0_adj*f1_rA(s_ax, iA) &
                        + dot_product(w1_adj, f2_r_rA(:, s_ax, iA)) &
                        + sum(w2_adj*f3_rr_rA(:, :, s_ax, iA))
                  call check(error, vjp(s_ax, iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do

            !* ----------------------------------- vjp_f1_rad ---------------------------------- *!
            ! Radius-ladder contraction; one scalar per active atom
            call prim%f3_rr_rad(rad1, rad2, rad3)
            call prim%vjp_f1_rad(w0_adj, w1_adj, w2_adj, vjp_rad)
            do iA = 1, prim%active_count()
               acc = w0_adj*rad1(iA) &
                     + dot_product(w1_adj, rad2(:, iA)) &
                     + sum(w2_adj*rad3(:, :, iA))
               call check(error, vjp_rad(iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do

            !* ----------------------- nuclear Hessian exchange symmetry ----------------------- *!
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
      real(wp) :: f0_on, f0_near, lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim)
      real(wp) :: f3_rrr(ndim, ndim, ndim), f4_rrrr(ndim, ndim, ndim, ndim)
      real(wp) :: f1_rA(ndim, nat), f2_r_rA(ndim, ndim, nat)
      real(wp) :: f3_rr_rA(ndim, ndim, ndim, nat)
      real(wp) :: f4_rrr_rA(ndim, ndim, ndim, ndim, nat)
      real(wp) :: f2_rArB(ndim, nat, ndim, nat)
      real(wp) :: f3_r_rArB(ndim, ndim, nat, ndim, nat)
      real(wp) :: f4_rr_rArB(ndim, ndim, ndim, nat, ndim, nat)
      real(wp) :: norm0, norm1_rA(ndim, nat)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat)
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
      call check(error, finite0(lsf0) .and. finite1(lsf1_r) .and. finite2(lsf2_rr), &
                 "f012_r is not finite on a nucleus")
      if (allocated(error)) return

      call prim%f3_rrr(lsf3_rrr=f3_rrr)
      call check(error, finite3(f3_rrr), "f3_rrr is not finite on a nucleus")
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
      call check(error, finite0(t0) .and. finite1(t1) .and. finite2(t2) .and. finite3(t3), &
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

   !> Finiteness of a matrix
   !>
   !> @param[in] x Matrix to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finite2(x) result(ok)
      !> Matrix to test
      real(wp), intent(in) :: x(:, :)
      ok = finiten(reshape(x, [size(x)]))
   end function finite2

   !> Finiteness of a rank-3 tensor
   !>
   !> @param[in] x Tensor to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finite3(x) result(ok)
      !> Tensor to test
      real(wp), intent(in) :: x(:, :, :)
      ok = finiten(reshape(x, [size(x)]))
   end function finite3

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
      real(wp), allocatable :: radii(:), points(:, :), shares(:), f1_rad(:), f3_rr_rad(:, :, :)
      integer :: icase, ipt, iblend, igamma, nat

      call get_test_structures(mols)
      do icase = 1, size(mols)
         nat = mols(icase)%nat
         call get_test_radii(mols(icase), radii)
         call get_test_points(mols(icase), points, n_points_svdw)
         allocate (shares(nat), f1_rad(nat), f3_rr_rad(ndim, ndim, nat))
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               call init_lsf(lsf, mols(icase), radii, 2, "svdw", blend_k=svdw_blend_k_values(iblend), &
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
      real(wp), allocatable :: f1_rad(:), f3_rr_rad(:, :, :)
      real(wp) :: p, eps
      integer :: icase, ipt, nat

      call get_test_structures(mols)
      do icase = 1, size(mols)
         nat = mols(icase)%nat
         if (nat < 2) cycle
         call get_test_radii(mols(icase), radii)
         call get_test_points(mols(icase), points, n_points_sses)
         allocate (shares(nat), slope(nat), f1_rad(nat), f3_rr_rad(ndim, ndim, nat))
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
   !*                        Surface-free exclusion certificate                          *!
   !* ================================================================================= *!

   !> Probe certified surface-free ball boundaries
   !>
   !> `exclusion_radius(S(x))` promises no zero inside `B(x, r)`
   !> Require unchanged sign across each claimed boundary
   !>
   !> Probe `0.999 r` to allow tight SvdW `r = |S|` claims
   !>
   !> Caller checks `n_certified > 0` to avoid vacuous success
   !>
   !> @param[inout] lsf          LSF bound to the structure the points came from
   !> @param[in]    points       Sampling points (3, npoint)
   !> @param[out]   n_certified  Points that claimed a positive radius
   !> @param[out]   error        Test failure
   subroutine probe_exclusion_balls(lsf, points, n_certified, error)
      !> LSF bound to the structure the points came from
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Sampling points
      real(wp), intent(in) :: points(:, :)
      !> Points that claimed a positive radius
      integer, intent(out) :: n_certified
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Fraction of the claimed radius actually probed
      real(wp), parameter :: probe_frac = 0.999_wp

      type(mctc_error), allocatable :: lsf_err
      real(wp) :: center(ndim), probe(ndim), dir(ndim)
      real(wp) :: lsf_center, lsf_probe, r
      integer :: ipoint, idir

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

      n_certified = 0

      do ipoint = 1, size(points, 2)
         center = points(:, ipoint)

         call lsf%prepare(center, lsf_err)
         if (allocated(lsf_err)) then
            call test_failed(error, "LSF prepare failed: "//lsf_err%message)
            return
         end if
         call lsf%f0(lsf_center)

         r = lsf%exclusion_radius(lsf_center)
         if (r <= 0.0_wp) cycle
         n_certified = n_certified + 1

         do idir = 1, ndir
            dir = dirs(:, idir)
            probe = center + probe_frac*r*dir/norm2(dir)

            call lsf%prepare(probe, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if
            call lsf%f0(lsf_probe)

            if (lsf_probe*lsf_center <= 0.0_wp) then
               call test_failed(error, &
                                "exclusion radius over-claims: S changes sign inside "// &
                                "the ball it certifies as surface-free")
               return
            end if
         end do
      end do
   end subroutine probe_exclusion_balls

   !> Check SvdW exclusion radius across blend sharpness
   !>
   !> Tight `r = |S(x)|` claim from 1-Lipschitz property
   subroutine test_svdw_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: iblend, n_certified

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call get_test_points(mol, points)

      do iblend = 1, n_svdw_blends
         call init_lsf(lsf, mol, radii, 0, "svdw", &
                       blend_k=svdw_blend_k_values(iblend))

         call probe_exclusion_balls(lsf, points, n_certified, error)
         if (allocated(error)) return

         if (n_certified == 0) then
            call test_failed(error, "SvdW certified no ball at all: the probe "// &
                             "loop never ran")
            return
         end if
      end do
   end subroutine test_svdw_exclusion_radius

   !> Negative blend coefficients invalidate 1-Lipschitz bound
   !> Require no certificate under extrapolated weights
   subroutine test_svdw_exclusion_radius_gate(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)

      call init_lsf(lsf, mol, radii, 0, "svdw", blend_2b=1.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 1.0_wp, thr=0.0_wp, &
                 message="non-negative blends should certify r = |S|")
      if (allocated(error)) return

      call init_lsf(lsf, mol, radii, 0, "svdw", blend_2b=-1.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="a negative blend coefficient must certify nothing")
   end subroutine test_svdw_exclusion_radius_gate

   !> Check CFC exclusion radius
   !>
   !> `S = -log PD`; claim `r = |S| / L` with cached Lipschitz bound
   !> Conservative `L >= sup ||grad S||` favors smaller, safe balls
   subroutine test_cfc_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, n_certified, n_total

      ! Sweep fixture sizes and radii for `(n - 1)**(1/m)` bound
      ! Length scale uses two smallest radii
      call get_test_structures(mols, 10)

      n_total = 0
      do icase = 1, size(mols)
         call get_test_radii(mols(icase), radii)
         call get_test_points(mols(icase), points)

         call init_lsf(lsf, mols(icase), radii, 0, "cfc")

         call probe_exclusion_balls(lsf, points, n_certified, error)
         if (allocated(error)) return
         n_total = n_total + n_certified
      end do

      if (n_total == 0) then
         call test_failed(error, "CFC certified no ball at all: the probe loop "// &
                          "never ran")
      end if
   end subroutine test_cfc_exclusion_radius

   !> CFC certificate parameter limits
   !>
   !> - Decaying exponents: `a1, a2 < 0`
   !> - Nonnegative pair term: `c >= 0`, even `m`; positive `PD`
   !> - Faster pair decay: `2|a2| > |a1|`
   !> Any violation disables certification
   subroutine test_cfc_exclusion_radius_gate(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      real(wp) :: r

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)

      ! LiH CPCM radii 2.46 and 3.87 Bohr; L = 8.4 / Bohr
      ! Certificate 0.119 Bohr per |S|; catches radius, exponent, or unit errors
      call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, m=4))
      call cfc%update(mol, radii)
      r = cfc%exclusion_radius(-1.0_wp)
      if (r < 0.05_wp .or. r > 0.2_wp) then
         call test_failed(error, "default CFC parameters must certify a radius "// &
                          "near 0.1 Bohr per unit of |S|")
         return
      end if

      ! Certificate linear in supplied value
      call check(error, cfc%exclusion_radius(-3.0_wp), 3.0_wp*r, thr=0.0_wp, &
                 message="the certified radius must scale with |S|")
      if (allocated(error)) return

      ! `2|a2| < |a1|`: pair term outlives atomic terms
      call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp, a2=-7.0_wp, c=5.0_wp, m=4))
      call cfc%update(mol, radii)
      call check(error, cfc%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="a pair term decaying slower than the atomic one "// &
                 "must certify nothing")
      if (allocated(error)) return

      ! Odd `m`: `u**m` may make `PD` negative
      call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, m=3))
      call cfc%update(mol, radii)
      call check(error, cfc%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="an odd pair power must certify nothing")
      if (allocated(error)) return

      ! Negative coupling can make `PD` negative
      call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp, a2=-9.0_wp, c=-5.0_wp, m=4))
      call cfc%update(mol, radii)
      call check(error, cfc%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="a negative pair coupling must certify nothing")
      if (allocated(error)) return

      ! Growing atomic term lies outside the bound
      call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=1.0_wp, a2=-9.0_wp, c=5.0_wp, m=4))
      call cfc%update(mol, radii)
      call check(error, cfc%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="a non-decaying atomic term must certify nothing")
   end subroutine test_cfc_exclusion_radius_gate

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
      real(wp) :: dr(1), gr(3, 1), hr(3, 3, 1), rr(1, 1)
      real(wp) :: point(3), expected(3, 3)
      integer :: i

      call sses_setup(lsf, mol, radii, 1, 4)
      point = mol%xyz(:, 1) + [0.2_wp, 0.3_wp, 0.4_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f3_rrr(f, g, h, t)
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
      call lsf%f4_rrrr(fourth)
      call check(error, maxval(abs(t)) + maxval(abs(fourth)) < 1.0e-14_wp)
      if (allocated(error)) return
      call lsf%f3_rr_rad(dr, gr, hr)
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
      real(wp) :: f, g(3), h(3, 3), fj, gj(3), hj(3, 3)
      real(wp) :: n1(3, 3), n2(3, 3, 3), n3(3, 3, 3, 3), j1(3, 4), j2(3, 3, 4), j3(3, 3, 3, 4)
      real(wp) :: r1(3), r2(3, 3), r3(3, 3, 3), k1(4), k2(3, 4), k3(3, 3, 4)
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
                    maxval(abs(n3 - j3(:, :, :, 1:3))))
         diff = max(diff, maxval(abs(r1 - k1(1:3))), maxval(abs(r2 - k2(:, 1:3))), &
                    maxval(abs(r3 - k3(:, :, 1:3))))
         call check(error, diff < thr*(1.0_wp + abs(f)), "distant atom changes the water field")
         if (allocated(error)) return
         diff = max(maxval(abs(j1(:, 4))), maxval(abs(j2(:, :, 4))), maxval(abs(j3(:, :, :, 4))), &
                    abs(k1(4)), maxval(abs(k2(:, 4))), maxval(abs(k3(:, :, 4))))
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
      type(moist_cavity_drop_swif_sigmoid_bump_type) :: critical
      type(structure_type) :: mol
      real(wp) :: centers(3, 2), point(3), a, y, f, g(3), h(3, 3), norm, tau
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
      call new_swif_sigmoid_bump(critical, param%w_0ls_from, param%w_0ls_to, &
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
      real(wp) :: point(3), shift(3), rotation(3, 3), f, g(3), h(3, 3), f2, g2(3), h2(3, 3)

      call sses_setup(lsf, mol, radii, 3, 2)
      point = [1.1_wp, 0.7_wp, 0.2_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f012_r(f, g, h)
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
      call lsf%f012_r(f2, g2, h2)
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
      real(wp) :: f, g(3), h(3, 3), third(3, 3, 3), fourth(3, 3, 3, 3)
      real(wp) :: f2, g2(3), h2(3, 3), t2(3, 3, 3), q2(3, 3, 3, 3)
      real(wp) :: nuclear(3, 3), mixed(3, 3, 3), high(3, 3, 3, 3), nuclear2(3, 3)
      real(wp) :: mixed2(3, 3, 3), high2(3, 3, 3, 3), radial(3), radial2(3)
      real(wp) :: gr(3, 3), gr2(3, 3), hr(3, 3, 3), hr2(3, 3, 3)
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
      integer :: i, j, k, a, b, s, t

      call sses_setup(lsf, mol, radii, 3, 3)
      point = [1.1_wp, 0.7_wp, 0.2_wp]
      v = reshape([0.2_wp, -0.1_wp, 0.3_wp, -0.2_wp, 0.4_wp, 0.1_wp, 0.3_wp, 0.2_wp, -0.4_wp], [3, 3])
      vr = [0.1_wp, -0.2_wp, 0.3_wp]
      call sses_prepare(lsf, point, error)
      if (allocated(error)) return
      call lsf%f3_rr_rA(n1, n2, n3)
      call lsf%f3_rr_rad(r1, r2, r3)
      call lsf%f4_rrr_rA(n4)
      call lsf%f2_rArB(nn)
      call lsf%f3_r_rArB(nn1)
      call lsf%f4_rr_rArB(nn2)
      call lsf%f2_rA_rad(rn)
      call lsf%f3_r_rA_rad(rn1)
      call lsf%f4_rr_rA_rad(rn2)
      call lsf%f2_radrad(rr)
      call lsf%f3_r_radrad(rr1)
      call lsf%f4_rr_radrad(rr2)
      call lsf%hvp_f1_rA(v, h1, vrad=vr)
      call lsf%hvp_f2_r_rA(v, h2, vrad=vr)
      call lsf%hvp_f3_rr_rA(v, h3, vrad=vr)
      call lsf%hvp_f1_rad(v, vr, hr1)
      call lsf%hvp_f2_r_rad(v, vr, hr2)
      call lsf%hvp_f3_rr_rad(v, vr, hr3)
      call lsf%tangent_f0(v, t0)
      call lsf%tangent_f1_r(v, t1)
      call lsf%tangent_f2_rr(v, t2)
      call lsf%tangent_f3_rrr(v, t3)
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
      param%proj_level = 2
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

end module test_cavity_drop_lsf
