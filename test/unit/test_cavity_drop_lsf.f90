!> Orchestrator-level unit tests for the DROP level set functions (LSF)
!>
!> Each common FD check is written once against the abstract
!> [[moist_cavity_drop_lsf_type]] and dispatched twice in the collector:
!> once for the SvdW concrete ([[moist_cavity_drop_lsf_svdw_type]]) and
!> once for the CFC concrete ([[moist_cavity_drop_lsf_cfc_type]]). The
!> SvdW dispatch additionally sweeps over the blend_k x gamma parameter
!> grid that exercises SvdW's body-order weights; CFC has no analogous
!> knob and runs with its compiled defaults only.
!>
!> Common checks (run for both concretes):
!>   * `sign_convention`     interior < 0, exterior > 0
!>   * `neighbor_cutoff`     screening-threshold sweep vs threshold=0 ref
!>   * `f1_r_fd`             grad vs 4-point central FD of f0
!>   * `f2_rr_fd`            Hessian vs FD of grad
!>   * `f3_rrr_fd`           third spatial deriv vs FD of Hessian
!>   * `f1_rA_fd`            nuclear grad vs FD of f0
!>   * `f2_r_rA_fd`          mixed deriv vs FD of grad
!>   * `f3_rr_rA_fd`         mixed third deriv vs FD of Hessian
!>
!> Radius-derivative checks, both concretes. Their radius channels rest on
!> opposite structures -- SvdW's `d/dR_a` is a diagonal rescaling of an atom's
!> own tensors, CFC's a genuine extra derivative of a two-center kernel -- so
!> the same FD says two different things:
!>   * `f1_rad_fd`           radius grad vs FD of f0
!>   * `f2_r_rad_fd`         mixed spatial-radius deriv vs FD of grad
!>   * `f3_rr_rad_fd`        mixed third radius deriv vs FD of Hessian
!>   * `radrad_fd`           `f{2,3,4}_*_radrad` vs a one-radius FD of the
!>                           `f*_rad` ladder, element by element
!>   * `ra_rad_fd`           `f{2,3,4}_*_rA_rad` vs a one-radius FD of the
!>                           `f*_rA` ladder. Also pins the *order* of the two
!>                           atom slots, which that block is not symmetric in
!>   * `hvp_rad_fd`          radius row of the joint position/radius HVP vs
!>                           FD of the `f*_rad` ladder along `(v, vr)`
!>   * `hvp_ra_joint_fd`     nuclear row of the same HVP, `vrad` supplied, vs
!>                           FD of the `f*_rA` ladder along `(v, vr)`
!>   * `radius_pairwise`     both HVP rows against the *uncontracted* blocks
!>                           `f*_radrad` and `f*_rA_rad`. No finite differences:
!>                           an exact contraction identity, so it pins the two
!>                           uncontracted families at roundoff where every other
!>                           radius test only reaches FD accuracy.
!>   * `empty_active`        every accessor above at a point where screening
!>                           rejects all atoms. The one check that inspects a
!>                           whole result buffer instead of its active prefix.
!> All the FD ones above take CFC's looser thresholds; see [[radius_fd_thr]].
!>
!> High-order extensions, registered once per concrete because the SvdW
!> dispatch sweeps a blend_k x gamma parameter grid that CFC has no analogue of
!> (its four shape parameters are compiled-in constants):
!>   * `f2_rArB`, `f3_r_rArB`                  pure/mixed nuclear seconds
!>   * `f4_rrrr`, `f4_rrr_rA`, `f4_rr_rArB`    fourth-order derivatives
!>   * `normalized_f1_rA`                      normalized LSF nuclear grad
!>   * `contracted_vs_pairwise`                every `tangent_*` / `hvp_*`
!>                                             against the uncontracted tensors
!>
!> SvdW-only:
!>   * `body_order_scaling`                    1b/2b/3b weight reduction
!>   * `point_on_nucleus_is_finite`            CFC has no on-nucleus singularity
!>   * `exclusion_radius_*`                    CFC certifies nothing
!>
!> Test fixtures come from the shared MB16-43/Heavy28/Amino20x4/But14diol/UPU23
!> palette in `test_helpers::get_test_structures`. Per-atom radii come from
!> the project's standard CPCM table via `get_test_radii`, and sampling
!> points come from `get_test_points`.
!>
!> Nuclear-derivative FDs bypass `lsf%update(mol, radii)` (which would
!> re-init the concrete's per-atom caches and reset `max_deriv`) by calling
!> the base's `lsf%set_centers(...)`, which only moves the atoms and rebuilds
!> the spatial sort and screening bounds.
module test_cavity_drop_lsf
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, &
                           center_at_origin, fd4_scalar
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_model_gems_utils, only: BuildSuperStructure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_cavity_drop_lsf

   integer, parameter :: ndim = 3

   !> Concrete LSF selector strings; used by [[init_lsf]]
   character(len=*), parameter :: kind_svdw = "svdw"
   character(len=*), parameter :: kind_cfc = "cfc"

   !> Tolerance of the Hessian-free consistency check. The two dispatch
   !> branches of the generated SvdW kernel share a value and a gradient but
   !> not a common-subexpression schedule, so they agree to roundoff, not bits.
   real(wp), parameter :: HESSFREE_ABS = 1.0e-14_wp
   real(wp), parameter :: HESSFREE_REL = 1.0e-12_wp

   !> Poison written into every result buffer by the empty-active-list check.
   real(wp), parameter :: EMPTY_POISON = -1.0e30_wp

   real(wp), parameter :: STEP_SIZE = 1.0e-3_wp
   real(wp), parameter :: ABS_THR = 2.0e-10_wp
   real(wp), parameter :: REL_THR = 1.0e-9_wp

   !> Finite-difference thresholds of the CFC-only high-order block.
   !>
   !> CFC's exponents (a1 = -15, a2 = -9) make the level set vary on a length
   !> scale of R/15 rather than SvdW's 3/k, so at the shared `STEP_SIZE` a
   !> fourth-order central difference carries correspondingly more truncation
   !> error. Measured, not guessed: the order-4 checks (`f4_rrrr`, `f4_rrr_rA`,
   !> `f4_rr_rArB`) fail at SvdW's `2e-10 / 1e-9` and pass at `1e-9 / 1e-8`, so
   !> the true stencil error sits in that decade; these values leave one decade
   !> of headroom above it for other platforms. The orders below four pass at
   !> the SvdW thresholds and are not relaxed by this being loose -- they are
   !> additionally pinned at roundoff by `cfc_contracted_vs_pairwise`, which
   !> compares two independent derivation ladders rather than a difference
   !> quotient.
   real(wp), parameter :: CFC_ABS_THR = 1.0e-8_wp
   real(wp), parameter :: CFC_REL_THR = 1.0e-7_wp

   !> Sampling points per structure in the finite-difference drivers.
   !>
   !> CFC gets fewer than SvdW purely because of cost: one CFC evaluation runs an
   !> `O(n_active^2)` pair kernel where SvdW runs `O(n_active)` power sums, which
   !> makes a CFC point roughly 300x a SvdW one. The two counts buy different
   !> amounts of the *same* thing -- the points are independent samples of one
   !> pointwise identity, so the fourth through seventh add far less than the
   !> first three. The diversity that actually matters (five chemistries, five
   !> sizes, and for SvdW the 5x2 blend/gamma sweep) is untouched.
   !>
   !> Safe to tune: `get_test_points` seeds its sampler from `mol%nat` alone, so a
   !> smaller count returns a strict *prefix* of the same sequence rather than a
   !> different draw. Shrinking it can only remove points, never move them.
   integer, parameter :: n_points_svdw = 7
   integer, parameter :: n_points_cfc = 1

   integer, parameter :: n_svdw_blends = 5
   integer, parameter :: n_svdw_gammas = 2
   real(wp), parameter :: svdw_blend_k_values(n_svdw_blends) = [0.1_wp, 0.5_wp, 1.0_wp, 2.0_wp, 3.0_wp]
   real(wp), parameter :: svdw_gamma_values(n_svdw_gammas) = [0.0_wp, 0.7_wp]

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
                  new_unittest("svdw_screening_separation", test_svdw_screening_separation), &
                  new_unittest("cfc_screening_separation", test_cfc_screening_separation), &
                  !> LSFs
                  new_unittest("svdw_f1_r_fd", test_svdw_f1_r_fd), &
                  new_unittest("cfc_f1_r_fd", test_cfc_f1_r_fd), &
                  new_unittest("svdw_f2_rr_fd", test_svdw_f2_rr_fd), &
                  new_unittest("cfc_f2_rr_fd", test_cfc_f2_rr_fd), &
                  new_unittest("svdw_f3_rrr_fd", test_svdw_f3_rrr_fd), &
                  new_unittest("cfc_f3_rrr_fd", test_cfc_f3_rrr_fd), &
                  new_unittest("svdw_f1_ra_fd", test_svdw_f1_rA_fd), &
                  new_unittest("cfc_f1_ra_fd", test_cfc_f1_rA_fd), &
                  new_unittest("svdw_f2_r_ra_fd", test_svdw_f2_r_rA_fd), &
                  new_unittest("cfc_f2_r_ra_fd", test_cfc_f2_r_rA_fd), &
                  new_unittest("svdw_f3_rr_ra_fd", test_svdw_f3_rr_rA_fd), &
                  new_unittest("cfc_f3_rr_ra_fd", test_cfc_f3_rr_rA_fd), &
                  !> Radius derivatives
                  new_unittest("svdw_f1_rad_fd", test_svdw_f1_rad_fd), &
                  new_unittest("cfc_f1_rad_fd", test_cfc_f1_rad_fd), &
                  new_unittest("svdw_f2_r_rad_fd", test_svdw_f2_r_rad_fd), &
                  new_unittest("cfc_f2_r_rad_fd", test_cfc_f2_r_rad_fd), &
                  new_unittest("svdw_f3_rr_rad_fd", test_svdw_f3_rr_rad_fd), &
                  new_unittest("cfc_f3_rr_rad_fd", test_cfc_f3_rr_rad_fd), &
                  new_unittest("svdw_radrad_fd", test_svdw_radrad_fd), &
                  new_unittest("cfc_radrad_fd", test_cfc_radrad_fd), &
                  new_unittest("svdw_ra_rad_fd", test_svdw_rA_rad_fd), &
                  new_unittest("cfc_ra_rad_fd", test_cfc_rA_rad_fd), &
                  new_unittest("svdw_hvp_rad_fd", test_svdw_hvp_rad_fd), &
                  new_unittest("cfc_hvp_rad_fd", test_cfc_hvp_rad_fd), &
                  new_unittest("svdw_hvp_ra_joint_fd", test_svdw_hvp_rA_joint_fd), &
                  new_unittest("cfc_hvp_ra_joint_fd", test_cfc_hvp_rA_joint_fd), &
                  new_unittest("svdw_radius_pairwise", test_svdw_radius_pairwise), &
                  new_unittest("cfc_radius_pairwise", test_cfc_radius_pairwise), &
                  new_unittest("svdw_empty_active", test_svdw_empty_active), &
                  new_unittest("cfc_empty_active", test_cfc_empty_active), &
                  !> Consistency checks (with+without Hessian)
                  new_unittest("svdw_f012_hessfree_value_grad", test_svdw_f012_hessfree), &
                  new_unittest("cfc_f012_hessfree_value_grad", test_cfc_f012_hessfree), &
                  !> SvdW-only
                  new_unittest("svdw_f2_rarb", test_svdw_f2_rArB), &
                  new_unittest("svdw_f3_r_rarb", test_svdw_f3_r_rArB), &
                  new_unittest("svdw_f4_rrrr", test_svdw_f4_rrrr), &
                  new_unittest("svdw_f4_rrr_ra", test_svdw_f4_rrr_rA), &
                  new_unittest("svdw_f4_rr_rarb", test_svdw_f4_rr_rArB), &
                  new_unittest("svdw_normalized_f1_ra", test_svdw_normalized_f1_rA), &
                  new_unittest("svdw_contracted_vs_pairwise", test_svdw_contracted), &
                  !> CFC-only
                  new_unittest("cfc_f2_rarb", test_cfc_f2_rArB), &
                  new_unittest("cfc_f3_r_rarb", test_cfc_f3_r_rArB), &
                  new_unittest("cfc_f4_rrrr", test_cfc_f4_rrrr), &
                  new_unittest("cfc_f4_rrr_ra", test_cfc_f4_rrr_rA), &
                  new_unittest("cfc_f4_rr_rarb", test_cfc_f4_rr_rArB), &
                  new_unittest("cfc_normalized_f1_ra", test_cfc_normalized_f1_rA), &
                  new_unittest("cfc_contracted_vs_pairwise", test_cfc_contracted), &
                  new_unittest("svdw_point_on_nucleus_is_finite", test_svdw_on_nucleus), &
                  new_unittest("svdw_body_order_scaling", test_svdw_body_order_scaling), &
                  new_unittest("svdw_exclusion_radius_never_overclaims", &
                               test_svdw_exclusion_radius), &
                  new_unittest("svdw_exclusion_radius_gated_on_blends", &
                               test_svdw_exclusion_radius_gate), &
                  new_unittest("cfc_supplies_no_exclusion_radius", &
                               test_cfc_exclusion_radius) &
                  ]
   end subroutine collect_cavity_drop_lsf

   !* ================================================================================= *!
   !*                              Local helpers                                        *!
   !* ================================================================================= *!

   !> Allocate a fresh polymorphic LSF of the requested concrete kind and
   !> bind it to the given molecule. Optional `blend_k` / `blend_2b` /
   !> `blend_3b` apply only when `kind == "svdw"` (CFC has no equivalent knob).
   !> Optional `screening_threshold` sets the inherited base-type field
   !> *before* `new`/`update` so the SSD system picks it up; default is
   !> 0 (no screening), matching the rest of the suite.
   !>
   !> @param[out] lsf                 polymorphic LSF allocatable
   !> @param[in]  mol                 molecular structure to bind
   !> @param[in]  radii               per-atom radii (size mol%nat)
   !> @param[in]  max_deriv           highest spatial derivative order to enable
   !> @param[in]  kind                "svdw" or "cfc"
   !> @param[in]  blend_k             optional svdw blend sharpness override
   !> @param[in]  blend_2b            optional svdw 2-body weight override
   !> @param[in]  blend_3b            optional svdw 3-body weight override
   !> @param[in]  screening_threshold optional screening cutoff (default 0)
   subroutine init_lsf(lsf, mol, radii, max_deriv, kind, blend_k, blend_2b, &
                       blend_3b, screening_threshold)
      class(moist_cavity_drop_lsf_type), allocatable, intent(out) :: lsf
      !> Molecular structure to bind to the LSF
      type(structure_type), intent(in) :: mol
      !> Per-atom radii (size mol%nat)
      real(wp), intent(in) :: radii(:)
      !> Highest spatial derivative order to enable
      integer, intent(in) :: max_deriv
      !> Concrete kind selector: "svdw" or "cfc"
      character(len=*), intent(in) :: kind
      !> Optional SvdW blend_k override (ignored for CFC)
      real(wp), intent(in), optional :: blend_k
      !> Optional SvdW blend_2b override (ignored for CFC)
      real(wp), intent(in), optional :: blend_2b
      !> Optional SvdW blend_3b override (ignored for CFC)
      real(wp), intent(in), optional :: blend_3b
      !> Optional screening threshold (default 0, no screening)
      real(wp), intent(in), optional :: screening_threshold

      type(moist_cavity_drop_lsf_svdw_type) :: tmp_svdw
      type(moist_cavity_drop_lsf_cfc_type)  :: tmp_cfc
      real(wp) :: thr

      thr = 0.0_wp
      if (present(screening_threshold)) thr = screening_threshold

      select case (kind)
      case (kind_svdw)
         tmp_svdw%screening_threshold = thr
         call tmp_svdw%new(blend_k=blend_k, blend_2b=blend_2b, blend_3b=blend_3b)
         call tmp_svdw%update(mol, radii)
         call tmp_svdw%set_max_deriv(max_deriv)
         allocate (lsf, source=tmp_svdw)
      case (kind_cfc)
         tmp_cfc%screening_threshold = thr
         call tmp_cfc%new()
         call tmp_cfc%update(mol, radii)
         call tmp_cfc%set_max_deriv(max_deriv)
         allocate (lsf, source=tmp_cfc)
      case default
         error stop "init_lsf: unknown kind '"//kind//"'"
      end select
   end subroutine init_lsf

   !> Refresh only the geometry on the underlying concrete (no full LSF
   !> re-init). Used by nuclear-derivative FDs to perturb atom positions
   !> without wiping `max_deriv` or other concrete caches.
   !>
   !> `set_centers` is part of the abstract base API, so no `select type`
   !> dispatch is needed; `radii` is accepted only to keep the call sites
   !> reading as "move these atoms, keep these radii".
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

   !> Prepare the LSF at a displaced geometry, failing the test if it cannot.
   !>
   !> The prepare at the *base* geometry is checked inline in every driver,
   !> because it sits next to the active-list identity check that follows it.
   !> The displaced prepares inside the FD loops have no such context, and an
   !> unreported failure there is the worst kind: `f0` and friends would return
   !> the values of the previous point, and the difference quotient built from
   !> them looks like an ordinary numerical disagreement rather than a broken
   !> setup. A `.true.` return is the caller's signal to return at once.
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

   !> SvdW dispatch for the sign-convention check.
   subroutine test_svdw_sign_convention(error)
      type(error_type), allocatable, intent(out) :: error
      call run_sign_convention(error, kind_svdw)
   end subroutine test_svdw_sign_convention

   !> CFC dispatch for the sign-convention check.
   subroutine test_cfc_sign_convention(error)
      type(error_type), allocatable, intent(out) :: error
      call run_sign_convention(error, kind_cfc)
   end subroutine test_cfc_sign_convention

   !> Verifies LSF is negative deep inside an atom (high PD) and positive
   !> outside (low PD). Uses LiH explicitly - a small heteronuclear
   !> diatomic with a clear inside/outside.
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

      ! Inside: lies on atom 1, well inside the cavity
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

      ! Outside: well beyond the molecular bounding box along +x
      outside_pt = [maxval(mol%xyz(1, :)) + 20.0_wp, 0.0_wp, 0.0_wp]
      call lsf%prepare(outside_pt, lsf_err)
      call lsf%f0(val_outside)
      if (val_outside <= 0.0_wp) then
         call test_failed(error, "Expected lsf0 > 0 well outside the cavity")
         return
      end if
   end subroutine run_sign_convention

   !> SvdW dispatch for the molecule-pair separation check.
   subroutine test_svdw_screening_separation(error)
      type(error_type), allocatable, intent(out) :: error
      call run_screening_separation(error, kind_svdw)
   end subroutine test_svdw_screening_separation

   !> CFC dispatch for the molecule-pair separation check.
   subroutine test_cfc_screening_separation(error)
      type(error_type), allocatable, intent(out) :: error
      call run_screening_separation(error, kind_cfc)
   end subroutine test_cfc_screening_separation

   !> Verify screening in the regime it is actually used: a stationary centers
   !> molecule with a copy of itself pulled away from it, sampling the LSF only
   !> at points near the centers molecule (never far out in vacuum).
   !>
   !> Unlike `run_neighbor_cutoff`, which marches the evaluation point outward
   !> until conditioning degrades, here the points stay near the stationary
   !> centers molecule (well-conditioned, val ~ 0) and screening is driven by
   !> the *separation* of the moving copy: when adjacent its atoms are within
   !> the SSD cutoff and contribute (screened == unscreened); as it recedes
   !> past the cutoff (~25-30 bohr for the production threshold) its atoms are
   !> dropped, but by then their contribution is ~threshold, so the screened
   !> LSF still tracks the unscreened reference within ~ntot^2 * X.
   !>
   !> The centers and moving molecules are the same structure (a homo-dimer
   !> separation), joined with [[BuildSuperStructure]] (centers first, moving
   !> copy appended), so the moving atoms carry user ids > `nc`; that lets the
   !> realism guard confirm the moving copy contributes when adjacent and is
   !> fully screened away when far - the contribute->screen transition that
   !> makes this a non-vacuous test of production-style screening.
   !>
   !> The SvdW blending weights are pinned to the legacy set
   !> (`svdw_legacy_blend_*`) rather than the shipped defaults: the bound
   !> `ntot**2 * thr` is derived for unit 2b/3b weights, and this test is about
   !> the SSD screening machinery, not about whichever blending the production
   !> cavity currently ships.
   subroutine run_screening_separation(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: center_mol, moving_mol, super_mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_ref, lsf_scr
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: thr, gap, shift, val_ref, val_scr, diff, bound, perc
      integer :: icase, ithr, istep, ip, nc, ntot, min_nact, n_mov_pts

      real(wp), parameter :: separation_thresholds(6) = [ &
                             1.0e-8_wp, 1.0E-9_wp, 1.0e-10_wp, 1.0e-11_wp, 1.0e-12_wp, 1.0e-13_wp]

      integer, parameter :: separation_n_steps = 100
      real(wp), parameter :: separation_gap_step = 0.25_wp
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols, 10)

      do icase = 1, size(mols)
         ! Center and moving molecules; the centers stays fixed, the other is translated
         center_mol = mols(icase)
         call center_at_origin(center_mol)
         nc = center_mol%nat
         call get_test_points(center_mol, points, 12)

         do ithr = 1, size(separation_thresholds)
            thr = separation_thresholds(ithr)

            do istep = 0, separation_n_steps
               gap = real(istep, wp)*separation_gap_step

               ! Create superstructure from two molecules
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

                  ! The three body term introduces the largest errors prop. to thr * natoms^2
                  ! This is because the three body term is inside ntot^2 terms; thus we use this as the thr
                  call check(error, val_scr, val_ref, thr=real(ntot, wp)**2*thr, &
                             more="screened lsf diverged from unscreened during separation (lower than natoms^2*thr)")
                  if (allocated(error)) return

               end do
               deallocate (lsf_ref, lsf_scr)
            end do
         end do
      end do

   end subroutine run_screening_separation

   !> True if any currently-active atom belongs to the moving copy (user id
   !> > `nc`, the centers molecule's atom count), i.e. the moving molecule
   !> still contributes after screening. Reflects the most recent `prepare`.
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

   !> Assert that the active list is the unscreened identity
   !>
   !> Every nuclear- and radius-derivative accessor returns *active*-indexed
   !> results: slot `i` belongs to `active_atom(i)`, and slots above
   !> `active_count()` are never written. The FD drivers below allocate their
   !> analytic buffers at `nat` and index them with user-space atom ids, which is
   !> only correct while `init_lsf` leaves screening off. Assert that rather than
   !> assume it -- a screened fixture would otherwise compare uninitialized
   !> memory against a finite difference instead of failing.
   !>
   !> Call once after the first `prepare` of each geometry; the active list does
   !> not change under the FD displacements used here.
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

   !> SvdW dispatch for the spatial gradient FD check.
   subroutine test_svdw_f1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_r_fd(error, kind_svdw)
   end subroutine test_svdw_f1_r_fd

   !> CFC dispatch for the spatial gradient FD check.
   subroutine test_cfc_f1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_r_fd(error, kind_cfc)
   end subroutine test_cfc_f1_r_fd

   !> grad vs 4-point central FD of f0.
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
                     numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, STEP_SIZE)
                  end do
                  do i = 1, ndim
                     call check(error, analytic(i), numeric(i), &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
               deallocate (lsf)
            end do
         end do
      end do
   end subroutine run_f1_r_fd

   !> SvdW dispatch for the spatial Hessian FD check.
   subroutine test_svdw_f2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rr_fd(error, kind_svdw)
   end subroutine test_svdw_f2_rr_fd

   !> CFC dispatch for the spatial Hessian FD check.
   subroutine test_cfc_f2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_rr_fd(error, kind_cfc)
   end subroutine test_cfc_f2_rr_fd

   !> Hessian vs 4-point FD of grad.
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
                        numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), STEP_SIZE)
                     end do
                  end do
                  do j = 1, ndim
                     do i = 1, ndim
                        call check(error, analytic(i, j), numeric(i, j), &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
               deallocate (lsf)
            end do
         end do
      end do
   end subroutine run_f2_rr_fd

   !> SvdW dispatch for the spatial third-derivative FD check.
   subroutine test_svdw_f3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rrr_fd(error, kind_svdw)
   end subroutine test_svdw_f3_rrr_fd

   !> CFC dispatch for the spatial third-derivative FD check.
   subroutine test_cfc_f3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rrr_fd(error, kind_cfc)
   end subroutine test_cfc_f3_rrr_fd

   !> Third spatial derivative vs FD of Hessian (both pulled from
   !> f3_rrr so the FD and analytic branches share the same
   !> internal code path; matters at -O3).
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
                           numeric(i, j, axis) = fd4_scalar( &
                                                 hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), eps)
                        end do
                     end do
                  end do
                  do k = 1, ndim
                     do j = 1, ndim
                        do i = 1, ndim
                           call check(error, analytic(i, j, k), numeric(i, j, k), &
                                      thr_abs=ABS_THR, thr_rel=REL_THR)
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

   !> SvdW dispatch for the nuclear gradient FD check.
   subroutine test_svdw_f1_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rA_fd(error, kind_svdw)
   end subroutine test_svdw_f1_rA_fd

   !> CFC dispatch for the nuclear gradient FD check.
   subroutine test_cfc_f1_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rA_fd(error, kind_cfc)
   end subroutine test_cfc_f1_rA_fd

   !> Nuclear-position gradient vs FD of f0.
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
                        numeric = fd4_scalar(f_pp, f_p, f_m, f_mm, eps)
                        call check(error, analytic(axis, atom), numeric, &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
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

   !> SvdW dispatch for the mixed spatial-nuclear FD check.
   subroutine test_svdw_f2_r_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rA_fd(error, kind_svdw)
   end subroutine test_svdw_f2_r_rA_fd

   !> CFC dispatch for the mixed spatial-nuclear FD check.
   subroutine test_cfc_f2_r_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rA_fd(error, kind_cfc)
   end subroutine test_cfc_f2_r_rA_fd

   !> Mixed spatial-nuclear second derivative vs FD of spatial grad.
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
                           numeric = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), eps)
                           call check(error, analytic(i, axis, atom), numeric, &
                                      thr_abs=ABS_THR, thr_rel=REL_THR)
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

   !> SvdW dispatch for the mixed third FD check.
   subroutine test_svdw_f3_rr_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rA_fd(error, kind_svdw)
   end subroutine test_svdw_f3_rr_rA_fd

   !> CFC dispatch for the mixed third FD check.
   subroutine test_cfc_f3_rr_rA_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rA_fd(error, kind_cfc)
   end subroutine test_cfc_f3_rr_rA_fd

   !> Mixed third derivative (Hess w.r.t. spatial coords, grad w.r.t.
   !> nuclei) vs FD of spatial Hessian.
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
                              numeric = fd4_scalar(hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), eps)
                              call check(error, analytic(i, j, axis, atom), numeric, &
                                         thr_abs=ABS_THR, thr_rel=REL_THR)
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
   ! `R_a` is the *radius* of an active atom, not its position, so `f3_rr_rad`
   ! and its lower orders carry no trailing derivative index: the rank equals the
   ! spatial order. These three checks differ from the `_rA` block above in one
   ! further respect. A nuclear FD can move atoms with `set_centers`, which only
   ! rebuilds the spatial sort; a radius FD cannot, because the radii are baked
   ! into every concrete's per-atom cache. `refresh_radii` therefore goes through
   ! the full `update` and restores `max_deriv` afterwards.
   !
   ! Both concretes are registered. Their radius channels are built on opposite
   ! structures -- SvdW's `d/dR_a` is a diagonal rescaling of the atom's own
   ! tensors, CFC's is a genuine extra derivative of a two-center kernel -- so
   ! running the same FD against both is worth more than running it twice.

   !> SvdW dispatch for the radius gradient FD check.
   subroutine test_svdw_f1_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rad_fd(error, kind_svdw)
   end subroutine test_svdw_f1_rad_fd

   !> CFC dispatch for the radius gradient FD check.
   subroutine test_cfc_f1_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f1_rad_fd(error, kind_cfc)
   end subroutine test_cfc_f1_rad_fd

   !> Radius gradient dS/dR_a vs FD of f0.
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
                     numeric = fd4_scalar(f_pp, f_p, f_m, f_mm, eps)
                     call check(error, analytic(atom), numeric, &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
               deallocate (analytic, lsf)
            end do
         end do
         deallocate (radii_local, dummy_2nd)
      end do
   end subroutine run_f1_rad_fd

   !> SvdW dispatch for the mixed spatial-radius FD check.
   subroutine test_svdw_f2_r_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rad_fd(error, kind_svdw)
   end subroutine test_svdw_f2_r_rad_fd

   !> CFC dispatch for the mixed spatial-radius FD check.
   subroutine test_cfc_f2_r_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f2_r_rad_fd(error, kind_cfc)
   end subroutine test_cfc_f2_r_rad_fd

   !> Mixed spatial-radius second derivative vs FD of the spatial gradient.
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
                        numeric = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), eps)
                        call check(error, analytic(i, atom), numeric, &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
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

   !> SvdW dispatch for the mixed third radius FD check.
   subroutine test_svdw_f3_rr_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rad_fd(error, kind_svdw)
   end subroutine test_svdw_f3_rr_rad_fd

   !> CFC dispatch for the mixed third radius FD check.
   subroutine test_cfc_f3_rr_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f3_rr_rad_fd(error, kind_cfc)
   end subroutine test_cfc_f3_rr_rad_fd

   !> FD thresholds of the deep radius checks: CFC gets the same headroom the
   !> CFC high-order block gets, and for the same reason.
   !>
   !> CFC's radius derivatives grow faster than its nuclear ones. Both bring
   !> down a factor `|a1| s_a / R_a ~ 7.5`, but `R_a` sits in a *denominator*, so
   !> each further `d/dR_a` also differentiates `1/R_a` and picks up the extra
   !> quotient terms; the nuclear derivative has no such chain. That is why
   !> `f3_rr_rA` clears the tight thresholds for CFC at this step size and
   !> `f3_rr_rad` does not.
   !>
   !> Measured, not guessed: at `STEP_SIZE` the three deep radius checks miss the
   !> tight thresholds by 2-3x, and at half the step all three clear them, which
   !> is the h**4 the stencil promises and not a defect. `f1_rad` and `f2_r_rad`
   !> stay on the tight thresholds -- they pass there, and a check that can be
   !> strict should be.
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
      case (kind_cfc)
         thr_abs = CFC_ABS_THR
         thr_rel = CFC_REL_THR
      case default
         thr_abs = ABS_THR
         thr_rel = REL_THR
      end select
   end subroutine radius_fd_thr

   !> Mixed third derivative d^3S/(dr^2 dR_a) vs FD of the spatial Hessian.
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
                           numeric = fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), &
                                                h_mm(i, j), eps)
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
   ! [[run_radius_contracted_vs_pairwise]] ties `f*_radrad` and `f*_rA_rad` to the
   ! two `hvp_*` rows exactly, but only through one fixed direction pair, and the
   ! `hvp_*` FD checks differentiate along that same pair. An error shared between
   ! the two families that cancels in that one contraction survives all three. The
   ! two drivers below close the hole by differencing element by element: a single
   ! radius is perturbed on its own, and the entire first-derivative ladder of the
   ! other slot is differenced against it.
   !
   ! Perturbing the radius of atom `b` and differencing
   !   * the `f*_rad` ladder gives the `f*_radrad` column `(.., :, b)`
   !   * the `f*_rA`  ladder gives the `f*_rA_rad` column `(.., :, b)`
   ! The second is what pins the *order* of the two atom slots. That family is not
   ! symmetric in them -- the first index carries a position derivative, the second
   ! a radius one -- so a transposed implementation can still satisfy the
   ! contraction identity and fails only here.

   !> SvdW dispatch for the two-radius block FD check.
   subroutine test_svdw_radrad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radrad_fd(error, kind_svdw)
   end subroutine test_svdw_radrad_fd

   !> CFC dispatch for the two-radius block FD check.
   subroutine test_cfc_radrad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radrad_fd(error, kind_cfc)
   end subroutine test_cfc_radrad_fd

   !> `f{2,3,4}_*_radrad` vs a one-radius FD of the whole `f*_rad` ladder.
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
                        numeric = fd4_scalar(l1(other, 1), l1(other, 2), &
                                             l1(other, 3), l1(other, 4), eps)
                        call check(error, rr2(other, atom), numeric, &
                                   thr_abs=thr_abs, thr_rel=thr_rel)
                        if (allocated(error)) return
                        do i = 1, ndim
                           numeric = fd4_scalar(l2(i, other, 1), l2(i, other, 2), &
                                                l2(i, other, 3), l2(i, other, 4), eps)
                           call check(error, rr3(i, other, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                           do j = 1, ndim
                              numeric = fd4_scalar(l3(i, j, other, 1), l3(i, j, other, 2), &
                                                   l3(i, j, other, 3), l3(i, j, other, 4), eps)
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

   !> SvdW dispatch for the nuclear-radius block FD check.
   subroutine test_svdw_rA_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rA_rad_fd(error, kind_svdw)
   end subroutine test_svdw_rA_rad_fd

   !> CFC dispatch for the nuclear-radius block FD check.
   subroutine test_cfc_rA_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rA_rad_fd(error, kind_cfc)
   end subroutine test_cfc_rA_rad_fd

   !> `f{2,3,4}_*_rA_rad` vs a one-radius FD of the whole `f*_rA` ladder.
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
                           numeric = fd4_scalar(n1(s_ax, other, 1), n1(s_ax, other, 2), &
                                                n1(s_ax, other, 3), n1(s_ax, other, 4), eps)
                           call check(error, mx2(s_ax, other, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                           do i = 1, ndim
                              numeric = fd4_scalar(n2(i, s_ax, other, 1), &
                                                   n2(i, s_ax, other, 2), &
                                                   n2(i, s_ax, other, 3), &
                                                   n2(i, s_ax, other, 4), eps)
                              call check(error, mx3(i, s_ax, other, atom), numeric, &
                                         thr_abs=thr_abs, thr_rel=thr_rel)
                              if (allocated(error)) return
                              do j = 1, ndim
                                 numeric = fd4_scalar(n3(i, j, s_ax, other, 1), &
                                                      n3(i, j, s_ax, other, 2), &
                                                      n3(i, j, s_ax, other, 3), &
                                                      n3(i, j, s_ax, other, 4), eps)
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
   ! With `R = R(X)` the contraction direction of the level-set HVP is the joint
   ! pair `(v_B, vr_B)`. The two checks below step the whole molecule along that
   ! joint direction at once -- centers by `eps*v`, radii by `eps*vr` -- and
   ! difference the corresponding first-derivative ladder. That is a genuine
   ! end-to-end test of the composition: it fails if either half of the
   ! contraction is missing, and it fails if the two halves are mixed up, because
   ! `v` and `vr` are chosen unrelated to each other.
   !
   ! Two rows are covered:
   !   * `hvp_f*_rad`  radius row, `d/dR_a` retained  -- FD of `f*_rad`
   !   * `hvp_f*_rA`   nuclear row with `vrad` passed -- FD of `f*_rA`
   ! Between them they exercise all four blocks of the joint Hessian, including
   ! the radius-radius coupling between different atoms, which no per-atom
   ! quantity can produce.

   !> Deterministic, non-symmetric joint direction field, normalized to max 1.
   !>
   !> The normalization is what makes the FD checks that use this field hold at
   !> the shared `ABS_THR`/`REL_THR`. The direction sets the *effective* step:
   !> stepping by `eps*v` with `max|v| = m` gives a four-point central difference
   !> a truncation error scaling like `(m*eps)**4`, so an unnormalized field whose
   !> entries grow with the atom index (m ~ 4 at nat = 16) costs two orders of
   !> magnitude of accuracy for nothing. Scaling to `max|v| = 1` is measured, not
   !> guessed: without it these checks land at ~1e-9 relative, just above
   !> `REL_THR`, and the miss grows with molecule size.
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
         ! Deliberately unrelated to `v`: a real radius model ties the two, and a
         ! fixture that respected that tie could hide a term mixing them.
         vrad(ia) = 0.19_wp*real(mod(ia, 4) + 1, wp) - 0.43_wp
      end do

      scal = max(maxval(abs(v)), maxval(abs(vrad)))
      if (scal > 0.0_wp) then
         v = v/scal
         vrad = vrad/scal
      end if
   end subroutine joint_direction

   !> SvdW dispatch for the uncontracted radius-block cross-check.
   subroutine test_svdw_radius_pairwise(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radius_contracted_vs_pairwise(error, kind_svdw)
   end subroutine test_svdw_radius_pairwise

   !> CFC dispatch for the uncontracted radius-block cross-check.
   subroutine test_cfc_radius_pairwise(error)
      type(error_type), allocatable, intent(out) :: error
      call run_radius_contracted_vs_pairwise(error, kind_cfc)
   end subroutine test_cfc_radius_pairwise

   !> Contract the uncontracted radius blocks and recover both `hvp_*` rows.
   !>
   !> No finite differences anywhere: with `R = R(X)` the joint Hessian-vector
   !> product is an exact contraction of the four uncontracted blocks, so this
   !> holds to roundoff and is a far sharper gate than the FD tests above. It is
   !> the only check that pins `f*_radrad` and `f*_rA_rad` at `REL_THR` rather
   !> than at the loosened radius-FD thresholds (see [[radius_fd_thr]]).
   !>
   !> Both rows are verified, which matters because the nuclear-radius block is *not*
   !> symmetric in its two atom indices:
   !>
   !>    hvp_f{n}_rad(A)  = sum_B [ v_B . f_rA_rad(:, B, A) + vr_B f_radrad(A, B) ]
   !>    hvp_f{n}_rA(A,s) = sum_B [ v_B . f_rArB(s, A, :, B) + vr_B f_rA_rad(s, A, B) ]
   !>
   !> The first reads the nuclear-radius block with the position index on B, the second
   !> with it on A. A transposed implementation passes neither.
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

            !* ------------------------- the radius row -------------------------- *!
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

            !* --- the position row, whose radius half is the nuclear-radius block again --- *!
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

   !> SvdW dispatch for the radius-row HVP FD check.
   subroutine test_svdw_hvp_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rad_fd(error, kind_svdw)
   end subroutine test_svdw_hvp_rad_fd

   !> CFC dispatch for the radius-row HVP FD check.
   subroutine test_cfc_hvp_rad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rad_fd(error, kind_cfc)
   end subroutine test_cfc_hvp_rad_fd

   !> Radius row `hvp_f*_rad` vs joint-direction FD of the `f*_rad` ladder.
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
                     numeric = fd4_scalar(l1(atom, 1), l1(atom, 2), l1(atom, 3), &
                                          l1(atom, 4), eps)
                     call check(error, a1(atom), numeric, &
                                thr_abs=thr_abs, thr_rel=thr_rel)
                     if (allocated(error)) return
                     do i = 1, ndim
                        numeric = fd4_scalar(l2(i, atom, 1), l2(i, atom, 2), &
                                             l2(i, atom, 3), l2(i, atom, 4), eps)
                        call check(error, a2(i, atom), numeric, &
                                   thr_abs=thr_abs, thr_rel=thr_rel)
                        if (allocated(error)) return
                        do j = 1, ndim
                           numeric = fd4_scalar(l3(i, j, atom, 1), l3(i, j, atom, 2), &
                                                l3(i, j, atom, 3), l3(i, j, atom, 4), eps)
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

   !> SvdW dispatch for the joint nuclear-row HVP FD check.
   subroutine test_svdw_hvp_rA_joint_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rA_joint_fd(error, kind_svdw)
   end subroutine test_svdw_hvp_rA_joint_fd

   !> CFC dispatch for the joint nuclear-row HVP FD check.
   subroutine test_cfc_hvp_rA_joint_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_hvp_rA_joint_fd(error, kind_cfc)
   end subroutine test_cfc_hvp_rA_joint_fd

   !> Nuclear row `hvp_f*_rA(v, res, vrad)` vs joint-direction FD of `f*_rA`.
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
                        numeric = fd4_scalar(n1(s_ax, atom, 1), n1(s_ax, atom, 2), &
                                             n1(s_ax, atom, 3), n1(s_ax, atom, 4), eps)
                        call check(error, h1(s_ax, atom), numeric, &
                                   thr_abs=thr_abs, thr_rel=thr_rel)
                        if (allocated(error)) return
                        do i = 1, ndim
                           numeric = fd4_scalar(n2(i, s_ax, atom, 1), n2(i, s_ax, atom, 2), &
                                                n2(i, s_ax, atom, 3), n2(i, s_ax, atom, 4), eps)
                           call check(error, h2(i, s_ax, atom), numeric, &
                                      thr_abs=thr_abs, thr_rel=thr_rel)
                           if (allocated(error)) return
                           do j = 1, ndim
                              numeric = fd4_scalar(n3(i, j, s_ax, atom, 1), &
                                                   n3(i, j, s_ax, atom, 2), &
                                                   n3(i, j, s_ax, atom, 3), &
                                                   n3(i, j, s_ax, atom, 4), eps)
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
   ! Screening is allowed to reject every atom -- a point far enough outside the
   ! molecule contributes nothing -- and that leaves the accessors with no owned
   ! slot to write. Their normal contract, "writes the first `active_count()`
   ! slots and leaves the rest alone", degenerates there, and since every result
   ! is `intent(out)` a bare early return would hand back a buffer the caller is
   ! not allowed to read. The check below is the only one in this file that
   ! inspects a *whole* result buffer rather than its active prefix, so it is
   ! also the only one that can see that difference.
   !
   ! Every buffer is poisoned first. Poison surviving the call is exactly the
   ! symptom: it means the accessor returned without defining its result.

   !> SvdW dispatch for the empty-active-list check.
   subroutine test_svdw_empty_active(error)
      type(error_type), allocatable, intent(out) :: error
      call run_empty_active(error, kind_svdw)
   end subroutine test_svdw_empty_active

   !> CFC dispatch for the empty-active-list check.
   subroutine test_cfc_empty_active(error)
      type(error_type), allocatable, intent(out) :: error
      call run_empty_active(error, kind_cfc)
   end subroutine test_cfc_empty_active

   !> Every radius-channel accessor must define its result when nothing is active.
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

      !* A nonzero threshold is what lets screening reject everything; at the
      !* default 0 every atom stays a candidate however far away the point is.
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

      !* Uncontracted radius ladder
      g1 = EMPTY_POISON; g2 = EMPTY_POISON; g3 = EMPTY_POISON
      call lsf%f3_rr_rad(g1, g2, g3)
      call check_all_zero(error, reshape(g1, [size(g1)]), "f1_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(g2, [size(g2)]), "f2_r_rad")
      if (allocated(error)) return
      call check_all_zero(error, reshape(g3, [size(g3)]), "f3_rr_rad")
      if (allocated(error)) return

      !* Uncontracted two-radius and nuclear-radius blocks
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

      !* Radius row of the joint HVP
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

      !* Nuclear row, both with and without a radius direction: the two take
      !* different branches inside CFC, and only one of them allocates.
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

   !> Fail unless every element of `a` is exactly zero.
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

   !> Rebind the LSF to a jointly perturbed geometry and radius vector.
   !>
   !> `set_centers` alone will not do: the radii move too, and every concrete
   !> bakes those into its per-atom cache at `update` time.
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

   !> Rebind the LSF to a perturbed radius vector.
   !>
   !> The nuclear FDs get away with `set_centers`, which only moves atoms and
   !> rebuilds the spatial sort. Radii cannot be perturbed that way: every
   !> concrete bakes them into its per-atom cache at `update` time, so a radius
   !> FD has to re-run the full `update`. That resets `max_deriv`, hence the
   !> explicit restore -- without it the first `f012_r` after a perturbation
   !> would abort in `require_deriv` rather than return a wrong number.
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

   !> SvdW dispatch for the Hessian-free value+gradient consistency check.
   subroutine test_svdw_f012_hessfree(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f012_hessfree(error, kind_svdw)
   end subroutine test_svdw_f012_hessfree

   !> CFC dispatch for the Hessian-free value+gradient consistency check.
   subroutine test_cfc_f012_hessfree(error)
      type(error_type), allocatable, intent(out) :: error
      call run_f012_hessfree(error, kind_cfc)
   end subroutine test_cfc_f012_hessfree

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
                             blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  call lsf%prepare(points(:, ipt), lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  ! Hessian skipped (no lsf2_rr) vs Hessian computed.
                  !
                  ! Agreement is to a few ULP, not bit-for-bit: the SvdW kernel
                  ! is code-generated with a separate common-subexpression set
                  ! per derivative order, so the order-1 and order-2 branches
                  ! reach the same value and gradient by different (equally
                  ! valid) associations. The measured spread is ~1e-14 relative;
                  ! what this test still pins is that skipping the Hessian does
                  ! not change the *answer*.
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

   !> Pick the (blend, gamma) sweep dimensions for the requested kind.
   !> SvdW iterates over the parameter grid; CFC runs a single configuration.
   !> Sampling points per structure for a kind-dispatched driver
   !>
   !> See the [[n_points_cfc]] note above for why the two kinds differ.
   pure function fd_points(kind) result(npts)
      !> Concrete kind selector
      character(len=*), intent(in) :: kind
      !> Number of sampling points to request from `get_test_points`
      integer :: npts
      if (kind == kind_cfc) then
         npts = n_points_cfc
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
      if (kind == kind_svdw) then
         nblend = n_svdw_blends
         ngamma = n_svdw_gammas
      else
         nblend = 1
         ngamma = 1
      end if
   end subroutine svdw_sweep_sizes

   !> i-th SvdW blend_k value, or huge() (treated as "absent") for CFC.
   !> Returned as a function so init_lsf's optional argument is only
   !> defined for the SvdW kind.
   pure function svdw_sweep_blend(kind, i) result(val)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: i
      real(wp) :: val
      if (kind == kind_svdw) then
         val = svdw_blend_k_values(i)
      else
         !* Placeholder; CFC call sites pass this via "optional" wrapping below.
         val = 0.0_wp
      end if
   end function svdw_sweep_blend

   !> i-th SvdW gamma (3-body) value, or 0 for CFC.
   pure function svdw_sweep_gamma(kind, i) result(val)
      character(len=*), intent(in) :: kind
      integer, intent(in) :: i
      real(wp) :: val
      if (kind == kind_svdw) then
         val = svdw_gamma_values(i)
      else
         val = 0.0_wp
      end if
   end function svdw_sweep_gamma

   !* ================================================================================= *!
   !*                  SvdW-only extension tests (kept as-is from primitives)           *!
   !* ================================================================================= *!

   !> SvdW nuclear-position second derivative (f2_rArB) vs FD of f1_rA.
   !>
   !> `f2_rArB` is *active-indexed*, like its `f3_r_rArB` and
   !> `f4_rr_rArB` siblings: `analytic(:, iA, :, iB)` is the block for the
   !> atoms `active_atom(iA)` / `active_atom(iB)`, and the result is sized
   !> `(3, n_active, 3, n_active)`. The FD reference below is built per user-space
   !> atom, so the test asserts the two agree slot by slot after translating the
   !> active index; with `screening_threshold = 0` every atom stays active and the
   !> translation is the identity, which the guard below pins down.
   subroutine test_svdw_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iblend, igamma, atomA, axisA, atomB, axisB
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :)
      real(wp), allocatable :: rA_fwd(:, :), rA_fwd2(:, :)
      real(wp), allocatable :: rA_bwd(:, :), rA_bwd2(:, :)
      real(wp), allocatable :: dummy_rr_rA(:, :, :, :)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         !* FD buffers are sized by mol%nat only and reused across the
         !* (iblend, igamma, ipt) sweep.
         if (allocated(rA_fwd)) deallocate (rA_fwd)
         if (allocated(rA_fwd2)) deallocate (rA_fwd2)
         if (allocated(rA_bwd)) deallocate (rA_bwd)
         if (allocated(rA_bwd2)) deallocate (rA_bwd2)
         allocate (rA_fwd(ndim, mol%nat), rA_fwd2(ndim, mol%nat))
         allocate (rA_bwd(ndim, mol%nat), rA_bwd2(ndim, mol%nat))
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         if (allocated(analytic)) deallocate (analytic)
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol%nat))
         allocate (analytic(ndim, mol%nat, ndim, mol%nat))
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(3)
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%set_centers(centers_base)
                  call prim%prepare(point, lsf_err)
                  if (allocated(lsf_err)) then
                     call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                     return
                  end if
                  call prim%f2_rArB(analytic)
                  !* Active-indexed result, so assert the shape and defer the
                  !* identity check to the shared helper before the atom loops
                  !* below index `analytic` with user-space ids.
                  call check(error, size(analytic, 2), prim%active_count())
                  if (allocated(error)) return
                  call check_identity_active(error, prim, mol%nat)
                  if (allocated(error)) return
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
                              numeric = fd4_scalar(rA_fwd2(axisA, atomA), rA_fwd(axisA, atomA), &
                                                   rA_bwd(axisA, atomA), rA_bwd2(axisA, atomA), STEP_SIZE)
                              call check(error, analytic(axisA, atomA, axisB, atomB), numeric, &
                                         thr_abs=ABS_THR, thr_rel=REL_THR)
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
   end subroutine test_svdw_f2_rArB

   !> SvdW mixed third (spatial grad x nuclear Hessian) vs FD of f2_r_rA.
   subroutine test_svdw_f3_r_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iblend, igamma, iA, iB, axisA, axisB, jdir
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :)
      real(wp), allocatable :: r_rA_fwd(:, :, :), r_rA_fwd2(:, :, :)
      real(wp), allocatable :: r_rA_bwd(:, :, :), r_rA_bwd2(:, :, :)
      real(wp), allocatable :: dummy_rA(:, :), dummy_rr_rA(:, :, :, :)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         !* f3_r_rArB declares lsf1_rA and lsf2_r_rA as intent(in)
         !* non-allocatable assumed-shape; passing unallocated allocatables
         !* is undefined behavior, so size the dummies up front. All
         !* buffers depend only on mol%nat, so allocate once per icase.
         if (allocated(dummy_rA)) deallocate (dummy_rA)
         if (allocated(r_rA_fwd)) deallocate (r_rA_fwd)
         if (allocated(r_rA_fwd2)) deallocate (r_rA_fwd2)
         if (allocated(r_rA_bwd)) deallocate (r_rA_bwd)
         if (allocated(r_rA_bwd2)) deallocate (r_rA_bwd2)
         allocate (dummy_rA(ndim, mol%nat))
         allocate (r_rA_fwd(ndim, ndim, mol%nat), r_rA_fwd2(ndim, ndim, mol%nat))
         allocate (r_rA_bwd(ndim, ndim, mol%nat), r_rA_bwd2(ndim, ndim, mol%nat))
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         if (allocated(analytic)) deallocate (analytic)
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol%nat))
         allocate (analytic(ndim, ndim, mol%nat, ndim, mol%nat))
         dummy_rA = 0.0_wp
         r_rA_fwd = 0.0_wp
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(3)
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
                                 numeric = fd4_scalar(r_rA_fwd2(jdir, axisA, iA), &
                                                      r_rA_fwd(jdir, axisA, iA), r_rA_bwd(jdir, axisA, iA), &
                                                      r_rA_bwd2(jdir, axisA, iA), STEP_SIZE)
                                 call check(error, analytic(jdir, axisA, iA, axisB, iB), numeric, &
                                            thr_abs=ABS_THR, thr_rel=REL_THR)
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
   end subroutine test_svdw_f3_r_rArB

   !> SvdW pure spatial fourth derivative vs FD of f3_rrr.
   subroutine test_svdw_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, iblend, igamma, axis, i, j, kk
      real(wp) :: point(ndim), work_point(ndim)
      real(wp) :: t3_fwd(ndim, ndim, ndim), t3_fwd2(ndim, ndim, ndim)
      real(wp) :: t3_bwd(ndim, ndim, ndim), t3_bwd2(ndim, ndim, ndim)
      real(wp) :: analytic(ndim, ndim, ndim, ndim)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(4)
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
                              numeric = fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), &
                                                   t3_bwd(i, j, kk), t3_bwd2(i, j, kk), STEP_SIZE)
                              call check(error, analytic(i, j, kk, axis), numeric, &
                                         thr_abs=ABS_THR, thr_rel=REL_THR)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_svdw_f4_rrrr

   !> SvdW mixed fourth (3 r-axes + 1 R-axis) vs FD of f3_rrr.
   subroutine test_svdw_f4_rrr_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iblend, igamma, atom, axis, i, j, kk
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :)
      real(wp) :: t3_fwd(ndim, ndim, ndim), t3_fwd2(ndim, ndim, ndim)
      real(wp) :: t3_bwd(ndim, ndim, ndim), t3_bwd2(ndim, ndim, ndim)
      real(wp) :: numeric
      type(structure_type) :: mol_shift
      integer, allocatable :: atomic_numbers(:)
      integer :: iat
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate (atomic_numbers)
         allocate (atomic_numbers(mol%nat))
         do iat = 1, mol%nat
            atomic_numbers(iat) = mol%num(mol%id(iat))
         end do
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(4)
               if (allocated(analytic)) deallocate (analytic)
               allocate (analytic(ndim, ndim, ndim, ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%update(mol, radii)
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
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rrr(lsf3_rrr=t3_fwd)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rrr(lsf3_rrr=t3_fwd2)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rrr(lsf3_rrr=t3_bwd)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rrr(lsf3_rrr=t3_bwd2)
                        do i = 1, ndim
                           do j = 1, ndim
                              do kk = 1, ndim
                                 numeric = fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), &
                                                      t3_bwd(i, j, kk), t3_bwd2(i, j, kk), STEP_SIZE)
                                 call check(error, analytic(i, j, kk, axis, atom), numeric, &
                                            thr_abs=ABS_THR, thr_rel=REL_THR)
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
   end subroutine test_svdw_f4_rrr_rA

   !> SvdW mixed fourth (2 r-axes + 2 R-axes) vs FD of f3_rr_rA.
   subroutine test_svdw_f4_rr_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iblend, igamma, iA, iB, axisA, axisB, j, kk
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :, :)
      real(wp), allocatable :: dummy_rA(:, :), dummy_r_rA(:, :, :)
      real(wp), allocatable :: rr_rA_fwd(:, :, :, :), rr_rA_fwd2(:, :, :, :)
      real(wp), allocatable :: rr_rA_bwd(:, :, :, :), rr_rA_bwd2(:, :, :, :)
      real(wp) :: numeric
      type(structure_type) :: mol_shift
      integer, allocatable :: atomic_numbers(:)
      integer :: iat
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate (atomic_numbers)
         allocate (atomic_numbers(mol%nat))
         do iat = 1, mol%nat
            atomic_numbers(iat) = mol%num(mol%id(iat))
         end do
         !* f4_rr_rArB declares lsf1_rA and lsf2_r_rA as intent(in)
         !* non-allocatable assumed-shape; passing unallocated allocatables
         !* is undefined behavior, so size the dummies to the active-atom
         !* count up front. Both depend only on mol%nat.
         if (allocated(dummy_rA)) deallocate (dummy_rA)
         if (allocated(dummy_r_rA)) deallocate (dummy_r_rA)
         allocate (dummy_rA(ndim, mol%nat), dummy_r_rA(ndim, ndim, mol%nat))
         dummy_rA = 0.0_wp
         dummy_r_rA = 0.0_wp
         if (allocated(analytic)) deallocate (analytic)
         if (allocated(rr_rA_fwd)) deallocate (rr_rA_fwd, rr_rA_fwd2, rr_rA_bwd, rr_rA_bwd2)
         allocate (analytic(ndim, ndim, ndim, mol%nat, ndim, mol%nat))
         allocate (rr_rA_fwd(ndim, ndim, ndim, mol%nat), rr_rA_fwd2(ndim, ndim, ndim, mol%nat))
         allocate (rr_rA_bwd(ndim, ndim, ndim, mol%nat), rr_rA_bwd2(ndim, ndim, ndim, mol%nat))
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(4)
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%update(mol, radii)
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
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_fwd)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) + 2.0_wp*STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_fwd2)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) - STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_bwd)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) - 2.0_wp*STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%f3_rr_rA(lsf3_rr_rA=rr_rA_bwd2)
                        do iA = 1, mol%nat
                           do axisA = 1, ndim
                              do j = 1, ndim
                                 do kk = 1, ndim
                                    numeric = fd4_scalar(rr_rA_fwd2(j, kk, axisA, iA), &
                                                         rr_rA_fwd(j, kk, axisA, iA), rr_rA_bwd(j, kk, axisA, iA), &
                                                         rr_rA_bwd2(j, kk, axisA, iA), STEP_SIZE)
                                    call check(error, analytic(j, kk, axisA, iA, axisB, iB), &
                                               numeric, thr_abs=ABS_THR, thr_rel=REL_THR)
                                    if (allocated(error)) return
                                 end do
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
   end subroutine test_svdw_f4_rr_rArB

   !> SvdW normalized-LSF nuclear gradient vs FD of normalized f0.
   subroutine test_svdw_normalized_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      real(wp), allocatable :: deriv_rA(:, :)
      integer  :: icase, ipt, iblend, igamma, atom, axis
      real(wp) :: point(ndim)
      real(wp) :: normalized_val
      real(wp) :: f_forward, f_backward, f_forward2, f_backward2, numeric
      type(structure_type) :: mol_shift
      integer, allocatable :: atomic_numbers(:)
      integer :: iat
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate (atomic_numbers)
         allocate (atomic_numbers(mol%nat))
         do iat = 1, mol%nat
            atomic_numbers(iat) = mol%num(mol%id(iat))
         end do
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(2)
               if (allocated(deriv_rA)) deallocate (deriv_rA)
               allocate (deriv_rA(ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%update(mol, radii)
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
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%normalized_f01_rA(f_forward)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%normalized_f01_rA(f_forward2)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%normalized_f01_rA(f_backward)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*STEP_SIZE
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%set_centers(centers_local)
                        call prim%prepare(point, lsf_err)
                        call prim%normalized_f01_rA(f_backward2)
                        numeric = fd4_scalar(f_forward2, f_forward, f_backward, f_backward2, STEP_SIZE)
                        call check(error, deriv_rA(axis, atom), numeric, &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_svdw_normalized_f1_rA

   !> SvdW body-order weight reduction sanity check.
   !>
   !> Verify that selecting blend_2b=1 (pure pair-mean) collapses to the
   !> average of two atom SSDs, and blend_3b=1 (pure triple-mean) collapses
   !> to the mean of three atom SSDs.
   !* ================================================================================= *!
   !*                    Direction-contracted nuclear derivatives                       *!
   !* ================================================================================= *!

   !> The `tangent_*` / `hvp_*` / `vjp_*` accessors must equal the explicit contraction
   !>
   !> All three families exist only so a caller never has to form the full nuclear
   !> tensors; that is worth nothing unless they agree with them. This contracts
   !> the full tensors by hand against the same direction field and compares.
   !> Cheap and total: it exercises every contracted binding against an
   !> independently derived sibling rather than against a finite difference.
   subroutine test_svdw_contracted(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         nat = mol%nat
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         !* A deterministic, non-symmetric direction field
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

         prim%screening_threshold = 0.0_wp
         call prim%new(blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                       blend_3b=svdw_legacy_blend_3b)
         call prim%update(mol, radii)
         call prim%set_max_deriv(4)

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

            !* --------------------------- tangent_* ---------------------------- *!
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

            !* ----------------------------- hvp_* ------------------------------ *!
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

            !* ---------------------------- vjp_f1_rA --------------------------- *!
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

            !* --------------------------- vjp_f1_rad --------------------------- *!
            call prim%vjp_f1_rad(w0_adj, w1_adj, w2_adj, vjp_rad)
            do iA = 1, prim%active_count()
               acc = w0_adj*rad1(iA) &
                     + dot_product(w1_adj, rad2(:, iA)) &
                     + sum(w2_adj*rad3(:, :, iA))
               call check(error, vjp_rad(iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_svdw_contracted

   !* ================================================================================= *!
   !*                              CFC-only extension tests                             *!
   !* ================================================================================= *!
   !
   ! The CFC twins of the SvdW-only block above. They are separate subroutines
   ! rather than shared dispatches because the SvdW originals sweep the
   ! `blend_k` x `gamma` parameter grid, which CFC has no analogue of: its four
   ! shape parameters are compiled-in constants. Everything else -- fixtures,
   ! step size, FD stencil, index conventions -- is identical, so the two blocks
   ! can be diffed against each other.
   !
   ! CFC's exponents (a1 = -15, a2 = -9) make its high-order derivatives an order
   ! of magnitude stiffer than SvdW's, so the fourth-order finite differences
   ! carry visibly more truncation error at the shared step size. That is what
   ! `CFC_ABS_THR` / `CFC_REL_THR` are for; they are FD-stencil headroom, not a
   ! statement about the analytic derivatives, which the exact
   ! `cfc_contracted_vs_pairwise` cross-check below pins at roundoff.

   !> CFC nuclear-position second derivative (f2_rArB) vs FD of f1_rA.
   subroutine test_cfc_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(rA_fwd)) deallocate (rA_fwd, rA_fwd2, rA_bwd, rA_bwd2)
         allocate (rA_fwd(ndim, mol%nat), rA_fwd2(ndim, mol%nat))
         allocate (rA_bwd(ndim, mol%nat), rA_bwd2(ndim, mol%nat))
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         if (allocated(analytic)) deallocate (analytic)
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol%nat))
         allocate (analytic(ndim, mol%nat, ndim, mol%nat))

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(3)
         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            call prim%set_centers(centers_base)
            call prim%prepare(point, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if
            call prim%f2_rArB(analytic)
            !* Active-indexed result: assert the shape and that the active list
            !* is the identity here, so the atom loops below may index
            !* `analytic` with user-space ids.
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
                        numeric = fd4_scalar(rA_fwd2(axisA, atomA), rA_fwd(axisA, atomA), &
                                             rA_bwd(axisA, atomA), rA_bwd2(axisA, atomA), STEP_SIZE)
                        call check(error, analytic(axisA, atomA, axisB, atomB), numeric, &
                                   thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_cfc_f2_rArB

   !> CFC mixed third (spatial grad x nuclear Hessian) vs FD of f2_r_rA.
   subroutine test_cfc_f3_r_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(r_rA_fwd)) deallocate (r_rA_fwd, r_rA_fwd2, r_rA_bwd, r_rA_bwd2)
         allocate (r_rA_fwd(ndim, ndim, mol%nat), r_rA_fwd2(ndim, ndim, mol%nat))
         allocate (r_rA_bwd(ndim, ndim, mol%nat), r_rA_bwd2(ndim, ndim, mol%nat))
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         if (allocated(analytic)) deallocate (analytic)
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol%nat))
         allocate (analytic(ndim, ndim, mol%nat, ndim, mol%nat))

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(3)
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
                           numeric = fd4_scalar(r_rA_fwd2(jdir, axisA, iA), &
                                                r_rA_fwd(jdir, axisA, iA), r_rA_bwd(jdir, axisA, iA), &
                                                r_rA_bwd2(jdir, axisA, iA), STEP_SIZE)
                           call check(error, analytic(jdir, axisA, iA, axisB, iB), numeric, &
                                      thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_cfc_f3_r_rArB

   !> CFC pure spatial fourth derivative vs FD of f3_rrr.
   subroutine test_cfc_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, axis, i, j, kk
      real(wp) :: point(ndim), work_point(ndim)
      real(wp) :: t3_fwd(ndim, ndim, ndim), t3_fwd2(ndim, ndim, ndim)
      real(wp) :: t3_bwd(ndim, ndim, ndim), t3_bwd2(ndim, ndim, ndim)
      real(wp) :: analytic(ndim, ndim, ndim, ndim)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)
         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(4)
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
                        numeric = fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), &
                                             t3_bwd(i, j, kk), t3_bwd2(i, j, kk), STEP_SIZE)
                        call check(error, analytic(i, j, kk, axis), numeric, &
                                   thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_cfc_f4_rrrr

   !> CFC mixed fourth (3 r-axes + 1 R-axis) vs FD of f3_rrr.
   subroutine test_cfc_f4_rrr_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, atom, axis, i, j, kk
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :)
      real(wp) :: t3_fwd(ndim, ndim, ndim), t3_fwd2(ndim, ndim, ndim)
      real(wp) :: t3_bwd(ndim, ndim, ndim), t3_bwd2(ndim, ndim, ndim)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(analytic)) deallocate (analytic)
         allocate (analytic(ndim, ndim, ndim, ndim, mol%nat))

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(4)
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
                           numeric = fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), &
                                                t3_bwd(i, j, kk), t3_bwd2(i, j, kk), STEP_SIZE)
                           call check(error, analytic(i, j, kk, axis, atom), numeric, &
                                      thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_cfc_f4_rrr_rA

   !> CFC mixed fourth (2 r-axes + 2 R-axes) vs FD of f3_rr_rA.
   subroutine test_cfc_f4_rr_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer  :: icase, ipt, iA, iB, axisA, axisB, j, kk
      real(wp) :: point(ndim)
      real(wp), allocatable :: analytic(:, :, :, :, :, :)
      real(wp), allocatable :: rr_rA_fwd(:, :, :, :), rr_rA_fwd2(:, :, :, :)
      real(wp), allocatable :: rr_rA_bwd(:, :, :, :), rr_rA_bwd2(:, :, :, :)
      real(wp) :: numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(analytic)) deallocate (analytic)
         if (allocated(rr_rA_fwd)) deallocate (rr_rA_fwd, rr_rA_fwd2, rr_rA_bwd, rr_rA_bwd2)
         allocate (analytic(ndim, ndim, ndim, mol%nat, ndim, mol%nat))
         allocate (rr_rA_fwd(ndim, ndim, ndim, mol%nat), rr_rA_fwd2(ndim, ndim, ndim, mol%nat))
         allocate (rr_rA_bwd(ndim, ndim, ndim, mol%nat), rr_rA_bwd2(ndim, ndim, ndim, mol%nat))

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(4)
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
                              numeric = fd4_scalar(rr_rA_fwd2(j, kk, axisA, iA), &
                                                   rr_rA_fwd(j, kk, axisA, iA), rr_rA_bwd(j, kk, axisA, iA), &
                                                   rr_rA_bwd2(j, kk, axisA, iA), STEP_SIZE)
                              call check(error, analytic(j, kk, axisA, iA, axisB, iB), &
                                         numeric, thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
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
   end subroutine test_cfc_f4_rr_rArB

   !> CFC normalized-LSF nuclear gradient vs FD of normalized f0.
   subroutine test_cfc_normalized_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      real(wp), allocatable :: deriv_rA(:, :)
      integer  :: icase, ipt, atom, axis
      real(wp) :: point(ndim)
      real(wp) :: normalized_val
      real(wp) :: f_forward, f_backward, f_forward2, f_backward2, numeric
      type(mctc_error), allocatable :: lsf_err

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)
         allocate (centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(deriv_rA)) deallocate (deriv_rA)
         allocate (deriv_rA(ndim, mol%nat))

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(2)
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
                  numeric = fd4_scalar(f_forward2, f_forward, f_backward, f_backward2, STEP_SIZE)
                  call check(error, deriv_rA(axis, atom), numeric, &
                             thr_abs=CFC_ABS_THR, thr_rel=CFC_REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_cfc_normalized_f1_rA

   !> The CFC `tangent_*` / `hvp_*` accessors must equal the explicit contraction
   !>
   !> The CFC twin of [[test_svdw_contracted]], and the highest-value single test
   !> of this block: it exercises all four `tangent_*` and all three `hvp_*`
   !> bindings against the uncontracted tensors they are supposed to summarise,
   !> which are themselves finite-difference tested above. The two sides are
   !> genuinely independent -- the contracted families come out of a different
   !> derivation ladder in the generated kernel (`tg`/`hv`) than the uncontracted
   !> ones (`qn`/`qq`), and a sign or index error in either shows up here at
   !> roundoff rather than at finite-difference resolution.
   subroutine test_cfc_contracted(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: prim
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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         nat = mol%nat
         call get_test_radii(mol, radii)
         call get_test_points(mol, points, n_points_cfc)

         !* A deterministic, non-symmetric direction field
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

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol, radii)
         call prim%set_max_deriv(4)

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

            !* --------------------------- tangent_* ---------------------------- *!
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

            !* ----------------------------- hvp_* ------------------------------ *!
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

            !* ---------------------------- vjp_f1_rA --------------------------- *!
            !* The reverse mirror of `tangent_*`: the jet indices are contracted
            !* away, the nuclear index survives. The weights are deliberately
            !* non-symmetric in `w2` -- a kernel that folded the Hessian slot
            !* into a symmetric half would pass a symmetric probe and fail here.
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

            !* --------------------------- vjp_f1_rad --------------------------- *!
            !* The same contraction against the radius ladder. A radius is a
            !* scalar parameter, so the surviving index is the atom alone and
            !* the result is one number per active slot.
            call prim%f3_rr_rad(rad1, rad2, rad3)
            call prim%vjp_f1_rad(w0_adj, w1_adj, w2_adj, vjp_rad)
            do iA = 1, prim%active_count()
               acc = w0_adj*rad1(iA) &
                     + dot_product(w1_adj, rad2(:, iA)) &
                     + sum(w2_adj*rad3(:, :, iA))
               call check(error, vjp_rad(iA), acc, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do

            !* ---------------- nuclear Hessian exchange symmetry ---------------- *!
            !* d2S/(dR_A dR_B) must equal d2S/(dR_B dR_A) under a simultaneous
            !* swap of both nuclear slots. The two sides come from *different*
            !* kernel calls (the `qq` cross block is built with A in the pair's
            !* `a` slot and B in its `b` slot, and vice versa), so this is a real
            !* check of the two-nucleus lift, not an identity.
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
   end subroutine test_cfc_contracted

   !* ================================================================================= *!
   !*                          A point sitting on a nucleus                             *!
   !* ================================================================================= *!

   !> No accessor may produce NaN or Inf at a point that lands on a nucleus
   !>
   !> The signed sphere distance `||r - R_A||` is not differentiable there, and
   !> the generated kernel divides by it from the first derivative onwards. The
   !> LSF intercepts that case (see `atom_tensors`): the atom keeps its value
   !> contribution and loses every derivative one, which is the convention the
   !> retired SSD fill used. This pins that the guard is actually in the path of
   !> *every* accessor -- a single unguarded call site would show up as a NaN.
   !>
   !> The value is additionally checked for continuity against a point a
   !> whisker off the nucleus: `f0` itself is continuous there, only its
   !> derivatives are not.
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
      call prim%new(blend_k=svdw_legacy_blend_k, blend_2b=svdw_legacy_blend_2b, &
                    blend_3b=svdw_legacy_blend_3b)
      call prim%update(mol, radii)
      call prim%set_max_deriv(4)

      !* Exactly on nucleus 2
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

      !* The value is continuous across the nucleus even though its slope is not
      point = centers(:, 2) + [nudge, 0.0_wp, 0.0_wp]
      call prim%prepare(point, lsf_err)
      call prim%f0(f0_near)
      call check(error, f0_on, f0_near, thr_abs=1.0e-6_wp, thr_rel=1.0e-6_wp, &
                 more="f0 jumps at a nucleus")
   end subroutine test_svdw_on_nucleus

   !> Is this scalar a finite number?
   !> @param[in] x Value to test
   !> @returns     `.true.` unless `x` is NaN or infinite
   pure logical function finite0(x) result(ok)
      !> Value to test
      real(wp), intent(in) :: x
      ok = (x == x) .and. abs(x) <= huge(1.0_wp)
   end function finite0

   !> Is every element of a rank-1 array finite?
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

   !> Is every element of a vector finite?
   !> @param[in] x Vector to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finite1(x) result(ok)
      !> Vector to test
      real(wp), intent(in) :: x(:)
      ok = finiten(x)
   end function finite1

   !> Is every element of a matrix finite?
   !> @param[in] x Matrix to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finite2(x) result(ok)
      !> Matrix to test
      real(wp), intent(in) :: x(:, :)
      ok = finiten(reshape(x, [size(x)]))
   end function finite2

   !> Is every element of a rank-3 tensor finite?
   !> @param[in] x Tensor to test
   !> @returns     `.true.` unless any element is NaN or infinite
   pure logical function finite3(x) result(ok)
      !> Tensor to test
      real(wp), intent(in) :: x(:, :, :)
      ok = finiten(reshape(x, [size(x)]))
   end function finite3

   !> Signed distance from a point to a sphere surface
   !>
   !> The reference expression the body-order reduction is checked against.
   !> Spelled out here rather than imported: the level set function no longer
   !> owns a signed-distance module, the kernel is symbolic in the screening
   !> factors instead.
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
      call prim%new(blend_k=2.4_wp, blend_1b=0.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp)
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
      call prim%new(blend_k=1.7_wp, blend_1b=0.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp)
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
   !*                        Surface-free exclusion certificate                          *!
   !* ================================================================================= *!

   !> The exclusion radius must never over-claim.
   !>
   !> `exclusion_radius(S(x))` promises that `S` has no zero inside `B(x, r)`.
   !> The certified branch search rests entirely on that promise -- a radius one
   !> per cent too large silently hides branches -- so this samples the ball
   !> around a spread of points and insists `S` keeps its sign throughout.
   !>
   !> The SvdW claim is `r = |S(x)|` from the 1-Lipschitz property, which is
   !> tight: points at distance exactly `r` may sit on the surface. Sampling
   !> stops just inside at `0.999 r`.
   subroutine test_svdw_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      !> Fraction of the claimed radius actually probed
      real(wp), parameter :: probe_frac = 0.999_wp

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      type(mctc_error), allocatable :: lsf_err
      real(wp) :: centre(ndim), probe(ndim), dir(ndim)
      real(wp) :: lsf_centre, lsf_probe, r, dir_norm
      integer :: ipoint, idir, iblend

      ! A deterministic spread of probe directions; no RNG, so a failure here
      ! reproduces exactly.
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

      call get_structure(mol, "MB16-43", "01")
      call get_test_radii(mol, radii)
      call get_test_points(mol, points)

      do iblend = 1, n_svdw_blends
         call init_lsf(lsf, mol, radii, 0, kind_svdw, &
                       blend_k=svdw_blend_k_values(iblend))

         do ipoint = 1, size(points, 2)
            centre = points(:, ipoint)

            call lsf%prepare(centre, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if
            call lsf%f0(lsf_centre)

            r = lsf%exclusion_radius(lsf_centre)
            if (r <= 0.0_wp) cycle

            do idir = 1, ndir
               dir = dirs(:, idir)
               dir_norm = norm2(dir)
               probe = centre + probe_frac*r*dir/dir_norm

               call lsf%prepare(probe, lsf_err)
               if (allocated(lsf_err)) then
                  call test_failed(error, "LSF prepare failed: "//lsf_err%message)
                  return
               end if
               call lsf%f0(lsf_probe)

               if (lsf_probe*lsf_centre <= 0.0_wp) then
                  call test_failed(error, &
                                   "exclusion radius over-claims: S changes sign inside "// &
                                   "the ball it certifies as surface-free")
                  return
               end if
            end do
         end do
      end do
   end subroutine test_svdw_exclusion_radius

   !> The 1-Lipschitz bound needs non-negative blend coefficients. With a
   !> negative one the weighted mean of unit vectors becomes an extrapolation,
   !> the bound is lost, and the LSF must decline to certify anything rather
   !> than hand back a radius it cannot stand behind.
   subroutine test_svdw_exclusion_radius_gate(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)

      call init_lsf(lsf, mol, radii, 0, kind_svdw, blend_2b=1.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 1.0_wp, thr=0.0_wp, &
                 message="non-negative blends should certify r = |S|")
      if (allocated(error)) return

      call init_lsf(lsf, mol, radii, 0, kind_svdw, blend_2b=-1.0_wp)
      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="a negative blend coefficient must certify nothing")
   end subroutine test_svdw_exclusion_radius_gate

   !> CFC has no derived gradient bound yet, so it inherits the base "no
   !> certificate" answer. Level 9 then runs out of budget instead of
   !> reporting a completeness it cannot back up.
   subroutine test_cfc_exclusion_radius(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 0, kind_cfc)

      call check(error, lsf%exclusion_radius(-1.0_wp), 0.0_wp, thr=0.0_wp, &
                 message="CFC must not claim an exclusion radius")
   end subroutine test_cfc_exclusion_radius

end module test_cavity_drop_lsf
