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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
                  !* Active-indexed result: assert the shape and that the active
                  !* list is the identity here, so the atom loops below may index
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

   !> The `tangent_*` / `hvp_*` accessors must equal the explicit contraction
   !>
   !> Both families exist only so a caller never has to form the O(n**2) nuclear
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
                                           f2_rArB, f3_r_rArB, f4_rr_rArB, h1, h2, h3)
         allocate (f1_rA(ndim, nat), f2_r_rA(ndim, ndim, nat))
         allocate (f3_rr_rA(ndim, ndim, ndim, nat))
         allocate (f4_rrr_rA(ndim, ndim, ndim, ndim, nat))
         allocate (f2_rArB(ndim, nat, ndim, nat))
         allocate (f3_r_rArB(ndim, ndim, nat, ndim, nat))
         allocate (f4_rr_rArB(ndim, ndim, ndim, nat, ndim, nat))
         allocate (h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat))

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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
         call get_test_points(mol, points)
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
                                           f2_rArB, f3_r_rArB, f4_rr_rArB, h1, h2, h3)
         allocate (f1_rA(ndim, nat), f2_r_rA(ndim, ndim, nat))
         allocate (f3_rr_rA(ndim, ndim, ndim, nat))
         allocate (f4_rrr_rA(ndim, ndim, ndim, ndim, nat))
         allocate (f2_rArB(ndim, nat, ndim, nat))
         allocate (f3_r_rArB(ndim, ndim, nat, ndim, nat))
         allocate (f4_rr_rArB(ndim, ndim, ndim, nat, ndim, nat))
         allocate (h1(ndim, nat), h2(ndim, ndim, nat), h3(ndim, ndim, ndim, nat))

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
