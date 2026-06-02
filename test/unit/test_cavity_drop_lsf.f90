!> Orchestrator-level unit tests for the DROP level-set functions (LSF)
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
!> SvdW-only extensions (registered with svdw_ prefix only):
!>   * `f2_rArB`, `f3_r_rArB`                  pure/mixed nuclear seconds
!>   * `f4_rrrr`, `f4_rrr_rA`, `f4_rr_rArB`    fourth-order derivatives
!>   * `pou_f1_r`, `pou_f2_rr`, `pou_f2_r_rA`  partition-of-unity tests
!>   * `normalized_f1_rA`                      normalized LSF nuclear grad
!>   * `body_order_scaling`                    1b/2b/3b weight reduction
!>
!> Test fixtures come from the shared MB16-43/Heavy28/Amino20x4/But14diol/UPU23
!> palette in `test_helpers::get_test_structures`. Per-atom radii come from
!> the project's standard CPCM table via `get_test_radii`, and sampling
!> points come from `get_test_points`.
!>
!> Nuclear-derivative FDs bypass `lsf%update(mol, radii)` (which would
!> re-init the SSD and may renumber atoms) by calling
!> `lsf%ssd_system%update(...)` directly via a `select type` shim. The
!> shim is needed because the SSD system is held by each concrete subtype
!> (CFC and SvdW use the same `moist_cavity_drop_lsf_svdw_ssd_type`) and
!> is not part of the abstract base API.
module test_cavity_drop_lsf
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, fd4_scalar
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_ssd, only: ssd0
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_cavity_drop_lsf

   integer, parameter :: ndim = 3

   !> Concrete LSF selector strings; used by [[init_lsf]]
   character(len=*), parameter :: kind_svdw = "svdw"
   character(len=*), parameter :: kind_cfc  = "cfc"

   real(wp), parameter :: fd_h = 5.0e-4_wp
   real(wp), parameter :: atol = 2.0e-10_wp
   real(wp), parameter :: rtol = 1.0e-9_wp

   integer,  parameter :: n_svdw_blends = 5
   integer,  parameter :: n_svdw_gammas = 2
   real(wp), parameter :: svdw_blend_k_values(n_svdw_blends) = [0.1_wp, 0.5_wp, 1.0_wp, 2.0_wp, 3.0_wp]
   real(wp), parameter :: svdw_gamma_values(n_svdw_gammas) = [0.0_wp, 0.7_wp]

   ! Screening-cutoff sweep
   real(wp), parameter :: screening_thresholds(3) = [ &
      1.0e-12_wp, 1.0e-13_wp, 1.0e-14_wp]
   ! Number of 1-bohr outward steps marched from the molecule's bounding-box edge
   integer,  parameter :: screening_n_steps = 15
   real(wp), parameter :: screening_step_bohr = 1.0_wp
   ! Stop the march once -log(PD0) exceeds this
   real(wp), parameter :: screening_val_ref_max = 5.0_wp
   ! Floating-point noise floor for the screening diff comparison
   real(wp), parameter :: screening_noise_floor = 1.0e-12_wp

contains

   subroutine collect_cavity_drop_lsf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         !> Common
         new_unittest("svdw_sign_convention", test_svdw_sign_convention), &
         new_unittest("cfc_sign_convention", test_cfc_sign_convention), &
         new_unittest("svdw_neighbor_cutoff", test_svdw_neighbor_cutoff), &
         new_unittest("cfc_neighbor_cutoff", test_cfc_neighbor_cutoff), &
         !> LSF
         new_unittest("svdw_f1_r_fd", test_svdw_f1_r_fd), &
         new_unittest("cfc_f1_r_fd", test_cfc_f1_r_fd), &
         new_unittest("svdw_f2_rr_fd", test_svdw_f2_rr_fd), &
         new_unittest("cfc_f2_rr_fd", test_cfc_f2_rr_fd), &
         new_unittest("svdw_f3_rrr_fd", test_svdw_f3_rrr_fd), &
         new_unittest("cfc_f3_rrr_fd", test_cfc_f3_rrr_fd), &
         new_unittest("svdw_f1_rA_fd", test_svdw_f1_rA_fd), &
         new_unittest("cfc_f1_rA_fd", test_cfc_f1_rA_fd), &
         new_unittest("svdw_f2_r_rA_fd", test_svdw_f2_r_rA_fd), &
         new_unittest("cfc_f2_r_rA_fd", test_cfc_f2_r_rA_fd), &
         new_unittest("svdw_f3_rr_rA_fd", test_svdw_f3_rr_rA_fd), &
         new_unittest("cfc_f3_rr_rA_fd", test_cfc_f3_rr_rA_fd), &
         !> SvdW-only
         new_unittest("svdw_f2_rArB", test_svdw_f2_rArB), &
         new_unittest("svdw_f3_r_rArB", test_svdw_f3_r_rArB), &
         new_unittest("svdw_f4_rrrr", test_svdw_f4_rrrr), &
         new_unittest("svdw_f4_rrr_rA", test_svdw_f4_rrr_rA), &
         new_unittest("svdw_f4_rr_rArB", test_svdw_f4_rr_rArB), &
         new_unittest("svdw_pou_f1_r", test_svdw_pou_f1_r), &
         new_unittest("svdw_pou_f2_rr", test_svdw_pou_f2_rr), &
         new_unittest("svdw_pou_f2_r_rA", test_svdw_pou_f2_r_rA), &
         new_unittest("svdw_normalized_f1_rA", test_svdw_normalized_f1_rA), &
         new_unittest("svdw_body_order_scaling", test_svdw_body_order_scaling) &
         ]
   end subroutine collect_cavity_drop_lsf

   !* ================================================================================= *!
   !*                              Local helpers                                        *!
   !* ================================================================================= *!

   !> Allocate a fresh polymorphic LSF of the requested concrete kind and
   !> bind it to the given molecule. Optional `blend_k` / `blend_3b`
   !> apply only when `kind == "svdw"` (CFC has no equivalent knob).
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
   !> @param[in]  blend_3b            optional svdw 3-body weight override
   !> @param[in]  screening_threshold optional screening cutoff (default 0)
   subroutine init_lsf(lsf, mol, radii, max_deriv, kind, blend_k, blend_3b, &
         screening_threshold)
      class(moist_cavity_drop_lsf_type), allocatable, intent(out) :: lsf
      !> Molecular structure to bind to the LSF
      type(structure_type), intent(in) :: mol
      !> Per-atom radii (size mol%nat)
      real(wp), intent(in) :: radii(:)
      !> Highest spatial derivative order to enable
      integer,  intent(in) :: max_deriv
      !> Concrete kind selector: "svdw" or "cfc"
      character(len=*), intent(in) :: kind
      !> Optional SvdW blend_k override (ignored for CFC)
      real(wp), intent(in), optional :: blend_k
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
         call tmp_svdw%new(blend_k=blend_k, blend_3b=blend_3b)
         call tmp_svdw%update(mol, radii)
         call tmp_svdw%set_max_deriv(max_deriv)
         allocate(lsf, source=tmp_svdw)
      case (kind_cfc)
         tmp_cfc%screening_threshold = thr
         call tmp_cfc%new()
         call tmp_cfc%update(mol, radii)
         call tmp_cfc%set_max_deriv(max_deriv)
         allocate(lsf, source=tmp_cfc)
      case default
         error stop "init_lsf: unknown kind '"//kind//"'"
      end select
   end subroutine init_lsf

   !> Refresh only the SSD geometry on the underlying concrete (no full
   !> LSF re-init). Used by nuclear-derivative FDs to perturb atom
   !> positions without wiping `max_deriv` or other concrete caches.
   !> Performs a `select type` dispatch because `ssd_system` lives on
   !> the concretes, not the abstract base.
   !>
   !> @param[inout] lsf      polymorphic LSF (must be allocated)
   !> @param[in]    centers  perturbed positions (3, mol%nat)
   !> @param[in]    radii    per-atom radii (size mol%nat)
   subroutine refresh_ssd(lsf, centers, radii)
      !> Polymorphic LSF whose SSD geometry is refreshed in place
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Perturbed atom positions (3, mol%nat)
      real(wp), intent(in) :: centers(:, :)
      !> Per-atom radii (size mol%nat)
      real(wp), intent(in) :: radii(:)

      select type (lsf)
      type is (moist_cavity_drop_lsf_svdw_type)
         call lsf%ssd_system%update(centers, radii)
      type is (moist_cavity_drop_lsf_cfc_type)
         call lsf%ssd_system%update(centers, radii)
      class default
         error stop "refresh_ssd: unknown concrete LSF kind"
      end select
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

      call get_structure(mol, "MB16-43", "LiH")
      call get_test_radii(mol, radii)
      call init_lsf(lsf, mol, radii, 0, kind)

      !* Inside: lies on atom 1, well inside the cavity.
      inside_pt = mol%xyz(:, 1) + [0.05_wp, 0.05_wp, 0.0_wp]
      call lsf%prepare(inside_pt)
      call lsf%f0_screened(val_inside)
      if (val_inside >= 0.0_wp) then
         call test_failed(error, "Expected lsf0 < 0 deep inside the cavity")
         return
      end if

      !* Outside: well beyond the molecular bounding box along +x.
      outside_pt = [maxval(mol%xyz(1, :)) + 20.0_wp, 0.0_wp, 0.0_wp]
      call lsf%prepare(outside_pt)
      call lsf%f0_screened(val_outside)
      if (val_outside <= 0.0_wp) then
         call test_failed(error, "Expected lsf0 > 0 well outside the cavity")
         return
      end if
   end subroutine run_sign_convention

   !> SvdW dispatch for the neighbour-cutoff check.
   subroutine test_svdw_neighbor_cutoff(error)
      type(error_type), allocatable, intent(out) :: error
      call run_neighbor_cutoff(error, kind_svdw)
   end subroutine test_svdw_neighbor_cutoff

   !> CFC dispatch for the neighbour-cutoff check.
   subroutine test_cfc_neighbor_cutoff(error)
      type(error_type), allocatable, intent(out) :: error
      call run_neighbor_cutoff(error, kind_cfc)
   end subroutine test_cfc_neighbor_cutoff

   !> Verify the screening contract: a screened LSF with threshold X
   !> differs from the unscreened (threshold=0) LSF by at most a bound
   !> proportional to X, in the well-conditioned LSF regime.
   !>
   !> The check sweeps three tight thresholds (1e-12, 1e-13, 1e-14) and
   !> marches outward from the molecule's bounding-box edge in 1-bohr
   !> steps, stopping per-threshold once either (a) every atom has
   !> been screened (n_active == 0) or (b) the unscreened LSF value
   !> exceeds `screening_val_ref_max`. Looser thresholds and points
   !> with larger val_ref were excluded because they push the kernel
   !> into the regime where the underlying Z is tiny and the kernel's
   !> -log(Z)/blend_k amplification turns sub-threshold absolute drops
   !> into LSF-value changes well above the threshold; that's a
   !> kernel-conditioning issue rather than a screening bug, and the
   !> cavity machinery never evaluates LSF in that regime.
   !>
   !> Bound rationale: an atom passes the SSD screening only if its
   !> weight at the point is >= X. For SvdW this directly bounds its
   !> kernel contribution by X. For CFC the atomic exponent (a1=-15)
   !> is stricter than the SSD exponent (screen_k=3), so atomic terms
   !> from dropped atoms decay much faster than X - but pair terms
   !> between a dropped atom and a kept atom can persist. Counting both
   !> atomic and pair channels gives a worst-case
   !> P_dropped <= n_dropped*nat*X, so the LSF diff is bounded by
   !> ~nat^2*X over the noise floor.
   subroutine run_neighbor_cutoff(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf_ref, lsf_scr
      real(wp), allocatable :: radii(:)
      real(wp) :: centroid(ndim), direction(ndim), point(ndim)
      real(wp) :: max_extent, thr, val_ref, val_scr, diff, bound
      integer :: icase, ithr, istep, n_active_scr

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)

         !* Reference LSF: no screening.
         call init_lsf(lsf_ref, mol, radii, 0, kind)

         !* March outward from the molecule along +x. The starting offset
         !* is the bounding-box max along +x so the first point sits near
         !* the molecular surface where screening just begins to engage.
         centroid = sum(mol%xyz, dim=2) / real(mol%nat, wp)
         max_extent = maxval(mol%xyz(1, :)) - centroid(1)
         direction = [1.0_wp, 0.0_wp, 0.0_wp]

         do ithr = 1, size(screening_thresholds)
            thr = screening_thresholds(ithr)
            call init_lsf(lsf_scr, mol, radii, 0, kind, screening_threshold=thr)

            do istep = 0, screening_n_steps
               point = centroid + (max_extent &
                  + real(istep, wp) * screening_step_bohr) * direction

               call lsf_ref%prepare(point)
               call lsf_ref%f0_screened(val_ref)
               call lsf_scr%prepare(point)
               call lsf_scr%f0_screened(val_scr)

               !* Two march-termination guards:
               !*   (1) all atoms screened: the kernel reports n_active=0
               !*       and returns val=0; the comparison is undefined.
               !*   (2) val_ref exceeds the well-conditioned LSF range:
               !*       past this point the underlying Z is small enough
               !*       that the kernel amplifies sub-threshold absolute
               !*       drops into LSF-value changes above the threshold.
               n_active_scr = lsf_scr%active_count()
               if (n_active_scr == 0) exit
               if (val_ref > screening_val_ref_max) exit

               bound = real(mol%nat, wp)**2 * thr + screening_noise_floor
               diff = abs(val_scr - val_ref)
               if (diff > bound) then
                  call test_failed(error, &
                     "screened lsf diverged from unscreened beyond threshold bound")
                  return
               end if
            end do

            deallocate(lsf_scr)
         end do

         deallocate(lsf_ref)
      end do
   end subroutine run_neighbor_cutoff

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
                  call lsf%prepare(point)
                  call lsf%f012_r_screened(lsf1_r=analytic)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*fd_h
                     call lsf%prepare(shifted); call lsf%f0_screened(f_pp)
                     shifted = point; shifted(axis) = point(axis) + fd_h
                     call lsf%prepare(shifted); call lsf%f0_screened(f_p)
                     shifted = point; shifted(axis) = point(axis) - fd_h
                     call lsf%prepare(shifted); call lsf%f0_screened(f_m)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*fd_h
                     call lsf%prepare(shifted); call lsf%f0_screened(f_mm)
                     numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, fd_h)
                  end do
                  do i = 1, ndim
                     call check(error, analytic(i), numeric(i), &
                        thr_abs=atol, thr_rel=rtol)
                     if (allocated(error)) return
                  end do
               end do
               deallocate(lsf)
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
                  call lsf%prepare(point)
                  call lsf%f012_r_screened(lsf2_rr=analytic)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*fd_h
                     call lsf%prepare(shifted); call lsf%f012_r_screened(lsf1_r=g_pp)
                     shifted = point; shifted(axis) = point(axis) + fd_h
                     call lsf%prepare(shifted); call lsf%f012_r_screened(lsf1_r=g_p)
                     shifted = point; shifted(axis) = point(axis) - fd_h
                     call lsf%prepare(shifted); call lsf%f012_r_screened(lsf1_r=g_m)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*fd_h
                     call lsf%prepare(shifted); call lsf%f012_r_screened(lsf1_r=g_mm)
                     do i = 1, ndim
                        numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), fd_h)
                     end do
                  end do
                  do j = 1, ndim
                     do i = 1, ndim
                        call check(error, analytic(i, j), numeric(i, j), &
                           thr_abs=atol, thr_rel=rtol)
                        if (allocated(error)) return
                     end do
                  end do
               end do
               deallocate(lsf)
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
   !> f3_rrr_screened so the FD and analytic branches share the same
   !> internal code path; matters at -O3).
   subroutine run_f3_rrr_fd(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, axis, i, j, k, iblend, igamma, nblend, ngamma
      real(wp), allocatable :: analytic(:, :, :), dummy_third(:, :, :)
      real(wp) :: numeric(ndim, ndim, ndim), point(ndim), shifted(ndim)
      real(wp) :: hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
      real(wp) :: eps

      eps = fd_h
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
                  call lsf%prepare(point)
                  call lsf%f3_rrr_screened(lsf3_rrr=analytic)
                  do axis = 1, ndim
                     shifted = point; shifted(axis) = point(axis) + 2.0_wp*eps
                     call lsf%prepare(shifted)
                     call lsf%f3_rrr_screened(lsf2_rr=hess_pp, lsf3_rrr=dummy_third)
                     shifted = point; shifted(axis) = point(axis) + eps
                     call lsf%prepare(shifted)
                     call lsf%f3_rrr_screened(lsf2_rr=hess_p, lsf3_rrr=dummy_third)
                     shifted = point; shifted(axis) = point(axis) - eps
                     call lsf%prepare(shifted)
                     call lsf%f3_rrr_screened(lsf2_rr=hess_m, lsf3_rrr=dummy_third)
                     shifted = point; shifted(axis) = point(axis) - 2.0_wp*eps
                     call lsf%prepare(shifted)
                     call lsf%f3_rrr_screened(lsf2_rr=hess_mm, lsf3_rrr=dummy_third)
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
                              thr_abs=atol, thr_rel=rtol)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
               deallocate(lsf)
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

      eps = fd_h
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                  blend_k=svdw_sweep_blend(kind, iblend), &
                  blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate(analytic(ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_ssd(lsf, centers_base, radii)
                  call lsf%prepare(point)
                  call lsf%f3_rr_rA_screened(lsf1_rA=analytic, lsf3_rr_rA=dummy_3rd)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f0_screened(f_pp)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f0_screened(f_p)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f0_screened(f_m)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f0_screened(f_mm)
                        numeric = fd4_scalar(f_pp, f_p, f_m, f_mm, eps)
                        call check(error, analytic(axis, atom), numeric, &
                           thr_abs=atol, thr_rel=rtol)
                        if (allocated(error)) return
                     end do
                  end do
               end do
               deallocate(analytic, lsf)
            end do
         end do
         deallocate(centers_base, centers_local)
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

      eps = fd_h
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                  blend_k=svdw_sweep_blend(kind, iblend), &
                  blend_3b=svdw_sweep_gamma(kind, igamma))
               allocate(analytic(ndim, ndim, mol%nat))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_ssd(lsf, centers_base, radii)
                  call lsf%prepare(point)
                  call lsf%f3_rr_rA_screened(lsf2_r_rA=analytic, lsf3_rr_rA=dummy_3rd)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf1_r=g_pp)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf1_r=g_p)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf1_r=g_m)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf1_r=g_mm)
                        do i = 1, ndim
                           numeric = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), eps)
                           call check(error, analytic(i, axis, atom), numeric, &
                              thr_abs=atol, thr_rel=rtol)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
               deallocate(analytic, lsf)
            end do
         end do
         deallocate(centers_base, centers_local)
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

      eps = fd_h
      call svdw_sweep_sizes(kind, nblend, ngamma)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         do iblend = 1, nblend
            do igamma = 1, ngamma
               call init_lsf(lsf, mol, radii, 3, kind, &
                  blend_k=svdw_sweep_blend(kind, iblend), &
                  blend_3b=svdw_sweep_gamma(kind, igamma))
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call refresh_ssd(lsf, centers_base, radii)
                  call lsf%prepare(point)
                  call lsf%f3_rr_rA_screened(lsf3_rr_rA=analytic)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf2_rr=hess_pp)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf2_rr=hess_p)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf2_rr=hess_m)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                        call refresh_ssd(lsf, centers_local, radii)
                        call lsf%prepare(point); call lsf%f012_r_screened(lsf2_rr=hess_mm)
                        do j = 1, ndim
                           do i = 1, ndim
                              numeric = fd4_scalar(hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), eps)
                              call check(error, analytic(i, j, axis, atom), numeric, &
                                 thr_abs=atol, thr_rel=rtol)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
               deallocate(lsf)
            end do
         end do
         deallocate(centers_base, centers_local)
      end do
   end subroutine run_f3_rr_rA_fd

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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         !* FD buffers are sized by mol%nat only and reused across the
         !* (iblend, igamma, ipt) sweep.
         if (allocated(rA_fwd))  deallocate(rA_fwd)
         if (allocated(rA_fwd2)) deallocate(rA_fwd2)
         if (allocated(rA_bwd))  deallocate(rA_bwd)
         if (allocated(rA_bwd2)) deallocate(rA_bwd2)
         allocate(rA_fwd(ndim, mol%nat), rA_fwd2(ndim, mol%nat))
         allocate(rA_bwd(ndim, mol%nat), rA_bwd2(ndim, mol%nat))
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(3)
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%ssd_system%update(centers_base, radii)
                  call prim%prepare(point)
                  call prim%f2_rArB_screened(analytic)
                  do atomB = 1, mol%nat
                     do axisB = 1, ndim
                        centers_local = centers_base
                        centers_local(axisB, atomB) = centers_local(axisB, atomB) + fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf1_rA=rA_fwd, lsf3_rr_rA=dummy_rr_rA)
                        centers_local = centers_base
                        centers_local(axisB, atomB) = centers_local(axisB, atomB) + 2.0_wp*fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf1_rA=rA_fwd2, lsf3_rr_rA=dummy_rr_rA)
                        centers_local = centers_base
                        centers_local(axisB, atomB) = centers_local(axisB, atomB) - fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf1_rA=rA_bwd, lsf3_rr_rA=dummy_rr_rA)
                        centers_local = centers_base
                        centers_local(axisB, atomB) = centers_local(axisB, atomB) - 2.0_wp*fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf1_rA=rA_bwd2, lsf3_rr_rA=dummy_rr_rA)
                        do atomA = 1, mol%nat
                           do axisA = 1, ndim
                              numeric = fd4_scalar(rA_fwd2(axisA, atomA), rA_fwd(axisA, atomA), &
                                 rA_bwd(axisA, atomA), rA_bwd2(axisA, atomA), fd_h)
                              call check(error, analytic(axisA, atomA, axisB, atomB), numeric, &
                                 thr_abs=atol, thr_rel=rtol)
                              if (allocated(error)) return
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate(centers_base, centers_local)
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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         !* f3_r_rArB_screened declares lsf1_rA and lsf2_r_rA as intent(in)
         !* non-allocatable assumed-shape; passing unallocated allocatables
         !* is undefined behavior, so size the dummies up front. All
         !* buffers depend only on mol%nat, so allocate once per icase.
         if (allocated(dummy_rA))  deallocate(dummy_rA)
         if (allocated(r_rA_fwd))  deallocate(r_rA_fwd)
         if (allocated(r_rA_fwd2)) deallocate(r_rA_fwd2)
         if (allocated(r_rA_bwd))  deallocate(r_rA_bwd)
         if (allocated(r_rA_bwd2)) deallocate(r_rA_bwd2)
         allocate(dummy_rA(ndim, mol%nat))
         allocate(r_rA_fwd(ndim, ndim, mol%nat), r_rA_fwd2(ndim, ndim, mol%nat))
         allocate(r_rA_bwd(ndim, ndim, mol%nat), r_rA_bwd2(ndim, ndim, mol%nat))
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
                  call prim%ssd_system%update(centers_base, radii)
                  call prim%prepare(point)
                  call prim%f3_r_rArB_screened(dummy_rA, r_rA_fwd, analytic)
                  do iB = 1, mol%nat
                     do axisB = 1, ndim
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) + fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf2_r_rA=r_rA_fwd, lsf3_rr_rA=dummy_rr_rA)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) + 2.0_wp*fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf2_r_rA=r_rA_fwd2, lsf3_rr_rA=dummy_rr_rA)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) - fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf2_r_rA=r_rA_bwd, lsf3_rr_rA=dummy_rr_rA)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) - 2.0_wp*fd_h
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf2_r_rA=r_rA_bwd2, lsf3_rr_rA=dummy_rr_rA)
                        do iA = 1, mol%nat
                           do axisA = 1, ndim
                              do jdir = 1, ndim
                                 numeric = fd4_scalar(r_rA_fwd2(jdir, axisA, iA), &
                                    r_rA_fwd(jdir, axisA, iA), r_rA_bwd(jdir, axisA, iA), &
                                    r_rA_bwd2(jdir, axisA, iA), fd_h)
                                 call check(error, analytic(jdir, axisA, iA, axisB, iB), numeric, &
                                    thr_abs=atol, thr_rel=rtol)
                                 if (allocated(error)) return
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate(centers_base, centers_local)
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
      real(wp), allocatable :: t3_fwd(:, :, :), t3_fwd2(:, :, :)
      real(wp), allocatable :: t3_bwd(:, :, :), t3_bwd2(:, :, :)
      real(wp), allocatable :: analytic(:, :, :, :)
      real(wp) :: numeric

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
                  call prim%prepare(point)
                  call prim%f4_rrrr_screened(analytic)
                  do axis = 1, ndim
                     work_point = point
                     work_point(axis) = point(axis) + fd_h
                     call prim%prepare(work_point)
                     call prim%f3_rrr_screened(lsf3_rrr=t3_fwd)
                     work_point = point
                     work_point(axis) = point(axis) + 2.0_wp*fd_h
                     call prim%prepare(work_point)
                     call prim%f3_rrr_screened(lsf3_rrr=t3_fwd2)
                     work_point = point
                     work_point(axis) = point(axis) - fd_h
                     call prim%prepare(work_point)
                     call prim%f3_rrr_screened(lsf3_rrr=t3_bwd)
                     work_point = point
                     work_point(axis) = point(axis) - 2.0_wp*fd_h
                     call prim%prepare(work_point)
                     call prim%f3_rrr_screened(lsf3_rrr=t3_bwd2)
                     do i = 1, ndim
                        do j = 1, ndim
                           do kk = 1, ndim
                              numeric = fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), &
                                 t3_bwd(i, j, kk), t3_bwd2(i, j, kk), fd_h)
                              call check(error, analytic(i, j, kk, axis), numeric, &
                                 thr_abs=atol, thr_rel=rtol)
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
      real(wp), allocatable :: t3_fwd(:, :, :), t3_fwd2(:, :, :)
      real(wp), allocatable :: t3_bwd(:, :, :), t3_bwd2(:, :, :)
      real(wp) :: numeric
      type(structure_type) :: mol_shift
      integer, allocatable :: atomic_numbers(:)
      integer :: iat

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate(atomic_numbers)
         allocate(atomic_numbers(mol%nat))
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
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%update(mol, radii)
                  call prim%ssd_system%update(centers_base, radii)
                  call prim%prepare(point)
                  call prim%f4_rrr_rA_screened(analytic)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rrr_screened(lsf3_rrr=t3_fwd)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rrr_screened(lsf3_rrr=t3_fwd2)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rrr_screened(lsf3_rrr=t3_bwd)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rrr_screened(lsf3_rrr=t3_bwd2)
                        do i = 1, ndim
                           do j = 1, ndim
                              do kk = 1, ndim
                                 numeric = fd4_scalar(t3_fwd2(i, j, kk), t3_fwd(i, j, kk), &
                                    t3_bwd(i, j, kk), t3_bwd2(i, j, kk), fd_h)
                                 call check(error, analytic(i, j, kk, axis, atom), numeric, &
                                    thr_abs=atol, thr_rel=rtol)
                                 if (allocated(error)) return
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate(centers_base, centers_local)
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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate(atomic_numbers)
         allocate(atomic_numbers(mol%nat))
         do iat = 1, mol%nat
            atomic_numbers(iat) = mol%num(mol%id(iat))
         end do
         !* f4_rr_rArB_screened declares lsf1_rA and lsf2_r_rA as intent(in)
         !* non-allocatable assumed-shape; passing unallocated allocatables
         !* is undefined behavior, so size the dummies to the active-atom
         !* count up front. Both depend only on mol%nat.
         if (allocated(dummy_rA))   deallocate(dummy_rA)
         if (allocated(dummy_r_rA)) deallocate(dummy_r_rA)
         allocate(dummy_rA(ndim, mol%nat), dummy_r_rA(ndim, ndim, mol%nat))
         dummy_rA = 0.0_wp
         dummy_r_rA = 0.0_wp
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
                  call prim%ssd_system%update(centers_base, radii)
                  call prim%prepare(point)
                  call prim%f4_rr_rArB_screened(dummy_rA, dummy_r_rA, analytic)
                  do iB = 1, mol%nat
                     do axisB = 1, ndim
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) + fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf3_rr_rA=rr_rA_fwd)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) + 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf3_rr_rA=rr_rA_fwd2)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) - fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf3_rr_rA=rr_rA_bwd)
                        centers_local = centers_base
                        centers_local(axisB, iB) = centers_local(axisB, iB) - 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%f3_rr_rA_screened(lsf3_rr_rA=rr_rA_bwd2)
                        do iA = 1, mol%nat
                           do axisA = 1, ndim
                              do j = 1, ndim
                                 do kk = 1, ndim
                                    numeric = fd4_scalar(rr_rA_fwd2(j, kk, axisA, iA), &
                                       rr_rA_fwd(j, kk, axisA, iA), rr_rA_bwd(j, kk, axisA, iA), &
                                       rr_rA_bwd2(j, kk, axisA, iA), fd_h)
                                    call check(error, analytic(j, kk, axisA, iA, axisB, iB), &
                                       numeric, thr_abs=atol, thr_rel=rtol)
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
         deallocate(centers_base, centers_local)
      end do
   end subroutine test_svdw_f4_rr_rArB

   !> SvdW partition-of-unity spatial gradient vs FD of weight.
   subroutine test_svdw_pou_f1_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, iblend, igamma, i, axis
      real(wp) :: point(ndim), work_point(ndim)
      real(wp) :: weight, weight_forward, weight_backward
      real(wp) :: weight_forward2, weight_backward2
      real(wp) :: dweight_r(ndim)
      real(wp) :: numeric

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
               call prim%set_max_deriv(2)
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  do i = 1, mol%nat
                     call prim%prepare(point)
                     call prim%pou_f012_r_screened(i, weight, dweight_r=dweight_r)
                     do axis = 1, ndim
                        work_point = point
                        work_point(axis) = point(axis) + fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight_forward)
                        work_point = point
                        work_point(axis) = point(axis) + 2.0_wp*fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight_forward2)
                        work_point = point
                        work_point(axis) = point(axis) - fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight_backward)
                        work_point = point
                        work_point(axis) = point(axis) - 2.0_wp*fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight_backward2)
                        numeric = fd4_scalar(weight_forward2, weight_forward, &
                           weight_backward, weight_backward2, fd_h)
                        call check(error, dweight_r(axis), numeric, &
                           thr_abs=atol, thr_rel=rtol)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_svdw_pou_f1_r

   !> SvdW partition-of-unity spatial Hessian vs FD of POU gradient.
   subroutine test_svdw_pou_f2_rr(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      integer  :: icase, ipt, iblend, igamma, i, axisA, axisB
      real(wp) :: point(ndim), work_point(ndim)
      real(wp) :: weight, d2weight_rr(ndim, ndim)
      real(wp) :: grad_forward(ndim), grad_backward(ndim)
      real(wp) :: grad_forward2(ndim), grad_backward2(ndim)
      real(wp) :: numeric

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
               call prim%set_max_deriv(2)
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  do i = 1, mol%nat
                     call prim%prepare(point)
                     call prim%pou_f012_r_screened(i, weight, d2weight_rr=d2weight_rr)
                     do axisB = 1, ndim
                        work_point = point
                        work_point(axisB) = point(axisB) + fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight, dweight_r=grad_forward)
                        work_point = point
                        work_point(axisB) = point(axisB) + 2.0_wp*fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight, dweight_r=grad_forward2)
                        work_point = point
                        work_point(axisB) = point(axisB) - fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight, dweight_r=grad_backward)
                        work_point = point
                        work_point(axisB) = point(axisB) - 2.0_wp*fd_h
                        call prim%prepare(work_point)
                        call prim%pou_f012_r_screened(i, weight, dweight_r=grad_backward2)
                        do axisA = 1, ndim
                           numeric = fd4_scalar(grad_forward2(axisA), grad_forward(axisA), &
                              grad_backward(axisA), grad_backward2(axisA), fd_h)
                           call check(error, d2weight_rr(axisA, axisB), numeric, &
                              thr_abs=atol, thr_rel=rtol)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_svdw_pou_f2_rr

   !> SvdW partition-of-unity mixed spatial-nuclear derivative vs FD of POU gradient.
   subroutine test_svdw_pou_f2_r_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      real(wp), allocatable :: d2weight_r_rA(:, :, :)
      integer  :: icase, ipt, iblend, igamma, atomA, axisA, axis
      real(wp) :: point(ndim)
      real(wp) :: weight_dummy
      real(wp) :: grad_forward(ndim), grad_backward(ndim)
      real(wp) :: grad_forward2(ndim), grad_backward2(ndim)
      real(wp) :: numeric
      type(structure_type) :: mol_shift
      integer, allocatable :: atomic_numbers(:)
      integer :: iat
      integer, parameter :: owner_id = 1

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate(atomic_numbers)
         allocate(atomic_numbers(mol%nat))
         do iat = 1, mol%nat
            atomic_numbers(iat) = mol%num(mol%id(iat))
         end do
         if (allocated(d2weight_r_rA)) deallocate(d2weight_r_rA)
         allocate(d2weight_r_rA(ndim, ndim, mol%nat))
         do iblend = 1, n_svdw_blends
            do igamma = 1, n_svdw_gammas
               prim%screening_threshold = 0.0_wp
               call prim%new(blend_k=svdw_blend_k_values(iblend), &
                             blend_3b=svdw_gamma_values(igamma))
               call prim%update(mol, radii)
               call prim%set_max_deriv(2)
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%update(mol, radii)
                  call prim%ssd_system%update(centers_base, radii)
                  call prim%prepare(point)
                  call prim%pou_f2_r_rA_screened(owner_id, d2weight_r_rA)
                  do atomA = 1, mol%nat
                     do axisA = 1, ndim
                        centers_local = centers_base
                        centers_local(axisA, atomA) = centers_local(axisA, atomA) + fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%pou_f012_r_screened(owner_id, weight_dummy, dweight_r=grad_forward)
                        centers_local = centers_base
                        centers_local(axisA, atomA) = centers_local(axisA, atomA) + 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%pou_f012_r_screened(owner_id, weight_dummy, dweight_r=grad_forward2)
                        centers_local = centers_base
                        centers_local(axisA, atomA) = centers_local(axisA, atomA) - fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%pou_f012_r_screened(owner_id, weight_dummy, dweight_r=grad_backward)
                        centers_local = centers_base
                        centers_local(axisA, atomA) = centers_local(axisA, atomA) - 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%pou_f012_r_screened(owner_id, weight_dummy, dweight_r=grad_backward2)
                        do axis = 1, ndim
                           numeric = fd4_scalar(grad_forward2(axis), grad_forward(axis), &
                              grad_backward(axis), grad_backward2(axis), fd_h)
                           call check(error, d2weight_r_rA(axis, axisA, atomA), numeric, &
                              thr_abs=atol, thr_rel=rtol)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate(centers_base, centers_local)
      end do
   end subroutine test_svdw_pou_f2_r_rA

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

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         allocate(centers_base(ndim, mol%nat), centers_local(ndim, mol%nat))
         centers_base = mol%xyz
         if (allocated(atomic_numbers)) deallocate(atomic_numbers)
         allocate(atomic_numbers(mol%nat))
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
               do ipt = 1, size(points, 2)
                  point = points(:, ipt)
                  call prim%update(mol, radii)
                  call prim%ssd_system%update(centers_base, radii)
                  call prim%prepare(point)
                  call prim%normalized_f01_rA_screened(normalized_val, deriv_rA=deriv_rA)
                  do atom = 1, mol%nat
                     do axis = 1, ndim
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%normalized_f01_rA_screened(f_forward)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%normalized_f01_rA_screened(f_forward2)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%normalized_f01_rA_screened(f_backward)
                        centers_local = centers_base
                        centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*fd_h
                        call new(mol_shift, atomic_numbers, centers_local)
                        call prim%update(mol_shift, radii)
                        call prim%ssd_system%update(centers_local, radii)
                        call prim%prepare(point)
                        call prim%normalized_f01_rA_screened(f_backward2)
                        numeric = fd4_scalar(f_forward2, f_forward, f_backward, f_backward2, fd_h)
                        call check(error, deriv_rA(axis, atom), numeric, &
                           thr_abs=atol, thr_rel=rtol)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
         deallocate(centers_base, centers_local)
      end do
   end subroutine test_svdw_normalized_f1_rA

   !> SvdW body-order weight reduction sanity check.
   !>
   !> Verify that selecting blend_2b=1 (pure pair-mean) collapses to the
   !> average of two atom SSDs, and blend_3b=1 (pure triple-mean) collapses
   !> to the mean of three atom SSDs.
   subroutine test_svdw_body_order_scaling(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      integer, parameter :: n2 = 2, n3 = 3
      integer :: atomic_numbers(3)
      real(wp) :: centers2(ndim, n2), radii2(n2), point2(ndim), d2(n2)
      real(wp) :: centers3(ndim, n3), radii3(n3), point3(ndim), d3(n3)
      real(wp) :: lsf0, expected

      atomic_numbers = 1

      centers2 = reshape([ &
         -1.20_wp,  0.10_wp,  0.00_wp, &
          1.10_wp, -0.20_wp,  0.30_wp], [ndim, n2])
      radii2 = [0.85_wp, 1.05_wp]
      point2 = [0.35_wp, 0.45_wp, -0.25_wp]

      call new(mol, atomic_numbers(:n2), centers2)
      prim%screening_threshold = 0.0_wp
      call prim%new(blend_k=2.4_wp, blend_1b=0.0_wp, blend_2b=1.0_wp, blend_3b=0.0_wp)
      call prim%update(mol, radii2)
      call prim%set_max_deriv(0)
      call prim%ssd_system%update(centers2, radii2)
      call prim%prepare(point2)
      call prim%f0_screened(lsf0)
      d2(1) = ssd0(point2, centers2(:, 1), radii2(1))
      d2(2) = ssd0(point2, centers2(:, 2), radii2(2))
      expected = 0.5_wp*sum(d2)
      call check(error, lsf0, expected, thr_abs=atol, thr_rel=rtol)
      if (allocated(error)) return

      centers3 = reshape([ &
         -1.40_wp,  0.20_wp, -0.10_wp, &
          1.25_wp, -0.30_wp,  0.35_wp, &
          0.15_wp,  1.10_wp, -0.45_wp], [ndim, n3])
      radii3 = [0.80_wp, 1.00_wp, 0.75_wp]
      point3 = [0.10_wp, 0.25_wp, 0.40_wp]

      call new(mol, atomic_numbers(:n3), centers3)
      prim%screening_threshold = 0.0_wp
      call prim%new(blend_k=1.7_wp, blend_1b=0.0_wp, blend_2b=0.0_wp, blend_3b=1.0_wp)
      call prim%update(mol, radii3)
      call prim%set_max_deriv(0)
      call prim%ssd_system%update(centers3, radii3)
      call prim%prepare(point3)
      call prim%f0_screened(lsf0)
      d3(1) = ssd0(point3, centers3(:, 1), radii3(1))
      d3(2) = ssd0(point3, centers3(:, 2), radii3(2))
      d3(3) = ssd0(point3, centers3(:, 3), radii3(3))
      expected = sum(d3)/3.0_wp
      call check(error, lsf0, expected, thr_abs=atol, thr_rel=rtol)
   end subroutine test_svdw_body_order_scaling


end module test_cavity_drop_lsf
