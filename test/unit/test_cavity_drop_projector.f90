!> Tests for the DROP projector: candidate filtering, optimality classification,
!> saddle escape, retraction and the level-9 certified branch search
!>
!> TODO: This will be reworked with the new projection refactor
module test_cavity_drop_projector
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_projector, only: drop_projector_type
   use moist_cavity_drop_types, only: projection_workspace_type
   use moist_context, only: moist_context_type, new_context
   use moist_math_cell_grid, only: moist_cell_grid_type
   use moist_radii, only: default_cpcm_radii
   use test_helpers, only: get_test_cross, fill_legacy_radii
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_projector

   !> Solver tolerance of the projector under test
   real(wp), parameter :: PROJ_TOL = 1.0e-14_wp
   integer, parameter :: PROJ_MAXITER = 100
   !> Displacement bound of every projection (Bohr)
   real(wp), parameter :: MAX_DISP = 15.0_wp
   !> Position tolerance of the exact sphere and mirror geometries (Bohr)
   real(wp), parameter :: POS_THR = 1.0e-10_wp

contains

   !> Collect all DROP projector tests
   subroutine collect_cavity_drop_projector(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("filter_candidates_dedupe_and_slack", test_filter_candidates), &
                  new_unittest("optimality_classifies_min_max_and_non_kkt", test_optimality_classes), &
                  new_unittest("optimality_degenerate_gradient", test_optimality_degenerate_gradient), &
                  new_unittest("retraction_lands_on_surface", test_retraction), &
                  new_unittest("refine_point_escapes_saddle", test_refine_saddle), &
                  new_unittest("project_point_single_atom_levels", test_project_single_levels), &
                  new_unittest("project_point_midpoint_branches", test_project_midpoint), &
                  new_unittest("octree_budget_exhaustion_aborts", test_octree_budget), &
                  new_unittest("bad_filter_no_candidates", test_bad_filter_none, should_fail=.true.), &
                  new_unittest("bad_filter_count", test_bad_filter_count, should_fail=.true.), &
                  new_unittest("bad_filter_dimension", test_bad_filter_dimension, should_fail=.true.), &
                  new_unittest("bad_filter_slack", test_bad_filter_slack, should_fail=.true.), &
                  new_unittest("bad_strict_saddle", test_bad_strict_saddle, should_fail=.true.), &
                  new_unittest("bad_flat_optimality", test_bad_flat_optimality, should_fail=.true.), &
                  new_unittest("bad_retraction_centre", test_bad_retraction_centre, should_fail=.true.), &
                  new_unittest("bad_retraction_budget", test_bad_retraction_budget, should_fail=.true.), &
                  new_unittest("bad_escape_budget", test_bad_escape_budget, should_fail=.true.), &
                  new_unittest("bad_escape_centre", test_bad_escape_centre, should_fail=.true.), &
                  new_unittest("bad_refine_newton_budget", test_bad_refine_budget, should_fail=.true.), &
                  new_unittest("bad_octree_budget", test_bad_octree_budget, should_fail=.true.) &
                  ]

   end subroutine collect_cavity_drop_projector

   !* ================================================================================= *!
   !*                                 Candidate filtering                               *!
   !* ================================================================================= *!

   !> Duplicates are merged and branches beyond the squared-distance slack carry no weight
   subroutine test_filter_candidates(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(moist_error_type), allocatable :: failed
      logical, allocatable :: keep(:)
      real(wp) :: cand(3, 5)
      real(wp), parameter :: anchor(3) = 0.0_wp

      !> Anchor at the origin: rho^2 = 1.0, 1.0 (duplicate), 1.4, 1.6, 4.0
      cand(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
      cand(:, 2) = [1.0_wp, 1.0e-12_wp, 0.0_wp]
      cand(:, 3) = [0.0_wp, sqrt(1.4_wp), 0.0_wp]
      cand(:, 4) = [0.0_wp, 0.0_wp, sqrt(1.6_wp)]
      cand(:, 5) = [-2.0_wp, 0.0_wp, 0.0_wp]
      proj%branch_sep_cut = 1.0e-10_wp
      proj%branch_rho2_slack = 0.5_wp

      call proj%filter_candidates(cand, 5, keep, failed, anchor)
      call check(error, .not. allocated(failed), "filtering succeeds")
      if (allocated(error)) return
      call check(error, all(keep .eqv. [.true., .false., .true., .false., .false.]), &
                 "duplicate merged, rho^2 slack of 0.5 keeps 1.4 and rejects 1.6 and 4.0")
      if (allocated(error)) return

      !> A wider slack admits the branch that was out of range, and a wide
      !> separation cut merges candidates that are distinct at the default
      call proj%filter_candidates(cand, 5, keep, failed, anchor, pm_sep_cut_inp=0.7_wp)
      call check(error, all(keep .eqv. [.true., .false., .true., .true., .false.]), &
                 "slack override admits rho^2 = 1.6")
      if (allocated(error)) return
      call proj%filter_candidates(cand, 5, keep, failed, anchor, pp_sep_cut_inp=1.7_wp)
      call check(error, all(keep .eqv. [.true., .false., .false., .false., .false.]), &
                 "separation override merges everything within 1.7 Bohr of a kept point")
      if (allocated(error)) return

      !> The closest candidate defines the slack reference, whatever the list order
      cand(:, 1:2) = cand(:, 4:5)
      cand(:, 3) = [1.0_wp, 0.0_wp, 0.0_wp]
      call proj%filter_candidates(cand(:, 1:3), 3, keep, failed, anchor)
      call check(error, all(keep .eqv. [.false., .false., .true.]), &
                 "slack is measured from the closest candidate, not the first")

   end subroutine test_filter_candidates

   !* ================================================================================= *!
   !*                         Optimality classification (sphere)                        *!
   !* ================================================================================= *!

   !> On a sphere about an external anchor, the near pole is the minimum, the far
   !> pole is a constrained maximum, and points off the surface or off the
   !> stationarity condition are no KKT points
   subroutine test_optimality_classes(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda
      logical :: saddle, kkt, converged

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]

      !> Near pole
      call on_surface(proj, [rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "near pole classified without error")
      if (allocated(error)) return
      call check(error, kkt .and. converged .and. .not. saddle, "near pole is a constrained minimum")
      if (allocated(error)) return

      !> Far pole: stationary, but the displacement is maximal along the surface
      call on_surface(proj, [-rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "far pole flagged through the saddle output, no error")
      if (allocated(error)) return
      call check(error, kkt .and. saddle .and. .not. converged, "far pole is a saddle-type KKT point")
      if (allocated(error)) return

      !> Without a way to hand the saddle back, strict mode must refuse it
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, is_kkt=kkt)
      call check(error, allocated(failed), "strict mode refuses a saddle that the caller cannot handle")
      if (allocated(error)) return
      proj%strict_minimum_required = .false.
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, is_kkt=kkt)
      call check(error, .not. allocated(failed), "non-strict mode accepts it")
      if (allocated(error)) return
      proj%strict_minimum_required = .true.

      !> Off the surface
      call proj%check_constrained_optimality([rs + 0.2_wp, 0.0_wp, 0.0_wp], anchor, 1, 1, lambda, &
                                             failed, saddle, kkt, converged)
      call check(error, .not. kkt .and. .not. converged, "a point 0.2 Bohr off the surface is no KKT point")
      if (allocated(error)) return

      !> On the surface, but not the closest point: the stationarity condition fails
      call on_surface(proj, [0.0_wp, rs, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, saddle, kkt, converged)
      call check(error, .not. kkt .and. .not. converged, "the equator is no KKT point")

   end subroutine test_optimality_classes

   !> At a symmetric pinch the LSF gradient vanishes and only feasibility decides
   subroutine test_optimality_degenerate_gradient(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, half, mid(3), s_mid, grad(3), anchor(3)
      logical :: kkt

      !> Reference single-atom radius to size the pair
      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      call proj%destroy()
      deallocate (radii)

      !> Two equal atoms far enough apart that the midpoint is outside the surface
      half = 2.2_wp*rs
      call new(mol, [6, 6], reshape([-half, 0.0_wp, 0.0_wp, half, 0.0_wp, 0.0_wp], [3, 2]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      mid = 0.0_wp
      anchor = [0.0_wp, 0.0_wp, 0.5_wp]
      call proj%compute_ssd(mid)
      call proj%lsf%f012_r(lsf0=s_mid, lsf1_r=grad)
      call check(error, norm2(grad) < 1.0e-10_wp, "mirror symmetry makes the gradient vanish at the midpoint")
      if (allocated(error)) return
      call check(error, abs(s_mid) > 1.0e-3_wp, "the midpoint is off the surface")
      if (allocated(error)) return

      call proj%check_constrained_optimality(mid, anchor, 1, 1, 0.0_wp, failed, is_kkt=kkt)
      call check(error, .not. allocated(failed), "degenerate gradient does not raise")
      if (allocated(error)) return
      call check(error, .not. kkt, "an infeasible point with vanishing gradient is no KKT point")
      if (allocated(error)) return

      !> Tune the separation until the midpoint is on the surface: the pinch
      call tune_pinch(half, rs, error)
      if (allocated(error)) return
      call proj%destroy()
      deallocate (radii)
      call new(mol, [6, 6], reshape([-half, 0.0_wp, 0.0_wp, half, 0.0_wp, 0.0_wp], [3, 2]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call proj%compute_ssd(mid)
      call proj%lsf%f012_r(lsf0=s_mid, lsf1_r=grad)
      call check(error, abs(s_mid) < 1.0e-12_wp .and. norm2(grad) < 1.0e-10_wp, "fixture sits on the pinch")
      if (allocated(error)) return
      call proj%check_constrained_optimality(mid, anchor, 1, 1, 0.0_wp, failed, is_kkt=kkt)
      call check(error, kkt, "a feasible point with vanishing gradient passes as KKT, unclassified")

   contains

      !> Bisect the half-separation on the midpoint value of the pair's own LSF
      subroutine tune_pinch(h, r0, err)
         real(wp), intent(inout) :: h
         real(wp), intent(in) :: r0
         type(error_type), allocatable, intent(out) :: err

         real(wp) :: lo, hi, s, s_lo
         integer :: it

         lo = 0.7_wp*r0
         hi = 2.2_wp*r0
         s_lo = midpoint_value(lo, err)
         if (allocated(err)) return
         do it = 1, 60
            h = 0.5_wp*(lo + hi)
            s = midpoint_value(h, err)
            if (allocated(err)) return
            if (s*s_lo > 0.0_wp) then
               lo = h
            else
               hi = h
            end if
         end do
         h = 0.5_wp*(lo + hi)
      end subroutine tune_pinch

      !> LSF value at the midpoint of a pair at half-separation `h`
      function midpoint_value(h, err) result(s)
         real(wp), intent(in) :: h
         type(error_type), allocatable, intent(out) :: err
         real(wp) :: s

         type(drop_projector_type) :: p
         type(structure_type) :: m
         real(wp), allocatable :: rad(:)
         real(wp) :: g(3)

         s = 0.0_wp
         call new(m, [6, 6], reshape([-h, 0.0_wp, 0.0_wp, h, 0.0_wp, 0.0_wp], [3, 2]))
         call build_projector(p, m, rad, err)
         if (allocated(err)) return
         call p%compute_ssd(mid)
         call p%lsf%f012_r(lsf0=s, lsf1_r=g)
      end function midpoint_value

   end subroutine test_optimality_degenerate_gradient

   !* ================================================================================= *!
   !*                            Retraction, escape, refinement                         *!
   !* ================================================================================= *!

   !> Retraction returns a point with |S| < 1e-12 and a multiplier that makes it stationary
   subroutine test_retraction(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda, s, g(3)
      logical :: kkt, converged

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]

      !> Radial offset retracts to the near pole, which is then a KKT minimum
      call proj%retract_to_surface([rs + 0.4_wp, 0.0_wp, 0.0_wp], anchor, 1, r, lambda, failed)
      call check(error, .not. allocated(failed), "retraction succeeds")
      if (allocated(error)) return
      call proj%compute_ssd(r)
      call proj%lsf%f012_r(lsf0=s, lsf1_r=g)
      call check(error, abs(s) < 1.0e-12_wp, "retracted point is on the surface")
      if (allocated(error)) return
      call check(error, norm2(r - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "the radial offset retracts radially")
      if (allocated(error)) return
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, is_kkt=kkt, is_converged=converged)
      call check(error, kkt .and. converged, "the retraction multiplier is the stationarity multiplier")
      if (allocated(error)) return

      !> A skew start lands on the sphere as well
      call proj%retract_to_surface([-1.3_wp, 2.1_wp, 0.7_wp], anchor, 1, r, lambda, failed)
      call check(error, .not. allocated(failed), "skew retraction succeeds")
      if (allocated(error)) return
      call check(error, abs(norm2(r) - rs) < POS_THR, "skew start lands on the sphere")
      if (allocated(error)) return

   end subroutine test_retraction

   !> Newton converges onto the far pole, the classification sees a saddle, and
   !> the refinement escapes and re-polishes to the near pole
   subroutine test_refine_saddle(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), xl(3), xu(3), seed(3), gp(3), lambda, rho, phi, nrm(3)

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      xl = anchor - MAX_DISP
      xu = anchor + MAX_DISP

      !> A seed next to the far pole
      seed = [-rs - 0.05_wp, 0.02_wp, 0.0_wp]
      call proj%refine_point(anchor, seed, xl, xu, MAX_DISP, 1, 1, failed, gp, lambda, rho, phi, nrm)
      call check(error, .not. allocated(failed), "refinement succeeds")
      if (allocated(error)) return
      call check(error, norm2(gp - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "refinement ends on the near pole")
      if (allocated(error)) return
      call check(error, abs(rho - 3.0_wp) < POS_THR, "reported displacement is the anchor height above the sphere")
      if (allocated(error)) return
      call check(error, norm2(nrm - [1.0_wp, 0.0_wp, 0.0_wp]) < POS_THR, "unit normal points along the axis")
      if (allocated(error)) return

      !> The same seed without the check stays on the maximum
      proj%verify_optimality = .false.
      call proj%refine_point(anchor, seed, xl, xu, MAX_DISP, 1, 1, failed, gp)
      call check(error, .not. allocated(failed), "unverified refinement succeeds")
      if (allocated(error)) return
      call check(error, norm2(gp - [-rs, 0.0_wp, 0.0_wp]) < POS_THR, "without the check the maximum is kept")

   end subroutine test_refine_saddle

   !* ================================================================================= *!
   !*                                 Top-level projection                              *!
   !* ================================================================================= *!

   !> Every level lands on the near pole of the sphere with one branch
   subroutine test_project_single_levels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(projection_workspace_type) :: work
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), rho(:)
      real(wp) :: rs, anchor(3)
      integer, parameter :: levels(4) = [2, 3, 7, 9]
      integer :: il, n

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      call work%init(error=failed)
      anchor = [rs + 3.0_wp, 1.0_wp, -0.5_wp]

      do il = 1, size(levels)
         call proj%project_point(anchor, pts, n, 1, 1, proj_level=levels(il), initial_guess=anchor, &
                                 error=failed, rho_out=rho, work=work)
         call check(error, .not. allocated(failed), "projection succeeds")
         if (allocated(error)) return
         call check(error, n, 1, "a convex surface has one branch")
         if (allocated(error)) return
         call check(error, abs(norm2(pts(:, 1)) - rs) < POS_THR, "projected point is on the sphere")
         if (allocated(error)) return
         call check(error, norm2(pts(:, 1)/rs - anchor/norm2(anchor)) < POS_THR, &
                    "projection is along the line to the centre")
         if (allocated(error)) return
         call check(error, abs(rho(1) - (norm2(anchor) - rs)) < POS_THR, "displacement is the height above the sphere")
         if (allocated(error)) return
      end do

   end subroutine test_project_single_levels

   !> Between two equal atoms the gradient vanishes by symmetry; the two mirror
   !> branches are found by the multistart and the certified search alike
   subroutine test_project_midpoint(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(projection_workspace_type) :: work
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), rho(:)
      real(wp) :: rs, half, anchor(3), s, g(3)
      integer, parameter :: levels(4) = [2, 3, 7, 9]
      integer :: il, n, i

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      call proj%destroy()
      deallocate (radii)

      half = 2.2_wp*rs
      call new(mol, [6, 6], reshape([-half, 0.0_wp, 0.0_wp, half, 0.0_wp, 0.0_wp], [3, 2]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call work%init(error=failed)
      anchor = 0.0_wp

      do il = 1, size(levels)
         call proj%project_point(anchor, pts, n, 1, 1, proj_level=levels(il), initial_guess=anchor, &
                                 error=failed, rho_out=rho, work=work)
         call check(error, .not. allocated(failed), "projection from the symmetric point succeeds")
         if (allocated(error)) return
         do i = 1, n
            call proj%compute_ssd(pts(:, i))
            call proj%lsf%f012_r(lsf0=s, lsf1_r=g)
            call check(error, abs(s) < 1.0e-12_wp, "every branch is on the surface")
            if (allocated(error)) return
         end do
         if (levels(il) >= 7) then
            call check(error, n, 2, "both mirror branches are found")
            if (allocated(error)) return
            call check(error, norm2(pts(:, 1) + pts(:, 2)) < POS_THR, "the branches are mirror images")
            if (allocated(error)) return
            call check(error, abs(rho(1) - rho(2)) < POS_THR, "and equally far from the anchor")
            if (allocated(error)) return
         end if
      end do

   end subroutine test_project_midpoint

   !> A search that runs out of boxes aborts the build instead of dropping the anchor
   subroutine test_octree_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(projection_workspace_type) :: work
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rs, anchor(3)
      integer :: n

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) return
      call surface_radius(proj, rs, error)
      if (allocated(error)) return
      call work%init(error=failed)
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]

      proj%octree_max_boxes = 16
      call proj%project_point(anchor, pts, n, 1, 1, proj_level=9, error=failed, work=work)
      call check(error, allocated(proj%abort_error), "an exhausted box budget is raised as a build-aborting error")
      if (allocated(error)) return
      call check(error, index(proj%abort_error%message, "budget") > 0, "and names the budget")
      if (allocated(error)) return
      call check(error, n, 0, "and yields no points")

   end subroutine test_octree_budget

   !* ================================================================================= *!
   !*                                    Error paths                                    *!
   !* ================================================================================= *!

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(moist_error_type), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> An empty candidate list is refused
   subroutine test_bad_filter_none(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(moist_error_type), allocatable :: failed
      logical, allocatable :: keep(:)
      real(wp) :: cand(3, 2)
      real(wp), parameter :: anchor(3) = 0.0_wp

      cand = 1.0_wp
      cand(:, 2) = 2.0_wp
      call proj%filter_candidates(cand(:, 1:0), 0, keep, failed, anchor)
      call expect_error(error, failed, "filter_candidates: no candidates")
   end subroutine test_bad_filter_none

   !> A candidate count that disagrees with the array is refused
   subroutine test_bad_filter_count(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(moist_error_type), allocatable :: failed
      logical, allocatable :: keep(:)
      real(wp) :: cand(3, 2)
      real(wp), parameter :: anchor(3) = 0.0_wp

      cand = 1.0_wp
      cand(:, 2) = 2.0_wp
      call proj%filter_candidates(cand, 1, keep, failed, anchor)
      call expect_error(error, failed, "inconsistent candidate count")
   end subroutine test_bad_filter_count

   !> Candidates that are not three-dimensional are refused
   subroutine test_bad_filter_dimension(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(moist_error_type), allocatable :: failed
      logical, allocatable :: keep(:)
      real(wp) :: cand(3, 2)
      real(wp), parameter :: anchor(3) = 0.0_wp

      cand = 1.0_wp
      cand(:, 2) = 2.0_wp
      call proj%filter_candidates(cand(1:2, :), 2, keep, failed, anchor)
      call expect_error(error, failed, "leading dimension 3")
   end subroutine test_bad_filter_dimension

   !> A negative slack leaves even the closest candidate outside it and is refused
   subroutine test_bad_filter_slack(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(moist_error_type), allocatable :: failed
      logical, allocatable :: keep(:)
      real(wp) :: cand(3, 2)
      real(wp), parameter :: anchor(3) = 0.0_wp

      cand = 1.0_wp
      cand(:, 2) = 2.0_wp
      call proj%filter_candidates(cand, 2, keep, failed, anchor, pm_sep_cut_inp=-1.0_wp)
      call expect_error(error, failed, "no unique candidates after filtering")
   end subroutine test_bad_filter_slack

   !> Strict mode refuses a saddle that the caller cannot handle
   subroutine test_bad_strict_saddle(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda
      logical :: kkt

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call on_surface(proj, [-rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, is_kkt=kkt)
      call expect_error(error, failed, "Saddle point detected but strict minimum required")
   end subroutine test_bad_strict_saddle

   !> A surface flat to third order is reported as unreliable
   subroutine test_bad_flat_optimality(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda
      logical :: saddle, kkt, converged

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = 0.0_wp
      call on_surface(proj, [rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call proj%check_constrained_optimality(r, anchor, 1, 1, lambda, failed, saddle, kkt, converged)
      call expect_error(error, failed, "Degenerate to third order")
   end subroutine test_bad_flat_optimality

   !> The centre has no gradient to retract along
   subroutine test_bad_retraction_centre(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call proj%retract_to_surface([0.0_wp, 0.0_wp, 0.0_wp], anchor, 1, r, lambda, failed)
      call expect_error(error, failed, "Retraction failed: degenerate LSF gradient")
   end subroutine test_bad_retraction_centre

   !> One Newton step cannot retract to the surface from 6 Bohr
   subroutine test_bad_retraction_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      proj%retract_max_iter = 1
      call proj%retract_to_surface([rs + 6.0_wp, 0.5_wp, 0.0_wp], anchor, 1, r, lambda, failed)
      call expect_error(error, failed, "Retraction to surface did not converge")
   end subroutine test_bad_retraction_budget

   !> An exhausted iteration budget stops the Riemannian escape
   subroutine test_bad_escape_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r0(3), r(3), lambda0, lambda

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call on_surface(proj, [-rs, 0.0_wp, 0.0_wp], anchor, r0, lambda0, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      proj%riemann_alpha_init = 1.0_wp
      proj%riemann_max_iter = 1
      call proj%riemannian_newton_escape(r0, anchor, 1, lambda0, r, lambda, failed)
      call expect_error(error, failed, "Riemannian Newton did not converge within max iterations")
   end subroutine test_bad_escape_budget

   !> At the centre there is no surface normal to define a tangent plane
   subroutine test_bad_escape_centre(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call proj%riemannian_newton_escape([0.0_wp, 0.0_wp, 0.0_wp], anchor, 1, 0.0_wp, r, lambda, failed)
      call expect_error(error, failed, "Riemannian escape: degenerate LSF gradient")
   end subroutine test_bad_escape_centre

   !> A starved Newton budget is an error, not a wrong point
   subroutine test_bad_refine_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), xl(3), xu(3), seed(3), gp(3)

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      xl = anchor - MAX_DISP
      xu = anchor + MAX_DISP
      seed = [-rs - 0.05_wp, 0.02_wp, 0.0_wp]
      proj%newton_max_iter = 1
      call proj%refine_point(anchor, seed, xl, xu, MAX_DISP, 1, 1, failed, gp)
      call expect_error(error, failed, "Maximum iterations reached")
   end subroutine test_bad_refine_budget

   !> An exhausted octree box budget is raised as a build-aborting error
   subroutine test_bad_octree_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_projector_type) :: proj
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      type(projection_workspace_type) :: work
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :)
      real(wp) :: rs, anchor(3)
      integer :: n

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      ! Setup failures return cleanly, so the test is reported as failed
      call build_projector(proj, mol, radii, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call surface_radius(proj, rs, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call work%init(error=failed)
      if (allocated(failed)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]

      proj%octree_max_boxes = 16
      call proj%project_point(anchor, pts, n, 1, 1, proj_level=9, error=failed, work=work)
      call expect_error(error, proj%abort_error, "budget")
   end subroutine test_bad_octree_budget

   !* ================================================================================= *!
   !*                                      Fixtures                                     *!
   !* ================================================================================= *!

   !> Build a bare projector for a structure, bound to its atoms and screening grid
   !>
   !> @param[inout] proj  Projector to initialize
   !> @param[in]    mol   Structure, must outlive the projector
   !> @param[out]   radii Legacy radii of the atoms
   !> @param[out]   error Error handle
   subroutine build_projector(proj, mol, radii, error, blend_k, branch_weight_s)
      !> Projector to initialize
      type(drop_projector_type), intent(inout) :: proj
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii
      real(wp), allocatable, intent(out) :: radii(:)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), intent(in), optional :: blend_k, branch_weight_s
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(moist_cell_grid_type) :: grid
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: r_eff(:)
      real(wp) :: k_blend
      integer :: i

      k_blend = 2.5_wp
      if (present(blend_k)) k_blend = blend_k
      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return
      call param%new(nleb=50, tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=7, error=failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k_blend, blend_3b=1.0_wp))
      if (present(branch_weight_s)) then
         param%branch_weight_s = branch_weight_s
         call param%compute_derived(failed)
      end if
      lsf%screening_threshold = param%screening_threshold

      allocate (r_eff(mol%nat))
      do i = 1, mol%nat
         r_eff(i) = radii(i) + lsf%neighbor_cutoff(radii(i))
      end do
      call grid%build(mol%xyz, r_eff, full_scan_below=100)

      call proj%init(param, lsf, branch_sep_cut=param%branch_sep_cut, tol=param%proj_tol, maxiter=param%proj_maxiter)
      call proj%init_primitives(mol, radii, grid)

   end subroutine build_projector

   !> Radius of the isolated-atom surface, found by projecting a far anchor onto it
   !>
   !> @param[inout] proj Single-atom projector at the origin
   !> @param[out]   rs   Radius of the zero level set
   !> @param[out]   error Error handle
   subroutine surface_radius(proj, rs, error)
      !> Single-atom projector
      type(drop_projector_type), intent(inout) :: proj
      !> Radius
      real(wp), intent(out) :: rs
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(projection_workspace_type) :: work
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :)
      integer :: n

      rs = 0.0_wp
      call work%init(error=failed)
      call proj%project_point([8.0_wp, 0.0_wp, 0.0_wp], pts, n, 1, 1, proj_level=2, error=failed, work=work)
      if (allocated(failed)) then
         call test_failed(error, "surface radius projection failed: "//failed%message)
         return
      end if
      rs = norm2(pts(:, 1))
      call check(error, n, 1, "single projection for the radius")
      if (allocated(error)) return
      call check(error, abs(pts(1, 1) - rs) < 1.0e-12_wp, "radius probe lies on the axis")

   end subroutine surface_radius

   !> A point on the surface near `guess`, with the multiplier that makes it stationary
   !>
   !> @param[inout] proj   Projector
   !> @param[in]    guess  Approximate position
   !> @param[in]    anchor Anchor of the projection
   !> @param[out]   r      Surface point
   !> @param[out]   lambda Stationarity multiplier of the anchor objective
   !> @param[out]   error  Error handle
   subroutine on_surface(proj, guess, anchor, r, lambda, error)
      !> Projector
      type(drop_projector_type), intent(inout) :: proj
      !> Approximate position
      real(wp), intent(in) :: guess(3)
      !> Anchor
      real(wp), intent(in) :: anchor(3)
      !> Surface point
      real(wp), intent(out) :: r(3)
      !> Multiplier
      real(wp), intent(out) :: lambda
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: failed

      call proj%retract_to_surface(guess, anchor, 1, r, lambda, failed)
      if (allocated(failed)) call test_failed(error, "retraction of the fixture point failed: "//failed%message)

   end subroutine on_surface

   !> Build the five-carbon cross at a soft blend with a branch-friendly softmax
   !>
   !> @param[inout] cav   Cavity to build
   !> @param[inout] ctx   Borrowed run context
   !> @param[in]    mol   Cross structure
   !> @param[in]    level Projection level
   !> @param[out]   error Error handle
   subroutine build_cross(cav, ctx, mol, level, error)
      !> Cavity to build
      type(cavity_type_drop), intent(inout) :: cav
      !> Borrowed run context
      type(moist_context_type), intent(inout), target :: ctx
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Projection level
      integer, intent(in) :: level
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_svdw_type) :: template
      type(moist_error_type), allocatable :: failed

      call template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.0_wp, blend_3b=1.0_wp))
      call new_cavity_drop(cav, radius_model=default_cpcm_radii(), lsf_model=template, error=failed, &
                           param=moist_cavity_drop_parameters_type(num_leb=14, do_fine=.false., tolerance=PROJ_TOL, &
                                                                   proj_maxiter=PROJ_MAXITER, proj_level=level, &
                                                                   wleb_prune_level=4), ctx=ctx)
      if (allocated(failed)) then
         call test_failed(error, "cavity setup failed: "//failed%message)
         return
      end if
      cav%param%branch_weight_s = 0.5_wp
      call cav%param%compute_derived(failed)
      if (allocated(failed)) then
         call test_failed(error, "derived parameters failed: "//failed%message)
         return
      end if
      call cav%branch_weight%init(0.5_wp)
      call cav%update(mol, error=failed)
      if (allocated(failed)) call test_failed(error, "cavity build failed: "//failed%message)

   end subroutine build_cross

end module test_cavity_drop_projector
