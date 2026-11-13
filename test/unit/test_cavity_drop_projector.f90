!> Tests for the per-point DROP projection pieces: candidate deduplication,
!> optimality classification, saddle escape, retraction, and the single-anchor
!> projection entry at every level, including the octree retry of a failed
!> local solve
module test_cavity_drop_projector
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_math_packed, only: packed_full_contract2
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_converged, point_accepted, point_failed_iterations
   use moist_cavity_drop_projection_newton, only: newton_config_type, newton_scratch_type, &
                                                 newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_config_type, polish_batch, &
                                                 check_constrained_optimality, riemannian_newton_escape, &
                                                 retract_to_surface
   use moist_cavity_drop_projection_finalize, only: dedupe_candidates, cp_jacobian_weight
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_projection_worker, only: projection_worker_type
   use test_helpers, only: fill_legacy_radii
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_projector

   !> Solver tolerance of the projection under test
   real(wp), parameter :: PROJ_TOL = 1.0e-14_wp
   integer, parameter :: PROJ_MAXITER = 100
   !> Position tolerance of the exact sphere and mirror geometries (Bohr)
   real(wp), parameter :: POS_THR = 1.0e-10_wp

   !> Bare level set fixture: one structure, its LSF and the projection settings
   type :: fixture_type
      !> Structure
      type(structure_type) :: mol
      !> Legacy radii of the atoms
      real(wp), allocatable :: radii(:)
      !> Level set bound to the structure
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Projection settings
      type(moist_cavity_drop_parameters_type) :: param
      !> Polish settings
      type(polish_config_type) :: pcfg
   end type fixture_type

contains

   !> Collect all DROP projector tests
   subroutine collect_cavity_drop_projector(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dedupe_candidates_and_slack", test_dedupe_candidates), &
                  new_unittest("optimality_classifies_min_max_and_non_kkt", test_optimality_classes), &
                  new_unittest("optimality_degenerate_gradient", test_optimality_degenerate_gradient), &
                  new_unittest("retraction_lands_on_surface", test_retraction), &
                  new_unittest("polish_escapes_saddle", test_refine_saddle), &
                  new_unittest("project_anchor_single_atom_levels", test_project_single_levels), &
                  new_unittest("project_anchor_midpoint_branches", test_project_midpoint), &
                  new_unittest("project_anchor_failed_local_octree", test_project_failed_octree), &
                  new_unittest("project_anchor_failed_local_nearest_foot", test_project_nearest_foot), &
                  new_unittest("project_anchor_nearest_foot_confirmed", test_project_nearest_confirmed), &
                  new_unittest("bad_strict_saddle", test_bad_strict_saddle, should_fail=.true.), &
                  new_unittest("bad_flat_optimality", test_bad_flat_optimality, should_fail=.true.), &
                  new_unittest("bad_retraction_centre", test_bad_retraction_centre, should_fail=.true.), &
                  new_unittest("bad_retraction_budget", test_bad_retraction_budget, should_fail=.true.), &
                  new_unittest("bad_escape_budget", test_bad_escape_budget, should_fail=.true.), &
                  new_unittest("bad_escape_centre", test_bad_escape_centre, should_fail=.true.), &
                  new_unittest("bad_newton_budget", test_bad_newton_budget, should_fail=.true.), &
                  new_unittest("bad_octree_budget", test_bad_octree_budget, should_fail=.true.), &
                  new_unittest("bad_inflection_strict", test_bad_inflection_strict, should_fail=.true.), &
                  new_unittest("bad_retraction_critical_point", test_bad_retraction_critical, should_fail=.true.), &
                  new_unittest("bad_worker_level", test_bad_worker_level, should_fail=.true.), &
                  new_unittest("bad_worker_level_2", test_bad_worker_level_2, should_fail=.true.), &
                  new_unittest("bad_nearest_count", test_bad_nearest_count, should_fail=.true.), &
                  new_unittest("bad_singular_tangent_jacobian", test_bad_singular_jacobian, should_fail=.true.) &
                  ]

   end subroutine collect_cavity_drop_projector

   !* ================================================================================= *!
   !*                                 Candidate filtering                               *!
   !* ================================================================================= *!

   !> Duplicates are merged and branches beyond the squared-distance slack carry no weight
   subroutine test_dedupe_candidates(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      logical :: keep(5)
      real(wp) :: cand(3, 5)
      real(wp), parameter :: anchor(3) = 0.0_wp

      !> Anchor at the origin: rho^2 = 1.0, 1.0 (duplicate), 1.4, 1.6, 4.0
      cand(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
      cand(:, 2) = [1.0_wp, 1.0e-12_wp, 0.0_wp]
      cand(:, 3) = [0.0_wp, sqrt(1.4_wp), 0.0_wp]
      cand(:, 4) = [0.0_wp, 0.0_wp, sqrt(1.6_wp)]
      cand(:, 5) = [-2.0_wp, 0.0_wp, 0.0_wp]

      call dedupe_candidates(cand, anchor, 1.0e-10_wp, 0.5_wp, keep)
      call check(error, all(keep .eqv. [.true., .false., .true., .false., .false.]), &
                 "duplicate merged, rho^2 slack of 0.5 keeps 1.4 and rejects 1.6 and 4.0")
      if (allocated(error)) return

      !> A wider slack admits the branch that was out of range, and a wide
      !> separation cut merges candidates that are distinct at the default
      call dedupe_candidates(cand, anchor, 1.0e-10_wp, 0.7_wp, keep)
      call check(error, all(keep .eqv. [.true., .false., .true., .true., .false.]), &
                 "slack override admits rho^2 = 1.6")
      if (allocated(error)) return
      call dedupe_candidates(cand, anchor, 1.7_wp, 0.5_wp, keep)
      call check(error, all(keep .eqv. [.true., .false., .false., .false., .false.]), &
                 "separation override merges everything within 1.7 Bohr of a kept point")
      if (allocated(error)) return

      !> The closest candidate defines the slack reference, whatever the list order
      call dedupe_candidates(cand(:, [4, 5, 1]), anchor, 1.0e-10_wp, 0.5_wp, keep(1:3))
      call check(error, all(keep(1:3) .eqv. [.false., .false., .true.]), &
                 "slack is measured from the closest candidate, not the first")

   end subroutine test_dedupe_candidates

   !* ================================================================================= *!
   !*                         Optimality classification (sphere)                        *!
   !* ================================================================================= *!

   !> On a sphere about an external anchor, the near pole is the minimum, the far
   !> pole is a constrained maximum, and points off the surface or off the
   !> stationarity condition are no KKT points
   subroutine test_optimality_classes(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(polish_config_type) :: loose
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda
      logical :: saddle, kkt, converged

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]

      !> Near pole
      call on_surface(fx, [rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "near pole classified without error")
      if (allocated(error)) return
      call check(error, kkt .and. converged .and. .not. saddle, "near pole is a constrained minimum")
      if (allocated(error)) return

      !> Far pole: stationary, but the displacement is maximal along the surface
      call on_surface(fx, [-rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "far pole flagged through the saddle output, no error")
      if (allocated(error)) return
      call check(error, kkt .and. saddle .and. .not. converged, "far pole is a saddle-type KKT point")
      if (allocated(error)) return

      !> Without a way to hand the saddle back, strict mode must refuse it
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, is_kkt=kkt)
      call check(error, allocated(failed), "strict mode refuses a saddle that the caller cannot handle")
      if (allocated(error)) return
      loose = fx%pcfg
      loose%strict_minimum_required = .false.
      call check_constrained_optimality(loose, fx%lsf, r, anchor, lambda, failed, is_kkt=kkt)
      call check(error, .not. allocated(failed), "non-strict mode accepts it")
      if (allocated(error)) return

      !> Off the surface
      call check_constrained_optimality(fx%pcfg, fx%lsf, [rs + 0.2_wp, 0.0_wp, 0.0_wp], anchor, lambda, &
                                        failed, saddle, kkt, converged)
      call check(error, .not. kkt .and. .not. converged, "a point 0.2 Bohr off the surface is no KKT point")
      if (allocated(error)) return

      !> On the surface, but not the closest point: the stationarity condition fails
      call on_surface(fx, [0.0_wp, rs, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, saddle, kkt, converged)
      call check(error, .not. kkt .and. .not. converged, "the equator is no KKT point")

   end subroutine test_optimality_classes

   !> At a symmetric pinch the LSF gradient vanishes and only feasibility decides
   subroutine test_optimality_degenerate_gradient(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, half, mid(3), s_mid, grad(3), anchor(3)
      logical :: kkt

      !> Reference single-atom radius to size the pair
      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return

      !> Two equal atoms far enough apart that the midpoint is outside the surface
      half = 2.2_wp*rs
      call build_fixture(fx, carbon_pair(half), error)
      if (allocated(error)) return
      mid = 0.0_wp
      anchor = [0.0_wp, 0.0_wp, 0.5_wp]
      call lsf_value(fx, mid, s_mid, grad, error)
      if (allocated(error)) return
      call check(error, norm2(grad) < 1.0e-10_wp, "mirror symmetry makes the gradient vanish at the midpoint")
      if (allocated(error)) return
      call check(error, abs(s_mid) > 1.0e-3_wp, "the midpoint is off the surface")
      if (allocated(error)) return

      call check_constrained_optimality(fx%pcfg, fx%lsf, mid, anchor, 0.0_wp, failed, is_kkt=kkt)
      call check(error, .not. allocated(failed), "degenerate gradient does not raise")
      if (allocated(error)) return
      call check(error, .not. kkt, "an infeasible point with vanishing gradient is no KKT point")
      if (allocated(error)) return

      !> Tune the separation until the midpoint is on the surface: the pinch
      call tune_pinch(half, rs, error)
      if (allocated(error)) return
      call build_fixture(fx, carbon_pair(half), error)
      if (allocated(error)) return
      call lsf_value(fx, mid, s_mid, grad, error)
      if (allocated(error)) return
      call check(error, abs(s_mid) < 1.0e-12_wp .and. norm2(grad) < 1.0e-10_wp, "fixture sits on the pinch")
      if (allocated(error)) return
      call check_constrained_optimality(fx%pcfg, fx%lsf, mid, anchor, 0.0_wp, failed, is_kkt=kkt)
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

         type(fixture_type) :: p
         real(wp) :: g(3)

         s = 0.0_wp
         call build_fixture(p, carbon_pair(h), err)
         if (allocated(err)) return
         call lsf_value(p, mid, s, g, err)
      end function midpoint_value

   end subroutine test_optimality_degenerate_gradient

   !* ================================================================================= *!
   !*                            Retraction, escape, refinement                         *!
   !* ================================================================================= *!

   !> Retraction returns a point with |S| < 1e-12 and a multiplier that makes it stationary
   subroutine test_retraction(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda, s, g(3)
      logical :: kkt, converged

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]

      !> Radial offset retracts to the near pole, which is then a KKT minimum
      call retract_to_surface(fx%pcfg, fx%lsf, [rs + 0.4_wp, 0.0_wp, 0.0_wp], anchor, r, lambda, failed)
      call check(error, .not. allocated(failed), "retraction succeeds")
      if (allocated(error)) return
      call lsf_value(fx, r, s, g, error)
      if (allocated(error)) return
      call check(error, abs(s) < 1.0e-12_wp, "retracted point is on the surface")
      if (allocated(error)) return
      call check(error, norm2(r - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "the radial offset retracts radially")
      if (allocated(error)) return
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, is_kkt=kkt, &
                                        is_converged=converged)
      call check(error, kkt .and. converged, "the retraction multiplier is the stationarity multiplier")
      if (allocated(error)) return

      !> A skew start lands on the sphere as well
      call retract_to_surface(fx%pcfg, fx%lsf, [-1.3_wp, 2.1_wp, 0.7_wp], anchor, r, lambda, failed)
      call check(error, .not. allocated(failed), "skew retraction succeeds")
      if (allocated(error)) return
      call check(error, abs(norm2(r) - rs) < POS_THR, "skew start lands on the sphere")

   end subroutine test_retraction

   !> Newton converges onto the far pole, the polish sees a saddle, escapes and
   !> re-polishes to the near pole
   subroutine test_refine_saddle(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_config_type) :: ncfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), seed(3)
      integer :: icol

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      ncfg = newton_config(fx)

      !> A seed on the axis behind the far pole: off the axis the merit line
      !> search already walks Newton away from the maximum
      seed = [-rs - 0.05_wp, 0.0_wp, 0.0_wp]
      call batch%init(4, failed)
      call batch%begin_run(1, anchor, 1)
      call batch%append(seed, 0.0_wp, icol)
      call newton_solve_batch(batch, fx%lsf, ncfg, scratch, failed)
      call check(error, .not. allocated(failed), "Newton succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges")
      if (allocated(error)) return

      !> Without the check the maximum is what Newton found
      call check(error, norm2(batch%x(:, 1) - [-rs, 0.0_wp, 0.0_wp]) < POS_THR, &
                 "Newton alone stops on the maximum")
      if (allocated(error)) return

      call polish_batch(batch, fx%lsf, fx%pcfg, ncfg, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_accepted, "escaped point verified")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "polish ends on the near pole")
      if (allocated(error)) return
      call check(error, abs(batch%rho(1) - 3.0_wp) < POS_THR, &
                 "reported displacement is the anchor height above the sphere")
      if (allocated(error)) return
      call check(error, norm2(batch%normal(:, 1) - [1.0_wp, 0.0_wp, 0.0_wp]) < POS_THR, &
                 "unit normal points along the axis")

   end subroutine test_refine_saddle

   !* ================================================================================= *!
   !*                                 Top-level projection                              *!
   !* ================================================================================= *!

   !> Every level lands on the near pole of the sphere with one branch
   subroutine test_project_single_levels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), rho(:)
      real(wp) :: rs, anchor(3)
      integer :: il, level, n
      !> Every projection level
      integer, parameter :: levels(3) = [1, 3, 4]

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 1.0_wp, -0.5_wp]

      do il = 1, size(levels)
         level = levels(il)
         call new_worker(fx, level, worker, error)
         if (allocated(error)) return
         call worker%project_anchor(anchor, 1, n, failed, feet=pts, rho=rho)
         call check(error, .not. allocated(failed), "projection succeeds")
         if (allocated(error)) return
         call check(error, n, 1, "a convex surface has one branch")
         if (allocated(error)) return
         call check(error, abs(norm2(pts(:, 1)) - rs) < POS_THR, "projected point is on the sphere")
         if (allocated(error)) return
         call check(error, norm2(pts(:, 1)/rs - anchor/norm2(anchor)) < POS_THR, &
                    "projection is along the line to the centre")
         if (allocated(error)) return
         call check(error, abs(rho(1) - (norm2(anchor) - rs)) < POS_THR, &
                    "displacement is the height above the sphere")
         if (allocated(error)) return
         if (level == 1) then
            ! Converged local foot and no global search: no seeds at all
            call check(error, worker%n_seeds, 0, "a converged anchor proposes no seeds")
            if (allocated(error)) return
            call check(error, worker%n_escalated, 0, "a converged anchor is not escalated")
            if (allocated(error)) return
         end if
      end do

   end subroutine test_project_single_levels

   !> Between two equal atoms the gradient vanishes by symmetry; the local solve
   !> fails and every level's octree retry finds the two mirror branches
   subroutine test_project_midpoint(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), rho(:)
      real(wp) :: rs, half, anchor(3), s, g(3)
      integer :: il, level, n, i
      !> Every projection level
      integer, parameter :: levels(3) = [1, 3, 4]

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return

      half = 2.2_wp*rs
      call build_fixture(fx, carbon_pair(half), error)
      if (allocated(error)) return
      anchor = 0.0_wp

      do il = 1, size(levels)
         level = levels(il)
         call new_worker(fx, level, worker, error)
         if (allocated(error)) return
         call worker%project_anchor(anchor, 1, n, failed, feet=pts, rho=rho)
         call check(error, .not. allocated(failed), "projection from the symmetric point succeeds")
         if (allocated(error)) return
         do i = 1, n
            call lsf_value(fx, pts(:, i), s, g, error)
            if (allocated(error)) return
            call check(error, abs(s) < 1.0e-12_wp, "every branch is on the surface")
            if (allocated(error)) return
         end do
         call check(error, n, 2, "both mirror branches are found")
         if (allocated(error)) return
         call check(error, norm2(pts(:, 1) + pts(:, 2)) < POS_THR, "the branches are mirror images")
         if (allocated(error)) return
         call check(error, abs(rho(1) - rho(2)) < POS_THR, "and equally far from the anchor")
         if (allocated(error)) return
      end do

   end subroutine test_project_midpoint

   !> A failed local solve at level 1 is retried from the coarse octree of
   !> level 3: same seeds, same feet
   !>
   !> - midpoint of a carbon pair, where the vanishing gradient stalls Newton
   subroutine test_project_failed_octree(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), ref(:, :)
      real(wp) :: rs, anchor(3)
      integer :: n, n_ref, seeds, seeds_ref
      logical :: lost

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      call build_fixture(fx, carbon_pair(2.2_wp*rs), error)
      if (allocated(error)) return
      anchor = 0.0_wp

      call new_worker(fx, 3, worker, error)
      if (allocated(error)) return
      call worker%project_anchor(anchor, 1, n_ref, failed, feet=ref, n_seeds=seeds_ref)
      call check(error, .not. allocated(failed), "level-3 reference projection succeeds")
      if (allocated(error)) return

      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      call worker%project_anchor(anchor, 1, n, failed, feet=pts, n_seeds=seeds, lost=lost)
      call check(error, .not. allocated(failed), "projection succeeds")
      if (allocated(error)) return
      call check(error, worker%n_escalated, 1, "the failed local solve is escalated")
      if (allocated(error)) return
      call check(error, seeds, seeds_ref, "the retry proposes the level-3 octree seeds")
      if (allocated(error)) return
      call check(error, .not. lost, "the anchor keeps a foot")
      if (allocated(error)) return
      call check(error, n, n_ref, "the level-3 branch count")
      if (allocated(error)) return
      call check(error, maxval(abs(pts(:, 1:n) - ref(:, 1:n_ref))) < POS_THR, "the level-3 feet")

   end subroutine test_project_failed_octree

   !> A failed local solve gets the foot closest to the anchor
   !>
   !> - centre of two small atoms on x and two large ones on y, where the
   !>   gradient vanishes by symmetry; the x feet are the closest, the flatter y
   !>   feet slightly farther and inside the window
   subroutine test_project_nearest_foot(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), rho(:)
      real(wp) :: rho_x, rho_y
      integer :: level, n

      call nearest_foot_fixture(fx, rho_x, rho_y, error)
      if (allocated(error)) return

      do level = 1, 3, 2
         call new_worker(fx, level, worker, error)
         if (allocated(error)) return
         call worker%project_anchor([0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, feet=pts, rho=rho)
         call check(error, .not. allocated(failed), "projection from the symmetric centre succeeds")
         if (allocated(error)) return
         call check(error, worker%n_escalated, 1, "the local solve fails")
         if (allocated(error)) return
         call check(error, n > 0, "the anchor keeps a foot")
         if (allocated(error)) return
         call check(error, abs(minval(rho(1:n)) - rho_x) < POS_THR, "the closest foot is the x pole")
         if (allocated(error)) return
      end do
   end subroutine test_project_nearest_foot

   !> The second search finds the closest foot that the nearest seed missed
   !>
   !> - one seed in the first round: the nearest crossing lies on the flatter y
   !>   sheet and converges to its farther pole
   !> - without the second search that farther foot is all the anchor keeps
   !> - the weight cut is disabled, so the raw branch sets are compared
   subroutine test_project_nearest_confirmed(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), rho(:)
      real(wp) :: rho_x, rho_y
      integer :: n

      call nearest_foot_fixture(fx, rho_x, rho_y, error)
      if (allocated(error)) return

      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      worker%nearest_count = 1
      worker%fcfg%wleb_cut = -1.0_wp
      call worker%project_anchor([0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, feet=pts, rho=rho)
      call check(error, .not. allocated(failed), "projection succeeds")
      if (allocated(error)) return
      call check(error, n > 0, "the anchor keeps a foot")
      if (allocated(error)) return
      call check(error, abs(minval(rho(1:n)) - rho_x) < POS_THR, "the second search finds the x pole")
      if (allocated(error)) return

      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      worker%nearest_count = 1
      worker%confirm_nearest = .false.
      worker%fcfg%wleb_cut = -1.0_wp
      call worker%project_anchor([0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, feet=pts, rho=rho)
      call check(error, .not. allocated(failed), "projection without the second search succeeds")
      if (allocated(error)) return
      call check(error, n, 1, "the first round alone keeps one foot")
      if (allocated(error)) return
      call check(error, abs(rho(1) - rho_y) < POS_THR, "the nearest seed alone ends on the farther y pole")
   end subroutine test_project_nearest_confirmed

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
      !> The exact refusal text
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (err%message == expected) call test_failed(error, err%message)
   end subroutine expect_error

   !> Strict mode refuses a saddle that the caller cannot handle
   subroutine test_bad_strict_saddle(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda
      logical :: kkt

      ! Setup failures return cleanly, so the test is reported as failed
      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call on_surface(fx, [-rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, is_kkt=kkt)
      call expect_error(error, failed, "Saddle point detected but strict minimum required")
   end subroutine test_bad_strict_saddle

   !> A surface flat to third order is reported as unreliable
   subroutine test_bad_flat_optimality(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda
      logical :: saddle, kkt, converged

      if (.not. sphere_setup(fx, rs)) return
      anchor = 0.0_wp
      call on_surface(fx, [rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed, saddle, kkt, converged)
      call expect_error(error, failed, "Degenerate to third order; classification unreliable")
   end subroutine test_bad_flat_optimality

   !> The centre has no gradient to retract along
   subroutine test_bad_retraction_centre(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call retract_to_surface(fx%pcfg, fx%lsf, [0.0_wp, 0.0_wp, 0.0_wp], anchor, r, lambda, failed)
      call expect_error(error, failed, "Retraction failed: degenerate LSF gradient")
   end subroutine test_bad_retraction_centre

   !> One Newton step cannot retract to the surface from 6 Bohr
   subroutine test_bad_retraction_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      fx%pcfg%retract_max_iter = 1
      call retract_to_surface(fx%pcfg, fx%lsf, [rs + 6.0_wp, 0.5_wp, 0.0_wp], anchor, r, lambda, failed)
      call expect_error(error, failed, "Retraction to surface did not converge")
   end subroutine test_bad_retraction_budget

   !> An exhausted iteration budget stops the Riemannian escape
   subroutine test_bad_escape_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r0(3), r(3), lambda0, lambda

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call on_surface(fx, [-rs, 0.0_wp, 0.0_wp], anchor, r0, lambda0, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      fx%pcfg%riemann_alpha_init = 1.0_wp
      fx%pcfg%riemann_max_iter = 1
      call riemannian_newton_escape(fx%pcfg, fx%lsf, r0, anchor, lambda0, r, lambda, failed)
      call expect_error(error, failed, "Riemannian Newton did not converge within max iterations")
   end subroutine test_bad_escape_budget

   !> At the centre there is no surface normal to define a tangent plane
   subroutine test_bad_escape_centre(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), r(3), lambda

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call riemannian_newton_escape(fx%pcfg, fx%lsf, [0.0_wp, 0.0_wp, 0.0_wp], anchor, 0.0_wp, r, lambda, failed)
      call expect_error(error, failed, "Riemannian escape: degenerate LSF gradient")
   end subroutine test_bad_escape_centre

   !> A starved Newton budget is a failed point, not a wrong one
   subroutine test_bad_newton_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_config_type) :: ncfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3)
      integer :: icol

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      ncfg = newton_config(fx)
      ncfg%max_iter = 1
      call batch%init(4, failed)
      if (allocated(failed)) return
      call batch%begin_run(1, anchor, 1)
      call batch%append([-rs - 0.05_wp, 0.02_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, fx%lsf, ncfg, scratch, failed)
      if (allocated(failed)) return
      if (batch%status(1) == point_failed_iterations) then
         call test_failed(error, "Newton iteration budget exhausted")
      end if
   end subroutine test_bad_newton_budget

   !> An exhausted octree box budget is raised as a build-aborting error
   subroutine test_bad_octree_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3)
      integer :: n

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      call new_worker(fx, 4, worker, error, max_boxes=16)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call worker%project_anchor(anchor, 1, n, failed)
      call expect_error(error, failed, "Octree branch search: box budget exhausted before the "// &
                        "search ball was certified. Either the level set supplies no exclusion "// &
                        "radius (pruning nothing), or seed_size is too small for the admissible radius.")
   end subroutine test_bad_octree_budget

   !> A nearest-foot round without seeds is refused
   subroutine test_bad_nearest_count(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      integer :: n

      if (.not. sphere_setup(fx, rs)) return
      call build_fixture(fx, carbon_pair(2.2_wp*rs), error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call new_worker(fx, 1, worker, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      worker%nearest_count = 0
      call worker%project_anchor([0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed)
      call expect_error(error, failed, "DROP projection: nearest_count must be positive")
   end subroutine test_bad_nearest_count

   !> A flat reduced Hessian with a nonzero third derivative is an inflection,
   !> refused when a strict minimum is required
   !>
   !> - neck of a carbon pair: the anchor sits on the normal at the focal
   !>   distance of the (negative) meridian curvature, so the meridian direction
   !>   is flat, and a surface of revolution has no mirror symmetry along it
   subroutine test_bad_inflection_strict(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: r(3), lambda0, p(3, 1), v(1), g(3, 1), h(6, 1), n(3), tc(3), tm(3)
      real(wp) :: km, w, lambda, anchor(3)

      call build_fixture(fx, carbon_pair(1.6_wp), error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call on_surface(fx, [0.7_wp, 2.1_wp, 0.9_wp], [0.0_wp, 0.0_wp, 0.0_wp], r, lambda0, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      p(:, 1) = r
      call fx%lsf%evaluate_batch(p, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=failed)
      if (allocated(failed)) return
      n = g(:, 1)/norm2(g(:, 1))
      ! Circumferential and meridian tangents of the pair's surface of revolution
      tc = [0.0_wp, -r(3), r(2)]/hypot(r(2), r(3))
      tm = [n(2)*tc(3) - n(3)*tc(2), n(3)*tc(1) - n(1)*tc(3), n(1)*tc(2) - n(2)*tc(1)]
      km = packed_full_contract2(h(:, 1), tm, 3)
      if (km >= 0.0_wp) return
      w = fx%pcfg%phi_alpha
      lambda = w/km
      anchor = r - (lambda/w)*g(:, 1)
      call check_constrained_optimality(fx%pcfg, fx%lsf, r, anchor, lambda, failed)
      call expect_error(error, failed, "Inflection point detected but strict minimum required")
   end subroutine test_bad_inflection_strict

   !> Retraction that lands on a critical point of the level set is refused
   !>
   !> - the midpoint of a carbon pair has a vanishing gradient by symmetry; the
   !>   separation is bisected until the midpoint lies on the zero level set
   subroutine test_bad_retraction_critical(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: h_in, h_out, h_mid, s_mid, g(3), r(3), lambda
      integer :: it

      h_in = 1.0_wp
      h_out = 6.0_wp
      if (midpoint_value(h_in) >= 0.0_wp .or. midpoint_value(h_out) <= 0.0_wp) return
      do it = 1, 200
         h_mid = 0.5_wp*(h_in + h_out)
         s_mid = midpoint_value(h_mid)
         if (abs(s_mid) < 1.0e-13_wp) exit
         if (s_mid < 0.0_wp) then
            h_in = h_mid
         else
            h_out = h_mid
         end if
      end do
      if (abs(s_mid) >= 1.0e-13_wp) return
      call build_fixture(fx, carbon_pair(h_mid), error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call retract_to_surface(fx%pcfg, fx%lsf, [0.0_wp, 0.0_wp, 0.0_wp], [3.0_wp, 0.0_wp, 0.0_wp], &
                              r, lambda, failed)
      call expect_error(error, failed, "Retraction failed: degenerate LSF gradient at converged point")

   contains

      !> Level set at the pair midpoint for half-separation `h`
      real(wp) function midpoint_value(h) result(s0)
         real(wp), intent(in) :: h
         type(fixture_type) :: fp
         type(error_type), allocatable :: e

         s0 = 0.0_wp
         call build_fixture(fp, carbon_pair(h), e)
         if (allocated(e)) then
            deallocate (e)
            return
         end if
         call lsf_value(fp, [0.0_wp, 0.0_wp, 0.0_wp], s0, g, e)
         if (allocated(e)) deallocate (e)
      end function midpoint_value
   end subroutine test_bad_retraction_critical

   !> A worker refuses projection level 5
   subroutine test_bad_worker_level(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs

      if (.not. sphere_setup(fx, rs)) return
      param = fx%param
      param%proj_level = 5
      call worker%init(param, fx%lsf, fx%mol%xyz, failed)
      call expect_error(error, failed, "DROP projection: level must be 1, 3 or 4")
   end subroutine test_bad_worker_level

   !> A worker refuses projection level 2
   subroutine test_bad_worker_level_2(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs

      if (.not. sphere_setup(fx, rs)) return
      param = fx%param
      param%proj_level = 2
      call worker%init(param, fx%lsf, fx%mol%xyz, failed)
      call expect_error(error, failed, "DROP projection: level must be 1, 3 or 4")
   end subroutine test_bad_worker_level_2

   !> A tangent block that vanishes while both switches stay on is refused
   !>
   !> - `w I - lambda Hess S = 0` on the tangent plane; the focal switch is
   !>   shifted so that it does not suppress the singular point
   subroutine test_bad_singular_jacobian(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc
      type(moist_error_type), allocatable :: failed
      real(wp) :: hess(6), jac, w_f

      call new_smoothstep_sigmoid(f_crit, 0.1_wp, 0.2_wp)
      call new_smoothstep_sigmoid(f_foc, -2.0_wp, -1.0_wp)
      !> Packed identity
      hess = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]
      call cp_jacobian_weight(1.0_wp, [0.0_wp, 0.0_wp, 1.0_wp], hess, 1.0_wp, [0.0_wp, 0.0_wp, 2.0_wp], &
                              [0.0_wp, 0.0_wp, 0.0_wp], f_crit, f_foc, 1.0e-12_wp, jac, w_f, failed)
      call expect_error(error, failed, "Tangent Jacobian matrix B is singular after switching")
   end subroutine test_bad_singular_jacobian

   !* ================================================================================= *!
   !*                                      Fixtures                                     *!
   !* ================================================================================= *!

   !> One carbon at the origin
   function single_atom() result(mol)
      !> Structure
      type(structure_type) :: mol

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
   end function single_atom

   !> Two carbons on the x axis at +-h
   !>
   !> @param[in] h Half-separation (Bohr)
   function carbon_pair(h) result(mol)
      !> Half-separation
      real(wp), intent(in) :: h
      !> Structure
      type(structure_type) :: mol

      call new(mol, [6, 6], reshape([-h, 0.0_wp, 0.0_wp, h, 0.0_wp, 0.0_wp], [3, 2]))
   end function carbon_pair

   !> Two small atoms on x and two large ones on y around the origin, and their poles
   !>
   !> - radii 2 and 4 Bohr at 5 and 7.03 Bohr: the gradient vanishes at the
   !>   origin, the x poles are closest, the y poles 0.03 Bohr farther
   !> - the pole distances come from bisection of the level set along the axes
   !>
   !> @param[out] fx    Fixture
   !> @param[out] rho_x Distance of the x poles from the origin
   !> @param[out] rho_y Distance of the y poles from the origin
   !> @param[out] error Error handle
   subroutine nearest_foot_fixture(fx, rho_x, rho_y, error)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Pole distances
      real(wp), intent(out) :: rho_x, rho_y
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      rho_x = 0.0_wp
      rho_y = 0.0_wp
      call new(mol, [6, 6, 6, 6], reshape([-5.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, 0.0_wp, &
                                           0.0_wp, -7.03_wp, 0.0_wp, 0.0_wp, 7.03_wp, 0.0_wp], [3, 4]))
      call build_fixture(fx, mol, error)
      if (allocated(error)) return
      fx%radii = [2.0_wp, 2.0_wp, 4.0_wp, 4.0_wp]
      call fx%lsf%update(fx%mol, fx%radii)
      call axis_root(fx, [1.0_wp, 0.0_wp, 0.0_wp], 5.0_wp, rho_x, error)
      if (allocated(error)) return
      call axis_root(fx, [0.0_wp, 1.0_wp, 0.0_wp], 7.03_wp, rho_y, error)
      if (allocated(error)) return
      call check(error, rho_y > rho_x + 0.01_wp, "the y poles are farther than the x poles")
   end subroutine nearest_foot_fixture

   !> First zero of the level set along a ray from the origin, by bisection
   !>
   !> @param[inout] fx    Fixture
   !> @param[in]    dir   Unit direction
   !> @param[in]    hi    Ray length to a point inside
   !> @param[out]   r     Distance of the zero
   !> @param[out]   error Bracket failure
   subroutine axis_root(fx, dir, hi, r, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Direction and bracket end
      real(wp), intent(in) :: dir(3), hi
      !> Root distance
      real(wp), intent(out) :: r
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: lo, up, s_lo, s, g(3)
      integer :: it

      lo = 0.0_wp
      up = hi
      r = 0.0_wp
      call lsf_value(fx, lo*dir, s_lo, g, error)
      if (allocated(error)) return
      call lsf_value(fx, up*dir, s, g, error)
      if (allocated(error)) return
      call check(error, s_lo*s < 0.0_wp, "the ray brackets the surface")
      if (allocated(error)) return
      do it = 1, 200
         r = 0.5_wp*(lo + up)
         call lsf_value(fx, r*dir, s, g, error)
         if (allocated(error)) return
         if (s*s_lo > 0.0_wp) then
            lo = r
         else
            up = r
         end if
         if (up - lo < 1.0e-14_wp) exit
      end do
      r = 0.5_wp*(lo + up)
   end subroutine axis_root

   !> Bind a level set and projection settings to a structure
   !>
   !> @param[out] fx      Fixture
   !> @param[in]  mol     Structure
   !> @param[out] error   Error handle
   !> @param[in]  blend_k Level set softness (default 2.5)
   subroutine build_fixture(fx, mol, error, blend_k)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Level set softness
      real(wp), intent(in), optional :: blend_k

      type(moist_error_type), allocatable :: failed
      real(wp) :: k_blend

      k_blend = 2.5_wp
      if (present(blend_k)) k_blend = blend_k
      fx%mol = mol
      call fill_legacy_radii(mol, fx%radii, error)
      if (allocated(error)) return
      call fx%param%new(nleb=50, tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=1, error=failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call fx%lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k_blend, blend_3b=1.0_wp))
      fx%lsf%screening_threshold = fx%param%screening_threshold
      call fx%lsf%update(mol, fx%radii)
      fx%pcfg%phi_alpha = fx%param%phi_alpha
   end subroutine build_fixture

   !> Newton settings of the fixture
   !>
   !> @param[in] fx Fixture
   function newton_config(fx) result(cfg)
      !> Fixture
      type(fixture_type), intent(in) :: fx
      !> Newton settings
      type(newton_config_type) :: cfg

      cfg%phi_alpha = fx%param%phi_alpha
      cfg%tol = fx%param%proj_tol
      cfg%max_iter = fx%param%proj_maxiter
   end function newton_config

   !> A projection worker at `level` on the fixture
   !>
   !> @param[in]  fx        Fixture
   !> @param[in]  level     Projection level
   !> @param[out] worker    Configured worker
   !> @param[out] error     Error handle
   !> @param[in]  max_boxes Octree box budget override
   subroutine new_worker(fx, level, worker, error, max_boxes)
      !> Fixture
      type(fixture_type), intent(in) :: fx
      !> Projection level
      integer, intent(in) :: level
      !> Configured worker
      type(projection_worker_type), intent(out) :: worker
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Octree box budget
      integer, intent(in), optional :: max_boxes

      type(moist_cavity_drop_parameters_type) :: param
      type(moist_error_type), allocatable :: failed

      param = fx%param
      param%proj_level = level
      if (present(max_boxes)) param%octree_max_boxes = max_boxes
      call worker%init(param, fx%lsf, fx%mol%xyz, failed)
      if (allocated(failed)) call test_failed(error, "worker setup failed: "//failed%message)
   end subroutine new_worker

   !> Radius of the isolated-atom surface, found by projecting a far anchor onto it
   !>
   !> @param[in]  fx    Single-atom fixture at the origin
   !> @param[out] rs    Radius of the zero level set
   !> @param[out] error Error handle
   subroutine surface_radius(fx, rs, error)
      !> Single-atom fixture
      type(fixture_type), intent(in) :: fx
      !> Radius
      real(wp), intent(out) :: rs
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :)
      integer :: n

      rs = 0.0_wp
      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      call worker%project_anchor([8.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, feet=pts)
      if (allocated(failed)) then
         call test_failed(error, "surface radius projection failed: "//failed%message)
         return
      end if
      call check(error, n, 1, "single projection for the radius")
      if (allocated(error)) return
      rs = norm2(pts(:, 1))
      call check(error, abs(pts(1, 1) - rs) < 1.0e-12_wp, "radius probe lies on the axis")
   end subroutine surface_radius

   !> Single-atom fixture and its radius for the error paths; `.false.` on a setup failure
   !>
   !> @param[out] fx Fixture
   !> @param[out] rs Surface radius
   logical function sphere_setup(fx, rs) result(ok)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Surface radius
      real(wp), intent(out) :: rs

      type(error_type), allocatable :: error

      rs = 0.0_wp
      call build_fixture(fx, single_atom(), error)
      ok = .not. allocated(error)
      if (.not. ok) return
      call surface_radius(fx, rs, error)
      ok = .not. allocated(error)
   end function sphere_setup

   !> Level set value and gradient at one point
   !>
   !> @param[inout] fx    Fixture
   !> @param[in]    x     Point
   !> @param[out]   s     Value
   !> @param[out]   g     Gradient
   !> @param[out]   error Error handle
   subroutine lsf_value(fx, x, s, g, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Point
      real(wp), intent(in) :: x(3)
      !> Value and gradient
      real(wp), intent(out) :: s, g(3)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: failed
      real(wp) :: p(3, 1), v(1), gg(3, 1)

      p(:, 1) = x
      call fx%lsf%evaluate_batch(p, 1, v, lsf1_r=gg, error=failed)
      s = v(1)
      g = gg(:, 1)
      if (allocated(failed)) call test_failed(error, failed%message)
   end subroutine lsf_value

   !> A point on the surface near `guess`, with the multiplier that makes it stationary
   !>
   !> @param[inout] fx     Fixture
   !> @param[in]    guess  Approximate position
   !> @param[in]    anchor Anchor of the projection
   !> @param[out]   r      Surface point
   !> @param[out]   lambda Stationarity multiplier of the anchor objective
   !> @param[out]   error  Error handle
   subroutine on_surface(fx, guess, anchor, r, lambda, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
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

      call retract_to_surface(fx%pcfg, fx%lsf, guess, anchor, r, lambda, failed)
      if (allocated(failed)) call test_failed(error, "retraction of the fixture point failed: "//failed%message)

   end subroutine on_surface

end module test_cavity_drop_projector
