!> Tests for the per-point DROP projection pieces: candidate deduplication, the
!> horizon gate, the focal gate at umbilics, optimality classification, saddle escape (sphere, necks,
!> multi-sphere clusters, the threshold band and inflection points), retraction,
!> and the single-anchor projection entry at every level, including the octree
!> retry of a failed local solve, the restart of an anchor whose seeds all fail,
!> and an anchor without surface in reach
module test_cavity_drop_projector
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, drop_solve_parameters_type, &
                                           drop_polish_parameters_type
   use moist_math_packed, only: packed_full_contract2
   use moist_math_linalg, only: cross_product, setup_tangent_frame
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_converged, point_accepted, point_failed_iterations
   use moist_cavity_drop_projection_newton, only: newton_scratch_type, newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_batch, &
                                                 check_constrained_optimality, riemannian_newton_escape, &
                                                 retract_to_surface
   use moist_cavity_drop_projection_finalize, only: dedupe_candidates, cp_jacobian_weight
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_projection_worker, only: projection_worker_type
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: fill_legacy_radii
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, to_string
   implicit none(type, external)
   private

   public :: collect_cavity_drop_projector

   !> Solver tolerance of the projection under test
   real(wp), parameter :: PROJ_TOL = 1.0e-14_wp
   !> Position tolerance of the exact sphere and mirror geometries (Bohr)
   real(wp), parameter :: POS_THR = 1.0e-10_wp
   !> Half-separation of the neck fixtures in carbon radii; the two atoms overlap
   !> into a waist, and near it the axis cylinder of the waist radius lies inside the surface
   real(wp), parameter :: NECK_HALF = 0.8_wp
   !> Position tolerance on the neck points; the waist of unequal atoms is found by a
   !> golden-section search that stalls near the square root of machine precision
   real(wp), parameter :: NECK_THR = 1.0e-6_wp
   !> Surface step of the curvature probes (Bohr)
   real(wp), parameter :: PROBE_STEP = 0.05_wp
   !> Smallest objective change a curvature probe must show, against roundoff and a
   !> first-order residual
   real(wp), parameter :: PROBE_RISE = 1.0e-5_wp

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
      !> Run context
      type(moist_context_type) :: ctx
   end type fixture_type

contains

   !> Collect all DROP projector tests
   subroutine collect_cavity_drop_projector(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dedupe_candidates_and_slack", test_dedupe_candidates), &
                  new_unittest("horizon_gate", test_horizon_gate), &
                  new_unittest("focal_gate_at_umbilic", test_focal_gate_at_umbilic), &
                  new_unittest("optimality_classifies_min_max_and_non_kkt", test_optimality_classes), &
                  new_unittest("optimality_degenerate_gradient", test_optimality_degenerate_gradient), &
                  new_unittest("retraction_lands_on_surface", test_retraction), &
                  new_unittest("polish_escapes_saddle", test_refine_saddle), &
                  new_unittest("polish_escapes_neck_ring_saddle", test_escape_neck_ring), &
                  new_unittest("polish_escapes_neck_meridian_saddle", test_escape_neck_meridian), &
                  new_unittest("polish_escapes_exact_neck_saddle", test_escape_exact_neck), &
                  new_unittest("escape_pocket_axis_maximum", test_escape_pocket_maximum), &
                  new_unittest("escape_pocket_pass_saddle", test_escape_pocket_pass), &
                  new_unittest("escape_cavity_saddles_in_batch", test_escape_cavity_batch), &
                  new_unittest("escape_bent_chain_axis_saddle", test_escape_bent_chain), &
                  new_unittest("saddle_band_check_and_escape_agree", test_saddle_band_agreement), &
                  new_unittest("saddle_threshold_scales_with_hessian_norm", test_saddle_threshold_scales), &
                  new_unittest("polish_keeps_point_flat_to_third_order", test_polish_flat_third_order), &
                  new_unittest("polish_never_accepts_inflection_in_place", test_polish_inflection_not_kept), &
                  new_unittest("polish_escapes_inflection", test_polish_inflection_escapes), &
                  new_unittest("project_anchor_single_atom_levels", test_project_single_levels), &
                  new_unittest("project_anchor_midpoint_branches", test_project_midpoint), &
                  new_unittest("project_anchor_failed_local_octree", test_project_failed_octree), &
                  new_unittest("project_anchor_failed_seeds_restart", test_project_restart), &
                  new_unittest("project_anchor_no_surface_in_reach", test_project_unseeded), &
                  new_unittest("project_anchor_failed_local_nearest_point", test_project_nearest_point), &
                  new_unittest("project_anchor_nearest_point_confirmed", test_project_nearest_confirmed), &
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

      !> Anchor at the origin: pdist^2 = 1.0, 1.0 (duplicate), 1.4, 1.6, 4.0
      cand(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
      cand(:, 2) = [1.0_wp, 1.0e-12_wp, 0.0_wp]
      cand(:, 3) = [0.0_wp, sqrt(1.4_wp), 0.0_wp]
      cand(:, 4) = [0.0_wp, 0.0_wp, sqrt(1.6_wp)]
      cand(:, 5) = [-2.0_wp, 0.0_wp, 0.0_wp]

      call dedupe_candidates(cand, anchor, 1.0e-10_wp, 0.5_wp, keep)
      call check(error, all(keep .eqv. [.true., .false., .true., .false., .false.]), &
                 "duplicate merged, pdist^2 slack of 0.5 keeps 1.4 and rejects 1.6 and 4.0")
      if (allocated(error)) return

      !> A wider slack admits the branch that was out of range, and a wide
      !> separation cut merges candidates that are distinct at the default
      call dedupe_candidates(cand, anchor, 1.0e-10_wp, 0.7_wp, keep)
      call check(error, all(keep .eqv. [.true., .false., .true., .true., .false.]), &
                 "slack override admits pdist^2 = 1.6")
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
   !*                                 Closest-point gate                                *!
   !* ================================================================================= *!

   !> The horizon switch closes the gate of a point grazing the anchor sphere
   !>
   !> - flat level set, so the tangent block is the identity and the area scaling is the
   !>   cosine between surface and sphere normal
   !> - the gate is zero below the band, one above it, strictly between inside and growing
   !>   with the cosine; the band is on the signed cosine, so a normal pointing into the
   !>   anchor sphere closes the gate
   !> - the band comes from the parameters; `w_0hor_on = 0` switches it off, which leaves a
   !>   grazing point fully open, with its small area scaling
   subroutine test_horizon_gate(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: C_OFF = 0.05_wp, C_ON = 0.15_wp
      real(wp), parameter :: cosines(5) = [0.01_wp, 0.08_wp, 0.12_wp, 0.3_wp, 1.0_wp]
      type(moist_math_smoothstep_sigmoid_type) :: f_hor
      real(wp) :: g(3), jac, w_f(5), w_flip
      integer :: i

      call horizon_switch(C_OFF, C_ON, f_hor, error)
      if (allocated(error)) return

      do i = 1, size(cosines)
         g = [sqrt(1.0_wp - cosines(i)**2), 0.0_wp, cosines(i)]
         call flat_gate(g, f_hor, jac, w_f(i), error)
         if (allocated(error)) return
         if (w_f(i) > 0.0_wp) then
            call check(error, jac, cosines(i), "an open gate carries the cosine as area scaling", thr=1.0e-12_wp)
            if (allocated(error)) return
         end if
         call flat_gate(-g, f_hor, jac, w_flip, error)
         if (allocated(error)) return
         call check(error, w_flip == 0.0_wp, "an inward normal closes the gate")
         if (allocated(error)) return
      end do

      call check(error, w_f(1) == 0.0_wp, "the gate is closed below the horizon band")
      if (allocated(error)) return
      call check(error, w_f(2) > 0.0_wp .and. w_f(2) < w_f(3) .and. w_f(3) < 1.0_wp, &
                 "the gate rises strictly inside the horizon band")
      if (allocated(error)) return
      call check(error, w_f(4) == 1.0_wp .and. w_f(5) == 1.0_wp, "the gate is fully open above the horizon band")
      if (allocated(error)) return

      !> Switched off through the parameters
      call horizon_switch(C_OFF, 0.0_wp, f_hor, error)
      if (allocated(error)) return
      g = [sqrt(1.0_wp - cosines(1)**2), 0.0_wp, cosines(1)]
      call flat_gate(g, f_hor, jac, w_f(1), error)
      if (allocated(error)) return
      call check(error, w_f(1) == 1.0_wp, "with the switch off a grazing point is fully open")
      if (allocated(error)) return
      call check(error, jac, cosines(1), "and keeps its small area scaling", thr=1.0e-12_wp)

   end subroutine test_horizon_gate

   !> The focal margin of an umbilic point carries no square-root noise
   !>
   !> - isotropic level-set Hessian `h I`, so the tangent block is `beta I` with
   !>   `beta = w - lambda h` for every normal; normals off the axes give a generic
   !>   tangent frame, whose rounding leaves `B` umbilic only to machine precision
   !> - `tr^2/4 - det` then cancels to noise of order `eps beta^2`, and its square root
   !>   moved the margin by up to `sqrt(eps) beta`, about 1e-8 (the 2026-10-05 bug, found
   !>   at the symmetry axes of H2 and PbH4); the margins lie on the focal ramp, where
   !>   that noise reaches the gate
   !> - critical and horizon switches open, so the gate is the focal switch at the exact
   !>   margin
   subroutine test_focal_gate_at_umbilic(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Objective weight and multiplier
      real(wp), parameter :: W = 1.0_wp, LAMBDA = 0.25_wp
      !> Isotropic Hessian scales; exact focal margins 0.28125, 0.296875 and 0.3125
      real(wp), parameter :: HSCALE(3) = [2.875_wp, 2.8125_wp, 2.75_wp]
      !> Focal ramp, steep around every margin (gate 0.18 to 0.73)
      real(wp), parameter :: FOC_OFF = 0.1_wp, FOC_ON = 0.5_wp
      !> Gate deviation allowed
      !>
      !> - measured 1.7e-15 with the sum of squares
      !> - 7.3e-8 with the cancelling discriminant, above the bound for 26 of the 72 cases
      real(wp), parameter :: UMBILIC_THR = 1.0e-13_wp
      !> Surface normals tried per margin
      integer, parameter :: NDIR = 24
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_hor
      type(moist_error_type), allocatable :: failed
      real(wp) :: n(3), hess(6), jac, w_f, expected, beta, theta, phi
      integer :: i, k

      call new_smoothstep_sigmoid(f_crit, 0.1_wp, 0.2_wp)
      call new_smoothstep_sigmoid(f_foc, FOC_OFF, FOC_ON)
      call horizon_switch(0.05_wp, 0.15_wp, f_hor, error)
      if (allocated(error)) return

      do k = 1, size(HSCALE)
         hess = [HSCALE(k), 0.0_wp, HSCALE(k), 0.0_wp, 0.0_wp, HSCALE(k)]
         beta = W - LAMBDA*HSCALE(k)
         expected = f_foc%f0(beta/W)
         call check(error, expected > 1.0e-2_wp .and. expected < 0.99_wp, &
                    "fixture margin "//to_string(beta)//" is off the focal ramp")
         if (allocated(error)) return
         do i = 1, NDIR
            theta = 0.2_wp + 0.11_wp*real(i, wp)
            phi = 0.7_wp + 0.53_wp*real(i, wp)
            n = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
            ! Anchor along the normal, so the horizon switch is open
            call cp_jacobian_weight(W, n, hess, LAMBDA, 2.0_wp*n, [0.0_wp, 0.0_wp, 0.0_wp], &
                                    f_crit, f_foc, f_hor, jac, w_f, failed)
            if (allocated(failed)) then
               call test_failed(error, "umbilic gate failed: "//failed%message)
               return
            end if
            call check(error, abs(w_f - expected) <= UMBILIC_THR, &
                       "focal gate at an umbilic, margin "//to_string(beta)//", normal "//to_string(i)// &
                       ": "//to_string(w_f)//" instead of "//to_string(expected))
            if (allocated(error)) return
         end do
      end do

   end subroutine test_focal_gate_at_umbilic

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
      type(drop_polish_parameters_type) :: loose
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
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "near pole classified without error")
      if (allocated(error)) return
      call check(error, kkt .and. converged .and. .not. saddle, "near pole is a constrained minimum")
      if (allocated(error)) return

      !> Far pole: stationary, but the displacement is maximal along the surface
      call on_surface(fx, [-rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "far pole flagged through the saddle output, no error")
      if (allocated(error)) return
      call check(error, kkt .and. saddle .and. .not. converged, "far pole is a saddle-type KKT point")
      if (allocated(error)) return

      !> Without a way to hand the saddle back, strict mode must refuse it
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, is_kkt=kkt)
      call check(error, allocated(failed), "strict mode refuses a saddle that the caller cannot handle")
      if (allocated(error)) return
      loose = fx%param%polish
      loose%strict_minimum_required = .false.
      call check_constrained_optimality(fx%ctx, loose, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, is_kkt=kkt)
      call check(error, .not. allocated(failed), "non-strict mode accepts it")
      if (allocated(error)) return

      !> Off the surface
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        [rs + 0.2_wp, 0.0_wp, 0.0_wp], anchor, lambda, &
                                        failed, saddle, kkt, converged)
      call check(error, .not. kkt .and. .not. converged, "a point 0.2 Bohr off the surface is no KKT point")
      if (allocated(error)) return

      !> On the surface, but not the closest point: the stationarity condition fails
      call on_surface(fx, [0.0_wp, rs, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, saddle, kkt, converged)
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

      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        mid, anchor, 0.0_wp, failed, is_kkt=kkt)
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
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        mid, anchor, 0.0_wp, failed, is_kkt=kkt)
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
      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                              [rs + 0.4_wp, 0.0_wp, 0.0_wp], anchor, r, lambda, failed)
      call check(error, .not. allocated(failed), "retraction succeeds")
      if (allocated(error)) return
      call lsf_value(fx, r, s, g, error)
      if (allocated(error)) return
      call check(error, abs(s) < 1.0e-12_wp, "retracted point is on the surface")
      if (allocated(error)) return
      call check(error, norm2(r - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "the radial offset retracts radially")
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, is_kkt=kkt, &
                                        is_converged=converged)
      call check(error, kkt .and. converged, "the retraction multiplier is the stationarity multiplier")
      if (allocated(error)) return

      !> A skew start lands on the sphere as well
      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                              [-1.3_wp, 2.1_wp, 0.7_wp], anchor, r, lambda, failed)
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
      type(drop_solve_parameters_type) :: ncfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), seed(3)
      integer :: icol

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      ncfg = fx%param%solve

      !> A seed on the axis behind the far pole: off the axis the merit line
      !> search already walks Newton away from the maximum
      seed = [-rs - 0.05_wp, 0.0_wp, 0.0_wp]
      call batch%init(4, failed)
      call batch%begin_run(1, anchor, 1)
      call batch%append(seed, 0.0_wp, icol)
      call newton_solve_batch(batch, fx%lsf, ncfg, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "Newton succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges")
      if (allocated(error)) return

      !> Without the check the maximum is what Newton found
      call check(error, norm2(batch%x(:, 1) - [-rs, 0.0_wp, 0.0_wp]) < POS_THR, &
                 "Newton alone stops on the maximum")
      if (allocated(error)) return

      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, ncfg, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_accepted, "escaped point verified")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "polish ends on the near pole")
      if (allocated(error)) return
      call check(error, abs(batch%pdist(1) - 3.0_wp) < POS_THR, &
                 "reported displacement is the anchor height above the sphere")
      if (allocated(error)) return
      call check(error, norm2(batch%normal(:, 1) - [1.0_wp, 0.0_wp, 0.0_wp]) < POS_THR, &
                 "unit normal points along the axis")

   end subroutine test_refine_saddle

   !* ================================================================================= *!
   !*                           Saddle escape on necks and clusters                     *!
   !* ================================================================================= *!

   !> Newton converges onto the far neck point of a dumbbell, a mixed-signature saddle
   !> of the anchor distance; the polish walks round the neck ring to the near point
   !>
   !> - anchor between the axis and the neck wall, in the plane of the waist; the waist
   !>   is the narrowest cross-section, so between the atom centres the surface lies
   !>   outside the cylinder of the waist radius, and the end caps beyond them are far
   !>   from the anchor: the near neck point is the global minimum, at the waist radius
   !>   minus the anchor height
   !> - far neck point: farthest point of the ring (negative curvature), waist profile
   !>   curving away (positive); the two signs are probed on the surface, not read from
   !>   the Hessian the classification uses
   !> - equal atoms (waist on the mirror plane) and unequal atoms (waist off the plane)
   subroutine test_escape_neck_ring(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: scales(2, 2) = reshape([1.0_wp, 1.0_wp, 0.8_wp, 1.2_wp], [2, 2])
      real(wp), parameter :: HEIGHT = 0.9_wp
      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: xm, rm, anchor(3), far(3), near(3)
      integer :: ic
      logical :: saddle, kkt, converged

      do ic = 1, 2
         call neck_fixture(fx, scales(:, ic), error)
         if (allocated(error)) return
         if (ic == 1) then
            xm = 0.0_wp
            call ring_radius(fx, xm, rm, error)
         else
            call neck_location(fx, xm, rm, error)
         end if
         if (allocated(error)) return
         anchor = [xm, HEIGHT, 0.0_wp]
         far = [xm, -rm, 0.0_wp]
         near = [xm, rm, 0.0_wp]

         !> A seed in the plane of the symmetry: off the plane the merit line search
         !> already walks Newton away from the saddle
         call solve_columns(fx, anchor, reshape([xm, -rm - 0.2_wp, 0.0_wp], [3, 1]), batch, scratch, error)
         if (allocated(error)) return
         call check(error, batch%status(1), point_converged, "Newton converges")
         if (allocated(error)) return
         call check(error, norm2(batch%x(:, 1) - far) < NECK_THR, "Newton alone stops on the far neck point")
         if (allocated(error)) return
         call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                           batch%x(:, 1), anchor, batch%lambda(1), failed, saddle, kkt, converged)
         call check(error, .not. allocated(failed), "far neck point classified without error")
         if (allocated(error)) return
         call check(error, kkt .and. saddle .and. .not. converged, "far neck point is a saddle-type KKT point")
         if (allocated(error)) return
         !> Falls round the ring, rises along the waist profile
         call expect_mixed(fx, anchor, batch%x(:, 1), [0.0_wp, 0.0_wp, 1.0_wp], error)
         if (allocated(error)) return

         call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
         call check(error, .not. allocated(failed), "polish succeeds")
         if (allocated(error)) return
         call check(error, batch%status(1), point_accepted, "escaped point verified")
         if (allocated(error)) return
         call check(error, norm2(batch%x(:, 1) - near) < NECK_THR, "polish ends on the near neck point")
         if (allocated(error)) return
         call check(error, abs(batch%pdist(1) - (rm - HEIGHT)) < 1.0e-9_wp, &
                    "displacement is the waist radius minus the anchor height")
         if (allocated(error)) return
      end do

   end subroutine test_escape_neck_ring

   !> A mirror-plane anchor above the neck: Newton stops on the neck point, a saddle with a
   !> maximum along the meridian, and the polish ends on one of the two mirror points
   !>
   !> - anchor beyond the focal point of the waist profile, so the neck point is farthest
   !>   along the meridian and nearest along the ring
   !> - the escape is a descent along the meridian, where both sides are mirror images;
   !>   only the pair is asserted: off the plane, lower than the saddle, a strict minimum
   !>   with a mirror twin of equal multiplier
   !> - same input, same point
   subroutine test_escape_neck_meridian(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: HEIGHT = 2.0_wp
      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rm, anchor(3), neck(3), point(3), twin(3), seed(3, 1)
      logical :: saddle, kkt, converged

      call neck_fixture(fx, [1.0_wp, 1.0_wp], error)
      if (allocated(error)) return
      call ring_radius(fx, 0.0_wp, rm, error)
      if (allocated(error)) return
      anchor = [0.0_wp, rm + HEIGHT, 0.0_wp]
      neck = [0.0_wp, rm, 0.0_wp]
      seed(:, 1) = [0.0_wp, rm + 0.2_wp, 0.0_wp]

      call solve_columns(fx, anchor, seed, batch, scratch, error)
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - neck) < POS_THR, "Newton alone stops on the neck point")
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        batch%x(:, 1), anchor, batch%lambda(1), failed, saddle, kkt, converged)
      call check(error, .not. allocated(failed), "neck point classified without error")
      if (allocated(error)) return
      call check(error, kkt .and. saddle .and. .not. converged, "neck point is a saddle-type KKT point")
      if (allocated(error)) return
      !> Falls along the meridian, rises round the ring
      call expect_mixed(fx, anchor, batch%x(:, 1), [1.0_wp, 0.0_wp, 0.0_wp], error)
      if (allocated(error)) return

      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_accepted, "escaped point verified")
      if (allocated(error)) return
      point = batch%x(:, 1)
      call check(error, abs(point(1)) > 0.1_wp, "the point left the mirror plane")
      if (allocated(error)) return
      call check(error, batch%pdist(1) < HEIGHT - 0.05_wp, "the point is closer than the neck point")
      if (allocated(error)) return
      call check_surface_minimum(fx, anchor, point, batch%normal(:, 1), error)
      if (allocated(error)) return

      !> The other mirror point is as good: same multiplier, strict minimum
      twin = [-point(1), point(2), point(3)]
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        twin, anchor, batch%lambda(1), failed, is_kkt=kkt, is_converged=converged)
      call check(error, .not. allocated(failed) .and. kkt .and. converged, "the mirror twin is a strict minimum too")
      if (allocated(error)) return

      !> Same input, same point
      call solve_columns(fx, anchor, seed, batch, scratch, error)
      if (allocated(error)) return
      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "second polish succeeds")
      if (allocated(error)) return
      call check(error, maxval(abs(batch%x(:, 1) - point)) == 0.0_wp, "the escape is deterministic")

   end subroutine test_escape_neck_meridian

   !> A start exactly on the saddle: the escape walk alone and the polish both leave it
   !>
   !> - the neck points sit on the axes of the dumbbell, so the saddles are known without
   !>   a solve, and the retraction multiplier is the exact stationarity multiplier there
   !> - ring saddle: the polish ends on the near neck point; meridian saddle: on a point
   !>   off the mirror plane, a strict minimum
   subroutine test_escape_exact_neck(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rm, anchor(3), saddle_pt(3), lambda, r_out(3), lambda_out, s_out, g_out(3)
      logical :: saddle, kkt, converged
      integer :: icase

      call neck_fixture(fx, [1.0_wp, 1.0_wp], error)
      if (allocated(error)) return
      call ring_radius(fx, 0.0_wp, rm, error)
      if (allocated(error)) return

      do icase = 1, 2
         if (icase == 1) then
            !> Far neck point under an anchor inside the waist
            anchor = [0.0_wp, 0.9_wp, 0.0_wp]
            call on_surface(fx, [0.0_wp, -rm, 0.0_wp], anchor, saddle_pt, lambda, error)
         else
            !> Neck point under an anchor above it
            anchor = [0.0_wp, rm + 2.0_wp, 0.0_wp]
            call on_surface(fx, [0.0_wp, rm, 0.0_wp], anchor, saddle_pt, lambda, error)
         end if
         if (allocated(error)) return
         call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                           saddle_pt, anchor, lambda, failed, saddle, kkt, converged)
         call check(error, .not. allocated(failed), "the start is classified without error")
         if (allocated(error)) return
         call check(error, kkt .and. saddle .and. .not. converged, "the start is a saddle-type KKT point")
         if (allocated(error)) return

         !> The walk alone: the descent from the saddle, on the surface, to convergence
         call riemannian_newton_escape(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                       saddle_pt, anchor, lambda, r_out, lambda_out, failed)
         call check(error, .not. allocated(failed), "the walk converges within its budget")
         if (allocated(error)) return
         call lsf_value(fx, r_out, s_out, g_out, error)
         if (allocated(error)) return
         call check(error, abs(s_out) < 1.0e-12_wp, "the walk stays on the surface")
         if (allocated(error)) return
         call check(error, sum((r_out - anchor)**2) < sum((saddle_pt - anchor)**2) - 0.5_wp, &
                    "the walk descends from the saddle")
         if (allocated(error)) return
         if (icase == 1) then
            call check(error, norm2(r_out - [0.0_wp, rm, 0.0_wp]) < 0.1_wp, "the ring walk heads for the near neck point")
            if (allocated(error)) return
         end if

         !> The polish from a seed on the saddle
         call solve_columns(fx, anchor, reshape(saddle_pt, [3, 1]), batch, scratch, error, lambdas=[lambda])
         if (allocated(error)) return
         call check(error, batch%status(1), point_converged, "Newton accepts the exact saddle")
         if (allocated(error)) return
         call check(error, norm2(batch%x(:, 1) - saddle_pt) < POS_THR, "and stays on it")
         if (allocated(error)) return
         call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
         call check(error, .not. allocated(failed), "polish succeeds")
         if (allocated(error)) return
         call check(error, batch%status(1), point_accepted, "escaped point verified")
         if (allocated(error)) return
         call check(error, norm2(batch%x(:, 1) - saddle_pt) > 0.1_wp, "the polish left the saddle")
         if (allocated(error)) return
         call check_surface_minimum(fx, anchor, batch%x(:, 1), batch%normal(:, 1), error)
         if (allocated(error)) return
         if (icase == 1) then
            call check(error, norm2(batch%x(:, 1) - [0.0_wp, rm, 0.0_wp]) < NECK_THR, "polish ends on the near neck point")
            if (allocated(error)) return
         end if
      end do

   end subroutine test_escape_exact_neck

   !> Escape from the isotropic maximum on the C3 axis of an equilateral triangle
   !>
   !> - three equal spheres, circumradius one sphere radius; the anchor sits on the axis above
   !>   the pocket between them
   !> - D3h: the axis point has equal reduced curvatures, a maximum of the distance; the three
   !>   points over the atoms are the symmetry-equivalent minima
   !> - Newton from the axis stays on it (mirror planes), the polish must walk off in a direction
   !>   the symmetry does not pick, and end on one of the three minima
   subroutine test_escape_pocket_maximum(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), mins(3, 3), start(3, 1)
      logical :: saddle, kkt
      integer :: k

      call cluster_fixture(fx, triangle_xyz(1.0_wp), rs, error)
      if (allocated(error)) return
      anchor = [0.0_wp, 0.0_wp, rs]
      call pocket_minima(fx, anchor, rs, mins, error)
      if (allocated(error)) return
      call check(error, abs(norm2(mins(:, 1) - anchor) - norm2(mins(:, 2) - anchor)) < POS_THR .and. &
                 abs(norm2(mins(:, 1) - anchor) - norm2(mins(:, 3) - anchor)) < POS_THR, &
                 "the three minima are equally far from the anchor")
      if (allocated(error)) return

      !> Start on the axis: Newton stops on the maximum
      start(:, 1) = [0.0_wp, 0.0_wp, 0.7_wp*rs]
      call solve_columns(fx, anchor, start, batch, scratch, error)
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges on the axis")
      if (allocated(error)) return
      call check(error, norm2(batch%x(1:2, 1)) < POS_THR, "Newton stays on the symmetry axis")
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, batch%x(:, 1), &
                                        anchor, batch%lambda(1), failed, is_saddle=saddle, is_kkt=kkt)
      call check(error, kkt .and. saddle, "the axis point is a saddle-type KKT point")
      if (allocated(error)) return

      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_accepted, "the escaped point is verified")
      if (allocated(error)) return
      call check(error, any([(norm2(batch%x(:, 1) - mins(:, k)) < POS_THR, k=1, 3)]), &
                 "the escape ends on one of the three equivalent minima")

   end subroutine test_escape_pocket_maximum

   !> Escape from a mixed saddle between two equivalent minima, on and off the saddle
   !>
   !> - same triangle and anchor; over each edge the pass between the two neighbouring minima is
   !>   a saddle with one negative and one positive reduced curvature
   !> - the mirror plane through the opposite atom holds it: Newton from there stays in the plane
   !>   and stops on the pass, the descent leaves the plane, to a flanking minimum
   !> - the second start sits exactly on the saddle Newton found
   subroutine test_escape_pocket_pass(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), mins(3, 3), start(3, 1), pass(3), lambda_pass
      logical :: saddle, kkt
      integer :: k, stage

      call cluster_fixture(fx, triangle_xyz(1.0_wp), rs, error)
      if (allocated(error)) return
      anchor = [0.0_wp, 0.0_wp, rs]
      call pocket_minima(fx, anchor, rs, mins, error)
      if (allocated(error)) return

      !> Atom 1 sits on +y, so the pass over the opposite edge is in the plane x = 0, at y < 0
      start(:, 1) = [0.0_wp, -0.1_wp*rs, 0.6_wp*rs]
      call solve_columns(fx, anchor, start, batch, scratch, error)
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges")
      if (allocated(error)) return
      pass = batch%x(:, 1)
      lambda_pass = batch%lambda(1)
      call check(error, abs(pass(1)) < POS_THR .and. pass(2) < 0.0_wp, "Newton stays in the mirror plane, over the edge")
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, pass, &
                                        anchor, lambda_pass, failed, is_saddle=saddle, is_kkt=kkt)
      call check(error, kkt .and. saddle, "the pass is a saddle-type KKT point")
      if (allocated(error)) return
      call expect_mixed(fx, anchor, pass, [1.0_wp, 0.0_wp, 0.0_wp], error)
      if (allocated(error)) return

      do stage = 1, 2
         if (stage == 2) then
            !> Start exactly on the saddle
            call solve_columns(fx, anchor, reshape(pass, [3, 1]), batch, scratch, error, lambdas=[lambda_pass])
            if (allocated(error)) return
            call check(error, norm2(batch%x(:, 1) - pass) < POS_THR, "Newton leaves the saddle where it is")
            if (allocated(error)) return
         end if
         call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
         call check(error, .not. allocated(failed), "polish succeeds")
         if (allocated(error)) return
         call check(error, batch%status(1), point_accepted, "the escaped point is verified")
         if (allocated(error)) return
         !> The negative curvature runs out of the plane, down to the minima next to the pass
         call check(error, any([(norm2(batch%x(:, 1) - mins(:, k)) < POS_THR, k=2, 3)]), &
                    "the escape ends on a minimum flanking the pass")
         if (allocated(error)) return
         call check(error, norm2(batch%x(:, 1) - mins(:, 1)) > 1.0_wp, "and not on the minimum on the far side")
         if (allocated(error)) return
      end do

   end subroutine test_escape_pocket_pass

   !> Several saddles in one batch, one batch of a cavity with four minima
   !>
   !> - four spheres on a tetrahedron, anchor at its centre inside the cavity the spheres leave
   !>   open; Td: four equivalent minima towards the spheres, saddles over the six edges
   !> - columns: two inner saddles on different C2 axes, the outer saddle on one of them (a long
   !>   walk through the gap between two spheres), and one start that is already a minimum
   !> - every saddle column must end on a minimum, the minimum column must come back untouched
   subroutine test_escape_cavity_batch(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), starts(3, 4), landed(3, 4), x_min(3)
      logical :: saddle, kkt
      integer :: j

      call cluster_fixture(fx, tetrahedron_xyz(1.28_wp), rs, error)
      if (allocated(error)) return
      anchor = 0.0_wp
      starts(:, 1) = [0.5_wp*rs, 0.0_wp, 0.0_wp]
      starts(:, 2) = [0.0_wp, 0.0_wp, -0.5_wp*rs]
      starts(:, 3) = [0.9_wp*rs, 0.0_wp, 0.0_wp]
      starts(:, 4) = [0.15_wp, 0.15_wp, 0.15_wp]*rs
      call solve_columns(fx, anchor, starts, batch, scratch, error)
      if (allocated(error)) return
      do j = 1, 4
         call check(error, batch%status(j), point_converged, "every column converges")
         if (allocated(error)) return
      end do
      landed = batch%x(:, 1:4)

      !> Newton alone: three saddles on the C2 axes and a minimum
      do j = 1, 4
         call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, batch%x(:, j), &
                                           anchor, batch%lambda(j), failed, is_saddle=saddle, is_kkt=kkt)
         call check(error, kkt .and. (saddle .eqv. j < 4), "columns 1-3 are saddles, column 4 a minimum")
         if (allocated(error)) return
      end do
      call check(error, norm2(landed(2:3, 1)) < POS_THR .and. norm2(landed(1:2, 2)) < POS_THR .and. &
                 norm2(landed(2:3, 3)) < POS_THR, "the saddles sit on the C2 axes")
      if (allocated(error)) return
      call check(error, landed(1, 3) > landed(1, 1) + 1.0_wp, "column 3 is the outer saddle, farther along the axis")
      if (allocated(error)) return
      !> Towards the neighbouring spheres the surface falls away, across the edge it rises
      call expect_mixed(fx, anchor, landed(:, 1), [0.0_wp, 1.0_wp, 1.0_wp], error)
      if (allocated(error)) return
      call expect_mixed(fx, anchor, landed(:, 2), [1.0_wp, -1.0_wp, 0.0_wp], error)
      if (allocated(error)) return

      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return

      do j = 1, 4
         call check(error, batch%status(j), point_accepted, "every column is verified")
         if (allocated(error)) return
         !> The four minima are the vertex directions: equal coordinates up to sign, positive product
         x_min = batch%x(:, j)
         call check(error, abs(abs(x_min(1)) - abs(x_min(2))) < POS_THR .and. &
                    abs(abs(x_min(1)) - abs(x_min(3))) < POS_THR .and. x_min(1)*x_min(2)*x_min(3) > 0.0_wp, &
                    "the column ends on one of the four vertex minima")
         if (allocated(error)) return
         call check(error, abs(norm2(x_min - anchor) - norm2(batch%x(:, 4) - anchor)) < POS_THR, &
                    "all minima are equally far from the centre")
         if (allocated(error)) return
      end do
      call check(error, all(batch%x(:, 4) == landed(:, 4)), "the minimum column is left untouched")
      if (allocated(error)) return

      !> The inner saddles drain into the two minima beside their edge
      call check(error, batch%x(1, 1) > 0.0_wp, "the saddle on +x ends on a minimum of an edge sphere")
      if (allocated(error)) return
      call check(error, batch%x(3, 2) < 0.0_wp, "the saddle on -z ends on a minimum of an edge sphere")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 3) - landed(:, 3)) > 2.0_wp, "the outer saddle is left over a long walk")

   end subroutine test_escape_cavity_batch

   !> Escape from the saddle on the axis of a bent chain, between two mirror minima
   !>
   !> - spheres at the apex and at both ends of a 120 degree V; the anchor lies on the axis in
   !>   the cleft, farther out than the apex surface
   !> - the axis point is a saddle, the two points on the arm spheres are mirror images
   subroutine test_escape_bent_chain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3), start(3, 1), twin(3, 1), mirror(3)
      logical :: saddle, kkt

      call cluster_fixture(fx, bent_chain_xyz(1.2_wp, 120.0_wp), rs, error)
      if (allocated(error)) return
      anchor = [0.0_wp, 0.0_wp, -1.7_wp*rs]
      twin(:, 1) = [0.3_wp*rs, 0.0_wp, -1.2_wp*rs]
      call solve_columns(fx, anchor, twin, batch, scratch, error)
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges on the arm")
      if (allocated(error)) return
      twin(:, 1) = batch%x(:, 1)
      mirror = [-twin(1, 1), twin(2, 1), twin(3, 1)]

      start(:, 1) = [0.0_wp, 0.0_wp, -1.1_wp*rs]
      call solve_columns(fx, anchor, start, batch, scratch, error)
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges on the axis")
      if (allocated(error)) return
      call check(error, norm2(batch%x(1:2, 1)) < POS_THR, "Newton stays on the symmetry axis")
      if (allocated(error)) return
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, batch%x(:, 1), &
                                        anchor, batch%lambda(1), failed, is_saddle=saddle, is_kkt=kkt)
      call check(error, kkt .and. saddle, "the axis point is a saddle-type KKT point")
      if (allocated(error)) return
      call expect_mixed(fx, anchor, batch%x(:, 1), [1.0_wp, 0.0_wp, 0.0_wp], error)
      if (allocated(error)) return

      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_accepted, "the escaped point is verified")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - twin(:, 1)) < POS_THR .or. norm2(batch%x(:, 1) - mirror) < POS_THR, &
                 "the escape ends on one of the two mirror minima")

   end subroutine test_escape_bent_chain

   !* ================================================================================= *!
   !*                        Saddle threshold and degenerate points                     *!
   !* ================================================================================= *!

   !> Check and escape walk call the same points saddles, across the threshold band
   !>
   !> - sphere, anchor on the axis at height `h`, point at the far pole: the reduced
   !>   eigenvalue is `-phi_alpha h/R` twice, so `h` tunes it continuously through zero
   !>   (`h < 0` makes the point the near pole, a minimum)
   !> - `saddle_eig_tol = 1e-2` puts the threshold at about 5e-3 (`phi_alpha = 0.5`), so the
   !>   band is +-5e-3: `h = 0.2` is clearly below it, `h = +-0.02` inside it (the
   !>   eigenvalue is resolved to 1e-16), `h = -0.2` above
   !> - a saddle is walked off, a point inside the band or above it is left where it is,
   !>   exactly; `h = 0.02` is the case a walk threshold other than the check's would split
   !> - the polish of the clear saddle ends on the near pole
   subroutine test_saddle_band_agreement(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      real(wp) :: rs, walked
      logical :: saddle
      integer :: ih
      !> Anchor heights above the centre and the expected classification
      real(wp), parameter :: heights(4) = [0.2_wp, 0.02_wp, -0.02_wp, -0.2_wp]
      logical, parameter :: is_saddle(4) = [.true., .false., .false., .false.]

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      fx%param%polish%saddle_eig_tol = 1.0e-2_wp

      do ih = 1, size(heights)
         call far_pole_regime(fx, rs, heights(ih), saddle, walked, error)
         if (allocated(error)) return
         call check(error, saddle .eqv. is_saddle(ih), "classification of the far pole follows the band")
         if (allocated(error)) return
         if (is_saddle(ih)) then
            ! Both poles are 2 R apart; leaving the saddle means leaving its neighbourhood
            call check(error, walked > 1.0_wp, "the walk leaves a saddle")
         else
            call check(error, walked == 0.0_wp, "the walk does not move a point the check keeps")
         end if
         if (allocated(error)) return
      end do

      !> The polish of the clear saddle ends on the near pole
      call polish_far_pole(fx, rs, heights(1), error)

   end subroutine test_saddle_band_agreement

   !> The threshold scales with `||H_L||`: the same agreement at three Hessian scales
   !>
   !> - default `saddle_eig_tol = 1e-8`, so the band is `1e-8 ||H_L||_F`: about
   !>   `5e-10` for `phi_alpha = 0.05`, `1e-9` for `phi_alpha = 0.1` (an absolute 1e-8
   !>   would be above both) and `1e-7` for `phi_alpha = 10` (an absolute 1e-8 would be
   !>   below it); `||H_L||` is `phi_alpha` on the sphere to within 1e-3
   !> - eigenvalues `-3e-8`, `-6e-10`, `-3e-11` at `phi_alpha = 0.05`; `-3e-8`, `-3e-9`,
   !>   `-3e-11` at `phi_alpha = 0.1`; `-3e-6`, `-3e-8`, `-3e-10` at `phi_alpha = 10`;
   !>   `h = mu R / phi_alpha` from the sphere relation
   !> - a saddle just below the band gains `mu s^2/2` over an escape step `s`: about
   !>   `7.5e-11` for `-6e-10` and the default step, so the walk must accept a decrease
   !>   on the scale of the band, not an absolute one
   subroutine test_saddle_threshold_scales(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      real(wp) :: rs, walked
      logical :: saddle
      integer :: ip, ih
      real(wp), parameter :: phis(3) = [0.05_wp, 0.1_wp, 10.0_wp]
      !> Eigenvalues per scale: clearly below the band, just below or near its upper side, deep inside
      real(wp), parameter :: mus(3, 3) = reshape([-3.0e-8_wp, -6.0e-10_wp, -3.0e-11_wp, &
                                                  -3.0e-8_wp, -3.0e-9_wp, -3.0e-11_wp, &
                                                  -3.0e-6_wp, -3.0e-8_wp, -3.0e-10_wp], [3, 3])
      !> Whether the eigenvalue lies below the band; the second one does not at 10
      logical, parameter :: is_saddle(3, 3) = reshape([.true., .true., .false., &
                                                       .true., .true., .false., &
                                                       .true., .false., .false.], [3, 3])

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return

      do ip = 1, size(phis)
         fx%param%phi_alpha = phis(ip)
         do ih = 1, 3
            call far_pole_regime(fx, rs, -mus(ih, ip)*rs/phis(ip), saddle, walked, error)
            if (allocated(error)) return
            call check(error, saddle .eqv. is_saddle(ih, ip), "classification follows the relative band")
            if (allocated(error)) return
            call check(error, (walked > 1.0_wp) .eqv. is_saddle(ih, ip), "the walk agrees with the check")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_saddle_threshold_scales

   !> A point flat to third order is kept by the polish, where Newton left it
   !>
   !> - anchor at the centre of the sphere: every point is stationary and flat
   !>   (the check reports it as unreliable, see `bad_flat_optimality`), the polish
   !>   accepts it and does not move it
   subroutine test_polish_flat_third_order(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      real(wp) :: rs, x(3)
      integer :: status

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return

      call polish_column(fx, [0.0_wp, 0.0_wp, 0.0_wp], [-rs - 0.05_wp, 0.0_wp, 0.0_wp], 0.0_wp, status, x, error)
      if (allocated(error)) return
      call check(error, status, point_accepted, "the unclassifiable point is kept")
      if (allocated(error)) return
      call check(error, norm2(x - [-rs, 0.0_wp, 0.0_wp]) < POS_THR, "and stays where it is")

   end subroutine test_polish_flat_third_order

   !> An inflection point is never accepted where it stands
   !>
   !> - flat meridian direction at the neck of a carbon pair, third derivative nonzero
   !>   (`eps = 0`), and the same neck slightly past the focal distance (`eps = 1e-3`,
   !>   eigenvalue `-5e-4`, about 3e4 times the band): a saddle for the check
   !> - the polish may escape it or give it up, but must not report it a minimum while it
   !>   stays: a saddle that is accepted in place is a wrong point with full weight
   subroutine test_polish_inflection_not_kept(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: r(3), anchor(3), lambda, x(3)
      integer :: ie, status
      logical :: saddle, kkt
      real(wp), parameter :: epss(2) = [0.0_wp, 1.0e-3_wp]

      do ie = 1, size(epss)
         call inflection_point(fx, epss(ie), r, anchor, lambda, error)
         if (allocated(error)) return
         call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                           r, anchor, lambda, failed, saddle, kkt)
         call check(error, kkt .and. saddle .and. .not. allocated(failed), "the check calls the neck a saddle")
         if (allocated(error)) return
         call polish_column(fx, anchor, r, lambda, status, x, error)
         if (allocated(error)) return
         call check(error, .not. (status == point_accepted .and. norm2(x - r) < POS_THR), &
                    "a saddle is not accepted in place")
         if (allocated(error)) return
      end do

   end subroutine test_polish_inflection_not_kept

   !> An inflection point and a near-inflection saddle are escaped to a minimum
   !>
   !> - same necks as `polish_never_accepts_inflection_in_place`; the first-order change
   !>   along the flat or weakly negative direction is zero, only the third-order term
   !>   tells the downhill side from the uphill one
   subroutine test_polish_inflection_escapes(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      real(wp) :: r(3), anchor(3), lambda, x(3)
      integer :: ie, status
      real(wp), parameter :: epss(2) = [0.0_wp, 1.0e-3_wp]

      do ie = 1, size(epss)
         call inflection_point(fx, epss(ie), r, anchor, lambda, error)
         if (allocated(error)) return
         call polish_column(fx, anchor, r, lambda, status, x, error)
         if (allocated(error)) return
         call check(error, status, point_accepted, "the escape reaches a minimum")
         if (allocated(error)) return
         call check(error, norm2(x - r) > 0.1_wp, "away from the inflection point")
         if (allocated(error)) return
         call check(error, sum((x - anchor)**2) < sum((r - anchor)**2), "closer to the anchor than the inflection point")
         if (allocated(error)) return
      end do

   end subroutine test_polish_inflection_escapes

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
      real(wp), allocatable :: pts(:, :), pdist(:)
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
         call worker%project_anchor(fx%ctx, anchor, 1, n, failed, points=pts, pdist=pdist)
         call check(error, .not. allocated(failed), "projection succeeds")
         if (allocated(error)) return
         call check(error, n, 1, "a convex surface has one branch")
         if (allocated(error)) return
         call check(error, abs(norm2(pts(:, 1)) - rs) < POS_THR, "projected point is on the sphere")
         if (allocated(error)) return
         call check(error, norm2(pts(:, 1)/rs - anchor/norm2(anchor)) < POS_THR, &
                    "projection is along the line to the centre")
         if (allocated(error)) return
         call check(error, abs(pdist(1) - (norm2(anchor) - rs)) < POS_THR, &
                    "displacement is the height above the sphere")
         if (allocated(error)) return
         if (level == 1) then
            ! Converged local point and no global search: no seeds at all
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
      real(wp), allocatable :: pts(:, :), pdist(:)
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
      ! The mirror point on the other atom is a backward projection (c = -1),
      ! which the signed horizon switch closes under defaults
      fx%param%w_0hor_on = 0.0_wp
      anchor = 0.0_wp

      do il = 1, size(levels)
         level = levels(il)
         call new_worker(fx, level, worker, error)
         if (allocated(error)) return
         call worker%project_anchor(fx%ctx, anchor, 1, n, failed, points=pts, pdist=pdist)
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
         call check(error, abs(pdist(1) - pdist(2)) < POS_THR, "and equally far from the anchor")
         if (allocated(error)) return
      end do

   end subroutine test_project_midpoint

   !> A failed local solve at level 1 is retried from the coarse octree of
   !> level 3: same seeds, same points
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
      call worker%project_anchor(fx%ctx, anchor, 1, n_ref, failed, points=ref, n_seeds=seeds_ref)
      call check(error, .not. allocated(failed), "level-3 reference projection succeeds")
      if (allocated(error)) return

      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      call worker%project_anchor(fx%ctx, anchor, 1, n, failed, points=pts, n_seeds=seeds, lost=lost)
      call check(error, .not. allocated(failed), "projection succeeds")
      if (allocated(error)) return
      call check(error, worker%n_escalated, 1, "the failed local solve is escalated")
      if (allocated(error)) return
      call check(error, seeds, seeds_ref, "the retry proposes the level-3 octree seeds")
      if (allocated(error)) return
      call check(error, .not. lost, "the anchor keeps a point")
      if (allocated(error)) return
      call check(error, n, n_ref, "the level-3 branch count")
      if (allocated(error)) return
      call check(error, maxval(abs(pts(:, 1:n) - ref(:, 1:n_ref))) < POS_THR, "the level-3 points")

   end subroutine test_project_failed_octree

   !> An anchor whose seeds all fail is solved again with the full system
   !>
   !> - anchor on the axis of an argon-pair neck: the surface is a ring around it, and a
   !>   long cheap stage (50 steps) carries every seed onto the axis
   !> - `h = 4.2` gives a neck radius of about 1.04 Bohr; the same holds from about 3.6
   !>   (radius 2.4) up to the pinch-off just above 4.33, so the value is not tuned
   !> - every ring point is equally far from the anchor, which is the focal point of the
   !>   ring and sees its normals tangent to the anchor sphere, so the points carry no
   !>   weight; the weight cut is disabled to see them
   !> - negative control: without the restart (seeds solved once, as a worker
   !>   whose restart keeps the cheap steps does) the anchor is lost with seeds
   subroutine test_project_restart(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      type(structure_type) :: mol
      real(wp), parameter :: h = 4.2_wp
      real(wp), allocatable :: pdist(:)
      real(wp) :: anchor(3), rm
      integer :: n, seeds
      logical :: lost, unseeded

      call new(mol, [18, 18], reshape([0.0_wp, 0.0_wp, -h, 0.0_wp, 0.0_wp, h], [3, 2]))
      call build_fixture(fx, mol, error, blend_k=5.5_wp)
      if (allocated(error)) return
      fx%param%branch%reach_on = 1.6_wp
      fx%param%branch%reach_off = 2.4_wp
      anchor = 0.0_wp
      ! The waist lies in the plane z = 0, which holds the ray along +y from the origin
      call ring_radius(fx, 0.0_wp, rm, error)
      if (allocated(error)) return

      fx%param%solve%max_cheap_steps = 50
      call new_worker(fx, 3, worker, error)
      if (allocated(error)) return
      worker%solve_restart = worker%solve_seeded
      call worker%project_anchor(fx%ctx, anchor, 1, n, failed, n_seeds=seeds, lost=lost, unseeded=unseeded)
      call check(error, .not. allocated(failed), "control: projection succeeds")
      if (allocated(error)) return
      call check(error, seeds > 0 .and. lost .and. .not. unseeded, "control: every seed of the anchor fails")
      if (allocated(error)) return

      call new_worker(fx, 3, worker, error)
      if (allocated(error)) return
      worker%finalize%wleb_cut = -1.0_wp
      call worker%project_anchor(fx%ctx, anchor, 1, n, failed, pdist=pdist, lost=lost, unseeded=unseeded)
      call check(error, .not. allocated(failed), "projection succeeds")
      if (allocated(error)) return
      call check(error, worker%n_restarted, 1, "the anchor is solved again")
      if (allocated(error)) return
      call check(error, .not. lost .and. .not. unseeded .and. n > 0, "and keeps an accepted point")
      if (allocated(error)) return
      call check(error, maxval(abs(pdist(1:n) - rm)) < POS_THR, "every point lies on the neck ring")

   end subroutine test_project_restart

   !> An anchor with no surface within the reach is lost without a seed
   !>
   !> - midpoint of a carbon pair, where the vanishing gradient stalls Newton;
   !>   far apart, no sign change lies inside the reach ball
   !> - negative control: the close pair has surface in reach and keeps points
   subroutine test_project_unseeded(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3)
      integer :: n, seeds
      logical :: lost, unseeded

      call build_fixture(fx, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(fx, rs, error)
      if (allocated(error)) return
      anchor = 0.0_wp

      call build_fixture(fx, carbon_pair(2.2_wp*rs), error)
      if (allocated(error)) return
      call new_worker(fx, 3, worker, error)
      if (allocated(error)) return
      call worker%project_anchor(fx%ctx, anchor, 1, n, failed, lost=lost, unseeded=unseeded)
      call check(error, .not. allocated(failed), "close pair: projection succeeds")
      if (allocated(error)) return
      call check(error, n > 0 .and. .not. lost .and. .not. unseeded, "close pair: the anchor keeps points")
      if (allocated(error)) return

      call build_fixture(fx, carbon_pair(rs + 1.5_wp*worker%reach_off), error)
      if (allocated(error)) return
      call new_worker(fx, 3, worker, error)
      if (allocated(error)) return
      call worker%project_anchor(fx%ctx, anchor, 1, n, failed, n_seeds=seeds, lost=lost, unseeded=unseeded)
      if (allocated(failed)) then
         call test_failed(error, "far pair: projection failed: "//failed%message)
         return
      end if
      call check(error, seeds, 0, "far pair: the reach ball holds no seed")
      if (allocated(error)) return
      call check(error, n == 0 .and. lost, "far pair: the anchor is lost")
      if (allocated(error)) return
      call check(error, unseeded, "far pair: and reported as without surface in reach")

   end subroutine test_project_unseeded

   !> A failed local solve gets the point closest to the anchor
   !>
   !> - centre of two small atoms on x and two large ones on y, where the
   !>   gradient vanishes by symmetry; the x points are the closest, the flatter y
   !>   points slightly farther and inside the window
   subroutine test_project_nearest_point(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), pdist(:)
      real(wp) :: pdist_x, pdist_y
      integer :: level, n

      call nearest_point_fixture(fx, pdist_x, pdist_y, error)
      if (allocated(error)) return

      do level = 1, 3, 2
         call new_worker(fx, level, worker, error)
         if (allocated(error)) return
         call worker%project_anchor(fx%ctx, [0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, points=pts, pdist=pdist)
         call check(error, .not. allocated(failed), "projection from the symmetric centre succeeds")
         if (allocated(error)) return
         call check(error, worker%n_escalated, 1, "the local solve fails")
         if (allocated(error)) return
         call check(error, n > 0, "the anchor keeps a point")
         if (allocated(error)) return
         call check(error, abs(minval(pdist(1:n)) - pdist_x) < POS_THR, "the closest point is the x pole")
         if (allocated(error)) return
      end do
   end subroutine test_project_nearest_point

   !> The second search finds the closest point that the nearest seed missed
   !>
   !> - one seed in the first round: the nearest crossing lies on the flatter y
   !>   sheet and converges to its farther pole
   !> - without the second search that farther point is all the anchor keeps
   !> - the weight cut is disabled, so the raw branch sets are compared
   subroutine test_project_nearest_confirmed(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_worker_type) :: worker
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: pts(:, :), pdist(:)
      real(wp) :: pdist_x, pdist_y
      integer :: n

      call nearest_point_fixture(fx, pdist_x, pdist_y, error)
      if (allocated(error)) return

      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      worker%nearest_count = 1
      worker%finalize%wleb_cut = -1.0_wp
      call worker%project_anchor(fx%ctx, [0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, points=pts, pdist=pdist)
      call check(error, .not. allocated(failed), "projection succeeds")
      if (allocated(error)) return
      call check(error, n > 0, "the anchor keeps a point")
      if (allocated(error)) return
      call check(error, abs(minval(pdist(1:n)) - pdist_x) < POS_THR, "the second search finds the x pole")
      if (allocated(error)) return

      call new_worker(fx, 1, worker, error)
      if (allocated(error)) return
      worker%nearest_count = 1
      worker%confirm_nearest = .false.
      worker%finalize%wleb_cut = -1.0_wp
      call worker%project_anchor(fx%ctx, [0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, points=pts, pdist=pdist)
      call check(error, .not. allocated(failed), "projection without the second search succeeds")
      if (allocated(error)) return
      call check(error, n, 1, "the first round alone keeps one point")
      if (allocated(error)) return
      call check(error, abs(pdist(1) - pdist_y) < POS_THR, "the nearest seed alone ends on the farther y pole")
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
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, is_kkt=kkt)
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
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, saddle, kkt, converged)
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
      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                              [0.0_wp, 0.0_wp, 0.0_wp], anchor, r, lambda, failed)
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
      fx%param%polish%retract_max_iter = 1
      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                              [rs + 6.0_wp, 0.5_wp, 0.0_wp], anchor, r, lambda, failed)
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
      fx%param%polish%escape_alpha_init = 1.0_wp
      fx%param%polish%escape_max_iter = 1
      call riemannian_newton_escape(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                    r0, anchor, lambda0, r, lambda, failed)
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
      call riemannian_newton_escape(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                    [0.0_wp, 0.0_wp, 0.0_wp], anchor, 0.0_wp, r, lambda, failed)
      call expect_error(error, failed, "Riemannian escape: degenerate LSF gradient")
   end subroutine test_bad_escape_centre

   !> A starved Newton budget is a failed point, not a wrong one
   subroutine test_bad_newton_budget(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: ncfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3)
      integer :: icol

      if (.not. sphere_setup(fx, rs)) return
      anchor = [rs + 3.0_wp, 0.0_wp, 0.0_wp]
      ncfg = fx%param%solve
      ncfg%max_iter = 1
      call batch%init(4, failed)
      if (allocated(failed)) return
      call batch%begin_run(1, anchor, 1)
      call batch%append([-rs - 0.05_wp, 0.02_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, fx%lsf, ncfg, fx%param%phi_alpha, scratch, failed)
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
      call worker%project_anchor(fx%ctx, anchor, 1, n, failed)
      call expect_error(error, failed, "Octree branch search: box budget exhausted before the "// &
                        "search ball was certified. Either the level set supplies no exclusion "// &
                        "radius (pruning nothing), or seed_size is too small for the admissible radius.")
   end subroutine test_bad_octree_budget

   !> A nearest-point round without seeds is refused
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
      call worker%project_anchor(fx%ctx, [0.0_wp, 0.0_wp, 0.0_wp], 1, n, failed)
      call expect_error(error, failed, "DROP projection: nearest_count must be positive")
   end subroutine test_bad_nearest_count

   !> A flat reduced Hessian with a nonzero third derivative is an inflection,
   !> refused when a strict minimum is required
   !>
   !> - neck of a carbon pair: the anchor sits on the normal at the focal
   !>   distance of the (negative) meridian curvature, so the meridian direction
   !>   is flat, and a surface of revolution has no mirror symmetry along it
   !>   (the fixture of `polish_escapes_inflection` with `eps = 0`)
   subroutine test_bad_inflection_strict(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(fixture_type) :: fx
      type(moist_error_type), allocatable :: failed
      real(wp) :: r(3), lambda, anchor(3)

      call inflection_point(fx, 0.0_wp, r, anchor, lambda, error)
      if (allocated(error)) then
         deallocate (error)
         return
      end if
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed)
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
      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                              [0.0_wp, 0.0_wp, 0.0_wp], [3.0_wp, 0.0_wp, 0.0_wp], &
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

      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_hor
      type(moist_error_type), allocatable :: failed
      real(wp) :: hess(6), jac, w_f

      call new_smoothstep_sigmoid(f_crit, 0.1_wp, 0.2_wp)
      call new_smoothstep_sigmoid(f_foc, -2.0_wp, -1.0_wp)
      call new_smoothstep_sigmoid(f_hor, -2.0_wp, -1.0_wp)
      !> Packed identity
      hess = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]
      call cp_jacobian_weight(1.0_wp, [0.0_wp, 0.0_wp, 1.0_wp], hess, 1.0_wp, [0.0_wp, 0.0_wp, 2.0_wp], &
                              [0.0_wp, 0.0_wp, 0.0_wp], f_crit, f_foc, f_hor, jac, w_f, failed)
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
   !> @param[out] fx      Fixture
   !> @param[out] pdist_x Distance of the x poles from the origin
   !> @param[out] pdist_y Distance of the y poles from the origin
   !> @param[out] error   Error handle
   subroutine nearest_point_fixture(fx, pdist_x, pdist_y, error)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Pole distances
      real(wp), intent(out) :: pdist_x, pdist_y
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      pdist_x = 0.0_wp
      pdist_y = 0.0_wp
      call new(mol, [6, 6, 6, 6], reshape([-5.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.0_wp, 0.0_wp, &
                                           0.0_wp, -7.03_wp, 0.0_wp, 0.0_wp, 7.03_wp, 0.0_wp], [3, 4]))
      call build_fixture(fx, mol, error)
      if (allocated(error)) return
      fx%radii = [2.0_wp, 2.0_wp, 4.0_wp, 4.0_wp]
      call fx%lsf%update(fx%mol, fx%radii)
      call axis_root(fx, [1.0_wp, 0.0_wp, 0.0_wp], 5.0_wp, pdist_x, error)
      if (allocated(error)) return
      call axis_root(fx, [0.0_wp, 1.0_wp, 0.0_wp], 7.03_wp, pdist_y, error)
      if (allocated(error)) return
      call check(error, pdist_y > pdist_x + 0.01_wp, "the y poles are farther than the x poles")
   end subroutine nearest_point_fixture

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
      call fx%param%new(nleb=50, tolerance=PROJ_TOL, proj_level=1, error=failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      ! Anchors of these fixtures sit anywhere, not on an atomic sphere: a reach beyond all of them
      fx%param%branch%reach_on = 10.0_wp
      fx%param%branch%reach_off = 12.0_wp
      call fx%param%compute_derived(failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call fx%lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k_blend, blend_3b=1.0_wp))
      fx%lsf%screening_threshold = fx%param%screening_threshold
      call fx%lsf%update(mol, fx%radii)
      call new_context(fx%ctx, verbosity=0, nthreads=1)
   end subroutine build_fixture

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
      if (present(max_boxes)) param%seed%max_boxes = max_boxes
      ! Tests set user-facing fields on the fixture; resolve what derives from them
      call param%compute_derived(failed)
      if (allocated(failed)) then
         call test_failed(error, "parameter setup failed: "//failed%message)
         return
      end if
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
      call worker%project_anchor(fx%ctx, [8.0_wp, 0.0_wp, 0.0_wp], 1, n, failed, points=pts)
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

      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, guess, anchor, r, lambda, failed)
      if (allocated(failed)) call test_failed(error, "retraction of the fixture point failed: "//failed%message)

   end subroutine on_surface

   !> Horizon switch of a parameter set with the given edges, as the projection builds it
   !>
   !> @param[in]  w_off Horizon off-edge
   !> @param[in]  w_on  Horizon on-edge; zero switches it off
   !> @param[out] f_hor Horizon switch
   !> @param[out] error Invalid edges
   subroutine horizon_switch(w_off, w_on, f_hor, error)
      !> Edges
      real(wp), intent(in) :: w_off, w_on
      !> Horizon switch
      type(moist_math_smoothstep_sigmoid_type), intent(out) :: f_hor
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_parameters_type) :: param
      type(moist_error_type), allocatable :: failed

      call param%new(w_0hor_off=w_off, w_0hor_on=w_on, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "parameter setup failed: "//failed%message)
         return
      end if
      call new_smoothstep_sigmoid(f_hor, param%hor_off, param%hor_on)
   end subroutine horizon_switch

   !> Area scaling and gate of a point on a flat level set, with the other switches open
   !>
   !> - unit objective weight and multiplier, anchor 2 Bohr above its owner on z, so the
   !>   sphere normal is +z
   !>
   !> @param[in]  g     Level set gradient (unit)
   !> @param[in]  f_hor Horizon switch
   !> @param[out] jac   Area scaling
   !> @param[out] w_f   Gate
   !> @param[out] error Singular tangent block
   subroutine flat_gate(g, f_hor, jac, w_f, error)
      !> Gradient
      real(wp), intent(in) :: g(3)
      !> Horizon switch
      type(moist_math_smoothstep_sigmoid_type), intent(in) :: f_hor
      !> Area scaling and gate
      real(wp), intent(out) :: jac, w_f
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc
      type(moist_error_type), allocatable :: failed
      real(wp), parameter :: flat_hess(6) = 0.0_wp

      call new_smoothstep_sigmoid(f_crit, 0.1_wp, 0.2_wp)
      call new_smoothstep_sigmoid(f_foc, -2.0_wp, -1.0_wp)
      call cp_jacobian_weight(1.0_wp, g, flat_hess, 1.0_wp, [0.0_wp, 0.0_wp, 2.0_wp], &
                              [0.0_wp, 0.0_wp, 0.0_wp], f_crit, f_foc, f_hor, jac, w_f, failed)
      if (allocated(failed)) call test_failed(error, failed%message)
   end subroutine flat_gate

   !> Carbons at `rs` times the positions `xyz`, `rs` the surface radius of one carbon
   !>
   !> @param[out] fx    Fixture
   !> @param[in]  xyz   Positions in units of the carbon surface radius (3, n)
   !> @param[out] rs    Surface radius of one carbon
   !> @param[out] error Error handle
   !> @param[in]  radii Atom radii in units of the carbon radius (default all 1)
   subroutine cluster_fixture(fx, xyz, rs, error, radii)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Positions in units of the radius
      real(wp), intent(in) :: xyz(:, :)
      !> Radius
      real(wp), intent(out) :: rs
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Radii in units of the carbon radius
      real(wp), intent(in), optional :: radii(:)

      type(fixture_type) :: single
      type(structure_type) :: mol
      integer, allocatable :: z(:)

      rs = 0.0_wp
      call build_fixture(single, single_atom(), error)
      if (allocated(error)) return
      call surface_radius(single, rs, error)
      if (allocated(error)) return
      allocate (z(size(xyz, 2)))
      z = 6
      call new(mol, z, rs*xyz)
      call build_fixture(fx, mol, error)
      if (allocated(error) .or. .not. present(radii)) return
      fx%radii = radii*fx%radii
      call fx%lsf%update(fx%mol, fx%radii)
   end subroutine cluster_fixture

   !> Two overlapping carbons on the x axis at +-`NECK_HALF` surface radii, a waist between them
   !>
   !> @param[out] fx    Fixture
   !> @param[in]  scale Radius of each atom in carbon radii
   !> @param[out] error Error handle
   subroutine neck_fixture(fx, scale, error)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Radius factors
      real(wp), intent(in) :: scale(2)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: rs

      call cluster_fixture(fx, reshape([-NECK_HALF, 0.0_wp, 0.0_wp, NECK_HALF, 0.0_wp, 0.0_wp], [3, 2]), &
                           rs, error, radii=scale)
   end subroutine neck_fixture

   !> Equilateral triangle in the xy plane around the z axis, first vertex on +y
   !>
   !> @param[in] rc Circumradius
   function triangle_xyz(rc) result(xyz)
      !> Circumradius
      real(wp), intent(in) :: rc
      !> Positions
      real(wp) :: xyz(3, 3)

      real(wp), parameter :: c30 = 0.86602540378443865_wp

      xyz(:, 1) = [0.0_wp, rc, 0.0_wp]
      xyz(:, 2) = [-c30*rc, -0.5_wp*rc, 0.0_wp]
      xyz(:, 3) = [c30*rc, -0.5_wp*rc, 0.0_wp]
   end function triangle_xyz

   !> Regular tetrahedron around the origin, vertices along (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
   !>
   !> @param[in] rc Circumradius
   function tetrahedron_xyz(rc) result(xyz)
      !> Circumradius
      real(wp), intent(in) :: rc
      !> Positions
      real(wp) :: xyz(3, 4)

      real(wp) :: t

      t = rc/sqrt(3.0_wp)
      xyz(:, 1) = t*[1.0_wp, 1.0_wp, 1.0_wp]
      xyz(:, 2) = t*[1.0_wp, -1.0_wp, -1.0_wp]
      xyz(:, 3) = t*[-1.0_wp, 1.0_wp, -1.0_wp]
      xyz(:, 4) = t*[-1.0_wp, -1.0_wp, 1.0_wp]
   end function tetrahedron_xyz

   !> Bent three-atom chain in the xz plane: apex at the origin, arms opening towards -z
   !>
   !> @param[in] arm   Arm length
   !> @param[in] angle Apex angle in degrees
   function bent_chain_xyz(arm, angle) result(xyz)
      !> Arm length
      real(wp), intent(in) :: arm
      !> Apex angle (degrees)
      real(wp), intent(in) :: angle
      !> Positions
      real(wp) :: xyz(3, 3)

      real(wp) :: half

      half = 0.5_wp*angle*pi/180.0_wp
      xyz(:, 1) = [-arm*sin(half), 0.0_wp, -arm*cos(half)]
      xyz(:, 2) = 0.0_wp
      xyz(:, 3) = [arm*sin(half), 0.0_wp, -arm*cos(half)]
   end function bent_chain_xyz

   !> Radius of the surface at axial position `x`, by bisection along +y from the axis
   !>
   !> @param[inout] fx    Fixture, axis point inside
   !> @param[in]    x     Axial position
   !> @param[out]   r     Distance of the surface from the axis
   !> @param[out]   error Axis point outside
   subroutine ring_radius(fx, x, r, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Axial position
      real(wp), intent(in) :: x
      !> Radius
      real(wp), intent(out) :: r
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: lo, up, s, g(3)
      integer :: it

      lo = 0.0_wp
      up = 20.0_wp
      r = 0.0_wp
      call lsf_value(fx, [x, lo, 0.0_wp], s, g, error)
      if (allocated(error)) return
      call check(error, s < 0.0_wp, "the axis point is inside")
      if (allocated(error)) return
      do it = 1, 60
         r = 0.5_wp*(lo + up)
         call lsf_value(fx, [x, r, 0.0_wp], s, g, error)
         if (allocated(error)) return
         if (s < 0.0_wp) then
            lo = r
         else
            up = r
         end if
      end do
      r = 0.5_wp*(lo + up)
   end subroutine ring_radius

   !> Axial position and radius of the waist, by golden-section search of the ring radius
   !>
   !> @param[inout] fx    Fixture
   !> @param[out]   xm    Axial position of the narrowest cross-section
   !> @param[out]   rm    Its radius
   !> @param[out]   error Error handle
   subroutine neck_location(fx, xm, rm, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Waist position and radius
      real(wp), intent(out) :: xm, rm
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: GOLD = 0.6180339887498949_wp
      real(wp) :: a, b, c, d, fc, fd
      integer :: it

      !> Between the atom centres the radius is unimodal: full at the centres, pinched between
      a = fx%mol%xyz(1, 1)
      b = fx%mol%xyz(1, 2)
      c = b - GOLD*(b - a)
      d = a + GOLD*(b - a)
      call ring_radius(fx, c, fc, error)
      if (allocated(error)) return
      call ring_radius(fx, d, fd, error)
      if (allocated(error)) return
      do it = 1, 60
         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b - GOLD*(b - a)
            call ring_radius(fx, c, fc, error)
         else
            a = c
            c = d
            fc = fd
            d = a + GOLD*(b - a)
            call ring_radius(fx, d, fd, error)
         end if
         if (allocated(error)) return
      end do
      xm = 0.5_wp*(a + b)
      call ring_radius(fx, xm, rm, error)
   end subroutine neck_location

   !> The three minima over the atoms of the triangle pocket, by Newton alone
   !>
   !> @param[inout] fx     Triangle fixture
   !> @param[in]    anchor Anchor on the C3 axis
   !> @param[in]    rs     Surface radius of one carbon
   !> @param[out]   mins   Minimum over each atom (3, 3)
   !> @param[out]   error  A start that does not converge to a minimum
   subroutine pocket_minima(fx, anchor, rs, mins, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Anchor and radius
      real(wp), intent(in) :: anchor(3), rs
      !> Minima
      real(wp), intent(out) :: mins(3, 3)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: over_atom(3, 1), xyz(3)
      logical :: saddle, kkt
      integer :: k

      mins = 0.0_wp
      do k = 1, 3
         xyz = fx%mol%xyz(:, k)
         over_atom(:, 1) = [0.25_wp*xyz(1:2)/norm2(xyz(1:2))*rs, 0.7_wp*rs]
         call solve_columns(fx, anchor, over_atom, batch, scratch, error)
         if (allocated(error)) return
         call check(error, batch%status(1), point_converged, "Newton converges over an atom")
         if (allocated(error)) return
         call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, batch%x(:, 1), &
                                           anchor, batch%lambda(1), failed, is_saddle=saddle, is_kkt=kkt)
         call check(error, kkt .and. .not. saddle, "the point over an atom is a minimum")
         if (allocated(error)) return
         mins(:, k) = batch%x(:, 1)
      end do
   end subroutine pocket_minima

   !> Newton from each start, as the columns of one run
   !>
   !> @param[inout] fx      Fixture
   !> @param[in]    anchor  Anchor of the run
   !> @param[in]    starts  Start points (3, n)
   !> @param[inout] batch   Batch, rebuilt
   !> @param[inout] scratch Newton scratch
   !> @param[out]   error   Error handle
   !> @param[in]    lambdas Start multipliers (default 0)
   subroutine solve_columns(fx, anchor, starts, batch, scratch, error, lambdas)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Anchor
      real(wp), intent(in) :: anchor(3)
      !> Start points
      real(wp), intent(in) :: starts(:, :)
      !> Batch
      type(projection_batch_type), intent(inout) :: batch
      !> Newton scratch
      type(newton_scratch_type), intent(inout) :: scratch
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Start multipliers
      real(wp), intent(in), optional :: lambdas(:)

      type(moist_error_type), allocatable :: failed
      real(wp) :: lambda
      integer :: j, icol

      call batch%init(max(4, size(starts, 2)), failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call batch%begin_run(1, anchor, 1)
      do j = 1, size(starts, 2)
         lambda = 0.0_wp
         if (present(lambdas)) lambda = lambdas(j)
         call batch%append(starts(:, j), lambda, icol)
      end do
      call newton_solve_batch(batch, fx%lsf, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      if (allocated(failed)) call test_failed(error, "Newton failed: "//failed%message)
   end subroutine solve_columns

   !> One column through Newton and the polish
   !>
   !> @param[inout] fx     Fixture
   !> @param[in]    anchor Anchor
   !> @param[in]    seed   Start of the column
   !> @param[in]    lambda Multiplier of the start
   !> @param[out]   status Column status after the polish
   !> @param[out]   x      Column position after the polish
   !> @param[out]   error  Error handle
   subroutine polish_column(fx, anchor, seed, lambda, status, x, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Anchor, start and its multiplier
      real(wp), intent(in) :: anchor(3), seed(3), lambda
      !> Status and position after the polish
      integer, intent(out) :: status
      real(wp), intent(out) :: x(3)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed

      status = point_failed_iterations
      x = 0.0_wp
      call solve_columns(fx, anchor, reshape(seed, [3, 1]), batch, scratch, error, lambdas=[lambda])
      if (allocated(error)) return
      call check(error, batch%status(1), point_converged, "Newton converges")
      if (allocated(error)) return
      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      status = batch%status(1)
      x = batch%x(:, 1)
   end subroutine polish_column

   !> Classification and escape walk at the far pole of a sphere
   !>
   !> - anchor on the axis at height `h`; the point is `-R` on the axis
   !>
   !> @param[inout] fx     Single-atom fixture, polish and `phi_alpha` as set by the caller
   !> @param[in]    rs     Sphere radius
   !> @param[in]    h      Anchor height above the centre
   !> @param[out]   saddle Check's verdict
   !> @param[out]   walked Distance the escape walk travelled from the point
   !> @param[out]   error  Error handle
   subroutine far_pole_regime(fx, rs, h, saddle, walked, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Sphere radius and anchor height
      real(wp), intent(in) :: rs, h
      !> Check's verdict
      logical, intent(out) :: saddle
      !> Distance walked
      real(wp), intent(out) :: walked
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: failed
      real(wp) :: anchor(3), r(3), lambda, r_out(3), lambda_out
      logical :: kkt

      saddle = .false.
      walked = 0.0_wp
      anchor = [h, 0.0_wp, 0.0_wp]
      call on_surface(fx, [-rs, 0.0_wp, 0.0_wp], anchor, r, lambda, error)
      if (allocated(error)) return
      ! A point inside the band raises the third-order error, not a failure of this fixture
      call check_constrained_optimality(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                        r, anchor, lambda, failed, saddle, kkt)
      call check(error, kkt, "the far pole is a KKT point")
      if (allocated(error)) return
      ! An exhausted budget on the way down is of no interest; the distance is
      call riemannian_newton_escape(fx%ctx, fx%param%polish, fx%param%phi_alpha, fx%lsf, &
                                    r, anchor, lambda, r_out, lambda_out, failed)
      walked = norm2(r_out - r)
   end subroutine far_pole_regime

   !> Polish of the far pole ends on the near pole: Newton stops on the maximum, the
   !> check finds the saddle, the walk and a second Newton pass reach the minimum
   !>
   !> @param[inout] fx    Single-atom fixture
   !> @param[in]    rs    Sphere radius
   !> @param[in]    h     Anchor height above the centre
   !> @param[out]   error Error handle
   subroutine polish_far_pole(fx, rs, h, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Sphere radius and anchor height
      real(wp), intent(in) :: rs, h
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(projection_batch_type) :: batch
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed

      call solve_columns(fx, [h, 0.0_wp, 0.0_wp], reshape([-rs - 0.05_wp, 0.0_wp, 0.0_wp], [3, 1]), &
                         batch, scratch, error)
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - [-rs, 0.0_wp, 0.0_wp]) < POS_THR, "Newton stops on the far pole")
      if (allocated(error)) return
      call polish_batch(fx%ctx, batch, fx%lsf, fx%param%polish, fx%param%solve, fx%param%phi_alpha, scratch, failed)
      call check(error, .not. allocated(failed), "polish succeeds")
      if (allocated(error)) return
      call check(error, batch%status(1), point_accepted, "the escaped point is verified")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - [rs, 0.0_wp, 0.0_wp]) < POS_THR, "polish ends on the near pole")
   end subroutine polish_far_pole

   !> Neck point of a carbon pair with an anchor at the focal distance of its meridian curvature
   !>
   !> - the reduced eigenvalue along the meridian is `-phi_alpha eps`, the circumferential
   !>   one is positive (the neck is saddle-shaped); a surface of revolution has no
   !>   mirror symmetry along the meridian, so the third derivative there is nonzero
   !>
   !> @param[out] fx     Fixture
   !> @param[in]  eps    Relative distance of the anchor beyond the focal point (0: on it)
   !> @param[out] r      Point
   !> @param[out] anchor Anchor
   !> @param[out] lambda Multiplier of the point
   !> @param[out] error  Error handle
   subroutine inflection_point(fx, eps, r, anchor, lambda, error)
      !> Fixture
      type(fixture_type), intent(out) :: fx
      !> Focal offset
      real(wp), intent(in) :: eps
      !> Point, anchor, multiplier
      real(wp), intent(out) :: r(3), anchor(3), lambda
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: failed
      real(wp) :: lambda0, p(3, 1), v(1), g(3, 1), h(6, 1), n(3), tc(3), tm(3), km, w

      r = 0.0_wp
      anchor = 0.0_wp
      lambda = 0.0_wp
      call build_fixture(fx, carbon_pair(1.6_wp), error)
      if (allocated(error)) return
      call on_surface(fx, [0.7_wp, 2.1_wp, 0.9_wp], [0.0_wp, 0.0_wp, 0.0_wp], r, lambda0, error)
      if (allocated(error)) return
      p(:, 1) = r
      call fx%lsf%evaluate_batch(p, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      n = g(:, 1)/norm2(g(:, 1))
      ! Circumferential and meridian tangents of the pair's surface of revolution
      tc = [0.0_wp, -r(3), r(2)]/hypot(r(2), r(3))
      tm = cross_product(n, tc)
      km = packed_full_contract2(h(:, 1), tm, 3)
      call check(error, km < 0.0_wp, "the neck curves away from the anchor along the meridian")
      if (allocated(error)) return
      w = fx%param%phi_alpha
      lambda = w*(1.0_wp + eps)/km
      anchor = r - (lambda/w)*g(:, 1)
   end subroutine inflection_point

   !> The reduced curvatures of a surface point have mixed signs
   !>
   !> - `d_neg` is the tangent direction, known from the symmetry, along which the distance
   !>   to the anchor must fall; the other tangent direction must raise it
   !> - probed by short retracted steps, independent of the Hessian the polish classifies with
   !>
   !> @param[inout] fx     Fixture
   !> @param[in]    anchor Anchor
   !> @param[in]    x      Surface point
   !> @param[in]    d_neg  Direction of the negative curvature (made tangent and normalised)
   !> @param[out]   error  Error handle
   subroutine expect_mixed(fx, anchor, x, d_neg, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Anchor and surface point
      real(wp), intent(in) :: anchor(3), x(3)
      !> Direction of the negative curvature
      real(wp), intent(in) :: d_neg(3)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: s, g(3), n(3), dn(3), dp(3), change(4)

      call lsf_value(fx, x, s, g, error)
      if (allocated(error)) return
      n = g/norm2(g)
      dn = d_neg - dot_product(d_neg, n)*n
      dn = dn/norm2(dn)
      dp = cross_product(n, dn)
      call step_change(fx, anchor, x, dn, PROBE_STEP, change(1), error)
      if (allocated(error)) return
      call step_change(fx, anchor, x, -dn, PROBE_STEP, change(2), error)
      if (allocated(error)) return
      call step_change(fx, anchor, x, dp, PROBE_STEP, change(3), error)
      if (allocated(error)) return
      call step_change(fx, anchor, x, -dp, PROBE_STEP, change(4), error)
      if (allocated(error)) return
      call check(error, all(change(1:2) < -PROBE_RISE), &
                 "the distance falls along the negative curvature direction, both ways")
      if (allocated(error)) return
      call check(error, all(change(3:4) > PROBE_RISE), "and rises along the other tangent direction, both ways")
   end subroutine expect_mixed

   !> A surface point is a strict local minimum of the anchor distance along eight tangent rays
   !>
   !> - rays 45 degrees apart; each retracted probe must rise by more than `PROBE_RISE`
   !>
   !> @param[inout] fx     Fixture
   !> @param[in]    anchor Anchor
   !> @param[in]    r      Surface point
   !> @param[in]    normal Unit surface normal at `r`
   !> @param[out]   error  Test failure when a probe is not higher than `r`
   subroutine check_surface_minimum(fx, anchor, r, normal, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Anchor, point and normal
      real(wp), intent(in) :: anchor(3), r(3), normal(3)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: t1(3), t2(3), d(3), rise, ang
      integer :: k

      call setup_tangent_frame(normal, t1, t2)
      do k = 0, 7
         ang = 0.25_wp*pi*real(k, wp)
         d = cos(ang)*t1 + sin(ang)*t2
         call step_change(fx, anchor, r, d, PROBE_STEP, rise, error)
         if (allocated(error)) return
         call check(error, rise > PROBE_RISE, "the objective rises along every probed surface ray")
         if (allocated(error)) return
      end do
   end subroutine check_surface_minimum

   !> Change of the anchor objective after a retracted surface step
   !>
   !> @param[inout] fx     Fixture
   !> @param[in]    anchor Anchor
   !> @param[in]    x      Surface point
   !> @param[in]    d      Step direction
   !> @param[in]    delta  Step length
   !> @param[out]   dphi   Change of the objective
   !> @param[out]   error  Retraction failure
   subroutine step_change(fx, anchor, x, d, delta, dphi, error)
      !> Fixture
      type(fixture_type), intent(inout) :: fx
      !> Anchor, point and direction
      real(wp), intent(in) :: anchor(3), x(3), d(3)
      !> Step length
      real(wp), intent(in) :: delta
      !> Change of the objective
      real(wp), intent(out) :: dphi
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: failed
      real(wp) :: r(3), lambda

      dphi = 0.0_wp
      call retract_to_surface(fx%param%polish, fx%param%phi_alpha, fx%lsf, x + delta*d, anchor, r, lambda, failed)
      if (allocated(failed)) then
         call test_failed(error, "retraction of the probe step failed: "//failed%message)
         return
      end if
      dphi = 0.5_wp*fx%param%phi_alpha*(sum((r - anchor)**2) - sum((x - anchor)**2))
   end subroutine step_change

end module test_cavity_drop_projector
