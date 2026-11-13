!> Regression tests for the batched DROP projection: branch discovery in
!> connected surface patches, cluster branch sets against an exhaustive
!> surface-leaf reference, and the batch, Newton, polish, seeder and finalize
!> building blocks
module test_cavity_drop_projection
   use mctc_env, only: wp, moist_error_type => error_type
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_radii, only: default_cpcm_radii
   use moist_math_linalg, only: setup_tangent_frame
   use moist_math_packed, only: unpack_sym_matrix
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, &
                                                point_captured, point_failed_unreached
   use moist_cavity_drop_projection_newton, only: newton_config_type, newton_scratch_type, &
                                                 newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_config_type, polish_batch
   use moist_cavity_drop_projection_finalize, only: dedupe_candidates, canonical_branch_order, &
                                                   finalize_config_type, finalize_anchor_type, &
                                                   finalize_scratch_type, projection_rows_type, &
                                                   finalize_anchor
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seed_request_type
   use moist_cavity_drop_projection_seed_octree, only: octree_seeder_type
   use moist_cavity_drop_projection_worker, only: projection_worker_type
   implicit none(type, external)
   private

   public :: collect_cavity_drop_projection

   !> Level set evaluator for the unpruned surface-coverage reference
   type :: brute_force_context_type
      !> Geometry evaluator, separate from the finder under test
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
   end type brute_force_context_type

   !> Context for the analytic probe: a union of spheres, as a signed distance
   !>
   !> `S(x) = min_i (||x - c_i|| - R_i)` is exactly 1-Lipschitz, so `|S(x)|` is
   !> a valid surface-free radius -- the same certificate the SvdW level set
   !> supplies, without pulling the whole cavity machinery into a solver test
   type :: sphere_union_context
      !> Sphere centers (3, nsphere)
      real(wp) :: center(3, 2) = 0.0_wp
      !> Sphere radii
      real(wp) :: radius(2) = 1.0_wp
      !> Number of spheres actually used
      integer :: nsphere = 1
   end type sphere_union_context

   !> SvdW level set that counts the points it evaluates
   type, extends(moist_cavity_drop_lsf_svdw_type) :: counting_svdw_type
      !> Points evaluated, every order
      integer :: n_points = 0
   contains
      !> Count, then evaluate as the SvdW level set
      procedure :: evaluate_batch => counting_evaluate_batch
   end type counting_svdw_type

contains

   !> Collect branch discovery regressions
   !>
   !> @param[out] testsuite Projection tests
   subroutine collect_cavity_drop_projection(testsuite)
      !> Registered tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("octree_connected_minima", test_connected_minima), &
         & new_unittest("octree_shallow_minimum", test_shallow_minimum), &
         & new_unittest("octree_search_window", test_search_window), &
         & new_unittest("octree_rigid_motion", test_rigid_motion), &
         & new_unittest("octree_zero_window_on_surface", test_zero_window), &
         & new_unittest("search_window_validation", test_window_validation), &
         & new_unittest("octree_cluster_dimer", test_dimer_branching), &
         & new_unittest("octree_cluster_cross", test_branching_xyz_cross), &
         & new_unittest("octree_cluster_octahedron", test_octahedral), &
         & new_unittest("octree_cluster_cube", test_cube), &
         & new_unittest("octree_cluster_pentagon", test_pentagonal), &
         & new_unittest("octree_cluster_tetrahedron", test_tetrahedral), &
         & new_unittest("brute_force_broad_window", test_brute_force_broad_window), &
         & new_unittest("stationary_box_bounds", test_stationary_box_bounds), &
         & new_unittest("batch_flush_accumulator", test_batch_accumulator), &
         & new_unittest("newton_sphere_foot", test_newton_sphere_foot), &
         & new_unittest("newton_stage_switch", test_newton_stage_switch), &
         & new_unittest("newton_far_seed_unreached", test_newton_far_seed_unreached), &
         & new_unittest("newton_far_foot", test_newton_far_foot), &
         & new_unittest("polish_escapes_far_pole", test_polish_far_pole), &
         & new_unittest("multistart_capture_stops_duplicates", test_capture_stops_duplicates), &
         & new_unittest("multistart_capture_keeps_twin_feet", test_capture_keeps_twin_feet), &
         & new_unittest("finalize_canonical_order", test_finalize_canonical_order), &
         & new_unittest("finalize_sphere_anchor", test_finalize_sphere_anchor), &
         & new_unittest("finalize_branch_weight_fold", test_finalize_branch_fold), &
         & new_unittest("octree_two_sphere_branches", test_octree_two_sphere), &
         & new_unittest("octree_single_sphere_one_branch", test_octree_single_sphere), &
         & new_unittest("octree_sign_change_bound_is_tight", test_octree_sign_change_bound), &
         & new_unittest("octree_certifies_empty_ball", test_octree_empty_ball), &
         & new_unittest("octree_surface_patch_coverage", test_octree_surface_patch_coverage), &
         & new_unittest("octree_optional_region_certificate", test_octree_region_certificate), &
         & new_unittest("octree_seed_anchor_within_trace_tolerance", test_octree_trace_on_surface), &
         & new_unittest("octree_far_verified_foot_radius", test_octree_far_foot_radius), &
         & new_unittest("octree_far_traced_radius", test_octree_far_traced_radius), &
         & new_unittest("octree_seeds_on_sign_change_segments", test_octree_sign_change_seeds), &
         & new_unittest("octree_sign_free_region_no_seeds", test_octree_sign_free_region), &
         & new_unittest("octree_refinement_thin_feature", test_octree_refinement_thin_feature), &
         & new_unittest("bad_octree_seed_size", test_bad_octree_seed_size, should_fail=.true.), &
         & new_unittest("bad_octree_budgets", test_bad_octree_budgets, should_fail=.true.), &
         & new_unittest("bad_octree_max_depth", test_bad_octree_max_depth, should_fail=.true.), &
         & new_unittest("bad_octree_rho_max", test_bad_octree_rho_max, should_fail=.true.), &
         & new_unittest("bad_octree_negative_slack", test_bad_octree_negative_slack, should_fail=.true.), &
         & new_unittest("bad_octree_depth_cap", test_bad_octree_depth_cap, should_fail=.true.), &
         & new_unittest("bad_octree_survivors", test_bad_octree_survivors, should_fail=.true.), &
         & new_unittest("bad_octree_box_budget", test_bad_octree_box_budget, should_fail=.true.), &
         & new_unittest("bad_octree_seed_before_init", test_bad_octree_seed_before_init, should_fail=.true.), &
         & new_unittest("bad_octree_run_before_init", test_bad_octree_run_before_init, should_fail=.true.), &
         & new_unittest("bad_octree_refine_levels", test_bad_octree_refine_levels, should_fail=.true.), &
         & new_unittest("octree_neighbours_outside_box_budget", test_octree_neighbours_outside_budget), &
         & new_unittest("octree_refinement_outside_box_budget", test_octree_refinement_outside_budget), &
         & new_unittest("worker_octree_refinement_by_level", test_worker_octree_refinement), &
         & new_unittest("bad_batch_target_size", test_bad_batch_target_size, should_fail=.true.), &
         & new_unittest("bad_newton_unreached_anchor", test_bad_newton_unreached_anchor, should_fail=.true.) &
         & ]

   end subroutine collect_cavity_drop_projection

   !> Check exclusions against analytic stationary points on spheres
   !>
   !> @param[out] error Invalid whole-box exclusion
   subroutine test_stationary_box_bounds(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Atomic level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Molecular geometry
      type(structure_type) :: mol
      !> Sphere centers, anchor, and box center
      real(wp) :: positions(3, 1), anchor(3), center(3), direction(3)
      !> Rigid translation testing cancellation in coordinate differences
      real(wp) :: shift(3)
      !> Deterministic directions and translations
      integer :: i, trial

      positions = 0.0_wp
      call new(mol, [1], positions)
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      call lsf%update(mol, [1.0_wp])
      anchor = [2.0_wp, 0.0_wp, 0.0_wp]
      center = [0.0_wp, 1.0_wp, 0.0_wp]
      call check(error, lsf%stationary_box_excluded(anchor, center, 0.1_wp), &
                 message="A transverse sphere patch was not excluded")
      if (allocated(error)) return
      center = [1.0_wp, 0.1_wp, 0.0_wp]
      call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.1_wp), &
                 message="Excluded a stationary point on the box boundary")
      if (allocated(error)) return

      do trial = 1, 2
         shift = 0.0_wp
         if (trial == 2) shift = [2.0_wp**38, -2.0_wp**37, 2.0_wp**36]
         positions(:, 1) = shift
         call lsf%set_centers(positions)
         do i = 1, 16
            direction = [cos(real(i, wp)), sin(real(i, wp)), 0.125_wp*real(i, wp)]
            direction = direction/norm2(direction)
            anchor = shift + 2.0_wp*direction
            center = shift + direction
            call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.001_wp), &
                       message="Excluded the analytic nearest sphere point")
            if (allocated(error)) return
            center = shift - direction
            call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.001_wp), &
                       message="Excluded the analytic farthest sphere point")
            if (allocated(error)) return
         end do
      end do

      positions = 0.0_wp
      call lsf%set_centers(positions)
      anchor = [2.0_wp, 0.0_wp, 0.0_wp]
      center = [0.0_wp, 1.0_wp, 0.0_wp]
      lsf%param%blend_2b = -1.0_wp
      call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.1_wp), &
                 message="Applied positive-mixture bound to a negative coefficient")
      if (allocated(error)) return
      lsf%param%blend_2b = 0.0_wp
      lsf%param%blend_k = ieee_value(0.0_wp, ieee_quiet_nan)
      call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.1_wp), &
                 message="Applied bound to nonfinite blend parameters")
      if (allocated(error)) return
      lsf%param%blend_k = 5.5_wp
      lsf%center_bounds_valid = .false.
      call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.1_wp), &
                 message="Applied bound without a valid atomic enclosure")
   end subroutine test_stationary_box_bounds

   !> Recover the farther minimum lost by discrete survivor clustering
   !>
   !> @param[out] error Failed branch comparison
   subroutine test_connected_minima(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Anchor on the lower atomic reference sphere, in Bohr
      real(wp), parameter :: anchor(3) = &
                             [1.4183440527675635_wp, 1.4183440527675635_wp, -1.1316559472324363_wp]
      !> Independent meridian stationary points, in Bohr
      real(wp), parameter :: expected(3, 2) = reshape([ &
                                                      0.8215380815242048_wp, 0.8215380815242048_wp, &
                                                      0.3809228455330751_wp, 1.4183515439720356_wp, &
                                                      1.4183515439720356_wp, -1.1316484571542558_wp], [3, 2])

      call check_h2_branches(anchor, expected, error)
   end subroutine test_connected_minima

   !> Resolve both minima near the neck of the connected surface
   !>
   !> @param[out] error Failed branch comparison
   subroutine test_shallow_minimum(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Anchor close to the neck, in Bohr
      real(wp), parameter :: anchor(3) = &
                             [0.8494343226178617_wp, 0.0_wp, -0.24488375832755738_wp]
      !> Independent meridian stationary points, in Bohr
      real(wp), parameter :: expected(3, 2) = reshape([ &
                                                      0.6839585893716769_wp, 0.0_wp, 0.15815551220714166_wp, &
                                                      0.8557635344337985_wp, 0.0_wp, -0.23042303802218192_wp], [3, 2])

      call check_h2_branches(anchor, expected, error)
   end subroutine test_shallow_minimum

   !> Retain distant strict minima under a window independent of branch weights
   !>
   !> @param[out] error Missing reference minimum
   subroutine test_brute_force_broad_window(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Reference anchor requiring an excess squared distance of 16.83158
      real(wp), parameter :: anchor(3) = &
                             [2.062783409097889_wp, 0.74423217486183_wp, -3.657313249769131_wp]
      !> Independent meridian minima from the 1202-point rule, anchor 831
      real(wp), parameter :: expected(3, 2) = reshape([ &
                                                      0.5813856599086734_wp, 0.2097582868850674_wp, &
                                                      0.1310076148269001_wp, 2.0627834091191093_wp, &
                                                      0.7442321748694862_wp, -3.657313249780524_wp], [3, 2])

      call check_h2_branches(anchor, expected, &
                             error, sigma=0.0025_wp, window=20.0_wp, require_pruning=.true.)
   end subroutine test_brute_force_broad_window

   !> Keep discovery independent of softmax parameters under an explicit window
   !>
   !> @param[out] error Incorrect search support
   subroutine test_search_window(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Reference anchor, in Bohr
      real(wp), parameter :: anchor(3) = &
                             [1.4183440527675635_wp, 1.4183440527675635_wp, -1.1316559472324363_wp]
      !> Independent minima, in Bohr
      real(wp), parameter :: expected(3, 2) = reshape([ &
                                                      0.8215380815242048_wp, 0.8215380815242048_wp, &
                                                      0.3809228455330751_wp, 1.4183515439720356_wp, &
                                                      1.4183515439720356_wp, -1.1316484571542558_wp], [3, 2])

      call check_h2_branches(anchor, expected, error, sigma=0.0025_wp, window=3.2_wp)
      if (allocated(error)) return
      call check_h2_branches(anchor, expected, error, sigma=0.25_wp, window=3.2_wp)
      if (allocated(error)) return
      call check_h2_branches(anchor, expected(:, 2:2), error, window=0.02_wp)
   end subroutine test_search_window

   !> Retain both minima after rotating and translating the whole geometry
   !>
   !> @param[out] error Missing transformed branch
   subroutine test_rigid_motion(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Reference anchor, in Bohr
      real(wp), parameter :: anchor(3) = &
                             [1.4183440527675635_wp, 1.4183440527675635_wp, -1.1316559472324363_wp]
      !> Independent minima, in Bohr
      real(wp), parameter :: expected(3, 2) = reshape([ &
                                                      0.8215380815242048_wp, 0.8215380815242048_wp, &
                                                      0.3809228455330751_wp, 1.4183515439720356_wp, &
                                                      1.4183515439720356_wp, -1.1316484571542558_wp], [3, 2])

      call check_h2_branches(anchor, expected, error, transform=.true.)
   end subroutine test_rigid_motion

   !> A surface anchor is the only branch in a zero-distance window
   !>
   !> @param[out] error Zero-radius search failed
   subroutine test_zero_window(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Projection settings
      type(moist_cavity_drop_parameters_type) :: param
      !> One-sphere level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Level-4 projection worker
      type(projection_worker_type) :: worker
      !> Molecular geometry
      type(structure_type) :: mol
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Atom at the origin and exact point on a unit sphere
      real(wp) :: positions(3, 1), anchor(3)
      !> Returned points
      real(wp), allocatable :: feet(:, :)
      !> Returned point count
      integer :: n_points

      positions = 0.0_wp
      anchor = [1.0_wp, 0.0_wp, 0.0_wp]
      call new(mol, [1], positions)
      param%branch_search_rho2_slack = 0.0_wp
      param%proj_level = 4
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      call lsf%update(mol, [1.0_wp])
      call worker%init(param, lsf, positions, native_error)
      if (.not. allocated(native_error)) call worker%project_anchor(anchor, 1, n_points, native_error, feet=feet)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, 1)
      if (allocated(error)) return
      call check(error, norm2(feet(:, 1) - anchor), 0.0_wp, thr=1.0e-14_wp)
   end subroutine test_zero_window

   !> Validate the search sentinel and prevent nonfinite search bounds
   !>
   !> @param[out] error Invalid configuration accepted or default changed
   subroutine test_window_validation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Parameter set under test
      type(moist_cavity_drop_parameters_type) :: param
      !> Native validation result
      type(moist_error_type), allocatable :: native_error

      call param%compute_derived(native_error)
      call check(error,.not. allocated(native_error))
      if (allocated(error)) return
      call check(error, param%branch_rho2_slack(), &
                 2.0_wp*param%branch_dphi_max/param%phi_alpha, thr=1.0e-14_wp)
      if (allocated(error)) return
      param%branch_search_rho2_slack = -0.5_wp
      call param%compute_derived(native_error)
      call check(error, allocated(native_error))
      if (allocated(error)) return
      param%branch_search_rho2_slack = ieee_value(0.0_wp, ieee_quiet_nan)
      call param%compute_derived(native_error)
      call check(error, allocated(native_error))
      if (allocated(error)) return
      call param%init_defaults()
      call param%compute_derived(native_error)
      call check(error,.not. allocated(native_error))
      if (allocated(error)) return
      call check(error, param%branch_search_rho2_slack, -1.0_wp, thr=1.0e-14_wp)
   end subroutine test_window_validation

   !> Compare the level-4 projection with independently enumerated minima
   !>
   !> Reference construction: solve d/dtheta ||y(theta) - anchor||^2 = 0 on
   !> the exact two-focus meridian, retaining positive tangent curvature
   !> No projected neighbour points enter the production solver
   !>
   !> - the worker runs with the weight cut disabled, so the raw branch set is
   !>   compared rather than the softmax survivors
   !>
   !> @param[in] anchor Initial point, in Bohr
   !> @param[in] expected Distinct strict minima, shape (3, n), in Bohr
   !> @param[out] error Missing minimum or solver failure
   !> @param[in] sigma Optional softmax scale
   !> @param[in] window Optional squared-distance search slack, in Bohr^2
   !> @param[in] transform Rotate and translate the geometry when true
   !> @param[in] require_pruning Require fewer seeds than the exhaustive reference
   subroutine check_h2_branches(anchor, expected, error, sigma, window, transform, require_pruning)
      !> Reference point
      real(wp), intent(in) :: anchor(3)
      !> Independent minima
      real(wp), intent(in) :: expected(:, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Weight model scale
      real(wp), intent(in), optional :: sigma
      !> Search support
      real(wp), intent(in), optional :: window
      !> Apply a rigid motion
      logical, intent(in), optional :: transform
      !> Check algorithmic work reduction on the broad-window fixture
      logical, intent(in), optional :: require_pruning
      !> Projection settings
      type(moist_cavity_drop_parameters_type) :: param
      !> Production level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Level-4 projection worker
      type(projection_worker_type) :: worker
      !> Molecular geometry
      type(structure_type) :: mol
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Atom positions
      real(wp) :: positions(3, 2)
      !> Transformed fixture and proper orthogonal rotation
      real(wp) :: point(3), targets(3, size(expected, 2)), rotation(3, 3), shift(3)
      !> Native hydrogen radius, in Bohr
      real(wp), parameter :: radius = 2.456643962006573_wp
      !> Returned branch count and matching indices
      integer :: n_points, i, j, n_octree
      !> Distance from one reference minimum to the nearest native branch
      real(wp) :: distance, rho_start
      !> Returned and unpruned surface-leaf reference points
      real(wp), allocatable :: feet(:, :), reference(:, :)
      !> Exhaustive surface-leaf count
      integer :: reference_seeds

      positions = 0.0_wp
      positions(3, :) = [-2.55_wp, 2.55_wp]
      point = anchor
      targets = expected
      if (present(transform)) then
         if (transform) then
            rotation = reshape([0.36_wp, 0.48_wp, -0.8_wp, &
                                -0.8_wp, 0.6_wp, 0.0_wp, &
                                0.48_wp, 0.64_wp, 0.6_wp], [3, 3])
            shift = [1.3_wp, -0.7_wp, 2.1_wp]
            positions = matmul(rotation, positions) + spread(shift, 2, 2)
            point = matmul(rotation, anchor) + shift
            targets = matmul(rotation, expected) + spread(shift, 2, size(expected, 2))
         end if
      end if
      call new(mol, [1, 1], positions)
      param%tolerance = 1.0e-12_wp
      param%branch_weight_s = 0.25_wp
      param%proj_level = 4
      if (present(sigma)) param%branch_weight_s = sigma
      if (present(window)) param%branch_search_rho2_slack = window
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      lsf%screening_threshold = param%screening_threshold
      call lsf%update(mol, [radius, radius])
      call worker%init(param, lsf, positions, native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      worker%fcfg%wleb_cut = -1.0_wp
      call worker%project_anchor(point, 1, n_points, native_error, feet=feet)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, size(expected, 2), message="Unexpected or duplicate branches")
      if (allocated(error)) return
      do i = 1, size(expected, 2)
         distance = huge(1.0_wp)
         do j = 1, n_points
            distance = min(distance, norm2(feet(:, j) - targets(:, i)))
         end do
         call check(error, distance < 1.0e-8_wp, message="Octree missed a strict minimum")
         if (allocated(error)) return
      end do

      rho_start = 0.0_wp
      n_octree = 0
      select type (seeder => worker%seeder)
      type is (octree_seeder_type)
         rho_start = seeder%dbg_rho_start
         n_octree = seeder%n_points
      end select
      call brute_force_branches(lsf, point, param, rho_start, reference, native_error, reference_seeds)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, size(reference, 2), size(expected, 2), message="Brute-force H2 count differs")
      if (allocated(error)) return
      call check_branch_sets(reference, feet(:, 1:n_points), error)
      if (allocated(error)) return
      if (present(require_pruning)) then
         if (require_pruning) then
            call check(error, 2*n_octree < reference_seeds, &
                       message="Broad H2 search still refines most surface leaves")
         end if
      end if
   end subroutine check_h2_branches

   !> Compare every unpruned strict minimum with the finder result
   !>
   !> @param[in] reference Brute-force strict minima
   !> @param[in] actual Finder points for the same anchor
   !> @param[out] error Missing reference branch
   subroutine check_branch_sets(reference, actual, error)
      !> Reference points
      real(wp), intent(in) :: reference(:, :)
      !> Finder points
      real(wp), intent(in) :: actual(:, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Point and matching indices
      integer :: i, j, matches
      !> Matching diagnostic
      character(len=256) :: message
      !> Smallest displacement to a finder point
      real(wp) :: closest

      do i = 1, size(reference, 2)
         matches = 0
         closest = huge(1.0_wp)
         do j = 1, size(actual, 2)
            closest = min(closest, norm2(reference(:, i) - actual(:, j)))
            if (norm2(reference(:, i) - actual(:, j)) < 1.0e-6_wp) matches = matches + 1
         end do
         write (message, "(a,i0,a,i0,a,es12.4,a,3es20.10)") &
            "Brute-force branch ", i, ": matches=", matches, ", nearest=", closest, &
            ", point=", reference(:, i)
         ! Raw suppressed candidates can coincide at this matching tolerance
         ! Admissible duplicates are checked separately by the geometry test
         call check(error, closest < 1.0e-6_wp, message=trim(message))
         if (allocated(error)) return
      end do
   end subroutine check_branch_sets

   !> Refine every surface leaf without stationary-point or root-region pruning
   !>
   !> Keeps the original exhaustive seeding path inside the tests
   !> Shared local refinement is numerical; this is not a completeness proof
   !> Nonisolated stationary sets are excluded by the strict-curvature check
   !>
   !> @param[in] lsf Level set of the geometry
   !> @param[in] anchor Reference point
   !> @param[in] param Projection settings
   !> @param[in] radius Initial search radius enclosing the admissible window
   !> @param[out] points Distinct strict minima
   !> @param[out] error Native evaluation or search failure
   !> @param[out] n_seeds Optional exhaustive seed count
   subroutine brute_force_branches(lsf, anchor, param, radius, points, error, n_seeds)
      !> Level set
      type(moist_cavity_drop_lsf_svdw_type), intent(in) :: lsf
      !> Reference point
      real(wp), intent(in) :: anchor(3)
      !> Projection settings
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      !> Surface-search radius
      real(wp), intent(in) :: radius
      !> Reference minima
      real(wp), allocatable, intent(out) :: points(:, :)
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Reference work count
      integer, intent(out), optional :: n_seeds
      !> Surface-only octree
      type(octree_seeder_type) :: tree
      !> Reference geometry callback
      type(brute_force_context_type) :: context
      !> Local refinement of every leaf
      type(projection_batch_type) :: batch
      type(newton_config_type) :: ncfg
      type(polish_config_type) :: pcfg
      type(newton_scratch_type) :: scratch
      !> Candidate points and multipliers
      real(wp), allocatable :: candidates(:, :), multipliers(:), gradients(:, :), hessians(:, :, :)
      !> Distance-window and deduplication mask
      logical, allocatable :: keep(:)
      !> Surface data and tangent Hessian
      real(wp) :: t1(3), t2(3), b11, b12, b22, beta, p(3, 1), v(1), gradient(3)
      !> Seed and accepted point counts
      integer :: i, n, n_strict, icol

      if (present(n_seeds)) n_seeds = 0
      context%lsf = lsf
      call tree%init(error, seed_size=param%octree_seed_size, max_boxes=param%octree_max_boxes, &
                     max_survivors=param%octree_max_survivors, max_depth=param%octree_max_depth)
      if (allocated(error)) return
      p(:, 1) = anchor
      call context%lsf%evaluate_batch(p, 0, v, error=error)
      if (allocated(error)) return
      call tree%run(anchor, v(1), radius, param%branch_rho2_slack(), brute_force_probe, context, error)
      if (allocated(error)) return
      if (present(n_seeds)) n_seeds = tree%n_points

      ncfg%phi_alpha = param%phi_alpha
      ncfg%tol = param%proj_tol
      ncfg%max_iter = param%proj_maxiter
      pcfg%phi_alpha = param%phi_alpha
      call batch%init(tree%n_points + 1, error)
      if (allocated(error)) return
      call batch%begin_run(1, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      do i = 1, tree%n_points
         call batch%append(tree%points(:, i), 0.0_wp, icol)
      end do
      call newton_solve_batch(batch, context%lsf, ncfg, scratch, error)
      if (allocated(error)) return
      call polish_batch(batch, context%lsf, pcfg, ncfg, scratch, error)
      if (allocated(error)) return

      allocate (candidates(3, batch%n), multipliers(batch%n), gradients(3, batch%n), hessians(3, 3, batch%n))
      n = 0
      do i = 1, batch%n
         if (batch%status(i) /= point_accepted) cycle
         if (abs(batch%S(i)) > 1.0e-9_wp) cycle
         n = n + 1
         candidates(:, n) = batch%x(:, i)
         multipliers(n) = batch%lambda(i)
         gradients(:, n) = batch%gradS(:, i)
         call unpack_sym_matrix(batch%hessS(:, i), hessians(:, :, n), 3)
      end do
      if (n == 0) then
         allocate (points(3, 0))
         return
      end if
      allocate (keep(n))
      call dedupe_candidates(candidates(:, 1:n), anchor, param%branch_sep_cut, param%branch_rho2_slack(), keep)
      n_strict = 0
      do i = 1, n
         if (.not. keep(i)) cycle
         gradient = gradients(:, i)
         if (norm2(gradient) < 1.0e-12_wp) cycle
         if (norm2(param%phi_alpha*(candidates(:, i) - anchor) - multipliers(i)*gradient) > 1.0e-8_wp) cycle
         call setup_tangent_frame(gradient/norm2(gradient), t1, t2)
         b11 = param%phi_alpha - multipliers(i)*dot_product(t1, matmul(hessians(:, :, i), t1))
         b12 = -multipliers(i)*dot_product(t1, matmul(hessians(:, :, i), t2))
         b22 = param%phi_alpha - multipliers(i)*dot_product(t2, matmul(hessians(:, :, i), t2))
         beta = 0.5_wp*(b11 + b22 - sqrt((b11 - b22)**2 + 4.0_wp*b12**2))
         if (beta <= 1.0e-8_wp) cycle
         n_strict = n_strict + 1
         candidates(:, n_strict) = candidates(:, i)
      end do
      points = candidates(:, 1:n_strict)
   end subroutine brute_force_branches

   !> Surface-only batched probe for the brute-force reference
   !>
   !> @param[in] points Evaluation points
   !> @param[out] lsf0 Surface values
   !> @param[out] radius Surface-free radii
   !> @param[out] reached Every point is evaluated; an unreached one is an error
   !> @param[inout] context Reference geometry
   !> @param[out] error Evaluation failure
   subroutine brute_force_probe(points, lsf0, radius, reached, context, error)
      !> Evaluation points
      real(wp), intent(in), contiguous :: points(:, :)
      !> Surface values
      real(wp), intent(out), contiguous :: lsf0(:)
      !> Surface-free radii
      real(wp), intent(out), contiguous :: radius(:)
      !> Evaluated points
      logical, intent(out), contiguous :: reached(:)
      !> Reference geometry
      class(*), intent(inout) :: context
      !> Evaluation failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Point index
      integer :: k

      select type (context)
      type is (brute_force_context_type)
         call context%lsf%evaluate_batch(points, 0, lsf0, error=error)
         if (allocated(error)) return
         do k = 1, size(points, 2)
            radius(k) = context%lsf%exclusion_radius(lsf0(k))
            reached(k) = .true.
         end do
      class default
         error stop "brute_force_probe: foreign context"
      end select
   end subroutine brute_force_probe

   !> Check the dimer fixture from the development branching set
   !>
   !> @param[out] error Projection regression failure
   subroutine test_dimer_branching(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Development fixture geometry
      type(structure_type) :: mol

      call new(mol, [6, 6], reshape([ &
                                    0.0000_wp, 0.000_wp, 0.00_wp, &
                                    9.6_wp, 0.000_wp, 0.00_wp], [3, 2]))

      call check_cluster_branches(error, "dimer_branching (C-C 9.6 bohr)", &
                                  mol, 50, 0.8_wp)
   end subroutine test_dimer_branching

   !> Check the cross fixture from the development branching set
   !>
   !> @param[out] error Projection regression failure
   subroutine test_branching_xyz_cross(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Development fixture geometry
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                             0.00_wp, 4.21_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, 4.22_wp, &
                                             0.00_wp, -4.18_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, -4.15_wp, &
                                             0.02_wp, 0.10_wp, -0.20_wp], [3, 5])*aatoau)

      call check_cluster_branches(error, "branching_xyz_cross (5C)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_branching_xyz_cross

   !> Check the octahedron fixture from the development branching set
   !>
   !> @param[out] error Projection regression failure
   subroutine test_octahedral(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Development fixture geometry
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6, 6], reshape([ &
                                                4.000_wp, 0.000_wp, 0.000_wp, &
                                                -4.001_wp, 0.001_wp, 0.000_wp, &
                                                0.000_wp, 4.000_wp, 0.002_wp, &
                                                0.000_wp, -4.000_wp, 0.000_wp, &
                                                0.000_wp, 0.001_wp, 4.001_wp, &
                                                0.001_wp, 0.000_wp, -4.000_wp], [3, 6]))

      call check_cluster_branches(error, "octahedral_6C (R=4 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_octahedral

   !> Check the cube fixture from the development branching set
   !>
   !> @param[out] error Projection regression failure
   subroutine test_cube(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Development fixture geometry
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6, 6, 6, 6], reshape([ &
                                                      3.000_wp, 3.000_wp, 3.000_wp, &
                                                      -3.000_wp, 3.000_wp, 3.000_wp, &
                                                      3.000_wp, -3.000_wp, 3.001_wp, &
                                                      -3.001_wp, -3.001_wp, 3.000_wp, &
                                                      3.000_wp, 3.000_wp, -3.000_wp, &
                                                      -3.001_wp, 3.000_wp, -3.000_wp, &
                                                      3.000_wp, -3.001_wp, -3.001_wp, &
                                                      -3.000_wp, -3.000_wp, -3.000_wp], [3, 8]))

      call check_cluster_branches(error, "cube_8C (a=6 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_cube

   !> Check the pentagon fixture from the development branching set
   !>
   !> @param[out] error Projection regression failure
   subroutine test_pentagonal(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Development fixture geometry
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                             5.000_wp, 0.000_wp, 0.001_wp, &
                                             1.545_wp, 4.755_wp, 0.000_wp, &
                                             -4.045_wp, 2.939_wp, 0.000_wp, &
                                             -4.046_wp, -2.939_wp, 0.001_wp, &
                                             1.544_wp, -4.755_wp, 0.000_wp], [3, 5]))

      call check_cluster_branches(error, "pentagonal_5C (R=5 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_pentagonal

   !> Check the tetrahedron fixture from the development branching set
   !>
   !> @param[out] error Projection regression failure
   subroutine test_tetrahedral(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Development fixture geometry
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6], reshape([ &
                                          3.000_wp, 3.000_wp, 3.000_wp, &
                                          -3.000_wp, -3.000_wp, 3.001_wp, &
                                          -3.000_wp, 3.001_wp, -3.000_wp, &
                                          3.002_wp, -3.001_wp, -3.001_wp], [3, 4]))

      call check_cluster_branches(error, "tetrahedral_4C (edge ~8.5 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_tetrahedral

   !> Require fine-octree recovery of the admissible minima the coarse octree finds
   !>
   !> Level 3 is a comparison set, not a completeness certificate; extra
   !> valid level-4 branches are allowed. Both run with the weight cut
   !> disabled, so raw branch sets are compared
   !>
   !> @param[out] error Missing, duplicate, or invalid branch
   !> @param[in] title Fixture label
   !> @param[in] mol Cluster geometry in bohr
   !> @param[in] nleb Reference sphere quadrature size
   !> @param[in] blend_k Level-set softness
   subroutine check_cluster_branches(error, title, mol, nleb, blend_k)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Fixture label
      character(*), intent(in) :: title
      !> Cluster geometry
      type(structure_type), intent(in) :: mol
      !> Points per reference sphere
      integer, intent(in) :: nleb
      !> Level-set softness
      real(wp), intent(in) :: blend_k
      !> Shared serial context
      type(moist_context_type), target :: ctx
      !> Anchor-phase cavity
      type(cavity_type_drop), allocatable :: cavity
      !> Coarse (level 3) and fine (level 4) octree projections
      type(projection_rows_type) :: comparison, octree
      !> Reference level set for the geometry checks
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Native construction error
      type(moist_error_type), allocatable :: native_error
      !> Candidate bounds per anchor
      integer, allocatable :: first(:), last(:)
      !> Point indices, matching count, anchor slot, relevant comparison rows
      integer :: i, j, aid, matches, nrelevant
      !> Branch comparison message
      character(len=256) :: message

      call new_context(ctx, verbosity=0, nthreads=1)
      allocate (cavity)
      call build_cluster_anchors(cavity, ctx, mol, nleb, blend_k, lsf, native_error)
      if (allocated(native_error)) then
         call test_failed(error, title//": anchors: "//native_error%message)
         return
      end if
      call project_cluster(cavity, 3, comparison, native_error)
      if (allocated(native_error)) then
         call test_failed(error, title//": level 3: "//native_error%message)
         return
      end if
      call project_cluster(cavity, 4, octree, native_error)
      if (allocated(native_error)) then
         call test_failed(error, title//": level 4: "//native_error%message)
         return
      end if
      call check_projection_geometry(cavity, lsf, comparison, error)
      if (allocated(error)) return
      call check_projection_geometry(cavity, lsf, octree, error)
      if (allocated(error)) return

      call row_bounds(octree, cavity%anchor%ngrid, first, last)
      nrelevant = 0
      do i = 1, comparison%n
         if (.not. relevant_row(cavity, comparison, i)) cycle
         nrelevant = nrelevant + 1
         aid = comparison%anchor(i)
         matches = 0
         if (first(aid) > 0) then
            do j = first(aid), last(aid)
               if (.not. relevant_row(cavity, octree, j)) cycle
               if (norm2(comparison%xyz(:, i) - octree%xyz(:, j)) <= 1.0e-6_wp) matches = matches + 1
            end do
         end if
         write (message, "(a,a,i0,a,i0)") title, ": anchor ", cavity%anchor%id(aid), &
            ", matching octree branches ", matches
         call check(error, matches, 1, message=trim(message))
         if (allocated(error)) return
      end do
      call check(error, nrelevant > 0, message="Comparison set contains no admissible minimum")
      if (allocated(error)) return
      call check_cluster_brute_force(cavity, lsf, octree, error)
      if (allocated(error)) return
      deallocate (cavity)
      call ctx%delete()
   end subroutine check_cluster_branches

   !> First and last row of every anchor in a row buffer; rows of one anchor are contiguous
   !>
   !> @param[in] rows Row buffer
   !> @param[in] nmax Anchor count
   !> @param[out] first First row per anchor (0: none)
   !> @param[out] last Last row per anchor
   subroutine row_bounds(rows, nmax, first, last)
      !> Row buffer
      type(projection_rows_type), intent(in) :: rows
      !> Anchor count
      integer, intent(in) :: nmax
      !> Row bounds per anchor
      integer, allocatable, intent(out) :: first(:), last(:)
      !> Row and anchor
      integer :: i, aid

      allocate (first(nmax), source=0)
      allocate (last(nmax), source=0)
      do i = 1, rows%n
         aid = rows%anchor(i)
         if (first(aid) == 0) first(aid) = i
         last(aid) = i
      end do
   end subroutine row_bounds

   !> Compare each cluster anchor against exhaustive surface-leaf refinement
   !>
   !> @param[in] cavity Anchor-phase cavity
   !> @param[in] lsf Level set of the geometry
   !> @param[in] rows Raw level-4 rows
   !> @param[out] error Missing reference branch
   subroutine check_cluster_brute_force(cavity, lsf, rows, error)
      !> Anchors and settings
      type(cavity_type_drop), intent(in) :: cavity
      !> Level set
      type(moist_cavity_drop_lsf_svdw_type), intent(in) :: lsf
      !> Raw projected points
      type(projection_rows_type), intent(in) :: rows
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Native reference failure
      type(moist_error_type), allocatable :: native_error
      !> Reference minima
      real(wp), allocatable :: reference(:, :)
      !> Anchor and a conservative radius from a known surface point
      real(wp) :: anchor(3), radius
      !> Row bounds and total reference count
      integer, allocatable :: first(:), last(:)
      integer :: aid, j, total

      call row_bounds(rows, cavity%anchor%ngrid, first, last)
      total = 0
      do aid = 1, cavity%anchor%ngrid
         if (first(aid) == 0) cycle
         anchor = cavity%anchor%xyz(:, aid)
         radius = huge(1.0_wp)
         do j = first(aid), last(aid)
            radius = min(radius, sum((rows%xyz(:, j) - anchor)**2))
         end do
         radius = sqrt(radius + cavity%param%branch_rho2_slack())
         call brute_force_branches(lsf, anchor, cavity%param, radius, reference, native_error)
         if (allocated(native_error)) then
            call test_failed(error, native_error%message)
            return
         end if
         total = total + size(reference, 2)
         call check_branch_sets(reference, rows%xyz(:, first(aid):last(aid)), error)
         if (allocated(error)) return
      end do
      call check(error, total > 0, message="Brute-force cluster reference is empty")
   end subroutine check_cluster_brute_force

   !> Build the anchor phase of a cluster cavity with the development fixtures' settings
   !>
   !> @param[out] cavity Cavity after the pre-filter
   !> @param[in,out] ctx Borrowed context
   !> @param[in] mol Cluster geometry in bohr
   !> @param[in] nleb Reference quadrature size
   !> @param[in] blend_k Level-set softness
   !> @param[out] lsf The cavity's level set, updated for the geometry
   !> @param[out] error Native construction failure
   subroutine build_cluster_anchors(cavity, ctx, mol, nleb, blend_k, lsf, error)
      !> Anchor-phase cavity
      type(cavity_type_drop), intent(out) :: cavity
      !> Borrowed construction context
      type(moist_context_type), intent(inout), target :: ctx
      !> Cluster geometry
      type(structure_type), intent(in) :: mol
      !> Points per reference sphere
      integer, intent(in) :: nleb
      !> Level-set softness
      real(wp), intent(in) :: blend_k
      !> Level set of the geometry
      type(moist_cavity_drop_lsf_svdw_type), intent(out) :: lsf
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Numerical settings shared by both searches
      type(moist_cavity_drop_parameters_type) :: param
      !> Anchors dropped by the pre-filter
      integer :: nremoved

      param%num_leb = nleb
      param%proj_maxiter = 500
      param%tolerance = 1.0e-12_wp
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k, blend_3b=1.0_wp))
      call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=lsf, &
                           param=param, error=error, ctx=ctx)
      if (allocated(error)) return
      cavity%nsph = mol%nat
      call cavity%radius_model%update(mol, error)
      if (allocated(error)) return
      cavity%radii = cavity%radius_model%f0
      cavity%mol = mol
      cavity%sphxyz = mol%xyz
      call cavity%lsf_model%update(mol, cavity%radii)
      call cavity%ensure_lebedev_cache(error)
      if (allocated(error)) return
      call cavity%anchor%fill(cavity%sphxyz, cavity%radii, cavity%ang_grid, cavity%ang_weight)
      call cavity%iswig%update(mol, cavity%radii, wleb_max=maxval(cavity%anchor%wleb0))
      call cavity%anchor%compute_switching(cavity%iswig, ctx%get_num_threads())
      call cavity%anchor%filter(cavity%param%wleb_cut, nremoved)
      lsf%screening_threshold = cavity%param%screening_threshold
      call lsf%update(mol, cavity%radii)
   end subroutine build_cluster_anchors

   !> Raw branch sets of every anchor at one level, with the fixture's Lebedev
   !> weights and switching values
   !>
   !> @param[in] cavity Anchor-phase cavity
   !> @param[in] level Projection level
   !> @param[out] rows Raw rows, weight cut disabled
   !> @param[out] error Native failure
   subroutine project_cluster(cavity, level, rows, error)
      !> Anchor-phase cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Projection level
      integer, intent(in) :: level
      !> Raw rows
      type(projection_rows_type), intent(out) :: rows
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Single worker for the whole anchor set
      type(projection_worker_type) :: worker
      !> Projection settings at this level
      type(moist_cavity_drop_parameters_type) :: param
      !> Per-anchor outputs
      integer, allocatable :: nrows(:), slots(:)
      logical, allocatable :: lost(:)
      integer :: i

      param = cavity%param
      param%proj_level = level
      call worker%init(param, cavity%lsf_model, cavity%mol%xyz, error)
      if (allocated(error)) return
      worker%fcfg%wleb_cut = -1.0_wp
      allocate (nrows(cavity%anchor%ngrid), lost(cavity%anchor%ngrid))
      slots = [(i, i=1, cavity%anchor%ngrid)]
      call worker%process_chunk(slots, cavity%anchor%xyz(:, 1:cavity%anchor%ngrid), &
                                cavity%anchor%owner(1:cavity%anchor%ngrid), &
                                cavity%anchor%wleb0(1:cavity%anchor%ngrid), &
                                cavity%anchor%f0(1:cavity%anchor%ngrid), nrows, lost, error)
      if (allocated(error)) return
      rows = worker%rows
   end subroutine project_cluster

   !> Whether a raw row has a nonnegligible geometric contribution
   !>
   !> @param[in] cavity Settings
   !> @param[in] rows Raw rows
   !> @param[in] i Row index
   pure function relevant_row(cavity, rows, i) result(relevant)
      !> Settings
      type(cavity_type_drop), intent(in) :: cavity
      !> Raw rows
      type(projection_rows_type), intent(in) :: rows
      !> Row index
      integer, intent(in) :: i
      !> Positive geometric contribution before outgoing branch fractions
      logical :: relevant

      relevant = abs(cavity%anchor%wleb0(rows%anchor(i))*cavity%anchor%f0(rows%anchor(i))* &
                     rows%cpjac(i)*rows%w_f(i)) > cavity%param%wleb_cut
   end function relevant_row

   !> Verify surface membership, stationarity, and strict local minimality
   !>
   !> Checks use the defining equations, independently of branch matching
   !>
   !> @param[in] cavity Anchors and settings
   !> @param[inout] lsf Level set evaluator
   !> @param[in] rows Raw rows
   !> @param[out] error Invalid or duplicate admissible representation
   subroutine check_projection_geometry(cavity, lsf, rows, error)
      !> Anchors and settings
      type(cavity_type_drop), intent(in) :: cavity
      !> Level set evaluator
      type(moist_cavity_drop_lsf_svdw_type), intent(inout) :: lsf
      !> Raw rows
      type(projection_rows_type), intent(in) :: rows
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Native evaluation failure
      type(moist_error_type), allocatable :: native_error
      !> Level-set value, gradient and Hessian
      real(wp) :: p(3, 1), v(1), g(3, 1), h(6, 1), gradient(3), hessian(3, 3)
      !> Tangent frame and constrained Hessian entries
      real(wp) :: tangent1(3), tangent2(3), b11, b12, b22, beta_min
      !> Stationarity residual and objective coefficient
      real(wp) :: residual(3), alpha
      !> Point and sibling indices
      integer :: i, j

      alpha = cavity%param%phi_alpha
      do i = 1, rows%n
         if (.not. relevant_row(cavity, rows, i)) cycle
         p(:, 1) = rows%xyz(:, i)
         call lsf%evaluate_batch(p, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=native_error)
         if (allocated(native_error)) then
            call test_failed(error, native_error%message)
            return
         end if
         gradient = g(:, 1)
         call unpack_sym_matrix(h(:, 1), hessian, 3)
         residual = alpha*(rows%xyz(:, i) - cavity%anchor%xyz(:, rows%anchor(i))) - rows%lambda(i)*gradient
         call check(error, abs(v(1)) < 1.0e-9_wp, message="Projected point is off the surface")
         if (allocated(error)) return
         call check(error, norm2(residual) < 1.0e-8_wp, message="Projected point is not stationary")
         if (allocated(error)) return
         call setup_tangent_frame(gradient/norm2(gradient), tangent1, tangent2)
         b11 = alpha - rows%lambda(i)*dot_product(tangent1, matmul(hessian, tangent1))
         b12 = -rows%lambda(i)*dot_product(tangent1, matmul(hessian, tangent2))
         b22 = alpha - rows%lambda(i)*dot_product(tangent2, matmul(hessian, tangent2))
         beta_min = 0.5_wp*(b11 + b22 - sqrt((b11 - b22)**2 + 4.0_wp*b12**2))
         call check(error, beta_min > 1.0e-8_wp, message="Admissible point is not a strict local minimum")
         if (allocated(error)) return
         do j = i + 1, rows%n
            if (rows%anchor(j) /= rows%anchor(i)) exit
            if (.not. relevant_row(cavity, rows, j)) cycle
            call check(error, norm2(rows%xyz(:, i) - rows%xyz(:, j)) > 1.0e-7_wp, &
                       message="Duplicate admissible branch")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_projection_geometry

   !* ================================================================================= *!
   !*                            Batched projection modules                             *!
   !* ================================================================================= *!

   !> Flush accumulator: anchors are never split, an oversize anchor fills a batch alone
   !>
   !> @param[out] error Accumulator violation
   subroutine test_batch_accumulator(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Batch under test
      type(projection_batch_type) :: batch
      !> Native failure
      type(moist_error_type), allocatable :: failed
      !> Points per anchor and the batch each anchor must land in
      integer, parameter :: sizes(6) = [2, 1, 3, 5, 1, 1]
      integer, parameter :: expected_batch(6) = [1, 1, 2, 3, 4, 4]
      !> Anchor, point, column, batch counters
      integer :: ia, k, icol, ibatch, cap, irun
      !> Recorded anchor jet
      real(wp) :: S, grad(3)

      call batch%init(0, failed)
      call check(error, allocated(failed), "a non-positive target is refused")
      if (allocated(error)) return
      call batch%init(4, failed)
      call check(error, .not. allocated(failed), "positive target accepted")
      if (allocated(error)) return

      ibatch = 1
      do ia = 1, size(sizes)
         if (batch%would_overflow(sizes(ia))) then
            call check_runs_contiguous(batch, error)
            if (allocated(error)) return
            call batch%clear()
            ibatch = ibatch + 1
         end if
         call check(error, ibatch, expected_batch(ia), "anchor lands in the expected flush")
         if (allocated(error)) return
         call batch%begin_run(ia, [real(ia, wp), 0.0_wp, 0.0_wp], 1)
         do k = 1, sizes(ia)
            call batch%append([real(ia, wp), real(k, wp), 0.0_wp], 0.0_wp, icol)
         end do
         call check(error, batch%run_len(batch%nrun), sizes(ia), "run holds all of its anchor's points")
         if (allocated(error)) return
         if (ia == 4) then
            call check(error, batch%n, 5, "oversize anchor accepted whole into an empty batch")
            if (allocated(error)) return
         end if
      end do
      call check_runs_contiguous(batch, error)
      if (allocated(error)) return

      call check(error, all(batch%status(1:batch%n) == point_pending), "appended points are pending")
      if (allocated(error)) return
      call check(error, all(ieee_is_nan(batch%S(1:batch%n))), "jets start as NaN, not zero")
      if (allocated(error)) return

      cap = batch%capacity
      call batch%clear()
      call check(error, batch%n == 0 .and. batch%nrun == 0 .and. batch%capacity == cap, &
                 "clear keeps the storage")
      if (allocated(error)) return

      call batch%begin_run(7, [0.0_wp, 0.0_wp, 0.0_wp], 2)
      call batch%append([0.0_wp, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      irun = batch%run_of(icol)
      call batch%set_anchor_jet(irun, 0.5_wp, [1.0_wp, 2.0_wp, 3.0_wp])
      call batch%anchor_jet(irun, S, grad)
      call check(error, S == 0.5_wp .and. all(grad == [1.0_wp, 2.0_wp, 3.0_wp]), "anchor jet round trip")
   end subroutine test_batch_accumulator

   !> Every run is a contiguous, gap-free column range
   !>
   !> @param[in]  batch Batch
   !> @param[out] error Broken run layout
   subroutine check_runs_contiguous(batch, error)
      !> Batch
      type(projection_batch_type), intent(in) :: batch
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Run and column
      integer :: irun, j, next

      next = 1
      do irun = 1, batch%nrun
         call check(error, batch%run_start(irun), next, "runs follow each other without gaps")
         if (allocated(error)) return
         do j = batch%run_start(irun), batch%run_start(irun) + batch%run_len(irun) - 1
            call check(error, batch%run_of(j), irun, "column tagged with its run")
            if (allocated(error)) return
         end do
         next = next + batch%run_len(irun)
      end do
      call check(error, next - 1, batch%n, "runs cover every column")
   end subroutine check_runs_contiguous

   !> Newton from the anchor lands on the analytic foot of a single-atom sphere
   !>
   !> @param[out] error Wrong foot or multiplier
   subroutine test_newton_sphere_foot(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: cfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchors, feet
      real(wp) :: rs, anchors(3, 4), foot(3), u(3), residual(3)
      !> Steps per stage
      integer :: n1(4), n2(4)
      !> Anchor and column
      integer :: k, icol

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return

      ! Outside on an axis, outside in a generic direction, inside, far outside
      anchors(:, 1) = [rs + 2.0_wp, 0.0_wp, 0.0_wp]
      anchors(:, 2) = (rs + 0.7_wp)*unit_vector([0.3_wp, -1.0_wp, 1.1_wp])
      anchors(:, 3) = 0.6_wp*rs*unit_vector([-0.4_wp, 0.5_wp, -0.2_wp])
      anchors(:, 4) = 2.5_wp*rs*unit_vector([0.2_wp, 0.9_wp, -0.6_wp])

      call batch%init(8, failed)
      do k = 1, 4
         call batch%begin_run(k, anchors(:, k), 1)
         call batch%append(anchors(:, k), 0.0_wp, icol)
      end do
      call newton_solve_batch(batch, lsf, cfg, scratch, failed, stage1_steps=n1, stage2_steps=n2)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if

      do k = 1, 4
         call check(error, batch%status(k), point_converged, "anchor converges")
         if (allocated(error)) return
         u = unit_vector(anchors(:, k))
         foot = rs*u
         call check(error, norm2(batch%x(:, k) - foot) < 1.0e-9_wp, "foot is the radial projection")
         if (allocated(error)) return
         call check(error, abs(batch%S(k)) < 1.0e-10_wp, "final jet is on the surface")
         if (allocated(error)) return
         residual = cfg%phi_alpha*(batch%x(:, k) - anchors(:, k)) - batch%lambda(k)*batch%gradS(:, k)
         call check(error, norm2(residual) < 1.0e-9_wp, "final multiplier is stationary")
         if (allocated(error)) return
         call check(error, abs(batch%rho(k) - abs(norm2(anchors(:, k)) - rs)) < 1.0e-9_wp, "rho")
         if (allocated(error)) return
         call check(error, abs(abs(dot_product(batch%normal(:, k), u)) - 1.0_wp) < 1.0e-10_wp, &
                    "normal is radial")
         if (allocated(error)) return
         call check(error, .not. any(ieee_is_nan(batch%hessS(:, k))), "final jet carries the Hessian")
         if (allocated(error)) return
         call check(error, batch%run_jet(k), "anchor jet recorded without an extra evaluation")
         if (allocated(error)) return
         call check(error, n1(k) + n2(k) <= 20, "converges in a handful of steps")
         if (allocated(error)) return
      end do
   end subroutine test_newton_sphere_foot

   !> A seed beyond the level set reach fails alone, its run's anchor seed still converges
   !>
   !> @param[out] error Far seed aborted the batch or was not marked unreached
   subroutine test_newton_far_seed_unreached(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: cfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, far seed and its value
      real(wp) :: rs, anchor(3), far(3, 1), v(1)
      !> Columns
      integer :: icol_anchor, icol_far

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return

      anchor = [rs + 2.0_wp, 0.0_wp, 0.0_wp]
      far(:, 1) = [80.0_wp, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(far, 0, v, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "far probe: "//failed%message)
         return
      end if
      call check(error, lsf%unreached(1), "fixture seed lies beyond the screening reach")
      if (allocated(error)) return

      call batch%init(8, failed)
      call batch%begin_run(1, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol_anchor)
      call batch%append(far(:, 1), 0.0_wp, icol_far)
      call newton_solve_batch(batch, lsf, cfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, batch%status(icol_far), point_failed_unreached, "far seed is unreached")
      if (allocated(error)) return
      call check(error, batch%status(icol_anchor), point_converged, "anchor seed converges")
      if (allocated(error)) return
      call check(error, ieee_is_nan(batch%S(icol_far)), "unreached jet is never stored")
   end subroutine test_newton_far_seed_unreached

   !> Newton from an anchor far outside the sphere walks the whole distance to its foot
   !>
   !> - the foot lies more than 15 Bohr from the anchor along x, many `max_step`
   !>   steps away
   !>
   !> @param[out] error Far foot not reached
   subroutine test_newton_far_foot(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: cfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, analytic foot and anchor value
      real(wp) :: rs, anchor(3, 1), foot(3), v(1)
      !> Column
      integer :: icol
      !> Absolute foot tolerance (Bohr)
      real(wp), parameter :: tol = 1.0e-9_wp

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return

      anchor(:, 1) = [rs + 18.0_wp, 0.0_wp, 0.0_wp]
      foot = [rs, 0.0_wp, 0.0_wp]
      call check(error, abs(anchor(1, 1) - foot(1)) > 15.0_wp, "fixture foot lies beyond 15 Bohr")
      if (allocated(error)) return
      call lsf%evaluate_batch(anchor, 0, v, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "anchor probe: "//failed%message)
         return
      end if
      call check(error, .not. lsf%unreached(1), "fixture anchor lies within the screening reach")
      if (allocated(error)) return

      call batch%init(4, failed)
      call batch%begin_run(1, anchor(:, 1), 1)
      call batch%append(anchor(:, 1), 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, cfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, batch%status(icol), point_converged, "far anchor converges")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, icol) - foot) <= tol .or. &
                 norm2(batch%x(:, icol) - foot) <= 10.0_wp*tol*rs, "foot is the radial projection")
      if (allocated(error)) return
      call check(error, abs(batch%rho(icol) - 18.0_wp) <= tol .or. &
                 abs(batch%rho(icol) - 18.0_wp) <= 10.0_wp*tol*18.0_wp, "rho is the full distance")
   end subroutine test_newton_far_foot

   !> The stage-1 exits: residual threshold, the cheap-step budget, and a zero
   !> budget
   !>
   !> @param[out] error Wrong stage accounting
   subroutine test_newton_stage_switch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: cfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, seed
      real(wp) :: rs, anchor(3), seed(3)
      !> Steps per stage
      integer :: n1(1), n2(1)
      !> Column
      integer :: icol

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      call batch%init(4, failed)

      !> Off-axis seed, near anchor: cheap steps contract, the residual rule switches
      !> (from the anchor itself the radial cheap step is already exact on a sphere)
      anchor = [rs + 1.0_wp, 0.0_wp, 0.0_wp]
      seed = rs*[cos(0.3_wp), sin(0.3_wp), 0.0_wp]
      cfg = newton_config_type(switch_fnorm=0.1_wp, max_cheap_steps=100)
      call batch%begin_run(1, anchor, 1)
      call batch%append(seed, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, cfg, scratch, failed, stage1_steps=n1, stage2_steps=n2)
      call check_stage_run(error, "residual switch", failed, batch, rs*unit_vector(anchor), n1(1), n2(1), 1, 20)
      if (allocated(error)) return

      !> Far anchor, seed off the axis, residual switch off: the budget switches
      call batch%clear()
      anchor = [3.0_wp*rs, 0.0_wp, 0.0_wp]
      seed = rs*[cos(0.5_wp), sin(0.5_wp), 0.0_wp]
      cfg = newton_config_type(switch_fnorm=0.0_wp, max_cheap_steps=3)
      call batch%begin_run(1, anchor, 1)
      call batch%append(seed, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, cfg, scratch, failed, stage1_steps=n1, stage2_steps=n2)
      call check_stage_run(error, "cheap-step budget", failed, batch, rs*unit_vector(anchor), n1(1), n2(1), 3, 3)
      if (allocated(error)) return

      !> No cheap step at all
      call batch%clear()
      anchor = [rs + 1.0_wp, 0.3_wp, -0.2_wp]
      cfg = newton_config_type(max_cheap_steps=0)
      call batch%begin_run(1, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, cfg, scratch, failed, stage1_steps=n1, stage2_steps=n2)
      call check_stage_run(error, "budget switch", failed, batch, rs*unit_vector(anchor), n1(1), n2(1), 0, 0)
   end subroutine test_newton_stage_switch

   !> Check one single-column stage run
   !>
   !> @param[out] error  Test failure
   !> @param[in]  label  Case label
   !> @param[in]  failed Native Newton failure
   !> @param[in]  batch  Solved batch
   !> @param[in]  foot   Analytic foot
   !> @param[in]  n1     Cheap steps taken
   !> @param[in]  n2     Full steps taken
   !> @param[in]  n1_min Fewest cheap steps expected
   !> @param[in]  n1_max Most cheap steps expected
   subroutine check_stage_run(error, label, failed, batch, foot, n1, n2, n1_min, n1_max)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Case label
      character(len=*), intent(in) :: label
      !> Native Newton failure
      type(moist_error_type), allocatable, intent(in) :: failed
      !> Solved batch
      type(projection_batch_type), intent(in) :: batch
      !> Analytic foot
      real(wp), intent(in) :: foot(3)
      !> Steps taken per stage and the cheap-step window
      integer, intent(in) :: n1, n2, n1_min, n1_max
      !> Message buffer
      character(len=160) :: msg

      if (allocated(failed)) then
         call test_failed(error, label//": "//failed%message)
         return
      end if
      write (msg, "(a,a,i0,a,i0)") label, ": cheap steps ", n1, ", full steps ", n2
      call check(error, batch%status(1), point_converged, trim(msg)//": converged")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - foot) < 1.0e-9_wp, trim(msg)//": analytic foot")
      if (allocated(error)) return
      call check(error, n1 >= n1_min .and. n1 <= n1_max, trim(msg)//": cheap steps in window")
      if (allocated(error)) return
      call check(error, n2 >= 1, trim(msg)//": the full system finishes")
   end subroutine check_stage_run

   !> A Newton point stuck on the far pole (a constrained maximum) is escaped
   !> and re-polished onto the near pole; the anchor's own foot is accepted as is
   !>
   !> @param[out] error Wrong classification or escape
   subroutine test_polish_far_pole(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: ncfg
      type(polish_config_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, near pole
      real(wp) :: rs, anchor(3), pole(3)
      !> Column
      integer :: icol

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      anchor = [rs + 2.0_wp, 0.0_wp, 0.0_wp]
      pole = [rs, 0.0_wp, 0.0_wp]

      call batch%init(4, failed)
      call batch%begin_run(1, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      call batch%append([-rs, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:2) == point_converged), "both starts converge")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 2) + pole) < 1.0e-9_wp, "second start sits on the far pole")
      if (allocated(error)) return

      call polish_batch(batch, lsf, pcfg, ncfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "polish: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:2) == point_accepted), "both end as verified minima")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - pole) < 1.0e-9_wp, "anchor foot untouched")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 2) - pole) < 1.0e-8_wp, "far pole escaped onto the near pole")
   end subroutine test_polish_far_pole

   !> The seeded pass stops starts that head into the accepted local foot
   !>
   !> - one sphere at level 3: every octree start of the anchor ends on the
   !>   local foot, so all of them are duplicates
   !> - the same worker with capture off is the reference: the same rows, more
   !>   level set evaluations, nothing captured
   !>
   !> @param[out] error Nothing captured, no evaluations saved, or different rows
   subroutine test_capture_stops_duplicates(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Projection settings
      type(moist_cavity_drop_parameters_type) :: param
      !> Counting one-sphere level set the worker clones
      type(counting_svdw_type) :: lsf
      !> Worker under test
      type(projection_worker_type) :: worker
      !> One-atom geometry
      type(structure_type) :: mol
      type(moist_error_type), allocatable :: failed
      !> Atom position, anchor
      real(wp) :: positions(3, 1), anchor(3)
      !> Rows with capture on (1) and off (2)
      real(wp), allocatable :: feet_on(:, :), feet_off(:, :), w_on(:), w_off(:)
      !> Rows, captured starts and evaluated points per run
      integer :: n_on, n_off, captured(2), evaluated(2), k

      positions = 0.0_wp
      call new(mol, [6], positions)
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, blend_3b=1.0_wp))
      lsf%screening_threshold = 1.0e-11_wp
      call lsf%update(mol, [3.0_wp])
      anchor = 4.0_wp*unit_vector([0.7_wp, -0.4_wp, 0.3_wp])

      param%num_leb = 302
      param%proj_level = 3
      call param%compute_derived(failed)
      if (.not. allocated(failed)) call worker%init(param, lsf, positions, failed)
      if (allocated(failed)) then
         call test_failed(error, "worker: "//failed%message)
         return
      end if

      call project(.true., 1, n_on, feet_on, w_on)
      if (allocated(error)) return
      call project(.false., 2, n_off, feet_off, w_off)
      if (allocated(error)) return

      call check(error, captured(1) > 0, "the seeded pass captures duplicate starts")
      if (allocated(error)) return
      call check(error, captured(2), 0, "nothing is captured with capture off")
      if (allocated(error)) return
      call check(error, evaluated(1) < evaluated(2), "captured starts save level set evaluations")
      if (allocated(error)) return
      call check(error, n_on >= 1 .and. n_on == n_off, "the same number of rows")
      if (allocated(error)) return
      do k = 1, n_on
         call check(error, same_value(feet_on(:, k), feet_off(:, k), 1.0e-9_wp), "the same feet")
         if (allocated(error)) return
         call check(error, same_value([w_on(k)], [w_off(k)], 1.0e-9_wp), "the same branch weights")
         if (allocated(error)) return
      end do

   contains

      !> Project the anchor with capture on or off and record the counts
      !>
      !> @param[in]  capture Capture duplicates in the seeded pass
      !> @param[in]  run     Slot of the counts
      !> @param[out] n       Rows
      !> @param[out] feet    Feet (3, n)
      !> @param[out] w       Branch weights
      subroutine project(capture, run, n, feet, w)
         logical, intent(in) :: capture
         integer, intent(in) :: run
         integer, intent(out) :: n
         real(wp), allocatable, intent(out) :: feet(:, :)
         real(wp), allocatable, intent(out) :: w(:)
         integer :: n0

         worker%ncfg_seeded%capture_duplicates = capture
         n0 = worker%n_captured
         call reset_count(worker%lsf)
         call worker%project_anchor(anchor, 1, n, failed, feet=feet, weights=w)
         if (allocated(failed)) then
            call test_failed(error, "projection: "//failed%message)
            return
         end if
         captured(run) = worker%n_captured - n0
         select type (clone => worker%lsf)
         type is (counting_svdw_type)
            evaluated(run) = clone%n_points
         class default
            call test_failed(error, "the worker's level set is not the counting clone")
         end select
      end subroutine project

   end subroutine test_capture_stops_duplicates

   !> Close twin feet both survive capture
   !>
   !> - two overlapping carbons; the anchor on their mirror plane lies just
   !>   beyond the focal point of the neck, so its two mirror feet are close
   !> - a run as the seeded pass builds it: the accepted foot on one side
   !>   first, then starts across the whole neck, the middle one on the saddle
   !> - with capture on, the run keeps both feet, a start of each side is
   !>   captured, and the feet are those of the run with capture off
   !>
   !> @param[out] error A twin lost or a foot moved
   subroutine test_capture_keeps_twin_feet(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Two-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Run of the accepted foot, run under test
      type(projection_batch_type) :: local, batch
      type(newton_config_type) :: ncfg
      type(polish_config_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      type(structure_type) :: mol
      !> Anchor, distinct feet with capture on (1) and off (2)
      real(wp) :: anchor(3), feet(3, 32, 2), mirror(3)
      !> Distinct feet and captured starts per run, column, start, match
      integer :: nfeet(2), captured(2), icol, k, run, i, m

      call new(mol, [6, 6], reshape([-2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, blend_3b=1.0_wp))
      lsf%screening_threshold = 1.0e-11_wp
      call lsf%update(mol, [3.0_wp, 3.0_wp])
      anchor = [0.0_wp, 3.67_wp, 0.0_wp]

      ! Accepted foot on the +x side
      call local%init(4, failed)
      call local%begin_run(1, anchor, 1)
      call local%append([0.6_wp, 2.6_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(local, lsf, ncfg, scratch, failed)
      if (.not. allocated(failed)) call polish_batch(local, lsf, pcfg, ncfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "local foot: "//failed%message)
         return
      end if
      call check(error, local%status(1), point_accepted, "local foot verified")
      if (allocated(error)) return

      call batch%init(32, failed)
      do run = 1, 2
         call batch%clear()
         call batch%begin_run(1, anchor, 1)
         call batch%append_copy(local, 1, icol)
         do k = -10, 10
            call batch%append([0.1_wp*k, 2.6_wp, 0.0_wp], 0.0_wp, icol)
         end do
         ncfg%capture_duplicates = run == 1
         call newton_solve_batch(batch, lsf, ncfg, scratch, failed)
         if (.not. allocated(failed)) call polish_batch(batch, lsf, pcfg, ncfg, scratch, failed)
         if (allocated(failed)) then
            call test_failed(error, "seeded run: "//failed%message)
            return
         end if
         captured(run) = count(batch%status(1:batch%n) == point_captured)
         call distinct_feet(batch, feet(:, :, run), nfeet(run))
      end do

      call check(error, captured(1) > 0, "the run captures duplicate starts")
      if (allocated(error)) return
      call check(error, nfeet(1), 2, "both twin feet kept")
      if (allocated(error)) return
      mirror = [-feet(1, 2, 1), feet(2, 2, 1), feet(3, 2, 1)]
      call check(error, same_value(feet(:, 1, 1), mirror, 1.0e-7_wp), "the two feet are mirror images")
      if (allocated(error)) return
      call check(error, norm2(feet(:, 1, 1) - feet(:, 2, 1)) < 0.1_wp, "the fixture's twins are close")
      if (allocated(error)) return
      do i = 1, nfeet(1)
         m = 0
         do k = 1, nfeet(2)
            if (same_value(feet(:, i, 1), feet(:, k, 2), 1.0e-7_wp)) m = k
         end do
         call check(error, m > 0, "every foot is a foot of the run with capture off")
         if (allocated(error)) return
      end do
      do k = 1, nfeet(2)
         m = 0
         do i = 1, nfeet(1)
            if (same_value(feet(:, i, 1), feet(:, k, 2), 1.0e-7_wp)) m = i
         end do
         call check(error, m > 0, "no foot of the run with capture off is lost")
         if (allocated(error)) return
      end do

   contains

      !> Accepted feet of the run, deduplicated as the finalize does
      !>
      !> @param[in]  b     Solved and polished run
      !> @param[out] pts   Distinct feet (3, :)
      !> @param[out] n     Distinct feet found
      subroutine distinct_feet(b, pts, n)
         type(projection_batch_type), intent(in) :: b
         real(wp), intent(out) :: pts(:, :)
         integer, intent(out) :: n
         real(wp), allocatable :: acc(:, :)
         logical, allocatable :: keep(:)
         integer :: j, na

         na = count(b%status(1:b%n) == point_accepted)
         allocate (acc(3, na), keep(na))
         na = 0
         do j = 1, b%n
            if (b%status(j) /= point_accepted) cycle
            na = na + 1
            acc(:, na) = b%x(:, j)
         end do
         call dedupe_candidates(acc, anchor, 1.0e-8_wp, huge(1.0_wp), keep)
         n = 0
         do j = 1, na
            if (.not. keep(j)) cycle
            n = n + 1
            pts(:, n) = acc(:, j)
         end do
      end subroutine distinct_feet

   end subroutine test_capture_keeps_twin_feet

   !> Whether two vectors agree within `tol`, or `10 tol` relative to the reference
   !>
   !> @param[in] a   Value
   !> @param[in] ref Reference
   !> @param[in] tol Absolute tolerance
   pure function same_value(a, ref, tol) result(same)
      !> Value
      real(wp), intent(in) :: a(:)
      !> Reference
      real(wp), intent(in) :: ref(:)
      !> Absolute tolerance
      real(wp), intent(in) :: tol
      !> Agreement
      logical :: same

      real(wp) :: diff

      diff = maxval(abs(a - ref))
      same = diff <= tol .or. diff <= 10.0_wp*tol*maxval(abs(ref))
   end function same_value

   !> Reset the evaluation count of a counting level set
   !>
   !> @param[inout] lsf Level set, a counting one
   subroutine reset_count(lsf)
      !> Level set, a counting one
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf

      select type (lsf)
      type is (counting_svdw_type)
         lsf%n_points = 0
      end select
   end subroutine reset_count

   !> Count the points, then evaluate as the SvdW level set
   !>
   !> @param[inout] self              Counting level set
   !> @param[in]    points            Evaluation points (3, n)
   !> @param[in]    order             Highest spatial derivative to write
   !> @param[out]   lsf0              Values (n)
   !> @param[out]   lsf1_r            Gradients (3, n)
   !> @param[out]   lsf2_rr_packed    Packed Hessians (6, n)
   !> @param[out]   lsf3_rrr          Packed third derivatives (10, n)
   !> @param[out]   lsf4_rrrr         Packed fourth derivatives (15, n)
   !> @param[in]    candidate_indices Shared atom list
   !> @param[out]   error             Evaluation failure
   subroutine counting_evaluate_batch(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                      lsf3_rrr, lsf4_rrrr, candidate_indices, error)
      class(counting_svdw_type), intent(inout) :: self
      real(wp), intent(in), contiguous :: points(:, :)
      integer, intent(in) :: order
      real(wp), intent(out), contiguous :: lsf0(:)
      real(wp), intent(out), optional, contiguous :: lsf1_r(:, :)
      real(wp), intent(out), optional, contiguous :: lsf2_rr_packed(:, :)
      real(wp), intent(out), optional, contiguous :: lsf3_rrr(:, :)
      real(wp), intent(out), optional, contiguous :: lsf4_rrrr(:, :)
      integer, intent(in), optional :: candidate_indices(:)
      type(moist_error_type), allocatable, intent(out) :: error

      self%n_points = self%n_points + size(points, 2)
      call self%moist_cavity_drop_lsf_svdw_type%evaluate_batch(points, order, lsf0, lsf1_r, &
                                                               lsf2_rr_packed, lsf3_rrr, lsf4_rrrr, &
                                                               candidate_indices, error)
   end subroutine counting_evaluate_batch

   !> The canonical branch order is lexicographic with a tie tolerance
   !>
   !> @param[out] error Wrong order
   subroutine test_finalize_canonical_order(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Feet and order
      real(wp) :: pts(3, 3)
      integer :: order(3)

      ! x ties within 1e-8 fall through to y
      pts(:, 1) = [1.0_wp, 2.0_wp, 0.0_wp]
      pts(:, 2) = [1.0_wp + 1.0e-10_wp, -1.0_wp, 0.0_wp]
      pts(:, 3) = [-3.0_wp, 5.0_wp, 0.0_wp]
      call canonical_branch_order(pts, order)
      call check(error, all(order == [3, 2, 1]), "lexicographic order with tie tolerance")
   end subroutine test_finalize_canonical_order

   !> Finalize on a sphere: duplicate feet merge into one row whose CP Jacobian is
   !> the radial area ratio `(rs / R)^2` of the concentric anchor sphere
   !>
   !> @param[out] error Wrong row
   subroutine test_finalize_sphere_anchor(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points, solver and finalize state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: ncfg
      type(polish_config_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(finalize_config_type) :: fcfg
      type(finalize_anchor_type) :: anchor
      type(finalize_scratch_type) :: fsc
      type(projection_rows_type) :: rows
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_wleb
      type(branch_weight_type) :: branch_weight
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor sphere radius, foot
      real(wp) :: rs, foot(3)
      real(wp), parameter :: r_anchor = 3.0_wp
      !> Columns and counts
      integer :: icol, nbranch, nrows

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      anchor%index = 5
      anchor%xyz = r_anchor*unit_vector([0.2_wp, 0.3_wp, 0.9_wp])
      anchor%owner_center = 0.0_wp
      anchor%wleb = 0.01_wp
      anchor%f = 1.0_wp
      foot = rs*unit_vector(anchor%xyz)

      ! The anchor and two nearby starts all reach the same foot
      call batch%init(4, failed)
      call batch%begin_run(anchor%index, anchor%xyz, 1)
      call batch%append(anchor%xyz, 0.0_wp, icol)
      call batch%append(anchor%xyz + [1.0e-3_wp, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call batch%append(foot, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, scratch, failed)
      if (.not. allocated(failed)) call polish_batch(batch, lsf, pcfg, ncfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "Newton/polish: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:3) == point_accepted), "all starts verified")
      if (allocated(error)) return

      call new_smoothstep_sigmoid(f_crit, 0.25_wp, 0.6_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_foc, 0.1_wp, 0.3_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_wleb, 1.0e-12_wp, 1.0e-10_wp)
      call branch_weight%init(0.0025_wp)
      fcfg%rho2_slack = 0.26_wp
      call finalize_anchor(fcfg, batch, 1, 3, anchor, f_crit, f_foc, f_wleb, branch_weight, &
                           rows, nbranch, nrows, fsc, failed)
      if (allocated(failed)) then
         call test_failed(error, "finalize: "//failed%message)
         return
      end if
      call check(error, nbranch == 1 .and. nrows == 1 .and. rows%n == 1, "duplicates merge into one row")
      if (allocated(error)) return
      call check(error, rows%anchor(1) == 5 .and. rows%branch(1) == 1, &
                 "row indexing")
      if (allocated(error)) return
      call check(error, norm2(rows%xyz(:, 1) - foot) < 1.0e-9_wp, "row foot")
      if (allocated(error)) return
      call check(error, abs(rows%cpjac(1) - (rs/r_anchor)**2) < 1.0e-8_wp, "radial area ratio")
      if (allocated(error)) return
      call check(error, rows%w_f(1) > 0.0_wp .and. rows%w_f(1) <= 1.0_wp .and. rows%wbranch(1) == 1.0_wp, &
                 "switching and branch weight")
      if (allocated(error)) return

      ! An anchor below the switching cutoff appends nothing
      anchor%f = 0.0_wp
      call finalize_anchor(fcfg, batch, 1, 3, anchor, f_crit, f_foc, f_wleb, branch_weight, &
                           rows, nbranch, nrows, fsc, failed)
      call check(error, nbranch == 1 .and. nrows == 0 .and. rows%n == 1, "rows below the cut are never appended")
   end subroutine test_finalize_sphere_anchor

   !> Finalize with two mirror branches folds each branch weight into the row
   !> weight: every row is its single-branch weight times its branch weight,
   !> and the rows of the anchor sum to the single-branch weight
   !>
   !> - midpoint of a carbon pair; the single-branch references finalize each
   !>   column alone, where the branch weight is one
   !> - negative control: both branch weights are below one, so a lost fold
   !>   fails the per-row check
   !>
   !> @param[out] error Wrong row weight
   subroutine test_finalize_branch_fold(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Two-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points, solver and finalize state
      type(projection_batch_type) :: batch
      type(newton_config_type) :: ncfg
      type(polish_config_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(finalize_config_type) :: fcfg
      type(finalize_anchor_type) :: anchor
      type(finalize_scratch_type) :: fsc
      type(projection_rows_type) :: rows, single
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_wleb
      type(branch_weight_type) :: branch_weight
      type(moist_error_type), allocatable :: failed
      type(structure_type) :: mol
      !> Single-atom radius, atom offset, reference weight
      real(wp) :: rs, h, w_ref
      !> Columns, counts, rows
      integer :: icol, nbranch, nrows, r, m, match

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      h = 2.2_wp*rs
      call new(mol, [6, 6], reshape([-h, 0.0_wp, 0.0_wp, h, 0.0_wp, 0.0_wp], [3, 2]))
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, blend_3b=1.0_wp))
      lsf%screening_threshold = 1.0e-11_wp
      call lsf%update(mol, [3.0_wp, 3.0_wp])

      anchor%index = 7
      anchor%xyz = 0.0_wp
      anchor%owner_center = [-h, 0.0_wp, 0.0_wp]
      anchor%wleb = 0.01_wp
      anchor%f = 1.0_wp

      ! One start near each atom's near pole
      call batch%init(4, failed)
      call batch%begin_run(anchor%index, anchor%xyz, 1)
      call batch%append([-(h - rs), 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call batch%append([h - rs, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, scratch, failed)
      if (.not. allocated(failed)) call polish_batch(batch, lsf, pcfg, ncfg, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "Newton/polish: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:2) == point_accepted), "both mirror feet verified")
      if (allocated(error)) return

      call new_smoothstep_sigmoid(f_crit, 0.25_wp, 0.6_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_foc, 0.1_wp, 0.3_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_wleb, 1.0e-12_wp, 1.0e-10_wp)
      call branch_weight%init(0.0025_wp)
      fcfg%rho2_slack = 0.26_wp

      call finalize_anchor(fcfg, batch, 1, 2, anchor, f_crit, f_foc, f_wleb, branch_weight, &
                           rows, nbranch, nrows, fsc, failed)
      if (.not. allocated(failed)) then
         call finalize_anchor(fcfg, batch, 1, 1, anchor, f_crit, f_foc, f_wleb, branch_weight, &
                              single, nbranch, nrows, fsc, failed)
      end if
      if (.not. allocated(failed)) then
         call finalize_anchor(fcfg, batch, 2, 2, anchor, f_crit, f_foc, f_wleb, branch_weight, &
                              single, nbranch, nrows, fsc, failed)
      end if
      if (allocated(failed)) then
         call test_failed(error, "finalize: "//failed%message)
         return
      end if
      call check(error, rows%n == 2 .and. single%n == 2, "two kept rows, one per single-branch run")
      if (allocated(error)) return
      call check(error, all(single%wbranch(1:2) == 1.0_wp), "a single branch has branch weight one")
      if (allocated(error)) return
      call check(error, all(rows%wbranch(1:2) < 1.0_wp) .and. abs(sum(rows%wbranch(1:2)) - 1.0_wp) < 1.0e-12_wp, &
                 "the kept branch weights share the anchor")
      if (allocated(error)) return

      do r = 1, rows%n
         match = 0
         do m = 1, single%n
            if (norm2(rows%xyz(:, r) - single%xyz(:, m)) < 1.0e-12_wp) match = m
         end do
         call check(error, match > 0, "every kept row has its single-branch reference")
         if (allocated(error)) return
         w_ref = single%wleb(match)*rows%wbranch(r)
         call check(error, abs(rows%wleb(r) - w_ref) < max(1.0e-15_wp, 1.0e-12_wp*abs(w_ref)), &
                    "row weight is the single-branch weight times the branch weight")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(rows%wleb(1:2)) - single%wleb(1)) < 1.0e-12_wp*single%wleb(1), &
                 "the kept rows sum to the anchor's single-branch weight")
   end subroutine test_finalize_branch_fold

   !> Single carbon at the origin and the radius of its zero level set
   !>
   !> - the radius comes from bisection of `S` along an axis, independent of
   !>   any projection code
   !>
   !> @param[out] lsf   Single-atom SvdW level set
   !> @param[out] rs    Surface radius
   !> @param[out] error Bracket failure
   subroutine sphere_lsf(lsf, rs, error)
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type), intent(out) :: lsf
      !> Surface radius
      real(wp), intent(out) :: rs
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> One-atom geometry
      type(structure_type) :: mol
      !> Bracket, probe point and value
      real(wp) :: lo, hi, p(3, 1), v(1), v_lo
      !> Native failure
      type(moist_error_type), allocatable :: failed
      !> Bisection step
      integer :: it

      call new(mol, [6], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, blend_3b=1.0_wp))
      lsf%screening_threshold = 1.0e-11_wp
      call lsf%update(mol, [3.0_wp])

      lo = 0.5_wp
      hi = 12.0_wp
      p(:, 1) = [lo, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(p, 0, v, error=failed)
      v_lo = v(1)
      p(:, 1) = [hi, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(p, 0, v, error=failed)
      rs = 0.0_wp
      call check(error, v_lo*v(1) < 0.0_wp, "sphere bracket")
      if (allocated(error)) return
      do it = 1, 200
         rs = 0.5_wp*(lo + hi)
         p(:, 1) = [rs, 0.0_wp, 0.0_wp]
         call lsf%evaluate_batch(p, 0, v, error=failed)
         if (v(1)*v_lo > 0.0_wp) then
            lo = rs
         else
            hi = rs
         end if
         if (hi - lo < 1.0e-14_wp) exit
      end do
      rs = 0.5_wp*(lo + hi)
   end subroutine sphere_lsf

   !> Unit vector along `v`
   !>
   !> @param[in] v Direction
   pure function unit_vector(v) result(u)
      !> Direction
      real(wp), intent(in) :: v(3)
      !> Unit vector
      real(wp) :: u(3)

      u = v/norm2(v)
   end function unit_vector

   !* ================================================================================= *!
   !*                  Octree surface coverage on an analytic level set                 *!
   !* ================================================================================= *!

   !> Probe callback: level set value and exact Lipschitz exclusion radius
   !>
   !> @param[in]  x       Evaluation point
   !> @param[out] lsf0    Level set value
   !> @param[out] radius  Surface-free radius around `x`
   !> @param[in]  context Sphere union definition
   subroutine sphere_union_probe(x, lsf0, radius, context)
      real(wp), intent(in) :: x(3)
      real(wp), intent(out) :: lsf0
      real(wp), intent(out) :: radius
      class(*), intent(in) :: context

      integer :: i
      real(wp) :: d

      lsf0 = huge(1.0_wp)
      radius = 0.0_wp

      select type (context)
      type is (sphere_union_context)
         do i = 1, context%nsphere
            d = norm2(x - context%center(:, i)) - context%radius(i)
            lsf0 = min(lsf0, d)
         end do
      class default
         error stop "sphere_union_probe: foreign context"
      end select

      radius = abs(lsf0)
   end subroutine sphere_union_probe

   !> Batched probe: the scalar sphere-union probe at every point
   !>
   !> @param[in]    points  Evaluation points (3, n)
   !> @param[out]   lsf0    Level set values
   !> @param[out]   radius  Surface-free radii
   !> @param[out]   reached Every point is evaluated
   !> @param[inout] context Sphere union definition
   !> @param[out]   error   Never set
   subroutine sphere_union_probe_batch(points, lsf0, radius, reached, context, error)
      real(wp), intent(in), contiguous :: points(:, :)
      real(wp), intent(out), contiguous :: lsf0(:)
      real(wp), intent(out), contiguous :: radius(:)
      logical, intent(out), contiguous :: reached(:)
      class(*), intent(inout) :: context
      type(moist_error_type), allocatable, intent(out) :: error

      integer :: k

      do k = 1, size(points, 2)
         call sphere_union_probe(points(:, k), lsf0(k), radius(k), context)
         reached(k) = .true.
      end do
   end subroutine sphere_union_probe_batch

   !> Batched probe of the sphere union that certifies no surface-free ball
   !>
   !> @param[in]    points  Evaluation points (3, n)
   !> @param[out]   lsf0    Level set values
   !> @param[out]   radius  Zero radii
   !> @param[out]   reached Every point is evaluated
   !> @param[inout] context Sphere union definition
   !> @param[out]   error   Never set
   subroutine sphere_union_probe_blind(points, lsf0, radius, reached, context, error)
      real(wp), intent(in), contiguous :: points(:, :)
      real(wp), intent(out), contiguous :: lsf0(:)
      real(wp), intent(out), contiguous :: radius(:)
      logical, intent(out), contiguous :: reached(:)
      class(*), intent(inout) :: context
      type(moist_error_type), allocatable, intent(out) :: error

      integer :: k

      call sphere_union_probe_batch(points, lsf0, radius, reached, context, error)
      do k = 1, size(points, 2)
         radius(k) = 0.0_wp
      end do
   end subroutine sphere_union_probe_blind

   !> Two equal spheres straddling the anchor: the closest-point projection is
   !> two-valued and the search must return both branches
   subroutine test_octree_two_sphere(error)
      type(error_type), allocatable, intent(out) :: error

      type(octree_seeder_type) :: octree
      type(sphere_union_context) :: ctx
      type(moist_error_type), allocatable :: solver_error
      real(wp) :: anchor(3), lsf0, radius
      integer :: i, n_left, n_right

      ! Spheres of radius 1 centered at x = -2 and x = +2. The anchor sits at
      ! the origin, equidistant from both surfaces at rho = 1
      ctx%nsphere = 2
      ctx%center(:, 1) = [-2.0_wp, 0.0_wp, 0.0_wp]
      ctx%center(:, 2) = [2.0_wp, 0.0_wp, 0.0_wp]
      ctx%radius = [1.0_wp, 1.0_wp]
      anchor = 0.0_wp

      call octree%init(seed_size=0.1_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree init: "//solver_error%message)
         return
      end if

      call sphere_union_probe(anchor, lsf0, radius, ctx)

      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=2.0_wp, &
                      rho2_slack=0.5_wp, probe=sphere_union_probe_batch, &
                      context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree run: "//solver_error%message)
         return
      end if

      call check(error, octree%n_points > 0, message="surface coverage produced no seeds")
      if (allocated(error)) return

      ! Both admissible surface patches have seeds near their surfaces
      n_left = 0
      n_right = 0
      do i = 1, octree%n_points
         if (octree%points(1, i) < 0.0_wp) n_left = n_left + 1
         if (octree%points(1, i) > 0.0_wp) n_right = n_right + 1
         call sphere_union_probe(octree%points(:, i), lsf0, radius, ctx)
         call check(error, radius <= sqrt(3.0_wp)*octree%seed_size/2.0_wp, &
                    message="seed lies farther than a leaf circumradius from the surface")
         if (allocated(error)) return
      end do
      call check(error, n_left > 0 .and. n_right > 0, &
                 message="one of the two surface patches has no seed")
   end subroutine test_octree_two_sphere

   !> A single sphere admits exactly one branch, and the search tightens its own
   !> radius down to that branch instead of certifying the ball it started from
   subroutine test_octree_single_sphere(error)
      type(error_type), allocatable, intent(out) :: error

      type(octree_seeder_type) :: octree
      type(sphere_union_context) :: ctx
      type(moist_error_type), allocatable :: solver_error
      real(wp) :: anchor(3), lsf0, radius

      ctx%nsphere = 1
      ctx%center(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      ctx%radius(1) = 2.0_wp
      anchor = [1.0_wp, 0.0_wp, 0.0_wp]

      call octree%init(seed_size=0.1_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree init: "//solver_error%message)
         return
      end if

      call sphere_union_probe(anchor, lsf0, radius, ctx)

      ! Start from a ball far wider than needed: the anchor is 1 Bohr inside a
      ! sphere of radius 2, so the only branch sits at rho = 1
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=5.0_wp, &
                      rho2_slack=0.25_wp, probe=sphere_union_probe_batch, &
                      context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree run: "//solver_error%message)
         return
      end if

      call check(error, octree%n_points > 0, message="surface coverage produced no seeds")
      if (allocated(error)) return

      ! rho_min = 1 with slack 0.25 gives sqrt(1.25) ~ 1.118, well below the 5.0
      ! the run started from: the sign-change bound did the tightening
      call check(error, octree%rho_max_final < 1.2_wp, &
                 message="search failed to tighten its own admissible radius")
   end subroutine test_octree_single_sphere

   !> The sign-change bound is refined by the probed center's own surface-free
   !> ball, which must leave it above the truth and should land close to it
   subroutine test_octree_sign_change_bound(error)
      type(error_type), allocatable, intent(out) :: error

      type(octree_seeder_type) :: octree
      type(sphere_union_context) :: ctx
      type(moist_error_type), allocatable :: solver_error
      real(wp) :: anchor(3), lsf0, radius, rho_true

      !> Squared-distance slack defining the admissible set
      real(wp), parameter :: slack = 0.25_wp

      ctx%nsphere = 1
      ctx%center(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      ctx%radius(1) = 1.0_wp
      anchor = [3.0_wp, 0.0_wp, 0.0_wp]

      call octree%init(seed_size=0.1_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree init: "//solver_error%message)
         return
      end if

      call sphere_union_probe(anchor, lsf0, radius, ctx)

      ! The anchor sits 2 Bohr outside the sphere, so rho_min is exactly 2 and
      ! no correct run may certify less than sqrt(2^2 + slack) = 2.0616. A run
      ! that bounds the crossing by the probed center alone, without taking off
      ! that center's surface-free radius, stops around 2.1035 instead
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=5.0_wp, &
                      rho2_slack=slack, probe=sphere_union_probe_batch, &
                      context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree run: "//solver_error%message)
         return
      end if

      call check(error, octree%n_points > 0, message="surface coverage produced no seeds")
      if (allocated(error)) return

      rho_true = sqrt(4.0_wp + slack)

      ! The epsilon absorbs a last-bit difference in how the bound is summed,
      ! not a real shortfall: the margin a working run keeps here is ~2e-3
      call check(error, octree%rho_max_final >= rho_true - 1.0e-10_wp, &
                 message="certified radius fell below the true admissible one")
      if (allocated(error)) return

      call check(error, octree%rho_max_final <= 1.01_wp*rho_true, &
                 message="sign-change bound is not being refined by the "// &
                 "exclusion radius")
   end subroutine test_octree_sign_change_bound

   !> A ball holding no surface is certified empty, with no seeds and without
   !> spending more than the root box
   subroutine test_octree_empty_ball(error)
      type(error_type), allocatable, intent(out) :: error

      type(octree_seeder_type) :: octree
      type(sphere_union_context) :: ctx
      type(moist_error_type), allocatable :: solver_error
      real(wp) :: anchor(3), lsf0, radius

      ctx%nsphere = 1
      ctx%center(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      ctx%radius(1) = 10.0_wp
      anchor = 0.0_wp

      call octree%init(seed_size=0.1_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree init: "//solver_error%message)
         return
      end if

      call sphere_union_probe(anchor, lsf0, radius, ctx)

      ! The surface is 10 Bohr away; a ball of radius 1 around the anchor holds
      ! none of it, and the root box alone proves that
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=1.0_wp, &
                      rho2_slack=0.25_wp, probe=sphere_union_probe_batch, &
                      context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, "octree run: "//solver_error%message)
         return
      end if

      call check(error, octree%n_points, 0, &
                 message="empty ball produced seeds")
      if (allocated(error)) return
      call check(error, octree%n_boxes_visited, 1, &
                 message="empty ball took more than the root box to certify")
   end subroutine test_octree_empty_ball

   !> Cover an analytically known spherical cap throughout its admissible window
   subroutine test_octree_surface_patch_coverage(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Search under test
      type(octree_seeder_type) :: octree
      !> Analytic sphere
      type(sphere_union_context) :: ctx
      !> Search error
      type(moist_error_type), allocatable :: solver_error
      !> Search anchor and level-set probe
      real(wp) :: anchor(3), lsf0, radius
      !> Exact cap samples and distance to the closest seed box
      real(wp) :: point(3), theta, phi, distance
      !> Sample and seed indices
      integer :: itheta, iphi, i

      ctx%nsphere = 1
      ctx%center(:, 1) = 0.0_wp
      ctx%radius(1) = 1.0_wp
      anchor = [2.0_wp, 0.0_wp, 0.0_wp]
      call sphere_union_probe(anchor, lsf0, radius, ctx)
      call octree%init(seed_size=0.1_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, solver_error%message)
         return
      end if
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=2.0_wp, &
                      rho2_slack=0.5_wp, probe=sphere_union_probe_batch, &
                      context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, solver_error%message)
         return
      end if

      ! rho^2 = 5 - 4*cos(theta); the admissible cap has cos(theta) >= 7/8
      do itheta = 0, 8
         theta = acos(0.875_wp)*real(itheta, wp)/8.0_wp
         do iphi = 0, 23
            phi = 2.0_wp*acos(-1.0_wp)*real(iphi, wp)/24.0_wp
            point = [cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi)]
            distance = huge(1.0_wp)
            do i = 1, octree%n_points
               distance = min(distance, maxval(abs(point - octree%points(:, i))))
            end do
            call check(error, distance <= 0.5_wp*octree%seed_size + 1.0e-12_wp, &
                       message="an admissible surface point has no seed box")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_octree_surface_patch_coverage

   !> Keep analytic sphere stationary points while discarding transverse boxes
   subroutine test_octree_region_certificate(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Reused tree
      type(octree_seeder_type) :: tree
      !> Unit sphere at the origin
      type(sphere_union_context) :: context
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Distance origin and probe values
      real(wp) :: anchor(3), value, radius, distance
      !> Unpruned count and nearest/farthest sphere point indices
      integer :: unpruned, sign_index, i

      anchor = [2.0_wp, 0.0_wp, 0.0_wp]
      call sphere_union_probe(anchor, value, radius, context)
      call tree%init(seed_size=0.1_wp, error=native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call tree%run(anchor, value, 4.0_wp, 20.0_wp, sphere_union_probe_batch, context, native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      unpruned = tree%n_points
      call tree%run(anchor, value, 4.0_wp, 20.0_wp, sphere_union_probe_batch, context, native_error, &
                    exclude_region=sphere_stationarity_excluded)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, tree%n_regions_excluded > 0 .and. 2*tree%n_points < unpruned, &
                 message="Optional certificate did not discard transverse boxes")
      if (allocated(error)) return
      do sign_index = -1, 1, 2
         distance = huge(1.0_wp)
         do i = 1, tree%n_points
            distance = min(distance, maxval(abs(tree%points(:, i) - [real(sign_index, wp), 0.0_wp, 0.0_wp])))
         end do
         call check(error, distance <= 0.5_wp*tree%seed_size, &
                    message="Optional certificate removed an analytic stationary point")
         if (allocated(error)) return
      end do
      call tree%run(anchor, value, 4.0_wp, 20.0_wp, sphere_union_probe_batch, context, native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, tree%n_points, unpruned, message="Reference search changed after a pruned run")
      if (allocated(error)) return
      call check(error, tree%n_regions_excluded, 0, message="Region exclusion count was not reset")
   end subroutine test_octree_region_certificate

   !> With no slack, an anchor within the sphere-trace tolerance of the surface
   !> is its own single seed, and the last search's tallies are cleared
   !>
   !> - the trace stops at the anchor, so the starting radius is zero; a
   !>   nonzero value at the anchor skips the exact-surface shortcut
   !>
   !> @param[out] error Missing seed or stale tally
   subroutine test_octree_trace_on_surface(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(octree_seeder_type) :: octree
      type(seed_request_type) :: request
      type(moist_error_type), allocatable :: failed
      !> Surface radius, probe point and value
      real(wp) :: rs, p(3, 1), v(1)

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      call octree%init(failed, seed_size=0.2_wp)
      if (allocated(failed)) then
         call test_failed(error, "octree init: "//failed%message)
         return
      end if

      ! A regular search first, so the tallies hold something
      p(:, 1) = [rs + 1.0_wp, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(p, 0, v, error=failed)
      request%anchor = p(:, 1)
      request%lsf0_anchor = v(1)
      request%rho2_slack = 0.25_wp
      if (.not. allocated(failed)) call octree%seed(lsf, request, failed)
      if (allocated(failed)) then
         call test_failed(error, "octree seed: "//failed%message)
         return
      end if
      call check(error, octree%n_boxes_visited > 0, "the regular search probes boxes")
      if (allocated(error)) return

      p(:, 1) = [rs + 1.0e-9_wp, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(p, 0, v, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "anchor value: "//failed%message)
         return
      end if
      call check(error, v(1) /= 0.0_wp .and. abs(v(1)) < 1.0e-7_wp, "anchor next to, not on, the surface")
      if (allocated(error)) return
      request%anchor = p(:, 1)
      request%lsf0_anchor = v(1)
      request%rho2_slack = 0.0_wp
      request%rho_start = -1.0_wp
      call octree%seed(lsf, request, failed)
      if (allocated(failed)) then
         call test_failed(error, "octree seed: "//failed%message)
         return
      end if
      call check(error, octree%n_points, 1, "one seed")
      if (allocated(error)) return
      call check(error, all(octree%points(:, 1) == request%anchor), "the seed is the anchor")
      if (allocated(error)) return
      call check(error, octree%n_boxes_visited == 0 .and. octree%rho_max_final == 0.0_wp, &
                 "the tallies of the previous search are cleared")
   end subroutine test_octree_trace_on_surface

   !> A verified foot beyond the fallback search radius bounds the search as it is
   !>
   !> @param[out] error Missing or misplaced seeds
   subroutine test_octree_far_foot_radius(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_octree_far_surface(error, from_foot=.true.)
   end subroutine test_octree_far_foot_radius

   !> A surface beyond the fallback search radius, reached by the trace, bounds the
   !> search as it is
   !>
   !> @param[out] error Missing or misplaced seeds
   subroutine test_octree_far_traced_radius(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_octree_far_surface(error, from_foot=.false.)
   end subroutine test_octree_far_traced_radius

   !> Seeds at a single sphere whose only surface lies beyond the fallback radius
   !>
   !> - the anchor sits `distance` outside the sphere on the x axis, so the
   !>   analytic foot is the radial one and the whole window lies beyond the
   !>   request's `rho_fallback`
   !> - the request carries the radius of that foot, `sqrt(distance^2 + slack)`,
   !>   or no radius, so the seeder traces the surface itself
   !> - every seed lies near the foot and at least one exists
   !>
   !> @param[out] error     Missing or misplaced seeds
   !> @param[in]  from_foot Whether the request carries the radius of the foot
   subroutine check_octree_far_surface(error, from_foot)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Whether the request carries the radius of the foot
      logical, intent(in) :: from_foot
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(octree_seeder_type) :: octree
      type(seed_request_type) :: request
      type(moist_error_type), allocatable :: failed
      !> Surface radius, anchor and its value, analytic foot
      real(wp) :: rs, p(3, 1), v(1), foot(3)
      !> Seed index
      integer :: i
      !> Anchor distance from the surface (Bohr)
      real(wp), parameter :: distance = 18.0_wp
      !> Squared-distance window (Bohr^2)
      real(wp), parameter :: slack = 0.25_wp
      !> Largest seed distance from the foot (Bohr), far below `distance`
      real(wp), parameter :: near = 1.0_wp

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      call octree%init(failed, seed_size=0.2_wp)
      if (allocated(failed)) then
         call test_failed(error, "octree init: "//failed%message)
         return
      end if

      p(:, 1) = [rs + distance, 0.0_wp, 0.0_wp]
      foot = [rs, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(p, 0, v, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "anchor value: "//failed%message)
         return
      end if
      call check(error, .not. lsf%unreached(1), "fixture anchor lies within the screening reach")
      if (allocated(error)) return

      request%anchor = p(:, 1)
      request%lsf0_anchor = v(1)
      request%rho2_slack = slack
      request%rho_start = -1.0_wp
      if (from_foot) request%rho_start = sqrt(distance**2 + slack)
      call check(error, distance > request%rho_fallback, "fixture surface lies beyond the fallback radius")
      if (allocated(error)) return

      call octree%seed(lsf, request, failed)
      if (allocated(failed)) then
         call test_failed(error, "octree seed: "//failed%message)
         return
      end if
      call check(error, octree%n_points > 0, "the far surface produced no seeds")
      if (allocated(error)) return
      do i = 1, octree%n_points
         call check(error, norm2(octree%points(:, i) - foot) <= near, "a seed lies far from the foot")
         if (allocated(error)) return
      end do
   end subroutine check_octree_far_surface

   !> Seeds sit on sign-change segments between face-neighbour leaf centers
   !>
   !> - a segment runs along one lattice axis, so a seed shares two coordinates
   !>   with the leaf-center lattice; a leaf-center seed shares all three
   !> - the seed lies between two centers on opposite sides of the surface, so
   !>   the exact distance of the sphere puts it within one leaf edge of it
   !>
   !> @param[out] error Seed off a sign-change segment
   subroutine test_octree_sign_change_seeds(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Search under test
      type(octree_seeder_type) :: octree
      !> Unit sphere at the origin
      type(sphere_union_context) :: ctx
      !> Search error
      type(moist_error_type), allocatable :: solver_error
      !> Anchor, probe value and radius, leaf edge, lattice coordinate
      real(wp) :: anchor(3), lsf0, radius, edge, q
      !> Seed, coordinate and lattice-coordinate counters
      integer :: i, c, n_lattice

      ctx%nsphere = 1
      ctx%center(:, 1) = 0.0_wp
      ctx%radius(1) = 1.0_wp
      anchor = [2.0_wp, 0.0_wp, 0.0_wp]
      call sphere_union_probe(anchor, lsf0, radius, ctx)
      call octree%init(seed_size=0.1_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, solver_error%message)
         return
      end if
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=2.0_wp, rho2_slack=0.5_wp, &
                      probe=sphere_union_probe_batch, context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, solver_error%message)
         return
      end if
      call check(error, octree%n_points > 0, message="the sphere cap produced no seeds")
      if (allocated(error)) return

      edge = 2.0_wp*octree%root_half/real(2**octree%dbg_leaf_depth, wp)
      do i = 1, octree%n_points
         call sphere_union_probe(octree%points(:, i), lsf0, radius, ctx)
         call check(error, abs(lsf0) <= edge, message="a seed lies farther than a leaf edge from the surface")
         if (allocated(error)) return
         ! Center offsets are `edge/2 (2 lat + 1 - 2^depth)` for integer `lat`
         n_lattice = 0
         do c = 1, 3
            q = 0.5_wp*((octree%points(c, i) - anchor(c))/(0.5_wp*edge) + real(2**octree%dbg_leaf_depth - 1, wp))
            if (abs(q - real(nint(q), wp)) < 1.0e-9_wp) n_lattice = n_lattice + 1
         end do
         call check(error, n_lattice, 2, message="a seed is not on a segment between face-neighbour centers")
         if (allocated(error)) return
      end do
   end subroutine test_octree_sign_change_seeds

   !> A ball without surface that nothing certifies yields leaves but no seeds
   !>
   !> - the blind probe frees no box, so every leaf survives, has no sign change
   !>   and is refined without finding one
   !>
   !> @param[out] error Seed from a sign-free leaf
   subroutine test_octree_sign_free_region(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Search under test
      type(octree_seeder_type) :: octree
      !> Sphere far larger than the search ball
      type(sphere_union_context) :: ctx
      !> Search error
      type(moist_error_type), allocatable :: solver_error
      !> Anchor and probe
      real(wp) :: anchor(3), lsf0, radius
      !> Surviving leaves
      integer :: n_leaves

      ctx%nsphere = 1
      ctx%center(:, 1) = 0.0_wp
      ctx%radius(1) = 10.0_wp
      anchor = 0.0_wp
      call sphere_union_probe(anchor, lsf0, radius, ctx)
      call octree%init(seed_size=0.5_wp, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, solver_error%message)
         return
      end if
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=1.0_wp, rho2_slack=0.25_wp, &
                      probe=sphere_union_probe_blind, context=ctx, error=solver_error)
      if (allocated(solver_error)) then
         call test_failed(error, solver_error%message)
         return
      end if
      n_leaves = octree%dbg_leaves(octree%dbg_leaf_depth)
      call check(error, n_leaves > 0, message="the uncertified ball kept no leaf")
      if (allocated(error)) return
      call check(error, octree%dbg_refined, n_leaves, message="a leaf without evidence was not refined")
      if (allocated(error)) return
      call check(error, octree%n_points, 0, message="a sign-free region produced seeds")
   end subroutine test_octree_sign_free_region

   !> A feature thinner than a leaf, between the leaf centers, is found by the child probes
   !>
   !> - a sphere of radius 0.05 around a child center of the leaf at
   !>   (0.25, 0.25, 0.25), edge 0.5: no leaf center is inside it, so no leaf
   !>   has evidence, and the crossing on the radial parent-child segment is
   !>   exact for the distance function
   !> - with the refinement disabled the search finds nothing
   !>
   !> @param[out] error Missed or spurious feature
   subroutine test_octree_refinement_thin_feature(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Search under test
      type(octree_seeder_type) :: octree
      !> Small sphere off every leaf center
      type(sphere_union_context) :: ctx
      !> Search error
      type(moist_error_type), allocatable :: solver_error
      !> Anchor and probe
      real(wp) :: anchor(3), lsf0, radius, closest
      !> Seed index and refinement setting
      integer :: i, levels

      ctx%nsphere = 1
      ctx%center(:, 1) = [0.375_wp, 0.375_wp, 0.375_wp]
      ctx%radius(1) = 0.05_wp
      anchor = 0.0_wp
      call sphere_union_probe(anchor, lsf0, radius, ctx)
      do levels = 1, 0, -1
         call octree%init(seed_size=0.5_wp, refine_levels=levels, error=solver_error)
         if (allocated(solver_error)) then
            call test_failed(error, solver_error%message)
            return
         end if
         call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=1.0_wp, rho2_slack=0.25_wp, &
                         probe=sphere_union_probe_batch, context=ctx, error=solver_error)
         if (allocated(solver_error)) then
            call test_failed(error, solver_error%message)
            return
         end if
         call check(error, octree%dbg_evidence, 0, message="a leaf center resolves the thin feature")
         if (allocated(error)) return
         if (levels == 0) then
            call check(error, octree%n_points, 0, message="the feature was found without refinement")
            return
         end if
         call check(error, octree%n_points, 1, message="refinement did not seed the thin feature once")
         if (allocated(error)) return
         closest = huge(1.0_wp)
         do i = 1, octree%n_points
            call sphere_union_probe(octree%points(:, i), lsf0, radius, ctx)
            closest = min(closest, abs(lsf0))
         end do
         call check(error, closest < 1.0e-12_wp, message="the child seed is off the feature surface")
         if (allocated(error)) return
      end do
   end subroutine test_octree_refinement_thin_feature

   !> Exact stationary-axis exclusion for the analytic sphere fixture
   !>
   !> @param[in] center     Box midpoint (3)
   !> @param[in] half_width Cube half-width
   !> @param[in] context    Single sphere centered on the anchor's x-axis
   function sphere_stationarity_excluded(center, half_width, context) result(excluded)
      !> Box midpoint
      real(wp), intent(in) :: center(3)
      !> Cube half-width
      real(wp), intent(in) :: half_width
      !> Analytic sphere fixture
      class(*), intent(in) :: context
      !> Box cannot intersect the stationary x-axis
      logical :: excluded

      excluded = .false.
      select type (context)
      type is (sphere_union_context)
         excluded = any(abs(center(2:3) - context%center(2:3, 1)) > half_width)
      end select
   end function sphere_stationarity_excluded

   !> Fail an expected-failure test only on the targeted refusal
   !>
   !> - test-drive inverts `should_fail` tests, so any other outcome, including
   !>   a different error, reports the test as failed
   !>
   !> @param[out] error    Test failure, set only on the expected error
   !> @param[in]  err      Library error, possibly unallocated
   !> @param[in]  expected The exact refusal text
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

   !> Run the octree on a 2 Bohr sphere from an anchor 1 Bohr off its center
   !>
   !> @param[inout] octree    Initialized seeder
   !> @param[in]    rho_max   Starting radius
   !> @param[in]    rho2_slack Squared-distance window
   !> @param[out]   failed    Native refusal
   subroutine run_sphere_octree(octree, rho_max, rho2_slack, failed)
      !> Initialized seeder
      type(octree_seeder_type), intent(inout) :: octree
      !> Starting radius and window
      real(wp), intent(in) :: rho_max, rho2_slack
      !> Native refusal
      type(moist_error_type), allocatable, intent(out) :: failed
      type(sphere_union_context) :: ctx
      real(wp) :: anchor(3), lsf0, radius

      ctx%nsphere = 1
      ctx%center(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      ctx%radius(1) = 2.0_wp
      anchor = [1.0_wp, 0.0_wp, 0.0_wp]
      call sphere_union_probe(anchor, lsf0, radius, ctx)
      call octree%run(anchor=anchor, lsf0_anchor=lsf0, rho_max=rho_max, rho2_slack=rho2_slack, &
                      probe=sphere_union_probe_batch, context=ctx, error=failed)
   end subroutine run_sphere_octree

   !> A nonpositive seed size is refused
   subroutine test_bad_octree_seed_size(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, seed_size=0.0_wp)
      call expect_error(error, failed, "Octree branch search: seed_size must be positive")
   end subroutine test_bad_octree_seed_size

   !> A nonpositive box or survivor budget is refused
   subroutine test_bad_octree_budgets(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, max_boxes=0)
      call expect_error(error, failed, "Octree branch search: box and survivor budgets must be positive")
   end subroutine test_bad_octree_budgets

   !> A depth limit outside the addressable lattice is refused
   subroutine test_bad_octree_max_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, max_depth=0)
      call expect_error(error, failed, "Octree branch search: max_depth out of range")
   end subroutine test_bad_octree_max_depth

   !> A nonpositive starting radius is refused
   subroutine test_bad_octree_rho_max(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, seed_size=0.1_wp)
      if (allocated(failed)) return
      call run_sphere_octree(octree, 0.0_wp, 0.25_wp, failed)
      call expect_error(error, failed, "Octree branch search: rho_max must be positive")
   end subroutine test_bad_octree_rho_max

   !> A negative squared-distance slack would poison the radius cap with a NaN,
   !> after which every comparison against it is false and nothing is bounded
   subroutine test_bad_octree_negative_slack(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, seed_size=0.1_wp)
      if (allocated(failed)) return
      call run_sphere_octree(octree, 5.0_wp, -0.25_wp, failed)
      call expect_error(error, failed, "Octree branch search: rho2_slack must not be negative")
   end subroutine test_bad_octree_negative_slack

   !> A depth cap that stops short of the requested seed size must be reported
   !>
   !> - continuing with coarser leaves would merge minima the configured
   !>   resolution separates while still claiming a completed enumeration; a
   !>   5 Bohr ball needs depth 7 to reach 0.1 Bohr leaves, three is not enough
   subroutine test_bad_octree_depth_cap(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, seed_size=0.1_wp, max_depth=3)
      if (allocated(failed)) return
      call run_sphere_octree(octree, 5.0_wp, 0.25_wp, failed)
      call expect_error(error, failed, "Octree branch search: max_depth (3) bottoms out at a leaf edge "// &
                        "of 1.2500 Bohr, coarser than the requested seed_size of .1000 Bohr. Raise "// &
                        "max_depth to at least 7 or loosen seed_size.")
   end subroutine test_bad_octree_depth_cap

   !> More surviving leaves than the survivor budget is refused
   subroutine test_bad_octree_survivors(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, seed_size=4.0_wp, max_survivors=1)
      if (allocated(failed)) return
      call run_sphere_octree(octree, 5.0_wp, 0.25_wp, failed)
      call expect_error(error, failed, "Octree branch search: survivor budget exhausted; "// &
                        "seed_size is too small for the admissible radius.")
   end subroutine test_bad_octree_survivors

   !> A next frontier beyond the box budget is refused while it is generated
   !>
   !> - the root contains the anchor and splits into eight children that all
   !>   pass the cap, nine boxes against a budget of eight
   subroutine test_bad_octree_box_budget(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, seed_size=0.1_wp, max_boxes=8)
      if (allocated(failed)) return
      call run_sphere_octree(octree, 5.0_wp, 0.25_wp, failed)
      call expect_error(error, failed, "Octree branch search: box budget exhausted before the "// &
                        "search ball was certified. Either the level set supplies no exclusion "// &
                        "radius (pruning nothing), or seed_size is too small for the admissible radius.")
   end subroutine test_bad_octree_box_budget

   !> Seeding with an uninitialized octree is refused
   subroutine test_bad_octree_seed_before_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_request_type) :: request
      type(moist_error_type), allocatable :: failed

      request%anchor = [1.0_wp, 0.0_wp, 0.0_wp]
      call octree%seed(lsf, request, failed)
      call expect_error(error, failed, "Octree branch search: run before init")
   end subroutine test_bad_octree_seed_before_init

   !> Running an uninitialized octree is refused
   subroutine test_bad_octree_run_before_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call run_sphere_octree(octree, 5.0_wp, 0.25_wp, failed)
      call expect_error(error, failed, "Octree branch search: run before init")
   end subroutine test_bad_octree_run_before_init

   !> A negative refinement depth is refused
   subroutine test_bad_octree_refine_levels(error)
      type(error_type), allocatable, intent(out) :: error
      type(octree_seeder_type) :: octree
      type(moist_error_type), allocatable :: failed

      call octree%init(failed, refine_levels=-1)
      call expect_error(error, failed, "Octree branch search: refine_levels must not be negative")
   end subroutine test_bad_octree_refine_levels

   !> Face-neighbour probes are not charged to the box budget
   !>
   !> - the budget is the tree descent of the same search, so every depth fits
   !>   exactly and the face neighbours lie beyond it
   !> - the seeds equal those of the unbudgeted search
   !>
   !> @param[out] error Refusal, or seeds differing from the unbudgeted search
   subroutine test_octree_neighbours_outside_budget(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Search under test
      type(octree_seeder_type) :: octree
      !> Search error
      type(moist_error_type), allocatable :: failed
      !> Seeds of the unbudgeted search
      real(wp), allocatable :: ref_points(:, :)
      !> Tree-descent boxes of the unbudgeted search
      integer :: descent

      call octree%init(failed, seed_size=0.5_wp)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call run_sphere_octree(octree, 5.0_wp, 0.25_wp, failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call check(error, octree%dbg_neighbours > 0, "the fixture probes no face neighbour")
      if (allocated(error)) return
      call check(error, octree%n_boxes_visited, sum(octree%dbg_probed), &
                 message="the budgeted count is not the tree descent")
      if (allocated(error)) return
      descent = octree%n_boxes_visited
      ref_points = octree%points(:, 1:octree%n_points)

      call octree%init(failed, seed_size=0.5_wp, max_boxes=descent)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call run_sphere_octree(octree, 5.0_wp, 0.25_wp, failed)
      if (allocated(failed)) then
         call test_failed(error, "budget of the tree descent: "//failed%message)
         return
      end if
      call check(error, octree%n_points, size(ref_points, 2), &
                 message="the budgeted search emits a different seed count")
      if (allocated(error)) return
      call check(error, all(octree%points(:, 1:octree%n_points) == ref_points), &
                 "the budgeted search emits different seeds")
   end subroutine test_octree_neighbours_outside_budget

   !> Refinement probes are not charged to the box budget
   !>
   !> - a blind search of a sign-free ball refines every leaf; the budget is
   !>   its tree descent, so the face neighbours and children lie beyond it
   !> - the seeds and the refinement equal those of the unbudgeted search
   !>
   !> @param[out] error Refusal, or a result differing from the unbudgeted search
   subroutine test_octree_refinement_outside_budget(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Search under test
      type(octree_seeder_type) :: octree
      !> Sphere far larger than the search ball
      type(sphere_union_context) :: ctx
      !> Search error
      type(moist_error_type), allocatable :: failed
      !> Anchor and probe
      real(wp) :: anchor(3), lsf0, radius
      !> Seeds of the unbudgeted search
      real(wp), allocatable :: ref_points(:, :)
      !> Tree-descent boxes and children probed by the unbudgeted search
      integer :: descent, children

      ctx%nsphere = 1
      ctx%center(:, 1) = 0.0_wp
      ctx%radius(1) = 10.0_wp
      anchor = 0.0_wp
      call sphere_union_probe(anchor, lsf0, radius, ctx)
      call octree%init(failed, seed_size=0.5_wp)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call octree%run(anchor, lsf0, 1.0_wp, 0.25_wp, sphere_union_probe_blind, ctx, failed)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call check(error, octree%dbg_children > 0, "the fixture probes no refinement child")
      if (allocated(error)) return
      descent = octree%n_boxes_visited
      children = octree%dbg_children
      ref_points = octree%points(:, 1:octree%n_points)

      call octree%init(failed, seed_size=0.5_wp, max_boxes=descent)
      if (allocated(failed)) then
         call test_failed(error, failed%message)
         return
      end if
      call octree%run(anchor, lsf0, 1.0_wp, 0.25_wp, sphere_union_probe_blind, ctx, failed)
      if (allocated(failed)) then
         call test_failed(error, "budget of the tree descent: "//failed%message)
         return
      end if
      call check(error, octree%dbg_children, children, &
                 message="the budgeted search probes a different number of children")
      if (allocated(error)) return
      call check(error, octree%n_points, size(ref_points, 2), &
                 message="the budgeted search emits a different seed count")
      if (allocated(error)) return
      call check(error, all(octree%points(:, 1:octree%n_points) == ref_points), &
                 "the budgeted search emits different seeds")
   end subroutine test_octree_refinement_outside_budget

   !> The worker refines evidence-free leaves in the fine octree only
   !>
   !> - the octree each level's worker builds (the fallback at level 1, the
   !>   seeder at levels 3, 4) runs a blind search of a sign-free ball, so every
   !>   surviving leaf lacks evidence
   !> - the coarse octree of levels 1, 3 keeps such leaves but probes no child;
   !>   the fine octree of level 4 probes their children
   !>
   !> @param[out] error Refinement differing from the level's rule
   subroutine test_worker_octree_refinement(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Projection settings
      type(moist_cavity_drop_parameters_type) :: param
      !> One-sphere level set the worker is built on
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Worker under test
      type(projection_worker_type) :: worker
      !> Molecular geometry
      type(structure_type) :: mol
      !> Sphere far larger than the search ball
      type(sphere_union_context) :: ctx
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Atom position, anchor and probe
      real(wp) :: positions(3, 1), anchor(3), lsf0, radius
      !> Projection level and its index in `levels`
      integer :: il, level
      !> Every projection level
      integer, parameter :: levels(3) = [1, 3, 4]
      !> Level label for the messages
      character(len=8) :: label

      positions = 0.0_wp
      call new(mol, [1], positions)
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      call lsf%update(mol, [1.0_wp])
      ctx%nsphere = 1
      ctx%center(:, 1) = 0.0_wp
      ctx%radius(1) = 10.0_wp
      anchor = 0.0_wp
      call sphere_union_probe(anchor, lsf0, radius, ctx)

      do il = 1, size(levels)
         level = levels(il)
         write (label, "(a,i0)") "level ", level
         param%proj_level = level
         call param%compute_derived(native_error)
         if (.not. allocated(native_error)) call worker%init(param, lsf, positions, native_error)
         if (allocated(native_error)) then
            call test_failed(error, trim(label)//": "//native_error%message)
            return
         end if
         if (level == 1) then
            call check(error, allocated(worker%fallback), message=trim(label)//": no fallback octree")
            if (allocated(error)) return
            call run_blind(worker%fallback)
         else
            call check(error, allocated(worker%seeder), message=trim(label)//": no seeder")
            if (allocated(error)) return
            call run_blind(worker%seeder)
         end if
         if (allocated(error)) return
      end do

   contains

      !> Run the blind search with `seeder` and check its refinement
      !>
      !> @param[inout] seeder Octree built by the worker
      subroutine run_blind(seeder)
         class(projection_seeder_type), intent(inout) :: seeder

         select type (seeder)
         type is (octree_seeder_type)
            call seeder%run(anchor, lsf0, 1.0_wp, 0.25_wp, sphere_union_probe_blind, ctx, native_error)
            if (allocated(native_error)) then
               call test_failed(error, trim(label)//": "//native_error%message)
               return
            end if
            call check(error, seeder%dbg_refined > 0, &
                       message=trim(label)//": the fixture keeps no leaf without evidence")
            if (allocated(error)) return
            if (level <= 3) then
               call check(error, seeder%dbg_children, 0, &
                          message=trim(label)//": the coarse octree probed refinement children")
            else
               call check(error, seeder%dbg_children > 0, &
                          message=trim(label)//": the fine octree probed no refinement child")
            end if
         class default
            call test_failed(error, trim(label)//": the worker's search is not an octree")
         end select
      end subroutine run_blind

   end subroutine test_worker_octree_refinement

   !> A batch with a nonpositive flush target is refused
   subroutine test_bad_batch_target_size(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_batch_type) :: batch
      type(moist_error_type), allocatable :: failed

      call batch%init(0, failed)
      call expect_error(error, failed, "Projection batch: target size must be positive")
   end subroutine test_bad_batch_target_size

   !> An anchor beyond the level set reach is refused, its anchor jet cannot be recorded
   subroutine test_bad_newton_unreached_anchor(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(projection_batch_type) :: batch
      type(newton_config_type) :: cfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3)
      integer :: icol

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      anchor = [80.0_wp, 0.0_wp, 0.0_wp]
      call batch%init(8, failed)
      call batch%begin_run(7, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, cfg, scratch, failed)
      call expect_error(error, failed, &
                        "Projection Newton: anchor 7 lies beyond the level set reach, no atom screens it")
   end subroutine test_bad_newton_unreached_anchor

end module test_cavity_drop_projection
