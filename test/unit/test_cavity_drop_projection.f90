!> Regression tests for the batched DROP projection: branch discovery in
!> connected surface patches, cluster branch sets against an exhaustive
!> surface-leaf reference, and the batch, Newton, polish, seeder and finalize
!> building blocks
module test_cavity_drop_projection
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, moist_error_type => error_type
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   use mctc_io, only: structure_type, new
   use mctc_io_constants, only: pi
   use mctc_io_convert, only: aatoau
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_radii, only: default_cpcm_radii
   use moist_math_linalg, only: setup_tangent_frame
   use moist_math_packed, only: unpack_sym_matrix
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use test_helpers, only: expect_error
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, drop_solve_parameters_type, &
                                           drop_polish_parameters_type, drop_finalize_parameters_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_lsf_sses, only: moist_cavity_drop_lsf_sses_type
   use moist_cavity_drop_projection_batch, only: projection_batch_type, &
                                                point_pending, point_converged, point_accepted, &
                                                point_captured, point_failed
   use moist_cavity_drop_projection_newton, only: newton_scratch_type, newton_solve_batch
   use moist_cavity_drop_projection_polish, only: polish_batch, retract_to_surface
   use moist_cavity_drop_projection_finalize, only: dedupe_candidates, canonical_branch_order, &
                                                   finalize_anchor_type, &
                                                   finalize_scratch_type, projection_rows_type, &
                                                   finalize_anchor
   use moist_math_smoothstep, only: moist_math_smoothstep_sigmoid_type, new_smoothstep_sigmoid
   use moist_cavity_drop_branching, only: branch_weight_type
   use moist_cavity_drop_projection_seed, only: projection_seeder_type, seeder_scratch_type, seed_request_type, &
                                                seed_build_type
   use moist_cavity_drop_projection_parameters, only: drop_seed_parameters_type
   use moist_cavity_drop_projection_worker, only: projection_worker_type
   implicit none(type, external)
   private

   public :: collect_cavity_drop_projection

   !> Level set evaluator for the unpruned surface-coverage reference
   !> SvdW level set that counts the points it evaluates
   type, extends(moist_cavity_drop_lsf_svdw_type) :: counting_svdw_type
      !> Points evaluated, every order
      integer :: n_points = 0
      !> Point whose evaluations with a Hessian are counted
      real(wp) :: probe(3) = huge(1.0_wp)
      !> Evaluations of `probe` at order 2 or higher
      integer :: n_probe_hess = 0
   contains
      !> Count, then evaluate as the SvdW level set
      procedure :: evaluate_batch => counting_evaluate_batch
   end type counting_svdw_type

   !> SvdW level set whose value at one point is NaN, gradient and reach kept:
   !> a host density that is not finite
   type, extends(moist_cavity_drop_lsf_svdw_type) :: nan_value_svdw_type
      !> Point whose value is replaced by NaN; non-finite points are answered alike
      real(wp) :: poisoned(3) = huge(1.0_wp)
   contains
      !> Evaluate as the SvdW level set, then replace the value at `poisoned`
      procedure :: evaluate_batch => nan_value_evaluate_batch
   end type nan_value_svdw_type

contains

   !> Collect branch discovery regressions
   !>
   !> @param[out] testsuite Projection tests
   subroutine collect_cavity_drop_projection(testsuite)
      !> Registered tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("h2_connected_minima", test_connected_minima), &
         & new_unittest("h2_shallow_minimum", test_shallow_minimum), &
         & new_unittest("h2_search_window", test_search_window), &
         & new_unittest("h2_rigid_motion", test_rigid_motion), &
         & new_unittest("h2_zero_window_on_surface", test_zero_window), &
         & new_unittest("cluster_dimer", test_dimer_branching), &
         & new_unittest("cluster_cross", test_branching_xyz_cross), &
         & new_unittest("cluster_octahedron", test_octahedral), &
         & new_unittest("cluster_cube", test_cube), &
         & new_unittest("cluster_pentagon", test_pentagonal), &
         & new_unittest("cluster_tetrahedron", test_tetrahedral), &
         & new_unittest("brute_force_broad_window", test_brute_force_broad_window), &
         & new_unittest("stationary_box_bounds", test_stationary_box_bounds), &
         & new_unittest("batch_flush_accumulator", test_batch_accumulator), &
         & new_unittest("newton_sphere_point", test_newton_sphere_point), &
         & new_unittest("newton_start_multiplier", test_newton_start_multiplier), &
         & new_unittest("newton_far_seed_unreached", test_newton_far_seed_unreached), &
         & new_unittest("newton_far_point", test_newton_far_point), &
         & new_unittest("newton_singular_system", test_newton_singular_system), &
         & new_unittest("newton_non_finite_value_fails_alone", test_newton_non_finite_value), &
         & new_unittest("polish_escapes_far_pole", test_polish_far_pole), &
         & new_unittest("polish_fails_non_finite_jet", test_polish_non_finite_jet), &
         & new_unittest("polish_one_sided_saddle_repeats_captured", test_polish_one_sided_saddle), &
         & new_unittest("multistart_capture_stops_duplicates", test_capture_stops_duplicates), &
         & new_unittest("multistart_capture_keeps_twin_points", test_capture_keeps_twin_points), &
         & new_unittest("multistart_capture_only_predictions", test_capture_only_predictions), &
         & new_unittest("finalize_canonical_order", test_finalize_canonical_order), &
         & new_unittest("finalize_canonical_order_chain", test_finalize_canonical_order_chain), &
         & new_unittest("finalize_sphere_anchor", test_finalize_sphere_anchor), &
         & new_unittest("finalize_branch_weight_fold", test_finalize_branch_fold), &
         & new_unittest("rmesh_build_thread_independent", test_rmesh_thread_independent), &
         & new_unittest("rmesh_level3_filter_and_level4_unfiltered", test_rmesh_level_filters), &
         & new_unittest("rmesh_shell_query_partitions_ball", test_rmesh_shell_query), &
         & new_unittest("rmesh_sses_build_thread_independent", test_rmesh_sses_thread_independent), &
         & new_unittest("sses_projection_thread_independent", test_sses_projection_thread_independent), &
         & new_unittest("level2_selective_rmesh_branches", test_level2_selective_rmesh), &
         & new_unittest("level2_seed_filters", test_level2_seed_filters), &
         & new_unittest("bad_rmesh_sample_size", test_bad_rmesh_sample_size, should_fail=.true.), &
         & new_unittest("bad_rmesh_refine_levels", test_bad_rmesh_refine_levels, should_fail=.true.), &
         & new_unittest("bad_rmesh_build_before_init", test_bad_rmesh_build_before_init, should_fail=.true.), &
         & new_unittest("bad_rmesh_seed_before_build", test_bad_rmesh_seed_before_build, should_fail=.true.), &
         & new_unittest("bad_rmesh_radii", test_bad_rmesh_radii, should_fail=.true.), &
         & new_unittest("bad_rmesh_reach", test_bad_rmesh_reach, should_fail=.true.), &
         & new_unittest("bad_rmesh_lattice_range", test_bad_rmesh_lattice_range, should_fail=.true.), &
         & new_unittest("bad_rmesh_sample_size_infinite", test_bad_rmesh_sample_size_infinite, should_fail=.true.), &
         & new_unittest("bad_rmesh_cone_scale", test_bad_rmesh_cone_scale, should_fail=.true.), &
         & new_unittest("bad_rmesh_cone_scale_infinite", test_bad_rmesh_cone_scale_infinite, should_fail=.true.), &
         & new_unittest("bad_rmesh_anchors_without_atoms", test_bad_rmesh_anchors_without_atoms, should_fail=.true.), &
         & new_unittest("rmesh_back_facing_needs_owner", test_rmesh_back_facing_needs_owner), &
         & new_unittest("bad_batch_target_size", test_bad_batch_target_size, should_fail=.true.), &
         & new_unittest("bad_newton_unreached_anchor", test_bad_newton_unreached_anchor, should_fail=.true.), &
         & new_unittest("worker_skipped_anchor_empty_outputs", test_worker_skipped_anchor) &
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
      real(wp) :: positions(3, 1), pos2(3, 2), anchor(3), center(3), direction(3)
      real(wp) :: atom_lo(3), atom_hi(3)
      !> Rigid translation testing cancellation in coordinate differences
      real(wp) :: shift(3)
      !> Deterministic directions and translations
      integer :: i, trial
      logical :: can_exclude

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
      if (allocated(error)) return

      ! Distant second atom encloses the anchor globally along x, while localized
      ! screening bounds still exclude the transverse patch near atom 1
      pos2 = reshape([0.0_wp, 0.0_wp, 0.0_wp, 40.0_wp, 0.0_wp, 0.0_wp], [3, 2])
      call new(mol, [1, 1], pos2)
      call lsf%update(mol, [1.0_wp, 1.0_wp])
      anchor = [2.0_wp, 0.0_wp, 0.0_wp]
      center = [0.0_wp, 1.0_wp, 0.0_wp]
      call lsf%stationary_region_bounds(anchor, center - 0.5_wp, center + 0.5_wp, &
                                        atom_lo, atom_hi, can_exclude)
      call check(error, can_exclude, message="Localized region bounds did not enable exclusion")
      if (allocated(error)) return
      call check(error, lsf%stationary_box_excluded(anchor, center, 0.1_wp, &
                                                    atom_lower=atom_lo, atom_upper=atom_hi), &
                 message="Localized region bounds did not exclude transverse patch")
      if (allocated(error)) return
      call check(error, lsf%stationary_box_excluded(anchor, center, 0.1_wp), &
                 message="Direct localized stationary_box_excluded did not exclude transverse patch")
      if (allocated(error)) return
      center = [1.0_wp, 0.0_wp, 0.0_wp]
      call check(error,.not. lsf%stationary_box_excluded(anchor, center, 0.1_wp, &
                                                         atom_lower=atom_lo, atom_upper=atom_hi), &
                 message="Localized bounds excluded analytic nearest point")
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

   !> Retain a distant strict minimum inside a wide branch window
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
                             error, dpdist_max=5.0_wp, require_pruning=.true.)
   end subroutine test_brute_force_broad_window

   !> Discovery follows the branch window: a wide one finds both minima, a narrow
   !> one only the near one
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

      call check_h2_branches(anchor, expected, error, dpdist_max=5.0_wp)
      if (allocated(error)) return
      call check_h2_branches(anchor, expected(:, 2:2), error, dpdist_max=0.01_wp)
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
   !> - the window radius behind it is zero, so the worker finalizes it without
   !>   a search
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
      !> Run context
      type(moist_context_type) :: ctx
      !> Molecular geometry
      type(structure_type) :: mol
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Atom at the origin and exact point on a unit sphere
      real(wp) :: positions(3, 1), anchor(3)
      !> Returned points
      real(wp), allocatable :: points(:, :)
      !> Returned point count
      integer :: n_points

      positions = 0.0_wp
      anchor = [1.0_wp, 0.0_wp, 0.0_wp]
      call new(mol, [1], positions)
      param%branch%dpdist_max = 0.0_wp
      param%proj_level = 4
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      call lsf%update(mol, [1.0_wp])
      call worker%init(param, lsf, positions, native_error)
      call new_context(ctx, verbosity=0, nthreads=1)
      if (.not. allocated(native_error)) call worker%project_anchor(ctx, anchor, 1, n_points, native_error, points=points)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, 1)
      if (allocated(error)) return
      call check(error, norm2(points(:, 1) - anchor), 0.0_wp, thr=1.0e-14_wp)
      if (allocated(error)) return
      call check(error, worker%n_seeded, 0, message="a zero window was searched")
   end subroutine test_zero_window

   !> A skipped anchor gives empty outputs from a fresh worker
   !>
   !> - a switching value below the weight cut skips the anchor, so the row
   !>   buffer of the fresh worker never grew
   !>
   !> @param[out] error Outputs missing or not empty
   subroutine test_worker_skipped_anchor(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Projection settings
      type(moist_cavity_drop_parameters_type) :: param
      !> One-sphere level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Fresh projection worker
      type(projection_worker_type) :: worker
      !> Run context
      type(moist_context_type) :: ctx
      !> Molecular geometry
      type(structure_type) :: mol
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Atom at the origin and anchor on its unit sphere
      real(wp) :: positions(3, 1), anchor(3)
      !> Returned rows
      real(wp), allocatable :: points(:, :), normals(:, :), lambda(:), pdist(:), phi(:), weights(:)
      !> Returned point count
      integer :: n_points

      positions = 0.0_wp
      anchor = [1.0_wp, 0.0_wp, 0.0_wp]
      call new(mol, [1], positions)
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      call lsf%update(mol, [1.0_wp])
      call worker%init(param, lsf, positions, native_error)
      call new_context(ctx, verbosity=0, nthreads=1)
      if (.not. allocated(native_error)) then
         call worker%project_anchor(ctx, anchor, 1, n_points, native_error, points=points, normals=normals, &
                                    lambda=lambda, pdist=pdist, phi=phi, weights=weights, f=0.0_wp)
      end if
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, 0, message="a skipped anchor kept rows")
      if (allocated(error)) return
      call check(error, allocated(points) .and. allocated(normals) .and. allocated(lambda) .and. allocated(pdist) &
                 .and. allocated(phi) .and. allocated(weights), message="an output of a skipped anchor is missing")
      if (allocated(error)) return
      call check(error, all(shape(points) == [3, 0]) .and. all(shape(normals) == [3, 0]), &
                 message="points or normals of a skipped anchor are not (3, 0)")
      if (allocated(error)) return
      call check(error, size(lambda) + size(pdist) + size(phi) + size(weights), 0, &
                 message="a row vector of a skipped anchor is not empty")
   end subroutine test_worker_skipped_anchor

   !> Compare the level-4 projection with independently enumerated minima
   !>
   !> Reference construction: solve d/dtheta ||y(theta) - anchor||^2 = 0 on
   !> the exact two-focus meridian, retaining positive tangent curvature
   !> No projected neighbour points enter the production solver
   !>
   !> - the worker runs with the weight cut disabled, so the raw branch set is
   !>   compared rather than the rows above the cut
   !>
   !> @param[in] anchor Initial point, in Bohr
   !> @param[in] expected Distinct strict minima, shape (3, n), in Bohr
   !> @param[out] error Missing minimum or solver failure
   !> @param[in] dpdist_max Optional branch window, in Bohr
   !> @param[in] transform Rotate and translate the geometry when true
   !> @param[in] require_pruning Require fewer seeds than the exhaustive reference
   subroutine check_h2_branches(anchor, expected, error, dpdist_max, transform, require_pruning)
      !> Reference point
      real(wp), intent(in) :: anchor(3)
      !> Independent minima
      real(wp), intent(in) :: expected(:, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Weight model scale
      real(wp), intent(in), optional :: dpdist_max
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
      !> Run context
      type(moist_context_type) :: ctx
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
      integer :: n_points, i, j
      !> Distance from one reference minimum to the nearest native branch
      real(wp) :: distance, pdist_start
      !> Returned and unpruned surface-leaf reference points
      real(wp), allocatable :: points(:, :), reference(:, :)
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
      param%branch%dpdist_max = 5.0_wp
      param%proj_level = 4
      ! The reach is no part of these searches: beyond every window they ask for
      param%branch%reach_on = 10.0_wp
      param%branch%reach_off = 12.0_wp
      if (present(dpdist_max)) param%branch%dpdist_max = dpdist_max
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
      worker%finalize%wleb_cut = -1.0_wp
      call new_context(ctx, verbosity=0, nthreads=1)
      call worker%project_anchor(ctx, point, 1, n_points, native_error, points=points)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, size(expected, 2), message="Unexpected or duplicate branches")
      if (allocated(error)) return
      do i = 1, size(expected, 2)
         distance = huge(1.0_wp)
         do j = 1, n_points
            distance = min(distance, norm2(points(:, j) - targets(:, i)))
         end do
         call check(error, distance < 1.0e-8_wp, message="Search missed a strict minimum")
         if (allocated(error)) return
      end do

      ! The brute force covers the largest ball the search asked for
      pdist_start = worker%pdist_search_max
      call brute_force_branches(lsf, point, param, pdist_start, reference, native_error, reference_seeds)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, size(reference, 2), size(expected, 2), message="Brute-force H2 count differs")
      if (allocated(error)) return
      call check_branch_sets(reference, points(:, 1:n_points), error)
      if (allocated(error)) return
      if (present(require_pruning)) then
         if (require_pruning) then
            call check(error, 2*worker%n_seeds < reference_seeds, &
                       message="Broad H2 search still solves most surface leaves")
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
   !> Every surface leaf of the ball ([[surface_leaves]]) seeds Newton; minima
   !> inside the ball only
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
      !> Geometry evaluator, separate from the finder under test
      type(moist_cavity_drop_lsf_svdw_type) :: probe
      !> Centers of the surface leaves
      real(wp), allocatable :: leaves(:, :)
      !> Local refinement of every leaf
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      !> Run context
      type(moist_context_type) :: ctx
      !> Candidate points and multipliers
      real(wp), allocatable :: candidates(:, :), multipliers(:), gradients(:, :), hessians(:, :, :)
      !> Distance-window and deduplication mask
      logical, allocatable :: keep(:)
      !> Surface data and tangent Hessian
      real(wp) :: t1(3), t2(3), b11, b12, b22, beta, p(3, 1), v(1), gradient(3)
      !> Seed and accepted point counts
      integer :: i, n, n_strict, icol
      !> Branch window of the settings
      type(branch_weight_type) :: bw

      if (present(n_seeds)) n_seeds = 0
      probe = lsf
      p(:, 1) = anchor
      call probe%evaluate_batch(p, 0, v, error=error)
      if (allocated(error)) return
      call bw%init(param%branch%dpdist_max, param%phi_alpha)
      call surface_leaves(probe, anchor, radius, param%seed%seed_size, leaves, error)
      if (allocated(error)) return
      if (present(n_seeds)) n_seeds = size(leaves, 2)

      ncfg%tol = param%solve%tol
      ncfg%max_iter = param%solve%max_iter
      call batch%init(size(leaves, 2) + 1, error)
      if (allocated(error)) return
      call batch%begin_run(1, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      do i = 1, size(leaves, 2)
         call batch%append(leaves(:, i), 0.0_wp, icol)
      end do
      call newton_solve_batch(batch, probe, ncfg, pcfg, param%phi_alpha, scratch, error)
      if (allocated(error)) return
      call new_context(ctx, verbosity=0, nthreads=1)
      call polish_batch(ctx, batch, probe, pcfg, ncfg, param%phi_alpha, scratch, error)
      if (allocated(error)) return

      allocate (candidates(3, batch%n), multipliers(batch%n), gradients(3, batch%n), hessians(3, 3, batch%n))
      n = 0
      do i = 1, batch%n
         if (batch%status(i) /= point_accepted) cycle
         if (abs(batch%S(i)) > 1.0e-9_wp) cycle
         ! Only the searched ball: a leaf touching it may converge beyond
         if (norm2(batch%x(:, i) - anchor) > radius) cycle
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
      call dedupe_candidates(candidates(:, 1:n), anchor, param%branch%sep_cut, bw%x_max**2, keep, &
                             pdist_slack=bw%x_max, pdist_soft=bw%soft)
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

   !> Centers of every leaf of edge at most `h` that may hold surface inside a ball
   !>
   !> - exhaustive box descent from the cube around the ball, independent of the
   !>   seeders: a box is dropped only when it misses the ball or the exclusion
   !>   radius at its center covers it
   !>
   !> @param[inout] lsf    Level set, updated for the geometry
   !> @param[in]    anchor Ball center
   !> @param[in]    radius Ball radius (Bohr)
   !> @param[in]    h      Largest leaf edge (Bohr)
   !> @param[out]   leaves Leaf centers (3, n)
   !> @param[out]   error  Level set failure
   subroutine surface_leaves(lsf, anchor, radius, h, leaves, error)
      type(moist_cavity_drop_lsf_svdw_type), intent(inout) :: lsf
      real(wp), intent(in) :: anchor(3)
      real(wp), intent(in) :: radius
      real(wp), intent(in) :: h
      real(wp), allocatable, intent(out) :: leaves(:, :)
      type(moist_error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: centers(:, :), next(:, :), values(:)
      real(wp) :: edge, half_diag, gap(3)
      integer :: i, n, nnext, ix, iy, iz
      logical :: keep

      edge = 2.0_wp*radius
      centers = reshape(anchor, [3, 1])
      allocate (leaves(3, 0))
      do while (size(centers, 2) > 0)
         allocate (values(size(centers, 2)))
         call lsf%evaluate_batch(centers, 0, values, error=error)
         if (allocated(error)) return
         half_diag = 0.5_wp*sqrt(3.0_wp)*edge
         allocate (next(3, 8*size(centers, 2)))
         nnext = 0
         n = 0
         do i = 1, size(centers, 2)
            gap = max(abs(centers(:, i) - anchor) - 0.5_wp*edge, 0.0_wp)
            keep = norm2(gap) <= radius .and. lsf%exclusion_radius(values(i)) < half_diag
            if (.not. keep) cycle
            if (edge <= h) then
               n = n + 1
               centers(:, n) = centers(:, i)
               cycle
            end if
            do iz = -1, 1, 2
               do iy = -1, 1, 2
                  do ix = -1, 1, 2
                     nnext = nnext + 1
                     next(:, nnext) = centers(:, i) + 0.25_wp*edge*real([ix, iy, iz], wp)
                  end do
               end do
            end do
         end do
         leaves = reshape([leaves, centers(:, 1:n)], [3, size(leaves, 2) + n])
         centers = next(:, 1:nnext)
         deallocate (next, values)
         edge = 0.5_wp*edge
      end do
   end subroutine surface_leaves

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
                                  mol, 50, 1.0_wp)
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

   !> Require level-4 recovery of the admissible minima level 3 finds
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
      !> Level-3 and level-4 projections
      type(projection_rows_type) :: comparison, fine
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
      call project_cluster(cavity, 4, fine, native_error)
      if (allocated(native_error)) then
         call test_failed(error, title//": level 4: "//native_error%message)
         return
      end if
      call check_projection_geometry(cavity, lsf, comparison, error)
      if (allocated(error)) return
      call check_projection_geometry(cavity, lsf, fine, error)
      if (allocated(error)) return

      call row_bounds(fine, cavity%anchor%ngrid, first, last)
      nrelevant = 0
      do i = 1, comparison%n
         if (.not. relevant_row(cavity, comparison, i)) cycle
         nrelevant = nrelevant + 1
         aid = comparison%anchor(i)
         matches = 0
         if (first(aid) > 0) then
            do j = first(aid), last(aid)
               if (.not. relevant_row(cavity, fine, j)) cycle
               if (norm2(comparison%xyz(:, i) - fine%xyz(:, j)) <= 1.0e-6_wp) matches = matches + 1
            end do
         end if
         write (message, "(a,a,i0,a,i0)") title, ": anchor ", cavity%anchor%id(aid), &
            ", matching level-4 branches ", matches
         call check(error, matches, 1, message=trim(message))
         if (allocated(error)) return
      end do
      call check(error, nrelevant > 0, message="Comparison set contains no admissible minimum")
      if (allocated(error)) return
      call check_cluster_brute_force(cavity, lsf, fine, error)
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
      !> Level-4 settings of the compared rows
      type(moist_cavity_drop_parameters_type) :: param
      !> Branch window of the settings
      type(branch_weight_type) :: bw

      param = cavity%param
      param%proj_level = 4
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call bw%init(param%branch%dpdist_max, param%phi_alpha)
      call row_bounds(rows, cavity%anchor%ngrid, first, last)
      total = 0
      do aid = 1, cavity%anchor%ngrid
         if (first(aid) == 0) cycle
         anchor = cavity%anchor%xyz(:, aid)
         radius = huge(1.0_wp)
         do j = first(aid), last(aid)
            radius = min(radius, sum((rows%xyz(:, j) - anchor)**2))
         end do
         ! No search goes beyond the reach
         radius = min(bw%pdist_max(sqrt(radius)), param%branch%reach_off)
         call brute_force_branches(lsf, anchor, param, radius, reference, native_error)
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
      call cavity%anchor%filter(cavity%param%finalize%wleb_cut, nremoved)
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
      call param%compute_derived(error)
      if (allocated(error)) return
      call worker%init(param, cavity%lsf_model, cavity%mol%xyz, error)
      if (allocated(error)) return
      worker%finalize%wleb_cut = -1.0_wp
      allocate (nrows(cavity%anchor%ngrid), lost(cavity%anchor%ngrid))
      slots = [(i, i=1, cavity%anchor%ngrid)]
      call worker%process_chunk(cavity%ctx, slots, cavity%anchor%xyz(:, 1:cavity%anchor%ngrid), &
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
                     rows%cpjac(i)*rows%w_f(i)) > cavity%param%finalize%wleb_cut
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

   !> Newton from the anchor lands on the analytic point of a single-atom sphere
   !>
   !> @param[out] error Wrong point or multiplier
   subroutine test_newton_sphere_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchors, points
      real(wp) :: rs, anchors(3, 4), point(3), u(3), residual(3)
      !> Steps per column
      integer :: nsteps(4)
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
      call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed, steps=nsteps)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if

      do k = 1, 4
         call check(error, batch%status(k), point_converged, "anchor converges")
         if (allocated(error)) return
         u = unit_vector(anchors(:, k))
         point = rs*u
         call check(error, norm2(batch%x(:, k) - point) < 1.0e-9_wp, "point is the radial projection")
         if (allocated(error)) return
         call check(error, abs(batch%S(k)) < 1.0e-10_wp, "final jet is on the surface")
         if (allocated(error)) return
         residual = default_phi_alpha()*(batch%x(:, k) - anchors(:, k)) - batch%lambda(k)*batch%gradS(:, k)
         call check(error, norm2(residual) < 1.0e-9_wp, "final multiplier is stationary")
         if (allocated(error)) return
         call check(error, abs(batch%pdist(k) - abs(norm2(anchors(:, k)) - rs)) < 1.0e-9_wp, "pdist")
         if (allocated(error)) return
         call check(error, abs(abs(dot_product(batch%normal(:, k), u)) - 1.0_wp) < 1.0e-10_wp, &
                    "normal is radial")
         if (allocated(error)) return
         call check(error, .not. any(ieee_is_nan(batch%hessS(:, k))), "final jet carries the Hessian")
         if (allocated(error)) return
         call check(error, batch%run_jet(k), "anchor jet recorded without an extra evaluation")
         if (allocated(error)) return
         call check(error, nsteps(k) <= 20, "converges in a handful of steps")
         if (allocated(error)) return
      end do
   end subroutine test_newton_sphere_point

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
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
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
      call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, batch%status(icol_far), point_failed, "far seed is unreached")
      if (allocated(error)) return
      call check(error, batch%status(icol_anchor), point_converged, "anchor seed converges")
      if (allocated(error)) return
      call check(error, ieee_is_nan(batch%S(icol_far)), "unreached jet is never stored")
   end subroutine test_newton_far_seed_unreached

   !> Newton from an anchor far outside the sphere walks the whole distance to its point
   !>
   !> - the point lies more than 15 Bohr from the anchor along x, many `max_step`
   !>   steps away
   !>
   !> @param[out] error Far point not reached
   subroutine test_newton_far_point(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, analytic point and anchor value
      real(wp) :: rs, anchor(3, 1), point(3), v(1)
      !> Column
      integer :: icol
      !> Absolute point tolerance (Bohr)
      real(wp), parameter :: tol = 1.0e-9_wp

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return

      anchor(:, 1) = [rs + 18.0_wp, 0.0_wp, 0.0_wp]
      point = [rs, 0.0_wp, 0.0_wp]
      call check(error, abs(anchor(1, 1) - point(1)) > 15.0_wp, "fixture point lies beyond 15 Bohr")
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
      call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, batch%status(icol), point_converged, "far anchor converges")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, icol) - point) <= tol .or. &
                 norm2(batch%x(:, icol) - point) <= 10.0_wp*tol*rs, "point is the radial projection")
      if (allocated(error)) return
      call check(error, abs(batch%pdist(icol) - 18.0_wp) <= tol .or. &
                 abs(batch%pdist(icol) - 18.0_wp) <= 10.0_wp*tol*18.0_wp, "pdist is the full distance")
   end subroutine test_newton_far_point

   !> The two starts of a column: at `lambda = 0` the first step needs no
   !> Hessian, at `lambda /= 0` the start jet is refreshed to order 2
   !>
   !> - one column per solve, on a counting level set that counts the order-2
   !>   evaluations at the seed: none at `lambda = 0`, the one refresh at
   !>   `lambda /= 0`; every later iterate lies off the seed
   !>
   !> @param[out] error Wrong point, step count or Hessian evaluations
   subroutine test_newton_start_multiplier(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Counting single-atom level set
      type(counting_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, seed, analytic point
      real(wp) :: rs, anchor(3), seed(3), point(3)
      !> Start multipliers
      real(wp), parameter :: lambdas(2) = [0.0_wp, -0.5_wp]
      !> Steps of the column
      integer :: nsteps(1)
      !> Column and start
      integer :: icol, k

      call sphere_lsf(lsf%moist_cavity_drop_lsf_svdw_type, rs, error)
      if (allocated(error)) return
      call batch%init(4, failed)

      ! Off-axis seed on the surface, once without and once with a multiplier
      anchor = [rs + 1.0_wp, 0.0_wp, 0.0_wp]
      seed = rs*[cos(0.3_wp), sin(0.3_wp), 0.0_wp]
      point = rs*unit_vector(anchor)
      lsf%probe = seed
      do k = 1, 2
         call batch%clear()
         call batch%begin_run(1, anchor, 1)
         call batch%append(seed, lambdas(k), icol)
         lsf%n_probe_hess = 0
         call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed, steps=nsteps)
         if (allocated(failed)) then
            call test_failed(error, "batched Newton: "//failed%message)
            return
         end if
         call check(error, batch%status(1), point_converged, "column converges")
         if (allocated(error)) return
         call check(error, norm2(batch%x(:, 1) - point) < 1.0e-9_wp, "analytic point")
         if (allocated(error)) return
         call check(error, nsteps(1) >= 1 .and. nsteps(1) <= 20, "steps counted")
         if (allocated(error)) return
         call check(error, .not. any(ieee_is_nan(batch%hessS(:, 1))), "final jet carries the Hessian")
         if (allocated(error)) return
         call check(error, lsf%n_probe_hess, k - 1, "order-2 evaluations at the seed: none at lambda = 0, one else")
         if (allocated(error)) return
      end do
   end subroutine test_newton_start_multiplier

   !> A singular KKT system: a singular tangent block takes the flat step, a
   !> regular one is a vanishing gradient
   !>
   !> - one sphere, the column on its surface on the anchor's axis, where both
   !>   tangent curvatures of `S` are `kappa`
   !> - start multiplier `w / kappa`: the tangent block `w - lambda kappa`
   !>   vanishes; the flat step resets the multiplier and Newton converges on the point
   !> - start multiplier `1e17`: the tangent block is regular, the gradient
   !>   (about one) is below the round-off of `lambda Hess S`, the system singular:
   !>   `point_failed` before any step
   !>
   !> @param[out] error Wrong status, point or step count
   subroutine test_newton_singular_system(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, start, its jet, tangent curvature, start multipliers
      real(wp) :: rs, anchor(3), start(3, 1), v(1), g(3, 1), h(6, 1), kappa, lambdas(2)
      !> Status each start must end with, steps of the column
      integer :: expected(2), nsteps(1)
      !> Column and case
      integer :: icol, k

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      anchor = [rs + 1.0_wp, 0.0_wp, 0.0_wp]
      start(:, 1) = [rs, 0.0_wp, 0.0_wp]
      call lsf%evaluate_batch(start, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "start jet: "//failed%message)
         return
      end if
      ! Packed yy: a tangent direction of the start
      kappa = h(3, 1)
      call check(error, kappa > 0.0_wp, "the level set curves along the tangents")
      if (allocated(error)) return

      lambdas = [default_phi_alpha()/kappa, 1.0e17_wp]
      expected = [point_converged, point_failed]
      call batch%init(4, failed)
      do k = 1, 2
         call batch%clear()
         call batch%begin_run(1, anchor, 1)
         call batch%append(start(:, 1), lambdas(k), icol)
         call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed, steps=nsteps)
         if (allocated(failed)) then
            call test_failed(error, "batched Newton: "//failed%message)
            return
         end if
         call check(error, batch%status(1), expected(k), "singular tangent block steps, regular one fails")
         if (allocated(error)) return
      end do
      call check(error, nsteps(1), 0, "the vanishing gradient fails before any step")
      if (allocated(error)) return
   end subroutine test_newton_singular_system

   !> A column whose level set value is not finite fails alone
   !>
   !> - no bundled level set hands out such a value inside its reach; a host
   !>   density callback can
   !> - the anchor column of the same run takes the steps it takes on its own:
   !>   the line search must not try it in place of the other column
   !>
   !> @param[out] error Non-finite column not failed, or the anchor column disturbed
   subroutine test_newton_non_finite_value(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set, its value NaN at the seed
      type(nan_value_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, seed
      real(wp) :: rs, anchor(3), seed(3)
      !> Steps per column, steps of the anchor column on its own
      integer :: nsteps(2), nsteps_alone
      !> Column and solve
      integer :: icol, k

      call sphere_lsf(lsf%moist_cavity_drop_lsf_svdw_type, rs, error)
      if (allocated(error)) return
      anchor = [rs + 2.0_wp, 0.0_wp, 0.0_wp]
      seed = (rs + 0.5_wp)*[cos(0.3_wp), sin(0.3_wp), 0.0_wp]
      call batch%init(4, failed)

      ! The anchor column alone, then followed by the seed with the NaN value
      nsteps_alone = -1
      do k = 1, 2
         call batch%clear()
         call batch%begin_run(1, anchor, 1)
         call batch%append(anchor, 0.0_wp, icol)
         if (k == 2) then
            call batch%append(seed, 0.0_wp, icol)
            lsf%poisoned = seed
         end if
         call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed, steps=nsteps)
         if (allocated(failed)) then
            call test_failed(error, "batched Newton: "//failed%message)
            return
         end if
         call check(error, batch%status(1), point_converged, "anchor column converges")
         if (allocated(error)) return
         if (k == 1) nsteps_alone = nsteps(1)
      end do
      call check(error, nsteps(1), nsteps_alone, "anchor column takes the steps it takes alone")
      if (allocated(error)) return
      call check(error, batch%status(2), point_failed, "non-finite column fails")
   end subroutine test_newton_non_finite_value

   !> A Newton point stuck on the far pole (a constrained maximum) is escaped
   !> and re-polished onto the near pole; the anchor's own point is accepted as is
   !>
   !> @param[out] error Wrong classification or escape
   subroutine test_polish_far_pole(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
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
      call new_context(ctx, verbosity=0, nthreads=1)

      call batch%init(4, failed)
      call batch%begin_run(1, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      call batch%append([-rs, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:2) == point_converged), "both starts converge")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 2) + pole) < 1.0e-9_wp, "second start sits on the far pole")
      if (allocated(error)) return

      call polish_batch(ctx, batch, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "polish: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:2) == point_accepted), "both end as verified minima")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - pole) < 1.0e-9_wp, "anchor point untouched")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 2) - pole) < 1.0e-8_wp, "far pole escaped onto the near pole")
   end subroutine test_polish_far_pole

   !> A converged column whose jet is not finite is never accepted
   !>
   !> - one sphere, the column on the anchor's own point, a minimum
   !> - its Hessian set to NaN, as a host density callback can hand one back
   !>   without an error
   !>
   !> @param[out] error The column was kept
   subroutine test_polish_non_finite_jet(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor
      real(wp) :: rs, anchor(3)
      !> Column
      integer :: icol

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      anchor = [rs + 2.0_wp, 0.0_wp, 0.0_wp]
      call new_context(ctx, verbosity=0, nthreads=1)

      call batch%init(4, failed)
      call batch%begin_run(1, anchor, 1)
      call batch%append([rs, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, batch%status(1), point_converged, "the start converges")
      if (allocated(error)) return

      batch%hessS(:, 1) = ieee_value(0.0_wp, ieee_quiet_nan)
      call polish_batch(ctx, batch, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "polish: "//failed%message)
         return
      end if
      call check(error, batch%status(1), point_failed, "a NaN Hessian fails the column")
   end subroutine test_polish_non_finite_jet

   !> Repeat visits of a saddle with one descending sign are duplicates, not failures
   !>
   !> - carbon pair, a neck point just off the mirror plane, the anchor a little
   !>   beyond the focal point of the meridian curvature: a weak saddle its third
   !>   order dominates, so one meridian sign descends and the other climbs, while
   !>   the two objectives agree within the symmetry tolerance of the saddle cache
   !> - four columns on the saddle with `capture_duplicates`: the first escapes to
   !>   a minimum, each repeat is its duplicate and stays where it is
   !>
   !> @param[out] error A repeat that failed, or no escape of the first visit
   subroutine test_polish_one_sided_saddle(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Two-atom level set and its geometry
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(structure_type) :: mol
      !> Points and solver state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Objective weight, saddle with its jet, normal, meridian tangent and curvature
      real(wp) :: w, saddle(3, 1), v(1), g(3, 1), h(6, 1), hess(3, 3), n(3), tm(3), km
      !> Multiplier and anchor that make the point a saddle
      real(wp) :: lambda, anchor(3)
      !> Column and visit
      integer :: icol, j
      !> Relative distance of the anchor beyond the focal point
      real(wp), parameter :: eps = 1.0e-5_wp

      call new(mol, [6, 6], reshape([-1.6_wp, 0.0_wp, 0.0_wp, 1.6_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, blend_3b=1.0_wp))
      lsf%screening_threshold = 1.0e-11_wp
      call lsf%update(mol, [3.0_wp, 3.0_wp])
      call new_context(ctx, verbosity=0, nthreads=1)
      w = default_phi_alpha()

      ! Neck point in the xy plane: meridian tangent is the normal turned by 90 degrees
      call retract_to_surface(pcfg, w, lsf, [0.02_wp, 2.5_wp, 0.0_wp], [0.0_wp, 0.0_wp, 0.0_wp], &
                              saddle(:, 1), lambda, failed)
      if (allocated(failed)) then
         call test_failed(error, "neck point: "//failed%message)
         return
      end if
      call lsf%evaluate_batch(saddle, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "neck jet: "//failed%message)
         return
      end if
      n = unit_vector(g(:, 1))
      tm = [n(2), -n(1), 0.0_wp]
      call unpack_sym_matrix(h(:, 1), hess, 3)
      km = dot_product(tm, matmul(hess, tm))
      call check(error, km < 0.0_wp, "the neck curves away from the anchor along the meridian")
      if (allocated(error)) return
      lambda = w*(1.0_wp + eps)/km
      anchor = saddle(:, 1) - (lambda/w)*g(:, 1)

      call batch%init(4, failed)
      call batch%begin_run(1, anchor, 1)
      do j = 1, 4
         call batch%append(saddle(:, 1), lambda, icol)
      end do
      call newton_solve_batch(batch, lsf, ncfg, pcfg, w, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "batched Newton: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:4) == point_converged), "every visit converges on the saddle")
      if (allocated(error)) return

      ncfg%capture_duplicates = .true.
      call polish_batch(ctx, batch, lsf, pcfg, ncfg, w, scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "polish: "//failed%message)
         return
      end if
      call check(error, batch%status(1), point_accepted, "the first visit escapes to a minimum")
      if (allocated(error)) return
      call check(error, norm2(batch%x(:, 1) - saddle(:, 1)) > 1.0e-3_wp, "away from the saddle")
      if (allocated(error)) return
      do j = 2, 4
         call check(error, batch%status(j), point_captured, "a repeat visit is a duplicate, not a failed saddle")
         if (allocated(error)) return
         call check(error, all(batch%x(:, j) == saddle(:, 1)), "and stays on the saddle")
         if (allocated(error)) return
      end do
   end subroutine test_polish_one_sided_saddle

   !> The seeded pass stops starts that head into the accepted local point
   !>
   !> - one sphere at level 3: every seed of the anchor ends on the
   !>   local point, so all of them are duplicates
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
      !> Run context
      type(moist_context_type) :: ctx
      !> One-atom geometry
      type(structure_type) :: mol
      type(moist_error_type), allocatable :: failed
      !> Atom position, anchor
      real(wp) :: positions(3, 1), anchor(3)
      !> Rows with capture on (1) and off (2)
      real(wp), allocatable :: points_on(:, :), points_off(:, :), w_on(:), w_off(:)
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
      call new_context(ctx, verbosity=0, nthreads=1)

      call project(.true., 1, n_on, points_on, w_on)
      if (allocated(error)) return
      call project(.false., 2, n_off, points_off, w_off)
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
         call check(error, same_value(points_on(:, k), points_off(:, k), 1.0e-9_wp), "the same points")
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
      !> @param[out] points    Points (3, n)
      !> @param[out] w       Branch weights
      subroutine project(capture, run, n, points, w)
         logical, intent(in) :: capture
         integer, intent(in) :: run
         integer, intent(out) :: n
         real(wp), allocatable, intent(out) :: points(:, :)
         real(wp), allocatable, intent(out) :: w(:)
         integer :: n0

         worker%solve_seeded%capture_duplicates = capture
         n0 = worker%n_captured
         call reset_count(worker%lsf)
         call worker%project_anchor(ctx, anchor, 1, n, failed, points=points, weights=w)
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

   !> Close twin points both survive capture
   !>
   !> - two overlapping carbons; the anchor on their mirror plane lies just
   !>   beyond the focal point of the neck, so its two mirror points are close
   !> - a run as the seeded pass builds it: the accepted point on one side
   !>   first, then starts across the whole neck, the middle one on the saddle
   !> - with capture on, the run keeps both points, a start of each side is
   !>   captured, and the points are those of the run with capture off
   !>
   !> @param[out] error A twin lost or a point moved
   subroutine test_capture_keeps_twin_points(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Two-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Run of the accepted point, run under test
      type(projection_batch_type) :: local, batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      type(structure_type) :: mol
      !> Anchor, distinct points with capture on (1) and off (2)
      real(wp) :: anchor(3), points(3, 32, 2), mirror(3)
      !> Distinct points and captured starts per run, column, start, match
      integer :: npoint(2), captured(2), icol, k, run, i, m

      call new(mol, [6, 6], reshape([-2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, blend_3b=1.0_wp))
      lsf%screening_threshold = 1.0e-11_wp
      call lsf%update(mol, [3.0_wp, 3.0_wp])
      anchor = [0.0_wp, 3.67_wp, 0.0_wp]
      call new_context(ctx, verbosity=0, nthreads=1)

      ! Accepted point on the +x side
      call local%init(4, failed)
      call local%begin_run(1, anchor, 1)
      call local%append([0.6_wp, 2.6_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(local, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
      if (.not. allocated(failed)) call polish_batch(ctx, local, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "local point: "//failed%message)
         return
      end if
      call check(error, local%status(1), point_accepted, "local point verified")
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
         call newton_solve_batch(batch, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
         if (.not. allocated(failed)) call polish_batch(ctx, batch, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
         if (allocated(failed)) then
            call test_failed(error, "seeded run: "//failed%message)
            return
         end if
         captured(run) = count(batch%status(1:batch%n) == point_captured)
         call distinct_points(batch, points(:, :, run), npoint(run))
      end do

      call check(error, captured(1) > 0, "the run captures duplicate starts")
      if (allocated(error)) return
      call check(error, npoint(1), 2, "both twin points kept")
      if (allocated(error)) return
      mirror = [-points(1, 2, 1), points(2, 2, 1), points(3, 2, 1)]
      call check(error, same_value(points(:, 1, 1), mirror, 1.0e-7_wp), "the two points are mirror images")
      if (allocated(error)) return
      call check(error, norm2(points(:, 1, 1) - points(:, 2, 1)) < 0.1_wp, "the fixture's twins are close")
      if (allocated(error)) return
      do i = 1, npoint(1)
         m = 0
         do k = 1, npoint(2)
            if (same_value(points(:, i, 1), points(:, k, 2), 1.0e-7_wp)) m = k
         end do
         call check(error, m > 0, "every point is a point of the run with capture off")
         if (allocated(error)) return
      end do
      do k = 1, npoint(2)
         m = 0
         do i = 1, npoint(1)
            if (same_value(points(:, i, 1), points(:, k, 2), 1.0e-7_wp)) m = i
         end do
         call check(error, m > 0, "no point of the run with capture off is lost")
         if (allocated(error)) return
      end do

   contains

      !> Accepted points of the run, deduplicated as the finalize does
      !>
      !> @param[in]  b     Solved and polished run
      !> @param[out] pts   Distinct points (3, :)
      !> @param[out] n     Distinct points found
      subroutine distinct_points(b, pts, n)
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
      end subroutine distinct_points

   end subroutine test_capture_keeps_twin_points

   !> Capture stops a column only on a Newton prediction
   !>
   !> - one sphere, anchor `0.05` Bohr above it; the run leads with the accepted
   !>   point `p`, the second column starts on the surface `1e-3` rad off it, where
   !>   the anchor's projection onto the tangent plane, `pdist 1e-3` from `p`, is
   !>   within `capture_distance` of it, while the start itself is not
   !> - `lambda = 0`: the flat step is the full system, a prediction; captured
   !> - multiplier of `p`: the full step is the prediction; captured
   !> - a large positive multiplier: the tangent block `w I - lambda Hess S` is
   !>   negative, the full step no merit descent; the flat step taken instead lands
   !>   on `p` as well, but it is a fallback direction: with a one-step budget the
   !>   column runs out of iterations instead
   !>
   !> @param[out] error A fallback step captured, or a prediction did not
   subroutine test_capture_only_predictions(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Single-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Run of the accepted point, run under test
      type(projection_batch_type) :: local, batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor, start, its jet, tangential curvature, start multipliers
      real(wp) :: rs, anchor(3), start(3, 1), v(1), g(3, 1), h(6, 1), kappa, lambdas(3)
      !> Status each start must end with
      integer :: expected(3)
      !> Columns and case
      integer :: icol, k
      !> Angle of the start off the axis (rad)
      real(wp), parameter :: theta = 1.0e-3_wp

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      anchor = [rs + 0.05_wp, 0.0_wp, 0.0_wp]
      call new_context(ctx, verbosity=0, nthreads=1)

      call local%init(4, failed)
      call local%begin_run(1, anchor, 1)
      call local%append(anchor, 0.0_wp, icol)
      call newton_solve_batch(local, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
      if (.not. allocated(failed)) call polish_batch(ctx, local, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "accepted point: "//failed%message)
         return
      end if
      call check(error, local%status(1), point_accepted, "accepted point verified")
      if (allocated(error)) return

      ! Off the axis in the xy-plane: z is a tangent direction of the start
      start(:, 1) = rs*[cos(theta), sin(theta), 0.0_wp]
      call check(error, norm2(start(:, 1) - local%x(:, 1)) > 1.0e-3_wp, "the start itself is no duplicate")
      if (allocated(error)) return
      call lsf%evaluate_batch(start, 2, v, lsf1_r=g, lsf2_rr_packed=h, error=failed)
      if (allocated(failed)) then
         call test_failed(error, "start jet: "//failed%message)
         return
      end if
      kappa = h(6, 1)
      call check(error, kappa > 0.0_wp, "the level set curves along the tangent z")
      if (allocated(error)) return

      lambdas = [0.0_wp, local%lambda(1), 4.0_wp*default_phi_alpha()/kappa]
      expected = [point_captured, point_captured, point_failed]
      ncfg%capture_duplicates = .true.
      ncfg%max_iter = 1
      call batch%init(4, failed)
      do k = 1, 3
         call batch%clear()
         call batch%begin_run(1, anchor, 1)
         call batch%append_copy(local, 1, icol)
         call batch%append(start(:, 1), lambdas(k), icol)
         call newton_solve_batch(batch, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
         if (allocated(failed)) then
            call test_failed(error, "batched Newton: "//failed%message)
            return
         end if
         call check(error, batch%status(1), point_accepted, "the accepted point stays accepted")
         if (allocated(error)) return
         call check(error, batch%status(2), expected(k), "capture only on a Newton prediction")
         if (allocated(error)) return
      end do
   end subroutine test_capture_only_predictions

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

      integer :: k

      self%n_points = self%n_points + size(points, 2)
      if (order >= 2) then
         do k = 1, size(points, 2)
            if (all(points(:, k) == self%probe)) self%n_probe_hess = self%n_probe_hess + 1
         end do
      end if
      call self%moist_cavity_drop_lsf_svdw_type%evaluate_batch(points, order, lsf0, lsf1_r, &
                                                               lsf2_rr_packed, lsf3_rrr, lsf4_rrrr, &
                                                               candidate_indices, error)
   end subroutine counting_evaluate_batch

   !> Evaluate as the SvdW level set, then replace the value at `poisoned` by NaN
   !>
   !> - a non-finite point is answered as `poisoned`: a NaN value and no error
   !>
   !> @param[inout] self              Level set with a NaN value at one point
   !> @param[in]    points            Evaluation points (3, n)
   !> @param[in]    order             Highest spatial derivative to write
   !> @param[out]   lsf0              Values (n)
   !> @param[out]   lsf1_r            Gradients (3, n)
   !> @param[out]   lsf2_rr_packed    Packed Hessians (6, n)
   !> @param[out]   lsf3_rrr          Packed third derivatives (10, n)
   !> @param[out]   lsf4_rrrr         Packed fourth derivatives (15, n)
   !> @param[in]    candidate_indices Shared atom list
   !> @param[out]   error             Evaluation failure
   subroutine nan_value_evaluate_batch(self, points, order, lsf0, lsf1_r, lsf2_rr_packed, &
                                       lsf3_rrr, lsf4_rrrr, candidate_indices, error)
      class(nan_value_svdw_type), intent(inout) :: self
      real(wp), intent(in), contiguous :: points(:, :)
      integer, intent(in) :: order
      real(wp), intent(out), contiguous :: lsf0(:)
      real(wp), intent(out), optional, contiguous :: lsf1_r(:, :)
      real(wp), intent(out), optional, contiguous :: lsf2_rr_packed(:, :)
      real(wp), intent(out), optional, contiguous :: lsf3_rrr(:, :)
      real(wp), intent(out), optional, contiguous :: lsf4_rrrr(:, :)
      integer, intent(in), optional :: candidate_indices(:)
      type(moist_error_type), allocatable, intent(out) :: error

      !> Points with every non-finite one moved to `poisoned`
      real(wp) :: finite_points(3, size(points, 2))
      integer :: k

      finite_points = points
      do k = 1, size(points, 2)
         if (.not. all(abs(points(:, k)) <= huge(1.0_wp))) finite_points(:, k) = self%poisoned
      end do
      call self%moist_cavity_drop_lsf_svdw_type%evaluate_batch(finite_points, order, lsf0, lsf1_r, &
                                                               lsf2_rr_packed, lsf3_rrr, lsf4_rrrr, &
                                                               candidate_indices, error)
      if (allocated(error)) return
      do k = 1, size(points, 2)
         if (all(finite_points(:, k) == self%poisoned)) lsf0(k) = ieee_value(0.0_wp, ieee_quiet_nan)
      end do
   end subroutine nan_value_evaluate_batch

   !> The canonical branch order is lexicographic with a tie tolerance
   !>
   !> @param[out] error Wrong order
   subroutine test_finalize_canonical_order(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Points and order
      real(wp) :: pts(3, 3)
      integer :: order(3)

      ! x ties within 1e-8 fall through to y
      pts(:, 1) = [1.0_wp, 2.0_wp, 0.0_wp]
      pts(:, 2) = [1.0_wp + 1.0e-10_wp, -1.0_wp, 0.0_wp]
      pts(:, 3) = [-3.0_wp, 5.0_wp, 0.0_wp]
      call canonical_branch_order(pts, order)
      call check(error, all(order == [3, 2, 1]), "lexicographic order with tie tolerance")
   end subroutine test_finalize_canonical_order

   !> The canonical order of a chain of near-ties does not depend on the input order
   !>
   !> - `a`, `b`, `c` step by 0.6e-8 in x, below the tie tolerance 1e-8, while
   !>   `a` and `c` are 1.2e-8 apart: a pairwise tie test orders `c < b < a < c`,
   !>   a cycle, and its sort result depends on the input order
   !> - the chain is one tie cluster, so y decides inside it; `d` lies far off in x
   !>   and sorts by x, whatever its y
   !> - every permutation of the input gives the same sorted points: c, b, a, d
   !>
   !> @param[out] error Order depends on the input order
   subroutine test_finalize_canonical_order_chain(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Tie tolerance of the canonical order (Bohr)
      real(wp), parameter :: tie = 1.0e-8_wp
      !> Fixture: a, b, c chained in x with y falling, d far off in x
      real(wp), parameter :: base(3, 4) = reshape([ &
         1.0_wp, 2.0_wp, 0.0_wp, &
         1.0_wp + 0.6e-8_wp, 1.0_wp, 0.0_wp, &
         1.0_wp + 1.2e-8_wp, 0.0_wp, 0.0_wp, &
         2.0_wp, -5.0_wp, 0.0_wp], [3, 4])
      !> Expected sorted fixture columns
      integer, parameter :: expected(4) = [3, 2, 1, 4]
      real(wp) :: pts(3, 4)
      integer :: perm(4), order(4), i1, i2, i3, i4, nperm, k
      !> Rendered input, result and expected orders
      character(len=:), allocatable :: in_text, got_text, want_text

      ! Non-vacuity: the fixture is a chain on which the pairwise tie is not transitive
      call check(error, abs(base(1, 2) - base(1, 1)) <= tie .and. abs(base(1, 3) - base(1, 2)) <= tie &
                 .and. abs(base(1, 3) - base(1, 1)) > tie, "fixture is no chain of near-ties in x")
      if (allocated(error)) return

      nperm = 0
      do i1 = 1, 4
         do i2 = 1, 4
            do i3 = 1, 4
               do i4 = 1, 4
                  perm = [i1, i2, i3, i4]
                  if (any(perm(2:) == i1) .or. any(perm(3:) == i2) .or. i3 == i4) cycle
                  nperm = nperm + 1
                  pts = base(:, perm)
                  call canonical_branch_order(pts, order)
                  do k = 1, 4
                     if (perm(order(k)) /= expected(k)) then
                        call format_ints(perm, in_text)
                        call format_ints(perm(order), got_text)
                        call format_ints(expected, want_text)
                        call test_failed(error, "input "//in_text//" sorts to "//got_text// &
                                         ", expected "//want_text//" for every input order")
                        return
                     end if
                  end do
               end do
            end do
         end do
      end do
      call check(error, nperm == 24, "not every permutation of the fixture was sorted")

   contains

      !> Integers as a bracketed list
      pure subroutine format_ints(v, s)
         !> Integers
         integer, intent(in) :: v(:)
         !> Bracketed list
         character(len=:), allocatable, intent(out) :: s
         character(len=12) :: buf
         integer :: i

         s = "["
         do i = 1, size(v)
            write (buf, "(i0)") v(i)
            s = s//trim(buf)
            if (i < size(v)) s = s//", "
         end do
         s = s//"]"
      end subroutine format_ints

   end subroutine test_finalize_canonical_order_chain

   !> Finalize on a sphere: duplicate points merge into one row whose CP Jacobian is
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
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
      type(newton_scratch_type) :: scratch
      type(drop_finalize_parameters_type) :: fcfg
      type(finalize_anchor_type) :: anchor
      type(finalize_scratch_type) :: fsc
      type(projection_rows_type) :: rows
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_hor, f_wleb
      type(branch_weight_type) :: branch_weight
      !> Reach gate open far beyond every point
      type(moist_math_smoothstep_sigmoid_type) :: open_reach
      type(moist_error_type), allocatable :: failed
      !> Sphere radius, anchor sphere radius, point
      real(wp) :: rs, point(3)
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
      point = rs*unit_vector(anchor%xyz)
      call new_context(ctx, verbosity=0, nthreads=1)

      ! The anchor and two nearby starts all reach the same point
      call batch%init(4, failed)
      call batch%begin_run(anchor%index, anchor%xyz, 1)
      call batch%append(anchor%xyz, 0.0_wp, icol)
      call batch%append(anchor%xyz + [1.0e-3_wp, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call batch%append(point, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
      if (.not. allocated(failed)) call polish_batch(ctx, batch, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "Newton/polish: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:3) == point_accepted), "all starts verified")
      if (allocated(error)) return

      call new_smoothstep_sigmoid(f_crit, 0.25_wp, 0.6_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_foc, 0.2_wp, 0.6_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_wleb, 1.0e-12_wp, 1.0e-10_wp)
      call new_smoothstep_sigmoid(f_hor, -2.0_wp, -1.0_wp)
      call new_smoothstep_sigmoid(open_reach, 2.0e2_wp, 1.0e2_wp)
      call branch_weight%init(0.2_wp, default_phi_alpha())
      call finalize_anchor(fcfg, default_phi_alpha(), batch, 1, 3, anchor, f_crit, f_foc, f_hor, &
                           f_wleb, open_reach, branch_weight, rows, nbranch, nrows, fsc, failed)
      if (allocated(failed)) then
         call test_failed(error, "finalize: "//failed%message)
         return
      end if
      call check(error, nbranch == 1 .and. nrows == 1 .and. rows%n == 1, "duplicates merge into one row")
      if (allocated(error)) return
      call check(error, rows%anchor(1) == 5 .and. rows%branch(1) == 1, &
                 "row indexing")
      if (allocated(error)) return
      call check(error, norm2(rows%xyz(:, 1) - point) < 1.0e-9_wp, "row point")
      if (allocated(error)) return
      call check(error, abs(rows%cpjac(1) - (rs/r_anchor)**2) < 1.0e-8_wp, "radial area ratio")
      if (allocated(error)) return
      call check(error, rows%w_f(1) > 0.0_wp .and. rows%w_f(1) <= 1.0_wp &
                 .and. abs(rows%wbranch(1) - rows%w_f(1)) < 1.0e-14_wp, &
                 "gate, and a single branch has it as branch weight")
      if (allocated(error)) return

      ! An anchor below the switching cutoff appends nothing
      anchor%f = 0.0_wp
      call finalize_anchor(fcfg, default_phi_alpha(), batch, 1, 3, anchor, f_crit, f_foc, f_hor, &
                           f_wleb, open_reach, branch_weight, rows, nbranch, nrows, fsc, failed)
      call check(error, nbranch == 1 .and. nrows == 0 .and. rows%n == 1, "rows below the cut are never appended")
   end subroutine test_finalize_sphere_anchor

   !> Finalize with two mirror branches stores each branch weight next to the row weight
   !> weight: every row is its single-branch weight times its branch weight,
   !> and the rows of the anchor sum to the single-branch weight
   !>
   !> - midpoint of a carbon pair; the single-branch references finalize each
   !>   column alone, where the branch weight is one
   !> - negative control: both branch weights are below one, so a lost fold
   !>   fails the per-row check
   !> - closed sibling: with the horizon switch on, the mirror point on the other
   !>   atom has a closed gate; it is no row and has no influence on the other
   !>
   !> @param[out] error Wrong row weight
   subroutine test_finalize_branch_fold(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Two-atom level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Points, solver and finalize state
      type(projection_batch_type) :: batch
      type(drop_solve_parameters_type) :: ncfg
      type(drop_polish_parameters_type) :: pcfg
      type(moist_context_type) :: ctx
      type(newton_scratch_type) :: scratch
      type(drop_finalize_parameters_type) :: fcfg
      type(finalize_anchor_type) :: anchor
      type(finalize_scratch_type) :: fsc
      type(projection_rows_type) :: rows, single, closed
      type(moist_math_smoothstep_sigmoid_type) :: f_crit, f_foc, f_hor, f_wleb
      type(branch_weight_type) :: branch_weight
      !> Reach gate open far beyond every point
      type(moist_math_smoothstep_sigmoid_type) :: open_reach
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
      call new_context(ctx, verbosity=0, nthreads=1)

      ! One start near each atom's near pole
      call batch%init(4, failed)
      call batch%begin_run(anchor%index, anchor%xyz, 1)
      call batch%append([-(h - rs), 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call batch%append([h - rs, 0.0_wp, 0.0_wp], 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, ncfg, pcfg, default_phi_alpha(), scratch, failed)
      if (.not. allocated(failed)) call polish_batch(ctx, batch, lsf, pcfg, ncfg, default_phi_alpha(), scratch, failed)
      if (allocated(failed)) then
         call test_failed(error, "Newton/polish: "//failed%message)
         return
      end if
      call check(error, all(batch%status(1:2) == point_accepted), "both mirror points verified")
      if (allocated(error)) return

      call new_smoothstep_sigmoid(f_crit, 0.25_wp, 0.6_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_foc, 0.2_wp, 0.6_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call new_smoothstep_sigmoid(f_wleb, 1.0e-12_wp, 1.0e-10_wp)
      call new_smoothstep_sigmoid(f_hor, -2.0_wp, -1.0_wp)
      call new_smoothstep_sigmoid(open_reach, 2.0e2_wp, 1.0e2_wp)
      call branch_weight%init(0.2_wp, default_phi_alpha())

      call finalize_anchor(fcfg, default_phi_alpha(), batch, 1, 2, anchor, f_crit, f_foc, f_hor, &
                           f_wleb, open_reach, branch_weight, rows, nbranch, nrows, fsc, failed)
      if (.not. allocated(failed)) then
         call finalize_anchor(fcfg, default_phi_alpha(), batch, 1, 1, anchor, f_crit, f_foc, f_hor, &
                              f_wleb, open_reach, branch_weight, single, nbranch, nrows, fsc, failed)
      end if
      if (.not. allocated(failed)) then
         call finalize_anchor(fcfg, default_phi_alpha(), batch, 2, 2, anchor, f_crit, f_foc, f_hor, &
                              f_wleb, open_reach, branch_weight, single, nbranch, nrows, fsc, failed)
      end if
      if (allocated(failed)) then
         call test_failed(error, "finalize: "//failed%message)
         return
      end if
      call check(error, rows%n == 2 .and. single%n == 2, "two kept rows, one per single-branch run")
      if (allocated(error)) return
      call check(error, all(abs(single%wbranch(1:2) - single%w_f(1:2)) < 1.0e-14_wp), &
                 "a single branch has its gate as branch weight")
      if (allocated(error)) return
      call check(error, all(single%w_f(1:2) == 1.0_wp), "the mirror points are fully open")
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
         w_ref = single%wleb(match)
         call check(error, abs(rows%wleb(r) - w_ref) < max(1.0e-15_wp, 1.0e-12_wp*abs(w_ref)), &
                    "row weight is the single-branch weight, without the branch weight")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(rows%wleb(1:2)*rows%wbranch(1:2)) - single%wleb(1)) &
                 < 1.0e-12_wp*single%wleb(1), &
                 "the branch-weighted rows sum to the anchor's single-branch weight")

      ! Closed sibling: under the horizon switch the mirror point on the other
      ! atom faces away from the anchor sphere; it is no row and the row left is
      ! the one of the open point finalized alone, bit for bit
      call new_smoothstep_sigmoid(f_hor, 0.0_wp, 0.15_wp, p_hi=0.8_wp, a_hi=1.6_wp, p_lo=0.8_wp, a_lo=1.6_wp)
      call closed%clear()
      call single%clear()
      call finalize_anchor(fcfg, default_phi_alpha(), batch, 1, 2, anchor, f_crit, f_foc, f_hor, &
                           f_wleb, open_reach, branch_weight, closed, nbranch, nrows, fsc, failed)
      if (.not. allocated(failed)) then
         call check(error, nbranch == 2 .and. nrows == 1 .and. closed%n == 1, &
                    "a closed sibling is a row or takes the open one out")
         if (allocated(error)) return
         call finalize_anchor(fcfg, default_phi_alpha(), batch, 1, 1, anchor, f_crit, f_foc, f_hor, &
                              f_wleb, open_reach, branch_weight, single, nbranch, nrows, fsc, failed)
      end if
      if (allocated(failed)) then
         call test_failed(error, "finalize: "//failed%message)
         return
      end if
      call check(error, single%n == 1 .and. all(closed%xyz(:, 1) == single%xyz(:, 1)) &
                 .and. closed%wbranch(1) == single%wbranch(1) .and. closed%wleb(1) == single%wleb(1) &
                 .and. closed%w_f(1) == single%w_f(1) .and. closed%branch(1) == 1, &
                 "the closed sibling changes the row of the open point")
      if (allocated(error)) return
      call check(error, closed%wbranch(1) == 1.0_wp, "the open point does not keep its whole weight")
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

   !> Objective weight of the compiled DROP defaults
   function default_phi_alpha() result(alpha)
      !> Objective weight
      real(wp) :: alpha
      !> Compiled defaults
      type(moist_cavity_drop_parameters_type) :: defaults

      alpha = defaults%phi_alpha
   end function default_phi_alpha

   !* ================================================================================= *!
   !*                                  Batch and Newton                                 *!
   !* ================================================================================= *!

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
      type(drop_solve_parameters_type) :: cfg
      type(drop_polish_parameters_type) :: pcfg
      type(newton_scratch_type) :: scratch
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs, anchor(3)
      integer :: icol
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      anchor = [80.0_wp, 0.0_wp, 0.0_wp]
      call batch%init(8, failed)
      call batch%begin_run(7, anchor, 1)
      call batch%append(anchor, 0.0_wp, icol)
      call newton_solve_batch(batch, lsf, cfg, pcfg, default_phi_alpha(), scratch, failed)
      call expect_error(error, failed, &
                        "Projection Newton: anchor 7 lies beyond the level set reach, no atom screens it")
   end subroutine test_bad_newton_unreached_anchor

   !* ================================================================================= *!
   !*                              Surface mesh seeder                                  *!
   !* ================================================================================= *!

   !> Surface mesh build input of a cluster cavity: every anchor seeded
   !>
   !> @param[in]  cavity   Anchor-phase cavity
   !> @param[in]  nthreads Build threads
   !> @param[out] geometry Build input
   subroutine cluster_build_input(cavity, nthreads, geometry)
      !> Anchor-phase cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Build threads
      integer, intent(in) :: nthreads
      !> Build input
      type(seed_build_type), intent(out) :: geometry

      geometry%centers = cavity%mol%xyz
      geometry%radii = cavity%radii
      geometry%anchors = cavity%anchor%xyz(:, 1:cavity%anchor%ngrid)
      geometry%reach = cavity%param%branch%reach_off
      geometry%nthreads = nthreads
   end subroutine cluster_build_input

   !> Built level-4 surface mesh seeder of a cluster cavity
   !>
   !> @param[in]  cavity   Anchor-phase cavity
   !> @param[in]  nthreads Build threads
   !> @param[out] seeder   Built seeder
   !> @param[out] error    Native failure
   subroutine build_cluster_rmesh(cavity, nthreads, seeder, error)
      !> Anchor-phase cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Build threads
      integer, intent(in) :: nthreads
      !> Built seeder
      type(projection_seeder_type), intent(out) :: seeder
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Level-4 settings
      type(moist_cavity_drop_parameters_type) :: param
      !> Build input
      type(seed_build_type) :: geometry

      param = cavity%param
      param%proj_level = 4
      call param%compute_derived(error)
      if (allocated(error)) return
      call init_mesh(seeder, error, sample_size=param%seed%seed_size, refine_levels=param%seed%refine_levels)
      if (allocated(error)) return
      call cluster_build_input(cavity, nthreads, geometry)
      call seeder%build(cavity%lsf_model, geometry, error)
   end subroutine build_cluster_rmesh

   !> The built mesh does not depend on the build's thread count
   !>
   !> @param[out] error Meshes differ
   subroutine test_rmesh_thread_independent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Shared serial context
      type(moist_context_type), target :: ctx
      !> Anchor-phase cavity
      type(cavity_type_drop), allocatable :: cavity
      !> Fixture geometry
      type(structure_type) :: mol
      !> Reference level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Serial and threaded builds
      type(projection_seeder_type) :: serial, threaded
      !> Native failure
      type(moist_error_type), allocatable :: native_error

      call new(mol, [6, 6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                        5.2_wp, 0.3_wp, 0.0_wp, &
                                        2.4_wp, 4.6_wp, 0.7_wp], [3, 3]))
      call new_context(ctx, verbosity=0, nthreads=1)
      allocate (cavity)
      call build_cluster_anchors(cavity, ctx, mol, 302, 0.8_wp, lsf, native_error)
      if (.not. allocated(native_error)) call build_cluster_rmesh(cavity, 1, serial, native_error)
      if (.not. allocated(native_error)) call build_cluster_rmesh(cavity, 3, threaded, native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, serial%nsample, threaded%nsample, message="Sample counts differ")
      if (allocated(error)) return
      call check(error, serial%nsample > 0, message="The mesh holds no sample")
      if (allocated(error)) return
      call check(error, all(serial%y == threaded%y) .and. all(serial%normal == threaded%normal) &
                 .and. all(serial%form == threaded%form), message="Samples differ between thread counts")
      if (allocated(error)) return
      call check(error, serial%cover == threaded%cover, message="Covering radii differ")
      deallocate (cavity)
      call ctx%delete()
   end subroutine test_rmesh_thread_independent

   !> The level-3 seed filters cut the seeds and keep every row of the unfiltered search
   !>
   !> - two meshes of the level-3 spacing: the level-3 one with its filters, and
   !>   one without, as level 4 keeps every sample
   !>
   !> @param[out] error Filter did not reduce seeds or lost a row
   subroutine test_rmesh_level_filters(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Shared serial context
      type(moist_context_type), target :: ctx
      !> Anchor-phase cavity
      type(cavity_type_drop), allocatable :: cavity
      !> Fixture geometry
      type(structure_type) :: mol
      !> Reference level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Level-3 meshes with and without the filters, and the workers for level 3 and 4
      type(projection_seeder_type) :: seeder, unfiltered
      type(seed_build_type) :: geometry
      type(projection_worker_type) :: w3, w4
      type(moist_cavity_drop_parameters_type) :: param
      type(moist_error_type), allocatable :: native_error
      integer, allocatable :: nrows(:), slots(:)
      logical, allocatable :: lost(:)
      integer :: i, n, ia, ib, matches

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                             0.00_wp, 4.21_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, 4.22_wp, &
                                             0.00_wp, -4.18_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, -4.15_wp, &
                                             0.02_wp, 0.10_wp, -0.20_wp], [3, 5])*aatoau)
      call new_context(ctx, verbosity=0, nthreads=1)
      allocate (cavity)
      call build_cluster_anchors(cavity, ctx, mol, 302, 1.0_wp, lsf, native_error)
      if (.not. allocated(native_error)) then
         param = cavity%param
         param%proj_level = 3
         call param%compute_derived(native_error)
      end if
      if (.not. allocated(native_error)) then
         call init_mesh(seeder, native_error, sample_size=param%seed%seed_size, refine_levels=param%seed%refine_levels, &
                          filter_cone=param%seed%filter_cone, filter_grazing=param%seed%filter_grazing, &
                          filter_back_facing=param%seed%filter_back_facing, filter_focal=param%seed%filter_focal, &
                          cone_scale=param%seed%cone_scale)
      end if
      if (.not. allocated(native_error)) then
         call init_mesh(unfiltered, native_error, sample_size=param%seed%seed_size, &
                              refine_levels=param%seed%refine_levels)
      end if
      if (.not. allocated(native_error)) then
         call cluster_build_input(cavity, 1, geometry)
         call seeder%build(cavity%lsf_model, geometry, native_error)
      end if
      if (.not. allocated(native_error)) call unfiltered%build(cavity%lsf_model, geometry, native_error)
      if (.not. allocated(native_error)) call w3%init(param, cavity%lsf_model, cavity%mol%xyz, native_error)
      if (.not. allocated(native_error)) then
         param%proj_level = 4
         call param%compute_derived(native_error)
      end if
      if (.not. allocated(native_error)) call w4%init(param, cavity%lsf_model, cavity%mol%xyz, native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if

      w3%finalize%wleb_cut = -1.0_wp
      w4%finalize%wleb_cut = -1.0_wp
      n = cavity%anchor%ngrid
      allocate (nrows(n), lost(n))
      slots = [(i, i=1, n)]
      call w3%process_chunk(cavity%ctx, slots, cavity%anchor%xyz(:, 1:n), cavity%anchor%owner(1:n), &
                            cavity%anchor%wleb0(1:n), cavity%anchor%f0(1:n), nrows, lost, native_error, &
                            seeder=seeder)
      if (.not. allocated(native_error)) then
         call w4%process_chunk(cavity%ctx, slots, cavity%anchor%xyz(:, 1:n), cavity%anchor%owner(1:n), &
                               cavity%anchor%wleb0(1:n), cavity%anchor%f0(1:n), nrows, lost, native_error, &
                               seeder=unfiltered)
      end if
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if

      call check(error, w3%n_seeds > 0 .and. 2*w3%n_seeds < w4%n_seeds, &
                 message="The filters did not cut level-3 seeds by >50% relative to the unfiltered search")
      if (allocated(error)) return
      call check_projection_geometry(cavity, lsf, w3%rows, error)
      if (allocated(error)) return
      do ia = 1, w4%rows%n
         if (.not. relevant_row(cavity, w4%rows, ia)) cycle
         matches = 0
         do ib = 1, w3%rows%n
            if (w3%rows%anchor(ib) /= w4%rows%anchor(ia)) cycle
            if (.not. relevant_row(cavity, w3%rows, ib)) cycle
            if (norm2(w4%rows%xyz(:, ia) - w3%rows%xyz(:, ib)) <= 1.0e-6_wp) matches = matches + 1
         end do
         call check(error, matches, 1, message="The filters dropped a relevant row of the unfiltered search")
         if (allocated(error)) return
      end do
      deallocate (cavity)
      call ctx%delete()
   end subroutine test_rmesh_level_filters

   !> Every anchor of a cluster cavity through one worker, the raw rows kept
   !>
   !> - the worker builds its own seeder for the chunk unless `seeder` is given
   !>
   !> @param[in]    cavity Anchor-phase cavity
   !> @param[in]    param  Projection settings, derived
   !> @param[inout] worker Worker, set up here
   !> @param[out]   nrows  Rows per anchor
   !> @param[out]   error  Native failure
   !> @param[in]    seeder Built seeder to search with
   subroutine project_cluster_with(cavity, param, worker, nrows, error, seeder)
      !> Anchor-phase cavity
      type(cavity_type_drop), intent(in) :: cavity
      !> Projection settings
      type(moist_cavity_drop_parameters_type), intent(in) :: param
      !> Worker
      type(projection_worker_type), intent(inout) :: worker
      !> Rows per anchor
      integer, allocatable, intent(out) :: nrows(:)
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Built seeder
      type(projection_seeder_type), intent(in), optional :: seeder
      logical, allocatable :: lost(:)
      integer :: n, k

      n = cavity%anchor%ngrid
      allocate (nrows(n), lost(n))
      call worker%init(param, cavity%lsf_model, cavity%mol%xyz, error)
      if (allocated(error)) return
      worker%finalize%wleb_cut = -1.0_wp
      call worker%process_chunk(cavity%ctx, [(k, k=1, n)], cavity%anchor%xyz(:, 1:n), cavity%anchor%owner(1:n), &
                                cavity%anchor%wleb0(1:n), cavity%anchor%f0(1:n), nrows, lost, error, seeder=seeder)
   end subroutine project_cluster_with

   !> Whether every relevant row of anchor `aid` in `a` has exactly one relevant
   !> row of that anchor in `b` within 1e-6 Bohr
   !>
   !> @param[in] cavity Settings and anchors
   !> @param[in] a      Rows to match
   !> @param[in] b      Candidate partners
   !> @param[in] aid    Anchor
   function rows_covered(cavity, a, b, aid) result(covered)
      !> Settings and anchors
      type(cavity_type_drop), intent(in) :: cavity
      !> Rows to match and candidate partners
      type(projection_rows_type), intent(in) :: a, b
      !> Anchor
      integer, intent(in) :: aid
      !> Every relevant row matched once
      logical :: covered
      integer :: ia, ib, matches

      covered = .true.
      do ia = 1, a%n
         if (a%anchor(ia) /= aid) cycle
         if (.not. relevant_row(cavity, a, ia)) cycle
         matches = 0
         do ib = 1, b%n
            if (b%anchor(ib) /= aid) cycle
            if (.not. relevant_row(cavity, b, ib)) cycle
            if (norm2(a%xyz(:, ia) - b%xyz(:, ib)) <= 1.0e-6_wp) matches = matches + 1
         end do
         covered = covered .and. matches == 1
      end do
   end function rows_covered

   !> Cross cluster of [[test_rmesh_level_filters]] with a soft level set, anchors at 302 points per atom
   !>
   !> @param[out]   cavity Anchor-phase cavity
   !> @param[inout] ctx    Borrowed context
   !> @param[out]   lsf    Level set of the geometry
   !> @param[out]   error  Native failure
   subroutine soft_cross_cluster(cavity, ctx, lsf, error)
      !> Anchor-phase cavity
      type(cavity_type_drop), intent(out) :: cavity
      !> Borrowed context
      type(moist_context_type), intent(inout), target :: ctx
      !> Level set of the geometry
      type(moist_cavity_drop_lsf_svdw_type), intent(out) :: lsf
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                             0.00_wp, 4.21_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, 4.22_wp, &
                                             0.00_wp, -4.18_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, -4.15_wp, &
                                             0.02_wp, 0.10_wp, -0.20_wp], [3, 5])*aatoau)
      call build_cluster_anchors(cavity, ctx, mol, 302, 5.5_wp, lsf, error)
   end subroutine soft_cross_cluster

   !> Level 2 finalizes a close fully open local point and otherwise searches as level 3
   !>
   !> - the cross cluster at level 2, and at level 3 with level 2's mesh spacing
   !>   and seed filters, so only the routing differs; each worker builds its own
   !>   mesh from its chunk's search balls
   !> - level 3 searches the anchors level 2 searches and those level 2
   !>   finalizes from a fully open local point within 0.5 Bohr
   !> - an anchor whose nearest point lies beyond 0.5 Bohr keeps the same
   !>   relevant rows at both levels; a nearer one keeps a subset at level 2
   !>
   !> @param[out] error Level 2 routed or searched differently
   subroutine test_level2_selective_rmesh(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Shared serial context
      type(moist_context_type), target :: ctx
      !> Anchor-phase cavity
      type(cavity_type_drop), allocatable :: cavity
      !> Reference level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Level-2 and level-3 workers and their settings
      type(projection_worker_type) :: w2, w3
      type(moist_cavity_drop_parameters_type) :: param2, param3
      type(moist_error_type), allocatable :: native_error
      integer, allocatable :: nrows2(:), nrows3(:), first(:), last(:)
      integer :: aid
      logical :: same
      character(len=64) :: message
      !> Level 2's distance below which a fully open local point is finalized (Bohr)
      real(wp), parameter :: close_pdist = 0.5_wp

      call new_context(ctx, verbosity=0, nthreads=1)
      allocate (cavity)
      call soft_cross_cluster(cavity, ctx, lsf, native_error)
      if (.not. allocated(native_error)) then
         param2 = cavity%param
         param2%proj_level = 2
         call param2%compute_derived(native_error)
      end if
      if (.not. allocated(native_error)) then
         param3 = cavity%param
         param3%proj_level = 3
         call param3%compute_derived(native_error)
      end if
      if (.not. allocated(native_error)) then
         param3%seed%seed_size = param2%seed%seed_size
         param3%seed%cone_scale = param2%seed%cone_scale
         param3%seed%filter_grazing = param2%seed%filter_grazing
         param3%seed%filter_back_facing = param2%seed%filter_back_facing
         param3%seed%filter_focal = param2%seed%filter_focal
         call project_cluster_with(cavity, param2, w2, nrows2, native_error)
      end if
      if (.not. allocated(native_error)) call project_cluster_with(cavity, param3, w3, nrows3, native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if

      call check(error, w2%n_window_skipped > 0 .and. w3%n_window_skipped == 0, &
                 message="Level 2 finalized no close fully open local point")
      if (allocated(error)) return
      call check(error, w3%n_seeded, w2%n_seeded + w2%n_window_skipped, &
                 message="Level 3 searched other anchors than level 2 and its close ones")
      if (allocated(error)) return
      call check(error, count(nrows2 > 1) > 0, message="Level 2 found no branched anchors on cross cluster")
      if (allocated(error)) return
      call check_projection_geometry(cavity, lsf, w2%rows, error)
      if (allocated(error)) return

      call row_bounds(w3%rows, cavity%anchor%ngrid, first, last)
      do aid = 1, cavity%anchor%ngrid
         same = rows_covered(cavity, w2%rows, w3%rows, aid)
         if (first(aid) > 0) then
            if (minval(w3%rows%pdist(first(aid):last(aid))) > close_pdist) then
               same = same .and. rows_covered(cavity, w3%rows, w2%rows, aid)
            end if
         end if
         write (message, "(a,i0)") "Level 2 and level 3 rows differ, anchor ", cavity%anchor%id(aid)
         call check(error, same, message=trim(message))
         if (allocated(error)) return
      end do
      deallocate (cavity)
      call ctx%delete()
   end subroutine test_level2_selective_rmesh

   !> Level 2's seed filters cut seeds and keep every relevant row
   !>
   !> - the cross cluster of [[test_level2_selective_rmesh]] at level 2, at
   !>   level 2 without its back-facing or its grazing filter, with the focal
   !>   filter alone, without any filter, with a wider tangent cone (2.0), and
   !>   with that wider cone alone
   !> - each filter of level 2 leaves fewer seeds to solve; the focal filter is
   !>   compared without the others, which drop its samples first here
   !> - level 2 keeps the relevant rows of the run with the wider cone alone
   !>
   !> @param[out] error A filter cut no seed, or a relevant row was lost
   subroutine test_level2_seed_filters(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Shared serial context
      type(moist_context_type), target :: ctx
      !> Anchor-phase cavity
      type(cavity_type_drop), allocatable :: cavity
      !> Reference level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      !> Workers of level 2, without its back-facing or grazing filter, with the
      !> focal filter alone, without any filter, with the wider cone, and with
      !> the wider cone alone
      type(projection_worker_type) :: w2, w_back, w_graze, w_focal, w_none, w_cone, w_cone_only
      type(moist_cavity_drop_parameters_type) :: param, variant
      type(moist_error_type), allocatable :: native_error
      integer, allocatable :: nrows(:)
      integer :: aid
      character(len=64) :: message

      call new_context(ctx, verbosity=0, nthreads=1)
      allocate (cavity)
      call soft_cross_cluster(cavity, ctx, lsf, native_error)
      if (.not. allocated(native_error)) then
         param = cavity%param
         param%proj_level = 2
         call param%compute_derived(native_error)
      end if
      if (.not. allocated(native_error)) then
         call check(error, param%seed%filter_cone .and. param%seed%filter_grazing .and. &
                    param%seed%filter_back_facing .and. param%seed%filter_focal .and. &
                    param%seed%cone_scale == 1.5_wp, message="Level 2 settings without its filters")
         if (allocated(error)) return
         call project_cluster_with(cavity, param, w2, nrows, native_error)
      end if
      if (.not. allocated(native_error)) then
         variant = param
         variant%seed%filter_back_facing = .false.
         call project_cluster_with(cavity, variant, w_back, nrows, native_error)
      end if
      if (.not. allocated(native_error)) then
         variant = param
         variant%seed%filter_grazing = .false.
         call project_cluster_with(cavity, variant, w_graze, nrows, native_error)
      end if
      if (.not. allocated(native_error)) then
         variant = param
         variant%seed%filter_cone = .false.
         variant%seed%filter_back_facing = .false.
         variant%seed%filter_grazing = .false.
         call project_cluster_with(cavity, variant, w_focal, nrows, native_error)
      end if
      if (.not. allocated(native_error)) then
         variant%seed%filter_focal = .false.
         call project_cluster_with(cavity, variant, w_none, nrows, native_error)
      end if
      if (.not. allocated(native_error)) then
         variant = param
         variant%seed%cone_scale = 2.0_wp
         call project_cluster_with(cavity, variant, w_cone, nrows, native_error)
      end if
      if (.not. allocated(native_error)) then
         variant%seed%filter_back_facing = .false.
         variant%seed%filter_grazing = .false.
         variant%seed%filter_focal = .false.
         call project_cluster_with(cavity, variant, w_cone_only, nrows, native_error)
      end if
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if

      call check(error, w2%n_seeds > 0 .and. w2%n_seeds < w_back%n_seeds, &
                 message="The back-facing filter cut no seed")
      if (allocated(error)) return
      call check(error, w2%n_seeds < w_graze%n_seeds, message="The grazing filter cut no seed")
      if (allocated(error)) return
      call check(error, w_focal%n_seeds < w_none%n_seeds, message="The focal filter cut no seed")
      if (allocated(error)) return
      call check(error, w2%n_seeds < w_cone%n_seeds, message="The narrower tangent cone cut no seed")
      if (allocated(error)) return
      do aid = 1, cavity%anchor%ngrid
         write (message, "(a,i0)") "The level-2 filters lost a relevant row, anchor ", cavity%anchor%id(aid)
         call check(error, rows_covered(cavity, w_cone_only%rows, w2%rows, aid), message=trim(message))
         if (allocated(error)) return
      end do
      deallocate (cavity)
      call ctx%delete()
   end subroutine test_level2_seed_filters

   !> A shell query and the ball inside it give the ball's seeds
   !>
   !> - for anchors around one sphere: the seeds of the ball `r`, as a multiset,
   !>   are those of the ball `r_in` and of the shell `(r_in, r]` together, and
   !>   no seed of the shell is one of the inner ball
   !> - one atom, so the level set's stationarity certificate takes the same
   !>   atom bounds for every ball and decides each bin alike
   !>
   !> @param[out] error Shell and inner ball do not partition the ball
   subroutine test_rmesh_shell_query(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(projection_seeder_type) :: seeder
      type(seeder_scratch_type) :: scratch
      type(seed_build_type) :: geometry
      type(seed_request_type) :: request
      type(moist_error_type), allocatable :: failed
      !> Seeds of the ball, the inner ball and the shell
      real(wp), allocatable :: ball(:, :), inner(:, :), shell(:, :)
      real(wp) :: rs
      logical :: shell_flag
      integer :: ia
      !> Ball and inner ball radii (Bohr)
      real(wp), parameter :: r_out = 1.4_wp, r_in = 0.6_wp

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      call init_mesh(seeder, failed, sample_size=0.2_wp, refine_levels=1)
      if (allocated(failed)) then
         call test_failed(error, "init: "//failed%message)
         return
      end if
      call single_atom_input(geometry)
      geometry%anchors = reshape([rs, 0.0_wp, 0.0_wp, 0.0_wp, rs + 0.4_wp, 0.0_wp, &
                                  0.0_wp, 0.0_wp, rs - 0.3_wp], [3, 3])
      call seeder%build(lsf, geometry, failed)
      if (allocated(failed)) then
         call test_failed(error, "build: "//failed%message)
         return
      end if
      do ia = 1, size(geometry%anchors, 2)
         request%anchor = geometry%anchors(:, ia)
         request%owner = 1
         request%pdist_start = r_out
         request%pdist_inner = -1.0_wp
         call query(ball, shell_flag)
         if (allocated(error)) return
         call check(error, .not. shell_flag, "a ball query is not a shell")
         if (allocated(error)) return
         request%pdist_start = r_in
         call query(inner, shell_flag)
         if (allocated(error)) return
         request%pdist_start = r_out
         request%pdist_inner = r_in
         call query(shell, shell_flag)
         if (allocated(error)) return
         call check(error, shell_flag, "a shell query reports its shell")
         if (allocated(error)) return
         call check(error, size(inner, 2) > 0 .and. size(shell, 2) > 0, "the fixture has seeds inside and outside")
         if (allocated(error)) return
         call check(error, size(ball, 2), size(inner, 2) + size(shell, 2), "the ball holds the inner ball and the shell")
         if (allocated(error)) return
         call check(error, same_multiset(ball, reshape([inner, shell], [3, size(inner, 2) + size(shell, 2)])), &
                    "the ball's seeds are the inner ball's and the shell's")
         if (allocated(error)) return
         call check(error, disjoint(inner, shell), "a seed of the inner ball is in the shell")
         if (allocated(error)) return
      end do

   contains

      !> Seeds of the current request
      !>
      !> @param[out] points    Seeds (3, n)
      !> @param[out] filtered  The scratch's `shell_filtered`
      subroutine query(points, filtered)
         real(wp), allocatable, intent(out) :: points(:, :)
         logical, intent(out) :: filtered

         call seeder%seed(scratch, lsf, request, failed)
         if (allocated(failed)) then
            call test_failed(error, "seed: "//failed%message)
            allocate (points(3, 0))
            filtered = .false.
            return
         end if
         points = scratch%points(:, 1:scratch%n_points)
         filtered = scratch%shell_filtered
      end subroutine query

      !> Whether two point lists hold the same points, bit for bit, with multiplicity
      !>
      !> @param[in] a First list
      !> @param[in] b Second list, as long as `a`
      logical function same_multiset(a, b)
         real(wp), intent(in) :: a(:, :), b(:, :)

         logical, allocatable :: used(:)
         integer :: i, j

         same_multiset = .false.
         allocate (used(size(b, 2)), source=.false.)
         do i = 1, size(a, 2)
            do j = 1, size(b, 2)
               if (used(j)) cycle
               if (all(a(:, i) == b(:, j))) exit
            end do
            if (j > size(b, 2)) return
            used(j) = .true.
         end do
         same_multiset = all(used)
      end function same_multiset

      !> Whether no point of `a` equals one of `b`
      !>
      !> @param[in] a First list
      !> @param[in] b Second list
      logical function disjoint(a, b)
         real(wp), intent(in) :: a(:, :), b(:, :)

         integer :: i, j

         disjoint = .false.
         do i = 1, size(a, 2)
            do j = 1, size(b, 2)
               if (all(a(:, i) == b(:, j))) return
            end do
         end do
         disjoint = .true.
      end function disjoint

   end subroutine test_rmesh_shell_query

   !> SSES geometry of three atoms with anchors on their spheres
   !>
   !> @param[out] lsf      SSES level set, updated for the geometry
   !> @param[out] mol      Geometry
   !> @param[out] geometry Build input: six anchors per atom, balls of 3 Bohr
   subroutine sses_cluster(lsf, mol, geometry)
      type(moist_cavity_drop_lsf_sses_type), intent(out) :: lsf
      type(structure_type), intent(out) :: mol
      type(seed_build_type), intent(out) :: geometry

      real(wp), parameter :: radius = 3.2_wp
      real(wp), parameter :: dirs(3, 6) = reshape([1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp, &
                                                   0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, &
                                                   0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, -1.0_wp], [3, 6])
      integer :: i, d

      call new(mol, [6, 6, 6], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                        5.2_wp, 0.3_wp, 0.0_wp, &
                                        2.4_wp, 4.6_wp, 0.7_wp], [3, 3]))
      call lsf%new()
      call lsf%update(mol, [radius, radius, radius])
      geometry%centers = mol%xyz
      geometry%radii = [radius, radius, radius]
      allocate (geometry%anchors(3, 18))
      do i = 1, 3
         do d = 1, 6
            geometry%anchors(:, 6*(i - 1) + d) = mol%xyz(:, i) + radius*dirs(:, d)
         end do
      end do
      geometry%reach = 5.67_wp
   end subroutine sses_cluster

   !> An SSES mesh does not depend on the build's thread count
   !>
   !> - SSES warm-starts its moment solve from the last point a clone evaluated;
   !>   every build slice starts without that history, so builds at 1 and 3
   !>   threads, and two at 3, give the same samples bit for bit
   !>
   !> @param[out] error Meshes differ
   subroutine test_rmesh_sses_thread_independent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      type(seed_build_type) :: geometry
      type(projection_seeder_type) :: serial, threaded
      type(moist_error_type), allocatable :: failed
      integer :: run

      call sses_cluster(lsf, mol, geometry)
      call init_mesh(serial, failed, sample_size=0.5_wp, refine_levels=0)
      if (.not. allocated(failed)) call serial%build(lsf, geometry, failed)
      if (allocated(failed)) then
         call test_failed(error, "serial build: "//failed%message)
         return
      end if
      call check(error, serial%nsample > 0, message="The mesh holds no sample")
      if (allocated(error)) return
      do run = 1, 2
         geometry%nthreads = 3
         call init_mesh(threaded, failed, sample_size=0.5_wp, refine_levels=0)
         if (.not. allocated(failed)) call threaded%build(lsf, geometry, failed)
         if (allocated(failed)) then
            call test_failed(error, "threaded build: "//failed%message)
            return
         end if
         call check(error, serial%nsample, threaded%nsample, message="Sample counts differ")
         if (allocated(error)) return
         call check(error, all(serial%y == threaded%y) .and. all(serial%normal == threaded%normal) &
                    .and. all(serial%form == threaded%form) .and. serial%cover == threaded%cover, &
                    message="Samples differ between thread counts")
         if (allocated(error)) return
      end do
   end subroutine test_rmesh_sses_thread_independent

   !> An SSES projection does not depend on the thread count
   !>
   !> - the level-3 projector of a three-atom SSES cavity at 1 and 3 threads,
   !>   twice at 3, in chunks of 8 anchors: the same rows bit for bit
   !>
   !> @param[out] error Rows differ
   subroutine test_sses_projection_thread_independent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_sses_type) :: lsf
      type(structure_type) :: mol
      type(seed_build_type) :: geometry
      type(moist_error_type), allocatable :: failed
      real(wp), allocatable :: ref(:, :), rows(:, :)
      integer :: run
      integer, parameter :: threads(3) = [1, 3, 3]

      call sses_cluster(lsf, mol, geometry)
      do run = 1, size(threads)
         call project(threads(run), rows)
         if (allocated(error)) return
         if (run == 1) then
            call check(error, size(rows, 2) > 0, message="The projection has no row")
            if (allocated(error)) return
            ref = rows
            cycle
         end if
         call check(error, size(rows, 2), size(ref, 2), message="Row counts differ between thread counts")
         if (allocated(error)) return
         call check(error, all(rows == ref), message="Rows differ between thread counts")
         if (allocated(error)) return
      end do

   contains

      !> Rows of the level-3 projector at `nthreads`
      !>
      !> @param[in]  nthreads Threads of the run context
      !> @param[out] xyz      Row positions (3, n)
      subroutine project(nthreads, xyz)
         integer, intent(in) :: nthreads
         real(wp), allocatable, intent(out) :: xyz(:, :)

         type(moist_context_type), target :: ctx
         type(cavity_type_drop), allocatable :: cavity
         type(moist_cavity_drop_parameters_type) :: param
         integer :: nremoved

         call new_context(ctx, verbosity=0, nthreads=nthreads)
         param%num_leb = 302
         param%proj_level = 3
         ! Many small chunks, which the threads take in any order
         param%batch_size = 8
         allocate (cavity)
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=lsf, param=param, &
                              error=failed, ctx=ctx)
         if (.not. allocated(failed)) then
            cavity%nsph = mol%nat
            call cavity%radius_model%update(mol, failed)
         end if
         if (.not. allocated(failed)) then
            cavity%radii = cavity%radius_model%f0
            cavity%mol = mol
            cavity%sphxyz = mol%xyz
            call cavity%lsf_model%update(mol, cavity%radii)
            call cavity%ensure_lebedev_cache(failed)
         end if
         if (.not. allocated(failed)) then
            call cavity%anchor%fill(cavity%sphxyz, cavity%radii, cavity%ang_grid, cavity%ang_weight)
            call cavity%iswig%update(mol, cavity%radii, wleb_max=maxval(cavity%anchor%wleb0))
            call cavity%anchor%compute_switching(cavity%iswig, ctx%get_num_threads())
            call cavity%anchor%filter(cavity%param%finalize%wleb_cut, nremoved)
            call cavity%project_all_points(failed)
         end if
         if (allocated(failed)) then
            call test_failed(error, failed%message)
            allocate (xyz(3, 0))
            return
         end if
         xyz = cavity%xyz(:, 1:cavity%ngrid)
         deallocate (cavity)
         call ctx%delete()
      end subroutine project

   end subroutine test_sses_projection_thread_independent

   !> One-atom build input around the origin
   !>
   !> @param[out] geometry Build input
   subroutine single_atom_input(geometry)
      !> Build input
      type(seed_build_type), intent(out) :: geometry

      allocate (geometry%centers(3, 1), source=0.0_wp)
      geometry%radii = [3.0_wp]
      geometry%anchors = reshape([3.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      geometry%reach = 5.67_wp
   end subroutine single_atom_input

   !> Initialize a mesh seeder from the default seed settings and the given overrides
   !>
   !> @param[inout] seeder             Seeder
   !> @param[out]   error              Refused settings
   !> @param[in]    sample_size        Sample spacing (Bohr)
   !> @param[in]    refine_levels      Levels below the sample spacing
   !> @param[in]    filter_cone        Tangent-cone filter
   !> @param[in]    filter_grazing     Grazing filter
   !> @param[in]    filter_back_facing Back-facing filter
   !> @param[in]    filter_focal       Focal filter
   !> @param[in]    cone_scale         Tangent-cone scale
   subroutine init_mesh(seeder, error, sample_size, refine_levels, filter_cone, filter_grazing, &
                        filter_back_facing, filter_focal, cone_scale)
      type(projection_seeder_type), intent(inout) :: seeder
      type(moist_error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: sample_size
      integer, intent(in), optional :: refine_levels
      logical, intent(in), optional :: filter_cone
      logical, intent(in), optional :: filter_grazing
      logical, intent(in), optional :: filter_back_facing
      logical, intent(in), optional :: filter_focal
      real(wp), intent(in), optional :: cone_scale

      type(drop_seed_parameters_type) :: seed

      if (present(sample_size)) seed%seed_size = sample_size
      if (present(refine_levels)) seed%refine_levels = refine_levels
      if (present(filter_cone)) seed%filter_cone = filter_cone
      if (present(filter_grazing)) seed%filter_grazing = filter_grazing
      if (present(filter_back_facing)) seed%filter_back_facing = filter_back_facing
      if (present(filter_focal)) seed%filter_focal = filter_focal
      if (present(cone_scale)) seed%cone_scale = cone_scale
      call seeder%init(seed, error)
   end subroutine init_mesh

   !> A nonpositive sample size is refused
   subroutine test_bad_rmesh_sample_size(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_error_type), allocatable :: failed

      call init_mesh(seeder, failed, sample_size=0.0_wp)
      call expect_error(error, failed, "Surface mesh seeder: sample size must be positive")
   end subroutine test_bad_rmesh_sample_size

   !> A negative refinement depth is refused
   subroutine test_bad_rmesh_refine_levels(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_error_type), allocatable :: failed

      call init_mesh(seeder, failed, refine_levels=-1)
      call expect_error(error, failed, "Surface mesh seeder: refine_levels must not be negative")
   end subroutine test_bad_rmesh_refine_levels

   !> A build before `init` is refused
   subroutine test_bad_rmesh_build_before_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_build_type) :: geometry
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      call single_atom_input(geometry)
      call seeder%build(lsf, geometry, failed)
      call expect_error(error, failed, "Surface mesh seeder: build before init")
   end subroutine test_bad_rmesh_build_before_init

   !> Seeding before the build is refused
   subroutine test_bad_rmesh_seed_before_build(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(seeder_scratch_type) :: scratch
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_request_type) :: request
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      call init_mesh(seeder, failed)
      if (allocated(failed)) return
      request%anchor = [rs, 0.0_wp, 0.0_wp]
      request%pdist_start = 1.0_wp
      call seeder%seed(scratch, lsf, request, failed)
      call expect_error(error, failed, "Surface mesh seeder: seed before build")
   end subroutine test_bad_rmesh_seed_before_build

   !> Radii that do not match the centers are refused
   subroutine test_bad_rmesh_radii(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_build_type) :: geometry
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      call init_mesh(seeder, failed)
      if (allocated(failed)) return
      call single_atom_input(geometry)
      geometry%radii = [3.0_wp, 3.0_wp]
      call seeder%build(lsf, geometry, failed)
      call expect_error(error, failed, "Surface mesh seeder: one radius per atom center required")
   end subroutine test_bad_rmesh_radii

   !> A nonpositive reach is refused
   subroutine test_bad_rmesh_reach(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_build_type) :: geometry
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      call init_mesh(seeder, failed)
      if (allocated(failed)) return
      call single_atom_input(geometry)
      geometry%reach = 0.0_wp
      call seeder%build(lsf, geometry, failed)
      call expect_error(error, failed, "Surface mesh seeder: reach must be positive")
   end subroutine test_bad_rmesh_reach

   !> A region beyond the lattice key range is refused
   subroutine test_bad_rmesh_lattice_range(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_build_type) :: geometry
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      call init_mesh(seeder, failed)
      if (allocated(failed)) return
      geometry%centers = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0e8_wp, 0.0_wp, 0.0_wp], [3, 2])
      geometry%radii = [3.0_wp, 3.0_wp]
      geometry%anchors = reshape([3.0_wp, 0.0_wp, 0.0_wp], [3, 1])
      geometry%reach = 5.67_wp
      call seeder%build(lsf, geometry, failed)
      call expect_error(error, failed, "Surface mesh seeder: the region does not fit the lattice keys; "// &
                        "raise the sample size or split the geometry")
   end subroutine test_bad_rmesh_lattice_range

   !> An infinite sample size is refused
   subroutine test_bad_rmesh_sample_size_infinite(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_error_type), allocatable :: failed

      call init_mesh(seeder, failed, sample_size=ieee_value(1.0_wp, ieee_positive_inf))
      call expect_error(error, failed, "Surface mesh seeder: sample size must be finite")
   end subroutine test_bad_rmesh_sample_size_infinite

   !> A negative tangent-cone scale is refused
   subroutine test_bad_rmesh_cone_scale(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_error_type), allocatable :: failed

      call init_mesh(seeder, failed, cone_scale=-1.5_wp)
      call expect_error(error, failed, "Surface mesh seeder: cone scale must be positive")
   end subroutine test_bad_rmesh_cone_scale

   !> An infinite tangent-cone scale is refused
   subroutine test_bad_rmesh_cone_scale_infinite(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_error_type), allocatable :: failed

      call init_mesh(seeder, failed, cone_scale=ieee_value(1.0_wp, ieee_positive_inf))
      call expect_error(error, failed, "Surface mesh seeder: cone scale must be finite")
   end subroutine test_bad_rmesh_cone_scale_infinite

   !> Anchors without an atom center are refused
   subroutine test_bad_rmesh_anchors_without_atoms(error)
      type(error_type), allocatable, intent(out) :: error
      type(projection_seeder_type) :: seeder
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(seed_build_type) :: geometry
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Fixture failure, discarded: a clean return fails the inverted test
      type(error_type), allocatable :: setup

      call sphere_lsf(lsf, rs, setup)
      if (allocated(setup)) then
         ! Discarded with stat 0, so the test-drive finalizer does not escalate
         setup%stat = 0
         return
      end if
      call init_mesh(seeder, failed)
      if (allocated(failed)) return
      call single_atom_input(geometry)
      deallocate (geometry%centers, geometry%radii)
      allocate (geometry%centers(3, 0), geometry%radii(0))
      call seeder%build(lsf, geometry, failed)
      call expect_error(error, failed, "Surface mesh seeder: anchors without an atom center")
   end subroutine test_bad_rmesh_anchors_without_atoms

   !> The back-facing filter is skipped for a request without an owner
   !>
   !> - one sphere and a ball that holds it whole: with the owner the filter
   !>   drops the far side, without an owner every sample is kept, as for an
   !>   anchor without an accepted local point
   !>
   !> @param[out] error A request without an owner was filtered
   subroutine test_rmesh_back_facing_needs_owner(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      type(projection_seeder_type) :: seeder
      type(seeder_scratch_type) :: scratch
      type(seed_build_type) :: geometry
      type(seed_request_type) :: request
      type(moist_error_type), allocatable :: failed
      real(wp) :: rs
      !> Seeds without the filter, with the owner, and without an owner
      integer :: n_all, n_owner, n_unknown

      call sphere_lsf(lsf, rs, error)
      if (allocated(error)) return
      call init_mesh(seeder, failed, filter_back_facing=.true.)
      if (allocated(failed)) then
         call test_failed(error, "init: "//failed%message)
         return
      end if
      call single_atom_input(geometry)
      geometry%anchors = reshape([rs + 0.4_wp, 0.0_wp, 0.0_wp], [3, 1])
      geometry%reach = 2.0_wp*rs + 1.0_wp
      call seeder%build(lsf, geometry, failed)
      if (allocated(failed)) then
         call test_failed(error, "build: "//failed%message)
         return
      end if

      ! The atom center is the origin, the default `owner_center`
      request%anchor = geometry%anchors(:, 1)
      request%pdist_start = geometry%reach
      request%owner = 1
      call count_seeds(n_all)
      if (allocated(error)) return
      request%has_accepted = .true.
      call count_seeds(n_owner)
      if (allocated(error)) return
      request%owner = 0
      call count_seeds(n_unknown)
      if (allocated(error)) return
      call check(error, n_owner < n_all, "the fixture has no back-facing sample to drop")
      if (allocated(error)) return
      call check(error, n_unknown, n_all, "a request without an owner was filtered")

   contains

      !> Seeds of the current request
      !>
      !> @param[out] n Seeds
      subroutine count_seeds(n)
         integer, intent(out) :: n

         n = 0
         call seeder%seed(scratch, lsf, request, failed)
         if (allocated(failed)) then
            call test_failed(error, "seed: "//failed%message)
            return
         end if
         n = scratch%n_points
      end subroutine count_seeds

   end subroutine test_rmesh_back_facing_needs_owner

   !* ================================================================================= *!
   !*                              Spherical chart seeder                               *!
   !* ================================================================================= *!


end module test_cavity_drop_projection
