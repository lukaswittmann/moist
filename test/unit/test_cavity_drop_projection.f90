!> Regression tests for distinct minima in connected DROP surface patches
module test_cavity_drop_projection
   use mctc_env, only: wp, moist_error_type => error_type
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_radii, only: default_cpcm_radii
   use moist_math_linalg, only: setup_tangent_frame
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_projector, only: drop_projector_type
   use moist_cavity_drop_types, only: projection_workspace_type
   use moist_math_cell_grid, only: moist_cell_grid_type
   use moist_math_solver_octree_branch, only: moist_math_octree_branch_type
   implicit none(type, external)
   private

   public :: collect_cavity_drop_projection

   !> Borrowed projector for the unpruned surface-coverage reference
   type :: brute_force_context_type
      !> Geometry evaluator, separate from the finder under test
      type(drop_projector_type), pointer :: projector => null()
   end type brute_force_context_type

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
         & new_unittest("stationary_box_bounds", test_stationary_box_bounds) &
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
      !> Production projector
      type(drop_projector_type), target :: projector
      !> Projection scratch
      type(projection_workspace_type) :: work
      !> Molecular geometry
      type(structure_type) :: mol
      !> Complete atom lookup
      type(moist_cell_grid_type) :: atom_grid
      !> Native error
      type(moist_error_type), allocatable :: native_error
      !> Atom at the origin and exact point on a unit sphere
      real(wp) :: positions(3, 1), anchor(3)
      !> Returned point count
      integer :: n_points

      positions = 0.0_wp
      anchor = [1.0_wp, 0.0_wp, 0.0_wp]
      call new(mol, [1], positions)
      param%branch_search_rho2_slack = 0.0_wp
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      call projector%init(param, lsf, param%branch_sep_cut, verbosity=0)
      call atom_grid%build(positions, [1.0_wp], full_scan_below=200)
      call projector%init_primitives(mol, [1.0_wp], atom_grid)
      call projector%project_point(anchor, n_points=n_points, owner=1, index=1, &
                                   proj_level=9, error=native_error, work=work)
      if (allocated(projector%abort_error)) then
         call test_failed(error, projector%abort_error%message)
         return
      end if
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, 1)
      if (allocated(error)) return
      call check(error, norm2(work%points(:, 1) - anchor), 0.0_wp, thr=1.0e-14_wp)
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

   !> Compare the production projector with independently enumerated minima
   !>
   !> Reference construction: solve d/dtheta ||y(theta) - anchor||^2 = 0 on
   !> the exact two-focus meridian, retaining positive tangent curvature
   !> No projected neighbour points enter the production solver
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
      !> Production projector
      type(drop_projector_type) :: projector
      !> Candidate storage
      type(projection_workspace_type) :: work
      !> Molecular geometry
      type(structure_type) :: mol
      !> Atom lookup required by the production projector
      type(moist_cell_grid_type) :: atom_grid
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
      real(wp) :: distance
      !> Unpruned surface-leaf reference points
      real(wp), allocatable :: reference(:, :)
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
      if (present(sigma)) param%branch_weight_s = sigma
      if (present(window)) param%branch_search_rho2_slack = window
      call param%compute_derived(native_error)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type())
      lsf%screening_threshold = param%screening_threshold
      call projector%init(param, lsf, param%branch_sep_cut, &
                          tol=param%proj_tol, maxiter=param%proj_maxiter, verbosity=0)
      call atom_grid%build(positions, [radius, radius], full_scan_below=200)
      call projector%init_primitives(mol, [radius, radius], atom_grid)
      call projector%project_point(point, n_points=n_points, owner=1, index=1, &
                                   proj_level=9, error=native_error, work=work)
      if (allocated(projector%abort_error)) then
         call test_failed(error, projector%abort_error%message)
         return
      end if
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, n_points, size(expected, 2), message="Unexpected or duplicate branches")
      if (allocated(error)) return
      call check(error, all(work%converged(1:n_points)), message="Unconverged branch returned")
      if (allocated(error)) return
      do i = 1, size(expected, 2)
         distance = huge(1.0_wp)
         do j = 1, n_points
            distance = min(distance, norm2(work%points(:, j) - targets(:, i)))
         end do
         call check(error, distance < 1.0e-8_wp, message="Octree missed a strict minimum")
         if (allocated(error)) return
      end do
      call brute_force_branches(projector, point, 1, param%phi_alpha, &
                                projector%octree%dbg_rho_start, reference, native_error, reference_seeds)
      if (allocated(native_error)) then
         call test_failed(error, native_error%message)
         return
      end if
      call check(error, size(reference, 2), size(expected, 2), message="Brute-force H2 count differs")
      if (allocated(error)) return
      call check_branch_sets(reference, work%points(:, 1:n_points), error)
      if (allocated(error)) return
      if (present(require_pruning)) then
         if (require_pruning) then
            call check(error, 2*projector%octree%n_seeds < reference_seeds, &
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
   !> @param[in,out] projector Geometry and local refinement implementation
   !> @param[in] anchor Reference point
   !> @param[in] owner Owner atom
   !> @param[in] alpha Quadratic projection coefficient
   !> @param[in] radius Initial search radius enclosing the admissible window
   !> @param[out] points Distinct strict minima
   !> @param[out] error Native evaluation or search failure
   !> @param[out] n_seeds Optional exhaustive seed count
   subroutine brute_force_branches(projector, anchor, owner, alpha, radius, points, error, n_seeds)
      !> Borrowed projector
      type(drop_projector_type), intent(inout), target :: projector
      !> Reference point
      real(wp), intent(in) :: anchor(3)
      !> Owner atom
      integer, intent(in) :: owner
      !> Objective coefficient
      real(wp), intent(in) :: alpha
      !> Surface-search radius
      real(wp), intent(in) :: radius
      !> Reference minima
      real(wp), allocatable, intent(out) :: points(:, :)
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Reference work count
      integer, intent(out), optional :: n_seeds
      !> Original surface-only octree
      type(moist_math_octree_branch_type) :: tree
      !> Borrowed geometry callback
      type(brute_force_context_type) :: context
      !> Failed individual refinements
      type(moist_error_type), allocatable :: seed_error
      !> Candidate points and multipliers
      real(wp), allocatable :: candidates(:, :), multipliers(:)
      !> Distance-window and deduplication mask
      logical, allocatable :: keep(:)
      !> Surface data and tangent Hessian
      real(wp) :: value, gradient(3), hessian(3, 3), t1(3), t2(3), b11, b12, b22, beta
      !> Seed, solution, and multiplier
      real(wp) :: seed(3), point(3), lambda, lsf_anchor
      !> Seed and accepted point counts
      integer :: i, n, n_strict

      if (present(n_seeds)) n_seeds = 0
      context%projector => projector
      call tree%init(seed_size=projector%octree_seed_size, &
                     max_boxes=projector%octree_max_boxes, &
                     max_survivors=projector%octree_max_survivors, &
                     max_depth=projector%octree_max_depth, error=error)
      if (allocated(error)) return
      call projector%lsf%set_max_deriv(0)
      projector%ssd_cache_valid = .false.
      call projector%compute_ssd(anchor)
      call projector%lsf%f0(lsf_anchor)
      call tree%run(anchor, lsf_anchor, radius, projector%branch_rho2_slack, &
                    brute_force_probe, context, error)
      call projector%lsf%set_max_deriv(2)
      projector%ssd_cache_valid = .false.
      if (allocated(error)) return
      if (allocated(projector%abort_error)) then
         error = projector%abort_error
         return
      end if
      if (present(n_seeds)) n_seeds = tree%n_seeds
      allocate (candidates(3, tree%n_seeds + 1), multipliers(tree%n_seeds + 1))
      n = 0
      do i = 0, tree%n_seeds
         seed = anchor
         if (i > 0) seed = tree%seeds(:, i)
         call projector%refine_point(anchor, seed, &
                                     anchor - projector%max_displacement_threshold, &
                                     anchor + projector%max_displacement_threshold, &
                                     projector%max_displacement_threshold, owner, 0, &
                                     seed_error, point, lambda_out=lambda)
         if (allocated(projector%abort_error)) then
            error = projector%abort_error
            return
         end if
         if (allocated(seed_error)) then
            deallocate (seed_error)
            cycle
         end if
         n = n + 1
         candidates(:, n) = point
         multipliers(n) = lambda
      end do
      if (n == 0) then
         allocate (points(3, 0))
         return
      end if
      call projector%filter_candidates(candidates(:, 1:n), n, keep, error, anchor)
      if (allocated(error)) return
      n_strict = 0
      do i = 1, n
         if (.not. keep(i)) cycle
         call projector%compute_ssd(candidates(:, i))
         call projector%lsf%f012_r(value, gradient, hessian)
         if (abs(value) > 1.0e-9_wp .or. norm2(gradient) < 1.0e-12_wp) cycle
         if (norm2(alpha*(candidates(:, i) - anchor) - multipliers(i)*gradient) > 1.0e-8_wp) cycle
         call setup_tangent_frame(gradient/norm2(gradient), t1, t2)
         b11 = alpha - multipliers(i)*dot_product(t1, matmul(hessian, t1))
         b12 = -multipliers(i)*dot_product(t1, matmul(hessian, t2))
         b22 = alpha - multipliers(i)*dot_product(t2, matmul(hessian, t2))
         beta = 0.5_wp*(b11 + b22 - sqrt((b11 - b22)**2 + 4.0_wp*b12**2))
         if (beta <= 1.0e-8_wp) cycle
         n_strict = n_strict + 1
         candidates(:, n_strict) = candidates(:, i)
      end do
      points = candidates(:, 1:n_strict)
   end subroutine brute_force_branches

   !> Surface-only callback for the brute-force reference
   !>
   !> @param[in] x Evaluation point
   !> @param[out] lsf0 Surface value
   !> @param[out] radius Surface-free radius
   !> @param[in] context Borrowed reference geometry
   subroutine brute_force_probe(x, lsf0, radius, context)
      !> Evaluation point
      real(wp), intent(in) :: x(3)
      !> Surface value
      real(wp), intent(out) :: lsf0
      !> Surface exclusion radius
      real(wp), intent(out) :: radius
      !> Borrowed callback state
      class(*), intent(in) :: context

      lsf0 = 0.0_wp
      radius = 0.0_wp
      select type (context)
      type is (brute_force_context_type)
         call context%projector%compute_ssd(x)
         call context%projector%lsf%f0(lsf0)
         radius = context%projector%lsf%exclusion_radius(lsf0)
      end select
   end subroutine brute_force_probe

   !> Check the dimer fixture from the deflation comparison
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

   !> Check the cross fixture from the deflation comparison
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

   !> Check the octahedron fixture from the deflation comparison
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

   !> Check the cube fixture from the deflation comparison
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

   !> Check the pentagon fixture from the deflation comparison
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

   !> Check the tetrahedron fixture from the deflation comparison
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

   !> Require octree recovery of admissible minima found by fine multistart
   !>
   !> Level 8 is a comparison set, not a completeness certificate
   !> Extra valid level-9 branches are allowed
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
      !> Comparison and octree projections before outgoing weighting
      type(cavity_type_drop), allocatable :: comparison, octree
      !> Native construction error
      type(moist_error_type), allocatable :: native_error
      !> Candidate bounds and comparison counts per anchor
      integer, allocatable :: first(:), last(:), counts(:)
      !> Point indices, matching count, and original anchor id
      integer :: i, j, aid, matches
      !> Branch comparison message
      character(len=256) :: message

      call new_context(ctx, verbosity=0, nthreads=1)
      allocate (comparison, octree)
      call build_cluster_projection(comparison, ctx, mol, nleb, blend_k, 8, native_error)
      if (allocated(native_error)) then
         call test_failed(error, title//": level 8: "//native_error%message)
         return
      end if
      call build_cluster_projection(octree, ctx, mol, nleb, blend_k, 9, native_error)
      if (allocated(native_error)) then
         call test_failed(error, title//": level 9: "//native_error%message)
         return
      end if
      call check(error, octree%nmax, comparison%nmax, message="Reference anchors differ")
      if (allocated(error)) return
      call check_projection_geometry(comparison, error)
      if (allocated(error)) return
      call check_projection_geometry(octree, error)
      if (allocated(error)) return

      allocate (first(mol%nat*nleb), source=0)
      allocate (last(mol%nat*nleb), source=0)
      allocate (counts(mol%nat*nleb), source=0)
      do i = 1, octree%ngrid
         aid = octree%anchor_id(i)
         if (first(aid) == 0) first(aid) = i
         last(aid) = i
      end do

      do i = 1, comparison%ngrid
         if (.not. relevant_projection(comparison, i)) cycle
         aid = comparison%anchor_id(i)
         counts(aid) = counts(aid) + 1
         matches = 0
         if (first(aid) > 0) then
            do j = first(aid), last(aid)
               if (.not. relevant_projection(octree, j)) cycle
               if (norm2(comparison%xyz(:, i) - octree%xyz(:, j)) <= 1.0e-6_wp) matches = matches + 1
            end do
         end if
         write (message, "(a,a,i0,a,i0)") title, ": anchor ", aid, ", matching octree branches ", matches
         call check(error, matches, 1, message=trim(message))
         if (allocated(error)) return
      end do
      call check(error, sum(counts) > 0, message="Comparison set contains no admissible minimum")
      if (allocated(error)) return
      call check_cluster_brute_force(octree, error)
      if (allocated(error)) return
      deallocate (comparison, octree)
      call ctx%delete()
   end subroutine check_cluster_branches

   !> Compare each cluster anchor against exhaustive surface-leaf refinement
   !>
   !> @param[in] cavity Raw finder result
   !> @param[out] error Missing reference branch
   subroutine check_cluster_brute_force(cavity, error)
      !> Raw projected points
      type(cavity_type_drop), intent(in) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Independent local-solver state
      type(drop_projector_type), target :: projector
      !> Native reference failure
      type(moist_error_type), allocatable :: native_error
      !> Reference minima
      real(wp), allocatable :: reference(:, :)
      !> Anchor and a conservative radius from a known surface point
      real(wp) :: anchor(3), radius
      !> Anchor-group bounds and total reference count
      integer :: first, last, j, total

      call projector%init(cavity%param, cavity%lsf_model, cavity%param%branch_sep_cut, &
                          tol=cavity%param%proj_tol, maxiter=cavity%param%proj_maxiter, verbosity=0)
      call projector%init_primitives(cavity%mol, cavity%radii, cavity%mol_cell_grid)
      first = 1
      total = 0
      do while (first <= cavity%ngrid)
         last = first
         do while (last < cavity%ngrid)
            if (cavity%anchor_id(last + 1) /= cavity%anchor_id(first)) exit
            last = last + 1
         end do
         anchor = cavity%anchorxyz(:, first)
         radius = huge(1.0_wp)
         do j = first, last
            if (.not. cavity%converged(j)) cycle
            radius = min(radius, sum((cavity%xyz(:, j) - anchor)**2))
         end do
         call check(error, radius < huge(1.0_wp), message="No surface point bounds the reference search")
         if (allocated(error)) return
         radius = sqrt(radius + cavity%param%branch_rho2_slack())
         call brute_force_branches(projector, anchor, cavity%owner(first), cavity%param%phi_alpha, &
                                   radius, reference, native_error)
         if (allocated(native_error)) then
            call test_failed(error, native_error%message)
            return
         end if
         total = total + size(reference, 2)
         call check_branch_sets(reference, cavity%xyz(:, first:last), error)
         if (allocated(error)) return
         first = last + 1
      end do
      call check(error, total > 0, message="Brute-force cluster reference is empty")
   end subroutine check_cluster_brute_force

   !> Build raw projected representations with the development fixtures' settings
   !>
   !> The outgoing softmax and final point filter are deliberately not applied
   !>
   !> @param[out] cavity Raw projected cavity and geometric admissibility factors
   !> @param[in,out] ctx Borrowed context
   !> @param[in] mol Cluster geometry in bohr
   !> @param[in] nleb Reference quadrature size
   !> @param[in] blend_k Level-set softness
   !> @param[in] level Projection solver level
   !> @param[out] error Native construction failure
   subroutine build_cluster_projection(cavity, ctx, mol, nleb, blend_k, level, error)
      !> Projection result
      type(cavity_type_drop), intent(out) :: cavity
      !> Borrowed construction context
      type(moist_context_type), intent(inout), target :: ctx
      !> Cluster geometry
      type(structure_type), intent(in) :: mol
      !> Points per reference sphere
      integer, intent(in) :: nleb
      !> Level-set softness
      real(wp), intent(in) :: blend_k
      !> Solver level
      integer, intent(in) :: level
      !> Native failure
      type(moist_error_type), allocatable, intent(out) :: error
      !> Numerical settings shared by both searches
      type(moist_cavity_drop_parameters_type) :: param
      !> Development fixture level set
      type(moist_cavity_drop_lsf_svdw_type) :: lsf

      param%num_leb = nleb
      param%proj_level = level
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
      call cavity%fill_arrays(error)
      if (allocated(error)) return
      call cavity%setup_mol_cell_grid(error)
      if (allocated(error)) return
      call cavity%compute_switching(error)
      if (allocated(error)) return
      call cavity%filter_arrays("Reference prefilter", error)
      if (allocated(error)) return
      cavity%nmax = cavity%ngrid
      call cavity%project_all_points(error)
      if (allocated(error)) return
      call cavity%compute_cp_jacobian_scaling(error)
   end subroutine build_cluster_projection

   !> Whether a raw representation has nonnegligible geometric contribution
   !>
   !> @param[in] cavity Raw projection and geometric factors
   !> @param[in] i Point index
   pure function relevant_projection(cavity, i) result(relevant)
      !> Projection result
      type(cavity_type_drop), intent(in) :: cavity
      !> Point index
      integer, intent(in) :: i
      !> Positive geometric contribution before outgoing branch fractions
      logical :: relevant

      relevant = cavity%converged(i) .and. &
                 abs(cavity%wleb(i)*cavity%f(i)*cavity%cpjac_scal0(i)*cavity%w_f0(i)) > cavity%param%wleb_cut
   end function relevant_projection

   !> Verify surface membership, stationarity, and strict local minimality
   !>
   !> Checks use the defining equations, independently of branch matching
   !>
   !> @param[in,out] cavity Projected points and level-set evaluator
   !> @param[out] error Invalid or duplicate admissible representation
   subroutine check_projection_geometry(cavity, error)
      !> Projection result and mutable evaluation cache
      type(cavity_type_drop), intent(inout) :: cavity
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Evaluation failure
      type(moist_error_type), allocatable :: native_error
      !> Surface value, gradient, and Hessian
      real(wp) :: value, gradient(3), hessian(3, 3)
      !> Tangent frame and constrained Hessian entries
      real(wp) :: tangent1(3), tangent2(3), b11, b12, b22, beta_min
      !> Stationarity residual and objective coefficient
      real(wp) :: residual(3), alpha
      !> Point and sibling indices
      integer :: i, j

      alpha = cavity%param%phi_alpha
      call cavity%lsf_model%set_max_deriv(2)
      do i = 1, cavity%ngrid
         if (.not. relevant_projection(cavity, i)) cycle
         call cavity%lsf_model%prepare(cavity%xyz(:, i), native_error)
         if (allocated(native_error)) then
            call test_failed(error, native_error%message)
            return
         end if
         call cavity%lsf_model%f012_r(value, gradient, hessian)
         residual = alpha*(cavity%xyz(:, i) - cavity%anchorxyz(:, i)) - cavity%lambda0(i)*gradient
         call check(error, abs(value) < 1.0e-9_wp, message="Projected point is off the surface")
         if (allocated(error)) return
         call check(error, norm2(residual) < 1.0e-8_wp, message="Projected point is not stationary")
         if (allocated(error)) return
         call setup_tangent_frame(gradient/norm2(gradient), tangent1, tangent2)
         b11 = alpha - cavity%lambda0(i)*dot_product(tangent1, matmul(hessian, tangent1))
         b12 = -cavity%lambda0(i)*dot_product(tangent1, matmul(hessian, tangent2))
         b22 = alpha - cavity%lambda0(i)*dot_product(tangent2, matmul(hessian, tangent2))
         beta_min = 0.5_wp*(b11 + b22 - sqrt((b11 - b22)**2 + 4.0_wp*b12**2))
         call check(error, beta_min > 1.0e-8_wp, message="Admissible point is not a strict local minimum")
         if (allocated(error)) return
         do j = i + 1, cavity%ngrid
            if (cavity%anchor_id(j) /= cavity%anchor_id(i)) exit
            if (.not. relevant_projection(cavity, j)) cycle
            call check(error, norm2(cavity%xyz(:, i) - cavity%xyz(:, j)) > 1.0e-7_wp, &
                       message="Duplicate admissible branch")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_projection_geometry

end module test_cavity_drop_projection
