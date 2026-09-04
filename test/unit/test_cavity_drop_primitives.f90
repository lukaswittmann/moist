module test_cavity_drop_primitives
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, fd4_scalar, &
                           fd4_offsets, rel_deviation
   use moist_utils_env, only: get_env
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_objective_phi, only: moist_cavity_drop_objective_phi_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_smooth_step_type, &
                                          new_swif_smooth_step, &
                                          moist_cavity_drop_swif_sigmoid_bump_type, &
                                          new_swif_sigmoid_bump
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig, new_iswig, &
                                         iswig_workspace_type
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
                                                   drop_seed_result_type, &
                                                   drop_seed_state_tangent_type, &
                                                   drop_seed_input_tangent_type, &
                                                   drop_seed_result_tangent_type, &
                                                   build_seed_state, apply_seed, &
                                                   apply_seed_tangent, &
                                                   switched_eigenvalue_response, &
                                                   switched_eigenvalue_curvature, &
                                                   seed_state_ok, seed_curv_disc_guard
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_getrf, only: lapack_getrf
   use moist_math_lapack_getrs, only: lapack_getrs
   use moist_math_lapack_kinds, only: lapack_ik
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_primitives

   integer, parameter :: ndim = 3
   integer, parameter :: owner_dummy = 1

   !> Finite-difference thresholds of the phi-primitive comparisons.
   !>
   !> Tightened to the project-wide `1e-10 / 1e-10` target on 2026-09-03 from
   !> `5e-9 / 5e-9`. Measured with the failure suppressed so the whole block is
   !> scanned rather than stopping at the first exceedance: the worst deviation
   !> anywhere in the phi ladder is `1.75e-10` absolute against a reference of
   !> `10.19` (relative `1.7e-11`, which the `max(thr_abs, thr_rel * expected)`
   !> rule passes), and `phi_f1_r` / `phi_f012_r` fail on a `1.53e-10` deviation
   !> against a reference of `-0.777`, where the absolute bound is the binding
   !> one. This block is therefore *at* `1e-10`, not orders from it -- one more
   !> decade of step refinement would likely close it.
   real(wp), parameter :: ABS_THR = 1.0e-10_wp
   real(wp), parameter :: REL_THR = 1.0e-10_wp
   real(wp), parameter :: STEP_SIZE = 1.0e-4_wp

   !* --------------------------- iSwiG switching fixtures -------------------------- *!

   !> Structures in the iSwiG probe set; must be a positive multiple of five
   integer, parameter :: n_iswig_structures = 5
   !> Owner atoms probed per structure
   integer, parameter :: n_iswig_owners = 3
   !> Probe offsets per owner
   integer, parameter :: n_iswig_dirs = 4
   !> Probes per structure and neighbour-list path
   integer, parameter :: n_iswig_probes = n_iswig_owners*n_iswig_dirs

   !> Born width of the nleb = 110 Lebedev grid (`parameters.f90:24`)
   real(wp), parameter :: iswig_swx = 4.900490_wp
   !> Representative largest raw Lebedev weight, as `setup.f90:35` passes it
   real(wp), parameter :: iswig_wleb_max = 1.25e-2_wp

   !> Tolerance of the golden comparison, measured as `rel_deviation`, `|a - b| / (1 + |b|)`
   real(wp), parameter :: iswig_golden_tol = 1.0e-12_wp

   !> Gaussian width used by the finite-difference tests
   real(wp), parameter :: iswig_fd_xi = 1.6_wp

   !> Step sizes of the iSwiG finite differences, coarse first
   real(wp), parameter :: iswig_fd_steps(2) = [4.0e-3_wp, 2.0e-3_wp]

   !> Tolerance of the block-versus-sparse comparison
   real(wp), parameter :: iswig_block_tol = 1.0e-13_wp

   !> Tolerances of the iSwiG finite-difference comparisons.
   !>
   !> Tightened to the project-wide `1e-10 / 1e-10` target on 2026-09-03 from
   !> `1e-8 / 1e-8`. Measured worst deviation over the whole block, failure
   !> suppressed: `3.48e-9` absolute at a reference of `2.214` in
   !> `iswig_swi1_ra_fd` and `3.85e-9` at `-1.306` in `iswig_swi2_rarb_fd`, so
   !> both fail by a factor of roughly 30. `iswig_fd_steps` is `4e-3 / 2e-3`,
   !> coarse enough that this is truncation of the difference rather than the
   !> arithmetic of the analytic row: the step, not the tolerance, is what would
   !> have to move.
   real(wp), parameter :: iswig_fd_abs = 1.0e-10_wp
   real(wp), parameter :: iswig_fd_rel = 1.0e-10_wp

   !> Angular offsets from the owner/neighbour sphere-intersection circle
   real(wp), parameter :: iswig_theta_off(n_iswig_dirs) = &
                          [-0.12_wp, -0.04_wp, 0.04_wp, 0.12_wp]

   !* --------------------------- bordered KKT fixtures ---------------------------- *!

   !> Bordered-KKT systems exercised by the factor tests
   integer, parameter :: kkt_ncase = 6
   !> Right-hand sides per system; the widest batch in the gradient path
   integer, parameter :: kkt_nrhs = 7

   !> Tolerance of the `gesv` comparison, measured as `rel_deviation`
   real(wp), parameter :: kkt_gesv_tol = 1.0e-12_wp

   !* ------------------------ seed-state tangent fixtures -------------------------- *!

   !> Central-difference steps of the seed-state tangent test
   real(wp), parameter :: seed_fd_steps(2) = [1.0e-5_wp, 1.0e-6_wp]

   !> Tolerance of the seed-state tangent comparison, measured as `rel_deviation`.
   !>
   !> Tightened to the project-wide `1e-10` target on 2026-09-03 from `1e-6`.
   !> Measured worst `rel_deviation` over the whole seed block, failure
   !> suppressed: `1.73e-8` on `dinv_J` at `h = 1e-5`, two orders above the
   !> bound. Nine of the ten seed tangent cases fail; the smallest failing
   !> deviation is `1.18e-10` on `dJ` at `h = 1e-6`, so the block spans the
   !> bound rather than sitting uniformly above it.
   real(wp), parameter :: seed_fd_tol = 1.0e-10_wp

   !> Tolerance for the seed-linearity identity, which is algebra rather than a
   !> finite difference: `apply_seed_tangent` is linear in the seed, so summing
   !> two single-entry seeds must reproduce the pair seed to rounding
   real(wp), parameter :: seed_linearity_tol = 1.0e-12_wp

   !> Smallest `|B12|` the fixture may reach
   real(wp), parameter :: seed_fd_b12_min = 1.0e-6_wp

   !> Smallest eigenvalue gap `sqrt_disc_B` the fixture may reach
   real(wp), parameter :: seed_fd_gap_min = 0.1_wp

   !> Smallest `proj_surf` the fixture may reach
   real(wp), parameter :: seed_fd_proj_min = 0.5_wp

   !> Smallest `disc_curv` the curvature fixtures may reach
   real(wp), parameter :: seed_fd_disc_min = 1.0e-3_wp

   !> Anti-vacuity floor on the finite-difference reference
   real(wp), parameter :: seed_fd_live_min = 1.0e-4_wp

   !* ------------------- switched-eigenvalue identity fixtures --------------------- *!

   !> Independent random draws of the switched-eigenvalue identity tests
   integer, parameter :: seed_identity_ndraw = 7

   !> Multiplier and modulus of the identity fixture's Lehmer generator; the
   !> product stays inside `int64`, so no draw can overflow
   integer(int64), parameter :: seed_identity_mult = 16807_int64
   integer(int64), parameter :: seed_identity_mod = 2147483647_int64

   !> First state of that generator, so every run draws the same numbers
   integer(int64), parameter :: seed_identity_seed = 20260903_int64

   !> Tolerance of the identity comparison, measured as `rel_deviation`. These
   !> are algebraic identities rather than finite differences, so the only
   !> admissible gap is accumulated roundoff
   real(wp), parameter :: seed_identity_tol = 1.0e-12_wp

   !> Anti-vacuity floor on the naive-assembly reference of a draw
   real(wp), parameter :: seed_identity_live_min = 1.0e-3_wp

contains

   !> Collect the phi primitive FD tests plus the switching nuclear gradient
   subroutine collect_cavity_drop_primitives(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("phi_f0", test_phi_f0), &
                  new_unittest("phi_f1_r", test_phi_f1_r), &
                  new_unittest("phi_f2_rr", test_phi_f2_rr), &
                  new_unittest("phi_f3_rrr", test_phi_f3_rrr), &
                  new_unittest("phi_f4_rrrr", test_phi_f4_rrrr), &
                  new_unittest("phi_f1_ra", test_phi_f1_rA), &
                  new_unittest("phi_f2_rarb", test_phi_f2_rArB), &
                  new_unittest("phi_f2_r_ra", test_phi_f2_r_rA), &
                  new_unittest("phi_f012_r", test_phi_f012_r), &
                  new_unittest("switching_f1_ra", test_switching_f1_rA), &
                  new_unittest("iswig_swi1_matches_golden", test_iswig_swi1_matches_golden), &
                  new_unittest("iswig_swi1_ra_fd", test_iswig_swi1_rA_fd), &
                  new_unittest("iswig_dxi_fd", test_iswig_dxi_fd), &
                  new_unittest("iswig_swi2_rarb_fd", test_iswig_swi2_rArB_fd), &
                  new_unittest("iswig_swi2_xi_channel_fd", test_iswig_swi2_xi_channel_fd), &
                  new_unittest("iswig_swi2_symmetry", test_iswig_swi2_symmetry), &
                  new_unittest("iswig_translation_null", test_iswig_translation_null), &
                  new_unittest("iswig_isolated_and_fallback", test_iswig_isolated_and_fallback), &
                  new_unittest("iswig_swi2_block_matches_sparse", &
                               test_iswig_swi2_block_matches_sparse), &
                  new_unittest("iswig_swi2_block_symmetry", test_iswig_swi2_block_symmetry), &
                  new_unittest("iswig_swi2_block_translation", &
                               test_iswig_swi2_block_translation), &
                  new_unittest("iswig_swi2_block_xi_fd", test_iswig_swi2_block_xi_fd), &
                  new_unittest("iswig_swi2_block_guarded", test_iswig_swi2_block_guarded), &
                  new_unittest("kkt_factor_matches_lu", test_kkt_factor_matches_lu), &
                  new_unittest("kkt_factor_matches_gesv", test_kkt_factor_matches_gesv), &
                  new_unittest("kkt_factor_reuse", test_kkt_factor_reuse), &
                  new_unittest("kkt_factor_singular", test_kkt_factor_singular), &
                  new_unittest("switched_eigenvalue_response_identity", &
                               test_switched_eigenvalue_response_identity), &
                  new_unittest("switched_eigenvalue_curvature_identity", &
                               test_switched_eigenvalue_curvature_identity), &
                  new_unittest("seed_state_tangent_plain", test_seed_state_tangent_plain), &
                  new_unittest("seed_state_tangent_prune", test_seed_state_tangent_prune), &
                  new_unittest("seed_state_tangent_curvature", &
                               test_seed_state_tangent_curvature), &
                  new_unittest("seed_state_tangent_curvature_prune", &
                               test_seed_state_tangent_curvature_prune), &
                  new_unittest("seed_result_tangent_plain", test_seed_result_tangent_plain), &
                  new_unittest("seed_result_tangent_prune", test_seed_result_tangent_prune), &
                  new_unittest("seed_result_tangent_curvature", &
                               test_seed_result_tangent_curvature), &
                  new_unittest("seed_result_asym_pair_plain", &
                               test_seed_result_asym_pair_plain), &
                  new_unittest("seed_result_asym_pair_curvature", &
                               test_seed_result_asym_pair_curvature), &
                  new_unittest("seed_result_tangent_curvature_prune", &
                               test_seed_result_tangent_curvature_prune) &
                  ]
   end subroutine collect_cavity_drop_primitives

   !> Test phi value against the direct quadratic expression.
   subroutine test_phi_f0(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim), analytic, reference
      integer :: icase, ipt, jpt

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jpt = 1, size(points, 2)
               anchor = points(:, jpt)
               analytic = phi%f0(point, anchor, owner_dummy)
               reference = 0.5_wp*param%phi_alpha*sum((point - anchor)**2)
               call check(error, analytic, reference, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_phi_f0

   !> Test phi point gradient against a finite difference of the value.
   subroutine test_phi_f1_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim), analytic(ndim), numeric(ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: work_point(ndim)
      integer :: icase, ipt, jpt, i, axis
      real(wp) :: h

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jpt = 1, size(points, 2)
               anchor = points(:, jpt)
               analytic = phi%f1_r(point, anchor, owner_dummy)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  f_pp = phi%f0(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  f_p = phi%f0(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  f_m = phi%f0(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  f_mm = phi%f0(work_point, anchor, owner_dummy)
                  numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
               end do

               do i = 1, ndim
                  call check(error, analytic(i), numeric(i), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_phi_f1_r

   !> Test phi point Hessian against a finite difference of the point gradient.
   subroutine test_phi_f2_rr(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim), analytic(ndim, ndim), numeric(ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: work_point(ndim)
      integer :: icase, ipt, jpt, i, j, axis
      real(wp) :: h

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jpt = 1, size(points, 2)
               anchor = points(:, jpt)
               analytic = phi%f2_rr(point, anchor, owner_dummy)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  g_pp = phi%f1_r(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  g_p = phi%f1_r(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  g_m = phi%f1_r(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  g_mm = phi%f1_r(work_point, anchor, owner_dummy)
                  do i = 1, ndim
                     numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
                  end do
               end do

               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_phi_f2_rr

   !> Test phi third point derivative against a finite difference of the Hessian.
   !> Mol-less pure-math test: phi = 0.5*alpha*(r-anchor)^2 has analytic
   !> derivatives that hold for any two distinct points, so the (point,
   !> anchor) pair is hard-coded rather than sourced from a molecule.
   subroutine test_phi_f3_rrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      real(wp) :: analytic(ndim, ndim, ndim), numeric(ndim, ndim, ndim)
      real(wp) :: hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
      real(wp) :: point(ndim), anchor(ndim), work_point(ndim)
      integer :: i, j, k, axis
      real(wp) :: h
      real(wp), parameter :: points(ndim, 2) = reshape([ &
                                                       -1.24_wp, 0.56_wp, 0.20_wp, &
                                                       0.60_wp, -0.90_wp, 0.80_wp], [ndim, 2])

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE
      point = points(:, 1)
      anchor = points(:, 2)
      analytic = phi%f3_rrr(point, anchor, owner_dummy)

      do axis = 1, ndim
         work_point = point
         work_point(axis) = point(axis) + 2.0_wp*h
         hess_pp = phi%f2_rr(work_point, anchor, owner_dummy)
         work_point = point
         work_point(axis) = point(axis) + h
         hess_p = phi%f2_rr(work_point, anchor, owner_dummy)
         work_point = point
         work_point(axis) = point(axis) - h
         hess_m = phi%f2_rr(work_point, anchor, owner_dummy)
         work_point = point
         work_point(axis) = point(axis) - 2.0_wp*h
         hess_mm = phi%f2_rr(work_point, anchor, owner_dummy)
         do i = 1, ndim
            do j = 1, ndim
               numeric(i, j, axis) = fd4_scalar( &
                                     hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), h)
            end do
         end do
      end do

      do i = 1, ndim
         do j = 1, ndim
            do k = 1, ndim
               call check(error, analytic(i, j, k), numeric(i, j, k), &
                          thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_phi_f3_rrr

   !> Test phi fourth point derivative against a finite difference of the third derivative.
   !> Mol-less pure-math test (see test_phi_f3_rrr).
   subroutine test_phi_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      real(wp) :: analytic(ndim, ndim, ndim, ndim), numeric(ndim, ndim, ndim, ndim)
      real(wp) :: third_pp(ndim, ndim, ndim), third_p(ndim, ndim, ndim)
      real(wp) :: third_m(ndim, ndim, ndim), third_mm(ndim, ndim, ndim)
      real(wp) :: point(ndim), anchor(ndim), work_point(ndim)
      integer :: i, j, k, m, axis
      real(wp) :: h
      real(wp), parameter :: points(ndim, 2) = reshape([ &
                                                       -1.24_wp, 0.56_wp, 0.20_wp, &
                                                       0.60_wp, -0.90_wp, 0.80_wp], [ndim, 2])

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE
      point = points(:, 1)
      anchor = points(:, 2)
      analytic = phi%f4_rrrr(point, anchor, owner_dummy)

      do axis = 1, ndim
         work_point = point
         work_point(axis) = point(axis) + 2.0_wp*h
         third_pp = phi%f3_rrr(work_point, anchor, owner_dummy)
         work_point = point
         work_point(axis) = point(axis) + h
         third_p = phi%f3_rrr(work_point, anchor, owner_dummy)
         work_point = point
         work_point(axis) = point(axis) - h
         third_m = phi%f3_rrr(work_point, anchor, owner_dummy)
         work_point = point
         work_point(axis) = point(axis) - 2.0_wp*h
         third_mm = phi%f3_rrr(work_point, anchor, owner_dummy)
         do i = 1, ndim
            do j = 1, ndim
               do k = 1, ndim
                  numeric(i, j, k, axis) = fd4_scalar( &
                                           third_pp(i, j, k), third_p(i, j, k), third_m(i, j, k), third_mm(i, j, k), h)
               end do
            end do
         end do
      end do

      do i = 1, ndim
         do j = 1, ndim
            do k = 1, ndim
               do m = 1, ndim
                  call check(error, analytic(i, j, k, m), numeric(i, j, k, m), &
                             thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_phi_f4_rrrr

   !> Test phi nuclear gradient against a finite difference of the anchor point.
   subroutine test_phi_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim)
      real(wp), allocatable :: analytic(:, :), numeric(:, :)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: work_anchor(ndim)
      integer :: icase, owner, i, j, axis
      real(wp) :: h

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         point = points(:, 1)
         anchor = points(:, 2)
         do owner = 1, mol%nat
            analytic = phi%f1_rA(point, anchor, owner)
            allocate (numeric(ndim, mol%nat), source=0.0_wp)
            do axis = 1, ndim
               work_anchor = anchor
               work_anchor(axis) = anchor(axis) + 2.0_wp*h
               f_pp = phi%f0(point, work_anchor, owner)
               work_anchor = anchor
               work_anchor(axis) = anchor(axis) + h
               f_p = phi%f0(point, work_anchor, owner)
               work_anchor = anchor
               work_anchor(axis) = anchor(axis) - h
               f_m = phi%f0(point, work_anchor, owner)
               work_anchor = anchor
               work_anchor(axis) = anchor(axis) - 2.0_wp*h
               f_mm = phi%f0(point, work_anchor, owner)
               numeric(axis, owner) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
            end do

            do i = 1, ndim
               do j = 1, mol%nat
                  call check(error, analytic(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
            deallocate (numeric)
         end do
      end do
   end subroutine test_phi_f1_rA

   !> Test phi nuclear Hessian against finite differences of the nuclear gradient.
   subroutine test_phi_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: analytic(:, :, :, :), numeric(:, :, :, :)
      real(wp), allocatable :: g_pp(:, :), g_p(:, :), g_m(:, :), g_mm(:, :)
      real(wp) :: work_anchor(ndim)
      integer :: icase, owner, i, j, atom_a, atom_b, axis
      real(wp) :: h

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         do owner = 1, mol%nat
            analytic = phi%f2_rArB(points(:, 1), points(:, 2), owner)
            allocate (numeric(ndim, ndim, mol%nat, mol%nat), source=0.0_wp)
            do axis = 1, ndim
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) + 2.0_wp*h
               g_pp = phi%f1_rA(points(:, 1), work_anchor, owner)
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) + h
               g_p = phi%f1_rA(points(:, 1), work_anchor, owner)
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) - h
               g_m = phi%f1_rA(points(:, 1), work_anchor, owner)
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) - 2.0_wp*h
               g_mm = phi%f1_rA(points(:, 1), work_anchor, owner)
               do i = 1, ndim
                  do atom_a = 1, mol%nat
                     numeric(i, axis, atom_a, owner) = fd4_scalar( &
                                                       g_pp(i, atom_a), g_p(i, atom_a), g_m(i, atom_a), g_mm(i, atom_a), h)
                  end do
               end do
            end do

            do i = 1, ndim
               do j = 1, ndim
                  do atom_a = 1, mol%nat
                     do atom_b = 1, mol%nat
                        call check(error, analytic(i, j, atom_a, atom_b), numeric(i, j, atom_a, atom_b), &
                                   thr_abs=ABS_THR, thr_rel=REL_THR)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
            deallocate (numeric)
         end do
      end do
   end subroutine test_phi_f2_rArB

   !> Test phi mixed point-nuclear Hessian against finite differences of point gradient.
   subroutine test_phi_f2_r_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: analytic(:, :, :), numeric(:, :, :)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: work_anchor(ndim)
      integer :: icase, owner, i, j, atom, axis
      real(wp) :: h

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         do owner = 1, mol%nat
            analytic = phi%f2_r_rA(points(:, 1), points(:, 2), owner)
            allocate (numeric(ndim, ndim, mol%nat), source=0.0_wp)
            do axis = 1, ndim
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) + 2.0_wp*h
               g_pp = phi%f1_r(points(:, 1), work_anchor, owner)
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) + h
               g_p = phi%f1_r(points(:, 1), work_anchor, owner)
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) - h
               g_m = phi%f1_r(points(:, 1), work_anchor, owner)
               work_anchor = points(:, 2)
               work_anchor(axis) = points(axis, 2) - 2.0_wp*h
               g_mm = phi%f1_r(points(:, 1), work_anchor, owner)
               do i = 1, ndim
                  numeric(i, axis, owner) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
               end do
            end do

            do i = 1, ndim
               do j = 1, ndim
                  do atom = 1, mol%nat
                     call check(error, analytic(i, j, atom), numeric(i, j, atom), &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do
            deallocate (numeric)
         end do
      end do
   end subroutine test_phi_f2_r_rA

   !> Test combined phi value, gradient, and Hessian against finite differences.
   subroutine test_phi_f012_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_objective_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim), work_point(ndim)
      real(wp) :: analytic_val, numeric_val
      real(wp) :: analytic_grad(ndim), numeric_grad(ndim)
      real(wp) :: analytic_hess(ndim, ndim), numeric_hess(ndim, ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      integer :: icase, ipt, jpt, i, j, axis
      real(wp) :: h

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)
         call phi%set_input(mol, radii)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jpt = 1, size(points, 2)
               anchor = points(:, jpt)
               call phi%f012_r(point, anchor, owner_dummy, analytic_val, analytic_grad, analytic_hess)
               numeric_val = phi%f0(point, anchor, owner_dummy)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  f_pp = phi%f0(work_point, anchor, owner_dummy)
                  g_pp = phi%f1_r(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  f_p = phi%f0(work_point, anchor, owner_dummy)
                  g_p = phi%f1_r(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  f_m = phi%f0(work_point, anchor, owner_dummy)
                  g_m = phi%f1_r(work_point, anchor, owner_dummy)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  f_mm = phi%f0(work_point, anchor, owner_dummy)
                  g_mm = phi%f1_r(work_point, anchor, owner_dummy)

                  numeric_grad(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
                  do i = 1, ndim
                     numeric_hess(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
                  end do
               end do

               call check(error, analytic_val, numeric_val, thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
               do i = 1, ndim
                  call check(error, analytic_grad(i), numeric_grad(i), &
                             thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic_hess(i, j), numeric_hess(i, j), &
                                thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_phi_f012_r

   !> Test switching function nuclear gradient via finite difference.
   !> Builds an LSF-svdw scaffold to obtain f0 and nuclear gradients,
   !> wraps them through the smooth-step switching function, and
   !> FD-checks against the analytic switching gradient.
   subroutine test_switching_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol_base, mol_shift
      type(moist_cavity_drop_lsf_svdw_type) :: prim
      type(moist_cavity_drop_swif_smooth_step_type) :: sw
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp), allocatable :: centers_base(:, :), centers_local(:, :)
      integer :: icase, ipt, atom, axis
      real(wp) :: point(ndim)
      real(wp) :: lsf0
      real(wp), allocatable :: lsf1(:, :)
      real(wp), allocatable :: analytic(:, :)
      real(wp), allocatable :: dummy_rr_rA(:, :, :, :)
      real(wp) :: f_pp, f_p, f_m, f_mm, numeric
      real(wp) :: lsf0_tmp
      real(wp) :: eps
      type(mctc_error), allocatable :: lsf_err

      eps = STEP_SIZE

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol_base = mols(icase)
         call get_test_radii(mol_base, radii)
         call get_test_points(mol_base, points)
         allocate (centers_base(ndim, mol_base%nat), centers_local(ndim, mol_base%nat))
         centers_base = mol_base%xyz
         call new_swif_smooth_step(sw, -0.5_wp, 0.5_wp)

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol_base, radii)
         call prim%set_max_deriv(3)

         if (allocated(lsf1)) deallocate (lsf1)
         if (allocated(dummy_rr_rA)) deallocate (dummy_rr_rA)
         allocate (lsf1(ndim, mol_base%nat))
         allocate (dummy_rr_rA(ndim, ndim, ndim, mol_base%nat))

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            call prim%update(mol_base, radii)
            call prim%set_centers(centers_base)
            call prim%prepare(point, lsf_err)
            if (allocated(lsf_err)) then
               call test_failed(error, "LSF prepare failed: "//lsf_err%message)
               return
            end if
            call prim%f0(lsf0)
            call prim%f3_rr_rA(lsf1_rA=lsf1, &
                               lsf3_rr_rA=dummy_rr_rA)
            analytic = sw%f1_rA(lsf0, lsf1)

            do atom = 1, mol_base%nat
               do axis = 1, ndim
                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%set_centers(centers_local)
                  call prim%prepare(point, lsf_err)
                  call prim%f0(lsf0_tmp)
                  f_pp = sw%f0(lsf0_tmp)

                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) + eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%set_centers(centers_local)
                  call prim%prepare(point, lsf_err)
                  call prim%f0(lsf0_tmp)
                  f_p = sw%f0(lsf0_tmp)

                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) - eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%set_centers(centers_local)
                  call prim%prepare(point, lsf_err)
                  call prim%f0(lsf0_tmp)
                  f_m = sw%f0(lsf0_tmp)

                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%set_centers(centers_local)
                  call prim%prepare(point, lsf_err)
                  call prim%f0(lsf0_tmp)
                  f_mm = sw%f0(lsf0_tmp)

                  numeric = fd4_scalar(f_pp, f_p, f_m, f_mm, eps)

                  call check(error, &
                             analytic(axis, atom), numeric, &
                             thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_switching_f1_rA

   !* ================================================================================= *!
   !*                          iSwiG switching function                                 *!
   !* ================================================================================= *!

   !> Build the iSwiG switching function of one probe case.
   !>
   !> `use_adj` selects the neighbour-list path: with a Lebedev weight the
   !> constructor builds the sorted adjacency list, without one it destroys it
   !> and every routine falls back to a full loop over all atoms.
   !>
   !> @param[in]  mol     Molecular structure
   !> @param[in]  radii   Atomic radii (bohr)
   !> @param[in]  use_adj Build the adjacency list (`.true.`) or the fallback
   !> @param[out] iswig   Configured switching function
   subroutine iswig_setup(mol, radii, use_adj, iswig)
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Whether to build the sorted adjacency list
      logical, intent(in) :: use_adj
      !> Configured switching function
      type(moist_cavity_drop_iswig), intent(out) :: iswig

      call new_iswig(iswig, iswig_swx)
      if (use_adj) then
         call iswig%update(mol, radii, wleb_max=iswig_wleb_max)
      else
         call iswig%update(mol, radii)
      end if
   end subroutine iswig_setup

   !> Build the probe list of one case: owner atoms, evaluation points and
   !> Lebedev weights.
   !>
   !> Every point is *anchor-like*, `xyz_owner + R_owner w`, because that is the
   !> only thing the centre-based neighbour list can screen correctly: the break
   !> threshold `R_owner + R_max + erf_cutoff/xi` assumes
   !> `|pos - xyz_owner| <= R_owner`. The fallback path carries no such
   !> restriction, but sharing the geometry lets the two paths be compared
   !> directly.
   !>
   !> `w` is placed on the circle where the owner's sphere cuts its nearest
   !> neighbour's, offset by [[iswig_theta_off]], so the probes straddle the
   !> switching transition rather than saturating at zero or one.
   !>
   !> @param[in]  mol    Molecular structure
   !> @param[in]  radii  Atomic radii (bohr)
   !> @param[out] owners Owner atom of each probe (n_iswig_probes)
   !> @param[out] pos    Evaluation point of each probe (3, n_iswig_probes)
   !> @param[out] wleb   Lebedev weight of each probe (n_iswig_probes)
   subroutine iswig_probes(mol, radii, owners, pos, wleb)
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Owner atom of each probe
      integer, allocatable, intent(out) :: owners(:)
      !> Evaluation point of each probe
      real(wp), allocatable, intent(out) :: pos(:, :)
      !> Lebedev weight of each probe
      real(wp), allocatable, intent(out) :: wleb(:)

      integer :: iprobe, iowner, idir, own, near, jat
      real(wp) :: dist, best, uvec(ndim), tvec(ndim), seed(ndim)
      real(wp) :: cos_c, theta, rad_own

      allocate (owners(n_iswig_probes), pos(ndim, n_iswig_probes), wleb(n_iswig_probes))

      do iprobe = 1, n_iswig_probes
         iowner = (iprobe - 1)/n_iswig_dirs + 1
         idir = mod(iprobe - 1, n_iswig_dirs) + 1
         own = mod(iowner - 1, mol%nat) + 1
         owners(iprobe) = own
         rad_own = radii(own)

         ! Nearest other nucleus: the one whose surface the owner's sphere is
         ! most likely to actually cut.
         near = 0
         best = huge(1.0_wp)
         do jat = 1, mol%nat
            if (jat == own) cycle
            dist = norm2(mol%xyz(:, jat) - mol%xyz(:, own))
            if (dist < best) then
               best = dist
               near = jat
            end if
         end do

         if (near == 0) then
            ! Single atom: nothing to cut, any direction will do.
            uvec = [0.0_wp, 0.0_wp, 1.0_wp]
            cos_c = 0.0_wp
         else
            uvec = (mol%xyz(:, near) - mol%xyz(:, own))/best
            ! |pos - xyz_near|**2 = R_own**2 + d**2 - 2 R_own d cos(theta), so the
            ! intersection circle sits at this cosine. Clamp when the spheres do
            ! not meet; the probe is then simply the closest approach.
            cos_c = (rad_own*rad_own + best*best - radii(near)*radii(near)) &
                    /(2.0_wp*rad_own*best)
            cos_c = max(-1.0_wp, min(1.0_wp, cos_c))
         end if

         ! Any unit vector orthogonal to uvec; the seed is deliberately off-axis
         ! so a symmetric fixture cannot make the cross product degenerate.
         seed = [0.3_wp, -0.7_wp, 1.1_wp]
         tvec = seed - dot_product(seed, uvec)*uvec
         if (norm2(tvec) < 1.0e-8_wp) then
            seed = [1.1_wp, 0.3_wp, -0.7_wp]
            tvec = seed - dot_product(seed, uvec)*uvec
         end if
         tvec = tvec/norm2(tvec)

         theta = acos(cos_c) + iswig_theta_off(idir)
         pos(:, iprobe) = mol%xyz(:, own) + rad_own*(cos(theta)*uvec + sin(theta)*tvec)

         ! Spread the widths over the grid's weight range; all stay at or below
         ! `iswig_wleb_max`, which is what the adjacency cutoff was built for.
         wleb(iprobe) = iswig_wleb_max*(0.3_wp + 0.175_wp*real(idir, wp))
      end do
   end subroutine iswig_probes

   !> Compare the switching value and gradient against the committed fixture.
   !>
   !> `test/unit/data/iswig_swi1_golden.txt` was produced by the hand-written
   !> implementation that predates the sparse rows and the generated pair
   !> kernel, so this is what pins the rewrite to the pre-existing mathematics.
   !> Every atom of every probe is compared, not only the nonzero ones: a
   !> spurious nonzero is exactly the failure mode a sparse rewrite invites.
   subroutine test_iswig_swi1_matches_golden(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), grad(:, :)
      integer, allocatable :: owners(:)
      character(len=:), allocatable :: path
      character(len=3) :: tag, want_tag
      character(len=2) :: quantity
      integer :: icase, iprobe, iatom, axis, ipath, unit, stat
      integer :: want_case, want_probe, want_atom, want_axis
      real(wp) :: xi, f0, reference, got
      logical :: use_adj

      path = get_env("MOIST_SOURCE_ROOT", default=".")//"/test/unit/data/iswig_swi1_golden.txt"
      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "cannot open the iSwiG golden fixture: "//path)
         return
      end if
      call skip_comments(unit)

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(grad)) deallocate (grad)
         allocate (grad(ndim, mol%nat))

         do ipath = 1, 2
            use_adj = ipath == 1
            tag = merge("adj", "fbk", use_adj)
            call iswig_setup(mol, radii, use_adj, iswig)
            call work%init(iswig)
            call iswig_probes(mol, radii, owners, pos, wleb)

            do iprobe = 1, n_iswig_probes
               xi = iswig%xi0(owners(iprobe), wleb(iprobe))
               f0 = iswig%swi0(pos(:, iprobe), owners(iprobe), xi)
               read (unit, *, iostat=stat) want_tag, want_case, want_probe, quantity, &
                  want_atom, want_axis, reference
               if (stat /= 0 .or. want_tag /= tag .or. want_case /= icase &
                   .or. want_probe /= iprobe .or. quantity /= "f0") then
                  call test_failed(error, "iSwiG golden fixture is out of step at f0")
                  close (unit)
                  return
               end if
               call check_golden(error, f0, reference, tag, icase, iprobe, "f0")
               if (allocated(error)) then
                  close (unit)
                  return
               end if

               call iswig%swi1_rA(pos(:, iprobe), owners(iprobe), xi, work, grad)
               do iatom = 1, mol%nat
                  do axis = 1, ndim
                     read (unit, *, iostat=stat) want_tag, want_case, want_probe, &
                        quantity, want_atom, want_axis, reference
                     if (stat /= 0 .or. want_atom /= iatom .or. want_axis /= axis &
                         .or. quantity /= "f1") then
                        call test_failed(error, "iSwiG golden fixture is out of step at f1")
                        close (unit)
                        return
                     end if
                     got = grad(axis, iatom)
                     call check_golden(error, got, reference, tag, icase, iprobe, "f1")
                     if (allocated(error)) then
                        close (unit)
                        return
                     end if
                  end do
               end do
            end do
            deallocate (owners, pos, wleb)
            call work%destroy()
         end do
      end do

      close (unit)
   end subroutine test_iswig_swi1_matches_golden

   !> Advance past the fixture's `#` header
   subroutine skip_comments(unit)
      !> Open fixture unit
      integer, intent(in) :: unit

      character(len=1) :: first
      integer :: stat

      do
         read (unit, "(a1)", iostat=stat) first
         if (stat /= 0) return
         if (first /= "#") then
            backspace (unit)
            return
         end if
      end do
   end subroutine skip_comments

   !> Report one golden mismatch with enough context to locate it
   subroutine check_golden(error, got, reference, tag, icase, iprobe, quantity)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Value produced now
      real(wp), intent(in) :: got
      !> Value from the fixture
      real(wp), intent(in) :: reference
      !> Neighbour-list path
      character(len=*), intent(in) :: tag
      !> Case and probe indices
      integer, intent(in) :: icase, iprobe
      !> Quantity name
      character(len=*), intent(in) :: quantity

      character(len=160) :: message

      if (rel_deviation(got, reference) <= iswig_golden_tol) return

      write (message, "(a,1x,a,1x,i0,1x,i0,a,es24.16,a,es24.16)") &
         "iSwiG golden mismatch", tag, icase, iprobe, " "//quantity//": got ", &
         got, " want ", reference
      call test_failed(error, trim(message))
   end subroutine check_golden

   !* --------------------------- finite-difference support ------------------------- *!

   !> Rebuild the switching function on a displaced geometry and evaluate it.
   !>
   !> Walks the joint parameter along `(v, vxi)`: every nucleus by `tau*v`, the
   !> Gaussian width by `tau*vxi`, and the surface point by `tau*v(:, owner)` -
   !> the last because the anchor moves rigidly with its owner, which is the
   !> convention the whole switching row is written in.
   !>
   !> @param[in]  mol     Molecular structure at `tau = 0`
   !> @param[in]  radii   Atomic radii (bohr)
   !> @param[in]  use_adj Whether to build the adjacency list
   !> @param[in]  owner   Owner atom index
   !> @param[in]  pos0    Surface point at `tau = 0`
   !> @param[in]  xi0     Gaussian width at `tau = 0`
   !> @param[in]  v       Nuclear direction (3, nat)
   !> @param[in]  vxi     Width direction
   !> @param[in]  tau     Step along the joint direction
   !> @param[out] f_val   Switching value
   !> @param[out] grad    Dense nuclear gradient (3, nat), optional
   !> @param[out] dxi     Width derivative, optional
   subroutine iswig_at(mol, radii, use_adj, owner, pos0, xi0, v, vxi, tau, f_val, grad, dxi)
      !> Molecular structure at `tau = 0`
      type(structure_type), intent(in) :: mol
      !> Atomic radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Whether to build the adjacency list
      logical, intent(in) :: use_adj
      !> Owner atom index
      integer, intent(in) :: owner
      !> Surface point at `tau = 0`
      real(wp), intent(in) :: pos0(ndim)
      !> Gaussian width at `tau = 0`
      real(wp), intent(in) :: xi0
      !> Nuclear direction
      real(wp), intent(in) :: v(:, :)
      !> Width direction
      real(wp), intent(in) :: vxi
      !> Step along the joint direction
      real(wp), intent(in) :: tau
      !> Switching value
      real(wp), intent(out) :: f_val
      !> Dense nuclear gradient
      real(wp), intent(out), optional :: grad(:, :)
      !> Width derivative
      real(wp), intent(out), optional :: dxi

      type(structure_type) :: moved
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp) :: pos(ndim), xi_t, dxi_local

      moved = mol
      moved%xyz = mol%xyz + tau*v
      pos = pos0 + tau*v(:, owner)
      xi_t = xi0 + tau*vxi

      call iswig_setup(moved, radii, use_adj, iswig)
      f_val = iswig%swi0(pos, owner, xi_t)
      if (present(grad)) then
         call work%init(iswig)
         call iswig%swi1_rA(pos, owner, xi_t, work, grad, dxi_local)
         call work%destroy()
         if (present(dxi)) dxi = dxi_local
      end if
   end subroutine iswig_at

   !> Dense expansion of the sparse contracted second derivative
   !>
   !> @param[in]  iswig  Switching function
   !> @param[in]  pos    Surface point
   !> @param[in]  owner  Owner atom index
   !> @param[in]  xi     Gaussian width
   !> @param[in]  v      Nuclear direction (3, nat)
   !> @param[in]  vxi    Width direction
   !> @param[out] hvp    sum_B v_B . d2 f / (d r_A d r_B), dense (3, nat)
   !> @param[out] dxi2   Directional derivative of the width derivative
   subroutine iswig_hvp_dense(iswig, pos, owner, xi, v, vxi, hvp, dxi2)
      !> Switching function
      type(moist_cavity_drop_iswig), intent(in) :: iswig
      !> Surface point
      real(wp), intent(in) :: pos(ndim)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width
      real(wp), intent(in) :: xi
      !> Nuclear direction
      real(wp), intent(in) :: v(:, :)
      !> Width direction
      real(wp), intent(in) :: vxi
      !> Contracted second derivative, dense
      real(wp), intent(out) :: hvp(:, :)
      !> Directional derivative of the width derivative
      real(wp), intent(out) :: dxi2

      type(iswig_workspace_type) :: work
      real(wp), allocatable :: rows2(:, :)
      real(wp) :: owner_row2(ndim), f_val
      integer :: jj

      call work%init(iswig)
      call iswig%swi_collect(pos, owner, xi, f_val, work)
      allocate (rows2(ndim, max(work%n_nb, 1)))
      call iswig%swi2_rArB_sparse(work, v, rows2, owner_row2, dxi2, vxi=vxi)

      hvp = 0.0_wp
      do jj = 1, work%n_nb
         hvp(:, work%idx(jj)) = rows2(:, jj)
      end do
      hvp(:, owner) = hvp(:, owner) + owner_row2
      call work%destroy()
   end subroutine iswig_hvp_dense

   !> Deterministic pseudo-random direction, seeded from the case indices
   !>
   !> @param[in]  nat  Atom count
   !> @param[in]  seed Seed
   !> @param[out] v    Direction (3, nat), entries in [-1, 1]
   subroutine iswig_direction(nat, seed, v)
      !> Atom count
      integer, intent(in) :: nat
      !> Seed
      integer, intent(in) :: seed
      !> Direction
      real(wp), allocatable, intent(out) :: v(:, :)

      integer :: iatom, axis

      ! A deterministic, obviously reproducible spread. Mutually incommensurate
      ! strides keep the entries from lining up into a symmetry direction of a
      ! symmetric fixture, which would silently weaken every check below.
      allocate (v(ndim, nat))
      do iatom = 1, nat
         do axis = 1, ndim
            v(axis, iatom) = sin(0.7_wp*real(11*seed + 13*iatom + 29*axis, wp))
         end do
      end do
   end subroutine iswig_direction

   !* ------------------------------- iSwiG derivative tests ------------------------ *!

   !> Directional nuclear gradient against finite differences of the value.
   subroutine test_iswig_swi1_rA_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), grad(:, :), v(:, :)
      integer, allocatable :: owners(:)
      integer :: icase, iprobe, ipath, istep, ioff
      real(wp) :: f_val, fv(4), analytic, numeric, err(2), h
      logical :: use_adj

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(grad)) deallocate (grad)
         allocate (grad(ndim, mol%nat))
         call iswig_direction(mol%nat, icase, v)
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            do iprobe = 1, n_iswig_probes
               call iswig_at(mol, radii, use_adj, owners(iprobe), pos(:, iprobe), &
                             iswig_fd_xi, v, 0.0_wp, 0.0_wp, f_val, grad)
               analytic = sum(grad*v)

               do istep = 1, 2
                  h = iswig_fd_steps(istep)
                  do ioff = 1, 4
                     call iswig_at(mol, radii, use_adj, owners(iprobe), pos(:, iprobe), &
                                   iswig_fd_xi, v, 0.0_wp, fd4_offsets(ioff)*h, fv(ioff))
                  end do
                  numeric = fd4_scalar(fv(1), fv(2), fv(3), fv(4), h)
                  err(istep) = abs(analytic - numeric)
               end do

               call check(error, analytic, numeric, &
                          thr_abs=iswig_fd_abs, thr_rel=iswig_fd_rel)
               if (allocated(error)) return
               call check_converges(error, err, "swi1_rA")
               if (allocated(error)) return
            end do
         end do
         deallocate (owners, pos, wleb, v)
      end do
   end subroutine test_iswig_swi1_rA_fd

   !> Width derivative against finite differences of the value.
   !>
   !> `dxi` has no other route to a test: both production callers multiply it by
   !> an identically-zero `xi1_rA`, so the golden fixture cannot see it either.
   subroutine test_iswig_dxi_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), grad(:, :), zero_v(:, :)
      integer, allocatable :: owners(:)
      integer :: icase, iprobe, ipath, istep, ioff
      real(wp) :: f_val, fv(4), analytic, numeric, err(2), h
      logical :: use_adj

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(grad)) deallocate (grad)
         if (allocated(zero_v)) deallocate (zero_v)
         allocate (grad(ndim, mol%nat))
         allocate (zero_v(ndim, mol%nat), source=0.0_wp)
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            do iprobe = 1, n_iswig_probes
               call iswig_at(mol, radii, use_adj, owners(iprobe), pos(:, iprobe), &
                             iswig_fd_xi, zero_v, 0.0_wp, 0.0_wp, f_val, grad, analytic)

               do istep = 1, 2
                  h = iswig_fd_steps(istep)
                  do ioff = 1, 4
                     call iswig_at(mol, radii, use_adj, owners(iprobe), pos(:, iprobe), &
                                   iswig_fd_xi, zero_v, 1.0_wp, fd4_offsets(ioff)*h, fv(ioff))
                  end do
                  numeric = fd4_scalar(fv(1), fv(2), fv(3), fv(4), h)
                  err(istep) = abs(analytic - numeric)
               end do

               call check(error, analytic, numeric, &
                          thr_abs=iswig_fd_abs, thr_rel=iswig_fd_rel)
               if (allocated(error)) return
               call check_converges(error, err, "dxi")
               if (allocated(error)) return
            end do
         end do
         deallocate (owners, pos, wleb)
      end do
   end subroutine test_iswig_dxi_fd

   !> Contracted second derivative against finite differences of the gradient.
   subroutine test_iswig_swi2_rArB_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_iswig_hvp_fd(error, with_xi=.false.)
   end subroutine test_iswig_swi2_rArB_fd

   !> The width channel of the contracted second derivative.
   !>
   !> Drives the joint direction with `v = 0` and `vxi = 1`, so what is under
   !> test is the mixed `d2 f / (d r_A d xi)` row and the `d2 f / d xi2` term -
   !> both of which are zero in every production path today, and both of which
   !> become live with geometry-dependent radii.
   subroutine test_iswig_swi2_xi_channel_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_iswig_hvp_fd(error, with_xi=.true.)
   end subroutine test_iswig_swi2_xi_channel_fd

   !> Shared driver of the two Hessian-vector finite-difference tests.
   !>
   !> @param[out] error   Error handle
   !> @param[in]  with_xi Drive the width channel instead of the nuclear one
   subroutine run_iswig_hvp_fd(error, with_xi)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Drive the width channel instead of the nuclear one
      logical, intent(in) :: with_xi

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), v(:, :)
      real(wp), allocatable :: hvp(:, :), gv(:, :, :), numeric(:, :)
      integer, allocatable :: owners(:)
      integer :: icase, iprobe, ipath, istep, ioff, iatom, axis
      real(wp) :: f_val, dxi2, dxi_v(4), dxi2_num, vxi, h, err(2), worst
      logical :: use_adj

      vxi = merge(1.0_wp, 0.0_wp, with_xi)

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(hvp)) deallocate (hvp, gv, numeric)
         allocate (hvp(ndim, mol%nat), numeric(ndim, mol%nat))
         allocate (gv(ndim, mol%nat, 4))
         call iswig_direction(mol%nat, icase + 977, v)
         if (with_xi) v = 0.0_wp
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)

            do iprobe = 1, n_iswig_probes
               call iswig_hvp_dense(iswig, pos(:, iprobe), owners(iprobe), iswig_fd_xi, &
                                    v, vxi, hvp, dxi2)

               do istep = 1, 2
                  h = iswig_fd_steps(istep)
                  do ioff = 1, 4
                     call iswig_at(mol, radii, use_adj, owners(iprobe), pos(:, iprobe), &
                                   iswig_fd_xi, v, vxi, fd4_offsets(ioff)*h, f_val, &
                                   gv(:, :, ioff), dxi_v(ioff))
                  end do
                  worst = 0.0_wp
                  do iatom = 1, mol%nat
                     do axis = 1, ndim
                        numeric(axis, iatom) = fd4_scalar(gv(axis, iatom, 1), &
                                                          gv(axis, iatom, 2), &
                                                          gv(axis, iatom, 3), &
                                                          gv(axis, iatom, 4), h)
                        worst = max(worst, abs(hvp(axis, iatom) - numeric(axis, iatom)))
                     end do
                  end do
                  dxi2_num = fd4_scalar(dxi_v(1), dxi_v(2), dxi_v(3), dxi_v(4), h)
                  err(istep) = max(worst, abs(dxi2 - dxi2_num))
               end do

               do iatom = 1, mol%nat
                  do axis = 1, ndim
                     call check(error, hvp(axis, iatom), numeric(axis, iatom), &
                                thr_abs=iswig_fd_abs, thr_rel=iswig_fd_rel)
                     if (allocated(error)) return
                  end do
               end do
               call check(error, dxi2, dxi2_num, thr_abs=iswig_fd_abs, thr_rel=iswig_fd_rel)
               if (allocated(error)) return
               call check_converges(error, err, "swi2_rArB")
               if (allocated(error)) return
            end do
         end do
         deallocate (owners, pos, wleb, v)
      end do
   end subroutine run_iswig_hvp_fd

   !> Assert that halving the step did not make the finite difference worse.
   !>
   !> A fourth-order stencil should improve by 16x until roundoff takes over, so
   !> this only rejects an error that *grows* - the signature of a derivative
   !> that is wrong rather than merely truncated.
   !>
   !> @param[out] error Error handle
   !> @param[in]  err   Errors at the coarse and fine steps
   !> @param[in]  what  Name for the failure message
   subroutine check_converges(error, err, what)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Errors at the coarse and fine steps
      real(wp), intent(in) :: err(2)
      !> Name for the failure message
      character(len=*), intent(in) :: what

      character(len=160) :: message

      ! Below this the comparison is measuring roundoff, not truncation.
      if (err(2) < 1.0e-11_wp) return
      if (err(2) <= err(1)) return

      write (message, "(a,a,a,es12.4,a,es12.4)") "iSwiG ", what, &
         " finite difference diverges: ", err(1), " -> ", err(2)
      call test_failed(error, trim(message))
   end subroutine check_converges

   !> Symmetry of the contracted second derivative, `v^T H u == u^T H v`.
   !>
   !> Over the joint `(positions, xi)` parameter, so the mixed row and column
   !> are under test too. Not exact: the two products run different arithmetic -
   !> `dr_k`, `s` and `dl` all differ per direction - so this is a tolerance
   !> check, unlike the null mode below.
   subroutine test_iswig_swi2_symmetry(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), v(:, :), u(:, :)
      real(wp), allocatable :: hu(:, :), hv(:, :)
      integer, allocatable :: owners(:)
      integer :: icase, iprobe, ipath
      real(wp) :: hu_xi, hv_xi, vhu, uhv, vxi, uxi
      logical :: use_adj

      vxi = 0.37_wp
      uxi = -0.62_wp

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(hu)) deallocate (hu, hv)
         allocate (hu(ndim, mol%nat), hv(ndim, mol%nat))
         call iswig_direction(mol%nat, icase + 31, v)
         call iswig_direction(mol%nat, icase + 617, u)
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)

            do iprobe = 1, n_iswig_probes
               call iswig_hvp_dense(iswig, pos(:, iprobe), owners(iprobe), iswig_fd_xi, &
                                    u, uxi, hu, hu_xi)
               call iswig_hvp_dense(iswig, pos(:, iprobe), owners(iprobe), iswig_fd_xi, &
                                    v, vxi, hv, hv_xi)
               vhu = sum(v*hu) + vxi*hu_xi
               uhv = sum(u*hv) + uxi*hv_xi
               call check(error, vhu, uhv, thr_abs=1.0e-10_wp, thr_rel=1.0e-10_wp)
               if (allocated(error)) return
            end do
         end do
         deallocate (owners, pos, wleb, v, u)
      end do
   end subroutine test_iswig_swi2_symmetry

   !> A uniform translation is a null mode of both derivative levels.
   !>
   !> The second derivative vanishes *exactly*: a uniform direction makes every
   !> `w_k = v_owner - v_k` identically zero, so `dr_k`, `s` and every `M_k w_k`
   !> are true zeros rather than cancelling sums. The first derivative only
   !> cancels to roundoff, because `sum_A grad_A . v` sums the same terms in a
   !> different order than the owner row accumulated them.
   subroutine test_iswig_translation_null(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), v(:, :), hvp(:, :), grad(:, :)
      integer, allocatable :: owners(:)
      integer :: icase, iprobe, ipath, iatom, axis
      real(wp) :: dxi2, dxi
      logical :: use_adj

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(hvp)) deallocate (hvp, grad, v)
         allocate (hvp(ndim, mol%nat), grad(ndim, mol%nat), v(ndim, mol%nat))
         v(1, :) = 0.3_wp
         v(2, :) = -0.7_wp
         v(3, :) = 1.1_wp
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)
            call work%init(iswig)

            do iprobe = 1, n_iswig_probes
               call iswig%swi1_rA(pos(:, iprobe), owners(iprobe), iswig_fd_xi, work, &
                                  grad, dxi)
               call check(error, sum(grad*v), 0.0_wp, thr_abs=1.0e-12_wp, thr_rel=0.0_wp)
               if (allocated(error)) return

               call iswig_hvp_dense(iswig, pos(:, iprobe), owners(iprobe), iswig_fd_xi, &
                                    v, 0.0_wp, hvp, dxi2)
               do iatom = 1, mol%nat
                  do axis = 1, ndim
                     call check(error, hvp(axis, iatom), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                                more="iSwiG translation hvp is not exactly zero")
                     if (allocated(error)) return
                  end do
               end do
               call check(error, dxi2, 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                          more="iSwiG translation dxi2 is not exactly zero")
               if (allocated(error)) return
            end do
            call work%destroy()
         end do
         deallocate (owners, pos, wleb)
      end do
   end subroutine test_iswig_translation_null

   !> The isolated-atom and `adj_list%n == 0` corners.
   !>
   !> Two atoms 50 bohr apart put the owner outside every cutoff, so the
   !> adjacency path caches no neighbours at all. That is the `n_nb == 0` corner
   !> the sparse rows have to survive: value one, no rows, and both derivative
   !> levels exactly zero.
   !>
   !> The same fixture also pins the two traversals against each other. They are
   !> not identical by construction - the fallback deliberately carries no
   !> per-atom erf cutoff, a pre-existing asymmetry this refactor preserves - so
   !> the comparison is to a tolerance, on a compact structure where the cutoff
   !> is not what decides the answer.
   subroutine test_iswig_isolated_and_fallback(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp) :: radii(2), pos(ndim), grad(ndim, 2), hvp(ndim, 2), v(ndim, 2)
      real(wp) :: rows(ndim, 2), owner_row(ndim), f_val, dxi, dxi2
      integer :: ipath, iatom, axis
      logical :: use_adj

      call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     50.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      radii = [3.0_wp, 3.0_wp]
      pos = [3.0_wp, 0.0_wp, 0.0_wp]
      v(:, 1) = [0.3_wp, -0.7_wp, 1.1_wp]
      v(:, 2) = [-0.9_wp, 0.4_wp, 0.2_wp]

      do ipath = 1, 2
         use_adj = ipath == 1
         call iswig_setup(mol, radii, use_adj, iswig)
         call work%init(iswig)
         call iswig%swi_collect(pos, 1, iswig_fd_xi, f_val, work)

         ! The far neighbour is out of the adjacency list entirely; on the
         ! fallback path it is cached but contributes an exactly unit factor.
         if (use_adj .and. work%n_nb /= 0) then
            call test_failed(error, "isolated owner still cached a neighbour")
            return
         end if
         call check(error, f_val, 1.0_wp, thr_abs=1.0e-14_wp, thr_rel=0.0_wp)
         if (allocated(error)) return

         call iswig%swi1_rA_sparse(work, rows, owner_row, dxi)
         call iswig%swi2_rArB_sparse(work, v, rows, owner_row, dxi2, vxi=0.5_wp)
         call check(error, dxi, 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                    more="iSwiG isolated dxi is not exactly zero")
         if (allocated(error)) return
         call check(error, dxi2, 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                    more="iSwiG isolated dxi2 is not exactly zero")
         if (allocated(error)) return

         call iswig%swi1_rA(pos, 1, iswig_fd_xi, work, grad)
         call iswig_hvp_dense(iswig, pos, 1, iswig_fd_xi, v, 0.5_wp, hvp, dxi2)
         do iatom = 1, 2
            do axis = 1, ndim
               call check(error, grad(axis, iatom), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                          more="iSwiG isolated gradient is not exactly zero")
               if (allocated(error)) return
               call check(error, hvp(axis, iatom), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                          more="iSwiG isolated hvp is not exactly zero")
               if (allocated(error)) return
            end do
         end do
         call work%destroy()
      end do

      call check_paths_agree(error)
   end subroutine test_iswig_isolated_and_fallback

   !> The two traversals must agree on a compact structure, where no cutoff is
   !> deciding the answer.
   subroutine check_paths_agree(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: adj, fbk
      type(iswig_workspace_type) :: work_a, work_f
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), ga(:, :), gf(:, :)
      integer, allocatable :: owners(:)
      integer :: icase, iprobe, iatom, axis
      real(wp) :: dxi_a, dxi_f

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(ga)) deallocate (ga, gf)
         allocate (ga(ndim, mol%nat), gf(ndim, mol%nat))
         call iswig_setup(mol, radii, .true., adj)
         call iswig_setup(mol, radii, .false., fbk)
         call work_a%init(adj)
         call work_f%init(fbk)
         call iswig_probes(mol, radii, owners, pos, wleb)

         do iprobe = 1, n_iswig_probes
            call adj%swi1_rA(pos(:, iprobe), owners(iprobe), iswig_fd_xi, work_a, ga, dxi_a)
            call fbk%swi1_rA(pos(:, iprobe), owners(iprobe), iswig_fd_xi, work_f, gf, dxi_f)
            call check(error, dxi_a, dxi_f, thr_abs=1.0e-10_wp, thr_rel=1.0e-10_wp)
            if (allocated(error)) return
            do iatom = 1, mol%nat
               do axis = 1, ndim
                  call check(error, ga(axis, iatom), gf(axis, iatom), &
                             thr_abs=1.0e-10_wp, thr_rel=1.0e-10_wp)
                  if (allocated(error)) return
               end do
            end do
         end do
         call work_a%destroy()
         call work_f%destroy()
         deallocate (owners, pos, wleb)
      end do
   end subroutine check_paths_agree

   !* --------------------- iSwiG local second-derivative block ---------------------- *!

   !> Assemble the local second-derivative block of one probe.
   !>
   !> The block routine allocates nothing, so the arrays are sized here from
   !> `work%capacity + 1`, which bounds the influence set whatever the traversal
   !> ends up caching. Only `1:n` of each is ever read back.
   !>
   !> @param[in]    iswig Switching function
   !> @param[in]    pos   Surface point
   !> @param[in]    owner Owner atom index
   !> @param[in]    xi    Gaussian width
   !> @param[inout] work  Caller-owned neighbour cache, sized by `init`
   !> @param[out]   n     Influence-set size
   !> @param[out]   idx   Atom ids of the influence set, owner first
   !> @param[out]   blk   Position-position block over the influence set
   !> @param[out]   mix   Mixed position-width rows
   !> @param[out]   d2xi  Second derivative w.r.t. the width
   subroutine iswig_block_at(iswig, pos, owner, xi, work, n, idx, blk, mix, d2xi)
      !> Switching function
      type(moist_cavity_drop_iswig), intent(in) :: iswig
      !> Surface point
      real(wp), intent(in) :: pos(ndim)
      !> Owner atom index
      integer, intent(in) :: owner
      !> Gaussian width
      real(wp), intent(in) :: xi
      !> Caller-owned neighbour cache
      type(iswig_workspace_type), intent(inout) :: work
      !> Influence-set size
      integer, intent(out) :: n
      !> Atom ids of the influence set
      integer, allocatable, intent(out) :: idx(:)
      !> Position-position block
      real(wp), allocatable, intent(out) :: blk(:, :, :, :)
      !> Mixed position-width rows
      real(wp), allocatable, intent(out) :: mix(:, :)
      !> Second derivative w.r.t. the width
      real(wp), intent(out) :: d2xi

      real(wp) :: f_val
      integer :: cap

      call iswig%swi_collect(pos, owner, xi, f_val, work)
      cap = work%capacity + 1
      allocate (idx(cap), blk(ndim, cap, ndim, cap), mix(ndim, cap))
      call iswig%swi2_rArB_block(work, n, idx, blk, mix, d2xi)
   end subroutine iswig_block_at

   !> Contract a local block over its influence set against a dense direction.
   !>
   !> The generic contraction, in the natural index order; the translation test
   !> below deliberately does not use it, and says why.
   !>
   !> @param[in]  n    Influence-set size
   !> @param[in]  idx  Atom ids of the influence set
   !> @param[in]  blk  Position-position block
   !> @param[in]  mix  Mixed position-width rows
   !> @param[in]  d2xi Second derivative w.r.t. the width
   !> @param[in]  v    Nuclear direction (3, nat)
   !> @param[in]  vxi  Width direction
   !> @param[out] hrow Contracted rows over the influence set (3, >= n)
   !> @param[out] dxi2 Directional derivative of the width derivative
   subroutine iswig_block_contract(n, idx, blk, mix, d2xi, v, vxi, hrow, dxi2)
      !> Influence-set size
      integer, intent(in) :: n
      !> Atom ids of the influence set
      integer, intent(in) :: idx(:)
      !> Position-position block
      real(wp), intent(in) :: blk(:, :, :, :)
      !> Mixed position-width rows
      real(wp), intent(in) :: mix(:, :)
      !> Second derivative w.r.t. the width
      real(wp), intent(in) :: d2xi
      !> Nuclear direction
      real(wp), intent(in) :: v(:, :)
      !> Width direction
      real(wp), intent(in) :: vxi
      !> Contracted rows over the influence set
      real(wp), intent(out) :: hrow(:, :)
      !> Directional derivative of the width derivative
      real(wp), intent(out) :: dxi2

      integer :: ia, ib

      dxi2 = d2xi*vxi
      do ia = 1, n
         hrow(:, ia) = mix(:, ia)*vxi
         do ib = 1, n
            hrow(:, ia) = hrow(:, ia) + matmul(blk(:, ia, :, ib), v(:, idx(ib)))
         end do
         dxi2 = dxi2 + dot_product(mix(:, ia), v(:, idx(ia)))
      end do
   end subroutine iswig_block_contract

   !> The block and the sparse contraction are the same object by two routes.
   !>
   !> `swi2_rArB_sparse` collapses the rank-one term into one scalar per
   !> direction; the block never forms that scalar and multiplies the stored
   !> outer product instead. Agreement is therefore tight but not bitwise - the
   !> two factor the same algebra differently, and the compiler is free to
   !> contract either into an FMA.
   subroutine test_iswig_swi2_block_matches_sparse(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:), v(:, :)
      real(wp), allocatable :: blk(:, :, :, :), mix(:, :), rows2(:, :), hrow(:, :)
      integer, allocatable :: owners(:), idx(:)
      integer :: icase, iprobe, ipath, idir, jj, n
      real(wp) :: owner_row2(ndim), d2xi, dxi2_ref, dxi2_blk
      logical :: use_adj

      !> Width directions probed alongside the nuclear ones
      real(wp), parameter :: vxi_set(3) = [0.0_wp, 0.83_wp, -1.4_wp]

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(rows2)) deallocate (rows2, hrow)
         allocate (rows2(ndim, mol%nat), hrow(ndim, mol%nat + 1))
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)
            call work%init(iswig)

            do iprobe = 1, n_iswig_probes
               do idir = 1, size(vxi_set)
                  call iswig_direction(mol%nat, icase + 53*idir, v)
                  call iswig_block_at(iswig, pos(:, iprobe), owners(iprobe), &
                                      iswig_fd_xi, work, n, idx, blk, mix, d2xi)
                  call iswig%swi2_rArB_sparse(work, v, rows2, owner_row2, dxi2_ref, &
                                              vxi=vxi_set(idir))
                  call iswig_block_contract(n, idx, blk, mix, d2xi, v, vxi_set(idir), &
                                            hrow, dxi2_blk)

                  call check_close(error, hrow(:, 1), owner_row2, "block owner row")
                  if (allocated(error)) return
                  do jj = 1, work%n_nb
                     call check_close(error, hrow(:, 1 + jj), rows2(:, jj), &
                                      "block neighbour row")
                     if (allocated(error)) return
                  end do
                  call check(error, dxi2_blk, dxi2_ref, &
                             thr_abs=iswig_block_tol, thr_rel=iswig_block_tol)
                  if (allocated(error)) return
               end do
            end do
            call work%destroy()
         end do
         deallocate (owners, pos, wleb)
      end do
   end subroutine test_iswig_swi2_block_matches_sparse

   !> The block is symmetric to the bit, not merely to roundoff.
   !>
   !> Only the upper triangle over `(A, B)` is computed and the rest is a copy,
   !> so `blk(i, A, j, B) == blk(j, B, i, A)` holds exactly - including inside
   !> the diagonal blocks, where the two entries would otherwise be the same
   !> product formed in the other order.
   subroutine test_iswig_swi2_block_symmetry(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:)
      real(wp), allocatable :: blk(:, :, :, :), mix(:, :)
      integer, allocatable :: owners(:), idx(:)
      integer :: icase, iprobe, ipath, ia, ib, i, j, n
      real(wp) :: d2xi
      logical :: use_adj
      character(len=160) :: message

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)
            call work%init(iswig)

            do iprobe = 1, n_iswig_probes
               call iswig_block_at(iswig, pos(:, iprobe), owners(iprobe), &
                                   iswig_fd_xi, work, n, idx, blk, mix, d2xi)
               do ib = 1, n
                  do ia = 1, n
                     do j = 1, ndim
                        do i = 1, ndim
                           if (blk(i, ia, j, ib) == blk(j, ib, i, ia)) cycle
                           write (message, "(a,4(1x,i0),a,es24.16,a,es24.16)") &
                              "iSwiG block is not symmetric at", i, ia, j, ib, &
                              ": ", blk(i, ia, j, ib), " vs ", blk(j, ib, i, ia)
                           call test_failed(error, trim(message))
                           return
                        end do
                     end do
                  end do
               end do
            end do
            call work%destroy()
         end do
         deallocate (owners, pos, wleb)
      end do
   end subroutine test_iswig_swi2_block_symmetry

   !> A uniform translation is a null mode of the block as well.
   !>
   !> The owner strip is the running negation of the neighbour blocks, so every
   !> neighbour row sums to the exact zero 3x3 - provided it is summed in that
   !> same order, neighbours first and the owner last, which is why the
   !> contraction here is written by hand rather than through
   !> [[iswig_block_contract]]. Summing the matrix-vector products instead would
   !> re-associate the cancellation and leave roundoff, exactly as
   !> [[test_iswig_translation_null]] observes for the first derivative.
   !>
   !> The owner's own row is the one place the block cannot be both bitwise
   !> symmetric and bitwise translation-invariant, so it is checked to
   !> roundoff; the routine's header explains which of the two it keeps.
   subroutine test_iswig_swi2_block_translation(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:)
      real(wp), allocatable :: blk(:, :, :, :), mix(:, :)
      real(wp), allocatable :: v(:, :), rows2(:, :)
      integer, allocatable :: owners(:), idx(:)
      integer :: icase, iprobe, ipath, ia, ib, jj, axis, n
      real(wp) :: d2xi, acc(ndim, ndim), macc(ndim), hvp(ndim)
      real(wp) :: owner_row2(ndim), dxi2_sp
      logical :: use_adj

      !> The uniform direction every atom is displaced along
      real(wp), parameter :: tvec(ndim) = [0.3_wp, -0.7_wp, 1.1_wp]

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call iswig_probes(mol, radii, owners, pos, wleb)

         ! The same uniform displacement, as a dense direction for the
         ! contracted routine to annihilate.
         if (allocated(v)) deallocate (v, rows2)
         allocate (v(ndim, mol%nat), rows2(ndim, mol%nat))
         do ia = 1, mol%nat
            v(:, ia) = tvec
         end do

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)
            call work%init(iswig)

            do iprobe = 1, n_iswig_probes
               call iswig_block_at(iswig, pos(:, iprobe), owners(iprobe), &
                                   iswig_fd_xi, work, n, idx, blk, mix, d2xi)

               do ia = 1, n
                  acc(:, :) = 0.0_wp
                  do ib = 2, n
                     acc(:, :) = acc(:, :) + blk(:, ia, :, ib)
                  end do
                  acc(:, :) = acc(:, :) + blk(:, ia, :, 1)
                  hvp = matmul(acc, tvec)
                  do axis = 1, ndim
                     if (ia == 1) then
                        call check(error, hvp(axis), 0.0_wp, &
                                   thr_abs=1.0e-12_wp, thr_rel=0.0_wp)
                     else
                        call check(error, hvp(axis), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                                   more="iSwiG block translation row is not exactly zero")
                     end if
                     if (allocated(error)) return
                  end do
               end do

               macc(:) = 0.0_wp
               do ia = 2, n
                  macc(:) = macc(:) + mix(:, ia)
               end do
               macc(:) = macc(:) + mix(:, 1)
               do axis = 1, ndim
                  call check(error, macc(axis), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                             more="iSwiG block translation mix is not exactly zero")
                  if (allocated(error)) return
               end do
               ! `swi2_rArB_sparse` has to annihilate the same translation,
               ! width channel included: `wvec` vanishes for every neighbour, so
               ! `s_sum` and `t_sum` are exact zeros and `vxi` is zero here by
               ! construction. Unlike the block, its owner row is the running
               ! negation, so it cancels exactly too.
               call iswig%swi2_rArB_sparse(work, v, rows2, owner_row2, dxi2_sp, &
                                           vxi=0.0_wp)
               do jj = 1, work%n_nb
                  do axis = 1, ndim
                     call check(error, rows2(axis, jj), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                                more="iSwiG sparse translation row is not exactly zero")
                     if (allocated(error)) return
                  end do
               end do
               do axis = 1, ndim
                  call check(error, owner_row2(axis), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                             more="iSwiG sparse translation owner row is not exactly zero")
                  if (allocated(error)) return
               end do
               call check(error, dxi2_sp, 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                          more="iSwiG sparse translation dxi2 is not exactly zero")
               if (allocated(error)) return
            end do
            call work%destroy()
         end do
         deallocate (owners, pos, wleb)
      end do
   end subroutine test_iswig_swi2_block_translation

   !> The width channel against finite differences of the first derivative.
   !>
   !> `mix` is `d2 f / (d r_A d xi)` and `d2xi` is `d2 f / d xi2`, so both are
   !> derivatives of [[iswig_swi_f1_rA_sparse]]'s output w.r.t. the width alone,
   !> at frozen geometry. The comparison runs in the dense per-atom layout
   !> because the cached neighbour set - and with it the compact index space -
   !> may change across the stencil.
   subroutine test_iswig_swi2_block_xi_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: radii(:), pos(:, :), wleb(:)
      real(wp), allocatable :: blk(:, :, :, :), mix(:, :)
      real(wp), allocatable :: dmix(:, :), gv(:, :, :), numeric(:, :)
      integer, allocatable :: owners(:), idx(:)
      integer :: icase, iprobe, ipath, istep, ioff, iatom, axis, ia, n
      real(wp) :: d2xi, dxi_v(4), dxi2_num, h, err(2), worst
      logical :: use_adj

      call get_test_structures(mols, n_iswig_structures)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         if (allocated(dmix)) deallocate (dmix, gv, numeric)
         allocate (dmix(ndim, mol%nat), numeric(ndim, mol%nat))
         allocate (gv(ndim, mol%nat, 4))
         call iswig_probes(mol, radii, owners, pos, wleb)

         do ipath = 1, 2
            use_adj = ipath == 1
            call iswig_setup(mol, radii, use_adj, iswig)
            call work%init(iswig)

            do iprobe = 1, n_iswig_probes
               call iswig_block_at(iswig, pos(:, iprobe), owners(iprobe), &
                                   iswig_fd_xi, work, n, idx, blk, mix, d2xi)
               dmix = 0.0_wp
               do ia = 1, n
                  dmix(:, idx(ia)) = mix(:, ia)
               end do

               do istep = 1, 2
                  h = iswig_fd_steps(istep)
                  do ioff = 1, 4
                     call iswig%swi1_rA(pos(:, iprobe), owners(iprobe), &
                                        iswig_fd_xi + fd4_offsets(ioff)*h, work, &
                                        gv(:, :, ioff), dxi_v(ioff))
                  end do
                  worst = 0.0_wp
                  do iatom = 1, mol%nat
                     do axis = 1, ndim
                        numeric(axis, iatom) = fd4_scalar(gv(axis, iatom, 1), &
                                                          gv(axis, iatom, 2), &
                                                          gv(axis, iatom, 3), &
                                                          gv(axis, iatom, 4), h)
                        worst = max(worst, abs(dmix(axis, iatom) - numeric(axis, iatom)))
                     end do
                  end do
                  dxi2_num = fd4_scalar(dxi_v(1), dxi_v(2), dxi_v(3), dxi_v(4), h)
                  err(istep) = max(worst, abs(d2xi - dxi2_num))
               end do

               do iatom = 1, mol%nat
                  do axis = 1, ndim
                     call check(error, dmix(axis, iatom), numeric(axis, iatom), &
                                thr_abs=iswig_fd_abs, thr_rel=iswig_fd_rel)
                     if (allocated(error)) return
                  end do
               end do
               call check(error, d2xi, dxi2_num, thr_abs=iswig_fd_abs, thr_rel=iswig_fd_rel)
               if (allocated(error)) return
               call check_converges(error, err, "swi2 block xi channel")
               if (allocated(error)) return
            end do
            call work%destroy()
         end do
         deallocate (owners, pos, wleb)
      end do
   end subroutine test_iswig_swi2_block_xi_fd

   !> A saturated or degenerate neighbour leaves exact zeros, not small ones.
   !>
   !> Three corners of the guarded pair interface, in one fixture family:
   !>
   !>   * an owner with no cached neighbour at all - the influence set is the
   !>     owner alone and the whole block vanishes;
   !>   * a neighbour 50 bohr away, which the fallback traversal caches with
   !>     every log-derivative coefficient exactly zero - its rows and columns
   !>     of the block must be true zeros while a live neighbour's are not;
   !>   * a point sitting exactly on a neighbour's centre, the slot
   !>     [[iswig_swi_collect]] zeroes wholesale because no row can be formed
   !>     for it.
   subroutine test_iswig_swi2_block_guarded(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_iswig) :: iswig
      type(iswig_workspace_type) :: work
      real(wp), allocatable :: blk(:, :, :, :), mix(:, :)
      integer, allocatable :: idx(:)
      real(wp) :: radii(3), pos(ndim), d2xi, biggest
      integer :: ipath, ia, i, j, n, far, live
      logical :: use_adj

      call new(mol, [8, 8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                        3.6_wp, 0.0_wp, 0.0_wp, &
                                        50.0_wp, 0.0_wp, 0.0_wp], [3, 3]))
      radii = [3.0_wp, 3.0_wp, 3.0_wp]
      ! On the owner's sphere and inside atom 2's switching shell, so atom 2 is
      ! live while atom 3 saturates.
      pos = [2.4_wp, 1.8_wp, 0.0_wp]

      do ipath = 1, 2
         use_adj = ipath == 1
         call iswig_setup(mol, radii, use_adj, iswig)
         call work%init(iswig)
         call iswig_block_at(iswig, pos, 1, iswig_fd_xi, work, n, idx, blk, mix, d2xi)

         far = influence_slot(n, idx, 3)
         live = influence_slot(n, idx, 2)
         if (live == 0) then
            call test_failed(error, "the live neighbour left the influence set")
            return
         end if

         ! The saturated neighbour, when it is cached at all, may not touch a
         ! single entry of the block.
         if (far > 0) then
            do ia = 1, n
               do j = 1, ndim
                  do i = 1, ndim
                     call check(error, blk(i, far, j, ia), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                                more="iSwiG saturated row is not exactly zero")
                     if (allocated(error)) return
                     call check(error, blk(i, ia, j, far), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                                more="iSwiG saturated column is not exactly zero")
                     if (allocated(error)) return
                  end do
               end do
            end do
            do i = 1, ndim
               call check(error, mix(i, far), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                          more="iSwiG saturated mix is not exactly zero")
               if (allocated(error)) return
            end do
         end if

         ! ... and the live one must not be zero, or the check above is vacuous.
         biggest = 0.0_wp
         do ia = 1, n
            do j = 1, ndim
               do i = 1, ndim
                  biggest = max(biggest, abs(blk(i, live, j, ia)))
               end do
            end do
         end do
         if (biggest == 0.0_wp) then
            call test_failed(error, "the live neighbour contributed nothing")
            return
         end if
         call work%destroy()
      end do

      ! An owner 50 bohr from anything: the adjacency traversal caches nothing,
      ! so the influence set is the owner alone and every output vanishes.
      call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     50.0_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call iswig_setup(mol, radii(:2), .true., iswig)
      call work%init(iswig)
      call iswig_block_at(iswig, [3.0_wp, 0.0_wp, 0.0_wp], 1, iswig_fd_xi, work, &
                          n, idx, blk, mix, d2xi)
      if (n /= 1) then
         call test_failed(error, "isolated owner still cached a neighbour")
         return
      end if
      call check_block_zero(error, n, blk, mix, d2xi, "isolated block")
      if (allocated(error)) return
      call work%destroy()

      ! The point sits exactly on the neighbour's centre: `collect` zeroes that
      ! slot wholesale, because every coefficient of it is a log-derivative.
      call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     3.6_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call iswig_setup(mol, [3.6_wp, 3.0_wp], .true., iswig)
      call work%init(iswig)
      call iswig_block_at(iswig, mol%xyz(:, 2), 1, iswig_fd_xi, work, &
                          n, idx, blk, mix, d2xi)
      if (n /= 2) then
         call test_failed(error, "the degenerate neighbour was not cached")
         return
      end if
      call check_block_zero(error, n, blk, mix, d2xi, "degenerate block")
      if (allocated(error)) return
      call work%destroy()
   end subroutine test_iswig_swi2_block_guarded

   !> Position of an atom in the influence set, or zero when it is absent
   !>
   !> @param[in] n     Influence-set size
   !> @param[in] idx   Atom ids of the influence set
   !> @param[in] iatom Atom to locate
   !> @return    slot  Index into the influence set, or zero
   pure function influence_slot(n, idx, iatom) result(slot)
      !> Influence-set size
      integer, intent(in) :: n
      !> Atom ids of the influence set
      integer, intent(in) :: idx(:)
      !> Atom to locate
      integer, intent(in) :: iatom
      !> Index into the influence set, or zero
      integer :: slot

      integer :: ia

      slot = 0
      do ia = 1, n
         if (idx(ia) == iatom) slot = ia
      end do
   end function influence_slot

   !> Assert that a whole block, its width rows and its width curvature vanish
   !>
   !> @param[out] error Error handle
   !> @param[in]  n     Influence-set size
   !> @param[in]  blk   Position-position block
   !> @param[in]  mix   Mixed position-width rows
   !> @param[in]  d2xi  Second derivative w.r.t. the width
   !> @param[in]  what  Name for the failure message
   subroutine check_block_zero(error, n, blk, mix, d2xi, what)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Influence-set size
      integer, intent(in) :: n
      !> Position-position block
      real(wp), intent(in) :: blk(:, :, :, :)
      !> Mixed position-width rows
      real(wp), intent(in) :: mix(:, :)
      !> Second derivative w.r.t. the width
      real(wp), intent(in) :: d2xi
      !> Name for the failure message
      character(len=*), intent(in) :: what

      integer :: ia, ib, i, j

      do ib = 1, n
         do ia = 1, n
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, blk(i, ia, j, ib), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                             more="iSwiG "//what//" is not exactly zero")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
      do ia = 1, n
         do i = 1, ndim
            call check(error, mix(i, ia), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                       more="iSwiG "//what//" is not exactly zero")
            if (allocated(error)) return
         end do
      end do
      call check(error, d2xi, 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp, &
                 more="iSwiG "//what//" is not exactly zero")
   end subroutine check_block_zero

   !> Compare two 3-vectors at the block's agreement tolerance
   !>
   !> @param[out] error Error handle
   !> @param[in]  got   Vector produced by the block route
   !> @param[in]  want  Vector produced by the sparse route
   !> @param[in]  what  Name for the failure message
   subroutine check_close(error, got, want, what)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Vector produced by the block route
      real(wp), intent(in) :: got(ndim)
      !> Vector produced by the sparse route
      real(wp), intent(in) :: want(ndim)
      !> Name for the failure message
      character(len=*), intent(in) :: what

      character(len=160) :: message
      integer :: axis

      do axis = 1, ndim
         if (abs(got(axis) - want(axis)) <= &
             max(iswig_block_tol, iswig_block_tol*abs(want(axis)))) cycle
         write (message, "(a,a,i0,a,es24.16,a,es24.16)") what, " mismatch on axis ", &
            axis, ": got ", got(axis), " want ", want(axis)
         call test_failed(error, trim(message))
         return
      end do
   end subroutine check_close

   !* ------------------------------ bordered KKT tests ----------------------------- *!

   !> Deterministic well-conditioned bordered-KKT system, seeded from the case index
   !>
   !> Every case but the last builds `H` symmetric, as the actual `phi_rr - lambda S_rr` is.
   !> The last is asymmetric: the routine under test takes a general `(3, 3)` to test for
   !> indexing errors.
   !>
   !> @param[in]  seed         Case index
   !> @param[out] H_lagrangian Lagrangian Hessian
   !> @param[out] lsf1_r       Level-set gradient
   !> @param[out] rhs          Right-hand side batch
   pure subroutine kkt_fixture(seed, H_lagrangian, lsf1_r, rhs)
      !> Case index
      integer, intent(in) :: seed
      !> Lagrangian Hessian
      real(wp), intent(out) :: H_lagrangian(3, 3)
      !> Level-set gradient
      real(wp), intent(out) :: lsf1_r(3)
      !> Right-hand side batch
      real(wp), intent(out) :: rhs(4, kkt_nrhs)

      integer :: i, j

      do i = 1, ndim
         do j = 1, i
            H_lagrangian(i, j) = sin(0.9_wp*real(7*seed + 5*i + 11*j, wp))
            H_lagrangian(j, i) = H_lagrangian(i, j)
         end do
         H_lagrangian(i, i) = H_lagrangian(i, i) + 3.0_wp
         lsf1_r(i) = 1.0_wp + 0.5_wp*sin(0.7_wp*real(13*seed + 17*i, wp))
      end do
      if (seed == kkt_ncase) then
         do i = 1, ndim
            do j = 1, i - 1
               H_lagrangian(i, j) = H_lagrangian(i, j) &
                                    + 0.5_wp*sin(0.3_wp*real(31*i + 37*j, wp))
            end do
         end do
      end if
      do j = 1, kkt_nrhs
         do i = 1, 4
            rhs(i, j) = sin(0.6_wp*real(3*seed + 19*i + 23*j, wp))
         end do
      end do
   end subroutine kkt_fixture

   !> Assemble the bordered matrix the way the pre-refactor helper did
   !>
   !> @param[in]  H_lagrangian Lagrangian Hessian
   !> @param[in]  lsf1_r       Level-set gradient
   !> @param[out] kkt_mat      Bordered matrix
   pure subroutine kkt_assemble(H_lagrangian, lsf1_r, kkt_mat)
      !> Lagrangian Hessian
      real(wp), intent(in) :: H_lagrangian(3, 3)
      !> Level-set gradient
      real(wp), intent(in) :: lsf1_r(3)
      !> Bordered matrix
      real(wp), intent(out) :: kkt_mat(4, 4)

      kkt_mat = 0.0_wp
      kkt_mat(1:3, 1:3) = H_lagrangian
      kkt_mat(1:3, 4) = -lsf1_r
      kkt_mat(4, 1:3) = lsf1_r
   end subroutine kkt_assemble

   !> The factor object is exactly a bordered assembly plus `getrf`/`getrs`
   subroutine test_kkt_factor_matches_lu(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: kkt
      type(mctc_error), allocatable :: kkt_error
      real(wp) :: H_lagrangian(ndim, ndim), lsf1_r(ndim)
      real(wp) :: rhs_fac(4, kkt_nrhs), rhs_ref(4, kkt_nrhs), kkt_mat(4, 4)
      integer(lapack_ik) :: ipiv(4), info_ref
      integer :: icase, i, j
      character(len=160) :: message

      do icase = 1, kkt_ncase
         call kkt_fixture(icase, H_lagrangian, lsf1_r, rhs_ref)
         rhs_fac = rhs_ref

         call kkt_assemble(H_lagrangian, lsf1_r, kkt_mat)
         call lapack_getrf(4_lapack_ik, 4_lapack_ik, kkt_mat, 4_lapack_ik, ipiv, info_ref)
         if (info_ref /= 0_lapack_ik) then
            write (message, "(a,i0,a,i0)") "KKT case ", icase, &
               ": the reference getrf failed on a nonsingular fixture, status ", info_ref
            call test_failed(error, trim(message))
            return
         end if
         call lapack_getrs("n", 4_lapack_ik, int(kkt_nrhs, lapack_ik), kkt_mat, &
                           4_lapack_ik, ipiv, rhs_ref, 4_lapack_ik, info_ref)

         call kkt%factor(H_lagrangian, lsf1_r, kkt_error)
         if (allocated(kkt_error)) then
            write (message, "(a,i0,2a)") "KKT case ", icase, &
               ": factorization of a nonsingular fixture failed: ", trim(kkt_error%message)
            call test_failed(error, trim(message))
            return
         end if
         call kkt%solve(rhs_fac, kkt_error)
         if (allocated(kkt_error)) then
            write (message, "(a,i0,2a)") "KKT case ", icase, &
               ": solve against valid factors failed: ", trim(kkt_error%message)
            call test_failed(error, trim(message))
            return
         end if

         do j = 1, kkt_nrhs
            do i = 1, 4
               if (rhs_fac(i, j) == rhs_ref(i, j)) cycle
               write (message, "(a,i0,a,i0,1x,i0,a,es24.16,a,es24.16)") "KKT case ", icase, &
                  ": factor solution differs from a direct getrf/getrs at", i, j, &
                  ": ", rhs_fac(i, j), " vs ", rhs_ref(i, j)
               call test_failed(error, trim(message))
               return
            end do
         end do
      end do
   end subroutine test_kkt_factor_matches_lu

   !> The factor/solve pair reproduces the one-shot `gesv` helper it replaced
   subroutine test_kkt_factor_matches_gesv(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: kkt
      type(mctc_error), allocatable :: kkt_error
      real(wp) :: H_lagrangian(ndim, ndim), lsf1_r(ndim)
      real(wp) :: rhs_fac(4, kkt_nrhs), rhs_ref(4, kkt_nrhs), kkt_mat(4, 4)
      real(wp) :: dev
      integer(lapack_ik) :: ipiv(4), info_ref
      integer :: icase, i, j
      character(len=160) :: message

      do icase = 1, kkt_ncase
         call kkt_fixture(icase, H_lagrangian, lsf1_r, rhs_ref)
         rhs_fac = rhs_ref

         call kkt_assemble(H_lagrangian, lsf1_r, kkt_mat)
         call lapack_gesv(4_lapack_ik, int(kkt_nrhs, lapack_ik), kkt_mat, &
                          4_lapack_ik, ipiv, rhs_ref, 4_lapack_ik, info_ref)

         call kkt%factor(H_lagrangian, lsf1_r, kkt_error)
         if (allocated(kkt_error)) then
            write (message, "(a,i0,2a)") "KKT case ", icase, &
               ": factorization of a nonsingular fixture failed: ", trim(kkt_error%message)
            call test_failed(error, trim(message))
            return
         end if
         call kkt%solve(rhs_fac, kkt_error)
         if (allocated(kkt_error)) then
            write (message, "(a,i0,2a)") "KKT case ", icase, &
               ": solve against valid factors failed: ", trim(kkt_error%message)
            call test_failed(error, trim(message))
            return
         end if
         if (info_ref /= 0_lapack_ik) then
            write (message, "(a,i0,a,i0)") "KKT case ", icase, &
               ": gesv failed on a nonsingular fixture, status ", info_ref
            call test_failed(error, trim(message))
            return
         end if
         do j = 1, kkt_nrhs
            do i = 1, 4
               dev = rel_deviation(rhs_fac(i, j), rhs_ref(i, j))
               if (dev <= kkt_gesv_tol) cycle
               write (message, "(a,i0,a,i0,1x,i0,a,es12.4,a,es24.16,a,es24.16)") &
                  "KKT case ", icase, ": factor solution differs from gesv at", i, j, &
                  " by ", dev, ": ", rhs_fac(i, j), " vs ", rhs_ref(i, j)
               call test_failed(error, trim(message))
               return
            end do
         end do
      end do
   end subroutine test_kkt_factor_matches_gesv

   !> Two solves against one factorization equal a single solve of the batch
   subroutine test_kkt_factor_reuse(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: kkt
      type(mctc_error), allocatable :: kkt_error
      real(wp) :: H_lagrangian(ndim, ndim), lsf1_r(ndim)
      real(wp) :: rhs_split(4, kkt_nrhs), rhs_ref(4, kkt_nrhs)
      integer :: icase, i, j
      character(len=160) :: message

      do icase = 1, kkt_ncase
         call kkt_fixture(icase, H_lagrangian, lsf1_r, rhs_ref)
         rhs_split = rhs_ref

         call kkt%factor(H_lagrangian, lsf1_r, kkt_error)
         if (allocated(kkt_error)) then
            write (message, "(a,i0,2a)") "KKT case ", icase, &
               ": factorization of a nonsingular fixture failed: ", trim(kkt_error%message)
            call test_failed(error, trim(message))
            return
         end if

         ! One factorization, one batch -- the reference.
         call kkt%solve(rhs_ref, kkt_error)
         ! The same factorization, reused across two later calls.
         call kkt%solve(rhs_split(:, 1:4), kkt_error)
         call kkt%solve(rhs_split(:, 5:kkt_nrhs), kkt_error)
         if (allocated(kkt_error)) then
            write (message, "(a,i0,2a)") "KKT case ", icase, &
               ": a solve against valid factors failed: ", trim(kkt_error%message)
            call test_failed(error, trim(message))
            return
         end if

         do j = 1, kkt_nrhs
            do i = 1, 4
               if (rhs_split(i, j) == rhs_ref(i, j)) cycle
               write (message, "(a,i0,a,i0,1x,i0,a,es24.16,a,es24.16)") "KKT case ", icase, &
                  ": reused factors differ from a single batch solve at", i, j, &
                  ": ", rhs_split(i, j), " vs ", rhs_ref(i, j)
               call test_failed(error, trim(message))
               return
            end do
         end do
      end do
   end subroutine test_kkt_factor_reuse

   !> A vanishing level-set gradient is reported, not silently solved
   subroutine test_kkt_factor_singular(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(drop_kkt_factor_type) :: kkt
      type(mctc_error), allocatable :: kkt_error
      real(wp) :: H_lagrangian(ndim, ndim), lsf1_r(ndim)
      real(wp) :: rhs_ref(4, kkt_nrhs), kkt_mat(4, 4)
      integer(lapack_ik) :: ipiv(4), info_ref

      call kkt_fixture(1, H_lagrangian, lsf1_r, rhs_ref)
      lsf1_r = 0.0_wp

      call kkt_assemble(H_lagrangian, lsf1_r, kkt_mat)
      call lapack_gesv(4_lapack_ik, int(kkt_nrhs, lapack_ik), kkt_mat, &
                       4_lapack_ik, ipiv, rhs_ref, 4_lapack_ik, info_ref)

      call kkt%factor(H_lagrangian, lsf1_r, kkt_error)

      if (.not. allocated(kkt_error)) then
         call test_failed(error, &
                          "KKT factorization reported success on a vanishing level-set gradient")
         return
      end if
      if (info_ref == 0_lapack_ik) then
         call test_failed(error, &
                          "gesv reported success on a vanishing level-set gradient")
      end if
   end subroutine test_kkt_factor_singular

   !* -------------------- switched-eigenvalue identity tests ------------------------ *!

   !> One draw of the identity fixture's Lehmer generator, mapped to `[-1, 1)`
   !>
   !> The generator state is an argument rather than a `save`d local because
   !> test-drive runs the tests of a suite in an `!$omp parallel do`.
   !>
   !> @param[inout] rng  Generator state
   !> @param[out]   val  Pseudo-random number in `[-1, 1)`
   pure subroutine identity_draw(rng, val)
      !> Generator state
      integer(int64), intent(inout) :: rng
      !> Pseudo-random number
      real(wp), intent(out) :: val

      rng = mod(seed_identity_mult*rng, seed_identity_mod)
      val = 2.0_wp*(real(rng, wp)/real(seed_identity_mod, wp)) - 1.0_wp
   end subroutine identity_draw

   !> `identity_draw` over a vector
   !>
   !> @param[inout] rng  Generator state
   !> @param[out]   vec  Pseudo-random vector
   pure subroutine identity_draw_vec(rng, vec)
      !> Generator state
      integer(int64), intent(inout) :: rng
      !> Pseudo-random vector
      real(wp), intent(out) :: vec(:)

      integer :: i

      do i = 1, size(vec)
         call identity_draw(rng, vec(i))
      end do
   end subroutine identity_draw_vec

   !> A pseudo-random symmetric matrix, as `M + transpose(M)`
   !>
   !> @param[inout] rng  Generator state
   !> @param[out]   mat  Pseudo-random symmetric matrix
   pure subroutine identity_draw_sym(rng, mat)
      !> Generator state
      integer(int64), intent(inout) :: rng
      !> Pseudo-random symmetric matrix
      real(wp), intent(out) :: mat(ndim, ndim)

      integer :: i, j
      real(wp) :: raw(ndim, ndim)

      do j = 1, ndim
         do i = 1, ndim
            call identity_draw(rng, raw(i, j))
         end do
      end do
      mat = raw + transpose(raw)
   end subroutine identity_draw_sym

   !> Build one draw of the switched-eigenvalue identity fixture
   !>
   !> `n_surf` is a random unit vector and `u_switch` is a random vector
   !> projected onto the plane orthogonal to it and normalised, which is the
   !> `n . u = 0` precondition the contraction rests on -- the same one
   !> `setup_tangent_frame` delivers in production.
   !>
   !> `A_mat` and all three of its perturbations are symmetrised, because
   !> `A = alpha I - lambda H` and every derivative of it is a combination of
   !> `H` and its own derivatives.
   !>
   !> The normal perturbations are deliberately left free. The identity never
   !> assumes them orthogonal to `n`, and in production they are not: `dn` is
   !> orthogonal to `n` only to first order, and `ddn` is not orthogonal at all.
   !> Constraining them here would hide any term that survives only off that
   !> constraint.
   !>
   !> @param[inout] rng       Generator state
   !> @param[out]   n_surf    Unit normal
   !> @param[out]   u_switch  Unit vector orthogonal to `n_surf`
   !> @param[out]   A_mat     Symmetric matrix
   !> @param[out]   dA_b      Symmetric `b` perturbation of `A_mat`
   !> @param[out]   dA_v      Symmetric `v` perturbation of `A_mat`
   !> @param[out]   ddA       Symmetric second-order perturbation of `A_mat`
   !> @param[out]   dn_b      Unconstrained `b` perturbation of the normal
   !> @param[out]   dn_v      Unconstrained `v` perturbation of the normal
   !> @param[out]   ddn       Unconstrained second-order perturbation of the normal
   pure subroutine seed_identity_draw(rng, n_surf, u_switch, A_mat, dA_b, dA_v, ddA, &
                                      dn_b, dn_v, ddn)
      !> Generator state
      integer(int64), intent(inout) :: rng
      !> Unit normal and the unit vector orthogonal to it
      real(wp), intent(out) :: n_surf(ndim), u_switch(ndim)
      !> Symmetric matrix and its three symmetric perturbations
      real(wp), intent(out) :: A_mat(ndim, ndim), dA_b(ndim, ndim)
      real(wp), intent(out) :: dA_v(ndim, ndim), ddA(ndim, ndim)
      !> Unconstrained perturbations of the normal
      real(wp), intent(out) :: dn_b(ndim), dn_v(ndim), ddn(ndim)

      real(wp) :: raw(ndim)

      ! A unit normal, and a unit tangent built by projecting the plane out of a
      ! second draw. Both draws are re-taken until they are far from degenerate,
      ! so the normalisations stay well conditioned
      do
         call identity_draw_vec(rng, raw)
         if (norm2(raw) > 0.3_wp) exit
      end do
      n_surf = raw/norm2(raw)
      do
         call identity_draw_vec(rng, raw)
         raw = raw - n_surf*dot_product(n_surf, raw)
         if (norm2(raw) > 0.3_wp) exit
      end do
      u_switch = raw/norm2(raw)

      call identity_draw_sym(rng, A_mat)
      call identity_draw_sym(rng, dA_b)
      call identity_draw_sym(rng, dA_v)
      call identity_draw_sym(rng, ddA)

      call identity_draw_vec(rng, dn_b)
      call identity_draw_vec(rng, dn_v)
      call identity_draw_vec(rng, ddn)
   end subroutine seed_identity_draw

   !> Assemble `dM_b` and `ddM` the obvious way, with `matmul` and outer products
   !>
   !> This is the naive route the contraction replaces, written out exactly as
   !> `apply_seed` and `apply_seed_tangent` used to write it:
   !>
   !>     P    = I - n n^T
   !>     dP_x = -(dn_x n^T + n dn_x^T)
   !>     ddP  = -(ddn n^T + dn_b dn_v^T + dn_v dn_b^T + n ddn^T)
   !>     dM_b = dP_b A P + P dA_b P + P A dP_b
   !>     ddM  = d_v(dM_b), all nine product-rule terms
   !>
   !> @param[in]  n_surf  Unit normal
   !> @param[in]  A_mat   Symmetric matrix
   !> @param[in]  dA_b    Symmetric `b` perturbation of `A_mat`
   !> @param[in]  dA_v    Symmetric `v` perturbation of `A_mat`
   !> @param[in]  ddA     Symmetric second-order perturbation of `A_mat`
   !> @param[in]  dn_b    `b` perturbation of the normal
   !> @param[in]  dn_v    `v` perturbation of the normal
   !> @param[in]  ddn     Second-order perturbation of the normal
   !> @param[out] dM_b    `d_b(P A P)`
   !> @param[out] ddM     `d_v(d_b(P A P))`
   pure subroutine seed_identity_assemble(n_surf, A_mat, dA_b, dA_v, ddA, &
                                          dn_b, dn_v, ddn, dM_b, ddM)
      !> Unit normal
      real(wp), intent(in) :: n_surf(ndim)
      !> Symmetric matrix and its three symmetric perturbations
      real(wp), intent(in) :: A_mat(ndim, ndim), dA_b(ndim, ndim)
      real(wp), intent(in) :: dA_v(ndim, ndim), ddA(ndim, ndim)
      !> Perturbations of the normal
      real(wp), intent(in) :: dn_b(ndim), dn_v(ndim), ddn(ndim)
      !> Naive `d_b(P A P)` and `d_v(d_b(P A P))`
      real(wp), intent(out) :: dM_b(ndim, ndim), ddM(ndim, ndim)

      real(wp) :: P(ndim, ndim), dP_b(ndim, ndim), dP_v(ndim, ndim), ddP(ndim, ndim)
      integer :: i

      P = -outer(n_surf, n_surf)
      do i = 1, ndim
         P(i, i) = P(i, i) + 1.0_wp
      end do
      dP_b = -(outer(dn_b, n_surf) + outer(n_surf, dn_b))
      dP_v = -(outer(dn_v, n_surf) + outer(n_surf, dn_v))
      ddP = -(outer(ddn, n_surf) + outer(dn_b, dn_v) &
              + outer(dn_v, dn_b) + outer(n_surf, ddn))

      dM_b = matmul(dP_b, matmul(A_mat, P)) &
             + matmul(P, matmul(dA_b, P)) &
             + matmul(P, matmul(A_mat, dP_b))

      ddM = matmul(ddP, matmul(A_mat, P)) &
            + matmul(dP_b, matmul(dA_v, P)) &
            + matmul(dP_b, matmul(A_mat, dP_v)) &
            + matmul(dP_v, matmul(dA_b, P)) &
            + matmul(P, matmul(ddA, P)) &
            + matmul(P, matmul(dA_b, dP_v)) &
            + matmul(dP_v, matmul(A_mat, dP_b)) &
            + matmul(P, matmul(dA_v, dP_b)) &
            + matmul(P, matmul(A_mat, ddP))
   end subroutine seed_identity_assemble

   !> Outer product of two 3-vectors
   !>
   !> @param[in] a  Left vector
   !> @param[in] b  Right vector
   pure function outer(a, b) result(mat)
      !> Left and right vectors
      real(wp), intent(in) :: a(ndim), b(ndim)
      !> Outer product `a b^T`
      real(wp) :: mat(ndim, ndim)

      integer :: i, j

      do j = 1, ndim
         do i = 1, ndim
            mat(i, j) = a(i)*b(j)
         end do
      end do
   end function outer

   !> Compare one contracted quantity against its naive assembly
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Draw description, carrying the draw index
   !> @param[in]    name   Quantity name, so a failure localises
   !> @param[in]    got    Contracted value
   !> @param[in]    want   Naive-assembly reference
   subroutine check_seed_identity(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Draw description
      character(len=*), intent(in) :: tag
      !> Quantity name
      character(len=*), intent(in) :: name
      !> Contracted value and its naive-assembly reference
      real(wp), intent(in) :: got, want

      real(wp) :: dev
      character(len=256) :: message

      if (allocated(error)) return
      dev = rel_deviation(got, want)
      if (dev <= seed_identity_tol) return
      write (message, "(5a,es12.4,a,es24.16,a,es24.16)") "seed identity [", tag, "] ", &
         name, " differs from its naive assembly by ", dev, ": ", got, " vs ", want
      call test_failed(error, trim(message))
   end subroutine check_seed_identity

   !> `check_seed_identity` over a vector, with the component in the message
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Draw description
   !> @param[in]    name   Quantity name
   !> @param[in]    got    Contracted value
   !> @param[in]    want   Naive-assembly reference
   subroutine check_seed_identity_vec(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Draw description
      character(len=*), intent(in) :: tag
      !> Quantity name
      character(len=*), intent(in) :: name
      !> Contracted value and its naive-assembly reference
      real(wp), intent(in) :: got(:), want(:)

      character(len=64) :: elem
      integer :: i

      do i = 1, size(got)
         write (elem, "(a,i0,a)") name//"(", i, ")"
         call check_seed_identity(error, tag, trim(elem), got(i), want(i))
         if (allocated(error)) return
      end do
   end subroutine check_seed_identity_vec

   !> `switched_eigenvalue_response` against the naive `P A P` assembly
   !>
   !> No finite differences are involved: both sides are exact functions of the
   !> same draw, so the only admissible gap is roundoff. This is what makes the
   !> contraction a separate routine -- it pins the algebra itself rather than
   !> its downstream consequences.
   !>
   !> @param[out] error  Error handle
   subroutine test_switched_eigenvalue_response_identity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: n_surf(ndim), u_switch(ndim), A_mat(ndim, ndim)
      real(wp) :: dA_b(ndim, ndim), dA_v(ndim, ndim), ddA(ndim, ndim)
      real(wp) :: dn_b(ndim), dn_v(ndim), ddn(ndim)
      real(wp) :: dM_b(ndim, ndim), ddM(ndim, ndim)
      real(wp) :: dlambda_switch, dM_u(ndim), want_vec(ndim), want_scalar
      integer(int64) :: rng
      integer :: idraw
      character(len=64) :: tag

      rng = seed_identity_seed
      do idraw = 1, seed_identity_ndraw
         write (tag, "(a,i0)") "response, draw ", idraw
         call seed_identity_draw(rng, n_surf, u_switch, A_mat, dA_b, dA_v, ddA, &
                                 dn_b, dn_v, ddn)
         call seed_identity_assemble(n_surf, A_mat, dA_b, dA_v, ddA, dn_b, dn_v, ddn, &
                                     dM_b, ddM)

         want_vec = matmul(dM_b, u_switch)
         want_scalar = dot_product(u_switch, want_vec)
         call check_seed_floor(error, trim(tag), "reference", "u . (dM u)", &
                               abs(want_scalar), seed_identity_live_min)
         call check_seed_floor(error, trim(tag), "reference", "dM u", &
                               maxval(abs(want_vec)), seed_identity_live_min)
         if (allocated(error)) return

         call switched_eigenvalue_response(n_surf, u_switch, A_mat, dA_b, dn_b, &
                                           dlambda_switch, dM_u)
         call check_seed_identity(error, trim(tag), "dlambda_switch", dlambda_switch, &
                                  want_scalar)
         call check_seed_identity_vec(error, trim(tag), "dM_u", dM_u, want_vec)
         if (allocated(error)) return

         ! The optional output must not change the scalar it is computed beside
         call switched_eigenvalue_response(n_surf, u_switch, A_mat, dA_b, dn_b, &
                                           dlambda_switch)
         call check_seed_identity(error, trim(tag), "dlambda_switch, no dM_u", &
                                  dlambda_switch, want_scalar)
         if (allocated(error)) return
      end do
   end subroutine test_switched_eigenvalue_response_identity

   !> `switched_eigenvalue_curvature` against the naive `P A P` assembly
   !>
   !> Same construction as the first-order identity test, one order up: the
   !> reference is `u . (ddM u)` with `ddM` assembled from all nine
   !> product-rule terms.
   !>
   !> @param[out] error  Error handle
   subroutine test_switched_eigenvalue_curvature_identity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: n_surf(ndim), u_switch(ndim), A_mat(ndim, ndim)
      real(wp) :: dA_b(ndim, ndim), dA_v(ndim, ndim), ddA(ndim, ndim)
      real(wp) :: dn_b(ndim), dn_v(ndim), ddn(ndim)
      real(wp) :: dM_b(ndim, ndim), ddM(ndim, ndim)
      real(wp) :: ddlambda_switch, want_scalar
      integer(int64) :: rng
      integer :: idraw
      character(len=64) :: tag

      rng = seed_identity_seed
      do idraw = 1, seed_identity_ndraw
         write (tag, "(a,i0)") "curvature, draw ", idraw
         call seed_identity_draw(rng, n_surf, u_switch, A_mat, dA_b, dA_v, ddA, &
                                 dn_b, dn_v, ddn)
         call seed_identity_assemble(n_surf, A_mat, dA_b, dA_v, ddA, dn_b, dn_v, ddn, &
                                     dM_b, ddM)

         want_scalar = dot_product(u_switch, matmul(ddM, u_switch))
         call check_seed_floor(error, trim(tag), "reference", "u . (ddM u)", &
                               abs(want_scalar), seed_identity_live_min)
         if (allocated(error)) return

         call switched_eigenvalue_curvature(n_surf, u_switch, A_mat, dA_b, dA_v, ddA, &
                                            dn_b, dn_v, ddn, ddlambda_switch)
         call check_seed_identity(error, trim(tag), "ddlambda_switch", ddlambda_switch, &
                                  want_scalar)
         if (allocated(error)) return
      end do
   end subroutine test_switched_eigenvalue_curvature_identity

   !* --------------------------- seed-state tangent tests --------------------------- *!

   !> Switching functions of the seed-state tangent fixture
   !>
   !> @param[out] f_crit  Critical-gradient switch
   !> @param[out] f_foc   Focusing switch
   !> @param[out] f_wleb  Lebedev-weight pruning switch
   subroutine seed_state_switches(f_crit, f_foc, f_wleb)
      !> Critical-gradient switch
      type(moist_cavity_drop_swif_sigmoid_bump_type), intent(out) :: f_crit
      !> Focusing switch
      type(moist_cavity_drop_swif_sigmoid_bump_type), intent(out) :: f_foc
      !> Lebedev-weight pruning switch
      type(moist_cavity_drop_swif_sigmoid_bump_type), intent(out) :: f_wleb

      ! The upper edge sits just above `g_norm = 0.9566`, so `f_crit0` lands near
      ! 0.5 rather than on the lower plateau. It multiplies the entire
      ! switched-eigenvalue channel of `res%dw_f`; with the wider [0.2, 2.5] the
      ! prefactor is 8.9e-4 and that channel is three orders below the
      ! gradient-norm channel beside it, far enough down that a wrong term in
      ! `u . (ddM u)` moves `dw_f` by less than the comparison tolerance.
      call new_swif_sigmoid_bump(f_crit, 0.2_wp, 1.7_wp)
      call new_swif_sigmoid_bump(f_foc, -3.0_wp, 3.0_wp)
      ! Narrow on purpose. The fixture sits at `|w_pre_i| ~ 0.21`, and a wide
      ! window such as [0.02, 1.5] puts that on the switch's upper plateau,
      ! where `wleb_prune_factor` and its tangent are both identically zero and
      ! every pruning assertion below becomes vacuous.
      call new_swif_sigmoid_bump(f_wleb, 0.05_wp, 0.40_wp)
   end subroutine seed_state_switches

   !> Fill the `Inputs` block of a seed state with the tangent-test fixture
   !>
   !> The numbers are chosen so that every guarded branch stays on its regular
   !> side: `|B12|` is far from `eig_2x2_symmetric`'s diagonal branch, the
   !> eigenvalue gap `sqrt_disc_B` is far from the degenerate one, the curvature
   !> discriminant is far above `seed_curv_disc_guard`, and `|w_pre_i|` lands
   !> inside the `f_wleb` transition window.
   !>
   !> `lsf3_rrr` is fully symmetric on purpose. `apply_seed` forms `dB12` with
   !> the `A`-symmetry shortcut `q1 . (A dq2) = (A q1) . dq2`, which needs
   !> `res%dH` -- and therefore `lsf3_rrr(:, :, k)` -- symmetric in its first two
   !> indices. A third derivative without that symmetry breaks the whole `B`
   !> chain and everything downstream of it for a reason that is a property of
   !> the fixture, not a defect of the tangent.
   !>
   !> @param[in]  want_curvature  Whether the curvature block is requested
   !> @param[out] state           Seed state with only the `Inputs` block filled
   pure subroutine seed_state_fixture(want_curvature, state)
      !> Whether the curvature block is requested
      logical, intent(in) :: want_curvature
      !> Seed state
      type(drop_seed_state_type), intent(out) :: state

      integer :: i, j, k
      real(wp) :: ri, rj, rk

      state%lsf1_r = [0.70_wp, -0.35_wp, 0.55_wp]
      state%lsf2_rr = reshape([0.90_wp, 0.21_wp, -0.13_wp, &
                               0.21_wp, 0.63_wp, 0.17_wp, &
                               -0.13_wp, 0.17_wp, 1.11_wp], [ndim, ndim])
      do k = 1, ndim
         rk = real(k, wp)
         do j = 1, ndim
            rj = real(j, wp)
            do i = 1, ndim
               ri = real(i, wp)
               state%lsf3_rrr(i, j, k) = 0.013_wp*(ri*rj + rj*rk + ri*rk) &
                                         + 0.007_wp*(ri + rj + rk) &
                                         + 0.004_wp*ri*rj*rk
            end do
         end do
      end do
      state%lambda_val = 0.37_wp
      state%alpha_coeff = 1.23_wp
      state%anchor = [1.30_wp, 0.40_wp, -0.90_wp]
      state%owner_xyz = [0.10_wp, -0.20_wp, 0.30_wp]
      state%anchor_wleb0 = 0.41_wp
      state%cpjac_scal0 = 0.83_wp
      state%w_f0 = 0.61_wp
      state%wbranch = 0.77_wp
      state%wleb = 0.29_wp
      state%xi0 = 1.90_wp
      state%want_curvature = want_curvature
   end subroutine seed_state_fixture

   !> Seed direction of the tangent test; `dlsf2_rr` is symmetric, as `d(grad^2 S)` is
   !>
   !> @param[out] dlsf1_r   Seed perturbation of `grad S` at fixed `r`
   !> @param[out] dlsf2_rr  Seed perturbation of `grad^2 S` at fixed `r`
   !> @param[out] dr        Induced motion of the projected point
   !> @param[out] dlambda   Induced change of the Lagrange multiplier
   pure subroutine seed_state_seed(dlsf1_r, dlsf2_rr, dr, dlambda)
      !> Seed perturbation of the level-set gradient
      real(wp), intent(out) :: dlsf1_r(ndim)
      !> Seed perturbation of the level-set Hessian
      real(wp), intent(out) :: dlsf2_rr(ndim, ndim)
      !> Induced motion of the projected point
      real(wp), intent(out) :: dr(ndim)
      !> Induced change of the multiplier
      real(wp), intent(out) :: dlambda

      dlsf1_r = [0.13_wp, 0.29_wp, -0.19_wp]
      dlsf2_rr = reshape([0.15_wp, 0.09_wp, 0.27_wp, &
                          0.09_wp, -0.22_wp, -0.18_wp, &
                          0.27_wp, -0.18_wp, 0.31_wp], [ndim, ndim])
      dr = [0.23_wp, -0.11_wp, 0.31_wp]
      dlambda = 0.19_wp
   end subroutine seed_state_seed

   !> Displace the fixture `Inputs` along the input tangent the seed induces
   !>
   !> `build_seed_state` maps the `Inputs` block to the `Derived` block, so its
   !> tangent is a directional derivative in those inputs alone. Five of them
   !> move: the level-set gradient and Hessian by `res%dg` and `res%dH`, the
   !> multiplier by `dlambda`, and -- because `cpjac_scal0` and `w_f0` are the
   !> primal `J` and `f` at this very point, which is exactly why `apply_seed`
   !> forms `dw_pre` as the product rule of those two -- by `res%dJ` and
   !> `res%dw_f`. Everything else is held fixed. `lsf3_rrr` is not displaced
   !> because `build_seed_state` never reads it.
   !>
   !> @param[in]  want_curvature  Whether the curvature block is requested
   !> @param[in]  res             Linear response at the base point
   !> @param[in]  dlambda         Induced change of the multiplier
   !> @param[in]  step            Signed displacement
   !> @param[out] state           Displaced seed state, `Inputs` block only
   pure subroutine displace_seed_state(want_curvature, res, dlambda, step, state)
      !> Whether the curvature block is requested
      logical, intent(in) :: want_curvature
      !> Linear response at the base point
      type(drop_seed_result_type), intent(in) :: res
      !> Induced change of the multiplier
      real(wp), intent(in) :: dlambda
      !> Signed displacement
      real(wp), intent(in) :: step
      !> Displaced seed state
      type(drop_seed_state_type), intent(out) :: state

      call seed_state_fixture(want_curvature, state)
      state%lsf1_r = state%lsf1_r + step*res%dg
      state%lsf2_rr = state%lsf2_rr + step*res%dH
      state%lambda_val = state%lambda_val + step*dlambda
      state%cpjac_scal0 = state%cpjac_scal0 + step*res%dJ
      state%w_f0 = state%w_f0 + step*res%dw_f
   end subroutine displace_seed_state

   !> Central difference of the whole derived block, field by field
   !>
   !> @param[in]  state_p   Derived block at `+h`
   !> @param[in]  state_m   Derived block at `-h`
   !> @param[in]  inv2h     `1 / (2 h)`
   !> @param[in]  u_switch  Switched eigenvector at the *base* point
   !> @param[out] fd        Finite-difference reference for the state tangent
   pure subroutine seed_state_central_difference(state_p, state_m, inv2h, u_switch, fd)
      !> Displaced states
      type(drop_seed_state_type), intent(in) :: state_p, state_m
      !> Reciprocal of twice the step
      real(wp), intent(in) :: inv2h
      !> Switched eigenvector at the base point, held fixed on purpose
      real(wp), intent(in) :: u_switch(ndim)
      !> Finite-difference reference
      type(drop_seed_state_tangent_type), intent(out) :: fd

      fd%dA_mat = (state_p%A_mat - state_m%A_mat)*inv2h
      fd%dq1 = (state_p%q1 - state_m%q1)*inv2h
      fd%dq2 = (state_p%q2 - state_m%q2)*inv2h
      fd%dAq1 = (state_p%Aq1 - state_m%Aq1)*inv2h
      fd%dAq2 = (state_p%Aq2 - state_m%Aq2)*inv2h
      fd%dB11 = (state_p%B11 - state_m%B11)*inv2h
      fd%dB12 = (state_p%B12 - state_m%B12)*inv2h
      fd%dB22 = (state_p%B22 - state_m%B22)*inv2h
      fd%ddet_B = (state_p%det_B - state_m%det_B)*inv2h
      fd%dBinv11 = (state_p%Binv11 - state_m%Binv11)*inv2h
      fd%dBinv12 = (state_p%Binv12 - state_m%Binv12)*inv2h
      fd%dBinv22 = (state_p%Binv22 - state_m%Binv22)*inv2h
      fd%dlambda_switch = (state_p%lambda_switch - state_m%lambda_switch)*inv2h
      ! `M = P A P` is not itself a derived field, so the reference assembles the
      ! tangent projector `P = I - n n^T` at either end from the stored normal.
      ! `dM_u` is `(dM) u` at the base eigenvector, not `d(M u)`, so the
      ! difference is taken on `M` alone and contracted with the fixed
      ! `u_switch` afterwards. Differencing the whole product would add `M du`,
      ! which is `O(1)` here because `du` leaves the `u` eigendirection
      block
         real(wp) :: P_p(ndim, ndim), P_m(ndim, ndim)
         integer :: i, j

         do j = 1, ndim
            do i = 1, ndim
               P_p(i, j) = -state_p%n_surf(i)*state_p%n_surf(j)
               P_m(i, j) = -state_m%n_surf(i)*state_m%n_surf(j)
            end do
            P_p(j, j) = P_p(j, j) + 1.0_wp
            P_m(j, j) = P_m(j, j) + 1.0_wp
         end do
         fd%dM_u = matmul((matmul(P_p, matmul(state_p%A_mat, P_p)) &
                           - matmul(P_m, matmul(state_m%A_mat, P_m)))*inv2h, u_switch)
      end block
      fd%dvmin_B = (state_p%vmin_B - state_m%vmin_B)*inv2h
      fd%du_switch = (state_p%u_switch - state_m%u_switch)*inv2h
      fd%dtau1 = (state_p%tau1 - state_m%tau1)*inv2h
      fd%dtau2 = (state_p%tau2 - state_m%tau2)*inv2h
      fd%dw1 = (state_p%w1 - state_m%w1)*inv2h
      fd%dw2 = (state_p%w2 - state_m%w2)*inv2h
      fd%dy1 = (state_p%y1 - state_m%y1)*inv2h
      fd%dy2 = (state_p%y2 - state_m%y2)*inv2h
      fd%dcross_vec = (state_p%cross_vec - state_m%cross_vec)*inv2h
      fd%dinv_J = (state_p%inv_J - state_m%inv_J)*inv2h
      fd%dn_dot_q1 = (state_p%n_dot_q1 - state_m%n_dot_q1)*inv2h
      fd%dproj_surf = (state_p%proj_surf - state_m%proj_surf)*inv2h
      fd%dv_norm_surf = (state_p%v_norm_surf - state_m%v_norm_surf)*inv2h
      fd%df_crit0 = (state_p%f_crit0 - state_m%f_crit0)*inv2h
      fd%df_crit_dS = (state_p%f_crit_dS - state_m%f_crit_dS)*inv2h
      fd%df_foc_f0 = (state_p%f_foc_f0 - state_m%f_foc_f0)*inv2h
      fd%df_foc_dS = (state_p%f_foc_dS - state_m%f_foc_dS)*inv2h
      fd%dwleb_prune_factor = (state_p%wleb_prune_factor - state_m%wleb_prune_factor)*inv2h
      fd%dg_norm_sq = (state_p%g_norm_sq - state_m%g_norm_sq)*inv2h
      fd%dHq1 = (state_p%Hq1 - state_m%Hq1)*inv2h
      fd%dHq2 = (state_p%Hq2 - state_m%Hq2)*inv2h
      fd%dS11 = (state_p%S11 - state_m%S11)*inv2h
      fd%dS12 = (state_p%S12 - state_m%S12)*inv2h
      fd%dS22 = (state_p%S22 - state_m%S22)*inv2h
      fd%dT_curv = (state_p%T_curv - state_m%T_curv)*inv2h
      fd%dKM_curv = (state_p%KM_curv - state_m%KM_curv)*inv2h
      fd%dhalf_diff = (state_p%half_diff - state_m%half_diff)*inv2h
      fd%ddisc_curv = (state_p%disc_curv - state_m%disc_curv)*inv2h
   end subroutine seed_state_central_difference

   !> Fail unless a scalar stays above a floor, naming the quantity
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Fixture and step description
   !> @param[in]    what   What the floor protects
   !> @param[in]    name   Quantity name
   !> @param[in]    val    Quantity value
   !> @param[in]    floor  Smallest admissible value
   subroutine check_seed_floor(error, tag, what, name, val, floor)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> What the floor protects
      character(len=*), intent(in) :: what
      !> Quantity name
      character(len=*), intent(in) :: name
      !> Quantity value
      real(wp), intent(in) :: val
      !> Smallest admissible value
      real(wp), intent(in) :: floor

      character(len=256) :: message

      if (allocated(error)) return
      if (val > floor) return
      write (message, "(7a,es24.16,a,es12.4)") "seed tangent [", tag, "] ", what, " ", name, &
         " collapsed to ", val, ", which must stay above ", floor
      call test_failed(error, trim(message))
   end subroutine check_seed_floor

   !> Fail unless the fixture is still the one the comparison assumes
   !>
   !> A flipped branch makes the finite difference meaningless rather than the
   !> tangent wrong, so these are checked at the base point and at both
   !> displacements before any field is compared.
   !>
   !> @param[inout] error           Error handle
   !> @param[in]    tag             Fixture and step description
   !> @param[in]    state           State to guard
   !> @param[in]    ref             Base-point state, for the discrete choices
   !> @param[in]    status          Status of `build_seed_state`
   !> @param[in]    use_wleb_prune  Whether Lebedev-weight pruning is active
   subroutine check_seed_guards(error, tag, state, ref, status, use_wleb_prune)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> State to guard and the base-point state
      type(drop_seed_state_type), intent(in) :: state, ref
      !> Status of `build_seed_state`
      integer, intent(in) :: status
      !> Whether Lebedev-weight pruning is active
      logical, intent(in) :: use_wleb_prune

      character(len=256) :: message

      if (allocated(error)) return
      if (status /= seed_state_ok) then
         write (message, "(3a,i0)") "seed tangent [", tag, &
            "] fixture is degenerate, build_seed_state returned status ", status
         call test_failed(error, trim(message))
         return
      end if
      if (state%min_axis /= ref%min_axis) then
         write (message, "(3a,i0,a,i0)") "seed tangent [", tag, &
            "] fixture flipped min_axis from ", ref%min_axis, " to ", state%min_axis
         call test_failed(error, trim(message))
         return
      end if
      ! Off `eig_2x2_symmetric`'s `|b| <= 1e-14` diagonal branch, and away from
      ! the degenerate eigenvalue that `dvmin_B` divides by
      call check_seed_floor(error, tag, "guard", "abs(B12)", abs(state%B12), seed_fd_b12_min)
      call check_seed_floor(error, tag, "guard", "sqrt_disc_B", state%sqrt_disc_B, &
                            seed_fd_gap_min)
      call check_seed_floor(error, tag, "guard", "proj_surf", state%proj_surf, &
                            seed_fd_proj_min)
      if (allocated(error)) return
      if (state%want_curvature) then
         call check_seed_floor(error, tag, "guard", "disc_curv", state%disc_curv, &
                               max(seed_fd_disc_min, seed_curv_disc_guard))
         if (allocated(error)) return
      end if
      if (use_wleb_prune) then
         if (sign(1.0_wp, state%w_pre_i) /= sign(1.0_wp, ref%w_pre_i)) then
            write (message, "(3a,es24.16,a,es24.16)") "seed tangent [", tag, &
               "] fixture flipped the sign of w_pre_i from ", ref%w_pre_i, " to ", state%w_pre_i
            call test_failed(error, trim(message))
            return
         end if
         ! The pruning chain is only live off the switch's plateaus
         call check_seed_floor(error, tag, "pruning channel", "abs(f_wleb_ds)", &
                               abs(state%f_wleb_ds), seed_fd_live_min)
         call check_seed_floor(error, tag, "pruning channel", "abs(f_wleb_d2S)", &
                               abs(state%f_wleb_d2S), seed_fd_live_min)
         call check_seed_floor(error, tag, "pruning channel", "abs(wleb_prune_factor)", &
                               abs(state%wleb_prune_factor), seed_fd_live_min)
         if (allocated(error)) return
      end if
   end subroutine check_seed_guards

   !> Compare one tangent component against its central difference
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Fixture and step description
   !> @param[in]    name   Field name, so a failure localises to one field
   !> @param[in]    got    Analytic tangent
   !> @param[in]    want   Central-difference reference
   subroutine check_seed_fd(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Field name
      character(len=*), intent(in) :: name
      !> Analytic tangent and its central-difference reference
      real(wp), intent(in) :: got, want

      real(wp) :: dev
      character(len=256) :: message

      if (allocated(error)) return
      dev = rel_deviation(got, want)
      if (dev <= seed_fd_tol) return
      write (message, "(5a,es12.4,a,es24.16,a,es24.16)") "seed tangent [", tag, "] ", name, &
         " differs from its central difference by ", dev, ": ", got, " vs ", want
      call test_failed(error, trim(message))
   end subroutine check_seed_fd

   !> `check_seed_fd` over a vector field, with the component in the message
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Fixture and step description
   !> @param[in]    name   Field name
   !> @param[in]    got    Analytic tangent
   !> @param[in]    want   Central-difference reference
   subroutine check_seed_fd_vec(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Field name
      character(len=*), intent(in) :: name
      !> Analytic tangent and its central-difference reference
      real(wp), intent(in) :: got(:), want(:)

      character(len=64) :: elem
      integer :: i

      do i = 1, size(got)
         write (elem, "(a,i0,a)") name//"(", i, ")"
         call check_seed_fd(error, tag, trim(elem), got(i), want(i))
         if (allocated(error)) return
      end do
   end subroutine check_seed_fd_vec

   !> `check_seed_fd` over a matrix field, with the element in the message
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Fixture and step description
   !> @param[in]    name   Field name
   !> @param[in]    got    Analytic tangent
   !> @param[in]    want   Central-difference reference
   subroutine check_seed_fd_mat(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Field name
      character(len=*), intent(in) :: name
      !> Analytic tangent and its central-difference reference
      real(wp), intent(in) :: got(:, :), want(:, :)

      character(len=64) :: elem
      integer :: i, j

      do j = 1, size(got, 2)
         do i = 1, size(got, 1)
            write (elem, "(a,i0,a,i0,a)") name//"(", i, ",", j, ")"
            call check_seed_fd(error, tag, trim(elem), got(i, j), want(i, j))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_seed_fd_mat

   !> Fail unless the channels under test actually carry a signal
   !>
   !> A tangent field that is zero for a trivial reason passes any tolerance, so
   !> the finite-difference reference of every field that is expected to move is
   !> required to be non-negligible before it is compared. Fields that are
   !> switched off in this configuration are excluded here and are still
   !> compared -- against a reference of zero, which is the assertion that they
   !> stay off.
   !>
   !> @param[inout] error           Error handle
   !> @param[in]    tag             Fixture and step description
   !> @param[in]    fd              Central-difference reference
   !> @param[in]    want_curvature  Whether the curvature block is requested
   !> @param[in]    use_wleb_prune  Whether Lebedev-weight pruning is active
   subroutine check_seed_tangent_live(error, tag, fd, want_curvature, use_wleb_prune)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Central-difference reference
      type(drop_seed_state_tangent_type), intent(in) :: fd
      !> Configuration flags
      logical, intent(in) :: want_curvature, use_wleb_prune

      if (allocated(error)) return
      call check_seed_floor(error, tag, "reference", "dA_mat", maxval(abs(fd%dA_mat)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dq1", maxval(abs(fd%dq1)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dq2", maxval(abs(fd%dq2)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dAq1", maxval(abs(fd%dAq1)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dAq2", maxval(abs(fd%dAq2)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dB11", abs(fd%dB11), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dB12", abs(fd%dB12), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dB22", abs(fd%dB22), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "ddet_B", abs(fd%ddet_B), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dBinv11", abs(fd%dBinv11), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dBinv12", abs(fd%dBinv12), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dBinv22", abs(fd%dBinv22), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dlambda_switch", &
                            abs(fd%dlambda_switch), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dM_u", maxval(abs(fd%dM_u)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dvmin_B", maxval(abs(fd%dvmin_B)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "du_switch", maxval(abs(fd%du_switch)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dtau1", maxval(abs(fd%dtau1)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dtau2", maxval(abs(fd%dtau2)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dw1", maxval(abs(fd%dw1)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dw2", maxval(abs(fd%dw2)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dy1", maxval(abs(fd%dy1)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dy2", maxval(abs(fd%dy2)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dcross_vec", &
                            maxval(abs(fd%dcross_vec)), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dinv_J", abs(fd%dinv_J), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dn_dot_q1", abs(fd%dn_dot_q1), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dproj_surf", abs(fd%dproj_surf), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dv_norm_surf", abs(fd%dv_norm_surf), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "df_crit0", abs(fd%df_crit0), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "df_crit_dS", abs(fd%df_crit_dS), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "df_foc_f0", abs(fd%df_foc_f0), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "df_foc_dS", abs(fd%df_foc_dS), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dg_norm_sq", abs(fd%dg_norm_sq), &
                            seed_fd_live_min)
      if (allocated(error)) return
      if (use_wleb_prune) then
         call check_seed_floor(error, tag, "reference", "dwleb_prune_factor", &
                               abs(fd%dwleb_prune_factor), seed_fd_live_min)
         if (allocated(error)) return
      end if
      if (want_curvature) then
         call check_seed_floor(error, tag, "reference", "dHq1", maxval(abs(fd%dHq1)), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dHq2", maxval(abs(fd%dHq2)), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dS11", abs(fd%dS11), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dS12", abs(fd%dS12), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dS22", abs(fd%dS22), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dT_curv", abs(fd%dT_curv), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dKM_curv", abs(fd%dKM_curv), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dhalf_diff", abs(fd%dhalf_diff), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "ddisc_curv", abs(fd%ddisc_curv), &
                               seed_fd_live_min)
         if (allocated(error)) return
      end if
   end subroutine check_seed_tangent_live

   !> Compare every component of `drop_seed_state_tangent_type` against the
   !> central difference of the corresponding derived field
   !>
   !> @param[inout] error   Error handle
   !> @param[in]    tag     Fixture and step description
   !> @param[in]    dstate  Analytic state tangent from `apply_seed`
   !> @param[in]    fd      Central-difference reference
   subroutine check_seed_tangent_fields(error, tag, dstate, fd)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Analytic state tangent and its central-difference reference
      type(drop_seed_state_tangent_type), intent(in) :: dstate, fd

      call check_seed_fd_mat(error, tag, "dA_mat", dstate%dA_mat, fd%dA_mat)
      call check_seed_fd_vec(error, tag, "dq1", dstate%dq1, fd%dq1)
      call check_seed_fd_vec(error, tag, "dq2", dstate%dq2, fd%dq2)
      ! The scratch inside `apply_seed` carries `dA . q` only; the state tangent
      ! must add `A . dq`, and this is the assertion that sees the difference
      call check_seed_fd_vec(error, tag, "dAq1", dstate%dAq1, fd%dAq1)
      call check_seed_fd_vec(error, tag, "dAq2", dstate%dAq2, fd%dAq2)
      call check_seed_fd(error, tag, "dB11", dstate%dB11, fd%dB11)
      call check_seed_fd(error, tag, "dB12", dstate%dB12, fd%dB12)
      call check_seed_fd(error, tag, "dB22", dstate%dB22, fd%dB22)
      call check_seed_fd(error, tag, "ddet_B", dstate%ddet_B, fd%ddet_B)
      call check_seed_fd(error, tag, "dBinv11", dstate%dBinv11, fd%dBinv11)
      call check_seed_fd(error, tag, "dBinv12", dstate%dBinv12, fd%dBinv12)
      call check_seed_fd(error, tag, "dBinv22", dstate%dBinv22, fd%dBinv22)
      ! Two independent routes to one scalar. `dstate%dlambda_switch` is the
      ! basis-invariant `u . (dM u)`, while the reference differences the
      ! primal `lambda_switch`, which `eig_2x2_symmetric` builds from
      ! `trace*trace - 4*det`. The cancellation floor of that construction
      ! scales as `eps/sqrt_disc_B**2`, which at the fixture's gap is orders
      ! below the finite-difference noise, so the ordinary tolerance covers it.
      call check_seed_fd(error, tag, "dlambda_switch", dstate%dlambda_switch, &
                         fd%dlambda_switch)
      ! The only direct assertion on the stored `dM_u`; the second-order path
      ! otherwise reaches it only through its contraction with `du_switch`
      call check_seed_fd_vec(error, tag, "dM_u", dstate%dM_u, fd%dM_u)
      call check_seed_fd_vec(error, tag, "dvmin_B", dstate%dvmin_B, fd%dvmin_B)
      call check_seed_fd_vec(error, tag, "du_switch", dstate%du_switch, fd%du_switch)
      call check_seed_fd_vec(error, tag, "dtau1", dstate%dtau1, fd%dtau1)
      call check_seed_fd_vec(error, tag, "dtau2", dstate%dtau2, fd%dtau2)
      call check_seed_fd_vec(error, tag, "dw1", dstate%dw1, fd%dw1)
      call check_seed_fd_vec(error, tag, "dw2", dstate%dw2, fd%dw2)
      call check_seed_fd_vec(error, tag, "dy1", dstate%dy1, fd%dy1)
      call check_seed_fd_vec(error, tag, "dy2", dstate%dy2, fd%dy2)
      call check_seed_fd_vec(error, tag, "dcross_vec", dstate%dcross_vec, fd%dcross_vec)
      call check_seed_fd(error, tag, "dinv_J", dstate%dinv_J, fd%dinv_J)
      call check_seed_fd(error, tag, "dn_dot_q1", dstate%dn_dot_q1, fd%dn_dot_q1)
      call check_seed_fd(error, tag, "dproj_surf", dstate%dproj_surf, fd%dproj_surf)
      call check_seed_fd(error, tag, "dv_norm_surf", dstate%dv_norm_surf, fd%dv_norm_surf)
      call check_seed_fd(error, tag, "df_crit0", dstate%df_crit0, fd%df_crit0)
      call check_seed_fd(error, tag, "df_crit_dS", dstate%df_crit_dS, fd%df_crit_dS)
      call check_seed_fd(error, tag, "df_foc_f0", dstate%df_foc_f0, fd%df_foc_f0)
      call check_seed_fd(error, tag, "df_foc_dS", dstate%df_foc_dS, fd%df_foc_dS)
      ! Zero on both sides when pruning is off, which is the assertion that the
      ! chain stays switched off
      call check_seed_fd(error, tag, "dwleb_prune_factor", dstate%dwleb_prune_factor, &
                         fd%dwleb_prune_factor)
      call check_seed_fd(error, tag, "dg_norm_sq", dstate%dg_norm_sq, fd%dg_norm_sq)
      ! Likewise zero on both sides without `want_curvature`, which is what the
      ! early return in `apply_seed` and the type's default initialisation claim
      call check_seed_fd_vec(error, tag, "dHq1", dstate%dHq1, fd%dHq1)
      call check_seed_fd_vec(error, tag, "dHq2", dstate%dHq2, fd%dHq2)
      call check_seed_fd(error, tag, "dS11", dstate%dS11, fd%dS11)
      call check_seed_fd(error, tag, "dS12", dstate%dS12, fd%dS12)
      call check_seed_fd(error, tag, "dS22", dstate%dS22, fd%dS22)
      call check_seed_fd(error, tag, "dT_curv", dstate%dT_curv, fd%dT_curv)
      call check_seed_fd(error, tag, "dKM_curv", dstate%dKM_curv, fd%dKM_curv)
      call check_seed_fd(error, tag, "dhalf_diff", dstate%dhalf_diff, fd%dhalf_diff)
      call check_seed_fd(error, tag, "ddisc_curv", dstate%ddisc_curv, fd%ddisc_curv)
   end subroutine check_seed_tangent_fields

   !> Central-difference the derived seed state along the seed-induced tangent
   !>
   !> `build_seed_state` is a map from the `Inputs` block to the `Derived`
   !> block, and `apply_seed`'s optional `dstate` claims to be its directional
   !> derivative along the input tangent the seed induces. This displaces the
   !> inputs by `+/- h` along that tangent, rebuilds the state at each end and
   !> compares `(F(+h) - F(-h))/(2h)` against every component of
   !> `drop_seed_state_tangent_type`, plus the two `drop_seed_result_type`
   !> channels that share the same displaced states.
   !>
   !> @param[out] error           Error handle
   !> @param[in]  want_curvature  Whether the curvature block is requested
   !> @param[in]  use_wleb_prune  Whether Lebedev-weight pruning is active
   subroutine run_seed_state_tangent(error, want_curvature, use_wleb_prune)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Configuration flags
      logical, intent(in) :: want_curvature, use_wleb_prune

      type(moist_cavity_drop_swif_sigmoid_bump_type) :: f_crit, f_foc, f_wleb
      type(drop_seed_state_type) :: state, state_p, state_m
      type(drop_seed_result_type) :: res
      type(drop_seed_state_tangent_type) :: dstate, fd
      real(wp) :: dlsf1_r(ndim), dlsf2_rr(ndim, ndim), dr(ndim), dlambda
      real(wp) :: h, inv2h, fd_gnorm, fd_jval
      integer :: status, status_p, status_m, istep
      character(len=64) :: tag

      call seed_state_switches(f_crit, f_foc, f_wleb)
      call seed_state_fixture(want_curvature, state)
      call seed_state_seed(dlsf1_r, dlsf2_rr, dr, dlambda)

      call build_seed_state(state, f_crit, f_foc, f_wleb, use_wleb_prune, status)
      call check_seed_guards(error, "base point", state, state, status, use_wleb_prune)
      if (allocated(error)) return

      call apply_seed(state, dlsf1_r, dlsf2_rr, dr, dlambda, res, dstate)

      do istep = 1, size(seed_fd_steps)
         h = seed_fd_steps(istep)
         inv2h = 0.5_wp/h
         write (tag, "(a,l1,a,l1,a,es9.2)") "curv ", want_curvature, ", prune ", &
            use_wleb_prune, ", h ", h

         call displace_seed_state(want_curvature, res, dlambda, h, state_p)
         call build_seed_state(state_p, f_crit, f_foc, f_wleb, use_wleb_prune, status_p)
         call check_seed_guards(error, trim(tag)//", +h", state_p, state, status_p, &
                                use_wleb_prune)
         if (allocated(error)) return

         call displace_seed_state(want_curvature, res, dlambda, -h, state_m)
         call build_seed_state(state_m, f_crit, f_foc, f_wleb, use_wleb_prune, status_m)
         call check_seed_guards(error, trim(tag)//", -h", state_m, state, status_m, &
                                use_wleb_prune)
         if (allocated(error)) return

         call seed_state_central_difference(state_p, state_m, inv2h, state%u_switch, fd)
         call check_seed_tangent_live(error, trim(tag), fd, want_curvature, use_wleb_prune)
         if (allocated(error)) return
         call check_seed_tangent_fields(error, trim(tag), dstate, fd)
         if (allocated(error)) return

         ! Two `drop_seed_result_type` channels the same displaced states pin
         fd_gnorm = (state_p%g_norm - state_m%g_norm)*inv2h
         call check_seed_fd(error, trim(tag), "res%d_gnorm", res%d_gnorm, fd_gnorm)
         fd_jval = (norm2(state_p%cross_vec) - norm2(state_m%cross_vec))*inv2h
         call check_seed_fd(error, trim(tag), "res%dJ", res%dJ, fd_jval)
         if (allocated(error)) return
      end do
   end subroutine run_seed_state_tangent

   !> Seed-state tangent without curvature and without Lebedev-weight pruning
   subroutine test_seed_state_tangent_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_state_tangent(error, .false., .false.)
   end subroutine test_seed_state_tangent_plain

   !> Seed-state tangent with Lebedev-weight pruning, without curvature
   subroutine test_seed_state_tangent_prune(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_state_tangent(error, .false., .true.)
   end subroutine test_seed_state_tangent_prune

   !> Seed-state tangent with curvature, without Lebedev-weight pruning
   subroutine test_seed_state_tangent_curvature(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_state_tangent(error, .true., .false.)
   end subroutine test_seed_state_tangent_curvature

   !> Seed-state tangent with curvature and Lebedev-weight pruning
   subroutine test_seed_state_tangent_curvature_prune(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_state_tangent(error, .true., .true.)
   end subroutine test_seed_state_tangent_curvature_prune


   !* -------------------------- seed-result tangent tests --------------------------- *!

   !> Second-direction tangent of the seed of [[seed_state_seed]]
   !>
   !> `ddlsf2_rr` is symmetric because `d(grad^2 S)` is, and because
   !> `apply_seed`'s `dB12` shortcut needs the seed Hessian symmetric at every
   !> displaced point, not only at the base point.
   !>
   !> @param[out] ddlsf1_r   Second-direction tangent of the seed `grad S`
   !> @param[out] ddlsf2_rr  Second-direction tangent of the seed `grad^2 S`
   !> @param[out] ddr        Second-direction tangent of the point motion
   !> @param[out] ddlambda   Second-direction tangent of the multiplier motion
   pure subroutine seed_result_seed_tangent(ddlsf1_r, ddlsf2_rr, ddr, ddlambda)
      !> Second-direction tangent of the seed gradient
      real(wp), intent(out) :: ddlsf1_r(ndim)
      !> Second-direction tangent of the seed Hessian
      real(wp), intent(out) :: ddlsf2_rr(ndim, ndim)
      !> Second-direction tangent of the point motion
      real(wp), intent(out) :: ddr(ndim)
      !> Second-direction tangent of the multiplier motion
      real(wp), intent(out) :: ddlambda

      ddlsf1_r = [-0.17_wp, 0.24_wp, 0.11_wp]
      ddlsf2_rr = reshape([0.21_wp, -0.14_wp, 0.08_wp, &
                           -0.14_wp, 0.17_wp, 0.25_wp, &
                           0.08_wp, 0.25_wp, -0.19_wp], [ndim, ndim])
      ddr = [-0.09_wp, 0.26_wp, 0.14_wp]
      ddlambda = -0.21_wp
   end subroutine seed_result_seed_tangent

   !> `Inputs` tangent of the second direction `v`
   !>
   !> Three of the ten fields are not free. `apply_seed` is the only producer
   !> of the `dstate_v`/`res_v` pair that `apply_seed_tangent` consumes, and it
   !> carries no seed channel for `anchor_wleb0`, and it builds
   !> `dwleb_prune_factor` out of `res%dJ` and `res%dw_f`. So `danchor_wleb0` is
   !> zero whenever pruning is active and free otherwise, because with pruning
   !> off `anchor_wleb0` never reaches the derived block at all; and
   !> `dcpjac_scal0`/`dw_f0` are left for the caller to fill once `res_v` exists.
   !>
   !> The jet fields are free but must carry the symmetry the primal assumes:
   !> `dlsf2_rr` in both indices and `dlsf3_rrr` in its first two, for exactly
   !> the reason `seed_state_fixture` documents.
   !>
   !> @param[in]  use_wleb_prune  Whether Lebedev-weight pruning is active
   !> @param[out] dinp_v          Input tangent along the second direction
   pure subroutine seed_result_input_tangent(use_wleb_prune, dinp_v)
      !> Whether Lebedev-weight pruning is active
      logical, intent(in) :: use_wleb_prune
      !> Input tangent along the second direction
      type(drop_seed_input_tangent_type), intent(out) :: dinp_v

      integer :: i, j, k
      real(wp) :: ri, rj, rk

      dinp_v%dlsf1_r = [0.31_wp, 0.17_wp, -0.24_wp]
      dinp_v%dlsf2_rr = reshape([0.19_wp, -0.11_wp, 0.23_wp, &
                                 -0.11_wp, 0.26_wp, 0.14_wp, &
                                 0.23_wp, 0.14_wp, -0.18_wp], [ndim, ndim])
      do k = 1, ndim
         rk = real(k, wp)
         do j = 1, ndim
            rj = real(j, wp)
            do i = 1, ndim
               ri = real(i, wp)
               dinp_v%dlsf3_rrr(i, j, k) = 0.009_wp*(ri*rj + rj*rk + ri*rk) &
                                           - 0.005_wp*(ri + rj + rk) &
                                           + 0.003_wp*ri*rj*rk
            end do
         end do
      end do
      dinp_v%dlambda_val = -0.27_wp
      if (use_wleb_prune) then
         dinp_v%danchor_wleb0 = 0.0_wp
      else
         dinp_v%danchor_wleb0 = 0.33_wp
      end if
      dinp_v%dcpjac_scal0 = 0.0_wp
      dinp_v%dw_f0 = 0.0_wp
      dinp_v%dwbranch = -0.13_wp
      dinp_v%dwleb = 0.21_wp
      dinp_v%dxi0 = -0.35_wp
   end subroutine seed_result_input_tangent

   !> Displace the `Inputs` block of the fixture linearly along `dinp_v`
   !>
   !> Unlike `displace_seed_state`, which moves only the five inputs the
   !> first-order state tangent depends on, this moves every input `dinp_v`
   !> carries: the result tangent reads `lsf3_rrr`, `wbranch`, `wleb` and `xi0`
   !> as well. `alpha_coeff` is not among them, because it is a fixed parameter
   !> with no tangent field; see `drop_seed_input_tangent_type`.
   !>
   !> @param[in]  want_curvature  Whether the curvature block is requested
   !> @param[in]  dinp_v          Input tangent along the second direction
   !> @param[in]  step            Signed displacement
   !> @param[out] state           Displaced seed state, `Inputs` block only
   pure subroutine displace_seed_inputs(want_curvature, dinp_v, step, state)
      !> Whether the curvature block is requested
      logical, intent(in) :: want_curvature
      !> Input tangent along the second direction
      type(drop_seed_input_tangent_type), intent(in) :: dinp_v
      !> Signed displacement
      real(wp), intent(in) :: step
      !> Displaced seed state
      type(drop_seed_state_type), intent(out) :: state

      call seed_state_fixture(want_curvature, state)
      state%lsf1_r = state%lsf1_r + step*dinp_v%dlsf1_r
      state%lsf2_rr = state%lsf2_rr + step*dinp_v%dlsf2_rr
      state%lsf3_rrr = state%lsf3_rrr + step*dinp_v%dlsf3_rrr
      state%lambda_val = state%lambda_val + step*dinp_v%dlambda_val
      state%anchor_wleb0 = state%anchor_wleb0 + step*dinp_v%danchor_wleb0
      state%cpjac_scal0 = state%cpjac_scal0 + step*dinp_v%dcpjac_scal0
      state%w_f0 = state%w_f0 + step*dinp_v%dw_f0
      state%wbranch = state%wbranch + step*dinp_v%dwbranch
      state%wleb = state%wleb + step*dinp_v%dwleb
      state%xi0 = state%xi0 + step*dinp_v%dxi0
   end subroutine displace_seed_inputs

   !> Central difference of the whole linear response, field by field
   !>
   !> @param[in]  res_p  Response at `+h`
   !> @param[in]  res_m  Response at `-h`
   !> @param[in]  inv2h  `1 / (2 h)`
   !> @param[out] fd     Finite-difference reference for the result tangent
   pure subroutine seed_result_central_difference(res_p, res_m, inv2h, fd)
      !> Displaced responses
      type(drop_seed_result_type), intent(in) :: res_p, res_m
      !> Reciprocal of twice the step
      real(wp), intent(in) :: inv2h
      !> Finite-difference reference
      type(drop_seed_result_tangent_type), intent(out) :: fd

      fd%dg = (res_p%dg - res_m%dg)*inv2h
      fd%dH = (res_p%dH - res_m%dH)*inv2h
      fd%dn_surf = (res_p%dn_surf - res_m%dn_surf)*inv2h
      fd%d_gnorm = (res_p%d_gnorm - res_m%d_gnorm)*inv2h
      fd%dJ = (res_p%dJ - res_m%dJ)*inv2h
      fd%dw_f = (res_p%dw_f - res_m%dw_f)*inv2h
      fd%dwleb = (res_p%dwleb - res_m%dwleb)*inv2h
      fd%dxi = (res_p%dxi - res_m%dxi)*inv2h
      fd%dk1 = (res_p%dk1 - res_m%dk1)*inv2h
      fd%dk2 = (res_p%dk2 - res_m%dk2)*inv2h
   end subroutine seed_result_central_difference

   !> Fail unless the result channels under test actually carry a signal
   !>
   !> Same anti-vacuity rule as `check_seed_tangent_live`: a component that is
   !> zero for a trivial reason passes any tolerance. `dk1`/`dk2` are excluded
   !> without `want_curvature` and are still compared, against a reference of
   !> zero, which is the assertion that the early return leaves them off.
   !>
   !> @param[inout] error           Error handle
   !> @param[in]    tag             Fixture and step description
   !> @param[in]    fd              Central-difference reference
   !> @param[in]    want_curvature  Whether the curvature block is requested
   subroutine check_seed_result_live(error, tag, fd, want_curvature)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Central-difference reference
      type(drop_seed_result_tangent_type), intent(in) :: fd
      !> Whether the curvature block is requested
      logical, intent(in) :: want_curvature

      if (allocated(error)) return
      call check_seed_floor(error, tag, "reference", "dg", maxval(abs(fd%dg)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dH", maxval(abs(fd%dH)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dn_surf", maxval(abs(fd%dn_surf)), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "d_gnorm", abs(fd%d_gnorm), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dJ", abs(fd%dJ), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dw_f", abs(fd%dw_f), seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dwleb", abs(fd%dwleb), &
                            seed_fd_live_min)
      call check_seed_floor(error, tag, "reference", "dxi", abs(fd%dxi), seed_fd_live_min)
      if (allocated(error)) return
      if (want_curvature) then
         call check_seed_floor(error, tag, "reference", "dk1", abs(fd%dk1), &
                               seed_fd_live_min)
         call check_seed_floor(error, tag, "reference", "dk2", abs(fd%dk2), &
                               seed_fd_live_min)
         if (allocated(error)) return
      end if
   end subroutine check_seed_result_live

   !> Compare every component of `drop_seed_result_tangent_type` against the
   !> central difference of the corresponding response channel
   !>
   !> @param[inout] error  Error handle
   !> @param[in]    tag    Fixture and step description
   !> @param[in]    dres   Analytic result tangent from `apply_seed_tangent`
   !> @param[in]    fd     Central-difference reference
   subroutine check_seed_result_fields(error, tag, dres, fd)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture and step description
      character(len=*), intent(in) :: tag
      !> Analytic result tangent and its central-difference reference
      type(drop_seed_result_tangent_type), intent(in) :: dres, fd

      call check_seed_fd_vec(error, tag, "dg", dres%dg, fd%dg)
      call check_seed_fd_mat(error, tag, "dH", dres%dH, fd%dH)
      call check_seed_fd_vec(error, tag, "dn_surf", dres%dn_surf, fd%dn_surf)
      call check_seed_fd(error, tag, "d_gnorm", dres%d_gnorm, fd%d_gnorm)
      call check_seed_fd(error, tag, "dJ", dres%dJ, fd%dJ)
      call check_seed_fd(error, tag, "dw_f", dres%dw_f, fd%dw_f)
      call check_seed_fd(error, tag, "dwleb", dres%dwleb, fd%dwleb)
      call check_seed_fd(error, tag, "dxi", dres%dxi, fd%dxi)
      ! Zero on both sides without `want_curvature`, which is what the early
      ! return in `apply_seed_tangent` and the type's default initialisation claim
      call check_seed_fd(error, tag, "dk1", dres%dk1, fd%dk1)
      call check_seed_fd(error, tag, "dk2", dres%dk2, fd%dk2)
   end subroutine check_seed_result_fields

   !> Central-difference the linear response along a second direction
   !>
   !> `apply_seed_tangent` claims to be `d/dv [ res_b(p, b) ]`, so this displaces
   !> both arguments of that map linearly -- every `Inputs` field by `dinp_v` and
   !> all four seed components by their own tangent -- rebuilds the state at each
   !> end, calls `apply_seed` exactly there, and compares
   !> `(res_b(+h) - res_b(-h))/(2h)` against every component of
   !> `drop_seed_result_tangent_type`. This is a first central difference, not a
   !> nested one: `res_b(t)` is exact at every `t`.
   !>
   !> @param[out] error           Error handle
   !> @param[in]  want_curvature  Whether the curvature block is requested
   !> @param[in]  use_wleb_prune  Whether Lebedev-weight pruning is active
   subroutine run_seed_result_tangent(error, want_curvature, use_wleb_prune)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Configuration flags
      logical, intent(in) :: want_curvature, use_wleb_prune

      type(moist_cavity_drop_swif_sigmoid_bump_type) :: f_crit, f_foc, f_wleb
      type(drop_seed_state_type) :: state, state_p, state_m
      type(drop_seed_result_type) :: res_b, res_v, res_p, res_m
      type(drop_seed_state_tangent_type) :: dstate_b, dstate_v
      type(drop_seed_input_tangent_type) :: dinp_v
      type(drop_seed_result_tangent_type) :: dres, fd
      real(wp) :: dlsf1_r(ndim), dlsf2_rr(ndim, ndim), dr(ndim), dlambda
      real(wp) :: ddlsf1_r(ndim), ddlsf2_rr(ndim, ndim), ddr(ndim), ddlambda
      real(wp) :: dlsf1_r_v(ndim), dlsf2_rr_v(ndim, ndim), dr_v(ndim), dlambda_v
      real(wp) :: h, inv2h
      integer :: status, status_p, status_m, istep, kaxis
      character(len=64) :: tag

      call seed_state_switches(f_crit, f_foc, f_wleb)
      call seed_state_fixture(want_curvature, state)
      call build_seed_state(state, f_crit, f_foc, f_wleb, use_wleb_prune, status)
      call check_seed_guards(error, "base point", state, state, status, use_wleb_prune)
      if (allocated(error)) return

      call seed_state_seed(dlsf1_r, dlsf2_rr, dr, dlambda)
      call seed_result_seed_tangent(ddlsf1_r, ddlsf2_rr, ddr, ddlambda)
      call seed_result_input_tangent(use_wleb_prune, dinp_v)

      ! `dinp_v` and the `v` seed are not independent. `dstate_v` is the tangent
      ! of the derived block along `v`, so the seed has to be the one that
      ! produces the chosen displacement of the inputs: `apply_seed` forms
      ! `res%dg = dlsf1_r + lsf2_rr . dr` and `res%dH = dlsf2_rr + lsf3 . dr`, and
      ! those totals are what `build_seed_state` sees move. Subtracting the
      ! point-motion term off the target leaves the jet part of the seed. `dr_v`
      ! is deliberately nonzero, so that coupling is exercised rather than
      ! trivialised away
      dr_v = [-0.17_wp, 0.28_wp, 0.13_wp]
      dlsf1_r_v = dinp_v%dlsf1_r - matmul(state%lsf2_rr, dr_v)
      dlsf2_rr_v = dinp_v%dlsf2_rr
      do kaxis = 1, ndim
         dlsf2_rr_v = dlsf2_rr_v - state%lsf3_rrr(:, :, kaxis)*dr_v(kaxis)
      end do
      dlambda_v = dinp_v%dlambda_val

      call apply_seed(state, dlsf1_r_v, dlsf2_rr_v, dr_v, dlambda_v, res_v, dstate_v)
      ! `cpjac_scal0` and `w_f0` are the primal `J` and `f` at this very point, so
      ! their displacement is fixed by the `v` response rather than free, exactly
      ! as `displace_seed_state` fixes it for the first-order test
      dinp_v%dcpjac_scal0 = res_v%dJ
      dinp_v%dw_f0 = res_v%dw_f

      call apply_seed(state, dlsf1_r, dlsf2_rr, dr, dlambda, res_b, dstate_b)
      call apply_seed_tangent(state, dstate_v, dinp_v, res_v, &
                              dlsf1_r, dlsf2_rr, dr, dlambda, &
                              ddlsf1_r, ddlsf2_rr, ddr, ddlambda, &
                              res_b, dstate_b, dres)

      do istep = 1, size(seed_fd_steps)
         h = seed_fd_steps(istep)
         inv2h = 0.5_wp/h
         write (tag, "(a,l1,a,l1,a,es9.2)") "result curv ", want_curvature, ", prune ", &
            use_wleb_prune, ", h ", h

         call displace_seed_inputs(want_curvature, dinp_v, h, state_p)
         call build_seed_state(state_p, f_crit, f_foc, f_wleb, use_wleb_prune, status_p)
         call check_seed_guards(error, trim(tag)//", +h", state_p, state, status_p, &
                                use_wleb_prune)
         if (allocated(error)) return
         call apply_seed(state_p, dlsf1_r + h*ddlsf1_r, dlsf2_rr + h*ddlsf2_rr, &
                         dr + h*ddr, dlambda + h*ddlambda, res_p)

         call displace_seed_inputs(want_curvature, dinp_v, -h, state_m)
         call build_seed_state(state_m, f_crit, f_foc, f_wleb, use_wleb_prune, status_m)
         call check_seed_guards(error, trim(tag)//", -h", state_m, state, status_m, &
                                use_wleb_prune)
         if (allocated(error)) return
         call apply_seed(state_m, dlsf1_r - h*ddlsf1_r, dlsf2_rr - h*ddlsf2_rr, &
                         dr - h*ddr, dlambda - h*ddlambda, res_m)

         call seed_result_central_difference(res_p, res_m, inv2h, fd)
         call check_seed_result_live(error, trim(tag), fd, want_curvature)
         if (allocated(error)) return
         call check_seed_result_fields(error, trim(tag), dres, fd)
         if (allocated(error)) return
      end do
   end subroutine run_seed_result_tangent

   !> Seed-result tangent without curvature and without Lebedev-weight pruning
   !* ================================================================================= *!
   !*                    Asymmetric jet-basis seeds and their rescue                    *!
   !* ================================================================================= *!
   !
   ! Every other fixture in this file seeds a *symmetric* Hessian, because that is
   ! what `apply_seed`'s `dB12` shortcut and `apply_seed_tangent`'s three
   ! factor-of-two collapses require. Production does not: `seed_jet_basis` emits
   ! the nine Hessian basis seeds as single-entry matrices, so the six
   ! off-diagonal ones are maximally asymmetric and each one's `dres` is wrong by
   ! order 100%, sign included.
   !
   ! The reason the shipped adjoint is nonetheless correct is linearity: the
   ! error cancels exactly between `E_ij` and `E_ji`. That cancellation is load
   ! bearing and, until these tests, entirely undefended -- the precondition
   ! block on `apply_seed_tangent` tells a driver to rely on it, and nothing
   ! checked that it holds.

   !> Sum two result tangents component by component
   pure subroutine seed_result_tangent_add(a, b, c)
      !> Terms
      type(drop_seed_result_tangent_type), intent(in) :: a, b
      !> Sum
      type(drop_seed_result_tangent_type), intent(out) :: c

      c%dg = a%dg + b%dg
      c%dH = a%dH + b%dH
      c%dn_surf = a%dn_surf + b%dn_surf
      c%d_gnorm = a%d_gnorm + b%d_gnorm
      c%dJ = a%dJ + b%dJ
      c%dw_f = a%dw_f + b%dw_f
      c%dwleb = a%dwleb + b%dwleb
      c%dxi = a%dxi + b%dxi
      c%dk1 = a%dk1 + b%dk1
      c%dk2 = a%dk2 + b%dk2
   end subroutine seed_result_tangent_add

   !> Fail unless two result tangents agree to `seed_linearity_tol`
   subroutine check_seed_result_same(error, tag, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture description
      character(len=*), intent(in) :: tag
      !> Quantities that must agree
      type(drop_seed_result_tangent_type), intent(in) :: got, want

      call check_seed_lin_vec(error, tag, "dg", got%dg, want%dg)
      call check_seed_lin_mat(error, tag, "dH", got%dH, want%dH)
      call check_seed_lin_vec(error, tag, "dn_surf", got%dn_surf, want%dn_surf)
      call check_seed_lin(error, tag, "d_gnorm", got%d_gnorm, want%d_gnorm)
      call check_seed_lin(error, tag, "dJ", got%dJ, want%dJ)
      call check_seed_lin(error, tag, "dw_f", got%dw_f, want%dw_f)
      call check_seed_lin(error, tag, "dwleb", got%dwleb, want%dwleb)
      call check_seed_lin(error, tag, "dxi", got%dxi, want%dxi)
      call check_seed_lin(error, tag, "dk1", got%dk1, want%dk1)
      call check_seed_lin(error, tag, "dk2", got%dk2, want%dk2)
   end subroutine check_seed_result_same

   !> Scalar arm of [[check_seed_result_same]]
   subroutine check_seed_lin(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture description and field name
      character(len=*), intent(in) :: tag, name
      !> Quantities that must agree
      real(wp), intent(in) :: got, want

      real(wp) :: dev
      character(len=256) :: message

      if (allocated(error)) return
      dev = rel_deviation(got, want)
      if (dev <= seed_linearity_tol) return
      write (message, "(5a,es12.4,a,es24.16,a,es24.16)") "seed linearity [", tag, "] ", &
         name, " breaks by ", dev, ": ", got, " vs ", want
      call test_failed(error, trim(message))
   end subroutine check_seed_lin

   !> Vector arm of [[check_seed_result_same]]
   subroutine check_seed_lin_vec(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture description and field name
      character(len=*), intent(in) :: tag, name
      !> Quantities that must agree
      real(wp), intent(in) :: got(:), want(:)

      character(len=32) :: label
      integer :: i

      do i = 1, size(got)
         if (allocated(error)) return
         write (label, "(a,a,i0,a)") name, "(", i, ")"
         call check_seed_lin(error, tag, trim(label), got(i), want(i))
      end do
   end subroutine check_seed_lin_vec

   !> Matrix arm of [[check_seed_result_same]]
   subroutine check_seed_lin_mat(error, tag, name, got, want)
      !> Error handle
      type(error_type), allocatable, intent(inout) :: error
      !> Fixture description and field name
      character(len=*), intent(in) :: tag, name
      !> Quantities that must agree
      real(wp), intent(in) :: got(:, :), want(:, :)

      character(len=32) :: label
      integer :: i, j

      do j = 1, size(got, 2)
         do i = 1, size(got, 1)
            if (allocated(error)) return
            write (label, "(a,a,i0,a,i0,a)") name, "(", i, ",", j, ")"
            call check_seed_lin(error, tag, trim(label), got(i, j), want(i, j))
         end do
      end do
   end subroutine check_seed_lin_mat

   !> Check that the asymmetric-seed error cancels over a transpose pair
   !>
   !> Two independent assertions, and they fail for different reasons:
   !>
   !>   1. **Linearity.** `dres(E12) + dres(E21)` must equal `dres(E12 + E21)`
   !>      computed directly, to rounding. This is algebra, not a finite
   !>      difference, so it is held to `seed_linearity_tol`. It would break if
   !>      anyone introduced a seed-dependent branch into the second-order chain.
   !>   2. **Correctness of the pair.** That same sum must equal the central
   !>      difference of `apply_seed` along `v` taken with the *pair* seed. This
   !>      is the statement a driver actually depends on: the per-seed values are
   !>      wrong, and their sum is right.
   !>
   !> The seed tangents are all zero on purpose. `seed_jet_basis` emits constant
   !> basis seeds, so `d_v(seed) = 0` is the production case, not a simplification.
   !>
   !> @param[out] error           Error handle
   !> @param[in]  want_curvature  Whether the curvature block is requested
   subroutine run_seed_result_asym_pair(error, want_curvature)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Whether the curvature block is requested
      logical, intent(in) :: want_curvature

      type(moist_cavity_drop_swif_sigmoid_bump_type) :: f_crit, f_foc, f_wleb
      type(drop_seed_state_type) :: state, state_p, state_m
      type(drop_seed_result_type) :: res_v, res_ij, res_ji, res_sym, res_p, res_m
      type(drop_seed_state_tangent_type) :: dstate_v, dstate_ij, dstate_ji, dstate_sym
      type(drop_seed_input_tangent_type) :: dinp_v
      type(drop_seed_result_tangent_type) :: dres_ij, dres_ji, dres_sym, dres_pair, fd
      real(wp) :: e_ij(ndim, ndim), e_ji(ndim, ndim), e_sym(ndim, ndim)
      real(wp) :: zero_r(ndim), zero_m(ndim, ndim)
      real(wp) :: dlsf1_r_v(ndim), dlsf2_rr_v(ndim, ndim), dr_v(ndim), dlambda_v
      real(wp) :: h, inv2h
      integer :: status, status_p, status_m, istep, kaxis
      character(len=64) :: tag

      zero_r = 0.0_wp
      zero_m = 0.0_wp

      ! Exactly what `seed_jet_basis` emits for the (1,2) and (2,1) basis
      ! elements: a single unit entry, nothing symmetrised
      e_ij = 0.0_wp
      e_ij(1, 2) = 1.0_wp
      e_ji = 0.0_wp
      e_ji(2, 1) = 1.0_wp
      e_sym = e_ij + e_ji

      call seed_state_switches(f_crit, f_foc, f_wleb)
      call seed_state_fixture(want_curvature, state)
      call build_seed_state(state, f_crit, f_foc, f_wleb, .false., status)
      call check_seed_guards(error, "asym base point", state, state, status, .false.)
      if (allocated(error)) return

      ! Same `v` construction as `run_seed_result_tangent`; only the `b` seed
      ! differs here
      call seed_result_input_tangent(.false., dinp_v)
      dr_v = [-0.17_wp, 0.28_wp, 0.13_wp]
      dlsf1_r_v = dinp_v%dlsf1_r - matmul(state%lsf2_rr, dr_v)
      dlsf2_rr_v = dinp_v%dlsf2_rr
      do kaxis = 1, ndim
         dlsf2_rr_v = dlsf2_rr_v - state%lsf3_rrr(:, :, kaxis)*dr_v(kaxis)
      end do
      dlambda_v = dinp_v%dlambda_val

      call apply_seed(state, dlsf1_r_v, dlsf2_rr_v, dr_v, dlambda_v, res_v, dstate_v)
      dinp_v%dcpjac_scal0 = res_v%dJ
      dinp_v%dw_f0 = res_v%dw_f

      call apply_seed(state, zero_r, e_ij, zero_r, 0.0_wp, res_ij, dstate_ij)
      call apply_seed(state, zero_r, e_ji, zero_r, 0.0_wp, res_ji, dstate_ji)
      call apply_seed(state, zero_r, e_sym, zero_r, 0.0_wp, res_sym, dstate_sym)

      call apply_seed_tangent(state, dstate_v, dinp_v, res_v, &
                              zero_r, e_ij, zero_r, 0.0_wp, &
                              zero_r, zero_m, zero_r, 0.0_wp, &
                              res_ij, dstate_ij, dres_ij)
      call apply_seed_tangent(state, dstate_v, dinp_v, res_v, &
                              zero_r, e_ji, zero_r, 0.0_wp, &
                              zero_r, zero_m, zero_r, 0.0_wp, &
                              res_ji, dstate_ji, dres_ji)
      call apply_seed_tangent(state, dstate_v, dinp_v, res_v, &
                              zero_r, e_sym, zero_r, 0.0_wp, &
                              zero_r, zero_m, zero_r, 0.0_wp, &
                              res_sym, dstate_sym, dres_sym)

      call seed_result_tangent_add(dres_ij, dres_ji, dres_pair)

      write (tag, "(a,l1)") "asym linearity, curv ", want_curvature
      call check_seed_result_same(error, trim(tag), dres_pair, dres_sym)
      if (allocated(error)) return

      do istep = 1, size(seed_fd_steps)
         h = seed_fd_steps(istep)
         inv2h = 0.5_wp/h
         write (tag, "(a,l1,a,es9.2)") "asym pair, curv ", want_curvature, ", h ", h

         call displace_seed_inputs(want_curvature, dinp_v, h, state_p)
         call build_seed_state(state_p, f_crit, f_foc, f_wleb, .false., status_p)
         call check_seed_guards(error, trim(tag)//", +h", state_p, state, status_p, .false.)
         if (allocated(error)) return
         call apply_seed(state_p, zero_r, e_sym, zero_r, 0.0_wp, res_p)

         call displace_seed_inputs(want_curvature, dinp_v, -h, state_m)
         call build_seed_state(state_m, f_crit, f_foc, f_wleb, .false., status_m)
         call check_seed_guards(error, trim(tag)//", -h", state_m, state, status_m, .false.)
         if (allocated(error)) return
         call apply_seed(state_m, zero_r, e_sym, zero_r, 0.0_wp, res_m)

         call seed_result_central_difference(res_p, res_m, inv2h, fd)
         call check_seed_result_fields(error, trim(tag), dres_pair, fd)
         if (allocated(error)) return
      end do
   end subroutine run_seed_result_asym_pair

   !> Non-curvature dispatch for the asymmetric-pair rescue check.
   subroutine test_seed_result_asym_pair_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      call run_seed_result_asym_pair(error, .false.)
   end subroutine test_seed_result_asym_pair_plain

   !> Curvature dispatch for the asymmetric-pair rescue check.
   subroutine test_seed_result_asym_pair_curvature(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      call run_seed_result_asym_pair(error, .true.)
   end subroutine test_seed_result_asym_pair_curvature

   subroutine test_seed_result_tangent_plain(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_result_tangent(error, .false., .false.)
   end subroutine test_seed_result_tangent_plain

   !> Seed-result tangent with Lebedev-weight pruning, without curvature
   subroutine test_seed_result_tangent_prune(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_result_tangent(error, .false., .true.)
   end subroutine test_seed_result_tangent_prune

   !> Seed-result tangent with curvature, without Lebedev-weight pruning
   subroutine test_seed_result_tangent_curvature(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_result_tangent(error, .true., .false.)
   end subroutine test_seed_result_tangent_curvature

   !> Seed-result tangent with curvature and Lebedev-weight pruning
   subroutine test_seed_result_tangent_curvature_prune(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_seed_result_tangent(error, .true., .true.)
   end subroutine test_seed_result_tangent_curvature_prune

end module test_cavity_drop_primitives
