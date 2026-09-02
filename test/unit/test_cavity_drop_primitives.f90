module test_cavity_drop_primitives
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
                                          new_swif_smooth_step
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig, new_iswig, &
                                         iswig_workspace_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none (type, external)
   private

   public :: collect_cavity_drop_primitives

   integer, parameter :: ndim = 3
   integer, parameter :: owner_dummy = 1

   real(wp), parameter :: ABS_THR = 5.0e-9_wp
   real(wp), parameter :: REL_THR = 5.0e-9_wp
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

   !> Tolerances of the iSwiG finite-difference comparisons
   real(wp), parameter :: iswig_fd_abs = 1.0e-8_wp
   real(wp), parameter :: iswig_fd_rel = 1.0e-8_wp

   !> Angular offsets from the owner/neighbour sphere-intersection circle
   real(wp), parameter :: iswig_theta_off(n_iswig_dirs) = &
                          [-0.12_wp, -0.04_wp, 0.04_wp, 0.12_wp]

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
                  new_unittest("iswig_swi2_block_guarded", test_iswig_swi2_block_guarded) &
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
         read (unit, '(a1)', iostat=stat) first
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

      write (message, '(a,1x,a,1x,i0,1x,i0,a,es24.16,a,es24.16)') &
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

      write (message, '(a,a,a,es12.4,a,es12.4)') "iSwiG ", what, &
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
                     call check_exact_zero(error, hvp(axis, iatom), "translation hvp")
                     if (allocated(error)) return
                  end do
               end do
               call check_exact_zero(error, dxi2, "translation dxi2")
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
         call check_exact_zero(error, dxi, "dxi")
         if (allocated(error)) return
         call check_exact_zero(error, dxi2, "dxi2")
         if (allocated(error)) return

         call iswig%swi1_rA(pos, 1, iswig_fd_xi, work, grad)
         call iswig_hvp_dense(iswig, pos, 1, iswig_fd_xi, v, 0.5_wp, hvp, dxi2)
         do iatom = 1, 2
            do axis = 1, ndim
               call check_exact_zero(error, grad(axis, iatom), "isolated gradient")
               if (allocated(error)) return
               call check_exact_zero(error, hvp(axis, iatom), "isolated hvp")
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
                           write (message, '(a,4(1x,i0),a,es24.16,a,es24.16)') &
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
                        call check_exact_zero(error, hvp(axis), "block translation row")
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
                  call check_exact_zero(error, macc(axis), "block translation mix")
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
                     call check_exact_zero(error, rows2(axis, jj), &
                                           "sparse translation row")
                     if (allocated(error)) return
                  end do
               end do
               do axis = 1, ndim
                  call check_exact_zero(error, owner_row2(axis), &
                                        "sparse translation owner row")
                  if (allocated(error)) return
               end do
               call check_exact_zero(error, dxi2_sp, "sparse translation dxi2")
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
                     call check_exact_zero(error, blk(i, far, j, ia), "saturated row")
                     if (allocated(error)) return
                     call check_exact_zero(error, blk(i, ia, j, far), "saturated column")
                     if (allocated(error)) return
                  end do
               end do
            end do
            do i = 1, ndim
               call check_exact_zero(error, mix(i, far), "saturated mix")
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
                  call check_exact_zero(error, blk(i, ia, j, ib), what)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
      do ia = 1, n
         do i = 1, ndim
            call check_exact_zero(error, mix(i, ia), what)
            if (allocated(error)) return
         end do
      end do
      call check_exact_zero(error, d2xi, what)
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
         write (message, '(a,a,i0,a,es24.16,a,es24.16)') what, " mismatch on axis ", &
            axis, ": got ", got(axis), " want ", want(axis)
         call test_failed(error, trim(message))
         return
      end do
   end subroutine check_close

   !> Assert a value is exactly zero, not merely small
   !>
   !> @param[out] error Error handle
   !> @param[in]  value Value that must be a true zero
   !> @param[in]  what  Name for the failure message
   subroutine check_exact_zero(error, value, what)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Value that must be a true zero
      real(wp), intent(in) :: value
      !> Name for the failure message
      character(len=*), intent(in) :: what

      character(len=160) :: message

      if (value == 0.0_wp) return

      write (message, '(a,a,a,es24.16)') "iSwiG ", what, " is not exactly zero: ", value
      call test_failed(error, trim(message))
   end subroutine check_exact_zero

end module test_cavity_drop_primitives
