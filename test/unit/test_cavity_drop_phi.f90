module test_cavity_drop_phi
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, fd4_scalar
   use moist_cavity_drop_phi, only: moist_cavity_drop_phi_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_math_packed, only: unpack_sym_matrix
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_cavity_drop_phi

   integer, parameter :: ndim = 3
   integer, parameter :: owner_dummy = 1

   real(wp), parameter :: ABS_THR = 5.0e-9_wp
   real(wp), parameter :: REL_THR = 5.0e-9_wp
   real(wp), parameter :: STEP_SIZE = 1.0e-4_wp

contains

   !> Collect the tests of the projection objective phi
   subroutine collect_cavity_drop_phi(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("phi_f0", test_phi_f0), &
                  new_unittest("phi_f1_r", test_phi_f1_r), &
                  new_unittest("phi_f2_rr", test_phi_f2_rr), &
                  new_unittest("phi_f3_rrr", test_phi_f3_rrr), &
                  new_unittest("phi_f4_rrrr", test_phi_f4_rrrr), &
                  new_unittest("phi_f2_r_ra", test_phi_f2_r_rA), &
                  new_unittest("phi_f012_r", test_phi_f012_r) &
                  ]
   end subroutine collect_cavity_drop_phi

   !> Test phi value against the direct quadratic expression
   subroutine test_phi_f0(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
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

   !> Test phi point gradient against a finite difference of the value
   subroutine test_phi_f1_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
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
                  call fd4_scalar(f_pp, f_p, f_m, f_mm, h, numeric(axis), error)
                  if (allocated(error)) return
               end do

               do i = 1, ndim
                  call check(error, analytic(i), numeric(i), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_phi_f1_r

   !> Test phi point Hessian against a finite difference of the point gradient
   subroutine test_phi_f2_rr(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim), analytic(ndim, ndim), numeric(ndim, ndim)
      real(wp) :: analytic_packed(6)
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
               analytic_packed = phi%f2_rr(point, anchor, owner_dummy)
               call unpack_sym_matrix(analytic_packed, analytic, ndim)

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
                     call fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h, numeric(i, axis), error)
                     if (allocated(error)) return
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

   !> Test phi third point derivative against a finite difference of the Hessian
   !> Mol-less pure-math test: phi = 0.5*alpha*(r-anchor)^2 has analytic
   !> derivatives that hold for any two distinct points, so the (point,
   !> anchor) pair is hard-coded rather than sourced from a molecule
   subroutine test_phi_f3_rrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      real(wp) :: analytic(ndim, ndim, ndim), numeric(ndim, ndim, ndim)
      real(wp) :: hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
      real(wp) :: hess_packed(6)
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
         hess_packed = phi%f2_rr(work_point, anchor, owner_dummy)
         call unpack_sym_matrix(hess_packed, hess_pp, ndim)
         work_point = point
         work_point(axis) = point(axis) + h
         hess_packed = phi%f2_rr(work_point, anchor, owner_dummy)
         call unpack_sym_matrix(hess_packed, hess_p, ndim)
         work_point = point
         work_point(axis) = point(axis) - h
         hess_packed = phi%f2_rr(work_point, anchor, owner_dummy)
         call unpack_sym_matrix(hess_packed, hess_m, ndim)
         work_point = point
         work_point(axis) = point(axis) - 2.0_wp*h
         hess_packed = phi%f2_rr(work_point, anchor, owner_dummy)
         call unpack_sym_matrix(hess_packed, hess_mm, ndim)
         do i = 1, ndim
            do j = 1, ndim
               call fd4_scalar(hess_pp(i, j), hess_p(i, j), hess_m(i, j), hess_mm(i, j), h, numeric(i, j, axis), error)
               if (allocated(error)) return
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

   !> Test phi fourth point derivative against a finite difference of the third derivative
   !> Mol-less pure-math test (see test_phi_f3_rrr)
   subroutine test_phi_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
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
                  call fd4_scalar(third_pp(i, j, k), third_p(i, j, k), third_m(i, j, k), third_mm(i, j, k), h, &
                                  numeric(i, j, k, axis), error)
                  if (allocated(error)) return
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

   !> Test phi mixed point-nuclear Hessian against finite differences of point gradient
   subroutine test_phi_f2_r_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
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
                  call fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h, numeric(i, axis, owner), error)
                  if (allocated(error)) return
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

   !> Test combined phi value, gradient, and Hessian against finite differences
   subroutine test_phi_f012_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cavity_drop_phi_type) :: phi
      type(moist_cavity_drop_parameters_type) :: param
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: point(ndim), anchor(ndim), work_point(ndim)
      real(wp) :: analytic_val, numeric_val
      real(wp) :: analytic_grad(ndim), numeric_grad(ndim)
      real(wp) :: analytic_hess(ndim, ndim), numeric_hess(ndim, ndim)
      real(wp) :: analytic_hess_packed(6)
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
               call phi%f012_r(point, anchor, owner_dummy, analytic_val, analytic_grad, &
                               analytic_hess_packed)
               call unpack_sym_matrix(analytic_hess_packed, analytic_hess, ndim)
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

                  call fd4_scalar(f_pp, f_p, f_m, f_mm, h, numeric_grad(axis), error)
                  if (allocated(error)) return
                  do i = 1, ndim
                     call fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h, numeric_hess(i, axis), error)
                     if (allocated(error)) return
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

end module test_cavity_drop_phi
