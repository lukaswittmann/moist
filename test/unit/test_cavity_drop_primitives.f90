module test_cavity_drop_primitives
   use mctc_env_accuracy, only: wp
   use mctc_io, only: structure_type
   use test_helpers, only: get_test_structures, get_test_radii, get_test_points, fd4_scalar
   use moist_cavity_drop_lsf_svdw_ssd, only: ssd0, ssd1_r, ssd2_rr, ssd3_rrr, ssd4_rrrr, &
                                    ssd012_r, ssd2_r_rA, ssd1_rA, ssd2_rArB
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_objective_phi, only: moist_cavity_drop_objective_phi_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_switching, only: moist_cavity_drop_smooth_step_swif, &
                                         new_smooth_step_swif
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_cavity_drop_primitives

   public :: ndim, atol, rtol

   integer, parameter :: ndim = 3
   integer, parameter :: owner_dummy = 1

   real(wp), parameter :: atol = 5.0e-9_wp
   real(wp), parameter :: rtol = 5.0e-9_wp
   real(wp), parameter :: fd_h = 1.0e-4_wp

contains

   !> Collect SSD/phi primitive FD tests plus the switching nuclear-grad
   subroutine collect_cavity_drop_primitives(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ssd_f1_r        ", test_ssd_f1_r), &
                  new_unittest("ssd_f2_rr       ", test_ssd_f2_rr), &
                  new_unittest("ssd_f3_rrr      ", test_ssd_f3_rrr), &
                  new_unittest("ssd_f4_rrrr     ", test_ssd_f4_rrrr), &
                  new_unittest("ssd_f1_rA       ", test_ssd_f1_rA), &
                  new_unittest("ssd_f2_rArB     ", test_ssd_f2_rArB), &
                  new_unittest("ssd_f2_r_rA     ", test_ssd_f2_r_rA), &
                  new_unittest("ssd_f012_r      ", test_ssd_f012_r), &
                  new_unittest("phi_f0          ", test_phi_f0), &
                  new_unittest("phi_f1_r        ", test_phi_f1_r), &
                  new_unittest("phi_f2_rr       ", test_phi_f2_rr), &
                  new_unittest("phi_f3_rrr      ", test_phi_f3_rrr), &
                  new_unittest("phi_f4_rrrr     ", test_phi_f4_rrrr), &
                  new_unittest("phi_f1_rA       ", test_phi_f1_rA), &
                  new_unittest("phi_f2_rArB     ", test_phi_f2_rArB), &
                  new_unittest("phi_f2_r_rA     ", test_phi_f2_r_rA), &
                  new_unittest("phi_f012_r      ", test_phi_f012_r), &
                  new_unittest("switching_f1_rA ", test_switching_f1_rA) &
                  ]
   end subroutine collect_cavity_drop_primitives

   subroutine test_ssd_f1_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim)
      real(wp) :: numeric(ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: work_point(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd1_r(point, center, radius)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  f_pp = ssd0(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  f_p = ssd0(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  f_m = ssd0(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  f_mm = ssd0(work_point, center, radius)
                  numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
               end do

               do i = 1, ndim
                  call check(error, analytic(i), numeric(i), &
                                                thr_abs=atol, thr_rel=rtol)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f1_r

   subroutine test_ssd_f2_rr(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, j, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim, ndim)
      real(wp) :: numeric(ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: work_point(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd2_rr(point, center, radius)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  g_pp = ssd1_r(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  g_p = ssd1_r(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  g_m = ssd1_r(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  g_mm = ssd1_r(work_point, center, radius)
                  do i = 1, ndim
                     numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
                  end do
               end do

               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic(i, j), numeric(i, j), &
                                                   thr_abs=atol, thr_rel=rtol)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f2_rr

   subroutine test_ssd_f3_rrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, j, k, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim, ndim, ndim)
      real(wp) :: numeric(ndim, ndim, ndim)
      real(wp) :: hess_pp(ndim, ndim), hess_p(ndim, ndim), hess_m(ndim, ndim), hess_mm(ndim, ndim)
      real(wp) :: work_point(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd3_rrr(point, center, radius)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  hess_pp = ssd2_rr(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  hess_p = ssd2_rr(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  hess_m = ssd2_rr(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  hess_mm = ssd2_rr(work_point, center, radius)
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
                                                      thr_abs=atol, thr_rel=rtol)
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f3_rrr

   subroutine test_ssd_f4_rrrr(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, j, k, m, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim, ndim, ndim, ndim)
      real(wp) :: numeric(ndim, ndim, ndim, ndim)
      real(wp) :: third_pp(ndim, ndim, ndim), third_p(ndim, ndim, ndim)
      real(wp) :: third_m(ndim, ndim, ndim), third_mm(ndim, ndim, ndim)
      real(wp) :: work_point(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd4_rrrr(point, center, radius)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  third_pp = ssd3_rrr(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  third_p = ssd3_rrr(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  third_m = ssd3_rrr(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  third_mm = ssd3_rrr(work_point, center, radius)
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
                              thr_abs=atol, thr_rel=rtol)
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f4_rrrr

   subroutine test_ssd_f2_r_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, j, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim, ndim)
      real(wp) :: numeric(ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: work_center(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd2_r_rA(point, center, radius)

               do axis = 1, ndim
                  work_center = center
                  work_center(axis) = center(axis) + 2.0_wp*h
                  g_pp = ssd1_r(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) + h
                  g_p = ssd1_r(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) - h
                  g_m = ssd1_r(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) - 2.0_wp*h
                  g_mm = ssd1_r(point, work_center, radius)
                  do i = 1, ndim
                     numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
                  end do
               end do

               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic(i, j), numeric(i, j), &
                                                   thr_abs=atol, thr_rel=rtol)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f2_r_rA

   subroutine test_ssd_f1_rA(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim)
      real(wp) :: numeric(ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: work_center(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd1_rA(point, center, radius)

               do axis = 1, ndim
                  work_center = center
                  work_center(axis) = center(axis) + 2.0_wp*h
                  f_pp = ssd0(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) + h
                  f_p = ssd0(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) - h
                  f_m = ssd0(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) - 2.0_wp*h
                  f_mm = ssd0(point, work_center, radius)
                  numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
               end do

               do i = 1, ndim
                  call check(error, analytic(i), numeric(i), &
                                                thr_abs=atol, thr_rel=rtol)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f1_rA

   subroutine test_ssd_f2_rArB(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, j, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic(ndim, ndim)
      real(wp) :: numeric(ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: work_center(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               analytic = ssd2_rArB(point, center, radius)

               do axis = 1, ndim
                  work_center = center
                  work_center(axis) = center(axis) + 2.0_wp*h
                  g_pp = ssd1_rA(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) + h
                  g_p = ssd1_rA(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) - h
                  g_m = ssd1_rA(point, work_center, radius)
                  work_center = center
                  work_center(axis) = center(axis) - 2.0_wp*h
                  g_mm = ssd1_rA(point, work_center, radius)
                  do i = 1, ndim
                     numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
                  end do
               end do

               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic(i, j), numeric(i, j), &
                                                   thr_abs=atol, thr_rel=rtol)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f2_rArB

   !> Test combined SSD value, gradient, and Hessian against finite differences.
   subroutine test_ssd_f012_r(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), allocatable :: mols(:)
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ipt, jc, i, j, axis
      real(wp) :: point(ndim), center(ndim), radius
      real(wp) :: analytic_val, numeric_val
      real(wp) :: analytic_grad(ndim), numeric_grad(ndim)
      real(wp) :: analytic_hess(ndim, ndim), numeric_hess(ndim, ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: work_point(ndim)
      real(wp) :: h

      h = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol = mols(icase)
         call get_test_radii(mol, radii)
         call get_test_points(mol, points)

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            do jc = 1, mol%nat
               center = mol%xyz(:, jc)
               radius = radii(jc)
               call ssd012_r(point, center, radius, analytic_val, analytic_grad, analytic_hess)
               numeric_val = ssd0(point, center, radius)

               do axis = 1, ndim
                  work_point = point
                  work_point(axis) = point(axis) + 2.0_wp*h
                  f_pp = ssd0(work_point, center, radius)
                  g_pp = ssd1_r(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) + h
                  f_p = ssd0(work_point, center, radius)
                  g_p = ssd1_r(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - h
                  f_m = ssd0(work_point, center, radius)
                  g_m = ssd1_r(work_point, center, radius)
                  work_point = point
                  work_point(axis) = point(axis) - 2.0_wp*h
                  f_mm = ssd0(work_point, center, radius)
                  g_mm = ssd1_r(work_point, center, radius)

                  numeric_grad(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
                  do i = 1, ndim
                     numeric_hess(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
                  end do
               end do

               call check(error, analytic_val, numeric_val, thr_abs=atol, thr_rel=rtol)
               if (allocated(error)) return
               do i = 1, ndim
                  call check(error, analytic_grad(i), numeric_grad(i), &
                                                thr_abs=atol, thr_rel=rtol)
                  if (allocated(error)) return
               end do
               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic_hess(i, j), numeric_hess(i, j), &
                                                   thr_abs=atol, thr_rel=rtol)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_ssd_f012_r

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
               call check(error, analytic, reference, thr_abs=atol, thr_rel=rtol)
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
      h = fd_h

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
                  call check(error, analytic(i), numeric(i), thr_abs=atol, thr_rel=rtol)
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
      h = fd_h

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
                     call check(error, analytic(i, j), numeric(i, j), thr_abs=atol, thr_rel=rtol)
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
         -1.24_wp,  0.56_wp,  0.20_wp, &
          0.60_wp, -0.90_wp,  0.80_wp], [ndim, 2])

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = fd_h
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
                                             thr_abs=atol, thr_rel=rtol)
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
         -1.24_wp,  0.56_wp,  0.20_wp, &
          0.60_wp, -0.90_wp,  0.80_wp], [ndim, 2])

      param%phi_alpha = 0.7_wp
      call phi%set_parameters(param)
      h = fd_h
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
                                                thr_abs=atol, thr_rel=rtol)
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
      h = fd_h

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
                  call check(error, analytic(i, j), numeric(i, j), thr_abs=atol, thr_rel=rtol)
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
      h = fd_h

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
                                                      thr_abs=atol, thr_rel=rtol)
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
      h = fd_h

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
                                                thr_abs=atol, thr_rel=rtol)
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
      h = fd_h

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

               call check(error, analytic_val, numeric_val, thr_abs=atol, thr_rel=rtol)
               if (allocated(error)) return
               do i = 1, ndim
                  call check(error, analytic_grad(i), numeric_grad(i), &
                                                thr_abs=atol, thr_rel=rtol)
                  if (allocated(error)) return
               end do
               do i = 1, ndim
                  do j = 1, ndim
                     call check(error, analytic_hess(i, j), numeric_hess(i, j), &
                                                   thr_abs=atol, thr_rel=rtol)
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
      type(moist_cavity_drop_smooth_step_swif) :: sw
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

      eps = fd_h

      call get_test_structures(mols)
      do icase = 1, size(mols)
         mol_base = mols(icase)
         call get_test_radii(mol_base, radii)
         call get_test_points(mol_base, points)
         allocate (centers_base(ndim, mol_base%nat), centers_local(ndim, mol_base%nat))
         centers_base = mol_base%xyz
         call new_smooth_step_swif(sw, -0.5_wp, 0.5_wp)

         prim%screening_threshold = 0.0_wp
         call prim%new()
         call prim%update(mol_base, radii)
         call prim%set_max_deriv(3)

         if (allocated(lsf1)) deallocate (lsf1)
         allocate (lsf1(ndim, mol_base%nat))

         do ipt = 1, size(points, 2)
            point = points(:, ipt)
            call prim%update(mol_base, radii)
            call prim%ssd_system%update(centers_base, radii)
            call prim%prepare(point)
            call prim%f0_screened(lsf0)
            call prim%f3_rr_rA_screened(lsf1_rA=lsf1, &
                                        lsf3_rr_rA=dummy_rr_rA)
            analytic = sw%f1_rA(lsf0, lsf1)

            do atom = 1, mol_base%nat
               do axis = 1, ndim
                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) + 2.0_wp*eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%ssd_system%update(centers_local, radii)
                  call prim%prepare(point)
                  call prim%f0_screened(lsf0_tmp)
                  f_pp = sw%f0(lsf0_tmp)

                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) + eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%ssd_system%update(centers_local, radii)
                  call prim%prepare(point)
                  call prim%f0_screened(lsf0_tmp)
                  f_p = sw%f0(lsf0_tmp)

                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) - eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%ssd_system%update(centers_local, radii)
                  call prim%prepare(point)
                  call prim%f0_screened(lsf0_tmp)
                  f_m = sw%f0(lsf0_tmp)

                  mol_shift = mol_base
                  centers_local = centers_base
                  centers_local(axis, atom) = centers_local(axis, atom) - 2.0_wp*eps
                  mol_shift%xyz = centers_local
                  call prim%update(mol_shift, radii)
                  call prim%ssd_system%update(centers_local, radii)
                  call prim%prepare(point)
                  call prim%f0_screened(lsf0_tmp)
                  f_mm = sw%f0(lsf0_tmp)

                  numeric = fd4_scalar(f_pp, f_p, f_m, f_mm, eps)

                  call check(error, &
                                                analytic(axis, atom), numeric, &
                                                thr_abs=atol, thr_rel=rtol)
                  if (allocated(error)) return
               end do
            end do
         end do
         deallocate (centers_base, centers_local)
      end do
   end subroutine test_switching_f1_rA

end module test_cavity_drop_primitives
