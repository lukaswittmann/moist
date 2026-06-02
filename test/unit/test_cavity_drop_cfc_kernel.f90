!> Kernel-level unit tests for the auto-generated CFC pseudo-density module.
!>
!> Validates the sympy-generated `cfc_atomic_term_eval`, `cfc_pair_term_eval`,
!> and `cfc_log_lift` routines independently of the orchestrator. Catches
!> regressions in the codegen script before they propagate into `cfc.f90`.
!>
!> Each derivative tensor is checked by 4-point central finite difference
!> against the next-lower-order tensor. Index symmetries of the symbolic
!> tensors are checked separately.
module test_cavity_drop_cfc_kernel
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_cfc_kernel, only: cfc_atomic_term_eval, &
      cfc_pair_term_eval, cfc_log_lift
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_cavity_drop_cfc_kernel

   integer, parameter :: ndim = 3

   !> Klamt-Diedenhofen 2018 reference parameter set
   real(wp), parameter :: a1_ref = -15.0_wp
   real(wp), parameter :: a2_ref = -9.0_wp
   real(wp), parameter :: c_ref  = 5.0_wp

   !> Auxiliary parameter sweeps (the kernel is parametric in a1, a2, c)
   real(wp), parameter :: a1_values(2) = [-15.0_wp, -8.0_wp]
   real(wp), parameter :: a2_values(2) = [-9.0_wp, -5.0_wp]
   real(wp), parameter :: c_values(2)  = [5.0_wp, 2.5_wp]

   !> Atomic-term test triples (d_a, R_a)
   integer,  parameter :: n_atom_pts = 6
   real(wp), parameter :: atom_pts(ndim, n_atom_pts) = reshape([ &
      0.60_wp,  0.20_wp, -0.10_wp, &
      1.10_wp, -0.40_wp,  0.30_wp, &
     -0.50_wp,  0.90_wp,  0.20_wp, &
      0.80_wp,  0.80_wp,  0.80_wp, &
     -1.30_wp, -0.40_wp,  0.70_wp, &
      0.30_wp, -1.20_wp,  0.50_wp], [ndim, n_atom_pts])
   real(wp), parameter :: atom_radii(n_atom_pts) = [ &
      0.90_wp, 1.10_wp, 1.30_wp, 0.85_wp, 1.50_wp, 1.00_wp]

   !> Pair-term test quadruples (d_a, d_b, R_a, R_b)
   integer,  parameter :: n_pair_pts = 5
   real(wp), parameter :: pair_d_a(ndim, n_pair_pts) = reshape([ &
      0.50_wp,  0.30_wp,  0.10_wp, &
      1.00_wp, -0.20_wp,  0.40_wp, &
     -0.40_wp,  0.80_wp,  0.20_wp, &
      0.70_wp, -0.50_wp, -0.30_wp, &
     -0.90_wp, -0.30_wp,  0.50_wp], [ndim, n_pair_pts])
   real(wp), parameter :: pair_d_b(ndim, n_pair_pts) = reshape([ &
     -0.40_wp,  0.20_wp,  0.30_wp, &
      0.30_wp,  0.90_wp, -0.10_wp, &
      0.60_wp, -0.20_wp,  0.80_wp, &
     -0.60_wp,  0.40_wp,  0.30_wp, &
      0.70_wp, -0.50_wp, -0.20_wp], [ndim, n_pair_pts])
   real(wp), parameter :: pair_R_a(n_pair_pts) = [0.90_wp, 1.10_wp, 1.30_wp, 0.95_wp, 1.20_wp]
   real(wp), parameter :: pair_R_b(n_pair_pts) = [1.20_wp, 0.85_wp, 1.05_wp, 1.40_wp, 0.95_wp]

   !> FD steps and thresholds
   real(wp), parameter :: STEP_SIZE   = 2.0e-4_wp
   real(wp), parameter :: ABS_THR = 5.0E-9_wp
   real(wp), parameter :: REL_THR = 5.0E-8_wp


contains

   subroutine collect_cavity_drop_cfc_kernel(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("atom_pd1_r_fd       ", test_atom_pd1_r_fd), &
         new_unittest("atom_pd2_rr_fd      ", test_atom_pd2_rr_fd), &
         new_unittest("atom_pd3_rrr_fd     ", test_atom_pd3_rrr_fd), &
         new_unittest("atom_pd2_symmetry   ", test_atom_pd2_symmetry), &
         new_unittest("atom_pd3_symmetry   ", test_atom_pd3_symmetry), &
         new_unittest("pair_pd1_a_fd       ", test_pair_pd1_a_fd), &
         new_unittest("pair_pd1_b_fd       ", test_pair_pd1_b_fd), &
         new_unittest("pair_pd2_aa_fd      ", test_pair_pd2_aa_fd), &
         new_unittest("pair_pd2_bb_fd      ", test_pair_pd2_bb_fd), &
         new_unittest("pair_pd2_ab_fd      ", test_pair_pd2_ab_fd), &
         new_unittest("pair_pd3_aaa_fd     ", test_pair_pd3_aaa_fd), &
         new_unittest("pair_pd3_bbb_fd     ", test_pair_pd3_bbb_fd), &
         new_unittest("pair_pd3_aab_fd     ", test_pair_pd3_aab_fd), &
         new_unittest("pair_pd3_abb_fd     ", test_pair_pd3_abb_fd), &
         new_unittest("pair_tensor_symm    ", test_pair_tensor_symmetries), &
         new_unittest("pair_swap_a_b       ", test_pair_swap_invariance), &
         new_unittest("log_lift_grad_fd    ", test_log_lift_grad_fd), &
         new_unittest("log_lift_hess_fd    ", test_log_lift_hess_fd), &
         new_unittest("log_lift_third_fd   ", test_log_lift_third_fd) &
         ]
   end subroutine collect_cavity_drop_cfc_kernel

   !* ================================================================================= *!
   !*                       Local single-call wrappers (zeroing)                        *!
   !* ================================================================================= *!

   !> Evaluate the atomic-term kernel at one point, returning a fresh
   !> tensor stack (initialised to zero so the kernel's accumulator
   !> semantics produce the single-point value).
   subroutine eval_atom(d_a, R_a, a1, max_deriv, pd0, pd1_r, pd2_rr, pd3_rrr)
      real(wp), intent(in)  :: d_a(ndim), R_a, a1
      integer,  intent(in)  :: max_deriv
      real(wp), intent(out) :: pd0
      real(wp), intent(out) :: pd1_r(ndim)
      real(wp), intent(out) :: pd2_rr(ndim, ndim)
      real(wp), intent(out) :: pd3_rrr(ndim, ndim, ndim)

      pd0 = 0.0_wp
      pd1_r = 0.0_wp
      pd2_rr = 0.0_wp
      pd3_rrr = 0.0_wp
      call cfc_atomic_term_eval(d_a, R_a, a1, max_deriv, pd0, pd1_r, pd2_rr, pd3_rrr)
   end subroutine eval_atom

   !> Evaluate the pair-term kernel at one (d_a, d_b) configuration.
   subroutine eval_pair(d_a, d_b, R_a, R_b, a2, c_par, max_deriv, &
                       pd0, pd1_a, pd1_b, pd2_aa, pd2_ab, pd2_bb, &
                       pd3_aaa, pd3_aab, pd3_abb, pd3_bbb)
      real(wp), intent(in)  :: d_a(ndim), d_b(ndim), R_a, R_b, a2, c_par
      integer,  intent(in)  :: max_deriv
      real(wp), intent(out) :: pd0
      real(wp), intent(out) :: pd1_a(ndim), pd1_b(ndim)
      real(wp), intent(out) :: pd2_aa(ndim, ndim), pd2_ab(ndim, ndim), pd2_bb(ndim, ndim)
      real(wp), intent(out) :: pd3_aaa(ndim, ndim, ndim), pd3_aab(ndim, ndim, ndim)
      real(wp), intent(out) :: pd3_abb(ndim, ndim, ndim), pd3_bbb(ndim, ndim, ndim)

      pd0 = 0.0_wp
      pd1_a = 0.0_wp;  pd1_b = 0.0_wp
      pd2_aa = 0.0_wp; pd2_ab = 0.0_wp; pd2_bb = 0.0_wp
      pd3_aaa = 0.0_wp; pd3_aab = 0.0_wp
      pd3_abb = 0.0_wp; pd3_bbb = 0.0_wp
      call cfc_pair_term_eval(d_a, d_b, R_a, R_b, a2, c_par, max_deriv, &
         pd0, pd1_a, pd1_b, pd2_aa, pd2_ab, pd2_bb, &
         pd3_aaa, pd3_aab, pd3_abb, pd3_bbb)
   end subroutine eval_pair

   !> 4-point central FD: f'(x) ~ (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12 h).
   pure real(wp) function fd4_scalar(fpp, fp, fm, fmm, h) result(df)
      real(wp), intent(in) :: fpp, fp, fm, fmm, h
      df = (-fpp + 8.0_wp*fp - 8.0_wp*fm + fmm) / (12.0_wp*h)
   end function fd4_scalar

   !* ================================================================================= *!
   !*                              Atomic-term FD tests                                 *!
   !* ================================================================================= *!

   !> pd1_r[i] = d/d(d_a,i) pd0
   subroutine test_atom_pd1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia1, axis, i
      real(wp) :: d_a(ndim), R_a, a1, h
      real(wp) :: pd0_a, pd1_a(ndim), pd2_a(ndim, ndim), pd3_a(ndim, ndim, ndim)
      real(wp) :: numeric(ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm, dummy1(ndim), dummy2(ndim, ndim), dummy3(ndim, ndim, ndim)
      real(wp) :: shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_atom_pts
         d_a = atom_pts(:, ipt)
         R_a = atom_radii(ipt)
         do ia1 = 1, size(a1_values)
            a1 = a1_values(ia1)
            call eval_atom(d_a, R_a, a1, 1, pd0_a, pd1_a, pd2_a, pd3_a)
            do axis = 1, ndim
               shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
               call eval_atom(shifted, R_a, a1, 0, f_pp, dummy1, dummy2, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) + h
               call eval_atom(shifted, R_a, a1, 0, f_p,  dummy1, dummy2, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) - h
               call eval_atom(shifted, R_a, a1, 0, f_m,  dummy1, dummy2, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
               call eval_atom(shifted, R_a, a1, 0, f_mm, dummy1, dummy2, dummy3)
               numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
            end do
            do i = 1, ndim
               call check(error, pd1_a(i), numeric(i), thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_atom_pd1_r_fd

   !> pd2_rr[i,j] = d/d(d_a,j) pd1_r[i]
   subroutine test_atom_pd2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia1, axis, i, j
      real(wp) :: d_a(ndim), R_a, a1, h
      real(wp) :: pd0_a, pd1_a(ndim), pd2_a(ndim, ndim), pd3_a(ndim, ndim, ndim)
      real(wp) :: numeric(ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: shifted(ndim)
      real(wp) :: dummy0, dummy2(ndim, ndim), dummy3(ndim, ndim, ndim)

      h = STEP_SIZE
      do ipt = 1, n_atom_pts
         d_a = atom_pts(:, ipt)
         R_a = atom_radii(ipt)
         do ia1 = 1, size(a1_values)
            a1 = a1_values(ia1)
            call eval_atom(d_a, R_a, a1, 2, pd0_a, pd1_a, pd2_a, pd3_a)
            do axis = 1, ndim
               shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
               call eval_atom(shifted, R_a, a1, 1, dummy0, g_pp, dummy2, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) + h
               call eval_atom(shifted, R_a, a1, 1, dummy0, g_p,  dummy2, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) - h
               call eval_atom(shifted, R_a, a1, 1, dummy0, g_m,  dummy2, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
               call eval_atom(shifted, R_a, a1, 1, dummy0, g_mm, dummy2, dummy3)
               do i = 1, ndim
                  numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
               end do
            end do
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd2_a(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_atom_pd2_rr_fd

   !> pd3_rrr[i,j,k] = d/d(d_a,k) pd2_rr[i,j]
   subroutine test_atom_pd3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia1, axis, i, j, k
      real(wp) :: d_a(ndim), R_a, a1, h
      real(wp) :: pd0_a, pd1_a(ndim), pd2_a(ndim, ndim), pd3_a(ndim, ndim, ndim)
      real(wp) :: numeric(ndim, ndim, ndim)
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: shifted(ndim)
      real(wp) :: dummy0, dummy1(ndim), dummy3(ndim, ndim, ndim)

      h = STEP_SIZE
      do ipt = 1, n_atom_pts
         d_a = atom_pts(:, ipt)
         R_a = atom_radii(ipt)
         do ia1 = 1, size(a1_values)
            a1 = a1_values(ia1)
            call eval_atom(d_a, R_a, a1, 3, pd0_a, pd1_a, pd2_a, pd3_a)
            do axis = 1, ndim
               shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
               call eval_atom(shifted, R_a, a1, 2, dummy0, dummy1, h_pp, dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) + h
               call eval_atom(shifted, R_a, a1, 2, dummy0, dummy1, h_p,  dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) - h
               call eval_atom(shifted, R_a, a1, 2, dummy0, dummy1, h_m,  dummy3)
               shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
               call eval_atom(shifted, R_a, a1, 2, dummy0, dummy1, h_mm, dummy3)
               do j = 1, ndim
                  do i = 1, ndim
                     numeric(i, j, axis) = fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), h_mm(i, j), h)
                  end do
               end do
            end do
            do k = 1, ndim
               do j = 1, ndim
                  do i = 1, ndim
                     call check(error, pd3_a(i, j, k), numeric(i, j, k), &
                        thr_abs=ABS_THR, thr_rel=REL_THR)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_atom_pd3_rrr_fd

   subroutine test_atom_pd2_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia1, i, j
      real(wp) :: pd0_a, pd1_a(ndim), pd2_a(ndim, ndim), pd3_a(ndim, ndim, ndim)

      do ipt = 1, n_atom_pts
         do ia1 = 1, size(a1_values)
            call eval_atom(atom_pts(:, ipt), atom_radii(ipt), a1_values(ia1), 2, &
               pd0_a, pd1_a, pd2_a, pd3_a)
            do j = 1, ndim
               do i = j + 1, ndim
                  call check(error, pd2_a(i, j), pd2_a(j, i), thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_atom_pd2_symmetry

   subroutine test_atom_pd3_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia1, i, j, k
      real(wp) :: pd0_a, pd1_a(ndim), pd2_a(ndim, ndim), pd3_a(ndim, ndim, ndim)
      real(wp) :: ref

      do ipt = 1, n_atom_pts
         do ia1 = 1, size(a1_values)
            call eval_atom(atom_pts(:, ipt), atom_radii(ipt), a1_values(ia1), 3, &
               pd0_a, pd1_a, pd2_a, pd3_a)
            do k = 1, ndim
               do j = 1, ndim
                  do i = 1, ndim
                     ref = pd3_a(i, j, k)
                     call check(error, pd3_a(i, k, j), ref, thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
                     if (allocated(error)) return
                     call check(error, pd3_a(j, i, k), ref, thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
                     if (allocated(error)) return
                     call check(error, pd3_a(j, k, i), ref, thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
                     if (allocated(error)) return
                     call check(error, pd3_a(k, i, j), ref, thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
                     if (allocated(error)) return
                     call check(error, pd3_a(k, j, i), ref, thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_atom_pd3_symmetry

   !* ================================================================================= *!
   !*                               Pair-term FD tests                                  *!
   !* ================================================================================= *!

   !> pd1_a[i] = d/d(d_a,i) pd0 (with d_b fixed)
   subroutine test_pair_pd1_a_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia2, ic, axis, i
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, a2, c_par, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: numeric(ndim), shifted(ndim)
      real(wp) :: t_pd1a(ndim), t_pd1b(ndim), t_pd2aa(ndim, ndim)
      real(wp) :: t_pd2ab(ndim, ndim), t_pd2bb(ndim, ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         do ia2 = 1, size(a2_values)
            do ic = 1, size(c_values)
               a2 = a2_values(ia2); c_par = c_values(ic)
               call eval_pair(d_a, d_b, R_a, R_b, a2, c_par, 1, pd0, pd1a, pd1b, &
                  pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
               do axis = 1, ndim
                  shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
                  call eval_pair(shifted, d_b, R_a, R_b, a2, c_par, 0, f_pp, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  shifted = d_a; shifted(axis) = d_a(axis) + h
                  call eval_pair(shifted, d_b, R_a, R_b, a2, c_par, 0, f_p, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  shifted = d_a; shifted(axis) = d_a(axis) - h
                  call eval_pair(shifted, d_b, R_a, R_b, a2, c_par, 0, f_m, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
                  call eval_pair(shifted, d_b, R_a, R_b, a2, c_par, 0, f_mm, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
               end do
               do i = 1, ndim
                  call check(error, pd1a(i), numeric(i), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_pd1_a_fd

   !> pd1_b[i] = d/d(d_b,i) pd0 (with d_a fixed)
   subroutine test_pair_pd1_b_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, ia2, ic, axis, i
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, a2, c_par, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm
      real(wp) :: numeric(ndim), shifted(ndim)
      real(wp) :: t_pd1a(ndim), t_pd1b(ndim), t_pd2aa(ndim, ndim)
      real(wp) :: t_pd2ab(ndim, ndim), t_pd2bb(ndim, ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         do ia2 = 1, size(a2_values)
            do ic = 1, size(c_values)
               a2 = a2_values(ia2); c_par = c_values(ic)
               call eval_pair(d_a, d_b, R_a, R_b, a2, c_par, 1, pd0, pd1a, pd1b, &
                  pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
               do axis = 1, ndim
                  shifted = d_b; shifted(axis) = d_b(axis) + 2.0_wp*h
                  call eval_pair(d_a, shifted, R_a, R_b, a2, c_par, 0, f_pp, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  shifted = d_b; shifted(axis) = d_b(axis) + h
                  call eval_pair(d_a, shifted, R_a, R_b, a2, c_par, 0, f_p, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  shifted = d_b; shifted(axis) = d_b(axis) - h
                  call eval_pair(d_a, shifted, R_a, R_b, a2, c_par, 0, f_m, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  shifted = d_b; shifted(axis) = d_b(axis) - 2.0_wp*h
                  call eval_pair(d_a, shifted, R_a, R_b, a2, c_par, 0, f_mm, &
                     t_pd1a, t_pd1b, t_pd2aa, t_pd2ab, t_pd2bb, &
                     t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
                  numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
               end do
               do i = 1, ndim
                  call check(error, pd1b(i), numeric(i), thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_pd1_b_fd

   !> pd2_aa[i,j] = d/d(d_a,j) pd1_a[i]
   subroutine test_pair_pd2_aa_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: dummy_b(ndim)
      real(wp) :: t_pd2aa(ndim, ndim), t_pd2ab(ndim, ndim), t_pd2bb(ndim, ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 2, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_pp, dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) + h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_p,  dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) - h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_m,  dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_mm, dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do i = 1, ndim
               numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
            end do
         end do
         do j = 1, ndim
            do i = 1, ndim
               call check(error, pd2aa(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_pair_pd2_aa_fd

   !> pd2_bb[i,j] = d/d(d_b,j) pd1_b[i]
   subroutine test_pair_pd2_bb_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: dummy_a(ndim)
      real(wp) :: t_pd2aa(ndim, ndim), t_pd2ab(ndim, ndim), t_pd2bb(ndim, ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 2, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_b; shifted(axis) = d_b(axis) + 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, dummy_a, g_pp, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) + h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, dummy_a, g_p, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, dummy_a, g_m, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, dummy_a, g_mm, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do i = 1, ndim
               numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
            end do
         end do
         do j = 1, ndim
            do i = 1, ndim
               call check(error, pd2bb(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_pair_pd2_bb_fd

   !> pd2_ab[i,j] = d/d(d_b,j) pd1_a[i] (NOT symmetric in i,j)
   subroutine test_pair_pd2_ab_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: dummy_b(ndim)
      real(wp) :: t_pd2aa(ndim, ndim), t_pd2ab(ndim, ndim), t_pd2bb(ndim, ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 2, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_b; shifted(axis) = d_b(axis) + 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_pp, dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) + h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_p,  dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_m,  dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 1, dummy_pd0, g_mm, dummy_b, &
               t_pd2aa, t_pd2ab, t_pd2bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do i = 1, ndim
               numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
            end do
         end do
         do j = 1, ndim
            do i = 1, ndim
               call check(error, pd2ab(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_pair_pd2_ab_fd

   !> pd3_aaa[i,j,k] = d/d(d_a,k) pd2_aa[i,j]
   subroutine test_pair_pd3_aaa_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j, k
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: dummy_ab(ndim, ndim), dummy_bb(ndim, ndim)
      real(wp) :: t_pd1a(ndim), t_pd1b(ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 3, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_pp, dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) + h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_p,  dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) - h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_m,  dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_mm, dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do j = 1, ndim
               do i = 1, ndim
                  numeric(i, j, axis) = fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), h_mm(i, j), h)
               end do
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3aaa(i, j, k), numeric(i, j, k), &
                     thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_pd3_aaa_fd

   !> pd3_bbb[i,j,k] = d/d(d_b,k) pd2_bb[i,j]
   subroutine test_pair_pd3_bbb_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j, k
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: dummy_aa(ndim, ndim), dummy_ab(ndim, ndim)
      real(wp) :: t_pd1a(ndim), t_pd1b(ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 3, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_b; shifted(axis) = d_b(axis) + 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_pp, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) + h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_p,  t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_m,  t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_mm, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do j = 1, ndim
               do i = 1, ndim
                  numeric(i, j, axis) = fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), h_mm(i, j), h)
               end do
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3bbb(i, j, k), numeric(i, j, k), &
                     thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_pd3_bbb_fd

   !> pd3_aab[i,j,k] = d/d(d_b,k) pd2_aa[i,j] (a-symmetric in i,j)
   subroutine test_pair_pd3_aab_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j, k
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: dummy_ab(ndim, ndim), dummy_bb(ndim, ndim)
      real(wp) :: t_pd1a(ndim), t_pd1b(ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 3, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_b; shifted(axis) = d_b(axis) + 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_pp, dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) + h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_p,  dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_m,  dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_b; shifted(axis) = d_b(axis) - 2.0_wp*h
            call eval_pair(d_a, shifted, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               h_mm, dummy_ab, dummy_bb, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do j = 1, ndim
               do i = 1, ndim
                  numeric(i, j, axis) = fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), h_mm(i, j), h)
               end do
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3aab(i, j, k), numeric(i, j, k), &
                     thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_pd3_aab_fd

   !> pd3_abb[i,j,k] = d/d(d_a,i) pd2_bb[j,k] (b-symmetric in j,k)
   subroutine test_pair_pd3_abb_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j, k
      real(wp) :: d_a(ndim), d_b(ndim), R_a, R_b, h
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: dummy_aa(ndim, ndim), dummy_ab(ndim, ndim)
      real(wp) :: t_pd1a(ndim), t_pd1b(ndim)
      real(wp) :: t_pd3aaa(ndim, ndim, ndim), t_pd3aab(ndim, ndim, ndim)
      real(wp) :: t_pd3abb(ndim, ndim, ndim), t_pd3bbb(ndim, ndim, ndim)
      real(wp) :: dummy_pd0, numeric(ndim, ndim, ndim), shifted(ndim)

      h = STEP_SIZE
      do ipt = 1, n_pair_pts
         d_a = pair_d_a(:, ipt); d_b = pair_d_b(:, ipt)
         R_a = pair_R_a(ipt);    R_b = pair_R_b(ipt)
         call eval_pair(d_a, d_b, R_a, R_b, a2_ref, c_ref, 3, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)
         do axis = 1, ndim
            shifted = d_a; shifted(axis) = d_a(axis) + 2.0_wp*h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_pp, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) + h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_p,  t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) - h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_m,  t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            shifted = d_a; shifted(axis) = d_a(axis) - 2.0_wp*h
            call eval_pair(shifted, d_b, R_a, R_b, a2_ref, c_ref, 2, dummy_pd0, t_pd1a, t_pd1b, &
               dummy_aa, dummy_ab, h_mm, t_pd3aaa, t_pd3aab, t_pd3abb, t_pd3bbb)
            do k = 1, ndim
               do j = 1, ndim
                  numeric(axis, j, k) = fd4_scalar(h_pp(j, k), h_p(j, k), h_m(j, k), h_mm(j, k), h)
               end do
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3abb(i, j, k), numeric(i, j, k), &
                     thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_pd3_abb_fd

   !> Verify the index-permutation symmetries of pd2_aa/pd2_bb (symmetric)
   !> and pd3_aaa/pd3_aab/pd3_abb/pd3_bbb (fully or partly symmetric).
   subroutine test_pair_tensor_symmetries(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, i, j, k
      real(wp) :: pd0, pd1a(ndim), pd1b(ndim)
      real(wp) :: pd2aa(ndim, ndim), pd2ab(ndim, ndim), pd2bb(ndim, ndim)
      real(wp) :: pd3aaa(ndim, ndim, ndim), pd3aab(ndim, ndim, ndim)
      real(wp) :: pd3abb(ndim, ndim, ndim), pd3bbb(ndim, ndim, ndim)
      real(wp) :: ref
      real(wp), parameter :: stol = 1.0e-13_wp

      do ipt = 1, n_pair_pts
         call eval_pair(pair_d_a(:, ipt), pair_d_b(:, ipt), pair_R_a(ipt), pair_R_b(ipt), &
            a2_ref, c_ref, 3, pd0, pd1a, pd1b, &
            pd2aa, pd2ab, pd2bb, pd3aaa, pd3aab, pd3abb, pd3bbb)

         do j = 1, ndim
            do i = j + 1, ndim
               call check(error, pd2aa(i, j), pd2aa(j, i), thr_abs=stol, thr_rel=stol)
               if (allocated(error)) return
               call check(error, pd2bb(i, j), pd2bb(j, i), thr_abs=stol, thr_rel=stol)
               if (allocated(error)) return
            end do
         end do

         !* pd3_aaa fully symmetric *!
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  ref = pd3aaa(i, j, k)
                  call check(error, pd3aaa(j, i, k), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3aaa(j, k, i), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3aaa(i, k, j), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3aaa(k, i, j), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3aaa(k, j, i), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  ref = pd3bbb(i, j, k)
                  call check(error, pd3bbb(j, i, k), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3bbb(j, k, i), ref, thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
               end do
            end do
         end do

         !* pd3_aab symmetric in first two indices *!
         do k = 1, ndim
            do j = 1, ndim
               do i = j + 1, ndim
                  call check(error, pd3aab(i, j, k), pd3aab(j, i, k), thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
               end do
            end do
         end do

         !* pd3_abb symmetric in last two indices *!
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3abb(i, j, k), pd3abb(i, k, j), thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_tensor_symmetries

   !> Swap-(a,b) invariance: swapping the two atoms in the pair kernel
   !> exchanges aa <-> bb tensors and transposes ab, with appropriate
   !> index re-labellings. Verifies the sympy expression respects the
   !> (a,b) symmetry of the underlying K_pair.
   subroutine test_pair_swap_invariance(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, i, j, k
      real(wp) :: pd0_1, pd1a_1(ndim), pd1b_1(ndim)
      real(wp) :: pd2aa_1(ndim, ndim), pd2ab_1(ndim, ndim), pd2bb_1(ndim, ndim)
      real(wp) :: pd3aaa_1(ndim, ndim, ndim), pd3aab_1(ndim, ndim, ndim)
      real(wp) :: pd3abb_1(ndim, ndim, ndim), pd3bbb_1(ndim, ndim, ndim)
      real(wp) :: pd0_2, pd1a_2(ndim), pd1b_2(ndim)
      real(wp) :: pd2aa_2(ndim, ndim), pd2ab_2(ndim, ndim), pd2bb_2(ndim, ndim)
      real(wp) :: pd3aaa_2(ndim, ndim, ndim), pd3aab_2(ndim, ndim, ndim)
      real(wp) :: pd3abb_2(ndim, ndim, ndim), pd3bbb_2(ndim, ndim, ndim)
      real(wp), parameter :: stol = 1.0e-12_wp

      do ipt = 1, n_pair_pts
         call eval_pair(pair_d_a(:, ipt), pair_d_b(:, ipt), pair_R_a(ipt), pair_R_b(ipt), &
            a2_ref, c_ref, 3, pd0_1, pd1a_1, pd1b_1, &
            pd2aa_1, pd2ab_1, pd2bb_1, pd3aaa_1, pd3aab_1, pd3abb_1, pd3bbb_1)
         call eval_pair(pair_d_b(:, ipt), pair_d_a(:, ipt), pair_R_b(ipt), pair_R_a(ipt), &
            a2_ref, c_ref, 3, pd0_2, pd1a_2, pd1b_2, &
            pd2aa_2, pd2ab_2, pd2bb_2, pd3aaa_2, pd3aab_2, pd3abb_2, pd3bbb_2)

         call check(error, pd0_1, pd0_2, thr_abs=stol, thr_rel=stol)
         if (allocated(error)) return
         do i = 1, ndim
            call check(error, pd1a_1(i), pd1b_2(i), thr_abs=stol, thr_rel=stol)
            if (allocated(error)) return
            call check(error, pd1b_1(i), pd1a_2(i), thr_abs=stol, thr_rel=stol)
            if (allocated(error)) return
         end do
         do j = 1, ndim
            do i = 1, ndim
               call check(error, pd2aa_1(i, j), pd2bb_2(i, j), thr_abs=stol, thr_rel=stol)
               if (allocated(error)) return
               call check(error, pd2bb_1(i, j), pd2aa_2(i, j), thr_abs=stol, thr_rel=stol)
               if (allocated(error)) return
               !* pd2_ab(i,j) under swap becomes pd2_ab(j,i): the a-index
               !* (now from the original b) lives in slot j of the swapped
               !* expression, and vice versa.
               call check(error, pd2ab_1(i, j), pd2ab_2(j, i), thr_abs=stol, thr_rel=stol)
               if (allocated(error)) return
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3aaa_1(i, j, k), pd3bbb_2(i, j, k), thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3bbb_1(i, j, k), pd3aaa_2(i, j, k), thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  !* pd3_aab(i,j,k) (a-symmetric (i,j), b-singleton k) swaps to
                  !* pd3_abb(k,i,j) (a-singleton k, b-symmetric (i,j)).
                  call check(error, pd3aab_1(i, j, k), pd3abb_2(k, i, j), thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3abb_1(i, j, k), pd3aab_2(j, k, i), thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_swap_invariance

   !* ================================================================================= *!
   !*                            log_lift FD tests                                      *!
   !* ================================================================================= *!

   !> Reference scalar field PD(r) = sum of two atom-centered Gaussians.
   !> Used as a known smooth positive PD to verify cfc_log_lift against
   !> finite differences of log(PD).
   pure subroutine ref_pd(r, pd0, pd1, pd2, pd3)
      real(wp), intent(in) :: r(ndim)
      real(wp), intent(out) :: pd0, pd1(ndim), pd2(ndim, ndim), pd3(ndim, ndim, ndim)
      real(wp), parameter :: c1(ndim) = [0.7_wp, -0.3_wp, 0.4_wp]
      real(wp), parameter :: c2(ndim) = [-0.5_wp, 0.6_wp, -0.2_wp]
      real(wp), parameter :: alpha1 = 1.3_wp
      real(wp), parameter :: alpha2 = 0.9_wp
      real(wp) :: d1(ndim), d2(ndim)
      real(wp) :: e1, e2
      integer  :: i, j, k

      d1 = r - c1; d2 = r - c2
      e1 = exp(-alpha1*dot_product(d1, d1))
      e2 = exp(-alpha2*dot_product(d2, d2))
      pd0 = e1 + e2

      do i = 1, ndim
         pd1(i) = -2.0_wp*alpha1*d1(i)*e1 - 2.0_wp*alpha2*d2(i)*e2
      end do

      do j = 1, ndim
         do i = 1, ndim
            pd2(i, j) = (4.0_wp*alpha1**2 * d1(i)*d1(j) &
                       - merge(2.0_wp*alpha1, 0.0_wp, i == j)) * e1 &
                      + (4.0_wp*alpha2**2 * d2(i)*d2(j) &
                       - merge(2.0_wp*alpha2, 0.0_wp, i == j)) * e2
         end do
      end do

      do k = 1, ndim
         do j = 1, ndim
            do i = 1, ndim
               pd3(i, j, k) = ( &
                   -8.0_wp*alpha1**3 * d1(i)*d1(j)*d1(k) &
                   + 4.0_wp*alpha1**2 * ( &
                        merge(d1(k), 0.0_wp, i == j) + &
                        merge(d1(j), 0.0_wp, i == k) + &
                        merge(d1(i), 0.0_wp, j == k))) * e1 &
                 + ( &
                   -8.0_wp*alpha2**3 * d2(i)*d2(j)*d2(k) &
                   + 4.0_wp*alpha2**2 * ( &
                        merge(d2(k), 0.0_wp, i == j) + &
                        merge(d2(j), 0.0_wp, i == k) + &
                        merge(d2(i), 0.0_wp, j == k))) * e2
            end do
         end do
      end do
   end subroutine ref_pd

   pure subroutine ref_log_pd_value(r, val)
      real(wp), intent(in) :: r(ndim)
      real(wp), intent(out) :: val
      real(wp) :: pd0, pd1(ndim), pd2(ndim, ndim), pd3(ndim, ndim, ndim)
      call ref_pd(r, pd0, pd1, pd2, pd3)
      val = log(pd0)
   end subroutine ref_log_pd_value

   subroutine test_log_lift_grad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i
      real(wp) :: r(ndim), pd0, pd1(ndim), pd2(ndim, ndim), pd3(ndim, ndim, ndim)
      real(wp) :: lpd0, lpd1(ndim), lpd2(ndim, ndim), lpd3(ndim, ndim, ndim)
      real(wp) :: f_pp, f_p, f_m, f_mm, h
      real(wp) :: shifted(ndim), numeric(ndim)

      h = STEP_SIZE
      do ipt = 1, n_atom_pts
         r = atom_pts(:, ipt)
         call ref_pd(r, pd0, pd1, pd2, pd3)
         call cfc_log_lift(pd0, pd1, pd2, pd3, 1, lpd0, lpd1, lpd2, lpd3)
         do axis = 1, ndim
            shifted = r; shifted(axis) = r(axis) + 2.0_wp*h
            call ref_log_pd_value(shifted, f_pp)
            shifted = r; shifted(axis) = r(axis) + h
            call ref_log_pd_value(shifted, f_p)
            shifted = r; shifted(axis) = r(axis) - h
            call ref_log_pd_value(shifted, f_m)
            shifted = r; shifted(axis) = r(axis) - 2.0_wp*h
            call ref_log_pd_value(shifted, f_mm)
            numeric(axis) = fd4_scalar(f_pp, f_p, f_m, f_mm, h)
         end do
         do i = 1, ndim
            call check(error, lpd1(i), numeric(i), thr_abs=ABS_THR, thr_rel=REL_THR)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_log_lift_grad_fd

   subroutine test_log_lift_hess_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j
      real(wp) :: r(ndim), pd0, pd1(ndim), pd2(ndim, ndim), pd3(ndim, ndim, ndim)
      real(wp) :: lpd0, lpd1(ndim), lpd2(ndim, ndim), lpd3(ndim, ndim, ndim)
      real(wp) :: pd0_s, pd1_s(ndim), pd2_s(ndim, ndim), pd3_s(ndim, ndim, ndim)
      real(wp) :: g_pp(ndim), g_p(ndim), g_m(ndim), g_mm(ndim)
      real(wp) :: lpd1_s(ndim), lpd2_s(ndim, ndim), lpd3_s(ndim, ndim, ndim), lpd0_s
      real(wp) :: h, shifted(ndim), numeric(ndim, ndim)

      h = STEP_SIZE
      do ipt = 1, n_atom_pts
         r = atom_pts(:, ipt)
         call ref_pd(r, pd0, pd1, pd2, pd3)
         call cfc_log_lift(pd0, pd1, pd2, pd3, 2, lpd0, lpd1, lpd2, lpd3)
         do axis = 1, ndim
            shifted = r; shifted(axis) = r(axis) + 2.0_wp*h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 1, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            g_pp = lpd1_s
            shifted = r; shifted(axis) = r(axis) + h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 1, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            g_p = lpd1_s
            shifted = r; shifted(axis) = r(axis) - h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 1, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            g_m = lpd1_s
            shifted = r; shifted(axis) = r(axis) - 2.0_wp*h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 1, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            g_mm = lpd1_s
            do i = 1, ndim
               numeric(i, axis) = fd4_scalar(g_pp(i), g_p(i), g_m(i), g_mm(i), h)
            end do
         end do
         do j = 1, ndim
            do i = 1, ndim
               call check(error, lpd2(i, j), numeric(i, j), thr_abs=ABS_THR, thr_rel=REL_THR)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_log_lift_hess_fd

   subroutine test_log_lift_third_fd(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ipt, axis, i, j, k
      real(wp) :: r(ndim), pd0, pd1(ndim), pd2(ndim, ndim), pd3(ndim, ndim, ndim)
      real(wp) :: lpd0, lpd1(ndim), lpd2(ndim, ndim), lpd3(ndim, ndim, ndim)
      real(wp) :: pd0_s, pd1_s(ndim), pd2_s(ndim, ndim), pd3_s(ndim, ndim, ndim)
      real(wp) :: lpd0_s, lpd1_s(ndim), lpd2_s(ndim, ndim), lpd3_s(ndim, ndim, ndim)
      real(wp) :: h_pp(ndim, ndim), h_p(ndim, ndim), h_m(ndim, ndim), h_mm(ndim, ndim)
      real(wp) :: h, shifted(ndim), numeric(ndim, ndim, ndim)

      h = STEP_SIZE
      do ipt = 1, n_atom_pts
         r = atom_pts(:, ipt)
         call ref_pd(r, pd0, pd1, pd2, pd3)
         call cfc_log_lift(pd0, pd1, pd2, pd3, 3, lpd0, lpd1, lpd2, lpd3)
         do axis = 1, ndim
            shifted = r; shifted(axis) = r(axis) + 2.0_wp*h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 2, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            h_pp = lpd2_s
            shifted = r; shifted(axis) = r(axis) + h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 2, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            h_p = lpd2_s
            shifted = r; shifted(axis) = r(axis) - h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 2, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            h_m = lpd2_s
            shifted = r; shifted(axis) = r(axis) - 2.0_wp*h
            call ref_pd(shifted, pd0_s, pd1_s, pd2_s, pd3_s)
            call cfc_log_lift(pd0_s, pd1_s, pd2_s, pd3_s, 2, lpd0_s, lpd1_s, lpd2_s, lpd3_s)
            h_mm = lpd2_s
            do j = 1, ndim
               do i = 1, ndim
                  numeric(i, j, axis) = fd4_scalar(h_pp(i, j), h_p(i, j), h_m(i, j), h_mm(i, j), h)
               end do
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, lpd3(i, j, k), numeric(i, j, k), &
                     thr_abs=ABS_THR, thr_rel=REL_THR)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_log_lift_third_fd

end module test_cavity_drop_cfc_kernel
