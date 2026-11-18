!> Kernel-level unit tests for the auto-generated CFC pseudo-density module
module test_cavity_drop_cfc
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use moist_api, only: vp_error
   use moist_cavity_drop_lsf_cfc_kernel, only: cfc_atomic_term_eval, &
                                               cfc_pair_term_eval, cfc_spatial_eval, &
                                               cfc_nuclear_eval, cfc_hessian_eval, &
                                               cfc_tangent_eval, cfc_hvp_eval
   use test_helpers, only: fd_probe_type, fd_derivative, fd_compare
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, to_string
   implicit none(type, external)
   private

   public :: collect_cavity_drop_cfc

   integer, parameter :: ndim = 3

   !> Fade band at infinity in sigma: every center on the unfaded kernel body
   real(wp), parameter :: no_fade_start = huge(1.0_wp), no_fade_width = 1.0_wp

   !> Klamt-Diedenhofen 2018 reference parameter set
   real(wp), parameter :: a1_ref = -15.0_wp
   real(wp), parameter :: a2_ref = -9.0_wp
   real(wp), parameter :: c_ref = 5.0_wp

   !> Auxiliary parameter sweeps (the kernel is parametric in a1, a2, c)
   real(wp), parameter :: a1_values(2) = [-15.0_wp, -8.0_wp]
   real(wp), parameter :: a2_values(2) = [-9.0_wp, -5.0_wp]
   real(wp), parameter :: c_values(2) = [5.0_wp, 2.5_wp]

   !> Atomic-term test triples (d_a, R_a)
   integer, parameter :: n_atom_pts = 6
   real(wp), parameter :: atom_pts(ndim, n_atom_pts) = reshape([ &
                                                               0.60_wp, 0.20_wp, -0.10_wp, &
                                                               1.10_wp, -0.40_wp, 0.30_wp, &
                                                               -0.50_wp, 0.90_wp, 0.20_wp, &
                                                               0.80_wp, 0.80_wp, 0.80_wp, &
                                                               -1.30_wp, -0.40_wp, 0.70_wp, &
                                                               0.30_wp, -1.20_wp, 0.50_wp], &
                                                               [ndim, n_atom_pts])
   real(wp), parameter :: atom_radii(n_atom_pts) = [ &
                          0.90_wp, 1.10_wp, 1.30_wp, 0.85_wp, 1.50_wp, 1.00_wp]

   !> Pair-term test quadruples (d_a, d_b, R_a, R_b)
   integer, parameter :: n_pair_pts = 5
   real(wp), parameter :: pair_d_a(ndim, n_pair_pts) = reshape([ &
                                                               0.50_wp, 0.30_wp, 0.10_wp, &
                                                               1.00_wp, -0.20_wp, 0.40_wp, &
                                                               -0.40_wp, 0.80_wp, 0.20_wp, &
                                                               0.70_wp, -0.50_wp, -0.30_wp, &
                                                               -0.90_wp, -0.30_wp, 0.50_wp], &
                                                               [ndim, n_pair_pts])
   real(wp), parameter :: pair_d_b(ndim, n_pair_pts) = reshape([ &
                                                               -0.40_wp, 0.20_wp, 0.30_wp, &
                                                               0.30_wp, 0.90_wp, -0.10_wp, &
                                                               0.60_wp, -0.20_wp, 0.80_wp, &
                                                               -0.60_wp, 0.40_wp, 0.30_wp, &
                                                               0.70_wp, -0.50_wp, -0.20_wp], &
                                                               [ndim, n_pair_pts])
   real(wp), parameter :: pair_R_a(n_pair_pts) = [0.90_wp, 1.10_wp, 1.30_wp, 0.95_wp, 1.20_wp]
   real(wp), parameter :: pair_R_b(n_pair_pts) = [1.20_wp, 0.85_wp, 1.05_wp, 1.40_wp, 0.95_wp]

   !> FD steps and thresholds
   real(wp), parameter :: STEP_SIZE = 2.0e-4_wp
   real(wp), parameter :: ABS_THR = 5.0E-9_wp
   real(wp), parameter :: REL_THR = 5.0E-8_wp

   !> Kernel term of an FD check
   integer, parameter :: term_atom = 1, term_pair = 2, term_lift = 3

   !> One kernel output with offset `d_a` or `d_b` displaced along one axis
   type, extends(fd_probe_type) :: cfc_kernel_probe_type
      !> Kernel term, `term_*`
      integer :: term = term_atom
      !> Sampled order, 0..2
      integer :: order = 0
      !> Pair: sample the tensor of center b (`pd1_b`, `pd2_bb`), else of center a
      logical :: side_b = .false.
      !> Pair: displace `d_b`, else `d_a`
      logical :: move_b = .false.
      !> Displaced axis
      integer :: axis = 1
      !> Offset to center a, or the lift point
      real(wp) :: d_a(ndim) = 0.0_wp
      !> Offset to center b
      real(wp) :: d_b(ndim) = 0.0_wp
      !> Radius of center a
      real(wp) :: R_a = 1.0_wp
      !> Radius of center b
      real(wp) :: R_b = 1.0_wp
      !> Exponent, `a1` (atom) or `a2` (pair)
      real(wp) :: expo = 0.0_wp
      !> Pair parameter `c`
      real(wp) :: c_par = 0.0_wp
   contains
      !> Sampled output at one displacement
      procedure :: sample => sample_cfc_kernel
   end type cfc_kernel_probe_type

contains

   subroutine collect_cavity_drop_cfc(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("atom_pd1_r_fd", test_atom_pd1_r_fd), &
                  new_unittest("atom_pd2_rr_fd", test_atom_pd2_rr_fd), &
                  new_unittest("atom_pd3_rrr_fd", test_atom_pd3_rrr_fd), &
                  new_unittest("atom_pd2_symmetry", test_atom_pd2_symmetry), &
                  new_unittest("atom_pd3_symmetry", test_atom_pd3_symmetry), &
                  new_unittest("pair_pd1_a_fd", test_pair_pd1_a_fd), &
                  new_unittest("pair_pd1_b_fd", test_pair_pd1_b_fd), &
                  new_unittest("pair_pd2_aa_fd", test_pair_pd2_aa_fd), &
                  new_unittest("pair_pd2_bb_fd", test_pair_pd2_bb_fd), &
                  new_unittest("pair_pd2_ab_fd", test_pair_pd2_ab_fd), &
                  new_unittest("pair_pd3_aaa_fd", test_pair_pd3_aaa_fd), &
                  new_unittest("pair_pd3_bbb_fd", test_pair_pd3_bbb_fd), &
                  new_unittest("pair_pd3_aab_fd", test_pair_pd3_aab_fd), &
                  new_unittest("pair_pd3_abb_fd", test_pair_pd3_abb_fd), &
                  new_unittest("pair_tensor_symm", test_pair_tensor_symmetries), &
                  new_unittest("pair_swap_a_b", test_pair_swap_invariance), &
                  new_unittest("log_lift_grad_fd", test_log_lift_grad_fd), &
                  new_unittest("log_lift_hess_fd", test_log_lift_hess_fd), &
                  new_unittest("log_lift_third_fd", test_log_lift_third_fd), &
                  new_unittest("lift_defines_every_output", test_lift_defines_every_output) &
                  ]
   end subroutine collect_cavity_drop_cfc

   !* ================================================================================= *!
   !*                       Local single-call wrappers (zeroing)                        *!
   !* ================================================================================= *!

   !> Evaluate the atomic-term kernel at one point, returning a fresh
   !> tensor stack (initialised to zero so the kernel's accumulator
   !> semantics produce the single-point value)
   subroutine eval_atom(d_a, R_a, a1, max_deriv, pd0, pd1_r, pd2_rr, pd3_rrr)
      real(wp), intent(in)  :: d_a(ndim), R_a, a1
      integer, intent(in)  :: max_deriv
      real(wp), intent(out) :: pd0
      real(wp), intent(out) :: pd1_r(ndim)
      real(wp), intent(out) :: pd2_rr(ndim, ndim)
      real(wp), intent(out) :: pd3_rrr(ndim, ndim, ndim)

      !> Order-4 accumulator of the kernel; this suite checks orders 0..3 only
      real(wp) :: pd4_rrrr(ndim, ndim, ndim, ndim)

      pd0 = 0.0_wp
      pd1_r = 0.0_wp
      pd2_rr = 0.0_wp
      pd3_rrr = 0.0_wp
      pd4_rrrr = 0.0_wp
      call cfc_atomic_term_eval(d_a, R_a, a1, no_fade_start, no_fade_width, max_deriv, &
                                pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr)
   end subroutine eval_atom

   !> Evaluate the pair-term kernel at one (d_a, d_b) configuration
   subroutine eval_pair(d_a, d_b, R_a, R_b, a2, c_par, max_deriv, &
                        pd0, pd1_a, pd1_b, pd2_aa, pd2_ab, pd2_bb, &
                        pd3_aaa, pd3_aab, pd3_abb, pd3_bbb)
      real(wp), intent(in)  :: d_a(ndim), d_b(ndim), R_a, R_b, a2, c_par
      integer, intent(in)  :: max_deriv
      real(wp), intent(out) :: pd0
      real(wp), intent(out) :: pd1_a(ndim), pd1_b(ndim)
      real(wp), intent(out) :: pd2_aa(ndim, ndim), pd2_ab(ndim, ndim), pd2_bb(ndim, ndim)
      real(wp), intent(out) :: pd3_aaa(ndim, ndim, ndim), pd3_aab(ndim, ndim, ndim)
      real(wp), intent(out) :: pd3_abb(ndim, ndim, ndim), pd3_bbb(ndim, ndim, ndim)

      !> Order-4 split accumulators; this suite checks total orders 0..3 only
      real(wp) :: pd4_aaaa(ndim, ndim, ndim, ndim), pd4_aaab(ndim, ndim, ndim, ndim)
      real(wp) :: pd4_aabb(ndim, ndim, ndim, ndim), pd4_abbb(ndim, ndim, ndim, ndim)
      real(wp) :: pd4_bbbb(ndim, ndim, ndim, ndim)

      pd0 = 0.0_wp
      pd1_a = 0.0_wp; pd1_b = 0.0_wp
      pd2_aa = 0.0_wp; pd2_ab = 0.0_wp; pd2_bb = 0.0_wp
      pd3_aaa = 0.0_wp; pd3_aab = 0.0_wp
      pd3_abb = 0.0_wp; pd3_bbb = 0.0_wp
      pd4_aaaa = 0.0_wp; pd4_aaab = 0.0_wp; pd4_aabb = 0.0_wp
      pd4_abbb = 0.0_wp; pd4_bbbb = 0.0_wp
      call cfc_pair_term_eval(d_a, d_b, R_a, R_b, a2, c_par, no_fade_start, no_fade_width, &
                              max_deriv, pd0, pd1_a, pd1_b, pd2_aa, pd2_ab, pd2_bb, &
                              pd3_aaa, pd3_aab, pd3_abb, pd3_bbb, &
                              pd4_aaaa, pd4_aaab, pd4_aabb, pd4_abbb, pd4_bbbb)
   end subroutine eval_pair

   !> Lift an accumulated pseudo-density stack to derivatives of `+log PD`
   !>
   !> The kernel's `cfc_spatial_eval` returns derivatives of the level set
   !> `S = -log PD` (the sign the LSF contract wants); the checks below are
   !> written against `+log PD`, so this wrapper negates. Orders above
   !> `max_deriv` come back as exact zeros from the kernel, not as garbage,
   !> which is what lets this wrapper negate all four unconditionally
   subroutine eval_log_lift(pd0, pd1_r, pd2_rr, pd3_rrr, max_deriv, &
                            lpd0, lpd1_r, lpd2_rr, lpd3_rrr)
      real(wp), intent(in)  :: pd0
      real(wp), intent(in)  :: pd1_r(ndim)
      real(wp), intent(in)  :: pd2_rr(ndim, ndim)
      real(wp), intent(in)  :: pd3_rrr(ndim, ndim, ndim)
      integer, intent(in)  :: max_deriv
      real(wp), intent(out) :: lpd0
      real(wp), intent(out) :: lpd1_r(ndim)
      real(wp), intent(out) :: lpd2_rr(ndim, ndim)
      real(wp), intent(out) :: lpd3_rrr(ndim, ndim, ndim)

      !> Order-4 input and output of the kernel, unused here
      real(wp) :: pd4_rrrr(ndim, ndim, ndim, ndim), f4_rrrr(ndim, ndim, ndim, ndim)

      pd4_rrrr = 0.0_wp
      call cfc_spatial_eval(pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr, max_deriv, &
                            lpd0, lpd1_r, lpd2_rr, lpd3_rrr, f4_rrrr)
      lpd0 = -lpd0
      lpd1_r = -lpd1_r
      lpd2_rr = -lpd2_rr
      lpd3_rrr = -lpd3_rrr
   end subroutine eval_log_lift

   !* ================================================================================= *!
   !*                         Finite-difference probe and checks                        *!
   !* ================================================================================= *!

   !> Analytic order `order + 1` of one kernel output against the FD of order `order`
   !>
   !> - one [[fd_derivative]] per displaced axis, one [[fd_compare]] per axis
   !>
   !> @param[out]    error  test failure, first in axis order
   !> @param[in,out] probe  kernel, term, sampled output and undisplaced arguments
   !> @param[in]     label  checked derivative and case
   subroutine check_kernel_fd(error, probe, label)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Probe
      type(cfc_kernel_probe_type), intent(inout) :: probe
      !> Checked derivative and case
      character(len=*), intent(in) :: label

      !> Analytic derivative [sampled size, axis]
      real(wp) :: analytic(ndim**probe%order, ndim)
      real(wp) :: numeric(ndim**probe%order)
      logical :: valid(ndim**probe%order)
      integer :: axis

      call kernel_analytic(probe, analytic)
      do axis = 1, ndim
         probe%axis = axis
         call fd_derivative(probe, STEP_SIZE, 4, numeric, valid, error)
         if (allocated(error)) return
         call fd_compare(error, analytic(:, axis), numeric, valid, ABS_THR, REL_THR, &
                         label//" axis "//to_string(axis))
         if (allocated(error)) return
      end do
   end subroutine check_kernel_fd

   !> Atomic-term FD over the test triples and `a1_values`
   !>
   !> @param[out] error  test failure
   !> @param[in]  name   checked derivative
   !> @param[in]  order  sampled order, 0..2
   subroutine run_atom_fd(error, name, order)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Checked derivative
      character(len=*), intent(in) :: name
      !> Sampled order
      integer, intent(in) :: order

      type(cfc_kernel_probe_type) :: probe
      integer :: ipt, ia1

      probe%term = term_atom
      probe%order = order
      do ipt = 1, n_atom_pts
         probe%d_a = atom_pts(:, ipt)
         probe%R_a = atom_radii(ipt)
         do ia1 = 1, size(a1_values)
            probe%expo = a1_values(ia1)
            call check_kernel_fd(error, probe, name//" vs FD, point "//to_string(ipt)// &
                                 " a1 "//to_string(ia1))
            if (allocated(error)) return
         end do
      end do
   end subroutine run_atom_fd

   !> Pair-term FD over the test quadruples, at the reference parameters or swept
   !>
   !> @param[out] error   test failure
   !> @param[in]  name    checked derivative
   !> @param[in]  order   sampled order, 0..2
   !> @param[in]  side_b  sample the tensor of center b (`pd1_b`, `pd2_bb`)
   !> @param[in]  move_b  displace `d_b`, else `d_a`
   !> @param[in]  sweep   `a2_values` x `c_values`, else `a2_ref`, `c_ref`
   subroutine run_pair_fd(error, name, order, side_b, move_b, sweep)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Checked derivative
      character(len=*), intent(in) :: name
      !> Sampled order
      integer, intent(in) :: order
      !> Tensor of center b
      logical, intent(in) :: side_b
      !> `d_b` displaced
      logical, intent(in) :: move_b
      !> Parameter sweep
      logical, intent(in) :: sweep

      type(cfc_kernel_probe_type) :: probe
      character(len=:), allocatable :: label
      integer :: ipt, ia2, ic

      probe%term = term_pair
      probe%order = order
      probe%side_b = side_b
      probe%move_b = move_b
      probe%expo = a2_ref
      probe%c_par = c_ref
      do ipt = 1, n_pair_pts
         probe%d_a = pair_d_a(:, ipt)
         probe%d_b = pair_d_b(:, ipt)
         probe%R_a = pair_R_a(ipt)
         probe%R_b = pair_R_b(ipt)
         label = name//" vs FD, point "//to_string(ipt)
         if (.not. sweep) then
            call check_kernel_fd(error, probe, label)
            if (allocated(error)) return
            cycle
         end if
         do ia2 = 1, size(a2_values)
            do ic = 1, size(c_values)
               probe%expo = a2_values(ia2)
               probe%c_par = c_values(ic)
               call check_kernel_fd(error, probe, label//" a2 "//to_string(ia2)//" c "//to_string(ic))
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine run_pair_fd

   !> Log-lift FD of [[ref_pd]] over the atomic test points
   !>
   !> @param[out] error  test failure
   !> @param[in]  name   checked derivative
   !> @param[in]  order  sampled order, 0..2; order 0 samples `log PD` directly
   subroutine run_lift_fd(error, name, order)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Checked derivative
      character(len=*), intent(in) :: name
      !> Sampled order
      integer, intent(in) :: order

      type(cfc_kernel_probe_type) :: probe
      integer :: ipt

      probe%term = term_lift
      probe%order = order
      do ipt = 1, n_atom_pts
         probe%d_a = atom_pts(:, ipt)
         call check_kernel_fd(error, probe, name//" vs FD, point "//to_string(ipt))
         if (allocated(error)) return
      end do
   end subroutine run_lift_fd

   !> Sampled kernel output with `d_a` or `d_b` displaced by `t` along the probe's axis
   !>
   !> @param[in,out] self   probe
   !> @param[in]     t      displacement
   !> @param[out]    vals   sampled tensor, column-major
   !> @param[out]    valid  every value present
   !> @param[out]    error  never set
   subroutine sample_cfc_kernel(self, t, vals, valid, error)
      !> Probe
      class(cfc_kernel_probe_type), intent(inout) :: self
      !> Displacement
      real(wp), intent(in) :: t
      !> Sampled values
      real(wp), intent(out) :: vals(:)
      !> Value present
      logical, intent(out) :: valid(:)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: d_a(ndim), d_b(ndim), t0, t1a(ndim), t1b(ndim)
      real(wp) :: t2aa(ndim, ndim), t2ab(ndim, ndim), t2bb(ndim, ndim)
      real(wp) :: t3aaa(ndim, ndim, ndim), t3aab(ndim, ndim, ndim)
      real(wp) :: t3abb(ndim, ndim, ndim), t3bbb(ndim, ndim, ndim)

      d_a = self%d_a
      d_b = self%d_b
      if (self%move_b) then
         d_b(self%axis) = self%d_b(self%axis) + t
      else
         d_a(self%axis) = self%d_a(self%axis) + t
      end if
      valid = .true.
      if (self%term == term_lift .and. self%order == 0) then
         call ref_log_pd_value(d_a, vals(1))
         return
      end if
      call kernel_eval(self, d_a, d_b, self%order, t0, t1a, t1b, t2aa, t2ab, t2bb, &
                       t3aaa, t3aab, t3abb, t3bbb)
      select case (self%order)
      case (0)
         vals(1) = t0
      case (1)
         if (self%side_b) then
            vals = t1b
         else
            vals = t1a
         end if
      case default
         if (self%side_b) then
            vals = reshape(t2bb, [ndim**2])
         else
            vals = reshape(t2aa, [ndim**2])
         end if
      end select
   end subroutine sample_cfc_kernel

   !> Analytic derivative of the probe's sampled output along each axis of the displaced offset
   !>
   !> @param[in]  probe     probe at its undisplaced arguments
   !> @param[out] analytic  derivative [sampled size, axis]
   subroutine kernel_analytic(probe, analytic)
      !> Probe
      type(cfc_kernel_probe_type), intent(in) :: probe
      !> Derivative
      real(wp), intent(out) :: analytic(:, :)

      real(wp) :: t0, t1a(ndim), t1b(ndim), t2aa(ndim, ndim), t2ab(ndim, ndim), t2bb(ndim, ndim)
      real(wp) :: t3aaa(ndim, ndim, ndim), t3aab(ndim, ndim, ndim)
      real(wp) :: t3abb(ndim, ndim, ndim), t3bbb(ndim, ndim, ndim)
      integer :: axis

      call kernel_eval(probe, probe%d_a, probe%d_b, probe%order + 1, t0, t1a, t1b, t2aa, t2ab, &
                       t2bb, t3aaa, t3aab, t3abb, t3bbb)
      select case (probe%order)
      case (0)
         if (probe%move_b) then
            analytic(1, :) = t1b
         else
            analytic(1, :) = t1a
         end if
      case (1)
         if (probe%side_b) then
            analytic = t2bb
         else if (probe%move_b) then
            analytic = t2ab
         else
            analytic = t2aa
         end if
      case default
         if (probe%side_b .and. .not. probe%move_b) then
            do axis = 1, ndim
               analytic(:, axis) = reshape(t3abb(axis, :, :), [ndim**2])
            end do
         else if (probe%side_b) then
            analytic = reshape(t3bbb, [ndim**2, ndim])
         else if (probe%move_b) then
            analytic = reshape(t3aab, [ndim**2, ndim])
         else
            analytic = reshape(t3aaa, [ndim**2, ndim])
         end if
      end select
   end subroutine kernel_analytic

   !> Kernel outputs of the probe's term up to `max_deriv`, pair layout
   !>
   !> - atom and lift fill the center-a tensors only, at offset or point `d_a`
   !>
   !> @param[in]  probe      kernel term and parameters
   !> @param[in]  d_a        offset to center a, or the lift point
   !> @param[in]  d_b        offset to center b (pair)
   !> @param[in]  max_deriv  highest order
   !> @param[out] t0         value
   !> @param[out] t1a        gradient in `d_a`
   !> @param[out] t1b        gradient in `d_b` (pair)
   !> @param[out] t2aa       second derivative in `d_a`
   !> @param[out] t2ab       mixed second derivative (pair)
   !> @param[out] t2bb       second derivative in `d_b` (pair)
   !> @param[out] t3aaa      third derivative in `d_a`
   !> @param[out] t3aab      `aab` third derivative (pair)
   !> @param[out] t3abb      `abb` third derivative (pair)
   !> @param[out] t3bbb      third derivative in `d_b` (pair)
   subroutine kernel_eval(probe, d_a, d_b, max_deriv, t0, t1a, t1b, t2aa, t2ab, t2bb, &
                          t3aaa, t3aab, t3abb, t3bbb)
      !> Kernel term and parameters
      class(cfc_kernel_probe_type), intent(in) :: probe
      !> Offset to center a, or the lift point
      real(wp), intent(in) :: d_a(ndim)
      !> Offset to center b
      real(wp), intent(in) :: d_b(ndim)
      !> Highest order
      integer, intent(in) :: max_deriv
      !> Value
      real(wp), intent(out) :: t0
      !> Gradient in `d_a`
      real(wp), intent(out) :: t1a(ndim)
      !> Gradient in `d_b`
      real(wp), intent(out) :: t1b(ndim)
      !> Second derivative in `d_a`
      real(wp), intent(out) :: t2aa(ndim, ndim)
      !> Mixed second derivative
      real(wp), intent(out) :: t2ab(ndim, ndim)
      !> Second derivative in `d_b`
      real(wp), intent(out) :: t2bb(ndim, ndim)
      !> Third derivative in `d_a`
      real(wp), intent(out) :: t3aaa(ndim, ndim, ndim)
      !> `aab` third derivative
      real(wp), intent(out) :: t3aab(ndim, ndim, ndim)
      !> `abb` third derivative
      real(wp), intent(out) :: t3abb(ndim, ndim, ndim)
      !> Third derivative in `d_b`
      real(wp), intent(out) :: t3bbb(ndim, ndim, ndim)

      real(wp) :: pd0, pd1(ndim), pd2(ndim, ndim), pd3(ndim, ndim, ndim)

      select case (probe%term)
      case (term_atom)
         call eval_atom(d_a, probe%R_a, probe%expo, max_deriv, t0, t1a, t2aa, t3aaa)
      case (term_pair)
         call eval_pair(d_a, d_b, probe%R_a, probe%R_b, probe%expo, probe%c_par, max_deriv, &
                        t0, t1a, t1b, t2aa, t2ab, t2bb, t3aaa, t3aab, t3abb, t3bbb)
      case default
         call ref_pd(d_a, pd0, pd1, pd2, pd3)
         call eval_log_lift(pd0, pd1, pd2, pd3, max_deriv, t0, t1a, t2aa, t3aaa)
      end select
   end subroutine kernel_eval

   !* ================================================================================= *!
   !*                              Atomic-term FD tests                                 *!
   !* ================================================================================= *!

   !> pd1_r[i] = d/d(d_a,i) pd0
   subroutine test_atom_pd1_r_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_atom_fd(error, "atom pd1_r", 0)
   end subroutine test_atom_pd1_r_fd

   !> pd2_rr[i,j] = d/d(d_a,j) pd1_r[i]
   subroutine test_atom_pd2_rr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_atom_fd(error, "atom pd2_rr", 1)
   end subroutine test_atom_pd2_rr_fd

   !> pd3_rrr[i,j,k] = d/d(d_a,k) pd2_rr[i,j]
   subroutine test_atom_pd3_rrr_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_atom_fd(error, "atom pd3_rrr", 2)
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
                  call check(error, pd2_a(i, j), pd2_a(j, i), &
                             thr_abs=1.0e-14_wp, thr_rel=1.0e-12_wp)
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
      call run_pair_fd(error, "pair pd1_a", 0, side_b=.false., move_b=.false., sweep=.true.)
   end subroutine test_pair_pd1_a_fd

   !> pd1_b[i] = d/d(d_b,i) pd0 (with d_a fixed)
   subroutine test_pair_pd1_b_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd1_b", 0, side_b=.false., move_b=.true., sweep=.true.)
   end subroutine test_pair_pd1_b_fd

   !> pd2_aa[i,j] = d/d(d_a,j) pd1_a[i]
   subroutine test_pair_pd2_aa_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd2_aa", 1, side_b=.false., move_b=.false., sweep=.false.)
   end subroutine test_pair_pd2_aa_fd

   !> pd2_bb[i,j] = d/d(d_b,j) pd1_b[i]
   subroutine test_pair_pd2_bb_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd2_bb", 1, side_b=.true., move_b=.true., sweep=.false.)
   end subroutine test_pair_pd2_bb_fd

   !> pd2_ab[i,j] = d/d(d_b,j) pd1_a[i] (NOT symmetric in i,j)
   subroutine test_pair_pd2_ab_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd2_ab", 1, side_b=.false., move_b=.true., sweep=.false.)
   end subroutine test_pair_pd2_ab_fd

   !> pd3_aaa[i,j,k] = d/d(d_a,k) pd2_aa[i,j]
   subroutine test_pair_pd3_aaa_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd3_aaa", 2, side_b=.false., move_b=.false., sweep=.false.)
   end subroutine test_pair_pd3_aaa_fd

   !> pd3_bbb[i,j,k] = d/d(d_b,k) pd2_bb[i,j]
   subroutine test_pair_pd3_bbb_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd3_bbb", 2, side_b=.true., move_b=.true., sweep=.false.)
   end subroutine test_pair_pd3_bbb_fd

   !> pd3_aab[i,j,k] = d/d(d_b,k) pd2_aa[i,j] (a-symmetric in i,j)
   subroutine test_pair_pd3_aab_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd3_aab", 2, side_b=.false., move_b=.true., sweep=.false.)
   end subroutine test_pair_pd3_aab_fd

   !> pd3_abb[i,j,k] = d/d(d_a,i) pd2_bb[j,k] (b-symmetric in j,k)
   subroutine test_pair_pd3_abb_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_pair_fd(error, "pair pd3_abb", 2, side_b=.true., move_b=.false., sweep=.false.)
   end subroutine test_pair_pd3_abb_fd

   !> Verify the index-permutation symmetries of pd2_aa/pd2_bb (symmetric)
   !> and pd3_aaa/pd3_aab/pd3_abb/pd3_bbb (fully or partly symmetric)
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
   !> (a,b) symmetry of the underlying K_pair
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
               !* expression, and vice versa
               call check(error, pd2ab_1(i, j), pd2ab_2(j, i), thr_abs=stol, thr_rel=stol)
               if (allocated(error)) return
            end do
         end do
         do k = 1, ndim
            do j = 1, ndim
               do i = 1, ndim
                  call check(error, pd3aaa_1(i, j, k), pd3bbb_2(i, j, k), &
                             thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3bbb_1(i, j, k), pd3aaa_2(i, j, k), &
                             thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  !* pd3_aab(i,j,k) (a-symmetric (i,j), b-singleton k) swaps to
                  !* pd3_abb(k,i,j) (a-singleton k, b-symmetric (i,j))
                  call check(error, pd3aab_1(i, j, k), pd3abb_2(k, i, j), &
                             thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
                  call check(error, pd3abb_1(i, j, k), pd3aab_2(j, k, i), &
                             thr_abs=stol, thr_rel=stol)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_pair_swap_invariance

   !* ================================================================================= *!
   !*                            log_lift FD tests                                      *!
   !* ================================================================================= *!

   !> Reference scalar field PD(r) = sum of two atom-centered Gaussians
   !> Used as a known smooth positive PD to verify the level-set lift against
   !> finite differences of log(PD)
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
            pd2(i, j) = (4.0_wp*alpha1**2*d1(i)*d1(j) &
                         - merge(2.0_wp*alpha1, 0.0_wp, i == j))*e1 &
                        + (4.0_wp*alpha2**2*d2(i)*d2(j) &
                           - merge(2.0_wp*alpha2, 0.0_wp, i == j))*e2
         end do
      end do

      do k = 1, ndim
         do j = 1, ndim
            do i = 1, ndim
               pd3(i, j, k) = ( &
                              -8.0_wp*alpha1**3*d1(i)*d1(j)*d1(k) &
                              + 4.0_wp*alpha1**2*( &
                              merge(d1(k), 0.0_wp, i == j) + &
                              merge(d1(j), 0.0_wp, i == k) + &
                              merge(d1(i), 0.0_wp, j == k)))*e1 &
                              + ( &
                              -8.0_wp*alpha2**3*d2(i)*d2(j)*d2(k) &
                              + 4.0_wp*alpha2**2*( &
                              merge(d2(k), 0.0_wp, i == j) + &
                              merge(d2(j), 0.0_wp, i == k) + &
                              merge(d2(i), 0.0_wp, j == k)))*e2
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

   !> log PD gradient of the lift vs FD of log PD
   subroutine test_log_lift_grad_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_lift_fd(error, "log lift grad", 0)
   end subroutine test_log_lift_grad_fd

   !> log PD Hessian of the lift vs FD of its gradient
   subroutine test_log_lift_hess_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_lift_fd(error, "log lift hess", 1)
   end subroutine test_log_lift_hess_fd

   !> log PD third derivative of the lift vs FD of its Hessian
   subroutine test_log_lift_third_fd(error)
      type(error_type), allocatable, intent(out) :: error
      call run_lift_fd(error, "log lift third", 2)
   end subroutine test_log_lift_third_fd

   !* ================================================================================= *!
   !*                  Every lift result is defined on every return path                *!
   !* ================================================================================= *!

   !> No lift routine may leave an `intent(out)` result undefined
   !>
   !> Each of the five lift routines dispatches on `max_deriv` and its lower
   !> branches assign only the orders they were asked for. Without an
   !> unconditional zeroing prologue the higher-order dummies would come back
   !> undefined -- a value the caller cannot even inspect safely, and one that a
   !> `-finit-real=snan` build cannot catch because these are dummy arguments,
   !> not locals
   !>
   !> This is the standing enforcement of that contract, and it is a real test,
   !> not an assertion of intent: every result is poisoned with a sentinel before
   !> the call, so a routine that stopped zeroing would hand the sentinel straight
   !> back. The orders *at or below* `max_deriv` are additionally required to be
   !> non-zero somewhere, so the check cannot pass by the routine having become a
   !> no-op that zeroes everything
   !>
   !> The `PD <= 0` early return (a fully screened-out point) is exercised too:
   !> it shares the same prologue, so all-zero is the correct answer there for
   !> every order
   subroutine test_lift_defines_every_output(error)
      type(error_type), allocatable, intent(out) :: error

      !> Value written into every result before the call; any of it surviving
      !> means the routine left that result undefined
      real(wp), parameter :: poison = -1.2345678901234e30_wp

      !> Aggregate pseudo-density stack; deliberately generic, not a real system
      real(wp) :: pd0, pd1_r(ndim), pd2_rr(ndim, ndim), pd3_rrr(ndim, ndim, ndim)
      real(wp) :: pd4_rrrr(ndim, ndim, ndim, ndim)
      !> Per-atom, contracted and two-nucleus input families
      real(wp) :: qa0(ndim), qa1_r(ndim, ndim), qa2_rr(ndim, ndim, ndim)
      real(wp) :: qa3_rrr(ndim, ndim, ndim, ndim)
      real(wp) :: qq0(ndim, ndim), qq1_r(ndim, ndim, ndim), qq2_rr(ndim, ndim, ndim, ndim)
      real(wp) :: tg0, tg1_r(ndim), tg2_rr(ndim, ndim), tg3_rrr(ndim, ndim, ndim)
      real(wp) :: hv0(ndim), hv1_r(ndim, ndim), hv2_rr(ndim, ndim, ndim)
      !> Lift results
      real(wp) :: f0, f1_r(ndim), f2_rr(ndim, ndim), f3_rrr(ndim, ndim, ndim)
      real(wp) :: f4_rrrr(ndim, ndim, ndim, ndim)
      real(wp) :: n1(ndim), n2(ndim, ndim), n3(ndim, ndim, ndim), n4(ndim, ndim, ndim, ndim)
      real(wp) :: e2(ndim, ndim), e3(ndim, ndim, ndim), e4(ndim, ndim, ndim, ndim)
      real(wp) :: t0, t1(ndim), t2(ndim, ndim), t3(ndim, ndim, ndim)
      real(wp) :: h1(ndim), h2(ndim, ndim), h3(ndim, ndim, ndim)
      !> Dispatch level and a spare index
      integer :: md, ipass

      call fill_ramp(pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr, 1)
      call fill_ramp(tg0, tg1_r, tg2_rr, tg3_rrr, pd4_rrrr, 2)
      qa0 = [0.7_wp, -0.3_wp, 0.4_wp]
      qa1_r = 0.21_wp
      qa2_rr = -0.13_wp
      qa3_rrr = 0.07_wp
      qq0 = 0.31_wp
      qq1_r = -0.19_wp
      qq2_rr = 0.11_wp
      hv0 = [0.5_wp, 0.2_wp, -0.6_wp]
      hv1_r = -0.23_wp
      hv2_rr = 0.17_wp

      ! Pass 1 runs the normal path; pass 2 drives PD to zero so the early
      ! return is covered by exactly the same assertions
      do ipass = 1, 2
         if (ipass == 2) pd0 = 0.0_wp

         do md = 0, 4
            f0 = poison; f1_r = poison; f2_rr = poison
            f3_rrr = poison; f4_rrrr = poison
            call cfc_spatial_eval(pd0, pd1_r, pd2_rr, pd3_rrr, pd4_rrrr, md, &
                                  f0, f1_r, f2_rr, f3_rrr, f4_rrrr)
            call no_poison(error, [f0], "cfc_spatial_eval f0", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(f1_r, [size(f1_r)]), "cfc_spatial_eval f1_r", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(f2_rr, [size(f2_rr)]), "cfc_spatial_eval f2_rr", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(f3_rrr, [size(f3_rrr)]), &
                           "cfc_spatial_eval f3_rrr", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(f4_rrrr, [size(f4_rrrr)]), &
                           "cfc_spatial_eval f4_rrrr", poison)
            if (allocated(error)) return
            if (ipass == 1 .and. md >= 4) then
               call check(error, any(f4_rrrr /= 0.0_wp), &
                          "cfc_spatial_eval returned an all-zero order-4 tensor at max_deriv 4")
               if (allocated(error)) return
            end if
            if (ipass == 1 .and. md < 4) then
               call check(error, all(f4_rrrr == 0.0_wp), &
                          "cfc_spatial_eval left order 4 non-zero below max_deriv 4")
               if (allocated(error)) return
            end if
         end do

         do md = 0, 3
            n1 = poison; n2 = poison; n3 = poison; n4 = poison
            call cfc_nuclear_eval(pd0, pd1_r, pd2_rr, pd3_rrr, qa0, qa1_r, qa2_rr, qa3_rrr, &
                                  md, n1, n2, n3, n4)
            call no_poison(error, reshape(n1, [size(n1)]), "cfc_nuclear_eval f1_rA", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(n2, [size(n2)]), "cfc_nuclear_eval f2_r_rA", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(n3, [size(n3)]), "cfc_nuclear_eval f3_rr_rA", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(n4, [size(n4)]), "cfc_nuclear_eval f4_rrr_rA", poison)
            if (allocated(error)) return

            t0 = poison; t1 = poison; t2 = poison; t3 = poison
            call cfc_tangent_eval(pd0, pd1_r, pd2_rr, pd3_rrr, tg0, tg1_r, tg2_rr, tg3_rrr, &
                                  md, t0, t1, t2, t3)
            call no_poison(error, [t0], "cfc_tangent_eval tangent_f0", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(t1, [size(t1)]), "cfc_tangent_eval tangent_f1_r", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(t2, [size(t2)]), "cfc_tangent_eval tangent_f2_rr", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(t3, [size(t3)]), &
                           "cfc_tangent_eval tangent_f3_rrr", poison)
            if (allocated(error)) return
         end do

         do md = 0, 2
            e2 = poison; e3 = poison; e4 = poison
            call cfc_hessian_eval(pd0, pd1_r, pd2_rr, qa0, qa1_r, qa2_rr, &
                                  qa0, qa1_r, qa2_rr, qq0, qq1_r, qq2_rr, md, e2, e3, e4)
            call no_poison(error, reshape(e2, [size(e2)]), "cfc_hessian_eval f2_rArB", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(e3, [size(e3)]), "cfc_hessian_eval f3_r_rArB", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(e4, [size(e4)]), "cfc_hessian_eval f4_rr_rArB", poison)
            if (allocated(error)) return

            h1 = poison; h2 = poison; h3 = poison
            call cfc_hvp_eval(pd0, pd1_r, pd2_rr, tg0, tg1_r, tg2_rr, qa0, qa1_r, qa2_rr, &
                              hv0, hv1_r, hv2_rr, md, h1, h2, h3)
            call no_poison(error, reshape(h1, [size(h1)]), "cfc_hvp_eval hvp_f1_rA", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(h2, [size(h2)]), "cfc_hvp_eval hvp_f2_r_rA", poison)
            if (allocated(error)) return
            call no_poison(error, reshape(h3, [size(h3)]), "cfc_hvp_eval hvp_f3_rr_rA", poison)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_lift_defines_every_output

   !> Fill a whole pseudo-density stack with a deterministic non-degenerate ramp
   !>
   !> @param[out] a0    Order-0 slot
   !> @param[out] a1    Order-1 slot
   !> @param[out] a2    Order-2 slot
   !> @param[out] a3    Order-3 slot
   !> @param[out] a4    Order-4 slot
   !> @param[in]  seed  Distinguishes two stacks that must not be identical
   subroutine fill_ramp(a0, a1, a2, a3, a4, seed)
      real(wp), intent(out) :: a0
      real(wp), intent(out) :: a1(ndim)
      real(wp), intent(out) :: a2(ndim, ndim)
      real(wp), intent(out) :: a3(ndim, ndim, ndim)
      real(wp), intent(out) :: a4(ndim, ndim, ndim, ndim)
      integer, intent(in) :: seed

      integer :: i, j, k, l

      a0 = 1.0_wp + 0.5_wp*real(seed, wp)
      do i = 1, ndim
         a1(i) = 0.11_wp*real(i*seed, wp) - 0.2_wp
         do j = 1, ndim
            a2(i, j) = 0.07_wp*real(i + j, wp) - 0.05_wp*real(seed, wp)
            do k = 1, ndim
               a3(i, j, k) = 0.013_wp*real(i + 2*j + 3*k + seed, wp)
               do l = 1, ndim
                  a4(i, j, k, l) = 0.003_wp*real(i*j + k*l + seed, wp)
               end do
            end do
         end do
      end do
   end subroutine fill_ramp

   !> Fail unless every element differs from the poison sentinel
   !>
   !> @param[inout] error testdrive failure
   !> @param[in]    vals  Flattened result array
   !> @param[in]    what  Result name, for the message
   !> @param[in]    p     Poison sentinel
   subroutine no_poison(error, vals, what, p)
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: vals(:)
      character(len=*), intent(in) :: what
      real(wp), intent(in) :: p

      call check(error, .not. any(vals == p), &
                 what//" was left undefined above max_deriv (poison survived)")
   end subroutine no_poison

end module test_cavity_drop_cfc
