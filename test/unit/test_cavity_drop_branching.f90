!> DROP branching tests; contains tests for the branch weights individually,
!> and full DROP nuclear gradients
!>
!> - branch weights are the pure kernel of the branch weight model, one
!>   anchor group at a time, no cavity
!> - [[test_cross_branching]] runs the per-point derivative check of the
!>   gradient suite on the five-carbon cross, whose anchors branch, so the
!>   branch post-pass is exercised
!> - [[test_tetra_branching]] runs the same check on a carbon tetrahedron
!>   whose branches the default gates keep
!> - [[test_branching_xyz_totals]] checks the total area and volume gradients
!>   of the same cross against finite differences
module test_cavity_drop_branching
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use test_helpers, only: get_test_cross, get_test_tetrahedron, get_test_cube, &
                           nudge_structure, get_test_forward_pair, get_test_sibling_gate
   use test_cavity_drop_fd, only: check_drop_row_gradient, check_drop_totals_gradient, drop_fd_param, drop_fd_svdw
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_branching, only: branch_weight_type, branch_soft_fraction
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_cavity_drop_lsf_sses, only: moist_cavity_drop_lsf_sses_type
   use moist_cavity_drop_lsf_sses_param, only: moist_cavity_drop_lsf_sses_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, drop_solve_parameters_type, &
                                           drop_branch_parameters_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   implicit none(type, external)
   private

   public :: collect_cavity_drop_branching

   integer, parameter :: ndim = 3

   !> FD step of the per-point checks, [[check_drop_row_gradient]]
   !>
   !> Measured on the steepest gate row of [[test_cross_branching]]
   !> (`d xi = -1744`, `d w_f = 21`): the O(h^6) truncation of the stencil is
   !> 1.6e-4 (9e-8 relative) at 2.5e-4 and drops below 3e-9 at 1.25e-4, the
   !> factor 64 of h^6; round-off grows as 1/h instead and is covered per point
   !> by the driver's stencil noise floor
   real(wp), parameter :: STEP_SIZE = 1.25E-4_wp

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   integer, parameter :: PROJ_MAXITER = 50
   !> Cheap Newton steps of the soft-level-set (k <= 1.5) FD fixtures: the
   !> library's short cheap stage solves their many duplicate seeds in full
   !> steps, about twice the time
   integer, parameter :: SOFT_CHEAP_STEPS = 50

   !> Branch window for [[test_branching_xyz_totals]] so we have a lot of branching
   real(wp), parameter :: TOTALS_BRANCH_DPDIST = 2.8_wp

   !> Branch weight threshold
   real(wp), parameter :: WEIGHT_THR = 1.0e-14_wp
   !> Objective weight of the branch-weight unit tests, `Phi = (ALPHA/2) pdist^2`
   real(wp), parameter :: ALPHA = 0.5_wp
   !> Lebedev points per atom of the sibling and horizon gradient checks
   integer, parameter :: SWITCH_NUM_LEB = 302
   !> Area that counts as a live contribution of a row in a switch band (bohr^2)
   real(wp), parameter :: LIVE_AREA = 1.0e-3_wp
   !> Out-of-plane shift of the forward pair off its mirror geometry (bohr)
   real(wp), parameter :: FORWARD_DZ = 0.06_wp
   !> FD step of the forward-pair totals (bohr)
   real(wp), parameter :: FORWARD_STEP = 2.5e-4_wp
   !> FD step of the horizon totals on the switching anchors (bohr)
   !>
   !> - their volume gradient is small, so the relative threshold is tight; at
   !>   1e-3 the O(h^4) truncation reaches it (1.07 of it, 17 at 2e-3)
   real(wp), parameter :: HORIZON_STEP = 5.0e-4_wp

contains

   !> Register the branch weight and branching gradient tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_cavity_drop_branching(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("branch_partition_of_unity", test_branch_partition_of_unity), &
                  new_unittest("branch_singletons", test_branch_singletons), &
                  new_unittest("branch_ties", test_branch_ties), &
                  new_unittest("branch_window", test_branch_window), &
                  new_unittest("branch_shift_invariance", test_branch_shift_invariance), &
                  new_unittest("branch_distance_window", test_branch_distance_window), &
                  new_unittest("branch_jacobian", test_branch_jacobian), &
                  new_unittest("branch_closed_no_influence", test_branch_closed_no_influence), &
                  new_unittest("branch_gate_strength", test_branch_gate_strength), &
                  new_unittest("branch_gate_sum", test_branch_gate_sum), &
                  new_unittest("branch_birth", test_branch_birth), &
                  new_unittest("closed_sibling_totals", test_closed_sibling_totals), &
                  new_unittest("reach_gate_totals", test_reach_gate_totals), &
                  new_unittest("search_radius_equivalence", test_search_radius_equivalence), &
                  new_unittest("scan_area_continuity", test_scan_area_continuity), &
                  new_unittest("front_gate_scan_continuity", test_front_gate_scan), &
                  new_unittest("row_cut_area_scaling", test_row_cut_scaling), &
                  new_unittest("horizon_totals", test_horizon_totals), &
                  new_unittest("forward_pair_shares", test_forward_pair_shares), &
                  new_unittest("forward_pair_totals", test_forward_pair_totals), &
                  new_unittest("forward_pair_rows", test_forward_pair_rows), &
                  new_unittest("forward_pair_swap_trapezoid", test_forward_pair_swap), &
                  new_unittest("sibling_gate_totals", test_sibling_gate_totals), &
                  new_unittest("cfc_branching_totals", test_cfc_branching_totals), &
                  new_unittest("sses_branching_totals", test_sses_branching_totals), &
                  new_unittest("focal_margin_alpha_invariance", test_focal_margin_alpha_invariance), &
                  new_unittest("cross_branching", test_cross_branching), &
                  new_unittest("tetra_branching", test_tetra_branching), &
                  new_unittest("branching_xyz_totals", test_branching_xyz_totals) &
                  ]

   end subroutine collect_cavity_drop_branching

   !* ================================================================================= *!
   !*                                   Branch weights                                  *!
   !* ================================================================================= *!

   !> With every point open the weights sum to one for every group size, and a
   !> pair follows the closed form `w_1 = 1/(1 + h)`, `w_2 = h/(1 + h)`
   subroutine test_branch_partition_of_unity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: x_max = 0.8_wp
      real(wp), parameter :: x5(5) = [1.30_wp, 1.05_wp, 1.62_wp, 1.11_wp, 1.44_wp]
      real(wp), parameter :: open5(5) = 1.0_wp
      real(wp) :: phi5(5), w(5), h
      integer :: n

      call bw%init(x_max, ALPHA)
      phi5 = phi_at(x5, x_max)

      do n = 1, 5
         call bw%weights(phi5(1:n), open5(1:n), w(1:n))
         call check(error, sum(w(1:n)), 1.0_wp, thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
         if (allocated(error)) return
         call check(error, all(w(1:n) > 0.0_wp) .and. all(w(1:n) <= 1.0_wp), &
                    "weights inside the window are in (0, 1]")
         if (allocated(error)) return
      end do

      !> Pair at a quarter of the window, against the defining exponentials
      h = ref_pair_factor(0.25_wp)
      call bw%weights(phi_at([2.0_wp, 2.0_wp + 0.25_wp*x_max], x_max), open5(1:2), w(1:2))
      call check(error, w(1), 1.0_wp/(1.0_wp + h), thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
      if (allocated(error)) return
      call check(error, w(2), h/(1.0_wp + h), thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)

   end subroutine test_branch_partition_of_unity

   !> A single branch has its gate as weight, whatever its objective
   subroutine test_branch_singletons(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: phi(4) = [0.0_wp, 5.0_wp, -2.0_wp, 100.0_wp]
      real(wp), parameter :: gate(4) = [1.0_wp, 0.375_wp, 0.0_wp, 1.0e-9_wp]
      real(wp) :: w(1), jac_phi(1, 1), jac_gate(1, 1)
      integer :: i

      call bw%init(1.0_wp, ALPHA)
      do i = 1, size(phi)
         call bw%weights_jac(phi(i:i), gate(i:i), w, jac_phi, jac_gate)
         call check(error, w(1), gate(i), thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
         if (allocated(error)) return
         call check(error, jac_phi(1, 1) == 0.0_wp, "a single branch does not depend on its objective")
         if (allocated(error)) return
         call check(error, jac_gate(1, 1), 1.0_wp, thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
         if (allocated(error)) return
      end do

   end subroutine test_branch_singletons

   !> Coinciding objectives share by the gates, and the weights are flat in the
   !> objectives there
   subroutine test_branch_ties(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: g = 0.3125_wp
      real(wp), parameter :: open4(4) = 1.0_wp
      real(wp) :: phi(4), w(4), jac(4, 4), jac_gate(4, 4), w3(3), w2(2)

      call bw%init(0.3_wp, ALPHA)
      phi = 0.75_wp
      call bw%weights_jac(phi, open4, w, jac, jac_gate)
      call check(error, all(w == 0.25_wp), "equal weights at a four-way tie")
      if (allocated(error)) return
      call check(error, all(jac == 0.0_wp), "vanishing derivative at a four-way tie")
      if (allocated(error)) return

      !> Two-way tie at the minimum next to a third sibling inside the window:
      !> the tied pair stays equal and the weights stay finite and normalized
      call bw%weights([0.5_wp, 0.5_wp, 0.6_wp], open4(1:3), w3)
      call check(error, w3(1) == w3(2) .and. w3(3) < w3(1) .and. w3(3) > 0.0_wp, &
                 "tied pair above the third sibling")
      if (allocated(error)) return
      call check(error, sum(w3), 1.0_wp, thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
      if (allocated(error)) return

      !> An open branch tied with a partly open one: `1/(1 + g)` and `g/(1 + g)`
      call bw%weights([0.5_wp, 0.5_wp], [1.0_wp, g], w2)
      call check(error, w2(1), 1.0_wp/(1.0_wp + g), thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
      if (allocated(error)) return
      call check(error, w2(2), g/(1.0_wp + g), thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)

   end subroutine test_branch_ties

   !> A sibling at or beyond the window of a fully open one has weight zero with
   !> zero derivative, and the other siblings do not see it at all; only the
   !> gate of the open one in front of it could uncover it
   subroutine test_branch_window(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: x_max = 0.3_wp
      real(wp), parameter :: open4(4) = 1.0_wp
      real(wp) :: phi(4), w(4), jac(4, 4), jac_gate(4, 4), w_in(3), jac_in(3, 3)
      real(wp) :: jac_gate_in(3, 3), w_far(4)

      call bw%init(x_max, ALPHA)
      phi = phi_at([1.00_wp, 1.12_wp, 1.21_wp, 1.00_wp + x_max], x_max)
      call bw%weights_jac(phi, open4, w, jac, jac_gate)
      call bw%weights_jac(phi(1:3), open4(1:3), w_in, jac_in, jac_gate_in)

      call check(error, w(4), 0.0_wp, thr_abs=0.0_wp, thr_rel=0.0_wp)
      if (allocated(error)) return
      call check(error, all(jac(4, :) == 0.0_wp) .and. all(jac(:, 4) == 0.0_wp), &
                 "a sibling at the window edge has no derivative coupling")
      if (allocated(error)) return
      call check(error, all(jac_gate(:, 4) == 0.0_wp), "the gate of a sibling at the window edge does not enter")
      if (allocated(error)) return
      call check(error, jac_gate(4, 1) < 0.0_wp .and. all(jac_gate(4, 2:3) == 0.0_wp), &
                 "only the gate of the open branch one window in front holds the far sibling at zero")
      if (allocated(error)) return
      call check(error, all(w(1:3) == w_in), "inside siblings unchanged by one at the edge")
      if (allocated(error)) return
      ! The gate of the open branch in front holds the far one back, so its column differs
      call check(error, all(jac(1:3, 1:3) == jac_in) .and. all(jac_gate(1:3, 2:3) == jac_gate_in(:, 2:3)), &
                 "inside Jacobians unchanged by one at the edge")
      if (allocated(error)) return

      !> Moving the outside sibling further away changes nothing
      phi(4) = 57.0_wp
      call bw%weights(phi, open4, w_far)
      call check(error, all(w_far == w), "weights independent of an outside sibling")
      if (allocated(error)) return

      !> A zero window gives everything to the closest branch
      call bw%init(0.0_wp, ALPHA)
      call bw%weights([1.0_wp, 1.0_wp + 1.0e-9_wp], open4(1:2), w(1:2))
      call check(error, w(1) == 1.0_wp .and. w(2) == 0.0_wp, "zero window is closest-branch only")

   end subroutine test_branch_window

   !> Only differences of the distance coordinates enter
   subroutine test_branch_shift_invariance(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: x_max = 0.5_wp
      real(wp), parameter :: x(3) = [1.25_wp, 1.5_wp, 1.375_wp]
      real(wp), parameter :: gate(3) = [0.75_wp, 1.0_wp, 0.5_wp]
      real(wp) :: w(3), w_shift(3), jac(3, 3), jac_gate(3, 3)

      call bw%init(x_max, ALPHA)
      call bw%weights_jac(phi_at(x, x_max), gate, w, jac, jac_gate)
      call bw%weights(phi_at(x + 16.0_wp, x_max), gate, w_shift)
      call check(error, maxval(abs(w_shift - w)) <= 1.0e-12_wp, "weights invariant to a common shift")
      if (allocated(error)) return
      !> `d phi / d x = ALPHA x` along a common shift of the coordinates
      call check(error, maxval(abs(matmul(jac, ALPHA*x))) <= 1.0e-13_wp, &
                 "a common shift has no first-order effect")

   end subroutine test_branch_shift_invariance

   !> The window is a projection distance, the same close to the anchor and far
   !> from it, and the weights stay smooth for a branch at the anchor itself
   subroutine test_branch_distance_window(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: x_max = 0.2_wp
      real(wp), parameter :: pdist_near(3) = [2.0_wp, 20.0_wp, 200.0_wp]
      real(wp), parameter :: open2(2) = 1.0_wp
      real(wp) :: w(2), w_half(3), jac(2, 2), jac_gate(2, 2)
      integer :: i

      call bw%init(x_max, ALPHA)
      do i = 1, size(pdist_near)
         call bw%weights(phi_of([pdist_near(i), pdist_near(i) + 1.01_wp*x_max]), open2, w)
         call check(error, w(1) == 1.0_wp .and. w(2) == 0.0_wp, &
                    "a branch one window behind an open one has no weight at any distance")
         if (allocated(error)) return
         call bw%weights(phi_of([pdist_near(i), pdist_near(i) + 0.5_wp*x_max]), open2, w)
         w_half(i) = w(2)
      end do
      call check(error, minval(w_half) > 0.1_wp, "a branch half a window behind shares")
      if (allocated(error)) return
      call check(error, maxval(w_half) - minval(w_half) <= 1.0e-2_wp, &
                 "the share at half a window does not depend on the distance")
      if (allocated(error)) return

      !> An anchor on the surface: the window is wider, at most `sqrt(2) x_max`,
      !> and the derivatives stay finite
      call bw%weights_jac(phi_of([0.0_wp, 1.42_wp*x_max]), open2, w, jac, jac_gate)
      call check(error, w(1) == 1.0_wp .and. w(2) == 0.0_wp, "window at the surface is bounded")
      if (allocated(error)) return
      call bw%weights_jac(phi_of([0.0_wp, 0.5_wp*x_max]), open2, w, jac, jac_gate)
      call check(error, w(2) > 0.0_wp .and. all(abs(jac) < 1.0e3_wp) .and. maxval(abs(jac)) > 0.0_wp, &
                 "finite objective derivatives for a branch at the anchor")

   end subroutine test_branch_distance_window

   !> Both Jacobians against central differences in the middle of the window
   !> with partly open gates, where the weights actually move, and with a
   !> branch behind the window of a partly open one; `weights_grad` is their
   !> contraction
   subroutine test_branch_jacobian(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: x_max = 0.4_wp, step = 1.0e-5_wp
      real(wp), parameter :: phi(5) = [0.90_wp, 1.10_wp, 1.02_wp, 1.19_wp, 2.6_wp]
      real(wp), parameter :: gate(5) = [0.6_wp, 0.9_wp, 0.35_wp, 0.8_wp, 0.7_wp]
      real(wp) :: w(5), jac(5, 5), jac_gate(5, 5), wp_(5), wm_(5), shifted(5), fd(5)
      real(wp) :: dphi(2, 5), dgate(2, 5), wg(5), dw(2, 5)
      integer :: l

      call bw%init(x_max, ALPHA)
      call bw%weights_jac(phi, gate, w, jac, jac_gate)
      call check(error, maxval(abs(jac)) > 0.1_wp .and. maxval(abs(jac_gate(1:4, 1:4))) > 0.1_wp, &
                 "the probe sits where the weights respond")
      if (allocated(error)) return
      call check(error, w(5) > 1.0e-3_wp .and. all(jac(:, 5) == 0.0_wp) .and. minval(abs(jac_gate(5, :))) > 1.0e-4_wp, &
                 "the far branch does not keep a weight set by the gates in front of it and not by its distance")
      if (allocated(error)) return

      do l = 1, 5
         shifted = phi
         shifted(l) = phi(l) + step
         call bw%weights(shifted, gate, wp_)
         shifted(l) = phi(l) - step
         call bw%weights(shifted, gate, wm_)
         fd = (wp_ - wm_)/(2.0_wp*step)
         call check(error, maxval(abs(fd - jac(:, l))) <= 1.0e-7_wp, &
                    "objective Jacobian column "//to_string(l)//" disagrees with finite differences")
         if (allocated(error)) return

         shifted = gate
         shifted(l) = gate(l) + step
         call bw%weights(phi, shifted, wp_)
         shifted(l) = gate(l) - step
         call bw%weights(phi, shifted, wm_)
         fd = (wp_ - wm_)/(2.0_wp*step)
         call check(error, maxval(abs(fd - jac_gate(:, l))) <= 1.0e-8_wp, &
                    "gate Jacobian column "//to_string(l)//" disagrees with finite differences")
         if (allocated(error)) return
      end do

      dphi(1, :) = [0.3_wp, -0.2_wp, 0.7_wp, 0.1_wp, 0.4_wp]
      dphi(2, :) = [-1.1_wp, 0.4_wp, 0.0_wp, 0.6_wp, -0.3_wp]
      dgate(1, :) = [0.2_wp, 0.5_wp, -0.3_wp, 0.0_wp, 0.6_wp]
      dgate(2, :) = [-0.4_wp, 0.1_wp, 0.9_wp, -0.7_wp, 0.2_wp]
      call bw%weights_grad(phi, gate, dphi, dgate, wg, dw)
      call check(error, maxval(abs(wg - w)) <= WEIGHT_THR, "weights_grad returns the weights")
      if (allocated(error)) return
      call check(error, maxval(abs(dw - matmul(dphi, transpose(jac)) - matmul(dgate, transpose(jac_gate)))) &
                 <= 1.0e-13_wp, "weights_grad is the contraction of the Jacobians")

   end subroutine test_branch_jacobian

   !> A branch with a closed gate has no weight and takes none, wherever it
   !> lies: closer than its siblings, tied with them, between, farther, or more
   !> than a window in front of all of them; values and Jacobians of the others
   !> are those of the group without it
   subroutine test_branch_closed_no_influence(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: phi(3) = [1.00_wp, 1.10_wp, 1.20_wp]
      real(wp), parameter :: gate(3) = [1.0_wp, 0.7_wp, 0.4_wp]
      !> Closer, tied, between, farther, and far in front of the group
      real(wp), parameter :: phi_ghost(5) = [0.81_wp, 1.0_wp, 1.15_wp, 1.4_wp, 0.2_wp]
      real(wp) :: w(4), jac(4, 4), jac_gate(4, 4), wref(3), jac_ref(3, 3), jac_gate_ref(3, 3)
      integer :: i

      call bw%init(0.4_wp, ALPHA)
      call bw%weights_jac(phi, gate, wref, jac_ref, jac_gate_ref)

      do i = 1, size(phi_ghost)
         call bw%weights_jac([phi, phi_ghost(i)], [gate, 0.0_wp], w, jac, jac_gate)
         call check(error, w(4) == 0.0_wp .and. all(abs(w(1:3) - wref) <= WEIGHT_THR), &
                    "a closed branch has a weight or takes one, case "//to_string(i))
         if (allocated(error)) return
         call check(error, maxval(abs(jac(1:3, 1:3) - jac_ref)) <= 1.0e-12_wp &
                    .and. maxval(abs(jac_gate(1:3, 1:3) - jac_gate_ref)) <= 1.0e-12_wp, &
                    "the siblings do not respond as without the closed branch, case "//to_string(i))
         if (allocated(error)) return
         call check(error, all(jac(:, 4) == 0.0_wp) .and. all(jac(4, :) == 0.0_wp), &
                    "the objective of a closed branch enters, case "//to_string(i))
         if (allocated(error)) return
      end do

   end subroutine test_branch_closed_no_influence

   !> The strength of a branch is its gate: a closed branch in front gives a
   !> sibling its whole gate, a fully open one a window in front none of it, and
   !> between the two the far one keeps the share `1 - g` of its gate
   subroutine test_branch_gate_strength(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: phi(2) = [1.00_wp, 5.00_wp], g_far = 0.8_wp
      real(wp) :: w(2), w_last
      integer :: i

      call bw%init(0.4_wp, ALPHA)
      call bw%weights(phi, [0.0_wp, g_far], w)
      call check(error, w(1) == 0.0_wp .and. abs(w(2) - g_far) <= WEIGHT_THR, &
                 "a closed branch in front does not leave the far one its gate")
      if (allocated(error)) return
      call bw%weights(phi, [1.0_wp, g_far], w)
      call check(error, w(1) == 1.0_wp .and. w(2) == 0.0_wp, "a fully open branch a window in front leaves a weight")
      if (allocated(error)) return

      w_last = g_far
      do i = 1, 9
         call bw%weights(phi, [0.1_wp*real(i, wp), g_far], w)
         call check(error, w(2) < w_last .and. w(2) > 0.0_wp .and. w(1) > 0.0_wp &
                    .and. w(1) <= 0.1_wp*real(i, wp), &
                    "the far weight does not fall as the gate in front opens, step "//to_string(i))
         if (allocated(error)) return
         call check(error, abs(w(1) - 0.1_wp*real(i, wp)) <= WEIGHT_THR &
                    .and. abs(w(2) - g_far*(1.0_wp - 0.1_wp*real(i, wp))) <= WEIGHT_THR, &
                    "the far branch does not keep the share the gate in front leaves, step "//to_string(i))
         if (allocated(error)) return
         w_last = w(2)
      end do

   end subroutine test_branch_gate_strength

   !> The weights never exceed the gates and sum to at most one; they sum to
   !> one when the closest branch is fully open, and nearly closed branches act
   !> alone
   subroutine test_branch_gate_sum(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: phi(3) = [1.00_wp, 1.10_wp, 1.20_wp]
      real(wp), parameter :: t = 1.0e-4_wp
      real(wp) :: w(3), gate(3), x(3)
      integer :: i, j, k

      call bw%init(0.4_wp, ALPHA)

      call bw%weights(phi, [1.0_wp, 0.5_wp, 0.6_wp], w)
      call check(error, sum(w), 1.0_wp, thr_abs=WEIGHT_THR, thr_rel=WEIGHT_THR)
      if (allocated(error)) return

      call bw%weights(phi, [0.2_wp, 0.5_wp, 0.6_wp], w)
      call check(error, sum(w) < 1.0_wp .and. all(w > 0.0_wp), "a partly open closest branch leaves a deficit")
      if (allocated(error)) return

      !> All nearly closed: each weight is its own gate to first order
      call bw%weights(phi, [t, 2.0_wp*t, 3.0_wp*t], w)
      call check(error, maxval(abs(w/[t, 2.0_wp*t, 3.0_wp*t] - 1.0_wp)) <= 10.0_wp*t, &
                 "nearly closed branches do not share")
      if (allocated(error)) return

      !> Bounds on a grid of distances and gates, inside and beyond the window
      do i = 0, 4
         do j = 0, 4
            do k = 0, 4
               x = [1.0_wp, 1.0_wp + 0.1_wp*real(i, wp), 1.0_wp + 0.25_wp*real(j, wp)]
               gate = [0.25_wp*real(k, wp), 1.0_wp - 0.2_wp*real(j, wp), 0.5_wp + 0.1_wp*real(i, wp)]
               call bw%weights(phi_at(x, 0.4_wp), gate, w)
               call check(error, all(w <= gate + WEIGHT_THR) .and. all(w >= 0.0_wp) &
                          .and. sum(w) <= 1.0_wp + WEIGHT_THR, &
                          "weight above its gate or total above one at "//to_string(i)//to_string(j)//to_string(k))
               if (allocated(error)) return
            end do
         end do
      end do

   end subroutine test_branch_gate_sum

   !> A branch born with a vanishing gate enters continuously, inside the
   !> window and a window or more in front of its siblings: at gate zero the
   !> weights are the ones without it, and they move in proportion to the gate
   subroutine test_branch_birth(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(branch_weight_type) :: bw
      real(wp), parameter :: x_max = 0.4_wp
      real(wp), parameter :: x(2) = [1.10_wp, 1.30_wp]
      real(wp), parameter :: x_born(4) = [0.95_wp, 1.20_wp, 1.45_wp, 0.40_wp]
      real(wp), parameter :: gate(2) = [0.9_wp, 1.0_wp]
      real(wp) :: w_ref(2), w(3), w_small(3), w_smaller(3)
      integer :: i

      call bw%init(x_max, ALPHA)
      call bw%weights(phi_at(x, x_max), gate, w_ref)
      do i = 1, size(x_born)
         call bw%weights(phi_at([x, x_born(i)], x_max), [gate, 0.0_wp], w)
         call check(error, maxval(abs(w(1:2) - w_ref)) <= WEIGHT_THR .and. w(3) == 0.0_wp, &
                    "a branch born at zero gate changes nothing")
         if (allocated(error)) return
         call bw%weights(phi_at([x, x_born(i)], x_max), [gate, 1.0e-6_wp], w_small)
         call bw%weights(phi_at([x, x_born(i)], x_max), [gate, 1.0e-7_wp], w_smaller)
         call check(error, maxval(abs(w_small(1:2) - w_ref)) <= 1.0e-5_wp, &
                    "a weak newborn branch moves the weights by more than its gate")
         if (allocated(error)) return
         call check(error, abs(maxval(abs(w_smaller(1:2) - w_ref)) - 0.1_wp*maxval(abs(w_small(1:2) - w_ref))) &
                    <= 1.0e-3_wp*maxval(abs(w_small(1:2) - w_ref)) + WEIGHT_THR, &
                    "the weights do not approach their values without the newborn branch linearly")
         if (allocated(error)) return
      end do

   end subroutine test_branch_birth

   !* ================================================================================= *!
   !*                   Gradients of partly open siblings and switch bands              *!
   !* ================================================================================= *!

   !> Total area and volume gradients against finite differences with many
   !> rows inside a wide band of the signed horizon switch, which closes where
   !> the area scaling vanishes; on the switching anchors
   subroutine test_horizon_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: HOR_OFF = 0.0_wp, HOR_ON = 0.6_wp
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: area_band, c, n_sph(3)
      integer :: i, ia

      call get_test_cube(3.8186186036_wp, mol)
      call build_switch_cavity(cavity, ctx, mol, 4.0_wp, 0.2_wp, HOR_OFF, HOR_ON, error)
      if (allocated(error)) return
      area_band = 0.0_wp
      do i = 1, cavity%ngrid
         ia = cavity%anchor_of_row(i)
         n_sph = cavity%anchor%xyz(:, ia) - mol%xyz(:, cavity%owner(i))
         c = dot_product(cavity%normal0(:, i), n_sph)/norm2(n_sph)
         if (c > HOR_OFF .and. c < HOR_ON) area_band = area_band + cavity%a(i)
      end do
      call check(error, area_band > LIVE_AREA, "fixture has no area inside the horizon band")
      if (allocated(error)) return
      call check_drop_totals_gradient(cavity, mol, error, step=HORIZON_STEP, subset=.true.)

   end subroutine test_horizon_totals

   !> Total area and volume gradients against finite differences on an argon
   !> cube whose anchors find second points with closed gates inside a 1 bohr window
   !>
   !> - a point with a closed gate is no row: every row has an area and a
   !>   weight of its own, and the gradients need no row of the closed points
   !> - gradients on the switching anchors
   subroutine test_closed_sibling_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx

      call get_test_cube(3.2731016603_wp, mol)
      call build_switch_cavity(cavity, ctx, mol, 4.0_wp, 1.0_wp, 0.0_wp, 0.15_wp, error)
      if (allocated(error)) return
      call check(error, cavity%n_extended > 0, "no anchor of the fixture is searched beyond its local point")
      if (allocated(error)) return
      call check(error, all(cavity%a(1:cavity%ngrid) > 0.0_wp) .and. all(cavity%wbranch(1:cavity%ngrid) > 0.0_wp) &
                 .and. all(cavity%w_f0(1:cavity%ngrid) > 0.0_wp), "a row without area, weight or gate")
      if (allocated(error)) return
      call check_drop_totals_gradient(cavity, mol, error, subset=.true.)

   end subroutine test_closed_sibling_totals

   !> Reach gate on a built cavity: no row at or beyond the closed edge, rows
   !> with area inside the band, and total area and volume gradients against
   !> finite differences
   !>
   !> - five-atom fixture at `blend_k = 1.5`, whose points lie up to several
   !>   bohr from their anchors; the reference build with a long reach holds
   !>   rows beyond the closed edge of the short one
   subroutine test_reach_gate_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Reach band of the check (bohr)
      real(wp), parameter :: R_ON = 0.9_wp, R_OFF = 1.5_wp
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity, far
      type(moist_context_type), target :: ctx, ctx_far
      real(wp) :: area_band, area_far
      integer :: i

      call get_test_sibling_gate(mol)
      call build_switch_cavity(far, ctx_far, mol, 1.5_wp, 0.2_wp, 0.0_wp, 0.15_wp, error, &
                               reach_on=8.0_wp, reach_off=10.0_wp)
      if (allocated(error)) return
      area_far = sum(far%a(1:far%ngrid), mask=far%pdist(1:far%ngrid) >= R_OFF)
      call check(error, area_far > LIVE_AREA, "the long-reach build has no area beyond the closed edge")
      if (allocated(error)) return

      call build_switch_cavity(cavity, ctx, mol, 1.5_wp, 0.2_wp, 0.0_wp, 0.15_wp, error, &
                               reach_on=R_ON, reach_off=R_OFF)
      if (allocated(error)) return
      call check(error, all(cavity%pdist(1:cavity%ngrid) < R_OFF), "a row at or beyond the closed edge of the reach")
      if (allocated(error)) return
      call check(error, cavity%pdist_search_max <= R_OFF, "an anchor was searched beyond the reach")
      if (allocated(error)) return
      area_band = 0.0_wp
      do i = 1, cavity%ngrid
         if (cavity%pdist(i) > R_ON) area_band = area_band + cavity%a(i)
      end do
      call check(error, area_band > LIVE_AREA, "fixture has no area inside the reach band")
      if (allocated(error)) return
      call check(error, sum(cavity%a(1:cavity%ngrid)) < sum(far%a(1:far%ngrid)) - area_far + LIVE_AREA, &
                 "the reach gate takes no area off")
      if (allocated(error)) return
      call check_drop_totals_gradient(cavity, mol, error, step=FORWARD_STEP)

   end subroutine test_reach_gate_totals

   !> The adaptive search radius gives the rows of a search of every anchor to
   !> the full reach: same points, weights and areas
   !>
   !> - the adaptive run searches most anchors only up to the window behind a
   !>   fully open point; what lies beyond has no weight and gives none
   !> - fixtures with points behind partly open ones, with forward branches and
   !>   with closed second points
   subroutine test_search_radius_equivalence(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Agreement of points (bohr) and of weights and areas (relative)
      real(wp), parameter :: POS_TOL = 1.0e-9_wp, REL_TOL = 1.0e-9_wp
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: adaptive, full
      type(moist_context_type), target :: ctx, ctx_full
      real(wp) :: blend_k, dpdist
      integer :: ifix
      character(len=:), allocatable :: label

      do ifix = 1, 3
         select case (ifix)
         case (1)
            call get_test_sibling_gate(mol)
            blend_k = 1.5_wp
            dpdist = 0.2_wp
            label = "sibling gate"
         case (2)
            call get_test_forward_pair(FORWARD_DZ, mol)
            blend_k = 2.0_wp
            dpdist = 0.2_wp
            label = "forward pair"
         case default
            call get_test_cube(3.2731016603_wp, mol)
            blend_k = 4.0_wp
            dpdist = 1.0_wp
            label = "cube"
         end select
         call build_switch_cavity(adaptive, ctx, mol, blend_k, dpdist, 0.0_wp, 0.15_wp, error)
         if (allocated(error)) return
         call build_switch_cavity(full, ctx_full, mol, blend_k, dpdist, 0.0_wp, 0.15_wp, error, &
                                  full_reach=.true.)
         if (allocated(error)) return
         call check(error, full%n_extended == full%n_searched .and. adaptive%n_extended < adaptive%n_searched &
                    .and. adaptive%n_extended > 0, label//": the two runs do not differ in their searches")
         if (allocated(error)) return
         call check(error, abs(full%pdist_search_max - full%param%branch%reach_off) <= 1.0e-12_wp, &
                    label//": the full search does not go to the reach")
         if (allocated(error)) return
         call check(error, adaptive%ngrid, full%ngrid, label//": row counts differ")
         if (allocated(error)) return
         call check(error, all(adaptive%anchor_of_row(1:full%ngrid) == full%anchor_of_row(1:full%ngrid)) &
                    .and. all(adaptive%branch(1:full%ngrid) == full%branch(1:full%ngrid)), &
                    label//": rows belong to different anchors or branches")
         if (allocated(error)) return
         call check(error, maxval(abs(adaptive%xyz(:, 1:full%ngrid) - full%xyz(:, 1:full%ngrid))) <= POS_TOL, &
                    label//": points differ by "// &
                    to_string(maxval(abs(adaptive%xyz(:, 1:full%ngrid) - full%xyz(:, 1:full%ngrid)))))
         if (allocated(error)) return
         call check(error, maxval(abs(adaptive%wbranch(1:full%ngrid) - full%wbranch(1:full%ngrid))) <= REL_TOL &
                    .and. maxval(abs(adaptive%a(1:full%ngrid) - full%a(1:full%ngrid)) &
                                 /max(full%a(1:full%ngrid), LIVE_AREA)) <= REL_TOL, &
                    label//": weights or areas differ")
         if (allocated(error)) return
         deallocate (adaptive, full)
      end do

   end subroutine test_search_radius_equivalence

   !> The total area stays continuous along a breathing tetrahedron of four
   !> argon atoms, across the births of points behind partly open closer ones
   !>
   !> - half edge from 2.8 to 3.2 bohr at `blend_k = 1.5`: anchors facing the
   !>   inside gain a second row across the cage, more than a window behind a
   !>   partly open closest row, and the number of branched anchors changes
   !>   several times
   !> - bound: the largest second difference of the area at a step of 0.01
   !>   bohr is below a quarter of the largest one at 0.04 bohr, taken around
   !>   the place of the latter; a jump keeps its size under refinement, a
   !>   smooth curve falls 16-fold
   !> - a reach of 6.5 to 7.5 bohr holds the rows across the cage and keeps
   !>   the searches of the scan short
   !> - no surface nucleates in front of a weighted point here: the void at the
   !>   center opens later (half edge about 3.5), when no row is left behind it
   subroutine test_scan_area_continuity(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NCOARSE = 11, NFINE = 9
      real(wp), parameter :: HALF0 = 2.8_wp, H_COARSE = 0.04_wp, H_FINE = 0.01_wp
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: coarse(NCOARSE), fine(NFINE), d2_coarse, d2_fine, d2
      integer :: nbranched(NCOARSE), i, iworst
      logical :: behind

      call argon_tetrahedron(HALF0, mol)
      call build_switch_cavity(cavity, ctx, mol, 1.5_wp, 0.2_wp, 0.0_wp, 0.15_wp, error, &
                               reach_on=6.5_wp, reach_off=7.5_wp)
      if (allocated(error)) return
      behind = .false.
      do i = 1, NCOARSE
         call argon_tetrahedron(HALF0 + real(i - 1, wp)*H_COARSE, mol)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         coarse(i) = sum(cavity%a(1:cavity%ngrid))
         nbranched(i) = count(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)
         behind = behind .or. weighted_behind_partial(cavity)
      end do
      call check(error, behind .and. maxval(nbranched) - minval(nbranched) >= 24, &
                 "the scan crosses no birth of a row behind a partly open one (branched anchors "// &
                 to_string(minval(nbranched))//" to "//to_string(maxval(nbranched))//")")
      if (allocated(error)) return

      d2_coarse = 0.0_wp
      iworst = 2
      do i = 2, NCOARSE - 1
         d2 = abs(coarse(i + 1) - 2.0_wp*coarse(i) + coarse(i - 1))
         if (d2 > d2_coarse) then
            d2_coarse = d2
            iworst = i
         end if
      end do
      call check(error, d2_coarse > 1.0e-2_wp, "the area does not bend along the scan: "//to_string(d2_coarse))
      if (allocated(error)) return

      ! One coarse interval to either side of the worst node, at the fine step
      do i = 1, NFINE
         call argon_tetrahedron(HALF0 + real(iworst - 2, wp)*H_COARSE + real(i - 1, wp)*H_FINE, mol)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         fine(i) = sum(cavity%a(1:cavity%ngrid))
      end do
      d2_fine = 0.0_wp
      do i = 2, NFINE - 1
         d2_fine = max(d2_fine, abs(fine(i + 1) - 2.0_wp*fine(i) + fine(i - 1)))
      end do
      call check(error, d2_fine <= 0.25_wp*d2_coarse, &
                 "second differences of the area "//to_string(d2_coarse)//" (step 0.04) and "//to_string(d2_fine)// &
                 " (step 0.01) do not fall as those of a continuous curve")

   end subroutine test_scan_area_continuity

   !> Area along a scan in which the gate of an anchor's nearest row opens while
   !> a farther row of the same anchor holds weight
   !>
   !> - Ar4 tetrahedron, SvdW `k = 3`, default switches: as the cube grows, a
   !>   point about 0.2 bohr from the anchor opens through its critical-gradient
   !>   gate in front of one about 3.7 bohr away, which gives its weight up
   !> - the second differences of the area fall at a four times finer step as
   !>   those of a continuous curve do; with a critical-gradient band of width
   !>   0.001 they do not (measured)
   subroutine test_front_gate_scan(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NCOARSE = 11, NFINE = 9
      real(wp), parameter :: HALF0 = 3.06_wp, H_COARSE = 0.016_wp, H_FINE = 0.004_wp
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cavity_error
      type(moist_cavity_drop_parameters_type) :: defaults
      real(wp) :: coarse(NCOARSE), fine(NFINE), d2_coarse, d2_fine, d2
      real(wp), allocatable :: gmin(:), gmax(:)
      integer :: i, iworst

      call argon_tetrahedron(HALF0, mol)
      call build_switch_cavity(cavity, ctx, mol, 3.0_wp, defaults%branch%dpdist_max, defaults%w_0hor_off, &
                               defaults%w_0hor_on, error)
      if (allocated(error)) return
      allocate (gmin(mol%nat*SWITCH_NUM_LEB), source=1.0_wp)
      allocate (gmax(mol%nat*SWITCH_NUM_LEB), source=0.0_wp)
      do i = 1, NCOARSE
         call argon_tetrahedron(HALF0 + real(i - 1, wp)*H_COARSE, mol)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         coarse(i) = sum(cavity%a(1:cavity%ngrid))
         call track_front_gates(cavity, gmin, gmax)
      end do
      call check(error, maxval(gmax - gmin) > 0.5_wp, &
                 "no front row opens by more than half in front of a weighted one: "//to_string(maxval(gmax - gmin)))
      if (allocated(error)) return

      d2_coarse = 0.0_wp
      iworst = 2
      do i = 2, NCOARSE - 1
         d2 = abs(coarse(i + 1) - 2.0_wp*coarse(i) + coarse(i - 1))
         if (d2 > d2_coarse) then
            d2_coarse = d2
            iworst = i
         end if
      end do
      call check(error, d2_coarse > 1.0e-1_wp, "the area does not bend along the scan: "//to_string(d2_coarse))
      if (allocated(error)) return

      ! One coarse interval to either side of the worst node, at the fine step
      do i = 1, NFINE
         call argon_tetrahedron(HALF0 + real(iworst - 2, wp)*H_COARSE + real(i - 1, wp)*H_FINE, mol)
         call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         fine(i) = sum(cavity%a(1:cavity%ngrid))
      end do
      d2_fine = 0.0_wp
      do i = 2, NFINE - 1
         d2_fine = max(d2_fine, abs(fine(i + 1) - 2.0_wp*fine(i) + fine(i - 1)))
      end do
      call check(error, d2_fine <= 0.25_wp*d2_coarse, &
                 "second differences of the area "//to_string(d2_coarse)//" (step 0.016) and "//to_string(d2_fine)// &
                 " (step 0.004) do not fall as those of a continuous curve")
   end subroutine test_front_gate_scan

   !> Range of the gates of partly open nearest rows with a weighted row behind
   !>
   !> - per anchor id, over the builds seen so far; a row behind lies more than
   !>   1.5 branch windows farther and has a branch weight above 1e-3
   !>
   !> @param[in]    cavity Built cavity
   !> @param[in,out] gmin   Smallest such gate per anchor id
   !> @param[in,out] gmax   Largest such gate per anchor id
   subroutine track_front_gates(cavity, gmin, gmax)
      type(cavity_type_drop), intent(in) :: cavity
      real(wp), intent(inout) :: gmin(:), gmax(:)

      integer :: ia, id, r, r0, n, near
      real(wp) :: g
      logical :: behind

      do ia = 1, cavity%anchor%ngrid
         r0 = cavity%offset(ia)
         n = cavity%offset(ia + 1) - r0
         if (n < 2) cycle
         near = r0 + minloc(cavity%pdist(r0 + 1:r0 + n), dim=1)
         g = cavity%w_f0(near)
         if (.not. (g > 0.0_wp .and. g < 1.0_wp)) cycle
         behind = .false.
         do r = r0 + 1, r0 + n
            behind = behind .or. (cavity%pdist(r) > cavity%pdist(near) + 1.5_wp*cavity%param%branch%dpdist_max &
                                  .and. cavity%wbranch(r) > 1.0e-3_wp)
         end do
         if (.not. behind) cycle
         id = cavity%anchor%id(ia)
         gmin(id) = min(gmin(id), g)
         gmax(id) = max(gmax(id), g)
      end do
   end subroutine track_front_gates

   !> A branch below the row cut moves the area by an amount that falls with the
   !> cut: taken out, its own area is below the cut and the shares of its
   !> siblings move by at most the cut times their Jacobian ratio
   !>
   !> - the fade-out window of `front_gate_scan_continuity`, built with the
   !>   default cut and with raised ones; anchors lose a branch to the cut and
   !>   keep another
   !> - the largest area difference to the default cut falls at least threefold
   !>   for a tenfold smaller cut
   subroutine test_row_cut_scaling(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NSCAN = 11, NAT = 4
      real(wp), parameter :: HALF0 = 3.10_wp, H_STEP = 0.008_wp
      real(wp), parameter :: CUTS(2) = [1.0e-4_wp, 1.0e-5_wp]
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cavity_error
      type(moist_cavity_drop_parameters_type) :: defaults
      real(wp) :: area_ref(NSCAN), shift(size(CUTS))
      integer, allocatable :: rows_ref(:, :), rows_cut(:)
      integer :: i, ic, nkept

      allocate (rows_ref(NAT*SWITCH_NUM_LEB, NSCAN), rows_cut(NAT*SWITCH_NUM_LEB))
      shift = 0.0_wp
      nkept = 0
      do ic = 0, size(CUTS)
         call argon_tetrahedron(HALF0, mol)
         call build_switch_cavity(cavity, ctx, mol, 3.0_wp, defaults%branch%dpdist_max, defaults%w_0hor_off, &
                                  defaults%w_0hor_on, error)
         if (allocated(error)) return
         if (ic > 0) cavity%param%finalize%wleb_cut = CUTS(ic)
         do i = 1, NSCAN
            call argon_tetrahedron(HALF0 + real(i - 1, wp)*H_STEP, mol)
            call cavity%update(mol, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            if (ic == 0) then
               area_ref(i) = sum(cavity%a(1:cavity%ngrid))
               call rows_per_anchor(cavity, rows_ref(:, i))
            else
               shift(ic) = max(shift(ic), abs(sum(cavity%a(1:cavity%ngrid)) - area_ref(i)))
               if (ic == 1) then
                  call rows_per_anchor(cavity, rows_cut)
                  nkept = nkept + count(rows_cut >= 1 .and. rows_cut < rows_ref(:, i))
               end if
            end if
         end do
         deallocate (cavity)
      end do
      call check(error, nkept > 0 .and. shift(1) > 1.0e-3_wp, &
                 "the raised cut takes no branch from an anchor that keeps another ("//to_string(nkept)// &
                 ") or does not move the area ("//to_string(shift(1))//")")
      if (allocated(error)) return
      call check(error, shift(2) <= shift(1)/3.0_wp, &
                 "area shift "//to_string(shift(1))//" at cut 1e-4 and "//to_string(shift(2))// &
                 " at 1e-5 does not fall with the cut")

   contains

      !> Rows per anchor id of a build
      subroutine rows_per_anchor(cavity, nrow)
         type(cavity_type_drop), intent(in) :: cavity
         integer, intent(out) :: nrow(:)

         integer :: ia

         nrow = 0
         do ia = 1, cavity%anchor%ngrid
            nrow(cavity%anchor%id(ia)) = cavity%offset(ia + 1) - cavity%offset(ia)
         end do
      end subroutine rows_per_anchor

   end subroutine test_row_cut_scaling

   !* ================================================================================= *!
   !*                        Branching on the other level sets                          *!
   !* ================================================================================= *!

   !> Total area and volume gradients against finite differences on a CFC cavity
   !> with branched anchors and partly open rows under the default switches, on
   !> the switching anchors
   subroutine test_cfc_branching_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(moist_cavity_drop_lsf_cfc_type) :: cfc

      call argon_tetrahedron(3.2_wp, mol)
      call cfc%new()
      call build_default_switch_cavity(cavity, ctx, mol, cfc, error)
      if (allocated(error)) return
      call require_branched_partial(cavity, 24, error)
      if (allocated(error)) return
      call check_drop_totals_gradient(cavity, mol, error, subset=.true.)
   end subroutine test_cfc_branching_totals

   !> Total area and volume gradients against finite differences on an SSES
   !> cavity with a 2.7 bohr probe, branched anchors and partly open rows, on
   !> the switching anchors
   subroutine test_sses_branching_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(moist_cavity_drop_lsf_sses_type) :: sses

      call argon_tetrahedron(3.6_wp, mol)
      call sses%new(moist_cavity_drop_lsf_sses_param_type(probe_radius=2.7_wp))
      call build_default_switch_cavity(cavity, ctx, mol, sses, error)
      if (allocated(error)) return
      call require_branched_partial(cavity, 12, error)
      if (allocated(error)) return
      call check_drop_totals_gradient(cavity, mol, error, step=2.0e-3_wp, subset=.true.)
   end subroutine test_sses_branching_totals

   !> Switch cavity on a given level set with the default switch bands
   !>
   !> @param[out]   cavity Built cavity
   !> @param[in,out] ctx    Run context borrowed by the cavity
   !> @param[in]    mol    Structure
   !> @param[in]    lsf    Level set
   !> @param[out]   error  Error handle
   subroutine build_default_switch_cavity(cavity, ctx, mol, lsf, error)
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      type(moist_context_type), target, intent(inout) :: ctx
      type(structure_type), intent(in) :: mol
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_parameters_type) :: defaults

      call build_switch_cavity(cavity, ctx, mol, 0.0_wp, defaults%branch%dpdist_max, defaults%w_0hor_off, &
                               defaults%w_0hor_on, error, lsf=lsf)
   end subroutine build_default_switch_cavity

   !> Require branched anchors and area on partly open rows
   !>
   !> @param[in]  cavity     Built cavity
   !> @param[in]  min_branch Fewest branched anchors
   !> @param[out] error      Error handle
   subroutine require_branched_partial(cavity, min_branch, error)
      type(cavity_type_drop), intent(in) :: cavity
      integer, intent(in) :: min_branch
      type(error_type), allocatable, intent(out) :: error

      integer :: nbranched
      real(wp) :: area_partial

      nbranched = count(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)
      area_partial = sum(cavity%a(1:cavity%ngrid), mask=cavity%w_f0(1:cavity%ngrid) > 0.0_wp &
                         .and. cavity%w_f0(1:cavity%ngrid) < 1.0_wp)
      call check(error, nbranched >= min_branch .and. area_partial > 1.0_wp, &
                 "fixture has "//to_string(nbranched)//" branched anchors and "//to_string(area_partial)// &
                 " bohr^2 on partly open rows")
   end subroutine require_branched_partial

   !> Four argon atoms on alternate corners of a cube, nudged off its symmetry
   !>
   !> @param[in]  half Half edge of the cube (bohr)
   !> @param[out] mol  Structure
   subroutine argon_tetrahedron(half, mol)
      real(wp), intent(in) :: half
      type(structure_type), intent(out) :: mol

      call new(mol, [18, 18, 18, 18], half*reshape([1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, &
                                                     -1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, -1.0_wp, 1.0_wp], [3, 4]))
      call nudge_structure(mol, 1.0e-3_wp)
   end subroutine argon_tetrahedron

   !* ================================================================================= *!
   !*                 Branching between forward-facing points, defaults                 *!
   !* ================================================================================= *!

   !> The forward pair branches under the default switches: one anchor shares
   !> between two forward-facing points clear of the other spheres, equally on
   !> the mirror-symmetric geometry and unequally, both open, off it
   subroutine test_forward_pair_shares(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: w(2)
      logical :: found

      call get_test_forward_pair(0.0_wp, mol)
      call build_switch_cavity(cavity, ctx, mol, 2.0_wp, 0.2_wp, 0.0_wp, 0.15_wp, error)
      if (allocated(error)) return
      call forward_branch(cavity, mol, w, found)
      call check(error, found, "no anchor shares between two forward-facing points")
      if (allocated(error)) return
      call check(error, w(1) > 0.01_wp .and. abs(w(1) - w(2)) <= 1.0e-10_wp, &
                 "mirror-image points do not share equally: "//to_string(w(1))//", "//to_string(w(2)))
      if (allocated(error)) return

      call get_test_forward_pair(FORWARD_DZ, mol)
      call build_switch_cavity(cavity, ctx, mol, 2.0_wp, 0.2_wp, 0.0_wp, 0.15_wp, error)
      if (allocated(error)) return
      call forward_branch(cavity, mol, w, found)
      call check(error, found, "no anchor shares between two forward-facing points off the mirror geometry")
      if (allocated(error)) return
      call check(error, minval(w) > 0.01_wp .and. abs(w(1) - w(2)) > 0.05_wp, &
                 "off the mirror geometry the shares are not unequal and open: "//to_string(w(1))//", "// &
                 to_string(w(2)))

   end subroutine test_forward_pair_shares

   !> Total area and volume gradients against finite differences on the
   !> forward pair, on and off the mirror geometry
   subroutine test_forward_pair_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: w(2), dz
      logical :: found
      integer :: i

      do i = 1, 2
         dz = merge(0.0_wp, FORWARD_DZ, i == 1)
         call get_test_forward_pair(dz, mol)
         call build_switch_cavity(cavity, ctx, mol, 2.0_wp, 0.2_wp, 0.0_wp, 0.15_wp, error)
         if (allocated(error)) return
         call forward_branch(cavity, mol, w, found)
         call check(error, found, "fixture lost its forward branches")
         if (allocated(error)) return
         ! The shares move fast with the out-of-plane shift: a 1e-3 step leaves
         ! an O(h^4) truncation of 4e-8, 2.5e-4 one of 3e-10
         call check_drop_totals_gradient(cavity, mol, error, step=FORWARD_STEP)
         if (allocated(error)) return
      end do

   end subroutine test_forward_pair_totals

   !> Per-point derivative suite (positions, widths, switching values, areas)
   !> on the forward pair off the mirror geometry, whose branches share unequally
   subroutine test_forward_pair_rows(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_parameters_type) :: param

      call get_test_forward_pair(FORWARD_DZ, mol)
      param = drop_fd_param()
      param%num_leb = SWITCH_NUM_LEB
      param%proj_level = 3
      call check_drop_row_gradient(error, mol, param=param, lsf=drop_fd_svdw(blend_k=2.0_wp, blend_3b=3.0_wp), &
                                   require_branching=.true., step=STEP_SIZE, noise_floor=.true.)

   end subroutine test_forward_pair_rows

   !> Across the mirror geometry the two forward points swap order; the area
   !> stays smooth there: the trapezoid residual of the analytic slope falls as
   !> the cube of the step
   subroutine test_forward_pair_swap(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NSAMPLE = 9
      real(wp), parameter :: DT = 5.0e-3_wp, T0 = -2.0e-2_wp
      !> Fall of the residual per halving at least, cubic: 8, linear: 2
      real(wp), parameter :: MIN_FALL = 4.0_wp
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: area(NSAMPLE), slope(NSAMPLE), res(3), w(2), wfirst(2)
      logical :: found
      integer :: i, istep

      do i = 1, NSAMPLE
         call get_test_forward_pair(T0 + real(i - 1, wp)*DT, mol)
         if (i == 1) then
            call build_switch_cavity(cavity, ctx, mol, 2.0_wp, 0.2_wp, 0.0_wp, 0.15_wp, error)
            if (allocated(error)) return
         else
            call cavity%update(mol, error=cavity_error)
         end if
         if (.not. allocated(cavity_error)) call cavity%get_gradient(cavity_error)
         if (allocated(cavity_error)) then
            call test_failed(error, cavity_error%message)
            return
         end if
         area(i) = sum(cavity%a(1:cavity%ngrid))
         slope(i) = cavity%A_tot1_rA(3, 3)
         if (i == 1 .or. i == NSAMPLE) then
            call forward_branch(cavity, mol, w, found)
            call check(error, found, "fixture lost its forward branches")
            if (allocated(error)) return
            if (i == 1) wfirst = w
         end if
      end do
      call check(error, (wfirst(1) - wfirst(2))*(w(1) - w(2)) < 0.0_wp, "the points do not swap order")
      if (allocated(error)) return

      do istep = 1, 3
         res(istep) = 0.0_wp
         do i = 1, NSAMPLE - 2**(istep - 1), 2**(istep - 1)
            res(istep) = max(res(istep), abs(area(i + 2**(istep - 1)) - area(i) &
                             - 0.5_wp*real(2**(istep - 1), wp)*DT*(slope(i) + slope(i + 2**(istep - 1)))))
         end do
      end do
      call check(error, res(3) >= MIN_FALL*res(2) .and. res(2) >= MIN_FALL*res(1), &
                 "trapezoid residuals "//to_string(res(3))//", "//to_string(res(2))//", "//to_string(res(1))// &
                 " do not fall as the cube of the step")

   end subroutine test_forward_pair_swap

   !> Total area and volume gradients against finite differences on the
   !> five-atom fixture whose anchors hold a partly open point next to a
   !> forward-facing row, under the default switches, on the switching anchors
   subroutine test_sibling_gate_totals(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: live

      call get_test_sibling_gate(mol)
      call build_switch_cavity(cavity, ctx, mol, 1.5_wp, 0.2_wp, 0.0_wp, 0.15_wp, error, &
                               cheap_steps=SOFT_CHEAP_STEPS)
      if (allocated(error)) return
      live = sibling_gate_response(cavity)
      call check(error, live > 1.0e-2_wp, &
                 "fixture has no partly open row that moves a sibling weight ("//to_string(live)//")")
      if (allocated(error)) return
      call check(error, weighted_behind_partial(cavity), &
                 "no anchor with a partly open closest row and a row with area more than a window behind it")
      if (allocated(error)) return
      call check_drop_totals_gradient(cavity, mol, error, step=FORWARD_STEP, subset=.true.)

   end subroutine test_sibling_gate_totals

   !> The gate, the branch weights and the areas do not depend on the objective
   !> weight `phi_alpha`: the focal switch acts on the dimensionless focal margin
   !>
   !> - sibling-gate fixture, whose anchors branch with partly open rows
   !> - the projections are the same points for every weight, so rows are
   !>   compared one to one
   !> - sensitivity: moving the focal band at fixed weight changes the gates, so
   !>   rows sit on the focal ramp and a switch argument that scaled with the
   !>   weight would show
   subroutine test_focal_margin_alpha_invariance(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: alphas(2) = [0.25_wp, 1.0_wp]
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: ref, cav
      type(moist_context_type), target :: ctx_ref, ctx
      integer :: i

      call get_test_sibling_gate(mol)
      call build(ref, ctx_ref, 0.5_wp, 1.0_wp, error)
      if (allocated(error)) return
      call check(error, any(ref%offset(2:ref%anchor%ngrid + 1) - ref%offset(1:ref%anchor%ngrid) > 1), &
                 "fixture does not branch")
      if (allocated(error)) return

      ! Measured differences 2e-14 to 4e-14 at a projection tolerance of 1e-14
      do i = 1, size(alphas)
         call build(cav, ctx, alphas(i), 1.0_wp, error)
         if (allocated(error)) return
         call check(error, cav%ngrid, ref%ngrid, "row count depends on the objective weight")
         if (allocated(error)) return
         call check(error, maxval(abs(cav%xyz - ref%xyz)) <= 1.0e-11_wp, &
                    "projected points depend on the objective weight")
         if (allocated(error)) return
         call check(error, maxval(abs(cav%w_f0 - ref%w_f0)) <= 1.0e-11_wp, &
                    "gate depends on the objective weight "//to_string(alphas(i))//": "// &
                    to_string(maxval(abs(cav%w_f0 - ref%w_f0))))
         if (allocated(error)) return
         call check(error, maxval(abs(cav%wbranch - ref%wbranch)) <= 1.0e-11_wp, &
                    "branch weights depend on the objective weight")
         if (allocated(error)) return
         call check(error, maxval(abs(cav%a - ref%a)) <= 1.0e-11_wp*maxval(ref%a), &
                    "areas depend on the objective weight")
         if (allocated(error)) return
      end do

      ! Sensitivity: the same weight with the focal band moved up by a quarter
      call build(cav, ctx, 0.5_wp, 1.25_wp, error)
      if (allocated(error)) return
      call check(error, abs(sum(cav%a) - sum(ref%a)) > LIVE_AREA, &
                 "fixture has no area on the ramp of the focal switch")

   contains

      !> Cavity of the fixture with objective weight `alpha` and the focal band scaled by `band`
      subroutine build(cavity, ctx_loc, alpha_loc, band, err)
         type(cavity_type_drop), allocatable, intent(out) :: cavity
         type(moist_context_type), target, intent(inout) :: ctx_loc
         real(wp), intent(in) :: alpha_loc, band
         type(error_type), allocatable, intent(out) :: err

         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         type(moist_cavity_drop_parameters_type) :: param
         type(mctc_error), allocatable :: cavity_error

         call new_context(ctx_loc, nthreads=0, verbosity=0)
         param = moist_cavity_drop_parameters_type(num_leb=SWITCH_NUM_LEB, tolerance=PROJ_TOL, &
                                                   solve=drop_solve_parameters_type(max_iter=PROJ_MAXITER), &
                                                   proj_level=3, phi_alpha=alpha_loc)
         param%w_0tra_off = band*param%w_0tra_off
         param%w_0tra_on = band*param%w_0tra_on
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.5_wp, blend_3b=3.0_wp))
         allocate (cavity)
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
                              error=cavity_error, param=param, ctx=ctx_loc)
         if (.not. allocated(cavity_error)) call cavity%update(mol, error=cavity_error)
         if (allocated(cavity_error)) call test_failed(err, cavity_error%message)
      end subroutine build

   end subroutine test_focal_margin_alpha_invariance

   !> Whether an anchor has a partly open closest row and a row with area more
   !> than the branch window behind it
   !>
   !> @param[in] cavity Built cavity
   !> @return    found  Such an anchor exists
   function weighted_behind_partial(cavity) result(found)
      type(cavity_type_drop), intent(in) :: cavity
      logical :: found

      integer :: ia, r, r0, n, near

      found = .false.
      do ia = 1, cavity%anchor%ngrid
         r0 = cavity%offset(ia)
         n = cavity%offset(ia + 1) - r0
         if (n < 2) cycle
         near = r0 + minloc(cavity%pdist(r0 + 1:r0 + n), dim=1)
         if (.not. cavity%w_f0(near) < 1.0_wp) cycle
         do r = r0 + 1, r0 + n
            found = found .or. (cavity%pdist(r) > cavity%pdist(near) + 1.5_wp*cavity%param%branch%dpdist_max &
                                .and. cavity%a(r) > LIVE_AREA)
         end do
      end do
   end function weighted_behind_partial

   !> Branch weights of the anchor with the largest area that shares between
   !> two rows facing forward (cosine above the horizon band) with an open
   !> reach gate, ordered by the height of the point
   !>
   !> @param[in]  cavity Built cavity
   !> @param[in]  mol    Structure of the build
   !> @param[out] w      Branch weights, lower point first
   !> @param[out] found  Whether such an anchor exists
   subroutine forward_branch(cavity, mol, w, found)
      type(cavity_type_drop), intent(in) :: cavity
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: w(2)
      logical, intent(out) :: found

      integer :: ia, r, r0
      real(wp) :: n_sph(3), c, area_best
      logical :: forward

      found = .false.
      w = 0.0_wp
      area_best = 0.0_wp
      do ia = 1, cavity%anchor%ngrid
         r0 = cavity%offset(ia)
         if (cavity%offset(ia + 1) - r0 /= 2) cycle
         forward = .true.
         do r = r0 + 1, r0 + 2
            n_sph = cavity%anchor%xyz(:, ia) - mol%xyz(:, cavity%owner(r))
            c = dot_product(cavity%normal0(:, r), n_sph)/norm2(n_sph)
            forward = forward .and. c > cavity%param%w_0hor_on .and. cavity%wbranch(r) > 0.0_wp &
                      .and. cavity%f_reach%f0(cavity%pdist(r)) == 1.0_wp
         end do
         if (.not. forward .or. sum(cavity%a(r0 + 1:r0 + 2)) <= area_best) cycle
         found = .true.
         area_best = sum(cavity%a(r0 + 1:r0 + 2))
         if (cavity%xyz(3, r0 + 1) <= cavity%xyz(3, r0 + 2)) then
            w = cavity%wbranch(r0 + 1:r0 + 2)
         else
            w = cavity%wbranch([r0 + 2, r0 + 1])
         end if
      end do
   end subroutine forward_branch

   !> Cavity for the switch checks: 302 points, level 3, SvdW unless another level set is given
   !>
   !> @param[out]   cavity   Built cavity
   !> @param[in,out] ctx      Run context borrowed by the cavity
   !> @param[in]    mol      Structure
   !> @param[in]    blend_k  SvdW blending sharpness
   !> @param[in]    dpdist   Branch window (bohr)
   !> @param[in]    hor_off  Horizon off-edge
   !> @param[in]    hor_on   Horizon on-edge
   !> @param[out]   error    Error handle
   !> @param[in]    blend_3b SvdW three-body weight (default 3)
   !> @param[in]    reach_on   Open edge of the reach gate (default: the library's)
   !> @param[in]    reach_off  Closed edge of the reach gate (default: the library's)
   !> @param[in]    full_reach Search every anchor to the reach (default off)
   !> @param[in]    cheap_steps Cheap Newton steps per point (default: the library's)
   !> @param[in]    lsf      Level set in place of SvdW; `blend_k` and `blend_3b` then unused
   subroutine build_switch_cavity(cavity, ctx, mol, blend_k, dpdist, hor_off, hor_on, error, blend_3b, &
                                  reach_on, reach_off, full_reach, cheap_steps, lsf)
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      type(moist_context_type), target, intent(inout) :: ctx
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: blend_k, dpdist, hor_off, hor_on
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: blend_3b
      real(wp), intent(in), optional :: reach_on, reach_off
      logical, intent(in), optional :: full_reach
      integer, intent(in), optional :: cheap_steps
      class(moist_cavity_drop_lsf_type), intent(in), optional :: lsf

      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      type(moist_cavity_drop_parameters_type) :: param
      type(mctc_error), allocatable :: cavity_error
      real(wp) :: gamma3

      gamma3 = 3.0_wp
      if (present(blend_3b)) gamma3 = blend_3b
      call new_context(ctx, nthreads=0, verbosity=0)
      param = moist_cavity_drop_parameters_type(num_leb=SWITCH_NUM_LEB, tolerance=PROJ_TOL, &
                                                solve=drop_solve_parameters_type(max_iter=PROJ_MAXITER), &
                                                proj_level=3, branch=drop_branch_parameters_type(dpdist_max=dpdist))
      param%w_0hor_off = hor_off
      param%w_0hor_on = hor_on
      if (present(reach_on)) param%branch%reach_on = reach_on
      if (present(reach_off)) param%branch%reach_off = reach_off
      if (present(full_reach)) param%branch%full_reach = full_reach
      if (present(cheap_steps)) param%solve%max_cheap_steps = cheap_steps
      allocate (cavity)
      if (present(lsf)) then
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=lsf, &
                              error=cavity_error, param=param, ctx=ctx)
      else
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k, blend_3b=gamma3))
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
                              error=cavity_error, param=param, ctx=ctx)
      end if
      if (.not. allocated(cavity_error)) call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
   end subroutine build_switch_cavity

   !> Largest response of a row weight to the gate of a partly open sibling row
   !>
   !> @param[in] cavity Built cavity
   !> @return    live   `max |d w_m / d g_l|`, `l /= m` rows of one anchor, `0 < g_l < 1`
   function sibling_gate_response(cavity) result(live)
      type(cavity_type_drop), intent(in) :: cavity
      real(wp) :: live

      integer :: ia, r0, n, m, l
      real(wp), allocatable :: w(:), jac_phi(:, :), jac_gate(:, :)

      live = 0.0_wp
      do ia = 1, cavity%anchor%ngrid
         r0 = cavity%offset(ia)
         n = cavity%offset(ia + 1) - r0
         if (n < 2) cycle
         allocate (w(n), jac_phi(n, n), jac_gate(n, n))
         call cavity%branch_weight%weights_jac(cavity%phi0(r0 + 1:r0 + n), cavity%w_f0(r0 + 1:r0 + n), &
                                               w, jac_phi, jac_gate)
         do l = 1, n
            if (.not. (cavity%w_f0(r0 + l) > 0.0_wp .and. cavity%w_f0(r0 + l) < 1.0_wp)) cycle
            do m = 1, n
               if (m /= l) live = max(live, abs(jac_gate(m, l)))
            end do
         end do
         deallocate (w, jac_phi, jac_gate)
      end do
   end function sibling_gate_response

   !* ================================================================================= *!
   !*                               Branch weight helpers                               *!
   !* ================================================================================= *!

   !> Objective of a branch at the projection distance `pdist`
   !>
   !> @param[in] pdist Projection distance
   !> @return    phi   Objective `(ALPHA/2) pdist^2`
   elemental function phi_of(pdist) result(phi)
      !> Projection distance
      real(wp), intent(in) :: pdist
      !> Objective
      real(wp) :: phi

      phi = 0.5_wp*ALPHA*pdist**2
   end function phi_of

   !> Objective of a branch at the softened distance coordinate `x`
   !>
   !> @param[in] x     Distance coordinate `sqrt(pdist^2 + eps^2)`
   !> @param[in] x_max Window, which sets `eps = branch_soft_fraction*x_max`
   !> @return    phi   Objective
   elemental function phi_at(x, x_max) result(phi)
      !> Distance coordinate
      real(wp), intent(in) :: x
      !> Window
      real(wp), intent(in) :: x_max
      !> Objective
      real(wp) :: phi

      phi = 0.5_wp*ALPHA*(x**2 - (branch_soft_fraction*x_max)**2)
   end function phi_at

   !> Pair factor from its defining exponentials, `eta(1-t)/(eta(t) + eta(1-t))`
   !>
   !> @param[in] t Coordinate difference over the window, in (0, 1)
   !> @return    h Pair factor
   pure function ref_pair_factor(t) result(h)
      !> Scaled difference
      real(wp), intent(in) :: t
      !> Pair factor
      real(wp) :: h

      h = exp(-1.0_wp/(1.0_wp - t))/(exp(-1.0_wp/t) + exp(-1.0_wp/(1.0_wp - t)))

   end function ref_pair_factor

   !* ================================================================================= *!
   !*                                Branching gradients                                *!
   !* ================================================================================= *!

   !> Full derivative suite on a geometry whose anchors actually branch
   !>
   !> The five-carbon cross has concave seams where the projection has several
   !> minima per anchor, so the projector returns siblings. Keeping them alive
   !> needs two departures from the production settings:
   !>
   !> - the siblings are backward projections (`c < 0`), which the signed
   !>   horizon switch closes; it is switched off
   !> - a branch window far wider than the production one, inside which no
   !>   second point of this fixture falls
   !>
   !> Exercises the branch post-pass, which adds the branch-weight response to
   !> `f1_rA`, `a_i1_rA` and `v1_rA` after the main loop; `wleb1_rA` and
   !> `xi1_rA` carry no share
   subroutine test_cross_branching(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_parameters_type) :: param

      call get_test_cross(mol)

      param = drop_fd_param()
      param%proj_level = 3
      param%branch%dpdist_max = 7.0_wp
      param%w_0hor_on = 0.0_wp
      param%solve%max_cheap_steps = SOFT_CHEAP_STEPS
      call check_drop_row_gradient(error, mol, param=param, lsf=drop_fd_svdw(blend_k=1.0_wp, blend_3b=1.0_wp), &
                                   require_branching=.true., step=STEP_SIZE, noise_floor=.true.)
   end subroutine test_cross_branching

   !> Full derivative suite on branches the default gates keep
   !>
   !> Four carbons around a hollow at the default SvdW blend: anchors facing
   !> the hollow keep a partly open closest point and an open one on its wall
   !> behind it, both inside the default reach and branch window (36 anchors
   !> branched at 110 points, second points at most 3.3 Bohr out)
   subroutine test_tetra_branching(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_cavity_drop_parameters_type) :: param

      call get_test_tetrahedron(mol)

      param = drop_fd_param()
      param%num_leb = 110
      param%proj_level = 3
      call check_drop_row_gradient(error, mol, param=param, lsf=drop_fd_svdw(blend_k=5.5_wp, blend_3b=3.0_wp), &
                                   require_branching=.true., step=STEP_SIZE)
   end subroutine test_tetra_branching

   !* ================================================================================= *!
   !*                                  Branching totals                                 *!
   !* ================================================================================= *!

   !> Total-area and total-volume gradient test for four carbons in a planar cross plus a central carbon
   !>
   !> Wider branch window, so seam anchors keep their siblings; fails if no anchor branches
   subroutine test_branching_xyz_totals(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, nthreads=0, verbosity=0)

      call get_test_cross(mol)

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=1.0_wp, &
            blend_3b=1.0_wp))
         ! The cross's siblings are backward projections (c < 0), which the
         ! signed horizon switch closes under defaults
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=110, do_fine=.true., &
            tolerance=PROJ_TOL, solve=drop_solve_parameters_type(max_iter=PROJ_MAXITER, &
            max_cheap_steps=SOFT_CHEAP_STEPS), proj_level=3, w_0hor_on=0.0_wp, &
            branch=drop_branch_parameters_type(dpdist_max=TOTALS_BRANCH_DPDIST)), ctx=ctx)
      end block
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
      if (allocated(error)) return

      !> The cavity must *borrow* the caller-owned run context, not copy it
      call check(error, associated(cavity%ctx, ctx), &
                 "DROP cavity does not borrow the caller-owned run context")
      if (allocated(error)) return
      call check(error, cavity%ctx%verbosity() == 0, &
                 "DROP cavity does not observe the context verbosity")
      if (allocated(error)) return

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) call test_failed(error, cavity_error%message)
      if (allocated(error)) return

      if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
         call test_failed(error, "Fixture requested branching but no anchor branched")
         return
      end if

      call check_drop_totals_gradient(cavity, mol, error)
   end subroutine test_branching_xyz_totals

end module test_cavity_drop_branching
