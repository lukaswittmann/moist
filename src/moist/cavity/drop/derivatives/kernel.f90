!> Per-grid point sensitivity kernel shared by the DROP adjoint paths
!>
!> Both DROP reverse-mode paths (electronic and nuclear) differentiate the same
!> per-point map: from a perturbation of the projected point `r`, the multiplier
!> `lambda` and the level-set jet `(S, grad S, grad^2 S)`, through the tangent
!> frame, the closest-point Jacobian `J`, the switching function(s), the Lebedev
!> weight and the Gaussian width, to the surface observables
!>
!> That map is linear in the perturbation, so a caller obtains the derivative
!> along any parameter by seeding it once per basis direction and summing
!> The two stages are:
!>
!>   1. [[build_seed_state]]: everything that depends only on the grid point
!>      evaluated once
!>   2. [[apply_seed]]: the linear response to one seed, evaluated once per
!>      basis direction
!>
!> Degeneracy is reported: [[build_seed_state]] returns a status code and
!> leaves the fatal-versus-skip decision to the caller
!>
!> Every seed of [[apply_seed]] holds the anchor's sphere frame rigid (its
!> normal `n_sph`); a rotation of that frame, the molecular-frame term of a
!> rotated anchor, is [[apply_sphere_rotation]]
!>
module moist_cavity_drop_derivatives_kernel
   use mctc_env_accuracy, only: wp
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric
   use moist_math_packed, only: packed_idx33, packed_contract2, packed_contract3
   use moist_math_smoothstep, only: moist_math_smoothstep_type
   use moist_cavity_drop_branching, only: branch_weight_type

   implicit none(type, external)
   private

   public :: drop_seed_state_type, drop_seed_result_type, drop_surface_weights_type
   public :: build_seed_state, apply_seed, apply_sphere_rotation, compute_branch_adj
   public :: seed_status_message
   public :: seed_state_ok, seed_state_singular_gradient
   public :: seed_state_singular_bmat, seed_state_singular_jacobian, seed_state_singular_kkt
   public :: seed_weight_tol, seed_det_b_guard, seed_curv_disc_guard

   !> Grid point is usable
   integer, parameter :: seed_state_ok = 0
   !> `|grad S|` vanished; the surface normal is undefined
   integer, parameter :: seed_state_singular_gradient = 1
   !> Tangent-restricted KKT matrix `B` is singular
   integer, parameter :: seed_state_singular_bmat = 2
   !> Closest-point Jacobian `J` vanished
   integer, parameter :: seed_state_singular_jacobian = 3
   !> Bordered KKT sensitivity system is singular
   integer, parameter :: seed_state_singular_kkt = 4

   !> Magnitude below which a weight or a norm counts as zero
   real(wp), parameter :: seed_weight_tol = 1.0e-30_wp
   !> Guard on `det(B)` before the tangent-restricted inverse is formed
   real(wp), parameter :: seed_det_b_guard = 1.0e-30_wp
   !> Below this discriminant the `k1`/`k2` split is treated as umbilic; the
   !> individual principal-curvature derivatives are ill-defined at `k1 = k2`,
   !> while the mean and Gaussian curvatures stay smooth
   real(wp), parameter :: seed_curv_disc_guard = 1.0e-10_wp
   !> Below this gap of `B` the switched eigenvalue is treated as degenerate
   real(wp), parameter :: seed_switch_gap_guard = 1.0e-10_wp

   !> Per-grid point forward state consumed by [[apply_seed]]
   !>
   !> Fields are filled in two stages:
   !>
   !> - the `Inputs` block is written by the caller before [[build_seed_state]] runs
   !> - everything below it is derived
   type :: drop_seed_state_type

      !* --------------------------------- Inputs --------------------------------- *!

      !> Level-set gradient and packed Hessian at the projected point
      !>
      !> - packed order of [[moist_math_packed]]
      real(wp) :: lsf1_r(3) = 0.0_wp, lsf2_rr(6) = 0.0_wp
      !> Packed third spatial derivative of the level set at the projected point
      real(wp) :: lsf3_rrr(10) = 0.0_wp
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val = 0.0_wp
      !> Objective coefficient `phi_alpha`
      real(wp) :: alpha_coeff = 0.0_wp
      !> Anchor position and its owner sphere center
      real(wp) :: anchor(3) = 0.0_wp, owner_xyz(3) = 0.0_wp
      !> Grid-level weights entering the `wleb` / `xi` chain
      real(wp) :: anchor_wleb0 = 0.0_wp, cpjac_scal0 = 0.0_wp
      real(wp) :: wleb = 0.0_wp, xi0 = 0.0_wp
      !> Reach gate of the projected point and its partial in that point; its
      !> partial in the anchor is the caller's
      real(wp) :: f_reach_f0 = 1.0_wp, f_reach_dy(3) = 0.0_wp
      !> Whether the caller needs the principal-curvature response
      logical :: want_curvature = .false.

      !* -------------------------------- Derived --------------------------------- *!

      !> Level-set gradient norm and its square
      real(wp) :: g_norm = 0.0_wp, g_norm_sq = 0.0_wp
      !> Tangent-restricted KKT matrix `A = alpha*I - lambda*H`
      real(wp) :: A_mat(3, 3) = 0.0_wp
      !> Outward normal and the surface tangent frame built from it
      real(wp) :: n_surf(3) = 0.0_wp, q1(3) = 0.0_wp, q2(3) = 0.0_wp
      !> `A` applied to the tangent frame
      real(wp) :: Aq1(3) = 0.0_wp, Aq2(3) = 0.0_wp
      !> `B = Q^T A Q` and its inverse
      real(wp) :: B11 = 0.0_wp, B12 = 0.0_wp, B22 = 0.0_wp, det_B = 0.0_wp
      real(wp) :: Binv11 = 0.0_wp, Binv12 = 0.0_wp, Binv22 = 0.0_wp
      !> Eigenvector of `B` for the switched eigenvalue, in the frame `[q1, q2]`
      real(wp) :: v_switch(2) = 0.0_wp
      !> Both eigenvalues of `B` coincide, `v_switch` is arbitrary
      logical :: degenerate_switch = .false.
      !> Sphere tangent frame and its projection into the surface tangent plane
      real(wp) :: t1_vec(3) = 0.0_wp, t2_vec(3) = 0.0_wp
      real(wp) :: tau1(2) = 0.0_wp, tau2(2) = 0.0_wp
      real(wp) :: w1(2) = 0.0_wp, w2(2) = 0.0_wp
      !> Lifted tangent vectors, their cross product and `1/J`
      real(wp) :: y1(3) = 0.0_wp, y2(3) = 0.0_wp
      real(wp) :: cross_vec(3) = 0.0_wp, inv_J = 0.0_wp
      !> Gram-Schmidt data for the tangent-frame derivative
      integer :: min_axis = 1
      real(wp) :: n_dot_q1 = 0.0_wp, proj_surf = 0.0_wp, v_norm_surf = 0.0_wp
      !> Switching-function values and slopes
      real(wp) :: f_crit0 = 0.0_wp, f_crit_dS = 0.0_wp
      real(wp) :: f_foc_f0 = 0.0_wp, f_foc_dS = 0.0_wp
      !> Horizon switch of the signed `n_surf . n_sph`, its value and its slope
      !> along `d n_surf`
      real(wp) :: f_hor_f0 = 1.0_wp, f_hor_dn(3) = 0.0_wp
      !> Unit sphere normal at the anchor
      real(wp) :: n_sph(3) = 0.0_wp
      !> Horizon cosine `n_surf . n_sph`
      real(wp) :: cos_hor = 0.0_wp
      !> Slope of the horizon switch in its cosine
      real(wp) :: f_hor_ds = 0.0_wp
      !> Prune switch of the geometric weight `anchor_wleb0*cpjac`, its value and
      !> its slope in the signed weight
      real(wp) :: wleb_prune_f0 = 1.0_wp, wleb_prune_dS = 0.0_wp
      !> `H q_a`, shape operator `S_ab = q_a^T H q_b/|g|` and its half gap
      !> `disc = hypot((S11-S22)/2, S12)` (only when `want_curvature`)
      real(wp) :: Hq1(3) = 0.0_wp, Hq2(3) = 0.0_wp
      real(wp) :: S11_curv = 0.0_wp, S12_curv = 0.0_wp, S22_curv = 0.0_wp
      real(wp) :: half_diff_curv = 0.0_wp, disc_curv = 0.0_wp

   end type drop_seed_state_type

   !> Linear response of the per-point map to one seed
   type :: drop_seed_result_type
      !> Total sensitivity of the level-set gradient and Hessian
      real(wp) :: dg(3) = 0.0_wp, dH(3, 3) = 0.0_wp
      !> Sensitivity of the outward normal and of `|grad S|`
      real(wp) :: dn_surf(3) = 0.0_wp, d_gnorm = 0.0_wp
      !> Sensitivity of the closest-point Jacobian and of the gate
      !> `f_crit * f_foc * f_hor * f_reach` (times the prune switch); the
      !> reach gate through the projected point only
      real(wp) :: dJ = 0.0_wp, dw_f = 0.0_wp
      !> Sensitivity of the Lebedev weight and the Gaussian width
      real(wp) :: dwleb = 0.0_wp, dxi = 0.0_wp
      !> Sensitivity of the principal curvatures (zero unless `want_curvature`)
      real(wp) :: dk1 = 0.0_wp, dk2 = 0.0_wp
   end type drop_seed_result_type

   !> Surface adjoints reduced to the channels the seed loop actually reads
   !>
   !> The area and integration-weight channels of a
   !> `cavity_adjoint_type` are *derived*: with
   !> `a_i = R_I^2 f_i wleb_i`, `w_i = wleb_i` and
   !> `xi_i = swx/(R_I sqrt(wleb_i))` we have `a = c f/xi^2` and `w = c/xi^2`,
   !> so both fold into the width channel through `d/dxi`. The area channel
   !> folds into the switching channel as well, through
   !> `da/df = R_I^2 wleb_i` -- but only for a parameter that moves `f`
   !>
   !> Holds the result of that folding, so the seed loops read one flat set of
   !> weights and never have to know which channel they came from
   type :: drop_surface_weights_type
      !> Gaussian-width adjoint, with the area and weight channels folded in
      real(wp), allocatable :: w_xi(:)
      !> Switching adjoint; carries the area fold only when it was requested
      real(wp), allocatable :: w_f(:)
      !> Projected-position adjoint
      real(wp), allocatable :: w_xyz(:, :)
      !> Outward-normal adjoint
      real(wp), allocatable :: w_n(:, :)
      !> Principal-curvature adjoints
      real(wp), allocatable :: w_k1(:), w_k2(:)
      !> Branch-objective and gate adjoints of the rows from the branch-weight
      !> reverse pass
      real(wp), allocatable :: branch_phi_adj(:), branch_gate_adj(:)
      !> Whether the normal channel carries anything
      logical :: have_wn = .false.
      !> Whether either curvature channel carries anything
      logical :: have_wk = .false.
   end type drop_surface_weights_type

contains

   !> Evaluate every per-grid point quantity the seed loop reuses
   !>
   !> The caller fills the `Inputs` block of `state` first
   !>
   !> - on a degenerate point `status` is set and the derived fields are left
   !>   incomplete
   !>
   !> @param[inout] state             Seed state; inputs read, derived fields written
   !> @param[in]    f_crit            Critical-gradient switching function
   !> @param[in]    f_foc             Focusing switching function
   !> @param[in]    f_hor             Horizon switching function
   !> @param[in]    f_wleb            Lebedev-weight pruning switching function
   !> @param[in]    use_wleb_prune    Whether Lebedev-weight pruning is active
   !> @param[out]   status            One of the `seed_state_*` codes
   subroutine build_seed_state(state, f_crit, f_foc, f_hor, f_wleb, use_wleb_prune, status)
      !> Seed state
      type(drop_seed_state_type), intent(inout) :: state
      !> Switching functions owned by the cavity
      class(moist_math_smoothstep_type), intent(in) :: f_crit, f_foc, f_hor, f_wleb
      !> Whether Lebedev-weight pruning is active
      logical, intent(in) :: use_wleb_prune
      !> Degeneracy status
      integer, intent(out) :: status

      !> Scratch
      real(wp) :: lambda_switch, beta_max
      real(wp) :: vmin_B(2), vmax_B(2)
      !> Lebedev pruning scratch
      real(wp) :: w_pre_i, f_wleb_s, f_wleb_ds
      !> Horizon scratch
      real(wp) :: n_sph(3), cos_hor, f_hor_ds
      !> Cartesian indices
      integer :: iaxis, jaxis

      status = seed_state_ok

      state%g_norm_sq = dot_product(state%lsf1_r, state%lsf1_r)
      state%g_norm = sqrt(state%g_norm_sq)
      if (state%g_norm <= seed_weight_tol) then
         status = seed_state_singular_gradient
         return
      end if
      call f_crit%eval(state%g_norm, state%f_crit0, state%f_crit_dS)

      ! A = alpha*I - lambda*H
      do jaxis = 1, 3
         do iaxis = 1, 3
            state%A_mat(iaxis, jaxis) = -state%lambda_val*state%lsf2_rr(packed_idx33(iaxis, jaxis))
         end do
      end do
      state%A_mat(1, 1) = state%A_mat(1, 1) + state%alpha_coeff
      state%A_mat(2, 2) = state%A_mat(2, 2) + state%alpha_coeff
      state%A_mat(3, 3) = state%A_mat(3, 3) + state%alpha_coeff

      state%n_surf = state%lsf1_r/state%g_norm

      call setup_tangent_frame(state%n_surf, state%q1, state%q2)
      state%Aq1 = matmul(state%A_mat, state%q1)
      state%Aq2 = matmul(state%A_mat, state%q2)
      state%B11 = dot_product(state%q1, state%Aq1)
      state%B12 = dot_product(state%q1, state%Aq2)
      state%B22 = dot_product(state%q2, state%Aq2)
      state%det_B = state%B11*state%B22 - state%B12*state%B12

      ! Gram-Schmidt data of the frame derivative
      state%min_axis = minloc(abs(state%n_surf), dim=1)
      state%n_dot_q1 = state%n_surf(state%min_axis)
      state%proj_surf = 1.0_wp - state%n_dot_q1**2
      state%v_norm_surf = sqrt(max(state%proj_surf, 1.0e-30_wp))

      call eig_2x2_symmetric(state%B11, state%B12, state%B22, lambda_switch, beta_max, &
                             vmin_B, vmax_B)
      state%v_switch = vmin_B
      state%degenerate_switch = beta_max - lambda_switch <= seed_switch_gap_guard
      ! Focal margin `beta_min/alpha`; the slope is kept in `beta_min`
      call f_foc%eval(lambda_switch/state%alpha_coeff, state%f_foc_f0, state%f_foc_dS)
      state%f_foc_dS = state%f_foc_dS/state%alpha_coeff

      if (abs(state%det_B) <= seed_det_b_guard) then
         status = seed_state_singular_bmat
         return
      end if
      state%Binv11 = state%B22/state%det_B
      state%Binv12 = -state%B12/state%det_B
      state%Binv22 = state%B11/state%det_B

      ! The sphere normal is rigid under every seed, so the cosine moves with n_surf alone
      n_sph = state%anchor - state%owner_xyz
      n_sph = n_sph/norm2(n_sph)
      cos_hor = dot_product(state%n_surf, n_sph)
      call f_hor%eval(cos_hor, state%f_hor_f0, f_hor_ds)
      state%f_hor_dn = f_hor_ds*n_sph
      state%n_sph = n_sph
      state%cos_hor = cos_hor
      state%f_hor_ds = f_hor_ds

      call setup_tangent_frame(state%anchor - state%owner_xyz, state%t1_vec, state%t2_vec)
      state%tau1(1) = dot_product(state%q1, state%t1_vec)
      state%tau1(2) = dot_product(state%q2, state%t1_vec)
      state%tau2(1) = dot_product(state%q1, state%t2_vec)
      state%tau2(2) = dot_product(state%q2, state%t2_vec)
      state%w1(1) = state%Binv11*state%tau1(1) + state%Binv12*state%tau1(2)
      state%w1(2) = state%Binv12*state%tau1(1) + state%Binv22*state%tau1(2)
      state%w2(1) = state%Binv11*state%tau2(1) + state%Binv12*state%tau2(2)
      state%w2(2) = state%Binv12*state%tau2(1) + state%Binv22*state%tau2(2)
      state%y1 = state%alpha_coeff*(state%w1(1)*state%q1 + state%w1(2)*state%q2)
      state%y2 = state%alpha_coeff*(state%w2(1)*state%q1 + state%w2(2)*state%q2)

      state%cross_vec(1) = state%y1(2)*state%y2(3) - state%y1(3)*state%y2(2)
      state%cross_vec(2) = state%y1(3)*state%y2(1) - state%y1(1)*state%y2(3)
      state%cross_vec(3) = state%y1(1)*state%y2(2) - state%y1(2)*state%y2(1)
      block
         real(wp) :: J_val
         J_val = sqrt(dot_product(state%cross_vec, state%cross_vec))
         if (J_val <= seed_weight_tol) then
            status = seed_state_singular_jacobian
            return
         end if
         state%inv_J = 1.0_wp/J_val
      end block

      ! The prune switch acts on |w_pre| and is a factor of the gate
      if (use_wleb_prune) then
         w_pre_i = state%anchor_wleb0*state%cpjac_scal0
         call f_wleb%eval(abs(w_pre_i), f_wleb_s, f_wleb_ds)
         state%wleb_prune_f0 = f_wleb_s
         state%wleb_prune_dS = sign(f_wleb_ds, w_pre_i)
      else
         state%wleb_prune_f0 = 1.0_wp
         state%wleb_prune_dS = 0.0_wp
      end if

      call curvature_state(state)

   contains

      !> Shape operator in the surface frame [q1, q2]; KM and the gap are
      !> rotation invariants, the gap a sum of squares (accurate near umbilics)
      subroutine curvature_state(st)
         type(drop_seed_state_type), intent(inout) :: st

         if (.not. st%want_curvature) return
         call packed_contract2(st%lsf2_rr, st%q1, st%Hq1, 3)
         call packed_contract2(st%lsf2_rr, st%q2, st%Hq2, 3)
         st%S11_curv = dot_product(st%q1, st%Hq1)/st%g_norm
         st%S12_curv = dot_product(st%q1, st%Hq2)/st%g_norm
         st%S22_curv = dot_product(st%q2, st%Hq2)/st%g_norm
         st%half_diff_curv = 0.5_wp*(st%S11_curv - st%S22_curv)
         st%disc_curv = hypot(st%half_diff_curv, st%S12_curv)
      end subroutine curvature_state

   end subroutine build_seed_state

   !> Propagate one seed through the per-grid point map
   !>
   !> The seed is the perturbation of the level-set jet at the *fixed* point
   !> (`dlsf1_r`, `dlsf2_rr`) together with the induced motion of the projected
   !> point and its multiplier (`dr`, `dlambda`), which the caller obtains from
   !> the bordered KKT system
   !>
   !> - a perturbation of the level-set *value* enters only through that system,
   !>   so it has no argument here
   !>
   !> @param[in]  state    Per-grid point forward state from [[build_seed_state]]
   !> @param[in]  dlsf1_r  Seed perturbation of `grad S` at fixed `r`
   !> @param[in]  dlsf2_rr Seed perturbation of `grad^2 S` at fixed `r`
   !> @param[in]  dr       Induced motion of the projected point
   !> @param[in]  dlambda  Induced change of the Lagrange multiplier
   !> @param[out] res      Linear response of the per-point map
   pure subroutine apply_seed(state, dlsf1_r, dlsf2_rr, dr, dlambda, res)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Seed perturbation of the level-set gradient and Hessian
      real(wp), intent(in) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Induced motion of the projected point and its multiplier
      real(wp), intent(in) :: dr(3), dlambda
      !> Linear response
      type(drop_seed_result_type), intent(out) :: res

      !> Sensitivity of the tangent-restricted KKT matrix
      real(wp) :: dA(3, 3)
      !> Tangent-frame derivative scratch
      real(wp) :: v_tmp(3), dq1(3), dq2(3)
      !> `B` and `B^-1` sensitivities
      real(wp) :: dAq1(3), dAq2(3), dB11, dB12, dB22, ddet_B
      real(wp) :: dBinv11, dBinv12, dBinv22
      !> Switched eigenvalue response and the symmetric part of `q1^T dA q2`
      real(wp) :: dlambda_switch, dB12_sym, dfoc
      !> Lifted tangent-vector sensitivities
      real(wp) :: dtau1(2), dtau2(2), dw1(2), dw2(2)
      real(wp) :: dy1(3), dy2(3), dcross(3)
      !> Lebedev-weight chain scratch
      real(wp) :: dw_pre
      !> Shape-operator, mean-curvature and half-gap sensitivities
      real(wp) :: dS11, dS12, dS22, dKM, d_disc
      !> `H dr` and `T dr`, packed (6)
      real(wp) :: Hdr(3), Tdr(6)
      !> Cartesian indices
      integer :: iaxis, jaxis

      call packed_contract2(state%lsf2_rr, dr, Hdr, 3)
      res%dg = dlsf1_r + Hdr
      res%dH = dlsf2_rr
      if (any(dr /= 0.0_wp)) then
         call packed_contract3(state%lsf3_rrr, dr, Tdr, 3)
         do jaxis = 1, 3
            do iaxis = 1, 3
               res%dH(iaxis, jaxis) = res%dH(iaxis, jaxis) + Tdr(packed_idx33(iaxis, jaxis))
            end do
         end do
      end if
      do jaxis = 1, 3
         do iaxis = 1, 3
            dA(iaxis, jaxis) = -dlambda*state%lsf2_rr(packed_idx33(iaxis, jaxis)) &
                               - state%lambda_val*res%dH(iaxis, jaxis)
         end do
      end do

      res%dn_surf = (res%dg - state%n_surf*dot_product(state%n_surf, res%dg))/state%g_norm

      ! dQ/dp: q1 comes from Gram-Schmidt of e_k against n, q2 = n x q1
      v_tmp = -res%dn_surf(state%min_axis)*state%n_surf &
              - state%n_dot_q1*res%dn_surf
      if (state%proj_surf > 1.0e-30_wp) then
         dq1 = (v_tmp - state%q1*dot_product(state%q1, v_tmp))/state%v_norm_surf
      else
         dq1 = 0.0_wp
      end if
      dq2(1) = res%dn_surf(2)*state%q1(3) - res%dn_surf(3)*state%q1(2) &
               + state%n_surf(2)*dq1(3) - state%n_surf(3)*dq1(2)
      dq2(2) = res%dn_surf(3)*state%q1(1) - res%dn_surf(1)*state%q1(3) &
               + state%n_surf(3)*dq1(1) - state%n_surf(1)*dq1(3)
      dq2(3) = res%dn_surf(1)*state%q1(2) - res%dn_surf(2)*state%q1(1) &
               + state%n_surf(1)*dq1(2) - state%n_surf(2)*dq1(1)

      ! dB/dp with B = Q^T A Q, using A-symmetry for the mixed term
      dAq1 = matmul(dA, state%q1)
      dAq2 = matmul(dA, state%q2)
      dB11 = 2.0_wp*dot_product(dq1, state%Aq1) + dot_product(state%q1, dAq1)
      dB12 = dot_product(dq1, state%Aq2) + dot_product(dq2, state%Aq1) &
             + dot_product(dAq1, state%q2)
      dB22 = 2.0_wp*dot_product(dq2, state%Aq2) + dot_product(state%q2, dAq2)

      ddet_B = dB11*state%B22 + state%B11*dB22 - 2.0_wp*state%B12*dB12
      dBinv11 = (dB22*state%det_B - state%B22*ddet_B)/(state%det_B*state%det_B)
      dBinv12 = (-dB12*state%det_B + state%B12*ddet_B)/(state%det_B*state%det_B)
      dBinv22 = (dB11*state%det_B - state%B11*ddet_B)/(state%det_B*state%det_B)

      ! Switched-eigenvalue response, d(lambda) = v^T dB v
      ! - With antisymmetric Q^T dQ and B v = lambda v, the in-plane rotation of the frame
      !   drops out, so this equals u^T d(P A P) u for u = Q v
      ! - The cross term takes the symmetric part of q1^T dA q2, exact for any seed
      dB12_sym = dot_product(dq1, state%Aq2) + dot_product(dq2, state%Aq1) &
                 + 0.5_wp*(dot_product(state%q1, dAq2) + dot_product(dAq1, state%q2))
      if (state%degenerate_switch) then
         dlambda_switch = 0.5_wp*(dB11 + dB22)
      else
         dlambda_switch = state%v_switch(1)*state%v_switch(1)*dB11 &
                          + 2.0_wp*state%v_switch(1)*state%v_switch(2)*dB12_sym &
                          + state%v_switch(2)*state%v_switch(2)*dB22
      end if

      ! The sphere tangent frame is rigid, so dt1 = dt2 = 0
      dtau1(1) = dot_product(dq1, state%t1_vec)
      dtau1(2) = dot_product(dq2, state%t1_vec)
      dtau2(1) = dot_product(dq1, state%t2_vec)
      dtau2(2) = dot_product(dq2, state%t2_vec)

      dw1(1) = dBinv11*state%tau1(1) + state%Binv11*dtau1(1) &
               + dBinv12*state%tau1(2) + state%Binv12*dtau1(2)
      dw1(2) = dBinv12*state%tau1(1) + state%Binv12*dtau1(1) &
               + dBinv22*state%tau1(2) + state%Binv22*dtau1(2)
      dw2(1) = dBinv11*state%tau2(1) + state%Binv11*dtau2(1) &
               + dBinv12*state%tau2(2) + state%Binv12*dtau2(2)
      dw2(2) = dBinv12*state%tau2(1) + state%Binv12*dtau2(1) &
               + dBinv22*state%tau2(2) + state%Binv22*dtau2(2)

      dy1 = state%alpha_coeff*(dw1(1)*state%q1 + state%w1(1)*dq1 &
                               + dw1(2)*state%q2 + state%w1(2)*dq2)
      dy2 = state%alpha_coeff*(dw2(1)*state%q1 + state%w2(1)*dq1 &
                               + dw2(2)*state%q2 + state%w2(2)*dq2)

      dcross(1) = dy1(2)*state%y2(3) - dy1(3)*state%y2(2) &
                  + state%y1(2)*dy2(3) - state%y1(3)*dy2(2)
      dcross(2) = dy1(3)*state%y2(1) - dy1(1)*state%y2(3) &
                  + state%y1(3)*dy2(1) - state%y1(1)*dy2(3)
      dcross(3) = dy1(1)*state%y2(2) - dy1(2)*state%y2(1) &
                  + state%y1(1)*dy2(2) - state%y1(2)*dy2(1)
      res%dJ = dot_product(state%cross_vec, dcross)*state%inv_J

      res%d_gnorm = dot_product(state%n_surf, res%dg)
      ! wleb = anchor wleb * J; the gate, its prune factor included, acts on
      ! the switching value instead
      dw_pre = state%anchor_wleb0*res%dJ
      res%dwleb = dw_pre
      dfoc = state%f_foc_dS*dlambda_switch
      res%dw_f = state%f_reach_f0*(state%wleb_prune_f0*(state%f_hor_f0*(state%f_foc_f0*state%f_crit_dS*res%d_gnorm &
                                                                        + state%f_crit0*dfoc) &
                                                        + state%f_crit0*state%f_foc_f0 &
                                                        *dot_product(state%f_hor_dn, res%dn_surf)) &
                                   + state%f_crit0*state%f_foc_f0*state%f_hor_f0*state%wleb_prune_dS*dw_pre) &
                 + state%wleb_prune_f0*state%f_crit0*state%f_foc_f0*state%f_hor_f0*dot_product(state%f_reach_dy, dr)
      if (state%wleb > seed_weight_tol) then
         res%dxi = -0.5_wp*state%xi0*res%dwleb/state%wleb
      else
         res%dxi = 0.0_wp
      end if

      if (.not. state%want_curvature) return

      associate (q1 => state%q1, q2 => state%q2, dH => res%dH)
         ! d(q_a^T H q_b) = dq_a.H q_b + dq_b.H q_a + q_a^T dH q_b (H symmetric)
         dS11 = (2.0_wp*dot_product(dq1, state%Hq1) + dot_product(q1, matmul(dH, q1))) &
                /state%g_norm - state%S11_curv*res%d_gnorm/state%g_norm
         dS12 = (dot_product(dq1, state%Hq2) + dot_product(dq2, state%Hq1) &
                 + dot_product(q1, matmul(dH, q2)))/state%g_norm &
                - state%S12_curv*res%d_gnorm/state%g_norm
         dS22 = (2.0_wp*dot_product(dq2, state%Hq2) + dot_product(q2, matmul(dH, q2))) &
                /state%g_norm - state%S22_curv*res%d_gnorm/state%g_norm
         dKM = 0.5_wp*(dS11 + dS22)

         ! disc = hypot(half_diff, S12)
         if (state%disc_curv > seed_curv_disc_guard) then
            d_disc = (state%half_diff_curv*0.5_wp*(dS11 - dS22) &
                      + state%S12_curv*dS12)/state%disc_curv
         else
            d_disc = 0.0_wp
         end if
         res%dk1 = dKM + d_disc
         res%dk2 = dKM - d_disc
      end associate

   end subroutine apply_seed

   !> Response of the per-point map to a rotation of the anchor's sphere frame
   !>
   !> The rotation `omega` turns the sphere normal by `omega x n_sph` at a
   !> fixed projected point, multiplier and level set; the motion of the anchor
   !> itself is an anchor seed of [[apply_seed]]
   !>
   !> - `n_sph` enters through `c = n_surf . n_sph` alone: the horizon switch,
   !>   and `J = alpha**2 |c|/|det B|` (`y1 x y2 = alpha**2/det(B) c n_surf`
   !>   for orthonormal right-handed frames), so `dJ = J dc/c`
   !> - the Lebedev weight, the width and the prune switch follow `J`
   !> - sets `dJ`, `dwleb`, `dxi` and `dw_f`; every other response is zero
   !>
   !> @param[in]  state Per-grid point forward state from [[build_seed_state]]
   !> @param[in]  omega Rotation (angular displacement) of the sphere frame
   !> @param[out] res   Linear response of the per-point map
   pure subroutine apply_sphere_rotation(state, omega, res)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Rotation of the sphere frame
      real(wp), intent(in) :: omega(3)
      !> Linear response
      type(drop_seed_result_type), intent(out) :: res

      !> Change of `n_surf . n_sph`, and of the geometric weight before the prune
      real(wp) :: dcos, dw_pre

      ! n_surf . (omega x n_sph) = omega . (n_sph x n_surf)
      dcos = omega(1)*(state%n_sph(2)*state%n_surf(3) - state%n_sph(3)*state%n_surf(2)) &
             + omega(2)*(state%n_sph(3)*state%n_surf(1) - state%n_sph(1)*state%n_surf(3)) &
             + omega(3)*(state%n_sph(1)*state%n_surf(2) - state%n_sph(2)*state%n_surf(1))
      if (abs(state%cos_hor) > seed_weight_tol .and. state%inv_J > 0.0_wp) then
         res%dJ = dcos/(state%cos_hor*state%inv_J)
      end if
      dw_pre = state%anchor_wleb0*res%dJ
      res%dwleb = dw_pre
      res%dw_f = state%f_reach_f0*state%f_crit0*state%f_foc_f0 &
                 *(state%wleb_prune_f0*state%f_hor_ds*dcos + state%f_hor_f0*state%wleb_prune_dS*dw_pre)
      if (state%wleb > seed_weight_tol) then
         res%dxi = -0.5_wp*state%xi0*res%dwleb/state%wleb
      end if
   end subroutine apply_sphere_rotation

   !> Reverse pass over the branch weights
   !>
   !> Within an anchor group the switching value carries the branch weight,
   !> `f_m = f_anchor * p_m`, a function of the objectives and gates of every
   !> row of the anchor
   !>
   !> - the switching channel handles `d(f_anchor)`
   !> - the pass converts the adjoint of `p_m` into the adjoints of the
   !>   objectives and gates through the transposed weight Jacobians; the seed
   !>   loops then couple them to the point motion and the gate response of
   !>   every row
   !>
   !> Groups are the CSR rows of one anchor, the set its weights were formed
   !> over; a single row has `p = g`
   !>
   !> @param[in]  offset          Row offsets per anchor (nmax + 1)
   !> @param[in]  phi0            Branch objective per grid point (ngrid)
   !> @param[in]  gate0           Gate per grid point (ngrid)
   !> @param[in]  f_anchor        Switching value per anchor (nmax)
   !> @param[in]  branch_weight   Branch weight model
   !> @param[in]  w_f             Adjoint of the switching value per grid point (ngrid)
   !> @param[out] branch_phi_adj  Adjoint of the branch objective Phi (ngrid)
   !> @param[out] branch_gate_adj Adjoint of the gate (ngrid)
   ! TODO: This is serial; it sh/could be parallelized
   pure subroutine compute_branch_adj(offset, phi0, gate0, f_anchor, branch_weight, w_f, &
                                      branch_phi_adj, branch_gate_adj)
      !> Row offsets per anchor
      integer, intent(in) :: offset(:)
      !> Branch objective per grid point
      real(wp), intent(in) :: phi0(:)
      !> Gate per grid point
      real(wp), intent(in) :: gate0(:)
      !> Switching value per anchor
      real(wp), intent(in) :: f_anchor(:)
      !> Branch weight model
      type(branch_weight_type), intent(in) :: branch_weight
      !> Adjoint of the switching value
      real(wp), intent(in) :: w_f(:)
      !> Adjoints of the rows
      real(wp), intent(out) :: branch_phi_adj(:), branch_gate_adj(:)

      !> Anchor and group bookkeeping
      integer :: ia, r0, n, m, l, max_group
      !> Weights and their Jacobians of one group
      real(wp), allocatable :: weights(:), jac_phi(:, :), jac_gate(:, :)
      !> Adjoint of one branch weight
      real(wp) :: adj_p

      branch_phi_adj = 0.0_wp
      branch_gate_adj = 0.0_wp
      if (size(branch_phi_adj) <= 0) return

      max_group = maxval(offset(2:) - offset(:size(offset) - 1))
      allocate (weights(max_group))
      allocate (jac_phi(max_group, max_group), jac_gate(max_group, max_group))

      do ia = 1, size(offset) - 1
         r0 = offset(ia)
         n = offset(ia + 1) - r0
         if (n <= 0) cycle
         if (n == 1) then
            branch_gate_adj(r0 + 1) = f_anchor(ia)*w_f(r0 + 1)
            cycle
         end if
         if (all(abs(w_f(r0 + 1:r0 + n)) <= seed_weight_tol)) cycle

         call branch_weight%weights_jac(phi0(r0 + 1:r0 + n), gate0(r0 + 1:r0 + n), &
                                        weights(1:n), jac_phi(1:n, 1:n), jac_gate(1:n, 1:n))
         do m = 1, n
            adj_p = f_anchor(ia)*w_f(r0 + m)
            do l = 1, n
               branch_phi_adj(r0 + l) = branch_phi_adj(r0 + l) + jac_phi(m, l)*adj_p
               branch_gate_adj(r0 + l) = branch_gate_adj(r0 + l) + jac_gate(m, l)*adj_p
            end do
         end do
      end do

   end subroutine compute_branch_adj

   !> Render a degeneracy status as a diagnostic message
   !>
   !> Callers prepend their own context and append the offending grid point
   !>
   !> @param[in] status  One of the `seed_state_*` codes
   !> @returns           Human-readable description of the degeneracy
   pure function seed_status_message(status) result(msg)
      !> Degeneracy status
      integer, intent(in) :: status
      !> Description
      character(len=:), allocatable :: msg

      select case (status)
      case (seed_state_singular_gradient)
         msg = "level set gradient vanishes"
      case (seed_state_singular_bmat)
         msg = "tangent Jacobian matrix B is singular after switching"
      case (seed_state_singular_jacobian)
         msg = "closest-point Jacobian vanishes"
      case (seed_state_singular_kkt)
         msg = "KKT sensitivity solve failed"
      case default
         msg = "unknown sensitivity-kernel failure"
      end select

   end function seed_status_message

end module moist_cavity_drop_derivatives_kernel
