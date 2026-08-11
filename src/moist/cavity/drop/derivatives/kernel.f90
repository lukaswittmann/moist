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
module moist_cavity_drop_derivatives_kernel
   use mctc_env_accuracy, only: wp
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif

   implicit none (type, external)
   private

   public :: drop_seed_state_type, drop_seed_result_type, drop_surface_weights_type
   public :: build_seed_state, apply_seed, compute_branch_phi_adj
   public :: seed_status_message
   public :: seed_state_ok, seed_state_singular_gradient
   public :: seed_state_singular_bmat, seed_state_singular_jacobian
   public :: seed_weight_tol, seed_det_b_guard, seed_curv_disc_guard

   !> Grid point is usable
   integer, parameter :: seed_state_ok = 0
   !> `|grad S|` vanished; the surface normal is undefined
   integer, parameter :: seed_state_singular_gradient = 1
   !> Tangent-restricted KKT matrix `B` is singular
   integer, parameter :: seed_state_singular_bmat = 2
   !> Closest-point Jacobian `J` vanished
   integer, parameter :: seed_state_singular_jacobian = 3

   !> Magnitude below which a weight or a norm counts as zero
   real(wp), parameter :: seed_weight_tol = 1.0e-30_wp
   !> Guard on `det(B)` before the tangent-restricted inverse is formed
   real(wp), parameter :: seed_det_b_guard = 1.0e-30_wp
   !> Below this discriminant the `k1`/`k2` split is treated as umbilic; the
   !> individual principal-curvature derivatives are ill-defined at `k1 = k2`,
   !> while the mean and Gaussian curvatures stay smooth
   real(wp), parameter :: seed_curv_disc_guard = 1.0e-10_wp

   !> Per-grid point forward state consumed by [[apply_seed]]
   !>
   !> Fields are filled in two stages. The `Inputs` block is written by the
   !> caller before [[build_seed_state]] runs; everything below it is derived.
   type :: drop_seed_state_type

      !* --------------------------------- Inputs --------------------------------- *!

      !> Level-set gradient and Hessian at the projected point
      real(wp) :: lsf1_r(3) = 0.0_wp, lsf2_rr(3, 3) = 0.0_wp
      !> Third spatial derivative of the level set at the projected point
      real(wp) :: lsf3_rrr(3, 3, 3) = 0.0_wp
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val = 0.0_wp
      !> Objective coefficient `phi_alpha`
      real(wp) :: alpha_coeff = 0.0_wp
      !> Anchor position and its owner sphere centre
      real(wp) :: anchor(3) = 0.0_wp, owner_xyz(3) = 0.0_wp
      !> Grid-level weights entering the `wleb` / `xi` chain
      real(wp) :: anchor_wleb0 = 0.0_wp, cpjac_scal0 = 0.0_wp, w_f0 = 0.0_wp
      real(wp) :: wbranch = 0.0_wp, wleb = 0.0_wp, xi0 = 0.0_wp
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
      !> Eigenvector of `B` for the switched eigenvalue, lifted to 3D
      real(wp) :: u_switch(3) = 0.0_wp
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
      !> Tangent projector and `A P`
      real(wp) :: P_tan(3, 3) = 0.0_wp, AP_tan(3, 3) = 0.0_wp
      !> Switching-function values and slopes
      real(wp) :: f_crit0 = 0.0_wp, f_crit_dS = 0.0_wp
      real(wp) :: f_foc_f0 = 0.0_wp, f_foc_dS = 0.0_wp
      !> Lebedev-weight pruning chain factor
      real(wp) :: wleb_prune_factor = 1.0_wp
      !> Shape-operator invariants (only when `want_curvature`)
      real(wp) :: Hn(3) = 0.0_wp, Cn(3) = 0.0_wp
      real(wp) :: T_curv = 0.0_wp, D_curv = 0.0_wp
      real(wp) :: KM_curv = 0.0_wp, disc_curv = 0.0_wp

   end type drop_seed_state_type

   !> Linear response of the per-point map to one seed
   type :: drop_seed_result_type
      !> Total sensitivity of the level-set gradient and Hessian
      real(wp) :: dg(3) = 0.0_wp, dH(3, 3) = 0.0_wp
      !> Sensitivity of the outward normal and of `|grad S|`
      real(wp) :: dn_surf(3) = 0.0_wp, d_gnorm = 0.0_wp
      !> Sensitivity of the closest-point Jacobian and the focusing switch
      real(wp) :: dJ = 0.0_wp, dw_f = 0.0_wp
      !> Sensitivity of the Lebedev weight and the Gaussian width
      real(wp) :: dwleb = 0.0_wp, dxi = 0.0_wp
      !> Sensitivity of the principal curvatures (zero unless `want_curvature`)
      real(wp) :: dk1 = 0.0_wp, dk2 = 0.0_wp
   end type drop_seed_result_type

   !> Surface adjoints reduced to the channels the seed loop actually reads
   !>
   !> The area and integration-weight channels of a
   !> `cavity_surface_adjoint_type` are *derived*: with
   !> `a_i = R_I^2 f_i wleb_i`, `w_i = wleb_i` and
   !> `xi_i = swx/(R_I sqrt(wleb_i))` we have `a = c f/xi^2` and `w = c/xi^2`,
   !> so both fold into the width channel through `d/dxi`. The area channel
   !> folds into the switching channel as well, through
   !> `da/df = R_I^2 wleb_i` -- but only for a parameter that moves `f`.
   !>
   !> This type holds the result of that folding, so the seed loops read one
   !> flat set of weights and never have to know which channel they came from.
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
      !> Branch-objective adjoint from the softmax reverse pass
      real(wp), allocatable :: branch_phi_adj(:)
      !> Whether the normal channel carries anything
      logical :: have_wn = .false.
      !> Whether either curvature channel carries anything
      logical :: have_wk = .false.
   end type drop_surface_weights_type

contains

   !> Evaluate every per-grid point quantity the seed loop reuses
   !>
   !> The caller fills the `Inputs` block of `state` first. On a degenerate
   !> point `status` is set and the derived fields are left incomplete.
   !>
   !> @param[inout] state             Seed state; inputs read, derived fields written
   !> @param[in]    f_crit            Critical-gradient switching function
   !> @param[in]    f_foc             Focusing switching function
   !> @param[in]    f_wleb            Lebedev-weight pruning switching function
   !> @param[in]    use_wleb_prune    Whether Lebedev-weight pruning is active
   !> @param[out]   status            One of the `seed_state_*` codes
   subroutine build_seed_state(state, f_crit, f_foc, f_wleb, use_wleb_prune, status)
      !> Seed state
      type(drop_seed_state_type), intent(inout) :: state
      !> Switching functions owned by the cavity
      class(moist_cavity_drop_swif), intent(in) :: f_crit, f_foc, f_wleb
      !> Whether Lebedev-weight pruning is active
      logical, intent(in) :: use_wleb_prune
      !> Degeneracy status
      integer, intent(out) :: status

      !> 2x2 eigen-decomposition scratch
      real(wp) :: tr_B, disc, sqrt_disc, beta1, beta2, lambda_switch
      real(wp) :: vmin_B(2), vmax_B(2)
      !> Lebedev pruning scratch
      real(wp) :: w_pre_i, f_wleb_s, f_wleb_ds
      !> Adjugate of the level-set Hessian
      real(wp) :: adjH(3, 3)
      !> Trace and normal contractions of the Hessian
      real(wp) :: trH, nHn, nCn

      status = seed_state_ok

      state%g_norm_sq = dot_product(state%lsf1_r, state%lsf1_r)
      state%g_norm = sqrt(state%g_norm_sq)
      if (state%g_norm <= seed_weight_tol) then
         status = seed_state_singular_gradient
         return
      end if
      call f_crit%eval(state%g_norm, state%f_crit0, state%f_crit_dS)

      ! A = alpha*I - lambda*H
      state%A_mat = -state%lambda_val*state%lsf2_rr
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
      tr_B = state%B11 + state%B22
      state%det_B = state%B11*state%B22 - state%B12*state%B12
      if (abs(state%det_B) <= seed_det_b_guard) then
         status = seed_state_singular_bmat
         return
      end if
      disc = max(0.25_wp*tr_B*tr_B - state%det_B, 0.0_wp)
      sqrt_disc = sqrt(disc)
      beta1 = 0.5_wp*tr_B + sqrt_disc
      beta2 = 0.5_wp*tr_B - sqrt_disc
      call eig_2x2_symmetric(state%B11, state%B12, state%B22, lambda_switch, beta1, &
                             vmin_B, vmax_B)
      state%u_switch = vmin_B(1)*state%q1 + vmin_B(2)*state%q2
      lambda_switch = beta2
      call f_foc%eval(lambda_switch, state%f_foc_f0, state%f_foc_dS)

      state%Binv11 = state%B22/state%det_B
      state%Binv12 = -state%B12/state%det_B
      state%Binv22 = state%B11/state%det_B

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

      state%min_axis = minloc(abs(state%n_surf), dim=1)
      state%n_dot_q1 = state%n_surf(state%min_axis)
      state%proj_surf = 1.0_wp - state%n_dot_q1**2
      state%v_norm_surf = sqrt(max(state%proj_surf, 1.0e-30_wp))

      state%P_tan(:, :) = -spread(state%n_surf, dim=2, ncopies=3) &
                          *spread(state%n_surf, dim=1, ncopies=3)
      state%P_tan(1, 1) = state%P_tan(1, 1) + 1.0_wp
      state%P_tan(2, 2) = state%P_tan(2, 2) + 1.0_wp
      state%P_tan(3, 3) = state%P_tan(3, 3) + 1.0_wp
      state%AP_tan = matmul(state%A_mat, state%P_tan)

      ! d(w_pre * S)/dp = (S + |w_pre|*S') * d(w_pre)/dp
      if (use_wleb_prune) then
         w_pre_i = state%anchor_wleb0*state%cpjac_scal0*state%w_f0
         call f_wleb%eval(abs(w_pre_i), f_wleb_s, f_wleb_ds)
         state%wleb_prune_factor = f_wleb_s + abs(w_pre_i)*f_wleb_ds
      else
         state%wleb_prune_factor = 1.0_wp
      end if

      ! Frame-independent shape-operator invariants of the level set:
      !   T = k1 + k2 = (tr H - n^T H n)/|g|
      !   D = k1 * k2 = (n^T adj(H) n)/|g|^2     (Goldman, implicit surfaces)
      !   k1,k2 = T/2 +/- sqrt((T/2)^2 - D)
      if (state%want_curvature) then
         associate (H => state%lsf2_rr, n => state%n_surf)
            state%Hn = matmul(H, n)
            trH = H(1, 1) + H(2, 2) + H(3, 3)
            nHn = dot_product(n, state%Hn)
            state%T_curv = (trH - nHn)/state%g_norm

            adjH(1, 1) = H(2, 2)*H(3, 3) - H(2, 3)*H(2, 3)
            adjH(2, 2) = H(1, 1)*H(3, 3) - H(1, 3)*H(1, 3)
            adjH(3, 3) = H(1, 1)*H(2, 2) - H(1, 2)*H(1, 2)
            adjH(1, 2) = H(1, 3)*H(2, 3) - H(1, 2)*H(3, 3)
            adjH(1, 3) = H(1, 2)*H(2, 3) - H(2, 2)*H(1, 3)
            adjH(2, 3) = H(1, 2)*H(1, 3) - H(1, 1)*H(2, 3)
            adjH(2, 1) = adjH(1, 2)
            adjH(3, 1) = adjH(1, 3)
            adjH(3, 2) = adjH(2, 3)

            state%Cn = matmul(adjH, n)
            nCn = dot_product(n, state%Cn)
            state%D_curv = nCn/state%g_norm_sq
            state%KM_curv = 0.5_wp*state%T_curv
            state%disc_curv = sqrt(max(state%KM_curv*state%KM_curv - state%D_curv, 0.0_wp))
         end associate
      end if

   end subroutine build_seed_state

   !> Propagate one seed through the per-grid point map
   !>
   !> The seed is the perturbation of the level-set jet at the *fixed* point
   !> (`dlsf1_r`, `dlsf2_rr`) together with the induced motion of the projected
   !> point and its multiplier (`dr`, `dlambda`), which the caller obtains from
   !> the bordered KKT system. A perturbation of the level-set *value* enters
   !> only through that system, so it has no argument here.
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
      !> Tangent projector sensitivities and the switched eigenvalue response
      real(wp) :: dP_tan(3, 3), dM_tan(3, 3), dlambda_switch
      !> Lifted tangent-vector sensitivities
      real(wp) :: dtau1(2), dtau2(2), dw1(2), dw2(2)
      real(wp) :: dy1(3), dy2(3), dcross(3)
      !> Lebedev-weight chain scratch
      real(wp) :: dw_pre
      !> Curvature-invariant sensitivities
      real(wp) :: dadjH(3, 3), dtrH, dnHn, dT, dnCn, dD, d_disc
      !> Cartesian index
      integer :: kaxis

      res%dg = dlsf1_r + matmul(state%lsf2_rr, dr)
      res%dH = dlsf2_rr
      do kaxis = 1, 3
         res%dH(:, :) = res%dH(:, :) + state%lsf3_rrr(:, :, kaxis)*dr(kaxis)
      end do
      dA = -dlambda*state%lsf2_rr - state%lambda_val*res%dH

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

      ! Basis-invariant switched-eigenvalue response via M = P A P
      dP_tan(:, :) = -(spread(res%dn_surf, dim=2, ncopies=3) &
                       *spread(state%n_surf, dim=1, ncopies=3) &
                       + spread(state%n_surf, dim=2, ncopies=3) &
                       *spread(res%dn_surf, dim=1, ncopies=3))
      dM_tan = matmul(dP_tan, state%AP_tan) &
               + matmul(state%P_tan, matmul(dA, state%P_tan)) &
               + matmul(state%P_tan, matmul(state%A_mat, dP_tan))
      dlambda_switch = dot_product(state%u_switch, matmul(dM_tan, state%u_switch))

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
      res%dw_f = state%f_foc_f0*state%f_crit_dS*res%d_gnorm &
                 + state%f_crit0*state%f_foc_dS*dlambda_switch

      dw_pre = state%anchor_wleb0*state%w_f0*res%dJ &
               + state%anchor_wleb0*state%cpjac_scal0*res%dw_f
      res%dwleb = state%wbranch*state%wleb_prune_factor*dw_pre

      if (state%wleb > seed_weight_tol) then
         res%dxi = -0.5_wp*state%xi0*res%dwleb/state%wleb
      else
         res%dxi = 0.0_wp
      end if

      if (.not. state%want_curvature) return

      associate (H => state%lsf2_rr, n => state%n_surf, dH => res%dH)
         ! d(n^T H n) = 2 (dn . H n) + n^T dH n, H symmetric
         dtrH = dH(1, 1) + dH(2, 2) + dH(3, 3)
         dnHn = 2.0_wp*dot_product(res%dn_surf, state%Hn) &
                + dot_product(n, matmul(dH, n))
         dT = (dtrH - dnHn)/state%g_norm - state%T_curv*res%d_gnorm/state%g_norm

         dadjH(1, 1) = dH(2, 2)*H(3, 3) + H(2, 2)*dH(3, 3) - 2.0_wp*H(2, 3)*dH(2, 3)
         dadjH(2, 2) = dH(1, 1)*H(3, 3) + H(1, 1)*dH(3, 3) - 2.0_wp*H(1, 3)*dH(1, 3)
         dadjH(3, 3) = dH(1, 1)*H(2, 2) + H(1, 1)*dH(2, 2) - 2.0_wp*H(1, 2)*dH(1, 2)
         dadjH(1, 2) = dH(1, 3)*H(2, 3) + H(1, 3)*dH(2, 3) &
                       - dH(1, 2)*H(3, 3) - H(1, 2)*dH(3, 3)
         dadjH(1, 3) = dH(1, 2)*H(2, 3) + H(1, 2)*dH(2, 3) &
                       - dH(2, 2)*H(1, 3) - H(2, 2)*dH(1, 3)
         dadjH(2, 3) = dH(1, 2)*H(1, 3) + H(1, 2)*dH(1, 3) &
                       - dH(1, 1)*H(2, 3) - H(1, 1)*dH(2, 3)
         dadjH(2, 1) = dadjH(1, 2)
         dadjH(3, 1) = dadjH(1, 3)
         dadjH(3, 2) = dadjH(2, 3)

         dnCn = 2.0_wp*dot_product(res%dn_surf, state%Cn) &
                + dot_product(n, matmul(dadjH, n))
         dD = dnCn/state%g_norm_sq - 2.0_wp*state%D_curv*res%d_gnorm/state%g_norm

         ! disc^2 = KM^2 - KG, so 2 disc d(disc) = (T/2) dT - dD
         if (state%disc_curv > seed_curv_disc_guard) then
            d_disc = (state%KM_curv*dT - dD)/(2.0_wp*state%disc_curv)
         else
            d_disc = 0.0_wp
         end if
         res%dk1 = 0.5_wp*dT + d_disc
         res%dk2 = 0.5_wp*dT - d_disc
      end associate

   end subroutine apply_seed

   !> Reverse pass over the branch-weight softmax
   !>
   !> Within an anchor group the Lebedev weight carries a softmax factor,
   !> `wleb_m = base_m * p_m`. The per-point seed loop handles `d(base_m)`;
   !> this pass converts the remaining width-induced adjoint `dL/dp_m` into
   !> `dL/dPhi_m`, which the seed loop then couples to the point motion.
   !>
   !> Groups are runs of equal `anchor_id`; points with `branch_count <= 1`
   !> carry no softmax factor and stay at zero.
   !>
   !> @param[in]  branch_count    Number of branches per grid point (ngrid)
   !> @param[in]  anchor_id       Anchor group id per grid point (ngrid)
   !> @param[in]  wbranch         Softmax branch weight per grid point (ngrid)
   !> @param[in]  wleb            Lebedev weight per grid point (ngrid)
   !> @param[in]  xi0             Gaussian width per grid point (ngrid)
   !> @param[in]  sigma_phi       Softmax temperature
   !> @param[in]  w_xi            Effective Gaussian-width adjoint (ngrid)
   !> @param[out] branch_phi_adj  Adjoint of the branch objective Phi (ngrid)
   pure subroutine compute_branch_phi_adj(branch_count, anchor_id, wbranch, wleb, xi0, &
                                          sigma_phi, w_xi, branch_phi_adj)
      !> Branch bookkeeping per grid point
      integer, intent(in) :: branch_count(:), anchor_id(:)
      !> Branch weight, Lebedev weight and Gaussian width per grid point
      real(wp), intent(in) :: wbranch(:), wleb(:), xi0(:)
      !> Softmax temperature
      real(wp), intent(in) :: sigma_phi
      !> Effective Gaussian-width adjoint
      real(wp), intent(in) :: w_xi(:)
      !> Adjoint of the branch objective
      real(wp), intent(out) :: branch_phi_adj(:)

      !> Grid extent and group bookkeeping
      integer :: ngrid, igroup_start, igroup_end, group_size, m_branch, im_grid
      !> Weight-adjoint scratch
      real(wp) :: adj_wleb, adj_branch, mean_adj_branch, factor_m

      branch_phi_adj = 0.0_wp
      ngrid = size(branch_count)
      if (ngrid <= 0) return
      if (.not. any(branch_count > 1)) return
      if (sigma_phi <= seed_weight_tol) return

      igroup_start = 1
      do while (igroup_start <= ngrid)
         if (branch_count(igroup_start) <= 1) then
            igroup_start = igroup_start + 1
            cycle
         end if

         ! Extend the group while anchor_id stays the same
         igroup_end = igroup_start
         do while (igroup_end < ngrid)
            if (anchor_id(igroup_end + 1) /= anchor_id(igroup_start)) exit
            igroup_end = igroup_end + 1
         end do
         group_size = igroup_end - igroup_start + 1

         mean_adj_branch = 0.0_wp
         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            adj_branch = 0.0_wp
            if (abs(w_xi(im_grid)) > seed_weight_tol &
                .and. wleb(im_grid) > seed_weight_tol &
                .and. wbranch(im_grid) > tiny(1.0_wp)) then
               adj_wleb = -0.5_wp*w_xi(im_grid)*xi0(im_grid)/wleb(im_grid)
               factor_m = wleb(im_grid)/wbranch(im_grid)
               adj_branch = adj_wleb*factor_m
            end if
            mean_adj_branch = mean_adj_branch + wbranch(im_grid)*adj_branch
            branch_phi_adj(im_grid) = adj_branch
         end do

         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            branch_phi_adj(im_grid) = -wbranch(im_grid) &
                                      *(branch_phi_adj(im_grid) - mean_adj_branch)/sigma_phi
         end do

         igroup_start = igroup_end + 1
      end do

   end subroutine compute_branch_phi_adj

   !> Render a degeneracy status as a diagnostic message
   !>
   !> Callers prepend their own context and append the offending grid point.
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
      case default
         msg = "unknown sensitivity-kernel failure"
      end select

   end function seed_status_message

end module moist_cavity_drop_derivatives_kernel
