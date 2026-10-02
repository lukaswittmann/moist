!> Per-grid point sensitivity kernel shared by the DROP adjoint paths
!>
!> - Per-point map shared by the electronic and nuclear reverse-mode paths
!> - Input: perturbation of the projected point `r`, the multiplier `lambda` and
!>   the level-set jet `(S, grad S, grad^2 S)`
!> - Chain: tangent frame, closest-point Jacobian `J`, switching functions,
!>   Lebedev weight, Gaussian width, surface observables
!> - Linear in the perturbation: seed once per basis direction and sum
!> - [[build_seed_state]]: seed-independent state, once per grid point
!> - [[apply_seed]]: linear response to one seed, once per basis direction
!> - Degeneracy returned as a status code; fatal or skip is the caller's decision
module moist_cavity_drop_derivatives_kernel
   use mctc_env_accuracy, only: wp
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric, &
      & eig_2x2_offdiag_tol
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_type

   implicit none(type, external)
   private

   public :: drop_seed_state_type, drop_seed_result_type, drop_surface_weights_type
   public :: drop_seed_state_tangent_type, drop_seed_input_tangent_type

   public :: build_seed_state, apply_seed, apply_seed_tangent, compute_branch_phi_adj
   public :: branch_point_adjoint, seed_contribution
   public :: next_branch_group, max_branch_group_size
   public :: switched_eigenvalue_response, switched_eigenvalue_curvature
   public :: seed_status_message
   public :: drop_n_host_jet, add_host_jet
   public :: seed_state_ok, seed_state_singular_gradient
   public :: seed_state_singular_bmat, seed_state_singular_jacobian
   public :: seed_weight_tol, seed_det_b_guard, seed_curv_disc_guard

   !> Grid point is usable
   integer, parameter :: seed_state_ok = 0
   !> `|grad S|` vanished, surface normal undefined
   integer, parameter :: seed_state_singular_gradient = 1
   !> Tangent-restricted KKT matrix `B` is singular
   integer, parameter :: seed_state_singular_bmat = 2
   !> Closest-point Jacobian `J` vanished
   integer, parameter :: seed_state_singular_jacobian = 3

   !> Length of one host jet tangent
   !>
   !> - Level set and spatial derivatives to third order, `1 + 3 + 9 + 27`
   !> - Full Cartesian tensors, packed by order in Fortran order
   !> - Packing of the second-order host exchange, see `coupling_tangent_type`
   integer, parameter :: drop_n_host_jet = 40

   !> Magnitude below which a weight or a norm counts as zero
   real(wp), parameter :: seed_weight_tol = 1.0e-30_wp
   !> Guard on `det(B)` before the tangent-restricted inverse
   real(wp), parameter :: seed_det_b_guard = 1.0e-30_wp
   !> Discriminant below which the `k1`/`k2` split counts as umbilic
   !>
   !> - `dk1`, `dk2` ill-defined at `k1 = k2`
   !> - Mean and Gaussian curvature stay smooth
   real(wp), parameter :: seed_curv_disc_guard = 1.0e-10_wp
   !> Eigenvalue gap of `B` below which the switched eigenvector counts as degenerate
   !>
   !> - `d(u_switch)` scales as `1/gap`
   !> - Eigenvector arbitrary at `gap = 0`, no derivative exists
   real(wp), parameter :: seed_eig_gap_guard = 1.0e-10_wp

   !> Per-grid point forward state consumed by [[apply_seed]]
   !>
   !> - `Inputs` block: written by the caller before [[build_seed_state]]
   !> - `Derived` block: written by [[build_seed_state]]
   type :: drop_seed_state_type

      !* ----------------------------------- Inputs ----------------------------------- *!

      !> Level-set gradient and Hessian at the projected point
      real(wp) :: lsf1_r(3) = 0.0_wp, lsf2_rr(3, 3) = 0.0_wp
      !> Third spatial derivative of the level set at the projected point
      real(wp) :: lsf3_rrr(3, 3, 3) = 0.0_wp
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val = 0.0_wp
      !> Objective coefficient `phi_alpha`
      real(wp) :: alpha_coeff = 0.0_wp
      !> Anchor position and its owner sphere center
      real(wp) :: anchor(3) = 0.0_wp, owner_xyz(3) = 0.0_wp
      !> Grid-level weights entering the `wleb` / `xi` chain
      real(wp) :: anchor_wleb0 = 0.0_wp, cpjac_scal0 = 0.0_wp, w_f0 = 0.0_wp
      !> Softmax branch weight, Lebedev weight and Gaussian width
      real(wp) :: wbranch = 0.0_wp, wleb = 0.0_wp, xi0 = 0.0_wp
      !> Whether the caller needs the principal-curvature response
      logical :: want_curvature = .false.

      !* ----------------------------------- Derived ---------------------------------- *!

      !> Level-set gradient norm and its square
      real(wp) :: g_norm = 0.0_wp, g_norm_sq = 0.0_wp
      !> Tangent-restricted KKT matrix `A = alpha*I - lambda*H`
      real(wp) :: A_mat(3, 3) = 0.0_wp
      !> Outward normal and the surface tangent frame built from it
      real(wp) :: n_surf(3) = 0.0_wp, q1(3) = 0.0_wp, q2(3) = 0.0_wp
      !> `A` applied to the tangent frame
      real(wp) :: Aq1(3) = 0.0_wp, Aq2(3) = 0.0_wp
      !> `B = Q^T A Q` and its determinant
      real(wp) :: B11 = 0.0_wp, B12 = 0.0_wp, B22 = 0.0_wp, det_B = 0.0_wp
      !> Inverse of `B`
      real(wp) :: Binv11 = 0.0_wp, Binv12 = 0.0_wp, Binv22 = 0.0_wp
      !> Eigenvector of `B` for the switched eigenvalue, lifted to 3D
      real(wp) :: u_switch(3) = 0.0_wp
      !> Whether `vmin_B` is `[B12, lambda - B11]`, not a canonical basis vector
      !>
      !> - `vmin_norm` meaningful only when true
      logical :: vmin_offdiag = .false.
      !> Switched (smaller) eigenvalue of `B` and the gap `sqrt(disc)` to the larger one
      real(wp) :: lambda_switch = 0.0_wp, sqrt_disc_B = 0.0_wp
      !> Normalised 2D eigenvector of `B` and the norm it was divided by
      real(wp) :: vmin_B(2) = 0.0_wp, vmin_norm = 0.0_wp
      !> Sphere tangent frame
      real(wp) :: t1_vec(3) = 0.0_wp, t2_vec(3) = 0.0_wp
      !> Sphere frame projected onto the surface tangent frame
      real(wp) :: tau1(2) = 0.0_wp, tau2(2) = 0.0_wp
      !> `B^-1` images of `tau1`, `tau2`
      real(wp) :: w1(2) = 0.0_wp, w2(2) = 0.0_wp
      !> Lifted tangent vectors
      real(wp) :: y1(3) = 0.0_wp, y2(3) = 0.0_wp
      !> Cross product of `y1`, `y2` and `1/J`
      real(wp) :: cross_vec(3) = 0.0_wp, inv_J = 0.0_wp
      !> Gram-Schmidt axis of the tangent-frame derivative
      integer :: min_axis = 1
      !> Gram-Schmidt data behind `q1`
      real(wp) :: n_dot_q1 = 0.0_wp, proj_surf = 0.0_wp, v_norm_surf = 0.0_wp
      !> Critical-gradient switch value and slope
      real(wp) :: f_crit0 = 0.0_wp, f_crit_dS = 0.0_wp
      !> Focusing switch value and slope
      real(wp) :: f_foc_f0 = 0.0_wp, f_foc_dS = 0.0_wp
      !> Switching-function curvatures for the second-order chain
      real(wp) :: f_crit_d2S = 0.0_wp, f_foc_d2S = 0.0_wp
      !> Lebedev-weight pruning chain factor
      real(wp) :: wleb_prune_factor = 1.0_wp
      !> Whether Lebedev-weight pruning was active at build time
      logical :: use_wleb_prune = .false.
      !> Pre-pruning weight product, pruning-switch slope and curvature at `|w_pre_i|`
      !>
      !> - All zero when pruning is off
      real(wp) :: w_pre_i = 0.0_wp, f_wleb_ds = 0.0_wp, f_wleb_d2S = 0.0_wp
      !> Level-set Hessian on the surface tangent frame, only when `want_curvature`
      real(wp) :: Hq1(3) = 0.0_wp, Hq2(3) = 0.0_wp
      !> Shape operator in that frame, only when `want_curvature`
      real(wp) :: S11 = 0.0_wp, S12 = 0.0_wp, S22 = 0.0_wp
      !> Shape-operator trace and mean curvature, only when `want_curvature`
      real(wp) :: T_curv = 0.0_wp, KM_curv = 0.0_wp
      !> Half-difference and gap `|k1 - k2|/2`, a sum of squares; only when `want_curvature`
      real(wp) :: half_diff = 0.0_wp, disc_curv = 0.0_wp

   end type drop_seed_state_type

   !> Linear response of the per-point map to one seed
   !>
   !> - Also the second-order response of [[apply_seed_tangent]]: its `dg` is the
   !>   derivative of `res%dg` along the second direction
   !> - Default initialisation is load bearing: the `want_curvature` early return
   !>   of [[apply_seed_tangent]] leaves `dk1`, `dk2` untouched
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

   !> Tangent of the derived block of [[drop_seed_state_type]] along one seed
   !>
   !> - Returned by [[apply_seed]] through its optional `dstate`, no second routine:
   !>   one copy of the chain keeps the first-order path bit-for-bit unchanged
   !> - One seed direction per instance
   !> - Only derived fields [[apply_seed]] reads that are not frozen
   !> - `dn_surf`, `d_gnorm`, `dH` live on [[drop_seed_result_type]]
   !> - No `dt1_vec`, `dt2_vec`: the sphere tangent frame is rigid at every order,
   !>   `anchor - owner_xyz = R_own * u_leb`
   !> - They turn nonzero only for geometry-dependent radii, not implemented
   !> - `min_axis`, `want_curvature` are frozen discrete choices
   !> - Default initialisation is load bearing: the `want_curvature` early return
   !>   of [[apply_seed]] leaves the curvature block untouched
   type :: drop_seed_state_tangent_type

      !* ---------------------------- KKT matrix and frame ---------------------------- *!

      !> Sensitivity of the tangent-restricted KKT matrix `A`
      real(wp) :: dA_mat(3, 3) = 0.0_wp
      !> Sensitivity of the surface tangent frame
      real(wp) :: dq1(3) = 0.0_wp, dq2(3) = 0.0_wp
      !> Sensitivity of `A` applied to the tangent frame, both product-rule terms
      real(wp) :: dAq1(3) = 0.0_wp, dAq2(3) = 0.0_wp
      !> Sensitivity of `B = Q^T A Q` and of its determinant
      real(wp) :: dB11 = 0.0_wp, dB12 = 0.0_wp, dB22 = 0.0_wp, ddet_B = 0.0_wp
      !> Sensitivity of `B^-1`
      real(wp) :: dBinv11 = 0.0_wp, dBinv12 = 0.0_wp, dBinv22 = 0.0_wp
      !> Sensitivity of the switched eigenvalue, from the basis-invariant route
      real(wp) :: dlambda_switch = 0.0_wp
      !> `(dM) u` of `M = P A P` at the base eigenvector, not `d(M u)`
      !>
      !> - Depends on the seed alone; stored to keep the `b` chain out of the
      !>   direction-pair loop
      !> - Paired with `du_switch` in [[apply_seed_tangent]]; `dM` itself is never read
      real(wp) :: dM_u(3) = 0.0_wp
      !> Sensitivity of the switched eigenvector, in the `B` basis and lifted
      real(wp) :: dvmin_B(2) = 0.0_wp, du_switch(3) = 0.0_wp

      !* ---------------------- Lifted tangents and the Jacobian ---------------------- *!

      !> Sensitivity of the sphere frame projected onto the surface frame
      real(wp) :: dtau1(2) = 0.0_wp, dtau2(2) = 0.0_wp
      !> Sensitivity of the `B^-1` images of those projections
      real(wp) :: dw1(2) = 0.0_wp, dw2(2) = 0.0_wp
      !> Sensitivity of the lifted tangent vectors and of their cross product
      real(wp) :: dy1(3) = 0.0_wp, dy2(3) = 0.0_wp, dcross_vec(3) = 0.0_wp
      !> Sensitivity of `1/J`
      real(wp) :: dinv_J = 0.0_wp

      !* -------------------------------- Gram-Schmidt -------------------------------- *!

      !> Sensitivity of the Gram-Schmidt data behind `q1`
      real(wp) :: dn_dot_q1 = 0.0_wp, dproj_surf = 0.0_wp, dv_norm_surf = 0.0_wp

      !* ---------------------------- Switching and weights --------------------------- *!

      !> Sensitivity of the critical-gradient switch value and slope
      real(wp) :: df_crit0 = 0.0_wp, df_crit_dS = 0.0_wp
      !> Sensitivity of the focusing switch value and slope
      real(wp) :: df_foc_f0 = 0.0_wp, df_foc_dS = 0.0_wp
      !> Sensitivity of the Lebedev-weight pruning factor, zero when pruning is off
      real(wp) :: dwleb_prune_factor = 0.0_wp
      !> Sensitivity of `|grad S|^2`
      real(wp) :: dg_norm_sq = 0.0_wp

      !* ---------------------------- Curvature invariants ---------------------------- *!

      !> Total `d(H q_a) = dH q_a + H dq_a`
      !>
      !> - Depends on the seed alone; stored for the direction-pair loop of
      !>   [[apply_seed_tangent]]
      real(wp) :: dHq1(3) = 0.0_wp, dHq2(3) = 0.0_wp
      !> Sensitivity of the shape-operator entries
      real(wp) :: dS11 = 0.0_wp, dS12 = 0.0_wp, dS22 = 0.0_wp
      !> Sensitivity of the shape-operator trace and mean curvature
      real(wp) :: dT_curv = 0.0_wp, dKM_curv = 0.0_wp
      !> Sensitivity of the half-difference and the eigenvalue gap
      real(wp) :: dhalf_diff = 0.0_wp, ddisc_curv = 0.0_wp

   end type drop_seed_state_tangent_type

   !> v-direction tangent of the `Inputs` block of [[drop_seed_state_type]]
   !>
   !> - Passed explicitly, not rebuilt in [[apply_seed_tangent]]: identities such as
   !>   `dlsf1_r = res_v%dg`, `dcpjac_scal0 = res_v%dJ` belong to the driver
   !> - Keeps the kernel testable along an arbitrary path
   !> - Only fields [[apply_seed]] reads
   !> - No tangent of `alpha_coeff` on purpose: `param%phi_alpha` is fixed and
   !>   [[apply_seed]] has no seed channel for it
   !> - Geometry-dependent `alpha` needs this field and a `dalpha I` term in `ddA`
   !>   of [[apply_seed_tangent]], both together
   type :: drop_seed_input_tangent_type
      !> Tangent of the level-set gradient and Hessian at the projected point
      real(wp) :: dlsf1_r(3) = 0.0_wp, dlsf2_rr(3, 3) = 0.0_wp
      !> Tangent of the level-set third derivative at the projected point
      real(wp) :: dlsf3_rrr(3, 3, 3) = 0.0_wp
      !> Tangent of the multiplier
      real(wp) :: dlambda_val = 0.0_wp
      !> Tangent of the grid-level weight scalars
      real(wp) :: danchor_wleb0 = 0.0_wp, dcpjac_scal0 = 0.0_wp, dw_f0 = 0.0_wp
      !> Tangent of the branch weight, Lebedev weight and Gaussian width
      real(wp) :: dwbranch = 0.0_wp, dwleb = 0.0_wp, dxi0 = 0.0_wp
   end type drop_seed_input_tangent_type

   !> Surface adjoints reduced to the channels the seed loop reads
   !>
   !> - Area and integration-weight channels are derived: `a_i = R_I^2 f_i wleb_i`,
   !>   `w_i = wleb_i`, `xi_i = swx/(R_I sqrt(wleb_i))`
   !> - Hence `a = c f/xi^2`, `w = c/xi^2`: both fold into the width channel
   !> - Area also folds into the switching channel through `da/df = R_I^2 wleb_i`,
   !>   only for a parameter that moves `f`
   type :: drop_surface_weights_type
      !> Gaussian-width adjoint, with the area and weight channels folded in
      real(wp), allocatable :: w_xi(:)
      !> Switching adjoint, with the area fold only when requested
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
   !> - `Inputs` block of `state` filled by the caller first
   !> - Degenerate point: `status` set, derived fields left incomplete
   !>
   !> @param[in,out] state           inputs read, derived fields written
   !> @param[in]     f_crit          critical-gradient switching function
   !> @param[in]     f_foc           focusing switching function
   !> @param[in]     f_wleb          Lebedev-weight pruning switching function
   !> @param[in]     use_wleb_prune  whether Lebedev-weight pruning is active
   !> @param[out]    status          one of the `seed_state_*` codes
   subroutine build_seed_state(state, f_crit, f_foc, f_wleb, use_wleb_prune, status)
      !> Seed state
      type(drop_seed_state_type), intent(inout) :: state
      !> Switching functions owned by the cavity
      class(moist_cavity_drop_swif_type), intent(in) :: f_crit, f_foc, f_wleb
      !> Whether Lebedev-weight pruning is active
      logical, intent(in) :: use_wleb_prune
      !> Degeneracy status
      integer, intent(out) :: status

      !> Scratch
      real(wp) :: lambda_switch, beta_max
      real(wp) :: vmin_B(2), vmax_B(2)
      !> Lebedev pruning scratch
      real(wp) :: w_pre_i, f_wleb_s, f_wleb_ds, f_wleb_d2s
      !> Off-diagonal entry of the unnormalised `vmin_B`
      real(wp) :: vmin_off

      status = seed_state_ok

      state%g_norm_sq = dot_product(state%lsf1_r, state%lsf1_r)
      state%g_norm = sqrt(state%g_norm_sq)
      if (state%g_norm <= seed_weight_tol) then
         status = seed_state_singular_gradient
         return
      end if
      call f_crit%eval(state%g_norm, state%f_crit0, state%f_crit_dS, state%f_crit_d2S)

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
      state%det_B = state%B11*state%B22 - state%B12*state%B12
      if (abs(state%det_B) <= seed_det_b_guard) then
         status = seed_state_singular_bmat
         return
      end if
      call eig_2x2_symmetric(state%B11, state%B12, state%B22, lambda_switch, beta_max, &
                             vmin_B, vmax_B)
      state%u_switch = vmin_B(1)*state%q1 + vmin_B(2)*state%q2
      ! Eigen data for the tangent; `vmin_norm` recomputed, not returned by the solver
      state%lambda_switch = lambda_switch
      ! Gap as `hypot(B11 - B22, 2 B12)`, not `beta_max - lambda_switch`
      ! - `trace^2 - 4 det` cancels to zero for gaps below `sqrt(eps)*|trace|`
      ! - The branch test below does not protect the division by the gap
      state%sqrt_disc_B = hypot(state%B11 - state%B22, 2.0_wp*state%B12)
      state%vmin_B = vmin_B
      state%vmin_offdiag = abs(state%B12) > eig_2x2_offdiag_tol
      if (state%vmin_offdiag) then
         vmin_off = lambda_switch - state%B11
         state%vmin_norm = sqrt(state%B12*state%B12 + vmin_off*vmin_off)
      else
         state%vmin_norm = 0.0_wp
      end if
      call f_foc%eval(lambda_switch, state%f_foc_f0, state%f_foc_dS, state%f_foc_d2S)

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

      ! d(w_pre * S)/dp = (S + |w_pre|*S') * d(w_pre)/dp
      state%use_wleb_prune = use_wleb_prune
      if (use_wleb_prune) then
         w_pre_i = state%anchor_wleb0*state%cpjac_scal0*state%w_f0
         call f_wleb%eval(abs(w_pre_i), f_wleb_s, f_wleb_ds, f_wleb_d2s)
         state%wleb_prune_factor = f_wleb_s + abs(w_pre_i)*f_wleb_ds
         state%w_pre_i = w_pre_i
         state%f_wleb_ds = f_wleb_ds
         state%f_wleb_d2S = f_wleb_d2s
      else
         state%wleb_prune_factor = 1.0_wp
         state%w_pre_i = 0.0_wp
         state%f_wleb_ds = 0.0_wp
         state%f_wleb_d2S = 0.0_wp
      end if

      ! Shape operator in the surface tangent frame Q = [q1, q2]
      !   S_ab = q_a^T H q_b / |g|,   a, b in {1, 2}
      !   k1,k2 = KM +/- sqrt(half_diff^2 + S12^2),   k1 >= k2
      !   KM = (S11 + S22)/2,   half_diff = (S11 - S22)/2
      ! - Discriminant as a sum of squares, not `sqrt(KM^2 - KG)`: keeps small gaps,
      !   [[apply_seed]] divides by it; same form as `properties.f90`
      ! - `KM`, `disc` invariant under a rotation of Q: any smooth orthonormal frame
      !   works when differentiated consistently, as `dq1`/`dq2` are
      if (state%want_curvature) then
         associate (H => state%lsf2_rr)
            state%Hq1 = matmul(H, state%q1)
            state%Hq2 = matmul(H, state%q2)
            state%S11 = dot_product(state%q1, state%Hq1)/state%g_norm
            state%S12 = dot_product(state%q1, state%Hq2)/state%g_norm
            state%S22 = dot_product(state%q2, state%Hq2)/state%g_norm

            state%T_curv = state%S11 + state%S22
            state%KM_curv = 0.5_wp*state%T_curv
            state%half_diff = 0.5_wp*(state%S11 - state%S22)
            state%disc_curv = sqrt(state%half_diff*state%half_diff &
                                   + state%S12*state%S12)
         end associate
      end if

   end subroutine build_seed_state

   !> Response of the switched eigenvalue to one seed, contracted
   !>
   !> - `dlambda = u . (dM u)`, `u` the base eigenvector
   !> - `M = P A P`, `P = I - n n^T`
   !> - Contracted form avoids the `dM` assembly and its `spread` heap temporaries
   !> - Precondition 1: `n . u = 0`, by construction in [[setup_tangent_frame]]
   !> - Precondition 2: `P u = u`, follows from 1
   !> - Precondition 3: `A` symmetric; the factor of two rests on
   !>   `u . (A n) = n . (A u)`
   !> - `dA` need not be symmetric: both outputs are exact for arbitrary `dA`
   !> - `dM u` feeds the caller's collapse `2 du_v . (dM_b u)`, which needs `dM_b`
   !>   and hence `dA_b` symmetric; see [[apply_seed_tangent]]
   !>
   !> With `a = dn . u`, so `dP u = -a n`
   !>
   !>     dlambda = u . (dA u) - 2 a (n . A u)
   !>     dM u    = -(dn (n . Au) + n (dn . Au))
   !>               + (dA u - n (n . dA u))
   !>               - a (An - n (n . An))
   !>
   !> @param[in]  n_surf          unit length
   !> @param[in]  u_switch        tangent to the surface
   !> @param[in]  A_mat           symmetric
   !> @param[in]  dA              need not be symmetric
   !> @param[in]  dn              sensitivity of the normal
   !> @param[out] dlambda_switch  `u . (dM u)`
   !> @param[out] dM_u            optional, for the second-order chain
   pure subroutine switched_eigenvalue_response(n_surf, u_switch, A_mat, dA, dn, &
                                                dlambda_switch, dM_u)
      !> Outward unit normal and the switched eigenvector
      real(wp), intent(in) :: n_surf(3), u_switch(3)
      !> Symmetric KKT matrix and its sensitivity, which need not be symmetric
      real(wp), intent(in) :: A_mat(3, 3), dA(3, 3)
      !> Sensitivity of the normal
      real(wp), intent(in) :: dn(3)
      !> Response of the switched eigenvalue
      real(wp), intent(out) :: dlambda_switch
      !> `dM u`, formed only when present
      real(wp), intent(out), optional :: dM_u(3)

      !> Matrix-vector products the identity is built from
      real(wp) :: Au(3), An(3), dA_u(3)
      !> `dn . u` and `n . Au`
      real(wp) :: a_dn, n_dot_Au

      Au = matmul(A_mat, u_switch)
      dA_u = matmul(dA, u_switch)
      a_dn = dot_product(dn, u_switch)
      n_dot_Au = dot_product(n_surf, Au)

      dlambda_switch = dot_product(u_switch, dA_u) - 2.0_wp*a_dn*n_dot_Au

      if (present(dM_u)) then
         An = matmul(A_mat, n_surf)
         dM_u = -(dn*n_dot_Au + n_surf*dot_product(dn, Au)) &
                + (dA_u - n_surf*dot_product(n_surf, dA_u)) &
                - a_dn*(An - n_surf*dot_product(n_surf, An))
      end if

   end subroutine switched_eigenvalue_response

   !> Curvature of the switched eigenvalue in two directions, contracted
   !>
   !> - Returns `u . (ddM u)` of `M = P A P`
   !> - `ddX` means `d_v(d_b X)`
   !> - Innermost `(b, v)` loop: contracted form avoids the `ddM` assembly and its
   !>   `spread` heap temporaries
   !> - Preconditions 1 to 3 of [[switched_eigenvalue_response]] carry over
   !> - Additionally `dA_b`, `dA_v` symmetric: the second bracket is four terms
   !>   collapsed through `u . (dA n) = n . (dA u)`
   !> - Symmetric `dA` is a condition on the seed's `dlsf2_rr`, broken by the
   !>   single-entry Hessian seeds of [[seed_jet_basis]]; see [[apply_seed_tangent]]
   !> - Only the symmetric part of `ddA` is read
   !> - Eigenvector cross terms `2 du_v . (dM_b u)` are added by the caller
   !>
   !> With `a_b = dn_b . u`, `a_v = dn_v . u` and `c = ddn . u`
   !>
   !>     u . (ddM u) = -2 [ a_v (dn_b . Au) + a_b (dn_v . Au) + c (n . Au) ]
   !>                   -2 [ a_b (n . dA_v u) + a_v (n . dA_b u) ]
   !>                   +2 a_b a_v (n . An)
   !>                   + u . (ddA u)
   !>
   !> @param[in]  n_surf           unit length
   !> @param[in]  u_switch         tangent to the surface
   !> @param[in]  A_mat            symmetric
   !> @param[in]  dA_b             along `b`, required symmetric
   !> @param[in]  dA_v             along `v`, required symmetric
   !> @param[in]  ddA              only its symmetric part is read
   !> @param[in]  dn_b             sensitivity of the normal along `b`
   !> @param[in]  dn_v             sensitivity of the normal along `v`
   !> @param[in]  ddn              second-order sensitivity of the normal
   !> @param[out] ddlambda_switch  `u . (ddM u)`
   pure subroutine switched_eigenvalue_curvature(n_surf, u_switch, A_mat, &
                                                 dA_b, dA_v, ddA, &
                                                 dn_b, dn_v, ddn, ddlambda_switch)
      !> Outward unit normal and the switched eigenvector
      real(wp), intent(in) :: n_surf(3), u_switch(3)
      !> Tangent-restricted KKT matrix, symmetric
      real(wp), intent(in) :: A_mat(3, 3)
      !> Symmetric first-order sensitivities of `A` and the second-order one
      real(wp), intent(in) :: dA_b(3, 3), dA_v(3, 3), ddA(3, 3)
      !> First- and second-order sensitivities of the normal
      real(wp), intent(in) :: dn_b(3), dn_v(3), ddn(3)
      !> `u . (ddM u)`
      real(wp), intent(out) :: ddlambda_switch

      !> Matrix-vector products the identity is built from
      real(wp) :: Au(3), An(3), dA_b_u(3), dA_v_u(3)
      !> `dn_b . u`, `dn_v . u`, `ddn . u` and `n . Au`
      real(wp) :: a_b, a_v, c_ddn, n_dot_Au

      Au = matmul(A_mat, u_switch)
      An = matmul(A_mat, n_surf)
      dA_b_u = matmul(dA_b, u_switch)
      dA_v_u = matmul(dA_v, u_switch)

      a_b = dot_product(dn_b, u_switch)
      a_v = dot_product(dn_v, u_switch)
      c_ddn = dot_product(ddn, u_switch)
      n_dot_Au = dot_product(n_surf, Au)

      ddlambda_switch = -2.0_wp*(a_v*dot_product(dn_b, Au) &
                                 + a_b*dot_product(dn_v, Au) &
                                 + c_ddn*n_dot_Au) &
                        - 2.0_wp*(a_b*dot_product(n_surf, dA_v_u) &
                                  + a_v*dot_product(n_surf, dA_b_u)) &
                        + 2.0_wp*a_b*a_v*dot_product(n_surf, An) &
                        + dot_product(u_switch, matmul(ddA, u_switch))

   end subroutine switched_eigenvalue_curvature

   !> Propagate one seed through the per-grid point map
   !>
   !> - Seed: jet perturbation at the fixed point (`dlsf1_r`, `dlsf2_rr`) plus the
   !>   induced `dr`, `dlambda` from the bordered KKT system
   !> - Level-set value perturbation enters only through that system, no argument
   !> - `dstate` fields computed only when it is present
   !>
   !> @param[in]  state     forward state from [[build_seed_state]]
   !> @param[in]  dlsf1_r   seed perturbation of `grad S` at fixed `r`
   !> @param[in]  dlsf2_rr  seed perturbation of `grad^2 S` at fixed `r`
   !> @param[in]  dr        induced motion of the projected point
   !> @param[in]  dlambda   induced change of the Lagrange multiplier
   !> @param[out] res       linear response of the per-point map
   !> @param[out] dstate    optional, tangent of the derived state along this seed
   pure subroutine apply_seed(state, dlsf1_r, dlsf2_rr, dr, dlambda, res, dstate)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Seed perturbation of the level-set gradient and Hessian
      real(wp), intent(in) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Induced motion of the projected point and its multiplier
      real(wp), intent(in) :: dr(3), dlambda
      !> Linear response
      type(drop_seed_result_type), intent(out) :: res
      !> Tangent of the derived block of `state`, computed only when present
      type(drop_seed_state_tangent_type), intent(out), optional :: dstate

      !> Sensitivity of the tangent-restricted KKT matrix
      real(wp) :: dA(3, 3)
      !> Tangent-frame derivative scratch
      real(wp) :: v_tmp(3), dq1(3), dq2(3)
      !> `B` and `B^-1` sensitivities
      real(wp) :: dAq1(3), dAq2(3), dB11, dB12, dB22, ddet_B
      real(wp) :: dBinv11, dBinv12, dBinv22
      !> Switched-eigenvalue response
      real(wp) :: dlambda_switch
      !> Lifted tangent-vector sensitivities
      real(wp) :: dtau1(2), dtau2(2), dw1(2), dw2(2)
      real(wp) :: dy1(3), dy2(3), dcross(3)
      !> Lebedev-weight chain scratch
      real(wp) :: dw_pre
      !> Curvature sensitivities; `dHq1_p`/`dHq2_p` are partial, `dH q_a` without `H dq_a`
      real(wp) :: dHq1_p(3), dHq2_p(3)
      real(wp) :: dN11, dN12, dN22, dS11, dS12, dS22
      real(wp) :: dT, dhalf_diff, d_disc
      !> Derivative of the `eig_2x2_symmetric` construction, only for `dstate`
      real(wp) :: dtrace_B, ddisc_B, dsqrt_disc_B, dlambda_switch_eig
      real(wp) :: dvmin_raw(2)
      !> Cartesian index
      integer :: kaxis

      res%dg = dlsf1_r + matmul(state%lsf2_rr, dr)
      res%dH = dlsf2_rr
      do kaxis = 1, 3
         res%dH(:, :) = res%dH(:, :) + state%lsf3_rrr(:, :, kaxis)*dr(kaxis)
      end do
      dA = -dlambda*state%lsf2_rr - state%lambda_val*res%dH

      res%dn_surf = (res%dg - state%n_surf*dot_product(state%n_surf, res%dg))/state%g_norm

      ! dQ/dp: q1 from Gram-Schmidt of e_k against n, q2 = n x q1
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

      ! dB/dp with B = Q^T A Q, A-symmetry used for the mixed term
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

      ! `u . (dM u)` contracted, see [[switched_eigenvalue_response]]
      ! - `dM_u` formed only for the state tangent
      if (present(dstate)) then
         call switched_eigenvalue_response(state%n_surf, state%u_switch, state%A_mat, &
                                           dA, res%dn_surf, dlambda_switch, dstate%dM_u)
      else
         call switched_eigenvalue_response(state%n_surf, state%u_switch, state%A_mat, &
                                           dA, res%dn_surf, dlambda_switch)
      end if

      if (present(dstate)) then
         dstate%dA_mat = dA
         dstate%dq1 = dq1
         dstate%dq2 = dq2
         ! Scratch `dAq1`/`dAq2` hold `dA q` only; the state tangent adds `A dq`
         dstate%dAq1 = dAq1 + matmul(state%A_mat, dq1)
         dstate%dAq2 = dAq2 + matmul(state%A_mat, dq2)
         dstate%dB11 = dB11
         dstate%dB12 = dB12
         dstate%dB22 = dB22
         dstate%ddet_B = ddet_B
         dstate%dBinv11 = dBinv11
         dstate%dBinv12 = dBinv12
         dstate%dBinv22 = dBinv22
         dstate%dlambda_switch = dlambda_switch

         ! d(v_min): `eig_2x2_symmetric` differentiated line by line, keeps its sign
         ! convention `v_min = [B12, lambda_min - B11]/norm`
         ! - Off the off-diagonal branch the primal is a canonical basis vector:
         !   derivative zero, `vmin_norm` not a norm the primal used
         ! - Vanishing gap: eigenvector arbitrary, zero is the guard's answer and not
         !   the derivative, as for `seed_curv_disc_guard`
         if (state%vmin_offdiag .and. state%sqrt_disc_B > seed_eig_gap_guard) then
            dtrace_B = dB11 + dB22
            ! d(disc) = d((B11 - B22)^2 + 4 B12^2), accurate at small gaps
            ddisc_B = 2.0_wp*(state%B11 - state%B22)*(dB11 - dB22) &
                      + 8.0_wp*state%B12*dB12
            dsqrt_disc_B = ddisc_B/(2.0_wp*state%sqrt_disc_B)
            dlambda_switch_eig = 0.5_wp*(dtrace_B - dsqrt_disc_B)
            dvmin_raw(1) = dB12
            dvmin_raw(2) = dlambda_switch_eig - dB11
            dstate%dvmin_B = (dvmin_raw &
                              - state%vmin_B*dot_product(state%vmin_B, dvmin_raw)) &
                             /state%vmin_norm
         else
            dstate%dvmin_B = 0.0_wp
         end if
         dstate%du_switch = dstate%dvmin_B(1)*state%q1 + state%vmin_B(1)*dq1 &
                            + dstate%dvmin_B(2)*state%q2 + state%vmin_B(2)*dq2

         dstate%dn_dot_q1 = res%dn_surf(state%min_axis)
         dstate%dproj_surf = -2.0_wp*state%n_dot_q1*dstate%dn_dot_q1
         ! Mirror of the primal `max(proj_surf, 1e-30)` clamp: norm constant once it bites
         if (state%proj_surf > 1.0e-30_wp) then
            dstate%dv_norm_surf = dstate%dproj_surf/(2.0_wp*state%v_norm_surf)
         else
            dstate%dv_norm_surf = 0.0_wp
         end if
      end if

      ! Sphere tangent frame is rigid, dt1 = dt2 = 0
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

      if (present(dstate)) then
         dstate%dtau1 = dtau1
         dstate%dtau2 = dtau2
         dstate%dw1 = dw1
         dstate%dw2 = dw2
         dstate%dy1 = dy1
         dstate%dy2 = dy2
         dstate%dcross_vec = dcross
         dstate%dinv_J = -res%dJ*state%inv_J*state%inv_J
         dstate%dg_norm_sq = 2.0_wp*state%g_norm*res%d_gnorm
         dstate%df_crit0 = state%f_crit_dS*res%d_gnorm
         dstate%df_crit_dS = state%f_crit_d2S*res%d_gnorm
         dstate%df_foc_f0 = state%f_foc_dS*dlambda_switch
         dstate%df_foc_dS = state%f_foc_d2S*dlambda_switch
         ! wleb_prune_factor = S(|w|) + |w| S'(|w|) with d|w| = sign(w_pre_i) * dw_pre
         if (state%use_wleb_prune) then
            dstate%dwleb_prune_factor = sign(1.0_wp, state%w_pre_i)*dw_pre &
                                        *(2.0_wp*state%f_wleb_ds &
                                          + abs(state%w_pre_i)*state%f_wleb_d2S)
         else
            dstate%dwleb_prune_factor = 0.0_wp
         end if
      end if

      if (.not. state%want_curvature) return

      associate (H => state%lsf2_rr, dH => res%dH)
         ! d(q_a^T H q_b) = dq_a . H q_b + dq_b . H q_a + q_a^T dH q_b
         ! - H-symmetry used for the mixed term, as in `dB12`
         dHq1_p = matmul(dH, state%q1)
         dHq2_p = matmul(dH, state%q2)
         dN11 = 2.0_wp*dot_product(dq1, state%Hq1) + dot_product(state%q1, dHq1_p)
         dN12 = dot_product(dq1, state%Hq2) + dot_product(dq2, state%Hq1) &
                + dot_product(state%q1, dHq2_p)
         dN22 = 2.0_wp*dot_product(dq2, state%Hq2) + dot_product(state%q2, dHq2_p)

         ! S_ab = N_ab/|g|
         dS11 = dN11/state%g_norm - state%S11*res%d_gnorm/state%g_norm
         dS12 = dN12/state%g_norm - state%S12*res%d_gnorm/state%g_norm
         dS22 = dN22/state%g_norm - state%S22*res%d_gnorm/state%g_norm

         dT = dS11 + dS22
         dhalf_diff = 0.5_wp*(dS11 - dS22)

         ! disc d(disc) = half_diff d(half_diff) + S12 dS12
         ! - Both products are O(disc); the quotient keeps the relative accuracy of `disc`
         if (state%disc_curv > seed_curv_disc_guard) then
            d_disc = (state%half_diff*dhalf_diff + state%S12*dS12)/state%disc_curv
         else
            d_disc = 0.0_wp
         end if
         res%dk1 = 0.5_wp*dT + d_disc
         res%dk2 = 0.5_wp*dT - d_disc

         if (present(dstate)) then
            ! Total d(H q_a), the form [[apply_seed_tangent]] contracts
            dstate%dHq1 = dHq1_p + matmul(H, dq1)
            dstate%dHq2 = dHq2_p + matmul(H, dq2)
            dstate%dS11 = dS11
            dstate%dS12 = dS12
            dstate%dS22 = dS22
            dstate%dT_curv = dT
            dstate%dKM_curv = 0.5_wp*dT
            dstate%dhalf_diff = dhalf_diff
            dstate%ddisc_curv = d_disc
         end if
      end associate

   end subroutine apply_seed

   !> Directional derivative of [[apply_seed]] along a second direction `v`
   !>
   !> - Returns `d/dv [ res_b ]` for the seed `b = (dlsf1_r, dlsf2_rr, dr, dlambda)`
   !> - Product rule on the code of [[apply_seed]], line by line; not a second
   !>   derivation of the geometry
   !> - PRECONDITION: the seed's `dlsf2_rr` is symmetric, or the caller sums the
   !>   transpose pair `E_ij + E_ji` before reading `dres`
   !> - Reason: three factor-of-two collapses need symmetric `dA_b`, namely
   !>   `2 du_v . (dM_b u)`, the second bracket of [[switched_eigenvalue_curvature]]
   !>   and `dBij = dqi . Aqj + qi . dAqj`
   !> - Broken by the off-diagonal Hessian seeds of [[seed_jet_basis]]: per-seed
   !>   `dres` is wrong at order 100 %, signs included
   !> - Equivalent to the pair sum: contract all nine seeds against a weight
   !>   symmetric in `(i, j)`; same rescue as the `dB12` shortcut of [[apply_seed]]
   !> - Not caught by the tests: every fixture seeds a symmetric Hessian
   !> - `res_b`, `dstate_b` come from one hoisted [[apply_seed]] call per seed and
   !>   are never recomputed here
   !> - `dstate_v`, `res_v` must be the true v-tangents of the derived state for the
   !>   input displacement `dinp_v`
   !> - [[apply_seed]] pairs are consistent only for `dinp_v%danchor_wleb0 = 0` and,
   !>   with pruning, `dinp_v%dcpjac_scal0 = res_v%dJ`, `dinp_v%dw_f0 = res_v%dw_f`;
   !>   the driver makes these identifications, other input tangents are free
   !> - `d_v(lsf3_rrr)` enters only through `dres%dH`; `lsf4_rrrr` is the driver's
   !> - No `dlsf1_r`, `dlsf2_rr` arguments: their coefficient is the identity, only
   !>   `ddlsf1_r`, `ddlsf2_rr` appear
   !> - TODO: retire the precondition by symmetrising the seed in
   !>   [[seed_jet_basis]] or on entry to [[apply_seed]]; too late here, since
   !>   `dstate_b%dM_u` is already built
   !> - Safe on the nuclear path (`f3_rr_rA` symmetric, only summation order
   !>   moves); undone since the LSF interface documents `w2` as a general `3x3`
   !>   and `w_lsf2` is caller-supplied through `api.f90`
   !>
   !> @param[in]  state      forward state from [[build_seed_state]]
   !> @param[in]  dstate_v   true v-tangent of the derived state
   !> @param[in]  dinp_v     input displacement along `v`
   !> @param[in]  res_v      response of direction `v`, consistent with `dstate_v`
   !> @param[in]  dr         induced point motion of seed `b`
   !> @param[in]  dlambda    induced multiplier change of seed `b`
   !> @param[in]  ddlsf1_r   v-tangent of the seed's `grad S` perturbation
   !> @param[in]  ddlsf2_rr  v-tangent of the seed's `grad^2 S` perturbation
   !> @param[in]  ddr        v-tangent of `dr`
   !> @param[in]  ddlambda   v-tangent of `dlambda`
   !> @param[in]  res_b      response of seed `b`, from [[apply_seed]]
   !> @param[in]  dstate_b   state tangent of seed `b`, from [[apply_seed]]
   !> @param[out] dres       second-order response
   pure subroutine apply_seed_tangent(state, dstate_v, dinp_v, res_v, dr, dlambda, &
                                      ddlsf1_r, ddlsf2_rr, ddr, ddlambda, &
                                      res_b, dstate_b, dres)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Tangent of the derived state along the second direction `v`
      type(drop_seed_state_tangent_type), intent(in) :: dstate_v
      !> Tangent of the state inputs along `v`
      type(drop_seed_input_tangent_type), intent(in) :: dinp_v
      !> Response of the `v` direction, carrying `dn_surf`, `d_gnorm` and `dH`
      type(drop_seed_result_type), intent(in) :: res_v
      !> Induced point motion and multiplier change of seed `b`
      real(wp), intent(in) :: dr(3), dlambda
      !> Tangent of that seed along `v`
      real(wp), intent(in) :: ddlsf1_r(3), ddlsf2_rr(3, 3), ddr(3), ddlambda
      !> Response of seed `b`, from [[apply_seed]]
      type(drop_seed_result_type), intent(in) :: res_b
      !> State tangent of seed `b`, from [[apply_seed]]
      type(drop_seed_state_tangent_type), intent(in) :: dstate_b
      !> Second-order response
      type(drop_seed_result_type), intent(out) :: dres

      !> Second-order sensitivity of the tangent-restricted KKT matrix
      real(wp) :: ddA(3, 3)
      !> Gram-Schmidt scratch of [[apply_seed]] for seed `b`, and its tangent
      real(wp) :: v_tmp_b(3), dv_tmp(3), ddq1(3), ddq2(3)
      !> Second-order `B` and `B^-1` sensitivities
      real(wp) :: ddAq1(3), ddAq2(3), ddB11, ddB12, ddB22, dddet_B
      real(wp) :: ddBinv11, ddBinv12, ddBinv22
      !> Second-order switched-eigenvalue response, and its `u . (ddM u)` part
      real(wp) :: ddlambda_switch, ddlambda_curv
      !> Second-order lifted tangent-vector sensitivities
      real(wp) :: ddtau1(2), ddtau2(2), ddw1(2), ddw2(2)
      real(wp) :: ddy1(3), ddy2(3), ddcross(3)
      !> Lebedev-weight and Jacobian chain scratch
      real(wp) :: dw_pre_b, ddw_pre, cross_dot_b
      !> Curvature scratch; `dH_b q_a` is the one partial `dstate%dHq_a` does not carry
      real(wp) :: dHq1_b(3), dHq2_b(3)
      real(wp) :: dN11_b, dN12_b, dN22_b, ddN11, ddN12, ddN22
      real(wp) :: ddS11, ddS12, ddS22, ddT, ddhalf_diff, dd_disc
      !> Cartesian index
      integer :: kaxis

      !* ================================ Level-set jet =============================== *!

      dres%dg = ddlsf1_r + matmul(dinp_v%dlsf2_rr, dr) + matmul(state%lsf2_rr, ddr)

      dres%dH = ddlsf2_rr
      do kaxis = 1, 3
         dres%dH(:, :) = dres%dH(:, :) + dinp_v%dlsf3_rrr(:, :, kaxis)*dr(kaxis) &
                         + state%lsf3_rrr(:, :, kaxis)*ddr(kaxis)
      end do

      ddA = -ddlambda*state%lsf2_rr - dlambda*dinp_v%dlsf2_rr &
            - dinp_v%dlambda_val*res_b%dH - state%lambda_val*dres%dH

      ! `d_gnorm` hoisted above its [[apply_seed]] position: `dn_surf` reuses `n . dg`
      dres%d_gnorm = dot_product(res_v%dn_surf, res_b%dg) &
                     + dot_product(state%n_surf, dres%dg)
      dres%dn_surf = (dres%dg - res_v%dn_surf*res_b%d_gnorm &
                      - state%n_surf*dres%d_gnorm)/state%g_norm &
                     - res_b%dn_surf*res_v%d_gnorm/state%g_norm

      !* ================================ Tangent frame =============================== *!

      ! `min_axis` is frozen, Gram-Schmidt axis fixed
      v_tmp_b = -res_b%dn_surf(state%min_axis)*state%n_surf &
                - state%n_dot_q1*res_b%dn_surf
      dv_tmp = -dres%dn_surf(state%min_axis)*state%n_surf &
               - res_b%dn_surf(state%min_axis)*res_v%dn_surf &
               - dstate_v%dn_dot_q1*res_b%dn_surf &
               - state%n_dot_q1*dres%dn_surf

      ! Guard rule: every branch tests the primal's condition on the primal's
      ! threshold, zero on the else
      ! - Raise a gate only in [[apply_seed]] and here together
      ! - Eigenvector and `use_wleb_prune` gates of [[apply_seed]] are absent: they
      !   guard only `dstate%dvmin_B` and `dstate%dwleb_prune_factor`, which arrive
      !   already guarded in `dstate_v`
      if (state%proj_surf > 1.0e-30_wp) then
         ddq1 = (dv_tmp - dstate_v%dq1*dot_product(state%q1, v_tmp_b) &
                 - state%q1*(dot_product(dstate_v%dq1, v_tmp_b) &
                             + dot_product(state%q1, dv_tmp)))/state%v_norm_surf &
                - dstate_b%dq1*dstate_v%dv_norm_surf/state%v_norm_surf
      else
         ddq1 = 0.0_wp
      end if

      ddq2(1) = dres%dn_surf(2)*state%q1(3) - dres%dn_surf(3)*state%q1(2) &
                + res_b%dn_surf(2)*dstate_v%dq1(3) - res_b%dn_surf(3)*dstate_v%dq1(2) &
                + res_v%dn_surf(2)*dstate_b%dq1(3) - res_v%dn_surf(3)*dstate_b%dq1(2) &
                + state%n_surf(2)*ddq1(3) - state%n_surf(3)*ddq1(2)
      ddq2(2) = dres%dn_surf(3)*state%q1(1) - dres%dn_surf(1)*state%q1(3) &
                + res_b%dn_surf(3)*dstate_v%dq1(1) - res_b%dn_surf(1)*dstate_v%dq1(3) &
                + res_v%dn_surf(3)*dstate_b%dq1(1) - res_v%dn_surf(1)*dstate_b%dq1(3) &
                + state%n_surf(3)*ddq1(1) - state%n_surf(1)*ddq1(3)
      ddq2(3) = dres%dn_surf(1)*state%q1(2) - dres%dn_surf(2)*state%q1(1) &
                + res_b%dn_surf(1)*dstate_v%dq1(2) - res_b%dn_surf(2)*dstate_v%dq1(1) &
                + res_v%dn_surf(1)*dstate_b%dq1(2) - res_v%dn_surf(2)*dstate_b%dq1(1) &
                + state%n_surf(1)*ddq1(2) - state%n_surf(2)*ddq1(1)

      !* ================================= KKT matrix ================================= *!

      ! Stored `dAq1`/`dAq2` are the full `dA q + A dq`, unlike the [[apply_seed]] scratch
      ! Form differentiated below, valid for symmetric `A` and the seed precondition
      !    dBij = dqi . Aqj + qi . dAqj
      ddAq1 = matmul(ddA, state%q1) + matmul(dstate_b%dA_mat, dstate_v%dq1) &
              + matmul(dstate_v%dA_mat, dstate_b%dq1) + matmul(state%A_mat, ddq1)
      ddAq2 = matmul(ddA, state%q2) + matmul(dstate_b%dA_mat, dstate_v%dq2) &
              + matmul(dstate_v%dA_mat, dstate_b%dq2) + matmul(state%A_mat, ddq2)

      ddB11 = dot_product(ddq1, state%Aq1) &
              + dot_product(dstate_b%dq1, dstate_v%dAq1) &
              + dot_product(dstate_v%dq1, dstate_b%dAq1) &
              + dot_product(state%q1, ddAq1)
      ddB12 = dot_product(ddq1, state%Aq2) &
              + dot_product(dstate_b%dq1, dstate_v%dAq2) &
              + dot_product(dstate_v%dq1, dstate_b%dAq2) &
              + dot_product(state%q1, ddAq2)
      ddB22 = dot_product(ddq2, state%Aq2) &
              + dot_product(dstate_b%dq2, dstate_v%dAq2) &
              + dot_product(dstate_v%dq2, dstate_b%dAq2) &
              + dot_product(state%q2, ddAq2)

      dddet_B = ddB11*state%B22 + dstate_b%dB11*dstate_v%dB22 &
                + dstate_v%dB11*dstate_b%dB22 + state%B11*ddB22 &
                - 2.0_wp*(dstate_v%dB12*dstate_b%dB12 + state%B12*ddB12)

      ! `dBinvXY = N/det^2`: outer quotient rule reuses the stored first-order value
      ddBinv11 = (ddB22*state%det_B + dstate_b%dB22*dstate_v%ddet_B &
                  - dstate_v%dB22*dstate_b%ddet_B - state%B22*dddet_B) &
                 /(state%det_B*state%det_B) &
                 - 2.0_wp*dstate_b%dBinv11*dstate_v%ddet_B/state%det_B
      ddBinv12 = (-ddB12*state%det_B - dstate_b%dB12*dstate_v%ddet_B &
                  + dstate_v%dB12*dstate_b%ddet_B + state%B12*dddet_B) &
                 /(state%det_B*state%det_B) &
                 - 2.0_wp*dstate_b%dBinv12*dstate_v%ddet_B/state%det_B
      ddBinv22 = (ddB11*state%det_B + dstate_b%dB11*dstate_v%ddet_B &
                  - dstate_v%dB11*dstate_b%ddet_B - state%B11*dddet_B) &
                 /(state%det_B*state%det_B) &
                 - 2.0_wp*dstate_b%dBinv22*dstate_v%ddet_B/state%det_B

      !* ========================= Basis-invariant eigenvalue ========================= *!

      ! `u . (ddM u)` contracted, see [[switched_eigenvalue_curvature]]
      ! - Both cross terms `dn_b dn_v^T`, `dn_v dn_b^T` of `d_v(dP_b)` are needed,
      !   as the `a_v (dn_b . Au)`, `a_b (dn_v . Au)` pair
      call switched_eigenvalue_curvature(state%n_surf, state%u_switch, state%A_mat, &
                                         dstate_b%dA_mat, dstate_v%dA_mat, ddA, &
                                         res_b%dn_surf, res_v%dn_surf, dres%dn_surf, &
                                         ddlambda_curv)

      ! Eigenvector terms of `d_v(u . dM_b u)` coincide for symmetric `dM_b`, not
      ! merely `M`: the seed precondition
      ddlambda_switch = 2.0_wp*dot_product(dstate_v%du_switch, dstate_b%dM_u) &
                        + ddlambda_curv

      !* ====================== Lifted tangents and the Jacobian ====================== *!

      ! Sphere tangent frame is rigid at every order, dt1 = dt2 = 0
      ddtau1(1) = dot_product(ddq1, state%t1_vec)
      ddtau1(2) = dot_product(ddq2, state%t1_vec)
      ddtau2(1) = dot_product(ddq1, state%t2_vec)
      ddtau2(2) = dot_product(ddq2, state%t2_vec)

      ddw1(1) = ddBinv11*state%tau1(1) + dstate_b%dBinv11*dstate_v%dtau1(1) &
                + dstate_v%dBinv11*dstate_b%dtau1(1) + state%Binv11*ddtau1(1) &
                + ddBinv12*state%tau1(2) + dstate_b%dBinv12*dstate_v%dtau1(2) &
                + dstate_v%dBinv12*dstate_b%dtau1(2) + state%Binv12*ddtau1(2)
      ddw1(2) = ddBinv12*state%tau1(1) + dstate_b%dBinv12*dstate_v%dtau1(1) &
                + dstate_v%dBinv12*dstate_b%dtau1(1) + state%Binv12*ddtau1(1) &
                + ddBinv22*state%tau1(2) + dstate_b%dBinv22*dstate_v%dtau1(2) &
                + dstate_v%dBinv22*dstate_b%dtau1(2) + state%Binv22*ddtau1(2)
      ddw2(1) = ddBinv11*state%tau2(1) + dstate_b%dBinv11*dstate_v%dtau2(1) &
                + dstate_v%dBinv11*dstate_b%dtau2(1) + state%Binv11*ddtau2(1) &
                + ddBinv12*state%tau2(2) + dstate_b%dBinv12*dstate_v%dtau2(2) &
                + dstate_v%dBinv12*dstate_b%dtau2(2) + state%Binv12*ddtau2(2)
      ddw2(2) = ddBinv12*state%tau2(1) + dstate_b%dBinv12*dstate_v%dtau2(1) &
                + dstate_v%dBinv12*dstate_b%dtau2(1) + state%Binv12*ddtau2(1) &
                + ddBinv22*state%tau2(2) + dstate_b%dBinv22*dstate_v%dtau2(2) &
                + dstate_v%dBinv22*dstate_b%dtau2(2) + state%Binv22*ddtau2(2)

      ! No `dalpha` term, see [[drop_seed_input_tangent_type]]
      ddy1 = state%alpha_coeff*(ddw1(1)*state%q1 + dstate_b%dw1(1)*dstate_v%dq1 &
                                + dstate_v%dw1(1)*dstate_b%dq1 + state%w1(1)*ddq1 &
                                + ddw1(2)*state%q2 + dstate_b%dw1(2)*dstate_v%dq2 &
                                + dstate_v%dw1(2)*dstate_b%dq2 + state%w1(2)*ddq2)
      ddy2 = state%alpha_coeff*(ddw2(1)*state%q1 + dstate_b%dw2(1)*dstate_v%dq1 &
                                + dstate_v%dw2(1)*dstate_b%dq1 + state%w2(1)*ddq1 &
                                + ddw2(2)*state%q2 + dstate_b%dw2(2)*dstate_v%dq2 &
                                + dstate_v%dw2(2)*dstate_b%dq2 + state%w2(2)*ddq2)

      ddcross(1) = ddy1(2)*state%y2(3) - ddy1(3)*state%y2(2) &
                   + dstate_b%dy1(2)*dstate_v%dy2(3) - dstate_b%dy1(3)*dstate_v%dy2(2) &
                   + dstate_v%dy1(2)*dstate_b%dy2(3) - dstate_v%dy1(3)*dstate_b%dy2(2) &
                   + state%y1(2)*ddy2(3) - state%y1(3)*ddy2(2)
      ddcross(2) = ddy1(3)*state%y2(1) - ddy1(1)*state%y2(3) &
                   + dstate_b%dy1(3)*dstate_v%dy2(1) - dstate_b%dy1(1)*dstate_v%dy2(3) &
                   + dstate_v%dy1(3)*dstate_b%dy2(1) - dstate_v%dy1(1)*dstate_b%dy2(3) &
                   + state%y1(3)*ddy2(1) - state%y1(1)*ddy2(3)
      ddcross(3) = ddy1(1)*state%y2(2) - ddy1(2)*state%y2(1) &
                   + dstate_b%dy1(1)*dstate_v%dy2(2) - dstate_b%dy1(2)*dstate_v%dy2(1) &
                   + dstate_v%dy1(1)*dstate_b%dy2(2) - dstate_v%dy1(2)*dstate_b%dy2(1) &
                   + state%y1(1)*ddy2(2) - state%y1(2)*ddy2(1)

      cross_dot_b = dot_product(state%cross_vec, dstate_b%dcross_vec)
      dres%dJ = (dot_product(dstate_v%dcross_vec, dstate_b%dcross_vec) &
                 + dot_product(state%cross_vec, ddcross))*state%inv_J &
                + cross_dot_b*dstate_v%dinv_J

      !* ============================ Switching and weights =========================== *!

      dres%dw_f = dstate_v%df_foc_f0*state%f_crit_dS*res_b%d_gnorm &
                  + state%f_foc_f0*dstate_v%df_crit_dS*res_b%d_gnorm &
                  + state%f_foc_f0*state%f_crit_dS*dres%d_gnorm &
                  + dstate_v%df_crit0*state%f_foc_dS*dstate_b%dlambda_switch &
                  + state%f_crit0*dstate_v%df_foc_dS*dstate_b%dlambda_switch &
                  + state%f_crit0*state%f_foc_dS*ddlambda_switch

      ! `dw_pre` of [[apply_seed]] is not stored; rebuilt from `res_b`
      dw_pre_b = state%anchor_wleb0*state%w_f0*res_b%dJ &
                 + state%anchor_wleb0*state%cpjac_scal0*res_b%dw_f
      ddw_pre = dinp_v%danchor_wleb0*state%w_f0*res_b%dJ &
                + state%anchor_wleb0*dinp_v%dw_f0*res_b%dJ &
                + state%anchor_wleb0*state%w_f0*dres%dJ &
                + dinp_v%danchor_wleb0*state%cpjac_scal0*res_b%dw_f &
                + state%anchor_wleb0*dinp_v%dcpjac_scal0*res_b%dw_f &
                + state%anchor_wleb0*state%cpjac_scal0*dres%dw_f
      dres%dwleb = dinp_v%dwbranch*state%wleb_prune_factor*dw_pre_b &
                   + state%wbranch*dstate_v%dwleb_prune_factor*dw_pre_b &
                   + state%wbranch*state%wleb_prune_factor*ddw_pre

      if (state%wleb > seed_weight_tol) then
         dres%dxi = -0.5_wp*(dinp_v%dxi0*res_b%dwleb + state%xi0*dres%dwleb)/state%wleb &
                    + 0.5_wp*state%xi0*res_b%dwleb*dinp_v%dwleb/(state%wleb*state%wleb)
      else
         dres%dxi = 0.0_wp
      end if

      ! `want_curvature` is frozen; default initialisation of `dres` zeroes `dk1`/`dk2`
      if (.not. state%want_curvature) return

      !* ============================ Curvature invariants ============================ *!

      associate (H => state%lsf2_rr, dH_b => res_b%dH, ddH => dres%dH)
         ! d(q_a . H q_b) = dq_a . H q_b + q_a . d(H q_b)
         ! - `d(H q_b)` is the stored `dstate_b%dHq_b`, no `q_a^T dH_b q_b` matvec
         ! - `ddN_ab` relies on symmetric `H`, `dH_b`, as `dN12` and `dB12` of
         !   [[apply_seed]] do
         dHq1_b = matmul(dH_b, state%q1)
         dHq2_b = matmul(dH_b, state%q2)

         dN11_b = dot_product(dstate_b%dq1, state%Hq1) &
                  + dot_product(state%q1, dstate_b%dHq1)
         dN12_b = dot_product(dstate_b%dq1, state%Hq2) &
                  + dot_product(state%q1, dstate_b%dHq2)
         dN22_b = dot_product(dstate_b%dq2, state%Hq2) &
                  + dot_product(state%q2, dstate_b%dHq2)

         ! d_v of the three-term form of [[apply_seed]], frame terms folded into `dHq`
         !   ddN_ab = ddq_a . Hq_b + ddq_b . Hq_a
         !          + dq_a^b . dHq_b^v + dq_b^b . dHq_a^v
         !          + dq_a^v . (dH_b q_b) + dq_b^v . (dH_b q_a)
         !          + q_a . (ddH q_b)
         ddN11 = 2.0_wp*dot_product(ddq1, state%Hq1) &
                 + 2.0_wp*dot_product(dstate_b%dq1, dstate_v%dHq1) &
                 + 2.0_wp*dot_product(dstate_v%dq1, dHq1_b) &
                 + dot_product(state%q1, matmul(ddH, state%q1))
         ddN12 = dot_product(ddq1, state%Hq2) + dot_product(ddq2, state%Hq1) &
                 + dot_product(dstate_b%dq1, dstate_v%dHq2) &
                 + dot_product(dstate_b%dq2, dstate_v%dHq1) &
                 + dot_product(dstate_v%dq1, dHq2_b) &
                 + dot_product(dstate_v%dq2, dHq1_b) &
                 + dot_product(state%q1, matmul(ddH, state%q2))
         ddN22 = 2.0_wp*dot_product(ddq2, state%Hq2) &
                 + 2.0_wp*dot_product(dstate_b%dq2, dstate_v%dHq2) &
                 + 2.0_wp*dot_product(dstate_v%dq2, dHq2_b) &
                 + dot_product(state%q2, matmul(ddH, state%q2))

         ! d_v of `dS_ab = dN_ab/|g| - S_ab d|g|/|g|`
         ddS11 = ddN11/state%g_norm &
                 - dN11_b*res_v%d_gnorm/state%g_norm_sq &
                 - dstate_v%dS11*res_b%d_gnorm/state%g_norm &
                 - state%S11*dres%d_gnorm/state%g_norm &
                 + state%S11*res_b%d_gnorm*res_v%d_gnorm/state%g_norm_sq
         ddS12 = ddN12/state%g_norm &
                 - dN12_b*res_v%d_gnorm/state%g_norm_sq &
                 - dstate_v%dS12*res_b%d_gnorm/state%g_norm &
                 - state%S12*dres%d_gnorm/state%g_norm &
                 + state%S12*res_b%d_gnorm*res_v%d_gnorm/state%g_norm_sq
         ddS22 = ddN22/state%g_norm &
                 - dN22_b*res_v%d_gnorm/state%g_norm_sq &
                 - dstate_v%dS22*res_b%d_gnorm/state%g_norm &
                 - state%S22*dres%d_gnorm/state%g_norm &
                 + state%S22*res_b%d_gnorm*res_v%d_gnorm/state%g_norm_sq

         ddT = ddS11 + ddS22
         ddhalf_diff = 0.5_wp*(ddS11 - ddS22)

         ! d_v of `disc d(disc) = half_diff d(half_diff) + S12 dS12`, solved for `dd_disc`
         ! - `- d_disc^v d_disc^b` comes off the left-hand side
         if (state%disc_curv > seed_curv_disc_guard) then
            dd_disc = (dstate_v%dhalf_diff*dstate_b%dhalf_diff &
                       + state%half_diff*ddhalf_diff &
                       + dstate_v%dS12*dstate_b%dS12 &
                       + state%S12*ddS12 &
                       - dstate_v%ddisc_curv*dstate_b%ddisc_curv)/state%disc_curv
         else
            dd_disc = 0.0_wp
         end if
         dres%dk1 = 0.5_wp*ddT + dd_disc
         dres%dk2 = 0.5_wp*ddT - dd_disc
      end associate

   end subroutine apply_seed_tangent

   !* ====================== Contiguous anchor-group iterator ====================== *!

   !> Advance to the next contiguous anchor group of the branch grid
   !>
   !> - Skip points while `branch_count <= 1`, then extend while `anchor_id` equals
   !>   that of the group's first point
   !> - Only the first point is gated on `branch_count`; later points of the same
   !>   `anchor_id` run join whatever their count
   !> - Both halves of that rule are load bearing
   !> - Contiguity comes from the stable `counting_argsort` of `projection.f90`,
   !>   not checked here; a split run is seen as two groups
   !> - Grid exhausted: `found` false, `first`/`last` an empty range
   !> - Caller starts `cursor` at 1 and never touches it inside its loop
   !> - Subroutine, not function: pure callers need the non-`intent(in)` `cursor`
   !>
   !> @param[in]     branch_count  branches per grid point (ngrid), sets the extent
   !> @param[in]     anchor_id     anchor group id per grid point (ngrid)
   !> @param[in,out] cursor        search position, advanced past the returned group
   !> @param[out]    first         first grid point of the group
   !> @param[out]    last          last grid point of the group
   !> @param[out]    found         `.false.` when no further group exists
   pure subroutine next_branch_group(branch_count, anchor_id, cursor, first, last, found)
      !> Branch bookkeeping per grid point
      integer, intent(in) :: branch_count(:), anchor_id(:)
      !> Search position, advanced past the returned group
      integer, intent(inout) :: cursor
      !> Bounds of the group
      integer, intent(out) :: first, last
      !> Whether a group was found
      logical, intent(out) :: found

      !> Grid extent
      integer :: ngrid

      found = .false.
      first = 1
      last = 0
      ngrid = size(branch_count)

      ! Skip singletons; only the group's head is gated
      do while (cursor <= ngrid)
         if (branch_count(cursor) > 1) exit
         cursor = cursor + 1
      end do
      if (cursor > ngrid) return

      ! Extend the group while anchor_id stays the same
      first = cursor
      last = first
      do while (last < ngrid)
         if (anchor_id(last + 1) /= anchor_id(first)) exit
         last = last + 1
      end do

      cursor = last + 1
      found = .true.

   end subroutine next_branch_group

   !> Largest number of grid points in one contiguous anchor group
   !>
   !> - Scratch bound for the `(:, nbranch)` group buffers of the branch post-pass
   !>   in `forward.f90` and `branch_stage` in `tangent_forward.f90`
   !> - Buffers are indexed by run length: `maxval(branch_count)` is the wrong bound
   !> - Taken from the [[next_branch_group]] walk the callers loop over: no group
   !>   they see is wider
   !> - Zero when no multi-branch group exists
   !>
   !> @param[in] branch_count  branches per grid point (ngrid), sets the extent
   !> @param[in] anchor_id     anchor group id per grid point (ngrid)
   pure function max_branch_group_size(branch_count, anchor_id) result(nmax)
      !> Branch bookkeeping per grid point
      integer, intent(in) :: branch_count(:), anchor_id(:)
      !> Maximum group extent
      integer :: nmax

      !> Group walk bookkeeping
      integer :: cursor, ifirst, ilast
      !> Whether a further group exists
      logical :: found

      nmax = 0
      cursor = 1
      do
         call next_branch_group(branch_count, anchor_id, cursor, ifirst, ilast, found)
         if (.not. found) exit
         nmax = max(nmax, ilast - ifirst + 1)
      end do

   end function max_branch_group_size

   !> Reverse pass over the branch-weight softmax
   !>
   !> - Within an anchor group `wleb_m = base_m * p_m`; seed loop handles `d(base_m)`
   !> - Converts the width-induced adjoint `dL/dp_m` into `dL/dPhi_m`, which the
   !>   seed loop couples to the point motion
   !> - Groups from [[next_branch_group]]; points with `branch_count <= 1` stay zero
   !> - TODO: serial, should/could be parallelized
   !>
   !> @param[in]  branch_count    number of branches per grid point (ngrid)
   !> @param[in]  anchor_id       anchor group id per grid point (ngrid)
   !> @param[in]  wbranch         softmax branch weight per grid point (ngrid)
   !> @param[in]  wleb            Lebedev weight per grid point (ngrid)
   !> @param[in]  xi0             Gaussian width per grid point (ngrid)
   !> @param[in]  sigma_phi       softmax temperature
   !> @param[in]  w_xi            effective Gaussian-width adjoint (ngrid)
   !> @param[out] branch_phi_adj  adjoint of the branch objective Phi (ngrid)
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
      integer :: ngrid, igroup_cursor, igroup_start, igroup_end, group_size
      integer :: m_branch, im_grid
      !> Whether a further anchor group exists
      logical :: have_group
      !> Weight-adjoint scratch
      real(wp) :: adj_branch, mean_adj_branch

      branch_phi_adj = 0.0_wp
      ngrid = size(branch_count)
      if (ngrid <= 0) return
      if (.not. any(branch_count > 1)) return
      if (sigma_phi <= seed_weight_tol) return

      igroup_cursor = 1
      do
         call next_branch_group(branch_count, anchor_id, igroup_cursor, &
                                igroup_start, igroup_end, have_group)
         if (.not. have_group) exit
         group_size = igroup_end - igroup_start + 1

         mean_adj_branch = 0.0_wp
         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            call branch_point_adjoint(w_xi(im_grid), wleb(im_grid), xi0(im_grid), &
                                      wbranch(im_grid), adj_branch)
            mean_adj_branch = mean_adj_branch + wbranch(im_grid)*adj_branch
            branch_phi_adj(im_grid) = adj_branch
         end do

         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            branch_phi_adj(im_grid) = -wbranch(im_grid) &
                                      *(branch_phi_adj(im_grid) - mean_adj_branch)/sigma_phi
         end do
      end do

   end subroutine compute_branch_phi_adj

   !> Branch adjoint of one point, optionally with its directional tangent
   !>
   !> - Shared by [[compute_branch_phi_adj]] (primal only) and
   !>   [[branch_phi_adj_tangent]]: one copy, same compiled statements, no ulp drift
   !> - Raw per-point adjoint, before the group reduction
   !>   `-wbranch*(adj - mean_adj)/sigma_phi` of [[compute_branch_phi_adj]]; not the
   !>   stored `branch_phi_adj`
   !> - `adj_wleb * factor` stays factored, not collapsed to
   !>   `-0.5 w_xi xi0/wbranch`, to match the derivation of the reverse pass
   !> - Tangent arguments all-or-none: `dw_xi`, `dwleb`, `dxi0`, `dwbranch`, `dadj_branch`
   !> - Tangent takes the primal's branch on the primal's condition, hard zero on
   !>   the else, no threshold of its own
   !>
   !> @param[in]  w_xi         folded width adjoint at the point
   !> @param[in]  wleb         Lebedev weight at the point
   !> @param[in]  xi0          Gaussian width at the point
   !> @param[in]  wbranch      softmax branch weight at the point
   !> @param[out] adj_branch   raw branch adjoint, zero when the gate is shut
   !> @param[in]  dw_xi        tangent of the folded width adjoint
   !> @param[in]  dwleb        tangent of the Lebedev weight
   !> @param[in]  dxi0         tangent of the Gaussian width
   !> @param[in]  dwbranch     tangent of the branch weight
   !> @param[out] dadj_branch  tangent of the branch adjoint, zero on the same gate
   pure subroutine branch_point_adjoint(w_xi, wleb, xi0, wbranch, adj_branch, &
                                        dw_xi, dwleb, dxi0, dwbranch, dadj_branch)
      !> Point values
      real(wp), intent(in) :: w_xi, wleb, xi0, wbranch
      !> Branch adjoint
      real(wp), intent(out) :: adj_branch
      !> Point tangents, present together with `dadj_branch` or not at all
      real(wp), intent(in), optional :: dw_xi, dwleb, dxi0, dwbranch
      !> Tangent of the branch adjoint; its presence selects the tangent half
      real(wp), intent(out), optional :: dadj_branch

      !> Weight-adjoint scratch and its tangent
      real(wp) :: adj_wleb, dadj_wleb, factor_m, dfactor_m

      adj_branch = 0.0_wp
      if (present(dadj_branch)) dadj_branch = 0.0_wp
      if (abs(w_xi) > seed_weight_tol &
          .and. wleb > seed_weight_tol &
          .and. wbranch > tiny(1.0_wp)) then
         adj_wleb = -0.5_wp*w_xi*xi0/wleb
         factor_m = wleb/wbranch
         adj_branch = adj_wleb*factor_m
         if (present(dadj_branch)) then
            dadj_wleb = -0.5_wp*(dw_xi*xi0 + w_xi*dxi0)/wleb &
                        + 0.5_wp*w_xi*xi0*dwleb/(wleb*wleb)
            dfactor_m = dwleb/wbranch - wleb*dwbranch/(wbranch*wbranch)
            dadj_branch = dadj_wleb*factor_m + adj_wleb*dfactor_m
         end if
      end if

   end subroutine branch_point_adjoint

   !> Adjoint contribution of one seed
   !>
   !> - Shared contraction of [[seed_jet_basis]], [[seed_anchor]] and the primal
   !>   seed loop of the fixed-adjoint Hessian: one copy, same compiled statements,
   !>   no ulp drift
   !> - Accumulation order is load bearing: position, width, branch, curvature
   !> - No switching term: `f_i` is an anchor-only iSwiG overlap, unchanged by a
   !>   level-set perturbation at fixed nuclei; anchor motion goes through the
   !>   caller's switching channel
   !> - Branch term gated on the stored `eff%branch_phi_adj`, allocated to zero by
   !>   [[prepare_surface_weights]]; [[seed_contribution_tangent]] uses the same gate
   !> - `branch_shift`: rigid-motion piece of an anchor seed, `-phi1_r(iaxis)` from
   !>   `phi = 0.5*alpha*|r - anchor|^2` at fixed `r`; omitted for a field seed
   !>
   !> @param[in] eff           folded surface adjoints
   !> @param[in] igrid         grid point
   !> @param[in] w_xyz_pt      effective position adjoint
   !> @param[in] dr            induced point motion of the seed
   !> @param[in] res           linear response of the seed
   !> @param[in] phi1_r        objective gradient at the projected point
   !> @param[in] branch_shift  rigid-motion shift of an anchor seed, omitted otherwise
   pure function seed_contribution(eff, igrid, w_xyz_pt, dr, res, phi1_r, branch_shift) &
      result(contribution)
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Effective position adjoint
      real(wp), intent(in) :: w_xyz_pt(3)
      !> Induced point motion
      real(wp), intent(in) :: dr(3)
      !> Linear response
      type(drop_seed_result_type), intent(in) :: res
      !> Objective gradient at the projected point
      real(wp), intent(in) :: phi1_r(3)
      !> Rigid-motion shift, present for an anchor seed only
      real(wp), intent(in), optional :: branch_shift
      !> Adjoint contribution
      real(wp) :: contribution

      !> Objective motion the branch adjoint contracts against
      real(wp) :: branch_dphi

      contribution = dot_product(w_xyz_pt, dr) + eff%w_xi(igrid)*res%dxi
      if (abs(eff%branch_phi_adj(igrid)) > seed_weight_tol) then
         branch_dphi = dot_product(phi1_r, dr)
         if (present(branch_shift)) branch_dphi = branch_dphi - branch_shift
         contribution = contribution + eff%branch_phi_adj(igrid)*branch_dphi
      end if
      if (eff%have_wk) then
         contribution = contribution + eff%w_k1(igrid)*res%dk1 + eff%w_k2(igrid)*res%dk2
      end if

   end function seed_contribution

   !> Add one host jet tangent onto a directional tangent of the jet
   !>
   !> - Host partial tangent at the fixed point, on top of the level set's own
   !>   nuclear tangent
   !> - Third order only where the second-order chain asks for it
   !>
   !> @param[in]     jet  host jet tangent, `drop_n_host_jet` entries
   !> @param[in,out] dv0  tangent of the value
   !> @param[in,out] dv1  tangent of the gradient
   !> @param[in,out] dv2  tangent of the Hessian
   !> @param[in,out] dv3  tangent of the third derivative, optional
   pure subroutine add_host_jet(jet, dv0, dv1, dv2, dv3)
      !> Host jet tangent
      real(wp), intent(in) :: jet(:)
      !> Tangents of the jet along the direction
      real(wp), intent(inout) :: dv0, dv1(3), dv2(3, 3)
      !> Tangent of the third derivative along the direction
      real(wp), intent(inout), optional :: dv3(3, 3, 3)

      dv0 = dv0 + jet(1)
      dv1 = dv1 + jet(2:4)
      dv2 = dv2 + reshape(jet(5:13), [3, 3])
      if (present(dv3)) dv3 = dv3 + reshape(jet(14:40), [3, 3, 3])
   end subroutine add_host_jet

   !> Diagnostic message for a degeneracy status
   !>
   !> - Callers prepend their own context and append the offending grid point
   !>
   !> @param[in] status  one of the `seed_state_*` codes
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
