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
!> TODO: Refactor types so that a clear intent in/out structure is visible
!>       this way, we will not have the risk of corruption during the run
module moist_cavity_drop_derivatives_kernel
   use mctc_env_accuracy, only: wp
   use moist_math_linalg, only: setup_tangent_frame, eig_2x2_symmetric, &
      & eig_2x2_offdiag_tol
   use moist_cavity_drop_switching, only: moist_cavity_drop_swif_type

   implicit none(type, external)
   private

   public :: drop_seed_state_type, drop_seed_result_type, drop_surface_weights_type
   public :: drop_seed_state_tangent_type, drop_seed_input_tangent_type
   public :: drop_seed_result_tangent_type
   public :: build_seed_state, apply_seed, apply_seed_tangent, compute_branch_phi_adj
   public :: switched_eigenvalue_response, switched_eigenvalue_curvature
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
   !> Below this eigenvalue gap of `B` the switched eigenvector is treated as
   !> degenerate. `d(u_switch)` scales as `1/gap`, and at `gap = 0` the
   !> eigenvector itself is arbitrary, so no derivative exists to compute
   real(wp), parameter :: seed_eig_gap_guard = 1.0e-10_wp

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
      !> Anchor position and its owner sphere center
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
      !> Whether `eig_2x2_symmetric` built `vmin_B` from `[B12, lambda - B11]`
      !> rather than falling back to a canonical basis vector; `vmin_norm` is
      !> meaningful only in the first case
      logical :: vmin_offdiag = .false.
      !> Switched (smaller) eigenvalue of `B` and the eigenvalue gap
      !> `lambda_max - lambda_min = sqrt(disc)`
      real(wp) :: lambda_switch = 0.0_wp, sqrt_disc_B = 0.0_wp
      !> Normalised 2D eigenvector of `B` and the norm it was divided by
      real(wp) :: vmin_B(2) = 0.0_wp, vmin_norm = 0.0_wp
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
      !> Switching-function curvatures, consumed by the second-order chain
      real(wp) :: f_crit_d2S = 0.0_wp, f_foc_d2S = 0.0_wp
      !> Lebedev-weight pruning chain factor
      real(wp) :: wleb_prune_factor = 1.0_wp
      !> Whether Lebedev-weight pruning was active when the state was built
      logical :: use_wleb_prune = .false.
      !> Pre-pruning weight product, and the pruning-switch slope and curvature
      !> at `|w_pre_i|`; all zero when pruning is off
      real(wp) :: w_pre_i = 0.0_wp, f_wleb_ds = 0.0_wp, f_wleb_d2S = 0.0_wp
      !> Level-set Hessian applied to the surface tangent frame, and the shape
      !> operator it induces there (only when `want_curvature`)
      real(wp) :: Hq1(3) = 0.0_wp, Hq2(3) = 0.0_wp
      real(wp) :: S11 = 0.0_wp, S12 = 0.0_wp, S22 = 0.0_wp
      !> Shape-operator invariants: trace, half-difference and the eigenvalue
      !> gap `|k1 - k2|/2`, the last as a sum of squares (only when
      !> `want_curvature`)
      real(wp) :: T_curv = 0.0_wp, KM_curv = 0.0_wp
      real(wp) :: half_diff = 0.0_wp, disc_curv = 0.0_wp

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

   !> Tangent of the derived block of [[drop_seed_state_type]] along one seed
   !>
   !> [[apply_seed]] forms this entire chain as local scratch on its way to
   !> [[drop_seed_result_type]] and then discards it. The second-order path
   !> needs it, so [[apply_seed]] hands it back through an optional argument
   !> rather than through a second routine: two copies of the same expression
   !> are free to contract differently under `-ffp-contract=fast`, and the
   !> shipped first-order path must stay bit-for-bit unchanged.
   !>
   !> One seed direction per instance. Only the derived fields that
   !> [[apply_seed]] reads and that are not frozen appear here:
   !>
   !>   * `dn_surf`, `d_gnorm` and `dH` live on [[drop_seed_result_type]] and
   !>     are not duplicated
   !>   * `dt1_vec` and `dt2_vec` are absent because they vanish identically.
   !>     The sphere tangent frame is rigid (see the comment on `dtau1` in
   !>     [[apply_seed]]): the anchor rides its owner sphere, so
   !>     `anchor - owner_xyz = R_own * u_leb` is invariant under every nuclear
   !>     direction, at every order. They become nonzero only once the radii
   !>     themselves are geometry dependent, which is not implemented
   !>   * `min_axis` and `want_curvature` are frozen discrete choices
   !>
   !> Every component is default initialised, which is load bearing rather than
   !> stylistic: the `want_curvature` early return in [[apply_seed]] leaves the
   !> curvature block untouched and relies on `intent(out)` default
   !> initialisation to zero it.
   type :: drop_seed_state_tangent_type

      !* -------------------------- KKT matrix and frame -------------------------- *!

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
      !> The vector `dM u` of the basis-invariant tangent operator `M = P A P`,
      !> taken at the *base* eigenvector: this is `(dM) u`, not `d(M u)`.
      !> Stored as the contracted vector rather than as the matrix because that
      !> is the only thing any consumer wants -- [[apply_seed_tangent]] pairs it
      !> with `du_switch` and nothing reads `dM` itself. It depends on the seed
      !> alone, so storing it also keeps the `b` chain out of the direction loop,
      !> which runs once per direction *pair*
      real(wp) :: dM_u(3) = 0.0_wp
      !> Sensitivity of the switched eigenvector, in the `B` basis and lifted
      real(wp) :: dvmin_B(2) = 0.0_wp, du_switch(3) = 0.0_wp

      !* -------------------- Lifted tangents and the Jacobian -------------------- *!

      !> Sensitivity of the sphere frame projected onto the surface frame
      real(wp) :: dtau1(2) = 0.0_wp, dtau2(2) = 0.0_wp
      !> Sensitivity of the `B^-1` images of those projections
      real(wp) :: dw1(2) = 0.0_wp, dw2(2) = 0.0_wp
      !> Sensitivity of the lifted tangent vectors and of their cross product
      real(wp) :: dy1(3) = 0.0_wp, dy2(3) = 0.0_wp, dcross_vec(3) = 0.0_wp
      !> Sensitivity of `1/J`
      real(wp) :: dinv_J = 0.0_wp

      !* ------------------------------ Gram-Schmidt ------------------------------ *!

      !> Sensitivity of the Gram-Schmidt data behind `q1`
      real(wp) :: dn_dot_q1 = 0.0_wp, dproj_surf = 0.0_wp, dv_norm_surf = 0.0_wp

      !* ------------------------- Switching and weights -------------------------- *!

      !> Sensitivity of the switching values and of their slopes
      real(wp) :: df_crit0 = 0.0_wp, df_crit_dS = 0.0_wp
      real(wp) :: df_foc_f0 = 0.0_wp, df_foc_dS = 0.0_wp
      !> Sensitivity of the Lebedev-weight pruning factor; identically zero when
      !> pruning is off
      real(wp) :: dwleb_prune_factor = 0.0_wp
      !> Sensitivity of `|grad S|^2`
      real(wp) :: dg_norm_sq = 0.0_wp

      !* -------------------------- Curvature invariants -------------------------- *!

      !> Sensitivity of the Hessian applied to the tangent frame, as the
      !> **total** `d(H q_a) = dH q_a + H dq_a`, matching what the old `dHn`
      !> carried. Stored rather than rebuilt for the same reason as `dM_u`: it
      !> depends on the seed alone, while [[apply_seed_tangent]] wants it once
      !> per direction *pair*
      real(wp) :: dHq1(3) = 0.0_wp, dHq2(3) = 0.0_wp
      !> Sensitivity of the shape-operator entries
      real(wp) :: dS11 = 0.0_wp, dS12 = 0.0_wp, dS22 = 0.0_wp
      !> Sensitivity of the shape-operator invariants
      real(wp) :: dT_curv = 0.0_wp, dKM_curv = 0.0_wp
      real(wp) :: dhalf_diff = 0.0_wp, ddisc_curv = 0.0_wp

   end type drop_seed_state_tangent_type

   !> v-direction tangent of the `Inputs` block of [[drop_seed_state_type]]
   !>
   !> Passed explicitly rather than reconstructed inside
   !> [[apply_seed_tangent]]. The driver knows that, for a physical direction,
   !> `dlsf1_r` is `res_v%dg`, `dcpjac_scal0` is `res_v%dJ` and so on -- but
   !> encoding those identities in the kernel would make it untestable along an
   !> arbitrary path, and would put the physics in the wrong place. Only the
   !> fields [[apply_seed]] actually reads appear
   !>
   !> `alpha_coeff` is the one real-valued input [[apply_seed]] reads that has no
   !> entry here, and that omission is a decision rather than an oversight. It is
   !> `param%phi_alpha`, a fixed parameter, and [[apply_seed]] carries no seed
   !> channel for it, so no consistent `(dstate_v, res_v)` pair can ever give it
   !> a nonzero tangent. Re-adding the field alone would not make a
   !> geometry-dependent `alpha` work either: `ddA` in [[apply_seed_tangent]]
   !> carries no `dalpha I` term, so the `y` chain would be right while the `A`
   !> tangent stayed wrong. Both have to move together
   type :: drop_seed_input_tangent_type
      !> Tangent of the level-set jet at the projected point
      real(wp) :: dlsf1_r(3) = 0.0_wp, dlsf2_rr(3, 3) = 0.0_wp
      real(wp) :: dlsf3_rrr(3, 3, 3) = 0.0_wp
      !> Tangent of the multiplier
      real(wp) :: dlambda_val = 0.0_wp
      !> Tangent of the grid-level weight scalars
      real(wp) :: danchor_wleb0 = 0.0_wp, dcpjac_scal0 = 0.0_wp, dw_f0 = 0.0_wp
      real(wp) :: dwbranch = 0.0_wp, dwleb = 0.0_wp, dxi0 = 0.0_wp
   end type drop_seed_input_tangent_type

   !> Tangent of [[drop_seed_result_type]] along a second direction
   !>
   !> Field names mirror [[drop_seed_result_type]] exactly: `dres%dg` is the
   !> directional derivative of `res%dg`, not a new quantity. Every component is
   !> default initialised, for the same reason as in
   !> [[drop_seed_state_tangent_type]]: the `want_curvature` early return in
   !> [[apply_seed_tangent]] leaves the curvature pair untouched
   type :: drop_seed_result_tangent_type
      real(wp) :: dg(3) = 0.0_wp, dH(3, 3) = 0.0_wp
      real(wp) :: dn_surf(3) = 0.0_wp, d_gnorm = 0.0_wp
      real(wp) :: dJ = 0.0_wp, dw_f = 0.0_wp
      real(wp) :: dwleb = 0.0_wp, dxi = 0.0_wp
      real(wp) :: dk1 = 0.0_wp, dk2 = 0.0_wp
   end type drop_seed_result_tangent_type

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
      ! Keep the eigen data the tangent needs; `vmin_norm` is the normalisation
      ! `eig_2x2_symmetric` divided by, recomputed here rather than returned
      state%lambda_switch = lambda_switch
      ! The eigenvalue gap as `hypot(B11 - B22, 2 B12)`, not the algebraically
      ! equal `beta_max - lambda_switch`. That form inherits
      ! `eig_2x2_symmetric`'s `disc = trace^2 - 4 det`, a cancellation that
      ! loses the entire gap once it drops below `sqrt(eps)*|trace|`: at
      ! `B = [[1.3, 1e-8], [1e-8, 1.3]]` it returns exactly zero while `|B12|`
      ! is six orders above that routine's own diagonal-branch threshold, so
      ! the branch test below is no protection for the division on it
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

      ! Shape operator of the level set in the surface tangent frame Q = [q1, q2]
      ! already built above for the closest-point Jacobian:
      !   S_ab = q_a^T H q_b / |g|,   a, b in {1, 2}
      ! whose eigenvalues are the principal curvatures k1 >= k2,
      !   k1,k2 = KM +/- sqrt(half_diff^2 + S12^2),
      ! with KM = (S11 + S22)/2 and half_diff = (S11 - S22)/2.
      !
      ! The discriminant is deliberately a sum of squares and not the
      ! algebraically equal invariant form `sqrt(KM^2 - KG)` with
      ! `KG = n^T adj(H) n / |g|^2`. `KM^2 - KG` is `((k1 - k2)/2)^2` written as
      ! a difference of two quantities of size `KM^2`, so it loses the gap
      ! exactly where the gap is small -- and [[apply_seed]] then divides by it.
      ! Same reason `sqrt_disc_B` above is a `hypot` and not `trace^2 - 4 det`,
      ! and same form `properties.f90` uses for the primal.
      !
      ! Both `KM` and `disc` are invariant under a rotation of Q, so any smooth
      ! orthonormal tangent frame gives the right derivative downstream; it need
      ! not be the one `compute_curvature` picks, only differentiated
      ! consistently, which `dq1`/`dq2` are.
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
   !> With `P = I - n n^T` the basis-invariant tangent operator is `M = P A P`,
   !> and the switched eigenvalue responds as `dlambda = u . (dM u)` with `u`
   !> the base eigenvector. Assembling `dM` to reach that one scalar costs four
   !> `spread`s and five matrix-matrix products, and gfortran never inlines
   !> `spread`, so every one of those allocates a heap temporary. Contracting
   !> the quadratic form first leaves two matrix-vector products and four dot
   !> products, and no temporaries at all.
   !>
   !> Three preconditions, all structural rather than incidental:
   !>
   !>   1. `n . u = 0`. `u_switch` is `vmin_B(1) q1 + vmin_B(2) q2` and
   !>      [[setup_tangent_frame]] builds `q1`, `q2` orthogonal to `n`, so this
   !>      holds to roundoff -- which is all a derivative needs
   !>   2. `P u = u`, which follows from 1
   !>   3. `A` is symmetric. `A = alpha I - lambda H` with `H` the level-set
   !>      Hessian, so this one is structural too. It is what the factor of two
   !>      in `dlambda` rests on, through `u . (A n) = n . (A u)`. `dA` is
   !>      *not* required to be symmetric anywhere in this routine: both
   !>      outputs are exact for an arbitrary `dA`
   !>
   !> From 1 and 2, `dP u = -(dn . u) n`, a scalar times the normal. Writing
   !> `a = dn . u`, the identity that replaces the assembly is
   !>
   !>     dlambda = u . (dA u) - 2 a (n . A u)
   !>     dM u    = -(dn (n . Au) + n (dn . Au))
   !>               + (dA u - n (n . dA u))
   !>               - a (An - n (n . An))
   !>
   !> `dM u` is handed back as one vector because the second-order chain
   !> contracts it once: the two eigenvector cross terms of `d_v(u . dM_b u)`
   !> collapse to `2 du_v . (dM_b u)`. That collapse needs `dM_b` symmetric --
   !> not merely `M` -- and `dM = dP A P + P dA P + P A dP` is symmetric only
   !> when `dA` is, which the nine Hessian basis seeds of [[seed_jet_basis]] are
   !> not. The collapse is the caller's, and so is the condition under which it
   !> holds; see the precondition block on [[apply_seed_tangent]].
   !>
   !> @param[in]  n_surf          Outward unit normal
   !> @param[in]  u_switch        Switched eigenvector, tangent to the surface
   !> @param[in]  A_mat           Tangent-restricted KKT matrix, symmetric
   !> @param[in]  dA              Sensitivity of `A`; need not be symmetric
   !> @param[in]  dn              Sensitivity of the normal
   !> @param[out] dlambda_switch  Response of the switched eigenvalue
   !> @param[out] dM_u            The vector `dM u`, for the second-order chain
   pure subroutine switched_eigenvalue_response(n_surf, u_switch, A_mat, dA, dn, &
                                                dlambda_switch, dM_u)
      !> Outward unit normal and the switched eigenvector
      real(wp), intent(in) :: n_surf(3), u_switch(3)
      !> Tangent-restricted KKT matrix, symmetric, and its sensitivity, which
      !> need not be
      real(wp), intent(in) :: A_mat(3, 3), dA(3, 3)
      !> Sensitivity of the normal
      real(wp), intent(in) :: dn(3)
      !> Response of the switched eigenvalue
      real(wp), intent(out) :: dlambda_switch
      !> The vector `dM u`; formed only when the caller asks for it
      real(wp), intent(out), optional :: dM_u(3)

      !> Matrix-vector products the identity is built from
      real(wp) :: Au(3), An(3), dA_u(3)
      !> `dn . u` and `n . Au`, each read more than once
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
   !> Returns the scalar `u . (ddM u)` of `M = P A P`, where `ddX` means
   !> `d_v(d_b X)`. Assembled the obvious way that scalar costs eight `spread`s
   !> and fourteen matrix-matrix products, and the assembly sits in the
   !> innermost `(b, v)` loop of the whole Hessian, where each un-inlined
   !> `spread` is one more heap allocation. Contracted it is five matrix-vector
   !> products and ten dot products.
   !>
   !> Preconditions 1 to 3 of [[switched_eigenvalue_response]] carry over
   !> unchanged: `n . u = 0`, hence `P u = u`, and `A` symmetric. Unlike that
   !> routine, this one needs `dA_b` and `dA_v` symmetric as well. The
   !> `-2 [ a_b (n . dA_v u) + a_v (n . dA_b u) ]` bracket below is four terms
   !> collapsed into two, and the collapse is `u . (dA n) = n . (dA u)`. `ddA`
   !> enters only as `u . (ddA u)`, so of it just the symmetric part is read.
   !>
   !> That is a real precondition, not a formality. `apply_seed` builds
   !> `dA = -dlambda lsf2_rr - lambda dH` with `dH = dlsf2_rr + lsf3_rrr . dr`;
   !> `lsf2_rr` and the `lsf3_rrr` contraction are symmetric, so the entire
   !> antisymmetric part of `dA` is `-lambda` times that of the seed's
   !> `dlsf2_rr`. Symmetric for a physical nuclear direction -- and false for
   !> the nine single-entry Hessian basis seeds of [[seed_jet_basis]]. See the
   !> precondition block on [[apply_seed_tangent]] for what a driver must do.
   !>
   !> The two eigenvector cross terms of `d_v(u . dM_b u)` coincide under the
   !> same condition and collapse into a factor of two. Those belong to the
   !> caller, as `2 du_v . (dM_b u)`; what this routine returns is the remaining
   !> term.
   !>
   !> With `a_b = dn_b . u`, `a_v = dn_v . u` and `c = ddn . u`,
   !>
   !>     u . (ddM u) = -2 [ a_v (dn_b . Au) + a_b (dn_v . Au) + c (n . Au) ]
   !>                   -2 [ a_b (n . dA_v u) + a_v (n . dA_b u) ]
   !>                   +2 a_b a_v (n . An)
   !>                   + u . (ddA u)
   !>
   !> @param[in]  n_surf           Outward unit normal
   !> @param[in]  u_switch         Switched eigenvector, tangent to the surface
   !> @param[in]  A_mat            Tangent-restricted KKT matrix, symmetric
   !> @param[in]  dA_b             Sensitivity of `A` along `b`, required symmetric
   !> @param[in]  dA_v             Sensitivity of `A` along `v`, required symmetric
   !> @param[in]  ddA              Second-order sensitivity of `A`; only its
   !>                              symmetric part is read
   !> @param[in]  dn_b             Sensitivity of the normal along `b`
   !> @param[in]  dn_v             Sensitivity of the normal along `v`
   !> @param[in]  ddn              Second-order sensitivity of the normal
   !> @param[out] ddlambda_switch  The scalar `u . (ddM u)`
   pure subroutine switched_eigenvalue_curvature(n_surf, u_switch, A_mat, &
                                                 dA_b, dA_v, ddA, &
                                                 dn_b, dn_v, ddn, ddlambda_switch)
      !> Outward unit normal and the switched eigenvector
      real(wp), intent(in) :: n_surf(3), u_switch(3)
      !> Tangent-restricted KKT matrix, symmetric
      real(wp), intent(in) :: A_mat(3, 3)
      !> First-order sensitivities of `A`, both required symmetric, and the
      !> second-order one, of which only the symmetric part is read
      real(wp), intent(in) :: dA_b(3, 3), dA_v(3, 3), ddA(3, 3)
      !> First- and second-order sensitivities of the normal
      real(wp), intent(in) :: dn_b(3), dn_v(3), ddn(3)
      !> The scalar `u . (ddM u)`
      real(wp), intent(out) :: ddlambda_switch

      !> Matrix-vector products the identity is built from
      real(wp) :: Au(3), An(3), dA_b_u(3), dA_v_u(3)
      !> `dn_b . u`, `dn_v . u`, `ddn . u` and `n . Au`, each read more than once
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
   !> @param[out] dstate   Tangent of the derived seed state along this seed
   pure subroutine apply_seed(state, dlsf1_r, dlsf2_rr, dr, dlambda, res, dstate)
      !> Per-grid point forward state
      type(drop_seed_state_type), intent(in) :: state
      !> Seed perturbation of the level-set gradient and Hessian
      real(wp), intent(in) :: dlsf1_r(3), dlsf2_rr(3, 3)
      !> Induced motion of the projected point and its multiplier
      real(wp), intent(in) :: dr(3), dlambda
      !> Linear response
      type(drop_seed_result_type), intent(out) :: res
      !> Tangent of the derived block of `state`; every field behind this
      !> argument is computed only when it is present
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
      !> Curvature sensitivities. `dHq1_p`/`dHq2_p` are the *partial* `dH q_a`;
      !> the total `d(H q_a)` that `dstate` carries adds `H dq_a`
      real(wp) :: dHq1_p(3), dHq2_p(3)
      real(wp) :: dN11, dN12, dN22, dS11, dS12, dS22
      real(wp) :: dT, dhalf_diff, d_disc
      !> Line-by-line derivative of the `eig_2x2_symmetric` construction,
      !> needed only for `dstate`
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

      ! `u . (dM u)` in contracted form, with `M = P A P` the basis-invariant
      ! tangent operator; see [[switched_eigenvalue_response]] for the identity
      ! and its three preconditions. The assembly it replaces built `dP` out of
      ! four `spread`s and `dM` out of five matrix-matrix products, none of them
      ! inlined, all of it to reach this one scalar. `dM_u` is the only piece of
      ! `dM` the second-order path reads, so it is formed only when the caller
      ! asks for the state tangent
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
         ! The scratch `dAq1`/`dAq2` above carry the `dA . q` term only, because
         ! the `dq` half rides separately in `dB11`/`dB12`/`dB22`. The state
         ! tangent needs the full product rule
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

         ! d(v_min) by differentiating the shipped `eig_2x2_symmetric`
         ! construction line by line, which reproduces its sign convention
         ! `v_min = [B12, lambda_min - B11]/norm` automatically.
         !
         ! The two guards are independent, and each covers a case the other
         ! does not. Off the off-diagonal branch the primal returns a canonical
         ! basis vector, which is piecewise constant, so the derivative is zero
         ! rather than this formula -- and there `vmin_norm` is not even the
         ! norm of a vector the primal used. A vanishing gap is the separate
         ! degeneracy: the eigenvector is arbitrary at `k1 = k2`, so zero is the
         ! guard's answer, the same one `seed_curv_disc_guard` gives, not a
         ! claim that the derivative is zero
         if (state%vmin_offdiag .and. state%sqrt_disc_B > seed_eig_gap_guard) then
            dtrace_B = dB11 + dB22
            ! `d(disc)` as `d((B11 - B22)^2 + 4 B12^2)`. Algebraically identical
            ! to the `trace^2 - 4 det` form, but that one loses relative
            ! accuracy exactly where the gap is small and this quantity matters
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
         ! Mirror the primal `max(proj_surf, 1e-30)` clamp: once it bites the
         ! norm is constant, exactly as the `dq1` branch above assumes
         if (state%proj_surf > 1.0e-30_wp) then
            dstate%dv_norm_surf = dstate%dproj_surf/(2.0_wp*state%v_norm_surf)
         else
            dstate%dv_norm_surf = 0.0_wp
         end if
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
         ! wleb_prune_factor = S(|w|) + |w| S'(|w|), so with
         ! d|w| = sign(w_pre_i) * dw_pre the chain collapses to one factor
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
         ! d(q_a^T H q_b) = dq_a . H q_b + dq_b . H q_a + q_a^T dH q_b, using the
         ! symmetry of H for the mixed term exactly as `dB12` above does
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

         ! disc^2 = half_diff^2 + S12^2, so
         ! disc d(disc) = half_diff d(half_diff) + S12 dS12 -- both products of
         ! quantities that are themselves O(disc), so the quotient keeps the
         ! relative accuracy `disc` was built with
         if (state%disc_curv > seed_curv_disc_guard) then
            d_disc = (state%half_diff*dhalf_diff + state%S12*dS12)/state%disc_curv
         else
            d_disc = 0.0_wp
         end if
         res%dk1 = 0.5_wp*dT + d_disc
         res%dk2 = 0.5_wp*dT - d_disc

         if (present(dstate)) then
            ! The *total* d(H q_a), the form [[apply_seed_tangent]] contracts
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
   !> [[apply_seed]] maps a seed `b = (dlsf1_r, dlsf2_rr, dr, dlambda)` to a
   !> linear response `res_b` that also depends on the per-point state. This
   !> routine returns `d/dv [ res_b ]`. It is the derivative *of that code*,
   !> obtained by applying the product rule to [[apply_seed]] line by line, and
   !> not a second derivation of the geometry: every `state%X` is replaced by the
   !> v-tangent of `X`, every first-order local by its own v-tangent, and every
   !> `res%dZ` by `dres%dZ`.
   !>
   !> --------------------------------------------------------------------------
   !>
   !> **PRECONDITION -- the seed's `dlsf2_rr` must be symmetric, or the caller
   !> must add the transpose seed before reading anything.**
   !>
   !> Three steps on the way to `dres` collapse a symmetric pair of terms into a
   !> factor of two: the `2 du_v . (dM_b u)` eigenvector term below, the
   !> `-2 [ a_b (n . dA_v u) + a_v (n . dA_b u) ]` bracket inside
   !> [[switched_eigenvalue_curvature]], and the `dBij = dqi . Aqj + qi . dAqj`
   !> restatement that the `ddB` lines differentiate. Each holds only for a
   !> symmetric `dA_b`, and [[apply_seed]] carries the seed's `dlsf2_rr` into
   !> `dA` scaled by `-lambda` and nothing else, so this is a condition on the
   !> *seed*, not on the geometry.
   !>
   !> [[seed_jet_basis]] is the producer that breaks it. Its nine Hessian basis
   !> seeds are single-entry matrices, `dlsf2_rr(iaxis, jaxis) = 1`, so the six
   !> off-diagonal ones are maximally asymmetric. Fed such a seed, this routine
   !> does **not** return the derivative of that seed's response, and not by a
   !> little: the per-seed error is of order 100 % and flips signs. Measured on
   !> `dres%dw_f`, analytic `-2.81e-3` against a finite difference of `+1.50e-3`.
   !>
   !> What a driver must do: never read a single off-diagonal Hessian seed's
   !> `dres`. Read only the sum over the transpose pair `E_ij + E_ji`, or --
   !> equivalently -- contract all nine seeds against a weight that is symmetric
   !> in `(i, j)`. This routine is linear in the seed, so `f(E_ij) + f(E_ji)` is
   !> `f(E_ij + E_ji)` and the error cancels to machine zero. It is the same
   !> rescue that [[apply_seed]]'s shipped `dB12` shortcut has always relied on,
   !> and it is why the shipped first-order adjoint is correct.
   !>
   !> The suite cannot catch a driver that gets this wrong. Every fixture seeds
   !> a symmetric Hessian by design, so the tests sit entirely inside the region
   !> where the collapse is valid; a per-seed driver will pass all of them.
   !>
   !> TODO: this precondition could be retired rather than documented, by
   !>       symmetrising the seed -- either where [[seed_jet_basis]] emits it or
   !>       on entry to [[apply_seed]]. Not on entry *here*: `dstate_b%dM_u` is
   !>       already built by then. It is known safe: the nuclear path contracts
   !>       these seeds against `f3_rr_rA`, which is symmetric in its two spatial
   !>       indices by equality of mixed partials and measures bit-for-bit
   !>       symmetric on SvdW and CFC, so `0.5*(w_ij + w_ji) == w_ij` and only the
   !>       per-seed summation order moves. It was left undone because the LSF
   !>       interface deliberately documents `w2` as a *general* `3x3` with no
   !>       symmetry assumption, and on the host path `w_lsf2` comes from the
   !>       caller through `api.f90`; symmetrising would quietly narrow that
   !>       contract to buy something this comment already provides
   !>
   !> --------------------------------------------------------------------------
   !>
   !> `res_b` and `dstate_b` are consumed, never recomputed. One
   !> `apply_seed(state, dlsf1_r, ..., res_b, dstate_b)` call is hoisted out of
   !> the direction loop: `b` ranges over the basis seeds while `v` ranges over
   !> `3N` nuclear directions, so rebuilding the `b` chain per direction is the
   !> dominant avoidable cost. It is also the anti-drift rule, since a second
   !> copy of a floating-point chain is free to contract differently under
   !> `-ffp-contract=fast`.
   !>
   !> Contract on the arguments: `dstate_v` and `res_v` must be the *true*
   !> v-tangents of the derived state for the same input displacement `dinp_v`.
   !> [[apply_seed]] is the only producer of such a pair today and it carries no
   !> seed channel for `anchor_wleb0`, and it builds `dwleb_prune_factor` out of
   !> `res%dJ` and `res%dw_f`. A pair it produces is therefore consistent only
   !> when `dinp_v%danchor_wleb0` vanishes and, with pruning active, when
   !> `dinp_v%dcpjac_scal0` and `dinp_v%dw_f0` are `res_v%dJ` and `res_v%dw_f`.
   !> The remaining input tangents are unconstrained. The routine itself makes
   !> none of these identifications; they belong to the driver.
   !>
   !> `d_v(lsf3_rrr)` enters only through `dres%dH`. `lsf4_rrrr` is the driver's
   !> concern and does not appear here.
   !>
   !> `dlsf1_r` and `dlsf2_rr` are deliberately unused, and the compiler warns
   !> about both. [[apply_seed]] is *linear* in the seed, so differentiating a
   !> term `c(p) * b` gives `dc * b + c * db`: a seed component survives here
   !> only if its coefficient is state dependent. Those two enter with the
   !> identity as coefficient, so only their own tangents `ddlsf1_r`/`ddlsf2_rr`
   !> appear, while `dr` and `dlambda` survive through `lsf2_rr` and `lsf3_rrr`.
   !> They are kept so the call site mirrors [[apply_seed]] argument for
   !> argument; nothing is missing.
   !>
   !> @param[in]  state      Per-grid point forward state from [[build_seed_state]]
   !> @param[in]  dstate_v   Tangent of the derived state along the second direction
   !> @param[in]  dinp_v     Tangent of the state inputs along the second direction
   !> @param[in]  res_v      Response of the second direction
   !> @param[in]  dlsf1_r    Seed perturbation of `grad S` at fixed `r`
   !> @param[in]  dlsf2_rr   Seed perturbation of `grad^2 S` at fixed `r`
   !> @param[in]  dr         Induced motion of the projected point
   !> @param[in]  dlambda    Induced change of the Lagrange multiplier
   !> @param[in]  ddlsf1_r   Second-direction tangent of `dlsf1_r`
   !> @param[in]  ddlsf2_rr  Second-direction tangent of `dlsf2_rr`
   !> @param[in]  ddr        Second-direction tangent of `dr`
   !> @param[in]  ddlambda   Second-direction tangent of `dlambda`
   !> @param[in]  res_b      Response of seed `b`, from [[apply_seed]]
   !> @param[in]  dstate_b   State tangent of seed `b`, from [[apply_seed]]
   !> @param[out] dres       Second-order response
   pure subroutine apply_seed_tangent(state, dstate_v, dinp_v, res_v, &
                                      dlsf1_r, dlsf2_rr, dr, dlambda, &
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
      !> The seed `b` being differentiated
      real(wp), intent(in) :: dlsf1_r(3), dlsf2_rr(3, 3), dr(3), dlambda
      !> Tangent of that seed along `v`
      real(wp), intent(in) :: ddlsf1_r(3), ddlsf2_rr(3, 3), ddr(3), ddlambda
      !> Response and state tangent of seed `b`, from [[apply_seed]]
      type(drop_seed_result_type), intent(in) :: res_b
      type(drop_seed_state_tangent_type), intent(in) :: dstate_b
      !> Second-order response
      type(drop_seed_result_tangent_type), intent(out) :: dres

      !> Second-order sensitivity of the tangent-restricted KKT matrix
      real(wp) :: ddA(3, 3)
      !> [[apply_seed]]'s Gram-Schmidt scratch for seed `b`, and its own tangent
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
      !> Curvature scratch, first order in `b` and second order. `dH_b q_a` is
      !> the one partial the total `dstate%dHq_a` does not already carry
      real(wp) :: dHq1_b(3), dHq2_b(3)
      real(wp) :: dN11_b, dN12_b, dN22_b, ddN11, ddN12, ddN22
      real(wp) :: ddS11, ddS12, ddS22, ddT, ddhalf_diff, dd_disc
      !> Cartesian index
      integer :: kaxis

      !*  =============================== Level-set jet ================================ *!

      dres%dg = ddlsf1_r + matmul(dinp_v%dlsf2_rr, dr) + matmul(state%lsf2_rr, ddr)

      dres%dH = ddlsf2_rr
      do kaxis = 1, 3
         dres%dH(:, :) = dres%dH(:, :) + dinp_v%dlsf3_rrr(:, :, kaxis)*dr(kaxis) &
                         + state%lsf3_rrr(:, :, kaxis)*ddr(kaxis)
      end do

      ddA = -ddlambda*state%lsf2_rr - dlambda*dinp_v%dlsf2_rr &
            - dinp_v%dlambda_val*res_b%dH - state%lambda_val*dres%dH

      ! `res%d_gnorm` is the same contraction `n . res%dg` that the `res%dn_surf`
      ! line subtracts, so it is hoisted above its position in [[apply_seed]]
      ! rather than formed twice
      dres%d_gnorm = dot_product(res_v%dn_surf, res_b%dg) &
                     + dot_product(state%n_surf, dres%dg)
      dres%dn_surf = (dres%dg - res_v%dn_surf*res_b%d_gnorm &
                      - state%n_surf*dres%d_gnorm)/state%g_norm &
                     - res_b%dn_surf*res_v%d_gnorm/state%g_norm

      !*  =============================== Tangent frame ================================ *!

      ! `min_axis` is a frozen discrete choice, so the Gram-Schmidt axis is fixed
      v_tmp_b = -res_b%dn_surf(state%min_axis)*state%n_surf &
                - state%n_dot_q1*res_b%dn_surf
      dv_tmp = -dres%dn_surf(state%min_axis)*state%n_surf &
               - res_b%dn_surf(state%min_axis)*res_v%dn_surf &
               - dstate_v%dn_dot_q1*res_b%dn_surf &
               - state%n_dot_q1*dres%dn_surf

      ! Guard consistency. Every branch in this routine tests exactly the
      ! primal's condition on exactly the primal's threshold and puts zero on the
      ! else, because there the primal set the quantity itself to zero. A tangent
      ! gated on a larger threshold than the value it differentiates would be
      ! wrong on the strip between the two, so if any of these gates has to be
      ! raised, [[apply_seed]] and [[apply_seed_tangent]] must be raised together.
      ! Two of [[apply_seed]]'s five gates do not appear here: the eigenvector
      ! chain (`vmin_offdiag .and. sqrt_disc_B > seed_eig_gap_guard`) guards only
      ! `dstate%dvmin_B`, and `use_wleb_prune` guards only
      ! `dstate%dwleb_prune_factor`. Neither lies on the path to
      ! [[drop_seed_result_type]]; both reach this routine already guarded, as
      ! components of `dstate_v`
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

      !*  ================================= KKT matrix ================================= *!

      ! `dstate%dAq1`/`dAq2` carry the FULL product rule `dA q + A dq`, unlike the
      ! partial scratch [[apply_seed]] feeds into its own `dB` lines, where the
      ! `dq` half rides separately. Restated in stored-field form, and using the
      ! symmetry of `A` together with the seed precondition on `dA`, all three
      ! entries share one shape,
      !    dBij = dqi . Aqj + qi . dAqj
      ! which is algebraically identical to the shipped lines and is what is
      ! differentiated below
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

      ! `dBinvXY` is `N/det^2`, so the outer quotient rule reuses the stored
      ! first-order value rather than rebuilding `N`
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

      !*  ========================= Basis-invariant eigenvalue ========================= *!

      ! `u . (ddM u)` in contracted form. The assembly it replaces built `ddP`
      ! out of eight `spread`s and `ddM` out of fourteen matrix-matrix products,
      ! all of it to reach this one scalar, in the innermost `(b, v)` loop; see
      ! [[switched_eigenvalue_curvature]] for the identity and its three
      ! preconditions. Both cross terms `dn_b dn_v^T` and `dn_v dn_b^T` of
      ! `d_v(dP_b)` are still in there, as the `a_v (dn_b . Au)` and
      ! `a_b (dn_v . Au)` pair; dropping either one remains the easiest error
      ! available here
      call switched_eigenvalue_curvature(state%n_surf, state%u_switch, state%A_mat, &
                                         dstate_b%dA_mat, dstate_v%dA_mat, ddA, &
                                         res_b%dn_surf, res_v%dn_surf, dres%dn_surf, &
                                         ddlambda_curv)

      ! The two eigenvector terms of `d_v(u . dM_b u)` coincide when `dM_b` is
      ! symmetric, which is what the seed precondition above buys. It is `dM_b`
      ! and not `M` that has to be symmetric here. `dstate_b%dM_u` is read
      ! rather than rebuilt because it depends on `b` alone; see its declaration
      ddlambda_switch = 2.0_wp*dot_product(dstate_v%du_switch, dstate_b%dM_u) &
                        + ddlambda_curv

      !*  ====================== Lifted tangents and the Jacobian ====================== *!

      ! The sphere tangent frame is rigid at every order, so `dt1 = dt2 = 0`
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

      ! `alpha_coeff` is a fixed parameter with no seed channel, so it
      ! contributes no `dalpha` term of its own here; see
      ! [[drop_seed_input_tangent_type]] for why that omission is deliberate
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

      !*  ============================ Switching and weights =========================== *!

      dres%dw_f = dstate_v%df_foc_f0*state%f_crit_dS*res_b%d_gnorm &
                  + state%f_foc_f0*dstate_v%df_crit_dS*res_b%d_gnorm &
                  + state%f_foc_f0*state%f_crit_dS*dres%d_gnorm &
                  + dstate_v%df_crit0*state%f_foc_dS*dstate_b%dlambda_switch &
                  + state%f_crit0*dstate_v%df_foc_dS*dstate_b%dlambda_switch &
                  + state%f_crit0*state%f_foc_dS*ddlambda_switch

      ! `dw_pre` is a local of [[apply_seed]] and is not stored. Rebuilding it
      ! here is two products off `res_b`, not a re-run of the `b` chain
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

      ! `want_curvature` is a frozen discrete choice, so the early return mirrors
      ! the primal's and `intent(out)` default initialisation zeroes `dk1`/`dk2`
      if (.not. state%want_curvature) return

      !*  ============================ Curvature invariants ============================ *!

      associate (H => state%lsf2_rr, dH_b => res_b%dH, ddH => dres%dH)
         ! `d(q_a . H q_b) = dq_a . H q_b + q_a . d(H q_b)`, and `d(H q_b)` is
         ! `dstate_b%dHq_b` already, so the `q_a^T dH_b q_b` matvec of
         ! [[apply_seed]]'s three-term form is not rebuilt here. `ddN_ab` below
         ! still leans on the symmetry of `H` and `dH_b`, exactly as
         ! [[apply_seed]]'s own `dN12` line and its `dB12` shortcut do
         dHq1_b = matmul(dH_b, state%q1)
         dHq2_b = matmul(dH_b, state%q2)

         dN11_b = dot_product(dstate_b%dq1, state%Hq1) &
                  + dot_product(state%q1, dstate_b%dHq1)
         dN12_b = dot_product(dstate_b%dq1, state%Hq2) &
                  + dot_product(state%q1, dstate_b%dHq2)
         dN22_b = dot_product(dstate_b%dq2, state%Hq2) &
                  + dot_product(state%q2, dstate_b%dHq2)

         ! d_v of [[apply_seed]]'s three-term form, with the two terms that
         ! differ only by which frame vector moves folded into the total `dHq`:
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

         ! d_v of `disc d(disc) = half_diff d(half_diff) + S12 dS12`, solved for
         ! `dd_disc`; the `- d_disc^v d_disc^b` term is the one that comes off
         ! the left-hand side
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
   ! TODO: This is serial; it sh/could be parallelized
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
