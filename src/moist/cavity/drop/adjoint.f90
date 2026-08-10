!> Reverse-mode surface -> level set adjoint contractions for the DROP cavity.
!>
!> These routines map per-tessera surface adjoint weights (Gaussian width,
!> integration weight, area, switch factor, projected position, and normal) onto
!> adjoint weights of the level set function value/gradient/Hessian.
!>
!> They provide the variational solvation-potential (Fock) response; the
!> nuclear-derivative contractions of the same quantities live in gradient.f90.
submodule(moist_cavity_drop) moist_cavity_drop_adjoint
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_linalg, only: eig_2x2_symmetric
   implicit none

contains

   !> Map accumulated surface adjoints into the generic potential container
   !>
   !> @param[inout] self      DROP cavity instance
   !> @param[in]    acc       Accumulated surface-observable adjoints
   !> @param[inout] potential Potential accumulator receiving the LSF channels
   !> @param[out]   error     Error object
   module subroutine get_surface_potential_drop(self, acc, potential, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Potential accumulator receiving the LSF channels
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      allocate (w0(self%ngrid), w1(3, self%ngrid), w2(3, 3, self%ngrid))
      call self%contract_surface_lsf_weights(acc, w0, w1, w2, error)
      if (allocated(error)) return

      if (.not. allocated(potential%w_lsf0)) then
         allocate (potential%w_lsf0(self%ngrid), source=0.0_wp)
         allocate (potential%w_lsf1(3, self%ngrid), source=0.0_wp)
         allocate (potential%w_lsf2(3, 3, self%ngrid), source=0.0_wp)
      else if (size(potential%w_lsf0) /= self%ngrid) then
         call fatal_error(error, "DROP surface potential grid-size mismatch")
         return
      end if
      potential%w_lsf0 = potential%w_lsf0 + w0
      potential%w_lsf1 = potential%w_lsf1 + w1
      potential%w_lsf2 = potential%w_lsf2 + w2

   end subroutine get_surface_potential_drop

   !> Contract per-grid surface adjoint weights to LSF value/gradient/Hessian adjoints
   !>
   !> Implements the DROP reverse-mode chain rule: a perturbation $p$
   !> of the level set function (value $S$, gradient $\nabla S$, Hessian $\nabla^2 S$) moves
   !> the projected surface point and every weight derived from it.
   !>
   !> This routine rewrites the upstream surface adjoint as the LSF-local adjoint,
   !> enforcing for every perturbation $p$ the identity
   !>
   !> $$
   !> \sum_i \Big[ w^{\xi}_i \, \frac{\partial \xi_i}{\partial p}
   !>            + \mathbf{w}^{\mathrm{xyz}}_i \cdot \frac{\partial \mathbf{r}_i}{\partial p} \Big]
   !> = \sum_i \Big[ w^{S}_i \, \frac{\partial S_i}{\partial p}
   !>            + \mathbf{w}^{S_r}_i \cdot \frac{\partial (\nabla S_i)}{\partial p}
   !>            + \sum_{a,b} w^{S_{rr}}_{ab,i} \, \frac{\partial (\nabla^2 S_i)_{ab}}{\partial p} \Big],
   !> $$
   !>
   !> where $w^{\xi}$ = `acc%w_xi`, $\mathbf{w}^{\mathrm{xyz}}$ = `acc%w_xyz`, and
   !> $w^{S}, \mathbf{w}^{S_r}, w^{S_{rr}}$ = `w_lsf0`, `w_lsf1`, `w_lsf2`. The
   !> outward-normal (`acc%w_n`) and principal-curvature (`acc%w_k1`, `acc%w_k2`)
   !> channels enter the same left-hand sum through their own $\partial/\partial p$
   !> sensitivities and are folded into the same `w_lsf` weights. The derived
   !> area (`acc%w_a`) and integration-weight (`acc%w_w`) channels are folded into
   !> the Gaussian-width channel via $a_i = c f_i/\xi_i^2$ and $w_i = c/\xi_i^2$.
   !>
   !> @param[in]  self    DROP cavity instance (must hold a projected grid)
   !> @param[in]  acc     Accumulated surface-observable adjoints
   !> @param[out] w_lsf0  Adjoint weights for LSF values S_i (ngrid)
   !> @param[out] w_lsf1  Adjoint weights for LSF gradients S_r_i (3, ngrid)
   !> @param[out] w_lsf2  Adjoint weights for LSF Hessians S_rr_i (3, 3, ngrid)
   !> @param[out] error   Error object, allocated on failure (KKT sensitivity solve)
   module subroutine contract_surface_lsf_weights(self, acc, w_lsf0, w_lsf1, w_lsf2, error)
      !> DROP cavity instance (must hold a projected grid)
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Adjoint weights for LSF values S_i (ngrid)
      real(wp), intent(out) :: w_lsf0(:)
      !> Adjoint weights for LSF gradients S_r_i (3, ngrid)
      real(wp), intent(out) :: w_lsf1(:, :)
      !> Adjoint weights for LSF Hessians S_rr_i (3, 3, ngrid)
      real(wp), intent(out) :: w_lsf2(:, :, :)
      !> Error object, allocated on failure (KKT sensitivity solve)
      type(error_type), allocatable, intent(out) :: error

      type(lsf_thread_slot) :: lsf_slot
      type(moist_cavity_drop_objective_phi_type) :: phi
      integer :: igrid, ibasis, iaxis, jaxis, kaxis, min_axis_idx
      integer :: igroup_start, igroup_end, group_size, im_grid, m_branch
      real(wp) :: point(3), anchor(3)
      integer :: owner_idx
      real(wp) :: lsf0, lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      real(wp) :: phi0, phi1_r(3), phi2_rr(3, 3)
      real(wp) :: lambda_val
      real(wp) :: kkt_mat_base(4, 4), kkt_mat(4, 4), kkt_rhs(4, 4)
      integer(lapack_ik) :: kkt_ipiv(4), kkt_info
      real(wp) :: dlsf0, dlsf1_r(3), dlsf2_rr(3, 3)
      real(wp) :: dr_dp(3), dlambda_dp, dg_dp(3), dH_dp(3, 3), dA_dp(3, 3)
      real(wp) :: alpha_coeff, g_norm, f_crit0, f_crit_dS, f_foc_f0, f_foc_dS
      real(wp) :: A_mat(3, 3), n_surf(3), q1(3), q2(3), Aq1(3), Aq2(3)
      real(wp) :: B11, B12, B22, tr_B, det_B, disc, sqrt_disc, beta1, beta2
      real(wp) :: lambda_switch, Binv11, Binv12, Binv22, vmin_B(2), vmax_B(2)
      real(wp) :: t1_vec(3), t2_vec(3), tau1(2), tau2(2), w1(2), w2(2)
      real(wp) :: y1(3), y2(3), cross_vec(3), J_val, inv_J, u_switch(3)
      real(wp) :: P_tan(3, 3), dP_tan(3, 3), AP_tan(3, 3), dM_tan(3, 3)
      real(wp) :: n_dot_q1_surf, proj_surf, v_norm_surf
      real(wp) :: dn_surf_dp(3), v_tmp(3), dq1_dp(3), dq2_dp(3)
      real(wp) :: dAq1(3), dAq2(3), dB11, dB12, dB22, ddet_B
      real(wp) :: dBinv11, dBinv12, dBinv22, dlambda_switch
      real(wp) :: dtau1(2), dtau2(2), dw1(2), dw2(2)
      real(wp) :: dy1_dp(3), dy2_dp(3), dcross_dp(3), dJ_dp
      real(wp) :: w_pre_i, f_wleb_s, f_wleb_ds, wleb_prune_factor
      real(wp) :: d_gnorm, dw_f_dp, dw_pre_dp, dwleb_dp, dxi_dp, contribution
      real(wp) :: sigma_phi, adj_wleb, adj_branch, mean_adj_branch, factor_m
      real(wp), allocatable :: branch_phi_adj(:)

      !> Effective primitive surface adjoints used by the DROP contraction
      real(wp), allocatable :: w_xi(:), w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp) :: normal_grad(3), w_xyz_local(3), nwn
      logical :: have_wn
      real(wp), parameter :: weight_tol = 1.0e-30_wp
      real(wp), parameter :: det_B_guard = 1.0e-30_wp

      !> Principal-curvature adjoint intermediates (frame-free invariants,
      !> matching the forward compute_curvature and the nuclear gradient).
      logical :: have_wk
      real(wp) :: g_norm_sq_c
      real(wp) :: Hn_curv(3), Cn_curv(3), adjH_curv(3, 3)
      real(wp) :: trH_curv, nHn_curv, T_curv, nCn_curv, D_curv, KM_curv, disc_curv
      real(wp) :: dadjH_curv(3, 3)
      real(wp) :: dtrH_c, dnHn_c, dT_c, dnCn_c, dD_c, d_disc_c, dk1_dp, dk2_dp
      real(wp), parameter :: curv_disc_guard = 1.0e-10_wp

      if (.not. acc%is_initialized()) then
         call fatal_error(error, "contract_surface_lsf_weights: accumulator is not initialized")
         return
      end if
      if (size(acc%w_xi) /= self%ngrid) then
         call fatal_error(error, "contract_surface_lsf_weights: accumulator grid size mismatch")
         return
      end if
      if (.not. allocated(self%xi0) .or. .not. allocated(self%a) .or. &
          .not. allocated(self%wleb)) then
         call fatal_error(error, "contract_surface_lsf_weights: cavity surface data are incomplete")
         return
      end if
      if (any(abs(self%xi0) <= weight_tol .and. &
              (abs(acc%w_a) > weight_tol .or. abs(acc%w_w) > weight_tol))) then
         call fatal_error(error, "contract_surface_lsf_weights: singular derived-weight conversion")
         return
      end if

      allocate (w_xi, source=acc%w_xi)
      allocate (w_xyz, source=acc%w_xyz)
      allocate (w_n, source=acc%w_n)
      allocate (w_k1, source=acc%w_k1)
      allocate (w_k2, source=acc%w_k2)
      where (abs(acc%w_a) > weight_tol)
         w_xi = w_xi - 2.0_wp*self%a*acc%w_a/self%xi0
      end where
      where (abs(acc%w_w) > weight_tol)
         w_xi = w_xi - 2.0_wp*self%wleb*acc%w_w/self%xi0
      end where

      allocate (lsf_slot%lsf, source=self%lsf_model)
      call lsf_slot%lsf%set_max_deriv(3)
      call phi%set_parameters(self%param)
      call phi%set_input(self%mol, self%radii)
      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (branch_phi_adj(self%ngrid), source=0.0_wp)

      w_lsf0 = 0.0_wp
      w_lsf1 = 0.0_wp
      w_lsf2 = 0.0_wp

      have_wn = any(abs(w_n) > weight_tol)
      have_wk = any(abs(w_k1) > weight_tol) .or. any(abs(w_k2) > weight_tol)

      ! Reverse pass for the branch-weight post-pass:
      ! wleb_m = base_m * p_m, where base_m contains anchor weight, cpjac,
      ! switches and pruning, and p_m is the branch softmax weight. The
      ! local grid loop below handles d(base_m)/dS; this pass converts the
      ! remaining xi-induced adjoint dL/dp_m to dL/dPhi_m
      if (self%ngrid > 0 .and. allocated(self%branch_count) &
          .and. any(self%branch_count(1:self%ngrid) > 1)) then
         igroup_start = 1
         do while (igroup_start <= self%ngrid)
            if (self%branch_count(igroup_start) <= 1) then
               igroup_start = igroup_start + 1
               cycle
            end if

            igroup_end = igroup_start
            do while (igroup_end < self%ngrid)
               if (self%anchor_id(igroup_end + 1) /= self%anchor_id(igroup_start)) exit
               igroup_end = igroup_end + 1
            end do
            group_size = igroup_end - igroup_start + 1

            sigma_phi = self%branch_weight%s
            if (sigma_phi <= weight_tol) then
               igroup_start = igroup_end + 1
               cycle
            end if

            mean_adj_branch = 0.0_wp
            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               adj_branch = 0.0_wp
               if (abs(w_xi(im_grid)) > weight_tol .and. self%wleb(im_grid) > weight_tol &
                   .and. self%wbranch(im_grid) > tiny(1.0_wp)) then
                  adj_wleb = -0.5_wp*w_xi(im_grid)*self%xi0(im_grid)/self%wleb(im_grid)
                  factor_m = self%wleb(im_grid)/self%wbranch(im_grid)
                  adj_branch = adj_wleb*factor_m
               end if
               mean_adj_branch = mean_adj_branch + self%wbranch(im_grid)*adj_branch
               branch_phi_adj(im_grid) = adj_branch
            end do

            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               branch_phi_adj(im_grid) = -self%wbranch(im_grid) &
                                         *(branch_phi_adj(im_grid) - mean_adj_branch)/sigma_phi
            end do

            igroup_start = igroup_end + 1
         end do
      end if

      do igrid = 1, self%ngrid
         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         call lsf_slot%lsf%prepare(point)
         call lsf_slot%lsf%f3_rrr_screened(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         call phi%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         kkt_mat_base = 0.0_wp
         kkt_mat_base(1:3, 1:3) = phi2_rr - lambda_val*lsf2_rr
         kkt_mat_base(1:3, 4) = -lsf1_r
         kkt_mat_base(4, 1:3) = lsf1_r

         alpha_coeff = self%param%phi_alpha
         g_norm = sqrt(dot_product(lsf1_r, lsf1_r))
         if (g_norm <= weight_tol) cycle
         call self%f_crit%eval(g_norm, f_crit0, f_crit_dS)

         A_mat = -lambda_val*lsf2_rr
         A_mat(1, 1) = A_mat(1, 1) + alpha_coeff
         A_mat(2, 2) = A_mat(2, 2) + alpha_coeff
         A_mat(3, 3) = A_mat(3, 3) + alpha_coeff

         n_surf = lsf1_r/g_norm

         ! Fold an optional outward-normal adjoint weight into the field channels:
         ! the direct grad-S contribution normal_grad = P_tan(w_n)/|grad S| enters
         ! w_lsf1 at the fixed projected point, and its point-motion coupling
         ! H @ normal_grad augments the effective position weight below.
         if (have_wn) then
            nwn = dot_product(n_surf, w_n(:, igrid))
            normal_grad = (w_n(:, igrid) - n_surf*nwn)/g_norm
            w_lsf1(:, igrid) = w_lsf1(:, igrid) + normal_grad
            w_xyz_local = w_xyz(:, igrid) + matmul(lsf2_rr, normal_grad)
         else
            w_xyz_local = w_xyz(:, igrid)
         end if

         call setup_tangent_frame(n_surf, q1, q2)
         Aq1 = matmul(A_mat, q1)
         Aq2 = matmul(A_mat, q2)
         B11 = dot_product(q1, Aq1)
         B12 = dot_product(q1, Aq2)
         B22 = dot_product(q2, Aq2)
         tr_B = B11 + B22
         det_B = B11*B22 - B12*B12
         if (abs(det_B) <= det_B_guard) cycle
         disc = max(0.25_wp*tr_B*tr_B - det_B, 0.0_wp)
         sqrt_disc = sqrt(disc)
         beta1 = 0.5_wp*tr_B + sqrt_disc
         beta2 = 0.5_wp*tr_B - sqrt_disc
         call eig_2x2_symmetric(B11, B12, B22, lambda_switch, beta1, vmin_B, vmax_B)
         u_switch = vmin_B(1)*q1 + vmin_B(2)*q2
         lambda_switch = beta2
         call self%f_foc%eval(lambda_switch, f_foc_f0, f_foc_dS)

         Binv11 = B22/det_B
         Binv12 = -B12/det_B
         Binv22 = B11/det_B
         call setup_tangent_frame(anchor - self%mol%xyz(:, owner_idx), t1_vec, t2_vec)
         tau1(1) = dot_product(q1, t1_vec)
         tau1(2) = dot_product(q2, t1_vec)
         tau2(1) = dot_product(q1, t2_vec)
         tau2(2) = dot_product(q2, t2_vec)
         w1(1) = Binv11*tau1(1) + Binv12*tau1(2)
         w1(2) = Binv12*tau1(1) + Binv22*tau1(2)
         w2(1) = Binv11*tau2(1) + Binv12*tau2(2)
         w2(2) = Binv12*tau2(1) + Binv22*tau2(2)
         y1 = alpha_coeff*(w1(1)*q1 + w1(2)*q2)
         y2 = alpha_coeff*(w2(1)*q1 + w2(2)*q2)
         cross_vec(1) = y1(2)*y2(3) - y1(3)*y2(2)
         cross_vec(2) = y1(3)*y2(1) - y1(1)*y2(3)
         cross_vec(3) = y1(1)*y2(2) - y1(2)*y2(1)
         J_val = sqrt(dot_product(cross_vec, cross_vec))
         if (J_val <= weight_tol) cycle
         inv_J = 1.0_wp/J_val
         min_axis_idx = minloc(abs(n_surf), dim=1)
         n_dot_q1_surf = n_surf(min_axis_idx)
         proj_surf = 1.0_wp - n_dot_q1_surf**2
         v_norm_surf = sqrt(max(proj_surf, 1.0e-30_wp))
         P_tan(:, :) = -spread(n_surf, dim=2, ncopies=3)*spread(n_surf, dim=1, ncopies=3)
         P_tan(1, 1) = P_tan(1, 1) + 1.0_wp
         P_tan(2, 2) = P_tan(2, 2) + 1.0_wp
         P_tan(3, 3) = P_tan(3, 3) + 1.0_wp
         AP_tan = matmul(A_mat, P_tan)

         if (self%param%wleb_prune_level > 0) then
            w_pre_i = self%anchor_wleb0(igrid)*self%cpjac_scal0(igrid)*self%w_f0(igrid)
            call self%f_wleb%eval(abs(w_pre_i), f_wleb_s, f_wleb_ds)
            wleb_prune_factor = f_wleb_s + abs(w_pre_i)*f_wleb_ds
         else
            wleb_prune_factor = 1.0_wp
         end if

         ! Grid-level principal-curvature invariants (frame independent), used
         ! when a curvature adjoint weight is supplied. Mirrors gradient.f90:
         !   T = k1 + k2 = (tr H - n^T H n)/|g|
         !   D = k1 * k2 = (n^T adj(H) n)/|g|^2
         !   k1,k2 = T/2 +/- sqrt((T/2)^2 - D)
         if (have_wk) then
            g_norm_sq_c = g_norm*g_norm
            Hn_curv = matmul(lsf2_rr, n_surf)
            trH_curv = lsf2_rr(1, 1) + lsf2_rr(2, 2) + lsf2_rr(3, 3)
            nHn_curv = dot_product(n_surf, Hn_curv)
            T_curv = (trH_curv - nHn_curv)/g_norm

            adjH_curv(1, 1) = lsf2_rr(2, 2)*lsf2_rr(3, 3) - lsf2_rr(2, 3)*lsf2_rr(2, 3)
            adjH_curv(2, 2) = lsf2_rr(1, 1)*lsf2_rr(3, 3) - lsf2_rr(1, 3)*lsf2_rr(1, 3)
            adjH_curv(3, 3) = lsf2_rr(1, 1)*lsf2_rr(2, 2) - lsf2_rr(1, 2)*lsf2_rr(1, 2)
            adjH_curv(1, 2) = lsf2_rr(1, 3)*lsf2_rr(2, 3) - lsf2_rr(1, 2)*lsf2_rr(3, 3)
            adjH_curv(1, 3) = lsf2_rr(1, 2)*lsf2_rr(2, 3) - lsf2_rr(2, 2)*lsf2_rr(1, 3)
            adjH_curv(2, 3) = lsf2_rr(1, 2)*lsf2_rr(1, 3) - lsf2_rr(1, 1)*lsf2_rr(2, 3)
            adjH_curv(2, 1) = adjH_curv(1, 2)
            adjH_curv(3, 1) = adjH_curv(1, 3)
            adjH_curv(3, 2) = adjH_curv(2, 3)

            Cn_curv = matmul(adjH_curv, n_surf)
            nCn_curv = dot_product(n_surf, Cn_curv)
            D_curv = nCn_curv/g_norm_sq_c
            KM_curv = 0.5_wp*T_curv
            disc_curv = sqrt(max(KM_curv*KM_curv - D_curv, 0.0_wp))
         end if

         ! KKT sensitivities for all 13 basis perturbations from one
         ! factorization: only the value (ibasis 1) and gradient (ibasis 2-4)
         ! perturbations enter the right-hand side; the nine Hessian
         ! perturbations have rhs = 0 and hence dr/dp = 0, dlambda/dp = 0.
         kkt_rhs = 0.0_wp
         kkt_rhs(4, 1) = -1.0_wp
         kkt_rhs(1, 2) = lambda_val
         kkt_rhs(2, 3) = lambda_val
         kkt_rhs(3, 4) = lambda_val
         kkt_mat = kkt_mat_base
         call lapack_gesv(4_lapack_ik, 4_lapack_ik, kkt_mat, 4_lapack_ik, &
                          kkt_ipiv, kkt_rhs, 4_lapack_ik, kkt_info)
         if (kkt_info /= 0_lapack_ik) then
            call fatal_error(error, &
                             "contract_surface_lsf_weights: KKT sensitivity solve failed")
            return
         end if

         do ibasis = 1, 13
            dlsf0 = 0.0_wp
            dlsf1_r = 0.0_wp
            dlsf2_rr = 0.0_wp
            if (ibasis == 1) then
               dlsf0 = 1.0_wp
            else if (ibasis <= 4) then
               dlsf1_r(ibasis - 1) = 1.0_wp
            else
               iaxis = (ibasis - 5)/3 + 1
               jaxis = mod(ibasis - 5, 3) + 1
               dlsf2_rr(iaxis, jaxis) = 1.0_wp
            end if

            if (ibasis <= 4) then
               dr_dp = kkt_rhs(1:3, ibasis)
               dlambda_dp = kkt_rhs(4, ibasis)
            else
               dr_dp = 0.0_wp
               dlambda_dp = 0.0_wp
            end if
            dg_dp = dlsf1_r + matmul(lsf2_rr, dr_dp)
            dH_dp = dlsf2_rr
            do kaxis = 1, 3
               dH_dp(:, :) = dH_dp(:, :) + lsf3_rrr(:, :, kaxis)*dr_dp(kaxis)
            end do
            dA_dp = -dlambda_dp*lsf2_rr - lambda_val*dH_dp
            dn_surf_dp = (dg_dp - n_surf*dot_product(n_surf, dg_dp))/g_norm
            v_tmp = -dn_surf_dp(min_axis_idx)*n_surf - n_dot_q1_surf*dn_surf_dp
            if (proj_surf > 1.0e-30_wp) then
               dq1_dp = (v_tmp - q1*dot_product(q1, v_tmp))/v_norm_surf
            else
               dq1_dp = 0.0_wp
            end if
            dq2_dp(1) = dn_surf_dp(2)*q1(3) - dn_surf_dp(3)*q1(2) &
                        + n_surf(2)*dq1_dp(3) - n_surf(3)*dq1_dp(2)
            dq2_dp(2) = dn_surf_dp(3)*q1(1) - dn_surf_dp(1)*q1(3) &
                        + n_surf(3)*dq1_dp(1) - n_surf(1)*dq1_dp(3)
            dq2_dp(3) = dn_surf_dp(1)*q1(2) - dn_surf_dp(2)*q1(1) &
                        + n_surf(1)*dq1_dp(2) - n_surf(2)*dq1_dp(1)
            dAq1 = matmul(dA_dp, q1)
            dAq2 = matmul(dA_dp, q2)
            dB11 = 2.0_wp*dot_product(dq1_dp, Aq1) + dot_product(q1, dAq1)
            dB12 = dot_product(dq1_dp, Aq2) + dot_product(dq2_dp, Aq1) &
                   + dot_product(dAq1, q2)
            dB22 = 2.0_wp*dot_product(dq2_dp, Aq2) + dot_product(q2, dAq2)
            ddet_B = dB11*B22 + B11*dB22 - 2.0_wp*B12*dB12
            dBinv11 = (dB22*det_B - B22*ddet_B)/(det_B*det_B)
            dBinv12 = (-dB12*det_B + B12*ddet_B)/(det_B*det_B)
            dBinv22 = (dB11*det_B - B11*ddet_B)/(det_B*det_B)
            dP_tan(:, :) = -(spread(dn_surf_dp, dim=2, ncopies=3)*spread(n_surf, dim=1, ncopies=3) &
                             + spread(n_surf, dim=2, ncopies=3)*spread(dn_surf_dp, dim=1, ncopies=3))
            dM_tan = matmul(dP_tan, AP_tan) &
                     + matmul(P_tan, matmul(dA_dp, P_tan)) &
                     + matmul(P_tan, matmul(A_mat, dP_tan))
            dlambda_switch = dot_product(u_switch, matmul(dM_tan, u_switch))
            dtau1(1) = dot_product(dq1_dp, t1_vec)
            dtau1(2) = dot_product(dq2_dp, t1_vec)
            dtau2(1) = dot_product(dq1_dp, t2_vec)
            dtau2(2) = dot_product(dq2_dp, t2_vec)
            dw1(1) = dBinv11*tau1(1) + Binv11*dtau1(1) &
                     + dBinv12*tau1(2) + Binv12*dtau1(2)
            dw1(2) = dBinv12*tau1(1) + Binv12*dtau1(1) &
                     + dBinv22*tau1(2) + Binv22*dtau1(2)
            dw2(1) = dBinv11*tau2(1) + Binv11*dtau2(1) &
                     + dBinv12*tau2(2) + Binv12*dtau2(2)
            dw2(2) = dBinv12*tau2(1) + Binv12*dtau2(1) &
                     + dBinv22*tau2(2) + Binv22*dtau2(2)
            dy1_dp = alpha_coeff*(dw1(1)*q1 + w1(1)*dq1_dp &
                                  + dw1(2)*q2 + w1(2)*dq2_dp)
            dy2_dp = alpha_coeff*(dw2(1)*q1 + w2(1)*dq1_dp &
                                  + dw2(2)*q2 + w2(2)*dq2_dp)
            dcross_dp(1) = dy1_dp(2)*y2(3) - dy1_dp(3)*y2(2) &
                           + y1(2)*dy2_dp(3) - y1(3)*dy2_dp(2)
            dcross_dp(2) = dy1_dp(3)*y2(1) - dy1_dp(1)*y2(3) &
                           + y1(3)*dy2_dp(1) - y1(1)*dy2_dp(3)
            dcross_dp(3) = dy1_dp(1)*y2(2) - dy1_dp(2)*y2(1) &
                           + y1(1)*dy2_dp(2) - y1(2)*dy2_dp(1)
            dJ_dp = dot_product(cross_vec, dcross_dp)*inv_J
            d_gnorm = dot_product(n_surf, dg_dp)
            dw_f_dp = f_foc_f0*f_crit_dS*d_gnorm + f_crit0*f_foc_dS*dlambda_switch
            dw_pre_dp = self%anchor_wleb0(igrid)*self%w_f0(igrid)*dJ_dp &
                        + self%anchor_wleb0(igrid)*self%cpjac_scal0(igrid)*dw_f_dp
            dwleb_dp = self%wbranch(igrid)*wleb_prune_factor*dw_pre_dp
            if (self%wleb(igrid) > weight_tol) then
               dxi_dp = -0.5_wp*self%xi0(igrid)*dwleb_dp/self%wleb(igrid)
            else
               dxi_dp = 0.0_wp
            end if

            ! self%f is an anchor-only iSwig overlap for electronic LSF
            ! perturbations, so df/dp = 0 and w_f does not contribute here
            contribution = dot_product(w_xyz_local, dr_dp) + w_xi(igrid)*dxi_dp
            if (abs(branch_phi_adj(igrid)) > weight_tol) then
               contribution = contribution + branch_phi_adj(igrid)*dot_product(phi1_r, dr_dp)
            end if

            ! Principal-curvature channel: fold w_k1/w_k2 through the same
            ! frame-free invariant kernel as gradient.f90, reusing this basis's
            ! total grad-S and Hessian sensitivities dg_dp, dH_dp.
            if (have_wk) then
               dtrH_c = dH_dp(1, 1) + dH_dp(2, 2) + dH_dp(3, 3)
               dnHn_c = 2.0_wp*dot_product(dn_surf_dp, Hn_curv) &
                        + dot_product(n_surf, matmul(dH_dp, n_surf))
               dT_c = (dtrH_c - dnHn_c)/g_norm - T_curv*d_gnorm/g_norm

               dadjH_curv(1, 1) = dH_dp(2, 2)*lsf2_rr(3, 3) + lsf2_rr(2, 2)*dH_dp(3, 3) &
                                  - 2.0_wp*lsf2_rr(2, 3)*dH_dp(2, 3)
               dadjH_curv(2, 2) = dH_dp(1, 1)*lsf2_rr(3, 3) + lsf2_rr(1, 1)*dH_dp(3, 3) &
                                  - 2.0_wp*lsf2_rr(1, 3)*dH_dp(1, 3)
               dadjH_curv(3, 3) = dH_dp(1, 1)*lsf2_rr(2, 2) + lsf2_rr(1, 1)*dH_dp(2, 2) &
                                  - 2.0_wp*lsf2_rr(1, 2)*dH_dp(1, 2)
               dadjH_curv(1, 2) = dH_dp(1, 3)*lsf2_rr(2, 3) + lsf2_rr(1, 3)*dH_dp(2, 3) &
                                  - dH_dp(1, 2)*lsf2_rr(3, 3) - lsf2_rr(1, 2)*dH_dp(3, 3)
               dadjH_curv(1, 3) = dH_dp(1, 2)*lsf2_rr(2, 3) + lsf2_rr(1, 2)*dH_dp(2, 3) &
                                  - dH_dp(2, 2)*lsf2_rr(1, 3) - lsf2_rr(2, 2)*dH_dp(1, 3)
               dadjH_curv(2, 3) = dH_dp(1, 2)*lsf2_rr(1, 3) + lsf2_rr(1, 2)*dH_dp(1, 3) &
                                  - dH_dp(1, 1)*lsf2_rr(2, 3) - lsf2_rr(1, 1)*dH_dp(2, 3)
               dadjH_curv(2, 1) = dadjH_curv(1, 2)
               dadjH_curv(3, 1) = dadjH_curv(1, 3)
               dadjH_curv(3, 2) = dadjH_curv(2, 3)

               dnCn_c = 2.0_wp*dot_product(dn_surf_dp, Cn_curv) &
                        + dot_product(n_surf, matmul(dadjH_curv, n_surf))
               dD_c = dnCn_c/g_norm_sq_c - 2.0_wp*D_curv*d_gnorm/g_norm

               if (disc_curv > curv_disc_guard) then
                  d_disc_c = (KM_curv*dT_c - dD_c)/(2.0_wp*disc_curv)
               else
                  d_disc_c = 0.0_wp
               end if
               dk1_dp = 0.5_wp*dT_c + d_disc_c
               dk2_dp = 0.5_wp*dT_c - d_disc_c

               contribution = contribution + w_k1(igrid)*dk1_dp + w_k2(igrid)*dk2_dp
            end if

            if (ibasis == 1) then
               w_lsf0(igrid) = w_lsf0(igrid) + contribution
            else if (ibasis <= 4) then
               w_lsf1(ibasis - 1, igrid) = w_lsf1(ibasis - 1, igrid) + contribution
            else
               iaxis = (ibasis - 5)/3 + 1
               jaxis = mod(ibasis - 5, 3) + 1
               w_lsf2(iaxis, jaxis, igrid) = w_lsf2(iaxis, jaxis, igrid) + contribution
            end if
         end do
      end do

   end subroutine contract_surface_lsf_weights

end submodule moist_cavity_drop_adjoint
