submodule(moist_cavity_drop) moist_cavity_drop_gradient
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_linalg, only: eig_2x2_symmetric
   implicit none

contains

   !> Contract per-grid surface adjoint weights to LSF value/gradient/Hessian adjoints
   !>
   !> Implements the DROP reverse-mode chain rule: a perturbation $p$
   !> of the level-set function (value $S$, gradient $\nabla S$, Hessian $\nabla^2 S$) moves
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
   !> where $w^{\xi}$ = `w_xi`, $\mathbf{w}^{\mathrm{xyz}}$ = `w_xyz`, and
   !> $w^{S}, \mathbf{w}^{S_r}, w^{S_{rr}}$ = `w_lsf0`, `w_lsf1`, `w_lsf2`.
   !>
   !> @param[in]  self    DROP cavity instance (must hold a projected grid)
   !> @param[in]  w_xi    Surface adjoint weights for Gaussian widths xi_i (ngrid)
   !> @param[in]  w_f     Surface adjoint weights for anchor switch factors f_i (ngrid); unused
   !> @param[in]  w_xyz   Surface adjoint weights for projected coordinates r_i (3, ngrid)
   !> @param[out] w_lsf0  Adjoint weights for LSF values S_i (ngrid)
   !> @param[out] w_lsf1  Adjoint weights for LSF gradients S_r_i (3, ngrid)
   !> @param[out] w_lsf2  Adjoint weights for LSF Hessians S_rr_i (3, 3, ngrid)
   !> @param[out] error   Error object, allocated on failure (KKT sensitivity solve)
   module subroutine contract_surface_lsf_weights(self, w_xi, w_f, w_xyz, w_lsf0, w_lsf1, w_lsf2, error)
      !> DROP cavity instance (must hold a projected grid)
      class(cavity_type_drop), intent(in) :: self
      !> Surface adjoint weights for Gaussian widths xi_i (ngrid)
      real(wp), intent(in) :: w_xi(:)
      !> Surface adjoint weights for anchor switch factors f_i (ngrid); unused (df/dp = 0)
      real(wp), intent(in) :: w_f(:)
      !> Surface adjoint weights for projected coordinates r_i (3, ngrid)
      real(wp), intent(in) :: w_xyz(:, :)
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
      real(wp) :: kkt_mat_base(4, 4), kkt_mat(4, 4), kkt_rhs(4, 1)
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
      real(wp), parameter :: weight_tol = 1.0e-30_wp
      real(wp), parameter :: det_B_guard = 1.0e-30_wp

      allocate (lsf_slot%lsf, source=self%lsf_model)
      call lsf_slot%lsf%set_max_deriv(3)
      call phi%set_parameters(self%param)
      call phi%set_input(self%mol, self%radii)
      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (branch_phi_adj(self%ngrid), source=0.0_wp)

      w_lsf0 = 0.0_wp
      w_lsf1 = 0.0_wp
      w_lsf2 = 0.0_wp

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

            kkt_rhs(1:3, 1) = lambda_val*dlsf1_r
            kkt_rhs(4, 1) = -dlsf0
            kkt_mat = kkt_mat_base
            call lapack_gesv(4_lapack_ik, 1_lapack_ik, kkt_mat, 4_lapack_ik, &
                             kkt_ipiv, kkt_rhs, 4_lapack_ik, kkt_info)
            if (kkt_info /= 0_lapack_ik) then
               call fatal_error(error, &
                                "contract_surface_lsf_weights: KKT sensitivity solve failed")
               return
            end if

            dr_dp = kkt_rhs(1:3, 1)
            dlambda_dp = kkt_rhs(4, 1)
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
            contribution = dot_product(w_xyz(:, igrid), dr_dp) + w_xi(igrid)*dxi_dp
            if (abs(branch_phi_adj(igrid)) > weight_tol) then
               contribution = contribution + branch_phi_adj(igrid) &
                              *dot_product(phi1_r, dr_dp)
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

      deallocate (lsf3_rrr, branch_phi_adj)

   end subroutine contract_surface_lsf_weights

   !> Compute first nuclear derivatives (gradients w.r.t. atomic coordinates r_A)
   !> of all DROP per-grid and per-sphere quantities.
   !>
   !> For every projected grid point a bordered-KKT sensitivity system is formed
   !> from the Lagrangian Hessian $\mathbf{H}_L = \phi_{\mathbf{rr}} - \lambda
   !> S_{\mathbf{rr}}$ and the LSF gradient $\mathbf{g} = S_{\mathbf r}$,
   !>
   !> $$
   !> \begin{bmatrix} \mathbf{H}_L & -\mathbf{g} \\ \mathbf{g}^{\top} & 0 \end{bmatrix}
   !> \begin{bmatrix} \partial \mathbf{r} / \partial r_A \\ \partial \lambda / \partial r_A \end{bmatrix}
   !> = \begin{bmatrix} -\big(\phi_{\mathbf r r_A} - \lambda\, S_{\mathbf r r_A}\big) \\ -S_{r_A} \end{bmatrix},
   !> $$
   !>
   !> where $\phi_{\mathbf r r_A}$ / $S_{\mathbf r r_A}$ are mixed spatial-nuclear
   !> second derivatives and $S_{r_A}$ is the explicit nuclear derivative of $S$.
   !>
   !> Solving it (one factorization, all active atoms and axes batched) yields the projected-point
   !> and multiplier derivatives `xyz1_rA`, `lambda1_rA`. These seed the chain rule for the surface normal,
   !> the closest-point Jacobian scaling `cpjac_scal1_rA`, the critical-gradient / focus switch `w_f1_rA`,
   !> the Lebedev weight `wleb1_rA`, the Gaussian width `xi1_rA`, the anchor switch factor `f1_rA`, and the
   !> per-point area / volume elements `a_i1_rA`, `v_i1_rA`, which accumulate into the per-sphere and total
   !> gradients `asph1_rA`, `vsph1_rA`, `A_tot1_rA`, `V_tot1_rA`. Optional point/normal derivatives
   !> (`rho1_rA`, `r_iI1_rA`, `normal1_rA`) are produced only when the corresponding request flag is set.
   !>
   !> @param[inout] self  DROP cavity instance; reads the projected grid and
   !>                     forward quantities, then allocates and fills the `*_rA`
   !>                     first-derivative arrays
   !> @param[out]   error Error object, allocated on failure (KKT sensitivity
   !>                     solve or singular tangent Jacobian)
   module subroutine compute_gradient_drop(self, error)

      !> DROP cavity instance; supplies the projected grid and receives the first-derivative (`*_rA`) arrays
      class(cavity_type_drop), intent(inout) :: self

      !> Error type
      type(error_type), allocatable, intent(out) :: error

      !> LSF thread slots
      type(lsf_thread_slot), allocatable :: lsf_threads(:)

      !> Phi thread slots
      type(moist_cavity_drop_objective_phi_type), allocatable :: phi_threads(:)

      !> Loop indices
      integer :: igrid, iatom, iaxis, jaxis, i, n_active
      integer, allocatable :: active_idx(:)

      !> OpenMP thread management
      integer :: nthreads, thread_slot
      integer :: timer_ref_thread
      logical :: abort_requested

      !> Timing flag
      logical :: do_timing

      !> Grid point data
      real(wp) :: point(3), anchor(3), rho_vec(3), rho_norm
      integer :: owner_idx

      !> LSF derivatives
      real(wp) :: lsf0, lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)

      !> Phi derivatives
      real(wp) :: phi0, phi1_r(3), phi2_rr(3, 3)
      real(wp), allocatable :: phi2_r_rA(:, :, :)

      !> KKT system
      real(wp) :: lambda_val
      real(wp) :: G_lagrangian(3), H_lagrangian(3, 3)
      real(wp) :: kkt_mat_base(4, 4), kkt_mat(4, 4), kkt_rhs(4, 1)
      real(wp) :: rhs_vec(4)
      integer(lapack_ik) :: kkt_ipiv(4), kkt_info
      real(wp), allocatable :: kkt_rhs_batch(:, :)

      !> swi: Rho derivatives
      real(wp) :: rho_unit(3), delta_matrix(3, 3)

      !> swi: POU derivatives
      real(wp) :: iswig_f0

      !> Point derivatives
      real(wp) :: r_iI0, r_iI_vec(3), r_iI_unit(3), r_iI_norm

      !> Volume derivatives
      real(wp) :: r_hat_dot_r, grad_r_hat_dot_r(3)

      !> Thread-local buffer for normal derivatives (used by volume gradient)
      real(wp), allocatable :: dn_dR_buf(:, :, :)

      !> Jacobian scaling derivatives
      real(wp) :: alpha_coeff, g_vec(3), g_norm_sq, g_norm
      real(wp) :: A_mat(3, 3)
      real(wp) :: t1_vec(3), t2_vec(3), y1(3), y2(3)
      real(wp) :: cross_vec(3), J_val, inv_J
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :)
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      real(wp) :: dA_dR(3, 3)
      real(wp) :: dg_dR(3)

      !> w_f switching function derivative intermediates
      real(wp) :: f_crit0, f_crit_dS, f_foc_f0, f_foc_dS, d_gnorm

      !> Lebedev weight switching intermediates
      real(wp) :: w_pre_i, f_wleb_s, f_wleb_ds, wleb_prune_factor, dw_pre_dR

      !> 2x2 tangent-restricted inverse and switch variables
      real(wp) :: n_surf(3), q1(3), q2(3)
      real(wp) :: Aq1(3), Aq2(3)
      real(wp) :: B11, B12, B22, tr_B, det_B, disc, sqrt_disc
      real(wp) :: beta1, beta2, lambda_switch
      real(wp), parameter :: det_B_guard = 1.0e-30_wp
      real(wp) :: Binv11, Binv12, Binv22
      real(wp) :: tau1(2), tau2(2), w1(2), w2(2)
      real(wp) :: vmin_B(2), vmax_B(2), u_switch(3)
      real(wp) :: P_tan(3, 3), dP_tan(3, 3), M_tan(3, 3), dM_tan(3, 3)
      real(wp) :: AP_tan(3, 3)

      !> Jacobian derivative intermediates
      real(wp) :: dn_dR(3)
      real(wp) :: dy1_dR(3), dy2_dR(3)
      real(wp) :: dcross_dR(3), dJ_dR
      real(wp) :: dlambda_val, dr_i_dR(3)
      real(wp) :: v_tmp(3)
      real(wp) :: dn_surf_dR(3), dq1_dR(3), dq2_dR(3)
      real(wp) :: dAq1(3), dAq2(3)
      real(wp) :: dB11, dB12, dB22, ddet_B
      real(wp) :: dBinv11, dBinv12, dBinv22
      real(wp) :: dtau1(2), dtau2(2), dw1(2), dw2(2), dlambda_switch
      integer :: min_axis_surf
      real(wp) :: proj_surf, v_norm_surf, n_dot_q1_surf

      ! Branch-weight post-pass state (serial, after main loop).
      ! Softmax weights_grad takes dphi in (nparam, nbranch) layout where
      ! nparam = 3 * nsph; we flatten (iatom, iaxis) -> (iatom - 1) * 3 + iaxis.
      integer :: igroup_start, igroup_end, group_size, m_branch, im_grid
      integer :: owner_m, k_param
      real(wp) :: pt_m(3), anch_m(3), phi1_r_m(3), dphi_m, factor_m
      real(wp) :: area_fac_m, rn_m, dwleb_branch, da_branch, dv_branch
      real(wp), allocatable :: branch_phi(:), branch_dphi(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)

      ! Thread-local xi buffer (avoids storing debug arrays)
      real(wp), allocatable :: anchor_xi_local(:, :)
      ! Scalar temps for atomic reduction (avoids gfortran aliasing issue)
      real(wp) :: ai_val, vi_val

      ! Thread-local accumulators (replace atomics for A_tot and V_tot)
      real(wp), allocatable :: A_tot_local(:, :), V_tot_local(:, :)

      ! Initialize thread-local primitives
      nthreads = max(1, omp_get_max_threads())
      timer_ref_thread = 1
      allocate (lsf_threads(nthreads))
      allocate (phi_threads(nthreads))
      do thread_slot = 1, nthreads
         allocate (lsf_threads(thread_slot)%lsf, source=self%lsf_model)
         ! Gradient uses third spatial derivatives (f3_rrr_screened, f3_rr_rA_screened);
         ! upgrade SSD storage so f3_rrr_arr is allocated before any %prepare call.
         call lsf_threads(thread_slot)%lsf%set_max_deriv(3)
         call phi_threads(thread_slot)%set_parameters(self%param)
         call phi_threads(thread_slot)%set_input(self%mol, self%radii)
      end do

      ! Timing proxy: only one thread contributes to section timers.
      do_timing = .true.

      ! Allocate gradient arrays
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      allocate (self%xyz1_rA(3, 3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%lambda1_rA)) deallocate (self%lambda1_rA)
      allocate (self%lambda1_rA(3, self%nsph, self%ngrid), source=0.0_wp)

      ! Optional derivative arrays (gated by request flags)
      if (self%request%rho) then
         if (allocated(self%rho1_rA)) deallocate (self%rho1_rA)
         allocate (self%rho1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      end if
      if (self%request%r_iI) then
         if (allocated(self%r_iI1_rA)) deallocate (self%r_iI1_rA)
         allocate (self%r_iI1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      end if
      if (self%request%normal) then
         if (allocated(self%normal1_rA)) deallocate (self%normal1_rA)
         allocate (self%normal1_rA(3, self%nsph, 3, self%ngrid), source=0.0_wp)
      end if
      if (allocated(self%a_i1_rA)) deallocate (self%a_i1_rA)
      allocate (self%a_i1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%v_i1_rA)) deallocate (self%v_i1_rA)
      allocate (self%v_i1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%A_tot1_rA)) deallocate (self%A_tot1_rA)
      allocate (self%A_tot1_rA(3, self%nsph), source=0.0_wp)
      if (allocated(self%asph1_rA)) deallocate (self%asph1_rA)
      allocate (self%asph1_rA(3, self%nsph, self%nsph), source=0.0_wp)
      if (allocated(self%V_tot1_rA)) deallocate (self%V_tot1_rA)
      allocate (self%V_tot1_rA(3, self%nsph), source=0.0_wp)
      if (allocated(self%vsph1_rA)) deallocate (self%vsph1_rA)
      allocate (self%vsph1_rA(3, self%nsph, self%nsph), source=0.0_wp)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      allocate (self%f1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%cpjac_scal1_rA)) deallocate (self%cpjac_scal1_rA)
      allocate (self%cpjac_scal1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%w_f1_rA)) deallocate (self%w_f1_rA)
      allocate (self%w_f1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%wleb1_rA)) deallocate (self%wleb1_rA)
      allocate (self%wleb1_rA(3, self%nsph, self%ngrid), source=0.0_wp)

      ! TODO: These are currently optional (just for testing)
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      allocate (self%xi1_rA(3, self%nsph, self%ngrid), source=0.0_wp)

      abort_requested = .false.

      ! Loop over all grid points (parallelized).
      !$omp parallel num_threads(nthreads) default(shared) private(thread_slot, do_timing, igrid, &
      !$omp& iatom, iaxis, jaxis, point, anchor, rho_vec, rho_norm, owner_idx, lsf0, lsf1_r, lsf2_rr, &
      !$omp& lsf1_rA, lsf2_r_rA, phi0, phi1_r, phi2_rr, phi2_r_rA, lambda_val, &
      !$omp& G_lagrangian, H_lagrangian, kkt_mat_base, kkt_mat, kkt_rhs, rhs_vec, &
      !$omp& kkt_ipiv, kkt_info, rho_unit, &
      !$omp& delta_matrix, r_iI_vec, r_iI_norm, r_hat_dot_r, &
      !$omp& grad_r_hat_dot_r, alpha_coeff, g_vec, g_norm_sq, g_norm, A_mat, &
      !$omp& t1_vec, t2_vec, y1, y2, cross_vec, &
      !$omp& J_val, inv_J, lsf3_rr_rA, lsf3_rrr, dA_dR, dg_dR, &
      !$omp& n_surf, q1, q2, Aq1, Aq2, B11, B12, B22, tr_B, det_B, disc, sqrt_disc, &
      !$omp& beta1, beta2, lambda_switch, Binv11, Binv12, Binv22, tau1, tau2, w1, w2, &
      !$omp& vmin_B, vmax_B, u_switch, P_tan, dP_tan, M_tan, dM_tan, &
      !$omp& dn_dR, dy1_dR, dy2_dR, dcross_dR, dJ_dR, dlambda_val, dr_i_dR, &
      !$omp& v_tmp, dn_surf_dR, dq1_dR, dq2_dR, dAq1, dAq2, &
      !$omp& dB11, dB12, dB22, ddet_B, dBinv11, dBinv12, dBinv22, &
      !$omp& dtau1, dtau2, dw1, dw2, &
      !$omp& dlambda_switch, &
      !$omp& min_axis_surf, proj_surf, v_norm_surf, n_dot_q1_surf, &
      !$omp& i, n_active, active_idx, &
      !$omp& f_crit0, f_crit_dS, f_foc_f0, f_foc_dS, d_gnorm, dn_dR_buf, anchor_xi_local, &
      !$omp& w_pre_i, f_wleb_s, f_wleb_ds, wleb_prune_factor, dw_pre_dR, &
      !$omp& ai_val, vi_val, kkt_rhs_batch, AP_tan, &
      !$omp& A_tot_local, V_tot_local)
      thread_slot = omp_get_thread_num() + 1
      do_timing = thread_slot == timer_ref_thread

      allocate (lsf3_rr_rA(3, 3, 3, self%nsph))
      allocate (lsf3_rrr(3, 3, 3))
      allocate (active_idx(self%nsph))
      allocate (dn_dR_buf(3, self%nsph, 3))
      allocate (lsf1_rA(3, self%nsph))
      allocate (lsf2_r_rA(3, 3, self%nsph))
      allocate (phi2_r_rA(3, 3, self%nsph))
      allocate (anchor_xi_local(3, self%nsph))
      allocate (kkt_rhs_batch(4, 3*self%nsph))
      allocate (A_tot_local(3, self%nsph), source=0.0_wp)
      allocate (V_tot_local(3, self%nsph), source=0.0_wp)

      !> The anchor_xi depends only on the nuclear geometry of the anchor system
      anchor_xi_local = 0.0_wp

      !$omp do schedule(dynamic)
      do igrid = 1, self%ngrid
         !$omp cancellation point do
         if (abort_requested) cycle

         !* -------------------------- Primitive derivatives -------------------------- *!
         if (do_timing) call self%timer%measure(20)

         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         ! Get phi derivatives
         call phi_threads(thread_slot)%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         ! Get cached phi derivatives
         phi2_r_rA = phi_threads(thread_slot)%f2_r_rA(point, anchor, owner_idx)

         ! Compute SSD on-the-fly for this point
         call lsf_threads(thread_slot)%lsf%prepare(point)

         ! Get nuclear and mixed derivatives
         call lsf_threads(thread_slot)%lsf%f3_rr_rA_screened(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         call lsf_threads(thread_slot)%lsf%f3_rrr_screened(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)

         ! Get LSF gradient magnitude
         g_norm_sq = dot_product(lsf1_r, lsf1_r)
         g_norm = sqrt(g_norm_sq)

         ! Compute Lagrangian gradient and hessian
         G_lagrangian = phi1_r - lambda_val*lsf1_r
         H_lagrangian = phi2_rr - lambda_val*lsf2_rr

         ! Compute r_i . n_i
         r_hat_dot_r = dot_product(point, self%normal0(:, igrid))

         if (do_timing) call self%timer%measure(20)

         !* ---------------------- r_i derivative ---------------------- *!
         if (do_timing) call self%timer%measure(21)

         ! Bordered KKT sensitivity solve:
         !   [H  -g] [dr/dR  ]   [b1]
         !   [g'  0] [dlambda] = [b4]
         ! where H = H_lagrangian and g = lsf1_r. Solve the full bordered
         ! system to avoid requiring H itself to be invertible.
         kkt_mat_base = 0.0_wp
         kkt_mat_base(1:3, 1:3) = H_lagrangian
         kkt_mat_base(1:3, 4) = -lsf1_r
         kkt_mat_base(4, 1:3) = lsf1_r

         ! Solve for each active atom and axis
         ! Screening: only active nodes have nonzero lsf1_rA / lsf2_r_rA;
         ! phi2_r_rA is nonzero only at owner_idx (which is always active).
         n_active = lsf_threads(thread_slot)%lsf%active_count()
         do i = 1, n_active
            active_idx(i) = lsf_threads(thread_slot)%lsf%active_atom(i)
         end do

         ! Pack all RHS into batch array
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               kkt_rhs_batch(1:3, (i - 1)*3 + iaxis) = -(phi2_r_rA(:, iaxis, iatom) &
                                                         - lambda_val*lsf2_r_rA(:, iaxis, iatom))
               kkt_rhs_batch(4, (i - 1)*3 + iaxis) = -lsf1_rA(iaxis, iatom)
            end do
         end do

         ! Single factorization + solve for all RHS
         kkt_mat = kkt_mat_base
         call lapack_gesv(4_lapack_ik, int(3*n_active, lapack_ik), kkt_mat, 4_lapack_ik, &
                          kkt_ipiv, kkt_rhs_batch, 4_lapack_ik, kkt_info)
         if (kkt_info /= 0_lapack_ik) then
            !$omp critical (compute_gradient_abort)
            if (.not. abort_requested) then
               abort_requested = .true.
               call fatal_error(error, "[Error] Bordered KKT sensitivity solve failed")
            end if
            !$omp end critical (compute_gradient_abort)
            !$omp cancel do
            cycle
         end if

         ! Unpack solutions
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               self%xyz1_rA(:, iaxis, iatom, igrid) = kkt_rhs_batch(1:3, (i - 1)*3 + iaxis)
               self%lambda1_rA(iaxis, iatom, igrid) = kkt_rhs_batch(4, (i - 1)*3 + iaxis)
            end do
         end do

         if (do_timing) call self%timer%measure(21)

         !* ---------------------- rho_i derivative (optional) ---------------------- *!
         if (allocated(self%rho1_rA)) then
            if (do_timing) call self%timer%measure(22)
            rho_vec = point - anchor
            rho_norm = sqrt(dot_product(rho_vec, rho_vec))

            if (rho_norm <= 1.0e-16_wp) then
               rho_unit = 0.0_wp
            else
               rho_unit = rho_vec/rho_norm
            end if

            ! Screening: xyz1_rA is zero for inactive atoms; delta_matrix only for owner (always active)
            do i = 1, n_active
               iatom = active_idx(i)

               ! Compute delta matrix (identity if iatom == owner_idx, zero otherwise)
               delta_matrix = 0.0_wp
               if (iatom == owner_idx) then
                  do iaxis = 1, 3
                     delta_matrix(iaxis, iaxis) = 1.0_wp
                  end do ! iaxis
               end if

               ! Compute rho derivative:rho_unit cdot ( dr/ dr_A - delta _A)
               do iaxis = 1, 3
                  self%rho1_rA(iaxis, iatom, igrid) = &
                     dot_product(rho_unit, self%xyz1_rA(:, iaxis, iatom, igrid) - delta_matrix(:, iaxis))
               end do ! iaxis
            end do ! i (active atoms)

            if (do_timing) call self%timer%measure(22)
         end if

         !* -------------------- r_iI derivative (optional) -------------------- *!
         if (allocated(self%r_iI1_rA)) then
            if (do_timing) call self%timer%measure(23)
            ! r_iI1_rA \equiv R_i when i \in I
            ! \frac{\mathbf r_{Ii}}{r_{iI}}^\top\cdot\left(\nabla_A \mathbf r_I
            ! - \nabla_A \mathbf r_i\right)

            ! r_iI1_rA \equiv \hat{r}_{Ii}^\top (\nabla_A r_I - \nabla_A r_i)
            r_iI_vec = self%mol%xyz(:, owner_idx) - point
            r_iI_norm = sqrt(dot_product(r_iI_vec, r_iI_vec))
            r_iI_vec = r_iI_vec/r_iI_norm

            ! Screening: xyz1_rA is zero for inactive; delta_matrix only for owner (always active)
            do i = 1, n_active
               iatom = active_idx(i)
               ! Reuse delta_matrix: identity for owner atom, zero otherwise
               delta_matrix = 0.0_wp
               if (iatom == owner_idx) then
                  do iaxis = 1, 3
                     delta_matrix(iaxis, iaxis) = 1.0_wp
                  end do ! iaxis
               end if

               do iaxis = 1, 3
                  self%r_iI1_rA(iaxis, iatom, igrid) = &
                     dot_product(r_iI_vec, delta_matrix(:, iaxis) - self%xyz1_rA(:, iaxis, iatom, igrid))
               end do ! iaxis
            end do ! i (active atoms)

            if (do_timing) call self%timer%measure(23)
         end if

         !* -------------------- surface normal derivative -------------------- *!
         ! n = grad(S) / ||grad(S)||
         ! dn/dr_A = (1/||grad(S)||) * [d(grad(S))/dr_A - n * (n^T * d(grad(S))/dr_A)]
         ! where d(grad(S))/dr_A = explicit + Hessian * dr/dr_A
         ! Always computed into thread-local buffer (needed by volume gradient).
         ! Stored to self%normal1_rA only when requested.
         if (do_timing) call self%timer%measure(30)

         dn_dR_buf = 0.0_wp

         ! Screening: lsf2_r_rA and xyz1_rA are zero for inactive atoms
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               ! Total derivative of grad(S) w.r.t. r_A:
               ! d(grad(S))/dr_A = ( d^2 S/ dr dr_A) + ( d^2 S/ dr^2 ) * ( dr/ dr_A)
               dg_dR = lsf2_r_rA(:, iaxis, iatom) &
                       + matmul(lsf2_rr, self%xyz1_rA(:, iaxis, iatom, igrid))

               ! Normal derivative:
               ! dn/dr_A = (1/||g||) * [dg/dr_A - n*(n^T*dg/dr_A)]
               dn_dR = (dg_dR - self%normal0(:, igrid) &
                        *dot_product(self%normal0(:, igrid), dg_dR))/g_norm

               ! Store in thread-local buffer for volume gradient
               dn_dR_buf(:, iatom, iaxis) = dn_dR

               ! Persist only when user requested normal derivatives
               if (allocated(self%normal1_rA)) &
                  self%normal1_rA(:, iatom, iaxis, igrid) = dn_dR
            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%timer%measure(30)

         !* -------------------- cpjac_scal derivative -------------------- *!
         ! 2x2 tangent-restricted approach matching projection.f90:
         !   n = g/|g|, Q = [q1,q2] from setup_tangent_frame(n)
         !   B = Q^T A Q (2x2), switch on beta2, Binv = B^{-1}
         !   tau_k = Q^T t_k, w_k = Binv*tau_k, y_k = alpha*Q*w_k
         !   J = |y1 x y2|
         ! Derivatives: dJ/dr_A via chain rule through all intermediates
         if (do_timing) call self%timer%measure(25)

         ! Coefficient alpha = 0.5 * w_a
         alpha_coeff = self%param%phi_alpha

         ! Get LSF gradient g and Hessian H at projected point
         g_vec = lsf1_r
         g_norm_sq = dot_product(g_vec, g_vec)
         g_norm = sqrt(g_norm_sq)

         ! Scalar derivative of w_f switching function w.r.t. ||g||
         call self%f_crit%eval(g_norm, f_crit0, f_crit_dS)

         ! Build A = alpha*I - lambda*H
         A_mat(1, 1) = alpha_coeff - lambda_val*lsf2_rr(1, 1)
         A_mat(1, 2) = -lambda_val*lsf2_rr(1, 2)
         A_mat(1, 3) = -lambda_val*lsf2_rr(1, 3)
         A_mat(2, 1) = -lambda_val*lsf2_rr(2, 1)
         A_mat(2, 2) = alpha_coeff - lambda_val*lsf2_rr(2, 2)
         A_mat(2, 3) = -lambda_val*lsf2_rr(2, 3)
         A_mat(3, 1) = -lambda_val*lsf2_rr(3, 1)
         A_mat(3, 2) = -lambda_val*lsf2_rr(3, 2)
         A_mat(3, 3) = alpha_coeff - lambda_val*lsf2_rr(3, 3)

         ! ---- Surface tangent frame Q = [q1, q2] from n = g/|g| ----
         n_surf = g_vec/g_norm
         call setup_tangent_frame(n_surf, q1, q2)

         ! ---- Precompute A*q1, A*q2 (reused in value and derivative) ----
         Aq1 = matmul(A_mat, q1)
         Aq2 = matmul(A_mat, q2)

         ! ---- Tangent-restricted KKT matrix B = Q^T A Q (2x2 symmetric) ----
         B11 = dot_product(q1, Aq1)
         B12 = dot_product(q1, Aq2)
         B22 = dot_product(q2, Aq2)

         ! ---- Analytic 2x2 eigenvalues of B ----
         tr_B = B11 + B22
         det_B = B11*B22 - B12*B12
         disc = 0.25_wp*tr_B*tr_B - det_B
         disc = max(disc, 0.0_wp)
         sqrt_disc = sqrt(disc)
         beta1 = 0.5_wp*tr_B + sqrt_disc
         beta2 = 0.5_wp*tr_B - sqrt_disc
         call eig_2x2_symmetric(B11, B12, B22, lambda_switch, beta1, vmin_B, vmax_B)
         u_switch = vmin_B(1)*q1 + vmin_B(2)*q2

         lambda_switch = beta2
         call self%f_foc%eval(lambda_switch, f_foc_f0, f_foc_dS)

         if (abs(det_B) <= det_B_guard) then
            !$omp critical (compute_gradient_abort)
            if (.not. abort_requested) then
               abort_requested = .true.
               call fatal_error(error, "[Error] Tangent Jacobian matrix B is singular after switching")
            end if
            !$omp end critical (compute_gradient_abort)
            !$omp cancel do
            cycle
         end if

         ! ---- Direct inverse of B ----
         Binv11 = B22/det_B
         Binv12 = -B12/det_B
         Binv22 = B11/det_B

         ! Recompute sphere tangent frame in-situ (avoids storing t_vec0)
         call setup_tangent_frame(anchor - self%mol%xyz(:, owner_idx), t1_vec, t2_vec)

         ! ---- Project sphere tangent vectors into surface tangent plane ----
         tau1(1) = dot_product(q1, t1_vec)
         tau1(2) = dot_product(q2, t1_vec)
         tau2(1) = dot_product(q1, t2_vec)
         tau2(2) = dot_product(q2, t2_vec)

         ! ---- Compute w_k = Binv * tau_k ----
         w1(1) = Binv11*tau1(1) + Binv12*tau1(2)
         w1(2) = Binv12*tau1(1) + Binv22*tau1(2)
         w2(1) = Binv11*tau2(1) + Binv12*tau2(2)
         w2(2) = Binv12*tau2(1) + Binv22*tau2(2)

         ! ---- Lift back to 3D: y_k = alpha * Q * w_k ----
         y1 = alpha_coeff*(w1(1)*q1 + w1(2)*q2)
         y2 = alpha_coeff*(w2(1)*q1 + w2(2)*q2)

         ! J = |y1 x y2|
         cross_vec(1) = y1(2)*y2(3) - y1(3)*y2(2)
         cross_vec(2) = y1(3)*y2(1) - y1(1)*y2(3)
         cross_vec(3) = y1(1)*y2(2) - y1(2)*y2(1)
         J_val = sqrt(dot_product(cross_vec, cross_vec))
         inv_J = 1.0_wp/J_val

         ! ---- Precompute Gram-Schmidt data for Q derivative ----
         ! Q = [q1, q2] is tangent frame from n_surf = g/|g|
         ! q1 built via Gram-Schmidt from e_{min_axis_surf} and n_surf
         min_axis_surf = minloc(abs(n_surf), dim=1)
         n_dot_q1_surf = n_surf(min_axis_surf)  ! = e_k . n_surf
         proj_surf = 1.0_wp - n_dot_q1_surf**2  ! = |v|^2
         v_norm_surf = sqrt(max(proj_surf, 1.0e-30_wp))

         ! Hoisted grid-level projector and product
         P_tan(:, :) = -spread(n_surf, dim=2, ncopies=3)*spread(n_surf, dim=1, ncopies=3)
         P_tan(1, 1) = P_tan(1, 1) + 1.0_wp
         P_tan(2, 2) = P_tan(2, 2) + 1.0_wp
         P_tan(3, 3) = P_tan(3, 3) + 1.0_wp
         AP_tan = matmul(A_mat, P_tan)

         ! Lebedev weight switching factor: d(w_pre * S)/dR = (S + |w_pre|*S') * d(w_pre)/dR
         if (self%param%wleb_prune_level > 0) then
            w_pre_i = self%anchor_wleb0(igrid)*self%cpjac_scal0(igrid)*self%w_f0(igrid)
            call self%f_wleb%eval(abs(w_pre_i), f_wleb_s, f_wleb_ds)
            wleb_prune_factor = f_wleb_s + abs(w_pre_i)*f_wleb_ds
         else
            wleb_prune_factor = 1.0_wp
         end if

         ! Loop over active atoms and axes to compute dJ/dr_A
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3

               ! Retrieve stored derivatives
               dlambda_val = self%lambda1_rA(iaxis, iatom, igrid)
               dr_i_dR = self%xyz1_rA(:, iaxis, iatom, igrid)

               ! dg/dr_A = explicit + Hessian * dr/dr_A
               dg_dR = lsf2_r_rA(:, iaxis, iatom) + matmul(lsf2_rr, dr_i_dR)

               ! dA/dr_A = -dlambda*H - lambda*(dH/dr_A)
               dA_dR = -dlambda_val*lsf2_rr
               do jaxis = 1, 3
                  dA_dR(:, :) = dA_dR(:, :) &
                                - lambda_val*lsf3_rrr(:, :, jaxis)*dr_i_dR(jaxis)
               end do
               dA_dR(:, :) = dA_dR(:, :) - lambda_val*lsf3_rr_rA(:, :, iaxis, iatom)

               ! ---- dn_surf/dr_A = (I - n*n^T) * dg/dr_A / |g| ----
               dn_surf_dR = (dg_dR - n_surf*dot_product(n_surf, dg_dR))/g_norm

               ! ---- dQ/dr_A: derivative of tangent frame Q = [q1, q2] ----
               ! q1 via Gram-Schmidt: v = e_k - (e_k.n)n, q1 = v/|v|
               ! dv/dr_A = -dn_surf(min_axis)*n - (e_k.n)*dn_surf
               v_tmp = -dn_surf_dR(min_axis_surf)*n_surf &
                       - n_dot_q1_surf*dn_surf_dR
               ! dq1/dr_A = (dv - q1*(q1.dv)) / |v|
               if (proj_surf > 1.0e-30_wp) then
                  dq1_dR = (v_tmp - q1*dot_product(q1, v_tmp))/v_norm_surf
               else
                  dq1_dR = 0.0_wp
               end if
               ! q2 = n x q1, dq2/dr_A = dn x q1 + n x dq1
               dq2_dR(1) = dn_surf_dR(2)*q1(3) - dn_surf_dR(3)*q1(2) &
                           + n_surf(2)*dq1_dR(3) - n_surf(3)*dq1_dR(2)
               dq2_dR(2) = dn_surf_dR(3)*q1(1) - dn_surf_dR(1)*q1(3) &
                           + n_surf(3)*dq1_dR(1) - n_surf(1)*dq1_dR(3)
               dq2_dR(3) = dn_surf_dR(1)*q1(2) - dn_surf_dR(2)*q1(1) &
                           + n_surf(1)*dq1_dR(2) - n_surf(2)*dq1_dR(1)

               ! ---- dB/dr_A: B = Q^T A Q ----
               ! Using precomputed Aq1, Aq2 and A-symmetry: q1^T*A*dq2 = Aq1^T*dq2
               dAq1 = matmul(dA_dR, q1)
               dAq2 = matmul(dA_dR, q2)
               dB11 = 2.0_wp*dot_product(dq1_dR, Aq1) + dot_product(q1, dAq1)
               dB12 = dot_product(dq1_dR, Aq2) + dot_product(dq2_dR, Aq1) &
                      + dot_product(dAq1, q2)
               dB22 = 2.0_wp*dot_product(dq2_dR, Aq2) + dot_product(q2, dAq2)

               ddet_B = dB11*B22 + B11*dB22 - 2.0_wp*B12*dB12
               dBinv11 = (dB22*det_B - B22*ddet_B)/(det_B*det_B)
               dBinv12 = (-dB12*det_B + B12*ddet_B)/(det_B*det_B)
               dBinv22 = (dB11*det_B - B11*ddet_B)/(det_B*det_B)

               ! ---- Basis-invariant d lambda_switch / dr_A ----
               ! Differentiate M = P A P; P_tan and AP_tan hoisted above atom loop.
               dP_tan(:, :) = -(spread(dn_surf_dR, dim=2, ncopies=3)*spread(n_surf, dim=1, ncopies=3) &
                                + spread(n_surf, dim=2, ncopies=3)*spread(dn_surf_dR, dim=1, ncopies=3))
               dM_tan = matmul(dP_tan, AP_tan) &
                        + matmul(P_tan, matmul(dA_dR, P_tan)) &
                        + matmul(P_tan, matmul(A_mat, dP_tan))
               dlambda_switch = dot_product(u_switch, matmul(dM_tan, u_switch))

               ! Sphere tangent frame is constant w.r.t. nuclear coordinates:
               ! n_sph = (anchor - R_I)/|...| has zero derivative for all atoms
               ! (anchor moves rigidly with owner, fixed for others).
               ! Therefore dt1/dr_A = dt2/dr_A = 0; t_vec1_rA stays at zero.

               ! ---- dtau_k/dr_A: tau_k = Q^T t_k (dt_k = 0) ----
               dtau1(1) = dot_product(dq1_dR, t1_vec)
               dtau1(2) = dot_product(dq2_dR, t1_vec)
               dtau2(1) = dot_product(dq1_dR, t2_vec)
               dtau2(2) = dot_product(dq2_dR, t2_vec)

               ! ---- dw_k/dr_A: w_k = Binv * tau_k ----
               dw1(1) = dBinv11*tau1(1) + Binv11*dtau1(1) &
                        + dBinv12*tau1(2) + Binv12*dtau1(2)
               dw1(2) = dBinv12*tau1(1) + Binv12*dtau1(1) &
                        + dBinv22*tau1(2) + Binv22*dtau1(2)
               dw2(1) = dBinv11*tau2(1) + Binv11*dtau2(1) &
                        + dBinv12*tau2(2) + Binv12*dtau2(2)
               dw2(2) = dBinv12*tau2(1) + Binv12*dtau2(1) &
                        + dBinv22*tau2(2) + Binv22*dtau2(2)

               ! ---- dy_k/dr_A: y_k = alpha * Q * w_k ----
               ! dy1 = alpha*(dw1(1)*q1 + w1(1)*dq1 + dw1(2)*q2 + w1(2)*dq2)
               dy1_dR = alpha_coeff*(dw1(1)*q1 + w1(1)*dq1_dR &
                                     + dw1(2)*q2 + w1(2)*dq2_dR)
               dy2_dR = alpha_coeff*(dw2(1)*q1 + w2(1)*dq1_dR &
                                     + dw2(2)*q2 + w2(2)*dq2_dR)

               ! d(y1 x y2)/dr_A = dy1 x y2 + y1 x dy2
               dcross_dR(1) = dy1_dR(2)*y2(3) - dy1_dR(3)*y2(2) &
                              + y1(2)*dy2_dR(3) - y1(3)*dy2_dR(2)
               dcross_dR(2) = dy1_dR(3)*y2(1) - dy1_dR(1)*y2(3) &
                              + y1(3)*dy2_dR(1) - y1(1)*dy2_dR(3)
               dcross_dR(3) = dy1_dR(1)*y2(2) - dy1_dR(2)*y2(1) &
                              + y1(1)*dy2_dR(2) - y1(2)*dy2_dR(1)

               ! dJ/dr_A = (y1 x y2) . d(y1 x y2)/dr_A / J
               dJ_dR = dot_product(cross_vec, dcross_dR)*inv_J

               self%cpjac_scal1_rA(iaxis, iatom, igrid) = dJ_dR

               ! d(w_f)/dr_A via chain rule: w_f'(||g||) * d(||g||)/dr_A
               d_gnorm = dot_product(n_surf, dg_dR)
               self%w_f1_rA(iaxis, iatom, igrid) = &
                  f_foc_f0*f_crit_dS*d_gnorm + f_crit0*f_foc_dS*dlambda_switch

               ! d(wleb)/dr_A = wbranch * (S + |w_pre|*S') * d(w_pre)/dr_A
               ! where w_pre = anchor_wleb * cpjac * w_f
               dw_pre_dR = self%anchor_wleb0(igrid)*self%w_f0(igrid)*dJ_dR &
                           + self%anchor_wleb0(igrid)*self%cpjac_scal0(igrid) &
                           *self%w_f1_rA(iaxis, iatom, igrid)
               self%wleb1_rA(iaxis, iatom, igrid) = &
                  self%wbranch(igrid)*wleb_prune_factor*dw_pre_dR

            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%timer%measure(25)

         !* ---------------------- xi_i derivative ---------------------- *!
         ! Compute derivative of Gaussian widths w.r.t. nuclear coordinates
         if (do_timing) call self%timer%measure(26)

         !> xi depends on cp_jac_scal (and derivative) *and* on anchor_xi (and derivative)
         if (allocated(self%xi1_rA)) &
            self%xi1_rA(:, :, igrid) = self%iswig%xi1_rA( &
                                       owner_idx, self%wleb(igrid), self%wleb1_rA(:, :, igrid), &
                                       active=active_idx(1:n_active))

         if (do_timing) call self%timer%measure(26)

         !* ---------------------- f_i derivative ---------------------- *!
         if (do_timing) call self%timer%measure(27)

         !> iswig switching derivatives: evaluated at anchor position
         !> Uses built-in sorted neighbor list for inner loops (early exit).
         self%f1_rA(:, :, igrid) = self%iswig%swi1_rA( &
                                   anchor, owner_idx, self%anchor_xi0(igrid), anchor_xi_local, &
                                   active=active_idx(1:n_active))

         if (do_timing) call self%timer%measure(27)

         !* ---------------------- a_i derivative ---------------------- *!
         if (do_timing) call self%timer%measure(28)
         ! Screening: f1_rA and wleb1_rA are zero for inactive atoms
         do i = 1, n_active
            iatom = active_idx(i)

            self%a_i1_rA(:, iatom, igrid) = self%radii(owner_idx)**2 &
                                            *(self%f1_rA(:, iatom, igrid)*self%wleb(igrid) &
                                              + self%f(igrid)*self%wleb1_rA(:, iatom, igrid))

         end do ! i (active atoms)

         if (do_timing) call self%timer%measure(28)

         !* ---------------------- v_i derivative ---------------------- *!
         if (do_timing) call self%timer%measure(29)

         ! New volume formula: v_i = (1/3) * a_i * (r_i . n_i)
         ! where n_i = dS / | dS| is the surface normal

         ! Screening: a_i1_rA, dn_dR_buf, xyz1_rA all zero for inactive atoms
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               ! Derivative of r_i . n_i:
               ! d(r_i . n_i)/dr_A = dr_i/dr_A . n_i + r_i . dn_i/dr_A
               ! Read normal derivative from thread-local buffer
               grad_r_hat_dot_r(iaxis) = &
                  +dot_product(self%xyz1_rA(:, iaxis, iatom, igrid), self%normal0(:, igrid)) &
                  + dot_product(point, dn_dR_buf(:, iatom, iaxis))

               ! Volume derivative: dv_i/dr_A = (1/3) * [da_i/dr_A * (r_i.n_i) + a_i * d(r_i.n_i)/dr_A]
               self%v_i1_rA(iaxis, iatom, igrid) = (1.0_wp/3.0_wp)*( &
                                                   self%a_i1_rA(iaxis, iatom, igrid)*r_hat_dot_r &
                                                   + self%a(igrid)*grad_r_hat_dot_r(iaxis))

            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%timer%measure(29)

         !* -------- Per-atom accumulation -------- *!
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               ai_val = self%a_i1_rA(iaxis, iatom, igrid)
               vi_val = self%v_i1_rA(iaxis, iatom, igrid)
               !$omp atomic
               self%asph1_rA(iaxis, owner_idx, iatom) = self%asph1_rA(iaxis, owner_idx, iatom) + ai_val
               A_tot_local(iaxis, iatom) = A_tot_local(iaxis, iatom) + ai_val
               !$omp atomic
               self%vsph1_rA(iaxis, owner_idx, iatom) = self%vsph1_rA(iaxis, owner_idx, iatom) + vi_val
               V_tot_local(iaxis, iatom) = V_tot_local(iaxis, iatom) + vi_val
            end do
         end do

      end do ! igrid
      !$omp end do

      !$omp critical (gradient_reduction)
      self%A_tot1_rA = self%A_tot1_rA + A_tot_local
      self%V_tot1_rA = self%V_tot1_rA + V_tot_local
      !$omp end critical (gradient_reduction)

      deallocate (lsf3_rrr, lsf3_rr_rA, active_idx, dn_dR_buf)
      deallocate (lsf1_rA, lsf2_r_rA, phi2_r_rA, anchor_xi_local)
      deallocate (kkt_rhs_batch, A_tot_local, V_tot_local)
      !$omp end parallel

      if (allocated(error)) return
      if (abort_requested) then
         call fatal_error(error, "Error: Gradient computation aborted. (unreachable in normal execution)")
         return
      end if

      !* ========================= Branch-weight post-pass ========================= *!
      if (do_timing) call self%timer%measure(31)

      ! TODO: This will (have to be) refactored; for now this is a slightly ugly solution (and not parallel)

      ! Assemble d(wbranch_m)/dr_A for every anchor group
      if (self%ngrid > 0 .and. any(self%branch_count(1:self%ngrid) > 1)) then
         allocate (branch_phi(maxval(self%branch_count(1:self%ngrid))), source=0.0_wp)
         allocate (branch_dphi(3*self%nsph, maxval(self%branch_count(1:self%ngrid))), &
                   source=0.0_wp)
         allocate (branch_weights(maxval(self%branch_count(1:self%ngrid))), source=0.0_wp)
         allocate (branch_dweights(3*self%nsph, maxval(self%branch_count(1:self%ngrid))), &
                   source=0.0_wp)

         igroup_start = 1
         do while (igroup_start <= self%ngrid)
            if (self%branch_count(igroup_start) <= 1) then
               igroup_start = igroup_start + 1
               cycle
            end if

            ! Extend group while anchor_id stays the same.
            igroup_end = igroup_start
            do while (igroup_end < self%ngrid)
               if (self%anchor_id(igroup_end + 1) /= self%anchor_id(igroup_start)) exit
               igroup_end = igroup_end + 1
            end do
            group_size = igroup_end - igroup_start + 1

            ! Gather phi and dphi for every branch.
            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               branch_phi(m_branch) = self%phi0(im_grid)

               owner_m = self%owner(im_grid)
               pt_m = self%xyz(:, im_grid)
               anch_m = self%anchorxyz(:, im_grid)
               ! phi = 0.5 * alpha * |r* - anch|^2, so d phi / d r = alpha * (r* - anch).
               phi1_r_m = self%param%phi_alpha*(pt_m - anch_m)

               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     ! Chain rule: phi1_r . d r*/d r_A^iaxis.
                     dphi_m = dot_product(phi1_r_m, &
                                          self%xyz1_rA(:, iaxis, iatom, im_grid))
                     ! Direct: anch moves rigidly with the owner atom, so
                     ! d phi / d R_owner = -alpha * (r* - anch) at fixed r*.
                     if (iatom == owner_m) then
                        dphi_m = dphi_m - phi1_r_m(iaxis)
                     end if
                     k_param = (iatom - 1)*3 + iaxis
                     branch_dphi(k_param, m_branch) = dphi_m
                  end do
               end do
            end do

            ! Softmax weights and their derivatives over the full branch set
            call self%branch_weight%weights_grad( &
               branch_phi(1:group_size), branch_dphi(:, 1:group_size), &
               weights=branch_weights(1:group_size), &
               dweights=branch_dweights(:, 1:group_size))

            ! Distribute d(wbranch)/dr_A back into each branch's wleb1_rA and
            ! propagate the branch-weight correction into the area and volume
            ! gradients (accumulated earlier)
            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               if (self%wbranch(im_grid) <= tiny(1.0_wp)) cycle
               owner_m = self%owner(im_grid)
               factor_m = self%wleb(im_grid)/self%wbranch(im_grid)
               area_fac_m = self%radii(owner_m)**2*self%f(im_grid)
               rn_m = dot_product(self%xyz(:, im_grid), self%normal0(:, im_grid))
               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     k_param = (iatom - 1)*3 + iaxis
                     dwleb_branch = factor_m*branch_dweights(k_param, m_branch)
                     self%wleb1_rA(iaxis, iatom, im_grid) = &
                        self%wleb1_rA(iaxis, iatom, im_grid) + dwleb_branch

                     da_branch = area_fac_m*dwleb_branch
                     dv_branch = (1.0_wp/3.0_wp)*da_branch*rn_m

                     self%a_i1_rA(iaxis, iatom, im_grid) = &
                        self%a_i1_rA(iaxis, iatom, im_grid) + da_branch
                     self%v_i1_rA(iaxis, iatom, im_grid) = &
                        self%v_i1_rA(iaxis, iatom, im_grid) + dv_branch

                     self%asph1_rA(iaxis, owner_m, iatom) = &
                        self%asph1_rA(iaxis, owner_m, iatom) + da_branch
                     self%vsph1_rA(iaxis, owner_m, iatom) = &
                        self%vsph1_rA(iaxis, owner_m, iatom) + dv_branch

                     self%A_tot1_rA(iaxis, iatom) = &
                        self%A_tot1_rA(iaxis, iatom) + da_branch
                     self%V_tot1_rA(iaxis, iatom) = &
                        self%V_tot1_rA(iaxis, iatom) + dv_branch
                  end do
               end do
            end do

            igroup_start = igroup_end + 1
         end do

         deallocate (branch_phi, branch_dphi, branch_weights, branch_dweights)
      end if

      if (do_timing) call self%timer%measure(31)

   end subroutine compute_gradient_drop

end submodule moist_cavity_drop_gradient
