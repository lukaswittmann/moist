!> Lgecacy forward-mode nuclear Jacobian of every DROP surface quantity
!>
!> This is a legacy duplicate of the new implementation as
!> [[moist_cavity_drop_derivatives_kernel:build_seed_state]] +
!> [[moist_cavity_drop_derivatives_kernel:apply_seed]], seeded with
!> `dlsf1_r = lsf2_r_rA(:,beta,A)` and `dlsf2_rr = lsf3_rr_rA(:,:,beta,A)`.
!>
!> This is left in the code base to allow for comparisons and reference
!> until possible bugs or inconsistencies in the reverse implementation
!> are resolved
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_forward
!$ use omp_lib, only: omp_get_thread_num
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_linalg, only: eig_2x2_symmetric
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_solve
   implicit none(type, external)

contains

   !> Compute first nuclear derivatives (gradients w.r.t. atomic coordinates r_A)
   !> of all DROP per-grid and per-sphere quantities
   !>
   !> For every projected grid point a bordered-KKT sensitivity system is formed
   !> from the Lagrangian Hessian $\mathbf{H}_L = \phi_{\mathbf{rr}} - \lambda
   !> S_{\mathbf{rr}}$ and the LSF gradient $\mathbf{g} = S_{\mathbf r}$,
   !>
   !> $$
   !> \begin{bmatrix} \mathbf{H}_L & -\mathbf{g} \\ \mathbf{g}^{\top} & 0 \end{bmatrix}
   !> \begin{bmatrix} \partial \mathbf{r} / \partial r_A \\ \partial \lambda / \partial r_A \end{bmatrix}
   !> =
   !> \begin{bmatrix} -\big(\phi_{\mathbf r r_A} - \lambda\, S_{\mathbf r r_A}\big) \\ -S_{r_A} \end{bmatrix},
   !> $$
   !>
   !> where $\phi_{\mathbf r r_A}$ / $S_{\mathbf r r_A}$ are mixed spatial-nuclear
   !> second derivatives and $S_{r_A}$ is the explicit nuclear derivative of $S$
   !>
   !> @param[inout] self  DROP cavity instance; reads the projected grid and
   !>                     forward quantities, then allocates and fills the `*_rA`
   !>                     first-derivative arrays
   !> @param[out]   error Error object, allocated on failure (KKT sensitivity
   !>                     solve or singular tangent Jacobian)
   module subroutine compute_gradient_drop(self, error, anchor_only)

      !> DROP cavity instance; supplies the projected grid and receives the first-derivative (`*_rA`) arrays
      class(cavity_type_drop), intent(inout) :: self

      !> Error type
      type(error_type), allocatable, intent(out) :: error

      !> Restrict each grid point's active atom to its owner (anchor motion only)
      !>
      !> For callback/isodensity LSFs the field nuclear derivatives are zero
      logical, intent(in), optional :: anchor_only

      !> Local copy of the anchor-only flag
      logical :: anchor_only_loc

      !> Per-thread LSF evaluators and projection objectives
      type(drop_worker_slots_type) :: slots

      !> Loop indices
      integer :: igrid, iatom, iaxis, jaxis, i, n_active
      integer, allocatable :: active_idx(:)

      !> OpenMP thread management
      integer :: thread_slot
      !> Thread whose timings stand in for the whole team, and the resulting
      !> per-thread gate. The timer is not thread-safe (see utils/timer.f90),
      !> so exactly one thread may touch it -- and only when the user asked for
      !> a detailed profile, since ~20 timer calls per grid point in the hot
      !> loop are not free.
      integer :: timer_ref_thread
      logical :: do_timing
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread LSF evaluation failure, handed to the latch
      type(error_type), allocatable :: lsf_error

      !> Pre-resolved timer handles for the per-grid point hot loop
      integer :: h_grad, h_prim, h_pos, h_disp, h_dist, h_norm, h_cpj, h_gw, &
                 h_sw, h_area, h_vol, h_bw

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
      real(wp) :: kkt_rhs(4, 1)
      real(wp) :: rhs_vec(4)
      integer(lapack_ik) :: kkt_info
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
      !> LSF active slot of entry `i` of `active_idx`, zero when that atom is
      !> screened away (only reachable in the anchor-only sweep, which walks
      !> every atom rather than the LSF's own active list)
      integer, allocatable :: lsf_slot(:)
      !> Nuclear partials of the atom entry currently being processed, gathered
      !> once from the active-indexed LSF outputs
      real(wp) :: s1_rA(3), s2_r_rA(3, 3), s3_rr_rA(3, 3, 3)
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

      !> Principal-curvature gradient intermediates (frame-free invariants).
      !> Grid-level (per igrid): shape-operator invariants of the LSF Hessian H;
      !> Hn = H n, adjH = adj(H), Cn = adj(H) n; T = k1+k2 = 2*KM, D = k1*k2 = KG.
      real(wp) :: Hn_curv(3), Cn_curv(3), adjH(3, 3)
      real(wp) :: trH_curv, nHn_curv, T_curv, nCn_curv, D_curv, KM_curv, disc_curv
      !> Per-(atom,axis): total nuclear derivative of H and adj(H) and the
      !> resulting curvature-invariant derivatives.
      real(wp) :: dH_curv(3, 3), dadjH(3, 3)
      real(wp) :: dtrH_c, dnHn_c, dT_c, dnCn_c, dD_c, d_disc_c
      !> Guard below which the k1/k2 split is treated as (near-)umbilic and the
      !> discriminant derivative is set to zero (individual k1/k2 derivatives are
      !> ill-defined at k1 = k2; KM and KG remain smooth).
      real(wp), parameter :: curv_disc_guard = 1.0e-10_wp

      ! Branch-weight post-pass state (serial, after main loop).
      ! Softmax weights_grad takes dphi in (nparam, nbranch) layout where
      ! nparam = 3 * nsph; we flatten (iatom, iaxis) -> (iatom - 1) * 3 + iaxis.
      integer :: igroup_start, igroup_end, group_size, m_branch, im_grid
      integer :: owner_m, k_param
      real(wp) :: pt_m(3), anch_m(3), phi1_r_m(3), dphi_m, factor_m
      real(wp) :: area_fac_m, rn_m, dwleb_branch, da_branch, dv_branch
      real(wp) :: xi_fac_m
      real(wp), allocatable :: branch_phi(:), branch_dphi(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)

      ! Thread-local xi buffer (avoids storing debug arrays)
      real(wp), allocatable :: anchor_xi_local(:, :)
      ! Scalar temps for atomic reduction (avoids gfortran aliasing issue)
      real(wp) :: ai_val, vi_val

      ! Thread-local accumulators (replace atomics for A_tot and V_tot)
      real(wp), allocatable :: A_tot_local(:, :), V_tot_local(:, :)

      anchor_only_loc = .false.
      if (present(anchor_only)) anchor_only_loc = anchor_only

      ! Initialize thread-local primitives
      ! Gradient uses third spatial derivatives (f3_rrr, f3_rr_rA);
      ! upgrade SSD storage so f3_rrr_arr is allocated before any %prepare call
      timer_ref_thread = 1
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)

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
      if (allocated(self%v1_rA)) deallocate (self%v1_rA)
      allocate (self%v1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
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
      if (allocated(self%xi1_rA)) deallocate (self%xi1_rA)
      allocate (self%xi1_rA(3, self%nsph, self%ngrid), source=0.0_wp)

      ! Principal-curvature gradients (diagnostic; gated by the same request
      ! flag as the forward compute_curvature). Mean/Gaussian curvature
      ! gradients are derived from these downstream, so they are not stored
      if (self%request%curvature) then
         if (allocated(self%k1_rA)) deallocate (self%k1_rA)
         allocate (self%k1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
         if (allocated(self%k2_rA)) deallocate (self%k2_rA)
         allocate (self%k2_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      end if

      call abort%reset()

      ! Pre-resolve the per-grid point timer handles under the "Gradients" node
      ! once, before the parallel region. Inside the loop the single timing
      ! thread uses these handles for pure-index start/stop with no lookup.
      ! The caller (get_gradient_drop) opens "Gradients" by name first, so the
      ! fine timers nest under whatever node is currently open -- regardless of
      ! how deep it sits. Fall back to a top-level "Gradients" if
      ! compute_gradient is invoked with nothing open.
      h_grad = self%ctx%timer%current()
      if (h_grad == 0) h_grad = self%ctx%timer%resolve("Gradients", 0, cat_gradient)
      h_prim = self%ctx%timer%resolve("Primitives", h_grad)
      h_pos = self%ctx%timer%resolve("Positions", h_grad)
      h_disp = self%ctx%timer%resolve("Displacement", h_grad)
      h_dist = self%ctx%timer%resolve("Distances", h_grad)
      h_norm = self%ctx%timer%resolve("Surface normal", h_grad)
      h_cpj = self%ctx%timer%resolve("CP Jacobian", h_grad)
      h_gw = self%ctx%timer%resolve("Gaussian widths", h_grad)
      h_sw = self%ctx%timer%resolve("Switching func.", h_grad)
      h_area = self%ctx%timer%resolve("Area", h_grad)
      h_vol = self%ctx%timer%resolve("Volume", h_grad)
      h_bw = self%ctx%timer%resolve("Branch weights", h_grad)

      ! Loop over all grid points (parallelized).
      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& iatom, iaxis, jaxis, point, anchor, rho_vec, rho_norm, owner_idx, lsf0, lsf1_r, lsf2_rr, &
      !$omp& lsf1_rA, lsf2_r_rA, phi0, phi1_r, phi2_rr, phi2_r_rA, lambda_val, &
      !$omp& G_lagrangian, H_lagrangian, kkt_rhs, rhs_vec, &
      !$omp& kkt_info, rho_unit, &
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
      !$omp& i, n_active, active_idx, lsf_slot, s1_rA, s2_r_rA, s3_rr_rA, &
      !$omp& f_crit0, f_crit_dS, f_foc_f0, f_foc_dS, d_gnorm, dn_dR_buf, anchor_xi_local, &
      !$omp& w_pre_i, f_wleb_s, f_wleb_ds, wleb_prune_factor, dw_pre_dR, &
      !$omp& ai_val, vi_val, kkt_rhs_batch, AP_tan, &
      !$omp& Hn_curv, Cn_curv, adjH, trH_curv, nHn_curv, T_curv, nCn_curv, &
      !$omp& D_curv, KM_curv, disc_curv, dH_curv, dadjH, dtrH_c, dnHn_c, &
      !$omp& dT_c, dnCn_c, dD_c, d_disc_c, &
      !$omp& A_tot_local, V_tot_local, lsf_error, do_timing)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1
      do_timing = thread_slot == timer_ref_thread .and. self%ctx%do_profile

      ! The LSF nuclear outputs are active-indexed and caller-owned; the buffers
      ! are sized to the atom count, which bounds `active_count()` from above.
      allocate (lsf3_rrr(3, 3, 3))
      allocate (lsf3_rr_rA(3, 3, 3, self%nsph))
      allocate (active_idx(self%nsph))
      allocate (lsf_slot(self%nsph))
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
         if (abort%requested) cycle

         !* -------------------------- Primitive derivatives -------------------------- *!
         if (do_timing) call self%ctx%timer%start(h_prim)

         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         ! Get phi derivatives
         call slots%phi(thread_slot)%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         ! Get cached phi derivatives
         phi2_r_rA = slots%phi(thread_slot)%f2_r_rA(point, anchor, owner_idx)

         ! Compute SSD on-the-fly for this point
         call slots%lsf(thread_slot)%lsf%prepare(point, lsf_error)
         if (allocated(lsf_error)) then
            call abort%latch_error(lsf_error, igrid)
            !$omp cancel do
            cycle
         end if

         ! Get nuclear and mixed derivatives
         call slots%lsf(thread_slot)%lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         call slots%lsf(thread_slot)%lsf%f3_rrr(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)

         ! Get LSF gradient magnitude
         g_norm_sq = dot_product(lsf1_r, lsf1_r)
         g_norm = sqrt(g_norm_sq)

         ! Compute Lagrangian gradient and hessian
         G_lagrangian = phi1_r - lambda_val*lsf1_r
         H_lagrangian = phi2_rr - lambda_val*lsf2_rr

         ! Compute r_i . n_i
         r_hat_dot_r = dot_product(point, self%normal0(:, igrid))

         if (do_timing) call self%ctx%timer%stop(h_prim)

         !* ---------------------- r_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_pos)

         ! Solve for each active atom and axis
         ! Screening: only active nodes have nonzero lsf1_rA / lsf2_r_rA;
         ! phi2_r_rA is nonzero only at owner_idx (which is always active)
         if (anchor_only_loc) then
            ! Walk every atom, not just the LSF-active ones; an atom screening
            ! dropped simply carries no LSF partials, which `lsf_slot` records
            n_active = self%nsph
            lsf_slot(1:n_active) = 0
            do i = 1, slots%lsf(thread_slot)%lsf%active_count()
               lsf_slot(slots%lsf(thread_slot)%lsf%active_atom(i)) = i
            end do
            do i = 1, n_active
               active_idx(i) = i
            end do
         else
            n_active = slots%lsf(thread_slot)%lsf%active_count()
            do i = 1, n_active
               active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
               lsf_slot(i) = i
            end do
         end if

         ! Pack all RHS into batch array
         do i = 1, n_active
            iatom = active_idx(i)
            call gather_lsf_partials(lsf_slot(i), lsf1_rA, lsf2_r_rA, lsf3_rr_rA, &
                                     s1_rA, s2_r_rA, s3_rr_rA)
            do iaxis = 1, 3
               kkt_rhs_batch(1:3, (i - 1)*3 + iaxis) = -(phi2_r_rA(:, iaxis, iatom) &
                                                         - lambda_val*s2_r_rA(:, iaxis))
               kkt_rhs_batch(4, (i - 1)*3 + iaxis) = -s1_rA(iaxis)
            end do
         end do

         ! Single factorization + solve for all RHS
         call drop_kkt_solve(H_lagrangian, lsf1_r, kkt_rhs_batch(:, 1:3*n_active), kkt_info)
         if (kkt_info /= 0_lapack_ik) then
            call abort%latch_message("[Error] Bordered KKT sensitivity solve failed", igrid)
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

         if (do_timing) call self%ctx%timer%stop(h_pos)

         !* ---------------------- rho_i derivative (optional) ---------------------- *!
         if (allocated(self%rho1_rA)) then
            if (do_timing) call self%ctx%timer%start(h_disp)
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

            if (do_timing) call self%ctx%timer%stop(h_disp)
         end if

         !* -------------------- r_iI derivative (optional) -------------------- *!
         if (allocated(self%r_iI1_rA)) then
            if (do_timing) call self%ctx%timer%start(h_dist)
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

            if (do_timing) call self%ctx%timer%stop(h_dist)
         end if

         !* -------------------- surface normal derivative -------------------- *!
         ! n = grad(S) / ||grad(S)||
         ! dn/dr_A = (1/||grad(S)||) * [d(grad(S))/dr_A - n * (n^T * d(grad(S))/dr_A)]
         ! where d(grad(S))/dr_A = explicit + Hessian * dr/dr_A
         ! Always computed into thread-local buffer (needed by volume gradient).
         ! Stored to self%normal1_rA only when requested.
         if (do_timing) call self%ctx%timer%start(h_norm)

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
               if (allocated(self%normal1_rA)) then
                  self%normal1_rA(:, iatom, iaxis, igrid) = dn_dR
               end if
            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%ctx%timer%stop(h_norm)

         !* -------------------- cpjac_scal derivative -------------------- *!
         ! 2x2 tangent-restricted approach matching projection.f90:
         !   n = g/|g|, Q = [q1,q2] from setup_tangent_frame(n)
         !   B = Q^T A Q (2x2), switch on beta2, Binv = B^{-1}
         !   tau_k = Q^T t_k, w_k = Binv*tau_k, y_k = alpha*Q*w_k
         !   J = |y1 x y2|
         ! Derivatives: dJ/dr_A via chain rule through all intermediates
         if (do_timing) call self%ctx%timer%start(h_cpj)

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

         ! Surface tangent frame Q = [q1, q2] from n = g/|g|
         n_surf = g_vec/g_norm
         call setup_tangent_frame(n_surf, q1, q2)

         ! Precompute A*q1, A*q2 (reused in value and derivative)
         Aq1 = matmul(A_mat, q1)
         Aq2 = matmul(A_mat, q2)

         ! Tangent-restricted KKT matrix B = Q^T A Q (2x2 symmetric)
         B11 = dot_product(q1, Aq1)
         B12 = dot_product(q1, Aq2)
         B22 = dot_product(q2, Aq2)

         ! Analytic 2x2 eigenvalues of B
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
            call abort%latch_message( &
               "[Error] Tangent Jacobian matrix B is singular after switching", igrid)
            !$omp cancel do
            cycle
         end if

         ! Direct inverse of B
         Binv11 = B22/det_B
         Binv12 = -B12/det_B
         Binv22 = B11/det_B

         ! Recompute sphere tangent frame in-situ (avoids storing t_vec0)
         call setup_tangent_frame(anchor - self%mol%xyz(:, owner_idx), t1_vec, t2_vec)

         ! Project sphere tangent vectors into surface tangent plane
         tau1(1) = dot_product(q1, t1_vec)
         tau1(2) = dot_product(q2, t1_vec)
         tau2(1) = dot_product(q1, t2_vec)
         tau2(2) = dot_product(q2, t2_vec)

         ! Compute w_k = Binv * tau_k
         w1(1) = Binv11*tau1(1) + Binv12*tau1(2)
         w1(2) = Binv12*tau1(1) + Binv22*tau1(2)
         w2(1) = Binv11*tau2(1) + Binv12*tau2(2)
         w2(2) = Binv12*tau2(1) + Binv22*tau2(2)

         ! Lift back to 3D: y_k = alpha * Q * w_k
         y1 = alpha_coeff*(w1(1)*q1 + w1(2)*q2)
         y2 = alpha_coeff*(w2(1)*q1 + w2(2)*q2)

         ! J = |y1 x y2|
         cross_vec(1) = y1(2)*y2(3) - y1(3)*y2(2)
         cross_vec(2) = y1(3)*y2(1) - y1(1)*y2(3)
         cross_vec(3) = y1(1)*y2(2) - y1(2)*y2(1)
         J_val = sqrt(dot_product(cross_vec, cross_vec))
         inv_J = 1.0_wp/J_val

         ! Precompute Gram-Schmidt data for Q derivative
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

         ! Grid-level principal-curvature invariants (frame independent)
         ! Shape operator of the level set: eigenvalues on the tangent plane are
         ! the principal curvatures k1 >= k2. Using invariants avoids the
         ! discontinuous tangent-frame choice of the forward compute_curvature:
         !   T = k1 + k2 = 2*KM = (tr H - n^T H n)/|g|
         !   D = k1 * k2 = KG   = (n^T adj(H) n)/|g|^2   (Goldman implicit-surface)
         !   k1,k2 = KM +/- sqrt(KM^2 - KG)
         if (allocated(self%k1_rA)) then
            Hn_curv = matmul(lsf2_rr, n_surf)
            trH_curv = lsf2_rr(1, 1) + lsf2_rr(2, 2) + lsf2_rr(3, 3)
            nHn_curv = dot_product(n_surf, Hn_curv)
            T_curv = (trH_curv - nHn_curv)/g_norm

            ! Adjugate (cofactor matrix) of the symmetric Hessian H = lsf2_rr
            adjH(1, 1) = lsf2_rr(2, 2)*lsf2_rr(3, 3) - lsf2_rr(2, 3)*lsf2_rr(2, 3)
            adjH(2, 2) = lsf2_rr(1, 1)*lsf2_rr(3, 3) - lsf2_rr(1, 3)*lsf2_rr(1, 3)
            adjH(3, 3) = lsf2_rr(1, 1)*lsf2_rr(2, 2) - lsf2_rr(1, 2)*lsf2_rr(1, 2)
            adjH(1, 2) = lsf2_rr(1, 3)*lsf2_rr(2, 3) - lsf2_rr(1, 2)*lsf2_rr(3, 3)
            adjH(1, 3) = lsf2_rr(1, 2)*lsf2_rr(2, 3) - lsf2_rr(2, 2)*lsf2_rr(1, 3)
            adjH(2, 3) = lsf2_rr(1, 2)*lsf2_rr(1, 3) - lsf2_rr(1, 1)*lsf2_rr(2, 3)
            adjH(2, 1) = adjH(1, 2)
            adjH(3, 1) = adjH(1, 3)
            adjH(3, 2) = adjH(2, 3)

            Cn_curv = matmul(adjH, n_surf)
            nCn_curv = dot_product(n_surf, Cn_curv)
            D_curv = nCn_curv/g_norm_sq
            KM_curv = 0.5_wp*T_curv
            disc_curv = sqrt(max(KM_curv*KM_curv - D_curv, 0.0_wp))
         end if

         ! Loop over active atoms and axes to compute dJ/dr_A
         do i = 1, n_active
            iatom = active_idx(i)
            call gather_lsf_partials(lsf_slot(i), lsf1_rA, lsf2_r_rA, lsf3_rr_rA, &
                                     s1_rA, s2_r_rA, s3_rr_rA)
            do iaxis = 1, 3

               ! Retrieve stored derivatives
               dlambda_val = self%lambda1_rA(iaxis, iatom, igrid)
               dr_i_dR = self%xyz1_rA(:, iaxis, iatom, igrid)

               ! dg/dr_A = explicit + Hessian * dr/dr_A
               dg_dR = s2_r_rA(:, iaxis) + matmul(lsf2_rr, dr_i_dR)

               ! dA/dr_A = -dlambda*H - lambda*(dH/dr_A)
               dA_dR = -dlambda_val*lsf2_rr
               do jaxis = 1, 3
                  dA_dR(:, :) = dA_dR(:, :) &
                                - lambda_val*lsf3_rrr(:, :, jaxis)*dr_i_dR(jaxis)
               end do
               dA_dR(:, :) = dA_dR(:, :) - lambda_val*s3_rr_rA(:, :, iaxis)

               ! dn_surf/dr_A = (I - n*n^T) * dg/dr_A / |g|
               dn_surf_dR = (dg_dR - n_surf*dot_product(n_surf, dg_dR))/g_norm

               ! dQ/dr_A: derivative of tangent frame Q = [q1, q2]
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

               ! dB/dr_A: B = Q^T A Q
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

               ! Basis-invariant d lambda_switch / dr_A
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

               ! dtau_k/dr_A: tau_k = Q^T t_k (dt_k = 0)
               dtau1(1) = dot_product(dq1_dR, t1_vec)
               dtau1(2) = dot_product(dq2_dR, t1_vec)
               dtau2(1) = dot_product(dq1_dR, t2_vec)
               dtau2(2) = dot_product(dq2_dR, t2_vec)

               ! dw_k/dr_A: w_k = Binv * tau_k
               dw1(1) = dBinv11*tau1(1) + Binv11*dtau1(1) &
                        + dBinv12*tau1(2) + Binv12*dtau1(2)
               dw1(2) = dBinv12*tau1(1) + Binv12*dtau1(1) &
                        + dBinv22*tau1(2) + Binv22*dtau1(2)
               dw2(1) = dBinv11*tau2(1) + Binv11*dtau2(1) &
                        + dBinv12*tau2(2) + Binv12*dtau2(2)
               dw2(2) = dBinv12*tau2(1) + Binv12*dtau2(1) &
                        + dBinv22*tau2(2) + Binv22*dtau2(2)

               ! dy_k/dr_A: y_k = alpha * Q * w_k
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

               ! Principal / mean / Gaussian curvature derivatives
               ! Total nuclear derivative of the LSF Hessian at the moving
               ! projected point:
               !   dH/dr_A = d^3S/(dr dr dR_A) + (d^3S/dr^3) . dr/dr_A
               if (allocated(self%k1_rA) .and. allocated(self%k2_rA)) then
                  dH_curv = s3_rr_rA(:, :, iaxis)
                  do jaxis = 1, 3
                     dH_curv(:, :) = dH_curv(:, :) &
                                     + lsf3_rrr(:, :, jaxis)*dr_i_dR(jaxis)
                  end do

                  ! dT = d[(tr H - n^T H n)/|g|]
                  ! d(n^T H n) = 2 (dn . H n) + n^T dH n   (H symmetric)
                  dtrH_c = dH_curv(1, 1) + dH_curv(2, 2) + dH_curv(3, 3)
                  dnHn_c = 2.0_wp*dot_product(dn_surf_dR, Hn_curv) &
                           + dot_product(n_surf, matmul(dH_curv, n_surf))
                  dT_c = (dtrH_c - dnHn_c)/g_norm - T_curv*d_gnorm/g_norm

                  ! d(adj H) via product rule (symmetric)
                  dadjH(1, 1) = dH_curv(2, 2)*lsf2_rr(3, 3) + lsf2_rr(2, 2)*dH_curv(3, 3) &
                                - 2.0_wp*lsf2_rr(2, 3)*dH_curv(2, 3)
                  dadjH(2, 2) = dH_curv(1, 1)*lsf2_rr(3, 3) + lsf2_rr(1, 1)*dH_curv(3, 3) &
                                - 2.0_wp*lsf2_rr(1, 3)*dH_curv(1, 3)
                  dadjH(3, 3) = dH_curv(1, 1)*lsf2_rr(2, 2) + lsf2_rr(1, 1)*dH_curv(2, 2) &
                                - 2.0_wp*lsf2_rr(1, 2)*dH_curv(1, 2)
                  dadjH(1, 2) = dH_curv(1, 3)*lsf2_rr(2, 3) + lsf2_rr(1, 3)*dH_curv(2, 3) &
                                - dH_curv(1, 2)*lsf2_rr(3, 3) - lsf2_rr(1, 2)*dH_curv(3, 3)
                  dadjH(1, 3) = dH_curv(1, 2)*lsf2_rr(2, 3) + lsf2_rr(1, 2)*dH_curv(2, 3) &
                                - dH_curv(2, 2)*lsf2_rr(1, 3) - lsf2_rr(2, 2)*dH_curv(1, 3)
                  dadjH(2, 3) = dH_curv(1, 2)*lsf2_rr(1, 3) + lsf2_rr(1, 2)*dH_curv(1, 3) &
                                - dH_curv(1, 1)*lsf2_rr(2, 3) - lsf2_rr(1, 1)*dH_curv(2, 3)
                  dadjH(2, 1) = dadjH(1, 2)
                  dadjH(3, 1) = dadjH(1, 3)
                  dadjH(3, 2) = dadjH(2, 3)

                  ! dD = d[(n^T adj(H) n)/|g|^2]
                  ! d(n^T adj n) = 2 (dn . adj(H) n) + n^T d(adj H) n
                  dnCn_c = 2.0_wp*dot_product(dn_surf_dR, Cn_curv) &
                           + dot_product(n_surf, matmul(dadjH, n_surf))
                  dD_c = dnCn_c/g_norm_sq - 2.0_wp*D_curv*d_gnorm/g_norm

                  ! Discriminant split: disc^2 = KM^2 - KG, so
                  ! 2 disc d(disc) = 2 KM dKM - dKG = (T/2) dT - dD.
                  if (disc_curv > curv_disc_guard) then
                     d_disc_c = (KM_curv*dT_c - dD_c)/(2.0_wp*disc_curv)
                  else
                     d_disc_c = 0.0_wp
                  end if

                  ! Principal-curvature gradients. The mean/Gaussian curvature
                  ! gradients are intentionally not stored: downstream they are
                  !   dKM = (k1_rA + k2_rA)/2,  dKG = k2*k1_rA + k1*k2_rA.
                  self%k1_rA(iaxis, iatom, igrid) = 0.5_wp*dT_c + d_disc_c
                  self%k2_rA(iaxis, iatom, igrid) = 0.5_wp*dT_c - d_disc_c
               end if

            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%ctx%timer%stop(h_cpj)

         !* ---------------------- xi_i derivative ---------------------- *!
         ! Compute derivative of Gaussian widths w.r.t. nuclear coordinates
         if (do_timing) call self%ctx%timer%start(h_gw)

         !> xi depends on cp_jac_scal (and derivative) *and* on anchor_xi (and derivative)
         if (allocated(self%xi1_rA)) then
            self%xi1_rA(:, :, igrid) = self%iswig%xi1_rA( &
                                       owner_idx, self%wleb(igrid), self%wleb1_rA(:, :, igrid), &
                                       active=active_idx(1:n_active))
         end if

         if (do_timing) call self%ctx%timer%stop(h_gw)

         !* ---------------------- f_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_sw)

         !> iswig switching derivatives: evaluated at anchor position
         !> Uses built-in sorted neighbor list for inner loops (early exit).
         self%f1_rA(:, :, igrid) = self%iswig%swi1_rA( &
                                   anchor, owner_idx, self%anchor_xi0(igrid), anchor_xi_local, &
                                   active=active_idx(1:n_active))

         if (do_timing) call self%ctx%timer%stop(h_sw)

         !* ---------------------- a_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_area)
         ! Screening: f1_rA and wleb1_rA are zero for inactive atoms
         do i = 1, n_active
            iatom = active_idx(i)

            self%a_i1_rA(:, iatom, igrid) = self%radii(owner_idx)**2 &
                                            *(self%f1_rA(:, iatom, igrid)*self%wleb(igrid) &
                                              + self%f(igrid)*self%wleb1_rA(:, iatom, igrid))

         end do ! i (active atoms)

         if (do_timing) call self%ctx%timer%stop(h_area)

         !* ---------------------- v_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_vol)

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
               self%v1_rA(iaxis, iatom, igrid) = (1.0_wp/3.0_wp)*( &
                                                 self%a_i1_rA(iaxis, iatom, igrid)*r_hat_dot_r &
                                                 + self%a(igrid)*grad_r_hat_dot_r(iaxis))

            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%ctx%timer%stop(h_vol)

         !* -------- Per-atom accumulation -------- *!
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               ai_val = self%a_i1_rA(iaxis, iatom, igrid)
               vi_val = self%v1_rA(iaxis, iatom, igrid)
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

      deallocate (lsf3_rrr, lsf3_rr_rA, active_idx, lsf_slot, dn_dR_buf)
      deallocate (lsf1_rA, lsf2_r_rA, phi2_r_rA, anchor_xi_local)
      deallocate (kkt_rhs_batch, A_tot_local, V_tot_local)
      !$omp end parallel

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call fatal_error(error, &
                             "Error: Gradient computation aborted. (unreachable in normal execution)")
         end if
         return
      end if

      !* ========================= Branch-weight post-pass ========================= *!
      call self%ctx%timer%start(h_bw)

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
               ! xi = swx/(R*sqrt(wleb)) -> branch weight of dxi/dwleb = -xi/(2*wleb)
               if (self%wleb(im_grid) > tiny(1.0_wp)) then
                  xi_fac_m = -0.5_wp*self%xi0(im_grid)/self%wleb(im_grid)
               else
                  xi_fac_m = 0.0_wp
               end if
               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     k_param = (iatom - 1)*3 + iaxis
                     dwleb_branch = factor_m*branch_dweights(k_param, m_branch)
                     self%wleb1_rA(iaxis, iatom, im_grid) = &
                        self%wleb1_rA(iaxis, iatom, im_grid) + dwleb_branch

                     if (allocated(self%xi1_rA)) then
                        self%xi1_rA(iaxis, iatom, im_grid) = &
                           self%xi1_rA(iaxis, iatom, im_grid) + xi_fac_m*dwleb_branch
                     end if

                     da_branch = area_fac_m*dwleb_branch
                     dv_branch = (1.0_wp/3.0_wp)*da_branch*rn_m

                     self%a_i1_rA(iaxis, iatom, im_grid) = &
                        self%a_i1_rA(iaxis, iatom, im_grid) + da_branch
                     self%v1_rA(iaxis, iatom, im_grid) = &
                        self%v1_rA(iaxis, iatom, im_grid) + dv_branch

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

      call self%ctx%timer%stop(h_bw)

   end subroutine compute_gradient_drop

   !> Compute anchor-only nuclear derivatives for the DROP cavity
   !>
   !> For a isodensity LSF the level set field has zero nuclear derivatives
   !> (in first order), so only the reference system motion derivatives
   !> are needed
   !>
   !> @param[inout] self   DROP cavity instance (must hold a projected grid)
   !> @param[out]   error  Error object, allocated on failure
   module subroutine compute_anchor_gradient(self, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Error object
      type(error_type), allocatable, intent(out) :: error

      call self%compute_gradient_drop(error, anchor_only=.true.)
   end subroutine compute_anchor_gradient

   !> Gather one atom's LSF nuclear partials out of the active-indexed outputs
   !>
   !> The LSF hands its nuclear tensors back indexed by active slot, so an atom
   !> that screening dropped has no slot at all. A missing slot is not an error:
   !> its partials are exactly zero, which is what the sweep over *all* atoms
   !> (`anchor_only`) needs and what the sweep over the active list never asks for.
   !>
   !> @param[in]  islot      LSF active slot, or zero when the atom is inactive
   !> @param[in]  lsf1_rA    Active-indexed dS/dR_A
   !> @param[in]  lsf2_r_rA  Active-indexed d^2S/(dr dR_A)
   !> @param[in]  lsf3_rr_rA Active-indexed d^3S/(dr^2 dR_A)
   !> @param[out] s1_rA      This atom's dS/dR_A [3]
   !> @param[out] s2_r_rA    This atom's d^2S/(dr dR_A) [3, 3]
   !> @param[out] s3_rr_rA   This atom's d^3S/(dr^2 dR_A) [3, 3, 3]
   pure subroutine gather_lsf_partials(islot, lsf1_rA, lsf2_r_rA, lsf3_rr_rA, &
                                       s1_rA, s2_r_rA, s3_rr_rA)
      !> LSF active slot, or zero
      integer, intent(in) :: islot
      !> Active-indexed nuclear gradient
      real(wp), intent(in) :: lsf1_rA(:, :)
      !> Active-indexed mixed second derivative
      real(wp), intent(in) :: lsf2_r_rA(:, :, :)
      !> Active-indexed mixed third derivative
      real(wp), intent(in) :: lsf3_rr_rA(:, :, :, :)
      !> Gathered nuclear gradient
      real(wp), intent(out) :: s1_rA(3)
      !> Gathered mixed second derivative
      real(wp), intent(out) :: s2_r_rA(3, 3)
      !> Gathered mixed third derivative
      real(wp), intent(out) :: s3_rr_rA(3, 3, 3)

      if (islot > 0) then
         s1_rA = lsf1_rA(:, islot)
         s2_r_rA = lsf2_r_rA(:, :, islot)
         s3_rr_rA = lsf3_rr_rA(:, :, :, islot)
      else
         s1_rA = 0.0_wp
         s2_r_rA = 0.0_wp
         s3_rr_rA = 0.0_wp
      end if
   end subroutine gather_lsf_partials

end submodule moist_cavity_drop_derivatives_forward
