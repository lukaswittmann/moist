!> Legacy forward-mode nuclear Jacobian of every DROP surface quantity
!>
!> - drives the shared per-point kernel, one seed per (atom, axis):
!>   [[moist_cavity_drop_derivatives_kernel:build_seed_state]] once per grid point,
!>   [[moist_cavity_drop_derivatives_kernel:apply_seed]] with
!>   `dlsf1_r = lsf2_r_rA(:,beta,A)`, `dlsf2_rr` the unpacked `lsf3_rr_rA(:,beta,A)`
!>   and the point motion from the bordered KKT solve
!> - own algebra only where the kernel has no output: point motion, displacement,
!>   distance, switching, area, volume and the branch-weight post-pass
!> - `3 * n_active` seeds per grid point and `O(N_sph * N_grid)` memory; kept as
!>   reference for the reverse paths, which contract the same kernel
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_forward
!$ use omp_lib, only: omp_get_thread_num
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_packed, only: packed_idx33
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, build_seed_state, apply_seed, seed_state_ok, &
      & seed_state_singular_kkt
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_solve, degenerate_point_error
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_type, only: add_volume_origin_gradient
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
      integer :: igrid, iatom, iaxis, jaxis, kaxis, i, n_active
      integer, allocatable :: active_idx(:)

      !> OpenMP thread management
      integer :: thread_slot
      !> Thread whose timings stand in for the whole team, and the resulting
      !> per-thread gate
      !>
      !> - the timer is not thread-safe (see utils/timer.f90), so exactly one
      !>   thread may touch it
      !> - only when the user asked for a detailed profile, since ~20 timer calls
      !>   per grid point in the hot loop are not free
      integer :: timer_ref_thread
      logical :: do_timing
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread LSF evaluation failure, handed to the latch
      type(error_type), allocatable :: lsf_error

      !> Pre-resolved timer handles for the per-grid point hot loop
      integer :: h_grad, h_prim, h_pos, h_disp, h_dist, h_seed, h_sw, h_area, h_vol, h_bw

      !> Grid point data
      real(wp) :: point(3), anchor(3), rho_vec(3), rho_norm
      integer :: owner_idx

      !> Level-set value; the spatial jet is written straight into the seed state
      real(wp) :: lsf0
      !> Active-indexed LSF nuclear partials, `lsf3_rr_rA` packed d^3S/(dr^2 dR_A)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :)

      !> Phi derivatives, packed Hessian
      real(wp) :: phi0, phi1_r(3), phi2_rr(6)
      real(wp), allocatable :: phi2_r_rA(:, :, :)

      !> KKT system, packed Lagrangian Hessian
      real(wp) :: lambda_val
      real(wp) :: H_lagrangian(6)
      integer(lapack_ik) :: kkt_info
      real(wp), allocatable :: kkt_rhs_batch(:, :)

      !> Per-point kernel state, its degeneracy status and the response to one seed
      type(drop_seed_state_type) :: state
      integer :: status
      type(drop_seed_result_type) :: res
      !> Whether the principal-curvature gradients are filled
      logical :: want_curvature

      !> swi: Rho derivatives
      real(wp) :: rho_unit(3), delta_matrix(3, 3)

      !> Point derivatives
      real(wp) :: r_iI_vec(3), r_iI_norm

      !> Volume derivatives
      real(wp) :: r_hat_dot_r, grad_r_hat_dot_r(3)
      !> Volume origin moving with the molecule
      real(wp) :: origin(3)

      !> Thread-local buffer for normal derivatives (used by volume gradient)
      real(wp), allocatable :: dn_dR_buf(:, :, :)

      !> LSF active slot of entry `i` of `active_idx`
      integer, allocatable :: lsf_slot(:)
      !> Nuclear partials of the atom entry currently being processed, s3_rr_rA packed
      real(wp) :: s1_rA(3), s2_r_rA(3, 3), s3_rr_rA(6, 3)
      !> Field seed of `grad^2 S` at fixed `r`: unpacked `s3_rr_rA(:, iaxis)`
      real(wp) :: dlsf2_rr(3, 3)

      ! Branch-weight post-pass state (serial, after main loop)
      ! Softmax weights_grad takes dphi in (nparam, nbranch) layout with
      ! nparam = 3 * nsph, flattening (iatom, iaxis) -> (iatom - 1) * 3 + iaxis
      integer :: igroup_start, igroup_end, group_size, m_branch, im_grid, ia_group, max_group
      integer :: owner_m, k_param
      real(wp) :: pt_m(3), anch_m(3), phi1_r_m(3), dphi_m, factor_m
      real(wp) :: area_fac_m, rn_m, dwleb_branch, da_branch, dv_branch
      real(wp) :: xi_fac_m
      real(wp), allocatable :: branch_phi(:), branch_dphi(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)

      ! Thread-local switching-gradient scratch
      type(iswig_workspace_type) :: iswig_work
      ! Scalar temps for atomic reduction (avoids gfortran aliasing issue)
      real(wp) :: ai_val, vi_val

      ! Thread-local accumulators (replace atomics for A_tot and V_tot)
      real(wp), allocatable :: A_tot_local(:, :), V_tot_local(:, :)

      ! Volume-origin term dv_i/dc and its per-sphere sums (shared and thread-local)
      real(wp) :: dv_dc(3)
      real(wp), allocatable :: dv_dc_sph(:, :), dv_dc_sph_local(:, :)
      integer :: isph

      call self%require_context(error)
      if (allocated(error)) return

      anchor_only_loc = .false.
      if (present(anchor_only)) anchor_only_loc = anchor_only

      ! Thread-local primitives; the gradient reads third spatial derivatives
      timer_ref_thread = 1
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)

      ! Allocate gradient arrays
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      allocate (self%xyz1_rA(3, 3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%lambda1_rA)) deallocate (self%lambda1_rA)
      allocate (self%lambda1_rA(3, self%nsph, self%ngrid), source=0.0_wp)

      ! Optional derivative arrays (gated by the property flags)
      if (self%param%do_rho) then
         if (allocated(self%rho1_rA)) deallocate (self%rho1_rA)
         allocate (self%rho1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      end if
      if (self%param%do_r_iI) then
         if (allocated(self%r_iI1_rA)) deallocate (self%r_iI1_rA)
         allocate (self%r_iI1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      end if
      if (self%param%do_normal) then
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

      ! Principal-curvature gradients, diagnostic and gated by the same property
      ! flag as the forward compute_curvature; mean/Gaussian curvature gradients
      ! are derived from these downstream, so they are not stored
      if (self%param%do_curvature) then
         if (allocated(self%k1_rA)) deallocate (self%k1_rA)
         allocate (self%k1_rA(3, self%nsph, self%ngrid), source=0.0_wp)
         if (allocated(self%k2_rA)) deallocate (self%k2_rA)
         allocate (self%k2_rA(3, self%nsph, self%ngrid), source=0.0_wp)
      end if

      want_curvature = allocated(self%k1_rA) .and. allocated(self%k2_rA)

      call abort%reset()

      origin = self%get_volume_origin()
      allocate (dv_dc_sph(3, self%nsph), source=0.0_wp)

      ! Pre-resolve the per-grid point timer handles under the "Gradients" node
      ! once, before the parallel region:
      !
      ! - inside the loop the single timing thread uses these handles for
      !   pure-index start/stop with no lookup
      ! - the caller (get_gradient_drop) opens "Gradients" by name first, so the
      !   fine timers nest under whatever node is currently open, regardless of
      !   how deep it sits
      ! - falls back to a top-level "Gradients" if compute_gradient is invoked
      !   with nothing open
      h_grad = self%ctx%timer%current()
      if (h_grad == 0) h_grad = self%ctx%timer%resolve("Gradients", 0, cat_gradient)
      h_prim = self%ctx%timer%resolve("Primitives", h_grad)
      h_pos = self%ctx%timer%resolve("Positions", h_grad)
      h_disp = self%ctx%timer%resolve("Displacement", h_grad)
      h_dist = self%ctx%timer%resolve("Distances", h_grad)
      h_seed = self%ctx%timer%resolve("Seed responses", h_grad)
      h_sw = self%ctx%timer%resolve("Switching func.", h_grad)
      h_area = self%ctx%timer%resolve("Area", h_grad)
      h_vol = self%ctx%timer%resolve("Volume", h_grad)
      h_bw = self%ctx%timer%resolve("Branch weights", h_grad)

      ! Loop over all grid points (parallelized)
      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& iatom, iaxis, jaxis, kaxis, point, anchor, rho_vec, rho_norm, owner_idx, lsf0, &
      !$omp& lsf1_rA, lsf2_r_rA, lsf3_rr_rA, phi0, phi1_r, phi2_rr, phi2_r_rA, lambda_val, &
      !$omp& H_lagrangian, kkt_info, kkt_rhs_batch, state, status, res, dlsf2_rr, &
      !$omp& rho_unit, delta_matrix, r_iI_vec, r_iI_norm, r_hat_dot_r, grad_r_hat_dot_r, &
      !$omp& i, n_active, active_idx, lsf_slot, s1_rA, s2_r_rA, s3_rr_rA, &
      !$omp& dn_dR_buf, iswig_work, ai_val, vi_val, &
      !$omp& A_tot_local, V_tot_local, dv_dc, dv_dc_sph_local, lsf_error, do_timing)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1
      do_timing = thread_slot == timer_ref_thread .and. self%ctx%do_profile

      ! The LSF nuclear outputs are active-indexed and caller-owned; the buffers
      ! are sized to the atom count, which bounds `active_count()` from above
      allocate (lsf3_rr_rA(6, 3, self%nsph))
      allocate (active_idx(self%nsph))
      allocate (lsf_slot(self%nsph))
      allocate (dn_dR_buf(3, self%nsph, 3))
      allocate (lsf1_rA(3, self%nsph))
      allocate (lsf2_r_rA(3, 3, self%nsph))
      allocate (phi2_r_rA(3, 3, self%nsph))
      call self%iswig%new_workspace(iswig_work)
      allocate (kkt_rhs_batch(4, 3*self%nsph))
      allocate (A_tot_local(3, self%nsph), source=0.0_wp)
      allocate (V_tot_local(3, self%nsph), source=0.0_wp)
      allocate (dv_dc_sph_local(3, self%nsph), source=0.0_wp)

      !$omp do schedule(dynamic)
      do igrid = 1, self%ngrid
         !$omp cancellation point do
         if (abort%requested) cycle

         !* -------------------------- Primitive derivatives -------------------------- *!
         if (do_timing) call self%ctx%timer%start(h_prim)

         point = self%xyz(:, igrid)
         anchor = self%anchor%xyz(:, self%anchor_of_row(igrid))
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
         call slots%lsf(thread_slot)%lsf%f3_rrr(lsf0, state%lsf1_r, state%lsf2_rr, state%lsf3_rrr)

         ! Per-point kernel state; a degenerate point is latched
         call fill_seed_state(self, igrid, want_curvature, state)
         state%lambda_val = lambda_val
         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)
         if (status /= seed_state_ok) then
            call abort%latch_status(status, igrid)
            !$omp cancel do
            cycle
         end if

         ! Compute Lagrangian hessian
         H_lagrangian = phi2_rr - lambda_val*state%lsf2_rr

         ! Compute (r_i - c) . n_i, c the volume origin
         r_hat_dot_r = dot_product(point - origin, self%normal0(:, igrid))

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
         call drop_kkt_solve(H_lagrangian, state%lsf1_r, kkt_rhs_batch(:, 1:3*n_active), kkt_info)
         if (kkt_info /= 0_lapack_ik) then
            call abort%latch_status(seed_state_singular_kkt, igrid)
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

         !* ------------------------------ Seed responses ----------------------------- *!
         ! One kernel seed per (atom, axis): the field moves by the LSF nuclear
         ! partials at fixed r, the point and its multiplier by the KKT solution
         if (do_timing) call self%ctx%timer%start(h_seed)

         dn_dR_buf = 0.0_wp

         ! Screening: the LSF partials are active-slot indexed, read through `lsf_slot`;
         ! an atom outside the point's reach carries an exactly vanishing field seed
         do i = 1, n_active
            iatom = active_idx(i)
            call gather_lsf_partials(lsf_slot(i), lsf1_rA, lsf2_r_rA, lsf3_rr_rA, &
                                     s1_rA, s2_r_rA, s3_rr_rA)
            do iaxis = 1, 3
               do kaxis = 1, 3
                  do jaxis = 1, 3
                     dlsf2_rr(jaxis, kaxis) = s3_rr_rA(packed_idx33(jaxis, kaxis), iaxis)
                  end do
               end do
               call apply_seed(state, s2_r_rA(:, iaxis), dlsf2_rr, &
                               self%xyz1_rA(:, iaxis, iatom, igrid), &
                               self%lambda1_rA(iaxis, iatom, igrid), res)

               ! Normal response, buffered for the volume gradient and
               ! persisted only when requested
               dn_dR_buf(:, iatom, iaxis) = res%dn_surf
               if (allocated(self%normal1_rA)) then
                  self%normal1_rA(:, iatom, iaxis, igrid) = res%dn_surf
               end if

               self%cpjac_scal1_rA(iaxis, iatom, igrid) = res%dJ
               self%w_f1_rA(iaxis, iatom, igrid) = res%dw_f
               ! Without the branch-weight response, added by the post-pass
               self%wleb1_rA(iaxis, iatom, igrid) = res%dwleb
               self%xi1_rA(iaxis, iatom, igrid) = res%dxi

               ! Principal-curvature gradients; the mean/Gaussian curvature
               ! gradients are intentionally not stored, being downstream
               !   dKM = (k1_rA + k2_rA)/2,  dKG = k2*k1_rA + k1*k2_rA
               if (want_curvature) then
                  self%k1_rA(iaxis, iatom, igrid) = res%dk1
                  self%k2_rA(iaxis, iatom, igrid) = res%dk2
               end if
            end do ! iaxis
         end do ! i (active atoms)

         if (do_timing) call self%ctx%timer%stop(h_seed)

         !* ---------------------- f_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_sw)

         !> iswig switching derivatives: evaluated at anchor position, at the
         !> anchor width, which has no nuclear dependence
         self%f1_rA(:, :, igrid) = 0.0_wp
         call self%iswig%add_f1_rA(anchor, owner_idx, self%anchor%xi0(self%anchor_of_row(igrid)), &
                                   1.0_wp, iswig_work, self%f1_rA(:, :, igrid))

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

         ! New volume formula: v_i = (1/3) * a_i * ((r_i - c) . n_i)
         ! where n_i = dS / | dS| is the surface normal and c the volume origin

         ! Screening: a_i1_rA, dn_dR_buf, xyz1_rA all zero for inactive atoms;
         ! the origin's own motion reaches every atom and is added below
         do i = 1, n_active
            iatom = active_idx(i)
            do iaxis = 1, 3
               ! Derivative of (r_i - c) . n_i at fixed c:
               ! dr_i/dr_A . n_i + (r_i - c) . dn_i/dr_A
               ! Read normal derivative from thread-local buffer
               grad_r_hat_dot_r(iaxis) = &
                  +dot_product(self%xyz1_rA(:, iaxis, iatom, igrid), self%normal0(:, igrid)) &
                  + dot_product(point - origin, dn_dR_buf(:, iatom, iaxis))

               ! Volume derivative: dv_i/dr_A = (1/3) * [da_i/dr_A * ((r_i-c).n_i) + a_i * d((r_i-c).n_i)/dr_A]
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

         !* -------- Volume-origin term -------- *!
         ! dv_i/dc = -a_i n_i/3 reaches every atom, LSF-active or not; the sphere
         ! and total sums are spread over the atoms once after the loop
         dv_dc = -self%a(igrid)*self%normal0(:, igrid)/3.0_wp
         call add_volume_origin_gradient(dv_dc, self%v1_rA(:, :, igrid))
         dv_dc_sph_local(:, owner_idx) = dv_dc_sph_local(:, owner_idx) + dv_dc

      end do ! igrid
      !$omp end do

      !$omp critical (gradient_reduction)
      self%A_tot1_rA = self%A_tot1_rA + A_tot_local
      self%V_tot1_rA = self%V_tot1_rA + V_tot_local
      dv_dc_sph = dv_dc_sph + dv_dc_sph_local
      !$omp end critical (gradient_reduction)

      deallocate (lsf3_rr_rA, active_idx, lsf_slot, dn_dR_buf)
      deallocate (lsf1_rA, lsf2_r_rA, phi2_r_rA)
      deallocate (kkt_rhs_batch, A_tot_local, V_tot_local, dv_dc_sph_local)
      !$omp end parallel

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("compute_gradient_drop", abort%status, abort%igrid, error)
         end if
         return
      end if

      ! Spread the sphere and total volume-origin terms over every atom; summed
      ! over the atoms they cancel the net force a nonzero sum_i a_i n_i leaves
      do isph = 1, self%nsph
         call add_volume_origin_gradient(dv_dc_sph(:, isph), self%vsph1_rA(:, isph, :))
      end do
      call add_volume_origin_gradient(sum(dv_dc_sph, dim=2), self%V_tot1_rA)

      !* ========================= Branch-weight post-pass ========================= *!
      call self%ctx%timer%start(h_bw)

      ! TODO: this will change with the branch weight refactor

      ! Assemble d(wbranch_m)/dr_A for every anchor group
      max_group = 0
      if (self%ngrid > 0) max_group = maxval(self%offset(2:self%anchor%ngrid + 1) - self%offset(1:self%anchor%ngrid))
      if (max_group > 1) then
         allocate (branch_phi(max_group), source=0.0_wp)
         allocate (branch_dphi(3*self%nsph, max_group), source=0.0_wp)
         allocate (branch_weights(max_group), source=0.0_wp)
         allocate (branch_dweights(3*self%nsph, max_group), source=0.0_wp)

         do ia_group = 1, self%anchor%ngrid
            ! An anchor's rows are one contiguous group
            igroup_start = self%offset(ia_group) + 1
            igroup_end = self%offset(ia_group + 1)
            group_size = igroup_end - igroup_start + 1
            if (group_size <= 1) cycle

            ! Gather phi and dphi for every branch
            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               branch_phi(m_branch) = self%phi0(im_grid)

               owner_m = self%owner(im_grid)
               pt_m = self%xyz(:, im_grid)
               anch_m = self%anchor%xyz(:, ia_group)
               ! phi = 0.5 * alpha * |r* - anch|^2, so d phi / d r = alpha * (r* - anch)
               phi1_r_m = self%param%phi_alpha*(pt_m - anch_m)

               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     ! Chain rule: phi1_r . d r*/d r_A^iaxis
                     dphi_m = dot_product(phi1_r_m, &
                                          self%xyz1_rA(:, iaxis, iatom, im_grid))
                     ! Direct: anch moves rigidly with the owner atom, so
                     ! d phi / d R_owner = -alpha * (r* - anch) at fixed r*
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
               rn_m = dot_product(self%xyz(:, im_grid) - origin, self%normal0(:, im_grid))
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
   !> @param[in]  islot      LSF active slot, or zero when the atom is inactive
   !> @param[in]  lsf1_rA    Active-indexed dS/dR_A
   !> @param[in]  lsf2_r_rA  Active-indexed d^2S/(dr dR_A)
   !> @param[in]  lsf3_rr_rA Active-indexed packed d^3S/(dr^2 dR_A)
   !> @param[out] s1_rA      dS/dR_A of this atom [3]
   !> @param[out] s2_r_rA    d^2S/(dr dR_A) of this atom [3, 3]
   !> @param[out] s3_rr_rA   Packed d^3S/(dr^2 dR_A) of this atom [6, 3]
   pure subroutine gather_lsf_partials(islot, lsf1_rA, lsf2_r_rA, lsf3_rr_rA, &
                                       s1_rA, s2_r_rA, s3_rr_rA)
      !> LSF active slot, or zero
      integer, intent(in) :: islot
      !> Active-indexed nuclear gradient
      real(wp), intent(in) :: lsf1_rA(:, :)
      !> Active-indexed mixed second derivative
      real(wp), intent(in) :: lsf2_r_rA(:, :, :)
      !> Active-indexed packed mixed third derivative
      real(wp), intent(in) :: lsf3_rr_rA(:, :, :)
      !> Gathered nuclear gradient
      real(wp), intent(out) :: s1_rA(3)
      !> Gathered mixed second derivative
      real(wp), intent(out) :: s2_r_rA(3, 3)
      !> Gathered packed mixed third derivative
      real(wp), intent(out) :: s3_rr_rA(:, :)

      if (islot > 0) then
         s1_rA = lsf1_rA(:, islot)
         s2_r_rA = lsf2_r_rA(:, :, islot)
         s3_rr_rA = lsf3_rr_rA(:, :, islot)
      else
         s1_rA = 0.0_wp
         s2_r_rA = 0.0_wp
         s3_rr_rA = 0.0_wp
      end if
   end subroutine gather_lsf_partials

end submodule moist_cavity_drop_derivatives_forward
