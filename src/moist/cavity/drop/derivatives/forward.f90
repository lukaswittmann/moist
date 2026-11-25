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
!> - atoms of a point: the LSF-active ones and the owner carry seeds (the anchor
!>   moves with its owner whether or not the level set reaches it); the
!>   switching neighbours of the anchor carry no seed and enter the
!>   switching, area and volume gradients alone
!> - molecular orientation: the anchors rotate with its sphere, which follows every
!>   atom; per point, the response to a unit rotation about each axis
!>   (anchor motion `omega x (anchor - c)` through the KKT system, the reach gate,
!>   the switching and the objective; the sphere normal through [[apply_sphere_rotation]])
!>   is spread over every atom ([[grid_frame_type]]) inside the point loop, while the
!>   point's arrays are in cache; the sphere and total area and volume take it once from
!>   the per-sphere sums, before the branch-weight post-pass, which then sees every frame term
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_forward
!$ use omp_lib, only: omp_get_thread_num
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_packed, only: packed_idx33
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, build_seed_state, apply_seed, apply_sphere_rotation, seed_state_ok, &
      & seed_state_singular_kkt
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_solve, degenerate_point_error
   use moist_cavity_drop_gaussian, only: iswig_workspace_type, iswig_pos_gradient
   use moist_cavity_type, only: add_volume_origin_gradient
   use moist_cavity_grid_frame, only: grid_frame_type
   use moist_math_linalg, only: cross_product
   implicit none(type, external)

   !> Per-point responses to a unit frame rotation, rows of the basis array
   !>
   !> Position (3), multiplier, normal (3), CP Jacobian, gate, Lebedev weight,
   !> width, principal curvatures (2), switching, area, volume, and the direct
   !> (fixed-point) change of the branch objective
   integer, parameter :: fb_xyz = 1, fb_lambda = 4, fb_normal = 5, fb_cpjac = 8, fb_w_f = 9, &
      & fb_wleb = 10, fb_xi = 11, fb_k1 = 12, fb_k2 = 13, fb_f = 14, fb_a = 15, fb_v = 16, &
      & fb_phi = 17, fb_n = 17

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
      class(cavity_type_drop), intent(inout), target :: self

      !> Error type
      type(error_type), allocatable, intent(out) :: error

      !> Frozen field: the level set's nuclear partials are dropped, each grid
      !> point moves with its owner's anchor only
      !>
      !> Equal to the full gradient for callback LSFs, whose field has no
      !> nuclear partials
      logical, intent(in), optional :: anchor_only

      !> Local copy of the anchor-only flag
      logical :: anchor_only_loc

      !> Per-thread LSF evaluators and projection objectives
      type(drop_worker_slots_type) :: slots

      !> Loop indices
      integer :: igrid, iatom, iaxis, jaxis, kaxis, i, n_active
      integer, allocatable :: active_idx(:)
      !> Entries of `active_idx` with the switching neighbours appended, and its membership mask
      integer :: n_row
      logical, allocatable :: in_row(:)

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
      integer :: h_grad, h_prim, h_pos, h_seed, h_sw, h_area, h_vol, h_bw

      !> Grid point data
      real(wp) :: point(3), anchor(3)
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
      ! weights_grad takes dphi in (nparam, nbranch) layout with
      ! nparam = 3 * nsph, flattening (iatom, iaxis) -> (iatom - 1) * 3 + iaxis
      integer :: igroup_start, igroup_end, group_size, m_branch, im_grid, ia_group, max_group
      integer :: owner_m, k_param, nbr
      real(wp) :: pt_m(3), anch_m(3), phi1_r_m(3), dphi_m
      real(wp) :: area_fac_m, rn_m, df_branch, da_branch, dv_branch
      real(wp), allocatable :: branch_phi(:), branch_gate(:)
      real(wp), allocatable :: branch_dphi(:, :), branch_dgate(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)
      !> Reach-gate partial of a row in the anchor, and the rest of its gate
      real(wp) :: reach_da(3), gate_rest

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

      !> Whether the anchors turn with a molecular frame
      logical :: frame
      !> Responses of every point to a unit rotation about each axis (fb_n, 3, ngrid)
      real(wp), allocatable :: frame_basis(:, :, :)
      !> Orientation Jacobian
      type(grid_frame_type) :: frame_jac
      !> Frame term of a point's normal with the axis before the atom, per
      !> thread (3, 3, nsph)
      real(wp), allocatable :: normal_buf(:, :, :)
      !> Area and volume responses per rotation summed per sphere (3, nsph)
      real(wp), allocatable :: frame_sum_a(:, :), frame_sum_v(:, :)
      !> KKT columns solved per point: the atom seeds, then the frame axes
      integer :: ncol
      !> Anchor offset from its center, anchor motion of one axis, switching
      !> gradient in the anchor
      real(wp) :: d_anchor(3), da_k(3), dpos(3)
      !> Response of the sphere normal's rotation
      type(drop_seed_result_type) :: res_rot
      !> Unit rotation about one axis, and the rotation's index of a row
      real(wp) :: omega(3)

      call self%require_context(error)
      if (allocated(error)) return

      anchor_only_loc = .false.
      if (present(anchor_only)) anchor_only_loc = anchor_only

      ! Thread-local primitives; the gradient reads third spatial derivatives
      timer_ref_thread = 1
      call slots%init(self%ctx, objective=.true.)

      ! Allocate gradient arrays
      if (allocated(self%xyz1_rA)) deallocate (self%xyz1_rA)
      allocate (self%xyz1_rA(3, 3, self%nsph, self%ngrid), source=0.0_wp)
      if (allocated(self%lambda1_rA)) deallocate (self%lambda1_rA)
      allocate (self%lambda1_rA(3, self%nsph, self%ngrid), source=0.0_wp)

      ! Optional derivative arrays (gated by the property flags)
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

      frame = .not. self%plain_grid()
      if (frame) then
         allocate (frame_basis(fb_n, 3, self%ngrid), source=0.0_wp)
         call self%get_frame(frame_jac)
      end if

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
      h_seed = self%ctx%timer%resolve("Seed responses", h_grad)
      h_sw = self%ctx%timer%resolve("Switching func.", h_grad)
      h_area = self%ctx%timer%resolve("Area", h_grad)
      h_vol = self%ctx%timer%resolve("Volume", h_grad)
      h_bw = self%ctx%timer%resolve("Branch weights", h_grad)

      ! Loop over all grid points (parallelized)
      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& iatom, iaxis, jaxis, kaxis, point, anchor, owner_idx, lsf0, &
      !$omp& lsf1_rA, lsf2_r_rA, lsf3_rr_rA, phi0, phi1_r, phi2_rr, phi2_r_rA, lambda_val, &
      !$omp& H_lagrangian, kkt_info, kkt_rhs_batch, state, status, res, dlsf2_rr, &
      !$omp& r_hat_dot_r, grad_r_hat_dot_r, &
      !$omp& i, n_active, active_idx, n_row, in_row, lsf_slot, s1_rA, s2_r_rA, s3_rr_rA, &
      !$omp& dn_dR_buf, iswig_work, ai_val, vi_val, &
      !$omp& A_tot_local, V_tot_local, dv_dc, dv_dc_sph_local, lsf_error, do_timing, &
      !$omp& reach_da, gate_rest, ncol, d_anchor, da_k, dpos, res_rot, omega, normal_buf)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1
      call slots%fill(thread_slot, self%lsf_model, 3, self%param, self%mol, self%radii)
      do_timing = thread_slot == timer_ref_thread .and. self%ctx%do_profile

      ! The LSF nuclear outputs are active-indexed and caller-owned; the buffers
      ! are sized to the atom count, which bounds `active_count()` from above
      allocate (lsf3_rr_rA(6, 3, self%nsph))
      allocate (active_idx(self%nsph))
      allocate (in_row(self%nsph), source=.false.)
      allocate (lsf_slot(self%nsph))
      allocate (dn_dR_buf(3, self%nsph, 3))
      allocate (lsf1_rA(3, self%nsph))
      allocate (lsf2_r_rA(3, 3, self%nsph))
      allocate (phi2_r_rA(3, 3, self%nsph))
      call self%iswig%new_workspace(iswig_work)
      allocate (kkt_rhs_batch(4, 3*self%nsph + 3))
      allocate (A_tot_local(3, self%nsph), source=0.0_wp)
      allocate (V_tot_local(3, self%nsph), source=0.0_wp)
      allocate (dv_dc_sph_local(3, self%nsph), source=0.0_wp)
      if (frame .and. allocated(self%normal1_rA)) allocate (normal_buf(3, 3, self%nsph))

      !$omp do schedule(dynamic)
      do igrid = 1, self%ngrid
         !$omp cancellation point do
         if (abort%is_set()) cycle
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

         ! Get nuclear and mixed derivatives; a frozen field has none
         if (.not. anchor_only_loc) call slots%lsf(thread_slot)%lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
         call slots%lsf(thread_slot)%lsf%f3_rrr(lsf0, state%lsf1_r, state%lsf2_rr, state%lsf3_rrr)

         ! Per-point kernel state; a degenerate point is latched
         call fill_seed_state(self, igrid, want_curvature, state, reach_da)
         state%lambda_val = lambda_val
         call build_seed_state(state, self%f_crit, self%f_foc, self%f_hor, self%f_wleb, &
                               self%param%finalize%wleb_prune, status)
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
         ! phi2_r_rA is nonzero only at owner_idx, which is seeded whether or
         ! not the level set reaches it
         if (anchor_only_loc) then
            ! Frozen field: only the owner moves the point, and it carries no
            ! field seed; the switching neighbours are appended below
            n_active = 1
            active_idx(1) = owner_idx
            lsf_slot(1) = 0
         else
            n_active = slots%lsf(thread_slot)%lsf%active_count()
            do i = 1, n_active
               active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
               lsf_slot(i) = i
            end do
            ! An owner outside the level set's reach, or a level set without
            ! active atoms, still moves the anchor
            if (.not. any(active_idx(1:n_active) == owner_idx)) then
               n_active = n_active + 1
               active_idx(n_active) = owner_idx
               lsf_slot(n_active) = 0
            end if
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

         ! Molecular frame: the anchor turns about each axis, d(anchor) =
         ! e_k x (anchor - c) at a fixed field, the anchor seed's RHS alpha*d(anchor)
         ncol = 3*n_active
         if (frame) then
            d_anchor = anchor - self%mol%xyz(:, owner_idx)
            do kaxis = 1, 3
               omega = 0.0_wp
               omega(kaxis) = 1.0_wp
               kkt_rhs_batch(1:3, ncol + kaxis) = self%param%phi_alpha*cross_product(omega, d_anchor)
               kkt_rhs_batch(4, ncol + kaxis) = 0.0_wp
            end do
            ncol = ncol + 3
         end if

         ! Single factorization + solve for all RHS
         call drop_kkt_solve(H_lagrangian, state%lsf1_r, kkt_rhs_batch(:, 1:ncol), kkt_info)
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

         ! The reach gate moves with the anchor at a fixed projected point; the
         ! owner need not be LSF-active
         gate_rest = state%wleb_prune_f0*state%f_crit0*state%f_foc_f0*state%f_hor_f0
         self%w_f1_rA(:, owner_idx, igrid) = self%w_f1_rA(:, owner_idx, igrid) + gate_rest*reach_da

         if (do_timing) call self%ctx%timer%stop(h_seed)

         !* ---------------------- f_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_sw)

         !> iswig switching derivatives: evaluated at anchor position, at the
         !> anchor width, which has no nuclear dependence
         !> `f = f_anchor * wbranch`; the response of the branch weight, the gate
         !> included, is added by the post-pass
         self%f1_rA(:, :, igrid) = 0.0_wp
         call self%iswig%add_f1_rA(anchor, owner_idx, self%anchor%xi0(self%anchor_of_row(igrid)), &
                                   self%wbranch(igrid), iswig_work, self%f1_rA(:, :, igrid))
         ! The anchor also turns with its sphere: the switching gradient in it
         if (frame) dpos = iswig_pos_gradient(iswig_work, self%anchor%xi0(self%anchor_of_row(igrid)), &
                                              self%wbranch(igrid))

         ! The switching neighbours of the anchor need not be LSF-active at the
         ! projected point: append them, seedless, for the area and volume
         n_row = n_active
         in_row(active_idx(1:n_active)) = .true.
         do i = 1, iswig_work%n
            iatom = iswig_work%k(i)
            if (in_row(iatom)) cycle
            in_row(iatom) = .true.
            n_row = n_row + 1
            active_idx(n_row) = iatom
         end do
         in_row(active_idx(1:n_row)) = .false.

         if (do_timing) call self%ctx%timer%stop(h_sw)

         !* ---------------------- a_i derivative ---------------------- *!
         if (do_timing) call self%ctx%timer%start(h_area)
         ! Screening: f1_rA is zero off the seeded atoms and the switching
         ! neighbours, wleb1_rA off the seeded atoms
         do i = 1, n_row
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

         ! Screening: dn_dR_buf and xyz1_rA are zero off the seeded atoms, a_i1_rA
         ! off them and the switching neighbours;
         ! the origin's own motion reaches every atom and is added below
         do i = 1, n_row
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
         do i = 1, n_row
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

         !* -------- Frame rotation, per axis -------- *!
         ! The anchor moves by da = e_k x (anchor - c) (point motion from the KKT
         ! columns after the atom seeds, reach gate, switching, objective) and the
         ! sphere normal turns with it
         if (frame) then
            do kaxis = 1, 3
               omega = 0.0_wp
               omega(kaxis) = 1.0_wp
               da_k = cross_product(omega, d_anchor)
               dlsf2_rr = 0.0_wp
               call apply_seed(state, [0.0_wp, 0.0_wp, 0.0_wp], dlsf2_rr, &
                               kkt_rhs_batch(1:3, 3*n_active + kaxis), kkt_rhs_batch(4, 3*n_active + kaxis), res)
               call apply_sphere_rotation(state, omega, res_rot)
               associate (fb => frame_basis(:, kaxis, igrid))
                  fb(fb_xyz:fb_xyz + 2) = kkt_rhs_batch(1:3, 3*n_active + kaxis)
                  fb(fb_lambda) = kkt_rhs_batch(4, 3*n_active + kaxis)
                  fb(fb_normal:fb_normal + 2) = res%dn_surf
                  fb(fb_cpjac) = res%dJ + res_rot%dJ
                  fb(fb_w_f) = res%dw_f + res_rot%dw_f + gate_rest*dot_product(reach_da, da_k)
                  fb(fb_wleb) = res%dwleb + res_rot%dwleb
                  fb(fb_xi) = res%dxi + res_rot%dxi
                  fb(fb_k1) = res%dk1
                  fb(fb_k2) = res%dk2
                  fb(fb_f) = dot_product(dpos, da_k)
                  fb(fb_a) = self%radii(owner_idx)**2*(fb(fb_f)*self%wleb(igrid) + self%f(igrid)*fb(fb_wleb))
                  ! d((r_i - c) . n_i) at fixed c, as in the volume gradient above
                  grad_r_hat_dot_r(kaxis) = dot_product(fb(fb_xyz:fb_xyz + 2), self%normal0(:, igrid)) &
                                            + dot_product(point - origin, res%dn_surf)
                  fb(fb_v) = (fb(fb_a)*r_hat_dot_r + self%a(igrid)*grad_r_hat_dot_r(kaxis))/3.0_wp
                  fb(fb_phi) = -dot_product(phi1_r, da_k)
               end associate
            end do

            ! Spread over every atom now, while the point's arrays are in cache
            associate (fb => frame_basis(:, :, igrid))
               call frame_jac%add_vector(self%comp(igrid), transpose(fb(fb_xyz:fb_xyz + 2, :)), &
                  & self%xyz1_rA(:, :, :, igrid))
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_lambda, :), self%lambda1_rA(:, :, igrid))
               if (allocated(self%normal1_rA)) then
                  normal_buf = 0.0_wp
                  call frame_jac%add_vector(self%comp(igrid), transpose(fb(fb_normal:fb_normal + 2, :)), normal_buf)
                  do iaxis = 1, 3
                     self%normal1_rA(:, :, iaxis, igrid) = self%normal1_rA(:, :, iaxis, igrid) &
                                                           + normal_buf(:, iaxis, :)
                  end do
               end if
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_cpjac, :), self%cpjac_scal1_rA(:, :, igrid))
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_w_f, :), self%w_f1_rA(:, :, igrid))
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_wleb, :), self%wleb1_rA(:, :, igrid))
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_xi, :), self%xi1_rA(:, :, igrid))
               if (want_curvature) then
                  call frame_jac%add_scalar(self%comp(igrid), fb(fb_k1, :), self%k1_rA(:, :, igrid))
                  call frame_jac%add_scalar(self%comp(igrid), fb(fb_k2, :), self%k2_rA(:, :, igrid))
               end if
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_f, :), self%f1_rA(:, :, igrid))
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_a, :), self%a_i1_rA(:, :, igrid))
               call frame_jac%add_scalar(self%comp(igrid), fb(fb_v, :), self%v1_rA(:, :, igrid))
            end associate
         end if

      end do ! igrid
      !$omp end do

      !$omp critical (gradient_reduction)
      self%A_tot1_rA = self%A_tot1_rA + A_tot_local
      self%V_tot1_rA = self%V_tot1_rA + V_tot_local
      dv_dc_sph = dv_dc_sph + dv_dc_sph_local
      !$omp end critical (gradient_reduction)

      deallocate (lsf3_rr_rA, active_idx, in_row, lsf_slot, dn_dR_buf)
      deallocate (lsf1_rA, lsf2_r_rA, phi2_r_rA)
      deallocate (kkt_rhs_batch, A_tot_local, V_tot_local, dv_dc_sph_local)
      if (allocated(normal_buf)) deallocate (normal_buf)
      call slots%release(thread_slot)
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

      ! Molecular frame: the sphere and total area and volume from the
      ! per-sphere sums, before the post-pass reads them
      if (frame) then
         allocate (frame_sum_a(3, self%grid%ncomp), frame_sum_v(3, self%grid%ncomp), source=0.0_wp)
         do igrid = 1, self%ngrid
            frame_sum_a(:, self%comp(igrid)) = frame_sum_a(:, self%comp(igrid)) + frame_basis(fb_a, :, igrid)
            frame_sum_v(:, self%comp(igrid)) = frame_sum_v(:, self%comp(igrid)) + frame_basis(fb_v, :, igrid)
         end do
         call frame_jac%add_spheres(frame_sum_a, self%asph1_rA)
         call frame_jac%add_spheres(frame_sum_v, self%vsph1_rA)
         call frame_jac%add_total(frame_sum_a, self%A_tot1_rA)
         call frame_jac%add_total(frame_sum_v, self%V_tot1_rA)
      end if

      !* ========================= Branch-weight post-pass ========================= *!
      call self%ctx%timer%start(h_bw)

      ! The branch weight of a row depends on the objectives and gates of
      ! every row of its anchor; a single row has `wbranch = w_f`

      ! Assemble d(wbranch_m)/dr_A for every anchor group
      max_group = 0
      if (self%ngrid > 0) then
         max_group = maxval(self%offset(2:self%anchor%ngrid + 1) - self%offset(1:self%anchor%ngrid))
      end if
      if (max_group > 0) then
         allocate (branch_phi(max_group), branch_gate(max_group), source=0.0_wp)
         allocate (branch_dphi(3*self%nsph, max_group), source=0.0_wp)
         allocate (branch_dgate(3*self%nsph, max_group), source=0.0_wp)
         allocate (branch_weights(max_group), source=0.0_wp)
         allocate (branch_dweights(3*self%nsph, max_group), source=0.0_wp)

         do ia_group = 1, self%anchor%ngrid
            ! An anchor's rows are one contiguous group
            igroup_start = self%offset(ia_group) + 1
            igroup_end = self%offset(ia_group + 1)
            group_size = igroup_end - igroup_start + 1
            if (group_size <= 0) cycle
            nbr = group_size

            ! Gather phi, gate and their derivatives for every row
            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               branch_phi(m_branch) = self%phi0(im_grid)
               branch_gate(m_branch) = self%w_f0(im_grid)
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
                     ! and turns with the owner's frame
                     if (frame) then
                        dphi_m = dphi_m + dot_product(frame_jac%gk(:, k_param, frame_jac%rot(self%comp(im_grid))), &
                                                      frame_basis(fb_phi, :, im_grid))
                     end if
                     branch_dphi(k_param, m_branch) = dphi_m
                     branch_dgate(k_param, m_branch) = self%w_f1_rA(iaxis, iatom, im_grid)
                  end do
               end do
            end do
            ! Branch weights and their derivatives over the anchor's group
            call self%branch_weight%weights_grad( &
               branch_phi(1:nbr), branch_gate(1:nbr), &
               branch_dphi(:, 1:nbr), branch_dgate(:, 1:nbr), &
               weights=branch_weights(1:nbr), &
               dweights=branch_dweights(:, 1:nbr))

            ! f = f_anchor * wbranch: add f_anchor * d(wbranch)/dr_A to the row's
            ! switching gradient and carry it into the area and volume gradients
            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               owner_m = self%owner(im_grid)
               area_fac_m = self%radii(owner_m)**2*self%wleb(im_grid)
               rn_m = dot_product(self%xyz(:, im_grid) - origin, self%normal0(:, im_grid))
               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     k_param = (iatom - 1)*3 + iaxis
                     df_branch = self%anchor%f0(ia_group)*branch_dweights(k_param, m_branch)
                     if (df_branch == 0.0_wp) cycle
                     self%f1_rA(iaxis, iatom, im_grid) = &
                        self%f1_rA(iaxis, iatom, im_grid) + df_branch

                     da_branch = area_fac_m*df_branch
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

         deallocate (branch_phi, branch_gate, branch_dphi, branch_dgate, &
                     branch_weights, branch_dweights)
      end if

      call self%ctx%timer%stop(h_bw)

   end subroutine compute_gradient_drop

   !> Compute anchor-only nuclear derivatives for the DROP cavity
   !>
   !> - frozen field: the level set's nuclear partials are dropped, each grid
   !>   point moves with its owner's anchor, the switching with its neighbours
   !> - equal to the full gradient for callback LSFs
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
