!> Reverse-mode surface -> level set adjoint contractions for the DROP cavity.
!>
!> These routines map per-tessera surface adjoint weights (Gaussian width,
!> integration weight, area, switch factor, projected position, and normal) onto
!> adjoint weights of the level set function value/gradient/Hessian.
!>
!> They provide the variational solvation-potential (Fock) response; the
!> nuclear-derivative contractions of the same quantities live in
!> nuclear_gradient.f90 (reverse mode) and gradient.f90 (forward mode).
!>
!> The per-grid-point sensitivity kernel is shared with the nuclear path and
!> lives in [[moist_cavity_drop_adjoint_kernel]].
submodule(moist_cavity_drop) moist_cavity_drop_adjoint
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_adjoint_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, build_seed_state, apply_seed, compute_branch_phi_adj, &
      & seed_state_ok, seed_weight_tol, seed_status_message
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
   !> The switching factor $f_i$ is an anchor-only iSwig overlap, so $\partial
   !> f_i/\partial p = 0$ for a level-set perturbation and `acc%w_f` does not
   !> contribute here. It does contribute to the nuclear gradient, which is why
   !> the area channel is folded into the width channel only.
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

      !> Level-set clone and objective used to rebuild the per-point jet
      type(lsf_thread_slot) :: lsf_slot
      type(moist_cavity_drop_objective_phi_type) :: phi
      !> Shared per-grid-point sensitivity kernel state and its response
      type(drop_seed_state_type) :: state
      type(drop_seed_result_type) :: res
      !> Grid, seed and Cartesian indices
      integer :: igrid, ibasis, iaxis, jaxis
      !> Degeneracy status returned by the kernel
      integer :: status
      !> Projected point, anchor and owner sphere
      real(wp) :: point(3), anchor(3)
      integer :: owner_idx
      !> Level-set jet at the projected point
      real(wp) :: lsf0, lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      !> Objective jet at the projected point
      real(wp) :: phi0, phi1_r(3), phi2_rr(3, 3)
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val
      !> Bordered KKT sensitivity system
      real(wp) :: kkt_mat_base(4, 4), kkt_mat(4, 4), kkt_rhs(4, 4)
      integer(lapack_ik) :: kkt_ipiv(4), kkt_info
      !> Seed perturbation of the level-set jet and the induced point motion
      real(wp) :: dlsf1_r(3), dlsf2_rr(3, 3)
      real(wp) :: dr_dp(3), dlambda_dp
      !> Accumulated adjoint contribution of one seed
      real(wp) :: contribution
      !> Softmax temperature and the branch objective adjoint
      real(wp) :: sigma_phi
      real(wp), allocatable :: branch_phi_adj(:)

      !> Effective primitive surface adjoints used by the DROP contraction
      real(wp), allocatable :: w_xi(:), w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp) :: normal_grad(3), w_xyz_local(3), nwn
      logical :: have_wn, have_wk

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
      if (any(abs(self%xi0) <= seed_weight_tol .and. &
              (abs(acc%w_a) > seed_weight_tol .or. abs(acc%w_w) > seed_weight_tol))) then
         call fatal_error(error, "contract_surface_lsf_weights: singular derived-weight conversion")
         return
      end if

      allocate (w_xi, source=acc%w_xi)
      allocate (w_xyz, source=acc%w_xyz)
      allocate (w_n, source=acc%w_n)
      allocate (w_k1, source=acc%w_k1)
      allocate (w_k2, source=acc%w_k2)
      where (abs(acc%w_a) > seed_weight_tol)
         w_xi = w_xi - 2.0_wp*self%a*acc%w_a/self%xi0
      end where
      where (abs(acc%w_w) > seed_weight_tol)
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

      have_wn = any(abs(w_n) > seed_weight_tol)
      have_wk = any(abs(w_k1) > seed_weight_tol) .or. any(abs(w_k2) > seed_weight_tol)

      ! Reverse pass for the branch-weight post-pass: converts the width-induced
      ! adjoint dL/dp_m into dL/dPhi_m, which the seed loop couples to dr/dp
      if (self%ngrid > 0 .and. allocated(self%branch_count)) then
         sigma_phi = self%branch_weight%s
         call compute_branch_phi_adj(self%branch_count(1:self%ngrid), &
                                     self%anchor_id(1:self%ngrid), &
                                     self%wbranch(1:self%ngrid), &
                                     self%wleb(1:self%ngrid), &
                                     self%xi0(1:self%ngrid), &
                                     sigma_phi, w_xi, branch_phi_adj)
      end if

      do igrid = 1, self%ngrid
         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         ! Serial loop, so an evaluation failure can be returned immediately
         ! (see the parallel loops for the general contract).
         call lsf_slot%lsf%prepare(point, error)
         if (allocated(error)) return
         call lsf_slot%lsf%f3_rrr_screened(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         call phi%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         kkt_mat_base = 0.0_wp
         kkt_mat_base(1:3, 1:3) = phi2_rr - lambda_val*lsf2_rr
         kkt_mat_base(1:3, 4) = -lsf1_r
         kkt_mat_base(4, 1:3) = lsf1_r

         state%lsf1_r = lsf1_r
         state%lsf2_rr = lsf2_rr
         state%lsf3_rrr = lsf3_rrr
         state%lambda_val = lambda_val
         state%alpha_coeff = self%param%phi_alpha
         state%anchor = anchor
         state%owner_xyz = self%mol%xyz(:, owner_idx)
         state%anchor_wleb0 = self%anchor_wleb0(igrid)
         state%cpjac_scal0 = self%cpjac_scal0(igrid)
         state%w_f0 = self%w_f0(igrid)
         state%wbranch = self%wbranch(igrid)
         state%wleb = self%wleb(igrid)
         state%xi0 = self%xi0(igrid)
         state%want_curvature = have_wk

         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)
         if (status /= seed_state_ok) then
            call degenerate_point_error(status, igrid, error)
            return
         end if

         ! Fold an optional outward-normal adjoint weight into the field channels:
         ! the direct grad-S contribution normal_grad = P_tan(w_n)/|grad S| enters
         ! w_lsf1 at the fixed projected point, and its point-motion coupling
         ! H @ normal_grad augments the effective position weight below.
         !
         ! This write happens only once the point is known to be usable, so a
         ! rejected point never leaves a half-contracted weight behind.
         if (have_wn) then
            nwn = dot_product(state%n_surf, w_n(:, igrid))
            normal_grad = (w_n(:, igrid) - state%n_surf*nwn)/state%g_norm
            w_lsf1(:, igrid) = w_lsf1(:, igrid) + normal_grad
            w_xyz_local = w_xyz(:, igrid) + matmul(lsf2_rr, normal_grad)
         else
            w_xyz_local = w_xyz(:, igrid)
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
            dlsf1_r = 0.0_wp
            dlsf2_rr = 0.0_wp
            if (ibasis <= 4) then
               if (ibasis > 1) dlsf1_r(ibasis - 1) = 1.0_wp
               dr_dp = kkt_rhs(1:3, ibasis)
               dlambda_dp = kkt_rhs(4, ibasis)
            else
               iaxis = (ibasis - 5)/3 + 1
               jaxis = mod(ibasis - 5, 3) + 1
               dlsf2_rr(iaxis, jaxis) = 1.0_wp
               dr_dp = 0.0_wp
               dlambda_dp = 0.0_wp
            end if

            call apply_seed(state, dlsf1_r, dlsf2_rr, dr_dp, dlambda_dp, res)

            ! self%f is an anchor-only iSwig overlap for electronic LSF
            ! perturbations, so df/dp = 0 and w_f does not contribute here
            contribution = dot_product(w_xyz_local, dr_dp) + w_xi(igrid)*res%dxi
            if (abs(branch_phi_adj(igrid)) > seed_weight_tol) then
               contribution = contribution + branch_phi_adj(igrid)*dot_product(phi1_r, dr_dp)
            end if
            if (have_wk) then
               contribution = contribution + w_k1(igrid)*res%dk1 + w_k2(igrid)*res%dk2
            end if

            if (ibasis == 1) then
               w_lsf0(igrid) = w_lsf0(igrid) + contribution
            else if (ibasis <= 4) then
               w_lsf1(ibasis - 1, igrid) = w_lsf1(ibasis - 1, igrid) + contribution
            else
               w_lsf2(iaxis, jaxis, igrid) = w_lsf2(iaxis, jaxis, igrid) + contribution
            end if
         end do
      end do

   end subroutine contract_surface_lsf_weights

   !> Report a degenerate projected point from the sensitivity kernel
   !>
   !> The forward gradient treats these geometries as fatal; the reverse paths
   !> match that policy so a Fock response or a nuclear gradient is never
   !> silently zeroed on a subset of the grid.
   !>
   !> @param[in]  status  One of the `seed_state_*` codes
   !> @param[in]  igrid   Offending grid point
   !> @param[out] error   Error object
   subroutine degenerate_point_error(status, igrid, error)
      !> Degeneracy status
      integer, intent(in) :: status
      !> Offending grid point
      integer, intent(in) :: igrid
      !> Error object
      type(error_type), allocatable, intent(out) :: error

      !> Rendered grid-point index
      character(len=32) :: idx

      write (idx, '(i0)') igrid
      call fatal_error(error, "[Error] contract_surface_lsf_weights: "// &
                       seed_status_message(status)//" at grid point "//trim(idx))

   end subroutine degenerate_point_error

end submodule moist_cavity_drop_adjoint
