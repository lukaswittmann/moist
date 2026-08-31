!> Reverse-mode surface -> level set adjoint contractions for the DROP cavity.
!>
!> Provides the variational cavity response (Fock) infrastucture
!>
!> These routines map per-point surface adjoint weights (Gaussian width,
!> integration weight, area, switch factor, projected position, and normal) onto
!> adjoint weights of the level set function value/gradient/Hessian
!>
!> The per-grid point sensitivity kernel is shared with the nuclear path
!> in [[moist_cavity_drop_derivatives_kernel]]
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_potential
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_surface_weights_type, build_seed_state, seed_state_ok
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_solve, seed_normal_channel, &
      & seed_jet_basis
   implicit none (type, external)

contains

   !> Map accumulated surface adjoints into the generic response container
   !>
   !> @param[inout] self     DROP cavity instance
   !> @param[in]    acc      Accumulated surface-observable adjoints
   !> @param[inout] response Response accumulator receiving the LSF channels
   !> @param[out]   error    Error object
   module subroutine get_surface_response_drop(self, acc, response, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Response accumulator receiving the LSF channels
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)

      allocate (w0(self%ngrid), w1(3, self%ngrid), w2(3, 3, self%ngrid))
      call self%contract_surface_lsf_weights(acc, w0, w1, w2, error)
      if (allocated(error)) return

      if (.not. allocated(response%lsf%w_value)) then
         allocate (response%lsf%w_value(self%ngrid), source=0.0_wp)
         allocate (response%lsf%w_gradient(3, self%ngrid), source=0.0_wp)
         allocate (response%lsf%w_hessian(3, 3, self%ngrid), source=0.0_wp)
      else if (size(response%lsf%w_value) /= self%ngrid) then
         call fatal_error(error, "DROP surface response grid-size mismatch")
         return
      end if
      response%lsf%w_value = response%lsf%w_value + w0
      response%lsf%w_gradient = response%lsf%w_gradient + w1
      response%lsf%w_hessian = response%lsf%w_hessian + w2

   end subroutine get_surface_response_drop

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
      !> Shared per-grid point sensitivity kernel state and its response
      type(drop_seed_state_type) :: state
      !> Grid, seed and Cartesian indices
      integer :: igrid
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
      real(wp) :: kkt_rhs(4, 4)
      integer(lapack_ik) :: kkt_info
      !> Point-local level-set adjoints built from the 13 jet seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Folded surface adjoints and the branch objective adjoint
      type(drop_surface_weights_type) :: eff
      real(wp) :: w_xyz_local(3)

      call check_surface_adjoint(self, acc, "contract_surface_lsf_weights", error)
      if (allocated(error)) return

      ! fold_switching = .false.: the electronic degrees of freedom leave the
      ! switching factor f untouched, so the area channel's da/df term is
      ! identically zero here. The nuclear path passes .true.
      call prepare_surface_weights(self, acc, .false., eff)

      allocate (lsf_slot%lsf, source=self%lsf_model)
      call lsf_slot%lsf%set_max_deriv(3)
      call phi%set_parameters(self%param)
      call phi%set_input(self%mol, self%radii)
      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)

      w_lsf0 = 0.0_wp
      w_lsf1 = 0.0_wp
      w_lsf2 = 0.0_wp

      do igrid = 1, self%ngrid
         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         ! Serial loop, so an evaluation failure can be returned immediately
         ! (see the parallel loops for the general contract).
         call lsf_slot%lsf%prepare(point, error)
         if (allocated(error)) return
         call lsf_slot%lsf%f3_rrr(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         call phi%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         state%lsf1_r = lsf1_r
         state%lsf2_rr = lsf2_rr
         state%lsf3_rrr = lsf3_rrr
         state%lambda_val = lambda_val
         call fill_seed_state(self, igrid, eff%have_wk, state)

         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)

         if (status /= seed_state_ok) cycle

         ! Fold an optional outward-normal adjoint weight into the field channels:
         ! the direct grad-S contribution normal_grad = P_tan(w_n)/|grad S| enters
         ! w_lsf1 at the fixed projected point, and its point-motion coupling
         ! H @ normal_grad augments the effective position weight below.
         !
         ! This write happens only once the point is known to be usable, so a
         ! rejected point never leaves a half-contracted weight behind.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_local)

         ! KKT sensitivities for all 13 basis perturbations from one
         ! factorization: only the value (ibasis 1) and gradient (ibasis 2-4)
         ! perturbations enter the right-hand side; the nine Hessian
         ! perturbations have rhs = 0 and hence dr/dp = 0, dlambda/dp = 0.
         ! Only the value (column 1) and gradient (columns 2-4) seeds move the
         ! point; the nine Hessian seeds have a zero right-hand side.
         kkt_rhs = 0.0_wp
         kkt_rhs(4, 1) = -1.0_wp
         kkt_rhs(1, 2) = lambda_val
         kkt_rhs(2, 3) = lambda_val
         kkt_rhs(3, 4) = lambda_val
         call drop_kkt_solve(phi2_rr - lambda_val*lsf2_rr, lsf1_r, kkt_rhs, kkt_info)
         if (kkt_info /= 0_lapack_ik) then
            call fatal_error(error, &
                             "contract_surface_lsf_weights: KKT sensitivity solve failed")
            return
         end if

         call seed_jet_basis(state, eff, igrid, phi1_r, kkt_rhs, w_xyz_local, &
                             w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

         w_lsf0(igrid) = w_lsf0(igrid) + w_lsf0_pt
         w_lsf1(:, igrid) = w_lsf1(:, igrid) + w_lsf1_pt
         w_lsf2(:, :, igrid) = w_lsf2(:, :, igrid) + w_lsf2_pt
      end do

   end subroutine contract_surface_lsf_weights

end submodule moist_cavity_drop_derivatives_potential
