!> Contract DROP surface adjoints into the density response
!>
!> Shares the sensitivity kernel with [[moist_cavity_drop_derivatives_kernel]]
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_potential
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_surface_weights_type, build_seed_state, seed_state_ok
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_solve, seed_normal_channel, &
      & seed_jet_basis
   implicit none(type, external)

contains

   !> Accumulate surface adjoints into the density response using dS/drho
   !>
   !> Converts level-set adjoints to density adjoints and adds them to the response;
   !> density-independent cavities contribute nothing
   !>
   !> @param[in,out] self     Requires an initialized level-set model
   !> @param[in]     acc      Surface weights on the current cavity grid
   !> @param[in,out] response Existing contributions are preserved
   !> @param[out]    error    Allocated on failure
   module subroutine get_surface_response_drop(self, acc, response, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(inout) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Response list receiving the density item
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: w0(:), w1(:, :), w2(:, :, :)
      !> Cavity density response
      type(density_response_type) :: item
      !> Density chain-rule factor dS/drho
      real(wp) :: factor

      if (.not. allocated(self%lsf_model)) then
         call fatal_error(error, "get_surface_response_drop: level-set model is not initialized")
         return
      end if
      if (.not. self%lsf_model%density_adjoint_factor(factor)) return

      allocate (w0(self%ngrid), w1(3, self%ngrid), w2(3, 3, self%ngrid))
      call self%contract_surface_lsf_weights(acc, w0, w1, w2, error)
      if (allocated(error)) return

      allocate (item%w_rho(self%ngrid), source=0.0_wp)
      allocate (item%w_grad_rho(3, self%ngrid), source=0.0_wp)
      allocate (item%w_hess_rho(3, 3, self%ngrid), source=0.0_wp)
      item%w_rho = item%w_rho + factor*w0
      item%w_grad_rho = item%w_grad_rho + factor*w1
      item%w_hess_rho = item%w_hess_rho + factor*w2
      call response_accumulate(response, item, error)

   end subroutine get_surface_response_drop

   !> Contract surface adjoints into LSF value, gradient, and Hessian weights
   !>
   !> Includes position, normal, and curvature channel; area and integration
   !> weights enter through the Gaussian width; anchor-only switching is fixed
   !>
   !> @param[in]  self   Requires a projected grid
   !> @param[in]  acc    Surface weights on the current cavity grid
   !> @param[out] w_lsf0 Shape (ngrid)
   !> @param[out] w_lsf1 Shape (3, ngrid)
   !> @param[out] w_lsf2 Shape (3, 3, ngrid)
   !> @param[out] error  Allocated on failure
   module subroutine contract_surface_lsf_weights(self, acc, w_lsf0, w_lsf1, w_lsf2, error)
      !> DROP cavity with a projected grid
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> LSF value adjoints (ngrid)
      real(wp), intent(out) :: w_lsf0(:)
      !> LSF gradient adjoints (3, ngrid)
      real(wp), intent(out) :: w_lsf1(:, :)
      !> LSF Hessian adjoints (3, 3, ngrid)
      real(wp), intent(out) :: w_lsf2(:, :, :)
      !> Error allocated on failure
      type(error_type), allocatable, intent(out) :: error

      !> Local level-set evaluator
      type(lsf_thread_slot) :: lsf_slot
      type(moist_cavity_drop_objective_phi_type) :: phi
      !> Point sensitivity state
      type(drop_seed_state_type) :: state
      !> Grid index
      integer :: igrid
      !> Kernel degeneracy status
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
      !> Local LSF adjoints from 13 jet seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Effective surface adjoints
      type(drop_surface_weights_type) :: eff
      real(wp) :: w_xyz_local(3)

      call check_surface_adjoint(self, acc, "contract_surface_lsf_weights", error)
      if (allocated(error)) return

      ! Density variations leave anchor-only switching fixed
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

         ! For usable points, fold normal adjoints into gradient and position weights
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_local)

         ! Solve value and gradient sensitivities together; Hessian seeds do not move the point
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
