!> Contract DROP surface adjoints into the density response
!>
!> Shares the sensitivity kernel with [[moist_cavity_drop_derivatives_kernel]]
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_potential
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type
   use moist_cavity_drop_derivatives_seeds, only: seed_normal_channel, seed_jet_basis
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

      call self%require_built(error)
      if (allocated(error)) return
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

      !> Level-set clone and objective used to rebuild the per-point jet; one
      !> slot, because this traversal is serial
      type(drop_worker_slots_type) :: slots
      !> Per-thread state of the grid point being opened
      type(drop_point_scratch_type) :: pt
      !> Failure latch of the prologue; read and cleared by this loop itself
      type(drop_abort_latch_type) :: abort
      !> Whether the shared prologue cleared the point
      logical :: point_ok
      !> Grid index
      integer :: igrid
      !> Point-local level-set adjoints built from the 13 jet seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Folded surface adjoints and the branch objective adjoint
      type(drop_surface_weights_type) :: eff
      real(wp) :: w_xyz_local(3)

      call self%require_context(error)
      if (allocated(error)) return
      call check_surface_adjoint(self, acc, "contract_surface_lsf_weights", error)
      if (allocated(error)) return

      ! Density variations leave anchor-only switching fixed
      call prepare_surface_weights(self, acc, .false., eff)

      ! One slot rather than the context's team: nothing below is parallel, and
      ! a clone carries the level set's screened-derivative cache with it
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii, &
                      nthreads=1)
      call pt%init(self%nsph, self%iswig, want_seed_batch=.true.)
      call abort%reset()

      w_lsf0 = 0.0_wp
      w_lsf1 = 0.0_wp
      w_lsf2 = 0.0_wp

      do igrid = 1, self%ngrid
         call drop_point_prologue(self, slots, 1, igrid, eff%have_wk, &
                                  "contract_surface_lsf_weights", abort, pt, point_ok)
         if (.not. point_ok) then
            if (allocated(abort%error)) then
               call move_alloc(abort%error, error)
               call pt%destroy()
               return
            end if
            call abort%reset()
            cycle
         end if

         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(pt%state, eff, igrid, pt%state%lsf2_rr, w_lsf1_pt, w_xyz_local)

         call seed_jet_basis(pt%state, eff, igrid, pt%phi1_r, pt%kkt_rhs, w_xyz_local, &
                             w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

         w_lsf0(igrid) = w_lsf0(igrid) + w_lsf0_pt
         w_lsf1(:, igrid) = w_lsf1(:, igrid) + w_lsf1_pt
         w_lsf2(:, :, igrid) = w_lsf2(:, :, igrid) + w_lsf2_pt
      end do

      call pt%destroy()

   end subroutine contract_surface_lsf_weights

end submodule moist_cavity_drop_derivatives_potential
