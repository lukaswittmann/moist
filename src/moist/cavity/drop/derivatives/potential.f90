!> Contract DROP surface adjoints into the density response
!>
!> Shares the sensitivity kernel with [[moist_cavity_drop_derivatives_kernel]]
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_potential
   !$ use omp_lib, only: omp_get_thread_num
   use moist_math_packed, only: unpack_sym_matrix
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_surface_weights_type, seed_state_ok
   use moist_cavity_drop_derivatives_seeds, only: seed_point_jet, degenerate_point_error
   use moist_cavity_drop_tile, only: surface_tile_width, surface_tile_candidates, &
      & surface_jet_tile
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
      type(cavity_adjoint_type), intent(in) :: acc
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
      type(cavity_adjoint_type), intent(in) :: acc
      !> LSF value adjoints (ngrid)
      real(wp), intent(out) :: w_lsf0(:)
      !> LSF gradient adjoints (3, ngrid)
      real(wp), intent(out) :: w_lsf1(:, :)
      !> LSF Hessian adjoints (3, 3, ngrid)
      real(wp), intent(out) :: w_lsf2(:, :, :)
      !> Error allocated on failure
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Thread slot
      integer :: thread_slot
      !> Point sensitivity state
      type(drop_seed_state_type) :: state
      !> Grid and tile-column index
      integer :: igrid, p
      !> Tile index, tile count, tile width, first grid point and size
      integer :: itile, ntile, tile, igrid0, np
      !> Shared candidate list of the tile
      integer, allocatable :: cand(:)
      integer :: nc
      !> Kernel degeneracy status
      integer :: status
      !> Owner sphere
      integer :: owner_idx
      !> Tile points and level-set jet, packed Hessian and third derivative
      real(wp), allocatable :: tile_xyz(:, :), tile_lsf0(:), tile_lsf1(:, :)
      real(wp), allocatable :: tile_lsf2(:, :), tile_lsf3(:, :)
      !> Objective gradient at the projected point
      real(wp) :: phi1_r(3)
      !> Bordered KKT sensitivity system
      real(wp) :: kkt_rhs(4, 4)
      !> Local LSF adjoints from 10 jet seeds, Hessian packed
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(6)
      !> Effective surface adjoints
      type(drop_surface_weights_type) :: eff
      real(wp) :: w_xyz_local(3)
      !> Timer handle
      integer :: h_weights

      call self%require_context(error)
      if (allocated(error)) return

      call check_surface_adjoint(self, acc, "contract_surface_lsf_weights", error)
      if (allocated(error)) return
      if (size(w_lsf0) /= self%ngrid .or. any(shape(w_lsf1) /= [3, self%ngrid]) &
          .or. any(shape(w_lsf2) /= [3, 3, self%ngrid])) then
         call fatal_error(error, "contract_surface_lsf_weights: weight shape mismatch")
         return
      end if

      w_lsf0 = 0.0_wp
      w_lsf1 = 0.0_wp
      w_lsf2 = 0.0_wp
      if (self%ngrid <= 0) return

      h_weights = self%ctx%timer%resolve("Surface LSF weights", self%ctx%timer%current(), &
                                         cat_gradient)
      call self%ctx%timer%start(h_weights)

      ! Density variations leave anchor-only switching fixed
      call prepare_surface_weights(self, acc, .false., eff)

      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      tile = surface_tile_width
      ntile = (self%ngrid + tile - 1)/tile

      call abort%reset()

      ! Tiles of `tile` points around one batched jet; outputs are per-point columns;
      ! threads share only the latch
      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, &
      !$omp& itile, igrid0, np, igrid, p, cand, nc, state, status, owner_idx, &
      !$omp& tile_xyz, tile_lsf0, tile_lsf1, tile_lsf2, tile_lsf3, &
      !$omp& phi1_r, kkt_rhs, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_local)
      thread_slot = 1
      !$ thread_slot = omp_get_thread_num() + 1

      allocate (cand(max(1, slots%lsf(thread_slot)%lsf%ncenters)))
      allocate (tile_xyz(3, tile), tile_lsf0(tile), tile_lsf1(3, tile), tile_lsf2(6, tile))
      allocate (tile_lsf3(10, tile))

      !$omp do schedule(dynamic, 1)
      do itile = 1, ntile
         if (abort%requested) cycle

         igrid0 = (itile - 1)*tile + 1
         np = min(tile, self%ngrid - igrid0 + 1)
         tile_xyz(:, 1:np) = self%xyz(:, igrid0:igrid0 + np - 1)

         call surface_tile_candidates(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand, nc)
         if (.not. surface_jet_tile(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand(1:nc), &
                                    igrid0, "contract_surface_lsf_weights", 3, tile_lsf0, tile_lsf1, &
                                    tile_lsf2, tile_lsf3, abort)) cycle

         do p = 1, np
            if (abort%requested) cycle

            igrid = igrid0 + p - 1
            owner_idx = self%owner(igrid)
            call fill_seed_state(self, igrid, eff%have_wk, state)

            ! A degenerate point or a singular KKT system voids the result
            call seed_point_jet(slots%phi(thread_slot), tile_xyz(:, p), owner_idx, tile_lsf1(:, p), &
                                tile_lsf2(:, p), tile_lsf3(:, p), self%lambda0(igrid), self%f_crit, &
                                self%f_foc, self%f_wleb, self%param%wleb_prune_level > 0, eff, &
                                igrid, state, kkt_rhs, phi1_r, w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, &
                                w_xyz_local, status)
            if (status /= seed_state_ok) then
               call abort%latch_status(status, igrid)
               cycle
            end if

            w_lsf0(igrid) = w_lsf0_pt
            w_lsf1(:, igrid) = w_lsf1_pt
            call unpack_sym_matrix(w_lsf2_pt, w_lsf2(:, :, igrid), 3)
         end do
      end do
      !$omp end do

      deallocate (cand, tile_xyz, tile_lsf0, tile_lsf1, tile_lsf2, tile_lsf3)
      !$omp end parallel

      call self%ctx%timer%stop(h_weights)

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("contract_surface_lsf_weights", abort%status, &
                                        abort%igrid, error)
         end if
      end if

   end subroutine contract_surface_lsf_weights

end submodule moist_cavity_drop_derivatives_potential
