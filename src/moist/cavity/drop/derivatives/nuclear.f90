!> Reverse-mode (z-vector) nuclear gradient for the DROP cavity
!>
!> Contracts an accumulated surface adjoint directly into `dE/dR_A`;
!> The per-point cost is independent of the system size, without persistent
!> cavity derivative arrays; the scratch is per tile (compact gradient rows)
!> and per thread (one `(3, nsph)` tile gradient, of which a tile touches and
!> clears only its own rows)
!>
!> Legacy forward path in forward.f90 builds the full forward Jacobian
!> of every surface observable: `3 * n_active` seeds per grid point
!> and `O(N_sph * N_grid)` memory
!>
!> A nuclear displacement reaches the surface through three routes:
!>
!>  1. Field: nuclear dependence of the level set
!>     - the per-point map is linear in its seed, so seeding the 10 components
!>       of the level-set jet gives adjoint weights `w_lsf0/1/2` contracting
!>       directly with the LSF's own nuclear partials
!>
!>  2. Anchor: anchor-owner sphere dependence, `d(anch)/dR_A = delta_{A,own}`
!>     - adds three seeds via the objective's mixed derivative of
!>       `-d^2 phi/dr dR = alpha*I`, affecting only the owner atom
!>
!>  3. Switching: the iSwiG `f_i` is a function of the nuclear coordinates and radii
!>     - contracted with a scalar weight
!>     - `f_i` belongs to the anchor, so the rows of one anchor share one sparse
!>       evaluation over the owner's neighbours
!>
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_nuclear
   !$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_surface_weights_type, seed_state_ok, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: seed_point_jet, seed_anchor, &
      & degenerate_point_error
   use moist_cavity_drop_tile, only: surface_tile_width, surface_tile_candidates, &
      & surface_jet_tile, surface_vjp_tile
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   implicit none(type, external)

   !> Compact nuclear-gradient rows of one surface tile
   type :: tile_rows_type
      !> Atoms the tile touched, in first-touch order
      integer, allocatable :: atoms(:)
      !> Their gradient rows (3, n)
      real(wp), allocatable :: rows(:, :)
   end type tile_rows_type

contains

   !> Contract a surface adjoint into the nuclear gradient
   !>
   !> Accumulates `dE/dR_A` for the energy whose surface adjoints `acc` holds
   !>
   !> - the result is *added* to `gradient`, so several cavities or several
   !>   passes can share one accumulator
   !> - the reverse-mode counterpart of running `compute_gradient_drop` and
   !>   contracting the resulting `*_rA` arrays
   !> - the two agree to round-off, which `test_cavity_drop_nuclear_adjoint`
   !>   asserts
   !>
   !> @param[in]    self     DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc      Accumulated surface-observable adjoints
   !> @param[inout] gradient Nuclear-gradient accumulator (3, nsph)
   !> @param[out]   error    Error object, allocated on failure
   module subroutine get_surface_gradient_drop(self, acc, gradient, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Per-tile gradient rows, summed in tile order after the region
      type(tile_rows_type), allocatable :: tile_rows(:)
      !> Gradient of the current tile (3, nsph), its touched atoms and their mask
      real(wp), allocatable :: tile_grad(:, :)
      integer, allocatable :: row_atoms(:)
      integer :: nrow
      logical, allocatable :: in_row(:)
      !> Tile-ordered sum of the tile rows
      real(wp), allocatable :: grad_sum(:, :)
      !> Thread bookkeeping
      integer :: thread_slot
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort

      !> Shared per-grid point sensitivity kernel state and its response
      type(drop_seed_state_type) :: state
      !> Degeneracy status
      integer :: status

      !> Grid, tile-column, candidate, atom and anchor indices
      integer :: igrid, p, k, iatom, ia, jrow
      !> Tile index, tile count, tile width, first grid point and size
      integer :: itile, ntile, tile, igrid0, np
      !> Shared candidate list of the tile
      integer, allocatable :: cand(:)
      integer :: nc

      !> Anchor and owner sphere
      real(wp) :: anchor(3)
      integer :: owner_idx
      !> Tile points and level-set jet, packed Hessian and third derivative
      real(wp), allocatable :: tile_xyz(:, :), tile_lsf0(:), tile_lsf1(:, :)
      real(wp), allocatable :: tile_lsf2(:, :), tile_lsf3(:, :)
      !> Tile field seeds, packed symmetric Hessian seed
      real(wp), allocatable :: tile_w0(:), tile_w1(:, :), tile_w2(:, :)
      !> Seed-contracted nuclear partials of the tile (3, nc, np), over a flat
      !> buffer sized for the largest candidate list
      real(wp), allocatable, target :: tile_res_buf(:)
      real(wp), pointer, contiguous :: tile_res(:, :, :)
      !> Objective gradient at the projected point
      real(wp) :: phi1_r(3)

      !> Bordered KKT system with four field seeds and three anchor seeds
      real(wp) :: kkt_rhs(4, 7)

      !> Point-local level-set adjoints from the 10 field seeds, Hessian packed
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(6)
      !> Effective position adjoint seen by every seed
      real(wp) :: w_xyz_local(3)
      !> Summed switching weight of an anchor's rows, and whether any row is live
      real(wp) :: w_switch
      logical :: have_switch
      !> iSwig switching-gradient scratch
      type(iswig_workspace_type) :: iswig_work

      !> Effective primitive surface adjoints
      type(drop_surface_weights_type) :: eff
      !> Timer handle
      integer :: h_sgrad

      call self%require_context(error)
      if (allocated(error)) return

      call check_surface_adjoint(self, acc, "get_surface_gradient_drop", error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%nsph])) then
         call fatal_error(error, "get_surface_gradient_drop: gradient shape mismatch")
         return
      end if
      if (self%ngrid <= 0) return

      h_sgrad = self%ctx%timer%resolve("Surface gradient", self%ctx%timer%current(), &
                                       cat_gradient)
      call self%ctx%timer%start(h_sgrad)

      !* -------------------------- Effective surface weights ------------------------- *!
      call prepare_surface_weights(self, acc, .true., eff)

      !* -------------------------------- Thread setup -------------------------------- *!
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      tile = surface_tile_width
      ntile = (self%ngrid + tile - 1)/tile
      allocate (tile_rows(ntile))

      call abort%reset()

      ! Tiles of `tile` points: one batched jet and one batched nuclear contraction per
      ! tile around the per-point KKT and seed work; every tile keeps its own gradient rows
      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, &
      !$omp& itile, igrid0, np, igrid, p, k, iatom, cand, nc, state, status, &
      !$omp& anchor, owner_idx, tile_xyz, tile_lsf0, tile_lsf1, tile_lsf2, tile_lsf3, &
      !$omp& tile_w0, tile_w1, tile_w2, tile_res_buf, tile_res, &
      !$omp& phi1_r, kkt_rhs, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_local, ia, jrow, w_switch, have_switch, &
      !$omp& iswig_work, tile_grad, row_atoms, nrow, in_row)
      thread_slot = 1
      !$ thread_slot = omp_get_thread_num() + 1

      allocate (cand(max(1, slots%lsf(thread_slot)%lsf%ncenters)))
      allocate (tile_xyz(3, tile), tile_lsf0(tile), tile_lsf1(3, tile), tile_lsf2(6, tile))
      allocate (tile_lsf3(10, tile))
      allocate (tile_w0(tile), tile_w1(3, tile), tile_w2(6, tile))
      allocate (tile_res_buf(3*size(cand)*tile))
      call self%iswig%new_workspace(iswig_work)
      allocate (tile_grad(3, self%nsph), source=0.0_wp)
      allocate (row_atoms(self%nsph))
      allocate (in_row(self%nsph), source=.false.)
      nrow = 0

      !$omp do schedule(dynamic, 1)
      do itile = 1, ntile
         if (abort%requested) cycle

         ! Clear the rows the previous tile of this thread touched
         do k = 1, nrow
            in_row(row_atoms(k)) = .false.
            tile_grad(:, row_atoms(k)) = 0.0_wp
         end do
         nrow = 0
         igrid0 = (itile - 1)*tile + 1
         np = min(tile, self%ngrid - igrid0 + 1)
         tile_xyz(:, 1:np) = self%xyz(:, igrid0:igrid0 + np - 1)

         ! The surface jet of the whole tile over one shared candidate list
         call surface_tile_candidates(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand, nc)
         if (.not. surface_jet_tile(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand(1:nc), &
                                    igrid0, "get_surface_gradient_drop", 3, tile_lsf0, tile_lsf1, &
                                    tile_lsf2, tile_lsf3, abort)) cycle

         do p = 1, np
            if (abort%requested) cycle

            igrid = igrid0 + p - 1
            anchor = self%anchorxyz(:, self%anchor_of_row(igrid))
            owner_idx = self%owner(igrid)
            call fill_seed_state(self, igrid, eff%have_wk, state)

            !* ------------------ Field seeds -> level-set adjoints ------------------- *!
            ! Columns 5-7 of the bordered KKT system are the anchor seeds:
            ! moving the owner rigidly leaves the field untouched and drives the system
            ! through -d^2 phi/dr dR_owner = +alpha*I
            kkt_rhs(:, 5:7) = 0.0_wp
            kkt_rhs(1, 5) = self%param%phi_alpha
            kkt_rhs(2, 6) = self%param%phi_alpha
            kkt_rhs(3, 7) = self%param%phi_alpha
            call seed_point_jet(slots%phi(thread_slot), tile_xyz(:, p), owner_idx, tile_lsf1(:, p), &
                                tile_lsf2(:, p), tile_lsf3(:, p), self%lambda0(igrid), self%f_crit, &
                                self%f_foc, self%f_wleb, self%param%wleb_prune_level > 0, eff, &
                                igrid, state, kkt_rhs, phi1_r, w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, &
                                w_xyz_local, status)
            if (status /= seed_state_ok) then
               call abort%latch_status(status, igrid)
               cycle
            end if

            ! Packed symmetric Hessian seed, the convention of the vjp
            tile_w0(p) = w_lsf0_pt
            tile_w1(:, p) = w_lsf1_pt
            tile_w2(:, p) = w_lsf2_pt

            !* ------------------------ Anchor channel (owner) ------------------------ *!
            call touch_row(owner_idx, in_row, row_atoms, nrow)
            call seed_anchor(state, eff, igrid, phi1_r, kkt_rhs, w_xyz_local, &
                             tile_grad(:, owner_idx))

            !* ----------------------- iSwig switching channel ------------------------ *!
            ! f_i depends on the nuclear geometry alone and is the same for every
            ! row of an anchor, whose rows are contiguous:
            ! the first row carries the summed weight of the run, evaluated once at fixed width
            ! (anchor_xi0 has no nuclear dependence)
            ia = self%anchor_of_row(igrid)
            if (igrid > 1) then
               if (self%anchor_of_row(igrid - 1) == ia) cycle
            end if
            w_switch = 0.0_wp
            have_switch = .false.
            do jrow = igrid, self%ngrid
               if (self%anchor_of_row(jrow) /= ia) exit
               if (abs(eff%w_f(jrow)) > seed_weight_tol) then
                  w_switch = w_switch + eff%w_f(jrow)
                  have_switch = .true.
               end if
            end do
            if (have_switch) then
               call self%iswig%add_f1_rA(anchor, owner_idx, self%anchor_xi0(ia), w_switch, &
                                         iswig_work, tile_grad)
               do k = 1, iswig_work%n
                  call touch_row(iswig_work%k(k), in_row, row_atoms, nrow)
               end do
            end if
         end do
         ! A latched point voids the result, so its tile is not contracted
         if (abort%requested) cycle

         !* -------------- Field channel: contract with nuclear partials -------------- *!
         ! The level set contracts the jet indices itself: `vjp_f1_rA_batch`
         ! returns the nuclear-gradient rows already weighted by the seeds, so
         ! the (3, 3, 3, n_active) mixed third derivative is never materialized
         tile_res(1:3, 1:nc, 1:np) => tile_res_buf(1:3*nc*np)
         if (.not. surface_vjp_tile(slots%lsf(thread_slot)%lsf, tile_xyz(:, 1:np), cand(1:nc), &
                                    tile_w0, tile_w1, tile_w2, igrid0, "get_surface_gradient_drop", &
                                    tile_res, abort)) cycle
         ! Column k belongs to candidate cand(k), in the LSF's candidate space
         do k = 1, nc
            call touch_row(slots%lsf(thread_slot)%lsf%cand_to_user(cand(k)), in_row, &
                           row_atoms, nrow)
         end do
         do p = 1, np
            do k = 1, nc
               iatom = slots%lsf(thread_slot)%lsf%cand_to_user(cand(k))
               tile_grad(:, iatom) = tile_grad(:, iatom) + tile_res(:, k, p)
            end do
         end do

         ! Keep the touched rows only; their order within a tile does not
         ! matter, every atom is summed over the tiles in tile order
         allocate (tile_rows(itile)%atoms(nrow), tile_rows(itile)%rows(3, nrow))
         do k = 1, nrow
            tile_rows(itile)%atoms(k) = row_atoms(k)
            tile_rows(itile)%rows(:, k) = tile_grad(:, row_atoms(k))
         end do
      end do
      !$omp end do

      deallocate (cand, tile_xyz, tile_lsf0, tile_lsf1, tile_lsf2, tile_lsf3)
      deallocate (tile_w0, tile_w1, tile_w2, tile_res_buf)
      deallocate (tile_grad, row_atoms, in_row)
      !$omp end parallel

      if (abort%requested) then
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("get_surface_gradient_drop", abort%status, abort%igrid, error)
         end if
         call self%ctx%timer%stop(h_sgrad)
         return
      end if

      ! Deterministic reduction: ascending tile order
      allocate (grad_sum(3, self%nsph), source=0.0_wp)
      do itile = 1, ntile
         do k = 1, size(tile_rows(itile)%atoms)
            iatom = tile_rows(itile)%atoms(k)
            grad_sum(:, iatom) = grad_sum(:, iatom) + tile_rows(itile)%rows(:, k)
         end do
      end do
      gradient = gradient + grad_sum

      call self%ctx%timer%stop(h_sgrad)

   end subroutine get_surface_gradient_drop

   !> Record that a tile wrote the gradient row of an atom
   !>
   !> @param[in]    iatom     Atom whose row was written
   !> @param[inout] in_row    Atoms recorded in this tile (nsph)
   !> @param[inout] row_atoms Recorded atoms, in first-touch order
   !> @param[inout] nrow      Number of recorded atoms
   pure subroutine touch_row(iatom, in_row, row_atoms, nrow)
      !> Atom whose row was written
      integer, intent(in) :: iatom
      !> Atoms recorded in this tile
      logical, intent(inout) :: in_row(:)
      !> Recorded atoms
      integer, intent(inout) :: row_atoms(:)
      !> Number of recorded atoms
      integer, intent(inout) :: nrow

      if (in_row(iatom)) return
      in_row(iatom) = .true.
      nrow = nrow + 1
      row_atoms(nrow) = iatom
   end subroutine touch_row

end submodule moist_cavity_drop_derivatives_nuclear
