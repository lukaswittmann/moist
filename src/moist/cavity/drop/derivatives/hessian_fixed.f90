!> Fixed-adjoint half of the DROP surface Hessian
!>
!> The nuclear gradient of `nuclear.f90` is `J^T omega` with `omega` the
!> accumulated surface adjoint. Its directional derivative splits into
!>
!>     d/dv [ J^T omega ]  =  (dJ^T/dv) omega  +  J^T (d omega/dv)
!>
!> and this submodule builds the first term: the second derivative of the
!> surface map contracted against adjoints held *fixed*. The adjoint-response
!> term belongs to the forward-tangent and weight-tangent passes and is not
!> computed here.
!>
!> Because the adjoints are fixed, the term is **direction free** -- the second
!> derivative of the geometry does not know which direction it is later
!> contracted with -- so the result is accumulated into a rank-4
!> `(3, nsph, 3, nsph)` object rather than into `ndir` gradient columns. Any
!> Hessian-vector product is a contraction of that; see the rank-4 note under
!> "Open questions" in `hessian.md` for the decision and its large-system cost.
!>
!>
!> Per grid point the three channels of the gradient are differentiated as:
!>
!>  1. Field: the row `vjp_f1_rA(w_lsf0, w_lsf1, w_lsf2)` moves both because
!>     the adjoint weights move and because the level-set jet is evaluated at
!>     the projected point, which rides along. [[drop_field_tangent]] owns that
!>     whole block; this routine owns the weight tangents feeding it
!>
!>  2. Anchor: the owner's three rigid directions, differentiated through the
!>     same seed chain
!>
!>  3. Switching: the iSwiG `f_i` depends on the nuclear geometry alone and its
!>     adjoint `eff%w_f` is fixed, so the entire channel is one weighted
!>     `swi2_rArB_block` scattered over the influence set -- with no loop over
!>     directions at all
!>
!> Channels 1 and 2 do have a direction loop, but a *local* one: the projected
!> point depends on the nuclei only through the level set (its active atoms)
!> and through the anchor (its owner sphere), so the loop runs over
!> `3 * |active union owner|` basis directions and the result is still the
!> direction-free block.
!>
!>
!> **Scope limit -- geometry-dependent effective weights are rejected.**
!> [[prepare_surface_weights]] folds the area and integration-weight channels
!> into `w_xi` and `w_f` through `self%a`, `self%wleb`, `self%xi0` and the
!> radii, and derives `branch_phi_adj` from the branch softmax; all of those are
!> geometry dependent, so differentiating them is the weight-tangent pass. Until
!> that pass exists this routine refuses an accumulator that would need it --
!> a live `w_a` or `w_w` channel, or a grid carrying multi-branch anchor groups
!> -- rather than returning a Hessian that is silently missing those terms.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_hessian_fixed
!$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, drop_seed_state_tangent_type, &
      & drop_seed_input_tangent_type, drop_seed_result_tangent_type, &
      & drop_surface_weights_type, build_seed_state, apply_seed, apply_seed_tangent, &
      & seed_state_ok, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type, drop_n_jet_seeds, &
      & seed_normal_channel, degenerate_point_error
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent, &
      & drop_field_tangent_work_type
   use moist_cavity_drop_derivatives_iswig_scatter, only: scatter_iswig_block
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Seeds pushed through the per-point kernel: the 13 level-set jet
   !> directions of [[seed_jet_basis]] followed by the three anchor directions
   !> of [[seed_anchor]], in that order
   integer, parameter :: n_point_seeds = drop_n_jet_seeds + 3

   !> Vanishing tangent of a basis seed. Every seed [[fill_seed_basis]] emits is
   !> a constant matrix, so its own derivative along a nuclear direction is zero
   !> and only the point motion it induces survives; these are the arguments
   !> that carries into [[apply_seed_tangent]]
   real(wp), parameter :: seed_dzero1(ndim) = 0.0_wp
   real(wp), parameter :: seed_dzero2(ndim, ndim) = 0.0_wp

contains

   !> Contract the second derivative of the surface map into a nuclear Hessian
   !>
   !> Accumulates `(dJ^T/dv) omega` for the energy whose surface adjoints `acc`
   !> holds. The result is *added* to `hessian`, so several cavities or several
   !> passes can share one accumulator, and the accumulator is left untouched
   !> when anything fails.
   !>
   !> Mirrors [[get_surface_gradient_drop]] throughout: same thread setup, same
   !> effective weights, same grid loop, same error latching and the same
   !> deterministic reduction. The differences are the derivative order
   !> (`max_deriv(4)`, for the fourth-order jet and mixed tensors), the local
   !> direction loop, and the rank-4 accumulator.
   !>
   !> @param[in]    self    DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc     Accumulated surface-observable adjoints, held fixed
   !> @param[inout] hessian Nuclear-Hessian accumulator (3, nsph, 3, nsph)
   !> @param[out]   error   Error object, allocated on failure
   module subroutine get_surface_hessian_fixed_drop(self, acc, hessian, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Per-thread Hessian buffers, summed deterministically after the region
      real(wp), allocatable :: hess_threads(:, :, :, :, :)
      !> Thread bookkeeping
      integer :: thread_slot, ithread
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure on its way to the latch
      type(error_type), allocatable :: worker_error

      !> Shared per-grid point sensitivity kernel state
      type(drop_seed_state_type) :: state
      !> Degeneracy status
      integer :: status

      !> Grid, atom, axis, seed and active-slot indices
      integer :: igrid, iatom, i, k, iaxis, ibasis
      integer :: n_active, ndir_atom, idir, dir_atom, dir_axis
      integer, allocatable :: active_idx(:), dir_atoms(:)

      !> Projected point, anchor and owner sphere
      real(wp) :: point(3), anchor(3)
      integer :: owner_idx
      !> Level-set jet at the projected point
      real(wp) :: lsf0, lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :), lsf4_rrrr(:, :, :, :)
      !> Objective jet at the projected point
      real(wp) :: phi0, phi1_r(3), phi2_rr(3, 3)
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val

      !> Bordered KKT system with four field seeds and three anchor seeds
      real(wp) :: kkt_rhs(4, 7)
      !> Factorization reused by every solve at this grid point
      type(drop_kkt_factor_type) :: kkt_fac

      !> The 16 basis seeds and the responses they produce
      real(wp) :: seed_dlsf1(3, n_point_seeds), seed_dlsf2(3, 3, n_point_seeds)
      real(wp) :: seed_x(4, n_point_seeds)
      type(drop_seed_result_type), allocatable :: res_seed(:)
      type(drop_seed_state_tangent_type), allocatable :: dstate_seed(:)

      !> Point-local level-set adjoint weights built from the 13 field seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Effective position adjoint seen by every seed, and the normal fold it
      !> is built from
      real(wp) :: w_xyz_local(3), normal_grad(3), nwn

      !> Nuclear tangent of the level-set jet at the *frozen* projected point
      real(wp) :: dv0, dv1(3), dv2(3, 3)
      real(wp), allocatable :: dv3(:, :, :)
      !> Right-hand side and solution of the directional projection response
      real(wp) :: rhs_v(4, 1), dr_v(3), dl_v
      !> Response, derived-state tangent and input tangent of one direction
      type(drop_seed_result_type) :: res_v
      type(drop_seed_state_tangent_type) :: dstate_v
      type(drop_seed_input_tangent_type) :: dinp_v
      !> Second-order response of one seed along one direction
      type(drop_seed_result_tangent_type) :: dres
      !> Tangent of the bordered KKT system and of its 16 right-hand sides
      real(wp) :: dH_lag(3, 3, 1), dg_tot(3, 1)
      real(wp) :: dseed_x(4, n_point_seeds)
      !> Tangents of the normal fold and of the level-set adjoint weights
      real(wp) :: dnormal_grad(3), dw_xyz_local(3)
      real(wp) :: dw_lsf0, dw_lsf1(3), dw_lsf2(3, 3)
      !> One seed's contribution to a weight tangent
      real(wp) :: contribution
      !> Basis direction and the field row it produces
      real(wp), allocatable :: vdir(:, :), field_row(:, :)
      !> Scratch of the field-contraction tangent, reused across points
      type(drop_field_tangent_work_type) :: ft_work

      !> iSwig neighbour cache and the local second-derivative block it feeds
      type(iswig_workspace_type) :: iswig_work
      real(wp), allocatable :: swi_blk(:, :, :, :), swi_mix(:, :)
      integer, allocatable :: swi_idx(:)
      real(wp) :: swi_f0, swi_d2xi
      integer :: swi_n

      !> Effective primitive surface adjoints
      type(drop_surface_weights_type) :: eff
      !> Timer handle
      integer :: h_shess

      call check_surface_adjoint(self, acc, "get_surface_hessian_fixed_drop", error)
      if (allocated(error)) return
      if (any(shape(hessian) /= [3, self%nsph, 3, self%nsph])) then
         call fatal_error(error, "get_surface_hessian_fixed_drop: hessian shape mismatch")
         return
      end if
      call check_frozen_weights(self, acc, error)
      if (allocated(error)) return
      if (self%ngrid <= 0) return

      h_shess = self%ctx%timer%resolve("Surface Hessian (fixed adjoint)", &
                                       self%ctx%timer%current(), cat_gradient)
      call self%ctx%timer%start(h_shess)

      !* -------------------------- Effective surface weights ------------------------- *!
      call prepare_surface_weights(self, acc, .true., eff)

      !* -------------------------------- Thread setup -------------------------------- *!
      ! Order 4 rather than the gradient path's 3: the field tangent reads
      ! `f4_rrrr` and `f4_rrr_rA`, and CFC asks for the highest *total* order.
      call slots%init(self%ctx, self%lsf_model, 4, self%param, self%mol, self%radii)
      allocate (hess_threads(3, self%nsph, 3, self%nsph, slots%nthreads), source=0.0_wp)

      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& iatom, i, k, iaxis, ibasis, n_active, ndir_atom, idir, dir_atom, dir_axis, &
      !$omp& active_idx, dir_atoms, state, status, &
      !$omp& point, anchor, owner_idx, lsf0, lsf1_r, lsf2_rr, lsf3_rrr, lsf4_rrrr, &
      !$omp& phi0, phi1_r, phi2_rr, lambda_val, kkt_rhs, kkt_fac, &
      !$omp& seed_dlsf1, seed_dlsf2, seed_x, res_seed, dstate_seed, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_local, normal_grad, nwn, &
      !$omp& dv0, dv1, dv2, dv3, rhs_v, dr_v, dl_v, res_v, dstate_v, dinp_v, dres, &
      !$omp& dH_lag, dg_tot, dseed_x, dnormal_grad, dw_xyz_local, &
      !$omp& dw_lsf0, dw_lsf1, dw_lsf2, contribution, vdir, field_row, ft_work, &
      !$omp& iswig_work, swi_blk, swi_mix, swi_idx, swi_f0, swi_d2xi, swi_n, worker_error)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (lsf4_rrrr(3, 3, 3, 3), source=0.0_wp)
      allocate (dv3(3, 3, 3), source=0.0_wp)
      allocate (active_idx(self%nsph))
      allocate (dir_atoms(self%nsph))
      allocate (vdir(3, self%nsph), source=0.0_wp)
      allocate (field_row(3, self%nsph), source=0.0_wp)
      allocate (res_seed(n_point_seeds))
      allocate (dstate_seed(n_point_seeds))
      call iswig_work%init(self%iswig)
      ! Grown on demand from `n_nb + 1` rather than sized to `nsph`: unlike the
      ! sparse rows of the gradient path a block is quadratic in the influence
      ! set, and a molecule-sized one would cost as much per thread as the
      ! rank-4 accumulator itself.
      allocate (swi_blk(3, 1, 3, 1), swi_mix(3, 1), swi_idx(1))

      !$omp do schedule(static, 8)
      do igrid = 1, self%ngrid
         if (abort%requested) cycle

         point = self%xyz(:, igrid)
         anchor = self%anchorxyz(:, igrid)
         owner_idx = self%owner(igrid)
         lambda_val = self%lambda0(igrid)

         call slots%lsf(thread_slot)%lsf%prepare(point, worker_error)

         ! The failure cannot be returned from inside this worksharing construct,
         ! so park it for the post-region promotion and let the flag drain the
         ! loop. The LSF's cached derivatives are substitutes; stop before
         ! reading them.
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            cycle
         end if

         call slots%lsf(thread_slot)%lsf%f3_rrr(lsf0, lsf1_r, lsf2_rr, lsf3_rrr)
         call slots%lsf(thread_slot)%lsf%f4_rrrr(lsf4_rrrr)
         call slots%phi(thread_slot)%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         state%lsf1_r = lsf1_r
         state%lsf2_rr = lsf2_rr
         state%lsf3_rrr = lsf3_rrr
         state%lambda_val = lambda_val
         call fill_seed_state(self, igrid, eff%have_wk, state)

         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)
         if (status /= seed_state_ok) then
            call abort%latch_status(status, igrid)
            cycle
         end if

         ! Outward-normal channel, as in the gradient path. `normal_grad` is
         ! kept here as well as folded, because the direction loop needs its
         ! own tangent and [[seed_normal_channel]] does not hand it back.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         normal_grad = 0.0_wp
         nwn = 0.0_wp
         if (eff%have_wn) then
            nwn = dot_product(state%n_surf, eff%w_n(:, igrid))
            normal_grad = (eff%w_n(:, igrid) - state%n_surf*nwn)/state%g_norm
         end if
         call seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_local)

         !* ------------------------ Bordered KKT sensitivities ----------------------- *!
         ! Columns 1-4 are the level-set value and gradient seeds; the nine
         ! Hessian seeds have a zero right-hand side. Columns 5-7 are the
         ! anchor seeds: moving the owner rigidly leaves the field untouched
         ! and drives the system through -d^2 phi/dr dR_owner = +alpha*I.
         kkt_rhs = 0.0_wp
         kkt_rhs(4, 1) = -1.0_wp
         kkt_rhs(1, 2) = lambda_val
         kkt_rhs(2, 3) = lambda_val
         kkt_rhs(3, 4) = lambda_val
         kkt_rhs(1, 5) = self%param%phi_alpha
         kkt_rhs(2, 6) = self%param%phi_alpha
         kkt_rhs(3, 7) = self%param%phi_alpha
         call kkt_fac%factor(phi2_rr - lambda_val*lsf2_rr, lsf1_r, worker_error)
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            cycle
         end if
         call kkt_fac%solve(kkt_rhs, worker_error)
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            cycle
         end if

         !* --------------------- Basis seeds and their responses --------------------- *!
         ! The 16 seeds are those of [[seed_jet_basis]] and [[seed_anchor]],
         ! collected into one array because the second-order chain needs each
         ! seed's `res` and `dstate` again inside the direction loop; rebuilding
         ! them per direction would be the dominant avoidable cost, and a second
         ! copy of a floating-point chain is free to contract differently.
         call fill_seed_basis(kkt_rhs, seed_dlsf1, seed_dlsf2, seed_x)
         do ibasis = 1, n_point_seeds
            call apply_seed(state, seed_dlsf1(:, ibasis), seed_dlsf2(:, :, ibasis), &
                            seed_x(1:3, ibasis), seed_x(4, ibasis), &
                            res_seed(ibasis), dstate_seed(ibasis))
            contribution = seed_contribution(eff, igrid, w_xyz_local, &
                                             seed_x(1:3, ibasis), res_seed(ibasis))
            call scatter_jet_weight(ibasis, contribution, w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)
         end do

         !* --------------------------- Local direction set --------------------------- *!
         ! The projected point, the jet and the anchor depend on the nuclei
         ! only through the level set's active atoms and through the owner, so
         ! every other column of this point's Hessian block is exactly zero.
         n_active = slots%lsf(thread_slot)%lsf%active_count()
         do i = 1, n_active
            active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
         end do
         ndir_atom = n_active
         dir_atoms(1:n_active) = active_idx(1:n_active)
         if (.not. any(dir_atoms(1:ndir_atom) == owner_idx)) then
            ndir_atom = ndir_atom + 1
            dir_atoms(ndir_atom) = owner_idx
         end if

         do idir = 1, 3*ndir_atom
            dir_atom = dir_atoms((idir - 1)/3 + 1)
            dir_axis = mod(idir - 1, 3) + 1
            vdir(dir_axis, dir_atom) = 1.0_wp

            !* ------------- Directional response of the projected point -------------- *!
            ! Tangent of the jet at the *frozen* point first; the projection
            ! then rides on it through the same bordered system the seeds use,
            ! with the anchor moving rigidly with its owner.
            call slots%lsf(thread_slot)%lsf%tangent_f0(vdir, dv0)
            call slots%lsf(thread_slot)%lsf%tangent_f1_r(vdir, dv1)
            call slots%lsf(thread_slot)%lsf%tangent_f2_rr(vdir, dv2)
            call slots%lsf(thread_slot)%lsf%tangent_f3_rrr(vdir, dv3)

            rhs_v = 0.0_wp
            rhs_v(1:3, 1) = lambda_val*dv1
            if (dir_atom == owner_idx) then
               rhs_v(dir_axis, 1) = rhs_v(dir_axis, 1) + self%param%phi_alpha
            end if
            rhs_v(4, 1) = -dv0
            call kkt_fac%solve(rhs_v, worker_error)
            if (allocated(worker_error)) then
               call abort%latch_error(worker_error, igrid)
               vdir(dir_axis, dir_atom) = 0.0_wp
               exit
            end if
            dr_v = rhs_v(1:3, 1)
            dl_v = rhs_v(4, 1)

            !* ----------------- Directional state and input tangents ----------------- *!
            ! A nuclear direction is just another seed of the same linear map,
            ! so `apply_seed` produces the whole forward tangent of this point.
            call apply_seed(state, dv1, dv2, dr_v, dl_v, res_v, dstate_v)

            ! The input tangents are *total*: `state%lsf1_r` is grad S at the
            ! projected point, so its v-tangent carries the point motion, which
            ! is exactly what `res_v%dg` and `res_v%dH` already are. The third
            ! spatial derivative has no such accessor and is folded by hand.
            dinp_v%dlsf1_r = res_v%dg
            dinp_v%dlsf2_rr = res_v%dH
            dinp_v%dlsf3_rrr = dv3
            do k = 1, 3
               dinp_v%dlsf3_rrr(:, :, :) = dinp_v%dlsf3_rrr(:, :, :) &
                                           + lsf4_rrrr(:, :, :, k)*dr_v(k)
            end do
            dinp_v%dlambda_val = dl_v
            dinp_v%dcpjac_scal0 = res_v%dJ
            dinp_v%dw_f0 = res_v%dw_f
            dinp_v%dwleb = res_v%dwleb
            dinp_v%dxi0 = res_v%dxi
            ! The anchor's Lebedev weight is a property of the rigid sphere and
            ! the branch weight is one for every group this routine admits, so
            ! both tangents vanish; see the scope limit in the module header.
            dinp_v%danchor_wleb0 = 0.0_wp
            dinp_v%dwbranch = 0.0_wp

            !* ---------------------- Tangent of the normal fold ---------------------- *!
            dnormal_grad = 0.0_wp
            dw_xyz_local = 0.0_wp
            if (eff%have_wn) then
               dnormal_grad = (-res_v%dn_surf*nwn &
                               - state%n_surf*dot_product(res_v%dn_surf, eff%w_n(:, igrid))) &
                              /state%g_norm &
                              - normal_grad*res_v%d_gnorm/state%g_norm
               dw_xyz_local = matmul(res_v%dH, normal_grad) + matmul(lsf2_rr, dnormal_grad)
            end if

            !* ----------------- Tangent of the seed right-hand sides ----------------- *!
            ! `K dx = db - dK x` on the same factors. Only the three gradient
            ! seeds carry a `db`: their right-hand side is the multiplier, and
            ! `alpha` and the unit value seed are constants.
            dH_lag(:, :, 1) = -dl_v*lsf2_rr - lambda_val*res_v%dH
            dg_tot(:, 1) = res_v%dg
            dseed_x = 0.0_wp
            do iaxis = 1, 3
               dseed_x(iaxis, 1 + iaxis) = dl_v
            end do
            call kkt_fac%solve_tangent(dH_lag, dg_tot, seed_x, dseed_x, worker_error)
            if (allocated(worker_error)) then
               call abort%latch_error(worker_error, igrid)
               vdir(dir_axis, dir_atom) = 0.0_wp
               exit
            end if

            !* ---------------------------- Weight tangents --------------------------- *!
            ! Every basis seed is a constant matrix, so its own tangent along a
            ! nuclear direction vanishes and only the induced point motion moves.
            dw_lsf0 = 0.0_wp
            dw_lsf1 = dnormal_grad
            dw_lsf2 = 0.0_wp
            do ibasis = 1, n_point_seeds
               call apply_seed_tangent(state, dstate_v, dinp_v, res_v, &
                                       seed_dlsf1(:, ibasis), seed_dlsf2(:, :, ibasis), &
                                       seed_x(1:3, ibasis), seed_x(4, ibasis), &
                                       seed_dzero1, seed_dzero2, &
                                       dseed_x(1:3, ibasis), dseed_x(4, ibasis), &
                                       res_seed(ibasis), dstate_seed(ibasis), dres)
               contribution = seed_contribution_tangent(eff, igrid, w_xyz_local, &
                                                        dw_xyz_local, seed_x(1:3, ibasis), &
                                                        dseed_x(1:3, ibasis), dres)
               if (ibasis > drop_n_jet_seeds) then
                  ! Anchor seed: the gradient row it feeds belongs to the owner.
                  iaxis = ibasis - drop_n_jet_seeds
                  hess_threads(iaxis, owner_idx, dir_axis, dir_atom, thread_slot) = &
                     hess_threads(iaxis, owner_idx, dir_axis, dir_atom, thread_slot) &
                     + contribution
               else
                  call scatter_jet_weight(ibasis, contribution, dw_lsf0, dw_lsf1, dw_lsf2)
               end if
            end do

            !* ----------------------------- Field channel ---------------------------- *!
            ! `dw_lsf2` is read only here, and only against `lsf3_rr_rA`, which
            ! is symmetric in its two spatial indices. That is what makes the
            ! nine asymmetric Hessian seeds legitimate: [[apply_seed_tangent]]
            ! is exact only on a symmetric seed, its per-seed error is
            ! antisymmetric in `(i, j)`, and an antisymmetric error contracted
            ! with a symmetric weight cancels to machine zero. Never read a
            ! single off-diagonal entry of `dw_lsf2` on its own.
            if (n_active > 0) then
               call drop_field_tangent(slots%lsf(thread_slot)%lsf, &
                                       w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, &
                                       dw_lsf0, dw_lsf1, dw_lsf2, dr_v, vdir, &
                                       ft_work, field_row)
               do i = 1, n_active
                  iatom = active_idx(i)
                  hess_threads(:, iatom, dir_axis, dir_atom, thread_slot) = &
                     hess_threads(:, iatom, dir_axis, dir_atom, thread_slot) &
                     + field_row(:, i)
               end do
            end if

            vdir(dir_axis, dir_atom) = 0.0_wp
         end do
         if (abort%requested) cycle

         !* ------------------------- iSwig switching channel ------------------------- *!
         ! `f_i` depends on the nuclear geometry alone and its adjoint is fixed,
         ! so the whole channel is one block over the influence set. The width
         ! rows `swi_mix` and `swi_d2xi` are not read: the width the switching
         ! factor is evaluated at is the *anchor* width, built from the rigid
         ! sphere's Lebedev weight, and so is nuclear-geometry independent --
         ! which is why the gradient path drops `swi_dxi` as well.
         if (abs(eff%w_f(igrid)) > seed_weight_tol) then
            call self%iswig%swi_collect(anchor, owner_idx, self%anchor_xi0(igrid), &
                                        swi_f0, iswig_work)
            if (size(swi_idx) < iswig_work%n_nb + 1) then
               deallocate (swi_blk, swi_mix, swi_idx)
               allocate (swi_blk(3, iswig_work%n_nb + 1, 3, iswig_work%n_nb + 1))
               allocate (swi_mix(3, iswig_work%n_nb + 1))
               allocate (swi_idx(iswig_work%n_nb + 1))
            end if
            call self%iswig%swi2_rArB_block(iswig_work, swi_n, swi_idx, swi_blk, &
                                            swi_mix, swi_d2xi)
            call scatter_iswig_block(swi_n, swi_idx, swi_blk, eff%w_f(igrid), &
                                     hess_threads(:, :, :, :, thread_slot))
         end if

      end do
      !$omp end do

      deallocate (lsf3_rrr, lsf4_rrrr, dv3, active_idx, dir_atoms, vdir, field_row)
      deallocate (res_seed, dstate_seed, swi_blk, swi_mix, swi_idx)
      call iswig_work%destroy()
      !$omp end parallel

      if (abort%requested) then
         ! An LSF failure or a KKT failure arrives as a ready-made error; a
         ! kernel degeneracy arrives as a status code that needs this routine's
         ! name to become a diagnostic.
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("get_surface_hessian_fixed_drop", abort%status, &
                                        abort%igrid, error)
         end if
         call self%ctx%timer%stop(h_shess)
         return
      end if

      ! Deterministic reduction: fixed thread order, independent of scheduling
      do ithread = 1, slots%nthreads
         hessian = hessian + hess_threads(:, :, :, :, ithread)
      end do

      call self%ctx%timer%stop(h_shess)

   end subroutine get_surface_hessian_fixed_drop

   !* ================================================================================= *!
   !*                              Scope of the fixed half                               *!
   !* ================================================================================= *!

   !> Reject an accumulator whose effective weights are geometry dependent
   !>
   !> [[prepare_surface_weights]] is the boundary this routine cannot cross on
   !> its own. Two of its outputs move with the nuclei:
   !>
   !>   * the area and integration-weight folds, which convert `w_a` and `w_w`
   !>     into `w_xi` and `w_f` through `self%a`, `self%wleb`, `self%xi0` and
   !>     the radii, and
   !>   * `branch_phi_adj`, which [[compute_branch_phi_adj]] derives from the
   !>     branch softmax -- and which is identically zero unless some anchor
   !>     group carries more than one branch.
   !>
   !> Differentiating either is the weight-tangent pass, not this one. Both are
   !> refused rather than silently omitted: a caller cannot tell a Hessian that
   !> is missing a term from one that is complete.
   !>
   !> The thresholds match [[prepare_surface_weights]]'s own, so a channel it
   !> would not fold is not rejected here either.
   !>
   !> @param[in]  self  DROP cavity instance
   !> @param[in]  acc   Accumulated surface-observable adjoints
   !> @param[out] error Error object, allocated when a weight would move
   subroutine check_frozen_weights(self, acc, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (any(abs(acc%w_a) > seed_weight_tol)) then
         call fatal_error(error, "get_surface_hessian_fixed_drop: the area channel folds"// &
                          " into geometry-dependent effective weights; its tangent is the"// &
                          " weight-tangent pass and is not implemented here")
         return
      end if
      if (any(abs(acc%w_w) > seed_weight_tol)) then
         call fatal_error(error, "get_surface_hessian_fixed_drop: the integration-weight"// &
                          " channel folds into geometry-dependent effective weights; its"// &
                          " tangent is the weight-tangent pass and is not implemented here")
         return
      end if
      if (allocated(self%branch_count)) then
         if (any(self%branch_count(1:self%ngrid) > 1)) then
            call fatal_error(error, "get_surface_hessian_fixed_drop: multi-branch anchor"// &
                             " groups make the branch adjoint geometry dependent; its"// &
                             " tangent is the weight-tangent pass and is not implemented here")
            return
         end if
      end if
   end subroutine check_frozen_weights

   !* ================================================================================= *!
   !*                            Basis seeds of one grid point                          *!
   !* ================================================================================= *!

   !> Lay out the 16 basis seeds in the order the second-order chain reads them
   !>
   !> Slots 1-13 are the level-set jet directions of [[seed_jet_basis]] -- the
   !> value, the three gradient components and the nine Hessian components, the
   !> last of which are single-entry matrices and therefore asymmetric -- and
   !> slots 14-16 the three anchor directions of [[seed_anchor]]. `x` collects
   !> the induced point motion and multiplier change of every seed in the
   !> `(4, nseed)` layout [[drop_kkt_solve_tangent]] expects; the nine Hessian
   !> seeds move neither, because the stationarity conditions see only `S` and
   !> `grad S`.
   !>
   !> @param[in]  kkt      Solved KKT sensitivities, columns 1-4 jet, 5-7 anchor
   !> @param[out] dlsf1_r  Gradient perturbation of each seed
   !> @param[out] dlsf2_rr Hessian perturbation of each seed
   !> @param[out] x        Induced point motion and multiplier change of each seed
   pure subroutine fill_seed_basis(kkt, dlsf1_r, dlsf2_rr, x)
      !> Solved KKT sensitivities
      real(wp), intent(in) :: kkt(:, :)
      !> Seed perturbations of the level-set jet
      real(wp), intent(out) :: dlsf1_r(ndim, n_point_seeds)
      real(wp), intent(out) :: dlsf2_rr(ndim, ndim, n_point_seeds)
      !> Induced point motion and multiplier change
      real(wp), intent(out) :: x(4, n_point_seeds)

      !> Seed and Cartesian indices
      integer :: ibasis, iaxis, jaxis

      dlsf1_r = 0.0_wp
      dlsf2_rr = 0.0_wp
      x = 0.0_wp

      do ibasis = 1, drop_n_jet_seeds
         if (ibasis <= 4) then
            if (ibasis > 1) dlsf1_r(ibasis - 1, ibasis) = 1.0_wp
            x(:, ibasis) = kkt(1:4, ibasis)
         else
            iaxis = (ibasis - 5)/3 + 1
            jaxis = mod(ibasis - 5, 3) + 1
            dlsf2_rr(iaxis, jaxis, ibasis) = 1.0_wp
         end if
      end do

      do iaxis = 1, ndim
         x(:, drop_n_jet_seeds + iaxis) = kkt(1:4, 4 + iaxis)
      end do
   end subroutine fill_seed_basis

   !> Place one jet seed's contribution in the level-set adjoint weights
   !>
   !> The same 13-slot layout [[fill_seed_basis]] writes, read back: slot 1 is
   !> the value, slots 2-4 the gradient and slots 5-13 the Hessian in row-major
   !> order. Anchor slots are not accepted; their contribution belongs to the
   !> owner's gradient row, not to a weight.
   !>
   !> @param[in]    ibasis Seed slot, `1 .. drop_n_jet_seeds`
   !> @param[in]    contrib Contribution of that seed
   !> @param[inout] w0     Level-set value adjoint
   !> @param[inout] w1     Level-set gradient adjoint
   !> @param[inout] w2     Level-set Hessian adjoint
   pure subroutine scatter_jet_weight(ibasis, contrib, w0, w1, w2)
      !> Seed slot
      integer, intent(in) :: ibasis
      !> Contribution of that seed
      real(wp), intent(in) :: contrib
      !> Level-set adjoints
      real(wp), intent(inout) :: w0, w1(ndim), w2(ndim, ndim)

      !> Cartesian indices of a Hessian slot
      integer :: iaxis, jaxis

      if (ibasis == 1) then
         w0 = w0 + contrib
      else if (ibasis <= 4) then
         w1(ibasis - 1) = w1(ibasis - 1) + contrib
      else if (ibasis <= drop_n_jet_seeds) then
         iaxis = (ibasis - 5)/3 + 1
         jaxis = mod(ibasis - 5, 3) + 1
         w2(iaxis, jaxis) = w2(iaxis, jaxis) + contrib
      end if
   end subroutine scatter_jet_weight

   !* ================================================================================= *!
   !*                        Adjoint contraction of one seed                            *!
   !* ================================================================================= *!

   !> Adjoint contribution of one seed
   !>
   !> The contraction [[seed_jet_basis]] and [[seed_anchor]] both perform, with
   !> the branch term left out: this submodule admits only grids on which
   !> `branch_phi_adj` vanishes identically, and [[check_frozen_weights]]
   !> enforces that. The switching channel is absent for the same reason it is
   !> absent there -- the iSwiG overlap is anchor-only, so a level-set
   !> perturbation at fixed nuclei leaves it alone, and the anchor's own motion
   !> is carried by the block of the switching channel instead.
   !>
   !> @param[in] eff      Folded surface adjoints
   !> @param[in] igrid    Grid point
   !> @param[in] w_xyz_pt Effective position adjoint
   !> @param[in] dr       Induced point motion of the seed
   !> @param[in] res      Linear response of the seed
   !> @return             Adjoint contribution
   pure function seed_contribution(eff, igrid, w_xyz_pt, dr, res) result(contribution)
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Effective position adjoint
      real(wp), intent(in) :: w_xyz_pt(ndim)
      !> Induced point motion
      real(wp), intent(in) :: dr(ndim)
      !> Linear response
      type(drop_seed_result_type), intent(in) :: res
      !> Adjoint contribution
      real(wp) :: contribution

      contribution = dot_product(w_xyz_pt, dr) + eff%w_xi(igrid)*res%dxi
      if (eff%have_wk) then
         contribution = contribution + eff%w_k1(igrid)*res%dk1 + eff%w_k2(igrid)*res%dk2
      end if
   end function seed_contribution

   !> Directional derivative of [[seed_contribution]]
   !>
   !> Term by term the product rule applied to the contraction above, with the
   !> surface adjoints themselves held fixed -- which is the whole premise of
   !> this half. The position adjoint still moves, because the normal fold
   !> inside it is built from the level-set gradient at the projected point.
   !>
   !> @param[in] eff       Folded surface adjoints
   !> @param[in] igrid     Grid point
   !> @param[in] w_xyz_pt  Effective position adjoint
   !> @param[in] dw_xyz_pt Tangent of the effective position adjoint
   !> @param[in] dr        Induced point motion of the seed
   !> @param[in] ddr       Tangent of that point motion
   !> @param[in] dres      Second-order response of the seed
   !> @return              Tangent of the adjoint contribution
   pure function seed_contribution_tangent(eff, igrid, w_xyz_pt, dw_xyz_pt, dr, ddr, dres) &
      result(contribution)
      !> Folded surface adjoints
      type(drop_surface_weights_type), intent(in) :: eff
      !> Grid point
      integer, intent(in) :: igrid
      !> Effective position adjoint and its tangent
      real(wp), intent(in) :: w_xyz_pt(ndim), dw_xyz_pt(ndim)
      !> Induced point motion and its tangent
      real(wp), intent(in) :: dr(ndim), ddr(ndim)
      !> Second-order response
      type(drop_seed_result_tangent_type), intent(in) :: dres
      !> Tangent of the adjoint contribution
      real(wp) :: contribution

      contribution = dot_product(dw_xyz_pt, dr) + dot_product(w_xyz_pt, ddr) &
                     + eff%w_xi(igrid)*dres%dxi
      if (eff%have_wk) then
         contribution = contribution + eff%w_k1(igrid)*dres%dk1 + eff%w_k2(igrid)*dres%dk2
      end if
   end function seed_contribution_tangent

end submodule moist_cavity_drop_derivatives_hessian_fixed
