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
!> Hessian-vector product is a contraction of that.
!>
!> The trade behind that choice: one traversal of the second-order chain serves
!> every direction, and the price is an accumulator *per thread*
!> (`hess_threads` below). Streaming into `ndir` gradient columns instead would
!> bound that memory, but the fixed half has no direction to specialise on, so
!> it would have to re-run the whole chain once per direction.
!>
!> The per-thread accumulator is *sparse* -- see [[drop_hess_sparse_type]] --
!> because a dense one is `9 nsph^2` doubles per thread and grows
!> quadratically: 1.15 GB at a thousand atoms on sixteen threads, long before
!> the grid costs anything comparable. What a grid point actually reaches is
!> its own neighbourhood, so the accumulator holds only the atom pairs the
!> traversal touches.
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
   use, intrinsic :: iso_fortran_env, only: int64
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, seed_status_message, &
      & drop_seed_result_type, drop_seed_state_tangent_type, &
      & drop_seed_input_tangent_type, drop_seed_result_tangent_type, &
      & drop_surface_weights_type, apply_seed, apply_seed_tangent, &
      & seed_weight_tol, seed_contribution
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type, drop_n_jet_seeds, &
      & seed_normal_channel, jet_seed_index, &
      & drop_jet_seed_value, drop_jet_seed_grad, drop_jet_seed_hess
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent_point, &
      & drop_field_tangent_dir, drop_field_tangent_work_type
   use moist_cavity_drop_gaussian_scatter, only: scatter_iswig_block_indexed
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

   !> Initial entry capacity and bucket count of a per-thread accumulator
   !>
   !> Neither is a bound: both grow on demand, and both are trimmed at birth to
   !> a molecule that cannot fill them. The bucket count stays a power of two --
   !> which is what makes the mask in [[hess_sparse_entry]] a legal index.
   integer, parameter :: hess_init_ent = 256
   integer, parameter :: hess_init_tab = 1024

   !> Load factor above which the bucket table doubles, as `num/den`
   integer, parameter :: hess_load_num = 7, hess_load_den = 10

   !> FNV-1a (32-bit) constants of [[hess_pair_hash]]; every product is masked
   !> back to 32 bits, so no intermediate leaves the signed 64-bit range
   integer(int64), parameter :: fnv_offset = 2166136261_int64
   integer(int64), parameter :: fnv_prime = 16777619_int64
   integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)

   !> Sparse per-thread accumulator of the fixed-adjoint Hessian
   !>
   !> One thread's share of `(3, nsph, 3, nsph)`, held as the `(3, 3)` blocks of
   !> the atom pairs it actually reaches. A grid point reaches the level set's
   !> active atoms, its anchor's owner and the switching influence set of that
   !> anchor -- all of them local neighbourhoods -- so the pair count grows like
   !> `nsph * k` and not like `nsph^2`.
   !>
   !> **This container reproduces the dense one to the bit, and that is a
   !> property of how it is used, not of the arithmetic.** Two orders carry it:
   !>
   !>  1. *Accumulation order.* Every contribution is added in place, in grid
   !>     order, exactly where the dense slab's `+=` stood. One pair has one
   !>     entry -- the hash guarantees that -- so an element receives its
   !>     contributions in the same sequence as before. Growth preserves it:
   !>     the entry arrays are copied by `move_alloc` so indices never move, and
   !>     the bucket table is rebuilt from the entries rather than the reverse.
   !>
   !>  2. *Merge order.* [[hess_sparse_reduce]] is called in the same fixed
   !>     `1 .. nthreads` sequence the dense reduction used, and adds each
   !>     thread's block to the destination in one add per element -- the same
   !>     add the whole-array `hessian = hessian + hess_threads(..., ithread)`
   !>     performed. A pair no thread touched is simply skipped, where the dense
   !>     path added an exact zero to it.
   !>
   !> The bucket table is open addressing with linear probing, the idiom
   !> [[timer_type]] already uses in this codebase, with entry indices as values
   !> and `0` for an empty slot.
   type :: drop_hess_sparse_type
      !> Entries in use
      integer :: nent = 0
      !> Pairs that exist at all, `nsph^2` clamped to the integer range
      !>
      !> The entry arrays double, so a container that ends up holding every
      !> pair would otherwise carry up to twice the dense slab. Capping the
      !> capacity here bounds the overshoot at the point where the two meet,
      !> and is exact: a distinct pair beyond this one cannot be asked for
      integer :: maxent = 0
      !> Row and column atom of each entry
      integer, allocatable :: pair_i(:), pair_j(:)
      !> Accumulated `(3, 3)` block of each entry, `(row axis, column axis, entry)`
      real(wp), allocatable :: blocks(:, :, :)
      !> Open-addressing bucket table; values are entry indices, `0` is empty
      integer, allocatable :: htab(:)
   end type drop_hess_sparse_type

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
      type(drop_hess_sparse_type), allocatable :: hess_threads(:)
      !> Thread bookkeeping
      integer :: thread_slot, ithread
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure on its way to the latch
      type(error_type), allocatable :: worker_error

      !> Shared per-grid point sensitivity kernel state
      type(drop_seed_state_type) :: state
      !> Whether the shared prologue cleared the point for this traversal
      logical :: point_ok

      !> Grid, atom, axis, seed and active-slot indices
      integer :: igrid, ient, i, k, iaxis, ibasis
      integer :: n_active, ndir_atom, idir, dir_atom, dir_axis
      integer :: dir_loc, owner_loc, ia, ib
      integer, allocatable :: active_idx(:), dir_atoms(:)
      !> Accumulator entry of every `(row, column)` pair this point can reach,
      !> resolved once per point and read by both the anchor and field channels
      integer, allocatable :: pair_ent(:, :)

      !> Anchor and owner sphere
      real(wp) :: anchor(3)
      integer :: owner_idx
      !> Level-set jet at the projected point
      real(wp) :: lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :), lsf4_rrrr(:, :, :, :)
      !> Objective jet at the projected point
      real(wp) :: phi1_r(3)
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
      real(wp), allocatable :: swi_blk(:, :, :, :)
      integer, allocatable :: swi_idx(:), swi_ent(:, :)
      real(wp) :: swi_f0
      integer :: swi_n

      !> Effective primitive surface adjoints
      type(drop_surface_weights_type) :: eff
      !> Timer handle
      integer :: h_shess
      !> Rendered grid index of a degenerate point
      character(len=32) :: idx

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
      allocate (hess_threads(slots%nthreads))

      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& ient, i, k, iaxis, ibasis, n_active, ndir_atom, idir, dir_atom, dir_axis, &
      !$omp& dir_loc, owner_loc, ia, ib, active_idx, dir_atoms, pair_ent, state, point_ok, &
      !$omp& anchor, owner_idx, lsf1_r, lsf2_rr, lsf3_rrr, lsf4_rrrr, &
      !$omp& phi1_r, lambda_val, kkt_rhs, kkt_fac, &
      !$omp& seed_dlsf1, seed_dlsf2, seed_x, res_seed, dstate_seed, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_local, normal_grad, nwn, &
      !$omp& dv0, dv1, dv2, dv3, rhs_v, dr_v, dl_v, res_v, dstate_v, dinp_v, dres, &
      !$omp& dH_lag, dg_tot, dseed_x, dnormal_grad, dw_xyz_local, &
      !$omp& dw_lsf0, dw_lsf1, dw_lsf2, contribution, vdir, field_row, ft_work, &
      !$omp& iswig_work, swi_blk, swi_idx, swi_ent, swi_f0, swi_n, worker_error)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (lsf4_rrrr(3, 3, 3, 3), source=0.0_wp)
      allocate (dv3(3, 3, 3), source=0.0_wp)
      allocate (active_idx(self%nsph))
      allocate (dir_atoms(self%nsph))
      ! Grown from the direction set it indexes, never sized to `nsph`, for the
      ! same reason `swi_ent` below is: a molecule-sized map is the quadratic
      ! again. Both start at one slot and are grown by the first point.
      allocate (pair_ent(1, 1))
      call hess_sparse_new(hess_threads(thread_slot), self%nsph)
      allocate (vdir(3, self%nsph), source=0.0_wp)
      allocate (field_row(3, self%nsph), source=0.0_wp)
      allocate (res_seed(n_point_seeds))
      allocate (dstate_seed(n_point_seeds))
      call iswig_work%init(self%iswig)
      ! Grown on demand from `n_nb + 1` rather than sized to `nsph`: unlike the
      ! sparse rows of the gradient path a block is quadratic in the influence
      ! set, and a molecule-sized one would put back per thread exactly the
      ! quadratic the sparse accumulator exists to avoid.
      allocate (swi_blk(3, 1, 3, 1), swi_idx(1), swi_ent(1, 1))

      !$omp do schedule(static, 8)
      do igrid = 1, self%ngrid
         !* -------------------------- Shared point prologue -------------------------- *!
         ! Point, jets, seed state and the solved jet and anchor seeds. Every
         ! failure path -- a latch already set, a refusing level set, a
         ! degenerate state, a singular bordered system -- has recorded itself.
         ! `lsf4_rrrr` is the one extra jet order this half needs: the field
         ! tangent reads it, and CFC asks for the highest *total* order.
         call drop_point_prologue(self, slots, thread_slot, igrid, eff%have_wk, &
                                  "get_surface_hessian_fixed_drop", abort, anchor, &
                                  owner_idx, lambda_val, lsf1_r, lsf2_rr, lsf3_rrr, &
                                  phi1_r, state, kkt_fac, point_ok, &
                                  lsf4_rrrr=lsf4_rrrr, kkt_rhs=kkt_rhs)
         if (.not. point_ok) cycle

         ! Outward-normal channel, as in the gradient path. The direction loop
         ! below needs `normal_grad` and `nwn` again to build their own
         ! tangents, so they are asked of the fold rather than recomputed here:
         ! a second copy of one floating-point chain is free to contract
         ! differently, and [[seed_normal_channel]] owns this one. Both come
         ! back zero when the channel is inactive.
         w_lsf0_pt = 0.0_wp
         w_lsf1_pt = 0.0_wp
         w_lsf2_pt = 0.0_wp
         call seed_normal_channel(state, eff, igrid, lsf2_rr, w_lsf1_pt, w_xyz_local, &
                                  normal_grad_pt=normal_grad, nwn_pt=nwn)

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
         end do

         ! The two bounds differ on purpose. The responses above are needed for
         ! all 16 seeds, because the direction loop reads every `res_seed` and
         ! `dstate_seed` again. The contraction below is a *jet* one: it builds
         ! the level-set adjoint weights, and only the 13 jet slots have a
         ! weight to land in. An anchor seed's contraction belongs to the
         ! owner's gradient row instead, which is the direction loop's
         ! `ibasis > drop_n_jet_seeds` branch -- and is why
         ! [[scatter_jet_weight]] deliberately has no anchor case.
         do ibasis = 1, drop_n_jet_seeds
            contribution = seed_contribution(eff, igrid, w_xyz_local, &
                                             seed_x(1:3, ibasis), res_seed(ibasis), phi1_r)
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
         owner_loc = ndir_atom
         do i = 1, ndir_atom
            if (dir_atoms(i) == owner_idx) then
               owner_loc = i
               exit
            end if
         end do

         ! Rows are the active atoms and the owner, columns the same set, so the
         ! whole square is reached and all of it is resolved here -- once per
         ! point rather than once per direction. Creating an entry commits
         ! nothing: an untouched one holds the exact zero the dense slab held.
         if (size(pair_ent, 1) < ndir_atom) then
            deallocate (pair_ent)
            allocate (pair_ent(ndir_atom, ndir_atom))
         end if
         do ib = 1, ndir_atom
            do ia = 1, ndir_atom
               call hess_sparse_entry(hess_threads(thread_slot), dir_atoms(ia), &
                                      dir_atoms(ib), pair_ent(ia, ib))
            end do
         end do

         ! `d^4S/(dr^3 dR_A)` is a function of the prepared point alone and is
         ! the most expensive accessor the field tangent calls, so its fill is
         ! the point half of [[drop_field_tangent]] and runs once here rather
         ! than `3 * ndir_atom` times inside the loop below. It has to run on
         ! every point that reaches the loop, not conditionally: the buffer is
         ! reused across points and a skipped fill would leave the previous
         ! point's tensor in it.
         call drop_field_tangent_point(slots%lsf(thread_slot)%lsf, ft_work)

         do idir = 1, 3*ndir_atom
            dir_loc = (idir - 1)/3 + 1
            dir_atom = dir_atoms(dir_loc)
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
            call kkt_fac%solve(rhs_v, "get_surface_hessian_fixed_drop", &
                               worker_error, igrid)
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
            call kkt_fac%solve_tangent(dH_lag, dg_tot, seed_x, dseed_x, &
                                       "get_surface_hessian_fixed_drop", &
                                       worker_error, igrid)
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
                  hess_threads(thread_slot)%blocks(iaxis, dir_axis, &
                                                   pair_ent(owner_loc, dir_loc)) = &
                     hess_threads(thread_slot)%blocks(iaxis, dir_axis, &
                                                      pair_ent(owner_loc, dir_loc)) &
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
               call drop_field_tangent_dir(slots%lsf(thread_slot)%lsf, &
                                           w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, &
                                           dw_lsf0, dw_lsf1, dw_lsf2, dr_v, vdir, &
                                           ft_work, field_row)
               do i = 1, n_active
                  ient = pair_ent(i, dir_loc)
                  hess_threads(thread_slot)%blocks(:, dir_axis, ient) = &
                     hess_threads(thread_slot)%blocks(:, dir_axis, ient) &
                     + field_row(:, i)
               end do
            end if

            vdir(dir_axis, dir_atom) = 0.0_wp
         end do
         if (abort%requested) cycle

         !* ------------------------- iSwig switching channel ------------------------- *!
         ! `f_i` depends on the nuclear geometry alone and its adjoint is fixed,
         ! so the whole channel is one block over the influence set. The width
         ! rows `mix` and `d2xi` are not asked for: the width the switching
         ! factor is evaluated at is the *anchor* width, built from the rigid
         ! sphere's Lebedev weight, and so is nuclear-geometry independent --
         ! which is why the gradient path drops `swi_dxi` as well. Left absent
         ! they cost nothing: [[iswig_swi_f2_rArB_block]] skips their pass.
         if (abs(eff%w_f(igrid)) > seed_weight_tol) then
            call self%iswig%swi_collect(anchor, owner_idx, self%anchor_xi0(igrid), &
                                        swi_f0, iswig_work)
            if (size(swi_idx) < iswig_work%n_nb + 1) then
               deallocate (swi_blk, swi_idx, swi_ent)
               allocate (swi_blk(3, iswig_work%n_nb + 1, 3, iswig_work%n_nb + 1))
               allocate (swi_idx(iswig_work%n_nb + 1))
               allocate (swi_ent(iswig_work%n_nb + 1, iswig_work%n_nb + 1))
            end if
            call self%iswig%swi2_rArB_block(iswig_work, swi_n, swi_idx, swi_blk)
            ! The influence set is not the direction set, so its pairs are
            ! resolved on their own. Any pair the two share resolves to the one
            ! entry the field channel already wrote, and lands after it.
            do ib = 1, swi_n
               do ia = 1, swi_n
                  call hess_sparse_entry(hess_threads(thread_slot), swi_idx(ia), &
                                         swi_idx(ib), swi_ent(ia, ib))
               end do
            end do
            call scatter_iswig_block_indexed(swi_n, swi_ent, swi_blk, eff%w_f(igrid), &
                                             hess_threads(thread_slot)%blocks)
         end if

      end do
      !$omp end do

      deallocate (lsf3_rrr, lsf4_rrrr, dv3, active_idx, dir_atoms, vdir, field_row)
      deallocate (res_seed, dstate_seed, swi_blk, swi_idx, swi_ent, pair_ent)
      call iswig_work%destroy()
      !$omp end parallel

      if (abort%requested) then
         ! An LSF failure or a KKT failure arrives as a ready-made error; a
         ! kernel degeneracy arrives as a status code that needs this routine's
         ! name to become a diagnostic.
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            write (idx, "(i0)") abort%igrid
            call fatal_error(error, "get_surface_hessian_fixed_drop: "// &
                             seed_status_message(abort%status)// &
                             " at grid point "//trim(idx))
         end if
         call self%ctx%timer%stop(h_shess)
         return
      end if

      ! Deterministic reduction: fixed thread order, independent of scheduling
      do ithread = 1, slots%nthreads
         call hess_sparse_reduce(hess_threads(ithread), hessian)
         call hess_sparse_destroy(hess_threads(ithread))
      end do

      call self%ctx%timer%stop(h_shess)

   end subroutine get_surface_hessian_fixed_drop


   !* ================================================================================= *!
   !*                       Sparse per-thread Hessian accumulator                       *!
   !* ================================================================================= *!

   !> Give one thread an empty accumulator
   !>
   !> Called from inside the parallel region on the thread that will own it, so
   !> the storage is first touched by the thread that writes it.
   !>
   !> @param[inout] self Accumulator of one thread
   !> @param[in]    nsph Atoms the pair indices range over
   subroutine hess_sparse_new(self, nsph)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self
      !> Atoms the pair indices range over
      integer, intent(in) :: nsph

      integer :: cap, ntab

      self%nent = 0
      self%maxent = int(min(int(nsph, int64)**2, int(huge(self%maxent), int64)))

      cap = max(1, min(hess_init_ent, self%maxent))
      allocate (self%pair_i(cap), self%pair_j(cap))
      allocate (self%blocks(ndim, ndim, cap))

      ! Halve the default table while the one below it would still hold `cap`
      ! entries under the load factor, so a molecule with a handful of pairs
      ! does not carry a kilobyte of empty buckets per thread
      ntab = hess_init_tab
      do while (ntab > 4 .and. cap*hess_load_den < (ntab/2)*hess_load_num)
         ntab = ntab/2
      end do
      allocate (self%htab(ntab), source=0)

   end subroutine hess_sparse_new

   !> Release one thread's accumulator
   !>
   !> @param[inout] self Accumulator of one thread
   subroutine hess_sparse_destroy(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      self%nent = 0
      self%maxent = 0
      if (allocated(self%pair_i)) deallocate (self%pair_i)
      if (allocated(self%pair_j)) deallocate (self%pair_j)
      if (allocated(self%blocks)) deallocate (self%blocks)
      if (allocated(self%htab)) deallocate (self%htab)

   end subroutine hess_sparse_destroy

   !> FNV-1a hash of an atom pair, masked to 32 bits
   !>
   !> Fed the two indices as whole words rather than byte by byte. Multiplying
   !> by an odd constant is a bijection modulo any power of two, so the low bits
   !> this hash is indexed by are a permutation of a contiguous run of atom ids
   !> -- which is exactly the run one grid point's neighbourhood hands it.
   !>
   !> @param[in] iatom Row atom
   !> @param[in] jatom Column atom
   !> @return          Hash value in `[0, 2^32)`
   pure function hess_pair_hash(iatom, jatom) result(h)
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Hash value
      integer(int64) :: h

      h = ieor(fnv_offset, iand(int(iatom, int64), mask32))
      h = iand(h*fnv_prime, mask32)
      h = ieor(h, iand(int(jatom, int64), mask32))
      h = iand(h*fnv_prime, mask32)

   end function hess_pair_hash

   !> Index of the `(3, 3)` block of one atom pair, creating it if it is new
   !>
   !> A new entry is zeroed here and nowhere else, so an entry that is created
   !> and then never written contributes an exact zero, and one that is created
   !> twice is impossible. Both are what makes the merge a single add per
   !> element per thread; see [[drop_hess_sparse_type]].
   !>
   !> @param[inout] self  Accumulator of one thread
   !> @param[in]    iatom Row atom
   !> @param[in]    jatom Column atom
   !> @param[out]   ient  Entry index of the pair
   subroutine hess_sparse_entry(self, iatom, jatom, ient)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self
      !> Row and column atom
      integer, intent(in) :: iatom, jatom
      !> Entry index of the pair
      integer, intent(out) :: ient

      integer :: slot

      slot = int(iand(hess_pair_hash(iatom, jatom), &
                      int(size(self%htab) - 1, int64)), kind(slot)) + 1
      do
         ient = self%htab(slot)
         if (ient == 0) exit
         if (self%pair_i(ient) == iatom .and. self%pair_j(ient) == jatom) return
         slot = slot + 1
         if (slot > size(self%htab)) slot = 1
      end do

      call hess_grow_entries(self)
      self%nent = self%nent + 1
      ient = self%nent
      self%pair_i(ient) = iatom
      self%pair_j(ient) = jatom
      self%blocks(:, :, ient) = 0.0_wp

      ! The probe above found this slot on the *current* table, so it is only
      ! valid while that table stands; a growth reinserts every entry instead.
      if (self%nent*hess_load_den >= size(self%htab)*hess_load_num) then
         call hess_grow_table(self)
      else
         self%htab(slot) = ient
      end if

   end subroutine hess_sparse_entry

   !> Double the entry arrays if the next entry would not fit
   !>
   !> `move_alloc` on a copy rather than a reallocation in place: entry indices
   !> are handed out to the caller and must not move, and the accumulated blocks
   !> are copied verbatim, never re-derived.
   !>
   !> @param[inout] self Accumulator of one thread
   subroutine hess_grow_entries(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      integer, allocatable :: new_i(:), new_j(:)
      real(wp), allocatable :: new_blocks(:, :, :)
      integer :: cap, new_cap

      cap = size(self%pair_i)
      if (self%nent < cap) return

      new_cap = min(2*cap, self%maxent)
      if (new_cap <= cap) then
         ! Unreachable: `maxent` counts every pair that exists, and each is
         ! entered once, so a full container is never asked for another
         error stop "moist DROP Hessian: the sparse accumulator was asked for more "// &
            "pairs than the molecule has"
      end if

      allocate (new_i(new_cap))
      new_i(1:cap) = self%pair_i(1:cap)
      call move_alloc(new_i, self%pair_i)

      allocate (new_j(new_cap))
      new_j(1:cap) = self%pair_j(1:cap)
      call move_alloc(new_j, self%pair_j)

      allocate (new_blocks(ndim, ndim, new_cap))
      new_blocks(:, :, 1:cap) = self%blocks(:, :, 1:cap)
      call move_alloc(new_blocks, self%blocks)

   end subroutine hess_grow_entries

   !> Double the bucket table and reinsert every entry
   !>
   !> Touches no block and no entry index -- only the pair-to-entry lookup is
   !> rebuilt, so neither the accumulated values nor the order they were
   !> accumulated in can move.
   !>
   !> @param[inout] self Accumulator of one thread
   subroutine hess_grow_table(self)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(inout) :: self

      integer, allocatable :: new_tab(:)
      integer :: ient, slot

      allocate (new_tab(size(self%htab)*2), source=0)
      call move_alloc(new_tab, self%htab)

      do ient = 1, self%nent
         slot = int(iand(hess_pair_hash(self%pair_i(ient), self%pair_j(ient)), &
                         int(size(self%htab) - 1, int64)), kind(slot)) + 1
         do
            if (self%htab(slot) == 0) then
               self%htab(slot) = ient
               exit
            end if
            slot = slot + 1
            if (slot > size(self%htab)) slot = 1
         end do
      end do

   end subroutine hess_grow_table

   !> Add one thread's accumulator to the caller's Hessian
   !>
   !> One add per element, on the elements this thread reached. Called in the
   !> fixed `1 .. nthreads` order, which is the second half of the
   !> bit-reproducibility argument in [[drop_hess_sparse_type]].
   !>
   !> @param[in]    self    Accumulator of one thread
   !> @param[inout] hessian Nuclear-Hessian accumulator (3, nsph, 3, nsph)
   subroutine hess_sparse_reduce(self, hessian)
      !> Accumulator of one thread
      type(drop_hess_sparse_type), intent(in) :: self
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)

      integer :: ient, iatom, jatom

      do ient = 1, self%nent
         iatom = self%pair_i(ient)
         jatom = self%pair_j(ient)
         hessian(:, iatom, :, jatom) = hessian(:, iatom, :, jatom) &
                                       + self%blocks(:, :, ient)
      end do

   end subroutine hess_sparse_reduce

   !* ================================================================================= *!
   !*                              Scope of the fixed half                              *!
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
   !> Which jet slot carries which direction is not decided here:
   !> [[jet_seed_index]] owns that layout and this routine only consumes it, so
   !> a change to the jet basis -- symmetrising the Hessian seeds, say -- does
   !> not have to be repeated in this file.
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

      !> Seed index, its kind and its Cartesian indices
      integer :: ibasis, slot_kind, iaxis, jaxis

      dlsf1_r = 0.0_wp
      dlsf2_rr = 0.0_wp
      x = 0.0_wp

      do ibasis = 1, drop_n_jet_seeds
         call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
         if (slot_kind == drop_jet_seed_hess) then
            dlsf2_rr(iaxis, jaxis, ibasis) = 1.0_wp
         else
            if (slot_kind == drop_jet_seed_grad) dlsf1_r(iaxis, ibasis) = 1.0_wp
            x(:, ibasis) = kkt(1:4, ibasis)
         end if
      end do

      do iaxis = 1, ndim
         x(:, drop_n_jet_seeds + iaxis) = kkt(1:4, 4 + iaxis)
      end do
   end subroutine fill_seed_basis

   !> Place one jet seed's contribution in the level-set adjoint weights
   !>
   !> The same layout [[fill_seed_basis]] writes, read back through the same
   !> classifier: [[jet_seed_index]] says which of the value, gradient and
   !> Hessian channels a slot belongs to, and this routine only places the
   !> number. Anchor slots are not accepted -- they classify as
   !> `drop_jet_seed_none` and land nowhere -- because their contribution
   !> belongs to the owner's gradient row, not to a weight.
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

      !> Kind of the slot and its Cartesian indices
      integer :: slot_kind, iaxis, jaxis

      call jet_seed_index(ibasis, slot_kind, iaxis, jaxis)
      select case (slot_kind)
      case (drop_jet_seed_value)
         w0 = w0 + contrib
      case (drop_jet_seed_grad)
         w1(iaxis) = w1(iaxis) + contrib
      case (drop_jet_seed_hess)
         w2(iaxis, jaxis) = w2(iaxis, jaxis) + contrib
      end select
   end subroutine scatter_jet_weight

   !* ================================================================================= *!
   !*                        Adjoint contraction of one seed                            *!
   !* ================================================================================= *!

   !> Directional derivative of [[seed_contribution]]
   !>
   !> Term by term the product rule applied to that contraction, with the
   !> surface adjoints themselves held fixed -- which is the whole premise of
   !> this half. The position adjoint still moves, because the normal fold
   !> inside it is built from the level-set gradient at the projected point.
   !>
   !> This one stays local, unlike its primal: this submodule is its only
   !> caller, so there is no second copy of the chain to share with. The branch
   !> term is absent for the reason the module header gives -- a grid carrying
   !> a multi-branch anchor group is rejected by [[check_frozen_weights]] --
   !> and the switching term for the reason [[seed_contribution]] gives.
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
