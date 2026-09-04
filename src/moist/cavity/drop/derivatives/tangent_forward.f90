!> Pass 1 of the DROP Hessian: forward tangent of the surface map
!>
!> Propagates a set of nuclear directions `v` through the same per-point map
!> the reverse path differentiates, and returns the directional derivatives of
!> the four grid scalars pass 2 consumes:
!>
!>     d_a, d_wleb, d_xi0, d_wbranch    (ngrid, ndir)
!>
!> This is the forward-mode dual of [[get_surface_gradient_drop]]: same thread
!> setup, same per-point `prepare` / jet / [[build_seed_state]] / KKT
!> factorization, same abort latch. What differs is the seeding. The reverse
!> path pushes 13 jet basis directions plus 3 anchor directions through
!> [[apply_seed]] and contracts the answers with a surface adjoint; this one
!> pushes exactly one seed per nuclear direction -- the direction's own image
!> in the level-set jet -- and keeps the response.
!>
!> Unlike the reverse path there is no cross-point reduction to make
!> deterministic: point `i` writes rows `(i, :)` of four arrays it shares with
!> nobody, so the grid loop needs neither per-thread buffers nor a fixed-order
!> sum. The one cross-point coupling in the scheme is the branch softmax, and
!> it runs outside the parallel region; see below.
!>
!> ## The three stages
!>
!>  1. **Grid loop (parallel).** Per point: the level-set jet and its
!>     directional nuclear tangents (`tangent_f0/f1_r/f2_rr` -- the LSF's own
!>     contracted accessors, so the active-slot index space never leaves the
!>     level set), the bordered KKT solve for `(dr, dlambda)` batched over all
!>     directions, [[apply_seed]] for the base Lebedev-weight motion, the
!>     sparse iSwiG rows for `d(f)`, and the branch objective's tangent
!>     `d(Phi)`.
!>  2. **Branch softmax (serial).** One [[branch_weight_type:weights_grad]]
!>     call per contiguous anchor group, with `nparam = ndir`, giving
!>     `d_wbranch` for every direction of every branch at once.
!>  3. **Assembly.** The branch motion is added to `d_wleb`, and `d_xi0` and
!>     `d_a` follow from it by the two identities below.
!>
!> ## The `wbranch` trap
!>
!> [[apply_seed]] does **not** return a complete `d(wleb)`. Its last line is
!>
!>     res%dwleb = state%wbranch * state%wleb_prune_factor * dw_pre
!>
!> with `wbranch` held fixed: in first-order reverse mode the branch weight's
!> own motion is not a term of this product, it is handled separately through
!> `branch_phi_adj` (see [[compute_branch_phi_adj]]). A forward tangent has no
!> such second channel and must put the term back,
!>
!>     d_wleb = res%dwleb + d_wbranch * wleb / wbranch
!>
!> and carry it on into `d_xi0` and `d_a`. Dropping it leaves every downstream
!> tangent wrong on multi-branch anchors only -- invisible to any test of
!> [[apply_seed]] itself, and invisible to a fixture that never branches. The
!> second-order kernel path anticipates the same correction: its
!> [[apply_seed_tangent]] takes `dinp_v%dwbranch` as an explicit input.
!>
!> ## Parallelisation of the branch stage
!>
!> The softmax couples every branch of one anchor group through its
!> normalisation, so a group split across two threads would silently corrupt
!> its reduction. Groups are runs of equal `anchor_id` -- contiguity is
!> guaranteed by the stable `counting_argsort` in `projection.f90` -- but the
!> grid loop above is chunked by grid point and knows nothing about them.
!> Stage 2 therefore runs **serially over groups**, outside the parallel
!> region, with all `ndir` directions batched into the one `weights_grad` call
!> the group needs. That is the cheapest of the three admissible choices here:
!> the walk touches only multi-branch groups, which are a small minority of the
!> grid, and it needs no group index to have been built.
!>
!> ## Two identities
!>
!> Both are exact, and both are used rather than recomputed:
!>
!>   * `xi0 = swx/(R sqrt(wleb))` with the *final* `wleb`
!>     (`projection.f90`, `compute_gaussians`), so
!>     `d_xi0 = -0.5 * xi0 * d_wleb / wleb`. [[apply_seed]] applies the same
!>     identity to its own partial `dwleb`; `res%dxi` is therefore the
!>     branch-frozen width tangent and is deliberately discarded here in favour
!>     of one application to the completed `d_wleb`.
!>   * `a = R^2 f wleb` (`properties.f90`, `compute_area_volume`), so
!>     `d_a = R^2 (wleb d_f + f d_wleb)`, where `f` is the **iSwiG** switching
!>     factor `self%f` -- not [[apply_seed]]'s `res%dw_f`, which is the
!>     `w_f = f_crit * f_foc` product that rides inside `wleb`.
!>
!> ## Redundancy of the outputs
!>
!> `d_wbranch` is reported separately *and* is already folded into `d_wleb`.
!> That is not a duplication to be optimized away: pass 2's
!> [[branch_point_adjoint]] differentiates `wleb/wbranch` and needs both halves
!> independently, and it is only consistent if `d_wleb` is the complete tangent.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_tangent_forward
!$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_seed_result_type, build_seed_state, apply_seed, seed_state_ok, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type, degenerate_point_error
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !> Forward tangent of the DROP surface map along a set of nuclear directions
   !>
   !> Every output is `(ngrid, ndir)` and is written in full: column `idir`
   !> holds the directional derivative along `dirs(:, :, idir)`.
   !>
   !> @param[in]  self      DROP cavity instance (must hold a projected grid)
   !> @param[in]  dirs      Nuclear directions `(3, nsph, ndir)`
   !> @param[out] d_a       Tangent of the area element `(ngrid, ndir)`
   !> @param[out] d_wleb    Tangent of the Lebedev weight `(ngrid, ndir)`
   !> @param[out] d_xi0     Tangent of the Gaussian width `(ngrid, ndir)`
   !> @param[out] d_wbranch Tangent of the branch weight `(ngrid, ndir)`
   !> @param[out] error     Error object, allocated on failure
   module subroutine get_surface_tangent_drop(self, dirs, d_a, d_wleb, d_xi0, &
                                              d_wbranch, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Directional tangents of the four grid scalars
      real(wp), intent(out) :: d_a(:, :), d_wleb(:, :), d_xi0(:, :), d_wbranch(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Thread bookkeeping
      integer :: thread_slot
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure on its way to the latch
      type(error_type), allocatable :: worker_error

      !> Shared per-grid point sensitivity kernel state and its response
      type(drop_seed_state_type) :: state
      type(drop_seed_result_type) :: res
      !> Degeneracy status
      integer :: status

      !> Grid, direction, sphere and Cartesian indices
      integer :: igrid, idir, ndir, jj, knb

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

      !> Directional nuclear tangents of the level-set jet at the *fixed* point
      real(wp) :: dlsf0
      real(wp), allocatable :: dlsf1_r(:, :), dlsf2_rr(:, :, :)
      !> Bordered KKT right-hand sides, one column per direction
      real(wp), allocatable :: kkt_rhs(:, :)
      !> Factorization reused by every direction at this grid point
      type(drop_kkt_factor_type) :: kkt_fac
      !> Induced motion of the projected point and of the multiplier
      real(wp) :: dr(3), dlambda

      !> iSwiG neighbour cache and the sparse switching rows it feeds
      type(iswig_workspace_type) :: iswig_work
      real(wp), allocatable :: swi_rows(:, :)
      real(wp) :: swi_owner_row(3), swi_f0, swi_dxi, df_dir

      !> Tangent of the branch objective, `d(Phi)` per point and direction
      real(wp), allocatable :: dphi(:, :)
      !> Softmax scratch of the branch stage
      real(wp), allocatable :: branch_phi(:), branch_dphi(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)
      !> Branch group bookkeeping
      integer :: igroup_start, igroup_end, group_size, m_branch, im_grid, nbranch_max

      !> Assembly scalars
      real(wp) :: wleb_i, wbranch_i, dwleb_i, r_own
      !> Timer handle
      integer :: h_stan

      !* ------------------------------- Shape guards --------------------------------- *!
      if (size(dirs, 1) /= ndim .or. size(dirs, 2) /= self%nsph) then
         call fatal_error(error, "get_surface_tangent_drop: dirs must be (3, nsph, ndir)")
         return
      end if
      ndir = size(dirs, 3)
      if (ndir <= 0) then
         call fatal_error(error, "get_surface_tangent_drop: no direction supplied")
         return
      end if
      if (size(d_a, 1) /= self%ngrid .or. size(d_a, 2) /= ndir .or. &
          size(d_wleb, 1) /= self%ngrid .or. size(d_wleb, 2) /= ndir .or. &
          size(d_xi0, 1) /= self%ngrid .or. size(d_xi0, 2) /= ndir .or. &
          size(d_wbranch, 1) /= self%ngrid .or. size(d_wbranch, 2) /= ndir) then
         call fatal_error(error, "get_surface_tangent_drop: every output must be"// &
                          " (ngrid, ndir)")
         return
      end if

      d_a = 0.0_wp
      d_wleb = 0.0_wp
      d_xi0 = 0.0_wp
      d_wbranch = 0.0_wp
      if (self%ngrid <= 0) return

      h_stan = self%ctx%timer%resolve("Surface tangent", self%ctx%timer%current(), &
                                      cat_gradient)
      call self%ctx%timer%start(h_stan)

      !* -------------------------------- Thread setup -------------------------------- *!
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      allocate (dphi(self%ngrid, ndir), source=0.0_wp)

      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& idir, jj, knb, state, res, status, &
      !$omp& point, anchor, owner_idx, lsf0, lsf1_r, lsf2_rr, lsf3_rrr, &
      !$omp& phi0, phi1_r, phi2_rr, lambda_val, &
      !$omp& dlsf0, dlsf1_r, dlsf2_rr, kkt_rhs, kkt_fac, dr, dlambda, &
      !$omp& iswig_work, swi_rows, swi_owner_row, swi_f0, swi_dxi, df_dir, &
      !$omp& worker_error)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (dlsf1_r(3, ndir), source=0.0_wp)
      allocate (dlsf2_rr(3, 3, ndir), source=0.0_wp)
      allocate (kkt_rhs(4, ndir), source=0.0_wp)
      call iswig_work%init(self%iswig)
      ! Sized to `nsph` rather than to the workspace capacity: `n_nb` is bounded
      ! by the atom count on either traversal, so this stays valid even if the
      ! workspace has to grow itself.
      allocate (swi_rows(3, self%nsph))

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
         call slots%phi(thread_slot)%f012_r(point, anchor, owner_idx, phi0, phi1_r, phi2_rr)

         state%lsf1_r = lsf1_r
         state%lsf2_rr = lsf2_rr
         state%lsf3_rrr = lsf3_rrr
         state%lambda_val = lambda_val
         call fill_seed_state(self, igrid, .false., state)

         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)
         if (status /= seed_state_ok) then
            call abort%latch_status(status, igrid)
            cycle
         end if

         !* --------------- Directional nuclear tangents of the jet ---------------- *!
         ! The level set contracts its own nuclear index: `tangent_*` returns
         ! `sum_B v_B . d(jet)/dR_B` at the fixed point, so the active-slot
         ! index space never has to be reconciled with the atom index space of
         ! `dirs` out here -- which is exactly where the two are easy to confuse.
         do idir = 1, ndir
            call slots%lsf(thread_slot)%lsf%tangent_f0(dirs(:, :, idir), dlsf0)
            call slots%lsf(thread_slot)%lsf%tangent_f1_r(dirs(:, :, idir), dlsf1_r(:, idir))
            call slots%lsf(thread_slot)%lsf%tangent_f2_rr(dirs(:, :, idir), dlsf2_rr(:, :, idir))

            ! Bordered right-hand side of the direction, with
            ! `d^2 phi/(dr dR_owner) = -alpha*I` (`objective_phi.f90`, `f2_r_rA`)
            ! and the anchor riding its owner rigidly.
            kkt_rhs(1:3, idir) = self%param%phi_alpha*dirs(:, owner_idx, idir) &
                                 + lambda_val*dlsf1_r(:, idir)
            kkt_rhs(4, idir) = -dlsf0
         end do

         !* ------------------------ Bordered KKT sensitivities ----------------------- *!
         ! One factorization, one batched solve: every direction shares the 4x4
         ! matrix of this grid point.
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

         !* ---------------------- Base Lebedev-weight response ----------------------- *!
         do idir = 1, ndir
            dr = kkt_rhs(1:3, idir)
            dlambda = kkt_rhs(4, idir)

            call apply_seed(state, dlsf1_r(:, idir), dlsf2_rr(:, :, idir), dr, dlambda, res)

            ! Branch-frozen half of `d(wleb)`; stage 3 completes it. `res%dxi`
            ! is the width tangent of exactly this incomplete weight and is not
            ! read at all.
            d_wleb(igrid, idir) = res%dwleb

            ! Tangent of the branch objective `Phi = 0.5 alpha |r* - anchor|^2`
            ! along the direction: the projected point moves by `dr`, the anchor
            ! rigidly with its owner.
            dphi(igrid, idir) = dot_product(phi1_r, dr - dirs(:, owner_idx, idir))
         end do

         !* ------------------------- iSwiG switching channel ------------------------- *!
         ! `f` is evaluated at the anchor with the anchor width, and
         ! `anchor_xi0` depends on the owner radius and the raw Lebedev weight
         ! alone, so it carries no nuclear tangent and `swi_dxi` is unused.
         ! Parked in `d_a` until stage 3 turns it into the area tangent.
         call self%iswig%swi_collect(anchor, owner_idx, self%anchor_xi0(igrid), &
                                     swi_f0, iswig_work)
         call self%iswig%swi1_rA_sparse(iswig_work, swi_rows, swi_owner_row, swi_dxi)
         do idir = 1, ndir
            df_dir = dot_product(swi_owner_row, dirs(:, owner_idx, idir))
            do jj = 1, iswig_work%n_nb
               knb = iswig_work%idx(jj)
               df_dir = df_dir + dot_product(swi_rows(:, jj), dirs(:, knb, idir))
            end do
            d_a(igrid, idir) = df_dir
         end do

      end do
      !$omp end do

      deallocate (lsf3_rrr, dlsf1_r, dlsf2_rr, kkt_rhs, swi_rows)
      call iswig_work%destroy()
      !$omp end parallel

      if (abort%requested) then
         ! An LSF failure or a KKT failure arrives as a ready-made error; a
         ! kernel degeneracy arrives as a status code that needs this routine's
         ! name to become a diagnostic.
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("get_surface_tangent_drop", abort%status, &
                                        abort%igrid, error)
         end if
         call self%ctx%timer%stop(h_stan)
         return
      end if

      !* ------------------------- Branch softmax (serial) ---------------------------- *!
      ! Serial over contiguous anchor groups, with every direction batched into
      ! the one `weights_grad` call per group; see the module header for why a
      ! group must never be split. Points outside a multi-branch group keep
      ! `d_wbranch = 0`, which is exact: their `wbranch` is the constant one.
      call branch_stage()

      !* -------------------------------- Assembly ------------------------------------ *!
      do idir = 1, ndir
         do igrid = 1, self%ngrid
            wleb_i = self%wleb(igrid)
            wbranch_i = self%wbranch(igrid)
            r_own = self%radii(self%owner(igrid))

            ! The wbranch trap: `apply_seed` froze the branch weight, so its
            ! motion is added back here. `wleb/wbranch` is the pre-branch weight
            ! the softmax multiplies, formed the way `forward.f90`'s branch
            ! post-pass and `compute_branch_phi_adj` both form it.
            dwleb_i = d_wleb(igrid, idir)
            if (wbranch_i > tiny(1.0_wp)) then
               dwleb_i = dwleb_i + (wleb_i/wbranch_i)*d_wbranch(igrid, idir)
            end if
            d_wleb(igrid, idir) = dwleb_i

            ! xi0 = swx/(R sqrt(wleb)); same guard as `apply_seed` and `iswig_xi0`
            if (wleb_i > seed_weight_tol) then
               d_xi0(igrid, idir) = -0.5_wp*self%xi0(igrid)*dwleb_i/wleb_i
            else
               d_xi0(igrid, idir) = 0.0_wp
            end if

            ! a = R^2 f wleb, with `d_a` still holding the parked `d(f)`
            d_a(igrid, idir) = r_own*r_own &
                               *(wleb_i*d_a(igrid, idir) + self%f(igrid)*dwleb_i)
         end do
      end do

      call self%ctx%timer%stop(h_stan)

   contains

      !> Differentiate the branch softmax over every contiguous anchor group
      !>
      !> Mirrors the group walk of [[compute_branch_phi_adj]] and of
      !> `forward.f90`'s branch post-pass: runs of equal `anchor_id` starting at
      !> a point with `branch_count > 1`. The softmax primitive takes its
      !> derivatives in `(nparam, nbranch)` layout, so passing `nparam = ndir`
      !> yields every direction of the group from one call.
      subroutine branch_stage()

         if (.not. allocated(self%branch_count) .or. .not. allocated(self%anchor_id)) return
         if (.not. any(self%branch_count(1:self%ngrid) > 1)) return

         nbranch_max = maxval(self%branch_count(1:self%ngrid))
         allocate (branch_phi(nbranch_max), source=0.0_wp)
         allocate (branch_dphi(ndir, nbranch_max), source=0.0_wp)
         allocate (branch_weights(nbranch_max), source=0.0_wp)
         allocate (branch_dweights(ndir, nbranch_max), source=0.0_wp)

         igroup_start = 1
         do while (igroup_start <= self%ngrid)
            if (self%branch_count(igroup_start) <= 1) then
               igroup_start = igroup_start + 1
               cycle
            end if

            ! Extend the group while anchor_id stays the same
            igroup_end = igroup_start
            do while (igroup_end < self%ngrid)
               if (self%anchor_id(igroup_end + 1) /= self%anchor_id(igroup_start)) exit
               igroup_end = igroup_end + 1
            end do
            group_size = igroup_end - igroup_start + 1

            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               branch_phi(m_branch) = self%phi0(im_grid)
               branch_dphi(:, m_branch) = dphi(im_grid, :)
            end do

            call self%branch_weight%weights_grad( &
               branch_phi(1:group_size), branch_dphi(:, 1:group_size), &
               weights=branch_weights(1:group_size), &
               dweights=branch_dweights(:, 1:group_size))

            do m_branch = 1, group_size
               im_grid = igroup_start + m_branch - 1
               d_wbranch(im_grid, :) = branch_dweights(:, m_branch)
            end do

            igroup_start = igroup_end + 1
         end do

         deallocate (branch_phi, branch_dphi, branch_weights, branch_dweights)

      end subroutine branch_stage

   end subroutine get_surface_tangent_drop

end submodule moist_cavity_drop_derivatives_tangent_forward
