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
!>  1. **Grid loop (parallel).** Per point: the level-set jet, its directional
!>     nuclear tangents along every direction -- from the mixed nuclear tensors
!>     materialised once ([[drop_field_jet_point]]) and contracted per direction
!>     ([[drop_field_jet_tangent]]), or, when the caller asks for it because it
!>     has a few directions, from the level set's contracted accessor
!>     `tangent_jet` per direction -- the bordered KKT solve for `(dr, dlambda)`
!>     batched over all
!>     directions, [[apply_seed]] for the base Lebedev-weight motion, the
!>     sparse iSwiG rows for `d(f)`, and the branch objective's tangent
!>     `d(Phi)`. The tensors are active-slot indexed, so each direction is
!>     gathered onto the point's active atoms first; that gather is the only
!>     place the two index spaces meet in this routine.
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
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_result_type, apply_seed, &
      & seed_weight_tol, next_branch_group, max_branch_group_size
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent_work_type, &
      & drop_field_jet_point, drop_field_jet_tangent
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !> Forward tangent of the four grid scalars the weight fold reads
   !>
   !> The DROP-internal form of pass 1. Every output is `(ngrid, ndir)` and is
   !> written in full: column `idir` holds the directional derivative along
   !> `dirs(:, :, idir)`. A thin wrapper over [[drop_surface_tangent_core]]
   !> that keeps the branch-weight tangent, which the generic
   !> [[cavity_surface_tangent_type]] has no channel for.
   !>
   !> @param[in]  self      DROP cavity instance (must hold a projected grid)
   !> @param[in]  dirs      Nuclear directions `(3, nsph, ndir)`
   !> @param[out] d_a       Tangent of the area element `(ngrid, ndir)`
   !> @param[out] d_wleb    Tangent of the Lebedev weight `(ngrid, ndir)`
   !> @param[out] d_xi0     Tangent of the Gaussian width `(ngrid, ndir)`
   !> @param[out] d_wbranch Tangent of the branch weight `(ngrid, ndir)`
   !> @param[out] error     Error object, allocated on failure
   !> @param[in]  contracted Take the jet tangents through the level set's
   !>                        contracted accessor per direction rather than off
   !>                        tensors materialised once per point; the caller's
   !>                        choice, see the core. Default `.false.`
   module subroutine get_surface_tangent_drop(self, dirs, d_a, d_wleb, d_xi0, &
                                              d_wbranch, error, contracted)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Directional tangents of the four grid scalars
      real(wp), intent(out) :: d_a(:, :), d_wleb(:, :), d_xi0(:, :), d_wbranch(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Whether the jet tangents come from the contracted accessor
      logical, intent(in), optional :: contracted

      !> Every channel of the tangent; only three are copied out
      type(cavity_surface_tangent_type) :: tangent
      !> Direction count
      integer :: ndir

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

      call tangent%init(self%ngrid, ndir, .false.)
      call drop_surface_tangent_core(self, dirs, .false., tangent, error, &
                                     d_wbranch=d_wbranch, contracted=contracted)
      if (allocated(error)) return
      d_a = tangent%d_a
      d_wleb = tangent%d_w
      d_xi0 = tangent%d_xi

   end subroutine get_surface_tangent_drop

   !> Forward tangent of every surface observable along a set of nuclear directions
   !>
   !> @param[in]    self    DROP cavity instance (must hold a projected grid)
   !> @param[in]    dirs    Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] tangent Surface tangent, initialised for `(ngrid, ndir)`
   !> @param[out]   error   Error object, allocated on failure
   module subroutine get_surface_tangent_full_drop(self, dirs, tangent, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Surface tangent
      type(cavity_surface_tangent_type), intent(inout) :: tangent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call drop_surface_tangent_core(self, dirs, tangent%have_curvature, tangent, error)

   end subroutine get_surface_tangent_full_drop

   !> Forward tangent of the DROP surface map along a set of nuclear directions
   !>
   !> Every channel of `tangent` is written in full: column `idir` holds the
   !> directional derivative along `dirs(:, :, idir)`. The curvature channels
   !> are filled only when `want_curvature` is set, which puts the curvature
   !> invariants into the seed state of every point; without it they are left
   !> zero and `tangent%have_curvature` says so. `d_wbranch`, when present,
   !> receives the softmax branch-weight tangent, which the generic tangent
   !> has no channel for and which pass 2 of the Hessian reads on its own.
   !>
   !> @param[in]    self           DROP cavity instance (must hold a projected grid)
   !> @param[in]    dirs           Nuclear directions `(3, nsph, ndir)`
   !> @param[in]    want_curvature Fill the curvature channels
   !> @param[inout] tangent        Surface tangent, initialised for `(ngrid, ndir)`
   !> @param[out]   error          Error object, allocated on failure
   !> @param[out]   d_wbranch      Tangent of the branch weight `(ngrid, ndir)`
   !> @param[in]    contracted     Take the jet tangents through the level set's
   !>                              contracted accessor per direction rather than
   !>                              off tensors materialised once per point; the
   !>                              caller's choice, see below. Default `.false.`
   module subroutine drop_surface_tangent_core(self, dirs, want_curvature, tangent, &
                                               error, d_wbranch, contracted)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Whether the curvature channels are filled
      logical, intent(in) :: want_curvature
      !> Surface tangent
      type(cavity_surface_tangent_type), intent(inout) :: tangent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Tangent of the softmax branch weight
      real(wp), intent(out), optional :: d_wbranch(:, :)
      !> Whether the jet tangents come from the contracted accessor
      logical, intent(in), optional :: contracted

      !> Branch-weight tangent, kept whether or not the caller asked for it:
      !> the assembly below reads it for every point
      real(wp), allocatable :: dwb(:, :)

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Thread bookkeeping
      integer :: thread_slot
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure on its way to the latch
      type(error_type), allocatable :: worker_error

      !> Per-thread state of the grid point being opened: the jets, the seed
      !> state, the bordered factorization and the buffers they are read into
      type(drop_point_scratch_type) :: pt
      !> Linear response of one seed
      type(drop_seed_result_type) :: res
      !> Whether the shared prologue cleared the point for this traversal
      logical :: point_ok

      !> Grid, direction, sphere, active-slot and Cartesian indices
      integer :: igrid, idir, ndir, jj, knb, i
      !> Whether the jet tensors are materialised per point
      logical :: materialise

      !> Mixed nuclear tensors of the level set at the point
      type(drop_field_tangent_work_type) :: ft_work
      !> One direction gathered onto the point's active slots
      real(wp), allocatable :: v_act(:, :)

      !> Directional nuclear tangents of the level-set jet at the *fixed* point
      real(wp) :: dlsf0
      real(wp), allocatable :: dlsf1_r(:, :), dlsf2_rr(:, :, :)
      !> Bordered KKT right-hand sides, one column per direction; this pass
      !> builds its own batch rather than the shared 16-seed one
      real(wp), allocatable :: dir_rhs(:, :)
      !> Induced motion of the projected point and of the multiplier
      real(wp) :: dr(3), dlambda

      !> Sparse switching rows of the owner sphere
      real(wp) :: swi_owner_row(3), swi_f0, swi_dxi, df_dir

      !> Tangent of the branch objective, `d(Phi)` per point and direction
      real(wp), allocatable :: dphi(:, :)
      !> Softmax scratch of the branch stage
      real(wp), allocatable :: branch_phi(:), branch_dphi(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)
      !> Branch group bookkeeping
      integer :: igroup_cursor, igroup_start, igroup_end, group_size
      integer :: m_branch, im_grid, nbranch_max
      !> Whether a further anchor group exists
      logical :: have_group

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
      if (.not. tangent%is_initialized()) then
         call fatal_error(error, "get_surface_tangent_drop: tangent is not initialized")
         return
      end if
      if (size(tangent%d_a, 1) /= self%ngrid .or. size(tangent%d_a, 2) /= ndir) then
         call fatal_error(error, "get_surface_tangent_drop: tangent must be initialized"// &
                          " for (ngrid, ndir)")
         return
      end if
      if (present(d_wbranch)) then
         if (size(d_wbranch, 1) /= self%ngrid .or. size(d_wbranch, 2) /= ndir) then
            call fatal_error(error, "get_surface_tangent_drop: d_wbranch must be"// &
                             " (ngrid, ndir)")
            return
         end if
         d_wbranch = 0.0_wp
      end if

      tangent%d_xi = 0.0_wp
      tangent%d_f = 0.0_wp
      tangent%d_a = 0.0_wp
      tangent%d_w = 0.0_wp
      tangent%d_xyz = 0.0_wp
      tangent%d_n = 0.0_wp
      tangent%d_k1 = 0.0_wp
      tangent%d_k2 = 0.0_wp
      tangent%have_curvature = want_curvature
      if (self%ngrid <= 0) return
      allocate (dwb(self%ngrid, ndir), source=0.0_wp)

      h_stan = self%ctx%timer%resolve("Surface tangent", self%ctx%timer%current(), &
                                      cat_gradient)
      call self%ctx%timer%start(h_stan)

      !* -------------------------------- Thread setup -------------------------------- *!
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      allocate (dphi(self%ngrid, ndir), source=0.0_wp)
      ! The jet tensors are materialised once per point and contracted with every
      ! direction, or the level set's own contracted accessor runs per direction;
      ! for a few directions the latter costs less than one fill, for many the
      ! former. The choice is the caller's and deliberately not a function of
      ! `ndir`: the Hessian traversal hands this routine one *block* of its
      ! directions at a time, and a choice keyed on the block size would let a
      ! direction's tangent depend on how the set happened to be blocked, by an
      ! ulp, which the exact chunking guarantee of that traversal forbids. Its
      ! per-direction mode asks for the contracted path; everything else, this
      ! routine's direct callers included, materialises.
      materialise = .true.
      if (present(contracted)) materialise = .not. contracted

      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& pt, point_ok, idir, jj, knb, i, res, ft_work, v_act, &
      !$omp& dlsf0, dlsf1_r, dlsf2_rr, dir_rhs, dr, dlambda, &
      !$omp& swi_owner_row, swi_f0, swi_dxi, df_dir, &
      !$omp& worker_error)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      call pt%init(self%nsph, self%iswig)
      allocate (dlsf1_r(3, ndir), source=0.0_wp)
      allocate (dlsf2_rr(3, 3, ndir), source=0.0_wp)
      allocate (dir_rhs(4, ndir), source=0.0_wp)
      allocate (v_act(3, self%nsph), source=0.0_wp)

      !$omp do schedule(static, 8)
      do igrid = 1, self%ngrid
         !* -------------------------- Shared point prologue -------------------------- *!
         ! Point, jets, seed state and the bordered factorization. Every failure
         ! path -- a latch already set, a refusing level set, a degenerate
         ! state, a singular bordered system -- has recorded itself. The
         ! standard 16-seed batch is *not* requested: this pass has one
         ! right-hand side per nuclear direction, and it cannot be built before
         ! the level set's directional tangents below are known.
         call drop_point_prologue(self, slots, thread_slot, igrid, want_curvature, &
                                  "get_surface_tangent_drop", abort, pt, point_ok)
         if (.not. point_ok) cycle

         !* ----------------- Directional nuclear tangents of the jet ----------------- *!
         ! `sum_B v_B . d(jet)/dR_B` at the fixed point, for every direction, as
         ! contractions of the point's mixed nuclear tensors. The tensors are
         ! formed once per point and unconditionally -- the buffer outlives the
         ! point -- and each direction is gathered onto the active slots first,
         ! which is the one place this routine touches the slot index space.
         pt%n_active = slots%lsf(thread_slot)%lsf%active_count()
         do i = 1, pt%n_active
            pt%active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
         end do
         if (materialise) call drop_field_jet_point(slots%lsf(thread_slot)%lsf, ft_work)
         do idir = 1, ndir
            if (materialise) then
               do i = 1, pt%n_active
                  v_act(:, i) = dirs(:, pt%active_idx(i), idir)
               end do
               call drop_field_jet_tangent(ft_work, pt%n_active, v_act, dlsf0, &
                                           dlsf1_r(:, idir), dlsf2_rr(:, :, idir))
            else
               call slots%lsf(thread_slot)%lsf%tangent_jet(dirs(:, :, idir), dlsf0, &
                                                            dlsf1_r(:, idir), dlsf2_rr(:, :, idir))
            end if

            ! Bordered right-hand side of the direction, with
            ! `d^2 phi/(dr dR_owner) = -alpha*I` (`objective_phi.f90`, `f2_r_rA`)
            ! and the anchor riding its owner rigidly.
            dir_rhs(1:3, idir) = self%param%phi_alpha*dirs(:, pt%owner_idx, idir) &
                                 + pt%lambda_val*dlsf1_r(:, idir)
            dir_rhs(4, idir) = -dlsf0
         end do

         !* ------------------------ Bordered KKT sensitivities ----------------------- *!
         ! One batched solve on the factorization the prologue already built:
         ! every direction shares the 4x4 matrix of this grid point.
         call pt%kkt_fac%solve(dir_rhs, "get_surface_tangent_drop", worker_error, igrid)
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            cycle
         end if

         !* ---------------------- Base Lebedev-weight response ----------------------- *!
         do idir = 1, ndir
            dr = dir_rhs(1:3, idir)
            dlambda = dir_rhs(4, idir)

            ! The projected point moves with the bordered solve
            tangent%d_xyz(:, igrid, idir) = dr

            call apply_seed(pt%state, dlsf1_r(:, idir), dlsf2_rr(:, :, idir), dr, dlambda, res)

            ! Branch-frozen half of `d(wleb)`; stage 3 completes it. `res%dxi`
            ! is the width tangent of exactly this incomplete weight and is not
            ! read at all.
            tangent%d_w(igrid, idir) = res%dwleb

            ! The outward normal `grad S/|grad S|` at the moving point, and the
            ! principal curvatures when the seed state carries them
            tangent%d_n(:, igrid, idir) = res%dn_surf
            if (want_curvature) then
               tangent%d_k1(igrid, idir) = res%dk1
               tangent%d_k2(igrid, idir) = res%dk2
            end if

            ! Tangent of the branch objective `Phi = 0.5 alpha |r* - anchor|^2`
            ! along the direction: the projected point moves by `dr`, the anchor
            ! rigidly with its owner.
            dphi(igrid, idir) = dot_product(pt%phi1_r, dr - dirs(:, pt%owner_idx, idir))
         end do

         !* ------------------------- iSwiG switching channel ------------------------- *!
         ! `f` is evaluated at the anchor with the anchor width, and
         ! `anchor_xi0` depends on the owner radius and the raw Lebedev weight
         ! alone, so it carries no nuclear tangent and `swi_dxi` is unused.
         call self%iswig%swi_collect(pt%anchor, pt%owner_idx, self%anchor_xi0(igrid), &
                                     swi_f0, pt%iswig_work)
         call self%iswig%swi1_rA_sparse(pt%iswig_work, pt%swi_rows, swi_owner_row, swi_dxi)
         do idir = 1, ndir
            df_dir = dot_product(swi_owner_row, dirs(:, pt%owner_idx, idir))
            do jj = 1, pt%iswig_work%n_nb
               knb = pt%iswig_work%idx(jj)
               df_dir = df_dir + dot_product(pt%swi_rows(:, jj), dirs(:, knb, idir))
            end do
            tangent%d_f(igrid, idir) = df_dir
         end do

      end do
      !$omp end do

      deallocate (dlsf1_r, dlsf2_rr, dir_rhs, v_act)
      call pt%destroy()
      !$omp end parallel

      if (abort%requested) then
         call abort%raise("get_surface_tangent_drop", error)
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
            dwleb_i = tangent%d_w(igrid, idir)
            if (wbranch_i > tiny(1.0_wp)) then
               dwleb_i = dwleb_i + (wleb_i/wbranch_i)*dwb(igrid, idir)
            end if
            tangent%d_w(igrid, idir) = dwleb_i

            ! xi0 = swx/(R sqrt(wleb)); same guard as `apply_seed` and `iswig_xi0`
            if (wleb_i > seed_weight_tol) then
               tangent%d_xi(igrid, idir) = -0.5_wp*self%xi0(igrid)*dwleb_i/wleb_i
            else
               tangent%d_xi(igrid, idir) = 0.0_wp
            end if

            ! a = R^2 f wleb, with `f` the iSwiG switching factor
            tangent%d_a(igrid, idir) = r_own*r_own &
                                       *(wleb_i*tangent%d_f(igrid, idir) + self%f(igrid)*dwleb_i)
         end do
      end do
      if (present(d_wbranch)) d_wbranch = dwb

      call self%ctx%timer%stop(h_stan)

   contains

      !> Differentiate the branch softmax over every contiguous anchor group
      !>
      !> Walks the groups with the shared [[next_branch_group]], exactly as
      !> [[compute_branch_phi_adj]] and `forward.f90`'s branch post-pass do:
      !> runs of equal `anchor_id` starting at a point with `branch_count > 1`.
      !> The softmax primitive takes its derivatives in `(nparam, nbranch)`
      !> layout, so passing `nparam = ndir` yields every direction of the group
      !> from one call.
      !>
      !> The scratch is sized by [[max_branch_group_size]] and not by
      !> `maxval(branch_count)`: it is indexed by the group's run length, which
      !> is what that function returns and what the walk below produces.
      subroutine branch_stage()

         if (.not. allocated(self%branch_count) .or. .not. allocated(self%anchor_id)) return
         if (.not. any(self%branch_count(1:self%ngrid) > 1)) return

         nbranch_max = max_branch_group_size(self%branch_count(1:self%ngrid), &
                                             self%anchor_id(1:self%ngrid))
         allocate (branch_phi(nbranch_max), source=0.0_wp)
         allocate (branch_dphi(ndir, nbranch_max), source=0.0_wp)
         allocate (branch_weights(nbranch_max), source=0.0_wp)
         allocate (branch_dweights(ndir, nbranch_max), source=0.0_wp)

         igroup_cursor = 1
         do
            call next_branch_group(self%branch_count(1:self%ngrid), &
                                   self%anchor_id(1:self%ngrid), igroup_cursor, &
                                   igroup_start, igroup_end, have_group)
            if (.not. have_group) exit
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
               dwb(im_grid, :) = branch_dweights(:, m_branch)
            end do
         end do

         deallocate (branch_phi, branch_dphi, branch_weights, branch_dweights)

      end subroutine branch_stage

   end subroutine drop_surface_tangent_core

end submodule moist_cavity_drop_derivatives_tangent_forward
