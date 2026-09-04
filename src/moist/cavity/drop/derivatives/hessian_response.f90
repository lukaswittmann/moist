!> Adjoint-response half of the DROP surface Hessian
!>
!> The nuclear gradient of `nuclear.f90` is `J^T omega`. Its directional
!> derivative splits into
!>
!>     d/dv [ J^T omega ]  =  (dJ^T/dv) omega  +  J^T (d omega/dv)
!>
!> [[get_surface_hessian_fixed_drop]] builds the first term, with the adjoints
!> held fixed and the surface map differentiated. This submodule builds the
!> second: the **primal map is held fixed and the adjoints move**.
!>
!> ## Why this is the gradient traversal again, and not new mathematics
!>
!> The raw host adjoints of a `cavity_surface_adjoint_type` never move -- they
!> belong to the energy expression, not to the geometry. What moves is
!> [[prepare_surface_weights]]'s *folding* of them: `w_xi` absorbs the area and
!> integration-weight channels through `a`, `wleb` and `xi0`, `w_f` absorbs the
!> area channel through `R^2 wleb`, and `branch_phi_adj` is derived from the
!> branch softmax. All three are functions of the geometry.
!>
!> Write the gradient as `G(R) = Phi(R) . eff(R)`, with `Phi` the per-point
!> primal map [[get_surface_gradient_drop]] contracts and `eff` the folded
!> weights. `G` is **linear** in `eff` -- every channel enters exactly once, in
!> the seed contraction, in the normal fold, in the field row and in the
!> switching row -- so
!>
!>     d/dv G  =  (dPhi/dv) . eff   +   Phi . (d eff/dv)
!>                ^ hessian_fixed.f90    ^ this submodule
!>
!> and the second term is nothing but `get_surface_gradient_drop` run with
!> `d(eff)` substituted for `eff`. That is what this routine does, seed for
!> seed, and it is why it reuses [[seed_jet_basis]] and [[seed_anchor]]
!> verbatim rather than restating their contractions.
!>
!> A consequence worth stating, because the second-order chain next door has to
!> worry about it and this one does not: the traversal here is **first order in
!> the seed chain**. It calls [[apply_seed]] (through the two seed routines) and
!> never [[apply_seed_tangent]], so the `PRECONDITION` on that routine -- the
!> nine single-entry, asymmetric Hessian jet seeds whose individual tangents are
!> wrong by order 100 % -- does not apply. The seeds still have to be read as a
!> contracted set rather than one at a time, and they are: their contributions
!> land in `w_lsf2_pt` and leave only through `vjp_f1_rA`, against the level
!> set's mixed third derivative, which is symmetric in its two spatial indices.
!> That is the same rescue the shipped first-order adjoint has always relied on.
!>
!> ## The three passes
!>
!>  1. **Forward tangent** -- [[get_surface_tangent_drop]] pushes every nuclear
!>     direction through the per-point map and returns `d_a`, `d_wleb`, `d_xi0`
!>     and `d_wbranch`, one column per direction.
!>  2. **Weight tangent** -- [[prepare_surface_weights_tangent]] turns those
!>     four into `d(eff)`. Only three channels move: `w_xi`, `w_f` and
!>     `branch_phi_adj`. `w_xyz`, `w_n`, `w_k1` and `w_k2` are `source=`-copies
!>     of the raw adjoints, so their tangent is identically zero and pass 2
!>     deliberately does not emit them.
!>  3. **Contraction** -- the gradient traversal, with `deff` in place of `eff`.
!>
!> Because pass 2 emits three channels and not seven, the `deff` objects this
!> routine builds carry `w_xi`, `w_f` and `branch_phi_adj` and nothing else.
!> Two consequences follow directly, and both are taken:
!>
!>   * the **normal channel is skipped outright**. With `d(w_xyz)` and `d(w_n)`
!>     identically zero, [[seed_normal_channel]]'s output is the zero vector at
!>     every point, so the effective position adjoint the seeds see is a hard
!>     zero rather than something to compute. `deff%w_xyz` is therefore not
!>     allocated at all, and the fold is not called;
!>   * the **curvature channel is skipped** for the same reason, through
!>     `deff%have_wk = .false.` -- which also means `fill_seed_state` is asked
!>     for `want_curvature = .false.`, exactly as the forward tangent asks.
!>
!> ## The group reduction, and where it is allowed to be parallel
!>
!> [[branch_phi_adj_tangent]] reduces one scalar over each contiguous anchor
!> group and is the only cross-point coupling in the scheme. A group split
!> across two threads corrupts its reduction silently. This routine therefore
!> runs the whole of pass 2 **serially, one call per direction**, outside any
!> parallel region: the primitive is itself serial over the entire grid, so
!> every group is seen whole by one call, and directions are independent of each
!> other by construction.
!>
!> The grid loop that follows is parallel over grid points, and that is safe
!> because it is not the reduction -- `dbranch_phi_adj` is a finished per-point
!> array by the time the loop starts, read exactly as
!> [[get_surface_gradient_drop]] reads the primal `branch_phi_adj`.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_hessian_response
!$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_state_type, &
      & drop_surface_weights_type, build_seed_state, seed_state_ok, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type, &
      & seed_jet_basis, seed_anchor, degenerate_point_error
   use moist_cavity_drop_derivatives_weights_tangent, only: prepare_surface_weights_tangent
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !> Contract the moving surface adjoints against the primal map
   !>
   !> Accumulates `J^T (d omega/dv)` for the energy whose raw surface adjoints
   !> `acc` holds, one gradient column per nuclear direction. The result is
   !> *added* to `hvp`, so the caller may already hold the fixed half, and the
   !> accumulator is left untouched when anything fails.
   !>
   !> Mirrors [[get_surface_gradient_drop]] throughout: same thread setup, same
   !> per-point jet and bordered KKT solve, same error latching and the same
   !> deterministic reduction. What differs is the weights -- `d(eff)` rather
   !> than `eff` -- and the direction loop those weights force inside the grid
   !> loop.
   !>
   !> @param[in]    self  DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc   Raw surface-observable adjoints, held fixed
   !> @param[in]    dirs  Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hvp   Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]   error Error object, allocated on failure
   module subroutine get_surface_hessian_response_drop(self, acc, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Hessian-vector accumulator
      real(wp), intent(inout) :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Per-thread accumulators, summed deterministically after the region
      real(wp), allocatable :: hvp_threads(:, :, :, :)
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

      !> Grid, direction, atom and active-slot indices
      integer :: igrid, idir, ndir, iatom, i, n_active
      integer, allocatable :: active_idx(:)

      !> Projected point, anchor and owner sphere
      real(wp) :: point(3), anchor(3)
      integer :: owner_idx
      !> Level-set jet at the projected point
      real(wp) :: lsf0, lsf1_r(3), lsf2_rr(3, 3)
      real(wp), allocatable :: lsf3_rrr(:, :, :)
      !> Jet-contracted nuclear partials of the level set, one column per active atom
      real(wp), allocatable :: vjp_pt(:, :)
      !> Objective jet at the projected point
      real(wp) :: phi0, phi1_r(3), phi2_rr(3, 3)
      !> Lagrange multiplier of the projection
      real(wp) :: lambda_val

      !> Bordered KKT system with four field seeds and three anchor seeds
      real(wp) :: kkt_rhs(4, 7)
      !> Factorization reused by every solve at this grid point
      type(drop_kkt_factor_type) :: kkt_fac

      !> Point-local level-set adjoint weights built from the 13 field seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Effective position adjoint seen by every seed; identically zero here
      real(wp) :: w_xyz_local(3)
      !> iSwig neighbour cache and the sparse switching rows it feeds
      type(iswig_workspace_type) :: iswig_work
      real(wp), allocatable :: swi_rows(:, :)
      real(wp) :: swi_owner_row(3), swi_f0, swi_dxi
      logical :: swi_live
      integer :: jj, knb

      !> Primal folded surface adjoints, and their tangent along each direction
      type(drop_surface_weights_type) :: eff
      type(drop_surface_weights_type), allocatable :: deff(:)
      !> Timer handle
      integer :: h_shres

      !* ------------------------------- Shape guards --------------------------------- *!
      call check_surface_adjoint(self, acc, "get_surface_hessian_response_drop", error)
      if (allocated(error)) return
      if (size(dirs, 1) /= ndim .or. size(dirs, 2) /= self%nsph) then
         call fatal_error(error, "get_surface_hessian_response_drop: dirs must be"// &
                          " (3, nsph, ndir)")
         return
      end if
      ndir = size(dirs, 3)
      if (ndir <= 0) then
         call fatal_error(error, "get_surface_hessian_response_drop: no direction supplied")
         return
      end if
      if (size(hvp, 1) /= ndim .or. size(hvp, 2) /= self%nsph .or. size(hvp, 3) /= ndir) then
         call fatal_error(error, "get_surface_hessian_response_drop: hvp must be"// &
                          " (3, nsph, ndir)")
         return
      end if
      if (self%ngrid <= 0) return

      h_shres = self%ctx%timer%resolve("Surface Hessian (adjoint response)", &
                                       self%ctx%timer%current(), cat_gradient)
      call self%ctx%timer%start(h_shres)

      !* ----------------------- Passes 1 and 2: the moving weights -------------------- *!
      call prepare_surface_weights(self, acc, .true., eff)
      call weight_tangents(self, acc, eff, dirs, deff, error)
      if (allocated(error)) then
         call self%ctx%timer%stop(h_shres)
         return
      end if

      !* -------------------------------- Thread setup -------------------------------- *!
      ! Order 3, as on the gradient path: the primal map is the one being
      ! contracted, so nothing above the third jet derivative is read.
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      allocate (hvp_threads(3, self%nsph, ndir, slots%nthreads), source=0.0_wp)

      call abort%reset()

      !$omp parallel num_threads(slots%nthreads) default(shared) private(thread_slot, igrid, &
      !$omp& idir, iatom, i, n_active, active_idx, state, status, &
      !$omp& point, anchor, owner_idx, lsf0, lsf1_r, lsf2_rr, lsf3_rrr, &
      !$omp& vjp_pt, phi0, phi1_r, phi2_rr, lambda_val, &
      !$omp& kkt_rhs, kkt_fac, &
      !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_local, &
      !$omp& iswig_work, swi_rows, swi_owner_row, swi_f0, swi_dxi, swi_live, &
      !$omp& jj, knb, worker_error)
      thread_slot = 1
!$    thread_slot = omp_get_thread_num() + 1

      allocate (lsf3_rrr(3, 3, 3), source=0.0_wp)
      allocate (vjp_pt(3, self%nsph), source=0.0_wp)
      allocate (active_idx(self%nsph))
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
         ! No curvature: `d(w_k1)` and `d(w_k2)` vanish identically whatever the
         ! host put in `acc`, so this half never reads a curvature response.
         call fill_seed_state(self, igrid, .false., state)

         call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                               self%param%wleb_prune_level > 0, status)
         if (status /= seed_state_ok) then
            call abort%latch_status(status, igrid)
            cycle
         end if

         ! Outward-normal channel: `d(w_xyz)` and `d(w_n)` are identically zero
         ! (`weights_tangent.f90`, module header -- a consumer of pass 2 "gains
         ! the right to skip the normal and curvature channels of the tangent
         ! contraction outright"), so the effective position adjoint every seed
         ! sees is the zero vector and [[seed_normal_channel]] is not called.
         w_xyz_local = 0.0_wp

         !* ------------------------ Bordered KKT sensitivities ----------------------- *!
         ! Columns 1-4 are the level-set value and gradient seeds; the nine
         ! Hessian seeds have a zero right-hand side. Columns 5-7 are the
         ! anchor seeds: moving the owner rigidly leaves the field untouched
         ! and drives the system through -d^2 phi/dr dR_owner = +alpha*I.
         ! Direction independent, so solved once for the whole direction loop.
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

         !* --------------------- Direction-independent point data -------------------- *!
         n_active = slots%lsf(thread_slot)%lsf%active_count()
         do i = 1, n_active
            active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
         end do

         ! The iSwiG rows depend on the geometry alone; only the scalar weight
         ! in front of them is per direction. Collected once when any direction
         ! carries a live switching tangent.
         swi_live = .false.
         do idir = 1, ndir
            if (abs(deff(idir)%w_f(igrid)) > seed_weight_tol) swi_live = .true.
         end do
         if (swi_live) then
            call self%iswig%swi_collect(anchor, owner_idx, self%anchor_xi0(igrid), &
                                        swi_f0, iswig_work)
            call self%iswig%swi1_rA_sparse(iswig_work, swi_rows, swi_owner_row, swi_dxi)
         end if

         !* ------------------------------ Direction loop ----------------------------- *!
         do idir = 1, ndir

            !* --------------- Field seeds -> level-set adjoint tangents -------------- *!
            w_lsf0_pt = 0.0_wp
            w_lsf1_pt = 0.0_wp
            w_lsf2_pt = 0.0_wp
            call seed_jet_basis(state, deff(idir), igrid, phi1_r, kkt_rhs, w_xyz_local, &
                                w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

            !* ------------- Field channel: contract with nuclear partials ------------ *!
            ! As on the gradient path, the level set contracts the jet indices
            ! itself, so the mixed third derivative is never materialized.
            call slots%lsf(thread_slot)%lsf%vjp_f1_rA(w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, vjp_pt)
            do i = 1, n_active
               iatom = active_idx(i)
               hvp_threads(:, iatom, idir, thread_slot) = &
                  hvp_threads(:, iatom, idir, thread_slot) + vjp_pt(:, i)
            end do

            !* ------------------------ Anchor channel (owner) ------------------------ *!
            call seed_anchor(state, deff(idir), igrid, phi1_r, kkt_rhs, w_xyz_local, &
                             hvp_threads(:, owner_idx, idir, thread_slot))

            !* ----------------------- iSwig switching channel ------------------------ *!
            if (abs(deff(idir)%w_f(igrid)) > seed_weight_tol) then
               do jj = 1, iswig_work%n_nb
                  knb = iswig_work%idx(jj)
                  hvp_threads(:, knb, idir, thread_slot) = &
                     hvp_threads(:, knb, idir, thread_slot) &
                     + deff(idir)%w_f(igrid)*swi_rows(:, jj)
               end do
               hvp_threads(:, owner_idx, idir, thread_slot) = &
                  hvp_threads(:, owner_idx, idir, thread_slot) &
                  + deff(idir)%w_f(igrid)*swi_owner_row
            end if

         end do

      end do
      !$omp end do

      deallocate (lsf3_rrr, vjp_pt, active_idx, swi_rows)
      call iswig_work%destroy()
      !$omp end parallel

      if (abort%requested) then
         ! An LSF failure or a KKT failure arrives as a ready-made error; a
         ! kernel degeneracy arrives as a status code that needs this routine's
         ! name to become a diagnostic.
         if (allocated(abort%error)) then
            call move_alloc(abort%error, error)
         else
            call degenerate_point_error("get_surface_hessian_response_drop", abort%status, &
                                        abort%igrid, error)
         end if
         call self%ctx%timer%stop(h_shres)
         return
      end if

      ! Deterministic reduction: fixed thread order, independent of scheduling
      do ithread = 1, slots%nthreads
         hvp = hvp + hvp_threads(:, :, :, ithread)
      end do

      call self%ctx%timer%stop(h_shres)

   end subroutine get_surface_hessian_response_drop

   !* ================================================================================= *!
   !*                          Passes 1 and 2: d(eff) per direction                     *!
   !* ================================================================================= *!

   !> Build the directional tangent of the folded surface weights
   !>
   !> Runs the forward tangent once for the whole batch and the weight tangent
   !> once per direction. The second loop is **serial on purpose**: it is the
   !> stage that carries [[branch_phi_adj_tangent]]'s group reduction, the only
   !> cross-point coupling in the scheme, and the primitive is serial over the
   !> whole grid, so every contiguous anchor group is seen whole by exactly one
   !> call. Directions are the admissible parallel axis if this ever needs one;
   !> grid points are not, and never will be.
   !>
   !> `deff(idir)` carries the three channels pass 2 emits and nothing else --
   !> see the module header for why the other five are absent rather than zero.
   !>
   !> @param[in]  self  DROP cavity instance
   !> @param[in]  acc   Raw surface adjoints, held fixed
   !> @param[in]  eff   Folded weights, as [[prepare_surface_weights]] returned them
   !> @param[in]  dirs  Nuclear directions `(3, nsph, ndir)`
   !> @param[out] deff  Tangent of the folded weights, one element per direction
   !> @param[out] error Error object, allocated on failure
   subroutine weight_tangents(self, acc, eff, dirs, deff, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded weights
      type(drop_surface_weights_type), intent(in) :: eff
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Tangent of the folded weights, one element per direction
      type(drop_surface_weights_type), allocatable, intent(out) :: deff(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Directional tangents of the four grid scalars pass 2 consumes
      real(wp), allocatable :: d_a(:, :), d_wleb(:, :), d_xi0(:, :), d_wbranch(:, :)
      !> Branch bookkeeping, defaulted when the cavity carries none
      integer, allocatable :: branch_count(:), anchor_id(:)
      !> Grid extent and direction index
      integer :: ngrid, ndir, idir

      ngrid = self%ngrid
      ndir = size(dirs, 3)

      !* --------------------------- Pass 1: forward tangent -------------------------- *!
      allocate (d_a(ngrid, ndir), d_wleb(ngrid, ndir), d_xi0(ngrid, ndir), &
                d_wbranch(ngrid, ndir))
      call self%get_surface_tangent(dirs, d_a, d_wleb, d_xi0, d_wbranch, error)
      if (allocated(error)) return

      ! `compute_branch_phi_adj` is skipped by the primal when the cavity holds
      ! no branch bookkeeping; a single-branch stand-in reproduces that early
      ! exit in the tangent instead of duplicating the guard.
      allocate (branch_count(ngrid), source=1)
      allocate (anchor_id(ngrid), source=0)
      if (allocated(self%branch_count)) branch_count = self%branch_count(1:ngrid)
      if (allocated(self%anchor_id)) anchor_id = self%anchor_id(1:ngrid)

      !* --------------------------- Pass 2: weight tangent --------------------------- *!
      allocate (deff(ndir))
      do idir = 1, ndir
         allocate (deff(idir)%w_xi(ngrid), deff(idir)%w_f(ngrid), &
                   deff(idir)%branch_phi_adj(ngrid))
         ! `have_wn` and `have_wk` keep their `.false.` default: the normal and
         ! curvature adjoints are copies of the fixed raw channels, so their
         ! tangent is identically zero.
         call prepare_surface_weights_tangent(acc, eff, .true., &
                                              self%a(1:ngrid), self%wleb(1:ngrid), &
                                              self%xi0(1:ngrid), self%wbranch(1:ngrid), &
                                              self%radii, self%owner(1:ngrid), &
                                              branch_count, anchor_id, &
                                              self%branch_weight%s, &
                                              d_a(:, idir), d_wleb(:, idir), &
                                              d_xi0(:, idir), d_wbranch(:, idir), &
                                              deff(idir)%w_xi, deff(idir)%w_f, &
                                              deff(idir)%branch_phi_adj, error)
         if (allocated(error)) return
      end do

   end subroutine weight_tangents

end submodule moist_cavity_drop_derivatives_hessian_response
