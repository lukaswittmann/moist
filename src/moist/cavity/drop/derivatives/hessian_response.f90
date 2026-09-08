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
!> seed, and it is why it reuses the shipped jet and anchor bases verbatim
!> rather than restating their contractions.
!>
!> It reuses them in halves, though, not whole. [[seed_jet_basis]] and
!> [[seed_anchor]] each seed the kernel and contract the response in one sweep,
!> which is right for a caller with one weight set per point; this one has
!> `ndir` of them and the seeding depends on none. So the traversal calls
!> [[seed_jet_basis_apply]] and [[seed_anchor_apply]] once per grid point and
!> [[seed_jet_basis_contract]] and [[seed_anchor_contract]] once per direction,
!> which is the same shape [[get_surface_hessian_fixed_drop]] already has.
!>
!> A consequence worth stating, because the second-order chain next door has to
!> worry about it and this one does not: the traversal here is **first order in
!> the seed chain**. It calls [[apply_seed]] (through the two apply halves) and
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
!>  1. **Forward tangent** -- [[get_surface_tangent_drop]] pushes each nuclear
!>     direction through the per-point map and returns `d_a`, `d_wleb`, `d_xi0`
!>     and `d_wbranch`, one column per direction.
!>  2. **Weight tangent** -- [[prepare_surface_weights_tangent]] turns those
!>     four into `d(eff)`. Only three channels move: `w_xi`, `w_f` and
!>     `branch_phi_adj`. `w_xyz`, `w_n`, `w_k1` and `w_k2` are `source=`-copies
!>     of the raw adjoints, so their tangent is identically zero and pass 2
!>     deliberately does not emit them.
!>  3. **Contraction** -- the gradient traversal, with `deff` in place of `eff`.
!>
!> All three run on a **block** of directions rather than on all of them at
!> once, because passes 1 and 2 are the only part of this scheme whose memory
!> grows with the direction count and nothing above bounds that count; see
!> `hvp_chunk_dirs` for the bound, the price and why the common case pays
!> nothing. Blocks are disjoint in the direction index and independent of each
!> other, so the blocking is invisible in the result down to the last bit.
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
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_result_type, &
      & drop_surface_weights_type, seed_weight_tol
   use moist_cavity_drop_derivatives_seeds, only: &
      & seed_jet_basis_apply, seed_jet_basis_contract, &
      & seed_anchor_apply, seed_anchor_contract, &
      & drop_n_jet_seeds, drop_n_anchor_seeds
   use moist_cavity_drop_derivatives_weights_tangent, only: prepare_surface_weights_tangent
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

   !> Nuclear directions carried by one block of the traversal
   !>
   !> Passes 1 and 2 are per direction and materialize seven `(ngrid, ndir)`
   !> grid arrays -- the four forward tangents `d_a`, `d_wleb`, `d_xi0`,
   !> `d_wbranch` and the three moving channels of `deff` -- and the thread
   !> accumulator adds a `(3, nsph, ndir, nthreads)` one. Nothing above this
   !> submodule bounds `ndir`: [[get_hessian_drop]] asks for `3 nsph`
   !> directions, so all eight grow *quadratically* in the system size and a
   !> medium molecule runs out of memory rather than running slowly.
   !>
   !> The traversal is therefore blocked: directions are processed
   !> `hvp_chunk_dirs` at a time and every array above is sized by the block
   !> rather than by `ndir`, which bounds the working set at
   !>
   !>     bytes  =  8 C (7 ngrid + 3 nsph nthreads),    C = min(ndir, chunk)
   !>
   !> and so makes it linear in the system size instead of quadratic.
   !>
   !> The price is that each block re-traverses the grid, and with it the 16
   !> direction-independent seeds of the point prologue:
   !>
   !>     apply_seed calls per grid point  =  16 ceil(ndir / C)
   !>
   !> against the `16 ndir` a per-direction seed chain would cost and the `16`
   !> of an unblocked traversal. At the default every `ndir <= 192` -- every
   !> Hessian-vector product of a system up to 64 atoms, and every explicit
   !> direction set a caller is likely to hand in -- is a single block and pays
   !> exactly `16`.
   integer, parameter :: hvp_chunk_dirs = 192

contains

   !> Contract the moving surface adjoints against the primal map
   !>
   !> Accumulates `J^T (d omega/dv)` for the energy whose raw surface adjoints
   !> `acc` holds, one gradient column per nuclear direction. The result is
   !> *added* to `hvp`, so the caller may already hold the fixed half, and the
   !> accumulator is left untouched when anything fails.
   !>
   !> This half is the one that needs **both** forms of the adjoints. `eff` is
   !> the primal fold, produced once by the caller in `hessian.f90` and shared
   !> with the fixed half; `acc` is what that fold was built from, and pass 2
   !> below differentiates the one out of the other. Neither substitutes for the
   !> other, which is why the folding is not repeated here.
   !>
   !> Mirrors [[get_surface_gradient_drop]] throughout: same thread setup, same
   !> per-point jet and bordered KKT solve, same error latching and the same
   !> deterministic reduction. What differs is the weights -- `d(eff)` rather
   !> than `eff` -- and the direction loop those weights force inside the grid
   !> loop.
   !>
   !> @param[in]    self  DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc   Raw surface-observable adjoints, held fixed
   !> @param[in]    eff   Folded surface adjoints of the base geometry
   !> @param[in]    dirs  Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hvp   Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]   error Error object, allocated on failure
   module subroutine get_surface_hessian_response_drop(self, acc, eff, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Raw surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded surface adjoints, as [[prepare_surface_weights]] returned them
      type(drop_surface_weights_type), intent(in) :: eff
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
      !> Accumulator as it was handed in, kept only while a block can still fail
      real(wp), allocatable :: hvp_entry(:, :, :)
      !> Thread bookkeeping
      integer :: thread_slot, ithread
      !> First failure seen anywhere in the parallel region
      type(drop_abort_latch_type) :: abort

      !> Per-thread state of the grid point being opened: the jets, the seed
      !> state, the bordered factorization, the solved seed batch and the
      !> buffers they are read into
      type(drop_point_scratch_type) :: pt
      !> Whether the shared prologue cleared the point for this traversal
      logical :: point_ok

      !> Grid, direction, atom and active-slot indices
      integer :: igrid, idir, ndir, iatom, i
      !> Direction block: its size, its bounds in `dirs` and its own extent
      integer :: nchunk, ilo, ihi, nblk

      !> Linear responses of the 16 seeds, and their induced point motion
      type(drop_seed_result_type) :: res_jet(drop_n_jet_seeds)
      type(drop_seed_result_type) :: res_anchor(drop_n_anchor_seeds)
      real(wp) :: x_jet(4, drop_n_jet_seeds), x_anchor(4, drop_n_anchor_seeds)

      !> Point-local level-set adjoint weights built from the 13 field seeds
      real(wp) :: w_lsf0_pt, w_lsf1_pt(3), w_lsf2_pt(3, 3)
      !> Effective position adjoint seen by every seed; identically zero here
      real(wp) :: w_xyz_local(3)
      !> Sparse switching rows of the owner sphere
      real(wp) :: swi_owner_row(3), swi_f0, swi_dxi
      logical :: swi_live
      integer :: jj, knb

      !> Tangent of the folded surface adjoints along each direction
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
      ! No branch guard here, and that asymmetry with the public accessors is
      ! deliberate rather than an omission. They refuse a multi-branch grid
      ! because the *composite* is short a second-order branch term; this half
      ! on its own is not. The folding of the branched weights is what it
      ! differentiates, and pass 2 ([[prepare_surface_weights_tangent]]) emits
      ! `d(branch_phi_adj)` as one of its three moving channels. Guarding here
      ! would reject grids this traversal handles correctly, and its own suite
      ! finite-differences one of them.
      if (self%ngrid <= 0) return

      h_shres = self%ctx%timer%resolve("Surface Hessian (adjoint response)", &
                                       self%ctx%timer%current(), cat_gradient)
      call self%ctx%timer%start(h_shres)

      !* -------------------------------- Thread setup -------------------------------- *!
      ! Order 3, as on the gradient path: the primal map is the one being
      ! contracted, so nothing above the third jet derivative is read.
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)

      !* ------------------------------ Direction blocks ------------------------------ *!
      ! See `hvp_chunk_dirs` for the memory bound and the recompute factor. A
      ! direction set no larger than one block -- the common case -- takes a
      ! single iteration and is exactly the unblocked traversal.
      !
      ! `hvp` is `intent(inout)` and this routine owes the caller an untouched
      ! accumulator when anything fails, which a multi-block run can no longer
      ! promise by construction: the reduction lands in `hvp` block by block, so
      ! that a direction's column sees the same additions in the same order
      ! whatever the blocking. A copy taken up front is what restores the
      ! promise, and it is taken only when a second block can actually fail.
      nchunk = min(ndir, hvp_chunk_dirs)
      allocate (hvp_threads(3, self%nsph, nchunk, slots%nthreads))
      if (nchunk < ndir) allocate (hvp_entry, source=hvp)

      do ilo = 1, ndir, nchunk
         ihi = min(ilo + nchunk - 1, ndir)
         nblk = ihi - ilo + 1

         !* ---------------------- Pass 2: the moving weights ------------------------- *!
         ! Serial over this block's directions, for the group-reduction reason
         ! [[weight_tangents]] documents; `deff` is indexed `1 .. nblk`, and the
         ! global direction index appears nowhere below this line.
         call weight_tangents(self, acc, eff, dirs(:, :, ilo:ihi), deff, error)
         if (allocated(error)) then
            if (allocated(hvp_entry)) hvp = hvp_entry
            call self%ctx%timer%stop(h_shres)
            return
         end if

         hvp_threads = 0.0_wp
         call abort%reset()

         !$omp parallel num_threads(slots%nthreads) default(shared) &
         !$omp& private(thread_slot, igrid, pt, point_ok, idir, iatom, i, &
         !$omp& res_jet, res_anchor, x_jet, x_anchor, &
         !$omp& w_lsf0_pt, w_lsf1_pt, w_lsf2_pt, w_xyz_local, &
         !$omp& swi_owner_row, swi_f0, swi_dxi, swi_live, &
         !$omp& jj, knb)
         thread_slot = 1
!$       thread_slot = omp_get_thread_num() + 1

         call pt%init(self%nsph, self%iswig, want_seed_batch=.true., want_vjp=.true.)

         !$omp do schedule(static, 8)
         do igrid = 1, self%ngrid
            !* ------------------------- Shared point prologue ------------------------ *!
            ! Point, jets, seed state and the bordered right-hand sides of the
            ! jet and anchor seeds. No curvature: `d(w_k1)` and `d(w_k2)` vanish
            ! identically whatever the host put in `acc`, so this half never
            ! reads a curvature response and asks [[fill_seed_state]] for none.
            call drop_point_prologue(self, slots, thread_slot, igrid, .false., &
                                     "get_surface_hessian_response_drop", abort, &
                                     pt, point_ok)
            if (.not. point_ok) cycle

            ! Outward-normal channel: `d(w_xyz)` and `d(w_n)` are identically
            ! zero (`weights_tangent.f90`, module header -- a consumer of pass 2
            ! "gains the right to skip the normal and curvature channels of the
            ! tangent contraction outright"), so the effective position adjoint
            ! every seed sees is the zero vector and [[seed_normal_channel]] is
            ! not called.
            w_xyz_local = 0.0_wp

            !* -------------------- Direction-independent point data ------------------ *!
            ! The 16 seeds of [[seed_jet_basis]] and [[seed_anchor]], applied
            ! once for the whole direction loop. Their construction and
            ! [[apply_seed]] -- the eigen decomposition, the tangent frame, the
            ! Jacobian -- read `state` and the bordered solve alone; a weight
            ! set reaches only the contraction. So the direction loop below runs
            ! the cheap half sixteen times and the expensive one not at all,
            ! which is the shape [[get_surface_hessian_fixed_drop]] already has
            ! and the reason this traversal costs 16 applies per grid point
            ! rather than `16 ndir` of them.
            call seed_jet_basis_apply(pt%state, pt%kkt_rhs, res_jet, x_jet)
            call seed_anchor_apply(pt%state, pt%kkt_rhs, res_anchor, x_anchor)

            pt%n_active = slots%lsf(thread_slot)%lsf%active_count()
            do i = 1, pt%n_active
               pt%active_idx(i) = slots%lsf(thread_slot)%lsf%active_atom(i)
            end do

            ! The iSwiG rows depend on the geometry alone; only the scalar
            ! weight in front of them is per direction. Collected once when any
            ! direction carries a live switching tangent, so the scan stops at
            ! the first one that does.
            swi_live = .false.
            do idir = 1, nblk
               if (abs(deff(idir)%w_f(igrid)) > seed_weight_tol) then
                  swi_live = .true.
                  exit
               end if
            end do
            if (swi_live) then
               call self%iswig%swi_collect(pt%anchor, pt%owner_idx, self%anchor_xi0(igrid), &
                                           swi_f0, pt%iswig_work)
               call self%iswig%swi1_rA_sparse(pt%iswig_work, pt%swi_rows, &
                                              swi_owner_row, swi_dxi)
            end if

            !* ---------------------------- Direction loop ---------------------------- *!
            do idir = 1, nblk

               !* ------------- Field seeds -> level-set adjoint tangents ------------- *!
               w_lsf0_pt = 0.0_wp
               w_lsf1_pt = 0.0_wp
               w_lsf2_pt = 0.0_wp
               call seed_jet_basis_contract(deff(idir), igrid, pt%phi1_r, w_xyz_local, &
                                            res_jet, x_jet, &
                                            w_lsf0_pt, w_lsf1_pt, w_lsf2_pt)

               !* ----------- Field channel: contract with nuclear partials ----------- *!
               ! As on the gradient path, the level set contracts the jet
               ! indices itself, so the mixed third derivative is never
               ! materialized.
               call slots%lsf(thread_slot)%lsf%vjp_f1_rA(w_lsf0_pt, w_lsf1_pt, &
                                                         w_lsf2_pt, pt%vjp_pt)
               do i = 1, pt%n_active
                  iatom = pt%active_idx(i)
                  hvp_threads(:, iatom, idir, thread_slot) = &
                     hvp_threads(:, iatom, idir, thread_slot) + pt%vjp_pt(:, i)
               end do

               !* --------------------- Anchor channel (owner) ------------------------ *!
               call seed_anchor_contract(deff(idir), igrid, pt%phi1_r, w_xyz_local, &
                                         res_anchor, x_anchor, &
                                         hvp_threads(:, pt%owner_idx, idir, thread_slot))

               !* -------------------- iSwig switching channel ------------------------ *!
               if (abs(deff(idir)%w_f(igrid)) > seed_weight_tol) then
                  do jj = 1, pt%iswig_work%n_nb
                     knb = pt%iswig_work%idx(jj)
                     hvp_threads(:, knb, idir, thread_slot) = &
                        hvp_threads(:, knb, idir, thread_slot) &
                        + deff(idir)%w_f(igrid)*pt%swi_rows(:, jj)
                  end do
                  hvp_threads(:, pt%owner_idx, idir, thread_slot) = &
                     hvp_threads(:, pt%owner_idx, idir, thread_slot) &
                     + deff(idir)%w_f(igrid)*swi_owner_row
               end if

            end do

         end do
         !$omp end do

         call pt%destroy()
         !$omp end parallel

         if (abort%requested) then
            call abort%raise("get_surface_hessian_response_drop", error)
            if (allocated(hvp_entry)) hvp = hvp_entry
            call self%ctx%timer%stop(h_shres)
            return
         end if

         ! Deterministic reduction: fixed thread order, independent of
         ! scheduling. One block's columns at a time, and blocks are disjoint in
         ! `idir`, so a direction sees exactly the same additions in the same
         ! order however the blocking falls.
         do ithread = 1, slots%nthreads
            hvp(:, :, ilo:ihi) = hvp(:, :, ilo:ihi) + hvp_threads(:, :, 1:nblk, ithread)
         end do

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
   !> `dirs` is one *block* of the caller's direction set and both extents are
   !> taken from it, so `deff` is indexed within the block and the seven
   !> `(ngrid, ndir)` arrays this routine holds are block sized. That is where
   !> the working set of the whole traversal is bounded; see `hvp_chunk_dirs`.
   !>
   !> @param[in]  self  DROP cavity instance
   !> @param[in]  acc   Raw surface adjoints, held fixed
   !> @param[in]  eff   Folded weights, as [[prepare_surface_weights]] returned them
   !> @param[in]  dirs  Nuclear directions of one block, `(3, nsph, ndir)`
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
