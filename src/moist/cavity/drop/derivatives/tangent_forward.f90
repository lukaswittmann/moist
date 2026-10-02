!> Pass 1 of the DROP Hessian: forward tangent of the surface map
!>
!> - Forward-mode dual of [[get_surface_gradient_drop]]: same thread setup,
!>   point prologue, KKT factorization and abort latch
!> - One seed per nuclear direction, its image in the level-set jet
!> - Outputs `d_a, d_wleb, d_xi0, d_wbranch`, each `(ngrid, ndir)`
!> - Point `i` writes only rows `(i, :)`: no thread buffers, no ordered sum
!> - Stage 1, grid loop (parallel): jet tangents, batched KKT solve for
!>   `(dr, dlambda)`, [[apply_seed]], iSwiG rows for `d(f)`, `d(Phi)`
!> - Jet tangents from [[drop_field_jet_point]] and
!>   [[drop_field_jet_tangent]], or from `tangent_jet` per direction
!> - Jet tensors are active-slot indexed; directions are gathered per point
!> - Stage 2, branch softmax (serial): one
!>   [[branch_weight_type:weights_grad]] call per anchor group, `nparam = ndir`
!> - Stage 2 stays serial: a softmax group is never split across threads
!> - Groups are runs of equal `anchor_id`, contiguous by the stable
!>   `counting_argsort` in `projection.f90`
!> - Stage 3, assembly: `d_wleb = res%dwleb + d_wbranch * wleb / wbranch`,
!>   then `d_xi0` and `d_a` from the completed `d_wleb`
!> - wbranch trap: [[apply_seed]] freezes `wbranch`; reverse mode carries its
!>   motion via [[compute_branch_phi_adj]], forward mode adds it in stage 3
!> - Trap shows on multi-branch anchors only; [[apply_seed_tangent]] takes
!>   the same term as `dinp_v%dwbranch`
!> - `xi0 = swx/(R sqrt(wleb))` (`compute_gaussians`) with the final `wleb`, so
!>   `d_xi0 = -0.5 * xi0 * d_wleb / wleb`; branch-frozen `res%dxi` is discarded
!> - `a = R^2 f wleb` (`compute_area_volume`), so
!>   `d_a = R^2 (wleb d_f + f d_wleb)`
!> - `f` is the iSwiG `self%f`, not `res%dw_f` (`w_f = f_crit * f_foc`, in `wleb`)
!> - `d_wbranch` is reported and also folded into `d_wleb`:
!>   [[branch_point_adjoint]] needs both, with `d_wleb` complete
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_tangent_forward
!$ use omp_lib, only: omp_get_thread_num
   use moist_cavity_drop_threads, only: drop_worker_slots_type, drop_abort_latch_type, &
      & drop_point_scratch_type
   use moist_cavity_drop_derivatives_kernel, only: drop_seed_result_type, apply_seed, &
      & seed_weight_tol, next_branch_group, max_branch_group_size, drop_n_host_jet, &
      & add_host_jet
   use moist_cavity_drop_derivatives_field_tangent, only: drop_field_tangent_work_type, &
      & drop_field_jet_point, drop_field_jet_tangent
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !> Forward tangent of the four grid scalars the weight fold reads
   !>
   !> - DROP-internal form of pass 1, wrapper over [[drop_surface_tangent_core]]
   !> - Keeps the branch-weight tangent, which [[cavity_surface_tangent_type]]
   !>   has no channel for
   !> - Outputs fully written; column `idir` is the tangent along
   !>   `dirs(:, :, idir)`
   !>
   !> @param[in]  self       DROP cavity holding a projected grid
   !> @param[in]  dirs       nuclear directions `(3, nsph, ndir)`
   !> @param[out] d_a        area element tangent `(ngrid, ndir)`
   !> @param[out] d_wleb     Lebedev weight tangent `(ngrid, ndir)`
   !> @param[out] d_xi0      Gaussian width tangent `(ngrid, ndir)`
   !> @param[out] d_wbranch  branch weight tangent `(ngrid, ndir)`
   !> @param[out] error      allocated on failure
   !> @param[in]  contracted jet tangents from the contracted accessor per
   !>                        direction, default `.false.`; see the core
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
      !> Jet tangents from the contracted accessor
      logical, intent(in), optional :: contracted

      !> Full tangent, three channels copied out
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

   !> Forward tangent of every surface observable along nuclear directions
   !>
   !> @param[in]     self    DROP cavity holding a projected grid
   !> @param[in]     dirs    nuclear directions `(3, nsph, ndir)`
   !> @param[in,out] tangent surface tangent, initialised for `(ngrid, ndir)`
   !> @param[out]    error   allocated on failure
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

   !> Forward tangent of the DROP surface map along nuclear directions
   !>
   !> - Every channel of `tangent` fully written; column `idir` is the tangent
   !>   along `dirs(:, :, idir)`
   !> - Curvature channels filled only with `want_curvature`, else zero;
   !>   `tangent%have_curvature` records which
   !> - `d_wbranch`: softmax branch-weight tangent, read by pass 2; the generic
   !>   tangent has no channel for it
   !> - `contracted` is the caller's choice and never a function of `ndir`: the
   !>   Hessian traversal passes blocks of directions, and a tangent must not
   !>   depend on the blocking
   !> - Contracted path suits few directions; only the traversal's
   !>   per-direction mode requests it
   !>
   !> @param[in]     self           DROP cavity holding a projected grid
   !> @param[in]     dirs           nuclear directions `(3, nsph, ndir)`
   !> @param[in]     want_curvature fill the curvature channels
   !> @param[in,out] tangent        surface tangent, initialised for `(ngrid, ndir)`
   !> @param[out]    error          allocated on failure
   !> @param[out]    d_wbranch      branch weight tangent `(ngrid, ndir)`, optional
   !> @param[in]     contracted     jet tangents from the contracted accessor per
   !>                               direction instead of per-point tensors,
   !>                               default `.false.`
   !> @param[in]     host_jets      host partial jet tangents at the fixed points,
   !>                               `(40, ngrid, ndir)` packed by spatial order,
   !>                               added to the level set's own; optional
   module subroutine drop_surface_tangent_core(self, dirs, want_curvature, tangent, &
                                               error, d_wbranch, contracted, host_jets)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Fill the curvature channels
      logical, intent(in) :: want_curvature
      !> Surface tangent
      type(cavity_surface_tangent_type), intent(inout) :: tangent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Tangent of the softmax branch weight
      real(wp), intent(out), optional :: d_wbranch(:, :)
      !> Jet tangents from the contracted accessor
      logical, intent(in), optional :: contracted
      !> Host partial jet tangents along every direction
      real(wp), intent(in), optional :: host_jets(:, :, :)

      !> Host jet tangents folded in
      logical :: have_jets

      !> Branch-weight tangent, always kept: the assembly reads it at every point
      real(wp), allocatable :: dwb(:, :)

      !> Per-thread level-set clones and objectives
      type(drop_worker_slots_type) :: slots
      !> Thread bookkeeping
      integer :: thread_slot
      !> First failure in the parallel region
      type(drop_abort_latch_type) :: abort
      !> Per-thread failure, handed to the latch
      type(error_type), allocatable :: worker_error

      !> Per-thread point state: jets, seed state, bordered factorization, buffers
      type(drop_point_scratch_type) :: pt
      !> Linear response of one seed
      type(drop_seed_result_type) :: res
      !> Point cleared by the shared prologue
      logical :: point_ok

      !> Grid, direction, sphere, active-slot and Cartesian indices
      integer :: igrid, idir, ndir, jj, knb, i
      !> Jet tensors materialised per point
      logical :: materialise

      !> Mixed nuclear tensors of the level set at the point
      type(drop_field_tangent_work_type) :: ft_work
      !> One direction gathered onto the point's active slots
      real(wp), allocatable :: v_act(:, :)

      !> Directional nuclear tangents of the level-set jet at the fixed point
      real(wp) :: dlsf0
      real(wp), allocatable :: dlsf1_r(:, :), dlsf2_rr(:, :, :)
      !> Bordered KKT right-hand sides, one column per direction
      real(wp), allocatable :: dir_rhs(:, :)
      !> Induced motion of the projected point and of the multiplier
      real(wp) :: dr(3), dlambda

      !> Sparse switching rows of the owner sphere
      real(wp) :: swi_owner_row(3), swi_f0, swi_dxi, df_dir

      !> Branch objective tangent `d(Phi)` `(ngrid, ndir)`
      real(wp), allocatable :: dphi(:, :)
      !> Softmax scratch of the branch stage
      real(wp), allocatable :: branch_phi(:), branch_dphi(:, :)
      real(wp), allocatable :: branch_weights(:), branch_dweights(:, :)
      !> Branch group bookkeeping
      integer :: igroup_cursor, igroup_start, igroup_end, group_size
      integer :: m_branch, im_grid, nbranch_max
      !> Further anchor group exists
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

      call tangent%zero()
      tangent%have_curvature = want_curvature
      if (self%ngrid <= 0) return
      allocate (dwb(self%ngrid, ndir), source=0.0_wp)

      h_stan = self%ctx%timer%resolve("Surface tangent", self%ctx%timer%current(), &
                                      cat_gradient)
      call self%ctx%timer%start(h_stan)

      !* -------------------------------- Thread setup -------------------------------- *!
      call slots%init(self%ctx, self%lsf_model, 3, self%param, self%mol, self%radii)
      allocate (dphi(self%ngrid, ndir), source=0.0_wp)
      ! Per-point tensors by default, contracted accessor on request
      ! Never keyed on `ndir`; see the procedure doc
      materialise = .true.
      if (present(contracted)) materialise = .not. contracted
      have_jets = .false.
      if (present(host_jets)) then
         if (any(shape(host_jets) /= [drop_n_host_jet, self%ngrid, ndir])) then
            call fatal_error(error, "get_surface_tangent_drop: host jet tangents must be"// &
                             " (40, ngrid, ndir)")
            call self%ctx%timer%stop(h_stan)
            return
         end if
         have_jets = .true.
      end if

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
         ! Point, jets, seed state, bordered factorization; failures latch themselves
         ! No 16-seed batch: one right-hand side per direction, built below
         call drop_point_prologue(self, slots, thread_slot, igrid, want_curvature, &
                                  "get_surface_tangent_drop", abort, pt, point_ok)
         if (.not. point_ok) cycle

         !* ----------------- Directional nuclear tangents of the jet ----------------- *!
         ! `sum_B v_B . d(jet)/dR_B` at the fixed point, per direction
         ! Tensor fill is unconditional: the buffer outlives the point
         ! Directions are gathered onto the prologue's active slots
         if (materialise) call drop_field_jet_point(slots%lsf(thread_slot)%lsf, ft_work)
         do idir = 1, ndir
            if (materialise) then
               do i = 1, pt%n_active
                  v_act(:, i) = dirs(:, pt%active_idx(i), idir)
               end do
               call drop_field_jet_tangent(ft_work, pt%n_active, v_act, dlsf0, &
                                           dlsf1_r(:, idir), dlsf2_rr(:, :, idir))
            else if (pt%n_active > 0) then
               call slots%lsf(thread_slot)%lsf%tangent_jet(dirs(:, :, idir), dlsf0, &
                                                            dlsf1_r(:, idir), dlsf2_rr(:, :, idir))
            else
               ! No active atoms: no nuclear partials, and the contracted
               ! accessor is the erroring default
               dlsf0 = 0.0_wp
               dlsf1_r(:, idir) = 0.0_wp
               dlsf2_rr(:, :, idir) = 0.0_wp
            end if
            ! Host partial tangents add to the level set's own; the whole
            ! tangent for a host-defined level set
            if (have_jets) then
               call add_host_jet(host_jets(:, igrid, idir), dlsf0, dlsf1_r(:, idir), &
                                 dlsf2_rr(:, :, idir))
            end if

            ! Bordered right-hand side, with `d^2 phi/(dr dR_owner) = -alpha*I`
            ! (`objective_phi.f90`, `f2_r_rA`); anchor rigid with its owner
            dir_rhs(1:3, idir) = self%param%phi_alpha*dirs(:, pt%owner_idx, idir) &
                                 + pt%state%lambda_val*dlsf1_r(:, idir)
            dir_rhs(4, idir) = -dlsf0
         end do

         !* ------------------------ Bordered KKT sensitivities ----------------------- *!
         ! Batched solve on the prologue's factorization; one 4x4 matrix per point
         call pt%kkt_fac%solve(dir_rhs, "get_surface_tangent_drop", worker_error, igrid)
         if (allocated(worker_error)) then
            call abort%latch_error(worker_error, igrid)
            cycle
         end if

         !* ---------------------- Base Lebedev-weight response ----------------------- *!
         do idir = 1, ndir
            dr = dir_rhs(1:3, idir)
            dlambda = dir_rhs(4, idir)

            ! Projected point moves with the bordered solve
            tangent%d_xyz(:, igrid, idir) = dr

            call apply_seed(pt%state, dlsf1_r(:, idir), dlsf2_rr(:, :, idir), dr, dlambda, res)

            ! Branch-frozen half of `d(wleb)`, completed in stage 3
            ! `res%dxi` belongs to this incomplete weight and is not read
            tangent%d_w(igrid, idir) = res%dwleb

            ! Outward normal `grad S/|grad S|` at the moving point; principal
            ! curvatures when the seed state carries them
            tangent%d_n(:, igrid, idir) = res%dn_surf
            if (want_curvature) then
               tangent%d_k1(igrid, idir) = res%dk1
               tangent%d_k2(igrid, idir) = res%dk2
            end if

            ! Tangent of `Phi = 0.5 alpha |r* - anchor|^2`: point moves by `dr`,
            ! anchor rigidly with its owner
            dphi(igrid, idir) = dot_product(pt%phi1_r, dr - dirs(:, pt%owner_idx, idir))
         end do

         !* ------------------------- iSwiG switching channel ------------------------- *!
         ! `f` at the anchor with the anchor width
         ! `anchor_xi0` depends on owner radius and raw Lebedev weight only:
         ! no nuclear tangent, `swi_dxi` unused
         call self%iswig%swi_collect(pt%state%anchor, pt%owner_idx, self%anchor_xi0(igrid), &
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
      ! Points outside a multi-branch group keep `d_wbranch = 0`, exact since
      ! their `wbranch` is the constant one
      call branch_stage()

      !* -------------------------------- Assembly ------------------------------------ *!
      do idir = 1, ndir
         do igrid = 1, self%ngrid
            wleb_i = self%wleb(igrid)
            wbranch_i = self%wbranch(igrid)
            r_own = self%radii(self%owner(igrid))

            ! wbranch trap: `apply_seed` froze the branch weight, add its motion
            ! `wleb/wbranch` is the pre-branch weight the softmax multiplies,
            ! formed as in `forward.f90` and `compute_branch_phi_adj`
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
      !> - Serial only; see the module header
      !> - Group walk by [[next_branch_group]], as in [[compute_branch_phi_adj]]
      !>   and the branch post-pass of `forward.f90`
      !> - Groups are runs of equal `anchor_id` starting at `branch_count > 1`
      !> - Softmax derivatives in `(nparam, nbranch)` layout; `nparam = ndir`
      !>   gives every direction of a group in one call
      !> - Scratch sized by [[max_branch_group_size]], not
      !>   `maxval(branch_count)`: it is indexed by the group run length
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
