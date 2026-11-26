!> Preparation shared by all DROP reverse-mode implementations
!>
!> For now, these are:
!> * `potential.f90`, the electronic path, which contracts a surface adjoint
!>   into level-set adjoint weights (cavity response fock), and
!> * `nuclear.f90`, the nuclear path, which contracts a surface adjoint
!>   against nuclear partial derivatives
!>
!> They share: validation of the accumulator, folding of the derived weight
!> channels, running the branch-weight reverse pass, and copying the cavity's
!> grid-level scalars into the kernel's seed state
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_weights
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type, &
      & drop_seed_state_type, compute_branch_adj, seed_weight_tol
   implicit none(type, external)

contains

   !> Reject a surface adjoint the cavity cannot contract
   !>
   !> Every message is prefixed with the caller's name, so a failure names the
   !> entry point the user actually called
   !>
   !> @param[in]  self    DROP cavity instance
   !> @param[in]  acc     Accumulated surface-observable adjoints
   !> @param[in]  context Calling routine, used to prefix the diagnostics
   !> @param[out] error   Error object, allocated on failure
   module subroutine check_surface_adjoint(self, acc, context, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_adjoint_type), intent(in) :: acc
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_built(error)
      if (allocated(error)) then
         error%message = context//": "//error%message
         return
      end if
      if (.not. allocated(self%lsf_model)) then
         call fatal_error(error, context//": level-set model is not initialized")
         return
      end if
      if (.not. acc%is_initialized()) then
         call fatal_error(error, context//": accumulator is not initialized")
         return
      end if
      if (size(acc%w_xi) /= self%ngrid) then
         call fatal_error(error, context//": accumulator grid size mismatch")
         return
      end if
      if (.not. allocated(self%xi0) .or. .not. allocated(self%a) .or. &
          .not. allocated(self%wleb) .or. .not. allocated(self%f)) then
         call fatal_error(error, context//": cavity surface data are incomplete")
         return
      end if
      ! The area and integration-weight channels are converted through 1/xi^2
      ! A vanishing width with a live area or weight adjoint has no finite
      ! conversion, and must not be silently dropped
      if (any(abs(self%xi0) <= seed_weight_tol .and. &
              (abs(acc%w_a) > seed_weight_tol .or. abs(acc%w_w) > seed_weight_tol))) then
         call fatal_error(error, context//": singular derived-weight conversion")
         return
      end if
   end subroutine check_surface_adjoint

   !> Fold the derived weight channels and run the branch reverse pass
   !>
   !> `fold_switching` is the *only* difference between the electronic and the
   !> nuclear preparation
   !>
   !> - the area channel folds into the switching channel through
   !>   `da/df = R_I^2 wleb_i`
   !> - a parameter leaving `f` fixed, an electronic degree of freedom, may skip
   !>   that term; a nuclear displacement may not
   !>
   !> Call [[check_surface_adjoint]] first; this routine assumes a valid
   !> accumulator and does not re-validate
   !>
   !> @param[in]  self           DROP cavity instance
   !> @param[in]  acc            Accumulated surface-observable adjoints
   !> @param[in]  fold_switching Whether to fold the area channel into `w_f`
   !> @param[out] eff            Folded weights and the branch adjoint
   module subroutine prepare_surface_weights(self, acc, fold_switching, eff)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_adjoint_type), intent(in) :: acc
      !> Whether the area channel also folds into the switching channel
      logical, intent(in) :: fold_switching
      !> Folded weights
      type(drop_surface_weights_type), intent(out) :: eff

      !> Grid index and grid size
      integer :: igrid, n
      !> Team size, anchor count, anchor block, its first/last anchor and rows
      integer :: nthreads, nanchor, nblock, iblock, ia0, ia1, r0, r1
      !> Anchors per block and the number of blocks holding one more
      integer :: nper, nextra
      !> Adjoint of the switching value with the area channel folded in
      real(wp), allocatable :: w_f_area(:)
      !> Switching adjoint with the area fold of one point
      real(wp) :: w_f_folded
      !> Whether the branch reverse pass runs, and the channel flags
      logical :: have_branch, have_wn, have_wk

      ! One parallel pass over the grid, every array first touched by the thread
      ! that fills it: the copies and folds are elementwise, so the result does
      ! not depend on the team size
      n = self%ngrid
      nthreads = self%ctx%get_num_threads()
      allocate (eff%w_xi(n), eff%w_f(n), eff%w_xyz(3, n), eff%w_n(3, n), eff%w_k1(n), eff%w_k2(n))
      allocate (eff%branch_phi_adj(n), eff%branch_gate_adj(n))
      have_branch = n > 0 .and. allocated(self%offset)
      if (have_branch) allocate (w_f_area(n))
      have_wn = .false.
      have_wk = .false.

      !$omp parallel do num_threads(nthreads) default(shared) private(igrid, w_f_folded) &
      !$omp& schedule(static) reduction(.or.:have_wn, have_wk)
      do igrid = 1, n
         eff%w_xyz(:, igrid) = acc%w_xyz(:, igrid)
         eff%w_n(:, igrid) = acc%w_n(:, igrid)
         eff%w_k1(igrid) = acc%w_k1(igrid)
         eff%w_k2(igrid) = acc%w_k2(igrid)
         have_wn = have_wn .or. any(abs(acc%w_n(:, igrid)) > seed_weight_tol)
         have_wk = have_wk .or. abs(acc%w_k1(igrid)) > seed_weight_tol &
                   .or. abs(acc%w_k2(igrid)) > seed_weight_tol

         ! a_i = c*f_i/xi_i^2 and w_i = c/xi_i^2, so da/dxi = -2a/xi and
         ! dw/dxi = -2w/xi; both land on the width channel
         eff%w_xi(igrid) = acc%w_xi(igrid)
         if (abs(acc%w_a(igrid)) > seed_weight_tol) then
            eff%w_xi(igrid) = eff%w_xi(igrid) &
                              - 2.0_wp*self%a(igrid)*acc%w_a(igrid)/self%xi0(igrid)
         end if
         if (abs(acc%w_w(igrid)) > seed_weight_tol) then
            eff%w_xi(igrid) = eff%w_xi(igrid) &
                              - 2.0_wp*self%wleb(igrid)*acc%w_w(igrid)/self%xi0(igrid)
         end if

         ! The area channel always reaches the branch weights p through
         ! a = c*f*wleb; it folds into `w_f` itself only with `fold_switching`
         w_f_folded = acc%w_f(igrid)
         if (abs(acc%w_a(igrid)) > seed_weight_tol) then
            w_f_folded = w_f_folded + acc%w_a(igrid)*self%radii(self%owner(igrid))**2*self%wleb(igrid)
         end if
         if (fold_switching) then
            eff%w_f(igrid) = w_f_folded
         else
            eff%w_f(igrid) = acc%w_f(igrid)
         end if
         if (have_branch) w_f_area(igrid) = w_f_folded
      end do
      !$omp end parallel do
      eff%have_wn = have_wn
      eff%have_wk = have_wk

      ! Reverse pass for the branch-weight post-pass: f = f_anchor*p, so the
      ! adjoint of f gives dL/dp_m and from it dL/dPhi_l and dL/dg_l of every
      ! row of the anchor, which the seed loops couple to the point motion and
      ! the gate response
      ! Anchors are independent: blocks of anchors run as rebased CSR slices,
      ! each writing only its own rows
      if (.not. have_branch) then
         eff%branch_phi_adj = 0.0_wp
         eff%branch_gate_adj = 0.0_wp
         return
      end if
      ! Several blocks per thread, handed out dynamically, as the rows per anchor vary;
      ! the split by quotient and remainder cannot overflow
      nanchor = size(self%offset) - 1
      nblock = 1
      if (nthreads > 1) nblock = max(1, min(nanchor, 4*nthreads))
      nper = nanchor/nblock
      nextra = mod(nanchor, nblock)
      !$omp parallel do num_threads(nthreads) default(shared) private(iblock, ia0, ia1, r0, r1) &
      !$omp& schedule(dynamic, 1)
      do iblock = 1, nblock
         ia0 = (iblock - 1)*nper + min(iblock - 1, nextra) + 1
         ia1 = iblock*nper + min(iblock, nextra)
         r0 = self%offset(ia0) + 1
         r1 = self%offset(ia1 + 1)
         call compute_branch_adj(self%offset(ia0:ia1 + 1) - self%offset(ia0), self%phi0(r0:r1), &
                                 self%w_f0(r0:r1), self%anchor%f0(ia0:ia1), self%branch_weight, &
                                 w_f_area(r0:r1), eff%branch_phi_adj(r0:r1), eff%branch_gate_adj(r0:r1))
      end do
      !$omp end parallel do
   end subroutine prepare_surface_weights

   !> Copy the cavity's grid-level scalars into a kernel seed state
   !>
   !> Only the fields the cavity owns are written
   !>
   !> - the level-set jet (`lsf1_r`, `lsf2_rr`, `lsf3_rrr`) and the multiplier
   !>   come from the caller's own LSF evaluation and are set by the caller
   !> - the reach gate and its partial in the projected point go into the state;
   !>   its partial in the anchor, which moves with the owner alone, is returned
   !>   when asked for
   !>
   !> @param[in]    self           DROP cavity instance
   !> @param[in]    igrid          Grid point to describe
   !> @param[in]    want_curvature Whether the curvature invariants are needed
   !> @param[inout] state          Seed state; its cavity-owned inputs are set
   !> @param[out]   reach_da       Reach-gate partial in the anchor
   module subroutine fill_seed_state(self, igrid, want_curvature, state, reach_da)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Grid point
      integer, intent(in) :: igrid
      !> Whether the curvature invariants are needed
      logical, intent(in) :: want_curvature
      !> Seed state
      type(drop_seed_state_type), intent(inout) :: state
      !> Reach-gate partial in the anchor
      real(wp), intent(out), optional :: reach_da(3)

      !> Path from the anchor to the point, its length, the gate's slope in it
      real(wp) :: v(3), pdist, fp

      state%alpha_coeff = self%param%phi_alpha
      state%anchor = self%anchor%xyz(:, self%anchor_of_row(igrid))
      state%owner_xyz = self%mol%xyz(:, self%owner(igrid))
      state%anchor_wleb0 = self%anchor%wleb0(self%anchor_of_row(igrid))
      state%cpjac_scal0 = self%cpjac_scal0(igrid)
      state%wleb = self%wleb(igrid)
      state%xi0 = self%xi0(igrid)
      state%want_curvature = want_curvature

      ! A closed gate stays closed to first order, every factor at zero being flat;
      ! the reach gate is flat wherever the path length has its cusp
      state%f_reach_f0 = 0.0_wp
      state%f_reach_dy = 0.0_wp
      if (self%w_f0(igrid) > 0.0_wp) then
         v = self%xyz(:, igrid) - state%anchor
         pdist = sqrt(dot_product(v, v))
         call self%f_reach%eval(pdist, state%f_reach_f0, fp)
         if (fp /= 0.0_wp) state%f_reach_dy = fp*v/pdist
      end if
      if (present(reach_da)) reach_da = -state%f_reach_dy
   end subroutine fill_seed_state

end submodule moist_cavity_drop_derivatives_weights
