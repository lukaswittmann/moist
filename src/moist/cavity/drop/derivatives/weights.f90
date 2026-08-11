!> Preparation shared by all DROP reverse-mode implementations
!>
!> For now, these are:
!> * `potential.f90`, the electronic path, which contracts a surface adjoint
!>   into level-set adjoint weights (cavity response fock), and
!> * `nuclear.f90`, the nuclear path, which contracts a surface adjoint
!>   against nuclear partial derivatives
!>
!> They share: validation of the accumulator, folding of the derived weight
!> channels, running the branch-softmax reverse pass, and copying the cavity's
!> grid-level scalars into the kernel's seed state
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_weights
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type, &
      & drop_seed_state_type, compute_branch_phi_adj, seed_weight_tol
   implicit none (type, external)

contains

   !> Reject a surface adjoint the cavity cannot contract
   !>
   !> Every message is prefixed with the caller's name, so a failure names the
   !> entry point the user actually called.
   !>
   !> @param[in]  self    DROP cavity instance
   !> @param[in]  acc     Accumulated surface-observable adjoints
   !> @param[in]  context Calling routine, used to prefix the diagnostics
   !> @param[out] error   Error object, allocated on failure
   module subroutine check_surface_adjoint(self, acc, context, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

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
      ! The area and integration-weight channels are converted through 1/xi^2.
      ! A vanishing width with a live area or weight adjoint has no finite
      ! conversion, and must not be silently dropped.
      if (any(abs(self%xi0) <= seed_weight_tol .and. &
              (abs(acc%w_a) > seed_weight_tol .or. abs(acc%w_w) > seed_weight_tol))) then
         call fatal_error(error, context//": singular derived-weight conversion")
         return
      end if
   end subroutine check_surface_adjoint

   !> Fold the derived weight channels and run the branch reverse pass
   !>
   !> `fold_switching` is the *only* difference between the electronic and the
   !> nuclear preparation. The area channel folds into the switching channel
   !> through `da/df = R_I^2 wleb_i`; a parameter that leaves `f` fixed -- an
   !> electronic degree of freedom -- may skip that term, a nuclear
   !> displacement may not.
   !>
   !> Call [[check_surface_adjoint]] first; this routine assumes a valid
   !> accumulator and does not re-validate.
   !>
   !> @param[in]  self           DROP cavity instance
   !> @param[in]  acc            Accumulated surface-observable adjoints
   !> @param[in]  fold_switching Whether to fold the area channel into `w_f`
   !> @param[out] eff            Folded weights and the branch adjoint
   module subroutine prepare_surface_weights(self, acc, fold_switching, eff)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Whether the area channel also folds into the switching channel
      logical, intent(in) :: fold_switching
      !> Folded weights
      type(drop_surface_weights_type), intent(out) :: eff

      !> Grid index
      integer :: igrid
      !> Softmax temperature of the branch weights
      real(wp) :: sigma_phi

      allocate (eff%w_xi, source=acc%w_xi)
      allocate (eff%w_f, source=acc%w_f)
      allocate (eff%w_xyz, source=acc%w_xyz)
      allocate (eff%w_n, source=acc%w_n)
      allocate (eff%w_k1, source=acc%w_k1)
      allocate (eff%w_k2, source=acc%w_k2)

      ! a_i = c*f_i/xi_i^2 and w_i = c/xi_i^2, so da/dxi = -2a/xi and
      ! dw/dxi = -2w/xi; both land on the width channel.
      where (abs(acc%w_a) > seed_weight_tol)
         eff%w_xi = eff%w_xi - 2.0_wp*self%a*acc%w_a/self%xi0
      end where
      where (abs(acc%w_w) > seed_weight_tol)
         eff%w_xi = eff%w_xi - 2.0_wp*self%wleb*acc%w_w/self%xi0
      end where

      if (fold_switching) then
         do igrid = 1, self%ngrid
            if (abs(acc%w_a(igrid)) > seed_weight_tol) then
               eff%w_f(igrid) = eff%w_f(igrid) &
                                + acc%w_a(igrid)*self%radii(self%owner(igrid))**2*self%wleb(igrid)
            end if
         end do
      end if

      eff%have_wn = any(abs(eff%w_n) > seed_weight_tol)
      eff%have_wk = any(abs(eff%w_k1) > seed_weight_tol) &
                    .or. any(abs(eff%w_k2) > seed_weight_tol)

      ! Reverse pass for the branch-weight post-pass: converts the width-induced
      ! adjoint dL/dp_m into dL/dPhi_m, which the seed loop couples to dr/dp.
      ! Runs after the folds, since it reads the folded width channel.
      allocate (eff%branch_phi_adj(self%ngrid), source=0.0_wp)
      if (self%ngrid > 0 .and. allocated(self%branch_count)) then
         sigma_phi = self%branch_weight%s
         call compute_branch_phi_adj(self%branch_count(1:self%ngrid), &
                                     self%anchor_id(1:self%ngrid), &
                                     self%wbranch(1:self%ngrid), &
                                     self%wleb(1:self%ngrid), &
                                     self%xi0(1:self%ngrid), &
                                     sigma_phi, eff%w_xi, eff%branch_phi_adj)
      end if
   end subroutine prepare_surface_weights

   !> Copy the cavity's grid-level scalars into a kernel seed state
   !>
   !> Only the fields the cavity owns are written. The level-set jet
   !> (`lsf1_r`, `lsf2_rr`, `lsf3_rrr`) and the multiplier come from the
   !> caller's own LSF evaluation and are set by the caller.
   !>
   !> @param[in]    self           DROP cavity instance
   !> @param[in]    igrid          Grid point to describe
   !> @param[in]    want_curvature Whether the curvature invariants are needed
   !> @param[inout] state          Seed state; its cavity-owned inputs are set
   module subroutine fill_seed_state(self, igrid, want_curvature, state)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Grid point
      integer, intent(in) :: igrid
      !> Whether the curvature invariants are needed
      logical, intent(in) :: want_curvature
      !> Seed state
      type(drop_seed_state_type), intent(inout) :: state

      state%alpha_coeff = self%param%phi_alpha
      state%anchor = self%anchorxyz(:, igrid)
      state%owner_xyz = self%mol%xyz(:, self%owner(igrid))
      state%anchor_wleb0 = self%anchor_wleb0(igrid)
      state%cpjac_scal0 = self%cpjac_scal0(igrid)
      state%w_f0 = self%w_f0(igrid)
      state%wbranch = self%wbranch(igrid)
      state%wleb = self%wleb(igrid)
      state%xi0 = self%xi0(igrid)
      state%want_curvature = want_curvature
   end subroutine fill_seed_state

end submodule moist_cavity_drop_derivatives_weights
