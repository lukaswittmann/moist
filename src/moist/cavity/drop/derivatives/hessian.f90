!> Public DROP Hessian accessors
!>
!> The nuclear gradient of `nuclear.f90` is `J^T omega`, and its directional
!> derivative splits into
!>
!>     d/dv [ J^T omega ]  =  (dJ^T/dv) omega  +  J^T (d omega/dv)
!>                            ^ hessian_fixed     ^ hessian_response
!>
!> Both halves are implemented and separately verified next door. This
!> submodule owns nothing of the mathematics; it owns the *composition* and the
!> two public entry points that composition makes reachable:
!>
!>   * [[get_surface_hessian_drop]], the Hessian-vector product, which is the
!>     primitive -- the response half is intrinsically per direction, so a
!>     direction set is what it must be asked for;
!>   * [[get_hessian_drop]], the dense `(3, nsph, 3, nsph)` block, which is a
!>     wrapper over the same two calls with the `3 nsph` Cartesian unit
!>     directions supplied **in one batch**. One batched call, not `3 nsph`
!>     serial ones: the response half rebuilds the whole per-point primal map
!>     once and loops directions inside it, so the grid traversal is paid once
!>     either way while `3 nsph` separate calls would pay it `3 nsph` times.
!>     The fixed half is direction free and already rank 4, so the dense path
!>     takes its block directly rather than contracting it against unit vectors
!>     and re-expanding the result.
!>
!>
!> ## The frozen surrogate, and why it is an identity rather than a workaround
!>
!> [[check_frozen_weights]] refuses any accumulator carrying a live `w_a` or
!> `w_w`: those channels fold into `w_xi` and `w_f` through `a`, `wleb`, `xi0`
!> and the radii, all geometry dependent, and differentiating that fold is the
!> response half's job -- so a fixed half handed such an accumulator would
!> silently omit a term. But an accumulator with live `w_a`/`w_w` is precisely
!> the one the response half exists for, and the composite needs both halves
!> driven by the same energy.
!>
!> The resolution is that the fixed half reads its accumulator through exactly
!> one door, [[prepare_surface_weights]]. Write the folded weights at the base
!> geometry as `E = eff(R_0, acc)`. [[freeze_surface_adjoint]] builds
!>
!>     acc_frozen:  w_xi := E%w_xi,  w_f := E%w_f,  w_a := 0,  w_w := 0
!>
!> with `w_xyz`, `w_n`, `w_k1` and `w_k2` copied unchanged. With both derived
!> channels zero, [[prepare_surface_weights]] folds nothing, so it returns
!> `E` again -- bit for bit, since the folds are guarded on
!> `|w_a| > seed_weight_tol` and are therefore not merely small but skipped.
!> The fixed half consequently sees the same effective weights it would have
!> seen from `acc`, and the substitution is an identity, not an approximation.
!> The guard is left exactly as it stands.
!>
!> `branch_phi_adj` is re-derived rather than copied, and it agrees for the
!> same reason: [[compute_branch_phi_adj]] is a function of the folded `w_xi`
!> and of grid quantities, and the folded `w_xi` is unchanged. It is moot in
!> practice, because [[check_frozen_weights]] also refuses a multi-branch grid
!> -- that refusal propagates out of both entry points unchanged, and is the
!> honest answer: on a branched grid the fixed half would be missing the
!> second-order branch term.
!>
!> The response half is driven by the **raw** `acc`, not by the surrogate: it
!> is the fold itself that moves, so freezing it there would zero the very term
!> it computes ([[test_frozen_response_is_zero]] asserts exactly that).
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_hessian
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type
   implicit none(type, external)

   !> Cartesian dimension
   integer, parameter :: ndim = 3

contains

   !* ================================================================================= *!
   !*                            Public Hessian accessors                               *!
   !* ================================================================================= *!

   !> Hessian-vector products of the DROP surface contribution
   !>
   !> Accumulates `d/dv [ J^T omega ]` for the energy whose surface adjoints
   !> `acc` holds, one gradient column per supplied nuclear direction. The
   !> result is *added* to `hvp`, and the accumulator is left untouched when
   !> anything fails -- both halves are formed in local buffers first, so a
   !> failure in the second one cannot leave the first one behind.
   !>
   !> @param[in]    self  DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc   Accumulated surface-observable adjoints
   !> @param[in]    dirs  Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hvp   Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]   error Error object, allocated on failure
   module subroutine get_surface_hessian_drop(self, acc, dirs, hvp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Hessian-vector accumulator
      real(wp), intent(inout) :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Direction-free fixed half, and the response half of every direction
      real(wp), allocatable :: hess_fixed(:, :, :, :), total(:, :, :)
      !> Extents and loop indices
      integer :: ndir, idir, iatom, iaxis

      !* ------------------------------- Shape guards --------------------------------- *!
      call check_surface_adjoint(self, acc, "get_surface_hessian_drop", error)
      if (allocated(error)) return
      if (size(dirs, 1) /= ndim .or. size(dirs, 2) /= self%nsph) then
         call fatal_error(error, "get_surface_hessian_drop: dirs must be (3, nsph, ndir)")
         return
      end if
      ndir = size(dirs, 3)
      if (ndir <= 0) then
         call fatal_error(error, "get_surface_hessian_drop: no direction supplied")
         return
      end if
      if (size(hvp, 1) /= ndim .or. size(hvp, 2) /= self%nsph .or. size(hvp, 3) /= ndir) then
         call fatal_error(error, "get_surface_hessian_drop: hvp must be (3, nsph, ndir)")
         return
      end if
      if (self%ngrid <= 0) return

      !* --------------------------------- Both halves --------------------------------- *!
      allocate (hess_fixed(ndim, self%nsph, ndim, self%nsph), source=0.0_wp)
      allocate (total(ndim, self%nsph, ndir), source=0.0_wp)

      call surface_hessian_halves(self, acc, dirs, hess_fixed, total, error)
      if (allocated(error)) return

      !* -------------------------- Contract the fixed half ----------------------------- *!
      ! `total` already holds the response half of every direction, so the
      ! contraction lands on top of it and the two are summed exactly once.
      do idir = 1, ndir
         do iatom = 1, self%nsph
            do iaxis = 1, ndim
               total(:, :, idir) = total(:, :, idir) &
                                   + hess_fixed(:, :, iaxis, iatom)*dirs(iaxis, iatom, idir)
            end do
         end do
      end do

      hvp = hvp + total
   end subroutine get_surface_hessian_drop

   !> Dense nuclear Hessian of the DROP surface contribution
   !>
   !> Accumulates the full `(3, nsph, 3, nsph)` block for the energy whose
   !> surface adjoints `acc` holds. The result is *added* to `hessian`, and the
   !> accumulator is left untouched when anything fails.
   !>
   !> Column `(beta, B)` of the block is the Hessian-vector product along the
   !> Cartesian unit direction `e_(beta, B)`, and that is how it is obtained:
   !> all `3 nsph` unit directions are handed to the response half in a single
   !> call. The fixed half is direction free, so its rank-4 block is added
   !> column for column with no contraction at all -- which is both cheaper and
   !> exact, and leaves this path bit-for-bit equal to
   !> [[get_surface_hessian_drop]] driven with the same directions.
   !>
   !> @param[in]    self    DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc     Accumulated surface-observable adjoints
   !> @param[inout] hessian Nuclear-Hessian accumulator `(3, nsph, 3, nsph)`
   !> @param[out]   error   Error object, allocated on failure
   module subroutine get_hessian_drop(self, acc, hessian, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Cartesian unit directions, one per nuclear degree of freedom
      real(wp), allocatable :: dirs(:, :, :)
      !> Direction-free fixed half, and the response half of every column
      real(wp), allocatable :: hess_fixed(:, :, :, :), resp(:, :, :)
      !> Extents and loop indices
      integer :: ndir, idir, iatom, iaxis

      !* ------------------------------- Shape guards --------------------------------- *!
      call check_surface_adjoint(self, acc, "get_hessian_drop", error)
      if (allocated(error)) return
      if (any(shape(hessian) /= [ndim, self%nsph, ndim, self%nsph])) then
         call fatal_error(error, "get_hessian_drop: hessian shape mismatch")
         return
      end if
      if (self%ngrid <= 0 .or. self%nsph <= 0) return

      !* ---------------------------- Cartesian unit directions ------------------------ *!
      ndir = ndim*self%nsph
      allocate (dirs(ndim, self%nsph, ndir), source=0.0_wp)
      do iatom = 1, self%nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do

      !* --------------------------------- Both halves --------------------------------- *!
      allocate (hess_fixed(ndim, self%nsph, ndim, self%nsph), source=0.0_wp)
      allocate (resp(ndim, self%nsph, ndir), source=0.0_wp)

      call surface_hessian_halves(self, acc, dirs, hess_fixed, resp, error)
      if (allocated(error)) return

      do iatom = 1, self%nsph
         do iaxis = 1, ndim
            idir = ndim*(iatom - 1) + iaxis
            hessian(:, :, iaxis, iatom) = hessian(:, :, iaxis, iatom) &
                                          + hess_fixed(:, :, iaxis, iatom) + resp(:, :, idir)
         end do
      end do
   end subroutine get_hessian_drop

   !* ================================================================================= *!
   !*                              Composition of the halves                            *!
   !* ================================================================================= *!

   !> Evaluate both halves of the surface Hessian into caller-owned buffers
   !>
   !> The single place the two halves meet. `hess_fixed` and `resp` are written
   !> by the halves themselves, which *add* to what they are given, so both are
   !> expected zeroed on entry and are the caller's staging buffers rather than
   !> its accumulators -- that is what keeps a public accumulator untouched
   !> when the second half fails.
   !>
   !> @param[in]    self       DROP cavity instance
   !> @param[in]    acc        Accumulated surface-observable adjoints
   !> @param[in]    dirs       Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hess_fixed Fixed-adjoint half `(3, nsph, 3, nsph)`
   !> @param[inout] resp       Adjoint-response half `(3, nsph, ndir)`
   !> @param[out]   error      Error object, allocated on failure
   subroutine surface_hessian_halves(self, acc, dirs, hess_fixed, resp, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Fixed-adjoint half
      real(wp), intent(inout) :: hess_fixed(:, :, :, :)
      !> Adjoint-response half
      real(wp), intent(inout) :: resp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surrogate accumulator with the base-geometry fold already applied
      type(cavity_surface_adjoint_type) :: acc_frozen

      call freeze_surface_adjoint(self, acc, acc_frozen, error)
      if (allocated(error)) return

      ! The fixed half, with the surrogate: see the module header for why the
      ! substitution is an identity. A multi-branch grid is refused here, and
      ! that refusal is passed on rather than absorbed.
      call self%get_surface_hessian_fixed(acc_frozen, hess_fixed, error)
      if (allocated(error)) return

      ! The response half, with the raw accumulator: the fold is what moves.
      call self%get_surface_hessian_response(acc, dirs, resp, error)
   end subroutine surface_hessian_halves

   !> Build the frozen surrogate the fixed half can accept
   !>
   !> Replaces the width and switching channels by their base-geometry *folded*
   !> values and drops the area and integration-weight channels, so that
   !> [[prepare_surface_weights]] reproduces those same folded weights at every
   !> geometry. The remaining four channels are `source=`-copies in
   !> [[prepare_surface_weights]] anyway and are carried over unchanged.
   !>
   !> @param[in]  self   DROP cavity instance
   !> @param[in]  acc    Accumulated surface-observable adjoints
   !> @param[out] frozen Surrogate accumulator
   !> @param[out] error  Error object, allocated on failure
   subroutine freeze_surface_adjoint(self, acc, frozen, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Surrogate accumulator
      type(cavity_surface_adjoint_type), intent(out) :: frozen
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Base-geometry folded weights
      type(drop_surface_weights_type) :: eff

      ! `fold_switching = .true.`, as on the nuclear path: a nuclear
      ! displacement moves `f`, so the area channel's `da/df` term is part of
      ! the effective switching adjoint the fixed half has to see.
      call prepare_surface_weights(self, acc, .true., eff)

      ! `init` zeroes every channel, which is where `w_a` and `w_w` are left.
      call frozen%init(self%ngrid)
      call frozen%add_surface_weights(error, w_xi=eff%w_xi, w_f=eff%w_f, &
                                      w_xyz=acc%w_xyz, w_n=acc%w_n, &
                                      w_k1=acc%w_k1, w_k2=acc%w_k2)
   end subroutine freeze_surface_adjoint

end submodule moist_cavity_drop_derivatives_hessian
