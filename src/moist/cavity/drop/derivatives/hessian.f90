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
!> ## Where the fold happens
!>
!> What the host accumulates is not what either half contracts.
!> [[prepare_surface_weights]] folds the area and integration-weight channels
!> into `w_xi` and `w_f` through `a`, `wleb`, `xi0` and the radii, and derives
!> `branch_phi_adj` from the branch softmax; everything it returns -- call it
!> `eff(R)` -- is a function of the geometry. That fold is performed **once**,
!> here in [[surface_hessian_halves]], and the result is handed to both halves:
!>
!>   * the fixed half takes `eff` and nothing else. Its term is
!>     `(dPhi/dv) . eff` with the weights held fixed, so the raw channels can
!>     tell it nothing it does not already read out of `eff`;
!>   * the response half takes `eff` **and** the raw `acc`. Its term is
!>     `Phi . (d eff/dv)`, and [[prepare_surface_weights_tangent]]
!>     differentiates the fold out of the raw channels and the primal `eff`
!>     together, so it genuinely needs both.
!>
!> Writing the gradient as `G(R) = Phi(R) . eff(R)`, the two halves are the two
!> terms of the product rule, and the split is exact because `G` is linear in
!> `eff` -- every channel enters exactly once.
!>
!> The seam sits at this level rather than inside the halves so that `eff` is
!> the *same object* in both of them by construction. Folding separately in
!> each half would make that an argument to be made rather than a fact, and it
!> is the argument that used to force a surrogate accumulator through the fixed
!> half's door.
!>
!>
!> ## What the composite still refuses
!>
!> A grid carrying a multi-branch anchor group. There the projected point is a
!> softmax over several anchors, and its second derivative carries a branch
!> term that **neither** half supplies: the fixed half's second-order chain
!> omits it (see [[seed_contribution_tangent]]), and the response half moves
!> `branch_phi_adj` without differentiating the branch geometry a second time.
!> Both public entry points below reject such a grid rather than return a
!> Hessian silently short a term.
!>
!> The refusal lives here and not in the halves on purpose. The response half
!> is correct on a branched grid when it is asked for on its own -- guarding
!> inside it would reject calls it handles -- and the fixed half is only ever
!> reachable through this composition.
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
      call check_single_branch(self, "get_surface_hessian_drop", error)
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

      !* --------------------------------- Both halves -------------------------------- *!
      allocate (hess_fixed(ndim, self%nsph, ndim, self%nsph), source=0.0_wp)
      allocate (total(ndim, self%nsph, ndir), source=0.0_wp)

      call surface_hessian_halves(self, acc, dirs, hess_fixed, total, error)
      if (allocated(error)) return

      !* --------------------------- Contract the fixed half -------------------------- *!
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
      call check_single_branch(self, "get_hessian_drop", error)
      if (allocated(error)) return
      if (any(shape(hessian) /= [ndim, self%nsph, ndim, self%nsph])) then
         call fatal_error(error, "get_hessian_drop: hessian shape mismatch")
         return
      end if
      if (self%ngrid <= 0 .or. self%nsph <= 0) return

      !* -------------------------- Cartesian unit directions ------------------------- *!
      ndir = ndim*self%nsph
      allocate (dirs(ndim, self%nsph, ndir), source=0.0_wp)
      do iatom = 1, self%nsph
         do iaxis = 1, ndim
            dirs(iaxis, iatom, ndim*(iatom - 1) + iaxis) = 1.0_wp
         end do
      end do

      !* --------------------------------- Both halves -------------------------------- *!
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
   !> The single place the two halves meet, and the single place the surface
   !> adjoints are folded. `hess_fixed` and `resp` are written by the halves
   !> themselves, which *add* to what they are given, so both are expected
   !> zeroed on entry and are the caller's staging buffers rather than its
   !> accumulators -- that is what keeps a public accumulator untouched when the
   !> second half fails.
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

      !> Folded surface adjoints of the base geometry, read by both halves
      type(drop_surface_weights_type) :: eff

      ! `fold_switching = .true.`, as on the nuclear path: a nuclear
      ! displacement moves `f`, so the area channel's `da/df` term is part of
      ! the effective switching adjoint both halves have to see.
      call prepare_surface_weights(self, acc, .true., eff)

      ! The fixed half: `(dPhi/dv) . eff`, the folded weights held fixed. It
      ! never sees the raw channels, because there is nothing in them it could
      ! use that `eff` does not already carry.
      call self%get_surface_hessian_fixed(eff, hess_fixed, error)
      if (allocated(error)) return

      ! The response half: `Phi . (d eff/dv)`. The fold is what moves here, so
      ! this one needs the raw channels it was built from as well as the
      ! primal `eff` it is differentiated around.
      call self%get_surface_hessian_response(acc, eff, dirs, resp, error)
   end subroutine surface_hessian_halves

   !* ================================================================================= *!
   !*                              Scope of the composite                               *!
   !* ================================================================================= *!

   !> Reject a grid whose anchor groups branch
   !>
   !> The one restriction the composite carries, and the module header has why:
   !> a branched anchor group puts a second-order branch term into `d/dv G` that
   !> neither half offers, so the sum of the two would be a Hessian missing a
   !> term with nothing to say so. `branch_count` is `1` everywhere on an
   !> unbranched grid and unallocated before the first projection, so both are
   !> accepted.
   !>
   !> Every message is prefixed with the caller's name, as
   !> [[check_surface_adjoint]] does, so a failure names the entry point the
   !> user actually called.
   !>
   !> @param[in]  self    DROP cavity instance
   !> @param[in]  context Calling routine, used to prefix the diagnostics
   !> @param[out] error   Error object, allocated on a branched grid
   subroutine check_single_branch(self, context, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%branch_count)) return
      if (self%ngrid <= 0) return
      if (any(self%branch_count(1:self%ngrid) > 1)) then
         call fatal_error(error, context//": multi-branch anchor groups carry a"// &
                          " second-order branch term that neither half of the surface"// &
                          " Hessian supplies")
      end if
   end subroutine check_single_branch

end submodule moist_cavity_drop_derivatives_hessian
