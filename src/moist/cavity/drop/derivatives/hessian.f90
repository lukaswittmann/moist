!> Public DROP Hessian accessors
!>
!> The nuclear gradient of `nuclear.f90` is `J^T omega`, and its directional
!> derivative splits into
!>
!>     d/dv [ J^T omega ]  =  (dJ^T/dv) omega  +  J^T (d omega/dv)
!>                            ^ fixed channel     ^ response channel
!>
!> Both halves are channels of the one traversal next door,
!> [[drop_hessian_traverse]]. This submodule owns nothing of the mathematics;
!> it owns the *composition* -- which form of the fixed channel to ask for --
!> and the two public entry points:
!>
!>   * [[get_surface_hessian_drop]], the Hessian-vector product along a set of
!>     directions;
!>   * [[get_hessian_drop]], the dense `(3, nsph, 3, nsph)` block.
!>
!>
!> ## Choosing the form of the fixed channel
!>
!> The response half is intrinsically per direction. The fixed half is offered
!> in two forms (see the traversal's header): the direction-free rank-4 block,
!> whose explicit nuclear motion is one level-set block per grid point and
!> whose second-order chain runs on a fixed 23-element basis per point, and
!> the per-direction column, built from the supplied directions alone. Per grid
!> point their costs are
!>
!>     rank-4:         the fills + one `vjp_f2_rArB` block + 23 chains
!>                     + two matrix products
!>     per direction:  the fills + the explicit motion (one `hvp_jet_rA` pass
!>                     per direction below `drop_hvp_apply_min` directions, one
!>                     applied block `vjp_f2_rArB_apply` from there on)
!>                     + ndir chains and contractions
!>
!> and the `O(n_active)` kernel passes dominate both, so the crossover is a
!> *number of directions*, not a fraction of the basis. Measured on the
!> polyalanine set (SvdW, 110-point Lebedev grids, one thread): the per-direction
!> form costs about 0.13 s + 0.022 s per direction (83 atoms) and 0.75 s +
!> 0.062 s per direction (163 atoms) -- the intercept is the point's fills and
!> the applied block, the slope one chain and its contractions -- against a
!> rank-4 fixed half of 0.52 s and 1.85 s, so the two meet at about 12 and 18
!> directions. [[hvp_fixed_mode]] therefore runs per direction up to
!> `drop_hvp_per_dir_max` directions and rank-4 beyond that: the rank-4 form is
!> never worse than the dense path, and the per-direction form is never asked
!> to do more than about one dense fixed half's worth of work. The bound is a
!> constant because the quadratic terms of the rank-4 form (the two matrix
!> products and the scatter) are still small at these sizes; on much larger
!> active sets the true crossover moves up and the rule errs towards rank-4,
!> which is bounded, rather than towards the unbounded per-direction cost.
!> Beyond the bound a Hessian-vector product also carries the rank-4 form's
!> memory: the dense `(3, nsph, 3, nsph)` staging block below and the
!> traversal's per-thread sparse accumulators.
!>
!> [[get_hessian_drop]] asks for the rank-4 form with all `3 nsph` unit
!> directions handed to the response half **in one batch**: the traversal
!> walks the grid once (per direction block) either way, where `3 nsph`
!> separate calls would walk it `3 nsph` times. The rank-4 block is added
!> column for column with no contraction at all, so this path is bit-for-bit
!> what [[get_surface_hessian_drop]] returns for the same unit directions.
!>
!>
!> ## Where the fold happens
!>
!> What the host accumulates is not what either half contracts.
!> [[prepare_surface_weights]] folds the area and integration-weight channels
!> into `w_xi` and `w_f` through `a`, `wleb`, `xi0` and the radii, and derives
!> `branch_phi_adj` from the branch softmax; everything it returns -- `eff(R)`
!> -- is a function of the geometry. That fold is performed **once**, in
!> [[surface_hessian_halves]], and the result is handed to both halves: the
!> fixed half takes `eff` alone (its term is `(dPhi/dv) . eff` with the weights
!> held fixed), the response half takes `eff` **and** the raw `acc` (its term
!> is `Phi . (d eff/dv)`, and [[prepare_surface_weights_tangent]]
!> differentiates the fold out of the raw channels and the primal `eff`
!> together). Writing the gradient as `G(R) = Phi(R) . eff(R)`, the two halves
!> are the two terms of the product rule, and the split is exact because `G`
!> is linear in `eff`. The seam sits here so that `eff` is the *same object* in
!> both halves by construction.
!>
!>
!> ## What the composite refuses
!>
!> A grid carrying a multi-branch anchor group. There the projected point is a
!> softmax over several anchors, and its second derivative carries a branch
!> term that **neither** half supplies: the fixed half's second-order chain
!> omits it (see [[seed_contribution_tangent]]), and the response half moves
!> `branch_phi_adj` without differentiating the branch geometry a second time.
!> Both public entry points reject such a grid rather than return a Hessian
!> silently short a term. The refusal lives here and not in the halves: the
!> response half is correct on a branched grid on its own.
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
   !> anything fails -- both halves are formed in local buffers first.
   !>
   !> With `omega_v` present the model's adjoint response is folded into the
   !> traversal's response channel, block by block, so the columns are those
   !> of the full model Hessian; without it they are the frozen-adjoint ones.
   !>
   !> @param[in]    self    DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc     Accumulated surface-observable adjoints
   !> @param[in]    dirs    Nuclear directions `(3, nsph, ndir)`
   !> @param[inout] hvp     Hessian-vector accumulator `(3, nsph, ndir)`
   !> @param[out]   error   Error object, allocated on failure
   !> @param[inout] omega_v Surface-adjoint response of the model, optional
   module subroutine get_surface_hessian_drop(self, acc, dirs, hvp, error, omega_v)
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
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v

      !> Direction-free fixed half (rank-4 form only), and the staged columns
      real(wp), allocatable :: hess_fixed(:, :, :, :), total(:, :, :)
      !> Form of the fixed channel
      integer :: fixed_mode
      !> Extents and loop indices
      integer :: ndir, idir, iatom, iaxis

      !* ------------------------------- Shape guards --------------------------------- *!
      call check_surface_adjoint(self, acc, "get_surface_hessian_drop", error)
      if (allocated(error)) return
      call check_single_branch(self, "get_surface_hessian_drop", error)
      if (allocated(error)) return
      call check_direction_set(self, dirs, hvp, "get_surface_hessian_drop", error)
      if (allocated(error)) return
      if (self%ngrid <= 0) return
      ndir = size(dirs, 3)

      !* --------------------------------- Both halves -------------------------------- *!
      fixed_mode = hvp_fixed_mode(ndir, self%nsph)
      allocate (total(ndim, self%nsph, ndir), source=0.0_wp)

      if (fixed_mode == drop_fixed_per_dir) then
         ! Both channels land their columns in `total` directly
         call surface_hessian_halves(self, acc, dirs, fixed_mode, "get_surface_hessian_drop", &
                                     total, error, omega_v=omega_v)
         if (allocated(error)) return
      else
         allocate (hess_fixed(ndim, self%nsph, ndim, self%nsph), source=0.0_wp)
         call surface_hessian_halves(self, acc, dirs, fixed_mode, "get_surface_hessian_drop", &
                                     total, error, hess_fixed=hess_fixed, omega_v=omega_v)
         if (allocated(error)) return

         !* -------------------------- Contract the fixed half ------------------------ *!
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
      end if

      hvp = hvp + total
   end subroutine get_surface_hessian_drop

   !> Dense nuclear Hessian of the DROP surface contribution
   !>
   !> Accumulates the full `(3, nsph, 3, nsph)` block for the energy whose
   !> surface adjoints `acc` holds. The result is *added* to `hessian`, and the
   !> accumulator is left untouched when anything fails.
   !>
   !> Column `(beta, B)` of the block is the Hessian-vector product along the
   !> Cartesian unit direction `e_(beta, B)`: all `3 nsph` unit directions are
   !> handed to the response half in a single call, and the rank-4 fixed half
   !> is added column for column with no contraction at all.
   !>
   !> @param[in]    self    DROP cavity instance (must hold a projected grid)
   !> @param[in]    acc     Accumulated surface-observable adjoints
   !> @param[inout] hessian Nuclear-Hessian accumulator `(3, nsph, 3, nsph)`
   !> @param[out]   error   Error object, allocated on failure
   !> @param[inout] omega_v Surface-adjoint response of the model, optional
   module subroutine get_hessian_drop(self, acc, hessian, error, omega_v)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v

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

      call surface_hessian_halves(self, acc, dirs, drop_fixed_rank4, "get_hessian_drop", &
                                  resp, error, hess_fixed=hess_fixed, omega_v=omega_v)
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

   !> Form of the fixed channel for a Hessian-vector product of `ndir` directions
   !>
   !> Per direction up to `drop_hvp_per_dir_max` directions and short of a full
   !> Cartesian basis, rank-4 otherwise; the module header has the cost model
   !> and the measurements the bound comes from.
   !>
   !> @param[in] ndir Directions asked for
   !> @param[in] nsph Spheres of the cavity
   !> @return         `drop_fixed_per_dir` or `drop_fixed_rank4`
   pure function hvp_fixed_mode(ndir, nsph) result(mode)
      !> Directions asked for
      integer, intent(in) :: ndir
      !> Spheres of the cavity
      integer, intent(in) :: nsph
      !> Form of the fixed channel
      integer :: mode

      if (ndir <= drop_hvp_per_dir_max .and. ndir < ndim*nsph) then
         mode = drop_fixed_per_dir
      else
         mode = drop_fixed_rank4
      end if
   end function hvp_fixed_mode

   !> Evaluate both halves of the surface Hessian into caller-owned buffers
   !>
   !> The single place the two halves meet, and the single place the surface
   !> adjoints are folded. Both buffers are *added* to by the traversal, so
   !> they are expected zeroed on entry and are the caller's staging buffers
   !> rather than its accumulators -- that is what keeps a public accumulator
   !> untouched when anything fails.
   !>
   !> `columns` receives the response half of every direction and, in the
   !> per-direction mode, the fixed half as well; `hess_fixed` receives the
   !> rank-4 fixed half and must be present in that mode.
   !>
   !> `context` is the public entry point the user actually called, threaded
   !> down so that a singular bordered system at some grid point names it.
   !>
   !> @param[in]    self       DROP cavity instance
   !> @param[in]    acc        Accumulated surface-observable adjoints
   !> @param[in]    dirs       Nuclear directions `(3, nsph, ndir)`
   !> @param[in]    fixed_mode Form of the fixed channel
   !> @param[in]    context    Calling routine, used to prefix the diagnostics
   !> @param[inout] columns    Per-direction half or halves `(3, nsph, ndir)`
   !> @param[out]   error      Error object, allocated on failure
   !> @param[inout] hess_fixed Rank-4 fixed half `(3, nsph, 3, nsph)`
   !> @param[inout] omega_v    Surface-adjoint response of the model, optional
   subroutine surface_hessian_halves(self, acc, dirs, fixed_mode, context, columns, error, &
                                     hess_fixed, omega_v)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Accumulated surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Form of the fixed channel
      integer, intent(in) :: fixed_mode
      !> Calling routine, so a failure names the entry point the user called
      character(len=*), intent(in) :: context
      !> Per-direction columns
      real(wp), intent(inout) :: columns(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Rank-4 fixed half
      real(wp), intent(inout), optional :: hess_fixed(:, :, :, :)
      !> Surface-adjoint response of the model
      class(surface_adjoint_response_type), intent(inout), optional :: omega_v

      !> Folded surface adjoints of the base geometry, read by both channels
      type(drop_surface_weights_type) :: eff

      ! `fold_switching = .true.`, as on the nuclear path: a nuclear
      ! displacement moves `f`, so the area channel's `da/df` term is part of
      ! the effective switching adjoint both halves have to see.
      call prepare_surface_weights(self, acc, .true., eff)

      ! The fixed channel is `(dPhi/dv) . eff` with the folded weights held
      ! fixed and never sees the raw ones. The response channel is
      ! `Phi . (d eff/dv)`: the fold is what moves there, so it needs the raw
      ! channels as well as the primal `eff`. Hence both objects go in, and one
      ! traversal serves the two.
      if (fixed_mode == drop_fixed_rank4) then
         call drop_hessian_traverse(self, eff, fixed_mode, .true., context, &
                                    acc=acc, dirs=dirs, hess_fixed=hess_fixed, hvp=columns, &
                                    error=error, omega_v=omega_v)
      else
         call drop_hessian_traverse(self, eff, fixed_mode, .true., context, &
                                    acc=acc, dirs=dirs, hvp=columns, error=error, &
                                    omega_v=omega_v)
      end if
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
