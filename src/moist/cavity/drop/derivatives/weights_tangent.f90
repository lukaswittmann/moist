!> Tangent of the DROP surface-weight folding
!>
!> Pass 2 of the DROP Hessian. [[prepare_surface_weights]]
!> (`derivatives/weights.f90`) folds the raw host adjoints of a
!> `cavity_surface_adjoint_type` into the effective weights of a
!> [[drop_surface_weights_type]]; this module differentiates that folding along
!> one direction, **with the raw adjoints held fixed**.
!>
!> That is the fixed-adjoint half of the composition: the folded weights are
!> geometry dependent even when the host adjoints are not, so
!> `d(sum_c eff_c dGamma_c/dx)` carries a `sum_c d(eff_c) dGamma_c/dx` term and
!> `d(eff)` is what produces it.
!>
!> ## What moves and what does not
!>
!> Only three of the eight fields of [[drop_surface_weights_type]] have a
!> nonzero tangent, and the other five are frozen for a reason worth stating
!> once rather than rediscovering:
!>
!>   * `w_xyz`, `w_n`, `w_k1` and `w_k2` are `source=`-copies of the raw
!>     adjoints (`weights.f90`, the `allocate` block). With the raw adjoints
!>     fixed their tangent is identically zero, so this module does not emit
!>     them at all -- a consumer must not contract them, and gains the right to
!>     skip the normal and curvature channels of the tangent contraction
!>     outright.
!>   * `have_wn` and `have_wk` are `any(abs(...) > seed_weight_tol)` over those
!>     same copies. They are functions of the fixed raw adjoints only, so they
!>     do not move with the geometry either; the caller reuses the primal
!>     `eff%have_wn` / `eff%have_wk` unchanged.
!>
!> What is left is `w_xi` (both derived-channel folds), `w_f` (the area fold,
!> only when `fold_switching`), and `branch_phi_adj`.
!>
!> ## What the producer must supply
!>
!> The folding reads `a`, `wleb`, `xi0` and `wbranch` from the cavity, so its
!> tangent needs those four per grid point per direction, and nothing else that
!> moves: `radii` is geometry independent for both shipped radius models
!> (`radii/type.f90`, `f1_rA` zeroed), and `owner`, `branch_count`, `anchor_id`
!> and `sigma_phi` are discrete or parametric.
!>
!> The four are *not* independent -- `a = R^2 f wleb` and
!> `xi0 = swx/(R sqrt(wleb))`, so `dwleb = -2 wleb dxi0/xi0` exactly (the same
!> identity [[apply_seed]] uses at `res%dxi`) -- but they are taken as four
!> independent arguments here on purpose, for the reason
!> [[drop_seed_input_tangent_type]] gives: baking the identities into the
!> primitive would make it untestable along an arbitrary linear path and would
!> put the physics in the driver's place. The identities are the *caller's*
!> contract, restated on the routine.
!>
!> ## Parallelisation -- read this before writing the driver
!>
!> [[branch_phi_adj_tangent]] reduces one scalar over each contiguous anchor
!> group, exactly as [[compute_branch_phi_adj]] does. It is the only
!> cross-point coupling in the whole scheme. **A driver may parallelise it over
!> anchor groups or over directions, never over grid points**: a group split
!> across two threads silently corrupts `mean_adj_branch` and its tangent, with
!> no error and no obviously wrong answer. Directions are the easy axis --
!> every direction is independent of every other. Contiguity of a group is
!> guaranteed by the stable `counting_argsort` at `projection.f90:542-583`, not
!> by anything this module checks.
!>
!> The primitive itself is serial over the whole grid, so a driver that calls
!> it once per direction is already safe.
!>
!> TODO: with geometry-dependent radii (`radii.md`) the `w_f` fold gains a
!>       `2 w_a R dR wleb` term, and this primitive then needs a `dradii`
!>       argument for it. Nothing else here changes: the `w_xi` fold reaches
!>       the radii only through `a`, `wleb` and `xi0`, whose `dR` content
!>       arrives with the tangents the caller already supplies.
module moist_cavity_drop_derivatives_weights_tangent
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type, seed_weight_tol

   implicit none(type, external)
   private

   public :: prepare_surface_weights_tangent, branch_phi_adj_tangent

contains

   !> Tangent of [[prepare_surface_weights]] along one direction
   !>
   !> Mirrors the primal statement for statement: the two derived-channel folds
   !> into `w_xi`, the conditional area fold into `w_f`, then the branch reverse
   !> pass, which reads the *folded* width channel and therefore has to run
   !> last, on `dw_xi` rather than on zero.
   !>
   !> ## Guards
   !>
   !> The primal folds behind `abs(acc%w_a(i)) > seed_weight_tol` and
   !> `abs(acc%w_w(i)) > seed_weight_tol`. Those conditions read the raw
   !> adjoints, which this pass holds fixed, so they are *exactly* constant
   !> along the direction: taking the primal's branch is not the usual
   !> piecewise-tangent convention here, it is an identity, and a point that
   !> skips a fold skips its tangent for every displacement. That is what makes
   !> the else-branch a hard zero rather than an approximation. It would not
   !> hold if the adjoints moved -- with a moving `w_a` the guard would be a
   !> genuine kink and the tangent only correct strictly inside a piece.
   !>
   !> The guards inside [[branch_phi_adj_tangent]] are the opposite case and are
   !> documented there.
   !>
   !> ## Argument contract
   !>
   !> `acc` must be the accumulator [[prepare_surface_weights]] was called with,
   !> `eff` the weights it returned for that `acc`, and `a`, `wleb`, `xi0`,
   !> `wbranch`, `radii`, `owner`, `branch_count`, `anchor_id` and `sigma_phi`
   !> the cavity state it read. Only `eff%w_xi` is read, but the whole object is
   !> taken so that the pairing is visible at the call site.
   !>
   !> For a physical direction the tangents are related by the identities in
   !> the module header -- `dwleb = -2 wleb dxi0/xi0` and
   !> `da = R^2 (wleb df + f dwleb)`. This routine does not impose them; a
   !> caller that supplies an inconsistent triple gets the exact tangent of the
   !> folding along that inconsistent path, which is what makes an arbitrary
   !> finite-difference path a valid test and an inconsistent driver a silent
   !> physics error.
   !>
   !> @param[in]  acc             Raw surface adjoints, held fixed
   !> @param[in]  eff             Folded weights, as [[prepare_surface_weights]] returned them
   !> @param[in]  fold_switching  Whether the primal folded the area channel into `w_f`
   !> @param[in]  a               Surface area per grid point (ngrid)
   !> @param[in]  wleb            Final Lebedev weight per grid point (ngrid)
   !> @param[in]  xi0             Gaussian width per grid point (ngrid)
   !> @param[in]  wbranch         Softmax branch weight per grid point (ngrid)
   !> @param[in]  radii           Sphere radii (nsph)
   !> @param[in]  owner           Owner sphere per grid point (ngrid)
   !> @param[in]  branch_count    Branches in the point's anchor group (ngrid)
   !> @param[in]  anchor_id       Anchor group id per grid point (ngrid)
   !> @param[in]  sigma_phi       Softmax temperature, `branch_weight%s`
   !> @param[in]  da              Directional tangent of `a` (ngrid)
   !> @param[in]  dwleb           Directional tangent of `wleb` (ngrid)
   !> @param[in]  dxi0            Directional tangent of `xi0` (ngrid)
   !> @param[in]  dwbranch        Directional tangent of `wbranch` (ngrid)
   !> @param[out] dw_xi           Tangent of the folded width channel (ngrid)
   !> @param[out] dw_f            Tangent of the folded switching channel (ngrid)
   !> @param[out] dbranch_phi_adj Tangent of the branch-objective adjoint (ngrid)
   !> @param[out] error           Error object, allocated on inconsistent shapes
   subroutine prepare_surface_weights_tangent(acc, eff, fold_switching, &
                                              a, wleb, xi0, wbranch, radii, owner, &
                                              branch_count, anchor_id, sigma_phi, &
                                              da, dwleb, dxi0, dwbranch, &
                                              dw_xi, dw_f, dbranch_phi_adj, error)
      !> Raw surface adjoints, held fixed
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded weights, as returned by [[prepare_surface_weights]] for `acc`
      type(drop_surface_weights_type), intent(in) :: eff
      !> Whether the area channel also folded into `w_f`
      logical, intent(in) :: fold_switching
      !> Cavity grid scalars the folding reads
      real(wp), contiguous, intent(in) :: a(:), wleb(:), xi0(:), wbranch(:), radii(:)
      !> Owner sphere and branch bookkeeping per grid point
      integer, contiguous, intent(in) :: owner(:), branch_count(:), anchor_id(:)
      !> Softmax temperature of the branch weights
      real(wp), intent(in) :: sigma_phi
      !> Directional tangents of the cavity grid scalars
      real(wp), contiguous, intent(in) :: da(:), dwleb(:), dxi0(:), dwbranch(:)
      !> Tangent of the folded weights
      real(wp), contiguous, intent(out) :: dw_xi(:), dw_f(:), dbranch_phi_adj(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grid extent and grid index
      integer :: ngrid, igrid
      !> Owner radius of the current grid point
      real(wp) :: r_own
      !> Rendered shapes; fixed length, so the message build is thread safe
      character(len=64) :: shapes

      if (.not. acc%is_initialized() .or. .not. allocated(eff%w_xi)) then
         call fatal_error(error, "Surface-weight tangent: the accumulator is not"// &
                          " initialized or the folded weights are unallocated")
         return
      end if
      ngrid = size(acc%w_xi)

      if (size(eff%w_xi) /= ngrid .or. &
          size(a) /= ngrid .or. size(wleb) /= ngrid .or. size(xi0) /= ngrid .or. &
          size(wbranch) /= ngrid .or. size(owner) /= ngrid .or. &
          size(branch_count) /= ngrid .or. size(anchor_id) /= ngrid .or. &
          size(da) /= ngrid .or. size(dwleb) /= ngrid .or. size(dxi0) /= ngrid .or. &
          size(dwbranch) /= ngrid .or. size(dw_xi) /= ngrid .or. size(dw_f) /= ngrid .or. &
          size(dbranch_phi_adj) /= ngrid) then
         write (shapes, "(2(a, i0))") "ngrid ", ngrid, ", smallest argument ", &
            min(size(eff%w_xi), size(a), size(wleb), size(xi0), size(wbranch), &
                size(owner), size(branch_count), size(anchor_id), size(da), &
                size(dwleb), size(dxi0), size(dwbranch), size(dw_xi), size(dw_f), &
                size(dbranch_phi_adj))
         call fatal_error(error, "Surface-weight tangent has inconsistent shapes ("// &
                          trim(shapes)//"; every grid argument must be ngrid)")
         return
      end if

      dw_xi = 0.0_wp
      dw_f = 0.0_wp

      ! d(-2 a w_a/xi) and d(-2 wleb w_w/xi) with w_a, w_w fixed. Written as
      ! (dX - X dxi/xi)/xi rather than dX/xi - X dxi/xi^2 so the second term
      ! never forms 1/xi^2 explicitly; check_surface_adjoint has already
      ! rejected a vanishing xi under a live area or weight adjoint.
      do igrid = 1, ngrid
         if (abs(acc%w_a(igrid)) > seed_weight_tol) then
            dw_xi(igrid) = dw_xi(igrid) - 2.0_wp*acc%w_a(igrid) &
                           *(da(igrid) - a(igrid)*dxi0(igrid)/xi0(igrid))/xi0(igrid)
         end if
         if (abs(acc%w_w(igrid)) > seed_weight_tol) then
            dw_xi(igrid) = dw_xi(igrid) - 2.0_wp*acc%w_w(igrid) &
                           *(dwleb(igrid) - wleb(igrid)*dxi0(igrid)/xi0(igrid))/xi0(igrid)
         end if
      end do

      ! d(w_a R^2 wleb). The radius is the only factor of the primal's area
      ! fold that this pass treats as frozen; see the module TODO.
      if (fold_switching) then
         do igrid = 1, ngrid
            if (abs(acc%w_a(igrid)) > seed_weight_tol) then
               r_own = radii(owner(igrid))
               dw_f(igrid) = dw_f(igrid) &
                             + acc%w_a(igrid)*r_own*r_own*dwleb(igrid)
            end if
         end do
      end if

      ! Runs last and on the folded `dw_xi`, because the primal branch pass runs
      ! last and on the folded `w_xi`.
      call branch_phi_adj_tangent(branch_count, anchor_id, wbranch, wleb, xi0, &
                                  sigma_phi, eff%w_xi, dwbranch, dwleb, dxi0, dw_xi, &
                                  dbranch_phi_adj)

   end subroutine prepare_surface_weights_tangent

   !> Tangent of the branch-weight reverse pass, [[compute_branch_phi_adj]]
   !>
   !> Same group walk, same early exits, same guard conditions; the only
   !> difference is that every quantity carries its directional derivative
   !> alongside it. Per contiguous anchor group the primal forms
   !>
   !> ```
   !>   adj_m  = adj_wleb_m * factor_m,  mean = sum_m wbranch_m adj_m
   !>   Phi_m  = -wbranch_m (adj_m - mean)/sigma_phi
   !> ```
   !>
   !> so the tangent needs `dmean` before it can write any point of the group.
   !> That is the coupling: a point's `dbranch_phi_adj` depends on every other
   !> point of its group, through `dmean`, even when its own `adj_m` is gated
   !> off. Hence the two passes, and hence the parallelisation constraint in the
   !> module header.
   !>
   !> ## Guards
   !>
   !> Unlike the folds in [[prepare_surface_weights_tangent]], the per-point
   !> gate here reads `w_xi`, `wleb` and `wbranch`, all of which *do* move with
   !> the geometry. So this is the ordinary piecewise case and it follows the
   !> project convention: take the primal's branch on the primal's condition,
   !> zero on the else, never a separate second-order threshold. The tangent is
   !> then exactly the tangent of the value it differentiates everywhere except
   !> on the measure-zero switching surface itself.
   !>
   !> The early exits (`sigma_phi <= seed_weight_tol`, no group with
   !> `branch_count > 1`) are frozen discrete choices in the same sense.
   !>
   !> @param[in]  branch_count    Branches per grid point (ngrid)
   !> @param[in]  anchor_id       Anchor group id per grid point (ngrid)
   !> @param[in]  wbranch         Softmax branch weight per grid point (ngrid)
   !> @param[in]  wleb            Lebedev weight per grid point (ngrid)
   !> @param[in]  xi0             Gaussian width per grid point (ngrid)
   !> @param[in]  sigma_phi       Softmax temperature
   !> @param[in]  w_xi            Folded width adjoint, `eff%w_xi` (ngrid)
   !> @param[in]  dwbranch        Directional tangent of `wbranch` (ngrid)
   !> @param[in]  dwleb           Directional tangent of `wleb` (ngrid)
   !> @param[in]  dxi0            Directional tangent of `xi0` (ngrid)
   !> @param[in]  dw_xi           Directional tangent of the folded `w_xi` (ngrid)
   !> @param[out] dbranch_phi_adj Tangent of the branch-objective adjoint (ngrid)
   pure subroutine branch_phi_adj_tangent(branch_count, anchor_id, wbranch, wleb, xi0, &
                                          sigma_phi, w_xi, dwbranch, dwleb, dxi0, dw_xi, &
                                          dbranch_phi_adj)
      !> Branch bookkeeping per grid point
      integer, intent(in) :: branch_count(:), anchor_id(:)
      !> Branch weight, Lebedev weight and Gaussian width per grid point
      real(wp), intent(in) :: wbranch(:), wleb(:), xi0(:)
      !> Softmax temperature
      real(wp), intent(in) :: sigma_phi
      !> Effective Gaussian-width adjoint
      real(wp), intent(in) :: w_xi(:)
      !> Directional tangents of the same quantities
      real(wp), intent(in) :: dwbranch(:), dwleb(:), dxi0(:), dw_xi(:)
      !> Tangent of the adjoint of the branch objective
      real(wp), intent(out) :: dbranch_phi_adj(:)

      !> Grid extent and group bookkeeping
      integer :: ngrid, igroup_start, igroup_end, group_size, m_branch, im_grid
      !> Per-point branch adjoint and its tangent
      real(wp) :: adj_branch, dadj_branch
      !> Group reduction and its tangent
      real(wp) :: mean_adj_branch, dmean_adj_branch

      dbranch_phi_adj = 0.0_wp
      ngrid = size(branch_count)
      if (ngrid <= 0) return
      if (.not. any(branch_count > 1)) return
      if (sigma_phi <= seed_weight_tol) return

      igroup_start = 1
      do while (igroup_start <= ngrid)
         if (branch_count(igroup_start) <= 1) then
            igroup_start = igroup_start + 1
            cycle
         end if

         ! Extend the group while anchor_id stays the same
         igroup_end = igroup_start
         do while (igroup_end < ngrid)
            if (anchor_id(igroup_end + 1) /= anchor_id(igroup_start)) exit
            igroup_end = igroup_end + 1
         end do
         group_size = igroup_end - igroup_start + 1

         mean_adj_branch = 0.0_wp
         dmean_adj_branch = 0.0_wp
         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            call branch_point_adjoint(w_xi(im_grid), wleb(im_grid), xi0(im_grid), &
                                      wbranch(im_grid), dw_xi(im_grid), dwleb(im_grid), &
                                      dxi0(im_grid), dwbranch(im_grid), &
                                      adj_branch, dadj_branch)
            mean_adj_branch = mean_adj_branch + wbranch(im_grid)*adj_branch
            dmean_adj_branch = dmean_adj_branch + dwbranch(im_grid)*adj_branch &
                               + wbranch(im_grid)*dadj_branch
         end do

         ! The per-point pair is recomputed rather than parked in a scratch
         ! array: the group is 2-4 points wide, so a buffer would cost an
         ! allocation or an automatic array bounded only by ngrid, to save four
         ! flops and two divisions.
         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            call branch_point_adjoint(w_xi(im_grid), wleb(im_grid), xi0(im_grid), &
                                      wbranch(im_grid), dw_xi(im_grid), dwleb(im_grid), &
                                      dxi0(im_grid), dwbranch(im_grid), &
                                      adj_branch, dadj_branch)
            dbranch_phi_adj(im_grid) = &
               -(dwbranch(im_grid)*(adj_branch - mean_adj_branch) &
                 + wbranch(im_grid)*(dadj_branch - dmean_adj_branch))/sigma_phi
         end do

         igroup_start = igroup_end + 1
      end do

   end subroutine branch_phi_adj_tangent

   !> One point's branch adjoint and its directional tangent
   !>
   !> The primal half is [[compute_branch_phi_adj]]'s inner block verbatim, so
   !> the two stay comparable line by line; the tangent is its product rule.
   !> `adj_wleb * factor` is kept factored rather than collapsed to
   !> `-0.5 w_xi xi0/wbranch` -- the `wleb` cancels algebraically -- because the
   !> factored form is the one the primal evaluates and the one a reviewer can
   !> diff against it.
   !>
   !> @param[in]  w_xi        Folded width adjoint at the point
   !> @param[in]  wleb        Lebedev weight at the point
   !> @param[in]  xi0         Gaussian width at the point
   !> @param[in]  wbranch     Softmax branch weight at the point
   !> @param[in]  dw_xi       Tangent of the folded width adjoint
   !> @param[in]  dwleb       Tangent of the Lebedev weight
   !> @param[in]  dxi0        Tangent of the Gaussian width
   !> @param[in]  dwbranch    Tangent of the branch weight
   !> @param[out] adj_branch  Branch adjoint, zero when the primal gate is shut
   !> @param[out] dadj_branch Its directional tangent, zero on the same gate
   pure subroutine branch_point_adjoint(w_xi, wleb, xi0, wbranch, &
                                        dw_xi, dwleb, dxi0, dwbranch, &
                                        adj_branch, dadj_branch)
      !> Point values
      real(wp), intent(in) :: w_xi, wleb, xi0, wbranch
      !> Point tangents
      real(wp), intent(in) :: dw_xi, dwleb, dxi0, dwbranch
      !> Branch adjoint and its tangent
      real(wp), intent(out) :: adj_branch, dadj_branch

      !> Weight-adjoint scratch and its tangent
      real(wp) :: adj_wleb, dadj_wleb, factor_m, dfactor_m

      adj_branch = 0.0_wp
      dadj_branch = 0.0_wp
      if (abs(w_xi) > seed_weight_tol &
          .and. wleb > seed_weight_tol &
          .and. wbranch > tiny(1.0_wp)) then
         adj_wleb = -0.5_wp*w_xi*xi0/wleb
         dadj_wleb = -0.5_wp*(dw_xi*xi0 + w_xi*dxi0)/wleb &
                     + 0.5_wp*w_xi*xi0*dwleb/(wleb*wleb)
         factor_m = wleb/wbranch
         dfactor_m = dwleb/wbranch - wleb*dwbranch/(wbranch*wbranch)
         adj_branch = adj_wleb*factor_m
         dadj_branch = dadj_wleb*factor_m + adj_wleb*dfactor_m
      end if

   end subroutine branch_point_adjoint

end module moist_cavity_drop_derivatives_weights_tangent
