!> Tangent of the DROP surface-weight folding
!>
!> - Pass 2 of the DROP Hessian: directional tangent of
!>   [[prepare_surface_weights]] with the raw adjoints held fixed
!> - Supplies the `sum_c d(eff_c) dGamma_c/dx` term of the fixed-adjoint half
!> - Nonzero tangents: `w_xi`, `w_f` (only with `fold_switching`) and
!>   `branch_phi_adj`
!> - `w_xyz`, `w_n`, `w_k1`, `w_k2` copy the raw adjoints: zero tangent, not
!>   emitted, consumers must not contract them
!> - `have_wn`, `have_wk` do not move; reuse the primal `eff` flags
!> - Producer supplies the tangents of `a`, `wleb`, `xi0`, `wbranch` per grid
!>   point and direction
!> - `radii` fixed for both shipped radius models (`f1_rA` zero); `owner`,
!>   `branch_count`, `anchor_id`, `sigma_phi` discrete or parametric
!> - Tangents taken as independent arguments, so any linear path is testable
!> - Caller's contract: `dwleb = -2 wleb dxi0/xi0`,
!>   `da = R^2 (wleb df + f dwleb)`
!> - [[branch_phi_adj_tangent]] reduces over contiguous anchor groups, the only
!>   cross-point coupling: parallelise over groups or directions, never over
!>   grid points; a split group corrupts the result without an error
!> - Group contiguity from the stable `counting_argsort` in `projection.f90`,
!>   unchecked here; shared walk is [[next_branch_group]]
!> - Primitive serial over the grid, so one call per direction is safe
!>
!> TODO: radius model with nonzero `f1_rA` adds `2 w_a R dR wleb` to the `w_f`
!>       fold and needs a `dradii` argument; the `w_xi` fold stays as it is
module moist_cavity_drop_derivatives_weights_tangent
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_drop_derivatives_kernel, only: drop_surface_weights_type, &
      & seed_weight_tol, branch_point_adjoint, next_branch_group

   implicit none(type, external)
   private

   public :: prepare_surface_weights_tangent, branch_phi_adj_tangent

contains

   !> Tangent of [[prepare_surface_weights]] along one direction
   !>
   !> - Order of the primal: two folds into `w_xi`, area fold into `w_f`, branch
   !>   reverse pass last on the folded `dw_xi`
   !> - Fold guards compare the fixed raw `w_a`, `w_w` with `seed_weight_tol`:
   !>   constant along the direction, so a skipped fold is an exact zero
   !> - Exact zero holds only while the raw adjoints are fixed
   !> - Guards of [[branch_phi_adj_tangent]] read moving quantities, see there
   !> - `acc`, `eff` and the cavity state must be those of the primal call
   !> - Only `eff%w_xi` is read
   !> - Identities between `da`, `dwleb`, `dxi0` (module header) are not
   !>   imposed: inconsistent tangents give the exact tangent along that path,
   !>   a silent physics error in a driver
   !>
   !> @param[in]  acc              raw surface adjoints, held fixed
   !> @param[in]  eff              folded weights from [[prepare_surface_weights]]
   !> @param[in]  fold_switching   whether the primal folded the area channel into `w_f`
   !> @param[in]  a                surface area per grid point `(ngrid)`
   !> @param[in]  wleb             final Lebedev weight per grid point `(ngrid)`
   !> @param[in]  xi0              Gaussian width per grid point `(ngrid)`
   !> @param[in]  wbranch          softmax branch weight per grid point `(ngrid)`
   !> @param[in]  radii            sphere radii `(nsph)`
   !> @param[in]  owner            owner sphere per grid point `(ngrid)`
   !> @param[in]  branch_count     branches in the point's anchor group `(ngrid)`
   !> @param[in]  anchor_id        anchor group id per grid point `(ngrid)`
   !> @param[in]  sigma_phi        softmax temperature, `branch_weight%s`
   !> @param[in]  da               directional tangent of `a` `(ngrid)`
   !> @param[in]  dwleb            directional tangent of `wleb` `(ngrid)`
   !> @param[in]  dxi0             directional tangent of `xi0` `(ngrid)`
   !> @param[in]  dwbranch         directional tangent of `wbranch` `(ngrid)`
   !> @param[out] dw_xi            tangent of the folded width channel `(ngrid)`
   !> @param[out] dw_f             tangent of the folded switching channel `(ngrid)`
   !> @param[out] dbranch_phi_adj  tangent of the branch-objective adjoint `(ngrid)`
   !> @param[out] error            error object, allocated on inconsistent shapes
   subroutine prepare_surface_weights_tangent(acc, eff, fold_switching, &
                                              a, wleb, xi0, wbranch, radii, owner, &
                                              branch_count, anchor_id, sigma_phi, &
                                              da, dwleb, dxi0, dwbranch, &
                                              dw_xi, dw_f, dbranch_phi_adj, error)
      !> Raw surface adjoints, held fixed
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Folded weights of [[prepare_surface_weights]] for `acc`
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
      !> Rendered shapes, fixed length for a thread-safe message build
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

      ! d(-2 a w_a/xi) and d(-2 wleb w_w/xi) at fixed w_a, w_w, written as
      ! (dX - X dxi/xi)/xi to avoid 1/xi^2; check_surface_adjoint rejects a
      ! vanishing xi under a live area or weight adjoint
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

      ! d(w_a R^2 wleb) with the radius frozen, see the module TODO
      if (fold_switching) then
         do igrid = 1, ngrid
            if (abs(acc%w_a(igrid)) > seed_weight_tol) then
               r_own = radii(owner(igrid))
               dw_f(igrid) = dw_f(igrid) &
                             + acc%w_a(igrid)*r_own*r_own*dwleb(igrid)
            end if
         end do
      end if

      ! Last and on the folded `dw_xi`, as the primal pass on the folded `w_xi`
      call branch_phi_adj_tangent(branch_count, anchor_id, wbranch, wleb, xi0, &
                                  sigma_phi, eff%w_xi, dwbranch, dwleb, dxi0, dw_xi, &
                                  dbranch_phi_adj)

   end subroutine prepare_surface_weights_tangent

   !> Tangent of the branch-weight reverse pass, [[compute_branch_phi_adj]]
   !>
   !> - Same group walk, early exits and guard conditions as the primal
   !> - Per contiguous anchor group: `adj_m = adj_wleb_m factor_m`,
   !>   `mean = sum_m wbranch_m adj_m`,
   !>   `Phi_m = -wbranch_m (adj_m - mean)/sigma_phi`
   !> - `dmean` couples every point of a group, gated-off points included:
   !>   two passes per group
   !> - Serial over a group, never split one across threads, see module header
   !> - Per-point gate reads the moving `w_xi`, `wleb`, `wbranch`: primal's
   !>   branch on the primal's condition, zero on the else, no separate
   !>   second-order threshold
   !> - Early exits `sigma_phi <= seed_weight_tol` and no `branch_count > 1`
   !>   are frozen discrete choices
   !>
   !> @param[in]  branch_count     branches per grid point `(ngrid)`
   !> @param[in]  anchor_id        anchor group id per grid point `(ngrid)`
   !> @param[in]  wbranch          softmax branch weight per grid point `(ngrid)`
   !> @param[in]  wleb             Lebedev weight per grid point `(ngrid)`
   !> @param[in]  xi0              Gaussian width per grid point `(ngrid)`
   !> @param[in]  sigma_phi        softmax temperature
   !> @param[in]  w_xi             folded width adjoint, `eff%w_xi` `(ngrid)`
   !> @param[in]  dwbranch         directional tangent of `wbranch` `(ngrid)`
   !> @param[in]  dwleb            directional tangent of `wleb` `(ngrid)`
   !> @param[in]  dxi0             directional tangent of `xi0` `(ngrid)`
   !> @param[in]  dw_xi            directional tangent of the folded `w_xi` `(ngrid)`
   !> @param[out] dbranch_phi_adj  tangent of the branch-objective adjoint `(ngrid)`
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
      integer :: ngrid, igroup_cursor, igroup_start, igroup_end, group_size
      integer :: m_branch, im_grid
      !> Whether a further anchor group exists
      logical :: have_group
      !> Per-point branch adjoint and its tangent
      real(wp) :: adj_branch, dadj_branch
      !> Group reduction and its tangent
      real(wp) :: mean_adj_branch, dmean_adj_branch

      dbranch_phi_adj = 0.0_wp
      ngrid = size(branch_count)
      if (ngrid <= 0) return
      if (.not. any(branch_count > 1)) return
      if (sigma_phi <= seed_weight_tol) return

      igroup_cursor = 1
      do
         call next_branch_group(branch_count, anchor_id, igroup_cursor, &
                                igroup_start, igroup_end, have_group)
         if (.not. have_group) exit
         group_size = igroup_end - igroup_start + 1

         mean_adj_branch = 0.0_wp
         dmean_adj_branch = 0.0_wp
         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            call branch_point_adjoint(w_xi(im_grid), wleb(im_grid), xi0(im_grid), &
                                      wbranch(im_grid), adj_branch, &
                                      dw_xi(im_grid), dwleb(im_grid), dxi0(im_grid), &
                                      dwbranch(im_grid), dadj_branch)
            mean_adj_branch = mean_adj_branch + wbranch(im_grid)*adj_branch
            dmean_adj_branch = dmean_adj_branch + dwbranch(im_grid)*adj_branch &
                               + wbranch(im_grid)*dadj_branch
         end do

         ! Per-point pair recomputed, no scratch array: groups are 2-4 points
         do m_branch = 1, group_size
            im_grid = igroup_start + m_branch - 1
            call branch_point_adjoint(w_xi(im_grid), wleb(im_grid), xi0(im_grid), &
                                      wbranch(im_grid), adj_branch, &
                                      dw_xi(im_grid), dwleb(im_grid), dxi0(im_grid), &
                                      dwbranch(im_grid), dadj_branch)
            dbranch_phi_adj(im_grid) = &
               -(dwbranch(im_grid)*(adj_branch - mean_adj_branch) &
                 + wbranch(im_grid)*(dadj_branch - dmean_adj_branch))/sigma_phi
         end do
      end do

   end subroutine branch_phi_adj_tangent

end module moist_cavity_drop_derivatives_weights_tangent
