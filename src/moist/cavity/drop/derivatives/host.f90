!> Analytic DROP surface derivatives for arbitrary host parameters
!>
!> - Host supplies the spatial level-set jet and its parameter derivatives at
!>   the fixed projected point
!> - Bordered projection system differentiated twice, then the seed kernels of
!>   the nuclear Hessian
!> - Density directions move the surface even with every nuclear direction zero
!> - Branch weight couples the points of an anchor group: held fixed per
!>   point, motion added over the grid by [[drop_host_branch_derivatives]]
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_host
   use moist_cavity_drop_derivatives_kernel, only: build_seed_state, apply_seed, &
      & apply_seed_tangent, drop_seed_result_type, drop_seed_state_tangent_type, &
      & drop_seed_input_tangent_type, seed_state_ok, seed_weight_tol, &
      & next_branch_group, max_branch_group_size
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type
   implicit none(type, external)
contains
   !> First and second host-parameter derivatives of one projected surface point
   !>
   !> - Jets are partial derivatives of the level set `S`, not the density, at
   !>   the fixed projected point
   !> - Cartesian tensors packed by spatial order, each in Fortran order
   !> - `jet`: spatial orders 0 to 4 at 1, 2:4, 5:13, 14:40, 41:121
   !> - `jet1`: first parameter derivative, spatial orders 0 to 3
   !> - `jet2`: second parameter derivative, spatial orders 0 to 2
   !> - Rows of `d1` and `d2`: `(x, y, z, xi, f)`
   !> - Row `f` from iSwiG on the rigid anchor sphere, moved by `dirs` alone
   !> - Branch weight held fixed: row `xi` of a multi-branch grid is completed
   !>   by [[drop_host_branch_derivatives]] once every point is known
   !> - `1 <= igrid <= ngrid` and `n >= 1`
   !> - `d1`, `d2` undefined on failure
   !>
   !> @param[in]  self   DROP cavity instance holding a projected grid
   !> @param[in]  igrid  grid point index
   !> @param[in]  dirs   nuclear/host directions `(3, nsph, n)`
   !> @param[in]  jet    level-set jet `(121)`
   !> @param[in]  jet1   first parameter derivative of the jet `(40, n)`
   !> @param[in]  jet2   second parameter derivative of the jet `(13, n, n)`
   !> @param[out] d1     first derivatives of `(x, y, z, xi, f)`, `(5, n)`
   !> @param[out] d2     second derivatives of `(x, y, z, xi, f)`, `(5, n, n)`
   !> @param[out] error  error object, allocated on failure
   module subroutine drop_host_point_derivatives(self, igrid, dirs, jet, jet1, jet2, d1, d2, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Grid point index
      integer, intent(in) :: igrid
      !> Directions, level-set jet and its first and second parameter derivatives
      real(wp), intent(in) :: dirs(:, :, :), jet(:), jet1(:, :), jet2(:, :, :)
      !> First and second derivatives of the surface point
      real(wp), intent(out) :: d1(:, :), d2(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Seed state of the grid point
      type(drop_seed_state_type) :: state
      !> Factored bordered projection system
      type(drop_kkt_factor_type) :: fac
      !> Point scratch of the iSwiG evaluation
      type(drop_point_scratch_type) :: pt
      !> First-order seed result per direction
      type(drop_seed_result_type), allocatable :: res(:)
      !> Seed-state tangent per direction
      type(drop_seed_state_tangent_type), allocatable :: ds(:)
      !> Seed input tangent of the current direction
      type(drop_seed_input_tangent_type) :: di
      !> Second-order seed result
      type(drop_seed_result_type) :: ddres
      !> First- and second-order bordered solutions; parameter derivatives of
      !> the level-set gradient, Hessian and third derivative
      real(wp), allocatable :: x(:, :), xx(:, :), dg(:, :), dh(:, :, :), dt(:, :, :, :)
      !> iSwiG second-derivative block over the influence set
      real(wp), allocatable :: swi2(:, :, :, :)
      !> Atom ids of the influence set, owner first
      integer, allocatable :: idx(:)
      !> Lagrangian Hessian and its tangent, second-order gradient and Hessian
      !> terms, fourth spatial derivative
      real(wp) :: hlag(3, 3), ddg(3), ddh(3, 3), dhlag(3, 3), fourth(3, 3, 3, 3)
      !> iSwiG value, owner gradient row and width derivative; second-order
      !> constraint term
      real(wp) :: f0, own(3), dxi, dpartial0
      !> Direction count, loop indices, owner sphere, seed status, influence-set size
      integer :: n, p, q, k, a, b, owner, status, nswi
      !> Routine name for diagnostics
      character(len=*), parameter :: context = 'drop_host_point_derivatives'

      n = size(dirs, 3)
      if (igrid < 1 .or. igrid > self%ngrid .or. n < 1) then
         call fatal_error(error, context//': invalid grid point or direction count')
         return
      end if
      if (any(shape(dirs) /= [3, self%nsph, n]) .or. size(jet) /= 121 .or. &
          any(shape(jet1) /= [40, n]) .or. any(shape(jet2) /= [13, n, n]) .or. &
          any(shape(d1) /= [5, n]) .or. any(shape(d2) /= [5, n, n])) then
         call fatal_error(error, context//': inconsistent host derivative shapes')
         return
      end if
      owner = self%owner(igrid)
      state%lsf1_r = jet(2:4)
      state%lsf2_rr = reshape(jet(5:13), [3, 3])
      state%lsf3_rrr = reshape(jet(14:40), [3, 3, 3])
      fourth = reshape(jet(41:121), [3, 3, 3, 3])
      state%lambda_val = self%lambda0(igrid)
      call fill_seed_state(self, igrid, .false., state)
      call build_seed_state(state, self%f_crit, self%f_foc, self%f_wleb, &
                            self%param%wleb_prune_level > 0, status)
      if (status /= seed_state_ok) then
         call fatal_error(error, context//': degenerate projection derivative state')
         return
      end if
      hlag = -state%lambda_val*state%lsf2_rr
      do k = 1, 3
         hlag(k, k) = hlag(k, k) + self%param%phi_alpha
      end do
      call fac%factor(hlag, state%lsf1_r, context, error, igrid)
      if (allocated(error)) return
      allocate (x(4, n), xx(4, n), dg(3, n), dh(3, 3, n), dt(3, 3, 3, n), res(n), ds(n))
      dg = jet1(2:4, :)
      dh = reshape(jet1(5:13, :), [3, 3, n])
      dt = reshape(jet1(14:40, :), [3, 3, 3, n])
      x(1:3, :) = state%lambda_val*dg + self%param%phi_alpha*dirs(:, owner, :)
      x(4, :) = -jet1(1, :)
      call fac%solve(x, context, error, igrid)
      if (allocated(error)) return
      d1 = 0.0_wp
      d2 = 0.0_wp
      do p = 1, n
         call apply_seed(state, dg(:, p), dh(:, :, p), x(1:3, p), x(4, p), res(p), ds(p))
         d1(1:3, p) = x(1:3, p)
         d1(4, p) = res(p)%dxi
      end do
      do q = 1, n
         di = drop_seed_input_tangent_type()
         di%dlsf1_r = res(q)%dg
         di%dlsf2_rr = res(q)%dH
         di%dlsf3_rrr = dt(:, :, :, q)
         do k = 1, 3
            di%dlsf3_rrr = di%dlsf3_rrr + fourth(:, :, :, k)*x(k, q)
         end do
         di%dlambda_val = x(4, q)
         di%dcpjac_scal0 = res(q)%dJ
         di%dw_f0 = res(q)%dw_f
         di%dwleb = res(q)%dwleb
         di%dxi0 = res(q)%dxi
         dhlag = -x(4, q)*state%lsf2_rr - state%lambda_val*res(q)%dH
         do p = 1, n
            ddg = jet2(2:4, p, q) + matmul(dh(:, :, p), x(1:3, q))
            dpartial0 = jet2(1, p, q) + dot_product(dg(:, p), x(1:3, q))
            xx(1:3, p) = x(4, q)*dg(:, p) + state%lambda_val*ddg &
               & - matmul(dhlag, x(1:3, p)) + res(q)%dg*x(4, p)
            xx(4, p) = -dpartial0 - dot_product(res(q)%dg, x(1:3, p))
         end do
         call fac%solve(xx, context, error, igrid)
         if (allocated(error)) return
         do p = 1, n
            ddg = jet2(2:4, p, q) + matmul(dh(:, :, p), x(1:3, q))
            ddh = reshape(jet2(5:13, p, q), [3, 3])
            do k = 1, 3
               ddh = ddh + dt(:, :, k, p)*x(k, q)
            end do
            call apply_seed_tangent(state, ds(q), di, res(q), &
               & x(1:3, p), x(4, p), ddg, ddh, xx(1:3, p), xx(4, p), res(p), ds(p), ddres)
            d2(1:3, p, q) = xx(1:3, p)
            d2(4, p, q) = ddres%dxi
         end do
      end do

      ! iSwiG on the rigid anchor sphere, independent of density
      call pt%init(self%nsph, self%iswig, .false., .false., .false.)
      call self%iswig%swi_collect(self%anchorxyz(:, igrid), owner, self%anchor_xi0(igrid), f0, pt%iswig_work)
      call self%iswig%swi1_rA_sparse(pt%iswig_work, pt%swi_rows, own, dxi)
      nswi = pt%iswig_work%n_nb + 1
      allocate (idx(nswi), swi2(3, nswi, 3, nswi))
      call self%iswig%swi2_rArB_block(pt%iswig_work, nswi, idx, swi2)
      do p = 1, n
         d1(5, p) = dot_product(own, dirs(:, owner, p))
         do k = 1, pt%iswig_work%n_nb
            d1(5, p) = d1(5, p) + dot_product(pt%swi_rows(:, k), dirs(:, pt%iswig_work%idx(k), p))
         end do
         do q = 1, n
            do a = 1, nswi
               do b = 1, nswi
                  d2(5, p, q) = d2(5, p, q) &
                     & + dot_product(dirs(:, idx(a), p), matmul(swi2(:, a, :, b), dirs(:, idx(b), q)))
               end do
            end do
         end do
      end do
      call pt%destroy()
   end subroutine drop_host_point_derivatives

   !> Add the branch-weight motion to the host derivatives of the whole grid
   !>
   !> - Completes row `xi` of [[drop_host_point_derivatives]], which holds the
   !>   softmax branch weight `p` of each point fixed
   !> - `xi = xi_frozen (p/p0)^(-1/2)`: the Lebedev weight is `p` times a
   !>   point-local factor and the width scales with its inverse root
   !> - First order: `-0.5 xi dp/p`
   !> - Second order: `-0.5 (dxi_p dp_q + dxi_q dp_p)/p`
   !>   `+ xi (0.75 dp_p dp_q/p^2 - 0.5 ddp_pq/p)`, `dxi` the frozen first order
   !> - `p` from the softmax over `Phi = 0.5 alpha |r - anchor|^2` of the group
   !> - Anchor rigid with its owner: `dPhi = phi1_r . (dr - v_owner)`,
   !>   `ddPhi = alpha (dr_p - v_p) . (dr_q - v_q) + phi1_r . ddr_pq`
   !> - Rows `x, y, z, f` do not read the branch weight and stay as they are
   !> - Guards as in `drop_surface_tangent_core`: no motion below the weight
   !>   thresholds
   !> - Single-branch grid: returns without touching `d1`, `d2`
   !> - Groups from [[next_branch_group]], serial
   !> - `d1`, `d2` untouched on failure
   !>
   !> @param[in]     self   DROP cavity instance holding a projected grid
   !> @param[in]     dirs   nuclear/host directions `(3, nsph, n)`
   !> @param[in,out] d1     first derivatives of `(x, y, z, xi, f)`, `(5, ngrid, n)`
   !> @param[in,out] d2     second derivatives of `(x, y, z, xi, f)`, `(5, ngrid, n, n)`
   !> @param[out]    error  error object, allocated on failure
   module subroutine drop_host_branch_derivatives(self, dirs, d1, d2, error)
      !> DROP cavity instance
      class(cavity_type_drop), intent(in) :: self
      !> Directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> First derivatives of the surface points
      real(wp), intent(inout) :: d1(:, :, :)
      !> Second derivatives of the surface points
      real(wp), intent(inout) :: d2(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Objective values and softmax weights of one group
      real(wp), allocatable :: phi(:), wgt(:)
      !> First derivatives of the objective and of the weights
      real(wp), allocatable :: dphi(:, :), dwgt(:, :)
      !> Second derivatives of the objective and of the weights
      real(wp), allocatable :: ddphi(:, :, :), ddwgt(:, :, :)
      !> Point motion relative to the rigid anchor, per direction
      real(wp), allocatable :: rel(:, :)
      !> Frozen first derivative of the width
      real(wp), allocatable :: dxi(:)
      !> Objective gradient at the projected point
      real(wp) :: phi1(3)
      !> Width and branch weight of the point
      real(wp) :: xi, wb
      !> Direction count, largest group, group walk bookkeeping
      integer :: n, nmax, cursor, first, last, nb
      !> Loop indices and grid point
      integer :: m, p, q, im
      !> Whether a further group exists
      logical :: found
      !> Routine name for diagnostics
      character(len=*), parameter :: context = 'drop_host_branch_derivatives'

      n = size(dirs, 3)
      if (n < 1 .or. self%ngrid < 1) then
         call fatal_error(error, context//': invalid grid or direction count')
         return
      end if
      if (any(shape(dirs) /= [3, self%nsph, n]) .or. any(shape(d1) /= [5, self%ngrid, n]) .or. &
          any(shape(d2) /= [5, self%ngrid, n, n])) then
         call fatal_error(error, context//': inconsistent host derivative shapes')
         return
      end if
      if (.not. allocated(self%branch_count) .or. .not. allocated(self%anchor_id)) return
      if (.not. any(self%branch_count(1:self%ngrid) > 1)) return

      nmax = max_branch_group_size(self%branch_count(1:self%ngrid), self%anchor_id(1:self%ngrid))
      allocate (phi(nmax), wgt(nmax), dphi(n, nmax), dwgt(n, nmax), ddphi(n, n, nmax), &
                ddwgt(n, n, nmax), rel(3, n), dxi(n))

      cursor = 1
      do
         call next_branch_group(self%branch_count(1:self%ngrid), self%anchor_id(1:self%ngrid), &
                                cursor, first, last, found)
         if (.not. found) exit
         nb = last - first + 1

         do m = 1, nb
            im = first + m - 1
            phi(m) = self%phi0(im)
            phi1 = self%param%phi_alpha*(self%xyz(:, im) - self%anchorxyz(:, im))
            rel = d1(1:3, im, :) - dirs(:, self%owner(im), :)
            do q = 1, n
               dphi(q, m) = dot_product(phi1, rel(:, q))
               do p = 1, n
                  ddphi(p, q, m) = self%param%phi_alpha*dot_product(rel(:, p), rel(:, q)) &
                     & + dot_product(phi1, d2(1:3, im, p, q))
               end do
            end do
         end do

         call self%branch_weight%weights_hess(phi(1:nb), dphi(:, 1:nb), ddphi(:, :, 1:nb), &
            & wgt(1:nb), dwgt(:, 1:nb), ddwgt(:, :, 1:nb))

         do m = 1, nb
            im = first + m - 1
            wb = self%wbranch(im)
            if (wb <= tiny(1.0_wp) .or. self%wleb(im) <= seed_weight_tol) cycle
            xi = self%xi0(im)
            dxi = d1(4, im, :)
            do q = 1, n
               do p = 1, n
                  d2(4, im, p, q) = d2(4, im, p, q) &
                     & - 0.5_wp*(dxi(p)*dwgt(q, m) + dxi(q)*dwgt(p, m))/wb &
                     & + xi*(0.75_wp*dwgt(p, m)*dwgt(q, m)/(wb*wb) - 0.5_wp*ddwgt(p, q, m)/wb)
               end do
            end do
            d1(4, im, :) = dxi - 0.5_wp*xi*dwgt(:, m)/wb
         end do
      end do
   end subroutine drop_host_branch_derivatives
end submodule moist_cavity_drop_derivatives_host
