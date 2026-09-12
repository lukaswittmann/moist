!> Analytic DROP surface derivatives for arbitrary host parameters.
!>
!> The host supplies the spatial level-set jet and its parameter derivatives
!> at the fixed projected point. Differentiate the bordered projection system
!> twice, then use the same seed kernels as the nuclear Hessian. In particular,
!> density directions move the surface even when every nuclear direction is zero.
submodule(moist_cavity_drop) moist_cavity_drop_derivatives_host
   use moist_cavity_drop_derivatives_kernel, only: build_seed_state, apply_seed, &
      & apply_seed_tangent, drop_seed_result_type, drop_seed_state_tangent_type, &
      & drop_seed_input_tangent_type, seed_state_ok
   use moist_cavity_drop_derivatives_seeds, only: drop_kkt_factor_type
   implicit none(type, external)
contains
   module subroutine drop_host_point_derivatives(self, igrid, dirs, jet, jet1, jet2, d1, d2, error)
      class(cavity_type_drop), intent(in) :: self
      integer, intent(in) :: igrid
      real(wp), intent(in) :: dirs(:, :, :), jet(:), jet1(:, :), jet2(:, :, :)
      real(wp), intent(out) :: d1(:, :), d2(:, :, :)
      type(error_type), allocatable, intent(out) :: error
      type(drop_seed_state_type) :: state
      type(drop_kkt_factor_type) :: fac
      type(drop_point_scratch_type) :: pt
      type(drop_seed_result_type), allocatable :: res(:)
      type(drop_seed_state_tangent_type), allocatable :: ds(:)
      type(drop_seed_input_tangent_type) :: di
      type(drop_seed_result_type) :: ddres
      real(wp), allocatable :: x(:, :), xx(:, :), dg(:, :), dh(:, :, :), dt(:, :, :, :)
      real(wp), allocatable :: swi2(:, :, :, :)
      integer, allocatable :: idx(:)
      real(wp) :: hlag(3, 3), ddg(3), ddh(3, 3), dhlag(3, 3), fourth(3, 3, 3, 3)
      real(wp) :: f0, own(3), dxi, dpartial0
      integer :: n, p, q, k, a, b, owner, status, nswi
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
      if (any(self%branch_count /= 1)) then
         call fatal_error(error, context//': Hessians require single-branch DROP projection')
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

      ! iSwiG is evaluated on the rigid anchor sphere, independent of density.
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
end submodule moist_cavity_drop_derivatives_host
