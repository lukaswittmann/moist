!> Second-order contractions of the Gaussian PCM interaction matrix
!>
!> The nuclear Hessian of a PCM energy needs two things from the matrix that
!> the gradient does not: the matrix moved along a surface tangent, applied to
!> the charges, and the response of the gradient's surface weights along that
!> tangent. Both are pair sums over the grid, so both fuse a whole batch of
!> directions into one pass over the pairs: the pair kernel is evaluated once
!> and contracted with every direction of the batch.
!>
!> The tangent enters through the three independent scalars of the kernel,
!> `(xi_i, xi_j, r2)`, with `dr2 = 2 (r_i - r_j) . (d_i - d_j)` for a pair
!> moved by `(d_i, d_j)`. The saturated regime keeps only the `r2` channel.
submodule(moist_model_component_pcm_amat) moist_model_component_pcm_amat_hessian
   use mctc_env, only: fatal_error
   use moist_model_component_pcm_amat_kernel, only: pcm_amat_far_grad, pcm_amat_far_hess, &
      pcm_amat_near_grad, pcm_amat_near_hess, pcm_amat_diag_grad, pcm_amat_diag_hess
   implicit none(type, external)

contains

   !> Validate the tangent batch against the surface and the charges
   !>
   !> @param[in]  ngrid   Number of grid points
   !> @param[in]  q       Charges (ngrid)
   !> @param[in]  d_xi    Width tangents (ngrid, ndir)
   !> @param[in]  d_f     Switching-factor tangents (ngrid, ndir)
   !> @param[in]  d_xyz   Position tangents (3, ngrid, ndir)
   !> @param[in]  context Calling routine, used to prefix the diagnostics
   !> @param[out] ndir    Number of directions of the batch
   !> @param[out] error   Error handling
   subroutine validate_tangent_batch(ngrid, q, d_xi, d_f, d_xyz, context, ndir, error)
      !> Number of grid points
      integer, intent(in) :: ngrid
      !> Charges
      real(wp), intent(in) :: q(:)
      !> Width and switching-factor tangents
      real(wp), intent(in) :: d_xi(:, :), d_f(:, :)
      !> Position tangents
      real(wp), intent(in) :: d_xyz(:, :, :)
      !> Calling routine
      character(len=*), intent(in) :: context
      !> Number of directions of the batch
      integer, intent(out) :: ndir
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ndir = size(d_xi, 2)
      if (size(q) /= ngrid) then
         call fatal_error(error, context//": charge vector does not match the surface")
      else if (ndir < 1) then
         call fatal_error(error, context//": no direction supplied")
      else if (size(d_xi, 1) /= ngrid .or. size(d_f, 1) /= ngrid .or. size(d_f, 2) /= ndir .or. &
               size(d_xyz, 1) /= 3 .or. size(d_xyz, 2) /= ngrid .or. size(d_xyz, 3) /= ndir) then
         call fatal_error(error, context//": surface tangent shape mismatch")
      end if
   end subroutine validate_tangent_batch

   !> Apply the matrix moved along a batch of surface tangents to the charges
   !>
   !> @param[in]  xi     Gaussian widths (ngrid)
   !> @param[in]  f      Gaussian switching factors (ngrid)
   !> @param[in]  xyz    Surface positions (3, ngrid)
   !> @param[in]  q      Charges the moved matrix is applied to (ngrid)
   !> @param[in]  d_xi   Width tangents (ngrid, ndir)
   !> @param[in]  d_f    Switching-factor tangents (ngrid, ndir)
   !> @param[in]  d_xyz  Position tangents (3, ngrid, ndir)
   !> @param[out] daq    `(dA . Gamma_v) q` per direction (ngrid, ndir)
   !> @param[out] error  Error handling
   module subroutine pcm_amat_tangent_apply(xi, f, xyz, q, d_xi, d_f, d_xyz, daq, error)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Charges the moved matrix is applied to
      real(wp), intent(in) :: q(:)
      !> Width and switching-factor tangents
      real(wp), intent(in) :: d_xi(:, :), d_f(:, :)
      !> Position tangents
      real(wp), intent(in) :: d_xyz(:, :, :)
      !> Moved matrix applied to the charges, per direction
      real(wp), intent(out) :: daq(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface-point, direction and axis indices with their extents
      integer :: i, j, v, iaxis, ngrid, ndir
      !> Row-point width, position, saturation bound and charge of the column point
      real(wp) :: xi_i, xyz_i(3), bound_i, q_j
      !> Pair displacement and squared separation
      real(wp) :: u(3), r2
      !> Kernel value and its three derivative channels
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      !> Diagonal kernel value and its two derivative channels
      real(wp) :: a_diag, a_diag_xi, a_diag_f
      !> Tangent of the squared separation over two
      real(wp) :: dd
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)
      !> Direction-major copies of the tangents and the result
      real(wp), allocatable :: dxi_t(:, :), df_t(:, :), dxyz_t(:, :, :), daq_t(:, :)
      !> Per-thread accumulator over the directions
      real(wp), allocatable :: acc(:)

      daq = 0.0_wp
      call validate_pcm_surface(xi, f, xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)
      call validate_tangent_batch(ngrid, q, d_xi, d_f, d_xyz, "pcm_amat_tangent_apply", &
                                  ndir, error)
      if (allocated(error)) return
      if (size(daq, 1) /= ngrid .or. size(daq, 2) /= ndir) then
         call fatal_error(error, "pcm_amat_tangent_apply: output shape mismatch")
         return
      end if

      allocate (bound(ngrid))
      call saturation_bounds(xi, bound)

      ! Direction-major layout: the inner loop runs over the directions of one
      ! pair, so the tangents of a point must be contiguous over the batch
      allocate (dxi_t(ndir, ngrid), df_t(ndir, ngrid), dxyz_t(ndir, 3, ngrid), daq_t(ndir, ngrid))
      do i = 1, ngrid
         dxi_t(:, i) = d_xi(i, :)
         df_t(:, i) = d_f(i, :)
         do iaxis = 1, 3
            dxyz_t(:, iaxis, i) = d_xyz(iaxis, i, :)
         end do
      end do

      !$omp parallel default(none) &
      !$omp shared(xi, f, xyz, q, bound, dxi_t, df_t, dxyz_t, daq_t, ngrid, ndir) &
      !$omp private(i, j, v, xi_i, xyz_i, bound_i, q_j, u, r2, a, a_xi_i, a_xi_j, a_r2, &
      !$omp         a_diag, a_diag_xi, a_diag_f, dd, acc)
      allocate (acc(ndir))
      !$omp do schedule(dynamic, 8)
      do i = 1, ngrid
         xi_i = xi(i)
         xyz_i = xyz(:, i)
         bound_i = bound(i)

         call pcm_amat_diag_grad(xi_i, f(i), a_diag, a_diag_xi, a_diag_f)
         do v = 1, ndir
            acc(v) = (a_diag_xi*dxi_t(v, i) + a_diag_f*df_t(v, i))*q(i)
         end do

         do j = 1, ngrid
            if (j == i) cycle
            u = xyz_i - xyz(:, j)
            r2 = max(u(1)*u(1) + u(2)*u(2) + u(3)*u(3), r2_floor)
            q_j = q(j)
            if (r2 >= bound_i + bound(j)) then
               call pcm_amat_far_grad(r2, a, a_r2)
               do v = 1, ndir
                  dd = u(1)*(dxyz_t(v, 1, i) - dxyz_t(v, 1, j)) &
                       + u(2)*(dxyz_t(v, 2, i) - dxyz_t(v, 2, j)) &
                       + u(3)*(dxyz_t(v, 3, i) - dxyz_t(v, 3, j))
                  acc(v) = acc(v) + 2.0_wp*a_r2*dd*q_j
               end do
            else
               call pcm_amat_near_grad(xi_i, xi(j), r2, a, a_xi_i, a_xi_j, a_r2)
               do v = 1, ndir
                  dd = u(1)*(dxyz_t(v, 1, i) - dxyz_t(v, 1, j)) &
                       + u(2)*(dxyz_t(v, 2, i) - dxyz_t(v, 2, j)) &
                       + u(3)*(dxyz_t(v, 3, i) - dxyz_t(v, 3, j))
                  acc(v) = acc(v) + (a_xi_i*dxi_t(v, i) + a_xi_j*dxi_t(v, j) &
                                     + 2.0_wp*a_r2*dd)*q_j
               end do
            end if
         end do

         daq_t(:, i) = acc
      end do
      !$omp end do
      !$omp end parallel

      do i = 1, ngrid
         daq(i, :) = daq_t(:, i)
      end do
   end subroutine pcm_amat_tangent_apply

   !> Response of the surface-variable weights along a batch of surface tangents
   !>
   !> The weights of [[pcm_amat_surface_weights]] at `q1 = q2 = q` are
   !> `w = d(q^T A q)/d(Gamma)` at fixed charges. Along a direction on which the
   !> surface moves by `Gamma_v` and the charges by `dq`, they respond by
   !>
   !>     dw = (d2(q^T A q)/dGamma2) Gamma_v + 2 w(dq, q)
   !>
   !> where `w(dq, q)` is the bilinear weight of the two charge vectors. Both
   !> terms are formed in the one pair pass. The charge response is contracted
   !> without a smallness gate on the charge product: a pair whose base
   !> product vanishes still carries a first-order response.
   !>
   !> @param[in]  xi      Gaussian widths (ngrid)
   !> @param[in]  f       Gaussian switching factors (ngrid)
   !> @param[in]  xyz     Surface positions (3, ngrid)
   !> @param[in]  q       Charges (ngrid)
   !> @param[in]  dq      Charge response per direction (ngrid, ndir)
   !> @param[in]  d_xi    Width tangents (ngrid, ndir)
   !> @param[in]  d_f     Switching-factor tangents (ngrid, ndir)
   !> @param[in]  d_xyz   Position tangents (3, ngrid, ndir)
   !> @param[out] dw_xi   Width-weight response (ngrid, ndir)
   !> @param[out] dw_f    Switching-factor-weight response (ngrid, ndir)
   !> @param[out] dw_xyz  Position-weight response (3, ngrid, ndir)
   !> @param[out] error   Error handling
   module subroutine pcm_amat_surface_weights_response(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, &
                                                       dw_xi, dw_f, dw_xyz, error)
      !> Gaussian widths
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors
      real(wp), intent(in) :: f(:)
      !> Surface positions
      real(wp), intent(in) :: xyz(:, :)
      !> Charges
      real(wp), intent(in) :: q(:)
      !> Charge response per direction
      real(wp), intent(in) :: dq(:, :)
      !> Width and switching-factor tangents
      real(wp), intent(in) :: d_xi(:, :), d_f(:, :)
      !> Position tangents
      real(wp), intent(in) :: d_xyz(:, :, :)
      !> Width-weight and switching-factor-weight responses
      real(wp), intent(out) :: dw_xi(:, :), dw_f(:, :)
      !> Position-weight response
      real(wp), intent(out) :: dw_xyz(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Surface-point, direction and axis indices with their extents
      integer :: i, j, v, iaxis, ngrid, ndir
      !> Row-point width, position, saturation bound and charges of the pair
      real(wp) :: xi_i, xyz_i(3), bound_i, q_i, q_j
      !> Base and response charge products of the pair
      real(wp) :: qsym, dqsym
      !> Pair displacement and squared separation
      real(wp) :: u(3), r2
      !> Kernel value, first derivatives, and second derivatives
      real(wp) :: a, a_xi_i, a_xi_j, a_r2
      real(wp) :: a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2
      !> Diagonal kernel value, first derivatives, and second derivatives
      real(wp) :: a_diag, a_diag_xi, a_diag_f, a_diag_xx, a_diag_xf, a_diag_ff
      !> Tangent of the squared separation, the moved kernel derivative, and
      !> the coefficient of the pair displacement
      real(wp) :: dr2, da_r2, coef
      !> Per-point saturation bounds
      real(wp), allocatable :: bound(:)
      !> Direction-major copies of the tangents and the charge response
      real(wp), allocatable :: dxi_t(:, :), df_t(:, :), dxyz_t(:, :, :), dq_t(:, :)
      !> Direction-major results
      real(wp), allocatable :: dw_xi_t(:, :), dw_f_t(:, :), dw_xyz_t(:, :, :)
      !> Per-thread accumulators over the directions
      real(wp), allocatable :: acc_xi(:), acc_xyz(:, :)

      dw_xi = 0.0_wp
      dw_f = 0.0_wp
      dw_xyz = 0.0_wp
      call validate_pcm_surface(xi, f, xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)
      call validate_tangent_batch(ngrid, q, d_xi, d_f, d_xyz, &
                                  "pcm_amat_surface_weights_response", ndir, error)
      if (allocated(error)) return
      if (size(dq, 1) /= ngrid .or. size(dq, 2) /= ndir) then
         call fatal_error(error, "pcm_amat_surface_weights_response: charge response shape mismatch")
         return
      end if
      if (size(dw_xi, 1) /= ngrid .or. size(dw_xi, 2) /= ndir .or. &
          size(dw_f, 1) /= ngrid .or. size(dw_f, 2) /= ndir .or. &
          size(dw_xyz, 1) /= 3 .or. size(dw_xyz, 2) /= ngrid .or. size(dw_xyz, 3) /= ndir) then
         call fatal_error(error, "pcm_amat_surface_weights_response: output shape mismatch")
         return
      end if

      allocate (bound(ngrid))
      call saturation_bounds(xi, bound)

      allocate (dxi_t(ndir, ngrid), df_t(ndir, ngrid), dxyz_t(ndir, 3, ngrid), dq_t(ndir, ngrid))
      allocate (dw_xi_t(ndir, ngrid), dw_f_t(ndir, ngrid), dw_xyz_t(ndir, 3, ngrid))
      do i = 1, ngrid
         dxi_t(:, i) = d_xi(i, :)
         df_t(:, i) = d_f(i, :)
         dq_t(:, i) = dq(i, :)
         do iaxis = 1, 3
            dxyz_t(:, iaxis, i) = d_xyz(iaxis, i, :)
         end do
      end do

      !$omp parallel default(none) &
      !$omp shared(xi, f, xyz, q, bound, dxi_t, df_t, dxyz_t, dq_t, dw_xi_t, dw_f_t, dw_xyz_t, &
      !$omp        ngrid, ndir) &
      !$omp private(i, j, v, xi_i, xyz_i, bound_i, q_i, q_j, qsym, dqsym, u, r2, &
      !$omp         a, a_xi_i, a_xi_j, a_r2, a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2, &
      !$omp         a_diag, a_diag_xi, a_diag_f, a_diag_xx, a_diag_xf, a_diag_ff, &
      !$omp         dr2, da_r2, coef, acc_xi, acc_xyz)
      allocate (acc_xi(ndir), acc_xyz(ndir, 3))
      !$omp do schedule(dynamic, 8)
      do i = 1, ngrid
         xi_i = xi(i)
         xyz_i = xyz(:, i)
         bound_i = bound(i)
         q_i = q(i)

         ! Self term: second derivatives along the tangent, charge response
         ! through the bilinear weight
         call pcm_amat_diag_hess(xi_i, f(i), a_diag, a_diag_xi, a_diag_f, &
                                 a_diag_xx, a_diag_xf, a_diag_ff)
         do v = 1, ndir
            acc_xi(v) = q_i*q_i*(a_diag_xx*dxi_t(v, i) + a_diag_xf*df_t(v, i)) &
                        + 2.0_wp*dq_t(v, i)*q_i*a_diag_xi
            dw_f_t(v, i) = q_i*q_i*(a_diag_xf*dxi_t(v, i) + a_diag_ff*df_t(v, i)) &
                           + 2.0_wp*dq_t(v, i)*q_i*a_diag_f
         end do
         acc_xyz = 0.0_wp

         do j = 1, ngrid
            if (j == i) cycle
            u = xyz_i - xyz(:, j)
            r2 = max(u(1)*u(1) + u(2)*u(2) + u(3)*u(3), r2_floor)
            q_j = q(j)
            qsym = 2.0_wp*q_i*q_j
            if (r2 >= bound_i + bound(j)) then
               ! Saturated: only the separation channel survives
               call pcm_amat_far_hess(r2, a, a_r2, a_r2r2)
               do v = 1, ndir
                  dr2 = 2.0_wp*(u(1)*(dxyz_t(v, 1, i) - dxyz_t(v, 1, j)) &
                                + u(2)*(dxyz_t(v, 2, i) - dxyz_t(v, 2, j)) &
                                + u(3)*(dxyz_t(v, 3, i) - dxyz_t(v, 3, j)))
                  dqsym = 2.0_wp*(dq_t(v, i)*q_j + dq_t(v, j)*q_i)
                  coef = 2.0_wp*(qsym*a_r2r2*dr2 + dqsym*a_r2)
                  acc_xyz(v, 1) = acc_xyz(v, 1) + coef*u(1) &
                                  + 2.0_wp*qsym*a_r2*(dxyz_t(v, 1, i) - dxyz_t(v, 1, j))
                  acc_xyz(v, 2) = acc_xyz(v, 2) + coef*u(2) &
                                  + 2.0_wp*qsym*a_r2*(dxyz_t(v, 2, i) - dxyz_t(v, 2, j))
                  acc_xyz(v, 3) = acc_xyz(v, 3) + coef*u(3) &
                                  + 2.0_wp*qsym*a_r2*(dxyz_t(v, 3, i) - dxyz_t(v, 3, j))
               end do
            else
               call pcm_amat_near_hess(xi_i, xi(j), r2, a, a_xi_i, a_xi_j, a_r2, &
                                       a_ii, a_ij, a_jj, a_ir2, a_jr2, a_r2r2)
               do v = 1, ndir
                  dr2 = 2.0_wp*(u(1)*(dxyz_t(v, 1, i) - dxyz_t(v, 1, j)) &
                                + u(2)*(dxyz_t(v, 2, i) - dxyz_t(v, 2, j)) &
                                + u(3)*(dxyz_t(v, 3, i) - dxyz_t(v, 3, j)))
                  dqsym = 2.0_wp*(dq_t(v, i)*q_j + dq_t(v, j)*q_i)
                  acc_xi(v) = acc_xi(v) &
                              + qsym*(a_ii*dxi_t(v, i) + a_ij*dxi_t(v, j) + a_ir2*dr2) &
                              + dqsym*a_xi_i
                  da_r2 = a_ir2*dxi_t(v, i) + a_jr2*dxi_t(v, j) + a_r2r2*dr2
                  coef = 2.0_wp*(qsym*da_r2 + dqsym*a_r2)
                  acc_xyz(v, 1) = acc_xyz(v, 1) + coef*u(1) &
                                  + 2.0_wp*qsym*a_r2*(dxyz_t(v, 1, i) - dxyz_t(v, 1, j))
                  acc_xyz(v, 2) = acc_xyz(v, 2) + coef*u(2) &
                                  + 2.0_wp*qsym*a_r2*(dxyz_t(v, 2, i) - dxyz_t(v, 2, j))
                  acc_xyz(v, 3) = acc_xyz(v, 3) + coef*u(3) &
                                  + 2.0_wp*qsym*a_r2*(dxyz_t(v, 3, i) - dxyz_t(v, 3, j))
               end do
            end if
         end do

         dw_xi_t(:, i) = acc_xi
         dw_xyz_t(:, :, i) = acc_xyz
      end do
      !$omp end do
      !$omp end parallel

      do i = 1, ngrid
         dw_xi(i, :) = dw_xi_t(:, i)
         dw_f(i, :) = dw_f_t(:, i)
         do iaxis = 1, 3
            dw_xyz(iaxis, i, :) = dw_xyz_t(:, iaxis, i)
         end do
      end do
   end subroutine pcm_amat_surface_weights_response

end submodule moist_model_component_pcm_amat_hessian
