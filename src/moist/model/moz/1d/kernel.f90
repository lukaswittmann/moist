!> Per-k kernels of the 1D UV equation
!>
!> - Block products h = W c P and transpose over vv-style k lanes
!>   - UV column (i - 1) nv + a: solute site i, solvent site a
!> - W(k): solute intramolecular correlation, packed symmetric (kb, nu(nu+1)/2)
!> - P(k): solvent response P_ba = w_ba + rho_b h_ba, full (kb, nv, nv)
!> - Charge-moment correction of P at small k
!> - Pointwise energy-response kernels shared with 3D: response_kernel.f90
!>
!> - Products follow h_uv = w_u c_uv (w_vv + rho h_vv) of the 1D UV equation
!> - Charge-moment correction adapted from moist_dev/pysofft (44874e7e),
!>   model/rism/1d/equation.f90
module moist_model_moz_1d_kernel
   use mctc_env, only: wp
   use moist_model_moz_response_kernel, only: uv_response_rhs, uv_response_sr, uv_response_lr
   use moist_math_packed_sym, only: packed_index_lut
   use moist_math_lapack_getrf, only: wrap_getrf
   use moist_math_lapack_getrs, only: wrap_getrs
   implicit none(type, external)
   private

   public :: uv_kblock, uv_block_product, uv_block_product_transpose, uv_charge_moment_correction
   public :: uv_block_intramolecular_adjoint
   public :: uv_response_rhs, uv_response_sr, uv_response_lr

   !> Grid points per block: the SIMD dimension of the block products
   integer, parameter :: uv_kblock = 32

contains

   !> h = W c P for one block of k lanes
   !>
   !> h(:, (i-1) nv + a) = sum_b [sum_j W_ij c(:, (j-1) nv + b)] P_ba
   !>
   !> @param[in]  lut  packed-index lookup for nu solute sites
   !> @param[in]  nv   solvent sites
   !> @param[in]  w    solute W(k), packed (uv_kblock, nu(nu+1)/2)
   !> @param[in]  c    c(k), (uv_kblock, nu nv)
   !> @param[in]  p    solvent response P(k), (uv_kblock, nv, nv)
   !> @param[out] h    h(k), (uv_kblock, nu nv)
   pure subroutine uv_block_product(lut, nv, w, c, p, h)
      !> Packed-index lookup for the solute
      type(packed_index_lut), intent(in) :: lut
      !> Solvent sites
      integer, intent(in) :: nv
      !> Solute intramolecular correlation
      real(wp), intent(in) :: w(:, :)
      !> Direct correlation
      real(wp), intent(in) :: c(:, :)
      !> Solvent response
      real(wp), intent(in) :: p(:, :, :)
      !> Total correlation
      real(wp), intent(out) :: h(:, :)
      real(wp) :: wc(size(c, 1), nv)
      integer :: i, j, a, b

      do i = 1, lut%ns
         wc = 0.0_wp
         do j = 1, lut%ns
            do b = 1, nv
               wc(:, b) = wc(:, b) + w(:, lut%idx(i, j))*c(:, (j - 1)*nv + b)
            end do
         end do
         do a = 1, nv
            h(:, (i - 1)*nv + a) = 0.0_wp
            do b = 1, nv
               h(:, (i - 1)*nv + a) = h(:, (i - 1)*nv + a) + wc(:, b)*p(:, b, a)
            end do
         end do
      end do
   end subroutine uv_block_product

   !> Transpose of uv_block_product in the Euclidean product: v = W^T u P^T
   !>
   !> - v(:, (j-1) nv + b) = sum_i W_ij sum_a u(:, (i-1) nv + a) P_ba
   !> - W symmetric
   !>
   !> @param[in]  lut  packed-index lookup for nu solute sites
   !> @param[in]  nv   solvent sites
   !> @param[in]  w    solute W(k), packed (uv_kblock, nu(nu+1)/2)
   !> @param[in]  u    direction, (uv_kblock, nu nv)
   !> @param[in]  p    solvent response P(k), (uv_kblock, nv, nv)
   !> @param[out] v    transpose product, (uv_kblock, nu nv)
   pure subroutine uv_block_product_transpose(lut, nv, w, u, p, v)
      !> Packed-index lookup for the solute
      type(packed_index_lut), intent(in) :: lut
      !> Solvent sites
      integer, intent(in) :: nv
      !> Solute intramolecular correlation
      real(wp), intent(in) :: w(:, :)
      !> Direction
      real(wp), intent(in) :: u(:, :)
      !> Solvent response
      real(wp), intent(in) :: p(:, :, :)
      !> Transpose product
      real(wp), intent(out) :: v(:, :)
      real(wp) :: up(size(u, 1), nv)
      integer :: i, j, a, b

      v = 0.0_wp
      do i = 1, lut%ns
         ! (u P^T) of solute row i
         up = 0.0_wp
         do b = 1, nv
            do a = 1, nv
               up(:, b) = up(:, b) + u(:, (i - 1)*nv + a)*p(:, b, a)
            end do
         end do
         do j = 1, lut%ns
            do b = 1, nv
               v(:, (j - 1)*nv + b) = v(:, (j - 1)*nv + b) + w(:, lut%idx(i, j))*up(:, b)
            end do
         end do
      end do
   end subroutine uv_block_product_transpose

   !> Pull back h = W c P into the packed symmetric solute W
   !>
   !> - Off-diagonal packed entries receive both W_ij and W_ji contributions
   !>
   !> @param[in]  lut  packed-index lookup for nu solute sites
   !> @param[in]  nv   solvent sites
   !> @param[in]  c    physical c(k), (kb, nu nv)
   !> @param[in]  p    fixed solvent response, (kb, nv, nv)
   !> @param[in]  u    adjoint of h(k), (kb, nu nv)
   !> @param[out] ww   adjoint of packed W(k), (kb, nu(nu+1)/2)
   pure subroutine uv_block_intramolecular_adjoint(lut, nv, c, p, u, ww)
      !> Packed-index lookup for nu solute sites
      type(packed_index_lut), intent(in) :: lut
      !> Solvent sites
      integer, intent(in) :: nv
      !> Physical c(k), solvent response P, adjoint of h(k)
      real(wp), intent(in) :: c(:, :), p(:, :, :), u(:, :)
      !> Adjoint of packed W(k), (kb, nu(nu+1)/2)
      real(wp), intent(out) :: ww(:, :)
      real(wp) :: up(size(u, 1), nv)
      integer :: i, j, a, b

      ww = 0.0_wp
      do i = 1, lut%ns
         up = 0.0_wp
         do b = 1, nv
            do a = 1, nv
               up(:, b) = up(:, b) + u(:, (i - 1)*nv + a)*p(:, b, a)
            end do
         end do
         do j = 1, lut%ns
            do b = 1, nv
               ww(:, lut%idx(i, j)) = ww(:, lut%idx(i, j)) + up(:, b)*c(:, (j - 1)*nv + b)
            end do
         end do
      end do
   end subroutine uv_block_intramolecular_adjoint

   !> Restore the k -> 0 charge moment of the solvent response
   !>
   !> - Charge moment S_a(k) = sum_b q_b P_ba(k) of a neutral solvent vanishes as k^2
   !>   - Assembled by cancellation from entries of order one
   !>   - A finite VV residual corrupts it at the lowest modes
   !>   - A net-charged solute multiplies it by the 4 pi/k^2 pole of u_lr
   !> - S_a(k)/k^2 fitted by a quadratic in k^2 on [k_fit_lo, k_fit_hi]
   !>   - Window widened up to max_k_fit_hi to reach five points
   !>   - Imposed on modes below k_correct, and on the lowest mode below the window
   !> - Fitted and analytical moments use the minimum symmetric correction
   !>   Delta P = (q delta^T + delta q^T)/Q - (q^T delta) q q^T/Q^2
   !>   - Q = q^T q
   !>   - Restores q^T P, preserves symmetry and the charge-neutral block
   !> - No change for an uncharged solvent or a grid without a usable window
   !>
   !> @param[in]    k        wave numbers, ascending, 1/bohr (npts)
   !> @param[in]    charge   solvent site charges, e (nv)
   !> @param[inout] p        solvent response, (npts, nv, nv)
   !> @param[out]   applied  whether the correction was applied
   !> @param[in]    charge_response  optional analytical P q/k^2, (npts, nv); replaces the fit
   subroutine uv_charge_moment_correction(k, charge, p, applied, charge_response)
      !> Wave numbers
      real(wp), intent(in) :: k(:)
      !> Solvent site charges
      real(wp), intent(in) :: charge(:)
      !> Solvent response
      real(wp), intent(inout) :: p(:, :, :)
      !> Whether the correction was applied
      logical, intent(out) :: applied
      !> Optional analytical P q/k^2, (npts, nv)
      real(wp), intent(in), optional :: charge_response(:, :)
      !> Modes below this wave number are corrected, 1/bohr
      real(wp), parameter :: k_correct = 0.05_wp
      !> Fit window, 1/bohr
      real(wp), parameter :: k_fit_lo = 0.05_wp, k_fit_hi = 0.25_wp
      !> Widest window for coarse grids, 1/bohr, and fewest fit points
      real(wp), parameter :: max_k_fit_hi = 0.5_wp
      integer, parameter :: min_fit_points = 5
      real(wp) :: qq, xs, xref, pred, amat(3, 3), xpow(3)
      real(wp), allocatable :: coeff(:, :), delta(:), moment(:)
      integer :: nv, npts, nlow, jlo, jhi, j, a, m, n, info, ipiv(3)

      applied = .false.
      nv = size(charge)
      npts = size(k)
      qq = sum(charge**2)
      if (qq <= 0.0_wp .or. nv < 1 .or. npts < 1) return

      if (present(charge_response)) then
         if (any(shape(charge_response) /= [npts, nv])) return
         allocate (delta(nv))
         do j = 1, npts
            if (k(j) < 0.0_wp .or. k(j) >= k_correct) cycle
            do a = 1, nv
               delta(a) = k(j)**2*charge_response(j, a) - dot_product(charge, p(j, :, a))
            end do
            call apply_symmetric_charge_delta(charge, delta, qq, p(j, :, :))
            applied = .true.
         end do
         return
      end if

      nlow = 0
      jlo = 0
      jhi = 0
      do j = 1, npts
         if (k(j) > 0.0_wp .and. k(j) < k_correct) nlow = j
         if (jlo == 0 .and. k(j) >= k_fit_lo) jlo = j
         if (k(j) <= k_fit_hi) jhi = j
      end do
      if (jlo < 1) return
      ! A grid with no mode below k_correct still has its lowest mode below the window
      if (nlow < 1 .and. k(1) > 0.0_wp .and. k(1) < k_fit_lo) nlow = 1
      if (jhi - jlo + 1 < min_fit_points) then
         jhi = min(npts, jlo + min_fit_points - 1)
         if (jhi - jlo + 1 < min_fit_points) return
         if (k(jhi) > max_k_fit_hi) return
      end if
      if (nlow < 1 .or. jlo <= nlow) return

      ! Least squares for S_a/k^2 = c1 + c2 x + c3 x^2, x = (k/k_fit_hi)^2
      xref = k_fit_hi**2
      allocate (coeff(3, nv), delta(nv), moment(nv))
      amat = 0.0_wp
      coeff = 0.0_wp
      do j = jlo, jhi
         xs = k(j)**2/xref
         xpow = [1.0_wp, xs, xs*xs]
         do a = 1, nv
            moment(a) = dot_product(charge, p(j, :, a))/k(j)**2
         end do
         do m = 1, 3
            do n = 1, 3
               amat(m, n) = amat(m, n) + xpow(m)*xpow(n)
            end do
            coeff(m, :) = coeff(m, :) + xpow(m)*moment
         end do
      end do
      call wrap_getrf(amat, ipiv, info)
      if (info /= 0) return
      call wrap_getrs(amat, coeff, ipiv, info)
      if (info /= 0) return

      do j = 1, nlow
         xs = k(j)**2/xref
         do a = 1, nv
            pred = coeff(1, a) + coeff(2, a)*xs + coeff(3, a)*xs*xs
            delta(a) = pred*k(j)**2 - dot_product(charge, p(j, :, a))
         end do
         call apply_symmetric_charge_delta(charge, delta, qq, p(j, :, :))
      end do
      applied = .true.
   end subroutine uv_charge_moment_correction

   !> Minimum Frobenius-norm symmetric update with q^T Delta P = delta^T
   !>
   !> - For neutral vectors u and v, u^T Delta P v = 0
   !>
   !> @param[in]    charge  solvent site charges, e (nv)
   !> @param[in]    delta   target charge-moment change delta, (nv)
   !> @param[in]    qq      Q = q^T q
   !> @param[inout] p       solvent response block at one k, (nv, nv)
   pure subroutine apply_symmetric_charge_delta(charge, delta, qq, p)
      !> Solvent site charges, e (nv); target change delta (nv); Q = q^T q
      real(wp), intent(in) :: charge(:), delta(:), qq
      !> Solvent response block at one k, (nv, nv)
      real(wp), intent(inout) :: p(:, :)
      real(wp) :: qdelta
      integer :: a, b

      qdelta = dot_product(charge, delta)/qq**2
      do a = 1, size(charge)
         do b = 1, size(charge)
            p(b, a) = p(b, a) + (charge(b)*delta(a) + delta(b)*charge(a))/qq - &
               & qdelta*charge(b)*charge(a)
         end do
      end do
   end subroutine apply_symmetric_charge_delta

end module moist_model_moz_1d_kernel
