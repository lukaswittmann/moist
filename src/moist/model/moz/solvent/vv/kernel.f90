!> Per-k kernels of the 1D VV equation
!>
!> - Intramolecular correlation w(k), packed (npts, npair), at any k >= 0
!> - Radial transform of r-space tables at arbitrary wave numbers
!> - Lane-major block layout (vv_kblock, ns, ns, nblock) of w(k)
!> - OZ solve h = (I - w c R)^-1 w c w per block or per point, R = diag(rho)
!>
!> The OZ solve is adapted from moist_dev/pysofft (44874e7e), model/rism/1d/equation.f90,
!> and keeps its association of the products
module moist_model_moz_solvent_vv_kernel
   use mctc_env, only: wp
   use moist_context, only: resolve_num_threads
   use moist_math_packed_sym, only: packed_index_lut
   use moist_math_lapack, only: getrf, getrs
   implicit none(type, external)
   private

   public :: vv_kblock, vv_intramolecular, vv_transform_at, vv_build_blocks, vv_solve_block, vv_solve_point
   public :: vv_intramolecular_adjoint

   !> Grid points per block: the SIMD dimension of the blocked solve and products
   integer, parameter :: vv_kblock = 32

contains

   !> Intramolecular correlation w_ij(k) = sin(k d_ij)/(k d_ij) of one rigid molecule
   !>
   !> Series 1 - x^2/6 below x = 1e-4, so k = 0 gives 1
   !>
   !> @param[in]  k    wave numbers >= 0, 1/bohr (npts)
   !> @param[in]  xyz  site coordinates, bohr (3, ns)
   !> @param[in]  lut  packed-index lookup for ns sites
   !> @param[out] wk   w(k), packed (npts, npair); 1 for coincident sites
   pure subroutine vv_intramolecular(k, xyz, lut, wk)
      !> Wave numbers, 1/bohr
      real(wp), intent(in) :: k(:)
      !> Site coordinates, bohr
      real(wp), intent(in) :: xyz(:, :)
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Packed w(k)
      real(wp), intent(out) :: wk(:, :)
      real(wp) :: d, x
      integer :: i, j, ip, ik

      do j = 1, size(xyz, 2)
         do i = 1, j
            ip = lut%idx(i, j)
            d = norm2(xyz(:, i) - xyz(:, j))
            if (d == 0.0_wp) then
               wk(:, ip) = 1.0_wp
            else
               do ik = 1, size(k)
                  x = d*k(ik)
                  if (abs(x) < 1.0e-4_wp) then
                     wk(ik, ip) = 1.0_wp - x*x/6.0_wp
                  else
                     wk(ik, ip) = sin(x)/x
                  end if
               end do
            end if
         end do
      end do
   end subroutine vv_intramolecular

   !> Contract an adjoint of packed w_ij(k) into site coordinates
   !>
   !> - Small-x branch differentiates the quadratic used by vv_intramolecular
   !> - Longer series nearby avoids cancellation in the derivative of sin(x)/x
   !>
   !> @param[in]     k         wave numbers, 1/bohr (npts)
   !> @param[in]     xyz       site coordinates, bohr (3, ns)
   !> @param[in]     lut       packed-index lookup for ns sites
   !> @param[in]     ww        adjoint of packed w(k), (npts, npair)
   !> @param[in,out] gradient  coordinate derivative accumulator (3, ns)
   pure subroutine vv_intramolecular_adjoint(k, xyz, lut, ww, gradient)
      !> Wave numbers (1/bohr), site coordinates (bohr) and adjoint of packed w(k)
      real(wp), intent(in) :: k(:), xyz(:, :), ww(:, :)
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Coordinate derivative accumulator (3, ns)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp) :: delta(3), pair_gradient(3), d, x, x2, factor, seed
      integer :: i, j, ip, ik

      do j = 2, size(xyz, 2)
         do i = 1, j - 1
            ip = lut%idx(i, j)
            delta = xyz(:, i) - xyz(:, j)
            d = norm2(delta)
            if (d == 0.0_wp) cycle
            seed = 0.0_wp
            do ik = 1, size(k)
               x = d*k(ik)
               if (abs(x) < 1.0e-4_wp) then
                  factor = -1.0_wp/3.0_wp
               else if (abs(x) < 1.0e-2_wp) then
                  x2 = x*x
                  factor = -1.0_wp/3.0_wp + x2*(1.0_wp/30.0_wp + &
                     & x2*(-1.0_wp/840.0_wp + x2/45360.0_wp))
               else
                  factor = (x*cos(x) - sin(x))/(x*x*x)
               end if
               seed = seed + ww(ik, ip)*k(ik)**2*factor
            end do
            pair_gradient = seed*delta
            gradient(:, i) = gradient(:, i) + pair_gradient
            gradient(:, j) = gradient(:, j) - pair_gradient
         end do
      end do
   end subroutine vv_intramolecular_adjoint

   !> Radial transform of r-space tables at arbitrary wave numbers
   !>
   !> - f(k) = 4 pi sum_i w_i r_i^2 f(r_i) j0(k r_i)
   !> - Quadrature the DST-IV pair evaluates at its nodes, here at any k
   !> - No interpolation, f(r) is zero beyond the last node
   !>
   !> @param[in]  r         radial nodes, bohr (npts)
   !> @param[in]  w         radial weights, bohr (npts)
   !> @param[in]  fr        r-space tables (npts, ncol)
   !> @param[in]  k         wave numbers >= 0, 1/bohr (nk)
   !> @param[out] fk        transforms (nk, ncol)
   !> @param[in]  nthreads  OpenMP team size; zero takes omp_get_max_threads
   subroutine vv_transform_at(r, w, fr, k, fk, nthreads)
      !> Radial nodes and weights
      real(wp), intent(in) :: r(:), w(:)
      !> r-space tables
      real(wp), intent(in) :: fr(:, :)
      !> Wave numbers
      real(wp), intent(in) :: k(:)
      !> Transforms
      real(wp), intent(out) :: fk(:, :)
      !> OpenMP team size
      integer, intent(in) :: nthreads
      real(wp), parameter :: four_pi = 4.0_wp*acos(-1.0_wp)
      real(wp), allocatable :: kernel(:)
      real(wp) :: x
      integer :: ik, i, nt

      nt = resolve_num_threads(nthreads)
      !$omp parallel num_threads(nt) default(none) shared(r, w, fr, k, fk) private(kernel, ik, i, x)
      allocate (kernel(size(r)))
      !$omp do schedule(static)
      do ik = 1, size(k)
         do i = 1, size(r)
            x = k(ik)*r(i)
            if (abs(x) < 1.0e-4_wp) then
               kernel(i) = four_pi*w(i)*r(i)**2*(1.0_wp - x*x/6.0_wp)
            else
               kernel(i) = four_pi*w(i)*r(i)**2*sin(x)/x
            end if
         end do
         fk(ik, :) = matmul(kernel, fr)
      end do
      !$omp end do
      !$omp end parallel
   end subroutine vv_transform_at

   !> w(k) in the lane-major block layout
   !>
   !> - Grid point ik sits in lane mod(ik - 1, vv_kblock) + 1
   !> - Block index of grid point ik is (ik - 1)/vv_kblock + 1
   !> - Lanes past the end of the grid are zero
   !>
   !> @param[in]  wk     w(k), packed (npts, npair)
   !> @param[in]  lut    packed-index lookup for ns sites
   !> @param[out] w_blk  w(k), (vv_kblock, ns, ns, nblock)
   pure subroutine vv_build_blocks(wk, lut, w_blk)
      !> Packed w(k)
      real(wp), intent(in) :: wk(:, :)
      !> Packed-index lookup
      type(packed_index_lut), intent(in) :: lut
      !> Blocked w(k)
      real(wp), allocatable, intent(out) :: w_blk(:, :, :, :)
      integer :: ik, kk, ib, is, js, npts

      npts = size(wk, 1)
      allocate (w_blk(vv_kblock, lut%ns, lut%ns, (npts - 1)/vv_kblock + 1), source=0.0_wp)
      do ik = 1, npts
         kk = mod(ik - 1, vv_kblock) + 1
         ib = (ik - 1)/vv_kblock + 1
         do js = 1, lut%ns
            do is = 1, lut%ns
               w_blk(kk, is, js, ib) = wk(ik, lut%idx(is, js))
            end do
         end do
      end do
   end subroutine vv_build_blocks

   !> OZ solve of one block of vv_kblock points at once, lane-major
   !>
   !> - h = (I - w (c R))^-1 w (c w), associated as in vv_solve_point:
   !>   near k = 0 the system is ill-conditioned, and regrouping as (w c) w moved
   !>   h(k) by 1e-10 relative there, enough to stall the water VV solve
   !> - LU without row interchanges; a lane where partial pivoting would swap
   !>   (LAPACK's first-maximum rule) or meets a zero pivot is flagged, continues
   !>   on a unit pivot, and its h is meaningless: re-solve it with vv_solve_point
   !>
   !> @param[in]  ns          number of sites
   !> @param[in]  w           w(k), (vv_kblock, ns, ns)
   !> @param[in]  ck          c(k), (vv_kblock, ns, ns)
   !> @param[in]  rho         site number densities (ns)
   !> @param[out] hk          h(k), (vv_kblock, ns, ns)
   !> @param[out] need_pivot  lanes to re-solve with pivoting
   pure subroutine vv_solve_block(ns, w, ck, rho, hk, need_pivot)
      !> Number of sites
      integer, intent(in) :: ns
      !> Intramolecular correlation
      real(wp), intent(in) :: w(vv_kblock, ns, ns)
      !> Direct correlation
      real(wp), intent(in) :: ck(vv_kblock, ns, ns)
      !> Site number densities
      real(wp), intent(in) :: rho(ns)
      !> Total correlation
      real(wp), intent(out) :: hk(vv_kblock, ns, ns)
      !> Lanes whose result must be re-solved with pivoting
      logical, intent(out) :: need_pivot(vv_kblock)
      !> c w, and I - w c R with its LU factors
      real(wp) :: x(vv_kblock, ns, ns), a(vv_kblock, ns, ns)
      !> Reciprocal pivots and the largest subdiagonal excess over the pivot
      real(wp) :: inv(vv_kblock), excess(vv_kblock)
      integer :: i, j, k

      a = 0.0_wp
      x = 0.0_wp
      do j = 1, ns
         do k = 1, ns
            do i = 1, ns
               a(:, i, j) = a(:, i, j) + w(:, i, k)*(ck(:, k, j)*rho(j))
               x(:, i, j) = x(:, i, j) + ck(:, i, k)*w(:, k, j)
            end do
         end do
      end do
      a = -a
      do j = 1, ns
         a(:, j, j) = a(:, j, j) + 1.0_wp
      end do
      hk = 0.0_wp
      do j = 1, ns
         do k = 1, ns
            do i = 1, ns
               hk(:, i, j) = hk(:, i, j) + w(:, i, k)*x(:, k, j)
            end do
         end do
      end do

      ! LU without interchanges, forward elimination applied to hk alongside
      need_pivot = .false.
      do k = 1, ns
         excess = -abs(a(:, k, k))
         do i = k + 1, ns
            excess = max(excess, abs(a(:, i, k)) - abs(a(:, k, k)))
         end do
         need_pivot = need_pivot .or. excess > 0.0_wp .or. a(:, k, k) == 0.0_wp
         where (need_pivot) a(:, k, k) = 1.0_wp
         inv = 1.0_wp/a(:, k, k)
         do i = k + 1, ns
            a(:, i, k) = a(:, i, k)*inv
         end do
         do j = k + 1, ns
            do i = k + 1, ns
               a(:, i, j) = a(:, i, j) - a(:, i, k)*a(:, k, j)
            end do
         end do
         do j = 1, ns
            do i = k + 1, ns
               hk(:, i, j) = hk(:, i, j) - a(:, i, k)*hk(:, k, j)
            end do
         end do
      end do
      ! Back substitution
      do j = 1, ns
         do k = ns, 1, -1
            hk(:, k, j) = hk(:, k, j)/a(:, k, k)
            do i = 1, k - 1
               hk(:, i, j) = hk(:, i, j) - a(:, i, k)*hk(:, k, j)
            end do
         end do
      end do
   end subroutine vv_solve_block

   !> OZ solve at one point with LAPACK partial pivoting
   !>
   !> h = (I - w (c R))^-1 w (c w)
   !>
   !> @param[in]  ns     number of sites
   !> @param[in]  w      w(k), (ns, ns)
   !> @param[in]  ck     c(k), (ns, ns)
   !> @param[in]  rho    site number densities (ns)
   !> @param[out] hk     h(k), (ns, ns)
   !> @param[out] info   0 on success, LAPACK status otherwise
   subroutine vv_solve_point(ns, w, ck, rho, hk, info)
      !> Number of sites
      integer, intent(in) :: ns
      !> Intramolecular correlation
      real(wp), intent(in) :: w(ns, ns)
      !> Direct correlation
      real(wp), intent(in) :: ck(ns, ns)
      !> Site number densities
      real(wp), intent(in) :: rho(ns)
      !> Total correlation
      real(wp), intent(out) :: hk(ns, ns)
      !> LAPACK status
      integer, intent(out) :: info
      real(wp) :: ck_rho(ns, ns), temp(ns, ns), iwcp(ns, ns)
      integer :: ipiv(ns), is, js

      do js = 1, ns
         ck_rho(:, js) = ck(:, js)*rho(js)
      end do
      iwcp = -matmul(w, ck_rho)
      do is = 1, ns
         iwcp(is, is) = iwcp(is, is) + 1.0_wp
      end do
      call getrf(iwcp, ipiv, info)
      if (info /= 0) return
      temp = matmul(ck, w)
      hk = matmul(w, temp)
      call getrs(iwcp, hk, ipiv, info)
   end subroutine vv_solve_point

end module moist_model_moz_solvent_vv_kernel
