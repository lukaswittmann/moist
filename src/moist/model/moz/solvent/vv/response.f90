!> Transfer a bulk RISM solution through its short-range direct correlation
!>
!> Interpolate only cbar(k), then solve (I - rho W C) chi = W, with the
!> Coulomb part of C restored analytically. For a rigid molecule,
!> W = 11^T + k^2 D. Forming D directly avoids cancellation in W q / k^2.
module moist_model_moz_solvent_vv_response
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo
   use moist_math_packed_sym, only: packed_index_lut, npair_from_ns
   use moist_math_lapack_getrf, only: wrap_getrf
   use moist_math_lapack_getrs, only: wrap_getrs
   use moist_model_moz_solvent_vv_kernel, only: vv_transform_at
   use moist_model_moz_solvent_theory, only: bulk_theory_type, rism_theory_type
   use moist_model_moz_solvent_theory_type, only: bulk_charges_neutral, check_bulk_charge_neutrality
   implicit none(type, external)
   private
   public :: vv_response_type, new_vv_response, vv_response_point

   !> Native transform and cubic spline of the bulk short-range cbar
   type :: vv_response_type
      private
      !> Owned snapshot of the theory used by the accepted bulk solve
      class(bulk_theory_type), allocatable :: theory
      type(packed_index_lut) :: lut
      real(wp), allocatable :: nodes(:), cbar(:, :), second(:, :), distance2(:, :), q(:)
      real(wp) :: rho = 0.0_wp, beta = 0.0_wp, alpha = 1.0_wp
   contains
      procedure :: at => response_at
      !> Long-wavelength screening 1 - 1/epsilon of the bulk OZ response
      procedure :: charge_screening => response_charge_screening
   end type vv_response_type

contains

   !> Transform cbar on its native pair and prepare a clamped cubic spline
   !>
   !> Include k = 0 and the band edge by radial quadrature. The derivative at
   !> zero is exactly zero; the derivative at the band edge is also evaluated
   !> from the radial table, so sampling never clamps or extrapolates cbar.
   !> @param[in] theory Optional fixed bulk theory; ordinary RISM when absent
   subroutine new_vv_response(self, rgrid, kgrid, cr_short, xyz, q, rho, beta, alpha, nthreads, error, theory)
      type(vv_response_type), intent(out) :: self
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      real(wp), intent(in) :: cr_short(:, :), xyz(:, :), q(:), rho, beta, alpha
      integer, intent(in) :: nthreads
      type(error_type), allocatable, intent(out) :: error
      class(bulk_theory_type), intent(in), optional :: theory
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      real(wp), allocatable :: ends(:, :), derivative(:), kernel(:), x(:)
      real(wp) :: band
      integer :: n, ns, i, j

      n = kgrid%npts
      ns = size(q)
      if (n < 1 .or. size(xyz, 2) /= ns .or. &
         & any(shape(cr_short) /= [rgrid%npts, npair_from_ns(ns)])) then
         call fatal_error(error, "Bulk OZ response needs a native radial pair and packed direct correlations")
         return
      end if
      call check_bulk_charge_neutrality(q, error)
      if (allocated(error)) return
      call self%lut%init(ns)
      self%q = q
      self%rho = rho
      self%beta = beta
      self%alpha = alpha
      if (present(theory)) then
         allocate (self%theory, source=theory)
      else
         allocate (rism_theory_type :: self%theory)
      end if
      allocate (self%distance2(ns, ns))
      do j = 1, ns
         do i = 1, ns
            self%distance2(i, j) = sum((xyz(:, i) - xyz(:, j))**2)
         end do
      end do
      band = real(n, wp)*kgrid%w(1)
      self%nodes = [0.0_wp, kgrid%r, band]
      allocate (self%cbar(n + 2, size(cr_short, 2)), ends(2, size(cr_short, 2)))
      call new_radial_trafo(trafo, rgrid, kgrid, error)
      if (.not. allocated(error)) call trafo%fbt_r2k_all(cr_short, self%cbar(2:n + 1, :), error)
      if (allocated(error)) return
      call vv_transform_at(rgrid%r, rgrid%w, cr_short, [0.0_wp, band], ends, nthreads)
      self%cbar(1, :) = ends(1, :)
      self%cbar(n + 2, :) = ends(2, :)
      x = band*rgrid%r
      kernel = (4.0_wp*acos(-1.0_wp)/band)*rgrid%w*rgrid%r**2*(cos(x) - sin(x)/x)
      derivative = matmul(kernel, cr_short)
      call spline_second(self%nodes, self%cbar, derivative, self%second)
   end subroutine new_vv_response

   !> Second derivatives of a cubic spline, clamped to f'(0) = 0 and f'(band)
   subroutine spline_second(x, f, derivative, second)
      real(wp), intent(in) :: x(:), f(:, :), derivative(:)
      real(wp), allocatable, intent(out) :: second(:, :)
      real(wp), allocatable :: h(:), diagonal(:), upper(:)
      real(wp) :: factor
      integer :: n, i

      n = size(x)
      h = x(2:) - x(:n - 1)
      allocate (diagonal(n), upper(n), second(n, size(f, 2)))
      diagonal(1) = 2.0_wp*h(1)
      upper(1) = h(1)
      second(1, :) = 6.0_wp*(f(2, :) - f(1, :))/h(1)
      do i = 2, n - 1
         diagonal(i) = 2.0_wp*(h(i - 1) + h(i))
         upper(i) = h(i)
         second(i, :) = 6.0_wp*((f(i + 1, :) - f(i, :))/h(i) - (f(i, :) - f(i - 1, :))/h(i - 1))
      end do
      diagonal(n) = 2.0_wp*h(n - 1)
      upper(n) = 0.0_wp
      second(n, :) = 6.0_wp*(derivative - (f(n, :) - f(n - 1, :))/h(n - 1))
      do i = 2, n
         factor = h(i - 1)/diagonal(i - 1)
         diagonal(i) = diagonal(i) - factor*upper(i - 1)
         second(i, :) = second(i, :) - factor*second(i - 1, :)
      end do
      second(n, :) = second(n, :)/diagonal(n)
      do i = n - 1, 1, -1
         second(i, :) = (second(i, :) - upper(i)*second(i + 1, :))/diagonal(i)
      end do
   end subroutine spline_second

   !> Interpolate cbar and reconstruct chi at each requested wave number
   !>
   !> charge_response is chi q / k^2, computed without subtracting nearly
   !> equal susceptibility entries; at zero it holds the analytical limit.
   subroutine response_at(self, k, chi, error, charge_response)
      class(vv_response_type), intent(in) :: self
      real(wp), intent(in) :: k(:)
      real(wp), allocatable, intent(out) :: chi(:, :)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable, intent(out), optional :: charge_response(:, :)
      real(wp), allocatable :: short(:), c(:, :), p(:, :), z(:)
      real(wp) :: h, a, b
      integer :: ik, lo, hi, mid, ns, i, j, info, n
      character(len=48) :: buf

      if (.not. allocated(self%nodes)) then
         call fatal_error(error, "Bulk OZ response has no current direct correlation transform")
         return
      end if
      if (.not. all(ieee_is_finite(k))) then
         call fatal_error(error, "Bulk OZ response wave numbers must be finite")
         return
      end if
      n = size(self%nodes)
      if (any(k < 0.0_wp) .or. any(k > self%nodes(n))) then
         call fatal_error(error, "Bulk OZ response: a wave number lies outside the VV grid's band")
         return
      end if
      ns = self%lut%ns
      allocate (chi(size(k), size(self%cbar, 2)), short(size(self%cbar, 2)), c(ns, ns), p(ns, ns), z(ns))
      if (present(charge_response)) allocate (charge_response(size(k), ns))
      do ik = 1, size(k)
         lo = 1
         hi = n
         do while (hi - lo > 1)
            mid = (hi + lo)/2
            if (k(ik) < self%nodes(mid)) then
               hi = mid
            else
               lo = mid
            end if
         end do
         h = self%nodes(hi) - self%nodes(lo)
         a = (self%nodes(hi) - k(ik))/h
         b = (k(ik) - self%nodes(lo))/h
         short = a*self%cbar(lo, :) + b*self%cbar(hi, :) + &
            & ((a**3 - a)*self%second(lo, :) + (b**3 - b)*self%second(hi, :))*h*h/6.0_wp
         do j = 1, ns
            do i = 1, ns
               c(i, j) = short(self%lut%idx(i, j))
            end do
         end do
         call vv_response_point(k(ik), c, self%distance2, self%q, self%rho, self%beta, self%alpha, p, z, info, &
            & self%theory)
         if (info /= 0) then
            write (buf, "(es16.8)") k(ik)
            call fatal_error(error, "Bulk OZ susceptibility reconstruction failed at k = "//trim(adjustl(buf)))
            return
         end if
         do j = 1, ns
            do i = 1, j
               chi(ik, self%lut%idx(i, j)) = p(i, j)
            end do
         end do
         if (present(charge_response)) charge_response(ik, :) = z
      end do
   end subroutine response_at

   !> chi here omits the molecular density: S = 4 pi beta rho lim(q^T chi q/k^2)
   !> Use the regularized charge limit; contracting chi(0) would give zero.
   subroutine response_charge_screening(self, screening, error)
      class(vv_response_type), intent(in) :: self
      real(wp), intent(out) :: screening
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: chi(:, :), z(:, :)

      screening = 0.0_wp
      call self%at([0.0_wp], chi, error, z)
      if (allocated(error)) return
      screening = 4.0_wp*acos(-1.0_wp)*self%beta*self%rho*dot_product(self%q, z(1, :))
      if (.not. ieee_is_finite(screening) .or. screening < -64.0_wp*epsilon(1.0_wp) .or. &
         & screening > 1.0_wp + 64.0_wp*epsilon(1.0_wp)) then
         call fatal_error(error, "Bulk OZ response does not define dielectric screening in [0,1]")
         return
      end if
      ! A vanishing dipole can leave roundoff in the charge projection.
      screening = max(0.0_wp, min(1.0_wp, screening))
   end subroutine response_charge_screening

   !> Stable molecular OZ solve with analytical form factors and Coulomb tail
   !>
   !> C = cbar - gamma q q^T / k^2, gamma = 4 pi beta exp(-k^2/(4 alpha^2)).
   !> W = omega + rho zeta and D = (omega - 11^T)/k^2 + rho zeta/k^2.
   !> For neutral bulk charges, W q = k^2 D q, so A = I - rho W cbar + rho gamma (D q) q^T
   !> stays finite at zero. Solve A chi = W and A z = D q together.
   !> At zero chi = 11^T / (1 - rho sum(cbar(0))); z is its charge limit.
   !> @param[out] info Zero on success; -4 for invalid bulk charges; otherwise the solve status
   !> @param[in] theory Optional fixed bulk theory; ordinary RISM when absent
   subroutine vv_response_point(k, cbar, distance2, q, rho, beta, alpha, chi, z, info, theory)
      real(wp), intent(in) :: k, cbar(:, :), distance2(:, :), q(:), rho, beta, alpha
      real(wp), intent(out) :: chi(:, :), z(:)
      integer, intent(out) :: info
      class(bulk_theory_type), intent(in), optional :: theory
      real(wp) :: w(size(q), size(q)), d(size(q), size(q)), a(size(q), size(q))
      real(wp) :: rhs(size(q), size(q) + 1), dq(size(q)), x2, x, k2, gamma
      real(wp) :: zeta(size(q), size(q)), zeta_over_k2(size(q), size(q))
      integer :: ns, i, j, ipiv(size(q))

      if (.not. bulk_charges_neutral(q)) then
         chi = 0.0_wp
         z = 0.0_wp
         info = -4
         return
      end if
      ns = size(q)
      k2 = k*k
      do j = 1, ns
         do i = 1, ns
            x2 = k2*distance2(i, j)
            if (x2 < 0.01_wp) then
               d(i, j) = distance2(i, j)*(-1.0_wp/6.0_wp + x2*(1.0_wp/120.0_wp + &
                  & x2*(-1.0_wp/5040.0_wp + x2*(1.0_wp/362880.0_wp - x2/39916800.0_wp))))
               w(i, j) = 1.0_wp + k2*d(i, j)
            else
               x = sqrt(x2)
               w(i, j) = sin(x)/x
               d(i, j) = (w(i, j) - 1.0_wp)/k2
            end if
         end do
      end do
      if (present(theory)) then
         call theory%correction_at(k, zeta, zeta_over_k2)
         w = w + rho*zeta
         d = d + rho*zeta_over_k2
      end if
      dq = matmul(d, q)
      gamma = 4.0_wp*acos(-1.0_wp)*beta*exp(-0.25_wp*k2/(alpha*alpha))
      a = -rho*matmul(w, cbar)
      do j = 1, ns
         a(:, j) = a(:, j) + rho*gamma*dq*q(j)
         a(j, j) = a(j, j) + 1.0_wp
      end do
      call wrap_getrf(a, ipiv, info)
      if (info /= 0) return
      rhs(:, :ns) = w
      rhs(:, ns + 1) = dq
      call wrap_getrs(a, rhs, ipiv, info)
      if (info /= 0) return
      chi = 0.5_wp*(rhs(:, :ns) + transpose(rhs(:, :ns)))
      if (k == 0.0_wp) chi = 1.0_wp/(1.0_wp - rho*sum(cbar))
      z = rhs(:, ns + 1)
      if (.not. all(ieee_is_finite(chi)) .or. .not. all(ieee_is_finite(z))) info = 1
   end subroutine vv_response_point

end module moist_model_moz_solvent_vv_response
