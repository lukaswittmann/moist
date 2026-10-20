!> Transfer a bulk RISM solution through its short-range direct correlation
!>
!> - Interpolate only cbar(k), then solve (I - rho W C) chi = W
!> - Coulomb part of C restored analytically
!> - Rigid molecule: W = 11^T + k^2 D
!> - Forming D directly avoids cancellation in W q / k^2
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
      !> Packed-index lookup for ns sites
      type(packed_index_lut) :: lut
      !> Spline data of cbar and site geometry
      !>
      !> - nodes: 0, k grid, band edge, 1/bohr (npts + 2)
      !> - cbar, second: values and spline second derivatives, packed (npts + 2, npair)
      !> - distance2: squared site distances, bohr^2 (ns, ns)
      !> - q: site charges, e (ns)
      real(wp), allocatable :: nodes(:), cbar(:, :), second(:, :), distance2(:, :), q(:)
      !> State point and splitting
      !>
      !> - rho: molecular number density, 1/bohr^3
      !> - beta: inverse thermal energy, 1/Hartree
      !> - alpha: Coulomb splitting parameter, 1/bohr
      real(wp) :: rho = 0.0_wp, beta = 0.0_wp, alpha = 1.0_wp
   contains
      !> Interpolate cbar and reconstruct chi at wave numbers
      procedure :: at => response_at
      !> Long-wavelength screening 1 - 1/epsilon of the bulk OZ response
      procedure :: charge_screening => response_charge_screening
   end type vv_response_type

contains

   !> Transform cbar on its native pair and prepare a clamped cubic spline
   !>
   !> - k = 0 and the band edge by radial quadrature
   !> - Derivative at zero is exactly zero
   !> - Derivative at the band edge is also evaluated from the radial table
   !> - Sampling never clamps or extrapolates cbar
   !>
   !> @param[out] self      response with spline of the bulk cbar
   !> @param[in]  rgrid     real-space radial grid of the solve
   !> @param[in]  kgrid     reciprocal grid of the pair
   !> @param[in]  cr_short  short-range direct correlation, packed (npts, npair)
   !> @param[in]  xyz       site coordinates, bohr (3, ns)
   !> @param[in]  q         site charges, e (ns)
   !> @param[in]  rho       molecular number density, 1/bohr^3
   !> @param[in]  beta      inverse thermal energy, 1/Hartree
   !> @param[in]  alpha     Coulomb splitting parameter, 1/bohr
   !> @param[in]  nthreads  OpenMP team size; zero takes omp_get_max_threads
   !> @param[out] error     mismatched shapes, non-neutral charges or transform failure
   !> @param[in]  theory    optional fixed bulk theory; ordinary RISM when absent
   subroutine new_vv_response(self, rgrid, kgrid, cr_short, xyz, q, rho, beta, alpha, nthreads, error, theory)
      !> Response to build
      type(vv_response_type), intent(out) :: self
      !> Real-space and reciprocal radial grids
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Short-range c(r) and site data, units in the @param lines
      real(wp), intent(in) :: cr_short(:, :), xyz(:, :), q(:), rho, beta, alpha
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Setup failure
      type(error_type), allocatable, intent(out) :: error
      !> Fixed bulk theory, ordinary RISM when absent
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
   !>
   !> @param[in]  x           nodes, increasing (n)
   !> @param[in]  f           node values (n, ncol)
   !> @param[in]  derivative  f' at the last node (ncol)
   !> @param[out] second      second derivatives (n, ncol)
   subroutine spline_second(x, f, derivative, second)
      !> Nodes, values and last-node derivative
      real(wp), intent(in) :: x(:), f(:, :), derivative(:)
      !> Second derivatives
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
   !> - chi omits the molecular density
   !> - charge_response is chi q / k^2, computed without subtracting nearly equal susceptibility entries
   !> - charge_response holds the analytical limit at k = 0
   !>
   !> @param[in]  self             response with spline of the bulk cbar
   !> @param[in]  k                wave numbers in [0, band edge], 1/bohr (nk)
   !> @param[out] chi              susceptibility chi(k), packed (nk, npair)
   !> @param[out] error            missing transform, nonfinite or out-of-band k, failed solve
   !> @param[out] charge_response  chi q / k^2 (nk, ns), optional
   subroutine response_at(self, k, chi, error, charge_response)
      !> Response to evaluate
      class(vv_response_type), intent(in) :: self
      !> Wave numbers, 1/bohr
      real(wp), intent(in) :: k(:)
      !> Packed susceptibility
      real(wp), allocatable, intent(out) :: chi(:, :)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error
      !> Charge response chi q / k^2
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

   !> Long-wavelength dielectric screening of the bulk OZ response
   !>
   !> - chi omits the molecular density: S = 4 pi beta rho lim(q^T chi q/k^2)
   !> - Uses the regularized charge limit; contracting chi(0) would give zero
   !>
   !> @param[in]  self       response with spline of the bulk cbar
   !> @param[out] screening  S, dimensionless, in [0, 1]
   !> @param[out] error      failed evaluation or S outside [0, 1]
   subroutine response_charge_screening(self, screening, error)
      !> Response to evaluate
      class(vv_response_type), intent(in) :: self
      !> Screening S
      real(wp), intent(out) :: screening
      !> Evaluation failure
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
      ! A vanishing dipole can leave roundoff in the charge projection
      screening = max(0.0_wp, min(1.0_wp, screening))
   end subroutine response_charge_screening

   !> Stable molecular OZ solve with analytical form factors and Coulomb tail
   !>
   !> - C = cbar - gamma q q^T / k^2
   !> - gamma = 4 pi beta exp(-k^2/(4 alpha^2))
   !> - W = omega + rho zeta, D = (omega - 11^T)/k^2 + rho zeta/k^2
   !> - Neutral bulk charges give W q = k^2 D q
   !> - A = I - rho W cbar + rho gamma (D q) q^T
   !> - A stays finite at zero
   !> - Solve A chi = W and A z = D q together
   !> - At zero chi = 11^T / (1 - rho sum(cbar(0))); z is its charge limit
   !>
   !> @param[in]  k          wave number, 1/bohr
   !> @param[in]  cbar       short-range direct correlation at k (ns, ns)
   !> @param[in]  distance2  squared site distances, bohr^2 (ns, ns)
   !> @param[in]  q          site charges, e (ns)
   !> @param[in]  rho        molecular number density, 1/bohr^3
   !> @param[in]  beta       inverse thermal energy, 1/Hartree
   !> @param[in]  alpha      Coulomb splitting parameter, 1/bohr
   !> @param[out] chi        susceptibility chi(k) (ns, ns)
   !> @param[out] z          charge limit A^-1 D q (ns)
   !> @param[out] info       zero on success; -4 for invalid bulk charges; otherwise the solve status
   !> @param[in]  theory     optional fixed bulk theory; ordinary RISM when absent
   subroutine vv_response_point(k, cbar, distance2, q, rho, beta, alpha, chi, z, info, theory)
      !> Wave number, 1/bohr
      real(wp), intent(in) :: k
      !> Direct correlation at k, squared site distances, site charges and state point
      real(wp), intent(in) :: cbar(:, :), distance2(:, :), q(:), rho, beta, alpha
      !> Susceptibility and charge limit
      real(wp), intent(out) :: chi(:, :), z(:)
      !> Status
      integer, intent(out) :: info
      !> Fixed bulk theory, ordinary RISM when absent
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
