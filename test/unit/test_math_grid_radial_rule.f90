!> Test suite for the reference quadrature rules on [-1, 1]
!>
!>   - Chebyshev-II: exact sqrt(1-x^2)*p(x) moments through deg(p) = 2n-1,
!>     closed-form weight sum, weights evaluated at the stored nodes
!>   - Midpoint: exact through degree 1, closed-form x^2 deficit, equal
!>     weights 2/n, uniform node spacing 2/n
!>   - Gauss-Legendre: exact through degree 2n-1, closed-form degree-2n error,
!>     literal low-order nodes
!>   - Affine interval scaling for every rule
module test_math_grid_radial_rule
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_type, &
      & moist_math_grid_radial_rule_chebyshev2_type, new_chebyshev2_rule, &
      & moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule, &
      & moist_math_grid_radial_rule_gauss_legendre_type, new_gauss_legendre_rule
   implicit none(type, external)
   private

   public :: collect_math_grid_radial_rule

   !> Rule selectors for `make_rule`
   integer, parameter :: rule_chebyshev2 = 1, rule_midpoint = 2, rule_gauss_legendre = 3

contains

   !> Collect all math_grid_radial_rule tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_radial_rule(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("chebyshev2_weighted_moments", test_chebyshev2_weighted_moments), &
         new_unittest("chebyshev2_weight_sum", test_chebyshev2_weight_sum), &
         new_unittest("midpoint_moments", test_midpoint_moments), &
         new_unittest("gauss_legendre_exactness", test_gauss_legendre_exactness), &
         new_unittest("gauss_legendre_degree_2n_error", test_gauss_legendre_degree_2n), &
         new_unittest("gauss_legendre_literal_nodes", test_gauss_legendre_literal), &
         new_unittest("scaled_interval", test_scaled_interval) &
         ]
   end subroutine collect_math_grid_radial_rule

   !> Build a rule of the requested kind behind the abstract type
   !>
   !> @param[in]  kind  rule selector
   !> @param[out] rule  allocated rule
   subroutine make_rule(kind, rule)
      !> Rule selector
      integer, intent(in) :: kind
      !> Allocated rule
      class(moist_math_grid_radial_rule_type), allocatable, intent(out) :: rule

      type(moist_math_grid_radial_rule_chebyshev2_type) :: cheb
      type(moist_math_grid_radial_rule_midpoint_type) :: mid
      type(moist_math_grid_radial_rule_gauss_legendre_type) :: gl

      select case (kind)
      case (rule_chebyshev2)
         call new_chebyshev2_rule(cheb)
         allocate (rule, source=cheb)
      case (rule_midpoint)
         call new_midpoint_rule(mid)
         allocate (rule, source=mid)
      case default
         call new_gauss_legendre_rule(gl)
         allocate (rule, source=gl)
      end select
   end subroutine make_rule

   !> Exact integral of x^j over [-1, 1]
   !>
   !> @param[in] j  monomial degree (j >= 0)
   pure function monomial_moment(j) result(m)
      !> Monomial degree
      integer, intent(in) :: j
      !> Integral of x^j over [-1, 1]
      real(wp) :: m

      if (mod(j, 2) == 1) then
         m = 0.0_wp
      else
         m = 2.0_wp/real(j + 1, wp)
      end if
   end function monomial_moment

   !> Chebyshev-II integrates sqrt(1-x^2)*x^j exactly for j <= 2n-1
   !>
   !> Reference B((j+1)/2, 3/2) = Gamma((j+1)/2)*Gamma(3/2)/Gamma(j/2 + 2) for even j
   subroutine test_chebyshev2_weighted_moments(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(8) = [1, 2, 3, 5, 8, 16, 50, 99]
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), w(:)
      real(wp) :: s, ref
      integer :: in, n, j

      call make_rule(rule_chebyshev2, rule)
      do in = 1, size(sizes)
         n = sizes(in)
         call rule%generate(n, x, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do j = 0, 2*n - 1
            s = sum(w*sqrt(1.0_wp - x*x)*x**j)
            if (mod(j, 2) == 1) then
               ref = 0.0_wp
            else
               ref = gamma(0.5_wp*real(j + 1, wp))*gamma(1.5_wp)/gamma(0.5_wp*real(j, wp) + 2.0_wp)
            end if
            call check(error, abs(s - ref) <= 1.0e-14_wp, &
               & "Chebyshev-II: weighted moment not exact through degree 2n-1")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_chebyshev2_weighted_moments

   !> Chebyshev-II dx weights sum to pi/(n+1)*cot(pi/(2(n+1)))
   !>
   !> sum_{i=1..n} sin(i*pi/(n+1)) = cot(pi/(2(n+1))); the sum tends to 2 as
   !> 2 - pi^2/(6(n+1)^2)
   !>
   !> - Weight pi/(n+1)*sqrt((1-x_i)*(1+x_i)) at stored rounded x_i, within 4 eps
   !> - Weights from sin(i*pi/(n+1)): deviations of 30 to 3000 eps for n = 50..400
   subroutine test_chebyshev2_weight_sum(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(6) = [1, 2, 7, 50, 99, 400]
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), w(:)
      real(wp) :: ref, h
      integer :: in, n

      call make_rule(rule_chebyshev2, rule)
      do in = 1, size(sizes)
         n = sizes(in)
         call rule%generate(n, x, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         h = pi/real(n + 1, wp)
         call check(error, all(abs(w/h/sqrt((1.0_wp - x)*(1.0_wp + x)) - 1.0_wp) &
            & <= 4.0_wp*epsilon(1.0_wp)), "Chebyshev-II: dx weights must use the stored nodes")
         if (allocated(error)) return
         ref = h/tan(0.5_wp*h)
         call check(error, abs(sum(w) - ref) <= 1.0e-14_wp*ref, &
            & "Chebyshev-II: weight sum differs from pi/(n+1)*cot(pi/(2(n+1)))")
         if (allocated(error)) return
      end do
   end subroutine test_chebyshev2_weight_sum

   !> Midpoint rule is exact through degree 1; the x^2 sum is 2/3 - 2/(3n^2)
   !>
   !> Weights each equal 2/n; nodes start at -1 + 1/n with spacing 2/n
   subroutine test_midpoint_moments(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(5) = [1, 2, 7, 64, 1000]
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), w(:)
      real(wp) :: ref2
      integer :: in, n

      call make_rule(rule_midpoint, rule)
      do in = 1, size(sizes)
         n = sizes(in)
         call rule%generate(n, x, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(abs(w - 2.0_wp/real(n, wp)) <= 4.0_wp*epsilon(1.0_wp)), &
            & "Midpoint: weights must each equal 2/n")
         if (allocated(error)) return
         call check(error, abs(x(1) - (-1.0_wp + 1.0_wp/real(n, wp))) <= 4.0_wp*epsilon(1.0_wp), &
            & "Midpoint: first node must be half a cell above -1")
         if (allocated(error)) return
         if (n > 1) then
            call check(error, all(abs(x(2:) - x(:n - 1) - 2.0_wp/real(n, wp)) &
               & <= 4.0_wp*epsilon(1.0_wp)), "Midpoint: node spacing must be 2/n")
            if (allocated(error)) return
         end if
         call check(error, abs(sum(w) - 2.0_wp) <= 1.0e-14_wp, "Midpoint: weights must sum to 2")
         if (allocated(error)) return
         call check(error, abs(sum(w*x)) <= 1.0e-14_wp, "Midpoint: first moment must vanish")
         if (allocated(error)) return
         ref2 = 2.0_wp/3.0_wp - 2.0_wp/(3.0_wp*real(n, wp)**2)
         call check(error, abs(sum(w*x*x) - ref2) <= 1.0e-14_wp, &
            & "Midpoint: x^2 sum differs from 2/3 - 2/(3n^2)")
         if (allocated(error)) return
      end do
   end subroutine test_midpoint_moments

   !> Gauss-Legendre integrates x^j exactly for j <= 2n-1
   subroutine test_gauss_legendre_exactness(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(17) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 20, 32, 64, 100]
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), w(:)
      integer :: in, n, j

      call make_rule(rule_gauss_legendre, rule)
      do in = 1, size(sizes)
         n = sizes(in)
         call rule%generate(n, x, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do j = 0, 2*n - 1
            call check(error, abs(sum(w*x**j) - monomial_moment(j)) <= 2.0e-14_wp, &
               & "Gauss-Legendre: monomial moment not exact through degree 2n-1")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_gauss_legendre_exactness

   !> Gauss-Legendre degree-2n error equals the closed form
   !>
   !> For f = x^(2n): integral - sum = 2^(2n+1)*(n!)^4/((2n+1)*((2n)!)^2), so
   !> the rule is exact through 2n-1 and no further
   subroutine test_gauss_legendre_degree_2n(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), w(:)
      real(wp) :: en, s
      integer :: n

      call make_rule(rule_gauss_legendre, rule)
      do n = 1, 10
         call rule%generate(n, x, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         en = exp(real(2*n + 1, wp)*log(2.0_wp) + 4.0_wp*log_gamma(real(n + 1, wp)) &
            & - log(real(2*n + 1, wp)) - 2.0_wp*log_gamma(real(2*n + 1, wp)))
         s = sum(w*x**(2*n))
         call check(error, abs((monomial_moment(2*n) - s) - en) <= 1.0e-14_wp, &
            & "Gauss-Legendre: degree-2n error differs from the closed form")
         if (allocated(error)) return
      end do
   end subroutine test_gauss_legendre_degree_2n

   !> Gauss-Legendre literal nodes and weights for n = 1, 2, 3
   subroutine test_gauss_legendre_literal(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), w(:)
      real(wp), parameter :: tol = 4.0_wp*epsilon(1.0_wp)

      call make_rule(rule_gauss_legendre, rule)

      call rule%generate(1, x, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, x(1) == 0.0_wp .and. abs(w(1) - 2.0_wp) <= tol, &
         & "Gauss-Legendre n=1: expected x = 0, w = 2")
      if (allocated(error)) return

      call rule%generate(2, x, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(abs(x - [-1.0_wp, 1.0_wp]/sqrt(3.0_wp)) <= tol) &
         & .and. all(abs(w - 1.0_wp) <= tol), "Gauss-Legendre n=2: expected x = +-1/sqrt(3), w = 1")
      if (allocated(error)) return

      call rule%generate(3, x, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(abs(x - [-sqrt(0.6_wp), 0.0_wp, sqrt(0.6_wp)]) <= tol) &
         & .and. all(abs(w - [5.0_wp, 8.0_wp, 5.0_wp]/9.0_wp) <= tol), &
         & "Gauss-Legendre n=3: expected x = 0, +-sqrt(3/5), w = 8/9, 5/9")
      if (allocated(error)) return
      call check(error, x(2) == 0.0_wp, "Gauss-Legendre n=3: center node must be exactly 0")
   end subroutine test_gauss_legendre_literal

   !> Every rule maps affinely onto [a, b]; Gauss-Legendre stays exact there
   subroutine test_scaled_interval(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: a = 0.5_wp, b = 3.0_wp
      integer, parameter :: n = 6
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x0(:), w0(:), x(:), w(:)
      real(wp) :: ref, tol
      integer :: kind, j

      tol = 4.0_wp*epsilon(1.0_wp)
      do kind = rule_chebyshev2, rule_gauss_legendre
         call make_rule(kind, rule)
         call rule%generate(n, x0, w0, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call rule%generate(n, x, w, merr, lower=a, upper=b)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(abs(x - (0.5_wp*(a + b) + 0.5_wp*(b - a)*x0)) <= tol*b) &
            & .and. all(abs(w - 0.5_wp*(b - a)*w0) <= tol*abs(w)), &
            & "Scaled rule: nodes or weights are not the affine image of the reference rule")
         if (allocated(error)) return
         call check(error, abs(sum(w) - 0.5_wp*(b - a)*sum(w0)) <= tol*sum(w), &
            & "Scaled rule: weight sum must scale by (b - a)/2")
         if (allocated(error)) return
      end do

      call make_rule(rule_gauss_legendre, rule)
      call rule%generate(n, x, w, merr, lower=a, upper=b)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do j = 0, 2*n - 1
         ref = (b**(j + 1) - a**(j + 1))/real(j + 1, wp)
         call check(error, abs(sum(w*x**j) - ref) <= 1.0e-14_wp*ref, &
            & "Scaled Gauss-Legendre: not exact through degree 2n-1 on [a, b]")
         if (allocated(error)) return
      end do

      call make_rule(rule_midpoint, rule)
      call rule%generate(n, x, w, merr, lower=a, upper=b)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      ref = 0.5_wp*(b*b - a*a)
      call check(error, abs(sum(w*x) - ref) <= 1.0e-14_wp*ref, &
         & "Scaled midpoint: not exact for degree 1 on [a, b]")
   end subroutine test_scaled_interval

end module test_math_grid_radial_rule
