!> Test suite for radial mappings, the radial recipe, and the concrete radial grid
!>
!>   - Mappings: analytic mapped integrals (polynomials, Gaussians,
!>     exponentials), element-dependent scales, domain and divergent-endpoint
!>     errors, forward map against transform, polymorphic copies
!>   - Radial recipe: lower and upper cutoffs with requested and retained
!>     counts, invalid recipes, volume integrals, element dependence
!>   - Uniform radial pair: nodes, weights, spacing, transform tag
module test_math_grid_radial
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, &
      & ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_data_atomicrad, only: covalent_rad
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_type, &
      & moist_math_grid_radial_rule_chebyshev2_type, new_chebyshev2_rule, &
      & moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule, &
      & moist_math_grid_radial_rule_gauss_legendre_type, new_gauss_legendre_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_type, check_mapping_input, &
      & moist_math_grid_radial_mapping_linear_type, new_linear_mapping, &
      & moist_math_grid_radial_mapping_becke_type, new_becke_mapping, &
      & moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping, &
      & moist_math_grid_radial_mapping_knowles_type, new_knowles_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, &
      & moist_math_grid_radial_recipe_type, new_radial_grid, new_uniform_radial_pair, &
      & transform_quadrature, transform_dst4
   implicit none(type, external)
   private

   public :: collect_math_grid_radial

   !> Rule selectors for `make_rule`
   integer, parameter :: rule_chebyshev2 = 1, rule_midpoint = 2, rule_gauss_legendre = 3

   !> Two ULP at unit scale
   real(wp), parameter :: two_ulp = 2.0_wp*epsilon(1.0_wp)

contains

   !> Collect all math_grid_radial tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_radial(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("mapping_validation_contracts", test_mapping_validation), &
         new_unittest("mapping_output_overflow", test_mapping_overflow), &
         new_unittest("linear_polynomial_exact", test_linear_polynomial), &
         new_unittest("linear_matches_scaled_rule", test_linear_scaled_rule), &
         new_unittest("becke_mapped_integrals", test_becke_integrals), &
         new_unittest("becke_element_scale", test_becke_element_scale), &
         new_unittest("becke_domain_errors", test_becke_errors), &
         new_unittest("handymod_finite_integrals", test_handymod_integrals), &
         new_unittest("handymod_endpoints", test_handymod_endpoints), &
         new_unittest("handymod_domain_errors", test_handymod_errors), &
         new_unittest("knowles_mapped_integrals", test_knowles_integrals), &
         new_unittest("knowles_mura_knowles_nodes", test_knowles_mura_knowles), &
         new_unittest("knowles_element_table", test_knowles_element_table), &
         new_unittest("knowles_domain_errors", test_knowles_errors), &
         new_unittest("mapping_map_matches_transform", test_map_matches_transform), &
         new_unittest("mapping_polymorphic_copy", test_mapping_copy), &
         new_unittest("recipe_cutoffs", test_recipe_cutoffs), &
         new_unittest("recipe_errors", test_recipe_errors), &
         new_unittest("recipe_volume_integrals", test_recipe_integrals), &
         new_unittest("recipe_element_dependence", test_recipe_element), &
         new_unittest("uniform_pair_contract", test_uniform_pair), &
         new_unittest("uniform_pair_errors", test_uniform_pair_errors) &
         ]
   end subroutine collect_math_grid_radial

   !> Constructor bounds, forward domains, and direct input validation
   subroutine test_mapping_validation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_linear_type) :: lin
      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(moist_math_grid_radial_mapping_knowles_type) :: kn
      class(moist_math_grid_radial_mapping_type), allocatable :: mapping
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      real(wp) :: nan, inf, rx, x, u, ref
      integer :: kind

      nan = ieee_value(0.0_wp, ieee_quiet_nan)
      inf = ieee_value(0.0_wp, ieee_positive_inf)
      call new_linear_mapping(lin, 0.0_wp, inf, merr)
      call check(error, allocated(merr), "Linear mapping: infinite upper bound")
      if (allocated(error)) return
      call new_linear_mapping(lin, -0.5_wp, 2.0_wp, merr)
      call check(error, allocated(merr), "Linear mapping: negative lower bound")
      if (allocated(error)) return
      call new_linear_mapping(lin, 2.0_wp, 2.0_wp, merr)
      call check(error, allocated(merr), "Linear mapping: equal bounds")
      if (allocated(error)) return
      call new_linear_mapping(lin, 3.0_wp, 2.0_wp, merr)
      call check(error, allocated(merr), "Linear mapping: reversed bounds")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 0.0_wp, inf, 2.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: infinite upper bound")
      if (allocated(error)) return
      ! Documented constructor requirement: a = 2^m*(1 - 2^m + rmax - rmin) finite
      call new_handymod_mapping(hm, 0.0_wp, 0.8_wp*huge(1.0_wp), 2.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: denominator scale overflow")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, scale=inf)
      call check(error, allocated(merr), "Knowles mapping: infinite fixed scale")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, scale=1.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      ! -log(1 - u) = u + u^2/2 + u^3/3 + ...; at u = 1.25e-13 the cubic term is below
      ! 1e-26 relative, so u*(1 + u/2) is within 0.44 eps of a 60-digit reference,
      ! while the naive -log(1 - u) is off by 4e11 eps
      x = -1.0_wp + 1.0e-4_wp
      u = (0.5_wp*(1.0_wp + x))**3
      ref = u*(1.0_wp + 0.5_wp*u)
      call kn%transform(0, [x], [1.0_wp], r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, abs(r(1) - ref) <= 4.0_wp*epsilon(1.0_wp)*ref, &
         & "Knowles mapping: small radius retains relative accuracy")
      if (allocated(error)) return
      rx = kn%map(0, x)
      call check(error, abs(rx - ref) <= 4.0_wp*epsilon(1.0_wp)*ref, &
         & "Knowles map: small radius retains relative accuracy")
      if (allocated(error)) return
      call check_mapping_input([nan], [1.0_wp], "Direct input", r, w, merr)
      call check(error, allocated(merr), "Mapping input: NaN node")
      if (allocated(error)) return
      do kind = 1, 5
         call make_mapping(kind, mapping, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         rx = mapping%map(6, -2.0_wp)
         call check(error, .not. (rx == rx), "Mapping map: lower outside domain must be NaN")
         if (allocated(error)) return
         rx = mapping%map(6, 2.0_wp)
         call check(error, .not. (rx == rx), "Mapping map: upper outside domain must be NaN")
         if (allocated(error)) return
         deallocate (mapping)
      end do
   end subroutine test_mapping_validation

   !> Each mapping rejects a finite request whose exact output exceeds huge
   !>
   !> Linear, HandyMod, and Knowles take w_ref = huge at a node with |dr/dx| > 1
   !> (1.25, 1.52, 1.46), so the exact weight is unrepresentable however it is
   !> evaluated; Becke maps x = 0.5 with p = 0.8*huge to r = 2.4*huge
   subroutine test_mapping_overflow(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_linear_type) :: lin
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(moist_math_grid_radial_mapping_knowles_type) :: kn
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      real(wp) :: big

      big = huge(1.0_wp)
      ! dr/dx = (3 - 0.5)/2 = 1.25
      call new_linear_mapping(lin, 0.5_wp, 3.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call lin%transform(0, [0.0_wp], [big], r, w, merr)
      call check(error, allocated(merr), "Linear mapping: non-finite mapped weight")
      if (allocated(error)) return
      call new_becke_mapping(becke, merr, scale=0.8_wp*big)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call becke%transform(0, [0.5_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: non-finite mapped output")
      if (allocated(error)) return
      ! d = 5.5, c = 4, a = 10: dr/dx(0) = 5.5*2*10/8.5^2 = 1.52
      call new_handymod_mapping(hm, 0.5_wp, 6.0_wp, 2.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call hm%transform(0, [0.0_wp], [big], r, w, merr)
      call check(error, allocated(merr), "HandyMod mapping: non-finite mapped weight")
      if (allocated(error)) return
      ! t = 0.75, u = t^3: dr/dx(0.5) = 3*t^2/(2*(1 - u)) = 1.46
      call new_knowles_mapping(kn, 3.0_wp, merr, scale=1.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call kn%transform(0, [0.5_wp], [big], r, w, merr)
      call check(error, allocated(merr), "Knowles mapping: non-finite mapped weight")
   end subroutine test_mapping_overflow

   ! --------------------------------------------------------------------------
   ! Integrands and helpers
   ! --------------------------------------------------------------------------

   !> Spherically symmetric unit Gaussian f(r) = exp(-r^2)
   !>
   !> @param[in] r  Radius (bohr)
   pure function gaussian_unit_radial(r) result(val)
      !> Radius (bohr)
      real(wp), intent(in) :: r
      !> exp(-r^2)
      real(wp) :: val

      val = exp(-r*r)
   end function gaussian_unit_radial

   !> Constant integrand f(r) = 1
   !>
   !> @param[in] r  Radius (bohr)
   pure function one_radial(r) result(val)
      !> Radius (bohr)
      real(wp), intent(in) :: r
      !> 1
      real(wp) :: val

      val = 1.0_wp + 0.0_wp*r
   end function one_radial

   !> Antiderivative of r^2 exp(-r): -exp(-r)*(r^2 + 2r + 2)
   !>
   !> @param[in] r  Radius (bohr)
   pure function r2_exp_antiderivative(r) result(val)
      !> Radius (bohr)
      real(wp), intent(in) :: r
      !> Antiderivative value
      real(wp) :: val

      val = -exp(-r)*(r*r + 2.0_wp*r + 2.0_wp)
   end function r2_exp_antiderivative

   !> Build a rule of the requested kind behind the abstract type
   !>
   !> @param[in]  kind  Rule selector
   !> @param[out] rule  Allocated rule
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

   !> Generate a canonical rule and apply a mapping for element z
   !>
   !> @param[in]  kind     Rule selector
   !> @param[in]  n        Number of nodes
   !> @param[in]  mapping  Radial mapping
   !> @param[in]  z        Atomic number
   !> @param[out] r        Radii, shape (n)
   !> @param[out] w        dr weights, shape (n)
   !> @param[out] merr     Rule or mapping error
   subroutine mapped_nodes(kind, n, mapping, z, r, w, merr)
      !> Rule selector
      integer, intent(in) :: kind
      !> Number of nodes
      integer, intent(in) :: n
      !> Radial mapping
      class(moist_math_grid_radial_mapping_type), intent(in) :: mapping
      !> Atomic number
      integer, intent(in) :: z
      !> Radii
      real(wp), allocatable, intent(out) :: r(:)
      !> dr weights
      real(wp), allocatable, intent(out) :: w(:)
      !> Rule or mapping error
      type(mctc_error), allocatable, intent(out) :: merr

      class(moist_math_grid_radial_rule_type), allocatable :: rule
      real(wp), allocatable :: x(:), wx(:)

      call make_rule(kind, rule)
      call rule%generate(n, x, wx, merr)
      if (allocated(merr)) return
      call mapping%transform(z, x, wx, r, w, merr)
   end subroutine mapped_nodes

   !> Fill a recipe with a rule of the given kind, a mapping, and a node count
   !>
   !> @param[out] recipe   Recipe without cutoffs
   !> @param[in]  kind     Rule selector
   !> @param[in]  mapping  Radial mapping, copied
   !> @param[in]  npts     Number of nodes
   subroutine make_recipe(recipe, kind, mapping, npts)
      !> Recipe
      type(moist_math_grid_radial_recipe_type), intent(out) :: recipe
      !> Rule selector
      integer, intent(in) :: kind
      !> Radial mapping
      class(moist_math_grid_radial_mapping_type), intent(in) :: mapping
      !> Number of nodes
      integer, intent(in) :: npts

      call make_rule(kind, recipe%rule)
      allocate (recipe%mapping, source=mapping)
      recipe%npts = npts
   end subroutine make_recipe

   !> Largest relative deviation max_i |a_i - b_i|/|b_i|
   !>
   !> @param[in] a  Values under test
   !> @param[in] b  Reference values, nonzero
   pure function max_rel_dev(a, b) result(dev)
      !> Values under test
      real(wp), intent(in) :: a(:)
      !> Reference values
      real(wp), intent(in) :: b(:)
      !> Largest relative deviation
      real(wp) :: dev

      dev = maxval(abs(a - b)/abs(b))
   end function max_rel_dev

   ! --------------------------------------------------------------------------
   ! Linear mapping
   ! --------------------------------------------------------------------------

   !> Linear mapping with Gauss-Legendre integrates r^j exactly on [a, b] for j <= 2n-1
   subroutine test_linear_polynomial(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: a = 0.5_wp, b = 3.0_wp
      integer, parameter :: n = 5
      type(moist_math_grid_radial_mapping_linear_type) :: lin
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      real(wp) :: ref
      integer :: j

      call new_linear_mapping(lin, a, b, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call mapped_nodes(rule_gauss_legendre, n, lin, 0, r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(r >= a .and. r <= b), "Linear mapping: radius outside [lower, upper]")
      if (allocated(error)) return
      do j = 0, 2*n - 1
         ref = (b**(j + 1) - a**(j + 1))/real(j + 1, wp)
         call check(error, abs(sum(w*r**j) - ref) <= 1.0e-14_wp*ref, &
            & "Linear mapping: r^j not exact through degree 2n-1")
         if (allocated(error)) return
      end do
   end subroutine test_linear_polynomial

   !> Linear mapping of the canonical rule equals the rule generated on [lower, upper]
   subroutine test_linear_scaled_rule(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: a = 0.25_wp, b = 7.5_wp
      integer, parameter :: n = 9
      type(moist_math_grid_radial_mapping_linear_type) :: lin
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:), xs(:), ws(:)
      integer :: kind

      call new_linear_mapping(lin, a, b, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do kind = rule_chebyshev2, rule_gauss_legendre
         call mapped_nodes(kind, n, lin, 0, r, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call make_rule(kind, rule)
         call rule%generate(n, xs, ws, merr, lower=a, upper=b)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(abs(r - xs) <= two_ulp*b) .and. all(abs(w - ws) <= two_ulp*abs(ws)), &
            & "Linear mapping: differs from the rule generated on [lower, upper]")
         if (allocated(error)) return
      end do
   end subroutine test_linear_scaled_rule

   ! --------------------------------------------------------------------------
   ! Becke mapping
   ! --------------------------------------------------------------------------

   !> Becke mapping reproduces Gaussian and exponential moments on [0, inf)
   !>
   !> integral r^2 exp(-r^2) dr = sqrt(pi)/4, integral r^2 exp(-2r) dr = 1/4,
   !> for several scales
   subroutine test_becke_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: scales(3) = [0.7_wp, 1.0_wp, 2.0_wp]
      integer, parameter :: kinds(2) = [rule_chebyshev2, rule_gauss_legendre]
      integer, parameter :: sizes(2) = [200, 99]
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      integer :: is, ik

      do is = 1, size(scales)
         call new_becke_mapping(becke, merr, scale=scales(is))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do ik = 1, size(kinds)
            call mapped_nodes(kinds(ik), sizes(ik), becke, 0, r, w, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            call check(error, all(r >= 0.0_wp) .and. all(w > 0.0_wp), &
               & "Becke mapping: negative radius or non-positive weight")
            if (allocated(error)) return
            call check(error, abs(sum(w*r*r*exp(-r*r)) - 0.25_wp*sqrt(pi)) <= 1.0e-13_wp, &
               & "Becke mapping: Gaussian moment deviates from sqrt(pi)/4")
            if (allocated(error)) return
            call check(error, abs(sum(w*r*r*exp(-2.0_wp*r)) - 0.25_wp) <= 1.0e-13_wp, &
               & "Becke mapping: exponential moment deviates from 1/4")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_becke_integrals

   !> Element-dependent Becke scale equals radius_factor*covalent_rad(z); fixed scale ignores z
   subroutine test_becke_element_scale(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: zs(3) = [1, 6, 79]
      integer, parameter :: n = 20
      type(moist_math_grid_radial_mapping_becke_type) :: by_element, fixed
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:), rf(:), wf(:), r1(:), w1(:)
      real(wp) :: ratio
      integer :: iz

      call new_becke_mapping(by_element, merr, radius_factor=0.5_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call mapped_nodes(rule_chebyshev2, n, by_element, 1, r1, w1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do iz = 1, size(zs)
         call new_becke_mapping(fixed, merr, scale=0.5_wp*covalent_rad(zs(iz)))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call mapped_nodes(rule_chebyshev2, n, by_element, zs(iz), r, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call mapped_nodes(rule_chebyshev2, n, fixed, 0, rf, wf, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, max_rel_dev(r, rf) <= two_ulp .and. max_rel_dev(w, wf) <= two_ulp, &
            & "Becke mapping: element scale differs from radius_factor*covalent_rad(z)")
         if (allocated(error)) return
         ratio = covalent_rad(zs(iz))/covalent_rad(1)
         call check(error, max_rel_dev(r, ratio*r1) <= 4.0_wp*epsilon(1.0_wp) &
            & .and. max_rel_dev(w, ratio*w1) <= 4.0_wp*epsilon(1.0_wp), &
            & "Becke mapping: radii and weights must scale with covalent_rad(z)")
         if (allocated(error)) return
      end do

      ! A fixed scale does not look at z, even outside the covalent-radius table
      call new_becke_mapping(fixed, merr, scale=1.3_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call mapped_nodes(rule_chebyshev2, n, fixed, 0, r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call mapped_nodes(rule_chebyshev2, n, fixed, 200, rf, wf, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(r == rf) .and. all(w == wf), "Becke mapping: fixed scale must ignore z")
      if (allocated(error)) return

      call mapped_nodes(rule_chebyshev2, n, by_element, 0, r, w, merr)
      call check(error, allocated(merr), "Becke mapping: z = 0 has no covalent radius")
      if (allocated(error)) return
      call mapped_nodes(rule_chebyshev2, n, by_element, 119, r, w, merr)
      call check(error, allocated(merr), "Becke mapping: z = 119 has no covalent radius")
   end subroutine test_becke_element_scale

   !> Becke constructor and domain errors; x = 1 is the divergent endpoint
   subroutine test_becke_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      real(wp) :: nan, inf, rx

      nan = ieee_value(nan, ieee_quiet_nan)
      inf = ieee_value(inf, ieee_positive_inf)

      call new_becke_mapping(becke, merr)
      call check(error, allocated(merr), "Becke mapping: needs a scale")
      if (allocated(error)) return
      call new_becke_mapping(becke, merr, scale=1.0_wp, radius_factor=0.5_wp)
      call check(error, allocated(merr), "Becke mapping: scale and radius_factor are exclusive")
      if (allocated(error)) return
      call new_becke_mapping(becke, merr, scale=0.0_wp)
      call check(error, allocated(merr), "Becke mapping: scale must be > 0")
      if (allocated(error)) return
      call new_becke_mapping(becke, merr, radius_factor=-0.5_wp)
      call check(error, allocated(merr), "Becke mapping: radius_factor must be > 0")
      if (allocated(error)) return
      call new_becke_mapping(becke, merr, radius_factor=nan)
      call check(error, allocated(merr), "Becke mapping: radius_factor must be finite")
      if (allocated(error)) return
      call new_becke_mapping(becke, merr, scale=inf)
      call check(error, allocated(merr), "Becke mapping: scale must be finite")
      if (allocated(error)) return

      call new_becke_mapping(becke, merr, scale=2.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call becke%transform(0, [0.0_wp, 1.0_wp], [1.0_wp, 1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: x = 1 must be rejected, not dropped")
      if (allocated(error)) return
      call becke%transform(0, [1.5_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: x > 1 must be rejected")
      if (allocated(error)) return
      call becke%transform(0, [-1.5_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: x < -1 must be rejected")
      if (allocated(error)) return
      call becke%transform(0, [nan], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: NaN node must be rejected")
      if (allocated(error)) return
      call becke%transform(0, [0.0_wp], [nan], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: NaN weight must be rejected")
      if (allocated(error)) return
      call becke%transform(0, [0.0_wp, 0.5_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Becke mapping: node and weight sizes must match")
      if (allocated(error)) return

      ! x = -1 is a valid node at the nucleus: r = 0, dr/dx = p/2
      call becke%transform(0, [-1.0_wp], [1.0_wp], r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, r(1) == 0.0_wp .and. abs(w(1) - 1.0_wp) <= two_ulp, &
         & "Becke mapping: expected r = 0 and dr/dx = p/2 at x = -1")
      if (allocated(error)) return

      rx = becke%map(0, 1.0_wp)
      call check(error, .not. ieee_is_finite(rx) .and. rx > 0.0_wp, "Becke map: expected +inf at x = 1")
      if (allocated(error)) return
      rx = becke%map(0, 2.0_wp)
      call check(error, .not. (rx == rx), "Becke map: expected NaN outside [-1, 1]")
   end subroutine test_becke_errors

   ! --------------------------------------------------------------------------
   ! HandyMod mapping
   ! --------------------------------------------------------------------------

   !> HandyMod with Gauss-Legendre reproduces finite-interval volume and exponential moments
   subroutine test_handymod_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: rmins(3) = [0.5_wp, 0.0_wp, 0.5_wp]
      real(wp), parameter :: rmaxs(3) = [6.0_wp, 5.0_wp, 6.0_wp]
      real(wp), parameter :: ms(3) = [2.0_wp, 2.0_wp, 1.0_wp]
      integer, parameter :: n = 40
      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      real(wp) :: a, b, ref
      integer :: ic

      do ic = 1, size(ms)
         a = rmins(ic)
         b = rmaxs(ic)
         call new_handymod_mapping(hm, a, b, ms(ic), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call mapped_nodes(rule_gauss_legendre, n, hm, 0, r, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(r >= a .and. r <= b), "HandyMod mapping: radius outside [rmin, rmax]")
         if (allocated(error)) return
         ref = (b**3 - a**3)/3.0_wp
         call check(error, abs(sum(w*r*r) - ref) <= 1.0e-13_wp*ref, &
            & "HandyMod mapping: shell volume moment deviates")
         if (allocated(error)) return
         ref = r2_exp_antiderivative(b) - r2_exp_antiderivative(a)
         call check(error, abs(sum(w*r*r*exp(-r)) - ref) <= 1.0e-13_wp, &
            & "HandyMod mapping: exponential moment deviates")
         if (allocated(error)) return
      end do
   end subroutine test_handymod_integrals

   !> HandyMod endpoints: x = -1 -> rmin, x = 1 -> rmax; x = -1 rejected for m < 1
   subroutine test_handymod_endpoints(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)

      call new_handymod_mapping(hm, 0.5_wp, 6.0_wp, 2.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call hm%transform(0, [-1.0_wp, 1.0_wp], [1.0_wp, 1.0_wp], r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, r(1) == 0.5_wp .and. abs(r(2) - 6.0_wp) <= 6.0_wp*two_ulp, &
         & "HandyMod mapping: endpoints must map to rmin and rmax")
      if (allocated(error)) return
      call check(error, w(1) == 0.0_wp .and. w(2) > 0.0_wp .and. ieee_is_finite(w(2)), &
         & "HandyMod mapping: m = 2 has dr/dx = 0 at x = -1 and finite at x = 1")
      if (allocated(error)) return

      call new_handymod_mapping(hm, 0.5_wp, 6.0_wp, 0.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call hm%transform(0, [-1.0_wp, 0.0_wp], [1.0_wp, 1.0_wp], r, w, merr)
      call check(error, allocated(merr), "HandyMod mapping: x = -1 must be rejected for m < 1")
      if (allocated(error)) return
      call hm%transform(0, [0.0_wp, 1.0_wp], [1.0_wp, 1.0_wp], r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(ieee_is_finite(w)) .and. abs(r(2) - 6.0_wp) <= 6.0_wp*two_ulp, &
         & "HandyMod mapping: m < 1 stays finite at x = 1")
   end subroutine test_handymod_endpoints

   !> HandyMod constructor and domain errors, including rmax - rmin <= 2^m - 1
   subroutine test_handymod_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      real(wp) :: nan, inf

      nan = ieee_value(nan, ieee_quiet_nan)
      inf = ieee_value(inf, ieee_positive_inf)

      call new_handymod_mapping(hm, 0.0_wp, 3.0_wp, 2.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: rmax - rmin = 2^m - 1 must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 0.0_wp, 3.0_wp + 1.0e-12_wp, 2.0_wp, merr)
      call check(error, .not. allocated(merr), "HandyMod mapping: rmax - rmin > 2^m - 1 is valid")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 0.0_wp, 5.0_wp, 0.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: m = 0 must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 0.0_wp, 5.0_wp, -1.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: m < 0 must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 0.0_wp, 5.0_wp, inf, merr)
      call check(error, allocated(merr), "HandyMod mapping: infinite m must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 0.0_wp, 5.0_wp, 2000.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: 2^m overflow must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, 5.0_wp, 5.0_wp, 1.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: rmax <= rmin must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, -0.5_wp, 5.0_wp, 1.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: rmin < 0 must be rejected")
      if (allocated(error)) return
      call new_handymod_mapping(hm, nan, 5.0_wp, 1.0_wp, merr)
      call check(error, allocated(merr), "HandyMod mapping: NaN rmin must be rejected")
      if (allocated(error)) return

      call new_handymod_mapping(hm, 0.0_wp, 5.0_wp, 2.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call hm%transform(0, [1.5_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "HandyMod mapping: x > 1 must be rejected")
      if (allocated(error)) return
      call hm%transform(0, [nan], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "HandyMod mapping: NaN node must be rejected")
      if (allocated(error)) return
      call hm%transform(0, [0.0_wp], [1.0_wp, 1.0_wp], r, w, merr)
      call check(error, allocated(merr), "HandyMod mapping: node and weight sizes must match")
   end subroutine test_handymod_errors

   ! --------------------------------------------------------------------------
   ! Knowles mapping
   ! --------------------------------------------------------------------------

   !> Knowles mapping reproduces exponential and Gaussian moments on [0, inf)
   !>
   !> integral r^2 exp(-r) dr = 2, integral r^2 exp(-r^2) dr = sqrt(pi)/4
   subroutine test_knowles_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: scales(3) = [3.0_wp, 5.0_wp, 7.0_wp]
      integer, parameter :: n = 200
      type(moist_math_grid_radial_mapping_knowles_type) :: kn
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:)
      integer :: is, kind

      do is = 1, size(scales)
         call new_knowles_mapping(kn, 3.0_wp, merr, scale=scales(is))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do kind = rule_chebyshev2, rule_gauss_legendre
            call mapped_nodes(kind, n, kn, 0, r, w, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            call check(error, all(r >= 0.0_wp) .and. all(w > 0.0_wp), &
               & "Knowles mapping: negative radius or non-positive weight")
            if (allocated(error)) return
            call check(error, abs(sum(w*r*r*exp(-r*r)) - 0.25_wp*sqrt(pi)) <= 1.0e-13_wp, &
               & "Knowles mapping: Gaussian moment deviates from sqrt(pi)/4")
            if (allocated(error)) return
            ! exp(-r) becomes ~(1 - t)^R near t = 1, so a small R leaves a weak
            ! endpoint singularity and algebraic convergence; checked for R >= 5
            if (scales(is) < 5.0_wp) cycle
            if (kind == rule_midpoint) then
               ! Midpoint converges algebraically on the exponential tail
               call check(error, abs(sum(w*r*r*exp(-r)) - 2.0_wp) <= 1.0e-7_wp, &
                  & "Knowles mapping: midpoint exponential moment deviates from 2")
            else
               call check(error, abs(sum(w*r*r*exp(-r)) - 2.0_wp) <= 1.0e-13_wp, &
                  & "Knowles mapping: exponential moment deviates from 2")
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_knowles_integrals

   !> Midpoint + Knowles(k = 3) is the Mura-Knowles recipe
   !>
   !> With t_i = (i - 1/2)/n: r_i = -R*log(1 - t_i^3) and
   !> w_i = 3*R*t_i^2/((1 - t_i^3)*n); the reference logarithm uses a series
   !> for small arguments so it stays accurate near the nucleus
   subroutine test_knowles_mura_knowles(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: big_r = 5.0_wp
      integer, parameter :: n = 64
      type(moist_math_grid_radial_mapping_knowles_type) :: kn
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:), r_ref(:), w_ref(:)
      real(wp) :: t, u, s
      integer :: i, j

      call new_knowles_mapping(kn, 3.0_wp, merr, scale=big_r)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call mapped_nodes(rule_midpoint, n, kn, 0, r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      allocate (r_ref(n), w_ref(n))
      do i = 1, n
         t = (real(i, wp) - 0.5_wp)/real(n, wp)
         u = t*t*t
         if (u < 0.05_wp) then
            s = 0.0_wp
            do j = 40, 1, -1
               s = s*u + 1.0_wp/real(j, wp)
            end do
            r_ref(i) = big_r*u*s
         else
            r_ref(i) = -big_r*log(1.0_wp - u)
         end if
         w_ref(i) = 3.0_wp*big_r*t*t/((1.0_wp - u)*real(n, wp))
      end do
      call check(error, max_rel_dev(r, r_ref) <= 1.0e-14_wp, &
         & "Knowles mapping: midpoint radii differ from -R*log(1 - t^3)")
      if (allocated(error)) return
      call check(error, max_rel_dev(w, w_ref) <= 1.0e-14_wp, &
         & "Knowles mapping: midpoint weights differ from 3*R*t^2/((1 - t^3)*n)")
   end subroutine test_knowles_mura_knowles

   !> Knowles element table: z selects R from the table; z outside the table is an error
   subroutine test_knowles_element_table(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: table(3) = [1.0_wp, 2.5_wp, 4.0_wp]
      integer, parameter :: n = 30
      type(moist_math_grid_radial_mapping_knowles_type) :: by_element, fixed
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:), rf(:), wf(:)
      integer :: z

      call new_knowles_mapping(by_element, 3.0_wp, merr, element_scale=table)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do z = 1, size(table)
         call new_knowles_mapping(fixed, 3.0_wp, merr, scale=table(z))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call mapped_nodes(rule_gauss_legendre, n, by_element, z, r, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call mapped_nodes(rule_gauss_legendre, n, fixed, 0, rf, wf, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, max_rel_dev(r, rf) <= two_ulp .and. max_rel_dev(w, wf) <= two_ulp, &
            & "Knowles mapping: element table entry differs from the same fixed scale")
         if (allocated(error)) return
      end do
      call mapped_nodes(rule_gauss_legendre, n, by_element, 4, r, w, merr)
      call check(error, allocated(merr), "Knowles mapping: z beyond the table must be rejected")
      if (allocated(error)) return
      call mapped_nodes(rule_gauss_legendre, n, by_element, 0, r, w, merr)
      call check(error, allocated(merr), "Knowles mapping: z = 0 must be rejected")
      if (allocated(error)) return
      call mapped_nodes(rule_gauss_legendre, n, fixed, 200, r, w, merr)
      call check(error, .not. allocated(merr), "Knowles mapping: fixed scale must ignore z")
   end subroutine test_knowles_element_table

   !> Knowles constructor and domain errors; x = 1 is the divergent endpoint
   subroutine test_knowles_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_knowles_type) :: kn
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: r(:), w(:), empty(:)
      real(wp) :: nan, rx

      nan = ieee_value(nan, ieee_quiet_nan)
      allocate (empty(0))

      call new_knowles_mapping(kn, 0.0_wp, merr, scale=1.0_wp)
      call check(error, allocated(merr), "Knowles mapping: k = 0 must be rejected")
      if (allocated(error)) return
      call new_knowles_mapping(kn, nan, merr, scale=1.0_wp)
      call check(error, allocated(merr), "Knowles mapping: NaN k must be rejected")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr)
      call check(error, allocated(merr), "Knowles mapping: needs a scale")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, scale=1.0_wp, element_scale=[1.0_wp])
      call check(error, allocated(merr), "Knowles mapping: scale and element_scale are exclusive")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, scale=-1.0_wp)
      call check(error, allocated(merr), "Knowles mapping: scale must be > 0")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, element_scale=[1.0_wp, 0.0_wp])
      call check(error, allocated(merr), "Knowles mapping: element scales must be > 0")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, element_scale=[nan])
      call check(error, allocated(merr), "Knowles mapping: element scales must be finite")
      if (allocated(error)) return
      call new_knowles_mapping(kn, 3.0_wp, merr, element_scale=empty)
      call check(error, allocated(merr), "Knowles mapping: empty element table must be rejected")
      if (allocated(error)) return

      call new_knowles_mapping(kn, 3.0_wp, merr, scale=5.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call kn%transform(0, [0.0_wp, 1.0_wp], [1.0_wp, 1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Knowles mapping: x = 1 must be rejected, not dropped")
      if (allocated(error)) return
      call kn%transform(0, [1.5_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Knowles mapping: x > 1 must be rejected")
      if (allocated(error)) return
      call kn%transform(0, [-1.0_wp], [1.0_wp], r, w, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, r(1) == 0.0_wp .and. w(1) == 0.0_wp, &
         & "Knowles mapping: k = 3 has r = 0 and dr/dx = 0 at x = -1")
      if (allocated(error)) return
      rx = kn%map(0, 1.0_wp)
      call check(error, .not. ieee_is_finite(rx) .and. rx > 0.0_wp, "Knowles map: expected +inf at x = 1")
      if (allocated(error)) return

      call new_knowles_mapping(kn, 0.5_wp, merr, scale=5.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call kn%transform(0, [-1.0_wp], [1.0_wp], r, w, merr)
      call check(error, allocated(merr), "Knowles mapping: x = -1 must be rejected for k < 1")
   end subroutine test_knowles_errors

   ! --------------------------------------------------------------------------
   ! All mappings
   ! --------------------------------------------------------------------------

   !> Build one of five configured mappings behind the abstract type
   !>
   !> @param[in]  kind     Mapping selector (1..5)
   !> @param[out] mapping  Allocated mapping
   !> @param[out] merr     Constructor error
   subroutine make_mapping(kind, mapping, merr)
      !> Mapping selector
      integer, intent(in) :: kind
      !> Allocated mapping
      class(moist_math_grid_radial_mapping_type), allocatable, intent(out) :: mapping
      !> Constructor error
      type(mctc_error), allocatable, intent(out) :: merr

      type(moist_math_grid_radial_mapping_linear_type) :: lin
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(moist_math_grid_radial_mapping_knowles_type) :: kn

      select case (kind)
      case (1)
         call new_linear_mapping(lin, 0.5_wp, 3.0_wp, merr)
         allocate (mapping, source=lin)
      case (2)
         call new_becke_mapping(becke, merr, scale=1.3_wp)
         allocate (mapping, source=becke)
      case (3)
         call new_becke_mapping(becke, merr, radius_factor=0.5_wp)
         allocate (mapping, source=becke)
      case (4)
         call new_handymod_mapping(hm, 0.5_wp, 6.0_wp, 2.0_wp, merr)
         allocate (mapping, source=hm)
      case default
         call new_knowles_mapping(kn, 3.0_wp, merr, element_scale=[1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp])
         allocate (mapping, source=kn)
      end select
   end subroutine make_mapping

   !> The forward map agrees with the radii from transform for every mapping
   subroutine test_map_matches_transform(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 16, z = 6
      class(moist_math_grid_radial_mapping_type), allocatable :: mapping
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), wx(:), r(:), w(:)
      integer :: kind, i

      call make_rule(rule_gauss_legendre, rule)
      call rule%generate(n, x, wx, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do kind = 1, 5
         call make_mapping(kind, mapping, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call mapping%transform(z, x, wx, r, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do i = 1, n
            call check(error, abs(mapping%map(z, x(i)) - r(i)) <= two_ulp*abs(r(i)), &
               & "Mapping: forward map differs from transform radii")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_map_matches_transform

   !> A polymorphic copy transforms like its source after the source is gone
   subroutine test_mapping_copy(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 12, z = 5
      class(moist_math_grid_radial_mapping_type), allocatable :: source, copy
      class(moist_math_grid_radial_rule_type), allocatable :: rule
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:), wx(:), r0(:), w0(:), r1(:), w1(:)
      integer :: kind

      call make_rule(rule_chebyshev2, rule)
      call rule%generate(n, x, wx, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      do kind = 1, 5
         call make_mapping(kind, source, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call source%transform(z, x, wx, r0, w0, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         allocate (copy, source=source)
         deallocate (source)
         call copy%transform(z, x, wx, r1, w1, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(r1 == r0) .and. all(w1 == w0), &
            & "Mapping: polymorphic copy transforms differently from its source")
         if (allocated(error)) return
         deallocate (copy)
      end do
   end subroutine test_mapping_copy

   ! --------------------------------------------------------------------------
   ! Radial recipe and grid
   ! --------------------------------------------------------------------------

   !> Cutoffs drop r < rcut_lower and r > rcut_upper; counts record the truncation
   !>
   !> Cutoffs placed exactly on nodes check the strict comparisons
   subroutine test_recipe_cutoffs(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 40, z = 6
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_recipe_type) :: recipe
      type(moist_math_grid_radial_type) :: full, cut
      type(mctc_error), allocatable :: merr
      logical, allocatable :: keep(:)
      real(wp) :: lo, hi

      call new_becke_mapping(becke, merr, scale=1.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_recipe(recipe, rule_chebyshev2, becke, n)
      call new_radial_grid(full, recipe, z, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, full%npts, n, "Radial grid: all nodes retained without cutoffs")
      if (allocated(error)) return
      call check(error, full%npts_requested, n, "Radial grid: requested count")
      if (allocated(error)) return
      call check(error, full%transform, transform_quadrature, "Radial grid: tag must be transform_quadrature")
      if (allocated(error)) return
      call check(error, full%spacing == 0.0_wp .and. size(full%r) == n .and. size(full%w) == n, &
         & "Radial grid: quadrature grid records no spacing and owns npts nodes")
      if (allocated(error)) return

      ! Chebyshev-II + Becke radii descend, so r(10) > r(30)
      hi = full%r(10)
      lo = full%r(30)

      allocate (recipe%rcut_lower, recipe%rcut_upper)
      recipe%rcut_lower = lo
      recipe%rcut_upper = hi
      call new_radial_grid(cut, recipe, z, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      keep = full%r >= lo .and. full%r <= hi
      call check(error, cut%npts, 21, "Radial grid: both cutoffs keep nodes 10..30")
      if (allocated(error)) return
      call check(error, cut%npts_requested, n, "Radial grid: requested count survives cutoffs")
      if (allocated(error)) return
      call check(error, all(cut%r == pack(full%r, keep)) .and. all(cut%w == pack(full%w, keep)), &
         & "Radial grid: retained nodes must be the uncut nodes inside the cutoffs, in order")
      if (allocated(error)) return
      call check(error, cut%r(1) == hi .and. cut%r(cut%npts) == lo, &
         & "Radial grid: nodes exactly on a cutoff are retained")
      if (allocated(error)) return

      deallocate (recipe%rcut_upper)
      call new_radial_grid(cut, recipe, z, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, cut%npts, 30, "Radial grid: lower cutoff keeps nodes 1..30")
      if (allocated(error)) return
      call check(error, all(cut%r >= lo) .and. cut%npts_requested == n, "Radial grid: lower cutoff bookkeeping")
      if (allocated(error)) return

      deallocate (recipe%rcut_lower)
      allocate (recipe%rcut_upper)
      recipe%rcut_upper = hi
      call new_radial_grid(cut, recipe, z, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, cut%npts, 31, "Radial grid: upper cutoff keeps nodes 10..40")
      if (allocated(error)) return
      call check(error, all(cut%r <= hi) .and. cut%npts_requested == n, "Radial grid: upper cutoff bookkeeping")
   end subroutine test_recipe_cutoffs

   !> Invalid recipes and failing components are errors
   subroutine test_recipe_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_recipe_type) :: recipe
      type(moist_math_grid_radial_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: nan

      nan = ieee_value(nan, ieee_quiet_nan)

      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: recipe without rule must be rejected")
      if (allocated(error)) return
      call make_rule(rule_chebyshev2, recipe%rule)
      recipe%npts = 10
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: recipe without mapping must be rejected")
      if (allocated(error)) return

      call new_becke_mapping(becke, merr, radius_factor=0.5_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      allocate (recipe%mapping, source=becke)
      deallocate (recipe%rule)
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: mapping alone does not replace a rule")
      if (allocated(error)) return
      call make_recipe(recipe, rule_chebyshev2, becke, 0)
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: npts = 0 must be rejected")
      if (allocated(error)) return

      recipe%npts = 10
      call new_radial_grid(grid, recipe, 0, merr)
      call check(error, allocated(merr), "Radial grid: mapping error must propagate")
      if (allocated(error)) return

      allocate (recipe%rcut_lower)
      recipe%rcut_lower = nan
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: NaN rcut_lower must be rejected")
      if (allocated(error)) return

      allocate (recipe%rcut_upper)
      recipe%rcut_lower = 2.0_wp
      recipe%rcut_upper = nan
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: NaN rcut_upper must be rejected")
      if (allocated(error)) return
      recipe%rcut_upper = ieee_value(nan, ieee_positive_inf)
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: infinite rcut_upper must be rejected")
      if (allocated(error)) return

      recipe%rcut_upper = 2.0_wp
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: rcut_lower >= rcut_upper must be rejected")
      if (allocated(error)) return

      deallocate (recipe%rcut_lower)
      recipe%rcut_upper = 1.0e-6_wp
      call new_radial_grid(grid, recipe, 6, merr)
      call check(error, allocated(merr), "Radial grid: cutoffs removing every node must be rejected")
   end subroutine test_recipe_errors

   !> Volume integrals through the grid apply 4*pi*r^2 to the dr weights
   !>
   !> integral exp(-r^2) dV = pi^(3/2); integral_{|r| < R} dV = 4*pi*R^3/3
   subroutine test_recipe_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_mapping_handymod_type) :: hm
      type(moist_math_grid_radial_mapping_knowles_type) :: kn
      type(moist_math_grid_radial_recipe_type) :: recipe
      type(moist_math_grid_radial_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:)
      real(wp) :: res, field_res, ref

      ref = pi*sqrt(pi)

      call new_becke_mapping(becke, merr, scale=1.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_recipe(recipe, rule_chebyshev2, becke, 400)
      call new_radial_grid(grid, recipe, 1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call grid%integrate(gaussian_unit_radial, res)
      call check(error, abs(res - ref) <= 1.0e-13_wp, "Chebyshev-II + Becke: Gaussian volume integral")
      if (allocated(error)) return
      f = exp(-grid%r*grid%r)
      call grid%integrate_field(f, field_res)
      call check(error, abs(field_res - ref) <= 1.0e-13_wp, "Chebyshev-II + Becke: integrate_field")
      if (allocated(error)) return

      call new_knowles_mapping(kn, 3.0_wp, merr, scale=5.0_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_recipe(recipe, rule_gauss_legendre, kn, 200)
      call new_radial_grid(grid, recipe, 1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call grid%integrate(gaussian_unit_radial, res)
      call check(error, abs(res - ref) <= 1.0e-13_wp, "Gauss-Legendre + Knowles: Gaussian volume integral")
      if (allocated(error)) return

      call new_handymod_mapping(hm, 0.0_wp, 5.0_wp, 2.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_recipe(recipe, rule_gauss_legendre, hm, 40)
      call new_radial_grid(grid, recipe, 1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call grid%integrate(one_radial, res)
      ref = 4.0_wp*pi*125.0_wp/3.0_wp
      call check(error, abs(res - ref) <= 1.0e-13_wp*ref, "Gauss-Legendre + HandyMod: ball volume")
   end subroutine test_recipe_integrals

   !> One recipe gives element-scaled grids; a copied recipe is independent
   subroutine test_recipe_element(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 30
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_recipe_type) :: recipe, copy
      type(moist_math_grid_radial_type) :: grid_h, grid_c
      type(mctc_error), allocatable :: merr
      real(wp) :: ratio

      call new_becke_mapping(becke, merr, radius_factor=0.5_wp)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_recipe(recipe, rule_chebyshev2, becke, n)
      copy = recipe
      recipe%npts = 5
      deallocate (recipe%mapping)

      call new_radial_grid(grid_h, copy, 1, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call new_radial_grid(grid_c, copy, 6, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, grid_h%npts == n .and. grid_c%npts == n, &
         & "Radial recipe: a copy must not see later changes to its source")
      if (allocated(error)) return
      ratio = covalent_rad(6)/covalent_rad(1)
      call check(error, max_rel_dev(grid_c%r, ratio*grid_h%r) <= 4.0_wp*epsilon(1.0_wp) &
         & .and. max_rel_dev(grid_c%w, ratio*grid_h%w) <= 4.0_wp*epsilon(1.0_wp), &
         & "Radial recipe: element grids must scale with covalent_rad(z)")
   end subroutine test_recipe_element

   ! --------------------------------------------------------------------------
   ! Uniform radial pair
   ! --------------------------------------------------------------------------

   !> Uniform pair: midpoint nodes, constant weights, DST-IV tags and spacings
   !>
   !> integral exp(-r^2) dV = pi^(3/2) on the r-grid; its transform
   !> pi^(3/2)*exp(-k^2/4) integrates to (2*pi)^3 on the k-grid
   subroutine test_uniform_pair(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 2000
      real(wp), parameter :: dr = 0.01_wp
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:)
      real(wp) :: dk, res, ref
      integer :: i

      call new_uniform_radial_pair(rgrid, kgrid, n, dr, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      dk = pi/(real(n, wp)*dr)
      call check(error, rgrid%transform == transform_dst4 .and. kgrid%transform == transform_dst4, &
         & "Uniform pair: both grids must be tagged transform_dst4")
      if (allocated(error)) return
      call check(error, rgrid%npts == n .and. kgrid%npts == n .and. rgrid%npts_requested == n &
         & .and. kgrid%npts_requested == n, "Uniform pair: node counts")
      if (allocated(error)) return
      call check(error, rgrid%spacing == dr .and. abs(kgrid%spacing - dk) <= two_ulp*dk, &
         & "Uniform pair: recorded spacings must be dr and pi/(npts*dr)")
      if (allocated(error)) return
      call check(error, all(rgrid%w == rgrid%spacing) .and. all(kgrid%w == kgrid%spacing), &
         & "Uniform pair: weights must equal the spacing")
      if (allocated(error)) return
      do i = 1, n
         call check(error, abs(rgrid%r(i) - (real(i, wp) - 0.5_wp)*dr) <= two_ulp*rgrid%r(i) &
            & .and. abs(kgrid%r(i) - (real(i, wp) - 0.5_wp)*dk) <= 4.0_wp*epsilon(1.0_wp)*kgrid%r(i), &
            & "Uniform pair: nodes must be (i - 1/2)*spacing")
         if (allocated(error)) return
      end do

      call rgrid%integrate(gaussian_unit_radial, res)
      call check(error, abs(res - pi*sqrt(pi)) <= 1.0e-13_wp, "Uniform pair: r-grid Gaussian volume integral")
      if (allocated(error)) return
      f = pi*sqrt(pi)*exp(-0.25_wp*kgrid%r*kgrid%r)
      call kgrid%integrate_field(f, res)
      ref = 8.0_wp*pi**3
      call check(error, abs(res - ref) <= 1.0e-13_wp*ref, "Uniform pair: k-grid transform volume integral")
   end subroutine test_uniform_pair

   !> Uniform pair rejects invalid counts and spacings
   subroutine test_uniform_pair_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(mctc_error), allocatable :: merr
      real(wp) :: nan, inf

      nan = ieee_value(nan, ieee_quiet_nan)
      inf = ieee_value(inf, ieee_positive_inf)

      call new_uniform_radial_pair(rgrid, kgrid, 0, 0.1_wp, merr)
      call check(error, allocated(merr), "Uniform pair: npts = 0 must be rejected")
      if (allocated(error)) return
      call new_uniform_radial_pair(rgrid, kgrid, 10, 0.0_wp, merr)
      call check(error, allocated(merr), "Uniform pair: dr = 0 must be rejected")
      if (allocated(error)) return
      call new_uniform_radial_pair(rgrid, kgrid, 10, -0.1_wp, merr)
      call check(error, allocated(merr), "Uniform pair: dr < 0 must be rejected")
      if (allocated(error)) return
      call new_uniform_radial_pair(rgrid, kgrid, 10, nan, merr)
      call check(error, allocated(merr), "Uniform pair: NaN dr must be rejected")
      if (allocated(error)) return
      call new_uniform_radial_pair(rgrid, kgrid, 10, inf, merr)
      call check(error, allocated(merr), "Uniform pair: infinite dr must be rejected")
   end subroutine test_uniform_pair_errors

end module test_math_grid_radial
