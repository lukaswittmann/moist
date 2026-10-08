!> Test suite for radial mappings, the radial recipe, and the concrete radial grid
!>
!>   - Mappings: analytic mapped integrals (polynomials, Gaussians,
!>     exponentials), element-dependent scales, Mura-Knowles nodes
!>   - Radial recipe: lower and upper cutoffs, volume integrals, element
!>     dependence
!>   - Uniform radial pair: r- and k-grid volume integrals
module test_math_grid_radial
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_data_atomicrad, only: covalent_rad
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_type, &
      & moist_math_grid_radial_rule_chebyshev2_type, new_chebyshev2_rule, &
      & moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule, &
      & moist_math_grid_radial_rule_gauss_legendre_type, new_gauss_legendre_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_type, &
      & moist_math_grid_radial_mapping_linear_type, new_linear_mapping, &
      & moist_math_grid_radial_mapping_becke_type, new_becke_mapping, &
      & moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping, &
      & moist_math_grid_radial_mapping_knowles_type, new_knowles_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, &
      & moist_math_grid_radial_recipe_type, new_radial_grid, new_uniform_radial_pair
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
         new_unittest("linear_polynomial_exact", test_linear_polynomial), &
         new_unittest("becke_mapped_integrals", test_becke_integrals), &
         new_unittest("becke_element_scale", test_becke_element_scale), &
         new_unittest("handymod_finite_integrals", test_handymod_integrals), &
         new_unittest("knowles_mapped_integrals", test_knowles_integrals), &
         new_unittest("knowles_mura_knowles_nodes", test_knowles_mura_knowles), &
         new_unittest("recipe_cutoffs", test_recipe_cutoffs), &
         new_unittest("recipe_volume_integrals", test_recipe_integrals), &
         new_unittest("recipe_element_dependence", test_recipe_element), &
         new_unittest("uniform_pair_contract", test_uniform_pair) &
         ]
   end subroutine collect_math_grid_radial

   !* ----------------------------- Integrands and helpers ---------------------------- *!

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

   !> Element-wise relative check |a_i - b_i| <= tol*|b_i| against a finite reference
   !>
   !> @param[out] error  Test failure
   !> @param[in]  a      Values under test
   !> @param[in]  b      Reference values, nonzero
   !> @param[in]  tol    Relative tolerance
   !> @param[in]  more   Failure context
   subroutine check_rel_dev(error, a, b, tol, more)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Values under test
      real(wp), intent(in) :: a(:)
      !> Reference values
      real(wp), intent(in) :: b(:)
      !> Relative tolerance
      real(wp), intent(in) :: tol
      !> Failure context
      character(len=*), intent(in) :: more

      integer :: i

      do i = 1, size(b)
         call check(error, ieee_is_finite(b(i)), "reference value is not finite")
         if (allocated(error)) return
         call check(error, a(i), b(i), thr=tol*abs(b(i)), more=more)
         if (allocated(error)) return
      end do
   end subroutine check_rel_dev

   !* --------------------------------- Linear mapping -------------------------------- *!

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

   !* --------------------------------- Becke mapping --------------------------------- *!

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
         call check_rel_dev(error, r, rf, two_ulp, &
            & "Becke mapping: element scale differs from radius_factor*covalent_rad(z)")
         if (allocated(error)) return
         call check_rel_dev(error, w, wf, two_ulp, &
            & "Becke mapping: element scale differs from radius_factor*covalent_rad(z)")
         if (allocated(error)) return
         ratio = covalent_rad(zs(iz))/covalent_rad(1)
         call check_rel_dev(error, r, ratio*r1, 4.0_wp*epsilon(1.0_wp), &
            & "Becke mapping: radii and weights must scale with covalent_rad(z)")
         if (allocated(error)) return
         call check_rel_dev(error, w, ratio*w1, 4.0_wp*epsilon(1.0_wp), &
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
   end subroutine test_becke_element_scale

   !* -------------------------------- HandyMod mapping ------------------------------- *!

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

   !* -------------------------------- Knowles mapping -------------------------------- *!

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
      call check_rel_dev(error, r, r_ref, 1.0e-14_wp, &
         & "Knowles mapping: midpoint radii differ from -R*log(1 - t^3)")
      if (allocated(error)) return
      call check_rel_dev(error, w, w_ref, 1.0e-14_wp, &
         & "Knowles mapping: midpoint weights differ from 3*R*t^2/((1 - t^3)*n)")
   end subroutine test_knowles_mura_knowles

   !* ----------------------------- Radial recipe and grid ---------------------------- *!

   !> Cutoffs drop r < rcut_lower and r > rcut_upper
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
      call check(error, all(cut%r >= lo), "Radial grid: lower cutoff drops r < rcut_lower")
      if (allocated(error)) return

      deallocate (recipe%rcut_lower)
      allocate (recipe%rcut_upper)
      recipe%rcut_upper = hi
      call new_radial_grid(cut, recipe, z, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(cut%r <= hi), "Radial grid: upper cutoff drops r > rcut_upper")
   end subroutine test_recipe_cutoffs

   !> Volume integrals through the grid apply 4*pi*r^2 to the dr weights
   !>
   !> integral exp(-r^2) dV = pi^(3/2); integral_{|r| < R} dV = 4*pi*R^3/3; a field of
   !> the wrong length is rejected, and the next valid call succeeds
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
      ! A short field is rejected; the valid call below still succeeds
      call grid%integrate_field(f(2:), field_res, merr)
      call check(error, allocated(merr), "integrate_field accepted a short field")
      if (allocated(error)) return
      call grid%integrate_field(f, field_res, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
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
      ratio = covalent_rad(6)/covalent_rad(1)
      call check_rel_dev(error, grid_c%r, ratio*grid_h%r, 4.0_wp*epsilon(1.0_wp), &
         & "Radial recipe: element grids must scale with covalent_rad(z)")
      if (allocated(error)) return
      call check_rel_dev(error, grid_c%w, ratio*grid_h%w, 4.0_wp*epsilon(1.0_wp), &
         & "Radial recipe: element grids must scale with covalent_rad(z)")
   end subroutine test_recipe_element

   !* ------------------------------ Uniform radial pair ------------------------------ *!

   !> Uniform pair: r- and k-grid volume integrals
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
      real(wp) :: res, ref

      call new_uniform_radial_pair(rgrid, kgrid, n, dr, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call rgrid%integrate(gaussian_unit_radial, res)
      call check(error, abs(res - pi*sqrt(pi)) <= 1.0e-13_wp, "Uniform pair: r-grid Gaussian volume integral")
      if (allocated(error)) return
      f = pi*sqrt(pi)*exp(-0.25_wp*kgrid%r*kgrid%r)
      call kgrid%integrate_field(f, res, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      ref = 8.0_wp*pi**3
      call check(error, abs(res - ref) <= 1.0e-13_wp*ref, "Uniform pair: k-grid transform volume integral")
   end subroutine test_uniform_pair

end module test_math_grid_radial
