!> Test suite for the atomic grid layer
!>
!> - Shell policies: sector and arc requests, edge membership, hard floors
!>   and caps, target fallback to the largest admissible rule
!> - Atomic composition: point layout, rules kept across cache growth, ball
!>   and shell volumes, Gaussian and exponential integrals across radial
!>   scales and elements
!> - Recipes: element overrides, first listing wins
module test_math_grid_atomic
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_helpers, only: check_moist_error
   use moist_math_grid_angular_lebedev, only: grid_size, lebedev_negative_weight_sizes, &
      & lebedev_degree_table
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_type, &
      & moist_math_grid_radial_rule_chebyshev2_type, new_chebyshev2_rule, &
      & moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule, &
      & moist_math_grid_radial_rule_gauss_legendre_type, new_gauss_legendre_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_type, &
      & moist_math_grid_radial_mapping_linear_type, new_linear_mapping, &
      & moist_math_grid_radial_mapping_becke_type, new_becke_mapping, &
      & moist_math_grid_radial_mapping_knowles_type, new_knowles_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_radial_grid
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, moist_math_grid_angular_request_type, &
      & moist_math_grid_angular_generator_lebedev_type, new_lebedev_generator, new_lebedev_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_shell_type, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy, &
      & moist_math_grid_atomic_shell_sector_type, new_sector_shell_policy, &
      & moist_math_grid_atomic_shell_arc_type, new_arc_shell_policy, &
      & moist_math_grid_atomic_recipe_type, moist_math_grid_atomic_recipe_override_type, &
      & element_override_index, get_element_recipe, default_element_recipes
   use moist_math_grid_atomic_grid, only: moist_math_grid_atomic_type, new_atomic_grid
   implicit none(type, external)
   private

   public :: collect_math_grid_atomic

   !> Rule selectors for `make_rule`
   integer, parameter :: rule_chebyshev2 = 1, rule_midpoint = 2, rule_gauss_legendre = 3

   !> Center of the off-center Gaussian (bohr)
   real(wp), parameter :: off_center(3) = [0.3_wp, -0.2_wp, 0.4_wp]

   !> Repeating requests exercise cache preservation after growth
   type, extends(moist_math_grid_atomic_shell_type) :: repeating_shell_type
   contains
      !> Cycle through every Lebedev degree, then revisit earlier requests
      procedure :: request => repeating_request
   end type repeating_shell_type

contains

   !> Collect all math_grid_atomic tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_atomic(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("shell_sector_edges", test_sector_edges), &
         new_unittest("shell_arc_target_and_bands", test_arc_request), &
         new_unittest("shell_floor_and_cap", test_floor_and_cap), &
         new_unittest("shell_target_fallback", test_target_fallback), &
         new_unittest("atomic_layout", test_atomic_layout), &
         new_unittest("atomic_cache_repeated_requests", test_cache_repeated_requests), &
         new_unittest("atomic_ball_and_shell_volume", test_ball_volume), &
         new_unittest("atomic_radial_integrals", test_radial_integrals), &
         new_unittest("atomic_off_center_gaussian", test_off_center_gaussian), &
         new_unittest("recipe_element_overrides", test_element_overrides) &
         ]
   end subroutine collect_math_grid_atomic

   !* ================================================================================= *!
   !*                                     Integrands                                    *!
   !* ================================================================================= *!

   !> Constant integrand f = 1
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function one_3d(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> 1
      real(wp) :: val

      val = 1.0_wp + 0.0_wp*r(1)
   end function one_3d

   !> Squared distance f = |r|^2
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function r2_3d(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> |r|^2
      real(wp) :: val

      val = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
   end function r2_3d

   !> Cartesian moment f = x^2
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function x2_3d(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> x^2
      real(wp) :: val

      val = r(1)*r(1)
   end function x2_3d

   !> Cartesian moment f = x^2 y^2
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function x2y2_3d(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> x^2 y^2
      real(wp) :: val

      val = r(1)*r(1)*r(2)*r(2)
   end function x2y2_3d

   !> Gaussian exp(-|r|^2/2), integral (2*pi)^(3/2)
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function gauss_half(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val

      val = exp(-0.5_wp*(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
   end function gauss_half

   !> Gaussian exp(-|r|^2), integral pi^(3/2)
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function gauss_one(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val

      val = exp(-(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
   end function gauss_one

   !> Gaussian exp(-4|r|^2), integral (pi/4)^(3/2)
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function gauss_four(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val

      val = exp(-4.0_wp*(r(1)*r(1) + r(2)*r(2) + r(3)*r(3)))
   end function gauss_four

   !> Exponential exp(-|r|), integral 8*pi
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function exp_one(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val

      val = exp(-norm2(r))
   end function exp_one

   !> Exponential exp(-2|r|), integral pi
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function exp_two(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val

      val = exp(-2.0_wp*norm2(r))
   end function exp_two

   !> Off-center Gaussian exp(-|r - d|^2), integral pi^(3/2)
   !>
   !> @param[in] r  Atom-relative point (bohr)
   pure function gauss_off_center(r) result(val)
      !> Atom-relative point (bohr)
      real(wp), intent(in) :: r(3)
      !> Function value
      real(wp) :: val

      real(wp) :: d(3)

      d = r - off_center
      val = exp(-(d(1)*d(1) + d(2)*d(2) + d(3)*d(3)))
   end function gauss_off_center

   !* ================================================================================= *!
   !*                                      Helpers                                      *!
   !* ================================================================================= *!

   !> Relative deviation |a - b|/|b|
   !>
   !> @param[in] a  Value under test
   !> @param[in] b  Reference value, nonzero
   pure function rel_dev(a, b) result(dev)
      !> Value under test
      real(wp), intent(in) :: a
      !> Reference value
      real(wp), intent(in) :: b
      !> Relative deviation
      real(wp) :: dev

      dev = abs(a - b)/abs(b)
   end function rel_dev

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

   !> Assemble an atomic recipe from its parts; every part is copied
   !>
   !> @param[out] recipe      New recipe
   !> @param[in]  kind        Rule selector
   !> @param[in]  npts        Number of radial nodes
   !> @param[in]  mapping     Radial mapping
   !> @param[in]  positive    Lebedev generator with positive weights only
   !> @param[in]  shells      Shell policy
   subroutine assemble_recipe(recipe, kind, npts, mapping, positive, shells)
      !> New recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Rule selector
      integer, intent(in) :: kind
      !> Number of radial nodes
      integer, intent(in) :: npts
      !> Radial mapping
      class(moist_math_grid_radial_mapping_type), intent(in) :: mapping
      !> Lebedev generator with positive weights only
      logical, intent(in) :: positive
      !> Shell policy
      class(moist_math_grid_atomic_shell_type), intent(in) :: shells

      type(moist_math_grid_angular_generator_lebedev_type) :: generator

      call make_rule(kind, recipe%radial%rule)
      recipe%radial%npts = npts
      allocate (recipe%radial%mapping, source=mapping)
      call new_lebedev_generator(generator, positive_weights_only=positive)
      allocate (recipe%angular, source=generator)
      allocate (recipe%shells, source=shells)
   end subroutine assemble_recipe

   !> Recipe with a linear radial mapping on [lower, upper] and Gauss-Legendre nodes
   !>
   !> @param[out] recipe  New recipe
   !> @param[in]  npts    Number of radial nodes
   !> @param[in]  lower   Inner radius (bohr)
   !> @param[in]  upper   Outer radius (bohr)
   !> @param[in]  shells  Shell policy
   !> @param[in]  positive  Lebedev generator with positive weights only
   !> @param[out] merr    Mapping error
   subroutine linear_recipe(recipe, npts, lower, upper, shells, positive, merr)
      !> New recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Number of radial nodes
      integer, intent(in) :: npts
      !> Inner radius (bohr)
      real(wp), intent(in) :: lower
      !> Outer radius (bohr)
      real(wp), intent(in) :: upper
      !> Shell policy
      class(moist_math_grid_atomic_shell_type), intent(in) :: shells
      !> Lebedev generator with positive weights only
      logical, intent(in) :: positive
      !> Mapping error
      type(mctc_error), allocatable, intent(out) :: merr

      type(moist_math_grid_radial_mapping_linear_type) :: lin

      call new_linear_mapping(lin, lower, upper, merr)
      if (allocated(merr)) return
      call assemble_recipe(recipe, rule_gauss_legendre, npts, lin, positive, shells)
   end subroutine linear_recipe

   !> Recipe with a Chebyshev-II rule and an element-scaled Becke mapping
   !>
   !> @param[out] recipe         New recipe
   !> @param[in]  npts           Number of radial nodes
   !> @param[in]  radius_factor  Becke scale as a multiple of the covalent radius
   !> @param[in]  shells         Shell policy
   !> @param[in]  positive       Lebedev generator with positive weights only
   !> @param[out] merr           Mapping error
   subroutine becke_recipe(recipe, npts, radius_factor, shells, positive, merr)
      !> New recipe
      type(moist_math_grid_atomic_recipe_type), intent(out) :: recipe
      !> Number of radial nodes
      integer, intent(in) :: npts
      !> Becke scale as a multiple of the covalent radius
      real(wp), intent(in) :: radius_factor
      !> Shell policy
      class(moist_math_grid_atomic_shell_type), intent(in) :: shells
      !> Lebedev generator with positive weights only
      logical, intent(in) :: positive
      !> Mapping error
      type(mctc_error), allocatable, intent(out) :: merr

      type(moist_math_grid_radial_mapping_becke_type) :: becke

      call new_becke_mapping(becke, merr, radius_factor=radius_factor)
      if (allocated(merr)) return
      call assemble_recipe(recipe, rule_chebyshev2, npts, becke, positive, shells)
   end subroutine becke_recipe

   !> Smallest positive-weight Lebedev size in [nlo, nhi] meeting a target, else the largest
   !>
   !> Independent restatement of the selection rule; 0 if the window is empty
   !>
   !> @param[in] target  Soft target point count
   !> @param[in] nlo     Floor
   !> @param[in] nhi     Cap
   pure function positive_size_for_target(target, nlo, nhi) result(npts)
      !> Soft target point count
      real(wp), intent(in) :: target
      !> Floor
      integer, intent(in) :: nlo
      !> Cap
      integer, intent(in) :: nhi
      !> Selected size, 0 if none
      integer :: npts

      integer :: i

      npts = 0
      do i = 1, size(grid_size)
         if (any(grid_size(i) == lebedev_negative_weight_sizes)) cycle
         if (grid_size(i) < nlo .or. grid_size(i) > nhi) cycle
         npts = grid_size(i)
         if (real(grid_size(i), wp) >= target) return
      end do
   end function positive_size_for_target

   !* ================================================================================= *!
   !*                                   Shell policies                                  *!
   !* ================================================================================= *!

   !> Sector policy: sector of r, shells on an edge in the inner sector
   subroutine test_sector_edges(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_shell_sector_type) :: pol, single
      type(moist_math_grid_angular_request_type) :: req
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: radii(:)
      integer, allocatable :: expected(:)
      character(len=120) :: msg
      integer :: ir

      call new_sector_shell_policy(pol, [1.0_wp, 2.5_wp], [11, 17, 23], merr, &
         & min_points=26, max_points=974)
      call check_moist_error(error, merr, "Sector policy")
      if (allocated(error)) return

      radii = [0.0_wp, 0.5_wp, 1.0_wp, nearest(1.0_wp, 2.0_wp), 2.5_wp, &
         & nearest(2.5_wp, 2.0_wp), 40.0_wp]
      expected = [11, 11, 11, 17, 17, 23, 23]
      do ir = 1, size(radii)
         req = pol%request(8, radii(ir))
         write (msg, "(a,es23.16,a,i0,a,i0)") "Sector policy at r = ", radii(ir), ": degree ", &
            & req%min_degree, ", expected ", expected(ir)
         call check(error, req%min_degree, expected(ir), trim(msg))
         if (allocated(error)) return
         call check(error, req%min_points == 26 .and. req%max_points == 974 .and. req%npts == 0 &
            & .and. req%target_points == 0.0_wp, "Sector policy: hard bounds not copied")
         if (allocated(error)) return
      end do

      ! No edges: one sector covering every radius
      call new_sector_shell_policy(single, [real(wp) ::], [29], merr)
      call check_moist_error(error, merr, "Sector policy without edges")
      if (allocated(error)) return
      do ir = 1, size(radii)
         req = single%request(1, radii(ir))
         call check(error, req%min_degree, 29, "Sector policy without edges: wrong degree")
         if (allocated(error)) return
      end do
   end subroutine test_sector_edges

   !> Arc policy: band spacing, edge membership, target spelling, default bounds
   subroutine test_arc_request(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: edges(2) = [1.0_wp, 3.0_wp]
      real(wp), parameter :: spacing(3) = [0.3_wp, 0.5_wp, 0.9_wp]
      type(moist_math_grid_atomic_shell_arc_type) :: pol, bounded
      type(moist_math_grid_angular_request_type) :: req
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: radii(:), band(:)
      real(wp) :: h
      character(len=120) :: msg
      integer :: ir

      call new_arc_shell_policy(pol, edges, spacing, merr)
      call check_moist_error(error, merr, "Arc policy")
      if (allocated(error)) return
      call check(error, pol%min_points == 110 .and. pol%max_points == 5810, &
         & "Arc policy: default floor and cap are not 110 and 5810")
      if (allocated(error)) return

      radii = [0.0_wp, 0.2_wp, 1.0_wp, nearest(1.0_wp, 2.0_wp), 2.0_wp, 3.0_wp, &
         & nearest(3.0_wp, 2.0_wp), 12.0_wp]
      band = [0.3_wp, 0.3_wp, 0.3_wp, 0.5_wp, 0.5_wp, 0.5_wp, 0.9_wp, 0.9_wp]
      do ir = 1, size(radii)
         req = pol%request(6, radii(ir))
         h = band(ir)
         write (msg, "(a,es23.16,a,es23.16)") "Arc policy at r = ", radii(ir), ": target ", &
            & req%target_points
         call check(error, req%target_points == 4.0_wp*pi*radii(ir)*radii(ir)/(h*h), trim(msg))
         if (allocated(error)) return
         call check(error, req%min_degree == 0 .and. req%npts == 0 .and. req%min_points == 110 &
            & .and. req%max_points == 5810, "Arc policy: request carries more than the target and bounds")
         if (allocated(error)) return
      end do

      call new_arc_shell_policy(bounded, edges, spacing, merr, min_points=26, max_points=302)
      call check_moist_error(error, merr, "Arc policy with bounds")
      if (allocated(error)) return
      req = bounded%request(1, 2.0_wp)
      call check(error, req%min_points == 26 .and. req%max_points == 302, &
         & "Arc policy: explicit floor and cap not copied")
   end subroutine test_arc_request

   !> Hard floors lift and caps bound the selected sizes of every policy
   subroutine test_floor_and_cap(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_shell_constant_type) :: cpol
      type(moist_math_grid_atomic_shell_sector_type) :: spol
      type(moist_math_grid_atomic_shell_arc_type) :: apol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: ish

      ! Degree 3 alone selects 6 points; the floor lifts it to 110
      call new_constant_shell_policy(cpol, 3, merr, min_points=110)
      call check_moist_error(error, merr, "Constant policy with floor")
      if (allocated(error)) return
      call linear_recipe(recipe, 4, 0.0_wp, 3.0_wp, cpol, .true., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 6, merr)
      call check_moist_error(error, merr, "Constant policy with floor")
      if (allocated(error)) return
      call check(error, all(grid%shell_npts == 110) .and. all(grid%shell_degree == 17), &
         & "Constant policy: floor 110 not applied")
      if (allocated(error)) return

      ! Inner sector degree 5 lifted to the 26-point floor, outer sector degree 41 (590)
      call new_sector_shell_policy(spol, [1.5_wp], [5, 41], merr, min_points=26)
      call check_moist_error(error, merr, "Sector policy with floor")
      if (allocated(error)) return
      call linear_recipe(recipe, 6, 0.0_wp, 3.0_wp, spol, .true., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 6, merr)
      call check_moist_error(error, merr, "Sector policy with floor")
      if (allocated(error)) return
      do ish = 1, grid%nshell
         if (grid%shell_r(ish) <= 1.5_wp) then
            call check(error, grid%shell_npts(ish), 26, "Sector policy: inner shell not at its floor")
         else
            call check(error, grid%shell_npts(ish), 590, "Sector policy: outer shell not at degree 41")
         end if
         if (allocated(error)) return
      end do

      ! Spacing far below reach: every shell outside the floor capped at 434
      call new_arc_shell_policy(apol, [real(wp) ::], [1.0e-3_wp], merr, min_points=6, max_points=434)
      call check_moist_error(error, merr, "Arc policy with cap")
      if (allocated(error)) return
      call linear_recipe(recipe, 5, 1.0_wp, 3.0_wp, apol, .true., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 6, merr)
      call check_moist_error(error, merr, "Arc policy with cap")
      if (allocated(error)) return
      call check(error, all(grid%shell_npts == 434), "Arc policy: cap 434 not applied")
   end subroutine test_floor_and_cap

   !> Unreachable targets fall back to the largest admissible rule
   subroutine test_target_fallback(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: h = 0.1_wp
      integer, parameter :: nlo = 110, nhi = 302
      type(moist_math_grid_atomic_shell_arc_type) :: apol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: r, target
      character(len=160) :: msg
      integer :: ish, expected, nshort, nmet

      call new_arc_shell_policy(apol, [real(wp) ::], [h], merr, min_points=nlo, max_points=nhi)
      call check_moist_error(error, merr, "Arc policy")
      if (allocated(error)) return
      call linear_recipe(recipe, 12, 0.0_wp, 6.0_wp, apol, .true., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 6, merr)
      call check_moist_error(error, merr, "Arc policy atomic grid")
      if (allocated(error)) return

      nshort = 0
      nmet = 0
      do ish = 1, grid%nshell
         r = grid%shell_r(ish)
         target = 4.0_wp*pi*r*r/(h*h)
         expected = positive_size_for_target(target, nlo, nhi)
         write (msg, "(a,es12.5,a,i0,a,i0)") "Arc fallback at r = ", r, ": selected ", &
            & grid%shell_npts(ish), ", expected ", expected
         call check(error, grid%shell_npts(ish), expected, trim(msg))
         if (allocated(error)) return

         if (target > real(nhi, wp)) then
            nshort = nshort + 1
         else
            nmet = nmet + 1
         end if
      end do
      call check(error, nshort > 0 .and. nmet > 0, "Arc fallback: case does not cover both regimes")
   end subroutine test_target_fallback

   !* ================================================================================= *!
   !*                                 Atomic composition                                *!
   !* ================================================================================= *!

   !> Point order, shell bookkeeping, and the separate radius and direction
   !>
   !> - Shells in the radial grid's order, angular points in table order
   !> - xyz is exactly r*u; each shell's weights sum to 4*pi*r^2*w_rad
   subroutine test_atomic_layout(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: z = 6, nrad = 20
      type(moist_math_grid_atomic_shell_sector_type) :: spol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: grid
      type(moist_math_grid_radial_type) :: radial
      type(moist_math_grid_angular_type) :: s2
      type(mctc_error), allocatable :: merr
      real(wp) :: r, shell_sum
      integer :: ish, j, ip, i0

      call new_sector_shell_policy(spol, [1.0_wp], [7, 17], merr)
      call check_moist_error(error, merr, "Sector policy")
      if (allocated(error)) return
      call becke_recipe(recipe, nrad, 0.5_wp, spol, .true., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, z, merr)
      call check_moist_error(error, merr, "Layout grid")
      if (allocated(error)) return
      call new_radial_grid(radial, recipe%radial, z, merr)
      call check_moist_error(error, merr, "Layout radial grid")
      if (allocated(error)) return

      call check(error, grid%z == z .and. grid%nshell == nrad .and. grid%nshell_requested == nrad, &
         & "Layout: element or shell counts not recorded")
      if (allocated(error)) return
      call check(error, all(grid%shell_r == radial%r) .and. all(grid%shell_w == radial%w), &
         & "Layout: shells not in the radial grid's order")
      if (allocated(error)) return
      call check(error, any(grid%shell_npts == 26) .and. any(grid%shell_npts == 110), &
         & "Layout: sector policy did not produce two angular sizes")
      if (allocated(error)) return
      call check(error, grid%npts, sum(grid%shell_npts), "Layout: npts is not the sum of shell sizes")
      if (allocated(error)) return
      call check(error, grid%shell_offset(1) == 1 .and. grid%shell_offset(nrad + 1) == grid%npts + 1, &
         & "Layout: shell offsets do not span the points")
      if (allocated(error)) return

      ip = 0
      do ish = 1, nrad
         i0 = grid%shell_offset(ish)
         call check(error, grid%shell_offset(ish + 1) - i0, grid%shell_npts(ish), &
            & "Layout: shell offsets disagree with shell sizes")
         if (allocated(error)) return
         call new_lebedev_grid(s2, merr, npts=grid%shell_npts(ish))
         call check_moist_error(error, merr, "Layout angular rule")
         if (allocated(error)) return
         call check(error, grid%shell_degree(ish), s2%degree, "Layout: shell degree not recorded")
         if (allocated(error)) return

         r = grid%shell_r(ish)
         shell_sum = 0.0_wp
         do j = 1, s2%npts
            ip = ip + 1
            call check(error, ip == i0 + j - 1 .and. grid%shell(ip) == ish, &
               & "Layout: point-to-shell index out of order")
            if (allocated(error)) return
            call check(error, all(grid%u(:, ip) == s2%points(:, j)), &
               & "Layout: directions not in table order")
            if (allocated(error)) return
            call check(error, all(grid%xyz(:, ip) == r*grid%u(:, ip)), &
               & "Layout: local point is not r*u")
            if (allocated(error)) return
            shell_sum = shell_sum + grid%w(ip)
         end do
         call check(error, rel_dev(shell_sum, 4.0_wp*pi*r*r*grid%shell_w(ish)) <= 1.0e-14_wp, &
            & "Layout: shell weights do not sum to 4*pi*r^2*w_rad")
         if (allocated(error)) return
      end do
   end subroutine test_atomic_layout

   !> Cycle through the Lebedev degrees in consecutive unit-radius bands
   !>
   !> @param[in] self  Shell policy
   !> @param[in] z     Atomic number
   !> @param[in] r     Shell radius (bohr)
   pure function repeating_request(self, z, r) result(request)
      !> Shell policy
      class(repeating_shell_type), intent(in) :: self
      !> Atomic number
      integer, intent(in) :: z
      !> Shell radius (bohr)
      real(wp), intent(in) :: r
      !> Angular request
      type(moist_math_grid_angular_request_type) :: request

      request%min_degree = lebedev_degree_table(mod(int(r), size(lebedev_degree_table)) + 1) + 0*z
      request%min_points = self%min_points
      request%max_points = self%max_points
   end function repeating_request

   !> Repeated requests retain their angular rules across cache growth
   !>
   !> - One shell per Lebedev degree (32 distinct requests and rules), then
   !>   every request again; the distinct-request and distinct-rule stores of
   !>   the atomic grid start at 8 entries, so they must grow, and the second
   !>   pass reads entries the growth copied
   !> - Covers initial store capacities up to 31; a larger start skips growth
   subroutine test_cache_repeated_requests(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(repeating_shell_type) :: shells
      type(moist_math_grid_radial_mapping_linear_type) :: mapping
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: grid
      type(moist_math_grid_angular_type) :: angular
      type(mctc_error), allocatable :: merr
      character(len=120) :: msg
      integer :: ish, nshell, degree

      ! Midpoint shells at r = ish - 0.5, so int(r) = ish - 1 selects the degree
      nshell = 2*size(lebedev_degree_table)
      call new_linear_mapping(mapping, 0.0_wp, real(nshell, wp), merr)
      call check_moist_error(error, merr, "Repeating request mapping")
      if (allocated(error)) return
      call assemble_recipe(recipe, rule_midpoint, nshell, mapping, .false., shells)
      call new_atomic_grid(grid, recipe, 6, merr)
      call check_moist_error(error, merr, "Repeating request grid")
      if (allocated(error)) return
      write (msg, "(a,i0,a,i0)") "Repeating request grid: nshell ", grid%nshell, ", expected ", nshell
      call check(error, grid%nshell == nshell, trim(msg))
      if (allocated(error)) return
      do ish = 1, nshell
         degree = lebedev_degree_table(mod(ish - 1, size(lebedev_degree_table)) + 1)
         call new_lebedev_grid(angular, merr, degree=degree)
         call check_moist_error(error, merr, "Repeating request reference")
         if (allocated(error)) return
         write (msg, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "Shell ", ish, ": degree ", &
            & grid%shell_degree(ish), ", expected ", angular%degree, "; npts ", &
            & grid%shell_npts(ish), ", expected ", angular%npts
         call check(error, grid%shell_degree(ish) == angular%degree &
            & .and. grid%shell_npts(ish) == angular%npts, &
            & "Repeated request lost its angular rule after cache growth: "//trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_cache_repeated_requests

   !> Ball and shell volumes and exact moments with a linear Gauss-Legendre radial grid
   !>
   !> A missing or duplicated r^2 or 4*pi changes every value by orders of magnitude
   subroutine test_ball_volume(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: rball = 2.5_wp, rinner = 1.0_wp, tol = 1.0e-14_wp
      type(moist_math_grid_atomic_shell_constant_type) :: cpol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: val, ref
      character(len=120) :: msg

      ! 4 Gauss-Legendre nodes: exact through r^7; degree 5 (14 points): exact through x^2 y^2
      call new_constant_shell_policy(cpol, 5, merr)
      call check_moist_error(error, merr, "Constant policy")
      if (allocated(error)) return
      call linear_recipe(recipe, 4, 0.0_wp, rball, cpol, .false., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 8, merr)
      call check_moist_error(error, merr, "Ball grid")
      if (allocated(error)) return

      call grid%integrate(one_3d, val)
      ref = 4.0_wp/3.0_wp*pi*rball**3
      write (msg, "(a,es23.16,a,es23.16)") "Ball volume ", val, ", expected ", ref
      call check(error, rel_dev(val, ref) <= tol, trim(msg))
      if (allocated(error)) return
      call grid%integrate(r2_3d, val)
      ref = 4.0_wp*pi*rball**5/5.0_wp
      write (msg, "(a,es23.16,a,es23.16)") "Ball |r|^2 moment ", val, ", expected ", ref
      call check(error, rel_dev(val, ref) <= tol, trim(msg))
      if (allocated(error)) return
      call grid%integrate(x2_3d, val)
      ref = 4.0_wp*pi*rball**5/15.0_wp
      write (msg, "(a,es23.16,a,es23.16)") "Ball x^2 moment ", val, ", expected ", ref
      call check(error, rel_dev(val, ref) <= tol, trim(msg))
      if (allocated(error)) return
      call grid%integrate(x2y2_3d, val)
      ref = 4.0_wp*pi*rball**7/105.0_wp
      write (msg, "(a,es23.16,a,es23.16)") "Ball x^2 y^2 moment ", val, ", expected ", ref
      call check(error, rel_dev(val, ref) <= tol, trim(msg))
      if (allocated(error)) return

      call linear_recipe(recipe, 4, rinner, rball, cpol, .false., merr)
      if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 8, merr)
      call check_moist_error(error, merr, "Shell grid")
      if (allocated(error)) return
      call grid%integrate(one_3d, val)
      ref = 4.0_wp/3.0_wp*pi*(rball**3 - rinner**3)
      write (msg, "(a,es23.16,a,es23.16)") "Spherical shell volume ", val, ", expected ", ref
      call check(error, rel_dev(val, ref) <= tol, trim(msg))
   end subroutine test_ball_volume

   !> Gaussian and exponential integrals across elements and radial scales
   !>
   !> - Per-element default recipes for H, C, Na, and Br (Chebyshev-II + Becke)
   !> - Fixed Becke scales 1 and 3 bohr, and midpoint + Knowles (k = 3, R = 5)
   subroutine test_radial_integrals(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: elements(4) = [1, 6, 11, 35]
      !> Relative quadrature error per element; measured 4.7e-9, 4.1e-12, 4.3e-11, 3.2e-12
      real(wp), parameter :: element_tol(4) = [2.0e-8_wp, 2.0e-11_wp, 2.0e-10_wp, 2.0e-11_wp]
      !> Relative quadrature error per scale case; measured 1.1e-11, 3.0e-10, 3.5e-7
      real(wp), parameter :: scale_tol(3) = [5.0e-11_wp, 2.0e-9_wp, 2.0e-6_wp]
      type(moist_math_grid_atomic_recipe_type) :: defaults, recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(moist_math_grid_atomic_shell_constant_type) :: cpol
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_mapping_knowles_type) :: knowles
      type(moist_math_grid_atomic_type) :: grid
      type(mctc_error), allocatable :: merr
      character(len=64) :: label
      integer :: iz, icase

      call default_element_recipes(defaults, overrides, merr)
      call check_moist_error(error, merr, "Default recipes")
      if (allocated(error)) return
      do iz = 1, size(elements)
         call get_element_recipe(recipe, defaults, elements(iz), overrides)
         call new_atomic_grid(grid, recipe, elements(iz), merr)
         write (label, "(a,i0)") "Default recipe, z = ", elements(iz)
         call check_moist_error(error, merr, trim(label))
         if (allocated(error)) return
         call check_radial_integrals(error, grid, trim(label), element_tol(iz))
         if (allocated(error)) return
      end do

      call new_constant_shell_policy(cpol, 11, merr)
      call check_moist_error(error, merr, "Constant policy")
      if (allocated(error)) return
      do icase = 1, 3
         select case (icase)
         case (1)
            call new_becke_mapping(becke, merr, scale=1.0_wp)
            if (.not. allocated(merr)) call assemble_recipe(recipe, rule_chebyshev2, 75, becke, .true., cpol)
            label = "Chebyshev-II + Becke, scale 1"
         case (2)
            call new_becke_mapping(becke, merr, scale=3.0_wp)
            if (.not. allocated(merr)) call assemble_recipe(recipe, rule_chebyshev2, 75, becke, .true., cpol)
            label = "Chebyshev-II + Becke, scale 3"
         case default
            call new_knowles_mapping(knowles, 3.0_wp, merr, scale=5.0_wp)
            if (.not. allocated(merr)) call assemble_recipe(recipe, rule_midpoint, 75, knowles, .true., cpol)
            label = "Midpoint + Knowles, R = 5"
         end select
         if (.not. allocated(merr)) call new_atomic_grid(grid, recipe, 8, merr)
         call check_moist_error(error, merr, trim(label))
         if (allocated(error)) return
         call check_radial_integrals(error, grid, trim(label), scale_tol(icase))
         if (allocated(error)) return
      end do
   end subroutine test_radial_integrals

   !> Check three Gaussian and two exponential integrals on one grid
   !> @param[in]  grid   Atomic grid
   !> @param[in]  label  Case description
   !> @param[in]  tol    Relative tolerance
   subroutine check_radial_integrals(error, grid, label, tol)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Atomic grid
      type(moist_math_grid_atomic_type), intent(in) :: grid
      !> Case description
      character(len=*), intent(in) :: label
      !> Relative tolerance
      real(wp), intent(in) :: tol

      real(wp) :: val(5), ref(5), dev
      character(len=160) :: msg
      integer :: i

      call grid%integrate(gauss_half, val(1))
      call grid%integrate(gauss_one, val(2))
      call grid%integrate(gauss_four, val(3))
      call grid%integrate(exp_one, val(4))
      call grid%integrate(exp_two, val(5))
      ref = [(2.0_wp*pi)**1.5_wp, pi**1.5_wp, (0.25_wp*pi)**1.5_wp, 8.0_wp*pi, pi]
      dev = 0.0_wp
      do i = 1, 5
         dev = max(dev, rel_dev(val(i), ref(i)))
      end do
      write (msg, "(a,a,es10.3)") label, ": largest relative integral error ", dev
      call check(error, dev <= tol, trim(msg))
   end subroutine check_radial_integrals

   !> Off-center Gaussian with the carbon default recipe
   !>
   !> - Couples radial and angular quadrature: exp(-|r - d|^2) with |d| = 0.54 bohr
   !> - Same accuracy as a callback and as samples tabulated on the points
   !> - A field of the wrong length is rejected, and the next valid call succeeds
   subroutine test_off_center_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Relative quadrature error; measured 3.9e-13
      real(wp), parameter :: tol = 2.0e-12_wp
      type(moist_math_grid_atomic_recipe_type) :: defaults, recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      type(moist_math_grid_atomic_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: samples(:)
      real(wp) :: val, dev
      integer :: i
      character(len=120) :: msg

      call default_element_recipes(defaults, overrides, merr)
      call check_moist_error(error, merr, "Default recipes")
      if (allocated(error)) return
      call get_element_recipe(recipe, defaults, 6, overrides)
      call new_atomic_grid(grid, recipe, 6, merr)
      call check_moist_error(error, merr, "Carbon grid")
      if (allocated(error)) return
      call grid%integrate(gauss_off_center, val)
      dev = rel_dev(val, pi**1.5_wp)
      write (msg, "(a,es10.3)") "Off-center Gaussian: relative error ", dev
      call check(error, dev <= tol, trim(msg))
      if (allocated(error)) return

      allocate (samples(grid%npts))
      do i = 1, grid%npts
         samples(i) = gauss_off_center(grid%xyz(:, i))
      end do
      ! A short field is rejected; the valid call below still succeeds
      call grid%integrate_field(samples(2:), val, merr)
      call check(error, allocated(merr), "integrate_field accepted a short field")
      if (allocated(error)) return
      call grid%integrate_field(samples, val, merr)
      call check_moist_error(error, merr, "Off-center Gaussian samples")
      if (allocated(error)) return
      dev = rel_dev(val, pi**1.5_wp)
      write (msg, "(a,es10.3)") "Off-center Gaussian samples: relative error ", dev
      call check(error, dev <= tol, trim(msg))
   end subroutine test_off_center_gaussian

   !* ================================================================================= *!
   !*                                      Recipes                                      *!
   !* ================================================================================= *!

   !> Override lookup: first listing wins, default otherwise, unallocated list is absent
   subroutine test_element_overrides(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_atomic_shell_constant_type) :: cpol
      type(moist_math_grid_atomic_recipe_type) :: default_recipe, selected
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:), none(:)
      type(moist_math_grid_atomic_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_constant_shell_policy(cpol, 11, merr)
      call check_moist_error(error, merr, "Constant policy")
      if (allocated(error)) return
      call becke_recipe(default_recipe, 14, 0.5_wp, cpol, .true., merr)
      call check_moist_error(error, merr, "Default recipe")
      if (allocated(error)) return
      allocate (overrides(2))
      overrides(1)%elements = [6, 7]
      call becke_recipe(overrides(1)%recipe, 10, 0.5_wp, cpol, .true., merr)
      call check_moist_error(error, merr, "First override")
      if (allocated(error)) return
      overrides(2)%elements = [7, 8]
      call becke_recipe(overrides(2)%recipe, 12, 0.5_wp, cpol, .true., merr)
      call check_moist_error(error, merr, "Second override")
      if (allocated(error)) return

      call check(error, element_override_index(6, overrides) == 1 &
         & .and. element_override_index(7, overrides) == 1 &
         & .and. element_override_index(8, overrides) == 2 &
         & .and. element_override_index(1, overrides) == 0, "Overrides: wrong override index")
      if (allocated(error)) return
      call check(error, element_override_index(7, none) == 0 .and. element_override_index(7) == 0, &
         & "Overrides: unallocated or absent list does not select the default")
      if (allocated(error)) return
      allocate (none(0))
      call check(error, element_override_index(7, none) == 0, "Overrides: empty list selects an override")
      if (allocated(error)) return

      call get_element_recipe(selected, default_recipe, 7, overrides)
      call check(error, selected%radial%npts, 10, "Overrides: first listing does not win")
      if (allocated(error)) return
      call get_element_recipe(selected, default_recipe, 8, overrides)
      call check(error, selected%radial%npts, 12, "Overrides: second override not selected")
      if (allocated(error)) return
      call get_element_recipe(selected, default_recipe, 1, overrides)
      call check(error, selected%radial%npts, 14, "Overrides: default not selected")
      if (allocated(error)) return
      deallocate (none)
      call get_element_recipe(selected, default_recipe, 7, none)
      call check(error, selected%radial%npts, 14, "Overrides: unallocated list does not select the default")
      if (allocated(error)) return

      call get_element_recipe(selected, default_recipe, 8, overrides)
      call new_atomic_grid(grid, selected, 8, merr)
      call check_moist_error(error, merr, "Override grid")
      if (allocated(error)) return
      call check(error, grid%nshell_requested == 12 .and. grid%npts == 12*50, &
         & "Overrides: atomic grid not built from the override")
   end subroutine test_element_overrides

end module test_math_grid_atomic
