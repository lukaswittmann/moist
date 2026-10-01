!> Test suite for the angular grid, the Lebedev generator, and the Lebedev tables
!>
!> - Angular grid: unit directions, weights summing to 4*pi, monomial moments,
!>   achieved degree checked with zonal harmonics on all 32 rules,
!>   `integrate` against `integrate_field`, idempotent destroy
!> - `new_lebedev_grid`: node and weight structure, selector validation,
!>   degree selection, low-order Cartesian moments on a spread of rules, and
!>   every one of the 32 rules at its full algebraic degree; these checks
!>   divide the solid-angle integrals by 4*pi and compare sphere averages
!> - Lebedev metadata: negative-weight list against the generated tables,
!>   `lebedev_order_from_num` filter and error messages, raw tables with
!>   weights summing to 1 and points on the unit sphere
!> - Selection: all-weights default, positive-only opt-in, exact requests,
!>   minimum degree, target with floors and caps, incompatible constraints
module test_math_grid_angular
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_grid_angular_lebedev, only: grid_size, lebedev_degree_table, &
      & lebedev_negative_weight_sizes, lebedev_has_negative_weights, lebedev_order_from_num, &
      & get_angular_grid
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, integrand_angular, &
      & moist_math_grid_angular_generator_type, moist_math_grid_angular_request_type, &
      & moist_math_grid_angular_generator_lebedev_type, new_lebedev_generator, new_lebedev_grid
   implicit none(type, external)
   private

   public :: collect_math_grid_angular

   !> Point counts of a representative spread of supported Lebedev rules
   integer, parameter :: npts_list(4) = [6, 74, 434, 590]
   !> Tolerance for the structural (sum-to-one, unit-vector) checks
   real(wp), parameter :: struct_tol = 1.0e-14_wp
   !> Tolerance for the quadrature exactness checks
   real(wp), parameter :: quad_tol = 1.0e-12_wp
   !> Number of supported Lebedev rules
   integer, parameter :: nrules = 32
   !> Point counts of a representative spread of rules, with and without negative weights
   integer, parameter :: npts_list_wide(7) = [6, 26, 74, 110, 434, 590, 5810]
   !> Tolerance of the unit-vector and integrator-agreement checks
   real(wp), parameter :: unit_tol = 1.0e-14_wp
   !> Relative tolerance of the total weight; measured at most 2.0e-14 (3890 points)
   real(wp), parameter :: sum_tol = 5.0e-14_wp
   !> Bound on a zonal moment the rule integrates exactly; measured at most 3.9e-15
   real(wp), parameter :: exact_tol = 1.0e-13_wp
   !> Lower bound on the zonal error one degree above the rule; measured at least 9.6e-3
   real(wp), parameter :: inexact_margin = 1.0e-3_wp
   !> Bound on a monomial moment the rule integrates exactly; measured at most 5.9e-14
   real(wp), parameter :: mono_exact_tol = 1.0e-12_wp
   !> Lower bound on the worst monomial error one degree above the rule; measured at least 4.2e-8
   real(wp), parameter :: mono_inexact_margin = 1.0e-9_wp
   !> Message of an exact request for the 74-point rule with the filter on
   character(len=*), parameter :: msg_neg74 = "Lebedev size 74 has negative weights and is "// &
      & "not allowed here"

contains

   !> Collect all math_grid_angular tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_angular(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("angular_grid_unit_directions", test_unit_directions), &
         new_unittest("angular_grid_weights_sum_to_four_pi", test_weights_four_pi), &
         new_unittest("angular_grid_monomial_moments", test_monomial_moments), &
         new_unittest("angular_grid_achieved_degree", test_achieved_degree), &
         new_unittest("angular_grid_integrate_matches_field", test_integrate_matches_field), &
         new_unittest("angular_grid_destroy_idempotent", test_destroy), &
         new_unittest("angular_lebedev_weights_sum_to_four_pi", test_weights_sum), &
         new_unittest("angular_lebedev_points_are_unit_vectors", test_unit_vectors), &
         new_unittest("angular_lebedev_requires_one_selector", test_selector_validation), &
         new_unittest("angular_lebedev_rejects_unsupported_npts", test_unsupported_npts), &
         new_unittest("angular_lebedev_degree_selection", test_degree_selection), &
         new_unittest("angular_lebedev_degree_out_of_range", test_degree_out_of_range), &
         new_unittest("angular_lebedev_polynomial_exactness", test_polynomial_exactness), &
         new_unittest("angular_lebedev_exact_at_full_degree", test_full_degree_exactness), &
         new_unittest("angular_integrate_matches_integrate_field", test_integrate_agreement), &
         new_unittest("angular_lebedev_destroy_idempotent", test_destroy_idempotent), &
         new_unittest("lebedev_negative_weight_list", test_negative_weight_list), &
         new_unittest("lebedev_order_from_num_filter", test_order_from_num_filter), &
         new_unittest("angular_weights_sum_to_one", test_angular_weights_sum), &
         new_unittest("angular_points_on_unit_sphere", test_angular_unit_sphere), &
         new_unittest("lebedev_default_keeps_all_weights", test_default_all_weights), &
         new_unittest("lebedev_positive_only_opt_in", test_positive_only), &
         new_unittest("lebedev_exact_request_errors", test_exact_request_errors), &
         new_unittest("lebedev_min_degree_selection", test_min_degree_selection), &
         new_unittest("lebedev_target_selection", test_target_selection), &
         new_unittest("lebedev_incompatible_constraints", test_incompatible_constraints), &
         new_unittest("lebedev_new_lebedev_grid", test_new_lebedev_grid), &
         new_unittest("lebedev_polymorphic_generator", test_polymorphic_generator) &
         ]
   end subroutine collect_math_grid_angular

   !* ================================================================================= *!
   !*                                      Helpers                                      *!
   !* ================================================================================= *!

   !> Run a selection that must succeed and check the selected size
   !>
   !> Also checks that the grid is consistently filled for that size
   !>
   !> @param[out] error      Test failure
   !> @param[in]  generator  Generator under test
   !> @param[in]  request    Angular request
   !> @param[in]  expected   Expected point count
   !> @param[in]  label      Case description for the failure message
   subroutine expect_npts(error, generator, request, expected, label)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Generator under test
      class(moist_math_grid_angular_generator_type), intent(in) :: generator
      !> Angular request
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Expected point count
      integer, intent(in) :: expected
      !> Case description
      character(len=*), intent(in) :: label

      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      character(len=200) :: msg
      integer :: order

      call generator%select(request, grid, merr)
      if (allocated(merr)) then
         call test_failed(error, label//": "//merr%message)
         return
      end if
      write (msg, "(a,a,i0,a,i0)") label, ": selected ", grid%npts, " points, expected ", expected
      call check(error, grid%npts, expected, trim(msg))
      if (allocated(error)) return

      order = findloc(grid_size, expected, dim=1)
      call check(error, grid%degree, lebedev_degree_table(order), label//": wrong achieved degree")
      if (allocated(error)) return
      call check(error, size(grid%points, 2) == expected .and. size(grid%weights) == expected, &
         & label//": grid arrays do not match npts")
   end subroutine expect_npts

   !> Run a selection that must fail and leave the grid empty
   !>
   !> @param[out] error      Test failure
   !> @param[in]  generator  Generator under test
   !> @param[in]  request    Angular request
   !> @param[in]  label      Case description for the failure message
   !> @param[in]  message    Exact expected error message, optional
   subroutine expect_failure(error, generator, request, label, message)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Generator under test
      class(moist_math_grid_angular_generator_type), intent(in) :: generator
      !> Angular request
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Case description
      character(len=*), intent(in) :: label
      !> Exact expected error message
      character(len=*), intent(in), optional :: message

      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      call generator%select(request, grid, merr)
      call check(error, allocated(merr), label//": request was accepted")
      if (allocated(error)) return
      call check(error, grid%npts == 0 .and. .not. allocated(grid%points) &
         & .and. .not. allocated(grid%weights), label//": rejected grid is not empty")
      if (allocated(error)) return
      if (present(message)) then
         call check(error, merr%message == message, label//": unexpected message: "//merr%message)
      end if
   end subroutine expect_failure

   !> Build the grid of an exact point count with the default generator
   !>
   !> @param[out] error  Test failure
   !> @param[in]  npts   Supported point count
   !> @param[out] grid   Generated grid
   subroutine build_exact(error, npts, grid)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Supported point count
      integer, intent(in) :: npts
      !> Generated grid
      type(moist_math_grid_angular_type), intent(out) :: grid

      type(mctc_error), allocatable :: merr

      call new_lebedev_grid(grid, merr, npts=npts)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine build_exact

   !> Exact integral of the monomial x**a*y**b*z**c over the unit sphere
   !>
   !> `2*Gamma(A)*Gamma(B)*Gamma(C)/Gamma(A+B+C)` with `A = (a+1)/2` etc.,
   !> zero if any exponent is odd
   !>
   !> @param[in] a  Exponent of x
   !> @param[in] b  Exponent of y
   !> @param[in] c  Exponent of z
   pure function monomial_integral(a, b, c) result(r)
      !> Exponent of x
      integer, intent(in) :: a
      !> Exponent of y
      integer, intent(in) :: b
      !> Exponent of z
      integer, intent(in) :: c
      !> Integral over the unit sphere
      real(wp) :: r

      real(wp) :: ha, hb, hc

      if (mod(a, 2) /= 0 .or. mod(b, 2) /= 0 .or. mod(c, 2) /= 0) then
         r = 0.0_wp
         return
      end if
      ha = 0.5_wp*real(a + 1, wp)
      hb = 0.5_wp*real(b + 1, wp)
      hc = 0.5_wp*real(c + 1, wp)
      r = 2.0_wp*gamma(ha)*gamma(hb)*gamma(hc)/gamma(ha + hb + hc)
   end function monomial_integral

   !> Compare the grid average of an integrand against its exact value
   !>
   !> The average is the solid-angle integral divided by 4*pi
   !>
   !> @param[in,out] error     Test error handling
   !> @param[in]     grid      Grid instance under test
   !> @param[in]     f         Angular integrand
   !> @param[in]     expected  Exact sphere average of f
   !> @param[in]     label     Name of the integrand for the failure message
   subroutine check_average(error, grid, f, expected, label)
      !> Test error handling
      type(error_type), allocatable, intent(inout) :: error
      !> Grid instance under test
      type(moist_math_grid_angular_type), intent(in) :: grid
      !> Angular integrand
      procedure(integrand_angular) :: f
      !> Exact sphere average of f
      real(wp), intent(in) :: expected
      !> Name of the integrand
      character(len=*), intent(in) :: label

      real(wp) :: val

      call grid%integrate(f, val)
      val = val/(4.0_wp*pi)
      call check(error, abs(val - expected) < quad_tol, &
         & "Angular quadrature is not exact for "//label)
   end subroutine check_average

   !> Constant integrand f(u) = 1
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_one(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = 1.0_wp + 0.0_wp*uvec(1)
   end function f_one

   !> Degree-1 integrand f(u) = x
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_x(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = uvec(1)
   end function f_x

   !> Degree-2 integrand f(u) = x*y
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_xy(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = uvec(1)*uvec(2)
   end function f_xy

   !> Degree-2 integrand f(u) = 3*z^2 - 1 (proportional to Y_20)
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_3z2m1(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = 3.0_wp*uvec(3)**2 - 1.0_wp
   end function f_3z2m1

   !> Degree-2 integrand f(u) = x^2, average 1/3
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_x2(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = uvec(1)**2
   end function f_x2

   !> Degree-4 integrand f(u) = x^4 + y^4 + z^4, average 3/5
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_quartic(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = uvec(1)**4 + uvec(2)**4 + uvec(3)**4
   end function f_quartic

   !> Mixed-parity integrand used to cross-check the two quadrature entries
   !>
   !> @param[in]  uvec  Cartesian unit vector
   pure function f_mixed(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = exp(uvec(1))*(1.0_wp + 2.0_wp*uvec(2)*uvec(3)) + uvec(3)**3
   end function f_mixed

   !> Smooth mixed-parity integrand
   !>
   !> @param[in] uvec  Unit vector
   pure function f_smooth(uvec) result(val)
      !> Unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = exp(0.3_wp*uvec(1) - 0.2_wp*uvec(2))*(1.0_wp + uvec(3)*uvec(3))
   end function f_smooth

   !* ================================================================================= *!
   !*                                   Angular grid                                    *!
   !* ================================================================================= *!

   !> Nodes are Cartesian unit vectors and arrays match npts
   !>
   !> @param[out] error  Test failure
   subroutine test_unit_directions(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      integer :: k, i

      do k = 1, size(npts_list_wide)
         call build_exact(error, npts_list_wide(k), grid)
         if (allocated(error)) return
         call check(error, grid%npts, npts_list_wide(k), "Unexpected number of angular nodes")
         if (allocated(error)) return
         call check(error, size(grid%points, 1) == 3 .and. size(grid%points, 2) == grid%npts &
            & .and. size(grid%weights) == grid%npts, "Angular grid arrays do not match npts")
         if (allocated(error)) return
         do i = 1, grid%npts
            call check(error, abs(norm2(grid%points(:, i)) - 1.0_wp) < unit_tol, &
               & "Angular node is not a unit vector")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_unit_directions

   !> Weights sum to 4*pi, directly and through both integrators
   !>
   !> @param[out] error  Test failure
   subroutine test_weights_four_pi(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      real(wp), allocatable :: ones(:)
      real(wp) :: total
      integer :: k

      do k = 1, size(npts_list_wide)
         call build_exact(error, npts_list_wide(k), grid)
         if (allocated(error)) return

         call check(error, abs(sum(grid%weights) - 4.0_wp*pi) <= sum_tol*4.0_wp*pi, &
            & "Angular weights do not sum to 4*pi")
         if (allocated(error)) return

         allocate (ones(grid%npts), source=1.0_wp)
         call grid%integrate_field(ones, total)
         deallocate (ones)
         call check(error, abs(total - 4.0_wp*pi) <= sum_tol*4.0_wp*pi, &
            & "integrate_field of a unit field is not 4*pi")
         if (allocated(error)) return

         call grid%integrate(f_one, total)
         call check(error, abs(total - 4.0_wp*pi) <= sum_tol*4.0_wp*pi, &
            & "integrate of a unit integrand is not 4*pi")
         if (allocated(error)) return
      end do
   end subroutine test_weights_four_pi

   !> Every monomial x**a*y**b*z**c up to the rule's degree is exact
   !>
   !> - Rules 6 to 266 points (degrees 3 to 27), including 74, 230, 266;
   !>   beyond them the monomial errors at degree d+1 approach round-off
   !> - At least one monomial of degree d+1 is not exact, so the degree is sharp
   !>
   !> @param[out] error  Test failure
   subroutine test_monomial_moments(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of rules checked, orders 1..13
      integer, parameter :: nchecked = 13
      type(moist_math_grid_angular_type) :: grid
      real(wp), allocatable :: px(:, :), py(:, :), pz(:, :), samples(:)
      real(wp) :: q, max_low, max_high
      integer :: order, d, a, b, c, n, i
      character(len=120) :: msg

      do order = 1, nchecked
         call build_exact(error, grid_size(order), grid)
         if (allocated(error)) return
         d = grid%degree
         n = grid%npts
         allocate (px(n, 0:d + 1), py(n, 0:d + 1), pz(n, 0:d + 1), samples(n))
         px(:, 0) = 1.0_wp
         py(:, 0) = 1.0_wp
         pz(:, 0) = 1.0_wp
         do a = 1, d + 1
            px(:, a) = px(:, a - 1)*grid%points(1, :)
            py(:, a) = py(:, a - 1)*grid%points(2, :)
            pz(:, a) = pz(:, a - 1)*grid%points(3, :)
         end do

         max_low = 0.0_wp
         max_high = 0.0_wp
         do a = 0, d + 1
            do b = 0, d + 1 - a
               do c = 0, d + 1 - a - b
                  do i = 1, n
                     samples(i) = px(i, a)*py(i, b)*pz(i, c)
                  end do
                  call grid%integrate_field(samples, q)
                  q = abs(q - monomial_integral(a, b, c))
                  if (a + b + c <= d) then
                     max_low = max(max_low, q)
                  else
                     max_high = max(max_high, q)
                  end if
               end do
            end do
         end do
         deallocate (px, py, pz, samples)

         write (msg, "(a,i0,a,i0,a,es10.3)") "Lebedev rule with ", n, &
            & " points misses a monomial of degree <= ", d, ": ", max_low
         call check(error, max_low <= mono_exact_tol, trim(msg))
         if (allocated(error)) return
         write (msg, "(a,i0,a,i0)") "Lebedev rule with ", n, &
            & " points is exact for every monomial of degree ", d + 1
         call check(error, max_high >= mono_inexact_margin, trim(msg))
         if (allocated(error)) return
      end do
   end subroutine test_monomial_moments

   !> The recorded degree is exact and sharp on all 32 rules
   !>
   !> - Zonal harmonics P_l(a.u) integrate to 0 for 1 <= l <= degree and not
   !>   for l = degree + 1, for two directions off the symmetry axes
   !> - Measured: at most 3.9e-15 below the degree, at least 9.6e-3 above
   !>
   !> @param[out] error  Test failure
   subroutine test_achieved_degree(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      real(wp), allocatable :: t(:), p0(:), p1(:), p2(:)
      real(wp) :: dirs(3, 2), q, max_low
      integer :: order, idir, l, n, d
      character(len=120) :: msg

      dirs(:, 1) = [0.3_wp, -0.5_wp, 0.8_wp]
      dirs(:, 2) = [0.71_wp, 0.23_wp, -0.41_wp]
      dirs(:, 1) = dirs(:, 1)/norm2(dirs(:, 1))
      dirs(:, 2) = dirs(:, 2)/norm2(dirs(:, 2))

      do order = 1, nrules
         call build_exact(error, grid_size(order), grid)
         if (allocated(error)) return
         call check(error, grid%degree, lebedev_degree_table(order), "Grid degree differs from the table")
         if (allocated(error)) return
         d = grid%degree
         n = grid%npts
         allocate (t(n), p0(n), p1(n), p2(n))
         do idir = 1, 2
            t(:) = matmul(dirs(:, idir), grid%points)
            p0(:) = 1.0_wp
            p1(:) = t
            call grid%integrate_field(p1, q)
            max_low = abs(q)
            do l = 1, d
               p2(:) = (real(2*l + 1, wp)*t*p1 - real(l, wp)*p0)/real(l + 1, wp)
               call grid%integrate_field(p2, q)
               if (l + 1 <= d) max_low = max(max_low, abs(q))
               p0(:) = p1
               p1(:) = p2
            end do
            ! q now holds the integral of P_(d+1)
            write (msg, "(a,i0,a,i0,a,es10.3)") "Lebedev rule with ", n, &
               & " points is not exact through degree ", d, ": ", max_low
            call check(error, max_low <= exact_tol, trim(msg))
            if (allocated(error)) return
            write (msg, "(a,i0,a,i0)") "Lebedev rule with ", n, &
               & " points is exact beyond its degree ", d
            call check(error, abs(q) >= inexact_margin, trim(msg))
            if (allocated(error)) return
         end do
         deallocate (t, p0, p1, p2)
      end do
   end subroutine test_achieved_degree

   !> `integrate` and `integrate_field` agree to round-off
   !>
   !> @param[out] error  Test failure
   subroutine test_integrate_matches_field(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      real(wp), allocatable :: samples(:)
      real(wp) :: from_field, from_analytic
      integer :: i

      call build_exact(error, 230, grid)
      if (allocated(error)) return

      allocate (samples(grid%npts))
      do i = 1, grid%npts
         samples(i) = f_smooth(grid%points(:, i))
      end do
      call grid%integrate_field(samples, from_field)
      call grid%integrate(f_smooth, from_analytic)

      call check(error, abs(from_field - from_analytic) <= unit_tol*abs(from_field), &
         & "integrate and integrate_field disagree")
   end subroutine test_integrate_matches_field

   !> `destroy` is idempotent and resets the grid
   !>
   !> @param[out] error  Test failure
   subroutine test_destroy(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid

      call build_exact(error, 110, grid)
      if (allocated(error)) return

      call grid%destroy()
      call grid%destroy()
      call check(error, grid%npts == 0 .and. grid%degree == 0 .and. .not. allocated(grid%points) &
         & .and. .not. allocated(grid%weights), "destroy did not reset the grid")
   end subroutine test_destroy

   !* ================================================================================= *!
   !*                        Lebedev grids from new_lebedev_grid                        *!
   !* ================================================================================= *!

   !> Weights must sum to 4*pi (solid angle of the unit sphere)
   !>
   !> Checked both directly and through `integrate_field` of a unit field,
   !> each divided by 4*pi
   !>
   !> @param[out] error  Test failure
   subroutine test_weights_sum(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: ones(:)
      real(wp) :: total
      integer :: k

      do k = 1, size(npts_list)
         call new_lebedev_grid(grid, merr, npts=npts_list(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if

         call check(error, grid%npts, npts_list(k), "Unexpected number of angular nodes")
         if (allocated(error)) return

         call check(error, abs(sum(grid%weights)/(4.0_wp*pi) - 1.0_wp) < struct_tol, &
            & "Angular weights do not sum to 4*pi")
         if (allocated(error)) return

         allocate (ones(grid%npts), source=1.0_wp)
         call grid%integrate_field(ones, total)
         deallocate (ones)
         call check(error, abs(total/(4.0_wp*pi) - 1.0_wp) < struct_tol, &
            & "integrate_field of a unit field is not 4*pi")
         if (allocated(error)) return
      end do
   end subroutine test_weights_sum

   !> Nodes must be Cartesian unit vectors
   !>
   !> @param[out] error  Test failure
   subroutine test_unit_vectors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: r2
      integer :: k, i

      do k = 1, size(npts_list)
         call new_lebedev_grid(grid, merr, npts=npts_list(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if

         call check(error, size(grid%points, 1), 3, "Angular points are not 3-vectors")
         if (allocated(error)) return

         do i = 1, grid%npts
            r2 = sum(grid%points(:, i)**2)
            call check(error, abs(r2 - 1.0_wp) < struct_tol, &
               & "Angular node is not a unit vector")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_unit_vectors

   !> Exactly one of npts= / degree= must be given
   !>
   !> @param[out] error  Test failure
   subroutine test_selector_validation(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_lebedev_grid(grid, merr)
      call check(error, allocated(merr), "Missing npts=/degree= was accepted")
      if (allocated(error)) return
      deallocate (merr)

      call new_lebedev_grid(grid, merr, npts=74, degree=13)
      call check(error, allocated(merr), "Both npts= and degree= were accepted")
      if (allocated(error)) return
      call check(error, grid%npts, 0, "Rejected grid was left initialized")
   end subroutine test_selector_validation

   !> An unsupported point count must be rejected
   !>
   !> @param[out] error  Test failure
   subroutine test_unsupported_npts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_lebedev_grid(grid, merr, npts=7)
      call check(error, allocated(merr), "Unsupported Lebedev size was accepted")
      if (allocated(error)) return
      call check(error, grid%npts, 0, "Rejected grid was left initialized")
   end subroutine test_unsupported_npts

   !> degree= selects the smallest rule reaching the requested exactness
   !>
   !> %degree reports a value at least as large as the request
   !>
   !> @param[out] error  Test failure
   subroutine test_degree_selection(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Requested exactness degrees
      integer, parameter :: request(5) = [0, 3, 10, 12, 131]
      !> Point count of the smallest rule reaching each request
      integer, parameter :: expect_npts(5) = [6, 6, 50, 74, 5810]
      !> Exactness degree of that rule
      integer, parameter :: expect_degree(5) = [3, 3, 11, 13, 131]
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, size(request)
         call new_lebedev_grid(grid, merr, degree=request(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if

         call check(error, grid%npts, expect_npts(k), &
            & "degree= did not select the smallest sufficient rule")
         if (allocated(error)) return

         call check(error, grid%degree, expect_degree(k), &
            & "Unexpected reported exactness degree")
         if (allocated(error)) return

         call check(error, grid%degree >= request(k), &
            & "Reported degree is below the requested degree")
         if (allocated(error)) return
      end do
   end subroutine test_degree_selection

   !> Degrees outside the supported range must be rejected
   !>
   !> @param[out] error  Test failure
   subroutine test_degree_out_of_range(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_lebedev_grid(grid, merr, degree=132)
      call check(error, allocated(merr), "Degree beyond the highest rule was accepted")
      if (allocated(error)) return
      deallocate (merr)

      call new_lebedev_grid(grid, merr, degree=-1)
      call check(error, allocated(merr), "Negative degree was accepted")
   end subroutine test_degree_out_of_range

   !> A rule of algebraic degree d reproduces exact sphere averages of polynomials
   !>
   !> In the Cartesian components, up to that degree; the solid-angle
   !> integral divided by 4*pi is the average, so
   !>   <1>=1, <x>=<x*y>=0, <3z^2-1>=0, <x^2>=1/3,
   !>   <x^4+y^4+z^4>=3/5
   !>
   !> @param[out] error  Test failure
   subroutine test_polynomial_exactness(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, size(npts_list)
         call new_lebedev_grid(grid, merr, npts=npts_list(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if

         call check_average(error, grid, f_one, 1.0_wp, "<1>")
         if (allocated(error)) return
         call check_average(error, grid, f_x, 0.0_wp, "<x>")
         if (allocated(error)) return
         call check_average(error, grid, f_xy, 0.0_wp, "<x*y>")
         if (allocated(error)) return
         call check_average(error, grid, f_3z2m1, 0.0_wp, "<3z^2-1>")
         if (allocated(error)) return
         call check_average(error, grid, f_x2, 1.0_wp/3.0_wp, "<x^2>")
         if (allocated(error)) return

         ! Degree 4; the 6-point rule is only exact to degree 3
         if (grid%degree >= 4) then
            call check_average(error, grid, f_quartic, 0.6_wp, "<x^4+y^4+z^4>")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_polynomial_exactness

   !> Every supported rule integrates a polynomial of its full degree exactly
   !>
   !> - Rules visited through `degree=`: request 0, then one above the
   !>   reported degree of the rule just built, until the request exceeds the
   !>   table; all 32 rules must be visited
   !> - Integrand `(a.u)^n` with `n` = degree - 1, the highest nontrivial
   !>   power: the degrees are odd, and odd powers average to zero on any
   !>   inversion-symmetric rule
   !> - Fixed direction `a` off every symmetry axis of the octahedral rules;
   !>   exact sphere average `1/(n+1)`, the solid-angle integral over 4*pi
   !> - Relative error; measured at most 3.6e-15 over all rules
   !>
   !> @param[out] error  Test failure
   subroutine test_full_degree_exactness(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of supported Lebedev rules
      integer, parameter :: nrules = 32
      !> Acceptance bound on the relative error
      real(wp), parameter :: thr = 2.0e-14_wp
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: samples(:)
      real(wp) :: dir(3), average, exact
      integer :: request, n, i, visited
      character(len=80) :: msg

      dir = [0.3_wp, -0.5_wp, 0.8_wp]
      dir = dir/norm2(dir)
      request = 0
      visited = 0
      do
         call new_lebedev_grid(grid, merr, degree=request)
         if (allocated(merr)) exit
         visited = visited + 1

         n = grid%degree - 1
         allocate (samples(grid%npts))
         do i = 1, grid%npts
            samples(i) = dot_product(dir, grid%points(:, i))**n
         end do
         call grid%integrate_field(samples, average)
         average = average/(4.0_wp*pi)
         deallocate (samples)

         exact = 1.0_wp/real(n + 1, wp)
         write (msg, "(a,i0,a,i0)") "Lebedev rule with ", grid%npts, &
            & " points is not exact at degree ", n
         call check(error, abs(average - exact) <= thr*exact, trim(msg))
         if (allocated(error)) return
         request = grid%degree + 1
      end do

      call check(error, visited, nrules, "degree= selection did not visit every Lebedev rule")
   end subroutine test_full_degree_exactness

   !> `integrate` and `integrate_field` must agree to round-off
   !>
   !> Analytic integrand vs. tabulated samples, on the same function;
   !> both compared as sphere averages
   !>
   !> @param[out] error  Test failure
   subroutine test_integrate_agreement(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: samples(:)
      real(wp) :: from_field, from_analytic
      integer :: i

      call new_lebedev_grid(grid, merr, npts=230)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      allocate (samples(grid%npts))
      do i = 1, grid%npts
         samples(i) = f_mixed(grid%points(:, i))
      end do

      call grid%integrate_field(samples, from_field)
      call grid%integrate(f_mixed, from_analytic)
      from_field = from_field/(4.0_wp*pi)
      from_analytic = from_analytic/(4.0_wp*pi)

      call check(error, abs(from_field - from_analytic) < struct_tol, &
         & "integrate and integrate_field disagree")
   end subroutine test_integrate_agreement

   !> `destroy` must be idempotent and reset the grid to its empty state
   !>
   !> @param[out] error  Test failure
   subroutine test_destroy_idempotent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_lebedev_grid(grid, merr, npts=110)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call grid%destroy()
      call grid%destroy()

      call check(error, grid%npts, 0, "destroy did not reset npts")
      if (allocated(error)) return
      call check(error, grid%degree, 0, "destroy did not reset the degree")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%points), "destroy left points allocated")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%weights), "destroy left weights allocated")
   end subroutine test_destroy_idempotent

   !* ================================================================================= *!
   !*                                 Lebedev metadata                                  *!
   !* ================================================================================= *!

   !> The negative-weight list matches the signs of the generated tables
   !>
   !> Measured minimum weights: 74 -> -2.96e-2, 230 -> -5.52e-2, 266 -> -2.52e-3
   !>
   !> @param[out] error  Test failure
   subroutine test_negative_weight_list(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: x(:, :), w(:)
      integer :: order, nneg
      character(len=80) :: msg

      nneg = 0
      do order = 1, nrules
         allocate (x(3, grid_size(order)), w(grid_size(order)))
         call get_angular_grid(order, x, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         if (any(w < 0.0_wp)) nneg = nneg + 1
         write (msg, "(a,i0,a)") "Negative-weight flag of the ", grid_size(order), &
            & "-point rule does not match its table"
         call check(error, any(w < 0.0_wp) .eqv. lebedev_has_negative_weights(order), trim(msg))
         if (allocated(error)) return
         call check(error, any(w < 0.0_wp) .eqv. &
            & any(grid_size(order) == lebedev_negative_weight_sizes), trim(msg))
         if (allocated(error)) return
         deallocate (x, w)
      end do
      call check(error, nneg, size(lebedev_negative_weight_sizes), &
         & "Negative-weight list has the wrong length")
   end subroutine test_negative_weight_list

   !> `lebedev_order_from_num`: filter, error messages, order 0 on error
   !>
   !> @param[out] error  Test failure
   subroutine test_order_from_num_filter(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      integer :: i, order

      do i = 1, nrules
         call lebedev_order_from_num(grid_size(i), order, merr)
         call check(error, .not. allocated(merr) .and. order == i, "Supported size was rejected")
         if (allocated(error)) return
         call lebedev_order_from_num(grid_size(i), order, merr, positive_weights_only=.false.)
         call check(error, .not. allocated(merr) .and. order == i, &
            & "Supported size was rejected without the filter")
         if (allocated(error)) return

         call lebedev_order_from_num(grid_size(i), order, merr, positive_weights_only=.true.)
         if (lebedev_has_negative_weights(i)) then
            call check(error, allocated(merr) .and. order == 0, &
               & "Size with negative weights passed the filter")
         else
            call check(error, .not. allocated(merr) .and. order == i, &
               & "Size with positive weights failed the filter")
         end if
         if (allocated(error)) return
      end do

      call lebedev_order_from_num(74, order, merr, positive_weights_only=.true.)
      call check(error, allocated(merr), "74 passed the filter")
      if (allocated(error)) return
      call check(error, merr%message == msg_neg74, "Unexpected filter message: "//merr%message)
      if (allocated(error)) return

      call lebedev_order_from_num(7, order, merr)
      call check(error, allocated(merr) .and. order == 0, "Unsupported size 7 was accepted")
      if (allocated(error)) return
      call check(error, index(merr%message, "Unsupported Lebedev size 7;") == 1 &
         & .and. index(merr%message, "6, 14, 26") > 0 .and. index(merr%message, "5294, 5810") > 0, &
         & "Unsupported-size message lacks the size or the supported list: "//merr%message)
   end subroutine test_order_from_num_filter

   !> Lebedev weights should sum to 1 on the unit sphere
   !>
   !> @param[out] error  Test failure
   subroutine test_angular_weights_sum(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: table_npts(3) = [74, 230, 434]
      integer, parameter :: orders(3)    = [ 6,  12,  16]
      real(wp), allocatable :: xyz(:, :), w(:)
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, 3
         allocate (xyz(3, table_npts(k)), w(table_npts(k)))
         call get_angular_grid(orders(k), xyz, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, abs(sum(w) - 1.0_wp) < 1.0e-12_wp, &
            & "Lebedev weights do not sum to 1")
         deallocate (xyz, w)
         if (allocated(error)) return
      end do
   end subroutine test_angular_weights_sum

   !> Lebedev points should lie exactly on the unit sphere
   !>
   !> @param[out] error  Test failure
   subroutine test_angular_unit_sphere(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: table_npts(3) = [74, 230, 434]
      integer, parameter :: orders(3)    = [ 6,  12,  16]
      real(wp), allocatable :: xyz(:, :), w(:)
      type(mctc_error), allocatable :: merr
      integer :: k, i
      real(wp) :: r2

      do k = 1, 3
         allocate (xyz(3, table_npts(k)), w(table_npts(k)))
         call get_angular_grid(orders(k), xyz, w, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         do i = 1, table_npts(k)
            r2 = xyz(1, i)**2 + xyz(2, i)**2 + xyz(3, i)**2
            call check(error, abs(r2 - 1.0_wp) < 1.0e-12_wp, &
               & "Lebedev point off unit sphere")
            if (allocated(error)) exit
         end do
         deallocate (xyz, w)
         if (allocated(error)) return
      end do
   end subroutine test_angular_unit_sphere

   !* ================================================================================= *!
   !*                                 Lebedev selection                                 *!
   !* ================================================================================= *!

   !> The default generator admits every rule, including negative weights
   !>
   !> @param[out] error  Test failure
   subroutine test_default_all_weights(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: gen
      type(moist_math_grid_angular_request_type) :: req
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: k

      call check(error, .not. gen%positive_weights_only, "Filter is on by default")
      if (allocated(error)) return

      do k = 1, size(lebedev_negative_weight_sizes)
         req = moist_math_grid_angular_request_type(npts=lebedev_negative_weight_sizes(k))
         call gen%select(req, grid, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, any(grid%weights < 0.0_wp), "Default exact request lost the negative weights")
         if (allocated(error)) return
      end do

      req = moist_math_grid_angular_request_type(min_degree=13)
      call expect_npts(error, gen, req, 74, "default, min_degree 13")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_degree=24)
      call expect_npts(error, gen, req, 230, "default, min_degree 24")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=70.0_wp)
      call expect_npts(error, gen, req, 74, "default, target 70")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=250.0_wp)
      call expect_npts(error, gen, req, 266, "default, target 250")
   end subroutine test_default_all_weights

   !> The opt-in filter rejects exact negative-weight requests and skips them otherwise
   !>
   !> @param[out] error  Test failure
   subroutine test_positive_only(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: gen
      type(moist_math_grid_angular_request_type) :: req
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: k, degree

      call new_lebedev_generator(gen, positive_weights_only=.true.)
      call check(error, gen%positive_weights_only, "Constructor ignored positive_weights_only")
      if (allocated(error)) return

      req = moist_math_grid_angular_request_type(npts=74)
      call expect_failure(error, gen, req, "filter, exact 74", msg_neg74)
      if (allocated(error)) return
      do k = 1, size(lebedev_negative_weight_sizes)
         req = moist_math_grid_angular_request_type(npts=lebedev_negative_weight_sizes(k))
         call expect_failure(error, gen, req, "filter, exact negative-weight size")
         if (allocated(error)) return
      end do

      req = moist_math_grid_angular_request_type(min_degree=13)
      call expect_npts(error, gen, req, 86, "filter, min_degree 13")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_degree=24)
      call expect_npts(error, gen, req, 302, "filter, min_degree 24")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=70.0_wp)
      call expect_npts(error, gen, req, 86, "filter, target 70")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=200.0_wp)
      call expect_npts(error, gen, req, 302, "filter, target 200")
      if (allocated(error)) return

      ! Every degree request under the filter yields positive weights only
      do degree = 0, lebedev_degree_table(nrules)
         req = moist_math_grid_angular_request_type(min_degree=degree)
         call gen%select(req, grid, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(grid%weights > 0.0_wp), "Filtered selection has a non-positive weight")
         if (allocated(error)) return
      end do
   end subroutine test_positive_only

   !> Exact requests are strict: no substitution, hard constraints still apply
   !>
   !> @param[out] error  Test failure
   subroutine test_exact_request_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: gen
      type(moist_math_grid_angular_request_type) :: req
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      req = moist_math_grid_angular_request_type(npts=7)
      call expect_failure(error, gen, req, "exact 7")
      if (allocated(error)) return
      call gen%select(req, grid, merr)
      call check(error, index(merr%message, "Unsupported Lebedev size 7;") == 1, &
         & "Exact request error does not name the size: "//merr%message)
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(npts=5811)
      call expect_failure(error, gen, req, "exact 5811")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(npts=-6)
      call expect_failure(error, gen, req, "exact -6")
      if (allocated(error)) return

      req = moist_math_grid_angular_request_type(npts=26, min_degree=9)
      call expect_failure(error, gen, req, "exact 26 below min_degree 9")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(npts=110, max_points=100)
      call expect_failure(error, gen, req, "exact 110 above max_points 100")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(npts=110, min_points=146)
      call expect_failure(error, gen, req, "exact 110 below min_points 146")
      if (allocated(error)) return

      ! Satisfied hard constraints pass; the soft target does not move an exact request
      req = moist_math_grid_angular_request_type(npts=110, min_degree=17, min_points=110, &
         & max_points=110, target_points=5000.0_wp)
      call expect_npts(error, gen, req, 110, "exact 110 with matching constraints")
   end subroutine test_exact_request_errors

   !> Minimum degree selects the smallest admissible rule reaching it
   !>
   !> @param[out] error  Test failure
   subroutine test_min_degree_selection(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Requested minimum degrees
      integer, parameter :: request(14) = [0, 3, 4, 10, 12, 13, 14, 24, 26, 28, 32, 36, 130, 131]
      !> Selected sizes, default generator
      integer, parameter :: expect_all(14) = [6, 6, 14, 50, 74, 74, 86, 230, 266, 302, 434, 590, &
         & 5810, 5810]
      !> Selected sizes, positive weights only
      integer, parameter :: expect_pos(14) = [6, 6, 14, 50, 86, 86, 86, 302, 302, 302, 434, 590, &
         & 5810, 5810]
      type(moist_math_grid_angular_generator_lebedev_type) :: gen_all, gen_pos
      type(moist_math_grid_angular_request_type) :: req
      integer :: k
      character(len=40) :: label

      call new_lebedev_generator(gen_pos, positive_weights_only=.true.)
      do k = 1, size(request)
         req = moist_math_grid_angular_request_type(min_degree=request(k))
         write (label, "(a,i0)") "min_degree ", request(k)
         call expect_npts(error, gen_all, req, expect_all(k), "default, "//trim(label))
         if (allocated(error)) return
         call expect_npts(error, gen_pos, req, expect_pos(k), "filter, "//trim(label))
         if (allocated(error)) return
      end do

      ! A hard degree wins over a smaller target and composes with a larger one
      req = moist_math_grid_angular_request_type(min_degree=17, target_points=50.0_wp)
      call expect_npts(error, gen_all, req, 110, "min_degree 17, target 50")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_degree=17, target_points=200.0_wp)
      call expect_npts(error, gen_all, req, 230, "default, min_degree 17, target 200")
      if (allocated(error)) return
      call expect_npts(error, gen_pos, req, 302, "filter, min_degree 17, target 200")
   end subroutine test_min_degree_selection

   !> Soft target with hard floors and caps
   !>
   !> @param[out] error  Test failure
   subroutine test_target_selection(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: gen_all, gen_pos
      type(moist_math_grid_angular_request_type) :: req

      call new_lebedev_generator(gen_pos, positive_weights_only=.true.)

      ! No target: smallest admissible rule
      req = moist_math_grid_angular_request_type()
      call expect_npts(error, gen_all, req, 6, "empty request")
      if (allocated(error)) return

      ! Comparison is size >= target
      req = moist_math_grid_angular_request_type(target_points=590.0_wp)
      call expect_npts(error, gen_all, req, 590, "target 590")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=nearest(590.0_wp, -1.0_wp))
      call expect_npts(error, gen_all, req, 590, "target just below 590")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=nearest(590.0_wp, 1.0_wp))
      call expect_npts(error, gen_all, req, 770, "target just above 590")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=800.0_wp)
      call expect_npts(error, gen_all, req, 974, "target 800")
      if (allocated(error)) return

      ! Unreachable target inside the cap: largest admissible rule
      req = moist_math_grid_angular_request_type(target_points=800.0_wp, max_points=600)
      call expect_npts(error, gen_all, req, 590, "default, target 800, cap 600")
      if (allocated(error)) return
      call expect_npts(error, gen_pos, req, 590, "filter, target 800, cap 600")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=1.0e9_wp)
      call expect_npts(error, gen_all, req, 5810, "target 1e9")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=ieee_value(1.0_wp, ieee_positive_inf))
      call expect_npts(error, gen_all, req, 5810, "target +inf")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=1000.0_wp, max_points=300)
      call expect_npts(error, gen_all, req, 266, "default, target 1000, cap 300")
      if (allocated(error)) return
      call expect_npts(error, gen_pos, req, 194, "filter, target 1000, cap 300")
      if (allocated(error)) return

      ! Floor above the target
      req = moist_math_grid_angular_request_type(target_points=10.0_wp, min_points=110)
      call expect_npts(error, gen_all, req, 110, "target 10, floor 110")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=10.0_wp, min_points=100)
      call expect_npts(error, gen_all, req, 110, "target 10, floor 100")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=10.0_wp, min_points=74)
      call expect_npts(error, gen_all, req, 74, "default, target 10, floor 74")
      if (allocated(error)) return
      call expect_npts(error, gen_pos, req, 86, "filter, target 10, floor 74")
      if (allocated(error)) return

      ! Floor and cap of today's banded rule
      req = moist_math_grid_angular_request_type(target_points=250.0_wp, min_points=110, max_points=302)
      call expect_npts(error, gen_all, req, 266, "default, target 250 in [110, 302]")
      if (allocated(error)) return
      call expect_npts(error, gen_pos, req, 302, "filter, target 250 in [110, 302]")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=1000.0_wp, min_points=110, max_points=290)
      call expect_npts(error, gen_pos, req, 194, "filter, target 1000 in [110, 290]")
   end subroutine test_target_selection

   !> Incompatible or out-of-range hard constraints are errors, never relaxed
   !>
   !> @param[out] error  Test failure
   subroutine test_incompatible_constraints(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: gen_all, gen_pos
      type(moist_math_grid_angular_request_type) :: req

      call new_lebedev_generator(gen_pos, positive_weights_only=.true.)

      ! Degree 41 needs 590 points
      req = moist_math_grid_angular_request_type(min_degree=41, max_points=434)
      call expect_failure(error, gen_all, req, "min_degree 41, cap 434")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_degree=41, max_points=590)
      call expect_npts(error, gen_all, req, 590, "min_degree 41, cap 590")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_degree=132)
      call expect_failure(error, gen_all, req, "min_degree 132")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_points=700, max_points=600)
      call expect_failure(error, gen_all, req, "floor 700 above cap 600")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_points=5811)
      call expect_failure(error, gen_all, req, "floor 5811")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(max_points=5)
      call expect_failure(error, gen_all, req, "cap 5")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_points=200, max_points=290)
      call expect_npts(error, gen_all, req, 230, "default, [200, 290]")
      if (allocated(error)) return
      call expect_failure(error, gen_pos, req, "filter, [200, 290]")
      if (allocated(error)) return

      ! Out-of-range fields
      req = moist_math_grid_angular_request_type(min_degree=-1)
      call expect_failure(error, gen_all, req, "min_degree -1")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(min_points=-1)
      call expect_failure(error, gen_all, req, "min_points -1")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(max_points=-1)
      call expect_failure(error, gen_all, req, "max_points -1")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=-1.0_wp)
      call expect_failure(error, gen_all, req, "target -1")
      if (allocated(error)) return
      req = moist_math_grid_angular_request_type(target_points=ieee_value(1.0_wp, ieee_quiet_nan))
      call expect_failure(error, gen_all, req, "target NaN")
   end subroutine test_incompatible_constraints

   !> `new_lebedev_grid`: one selector, exact count or degree, filter pass-through
   !>
   !> @param[out] error  Test failure
   subroutine test_new_lebedev_grid(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_lebedev_grid(grid, merr)
      call check(error, allocated(merr) .and. grid%npts == 0, "No selector was accepted")
      if (allocated(error)) return
      call new_lebedev_grid(grid, merr, npts=74, degree=13)
      call check(error, allocated(merr) .and. grid%npts == 0, "Both selectors were accepted")
      if (allocated(error)) return
      call new_lebedev_grid(grid, merr, npts=0)
      call check(error, allocated(merr) .and. grid%npts == 0, "npts=0 was accepted")
      if (allocated(error)) return
      call new_lebedev_grid(grid, merr, degree=-1)
      call check(error, allocated(merr) .and. grid%npts == 0, "Negative degree was accepted")
      if (allocated(error)) return
      call new_lebedev_grid(grid, merr, degree=132)
      call check(error, allocated(merr) .and. grid%npts == 0, "Degree 132 was accepted")
      if (allocated(error)) return

      call new_lebedev_grid(grid, merr, npts=74, positive_weights_only=.true.)
      call check(error, allocated(merr), "Filtered exact 74 was accepted")
      if (allocated(error)) return
      call check(error, merr%message == msg_neg74, "Unexpected filter message: "//merr%message)
      if (allocated(error)) return

      call new_lebedev_grid(grid, merr, npts=74)
      call check(error, .not. allocated(merr) .and. grid%npts == 74 .and. grid%degree == 13, &
         & "Exact 74 failed")
      if (allocated(error)) return
      call new_lebedev_grid(grid, merr, degree=13)
      call check(error, .not. allocated(merr) .and. grid%npts == 74, "degree 13 did not select 74")
      if (allocated(error)) return
      call new_lebedev_grid(grid, merr, degree=13, positive_weights_only=.true.)
      call check(error, .not. allocated(merr) .and. grid%npts == 86 .and. grid%degree == 15, &
         & "Filtered degree 13 did not select 86")
   end subroutine test_new_lebedev_grid

   !> Selection dispatches through the abstract generator
   !>
   !> @param[out] error  Test failure
   subroutine test_polymorphic_generator(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: leb
      class(moist_math_grid_angular_generator_type), allocatable :: gen, copy
      type(moist_math_grid_angular_request_type) :: req

      call new_lebedev_generator(leb, positive_weights_only=.true.)
      allocate (gen, source=leb)
      allocate (copy, source=gen)
      deallocate (gen)

      req = moist_math_grid_angular_request_type(min_degree=13)
      call expect_npts(error, copy, req, 86, "polymorphic copy, min_degree 13")
      if (allocated(error)) return
      select type (copy)
      type is (moist_math_grid_angular_generator_lebedev_type)
         call check(error, copy%positive_weights_only, "Copy lost positive_weights_only")
      class default
         call test_failed(error, "Copy has the wrong dynamic type")
      end select
   end subroutine test_polymorphic_generator

end module test_math_grid_angular
