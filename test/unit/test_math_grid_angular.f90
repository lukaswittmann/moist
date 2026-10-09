!> Test suite for the angular grid and the Lebedev selection
!>
!> - Angular grid: unit directions, weights summing to 4*pi, monomial moments,
!>   achieved degree checked with zonal harmonics on all 32 rules
!> - Metadata: the negative-weight list matches the signs of the tables
!> - Selection: positive-only opt-in, minimum degree, target with floors and caps
module test_math_grid_angular
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_grid_angular_lebedev, only: grid_size, lebedev_degree_table, &
      & lebedev_negative_weight_sizes, lebedev_has_negative_weights, get_angular_grid
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, &
      & moist_math_grid_angular_generator_type, moist_math_grid_angular_request_type, &
      & moist_math_grid_angular_generator_lebedev_type, new_lebedev_generator, new_lebedev_grid
   implicit none(type, external)
   private

   public :: collect_math_grid_angular

   !> Number of supported Lebedev rules
   integer, parameter :: nrules = 32
   !> Point counts of a representative spread of rules, with and without negative weights
   integer, parameter :: npts_list_wide(7) = [6, 26, 74, 110, 434, 590, 5810]
   !> Tolerance of the unit-vector checks
   real(wp), parameter :: unit_tol = 1.0e-14_wp
   !> Relative tolerance of the total weight; measured at most 2.0e-14 (3890 points)
   real(wp), parameter :: sum_tol = 5.0e-14_wp
   !> Bound on a zonal moment the rule integrates exactly; measured at most 3.9e-15
   real(wp), parameter :: exact_tol = 1.0e-13_wp
   !> Lower bound on the zonal error one degree above the rule; measured at least 9.6e-3
   real(wp), parameter :: inexact_margin = 1.0e-3_wp
   !> Bound on a monomial moment the rule integrates exactly; measured at most 5.9e-14
   real(wp), parameter :: mono_exact_tol = 1.0e-12_wp
   !> Worst monomial error lower bound one degree above the rule; measured at least 4.2e-8
   real(wp), parameter :: mono_inexact_margin = 1.0e-9_wp

contains

   !> Collect all math_grid_angular tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_angular(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("angular_grid_unit_directions", test_unit_directions), &
         new_unittest("angular_grid_weights_sum_to_four_pi", test_weights_four_pi), &
         new_unittest("angular_grid_monomial_moments", test_monomial_moments), &
         new_unittest("angular_grid_achieved_degree", test_achieved_degree), &
         new_unittest("lebedev_negative_weight_list", test_negative_weight_list), &
         new_unittest("lebedev_positive_only_opt_in", test_positive_only), &
         new_unittest("lebedev_min_degree_selection", test_min_degree_selection), &
         new_unittest("lebedev_target_selection", test_target_selection) &
         ]
   end subroutine collect_math_grid_angular

   !* ------------------------------------ Helpers ------------------------------------ *!

   !> Run a selection that must succeed and check the selected size
   !>
   !> Also checks that the grid is consistently filled for that size
   !>
   !> @param[out] error      test failure
   !> @param[in]  generator  generator under test
   !> @param[in]  request    angular request
   !> @param[in]  expected   expected point count
   !> @param[in]  label      case description for the failure message
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

   !> Build the grid of an exact point count with the default generator
   !>
   !> @param[in]  npts  supported point count
   !> @param[out] grid  generated grid
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
   !> @param[in] a  exponent of x
   !> @param[in] b  exponent of y
   !> @param[in] c  exponent of z
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

   !> Quadratic integrand f(u) = 1 + x/2 + 3 z**2; its sphere integral is 8*pi
   !>
   !> f(1, 0, 0) = 1.5 differs from the spherical average 2, so a callback
   !> evaluated only at the first Lebedev node cannot pass
   !>
   !> @param[in] uvec  Cartesian unit vector
   pure function f_quadratic(uvec) result(val)
      !> Cartesian unit vector
      real(wp), intent(in) :: uvec(3)
      !> Function value
      real(wp) :: val

      val = 1.0_wp + 0.5_wp*uvec(1) + 3.0_wp*uvec(3)**2
   end function f_quadratic

   !* ---------------------------------- Angular grid --------------------------------- *!

   !> Nodes are Cartesian unit vectors and arrays match npts
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

   !> Weights sum to 4*pi; both integrators agree on a quadratic integrand
   !>
   !> - A field of the wrong length is rejected, and the next valid call succeeds
   !> - Every rule integrates 1 + x/2 + 3 z**2 exactly (8*pi), as a callback and
   !>   as tabulated samples
   subroutine test_weights_four_pi(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: samples(:)
      real(wp) :: total
      integer :: k, i

      do k = 1, size(npts_list_wide)
         call build_exact(error, npts_list_wide(k), grid)
         if (allocated(error)) return

         call check(error, abs(sum(grid%weights) - 4.0_wp*pi) <= sum_tol*4.0_wp*pi, &
            & "Angular weights do not sum to 4*pi")
         if (allocated(error)) return

         allocate (samples(grid%npts))
         do i = 1, grid%npts
            samples(i) = f_quadratic(grid%points(:, i))
         end do
         ! A short field is rejected; the valid call below still succeeds
         call grid%integrate_field(samples(2:), total, merr)
         call check(error, allocated(merr), "integrate_field accepted a short field")
         if (allocated(error)) return
         call grid%integrate_field(samples, total, merr)
         deallocate (samples)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, abs(total - 8.0_wp*pi) <= sum_tol*8.0_wp*pi, &
            & "integrate_field of the quadratic samples is not 8*pi")
         if (allocated(error)) return

         call grid%integrate(f_quadratic, total)
         call check(error, abs(total - 8.0_wp*pi) <= sum_tol*8.0_wp*pi, &
            & "integrate of the quadratic integrand is not 8*pi")
         if (allocated(error)) return
      end do
   end subroutine test_weights_four_pi

   !> Every monomial x**a*y**b*z**c up to the rule's degree is exact
   !>
   !> - Rules 6 to 266 points (degrees 3 to 27), including 74, 230, 266;
   !>   beyond them the monomial errors at degree d+1 approach round-off
   !> - At least one monomial of degree d+1 is not exact, so the degree is sharp
   subroutine test_monomial_moments(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of rules checked, orders 1..13
      integer, parameter :: nchecked = 13
      type(moist_math_grid_angular_type) :: grid
      real(wp), allocatable :: px(:, :), py(:, :), pz(:, :), samples(:)
      real(wp) :: q, max_low, max_high
      type(mctc_error), allocatable :: merr
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
                  call grid%integrate_field(samples, q, merr)
                  if (allocated(merr)) then
                     call test_failed(error, merr%message)
                     return
                  end if
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

   !> Check exactness and sharpness of the recorded degree on all 32 rules
   !>
   !> - Zonal harmonics P_l(a.u) integrate to 0 for 1 <= l <= degree and not
   !>   for l = degree + 1, for two directions off the symmetry axes
   !> - Measured: at most 3.9e-15 below the degree, at least 9.6e-3 above
   subroutine test_achieved_degree(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_type) :: grid
      real(wp), allocatable :: t(:), p0(:), p1(:), p2(:)
      real(wp) :: dirs(3, 2), q, max_low
      type(mctc_error), allocatable :: merr
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
            call grid%integrate_field(p1, q, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            max_low = abs(q)
            do l = 1, d
               p2(:) = (real(2*l + 1, wp)*t*p1 - real(l, wp)*p0)/real(l + 1, wp)
               call grid%integrate_field(p2, q, merr)
               if (allocated(merr)) then
                  call test_failed(error, merr%message)
                  return
               end if
               if (l + 1 <= d) max_low = max(max_low, abs(q))
               p0(:) = p1
               p1(:) = p2
            end do
            ! Integral of P_(d+1) in q
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

   !* -------------------------------- Lebedev metadata ------------------------------- *!

   !> Check the negative-weight list against the generated table signs
   !>
   !> Measured minimum weights: 74 -> -2.96e-2, 230 -> -5.52e-2, 266 -> -2.52e-3
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

   !* ------------------------------- Lebedev selection ------------------------------- *!

   !> Check exclusion of negative-weight rules by the opt-in filter
   subroutine test_positive_only(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_angular_generator_lebedev_type) :: gen
      type(moist_math_grid_angular_request_type) :: req
      type(moist_math_grid_angular_type) :: grid
      type(mctc_error), allocatable :: merr
      integer :: degree

      call new_lebedev_generator(gen, positive_weights_only=.true.)

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

   !> Minimum degree selects the smallest admissible rule reaching it
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

end module test_math_grid_angular
