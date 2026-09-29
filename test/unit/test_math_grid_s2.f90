!> Test suite for the moist math/grid/s2 submodule
!>
!> - Lebedev node and weight structure, selector validation, degree selection
!> - Polynomial exactness: low-order Cartesian moments on a spread of rules,
!>   and every one of the 32 rules at its full algebraic degree
!> - `integrate` vs `integrate_field`, no trafo, idempotent destroy
module test_math_grid_s2
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_grid_s2_base, only: moist_math_grid_s2_type, &
      & moist_math_grid_s2_trafo_type, integrand_s2
   use moist_math_grid_s2_lebedev, only: moist_math_grid_s2_lebedev_type, new_s2_grid_lebedev
   implicit none(type, external)
   private

   public :: collect_math_grid_s2

   !> Point counts of a representative spread of supported Lebedev rules
   integer, parameter :: npts_list(4) = [6, 74, 434, 590]
   !> Tolerance for the structural (sum-to-one, unit-vector) checks
   real(wp), parameter :: struct_tol = 1.0e-14_wp
   !> Tolerance for the quadrature exactness checks
   real(wp), parameter :: quad_tol = 1.0e-12_wp

contains

   !> Collect all math_grid_s2 tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_s2(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("s2_lebedev_weights_sum_to_one", test_weights_sum), &
         new_unittest("s2_lebedev_points_are_unit_vectors", test_unit_vectors), &
         new_unittest("s2_lebedev_requires_one_selector", test_selector_validation), &
         new_unittest("s2_lebedev_rejects_unsupported_npts", test_unsupported_npts), &
         new_unittest("s2_lebedev_degree_selection", test_degree_selection), &
         new_unittest("s2_lebedev_degree_out_of_range", test_degree_out_of_range), &
         new_unittest("s2_lebedev_polynomial_exactness", test_polynomial_exactness), &
         new_unittest("s2_lebedev_exact_at_full_degree", test_full_degree_exactness), &
         new_unittest("s2_integrate_matches_integrate_field", test_integrate_agreement), &
         new_unittest("s2_lebedev_has_no_trafo", test_no_trafo), &
         new_unittest("s2_lebedev_destroy_idempotent", test_destroy_idempotent) &
         ]
   end subroutine collect_math_grid_s2

   !> Weights must sum to 1 (Laikov convention)
   !>
   !> Checked both directly and through the base `integrate_field` of a unit field
   !>
   !> @param[out] error  Test failure
   subroutine test_weights_sum(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_s2_lebedev_type), target :: grid
      class(moist_math_grid_s2_type), pointer :: base
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: ones(:)
      real(wp) :: total
      integer :: k

      do k = 1, size(npts_list)
         call new_s2_grid_lebedev(grid, merr, npts=npts_list(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         base => grid

         call check(error, grid%npts, npts_list(k), "Unexpected number of S2 nodes")
         if (allocated(error)) return

         call check(error, abs(sum(grid%weights) - 1.0_wp) < struct_tol, &
            & "S2 weights do not sum to 1")
         if (allocated(error)) return

         allocate (ones(grid%npts), source=1.0_wp)
         call base%integrate_field(ones, total)
         deallocate (ones)
         call check(error, abs(total - 1.0_wp) < struct_tol, &
            & "integrate_field of a unit field is not 1")
         if (allocated(error)) return
      end do
   end subroutine test_weights_sum

   !> Nodes must be Cartesian unit vectors
   !>
   !> @param[out] error  Test failure
   subroutine test_unit_vectors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_s2_lebedev_type) :: grid
      type(mctc_error), allocatable :: merr
      real(wp) :: r2
      integer :: k, i

      do k = 1, size(npts_list)
         call new_s2_grid_lebedev(grid, merr, npts=npts_list(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if

         call check(error, size(grid%points, 1), 3, "S2 points are not 3-vectors")
         if (allocated(error)) return

         do i = 1, grid%npts
            r2 = sum(grid%points(:, i)**2)
            call check(error, abs(r2 - 1.0_wp) < struct_tol, &
               & "S2 node is not a unit vector")
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
      type(moist_math_grid_s2_lebedev_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_s2_grid_lebedev(grid, merr)
      call check(error, allocated(merr), "Missing npts=/degree= was accepted")
      if (allocated(error)) return
      deallocate (merr)

      call new_s2_grid_lebedev(grid, merr, npts=74, degree=13)
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
      type(moist_math_grid_s2_lebedev_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_s2_grid_lebedev(grid, merr, npts=7)
      call check(error, allocated(merr), "Unsupported Lebedev size was accepted")
      if (allocated(error)) return
      call check(error, grid%npts, 0, "Rejected grid was left initialized")
   end subroutine test_unsupported_npts

   !> degree= selects the smallest rule reaching the requested exactness
   !>
   !> %degree() reports a value at least as large as the request
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
      type(moist_math_grid_s2_lebedev_type), target :: grid
      class(moist_math_grid_s2_type), pointer :: base
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, size(request)
         call new_s2_grid_lebedev(grid, merr, degree=request(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         base => grid

         call check(error, grid%npts, expect_npts(k), &
            & "degree= did not select the smallest sufficient rule")
         if (allocated(error)) return

         call check(error, base%degree(), expect_degree(k), &
            & "Unexpected reported exactness degree")
         if (allocated(error)) return

         call check(error, base%degree() >= request(k), &
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
      type(moist_math_grid_s2_lebedev_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_s2_grid_lebedev(grid, merr, degree=132)
      call check(error, allocated(merr), "Degree beyond the highest rule was accepted")
      if (allocated(error)) return
      deallocate (merr)

      call new_s2_grid_lebedev(grid, merr, degree=-1)
      call check(error, allocated(merr), "Negative degree was accepted")
   end subroutine test_degree_out_of_range

   !> A rule of algebraic degree d reproduces exact sphere averages of polynomials
   !>
   !> In the Cartesian components, up to that degree; with weights summing
   !> to 1 the quadrature is the average, so
   !>   <1>=1, <x>=<x*y>=0, <3z^2-1>=0, <x^2>=1/3,
   !>   <x^4+y^4+z^4>=3/5
   !>
   !> @param[out] error  Test failure
   subroutine test_polynomial_exactness(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_s2_lebedev_type), target :: grid
      class(moist_math_grid_s2_type), pointer :: base
      type(mctc_error), allocatable :: merr
      integer :: k

      do k = 1, size(npts_list)
         call new_s2_grid_lebedev(grid, merr, npts=npts_list(k))
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         base => grid

         call check_average(error, base, f_one, 1.0_wp, "<1>")
         if (allocated(error)) return
         call check_average(error, base, f_x, 0.0_wp, "<x>")
         if (allocated(error)) return
         call check_average(error, base, f_xy, 0.0_wp, "<x*y>")
         if (allocated(error)) return
         call check_average(error, base, f_3z2m1, 0.0_wp, "<3z^2-1>")
         if (allocated(error)) return
         call check_average(error, base, f_x2, 1.0_wp/3.0_wp, "<x^2>")
         if (allocated(error)) return

         ! Degree 4; the 6-point rule is only exact to degree 3
         if (base%degree() >= 4) then
            call check_average(error, base, f_quartic, 0.6_wp, "<x^4+y^4+z^4>")
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
   !>   exact sphere average `1/(n+1)`
   !> - Relative error; measured at most 2.8e-15 over all rules
   !>
   !> @param[out] error  Test failure
   subroutine test_full_degree_exactness(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of supported Lebedev rules
      integer, parameter :: nrules = 32
      !> Acceptance bound on the relative error
      real(wp), parameter :: thr = 2.0e-14_wp
      type(moist_math_grid_s2_lebedev_type), target :: grid
      class(moist_math_grid_s2_type), pointer :: base
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
         call new_s2_grid_lebedev(grid, merr, degree=request)
         if (allocated(merr)) exit
         base => grid
         visited = visited + 1

         n = base%degree() - 1
         allocate (samples(grid%npts))
         do i = 1, grid%npts
            samples(i) = dot_product(dir, grid%points(:, i))**n
         end do
         call base%integrate_field(samples, average)
         deallocate (samples)

         exact = 1.0_wp/real(n + 1, wp)
         write (msg, "(a,i0,a,i0)") "Lebedev rule with ", grid%npts, &
            & " points is not exact at degree ", n
         call check(error, abs(average - exact) <= thr*exact, trim(msg))
         if (allocated(error)) return
         request = base%degree() + 1
      end do

      call check(error, visited, nrules, "degree= selection did not visit every Lebedev rule")
   end subroutine test_full_degree_exactness

   !> `integrate` and `integrate_field` must agree to round-off
   !>
   !> Analytic integrand vs. tabulated samples, on the same function
   !>
   !> @param[out] error  Test failure
   subroutine test_integrate_agreement(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_s2_lebedev_type), target :: grid
      class(moist_math_grid_s2_type), pointer :: base
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: samples(:)
      real(wp) :: from_field, from_analytic
      integer :: i

      call new_s2_grid_lebedev(grid, merr, npts=230)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      base => grid

      allocate (samples(grid%npts))
      do i = 1, grid%npts
         samples(i) = f_mixed(grid%points(:, i))
      end do

      call base%integrate_field(samples, from_field)
      call base%integrate(f_mixed, from_analytic)

      call check(error, abs(from_field - from_analytic) < struct_tol, &
         & "integrate and integrate_field disagree")
   end subroutine test_integrate_agreement

   !> Lebedev is an integration-only scheme, with no spherical harmonic trafo
   !>
   !> The inherited default `new_trafo` must report an error and leave
   !> `trafo` unallocated
   !>
   !> @param[out] error  Test failure
   subroutine test_no_trafo(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_s2_lebedev_type), target :: grid
      class(moist_math_grid_s2_type), pointer :: base
      class(moist_math_grid_s2_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr

      call new_s2_grid_lebedev(grid, merr, npts=74)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      base => grid

      call base%new_trafo(trafo, merr)
      call check(error, allocated(merr), "Lebedev grid claimed to have a transform")
      if (allocated(error)) return
      call check(error, .not. allocated(trafo), "new_trafo allocated a trafo on error")
   end subroutine test_no_trafo

   !> `destroy` must be idempotent and reset the grid to its empty state
   !>
   !> @param[out] error  Test failure
   subroutine test_destroy_idempotent(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_s2_lebedev_type) :: grid
      type(mctc_error), allocatable :: merr

      call new_s2_grid_lebedev(grid, merr, npts=110)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if

      call grid%destroy()
      call grid%destroy()

      call check(error, grid%npts, 0, "destroy did not reset npts")
      if (allocated(error)) return
      call check(error, grid%degree(), 0, "destroy did not reset the degree")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%points), "destroy left points allocated")
      if (allocated(error)) return
      call check(error, .not. allocated(grid%weights), "destroy left weights allocated")
   end subroutine test_destroy_idempotent

   !> Compare the grid average of an integrand against its exact value
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
      class(moist_math_grid_s2_type), intent(in) :: grid
      !> Angular integrand
      procedure(integrand_s2) :: f
      !> Exact sphere average of f
      real(wp), intent(in) :: expected
      !> Name of the integrand
      character(len=*), intent(in) :: label

      real(wp) :: val

      call grid%integrate(f, val)
      call check(error, abs(val - expected) < quad_tol, &
         & "S2 quadrature is not exact for "//label)
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

end module test_math_grid_s2
