!> Concrete single-center atomic grid
!>
!> - Atom-relative points and `r^2 dr dOmega` volume weights for one element,
!>   built from an atomic recipe; no molecule involved
!> - Shells in the radial grid's order, angular points in table order within
!>   each shell
!> - Shell radius and unit direction are kept separately, so a translated
!>   point can be formed as `R_A + r*u` in one expression
!> - Achieved angular count, degree, and spacing recorded per shell
module moist_math_grid_atomic_grid
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, new_radial_grid
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, &
      & moist_math_grid_angular_request_type
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type
   implicit none(type, external)
   private

   public :: moist_math_grid_atomic_type
   public :: new_atomic_grid
   public :: integrand_atomic

   !> Initial capacity of the distinct-request and distinct-rule stores
   integer, parameter :: initial_capacity = 8

   !> Scalar-valued integrand f(x) at an atom-relative point, used by `integrate`
   abstract interface
      !> Integrand at an atom-relative point
      !>
      !> @param[in] r  Atom-relative point (bohr)
      pure function integrand_atomic(r) result(val)
         import :: wp
         implicit none(type, external)
         !> Atom-relative point (bohr)
         real(wp), intent(in) :: r(3)
         !> Function value at r
         real(wp) :: val
      end function integrand_atomic
   end interface

   !> Single-center atomic grid
   !>
   !> Point i lies on shell `shell(i)` at radius `shell_r(shell(i))` in the
   !> unit direction `u(:, i)`; the points of shell s are
   !> `shell_offset(s)` to `shell_offset(s + 1) - 1`
   type :: moist_math_grid_atomic_type
      !> Atomic number the grid was built for
      integer :: z = 0
      !> Number of points
      integer :: npts = 0
      !> Number of retained radial shells
      integer :: nshell = 0
      !> Number of radial shells requested before cutoffs
      integer :: nshell_requested = 0
      !> Shell radii (bohr), shape (nshell), in the radial grid's order
      real(wp), allocatable :: shell_r(:)
      !> Radial `dr` weights (bohr), shape (nshell), without r^2 or 4*pi
      real(wp), allocatable :: shell_w(:)
      !> Achieved angular point count per shell, shape (nshell)
      integer, allocatable :: shell_npts(:)
      !> Achieved angular exactness degree per shell, shape (nshell)
      integer, allocatable :: shell_degree(:)
      !> Area-based spacing estimate `sqrt(4*pi*r**2/N)` per shell (bohr), shape (nshell)
      real(wp), allocatable :: shell_spacing(:)
      !> First point index of each shell, shape (nshell + 1); last entry is npts + 1
      integer, allocatable :: shell_offset(:)
      !> Shell index of each point, shape (npts)
      integer, allocatable :: shell(:)
      !> Unit direction of each point, shape (3, npts)
      real(wp), allocatable :: u(:, :)
      !> Atom-relative points `r*u` (bohr), shape (3, npts)
      real(wp), allocatable :: xyz(:, :)
      !> Volume weights `r**2*w_rad*w_ang` (bohr^3), shape (npts)
      real(wp), allocatable :: w(:)
   contains
      !> Volume quadrature of a field already tabulated on the points
      procedure :: integrate_field => atomic_integrate_field
      !> Volume quadrature of an analytic integrand sampled at the local points
      procedure :: integrate => atomic_integrate
      !> Release all storage (idempotent)
      procedure :: destroy => atomic_destroy
   end type moist_math_grid_atomic_type

contains

   !> Build the atomic grid of element z from a recipe
   !>
   !> - Radial grid from the radial recipe (after its cutoffs)
   !> - One angular request per shell from the shell policy, served by the
   !>   angular generator; identical requests are selected once and
   !>   identical rules are stored once
   !> - `w = r**2*w_rad*w_ang`; the angular weights already carry 4*pi
   !>
   !> @param[out] grid    New atomic grid; empty on error
   !> @param[in]  recipe  Atomic recipe
   !> @param[in]  z       Atomic number
   !> @param[out] error   Set on an incomplete or invalid recipe, a radial or
   !>                     angular failure, or an unrepresentable point count
   subroutine new_atomic_grid(grid, recipe, z, error)
      !> New atomic grid
      type(moist_math_grid_atomic_type), intent(out) :: grid
      !> Atomic recipe
      type(moist_math_grid_atomic_recipe_type), intent(in) :: recipe
      !> Atomic number
      integer, intent(in) :: z
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Radial grid of element z after cutoffs
      type(moist_math_grid_radial_type) :: radial
      !> Distinct angular rules, shape (nrule)
      type(moist_math_grid_angular_type), allocatable :: rules(:)
      !> Rule index of each shell, shape (nshell)
      integer, allocatable :: shell_rule(:)
      integer :: nshell, nrule, ish, ip, j, n, stat
      integer(int64) :: total
      real(wp) :: r, wr

      if (.not. allocated(recipe%angular)) then
         call fatal_error(error, "Atomic grid: recipe has no angular generator")
         return
      end if
      if (.not. allocated(recipe%shells)) then
         call fatal_error(error, "Atomic grid: recipe has no shell policy")
         return
      end if
      call recipe%shells%validate(error)
      if (allocated(error)) return

      call new_radial_grid(radial, recipe%radial, z, error)
      if (allocated(error)) return
      nshell = radial%npts

      call select_shell_rules(recipe, z, radial%r, rules, nrule, shell_rule, error)
      if (allocated(error)) return

      total = 0_int64
      do ish = 1, nshell
         total = total + int(rules(shell_rule(ish))%npts, int64)
      end do
      if (total > int(huge(0), int64)) then
         call fatal_error(error, "Atomic grid: point count exceeds the representable integer range")
         return
      end if
      n = int(total)

      allocate (grid%shell_r(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_w(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_npts(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_degree(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_spacing(nshell), stat=stat)
      if (stat == 0) allocate (grid%shell_offset(nshell + 1), stat=stat)
      if (stat == 0) allocate (grid%shell(n), stat=stat)
      if (stat == 0) allocate (grid%u(3, n), stat=stat)
      if (stat == 0) allocate (grid%xyz(3, n), stat=stat)
      if (stat == 0) allocate (grid%w(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Atomic grid: failed to allocate the grid points")
         call grid%destroy()
         return
      end if

      ip = 0
      grid%shell_offset(1) = 1
      do ish = 1, nshell
         associate (rule => rules(shell_rule(ish)))
            r = radial%r(ish)
            wr = radial%w(ish)
            grid%shell_r(ish) = r
            grid%shell_w(ish) = wr
            grid%shell_npts(ish) = rule%npts
            grid%shell_degree(ish) = rule%degree
            grid%shell_spacing(ish) = sqrt(4.0_wp*pi*r*r/real(rule%npts, wp))
            do j = 1, rule%npts
               ip = ip + 1
               grid%shell(ip) = ish
               grid%u(:, ip) = rule%points(:, j)
               grid%xyz(1, ip) = r*rule%points(1, j)
               grid%xyz(2, ip) = r*rule%points(2, j)
               grid%xyz(3, ip) = r*rule%points(3, j)
               grid%w(ip) = r*r*wr*rule%weights(j)
            end do
         end associate
         grid%shell_offset(ish + 1) = ip + 1
      end do

      grid%z = z
      grid%npts = n
      grid%nshell = nshell
      grid%nshell_requested = radial%npts_requested
   end subroutine new_atomic_grid

   !> Select the angular rule of every shell
   !>
   !> Requests equal in every field reuse the earlier selection, and rules
   !> equal in size, points, and weights are stored once
   !>
   !> @param[in]  recipe      Atomic recipe with generator and shell policy
   !> @param[in]  z           Atomic number
   !> @param[in]  radii       Shell radii (bohr), shape (nshell)
   !> @param[out] rules       Distinct rules, shape (nrule) or larger
   !> @param[out] nrule       Number of distinct rules
   !> @param[out] shell_rule  Rule index of each shell, shape (nshell)
   !> @param[out] error       Set if the generator cannot serve a request
   subroutine select_shell_rules(recipe, z, radii, rules, nrule, shell_rule, error)
      !> Atomic recipe
      type(moist_math_grid_atomic_recipe_type), intent(in) :: recipe
      !> Atomic number
      integer, intent(in) :: z
      !> Shell radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Distinct rules
      type(moist_math_grid_angular_type), allocatable, intent(out) :: rules(:)
      !> Number of distinct rules
      integer, intent(out) :: nrule
      !> Rule index of each shell
      integer, allocatable, intent(out) :: shell_rule(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Distinct requests seen so far, shape (nseen) or larger
      type(moist_math_grid_angular_request_type), allocatable :: seen(:)
      !> Rule index of each distinct request
      integer, allocatable :: seen_rule(:)
      type(moist_math_grid_angular_request_type) :: request
      type(moist_math_grid_angular_type) :: candidate
      type(error_type), allocatable :: select_error
      character(len=128) :: msg
      integer :: nshell, nseen, ish, k, irule

      nshell = size(radii)
      allocate (shell_rule(nshell))
      ! Distinct requests and rules are few, so their stores grow on demand;
      ! sized per shell they would cost one rule per shell before the
      ! caller's point-count guard
      allocate (rules(initial_capacity), seen(initial_capacity), seen_rule(initial_capacity))
      nrule = 0
      nseen = 0
      do ish = 1, nshell
         request = recipe%shells%request(z, radii(ish))

         irule = 0
         do k = 1, nseen
            if (same_request(seen(k), request)) then
               irule = seen_rule(k)
               exit
            end if
         end do

         if (irule == 0) then
            call recipe%angular%select(request, candidate, select_error)
            if (allocated(select_error)) then
               write (msg, "(a,i0,a,i0,a,es12.5,a)") "Atomic grid for z = ", z, ", shell ", ish, &
                  & " (r = ", radii(ish), " bohr): "
               call fatal_error(error, trim(msg)//" "//select_error%message)
               return
            end if
            do k = 1, nrule
               if (same_rule(rules(k), candidate)) then
                  irule = k
                  exit
               end if
            end do
            if (irule == 0) then
               if (nrule == size(rules)) call grow_rules(rules, nrule)
               nrule = nrule + 1
               irule = nrule
               rules(irule) = candidate
            end if
            if (nseen == size(seen)) call grow_requests(seen, seen_rule, nseen)
            nseen = nseen + 1
            seen(nseen) = request
            seen_rule(nseen) = irule
         end if

         shell_rule(ish) = irule
      end do
   end subroutine select_shell_rules

   !> Double the capacity of the distinct-rule store
   !>
   !> @param[in,out] rules  Distinct rules; the first n entries are kept
   !> @param[in]     n      Number of entries in use
   subroutine grow_rules(rules, n)
      !> Distinct rules
      type(moist_math_grid_angular_type), allocatable, intent(inout) :: rules(:)
      !> Number of entries in use
      integer, intent(in) :: n

      type(moist_math_grid_angular_type), allocatable :: grown(:)
      integer :: k

      allocate (grown(2*size(rules)))
      do k = 1, n
         grown(k) = rules(k)
      end do
      call move_alloc(grown, rules)
   end subroutine grow_rules

   !> Double the capacity of the distinct-request store
   !>
   !> @param[in,out] seen       Distinct requests; the first n entries are kept
   !> @param[in,out] seen_rule  Rule index of each distinct request
   !> @param[in]     n          Number of entries in use
   subroutine grow_requests(seen, seen_rule, n)
      !> Distinct requests
      type(moist_math_grid_angular_request_type), allocatable, intent(inout) :: seen(:)
      !> Rule index of each distinct request
      integer, allocatable, intent(inout) :: seen_rule(:)
      !> Number of entries in use
      integer, intent(in) :: n

      type(moist_math_grid_angular_request_type), allocatable :: grown(:)
      integer, allocatable :: grown_rule(:)

      allocate (grown(2*size(seen)), grown_rule(2*size(seen)))
      grown(1:n) = seen(1:n)
      grown_rule(1:n) = seen_rule(1:n)
      call move_alloc(grown, seen)
      call move_alloc(grown_rule, seen_rule)
   end subroutine grow_requests

   !> Two angular requests with equal fields
   !>
   !> @param[in] a  First request
   !> @param[in] b  Second request
   pure function same_request(a, b) result(same)
      !> First request
      type(moist_math_grid_angular_request_type), intent(in) :: a
      !> Second request
      type(moist_math_grid_angular_request_type), intent(in) :: b
      !> True if every field is equal
      logical :: same

      same = a%npts == b%npts .and. a%min_degree == b%min_degree &
         & .and. a%min_points == b%min_points .and. a%max_points == b%max_points &
         & .and. a%target_points == b%target_points
   end function same_request

   !> Two angular rules with equal size, degree, points, and weights
   !>
   !> @param[in] a  First rule
   !> @param[in] b  Second rule
   pure function same_rule(a, b) result(same)
      !> First rule
      type(moist_math_grid_angular_type), intent(in) :: a
      !> Second rule
      type(moist_math_grid_angular_type), intent(in) :: b
      !> True if both rules hold the same nodes and weights
      logical :: same

      same = a%npts == b%npts .and. a%degree == b%degree
      if (.not. same) return
      same = all(a%points == b%points) .and. all(a%weights == b%weights)
   end function same_rule

   !> Volume quadrature of tabulated values: `result = sum_i w(i)*f(i)`
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-point field values, shape (npts); any other length is an error
   !> @param[out] result  Quadrature result
   !> @param[out] error   Error handling
   subroutine atomic_integrate_field(self, f, result, error)
      !> Grid instance
      class(moist_math_grid_atomic_type), intent(in) :: self
      !> Per-point field values
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i

      result = 0.0_wp
      if (size(f) /= self%npts) then
         call fatal_error(error, "atomic grid: integrate_field needs one value per point")
         return
      end if
      do i = 1, self%npts
         result = result + self%w(i)*f(i)
      end do
   end subroutine atomic_integrate_field

   !> Volume quadrature of an analytic integrand: `result = sum_i w(i)*f(xyz(:, i))`
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Integrand at an atom-relative point
   !> @param[out] result  Quadrature result
   subroutine atomic_integrate(self, f, result)
      !> Grid instance
      class(moist_math_grid_atomic_type), intent(in) :: self
      !> Integrand at an atom-relative point
      procedure(integrand_atomic) :: f
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + self%w(i)*f(self%xyz(:, i))
      end do
   end subroutine atomic_integrate

   !> Release all grid storage; idempotent
   !>
   !> @param[in,out] self  Grid instance
   pure subroutine atomic_destroy(self)
      !> Grid instance
      class(moist_math_grid_atomic_type), intent(inout) :: self

      if (allocated(self%shell_r)) deallocate (self%shell_r)
      if (allocated(self%shell_w)) deallocate (self%shell_w)
      if (allocated(self%shell_npts)) deallocate (self%shell_npts)
      if (allocated(self%shell_degree)) deallocate (self%shell_degree)
      if (allocated(self%shell_spacing)) deallocate (self%shell_spacing)
      if (allocated(self%shell_offset)) deallocate (self%shell_offset)
      if (allocated(self%shell)) deallocate (self%shell)
      if (allocated(self%u)) deallocate (self%u)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%w)) deallocate (self%w)
      self%z = 0
      self%npts = 0
      self%nshell = 0
      self%nshell_requested = 0
   end subroutine atomic_destroy

end module moist_math_grid_atomic_grid
