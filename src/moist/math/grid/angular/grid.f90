!> Angular grids on the unit sphere and the generators that build them
!>
!> - Angular grid: unit directions and solid-angle (`dOmega`) weights summing
!>   to 4*pi; records the algebraic exactness degree the rule achieves; one
!>   type for every angular scheme, built by an angular generator
!> - Angular request and abstract generator: a generator selects and builds an
!>   admissible angular rule for a request; the scheme-specific selection
!>   lives in each concrete generator
!> - Lebedev-Laikov generator: selects among the 32 Lebedev-Laikov rules for
!>   an angular request, optionally excluding the rules with negative weights;
!>   weights scaled from the Laikov tables (sum 1) to dOmega (sum 4*pi)
module moist_math_grid_angular_grid
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use moist_math_grid_angular_lebedev, only: grid_size, lebedev_degree_table, &
      & lebedev_has_negative_weights, lebedev_order_from_num, get_angular_grid
   implicit none(type, external)
   private

   public :: moist_math_grid_angular_type
   public :: integrand_angular
   public :: moist_math_grid_angular_generator_type
   public :: moist_math_grid_angular_request_type
   public :: validate_angular_request
   public :: moist_math_grid_angular_generator_lebedev_type
   public :: new_lebedev_generator
   public :: new_lebedev_grid

   !> Scalar-valued angular integrand f(u) on the unit sphere, used by `integrate`
   abstract interface
      !> Angular integrand f(u) on the unit sphere
      !>
      !> @param[in] uvec  Cartesian unit vector on the unit sphere
      pure function integrand_angular(uvec) result(val)
         import :: wp
         implicit none(type, external)
         !> Cartesian unit vector on the unit sphere
         real(wp), intent(in) :: uvec(3)
         !> Function value at uvec
         real(wp) :: val
      end function integrand_angular
   end interface

   !> Unit-sphere grid: directions, solid-angle weights, and exactness degree
   type :: moist_math_grid_angular_type
      !> Number of angular nodes
      integer :: npts = 0
      !> Cartesian unit vectors of the nodes, shape (3, npts)
      real(wp), allocatable :: points(:, :)
      !> Solid-angle (dOmega) weights, shape (npts), summing to 4*pi
      real(wp), allocatable :: weights(:)
      !> Algebraic exactness degree achieved by the rule; 0 means none
      integer :: degree = 0
   contains
      !> Quadrature of a field already tabulated on the nodes
      procedure :: integrate_field => angular_integrate_field
      !> Quadrature of an analytic angular integrand sampled at the nodes
      procedure :: integrate => angular_integrate
      !> Release all storage (idempotent)
      procedure :: destroy => angular_destroy
   end type moist_math_grid_angular_type

   !> Angular resolution request
   !>
   !> - `npts > 0`: exact point count; the rule must exist, no substitution
   !> - Otherwise the smallest admissible rule meeting `target_points`, or the
   !>   largest admissible one within `max_points` if the target is unreachable
   !> - `min_degree`, `min_points`, `max_points` are hard constraints
   type :: moist_math_grid_angular_request_type
      !> Exact point count; 0 = none
      integer :: npts = 0
      !> Minimum algebraic exactness degree (hard)
      integer :: min_degree = 0
      !> Minimum point count (hard floor)
      integer :: min_points = 0
      !> Maximum point count (hard cap)
      integer :: max_points = huge(0)
      !> Desired point count (soft); 0 = none
      real(wp) :: target_points = 0.0_wp
   end type moist_math_grid_angular_request_type

   !> Abstract angular generator: select and build a rule for a request
   type, abstract :: moist_math_grid_angular_generator_type
   contains
      !> Select the rule satisfying a request and build its grid
      procedure(angular_select_i), deferred :: select
   end type moist_math_grid_angular_generator_type

   !> Deferred selection shared by every angular generator
   abstract interface
      !> Select the rule satisfying a request and build its grid
      !>
      !> @param[in]  self     Generator instance
      !> @param[in]  request  Angular resolution request
      !> @param[out] grid     Generated grid; empty on error
      !> @param[out] error    Set if no rule satisfies the request
      subroutine angular_select_i(self, request, grid, error)
         import :: moist_math_grid_angular_generator_type, moist_math_grid_angular_request_type, &
            & moist_math_grid_angular_type, error_type
         implicit none(type, external)
         !> Generator instance
         class(moist_math_grid_angular_generator_type), intent(in) :: self
         !> Angular resolution request
         type(moist_math_grid_angular_request_type), intent(in) :: request
         !> Generated grid; empty on error
         type(moist_math_grid_angular_type), intent(out) :: grid
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine angular_select_i
   end interface

   !> Lebedev-Laikov angular generator
   type, extends(moist_math_grid_angular_generator_type) :: moist_math_grid_angular_generator_lebedev_type
      !> Exclude the rules with negative weights (74, 230, 266 points)
      logical :: positive_weights_only = .false.
   contains
      !> Select the rule satisfying a request and build its grid
      procedure :: select => lebedev_select
   end type moist_math_grid_angular_generator_lebedev_type

contains

   !> Field quadrature: `result = sum_i weights(i)*f(i)`
   !>
   !> Solid-angle integral of f over the unit sphere; divide by 4*pi for the
   !> average
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-node field values, shape (npts); any other length is an error
   !> @param[out] result  Integral of f over the unit sphere
   !> @param[out] error   Error handling
   subroutine angular_integrate_field(self, f, result, error)
      !> Grid instance
      class(moist_math_grid_angular_type), intent(in) :: self
      !> Per-node field values, shape (npts)
      real(wp), intent(in) :: f(:)
      !> Integral of f over the unit sphere
      real(wp), intent(out) :: result
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Node index
      integer :: i

      result = 0.0_wp
      if (size(f) /= self%npts) then
         call fatal_error(error, "angular grid: integrate_field needs one value per node")
         return
      end if
      do i = 1, self%npts
         result = result + self%weights(i)*f(i)
      end do
   end subroutine angular_integrate_field

   !> Analytic quadrature: `result = sum_i weights(i)*f(u_i)`
   !>
   !> Samples the integrand at the unit vectors the grid owns
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Angular integrand f(u)
   !> @param[out] result  Integral of f over the unit sphere
   subroutine angular_integrate(self, f, result)
      !> Grid instance
      class(moist_math_grid_angular_type), intent(in) :: self
      !> Angular integrand f(u)
      procedure(integrand_angular) :: f
      !> Integral of f over the unit sphere
      real(wp), intent(out) :: result

      !> Node index
      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + self%weights(i)*f(self%points(:, i))
      end do
   end subroutine angular_integrate

   !> Release all grid storage; idempotent
   !>
   !> @param[in,out] self  Grid instance
   pure subroutine angular_destroy(self)
      !> Grid instance
      class(moist_math_grid_angular_type), intent(inout) :: self

      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%weights)) deallocate (self%weights)
      self%npts = 0
      self%degree = 0
   end subroutine angular_destroy

   !> Reject a request whose fields are out of range
   !>
   !> - Counts and degree must be non-negative
   !> - `target_points` must be non-negative and not NaN; +inf is an
   !>   unreachable target
   !> - Mutually incompatible hard constraints are left to the generator
   !>
   !> @param[in]  request  Angular resolution request
   !> @param[out] error    Set for an out-of-range field
   subroutine validate_angular_request(request, error)
      !> Angular resolution request
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (request%npts < 0) then
         call fatal_error(error, "Angular request: npts must be non-negative")
         return
      end if
      if (request%min_degree < 0) then
         call fatal_error(error, "Angular request: min_degree must be non-negative")
         return
      end if
      if (request%min_points < 0) then
         call fatal_error(error, "Angular request: min_points must be non-negative")
         return
      end if
      if (request%max_points < 0) then
         call fatal_error(error, "Angular request: max_points must be non-negative")
         return
      end if
      ! Tested apart so a NaN is never compared
      if (ieee_is_nan(request%target_points)) then
         call fatal_error(error, "Angular request: target_points must not be NaN")
         return
      end if
      if (request%target_points < 0.0_wp) then
         call fatal_error(error, "Angular request: target_points must be non-negative")
         return
      end if
   end subroutine validate_angular_request

   !> Create a Lebedev-Laikov angular generator
   !>
   !> @param[out] self                   Generator instance
   !> @param[in]  positive_weights_only  Exclude rules with negative weights, default false
   pure subroutine new_lebedev_generator(self, positive_weights_only)
      !> Generator instance
      type(moist_math_grid_angular_generator_lebedev_type), intent(out) :: self
      !> Exclude rules with negative weights
      logical, intent(in), optional :: positive_weights_only

      if (present(positive_weights_only)) self%positive_weights_only = positive_weights_only
   end subroutine new_lebedev_generator

   !> Lebedev-Laikov grid for an exact point count or a minimum degree
   !>
   !> Exactly one of `npts` and `degree` must be present
   !>
   !> @param[out] grid                   Generated grid; empty on error
   !> @param[out] error                  Set for an invalid or unsatisfiable request
   !> @param[in]  npts                   Exact point count, one of the supported sizes
   !> @param[in]  degree                 Minimum exactness degree; the smallest admissible rule is used
   !> @param[in]  positive_weights_only  Exclude rules with negative weights, default false
   subroutine new_lebedev_grid(grid, error, npts, degree, positive_weights_only)
      !> Generated grid; empty on error
      type(moist_math_grid_angular_type), intent(out) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Exact point count
      integer, intent(in), optional :: npts
      !> Minimum exactness degree
      integer, intent(in), optional :: degree
      !> Exclude rules with negative weights
      logical, intent(in), optional :: positive_weights_only

      !> Generator with the requested filter
      type(moist_math_grid_angular_generator_lebedev_type) :: generator
      !> Request built from the arguments
      type(moist_math_grid_angular_request_type) :: request

      if (present(npts) .eqv. present(degree)) then
         call fatal_error(error, "Lebedev grid needs exactly one of npts and degree")
         return
      end if

      if (present(npts)) then
         ! npts = 0 would mean "no exact count" in the request
         if (npts <= 0) then
            call fatal_error(error, "Lebedev grid: npts must be positive")
            return
         end if
         request%npts = npts
      else
         request%min_degree = degree
      end if

      call new_lebedev_generator(generator, positive_weights_only)
      call generator%select(request, grid, error)
   end subroutine new_lebedev_grid

   !> Select the Lebedev rule satisfying a request and build its grid
   !>
   !> - `npts > 0`: that exact rule; unsupported sizes, filtered sizes, and
   !>   violated hard constraints are errors, never substituted
   !> - Otherwise among the admissible rules (filter, `min_degree`,
   !>   `min_points`, `max_points`): the smallest with at least
   !>   `target_points` points, else the largest admissible one
   !> - No admissible rule is an error; the degree is never relaxed
   !>
   !> @param[in]  self     Generator instance
   !> @param[in]  request  Angular resolution request
   !> @param[out] grid     Generated grid; empty on error
   !> @param[out] error    Set if no rule satisfies the request
   subroutine lebedev_select(self, request, grid, error)
      !> Generator instance
      class(moist_math_grid_angular_generator_lebedev_type), intent(in) :: self
      !> Angular resolution request
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Generated grid; empty on error
      type(moist_math_grid_angular_type), intent(out) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Order index of the selected rule
      integer :: order

      call validate_angular_request(request, error)
      if (allocated(error)) return

      if (request%npts > 0) then
         call lebedev_exact_order(self, request, order, error)
      else
         call lebedev_select_order(self, request, order, error)
      end if
      if (allocated(error)) return

      call lebedev_build(order, grid, error)
   end subroutine lebedev_select

   !> Order index of an exact point-count request
   !>
   !> @param[in]  self     Generator instance
   !> @param[in]  request  Request with `npts > 0`
   !> @param[out] order    Order index into `grid_size`, 0 on error
   !> @param[out] error    Set for an unsupported, filtered, or constrained-out size
   subroutine lebedev_exact_order(self, request, order, error)
      !> Generator instance
      class(moist_math_grid_angular_generator_lebedev_type), intent(in) :: self
      !> Request with `npts > 0`
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Order index into `grid_size`
      integer, intent(out) :: order
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Error message buffer
      character(len=256) :: msg

      call lebedev_order_from_num(request%npts, order, error, &
         & positive_weights_only=self%positive_weights_only)
      if (allocated(error)) return

      if (lebedev_degree_table(order) < request%min_degree &
         & .or. grid_size(order) < request%min_points &
         & .or. grid_size(order) > request%max_points) then
         write (msg, "(a,i0,a,i0,a,i0,a,i0,a,i0)") "Lebedev rule with ", grid_size(order), &
            & " points (degree ", lebedev_degree_table(order), &
            & ") violates the hard constraints min_degree=", request%min_degree, &
            & ", min_points=", request%min_points, ", max_points=", request%max_points
         order = 0
         call fatal_error(error, trim(msg))
         return
      end if
   end subroutine lebedev_exact_order

   !> Order index selected by degree, count bounds, and target
   !>
   !> Walks the rules in ascending size, so the first admissible rule meeting
   !> the target is the smallest, and the last admissible one seen is the
   !> largest within the cap
   !>
   !> @param[in]  self     Generator instance
   !> @param[in]  request  Request with `npts = 0`
   !> @param[out] order    Order index into `grid_size`, 0 on error
   !> @param[out] error    Set if no rule satisfies the hard constraints
   subroutine lebedev_select_order(self, request, order, error)
      !> Generator instance
      class(moist_math_grid_angular_generator_lebedev_type), intent(in) :: self
      !> Request with `npts = 0`
      type(moist_math_grid_angular_request_type), intent(in) :: request
      !> Order index into `grid_size`
      integer, intent(out) :: order
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Loop index over the rules
      integer :: i
      !> Error message buffer
      character(len=256) :: msg

      order = 0
      do i = 1, size(grid_size)
         if (self%positive_weights_only .and. lebedev_has_negative_weights(i)) cycle
         if (lebedev_degree_table(i) < request%min_degree) cycle
         if (grid_size(i) < request%min_points) cycle
         if (grid_size(i) > request%max_points) exit
         order = i
         if (real(grid_size(i), wp) >= request%target_points) return
      end do

      if (order == 0) then
         write (msg, "(a,i0,a,i0,a,i0,a,l1)") &
            & "No Lebedev rule satisfies the hard constraints min_degree=", request%min_degree, &
            & ", min_points=", request%min_points, ", max_points=", request%max_points, &
            & ", positive_weights_only=", self%positive_weights_only
         call fatal_error(error, trim(msg))
      end if
   end subroutine lebedev_select_order

   !> Build the angular grid of one Lebedev rule
   !>
   !> Weights are `w_table*(4.0_wp*pi)`, the Laikov table scaled to dOmega
   !>
   !> @param[in]  order  Order index into `grid_size`
   !> @param[out] grid   Generated grid; empty on error
   !> @param[out] error  Set if the table generation fails
   subroutine lebedev_build(order, grid, error)
      !> Order index into `grid_size`
      integer, intent(in) :: order
      !> Generated grid; empty on error
      type(moist_math_grid_angular_type), intent(out) :: grid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Laikov-normalized table weights, summing to 1
      real(wp), allocatable :: w_table(:)
      !> Number of nodes of the rule
      integer :: n

      n = grid_size(order)
      allocate (grid%points(3, n), w_table(n))
      call get_angular_grid(order, grid%points, w_table, error)
      if (allocated(error)) then
         call grid%destroy()
         return
      end if

      grid%weights = w_table*(4.0_wp*pi)
      grid%npts = n
      grid%degree = lebedev_degree_table(order)
   end subroutine lebedev_build

end module moist_math_grid_angular_grid
