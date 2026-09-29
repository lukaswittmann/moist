!> Abstract base types for moist S2 (unit sphere) angular grids
module moist_math_grid_s2_base
   use mctc_env, only: wp, error_type, fatal_error
   implicit none(type, external)
   private

   public :: moist_math_grid_s2_type
   public :: moist_math_grid_s2_trafo_type
   public :: integrand_s2

   !> Scalar-valued angular integrand f(u) on the unit sphere, used by `integrate`
   abstract interface
      !> Angular integrand f(u) on the unit sphere
      !>
      !> @param[in] uvec  Cartesian unit vector on S2
      pure function integrand_s2(uvec) result(val)
         import :: wp
         implicit none(type, external)
         !> Cartesian unit vector on S2
         real(wp), intent(in) :: uvec(3)
         !> Function value at uvec
         real(wp) :: val
      end function integrand_s2
   end interface

   !> Abstract immutable unit-sphere grid: nodes, weights, and S2 quadrature
   type, abstract :: moist_math_grid_s2_type
      !> Number of angular nodes
      integer :: npts = 0
      !> Cartesian unit vectors of the nodes, shape (3, npts)
      real(wp), allocatable :: points(:, :)
      !> Quadrature weights, shape (npts), normalised to sum to 1
      real(wp), allocatable :: weights(:)
   contains
      !> Quadrature weight of node i (1..npts)
      procedure(s2_measure_i), deferred :: measure
      !> Release all storage (idempotent)
      procedure(s2_destroy_i), deferred :: destroy
      !> Algebraic exactness degree of the rule; 0 means "no exactness degree"
      procedure :: degree => s2_degree_default
      !> Quadrature of a field already tabulated on the nodes
      procedure :: integrate_field => s2_integrate_field_default
      !> Quadrature of an analytic angular integrand sampled at the nodes
      procedure :: integrate => s2_integrate_default
      !> Allocate the spherical-harmonic transform matching this grid
      procedure :: new_trafo => s2_new_trafo_default
   end type moist_math_grid_s2_type

   !> Abstract marker for a spherical-harmonic transform engine bound to an S2 grid
   type, abstract :: moist_math_grid_s2_trafo_type
   end type moist_math_grid_s2_trafo_type

   !> Deferred grid operations shared by every concrete S2 grid
   abstract interface
      !> Quadrature weight (fraction of the unit sphere) of node i
      !>
      !> @param[in] self  Grid instance
      !> @param[in] i     Grid node index
      pure function s2_measure_i(self, i) result(w)
         import :: moist_math_grid_s2_type, wp
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_s2_type), intent(in) :: self
         !> Grid node index
         integer, intent(in) :: i
         !> Weight of node i
         real(wp) :: w
      end function s2_measure_i

      !> Release all grid storage, idempotent
      !>
      !> Not pure: concrete grids may tear down external handles
      !>
      !> @param[in,out] self  Grid instance
      subroutine s2_destroy_i(self)
         import :: moist_math_grid_s2_type
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_s2_type), intent(inout) :: self
      end subroutine s2_destroy_i
   end interface

contains

   !> Default algebraic exactness degree
   !>
   !> 0, meaning the scheme has none (Fibonacci, quasi-Monte-Carlo, ...)
   !>
   !> @param[in]  self  Grid instance
   pure function s2_degree_default(self) result(d)
      !> Grid instance
      class(moist_math_grid_s2_type), intent(in) :: self
      !> Exactness degree, 0 if the scheme has none
      integer :: d

      d = 0
   end function s2_degree_default

   !> Default field quadrature: `result = sum_i measure(i) * f(i)`
   !>
   !> Consumes discrete field data already tabulated on the angular nodes
   !> With weights summing to 1 this is the *average* of f over the sphere;
   !> multiply by 4*pi for the solid-angle integral
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-node field values (length npts)
   !> @param[out] result  Quadrature result
   pure subroutine s2_integrate_field_default(self, f, result)
      !> Grid instance
      class(moist_math_grid_s2_type), intent(in) :: self
      !> Per-node field values (length npts)
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + self%measure(i)*f(i)
      end do
   end subroutine s2_integrate_field_default

   !> Default analytic quadrature: `result = sum_i measure(i) * f(u_i)`
   !>
   !> Samples the integrand at the unit vectors the grid owns
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Angular integrand f(u)
   !> @param[out] result  Quadrature result
   subroutine s2_integrate_default(self, f, result)
      !> Grid instance
      class(moist_math_grid_s2_type), intent(in) :: self
      !> Angular integrand f(u)
      procedure(integrand_s2) :: f
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + self%measure(i)*f(self%points(:, i))
      end do
   end subroutine s2_integrate_default

   !> Default transform factory: no spherical harmonic transform available
   !>
   !> Leaves `trafo` unallocated; pure integration schemes (Lebedev) inherit
   !> the default, sampling schemes override it; `self` has the target
   !> attribute so an overriding scheme may bind a pointer to its grid
   !>
   !> @param[in]  self   Grid instance
   !> @param[out] trafo  Left unallocated by this default
   !> @param[out] error  Always set by this default
   subroutine s2_new_trafo_default(self, trafo, error)
      !> Grid instance
      class(moist_math_grid_s2_type), intent(in), target :: self
      !> Left unallocated by this default
      class(moist_math_grid_s2_trafo_type), allocatable, intent(out) :: trafo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "This S2 grid has no spherical harmonic transform")
   end subroutine s2_new_trafo_default

end module moist_math_grid_s2_base
