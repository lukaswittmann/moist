!> Abstract base types for moist 1D radial grids
module moist_math_grid_1d_base
   use mctc_env, only: wp, error_type
   implicit none(type, external)
   private

   public :: moist_math_grid_1d_type
   public :: moist_math_grid_1d_trafo_type
   public :: integrand_radial

   !> Scalar-valued radial integrand f(r) used by `integrate`
   abstract interface
      !> Radial integrand f(r)
      !>
      !> @param[in] r  Radius (bohr)
      pure function integrand_radial(r) result(val)
         import :: wp
         implicit none(type, external)
         !> Radius (bohr)
         real(wp), intent(in) :: r
         !> Function value at r
         real(wp) :: val
      end function integrand_radial
   end interface

   !> Abstract immutable radial grid: nodes, weights, and 3D-radial quadrature
   type, abstract :: moist_math_grid_1d_type
      ! TODO: Move radial nodes, weights and reciprocal geometry into a radial domain
      ! Keep solver radial samples distinct from host atomic evaluation sites
      ! Defer the 1D ownership refactor; retain the current interface for now
      !> Number of radial points
      integer :: npts = 0
      !> r-space nodes (bohr)
      real(wp), allocatable :: r(:)
      !> k-space nodes (1/bohr)
      real(wp), allocatable :: k(:)
   contains
      !> 3D-radial quadrature weight (4*pi*r^2 dr) of grid point i
      procedure(radial_measure_i), deferred :: measure
      !> Allocate the transform engine matching this grid, bound to it
      procedure(radial_new_trafo_i), deferred :: new_trafo
      !> Release all storage (idempotent)
      procedure(radial_destroy_i), deferred :: destroy
      !> Quadrature of a field already tabulated on the nodes
      procedure :: integrate_field => radial_integrate_field_default
      !> Quadrature of an analytic radial integrand sampled at the nodes
      procedure :: integrate => radial_integrate_default
   end type moist_math_grid_1d_type

   !> Abstract per-thread radial Fourier-Bessel transform engine
   !>
   !> Bound to one grid (held by the concrete type as a typed pointer) and
   !> owns the mutable scratch a transform needs; not shared between threads,
   !> create one per thread; must not outlive the grid it is bound to
   !>
   !> Two call shapes are offered for each direction so the caller can pick the
   !> faster one case-by-case:
   !>
   !>   * scalar (`fbt_*`) -- transform a single field; loop / OpenMP across
   !>     several fields at the call site
   !>   * batched (`fbt_*_all`) -- transform every column of a 2D field in one
   !>     call; the base provides a default that simply loops the scalar form;
   !>     concrete trafos override it with a genuinely batched implementation
   !>     (a batched DST, a GEMM, ...) when one is available
   type, abstract :: moist_math_grid_1d_trafo_type
   contains
      !> Forward Fourier-Bessel transform (r -> k)
      procedure(radial_trafo_fbt_i), deferred :: fbt_r2k
      !> Backward Fourier-Bessel transform (k -> r)
      procedure(radial_trafo_fbt_i), deferred :: fbt_k2r
      !> Euclidean adjoint (transpose) of the forward transform (k -> r)
      procedure(radial_trafo_fbt_i), deferred :: fbt_r2k_adj
      !> Euclidean adjoint (transpose) of the backward transform (r -> k)
      procedure(radial_trafo_fbt_i), deferred :: fbt_k2r_adj
      !> Batched forward transform over the columns of a 2D field (r -> k)
      procedure :: fbt_r2k_all => radial_trafo_fbt_r2k_all_default
      !> Batched backward transform over the columns of a 2D field (k -> r)
      procedure :: fbt_k2r_all => radial_trafo_fbt_k2r_all_default
      !> Batched adjoint of the forward transform (k -> r)
      procedure :: fbt_r2k_adj_all => radial_trafo_fbt_r2k_adj_all_default
      !> Batched adjoint of the backward transform (r -> k)
      procedure :: fbt_k2r_adj_all => radial_trafo_fbt_k2r_adj_all_default
      !> Release scratch and detach from the grid (idempotent)
      procedure(radial_trafo_destroy_i), deferred :: destroy
   end type moist_math_grid_1d_trafo_type

   !> Deferred grid and trafo operations shared by every concrete 1D grid
   abstract interface
      !> 3D-radial integration weight (volume element) of grid point i
      !>
      !> @param[in] self  Grid instance
      !> @param[in] i     Grid point index
      pure function radial_measure_i(self, i) result(w)
         import :: moist_math_grid_1d_type, wp
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_1d_type), intent(in) :: self
         !> Grid point index
         integer, intent(in) :: i
         !> Weight of point i
         real(wp) :: w
      end function radial_measure_i

      !> Allocate the transform engine matching this grid and bind it to `self`
      !>
      !> `self` must have the target attribute and outlive the trafo
      !>
      !> @param[in]  self   Grid instance, target
      !> @param[out] trafo  Allocated, bound transform engine
      !> @param[out] error  Error handling
      subroutine radial_new_trafo_i(self, trafo, error)
         import :: moist_math_grid_1d_type, moist_math_grid_1d_trafo_type, error_type
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_1d_type), intent(in), target :: self
         !> Allocated, bound transform engine
         class(moist_math_grid_1d_trafo_type), allocatable, intent(out) :: trafo
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine radial_new_trafo_i

      !> Release all grid storage, idempotent
      !>
      !> Not pure: concrete grids may tear down transform plans or other
      !> external handles
      !>
      !> @param[in,out] self  Grid instance
      subroutine radial_destroy_i(self)
         import :: moist_math_grid_1d_type
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_1d_type), intent(inout) :: self
      end subroutine radial_destroy_i

      !> Apply a spherical Fourier-Bessel transform, mapping `f_in` to `f_out`
      !>
      !> Forward, backward, or either Euclidean adjoint, depending on the
      !> bound procedure; `self` is mutated (it owns the scratch) but is
      !> thread-local, so this is safe to call concurrently from distinct
      !> threads on distinct trafos; a backend failure is reported through
      !> `error`
      !>
      !> @param[in,out] self   Trafo instance (thread-local)
      !> @param[in]     f_in   Input field
      !> @param[out]    f_out  Output field
      !> @param[out]    error  Error handling
      subroutine radial_trafo_fbt_i(self, f_in, f_out, error)
         import :: moist_math_grid_1d_trafo_type, wp, error_type
         implicit none(type, external)
         !> Trafo instance (thread-local)
         class(moist_math_grid_1d_trafo_type), intent(inout) :: self
         !> Input field
         real(wp), intent(in) :: f_in(:)
         !> Output field
         real(wp), intent(out) :: f_out(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine radial_trafo_fbt_i

      !> Release trafo scratch and detach from the grid; idempotent
      !>
      !> @param[in,out] self  Trafo instance
      subroutine radial_trafo_destroy_i(self)
         import :: moist_math_grid_1d_trafo_type
         implicit none(type, external)
         !> Trafo instance
         class(moist_math_grid_1d_trafo_type), intent(inout) :: self
      end subroutine radial_trafo_destroy_i
   end interface

contains

   !> Default field quadrature: `result = sum_i measure(i) * f(i)`
   !>
   !> Consumes discrete field data already tabulated on the nodes (e.g.
   !> correlation functions from a converged RISM solve)
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-node field values (length npts)
   !> @param[out] result  Quadrature result
   pure subroutine radial_integrate_field_default(self, f, result)
      !> Grid instance
      class(moist_math_grid_1d_type), intent(in) :: self
      !> Per-node field values (length npts)
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + self%measure(i)*f(i)
      end do
   end subroutine radial_integrate_field_default

   !> Default analytic quadrature: `result = sum_i measure(i) * f(r_i)`
   !>
   !> Samples the integrand at the radial nodes the grid owns
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Radial integrand f(r)
   !> @param[out] result  Quadrature result
   subroutine radial_integrate_default(self, f, result)
      !> Grid instance
      class(moist_math_grid_1d_type), intent(in) :: self
      !> Radial integrand f(r)
      procedure(integrand_radial) :: f
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%npts
         result = result + self%measure(i)*f(self%r(i))
      end do
   end subroutine radial_integrate_default

   !> Default batched forward transform: scalar forward transform per column
   !>
   !> Applies the scalar forward transform to each column of `f_in`
   !> independently; concrete trafos override this with a genuinely batched
   !> implementation (e.g. a batched DST or a GEMM); stops at the first
   !> column whose scalar transform reports an error
   !>
   !> @param[in,out] self   Trafo instance (thread-local)
   !> @param[in]     f_in   Input fields, one per column
   !> @param[out]    f_out  Output fields, one per column
   !> @param[out]    error  Set if a column's scalar transform failed
   subroutine radial_trafo_fbt_r2k_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_1d_trafo_type), intent(inout) :: self
      !> Input fields, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> Output fields, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      do j = 1, size(f_in, 2)
         call self%fbt_r2k(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_r2k_all_default

   !> Default batched backward transform: scalar backward transform per column
   !>
   !> Stops at the first column whose scalar transform reports an error
   !>
   !> @param[in,out] self   Trafo instance (thread-local)
   !> @param[in]     f_in   Input fields, one per column
   !> @param[out]    f_out  Output fields, one per column
   !> @param[out]    error  Set if a column's scalar transform failed
   subroutine radial_trafo_fbt_k2r_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_1d_trafo_type), intent(inout) :: self
      !> Input fields, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> Output fields, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      do j = 1, size(f_in, 2)
         call self%fbt_k2r(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_k2r_all_default

   !> Default batched forward adjoint: scalar forward adjoint per column
   !>
   !> Stops at the first column whose scalar transform reports an error
   !>
   !> @param[in,out] self   Trafo instance (thread-local)
   !> @param[in]     f_in   Input fields, one per column
   !> @param[out]    f_out  Output fields, one per column
   !> @param[out]    error  Set if a column's scalar transform failed
   subroutine radial_trafo_fbt_r2k_adj_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_1d_trafo_type), intent(inout) :: self
      !> Input fields, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> Output fields, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      do j = 1, size(f_in, 2)
         call self%fbt_r2k_adj(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_r2k_adj_all_default

   !> Default batched backward adjoint: scalar backward adjoint per column
   !>
   !> Stops at the first column whose scalar transform reports an error
   !>
   !> @param[in,out] self   Trafo instance (thread-local)
   !> @param[in]     f_in   Input fields, one per column
   !> @param[out]    f_out  Output fields, one per column
   !> @param[out]    error  Set if a column's scalar transform failed
   subroutine radial_trafo_fbt_k2r_adj_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_1d_trafo_type), intent(inout) :: self
      !> Input fields, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> Output fields, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      do j = 1, size(f_in, 2)
         call self%fbt_k2r_adj(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_k2r_adj_all_default

end module moist_math_grid_1d_base
