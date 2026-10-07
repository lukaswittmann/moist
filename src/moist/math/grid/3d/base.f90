!> Abstract volumetric domain and separate 3D Fourier transform engine
module moist_math_grid_3d_base
   use mctc_env, only: wp, error_type, fatal_error
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_io, only: structure_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
!$ use omp_lib, only: omp_get_max_threads, omp_in_parallel
   implicit none(type, external)
   private

   public :: moist_math_grid_3d_type
   public :: moist_math_grid_3d_trafo_type
   public :: integrand_3d
   public :: check_trafo_blocks

   !> Scalar-valued 3D integrand used by `integrate`
   abstract interface
      !> Scalar field evaluated at a spatial point
      !>
      !> @param[in] r Point in bohr
      pure function integrand_3d(r) result(val)
         import :: wp
         implicit none(type, external)
         !> Point in bohr
         real(wp), intent(in) :: r(3)
         !> Function value at r
         real(wp) :: val
      end function integrand_3d
   end interface

   !> Abstract 3D domain with quadrature and a transform factory
   !>
   !> - Real-space points `point(i)`, i = 1..ngrid, and reciprocal-space
   !>   points `kpoint(j)`, j = 1..npts_k, addressed per index
   !> - Field-construction code (potentials, susceptibility tabulation) stays
   !>   grid-agnostic: no dependence on the concrete grid's `x/y/z` or
   !>   `kx/ky/kz` layout
   !> - `kpoint(j)` ordering matches the reciprocal-space field layout
   !>   produced by the grid's transform engine
   type, abstract :: moist_math_grid_3d_type
      !> Number of real-space points
      integer :: ngrid = 0
      !> Number of solute atoms
      integer :: natom = 0
      !> Real-space coordinates, bohr (3, ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Volume integration weights, bohr**3
      real(wp), allocatable :: w(:)
      !> Owning atom, one based; zero for unowned points
      integer, allocatable :: owner(:)
      !> Optional Gaussian exponents, inverse bohr
      real(wp), allocatable :: xi0(:)
      !> Number of reciprocal-space (k-space) points
      integer :: npts_k = 0
      !> Thread count for grid builds and transforms; 0 takes omp_get_max_threads
      integer :: nthreads = 0
   contains
      !> Rebuild geometry for the molecule
      procedure(update_grid), deferred :: update
      !> Contract volume adjoints into the nuclear gradient
      procedure :: get_volume_gradient => get_grid_volume_gradient_default
      !> Contract fixed volume adjoints into a nuclear Hessian-vector product
      procedure :: get_volume_hessian_vector => get_grid_volume_hessian_vector_default
      !> Contract fixed volume adjoints into the full nuclear Hessian
      procedure :: get_volume_hessian => get_grid_volume_hessian
      !> Validate state, adjoint channels and the nuclear tensors of one contraction
      procedure :: check_volume_adjoint => check_grid_volume_adjoint
      !> Whether Gaussian widths follow molecular geometry
      procedure :: has_geometry_dependent_xi0 => grid_xi0_dependent_default
      !> Diagnostic grid name
      procedure :: kind_name => grid_kind_name_default
      !> Worker count for grid builds and transforms; 1 inside an OpenMP region
      procedure :: team_size => grid_team_size
      !> Real-space coordinate (bohr) of grid point i (1..ngrid)
      procedure :: point => grid_point
      !> Reciprocal-space coordinate (1/bohr) of k-point j (1..npts_k)
      procedure(grid_kpoint_i), deferred :: kpoint
      !> Weighted volume integral of a scalar integrand sampled at the points
      procedure :: integrate => grid_integrate
      !> Integration weight (volume element) of grid point i
      procedure :: measure => grid_measure
      !> Allocate the transform engine matching this grid, bound to it
      procedure(grid_new_trafo_i), deferred :: new_trafo
      !> Release all storage (idempotent)
      procedure(grid_destroy_i), deferred :: destroy
      !> Weighted volume integral of a field already tabulated on the points
      procedure :: integrate_field => grid_integrate_field_default
      !> Common storage cleanup used by concrete destroy methods
      procedure :: clear_geometry => grid_clear_geometry
   end type moist_math_grid_3d_type

   !> Abstract 3D Fourier transform engine
   !>
   !> - Bound to one grid, held by the concrete type as a typed pointer; owns
   !>   whatever immutable transform plans/scratch the backend needs
   !> - Built once per problem (`new_trafo` + `prepare`), reused across solver
   !>   iterations; must not outlive the grid it is bound to
   !> - Batched `fft_r2k`/`fft_k2r` each issued from a single thread (the
   !>   backend threads the work internally), one instance shared rather than
   !>   per-thread
   type, abstract :: moist_math_grid_3d_trafo_type
   contains
      !> Forward Fourier transform (real space -> reciprocal space), all sites
      procedure(grid_trafo_fft_r2k_i), deferred :: fft_r2k
      !> Backward Fourier transform (reciprocal space -> real space), all sites
      procedure(grid_trafo_fft_k2r_i), deferred :: fft_k2r
      !> Release scratch and detach from the grid (idempotent)
      procedure(grid_trafo_destroy_i), deferred :: destroy
      !> Build the backend plans for `ntrans` simultaneous transforms;
      !> default is a no-op (stateless backends such as the Cartesian FFT
      !> need nothing)
      procedure :: prepare => grid_trafo_prepare_default
   end type moist_math_grid_3d_trafo_type

   !> Deferred grid and transform contracts
   abstract interface

      !> Rebuild grid geometry for the molecular structure
      !>
      !> @param[in,out] self Domain instance
      !> @param[in] mol Molecular structure
      !> @param[out] error Error handling
      subroutine update_grid(self, mol, error)
         import :: moist_math_grid_3d_type, structure_type, error_type
         implicit none(type, external)
         !> Domain instance
         class(moist_math_grid_3d_type), intent(inout) :: self
         !> Molecular structure
         type(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_grid

      !> Reciprocal-space coordinate (1/bohr) of k-point j (1..npts_k)
      !>
      !> @param[in] self Grid instance
      !> @param[in] j k-point index
      pure function grid_kpoint_i(self, j) result(k)
         import :: moist_math_grid_3d_type, wp
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_3d_type), intent(in) :: self
         !> k-point index
         integer, intent(in) :: j
         !> Coordinate of k-point j
         real(wp) :: k(3)
      end function grid_kpoint_i

      !> Allocate the transform engine matching this grid and bind it to `self`
      !>
      !> `self` must have the target attribute and outlive the trafo
      !>
      !> @param[in] self Grid instance
      !> @param[out] trafo Allocated, bound transform engine
      !> @param[out] error Error handling
      subroutine grid_new_trafo_i(self, trafo, error)
         import :: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type, error_type
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_3d_type), intent(in), target :: self
         !> Allocated, bound transform engine
         class(moist_math_grid_3d_trafo_type), allocatable, intent(out) :: trafo
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine grid_new_trafo_i

      !> Release all grid storage; idempotent
      !>
      !> Not pure: concrete grids may tear down transform plans or other
      !> external handles
      !>
      !> @param[in,out] self Grid instance
      subroutine grid_destroy_i(self)
         import :: moist_math_grid_3d_type
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_3d_type), intent(inout) :: self
      end subroutine grid_destroy_i

      !> Forward FT of all `nv` sites at once, real-space to reciprocal-space
      !>
      !> Real-space block `f_r(ngrid, nv)` -> reciprocal-space block
      !> `f_k(npts_k, nv)`; `self` may be mutated (it owns scratch); issue
      !> from a single thread; on error `f_k` is unspecified
      !>
      !> @param[in,out] self Trafo instance
      !> @param[in,out] f_r Real-space field block, shape (ngrid, nv)
      !> @param[out] f_k Reciprocal-space field block, shape (npts_k, nv)
      !> @param[out] error Error handling
      subroutine grid_trafo_fft_r2k_i(self, f_r, f_k, error)
         import :: moist_math_grid_3d_trafo_type, wp, error_type
         implicit none(type, external)
         !> Trafo instance
         class(moist_math_grid_3d_trafo_type), intent(inout) :: self
         !> Real-space field block, shape (ngrid, nv)
         real(wp), intent(inout), contiguous, target :: f_r(:, :)
         !> Reciprocal-space field block, shape (npts_k, nv)
         complex(wp), intent(out), contiguous, target :: f_k(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine grid_trafo_fft_r2k_i

      !> Backward FT of all `nv` sites at once: reciprocal-space block
      !>
      !> `f_k(npts_k, nv)` -> real-space block `f_r(ngrid, nv)`.  On error
      !> `f_r` is unspecified
      !>
      !> @param[in,out] self Trafo instance
      !> @param[in,out] f_k reciprocal-space field, shape (npts_k, nv); used as scratch, destroyed on return
      !> @param[out] f_r Real-space field block, shape (ngrid, nv)
      !> @param[out] error Error handling
      subroutine grid_trafo_fft_k2r_i(self, f_k, f_r, error)
         import :: moist_math_grid_3d_trafo_type, wp, error_type
         implicit none(type, external)
         !> Trafo instance
         class(moist_math_grid_3d_trafo_type), intent(inout) :: self
         !> Reciprocal-space field block, shape (npts_k, nv); used as transform scratch, treat as destroyed on return
         complex(wp), intent(inout), contiguous, target :: f_k(:, :)
         !> Real-space field block, shape (ngrid, nv)
         real(wp), intent(out), contiguous, target :: f_r(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine grid_trafo_fft_k2r_i

      !> Release trafo scratch and detach from the grid; idempotent
      !>
      !> @param[in,out] self Trafo instance
      subroutine grid_trafo_destroy_i(self)
         import :: moist_math_grid_3d_trafo_type
         implicit none(type, external)
         !> Trafo instance
         class(moist_math_grid_3d_trafo_type), intent(inout) :: self
      end subroutine grid_trafo_destroy_i
   end interface

contains

   !> Validate realized state, adjoint channels and the nuclear tensors of one contraction
   !>
   !> - Every call checks the update state and the adjoint channels
   !> - `direction` or `hessian` marks a curvature contraction, which rejects
   !>   geometry-dependent Gaussian-width adjoints
   !>
   !> @param[in] self Realized grid
   !> @param[in] acc Volume-observable adjoints
   !> @param[out] error Invalid state, channels or shapes
   !> @param[in] vector Optional nuclear gradient or Hessian-vector accumulator (3, natom)
   !> @param[in] direction Optional nuclear displacement direction (3, natom)
   !> @param[in] hessian Optional nuclear Hessian accumulator (3, natom, 3, natom)
   subroutine check_grid_volume_adjoint(self, acc, error, vector, direction, hessian)
      !> Realized grid
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Volume adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Invalid input
      type(error_type), allocatable, intent(out) :: error
      !> Nuclear gradient or Hessian-vector accumulator
      real(wp), intent(in), optional :: vector(:, :)
      !> Nuclear displacement direction
      real(wp), intent(in), optional :: direction(:, :)
      !> Nuclear Hessian accumulator
      real(wp), intent(in), optional :: hessian(:, :, :, :)

      ! Updates zero natom on entry and commit it only on success
      if (self%natom < 1 .or. self%ngrid < 1) then
         call fatal_error(error, "volume adjoint: update the grid successfully first")
         return
      end if
      if (.not. acc%is_initialized()) then
         call fatal_error(error, "volume adjoint: accumulator is not initialized")
         return
      end if
      if (acc%size() /= self%ngrid) then
         call fatal_error(error, "volume adjoint: grid size mismatch")
         return
      end if
      if (.not. all(ieee_is_finite(acc%w_xyz)) .or. .not. all(ieee_is_finite(acc%w_w)) &
         & .or. .not. all(ieee_is_finite(acc%w_xi))) then
         call fatal_error(error, "volume adjoint: adjoints must be finite")
         return
      end if
      if (.not. allocated(self%xi0)) then
         if (any(acc%w_xi /= 0.0_wp)) then
            call fatal_error(error, "volume adjoint: point grids have no Gaussian-width channel")
            return
         end if
      end if
      if (present(vector)) then
         if (any(shape(vector) /= [3, self%natom])) then
            call fatal_error(error, "volume adjoint: nuclear accumulator shape mismatch")
            return
         end if
      end if
      if (present(direction)) then
         if (any(shape(direction) /= [3, self%natom])) then
            call fatal_error(error, "volume Hessian: nuclear-direction shape mismatch")
            return
         end if
         if (.not. all(ieee_is_finite(direction))) then
            call fatal_error(error, "volume Hessian: nuclear direction must be finite")
            return
         end if
      end if
      if (present(hessian)) then
         if (any(shape(hessian) /= [3, self%natom, 3, self%natom])) then
            call fatal_error(error, "volume Hessian: nuclear-Hessian shape mismatch")
            return
         end if
      end if
      if (present(direction) .or. present(hessian)) then
         if (self%has_geometry_dependent_xi0()) then
            if (any(acc%w_xi /= 0.0_wp)) then
               call fatal_error(error, "volume Hessian: geometry-dependent Gaussian-width curvature is unsupported")
               return
            end if
         end if
      end if
   end subroutine check_grid_volume_adjoint

   !> Assemble a full nuclear Hessian from exact directional contractions
   !>
   !> @param[in] self Successfully updated grid
   !> @param[in] acc Fixed volume adjoints
   !> @param[in,out] hessian Nuclear Hessian accumulator (3, natom, 3, natom)
   !> @param[out] error Invalid input or allocation failure
   subroutine get_grid_volume_hessian(self, acc, hessian, error)
      !> Updated grid
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Fixed volume adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Nuclear Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Invalid input or allocation failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: local(:, :, :, :), direction(:, :)
      integer :: a, c, stat

      call self%check_volume_adjoint(acc, error, hessian=hessian)
      if (allocated(error)) return
      allocate (local(3, self%natom, 3, self%natom), stat=stat)
      if (stat == 0) allocate (direction(3, self%natom), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "volume Hessian: cannot allocate contraction scratch")
         return
      end if
      local = 0.0_wp
      do a = 1, self%natom
         do c = 1, 3
            direction = 0.0_wp
            direction(c, a) = 1.0_wp
            call self%get_volume_hessian_vector(acc, direction, local(:, :, c, a), error)
            if (allocated(error)) return
         end do
      end do
      hessian = hessian + local
   end subroutine get_grid_volume_hessian

   !> Default error for grids without a second nuclear response
   !>
   !> @param[in] self Grid instance
   !> @param[in] acc Fixed volume adjoints
   !> @param[in] direction Nuclear displacement direction
   !> @param[in,out] hessian_vector Nuclear Hessian-vector accumulator
   !> @param[out] error Unsupported grid
   subroutine get_grid_volume_hessian_vector_default(self, acc, direction, hessian_vector, error)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Fixed volume adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Nuclear direction
      real(wp), intent(in) :: direction(:, :)
      !> Nuclear Hessian-vector accumulator
      real(wp), intent(inout) :: hessian_vector(:, :)
      !> Unsupported grid
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, self%kind_name()//": volume Hessian is not implemented")
   end subroutine get_grid_volume_hessian_vector_default

   !> Weighted volume integral sampled on the domain coordinates
   !>
   !> @param[in] self Grid instance
   !> @param[in] f Integrand
   !> @param[out] result Volume integral
   subroutine grid_integrate(self, f, result)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Integrand
      procedure(integrand_3d) :: f
      !> Volume integral
      real(wp), intent(out) :: result
      !> Point index
      integer :: i

      result = 0.0_wp
      do i = 1, self%ngrid
         result = result + self%w(i)*f(self%xyz(:, i))
      end do
   end subroutine grid_integrate

   !> Real-space coordinate from the authoritative domain geometry
   !>
   !> @param[in] self Grid instance
   !> @param[in] i Point index
   pure function grid_point(self, i) result(r)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Point index
      integer, intent(in) :: i
      !> Point coordinates, bohr
      real(wp) :: r(3)
      r = self%xyz(:, i)
   end function grid_point

   !> Volume weight from the authoritative domain discretization
   !>
   !> @param[in] self Grid instance
   !> @param[in] i Point index
   pure function grid_measure(self, i) result(w)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Point index
      integer, intent(in) :: i
      !> Integration weight, bohr**3
      real(wp) :: w
      w = self%w(i)
   end function grid_measure

   !> Release common discretization storage
   !>
   !> @param[in,out] self Grid instance
   subroutine grid_clear_geometry(self)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(inout) :: self
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%w)) deallocate (self%w)
      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%xi0)) deallocate (self%xi0)
      self%ngrid = 0
      self%natom = 0
      self%npts_k = 0
   end subroutine grid_clear_geometry

   !> Default transform preparation: a no-op
   !>
   !> Cartesian FFT execution needs no plans and uses this default
   !> Molecular NUFFT engines override it to allocate batch-sized plans
   !>
   !> @param[in,out] self    Trafo instance
   !> @param[in]    ntrans  Number of simultaneous transforms to plan for
   !> @param[out]   error   Propagated error (never raised by the default)
   subroutine grid_trafo_prepare_default(self, ntrans, error)
      !> Trafo instance
      class(moist_math_grid_3d_trafo_type), intent(inout) :: self
      !> Number of simultaneous transforms (unused by the stateless default)
      integer, intent(in) :: ntrans

      !> Error handling (left unallocated: the default never fails)
      type(error_type), allocatable, intent(out) :: error

      ! Stateless default: nothing to build; self/ntrans deliberately unused,
      ! error returns unallocated
   end subroutine grid_trafo_prepare_default

   !> Validate the field blocks of one batched transform call
   !>
   !> - Grid built: `ngrid >= 1` and `npts_k >= 1`; a destroyed grid keeps its
   !>   backend extents, so empty blocks must not reach native code
   !> - Real-space block `(ngrid, nv)`, reciprocal-space block `(npts_k, nv)`
   !> - Both blocks share the batch width `nv >= 1`
   !> - Run before any backend call: native code trusts these extents
   !>
   !> @param[in]  ngrid    Real-space point count of the bound grid
   !> @param[in]  npts_k   Reciprocal-space point count of the bound grid
   !> @param[in]  shape_r  Shape of the caller's real-space block
   !> @param[in]  shape_k  Shape of the caller's reciprocal-space block
   !> @param[out] error    Set when the blocks do not fit the grid
   subroutine check_trafo_blocks(ngrid, npts_k, shape_r, shape_k, error)
      !> Real-space point count
      integer, intent(in) :: ngrid
      !> Reciprocal-space point count
      integer, intent(in) :: npts_k
      !> Shape of the real-space block
      integer, intent(in) :: shape_r(2)
      !> Shape of the reciprocal-space block
      integer, intent(in) :: shape_k(2)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (ngrid < 1 .or. npts_k < 1) then
         call fatal_error(error, "grid trafo: grid has no points; update the grid first")
      else if (shape_r(1) /= ngrid .or. shape_k(1) /= npts_k) then
         call fatal_error(error, "grid trafo: field block size does not match the grid")
      else if (shape_r(2) /= shape_k(2)) then
         call fatal_error(error, "grid trafo: real- and reciprocal-space blocks differ in batch width")
      else if (shape_r(2) < 1) then
         call fatal_error(error, "grid trafo: batch needs at least one column")
      end if
   end subroutine check_trafo_blocks

   !> Default field quadrature: `result = sum_i measure(i) * f(i)`
   !>
   !> - Unlike `integrate`, which samples an analytic function at the grid
   !>   points, consumes discrete field data already tabulated on the grid
   !>   (e.g. correlation functions from a converged RISM solve)
   !> - Generic version walks every point and queries its weight; concrete
   !>   grids with a uniform (or otherwise cheap) measure should override it
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-point field values (length ngrid)
   !> @param[out] result  Quadrature result
   pure subroutine grid_integrate_field_default(self, f, result)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Per-point field values (length ngrid)
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result

      integer :: i

      result = 0.0_wp
      do i = 1, self%ngrid
         result = result + self%measure(i)*f(i)
      end do
   end subroutine grid_integrate_field_default

   !> Report an unavailable volume geometry gradient
   !>
   !> @param[in] self Domain instance
   !> @param[in] acc Accumulated adjoints
   !> @param[in,out] gradient Nuclear-gradient accumulator
   !> @param[out] error Error handling
   subroutine get_grid_volume_gradient_default(self, acc, gradient, error)
      !> Domain instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Accumulated adjoints
      type(volume_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "The "//self%kind_name()// &
         & " domain does not provide a reverse-mode geometry gradient")

   end subroutine get_grid_volume_gradient_default

   !> Default Gaussian widths are geometry-independent
   !>
   !> @param[in] self Domain instance
   function grid_xi0_dependent_default(self) result(xi0_dependent)
      !> Domain instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Whether `xi0` depends on the nuclear positions
      logical :: xi0_dependent

      xi0_dependent = .false.

   end function grid_xi0_dependent_default

   !> Diagnostic name for a volume grid
   !>
   !> @param[in] self Domain instance
   function grid_kind_name_default(self) result(name)
      !> Domain instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Kind name
      character(len=:), allocatable :: name

      name = "generic"

   end function grid_kind_name_default

   !> Worker count for grid builds and transforms
   !>
   !> - 1 inside an OpenMP region: the C++ backends do not nest
   !> - otherwise `nthreads`, or `omp_get_max_threads` when unset
   !>
   !> @param[in] self Domain instance
   function grid_team_size(self) result(nt)
      !> Domain instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Worker count (>= 1)
      integer :: nt

      nt = 1
!$    if (.not. omp_in_parallel()) then
!$       nt = max(1, omp_get_max_threads())
!$       if (self%nthreads > 0) nt = self%nthreads
!$    end if

   end function grid_team_size

end module moist_math_grid_3d_base
