!> Abstract volumetric domain and separate 3D Fourier transform engine
module moist_math_grid_3d_base
   use mctc_env, only: wp, error_type, fatal_error
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_io, only: structure_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_integration, only: integrate_weighted_values
   use moist_math_grid_3d_kernel_integration, only: compensated_block_sum
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
   !> - Transform phases measured from `get_kref()`:
   !>   `f_k = sum_j w_j f(r_j) exp(-i k.(r_j - kref))`; analytic reciprocal
   !>   tables use the same reference
   type, abstract :: moist_math_grid_3d_type
      !> Number of real-space points
      integer :: ngrid = 0
      !> Number of solute atoms
      integer :: natom = 0
      !> Real-space coordinates, bohr (3, ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Final volume integration weights, bohr**3 (ngrid)
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
      !> Check the configured settings before any geometry is built
      procedure(validate_grid_i), deferred :: validate
      !> Common update preamble: mark unusable, validate, check the structure
      procedure, non_overridable :: check_update_input => grid_check_update_input
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
      !> Phase reference of the transforms (bohr)
      procedure(grid_get_kref_i), deferred :: get_kref
      !> Weighted volume integral of a scalar integrand sampled at the points
      procedure :: integrate => grid_integrate
      !> Allocate the transform engine matching this grid, bound to it
      procedure(grid_new_trafo_i), deferred :: new_trafo
      !> Release all storage (idempotent)
      procedure(grid_destroy_i), deferred :: destroy
      !> Weighted volume integral of a field already tabulated on the points
      procedure :: integrate_field => grid_integrate_field_default
      !> Compensated contribution of a field slice; merge its two components in block order
      procedure :: integrate_field_block => grid_integrate_field_block
      !> Common storage cleanup used by concrete destroy methods
      procedure :: clear_geometry => grid_clear_geometry
   end type moist_math_grid_3d_type

   !> Abstract 3D Fourier transform engine
   !>
   !> - Bound to one grid, held by the concrete type as a typed pointer; owns
   !>   the transform plans and mutable scratch the backend needs
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
      !> Transpose of fft_r2k under Re(sum(conjg(w_k)*df_k)) and sum(w_r*df_r)
      !>
      !> - Reciprocal input may be destroyed; includes quadrature and packed-spectrum weights
      !> - Optional; the default reports an unsupported operation
      procedure :: fft_r2k_adj => grid_trafo_fft_r2k_adj_default
      !> Transpose of fft_k2r with the same Euclidean real/complex pairing
      !>
      !> - Real input may be destroyed; includes inverse-transform normalization
      !> - Optional; the default reports an unsupported operation
      procedure :: fft_k2r_adj => grid_trafo_fft_k2r_adj_default
      !> Geometry pullback of fft_r2k at fixed field values and reciprocal lattice
      !>
      !> - Optional; the default reports an unsupported operation
      procedure :: accumulate_r2k_geometry => grid_trafo_r2k_geometry_default
      !> Geometry pullback of fft_k2r at fixed field values and reciprocal lattice
      !>
      !> - Optional; the default reports an unsupported operation
      procedure :: accumulate_k2r_geometry => grid_trafo_k2r_geometry_default
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

      !> Check the configured settings before any geometry is built
      !>
      !> @param[in] self Domain instance
      !> @param[out] error Invalid setting
      subroutine validate_grid_i(self, error)
         import :: moist_math_grid_3d_type, error_type
         implicit none(type, external)
         !> Domain instance
         class(moist_math_grid_3d_type), intent(in) :: self
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine validate_grid_i

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

      !> Phase reference of the transforms (bohr)
      !>
      !> @param[in] self Grid instance
      pure function grid_get_kref_i(self) result(kref)
         import :: moist_math_grid_3d_type, wp
         implicit none(type, external)
         !> Grid instance
         class(moist_math_grid_3d_type), intent(in) :: self
         !> Real-space phase origin (bohr)
         real(wp) :: kref(3)
      end function grid_get_kref_i

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
         & .or. .not. all(ieee_is_finite(acc%w_xi)) .or. .not. all(ieee_is_finite(acc%w_kref))) then
         call fatal_error(error, "volume adjoint: adjoints must be finite")
         return
      end if
      if (.not. allocated(self%xi0)) then
         if (any(acc%w_xi /= 0.0_wp)) then
            call fatal_error(error, "volume adjoint: point grids have no Gaussian-width channel")
            return
         end if
      end if
      if (self%npts_k < 1) then
         if (any(acc%w_kref /= 0.0_wp)) then
            call fatal_error(error, "volume adjoint: grids without reciprocal points have no phase-reference channel")
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
      real(wp) :: correction, term, subtotal

      result = 0.0_wp
      correction = 0.0_wp
      ! Preserve small contributions after the accumulated integral grows
      do i = 1, self%ngrid
         term = self%w(i)*f(self%xyz(:, i)) - correction
         subtotal = result + term
         correction = (subtotal - result) - term
         result = subtotal
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

   !> Common update preamble shared by every concrete grid
   !>
   !> - Zero `natom` first: any failure from here on leaves the grid unusable
   !> - Validate the settings, then require at least one atom and finite coordinates
   !> - Concrete updates add their own geometry checks afterwards
   !>
   !> @param[in,out] self Domain instance
   !> @param[in] mol Requested solute structure
   !> @param[out] error Invalid settings or structure
   subroutine grid_check_update_input(self, mol, error)
      !> Domain instance
      class(moist_math_grid_3d_type), intent(inout) :: self
      !> Requested solute structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      self%natom = 0
      call self%validate(error)
      if (allocated(error)) return
      if (mol%nat < 1) then
         call fatal_error(error, self%kind_name()//" domain: at least one solute atom is required")
         return
      end if
      if (.not. all(ieee_is_finite(mol%xyz))) then
         call fatal_error(error, self%kind_name()//" domain: solute coordinates must be finite")
      end if
   end subroutine grid_check_update_input

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

   !> Report an unsupported forward-transform transpose
   !>
   !> @param[in,out] self Transform engine
   !> @param[in,out] f_k Complex output seed (npts_k, nv)
   !> @param[out] f_r Real input sensitivity (ngrid, nv); unspecified on error
   !> @param[out] error Unsupported operation
   subroutine grid_trafo_fft_r2k_adj_default(self, f_k, f_r, error)
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: self
      !> Complex output seed (npts_k, nv)
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      !> Real input sensitivity (ngrid, nv); unspecified on error
      real(wp), intent(out), contiguous, target :: f_r(:, :)
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "3D transform does not support fft_r2k_adj")
   end subroutine grid_trafo_fft_r2k_adj_default

   !> Report an unsupported backward-transform transpose
   !>
   !> @param[in,out] self Transform engine
   !> @param[in,out] f_r Real output seed (ngrid, nv)
   !> @param[out] f_k Complex input sensitivity (npts_k, nv); unspecified on error
   !> @param[out] error Unsupported operation
   subroutine grid_trafo_fft_k2r_adj_default(self, f_r, f_k, error)
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: self
      !> Real output seed (ngrid, nv)
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      !> Complex input sensitivity (npts_k, nv); unspecified on error
      complex(wp), intent(out), contiguous, target :: f_k(:, :)
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "3D transform does not support fft_k2r_adj")
   end subroutine grid_trafo_fft_k2r_adj_default

   !> Report unsupported forward-transform geometry derivatives
   !>
   !> @param[in,out] self Transform engine
   !> @param[in] f_r Fixed primal real field (ngrid, nv)
   !> @param[in] w_k Complex output seed under Re(sum(conjg(w_k)*dF)) (npts_k, nv)
   !> @param[in,out] acc Geometry sensitivities, unchanged on error
   !> @param[out] error Unsupported operation
   subroutine grid_trafo_r2k_geometry_default(self, f_r, w_k, acc, error)
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: self
      !> Fixed primal real field (ngrid, nv)
      real(wp), intent(in) :: f_r(:, :)
      !> Complex output seed (npts_k, nv)
      complex(wp), intent(in) :: w_k(:, :)
      !> Geometry sensitivities, unchanged on error
      type(volume_adjoint_type), intent(inout) :: acc
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "3D transform does not support accumulate_r2k_geometry")
   end subroutine grid_trafo_r2k_geometry_default

   !> Report unsupported backward-transform geometry derivatives
   !>
   !> @param[in,out] self Transform engine
   !> @param[in] f_k Fixed primal reciprocal field (npts_k, nv)
   !> @param[in] w_r Real output seed under sum(w_r*df) (ngrid, nv)
   !> @param[in,out] acc Geometry sensitivities, unchanged on error
   !> @param[out] error Unsupported operation
   subroutine grid_trafo_k2r_geometry_default(self, f_k, w_r, acc, error)
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: self
      !> Fixed primal reciprocal field (npts_k, nv)
      complex(wp), intent(in) :: f_k(:, :)
      !> Real output seed (ngrid, nv)
      real(wp), intent(in) :: w_r(:, :)
      !> Geometry sensitivities, unchanged on error
      type(volume_adjoint_type), intent(inout) :: acc
      !> Unsupported operation
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "3D transform does not support accumulate_k2r_geometry")
   end subroutine grid_trafo_k2r_geometry_default

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

   !> Default field quadrature: `result = sum_i w(i) * f(i)`
   !>
   !> - Unlike `integrate`, which samples an analytic function at the grid
   !>   points, consumes discrete field data already tabulated on the grid
   !>   (e.g. correlation functions from a converged RISM solve)
   !> - Final weights from w(:), with fixed compensated blocks for large fields
   !>
   !> @param[in]  self    Grid instance
   !> @param[in]  f       Per-point field values (length ngrid); any other length is an error
   !> @param[out] result  Quadrature result
   !> @param[out] error   Error handling
   subroutine grid_integrate_field_default(self, f, result, error)
      !> Grid instance
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Per-point field values (length ngrid)
      real(wp), intent(in) :: f(:)
      !> Quadrature result
      real(wp), intent(out) :: result
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      result = 0.0_wp
      if (size(f) /= self%ngrid) then
         call fatal_error(error, "3d grid: integrate_field needs one value per grid point")
         return
      end if
      if (size(f) == 0) return
      result = integrate_weighted_values(f, self%team_size(), self%w)
   end subroutine grid_integrate_field_default

   !> Weighted contribution of f at grid points first:first+size(f)-1
   !>
   !> @param[in] self Grid providing final volume weights
   !> @param[in] f Bounded field slice
   !> @param[in] first One-based absolute grid index of f(1)
   !> @param[out] partial Neumaier total and correction, both in integral units
   !> @param[out] error Slice outside the grid
   subroutine grid_integrate_field_block(self, f, first, partial, error)
      !> Grid providing final volume weights
      class(moist_math_grid_3d_type), intent(in) :: self
      !> Bounded field slice
      real(wp), intent(in) :: f(:)
      !> One-based absolute grid index of f(1)
      integer, intent(in) :: first
      !> Neumaier total and correction, both in integral units
      real(wp), intent(out) :: partial(2)
      !> Slice outside the grid
      type(error_type), allocatable, intent(out) :: error

      partial = 0.0_wp
      if (first < 1 .or. first > self%ngrid + 1) then
         call fatal_error(error, "3d grid: field block starts outside the grid")
         return
      end if
      if (size(f) > self%ngrid - first + 1) then
         call fatal_error(error, "3d grid: field block extends beyond the grid")
         return
      end if
      if (size(f) == 0) return
      call compensated_block_sum(f, 1, size(f), partial(1), partial(2), self%w(first:first + size(f) - 1))
   end subroutine grid_integrate_field_block

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
