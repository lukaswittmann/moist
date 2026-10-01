!> Radial Fourier-Bessel transforms between paired grids
!>
!> Copied nodes, weights, operator, and scratch; no grid pointers
!>
!> One instance per thread; independent clones via `allocate(t, source=template)`
!>
!> Call shapes for each direction:
!>
!>   * scalar `fbt_*`: one field; caller loops or parallelizes
!>   * batched `fbt_*_all`: all columns of a 2D field
!>     Base loops scalar calls; concretes use batched DST or GEMM
!>
!> DST-IV on a uniform radial pair:
!>
!> - Nodes: r_i = (i - 0.5)*dr, k_j = (j - 0.5)*dk,
!>   dk = pi/(npts*dr)
!> - Forward: diag(2*pi*dr/k) . S . diag(r)
!> - Backward: diag(dk/(4*pi^2)/r) . S . diag(k)
!> - S: unnormalised DST-IV (RODFT11); k2r(r2k(f)) = f to roundoff
!> - pi = acos(-1), matching `new_uniform_radial_pair` and exact dk check
!>
!> General quadrature on any radial pair: dense spherical Fourier kernel
!> Stored dr/dk weights w_r and w_k:
!>
!>     F(k_j) = 4*pi * sum_i f(r_i) r_i^2 sinc(k_j r_i) w_r,i
!>     f(r_i) = 1/(2*pi^2) * sum_j F(k_j) k_j^2 sinc(k_j r_i) w_k,j
!>
!> - sinc(z) = sin(z)/z; sinc(0) = 1
!> - Adjoints: matrix transposes; no exact discrete inversion
!> - Accuracy limited to k nodes resolved by r spacing; no band-limit check
!> - Chebyshev-II + Becke outer k nodes remain unresolved
!> - For r, k <= 10: ~1e-11 at 128 nodes, roundoff at 200
!>   (Gaussian exponent 0.7; p_r = 1, p_k = 1.5)
!> - O(nr*nk) work and storage for both dense matrices
!>
!> `new_radial_trafo` dispatch by grid-pair tags:
!>
!> - Two `transform_dst4` tags, equal npts and spacings: DST-IV
!> - Two `transform_quadrature` tags: general quadrature
!> - Mixed, mismatched, unset (0), or unknown tags: error
!>
!> Dispatch uses tags and recorded DST-IV geometry; ignores node values
module moist_math_grid_radial_trafo
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use moist_math_fft_dst4, only: dst4_plan_type, dst4_work_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, transform_dst4, &
      & transform_quadrature
   implicit none(type, external)
   private

   public :: moist_math_grid_radial_trafo_type
   public :: check_trafo_sizes
   public :: check_trafo_shape
   public :: moist_math_grid_radial_trafo_dst4_type
   public :: new_dst4_trafo
   public :: moist_math_grid_radial_trafo_quadrature_type
   public :: new_quadrature_trafo
   public :: new_radial_trafo

   !> Abstract thread-local radial Fourier-Bessel transform
   type, abstract :: moist_math_grid_radial_trafo_type
      !> Number of r-space nodes
      integer :: nr = 0
      !> Number of k-space nodes
      integer :: nk = 0
   contains
      !> Implementation tag: transform_dst4 or transform_quadrature
      procedure(radial_trafo_implementation_i), deferred :: implementation
      !> Forward Fourier-Bessel transform (r -> k)
      procedure(radial_trafo_fbt_i), deferred :: fbt_r2k
      !> Backward Fourier-Bessel transform (k -> r)
      procedure(radial_trafo_fbt_i), deferred :: fbt_k2r
      !> Forward Euclidean adjoint (k -> r)
      procedure(radial_trafo_fbt_i), deferred :: fbt_r2k_adj
      !> Backward Euclidean adjoint (r -> k)
      procedure(radial_trafo_fbt_i), deferred :: fbt_k2r_adj
      !> Batched forward transform of 2D columns (r -> k)
      procedure :: fbt_r2k_all => radial_trafo_fbt_r2k_all_default
      !> Batched backward transform of 2D columns (k -> r)
      procedure :: fbt_k2r_all => radial_trafo_fbt_k2r_all_default
      !> Batched forward adjoint (k -> r)
      procedure :: fbt_r2k_adj_all => radial_trafo_fbt_r2k_adj_all_default
      !> Batched backward adjoint (r -> k)
      procedure :: fbt_k2r_adj_all => radial_trafo_fbt_k2r_adj_all_default
   end type moist_math_grid_radial_trafo_type

   !> Deferred radial-transform operations
   abstract interface
      !> Implementation transform tag
      !>
      !> @param[in] self  trafo instance
      pure function radial_trafo_implementation_i(self) result(tag)
         import :: moist_math_grid_radial_trafo_type
         implicit none(type, external)
         !> Trafo instance
         class(moist_math_grid_radial_trafo_type), intent(in) :: self
         !> transform_dst4 or transform_quadrature
         integer :: tag
      end function radial_trafo_implementation_i

      !> Transform `f_in` to `f_out`
      !>
      !> Forward, backward, or Euclidean adjoint by binding
      !> Thread-local `self` owns scratch
      !>
      !> @param[in,out] self   trafo instance (thread-local)
      !> @param[in]     f_in   input field, size nr (r-space) or nk (k-space)
      !> @param[out]    f_out  output field, size nk (k-space) or nr (r-space)
      !> @param[out]    error  set on a size mismatch or a backend failure
      subroutine radial_trafo_fbt_i(self, f_in, f_out, error)
         import :: moist_math_grid_radial_trafo_type, wp, error_type
         implicit none(type, external)
         !> Trafo instance (thread-local)
         class(moist_math_grid_radial_trafo_type), intent(inout) :: self
         !> Input field
         real(wp), intent(in) :: f_in(:)
         !> Output field
         real(wp), intent(out) :: f_out(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine radial_trafo_fbt_i
   end interface

   !> DST-IV Fourier-Bessel transform
   type, extends(moist_math_grid_radial_trafo_type) :: moist_math_grid_radial_trafo_dst4_type
      !> r-space spacing (bohr), read-only for callers
      real(wp) :: dr = 0.0_wp
      !> k-space spacing (1/bohr), read-only for callers
      real(wp) :: dk = 0.0_wp
      !> pi = acos(-1) for DST-IV prefactors
      real(wp) :: pi = 0.0_wp
      !> r-space nodes (bohr), shape (npts)
      real(wp), allocatable :: r(:)
      !> k-space nodes (1/bohr), shape (npts)
      real(wp), allocatable :: k(:)
      !> Precomputed 1/r, shape (npts)
      real(wp), allocatable :: inv_r(:)
      !> Precomputed 1/k, shape (npts)
      real(wp), allocatable :: inv_k(:)
      !> Single-column DST-IV plan (geometry only)
      type(dst4_plan_type) :: dst4
      !> Scratch for the single-column plan
      type(dst4_work_type) :: dst4_work
      !> Single-column work, shape (npts)
      real(wp), allocatable :: work(:)
      !> Single-column output, shape (npts)
      real(wp), allocatable :: tmp(:)
      !> DST-IV plan rebuilt for each batch width
      type(dst4_plan_type) :: dst4_all
      !> Scratch for the batched plan
      type(dst4_work_type) :: dst4_work_all
      !> Planned batch width (0 = none)
      integer :: nbatch = 0
      !> Batched work buffer, shape (npts, nbatch)
      real(wp), allocatable :: work_all(:, :)
      !> Batched output buffer, shape (npts, nbatch)
      real(wp), allocatable :: tmp_all(:, :)
   contains
      !> Transform tag: transform_dst4
      procedure :: implementation => dst4_trafo_implementation
      !> Forward transform (r -> k)
      procedure :: fbt_r2k => dst4_trafo_fbt_r2k
      !> Backward transform (k -> r)
      procedure :: fbt_k2r => dst4_trafo_fbt_k2r
      !> Forward adjoint (k -> r)
      procedure :: fbt_r2k_adj => dst4_trafo_fbt_r2k_adj
      !> Backward adjoint (r -> k)
      procedure :: fbt_k2r_adj => dst4_trafo_fbt_k2r_adj
      !> Batched forward transform (r -> k)
      procedure :: fbt_r2k_all => dst4_trafo_fbt_r2k_all
      !> Batched backward transform (k -> r)
      procedure :: fbt_k2r_all => dst4_trafo_fbt_k2r_all
      !> Batched forward adjoint (k -> r)
      procedure :: fbt_r2k_adj_all => dst4_trafo_fbt_r2k_adj_all
      !> Batched backward adjoint (r -> k)
      procedure :: fbt_k2r_adj_all => dst4_trafo_fbt_k2r_adj_all
   end type moist_math_grid_radial_trafo_dst4_type

   !> 4*pi, the forward-transform normalisation
   real(wp), parameter :: four_pi = 4.0_wp*pi
   !> 1/(2*pi^2), the backward-transform normalisation
   real(wp), parameter :: inv_two_pi2 = 1.0_wp/(2.0_wp*pi*pi)

   !> Dense Fourier-Bessel transform
   type, extends(moist_math_grid_radial_trafo_type) :: moist_math_grid_radial_trafo_quadrature_type
      !> Forward matrix, 4*pi*r_i^2*w_r,i*sinc(k_j r_i), shape (nk, nr)
      real(wp), allocatable :: fwd(:, :)
      !> Backward matrix, k_j^2*w_k,j*sinc(k_j r_i)/(2*pi^2), shape (nr, nk)
      real(wp), allocatable :: bwd(:, :)
   contains
      !> Transform tag: transform_quadrature
      procedure :: implementation => quadrature_trafo_implementation
      !> Forward transform (r -> k)
      procedure :: fbt_r2k => quadrature_trafo_fbt_r2k
      !> Backward transform (k -> r)
      procedure :: fbt_k2r => quadrature_trafo_fbt_k2r
      !> Forward adjoint (k -> r)
      procedure :: fbt_r2k_adj => quadrature_trafo_fbt_r2k_adj
      !> Backward adjoint (r -> k)
      procedure :: fbt_k2r_adj => quadrature_trafo_fbt_k2r_adj
      !> Batched forward transform (r -> k)
      procedure :: fbt_r2k_all => quadrature_trafo_fbt_r2k_all
      !> Batched backward transform (k -> r)
      procedure :: fbt_k2r_all => quadrature_trafo_fbt_k2r_all
      !> Batched forward adjoint (k -> r)
      procedure :: fbt_r2k_adj_all => quadrature_trafo_fbt_r2k_adj_all
      !> Batched backward adjoint (r -> k)
      procedure :: fbt_k2r_adj_all => quadrature_trafo_fbt_k2r_adj_all
   end type moist_math_grid_radial_trafo_quadrature_type

contains

   !> Check transform input and output lengths
   !>
   !> @param[in]  n_in      input length
   !> @param[in]  n_out     output length
   !> @param[in]  want_in   expected input length
   !> @param[in]  want_out  expected output length
   !> @param[in]  label     procedure name for error
   !> @param[out] error     size mismatch
   subroutine check_trafo_sizes(n_in, n_out, want_in, want_out, label, error)
      !> Size of the input field
      integer, intent(in) :: n_in
      !> Size of the output field
      integer, intent(in) :: n_out
      !> Expected input size
      integer, intent(in) :: want_in
      !> Expected output size
      integer, intent(in) :: want_out
      !> Procedure name used in the message
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (n_in /= want_in) then
         call fatal_error(error, "Radial trafo "//label//": input size does not match the grid")
         return
      end if
      if (n_out /= want_out) then
         call fatal_error(error, "Radial trafo "//label//": output size does not match the grid")
      end if
   end subroutine check_trafo_sizes

   !> Check batched transform shapes
   !>
   !> @param[in]  shape_in   input shape (n_in, nbatch)
   !> @param[in]  shape_out  output shape (n_out, nbatch)
   !> @param[in]  want_in    expected input column length
   !> @param[in]  want_out   expected output column length
   !> @param[in]  label      procedure name for error
   !> @param[out] error      column-length or batch-width mismatch
   subroutine check_trafo_shape(shape_in, shape_out, want_in, want_out, label, error)
      !> Shape of the input fields
      integer, intent(in) :: shape_in(2)
      !> Shape of the output fields
      integer, intent(in) :: shape_out(2)
      !> Expected input column length
      integer, intent(in) :: want_in
      !> Expected output column length
      integer, intent(in) :: want_out
      !> Procedure name used in the message
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (shape_in(2) /= shape_out(2)) then
         call fatal_error(error, "Radial trafo "//label//": input and output batch widths differ")
         return
      end if
      call check_trafo_sizes(shape_in(1), shape_out(1), want_in, want_out, label, error)
   end subroutine check_trafo_shape

   !> Apply scalar forward transform per column
   !>
   !> Stop at first failed column
   !>
   !> @param[in,out] self   trafo instance (thread-local)
   !> @param[in]     f_in   r-space fields, shape (nr, nbatch)
   !> @param[out]    f_out  k-space fields, shape (nk, nbatch)
   !> @param[out]    error  set on a size mismatch or a failed column
   subroutine radial_trafo_fbt_r2k_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_radial_trafo_type), intent(inout) :: self
      !> r-space fields, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> k-space fields, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      call check_trafo_shape(shape(f_in), shape(f_out), self%nr, self%nk, "fbt_r2k_all", error)
      if (allocated(error)) return
      do j = 1, size(f_in, 2)
         call self%fbt_r2k(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_r2k_all_default

   !> Apply scalar backward transform per column
   !>
   !> Stop at first failed column
   !>
   !> @param[in,out] self   trafo instance (thread-local)
   !> @param[in]     f_in   k-space fields, shape (nk, nbatch)
   !> @param[out]    f_out  r-space fields, shape (nr, nbatch)
   !> @param[out]    error  set on a size mismatch or a failed column
   subroutine radial_trafo_fbt_k2r_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_radial_trafo_type), intent(inout) :: self
      !> k-space fields, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> r-space fields, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      call check_trafo_shape(shape(f_in), shape(f_out), self%nk, self%nr, "fbt_k2r_all", error)
      if (allocated(error)) return
      do j = 1, size(f_in, 2)
         call self%fbt_k2r(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_k2r_all_default

   !> Apply scalar forward adjoint per column
   !>
   !> Stop at first failed column
   !>
   !> @param[in,out] self   trafo instance (thread-local)
   !> @param[in]     f_in   k-space vectors, shape (nk, nbatch)
   !> @param[out]    f_out  r-space vectors, shape (nr, nbatch)
   !> @param[out]    error  set on a size mismatch or a failed column
   subroutine radial_trafo_fbt_r2k_adj_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_radial_trafo_type), intent(inout) :: self
      !> k-space vectors, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> r-space vectors, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      call check_trafo_shape(shape(f_in), shape(f_out), self%nk, self%nr, "fbt_r2k_adj_all", error)
      if (allocated(error)) return
      do j = 1, size(f_in, 2)
         call self%fbt_r2k_adj(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_r2k_adj_all_default

   !> Apply scalar backward adjoint per column
   !>
   !> Stop at first failed column
   !>
   !> @param[in,out] self   trafo instance (thread-local)
   !> @param[in]     f_in   r-space vectors, shape (nr, nbatch)
   !> @param[out]    f_out  k-space vectors, shape (nk, nbatch)
   !> @param[out]    error  set on a size mismatch or a failed column
   subroutine radial_trafo_fbt_k2r_adj_all_default(self, f_in, f_out, error)
      !> Trafo instance (thread-local)
      class(moist_math_grid_radial_trafo_type), intent(inout) :: self
      !> r-space vectors, one per column
      real(wp), intent(in) :: f_in(:, :)
      !> k-space vectors, one per column
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: j

      call check_trafo_shape(shape(f_in), shape(f_out), self%nr, self%nk, "fbt_k2r_adj_all", error)
      if (allocated(error)) return
      do j = 1, size(f_in, 2)
         call self%fbt_k2r_adj(f_in(:, j), f_out(:, j), error)
         if (allocated(error)) return
      end do
   end subroutine radial_trafo_fbt_k2r_adj_all_default

   !> Build a DST-IV transform from a uniform radial pair
   !>
   !> Require two transform_dst4 tags, equal npts, and exact
   !> dk = acos(-1)/(npts*dr) from `new_uniform_radial_pair`
   !> Copy nodes; retain recorded spacings
   !>
   !> @param[out] self   new trafo
   !> @param[in]  rgrid  r-space grid, tagged transform_dst4
   !> @param[in]  kgrid  k-space grid, tagged transform_dst4
   !> @param[out] error  set on a tag, size, or spacing mismatch, or allocation failure
   subroutine new_dst4_trafo(self, rgrid, kgrid, error)
      !> New trafo
      type(moist_math_grid_radial_trafo_dst4_type), intent(out) :: self
      !> r-space grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i, npts, stat
      real(wp) :: pi_dst

      if (rgrid%transform /= transform_dst4 .or. kgrid%transform /= transform_dst4) then
         call fatal_error(error, "DST-IV radial trafo: both grids must be tagged transform_dst4")
         return
      end if
      if (rgrid%npts /= kgrid%npts) then
         call fatal_error(error, "DST-IV radial trafo: r and k grids differ in npts")
         return
      end if
      npts = rgrid%npts
      if (npts < 1) then
         call fatal_error(error, "DST-IV radial trafo: grids have no nodes")
         return
      end if
      if (.not. (allocated(rgrid%r) .and. allocated(kgrid%r))) then
         call fatal_error(error, "DST-IV radial trafo: grid nodes are not allocated")
         return
      end if
      if (size(rgrid%r) /= npts .or. size(kgrid%r) /= npts) then
         call fatal_error(error, "DST-IV radial trafo: node arrays do not match npts")
         return
      end if
      if (.not. (ieee_is_finite(rgrid%spacing) .and. rgrid%spacing > 0.0_wp)) then
         call fatal_error(error, "DST-IV radial trafo: r grid has no positive finite spacing")
         return
      end if
      ! Match pair constructor expression for exact equality
      pi_dst = acos(-1.0_wp)
      if (.not. kgrid%spacing == pi_dst/(real(npts, wp)*rgrid%spacing)) then
         call fatal_error(error, "DST-IV radial trafo: k spacing does not match pi/(npts*dr)")
         return
      end if

      self%nr = npts
      self%nk = npts
      self%pi = pi_dst
      self%dr = rgrid%spacing
      self%dk = kgrid%spacing

      allocate (self%r(npts), stat=stat)
      if (stat == 0) allocate (self%k(npts), stat=stat)
      if (stat == 0) allocate (self%inv_r(npts), stat=stat)
      if (stat == 0) allocate (self%inv_k(npts), stat=stat)
      if (stat == 0) allocate (self%work(npts), stat=stat)
      if (stat == 0) allocate (self%tmp(npts), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "DST-IV radial trafo: failed to allocate nodes and buffers")
         return
      end if
      do i = 1, npts
         self%r(i) = rgrid%r(i)
         self%k(i) = kgrid%r(i)
         self%inv_r(i) = 1.0_wp/self%r(i)
         self%inv_k(i) = 1.0_wp/self%k(i)
      end do

      call self%dst4%init(npts, 1)
      call self%dst4%new_work(self%dst4_work)
   end subroutine new_dst4_trafo

   !> Implementation transform tag
   !>
   !> @param[in] self  trafo instance
   pure function dst4_trafo_implementation(self) result(tag)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(in) :: self
      !> transform_dst4
      integer :: tag

      tag = transform_dst4
   end function dst4_trafo_implementation

   !> Forward Fourier-Bessel transform r -> k
   !>
   !> @param[in,out] self   trafo instance (owns the scratch)
   !> @param[in]     f_in   function in r-space, size npts
   !> @param[out]    f_out  function in k-space, size npts
   !> @param[out]    error  set on a size mismatch or a backend failure
   subroutine dst4_trafo_fbt_r2k(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> Function in r-space
      real(wp), intent(in) :: f_in(:)
      !> Function in k-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac

      call check_trafo_sizes(size(f_in), size(f_out), self%nr, self%nk, "fbt_r2k", error)
      if (allocated(error)) return
      prefac = 2.0_wp*self%pi*self%dr
      self%work = f_in*self%r
      call self%dst4%execute(self%dst4_work, self%work, self%tmp, error)
      if (allocated(error)) return
      f_out = prefac*self%tmp*self%inv_k
   end subroutine dst4_trafo_fbt_r2k

   !> Backward Fourier-Bessel transform k -> r
   !>
   !> @param[in,out] self   trafo instance (owns the scratch)
   !> @param[in]     f_in   function in k-space, size npts
   !> @param[out]    f_out  function in r-space, size npts
   !> @param[out]    error  set on a size mismatch or a backend failure
   subroutine dst4_trafo_fbt_k2r(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> Function in k-space
      real(wp), intent(in) :: f_in(:)
      !> Function in r-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac

      call check_trafo_sizes(size(f_in), size(f_out), self%nk, self%nr, "fbt_k2r", error)
      if (allocated(error)) return
      prefac = self%dk/(4.0_wp*self%pi*self%pi)
      self%work = f_in*self%k
      call self%dst4%execute(self%dst4_work, self%work, self%tmp, error)
      if (allocated(error)) return
      f_out = prefac*self%tmp*self%inv_r
   end subroutine dst4_trafo_fbt_k2r

   !> Adjoint of the forward transform, y_r = (FBT_rk)^T x_k
   !>
   !> Symmetric S in FBT_rk = diag(2*pi*dr/k) . S . diag(r)
   !> Transpose swaps r and k diagonals around the same DST
   !>
   !> @param[in,out] self   trafo instance (owns the scratch)
   !> @param[in]     f_in   vector in k-space, size npts
   !> @param[out]    f_out  vector in r-space, size npts
   !> @param[out]    error  set on a size mismatch or a backend failure
   subroutine dst4_trafo_fbt_r2k_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> Vector in k-space
      real(wp), intent(in) :: f_in(:)
      !> Vector in r-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac

      call check_trafo_sizes(size(f_in), size(f_out), self%nk, self%nr, "fbt_r2k_adj", error)
      if (allocated(error)) return
      prefac = 2.0_wp*self%pi*self%dr
      self%work = prefac*f_in*self%inv_k
      call self%dst4%execute(self%dst4_work, self%work, self%tmp, error)
      if (allocated(error)) return
      f_out = self%r*self%tmp
   end subroutine dst4_trafo_fbt_r2k_adj

   !> Adjoint of the backward transform, y_k = (FBT_kr)^T x_r
   !>
   !> FBT_kr = diag(dk/(4*pi^2)/r) . S . diag(k)
   !> Transpose swaps r and k diagonals around the same DST
   !>
   !> @param[in,out] self   trafo instance (owns the scratch)
   !> @param[in]     f_in   vector in r-space, size npts
   !> @param[out]    f_out  vector in k-space, size npts
   !> @param[out]    error  set on a size mismatch or a backend failure
   subroutine dst4_trafo_fbt_k2r_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> Vector in r-space
      real(wp), intent(in) :: f_in(:)
      !> Vector in k-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac

      call check_trafo_sizes(size(f_in), size(f_out), self%nr, self%nk, "fbt_k2r_adj", error)
      if (allocated(error)) return
      prefac = self%dk/(4.0_wp*self%pi*self%pi)
      self%work = prefac*f_in*self%inv_r
      call self%dst4%execute(self%dst4_work, self%work, self%tmp, error)
      if (allocated(error)) return
      f_out = self%k*self%tmp
   end subroutine dst4_trafo_fbt_k2r_adj

   !> Size batched plan and buffers for `nbatch` columns
   !>
   !> Rebuild only on width change; one plan serves four directions
   !> Geometry-only plan; buffer allocation is the failure point
   !>
   !> @param[in,out] self    trafo instance (owns the batched plan and buffers)
   !> @param[in]     nbatch  number of columns, >= 1
   !> @param[out]    error   set if the batched buffers cannot be allocated
   subroutine dst4_trafo_ensure_all(self, nbatch, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> Number of columns
      integer, intent(in) :: nbatch
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: stat

      if (self%nbatch == nbatch) return

      call self%dst4_all%destroy()
      call self%dst4_work_all%destroy()
      self%nbatch = 0
      if (allocated(self%work_all)) deallocate (self%work_all)
      if (allocated(self%tmp_all)) deallocate (self%tmp_all)
      allocate (self%work_all(self%nr, nbatch), stat=stat)
      if (stat == 0) allocate (self%tmp_all(self%nr, nbatch), stat=stat)
      if (stat /= 0) then
         if (allocated(self%work_all)) deallocate (self%work_all)
         call fatal_error(error, "DST-IV radial trafo: batched buffer allocation failed")
         return
      end if

      call self%dst4_all%init(self%nr, nbatch)
      call self%dst4_all%new_work(self%dst4_work_all)
      self%nbatch = nbatch
   end subroutine dst4_trafo_ensure_all

   !> Batched forward transform r -> k over all columns
   !>
   !> Empty batch leaves cached buffers untouched
   !>
   !> @param[in,out] self   trafo instance (owns the batched plan and buffers)
   !> @param[in]     f_in   r-space fields, shape (npts, nbatch)
   !> @param[out]    f_out  k-space fields, shape (npts, nbatch)
   !> @param[out]    error  set on a shape mismatch or a backend failure
   subroutine dst4_trafo_fbt_r2k_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> r-space fields
      real(wp), intent(in) :: f_in(:, :)
      !> k-space fields
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac
      integer :: j, nb

      call check_trafo_shape(shape(f_in), shape(f_out), self%nr, self%nk, "fbt_r2k_all", error)
      if (allocated(error)) return
      nb = size(f_in, 2)
      if (nb == 0) return
      call dst4_trafo_ensure_all(self, nb, error)
      if (allocated(error)) return
      prefac = 2.0_wp*self%pi*self%dr
      do j = 1, nb
         self%work_all(:, j) = f_in(:, j)*self%r
      end do
      call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
      if (allocated(error)) return
      do j = 1, nb
         f_out(:, j) = prefac*self%tmp_all(:, j)*self%inv_k
      end do
   end subroutine dst4_trafo_fbt_r2k_all

   !> Batched backward transform k -> r over all columns
   !>
   !> @param[in,out] self   trafo instance (owns the batched plan and buffers)
   !> @param[in]     f_in   k-space fields, shape (npts, nbatch)
   !> @param[out]    f_out  r-space fields, shape (npts, nbatch)
   !> @param[out]    error  set on a shape mismatch or a backend failure
   subroutine dst4_trafo_fbt_k2r_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> k-space fields
      real(wp), intent(in) :: f_in(:, :)
      !> r-space fields
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac
      integer :: j, nb

      call check_trafo_shape(shape(f_in), shape(f_out), self%nk, self%nr, "fbt_k2r_all", error)
      if (allocated(error)) return
      nb = size(f_in, 2)
      if (nb == 0) return
      call dst4_trafo_ensure_all(self, nb, error)
      if (allocated(error)) return
      prefac = self%dk/(4.0_wp*self%pi*self%pi)
      do j = 1, nb
         self%work_all(:, j) = f_in(:, j)*self%k
      end do
      call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
      if (allocated(error)) return
      do j = 1, nb
         f_out(:, j) = prefac*self%tmp_all(:, j)*self%inv_r
      end do
   end subroutine dst4_trafo_fbt_k2r_all

   !> Batched adjoint of the forward transform over all columns
   !>
   !> @param[in,out] self   trafo instance (owns the batched plan and buffers)
   !> @param[in]     f_in   k-space vectors, shape (npts, nbatch)
   !> @param[out]    f_out  r-space vectors, shape (npts, nbatch)
   !> @param[out]    error  set on a shape mismatch or a backend failure
   subroutine dst4_trafo_fbt_r2k_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> k-space vectors
      real(wp), intent(in) :: f_in(:, :)
      !> r-space vectors
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac
      integer :: j, nb

      call check_trafo_shape(shape(f_in), shape(f_out), self%nk, self%nr, "fbt_r2k_adj_all", error)
      if (allocated(error)) return
      nb = size(f_in, 2)
      if (nb == 0) return
      call dst4_trafo_ensure_all(self, nb, error)
      if (allocated(error)) return
      prefac = 2.0_wp*self%pi*self%dr
      do j = 1, nb
         self%work_all(:, j) = prefac*f_in(:, j)*self%inv_k
      end do
      call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
      if (allocated(error)) return
      do j = 1, nb
         f_out(:, j) = self%r*self%tmp_all(:, j)
      end do
   end subroutine dst4_trafo_fbt_r2k_adj_all

   !> Batched adjoint of the backward transform over all columns
   !>
   !> @param[in,out] self   trafo instance (owns the batched plan and buffers)
   !> @param[in]     f_in   r-space vectors, shape (npts, nbatch)
   !> @param[out]    f_out  k-space vectors, shape (npts, nbatch)
   !> @param[out]    error  set on a shape mismatch or a backend failure
   subroutine dst4_trafo_fbt_k2r_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_dst4_type), intent(inout) :: self
      !> r-space vectors
      real(wp), intent(in) :: f_in(:, :)
      !> k-space vectors
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: prefac
      integer :: j, nb

      call check_trafo_shape(shape(f_in), shape(f_out), self%nr, self%nk, "fbt_k2r_adj_all", error)
      if (allocated(error)) return
      nb = size(f_in, 2)
      if (nb == 0) return
      call dst4_trafo_ensure_all(self, nb, error)
      if (allocated(error)) return
      prefac = self%dk/(4.0_wp*self%pi*self%pi)
      do j = 1, nb
         self%work_all(:, j) = prefac*f_in(:, j)*self%inv_r
      end do
      call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
      if (allocated(error)) return
      do j = 1, nb
         f_out(:, j) = self%k*self%tmp_all(:, j)
      end do
   end subroutine dst4_trafo_fbt_k2r_adj_all

   !> Build a general-quadrature transform on a radial pair
   !>
   !> Accept any built pair, including unequal node counts and tags
   !>
   !> @param[out] self   new trafo
   !> @param[in]  rgrid  r-space grid, nodes in bohr, dr weights
   !> @param[in]  kgrid  k-space grid, nodes in 1/bohr, dk weights
   !> @param[out] error  set on an empty or inconsistent grid, or allocation failure
   subroutine new_quadrature_trafo(self, rgrid, kgrid, error)
      !> New trafo
      type(moist_math_grid_radial_trafo_quadrature_type), intent(out) :: self
      !> r-space grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: ar(:), bk(:)
      real(wp) :: s
      integer :: i, j, nr, nk, stat

      call check_grid(rgrid, "r", error)
      if (allocated(error)) return
      call check_grid(kgrid, "k", error)
      if (allocated(error)) return
      nr = rgrid%npts
      nk = kgrid%npts

      allocate (self%fwd(nk, nr), stat=stat)
      if (stat == 0) allocate (self%bwd(nr, nk), stat=stat)
      if (stat == 0) allocate (ar(nr), stat=stat)
      if (stat == 0) allocate (bk(nk), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Quadrature radial trafo: failed to allocate transform matrices")
         return
      end if
      self%nr = nr
      self%nk = nk

      do i = 1, nr
         ar(i) = four_pi*rgrid%r(i)*rgrid%r(i)*rgrid%w(i)
      end do
      do j = 1, nk
         bk(j) = inv_two_pi2*kgrid%r(j)*kgrid%r(j)*kgrid%w(j)
      end do
      do i = 1, nr
         do j = 1, nk
            s = sinc(kgrid%r(j)*rgrid%r(i))
            self%fwd(j, i) = ar(i)*s
            self%bwd(i, j) = bk(j)*s
         end do
      end do
   end subroutine new_quadrature_trafo

   !> Check built grid and npts-sized arrays
   !>
   !> @param[in]  grid   radial grid
   !> @param[in]  label  "r" or "k", used in the message
   !> @param[out] error  set on an unbuilt or inconsistent grid
   subroutine check_grid(grid, label, error)
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Grid label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (grid%npts < 1) then
         call fatal_error(error, "Quadrature radial trafo: "//label//" grid has no nodes")
         return
      end if
      if (.not. (allocated(grid%r) .and. allocated(grid%w))) then
         call fatal_error(error, "Quadrature radial trafo: "//label//" grid arrays are not allocated")
         return
      end if
      if (size(grid%r) /= grid%npts .or. size(grid%w) /= grid%npts) then
         call fatal_error(error, "Quadrature radial trafo: "//label//" grid arrays do not match npts")
      end if
   end subroutine check_grid

   !> sin(z)/z with sinc(0) = 1
   !>
   !> @param[in] z  argument
   elemental function sinc(z) result(s)
      !> Argument
      real(wp), intent(in) :: z
      !> sin(z)/z
      real(wp) :: s

      if (z == 0.0_wp) then
         s = 1.0_wp
      else
         s = sin(z)/z
      end if
   end function sinc

   !> Transform tag of this implementation
   !>
   !> @param[in] self  trafo instance
   pure function quadrature_trafo_implementation(self) result(tag)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(in) :: self
      !> transform_quadrature
      integer :: tag

      tag = transform_quadrature
   end function quadrature_trafo_implementation

   !> Forward transform r -> k, F = M f
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   function in r-space, size nr
   !> @param[out]    f_out  function in k-space, size nk
   !> @param[out]    error  set on a size mismatch
   subroutine quadrature_trafo_fbt_r2k(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> Function in r-space
      real(wp), intent(in) :: f_in(:)
      !> Function in k-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_sizes(size(f_in), size(f_out), self%nr, self%nk, "fbt_r2k", error)
      if (allocated(error)) return
      f_out = matmul(self%fwd, f_in)
   end subroutine quadrature_trafo_fbt_r2k

   !> Backward transform k -> r, f = B F
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   function in k-space, size nk
   !> @param[out]    f_out  function in r-space, size nr
   !> @param[out]    error  set on a size mismatch
   subroutine quadrature_trafo_fbt_k2r(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> Function in k-space
      real(wp), intent(in) :: f_in(:)
      !> Function in r-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_sizes(size(f_in), size(f_out), self%nk, self%nr, "fbt_k2r", error)
      if (allocated(error)) return
      f_out = matmul(self%bwd, f_in)
   end subroutine quadrature_trafo_fbt_k2r

   !> Adjoint of the forward transform, y_r = M^T x_k
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   vector in k-space, size nk
   !> @param[out]    f_out  vector in r-space, size nr
   !> @param[out]    error  set on a size mismatch
   subroutine quadrature_trafo_fbt_r2k_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> Vector in k-space
      real(wp), intent(in) :: f_in(:)
      !> Vector in r-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_sizes(size(f_in), size(f_out), self%nk, self%nr, "fbt_r2k_adj", error)
      if (allocated(error)) return
      f_out = matmul(f_in, self%fwd)
   end subroutine quadrature_trafo_fbt_r2k_adj

   !> Adjoint of the backward transform, y_k = B^T x_r
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   vector in r-space, size nr
   !> @param[out]    f_out  vector in k-space, size nk
   !> @param[out]    error  set on a size mismatch
   subroutine quadrature_trafo_fbt_k2r_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> Vector in r-space
      real(wp), intent(in) :: f_in(:)
      !> Vector in k-space
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_sizes(size(f_in), size(f_out), self%nr, self%nk, "fbt_k2r_adj", error)
      if (allocated(error)) return
      f_out = matmul(f_in, self%bwd)
   end subroutine quadrature_trafo_fbt_k2r_adj

   !> Batched forward transform as one matrix product, F = M f
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   r-space fields, shape (nr, nbatch)
   !> @param[out]    f_out  k-space fields, shape (nk, nbatch)
   !> @param[out]    error  set on a shape mismatch
   subroutine quadrature_trafo_fbt_r2k_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> r-space fields
      real(wp), intent(in) :: f_in(:, :)
      !> k-space fields
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_shape(shape(f_in), shape(f_out), self%nr, self%nk, "fbt_r2k_all", error)
      if (allocated(error)) return
      f_out = matmul(self%fwd, f_in)
   end subroutine quadrature_trafo_fbt_r2k_all

   !> Batched backward transform as one matrix product, f = B F
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   k-space fields, shape (nk, nbatch)
   !> @param[out]    f_out  r-space fields, shape (nr, nbatch)
   !> @param[out]    error  set on a shape mismatch
   subroutine quadrature_trafo_fbt_k2r_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> k-space fields
      real(wp), intent(in) :: f_in(:, :)
      !> r-space fields
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_shape(shape(f_in), shape(f_out), self%nk, self%nr, "fbt_k2r_all", error)
      if (allocated(error)) return
      f_out = matmul(self%bwd, f_in)
   end subroutine quadrature_trafo_fbt_k2r_all

   !> Batched adjoint of the forward transform, y_r = M^T x_k
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   k-space vectors, shape (nk, nbatch)
   !> @param[out]    f_out  r-space vectors, shape (nr, nbatch)
   !> @param[out]    error  set on a shape mismatch
   subroutine quadrature_trafo_fbt_r2k_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> k-space vectors
      real(wp), intent(in) :: f_in(:, :)
      !> r-space vectors
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_shape(shape(f_in), shape(f_out), self%nk, self%nr, "fbt_r2k_adj_all", error)
      if (allocated(error)) return
      f_out = matmul(transpose(self%fwd), f_in)
   end subroutine quadrature_trafo_fbt_r2k_adj_all

   !> Batched adjoint of the backward transform, y_k = B^T x_r
   !>
   !> @param[in,out] self   trafo instance
   !> @param[in]     f_in   r-space vectors, shape (nr, nbatch)
   !> @param[out]    f_out  k-space vectors, shape (nk, nbatch)
   !> @param[out]    error  set on a shape mismatch
   subroutine quadrature_trafo_fbt_k2r_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_quadrature_type), intent(inout) :: self
      !> r-space vectors
      real(wp), intent(in) :: f_in(:, :)
      !> k-space vectors
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_trafo_shape(shape(f_in), shape(f_out), self%nr, self%nk, "fbt_k2r_adj_all", error)
      if (allocated(error)) return
      f_out = matmul(transpose(self%bwd), f_in)
   end subroutine quadrature_trafo_fbt_k2r_adj_all

   !> Create a radial transform from grid-pair tags
   !>
   !> - Matching transform_dst4 tags, npts, and recorded spacing: DST-IV
   !> - Two transform_quadrature tags: general quadrature
   !> - Other combinations: error
   !> Copies needed grid data; source grids may be released
   !>
   !> Quadrature accuracy limited to k nodes resolved by r spacing
   !> No band-limit check
   !>
   !> @param[out] trafo  new trafo, one per thread; clone with allocate(source=)
   !> @param[in]  rgrid  r-space grid
   !> @param[in]  kgrid  k-space grid
   !> @param[out] error  set on an unset, unknown, or mixed tag pair, a
   !>                    mismatched DST-IV pair, or a construction failure
   subroutine new_radial_trafo(trafo, rgrid, kgrid, error)
      !> New trafo
      class(moist_math_grid_radial_trafo_type), allocatable, intent(out) :: trafo
      !> r-space grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_trafo_dst4_type), allocatable :: dst4
      type(moist_math_grid_radial_trafo_quadrature_type), allocatable :: quad

      call check_tag(rgrid%transform, "r", error)
      if (allocated(error)) return
      call check_tag(kgrid%transform, "k", error)
      if (allocated(error)) return

      if (rgrid%transform == transform_dst4 .and. kgrid%transform == transform_dst4) then
         allocate (dst4)
         call new_dst4_trafo(dst4, rgrid, kgrid, error)
         if (allocated(error)) return
         call move_alloc(dst4, trafo)
      else if (rgrid%transform == transform_quadrature &
         & .and. kgrid%transform == transform_quadrature) then
         allocate (quad)
         call new_quadrature_trafo(quad, rgrid, kgrid, error)
         if (allocated(error)) return
         call move_alloc(quad, trafo)
      else
         call fatal_error(error, "Radial trafo: exactly one grid is tagged for DST-IV; "// &
            & "build both with new_uniform_radial_pair or both with new_radial_grid")
      end if
   end subroutine new_radial_trafo

   !> Reject unset or unknown transform tags
   !>
   !> @param[in]  tag    transform tag of the grid
   !> @param[in]  label  "r" or "k", used in the message
   !> @param[out] error  set on an unset or unknown tag
   subroutine check_tag(tag, label, error)
      !> Transform tag
      integer, intent(in) :: tag
      !> Grid label
      character(len=*), intent(in) :: label
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select case (tag)
      case (transform_quadrature, transform_dst4)
      case (0)
         call fatal_error(error, "Radial trafo: "//label//" grid has no transform tag (grid not built)")
      case default
         call fatal_error(error, "Radial trafo: "//label//" grid has an unknown transform tag")
      end select
   end subroutine check_tag

end module moist_math_grid_radial_trafo
