!> Non-uniform Chebyshev-radii radial grid and its transform engine
module moist_math_grid_1d_chebyshev
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use moist_math_grid_1d_base, only: moist_math_grid_1d_type, &
      & moist_math_grid_1d_trafo_type
   use moist_math_quadrature_chebyshev, only: chebyshev2_radii
   implicit none(type, external)
   private

   public :: moist_math_grid_1d_chebyshev_type, new_chebyshev_radial_grid
   public :: moist_math_grid_1d_chebyshev_trafo_type

   !> 4*pi, the surface area of the unit sphere times the radial Jacobian
   real(wp), parameter :: four_pi = 4.0_wp*pi
   !> 1/(2*pi^2), the inverse-transform normalisation
   real(wp), parameter :: inv_two_pi2 = 1.0_wp/(2.0_wp*pi*pi)

   !> Immutable non-uniform radial grid on Chebyshev-2 nodes in r- and k-space
   !>
   !> Holds the dense forward/backward transform matrices, so it carries no
   !> per-call mutable state and is safe to share read-only across threads
   type, extends(moist_math_grid_1d_type) :: moist_math_grid_1d_chebyshev_type
      !> Number of k-space nodes (npts counts the r-space nodes)
      integer :: nk = 0
      !> r-space quadrature weights, carrying r^2 dr (length npts)
      real(wp), allocatable :: wr(:)
      !> Dense forward transform matrix M, shape (nk, npts)
      real(wp), allocatable :: fwd(:, :)
      !> Dense backward transform matrix B, shape (npts, nk)
      real(wp), allocatable :: bwd(:, :)
   contains
      procedure :: measure => chebyshev_measure
      procedure :: new_trafo => chebyshev_new_trafo
      procedure :: destroy => chebyshev_grid_destroy
   end type moist_math_grid_1d_chebyshev_type

   !> Transform engine for a Chebyshev radial grid
   !>
   !> Stateless apart from the grid pointer (the dense matrices live on the
   !> grid), so it needs no scratch and a single instance is safe to share --
   !> but one per thread is fine too
   type, extends(moist_math_grid_1d_trafo_type) :: moist_math_grid_1d_chebyshev_trafo_type
      !> Grid this trafo transforms on (not owned; must outlive the trafo)
      class(moist_math_grid_1d_chebyshev_type), pointer :: grid => null()
   contains
      procedure :: fbt_r2k => chebyshev_trafo_fbt_r2k
      procedure :: fbt_k2r => chebyshev_trafo_fbt_k2r
      procedure :: fbt_r2k_adj => chebyshev_trafo_fbt_r2k_adj
      procedure :: fbt_k2r_adj => chebyshev_trafo_fbt_k2r_adj
      procedure :: fbt_r2k_all => chebyshev_trafo_fbt_r2k_all
      procedure :: fbt_k2r_all => chebyshev_trafo_fbt_k2r_all
      procedure :: fbt_r2k_adj_all => chebyshev_trafo_fbt_r2k_adj_all
      procedure :: fbt_k2r_adj_all => chebyshev_trafo_fbt_k2r_adj_all
      procedure :: destroy => chebyshev_trafo_destroy
   end type moist_math_grid_1d_chebyshev_trafo_type

contains

   !> Initialize a Chebyshev radial grid and assemble its transform matrices
   !>
   !> @param[out] grid  Initialized grid object
   !> @param[in]  nr  Number of r-space nodes (nr >= 1)
   !> @param[in]  p_r  r-space Chebyshev scale (bohr), finite and positive
   !> @param[in]  nk  Number of k-space nodes (nk >= 1)
   !> @param[in]  p_k  k-space Chebyshev scale (1/bohr), finite and positive
   !> @param[out] error  Set on invalid parameters or allocation failure
   subroutine new_chebyshev_radial_grid(grid, nr, p_r, nk, p_k, error)
      !> Initialized grid object
      type(moist_math_grid_1d_chebyshev_type), intent(out) :: grid
      !> Number of r-space nodes
      integer, intent(in) :: nr
      !> r-space Chebyshev scale (bohr)
      real(wp), intent(in) :: p_r
      !> Number of k-space nodes
      integer, intent(in) :: nk
      !> k-space Chebyshev scale (1/bohr)
      real(wp), intent(in) :: p_k
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i, j, stat
      real(wp), allocatable :: wk(:)

      if (nr < 1 .or. nk < 1) then
         call fatal_error(error, "Chebyshev radial grid needs at least one r- and one k-space node")
         return
      end if
      if (.not. (ieee_is_finite(p_r) .and. ieee_is_finite(p_k))) then
         call fatal_error(error, "Chebyshev radial grid needs finite scales")
         return
      end if
      if (p_r <= 0.0_wp .or. p_k <= 0.0_wp) then
         call fatal_error(error, "Chebyshev radial grid needs positive scales")
         return
      end if

      grid%npts = nr
      grid%nk = nk

      allocate (grid%r(nr), grid%wr(nr), grid%k(nk), wk(nk), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate Chebyshev radial grid nodes")
         return
      end if
      call chebyshev2_radii(nr, p_r, grid%r, grid%wr)
      call chebyshev2_radii(nk, p_k, grid%k, wk)

      allocate (grid%fwd(nk, nr), grid%bwd(nr, nk), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate Chebyshev transform matrices")
         return
      end if
      do i = 1, nr
         do j = 1, nk
            grid%fwd(j, i) = four_pi/grid%k(j)*grid%wr(i)/grid%r(i) &
                             *sin(grid%k(j)*grid%r(i))
            grid%bwd(i, j) = inv_two_pi2/grid%r(i)*wk(j)/grid%k(j) &
                             *sin(grid%k(j)*grid%r(i))
         end do
      end do
   end subroutine new_chebyshev_radial_grid

   !> 3D-radial integration weight of grid point i
   !>
   !> 4*pi*wr_i (the Chebyshev weight already carries r^2 dr)
   !>
   !> @param[in] self  Grid instance
   !> @param[in] i     Grid point index
   pure function chebyshev_measure(self, i) result(w)
      !> Grid instance
      class(moist_math_grid_1d_chebyshev_type), intent(in) :: self
      !> Grid point index
      integer, intent(in) :: i
      !> Weight of point i
      real(wp) :: w

      w = four_pi*self%wr(i)
   end function chebyshev_measure

   !> Release all grid storage; idempotent
   !>
   !> @param[in,out] self  Grid instance
   subroutine chebyshev_grid_destroy(self)
      !> Grid instance
      class(moist_math_grid_1d_chebyshev_type), intent(inout) :: self

      if (allocated(self%r)) deallocate (self%r)
      if (allocated(self%k)) deallocate (self%k)
      if (allocated(self%wr)) deallocate (self%wr)
      if (allocated(self%fwd)) deallocate (self%fwd)
      if (allocated(self%bwd)) deallocate (self%bwd)
      self%npts = 0
      self%nk = 0
   end subroutine chebyshev_grid_destroy

   !> Allocate a `moist_math_grid_1d_chebyshev_trafo_type` bound to `self`
   !>
   !> Dispatched via the base `new_trafo`; the grid owns the dense transform
   !> matrices, so the trafo needs no scratch; the grid must outlive the trafo
   !>
   !> @param[in]  self  Grid instance (must have the target attribute)
   !> @param[out] trafo  Allocated, bound transform engine
   !> @param[out] error  Set on allocation failure
   subroutine chebyshev_new_trafo(self, trafo, error)
      !> Grid instance
      class(moist_math_grid_1d_chebyshev_type), intent(in), target :: self
      !> Allocated transform engine
      class(moist_math_grid_1d_trafo_type), allocatable, intent(out) :: trafo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_1d_chebyshev_trafo_type), allocatable :: t
      integer :: stat

      allocate (t, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate Chebyshev radial trafo")
         return
      end if
      t%grid => self
      call move_alloc(t, trafo)
   end subroutine chebyshev_new_trafo

   !> Detach the trafo from its grid, idempotent
   !>
   !> @param[in,out] self  Trafo instance
   subroutine chebyshev_trafo_destroy(self)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self

      self%grid => null()
   end subroutine chebyshev_trafo_destroy

   !> Forward Fourier-Bessel transform r -> k via the dense matrix M
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  Function in r-space (size npts)
   !> @param[out]    f_out  Function in k-space (size nk)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_r2k(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> Function in r-space (size npts)
      real(wp), intent(in) :: f_in(:)
      !> Function in k-space (size nk)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(self%grid%fwd, f_in)
   end subroutine chebyshev_trafo_fbt_r2k

   !> Backward Fourier-Bessel transform k -> r via the dense matrix B
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  Function in k-space (size nk)
   !> @param[out]    f_out  Function in r-space (size npts)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_k2r(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> Function in k-space (size nk)
      real(wp), intent(in) :: f_in(:)
      !> Function in r-space (size npts)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(self%grid%bwd, f_in)
   end subroutine chebyshev_trafo_fbt_k2r

   !> Adjoint of the forward transform: y_r = M^T x_k (k -> r)
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  Input vector in k-space (size nk)
   !> @param[out]    f_out  Output vector in r-space (size npts)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_r2k_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> Input vector in k-space (size nk)
      real(wp), intent(in) :: f_in(:)
      !> Output vector in r-space (size npts)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(f_in, self%grid%fwd)
   end subroutine chebyshev_trafo_fbt_r2k_adj

   !> Adjoint of the backward transform: y_k = B^T x_r (r -> k)
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  Input vector in r-space (size npts)
   !> @param[out]    f_out  Output vector in k-space (size nk)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_k2r_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> Input vector in r-space (size npts)
      real(wp), intent(in) :: f_in(:)
      !> Output vector in k-space (size nk)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(f_in, self%grid%bwd)
   end subroutine chebyshev_trafo_fbt_k2r_adj

   !> Batched forward transform r -> k as a single GEMM: F_k = M F_r
   !>
   !> Each column of `f_in`/`f_out` is an independent field
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  r-space fields, shape (npts, nbatch)
   !> @param[out]    f_out  k-space fields, shape (nk, nbatch)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_r2k_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> r-space fields, shape (npts, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> k-space fields, shape (nk, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(self%grid%fwd, f_in)
   end subroutine chebyshev_trafo_fbt_r2k_all

   !> Batched backward transform k -> r as a single GEMM: F_r = B F_k
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  k-space fields, shape (nk, nbatch)
   !> @param[out]    f_out  r-space fields, shape (npts, nbatch)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_k2r_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> k-space fields, shape (nk, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> r-space fields, shape (npts, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(self%grid%bwd, f_in)
   end subroutine chebyshev_trafo_fbt_k2r_all

   !> Batched adjoint of the forward transform: F_r = M^T F_k (single GEMM)
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  k-space fields, shape (nk, nbatch)
   !> @param[out]    f_out  r-space fields, shape (npts, nbatch)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_r2k_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> k-space fields, shape (nk, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> r-space fields, shape (npts, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(transpose(self%grid%fwd), f_in)
   end subroutine chebyshev_trafo_fbt_r2k_adj_all

   !> Batched adjoint of the backward transform: F_k = B^T F_r (single GEMM)
   !>
   !> @param[in,out] self  Trafo instance
   !> @param[in]     f_in  r-space fields, shape (npts, nbatch)
   !> @param[out]    f_out  k-space fields, shape (nk, nbatch)
   !> @param[out]    error  Never set; a dense-matrix transform cannot fail
   subroutine chebyshev_trafo_fbt_k2r_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_chebyshev_trafo_type), intent(inout) :: self
      !> r-space fields, shape (npts, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> k-space fields, shape (nk, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      f_out = matmul(transpose(self%grid%bwd), f_in)
   end subroutine chebyshev_trafo_fbt_k2r_adj_all

end module moist_math_grid_1d_chebyshev
