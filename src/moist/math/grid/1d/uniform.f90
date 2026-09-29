!> Uniform (equidistant) radial grid and its Fourier-Bessel transform engine
module moist_math_grid_1d_uniform
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_1d_base, only: moist_math_grid_1d_type, &
      & moist_math_grid_1d_trafo_type
   use moist_math_fft_dst4, only: dst4_plan_type, dst4_work_type
   implicit none(type, external)
   private

   public :: moist_math_grid_1d_uniform_type, new_uniform_radial_grid
   public :: moist_math_grid_1d_uniform_trafo_type

   !> Immutable uniform radial grid for 1D-RISM calculations
   !>
   !> Holds the equidistant nodes, precomputed reciprocals, and the cached
   !> DST-IV plan; contains no per-call mutable state, so it is safe to share
   !> read-only across threads
   type, extends(moist_math_grid_1d_type) :: moist_math_grid_1d_uniform_type
      !> r-space spacing (bohr)
      real(wp) :: dr = 0.0_wp
      !> k-space spacing (inverse bohr)
      real(wp) :: dk = 0.0_wp
      !> Precomputed pi = acos(-1.0_wp)
      real(wp) :: pi = 0.0_wp
      !> Precomputed 1/r
      real(wp), allocatable :: inv_r(:)
      !> Precomputed 1/k
      real(wp), allocatable :: inv_k(:)
      !> Shared single-column DST-IV plan, executed with per-thread arrays
      !>
      !> Forward and inverse Fourier-Bessel transforms differ only in the r/k
      !> weight diagonals applied around it, so one plan serves both
      type(dst4_plan_type) :: dst4
   contains
      procedure :: measure => uniform_measure
      procedure :: new_trafo => uniform_new_trafo
      procedure :: destroy => uniform_grid_destroy
   end type moist_math_grid_1d_uniform_type

   !> Per-thread Fourier-Bessel transform engine for a uniform radial grid
   !>
   !> Owns the work/output buffers and a pointer to its (immutable) grid
   type, extends(moist_math_grid_1d_trafo_type) :: moist_math_grid_1d_uniform_trafo_type
      !> Grid this trafo transforms on (not owned; must outlive the trafo)
      class(moist_math_grid_1d_uniform_type), pointer :: grid => null()
      !> Work array for the single-column transforms
      real(wp), allocatable :: work(:)
      !> Output array for the single-column transforms
      real(wp), allocatable :: tmp(:)
      !> Per-thread scratch for the grid's single-column DST-IV plan
      type(dst4_work_type) :: dst4_work
      !> Cached batched DST-IV plan (built lazily)
      type(dst4_plan_type) :: dst4_all
      !> Per-thread scratch for the batched DST-IV plan
      type(dst4_work_type) :: dst4_work_all
      !> Batch count the cached plan/buffers were built for (0 = none)
      integer :: nbatch = 0
      !> Batched work buffer, shape (npts, nbatch)
      real(wp), allocatable :: work_all(:, :)
      !> Batched output buffer, shape (npts, nbatch)
      real(wp), allocatable :: tmp_all(:, :)
   contains
      procedure :: fbt_r2k => uniform_trafo_fbt_r2k
      procedure :: fbt_k2r => uniform_trafo_fbt_k2r
      procedure :: fbt_r2k_adj => uniform_trafo_fbt_r2k_adj
      procedure :: fbt_k2r_adj => uniform_trafo_fbt_k2r_adj
      procedure :: fbt_r2k_all => uniform_trafo_fbt_r2k_all
      procedure :: fbt_k2r_all => uniform_trafo_fbt_k2r_all
      procedure :: fbt_r2k_adj_all => uniform_trafo_fbt_r2k_adj_all
      procedure :: fbt_k2r_adj_all => uniform_trafo_fbt_k2r_adj_all
      procedure :: destroy => uniform_trafo_destroy
   end type moist_math_grid_1d_uniform_trafo_type

contains

   !> Initialize a uniform radial grid for 1D-RISM
   !>
   !> Records the DST-IV geometry; the transform is later executed against the
   !> per-thread trafo buffers
   !>
   !> @param[out] grid  Initialized grid object
   !> @param[in]  npts  Number of grid points
   !> @param[in]  dr  r-space spacing in bohr, finite and positive
   !> @param[out] error  Set on invalid parameters or allocation failure
   subroutine new_uniform_radial_grid(grid, npts, dr, error)
      !> Initialized grid object
      type(moist_math_grid_1d_uniform_type), intent(out) :: grid
      !> Number of grid points
      integer, intent(in) :: npts
      !> r-space spacing in bohr
      real(wp), intent(in) :: dr
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i, stat

      if (npts < 1) then
         call fatal_error(error, "Uniform radial grid needs at least one point")
         return
      end if
      if (.not. ieee_is_finite(dr)) then
         call fatal_error(error, "Uniform radial grid needs a finite spacing")
         return
      end if
      if (dr <= 0.0_wp) then
         call fatal_error(error, "Uniform radial grid needs a positive spacing")
         return
      end if

      grid%npts = npts
      grid%dr = dr
      grid%pi = acos(-1.0_wp)
      grid%dk = grid%pi/(real(npts, wp)*dr)

      allocate (grid%r(npts), grid%k(npts), grid%inv_r(npts), grid%inv_k(npts), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate uniform radial grid nodes")
         return
      end if
      do i = 1, npts
         grid%r(i) = (real(i, wp) - 0.5_wp)*dr
         grid%k(i) = (real(i, wp) - 0.5_wp)*grid%dk
         grid%inv_r(i) = 1.0_wp/grid%r(i)
         grid%inv_k(i) = 1.0_wp/grid%k(i)
      end do

      call grid%dst4%init(npts, 1)
   end subroutine new_uniform_radial_grid

   !> 3D-radial integration weight of grid point i
   !>
   !> 4*pi*r_i^2*dr, the midpoint rule on the equidistant nodes
   !>
   !> @param[in] self  Grid instance
   !> @param[in] i     Grid point index
   pure function uniform_measure(self, i) result(w)
      !> Grid instance
      class(moist_math_grid_1d_uniform_type), intent(in) :: self
      !> Grid point index
      integer, intent(in) :: i
      !> Weight of point i
      real(wp) :: w

      w = 4.0_wp*self%pi*self%r(i)*self%r(i)*self%dr
   end function uniform_measure

   !> Deallocate grid arrays and destroy the DST-IV plan, idempotent
   !>
   !> @param[in,out] self  Grid instance
   subroutine uniform_grid_destroy(self)
      !> Grid instance
      class(moist_math_grid_1d_uniform_type), intent(inout) :: self

      call self%dst4%destroy()
      if (allocated(self%r)) deallocate (self%r)
      if (allocated(self%k)) deallocate (self%k)
      if (allocated(self%inv_r)) deallocate (self%inv_r)
      if (allocated(self%inv_k)) deallocate (self%inv_k)
      self%npts = 0
      self%dr = 0.0_wp
      self%dk = 0.0_wp
   end subroutine uniform_grid_destroy

   !> Allocate a `moist_math_grid_1d_uniform_trafo_type` bound to `self`
   !>
   !> Dispatched polymorphically via the base `new_trafo` binding
   !>
   !> @param[in]  self  Grid instance (must have the target attribute)
   !> @param[out] trafo  Allocated, bound transform engine
   !> @param[out] error  Set on allocation failure
   subroutine uniform_new_trafo(self, trafo, error)
      !> Grid instance
      class(moist_math_grid_1d_uniform_type), intent(in), target :: self
      !> Allocated transform engine
      class(moist_math_grid_1d_trafo_type), allocatable, intent(out) :: trafo
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_1d_uniform_trafo_type), allocatable :: t
      integer :: stat

      allocate (t, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate uniform radial trafo")
         return
      end if
      t%grid => self
      allocate (t%work(self%npts), t%tmp(self%npts), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to allocate uniform radial trafo buffers")
         return
      end if
      call self%dst4%new_work(t%dst4_work)
      call move_alloc(t, trafo)
   end subroutine uniform_new_trafo

   !> Release the trafo buffers and detach from the grid, idempotent
   !>
   !> @param[in,out] self  Trafo instance
   subroutine uniform_trafo_destroy(self)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self

      if (allocated(self%work)) deallocate (self%work)
      if (allocated(self%tmp)) deallocate (self%tmp)
      call self%dst4_work%destroy()
      call self%dst4_all%destroy()
      call self%dst4_work_all%destroy()
      if (allocated(self%work_all)) deallocate (self%work_all)
      if (allocated(self%tmp_all)) deallocate (self%tmp_all)
      self%nbatch = 0
      self%grid => null()
   end subroutine uniform_trafo_destroy

   !> Forward Fourier-Bessel transform r -> k (DST-IV)
   !>
   !> @param[in,out] self  Trafo instance (owns the scratch buffers)
   !> @param[in]     f_in  Function in r-space (size npts)
   !> @param[out]    f_out  Function in k-space (size npts)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_r2k(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> Function in r-space (size npts)
      real(wp), intent(in) :: f_in(:)
      !> Function in k-space (size npts)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac

      associate (g => self%grid)
         prefac = 2.0_wp*g%pi*g%dr
         self%work = f_in*g%r
         call g%dst4%execute(self%dst4_work, self%work, self%tmp, error)
         if (allocated(error)) return
         f_out = prefac*self%tmp*g%inv_k
      end associate
   end subroutine uniform_trafo_fbt_r2k

   !> Inverse Fourier-Bessel transform k -> r (DST-IV)
   !>
   !> @param[in,out] self  Trafo instance (owns the scratch buffers)
   !> @param[in]     f_in  Function in k-space (size npts)
   !> @param[out]    f_out  Function in r-space (size npts)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_k2r(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> Function in k-space (size npts)
      real(wp), intent(in) :: f_in(:)
      !> Function in r-space (size npts)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac

      associate (g => self%grid)
         prefac = g%dk/(4.0_wp*g%pi*g%pi)
         self%work = f_in*g%k
         call g%dst4%execute(self%dst4_work, self%work, self%tmp, error)
         if (allocated(error)) return
         f_out = prefac*self%tmp*g%inv_r
      end associate
   end subroutine uniform_trafo_fbt_k2r

   !> Adjoint (transpose) of the forward transform: y_r = (FBT_rk)^T x_k
   !>
   !> FBT_rk = diag(2*pi*dr*inv_k) . S . diag(r), with S = symmetric DST-IV
   !> (RODFT11); its transpose is diag(r) . S . diag(2*pi*dr*inv_k), i.e. the
   !> same DST call with the r- and k-weight diagonals swapped; the
   !> unnormalized S appears once in the forward and once here, so
   !> <FBT_rk a, b> = <a, y_r> holds to machine precision
   !>
   !> @param[in,out] self  Trafo instance (owns the scratch buffers)
   !> @param[in]     f_in  Input vector in k-space (size npts)
   !> @param[out]    f_out  Output vector in r-space (size npts)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_r2k_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> Input vector in k-space (size npts)
      real(wp), intent(in) :: f_in(:)
      !> Output vector in r-space (size npts)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac

      associate (g => self%grid)
         prefac = 2.0_wp*g%pi*g%dr
         self%work = prefac*f_in*g%inv_k
         call g%dst4%execute(self%dst4_work, self%work, self%tmp, error)
         if (allocated(error)) return
         f_out = g%r*self%tmp
      end associate
   end subroutine uniform_trafo_fbt_r2k_adj

   !> Adjoint (transpose) of the inverse transform: y_k = (FBT_kr)^T x_r
   !>
   !> FBT_kr = diag(dk/(4*pi^2)*inv_r) . S . diag(k); its transpose is
   !> diag(k) . S . diag(dk/(4*pi^2)*inv_r): the same DST call with the r/k
   !> weight diagonals swapped
   !>
   !> @param[in,out] self  Trafo instance (owns the scratch buffers)
   !> @param[in]     f_in  Input vector in r-space (size npts)
   !> @param[out]    f_out  Output vector in k-space (size npts)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_k2r_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> Input vector in r-space (size npts)
      real(wp), intent(in) :: f_in(:)
      !> Output vector in k-space (size npts)
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac

      associate (g => self%grid)
         prefac = g%dk/(4.0_wp*g%pi*g%pi)
         self%work = prefac*f_in*g%inv_r
         call g%dst4%execute(self%dst4_work, self%work, self%tmp, error)
         if (allocated(error)) return
         f_out = g%k*self%tmp
      end associate
   end subroutine uniform_trafo_fbt_k2r_adj

   !> Ensure the cached batched DST-IV plan and buffers fit `nbatch` columns
   !>
   !> Rebuilds them when the batch count changes; one plan serves all four
   !> directions (DST-IV is its own inverse), so the r/k weight diagonals are
   !> applied separately by the callers
   !>
   !> The backend is plan-free, so the plan is only a geometry record and this
   !> cannot fail; only the later `execute` call against ducc0 can, which is
   !> why this helper itself takes no error channel
   !>
   !> @param[in,out] self  Trafo instance (owns the cached plan/buffers)
   !> @param[in]     nbatch  Number of columns to transform together
   subroutine uniform_trafo_ensure_all(self, nbatch)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> Number of columns to transform together
      integer, intent(in) :: nbatch

      integer :: npts

      if (self%nbatch == nbatch) return

      npts = self%grid%npts
      call self%dst4_all%destroy()
      call self%dst4_work_all%destroy()
      if (allocated(self%work_all)) deallocate (self%work_all)
      if (allocated(self%tmp_all)) deallocate (self%tmp_all)
      allocate (self%work_all(npts, nbatch), self%tmp_all(npts, nbatch))

      call self%dst4_all%init(npts, nbatch)
      call self%dst4_all%new_work(self%dst4_work_all)
      self%nbatch = nbatch
   end subroutine uniform_trafo_ensure_all

   !> Batched forward Fourier-Bessel transform r -> k (DST-IV) over all columns
   !>
   !> An empty batch (`nbatch == 0`) returns immediately: `uniform_trafo_ensure_all`
   !> would otherwise treat a first call with `nbatch == 0` as already matching
   !> its initial cached width of 0 and skip allocating `work_all`/`tmp_all`,
   !> which are then passed unallocated to the non-allocatable `execute` dummies
   !>
   !> @param[in,out] self  Trafo instance (owns the cached plan/buffers)
   !> @param[in]     f_in  r-space fields, shape (npts, nbatch)
   !> @param[out]    f_out  k-space fields, shape (npts, nbatch)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_r2k_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> r-space fields, shape (npts, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> k-space fields, shape (npts, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac
      integer :: j, nb

      nb = size(f_in, 2)
      if (nb == 0) return
      call uniform_trafo_ensure_all(self, nb)
      associate (g => self%grid)
         prefac = 2.0_wp*g%pi*g%dr
         do j = 1, nb
            self%work_all(:, j) = f_in(:, j)*g%r
         end do
         call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
         if (allocated(error)) return
         do j = 1, nb
            f_out(:, j) = prefac*self%tmp_all(:, j)*g%inv_k
         end do
      end associate
   end subroutine uniform_trafo_fbt_r2k_all

   !> Batched inverse Fourier-Bessel transform k -> r (DST-IV) over all columns
   !>
   !> An empty batch returns immediately; see `uniform_trafo_fbt_r2k_all`
   !>
   !> @param[in,out] self  Trafo instance (owns the cached plan/buffers)
   !> @param[in]     f_in  k-space fields, shape (npts, nbatch)
   !> @param[out]    f_out  r-space fields, shape (npts, nbatch)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_k2r_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> k-space fields, shape (npts, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> r-space fields, shape (npts, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac
      integer :: j, nb

      nb = size(f_in, 2)
      if (nb == 0) return
      call uniform_trafo_ensure_all(self, nb)
      associate (g => self%grid)
         prefac = g%dk/(4.0_wp*g%pi*g%pi)
         do j = 1, nb
            self%work_all(:, j) = f_in(:, j)*g%k
         end do
         call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
         if (allocated(error)) return
         do j = 1, nb
            f_out(:, j) = prefac*self%tmp_all(:, j)*g%inv_r
         end do
      end associate
   end subroutine uniform_trafo_fbt_k2r_all

   !> Batched adjoint of the forward transform: y_r = (FBT_rk)^T x_k over all columns
   !>
   !> The r/k weight diagonals are swapped versus the forward call; an empty
   !> batch returns immediately; see `uniform_trafo_fbt_r2k_all`
   !>
   !> @param[in,out] self  Trafo instance (owns the cached plan/buffers)
   !> @param[in]     f_in  k-space vectors, shape (npts, nbatch)
   !> @param[out]    f_out  r-space vectors, shape (npts, nbatch)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_r2k_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> k-space vectors, shape (npts, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> r-space vectors, shape (npts, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac
      integer :: j, nb

      nb = size(f_in, 2)
      if (nb == 0) return
      call uniform_trafo_ensure_all(self, nb)
      associate (g => self%grid)
         prefac = 2.0_wp*g%pi*g%dr
         do j = 1, nb
            self%work_all(:, j) = prefac*f_in(:, j)*g%inv_k
         end do
         call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
         if (allocated(error)) return
         do j = 1, nb
            f_out(:, j) = g%r*self%tmp_all(:, j)
         end do
      end associate
   end subroutine uniform_trafo_fbt_r2k_adj_all

   !> Batched adjoint of the inverse transform: y_k = (FBT_kr)^T x_r over all columns
   !>
   !> The r/k weight diagonals are swapped versus the backward call; an empty
   !> batch returns immediately; see `uniform_trafo_fbt_r2k_all`
   !>
   !> @param[in,out] self  Trafo instance (owns the cached plan/buffers)
   !> @param[in]     f_in  r-space vectors, shape (npts, nbatch)
   !> @param[out]    f_out  k-space vectors, shape (npts, nbatch)
   !> @param[out]    error  Set if the DST-IV backend failed
   subroutine uniform_trafo_fbt_k2r_adj_all(self, f_in, f_out, error)
      !> Trafo instance
      class(moist_math_grid_1d_uniform_trafo_type), intent(inout) :: self
      !> r-space vectors, shape (npts, nbatch)
      real(wp), intent(in) :: f_in(:, :)
      !> k-space vectors, shape (npts, nbatch)
      real(wp), intent(out) :: f_out(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: prefac
      integer :: j, nb

      nb = size(f_in, 2)
      if (nb == 0) return
      call uniform_trafo_ensure_all(self, nb)
      associate (g => self%grid)
         prefac = g%dk/(4.0_wp*g%pi*g%pi)
         do j = 1, nb
            self%work_all(:, j) = prefac*f_in(:, j)*g%inv_r
         end do
         call self%dst4_all%execute(self%dst4_work_all, self%work_all, self%tmp_all, error)
         if (allocated(error)) return
         do j = 1, nb
            f_out(:, j) = g%k*self%tmp_all(:, j)
         end do
      end associate
   end subroutine uniform_trafo_fbt_k2r_adj_all

end module moist_math_grid_1d_uniform
