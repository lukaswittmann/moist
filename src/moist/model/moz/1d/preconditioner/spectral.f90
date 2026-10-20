!> Radial spectral small-k filter for 1D MOZ Krylov vectors
!>
!> - Applies `w(k) = w0 + (1 - w0) k^2/(k^2 + k0^2)` per site pair
!> - Through the radial Fourier-Bessel transform
!> - Also smooths the JFNK Picard fallback step
!> - Adapted from moist_dev/pysofft (44874e7e) for the uniform DST-IV radial pair
!> - Filter copies the pair's nodes and transform
!> - Columns distributed over the model's OpenMP team, one transform clone per thread
module moist_model_moz_1d_preconditioner_spectral
   !$ use omp_lib, only: omp_get_thread_num
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: resolve_num_threads
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, transform_dst4
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, new_radial_trafo, check_trafo_shape
   use moist_model_moz_solver_core_preconditioner, only: moist_model_moz_solver_preconditioner_type
   implicit none(type, external)
   private
   public :: moist_model_moz_1d_preconditioner_spectral_type

   !> Transform scratch owned by one OpenMP thread
   type :: trafo_slot_type
      !> Independent scalar transform, without batched scratch
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
   end type trafo_slot_type

   !> Radial filter owning copies of the uniform radial pair and its transform
   type, extends(moist_model_moz_solver_preconditioner_type) :: moist_model_moz_1d_preconditioner_spectral_type
      !> OpenMP team size; zero takes the current runtime setting
      integer :: nthreads = 0
      !> Fourier-Bessel transform of the radial pair; absent until `setup`
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      !> Per-thread scalar transforms; retained between applications
      type(trafo_slot_type), allocatable, private :: slots(:)
      !> Reused k-space field, (npts, npair)
      real(wp), allocatable, private :: vk_work(:, :)
      !> r-space nodes (npts), bohr
      real(wp), allocatable :: r(:)
      !> r-space quadrature weights (npts), bohr
      real(wp), allocatable :: w(:)
      !> k-space nodes (npts), 1/bohr
      real(wp), allocatable :: k(:)
      !> Uniform r spacing, bohr
      real(wp) :: dr = 0.0_wp
      !> Filter w0 + (1 - w0) k^2/(k^2 + k0^2) at the k nodes (npts)
      real(wp), allocatable :: kfilter(:)
      !> Characteristic wavenumber, 1/bohr
      real(wp) :: k0 = 1.0_wp
      !> Weight at k = 0; zero suppresses, one leaves the vector unchanged
      real(wp) :: w0 = 0.1_wp
   contains
      !> Bind the radial pair and build the filter
      procedure :: setup => spectral_setup
      !> Whether the filter is built
      procedure :: ready => spectral_ready
      !> Filter a flat vector in place
      procedure :: apply => spectral_apply
      !> Status for the solver notes
      procedure :: describe => spectral_describe
      !> Release the filter and transform
      procedure :: release => spectral_release
      !> Transform columns r -> k and apply the filter
      procedure :: forward => spectral_forward
      !> Transform columns k -> r
      procedure :: backward => spectral_backward
      !> Grow the per-thread scalar transforms
      procedure, private :: ensure_slots => spectral_ensure_slots
      !> Transform columns using the requested team size
      procedure, private :: transform_columns => spectral_transform_columns
   end type moist_model_moz_1d_preconditioner_spectral_type

contains

   !> Copy a uniform DST-IV radial pair and precompute the filter
   !>
   !> Only the uniform pair from new_uniform_radial_pair is accepted: its
   !> transform is exactly invertible and its nodes ascend uniformly, as the
   !> k-space inverse assumes
   !>
   !> @param[inout] self  Filter; a previous setup is replaced, also on failure
   !> @param[in]    rgrid r-space radial grid
   !> @param[in]    kgrid k-space radial grid paired with rgrid
   !> @param[in]    k0    Characteristic wavenumber, 1/bohr; positive
   !> @param[in]    w0    Weight at k = 0, in [0, 1]
   !> @param[out]   error Non-uniform or mismatched pair, or invalid filter parameters
   subroutine spectral_setup(self, rgrid, kgrid, k0, w0, error)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      !> r-space and k-space radial grids
      type(moist_math_grid_radial_type), intent(in) :: rgrid, kgrid
      !> Filter parameters
      real(wp), intent(in) :: k0, w0
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%release()
      if (.not. (k0 > 0.0_wp .and. k0 <= huge(1.0_wp)) .or. .not. (w0 >= 0.0_wp .and. w0 <= 1.0_wp)) then
         call fatal_error(error, "Spectral filter needs k0 > 0 and w0 in [0, 1]")
         return
      end if
      if (.not. allocated(rgrid%r) .or. .not. allocated(rgrid%w) .or. .not. allocated(kgrid%r)) then
         call fatal_error(error, "Spectral filter needs built radial grids")
         return
      end if
      if (size(rgrid%r) /= size(kgrid%r) .or. size(rgrid%r) == 0) then
         call fatal_error(error, "Spectral filter needs paired radial grids of equal nonzero size")
         return
      end if
      if (rgrid%transform /= transform_dst4 .or. kgrid%transform /= transform_dst4 .or. &
         & .not. rgrid%spacing > 0.0_wp) then
         call fatal_error(error, "Spectral filter needs the uniform DST-IV radial pair")
         return
      end if
      call new_radial_trafo(self%trafo, rgrid, kgrid, error)
      if (allocated(error)) then
         if (allocated(self%trafo)) deallocate (self%trafo)
         return
      end if
      ! Clone before a serial call can allocate field-sized batched scratch
      allocate (self%slots(1))
      allocate (self%slots(1)%trafo, source=self%trafo)
      self%k0 = k0
      self%w0 = w0
      self%dr = rgrid%spacing
      self%r = rgrid%r
      self%w = rgrid%w
      self%k = kgrid%r
      self%kfilter = w0 + (1.0_wp - w0)*self%k**2/(self%k**2 + k0**2)
   end subroutine spectral_setup

   !> Whether the transform and filter are built
   !>
   !> @param[in]  self   filter
   !> @param[in]  rmse   residual RMS of the previous Newton step; unused
   logical function spectral_ready(self, rmse) result(usable)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(in) :: self
      !> Residual RMS of the previous Newton step
      real(wp), intent(in) :: rmse
      usable = allocated(self%trafo) .and. allocated(self%kfilter)
   end function spectral_ready

   !> Status for the solver notes: the filter, or that it is not built
   !>
   !> @param[in] self Filter
   !> @param[in] rmse Residual RMS of the previous Newton step; unused
   function spectral_describe(self, rmse) result(text)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(in) :: self
      !> Residual RMS of the previous Newton step
      real(wp), intent(in) :: rmse
      !> Status text
      character(len=:), allocatable :: text
      text = "spectral"
      if (.not. self%ready(rmse)) text = text//", not built"
   end function spectral_describe

   !> Transform the columns of an r-space field to k-space and apply the filter
   !>
   !> @param[inout] self  Filter with a built transform
   !> @param[in]    vr    r-space field (npts, npair)
   !> @param[out]   vk    Filtered k-space field (npts, npair)
   !> @param[out]   error Unbuilt filter or transform failure
   subroutine spectral_forward(self, vr, vk, error)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      !> r-space field (npts, npair)
      real(wp), intent(in) :: vr(:, :)
      !> Filtered k-space field (npts, npair)
      real(wp), intent(out) :: vk(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%ready(0.0_wp)) then
         call fatal_error(error, "Spectral filter is not built")
         return
      end if
      call self%transform_columns(.true., vr, vk, error)
   end subroutine spectral_forward

   !> Transform the columns of a k-space field back to r-space
   !>
   !> @param[inout] self  Filter with a built transform
   !> @param[in]    vk    k-space field (npts, npair)
   !> @param[out]   vr    r-space field (npts, npair)
   !> @param[out]   error Unbuilt filter or transform failure
   subroutine spectral_backward(self, vk, vr, error)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      !> k-space field (npts, npair)
      real(wp), intent(in) :: vk(:, :)
      !> r-space field (npts, npair)
      real(wp), intent(out) :: vr(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%ready(0.0_wp)) then
         call fatal_error(error, "Spectral filter is not built")
         return
      end if
      call self%transform_columns(.false., vk, vr, error)
   end subroutine spectral_backward

   !> Filter a flat vector in place
   !>
   !> - Per site pair: r -> k, multiply, k -> r
   !> - No-op until the filter is built
   !> - Vector size not a multiple of npts leaves the vector unchanged
   !> - Forward transform failure leaves the vector unchanged
   !> - One reused field-sized k-space buffer; the back transform writes into the vector
   !>
   !> @param[inout] self  Filter and its transform scratch
   !> @param[inout] vflat Flat (npts*npair) vector, filtered in place
   subroutine spectral_apply(self, vflat)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      !> Flat vector, modified in place
      real(wp), contiguous, target, intent(inout) :: vflat(:)
      !> Pair-major view of vflat (npts, npair)
      real(wp), pointer :: v(:, :)
      type(error_type), allocatable :: error
      integer :: npts, npair

      if (.not. self%ready(0.0_wp)) return
      npts = size(self%kfilter)
      if (mod(size(vflat), npts) /= 0) return
      npair = size(vflat)/npts
      if (npair == 0) return
      v(1:npts, 1:npair) => vflat
      if (allocated(self%vk_work)) then
         if (any(shape(self%vk_work) /= [npts, npair])) deallocate (self%vk_work)
      end if
      if (.not. allocated(self%vk_work)) allocate (self%vk_work(npts, npair))
      call self%forward(v, self%vk_work, error)
      if (allocated(error)) return
      call self%backward(self%vk_work, v, error)
   end subroutine spectral_apply

   !> Grow the scalar transform clones without copying batched scratch
   !>
   !> @param[inout] self Filter after setup
   !> @param[in]    nt   Number of transform slots required
   subroutine spectral_ensure_slots(self, nt)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      !> Number of slots required
      integer, intent(in) :: nt
      type(trafo_slot_type), allocatable :: slots(:)
      !> Separate copy source; ifx rejects SOURCE= from the array being allocated
      class(moist_math_grid_radial_trafo_type), allocatable :: template
      integer :: it, n

      n = size(self%slots)
      if (n >= nt) return
      allocate (slots(nt))
      do it = 1, n
         call move_alloc(self%slots(it)%trafo, slots(it)%trafo)
      end do
      if (nt > n) allocate (template, source=slots(1)%trafo)
      do it = n + 1, nt
         allocate (slots(it)%trafo, source=template)
      end do
      call move_alloc(slots, self%slots)
   end subroutine spectral_ensure_slots

   !> Transform site-pair columns and apply the k filter in the forward direction
   !>
   !> @param[inout] self  Filter after setup; owns independent transform scratch
   !> @param[in]    to_k  Forward transform and filtering when true; backward otherwise
   !> @param[in]    fin   Input field (npts, npair)
   !> @param[out]   fout  Output field (npts, npair)
   !> @param[out]   error Shape mismatch or transform failure
   subroutine spectral_transform_columns(self, to_k, fin, fout, error)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      !> Forward transform and filtering
      logical, intent(in) :: to_k
      !> Input field
      real(wp), intent(in) :: fin(:, :)
      !> Output field
      real(wp), intent(out) :: fout(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable :: terr
      integer :: nt, it, ip, trafo_flag

      call check_trafo_shape(shape(fin), shape(fout), size(self%kfilter), size(self%kfilter), "spectral", error)
      if (allocated(error)) return
      if (size(fin, 2) == 0) return
      nt = min(size(fin, 2), resolve_num_threads(self%nthreads))
      if (nt == 1) then
         if (to_k) then
            call self%trafo%fbt_r2k_all(fin, fout, error)
            if (allocated(error)) return
            do ip = 1, size(fout, 2)
               fout(:, ip) = fout(:, ip)*self%kfilter
            end do
         else
            call self%trafo%fbt_k2r_all(fin, fout, error)
         end if
         return
      end if

      call self%ensure_slots(nt)
      trafo_flag = 0
      !$omp parallel num_threads(nt) default(none) shared(self, to_k, fin, fout) &
      !$omp private(it, ip, terr) reduction(max:trafo_flag)
      it = 1
      !$ it = omp_get_thread_num() + 1
      !$omp do schedule(static)
      do ip = 1, size(fin, 2)
         if (to_k) then
            call self%slots(it)%trafo%fbt_r2k(fin(:, ip), fout(:, ip), terr)
         else
            call self%slots(it)%trafo%fbt_k2r(fin(:, ip), fout(:, ip), terr)
         end if
         if (allocated(terr)) then
            trafo_flag = 1
            deallocate (terr)
         else if (to_k) then
            fout(:, ip) = fout(:, ip)*self%kfilter
         end if
      end do
      !$omp end do
      !$omp end parallel
      if (trafo_flag /= 0) call fatal_error(error, "Radial transform failed in the spectral filter")
   end subroutine spectral_transform_columns

   !> Release the filter and transform, safe on repeated calls
   !>
   !> @param[inout] self Filter
   subroutine spectral_release(self)
      !> Filter
      class(moist_model_moz_1d_preconditioner_spectral_type), intent(inout) :: self
      if (allocated(self%trafo)) deallocate (self%trafo)
      if (allocated(self%slots)) deallocate (self%slots)
      if (allocated(self%vk_work)) deallocate (self%vk_work)
      if (allocated(self%r)) deallocate (self%r)
      if (allocated(self%w)) deallocate (self%w)
      if (allocated(self%k)) deallocate (self%k)
      if (allocated(self%kfilter)) deallocate (self%kfilter)
      self%dr = 0.0_wp
   end subroutine spectral_release

end module moist_model_moz_1d_preconditioner_spectral
