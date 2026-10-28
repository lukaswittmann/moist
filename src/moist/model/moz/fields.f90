!> Threaded field helpers shared by the MOZ drivers
!>
!> - Fields are (npts, ncol) arrays; solver vectors hold them flat, column after column
!> - Iterate is feasible when c and t are finite and 1 + c + t >= g_floor
module moist_model_moz_fields
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: g_floor, state_feasible, scale_flat, subtract_flat, unflatten, residual_difference, copy_field
   public :: field_num_threads

   !> Lower bound on 1 + h; allows hard-core underflow and roundoff near zero
   real(wp), parameter :: g_floor = -1.0e-12_wp
   !> Use the same small-field cutoff as the closure and Krylov vector kernels
   integer, parameter :: parallel_field_size = 65536

contains

   !> Keep small radial fields on the same serial team as their Krylov vector work
   !>
   !> @param[in] nfield  number of field values
   !> @param[in] nt      resolved OpenMP team size
   pure integer function field_num_threads(nfield, nt) result(nthreads)
      !> Number of field values
      integer, intent(in) :: nfield
      !> Resolved team size
      integer, intent(in) :: nt
      nthreads = merge(nt, 1, nfield >= parallel_field_size)
   end function field_num_threads

   !> Finite c and t with 1 + h = 1 + c + t above the floor
   !>
   !> @param[in] c   direct correlation, (npts, ncol)
   !> @param[in] t   indirect correlation of the same shape
   !> @param[in] nt  OpenMP team size
   logical function state_feasible(c, t, nt) result(ok)
      !> Direct and indirect correlation
      real(wp), intent(in) :: c(:, :), t(:, :)
      !> Team size
      integer, intent(in) :: nt
      integer :: i, col
      ok = .true.
      if (size(c) < parallel_field_size) then
         ok = all(ieee_is_finite(c) .and. ieee_is_finite(t) .and. 1.0_wp + c + t >= g_floor)
         return
      end if
      !$omp parallel do collapse(2) schedule(static) num_threads(nt) reduction(.and.:ok)
      do col = 1, size(c, 2)
         do i = 1, size(c, 1)
            ok = ok .and. ieee_is_finite(c(i, col)) .and. ieee_is_finite(t(i, col)) .and. &
               & 1.0_wp + c(i, col) + t(i, col) >= g_floor
         end do
      end do
      !$omp end parallel do
   end function state_feasible

   !> c' v on the field layout: dc(i, col) = dcdt(i, col) v(i + (col - 1) npts)
   !>
   !> @param[in]  dcdt  closure derivative (npts, ncol)
   !> @param[in]  v     flat direction (npts ncol)
   !> @param[out] dc    product (npts, ncol)
   !> @param[in]  nt    OpenMP team size
   subroutine scale_flat(dcdt, v, dc, nt)
      !> Closure derivative
      real(wp), intent(in) :: dcdt(:, :)
      !> Flat direction
      real(wp), intent(in) :: v(:)
      !> Product
      real(wp), intent(out) :: dc(:, :)
      !> Team size
      integer, intent(in) :: nt
      integer :: i, col, npts
      npts = size(dc, 1)
      if (size(dc) < parallel_field_size) then
         do col = 1, size(dc, 2)
            dc(:, col) = dcdt(:, col)*v(1 + (col - 1)*npts:col*npts)
         end do
         return
      end if
      !$omp parallel do collapse(2) schedule(static) num_threads(nt)
      do col = 1, size(dc, 2)
         do i = 1, npts
            dc(i, col) = dcdt(i, col)*v(i + (col - 1)*npts)
         end do
      end do
      !$omp end parallel do
   end subroutine scale_flat

   !> Flat difference out = v - field or v - dcdt*field, with finite flag
   !>
   !> @param[in]  v       flat vector (npts ncol)
   !> @param[in]  field   field (npts, ncol)
   !> @param[out] out     flat difference (npts ncol)
   !> @param[in]  nt      OpenMP team size
   !> @param[out] finite  whether every value is finite
   !> @param[in]  dcdt    optional scaling of field (npts, ncol)
   subroutine subtract_flat(v, field, out, nt, finite, dcdt)
      !> Flat vector
      real(wp), intent(in) :: v(:)
      !> Field
      real(wp), intent(in) :: field(:, :)
      !> Flat difference
      real(wp), intent(out) :: out(:)
      !> Team size
      integer, intent(in) :: nt
      !> Finite flag
      logical, intent(out) :: finite
      !> Optional scaling
      real(wp), intent(in), optional :: dcdt(:, :)
      integer :: i, col, npts, k
      npts = size(field, 1)
      finite = .true.
      if (size(field) < parallel_field_size) then
         do col = 1, size(field, 2)
            k = (col - 1)*npts
            if (present(dcdt)) then
               out(k + 1:k + npts) = v(k + 1:k + npts) - dcdt(:, col)*field(:, col)
            else
               out(k + 1:k + npts) = v(k + 1:k + npts) - field(:, col)
            end if
            finite = finite .and. all(ieee_is_finite(out(k + 1:k + npts)))
         end do
         return
      end if
      if (present(dcdt)) then
         !$omp parallel do collapse(2) schedule(static) num_threads(nt) private(k) reduction(.and.:finite)
         do col = 1, size(field, 2)
            do i = 1, npts
               k = i + (col - 1)*npts
               out(k) = v(k) - dcdt(i, col)*field(i, col)
               finite = finite .and. ieee_is_finite(out(k))
            end do
         end do
         !$omp end parallel do
      else
         !$omp parallel do collapse(2) schedule(static) num_threads(nt) private(k) reduction(.and.:finite)
         do col = 1, size(field, 2)
            do i = 1, npts
               k = i + (col - 1)*npts
               out(k) = v(k) - field(i, col)
               finite = finite .and. ieee_is_finite(out(k))
            end do
         end do
         !$omp end parallel do
      end if
   end subroutine subtract_flat

   !> Field from a flat vector
   !>
   !> @param[in]  v      flat vector (npts ncol)
   !> @param[out] field  field (npts, ncol)
   !> @param[in]  nt     OpenMP team size
   subroutine unflatten(v, field, nt)
      !> Flat vector
      real(wp), intent(in) :: v(:)
      !> Field
      real(wp), intent(out) :: field(:, :)
      !> Team size
      integer, intent(in) :: nt
      integer :: i, col, npts
      npts = size(field, 1)
      if (size(field) < parallel_field_size) then
         do col = 1, size(field, 2)
            field(:, col) = v(1 + (col - 1)*npts:col*npts)
         end do
         return
      end if
      !$omp parallel do collapse(2) schedule(static) num_threads(nt)
      do col = 1, size(field, 2)
         do i = 1, npts
            field(i, col) = v(i + (col - 1)*npts)
         end do
      end do
      !$omp end parallel do
   end subroutine unflatten

   !> In place r = gamma - r; whether it is finite
   !>
   !> @param[in]    gamma   iterate (npts, ncol)
   !> @param[inout] r       map on input, residual on output (npts, ncol)
   !> @param[in]    nt      OpenMP team size
   !> @param[out]   finite  whether every value is finite
   subroutine residual_difference(gamma, r, nt, finite)
      !> Iterate
      real(wp), intent(in) :: gamma(:, :)
      !> Map, then residual
      real(wp), intent(inout) :: r(:, :)
      !> Team size
      integer, intent(in) :: nt
      !> Finite flag
      logical, intent(out) :: finite
      integer :: i, col
      finite = .true.
      if (size(r) < parallel_field_size) then
         r = gamma - r
         finite = all(ieee_is_finite(r))
         return
      end if
      !$omp parallel do collapse(2) schedule(static) num_threads(nt) reduction(.and.:finite)
      do col = 1, size(r, 2)
         do i = 1, size(r, 1)
            r(i, col) = gamma(i, col) - r(i, col)
            finite = finite .and. ieee_is_finite(r(i, col))
         end do
      end do
      !$omp end parallel do
   end subroutine residual_difference

   !> Copy into an existing field, sharing large staging passes across the team
   !>
   !> - Caller checks matching shapes
   !> - No allocation or reshape temporary
   !>
   !> @param[in]  source  field (npts, ncol)
   !> @param[out] field   preallocated destination of the same shape
   !> @param[in]  nt      OpenMP team size
   subroutine copy_field(source, field, nt)
      !> Source field
      real(wp), intent(in) :: source(:, :)
      !> Preallocated destination of the same shape
      real(wp), intent(out) :: field(:, :)
      !> Team size
      integer, intent(in) :: nt
      integer :: i, col

      if (nt == 1 .or. size(field) < parallel_field_size) then
         field = source
         return
      end if
      !$omp parallel do collapse(2) schedule(static) num_threads(nt)
      do col = 1, size(field, 2)
         do i = 1, size(field, 1)
            field(i, col) = source(i, col)
         end do
      end do
      !$omp end parallel do
   end subroutine copy_field

end module moist_model_moz_fields
