!> Memory utilities for array reallocation
module moist_utils_mem
   use mctc_env, only: wp
   implicit none
   private

   public :: grow_array
   public :: filter_array

   !> Generic interface to grow and copy arrays of various types
   interface grow_array
      module procedure :: grow_array_real_1d
      module procedure :: grow_array_real_2d
      module procedure :: grow_array_int_1d
      module procedure :: grow_array_logical_1d
   end interface grow_array

   !> Generic interface to filter (compact) arrays in-place using a logical mask
   interface filter_array
      module procedure :: filter_array_real_1d
      module procedure :: filter_array_real_2d
      module procedure :: filter_array_int_1d
      module procedure :: filter_array_logical_1d
   end interface filter_array

contains

   !> Reallocate 1D real array to new capacity, preserving existing data
   !>
   !> @param[inout] array      Array to reallocate
   !> @param[in]    new_size   New array size (must be >= current size)
   !> @param[in]    fill_value Value for new elements (default 0.0)
   subroutine grow_array_real_1d(array, new_size, fill_value)
      !> Array to reallocate
      real(wp), allocatable, intent(inout) :: array(:)
      !> New array size
      integer, intent(in) :: new_size
      !> Fill value for new elements
      real(wp), intent(in), optional :: fill_value

      real(wp) :: fill
      real(wp), allocatable :: tmp(:)
      integer :: old_size

      old_size = 0
      if (allocated(array)) old_size = size(array)

      !> Early return if array already has the requested size
      if (old_size == new_size) return

      !> Error if attempting to shrink array
      if (old_size > new_size) then
         ! TODO: Proper error propagration
         error stop "grow_array_real_1d: Cannot shrink array"
      end if

      fill = 0.0_wp
      if (present(fill_value)) fill = fill_value

      allocate (tmp(new_size))
      if (old_size > 0) tmp(1:old_size) = array(1:old_size)
      tmp(old_size + 1:new_size) = fill
      call move_alloc(tmp, array)
   end subroutine grow_array_real_1d

   !> Reallocate 2D real array to new capacity, preserving existing data
   !>
   !> @param[inout] array      Array to reallocate
   !> @param[in]    dim1       First dimension size (must match existing if allocated)
   !> @param[in]    dim2       Second dimension size (must be >= current size)
   !> @param[in]    fill_value Value for new elements (default 0.0)
   subroutine grow_array_real_2d(array, dim1, dim2, fill_value)
      !> Array to reallocate
      real(wp), allocatable, intent(inout) :: array(:, :)
      !> First dimension size
      integer, intent(in) :: dim1
      !> Second dimension size
      integer, intent(in) :: dim2
      !> Fill value for new elements
      real(wp), intent(in), optional :: fill_value

      real(wp) :: fill
      real(wp), allocatable :: tmp(:, :)
      integer :: old_dim1, old_dim2

      old_dim1 = 0
      old_dim2 = 0
      if (allocated(array)) then
         old_dim1 = size(array, dim=1)
         old_dim2 = size(array, dim=2)
      else
         old_dim1 = dim1
      end if

      !> Enforce that first dimension remains constant
      if (old_dim1 > 0 .and. old_dim1 /= dim1) then
         ! TODO: Proper error propagration
         error stop "grow_array_real_2d: Cannot change first dimension"
      end if

      !> Early return if array already has the requested size
      if (old_dim2 == dim2) return

      !> Error if attempting to shrink array
      if (old_dim2 > dim2) then
         ! TODO: Proper error propagration
         error stop "grow_array_real_2d: Cannot shrink array"
      end if

      fill = 0.0_wp
      if (present(fill_value)) fill = fill_value

      allocate (tmp(dim1, dim2))
      if (old_dim2 > 0) tmp(:, 1:old_dim2) = array(:, 1:old_dim2)
      tmp(:, old_dim2 + 1:dim2) = fill
      call move_alloc(tmp, array)
   end subroutine grow_array_real_2d

   !> Reallocate 1D integer array to new capacity, preserving existing data
   !>
   !> @param[inout] array      Array to reallocate
   !> @param[in]    new_size   New array size (must be >= current size)
   !> @param[in]    fill_value Value for new elements (default 0)
   subroutine grow_array_int_1d(array, new_size, fill_value)
      !> Array to reallocate
      integer, allocatable, intent(inout) :: array(:)
      !> New array size
      integer, intent(in) :: new_size
      !> Fill value for new elements
      integer, intent(in), optional :: fill_value

      integer :: fill
      integer, allocatable :: tmp(:)
      integer :: old_size

      old_size = 0
      if (allocated(array)) old_size = size(array)

      !> Early return if array already has the requested size
      if (old_size == new_size) return

      !> Error if attempting to shrink array
      if (old_size > new_size) then
         ! TODO: Proper error propagration
         error stop "grow_array_int_1d: Cannot shrink array"
      end if

      fill = 0
      if (present(fill_value)) fill = fill_value

      allocate (tmp(new_size))
      if (old_size > 0) tmp(1:old_size) = array(1:old_size)
      tmp(old_size + 1:new_size) = fill
      call move_alloc(tmp, array)
   end subroutine grow_array_int_1d

   !> Reallocate 1D logical array to new capacity, preserving existing data
   !>
   !> @param[inout] array      Array to reallocate
   !> @param[in]    new_size   New array size (must be >= current size)
   !> @param[in]    fill_value Value for new elements (default .false.)
   subroutine grow_array_logical_1d(array, new_size, fill_value)
      !> Array to reallocate
      logical, allocatable, intent(inout) :: array(:)
      !> New array size
      integer, intent(in) :: new_size
      !> Fill value for new elements
      logical, intent(in), optional :: fill_value

      logical :: fill
      logical, allocatable :: tmp(:)
      integer :: old_size

      old_size = 0
      if (allocated(array)) old_size = size(array)

      !> Early return if array already has the requested size
      if (old_size == new_size) return

      !> Error if attempting to shrink array
      if (old_size > new_size) then
         ! TODO: Proper error propagration
         error stop "grow_array_logical_1d: Cannot shrink array"
      end if

      fill = .false.
      if (present(fill_value)) fill = fill_value

      allocate (tmp(new_size))
      if (old_size > 0) tmp(1:old_size) = array(1:old_size)
      tmp(old_size + 1:new_size) = fill
      call move_alloc(tmp, array)
   end subroutine grow_array_logical_1d

   !> Filter allocated 1D real array in-place.
   !>
   !> Compacts arr(1:n) to keep only elements where keep is .true.
   !> Safely skips unallocated arrays.
   !>
   !> @param[inout] arr    Allocatable 1D real array to filter
   !> @param[in]    n      Number of elements to consider
   !> @param[in]    keep   Logical mask (size n)
   !> @param[in]    nvalid Number of true values in keep
   subroutine filter_array_real_1d(arr, n, keep, nvalid)
      !> Array to filter
      real(wp), allocatable, intent(inout) :: arr(:)
      !> Number of elements to consider
      integer, intent(in) :: n, nvalid
      !> Logical mask
      logical, intent(in) :: keep(:)
      real(wp), allocatable :: tmp(:)
      integer :: i, j

      if (.not. allocated(arr)) return
      allocate (tmp(nvalid))
      j = 0
      do i = 1, n
         if (keep(i)) then
            j = j + 1
            tmp(j) = arr(i)
         end if
      end do
      call move_alloc(tmp, arr)
   end subroutine filter_array_real_1d

   !> Filter allocated 2D real array (dim1, n) in-place.
   !>
   !> Compacts arr(:, 1:n) to keep only columns where keep is .true.
   !> Processes each row independently. Safely skips unallocated arrays.
   !>
   !> @param[inout] arr    Allocatable 2D real array to filter
   !> @param[in]    n      Number of columns to consider
   !> @param[in]    keep   Logical mask (size n)
   !> @param[in]    nvalid Number of true values in keep
   subroutine filter_array_real_2d(arr, n, keep, nvalid)
      !> Array to filter
      real(wp), allocatable, intent(inout) :: arr(:, :)
      !> Number of columns to consider
      integer, intent(in) :: n, nvalid
      !> Logical mask
      logical, intent(in) :: keep(:)
      real(wp), allocatable :: tmp(:, :)
      integer :: i, j, dim1

      if (.not. allocated(arr)) return
      dim1 = size(arr, 1)
      allocate (tmp(dim1, nvalid))
      j = 0
      do i = 1, n
         if (keep(i)) then
            j = j + 1
            tmp(:, j) = arr(:, i)
         end if
      end do
      call move_alloc(tmp, arr)
   end subroutine filter_array_real_2d

   !> Filter allocated 1D integer array in-place.
   !>
   !> Compacts arr(1:n) to keep only elements where keep is .true.
   !> Safely skips unallocated arrays.
   !>
   !> @param[inout] arr    Allocatable 1D integer array to filter
   !> @param[in]    n      Number of elements to consider
   !> @param[in]    keep   Logical mask (size n)
   !> @param[in]    nvalid Number of true values in keep
   subroutine filter_array_int_1d(arr, n, keep, nvalid)
      !> Array to filter
      integer, allocatable, intent(inout) :: arr(:)
      !> Number of elements to consider
      integer, intent(in) :: n, nvalid
      !> Logical mask
      logical, intent(in) :: keep(:)
      integer, allocatable :: tmp(:)
      integer :: i, j

      if (.not. allocated(arr)) return
      allocate (tmp(nvalid))
      j = 0
      do i = 1, n
         if (keep(i)) then
            j = j + 1
            tmp(j) = arr(i)
         end if
      end do
      call move_alloc(tmp, arr)
   end subroutine filter_array_int_1d

   !> Filter allocated 1D logical array in-place.
   !>
   !> Compacts arr(1:n) to keep only elements where keep is .true.
   !> Safely skips unallocated arrays.
   !>
   !> @param[inout] arr    Allocatable 1D logical array to filter
   !> @param[in]    n      Number of elements to consider
   !> @param[in]    keep   Logical mask (size n)
   !> @param[in]    nvalid Number of true values in keep
   subroutine filter_array_logical_1d(arr, n, keep, nvalid)
      !> Array to filter
      logical, allocatable, intent(inout) :: arr(:)
      !> Number of elements to consider
      integer, intent(in) :: n, nvalid
      !> Logical mask
      logical, intent(in) :: keep(:)
      logical, allocatable :: tmp(:)
      integer :: i, j

      if (.not. allocated(arr)) return
      allocate (tmp(nvalid))
      j = 0
      do i = 1, n
         if (keep(i)) then
            j = j + 1
            tmp(j) = arr(i)
         end if
      end do
      call move_alloc(tmp, arr)
   end subroutine filter_array_logical_1d

end module moist_utils_mem
