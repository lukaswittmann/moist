!> Binary container round trips, groups, attributes and refused operations
module test_utils_binio
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use mctc_env, only: wp, moist_error => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_utils_binio, only: binary_file
   implicit none(type, external)
   private

   public :: collect_utils_binio

   !> Tolerance of real round trips
   real(wp), parameter :: thr = 1.0e-12_wp

contains

   !> Collect all binary container tests
   subroutine collect_utils_binio(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("create_file", test_create_file), &
         new_unittest("write_read_scalar", test_write_read_scalar), &
         new_unittest("write_read_array1d", test_write_read_array1d), &
         new_unittest("write_read_array2d", test_write_read_array2d), &
         new_unittest("write_read_array3d", test_write_read_array3d), &
         new_unittest("write_read_string", test_write_read_string), &
         new_unittest("group_operations", test_group_operations), &
         new_unittest("integer_arrays", test_integer_arrays), &
         new_unittest("attributes", test_attributes), &
         new_unittest("delete", test_delete), &
         new_unittest("open_unopened_handle", test_open_unopened_handle), &
         new_unittest("nested_open_group_without_close", test_nested_open_group_without_close), &
         new_unittest("delete_missing_path", test_delete_missing_path), &
         new_unittest("delete_non_empty_group", test_delete_non_empty_group), &
         new_unittest("dataset_overwrite_same_path", test_dataset_overwrite_same_path), &
         new_unittest("type_mismatch_read", test_type_mismatch_read), &
         new_unittest("rank_mismatch_read", test_rank_mismatch_read), &
         new_unittest("read_only_blocks_mutation", test_read_only_blocks_mutation), &
         new_unittest("adda_missing_path", test_adda_missing_path), &
         new_unittest("group_open_close", test_group_open_close), &
         new_unittest("close_open_group", test_close_open_group), &
         new_unittest("get_shape", test_get_shape), &
         new_unittest("malformed_index", test_malformed_index) &
         ]

   end subroutine collect_utils_binio

   !> Delete a test file if it exists
   !>
   !> @param[in] filename Path of the file
   subroutine remove_file(filename)
      !> Path of the file
      character(*), intent(in) :: filename
      !> Unit used for the deletion
      integer :: unit
      !> I/O status code
      integer :: stat
      !> Whether the file exists
      logical :: exists

      inquire (file=filename, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=filename, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete", iostat=stat)
   end subroutine remove_file

   !> Fail on a library error, then close the container and remove its file
   !>
   !> @param[inout] error    Test failure
   !> @param[in]    ioerr    Library error, possibly unallocated
   !> @param[inout] binf     Container under test
   !> @param[in]    filename Path of the container
   subroutine require_ok(error, ioerr, binf, filename)
      !> Test failure
      type(error_type), allocatable, intent(inout) :: error
      !> Library error, possibly unallocated
      type(moist_error), allocatable, intent(in) :: ioerr
      !> Container under test
      type(binary_file), intent(inout) :: binf
      !> Path of the container
      character(*), intent(in) :: filename

      if (.not. allocated(ioerr)) return
      call test_failed(error, ioerr%message)
      call binf%close()
      call remove_file(filename)
   end subroutine require_ok

   !> Fail unless `ioerr` is set and names `expected`
   !>
   !> @param[inout] error    Test failure
   !> @param[in]    ioerr    Library error, possibly unallocated
   !> @param[in]    expected Substring the message must contain
   subroutine require_error(error, ioerr, expected)
      !> Test failure
      type(error_type), allocatable, intent(inout) :: error
      !> Library error, possibly unallocated
      type(moist_error), allocatable, intent(in) :: ioerr
      !> Substring the message must contain
      character(*), intent(in) :: expected

      if (.not. allocated(ioerr)) then
         call test_failed(error, "expected an error containing '"//expected//"'")
      else if (index(ioerr%message, expected) == 0) then
         call test_failed(error, "message '"//ioerr%message//"' lacks '"//expected//"'")
      end if
   end subroutine require_error

   !> Opening a new container creates a file that opens again
   subroutine test_create_file(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-create.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      logical :: opened, exists, removed

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      opened = binf%is_open
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      inquire (file=filename, exist=exists)

      ! An empty container still carries a committed index
      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call remove_file(filename)
      inquire (file=filename, exist=removed)

      call check(error, opened, "handle not open after creating the file")
      if (allocated(error)) return
      call check(error, .not. binf%is_open, "handle still open after closing")
      if (allocated(error)) return
      call check(error, exists, "binary file was not created")
      if (allocated(error)) return
      call check(error, .not. removed, "binary file was not deleted")
   end subroutine test_create_file

   !> Real and integer scalars survive a round trip
   subroutine test_write_read_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-scalar.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp) :: value_in, value_out
      integer :: int_in, int_out

      value_in = 3.14159265358979_wp
      int_in = 42

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/real_scalar", value_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/int_scalar", int_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/real_scalar", value_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/int_scalar", int_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, abs(value_out - value_in) < thr, "real scalar read/write mismatch")
      if (allocated(error)) return
      call check(error, int_out == int_in, "integer scalar read/write mismatch")
   end subroutine test_write_read_scalar

   !> A real vector survives a round trip
   subroutine test_write_read_array1d(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-array1d.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp), allocatable :: arr_in(:), arr_out(:)
      integer :: i

      allocate (arr_in(10))
      do i = 1, 10
         arr_in(i) = real(i, wp)*0.1_wp
      end do

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/array_1d", arr_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/array_1d", arr_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, size(arr_out) == size(arr_in), "array 1D size mismatch")
      if (allocated(error)) return
      call check(error, all(abs(arr_out - arr_in) < thr), "array 1D read/write mismatch")
   end subroutine test_write_read_array1d

   !> A real matrix survives a round trip with its shape
   subroutine test_write_read_array2d(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-array2d.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp), allocatable :: arr_in(:, :), arr_out(:, :)
      integer :: i, j

      allocate (arr_in(3, 4))
      do j = 1, 4
         do i = 1, 3
            arr_in(i, j) = real(i + j*10, wp)
         end do
      end do

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/array_2d", arr_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/array_2d", arr_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, all(shape(arr_out) == shape(arr_in)), "array 2D shape mismatch")
      if (allocated(error)) return
      call check(error, all(abs(arr_out - arr_in) < thr), "array 2D read/write mismatch")
   end subroutine test_write_read_array2d

   !> Real and integer rank-3 arrays survive a round trip with their shapes
   subroutine test_write_read_array3d(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-array3d.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp), allocatable :: real_in(:, :, :), real_out(:, :, :)
      integer, allocatable :: int_in(:, :, :), int_out(:, :, :)
      integer :: i, j, k

      allocate (real_in(2, 3, 4))
      do k = 1, 4
         do j = 1, 3
            do i = 1, 2
               real_in(i, j, k) = real(i + j*10 + k*100, wp)
            end do
         end do
      end do

      allocate (int_in(2, 2, 2))
      do k = 1, 2
         do j = 1, 2
            do i = 1, 2
               int_in(i, j, k) = i + j*10 + k*100
            end do
         end do
      end do

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/real_3d", real_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/int_3d", int_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/real_3d", real_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/int_3d", int_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, all(shape(real_out) == shape(real_in)), "real 3D array shape mismatch")
      if (allocated(error)) return
      call check(error, all(abs(real_out - real_in) < thr), "real 3D array read/write mismatch")
      if (allocated(error)) return
      call check(error, all(shape(int_out) == shape(int_in)), "integer 3D array shape mismatch")
      if (allocated(error)) return
      call check(error, all(int_out == int_in), "integer 3D array read/write mismatch")
   end subroutine test_write_read_array3d

   !> A string survives a round trip with its length
   subroutine test_write_read_string(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-string.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      character(:), allocatable :: str_in, str_out

      str_in = "Hello, moist!"

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/string_data", str_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/string_data", str_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, len(str_out) == len(str_in), "string length mismatch")
      if (allocated(error)) return
      call check(error, str_out == str_in, "string read/write mismatch")
   end subroutine test_write_read_string

   !> Writing a nested dataset creates its parent groups
   subroutine test_group_operations(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-group.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp) :: value_in, value_out
      logical :: group_exists, subgroup_exists, data_exists

      value_in = 42.0696969696_wp

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/group1/subgroup/data", value_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/group1", group_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/group1/subgroup", subgroup_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/group1/subgroup/data", data_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/group1/subgroup/data", value_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, group_exists, "group /group1 should exist")
      if (allocated(error)) return
      call check(error, subgroup_exists, "group /group1/subgroup should exist")
      if (allocated(error)) return
      call check(error, data_exists, "dataset /group1/subgroup/data should exist")
      if (allocated(error)) return
      call check(error, abs(value_out - value_in) < thr, "nested data read/write mismatch")
   end subroutine test_group_operations

   !> Integer arrays of rank 1 to 3 survive a round trip with their shapes
   subroutine test_integer_arrays(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-int-arrays.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      integer, allocatable :: arr1d_in(:), arr1d_out(:)
      integer, allocatable :: arr2d_in(:, :), arr2d_out(:, :)
      integer, allocatable :: arr3d_in(:, :, :), arr3d_out(:, :, :)
      integer :: i, j, k

      allocate (arr1d_in(5))
      do i = 1, 5
         arr1d_in(i) = i*10
      end do

      allocate (arr2d_in(3, 4))
      do j = 1, 4
         do i = 1, 3
            arr2d_in(i, j) = i + j*10
         end do
      end do

      allocate (arr3d_in(2, 2, 3))
      do k = 1, 3
         do j = 1, 2
            do i = 1, 2
               arr3d_in(i, j, k) = i + j*10 + k*100
            end do
         end do
      end do

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/int_1d", arr1d_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/int_2d", arr2d_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/int_3d", arr3d_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/int_1d", arr1d_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/int_2d", arr2d_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/int_3d", arr3d_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, size(arr1d_out) == size(arr1d_in), "integer 1D array size mismatch")
      if (allocated(error)) return
      call check(error, all(arr1d_out == arr1d_in), "integer 1D array read/write mismatch")
      if (allocated(error)) return
      call check(error, all(shape(arr2d_out) == shape(arr2d_in)), "integer 2D array shape mismatch")
      if (allocated(error)) return
      call check(error, all(arr2d_out == arr2d_in), "integer 2D array read/write mismatch")
      if (allocated(error)) return
      call check(error, all(shape(arr3d_out) == shape(arr3d_in)), "integer 3D array shape mismatch")
      if (allocated(error)) return
      call check(error, all(arr3d_out == arr3d_in), "integer 3D array read/write mismatch")
   end subroutine test_integer_arrays

   !> String attributes survive a round trip and a repeated name replaces the value
   subroutine test_attributes(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-attrs.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      character(:), allocatable :: units, description

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/dataset", 1.0_wp, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%adda("/dataset", "units", "bohr", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%adda("/dataset", "description", "Test attribute", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%adda("/dataset", "units", "angstrom", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%geta("/dataset", "units", units, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%geta("/dataset", "description", description, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, units == "angstrom", "string attribute 'units' mismatch")
      if (allocated(error)) return
      call check(error, description == "Test attribute", "string attribute 'description' mismatch")
   end subroutine test_attributes

   !> Deleting a dataset keeps its siblings and persists after reopening
   subroutine test_delete(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-delete.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp) :: value_in, value_out
      logical :: before, after, kept, reopened

      value_in = 123.456_wp

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/to_delete", value_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/to_keep", value_in*2.0_wp, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="rw", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/to_delete", before, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%delete("/to_delete", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/to_delete", after, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/to_keep", kept, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      ! The updated index replaces the original one on disk
      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/to_delete", reopened, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/to_keep", value_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, before, "dataset /to_delete should exist before delete")
      if (allocated(error)) return
      call check(error, .not. after, "dataset /to_delete should not exist after delete")
      if (allocated(error)) return
      call check(error, kept, "dataset /to_keep should still exist")
      if (allocated(error)) return
      call check(error, .not. reopened, "deletion of /to_delete was not committed")
      if (allocated(error)) return
      call check(error, abs(value_out - 2.0_wp*value_in) < thr, "dataset /to_keep changed by the update")
   end subroutine test_delete

   !> Adding to a handle that was never opened is refused
   subroutine test_open_unopened_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr

      call binf%add("/dataset", 1.0_wp, error=ioerr)
      call require_error(error, ioerr, "closed binary file")
      if (allocated(error)) return
      call check(error, .not. binf%is_open, "refused call opened the handle")
   end subroutine test_open_unopened_handle

   !> Opening a second group before closing the first is refused
   subroutine test_nested_open_group_without_close(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-nested-group.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, nested_err
      logical :: still_open, child_exists

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/parent/child", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%open_group("/parent", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open_group("child", error=nested_err)
      ! The first group stays open and relative names still resolve inside it
      still_open = binf%group_open
      call binf%exist("child", child_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%close_group(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, nested_err, "nested group")
      if (allocated(error)) return
      call check(error, still_open, "refused nested open closed the current group")
      if (allocated(error)) return
      call check(error, child_exists, "refused nested open moved the current group")
   end subroutine test_nested_open_group_without_close

   !> Deleting a missing path is refused and leaves the index intact
   subroutine test_delete_missing_path(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-delete-missing.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, delete_err
      logical :: present_kept

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/present", 1.0_wp, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%delete("/missing", error=delete_err)
      call binf%exist("/present", present_kept, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, delete_err, "Failed to delete")
      if (allocated(error)) return
      call check(error, present_kept, "refused delete removed another dataset")
   end subroutine test_delete_missing_path

   !> Deleting a group removes everything below it
   subroutine test_delete_non_empty_group(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-delete-group.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      logical :: group_exists, data_exists

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/group/data", 1.0_wp, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="rw", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%delete("/group", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/group", group_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/group/data", data_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, .not. group_exists, "deleting a non-empty group should remove the group")
      if (allocated(error)) return
      call check(error, .not. data_exists, "deleting a group should remove its datasets")
   end subroutine test_delete_non_empty_group

   !> Creating a dataset at an existing path is refused and keeps the original
   subroutine test_dataset_overwrite_same_path(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-overwrite.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, overwrite_err
      real(wp) :: value_in, value_out

      value_in = 2.5_wp

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/dataset", value_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/dataset", value_in*2.0_wp, error=overwrite_err)
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/dataset", value_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, overwrite_err, "path exists")
      if (allocated(error)) return
      call check(error, abs(value_out - value_in) < thr, "refused overwrite changed the dataset")
   end subroutine test_dataset_overwrite_same_path

   !> Reading a string dataset as an integer is refused
   subroutine test_type_mismatch_read(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-type-mismatch.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, mismatch_err
      integer :: int_out

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/string_data", "mismatch", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      int_out = -1
      call binf%get("/string_data", int_out, error=mismatch_err)
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, mismatch_err, "type mismatch")
      if (allocated(error)) return
      call check(error, int_out, 0)
   end subroutine test_type_mismatch_read

   !> Reading a rank-2 dataset as a scalar is refused
   subroutine test_rank_mismatch_read(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-rank-mismatch.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, mismatch_err
      integer, allocatable :: arr2d(:, :)
      integer :: int_out
      integer :: i, j

      allocate (arr2d(2, 3))
      do j = 1, 3
         do i = 1, 2
            arr2d(i, j) = i + j*10
         end do
      end do

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/matrix", arr2d, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      int_out = -1
      call binf%get("/matrix", int_out, error=mismatch_err)
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, mismatch_err, "rank mismatch")
      if (allocated(error)) return
      call check(error, int_out, 0)
   end subroutine test_rank_mismatch_read

   !> A container opened for reading refuses new datasets
   subroutine test_read_only_blocks_mutation(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-read-only.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, mutation_err
      logical :: blocked_exists, baseline_exists

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/baseline", 1.0_wp, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/blocked", 2.0_wp, error=mutation_err)
      call binf%exist("/blocked", blocked_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%exist("/baseline", baseline_exists, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, mutation_err, "read-only")
      if (allocated(error)) return
      call check(error, .not. blocked_exists, "refused dataset entered the index")
      if (allocated(error)) return
      call check(error, baseline_exists, "refused dataset removed another dataset")
   end subroutine test_read_only_blocks_mutation

   !> Adding an attribute to a missing path is refused
   subroutine test_adda_missing_path(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-adda-missing.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, adda_err

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/dataset", 1.0_wp, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%adda("/dataset/missing", "units", "angstrom", error=adda_err)
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call require_error(error, adda_err, "missing path")
   end subroutine test_adda_missing_path

   !> Relative names written inside an open group land below it
   subroutine test_group_open_close(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-group-open-close.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp) :: value_in, value_out
      integer :: int_in, int_out

      value_in = 9.876_wp
      int_in = 99

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/mygroup", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%open_group("/mygroup", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("data_real", value_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("data_int", int_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close_group(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/mygroup/data_real", value_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/mygroup/data_int", int_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, abs(value_out - value_in) < thr, "group data real read/write mismatch")
      if (allocated(error)) return
      call check(error, int_out == int_in, "group data int read/write mismatch")
   end subroutine test_group_open_close

   !> Closing the file while a group is open still commits its contents
   subroutine test_close_open_group(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-close-open-group.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      real(wp) :: value_in, value_out

      value_in = 9.123_wp

      call binf%open(filename, status="new", action="w", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("/group", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%open_group("/group", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%add("value", value_in, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get("/group/value", value_out, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, abs(value_out - value_in) < thr, &
                 "file content should remain readable after closing with an open group")
   end subroutine test_close_open_group

   !> Extents of datasets come from the index; groups have none
   subroutine test_get_shape(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-shape.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr, group_err
      real(wp) :: matrix(3, 4)
      integer, allocatable :: matrix_dims(:), scalar_dims(:), string_dims(:), group_dims(:)

      matrix = 1.0_wp
      call binf%open(filename, status="new", action="w", error=ioerr)
      if (.not. allocated(ioerr)) call binf%add("/data/matrix", matrix, error=ioerr)
      if (.not. allocated(ioerr)) call binf%add("/data/count", 7, error=ioerr)
      if (.not. allocated(ioerr)) call binf%add("/data/name", "water", error=ioerr)
      if (.not. allocated(ioerr)) call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return

      call binf%open(filename, status="old", action="r", error=ioerr)
      if (.not. allocated(ioerr)) call binf%get_shape("/data/matrix", matrix_dims, error=ioerr)
      if (.not. allocated(ioerr)) call binf%get_shape("/data/count", scalar_dims, error=ioerr)
      if (.not. allocated(ioerr)) call binf%get_shape("/data/name", string_dims, error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call binf%get_shape("/data", group_dims, error=group_err)
      call binf%close(error=ioerr)
      call require_ok(error, ioerr, binf, filename)
      if (allocated(error)) return
      call remove_file(filename)

      call check(error, size(matrix_dims) == 2, "matrix rank not reported")
      if (allocated(error)) return
      call check(error, all(matrix_dims == [3, 4]), "matrix extents not reported")
      if (allocated(error)) return
      call check(error, size(scalar_dims) == 0 .and. size(string_dims) == 0, "scalars reported extents")
      if (allocated(error)) return
      call require_error(error, group_err, "Not a dataset")
   end subroutine test_get_shape

   !> Index records the file cannot hold are refused on opening
   !>
   !> - Huge extents, a huge record count and a payload past the index each fail
   !>   with an error instead of an allocation
   subroutine test_malformed_index(error)
      type(error_type), allocatable, intent(out) :: error
      character(*), parameter :: filename = "moist-test-binio-malformed.mbin"
      type(binary_file) :: binf
      type(moist_error), allocatable :: ioerr
      integer(int64) :: index_offset
      integer :: variant, unit

      do variant = 1, 3
         call binf%open(filename, status="new", action="w", error=ioerr)
         if (.not. allocated(ioerr)) call binf%add("/x", [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp], error=ioerr)
         if (.not. allocated(ioerr)) call binf%close(error=ioerr)
         call require_ok(error, ioerr, binf, filename)
         if (allocated(error)) return

         ! Header: index offset at byte 21; index: count, then the record of /x with
         ! type and rank (8 bytes), dims (24 bytes), offset and size
         open (newunit=unit, file=filename, status="old", access="stream", form="unformatted", &
               action="readwrite")
         read (unit, pos=21) index_offset
         select case (variant)
         case (1)
            write (unit, pos=index_offset + 12_int64) 2_int64**40
         case (2)
            write (unit, pos=index_offset) huge(1_int32)
         case (3)
            write (unit, pos=index_offset + 36_int64) index_offset
         case default
            call require_error(error, ioerr, "invalid test variant")
            if (allocated(error)) return
         end select
         close (unit)

         call binf%open(filename, status="old", action="r", error=ioerr)
         call binf%close()
         call remove_file(filename)
         call require_error(error, ioerr, "Corrupt index")
         if (allocated(error)) return
      end do
   end subroutine test_malformed_index

end module test_utils_binio
