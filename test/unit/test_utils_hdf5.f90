module test_utils_hdf5
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, &
      & skip_test
   use mctc_env_error, only: moist_error_type => error_type
   use moist_utils_hdf5io, only: hdf5_file
   implicit none(type, external)
   private
   public :: collect_utils_hdf5

contains

   subroutine collect_utils_hdf5(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("CreateFile", test_create_file), &
         & new_unittest("WriteReadScalar", test_write_read_scalar), &
         & new_unittest("WriteReadArray1D", test_write_read_array1d), &
         & new_unittest("WriteReadArray2D", test_write_read_array2d), &
         & new_unittest("WriteReadArray3D", test_write_read_array3d), &
         & new_unittest("WriteReadString", test_write_read_string), &
         & new_unittest("GroupOperations", test_group_operations), &
         & new_unittest("IntegerArrays", test_integer_arrays), &
         & new_unittest("Attributes", test_attributes), &
         & new_unittest("Delete", test_delete), &
         & new_unittest("OpenUnopenedHandle", test_open_unopened_handle, &
         &    should_fail=.true.), &
         & new_unittest("NestedOpenGroupWithoutClose", &
         &    test_nested_open_group_without_close, should_fail=.true.), &
         & new_unittest("DeleteMissingPath", test_delete_missing_path, &
         &    should_fail=.true.), &
         & new_unittest("DeleteNonEmptyGroup", &
         &    test_delete_non_empty_group), &
         & new_unittest("DatasetOverwriteSamePath", &
         &    test_dataset_overwrite_same_path, should_fail=.true.), &
         & new_unittest("TypeMismatchRead", test_type_mismatch_read, &
         &    should_fail=.true.), &
         & new_unittest("RankMismatchRead", test_rank_mismatch_read, &
         &    should_fail=.true.), &
         & new_unittest("ReadOnlyBlocksMutation", &
         &    test_read_only_blocks_mutation, should_fail=.true.), &
         & new_unittest("AddaMissingPath", test_adda_missing_path, should_fail=.true.), &
         & new_unittest("GroupOpenClose", test_group_open_close), &
         & new_unittest("CloseOpenGroup", test_close_with_open_group) &
         & ]
   end subroutine collect_utils_hdf5

   subroutine cleanup_test_file(filename)
      character(*), intent(in) :: filename

      integer :: ierr
      integer :: unit
      logical :: exists

      inquire (file=filename, exist=exists)
      if (.not. exists) return

      open (newunit=unit, file=filename, status="old", iostat=ierr)
      if (ierr /= 0) return
      close (unit, status="delete")
   end subroutine cleanup_test_file

   subroutine test_create_file(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      logical :: exists
      type(moist_error_type), allocatable :: h5err

      call h5f%open("test_create.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      inquire (file="test_create.h5", exist=exists)
      call check(error, exists, more="HDF5 file was not created")
      if (allocated(error)) return

      open (unit=99, file="test_create.h5", status="old")
      close (unit=99, status="delete")

      inquire (file="test_create.h5", exist=exists)
      call check(error,.not. exists, more="HDF5 file was not deleted")
      if (allocated(error)) return
   end subroutine test_create_file

   subroutine test_write_read_scalar(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp) :: value_in, value_out
      integer :: int_in, int_out
      type(moist_error_type), allocatable :: h5err

      value_in = 3.14159265358979_wp
      int_in = 42

      call h5f%open("test_scalar.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/real_scalar", value_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/int_scalar", int_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_scalar.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/real_scalar", value_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/int_scalar", int_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, abs(value_out - value_in) < 1.0e-12_wp, &
         & more="Real scalar read/write mismatch")
      if (allocated(error)) return

      call check(error, int_out == int_in, &
         & more="Integer scalar read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_scalar.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_write_read_scalar

   subroutine test_write_read_array1d(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp), allocatable :: arr_in(:), arr_out(:)
      integer :: i
      type(moist_error_type), allocatable :: h5err

      allocate (arr_in(10))
      do i = 1, 10
         arr_in(i) = real(i, wp)*0.1_wp
      end do

      call h5f%open("test_array1d.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/array_1d", arr_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_array1d.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/array_1d", arr_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, size(arr_out) == size(arr_in), &
         & more="Array 1D size mismatch")
      if (allocated(error)) return

      call check(error, all(abs(arr_out - arr_in) < 1.0e-12_wp), &
         & more="Array 1D read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_array1d.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_write_read_array1d

   subroutine test_write_read_array2d(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp), allocatable :: arr_in(:, :), arr_out(:, :)
      integer :: i, j
      type(moist_error_type), allocatable :: h5err

      allocate (arr_in(3, 4))
      do j = 1, 4
         do i = 1, 3
            arr_in(i, j) = real(i + j*10, wp)
         end do
      end do

      call h5f%open("test_array2d.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/array_2d", arr_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_array2d.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/array_2d", arr_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, size(arr_out, 1) == size(arr_in, 1) .and. &
         & size(arr_out, 2) == size(arr_in, 2), &
         & more="Array 2D shape mismatch")
      if (allocated(error)) return

      call check(error, all(abs(arr_out - arr_in) < 1.0e-12_wp), &
         & more="Array 2D read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_array2d.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_write_read_array2d

   subroutine test_write_read_string(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      character(:), allocatable :: str_in, str_out
      type(moist_error_type), allocatable :: h5err

      str_in = "Hello, HDF5!"

      call h5f%open("test_string.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/string_data", str_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_string.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/string_data", str_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, str_out == str_in, &
         & more="String read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_string.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_write_read_string

   subroutine test_group_operations(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp) :: value_in, value_out
      logical :: exists
      type(moist_error_type), allocatable :: h5err

      value_in = 42.0696969696_wp

      call h5f%open("test_group.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/group1/subgroup/data", value_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_group.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%exist("/group1", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, exists, more="Group /group1 should exist")
      if (allocated(error)) return

      call h5f%exist("/group1/subgroup", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, exists, more="Group /group1/subgroup should exist")
      if (allocated(error)) return

      call h5f%exist("/group1/subgroup/data", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, exists, more="Dataset /group1/subgroup/data should exist")
      if (allocated(error)) return

      call h5f%get("/group1/subgroup/data", value_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, abs(value_out - value_in) < 1.0e-12_wp, &
         & more="Nested data read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_group.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_group_operations

   subroutine test_write_read_array3d(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp), allocatable :: real_in(:, :, :), real_out(:, :, :)
      integer, allocatable :: int_in(:, :, :), int_out(:, :, :)
      integer :: i, j, k
      type(moist_error_type), allocatable :: h5err

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

      call h5f%open("test_array3d.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/real_3d", real_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/int_3d", int_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_array3d.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/real_3d", real_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/int_3d", int_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, size(real_out, 1) == size(real_in, 1) .and. &
         & size(real_out, 2) == size(real_in, 2) .and. &
         & size(real_out, 3) == size(real_in, 3), &
         & more="Real 3D array shape mismatch")
      if (allocated(error)) return

      call check(error, all(abs(real_out - real_in) < 1.0e-12_wp), &
         & more="Real 3D array read/write mismatch")
      if (allocated(error)) return

      call check(error, size(int_out, 1) == size(int_in, 1) .and. &
         & size(int_out, 2) == size(int_in, 2) .and. &
         & size(int_out, 3) == size(int_in, 3), &
         & more="Integer 3D array shape mismatch")
      if (allocated(error)) return

      call check(error, all(int_out == int_in), &
         & more="Integer 3D array read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_array3d.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_write_read_array3d

   subroutine test_integer_arrays(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      integer, allocatable :: arr1d_in(:), arr1d_out(:)
      integer, allocatable :: arr2d_in(:, :), arr2d_out(:, :)
      integer, allocatable :: arr3d_in(:, :, :), arr3d_out(:, :, :)
      integer :: i, j, k
      type(moist_error_type), allocatable :: h5err

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

      call h5f%open("test_int_arrays.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/int_1d", arr1d_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/int_2d", arr2d_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/int_3d", arr3d_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_int_arrays.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/int_1d", arr1d_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/int_2d", arr2d_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/int_3d", arr3d_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, size(arr1d_out) == size(arr1d_in), &
         & more="Integer 1D array size mismatch")
      if (allocated(error)) return
      call check(error, all(arr1d_out == arr1d_in), &
         & more="Integer 1D array read/write mismatch")
      if (allocated(error)) return

      call check(error, size(arr2d_out, 1) == size(arr2d_in, 1) .and. &
         & size(arr2d_out, 2) == size(arr2d_in, 2), &
         & more="Integer 2D array shape mismatch")
      if (allocated(error)) return
      call check(error, all(arr2d_out == arr2d_in), &
         & more="Integer 2D array read/write mismatch")
      if (allocated(error)) return

      call check(error, size(arr3d_out, 1) == size(arr3d_in, 1) .and. &
         & size(arr3d_out, 2) == size(arr3d_in, 2) .and. &
         & size(arr3d_out, 3) == size(arr3d_in, 3), &
         & more="Integer 3D array shape mismatch")
      if (allocated(error)) return
      call check(error, all(arr3d_out == arr3d_in), &
         & more="Integer 3D array read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_int_arrays.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_integer_arrays

   subroutine test_attributes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      character(:), allocatable :: attr_value
      type(moist_error_type), allocatable :: h5err

      call h5f%open("test_attrs.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/dataset", 1.0_wp, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%adda("/dataset", "units", "angstrom", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%adda("/dataset", "description", "Test attribute", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_attrs.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%geta("/dataset", "units", attr_value, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, attr_value == "angstrom", &
         & more="String attribute 'units' mismatch")
      if (allocated(error)) return

      call h5f%geta("/dataset", "description", attr_value, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, attr_value == "Test attribute", &
         & more="String attribute 'description' mismatch")
      if (allocated(error)) return

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      open (unit=99, file="test_attrs.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_attributes

   subroutine test_adda_missing_path(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      logical :: adda_failed
      type(moist_error_type), allocatable :: h5err

      call h5f%open("test_adda_missing_path.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/dataset", 1.0_wp, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%adda("/dataset/missing", "units", "angstrom", error=h5err)
      adda_failed = allocated(h5err)

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      open (unit=99, file="test_adda_missing_path.h5", status="old")
      close (unit=99, status="delete")

      call check(error,.not. adda_failed, &
         & more="adda should fail when the target path does not exist")
      if (allocated(error)) return

   end subroutine test_adda_missing_path

   subroutine test_delete(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp) :: value_in
      logical :: exists
      type(moist_error_type), allocatable :: h5err

      value_in = 123.456_wp

      call h5f%open("test_delete.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/to_delete", value_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/to_keep", value_in*2.0_wp, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_delete.h5", status="old", action="rw", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%exist("/to_delete", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, exists, more="Dataset /to_delete should exist before delete")
      if (allocated(error)) return

      call h5f%delete("/to_delete", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%exist("/to_delete", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error,.not. exists, more="Dataset /to_delete should not exist after delete")
      if (allocated(error)) return

      call h5f%exist("/to_keep", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call check(error, exists, more="Dataset /to_keep should still exist")
      if (allocated(error)) return

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      open (unit=99, file="test_delete.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_delete

   subroutine test_open_unopened_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      type(moist_error_type), allocatable :: h5err
      logical :: failed

      call h5f%add("/dataset", 1.0_wp, error=h5err)
      failed = allocated(h5err)

      call check(error,.not. failed, &
         & more="Adding to an unopened HDF5 handle should fail")
      if (allocated(error)) return
   end subroutine test_open_unopened_handle

   subroutine test_nested_open_group_without_close(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      type(moist_error_type), allocatable :: h5err
      logical :: nested_failed
      character(*), parameter :: filename = "test_nested_open_group.h5"

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/parent/child", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%open_group("/parent", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%open_group("child", error=h5err)
      nested_failed = allocated(h5err)

      call h5f%close_group(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. nested_failed, &
         & more="Opening a second group without closing the first should fail")
      if (allocated(error)) return
   end subroutine test_nested_open_group_without_close

   subroutine test_delete_missing_path(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      type(moist_error_type), allocatable :: h5err
      logical :: delete_failed
      character(*), parameter :: filename = "test_delete_missing_path.h5"

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/present", 1.0_wp, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%delete("/missing", error=h5err)
      delete_failed = allocated(h5err)

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. delete_failed, &
         & more="Deleting a missing path should fail")
      if (allocated(error)) return
   end subroutine test_delete_missing_path

   subroutine test_delete_non_empty_group(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      type(moist_error_type), allocatable :: h5err
      logical :: exists
      character(*), parameter :: filename = "test_delete_non_empty_group.h5"

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/group/data", 1.0_wp, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%open(filename, status="old", action="rw", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%delete("/group", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%exist("/group", exists, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. exists, &
         & more="Deleting a non-empty group should remove the group")
      if (allocated(error)) return
   end subroutine test_delete_non_empty_group

   subroutine test_dataset_overwrite_same_path(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp) :: value_in
      type(moist_error_type), allocatable :: h5err
      logical :: overwrite_failed
      character(*), parameter :: filename = "test_overwrite_same_path.h5"

      value_in = 2.5_wp

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/dataset", value_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/dataset", value_in*2.0_wp, error=h5err)
      overwrite_failed = allocated(h5err)

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. overwrite_failed, &
         & more="Creating a dataset at an existing path should fail")
      if (allocated(error)) return
   end subroutine test_dataset_overwrite_same_path

   subroutine test_type_mismatch_read(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      character(:), allocatable :: str_in
      integer :: int_out
      type(moist_error_type), allocatable :: h5err
      logical :: mismatch_failed
      character(*), parameter :: filename = "test_type_mismatch_read.h5"

      str_in = "mismatch"

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/string_data", str_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%open(filename, status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%get("/string_data", int_out, error=h5err)
      mismatch_failed = allocated(h5err)

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. mismatch_failed, &
         & more="Reading a string dataset as an integer should fail")
      if (allocated(error)) return
   end subroutine test_type_mismatch_read

   subroutine test_rank_mismatch_read(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      integer, allocatable :: arr2d(:, :)
      integer :: int_out
      integer :: i, j
      type(moist_error_type), allocatable :: h5err
      logical :: mismatch_failed
      character(*), parameter :: filename = "test_rank_mismatch_read.h5"

      allocate (arr2d(2, 3))
      do j = 1, 3
         do i = 1, 2
            arr2d(i, j) = i + j*10
         end do
      end do

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/matrix", arr2d, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%open(filename, status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%get("/matrix", int_out, error=h5err)
      mismatch_failed = allocated(h5err)

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. mismatch_failed, &
         & more="Reading a rank-2 dataset as scalar should fail")
      if (allocated(error)) return
   end subroutine test_rank_mismatch_read

   subroutine test_read_only_blocks_mutation(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      type(moist_error_type), allocatable :: h5err
      logical :: mutation_failed
      character(*), parameter :: filename = "test_read_only_blocks_mutation.h5"

      call h5f%open(filename, status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/baseline", 1.0_wp, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call h5f%close(error=h5err)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call h5f%open(filename, status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if
      call h5f%add("/blocked", 2.0_wp, error=h5err)
      mutation_failed = allocated(h5err)

      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         call cleanup_test_file(filename)
         return
      end if

      call cleanup_test_file(filename)

      call check(error,.not. mutation_failed, &
         & more="Read-only mode should block dataset creation")
      if (allocated(error)) return
   end subroutine test_read_only_blocks_mutation

   subroutine test_group_open_close(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp) :: value_in, value_out
      integer :: int_in, int_out
      type(moist_error_type), allocatable :: h5err

      value_in = 9.876_wp
      int_in = 99

      call h5f%open("test_group_oc.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/mygroup", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%open_group("/mygroup", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("data_real", value_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("data_int", int_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close_group(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_group_oc.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/mygroup/data_real", value_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/mygroup/data_int", int_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call check(error, abs(value_out - value_in) < 1.0e-12_wp, &
         & more="Group data real read/write mismatch")
      if (allocated(error)) return

      call check(error, int_out == int_in, &
         & more="Group data int read/write mismatch")
      if (allocated(error)) return

      open (unit=99, file="test_group_oc.h5", status="old")
      close (unit=99, status="delete")
   end subroutine test_group_open_close

   subroutine test_close_with_open_group(error)
      type(error_type), allocatable, intent(out) :: error
      type(hdf5_file) :: h5f
      real(wp) :: value_in, value_out
      type(moist_error_type), allocatable :: h5err

      value_in = 9.123_wp

      call h5f%open("test_close_with_open_group.h5", status="new", action="w", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("/group", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%open_group("/group", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%add("value", value_in, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      call h5f%open("test_close_with_open_group.h5", status="old", action="r", error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%get("/group/value", value_out, error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if
      call h5f%close(error=h5err)
      if (allocated(h5err)) then
         call test_failed(error, h5err%message)
         return
      end if

      open (unit=99, file="test_close_with_open_group.h5", status="old")
      close (unit=99, status="delete")

      call check(error, abs(value_out - value_in) < 1.0e-12_wp, &
         & more="File content should remain readable after closing with an open group")
      if (allocated(error)) return

   end subroutine test_close_with_open_group

end module test_utils_hdf5
