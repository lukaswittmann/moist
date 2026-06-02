!> HDF5 file I/O utilities for reading and writing scientific data.
!>
!> This module provides a high-level object-oriented interface to HDF5,
!> supporting scalar and multidimensional arrays of integers, reals, and strings.
!> The [[hdf5_file]] type encapsulates file handles and provides generic
!> bindings for type-transparent read/write operations.
!>
!> The HDF5 Fortran interface is initialized on first use and never explicitly
!> finalized. This is safe because the operating system reclaims all resources
!> on process exit, and calling h5close_f while other library users exist can
!> cause crashes. The module is thread-safe for initialization when OpenMP is
!> enabled.
!>
!> When built without HDF5 support, all operations return an error indicating
!> that HDF5 is not available.
module moist_utils_hdf5io
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
#ifdef WITH_HDF5
   use HDF5
   use H5LT
#endif
   implicit none(type, external)

   public :: hdf5_file

   private

   !> Flag indicating whether HDF5 Fortran interface has been initialized.
   !> This is never reset to false - the OS handles cleanup on process exit.
   logical, save :: hdf5_initialized = .false.

   !> HDF5 file handle wrapper for reading and writing datasets.
   !>
   !> Provides methods to open/close files, manage groups, and read/write
   !> datasets of various types (integers, reals, strings) and ranks (0-3D).
   type :: hdf5_file

      !> Path to the HDF5 file
      character(:), allocatable :: filename
      !> Flag indicating whether the file is currently open
      logical :: is_open = .false.
      !> Flag indicating whether a group is currently open
      logical :: group_open = .false.
#ifdef WITH_HDF5
      !> Current location identifier (file or group)
      integer(HID_T), private :: lid = -1
      !> Group identifier when a group is open
      integer(HID_T), private :: gid = -1
      !> Saved location identifier before entering a group
      integer(HID_T), private :: glid = -1
#else
      !> Dummy location identifier (HDF5 not available)
      integer, private :: lid = -1
      !> Dummy group identifier (HDF5 not available)
      integer, private :: gid = -1
      !> Dummy saved location identifier (HDF5 not available)
      integer, private :: glid = -1
#endif

   contains
      !> Open an HDF5 file for reading or writing
      procedure :: open => hdf_open_file
      !> Close the HDF5 file
      procedure :: close => hdf_close_file

      !> Open a group within the file
      procedure :: open_group => hdf_open_group
      !> Close the current group
      procedure :: close_group => hdf_close_group

      !> Check if a dataset or group exists
      procedure :: exist => hdf_exist

      !> Delete a dataset or group
      procedure :: delete => hdf_delete

      !> Add a group, dataset, or attribute
      generic :: add => hdf_add_group, &
         hdf_add_int, &
         hdf_add_int1d, &
         hdf_add_int2d, &
         hdf_add_int3d, &
         hdf_add_real, &
         hdf_add_real1d, &
         hdf_add_real2d, &
         hdf_add_real3d, &
         hdf_add_string

      !> Retrieve a dataset value
      generic :: get => hdf_get_int, &
         hdf_get_int1d, &
         hdf_get_int2d, &
         hdf_get_int3d, &
         hdf_get_real, &
         hdf_get_real1d, &
         hdf_get_real2d, &
         hdf_get_real3d, &
         hdf_get_string

      !> Add a string attribute to a path
      generic :: adda => hdf_adda_string

      !> Get a string attribute from a path
      generic :: geta => hdf_geta_string

      !> Get the current HDF5 location identifier (file or group)
      procedure :: get_location_id => hdf_get_location_id

      procedure, private :: hdf_add_group
      procedure, private :: hdf_ensure_parent_groups
      procedure, private :: hdf_add_int
      procedure, private :: hdf_get_int
      procedure, private :: hdf_add_int1d
      procedure, private :: hdf_get_int1d
      procedure, private :: hdf_add_int2d
      procedure, private :: hdf_get_int2d
      procedure, private :: hdf_add_int3d
      procedure, private :: hdf_get_int3d
      procedure, private :: hdf_add_real
      procedure, private :: hdf_get_real
      procedure, private :: hdf_add_real1d
      procedure, private :: hdf_get_real1d
      procedure, private :: hdf_add_real2d
      procedure, private :: hdf_get_real2d
      procedure, private :: hdf_add_real3d
      procedure, private :: hdf_get_real3d
      procedure, private :: hdf_get_string
      procedure, private :: hdf_add_string
      procedure, private :: hdf_adda_string
      procedure, private :: hdf_geta_string
   end type hdf5_file

contains

   !> Get the current HDF5 location identifier (file or group).
   !> @param[in] self HDF5 file instance
   !> @return Location identifier, or -1 if HDF5 not available
   function hdf_get_location_id(self) result(loc_id)
      class(hdf5_file), intent(in) :: self
      integer :: loc_id

      loc_id = self%lid

   end function hdf_get_location_id

   !> Ensure the HDF5 Fortran interface is initialized.
   !> Thread-safe initialization using OpenMP critical section.
   !> @param[out] error Error handler
   subroutine hdf5_ensure_initialized(error)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr
      !> Flag indicating initialization failed
      logical :: init_failed

      init_failed = .false.

      !$omp critical(hdf5_init)
      if (.not. hdf5_initialized) then
         call h5open_f(ierr)
         if (ierr /= 0) then
            init_failed = .true.
         else
            hdf5_initialized = .true.
         end if
      end if
      !$omp end critical(hdf5_init)

      if (init_failed) then
         if (present(error)) call fatal_error(error, &
                                              "Unable to initialize HDF5 Fortran interface")
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf5_ensure_initialized

   !> Open an HDF5 file for reading or writing.
   !> @param[in,out] self     HDF5 file instance
   !> @param[in]     filename Path to the HDF5 file
   !> @param[in]     status   File status: 'old', 'new', or 'replace' (default: 'old')
   !> @param[in]     action   Access mode: 'read'/'r', 'write'/'w', or 'readwrite'/'rw' (default: 'rw')
   !> @param[out]    error    Error handler
   subroutine hdf_open_file(self, filename, status, action, error)
      class(hdf5_file), intent(inout) :: self
      character(*), intent(in) :: filename
      character(*), intent(in), optional :: status
      character(*), intent(in), optional :: action
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Normalized status string
      character(:), allocatable :: lstatus
      !> Normalized action string
      character(:), allocatable :: laction
      !> HDF5 error code
      integer :: ierr

      call hdf5_ensure_initialized(error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      self%filename = filename

      lstatus = 'old'
      if (present(status)) lstatus = to_lower(status)

      laction = 'rw'
      if (present(action)) laction = to_lower(action)

      select case (lstatus)
      case ('old')
         select case (laction)
         case ('read', 'r')
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, self%lid, ierr)
         case ('write', 'readwrite', 'w', 'rw')
            call h5fopen_f(filename, H5F_ACC_RDWR_F, self%lid, ierr)
         case default
            if (present(error)) call fatal_error(error, "Unsupported action: "//laction)
            return
         end select
         if (ierr /= 0) then
            if (present(error)) call fatal_error(error, "Failed to open HDF5 file: "//filename)
            return
         end if
      case ('new', 'replace')
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, self%lid, ierr)
         if (ierr /= 0) then
            if (present(error)) call fatal_error(error, "Failed to create HDF5 file: "//filename)
            return
         end if
      case default
         if (present(error)) call fatal_error(error, "Unsupported status: "//lstatus)
         return
      end select

      self%is_open = .true.
      self%group_open = .false.
      self%gid = -1
      self%glid = -1
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_open_file

   !> Close the HDF5 file.
   !> @param[in,out] self  HDF5 file instance
   !> @param[out]    error Error handler
   subroutine hdf_close_file(self, error)
      class(hdf5_file), intent(inout) :: self
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr
      !> File identifier to close
      integer(HID_T) :: fid

      if (.not. self%is_open) return

      fid = self%lid
      if (self%group_open) then
         call h5gclose_f(self%gid, ierr)
         if (ierr /= 0) then
            if (present(error)) call fatal_error(error, &
                                                 "Unable to close open HDF5 group before closing file: "//self%filename)
            return
         end if
         fid = self%glid
         self%group_open = .false.
         self%gid = -1
         self%glid = -1
         self%lid = fid
      end if

      call h5fclose_f(fid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Unable to close HDF5 file: "//self%filename)
         return
      end if

      self%is_open = .false.
      self%lid = -1
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_close_file

   !> Open a group within the HDF5 file for subsequent operations.
   !> @param[in,out] self  HDF5 file instance
   !> @param[in]     gname Name of the group to open
   !> @param[out]    error Error handler
   subroutine hdf_open_group(self, gname, error)
      class(hdf5_file), intent(inout) :: self
      character(*), intent(in) :: gname
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "open a group", error)) return
      if (self%group_open) then
         if (present(error)) call fatal_error(error, &
                                              "Cannot open a nested HDF5 group before closing the current group")
         return
      end if

      call h5gopen_f(self%lid, gname, self%gid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to open group: "//gname)
         return
      end if
      self%glid = self%lid
      self%lid = self%gid
      self%group_open = .true.
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_open_group

   !> Close the currently open group and restore the previous location.
   !> @param[in,out] self  HDF5 file instance
   !> @param[out]    error Error handler
   subroutine hdf_close_group(self, error)
      class(hdf5_file), intent(inout) :: self
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "close a group", error)) return
      if (.not. self%group_open) return

      call h5gclose_f(self%gid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to close group")
         return
      end if
      self%lid = self%glid
      self%gid = -1
      self%glid = -1
      self%group_open = .false.
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_close_group

   !> Ensure the HDF5 file is open before using its handle.
   !> @param[in]  self    HDF5 file instance
   !> @param[in]  action  Description of the attempted operation
   !> @param[out] error   Error handler
   !> @return True when the file is open
   logical function hdf_require_open(self, action, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: action
      type(error_type), allocatable, intent(out), optional :: error

      hdf_require_open = self%is_open
      if (.not. hdf_require_open) then
         if (present(error)) call fatal_error(error, &
                                              "Cannot "//action//" on a closed HDF5 file")
      end if

   end function hdf_require_open

   !> Ensure dataset rank matches expected rank.
   !> @param[in]  self          HDF5 file instance
   !> @param[in]  dname         Dataset name
   !> @param[in]  expected_rank Expected dataset rank
   !> @param[out] error         Error handler
   !> @return True when rank matches expected rank
   logical function hdf_expect_rank(self, dname, expected_rank, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: expected_rank
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset identifier
      integer(HID_T) :: did
      !> Dataspace identifier
      integer(HID_T) :: sid
      !> HDF5 error code
      integer :: ierr
      !> Actual dataset rank
      integer :: actual_rank

      hdf_expect_rank = .false.

      call h5dopen_f(self%lid, dname, did, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to open dataset: "//dname)
         return
      end if

      call h5dget_space_f(did, sid, ierr)
      if (ierr /= 0) then
         call h5dclose_f(did, ierr)
         if (present(error)) call fatal_error(error, "Failed to get dataspace for: "//dname)
         return
      end if

      call h5sget_simple_extent_ndims_f(sid, actual_rank, ierr)
      if (ierr /= 0) then
         call h5sclose_f(sid, ierr)
         call h5dclose_f(did, ierr)
         if (present(error)) call fatal_error(error, "Failed to get dataset rank for: "//dname)
         return
      end if

      call h5sclose_f(sid, ierr)
      call h5dclose_f(did, ierr)

      if (actual_rank /= expected_rank) then
         if (present(error)) call fatal_error(error, "Dataset rank mismatch for: "//dname)
         return
      end if

      hdf_expect_rank = .true.
#else
      hdf_expect_rank = .false.
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end function hdf_expect_rank

   !> Check whether a dataset or group exists in the file.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  name  Name of the dataset or group to check
   !> @param[out] exist True if the object exists
   !> @param[out] error Error handler
   subroutine hdf_exist(self, name, exist, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: name
      logical, intent(out) :: exist
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "check existence", error)) return

      call h5lexists_f(self%lid, name, exist, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to check existence of: "//name)
         return
      end if
#else
      exist = .false.
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_exist

   !> Delete a dataset or group from the file.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  name  Name of the dataset or group to delete
   !> @param[out] error Error handler
   subroutine hdf_delete(self, name, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: name
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "delete", error)) return

      call h5ldelete_f(self%lid, name, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to delete: "//name)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_delete

   !> Ensure parent groups exist for a given path (internal helper).
   !> Creates only intermediate groups, not the final path component.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  path  Full path (e.g., '/data/results/dataset')
   !> @param[out] error Error handler
   subroutine hdf_ensure_parent_groups(self, path, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: path
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Temporary group identifier
      integer(HID_T) :: gid
      !> HDF5 error code
      integer :: ierr
      !> Start position in path parsing
      integer :: sp
      !> End position in path parsing
      integer :: ep
      !> String length of path
      integer :: sl
      !> Flag indicating group existence
      logical :: gexist

      if (.not. hdf_require_open(self, "create a group", error)) return

      sl = len(path)
      sp = 1

      do
         ep = index(path(sp + 1:sl), "/")
         if (ep == 0) exit
         sp = sp + ep
         call h5lexists_f(self%lid, path(1:sp - 1), gexist, ierr)
         if (ierr /= 0) then
            if (present(error)) then
               call fatal_error(error, "Failed to check group existence: "// &
                                path(1:sp - 1))
            end if
            return
         end if
         if (.not. gexist) then
            call h5gcreate_f(self%lid, path(1:sp - 1), gid, ierr)
            if (ierr /= 0) then
               if (present(error)) call fatal_error(error, &
                                                    "Failed to create group: "//path(1:sp - 1))
               return
            end if
            call h5gclose_f(gid, ierr)
            if (ierr /= 0) then
               if (present(error)) then
                  call fatal_error(error, "Failed to close group: "// &
                                   path(1:sp - 1))
               end if
               return
            end if
         end if
      end do
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_ensure_parent_groups

   !> Create a group, creating parent groups as needed.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  gname Full path of the group to create (e.g., '/data/results')
   !> @param[out] error Error handler
   subroutine hdf_add_group(self, gname, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: gname
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Temporary group identifier
      integer(HID_T) :: gid
      !> HDF5 error code
      integer :: ierr
      !> String length of group name
      integer :: sl
      !> Flag indicating group existence
      logical :: gexist

      if (.not. hdf_require_open(self, "add a group", error)) return

      call self%hdf_ensure_parent_groups(gname, error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      sl = len(gname)
      if (gname(1:1) == "/" .and. sl > 1) then
         call h5lexists_f(self%lid, gname, gexist, ierr)
         if (ierr /= 0) then
            if (present(error)) then
               call fatal_error(error, "Failed to check final group existence: "//gname)
            end if
            return
         end if
         if (.not. gexist) then
            call h5gcreate_f(self%lid, gname, gid, ierr)
            if (ierr /= 0) then
               if (present(error)) call fatal_error(error, &
                                                    "Failed to create final group: "//gname)
               return
            end if
            call h5gclose_f(gid, ierr)
            if (ierr /= 0) then
               if (present(error)) then
                  call fatal_error(error, "Failed to close final group: "//gname)
               end if
               return
            end if
         end if
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_group

   !> Add a string attribute to a group or dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  path  Path to the group or dataset
   !> @param[in]  name  Attribute name
   !> @param[in]  value Attribute value
   !> @param[out] error Error handler
   subroutine hdf_adda_string(self, path, name, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: path
      character(*), intent(in) :: name
      character(*), intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr
      !> Flag indicating path exists
      logical :: path_exists

      if (.not. hdf_require_open(self, "add an attribute", error)) return

      call self%exist(path, path_exists, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if
      if (.not. path_exists) then
         if (present(error)) call fatal_error(error, &
                                              "Failed to set attribute on missing path: "//path)
         return
      end if
      call h5ltset_attribute_string_f(self%lid, path, name, value, ierr)
      if (ierr /= 0) then
         if (present(error)) then
            call fatal_error(error, &
                             "Failed to set attribute: "//name//" in "//path)
         end if
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_adda_string

   !> Get a string attribute from a group or dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  path  Path to the group or dataset
   !> @param[in]  name  Attribute name
   !> @param[out] value Allocatable string to receive the attribute value
   !> @param[out] error Error handler
   subroutine hdf_geta_string(self, path, name, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: path
      character(*), intent(in) :: name
      character(:), intent(out), allocatable :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr
      !> Attribute buffer
      character(len=1024) :: buffer

      if (.not. hdf_require_open(self, "get an attribute", error)) return

      call h5ltget_attribute_string_f(self%lid, path, name, buffer, ierr)
      if (ierr /= 0) then
         if (present(error)) then
            call fatal_error(error, &
                             "Failed to get attribute: "//name//" in "//path)
         end if
         return
      end if

      value = trim(buffer)
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_geta_string

   !> Write a scalar integer dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value Scalar integer value
   !> @param[out] error Error handler
   subroutine hdf_add_int(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataspace identifier
      integer(HID_T) :: sid
      !> Dataset identifier
      integer(HID_T) :: did
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5screate_f(H5S_SCALAR_F, sid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, &
                                              "Failed to create scalar dataspace for: "//dname)
         return
      end if

      call h5dcreate_f(self%lid, dname, &
                       h5kind_to_type(kind(value), H5_INTEGER_KIND), sid, did, ierr)
      if (ierr /= 0) then
         call h5sclose_f(sid, ierr)
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if

      call h5dwrite_f(did, h5kind_to_type(kind(value), H5_INTEGER_KIND), &
                      value, int(shape(value), HSIZE_T), ierr)
      if (ierr /= 0) then
         call h5dclose_f(did, ierr)
         call h5sclose_f(sid, ierr)
         if (present(error)) call fatal_error(error, "Failed to write dataset: "//dname)
         return
      end if

      call h5dclose_f(did, ierr)
      if (ierr /= 0) then
         call h5sclose_f(sid, ierr)
         if (present(error)) call fatal_error(error, "Failed to close dataset: "//dname)
         return
      end if

      call h5sclose_f(sid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to close dataspace for: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_int

   !> Write a 1D integer array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value 1D integer array
   !> @param[out] error Error handler
   subroutine hdf_add_int1d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5ltmake_dataset_f(self%lid, dname, &
                              rank(value), int(shape(value), HSIZE_T), &
                              h5kind_to_type(kind(value), H5_INTEGER_KIND), value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_int1d

   !> Write a 2D integer array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value 2D integer array
   !> @param[out] error Error handler
   subroutine hdf_add_int2d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5ltmake_dataset_f(self%lid, dname, &
                              rank(value), int(shape(value), HSIZE_T), &
                              h5kind_to_type(kind(value), H5_INTEGER_KIND), value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_int2d

   !> Write a 3D integer array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value 3D integer array
   !> @param[out] error Error handler
   subroutine hdf_add_int3d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5ltmake_dataset_f(self%lid, dname, &
                              rank(value), int(shape(value), HSIZE_T), &
                              h5kind_to_type(kind(value), H5_INTEGER_KIND), value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_int3d

   !> Write a scalar real dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value Scalar real value
   !> @param[out] error Error handler
   subroutine hdf_add_real(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataspace identifier
      integer(HID_T) :: sid
      !> Dataset identifier
      integer(HID_T) :: did
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5screate_f(H5S_SCALAR_F, sid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, &
                                              "Failed to create scalar dataspace for: "//dname)
         return
      end if

      call h5dcreate_f(self%lid, dname, &
                       h5kind_to_type(kind(value), H5_REAL_KIND), sid, did, ierr)
      if (ierr /= 0) then
         call h5sclose_f(sid, ierr)
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if

      call h5dwrite_f(did, h5kind_to_type(kind(value), H5_REAL_KIND), &
                      value, int(shape(value), HSIZE_T), ierr)
      if (ierr /= 0) then
         call h5dclose_f(did, ierr)
         call h5sclose_f(sid, ierr)
         if (present(error)) call fatal_error(error, "Failed to write dataset: "//dname)
         return
      end if

      call h5dclose_f(did, ierr)
      if (ierr /= 0) then
         call h5sclose_f(sid, ierr)
         if (present(error)) call fatal_error(error, "Failed to close dataset: "//dname)
         return
      end if

      call h5sclose_f(sid, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to close dataspace for: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_real

   !> Write a 1D real array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value 1D real array
   !> @param[out] error Error handler
   subroutine hdf_add_real1d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5ltmake_dataset_f(self%lid, dname, &
                              rank(value), int(shape(value), HSIZE_T), &
                              h5kind_to_type(kind(value), H5_REAL_KIND), value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_real1d

   !> Write a 2D real array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value 2D real array
   !> @param[out] error Error handler
   subroutine hdf_add_real2d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5ltmake_dataset_f(self%lid, dname, &
                              rank(value), int(shape(value), HSIZE_T), &
                              h5kind_to_type(kind(value), H5_REAL_KIND), value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_real2d

   !> Write a 3D real array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value 3D real array
   !> @param[out] error Error handler
   subroutine hdf_add_real3d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      call h5ltmake_dataset_f(self%lid, dname, &
                              rank(value), int(shape(value), HSIZE_T), &
                              h5kind_to_type(kind(value), H5_REAL_KIND), value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_real3d

   !> Write a string dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[in]  value String value
   !> @param[out] error Error handler
   subroutine hdf_add_string(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      character(*), intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "add a string dataset", error)) return

      call self%hdf_ensure_parent_groups(dname, error=error)
      if (present(error)) then
         if (allocated(error)) return
      end if
      call h5ltmake_dataset_string_f(self%lid, dname, value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to create string dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_add_string

   !> Read a string dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable string to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_string(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      character(:), intent(out), allocatable :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(1)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a string dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 0, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (character(dsize - 1) :: value)
      call h5ltread_dataset_string_f(self%lid, dname, value, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read string dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_string

   !> Read a scalar integer dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Scalar integer to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_int(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out) :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset identifier
      integer(HID_T) :: set_id
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 0, error)) return

      call h5dopen_f(self%lid, dname, set_id, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to open dataset: "//dname)
         return
      end if

      call h5dread_f(set_id, h5kind_to_type(kind(value), H5_INTEGER_KIND), &
                     value, int(shape(value), HSIZE_T), ierr)
      if (ierr /= 0) then
         call h5dclose_f(set_id, ierr)
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if

      call h5dclose_f(set_id, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to close dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_int

   !> Read a 1D integer array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 1D integer array to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_int1d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out), allocatable :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(1)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 1, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (value(dims(1)))

      call h5ltread_dataset_f(self%lid, dname, &
                              h5kind_to_type(kind(value), H5_INTEGER_KIND), value, dims, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_int1d

   !> Read a 2D integer array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 2D integer array to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_int2d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out), allocatable :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(2)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 2, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (value(dims(1), dims(2)))

      call h5ltread_dataset_f(self%lid, dname, &
                              h5kind_to_type(kind(value), H5_INTEGER_KIND), value, dims, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_int2d

   !> Read a 3D integer array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 3D integer array to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_int3d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out), allocatable :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(3)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 3, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (value(dims(1), dims(2), dims(3)))

      call h5ltread_dataset_f(self%lid, dname, &
                              h5kind_to_type(kind(value), H5_INTEGER_KIND), value, dims, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_int3d

   !> Read a scalar real dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Scalar real to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_real(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset identifier
      integer(HID_T) :: set_id
      !> HDF5 error code
      integer :: ierr

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 0, error)) return

      call h5dopen_f(self%lid, dname, set_id, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to open dataset: "//dname)
         return
      end if

      call h5dread_f(set_id, h5kind_to_type(kind(value), H5_REAL_KIND), &
                     value, int(shape(value), HSIZE_T), ierr)
      if (ierr /= 0) then
         call h5dclose_f(set_id, ierr)
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if

      call h5dclose_f(set_id, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to close dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_real

   !> Read a 1D real array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 1D real array to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_real1d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out), allocatable :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(1)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 1, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (value(dims(1)))

      call h5ltread_dataset_f(self%lid, dname, &
                              h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_real1d

   !> Read a 2D real array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 2D real array to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_real2d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out), allocatable :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(2)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 2, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (value(dims(1), dims(2)))

      call h5ltread_dataset_f(self%lid, dname, &
                              h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_real2d

   !> Read a 3D real array dataset.
   !> @param[in]  self  HDF5 file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 3D real array to receive the data
   !> @param[out] error Error handler
   subroutine hdf_get_real3d(self, dname, value, error)
      class(hdf5_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out), allocatable :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

#ifdef WITH_HDF5
      !> Dataset dimensions
      integer(HSIZE_T) :: dims(3)
      !> Dataset size in bytes
      integer(SIZE_T) :: dsize
      !> HDF5 error code
      integer :: ierr
      !> Dataset type
      integer :: dtype

      if (.not. hdf_require_open(self, "get a dataset", error)) return
      if (.not. hdf_expect_rank(self, dname, 3, error)) return

      call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to get dataset info for: "//dname)
         return
      end if

      allocate (value(dims(1), dims(2), dims(3)))

      call h5ltread_dataset_f(self%lid, dname, &
                              h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
      if (ierr /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
#else
      if (present(error)) call fatal_error(error, "moist is built without HDF5 support")
#endif

   end subroutine hdf_get_real3d

   !> Convert a string to lowercase.
   !> @param[in] str Input string
   !> @return Lowercase version of the input string
   elemental function to_lower(str)
      character(*), intent(in) :: str
      character(len(str)) :: to_lower
      !> Lowercase alphabet characters
      character(*), parameter :: lower = "abcdefghijklmnopqrstuvwxyz", &
                                 !> Uppercase alphabet characters
                                 upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      !> Loop indices
      integer :: i, j

      to_lower = str

      do concurrent(i=1:len(str))
         j = index(upper, str(i:i))
         if (j > 0) to_lower(i:i) = lower(j:j)
      end do

   end function to_lower

end module moist_utils_hdf5io
