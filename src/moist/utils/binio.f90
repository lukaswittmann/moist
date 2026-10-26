!> Self-describing binary container files built on intrinsic stream I/O
!>
!> This module provides a high-level object-oriented interface to a compact
!> binary container format, supporting scalar and multidimensional arrays of
!> integers, reals, and strings addressed by hierarchical path names such as
!> `/grid/r`.  The [[binary_file]] type encapsulates the unit and the in-memory
!> index and provides generic bindings for type-transparent read/write
!> operations
!>
!> Payloads are written with unformatted stream `write` statements, so a whole
!> array leaves the program in a single contiguous transfer with no conversion
!> and no chunking layer in between
!>
!> File layout:
!>
!> ```
!> [header, 28 bytes]  magic(8) endian(int32) version(int32) real_bytes(int32)
!>                     index_offset(int64)
!> [payloads]          raw dataset bytes in creation order, contiguous
!> [index]             entry count(int32) followed by one record per entry:
!>                     type(int32) rank(int32) dims(3 x int64) offset(int64)
!>                     nbytes(int64) path(int32 length + chars)
!>                     attributes(int32 count + name/value length-prefixed pairs)
!> ```
!>
module moist_utils_binio
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io_utils, only: to_lower
   implicit none(type, external)

   public :: binary_file

   private

   !> Magic bytes identifying a moist binary container
   character(len=*), parameter :: file_magic = "MOISTBIN"

   !> Container format version written into the header
   integer(int32), parameter :: format_version = 1_int32

   !> Byte-order probe: read back unchanged only on a matching platform
   integer(int32), parameter :: endian_probe = 16909060_int32

   !> Entry holds a group rather than data
   integer(int32), parameter :: type_group = 0_int32
   !> Entry holds default integer data
   integer(int32), parameter :: type_int = 1_int32
   !> Entry holds real(wp) data
   integer(int32), parameter :: type_real = 2_int32
   !> Entry holds character data
   integer(int32), parameter :: type_string = 3_int32

   !> Size of one stored integer element in bytes
   integer(int64), parameter :: int_bytes = storage_size(0_int32)/8
   !> Size of one stored real element in bytes
   integer(int64), parameter :: real_bytes = storage_size(1.0_wp)/8

   !> Total size of the file header in bytes
   integer(int64), parameter :: header_bytes = 28_int64
   !> Stream position of the index offset field inside the header
   integer(int64), parameter :: index_offset_pos = 21_int64
   !> Maximum supported dataset rank
   integer, parameter :: max_rank = 3
   !> Size of an index record with an empty path and no attributes in bytes
   integer(int64), parameter :: min_entry_bytes = 56_int64

   !> Named string attribute attached to a group or dataset
   type :: attribute_type
      !> Attribute name
      character(:), allocatable :: name
      !> Attribute value
      character(:), allocatable :: value
   end type attribute_type

   !> Index record describing one group or dataset in the container
   type :: entry_type
      !> Absolute path of the entry, e.g. '/grid/r'
      character(:), allocatable :: path
      !> Payload type, one of the `type_*` constants
      integer(int32) :: dtype = type_group
      !> Number of dimensions (0 for scalars, groups and strings)
      integer(int32) :: rank = 0_int32
      !> Extent along each dimension, in Fortran (column-major) order
      integer(int64) :: dims(max_rank) = 0_int64
      !> Stream position of the first payload byte (1-based)
      integer(int64) :: offset = 0_int64
      !> Payload size in bytes
      integer(int64) :: nbytes = 0_int64
      !> String attributes attached to this entry
      type(attribute_type), allocatable :: attrs(:)
   end type entry_type

   !> Binary container handle for reading and writing datasets
   !>
   !> Provides methods to open/close files, manage groups, and read/write
   !> datasets of various types (integers, reals, strings) and ranks (0-3D)
   type :: binary_file

      !> Path to the container file
      character(:), allocatable :: filename
      !> Flag indicating whether the file is currently open
      logical :: is_open = .false.
      !> Flag indicating whether a group is currently open
      logical :: group_open = .false.
      !> Connected unit number, -1 when closed
      integer, private :: unit = -1
      !> Flag indicating whether the file may be modified
      logical, private :: writable = .false.
      !> Path prefix applied to relative names while a group is open
      character(:), allocatable, private :: prefix
      !> Index of groups and datasets held in memory while the file is open
      type(entry_type), allocatable, private :: entries(:)
      !> Number of used slots in `entries`
      integer, private :: nentries = 0
      !> Stream position where the next payload is written
      integer(int64), private :: next_pos = 0_int64

   contains
      !> Open a container file for reading or writing
      procedure :: open => bin_open_file
      !> Close the container file
      procedure :: close => bin_close_file

      !> Open a group within the file
      procedure :: open_group => bin_open_group
      !> Close the current group
      procedure :: close_group => bin_close_group

      !> Check if a dataset or group exists
      procedure :: exist => bin_exist

      !> Extents of a dataset, without reading its payload
      procedure :: get_shape => bin_get_shape

      !> Delete a dataset or group
      procedure :: delete => bin_delete

      !> Add a group, dataset, or attribute
      generic :: add => bin_add_group, &
         bin_add_int, &
         bin_add_int1d, &
         bin_add_int2d, &
         bin_add_int3d, &
         bin_add_real, &
         bin_add_real1d, &
         bin_add_real2d, &
         bin_add_real3d, &
         bin_add_string

      !> Retrieve a dataset value
      generic :: get => bin_get_int, &
         bin_get_int1d, &
         bin_get_int2d, &
         bin_get_int3d, &
         bin_get_real, &
         bin_get_real1d, &
         bin_get_real2d, &
         bin_get_real3d, &
         bin_get_string

      !> Add a string attribute to a path
      generic :: adda => bin_adda_string

      !> Get a string attribute from a path
      generic :: geta => bin_geta_string

      procedure, private :: bin_add_group
      procedure, private :: bin_ensure_parent_groups
      procedure, private :: bin_add_int
      procedure, private :: bin_get_int
      procedure, private :: bin_add_int1d
      procedure, private :: bin_get_int1d
      procedure, private :: bin_add_int2d
      procedure, private :: bin_get_int2d
      procedure, private :: bin_add_int3d
      procedure, private :: bin_get_int3d
      procedure, private :: bin_add_real
      procedure, private :: bin_get_real
      procedure, private :: bin_add_real1d
      procedure, private :: bin_get_real1d
      procedure, private :: bin_add_real2d
      procedure, private :: bin_get_real2d
      procedure, private :: bin_add_real3d
      procedure, private :: bin_get_real3d
      procedure, private :: bin_get_string
      procedure, private :: bin_add_string
      procedure, private :: bin_adda_string
      procedure, private :: bin_geta_string
   end type binary_file

contains

   !> Open a container file for reading or writing
   !>
   !> @param[in,out] self     Container file instance
   !> @param[in]     filename Path to the container file
   !> @param[in]     status   File status: "old", or "new"/"replace" (both
   !>                         truncate an existing file); default "old"
   !> @param[in]     action   Access mode: "read"/"r", "write"/"w", or
   !>                         "readwrite"/"rw" (default: "readwrite")
   !> @param[out]    error    Error handler
   subroutine bin_open_file(self, filename, status, action, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: filename
      character(*), intent(in), optional :: status
      character(*), intent(in), optional :: action
      type(error_type), allocatable, intent(out), optional :: error

      !> Normalized status string
      character(:), allocatable :: lstatus
      !> Normalized action string
      character(:), allocatable :: laction
      !> I/O status code
      integer :: stat
      !> Unit connected to the file
      integer :: unit

      lstatus = "old"
      if (present(status)) lstatus = to_lower(status)

      laction = "rw"
      if (present(action)) laction = to_lower(action)

      select case (lstatus)
      case ("old")
         select case (laction)
         case ("read", "r")
            self%writable = .false.
         case ("write", "readwrite", "w", "rw")
            self%writable = .true.
         case default
            if (present(error)) call fatal_error(error, "Unsupported action: "//laction)
            return
         end select
      case ("new", "replace")
         self%writable = .true.
      case default
         if (present(error)) call fatal_error(error, "Unsupported status: "//lstatus)
         return
      end select

      if (self%writable) then
         open (newunit=unit, file=filename, status=merge("replace", "old    ", lstatus /= "old"), &
               access="stream", form="unformatted", action="readwrite", iostat=stat)
      else
         open (newunit=unit, file=filename, status="old", access="stream", &
               form="unformatted", action="read", iostat=stat)
      end if
      if (stat /= 0) then
         if (present(error)) then
            if (lstatus == "old") then
               call fatal_error(error, "Failed to open binary file: "//filename)
            else
               call fatal_error(error, "Failed to create binary file: "//filename)
            end if
         end if
         return
      end if

      self%unit = unit
      self%filename = filename
      self%prefix = ""
      self%nentries = 0
      if (allocated(self%entries)) deallocate (self%entries)
      allocate (self%entries(16))

      if (lstatus == "old") then
         call bin_read_index(self, error)
      else
         call bin_write_header(self, 0_int64, error)
         self%next_pos = header_bytes + 1_int64
      end if

      if (present(error)) then
         if (allocated(error)) then
            !> The failure is already reported; drop the connection and leave
            !> the handle closed
            close (self%unit, iostat=stat)
            self%unit = -1
            return
         end if
      end if

      self%is_open = .true.
      self%group_open = .false.

   end subroutine bin_open_file

   !> Close the container file, committing the index when it is writable
   !>
   !> @param[in,out] self  Container file instance
   !> @param[out]    error Error handler
   subroutine bin_close_file(self, error)
      class(binary_file), intent(inout) :: self
      type(error_type), allocatable, intent(out), optional :: error

      !> I/O status code
      integer :: stat

      if (.not. self%is_open) return

      if (self%writable) call bin_write_index(self, error)

      close (self%unit, iostat=stat)
      if (stat /= 0 .and. present(error)) then
         !> An index-write failure is the more useful diagnostic, so it wins
         if (.not. allocated(error)) then
            call fatal_error(error, "Unable to close binary file: "//self%filename)
         end if
      end if

      call bin_reset(self)

   end subroutine bin_close_file

   !> Drop all state associated with an open file
   !>
   !> @param[in,out] self Container file instance
   subroutine bin_reset(self)
      class(binary_file), intent(inout) :: self

      self%is_open = .false.
      self%group_open = .false.
      self%writable = .false.
      self%unit = -1
      self%prefix = ""
      self%nentries = 0
      self%next_pos = 0_int64
      if (allocated(self%entries)) deallocate (self%entries)

   end subroutine bin_reset

   !> Write the fixed-size file header
   !>
   !> @param[in]  self         Container file instance
   !> @param[in]  index_offset Stream position of the index block, 0 while incomplete
   !> @param[out] error        Error handler
   subroutine bin_write_header(self, index_offset, error)
      class(binary_file), intent(in) :: self
      integer(int64), intent(in) :: index_offset
      type(error_type), allocatable, intent(out), optional :: error

      !> I/O status code
      integer :: stat

      write (self%unit, pos=1, iostat=stat) file_magic, endian_probe, format_version, &
         int(real_bytes, int32), index_offset
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to write header of: "//self%filename)
      end if

   end subroutine bin_write_header

   !> Read and validate the header, then load the index into memory
   !>
   !> @param[in,out] self  Container file instance
   !> @param[out]    error Error handler
   subroutine bin_read_index(self, error)
      class(binary_file), intent(inout) :: self
      type(error_type), allocatable, intent(out), optional :: error

      !> Magic bytes read back from the file
      character(len=len(file_magic)) :: magic
      !> Byte-order probe read back from the file
      integer(int32) :: endian
      !> Format version read back from the file
      integer(int32) :: version
      !> Stored real width in bytes
      integer(int32) :: stored_real_bytes
      !> Stream position of the index block
      integer(int64) :: index_offset
      !> Size of the file in bytes
      integer(int64) :: file_bytes
      !> Bytes from the index block to the end of the file
      integer(int64) :: index_bytes
      !> Number of index records
      integer(int32) :: count
      !> Loop index over entries
      integer :: i
      !> I/O status code
      integer :: stat

      read (self%unit, pos=1, iostat=stat) magic, endian, version, stored_real_bytes, index_offset
      if (stat /= 0 .or. magic /= file_magic) then
         if (present(error)) call fatal_error(error, "Not a moist binary file: "//self%filename)
         return
      end if
      if (endian /= endian_probe) then
         if (present(error)) call fatal_error(error, "Byte order mismatch in: "//self%filename)
         return
      end if
      if (version /= format_version) then
         if (present(error)) then
            call fatal_error(error, "Unsupported format version in: "//self%filename)
         end if
         return
      end if
      if (int(stored_real_bytes, int64) /= real_bytes) then
         if (present(error)) call fatal_error(error, "Real kind mismatch in: "//self%filename)
         return
      end if
      if (index_offset <= 0_int64) then
         if (present(error)) call fatal_error(error, "Incomplete binary file: "//self%filename)
         return
      end if

      inquire (unit=self%unit, size=file_bytes, iostat=stat)
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to size binary file: "//self%filename)
         return
      end if
      if (index_offset <= header_bytes .or. index_offset > file_bytes - int_bytes + 1_int64) then
         if (present(error)) call fatal_error(error, "Corrupt index in: "//self%filename)
         return
      end if
      index_bytes = file_bytes - index_offset + 1_int64

      read (self%unit, pos=index_offset, iostat=stat) count
      if (stat /= 0 .or. count < 0_int32 .or. &
          int(count, int64)*min_entry_bytes > index_bytes - int_bytes) then
         if (present(error)) call fatal_error(error, "Corrupt index in: "//self%filename)
         return
      end if

      if (allocated(self%entries)) deallocate (self%entries)
      allocate (self%entries(max(int(count), 16)))
      self%nentries = int(count)

      do i = 1, self%nentries
         call bin_read_entry(self%unit, index_bytes, self%entries(i), stat)
         if (stat == 0) then
            if (.not. valid_entry(self%entries(i), index_offset)) stat = -1
         end if
         if (stat /= 0) then
            self%nentries = 0
            if (present(error)) call fatal_error(error, "Corrupt index in: "//self%filename)
            return
         end if
      end do

      self%next_pos = file_bytes + 1_int64

   end subroutine bin_read_index

   !> Whether an index record describes a payload the file can hold
   !>
   !> - Groups carry no payload; strings are rank 0 with one byte per character
   !> - Arrays hold the product of their extents in elements of their type
   !> - Payloads lie between the header and the index
   !>
   !> @param[in] ent          Index record
   !> @param[in] index_offset Stream position of the index block
   !> @return True when sizes, extents and position agree
   pure logical function valid_entry(ent, index_offset)
      type(entry_type), intent(in) :: ent
      integer(int64), intent(in) :: index_offset

      !> Elements of the payload
      integer(int64) :: count
      !> Bytes per element
      integer(int64) :: width
      !> Loop index over dimensions
      integer :: i

      valid_entry = .false.
      if (ent%rank < 0_int32 .or. ent%rank > int(max_rank, int32)) return
      if (any(ent%dims < 0_int64) .or. any(ent%dims(ent%rank + 1:) /= 0_int64)) return
      if (ent%nbytes < 0_int64) return

      select case (ent%dtype)
      case (type_group)
         valid_entry = ent%rank == 0_int32 .and. ent%nbytes == 0_int64
         return
      case (type_string)
         if (ent%rank /= 0_int32 .or. ent%nbytes > int(huge(0), int64)) return
         count = ent%nbytes
         width = 1_int64
      case (type_int)
         width = int_bytes
      case (type_real)
         width = real_bytes
      case default
         return
      end select

      if (ent%dtype /= type_string) then
         count = 1_int64
         do i = 1, int(ent%rank)
            if (ent%dims(i) > 0_int64 .and. count > huge(count)/ent%dims(i)) return
            count = count*ent%dims(i)
         end do
      end if
      if (count > huge(count)/width .or. ent%nbytes /= count*width) return

      if (ent%nbytes > 0_int64) then
         if (ent%offset <= header_bytes .or. ent%offset > index_offset - ent%nbytes) return
      end if
      valid_entry = .true.

   end function valid_entry

   !> Read one index record from the current stream position
   !>
   !> Lengths and counts beyond the index block are refused before allocating
   !>
   !> @param[in]  unit        Connected unit
   !> @param[in]  index_bytes Bytes from the index block to the end of the file
   !> @param[out] ent         Index record
   !> @param[out] stat        I/O status code, non-zero on failure
   subroutine bin_read_entry(unit, index_bytes, ent, stat)
      integer, intent(in) :: unit
      integer(int64), intent(in) :: index_bytes
      type(entry_type), intent(out) :: ent
      integer, intent(out) :: stat

      !> Number of attributes on the entry
      integer(int32) :: nattr
      !> Loop index over attributes
      integer :: i

      read (unit, iostat=stat) ent%dtype, ent%rank, ent%dims, ent%offset, ent%nbytes
      if (stat /= 0) return

      call read_string(ent%path)
      if (stat /= 0) return

      read (unit, iostat=stat) nattr
      if (stat /= 0) return
      if (nattr < 0_int32 .or. 2_int64*int_bytes*int(nattr, int64) > index_bytes) then
         stat = -1
         return
      end if
      allocate (ent%attrs(int(nattr)))
      do i = 1, int(nattr)
         call read_string(ent%attrs(i)%name)
         if (stat /= 0) return
         call read_string(ent%attrs(i)%value)
         if (stat /= 0) return
      end do

   contains

      !> Read one length-prefixed string of the record
      !>
      !> @param[out] str String read from the stream
      subroutine read_string(str)
         character(:), allocatable, intent(out) :: str

         !> Length of the string
         integer(int32) :: slen

         read (unit, iostat=stat) slen
         if (stat /= 0) return
         if (slen < 0_int32 .or. int(slen, int64) > index_bytes) then
            stat = -1
            return
         end if
         allocate (character(int(slen)) :: str)
         if (slen > 0_int32) read (unit, iostat=stat) str
      end subroutine read_string

   end subroutine bin_read_entry

   !> Append the index at the end of the file and commit it in the header
   !>
   !> @param[in,out] self  Container file instance
   !> @param[out]    error Error handler
   subroutine bin_write_index(self, error)
      class(binary_file), intent(in) :: self
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the index block
      integer(int64) :: index_offset
      !> Loop indices over entries and attributes
      integer :: i, j
      !> I/O status code
      integer :: stat

      index_offset = self%next_pos

      write (self%unit, pos=index_offset, iostat=stat) int(self%nentries, int32)
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to write index of: "//self%filename)
         return
      end if

      do i = 1, self%nentries
         associate (ent => self%entries(i))
            write (self%unit, iostat=stat) ent%dtype, ent%rank, ent%dims, ent%offset, ent%nbytes
            if (stat == 0) write (self%unit, iostat=stat) int(len(ent%path), int32), ent%path
            if (stat == 0) write (self%unit, iostat=stat) int(size(ent%attrs), int32)
            do j = 1, size(ent%attrs)
               if (stat /= 0) exit
               write (self%unit, iostat=stat) &
                  int(len(ent%attrs(j)%name), int32), ent%attrs(j)%name, &
                  int(len(ent%attrs(j)%value), int32), ent%attrs(j)%value
            end do
         end associate
         if (stat /= 0) then
            if (present(error)) call fatal_error(error, "Failed to write index of: "//self%filename)
            return
         end if
      end do

      call bin_write_header(self, index_offset, error)

   end subroutine bin_write_index

   !> Ensure the file is open before using its handle
   !>
   !> @param[in]  self    Container file instance
   !> @param[in]  action  Description of the attempted operation
   !> @param[out] error   Error handler
   !> @return True when the file is open
   logical function bin_require_open(self, action, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: action
      type(error_type), allocatable, intent(out), optional :: error

      bin_require_open = self%is_open
      if (.not. bin_require_open) then
         if (present(error)) call fatal_error(error, &
                                              "Cannot "//action//" on a closed binary file")
      end if

   end function bin_require_open

   !> Ensure the file is open and writable before modifying it
   !>
   !> @param[in]  self    Container file instance
   !> @param[in]  action  Description of the attempted operation
   !> @param[out] error   Error handler
   !> @return True when the file may be modified
   logical function bin_require_writable(self, action, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: action
      type(error_type), allocatable, intent(out), optional :: error

      bin_require_writable = .false.
      if (.not. bin_require_open(self, action, error)) return
      bin_require_writable = self%writable
      if (.not. bin_require_writable) then
         if (present(error)) call fatal_error(error, &
                                              "Cannot "//action//" on a read-only binary file")
      end if

   end function bin_require_writable

   !> Expand a possibly relative name into an absolute path
   !>
   !> Relative names are resolved against the currently open group
   !>
   !> @param[in] self Container file instance
   !> @param[in] name Dataset, group or attribute path
   !> @return Absolute path with a leading slash and no trailing slash
   function resolve_path(self, name) result(path)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: name
      character(:), allocatable :: path

      !> Last significant character of the name
      integer :: last

      last = len_trim(name)
      do while (last > 1)
         if (name(last:last) /= "/") exit
         last = last - 1
      end do

      if (last < 1) then
         path = self%prefix
         if (len(path) == 0) path = "/"
         return
      end if

      if (name(1:1) == "/") then
         path = name(1:last)
      else
         path = self%prefix//"/"//name(1:last)
      end if

   end function resolve_path

   !> Locate an index entry by absolute path
   !>
   !> @param[in] self Container file instance
   !> @param[in] path Absolute path to look up
   !> @return Index into `entries`, or 0 when the path is unknown
   integer function find_entry(self, path)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: path

      !> Loop index over entries
      integer :: i

      find_entry = 0
      do i = 1, self%nentries
         if (self%entries(i)%path == path) then
            find_entry = i
            return
         end if
      end do

   end function find_entry

   !> Append an entry to the in-memory index, growing it as needed
   !>
   !> @param[in,out] self Container file instance
   !> @param[in]     ent  Entry to append
   subroutine push_entry(self, ent)
      class(binary_file), intent(inout) :: self
      type(entry_type), intent(in) :: ent

      !> Scratch storage while the index is grown
      type(entry_type), allocatable :: tmp(:)

      if (self%nentries >= size(self%entries)) then
         allocate (tmp(max(2*size(self%entries), 16)))
         tmp(1:self%nentries) = self%entries(1:self%nentries)
         call move_alloc(tmp, self%entries)
      end if

      self%nentries = self%nentries + 1
      self%entries(self%nentries) = ent

   end subroutine push_entry

   !> Open a group within the file for subsequent operations
   !>
   !> Relative dataset names are then resolved inside this group
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     gname Name of the group to open
   !> @param[out]    error Error handler
   subroutine bin_open_group(self, gname, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: gname
      type(error_type), allocatable, intent(out), optional :: error

      !> Absolute path of the group
      character(:), allocatable :: path
      !> Index of the group entry
      integer :: idx

      if (.not. bin_require_open(self, "open a group", error)) return
      if (self%group_open) then
         if (present(error)) then
            call fatal_error(error, &
                             "Cannot open a nested group before closing the current group")
         end if
         return
      end if

      path = resolve_path(self, gname)
      idx = find_entry(self, path)
      if (idx == 0) then
         if (present(error)) call fatal_error(error, "Failed to open group: "//gname)
         return
      end if
      if (self%entries(idx)%dtype /= type_group) then
         if (present(error)) call fatal_error(error, "Not a group: "//gname)
         return
      end if

      self%prefix = path
      self%group_open = .true.

   end subroutine bin_open_group

   !> Close the currently open group and restore the previous location
   !>
   !> @param[in,out] self  Container file instance
   !> @param[out]    error Error handler
   subroutine bin_close_group(self, error)
      class(binary_file), intent(inout) :: self
      type(error_type), allocatable, intent(out), optional :: error

      if (.not. bin_require_open(self, "close a group", error)) return
      if (.not. self%group_open) return

      self%prefix = ""
      self%group_open = .false.

   end subroutine bin_close_group

   !> Check whether a dataset or group exists in the file
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  name  Name of the dataset or group to check
   !> @param[out] exist True if the object exists
   !> @param[out] error Error handler
   subroutine bin_exist(self, name, exist, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: name
      logical, intent(out) :: exist
      type(error_type), allocatable, intent(out), optional :: error

      exist = .false.
      if (.not. bin_require_open(self, "check existence", error)) return

      exist = find_entry(self, resolve_path(self, name)) > 0

   end subroutine bin_exist

   !> Read the extents of a dataset from the index, without its payload
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] dims  Extent along each dimension; empty for scalars and strings
   !> @param[out] error Error handler
   subroutine bin_get_shape(self, dname, dims, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, allocatable, intent(out) :: dims(:)
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the dataset entry
      integer :: idx

      if (.not. bin_require_open(self, "get a dataset shape", error)) return

      idx = find_entry(self, resolve_path(self, dname))
      if (idx == 0) then
         if (present(error)) call fatal_error(error, "Failed to open dataset: "//dname)
         return
      end if

      associate (ent => self%entries(idx))
         if (ent%dtype == type_group) then
            if (present(error)) call fatal_error(error, "Not a dataset: "//dname)
         else if (any(ent%dims > int(huge(0), int64))) then
            if (present(error)) call fatal_error(error, "Dataset extent exceeds default integers: "//dname)
         else
            dims = int(ent%dims(1:ent%rank))
         end if
      end associate

   end subroutine bin_get_shape

   !> Delete a dataset or group from the file
   !>
   !> Deleting a group also removes everything below it.  The payload bytes stay
   !> in the file but become unreachable, since readers only follow the index
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     name  Name of the dataset or group to delete
   !> @param[out]    error Error handler
   subroutine bin_delete(self, name, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: name
      type(error_type), allocatable, intent(out), optional :: error

      !> Absolute path of the entry to remove
      character(:), allocatable :: path
      !> Read and write cursors while compacting the index
      integer :: i, keep

      if (.not. bin_require_writable(self, "delete", error)) return

      path = resolve_path(self, name)
      if (find_entry(self, path) == 0) then
         if (present(error)) call fatal_error(error, "Failed to delete: "//name)
         return
      end if

      keep = 0
      do i = 1, self%nentries
         associate (entry_path => self%entries(i)%path)
            if (entry_path == path) cycle
            if (len(entry_path) > len(path)) then
               if (entry_path(1:len(path)) == path .and. &
                   entry_path(len(path) + 1:len(path) + 1) == "/") cycle
            end if
         end associate
         keep = keep + 1
         if (keep /= i) self%entries(keep) = self%entries(i)
      end do
      self%nentries = keep

   end subroutine bin_delete

   !> Ensure parent groups exist for a given path (internal helper)
   !>
   !> Creates only intermediate groups, not the final path component
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     path  Absolute path (e.g., '/data/results/dataset')
   !> @param[out]    error Error handler
   subroutine bin_ensure_parent_groups(self, path, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: path
      type(error_type), allocatable, intent(out), optional :: error

      !> Position of the current separator
      integer :: sp
      !> Offset of the next separator
      integer :: ep
      !> Index of an existing entry
      integer :: idx

      if (.not. bin_require_writable(self, "create a group", error)) return

      sp = 1
      do
         ep = index(path(sp + 1:), "/")
         if (ep == 0) exit
         sp = sp + ep
         idx = find_entry(self, path(1:sp - 1))
         if (idx > 0) then
            if (self%entries(idx)%dtype /= type_group) then
               if (present(error)) then
                  call fatal_error(error, "Path component is not a group: "//path(1:sp - 1))
               end if
               return
            end if
            cycle
         end if
         call push_entry(self, entry_type(path=path(1:sp - 1), attrs=[attribute_type ::]))
      end do

   end subroutine bin_ensure_parent_groups

   !> Create a group, creating parent groups as needed
   !>
   !> Creating a group that already exists is not an error
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     gname Path of the group to create (e.g., '/data/results')
   !> @param[out]    error Error handler
   subroutine bin_add_group(self, gname, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: gname
      type(error_type), allocatable, intent(out), optional :: error

      !> Absolute path of the group
      character(:), allocatable :: path
      !> Index of an existing entry
      integer :: idx

      if (.not. bin_require_writable(self, "add a group", error)) return

      path = resolve_path(self, gname)
      call self%bin_ensure_parent_groups(path, error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      if (path == "/") return

      idx = find_entry(self, path)
      if (idx > 0) then
         if (self%entries(idx)%dtype /= type_group) then
            if (present(error)) call fatal_error(error, "Path already exists as a dataset: "//gname)
         end if
         return
      end if

      call push_entry(self, entry_type(path=path, attrs=[attribute_type ::]))

   end subroutine bin_add_group

   !> Reserve an index slot and stream position for a new dataset
   !>
   !> @param[in,out] self   Container file instance
   !> @param[in]     dname  Dataset name as given by the caller
   !> @param[in]     dtype  Payload type code
   !> @param[in]     dims   Extent along each dimension, in Fortran order
   !> @param[in]     nbytes Payload size in bytes
   !> @param[out]    pos    Stream position to write the payload at
   !> @param[out]    ok     True when the dataset may be written
   !> @param[out]    error  Error handler
   subroutine begin_dataset(self, dname, dtype, dims, nbytes, pos, ok, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      integer(int32), intent(in) :: dtype
      integer(int64), intent(in) :: dims(:)
      integer(int64), intent(in) :: nbytes
      integer(int64), intent(out) :: pos
      logical, intent(out) :: ok
      type(error_type), allocatable, intent(out), optional :: error

      !> Entry describing the new dataset
      type(entry_type) :: ent
      !> Absolute path of the dataset
      character(:), allocatable :: path

      ok = .false.
      pos = 0_int64
      if (.not. bin_require_writable(self, "add a dataset", error)) return

      path = resolve_path(self, dname)
      call self%bin_ensure_parent_groups(path, error)
      if (present(error)) then
         if (allocated(error)) return
      end if

      if (find_entry(self, path) > 0) then
         if (present(error)) then
            call fatal_error(error, "Failed to create dataset, path exists: "//dname)
         end if
         return
      end if

      pos = self%next_pos

      ent%path = path
      ent%dtype = dtype
      ent%rank = int(size(dims), int32)
      ent%dims = 0_int64
      ent%dims(1:size(dims)) = dims
      ent%offset = pos
      ent%nbytes = nbytes
      allocate (ent%attrs(0))
      call push_entry(self, ent)

      ok = .true.

   end subroutine begin_dataset

   !> Commit a dataset payload after it has been written
   !>
   !> On a failed transfer the reserved index entry is dropped again
   !>
   !> @param[in,out] self   Container file instance
   !> @param[in]     dname  Dataset name as given by the caller
   !> @param[in]     nbytes Payload size in bytes
   !> @param[in]     stat   I/O status of the payload write
   !> @param[out]    error  Error handler
   subroutine end_dataset(self, dname, nbytes, stat, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      integer(int64), intent(in) :: nbytes
      integer, intent(in) :: stat
      type(error_type), allocatable, intent(out), optional :: error

      if (stat /= 0) then
         self%nentries = self%nentries - 1
         if (present(error)) call fatal_error(error, "Failed to write dataset: "//dname)
         return
      end if

      self%next_pos = self%next_pos + nbytes

   end subroutine end_dataset

   !> Locate a dataset and validate its type and rank before reading
   !>
   !> @param[in]  self          Container file instance
   !> @param[in]  dname         Dataset name as given by the caller
   !> @param[in]  dtype         Expected payload type code
   !> @param[in]  expected_rank Expected dataset rank
   !> @param[out] idx           Index of the dataset entry, 0 on failure
   !> @param[out] error         Error handler
   subroutine locate_dataset(self, dname, dtype, expected_rank, idx, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer(int32), intent(in) :: dtype
      integer, intent(in) :: expected_rank
      integer, intent(out) :: idx
      type(error_type), allocatable, intent(out), optional :: error

      idx = 0
      if (.not. bin_require_open(self, "get a dataset", error)) return

      idx = find_entry(self, resolve_path(self, dname))
      if (idx == 0) then
         if (present(error)) call fatal_error(error, "Failed to open dataset: "//dname)
         return
      end if

      if (self%entries(idx)%dtype /= dtype) then
         idx = 0
         if (present(error)) call fatal_error(error, "Dataset type mismatch for: "//dname)
         return
      end if

      if (int(self%entries(idx)%rank) /= expected_rank) then
         idx = 0
         if (present(error)) call fatal_error(error, "Dataset rank mismatch for: "//dname)
         return
      end if

   end subroutine locate_dataset

   !> Add a string attribute to a group or dataset
   !>
   !> An existing attribute of the same name is replaced
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     path  Path to the group or dataset
   !> @param[in]     name  Attribute name
   !> @param[in]     value Attribute value
   !> @param[out]    error Error handler
   subroutine bin_adda_string(self, path, name, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: path
      character(*), intent(in) :: name
      character(*), intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Scratch storage while the attribute list is grown
      type(attribute_type), allocatable :: tmp(:)
      !> Index of the target entry
      integer :: idx
      !> Loop index over attributes
      integer :: i
      !> Number of attributes already present
      integer :: nattr

      if (.not. bin_require_writable(self, "add an attribute", error)) return

      idx = find_entry(self, resolve_path(self, path))
      if (idx == 0) then
         if (present(error)) call fatal_error(error, &
                                              "Failed to set attribute on missing path: "//path)
         return
      end if

      nattr = size(self%entries(idx)%attrs)
      do i = 1, nattr
         if (self%entries(idx)%attrs(i)%name == name) then
            self%entries(idx)%attrs(i)%value = value
            return
         end if
      end do

      allocate (tmp(nattr + 1))
      tmp(1:nattr) = self%entries(idx)%attrs
      tmp(nattr + 1)%name = name
      tmp(nattr + 1)%value = value
      call move_alloc(tmp, self%entries(idx)%attrs)

   end subroutine bin_adda_string

   !> Get a string attribute from a group or dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  path  Path to the group or dataset
   !> @param[in]  name  Attribute name
   !> @param[out] value Allocatable string to receive the attribute value
   !> @param[out] error Error handler
   subroutine bin_geta_string(self, path, name, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: path
      character(*), intent(in) :: name
      character(:), intent(out), allocatable :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the target entry
      integer :: idx
      !> Loop index over attributes
      integer :: i

      if (.not. bin_require_open(self, "get an attribute", error)) return

      idx = find_entry(self, resolve_path(self, path))
      if (idx > 0) then
         do i = 1, size(self%entries(idx)%attrs)
            if (self%entries(idx)%attrs(i)%name == name) then
               value = self%entries(idx)%attrs(i)%value
               return
            end if
         end do
      end if

      if (present(error)) call fatal_error(error, &
                                           "Failed to get attribute: "//name//" in "//path)

   end subroutine bin_geta_string

   !> Write a scalar integer dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value Scalar integer value
   !> @param[out]    error Error handler
   subroutine bin_add_int(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = int_bytes
      call begin_dataset(self, dname, type_int, [integer(int64) ::], nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) int(value, int32)
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_int

   !> Write a 1D integer array dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value 1D integer array
   !> @param[out]    error Error handler
   subroutine bin_add_int1d(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = int_bytes*size(value, kind=int64)
      call begin_dataset(self, dname, type_int, shape(value, kind=int64), nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) int(value, int32)
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_int1d

   !> Write a 2D integer array dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value 2D integer array
   !> @param[out]    error Error handler
   subroutine bin_add_int2d(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = int_bytes*size(value, kind=int64)
      call begin_dataset(self, dname, type_int, shape(value, kind=int64), nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) int(value, int32)
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_int2d

   !> Write a 3D integer array dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value 3D integer array
   !> @param[out]    error Error handler
   subroutine bin_add_int3d(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      integer, intent(in) :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = int_bytes*size(value, kind=int64)
      call begin_dataset(self, dname, type_int, shape(value, kind=int64), nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) int(value, int32)
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_int3d

   !> Write a scalar real dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value Scalar real value
   !> @param[out]    error Error handler
   subroutine bin_add_real(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = real_bytes
      call begin_dataset(self, dname, type_real, [integer(int64) ::], nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) value
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_real

   !> Write a 1D real array dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value 1D real array
   !> @param[out]    error Error handler
   subroutine bin_add_real1d(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = real_bytes*size(value, kind=int64)
      call begin_dataset(self, dname, type_real, shape(value, kind=int64), nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) value
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_real1d

   !> Write a 2D real array dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value 2D real array
   !> @param[out]    error Error handler
   subroutine bin_add_real2d(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = real_bytes*size(value, kind=int64)
      call begin_dataset(self, dname, type_real, shape(value, kind=int64), nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) value
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_real2d

   !> Write a 3D real array dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value 3D real array
   !> @param[out]    error Error handler
   subroutine bin_add_real3d(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      real(wp), intent(in) :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = real_bytes*size(value, kind=int64)
      call begin_dataset(self, dname, type_real, shape(value, kind=int64), nbytes, pos, ok, error)
      if (.not. ok) return

      write (self%unit, pos=pos, iostat=stat) value
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_real3d

   !> Write a string dataset
   !>
   !> @param[in,out] self  Container file instance
   !> @param[in]     dname Dataset name
   !> @param[in]     value String value
   !> @param[out]    error Error handler
   subroutine bin_add_string(self, dname, value, error)
      class(binary_file), intent(inout) :: self
      character(*), intent(in) :: dname
      character(*), intent(in) :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Stream position of the payload
      integer(int64) :: pos
      !> Payload size in bytes
      integer(int64) :: nbytes
      !> True when the slot was reserved
      logical :: ok
      !> I/O status code
      integer :: stat

      nbytes = len(value, kind=int64)
      call begin_dataset(self, dname, type_string, [integer(int64) ::], nbytes, pos, ok, error)
      if (.not. ok) return

      stat = 0
      if (nbytes > 0_int64) write (self%unit, pos=pos, iostat=stat) value
      call end_dataset(self, dname, nbytes, stat, error)

   end subroutine bin_add_string

   !> Read a string dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable string to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_string(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      character(:), intent(out), allocatable :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_string, 0, idx, error)
      if (idx == 0) return

      allocate (character(int(self%entries(idx)%nbytes)) :: value, stat=stat)
      if (stat == 0 .and. self%entries(idx)%nbytes > 0_int64) then
         read (self%unit, pos=self%entries(idx)%offset, iostat=stat) value
      end if
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
      end if

   end subroutine bin_get_string

   !> Read a scalar integer dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Scalar integer to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_int(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out) :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Storage matching the on-disk element width
      integer(int32) :: buffer
      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      value = 0
      call locate_dataset(self, dname, type_int, 0, idx, error)
      if (idx == 0) return

      read (self%unit, pos=self%entries(idx)%offset, iostat=stat) buffer
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
      value = int(buffer)

   end subroutine bin_get_int

   !> Read a 1D integer array dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 1D integer array to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_int1d(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out), allocatable :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

      !> Storage matching the on-disk element width
      integer(int32), allocatable :: buffer(:)
      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_int, 1, idx, error)
      if (idx == 0) return

      associate (dims => self%entries(idx)%dims)
         allocate (buffer(dims(1)), stat=stat)
      end associate
      if (stat == 0) read (self%unit, pos=self%entries(idx)%offset, iostat=stat) buffer
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
      value = int(buffer)

   end subroutine bin_get_int1d

   !> Read a 2D integer array dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 2D integer array to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_int2d(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out), allocatable :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Storage matching the on-disk element width
      integer(int32), allocatable :: buffer(:, :)
      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_int, 2, idx, error)
      if (idx == 0) return

      associate (dims => self%entries(idx)%dims)
         allocate (buffer(dims(1), dims(2)), stat=stat)
      end associate
      if (stat == 0) read (self%unit, pos=self%entries(idx)%offset, iostat=stat) buffer
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
      value = int(buffer)

   end subroutine bin_get_int2d

   !> Read a 3D integer array dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 3D integer array to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_int3d(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      integer, intent(out), allocatable :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Storage matching the on-disk element width
      integer(int32), allocatable :: buffer(:, :, :)
      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_int, 3, idx, error)
      if (idx == 0) return

      associate (dims => self%entries(idx)%dims)
         allocate (buffer(dims(1), dims(2), dims(3)), stat=stat)
      end associate
      if (stat == 0) read (self%unit, pos=self%entries(idx)%offset, iostat=stat) buffer
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
         return
      end if
      value = int(buffer)

   end subroutine bin_get_int3d

   !> Read a scalar real dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Scalar real to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_real(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out) :: value
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      value = 0.0_wp
      call locate_dataset(self, dname, type_real, 0, idx, error)
      if (idx == 0) return

      read (self%unit, pos=self%entries(idx)%offset, iostat=stat) value
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
      end if

   end subroutine bin_get_real

   !> Read a 1D real array dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 1D real array to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_real1d(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out), allocatable :: value(:)
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_real, 1, idx, error)
      if (idx == 0) return

      associate (dims => self%entries(idx)%dims)
         allocate (value(dims(1)), stat=stat)
      end associate
      if (stat == 0) read (self%unit, pos=self%entries(idx)%offset, iostat=stat) value
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
      end if

   end subroutine bin_get_real1d

   !> Read a 2D real array dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 2D real array to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_real2d(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out), allocatable :: value(:, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_real, 2, idx, error)
      if (idx == 0) return

      associate (dims => self%entries(idx)%dims)
         allocate (value(dims(1), dims(2)), stat=stat)
      end associate
      if (stat == 0) read (self%unit, pos=self%entries(idx)%offset, iostat=stat) value
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
      end if

   end subroutine bin_get_real2d

   !> Read a 3D real array dataset
   !>
   !> @param[in]  self  Container file instance
   !> @param[in]  dname Dataset name
   !> @param[out] value Allocatable 3D real array to receive the data
   !> @param[out] error Error handler
   subroutine bin_get_real3d(self, dname, value, error)
      class(binary_file), intent(in) :: self
      character(*), intent(in) :: dname
      real(wp), intent(out), allocatable :: value(:, :, :)
      type(error_type), allocatable, intent(out), optional :: error

      !> Index of the dataset entry
      integer :: idx
      !> I/O status code
      integer :: stat

      call locate_dataset(self, dname, type_real, 3, idx, error)
      if (idx == 0) return

      associate (dims => self%entries(idx)%dims)
         allocate (value(dims(1), dims(2), dims(3)), stat=stat)
      end associate
      if (stat == 0) read (self%unit, pos=self%entries(idx)%offset, iostat=stat) value
      if (stat /= 0) then
         if (present(error)) call fatal_error(error, "Failed to read dataset: "//dname)
      end if

   end subroutine bin_get_real3d

end module moist_utils_binio
