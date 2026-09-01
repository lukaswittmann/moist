!> Self-describing table of the per-point and per-sphere results a cavity holds
module moist_cavity_fields
   use mctc_env, only: wp, error_type, fatal_error

   implicit none (type, external)
   private

   public :: cavity_field_info_type, cavity_field_query_type
   public :: cavity_field_real, cavity_field_int, cavity_field_bool
   public :: cavity_field_max_rank

   !> Element type tags. Mirrored by the `MOIST_FIELD_*` macros in moist.h;
   !> the values are part of the public C contract and must not be renumbered.
   integer, parameter :: cavity_field_real = 1
   integer, parameter :: cavity_field_int = 2
   integer, parameter :: cavity_field_bool = 3

   !> Highest array rank a field can have
   integer, parameter :: cavity_field_max_rank = 2

   !> Shape and type of one readable cavity field
   type :: cavity_field_info_type
      !> Lookup key, always the cavity's own name for the array
      character(len=:), allocatable :: name
      !> One-line description of what the field holds
      character(len=:), allocatable :: about
      !> One of `cavity_field_real`, `cavity_field_int`, `cavity_field_bool`
      integer :: dtype = cavity_field_real
      !> Array rank; 0 marks a scalar
      integer :: rank = 0
      !> Extent of each dimension, fastest-varying first; unused entries are 1
      integer :: dims(cavity_field_max_rank) = 1
   contains
      !> Number of elements the field writes
      procedure :: count => field_element_count
   end type cavity_field_info_type

   !> Walker handed to `cavity_type%list_fields`
   !>
   !> In enumeration mode it collects a descriptor per declared field.  In fetch
   !> mode it copies out the one field whose name matches `want` and ignores the
   !> rest, so listing stays cheap and reading copies exactly one array.
   type :: cavity_field_query_type
      !> Name being fetched; unallocated while enumerating
      character(len=:), allocatable :: want
      !> Descriptors of every declared field (enumeration mode)
      type(cavity_field_info_type), allocatable :: info(:)
      !> Number of entries filled in `info`
      integer :: nfield = 0
      !> Whether the fetched name was declared by the cavity
      logical :: found = .false.
      !> Shape and type of the fetched field
      type(cavity_field_info_type) :: hit
      !> Fetched payload, flattened in Fortran order; one of the three is
      !> allocated according to `hit%dtype`
      real(wp), allocatable :: rvals(:)
      integer, allocatable :: ivals(:)
      logical, allocatable :: lvals(:)
   contains
      !> Switch to enumeration mode and discard any earlier result
      procedure :: enumerate => query_enumerate
      !> Switch to fetch mode for one field name
      procedure :: fetch => query_fetch
      !> Index of a name in `info`, or 0 when it was not declared
      procedure :: index_of => query_index_of

      !> Declare a per-point or per-sphere real array
      procedure :: add_real
      !> Declare a rank-2 real array, e.g. (3, ngrid)
      procedure :: add_real2
      !> Declare an integer array
      procedure :: add_int
      !> Declare a logical array
      procedure :: add_bool
      !> Declare an allocatable real scalar
      procedure :: add_real_scalar
      !> Declare an allocatable integer scalar
      procedure :: add_int_scalar
      !> Declare a plain integer scalar that is always present
      procedure :: add_int_value

      !> Record a declaration; reports whether the payload is wanted
      procedure, private :: record
   end type cavity_field_query_type

contains

   !* ================================================================================= *!
   !*                                    Descriptors                                    *!
   !* ================================================================================= *!

   !> Number of elements a field writes
   !>
   !> @param[in] self  Field descriptor
   !> @return          Product of the declared extents; 1 for a scalar
   pure integer function field_element_count(self) result(count)
      class(cavity_field_info_type), intent(in) :: self
      integer :: idim

      count = 1
      do idim = 1, self%rank
         count = count*self%dims(idim)
      end do

   end function field_element_count

   !* ================================================================================= *!
   !*                                     Query modes                                   *!
   !* ================================================================================= *!

   !> Prepare the query to collect a descriptor for every declared field
   !>
   !> @param[inout] self  Query walker
   subroutine query_enumerate(self)
      class(cavity_field_query_type), intent(inout) :: self

      call reset_query(self)
      if (allocated(self%info)) deallocate (self%info)
      allocate (self%info(0))
      self%nfield = 0

   end subroutine query_enumerate

   !> Prepare the query to copy out a single named field
   !>
   !> @param[inout] self  Query walker
   !> @param[in]    name  Field to fetch
   subroutine query_fetch(self, name)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field to fetch
      character(len=*), intent(in) :: name

      call reset_query(self)
      self%want = name

   end subroutine query_fetch

   !> Clear everything a previous pass left behind
   !>
   !> @param[inout] self  Query walker
   subroutine reset_query(self)
      class(cavity_field_query_type), intent(inout) :: self

      if (allocated(self%want)) deallocate (self%want)
      if (allocated(self%rvals)) deallocate (self%rvals)
      if (allocated(self%ivals)) deallocate (self%ivals)
      if (allocated(self%lvals)) deallocate (self%lvals)
      self%found = .false.

   end subroutine reset_query

   !> Position of a field in the enumerated list
   !>
   !> @param[in] self  Query walker holding an enumeration result
   !> @param[in] name  Field to look for
   !> @return          1-based index, or 0 when the cavity does not have it
   pure integer function query_index_of(self, name) result(pos)
      class(cavity_field_query_type), intent(in) :: self
      !> Field to look for
      character(len=*), intent(in) :: name
      integer :: ifield

      pos = 0
      if (.not. allocated(self%info)) return
      do ifield = 1, self%nfield
         if (self%info(ifield)%name == name) then
            pos = ifield
            return
         end if
      end do

   end function query_index_of

   !* ================================================================================= *!
   !*                                 Field declaration                                 *!
   !* ================================================================================= *!

   !> Record one declared field and report whether its payload is wanted
   !>
   !> @param[inout] self   Query walker
   !> @param[in]    name   Field name
   !> @param[in]    about  One-line description
   !> @param[in]    dtype  Element type tag
   !> @param[in]    rank   Array rank, 0 for a scalar
   !> @param[in]    dims   Extents, fastest-varying first
   !> @return              True when the caller should copy the payload
   logical function record(self, name, about, dtype, rank, dims) result(wanted)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Element type tag
      integer, intent(in) :: dtype
      !> Array rank, 0 for a scalar
      integer, intent(in) :: rank
      !> Extents, fastest-varying first
      integer, intent(in) :: dims(:)

      type(cavity_field_info_type) :: info
      type(cavity_field_info_type), allocatable :: grown(:)

      wanted = .false.

      !> Fetching only ever wants one name, so reject the others before paying
      !> for a descriptor.
      if (allocated(self%want)) then
         if (name /= self%want) return
      end if

      info%name = name
      info%about = about
      info%dtype = dtype
      info%rank = rank
      info%dims = 1
      if (rank > 0) info%dims(:rank) = dims(:rank)

      if (allocated(self%want)) then
         self%hit = info
         self%found = .true.
         wanted = .true.
         return
      end if

      if (.not. allocated(self%info)) allocate (self%info(0))
      if (self%nfield >= size(self%info)) then
         allocate (grown(max(16, 2*size(self%info))))
         grown(:self%nfield) = self%info(:self%nfield)
         call move_alloc(grown, self%info)
      end if
      self%nfield = self%nfield + 1
      self%info(self%nfield) = info

   end function record

   !> Declare a real array indexed by grid point or sphere
   !>
   !> @param[inout] self    Query walker
   !> @param[in]    name    Field name
   !> @param[in]    about   One-line description
   !> @param[in]    values  Backing array; an unallocated array is not declared
   subroutine add_real(self, name, about, values)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Backing array
      real(wp), allocatable, intent(in) :: values(:)

      if (.not. allocated(values)) return
      if (.not. self%record(name, about, cavity_field_real, 1, [size(values)])) return
      self%rvals = values

   end subroutine add_real

   !> Declare a rank-2 real array, stored and handed out in Fortran order
   !>
   !> @param[inout] self    Query walker
   !> @param[in]    name    Field name
   !> @param[in]    about   One-line description
   !> @param[in]    values  Backing array; an unallocated array is not declared
   subroutine add_real2(self, name, about, values)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Backing array
      real(wp), allocatable, intent(in) :: values(:, :)

      if (.not. allocated(values)) return
      if (.not. self%record(name, about, cavity_field_real, 2, shape(values))) return
      self%rvals = reshape(values, [size(values)])

   end subroutine add_real2

   !> Declare an integer array
   !>
   !> `zero_based` marks a Fortran 1-based index that the C API hands out
   !> 0-based; the shift is applied here so every consumer sees one convention.
   !>
   !> @param[inout] self        Query walker
   !> @param[in]    name        Field name
   !> @param[in]    about       One-line description
   !> @param[in]    values      Backing array; an unallocated array is not declared
   !> @param[in]    zero_based  Subtract one on read (default false)
   subroutine add_int(self, name, about, values, zero_based)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Backing array
      integer, allocatable, intent(in) :: values(:)
      !> Subtract one on read
      logical, intent(in), optional :: zero_based

      logical :: shift

      shift = .false.
      if (present(zero_based)) shift = zero_based

      if (.not. allocated(values)) return
      if (.not. self%record(name, about, cavity_field_int, 1, [size(values)])) return
      if (shift) then
         self%ivals = values - 1
      else
         self%ivals = values
      end if

   end subroutine add_int

   !> Declare a logical array
   !>
   !> @param[inout] self    Query walker
   !> @param[in]    name    Field name
   !> @param[in]    about   One-line description
   !> @param[in]    values  Backing array; an unallocated array is not declared
   subroutine add_bool(self, name, about, values)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Backing array
      logical, allocatable, intent(in) :: values(:)

      if (.not. allocated(values)) return
      if (.not. self%record(name, about, cavity_field_bool, 1, [size(values)])) return
      self%lvals = values

   end subroutine add_bool

   !> Declare an allocatable real scalar
   !>
   !> @param[inout] self   Query walker
   !> @param[in]    name   Field name
   !> @param[in]    about  One-line description
   !> @param[in]    value  Backing scalar; an unallocated scalar is not declared
   subroutine add_real_scalar(self, name, about, value)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Backing scalar
      real(wp), allocatable, intent(in) :: value

      if (.not. allocated(value)) return
      if (.not. self%record(name, about, cavity_field_real, 0, [integer ::])) return
      self%rvals = [value]

   end subroutine add_real_scalar

   !> Declare an allocatable integer scalar
   !>
   !> @param[inout] self   Query walker
   !> @param[in]    name   Field name
   !> @param[in]    about  One-line description
   !> @param[in]    value  Backing scalar; an unallocated scalar is not declared
   subroutine add_int_scalar(self, name, about, value)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Backing scalar
      integer, allocatable, intent(in) :: value

      if (.not. allocated(value)) return
      if (.not. self%record(name, about, cavity_field_int, 0, [integer ::])) return
      self%ivals = [value]

   end subroutine add_int_scalar

   !> Declare a plain integer scalar that the cavity always carries
   !>
   !> @param[inout] self   Query walker
   !> @param[in]    name   Field name
   !> @param[in]    about  One-line description
   !> @param[in]    value  Value to report
   subroutine add_int_value(self, name, about, value)
      class(cavity_field_query_type), intent(inout) :: self
      !> Field name
      character(len=*), intent(in) :: name
      !> One-line description
      character(len=*), intent(in) :: about
      !> Value to report
      integer, intent(in) :: value

      if (.not. self%record(name, about, cavity_field_int, 0, [integer ::])) return
      self%ivals = [value]

   end subroutine add_int_value

end module moist_cavity_fields
