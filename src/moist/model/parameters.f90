!> Shared JSON/TOML file input and output, with JSON printing for native
!> parameter objects
!>
!> Derived types declare their fields in register_entries and defaults in
!> init_defaults
module moist_model_parameters
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use jonquil, only: json_load, json_dump
   use tomlf, only: toml_table, toml_value, toml_array, toml_error, &
      toml_load, toml_dump, get_value, set_value, len
   use tomlf_type, only: cast_to_table
   implicit none(type, external)
   private
   public :: moist_model_parameters_type

   !> Base for parameter values; no references to caller-owned fields are retained
   type, abstract :: moist_model_parameters_type
      private
      !> Temporary document, present only during a file or print operation
      class(toml_value), allocatable :: document
      !> Whether registered fields are being read from the document
      logical :: reading = .false.
      !> First field-conversion error during the current operation
      type(error_type), allocatable :: failure
   contains
      !> Declare the parameter fields through register_* calls
      procedure(register_entries_ifc), deferred :: register_entries
      !> Restore compiled defaults
      procedure(init_defaults_ifc), deferred :: init_defaults
      !> Validate values and recompute derived settings
      procedure :: validate => validate_parameters
      !> Read JSON or TOML with defaults for omitted fields
      procedure :: read_file
      !> Write JSON or TOML, selected by the file extension
      procedure :: write_file
      !> Print parameter values as formatted JSON
      procedure :: print_parameters
      procedure :: register_real_scalar
      procedure :: register_int_scalar
      procedure :: register_logical
      procedure :: register_real_vector
      procedure :: register_string
      procedure :: register_alloc_string
      procedure, private :: field_parent
      procedure, private :: collect_document
      procedure, private :: clear_document
   end type moist_model_parameters_type

   abstract interface
      subroutine init_defaults_ifc(self)
         import :: moist_model_parameters_type
         implicit none(type, external)
         class(moist_model_parameters_type), intent(inout) :: self
      end subroutine init_defaults_ifc
      subroutine register_entries_ifc(self)
         import :: moist_model_parameters_type
         implicit none(type, external)
         class(moist_model_parameters_type), intent(inout), target :: self
      end subroutine register_entries_ifc
   end interface

contains

   !> Select the parser or serializer from a case-insensitive file extension
   !>
   !> @param[in] filepath Input or output path
   !> @param[out] use_toml True for TOML, false for JSON
   !> @param[out] error Unsupported or missing extension
   subroutine resolve_file_format(filepath, use_toml, error)
      !> Parameter file path
      character(len=*), intent(in) :: filepath
      !> Selected format
      logical, intent(out) :: use_toml
      !> Format-selection error
      type(error_type), allocatable, intent(out) :: error
      !> Last extension separator
      integer :: dot

      use_toml = .false.
      dot = index(trim(filepath), ".", back=.true.)
      if (dot > 0) then
         select case (to_lower(trim(filepath(dot:))))
         case (".json")
            return
         case (".toml")
            use_toml = .true.
            return
         end select
      end if
      call fatal_error(error, "Parameter file extension must be .json or .toml: "//trim(filepath))
   end subroutine resolve_file_format

   !> Default validation for parameter sets with no derived values
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error Validation error
   subroutine validate_parameters(self, error)
      class(moist_model_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
   end subroutine validate_parameters

   !> Discard the temporary document after each operation
   !>
   !> @param[inout] self Parameter values
   subroutine clear_document(self)
      class(moist_model_parameters_type), intent(inout) :: self
      if (allocated(self%document)) deallocate(self%document)
      if (allocated(self%failure)) deallocate(self%failure)
      self%reading = .false.
   end subroutine clear_document

   !> Read JSON or TOML, apply defaults for missing fields, then validate
   !>
   !> @param[inout] self Parameter values
   !> @param[in] filepath Input path ending in .json or .toml
   !> @param[out] error File, conversion, or validation error
   subroutine read_file(self, filepath, error)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: filepath
      type(error_type), allocatable, intent(out) :: error
      type(toml_error), allocatable :: format_error
      type(toml_table), pointer :: root
      !> TOML parser result; moved into the shared document
      type(toml_table), allocatable :: table
      !> Selected input format
      logical :: use_toml

      call self%clear_document()
      call resolve_file_format(filepath, use_toml, error)
      if (allocated(error)) return
      if (use_toml) then
         call toml_load(table, filepath, error=format_error)
         call move_alloc(table, self%document)
      else
         call json_load(self%document, filepath, error=format_error)
      end if
      if (allocated(format_error)) then
         call fatal_error(error, "Cannot read parameter file: "//format_error%message)
      else
         root => cast_to_table(self%document)
         if (.not. associated(root)) then
            call fatal_error(error, "Parameter file must contain an object or table")
         else
            call self%init_defaults()
            self%reading = .true.
            call self%register_entries()
            call move_alloc(self%failure, error)
         end if
      end if
      call self%clear_document()
      if (.not. allocated(error)) call self%validate(error)
   end subroutine read_file

   !> Collect parameter values into a temporary document
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error Allocation or conversion error
   subroutine collect_document(self, error)
      class(moist_model_parameters_type), intent(inout), target :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: stat

      call self%clear_document()
      allocate(toml_table :: self%document, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate parameter document")
         return
      end if
      call self%register_entries()
      call move_alloc(self%failure, error)
   end subroutine collect_document

   !> Write JSON or TOML, replacing an existing file
   !>
   !> @param[inout] self Parameter values
   !> @param[in] filepath Output path ending in .json or .toml
   !> @param[out] error File or serialization error
   subroutine write_file(self, filepath, error)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: filepath
      type(error_type), allocatable, intent(out) :: error
      type(toml_error), allocatable :: format_error
      integer :: unit, stat
      character(len=512) :: message
      !> Selected output format
      logical :: use_toml

      call self%clear_document()
      call resolve_file_format(filepath, use_toml, error)
      if (allocated(error)) return
      call self%collect_document(error)
      if (.not. allocated(error)) then
         open(newunit=unit, file=filepath, status="replace", action="write", iostat=stat, iomsg=message)
         if (stat /= 0) then
            call fatal_error(error, "Cannot write parameter file: "//trim(message))
         else
            if (use_toml) then
               call toml_dump(self%document, unit, format_error)
            else
               call json_dump(self%document, unit, format_error)
            end if
            if (allocated(format_error)) call fatal_error(error, format_error%message)
            close(unit, iostat=stat, iomsg=message)
            if (stat /= 0 .and. .not. allocated(error)) call fatal_error(error, trim(message))
         end if
      end if
      call self%clear_document()
   end subroutine write_file

   !> Print parameter values as formatted JSON
   !>
   !> @param[inout] self Parameter values
   !> @param[out] error Serialization or output error
   !> @param[in] unit Output unit; defaults to standard output
   subroutine print_parameters(self, error, unit)
      class(moist_model_parameters_type), intent(inout), target :: self
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: unit
      type(toml_error), allocatable :: format_error
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      call self%collect_document(error)
      if (.not. allocated(error)) then
         call json_dump(self%document, output, format_error)
         if (allocated(format_error)) call fatal_error(error, format_error%message)
      end if
      call self%clear_document()
   end subroutine print_parameters

   !> Resolve a dotted key; missing input fields keep their default values
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[out] parent Parent object, or null when no field is available
   !> @param[out] leaf Final key segment
   subroutine field_parent(self, key, parent, leaf)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      type(toml_table), pointer, intent(out) :: parent
      character(len=:), allocatable, intent(out) :: leaf
      type(toml_table), pointer :: child
      class(toml_value), pointer :: value
      integer :: first, dot, stat

      nullify(parent)
      leaf = key
      if (allocated(self%failure) .or. .not. allocated(self%document)) return
      parent => cast_to_table(self%document)
      first = 1
      do
         dot = index(key(first:), ".")
         if (dot == 0) exit
         dot = dot + first - 1
         if (self%reading) then
            call parent%get(key(first:dot-1), value)
            if (.not. associated(value)) then
               nullify(parent)
               return
            end if
         end if
         call get_value(parent, key(first:dot-1), child, requested=.not. self%reading, stat=stat)
         if (stat /= 0 .or. .not. associated(child)) then
            call fatal_error(self%failure, "Invalid parameter group: "//key)
            nullify(parent)
            return
         end if
         parent => child
         first = dot + 1
      end do
      leaf = key(first:)
      if (self%reading) then
         call parent%get(leaf, value)
         if (.not. associated(value)) nullify(parent)
      end if
   end subroutine field_parent

   !> Read or write one real scalar field
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[inout] value Parameter value
   subroutine register_real_scalar(self, key, value)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      real(wp), intent(inout) :: value
      type(toml_table), pointer :: parent
      character(len=:), allocatable :: leaf
      real(wp) :: parsed
      integer :: stat

      call self%field_parent(key, parent, leaf)
      if (.not. associated(parent)) return
      if (self%reading) then
         call get_value(parent, leaf, parsed, stat=stat)
         if (stat == 0) value = parsed
      else
         call set_value(parent, leaf, value, stat=stat)
      end if
      if (stat /= 0) call fatal_error(self%failure, "Invalid parameter: "//key)
   end subroutine register_real_scalar

   !> Read or write one int scalar field
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[inout] value Parameter value
   subroutine register_int_scalar(self, key, value)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      integer, intent(inout) :: value
      type(toml_table), pointer :: parent
      character(len=:), allocatable :: leaf
      integer :: parsed
      integer :: stat

      call self%field_parent(key, parent, leaf)
      if (.not. associated(parent)) return
      if (self%reading) then
         call get_value(parent, leaf, parsed, stat=stat)
         if (stat == 0) value = parsed
      else
         call set_value(parent, leaf, value, stat=stat)
      end if
      if (stat /= 0) call fatal_error(self%failure, "Invalid parameter: "//key)
   end subroutine register_int_scalar

   !> Read or write one logical field
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[inout] value Parameter value
   subroutine register_logical(self, key, value)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      logical, intent(inout) :: value
      type(toml_table), pointer :: parent
      character(len=:), allocatable :: leaf
      logical :: parsed
      integer :: stat

      call self%field_parent(key, parent, leaf)
      if (.not. associated(parent)) return
      if (self%reading) then
         call get_value(parent, leaf, parsed, stat=stat)
         if (stat == 0) value = parsed
      else
         call set_value(parent, leaf, value, stat=stat)
      end if
      if (stat /= 0) call fatal_error(self%failure, "Invalid parameter: "//key)
   end subroutine register_logical

   !> Read or write a fixed-size real vector
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[inout] values Parameter values
   subroutine register_real_vector(self, key, values)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      real(wp), intent(inout) :: values(:)
      type(toml_table), pointer :: parent
      type(toml_array), pointer :: array
      character(len=:), allocatable :: leaf
      real(wp) :: parsed
      integer :: stat, i

      call self%field_parent(key, parent, leaf)
      if (.not. associated(parent)) return
      call get_value(parent, leaf, array, requested=.not. self%reading, stat=stat)
      if (stat == 0 .and. associated(array)) then
         if (self%reading) then
            if (len(array) /= size(values)) stat = 1
         end if
         if (stat == 0) then
            do i = 1, size(values)
               if (self%reading) then
                  call get_value(array, i, parsed, stat=stat)
                  if (stat == 0) values(i) = parsed
               else
                  call set_value(array, i, values(i), stat=stat)
               end if
               if (stat /= 0) exit
            end do
         end if
      else
         stat = 1
      end if
      if (stat /= 0) call fatal_error(self%failure, "Invalid parameter vector: "//key)
   end subroutine register_real_vector

   !> Read or write an allocatable string; absent values remain unallocated
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[inout] value Parameter value
   subroutine register_alloc_string(self, key, value)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      character(len=:), allocatable, intent(inout) :: value
      type(toml_table), pointer :: parent
      character(len=:), allocatable :: leaf, parsed
      integer :: stat

      if (.not. self%reading .and. .not. allocated(value)) return
      call self%field_parent(key, parent, leaf)
      if (.not. associated(parent)) return
      if (self%reading) then
         call get_value(parent, leaf, parsed, stat=stat)
         if (stat == 0) call move_alloc(parsed, value)
      else
         call set_value(parent, leaf, value, stat=stat)
      end if
      if (stat /= 0) call fatal_error(self%failure, "Invalid parameter string: "//key)
   end subroutine register_alloc_string

   !> Read or write a fixed-length string, rejecting oversized input
   !>
   !> @param[inout] self Active document
   !> @param[in] key Dotted field name
   !> @param[inout] value Parameter value
   subroutine register_string(self, key, value)
      class(moist_model_parameters_type), intent(inout), target :: self
      character(len=*), intent(in) :: key
      character(len=*), intent(inout) :: value
      character(len=:), allocatable :: text

      text = trim(value)
      call self%register_alloc_string(key, text)
      if (allocated(self%failure)) return
      if (len(text) > len(value)) then
         call fatal_error(self%failure, "Parameter string is too long: "//key)
         return
      end if
      value = text
   end subroutine register_string

end module moist_model_parameters
