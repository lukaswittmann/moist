!> Infrastructure for declarative parameter tables with JSON IO
!>
!> ## Idea
!> Define each parameter once (key, target variable, default) and let the base class
!> iterate a polymorphic list of parameter bindings to handle read/write/print
!>
!> ## Basic Usage
!> 1. Extend `moist_model_parameters_type`
!> 2. In `init_defaults`, set all default values
!> 3. In `register_entries`, call `register_*` helpers for every parameter
!> 4. Call `read_file`, `write_file`, or `print_parameters` from child code
!>
!> ## Supported Parameter Types
!> - Real scalars (`real(wp)`)
!> - Integer scalars (`integer`)
!> - Logical scalars (`logical`)
!> - Character strings (`character(len=*)`)
!> - Real vectors with fixed size
module moist_model_parameters
   use jonquil, only: json_object, json_value, json_array, json_error, json_load, &
      & get_value, cast_to_object, cast_to_array, len
   use mctc_env_accuracy, only: wp
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   use, intrinsic :: iso_fortran_env, only: output_unit

   implicit none(type, external)
   private

   !> Base interface for a single parameter entry
   type, abstract :: parameter_binding
      !> JSON key identifying the parameter
      character(len=:), allocatable :: key
   contains
      !> Read the parameter from a JSON object
      procedure(read_json_ifc), deferred :: read_json
      !> Write the parameter to a JSON object
      procedure(write_json_ifc), deferred :: write_json
      procedure(print_value_ifc), deferred :: print_value
   end type parameter_binding

   !> Wrapper node stored inside the list
   type :: parameter_node
      !> Actual polymorphic binding
      class(parameter_binding), allocatable :: item
   end type parameter_node

   !> Binding for a real scalar parameter
   type, extends(parameter_binding) :: real_scalar_parameter
      !> Pointer to the target variable
      real(wp), pointer :: value => null()
   contains
      procedure :: init => init_real_scalar
      procedure :: read_json => read_real_scalar
      procedure :: write_json => write_real_scalar
      procedure :: print_value => print_real_scalar
   end type real_scalar_parameter

   !> Binding for an integer scalar parameter
   type, extends(parameter_binding) :: int_scalar_parameter
      !> Pointer to the target variable
      integer, pointer :: value => null()
   contains
      procedure :: init => init_int_scalar
      procedure :: read_json => read_int_scalar
      procedure :: write_json => write_int_scalar
      procedure :: print_value => print_int_scalar
   end type int_scalar_parameter

   !> Binding for a real vector parameter
   type, extends(parameter_binding) :: real_vector_parameter
      !> Pointer to the target array
      real(wp), pointer :: value(:) => null()
   contains
      procedure :: init => init_real_vector
      procedure :: read_json => read_real_vector
      procedure :: write_json => write_real_vector
      procedure :: print_value => print_real_vector
   end type real_vector_parameter

   !> Binding for a character string parameter
   type, extends(parameter_binding) :: string_parameter
      !> Pointer to the target character variable
      character(len=:), pointer :: value => null()
   contains
      procedure :: init => init_string
      procedure :: read_json => read_string
      procedure :: write_json => write_string
      procedure :: print_value => print_string
   end type string_parameter

   !> Binding for a logical scalar parameter
   type, extends(parameter_binding) :: logical_parameter
      !> Pointer to the target logical variable
      logical, pointer :: value => null()
   contains
      procedure :: init => init_logical
      procedure :: read_json => read_logical
      procedure :: write_json => write_logical
      procedure :: print_value => print_logical
   end type logical_parameter

   public :: moist_model_parameters_type

   !> Abstract base class for model-specific parameter sets
   type, abstract :: moist_model_parameters_type
      private
      !> Registered parameter list
      type(parameter_node), allocatable :: params(:)
      !> Flag indicating whether entries were registered
      logical :: entries_registered = .false.
   contains
      procedure(register_entries_ifc), deferred :: register_entries
      procedure(init_defaults_ifc), deferred :: init_defaults
      procedure                                :: read_file
      procedure                                :: write_file
      procedure                                :: print_parameters
      procedure, pass :: register_real_scalar
      procedure, pass :: register_int_scalar
      procedure, pass :: register_real_vector
      procedure, pass :: register_string
      procedure, pass :: register_logical
      procedure, private :: ensure_entries
      procedure, private :: push_parameter
      procedure, private :: read_entries
      procedure, private :: write_entries
      procedure, private :: print_entries
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
         class(moist_model_parameters_type), intent(inout) :: self
      end subroutine register_entries_ifc

   end interface

   abstract interface
      subroutine read_json_ifc(self, json)
         import :: parameter_binding, json_object
         implicit none(type, external)
         class(parameter_binding), intent(inout) :: self
         type(json_object), pointer :: json
      end subroutine read_json_ifc

      subroutine write_json_ifc(self, unit, is_last, indent)
         import :: parameter_binding
         implicit none(type, external)
         class(parameter_binding), intent(in) :: self
         integer, intent(in) :: unit
         logical, intent(in) :: is_last
         integer, intent(in) :: indent
      end subroutine write_json_ifc

      subroutine print_value_ifc(self, pp)
         import :: parameter_binding, prettyprinter
         implicit none(type, external)
         class(parameter_binding), intent(in) :: self
         type(prettyprinter), intent(inout) :: pp
      end subroutine print_value_ifc
   end interface

contains

   !> Get group path from a dotted key (everything before the last dot).
   !> e.g. "solver.gmres.restart" -> "solver.gmres", "grid.nr" -> "grid", "theory" -> ""
   pure function key_group(key) result(grp)
      character(*), intent(in) :: key
      character(:), allocatable :: grp
      integer :: pos
      pos = index(key, '.', back=.true.)
      if (pos > 0) then
         grp = key(1:pos - 1)
      else
         grp = ''
      end if
   end function key_group

   !> Get leaf name from a dotted key (everything after the last dot).
   !> e.g. "solver.gmres.restart" -> "restart", "grid.nr" -> "nr", "theory" -> "theory"
   pure function key_leaf(key) result(leaf)
      character(*), intent(in) :: key
      character(:), allocatable :: leaf
      integer :: pos
      pos = index(key, '.', back=.true.)
      if (pos > 0) then
         leaf = key(pos + 1:)
      else
         leaf = trim(key)
      end if
   end function key_leaf

   !> Count dot-separated segments in a group path.
   !> Returns 0 for empty string, 1 for "grid", 2 for "solver.gmres", etc.
   pure function path_depth(path) result(d)
      character(*), intent(in) :: path
      integer :: d, i

      if (len_trim(path) == 0) then
         d = 0
         return
      end if
      d = 1
      do i = 1, len_trim(path)
         if (path(i:i) == '.') d = d + 1
      end do
   end function path_depth

   !> Get the k-th segment (1-based) from a dot-separated path.
   !> e.g. path_segment("solver.gmres", 1) = "solver", path_segment("solver.gmres", 2) = "gmres"
   pure function path_segment(path, k) result(seg)
      character(*), intent(in) :: path
      integer, intent(in) :: k
      character(:), allocatable :: seg
      integer :: start, count, i

      start = 1
      count = 1
      do i = 1, len_trim(path)
         if (path(i:i) == '.') then
            if (count == k) then
               seg = path(start:i - 1)
               return
            end if
            start = i + 1
            count = count + 1
         end if
      end do
      if (count == k) then
         seg = path(start:len_trim(path))
      else
         seg = ''
      end if
   end function path_segment

   !> Check whether two group paths share the same first d segments.
   !> Always returns .true. when d == 0.
   pure function paths_share_prefix(a, b, d) result(shared)
      character(*), intent(in) :: a, b
      integer, intent(in) :: d
      logical :: shared
      integer :: k

      shared = .true.
      do k = 1, d
         if (path_segment(a, k) /= path_segment(b, k)) then
            shared = .false.
            return
         end if
      end do
   end function paths_share_prefix

   !> Navigate a dotted key path to the parent JSON table and leaf key.
   !> Descends through all intermediate dots, e.g. "solver.gmres.restart"
   !> traverses root -> solver -> gmres, returning parent => gmres, leaf = "restart".
   subroutine resolve_json_path(root, full_key, parent, leaf)
      type(json_object), pointer :: root
      character(*), intent(in) :: full_key
      type(json_object), pointer, intent(out) :: parent
      character(:), allocatable, intent(out) :: leaf
      type(json_object), pointer :: current, child
      integer :: pos, start

      current => root
      start = 1

      do
         pos = index(full_key(start:), '.')
         if (pos == 0) then
            parent => current
            leaf = trim(full_key(start:))
            return
         end if
         pos = start + pos - 1

         call get_value(current, full_key(start:pos - 1), child, requested=.false.)
         if (.not. associated(child)) then
            nullify (parent)
            leaf = trim(full_key(start:))
            return
         end if
         current => child
         start = pos + 1
      end do
   end subroutine resolve_json_path

   !> Read all registered parameters from a JSON file
   subroutine read_file(self, filepath)
      use jonquil, only: json_value, cast_to_object
      class(moist_model_parameters_type), intent(inout) :: self
      !> Path to the JSON parameter file
      character(len=*), intent(in) :: filepath
      !> JSON document handle
      class(json_value), allocatable :: json_val
      type(json_object), pointer :: json
      type(json_error), allocatable :: error

      call self%ensure_entries()
      call self%init_defaults()

      call json_load(json_val, trim(filepath), error=error)
      if (allocated(error)) then
         write (output_unit, '(a,a)') '[Error] Failed to load parameter file: ', trim(filepath)
         write (output_unit, '(a)') error%message
         stop
      end if

      json => cast_to_object(json_val)
      if (.not. associated(json)) then
         write (output_unit, '(a)') '[Error] JSON root is not an object'
         stop
      end if

      call self%read_entries(json)
   end subroutine read_file

   !> Write all registered parameters to a JSON file
   subroutine write_file(self, filepath)
      class(moist_model_parameters_type), intent(inout) :: self
      !> Path to the JSON parameter file
      character(len=*), intent(in) :: filepath
      integer :: unit

      call self%ensure_entries()

      open (newunit=unit, file=trim(filepath), status='replace', action='write')
      write (unit, '(a)') '{'
      call self%write_entries(unit)
      write (unit, '(a)') '}'
      close (unit)
   end subroutine write_file

   !> Print registered parameters in a human-readable form
   subroutine print_parameters(self)
      class(moist_model_parameters_type), intent(inout) :: self

      call self%ensure_entries()
      call self%print_entries()
   end subroutine print_parameters

   !> Ensure that `register_entries` has been executed
   subroutine ensure_entries(self)
      class(moist_model_parameters_type), intent(inout) :: self
      if (.not. self%entries_registered) then
         call self%register_entries()
         self%entries_registered = .true.
      end if
   end subroutine ensure_entries

   !> Append a parameter binding to the internal list
   subroutine push_parameter(self, binding)
      class(moist_model_parameters_type), intent(inout) :: self
      !> Parameter binding to append
      class(parameter_binding), intent(in) :: binding
      type(parameter_node), allocatable :: tmp(:)
      integer :: n, new_size

      if (.not. allocated(self%params)) then
         allocate (self%params(1))
         new_size = 1
      else
         n = size(self%params)
         allocate (tmp(n + 1))
         tmp(1:n) = self%params
         call move_alloc(tmp, self%params)
         new_size = n + 1
      end if

      allocate (self%params(new_size)%item, source=binding)
   end subroutine push_parameter

   !> Iterate over the list and invoke each entry's read method
   subroutine read_entries(self, json)
      class(moist_model_parameters_type), intent(inout) :: self
      !> JSON document handle
      type(json_object), pointer        :: json
      integer :: i

      if (.not. allocated(self%params)) return
      do i = 1, size(self%params)
         call self%params(i)%item%read_json(json)
      end do
   end subroutine read_entries

   !> Iterate over the list and invoke each entry's write method.
   !> Supports arbitrarily nested dotted-key groups (e.g. "solver.gmres.restart").
   subroutine write_entries(self, unit)
      class(moist_model_parameters_type), intent(in) :: self
      !> Fortran output unit
      integer, intent(in) :: unit
      integer :: i, n, d, curr_depth, prev_depth, cd
      character(:), allocatable :: grp, nxt_grp, prev_grp
      logical :: has_next, val_is_last

      if (.not. allocated(self%params)) return
      n = size(self%params)
      prev_grp = ''
      prev_depth = 0

      do i = 1, n
         grp = key_group(self%params(i)%item%key)
         curr_depth = path_depth(grp)
         has_next = (i < n)
         if (has_next) then
            nxt_grp = key_group(self%params(i + 1)%item%key)
         else
            nxt_grp = ''
         end if

         !> Find common prefix depth between previous and current group paths
         cd = 0
         do d = 1, min(prev_depth, curr_depth)
            if (path_segment(prev_grp, d) /= path_segment(grp, d)) exit
            cd = d
         end do

         !> Close groups that are no longer in common
         do d = prev_depth, cd + 1, -1
            if (paths_share_prefix(prev_grp, grp, d - 1)) then
               write (unit, '(a,a)') repeat('  ', d), '},'
            else
               write (unit, '(a,a)') repeat('  ', d), '}'
            end if
         end do

         !> Open new groups
         do d = cd + 1, curr_depth
            write (unit, '(a,3a)') repeat('  ', d), '"', path_segment(grp, d), '": {'
         end do

         !> Write the value; is_last means no trailing comma
         val_is_last = .not. (has_next .and. paths_share_prefix(grp, nxt_grp, curr_depth))
         call self%params(i)%item%write_json(unit, val_is_last, 2*(curr_depth + 1))

         prev_grp = grp
         prev_depth = curr_depth
      end do

      !> Close all remaining open groups after the last entry
      do d = prev_depth, 1, -1
         write (unit, '(a,a)') repeat('  ', d), '}'
      end do
   end subroutine write_entries

   !> Iterate over the list and invoke each entry's print method.
   !> Supports arbitrarily nested dotted-key groups with push/pop indentation.
   subroutine print_entries(self)
      class(moist_model_parameters_type), intent(in) :: self
      type(prettyprinter) :: pp
      integer :: i, d, curr_depth, prev_depth, cd
      character(:), allocatable :: grp, prev_grp

      if (.not. allocated(self%params)) return

      pp = new_prettyprinter(unit=output_unit)
      call pp%push('Model Parameters:')

      prev_grp = ''
      prev_depth = 0

      do i = 1, size(self%params)
         grp = key_group(self%params(i)%item%key)
         curr_depth = path_depth(grp)

         !> Find common prefix depth
         cd = 0
         do d = 1, min(prev_depth, curr_depth)
            if (path_segment(prev_grp, d) /= path_segment(grp, d)) exit
            cd = d
         end do

         !> Pop groups that are no longer in common
         do d = prev_depth, cd + 1, -1
            call pp%pop()
         end do

         !> Push new groups
         do d = cd + 1, curr_depth
            call pp%blank()
            call pp%push(path_segment(grp, d)//':')
         end do

         call self%params(i)%item%print_value(pp)
         prev_grp = grp
         prev_depth = curr_depth
      end do

      !> Pop remaining groups
      do d = prev_depth, 1, -1
         call pp%pop()
      end do

      !> End with blank line
      call pp%blank()

   end subroutine print_entries

   !> Register a real scalar parameter
   subroutine register_real_scalar(self, key, value)
      class(moist_model_parameters_type), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target real variable
      real(wp), target, intent(inout) :: value
      type(real_scalar_parameter) :: binding

      call binding%init(key, value)
      call self%push_parameter(binding)
   end subroutine register_real_scalar

   !> Register an integer scalar parameter
   subroutine register_int_scalar(self, key, value)
      class(moist_model_parameters_type), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target integer variable
      integer, target, intent(inout) :: value
      type(int_scalar_parameter) :: binding

      call binding%init(key, value)
      call self%push_parameter(binding)
   end subroutine register_int_scalar

   !> Register a real vector parameter (fixed-size array)
   subroutine register_real_vector(self, key, values)
      class(moist_model_parameters_type), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target real array
      real(wp), target, intent(inout) :: values(:)
      type(real_vector_parameter) :: binding

      call binding%init(key, values)
      call self%push_parameter(binding)
   end subroutine register_real_vector

   !> Register a character string parameter
   subroutine register_string(self, key, value)
      class(moist_model_parameters_type), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target character variable
      character(len=*), target, intent(inout) :: value
      type(string_parameter) :: binding

      call binding%init(key, value)
      call self%push_parameter(binding)
   end subroutine register_string

   !> Register a logical scalar parameter
   subroutine register_logical(self, key, value)
      class(moist_model_parameters_type), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target logical variable
      logical, target, intent(inout) :: value
      type(logical_parameter) :: binding

      call binding%init(key, value)
      call self%push_parameter(binding)
   end subroutine register_logical

   !> Associate a key and target variable for a real scalar binding
   subroutine init_real_scalar(self, key, value)
      class(real_scalar_parameter), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target real variable
      real(wp), target, intent(inout) :: value

      self%key = trim(key)
      self%value => value
   end subroutine init_real_scalar

   !> Read a real scalar from JSON
   subroutine read_real_scalar(self, json)
      class(real_scalar_parameter), intent(inout) :: self
      type(json_object), pointer :: json
      type(json_object), pointer :: parent
      character(:), allocatable :: leaf
      real(wp) :: tmp

      call resolve_json_path(json, self%key, parent, leaf)
      if (.not. associated(parent)) return
      tmp = self%value
      call get_value(parent, leaf, tmp)
      self%value = tmp
   end subroutine read_real_scalar

   !> Write a real scalar to JSON
   subroutine write_real_scalar(self, unit, is_last, indent)
      class(real_scalar_parameter), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in) :: is_last
      integer, intent(in) :: indent

      write (unit, '(a,3a,g0.16,a)') repeat(' ', indent), &
         '"', key_leaf(self%key), '": ', self%value, merge(' ', ',', is_last)
   end subroutine write_real_scalar

   !> Print a real scalar value
   subroutine print_real_scalar(self, pp)
      class(real_scalar_parameter), intent(in) :: self
      type(prettyprinter), intent(inout) :: pp

      call pp%kv(key_leaf(self%key), self%value)
   end subroutine print_real_scalar

   !> Associate a key and target variable for an integer scalar binding
   subroutine init_int_scalar(self, key, value)
      class(int_scalar_parameter), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target integer variable
      integer, target, intent(inout) :: value

      self%key = trim(key)
      self%value => value
   end subroutine init_int_scalar

   !> Read an integer scalar from JSON
   subroutine read_int_scalar(self, json)
      class(int_scalar_parameter), intent(inout) :: self
      type(json_object), pointer :: json
      type(json_object), pointer :: parent
      character(:), allocatable :: leaf
      integer :: tmp

      call resolve_json_path(json, self%key, parent, leaf)
      if (.not. associated(parent)) return
      tmp = self%value
      call get_value(parent, leaf, tmp)
      self%value = tmp
   end subroutine read_int_scalar

   !> Write an integer scalar to JSON
   subroutine write_int_scalar(self, unit, is_last, indent)
      class(int_scalar_parameter), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in) :: is_last
      integer, intent(in) :: indent

      write (unit, '(a,3a,i0,a)') repeat(' ', indent), &
         '"', key_leaf(self%key), '": ', self%value, merge(' ', ',', is_last)
   end subroutine write_int_scalar

   !> Print an integer scalar value
   subroutine print_int_scalar(self, pp)
      class(int_scalar_parameter), intent(in) :: self
      type(prettyprinter), intent(inout) :: pp

      call pp%kv(key_leaf(self%key), self%value)
   end subroutine print_int_scalar

   !> Associate a key and target array for a real vector binding
   subroutine init_real_vector(self, key, value)
      class(real_vector_parameter), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target real array
      real(wp), target, intent(inout) :: value(:)

      self%key = trim(key)
      self%value => value
   end subroutine init_real_vector

   !> Read a real vector from JSON (requires fixed-size target)
   subroutine read_real_vector(self, json)
      class(real_vector_parameter), intent(inout) :: self
      type(json_object), pointer :: json
      type(json_object), pointer :: parent
      character(:), allocatable :: leaf
      class(json_value), pointer :: val
      type(json_array), pointer :: arr
      integer :: n, i, arr_len
      real(wp) :: elem

      n = size(self%value)
      if (n <= 0) return

      call resolve_json_path(json, self%key, parent, leaf)
      if (.not. associated(parent)) return

      call parent%get(leaf, val)
      if (.not. associated(val)) return

      arr => cast_to_array(val)
      if (.not. associated(arr)) return

      arr_len = len(arr)
      if (arr_len /= n) then
         write (output_unit, '(a,i0,a,i0)') '[Warning] Array size mismatch for '// &
            trim(self%key)//': expected ', n, ', got ', arr_len
         return
      end if

      do i = 1, min(n, arr_len)
         call get_value(arr, i, elem)
         self%value(i) = elem
      end do
   end subroutine read_real_vector

   !> Write a real vector to JSON
   subroutine write_real_vector(self, unit, is_last, indent)
      class(real_vector_parameter), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in) :: is_last
      integer, intent(in) :: indent
      integer :: n, i

      n = size(self%value)
      if (n <= 0) then
         write (unit, '(a,3a,a)') repeat(' ', indent), &
            '"', key_leaf(self%key), '": []', merge(' ', ',', is_last)
         return
      end if

      write (unit, '(a,3a)', advance='no') repeat(' ', indent), &
         '"', key_leaf(self%key), '": ['
      do i = 1, n
         if (i > 1) write (unit, '(a)', advance='no') ', '
         write (unit, '(g0.16)', advance='no') self%value(i)
      end do
      write (unit, '(a,a)') ']', merge(' ', ',', is_last)
   end subroutine write_real_vector

   !> Print a real vector as a compact 3-wide table
   subroutine print_real_vector(self, pp)
      class(real_vector_parameter), intent(in) :: self
      type(prettyprinter), intent(inout) :: pp
      type(prettylistprinter) :: plp
      integer :: n, nrows, row, col, idx

      n = size(self%value)
      if (n <= 0) then
         call pp%kv(key_leaf(self%key), '[]')
         return
      end if

      call pp%blank()
      plp = new_prettylistprinter( &
            widths=[6, 16, 6, 16, 6, 16], &
            headers=[character(len=16) :: '#', key_leaf(self%key), &
                     '#', key_leaf(self%key), '#', key_leaf(self%key)], &
            unit=pp%iu, &
            offset=pp%indent + 1)
      call plp%print_header()
      call plp%separator()
      nrows = (n + 2)/3
      do row = 1, nrows
         call plp%begin_row()
         do col = 0, 2
            idx = row + col*nrows
            if (idx <= n) then
               call plp%add(idx)
               call plp%add(self%value(idx))
            else
               call plp%skip()
               call plp%skip()
            end if
         end do
         call plp%end_row()
      end do
      call plp%blank()
   end subroutine print_real_vector

   !> Associate a key and target variable for a string binding
   subroutine init_string(self, key, value)
      class(string_parameter), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target character variable
      character(len=*), target, intent(inout) :: value

      self%key = trim(key)
      self%value => value
   end subroutine init_string

   !> Read a string from JSON
   subroutine read_string(self, json)
      class(string_parameter), intent(inout) :: self
      type(json_object), pointer :: json
      type(json_object), pointer :: parent
      character(:), allocatable :: leaf
      character(:), allocatable :: tmp

      call resolve_json_path(json, self%key, parent, leaf)
      if (.not. associated(parent)) return
      call get_value(parent, leaf, tmp)
      if (allocated(tmp)) self%value = tmp
   end subroutine read_string

   !> Write a string to JSON
   subroutine write_string(self, unit, is_last, indent)
      class(string_parameter), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in) :: is_last
      integer, intent(in) :: indent

      write (unit, '(a,3a,3a,a)') repeat(' ', indent), &
         '"', key_leaf(self%key), '": "', trim(self%value), '"', &
         merge(' ', ',', is_last)
   end subroutine write_string

   !> Print a string value
   subroutine print_string(self, pp)
      class(string_parameter), intent(in) :: self
      type(prettyprinter), intent(inout) :: pp

      call pp%kv(key_leaf(self%key), trim(self%value))
   end subroutine print_string

   !> Associate a key and target variable for a logical binding
   subroutine init_logical(self, key, value)
      class(logical_parameter), intent(inout) :: self
      !> JSON key name
      character(len=*), intent(in) :: key
      !> Target logical variable
      logical, target, intent(inout) :: value

      self%key = trim(key)
      self%value => value
   end subroutine init_logical

   !> Read a logical from JSON
   subroutine read_logical(self, json)
      class(logical_parameter), intent(inout) :: self
      type(json_object), pointer :: json
      type(json_object), pointer :: parent
      character(:), allocatable :: leaf
      logical :: tmp

      call resolve_json_path(json, self%key, parent, leaf)
      if (.not. associated(parent)) return
      tmp = self%value
      call get_value(parent, leaf, tmp)
      self%value = tmp
   end subroutine read_logical

   !> Write a logical to JSON
   subroutine write_logical(self, unit, is_last, indent)
      class(logical_parameter), intent(in) :: self
      integer, intent(in) :: unit
      logical, intent(in) :: is_last
      integer, intent(in) :: indent

      write (unit, '(a,3a,a,a)') repeat(' ', indent), &
         '"', key_leaf(self%key), '": ', &
         merge('true ', 'false', self%value), merge(' ', ',', is_last)
   end subroutine write_logical

   !> Print a logical value
   subroutine print_logical(self, pp)
      class(logical_parameter), intent(in) :: self
      type(prettyprinter), intent(inout) :: pp

      call pp%kv(key_leaf(self%key), merge('true ', 'false', self%value))
   end subroutine print_logical

end module moist_model_parameters
