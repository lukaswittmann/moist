!> Owned thermodynamic items with cached forward values and reverse sensitivities
module moist_model_moz_thermodynamics_list
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_moz_thermodynamics_type, only: thermodynamic_quantity_type, thermo_input_type
   implicit none(type, external)
   private
   public :: thermo_list_type

   !> List entry: owned item with cached value and sensitivity
   type :: thermo_slot_type
      !> Owned copy of the item
      class(thermodynamic_quantity_type), allocatable :: item
      !> Slot indices of the dependencies (n_dep), set by compile
      integer, allocatable :: dependencies(:)
      !> Cached forward value, in the item unit, and sensitivity dE/dvalue, Hartree per item unit
      real(wp) :: value = 0.0_wp, adjoint = 0.0_wp
      !> Selected as additive energy term; default false
      logical :: selected = .false.
   end type thermo_slot_type

   !> Owned thermodynamic items with cached forward and reverse results
   !>
   !> - Each input state is evaluated once in dependency order, diagnostics included
   !> - Reverse mode reuses that cache and seeds only the base energy and selected corrections
   type :: thermo_list_type
      private
      !> Item slots (n_items), in registration order
      type(thermo_slot_type), allocatable :: slots(:)
      !> Slot indices (n_items) in dependency order, set by compile
      integer, allocatable :: order(:)
      !> Normalized name of the uncorrected energy item
      character(len=:), allocatable :: base_energy
      !> Forward cache valid, reverse sensitivities valid; both default false
      logical :: evaluated = .false., differentiated = .false.
   contains
      procedure :: add => list_add
      procedure :: add_input => list_add_input
      procedure :: select => list_select
      procedure :: set_base_energy => list_set_base_energy
      procedure :: set_input => list_set_input
      procedure :: evaluate => list_evaluate
      procedure :: energy => list_energy
      procedure :: accumulate_adjoint => list_accumulate_adjoint
      procedure :: get_value => list_get_value
      procedure :: get_adjoint => list_get_adjoint
      procedure :: n_items => list_n_items
      procedure :: get_name => list_get_name
      procedure, private :: find => list_find
      procedure, private :: compile => list_compile
   end type thermo_list_type

contains

   !> Append an owned copy of an item
   !>
   !> - Registration order need not be dependency order
   !> - Name, group and dependency names are stored in lower case
   !> - Error for a missing, empty or duplicate name, a non-Hartree correction or a conflicting selection
   !> - Invalidates cached values and sensitivities
   !>
   !> @param[inout] self      item list
   !> @param[in]    item      item to copy; name set, unit Hartree if a correction
   !> @param[out]   error     error handler, allocated on failure
   !> @param[in]    selected  select as energy term, corrections only; default false
   subroutine list_add(self, item, error, selected)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Item to copy
      class(thermodynamic_quantity_type), intent(in) :: item
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      !> Select as energy term, corrections only; default false
      logical, optional, intent(in) :: selected
      type(thermo_slot_type), allocatable :: slots(:)
      integer :: n, i

      if (.not. allocated(item%name)) then
         call fatal_error(error, "A thermodynamic item needs a name")
         return
      end if
      if (len_trim(item%name) == 0 .or. self%find(item%name) /= 0) then
         call fatal_error(error, "Empty or duplicate thermodynamic item: "//item%name)
         return
      end if
      if (item%correction) then
         if (.not. allocated(item%unit)) then
            call fatal_error(error, "An energy correction must declare its Hartree unit")
            return
         end if
         if (canonical(item%unit) /= "hartree") then
            call fatal_error(error, "Energy corrections must return Hartree: "//item%name)
            return
         end if
      end if
      n = self%n_items()
      allocate (slots(n + 1))
      if (n > 0) slots(:n) = self%slots
      allocate (slots(n + 1)%item, source=item)
      associate (owned => slots(n + 1)%item)
         owned%name = canonical(item%name)
         if (allocated(owned%exclusive_group)) owned%exclusive_group = canonical(owned%exclusive_group)
         if (allocated(owned%dependencies)) then
            do i = 1, size(owned%dependencies)
               owned%dependencies(i) = canonical(owned%dependencies(i))
            end do
         else
            allocate (character(len=1) :: owned%dependencies(0))
         end if
      end associate
      ! Selection is checked on a temporary list so a rejection does not mutate self
      if (present(selected)) then
         if (selected) then
            if (.not. item%correction) then
               call fatal_error(error, "Only corrections may be selected: "//item%name)
               return
            end if
            do i = 1, n
               if (.not. slots(i)%selected) cycle
               if (same_group(slots(i)%item, slots(n + 1)%item)) then
                  call fatal_error(error, "Mutually exclusive corrections: "//slots(i)%item%name//" and "//item%name)
                  return
               end if
            end do
            slots(n + 1)%selected = .true.
         end if
      end if
      call move_alloc(slots, self%slots)
      if (allocated(self%order)) deallocate (self%order)
      self%evaluated = .false.
      self%differentiated = .false.
   end subroutine list_add

   !> Append a scalar input leaf
   !>
   !> - Invalidates cached values and sensitivities
   !>
   !> @param[inout] self   item list
   !> @param[in]    name   input name, case-insensitive
   !> @param[in]    value  scalar value, in unit
   !> @param[out]   error  error handler, allocated on failure
   !> @param[in]    unit   unit label; default unset
   subroutine list_add_input(self, name, value, error, unit)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Input name, case-insensitive
      character(len=*), intent(in) :: name
      !> Scalar value, in unit
      real(wp), intent(in) :: value
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      !> Unit label; default unset
      character(len=*), optional, intent(in) :: unit
      type(thermo_input_type) :: item
      item%name = name
      if (present(unit)) item%unit = unit
      item%scalar = value
      call self%add(item, error)
   end subroutine list_add_input

   !> Select or deselect a correction as an energy term
   !>
   !> - Reuses forward values, invalidates reverse sensitivities
   !> - Error for an unknown or non-correction item, or a second member of an exclusive group
   !>
   !> @param[inout] self      item list
   !> @param[in]    name      item name, case-insensitive
   !> @param[in]    selected  new selection state
   !> @param[out]   error     error handler, allocated on failure
   subroutine list_select(self, name, selected, error)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Item name, case-insensitive
      character(len=*), intent(in) :: name
      !> New selection state
      logical, intent(in) :: selected
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i, index
      index = self%find(name)
      if (index == 0) then
         call fatal_error(error, "Unknown thermodynamic item: "//name)
         return
      end if
      if (selected) then
         if (.not. self%slots(index)%item%correction) then
            call fatal_error(error, "Only corrections may be selected: "//name)
            return
         end if
         do i = 1, self%n_items()
            if (i == index .or. .not. self%slots(i)%selected) cycle
            if (same_group(self%slots(i)%item, self%slots(index)%item)) then
               call fatal_error(error, "Mutually exclusive corrections: "//self%slots(i)%item%name//" and "//name)
               return
            end if
         end do
      end if
      self%slots(index)%selected = selected
      self%differentiated = .false.
   end subroutine list_select

   !> Set the uncorrected energy root
   !>
   !> - Separate root, never a correction; must declare the Hartree unit
   !> - Invalidates reverse sensitivities
   !>
   !> @param[inout] self   item list
   !> @param[in]    name   base energy item name, case-insensitive
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_set_base_energy(self, name, error)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Base energy item name, case-insensitive
      character(len=*), intent(in) :: name
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      integer :: index
      index = self%find(name)
      if (index == 0) then
         call fatal_error(error, "Unknown base energy: "//name)
         return
      end if
      if (self%slots(index)%item%correction) then
         call fatal_error(error, "A correction cannot also be the base energy")
         return
      end if
      if (.not. allocated(self%slots(index)%item%unit)) then
         call fatal_error(error, "The base energy must declare its Hartree unit")
         return
      end if
      if (canonical(self%slots(index)%item%unit) /= "hartree") then
         call fatal_error(error, "The base energy must return Hartree")
         return
      end if
      self%base_energy = canonical(name)
      self%differentiated = .false.
   end subroutine list_set_base_energy

   !> Replace a scalar input value
   !>
   !> - Invalidates cached values and sensitivities
   !> - Error for an unknown name or a non-input item
   !>
   !> @param[inout] self   item list
   !> @param[in]    name   input name, case-insensitive
   !> @param[in]    value  new scalar value, in the input unit
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_set_input(self, name, value, error)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Input name, case-insensitive
      character(len=*), intent(in) :: name
      !> New scalar value, in the input unit
      real(wp), intent(in) :: value
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      integer :: index
      index = self%find(name)
      if (index == 0) then
         call fatal_error(error, "Unknown thermodynamic input: "//name)
         return
      end if
      select type (item => self%slots(index)%item)
      type is (thermo_input_type)
         item%scalar = value
      class default
         call fatal_error(error, "Thermodynamic item is not an input: "//name)
         return
      end select
      self%evaluated = .false.
      self%differentiated = .false.
   end subroutine list_set_input

   !> Resolve dependencies into slot indices and a topological order
   !>
   !> - No-op if the order exists
   !> - Error for a missing dependency or a cycle
   !>
   !> @param[inout] self   item list
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_compile(self, error)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      logical, allocatable :: done(:)
      integer, allocatable :: order(:)
      integer :: i, j, index, n, count, previous
      if (allocated(self%order)) return
      n = self%n_items()
      allocate (done(n), order(n))
      done = .false.
      do i = 1, n
         if (allocated(self%slots(i)%dependencies)) deallocate (self%slots(i)%dependencies)
         allocate (self%slots(i)%dependencies(size(self%slots(i)%item%dependencies)))
         do j = 1, size(self%slots(i)%dependencies)
            index = self%find(self%slots(i)%item%dependencies(j))
            if (index == 0) then
               call fatal_error(error, "Missing dependency of "//self%slots(i)%item%name//": "// &
                  & trim(self%slots(i)%item%dependencies(j)))
               return
            end if
            self%slots(i)%dependencies(j) = index
         end do
      end do
      count = 0
      do while (count < n)
         previous = count
         do i = 1, n
            if (done(i)) cycle
            if (.not. all(done(self%slots(i)%dependencies))) cycle
            count = count + 1
            order(count) = i
            done(i) = .true.
         end do
         if (count == previous) then
            call fatal_error(error, "Cycle in thermodynamic item dependencies")
            return
         end if
      end do
      call move_alloc(order, self%order)
   end subroutine list_compile

   !> Evaluate every item in dependency order
   !>
   !> - Reuses values until inputs change
   !> - Error for a failing item or a non-finite value
   !>
   !> @param[inout] self   item list
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_evaluate(self, error)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: values(:)
      integer :: i, index
      if (self%evaluated) return
      self%differentiated = .false.
      call self%compile(error)
      if (allocated(error)) return
      do i = 1, size(self%order)
         index = self%order(i)
         values = self%slots(self%slots(index)%dependencies)%value
         call self%slots(index)%item%evaluate(values, self%slots(index)%value, error)
         if (allocated(error)) return
         if (.not. ieee_is_finite(self%slots(index)%value)) then
            call fatal_error(error, "Non-finite thermodynamic item: "//self%slots(index)%item%name)
            return
         end if
      end do
      self%evaluated = .true.
   end subroutine list_evaluate

   !> Sum base energy and selected corrections only, in Hartree
   !>
   !> - Requires evaluate; zero base energy if none is set
   !>
   !> @param[in]    self    item list
   !> @param[out]   energy  total energy, Hartree
   !> @param[out]   error   error handler, allocated on failure
   subroutine list_energy(self, energy, error)
      !> Item list
      class(thermo_list_type), intent(in) :: self
      !> Total energy, Hartree
      real(wp), intent(out) :: energy
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      energy = 0.0_wp
      if (.not. self%evaluated) then
         call fatal_error(error, "Evaluate thermodynamic items before reading their energy")
         return
      end if
      if (allocated(self%base_energy)) energy = self%slots(self%find(self%base_energy))%value
      do i = 1, self%n_items()
         if (self%slots(i)%selected) energy = energy + self%slots(i)%value
      end do
   end subroutine list_energy

   !> Seed base energy and selected corrections, then accumulate in reverse order
   !>
   !> - Requires evaluate
   !> - Error for non-finite sensitivities
   !>
   !> @param[inout] self   item list
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_accumulate_adjoint(self, error)
      !> Item list
      class(thermo_list_type), intent(inout) :: self
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: values(:), adjoints(:)
      integer :: i, j, index, dependency
      self%differentiated = .false.
      if (.not. self%evaluated) then
         call fatal_error(error, "Evaluate thermodynamic items before differentiating them")
         return
      end if
      if (self%n_items() == 0) then
         self%differentiated = .true.
         return
      end if
      self%slots%adjoint = 0.0_wp
      if (allocated(self%base_energy)) self%slots(self%find(self%base_energy))%adjoint = 1.0_wp
      do i = 1, self%n_items()
         if (self%slots(i)%selected) self%slots(i)%adjoint = self%slots(i)%adjoint + 1.0_wp
      end do
      do i = size(self%order), 1, -1
         index = self%order(i)
         if (self%slots(index)%adjoint == 0.0_wp) cycle
         values = self%slots(self%slots(index)%dependencies)%value
         if (allocated(adjoints)) deallocate (adjoints)
         allocate (adjoints(size(values)), source=0.0_wp)
         call self%slots(index)%item%accumulate_adjoint(values, self%slots(index)%value, &
            & self%slots(index)%adjoint, adjoints)
         do j = 1, size(adjoints)
            dependency = self%slots(index)%dependencies(j)
            self%slots(dependency)%adjoint = self%slots(dependency)%adjoint + adjoints(j)
         end do
      end do
      if (.not. all(ieee_is_finite(self%slots%adjoint))) then
         call fatal_error(error, "Non-finite thermodynamic sensitivities")
         return
      end if
      self%differentiated = .true.
   end subroutine list_accumulate_adjoint

   !> Read a cached scalar in its declared unit
   !>
   !> - Requires a successful evaluate
   !>
   !> @param[in]    self   item list
   !> @param[in]    name   item name, case-insensitive
   !> @param[out]   value  cached value, in the item unit; zero on error
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_get_value(self, name, value, error)
      !> Item list
      class(thermo_list_type), intent(in) :: self
      !> Item name, case-insensitive
      character(len=*), intent(in) :: name
      !> Cached value, in the item unit; zero on error
      real(wp), intent(out) :: value
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      integer :: index
      value = 0.0_wp
      index = self%find(name)
      if (index == 0 .or. .not. self%evaluated) then
         call fatal_error(error, "No evaluated thermodynamic item: "//name)
         return
      end if
      value = self%slots(index)%value
   end subroutine list_get_value

   !> Read dE/dvalue in Hartree per item unit
   !>
   !> - Requires accumulate_adjoint
   !>
   !> @param[in]    self   item list
   !> @param[in]    name   item name, case-insensitive
   !> @param[out]   value  sensitivity, Hartree per item unit; zero on error
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_get_adjoint(self, name, value, error)
      !> Item list
      class(thermo_list_type), intent(in) :: self
      !> Item name, case-insensitive
      character(len=*), intent(in) :: name
      !> Sensitivity, Hartree per item unit; zero on error
      real(wp), intent(out) :: value
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      integer :: index
      value = 0.0_wp
      index = self%find(name)
      if (index == 0 .or. .not. self%differentiated) then
         call fatal_error(error, "No accumulated thermodynamic sensitivity: "//name)
         return
      end if
      value = self%slots(index)%adjoint
   end subroutine list_get_adjoint

   !> Number of items in the list
   !>
   !> @param[in]    self  item list
   pure function list_n_items(self) result(n)
      !> Item list
      class(thermo_list_type), intent(in) :: self
      !> Item count
      integer :: n
      n = 0
      if (allocated(self%slots)) n = size(self%slots)
   end function list_n_items

   !> Name of an item by registration index
   !>
   !> @param[in]    self   item list
   !> @param[in]    index  item index, 1 to n_items
   !> @param[out]   name   normalized item name
   !> @param[out]   error  error handler, allocated on failure
   subroutine list_get_name(self, index, name, error)
      !> Item list
      class(thermo_list_type), intent(in) :: self
      !> Item index, 1 to n_items
      integer, intent(in) :: index
      !> Normalized item name
      character(len=:), allocatable, intent(out) :: name
      !> Error handler, allocated on failure
      type(error_type), allocatable, intent(out) :: error
      if (index < 1 .or. index > self%n_items()) then
         call fatal_error(error, "Thermodynamic item index is out of range")
         return
      end if
      name = self%slots(index)%item%name
   end subroutine list_get_name

   !> Slot index of an item by name, zero if absent
   !>
   !> @param[in]    self  item list
   !> @param[in]    name  item name, case-insensitive
   pure function list_find(self, name) result(index)
      !> Item list
      class(thermo_list_type), intent(in) :: self
      !> Item name, case-insensitive
      character(len=*), intent(in) :: name
      !> Slot index, zero if not found
      integer :: index, i
      character(len=:), allocatable :: normalized
      index = 0
      normalized = canonical(name)
      do i = 1, self%n_items()
         if (self%slots(i)%item%name == normalized) then
            index = i
            return
         end if
      end do
   end function list_find

   !> Normalize a name: trim, adjust left, lower case
   !>
   !> @param[in]    name  name to normalize
   pure function canonical(name) result(normalized)
      !> Name to normalize
      character(len=*), intent(in) :: name
      ! Automatic character results keep name normalization reentrant on GCC 13
      !> Normalized name
      character(len=len_trim(adjustl(name))) :: normalized
      normalized = to_lower(trim(adjustl(name)))
   end function canonical

   !> Test whether two items share a nonempty exclusive group
   !>
   !> - Group labels compared case-insensitive
   !>
   !> @param[in]    a  first item
   !> @param[in]    b  second item
   pure function same_group(a, b) result(same)
      !> Items to compare
      class(thermodynamic_quantity_type), intent(in) :: a, b
      !> True if both items belong to the same nonempty group
      logical :: same
      same = .false.
      if (.not. allocated(a%exclusive_group) .or. .not. allocated(b%exclusive_group)) return
      if (len_trim(a%exclusive_group) == 0) return
      same = canonical(a%exclusive_group) == canonical(b%exclusive_group)
   end function same_group

end module moist_model_moz_thermodynamics_list
