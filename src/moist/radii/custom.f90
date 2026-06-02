module moist_radii_custom
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   implicit none
   private

   public :: custom_radius_type
   public :: new_custom_radii_atoms
   public :: new_custom_radii_elements

   !> Custom radii model with user-supplied, geometry-invariant radii.
   type, extends(radius_type) :: custom_radius_type
      !> If true, radii are supplied per atom in molecular order.
      logical :: has_atom_radii = .false.
      !> If true, radii are supplied per element by atomic number.
      logical :: has_element_radii = .false.
      !> Stored per-atom radii.
      real(wp), allocatable :: atom_radii(:)
      !> Lookup table for per-element radii indexed by atomic number.
      real(wp), allocatable :: element_radii(:)
   contains
      !> Update cached custom radii and zero derivatives.
      procedure :: update => update_custom_radii
      !> Print custom radii model status.
      procedure :: print => print_custom_radii
   end type custom_radius_type

contains

   !> Build a custom model from per-atom radii.
   !> @param[in]  radii  per-atom radii (bohr)
   !> @param[out] self       custom radii model
   !> @param[out] error      error on invalid radii input
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_custom_radii_atoms(radii, self, error, verbosity)
      !> Per-atom radii (bohr)
      real(wp), intent(in) :: radii(:)
      !> Custom radii model
      type(custom_radius_type), intent(out) :: self
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Optional print level
      integer, intent(in), optional :: verbosity

      integer :: iat

      self%nat = 0
      self%has_atom_radii = .false.
      self%has_element_radii = .false.
      if (present(verbosity)) then
         self%verbosity = verbosity
      else
         self%verbosity = 0
      end if

      if (size(radii) < 1) then
         call fatal_error(error, "new_custom_radii_atoms: radii list must not be empty")
         return
      end if

      do iat = 1, size(radii)
         if (radii(iat) <= 0.0_wp) then
            call fatal_error(error, "new_custom_radii_atoms: all radii must be positive")
            return
         end if
      end do

      allocate (self%atom_radii(size(radii)))
      self%atom_radii = radii
      self%has_atom_radii = .true.
   end subroutine new_custom_radii_atoms

   !> Build a custom model from per-element radii.
   !> @param[in]  atomic_numbers  atomic numbers for supplied radii
   !> @param[in]  radii           radii matching atomic_numbers (bohr)
   !> @param[out] self            custom radii model
   !> @param[out] error           error on invalid input
   !> @param[in]  verbosity       optional print level for diagnostics
   subroutine new_custom_radii_elements(atomic_numbers, radii, self, error, verbosity)
      !> Atomic numbers for supplied radii
      integer, intent(in) :: atomic_numbers(:)
      !> Radii values matching atomic_numbers (bohr)
      real(wp), intent(in) :: radii(:)
      !> Custom radii model
      type(custom_radius_type), intent(out) :: self
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Optional print level
      integer, intent(in), optional :: verbosity

      integer :: iat
      integer :: z
      integer :: max_z
      character(len=128) :: msg

      self%nat = 0
      self%has_atom_radii = .false.
      self%has_element_radii = .false.
      if (present(verbosity)) then
         self%verbosity = verbosity
      else
         self%verbosity = 0
      end if

      if (size(atomic_numbers) < 1) then
         call fatal_error(error, "new_custom_radii_elements: atomic_numbers must not be empty")
         return
      end if
      if (size(atomic_numbers) /= size(radii)) then
         write (msg, '(a,i0,a,i0,a)') "new_custom_radii_elements: atomic_numbers size (", &
            size(atomic_numbers), ") does not match radii size (", size(radii), ")."
         call fatal_error(error, trim(msg))
         return
      end if

      max_z = 0
      do iat = 1, size(atomic_numbers)
         z = atomic_numbers(iat)
         if (z < 1) then
            call fatal_error(error, "new_custom_radii_elements: atomic numbers must be >= 1")
            return
         end if
         if (radii(iat) <= 0.0_wp) then
            call fatal_error(error, "new_custom_radii_elements: all radii must be positive")
            return
         end if
         max_z = max(max_z, z)
      end do

      allocate (self%element_radii(max_z), source=-1.0_wp)
      do iat = 1, size(atomic_numbers)
         z = atomic_numbers(iat)
         if (self%element_radii(z) > 0.0_wp) then
            call fatal_error(error, "new_custom_radii_elements: duplicate atomic number in input")
            return
         end if
         self%element_radii(z) = radii(iat)
      end do

      self%has_element_radii = .true.
   end subroutine new_custom_radii_elements

   !> Print custom radii model status.
   !> @param[in] self  custom radii model
   !> @param[in] unit  optional output unit
   subroutine print_custom_radii(self, unit)
      !> Custom radii model
      class(custom_radius_type), intent(in) :: self
      !> Optional output unit
      integer, intent(in), optional :: unit

      integer :: iu

      iu = output_unit
      if (present(unit)) iu = unit

      write (iu, '(a)') "Custom radii model:"
      write (iu, '(a,l1)') "  has_atom_radii: ", self%has_atom_radii
      write (iu, '(a,l1)') "  has_element_radii: ", self%has_element_radii
      write (iu, '(a,i0)') "  cached nat: ", self%nat

      if (self%has_atom_radii .and. allocated(self%atom_radii)) then
         write (iu, '(a,i0)') "  number of atom radii: ", size(self%atom_radii)
      end if
      if (self%has_element_radii .and. allocated(self%element_radii)) then
         write (iu, '(a,i0)') "  max atomic number in element radii table: ", size(self%element_radii)
      end if
   end subroutine print_custom_radii

   !> Update cached radii for a molecular structure.
   !> @param[inout] self   custom radii model
   !> @param[in]    mol    molecular structure
   !> @param[out]   error  error handle on invalid setup/input
   subroutine update_custom_radii(self, mol, error)
      !> Custom radii model
      class(custom_radius_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer :: iat
      integer :: z
      character(len=128) :: msg

      self%nat = mol%nat
      if (self%nat < 1) then
         write (msg, '(a,i0)') "Invalid atom count in structure: ", self%nat
         call fatal_error(error, trim(msg))
         return
      end if

      if (allocated(self%f0)) deallocate (self%f0)
      if (allocated(self%f1_rA)) deallocate (self%f1_rA)
      allocate (self%f0(self%nat), self%f1_rA(3, self%nat, self%nat))
      self%f1_rA = 0.0_wp

      if (self%has_atom_radii .eqv. self%has_element_radii) then
         call fatal_error(error, "Custom radii model must define exactly one of atom or element radii")
         return
      end if

      if (self%has_atom_radii) then
         if (.not. allocated(self%atom_radii)) then
            call fatal_error(error, "Custom atom radii are not allocated")
            return
         end if
         if (size(self%atom_radii) /= self%nat) then
            write (msg, '(a,i0,a,i0,a)') "Custom atom radii size (", &
               size(self%atom_radii), ") does not match mol%nat (", self%nat, ")."
            call fatal_error(error, trim(msg))
            return
         end if
         self%f0 = self%atom_radii
         return
      end if

      if (.not. allocated(self%element_radii)) then
         call fatal_error(error, "Custom element radii are not allocated")
         return
      end if

      do iat = 1, self%nat
         z = mol%num(mol%id(iat))
         if (z < 1 .or. z > size(self%element_radii)) then
            write (msg, '(a,i0,a)') "No custom radius provided for atomic number ", z, "."
            call fatal_error(error, trim(msg))
            return
         end if
         if (self%element_radii(z) <= 0.0_wp) then
            write (msg, '(a,i0,a)') "No custom radius provided for atomic number ", z, "."
            call fatal_error(error, trim(msg))
            return
         end if
         self%f0(iat) = self%element_radii(z)
      end do
   end subroutine update_custom_radii

end module moist_radii_custom
