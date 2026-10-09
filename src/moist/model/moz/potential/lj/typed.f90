!> Lennard-Jones potential term with parameters by atom type
!>
!> - Structure typing and parameter lookup at build; force field chosen at construction
!> - `lj_set_gaff`: approximate types from element and geometry, GAFF2 Version 2.2.30 parameters
!> - `lj_set_oplsaa`: environment rules from explicit bond graph with all hydrogens
!> - OPLS-AA rule labels for compatible imported LJ keys
!> - OPLS-AA scope: LJ parameters; separate charges and bonded terms
!> - OPLS-AA geometric mixing on interaction
!> - Untyped atoms or missing parameters: atoms left for later terms
!> - Fallback examples: UFF element table, solvent LJ for ow/hw water placeholders
!> - Atoms uncovered by every term: error at set build
!> - Unprocessable structure, such as OPLS-AA without bond graph: error
module moist_model_moz_potential_lj_typed
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_data_atomicrad, only: max_elem
   use moist_model_moz_potential_base, only: potential_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   use moist_model_moz_potential_lj_gaff_typing, only: generate_gaff_atomtypes
   use moist_model_moz_potential_lj_gaff_parameters, only: lookup_gaff_lj
   use moist_model_moz_potential_lj_oplsaa_parameters, only: lookup_oplsaa_lj, oplsaa_key_len
   use moist_model_moz_potential_lj_oplsaa_rules, only: oplsaa_type_len
   use moist_model_moz_potential_lj_oplsaa_typing, only: generate_oplsaa_atomtypes
   use moist_model_moz_potential_utils, only: atom_label, atom_label_len
   implicit none(type, external)
   private

   public :: lj_typed_type, lj_set_gaff, lj_set_oplsaa, lj_atomtype_len
   public :: lj_gaff, lj_oplsaa

   !> GAFF typing with GAFF2 parameters
   integer, parameter :: lj_set_gaff = 1
   !> Native explicit-graph OPLS-AA LJ typing
   integer, parameter :: lj_set_oplsaa = 2

   !> Length of an atom-type label
   integer, parameter :: lj_atomtype_len = max(8, oplsaa_type_len)

   !> Lennard-Jones term by atom type
   type, extends(potential_type) :: lj_typed_type
      !> Force field, `lj_set_gaff` or `lj_set_oplsaa`; 0 until constructed
      integer :: set = 0
      !> Atom-type label per atom after build
      !>
      !> - GAFF: "type/GAFF" for covered atoms (GAFF2 table), bare type when uncovered
      !> - OPLS-AA: rule label, blank when uncovered
      character(len=lj_atomtype_len), allocatable :: atomtype(:)
      !> Source-prefixed LJ key per atom after OPLS-AA typing
      character(len=oplsaa_key_len), allocatable :: parameter_key(:)
   contains
      procedure :: name => lj_typed_name
      procedure :: build => lj_typed_build
   end type lj_typed_type

   !> GAFF2 by approximate GAFF typing
   type(lj_typed_type), parameter :: lj_gaff = lj_typed_type(set=lj_set_gaff)
   !> OPLS-AA by environment rules on the bond graph
   type(lj_typed_type), parameter :: lj_oplsaa = lj_typed_type(set=lj_set_oplsaa)

contains

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function lj_typed_name(self) result(name)
      !> Term
      class(lj_typed_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "lj_typed"
   end function lj_typed_name

   !> Type the structure and look up the parameters
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side, coordinates in bohr
   !> @param[out]    error       unknown force field, atomic number outside the tables or parameter failure
   !> @param[in]     solvent_id  ignored; typing reads only the structure
   subroutine lj_typed_build(self, mol, error, solvent_id)
      !> Term
      class(lj_typed_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Ignored
      integer, intent(in), optional :: solvent_id
      character(len=16) :: label

      if (allocated(self%pair)) deallocate (self%pair)
      if (allocated(self%atomtype)) deallocate (self%atomtype)
      if (allocated(self%parameter_key)) deallocate (self%parameter_key)
      select case (self%set)
      case (lj_set_gaff)
         call build_gaff(self, mol, error)
      case (lj_set_oplsaa)
         call build_oplsaa(self, mol, error)
      case default
         write (label, "(i0)") self%set
         call fatal_error(error, "Unknown typed Lennard-Jones force field "//trim(label)// &
            & "; use lj_gaff or lj_oplsaa")
      end select
   end subroutine lj_typed_build

   !> Type GAFF atoms and look up GAFF2/TM/UFF parameters
   !>
   !> @param[in,out] self   term, pair and labels unallocated
   !> @param[in]     mol    structure of this side, coordinates in bohr
   !> @param[out]    error  atomic number outside the tables or parameter failure
   subroutine build_gaff(self, mol, error)
      !> Term
      type(lj_typed_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=2), allocatable :: gaff(:)
      character(len=atom_label_len) :: label
      real(wp), allocatable :: sigma(:), epsilon(:)
      logical, allocatable :: covered(:)
      logical :: found
      integer :: iat, z

      if (mol%nat < 1) then
         call fatal_error(error, "GAFF Lennard-Jones typing needs at least one atom")
         return
      end if
      do iat = 1, mol%nat
         z = mol%num(mol%id(iat))
         if (z < 1 .or. z > max_elem) then
            call atom_label(mol, iat, label)
            call fatal_error(error, "GAFF Lennard-Jones typing cannot handle atom "//trim(label))
            return
         end if
      end do

      call generate_gaff_atomtypes(mol, gaff)

      allocate (sigma(mol%nat), epsilon(mol%nat), covered(mol%nat), self%atomtype(mol%nat))
      do iat = 1, mol%nat
         call lookup_gaff_lj(gaff(iat), sigma(iat), epsilon(iat), found)
         covered(iat) = found .and. sigma(iat) > 0.0_wp .and. epsilon(iat) > 0.0_wp
         if (covered(iat)) then
            self%atomtype(iat) = trim(gaff(iat))//"/GAFF"
         else
            self%atomtype(iat) = gaff(iat)
         end if
      end do

      allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error, covered, self%atomtype)
      if (allocated(error)) deallocate (self%pair)
   end subroutine build_gaff

   !> Type the bond graph and look up compatible imported OPLS-AA LJ rows
   !>
   !> - Unresolved environment or missing parameters: blank key and label, atoms left for later terms
   !>
   !> @param[in,out] self   term, pair and labels unallocated
   !> @param[in]     mol    explicit hydrogens, bonds and bond orders in atom order
   !> @param[out]    error  invalid graph, ambiguous rules or inconsistent parameters
   subroutine build_oplsaa(self, mol, error)
      !> Term
      type(lj_typed_type), intent(inout) :: self
      !> Molecular graph
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=oplsaa_type_len), allocatable :: labels(:)
      real(wp), allocatable :: sigma(:), epsilon(:)
      logical, allocatable :: covered(:)
      logical :: found
      integer :: iat, number

      call generate_oplsaa_atomtypes(mol, self%parameter_key, error, labels)
      if (allocated(error)) return
      allocate (sigma(mol%nat), epsilon(mol%nat), covered(mol%nat))
      sigma = 0.0_wp
      epsilon = 0.0_wp
      do iat = 1, mol%nat
         covered(iat) = len_trim(self%parameter_key(iat)) > 0
         if (.not. covered(iat)) cycle
         call lookup_oplsaa_lj(self%parameter_key(iat), sigma(iat), epsilon(iat), found, number)
         if (.not. found .or. number /= mol%num(mol%id(iat))) then
            call fatal_error(error, "Inconsistent OPLS-AA LJ mapping for "//trim(labels(iat)))
            deallocate (self%parameter_key)
            return
         end if
      end do
      allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error, covered, labels)
      if (allocated(error)) then
         deallocate (self%pair, self%parameter_key)
         return
      end if
      self%atomtype = labels
   end subroutine build_oplsaa

end module moist_model_moz_potential_lj_typed
