!> Lennard-Jones potential term with parameters by element
!>
!> - Sigma and epsilon lookup by atomic number; element table chosen at construction
!> - `lj_set_uff`: UFF, H through Lr (Rappe et al. 1992)
!> - `lj_set_dreiding`: DREIDING explicit-atom rows (Mayo et al. 1990), ordinary hydrogen LJ
!> - Hydrogen-bond term outside the DREIDING LJ term
!> - `lj_set_tm`: fourteen transition metals, former GAFF element fallback
!> - Elements absent from table: atoms left for later terms
!> - Mixing rule chosen on solvent
module moist_model_moz_potential_lj_element
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_symbols, only: to_symbol
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   use moist_model_moz_potential_lj_uff_parameters, only: lookup_uff_lj
   use moist_model_moz_potential_lj_dreiding_parameters, only: lookup_dreiding_lj
   use moist_model_moz_potential_lj_tm_parameters, only: lookup_tm_lj
   implicit none(type, external)
   private

   public :: lj_element_type, lj_set_uff, lj_set_dreiding, lj_set_tm
   public :: lj_uff, lj_dreiding, lj_tm

   !> UFF element table
   integer, parameter :: lj_set_uff = 1
   !> DREIDING element table
   integer, parameter :: lj_set_dreiding = 2
   !> Transition-metal element table
   integer, parameter :: lj_set_tm = 3

   !> Lennard-Jones term by element
   type, extends(potential_term_type) :: lj_element_type
      !> Element table, `lj_set_uff`, `lj_set_dreiding` or `lj_set_tm`; 0 until constructed
      integer :: set = 0
      !> "El/table" label per atom after update, e.g. "O/UFF"; blank for an uncovered atom
      character(len=12), allocatable :: atomtype(:)
   contains
      procedure :: name => lj_element_name
      procedure :: update => lj_element_update
   end type lj_element_type

   !> UFF element rows
   type(lj_element_type), parameter :: lj_uff = lj_element_type(set=lj_set_uff)
   !> DREIDING element rows
   type(lj_element_type), parameter :: lj_dreiding = lj_element_type(set=lj_set_dreiding)
   !> Transition-metal element rows
   type(lj_element_type), parameter :: lj_tm = lj_element_type(set=lj_set_tm)

contains

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function lj_element_name(self) result(name)
      !> Term
      class(lj_element_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "lj_element"
   end function lj_element_name

   !> Look up the parameters of every atom by element
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       unconstructed term, empty structure or parameter failure
   !> @param[in]     solvent_id  ignored; the lookup uses the elements
   subroutine lj_element_update(self, mol, error, solvent_id)
      !> Term
      class(lj_element_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Ignored
      integer, intent(in), optional :: solvent_id
      real(wp), allocatable :: sigma(:), epsilon(:)
      logical, allocatable :: covered(:)
      character(len=8) :: source
      character(len=16) :: label
      integer :: iat, number

      if (allocated(self%pair)) deallocate (self%pair)
      if (allocated(self%atomtype)) deallocate (self%atomtype)
      select case (self%set)
      case (lj_set_uff)
         source = "UFF"
      case (lj_set_dreiding)
         source = "DREIDING"
      case (lj_set_tm)
         source = "TM"
      case default
         write (label, "(i0)") self%set
         call fatal_error(error, "Unknown Lennard-Jones element table "//trim(label)// &
            & "; use lj_uff, lj_dreiding or lj_tm")
         return
      end select
      if (mol%nat < 1) then
         call fatal_error(error, "Lennard-Jones element term needs at least one atom")
         return
      end if
      allocate (sigma(mol%nat), epsilon(mol%nat), covered(mol%nat), self%atomtype(mol%nat))
      self%atomtype = ""
      do iat = 1, mol%nat
         number = mol%num(mol%id(iat))
         select case (self%set)
         case (lj_set_uff)
            call lookup_uff_lj(number, sigma(iat), epsilon(iat), covered(iat))
         case (lj_set_dreiding)
            call lookup_dreiding_lj(number, sigma(iat), epsilon(iat), covered(iat))
         case default
            call lookup_tm_lj(number, sigma(iat), epsilon(iat), covered(iat))
         end select
         if (covered(iat)) self%atomtype(iat) = trim(to_symbol(number))//"/"//trim(source)
      end do
      allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error, covered, self%atomtype)
      if (allocated(error)) deallocate (self%pair)
   end subroutine lj_element_update

end module moist_model_moz_potential_lj_element
