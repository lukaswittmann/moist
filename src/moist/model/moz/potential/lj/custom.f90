!> Lennard-Jones potential term with explicit per-atom parameters
!>
!> Sigma (bohr) and epsilon (Hartree) given per atom of the structure, in its
!> atom order; an optional coverage mask leaves atoms for a later term of the
!> potential set. The parameters are validated and stored at construction,
!> `build` only checks them against the structure
module moist_model_moz_potential_lj_custom
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_model_moz_potential_base, only: potential_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   implicit none(type, external)
   private

   public :: custom_lj_type, new_custom_lj

   !> "Type/Src" label of every atom with explicit parameters
   character(len=*), parameter :: custom_label = "custom"

   !> Lennard-Jones term with explicit parameters
   type, extends(potential_type) :: custom_lj_type
   contains
      procedure :: name => custom_lj_name
      procedure :: build => custom_lj_build
   end type custom_lj_type

contains

   !> Create a term from explicit per-atom parameters
   !>
   !> @param[out] self Term
   !> @param[in] sigma Sigma per atom, bohr
   !> @param[in] epsilon Epsilon per atom, Hartree
   !> @param[out] error Size mismatch, invalid covered parameter or no covered atom
   !> @param[in] covered Optional coverage mask; every atom when absent
   subroutine new_custom_lj(self, sigma, epsilon, error, covered)
      !> Term
      type(custom_lj_type), intent(out) :: self
      !> Sigma per atom, bohr
      real(wp), intent(in) :: sigma(:)
      !> Epsilon per atom, Hartree
      real(wp), intent(in) :: epsilon(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional coverage mask; every atom when absent
      logical, intent(in), optional :: covered(:)
      integer :: iat
      allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error, covered, [(custom_label, iat=1, size(sigma))])
      if (allocated(error)) then
         deallocate (self%pair)
         return
      end if
      if (.not. any(self%pair%covered)) then
         deallocate (self%pair)
         call fatal_error(error, "Custom Lennard-Jones term covers no atom")
      end if
   end subroutine new_custom_lj

   !> Diagnostic name
   !>
   !> @param[in] self Term
   pure function custom_lj_name(self) result(name)
      !> Term
      class(custom_lj_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "custom_lj"
   end function custom_lj_name

   !> Check the stored parameters against the structure
   !>
   !> @param[in,out] self Term
   !> @param[in] mol Structure of this side
   !> @param[out] error Unconstructed term or parameter count differing from the atom count
   !> @param[in] solvent_id Ignored; the parameters are given per atom
   subroutine custom_lj_build(self, mol, error, solvent_id)
      !> Term
      class(custom_lj_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Ignored
      integer, intent(in), optional :: solvent_id
      character(len=24) :: label
      if (.not. allocated(self%pair)) then
         call fatal_error(error, "Custom Lennard-Jones term has no parameters; construct it with new_custom_lj")
         return
      end if
      if (size(self%pair%sigma) /= mol%nat) then
         write (label, '(i0, "/", i0)') size(self%pair%sigma), mol%nat
         call fatal_error(error, "Custom Lennard-Jones parameters do not match the atoms of the structure "// &
            & "(parameters/atoms "//trim(label)//")")
      end if
   end subroutine custom_lj_build

end module moist_model_moz_potential_lj_custom
