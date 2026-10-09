!> Host point multipoles of the MOZ potential layer
!>
!> - Per-atom host charges, dipoles and quadrupoles through the `atomic_multipole` request
!> - Charges as tail charges
!> - Stub; refusal at `update`
module moist_model_moz_potential_electrostatic_multipole
   use mctc_env, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len
   implicit none(type, external)
   private

   public :: multipole_type

   !> Potential term of host point multipoles
   type, extends(potential_term_type) :: multipole_type
   contains
      procedure :: name => multipole_name
      procedure :: update => multipole_update
      !> Host-fed by construction, so a solvent potential refuses it at `add`
      procedure :: is_host_fed => multipole_is_host_fed
   end type multipole_type

contains

   !> Host-fed by construction
   !>
   !> @param[in] self  term
   pure function multipole_is_host_fed(self) result(host_fed)
      !> Term
      class(multipole_type), intent(in) :: self
      !> Host-fed flag
      logical :: host_fed
      host_fed = .true.
   end function multipole_is_host_fed

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function multipole_name(self) result(name)
      !> Term
      class(multipole_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "multipole"
   end function multipole_name

   !> Refuse construction of the host multipole stub
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         structure of this side
   !> @param[out]    error       always set
   !> @param[in]     solvent_id  unused
   subroutine multipole_update(self, mol, error, solvent_id)
      !> Term
      class(multipole_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Unused
      integer, intent(in), optional :: solvent_id
      call fatal_error(error, "Potential term 'multipole' is not implemented yet")
   end subroutine multipole_update

end module moist_model_moz_potential_electrostatic_multipole
