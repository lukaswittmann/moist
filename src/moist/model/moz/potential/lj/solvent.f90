!> Lennard-Jones potential term for an explicitly selected solvent model
!>
!> - `lj_set_spce`: standard SPC/E, with inactive hydrogen LJ
!> - `lj_set_cspce`: cSPC/E, with nonzero hydrogen LJ
!> - One neutral, closed-shell H2O molecule; sites by element, independent of GAFF typing
!> - Positive solvent ID: water required; absent or zero ID: custom water permitted
!> - Separate charge and mixing-rule selection
module moist_model_moz_potential_lj_solvent
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_model_moz_potential_term, only: potential_term_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   use moist_model_moz_potential_lj_spce_parameters, only: lookup_spce_lj
   use moist_model_moz_potential_lj_cspce_parameters, only: lookup_cspce_lj
   use moist_model_moz_solvent_water, only: is_water
   implicit none(type, external)
   private

   public :: lj_solvent_type, lj_set_spce, lj_set_cspce, lj_spce, lj_cspce

   !> Standard SPC/E water LJ model
   integer, parameter :: lj_set_spce = 1
   !> cSPC/E water LJ model, with the Luchko hydrogen modification
   integer, parameter :: lj_set_cspce = 2

   !> Lennard-Jones term by solvent model
   type, extends(potential_term_type) :: lj_solvent_type
      !> Solvent model, `lj_set_spce` or `lj_set_cspce`; 0 until constructed
      integer :: set = 0
      !> Water site per atom after update, ow/hw followed by SPCE or CSPCE
      character(len=8), allocatable :: atomtype(:)
   contains
      procedure :: name => lj_solvent_name
      procedure :: update => lj_solvent_update
   end type lj_solvent_type

   ! ifx workaround: protected, not parameter; ifx 2025.1-2025.3 fails copies of derived-type
   ! constants with allocatable components (forrtl 189); TODO: revert to parameter once ifx is fixed
   !> Standard SPC/E water sites
   type(lj_solvent_type), protected :: lj_spce = lj_solvent_type(set=lj_set_spce)
   !> cSPC/E water sites
   type(lj_solvent_type), protected :: lj_cspce = lj_solvent_type(set=lj_set_cspce)

contains

   !> Diagnostic name
   !>
   !> @param[in] self  term
   pure function lj_solvent_name(self) result(name)
      !> Term
      class(lj_solvent_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "lj_solvent"
   end function lj_solvent_name

   !> Validate the solvent and assign LJ parameters to its sites
   !>
   !> @param[in,out] self        term
   !> @param[in]     mol         one neutral, closed-shell water molecule, in any atom order
   !> @param[out]    error       unconstructed term, incompatible solvent or parameter failure
   !> @param[in]     solvent_id  water ID, or zero for custom water; optional
   subroutine lj_solvent_update(self, mol, error, solvent_id)
      !> Term
      class(lj_solvent_type), intent(inout) :: self
      !> Solvent structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Solvent ID
      integer, intent(in), optional :: solvent_id
      real(wp) :: sigma(3), epsilon(3)
      logical :: covered(3)
      integer :: numbers(3), water_id, iat
      character(len=16) :: label
      character(len=5) :: model_label

      if (allocated(self%pair)) deallocate (self%pair)
      if (allocated(self%atomtype)) deallocate (self%atomtype)
      select case (self%set)
      case (lj_set_spce)
         model_label = "SPCE"
      case (lj_set_cspce)
         model_label = "CSPCE"
      case default
         write (label, "(i0)") self%set
         call fatal_error(error, "Unknown Lennard-Jones solvent model "//trim(label)// &
            & "; use lj_spce or lj_cspce")
         return
      end select
      if (.not. allocated(mol%id) .or. .not. allocated(mol%num)) then
         call fatal_error(error, "SPC/E-style Lennard-Jones solvent model needs initialized atom identities")
         return
      end if
      if (size(mol%id) /= mol%nat .or. any(mol%id < 1) .or. any(mol%id > size(mol%num))) then
         call fatal_error(error, "SPC/E-style Lennard-Jones solvent model has invalid atom IDs")
         return
      end if
      ! Non-water, charged or open-shell structure: atoms left for later terms
      if (.not. is_water(mol, solvent_id)) then
         allocate (self%atomtype(mol%nat), source=repeat(" ", len(self%atomtype)))
         allocate (self%pair)
         call new_lj_12_6(self%pair, [(0.0_wp, iat=1, mol%nat)], [(0.0_wp, iat=1, mol%nat)], error, &
            & [(.false., iat=1, mol%nat)])
         if (allocated(error)) deallocate (self%pair, self%atomtype)
         return
      end if
      numbers = mol%num(mol%id)

      allocate (self%atomtype(mol%nat))
      do iat = 1, mol%nat
         if (self%set == lj_set_spce) then
            call lookup_spce_lj(numbers(iat), sigma(iat), epsilon(iat), covered(iat))
         else
            call lookup_cspce_lj(numbers(iat), sigma(iat), epsilon(iat), covered(iat))
         end if
         if (numbers(iat) == 8) then
            self%atomtype(iat) = "ow/"//trim(model_label)
         else
            self%atomtype(iat) = "hw/"//trim(model_label)
         end if
      end do
      allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error, covered, self%atomtype)
      if (allocated(error)) then
         deallocate (self%pair, self%atomtype)
         return
      end if
   end subroutine lj_solvent_update

end module moist_model_moz_potential_lj_solvent
