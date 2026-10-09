!> Lennard-Jones potential term for an explicitly selected solvent model
!>
!> `lj_set_spce` preserves the legacy SPC/E-style water parameters, including
!> nonzero hydrogen LJ. The model requires one neutral, closed-shell H2O
!> molecule and assigns sites by element, independently of GAFF typing
!>
!> A positive solvent ID must name water; an absent or zero ID permits a
!> custom water structure. Charges and the mixing rule are selected separately
module moist_model_moz_potential_lj_solvent
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_model_moz_potential_base, only: potential_type, potential_name_len
   use moist_model_moz_potential_lj_base, only: new_lj_12_6
   use moist_model_moz_potential_lj_spce_parameters, only: lookup_spce_lj
   use moist_model_moz_potential_utils, only: is_water
   implicit none(type, external)
   private

   public :: lj_solvent_type, lj_set_spce, lj_spce

   !> Legacy SPC/E-style water LJ model
   integer, parameter :: lj_set_spce = 1

   !> Lennard-Jones term by solvent model
   type, extends(potential_type) :: lj_solvent_type
      !> Solvent model, `lj_set_spce`; 0 until constructed
      integer :: set = 0
      !> Water site per atom after build, "ow/SPCE" or "hw/SPCE"
      character(len=8), allocatable :: atomtype(:)
   contains
      procedure :: name => lj_solvent_name
      procedure :: build => lj_solvent_build
   end type lj_solvent_type

   !> SPC/E-style water sites
   type(lj_solvent_type), parameter :: lj_spce = lj_solvent_type(set=lj_set_spce)

contains

   !> Diagnostic name
   !>
   !> @param[in] self Term
   pure function lj_solvent_name(self) result(name)
      !> Term
      class(lj_solvent_type), intent(in) :: self
      !> Name
      character(len=potential_name_len) :: name
      name = "lj_solvent"
   end function lj_solvent_name

   !> Validate the solvent and assign LJ parameters to its sites
   !>
   !> @param[in,out] self Term
   !> @param[in] mol One neutral, closed-shell water molecule, in any atom order
   !> @param[out] error Unconstructed term, incompatible solvent or parameter failure
   !> @param[in] solvent_id Water ID, or zero for custom water; optional
   subroutine lj_solvent_build(self, mol, error, solvent_id)
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

      if (allocated(self%pair)) deallocate (self%pair)
      if (allocated(self%atomtype)) deallocate (self%atomtype)
      if (self%set /= lj_set_spce) then
         write (label, "(i0)") self%set
         call fatal_error(error, "Unknown Lennard-Jones solvent model "//trim(label)// &
            & "; use lj_spce")
         return
      end if
      if (.not. allocated(mol%id) .or. .not. allocated(mol%num)) then
         call fatal_error(error, "SPC/E-style Lennard-Jones solvent model needs initialized atom identities")
         return
      end if
      if (size(mol%id) /= mol%nat .or. any(mol%id < 1) .or. any(mol%id > size(mol%num))) then
         call fatal_error(error, "SPC/E-style Lennard-Jones solvent model has invalid atom IDs")
         return
      end if
      ! Anything but neutral, closed-shell H2O is left uncovered for a later term
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
         call lookup_spce_lj(numbers(iat), sigma(iat), epsilon(iat), covered(iat))
         if (numbers(iat) == 8) then
            self%atomtype(iat) = "ow/SPCE"
         else
            self%atomtype(iat) = "hw/SPCE"
         end if
      end do
      allocate (self%pair)
      call new_lj_12_6(self%pair, sigma, epsilon, error, covered, self%atomtype)
      if (allocated(error)) then
         deallocate (self%pair, self%atomtype)
         return
      end if
   end subroutine lj_solvent_build

end module moist_model_moz_potential_lj_solvent
