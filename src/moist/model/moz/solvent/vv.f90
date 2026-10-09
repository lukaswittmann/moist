!> Bulk solvent of the MOZ models
!>
!> Not a solvation model: holds the solvent system (solvent id, temperature,
!> number density, structure) and the solvent potential, which refuses
!> host-data terms and carries the mixing rule and Ng split of the VV
!> tables. One site per atom. A custom solvent is a hand-filled
!> `solvation_system_type` (solvent id 0, temperature, number density and
!> structure). Will own the 1D VV solve and the solvent susceptibility; the UV
!> models copy it
!>
!>    call solvent%potential%add(lj_spce, error)
!>    call solvent%potential%add(coulomb_spce, error)
!>    call solvent%update(error)
!>    call solvent%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, error)
!>
!> The solvent does not watch its potential: update it again after a change
module moist_model_moz_solvent_vv
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_data_solvents, only: solvation_system_type
   use moist_model_moz_potential, only: moz_potential_type
   implicit none(type, external)
   private

   public :: solvent_vv_type, new_vv_solvent

   !> Bulk solvent: system and potential
   type :: solvent_vv_type
      !> Run context, borrowed for the lifetime of the solvent
      type(moist_context_type), pointer :: ctx => null()
      !> Solvent system: id, temperature, number density and structure (bohr)
      type(solvation_system_type), private :: system
      !> Solvent potential; refuses host-data terms. Second side of the UV tables
      type(moz_potential_type) :: potential
      !> Whether the solvent was constructed
      logical, private :: constructed = .false.
      !> Whether the VV solution is available
      logical, private :: solved = .false.
   contains
      !> Update the potential from the solvent structure
      procedure :: update => vv_update
      !> Solve the 1D VV equations; pending
      procedure :: solve => vv_solve
      !> Whether the VV solution is available
      procedure :: is_solved => vv_is_solved
      !> Number of solvent sites
      procedure :: nsite => vv_nsite
      !> Copy of the solvent structure
      procedure :: structure => vv_structure
      !> Solvent id; 0 for a custom solvent
      procedure :: solvent_id => vv_solvent_id
   end type solvent_vv_type

contains

   !> Construct a bulk solvent from a solvation system
   !>
   !> Reads the solvent id, temperature, number density
   !> (`solvent_number_density_au`) and structure (`solv_mol`, bohr) of the
   !> system; a table entry comes from `new_solvation_system`
   !>
   !> @param[out] self Solvent
   !> @param[in] ctx Run context, borrowed for the lifetime of the solvent
   !> @param[in] system Solvent system
   !> @param[out] error Missing structure or invalid state point
   subroutine new_vv_solvent(self, ctx, system, error)
      !> Solvent
      type(solvent_vv_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the solvent
      type(moist_context_type), intent(in), target :: ctx
      !> Solvent system
      type(solvation_system_type), intent(in) :: system
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. allocated(system%solv_mol)) then
         call fatal_error(error, "1D VV solvent system has no structure (solv_mol)")
         return
      end if
      if (system%solv_mol%nat < 1) then
         call fatal_error(error, "1D VV solvent requires at least one atom")
         return
      end if
      if (.not. ieee_is_finite(system%temperature)) then
         call fatal_error(error, "1D VV solvent temperature must be finite")
         return
      end if
      if (system%temperature <= 0.0_wp) then
         call fatal_error(error, "1D VV solvent temperature must be positive")
         return
      end if
      if (.not. ieee_is_finite(system%solvent_number_density_au)) then
         call fatal_error(error, "1D VV solvent number density must be finite")
         return
      end if
      if (system%solvent_number_density_au <= 0.0_wp) then
         call fatal_error(error, "1D VV solvent number density must be positive")
         return
      end if
      self%system = system
      call self%potential%refuse_host_fed_terms()
      self%ctx => ctx
      self%constructed = .true.
   end subroutine new_vv_solvent

   !> Update the potential from the solvent structure
   !>
   !> Terms keyed by solvent receive the solvent id; none for a custom solvent
   !> (id 0). Drops any VV solution
   !>
   !> @param[in,out] self Constructed solvent
   !> @param[out] error Unconstructed solvent, no terms or term failure
   subroutine vv_update(self, error)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%solved = .false.
      if (.not. self%constructed) then
         call fatal_error(error, "Construct the 1D VV solvent before update")
         return
      end if
      if (self%system%solvent_id > 0) then
         call self%potential%update(self%system%solv_mol, error, self%system%solvent_id)
      else
         call self%potential%update(self%system%solv_mol, error)
      end if
      if (allocated(error)) return
      if (.not. associated(self%ctx)) return
      if (self%ctx%writes(2)) then
         call self%ctx%message("Solvent potential parameters", 2)
         call self%potential%print_table(self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
   end subroutine vv_update

   !> Solve the 1D VV equations; pending
   !>
   !> @param[in,out] self Built solvent
   !> @param[out] error Always set
   subroutine vv_solve(self, error)
      !> Solvent
      class(solvent_vv_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%solved = .false.
      call fatal_error(error, "1D VV solve is pending")
   end subroutine vv_solve

   !> Whether the VV solution is available
   !>
   !> @param[in] self Solvent
   pure function vv_is_solved(self) result(solved)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solution flag
      logical :: solved
      solved = self%solved
   end function vv_is_solved

   !> Number of solvent sites, one per atom; 0 before construction
   !>
   !> @param[in] self Solvent
   pure function vv_nsite(self) result(nsite)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Site count
      integer :: nsite
      nsite = 0
      if (self%constructed) nsite = self%system%solv_mol%nat
   end function vv_nsite

   !> Copy of the solvent structure, bohr
   !>
   !> @param[in] self Constructed solvent
   function vv_structure(self) result(mol)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solvent structure; empty before construction
      type(structure_type) :: mol
      if (self%constructed) mol = self%system%solv_mol
   end function vv_structure

   !> Solvent id; 0 for a custom solvent or before construction
   !>
   !> @param[in] self Solvent
   pure function vv_solvent_id(self) result(id)
      !> Solvent
      class(solvent_vv_type), intent(in) :: self
      !> Solvent id
      integer :: id
      id = 0
      if (self%constructed) id = self%system%solvent_id
   end function vv_solvent_id

end module moist_model_moz_solvent_vv
