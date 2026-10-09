!> One-dimensional MOZ model state
!>
!> - Owned radial grid pair `rgrid`, `kgrid`, independent of the solute structure
!> - Solute potential `potential` (public) and a private copy of the updated bulk solvent
!> - Radial theory pending; the tables come from `potential%compute(rgrid, kgrid, ..., solvent=...)`
!> - Coupling scope of solute term i: i
!> - Monopole term owns its host-charge requirements
!> - The model does not watch its potential: update it again after a change
module moist_model_moz_1d_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_begin_registration, moist_phase_energy, &
      & moist_phase_response, moist_phase_gradient, coupling_snapshot, coupling_check_mandatory
   use moist_channels_response, only: response_type
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_kernel_table_1d, only: require_radial_grids
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   implicit none(type, external)
   private
   public :: model_moz_1d_type, new_moz_1d_model

   !> Direct 1D MOZ model; radial theory is pending
   type, extends(solvation_model_type) :: model_moz_1d_type
      !> Solute potential; add terms before `update`
      type(moz_potential_type) :: potential
      !> Owned real-space radial grid, bohr
      type(moist_math_grid_radial_type) :: rgrid
      !> Owned reciprocal radial grid, 1/bohr
      type(moist_math_grid_radial_type) :: kgrid
      !> Copy of the updated bulk solvent
      type(solvent_vv_type), allocatable, private :: solvent
   contains
      final :: destroy_moz_1d_model
      procedure :: update => moz_1d_update
      procedure :: get_energy => moz_1d_energy
      procedure :: get_response => moz_1d_response
      procedure :: get_gradient => moz_1d_gradient
      procedure :: atom_count => moz_1d_atom_count
      !> Declare the coupling requests of every solute term in its own scope
      procedure :: declare_pass => moz_1d_declare_pass
      !> Publish the radial grid as named fields
      procedure :: list_fields => moz_1d_list_fields
      !> Whether the model holds a solvent
      procedure :: has_solvent => moz_1d_has_solvent
      !> Number of solvent sites; 0 before construction
      procedure :: solvent_nsite => moz_1d_solvent_nsite
   end type model_moz_1d_type

contains

   !> Copy the radial grid pair and an updated bulk solvent
   !>
   !> @param[out] self     model
   !> @param[in]  ctx      run context, borrowed for the lifetime of the model
   !> @param[in]  rgrid    real-space radial grid, nodes r > 0; e.g. from `new_uniform_radial_pair`
   !> @param[in]  kgrid    reciprocal radial grid of the pair
   !> @param[in]  solvent  bulk solvent with an updated potential, copied into the model
   !> @param[out] error    unbuilt grid, node at r <= 0, solvent not updated or copy failure
   subroutine new_moz_1d_model(self, ctx, rgrid, kgrid, solvent, error)
      !> Model
      type(model_moz_1d_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the model
      type(moist_context_type), intent(in), target :: ctx
      !> Real-space radial grid, bohr
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> Reciprocal radial grid, 1/bohr
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Bulk solvent with an updated potential
      type(solvent_vv_type), intent(in) :: solvent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: stat
      call require_radial_grids(rgrid, kgrid, error)
      if (allocated(error)) return
      if (.not. solvent%potential%is_updated()) then
         call fatal_error(error, "Update the 1D VV solvent before constructing the 1D MOZ model")
         return
      end if
      ! TODO: Require solved solvent after VV solve implementation
      allocate (self%solvent, source=solvent, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy the 1D VV solvent")
         return
      end if
      self%rgrid = rgrid
      self%kgrid = kgrid
      self%ctx => ctx
   end subroutine new_moz_1d_model

   !> Update the solute potential from the structure on the radial grid and print it at verbosity 2
   !>
   !> @param[in,out] self   model
   !> @param[in]     mol    solute structure
   !> @param[out]    error  unconstructed model, empty structure, no terms, term failure or a term refusing the radial grid
   subroutine moz_1d_update(self, mol, error)
      !> Model
      class(model_moz_1d_type), intent(inout) :: self
      !> Solute structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
      if (.not. associated(self%ctx) .or. .not. allocated(self%solvent)) then
         call fatal_error(error, "Construct the 1D MOZ model before update")
         return
      end if
      call self%potential%update(mol, self%rgrid, error)
      if (allocated(error)) return
      if (self%ctx%writes(2)) then
         call self%ctx%message("Solute potential parameters", 2)
         call self%potential%print_table(self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
      self%updated = .true.
   end subroutine moz_1d_update

   !> Site count of the latest valid structure, from the solute potential
   !>
   !> @param[in] self  model
   function moz_1d_atom_count(self) result(natom)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Site count
      integer :: natom
      natom = self%potential%natom()
   end function moz_1d_atom_count

   !> Declare the coupling requests of every solute term in its own scope
   !>
   !> - Coupling scope of term i: i
   !> - No requests from parameter terms
   !> - Monopole term declares host charges for tables
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling to declare
   !> @param[out]    error     failed declaration
   subroutine moz_1d_declare_pass(self, coupling, error)
      !> Updated model
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call coupling_begin_registration(coupling)
      call self%potential%declare(self%rgrid, coupling, error)
      if (allocated(error)) return
      call coupling_snapshot(coupling, ngrid=self%rgrid%npts, natom=self%potential%natom())
   end subroutine moz_1d_declare_pass

   !> Publish the radial grid as named fields, the model's evaluation domain
   !>
   !> - Every solute site carries the same radii
   !>
   !> @param[in]     self   model
   !> @param[in,out] query  field walker
   subroutine moz_1d_list_fields(self, query)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
      if (.not. allocated(self%rgrid%r)) return
      call query%add_int_value("ngrid", "Number of radial grid points", self%rgrid%npts)
      call query%add_int_value("natom", "Number of solute atoms", self%potential%natom())
      call query%add_real("r", "Radii from each solute site, bohr (ngrid)", self%rgrid%r)
   end subroutine moz_1d_list_fields

   !> Whether the model holds a solvent
   !>
   !> @param[in] self  model
   pure function moz_1d_has_solvent(self) result(has)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Presence
      logical :: has
      has = allocated(self%solvent)
   end function moz_1d_has_solvent

   !> Number of solvent sites; 0 before construction
   !>
   !> @param[in] self  model
   pure function moz_1d_solvent_nsite(self) result(nsite)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Site count
      integer :: nsite
      nsite = 0
      if (allocated(self%solvent)) nsite = self%solvent%nsite()
   end function moz_1d_solvent_nsite

   !> Report pending 1D MOZ energy theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] energy Energy
   !> @param[out] error Error
   subroutine moz_1d_energy(self, coupling, energy, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Energy
      real(wp), intent(inout) :: energy
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, error)
      if (allocated(error)) return
      call fatal_error(error, "1D MOZ energy is not implemented")
   end subroutine moz_1d_energy

   !> Report pending 1D MOZ response theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[out] error Error
   subroutine moz_1d_response(self, coupling, response, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_response, error)
      if (allocated(error)) return
      call fatal_error(error, "1D MOZ response is not implemented")
   end subroutine moz_1d_response

   !> Report pending 1D MOZ gradient theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[in,out] gradient Gradient
   !> @param[out] error Error
   subroutine moz_1d_gradient(self, coupling, response, gradient, error)
      !> Self
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling
      class(coupling_type), intent(inout), target :: coupling
      !> Response
      type(response_type), intent(inout) :: response
      !> Gradient
      real(wp), intent(inout) :: gradient(:, :)
      !> Error
      type(error_type), allocatable, intent(out) :: error
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, error)
      if (allocated(error)) return
      call fatal_error(error, "1D MOZ gradient is not implemented")
   end subroutine moz_1d_gradient

   !> Release the coupling registry
   !>
   !> @param[in,out] self Model being destroyed
   subroutine destroy_moz_1d_model(self)
      !> Model being destroyed
      type(model_moz_1d_type), intent(inout) :: self
      call self%clear_couplings()
   end subroutine destroy_moz_1d_model

end module moist_model_moz_1d_type
