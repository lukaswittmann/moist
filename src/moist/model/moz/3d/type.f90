!> Three-dimensional MOZ model state
!>
!> - Owned spatial grid copy
!> - Solute potential `potential` (public) and a private copy of the updated bulk solvent
!> - Correlation theory pending; the tables come from `potential%compute(grid, ..., solvent=...)`
!> - Coupling scope of solute term i: i
!> - Point grids only, by design: Gaussian widths give the tables nothing and cost more
!> - The model does not watch its potential: update it again after a change
module moist_model_moz_3d_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_begin_registration, moist_phase_energy, &
      & moist_phase_response, moist_phase_gradient, coupling_snapshot, coupling_check_mandatory
   use moist_channels_response, only: response_type
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_solvent_vv, only: solvent_vv_type
   implicit none(type, external)
   private
   public :: model_moz_3d_type, new_moz_3d_model

   !> Refusal of a grid with Gaussian widths
   character(len=*), parameter :: gaussian_grid_refused = &
      & "MOZ 3D model supports point grids only; Gaussian-width grids are not supported yet"

   !> Direct 3D MOZ model; correlation theory is pending
   type, extends(solvation_model_type) :: model_moz_3d_type
      !> Solute potential; add terms before `update`
      type(moz_potential_type) :: potential
      !> Owned spatial grid
      class(moist_math_grid_3d_type), allocatable :: grid
      !> Copy of the updated bulk solvent
      type(solvent_vv_type), allocatable, private :: solvent
   contains
      final :: destroy_moz_3d_model
      procedure :: update => moz_3d_update
      procedure :: get_energy => moz_3d_energy
      procedure :: get_response => moz_3d_response
      procedure :: get_gradient => moz_3d_gradient
      procedure :: atom_count => moz_3d_atom_count
      !> Declare the coupling requests of every solute term in its own scope
      procedure :: declare_pass => moz_3d_declare_pass
      !> Publish the grid geometry as named fields
      procedure :: list_fields => moz_3d_list_fields
      !> Whether the model holds a solvent
      procedure :: has_solvent => moz_3d_has_solvent
      !> Number of solvent sites; 0 before construction
      procedure :: solvent_nsite => moz_3d_solvent_nsite
   end type model_moz_3d_type

contains

   !> Copy the configured spatial grid and an updated bulk solvent
   !>
   !> - Gaussian widths already present: refusal at construction
   !> - Geometry-dependent Gaussian widths: refusal at `update`
   !>
   !> @param[out] self     model
   !> @param[in]  ctx      run context, borrowed for the lifetime of the model
   !> @param[in]  grid     spatial point grid template
   !> @param[in]  solvent  bulk solvent with an updated potential, copied into the model
   !> @param[out] error    gaussian-width grid, solvent not updated or allocation failure
   subroutine new_moz_3d_model(self, ctx, grid, solvent, error)
      !> Model
      type(model_moz_3d_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the model
      type(moist_context_type), intent(in), target :: ctx
      !> Spatial point grid template, copied into the model
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Bulk solvent with an updated potential
      type(solvent_vv_type), intent(in) :: solvent
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Allocation status of the copies
      integer :: stat
      if (allocated(grid%xi0)) then
         call fatal_error(error, gaussian_grid_refused)
         return
      end if
      if (.not. solvent%potential%is_updated()) then
         call fatal_error(error, "Update the 1D VV solvent before constructing the 3D MOZ model")
         return
      end if
      ! TODO: Require solved solvent after VV solve implementation
      ! TODO: share the solvent check and copy with the 1D model
      allocate (self%solvent, source=solvent, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy the 1D VV solvent")
         return
      end if
      self%ctx => ctx
      allocate (self%grid, source=grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy 3D MOZ grid")
         return
      end if
      self%grid%nthreads = ctx%get_num_threads()
   end subroutine new_moz_3d_model

   !> Update the owned grid geometry, then the solute potential on it; print it at verbosity 2
   !>
   !> - Grid-dependent term data (EC fit) built from the updated grid
   !>
   !> @param[in,out] self   model
   !> @param[in]     mol    solute structure
   !> @param[out]    error  grid error, Gaussian-width grid, empty structure, no terms or term failure
   subroutine moz_3d_update(self, mol, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Solute structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
      if (.not. allocated(self%grid) .or. .not. allocated(self%solvent)) then
         call fatal_error(error, "Construct the 3D MOZ model before update")
         return
      end if
      call self%grid%update(mol, error)
      if (allocated(error)) return
      if (allocated(self%grid%xi0)) then
         call fatal_error(error, gaussian_grid_refused)
         return
      end if
      call self%potential%update(mol, self%grid, error)
      if (allocated(error)) return
      if (self%ctx%writes(2)) then
         call self%ctx%message("Solute potential parameters", 2)
         call self%potential%print_table(self%ctx%unit, error)
         if (allocated(error)) return
         call self%ctx%message("", 2)
      end if
      self%updated = .true.
   end subroutine moz_3d_update

   !> Site count of the latest valid structure, from the solute potential
   !>
   !> @param[in] self  model
   function moz_3d_atom_count(self) result(natom)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Site count
      integer :: natom
      natom = self%potential%natom()
   end function moz_3d_atom_count

   !> Declare the coupling requests of every solute term in its own scope
   !>
   !> - Coupling scope of term i: i
   !> - No requests from parameter terms
   !> - EC term declares point-potential probes at grid points
   !>
   !> @param[in,out] self      updated model
   !> @param[in,out] coupling  coupling to declare
   !> @param[out]    error     failed declaration
   subroutine moz_3d_declare_pass(self, coupling, error)
      !> Updated model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call coupling_begin_registration(coupling)
      call self%potential%declare(self%grid, coupling, error)
      if (allocated(error)) return
      call coupling_snapshot(coupling, ngrid=self%grid%ngrid, natom=self%grid%natom)
   end subroutine moz_3d_declare_pass

   !> Whether the model holds a solvent
   !>
   !> @param[in] self  model
   pure function moz_3d_has_solvent(self) result(has)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Presence
      logical :: has
      has = allocated(self%solvent)
   end function moz_3d_has_solvent

   !> Number of solvent sites; 0 before construction
   !>
   !> @param[in] self  model
   pure function moz_3d_solvent_nsite(self) result(nsite)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Site count
      integer :: nsite
      nsite = 0
      if (allocated(self%solvent)) nsite = self%solvent%nsite()
   end function moz_3d_solvent_nsite

   !> Publish the grid geometry as named fields, the model's evaluation domain
   !>
   !> @param[in]     self   model
   !> @param[in,out] query  field walker
   subroutine moz_3d_list_fields(self, query)
      !> Model
      class(model_moz_3d_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
      if (.not. allocated(self%grid)) return
      call query%add_int_value("ngrid", "Number of volume grid points", self%grid%ngrid)
      call query%add_int_value("natom", "Number of solute atoms", self%grid%natom)
      call query%add_real2("xyz", "Grid coordinates, bohr (3, ngrid)", self%grid%xyz)
      call query%add_real("w", "Integration weights, bohr**3 (ngrid)", self%grid%w)
      call query%add_int("owner", "Atom owner, zero based; -1 for unowned points", &
         & self%grid%owner, zero_based=.true.)
   end subroutine moz_3d_list_fields

   !> Report pending 3D MOZ energy theory
   !>
   !> @param[in,out] self      model
   !> @param[in,out] coupling  coupling owned by the model
   !> @param[in,out] energy    energy, untouched
   !> @param[out]    error     foreign coupling, model not updated, missing answers or pending theory
   subroutine moz_3d_energy(self, coupling, energy, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling owned by the model
      class(coupling_type), intent(inout), target :: coupling
      !> Energy, untouched
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, error)
      if (allocated(error)) return
      call fatal_error(error, "3D MOZ energy is not implemented")
   end subroutine moz_3d_energy

   !> Report pending 3D MOZ response theory
   !>
   !> @param[in,out] self      model
   !> @param[in,out] coupling  coupling owned by the model
   !> @param[in,out] response  response, untouched
   !> @param[out]    error     foreign coupling, model not updated, missing answers or pending theory
   subroutine moz_3d_response(self, coupling, response, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling owned by the model
      class(coupling_type), intent(inout), target :: coupling
      !> Response, untouched
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_response, error)
      if (allocated(error)) return
      call fatal_error(error, "3D MOZ response is not implemented")
   end subroutine moz_3d_response

   !> Report pending 3D MOZ gradient theory
   !>
   !> @param[in,out] self      model
   !> @param[in,out] coupling  coupling owned by the model
   !> @param[in,out] response  response, untouched
   !> @param[in,out] gradient  nuclear gradient, untouched
   !> @param[out]    error     foreign coupling, model not updated, missing answers or pending theory
   subroutine moz_3d_gradient(self, coupling, response, gradient, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling owned by the model
      class(coupling_type), intent(inout), target :: coupling
      !> Response, untouched
      type(response_type), intent(inout) :: response
      !> Nuclear gradient, untouched
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, error)
      if (allocated(error)) return
      call fatal_error(error, "3D MOZ gradient is not implemented")
   end subroutine moz_3d_gradient

   !> Release the coupling registry
   !>
   !> @param[in,out] self  model being destroyed
   subroutine destroy_moz_3d_model(self)
      !> Model being destroyed
      type(model_moz_3d_type), intent(inout) :: self
      call self%clear_couplings()
   end subroutine destroy_moz_3d_model

end module moist_model_moz_3d_type
