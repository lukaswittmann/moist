!> Three-dimensional MOZ model state
module moist_model_moz_3d_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, point_potential_request_type, &
      & gaussian_potential_request_type, atomic_charge_request_type, &
      & atomic_multipole_request_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_begin_registration, coupling_set_scope, &
      & coupling_register, request_require, coupling_snapshot, coupling_check_mandatory
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: model_moz_3d_type, new_moz_3d_model

   !> Direct 3D MOZ model; correlation theory is pending
   type, extends(solvation_model_type) :: model_moz_3d_type
      !> Owned spatial grid
      class(moist_math_grid_3d_type), allocatable :: grid
      ! TODO: replace the coupling_mode string with a typed selector set at construction
      !> Electrostatic coupling source: "qat", "ec" or "multipoles" (3D only)
      character(len=16) :: coupling_mode = "ec"
   contains
      final :: destroy_moz_3d_model
      procedure :: update => moz_3d_update
      procedure :: get_energy => moz_3d_energy
      procedure :: get_response => moz_3d_response
      procedure :: get_gradient => moz_3d_gradient
      procedure :: atom_count => moz_3d_atom_count
      !> Declare the requests of the electrostatic coupling source: point
      !> charges ("qat"), point multipoles ("multipoles") or the grid potential
      !> plus tail charges ("ec")
      procedure :: declare_pass => moz_3d_declare_pass
      !> Publish the grid geometry as named fields
      procedure :: list_fields => moz_3d_list_fields
   end type model_moz_3d_type

contains

   !> Copy the configured spatial grid
   !>
   !> @param[out] self Model
   !> @param[in] grid Spatial grid template
   !> @param[in] ctx Run context, borrowed for the lifetime of the model
   !> @param[out] error Allocation failure
   subroutine new_moz_3d_model(self, grid, ctx, error)
      !> Model
      type(model_moz_3d_type), intent(out) :: self
      !> Spatial grid template, copied into the model
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Run context, borrowed for the lifetime of the model
      type(moist_context_type), intent(in), target :: ctx
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Allocation status of the grid copy
      integer :: stat
      allocate (self%grid, source=grid, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy 3D MOZ grid")
         return
      end if
      self%grid%nthreads = ctx%get_num_threads()
      self%ctx => ctx
   end subroutine new_moz_3d_model

   !> Update owned grid geometry
   !>
   !> @param[in,out] self Model
   !> @param[in] mol Solute structure
   !> @param[out] error Grid error
   subroutine moz_3d_update(self, mol, error)
      !> Model
      class(model_moz_3d_type), intent(inout) :: self
      !> Solute structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
      if (.not. allocated(self%grid)) then
         call fatal_error(error, "Construct the 3D MOZ model before update")
         return
      end if
      call self%grid%update(mol, error)
      if (allocated(error)) return
      self%updated = .true.
   end subroutine moz_3d_update

   !> Number of solute atoms in the grid
   !>
   !> @param[in] self Self
   function moz_3d_atom_count(self) result(natom)
      !> Self
      class(model_moz_3d_type), intent(in) :: self
      !> Natom
      integer :: natom
      natom = 0
      if (allocated(self%grid)) natom = self%grid%natom
   end function moz_3d_atom_count

   !> Declare the requests of the electrostatic coupling source
   !>
   !> - "qat": the solute's partial charges, `q` in every phase
   !> - "multipoles": the solute's point multipoles, `q`, `mu` and `theta` in
   !>   every phase; `q` doubles as the tail charges
   !> - "ec": the host potential on the grid, `phi` in every phase and
   !>   `dphi_dr` for the gradient, plus the partial charges, `q` in every
   !>   phase, as the tail charges of the Ng split
   !>
   !> A Cartesian grid always realizes Gaussian widths; a molecular grid does
   !> so only when configured for it, otherwise it reports bare point values
   !> Either way, `allocated(grid%xi0)` after `update` selects the "ec"
   !> potential: present widths mean a Gaussian-probed potential, which also
   !> needs `dphi_dxi` for the gradient when the widths follow the geometry;
   !> their absence a bare point potential
   !>
   !> @param[in,out] self Updated model
   !> @param[in,out] coupling Coupling to declare
   !> @param[out] error Unknown coupling source or failed declaration
   subroutine moz_3d_declare_pass(self, coupling, error)
      !> Updated model
      class(model_moz_3d_type), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Phases every source requires its outputs in
      integer, parameter :: phases(3) = [moist_phase_energy, moist_phase_response, moist_phase_gradient]
      !> Per-atom partial charge request
      type(atomic_charge_request_type) :: charges
      !> Per-atom point multipole request
      type(atomic_multipole_request_type) :: multipoles
      integer :: ip

      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 0)
      select case (self%coupling_mode)
      case ("qat")
         do ip = 1, size(phases)
            call request_require(charges, phases(ip), "q", error)
            if (allocated(error)) return
         end do
         call coupling_register(coupling, "charges", charges, error)
      case ("multipoles")
         ! TODO: this is what gfn2 and gxtb will use
         do ip = 1, size(phases)
            call request_require(multipoles, phases(ip), "q", error)
            if (allocated(error)) return
            call request_require(multipoles, phases(ip), "mu", error)
            if (allocated(error)) return
            call request_require(multipoles, phases(ip), "theta", error)
            if (allocated(error)) return
         end do
         call coupling_register(coupling, "multipoles", multipoles, error)
      case ("ec")
         ! TODO: I added the gaussian one here; but as the point one works we
         !       should use that (integrals are much faster)
         if (allocated(self%grid%xi0)) then
            block
               !> Gaussian-probed potential request
               type(gaussian_potential_request_type) :: potential
               do ip = 1, size(phases)
                  call request_require(potential, phases(ip), "phi", error)
                  if (allocated(error)) return
               end do
               call request_require(potential, moist_phase_gradient, "dphi_dr", error)
               if (allocated(error)) return
               call request_require(potential, moist_phase_gradient, "dphi_dxi", &
                  & self%grid%has_geometry_dependent_xi0(), error)
               if (allocated(error)) return
               call coupling_register(coupling, "potential", potential, error)
            end block
         else
            block
               !> Bare point potential request
               type(point_potential_request_type) :: potential
               do ip = 1, size(phases)
                  call request_require(potential, phases(ip), "phi", error)
                  if (allocated(error)) return
               end do
               call request_require(potential, moist_phase_gradient, "dphi_dr", error)
               if (allocated(error)) return
               call coupling_register(coupling, "potential", potential, error)
            end block
         end if
         if (allocated(error)) return
         ! The tail charges of the Ng split
         do ip = 1, size(phases)
            call request_require(charges, phases(ip), "q", error)
            if (allocated(error)) return
         end do
         call coupling_register(coupling, "charges", charges, error)
      case default
         call fatal_error(error, "3D MOZ coupling_mode '"//trim(self%coupling_mode)//"' is not supported")
      end select
      if (allocated(error)) return
      call coupling_snapshot(coupling, ngrid=self%grid%ngrid, natom=self%grid%natom)

   end subroutine moz_3d_declare_pass

   !> Publish the grid geometry as named fields, the model's evaluation domain
   !>
   !> @param[in] self Self
   !> @param[in,out] query Field walker
   subroutine moz_3d_list_fields(self, query)
      !> Self
      class(model_moz_3d_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
      if (.not. allocated(self%grid)) return
      call query%add_int_value("ngrid", "Number of volume grid points", self%grid%ngrid)
      call query%add_int_value("natom", "Number of solute atoms", self%grid%natom)
      call query%add_real2("xyz", "Grid coordinates, bohr (3, ngrid)", self%grid%xyz)
      call query%add_real("w", "Integration weights, bohr**3 (ngrid)", self%grid%w)
      call query%add_real("xi0", "Gaussian exponents, inverse bohr (ngrid)", self%grid%xi0)
      call query%add_int("owner", "Atom owner, zero based; -1 for unowned points", &
         & self%grid%owner, zero_based=.true.)
   end subroutine moz_3d_list_fields

   !> Report pending 3D MOZ energy theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] energy Energy
   !> @param[out] error Error
   subroutine moz_3d_energy(self, coupling, energy, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
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
      call fatal_error(error, "3D MOZ energy is not implemented")
   end subroutine moz_3d_energy

   !> Report pending 3D MOZ response theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[out] error Error
   subroutine moz_3d_response(self, coupling, response, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
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
      call fatal_error(error, "3D MOZ response is not implemented")
   end subroutine moz_3d_response

   !> Report pending 3D MOZ gradient theory
   !>
   !> @param[in,out] self Self
   !> @param[in,out] coupling Coupling
   !> @param[in,out] response Response
   !> @param[in,out] gradient Gradient
   !> @param[out] error Error
   subroutine moz_3d_gradient(self, coupling, response, gradient, error)
      !> Self
      class(model_moz_3d_type), intent(inout) :: self
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
      call fatal_error(error, "3D MOZ gradient is not implemented")
   end subroutine moz_3d_gradient

   !> Release the coupling registry
   !>
   !> @param[in,out] self Model being destroyed
   subroutine destroy_moz_3d_model(self)
      !> Model being destroyed
      type(model_moz_3d_type), intent(inout) :: self
      call self%clear_couplings()
   end subroutine destroy_moz_3d_model

end module moist_model_moz_3d_type
