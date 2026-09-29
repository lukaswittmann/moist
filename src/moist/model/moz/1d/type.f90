!> One-dimensional MOZ model state
module moist_model_moz_1d_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_model_type, only: solvation_model_type
   use moist_channels_coupling, only: coupling_type, atomic_charge_request_type, &
      & moist_phase_energy, moist_phase_response, moist_phase_gradient, &
      & coupling_begin_registration, coupling_set_scope, coupling_register, &
      & request_require, coupling_snapshot, coupling_check_mandatory
   use moist_channels_response, only: response_type
   implicit none(type, external)
   private
   public :: model_moz_1d_type, new_moz_1d_model

   !> Direct 1D MOZ model; radial theory is pending
   type, extends(solvation_model_type) :: model_moz_1d_type
      private
      ! TODO: replace the coupling_mode string with a typed selector set at construction
      !> Electrostatic coupling source: "qat", "ec" or "multipoles" (3D only)
      character(len=16), public :: coupling_mode = "qat"
      !> Number of solute sites in the latest structure
      integer :: natom = 0
   contains
      final :: destroy_moz_1d_model
      procedure :: update => moz_1d_update
      procedure :: get_energy => moz_1d_energy
      procedure :: get_response => moz_1d_response
      procedure :: get_gradient => moz_1d_gradient
      procedure :: atom_count => moz_1d_atom_count
      !> Declare the requests of the electrostatic coupling source
      procedure :: declare_pass => moz_1d_declare_pass
   end type model_moz_1d_type

contains

   !> Attach a borrowed run context
   !>
   !> @param[out] self Model
   !> @param[in] ctx Run context, borrowed for the lifetime of the model
   !> @param[out] error Error handling
   subroutine new_moz_1d_model(self, ctx, error)
      !> Model
      type(model_moz_1d_type), intent(out) :: self
      !> Run context, borrowed for the lifetime of the model
      type(moist_context_type), intent(in), target :: ctx
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      self%ctx => ctx
   end subroutine new_moz_1d_model

   !> Record the solute site count
   !>
   !> @param[in,out] self Model
   !> @param[in] mol Solute structure
   !> @param[out] error Invalid structure
   subroutine moz_1d_update(self, mol, error)
      !> Model
      class(model_moz_1d_type), intent(inout) :: self
      !> Solute structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
      self%natom = 0
      if (.not. associated(self%ctx)) then
         call fatal_error(error, "Construct the 1D MOZ model before update")
         return
      end if
      if (mol%nat < 1) then
         call fatal_error(error, "1D MOZ requires at least one solute atom")
         return
      end if
      self%natom = mol%nat
      self%updated = .true.
   end subroutine moz_1d_update

   !> Number of sites in the latest valid structure
   !>
   !> @param[in] self Model
   function moz_1d_atom_count(self) result(natom)
      !> Model
      class(model_moz_1d_type), intent(in) :: self
      !> Site count
      integer :: natom
      natom = self%natom
   end function moz_1d_atom_count

   !> Declare the requests of the electrostatic coupling source
   !>
   !> - "qat": the solute's partial charges, `q` in every phase
   !> - "ec": pending the radial grid of the 1D model, refused by name
   !> - "multipoles": refused, the site-site potentials of 1D MOZ are spherical
   !>
   !> @param[in,out] self Updated model
   !> @param[in,out] coupling Coupling to declare
   !> @param[out] error Unknown or unsupported coupling source
   subroutine moz_1d_declare_pass(self, coupling, error)
      !> Updated model
      class(model_moz_1d_type), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-atom partial charge request
      type(atomic_charge_request_type) :: charges

      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 0)
      select case (self%coupling_mode)
      case ("qat")
         call request_require(charges, moist_phase_energy, "q", error)
         if (allocated(error)) return
         call request_require(charges, moist_phase_response, "q", error)
         if (allocated(error)) return
         call request_require(charges, moist_phase_gradient, "q", error)
         if (allocated(error)) return
         call coupling_register(coupling, "charges", charges, error)
         if (allocated(error)) return
      case ("ec")
         ! TODO: once the 1D model owns its radial grid, register `radial_potential`
         ! (phi in energy, response, gradient) and `atomic_charges` (q in all three
         ! phases, the tail charges of the Ng split), publish ngrid/natom/r through
         ! a `list_fields` override and snapshot ngrid from the radial grid and natom
         call fatal_error(error, "1D MOZ EC coupling requires the radial grid (pending grid refactor)")
         return
      case ("multipoles")
         call fatal_error(error, "1D MOZ coupling_mode 'multipoles' is not supported: "// &
            & "1D site-site potentials are spherical, multipoles are 3D only")
         return
      case default
         call fatal_error(error, "1D MOZ coupling_mode '"//trim(self%coupling_mode)//"' is not supported")
         return
      end select
      call coupling_snapshot(coupling, natom=self%natom)

   end subroutine moz_1d_declare_pass

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
