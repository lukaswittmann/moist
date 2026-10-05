!> Domain-independent solvation model interface and shared coupling lifecycle
!>
!> Every family (continuum, MOZ 1D, MOZ 3D, ...) shares one host coupling
!> protocol: mint a coupling, stage a phase, walk it and answer by name
!>
!> - owns that protocol as concrete base procedures
!> - a family supplies only `update`, `get_energy`, `get_response`,
!>   `get_gradient`, `atom_count` and the one coupling hook, `declare_pass`,
!>   that declares its own requests and snapshots the coupling extents
!> - a family with an evaluation domain overrides `list_fields` to publish
!>   its named arrays; the base publishes none
module moist_model_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_channels_response, only: response_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_registry_type, &
      & moist_phase_energy, moist_phase_response, moist_phase_gradient, &
      & coupling_arm, coupling_invalidate

   implicit none(type, external)
   private

   public :: solvation_model_type

   !> Abstract base solvation model
   type, abstract :: solvation_model_type
      !> Borrowed run context (verbosity/debug/timer)
      !>
      !> Set at construction; owned by the caller, never allocated or freed here
      type(moist_context_type), pointer :: ctx => null()
      !> Whether the latest model update completed successfully
      logical :: updated = .false.
      !> Registry of model-owned host couplings
      type(coupling_registry_type), private :: couplings
   contains

      procedure(update_model), deferred :: update
      procedure(get_model_energy), deferred :: get_energy
      procedure(get_model_response), deferred :: get_response
      procedure(get_model_gradient), deferred :: get_gradient
      !> Number of atoms in the family's owned geometry
      procedure(model_atom_count_i), deferred :: atom_count
      !> Declare this family's coupling requests and record the extents
      procedure(declare_model_pass), deferred :: declare_pass
      !> Declare the named fields of the evaluation domain; none by default
      procedure :: list_fields => model_list_fields

      !> Whether the latest update completed
      procedure :: is_updated => model_is_updated
      !> Invalidate cached results and every model-owned coupling
      procedure :: invalidate => invalidate_model
      !> Release the coupling registry; called from each family's finalizer
      procedure :: clear_couplings => model_clear_couplings
      !> Require a successfully updated model
      procedure :: require_updated => model_require_updated
      !> Require a coupling minted by this model
      procedure :: require_owned => model_require_owned
      !> Build a model-owned coupling and declare its requests
      procedure :: new_coupling => model_new_coupling
      procedure :: release_coupling => model_release_coupling
      !> Stage one phase of the coupling (declare, snapshot the grid size, arm)
      procedure, private :: stage => model_stage_coupling
      !> Stage the energy phase
      procedure :: prepare_energy => model_prepare_energy
      !> Stage the response phase
      procedure :: prepare_response => model_prepare_response
      !> Stage the gradient phase
      procedure :: prepare_gradient => model_prepare_gradient

   end type solvation_model_type

   !> Deferred procedure contracts
   abstract interface

      !> Update the solvation model with the current molecular structure
      !>
      !> - calculates all structure-dependent properties
      !>
      !> @param[in,out] self Instance of the solvation model
      !> @param[in] mol Molecular structure data
      !> @param[out] error Error handling
      subroutine update_model(self, mol, error)
         import solvation_model_type, structure_type, error_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_model

      !> Evaluate the solvation energy
      !>
      !> @param[in,out] self Instance of the solvation model
      !> @param[in,out] coupling Wavefunction data
      !> @param[in,out] energy Solvation energy
      !> @param[out] error Error handling
      subroutine get_model_energy(self, coupling, energy, error)
         import solvation_model_type, structure_type, wp, error_type, coupling_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(inout), target :: coupling
         !> Solvation energy
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_energy

      !> Get the solvation response (only for self-consistent models)
      !>
      !> @param[in,out] self Instance of the solvation model
      !> @param[in,out] coupling Wavefunction data
      !> @param[in,out] response Solvation response for the component
      !> @param[out] error Error handling
      subroutine get_model_response(self, coupling, response, error)
         import solvation_model_type, structure_type, wp, error_type, response_type, coupling_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(inout), target :: coupling
         !> Solvation response for the component
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_response

      !> Get the solvation energy gradient and the host part of the phase
      !>
      !> `response` is cleared once the request is accepted and returns what the
      !> host contracts with its own geometry derivatives (potential adjoint,
      !> Gaussian amplitudes); a rejected request, including an unimplemented
      !> theory, leaves it untouched
      !>
      !> @param[in,out] self Instance of the solvation model
      !> @param[in,out] coupling Wavefunction data
      !> @param[in,out] response Host part of the gradient phase
      !> @param[in,out] gradient Solvation gradient
      !> @param[out] error Error handling
      subroutine get_model_gradient(self, coupling, response, gradient, error)
         import solvation_model_type, structure_type, wp, error_type, response_type, coupling_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(inout), target :: coupling
         !> Host part of the gradient phase
         type(response_type), intent(inout) :: response
         !> Solvation gradient
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_gradient

      !> Number of atoms in the family's owned geometry, zero before construction
      !>
      !> @param[in] self Instance of the solvation model
      function model_atom_count_i(self) result(nat)
         import solvation_model_type
         implicit none(type, external)
         !> Instance of the solvation model
         class(solvation_model_type), intent(in) :: self
         !> Atom count
         integer :: nat
      end function model_atom_count_i

      !> Declare this family's coupling requests, then record the extents
      !>
      !> Called on an updated model; ends with `coupling_snapshot`, giving the
      !> counts the family's own evaluation domain has (`ngrid`, `natom`)
      !>
      !> @param[in,out] self Updated solvation model
      !> @param[in,out] coupling Coupling to declare
      !> @param[out] error Error handling
      subroutine declare_model_pass(self, coupling, error)
         import solvation_model_type, coupling_type, error_type
         implicit none(type, external)
         !> Updated solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Coupling to declare
         type(coupling_type), intent(inout) :: coupling
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine declare_model_pass

   end interface

contains

   !> Whether the latest update completed
   !>
   !> @param[in] self Model
   function model_is_updated(self) result(updated)
      !> Model
      class(solvation_model_type), intent(in) :: self
      !> Update status
      logical :: updated
      updated = self%updated
   end function model_is_updated

   !> Declare the named fields of the model's evaluation domain
   !>
   !> Default: none; a family with a cavity or a grid overrides it
   !>
   !> @param[in] self Model
   !> @param[in,out] query Field walker
   subroutine model_list_fields(self, query)
      !> Model
      class(solvation_model_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
   end subroutine model_list_fields

   !> Invalidate cached results and host answers
   !>
   !> @param[in,out] self Model
   subroutine invalidate_model(self)
      !> Model
      class(solvation_model_type), intent(inout) :: self
      self%updated = .false.
      call self%couplings%invalidate()
   end subroutine invalidate_model

   !> Release every coupling minted by this model
   !>
   !> Called from each family's own `final` procedure: a `final` subroutine
   !> takes a non-polymorphic dummy, which an abstract type cannot provide, so
   !> the registry stays private here and each concrete family destroys
   !> through this public hook instead of its own `final` touching it directly
   !>
   !> @param[in,out] self Model being destroyed
   subroutine model_clear_couplings(self)
      !> Model being destroyed
      class(solvation_model_type), intent(inout) :: self
      call self%couplings%clear()
   end subroutine model_clear_couplings

   !> Require a successfully updated model
   !>
   !> @param[in] self Model
   !> @param[out] error Error handling
   subroutine model_require_updated(self, error)
      !> Model
      class(solvation_model_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%updated) call fatal_error(error, "Solvation model must be updated first")
   end subroutine model_require_updated

   !> Require a coupling minted by this model
   !>
   !> @param[in] self Model
   !> @param[in] coupling Coupling to check
   !> @param[out] error Error handling
   subroutine model_require_owned(self, coupling, error)
      !> Model
      class(solvation_model_type), intent(in), target :: self
      !> Coupling to check
      class(coupling_type), intent(in), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%couplings%owns(coupling)) then
         call fatal_error(error, "Coupling belongs to a different model")
      end if
   end subroutine model_require_owned

   !> Build the host coupling of an updated model
   !>
   !> - multiple couplings may coexist; release unused ones with
   !>   `release_coupling`
   !>
   !> @param[in,out] self Instance
   !> @param[out] coupling Host coupling
   !> @param[out] error Error handling
   subroutine model_new_coupling(self, coupling, error)
      !> Updated model
      class(solvation_model_type), intent(inout), target :: self
      !> Coupling to build
      type(coupling_type), pointer, intent(out) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      nullify (coupling)
      call self%require_updated(error)
      if (allocated(error)) return
      call self%couplings%mint(coupling, error)
      if (allocated(error)) return
      call self%declare_pass(coupling, error)
      if (allocated(error)) call self%release_coupling(coupling)

   end subroutine model_new_coupling

   !> Release a coupling before destroying its parent model
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   subroutine model_release_coupling(self, coupling)
      !> Owning model
      class(solvation_model_type), intent(inout), target :: self
      !> Coupling to release; other aliases must no longer be used
      type(coupling_type), pointer, intent(inout) :: coupling
      call self%couplings%release(coupling)
   end subroutine model_release_coupling

   !> Stage per-output requirements and preserve still-valid raw answers
   !>
   !> - energy staging starts a new host evaluation; response and gradient
   !>   staging reuse outputs until geometry or declared scientific inputs change
   !> - every staging starts a new host walk: `next()` begins at the first request
   !>
   !> @param[in,out] self     Updated model
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[in]    phase    Phase index, `moist_phase_energy` and so on
   !> @param[out]   error    Foreign coupling, invalid phase or failed declaration
   subroutine model_stage_coupling(self, coupling, phase, error)
      !> Updated model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Phase index
      integer, intent(in) :: phase
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_updated(error)
      if (allocated(error)) return
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      if (phase == moist_phase_energy) call coupling_invalidate(coupling)
      call self%declare_pass(coupling, error)
      if (allocated(error)) return
      call coupling_arm(coupling, phase, error)

   end subroutine model_stage_coupling

   !> Stage the energy phase of the coupling
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   !> @param[out] error Error handling
   subroutine model_prepare_energy(self, coupling, error)
      !> Updated model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_energy, error)

   end subroutine model_prepare_energy

   !> Stage the response phase of the coupling
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   !> @param[out] error Error handling
   subroutine model_prepare_response(self, coupling, error)
      !> Updated model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_response, error)

   end subroutine model_prepare_response

   !> Stage the gradient phase of the coupling
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   !> @param[out] error Error handling
   subroutine model_prepare_gradient(self, coupling, error)
      !> Updated model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_gradient, error)

   end subroutine model_prepare_gradient

end module moist_model_type
