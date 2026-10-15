!> Domain-independent solvation model interface and shared coupling lifecycle
!>
!> Every family (continuum, MOZ 1D, MOZ 3D, ...) shares one host coupling
!> protocol: mint a coupling, stage a phase, walk it and answer by name
!>
!> - owns that protocol as concrete base procedures: `prepare_build(coupling)`
!>   stages the build phase, `build(coupling)` reads its answers
!> - a family supplies `update`, `get_energy`, `get_response`, `get_gradient`,
!>   `atom_count` and the one coupling hook, `declare_pass`, that declares its
!>   own requests and snapshots the coupling extents; one whose update waits
!>   for host answers overrides `defers_build` and `build`, calling
!>   `model_build` first
!> - `updated` and `built` record the latest successful update and build;
!>   `new_coupling` requires an updated model, the staging calls a built one
!> - a family with an evaluation domain overrides `list_fields` to publish
!>   its named arrays; the base publishes none
module moist_model_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_channels_response, only: response_type
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_registry_type, &
      & moist_phase_energy, moist_phase_response, moist_phase_gradient, moist_phase_build, &
      & coupling_arm, coupling_invalidate, coupling_check_mandatory

   implicit none(type, external)
   private

   !> Family overrides can call the base invalidation with their polymorphic self
   public :: solvation_model_type, model_build, invalidate_model, copy_model_state

   !> Abstract base solvation model
   type, abstract :: solvation_model_type
      !> Borrowed run context (verbosity/debug/timer)
      !>
      !> Set at construction; owned by the caller, never allocated or freed here
      type(moist_context_type), pointer :: ctx => null()
      !> Whether the latest model update completed successfully
      logical :: updated = .false.
      !> Whether the model is complete for evaluation: built by the latest
      !> update, or by the latest `build` when the update waits for host answers
      logical :: built = .false.
      !> Registry of model-owned host couplings
      type(coupling_registry_type), private :: couplings
   contains

      !> Checked copy into the same concrete family; unsupported unless overridden
      procedure :: copy_from => model_copy_from_default
      !> Dispatch copying through the model interface, including polymorphic variables
      generic :: assignment(=) => assign
      !> GCC 13 misdispatches inherited defined assignment when this binding is non_overridable
      procedure, private :: assign => model_assign

      !* ------------------------------ Update and build ------------------------------ *!

      procedure(update_model), deferred :: update
      !> Whether the latest update left part of the model to `build`; never by default
      procedure :: defers_build => model_defers_build_default
      !> Stage the build phase: clear and declare the coupling
      procedure, non_overridable :: prepare_build => model_prepare_build
      !> Read the build-phase answers of a staged coupling and complete the model
      procedure :: build => model_build

      !* --------------------------------- Evaluation --------------------------------- *!

      procedure(get_model_energy), deferred :: get_energy
      procedure(get_model_response), deferred :: get_response
      procedure(get_model_gradient), deferred :: get_gradient

      !* ----------------------------- Evaluation domain ------------------------------ *!

      !> Number of atoms in the family's owned geometry
      procedure(model_atom_count_i), deferred :: atom_count
      !> Declare this family's coupling requests and record the extents
      procedure(declare_model_pass), deferred :: declare_pass
      !> Declare the named fields of the evaluation domain; none by default
      procedure :: list_fields => model_list_fields

      !* ------------------------------- Host couplings ------------------------------- *!

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

      !* ------------------------------ State and guards ------------------------------ *!

      !> Whether the latest update completed
      procedure :: is_updated => model_is_updated
      !> Whether the model is built for evaluation
      procedure :: is_built => model_is_built
      !> Invalidate cached results and every model-owned coupling
      procedure :: invalidate => invalidate_model
      !> Release the coupling registry; called from each family's finalizer
      procedure :: clear_couplings => model_clear_couplings
      !> Require a successfully updated model
      procedure :: require_updated => model_require_updated
      !> Require a built model
      procedure :: require_built => model_require_built
      !> Require a coupling minted by this model
      procedure :: require_owned => model_require_owned

   end type solvation_model_type

   !> Deferred procedure contracts
   abstract interface

      !> Update the solvation model with the current molecular structure
      !>
      !> - calculates all structure-dependent properties that need no host state
      !> - invalidates first and sets `updated` and `built` on success
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
         type(structure_type), intent(in) :: mol
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
      !> - before its domain exists only the requests of the build phase
      !> - ends with `coupling_snapshot`, giving the counts the family's own
      !>   evaluation domain has (the continuum's from its cavity)
      !>
      !> @param[in,out] self Solvation model
      !> @param[in,out] coupling Coupling to declare
      !> @param[out] error Error handling
      subroutine declare_model_pass(self, coupling, error)
         import solvation_model_type, coupling_type, error_type
         implicit none(type, external)
         !> Solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Coupling to declare
         type(coupling_type), intent(inout) :: coupling
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine declare_model_pass

   end interface

contains

   !> Assign through the family's checked copy hook
   !>
   !> @param[inout] self Allocated destination of the same concrete type as source
   !> @param[in] source Model whose run context remains borrowed
   impure elemental subroutine model_assign(self, source)
      !> Destination model
      class(solvation_model_type), intent(inout) :: self
      !> Source model
      class(solvation_model_type), intent(in) :: source
      !> Assignment has no error argument; checked callers use copy_from directly
      type(error_type), allocatable :: error

      call self%copy_from(source, error)
      if (allocated(error)) error stop error%message
   end subroutine model_assign

   !> Report copying unsupported by a custom model family
   !>
   !> @param[inout] self Destination, unchanged on error
   !> @param[in] source Source, unchanged
   !> @param[out] error Unsupported copy operation
   subroutine model_copy_from_default(self, source, error)
      !> Destination model
      class(solvation_model_type), intent(inout) :: self
      !> Source model
      class(solvation_model_type), intent(in) :: source
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "This solvation model does not implement copying")
   end subroutine model_copy_from_default

   !> Prepare a same-family copy with a fresh coupling registry and borrowed context
   !>
   !> @param[inout] self Destination; invalidated before the family replaces its owned data
   !> @param[in] source Source model
   !> @param[out] error Incompatible concrete types; destination unchanged
   !> @return Whether family data need copying; false for self-copy or an error
   logical function copy_model_state(self, source, error) result(copy_data)
      !> Destination and source models
      class(solvation_model_type), intent(inout), target :: self
      class(solvation_model_type), intent(in), target :: source
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Borrowed pointer used only to compare model identity
      class(solvation_model_type), pointer :: destination

      copy_data = .false.
      destination => self
      if (associated(destination, source)) return
      if (.not. same_type_as(self, source)) then
         call fatal_error(error, "Model copies require the same concrete type")
         return
      end if
      call self%clear_couplings()
      call self%invalidate()
      self%ctx => source%ctx
      self%updated = source%updated
      self%built = source%built
      copy_data = .true.
   end function copy_model_state

   !* ================================================================================= *!
   !*                                 Update and build                                  *!
   !* ================================================================================= *!

   !> Default: `update` builds the whole model
   !>
   !> @param[in] self Model
   function model_defers_build_default(self) result(defers)
      !> Model
      class(solvation_model_type), intent(in) :: self
      !> Whether `update` leaves part of the model to `build`
      logical :: defers

      defers = .false.

   end function model_defers_build_default

   !> Stage the build phase of a coupling
   !>
   !> - clears every answer of the coupling; the atom count comes from the
   !>   structure of the latest `update`
   !> - a family whose domain waits for host answers is not built until the
   !>   next `build` succeeds, and every coupling's answers over it go; one
   !>   that `update` built completely keeps everything
   !>
   !> @param[in,out] self     Updated model
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[out]    error    Model not updated, foreign coupling or failed declaration
   subroutine model_prepare_build(self, coupling, error)
      !> Updated model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_updated(error)
      if (allocated(error)) return
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call coupling_invalidate(coupling)
      if (self%defers_build()) then
         self%built = .false.
         call self%couplings%invalidate()
      end if
      call self%declare_pass(coupling, error)
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_build, error)

   end subroutine model_prepare_build

   !> Read the build-phase answers of a staged coupling
   !>
   !> - the coupling must be staged by `prepare_build` and fully answered; a
   !>   missing output is reported by name
   !> - its build-phase answers stay for the later phases; every other
   !>   coupling of the model is invalidated
   !> - `built` is false from entry; on a model that `update` built completely
   !>   the walk is empty, nothing is rebuilt, the model is built again and
   !>   the coupling re-declared
   !> - a family that `defers_build` overrides `build`, calls this first and
   !>   then completes itself: built, the coupling re-declared and so
   !>   unstaged, only at the end of the override
   !>
   !> @param[in,out] self     Updated model
   !> @param[in,out] coupling Coupling staged by `prepare_build`
   !> @param[out]    error    Error handling
   subroutine model_build(self, coupling, error)
      !> Updated model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling staged by `prepare_build`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      self%built = .false.
      call coupling_check_mandatory(coupling, moist_phase_build, error)
      if (allocated(error)) return
      call self%couplings%invalidate(except=coupling)
      if (self%defers_build()) return
      self%built = .true.
      call self%declare_pass(coupling, error)

   end subroutine model_build

   !* ================================================================================= *!
   !*                                 Evaluation domain                                 *!
   !* ================================================================================= *!

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

   !* ================================================================================= *!
   !*                                  Host couplings                                   *!
   !* ================================================================================= *!

   !> Build the host coupling of an updated model
   !>
   !> - before a deferred domain is built the snapshot has no point count,
   !>   so grid outputs are refused until then
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
   !> - energy staging starts a new host evaluation and drops every answer but
   !>   those of the build phase; response and gradient staging reuse outputs
   !>   until geometry or declared scientific inputs change
   !> - every staging starts a new host walk: `next()` begins at the first request
   !>
   !> @param[in,out] self     Built model
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[in]    phase    Phase index, `moist_phase_energy` and so on
   !> @param[out]   error    Model not built, foreign coupling, invalid phase or failed declaration
   subroutine model_stage_coupling(self, coupling, phase, error)
      !> Built model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Phase index
      integer, intent(in) :: phase
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_built(error)
      if (allocated(error)) return
      call self%require_owned(coupling, error)
      if (allocated(error)) return
      if (phase == moist_phase_energy) call coupling_invalidate(coupling, keep_build=.true.)
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
      !> Built model
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
      !> Built model
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
      !> Built model
      class(solvation_model_type), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_gradient, error)

   end subroutine model_prepare_gradient

   !* ================================================================================= *!
   !*                                 State and guards                                  *!
   !* ================================================================================= *!

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

   !> Whether the model is built for evaluation
   !>
   !> @param[in] self Model
   function model_is_built(self) result(built)
      !> Model
      class(solvation_model_type), intent(in) :: self
      !> Build status
      logical :: built
      built = self%built
   end function model_is_built

   !> Invalidate cached results and host answers
   !>
   !> @param[in,out] self Model
   subroutine invalidate_model(self)
      !> Model
      class(solvation_model_type), intent(inout) :: self
      self%updated = .false.
      self%built = .false.
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

   !> Require a built model
   !>
   !> @param[in] self Model
   !> @param[out] error Error handling
   subroutine model_require_built(self, error)
      !> Model
      class(solvation_model_type), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%built) call fatal_error(error, "Solvation model must be built first")
   end subroutine model_require_built

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

end module moist_model_type
