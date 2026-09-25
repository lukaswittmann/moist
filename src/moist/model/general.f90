!> General list-driven solvation model
module moist_model_general
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_model_type, only: solvation_model_type, solvation_model_component_type
   use moist_channels_response, only: response_type, response_clear
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & coupling_registry_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_begin_registration, coupling_set_scope, &
      & coupling_snapshot, coupling_arm, coupling_invalidate, coupling_check_mandatory, &
      & coupling_make_view, coupling_close_view
   use moist_cavity_drop, only: cavity_type_drop
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type

   implicit none(type, external)
   private

   public :: solvation_model_general, new_model_general

   !> Owning box that makes heterogeneous components storable in one array
   type :: solvation_component_slot
      !> Concrete component owned by this slot
      class(solvation_model_component_type), allocatable :: item
   end type solvation_component_slot

   !> General solvation model with one cavity and an ordered component list
   type, extends(solvation_model_type) :: solvation_model_general
      !> Authoritative cavity shared by all components
      class(cavity_type), allocatable :: cavity
      !> Ordered heterogeneous component collection
      type(solvation_component_slot), allocatable :: components(:)
      !> Whether the latest model update completed successfully
      logical :: updated = .false.
      !> Registry of model-owned host calculations
      type(coupling_registry_type), private :: couplings
      !> Component declarations are frozen after the first attempted update
      logical :: configured = .false.
      !> Force the legacy forward nuclear-gradient path
      logical :: force_forward_gradient = .false.
   contains
      final :: destroy_model
      procedure :: set_isodensity_density
      procedure :: add_component
      procedure :: update => general_update
      procedure :: get_energy => general_get_energy
      procedure :: get_response => general_get_response
      procedure :: get_gradient => general_get_gradient
      !> Build a model-owned coupling and declare its requests
      procedure :: new_coupling => general_new_coupling
      procedure :: release_coupling => general_release_coupling
      !> Stage one phase of the coupling (declare, record the grid size, arm)
      procedure, private :: stage => general_stage_coupling
      !> Stage the energy phase
      procedure :: prepare_energy => general_prepare_energy
      !> Stage the response phase
      procedure :: prepare_response => general_prepare_response
      !> Stage the gradient phase
      procedure :: prepare_gradient => general_prepare_gradient
   end type solvation_model_general

contains

   !> Install internal isodensity data and invalidate all model results
   subroutine set_isodensity_density(self, density, error)
      !> Model owning the cavity
      class(solvation_model_general), intent(inout) :: self
      !> Cartesian-monomial density matrix
      real(wp), intent(in) :: density(:, :)
      !> Invalid cavity or density
      type(error_type), allocatable, intent(out) :: error

      self%updated = .false.
      call self%couplings%invalidate()
      if (allocated(self%cavity)) then
         select type (cavity => self%cavity)
         type is (cavity_type_drop)
            if (allocated(cavity%lsf_model)) then
               select type (lsf => cavity%lsf_model)
               type is (moist_cavity_drop_lsf_isodensity_internal_type)
                  call lsf%set_density(density, error)
                  return
               end select
            end if
         end select
      end if
      call fatal_error(error, "Model requires an internal isodensity cavity")
   end subroutine set_isodensity_density

   !> Construct an empty general model around an owned copy of a cavity
   subroutine new_model_general(self, cavity, ctx, error)
      !> General model
      type(solvation_model_general), intent(out) :: self
      !> Cavity configuration to copy
      class(cavity_type), intent(in) :: cavity
      !> Shared run context, which must outlive the model
      type(moist_context_type), intent(in), target :: ctx
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Allocation status of the cavity copy
      integer :: stat

      allocate (self%cavity, source=cavity, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Failed to copy model cavity")
         return
      end if
      self%ctx => ctx
      self%cavity%ctx => ctx
      allocate (self%components(0))
      self%updated = .false.
      self%configured = .false.

   end subroutine new_model_general

   !> Append an owned copy of a component before the first update
   subroutine add_component(self, component, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Component to copy
      class(solvation_model_component_type), intent(in) :: component
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grown slot array
      type(solvation_component_slot), allocatable :: grown(:)
      !> Slot index and old slot count
      integer :: i, n

      if (self%configured) then
         call fatal_error(error, "Components cannot be added after model update")
         return
      end if
      if (.not. allocated(self%components)) allocate (self%components(0))
      n = size(self%components)
      allocate (grown(n + 1))
      do i = 1, n
         call move_alloc(self%components(i)%item, grown(i)%item)
      end do
      allocate (grown(n + 1)%item, source=component)
      grown(n + 1)%item%ctx => self%ctx
      call move_alloc(grown, self%components)

   end subroutine add_component

   !> Update the authoritative cavity and every component
   subroutine general_update(self, mol, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Molecular structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      self%updated = .false.
      self%configured = .true.
      call self%couplings%invalidate()
      if (.not. allocated(self%cavity)) then
         call fatal_error(error, "General model has no cavity")
         return
      end if
      call self%cavity%update(mol, error)
      if (allocated(error)) return
      do i = 1, size(self%components)
         call self%components(i)%item%update(mol, self%cavity, error)
         if (allocated(error)) return
      end do
      self%updated = .true.

   end subroutine general_update

   !> Accumulate the energy of every component
   !>
   !> Requires a coupling staged by `prepare_energy`; reports any other staged
   !> phase, then every missing output of the energy phase, by name
   subroutine general_get_energy(self, coupling, energy, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(inout), target :: coupling
      !> Energy accumulator
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional energy accumulator
      real(wp) :: local_energy
      !> Component index
      integer :: i
      !> Component-local borrowed read interface
      type(coupling_view_type) :: view

      if (.not. self%couplings%owns(coupling)) then
         call fatal_error(error, "Coupling belongs to a different model")
         return
      end if
      call require_updated(self, error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_energy, error)
      if (allocated(error)) return
      local_energy = 0.0_wp
      do i = 1, size(self%components)
         call coupling_make_view(coupling, i, view)
         call self%components(i)%item%get_energy(view, self%cavity, local_energy, error)
         call coupling_close_view(coupling)
         if (allocated(error)) return
      end do
      energy = energy + local_energy

   end subroutine general_get_energy

   !> Assemble the host part of the response phase
   !>
   !> `response` is cleared after input validation and returns the potential
   !> adjoint, the Gaussian amplitudes and, for a cavity with field-dependent
   !> geometry, the density weights
   !>
   !> - components accumulate into it within this call
   !> - the coupling must be staged by `prepare_response`; a stale mandatory
   !>   request of the response phase is then reported by name
   subroutine general_get_response(self, coupling, response, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(inout), target :: coupling
      !> Response list, cleared after input validation
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Shared surface accumulator
      type(cavity_surface_adjoint_type) :: acc
      !> Component index
      integer :: i
      !> Component-local borrowed read interface
      type(coupling_view_type) :: view

      if (.not. self%couplings%owns(coupling)) then
         call fatal_error(error, "Coupling belongs to a different model")
         return
      end if
      call require_updated(self, error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_response, error)
      if (allocated(error)) return
      call response_clear(response)
      call acc%init(self%cavity%ngrid)
      do i = 1, size(self%components)
         call coupling_make_view(coupling, i, view)
         call self%components(i)%item%get_response(view, self%cavity, response, error)
         call coupling_close_view(coupling)
         if (allocated(error)) return
         call coupling_make_view(coupling, i, view)
         call self%components(i)%item%get_surface_weights(view, self%cavity, acc, error)
         call coupling_close_view(coupling)
         if (allocated(error)) return
      end do
      call self%cavity%get_surface_response(acc, response, error)

   end subroutine general_get_response

   !> Accumulate the nuclear gradient of every component
   !>
   !> `response` is cleared after input validation and returns the host part
   !> of the gradient phase: the potential adjoint and the Gaussian amplitudes,
   !> which the host contracts with its own geometry derivatives
   !>
   !> - in reverse mode (the default) the components never see `get_gradient`,
   !>   so that part is collected through their `get_response` after the
   !>   surface contraction
   !> - in forward mode every component emits it from its own `get_gradient`
   !> - the two branches are exclusive, so nothing is counted twice
   !> - the coupling must be staged by `prepare_gradient`; a stale mandatory
   !>   request of the gradient phase is then reported by name
   !> - the components' own checks run against the armed phase, so the
   !>   reverse-path `get_response` calls do not demand the response-phase
   !>   requests again
   subroutine general_get_gradient(self, coupling, response, gradient, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(inout), target :: coupling
      !> Response list, cleared after input validation
      type(response_type), intent(inout) :: response
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional local gradient
      real(wp), allocatable :: local(:, :)
      !> Component index
      integer :: i
      !> Component-local borrowed read interface
      type(coupling_view_type) :: view

      if (.not. self%couplings%owns(coupling)) then
         call fatal_error(error, "Coupling belongs to a different model")
         return
      end if
      call require_updated(self, error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%cavity%nsph])) then
         call fatal_error(error, "General-model gradient shape mismatch")
         return
      end if
      call response_clear(response)
      allocate (local(3, self%cavity%nsph), source=0.0_wp)

      if (.not. self%force_forward_gradient) then
         ! Reverse mode: every component states its surface adjoint, the cavity
         ! contracts the lot once, and nothing builds a nuclear Jacobian
         block
            type(cavity_surface_adjoint_type) :: acc

            call acc%init(self%cavity%ngrid)
            do i = 1, size(self%components)
               call coupling_make_view(coupling, i, view)
               call self%components(i)%item%get_direct_gradient(view, self%cavity, &
                                                                local, error)
               call coupling_close_view(coupling)
               if (allocated(error)) return
               call coupling_make_view(coupling, i, view)
               call self%components(i)%item%get_gradient_surface_weights(view, &
                                                                         self%cavity, acc, error)
               call coupling_close_view(coupling)
               if (allocated(error)) return
            end do
            call self%cavity%get_surface_gradient(acc, local, error)
            if (allocated(error)) return
         end block
         ! Host part of the phase: the already solved charges and the amplitudes
         do i = 1, size(self%components)
            call coupling_make_view(coupling, i, view)
            call self%components(i)%item%get_response(view, self%cavity, response, error)
            call coupling_close_view(coupling)
            if (allocated(error)) return
         end do
      else
         do i = 1, size(self%components)
            call coupling_make_view(coupling, i, view)
            call self%components(i)%item%get_gradient(view, self%cavity, response, local, &
                                                      error)
            call coupling_close_view(coupling)
            if (allocated(error)) return
         end do
      end if

      gradient = gradient + local

   end subroutine general_get_gradient

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Declare cavity and component requests, then record the grid size
   !>
   !> Calculations with matching inputs are shared between components
   !>
   !> @param[in,out] self     Updated general model
   !> @param[in,out] coupling Coupling to declare
   !> @param[out]   error    Error handling
   subroutine declare_general_pass(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 0)
      call self%cavity%declare_coupling(coupling, error)
      if (allocated(error)) return
      do i = 1, size(self%components)
         call coupling_set_scope(coupling, i)
         call self%components(i)%item%declare_coupling(self%cavity, coupling, error)
         if (allocated(error)) return
      end do
      call coupling_snapshot(coupling, self%cavity%ngrid)

   end subroutine declare_general_pass

   !> Build the host coupling of an updated model
   !>
   !> - multiple couplings may coexist; release unused ones with
   !>   `release_coupling`
   subroutine general_new_coupling(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout), target :: self
      !> Coupling to build
      type(coupling_type), pointer, intent(out) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      nullify (coupling)
      call require_updated(self, error)
      if (allocated(error)) return
      call self%couplings%mint(coupling, error)
      if (allocated(error)) return
      call declare_general_pass(self, coupling, error)
      if (allocated(error)) call self%release_coupling(coupling)

   end subroutine general_new_coupling

   !> Release a coupling before destroying its parent model
   subroutine general_release_coupling(self, coupling)
      !> Owning model
      class(solvation_model_general), intent(inout), target :: self
      !> Coupling to release; other aliases must no longer be used
      type(coupling_type), pointer, intent(inout) :: coupling
      call self%couplings%release(coupling)
   end subroutine general_release_coupling

   !> Stage per-output requirements and preserve still-valid raw answers
   !>
   !> - energy staging starts a new host evaluation; response and gradient
   !>   staging reuse outputs until geometry or declared scientific inputs change
   !> - every staging starts a new host walk: `next()` begins at the first request
   !>
   !> @param[in,out] self     Updated general model
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[in]    phase    Phase index, `moist_phase_energy` and so on
   !> @param[out]   error    Foreign coupling, invalid phase or failed declaration
   subroutine general_stage_coupling(self, coupling, phase, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Phase index
      integer, intent(in) :: phase
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call require_updated(self, error)
      if (allocated(error)) return
      if (.not. self%couplings%owns(coupling)) then
         call fatal_error(error, "Coupling belongs to a different model")
         return
      end if
      if (phase == moist_phase_energy) call coupling_invalidate(coupling)
      call declare_general_pass(self, coupling, error)
      if (allocated(error)) return
      call coupling_arm(coupling, phase, error)

   end subroutine general_stage_coupling

   !> Stage the energy phase of the coupling
   subroutine general_prepare_energy(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_energy, error)

   end subroutine general_prepare_energy

   !> Stage the response phase of the coupling
   subroutine general_prepare_response(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_response, error)

   end subroutine general_prepare_response

   !> Stage the gradient phase of the coupling
   subroutine general_prepare_gradient(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(coupling, moist_phase_gradient, error)

   end subroutine general_prepare_gradient

   !> Require a successfully updated model
   subroutine require_updated(self, error)
      !> General model
      class(solvation_model_general), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%updated) call fatal_error(error, "General model must be updated first")

   end subroutine require_updated

   !> Release collections owned by the model; borrowed handles must not outlive it
   subroutine destroy_model(self)
      !> Model being destroyed
      type(solvation_model_general), intent(inout) :: self
      call self%couplings%clear()
   end subroutine destroy_model
end module moist_model_general
