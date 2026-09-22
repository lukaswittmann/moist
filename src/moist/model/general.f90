!> General list-driven solvation model
module moist_model_general
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_type, only: solvation_model_type, solvation_model_component_type, cavity_type, &
      & snapshot_cavity_coupling, resolve_coupling_phase
   use moist_channels_response, only: response_type
   use moist_channels_request, only: coupling_type, coupling_view_type, &
      & coupling_registry_type, moist_phase_none, moist_phase_energy, &
            & moist_phase_response, moist_phase_gradient
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
      !> Build a model-owned coupling and its initial grid snapshot
      procedure :: new_coupling => general_new_coupling
      procedure :: release_coupling => general_release_coupling
      !> Stage one phase of the coupling (snapshot, arm, declare scientific inputs)
      procedure :: update_coupling => general_update_coupling
      !> Stage the energy phase (alias of `update_coupling(energy=.true.)`)
      procedure :: prepare_energy => general_prepare_energy
      !> Stage the response phase (alias of `update_coupling(response=.true.)`)
      procedure :: prepare_response => general_prepare_response
      !> Stage the gradient phase (alias of `update_coupling(gradient=.true.)`)
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
   !> Requires a coupling staged by `prepare_energy` and reports any other
   !> staged phase by name; then a stale mandatory request of the energy phase,
   !> or a failure latched by a mis-shaped `set`, by name
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
      call require_staged(coupling, moist_phase_energy, "get_energy", error)
      if (allocated(error)) return
      call coupling%check_mandatory(moist_phase_energy, error)
      if (allocated(error)) return
      local_energy = 0.0_wp
      do i = 1, size(self%components)
         call coupling%make_view(i, view)
         call self%components(i)%item%get_energy(view, self%cavity, local_energy, error)
         call coupling%close_view()
         if (allocated(error)) return
      end do
      energy = energy + local_energy

   end subroutine general_get_energy

   !> Assemble the host part of the response phase
   !>
   !> `response` is cleared after input validation and returns the surface charge, the
   !> Gaussian amplitudes and, for a cavity with field-dependent geometry, the
   !> density weights; components accumulate into it within this call. The
   !> coupling must be staged by `prepare_response`; a stale mandatory request
   !> of the response phase is then reported by name
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
      call require_staged(coupling, moist_phase_response, "get_response", error)
      if (allocated(error)) return
      call coupling%check_mandatory(moist_phase_response, error)
      if (allocated(error)) return
      call response%clear()
      call acc%init(self%cavity%ngrid)
      do i = 1, size(self%components)
         call coupling%make_view(i, view)
         call self%components(i)%item%get_response(view, self%cavity, response, error)
         call coupling%close_view()
         if (allocated(error)) return
         call coupling%make_view(i, view)
         call self%components(i)%item%get_surface_weights(view, self%cavity, acc, error)
         call coupling%close_view()
         if (allocated(error)) return
      end do
      call self%cavity%get_surface_response(acc, response, error)

   end subroutine general_get_response

   !> Accumulate the nuclear gradient of every component
   !>
   !> `response` is cleared after input validation and returns the host part of the gradient
   !> phase: the surface charge and the Gaussian amplitudes, which the host
   !> contracts with its own geometry derivatives. In reverse mode (the
   !> default) the components never see `get_gradient`, so that part is
   !> collected through their `get_response` after the surface contraction;
   !> in forward mode every component emits it from its own `get_gradient`
   !> The two branches are exclusive, so nothing is counted twice. The coupling
   !> must be staged by `prepare_gradient`; a stale mandatory request of the
   !> gradient phase is then reported by name. The components' own checks run
   !> against the armed phase, so the reverse-path `get_response` calls do not
   !> demand the response-phase requests again
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
      call require_staged(coupling, moist_phase_gradient, "get_gradient", error)
      if (allocated(error)) return
      call coupling%check_mandatory(moist_phase_gradient, error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%cavity%nsph])) then
         call fatal_error(error, "General-model gradient shape mismatch")
         return
      end if
      call response%clear()
      allocate (local(3, self%cavity%nsph), source=0.0_wp)

      if (.not. self%force_forward_gradient) then
         ! Reverse mode: every component states its surface adjoint, the cavity
         ! contracts the lot once. Nothing builds a nuclear Jacobian
         block
            type(cavity_surface_adjoint_type) :: acc

            call acc%init(self%cavity%ngrid)
            do i = 1, size(self%components)
               call coupling%make_view(i, view)
               call self%components(i)%item%get_direct_gradient(view, self%cavity, &
                                                                local, error)
               call coupling%close_view()
               if (allocated(error)) return
               call coupling%make_view(i, view)
               call self%components(i)%item%get_gradient_surface_weights(view, &
                                                                         self%cavity, acc, error)
               call coupling%close_view()
               if (allocated(error)) return
            end do
            call self%cavity%get_surface_gradient(acc, local, error)
            if (allocated(error)) return
         end block
         ! Host part of the phase: the already solved charges and the amplitudes
         do i = 1, size(self%components)
            call coupling%make_view(i, view)
            call self%components(i)%item%get_response(view, self%cavity, response, error)
            call coupling%close_view()
            if (allocated(error)) return
         end do
      else
         do i = 1, size(self%components)
            call coupling%make_view(i, view)
            call self%components(i)%item%get_gradient(view, self%cavity, response, local, &
                                                      error)
            call coupling%close_view()
            if (allocated(error)) return
         end do
      end if

      gradient = gradient + local

   end subroutine general_get_gradient

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Build the host coupling of an updated model
   !>
   !> Declare cavity and component requests, sharing calculations with matching
   !> inputs, then snapshot the grid. Multiple couplings may coexist; release
   !> unused ones with `release_coupling`
   subroutine general_new_coupling(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout), target :: self
      !> Coupling to build
      type(coupling_type), pointer, intent(out) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      nullify (coupling)
      call require_updated(self, error)
      if (allocated(error)) return
      call self%couplings%mint(coupling, error)
      if (allocated(error)) return
      call coupling%begin_registration()
      call coupling%set_scope(0)
      call self%cavity%declare_coupling(coupling, error)
      if (allocated(error)) then
         call self%release_coupling(coupling)
         return
      end if
      do i = 1, size(self%components)
         call coupling%set_scope(i)
         call self%components(i)%item%declare_coupling(self%cavity, coupling, error)
         if (allocated(error)) then
            call self%release_coupling(coupling)
            return
         end if
      end do
      call snapshot_cavity_coupling(self%cavity, coupling, error)
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
   !> Energy staging starts a new host evaluation. Response and gradient staging
   !> reuse outputs until geometry or declared scientific inputs change
   subroutine general_update_coupling(self, coupling, error, energy, response, gradient)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Stage the energy phase
      logical, intent(in), optional :: energy
      !> Stage the response phase
      logical, intent(in), optional :: response
      !> Stage the gradient phase
      logical, intent(in), optional :: gradient

      !> Resolved phase index and component index
      integer :: phase, i

      call require_updated(self, error)
      if (allocated(error)) return
      call resolve_coupling_phase(energy, response, gradient, phase, error)
      if (allocated(error)) return
      if (.not. self%couplings%owns(coupling)) then
         call fatal_error(error, "Coupling belongs to a different model")
         return
      end if
      if (phase == moist_phase_energy) call coupling%invalidate()
      call coupling%begin_registration()
      call coupling%set_scope(0)
      call self%cavity%declare_coupling(coupling, error)
      if (allocated(error)) return
      do i = 1, size(self%components)
         call coupling%set_scope(i)
         call self%components(i)%item%declare_coupling(self%cavity, coupling, error)
         if (allocated(error)) return
      end do
      call snapshot_cavity_coupling(self%cavity, coupling, error)
      if (allocated(error)) return
      call coupling%arm(phase, error)

   end subroutine general_update_coupling

   !> Stage the energy phase of the coupling
   subroutine general_prepare_energy(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%update_coupling(coupling, error, energy=.true.)

   end subroutine general_prepare_energy

   !> Stage the response phase of the coupling
   subroutine general_prepare_response(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%update_coupling(coupling, error, response=.true.)

   end subroutine general_prepare_response

   !> Stage the gradient phase of the coupling
   subroutine general_prepare_gradient(self, coupling, error)
      !> Updated general model
      class(solvation_model_general), intent(inout) :: self
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%update_coupling(coupling, error, gradient=.true.)

   end subroutine general_prepare_gradient

   !> Require a successfully updated model
   subroutine require_updated(self, error)
      !> General model
      class(solvation_model_general), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%updated) call fatal_error(error, "General model must be updated first")

   end subroutine require_updated

   !> Require a coupling staged for the phase the accessor reads
   !>
   !> `check_mandatory` only asks whether the requests of a phase carry an
   !> answer, which a coupling staged for a neighbouring phase can satisfy by
   !> accident: the answers of the response phase are still fresh when the
   !> gradient phase is read. Which phase was staged is the more fundamental
   !> fact, so it is checked first and reported naming both phases
   subroutine require_staged(coupling, phase, accessor, error)
      !> Coupling handed to the accessor
      class(coupling_type), intent(in) :: coupling
      !> Phase index the accessor reads
      integer, intent(in) :: phase
      !> Name of the accessor
      character(len=*), intent(in) :: accessor
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (coupling%phase == phase) return
      if (coupling%phase == moist_phase_none) then
         call fatal_error(error, accessor//" requires a coupling staged by prepare_"// &
            & trim(phase_name(phase))//"; this coupling has never been staged")
         return
      end if
      call fatal_error(error, accessor//" requires a coupling staged by prepare_"// &
         & trim(phase_name(phase))//"; this coupling is staged for the "// &
         & trim(phase_name(coupling%phase))//" phase")

   end subroutine require_staged

   !> Name of a phase for diagnostics, e.g. "gradient"
   !>
   !> Fixed length, so that no deferred-length temporary is created for the
   !> result; callers `trim` it
   !>
   pure function phase_name(phase) result(name)
      !> Phase index
      integer, intent(in) :: phase
      !> Phase name, blank padded
      character(len=8) :: name

      select case (phase)
      case (moist_phase_energy)
         name = "energy"
      case (moist_phase_response)
         name = "response"
      case (moist_phase_gradient)
         name = "gradient"
      case default
         name = "unknown"
      end select

   end function phase_name

   !> Release collections owned by the model; borrowed handles must not outlive it
   subroutine destroy_model(self)
      !> Model being destroyed
      type(solvation_model_general), intent(inout) :: self
      call self%couplings%clear()
   end subroutine destroy_model
end module moist_model_general
