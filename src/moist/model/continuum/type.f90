!> Typed continuum model with direct accumulation and family-owned geometry
module moist_model_continuum_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_type, only: model_continuum_component_type
   use moist_cavity_drop, only: cavity_type_drop
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_model_type, only: solvation_model_type
   use moist_channels_response, only: response_type, response_clear
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_begin_registration, coupling_set_scope, &
      & coupling_snapshot, coupling_check_mandatory, &
      & coupling_make_view, coupling_close_view

   implicit none(type, external)
   private

   public :: model_continuum_type, new_continuum_model

   !> Owning box that makes heterogeneous components storable in one array
   type :: solvation_component_slot
      !> Concrete component owned by this slot
      class(model_continuum_component_type), allocatable :: item
   end type solvation_component_slot

   !> Typed continuum model with one cavity and an ordered component list
   type, extends(solvation_model_type) :: model_continuum_type
      !> Owned geometry shared by this family's components
      class(cavity_type), allocatable :: cavity
      !> Ordered heterogeneous component collection
      type(solvation_component_slot), allocatable :: components(:)
      !> Component declarations are frozen after the first attempted update
      logical :: configured = .false.
      !> Force the legacy forward nuclear-gradient path
      logical :: force_forward_gradient = .false.
   contains
      final :: destroy_model
      procedure :: atom_count => model_atom_count
      procedure :: use_forward_gradient => model_use_forward_gradient
      procedure :: set_isodensity_density => set_cavity_density
      procedure :: add_component
      procedure :: update => continuum_update
      procedure :: get_energy => continuum_get_energy
      procedure :: get_response => continuum_get_response
      procedure :: get_gradient => continuum_get_gradient
      !> Declare the cavity and component requests of one coupling pass
      procedure :: declare_pass => declare_continuum_pass
   end type model_continuum_type

contains

   !> Construct an empty continuum model around an owned copy of a cavity
   !>
   !> The geometry copy must not carry state bound to its own address
   !>
   !> @param[out] self Instance
   !> @param[in] cavity Live cavity
   !> @param[in] ctx Borrowed run context
   !> @param[out] error Error handling
   subroutine new_continuum_model(self, cavity, ctx, error)
      !> Continuum model
      type(model_continuum_type), intent(out) :: self
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

   end subroutine new_continuum_model

   !> Append an owned copy of a component before the first update
   !>
   !> @param[in,out] self Instance
   !> @param[in] component Compatible component to copy
   !> @param[out] error Error handling
   subroutine add_component(self, component, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Component to copy
      class(model_continuum_component_type), intent(in) :: component
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Grown slot array
      type(solvation_component_slot), allocatable :: grown(:)
      !> Slot index and old slot count
      integer :: i, n

      if (.not. allocated(self%cavity)) then
         call fatal_error(error, "Construct the model before adding components")
         return
      end if
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

   !> Update the cavity and every component
   !>
   !> @param[in,out] self Instance
   !> @param[in] mol Molecular structure
   !> @param[out] error Error handling
   subroutine continuum_update(self, mol, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Molecular structure
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      call self%invalidate()
      self%configured = .true.
      if (.not. allocated(self%cavity)) then
         call fatal_error(error, "Continuum model has no cavity")
         return
      end if
      call self%cavity%update(mol, error)
      if (allocated(error)) return
      do i = 1, size(self%components)
         call self%components(i)%item%update(mol, self%cavity, error)
         if (allocated(error)) return
      end do
      self%updated = .true.

   end subroutine continuum_update

   !> Accumulate the energy of every component
   !>
   !> Requires a coupling staged by `prepare_energy`; reports any other staged
   !> phase, then every missing output of the energy phase, by name
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   !> @param[in,out] energy Energy accumulator
   !> @param[out] error Error handling
   subroutine continuum_get_energy(self, coupling, energy, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
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

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
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

   end subroutine continuum_get_energy

   !> Assemble the host part of the response phase
   !>
   !> `response` is cleared after input validation and returns the potential
   !> adjoint, the Gaussian amplitudes and, for a cavity with field-dependent
   !> geometry, the density weights
   !>
   !> - components accumulate into it within this call
   !> - the coupling must be staged by `prepare_response`; a stale mandatory
   !>   request of the response phase is then reported by name
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   !> @param[in,out] response Host response accumulator
   !> @param[out] error Error handling
   subroutine continuum_get_response(self, coupling, response, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(inout), target :: coupling
      !> Response list, cleared after input validation
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Shared cavity adjoint accumulator
      type(cavity_surface_adjoint_type) :: acc
      !> Component index
      integer :: i
      !> Component-local borrowed read interface
      type(coupling_view_type) :: view

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
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

   end subroutine continuum_get_response

   !> Accumulate the nuclear gradient of every component
   !>
   !> `response` is cleared after input validation and returns the host part
   !> of the gradient phase: the potential adjoint and the Gaussian amplitudes,
   !> which the host contracts with its own geometry derivatives
   !>
   !> - in reverse mode (the default) the components never see `get_gradient`,
   !>   so that part is collected through their `get_response` after the
   !>   cavity contraction
   !> - in forward mode every component emits it from its own `get_gradient`
   !> - the two branches are exclusive, so nothing is counted twice
   !> - the coupling must be staged by `prepare_gradient`; a stale mandatory
   !>   request of the gradient phase is then reported by name
   !> - the components' own checks run against the armed phase, so the
   !>   reverse-path `get_response` calls do not demand the response-phase
   !>   requests again
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Host coupling
   !> @param[in,out] response Host response accumulator
   !> @param[in,out] gradient Nuclear gradient accumulator
   !> @param[out] error Error handling
   subroutine continuum_get_gradient(self, coupling, response, gradient, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
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

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call self%require_updated(error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_gradient, error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%cavity%nsph])) then
         call fatal_error(error, "Continuum-model gradient shape mismatch")
         return
      end if
      call response_clear(response)
      allocate (local(3, self%cavity%nsph), source=0.0_wp)

      if (.not. self%force_forward_gradient) then
         ! Reverse mode: every component states its geometry adjoint, the
         ! cavity contracts the lot once, and nothing builds a nuclear Jacobian
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

   end subroutine continuum_get_gradient

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Declare cavity and component requests, then record the extents
   !>
   !> Calculations with matching inputs are shared between components
   !>
   !> @param[in,out] self     Updated continuum model
   !> @param[in,out] coupling Coupling to declare
   !> @param[out]   error    Error handling
   subroutine declare_continuum_pass(self, coupling, error)
      !> Updated continuum model
      class(model_continuum_type), intent(inout) :: self
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
      call coupling_snapshot(coupling, self%cavity%ngrid, self%cavity%nsph)

   end subroutine declare_continuum_pass

   !> Release collections owned by the model; borrowed handles must not outlive it
   !>
   !> @param[in,out] self Instance
   subroutine destroy_model(self)
      !> Model being destroyed
      type(model_continuum_type), intent(inout) :: self
      call self%clear_couplings()
   end subroutine destroy_model

   !> Number of atoms in the owned geometry
   !>
   !> @param[in] self Model
   function model_atom_count(self) result(nat)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> Atom count
      integer :: nat
      nat = 0
      if (allocated(self%cavity)) nat = self%cavity%nsph
   end function model_atom_count

   !> Select the forward gradient reference path
   !>
   !> @param[in,out] self Model
   !> @param[in] enabled Whether to use forward geometry derivatives
   subroutine model_use_forward_gradient(self, enabled)
      !> Model
      class(model_continuum_type), intent(inout) :: self
      !> Whether to use forward geometry derivatives
      logical, intent(in) :: enabled
      self%force_forward_gradient = enabled
   end subroutine model_use_forward_gradient

   !> Install density data and invalidate the model
   !>
   !> @param[in,out] self Model owning the cavity
   !> @param[in] density Cartesian-monomial density matrix
   !> @param[out] error Incompatible cavity or density
   subroutine set_cavity_density(self, density, error)
      !> Model owning the cavity
      class(model_continuum_type), intent(inout), target :: self
      !> Cartesian-monomial density matrix
      real(wp), intent(in) :: density(:, :)
      !> Incompatible cavity or density
      type(error_type), allocatable, intent(out) :: error
      call self%invalidate()
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
   end subroutine set_cavity_density

end module moist_model_continuum_type
