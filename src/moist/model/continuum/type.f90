!> Typed continuum model with direct accumulation and family-owned geometry
module moist_model_continuum_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_utils_timer, only: cat_build, cat_hessian
   use moist_cavity_type, only: cavity_type, surface_adjoint_response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type, &
      & cavity_hessian_exchange_type
   use moist_model_continuum_component_type, only: model_continuum_component_type
   use moist_model_type, only: solvation_model_type, model_build, invalidate_model, &
      & model_release_coupling
   use moist_channels_response, only: response_type, response_clear
   use moist_channels_fields, only: field_query_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, moist_phase_hvp_geometry, moist_phase_hvp, &
      & coupling_begin_registration, coupling_set_scope, &
      & coupling_snapshot, coupling_arm, coupling_check_mandatory, &
      & coupling_set_hvp_dirs, coupling_set_hvp_surface_tangent, coupling_get_hvp_dirs, &
      & coupling_make_view, coupling_close_view

   implicit none(type, external)
   private

   public :: model_continuum_type, new_continuum_model

   !> Owning box that makes heterogeneous components storable in one array
   type :: solvation_component_slot
      !> Concrete component owned by this slot
      class(model_continuum_component_type), allocatable :: item
      !> Energy of the latest evaluation; unallocated until computed
      real(wp), allocatable :: energy
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
      !> Cached Stage-1 surface/LSF tangent exchange for the active HVP batch
      type(cavity_hessian_exchange_type), private :: hvp_exchange
      !> Borrowed pointer to the coupling owning the active HVP batch
      type(coupling_type), pointer, private :: hvp_coupling => null()
   contains
      final :: destroy_model
      procedure :: atom_count => model_atom_count
      procedure :: use_forward_gradient => model_use_forward_gradient
      procedure :: add_component
      !> Number of components, available before the first update
      procedure :: component_count => continuum_component_count
      !> Name of one component, available before the first update
      procedure :: component_name => continuum_component_name
      !> One-line description of one component, available before the first update
      procedure :: component_description => continuum_component_description
      !> Publish the results of one component, e.g. its latest energy
      procedure :: list_component_fields => continuum_list_component_fields
      procedure :: update => continuum_update
      !> Whether the cavity waits for the build phase
      procedure :: defers_build => continuum_defers_build
      !> Read the build phase, then complete a cavity that waited for it and the components
      procedure :: build => continuum_build
      procedure :: invalidate => continuum_invalidate
      procedure :: release_coupling => continuum_release_coupling
      procedure :: get_energy => continuum_get_energy
      procedure :: get_response => continuum_get_response
      procedure :: get_gradient => continuum_get_gradient
      procedure :: prepare_hvp_geometry => continuum_prepare_hvp_geometry
      procedure :: prepare_hvp => continuum_prepare_hvp
      procedure :: get_hvp => continuum_get_hvp
      procedure :: get_hessian => continuum_get_hessian
      !> Declare the cavity and component requests of one coupling pass
      procedure :: declare_pass => declare_continuum_pass
      !> Publish the fields of the owned cavity
      procedure :: list_fields => continuum_list_fields
      !> Print the settings of the cavity and of every component
      procedure :: print_parameters => continuum_print_parameters
      !> Print the radii of the current structure
      procedure :: print_radii => continuum_print_radii
   end type model_continuum_type

   !> Adjoint response adapter binding the model's components and coupling into cavity Hessian traversal
   type, extends(surface_adjoint_response_type) :: continuum_adjoint_response_type
      type(solvation_component_slot), pointer :: components(:) => null()
      type(coupling_type), pointer :: coupling => null()
   contains
      procedure :: apply => continuum_adjoint_response_apply
   end type continuum_adjoint_response_type

contains

   !> Construct an empty continuum model around an owned copy of a cavity
   !>
   !> The geometry copy must not carry state bound to its own address
   !>
   !> - a cavity with its own context keeps it; one without runs on `ctx`
   !>
   !> @param[out] self Instance
   !> @param[in] ctx Borrowed run context
   !> @param[in] cavity Live cavity
   !> @param[out] error Error handling
   subroutine new_continuum_model(self, ctx, cavity, error)
      !> Continuum model
      type(model_continuum_type), intent(out) :: self
      !> Shared run context, which must outlive the model
      type(moist_context_type), intent(in), target :: ctx
      !> Cavity configuration to copy
      class(cavity_type), intent(in) :: cavity
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
      if (.not. associated(self%cavity%ctx)) self%cavity%ctx => ctx
      allocate (self%components(0))
      self%updated = .false.
      self%built = .false.
      self%configured = .false.

   end subroutine new_continuum_model

   !> Append an owned copy of a component before the first update
   !>
   !> - a component with its own context keeps it; one without runs on the model's
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
         call move_alloc(self%components(i)%energy, grown(i)%energy)
      end do
      allocate (grown(n + 1)%item, source=component)
      if (.not. associated(grown(n + 1)%item%ctx)) grown(n + 1)%item%ctx => self%ctx
      call move_alloc(grown, self%components)

   end subroutine add_component

   !> Drop the stored energy of every component
   !>
   !> @param[in,out] self Instance
   subroutine clear_component_energies(self)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self

      !> Component index
      integer :: i

      if (.not. allocated(self%components)) return
      do i = 1, size(self%components)
         if (allocated(self%components(i)%energy)) deallocate (self%components(i)%energy)
      end do

   end subroutine clear_component_energies

   !> Update the cavity and, unless it waits for the build phase, every component
   !>
   !> @param[in,out] self Instance
   !> @param[in] mol Molecular structure
   !> @param[out] error Error handling
   subroutine continuum_update(self, mol, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      call self%invalidate()
      call clear_component_energies(self)
      self%configured = .true.
      if (.not. allocated(self%cavity)) then
         call fatal_error(error, "Continuum model has no cavity")
         return
      end if
      do i = 1, size(self%components)
         call self%components(i)%item%configure_cavity(self%cavity, error)
         if (allocated(error)) return
      end do
      call self%cavity%update(mol, error)
      if (allocated(error)) return
      ! A cavity waiting for the build phase leaves the components to it
      if (.not. self%cavity%defers_build()) then
         call update_components(self, error)
         if (allocated(error)) return
      end if
      self%updated = .true.
      self%built = .not. self%defers_build()

   end subroutine continuum_update

   !> Whether the cavity waits for the build phase
   !>
   !> @param[in] self Instance
   function continuum_defers_build(self) result(defers)
      !> Continuum model
      class(model_continuum_type), intent(in) :: self
      !> Whether `update` leaves the surface and the components to `build`
      logical :: defers

      defers = .false.
      if (allocated(self%cavity)) defers = self%cavity%defers_build()

   end function continuum_defers_build

   !> Read the build phase, then complete a cavity that waited for it and
   !> every component
   !>
   !> - the base protocol first; a model that `update` built completely
   !>   stays built
   !> - otherwise built only once the cavity and every component are, then
   !>   the coupling records the new point count
   !> - timed as `Build` after the protocol, also when `update` built
   !>   everything
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Coupling staged by `prepare_build`; the cavity reads its own scope
   !> @param[out] error Error handling
   subroutine continuum_build(self, coupling, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Coupling staged by `prepare_build`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call model_build(self, coupling, error)
      if (allocated(error)) return
      call self%ctx%timer%start("Build", category=cat_build)
      call complete_build(self, coupling, error)
      call self%ctx%timer%stop("Build")

   end subroutine continuum_build

   !> Complete a cavity that waited for the build phase, then every component
   !>
   !> @param[in,out] self Instance
   !> @param[in,out] coupling Coupling staged by `prepare_build`; the cavity reads its own scope
   !> @param[out] error Error handling
   subroutine complete_build(self, coupling, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Coupling staged by `prepare_build`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ! Built completely by `update`, components included
      if (.not. self%defers_build()) return
      call self%hvp_exchange%destroy()
      nullify (self%hvp_coupling)
      call clear_component_energies(self)
      call self%cavity%build(coupling, error)
      if (allocated(error)) return
      call update_components(self, error)
      if (allocated(error)) return
      self%built = .true.
      ! The answers of the build phase survive the new point count
      call self%declare_pass(coupling, error)

   end subroutine complete_build

   !> Update every component on the built cavity and its structure
   !>
   !> @param[in,out] self Instance
   !> @param[out] error Error handling
   subroutine update_components(self, error)
      !> Continuum model
      class(model_continuum_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      ! `cavity%mol` is passed beside `cavity`: components must not modify it
      do i = 1, size(self%components)
         call self%components(i)%item%update(self%cavity%mol, self%cavity, error)
         if (allocated(error)) return
      end do

   end subroutine update_components

   !> Accumulate the energy of every component
   !>
   !> Requires a coupling staged by `prepare_energy`; reports any other staged
   !> phase, then every missing output of the energy phase, by name
   !>
   !> - once the inputs are accepted, every component energy is cleared, then
   !>   each component stores its own on success
   !> - stops at the first failing component: earlier energies stay readable,
   !>   the failing and later ones stay unavailable, `energy` is unchanged
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
      !> Energy of one component, accumulated from zero
      real(wp) :: component_energy
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
      call clear_component_energies(self)
      local_energy = 0.0_wp
      do i = 1, size(self%components)
         component_energy = 0.0_wp
         call coupling_make_view(coupling, i, view)
         call self%components(i)%item%get_energy(view, self%cavity, component_energy, error)
         call coupling_close_view(coupling)
         if (allocated(error)) return
         self%components(i)%energy = component_energy
         local_energy = local_energy + component_energy
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
   !> - in forward mode every component emits it from `get_forward_gradient`
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
            call self%components(i)%item%get_forward_gradient(view, self%cavity, response, local, &
                                                              error)
            call coupling_close_view(coupling)
            if (allocated(error)) return
         end do
      end if

      gradient = gradient + local

   end subroutine continuum_get_gradient

   !* ================================================================================= *!
   !*                                Second derivatives                                 *!
   !* ================================================================================= *!

   !> Clear cached HVP state alongside model invalidation
   subroutine continuum_invalidate(self)
      class(model_continuum_type), intent(inout) :: self

      call self%hvp_exchange%destroy()
      nullify (self%hvp_coupling)
      call invalidate_model(self)
   end subroutine continuum_invalidate

   !> Clear cached HVP state if its coupling is released
   subroutine continuum_release_coupling(self, coupling)
      class(model_continuum_type), intent(inout), target :: self
      type(coupling_type), pointer, intent(inout) :: coupling

      if (associated(coupling)) then
         if (associated(self%hvp_coupling, coupling)) then
            call self%hvp_exchange%destroy()
            nullify (self%hvp_coupling)
         end if
      end if
      call model_release_coupling(self, coupling)
   end subroutine continuum_release_coupling

   !> Check that the cavity and every component support reverse-mode nuclear Hessians
   subroutine require_hessian_ready(self, error)
      class(model_continuum_type), intent(in) :: self
      type(error_type), allocatable, intent(out) :: error

      character(len=:), allocatable :: label
      integer :: i

      call self%require_built(error)
      if (allocated(error)) return
      if (self%force_forward_gradient) then
         call fatal_error(error, "Continuum-model Hessian is the derivative of the reverse-mode"// &
                          " gradient; the forward gradient path has none")
         return
      end if
      if (.not. self%cavity%supports_hessian()) then
         call fatal_error(error, "This cavity does not support nuclear Hessians")
         return
      end if
      do i = 1, size(self%components)
         if (.not. self%components(i)%item%supports_hessian()) then
            label = "This"
            if (allocated(self%components(i)%item%name)) then
               label = trim(self%components(i)%item%name)
            end if
            call fatal_error(error, label//" component does not support nuclear Hessians")
            return
         end if
      end do
   end subroutine require_hessian_ready

   !> Accumulate primal reverse-mode gradient surface weights across all components
   subroutine accumulate_gradient_surface_weights(self, coupling, acc, error)
      class(model_continuum_type), intent(inout) :: self
      type(coupling_type), intent(inout), target :: coupling
      type(cavity_surface_adjoint_type), intent(out) :: acc
      type(error_type), allocatable, intent(out) :: error

      type(coupling_view_type) :: view
      integer :: i

      call acc%init(self%cavity%ngrid)
      do i = 1, size(self%components)
         call coupling_make_view(coupling, i, view)
         call self%components(i)%item%get_gradient_surface_weights(view, self%cavity, acc, error)
         call coupling_close_view(coupling)
         if (allocated(error)) return
      end do
   end subroutine accumulate_gradient_surface_weights

   !> Stage HVP stage 1: host level-set/density tangents at fixed base surface
   subroutine continuum_prepare_hvp_geometry(self, coupling, dirs, error)
      class(model_continuum_type), intent(inout), target :: self
      type(coupling_type), intent(inout), target :: coupling
      real(wp), intent(in) :: dirs(:, :, :)
      type(error_type), allocatable, intent(out) :: error

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call require_hessian_ready(self, error)
      if (allocated(error)) return
      if (size(dirs, 1) /= 3 .or. size(dirs, 2) /= self%cavity%nsph .or. size(dirs, 3) < 1) then
         call fatal_error(error, "Continuum-model HVP direction shape mismatch")
         return
      end if

      call self%hvp_exchange%destroy()
      nullify (self%hvp_coupling)
      call coupling_set_hvp_dirs(coupling, dirs)
      call self%declare_pass(coupling, error)
      if (allocated(error)) return
      call coupling_arm(coupling, moist_phase_hvp_geometry, error)
      if (allocated(error)) return
      self%hvp_coupling => coupling
   end subroutine continuum_prepare_hvp_geometry

   !> Stage HVP stage 2: compute surface tangents and stage potential/trace tangents
   subroutine continuum_prepare_hvp(self, coupling, error)
      class(model_continuum_type), intent(inout), target :: self
      type(coupling_type), intent(inout), target :: coupling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: dirs(:, :, :)
      logical :: want_curv
      integer :: ngrid, ndir, i

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call require_hessian_ready(self, error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_hvp_geometry, error)
      if (allocated(error)) return
      if (.not. associated(self%hvp_coupling, coupling)) then
         call fatal_error(error, "prepare_hvp requires the coupling staged by prepare_hvp_geometry")
         return
      end if
      call coupling_get_hvp_dirs(coupling, dirs)
      if (.not. allocated(dirs)) then
         call fatal_error(error, "prepare_hvp: missing HVP directions on coupling")
         return
      end if

      call self%cavity%read_hvp_geometry(coupling, self%hvp_exchange, error)
      if (allocated(error)) return

      ngrid = self%cavity%ngrid
      ndir = size(dirs, 3)
      want_curv = allocated(self%cavity%k1) .and. allocated(self%cavity%k2)
      call self%hvp_exchange%tangent%init(ngrid, ndir, want_curv)
      if (allocated(self%hvp_exchange%d_wbranch)) deallocate (self%hvp_exchange%d_wbranch)
      allocate (self%hvp_exchange%d_wbranch(ngrid, ndir), source=0.0_wp)
      if (ngrid > 0) then
         if (self%hvp_exchange%have_host_jets) then
            call self%cavity%get_surface_tangent(dirs, self%hvp_exchange%tangent, error, &
               & d_wbranch=self%hvp_exchange%d_wbranch, host_jets=self%hvp_exchange%host_jets)
         else
            call self%cavity%get_surface_tangent(dirs, self%hvp_exchange%tangent, error, &
               & d_wbranch=self%hvp_exchange%d_wbranch)
         end if
         if (allocated(error)) then
            call self%hvp_exchange%destroy()
            nullify (self%hvp_coupling)
            return
         end if
      end if
      self%hvp_exchange%have_tangent = .true.

      if (allocated(self%cavity%xi0)) then
         call coupling_set_hvp_surface_tangent(coupling, self%hvp_exchange%tangent%d_xyz, &
            & self%hvp_exchange%tangent%d_xi)
      else
         call coupling_set_hvp_surface_tangent(coupling, self%hvp_exchange%tangent%d_xyz)
      end if
      call self%declare_pass(coupling, error)
      if (allocated(error)) return
      do i = 1, size(self%components)
         call self%components(i)%item%reset_hvp(ngrid, ndir)
      end do
      call coupling_arm(coupling, moist_phase_hvp, error)
   end subroutine continuum_prepare_hvp

   !> Evaluate nuclear Hessian-vector products and second-order response tangents
   subroutine continuum_get_hvp(self, coupling, response, hvp, error)
      class(model_continuum_type), intent(inout), target :: self
      type(coupling_type), intent(inout), target :: coupling
      type(response_type), intent(inout) :: response
      real(wp), intent(inout) :: hvp(:, :, :)
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: dirs(:, :, :), local(:, :, :)
      type(response_type) :: local_resp
      type(cavity_surface_adjoint_type) :: acc
      type(continuum_adjoint_response_type) :: omega
      type(coupling_view_type) :: view
      integer :: ngrid, ndir, i, d0

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call require_hessian_ready(self, error)
      if (allocated(error)) return
      call coupling_check_mandatory(coupling, moist_phase_hvp, error)
      if (allocated(error)) return
      if (.not. associated(self%hvp_coupling, coupling) .or. .not. self%hvp_exchange%have_tangent) then
         call fatal_error(error, "get_hvp requires a coupling prepared by prepare_hvp")
         return
      end if
      call coupling_get_hvp_dirs(coupling, dirs)
      if (.not. allocated(dirs)) then
         call fatal_error(error, "get_hvp: missing HVP directions on coupling")
         return
      end if
      if (any(shape(hvp) /= shape(dirs))) then
         call fatal_error(error, "Continuum-model Hessian-vector product shape mismatch")
         return
      end if

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("Model HVP", category=cat_hessian)

      ngrid = self%cavity%ngrid
      ndir = size(dirs, 3)
      do i = 1, size(self%components)
         call self%components(i)%item%reset_hvp(ngrid, ndir)
      end do
      call self%hvp_exchange%init_outputs(ngrid, ndir, self%cavity%has_density_response())

      call accumulate_gradient_surface_weights(self, coupling, acc, error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if

      omega%components => self%components
      omega%coupling => coupling
      allocate (local, mold=hvp)
      local = 0.0_wp
      call self%cavity%get_surface_hessian(acc, dirs, local, error, omega_v=omega, &
         & exchange=self%hvp_exchange)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if

      if (ngrid > 0) then
         do i = 1, size(self%components)
            call coupling_make_view(coupling, i, view)
            call self%components(i)%item%get_trace_response_tangent(view, self%cavity, &
               & self%hvp_exchange%tangent, local_resp, error)
            call coupling_close_view(coupling)
            if (allocated(error)) then
               call self%ctx%timer%unwind(d0)
               return
            end if
         end do
         call self%cavity%get_surface_response_tangent(self%hvp_exchange, local_resp, error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
      end if

      call self%ctx%timer%stop("Model HVP")

      call response_clear(response)
      response = local_resp
      hvp = hvp + local
   end subroutine continuum_get_hvp

   !> Evaluate the dense nuclear Hessian of a host-independent continuum model
   subroutine continuum_get_hessian(self, coupling, hessian, error)
      class(model_continuum_type), intent(inout), target :: self
      type(coupling_type), intent(inout), target :: coupling
      real(wp), intent(inout) :: hessian(:, :, :, :)
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: local(:, :, :, :)
      type(cavity_surface_adjoint_type) :: acc
      type(continuum_adjoint_response_type) :: omega
      integer :: i, d0

      call self%require_owned(coupling, error)
      if (allocated(error)) return
      call require_hessian_ready(self, error)
      if (allocated(error)) return
      if (any(shape(hessian) /= [3, self%cavity%nsph, 3, self%cavity%nsph])) then
         call fatal_error(error, "Continuum-model Hessian shape mismatch")
         return
      end if
      if (self%cavity%has_density_response()) then
         call fatal_error(error, "Density-dependent cavities require the two-stage HVP"// &
                          " protocol (prepare_hvp_geometry / prepare_hvp / get_hvp)")
         return
      end if
      do i = 1, size(self%components)
         if (self%components(i)%item%needs_hvp_coupling()) then
            call fatal_error(error, trim(self%components(i)%item%name)// &
                             " requires the two-stage HVP protocol"// &
                             " (prepare_hvp_geometry / prepare_hvp / get_hvp)")
            return
         end if
      end do

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("Model Hessian", category=cat_hessian)

      call accumulate_gradient_surface_weights(self, coupling, acc, error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      omega%components => self%components
      omega%coupling => coupling
      allocate (local(3, self%cavity%nsph, 3, self%cavity%nsph), source=0.0_wp)
      call self%cavity%get_hessian(acc, local, error, omega_v=omega)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if

      call self%ctx%timer%stop("Model Hessian")

      hessian = hessian + local
   end subroutine continuum_get_hessian

   !> Apply component surface-weight tangents and direct nuclear HVP for one direction block
   subroutine continuum_adjoint_response_apply(self, cavity, dirs, tangent, dacc, hvp_direct, &
                                               error, first)
      class(continuum_adjoint_response_type), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      real(wp), intent(in) :: dirs(:, :, :)
      type(cavity_surface_tangent_type), intent(in) :: tangent
      type(cavity_surface_adjoint_type), intent(inout) :: dacc(:)
      real(wp), intent(inout) :: hvp_direct(:, :, :)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: first

      type(coupling_view_type) :: view
      integer :: i

      if (.not. associated(self%components) .or. .not. associated(self%coupling)) then
         call fatal_error(error, "Continuum-model adjoint response is not bound")
         return
      end if
      do i = 1, size(self%components)
         call coupling_make_view(self%coupling, i, view)
         call self%components(i)%item%get_hessian_surface_weights(view, cavity, dirs, tangent, &
                                                                  dacc, error, first=first)
         call coupling_close_view(self%coupling)
         if (allocated(error)) return
         call coupling_make_view(self%coupling, i, view)
         call self%components(i)%item%get_direct_hessian(view, cavity, dirs, tangent, &
                                                         hvp_direct, error, first=first)
         call coupling_close_view(self%coupling)
         if (allocated(error)) return
      end do
   end subroutine continuum_adjoint_response_apply

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Declare cavity and component requests, then record the extents
   !>
   !> - calculations with matching inputs are shared between components
   !> - until the model and its cavity are built only the cavity declares:
   !>   the requests of the build phase belong to its level set
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
      !> Counts the snapshot records; unallocated ones are absent
      integer, allocatable :: ngrid, natom

      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 0)
      call self%cavity%declare_coupling(coupling, error)
      if (allocated(error)) return
      if (self%is_built() .and. self%cavity%is_built()) then
         do i = 1, size(self%components)
            call coupling_set_scope(coupling, i)
            call self%components(i)%item%declare_coupling(self%cavity, coupling, error)
            if (allocated(error)) return
         end do
      end if
      call self%cavity%snapshot_counts(ngrid, natom)
      call coupling_snapshot(coupling, ngrid=ngrid, natom=natom)

   end subroutine declare_continuum_pass

   !> Release collections owned by the model; borrowed handles must not outlive it
   !>
   !> @param[in,out] self Instance
   subroutine destroy_model(self)
      !> Model being destroyed
      type(model_continuum_type), intent(inout) :: self
      call self%hvp_exchange%destroy()
      nullify (self%hvp_coupling)
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

   !> Publish the fields of the owned cavity, its evaluation domain
   !>
   !> @param[in] self Model
   !> @param[in,out] query Field walker
   subroutine continuum_list_fields(self, query)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query
      if (allocated(self%cavity)) call self%cavity%list_fields(query)
   end subroutine continuum_list_fields

   !> Number of components in the order they were added
   !>
   !> @param[in] self Model
   function continuum_component_count(self) result(ncomponent)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> Component count
      integer :: ncomponent
      ncomponent = 0
      if (allocated(self%components)) ncomponent = size(self%components)
   end function continuum_component_count

   !> Name of one component; repeated names are told apart by index
   !>
   !> A subroutine: gfortran shares the length of a deferred-length function
   !> result between OpenMP threads
   !>
   !> @param[in] self Model
   !> @param[in] index 1-based component index
   !> @param[out] name Component name, unallocated for an index outside the list
   subroutine continuum_component_name(self, index, name)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> 1-based component index
      integer, intent(in) :: index
      !> Component name
      character(len=:), allocatable, intent(out) :: name

      if (index < 1 .or. index > self%component_count()) return
      if (allocated(self%components(index)%item%name)) then
         name = self%components(index)%item%name
      else
         name = ""
      end if

   end subroutine continuum_component_name

   !> One-line description of one component
   !>
   !> A subroutine for the same reason as `continuum_component_name`
   !>
   !> @param[in] self Model
   !> @param[in] index 1-based component index
   !> @param[out] description Component description, unallocated for an index outside the list
   subroutine continuum_component_description(self, index, description)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> 1-based component index
      integer, intent(in) :: index
      !> Component description
      character(len=:), allocatable, intent(out) :: description

      if (index < 1 .or. index > self%component_count()) return
      if (allocated(self%components(index)%item%description)) then
         description = self%components(index)%item%description
      else
         description = ""
      end if

   end subroutine continuum_component_description

   !> Print the settings of the cavity and of every component
   !>
   !> - the cavity section first, then the components in list order, each
   !>   headed by its 1-based index
   !> - available before the first update
   !>
   !> @param[in] self Model
   !> @param[in] unit Output unit; defaults to the run context's unit
   !> @param[in] level Verbosity level; settings above it are left out, nothing
   !>                  prints at 0 and 1 prints the normal settings (default: each part's context
   !>                  verbosity)
   subroutine continuum_print_parameters(self, unit, level)
      !> Model
      class(model_continuum_type), intent(in), target :: self
      !> Output unit
      integer, intent(in), optional :: unit
      !> Verbosity level of the printout
      integer, intent(in), optional :: level

      !> Effective output unit
      integer :: iu
      !> Component index
      integer :: i

      iu = output_unit
      if (associated(self%ctx)) iu = self%ctx%unit
      if (present(unit)) iu = unit
      if (allocated(self%cavity)) call self%cavity%print_parameters(iu, level)
      do i = 1, self%component_count()
         call self%components(i)%item%print_parameters(iu, i, level)
      end do

   end subroutine continuum_print_parameters

   !> Print the radii of the current structure, the cavity's radius model
   !>
   !> - nothing before a successful update; a cavity waiting for the build
   !>   phase already has its radii
   !> - the cavity follows its own context, like `print_parameters`
   !>
   !> @param[in] self Model
   !> @param[in] unit Output unit; defaults to the run context's unit
   !> @param[in] level Verbosity level; nothing prints at 0 (default: the
   !>                  cavity's context verbosity)
   subroutine continuum_print_radii(self, unit, level)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> Output unit
      integer, intent(in), optional :: unit
      !> Verbosity level of the printout
      integer, intent(in), optional :: level

      !> Effective output unit
      integer :: iu

      if (.not. self%updated .or. .not. allocated(self%cavity)) return
      iu = output_unit
      if (associated(self%ctx)) iu = self%ctx%unit
      if (present(unit)) iu = unit
      call self%cavity%print_radii(iu, level)

   end subroutine continuum_print_radii

   !> Publish the results of one component
   !>
   !> - `energy`, the scalar of the latest evaluation; not listed before it,
   !>   after an update and when this component failed or was not reached
   !> - then the component's own results (`list_fields`), under the same rule
   !> - an index outside the list declares nothing
   !>
   !> @param[in] self Model
   !> @param[in] index 1-based component index
   !> @param[in,out] query Field walker
   subroutine continuum_list_component_fields(self, index, query)
      !> Model
      class(model_continuum_type), intent(in) :: self
      !> 1-based component index
      integer, intent(in) :: index
      !> Field walker
      type(field_query_type), intent(inout) :: query

      if (index < 1 .or. index > self%component_count()) return
      call query%add_real_scalar("energy", "Energy contribution of the latest evaluation, in Hartree", &
         & self%components(index)%energy)
      if (allocated(self%components(index)%energy)) call self%components(index)%item%list_fields(query)

   end subroutine continuum_list_component_fields

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

end module moist_model_continuum_type
