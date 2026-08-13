!> General list-driven solvation model
module moist_model_general
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_type, only: solvation_model_type, solvation_model_component_type, cavity_type
   use moist_channels, only: coupling_type, response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type

   implicit none
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
      !> Force the legacy forward nuclear-gradient path
      logical :: force_forward_gradient = .false.
   contains
      procedure :: add_component
      procedure :: update => general_update
      procedure :: get_trace_response => general_get_trace_response
      procedure :: get_energy => general_get_energy
      procedure :: get_response => general_get_response
      procedure :: get_gradient => general_get_gradient
   end type solvation_model_general

contains

   !> Construct an empty general model around an owned copy of a cavity
   !>
   !> @param[out] self   General model
   !> @param[in]  cavity Cavity configuration to copy
   !> @param[in]  ctx    Shared run context, which must outlive the model
   !> @param[out] error  Error handling
   subroutine new_model_general(self, cavity, ctx, error)
      !> General model
      type(solvation_model_general), intent(out) :: self
      !> Cavity configuration to copy
      class(cavity_type), intent(in) :: cavity
      !> Shared run context
      type(moist_context_type), intent(in), target :: ctx
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      allocate (self%cavity, source=cavity)
      self%ctx => ctx
      self%cavity%ctx => ctx
      allocate (self%components(0))
      self%updated = .false.
      if (.not. allocated(self%cavity)) call fatal_error(error, "Failed to copy model cavity")

   end subroutine new_model_general

   !> Append an owned copy of a component before the first update
   !>
   !> @param[inout] self      General model
   !> @param[in]    component Component to copy into the model
   !> @param[out]   error     Error handling
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

      if (self%updated) then
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
      call move_alloc(grown, self%components)

   end subroutine add_component

   !> Update the authoritative cavity and every component
   !>
   !> @param[inout] self  General model
   !> @param[in]    mol   Molecular structure
   !> @param[out]   error Error handling
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

   !> Return direct host-trace adjoints before charge-dependent host response
   !>
   !> @param[inout] self      General model
   !> @param[in]    coupling  Host coupling data
   !> @param[inout] response  Response accumulator
   !> @param[out]   error     Error handling
   subroutine general_get_trace_response(self, coupling, response, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Potential accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional local response
      type(response_type) :: local
      !> Component index
      integer :: i

      call require_updated(self, error)
      if (allocated(error)) return
      do i = 1, size(self%components)
         call self%components(i)%item%get_trace_response(coupling, self%cavity, local, error)
         if (allocated(error)) return
      end do
      call add_response(response, local, error)

   end subroutine general_get_trace_response

   !> Accumulate the energy of every component
   !>
   !> @param[inout] self     General model
   !> @param[in]    coupling Host coupling data
   !> @param[inout] energy   Energy accumulator
   !> @param[out]   error    Error handling
   subroutine general_get_energy(self, coupling, energy, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Energy accumulator
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional energy accumulator
      real(wp) :: local_energy
      !> Component index
      integer :: i

      call require_updated(self, error)
      if (allocated(error)) return
      local_energy = 0.0_wp
      do i = 1, size(self%components)
         call self%components(i)%item%get_energy(coupling, self%cavity, local_energy, error)
         if (allocated(error)) return
      end do
      energy = energy + local_energy

   end subroutine general_get_energy

   !> Assemble direct trace and cavity-response channels.
   !> @param[inout] self      General model
   !> @param[in]    coupling  Host coupling data
   !> @param[inout] response  Response accumulator
   !> @param[out]   error     Error handling
   subroutine general_get_response(self, coupling, response, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Potential accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional local response
      type(response_type) :: local
      !> Shared surface accumulator
      type(cavity_surface_adjoint_type) :: acc
      !> Component index
      integer :: i

      call require_updated(self, error)
      if (allocated(error)) return
      call acc%init(self%cavity%ngrid)
      do i = 1, size(self%components)
         call self%components(i)%item%get_response(coupling, self%cavity, local, error)
         if (allocated(error)) return
         call self%components(i)%item%get_surface_weights(coupling, self%cavity, acc, error)
         if (allocated(error)) return
      end do
      call self%cavity%get_surface_response(acc, local, error)
      if (allocated(error)) return
      call add_response(response, local, error)

   end subroutine general_get_response

   !> Accumulate the nuclear gradient of every component
   !>
   !> @param[inout] self     General model
   !> @param[in]    coupling Host coupling data
   !> @param[inout] gradient Nuclear-gradient accumulator
   !> @param[out]   error    Error handling
   subroutine general_get_gradient(self, coupling, gradient, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional local gradient
      real(wp), allocatable :: local(:, :)
      !> Component index
      integer :: i

      call require_updated(self, error)
      if (allocated(error)) return
      if (any(shape(gradient) /= [3, self%cavity%nsph])) then
         call fatal_error(error, "General-model gradient shape mismatch")
         return
      end if
      allocate (local(3, self%cavity%nsph), source=0.0_wp)

      if ( .not. self%force_forward_gradient) then
         ! Reverse mode: every component states its surface adjoint, the cavity
         ! contracts the lot once. Nothing builds a nuclear Jacobian.
         block
            type(cavity_surface_adjoint_type) :: acc

            call acc%init(self%cavity%ngrid)
            do i = 1, size(self%components)
               call self%components(i)%item%get_direct_gradient(coupling, self%cavity, &
                                                                local, error)
               if (allocated(error)) return
               call self%components(i)%item%get_gradient_surface_weights(coupling, &
                                                                         self%cavity, acc, error)
               if (allocated(error)) return
            end do
            call self%cavity%get_surface_gradient(acc, local, error)
            if (allocated(error)) return
         end block
      else
         do i = 1, size(self%components)
            call self%components(i)%item%get_gradient(coupling, self%cavity, local, error)
            if (allocated(error)) return
         end do
      end if

      gradient = gradient + local

   end subroutine general_get_gradient

   !> Require a successfully updated model.
   !> @param[in]  self  General model
   !> @param[out] error Error handling
   subroutine require_updated(self, error)
      !> General model
      class(solvation_model_general), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%updated) call fatal_error(error, "General model must be updated first")

   end subroutine require_updated

   !> Add every allocated channel from one response to another
   !>
   !> @param[inout] target Destination accumulator
   !> @param[in]    source Source contribution
   !> @param[out]   error  Error handling
   subroutine add_response(target, source, error)
      !> Destination accumulator
      type(response_type), intent(inout) :: target
      !> Source contribution
      type(response_type), intent(in) :: source
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call add_vector(target%electrostatics%surface_charge, &
         & source%electrostatics%surface_charge, "electrostatics%surface_charge", error)
      if (allocated(error)) return
      call add_vector(target%gostshyp%w_overlap, source%gostshyp%w_overlap, &
         & "gostshyp%w_overlap", error)
      if (allocated(error)) return
      call add_vector(target%gostshyp%w_normal_deriv, source%gostshyp%w_normal_deriv, &
         & "gostshyp%w_normal_deriv", error)
      if (allocated(error)) return
      call add_vector(target%lsf%w_value, source%lsf%w_value, "lsf%w_value", error)
      if (allocated(error)) return
      call add_matrix(target%lsf%w_gradient, source%lsf%w_gradient, "lsf%w_gradient", error)
      if (allocated(error)) return
      call add_tensor3(target%lsf%w_hessian, source%lsf%w_hessian, "lsf%w_hessian", error)

   end subroutine add_response

   !> Add an allocated vector contribution to a response channel
   !>
   !> @param[inout] target Destination accumulator
   !> @param[in]    source Source contribution
   !> @param[in]    name   Channel name used in diagnostics
   !> @param[out]   error  Error handling
   subroutine add_vector(target, source, name, error)
      !> Destination vector
      real(wp), allocatable, intent(inout) :: target(:)
      !> Source vector
      real(wp), allocatable, intent(in) :: source(:)
      !> Channel name
      character(len=*), intent(in) :: name
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(source)) return
      if (.not. allocated(target)) then
         allocate (target, source=source)
      else if (size(target) /= size(source)) then
         call fatal_error(error, "Potential shape mismatch for "//name)
      else
         target = target + source
      end if

   end subroutine add_vector

   !> Add an allocated matrix contribution to a response channel
   !>
   !> @param[inout] target Destination accumulator
   !> @param[in]    source Source contribution
   !> @param[in]    name   Channel name used in diagnostics
   !> @param[out]   error  Error handling
   subroutine add_matrix(target, source, name, error)
      !> Destination matrix
      real(wp), allocatable, intent(inout) :: target(:, :)
      !> Source matrix
      real(wp), allocatable, intent(in) :: source(:, :)
      !> Channel name
      character(len=*), intent(in) :: name
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(source)) return
      if (.not. allocated(target)) then
         allocate (target, source=source)
      else if (any(shape(target) /= shape(source))) then
         call fatal_error(error, "Potential shape mismatch for "//name)
      else
         target = target + source
      end if

   end subroutine add_matrix

   !> Add an allocated rank-three contribution to a response channel
   !>
   !> @param[inout] target Destination accumulator
   !> @param[in]    source Source contribution
   !> @param[in]    name   Channel name used in diagnostics
   !> @param[out]   error  Error handling
   subroutine add_tensor3(target, source, name, error)
      !> Destination tensor
      real(wp), allocatable, intent(inout) :: target(:, :, :)
      !> Source tensor
      real(wp), allocatable, intent(in) :: source(:, :, :)
      !> Channel name
      character(len=*), intent(in) :: name
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(source)) return
      if (.not. allocated(target)) then
         allocate (target, source=source)
      else if (any(shape(target) /= shape(source))) then
         call fatal_error(error, "Potential shape mismatch for "//name)
      else
         target = target + source
      end if

   end subroutine add_tensor3

end module moist_model_general
