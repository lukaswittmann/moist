!> General list-driven solvation model
module moist_model_general
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_type, only: solvation_model_type, solvation_model_component_type, cavity_type, &
      & surface_adjoint_response_type
   use moist_channels, only: coupling_type, response_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_cavity_surface_tangent, only: cavity_surface_tangent_type

   implicit none
   private

   public :: solvation_model_general, new_model_general

   !> Owning box that makes heterogeneous components storable in one array
   type :: solvation_component_slot
      !> Concrete component owned by this slot
      class(solvation_model_component_type), allocatable :: item
   end type solvation_component_slot

   !> Second-order surface channel of the general model
   !>
   !> Handed to the cavity's Hessian entry points, which call it back once per
   !> block of nuclear directions with the surface tangent of that block. It
   !> walks the component list exactly as the gradient does: every component
   !> states the response of its gradient-path surface adjoints, and any
   !> non-surface Hessian columns, and the cavity contracts the lot once.
   type, extends(surface_adjoint_response_type) :: general_adjoint_response_type
      !> The model's components, borrowed for the duration of one Hessian call
      type(solvation_component_slot), pointer :: components(:) => null()
      !> Host coupling data, borrowed likewise
      class(coupling_type), pointer :: coupling => null()
   contains
      procedure :: apply => general_adjoint_response_apply
   end type general_adjoint_response_type

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
      procedure :: get_hessian => general_get_hessian
      procedure :: get_hvp => general_get_hvp
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
      grown(n + 1)%item%ctx => self%ctx
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

   !> Dense nuclear Hessian of every component
   !>
   !> The derivative of [[general_get_gradient]] at fixed host coupling data:
   !> the same surface adjoints are accumulated, and the cavity differentiates
   !> their contraction with the components' adjoint response folded in through
   !> [[general_adjoint_response_type]]. The result is *added* to `hessian`, and
   !> the accumulator is left untouched when anything fails.
   !>
   !> @param[inout] self     General model
   !> @param[in]    coupling Host coupling data
   !> @param[inout] hessian  Nuclear-Hessian accumulator (3, nat, 3, nat)
   !> @param[out]   error    Error handling
   subroutine general_get_hessian(self, coupling, hessian, error)
      !> General model
      class(solvation_model_general), intent(inout), target :: self
      !> Host coupling data
      class(coupling_type), intent(in), target :: coupling
      !> Nuclear-Hessian accumulator
      real(wp), intent(inout) :: hessian(:, :, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional local Hessian
      real(wp), allocatable :: local(:, :, :, :)
      !> Shared surface accumulator and the components' response to it
      type(cavity_surface_adjoint_type) :: acc
      type(general_adjoint_response_type) :: omega

      call require_hessian_ready(self, error)
      if (allocated(error)) return
      if (any(shape(hessian) /= [3, self%cavity%nsph, 3, self%cavity%nsph])) then
         call fatal_error(error, "General-model Hessian shape mismatch")
         return
      end if

      call gradient_surface_weights(self, coupling, acc, error)
      if (allocated(error)) return
      omega%components => self%components
      omega%coupling => coupling

      allocate (local(3, self%cavity%nsph, 3, self%cavity%nsph), source=0.0_wp)
      call self%cavity%get_hessian(acc, local, error, omega_v=omega)
      if (allocated(error)) return

      hessian = hessian + local

   end subroutine general_get_hessian

   !> Nuclear Hessian-vector products of every component
   !>
   !> One column of [[general_get_hessian]] per supplied direction, without
   !> the dense block: the cavity runs its per-direction channels and the
   !> components' adjoint response for exactly these directions. The result is
   !> *added* to `hvp`, and the accumulator is left untouched when anything
   !> fails.
   !>
   !> @param[inout] self     General model
   !> @param[in]    coupling Host coupling data
   !> @param[in]    dirs     Nuclear directions (3, nat, ndir)
   !> @param[inout] hvp      Hessian-vector accumulator (3, nat, ndir)
   !> @param[out]   error    Error handling
   subroutine general_get_hvp(self, coupling, dirs, hvp, error)
      !> General model
      class(solvation_model_general), intent(inout), target :: self
      !> Host coupling data
      class(coupling_type), intent(in), target :: coupling
      !> Nuclear directions
      real(wp), intent(in) :: dirs(:, :, :)
      !> Hessian-vector accumulator
      real(wp), intent(inout) :: hvp(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Transactional local columns
      real(wp), allocatable :: local(:, :, :)
      !> Shared surface accumulator and the components' response to it
      type(cavity_surface_adjoint_type) :: acc
      type(general_adjoint_response_type) :: omega

      call require_hessian_ready(self, error)
      if (allocated(error)) return
      if (size(dirs, 1) /= 3 .or. size(dirs, 2) /= self%cavity%nsph .or. size(dirs, 3) < 1) then
         call fatal_error(error, "General-model Hessian-vector product direction shape mismatch")
         return
      end if
      if (any(shape(hvp) /= shape(dirs))) then
         call fatal_error(error, "General-model Hessian-vector product shape mismatch")
         return
      end if

      call gradient_surface_weights(self, coupling, acc, error)
      if (allocated(error)) return
      omega%components => self%components
      omega%coupling => coupling

      allocate (local, mold=hvp)
      local = 0.0_wp
      call self%cavity%get_surface_hessian(acc, dirs, local, error, omega_v=omega)
      if (allocated(error)) return

      hvp = hvp + local

   end subroutine general_get_hvp

   !> The surface adjoints the nuclear gradient contracts, from every component
   !>
   !> The Hessian differentiates exactly this accumulation, so it is built by
   !> the same hooks the reverse-mode gradient uses and by no others.
   !>
   !> @param[inout] self     General model
   !> @param[in]    coupling Host coupling data
   !> @param[out]   acc      Surface-adjoint accumulator
   !> @param[out]   error    Error handling
   subroutine gradient_surface_weights(self, coupling, acc, error)
      !> General model
      class(solvation_model_general), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Surface-adjoint accumulator
      type(cavity_surface_adjoint_type), intent(out) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      call acc%init(self%cavity%ngrid)
      do i = 1, size(self%components)
         call self%components(i)%item%get_gradient_surface_weights(coupling, self%cavity, &
                                                                   acc, error)
         if (allocated(error)) return
      end do

   end subroutine gradient_surface_weights

   !> Require a model whose Hessian is defined
   !>
   !> The Hessian is the derivative of the reverse-mode gradient; the legacy
   !> forward path has no second-order counterpart.
   !>
   !> @param[in]  self  General model
   !> @param[out] error Error handling
   subroutine require_hessian_ready(self, error)
      !> General model
      class(solvation_model_general), intent(in) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call require_updated(self, error)
      if (allocated(error)) return
      if (self%force_forward_gradient) then
         call fatal_error(error, "General-model Hessian is the derivative of the reverse-mode"// &
                          " gradient; the forward gradient path has none")
         return
      end if

   end subroutine require_hessian_ready

   !> Adjoint response of every component along one block of directions
   !>
   !> @param[inout] self       Response object
   !> @param[in]    cavity     Cavity the tangent was taken on
   !> @param[in]    dirs       Nuclear directions of the block (3, nsph, nblk)
   !> @param[in]    tangent    Surface tangent of the block
   !> @param[inout] dacc       Surface-adjoint response per direction (nblk)
   !> @param[inout] hvp_direct Non-surface Hessian columns of the block
   !> @param[out]   error      Error handling
   subroutine general_adjoint_response_apply(self, cavity, dirs, tangent, dacc, hvp_direct, &
                                             error)
      !> Response object
      class(general_adjoint_response_type), intent(inout) :: self
      !> Cavity the tangent was taken on
      class(cavity_type), intent(in) :: cavity
      !> Nuclear directions of the block
      real(wp), intent(in) :: dirs(:, :, :)
      !> Surface tangent of the block
      type(cavity_surface_tangent_type), intent(in) :: tangent
      !> Surface-adjoint response per direction
      type(cavity_surface_adjoint_type), intent(inout) :: dacc(:)
      !> Non-surface Hessian columns of the block
      real(wp), intent(inout) :: hvp_direct(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Component index
      integer :: i

      if (.not. associated(self%components) .or. .not. associated(self%coupling)) then
         call fatal_error(error, "General-model adjoint response was not bound to a model")
         return
      end if

      do i = 1, size(self%components)
         call self%components(i)%item%get_hessian_surface_weights(self%coupling, cavity, dirs, &
                                                                  tangent, dacc, error)
         if (allocated(error)) return
         call self%components(i)%item%get_direct_hessian(self%coupling, cavity, dirs, tangent, &
                                                         hvp_direct, error)
         if (allocated(error)) return
      end do

   end subroutine general_adjoint_response_apply

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
