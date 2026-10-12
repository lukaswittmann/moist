!> Shared host integration for single- and double-layer PCM models
module moist_model_continuum_component_pcm_layer
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & coupling_request_type, point_potential_request_type, gaussian_potential_request_type, &
      & moist_phase_energy, moist_phase_response, moist_phase_gradient, &
      & coupling_register, request_require
   use moist_channels_response, only: response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, solver_type
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver, pcm_dense_solver, &
      & pcm_gmres_solver, pcm_block_jacobi
   implicit none(type, external)
   private

   public :: pcm_layer_base, pcm_layer_feps

   !> Common host lifecycle; each extension owns one generated model state
   type, abstract, extends(model_continuum_component_pcm) :: pcm_layer_base
      !> Symmetrized potential adjoint
      real(wp), allocatable :: qsym(:)
      !> Normal-derivative adjoint, for DPCM and SLIC
      real(wp), allocatable :: qn(:)
      !> Host normal-derivative trace
      real(wp), allocatable :: vn(:)
      !> Energy returned by the generated solve
      real(wp) :: kernel_energy = 0.0_wp
      !> True when the cavity has Gaussian widths and switching factors
      logical :: smeared = .false.
      !> This model reads the normal derivative of the host potential
      logical :: needs_normal = .false.
      !> Whether the generated operator needs principal curvatures
      logical :: needs_curvature = .true.
      !> Whether Gaussian kernels consume cavity switching factors
      logical :: needs_switching = .true.
   contains
      procedure :: configure_cavity => layer_configure_cavity
      procedure :: update => layer_update
      procedure :: ensure_charges => layer_ensure_charges
      procedure :: potential_adjoint => layer_potential_adjoint
      procedure :: normal_trace_adjoint => layer_normal_trace_adjoint
      procedure :: uses_normal_trace => layer_uses_normal_trace
      procedure :: get_energy => layer_get_energy
      procedure :: get_surface_weights => layer_get_surface_weights
      procedure :: get_gradient_surface_weights => layer_get_gradient_surface_weights
      procedure :: get_gradient => layer_get_gradient
      procedure :: declare_coupling => layer_declare_coupling
      procedure :: get_host_surface_weights => layer_get_host_surface_weights
      procedure(assemble_layer_iface), deferred :: assemble_layer
      procedure(solve_layer_iface), deferred :: solve_layer
      procedure(layer_weights_iface), deferred :: layer_weights
   end type pcm_layer_base

   abstract interface
      !> Assemble a generated model operator
      subroutine assemble_layer_iface(self, cavity, solver, error)
         import :: pcm_layer_base, cavity_type, pcm_solver, error_type
         implicit none(type, external)
         class(pcm_layer_base), intent(inout) :: self
         class(cavity_type), intent(in) :: cavity
         class(pcm_solver), intent(in) :: solver
         type(error_type), allocatable, intent(out) :: error
      end subroutine assemble_layer_iface

      !> Solve the generated model for current host traces
      subroutine solve_layer_iface(self, error)
         import :: pcm_layer_base, error_type
         implicit none(type, external)
         class(pcm_layer_base), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_layer_iface

      !> Add fixed-trace geometry weights from the generated kernel
      subroutine layer_weights_iface(self, cavity, acc, error)
         import :: pcm_layer_base, cavity_type, cavity_surface_adjoint_type, error_type
         implicit none(type, external)
         class(pcm_layer_base), intent(inout) :: self
         class(cavity_type), intent(in) :: cavity
         class(cavity_surface_adjoint_type), intent(inout) :: acc
         type(error_type), allocatable, intent(out) :: error
      end subroutine layer_weights_iface
   end interface

contains

   !> Local-layer dielectric factor, including the conductor limit
   pure function pcm_layer_feps(epsilon) result(f)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      real(wp), intent(in) :: epsilon
      real(wp) :: f
      if (ieee_is_finite(epsilon)) then
         f = (epsilon - 1.0_wp)/(epsilon + 1.0_wp)
      else
         f = 1.0_wp
      end if
   end function pcm_layer_feps

   !> Whether this concrete layer model reads the normal trace
   pure function layer_uses_normal_trace(self) result(uses)
      class(pcm_layer_base), intent(in) :: self
      logical :: uses
      uses = self%needs_normal
   end function layer_uses_normal_trace

   !> Register the potential primitive matching the cavity discretization
   subroutine layer_declare_coupling(self, cavity, coupling, error)
      class(pcm_layer_base), intent(in) :: self
      class(cavity_type), intent(in) :: cavity
      type(coupling_type), intent(inout) :: coupling
      type(error_type), allocatable, intent(out) :: error
      class(coupling_request_type), allocatable :: potential
      logical :: normal, field, gaussian

      normal = self%uses_normal_trace()
      field = cavity%has_field_dependent_geometry()
      gaussian = allocated(cavity%xi0)
      if (gaussian) then
         allocate (gaussian_potential_request_type :: potential)
      else
         allocate (point_potential_request_type :: potential)
      end if
      call request_require(potential, moist_phase_energy, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_energy, "dphi_dr", normal, error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "dphi_dr", field .or. normal, error)
      if (allocated(error)) return
      if (gaussian) then
         call request_require(potential, moist_phase_response, "dphi_dxi", field, error)
         if (allocated(error)) return
      end if
      call request_require(potential, moist_phase_response, "d2phi_dr2", normal .and. field, error)
      if (allocated(error)) return
      if (gaussian) then
         call request_require(potential, moist_phase_response, "d2phi_dr_dxi", normal .and. field, error)
         if (allocated(error)) return
      end if
      call request_require(potential, moist_phase_gradient, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "dphi_dr", error)
      if (allocated(error)) return
      if (gaussian) then
         call request_require(potential, moist_phase_gradient, "dphi_dxi", error)
         if (allocated(error)) return
      end if
      call request_require(potential, moist_phase_gradient, "d2phi_dr2", normal, error)
      if (allocated(error)) return
      if (gaussian) then
         call request_require(potential, moist_phase_gradient, "d2phi_dr_dxi", normal, error)
         if (allocated(error)) return
      end if
      call coupling_register(coupling, "potential", potential, error)
   end subroutine layer_declare_coupling

   !> Request normals and principal curvatures before cavity update
   subroutine layer_configure_cavity(self, cavity, error)
      class(pcm_layer_base), intent(inout) :: self
      class(cavity_type), intent(inout) :: cavity
      type(error_type), allocatable, intent(out) :: error
      if (self%needs_curvature) then
         call cavity%configure(normal=.true., curvature=.true., error=error)
      else
         call cavity%configure(normal=.true., error=error)
      end if
   end subroutine layer_configure_cavity

   !> Assemble the generated operators on the current cavity
   subroutine layer_update(self, mol, cavity, error)
      class(pcm_layer_base), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      class(cavity_type), intent(inout) :: cavity
      type(error_type), allocatable, intent(out) :: error
      class(pcm_solver), allocatable :: solver
      integer :: n

      call self%require_context(error)
      if (allocated(error)) return
      self%mol_solu = mol
      self%charges_valid = .false.
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%vn)) deallocate (self%vn)
      if (allocated(self%q)) deallocate (self%q)
      if (allocated(self%qsym)) deallocate (self%qsym)
      if (allocated(self%qn)) deallocate (self%qn)
      n = cavity%ngrid
      if (.not. allocated(cavity%xyz) .or. .not. allocated(cavity%a) .or. &
          & .not. allocated(cavity%normal0)) then
         call fatal_error(error, "Layer PCM requires positions, areas, and normals")
         return
      end if
      if (self%needs_curvature) then
         if (.not. allocated(cavity%k1) .or. .not. allocated(cavity%k2)) then
            call fatal_error(error, "Layer PCM requires principal curvatures")
            return
         end if
      end if
      if (size(cavity%xyz, 2) < n .or. size(cavity%normal0, 2) < n .or. &
          & size(cavity%a) < n) then
         call fatal_error(error, "Layer PCM cavity geometry is stale")
         return
      end if
      if (self%needs_curvature) then
         if (size(cavity%k1) < n .or. size(cavity%k2) < n) then
            call fatal_error(error, "Layer PCM cavity curvature is stale")
            return
         end if
      end if
      if (self%param%solver == solver_type%iterative) then
         allocate (pcm_gmres_solver :: solver)
         select type (solver)
         type is (pcm_gmres_solver)
            solver%tol = self%param%solver_tol
            solver%max_iter = self%param%solver_maxiter
            solver%restart = self%param%solver_restart
            solver%storage = self%param%storage
            solver%memory_budget_gib = self%param%memory_budget_gib
            allocate (solver%preconditioner, source=pcm_block_jacobi(block_size=256))
         end select
      else
         allocate (pcm_dense_solver :: solver)
      end if
      solver%nonlinear_tol = self%param%newton_tol
      solver%nonlinear_max_iter = self%param%newton_maxiter
      self%smeared = allocated(cavity%xi0)
      if (self%smeared .and. self%needs_switching) then
         if (.not. allocated(cavity%f)) then
            call fatal_error(error, "Layer PCM Gaussian surface requires switching factors")
            return
         end if
      end if
      call self%assemble_layer(cavity, solver, error)
      if (allocated(error)) return
      allocate (self%q(n), self%qsym(n), source=0.0_wp)
      if (self%uses_normal_trace()) allocate (self%qn(n), source=0.0_wp)
   end subroutine layer_update

   !> Read host traces and run the generated solve when they changed
   subroutine layer_ensure_charges(self, coupling, cavity, error)
      class(pcm_layer_base), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: phi(:), vn(:)
      logical :: same

      if (.not. allocated(self%q)) then
         call fatal_error(error, "Layer PCM must be updated before solving")
         return
      end if
      call coupling%read("potential", "phi", phi, error)
      if (allocated(error)) return
      if (size(phi) /= cavity%ngrid) then
         call fatal_error(error, "Layer PCM potential grid mismatch")
         return
      end if
      if (self%uses_normal_trace()) then
         call self%normal_trace(coupling, cavity, vn, error)
         if (allocated(error)) return
      end if
      same = .false.
      if (allocated(self%phi) .and. self%charges_valid) then
         if (size(self%phi) == size(phi)) then
            same = all(self%phi == phi)
            if (self%uses_normal_trace()) then
               if (allocated(self%vn)) then
                  if (size(self%vn) == size(vn)) then
                     same = same .and. all(self%vn == vn)
                  else
                     same = .false.
                  end if
               else
                  same = .false.
               end if
            end if
         end if
      end if
      if (same) return
      self%phi = phi
      if (self%uses_normal_trace()) self%vn = vn
      self%charges_valid = .false.
      call self%solve_layer(error)
      if (allocated(error)) return
      self%charges_valid = .true.
   end subroutine layer_ensure_charges

   !> Return dE/dphi from the generated adjoint solve
   subroutine layer_potential_adjoint(self, w_phi, error)
      class(pcm_layer_base), intent(in) :: self
      real(wp), allocatable, intent(out) :: w_phi(:)
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%charges_valid) then
         call fatal_error(error, "Layer PCM adjoint is unavailable")
         return
      end if
      allocate (w_phi, source=self%qsym)
   end subroutine layer_potential_adjoint

   !> Return dE/d(n dot grad phi), zero for potential-only models
   subroutine layer_normal_trace_adjoint(self, w_vn, error)
      class(pcm_layer_base), intent(in) :: self
      real(wp), allocatable, intent(out) :: w_vn(:)
      type(error_type), allocatable, intent(out) :: error
      if (.not. self%charges_valid) then
         call fatal_error(error, "Layer PCM normal adjoint is unavailable")
         return
      end if
      if (allocated(self%qn)) then
         allocate (w_vn, source=self%qn)
      else
         allocate (w_vn(size(self%q)), source=0.0_wp)
      end if
   end subroutine layer_normal_trace_adjoint

   !> Add the generated polarization energy
   subroutine layer_get_energy(self, coupling, cavity, energy, error)
      class(pcm_layer_base), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(inout) :: cavity
      real(wp), intent(inout) :: energy
      type(error_type), allocatable, intent(out) :: error
      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      energy = energy + self%kernel_energy
   end subroutine layer_get_energy

   !> Contract host potential derivatives for point or Gaussian surfaces
   subroutine layer_get_host_surface_weights(self, coupling, acc, ngrid, error)
      class(pcm_layer_base), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      integer, intent(in) :: ngrid
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: grad_phi(:, :), dphi_dxi(:), w_xyz(:, :), w_xi(:)
      integer :: i

      call coupling%read("potential", "dphi_dr", grad_phi, error)
      if (allocated(error)) return
      if (any(shape(grad_phi) /= [3, ngrid])) then
         call fatal_error(error, "Layer PCM host position derivative shape mismatch")
         return
      end if
      allocate (w_xyz(3, ngrid))
      do i = 1, ngrid
         w_xyz(:, i) = self%qsym(i)*grad_phi(:, i)
      end do
      if (self%smeared) then
         call coupling%read("potential", "dphi_dxi", dphi_dxi, error)
         if (allocated(error)) return
         if (size(dphi_dxi) /= ngrid) then
            call fatal_error(error, "Layer PCM host width derivative shape mismatch")
            return
         end if
         w_xi = self%qsym*dphi_dxi
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_xi=w_xi)
      else
         call acc%add_surface_weights(error, w_xyz=w_xyz)
      end if
   end subroutine layer_get_host_surface_weights

   !> Add generated geometry weights and response of host raw traces
   subroutine layer_get_surface_weights(self, coupling, cavity, acc, error)
      class(pcm_layer_base), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      call self%layer_weights(cavity, acc, error)
      if (allocated(error)) return
      if (cavity%has_field_dependent_geometry()) then
         call self%get_host_surface_weights(coupling, acc, cavity%ngrid, error)
         if (allocated(error)) return
         call layer_add_host_normal_weights(self, coupling, cavity, acc, error)
      end if
   end subroutine layer_get_surface_weights

   !> Add all geometry weights needed for the nuclear gradient
   subroutine layer_get_gradient_surface_weights(self, coupling, cavity, acc, error)
      class(pcm_layer_base), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      call self%get_surface_weights(coupling, cavity, acc, error)
      if (allocated(error)) return
      if (.not. cavity%has_field_dependent_geometry()) then
         call self%get_host_surface_weights(coupling, acc, cavity%ngrid, error)
         if (allocated(error)) return
         call layer_add_host_normal_weights(self, coupling, cavity, acc, error)
      end if
   end subroutine layer_get_gradient_surface_weights

   !> Forward gradient path equivalent to the reverse surface contraction
   subroutine layer_get_gradient(self, coupling, cavity, response, gradient, error)
      class(pcm_layer_base), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(inout) :: cavity
      type(response_type), intent(inout) :: response
      real(wp), intent(inout) :: gradient(:, :)
      type(error_type), allocatable, intent(out) :: error
      type(cavity_surface_adjoint_type) :: acc
      call self%get_trace_response(coupling, cavity, response, error)
      if (allocated(error)) return
      call self%get_direct_gradient(coupling, cavity, gradient, error)
      if (allocated(error)) return
      call acc%init(cavity%ngrid)
      call self%get_gradient_surface_weights(coupling, cavity, acc, error)
      if (allocated(error)) return
      call cavity%get_surface_gradient(acc, gradient, error)
   end subroutine layer_get_gradient

   !> Add derivatives of vn = n dot grad(phi) with respect to the surface
   subroutine layer_add_host_normal_weights(self, coupling, cavity, acc, error)
      class(pcm_layer_base), intent(in) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: grad_phi(:, :), hess_phi(:, :, :), mixed(:, :)
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_xi(:)
      integer :: i, n

      if (.not. allocated(self%qn)) return
      n = cavity%ngrid
      call coupling%read("potential", "dphi_dr", grad_phi, error)
      if (allocated(error)) return
      call coupling%read("potential", "d2phi_dr2", hess_phi, error)
      if (allocated(error)) return
      if (any(shape(grad_phi) /= [3, n]) .or. any(shape(hess_phi) /= [3, 3, n])) then
         call fatal_error(error, "Layer PCM host potential derivative shape mismatch")
         return
      end if
      allocate (w_xyz(3, n), w_n(3, n))
      do i = 1, n
         w_xyz(:, i) = self%qn(i)*matmul(transpose(hess_phi(:, :, i)), cavity%normal0(:, i))
         w_n(:, i) = self%qn(i)*grad_phi(:, i)
      end do
      if (self%smeared) then
         call coupling%read("potential", "d2phi_dr_dxi", mixed, error)
         if (allocated(error)) return
         if (any(shape(mixed) /= [3, n])) then
            call fatal_error(error, "Layer PCM host mixed derivative shape mismatch")
            return
         end if
         allocate (w_xi(n))
         do i = 1, n
            w_xi(i) = self%qn(i)*dot_product(cavity%normal0(:, i), mixed(:, i))
         end do
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_xi=w_xi)
      else
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n)
      end if
   end subroutine layer_add_host_normal_weights

end module moist_model_continuum_component_pcm_layer
