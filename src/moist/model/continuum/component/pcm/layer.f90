!> Shared host integration for single- and double-layer PCM models
module moist_model_continuum_component_pcm_layer
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & coupling_request_type, point_potential_request_type, gaussian_potential_request_type, &
      & nuclear_charge_request_type, moist_phase_energy, moist_phase_response, moist_phase_gradient, &
      & coupling_register, request_require
   use moist_channels_response, only: response_type
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, solver_type, &
      & pcm_tangent_block, pcm_tangent_require
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver, pcm_dense_solver, &
      & pcm_gmres_solver, pcm_block_jacobi
   implicit none(type, external)
   private

   public :: pcm_layer_base, pcm_layer_feps, layer_configure_cavity

   !> Common host lifecycle; each extension owns one generated model state
   type, abstract, extends(model_continuum_component_pcm) :: pcm_layer_base
      !> Symmetrized potential adjoint
      real(wp), allocatable :: qsym(:)
      !> Normal-derivative adjoint, for DPCM and SLIC
      real(wp), allocatable :: qn(:)
      !> Host normal-derivative trace
      real(wp), allocatable :: vn(:)
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
      procedure :: get_surface_weights => layer_get_surface_weights
      procedure :: get_gradient_surface_weights => layer_get_gradient_surface_weights
      procedure :: get_gradient => layer_get_gradient
      procedure :: declare_coupling => layer_declare_coupling
      procedure :: get_host_surface_weights => layer_get_host_surface_weights
      procedure :: surface_weight_tangent => layer_surface_weight_tangent
      procedure(assemble_layer_iface), deferred :: assemble_layer
      procedure(solve_layer_iface), deferred :: solve_layer
      procedure(layer_weights_iface), deferred :: layer_weights
      procedure(layer_tangent_iface), deferred :: layer_tangent
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

      !> Directional derivative of the generated weights and trace adjoints
      !>
      !> - Called after `layer_weights` on the current solution
      !> - One direction block (nblk = 1); channels the kernel lacks are ignored
      !>   on input and stay zero on output
      !>
      !> @param[in,out] self   Layer component with current charges
      !> @param[in]     dr     Position direction (3, ngrid, 1)
      !> @param[in]     dn     Normal direction (3, ngrid, 1)
      !> @param[in]     dk1    First principal curvature direction (ngrid, 1)
      !> @param[in]     dk2    Second principal curvature direction (ngrid, 1)
      !> @param[in]     dxi    Gaussian width direction (ngrid, 1)
      !> @param[in]     df     Switching factor direction (ngrid, 1)
      !> @param[in]     da     Area direction (ngrid, 1)
      !> @param[in]     dv     Potential-trace direction (ngrid, 1)
      !> @param[in]     dvn    Normal-trace direction (ngrid, 1)
      !> @param[in,out] hvp_r  Position tangent (3, ngrid, 1), zero on entry
      !> @param[in,out] hvp_n  Normal tangent (3, ngrid, 1), zero on entry
      !> @param[in,out] hvp_k1 First principal curvature tangent (ngrid, 1), zero on entry
      !> @param[in,out] hvp_k2 Second principal curvature tangent (ngrid, 1), zero on entry
      !> @param[in,out] hvp_xi Gaussian width tangent (ngrid, 1), zero on entry
      !> @param[in,out] hvp_f  Switching factor tangent (ngrid, 1), zero on entry
      !> @param[in,out] hvp_a  Area tangent (ngrid, 1), zero on entry
      !> @param[in,out] hvp_v  Potential adjoint tangent (ngrid, 1), zero on entry
      !> @param[in,out] hvp_vn Normal-trace adjoint tangent (ngrid, 1), zero on entry
      !> @param[out]    error  Allocated on failure
      subroutine layer_tangent_iface(self, dr, dn, dk1, dk2, dxi, df, da, dv, dvn, hvp_r, hvp_n, &
            & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, hvp_vn, error)
         import :: pcm_layer_base, error_type, wp
         implicit none(type, external)
         !> Layer component with current charges
         class(pcm_layer_base), intent(inout) :: self
         !> Position direction
         real(wp), contiguous, intent(in) :: dr(:, :, :)
         !> Normal direction
         real(wp), contiguous, intent(in) :: dn(:, :, :)
         !> First principal curvature direction
         real(wp), contiguous, intent(in) :: dk1(:, :)
         !> Second principal curvature direction
         real(wp), contiguous, intent(in) :: dk2(:, :)
         !> Gaussian width direction
         real(wp), contiguous, intent(in) :: dxi(:, :)
         !> Switching factor direction
         real(wp), contiguous, intent(in) :: df(:, :)
         !> Area direction
         real(wp), contiguous, intent(in) :: da(:, :)
         !> Potential-trace direction
         real(wp), contiguous, intent(in) :: dv(:, :)
         !> Normal-trace direction
         real(wp), contiguous, intent(in) :: dvn(:, :)
         !> Position tangent
         real(wp), contiguous, intent(inout) :: hvp_r(:, :, :)
         !> Normal tangent
         real(wp), contiguous, intent(inout) :: hvp_n(:, :, :)
         !> First principal curvature tangent
         real(wp), contiguous, intent(inout) :: hvp_k1(:, :)
         !> Second principal curvature tangent
         real(wp), contiguous, intent(inout) :: hvp_k2(:, :)
         !> Gaussian width tangent
         real(wp), contiguous, intent(inout) :: hvp_xi(:, :)
         !> Switching factor tangent
         real(wp), contiguous, intent(inout) :: hvp_f(:, :)
         !> Area tangent
         real(wp), contiguous, intent(inout) :: hvp_a(:, :)
         !> Potential adjoint tangent
         real(wp), contiguous, intent(inout) :: hvp_v(:, :)
         !> Normal-trace adjoint tangent
         real(wp), contiguous, intent(inout) :: hvp_vn(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine layer_tangent_iface
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
      type(nuclear_charge_request_type) :: charges
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
      if (allocated(error)) return
      call request_require(charges, moist_phase_gradient, "z", error)
      if (allocated(error)) return
      call coupling_register(coupling, "nuclear_charges", charges, error)
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
            ! Typed allocate; ifx seems to have issues here sometimes
            allocate (pcm_block_jacobi :: solver%preconditioner)
            select type (preconditioner => solver%preconditioner)
            type is (pcm_block_jacobi)
               preconditioner%block_size = 256
            end select
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

   !> Directional derivative of the generated weights and trace adjoints
   !>
   !> - Surface channels the concrete kernel reads; xi and f on Gaussian surfaces only
   !> - `dw_vn` is allocated for normal-trace models when requested
   !>
   !> @param[in,out] self      Layer component with current charges
   !> @param[in]     cavity    Cavity the operators were assembled on
   !> @param[in]     direction Surface direction
   !> @param[in]     dphi      Potential-trace direction (ngrid)
   !> @param[in,out] tangent   Initialized accumulator, receives the weight tangent
   !> @param[out]    dw_phi    Tangent of the symmetrized potential adjoint (ngrid)
   !> @param[out]    error     Allocated on failure
   !> @param[in]     dvn       Normal-trace direction (ngrid); absent is zero
   !> @param[out]    dw_vn     Tangent of the normal-trace adjoint (ngrid)
   subroutine layer_surface_weight_tangent(self, cavity, direction, dphi, tangent, dw_phi, &
         & error, dvn, dw_vn)
      !> Layer component with current charges
      class(pcm_layer_base), intent(inout) :: self
      !> Cavity the operators were assembled on
      class(cavity_type), intent(in) :: cavity
      !> Surface direction
      class(cavity_surface_adjoint_type), intent(in) :: direction
      !> Potential-trace direction
      real(wp), intent(in) :: dphi(:)
      !> Accumulator of the weight tangent
      class(cavity_surface_adjoint_type), intent(inout) :: tangent
      !> Tangent of the potential adjoint
      real(wp), allocatable, intent(out) :: dw_phi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Normal-trace direction
      real(wp), intent(in), optional :: dvn(:)
      !> Tangent of the normal-trace adjoint
      real(wp), allocatable, intent(out), optional :: dw_vn(:)
      !> Kernel direction blocks
      real(wp), allocatable :: dr(:, :, :), dn(:, :, :), dk1(:, :), dk2(:, :), dxi(:, :), &
         & df(:, :), da(:, :), dv(:, :), dvn_block(:, :)
      !> Kernel tangent blocks
      real(wp), allocatable :: hvp_r(:, :, :), hvp_n(:, :, :), hvp_k1(:, :), hvp_k2(:, :), &
         & hvp_xi(:, :), hvp_f(:, :), hvp_a(:, :), hvp_v(:, :), hvp_vn(:, :)
      !> Kernel weights, recomputed because the tangent requires them
      type(cavity_surface_adjoint_type) :: scratch
      !> Number of grid points
      integer :: ngrid

      call pcm_tangent_require(self, cavity, dphi, error, dvn)
      if (allocated(error)) return
      ngrid = cavity%ngrid
      call pcm_tangent_block(direction%w_xyz, ngrid, "w_xyz", dr, error)
      if (allocated(error)) return
      call pcm_tangent_block(direction%w_n, ngrid, "w_n", dn, error)
      if (allocated(error)) return
      call pcm_tangent_block(direction%w_k1, ngrid, "w_k1", dk1, error)
      if (allocated(error)) return
      call pcm_tangent_block(direction%w_k2, ngrid, "w_k2", dk2, error)
      if (allocated(error)) return
      call pcm_tangent_block(direction%w_xi, ngrid, "w_xi", dxi, error)
      if (allocated(error)) return
      call pcm_tangent_block(direction%w_f, ngrid, "w_f", df, error)
      if (allocated(error)) return
      call pcm_tangent_block(direction%w_a, ngrid, "w_a", da, error)
      if (allocated(error)) return
      allocate (dv(ngrid, 1), dvn_block(ngrid, 1), source=0.0_wp)
      dv(:, 1) = dphi
      if (present(dvn)) dvn_block(:, 1) = dvn
      allocate (hvp_r(3, ngrid, 1), hvp_n(3, ngrid, 1), source=0.0_wp)
      allocate (hvp_k1(ngrid, 1), hvp_k2(ngrid, 1), hvp_xi(ngrid, 1), hvp_f(ngrid, 1), &
         & hvp_a(ngrid, 1), hvp_v(ngrid, 1), hvp_vn(ngrid, 1), source=0.0_wp)

      call scratch%init(ngrid)
      call self%layer_weights(cavity, scratch, error)
      if (allocated(error)) return
      call self%layer_tangent(dr, dn, dk1, dk2, dxi, df, da, dv, dvn_block, hvp_r, hvp_n, &
         & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, hvp_vn, error)
      if (allocated(error)) return
      if (self%smeared) then
         call tangent%add_surface_weights(error, w_xyz=hvp_r(:, :, 1), w_n=hvp_n(:, :, 1), &
            & w_k1=hvp_k1(:, 1), w_k2=hvp_k2(:, 1), w_xi=hvp_xi(:, 1), w_f=hvp_f(:, 1), &
            & w_a=hvp_a(:, 1))
      else
         call tangent%add_surface_weights(error, w_xyz=hvp_r(:, :, 1), w_n=hvp_n(:, :, 1), &
            & w_k1=hvp_k1(:, 1), w_k2=hvp_k2(:, 1), w_a=hvp_a(:, 1))
      end if
      if (allocated(error)) return
      dw_phi = hvp_v(:, 1)
      if (present(dw_vn) .and. self%uses_normal_trace()) dw_vn = hvp_vn(:, 1)
   end subroutine layer_surface_weight_tangent

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
