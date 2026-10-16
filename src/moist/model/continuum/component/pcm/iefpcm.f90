!> IEF-PCM component backed by generated Gaussian and point kernels
module moist_model_continuum_component_pcm_iefpcm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_layer, only: pcm_layer_base, pcm_layer_feps
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver
   use pcmgen_gauss_iefpcm, only: iefpcm_complex_state, &
      & iefpcm_assemble, iefpcm_medium, iefpcm_solve, iefpcm_weights, iefpcm_tangent
   use pcmgen_point_iefpcm, only: iefpcm_point_complex_state, &
      & iefpcm_point_assemble, iefpcm_point_medium, iefpcm_point_solve, &
      & iefpcm_point_weights, iefpcm_point_tangent
   implicit none(type, external)
   private

   public :: model_continuum_component_iefpcm, new_component_iefpcm

   !> Integral-equation-formalism PCM component
   !>
   !> - Physical input: solvent relative permittivity epsilon
   !> - Cavity positions, normals, principal curvatures and grid point areas
   !> - Gaussian charge widths and switching factors when available
   !> - Host electrostatic potential at the surface points
   !> - Host nuclear charges for gradients, with ECP core electrons removed
   type, extends(pcm_layer_base) :: model_continuum_component_iefpcm
      !> Gaussian surface kernel state
      type(iefpcm_complex_state) :: gauss
      !> Point surface kernel state
      type(iefpcm_point_complex_state) :: point
   contains
      procedure :: assemble_layer => iefpcm_assemble_layer
      procedure :: solve_layer => iefpcm_solve_layer
      procedure :: layer_weights => iefpcm_layer_weights
      procedure :: layer_tangent => iefpcm_layer_tangent
   end type model_continuum_component_iefpcm

contains

   !> Construct IEF-PCM from physical and numerical parameters
   !>
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self     initialized IEF-PCM component
   !> @param[in]  epsilon  dimensionless solvent relative permittivity, >= 1
   !> @param[out] error    allocated on invalid physical or numerical parameters
   !> @param[in]  param    optional solver controls, compiled defaults when omitted
   !> @param[in]  ctx      optional borrowed context, must outlive the component
   subroutine new_component_iefpcm(self, epsilon, error, param, ctx)
      !> IEF-PCM component
      type(model_continuum_component_iefpcm), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Solvent relative permittivity
      real(wp), intent(in) :: epsilon
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear solver and operator-storage controls
      type(moist_pcm_parameters_type), intent(in), optional :: param
      type(moist_pcm_parameters_type) :: settings

      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, "IEF-PCM epsilon must be >= 1")
         return
      end if
      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      if (present(ctx)) self%ctx => ctx
      self%epsilon = epsilon
      self%feps = pcm_layer_feps(epsilon)
      self%param = settings
      self%name = "IEF-PCM"
      self%description = "Integral-equation formulation of the polarizable continuum model"
   end subroutine new_component_iefpcm

   !> Assemble the IEF-PCM operators
   subroutine iefpcm_assemble_layer(self, cavity, solver, error)
      class(model_continuum_component_iefpcm), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(pcm_solver), intent(in) :: solver
      type(error_type), allocatable, intent(out) :: error
      integer :: ngrid

      ngrid = cavity%ngrid
      associate (r => cavity%xyz(:, 1:ngrid), n => cavity%normal0(:, 1:ngrid), &
         & k1 => cavity%k1(1:ngrid), k2 => cavity%k2(1:ngrid), &
         & a => cavity%a(1:ngrid))
         if (self%smeared) then
            associate (xi => cavity%xi0(1:ngrid), f => cavity%f(1:ngrid))
               call iefpcm_assemble(self%gauss%iefpcm_state, r, n, k1, k2, xi, f, a, error, solver)
               if (allocated(error)) return
               call iefpcm_medium(self%gauss%iefpcm_state, self%epsilon, error)
            end associate
         else
            call iefpcm_point_assemble(self%point%iefpcm_point_state, r, n, k1, k2, a, error, solver)
            if (allocated(error)) return
            call iefpcm_point_medium(self%point%iefpcm_point_state, self%epsilon, error)
         end if
      end associate
   end subroutine iefpcm_assemble_layer

   !> Solve the current IEF-PCM host traces
   subroutine iefpcm_solve_layer(self, error)
      class(model_continuum_component_iefpcm), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      if (self%smeared) then
         call iefpcm_solve(self%gauss%iefpcm_state, self%phi, self%kernel_energy, self%q, self%qsym, error)
      else
         call iefpcm_point_solve(self%point%iefpcm_point_state, self%phi, self%kernel_energy, self%q, self%qsym, error)
      end if
   end subroutine iefpcm_solve_layer

   !> Accumulate IEF-PCM operator surface weights
   subroutine iefpcm_layer_weights(self, cavity, acc, error)
      class(model_continuum_component_iefpcm), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp), allocatable :: w_a(:), w_xi(:), w_f(:)
      real(wp) :: w_eps
      integer :: ngrid

      ngrid = cavity%ngrid
      allocate (w_xyz(3, ngrid), w_n(3, ngrid), w_k1(ngrid), w_k2(ngrid), w_a(ngrid))
      if (self%smeared) then
         allocate (w_xi(ngrid), w_f(ngrid), source=0.0_wp)
         call iefpcm_weights(self%gauss%iefpcm_state, w_xyz, w_n, w_k1, w_k2, w_xi, w_f, w_a, w_eps, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a, w_xi=w_xi, w_f=w_f)
      else
         call iefpcm_point_weights(self%point%iefpcm_point_state, w_xyz, w_n, w_k1, w_k2, w_a, w_eps, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a)
      end if
   end subroutine iefpcm_layer_weights

   !> Directional derivative of the IEF-PCM weights and trace adjoints
   !>
   !> - no normal trace: `dvn` is ignored and `hvp_vn` stays zero
   !>
   !> @param[in,out] self   IEF-PCM component with current weights
   !> @param[in]     dr     Position direction (3, ngrid, 1)
   !> @param[in]     dn     Normal direction (3, ngrid, 1)
   !> @param[in]     dk1    First principal curvature direction (ngrid, 1)
   !> @param[in]     dk2    Second principal curvature direction (ngrid, 1)
   !> @param[in]     dxi    Gaussian width direction (ngrid, 1)
   !> @param[in]     df     Switching factor direction (ngrid, 1)
   !> @param[in]     da     Area direction (ngrid, 1)
   !> @param[in]     dv     Potential-trace direction (ngrid, 1)
   !> @param[in]     dvn    Normal-trace direction (ngrid, 1)
   !> @param[in,out] hvp_r  Position tangent (3, ngrid, 1)
   !> @param[in,out] hvp_n  Normal tangent (3, ngrid, 1)
   !> @param[in,out] hvp_k1 First principal curvature tangent (ngrid, 1)
   !> @param[in,out] hvp_k2 Second principal curvature tangent (ngrid, 1)
   !> @param[in,out] hvp_xi Gaussian width tangent (ngrid, 1)
   !> @param[in,out] hvp_f  Switching factor tangent (ngrid, 1)
   !> @param[in,out] hvp_a  Area tangent (ngrid, 1)
   !> @param[in,out] hvp_v  Potential adjoint tangent (ngrid, 1)
   !> @param[in,out] hvp_vn Normal-trace adjoint tangent (ngrid, 1)
   !> @param[out]    error  Allocated on failure
   subroutine iefpcm_layer_tangent(self, dr, dn, dk1, dk2, dxi, df, da, dv, dvn, hvp_r, hvp_n, &
         & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, hvp_vn, error)
      !> IEF-PCM component with current weights
      class(model_continuum_component_iefpcm), intent(inout) :: self
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

      if (self%smeared) then
         call iefpcm_tangent(self%gauss%iefpcm_state, dr, dn, dk1, dk2, dxi, df, da, dv, hvp_r, hvp_n, &
            & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, error)
      else
         call iefpcm_point_tangent(self%point%iefpcm_point_state, dr, dn, dk1, dk2, da, dv, hvp_r, hvp_n, &
            & hvp_k1, hvp_k2, hvp_a, hvp_v, error)
      end if
   end subroutine iefpcm_layer_tangent

end module moist_model_continuum_component_pcm_iefpcm
