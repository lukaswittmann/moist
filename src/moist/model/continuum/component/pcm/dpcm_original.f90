!> DPCM original component backed by generated Gaussian and point kernels
module moist_model_continuum_component_pcm_dpcm_original
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use mctc_env, only: wp, error_type, fatal_error
   use moist_channels_response, only: parameter_adjoint_response_type
   use moist_context, only: moist_context_type
   use moist_output_citations, only: citation
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_layer, only: pcm_layer_base, pcm_layer_feps
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver
   use pcmgen_gauss_dpcm_original, only: dpcm_original_complex_state, &
      & dpcm_original_assemble, dpcm_original_medium, dpcm_original_solve, dpcm_original_weights, dpcm_original_param_weights, dpcm_original_tangent
   use pcmgen_point_dpcm_original, only: dpcm_original_point_complex_state, &
      & dpcm_original_point_assemble, dpcm_original_point_medium, dpcm_original_point_solve, &
      & dpcm_original_point_weights, dpcm_original_point_param_weights, dpcm_original_point_tangent
   implicit none(type, external)
   private

   public :: model_continuum_component_dpcm_original, new_component_dpcm_original

   !> Original DPCM component
   !>
   !> - Physical input: solvent relative permittivity epsilon
   !> - Cavity positions, normals, principal curvatures and grid point areas
   !> - Gaussian charge widths when available
   !> - Host potential and its spatial derivative at the surface points
   !> - Host nuclear charges for gradients, with ECP core electrons removed
   type, extends(pcm_layer_base) :: model_continuum_component_dpcm_original
      !> Gaussian surface kernel state
      type(dpcm_original_complex_state) :: gauss
      !> Point surface kernel state
      type(dpcm_original_point_complex_state) :: point
   contains
      procedure :: assemble_layer => dpcm_original_assemble_layer
      procedure :: solve_layer => dpcm_original_solve_layer
      procedure :: layer_weights => dpcm_original_layer_weights
      procedure :: get_parameter_adjoint => dpcm_original_get_parameter_adjoint
      procedure :: layer_tangent => dpcm_original_layer_tangent
   end type model_continuum_component_dpcm_original

contains

   !> Construct DPCM original from physical and numerical parameters
   !>
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self     initialized original DPCM component
   !> @param[in]  epsilon  dimensionless solvent relative permittivity, >= 1
   !> @param[out] error    allocated on invalid physical or numerical parameters
   !> @param[in]  param    optional solver controls, compiled defaults when omitted
   !> @param[in]  ctx      optional borrowed context, must outlive the component
   subroutine new_component_dpcm_original(self, epsilon, error, param, ctx)
      !> Original DPCM component
      type(model_continuum_component_dpcm_original), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Solvent relative permittivity
      real(wp), intent(in) :: epsilon
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear solver and operator-storage controls
      type(moist_pcm_parameters_type), intent(in), optional :: param
      type(moist_pcm_parameters_type) :: settings

      if (epsilon < 1.0_wp .or. ieee_is_nan(epsilon)) then
         call fatal_error(error, "DPCM original epsilon must be >= 1")
         return
      end if
      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      if (present(ctx)) self%ctx => ctx
      self%epsilon = epsilon
      self%feps = pcm_layer_feps(epsilon)
      self%param = settings
      self%name = "DPCM original"
      self%description = "Dielectric polarizable continuum in the published Miertus-Scrocco-Tomasi form"
      self%needs_switching = .false.
      self%needs_normal = .true.
   end subroutine new_component_dpcm_original

   !> Assemble the DPCM original operators
   subroutine dpcm_original_assemble_layer(self, cavity, solver, error)
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(pcm_solver), intent(in) :: solver
      type(error_type), allocatable, intent(out) :: error
      integer :: ngrid

      call self%ctx%cite(citation%dpcm_original)
      ngrid = cavity%ngrid
      associate (r => cavity%xyz(:, 1:ngrid), n => cavity%normal0(:, 1:ngrid), &
         & k1 => cavity%k1(1:ngrid), k2 => cavity%k2(1:ngrid), &
         & a => cavity%a(1:ngrid))
         if (self%smeared) then
            associate (xi => cavity%xi0(1:ngrid))
               call dpcm_original_assemble(self%gauss%dpcm_original_state, r, n, k1, k2, xi, a, error, solver)
               if (allocated(error)) return
               call dpcm_original_medium(self%gauss%dpcm_original_state, self%epsilon, error)
            end associate
         else
            call dpcm_original_point_assemble(self%point%dpcm_original_point_state, r, n, k1, k2, a, error, solver)
            if (allocated(error)) return
            call dpcm_original_point_medium(self%point%dpcm_original_point_state, self%epsilon, error)
         end if
      end associate
   end subroutine dpcm_original_assemble_layer

   !> Solve the current DPCM original host traces
   subroutine dpcm_original_solve_layer(self, error)
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      if (self%smeared) then
         call dpcm_original_solve(self%gauss%dpcm_original_state, self%phi, self%vn, &
            & self%kernel_energy, self%q, self%qsym, self%qn, error)
      else
         call dpcm_original_point_solve(self%point%dpcm_original_point_state, self%phi, self%vn, &
            & self%kernel_energy, self%q, self%qsym, self%qn, error)
      end if
   end subroutine dpcm_original_solve_layer

   !> Accumulate DPCM original operator surface weights
   subroutine dpcm_original_layer_weights(self, cavity, acc, error)
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp), allocatable :: w_a(:), w_xi(:), w_f(:)
      integer :: ngrid

      ngrid = cavity%ngrid
      allocate (w_xyz(3, ngrid), w_n(3, ngrid), w_k1(ngrid), w_k2(ngrid), w_a(ngrid))
      if (self%smeared) then
         allocate (w_xi(ngrid), source=0.0_wp)
         call dpcm_original_weights(self%gauss%dpcm_original_state, w_xyz, w_n, w_k1, w_k2, w_xi, w_a, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a, w_xi=w_xi)
      else
         call dpcm_original_point_weights(self%point%dpcm_original_point_state, w_xyz, w_n, w_k1, w_k2, w_a, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a)
      end if
   end subroutine dpcm_original_layer_weights

   !> Directional derivative of the original DPCM weights and trace adjoints
   !>
   !> - no switching factors: `df` is ignored and `hvp_f` stays zero
   !>
   !> @param[in,out] self   original DPCM component with current weights
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
   subroutine dpcm_original_layer_tangent(self, dr, dn, dk1, dk2, dxi, df, da, dv, dvn, hvp_r, hvp_n, &
         & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, hvp_vn, error)
      !> original DPCM component with current weights
      class(model_continuum_component_dpcm_original), intent(inout) :: self
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
         call dpcm_original_tangent(self%gauss%dpcm_original_state, dr, dn, dk1, dk2, dxi, da, dv, dvn, hvp_r, &
            & hvp_n, hvp_k1, hvp_k2, hvp_xi, hvp_a, hvp_v, hvp_vn, error)
      else
         call dpcm_original_point_tangent(self%point%dpcm_original_point_state, dr, dn, dk1, dk2, da, dv, dvn, &
            & hvp_r, hvp_n, hvp_k1, hvp_k2, hvp_a, hvp_v, hvp_vn, error)
      end if
   end subroutine dpcm_original_layer_tangent

   !> Derivatives of the energy with respect to the physical parameters
   !>
   !> @param[in,out] self  dpcm_original component with current charges
   !> @param[out]    item  parameter derivatives
   !> @param[out]    error allocated on failure
   subroutine dpcm_original_get_parameter_adjoint(self, item, error)
      !> Component with current charges
      class(model_continuum_component_dpcm_original), intent(inout) :: self
      !> Parameter derivatives
      type(parameter_adjoint_response_type), intent(out) :: item
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: w_eps

      if (self%smeared) then
         call dpcm_original_param_weights(self%gauss%dpcm_original_state, w_eps, error)
      else
         call dpcm_original_point_param_weights(self%point%dpcm_original_point_state, w_eps, error)
      end if
      if (allocated(error)) return
      allocate (item%dg_deps, source=0.0_wp)
      if (ieee_is_finite(self%epsilon)) item%dg_deps = w_eps

   end subroutine dpcm_original_get_parameter_adjoint

end module moist_model_continuum_component_pcm_dpcm_original
