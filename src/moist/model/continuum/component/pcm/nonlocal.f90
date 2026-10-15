!> Nonlocal electrostatics component backed by generated boundary operators
module moist_model_continuum_component_pcm_nonlocal
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_utils_prettyprint, only: prettyprinter
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_layer, only: pcm_layer_base, pcm_layer_feps
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver
   use pcmgen_gauss_nonlocal, only: nonlocal_pcm_state, &
      & nonlocal_pcm_assemble, nonlocal_pcm_medium, nonlocal_pcm_solve, nonlocal_pcm_weights, nonlocal_pcm_tangent
   use pcmgen_point_nonlocal, only: nonlocal_pcm_point_state, &
      & nonlocal_pcm_point_assemble, nonlocal_pcm_point_medium, nonlocal_pcm_point_solve, &
      & nonlocal_pcm_point_weights, nonlocal_pcm_point_tangent
   implicit none(type, external)
   private

   public :: model_continuum_component_nonlocal, new_component_nonlocal

   !> Three-field nonlocal boundary system
   !>
   !> - Physical inputs: eps_solvent, eps_inf, lambda and optional eps_solute
   !> - Cavity positions, normals and grid point areas
   !> - Gaussian charge widths and switching factors when available
   !> - Host potential and its outward normal derivative at the surface points
   !> - Host nuclear charges for gradients, with ECP core electrons removed
   type, extends(pcm_layer_base) :: model_continuum_component_nonlocal
      !> High-frequency solvent relative permittivity, dimensionless
      real(wp) :: eps_inf = 1.0_wp
      !> Relative permittivity inside the cavity, dimensionless
      real(wp) :: eps_solute = 1.0_wp
      !> Solvent polarization correlation length, bohr
      real(wp) :: lambda = 1.0_wp
      !> Interior potential trace times 4*pi, shape (ngrid)
      real(wp), allocatable :: u(:)
      !> Exterior auxiliary trace times 4*pi, shape (ngrid)
      real(wp), allocatable :: w_ext(:)
      !> Gaussian surface kernel state
      type(nonlocal_pcm_state) :: gauss
      !> Point surface kernel state
      type(nonlocal_pcm_point_state) :: point
   contains
      procedure :: print_inputs => nonlocal_print_inputs
      procedure :: assemble_layer => nonlocal_assemble_layer
      procedure :: solve_layer => nonlocal_solve_layer
      procedure :: layer_weights => nonlocal_layer_weights
      procedure :: layer_tangent => nonlocal_layer_tangent
      procedure :: is_inert => nonlocal_is_inert
   end type model_continuum_component_nonlocal

contains

   !> Construct nonlocal electrostatics from permittivities and correlation length
   !>
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self         initialized nonlocal PCM component
   !> @param[in]  eps_solvent  dimensionless static solvent relative permittivity, finite and >= 1
   !> @param[in]  eps_inf      dimensionless high-frequency solvent relative permittivity, finite and >= 1
   !> @param[in]  lambda       solvent polarization correlation length in bohr, finite and > 0
   !> @param[out] error        allocated on invalid physical or numerical parameters
   !> @param[in]  param        optional solver controls, compiled defaults when omitted
   !> @param[in]  eps_solute   dimensionless relative permittivity inside the cavity, finite and >= 1; default 1
   !> @param[in]  ctx          optional borrowed context, must outlive the component
   subroutine new_component_nonlocal(self, eps_solvent, eps_inf, lambda, error, param, eps_solute, ctx)
      !> Nonlocal PCM component
      type(model_continuum_component_nonlocal), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Static solvent relative permittivity
      real(wp), intent(in) :: eps_solvent
      !> High-frequency solvent relative permittivity
      real(wp), intent(in) :: eps_inf
      !> Solvent polarization correlation length
      real(wp), intent(in) :: lambda
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear solver and operator-storage controls
      type(moist_pcm_parameters_type), intent(in), optional :: param
      !> Relative permittivity inside the cavity
      real(wp), intent(in), optional :: eps_solute
      type(moist_pcm_parameters_type) :: settings

      if (.not. ieee_is_finite(eps_solvent) .or. eps_solvent < 1.0_wp .or. &
          & .not. ieee_is_finite(eps_inf) .or. eps_inf < 1.0_wp .or. &
          & .not. ieee_is_finite(lambda) .or. lambda <= 0.0_wp) then
         call fatal_error(error, "Nonlocal PCM requires permittivities >= 1 and positive lambda")
         return
      end if
      if (present(eps_solute)) then
         if (.not. ieee_is_finite(eps_solute) .or. eps_solute < 1.0_wp) then
            call fatal_error(error, "Nonlocal PCM solute permittivity must be >= 1")
            return
         end if
         self%eps_solute = eps_solute
      end if
      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      if (present(ctx)) self%ctx => ctx
      self%epsilon = eps_solvent
      self%feps = pcm_layer_feps(eps_solvent)
      self%eps_inf = eps_inf
      self%lambda = lambda
      self%param = settings
      self%needs_normal = .true.
      self%needs_curvature = .false.
      self%name = "Nonlocal electrostatics"
      self%description = "Three-field nonlocal dielectric model with solvent polarization correlations"
   end subroutine new_component_nonlocal

   !> Print the physical inputs
   !>
   !> @param[in]    self Component instance
   !> @param[inout] pp   Pretty printer inside the component section
   subroutine nonlocal_print_inputs(self, pp)
      !> Component instance
      class(model_continuum_component_nonlocal), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp

      call pp%kv("Epsilon", self%epsilon)
      call pp%kv("Epsilon (high-frequency)", self%eps_inf)
      call pp%kv("Epsilon (solute)", self%eps_solute)
      call pp%kv("Correlation length", self%lambda, "bohr")

   end subroutine nonlocal_print_inputs

   !> Assemble Gaussian or point boundary operators
   subroutine nonlocal_assemble_layer(self, cavity, solver, error)
      class(model_continuum_component_nonlocal), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(pcm_solver), intent(in) :: solver
      type(error_type), allocatable, intent(out) :: error
      integer :: ngrid

      ngrid = cavity%ngrid
      associate (r => cavity%xyz(:, 1:ngrid), n => cavity%normal0(:, 1:ngrid), &
         & a => cavity%a(1:ngrid))
         if (self%smeared) then
            associate (xi => cavity%xi0(1:ngrid), f => cavity%f(1:ngrid))
               call nonlocal_pcm_assemble(self%gauss, r, n, xi, f, a, error, solver)
            end associate
            if (allocated(error)) return
            call nonlocal_pcm_medium(self%gauss, self%epsilon, self%eps_inf, &
               & self%lambda, self%eps_solute, error)
         else
            call nonlocal_pcm_point_assemble(self%point, r, n, a, error, solver)
            if (allocated(error)) return
            call nonlocal_pcm_point_medium(self%point, self%epsilon, self%eps_inf, &
               & self%lambda, self%eps_solute, error)
         end if
      end associate
   end subroutine nonlocal_assemble_layer

   !> Solve the three fields and retain both host-trace adjoints
   subroutine nonlocal_solve_layer(self, error)
      class(model_continuum_component_nonlocal), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      integer :: ngrid

      ngrid = size(self%phi)
      if (allocated(self%u)) deallocate (self%u)
      if (allocated(self%w_ext)) deallocate (self%w_ext)
      allocate (self%u(ngrid), self%w_ext(ngrid))
      if (self%smeared) then
         call nonlocal_pcm_solve(self%gauss, self%phi, self%vn, self%kernel_energy, &
            & self%u, self%q, self%w_ext, self%qsym, self%qn, error)
      else
         call nonlocal_pcm_point_solve(self%point, self%phi, self%vn, self%kernel_energy, &
            & self%u, self%q, self%w_ext, self%qsym, self%qn, error)
      end if
   end subroutine nonlocal_solve_layer

   !> Accumulate the operator response at fixed host traces
   subroutine nonlocal_layer_weights(self, cavity, acc, error)
      class(model_continuum_component_nonlocal), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_xi(:), w_f(:), w_a(:)
      real(wp) :: w_eps, w_eps_inf, w_lambda, w_eps_solute
      integer :: ngrid

      ngrid = cavity%ngrid
      allocate (w_xyz(3, ngrid), w_n(3, ngrid), w_a(ngrid))
      if (self%smeared) then
         allocate (w_xi(ngrid), w_f(ngrid))
         call nonlocal_pcm_weights(self%gauss, w_xyz, w_n, w_xi, w_f, w_a, w_eps, &
            & w_eps_inf, w_lambda, w_eps_solute, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_xi=w_xi, w_f=w_f, w_a=w_a)
      else
         call nonlocal_pcm_point_weights(self%point, w_xyz, w_n, w_a, w_eps, &
            & w_eps_inf, w_lambda, w_eps_solute, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_a=w_a)
      end if
   end subroutine nonlocal_layer_weights

   !> Directional derivative of the nonlocal weights and trace adjoints
   !>
   !> - no curvature: `dk1` and `dk2` are ignored
   !> - point surface: no widths or switching factors, `dxi` and `df` are ignored
   !>
   !> @param[in,out] self   nonlocal component with current weights
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
   subroutine nonlocal_layer_tangent(self, dr, dn, dk1, dk2, dxi, df, da, dv, dvn, hvp_r, hvp_n, &
         & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, hvp_vn, error)
      !> nonlocal component with current weights
      class(model_continuum_component_nonlocal), intent(inout) :: self
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
         call nonlocal_pcm_tangent(self%gauss, dr, dn, dxi, df, da, dv, dvn, hvp_r, hvp_n, hvp_xi, &
            & hvp_f, hvp_a, hvp_v, hvp_vn, error)
      else
         call nonlocal_pcm_point_tangent(self%point, dr, dn, da, dv, dvn, hvp_r, hvp_n, hvp_a, &
            & hvp_v, hvp_vn, error)
      end if
   end subroutine nonlocal_layer_tangent

   !> Vacuum bulk remains active when the other permittivities differ
   pure function nonlocal_is_inert(self) result(inert)
      class(model_continuum_component_nonlocal), intent(in) :: self
      logical :: inert
      inert = self%epsilon == 1.0_wp .and. self%eps_inf == 1.0_wp .and. self%eps_solute == 1.0_wp
   end function nonlocal_is_inert

end module moist_model_continuum_component_pcm_nonlocal
