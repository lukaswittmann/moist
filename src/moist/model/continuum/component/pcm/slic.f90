!> SLIC component backed by generated Gaussian and point kernels
module moist_model_continuum_component_pcm_slic
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_layer, only: pcm_layer_base, pcm_layer_feps
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver
   use pcmgen_gauss_slic, only: slic_state, &
      & slic_assemble, slic_medium, slic_solve, slic_weights
   use pcmgen_point_slic, only: slic_point_state, &
      & slic_point_assemble, slic_point_medium, slic_point_solve, &
      & slic_point_weights
   implicit none(type, external)
   private

   public :: model_continuum_component_slic, new_component_slic

   !> Solvation-layer interface condition component
   !>
   !> - Physical inputs: epsilon, alpha, beta, gamma; optional mu and phi_static
   !> - Perturbation `h(E_n) = alpha*tanh(beta*E_n - gamma) + mu`
   !> - `E_n`: total outward normal electric field in atomic units
   !> - Static energy contribution `phi_static*Q`, with net solute charge `Q` in elementary charges
   !> - Cavity positions, normals, principal curvatures and grid point areas
   !> - Gaussian charge widths when available
   !> - Host potential and its outward normal derivative at the surface points
   type, extends(pcm_layer_base) :: model_continuum_component_slic
      !> Gaussian surface kernel state
      type(slic_state) :: gauss
      !> Point surface kernel state
      type(slic_point_state) :: point
      !> Amplitude of the solvation-layer perturbation h, dimensionless
      real(wp) :: alpha = 0.0_wp
      !> Transition steepness, inverse atomic electric-field units
      real(wp) :: beta = 0.0_wp
      !> Offset of the tanh argument `beta*E_n - gamma`, dimensionless
      real(wp) :: gamma = 0.0_wp
      !> Additive shift of the perturbation h, dimensionless
      real(wp) :: mu = 0.0_wp
      !> Static potential offset, hartree per elementary charge
      real(wp) :: phi_static = 0.0_wp
   contains
      procedure :: assemble_layer => slic_assemble_layer
      procedure :: solve_layer => slic_solve_layer
      procedure :: layer_weights => slic_layer_weights
   end type model_continuum_component_slic

contains

   !> Construct SLIC from physical and numerical parameters
   !>
   !> - Flux constraint: `2*feps*(mu + abs(alpha)) < 1`
   !> - Dielectric factor `feps = (epsilon - 1)/(epsilon + 1)`, conductor limit 1
   !> - Numerical constraint: `param%newton_tol < 1`
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self        initialized SLIC component
   !> @param[in]  epsilon     dimensionless solvent relative permittivity, >= 1
   !> @param[in]  alpha       dimensionless perturbation amplitude, finite
   !> @param[in]  beta        inverse atomic electric-field units, finite
   !> @param[in]  gamma       dimensionless transition offset, finite
   !> @param[out] error       allocated on invalid physical or numerical parameters
   !> @param[in]  param       optional linear/nonlinear solver controls, compiled defaults when omitted
   !> @param[in]  mu          dimensionless additive shift, finite; default `alpha*tanh(gamma)`
   !> @param[in]  phi_static  static potential in hartree per elementary charge, finite; default 0
   !> @param[in]  ctx         optional borrowed context, must outlive the component
   subroutine new_component_slic(self, epsilon, alpha, beta, gamma, error, param, mu, phi_static, ctx)
      !> SLIC component
      type(model_continuum_component_slic), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Solvent relative permittivity
      real(wp), intent(in) :: epsilon
      !> Amplitude of the solvation-layer perturbation h
      real(wp), intent(in) :: alpha
      !> Transition steepness in the normal electric field
      real(wp), intent(in) :: beta
      !> Offset of the tanh argument `beta*E_n - gamma`
      real(wp), intent(in) :: gamma
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear and nonlinear solver controls
      type(moist_pcm_parameters_type), intent(in), optional :: param
      !> Additive shift of the perturbation h
      real(wp), intent(in), optional :: mu
      !> Static potential offset
      real(wp), intent(in), optional :: phi_static
      type(moist_pcm_parameters_type) :: settings

      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, "SLIC epsilon must be >= 1")
         return
      end if
      if (.not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(beta) .or. &
          & .not. ieee_is_finite(gamma)) then
         call fatal_error(error, "SLIC parameters must be finite")
         return
      end if
      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      if (present(ctx)) self%ctx => ctx
      self%epsilon = epsilon
      self%feps = pcm_layer_feps(epsilon)
      self%param = settings
      self%name = "SLIC"
      self%needs_switching = .false.
      self%needs_normal = .true.
      self%alpha = alpha
      self%beta = beta
      self%gamma = gamma
      self%mu = alpha*tanh(gamma)
      if (present(mu)) self%mu = mu
      if (present(phi_static)) self%phi_static = phi_static
      if (.not. ieee_is_finite(self%mu) .or. .not. ieee_is_finite(self%phi_static)) then
         call fatal_error(error, "SLIC mu and static potential must be finite")
         return
      end if
      if (.not. 2.0_wp*self%feps*(self%mu + abs(alpha)) < 1.0_wp) then
         call fatal_error(error, "SLIC flux condition is singular for these parameters")
         return
      end if
      if (self%param%newton_tol >= 1.0_wp) then
         call fatal_error(error, "SLIC Newton tolerance must be below one")
         return
      end if
   end subroutine new_component_slic

   !> Assemble the SLIC operators
   subroutine slic_assemble_layer(self, cavity, solver, error)
      class(model_continuum_component_slic), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(pcm_solver), intent(in) :: solver
      type(error_type), allocatable, intent(out) :: error
      integer :: ngrid

      ngrid = cavity%ngrid
      associate (r => cavity%xyz(:, 1:ngrid), n => cavity%normal0(:, 1:ngrid), &
         & k1 => cavity%k1(1:ngrid), k2 => cavity%k2(1:ngrid), &
         & a => cavity%a(1:ngrid))
         if (self%smeared) then
            associate (xi => cavity%xi0(1:ngrid))
               call slic_assemble(self%gauss, r, n, k1, k2, xi, a, error, solver)
               if (allocated(error)) return
               call slic_medium(self%gauss, self%epsilon, self%alpha, self%beta, self%gamma, self%mu, error)
            end associate
         else
            call slic_point_assemble(self%point, r, n, k1, k2, a, error, solver)
            if (allocated(error)) return
            call slic_point_medium(self%point, self%epsilon, self%alpha, self%beta, self%gamma, self%mu, error)
         end if
      end associate
   end subroutine slic_assemble_layer

   !> Solve the current SLIC host traces
   subroutine slic_solve_layer(self, error)
      class(model_continuum_component_slic), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      if (self%smeared) then
         call slic_solve(self%gauss, self%phi, self%vn, self%kernel_energy, self%q, self%qsym, self%qn, error)
      else
         call slic_point_solve(self%point, self%phi, self%vn, self%kernel_energy, self%q, self%qsym, self%qn, error)
      end if
      if (.not. allocated(error)) self%kernel_energy = self%kernel_energy + self%phi_static*self%mol_solu%charge
   end subroutine slic_solve_layer

   !> Accumulate SLIC operator surface weights
   subroutine slic_layer_weights(self, cavity, acc, error)
      class(model_continuum_component_slic), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp), allocatable :: w_a(:), w_xi(:), w_f(:)
      real(wp) :: w_eps
      real(wp) :: w_alpha, w_beta, w_gamma, w_mu
      integer :: ngrid

      ngrid = cavity%ngrid
      allocate (w_xyz(3, ngrid), w_n(3, ngrid), w_k1(ngrid), w_k2(ngrid), w_a(ngrid))
      if (self%smeared) then
         allocate (w_xi(ngrid), source=0.0_wp)
         call slic_weights(self%gauss, w_xyz, w_n, w_k1, w_k2, w_xi, w_a, w_eps, w_alpha, w_beta, w_gamma, w_mu, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a, w_xi=w_xi)
      else
         call slic_point_weights(self%point, w_xyz, w_n, w_k1, w_k2, w_a, w_eps, w_alpha, w_beta, w_gamma, w_mu, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a)
      end if
   end subroutine slic_layer_weights

end module moist_model_continuum_component_pcm_slic
