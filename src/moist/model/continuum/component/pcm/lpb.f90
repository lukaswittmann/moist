!> LPB component backed by generated Gaussian and point kernels
module moist_model_continuum_component_pcm_lpb
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_layer, only: pcm_layer_base, pcm_layer_feps
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver
   use pcmgen_gauss_lpb, only: lpb_complex_state, &
      & lpb_assemble, lpb_medium, lpb_solve, lpb_weights, lpb_tangent
   use pcmgen_point_lpb, only: lpb_point_complex_state, &
      & lpb_point_assemble, lpb_point_medium, lpb_point_solve, &
      & lpb_point_weights, lpb_point_tangent
   implicit none(type, external)
   private

   public :: model_continuum_component_lpb, new_component_lpb
   public :: debye_kappa

   !> Linearized Poisson-Boltzmann component
   !>
   !> - Physical inputs: solvent relative permittivity epsilon and inverse Debye length kappa
   !> - Cavity positions, normals, principal curvatures and grid point areas
   !> - Gaussian charge widths and switching factors when available
   !> - Host electrostatic potential at the surface points
   !> - Host nuclear charges for gradients, with ECP core electrons removed
   type, extends(pcm_layer_base) :: model_continuum_component_lpb
      !> Gaussian surface kernel state
      type(lpb_complex_state) :: gauss
      !> Point surface kernel state
      type(lpb_point_complex_state) :: point
      !> Inverse Debye screening length, bohr^-1
      real(wp) :: kappa = 0.0_wp
   contains
      procedure :: assemble_layer => lpb_assemble_layer
      procedure :: solve_layer => lpb_solve_layer
      procedure :: layer_weights => lpb_layer_weights
      procedure :: layer_tangent => lpb_layer_tangent
      procedure :: is_inert => lpb_is_inert
   end type model_continuum_component_lpb

contains

   !> Construct LPB from physical and numerical parameters
   !>
   !> - Kappa from ionic strength and temperature: debye_kappa
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self     initialized LPB component
   !> @param[in]  epsilon  dimensionless solvent relative permittivity, >= 1
   !> @param[in]  kappa    inverse Debye length in bohr^-1, finite and >= 0
   !> @param[out] error    allocated on invalid physical or numerical parameters
   !> @param[in]  param    optional solver controls, compiled defaults when omitted
   !> @param[in]  ctx      optional borrowed context, must outlive the component
   subroutine new_component_lpb(self, epsilon, kappa, error, param, ctx)
      !> LPB component
      type(model_continuum_component_lpb), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Solvent relative permittivity
      real(wp), intent(in) :: epsilon
      !> Inverse Debye screening length
      real(wp), intent(in) :: kappa
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear solver and operator-storage controls
      type(moist_pcm_parameters_type), intent(in), optional :: param
      type(moist_pcm_parameters_type) :: settings

      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, "LPB epsilon must be >= 1")
         return
      end if
      if (.not. ieee_is_finite(kappa) .or. kappa < 0.0_wp) then
         call fatal_error(error, "LPB kappa must be finite and nonnegative")
         return
      end if
      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      if (present(ctx)) self%ctx => ctx
      self%epsilon = epsilon
      self%feps = pcm_layer_feps(epsilon)
      self%param = settings
      self%name = "LPB"
      self%kappa = kappa
   end subroutine new_component_lpb

   !> Assemble the LPB operators
   subroutine lpb_assemble_layer(self, cavity, solver, error)
      class(model_continuum_component_lpb), intent(inout) :: self
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
               call lpb_assemble(self%gauss%lpb_state, r, n, k1, k2, xi, f, a, error, solver)
               if (allocated(error)) return
               call lpb_medium(self%gauss%lpb_state, self%epsilon, self%kappa, error)
            end associate
         else
            call lpb_point_assemble(self%point%lpb_point_state, r, n, k1, k2, a, error, solver)
            if (allocated(error)) return
            call lpb_point_medium(self%point%lpb_point_state, self%epsilon, self%kappa, error)
         end if
      end associate
   end subroutine lpb_assemble_layer

   !> Solve the current LPB host traces
   subroutine lpb_solve_layer(self, error)
      class(model_continuum_component_lpb), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      if (self%smeared) then
         call lpb_solve(self%gauss%lpb_state, self%phi, self%kernel_energy, self%q, self%qsym, error)
      else
         call lpb_point_solve(self%point%lpb_point_state, self%phi, self%kernel_energy, self%q, self%qsym, error)
      end if
   end subroutine lpb_solve_layer

   !> Accumulate LPB operator surface weights
   subroutine lpb_layer_weights(self, cavity, acc, error)
      class(model_continuum_component_lpb), intent(inout) :: self
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_n(:, :), w_k1(:), w_k2(:)
      real(wp), allocatable :: w_a(:), w_xi(:), w_f(:)
      real(wp) :: w_eps
      real(wp) :: w_kappa
      integer :: ngrid

      ngrid = cavity%ngrid
      allocate (w_xyz(3, ngrid), w_n(3, ngrid), w_k1(ngrid), w_k2(ngrid), w_a(ngrid))
      if (self%smeared) then
         allocate (w_xi(ngrid), w_f(ngrid), source=0.0_wp)
         call lpb_weights(self%gauss%lpb_state, w_xyz, w_n, w_k1, w_k2, w_xi, w_f, w_a, w_eps, w_kappa, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a, w_xi=w_xi, w_f=w_f)
      else
         call lpb_point_weights(self%point%lpb_point_state, w_xyz, w_n, w_k1, w_k2, w_a, w_eps, w_kappa, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_n=w_n, w_k1=w_k1, &
            & w_k2=w_k2, w_a=w_a)
      end if
   end subroutine lpb_layer_weights

   !> Directional derivative of the LPB weights and trace adjoints
   !>
   !> - no normal trace: `dvn` is ignored and `hvp_vn` stays zero
   !>
   !> @param[in,out] self   LPB component with current weights
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
   subroutine lpb_layer_tangent(self, dr, dn, dk1, dk2, dxi, df, da, dv, dvn, hvp_r, hvp_n, &
         & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, hvp_vn, error)
      !> LPB component with current weights
      class(model_continuum_component_lpb), intent(inout) :: self
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
         call lpb_tangent(self%gauss%lpb_state, dr, dn, dk1, dk2, dxi, df, da, dv, hvp_r, hvp_n, &
            & hvp_k1, hvp_k2, hvp_xi, hvp_f, hvp_a, hvp_v, error)
      else
         call lpb_point_tangent(self%point%lpb_point_state, dr, dn, dk1, dk2, da, dv, hvp_r, hvp_n, &
            & hvp_k1, hvp_k2, hvp_a, hvp_v, error)
      end if
   end subroutine lpb_layer_tangent

   !> Screening can polarize even when epsilon is one
   pure function lpb_is_inert(self) result(inert)
      class(model_continuum_component_lpb), intent(in) :: self
      logical :: inert
      inert = self%feps == 0.0_wp .and. self%kappa == 0.0_wp
   end function lpb_is_inert

   !> Inverse Debye screening length of a 1:1 electrolyte, in inverse bohr
   pure function debye_kappa(epsilon, ionic_strength, temperature) result(kappa)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use mctc_io_constants, only: pi, codata
      use mctc_io_convert, only: autoaa, autokj
      real(wp), intent(in) :: epsilon, ionic_strength, temperature
      real(wp) :: kappa, density, kt

      if (.not. (epsilon > 0.0_wp .and. ionic_strength >= 0.0_wp .and. temperature > 0.0_wp)) then
         kappa = ieee_value(kappa, ieee_quiet_nan)
         return
      end if
      density = ionic_strength*1.0e3_wp*codata%NA*(autoaa*1.0e-10_wp)**3
      kt = codata%kb*temperature*codata%NA*1.0e-3_wp/autokj
      kappa = sqrt(8.0_wp*pi*density/(epsilon*kt))
   end function debye_kappa

end module moist_model_continuum_component_pcm_lpb
