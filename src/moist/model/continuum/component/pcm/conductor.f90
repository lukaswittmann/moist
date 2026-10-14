!> Generated conductor kernel adapter for CPCM and COSMO
module moist_model_continuum_component_pcm_conductor
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use moist_context, only: moist_context_type
   use moist_utils_timer, only: cat_setup, cat_solve
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & coupling_request_type, gaussian_potential_request_type, &
      & point_potential_request_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_register, request_require
   use moist_model_continuum_component_pcm_type, only: model_continuum_component_pcm, &
      & solver_type, moist_pcm_parameters_type
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver, pcm_dense_solver, &
      & pcm_gmres_solver
   use pcmgen_gauss_conductor, only: conductor_state, &
      & conductor_assemble, conductor_medium, conductor_solve, conductor_weights
   use pcmgen_point_conductor, only: conductor_point_state, &
      & conductor_point_assemble, conductor_point_medium, conductor_point_solve, &
      & conductor_point_weights
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none(type, external)
   private

   public :: model_continuum_component_pcm_conductor
   public :: model_continuum_component_cpcm, new_component_cpcm
   public :: model_continuum_component_cosmo, new_component_cosmo

   !> Conductor-like PCM backed by the generated Gaussian or point kernel
   type, abstract, extends(model_continuum_component_pcm) :: model_continuum_component_pcm_conductor
      !> Gaussian surface kernel state
      type(conductor_state) :: gauss
      !> Point surface kernel state
      type(conductor_point_state) :: point
      !> Whether the assembled surface uses Gaussian charges
      logical :: smeared = .true.
   contains
      procedure :: update => conductor_update
      procedure :: ensure_charges => conductor_ensure_charges
      procedure :: get_surface_weights => conductor_get_surface_weights
      procedure :: get_gradient_surface_weights => conductor_get_gradient_surface_weights
      procedure :: get_host_surface_weights => conductor_get_host_surface_weights
      procedure :: declare_coupling => conductor_declare_coupling
   end type model_continuum_component_pcm_conductor

   !> CPCM with dielectric factor (epsilon - 1)/epsilon
   !>
   !> - Physical input: solvent relative permittivity epsilon
   !> - Cavity positions
   !> - Gaussian widths and switching factors, or point-charge grid point areas
   !> - Host electrostatic potential at the surface points
   type, extends(model_continuum_component_pcm_conductor) :: model_continuum_component_cpcm
   end type model_continuum_component_cpcm

   !> COSMO with dielectric factor (epsilon - 1)/(epsilon + 1/2)
   !>
   !> - Physical input: solvent relative permittivity epsilon
   !> - Cavity positions
   !> - Gaussian widths and switching factors, or point-charge grid point areas
   !> - Host electrostatic potential at the surface points
   type, extends(model_continuum_component_pcm_conductor) :: model_continuum_component_cosmo
   end type model_continuum_component_cosmo

contains

   !> Construct CPCM from physical and numerical parameters
   !>
   !> - Positive infinite epsilon: conductor limit
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self     initialized CPCM component
   !> @param[in]  epsilon  dimensionless solvent relative permittivity, >= 1
   !> @param[out] error    allocated on invalid physical or numerical parameters
   !> @param[in]  param    optional solver controls, compiled defaults when omitted
   !> @param[in]  ctx      optional borrowed context, must outlive the component
   subroutine new_component_cpcm(self, epsilon, error, param, ctx)
      !> CPCM component
      type(model_continuum_component_cpcm), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Solvent relative permittivity
      real(wp), intent(in) :: epsilon
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear solver and operator-storage controls
      type(moist_pcm_parameters_type), intent(in), optional :: param

      if (present(ctx)) self%ctx => ctx
      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, "[new_component_cpcm] Dielectric constant must be >= 1")
         return
      end if
      self%epsilon = epsilon
      if (ieee_is_finite(epsilon)) then
         self%feps = (epsilon - 1.0_wp)/epsilon
      else
         self%feps = 1.0_wp
      end if
      call configure_conductor(self, param, error)
      if (allocated(error)) return
      self%name = "CPCM"
      self%description = "Conductor-like polarizable continuum, f(eps) = (eps - 1)/eps"
   end subroutine new_component_cpcm

   !> Construct COSMO from physical and numerical parameters
   !>
   !> - Positive infinite epsilon: conductor limit
   !> - Standalone computation requires ctx; model-owned components inherit the model context
   !>
   !> @param[out] self     initialized COSMO component
   !> @param[in]  epsilon  dimensionless solvent relative permittivity, >= 1
   !> @param[out] error    allocated on invalid physical or numerical parameters
   !> @param[in]  param    optional solver controls, compiled defaults when omitted
   !> @param[in]  ctx      optional borrowed context, must outlive the component
   subroutine new_component_cosmo(self, epsilon, error, param, ctx)
      !> COSMO component
      type(model_continuum_component_cosmo), intent(out) :: self
      !> Borrowed run context; model context when omitted
      type(moist_context_type), intent(in), target, optional :: ctx
      !> Solvent relative permittivity
      real(wp), intent(in) :: epsilon
      !> Construction error
      type(error_type), allocatable, intent(out) :: error
      !> Linear solver and operator-storage controls
      type(moist_pcm_parameters_type), intent(in), optional :: param

      if (present(ctx)) self%ctx => ctx
      if (epsilon < 1.0_wp .or. epsilon /= epsilon) then
         call fatal_error(error, "[new_component_cosmo] Dielectric constant must be >= 1")
         return
      end if
      self%epsilon = epsilon
      if (ieee_is_finite(epsilon)) then
         self%feps = (epsilon - 1.0_wp)/(epsilon + 0.5_wp)
      else
         self%feps = 1.0_wp
      end if
      call configure_conductor(self, param, error)
      if (allocated(error)) return
      self%name = "COSMO"
      self%description = "Conductor-like screening model, f(eps) = (eps - 1)/(eps + 1/2)"
   end subroutine new_component_cosmo

   !> Copy validated numerical settings to a conductor component
   subroutine configure_conductor(self, param, error)
      class(model_continuum_component_pcm_conductor), intent(inout) :: self
      type(moist_pcm_parameters_type), intent(in), optional :: param
      type(error_type), allocatable, intent(out) :: error
      type(moist_pcm_parameters_type) :: settings

      if (present(param)) settings = param
      call settings%validate(error)
      if (allocated(error)) return
      self%param = settings
   end subroutine configure_conductor

   !> Register the potential primitive matching the conductor discretization
   subroutine conductor_declare_coupling(self, cavity, coupling, error)
      class(model_continuum_component_pcm_conductor), intent(in) :: self
      class(cavity_type), intent(in) :: cavity
      type(coupling_type), intent(inout) :: coupling
      type(error_type), allocatable, intent(out) :: error
      class(coupling_request_type), allocatable :: potential
      logical :: gaussian, field

      gaussian = allocated(cavity%xi0)
      field = cavity%has_field_dependent_geometry()
      if (gaussian) then
         allocate (gaussian_potential_request_type :: potential)
      else
         allocate (point_potential_request_type :: potential)
      end if
      call request_require(potential, moist_phase_energy, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "dphi_dr", field, error)
      if (allocated(error)) return
      if (gaussian) then
         call request_require(potential, moist_phase_response, "dphi_dxi", field, error)
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
      call coupling_register(coupling, "potential", potential, error)
   end subroutine conductor_declare_coupling

   !> Assemble the conductor operator for the current cavity
   subroutine conductor_update(self, mol, cavity, error)
      class(model_continuum_component_pcm_conductor), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      class(cavity_type), intent(inout) :: cavity
      type(error_type), allocatable, intent(out) :: error
      class(pcm_solver), allocatable :: solver
      integer :: d0

      call self%require_context(error)
      if (allocated(error)) return
      self%mol_solu = mol
      self%charges_valid = .false.
      if (allocated(self%phi)) deallocate (self%phi)
      if (allocated(self%q)) deallocate (self%q)
      if (.not. allocated(cavity%xyz)) then
         call fatal_error(error, "Conductor PCM requires surface positions")
         return
      end if

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("PCM setup", category=cat_setup)

      if (self%param%solver == solver_type%iterative) then
         allocate (pcm_gmres_solver :: solver)
         select type (solver)
         type is (pcm_gmres_solver)
            solver%tol = self%param%solver_tol
            solver%max_iter = self%param%solver_maxiter
            solver%restart = self%param%solver_restart
            solver%storage = self%param%storage
            solver%memory_budget_gib = self%param%memory_budget_gib
         end select
      else
         allocate (pcm_dense_solver :: solver)
      end if
      solver%nonlinear_tol = self%param%newton_tol
      solver%nonlinear_max_iter = self%param%newton_maxiter

      self%smeared = allocated(cavity%xi0) .and. allocated(cavity%f)
      if (self%smeared) then
         call self%ctx%timer%start("Interaction matrix")
         call conductor_assemble(self%gauss, cavity%xyz, cavity%xi0, cavity%f, error, solver)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Interaction matrix")
         call conductor_medium(self%gauss, self%feps, error)
      else if (allocated(cavity%a)) then
         call self%ctx%timer%start("Interaction matrix")
         call conductor_point_assemble(self%point, cavity%xyz, cavity%a, error, solver)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Interaction matrix")
         call conductor_point_medium(self%point, self%feps, error)
      else
         call fatal_error(error, "Conductor PCM requires Gaussian widths or surface areas")
      end if
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      allocate (self%q(cavity%ngrid), source=0.0_wp)
      call self%ctx%timer%stop("PCM setup")
   end subroutine conductor_update

   !> Solve the generated conductor system when the host potential changes
   subroutine conductor_ensure_charges(self, coupling, cavity, error)
      class(model_continuum_component_pcm_conductor), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: phi(:), grad_v(:)
      real(wp) :: energy
      integer :: d0

      if (.not. allocated(self%q)) then
         call fatal_error(error, "Conductor PCM must be updated before solving")
         return
      end if
      call coupling%read("potential", "phi", phi, error)
      if (allocated(error)) return
      if (size(phi) /= cavity%ngrid) then
         call fatal_error(error, "Conductor PCM potential size mismatch")
         return
      end if
      if (allocated(self%phi)) then
         if (size(self%phi) == size(phi)) then
            if (all(self%phi == phi) .and. self%charges_valid) return
         end if
      end if
      self%phi = phi
      self%charges_valid = .false.
      allocate (grad_v(cavity%ngrid))
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("PCM solve", category=cat_solve)
      if (self%smeared) then
         call conductor_solve(self%gauss, self%phi, energy, self%q, grad_v, error)
      else
         call conductor_point_solve(self%point, self%phi, energy, self%q, grad_v, error)
      end if
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if
      call self%ctx%timer%stop("PCM solve")
      self%charges_valid = .true.
   end subroutine conductor_ensure_charges

   !> Add generated operator weights to the cavity adjoint
   subroutine conductor_get_surface_weights(self, coupling, cavity, acc, error)
      class(model_continuum_component_pcm_conductor), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_a(:), w_xi(:), w_f(:)
      real(wp) :: w_feps

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      if (self%is_inert()) return
      allocate (w_xyz(3, cavity%ngrid))
      if (self%smeared) then
         allocate (w_xi(cavity%ngrid), w_f(cavity%ngrid))
         call conductor_weights(self%gauss, w_xyz, w_xi, w_f, w_feps, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_xi=w_xi, w_f=w_f)
      else
         allocate (w_a(cavity%ngrid))
         call conductor_point_weights(self%point, w_xyz, w_a, w_feps, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_a=w_a)
      end if
      if (allocated(error)) return
      if (cavity%has_field_dependent_geometry()) then
         call self%get_host_surface_weights(coupling, acc, cavity%ngrid, error)
      end if
   end subroutine conductor_get_surface_weights

   !> Include host trace geometry weights in the nuclear gradient
   subroutine conductor_get_gradient_surface_weights(self, coupling, cavity, acc, error)
      class(model_continuum_component_pcm_conductor), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error

      call self%get_surface_weights(coupling, cavity, acc, error)
      if (allocated(error)) return
      if (.not. cavity%has_field_dependent_geometry()) then
         call self%get_host_surface_weights(coupling, acc, cavity%ngrid, error)
      end if
   end subroutine conductor_get_gradient_surface_weights

   !> Contract host potential derivatives for Gaussian or point charges
   subroutine conductor_get_host_surface_weights(self, coupling, acc, ngrid, error)
      class(model_continuum_component_pcm_conductor), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      integer, intent(in) :: ngrid
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_phi(:), w_xyz(:, :), w_xi(:)
      integer :: i

      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call coupling%read("potential", "dphi_dr", w_xyz, error)
      if (allocated(error)) return
      if (any(shape(w_xyz) /= [3, ngrid])) then
         call fatal_error(error, "PCM spatial derivative grid mismatch")
         return
      end if
      do i = 1, ngrid
         w_xyz(:, i) = w_phi(i)*w_xyz(:, i)
      end do
      if (self%smeared) then
         call coupling%read("potential", "dphi_dxi", w_xi, error)
         if (allocated(error)) return
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_xi=w_phi*w_xi)
      else
         call acc%add_surface_weights(error, w_xyz=w_xyz)
      end if
   end subroutine conductor_get_host_surface_weights

end module moist_model_continuum_component_pcm_conductor
