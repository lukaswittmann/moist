!> PCM (Polarizable Continuum Model) abstract base type
!>
!> The abstract PCM base component, extended by the concrete implementations
!> (CPCM, COSMO, ...) in their respective modules
module moist_model_continuum_component_pcm_type
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env, only: wp, fatal_error
   use mctc_env_error, only: error_type
   use mctc_io, only: structure_type
   use moist_cavity_type, only: cavity_type
   use moist_model_continuum_component_type, only: model_continuum_component_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type, &
      & response_accumulate
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & gaussian_potential_request_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_register, request_require
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_electrostatics, only: &
      & pcm_electrostatic_direct_gradient, &
      & pcm_electrostatic_normal_direct_gradient
   use moist_utils_timer, only: cat_energy
   use moist_model_continuum_component_pcm_solvers, only: pcm_storage_auto, pcm_storage_dense, &
      & pcm_storage_matrix_free
   use moist_utils_prettyprint, only: prettyprinter
   implicit none(type, external)
   private

   public :: model_continuum_component_pcm
   public :: pcm_solver_type
   public :: solver_type

   !> PCM linear solver choices
   type :: pcm_solver_type
      !> Dense factorization
      integer :: dense = 1
      !> Restarted GMRES with block-Jacobi preconditioning
      integer :: iterative = 2
   end type pcm_solver_type

   !> Global instance of the solver enum
   type(pcm_solver_type), parameter :: solver_type = pcm_solver_type()

   public :: moist_pcm_parameters_type

   !> PCM solver configuration
   type, extends(moist_model_parameters_type) :: moist_pcm_parameters_type
      !> Linear solver selection
      integer :: solver = solver_type%dense
      !> Iterative solver convergence tolerance
      real(wp) :: solver_tol = 1.0e-8_wp
      !> Maximum iterative solver steps
      integer :: solver_maxiter = 1000
      !> Krylov dimension before GMRES restarts
      integer :: solver_restart = 100
      !> Operator storage choice for iterative solves
      integer :: storage = pcm_storage_auto
      !> Memory budget in GiB for automatic operator storage
      integer :: memory_budget_gib = 4
      !> Relative nonlinear residual threshold
      real(wp) :: newton_tol = 1.0e-10_wp
      !> Maximum Newton steps
      integer :: newton_maxiter = 50
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
      !> Check solver controls
      procedure :: validate => validate_pcm_parameters
   end type moist_pcm_parameters_type

   !> Abstract PCM base component
   !>
   !> Common infrastructure for PCM-family methods (CPCM, COSMO, IEF-PCM)
   !>
   !> Owns the host coupling and gradient lifecycle shared by generated PCM models
   type, abstract, extends(model_continuum_component_type) :: model_continuum_component_pcm

      !> Dielectric constant of the solvent
      real(wp) :: epsilon

      !> Dielectric scaling factor f( epsilon ) - variant-specific formula
      !>
      !> CPCM: f epsilon = ( epsilon -1)/ epsilon
      !> COSMO: f epsilon = ( epsilon -1)/( epsilon +0.5)
      real(wp) :: feps

      !> Surface charges on cavity grid points (ngrid)
      real(wp), allocatable :: q(:)

      !> Numerical settings shared by computation and parameter introspection
      type(moist_pcm_parameters_type) :: param

      !> Molecular electrostatic potential at cavity grid points (ngrid)
      !>
      !> Copied from the coupling's potential request by `ensure_charges`;
      !> never computed here
      real(wp), allocatable :: phi(:)

      !> Whether charges belong to the current cavity and host potential
      logical :: charges_valid = .false.

   contains

      !> Request normals when the model uses a normal potential trace
      procedure :: configure_cavity => pcm_component_configure_cavity

      !> Whether this component reads n dot grad(phi)
      procedure :: uses_normal_trace => pcm_component_uses_normal_trace

      !> Whether the current medium has no electrostatic response
      procedure :: is_inert => pcm_component_is_inert

      !> Adjoint of n dot grad(phi), zero for conductor models
      procedure :: normal_trace_adjoint => pcm_component_normal_trace_adjoint

      !> Read n dot grad(phi) from raw host potential derivatives
      procedure :: normal_trace => pcm_component_normal_trace

      !> Update a concrete PCM component
      procedure :: update => pcm_component_update

      !> Compute PCM solvation energy
      procedure :: get_energy => pcm_component_get_energy

      !> Compute PCM reaction potential
      procedure :: get_response => pcm_component_get_response

      !> Compute the potential adjoint `dE/dphi` (the surface charges for a stationary PCM)
      procedure :: get_trace_response => pcm_component_get_trace_response

      !> Compute PCM gradient with respect to nuclear coordinates
      procedure :: get_gradient => pcm_component_get_gradient

      !> Contract the explicit nuclear Jacobians of Gaussian surface geometry
      procedure :: get_forward_gradient => pcm_component_get_forward_gradient

      !> Accumulate the PCM cavity surface adjoint weights
      procedure :: get_surface_weights => pcm_component_get_surface_weights

      !> Accumulate the surface adjoints the reverse-mode gradient consumes
      procedure :: get_gradient_surface_weights => pcm_component_get_gradient_surface_weights

      !> Nuclear gradient contributions that bypass the cavity surface
      procedure :: get_direct_gradient => pcm_component_get_direct_gradient

      !> Nuclear charges of the solute, the moving sources of the potential
      procedure :: nuclear_charges => pcm_component_nuclear_charges

      !> Fold the host's direct trace-geometry weights into the accumulator
      procedure :: get_host_surface_weights => pcm_component_get_host_surface_weights

      !> Declare the potential and the host surface weights
      procedure :: declare_coupling => pcm_component_declare_coupling

      !> Solve for the surface charges unless they are already current
      procedure :: ensure_charges => pcm_ensure_charges

      !> Potential adjoint `w_phi = dE/dphi` of the current charges
      procedure :: potential_adjoint => pcm_component_potential_adjoint

      !> Linear solver settings
      procedure :: parameters => pcm_component_parameters

      !> Print the dielectric constant and its scaling factor
      procedure :: print_inputs => pcm_component_print_inputs

   end type model_continuum_component_pcm

contains

   !> Request normals before building a cavity for a normal-trace model
   subroutine pcm_component_configure_cavity(self, cavity, error)
      class(model_continuum_component_pcm), intent(inout) :: self
      class(cavity_type), intent(inout) :: cavity
      type(error_type), allocatable, intent(out) :: error

      if (self%uses_normal_trace()) call cavity%configure(normal=.true., error=error)
   end subroutine pcm_component_configure_cavity

   !> Conductor models read only the potential value
   pure function pcm_component_uses_normal_trace(self) result(uses)
      class(model_continuum_component_pcm), intent(in) :: self
      logical :: uses

      uses = .false.
   end function pcm_component_uses_normal_trace

   !> A conductor-like model is inert at zero dielectric factor
   pure function pcm_component_is_inert(self) result(inert)
      class(model_continuum_component_pcm), intent(in) :: self
      logical :: inert

      inert = self%feps == 0.0_wp
   end function pcm_component_is_inert

   !> Default zero adjoint of the normal potential trace
   subroutine pcm_component_normal_trace_adjoint(self, w_vn, error)
      class(model_continuum_component_pcm), intent(in) :: self
      real(wp), allocatable, intent(out) :: w_vn(:)
      type(error_type), allocatable, intent(out) :: error

      if (.not. allocated(self%q)) then
         call fatal_error(error, "PCM charges are unavailable")
         return
      end if
      allocate (w_vn(size(self%q)), source=0.0_wp)
   end subroutine pcm_component_normal_trace_adjoint

   !> Form n dot grad(phi) from the host's raw spatial derivative
   subroutine pcm_component_normal_trace(self, coupling, cavity, vn, error)
      class(model_continuum_component_pcm), intent(in) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      real(wp), allocatable, intent(out) :: vn(:)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: grad_phi(:, :)
      integer :: i

      if (.not. allocated(cavity%normal0)) then
         call fatal_error(error, "PCM normal trace requires cavity normals")
         return
      end if
      call coupling%read("potential", "dphi_dr", grad_phi, error)
      if (allocated(error)) return
      if (any(shape(grad_phi) /= [3, cavity%ngrid])) then
         call fatal_error(error, "PCM normal trace derivative shape mismatch")
         return
      end if
      allocate (vn(cavity%ngrid))
      do i = 1, cavity%ngrid
         vn(i) = dot_product(cavity%normal0(:, i), grad_phi(:, i))
      end do
   end subroutine pcm_component_normal_trace

   !> Validate PCM solver settings
   !>
   !> @param[inout] self Solver settings
   !> @param[out] error Invalid setting
   subroutine validate_pcm_parameters(self, error)
      class(moist_pcm_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      if (self%solver < solver_type%dense .or. self%solver > solver_type%iterative .or. &
          self%solver_tol <= 0.0_wp .or. self%solver_tol /= self%solver_tol .or. &
          self%solver_maxiter < 1 .or. self%solver_restart < 1 .or. &
          (self%storage /= pcm_storage_auto .and. self%storage /= pcm_storage_dense .and. &
           self%storage /= pcm_storage_matrix_free) .or. &
          self%memory_budget_gib < 1 .or. self%newton_tol <= 0.0_wp .or. &
          self%newton_tol /= self%newton_tol .or. self%newton_maxiter < 1) then
         call fatal_error(error, "Invalid PCM solver parameters")
      end if
   end subroutine validate_pcm_parameters

   !> Restore compiled parameter defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      class(moist_pcm_parameters_type), intent(inout) :: self
      type(moist_pcm_parameters_type) :: defaults

      self%solver = defaults%solver
      self%solver_tol = defaults%solver_tol
      self%solver_maxiter = defaults%solver_maxiter
      self%solver_restart = defaults%solver_restart
      self%storage = defaults%storage
      self%memory_budget_gib = defaults%memory_budget_gib
      self%newton_tol = defaults%newton_tol
      self%newton_maxiter = defaults%newton_maxiter
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_pcm_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("solver", self%solver)
      call self%register_real_scalar("solver_tol", self%solver_tol)
      call self%register_int_scalar("solver_maxiter", self%solver_maxiter)
      call self%register_int_scalar("solver_restart", self%solver_restart)
      call self%register_int_scalar("storage", self%storage)
      call self%register_int_scalar("memory_budget_gib", self%memory_budget_gib)
      call self%register_real_scalar("newton_tol", self%newton_tol)
      call self%register_int_scalar("newton_maxiter", self%newton_maxiter)
   end subroutine register_parameter_entries

   !> Linear solver settings of the component
   !>
   !> @param[in] self PCM component instance
   function pcm_component_parameters(self) result(param)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(in), target :: self
      !> Solver settings, valid while the component is
      class(moist_model_parameters_type), pointer :: param

      param => self%param

   end function pcm_component_parameters

   !> Print the dielectric constant and its scaling factor
   !>
   !> @param[in]    self PCM component instance
   !> @param[inout] pp   Pretty printer inside the component section
   subroutine pcm_component_print_inputs(self, pp)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp

      call pp%kv("Epsilon", self%epsilon)
      call pp%kv("f(eps)", self%feps)

   end subroutine pcm_component_print_inputs

   !> Report a concrete PCM implementation that omitted its update hook
   subroutine pcm_component_update(self, mol, cavity, error)
      class(model_continuum_component_pcm), intent(inout) :: self
      type(structure_type), intent(in) :: mol
      class(cavity_type), intent(inout) :: cavity
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "PCM component does not implement update")
   end subroutine pcm_component_update

   !> Report a concrete PCM implementation that omitted its solve hook
   subroutine pcm_ensure_charges(self, coupling, cavity, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Host data carrying the molecular potential trace
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "PCM component does not implement charge solving")
   end subroutine pcm_ensure_charges

   !> Potential adjoint `w_phi = dE/dphi` of the current charges
   !>
   !> Every chain-rule contraction through the host potential (Fock weights,
   !> the direct nuclear term, the surface-position and width weights) uses
   !> this vector, never the charges themselves. For a stationary PCM (CPCM,
   !> COSMO) it is the surface charge `q`; a non-symmetric response matrix
   !> (IEF-PCM, SS(V)PE) overrides it with the symmetrized adjoint
   !>
   !> @param[in]  self  PCM component with current surface charges
   !> @param[out] w_phi Potential adjoint (ngrid)
   !> @param[out] error Charges not current
   subroutine pcm_component_potential_adjoint(self, w_phi, error)
      !> PCM component with current surface charges
      class(model_continuum_component_pcm), intent(in) :: self
      !> Potential adjoint (ngrid)
      real(wp), allocatable, intent(out) :: w_phi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%charges_valid .or. .not. allocated(self%q)) then
         call fatal_error(error, "pcm_component_potential_adjoint: "// &
            & "surface charges are unavailable - call ensure_charges first")
         return
      end if
      w_phi = self%q

   end subroutine pcm_component_potential_adjoint

   !> Compute PCM solvation energy
   !>
   !> E_solv = 0.5 * dot(q, phi)
   !>
   !> The molecular potential phi comes from the coupling's potential request
   !>
   !> A missing potential fails its read by name; completeness is the model's
   !> entry check
   !>
   !> @param[in,out] self PCM component instance
   !> @param[in] coupling Host data carrying the molecular potential trace
   !> @param[in,out] cavity Live cavity owned by the orchestrating model
   !> @param[in,out] energy Solvation energy (inout to allow accumulation)
   !> @param[out] error Error handling
   subroutine pcm_component_get_energy(self, coupling, cavity, energy, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Host data carrying the molecular potential trace
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Solvation energy (inout to allow accumulation)
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: e_pcm
      !> Timer depth on entry, restored on every early return
      integer :: d0

      call self%require_context(error)
      if (allocated(error)) return

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("PCM energy", category=cat_energy)

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) then
         call self%ctx%timer%unwind(d0)
         return
      end if

      ! Compute energy: E = 0.5 * dot(q, phi)
      e_pcm = 0.5_wp*dot_product(self%q, self%phi)

      ! Accumulate into energy
      energy = energy + e_pcm

      call self%ctx%timer%stop("PCM energy")

   end subroutine pcm_component_get_energy

   !> Compute the PCM potential adjoint
   !>
   !>    dE/dphi_i = q_i
   !> by stationarity, returned as the `potential_adjoint_response_type` item; the
   !> host contracts it with its own potential integrals to build the Fock
   !> contribution F_uv += sum_i q_i V_uv(r_i)
   !>
   !> @param[in,out] self PCM component instance
   !> @param[in] coupling Wavefunction data
   !> @param[in,out] cavity Live cavity owned by the orchestrating model
   !> @param[in,out] response Solvation response
   !> @param[out] error Error handling
   subroutine pcm_component_get_response(self, coupling, cavity, response, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Solvation response
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%require_context(error)
      if (allocated(error)) return

      call self%get_trace_response(coupling, cavity, response, error)

   end subroutine pcm_component_get_response

   !> Accumulate the direct electrostatic trace adjoint
   !>
   !> @param[in,out] self      PCM component
   !> @param[in]    coupling  Host coupling data
   !> @param[in,out] cavity    Live model cavity
   !> @param[in,out] response  Direct trace-response accumulator
   !> @param[out]   error     Error handling
   subroutine pcm_component_get_trace_response(self, coupling, cavity, response, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Direct trace-response accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Potential adjoint item of this component
      type(potential_adjoint_response_type) :: item
      real(wp), allocatable :: w_vn(:)
      integer :: i

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      call self%potential_adjoint(item%w_phi, error)
      if (allocated(error)) return
      if (self%uses_normal_trace()) then
         if (.not. allocated(cavity%normal0)) then
            call fatal_error(error, "PCM normal trace requires cavity normals")
            return
         end if
         call self%normal_trace_adjoint(w_vn, error)
         if (allocated(error)) return
         if (size(w_vn) /= cavity%ngrid) then
            call fatal_error(error, "PCM normal adjoint size mismatch")
            return
         end if
         allocate (item%w_dphi_dr(3, cavity%ngrid))
         do i = 1, cavity%ngrid
            item%w_dphi_dr(:, i) = w_vn(i)*cavity%normal0(:, i)
         end do
      end if
      call response_accumulate(response, item, error)

   end subroutine pcm_component_get_trace_response

   !> Contract the reverse cavity adjoint and fixed-surface nuclear terms
   !>
   !> @param[in,out] self PCM component
   !> @param[in] coupling Scoped host data
   !> @param[in,out] cavity Live model cavity
   !> @param[in,out] response Host trace-response accumulator
   !> @param[in,out] gradient Nuclear gradient (3, nat), unchanged on failure
   !> @param[out] error Allocated on failure
   subroutine pcm_component_get_gradient(self, coupling, cavity, response, gradient, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Scoped host data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host trace-response accumulator
      type(response_type), intent(inout) :: response
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call pcm_component_gradient(self, coupling, cavity, response, gradient, error, forward=.false.)
   end subroutine pcm_component_get_gradient

   !> Contract explicit Gaussian surface Jacobians and fixed-surface terms
   !>
   !> @param[in,out] self PCM component
   !> @param[in] coupling Scoped host data
   !> @param[in,out] cavity Live model cavity
   !> @param[in,out] response Host trace-response accumulator
   !> @param[in,out] gradient Nuclear gradient (3, nat), unchanged on failure
   !> @param[out] error Allocated on failure
   subroutine pcm_component_get_forward_gradient(self, coupling, cavity, response, gradient, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Scoped host data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host trace-response accumulator
      type(response_type), intent(inout) :: response
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call pcm_component_gradient(self, coupling, cavity, response, gradient, error, forward=.true.)
   end subroutine pcm_component_get_forward_gradient

   !> Assemble a PCM gradient locally and commit only after success
   !>
   !> @param[in,out] self PCM component
   !> @param[in] coupling Scoped host data
   !> @param[in,out] cavity Live model cavity
   !> @param[in,out] response Host trace-response accumulator
   !> @param[in,out] gradient Nuclear gradient (3, nat), unchanged on failure
   !> @param[out] error Allocated on failure
   !> @param[in] forward Select explicit Jacobians instead of the reverse cavity pass
   subroutine pcm_component_gradient(self, coupling, cavity, response, gradient, error, forward)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Scoped host data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Host trace-response accumulator
      type(response_type), intent(inout) :: response
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Select explicit cavity Jacobians
      logical, intent(in) :: forward
      !> Component surface adjoint
      type(cavity_surface_adjoint_type) :: acc
      !> Local gradient contribution, committed on success
      real(wp), allocatable :: contribution(:, :)
      integer :: stat

      call self%require_context(error)
      if (allocated(error)) return
      call self%get_trace_response(coupling, cavity, response, error)
      if (allocated(error)) return
      allocate (contribution(size(gradient, 1), size(gradient, 2)), source=0.0_wp, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "[pcm_component_get_gradient] gradient allocation failed")
         return
      end if
      call self%get_direct_gradient(coupling, cavity, contribution, error)
      if (allocated(error)) return
      call acc%init(cavity%ngrid)
      call self%get_gradient_surface_weights(coupling, cavity, acc, error)
      if (allocated(error)) return
      if (forward) then
         call pcm_forward_surface_gradient(cavity, acc, contribution, error)
      else
         call cavity%get_surface_gradient(acc, contribution, error)
      end if
      if (allocated(error)) return
      gradient = gradient + contribution
   end subroutine pcm_component_gradient

   !> Contract the explicit position, width and switching Jacobians
   !>
   !> Other surface channels require derivatives outside the base cavity
   !> interface and are rejected here instead of falling back to reverse mode
   !>
   !> @param[in,out] cavity Live cavity, builds forward derivatives on demand
   !> @param[in] acc Surface adjoint (ngrid)
   !> @param[in,out] gradient Nuclear-gradient accumulator (3, nsph)
   !> @param[out] error Allocated on missing or inconsistent derivatives
   subroutine pcm_forward_surface_gradient(cavity, acc, gradient, error)
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Surface adjoint
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Surface point, atom and displacement axis
      integer :: igrid, iatom, iaxis

      if (any(acc%w_a /= 0.0_wp) .or. any(acc%w_w /= 0.0_wp) .or. &
          any(acc%w_n /= 0.0_wp) .or. any(acc%w_k1 /= 0.0_wp) .or. any(acc%w_k2 /= 0.0_wp)) then
         call fatal_error(error, "PCM forward gradient supports position, width and switching channels only")
         return
      end if
      if (any(shape(gradient) /= [3, cavity%nsph])) then
         call fatal_error(error, "PCM forward gradient shape mismatch")
         return
      end if
      if (.not. allocated(cavity%xi1_rA) .or. .not. allocated(cavity%f1_rA) .or. &
          .not. allocated(cavity%xyz1_rA)) then
         call cavity%get_gradient(error)
         if (allocated(error)) return
      end if
      if (.not. allocated(cavity%xi1_rA) .or. .not. allocated(cavity%f1_rA) .or. &
          .not. allocated(cavity%xyz1_rA)) then
         call fatal_error(error, "Cavity does not provide PCM forward nuclear derivatives")
         return
      end if
      if (any(shape(cavity%xi1_rA) /= [3, cavity%nsph, cavity%ngrid]) .or. &
          any(shape(cavity%f1_rA) /= [3, cavity%nsph, cavity%ngrid]) .or. &
          any(shape(cavity%xyz1_rA) /= [3, 3, cavity%nsph, cavity%ngrid])) then
         call fatal_error(error, "PCM forward surface-derivative shape mismatch")
         return
      end if

      do igrid = 1, cavity%ngrid
         gradient = gradient + acc%w_xi(igrid)*cavity%xi1_rA(:, :, igrid) &
            & + acc%w_f(igrid)*cavity%f1_rA(:, :, igrid)
         do iatom = 1, cavity%nsph
            do iaxis = 1, 3
               gradient(iaxis, iatom) = gradient(iaxis, iatom) &
                  & + dot_product(acc%w_xyz(:, igrid), cavity%xyz1_rA(:, iaxis, iatom, igrid))
            end do
         end do
      end do
   end subroutine pcm_forward_surface_gradient

   !> Nuclear charges of the solute
   !>
   !> The moving sources of the host potential are the solute nuclei, which
   !> enter only the direct nuclear gradient at fixed surface, weighted by the
   !> potential adjoint: every
   !> host quantity is a total, so nothing else in the component sees them
   !>
   !> The host's charges `z_eff` (atomic numbers minus ECP core electrons) when
   !> the structure carries them, else the atomic numbers
   !>
   !> @param[in]  self PCM component instance
   !> @param[out] za   Nuclear charges (nat)
   subroutine pcm_component_nuclear_charges(self, za)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(in) :: self
      !> Nuclear charges
      real(wp), allocatable, intent(out) :: za(:)

      !> Number of solute atoms and atom index
      integer :: nat, iatom

      nat = self%mol_solu%nat
      allocate (za(nat))
      if (allocated(self%mol_solu%z_eff)) then
         if (size(self%mol_solu%z_eff) == nat) then
            za(:) = real(self%mol_solu%z_eff, wp)
            return
         end if
      end if
      do iatom = 1, nat
         za(iatom) = real(self%mol_solu%num(self%mol_solu%id(iatom)), wp)
      end do

   end subroutine pcm_component_nuclear_charges

   !> Contract the raw spatial derivative with this component's potential adjoint
   !>
   !> The nuclear and electronic contributions share the host's potential convention
   !>
   !> @param[in]  coupling QM coupling data
   !> @param[in]  w_phi    Potential adjoint of this component (ngrid)
   !> @param[out] w_xyz    Host total surface-position weight (3, ngrid)
   !> @param[out] error    Error handling
   subroutine read_host_position_weight(coupling, w_phi, w_xyz, error)
      !> Scoped raw host primitives
      class(coupling_view_type), intent(in) :: coupling
      !> Potential adjoint of this component
      real(wp), intent(in) :: w_phi(:)
      !> Contracted position weights
      real(wp), allocatable, intent(out) :: w_xyz(:, :)
      !> Missing spatial derivative
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      call coupling%read("potential", "dphi_dr", w_xyz, error)
      if (allocated(error)) return
      do i = 1, size(w_phi)
         w_xyz(:, i) = w_phi(i)*w_xyz(:, i)
      end do
   end subroutine read_host_position_weight

   !> Add the operator and host trace surface weights for the gradient
   subroutine pcm_component_get_gradient_surface_weights(self, coupling, cavity, acc, error)
      !> PCM component with its own solved charges
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Raw host primitives in this component's scope
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error

      call self%get_surface_weights(coupling, cavity, acc, error)
      if (allocated(error)) return
      if (.not. cavity%has_field_dependent_geometry()) then
         call self%get_host_surface_weights(coupling, acc, cavity%ngrid, error)
      end if
   end subroutine pcm_component_get_gradient_surface_weights

   !> Nuclear gradient of the PCM electrostatics at fixed surface
   !>
   !> The solute nuclei move under the fixed surface, weighted by the potential
   !> adjoint (the surface charges of a stationary PCM); this term does
   !> not reach the energy through any cavity surface quantity, so it stays
   !> with the component instead of going through the cavity contraction
   !>
   !> @param[in,out] self     PCM component instance
   !> @param[in]    coupling QM coupling data
   !> @param[in,out] cavity   Live cavity
   !> @param[in,out] gradient Nuclear-gradient accumulator
   !> @param[out]   error    Error handling
   subroutine pcm_component_get_direct_gradient(self, coupling, cavity, gradient, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> QM coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity
      class(cavity_type), intent(inout) :: cavity
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Direct term and the nuclear charges
      real(wp), allocatable :: grad_direct(:, :), grad_normal(:, :), za(:)
      !> Potential adjoint of this component
      real(wp), allocatable :: w_phi(:)
      real(wp), allocatable :: w_vn(:), w_grad(:, :)
      integer :: i
      !> Solute atom count and grid size
      integer :: nat, ngrid

      nat = self%mol_solu%nat
      ngrid = cavity%ngrid
      if (size(gradient, 1) /= 3 .or. size(gradient, 2) /= nat) then
         call fatal_error(error, "[pcm_component_get_direct_gradient] gradient shape mismatch")
         return
      end if

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      if (self%is_inert()) return

      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call self%nuclear_charges(za)
      allocate (grad_direct(3, nat))
      call pcm_electrostatic_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
         & w_phi, za, grad_direct, error, xi=cavity%xi0)
      if (allocated(error)) return

      if (self%uses_normal_trace()) then
         if (.not. allocated(cavity%normal0)) then
            call fatal_error(error, "PCM normal trace requires cavity normals")
            return
         end if
         call self%normal_trace_adjoint(w_vn, error)
         if (allocated(error)) return
         if (size(w_vn) /= ngrid) then
            call fatal_error(error, "PCM normal adjoint size mismatch")
            return
         end if
         allocate (w_grad(3, ngrid), grad_normal(3, nat))
         do i = 1, ngrid
            w_grad(:, i) = w_vn(i)*cavity%normal0(:, i)
         end do
         call pcm_electrostatic_normal_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
            & w_grad, za, grad_normal, error, xi=cavity%xi0)
         if (allocated(error)) return
         grad_direct = grad_direct + grad_normal
      end if
      gradient = gradient + grad_direct

   end subroutine pcm_component_get_direct_gradient

   !> Report a concrete PCM implementation without surface derivatives
   subroutine pcm_component_get_surface_weights(self, coupling, cavity, acc, error)
      class(model_continuum_component_pcm), intent(inout) :: self
      class(coupling_view_type), intent(in) :: coupling
      class(cavity_type), intent(in) :: cavity
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "PCM component does not implement surface weights")
   end subroutine pcm_component_get_surface_weights

   !> Contract raw potential derivatives into position and Gaussian-width adjoints
   !>
   !> @param[in,out] self     PCM component instance
   !> @param[in]    coupling QM coupling data carrying the host weights
   !> @param[in,out] acc      Accumulated cavity surface adjoints
   !> @param[in]    ngrid    Expected electrostatic grid size
   !> @param[out]   error    Error handling
   subroutine pcm_component_get_host_surface_weights(self, coupling, acc, ngrid, error)
      !> PCM component with its own solved charges
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Raw host primitives
      class(coupling_view_type), intent(in) :: coupling
      !> Model-total surface adjoints
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Expected point count
      integer, intent(in) :: ngrid
      !> Missing or mismatched derivative
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xyz(:, :), w_xi(:)
      !> Potential adjoint of this component
      real(wp), allocatable :: w_phi(:)
      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call read_host_position_weight(coupling, w_phi, w_xyz, error)
      if (allocated(error)) return
      if (size(w_xyz, 2) /= ngrid) then
         call fatal_error(error, "PCM spatial derivative grid mismatch")
         return
      end if
      call coupling%read("potential", "dphi_dxi", w_xi, error)
      if (allocated(error)) return
      call acc%add_surface_weights(error, w_xyz=w_xyz, w_xi=w_phi*w_xi)
   end subroutine pcm_component_get_host_surface_weights

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Declare Gaussian potential outputs independently for each phase
   !>
   !> Spatial and width derivatives are needed in the response only when the
   !> level set follows the density, and in every nuclear gradient
   !>
   !> @param[in]    self     PCM component instance
   !> @param[in]    cavity   Cavity the model is built on
   !> @param[in,out] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine pcm_component_declare_coupling(self, cavity, coupling, error)
      !> PCM component
      class(model_continuum_component_pcm), intent(in) :: self
      !> Live cavity used to select phase requirements
      class(cavity_type), intent(in) :: cavity
      !> Component-scoped registration context
      type(coupling_type), intent(inout) :: coupling
      !> Registration error
      type(error_type), allocatable, intent(out) :: error
      type(gaussian_potential_request_type) :: potential
      logical :: need_normal

      need_normal = self%uses_normal_trace()
      call request_require(potential, moist_phase_energy, "phi", error)
      if (allocated(error)) return
      if (need_normal) then
         call request_require(potential, moist_phase_energy, "dphi_dr", error)
         if (allocated(error)) return
      end if
      call request_require(potential, moist_phase_response, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "dphi_dr", &
         & cavity%has_field_dependent_geometry() .or. need_normal, error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "dphi_dxi", cavity%has_field_dependent_geometry(), error)
      if (allocated(error)) return
      if (need_normal) then
         call request_require(potential, moist_phase_response, "d2phi_dr2", &
            & cavity%has_field_dependent_geometry(), error)
         if (allocated(error)) return
         call request_require(potential, moist_phase_response, "d2phi_dr_dxi", &
            & cavity%has_field_dependent_geometry(), error)
         if (allocated(error)) return
      end if
      call request_require(potential, moist_phase_gradient, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "dphi_dr", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "dphi_dxi", error)
      if (allocated(error)) return
      if (need_normal) then
         call request_require(potential, moist_phase_gradient, "d2phi_dr2", error)
         if (allocated(error)) return
         call request_require(potential, moist_phase_gradient, "d2phi_dr_dxi", error)
         if (allocated(error)) return
      end if
      call coupling_register(coupling, "potential", potential, error)
   end subroutine pcm_component_declare_coupling

end module moist_model_continuum_component_pcm_type
