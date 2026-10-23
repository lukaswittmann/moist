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
      & parameter_adjoint_response_type, response_accumulate
   use moist_channels_coupling, only: coupling_type, coupling_view_type, coupling_request_type, &
      & point_potential_request_type, gaussian_potential_request_type, nuclear_charge_request_type, &
      & moist_phase_energy, moist_phase_response, moist_phase_gradient, coupling_register, request_require
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_continuum_component_pcm_electrostatics, only: &
      & pcm_electrostatic_direct_gradient, &
      & pcm_electrostatic_normal_direct_gradient
   use moist_utils_timer, only: cat_energy
   use moist_model_continuum_component_pcm_solvers, only: pcm_storage_dense, &
      & pcm_storage_matrix_free
   use moist_utils_prettyprint, only: prettyprinter
   use moist_channels_fields, only: field_query_type
   implicit none(type, external)
   private

   public :: model_continuum_component_pcm
   public :: pcm_solver_type
   public :: solver_type
   public :: pcm_tangent_block, pcm_tangent_require
   public :: pcm_component_declare_coupling
   public :: pcm_component_list_fields

   !> One-block kernel direction of a surface channel, zero when unallocated
   interface pcm_tangent_block
      module procedure pcm_tangent_block_scalar
      module procedure pcm_tangent_block_vector
   end interface pcm_tangent_block

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
      integer :: solver = solver_type%iterative
      !> Iterative solver convergence tolerance
      real(wp) :: solver_tol = 1.0e-8_wp
      !> Maximum iterative solver steps
      integer :: solver_maxiter = 1000
      !> Krylov dimension before GMRES restarts
      integer :: solver_restart = 100
      !> Operator storage choice for iterative solves
      integer :: storage = pcm_storage_dense
      !> Relative nonlinear residual threshold
      real(wp) :: newton_tol = 1.0e-10_wp
      !> Maximum Newton steps
      integer :: newton_maxiter = 50
      !> Add the parameter derivatives of the energy to `get_response`, never to gradient reads
      logical :: parameter_derivatives = .false.
      !> Couple Gaussian surfaces through point potential integrals; no effect on point surfaces
      logical :: point_potential = .false.
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
      !> Check solver controls
      procedure :: validate => validate_pcm_parameters
      !> Print the solver and the operator storage by name
      procedure :: print_choices => print_pcm_choices
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

      !> Energy returned by the generated solve
      !>
      !> Residual-corrected for inexact solves, consistent with the gradient weights
      real(wp) :: kernel_energy = 0.0_wp

      !> Numerical settings shared by computation and parameter introspection
      type(moist_pcm_parameters_type) :: param

      !> Molecular electrostatic potential at cavity grid points (ngrid)
      !>
      !> Copied from the coupling's potential request by `ensure_charges`;
      !> never computed here
      real(wp), allocatable :: phi(:)

      !> Whether charges belong to the current cavity and host potential
      logical :: charges_valid = .false.

      !> Whether the host answers Gaussian potential integrals, set by `update`
      !>
      !> Read by the host-weight routines, which get no cavity; the operator
      !> discretization is the separate concern of the concrete models
      logical :: gaussian_potential = .false.

   contains

      !> Request normals when the model uses a normal potential trace
      procedure :: configure_cavity => pcm_component_configure_cavity

      !> Whether this component reads n dot grad(phi)
      procedure :: uses_normal_trace => pcm_component_uses_normal_trace

      !> Whether the host potential uses Gaussian integrals on a cavity
      procedure :: uses_gaussian_potential => pcm_component_uses_gaussian_potential

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

      !> Derivatives of the energy with respect to the physical parameters
      procedure :: get_parameter_adjoint => pcm_component_get_parameter_adjoint

      !> Compute the potential adjoint `dE/dphi` and, for normal-trace models, `dE/dvn`
      procedure :: get_trace_response => pcm_component_get_trace_response

      !> Compute PCM gradient with respect to nuclear coordinates
      procedure :: get_gradient => pcm_component_get_gradient

      !> Contract the explicit nuclear Jacobians of Gaussian surface geometry
      procedure :: get_forward_gradient => pcm_component_get_forward_gradient

      !> Accumulate the PCM cavity surface adjoint weights
      procedure :: get_surface_weights => pcm_component_get_surface_weights

      !> Accumulate the surface adjoints the reverse-mode gradient consumes
      procedure :: get_gradient_surface_weights => pcm_component_get_gradient_surface_weights

      !> Directional derivative of the operator surface weights and trace adjoints
      procedure :: surface_weight_tangent => pcm_component_surface_weight_tangent

      !> Nuclear gradient contributions that bypass the cavity surface
      procedure :: get_direct_gradient => pcm_component_get_direct_gradient

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

      !> Publish the potential adjoint of the current charges
      procedure :: list_fields => pcm_component_list_fields

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

   !> Whether the host potential uses Gaussian integrals on a cavity
   !>
   !> - Gaussian surface: Gaussian integrals unless `point_potential` is set
   !> - Point surface: always point integrals
   !>
   !> @param[in] self   PCM component
   !> @param[in] cavity Cavity the component is built on
   pure function pcm_component_uses_gaussian_potential(self, cavity) result(uses)
      !> PCM component
      class(model_continuum_component_pcm), intent(in) :: self
      !> Cavity the component is built on
      class(cavity_type), intent(in) :: cavity
      !> Gaussian potential integrals
      logical :: uses

      uses = allocated(cavity%xi0) .and. .not. self%param%point_potential
   end function pcm_component_uses_gaussian_potential

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

   !> Form the normal trace `vn = n dot grad(phi)` from the host's raw spatial derivative
   !>
   !> - Outward unit normal `n`, pointing into the solvent
   !> - Host potential only, nuclear plus electronic, without the reaction field
   !> - `vn = -E_n`, minus the outward normal component of the solute electric field
   !> - Gaussian surfaces: derivative of the smeared potential, so the smeared field
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
          (self%storage /= pcm_storage_dense .and. &
           self%storage /= pcm_storage_matrix_free) .or. &
          self%newton_tol <= 0.0_wp .or. &
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
      self%newton_tol = defaults%newton_tol
      self%newton_maxiter = defaults%newton_maxiter
      self%parameter_derivatives = defaults%parameter_derivatives
      self%point_potential = defaults%point_potential
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
      call self%register_real_scalar("newton_tol", self%newton_tol)
      call self%register_int_scalar("newton_maxiter", self%newton_maxiter)
      call self%register_logical("parameter_derivatives", self%parameter_derivatives)
      call self%register_logical("point_potential", self%point_potential)
   end subroutine register_parameter_entries

   !> Print the linear solver and, for the iterative one, the operator storage by name
   !>
   !> @param[in] self Parameter values
   !> @param[in,out] pp Pretty printer inside the component section
   subroutine print_pcm_choices(self, pp)
      !> Parameter values
      class(moist_pcm_parameters_type), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp
      !> Name of the current choice
      character(len=:), allocatable :: label

      if (self%solver == solver_type%dense) then
         label = "Dense"
      else if (self%solver == solver_type%iterative) then
         label = "Iterative (GMRES)"
      else
         label = "unknown"
      end if
      call pp%kv("Solver", self%solver, label)
      if (self%solver /= solver_type%iterative) return

      if (self%storage == pcm_storage_dense) then
         label = "Dense"
      else if (self%storage == pcm_storage_matrix_free) then
         label = "Matrix-free"
      else
         label = "unknown"
      end if
      call pp%kv("Storage", self%storage, label)
   end subroutine print_pcm_choices

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

   !> Print the dielectric constant
   !>
   !> @param[in]    self PCM component instance
   !> @param[inout] pp   Pretty printer inside the component section
   subroutine pcm_component_print_inputs(self, pp)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp

      call pp%kv("Epsilon", self%epsilon)

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
   !> this vector, never the charges themselves. It equals the surface charge
   !> `q` only for stationary potential-only models (CPCM, COSMO); IEF-PCM,
   !> SS(V)PE and LPB return the symmetrized adjoint, and models that also consume
   !> the normal trace (DPCM, DPCM-original, SLIC, nonlocal) return their own
   !> `dE/dphi`, with `dE/dvn` carried separately by `normal_trace_adjoint`
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

   !> Publish `potential_adjoint`, this component's `w_phi = dE/dphi` (ngrid)
   !>
   !> The share of the model's potential adjoint response item; not listed
   !> while the charges are not current
   !>
   !> @param[in]    self  PCM component
   !> @param[inout] query Field walker
   subroutine pcm_component_list_fields(self, query)
      !> PCM component
      class(model_continuum_component_pcm), intent(in) :: self
      !> Field walker
      type(field_query_type), intent(inout) :: query

      !> Potential adjoint of the current charges
      real(wp), allocatable :: w_phi(:)
      !> Charges not current
      type(error_type), allocatable :: error

      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call query%add_real("potential_adjoint", &
         & "Potential adjoint dE/dphi of this component, e (ngrid)", w_phi)

   end subroutine pcm_component_list_fields

   !> Compute PCM solvation energy
   !>
   !> Returns the energy of the generated solve, e.g. 0.5*dot(q, phi) for an
   !> exact conductor solve; residual-corrected for iterative solves
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

      energy = energy + self%kernel_energy

      call self%ctx%timer%stop("PCM energy")

   end subroutine pcm_component_get_energy

   !> Compute the PCM potential adjoint
   !>
   !> Returns the `potential_adjoint_response_type` item with `w_phi = dE/dphi`
   !> and, for models that consume the normal trace `vn`, the normal-trace
   !> adjoint `dE/dvn` as `w_dphi_dr`; the host contracts `w_phi` with its own
   !> potential integrals, `F_uv += sum_i w_phi_i V_uv(r_i)`, and `w_dphi_dr`
   !> with their position derivatives. `w_phi` equals the surface charge `q`
   !> only for stationary potential-only models (CPCM, COSMO)
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

      !> Parameter derivative item of this component
      type(parameter_adjoint_response_type) :: pitem

      call self%require_context(error)
      if (allocated(error)) return

      call self%get_trace_response(coupling, cavity, response, error)
      if (allocated(error)) return
      ! Parameter item on the response phase only; the reverse gradient path
      ! reuses this routine and must neither carry nor compute it
      if (self%param%parameter_derivatives .and. coupling%phase == moist_phase_response) then
         call self%get_parameter_adjoint(pitem, error)
         if (allocated(error)) return
         call response_accumulate(response, pitem, error)
      end if

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

   !> Derivatives of the energy with respect to the physical parameters
   !>
   !> Overridden by concrete models; runs only when the parameter
   !> `parameter_derivatives` is set and only in `get_response` (response phase), on the
   !> solution of the current traces
   !>
   !> @param[in,out] self  PCM component with current charges
   !> @param[out]    item  parameter derivatives
   !> @param[out]    error error handling
   subroutine pcm_component_get_parameter_adjoint(self, item, error)
      !> PCM component instance
      class(model_continuum_component_pcm), intent(inout) :: self
      !> Parameter derivatives
      type(parameter_adjoint_response_type), intent(out) :: item
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "Parameter derivatives are not available for "//trim(self%name))

   end subroutine pcm_component_get_parameter_adjoint

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

   !> Nuclear charges of the solute, answered by the host
   !>
   !> The moving sources of the host potential are the solute nuclei, which
   !> enter only the direct nuclear gradient at fixed surface, weighted by the
   !> potential adjoint: every host quantity is a total, so nothing else in the
   !> component sees them
   !>
   !> The charges of the host Hamiltonian, `nuclear_charges` `z`: the atomic
   !> numbers minus the core electrons an ECP replaces
   !>
   !> @param[in]  coupling Scoped raw host primitives
   !> @param[in]  nat      Number of solute atoms
   !> @param[out] za       Nuclear charges (nat)
   !> @param[out] error    Missing answer or atom-count mismatch
   subroutine read_host_nuclear_charges(coupling, nat, za, error)
      !> Scoped raw host primitives
      class(coupling_view_type), intent(in) :: coupling
      !> Number of solute atoms
      integer, intent(in) :: nat
      !> Nuclear charges
      real(wp), allocatable, intent(out) :: za(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call coupling%read("nuclear_charges", "z", za, error)
      if (allocated(error)) return
      if (size(za) /= nat) then
         call fatal_error(error, "PCM nuclear charges do not match the solute atom count")
      end if

   end subroutine read_host_nuclear_charges

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
   !> adjoint `w_phi` (equal to `q` only for stationary potential-only models);
   !> this term does not reach the energy through any cavity surface quantity,
   !> so it stays with the component instead of going through the cavity
   !> contraction. The nuclear potential matches the host request, Gaussian or
   !> point
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
      !> Gaussian host potential
      logical :: gaussian

      nat = self%mol_solu%nat
      ngrid = cavity%ngrid
      if (size(gradient, 1) /= 3 .or. size(gradient, 2) /= nat) then
         call fatal_error(error, "[pcm_component_get_direct_gradient] gradient shape mismatch")
         return
      end if

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return

      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call read_host_nuclear_charges(coupling, nat, za, error)
      if (allocated(error)) return
      gaussian = self%uses_gaussian_potential(cavity)
      allocate (grad_direct(3, nat))
      if (gaussian) then
         call pcm_electrostatic_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
            & w_phi, za, grad_direct, error, xi=cavity%xi0)
      else
         call pcm_electrostatic_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
            & w_phi, za, grad_direct, error)
      end if
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
         if (gaussian) then
            call pcm_electrostatic_normal_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
               & w_grad, za, grad_normal, error, xi=cavity%xi0)
         else
            call pcm_electrostatic_normal_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
               & w_grad, za, grad_normal, error)
         end if
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

   !> Report a concrete PCM implementation without a surface weight tangent
   !>
   !> Contract of the overrides, first piece of the nuclear Hessian:
   !> - Derivative along (surface fields, dphi, dvn) of the operator weights
   !>   `get_surface_weights` adds on a field-independent cavity (no host terms)
   !>   and of `potential_adjoint` and `normal_trace_adjoint`
   !> - Needs charges solved on the current traces
   !> - Direction channels as in the accumulator (w_xyz, w_n, w_k1, w_k2, w_xi,
   !>   w_f, w_a); an unallocated channel is a zero direction, a channel the
   !>   model does not read is ignored
   !>
   !> @param[in,out] self      PCM component with current charges
   !> @param[in]     cavity    Cavity the operators were assembled on
   !> @param[in]     direction Surface direction, one entry per grid point and channel
   !> @param[in]     dphi      Potential-trace direction (ngrid)
   !> @param[in,out] tangent   Initialized accumulator, receives the weight tangent
   !> @param[out]    dw_phi    Tangent of `potential_adjoint` (ngrid)
   !> @param[out]    error     Allocated on failure
   !> @param[in]     dvn       Normal-trace direction (ngrid); absent is zero
   !> @param[out]    dw_vn     Tangent of `normal_trace_adjoint` (ngrid), normal-trace models only
   subroutine pcm_component_surface_weight_tangent(self, cavity, direction, dphi, tangent, &
         & dw_phi, error, dvn, dw_vn)
      !> PCM component with current charges
      class(model_continuum_component_pcm), intent(inout) :: self
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

      call fatal_error(error, "PCM component does not implement the surface weight tangent")
   end subroutine pcm_component_surface_weight_tangent

   !> Check the preconditions of a surface weight tangent
   !>
   !> @param[in]  self   PCM component
   !> @param[in]  cavity Cavity the operators were assembled on
   !> @param[in]  dphi   Potential-trace direction (ngrid)
   !> @param[out] error  Allocated on stale charges or a grid mismatch
   !> @param[in]  dvn    Optional normal-trace direction (ngrid)
   subroutine pcm_tangent_require(self, cavity, dphi, error, dvn)
      !> PCM component
      class(model_continuum_component_pcm), intent(in) :: self
      !> Cavity the operators were assembled on
      class(cavity_type), intent(in) :: cavity
      !> Potential-trace direction
      real(wp), intent(in) :: dphi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Normal-trace direction
      real(wp), intent(in), optional :: dvn(:)

      if (.not. self%charges_valid .or. .not. allocated(self%q)) then
         call fatal_error(error, "PCM tangent needs charges solved on the current traces")
         return
      end if
      if (cavity%ngrid /= size(self%q)) then
         call fatal_error(error, "PCM tangent cavity grid does not match the operators")
         return
      end if
      if (size(dphi) /= cavity%ngrid) then
         call fatal_error(error, "PCM tangent potential direction does not match the cavity grid")
         return
      end if
      if (present(dvn)) then
         if (size(dvn) /= cavity%ngrid) then
            call fatal_error(error, "PCM tangent normal-trace direction does not match the cavity grid")
         end if
      end if
   end subroutine pcm_tangent_require

   !> One-block kernel direction of a scalar surface channel
   !>
   !> @param[in]  channel Direction channel (ngrid), unallocated for zero
   !> @param[in]  ngrid   Number of grid points
   !> @param[in]  name    Channel name for the error message
   !> @param[out] block   Kernel direction block (ngrid, 1)
   !> @param[out] error   Allocated on a size mismatch
   subroutine pcm_tangent_block_scalar(channel, ngrid, name, block, error)
      !> Direction channel
      real(wp), allocatable, intent(in) :: channel(:)
      !> Number of grid points
      integer, intent(in) :: ngrid
      !> Channel name
      character(len=*), intent(in) :: name
      !> Kernel direction block
      real(wp), allocatable, intent(out) :: block(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      allocate (block(ngrid, 1), source=0.0_wp)
      if (.not. allocated(channel)) return
      if (size(channel) /= ngrid) then
         call fatal_error(error, "PCM tangent direction "//name//" does not match the cavity grid")
         return
      end if
      block(:, 1) = channel
   end subroutine pcm_tangent_block_scalar

   !> One-block kernel direction of a vector surface channel
   !>
   !> @param[in]  channel Direction channel (3, ngrid), unallocated for zero
   !> @param[in]  ngrid   Number of grid points
   !> @param[in]  name    Channel name for the error message
   !> @param[out] block   Kernel direction block (3, ngrid, 1)
   !> @param[out] error   Allocated on a shape mismatch
   subroutine pcm_tangent_block_vector(channel, ngrid, name, block, error)
      !> Direction channel
      real(wp), allocatable, intent(in) :: channel(:, :)
      !> Number of grid points
      integer, intent(in) :: ngrid
      !> Channel name
      character(len=*), intent(in) :: name
      !> Kernel direction block
      real(wp), allocatable, intent(out) :: block(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      allocate (block(3, ngrid, 1), source=0.0_wp)
      if (.not. allocated(channel)) return
      if (any(shape(channel) /= [3, ngrid])) then
         call fatal_error(error, "PCM tangent direction "//name//" does not match the cavity grid")
         return
      end if
      block(:, :, 1) = channel
   end subroutine pcm_tangent_block_vector

   !> Contract raw potential derivatives into position and Gaussian-width adjoints
   !>
   !> The width adjoint only for Gaussian potential integrals; point integrals
   !> do not depend on the widths
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
      !> Raw host derivatives, contracted in place
      real(wp), allocatable :: w_xyz(:, :), w_xi(:)
      !> Potential adjoint of this component
      real(wp), allocatable :: w_phi(:)
      integer :: i

      call coupling%read("potential", "dphi_dr", w_xyz, error)
      if (allocated(error)) return
      if (any(shape(w_xyz) /= [3, ngrid])) then
         call fatal_error(error, "PCM spatial derivative grid mismatch")
         return
      end if
      if (self%gaussian_potential) then
         call coupling%read("potential", "dphi_dxi", w_xi, error)
         if (allocated(error)) return
         if (size(w_xi) /= ngrid) then
            call fatal_error(error, "PCM width derivative grid mismatch")
            return
         end if
      end if
      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      if (size(w_phi) /= ngrid) then
         call fatal_error(error, "PCM spatial derivative grid mismatch")
         return
      end if
      do i = 1, ngrid
         w_xyz(:, i) = w_phi(i)*w_xyz(:, i)
      end do
      if (self%gaussian_potential) then
         call acc%add_surface_weights(error, w_xyz=w_xyz, w_xi=w_phi*w_xi)
      else
         call acc%add_surface_weights(error, w_xyz=w_xyz)
      end if
   end subroutine pcm_component_get_host_surface_weights

   !* ================================================================================= *!
   !*                            Coupling declaration and staging                       *!
   !* ================================================================================= *!

   !> Declare the host potential outputs independently for each phase
   !>
   !> - `gaussian_potential` when `uses_gaussian_potential`, else `point_potential`
   !> - Spatial derivatives (and width derivatives of Gaussian integrals) in the
   !>   response only when the level set follows the density, and in every
   !>   nuclear gradient; second derivatives only for normal-trace models
   !> - The nuclear charges only in the gradient
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
      class(coupling_request_type), allocatable :: potential
      type(nuclear_charge_request_type) :: charges
      logical :: normal, field, gaussian

      normal = self%uses_normal_trace()
      field = cavity%has_field_dependent_geometry()
      gaussian = self%uses_gaussian_potential(cavity)
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
      call request_require(potential, moist_phase_response, "d2phi_dr2", normal .and. field, error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "dphi_dr", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "d2phi_dr2", normal, error)
      if (allocated(error)) return
      if (gaussian) then
         call request_require(potential, moist_phase_response, "dphi_dxi", field, error)
         if (allocated(error)) return
         call request_require(potential, moist_phase_response, "d2phi_dr_dxi", normal .and. field, error)
         if (allocated(error)) return
         call request_require(potential, moist_phase_gradient, "dphi_dxi", error)
         if (allocated(error)) return
         call request_require(potential, moist_phase_gradient, "d2phi_dr_dxi", normal, error)
         if (allocated(error)) return
      end if
      call coupling_register(coupling, "potential", potential, error)
      if (allocated(error)) return
      call request_require(charges, moist_phase_gradient, "z", error)
      if (allocated(error)) return
      call coupling_register(coupling, "nuclear_charges", charges, error)
   end subroutine pcm_component_declare_coupling

end module moist_model_continuum_component_pcm_type
