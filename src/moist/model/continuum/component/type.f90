!> Typed cavity component interface and standalone coupling lifecycle
module moist_model_continuum_component_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type
   use mctc_io_codata2018, only: atomic_unit_of_mass, atomic_unit_of_time, Bohr_radius
   use moist_context, only: moist_context_type
   use moist_cavity_type, only: cavity_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_channels_response, only: response_type
   use moist_channels_coupling, only: coupling_type, coupling_view_type, moist_phase_energy, &
      & moist_phase_response, moist_phase_gradient, coupling_begin_registration, &
      & coupling_set_scope, coupling_snapshot, coupling_arm, coupling_invalidate

   implicit none(type, external)
   private

   public :: model_continuum_component_type
   public :: autogpa

   !> Pressure Eh/bohr^3 -> GPa: the atomic unit of pressure m_e/(a_0 t_au^2) in Pa, times 1e-9
   real(wp), parameter :: autogpa = atomic_unit_of_mass/(Bohr_radius*atomic_unit_of_time**2)*1.0e-9_wp

   !> Abstract solvation model component
   type, abstract :: model_continuum_component_type
      !> Borrowed run context (verbosity/debug/timer)
      !>
      !> Set at construction; owned by the caller, never allocated or freed here
      type(moist_context_type), pointer :: ctx => null()
      !> Name of the component
      character(len=:), allocatable :: name
      !> One-line description of the component, set by its constructor
      character(len=:), allocatable :: description
      !> Molecular structure data for the component
      type(structure_type) :: mol_solu
      !> Linear scale factor applied to this contribution
      !>
      !> - the component multiplies its energy, solvation response and
      !>   surface/level set response by this constant, so the contribution
      !>   stays variational
      !> - 1.0 leaves it unchanged, 0.0 disables it
      real(wp) :: scale = 1.0_wp
      !> Error handling
      type(error_type), allocatable :: error
   contains

      procedure(update_component), deferred :: update
      procedure(get_component_energy), deferred :: get_energy
      procedure(get_component_response), deferred :: get_response
      procedure(get_component_gradient), deferred :: get_gradient
      !> Internal emission hook: accumulate the response items that depend only
      !> on the potential adjoint; called by a component's own `get_response` and
      !> `get_gradient`, while the public phase accessors are `get_energy`,
      !> `get_response` and `get_gradient`
      procedure :: get_trace_response => get_component_trace_response_default
      !> Accumulate component-specific surface adjoint weights
      procedure :: get_surface_weights => get_component_surface_weights_default
      !> Accumulate the host's direct trace-geometry surface adjoint weights
      procedure :: get_host_surface_weights => get_component_host_surface_weights_default
      !> Accumulate the surface adjoint weights the *nuclear gradient* needs
      procedure :: get_gradient_surface_weights => get_component_gradient_surface_weights_default
      !> Accumulate nuclear-gradient terms that do not flow through the surface
      procedure :: get_direct_gradient => get_component_direct_gradient_default
      !> Declare the host requests this component reads (none by default)
      procedure :: declare_coupling => declare_component_coupling_default
      !> Build the coupling of a bare component driven without a model
      procedure :: new_coupling => new_component_coupling
      !> Stage one phase of a bare component's coupling
      procedure, private :: stage => stage_component_coupling
      !> Stage the energy phase
      procedure :: prepare_energy => prepare_component_energy
      !> Stage the response phase
      procedure :: prepare_response => prepare_component_response
      !> Stage the gradient phase
      procedure :: prepare_gradient => prepare_component_gradient
      !> Solvent-independent settings; null for a component without any
      procedure :: parameters => component_parameters_default
      !> Print the solvent and physical inputs (none by default)
      procedure :: print_inputs => print_component_inputs_default
      !> Print name, description, scale, inputs and settings
      procedure :: print_parameters => print_component_parameters

   end type model_continuum_component_type

   !> Deferred procedure contracts
   abstract interface

      !> Update the solvation model component with the current molecular structure
      !>
      !> @param[in,out] self Instance
      !> @param[in] mol Molecular structure
      !> @param[in,out] cavity Live cavity
      !> @param[out] error Error handling
      subroutine update_component(self, mol, cavity, error)
         import model_continuum_component_type, structure_type, cavity_type, error_type
         implicit none(type, external)
         !> Instance of the solvation model component
         class(model_continuum_component_type), intent(inout) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Cavity the model is built on
         class(cavity_type), intent(inout) :: cavity
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_component

      !> Evaluate the solvation energy for the component
      !>
      !> @param[in,out] self Instance
      !> @param[in] coupling Host coupling
      !> @param[in,out] cavity Live cavity
      !> @param[in,out] energy Energy accumulator
      !> @param[out] error Error handling
      subroutine get_component_energy(self, coupling, cavity, energy, error)
         import model_continuum_component_type, cavity_type, wp, coupling_view_type, error_type
         implicit none(type, external)
         !> Instance of the solvation model component
         class(model_continuum_component_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_view_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> solvation energy for the component
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_energy

      !> Get the solvation response for the component
      !>
      !> @param[in,out] self Instance
      !> @param[in] coupling Host coupling
      !> @param[in,out] cavity Live cavity
      !> @param[in,out] response Host response accumulator
      !> @param[out] error Error handling
      subroutine get_component_response(self, coupling, cavity, response, error)
         import model_continuum_component_type, cavity_type, response_type, coupling_view_type, error_type
         implicit none(type, external)
         !> Instance of the solvation model component
         class(model_continuum_component_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_view_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> Solvation response for the component
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_response

      !> Get the solvation energy gradient for the component
      !>
      !> Accumulates into `response` the host part of the gradient phase (what
      !> the host contracts with its own geometry derivatives)
      !>
      !> @param[in,out] self Instance
      !> @param[in] coupling Host coupling
      !> @param[in,out] cavity Live cavity
      !> @param[in,out] response Host response accumulator
      !> @param[in,out] gradient Nuclear gradient accumulator
      !> @param[out] error Error handling
      subroutine get_component_gradient(self, coupling, cavity, response, gradient, error)
         import model_continuum_component_type, cavity_type, wp, response_type, coupling_view_type, &
            & error_type
         implicit none(type, external)
         !> Instance of the solvation model component
         class(model_continuum_component_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_view_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> Host part of the gradient phase for the component
         type(response_type), intent(inout) :: response
         !> Solvation gradient for the component
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_gradient

   end interface

contains

   !* ================================================================================= *!
   !*                                 Component hooks                                 *!
   !* ================================================================================= *!

   !> Default no-op direct trace-response hook
   !>
   !> @param[in,out] self      Solvation component
   !> @param[in]    coupling  Host coupling data
   !> @param[in,out] cavity    Live model cavity
   !> @param[in,out] response  Direct trace-response accumulator
   !> @param[out]   error     Error object
   subroutine get_component_trace_response_default(self, coupling, cavity, response, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Host coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Direct trace-response accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_trace_response_default

   !> Default no-op surface-weight hook for components without cavity response
   !>
   !> @param[in,out] self    Solvation component
   !> @param[in]    coupling     Wavefunction data
   !> @param[in]    cavity  Cavity data
   !> @param[in,out] acc     Cavity-specific surface-adjoint accumulator
   !> @param[out]   error   Error object
   subroutine get_component_surface_weights_default(self, coupling, cavity, acc, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Cavity data
      class(cavity_type), intent(in) :: cavity
      !> Surface-adjoint accumulator (a surface accumulator on a cavity)
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_surface_weights_default

   !> Default no-op host trace-geometry weight hook
   !>
   !> Surface traces built from the host's QM integrals (potential, normal
   !> derivative, ...) carry a surface dependence moist cannot differentiate, so
   !> the host supplies dE/d(xi, f, r, n) at fixed operator through the
   !> surface-weight requests of the coupling
   !>
   !> Components with such a trace override this hook to add those channels to
   !> the shared surface-adjoint accumulator; the rest inherit the no-op
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    coupling Wavefunction data carrying the host weights
   !> @param[in,out] acc      Surface-adjoint accumulator
   !> @param[in]    ngrid    Expected grid size of the component's cavity
   !> @param[out]   error    Error object
   subroutine get_component_host_surface_weights_default(self, coupling, acc, ngrid, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Surface-adjoint accumulator (a surface accumulator on a cavity)
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Expected grid size of the component's cavity
      integer, intent(in) :: ngrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_host_surface_weights_default

   !> Default gradient-side surface weights: the same ones the response uses
   !>
   !> For most components the surface adjoint of the energy is one object, so
   !> the reverse-mode nuclear gradient can reuse `get_surface_weights`
   !> verbatim; a component whose gradient legitimately consumes a different
   !> set of host channels overrides this (see `model_continuum_component_pcm`)
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    coupling Wavefunction data
   !> @param[in]    cavity   Cavity data
   !> @param[in,out] acc      Cavity-specific surface-adjoint accumulator
   !> @param[out]   error    Error object
   subroutine get_component_gradient_surface_weights_default(self, coupling, cavity, acc, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Cavity data
      class(cavity_type), intent(in) :: cavity
      !> Surface-adjoint accumulator (a surface accumulator on a cavity)
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%get_surface_weights(coupling, cavity, acc, error)

   end subroutine get_component_gradient_surface_weights_default

   !> Default no-op hook for nuclear-gradient terms outside the surface
   !>
   !> Used by the reverse-mode gradient path for contributions that do not
   !> reach the energy through a cavity surface quantity -- for PCM, the
   !> solute nuclei moving under fixed surface charges
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    coupling Wavefunction data
   !> @param[in,out] cavity   Cavity data
   !> @param[in,out] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error object
   subroutine get_component_direct_gradient_default(self, coupling, cavity, gradient, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Cavity data
      class(cavity_type), intent(inout) :: cavity
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_direct_gradient_default

   !> Default component-side declaration: the component reads no host request
   !>
   !> @param[in]    self     Solvation component
   !> @param[in]    cavity   Cavity the model is built on
   !> @param[in,out] coupling Coupling being declared
   !> @param[out]   error    Error handling
   subroutine declare_component_coupling_default(self, cavity, coupling, error)
      !> Solvation component
      class(model_continuum_component_type), intent(in) :: self
      !> Cavity the model is built on
      class(cavity_type), intent(in) :: cavity
      !> Coupling being declared
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine declare_component_coupling_default

   !* ================================================================================= *!
   !*                             Bare-component coupling                             *!
   !* ================================================================================= *!

   !> Declare the cavity's and the component's requests, then record the grid size
   !>
   !> @param[in]    self     Solvation component
   !> @param[in]    cavity   Updated cavity
   !> @param[in,out] coupling Coupling to declare
   !> @param[out]   error    Error handling
   subroutine declare_component_pass(self, cavity, coupling, error)
      !> Solvation component
      class(model_continuum_component_type), intent(in) :: self
      !> Updated cavity
      class(cavity_type), intent(in) :: cavity
      !> Coupling to declare
      type(coupling_type), intent(inout) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call coupling_begin_registration(coupling)
      call coupling_set_scope(coupling, 0)
      call cavity%declare_coupling(coupling, error)
      if (allocated(error)) return
      call coupling_set_scope(coupling, 1)
      call self%declare_coupling(cavity, coupling, error)
      if (allocated(error)) return
      call coupling_snapshot(coupling, ngrid=cavity%ngrid, natom=cavity%nsph)

   end subroutine declare_component_pass

   !> Build the coupling of a bare component driven without a model
   !>
   !> The request list never grows after this call; `prepare_*` re-declares
   !> the same requests with the phase's requirements
   !>
   !> @param[in]  self     Solvation component
   !> @param[in]  cavity   Updated cavity the component was updated on
   !> @param[out] coupling Coupling to build
   !> @param[out] error    Error handling
   subroutine new_component_coupling(self, cavity, coupling, error)
      !> Solvation component
      class(model_continuum_component_type), intent(in) :: self
      !> Updated cavity
      class(cavity_type), intent(in) :: cavity
      !> Coupling to build
      type(coupling_type), intent(out) :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call declare_component_pass(self, cavity, coupling, error)

   end subroutine new_component_coupling

   !> Stage one phase of a bare component's coupling
   !>
   !> - energy starts a new host evaluation; later phases reuse valid raw
   !>   answers
   !> - the coupling keeps no grid copy, so after changing the cavity the
   !>   caller must stage the energy phase before any later one
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    cavity   Live cavity
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[in]    phase    Phase index, `moist_phase_energy` and so on
   !> @param[out]   error    Invalid phase or failed declaration
   subroutine stage_component_coupling(self, cavity, coupling, phase, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Live cavity
      class(cavity_type), intent(in) :: cavity
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Phase index
      integer, intent(in) :: phase
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (phase == moist_phase_energy) call coupling_invalidate(coupling)
      call declare_component_pass(self, cavity, coupling, error)
      if (allocated(error)) return
      call coupling_arm(coupling, phase, error)

   end subroutine stage_component_coupling

   !> Stage the energy phase of a bare component's coupling
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    cavity   Live cavity
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[out]   error    Error handling
   subroutine prepare_component_energy(self, cavity, coupling, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Live cavity
      class(cavity_type), intent(in) :: cavity
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(cavity, coupling, moist_phase_energy, error)

   end subroutine prepare_component_energy

   !> Stage the response phase of a bare component's coupling
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    cavity   Live cavity
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[out]   error    Error handling
   subroutine prepare_component_response(self, cavity, coupling, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Live cavity
      class(cavity_type), intent(in) :: cavity
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(cavity, coupling, moist_phase_response, error)

   end subroutine prepare_component_response

   !> Stage the gradient phase of a bare component's coupling
   !>
   !> @param[in,out] self     Solvation component
   !> @param[in]    cavity   Live cavity
   !> @param[in,out] coupling Coupling built by `new_coupling`
   !> @param[out]   error    Error handling
   subroutine prepare_component_gradient(self, cavity, coupling, error)
      !> Solvation component
      class(model_continuum_component_type), intent(inout) :: self
      !> Live cavity
      class(cavity_type), intent(in) :: cavity
      !> Coupling built by `new_coupling`
      type(coupling_type), intent(inout), target :: coupling
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%stage(cavity, coupling, moist_phase_gradient, error)

   end subroutine prepare_component_gradient

   !* ================================================================================= *!
   !*                                    Printout                                     *!
   !* ================================================================================= *!

   !> Default settings hook: the component has no parameter object
   !>
   !> @param[in] self Solvation component
   function component_parameters_default(self) result(param)
      !> Solvation component
      class(model_continuum_component_type), intent(in), target :: self
      !> Settings of the component, never associated here
      class(moist_model_parameters_type), pointer :: param

      param => null()

   end function component_parameters_default

   !> Default no-op input printer for a component without solvent inputs
   !>
   !> @param[in]    self Solvation component
   !> @param[inout] pp   Pretty printer inside the component section
   subroutine print_component_inputs_default(self, pp)
      !> Solvation component
      class(model_continuum_component_type), intent(in) :: self
      !> Pretty printer inside the component section
      type(prettyprinter), intent(inout) :: pp

   end subroutine print_component_inputs_default

   !> Print one component section
   !>
   !> - the header is the name, or `Component <index> (<name>)` inside a model
   !> - then the description, `scale` when it is not 1, the solvent inputs and
   !>   the registered settings of `parameters()`
   !>
   !> @param[in] self  Solvation component
   !> @param[in] unit  Output unit; defaults to the run context's unit
   !> @param[in] index 1-based position in a model, shown in the header
   subroutine print_component_parameters(self, unit, index)
      !> Solvation component
      class(model_continuum_component_type), intent(in), target :: self
      !> Output unit
      integer, intent(in), optional :: unit
      !> 1-based position in a model
      integer, intent(in), optional :: index

      !> Section printer
      type(prettyprinter) :: pp
      !> Registered settings, when the component has any
      class(moist_model_parameters_type), pointer :: param
      !> Component name and section header
      character(len=:), allocatable :: name, header
      !> Index as text
      character(len=12) :: label
      !> Effective output unit
      integer :: iu

      iu = output_unit
      if (associated(self%ctx)) iu = self%ctx%unit
      if (present(unit)) iu = unit
      name = ""
      if (allocated(self%name)) name = self%name
      if (present(index)) then
         write (label, "(i0)") index
         header = "Component "//trim(label)
         if (len(name) > 0) header = header//" ("//name//")"
      else
         header = name
      end if

      pp = new_prettyprinter(unit=iu)
      call pp%push(header//":")
      if (allocated(self%description)) call pp%section(self%description)
      if (self%scale /= 1.0_wp) call pp%kv("scale", self%scale)
      call self%print_inputs(pp)
      param => self%parameters()
      if (associated(param)) call param%print_table(pp)
      call pp%pop()
      call pp%blank()

   end subroutine print_component_parameters

end module moist_model_continuum_component_type
