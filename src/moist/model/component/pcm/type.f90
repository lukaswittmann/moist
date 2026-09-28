!> PCM (Polarizable Continuum Model) abstract base type
!>
!> The abstract PCM base component, extended by the concrete implementations
!> (CPCM, COSMO, ...) in their respective modules
module moist_model_component_pcm_type
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env, only: wp, fatal_error
   use mctc_env_error, only: error_type
   use mctc_io, only: structure_type
   use moist_cavity_type, only: cavity_type
   use moist_model_type, only: solvation_model_component_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type, &
      & response_accumulate
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & gaussian_potential_request_type, moist_phase_energy, moist_phase_response, &
      & moist_phase_gradient, coupling_register, request_require
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
      & pcm_amat_surface_weights, pcm_amat_nuclear_gradient
   use moist_model_component_pcm_electrostatics, only: &
      & pcm_electrostatic_nuclear_gradient, pcm_electrostatic_direct_gradient
   use moist_utils_timer, only: cat_setup, cat_energy, cat_solve
   implicit none(type, external)
   private

   public :: solvation_model_component_pcm
   public :: pcm_solver_type
   public :: solver_type

   !> Enumerator for PCM solver types
   type :: pcm_solver_type
      !> Matrix inversion
      integer :: inversion = 1
      !> LU factorization (LAPACK GETRF+GETRS)
      integer :: lu = 2
      !> Cholesky factorization (SPD matrices)
      integer :: cholesky = 3
      !> Iterative (CG with preconditioner)
      integer :: iterative = 4
   end type pcm_solver_type

   !> Global instance of the solver enum
   type(pcm_solver_type), parameter :: solver_type = pcm_solver_type()

   public :: moist_pcm_parameters_type

   !> PCM solver configuration
   type, extends(moist_model_parameters_type) :: moist_pcm_parameters_type
      !> Linear solver selection
      integer :: solver = solver_type%cholesky
      !> Iterative solver convergence tolerance
      real(wp) :: solver_tol = 1.0e-10_wp
      !> Maximum iterative solver steps
      integer :: solver_maxiter = 50
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
   !> Matrix assembly uses the generic Gaussian surface kernel (assemble_pcm_amat)
   !>
   !> Wraps general moist solvers for the linear system solution
   type, abstract, extends(solvation_model_component_type) :: solvation_model_component_pcm

      !> Dielectric constant of the solvent
      real(wp) :: epsilon

      !> Dielectric scaling factor f( epsilon ) - variant-specific formula
      !>
      !> CPCM: f epsilon = ( epsilon -1)/ epsilon
      !> COSMO: f epsilon = ( epsilon -1)/( epsilon +0.5)
      real(wp) :: feps

      !> Surface charges on cavity grid points (ngrid)
      real(wp), allocatable :: q(:)

      !> PCM interaction matrix A (ngrid, ngrid)
      real(wp), allocatable :: amat(:, :)


      ! TODO: remove these parameters here and use the parameter type

      !> Solver type identifier
      integer :: solver = solver_type%lu

      !> Use external matrix (bypasses assembly if .true.)
      logical :: use_external_matrix = .false.

      !> Convergence tolerance for iterative solvers
      real(wp) :: solver_tol = 1.0e-10_wp

      !> Maximum iterations for iterative solvers
      integer :: solver_maxiter = 50

      !> Molecular electrostatic potential at cavity grid points (ngrid)
      !>
      !> Copied from the coupling's potential request by `ensure_charges`;
      !> never computed here
      real(wp), allocatable :: phi(:)

      !> Whether self%q holds the charges belonging to the current matrix and phi
      logical :: charges_valid = .false.

   contains

      !> Update PCM component: assembles matrix and prepares for charge solution
      procedure :: update => pcm_component_update

      !> Compute PCM solvation energy
      procedure :: get_energy => pcm_component_get_energy

      !> Compute PCM reaction potential
      procedure :: get_response => pcm_component_get_response

      !> Compute the potential adjoint `dE/dphi` (the surface charges for a stationary PCM)
      procedure :: get_trace_response => pcm_component_get_trace_response

      !> Compute PCM gradient with respect to nuclear coordinates
      procedure :: get_gradient => pcm_component_get_gradient

      !> Accumulate the PCM cavity surface adjoint weights
      procedure :: get_surface_weights => pcm_component_get_surface_weights

      !> Accumulate the surface adjoints the reverse-mode gradient consumes
      procedure :: get_gradient_surface_weights => pcm_component_get_gradient_surface_weights

      !> Nuclear gradient contributions that bypass the cavity surface
      procedure :: get_direct_gradient => pcm_component_get_direct_gradient

      !> Nuclear charges of the solute, the moving sources of the potential
      procedure :: nuclear_charges => pcm_component_nuclear_charges

      !> Contract the current PCM charges to Gaussian-surface weights
      procedure :: amat_surface_weights => pcm_component_amat_surface_weights

      !> Contract the current PCM charges to the A-matrix nuclear gradient
      procedure :: amat_nuclear_gradient => pcm_component_amat_nuclear_gradient

      !> Fold the host's direct trace-geometry weights into the accumulator
      procedure :: get_host_surface_weights => pcm_component_get_host_surface_weights

      !> Declare the potential and the host surface weights
      procedure :: declare_coupling => pcm_component_declare_coupling

      !> Solve for the surface charges unless they are already current
      procedure :: ensure_charges => pcm_ensure_charges

      !> Potential adjoint `w_phi = dE/dphi` of the current charges
      procedure :: potential_adjoint => pcm_component_potential_adjoint

      !> Set external matrix (bypasses internal assembly)
      procedure :: set_external_matrix => pcm_set_external_matrix

      !> Solve the PCM linear system A . q = rhs using selected solver
      procedure :: solve_system => pcm_solve_system

   end type solvation_model_component_pcm

contains

   !> Validate PCM solver settings
   !>
   !> @param[inout] self Solver settings
   !> @param[out] error Invalid setting
   subroutine validate_pcm_parameters(self, error)
      class(moist_pcm_parameters_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      if (self%solver < solver_type%inversion .or. self%solver > solver_type%iterative .or. &
          self%solver_tol <= 0.0_wp .or. self%solver_tol /= self%solver_tol .or. &
          self%solver_maxiter < 1) then
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
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_pcm_parameters_type), intent(inout), target :: self

      call self%register_int_scalar("solver", self%solver)
      call self%register_real_scalar("solver_tol", self%solver_tol)
      call self%register_int_scalar("solver_maxiter", self%solver_maxiter)
   end subroutine register_parameter_entries

   !> Update PCM base component
   !>
   !> Stores references to mol/cavity and assembles the PCM matrix
   !>
   !> (unless using an external matrix)
   !>
   !> @param[in,out] self PCM component instance
   !> @param[in] mol Molecular structure data
   !> @param[in,out] cavity Cavity type data
   !> @param[out] error Error handling
   subroutine pcm_component_update(self, mol, cavity, error)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(inout) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Cavity type data
      class(cavity_type), intent(inout) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: ngrid
      !> Timer depth on entry, restored on every early return
      integer :: d0

      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("PCM setup", category=cat_setup)

      ! Store references (the cavity is owned by the orchestrating model)
      self%mol_solu = mol

      ngrid = cavity%ngrid

      ! Drop the host potential so a new geometry always gets a fresh one,
      ! even if ngrid happens to be the same
      if (allocated(self%phi)) deallocate (self%phi)

      ! Allocate charge array
      if (allocated(self%q)) then
         if (size(self%q) /= ngrid) deallocate (self%q)
      end if
      if (.not. allocated(self%q)) then
         allocate (self%q(ngrid))
      end if
      self%q(:) = 0.0_wp

      ! Assemble or use external matrix
      if (.not. self%use_external_matrix) then
         if (allocated(self%amat)) then
            if (size(self%amat, 1) /= ngrid .or. size(self%amat, 2) /= ngrid) then
               deallocate (self%amat)
            end if
         end if
         if (.not. allocated(self%amat)) then
            allocate (self%amat(ngrid, ngrid))
         end if

         if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f)) then
            call fatal_error(error, &
               & "[pcm_component_update] Cavity does not provide a Gaussian PCM surface")
            call self%ctx%timer%unwind(d0)
            return
         end if

         ! Generic Gaussian surface-charge interaction matrix
         call self%ctx%timer%start("Interaction matrix")
         call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, self%amat, error)
         if (allocated(error)) then
            call self%ctx%timer%unwind(d0)
            return
         end if
         call self%ctx%timer%stop("Interaction matrix")
      else
         if (.not. allocated(self%amat)) then
            call fatal_error(error, &
               & "[pcm_component_update] External PCM matrix requested but not allocated")
            call self%ctx%timer%unwind(d0)
            return
         end if
         if (size(self%amat, 1) /= ngrid .or. size(self%amat, 2) /= ngrid) then
            call fatal_error(error, &
               & "[pcm_component_update] External PCM matrix dimension mismatch")
            call self%ctx%timer%unwind(d0)
            return
         end if
      end if

      ! Note: charge solving happens on demand in ensure_charges, once the
      ! wavefunction data (electrostatic potential phi ) is available; a new
      ! geometry means a new matrix, so any cached charges are stale
      self%charges_valid = .false.

      call self%ctx%timer%stop("PCM setup")

   end subroutine pcm_component_update

   !> Solve for the induced surface charges unless they are already current
   !>
   !> The molecular potential is read off the coupling's potential request
   !>
   !> (mandatory: a stale answer is reported by name, never read as zero)
   !> A host may change the potential without touching the geometry (the
   !> ordinary SCF pattern), so update() has had no chance to clear the cache:
   !> cached charges stay valid only while phi is unchanged, which is a
   !> component-internal solve cache and deliberately not the request's
   !> staleness flag
   !>
   !> @param[in,out] self     PCM component instance
   !> @param[in]    coupling Host data carrying the molecular potential trace
   !> @param[in]    cavity   Live cavity used to assemble the current matrix
   !> @param[out]   error    Error handling
   subroutine pcm_ensure_charges(self, coupling, cavity, error)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(inout) :: self
      !> Host data carrying the molecular potential trace
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity used to assemble the current matrix
      class(cavity_type), intent(in) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Number of cavity grid points
      integer :: ngrid
      !> Right-hand side of the PCM linear system
      real(wp), allocatable :: rhs(:)
      !> Potential answered by the host
      real(wp), allocatable :: phi(:)
      !> The potential request

      ! Check that update() has prepared the matrix data
      if (.not. allocated(self%amat)) then
         call fatal_error(error, &
            & "[pcm_ensure_charges] PCM matrix not allocated - call update() first")
         return
      end if

      ngrid = cavity%ngrid

      call coupling%read("potential", "phi", phi, error)
      if (allocated(error)) return
      if (size(phi) /= ngrid) then
         call fatal_error(error, &
            & "[pcm_ensure_charges] External potential size mismatch")
         return
      end if
      if (.not. allocated(self%phi)) then
         self%charges_valid = .false.
      else if (size(self%phi) /= ngrid) then
         self%charges_valid = .false.
      else if (any(self%phi /= phi)) then
         self%charges_valid = .false.
      end if
      self%phi = phi

      if (self%charges_valid) return

      ! Build RHS: b = -f(eps) * phi
      allocate (rhs(ngrid))
      rhs(:) = -self%feps*self%phi(:)

      ! Solve for charges: A*q = b
      call self%solve_system(self%amat, rhs, self%q, error)
      if (allocated(error)) return

      self%charges_valid = .true.

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
      class(solvation_model_component_pcm), intent(in) :: self
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
   !> Like every accessor, this starts by reporting a stale mandatory request
   !> of the armed phase (or a failure latched by a mis-shaped `set`) by name
   !>
   !> @param[in,out] self PCM component instance
   !> @param[in] coupling Host data carrying the molecular potential trace
   !> @param[in,out] cavity Live cavity owned by the orchestrating model
   !> @param[in,out] energy Solvation energy (inout to allow accumulation)
   !> @param[out] error Error handling
   subroutine pcm_component_get_energy(self, coupling, cavity, energy, error)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(inout) :: self
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

      call coupling%check_mandatory(coupling%phase, error)
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
      class(solvation_model_component_pcm), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Solvation response
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call coupling%check_mandatory(coupling%phase, error)
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
      class(solvation_model_component_pcm), intent(inout) :: self
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

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      call self%potential_adjoint(item%w_phi, error)
      if (allocated(error)) return
      call response_accumulate(response, item, error)

   end subroutine pcm_component_get_trace_response

   !> Compute the PCM energy gradient with respect to nuclear coordinates
   !>
   !> For `A q = -f(eps) phi` and `E = 1/2 q^T phi`, the stationary
   !> derivative is
   !>    dE/dR_A = q^T dphi/dR_A
   !>
   !>              + 1/(2 f(eps)) q^T (dA/dR_A) q
   !>
   !> The A-matrix contribution is obtained from [[pcm_component_amat_nuclear_gradient]]
   !>
   !> The potential adjoint the host contracts with its own `dphi/dR` is
   !> accumulated into `response` from the already solved charges
   !>
   !> @param[in,out] self     PCM component instance
   !> @param[in]    coupling Wavefunction and electrostatic coupling data
   !> @param[in,out] cavity   Live cavity owned by the orchestrating model
   !> @param[in,out] response Host part of the gradient phase (potential adjoint)
   !> @param[in,out] gradient Solvation gradient accumulator
   !> @param[out]   error    Error handling
   subroutine pcm_component_get_gradient(self, coupling, cavity, response, gradient, error)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(inout) :: self
      !> Wavefunction data
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Host part of the gradient phase
      type(response_type), intent(inout) :: response
      !> Solvation gradient (3, nat)
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Raw A-matrix, electrostatic, and host-width gradient contributions
      real(wp), allocatable :: grad_amat(:, :), grad_electrostatic(:, :), grad_width(:, :)
      !> Host total surface-position weight `dE_host/dr_i` at the cavity points
      real(wp), allocatable :: w_xyz(:, :)
      !> Nuclear charges, the moving sources of the potential
      real(wp), allocatable :: za(:)
      !> Number of solute atoms and cavity grid points
      integer :: nat, ngrid
      !> Grid point index
      integer :: igrid
      !> Raw inverse-length derivative contracted with this component's adjoint
      real(wp), allocatable :: w_xi(:)
      !> Potential adjoint of this component
      real(wp), allocatable :: w_phi(:)

      call coupling%check_mandatory(coupling%phase, error)
      if (allocated(error)) return

      nat = self%mol_solu%nat
      ngrid = cavity%ngrid
      if (size(gradient, 1) /= 3 .or. size(gradient, 2) /= nat) then
         call fatal_error(error, "[pcm_component_get_gradient] gradient shape mismatch")
         return
      end if
      if (cavity%nsph /= nat) then
         call fatal_error(error, &
            & "[pcm_component_get_gradient] cavity sphere count does not match the solute")
         return
      end if

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return

      ! The host contracts the solved charges with its own dA/dR; they are
      ! published even when they vanish identically below
      call self%get_trace_response(coupling, cavity, response, error)
      if (allocated(error)) return

      ! At eps == 1 the right-hand side and surface charges vanish, so the
      ! polarization energy and all of its derivatives are exactly zero
      if (self%feps == 0.0_wp) return

      if (.not. allocated(cavity%xi1_rA) .or. &
          .not. allocated(cavity%f1_rA) .or. &
          .not. allocated(cavity%xyz1_rA)) then
         call cavity%get_gradient(error)
         if (allocated(error)) return
      end if

      allocate (grad_amat(3, nat), grad_electrostatic(3, nat))
      call self%amat_nuclear_gradient(cavity, grad_amat, error)
      if (allocated(error)) return

      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call read_host_position_weight(coupling, w_phi, w_xyz, error)
      if (allocated(error)) return
      call self%nuclear_charges(za)

      call pcm_electrostatic_nuclear_gradient(cavity%xyz, &
         & self%mol_solu%xyz, cavity%xyz1_rA, w_phi, w_xyz, &
         & za, grad_electrostatic, error, xi=cavity%xi0)
      if (allocated(error)) return

      call coupling%read("potential", "dphi_dxi", w_xi, error)
      if (allocated(error)) return
      w_xi = w_phi*w_xi
      allocate (grad_width(3, nat), source=0.0_wp)
      do igrid = 1, ngrid
         grad_width = grad_width + w_xi(igrid)*cavity%xi1_rA(:, :, igrid)
      end do

      gradient = gradient + 0.5_wp*grad_amat/self%feps + &
         & grad_electrostatic + grad_width

   end subroutine pcm_component_get_gradient

   !> Nuclear charges of the solute
   !>
   !> The moving sources of the host potential are the solute nuclei, which
   !> enter only the direct nuclear gradient at fixed surface, weighted by the
   !> potential adjoint: every
   !> host quantity is a total, so nothing else in the component sees them
   !>
   !> @param[in]  self PCM component instance
   !> @param[out] za   Nuclear charges (nat)
   subroutine pcm_component_nuclear_charges(self, za)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(in) :: self
      !> Nuclear charges
      real(wp), allocatable, intent(out) :: za(:)

      !> Number of solute atoms and atom index
      integer :: nat, iatom

      nat = self%mol_solu%nat
      allocate (za(nat))
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

   !> Accumulate PCM matrix and raw-potential adjoints for the nuclear gradient
   !>
   !> Each component contracts its own charges once before the cavity reverse pass
   !>
   !> @param[in,out] self     PCM component instance
   !> @param[in]    coupling QM coupling data
   !> @param[in]    cavity   Cavity the PCM matrix was assembled on
   !> @param[in,out] acc      Accumulated cavity surface adjoints
   !> @param[out]   error    Error handling
   subroutine pcm_component_get_gradient_surface_weights(self, coupling, cavity, acc, error)
      !> PCM component with its own solved charges
      class(solvation_model_component_pcm), intent(inout) :: self
      !> Raw host primitives in this component's scope
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity
      class(cavity_type), intent(in) :: cavity
      !> Model-total surface adjoints
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :), host_xyz(:, :), host_xi(:)
      !> Potential adjoint of this component
      real(wp), allocatable :: w_phi(:)
      real(wp) :: prefactor
      integer :: ngrid
      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return
      if (self%feps == 0.0_wp) return
      ngrid = cavity%ngrid
      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      call self%amat_surface_weights(cavity, w_xi, w_f, w_xyz, error)
      if (allocated(error)) return
      prefactor = 0.5_wp/self%feps
      call acc%add_surface_weights(error, w_xi=prefactor*w_xi, w_f=prefactor*w_f, &
         & w_xyz=prefactor*w_xyz)
      if (allocated(error)) return
      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call read_host_position_weight(coupling, w_phi, host_xyz, error)
      if (allocated(error)) return
      call coupling%read("potential", "dphi_dxi", host_xi, error)
      if (allocated(error)) return
      call acc%add_surface_weights(error, w_xyz=host_xyz, w_xi=w_phi*host_xi)
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
      class(solvation_model_component_pcm), intent(inout) :: self
      !> QM coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Live cavity
      class(cavity_type), intent(inout) :: cavity
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Direct term and the nuclear charges
      real(wp), allocatable :: grad_direct(:, :), za(:)
      !> Potential adjoint of this component
      real(wp), allocatable :: w_phi(:)
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
      if (self%feps == 0.0_wp) return

      call self%potential_adjoint(w_phi, error)
      if (allocated(error)) return
      call self%nuclear_charges(za)
      allocate (grad_direct(3, nat))
      call pcm_electrostatic_direct_gradient(cavity%xyz, self%mol_solu%xyz, &
         & w_phi, za, grad_direct, error, xi=cavity%xi0)
      if (allocated(error)) return

      gradient = gradient + grad_direct

   end subroutine pcm_component_get_direct_gradient

   !> Contract the current PCM charges to Gaussian-surface A-matrix weights
   !>
   !> Computes the derivatives of `q^T A q` wrt. cavity quantities
   !>
   !> Surface charges must have been produced by [[pcm_ensure_charges]] for that update
   !>
   !> @param[in]  self    PCM component with current surface charges
   !> @param[in]  cavity  Live cavity carrying the Gaussian PCM surface
   !> @param[out] w_xi    Gaussian-width weights
   !> @param[out] w_f     Switching-factor weights
   !> @param[out] w_xyz   Grid point-position weights
   !> @param[out] error   Error handling
   subroutine pcm_component_amat_surface_weights(self, cavity, w_xi, w_f, w_xyz, error)
      !> PCM component with current surface charges
      class(solvation_model_component_pcm), intent(in) :: self
      !> Live cavity carrying the Gaussian PCM surface
      class(cavity_type), intent(in) :: cavity
      !> Gaussian-width weights
      real(wp), intent(out) :: w_xi(:)
      !> Switching-factor weights
      real(wp), intent(out) :: w_f(:)
      !> Grid point-position weights
      real(wp), intent(out) :: w_xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%charges_valid .or. .not. allocated(self%q)) then
         call fatal_error(error, "[pcm_component_amat_surface_weights] "// &
            & "surface charges are unavailable - call ensure_charges first")
         return
      end if
      if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) .or. &
          .not. allocated(cavity%xyz)) then
         call fatal_error(error, "[pcm_component_amat_surface_weights] "// &
            & "cavity does not provide a Gaussian PCM surface")
         return
      end if

      call pcm_amat_surface_weights(cavity%xi0, cavity%f, cavity%xyz, self%q, &
         & self%q, w_xi, w_f, w_xyz, error)

   end subroutine pcm_component_amat_surface_weights

   !> Contract the current PCM charges to the A-matrix nuclear gradient
   !>
   !> Computes `q^T (dA/dR_A) q` from the Gaussian-surface derivative arrays on
   !> the supplied live cavity -- the A-matrix contribution only, without
   !> dielectric scaling or the host electrostatic-potential contribution
   !>
   !> @param[in]  self     PCM component with current surface charges
   !> @param[in]  cavity   Live cavity carrying Gaussian-surface derivatives
   !> @param[out] grad_rA  A-matrix nuclear gradient
   !> @param[out] error    Error handling
   subroutine pcm_component_amat_nuclear_gradient(self, cavity, grad_rA, error)
      !> PCM component with current surface charges
      class(solvation_model_component_pcm), intent(in) :: self
      !> Live cavity carrying Gaussian-surface derivatives
      class(cavity_type), intent(in) :: cavity
      !> A-matrix nuclear gradient
      real(wp), intent(out) :: grad_rA(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Gaussian-width weights
      real(wp), allocatable :: w_xi(:)
      !> Switching-factor weights
      real(wp), allocatable :: w_f(:)
      !> Grid point-position weights
      real(wp), allocatable :: w_xyz(:, :)

      if (.not. allocated(cavity%xi1_rA) .or. .not. allocated(cavity%f1_rA) .or. &
          .not. allocated(cavity%xyz1_rA)) then
         call fatal_error(error, "[pcm_component_amat_nuclear_gradient] "// &
            & "Gaussian PCM derivatives are unavailable - compute the cavity gradient first")
         return
      end if

      allocate (w_xi(cavity%ngrid), w_f(cavity%ngrid), w_xyz(3, cavity%ngrid))
      call self%amat_surface_weights(cavity, w_xi, w_f, w_xyz, error)
      if (allocated(error)) return

      call pcm_amat_nuclear_gradient(cavity%xi1_rA, cavity%f1_rA, &
         & cavity%xyz1_rA, w_xi, w_f, w_xyz, grad_rA, error)

   end subroutine pcm_component_amat_nuclear_gradient

   !> Accumulate the PCM cavity surface adjoint weights
   !>
   !> With A q = -f(eps) phi the energy is
   !>
   !>    E = 1/2 q^T phi = -q^T A q/(2 f(eps))
   !>
   !> Differentiating at fixed host potential and eliminating dq/dp through the
   !> linear system gives
   !>    dE/dp = 1/(2 f(eps)) * q^T (dA/dp) q  +  q^T (dphi/dp)
   !>
   !> The caller contracts the finished accumulator once via
   !>
   !> `cavity%contract_surface_lsf_weights` to obtain the level-set adjoints
   !>
   !> @param[in,out] self     PCM component instance
   !> @param[in]    coupling QM coupling data
   !> @param[in]    cavity   Cavity the PCM matrix was assembled on
   !> @param[in,out] acc      Accumulated cavity surface adjoints
   !> @param[out]   error    Error handling
   subroutine pcm_component_get_surface_weights(self, coupling, cavity, acc, error)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(inout) :: self
      !> QM coupling data
      class(coupling_view_type), intent(in) :: coupling
      !> Cavity the PCM matrix was assembled on
      class(cavity_type), intent(in) :: cavity
      !> Accumulated cavity surface adjoints
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-grid point adjoints of q^T A q w.r.t. xi_i, f_i and r_i
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      !> Number of cavity grid points
      integer :: ngrid
      !> The 1/(2 f(eps)) response prefactor
      real(wp) :: prefactor

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      if (size(self%q) /= ngrid) then
         call fatal_error(error, "[pcm_component_get_surface_weights] "// &
            & "surface charges do not match the cavity grid")
         return
      end if
      if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f)) then
         call fatal_error(error, &
            & "[pcm_component_get_surface_weights] Cavity does not provide a Gaussian PCM surface")
         return
      end if

      ! eps == 1: q == 0, so both operator and host responses vanish
      if (self%feps == 0.0_wp) return

      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      call self%amat_surface_weights(cavity, w_xi, w_f, w_xyz, error)
      if (allocated(error)) return

      prefactor = 0.5_wp/self%feps
      call acc%add_surface_weights(error, w_xi=prefactor*w_xi, &
         & w_f=prefactor*w_f, w_xyz=prefactor*w_xyz)
      if (allocated(error)) return
      if (cavity%has_field_dependent_geometry()) then
         call self%get_host_surface_weights(coupling, acc, ngrid, error)
      end if

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
      class(solvation_model_component_pcm), intent(inout) :: self
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
      class(solvation_model_component_pcm), intent(in) :: self
      !> Live cavity used to select phase requirements
      class(cavity_type), intent(in) :: cavity
      !> Component-scoped registration context
      type(coupling_type), intent(inout) :: coupling
      !> Registration error
      type(error_type), allocatable, intent(out) :: error
      type(gaussian_potential_request_type) :: potential
      call request_require(potential, moist_phase_energy, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "dphi_dr", cavity%has_field_dependent_geometry(), error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_response, "dphi_dxi", cavity%has_field_dependent_geometry(), error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "phi", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "dphi_dr", error)
      if (allocated(error)) return
      call request_require(potential, moist_phase_gradient, "dphi_dxi", error)
      if (allocated(error)) return
      call coupling_register(coupling, "potential", potential, error)
   end subroutine pcm_component_declare_coupling

   !> Set external matrix (bypasses internal assembly)
   !>
   !> Allows user to provide a pre-computed PCM matrix
   !>
   !> @param[in,out] self PCM component instance
   !> @param[in] amat External matrix (ngrid, ngrid)
   subroutine pcm_set_external_matrix(self, amat)
      !> PCM component instance
      class(solvation_model_component_pcm), intent(inout) :: self
      !> External matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)

      self%use_external_matrix = .true.
      if (allocated(self%amat)) deallocate (self%amat)
      allocate (self%amat, source=amat)
      self%charges_valid = .false.

   end subroutine pcm_set_external_matrix

   !> Solve the PCM linear system A*q = rhs
   !>
   !> Dispatches to appropriate solver based on self%solver setting
   !>
   !> @param[in] self PCM component instance
   !> @param[in] amat System matrix (ngrid, ngrid)
   !> @param[in] rhs Right-hand side (ngrid)
   !> @param[out] q Solution vector - surface charges (ngrid)
   !> @param[out] error Error handling
   subroutine pcm_solve_system(self, amat, rhs, q, error)
      use moist_model_component_pcm_solvers, only: solve_pcm_lu, &
         & solve_pcm_cholesky, solve_pcm_iterative, solve_pcm_inversion
      !> PCM component instance
      class(solvation_model_component_pcm), intent(in) :: self
      !> System matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)
      !> Right-hand side (ngrid)
      real(wp), intent(in) :: rhs(:)
      !> Solution vector - surface charges (ngrid)
      real(wp), intent(out) :: q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Timer depth on entry, restored on every early return
      integer :: d0

      !> The timer lives on the borrowed context, so it is writable even though
      !>
      !> self is intent(in) -- only the pointer association would be fixed
      d0 = self%ctx%timer%current_depth()
      call self%ctx%timer%start("PCM solve", category=cat_solve)

      select case (self%solver)
      case (solver_type%lu)
         call solve_pcm_lu(amat, rhs, q, error, unit=self%ctx%unit)

      case (solver_type%cholesky)
         call solve_pcm_cholesky(amat, rhs, q, error)

      case (solver_type%iterative)
         call solve_pcm_iterative(amat, rhs, q, self%solver_tol, &
            & self%solver_maxiter, error)

      case (solver_type%inversion)
         call solve_pcm_inversion(amat, rhs, q, error)

      case default
         call fatal_error(error, "[pcm_solve_system] Unknown solver type")
         call self%ctx%timer%unwind(d0)
         return
      end select

      call self%ctx%timer%stop("PCM solve")

   end subroutine pcm_solve_system

end module moist_model_component_pcm_type
