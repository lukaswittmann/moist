!> PCM (Polarizable Continuum Model) abstract base type
!> This module defines the abstract PCM base component that is extended by
!> concrete implementations (CPCM, COSMO, ...) in their respective modules.
module moist_model_component_pcm_type
   use mctc_env, only: wp, fatal_error
   use mctc_env_error, only: error_type
   use mctc_io, only: structure_type
   use moist_type, only: solvation_model_component, cavity_type, &
      & potential_type, coupling_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
      & pcm_amat_surface_weights, pcm_amat_nuclear_gradient
   use moist_model_component_pcm_electrostatics, only: &
      & pcm_electrostatic_nuclear_gradient
   use moist_utils_timer, only: cat_setup, cat_energy, cat_solve
   implicit none (type, external)
   private

   public :: pcm_base
   public :: pcm_solver_type, pcm_potential_source
   public :: solver_type, potential_source

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

   !> Enumerator for electrostatic potential source
   type :: pcm_potential_source
      !> Compute from atomic point charges
      integer :: charges = 1
      !> Provided externally via input_potential (QM coupling)
      integer :: external = 2
   end type pcm_potential_source

   !> Global instances for solver/potential enums
   type(pcm_solver_type), parameter :: solver_type = pcm_solver_type()
   type(pcm_potential_source), parameter :: potential_source = pcm_potential_source()

   !> Abstract PCM base component
   !> Provides common infrastructure for PCM-family methods (CPCM, COSMO, IEF-PCM).
   !> Matrix assembly uses the generic Gaussian surface kernel (assemble_pcm_amat).
   !> Wraps general moist solvers for the linear system solution.
   type, abstract, extends(solvation_model_component) :: pcm_base

      !> Dielectric constant of the solvent
      real(wp) :: epsilon

      !> Dielectric scaling factor f( epsilon ) - variant-specific formula
      !> CPCM: f epsilon = ( epsilon -1)/ epsilon
      !> COSMO: f epsilon = ( epsilon -1)/( epsilon +0.5)
      real(wp) :: feps

      !> Surface charges on cavity grid points (ngrid)
      real(wp), allocatable :: q(:)

      !> PCM interaction matrix A (ngrid, ngrid)
      real(wp), allocatable :: amat(:, :)

      !> Solver type identifier
      integer :: solver = solver_type%lu

      !> Use external matrix (bypasses assembly if .true.)
      logical :: use_external_matrix = .false.

      !> Convergence tolerance for iterative solvers
      real(wp) :: solver_tol = 1.0e-10_wp

      !> Maximum iterations for iterative solvers
      integer :: solver_maxiter = 1000

      !> Electrostatic potential at cavity grid points (ngrid)
      !> Set externally via input_potential when using external potential source.
      real(wp), allocatable :: phi(:)

      !> Potential source strategy
      integer :: phi_source = potential_source%charges

      !> Whether self%q holds the charges belonging to the current matrix and phi
      logical :: charges_valid = .false.

   contains

      !> Update PCM component: assembles matrix and prepares for charge solution
      procedure :: update => pcm_base_update

      !> Compute PCM solvation energy
      procedure :: get_energy => pcm_base_get_energy

      !> Compute PCM reaction potential
      procedure :: get_potential => pcm_base_get_potential

      !> Compute the direct electrostatic trace adjoint (the surface charges)
      procedure :: get_trace_potential => pcm_base_get_trace_potential

      !> Compute PCM gradient with respect to nuclear coordinates
      procedure :: get_gradient => pcm_base_get_gradient

      !> Accumulate the PCM cavity surface adjoint weights
      procedure :: get_surface_weights => pcm_base_get_surface_weights

      !> Contract the current PCM charges to Gaussian-surface weights
      procedure :: amat_surface_weights => pcm_base_amat_surface_weights

      !> Contract the current PCM charges to the A-matrix nuclear gradient
      procedure :: amat_nuclear_gradient => pcm_base_amat_nuclear_gradient

      !> Fold the host's direct trace-geometry weights into the accumulator
      procedure :: get_host_surface_weights => pcm_base_get_host_surface_weights

      !> Solve for the surface charges unless they are already current
      procedure :: ensure_charges => pcm_ensure_charges

      !> Set external matrix (bypasses internal assembly)
      procedure :: set_external_matrix => pcm_set_external_matrix

      !> Solve the PCM linear system A . q = rhs using selected solver
      procedure :: solve_system => pcm_solve_system

      !> Set external electrostatic potential at cavity grid points
      procedure :: input_potential => pcm_input_potential

   end type pcm_base

contains

   !> Update PCM base component
   !> Stores references to mol/cavity and assembles the PCM matrix
   !> (unless using an external matrix).
   subroutine pcm_base_update(self, mol, cavity, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
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
               & "[pcm_base_update] Cavity does not provide a Gaussian PCM surface")
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
               & "[pcm_base_update] External PCM matrix requested but not allocated")
            call self%ctx%timer%unwind(d0)
            return
         end if
         if (size(self%amat, 1) /= ngrid .or. size(self%amat, 2) /= ngrid) then
            call fatal_error(error, &
               & "[pcm_base_update] External PCM matrix dimension mismatch")
            call self%ctx%timer%unwind(d0)
            return
         end if
      end if

      ! Note: charge solving happens on demand in ensure_charges, once the
      ! wavefunction data (electrostatic potential phi ) is available. A new
      ! geometry means a new matrix, so any cached charges are stale.
      self%charges_valid = .false.

      call self%ctx%timer%stop("PCM setup")

   end subroutine pcm_base_update

   !> Solve for the induced surface charges unless they are already current
   !>
   !> @param[inout] self     PCM component instance
   !> @param[in]    coupling Wavefunction data (used only for the internal phi)
   !> @param[in]    cavity   Live cavity used to assemble the current matrix
   !> @param[out]   error    Error handling
   subroutine pcm_ensure_charges(self, coupling, cavity, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data (used only when no external potential is set)
      class(coupling_type), intent(in) :: coupling
      !> Live cavity used to assemble the current matrix
      class(cavity_type), intent(in) :: cavity
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Number of cavity grid points
      integer :: ngrid
      !> Right-hand side of the PCM linear system
      real(wp), allocatable :: rhs(:)
      !> Freshly evaluated potential, compared against the cache before storing
      real(wp), allocatable :: phi_new(:)

      ! Check that update() has prepared the matrix data
      if (.not. allocated(self%amat)) then
         call fatal_error(error, &
            & "[pcm_ensure_charges] PCM matrix not allocated - call update() first")
         return
      end if

      ngrid = cavity%ngrid

      ! Obtain electrostatic potential at cavity grid points
      select case (self%phi_source)
      case (potential_source%charges)
         ! Compute internally from atomic point charges
         allocate (phi_new(ngrid))
         call self%ctx%timer%start("Molecular potential")
         call compute_point_charge_potential(self%mol_solu, cavity, coupling, phi_new)
         call self%ctx%timer%stop("Molecular potential")
         ! A host may change the atomic charges without touching the geometry
         ! (the ordinary SCF pattern), so update() has had no chance to clear
         ! the cache. Cached charges stay valid only while phi is unchanged.
         if (.not. allocated(self%phi)) then
            self%charges_valid = .false.
         else if (size(self%phi) /= ngrid) then
            self%charges_valid = .false.
         else if (any(self%phi /= phi_new)) then
            self%charges_valid = .false.
         end if
         call move_alloc(phi_new, self%phi)

      case (potential_source%external)
         ! Prefer the coupling trace. The explicit input_potential method remains
         ! available as a compatibility adapter for direct component users.
         if (allocated(coupling%elstat_umol)) then
            if (size(coupling%elstat_umol) /= ngrid) then
               call fatal_error(error, &
                  & "[pcm_ensure_charges] External potential size mismatch")
               return
            end if
            if (.not. allocated(self%phi)) then
               self%charges_valid = .false.
            else if (size(self%phi) /= ngrid) then
               self%charges_valid = .false.
            else if (any(self%phi /= coupling%elstat_umol)) then
               self%charges_valid = .false.
            end if
            self%phi = coupling%elstat_umol
         end if
         if (.not. allocated(self%phi)) then
            call fatal_error(error, &
               & "[pcm_ensure_charges] External potential source selected "// &
               & "but phi not supplied - call input_potential() first")
            return
         end if

      case default
         call fatal_error(error, "[pcm_ensure_charges] Unknown potential source")
         return
      end select

      if (self%charges_valid) return

      ! Build RHS: b = -f(eps) * phi
      allocate (rhs(ngrid))
      rhs(:) = -self%feps*self%phi(:)

      ! Solve for charges: A*q = b
      call self%solve_system(self%amat, rhs, self%q, error)
      if (allocated(error)) return

      self%charges_valid = .true.

   end subroutine pcm_ensure_charges

   !> Compute PCM solvation energy
   !> E_solv = 0.5 * dot(q, phi)
   !> If an external potential was provided via input_potential, uses that;
   !> otherwise computes the potential internally from atomic point charges.
   subroutine pcm_base_get_energy(self, coupling, cavity, energy, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data (used only when no external potential is set)
      class(coupling_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Solvation energy (inout to allow accumulation)
      real(wp), intent(inout) :: energy
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: e_pcm
      !> Timer depth on entry, restored on every early return
      integer :: d0

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

   end subroutine pcm_base_get_energy

   !> Compute the PCM reaction potential channel
   !>
   !>    dE/dphi_i = q_i
   !>
   !> which is returned in `potential%w_elstat_umol`; the host contracts it with
   !> its own potential integrals to build the Fock contribution
   !> F_uv += sum_i q_i V_uv(r_i)
   subroutine pcm_base_get_potential(self, coupling, cavity, potential, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Solvation potential
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%get_trace_potential(coupling, cavity, potential, error)

   end subroutine pcm_base_get_potential

   !> Accumulate the direct electrostatic trace adjoint
   !>
   !> @param[inout] self      PCM component
   !> @param[in]    coupling  Host coupling data
   !> @param[inout] cavity    Live model cavity
   !> @param[inout] potential Direct trace-potential accumulator
   !> @param[out]   error     Error handling
   subroutine pcm_base_get_trace_potential(self, coupling, cavity, potential, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Direct trace-potential accumulator
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Number of cavity grid points
      integer :: ngrid

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return

      ngrid = cavity%ngrid

      if (allocated(potential%w_elstat_umol)) then
         if (size(potential%w_elstat_umol) /= ngrid) deallocate (potential%w_elstat_umol)
      end if
      if (.not. allocated(potential%w_elstat_umol)) then
         allocate (potential%w_elstat_umol(ngrid), source=0.0_wp)
      end if

      potential%w_elstat_umol(:) = potential%w_elstat_umol(:) + self%q(:)

   end subroutine pcm_base_get_trace_potential

   !> Compute the PCM energy gradient with respect to nuclear coordinates
   !>
   !> For `A q = -f(eps) phi` and `E = 1/2 q^T phi`, the stationary
   !> derivative is
   !>
   !>    dE/dR_A = q^T dphi/dR_A
   !>              + 1/(2 f(eps)) q^T (dA/dR_A) q.
   !>
   !> The A-matrix contribution is obtained from [[pcm_base_amat_nuclear_gradient]]
   !>
   !> @param[inout] self     PCM component instance
   !> @param[in]    coupling Wavefunction and electrostatic coupling data
   !> @param[inout] cavity   Live cavity owned by the orchestrating model
   !> @param[inout] gradient Solvation gradient accumulator
   !> @param[out]   error    Error handling
   subroutine pcm_base_get_gradient(self, coupling, cavity, gradient, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Live cavity owned by the orchestrating model
      class(cavity_type), intent(inout) :: cavity
      !> Solvation gradient (3, nat)
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Raw A-matrix, electrostatic, and host-width gradient contributions
      real(wp), allocatable :: grad_amat(:, :), grad_electrostatic(:, :), grad_width(:, :)
      !> Charge-weighted electronic field at the cavity points
      real(wp), allocatable :: qefield(:, :)
      !> Moving point charges entering the molecular potential
      real(wp), allocatable :: source_charge(:)
      !> Number of solute atoms and cavity grid points
      integer :: nat, ngrid
      !> Atom and grid-point indices
      integer :: iatom, igrid

      nat = self%mol_solu%nat
      ngrid = cavity%ngrid
      if (size(gradient, 1) /= 3 .or. size(gradient, 2) /= nat) then
         call fatal_error(error, "[pcm_base_get_gradient] gradient shape mismatch")
         return
      end if
      if (cavity%nsph /= nat) then
         call fatal_error(error, &
            & "[pcm_base_get_gradient] cavity sphere count does not match the solute")
         return
      end if

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return

      ! At eps == 1 the right-hand side and surface charges vanish, so the
      ! polarization energy and all of its derivatives are exactly zero.
      if (self%feps == 0.0_wp) return

      if (.not. allocated(cavity%xi1_rA) .or. &
          .not. allocated(cavity%f1_rA) .or. &
          .not. allocated(cavity%xyz1_rA)) then
         call cavity%get_gradient()
      end if
      if (allocated(cavity%error)) then
         allocate (error, source=cavity%error)
         return
      end if

      allocate (grad_amat(3, nat), grad_electrostatic(3, nat))
      call self%amat_nuclear_gradient(cavity, grad_amat, error)
      if (allocated(error)) return

      allocate (qefield(3, ngrid), source=0.0_wp)
      allocate (source_charge(nat))
      select case (self%phi_source)
      case (potential_source%charges)
         if (.not. allocated(coupling%qat)) then
            call fatal_error(error, &
               & "[pcm_base_get_gradient] internal potential requires atomic charges")
            return
         end if
         if (size(coupling%qat, 1) /= nat .or. size(coupling%qat, 2) < 1) then
            call fatal_error(error, &
               & "[pcm_base_get_gradient] atomic charge shape mismatch")
            return
         end if
         do iatom = 1, nat
            source_charge(iatom) = sum(coupling%qat(iatom, :))
         end do

      case (potential_source%external)
         do iatom = 1, nat
            source_charge(iatom) = &
               & real(self%mol_solu%num(self%mol_solu%id(iatom)), wp)
         end do
         if (allocated(coupling%elstat_qefield)) then
            if (size(coupling%elstat_qefield, 1) /= 3 .or. &
                size(coupling%elstat_qefield, 2) /= ngrid) then
               call fatal_error(error, &
                  & "[pcm_base_get_gradient] electronic field shape mismatch")
               return
            end if
            qefield = coupling%elstat_qefield
         end if

      case default
         call fatal_error(error, "[pcm_base_get_gradient] unknown potential source")
         return
      end select

      call pcm_electrostatic_nuclear_gradient(cavity%xyz, &
         & self%mol_solu%xyz, cavity%xyz1_rA, self%q, qefield, &
         & source_charge, grad_electrostatic, error)
      if (allocated(error)) return

      allocate (grad_width(3, nat), source=0.0_wp)
      if (allocated(coupling%elstat_w_xi)) then
         if (size(coupling%elstat_w_xi) /= ngrid) then
            call fatal_error(error, &
               & "[pcm_base_get_gradient] host Gaussian-width weights have wrong size")
            return
         end if
         do igrid = 1, ngrid
            grad_width = grad_width + &
               & coupling%elstat_w_xi(igrid)*cavity%xi1_rA(:, :, igrid)
         end do
      end if

      gradient = gradient + 0.5_wp*grad_amat/self%feps + &
         & grad_electrostatic + grad_width

   end subroutine pcm_base_get_gradient

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
   !> @param[out] w_xyz   Tessera-position weights
   !> @param[out] error   Error handling
   subroutine pcm_base_amat_surface_weights(self, cavity, w_xi, w_f, w_xyz, error)
      !> PCM component with current surface charges
      class(pcm_base), intent(in) :: self
      !> Live cavity carrying the Gaussian PCM surface
      class(cavity_type), intent(in) :: cavity
      !> Gaussian-width weights
      real(wp), intent(out) :: w_xi(:)
      !> Switching-factor weights
      real(wp), intent(out) :: w_f(:)
      !> Tessera-position weights
      real(wp), intent(out) :: w_xyz(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      if (.not. self%charges_valid .or. .not. allocated(self%q)) then
         call fatal_error(error, "[pcm_base_amat_surface_weights] "// &
            & "surface charges are unavailable - call ensure_charges first")
         return
      end if
      if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f) .or. &
          .not. allocated(cavity%xyz)) then
         call fatal_error(error, "[pcm_base_amat_surface_weights] "// &
            & "cavity does not provide a Gaussian PCM surface")
         return
      end if

      call pcm_amat_surface_weights(cavity%xi0, cavity%f, cavity%xyz, self%q, &
         & self%q, w_xi, w_f, w_xyz, error)

   end subroutine pcm_base_amat_surface_weights

   !> Contract the current PCM charges to the A-matrix nuclear gradient.
   !>
   !> Computes `q^T (dA/dR_A) q` from the Gaussian-surface derivative arrays on
   !> the supplied live cavity. This is the A-matrix contribution only, without
   !> dielectric scaling or the host electrostatic-potential contribution.
   !>
   !> @param[in]  self     PCM component with current surface charges
   !> @param[in]  cavity   Live cavity carrying Gaussian-surface derivatives
   !> @param[out] grad_rA  A-matrix nuclear gradient
   !> @param[out] error    Error handling
   subroutine pcm_base_amat_nuclear_gradient(self, cavity, grad_rA, error)
      !> PCM component with current surface charges
      class(pcm_base), intent(in) :: self
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
      !> Tessera-position weights
      real(wp), allocatable :: w_xyz(:, :)

      if (.not. allocated(cavity%xi1_rA) .or. .not. allocated(cavity%f1_rA) .or. &
          .not. allocated(cavity%xyz1_rA)) then
         call fatal_error(error, "[pcm_base_amat_nuclear_gradient] "// &
            & "Gaussian PCM derivatives are unavailable - compute the cavity gradient first")
         return
      end if

      allocate (w_xi(cavity%ngrid), w_f(cavity%ngrid), w_xyz(3, cavity%ngrid))
      call self%amat_surface_weights(cavity, w_xi, w_f, w_xyz, error)
      if (allocated(error)) return

      call pcm_amat_nuclear_gradient(cavity%xi1_rA, cavity%f1_rA, &
         & cavity%xyz1_rA, w_xi, w_f, w_xyz, grad_rA, error)

   end subroutine pcm_base_amat_nuclear_gradient

   !> Accumulate the PCM cavity surface adjoint weights
   !>
   !> With A q = -f(eps) phi the energy is E = 1/2 q^T phi = -q^T A q/(2 f(eps)),
   !> so differentiating at fixed host potential and eliminating dq/dp through the
   !> linear system gives
   !>
   !>    dE/dp = 1/(2 f(eps)) * q^T (dA/dp) q  +  q^T (dphi/dp)
   !>
   !> The caller contracts the finished accumulator once via
   !> `cavity%contract_surface_lsf_weights` to obtain the level-set adjoints
   !>
   !> @param[inout] self     PCM component instance
   !> @param[in]    coupling QM coupling data
   !> @param[in]    cavity   Cavity the PCM matrix was assembled on
   !> @param[inout] acc      Accumulated cavity surface adjoints
   !> @param[out]   error    Error handling
   subroutine pcm_base_get_surface_weights(self, coupling, cavity, acc, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> QM coupling data
      class(coupling_type), intent(in) :: coupling
      !> Cavity the PCM matrix was assembled on
      class(cavity_type), intent(in) :: cavity
      !> Accumulated cavity surface adjoints
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Per-tessera adjoints of q^T A q w.r.t. xi_i, f_i and r_i
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      !> Number of cavity grid points
      integer :: ngrid
      !> The 1/(2 f(eps)) response prefactor
      real(wp) :: prefactor

      call self%ensure_charges(coupling, cavity, error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      if (size(self%q) /= ngrid) then
         call fatal_error(error, "[pcm_base_get_surface_weights] "// &
            & "surface charges do not match the cavity grid")
         return
      end if
      if (.not. allocated(cavity%xi0) .or. .not. allocated(cavity%f)) then
         call fatal_error(error, &
            & "[pcm_base_get_surface_weights] Cavity does not provide a Gaussian PCM surface")
         return
      end if

      ! eps == 1: q == 0, so both operator and host responses vanish.
      if (self%feps == 0.0_wp) return

      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      call self%amat_surface_weights(cavity, w_xi, w_f, w_xyz, error)
      if (allocated(error)) return

      prefactor = 0.5_wp/self%feps
      call acc%add_surface_weights(error, w_xi=prefactor*w_xi, &
         & w_f=prefactor*w_f, w_xyz=prefactor*w_xyz)
      if (allocated(error)) return
      call self%get_host_surface_weights(coupling, acc, ngrid, error)

   end subroutine pcm_base_get_surface_weights

   !> Fold the host's direct trace-geometry weights into the surface accumulator
   !>
   !> The PCM energy depends on the cavity twice: through the matrix A the
   !> component builds itself (handled in [[pcm_base_get_surface_weights]]) and
   !> through the potential trace phi the *host* evaluates at the grid points
   !>
   !> Host supplies dE/d(xi, f, r, n) at fixed operator in `coupling%elstat_w_*`,
   !> already contracted with the surface charges
   !>
   !> Unallocated channels are treated as zero, so a host with only a position
   !> response may supply just that one; allocated channels must match `ngrid`
   !>
   !> @param[inout] self     PCM component instance
   !> @param[in]    coupling QM coupling data carrying the host weights
   !> @param[inout] acc      Accumulated cavity surface adjoints
   !> @param[in]    ngrid    Expected electrostatic grid size
   !> @param[out]   error    Error handling
   subroutine pcm_base_get_host_surface_weights(self, coupling, acc, ngrid, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> QM coupling data
      class(coupling_type), intent(in) :: coupling
      !> Accumulated cavity surface adjoints
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Expected electrostatic grid size
      integer, intent(in) :: ngrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      block
         real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :), w_n(:, :)

         allocate (w_xi(ngrid), source=0.0_wp)
         allocate (w_f(ngrid), source=0.0_wp)
         allocate (w_xyz(3, ngrid), source=0.0_wp)
         allocate (w_n(3, ngrid), source=0.0_wp)

         if (allocated(coupling%elstat_w_xi)) then
            if (size(coupling%elstat_w_xi) /= ngrid) then
               call fatal_error(error, "Host electrostatic xi weights have wrong size")
               return
            end if
            w_xi = coupling%elstat_w_xi
         end if
         if (allocated(coupling%elstat_w_f)) then
            if (size(coupling%elstat_w_f) /= ngrid) then
               call fatal_error(error, "Host electrostatic f weights have wrong size")
               return
            end if
            w_f = coupling%elstat_w_f
         end if
         if (allocated(coupling%elstat_w_xyz)) then
            if (size(coupling%elstat_w_xyz, 1) /= 3 .or. &
                size(coupling%elstat_w_xyz, 2) /= ngrid) then
               call fatal_error(error, "Host electrostatic xyz weights have wrong shape")
               return
            end if
            w_xyz = coupling%elstat_w_xyz
         end if
         if (allocated(coupling%elstat_w_n)) then
            if (size(coupling%elstat_w_n, 1) /= 3 .or. &
                size(coupling%elstat_w_n, 2) /= ngrid) then
               call fatal_error(error, "Host electrostatic normal weights have wrong shape")
               return
            end if
            w_n = coupling%elstat_w_n
         end if
         call acc%add_surface_weights(error, w_xi=w_xi, w_f=w_f, &
            & w_xyz=w_xyz, w_n=w_n)
      end block

   end subroutine pcm_base_get_host_surface_weights

   !> Set external matrix (bypasses internal assembly)
   !> Allows user to provide a pre-computed PCM matrix.
   subroutine pcm_set_external_matrix(self, amat)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> External matrix (ngrid, ngrid)
      real(wp), intent(in) :: amat(:, :)

      self%use_external_matrix = .true.
      if (allocated(self%amat)) deallocate (self%amat)
      allocate (self%amat, source=amat)
      self%charges_valid = .false.

   end subroutine pcm_set_external_matrix

   !> Solve the PCM linear system A*q = rhs
   !> Dispatches to appropriate solver based on self%solver setting.
   subroutine pcm_solve_system(self, amat, rhs, q, error)
      use moist_model_component_pcm_solvers, only: solve_pcm_lu, &
         & solve_pcm_cholesky, solve_pcm_iterative, solve_pcm_inversion
      !> PCM component instance
      class(pcm_base), intent(in) :: self
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
      !> self is intent(in) -- only the pointer association would be fixed.
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

   !> Set external electrostatic potential at cavity grid points.
   !> Call this before get_energy to provide the potential computed by a QM
   !> code (e.g. from AO integrals: phi_i = sum_uv P_uv V_uv(r_i)).
   !> The surface charges q solved by get_energy are then available on
   !> self%q for the caller to build its Fock matrix contribution.
   !> @param[in] phi Electrostatic potential at grid points (ngrid)
   subroutine pcm_input_potential(self, phi)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Electrostatic potential at cavity grid points (ngrid)
      real(wp), intent(in) :: phi(:)

      if (allocated(self%phi)) deallocate (self%phi)
      allocate (self%phi, source=phi)
      self%charges_valid = .false.

   end subroutine pcm_input_potential

   !> Electrostatic potential of the solute atomic partial charges on the grid.
   !>
   !> This is the monopole approximation to the solute potential,
   !>
   !>    phi(i) = sum_j q_j / |r_i - R_j|,
   !>
   !> summing the Coulomb potential of the atomic partial charges carried by
   !> `coupling` (over spin channels where present). It does *not* include the
   !> continuous electron density, so it is only as good as the underlying
   !> charge model. Hosts able to supply the exact trace should select
   !> `potential_source%external` instead.
   !>
   !> @param[in]  mol      Molecular structure supplying the atom positions
   !> @param[in]  cavity   Cavity supplying the surface grid points
   !> @param[in]  coupling Wavefunction data carrying the atomic partial charges
   !> @param[out] phi      Potential at each grid point (ngrid)
   subroutine compute_point_charge_potential(mol, cavity, coupling, phi)
      !> Molecular structure supplying the atom positions
      type(structure_type), intent(in) :: mol
      !> Cavity supplying the surface grid points
      class(cavity_type), intent(in) :: cavity
      !> Wavefunction data carrying the atomic partial charges
      class(coupling_type), intent(in) :: coupling
      !> Potential at each grid point (ngrid)
      real(wp), intent(out) :: phi(:)

      !> Grid-point and atom indices with their extents
      integer :: i, j, nat, ngrid
      !> Separation vector, its length, and the atomic partial charge
      real(wp) :: r_vec(3), r_dist, q_atom
      !> Separation below which the singular self-term is skipped
      real(wp), parameter :: min_dist = 1.0e-10_wp

      nat = mol%nat
      ngrid = cavity%ngrid

      phi(:) = 0.0_wp

      ! Compute Coulomb potential from atomic charges
      !  phi _i = Sigma _j q_j / |r_i - R_j | (in atomic units)
      do i = 1, ngrid
         do j = 1, nat
            r_vec(:) = cavity%xyz(:, i) - mol%xyz(:, j)
            r_dist = sqrt(sum(r_vec**2))

            ! Avoid singularities
            if (r_dist < min_dist) cycle

            ! Get atomic charge (sum over spin channels if present)
            if (size(coupling%qat, 2) == 1) then
               q_atom = coupling%qat(j, 1)
            else
               q_atom = sum(coupling%qat(j, :))
            end if

            ! Accumulate Coulomb potential
            phi(i) = phi(i) + q_atom/r_dist
         end do
      end do

   end subroutine compute_point_charge_potential

end module moist_model_component_pcm_type
