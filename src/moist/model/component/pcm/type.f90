!> PCM (Polarizable Continuum Model) abstract base type
!> This module defines the abstract PCM base component that is extended by
!> concrete implementations (CPCM, COSMO, IEF-PCM) in their respective modules.
module moist_model_component_pcm_type
   use mctc_env, only: wp, fatal_error
   use mctc_env_error, only: error_type
   use mctc_io, only: structure_type
   use moist_type, only: solvation_model_component, cavity_type, &
      & potential_type, wavefunction_type
   implicit none
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
   !> Matrix assembly is delegated to the cavity type (cavity%get_amat).
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

   contains

      !> Update PCM component: assembles matrix and prepares for charge solution
      procedure :: update => pcm_base_update

      !> Compute PCM solvation energy
      procedure :: get_energy => pcm_base_get_energy

      !> Compute PCM reaction potential
      procedure :: get_potential => pcm_base_get_potential

      !> Compute PCM gradient with respect to nuclear coordinates
      procedure :: get_gradient => pcm_base_get_gradient

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

      ! Store references
      self%mol_solu = mol
      self%cavity = cavity

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

         ! Delegate matrix assembly to the cavity type
         call cavity%get_amat(self%amat, error)
         if (allocated(error)) return
      else
         if (.not. allocated(self%amat)) then
            call fatal_error(error, &
               & "[pcm_base_update] External PCM matrix requested but not allocated")
            return
         end if
         if (size(self%amat, 1) /= ngrid .or. size(self%amat, 2) /= ngrid) then
            call fatal_error(error, &
               & "[pcm_base_update] External PCM matrix dimension mismatch")
            return
         end if
      end if

      ! Note: charge solving happens in get_energy/get_potential when
      ! wavefunction data (electrostatic potential phi ) is available.

   end subroutine pcm_base_update

   !> Compute PCM solvation energy
   !> E_solv = 0.5 * dot(q, phi)
   !> If an external potential was provided via input_potential, uses that;
   !> otherwise computes the potential internally from atomic point charges.
   subroutine pcm_base_get_energy(self, wfn, energy, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data (used only when no external potential is set)
      type(wavefunction_type), intent(in) :: wfn
      !> Solvation energy (inout to allow accumulation)
      real(wp), intent(inout) :: energy(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: ngrid
      real(wp) :: e_pcm
      real(wp), allocatable :: rhs(:)

      ! Check that update() has prepared cavity and matrix data
      if (.not. allocated(self%cavity) .or. .not. allocated(self%amat)) then
         call fatal_error(error, &
            & "[pcm_base_get_energy] PCM matrix not allocated - call update() first")
         return
      end if

      ngrid = self%cavity%ngrid

      ! Obtain electrostatic potential at cavity grid points
      select case (self%phi_source)
      case (potential_source%charges)
         ! Compute internally from atomic point charges
         if (allocated(self%phi)) deallocate (self%phi)
         allocate (self%phi(ngrid))
         call compute_molecular_potential(self%mol_solu, self%cavity, wfn, self%phi)

      case (potential_source%external)
         ! Use externally provided potential (set via input_potential)
         if (.not. allocated(self%phi)) then
            call fatal_error(error, &
               & "[pcm_base_get_energy] External potential source selected "// &
               & "but phi not set - call input_potential() first")
            return
         end if

      case default
         call fatal_error(error, "[pcm_base_get_energy] Unknown potential source")
         return
      end select

      ! Build RHS: b = -f(eps) * phi
      allocate (rhs(ngrid))
      rhs(:) = -self%feps*self%phi(:)

      ! Solve for charges: A*q = b
      call self%solve_system(self%amat, rhs, self%q, error)
      if (allocated(error)) return

      ! Compute energy: E = 0.5 * dot(q, phi)
      e_pcm = 0.5_wp*dot_product(self%q, self%phi)

      ! Accumulate into energy array
      energy(1) = energy(1) + e_pcm

   end subroutine pcm_base_get_energy

   !> Compute PCM reaction potential
   !> Provides the potential from induced surface charges back to the wavefunction.
   subroutine pcm_base_get_potential(self, wfn, potential, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data
      type(wavefunction_type), intent(in) :: wfn
      !> Solvation potential
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ! Placeholder: implement potential computation
      ! This would compute the reaction field potential at atomic centers
      ! or orbital basis functions from the surface charges q.

      call fatal_error(error, "[pcm_base_get_potential] Not implemented yet")
      return

   end subroutine pcm_base_get_potential

   !> Compute PCM gradient with respect to nuclear coordinates
   !>  dE = 1/2 Sigma _i q_i dphi _i + cavity geometric derivatives
   subroutine pcm_base_get_gradient(self, wfn, gradient, error)
      !> PCM component instance
      class(pcm_base), intent(inout) :: self
      !> Wavefunction data
      type(wavefunction_type), intent(in) :: wfn
      !> Solvation gradient (3, nat)
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ! Placeholder: implement gradient computation
      ! Requires:
      ! 1. dphi _i : gradient of molecular potential at each surface point
      ! 2. Geometric derivatives of cavity (area, normal vectors)
      ! 3. Chain rule for surface charge derivatives

      call fatal_error(error, "[pcm_base_get_gradient] Not implemented yet")
      return

   end subroutine pcm_base_get_gradient

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

      select case (self%solver)
      case (solver_type%lu)
         call solve_pcm_lu(amat, rhs, q, error)

      case (solver_type%cholesky)
         call solve_pcm_cholesky(amat, rhs, q, error)

      case (solver_type%iterative)
         call solve_pcm_iterative(amat, rhs, q, self%solver_tol, &
            & self%solver_maxiter, error)

      case (solver_type%inversion)
         call solve_pcm_inversion(amat, rhs, q, error)

      case default
         call fatal_error(error, "[pcm_solve_system] Unknown solver type")
         return
      end select

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

   end subroutine pcm_input_potential

   !> Compute molecular electrostatic potential at cavity grid points
   !>  phi _i = Sigma _j q_j / |r_i - R_j |
   !> Uses atomic point charges to compute Coulomb potential at each cavity point.
   subroutine compute_molecular_potential(mol, cavity, wfn, phi)
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Cavity with grid points
      class(cavity_type), intent(in) :: cavity
      !> Wavefunction data (contains atomic charges)
      type(wavefunction_type), intent(in) :: wfn
      !> Output: potential at each grid point (ngrid)
      real(wp), intent(out) :: phi(:)

      integer :: i, j, nat, ngrid
      real(wp) :: r_vec(3), r_dist, q_atom
      real(wp), parameter :: min_dist = 1.0e-10_wp  ! Avoid division by zero

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
            if (size(wfn%qat, 2) == 1) then
               q_atom = wfn%qat(j, 1)
            else
               q_atom = sum(wfn%qat(j, :))
            end if

            ! Accumulate Coulomb potential
            phi(i) = phi(i) + q_atom/r_dist
         end do
      end do

   end subroutine compute_molecular_potential

end module moist_model_component_pcm_type
