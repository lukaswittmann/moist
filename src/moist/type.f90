!> Definition of the abstract base solvation model
module moist_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   use moist_context, only: moist_context_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none
   private

   public :: cavity_type
   public :: potential_type
   public :: cavity_surface_adjoint_type
   public :: coupling_type
   public :: solvation_model, solvation_model_component
   public :: solver_base_type
   public :: write_cavity_xyz_debug
   public :: write_cavity_csv_debug
   public :: write_cavity_pqr_debug

   !> Abstract base type containing minimal cavity/surface information
   !>
   !> Cavities within moist are per default discretized using Gaussians
   type, abstract :: cavity_type
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller. Never allocated/freed by the cavity.
      type(moist_context_type), pointer :: ctx => null()

      !> Sphere radii, bohr (nat)
      real(wp), allocatable :: radii(:)
      !> Radii model used to update cached radii.
      class(radius_type), allocatable :: radius_model

      !> Number of atomic spheres
      integer :: nsph = 0

      !> Area per sphere (nsph)
      real(wp), allocatable :: asph(:)

      !> Number of cavity points
      integer :: ngrid = 0

      !> Cartesian coordinates of atomic sphere centers (3, nsph)
      real(wp), allocatable :: sphxyz(:, :)

      !> Owner of each grid point (ngrid)
      integer, allocatable :: owner(:)

      !> Cartesian coords of points (3,ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Point area, bohr**2 (ngrid)
      real(wp), allocatable :: a(:)
      !> Gaussian switching factor of each surface point (ngrid)
      real(wp), allocatable :: f(:)
      !> Gaussian width of each surface point (ngrid)
      real(wp), allocatable :: xi0(:)
      !> Outward unit normal of each surface point (3, ngrid)
      real(wp), allocatable :: normal0(:, :)
      !> Point volume element, bohr**3 (ngrid)
      !>
      !> Divergence-theorem partition of the enclosed volume,
      !> v_i = a_i (r_i . n_i)/3, so that `total_volume` is `sum(v)`.
      real(wp), allocatable :: v(:)

      !> Total surface area, bohr**2
      real(wp), allocatable :: total_area
      !> Total cavity volume, bohr**3
      real(wp), allocatable :: total_volume

      !> Nuclear derivatives of surface Gaussian widths (3, nsph, ngrid)
      real(wp), allocatable :: xi1_rA(:, :, :)
      !> Nuclear derivatives of Gaussian switching factors (3, nsph, ngrid)
      real(wp), allocatable :: f1_rA(:, :, :)
      !> Nuclear derivatives of surface positions (3, 3, nsph, ngrid)
      real(wp), allocatable :: xyz1_rA(:, :, :, :)
      !> Nuclear derivatives of surface point volumes (3, nsph, ngrid).
      !> Contracting over the grid gives the total-volume nuclear gradient.
      real(wp), allocatable :: v1_rA(:, :, :)

      !> Error handling
      type(error_type), allocatable :: error
   contains
      procedure(update_cavity), deferred :: update
      procedure(get_cavity_gradient), deferred :: get_gradient
      !> Map accumulated surface-observable adjoints to host potential channels
      procedure :: get_surface_potential => get_cavity_surface_potential_default
      !> Contract accumulated surface-observable adjoints into the nuclear gradient
      procedure :: get_surface_gradient => get_cavity_surface_gradient_default
      !> Write grid to XYZ file for visualization
      procedure :: write_xyz_debug => write_cavity_xyz_debug
      !> Write grid to CSV file for visualization
      procedure :: write_csv_debug => write_cavity_csv_debug
      !> Write grid to PQR file for visualization
      procedure :: write_pqr_debug => write_cavity_pqr_debug
      !> Find disconnected cavities/islands in grid
      procedure :: find_disconnected_cavities => find_disconnected_cavities_base
      !> Print basic cavity information
      procedure :: print => print_cavity_info
   end type cavity_type

   ! Abstract interfaces for deferred procedures
   abstract interface

      subroutine update_cavity(self, mol, error)
         import :: cavity_type, structure_type, wp, error_type
         class(cavity_type), intent(inout) :: self
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_cavity

      subroutine get_cavity_gradient(self)
         import :: cavity_type
         class(cavity_type), intent(inout) :: self
      end subroutine get_cavity_gradient

   end interface

   !> Type for potential data (idential to tblite)
   type :: potential_type
      !> Atom-resolved charge-dependent potential
      real(wp), allocatable :: vat(:, :)
      !> Shell-resolved charge-dependent potential
      real(wp), allocatable :: vsh(:, :)
      !> Orbital-resolved charge-dependent potential
      real(wp), allocatable :: vao(:, :)
      !> Atom-resolved dipolar potential
      real(wp), allocatable :: vdp(:, :, :)
      !> Atom-resolved quadrupolar potential
      real(wp), allocatable :: vqp(:, :, :)
      !> Adjoint weights for cavity level set values (ngrid)
      real(wp), allocatable :: w_lsf0(:)
      !> Adjoint weights for cavity level set gradients (3, ngrid)
      real(wp), allocatable :: w_lsf1(:, :)
      !> Adjoint weights for cavity level set Hessians (3, 3, ngrid)
      real(wp), allocatable :: w_lsf2(:, :, :)
      !> Adjoint weights for the solute density at cavity points (ngrid);
      !> dE/drho_i; the host contracts this with its density evaluator to build
      !> the Fock contribution F_uv += sum_i w_rho(i) phi_u(r_i) phi_v(r_i)
      real(wp), allocatable :: w_rho(:)
      !> Nonlocal electrostatic adjoint dE/dumol on the electrostatic grid
      real(wp), allocatable :: w_elstat_umol(:)
      !> Nonlocal electrostatic adjoint dE/dqmol on the electrostatic grid
      real(wp), allocatable :: w_elstat_qmol(:)
      !> Inner ("elstat") cavity level set adjoint weights for electrostatics.
      !> Same meaning as w_lsf0/1/2 but on the inner cavity grid (ngrid_elstat),
      !> which differs from the outer grid.
      !> The two grids cannot be summed inside moist: the host contracts each
      !> set against the level set basis on its own cavity grid, then sums the
      !> resulting dE/de_c contributions. Only allocated when the model owns an
      !> inner cavity (internal level set path with LevelSet.g_iso_elstat > 0).
      real(wp), allocatable :: w_lsf0_elstat(:)
      real(wp), allocatable :: w_lsf1_elstat(:, :)
      real(wp), allocatable :: w_lsf2_elstat(:, :, :)
   end type potential_type

   !> QM/solute data supplied to solvation models for one coupling step.
   type :: coupling_type
      !> Number of electrons for each atom (nat, spin)
      real(wp), allocatable :: qat(:, :)
      !> Number of electrons for each shell (nsh, spin)
      real(wp), allocatable :: qsh(:, :)
      !> Atomic dipole moments for each atom (3, nat, spin)
      real(wp), allocatable :: dpat(:, :, :)
      !> Atomic quadrupole moments for each atom (5, nat, spin)
      real(wp), allocatable :: qpat(:, :, :)
      !> Molecular potential trace on the electrostatic cavity (ngrid)
      real(wp), allocatable :: elstat_umol(:)
      !> Molecular outward normal-derivative trace on the cavity (ngrid)
      real(wp), allocatable :: elstat_qmol(:)
      !> Direct host trace-geometry weights for Gaussian widths (ngrid)
      real(wp), allocatable :: elstat_w_xi(:)
      !> Direct host trace-geometry weights for switch factors (ngrid)
      real(wp), allocatable :: elstat_w_f(:)
      !> Direct host trace-geometry weights for positions (3, ngrid)
      real(wp), allocatable :: elstat_w_xyz(:, :)
      !> Direct host trace-geometry weights for normals (3, ngrid)
      real(wp), allocatable :: elstat_w_n(:, :)
      !> Charge-weighted electronic field on the electrostatic cavity (3, ngrid)
      !  (distinct from `elstat_w_xyz`, which is an adjoint with respect to surface positions)
      real(wp), allocatable :: elstat_qefield(:, :)
      !> Solute electron density at cavity points, native cavity order (ngrid)
      real(wp), allocatable :: rho_solute(:)
      !> Spatial gradient of the solute density at cavity points (3, ngrid)
      real(wp), allocatable :: rho_solute_gradient(:, :)
   contains
      !> Clear all supplied coupling arrays.
      procedure :: clear => clear_coupling
   end type coupling_type

   !> Abstract base solvation model
   type, abstract :: solvation_model
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller. Never allocated/freed by the model.
      type(moist_context_type), pointer :: ctx => null()

   contains

      procedure(update_model), deferred :: update
      procedure(get_model_energy), deferred :: get_energy
      procedure(get_model_potential), deferred :: get_potential
      procedure(get_model_gradient), deferred :: get_gradient

   end type solvation_model

   abstract interface

      !> Update the solvation model with the current molecular structure
      !> Calculate all structure-dependent properties
      subroutine update_model(self, mol, error)
         import solvation_model, structure_type, error_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_model

      !> Evaluate the solvation energy
      subroutine get_model_energy(self, coupling, energy, error)
         import solvation_model, structure_type, wp, error_type, coupling_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Solvation energy
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_energy

      !> Get the solvation potential (only for self-consistent models)
      subroutine get_model_potential(self, coupling, potential, error)
         import solvation_model, structure_type, wp, error_type, potential_type, coupling_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Solvation potential for the component
         type(potential_type), intent(inout) :: potential
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_potential

      !> Get the solvation energy gradient
      subroutine get_model_gradient(self, coupling, gradient, error)
         import solvation_model, structure_type, wp, error_type, coupling_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Solvation gradient
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_gradient

   end interface

   !> Abstract solvation model component
   type, abstract :: solvation_model_component
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller. Never allocated/freed by the component.
      type(moist_context_type), pointer :: ctx => null()
      !> Name of the component
      character(len=:), allocatable :: name
      !> Molecular structure data for the component
      type(structure_type) :: mol_solu
      !> Linear scale factor applied to this contribution.  The component
      !> multiplies its energy, solvation potential, and surface/level set
      !> response by this constant so the contribution stays variational: 1.0
      !> leaves it unchanged, 0.0 disables it.
      real(wp) :: scale = 1.0_wp
      !> Error handling
      type(error_type), allocatable :: error
   contains

      procedure(update_component), deferred :: update
      procedure(get_component_energy), deferred :: get_energy
      procedure(get_component_potential), deferred :: get_potential
      procedure(get_component_gradient), deferred :: get_gradient
      !> Accumulate direct host-trace adjoints needed before the host can build
      !> its charge-dependent response quantities.
      procedure :: get_trace_potential => get_component_trace_potential_default
      !> Accumulate component-specific surface adjoint weights.
      procedure :: get_surface_weights => get_component_surface_weights_default
      !> Accumulate the host's direct trace-geometry surface adjoint weights.
      procedure :: get_host_surface_weights => get_component_host_surface_weights_default
      !> Accumulate the surface adjoint weights the *nuclear gradient* needs.
      procedure :: get_gradient_surface_weights => get_component_gradient_surface_weights_default
      !> Accumulate nuclear-gradient terms that do not flow through the surface.
      procedure :: get_direct_gradient => get_component_direct_gradient_default

   end type solvation_model_component

   abstract interface

      !> Update the solvation model component with the current molecular structure
      subroutine update_component(self, mol, cavity, error)
         import solvation_model_component, structure_type, cavity_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Cavity type data
         class(cavity_type), intent(inout) :: cavity
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_component

      !> Evaluate the solvation energy for the component
      subroutine get_component_energy(self, coupling, cavity, energy, error)
         import solvation_model_component, cavity_type, wp, coupling_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> solvation energy for the component
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_energy

      !> Get the solvation potential for the component
      subroutine get_component_potential(self, coupling, cavity, potential, error)
         import solvation_model_component, cavity_type, potential_type, coupling_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> Solvation potential for the component
         type(potential_type), intent(inout) :: potential
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_potential

      !> Get the solvation energy gradient for the component
      subroutine get_component_gradient(self, coupling, cavity, gradient, error)
         import solvation_model_component, cavity_type, wp, coupling_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> Solvation gradient for the component
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_gradient

   end interface

   !> Abstract base type for nonlinear solvers
   !>
   !> Note: initialize() is not part of the abstract interface because
   !> different solver types (Newton for equations, SLSQP for optimization)
   !> require different initialization parameters. Each concrete solver
   !> provides its own initialize() method.
   type, abstract :: solver_base_type
   contains
      !> Solve the problem starting from initial guess
      procedure(solve_solver), deferred :: solve

      !> Clean up solver resources
      procedure(destroy_solver), deferred :: destroy
   end type solver_base_type

   abstract interface
      !> Solve the system/optimization problem
      subroutine solve_solver(self, x, error)
         import :: solver_base_type, wp, error_type
         class(solver_base_type), intent(inout), target :: self
         real(wp), dimension(:), intent(inout) :: x  !> Initial guess in, solution out
         type(error_type), allocatable, intent(out) :: error
      end subroutine solve_solver

      !> Destroy solver and free resources
      subroutine destroy_solver(self)
         import :: solver_base_type
         class(solver_base_type), intent(inout), target :: self
      end subroutine destroy_solver
   end interface

contains

   !> Output unit for a cavity: the borrowed run context's unit when one is
   !> attached, otherwise the standard output unit. The base procedures below
   !> are reachable on a cavity that never went through a constructor, so the
   !> association has to be checked rather than assumed.
   pure function cavity_unit(self) result(iunit)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Unit to write to
      integer :: iunit

      iunit = output_unit
      if (associated(self%ctx)) iunit = self%ctx%unit

   end function cavity_unit

   !> Default surface-potential hook for cavities without field-dependent geometry
   !>
   !> @param[inout] self      Cavity instance, unchanged
   !> @param[in]    acc       Surface-observable adjoints, unused
   !> @param[inout] potential Potential accumulator, unchanged
   !> @param[out]   error     Error handling
   subroutine get_cavity_surface_potential_default(self, acc, potential, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Potential accumulator
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_cavity_surface_potential_default

   !> Default reverse-mode nuclear-gradient hook
   !>
   !> Cavities that do not implement the surface-adjoint contraction must be
   !> reached through the forward path instead. Returning silently here would
   !> hand back a zero gradient, so this errors
   !>
   !> @param[in]    self     Cavity instance
   !> @param[in]    acc      Surface-observable adjoints, unused
   !> @param[inout] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error handling
   subroutine get_cavity_surface_gradient_default(self, acc, gradient, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call fatal_error(error, "This cavity does not provide a reverse-mode surface gradient")

   end subroutine get_cavity_surface_gradient_default

   !> Clear all arrays supplied for one QM-solvation coupling step.
   subroutine clear_coupling(self)
      !> Coupling data to clear
      class(coupling_type), intent(inout) :: self

      if (allocated(self%qat)) deallocate (self%qat)
      if (allocated(self%qsh)) deallocate (self%qsh)
      if (allocated(self%dpat)) deallocate (self%dpat)
      if (allocated(self%qpat)) deallocate (self%qpat)
      if (allocated(self%elstat_umol)) deallocate (self%elstat_umol)
      if (allocated(self%elstat_qmol)) deallocate (self%elstat_qmol)
      if (allocated(self%elstat_w_xi)) deallocate (self%elstat_w_xi)
      if (allocated(self%elstat_w_f)) deallocate (self%elstat_w_f)
      if (allocated(self%elstat_w_xyz)) deallocate (self%elstat_w_xyz)
      if (allocated(self%elstat_w_n)) deallocate (self%elstat_w_n)
      if (allocated(self%elstat_qefield)) deallocate (self%elstat_qefield)
      if (allocated(self%rho_solute)) deallocate (self%rho_solute)
      if (allocated(self%rho_solute_gradient)) deallocate (self%rho_solute_gradient)
   end subroutine clear_coupling

   !> Default no-op direct trace-potential hook
   !>
   !> @param[inout] self      Solvation component
   !> @param[in]    coupling  Host coupling data
   !> @param[inout] cavity    Live model cavity
   !> @param[inout] potential Direct trace-potential accumulator
   !> @param[out]   error     Error object
   subroutine get_component_trace_potential_default(self, coupling, cavity, potential, error)
      !> Solvation component
      class(solvation_model_component), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Direct trace-potential accumulator
      type(potential_type), intent(inout) :: potential
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_trace_potential_default

   !> Default no-op surface-weight hook for components without cavity response.
   !> @param[inout] self    Solvation component
   !> @param[in]    coupling     Wavefunction data
   !> @param[in]    cavity  Cavity data
   !> @param[inout] acc     Cavity-specific surface-adjoint accumulator
   !> @param[out]   error   Error object
   subroutine get_component_surface_weights_default(self, coupling, cavity, acc, error)
      !> Solvation component
      class(solvation_model_component), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Cavity data
      class(cavity_type), intent(in) :: cavity
      !> Cavity-specific surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_surface_weights_default

   !> Default no-op host trace-geometry weight hook
   !>
   !> Surface traces built from the host's QM integrals (potential, normal
   !> derivative, ...) carry a surface dependence moist cannot differentiate, so
   !> the host supplies dE/d(xi, f, r, n) at fixed operator in `coupling%elstat_w_*`.
   !>
   !> Components with such a trace override this hook to add those channels to
   !> the shared surface-adjoint accumulator; the rest inherit the no-op.
   !>
   !> @param[inout] self     Solvation component
   !> @param[in]    coupling Wavefunction data carrying the host weights
   !> @param[inout] acc      Surface-adjoint accumulator
   !> @param[in]    ngrid    Expected grid size of the component's cavity
   !> @param[out]   error    Error object
   subroutine get_component_host_surface_weights_default(self, coupling, acc, ngrid, error)
      !> Solvation component
      class(solvation_model_component), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Cavity-specific surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Expected grid size of the component's cavity
      integer, intent(in) :: ngrid
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_host_surface_weights_default

   !> Default gradient-side surface weights: the same ones the potential uses
   !>
   !> For most components the surface adjoint of the energy is one object, so
   !> the reverse-mode nuclear gradient can reuse `get_surface_weights`
   !> verbatim. A component whose gradient legitimately consumes a different
   !> set of host channels overrides this (see `pcm_base`).
   !>
   !> @param[inout] self     Solvation component
   !> @param[in]    coupling Wavefunction data
   !> @param[in]    cavity   Cavity data
   !> @param[inout] acc      Cavity-specific surface-adjoint accumulator
   !> @param[out]   error    Error object
   subroutine get_component_gradient_surface_weights_default(self, coupling, cavity, acc, error)
      !> Solvation component
      class(solvation_model_component), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Cavity data
      class(cavity_type), intent(in) :: cavity
      !> Cavity-specific surface-adjoint accumulator
      class(cavity_surface_adjoint_type), intent(inout) :: acc
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call self%get_surface_weights(coupling, cavity, acc, error)

   end subroutine get_component_gradient_surface_weights_default

   !> Default no-op hook for nuclear-gradient terms outside the surface
   !>
   !> Used by the reverse-mode gradient path for contributions that do not
   !> reach the energy through a cavity surface quantity -- for PCM, the
   !> solute nuclei moving under fixed surface charges.
   !>
   !> @param[inout] self     Solvation component
   !> @param[in]    coupling Wavefunction data
   !> @param[in]    cavity   Cavity data
   !> @param[inout] gradient Nuclear-gradient accumulator, unchanged
   !> @param[out]   error    Error object
   subroutine get_component_direct_gradient_default(self, coupling, cavity, gradient, error)
      !> Solvation component
      class(solvation_model_component), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Cavity data
      class(cavity_type), intent(inout) :: cavity
      !> Nuclear-gradient accumulator
      real(wp), intent(inout) :: gradient(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_direct_gradient_default

   !> Write grid points to an XYZ file as helium atoms (debug visualization)
   subroutine write_cavity_xyz_debug(self, filename, error)
      use mctc_io_convert, only: autoaa
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, 'write_xyz_debug: cavity grid not allocated')
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, 'write_xyz_debug: no grid points to write')
         return
      end if

      open (file=filename, newunit=unit, status='replace', action='write', iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, 'Could not open XYZ file for writing: '//trim(filename))
         return
      end if

      write (unit, '(i0)') self%ngrid
      write (unit, '(a)') 'drop cavity grid points as He (Angstrom)'
      do i = 1, self%ngrid
         write (unit, '(a2,1x,3f16.8)') 'He', &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa
      end do
      close (unit)

      write (cavity_unit(self), '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

   end subroutine write_cavity_xyz_debug

   !> Write grid points to a CSV file (debug visualization)
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: stat, i, unit

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, 'write_csv_debug: cavity grid not allocated')
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, 'write_csv_debug: point areas not allocated')
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, 'write_csv_debug: point owners not allocated')
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, 'write_csv_debug: no grid points to write')
         return
      end if

      open (file=filename, newunit=unit, status='replace', action='write', iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, 'Could not open CSV file for writing: '//trim(filename))
         return
      end if

      write (unit, '(a)') 'ngrid,x,y,z,owner,area'

      do i = 1, self%ngrid
         write (unit, '(i0,7('','',g0))') i, &
            self%xyz(1, i), self%xyz(2, i), self%xyz(3, i), &
            self%owner(i), self%a(i)
      end do
      close (unit)

      write (cavity_unit(self), '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

   end subroutine write_cavity_csv_debug

   !> Write grid points to a PQR file (debug visualization)
   !>
   !> Grid points are written as HETATM records with:
   !> - positions converted from bohr to Angstrom
   !> - charge set to 0.0
   !> - radius set to the final adapted integration weight `a(i)` (area element
   !>   with switching function applied), also converted to Angstrom
   !>
   !> @param[in]  self      Cavity instance
   !> @param[in]  filename  Output PQR file path
   subroutine write_cavity_pqr_debug(self, filename, error)
      use mctc_io_convert, only: autoaa
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, i

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, 'write_pqr_debug: cavity grid not allocated')
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, 'write_pqr_debug: point areas not allocated')
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, 'write_pqr_debug: point owners not allocated')
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, 'write_pqr_debug: no grid points to write')
         return
      end if

      open (file=filename, newunit=unit, status='replace', action='write', iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, 'Could not open PQR file for writing: '//trim(filename))
         return
      end if

      do i = 1, self%ngrid
         write (unit, '(a6,i5,1x,a4,a1,a3,1x,a1,i4,4x,3f8.3,f8.4,f7.4)') &
            'HETATM', i, 'GP  ', ' ', 'GRD', 'A', self%owner(i), &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa, &
            0.0_wp, &
            (sqrt(self%a(i)/(2.0_wp*pi))*autoaa + 0.0001_wp)
      end do
      write (unit, '(a)') 'END'
      close (unit)

      write (cavity_unit(self), '(a,1x,a)') '[Info] Wrote cavity PQR to', trim(filename)

   end subroutine write_cavity_pqr_debug
   

   !> Print basic cavity information (grid points, total area, total volume)
   subroutine print_cavity_info(self, unit)
      class(cavity_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer :: iunit
      type(prettyprinter) :: pp

      iunit = cavity_unit(self)
      if (present(unit)) iunit = unit

      if (.not. allocated(self%total_area) .or. .not. allocated(self%total_volume)) then
         write (iunit, '(a)') '[Warning] Cavity not fully initialized'
         return
      end if

      pp = new_prettyprinter(unit=iunit, fmt_len=20)

      call pp%blank()
      call pp%push('Results:')
      call pp%kv('Cavity points', self%ngrid)
      call pp%kv('Total area', self%total_area, 'bohr^2')
      call pp%kv('Total volume', self%total_volume, 'bohr^3')
      call pp%pop()
      call pp%blank()

   end subroutine print_cavity_info


   !> Find disconnected grid points / cavities / islands
   subroutine find_disconnected_cavities_base(self, disconnection_thrs, verbose_inp, error)
      class(cavity_type), intent(inout) :: self
      real(wp), intent(in), optional :: disconnection_thrs
      integer, intent(in), optional :: verbose_inp
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: cell_size, cell_size2, spacing_est, spacing_guess
      real(wp) :: min_xyz(3), max_xyz(3), bbox(3), volume
      real(wp) :: dx, dy, dz, dist2, min_dist2
      integer :: nx, ny, nz, ncell
      integer :: i, ix, iy, iz, lin, neighbour, qhead, qtail, current
      integer :: nxi, nyi, nzi, comp, search_rad
      integer :: nspacing_count
      integer :: alloc_stat
      integer, allocatable :: head(:), next(:), queue(:)
      integer, allocatable :: cell_ix(:), cell_iy(:), cell_iz(:)
      integer, allocatable :: comp_sizes(:)
      logical, allocatable :: visited(:)
      real(wp) :: thrs
      integer :: verbose
      integer :: iunit

      ! Set threshold (default 4.0 if not provided)
      if (present(disconnection_thrs)) then
         thrs = disconnection_thrs
      else
         thrs = 4.0_wp
      end if

      !> Silent unless the caller asks for output. This deliberately does not
      !> follow the context verbosity: the only caller passes no verbose_inp,
      !> and defaulting from the context would start printing island tables on
      !> every verbosity-2 run.
      verbose = 0
      if (present(verbose_inp)) verbose = verbose_inp

      iunit = cavity_unit(self)

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, 'find_disconnected_cavities: grid not allocated')
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, 'find_disconnected_cavities: no grid points')
         return
      end if

      min_xyz = [minval(self%xyz(1, :)), minval(self%xyz(2, :)), minval(self%xyz(3, :))]
      max_xyz = [maxval(self%xyz(1, :)), maxval(self%xyz(2, :)), maxval(self%xyz(3, :))]
      bbox = max_xyz - min_xyz
      volume = max(1.0e-12_wp, bbox(1)*bbox(2)*bbox(3))
      spacing_guess = max(1.0e-6_wp, (volume/real(self%ngrid, wp))**(1.0_wp/3.0_wp))

      ! First pass: estimate average nearest-neighbour spacing with a coarse grid.
      cell_size = spacing_guess
      nx = max(1, int(bbox(1)/cell_size) + 1)
      ny = max(1, int(bbox(2)/cell_size) + 1)
      nz = max(1, int(bbox(3)/cell_size) + 1)
      ncell = max(1, nx*ny*nz)

      allocate (head(ncell), source=0, stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for head')
         return
      end if
      allocate (next(self%ngrid), source=0, stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for next')
         return
      end if
      allocate (cell_ix(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for cell_ix')
         return
      end if
      allocate (cell_iy(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for cell_iy')
         return
      end if
      allocate (cell_iz(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for cell_iz')
         return
      end if

      do i = 1, self%ngrid
         ix = min(nx, max(1, int((self%xyz(1, i) - min_xyz(1))/cell_size) + 1))
         iy = min(ny, max(1, int((self%xyz(2, i) - min_xyz(2))/cell_size) + 1))
         iz = min(nz, max(1, int((self%xyz(3, i) - min_xyz(3))/cell_size) + 1))

         cell_ix(i) = ix
         cell_iy(i) = iy
         cell_iz(i) = iz
         lin = ix + nx*(iy - 1 + ny*(iz - 1))
         next(i) = head(lin)
         head(lin) = i
      end do

      ! Accumulate nearest-neighbour distances to form a characteristic spacing.
      spacing_est = 0.0_wp
      nspacing_count = 0
      do i = 1, self%ngrid
         ix = cell_ix(i); iy = cell_iy(i); iz = cell_iz(i)
         min_dist2 = huge(1.0_wp)
         do search_rad = 0, 2
            do nzi = max(1, iz - search_rad), min(nz, iz + search_rad)
               do nyi = max(1, iy - search_rad), min(ny, iy + search_rad)
                  do nxi = max(1, ix - search_rad), min(nx, ix + search_rad)
                     lin = nxi + nx*(nyi - 1 + ny*(nzi - 1))
                     neighbour = head(lin)
                     do while (neighbour /= 0)
                        if (neighbour /= i) then
                           dx = self%xyz(1, neighbour) - self%xyz(1, i)
                           dy = self%xyz(2, neighbour) - self%xyz(2, i)
                           dz = self%xyz(3, neighbour) - self%xyz(3, i)
                           dist2 = dx*dx + dy*dy + dz*dz
                           if (dist2 < min_dist2) min_dist2 = dist2
                        end if
                        neighbour = next(neighbour)
                     end do
                  end do
               end do
            end do
            if (min_dist2 < huge(1.0_wp)) exit
         end do

         if (min_dist2 < huge(1.0_wp)) then
            spacing_est = spacing_est + sqrt(min_dist2)
            nspacing_count = nspacing_count + 1
         end if
      end do

      if (nspacing_count == 0 .or. spacing_est <= 0.0_wp) then
         call fatal_error(error, 'find_disconnected_cavities: could not estimate grid spacing')
         deallocate (head, next, cell_ix, cell_iy, cell_iz)
         return
      end if

      if (verbose > 1) then
         write (iunit, '(a,i0,a,1x,f7.4,1x,a)') '[Info] Estimated average grid spacing: ', &
            nspacing_count, ' points, ', spacing_est/real(nspacing_count, wp), 'bohr'
      end if

      spacing_est = spacing_est/real(nspacing_count, wp)
      ! Use dimensionless threshold to derive connectivity radius.
      cell_size = thrs*spacing_est
      cell_size2 = cell_size*cell_size

      if (verbose > 1) then
         write (iunit, '(a,1x,f7.4,1x,a)') '[Info] Using cell size for connectivity search:', &
            cell_size, 'bohr'
      end if

      ! Rebuild grid for connectivity search with final cell size.
      nx = max(1, int((bbox(1)/cell_size)) + 1)
      ny = max(1, int((bbox(2)/cell_size)) + 1)
      nz = max(1, int((bbox(3)/cell_size)) + 1)
      ncell = max(1, nx*ny*nz)

      if (size(head) /= ncell) then
         deallocate (head)
         allocate (head(ncell), source=0, stat=alloc_stat)
         if (alloc_stat /= 0) then
            call fatal_error(error, 'find_disconnected_cavities: allocation failed for head resize')
            return
         end if
      else
         head = 0
      end if
      next = 0
      do i = 1, self%ngrid
         ix = min(nx, max(1, int((self%xyz(1, i) - min_xyz(1))/cell_size) + 1))
         iy = min(ny, max(1, int((self%xyz(2, i) - min_xyz(2))/cell_size) + 1))
         iz = min(nz, max(1, int((self%xyz(3, i) - min_xyz(3))/cell_size) + 1))

         cell_ix(i) = ix
         cell_iy(i) = iy
         cell_iz(i) = iz
         lin = ix + nx*(iy - 1 + ny*(iz - 1))
         next(i) = head(lin)
         head(lin) = i
      end do

      if (allocated(queue)) deallocate (queue)
      if (allocated(comp_sizes)) deallocate (comp_sizes)
      if (allocated(visited)) deallocate (visited)
      allocate (queue(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for queue')
         return
      end if
      allocate (comp_sizes(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for comp_sizes')
         return
      end if
      allocate (visited(self%ngrid), source=.false., stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, 'find_disconnected_cavities: allocation failed for visited')
         return
      end if

      ! BFS to label connected components using distance-limited neighbours.
      comp = 0
      do i = 1, self%ngrid
         if (visited(i)) cycle

         comp = comp + 1
         comp_sizes(comp) = 0

         qhead = 1
         qtail = 1
         queue(1) = i
         visited(i) = .true.

         do while (qhead <= qtail)
            current = queue(qhead)
            qhead = qhead + 1
            comp_sizes(comp) = comp_sizes(comp) + 1

            ix = cell_ix(current)
            iy = cell_iy(current)
            iz = cell_iz(current)

            do nzi = max(1, iz - 1), min(nz, iz + 1)
               do nyi = max(1, iy - 1), min(ny, iy + 1)
                  do nxi = max(1, ix - 1), min(nx, ix + 1)
                     lin = nxi + nx*(nyi - 1 + ny*(nzi - 1))
                     neighbour = head(lin)

                     do while (neighbour /= 0)
                        if (.not. visited(neighbour)) then
                           dx = self%xyz(1, neighbour) - self%xyz(1, current)
                           dy = self%xyz(2, neighbour) - self%xyz(2, current)
                           dz = self%xyz(3, neighbour) - self%xyz(3, current)
                           dist2 = dx*dx + dy*dy + dz*dz

                           if (dist2 <= cell_size2) then
                              visited(neighbour) = .true.
                              qtail = qtail + 1
                              queue(qtail) = neighbour
                           end if
                        end if
                        neighbour = next(neighbour)
                     end do
                  end do
               end do
            end do
         end do
      end do

      if (verbose > 1) then
         if (comp == 1) then
            write (iunit, '(a)') '[Info] No disconnected cavities found.'
            return
         else
            write (iunit, '(a,i3)') '[Info] Disconnected cavities found:', comp
            write (iunit, '(1x,a10,a10,a10)') 'id', 'npoints', '%'
            write (iunit, '(1x,a10,a10,a10)') '---------', '---------', '---------'
            do i = 1, comp
               write (iunit, '(1x,i10,i10,f10.2)') i, comp_sizes(i), &
                  real(comp_sizes(i), wp)/real(self%ngrid, wp)*100.0_wp
            end do
         end if
      end if

   end subroutine find_disconnected_cavities_base

end module moist_type
