!> Definition of the abstract base solvation model
module moist_type
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   use moist_context, only: moist_context_type
   use moist_cavity_surface_adjoint, only: cavity_surface_adjoint_type
   use moist_channels, only: coupling_type, response_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none
   private

   public :: cavity_type
   public :: cavity_surface_adjoint_type
   public :: solvation_model_type, solvation_model_component_type
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
   contains
      procedure(update_cavity), deferred :: update
      procedure(get_cavity_gradient), deferred :: get_gradient
      !> Map accumulated surface-observable adjoints to host response channels
      procedure :: get_surface_response => get_cavity_surface_response_default
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

      subroutine get_cavity_gradient(self, error)
         import :: cavity_type, error_type
         class(cavity_type), intent(inout) :: self
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_cavity_gradient

   end interface

   !> Abstract base solvation model
   type, abstract :: solvation_model_type
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller. Never allocated/freed by the model.
      type(moist_context_type), pointer :: ctx => null()

   contains

      procedure(update_model), deferred :: update
      procedure(get_model_energy), deferred :: get_energy
      procedure(get_model_response), deferred :: get_response
      procedure(get_model_gradient), deferred :: get_gradient

   end type solvation_model_type

   abstract interface

      !> Update the solvation model with the current molecular structure
      !> Calculate all structure-dependent properties
      subroutine update_model(self, mol, error)
         import solvation_model_type, structure_type, error_type
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_model

      !> Evaluate the solvation energy
      subroutine get_model_energy(self, coupling, energy, error)
         import solvation_model_type, structure_type, wp, error_type, coupling_type
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Solvation energy
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_energy

      !> Get the solvation response (only for self-consistent models)
      subroutine get_model_response(self, coupling, response, error)
         import solvation_model_type, structure_type, wp, error_type, response_type, coupling_type
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Solvation response for the component
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_response

      !> Get the solvation energy gradient
      subroutine get_model_gradient(self, coupling, gradient, error)
         import solvation_model_type, structure_type, wp, error_type, coupling_type
         !> Instance of the solvation model
         class(solvation_model_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Solvation gradient
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_gradient

   end interface

   !> Abstract solvation model component
   type, abstract :: solvation_model_component_type
      !> Borrowed run context (verbosity/debug/timer); set at construction,
      !> owned by the top-level caller. Never allocated/freed by the component.
      type(moist_context_type), pointer :: ctx => null()
      !> Name of the component
      character(len=:), allocatable :: name
      !> Molecular structure data for the component
      type(structure_type) :: mol_solu
      !> Linear scale factor applied to this contribution.  The component
      !> multiplies its energy, solvation response, and surface/level set
      !> response by this constant so the contribution stays variational: 1.0
      !> leaves it unchanged, 0.0 disables it.
      real(wp) :: scale = 1.0_wp
      !> Error handling
      type(error_type), allocatable :: error
   contains

      procedure(update_component), deferred :: update
      procedure(get_component_energy), deferred :: get_energy
      procedure(get_component_response), deferred :: get_response
      procedure(get_component_gradient), deferred :: get_gradient
      !> Accumulate direct host-trace adjoints needed before the host can build
      !> its charge-dependent response quantities.
      procedure :: get_trace_response => get_component_trace_response_default
      !> Accumulate component-specific surface adjoint weights.
      procedure :: get_surface_weights => get_component_surface_weights_default
      !> Accumulate the host's direct trace-geometry surface adjoint weights.
      procedure :: get_host_surface_weights => get_component_host_surface_weights_default
      !> Accumulate the surface adjoint weights the *nuclear gradient* needs.
      procedure :: get_gradient_surface_weights => get_component_gradient_surface_weights_default
      !> Accumulate nuclear-gradient terms that do not flow through the surface.
      procedure :: get_direct_gradient => get_component_direct_gradient_default

   end type solvation_model_component_type

   abstract interface

      !> Update the solvation model component with the current molecular structure
      subroutine update_component(self, mol, cavity, error)
         import solvation_model_component_type, structure_type, cavity_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component_type), intent(inout) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Cavity type data
         class(cavity_type), intent(inout) :: cavity
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine update_component

      !> Evaluate the solvation energy for the component
      subroutine get_component_energy(self, coupling, cavity, energy, error)
         import solvation_model_component_type, cavity_type, wp, coupling_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> solvation energy for the component
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_energy

      !> Get the solvation response for the component
      subroutine get_component_response(self, coupling, cavity, response, error)
         import solvation_model_component_type, cavity_type, response_type, coupling_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component_type), intent(inout) :: self
         !> Wavefunction data
         class(coupling_type), intent(in) :: coupling
         !> Live cavity owned by the orchestrating model
         class(cavity_type), intent(inout) :: cavity
         !> Solvation response for the component
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_response

      !> Get the solvation energy gradient for the component
      subroutine get_component_gradient(self, coupling, cavity, gradient, error)
         import solvation_model_component_type, cavity_type, wp, coupling_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component_type), intent(inout) :: self
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

   !> Default surface-response hook for cavities without field-dependent geometry
   !>
   !> @param[inout] self      Cavity instance, unchanged
   !> @param[in]    acc       Surface-observable adjoints, unused
   !> @param[inout] response  Response accumulator, unchanged
   !> @param[out]   error     Error handling
   subroutine get_cavity_surface_response_default(self, acc, response, error)
      !> Cavity instance
      class(cavity_type), intent(inout) :: self
      !> Surface-observable adjoints
      type(cavity_surface_adjoint_type), intent(in) :: acc
      !> Response accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_cavity_surface_response_default

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

   !> Default no-op direct trace-response hook
   !>
   !> @param[inout] self      Solvation component
   !> @param[in]    coupling  Host coupling data
   !> @param[inout] cavity    Live model cavity
   !> @param[inout] response  Direct trace-response accumulator
   !> @param[out]   error     Error object
   subroutine get_component_trace_response_default(self, coupling, cavity, response, error)
      !> Solvation component
      class(solvation_model_component_type), intent(inout) :: self
      !> Host coupling data
      class(coupling_type), intent(in) :: coupling
      !> Live model cavity
      class(cavity_type), intent(inout) :: cavity
      !> Direct trace-response accumulator
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

   end subroutine get_component_trace_response_default

   !> Default no-op surface-weight hook for components without cavity response.
   !> @param[inout] self    Solvation component
   !> @param[in]    coupling     Wavefunction data
   !> @param[in]    cavity  Cavity data
   !> @param[inout] acc     Cavity-specific surface-adjoint accumulator
   !> @param[out]   error   Error object
   subroutine get_component_surface_weights_default(self, coupling, cavity, acc, error)
      !> Solvation component
      class(solvation_model_component_type), intent(inout) :: self
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
   !> the host supplies dE/d(xi, f, r, n) at fixed operator in `coupling%electrostatics%w_*`.
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
      class(solvation_model_component_type), intent(inout) :: self
      !> Wavefunction data
      class(coupling_type), intent(in) :: coupling
      !> Cavity-specific surface-adjoint accumulator
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
   !> verbatim. A component whose gradient legitimately consumes a different
   !> set of host channels overrides this (see `solvation_model_component_pcm`).
   !>
   !> @param[inout] self     Solvation component
   !> @param[in]    coupling Wavefunction data
   !> @param[in]    cavity   Cavity data
   !> @param[inout] acc      Cavity-specific surface-adjoint accumulator
   !> @param[out]   error    Error object
   subroutine get_component_gradient_surface_weights_default(self, coupling, cavity, acc, error)
      !> Solvation component
      class(solvation_model_component_type), intent(inout) :: self
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
      class(solvation_model_component_type), intent(inout) :: self
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
         call fatal_error(error, "write_xyz_debug: cavity grid not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_xyz_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open XYZ file for writing: "//trim(filename))
         return
      end if

      write (unit, "(i0)") self%ngrid
      write (unit, "(a)") "drop cavity grid points as He (Angstrom)"
      do i = 1, self%ngrid
         write (unit, "(a2,1x,3f16.8)") "He", &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa
      end do
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

   end subroutine write_cavity_xyz_debug

   !> Write grid points to a CSV file (debug visualization)
   subroutine write_cavity_csv_debug(self, filename, error)
      class(cavity_type), intent(in) :: self
      character(len=*), intent(in) :: filename
      type(error_type), allocatable, intent(out) :: error

      integer :: stat, i, unit

      if (.not. allocated(self%xyz)) then
         call fatal_error(error, "write_csv_debug: cavity grid not allocated")
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, "write_csv_debug: point areas not allocated")
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, "write_csv_debug: point owners not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_csv_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open CSV file for writing: "//trim(filename))
         return
      end if

      write (unit, "(a)") "ngrid,x,y,z,owner,area"

      do i = 1, self%ngrid
         write (unit, '(i0,7('','',g0))') i, &
            self%xyz(1, i), self%xyz(2, i), self%xyz(3, i), &
            self%owner(i), self%a(i)
      end do
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity grid to", trim(filename)

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
         call fatal_error(error, "write_pqr_debug: cavity grid not allocated")
         return
      end if
      if (.not. allocated(self%a)) then
         call fatal_error(error, "write_pqr_debug: point areas not allocated")
         return
      end if
      if (.not. allocated(self%owner)) then
         call fatal_error(error, "write_pqr_debug: point owners not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "write_pqr_debug: no grid points to write")
         return
      end if

      open (file=filename, newunit=unit, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Could not open PQR file for writing: "//trim(filename))
         return
      end if

      do i = 1, self%ngrid
         write (unit, "(a6,i5,1x,a4,a1,a3,1x,a1,i4,4x,3f8.3,f8.4,f7.4)") &
            "HETATM", i, "GP  ", " ", "GRD", "A", self%owner(i), &
            self%xyz(1, i)*autoaa, &
            self%xyz(2, i)*autoaa, &
            self%xyz(3, i)*autoaa, &
            0.0_wp, &
            (sqrt(self%a(i)/(2.0_wp*pi))*autoaa + 0.0001_wp)
      end do
      write (unit, "(a)") "END"
      close (unit)

      write (cavity_unit(self), "(a,1x,a)") "[Info] Wrote cavity PQR to", trim(filename)

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
         write (iunit, "(a)") "[Warning] Cavity not fully initialized"
         return
      end if

      pp = new_prettyprinter(unit=iunit, fmt_len=20)

      call pp%blank()
      call pp%push("Results:")
      call pp%kv("Cavity points", self%ngrid)
      call pp%kv("Total area", self%total_area, "bohr^2")
      call pp%kv("Total volume", self%total_volume, "bohr^3")
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
         call fatal_error(error, "find_disconnected_cavities: grid not allocated")
         return
      end if
      if (self%ngrid <= 0) then
         call fatal_error(error, "find_disconnected_cavities: no grid points")
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
         call fatal_error(error, "find_disconnected_cavities: allocation failed for head")
         return
      end if
      allocate (next(self%ngrid), source=0, stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, "find_disconnected_cavities: allocation failed for next")
         return
      end if
      allocate (cell_ix(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, "find_disconnected_cavities: allocation failed for cell_ix")
         return
      end if
      allocate (cell_iy(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, "find_disconnected_cavities: allocation failed for cell_iy")
         return
      end if
      allocate (cell_iz(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, "find_disconnected_cavities: allocation failed for cell_iz")
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
         call fatal_error(error, "find_disconnected_cavities: could not estimate grid spacing")
         deallocate (head, next, cell_ix, cell_iy, cell_iz)
         return
      end if

      if (verbose > 1) then
         write (iunit, "(a,i0,a,1x,f7.4,1x,a)") "[Info] Estimated average grid spacing: ", &
            nspacing_count, " points, ", spacing_est/real(nspacing_count, wp), "bohr"
      end if

      spacing_est = spacing_est/real(nspacing_count, wp)
      ! Use dimensionless threshold to derive connectivity radius.
      cell_size = thrs*spacing_est
      cell_size2 = cell_size*cell_size

      if (verbose > 1) then
         write (iunit, "(a,1x,f7.4,1x,a)") "[Info] Using cell size for connectivity search:", &
            cell_size, "bohr"
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
            call fatal_error(error, "find_disconnected_cavities: allocation failed for head resize")
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
         call fatal_error(error, "find_disconnected_cavities: allocation failed for queue")
         return
      end if
      allocate (comp_sizes(self%ngrid), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, "find_disconnected_cavities: allocation failed for comp_sizes")
         return
      end if
      allocate (visited(self%ngrid), source=.false., stat=alloc_stat)
      if (alloc_stat /= 0) then
         call fatal_error(error, "find_disconnected_cavities: allocation failed for visited")
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
            write (iunit, "(a)") "[Info] No disconnected cavities found."
            return
         else
            write (iunit, "(a,i3)") "[Info] Disconnected cavities found:", comp
            write (iunit, "(1x,a10,a10,a10)") "id", "npoints", "%"
            write (iunit, "(1x,a10,a10,a10)") "---------", "---------", "---------"
            do i = 1, comp
               write (iunit, "(1x,i10,i10,f10.2)") i, comp_sizes(i), &
                  real(comp_sizes(i), wp)/real(self%ngrid, wp)*100.0_wp
            end do
         end if
      end if

   end subroutine find_disconnected_cavities_base

end module moist_type
