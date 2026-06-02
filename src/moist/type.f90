
!> Definition of the abstract base solvation model
module moist_type
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_constants, only: pi
   use mctc_io, only: structure_type
   use moist_radius_type, only: radius_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter

   implicit none
   private

   public :: cavity_type
   public :: potential_type
   public :: solvation_model, solvation_model_component
   public :: solver_base_type
   public :: write_cavity_xyz_debug
   public :: write_cavity_csv_debug
   public :: write_cavity_pqr_debug

   !> Abstract base type containing minimal cavity/surface information
   type, abstract :: cavity_type
      !> Area per sphere (nsph)
      real(wp), allocatable :: asph(:)
      !> Total surface area, bohr**2
      real(wp), allocatable :: total_area
      !> Total cavity volume, bohr**3
      real(wp), allocatable :: total_volume
      !> Sphere radii, bohr (nat)
      real(wp), allocatable :: radii(:)
      !> Radii model used to update cached radii.
      class(radius_type), allocatable :: radius_model
      !> Charge of each atom (nat)
      real(wp), allocatable :: qat(:)
      !> Gaussian width of each atom (nat) from eeqbceps
      real(wp), allocatable :: aat(:)
      !> Number of cavity points
      integer :: ngrid
      !> Cartesian coords of points (3,ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Point area, bohr**2 (ngrid)
      real(wp), allocatable :: a(:)
      !> Owner of each grid point (ngrid)
      integer, allocatable :: owner(:)
      !> Electrostatic potential at gridpoints (ngrid)
      real(wp), allocatable :: phi(:)
      !> Convergence flag for each grid point (ngrid)
      logical, allocatable :: converged(:)
      !> Error handling
      type(error_type), allocatable :: error
   contains
      procedure(update_cavity), deferred :: update
      procedure(get_cavity_gradient), deferred :: get_gradient
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
      !> Print detailed cavity information including atomic contributions
      procedure :: print_fine => print_cavity_fine
      !> Assemble PCM interaction matrix A (ngrid, ngrid)
      !> Default: collocation BEM using Coulomb kernel
      procedure :: get_amat => cavity_get_amat_collocation
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

   ! TODO: These will be changed soon

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
   end type potential_type

   !> Type wavefunction data (idential to tblite)
   type, public :: wavefunction_type
      !> Number of electrons for each atom, shape: [nat, spin]
      real(wp), allocatable :: qat(:, :)
      !> Number of electrons for each shell, shape: [nsh, spin]
      real(wp), allocatable :: qsh(:, :)
      !> Atomic dipole moments for each atom, shape: [3, nat, spin]
      real(wp), allocatable :: dpat(:, :, :)
      !> Atomic quadrupole moments for each atom, shape: [5, nat, spin]
      real(wp), allocatable :: qpat(:, :, :)
   end type wavefunction_type

   !> Abstract base solvation model
   type, abstract :: solvation_model

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
      subroutine get_model_energy(self, wfn, energy, error)
         import solvation_model, structure_type, wp, error_type, wavefunction_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Solvation energy
         real(wp), intent(inout) :: energy
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_energy

      !> Get the solvation potential (only for self-consistent models)
      subroutine get_model_potential(self, wfn, potential, error)
         import solvation_model, structure_type, wp, error_type, potential_type, wavefunction_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Solvation potential for the component
         type(potential_type), intent(inout) :: potential
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_potential

      !> Get the solvation energy gradient
      subroutine get_model_gradient(self, wfn, gradient, error)
         import solvation_model, structure_type, wp, error_type, wavefunction_type
         !> Instance of the solvation model
         class(solvation_model), intent(inout) :: self
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Solvation gradient
         real(wp), intent(inout) :: gradient(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_model_gradient

   end interface

   !> Abstract solvation model component
   type, abstract :: solvation_model_component
      !> Name of the component
      character(len=:), allocatable :: name
      !> Abstract cavity type
      class(cavity_type), allocatable :: cavity
      !> Molecular structure data for the component
      type(structure_type) :: mol_solu
      !> Print level for debugging
      integer :: verbosity = 2
      !> Error handling
      type(error_type), allocatable :: error
   contains

      procedure(update_component), deferred :: update
      procedure(get_component_energy), deferred :: get_energy
      procedure(get_component_potential), deferred :: get_potential
      procedure(get_component_gradient), deferred :: get_gradient

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
      subroutine get_component_energy(self, wfn, energy, error)
         import solvation_model_component, wp, wavefunction_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> solvation energy for the component
         real(wp), intent(inout) :: energy(:)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_energy

      !> Get the solvation potential for the component
      subroutine get_component_potential(self, wfn, potential, error)
         import solvation_model_component, potential_type, wavefunction_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Solvation potential for the component
         type(potential_type), intent(inout) :: potential
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine get_component_potential

      !> Get the solvation energy gradient for the component
      subroutine get_component_gradient(self, wfn, gradient, error)
         import solvation_model_component, wp, wavefunction_type, error_type
         !> Instance of the solvation model component
         class(solvation_model_component), intent(inout) :: self
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
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

   !> Write grid points to an XYZ file as helium atoms (debug visualization)
   subroutine write_cavity_xyz_debug(self, filename, error)
      use iso_fortran_env, only: output_unit
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
         if (allocated(self%converged)) then
            if (i <= size(self%converged) .and. self%converged(i)) then
               write (unit, '(a2,1x,3f16.8)') 'He', &
                  self%xyz(1, i)*autoaa, &
                  self%xyz(2, i)*autoaa, &
                  self%xyz(3, i)*autoaa
            else
               write (unit, '(a2,1x,3f16.8)') 'Xe', &
                  self%xyz(1, i)*autoaa, &
                  self%xyz(2, i)*autoaa, &
                  self%xyz(3, i)*autoaa
            end if
         else
            write (unit, '(a2,1x,3f16.8)') 'He', &
               self%xyz(1, i)*autoaa, &
               self%xyz(2, i)*autoaa, &
               self%xyz(3, i)*autoaa
         end if
      end do
      close (unit)

      write (output_unit, '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

   end subroutine write_cavity_xyz_debug

   !> Write grid points to a CSV file (debug visualization)
   subroutine write_cavity_csv_debug(self, filename, error)
      use iso_fortran_env, only: output_unit
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

      write (output_unit, '(a,1x,a)') '[Info] Wrote cavity grid to', trim(filename)

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
      use iso_fortran_env, only: output_unit
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
            (sqrt(self%a(i)/(2.0_wp*pi))*autoaa*0.66_wp + 0.0001_wp)
      end do
      write (unit, '(a)') 'END'
      close (unit)

      write (output_unit, '(a,1x,a)') '[Info] Wrote cavity PQR to', trim(filename)

   end subroutine write_cavity_pqr_debug

   !> Print basic cavity information (grid points, total area, total volume)
   subroutine print_cavity_info(self, unit)
      use iso_fortran_env, only: output_unit
      class(cavity_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer :: iunit
      type(prettyprinter) :: pp

      iunit = output_unit
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

   !> Print detailed cavity information including atomic surface areas
   subroutine print_cavity_fine(self, unit)
      use iso_fortran_env, only: output_unit
      class(cavity_type), intent(in) :: self
      integer, intent(in), optional :: unit
      integer :: iunit, i
      character(32) :: atom_label
      type(prettyprinter) :: pp

      iunit = output_unit
      if (present(unit)) iunit = unit

      pp = new_prettyprinter(unit=iunit, col_value=42, indent_step=2, fmt_len=16)

      ! Print atomic contributions if available
      if (allocated(self%asph)) then
         call pp%push('Atomic surface areas:')
         do i = 1, size(self%asph)
            write (atom_label, '("Atom ",i0)') i
            call pp%kv(trim(atom_label)//' area', self%asph(i), 'bohr^2')
            call pp%kv(trim(atom_label)//' share', 100.0_wp*self%asph(i)/self%total_area, '%')
         end do
         call pp%pop()
      else
         call pp%blank()
         write (iunit, '(a)') '[Info] Atomic surface areas not available'
      end if
      call pp%blank()

   end subroutine print_cavity_fine

   !> Find disconnected grid points / cavities / islands
   subroutine find_disconnected_cavities_base(self, disconnection_thrs, verbose_inp, error)
      use iso_fortran_env, only: output_unit
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

      ! Set threshold (default 4.0 if not provided)
      if (present(disconnection_thrs)) then
         thrs = disconnection_thrs
      else
         thrs = 4.0_wp
      end if

      verbose = 0
      if (present(verbose_inp)) verbose = verbose_inp

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
         write (output_unit, '(a,i0,a,1x,f7.4,1x,a)') '[Info] Estimated average grid spacing: ', &
            nspacing_count, ' points, ', spacing_est/real(nspacing_count, wp), 'bohr'
      end if

      spacing_est = spacing_est/real(nspacing_count, wp)
      ! Use dimensionless threshold to derive connectivity radius.
      cell_size = thrs*spacing_est
      cell_size2 = cell_size*cell_size

      if (verbose > 1) then
         write (output_unit, '(a,1x,f7.4,1x,a)') '[Info] Using cell size for connectivity search:', &
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
            write (output_unit, '(a)') '[Info] No disconnected cavities found.'
            return
         else
            write (output_unit, '(a,i3)') '[Info] Disconnected cavities found:', comp
            write (output_unit, '(1x,a10,a10,a10)') 'id', 'npoints', '%'
            write (output_unit, '(1x,a10,a10,a10)') '---------', '---------', '---------'
            do i = 1, comp
               write (output_unit, '(1x,i10,i10,f10.2)') i, comp_sizes(i), &
                  real(comp_sizes(i), wp)/real(self%ngrid, wp)*100.0_wp
            end do
         end if
      end if

   end subroutine find_disconnected_cavities_base

   !> Assemble PCM interaction matrix using collocation BEM (default).
   !> Builds the Coulomb interaction matrix between cavity grid points:
   !>   Diagonal:     Aii = 2pi / a_i   (self-potential correction)
   !>   Off-diagonal: Aij = -a_j / (4pi |r_i - r_j|)
   !> @param[out] amat  Interaction matrix (ngrid, ngrid)
   !> @param[out] error Error handling
   subroutine cavity_get_amat_collocation(self, amat, error)
      !> Cavity instance
      class(cavity_type), intent(in) :: self
      !> Output: assembled matrix (ngrid, ngrid)
      real(wp), intent(out) :: amat(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i, j, ngrid
      real(wp) :: r_vec(3), r_dist
      real(wp), parameter :: fourpi = 4.0_wp*pi

      ngrid = self%ngrid

      ! Check dimensions
      if (size(amat, 1) /= ngrid .or. size(amat, 2) /= ngrid) then
         call fatal_error(error, &
            & "[cavity_get_amat_collocation] Matrix dimension mismatch")
         return
      end if

      ! Standard collocation BEM interaction matrix
      ! Off-diagonal: Coulomb interaction weighted by area element
      ! Diagonal: self-potential correction for point charges on surface
      do i = 1, ngrid
         do j = 1, ngrid
            if (i == j) then
               amat(i, j) = (2.0_wp*pi)/self%a(i)
            else
               r_vec(:) = self%xyz(:, i) - self%xyz(:, j)
               r_dist = sqrt(sum(r_vec**2))
               amat(i, j) = -self%a(j)/(fourpi*r_dist)
            end if
         end do
      end do

   end subroutine cavity_get_amat_collocation

end module moist_type
