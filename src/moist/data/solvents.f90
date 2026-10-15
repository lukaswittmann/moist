
module moist_data_solvents
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new_structure
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use mctc_io_convert, only: autokcal, aatoau
   use mctc_io_codata2018, only: Avogadro_constant, Bohr_radius
   use mctc_io_codata2018, only: Hartree_energy, atomic_unit_of_mass
   use mctc_env_error, only: error_type, fatal_error
   use moist_data_mass, only: get_mass
   use mctc_io_utils, only: to_lower
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)

   integer, parameter, public :: max_solvents = 180

   !> Longest canonical solvent name, MOIST_NAME_MAX in the C API
   integer, parameter, public :: max_solvent_name_len = 32

   public :: solvent_system_type, new_solvent_system
   public :: get_solvent_id
   public :: solvent_descriptors_type, get_solvent_descriptors, get_solvent_descriptors_by_id
   public :: solvent_multipole_data_type, get_solvent_charges, get_solvent_multipoles

   private

   !> Bulk descriptors of one solvent table entry
   !>
   !> Table values at the temperatures of solvents.inc; a quiet NaN marks a
   !> value the table lacks and carries into the values derived from it
   type :: solvent_descriptors_type
      !> Table ID, frozen: new solvents are appended
      integer :: id = 0
      !> Canonical name, lowercase
      character(:), allocatable :: name
      !> Static relative permittivity
      real(wp) :: epsilon = 0.0_wp
      !> Refractive index at the sodium D line
      real(wp) :: refractive_index = 0.0_wp
      !> Abraham hydrogen-bond acidity
      real(wp) :: abraham_alpha = 0.0_wp
      !> Abraham hydrogen-bond basicity
      real(wp) :: abraham_beta = 0.0_wp
      !> Macroscopic surface tension in N/m
      real(wp) :: surface_tension_si = 0.0_wp
      !> Macroscopic surface tension in Eh/bohr**2
      real(wp) :: surface_tension_au = 0.0_wp
      !> Mass density in kg/m**3
      real(wp) :: mass_density_si = 0.0_wp
      !> Mass density in me/bohr**3
      real(wp) :: mass_density_au = 0.0_wp
      !> Molar mass in kg/mol
      real(wp) :: molar_mass_si = 0.0_wp
      !> Mass of one molecule in me
      real(wp) :: molecular_mass_au = 0.0_wp
      !> Number density in 1/m**3
      real(wp) :: number_density_si = 0.0_wp
      !> Number density in 1/bohr**3
      real(wp) :: number_density_au = 0.0_wp
      !> Volume per molecule in m**3
      real(wp) :: molecular_volume_si = 0.0_wp
      !> Volume per molecule in bohr**3
      real(wp) :: molecular_volume_au = 0.0_wp
   contains
      procedure :: print => print_solvent_descriptors
   end type solvent_descriptors_type

   !> One solvent at a thermodynamic state, with its geometry
   type :: solvent_system_type
      !> Table values and derived bulk properties
      type(solvent_descriptors_type) :: solvent
      !> Temperature in K
      real(wp) :: temperature
      !> Pressure in Pa
      real(wp) :: pressure_si
      !> Pressure in Eh/bohr**3
      real(wp) :: pressure_au
      !> Solvent geometry in bohr
      type(structure_type), allocatable :: solv_mol
   contains
      procedure :: print => print_solvent_system
      procedure :: get_charges => get_system_charges
      procedure :: get_multipoles => get_system_multipoles
   end type solvent_system_type

   !> Atom-centered Cartesian multipoles in e*bohr^rank
   type :: solvent_multipole_data_type
      real(wp), allocatable :: monopole(:)
      real(wp), allocatable :: dipole(:,:)
      real(wp), allocatable :: quadrupole(:,:)
      real(wp), allocatable :: octupole(:,:)
   end type solvent_multipole_data_type

contains

   !> Get the solvent ID from an alias
   subroutine get_solvent_id(alias, solvent_id, error)

      !> Solvent alias (case insensitive)
      character(len=*), intent(in) :: alias

      !> Solvent ID
      integer, intent(out) :: solvent_id

      !> Error output for unknown aliases
      type(error_type), allocatable, intent(out) :: error

      !> Iterables
      integer :: i, j

      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)
      character(:), allocatable :: query

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho

      !> Get basic solvent information
      include "solvents.inc"

      solvent_id = 0

      ! Reject blank queries, which would otherwise match the blank alias padding
      if (len_trim(alias) == 0) then
         call fatal_error(error, message="Empty solvent alias", stat=1)
         return
      end if

      ! Normalise
      query = trim(adjustl(to_lower(alias)))

      ! Search for the alias in the alias list, normalised like the query
      do i = 1, max_solvents
         do j = 1, 10
            if (len_trim(alias_list(j, i)) == 0) cycle
            if (query == trim(adjustl(to_lower(alias_list(j, i))))) then
               solvent_id = id_list(i)
               return
            end if
         end do
      end do

      ! If we reach here, the alias was not found
      call fatal_error(error, message="Unknown solvent: "//trim(alias), stat=1)

   end subroutine get_solvent_id

   !> Include subroutine to get the solvent geometry
   !>
   !> The routine is of the form
   !>
   !>   subroutine get_solvent_geometry(solvent_id, mol)
   !>     integer, intent(in) :: solvent_id
   !>     type(structure_type), intent(out) :: mol
   !>   end subroutine get_solvent_geometry
   !>
   include "solventgeometries.inc"

   !> Atomic charges in e
   include "solventcharges.inc"

   !> Signed, atom-centered Cartesian moments in e*bohr^rank
   include "solventmultipoles.inc"

   !> Get atomic charges for a gas, solvent, or conductor environment
   subroutine get_system_charges(self, environment, model, charges, error)
      class(solvent_system_type), intent(in) :: self
      character(len=*), intent(in) :: environment
      character(len=*), intent(in) :: model
      real(wp), allocatable, intent(out) :: charges(:)
      type(error_type), allocatable, intent(out) :: error

      call get_solvent_charges(self%solvent%id, environment, model, charges, error)
   end subroutine get_system_charges

   !> Get atom-centered multipoles for a gas, solvent, or conductor environment
   subroutine get_system_multipoles(self, environment, model, multipoles, error)
      class(solvent_system_type), intent(in) :: self
      character(len=*), intent(in) :: environment
      character(len=*), intent(in) :: model
      type(solvent_multipole_data_type), intent(out) :: multipoles
      type(error_type), allocatable, intent(out) :: error

      call get_solvent_multipoles(self%solvent%id, environment, model, multipoles, error)
   end subroutine get_system_multipoles


   !> Bulk descriptors of a solvent named by any of its aliases
   !>
   !> @param[in]  alias       solvent name or alias; case and surrounding blanks ignored
   !> @param[out] descriptors table values and derived bulk properties
   !> @param[out] error       blank or unknown alias
   subroutine get_solvent_descriptors(alias, descriptors, error)
      !> Solvent name or alias
      character(len=*), intent(in) :: alias
      !> Table values and derived bulk properties
      type(solvent_descriptors_type), intent(out) :: descriptors
      !> Blank or unknown alias
      type(error_type), allocatable, intent(out) :: error
      !> Matched table ID
      integer :: solvent_id

      call get_solvent_id(alias, solvent_id, error)
      if (allocated(error)) return
      call get_solvent_descriptors_by_id(solvent_id, descriptors, error)
   end subroutine get_solvent_descriptors

   !> Bulk descriptors of a solvent table ID
   !>
   !> @param[in]  solvent_id  table ID
   !> @param[out] descriptors table values and derived bulk properties
   !> @param[out] error       unknown ID
   subroutine get_solvent_descriptors_by_id(solvent_id, descriptors, error)
      !> Table ID
      integer, intent(in) :: solvent_id
      !> Table values and derived bulk properties
      type(solvent_descriptors_type), intent(out) :: descriptors
      !> Unknown ID
      type(error_type), allocatable, intent(out) :: error
      !> Solvent geometry, only needed for the molar mass
      type(structure_type), allocatable :: mol

      call load_solvent_entry(solvent_id, descriptors, mol, error)
   end subroutine get_solvent_descriptors_by_id

   !> Read one table entry, build its geometry and derive the bulk properties
   !>
   !> Shared by the descriptor getters and new_solvent_system; reads only
   !> constant tables, so concurrent calls are safe. Without a stored
   !> geometry the values derived from it are NaN and mol stays unallocated
   !>
   !> @param[in]  solvent_id  table ID
   !> @param[out] descriptors table values and derived bulk properties
   !> @param[out] mol         solvent geometry in bohr, if stored
   !> @param[out] error       unknown ID
   subroutine load_solvent_entry(solvent_id, descriptors, mol, error)
      !> Table ID
      integer, intent(in) :: solvent_id
      !> Table values and derived bulk properties
      type(solvent_descriptors_type), intent(out) :: descriptors
      !> Solvent geometry in bohr, if stored
      type(structure_type), allocatable, intent(out) :: mol
      !> Unknown ID
      type(error_type), allocatable, intent(out) :: error
      !> Iterables
      integer :: i
      !> Whether solvent_id matched an entry in the table
      logical :: found
      !> Per-atom mass accumulated into the molar mass
      real(wp) :: atomic_mass
      !> Marker of values the table lacks
      real(wp) :: missing
      character(len=64) :: id_msg

      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho

      !> Get basic solvent information
      include "solvents.inc"
      found = .false.
      do i = 1, max_solvents
         if (solvent_id == id_list(i)) then
            found = .true.
            descriptors%id = id_list(i)
            descriptors%name = trim(to_lower(name_list(i)))
            descriptors%epsilon = eps(i)
            descriptors%refractive_index = refr(i)
            descriptors%abraham_alpha = A(i)
            descriptors%abraham_beta = B(i)
            descriptors%surface_tension_si = g(i)*0.001_wp
            descriptors%mass_density_si = rho(i)
            exit
         end if
      end do

      ! Error for unknown solvent ID
      if (.not. found) then
         write (id_msg, "(a,i0)") "Unknown solvent ID: ", solvent_id
         call fatal_error(error, trim(id_msg))
         return
      end if

      ! Geometry-free values in atomic units
      descriptors%surface_tension_au = descriptors%surface_tension_si*(Bohr_radius**2)/Hartree_energy
      descriptors%mass_density_au = descriptors%mass_density_si*(Bohr_radius**3)/atomic_unit_of_mass

      call get_solvent_geometry(solvent_id, mol, error)
      if (allocated(error)) then
         ! No stored geometry: the values derived from it are missing
         deallocate (error)
         if (allocated(mol)) deallocate (mol)
         missing = ieee_value(missing, ieee_quiet_nan)
         descriptors%molar_mass_si = missing
         descriptors%molecular_mass_au = missing
         descriptors%number_density_si = missing
         descriptors%number_density_au = missing
         descriptors%molecular_volume_si = missing
         descriptors%molecular_volume_au = missing
         return
      end if

      ! Convert coordinates to atomic units
      mol%xyz = mol%xyz*aatoau

      !> Compute the atomic mass of the solvent
      descriptors%molar_mass_si = 0.0_wp
      do i = 1, mol%nat
         call get_mass(mol%num(mol%id(i)), atomic_mass, error)
         if (allocated(error)) return
         descriptors%molar_mass_si = descriptors%molar_mass_si + atomic_mass*0.001_wp
      end do

      ! Mass: kg/mol -> me
      descriptors%molecular_mass_au = descriptors%molar_mass_si/Avogadro_constant/atomic_unit_of_mass

      ! Density: kg/m^3 to 1/m^3
      descriptors%number_density_si = Avogadro_constant*descriptors%mass_density_si/descriptors%molar_mass_si

      ! Solvent number density: 1/m^3 -> 1/bohr^3
      descriptors%number_density_au = descriptors%number_density_si*Bohr_radius**3

      ! Solvent molecular volume
      descriptors%molecular_volume_si = descriptors%molar_mass_si/ &
                                        descriptors%mass_density_si/Avogadro_constant

      ! Solvent molecular volume: m^3 -> bohr^3
      descriptors%molecular_volume_au = descriptors%molecular_volume_si/(Bohr_radius**3)

   end subroutine load_solvent_entry

   !> Initialize a solvent system from solvent data
   subroutine new_solvent_system( &
      self, &
      solvent_id, &
      temperature, &
      pressure_si, &
      error &
      )

      !> Solvent system
      class(solvent_system_type), intent(out) :: self

      !> Solvent ID
      integer, intent(in) :: solvent_id

      !> Optional temperature and pressure (in SI units)
      real(wp), intent(in), optional :: temperature, pressure_si

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=64) :: id_msg

      ! Default values
      if (.not. present(temperature)) then
         self%temperature = 298.15_wp ! Default temperature: 25 degrees Celsius
      else
         if (temperature <= 0.0_wp) then
            call fatal_error(error, "Temperature must be positive.")
            return
         else
            self%temperature = temperature
         end if
      end if

      if (.not. present(pressure_si)) then
         self%pressure_si = 101325.0_wp ! Default pressure: 1 atm in Pa
      else
         if (pressure_si < 0.0_wp) then
            call fatal_error(error, "Pressure must be non-negative.")
            return
         else
            self%pressure_si = pressure_si
         end if
      end if

      call load_solvent_entry(solvent_id, self%solvent, self%solv_mol, error)
      if (allocated(error)) return
      if (.not. allocated(self%solv_mol)) then
         write (id_msg, "(a,i0,a)") " (ID ", self%solvent%id, ")"
         call fatal_error(error, "No geometry available for solvent '"//self%solvent%name//"'"//trim(id_msg))
         return
      end if

      ! Pressure: Pa (kg/m/s**2) -> Eh/bohr^3
      self%pressure_au = self%pressure_si/atomic_unit_of_mass*Bohr_radius*2.4188843E-17_wp**2

   end subroutine new_solvent_system

   !> Print the solvent system properties
   subroutine print_solvent_system(self)
      class(solvent_system_type), intent(in) :: self
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit, col_value=30, indent_step=2, fmt_len=16)

      call pp%blank()
      call pp%push("System properties:")
      call pp%kv("Temperature", self%temperature, "K")
      call pp%kv2("Pressure", self%pressure_si, "Pa", self%pressure_au, "au")
      call pp%pop()

      call pp%blank()
      if (allocated(self%solv_mol)) then
         call write_solvent_properties(pp, self%solvent, self%solv_mol%nat)
      else
         call write_solvent_properties(pp, self%solvent)
      end if
      call pp%blank()

   end subroutine print_solvent_system

   !> Print solvent descriptors as a "Solvent properties" block
   !>
   !> @param[in] self descriptors to print
   !> @param[in] unit output unit, standard output by default
   subroutine print_solvent_descriptors(self, unit)
      !> Descriptors to print
      class(solvent_descriptors_type), intent(in) :: self
      !> Output unit, standard output by default
      integer, intent(in), optional :: unit
      !> Block layout shared with print_solvent_system
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit, col_value=30, indent_step=2, fmt_len=16)
      if (present(unit)) call pp%set_unit(unit)
      call write_solvent_properties(pp, self)
      flush (pp%iu)
   end subroutine print_solvent_descriptors

   !> Write the "Solvent properties" block
   !>
   !> @param[inout] pp          pretty printer, unit and layout
   !> @param[in]    descriptors descriptors to print
   !> @param[in]    nat         atoms of the solvent geometry, omitted when not built
   subroutine write_solvent_properties(pp, descriptors, nat)
      !> Pretty printer, unit and layout
      type(prettyprinter), intent(inout) :: pp
      !> Descriptors to print
      type(solvent_descriptors_type), intent(in) :: descriptors
      !> Atoms of the solvent geometry
      integer, intent(in), optional :: nat

      call pp%push("Solvent properties:")
      if (allocated(descriptors%name)) call pp%kv("Name", trim(descriptors%name))
      call pp%kv("ID", descriptors%id)
      if (present(nat)) call pp%kv("Number of atoms", nat)
      call pp%kv2("Mass", descriptors%molar_mass_si, "kg/mol", descriptors%molecular_mass_au, "me")
      call pp%kv2("Mass density", descriptors%mass_density_si, "kg/m**3", &
                  descriptors%mass_density_au, "me/bohr**3")
      call pp%kv2("Number density", descriptors%number_density_si, "1/m**3", &
                  descriptors%number_density_au, "1/bohr**3")
      call pp%kv2("Molecular volume", descriptors%molecular_volume_si, "m**3", &
                  descriptors%molecular_volume_au, "bohr**3")
      call pp%kv2("Surface tension", descriptors%surface_tension_si, "N/m", &
                  descriptors%surface_tension_au, "Eh/bohr**2")
      call pp%kv("Rel. permitivity", descriptors%epsilon, "eps/eps0")
      call pp%kv("Refractive index", descriptors%refractive_index, "c/c0")
      call pp%kv("HB acidity", descriptors%abraham_alpha)
      call pp%kv("HB basicity", descriptors%abraham_beta)
      call pp%pop()
   end subroutine write_solvent_properties

end module moist_data_solvents
