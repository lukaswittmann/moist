
module moist_data_solvents
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new_structure
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_io_convert, only: autokcal, aatoau
   use mctc_io_codata2018, only: Avogadro_constant, Bohr_radius
   use mctc_io_codata2018, only: Hartree_energy, atomic_unit_of_mass
   use mctc_env_error, only: error_type, fatal_error
   use moist_data_mass, only: get_mass
   use mctc_io_utils, only: to_lower
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)

   integer, parameter, public :: max_solvents = 180

   public :: solvation_system_type, new_solvation_system
   public :: get_solvent_id

   private

   !> Solvent and solute properties and molecular structures
   type :: solvation_system_type

      integer :: solvent_id
      character(:), allocatable  :: solvent_name

      !> System info
      real(wp) :: temperature                 ! Temperature in Kelvin
      real(wp) :: pressure_si                 ! Pressure in Pa
      real(wp) :: pressure_au                 ! Pressure in atomic units

      !> Solvent properties
      real(wp) :: solvent_epsilon             ! Dielectric constant
      real(wp) :: solvent_refractive_index    ! Refractive index
      real(wp) :: solvent_alpha               ! Abrahams HB acidity
      real(wp) :: solvent_beta                ! Abrahams HB basicity
      real(wp) :: solvent_surface_tension_si  ! Surface tension in SI units (N/m)
      real(wp) :: solvent_surface_tension_au  ! Surface tension in atomic units
      real(wp) :: solvent_mass_density_si     ! Mass density in kg/m^3
      real(wp) :: solvent_mass_density_au     ! Density in atomic units
      real(wp) :: solvent_number_density_si   ! Solvent number density in 1/m^3
      real(wp) :: solvent_number_density_au   ! Solvent number density in atomic units (1/bohr^3)
      real(wp) :: solvent_molecular_volume_si ! Volume per solvent molecule in m^3
      real(wp) :: solvent_molecular_volume_au ! Volume per solvent molecule in atomic units
      real(wp) :: solvent_molar_mass_si       ! Molar mass of solvent in kg/mol
      real(wp) :: solvent_mass_au             ! Mass of solvent in atomic units (AU)

      !> Solute properties (that do *not* depend on the geometry)
      real(wp) :: solute_molar_mass_si ! Molar mass of solute in kg/mol
      real(wp) :: solute_mass_au       ! Mass of solute in atomic units (AU)

      !> Solvent geometry
      type(structure_type), allocatable :: solv_mol ! Geometry of the solvent molecule
      type(structure_type), allocatable :: solu_mol ! Geometry of the solute molecule

   contains

      procedure :: print => print_solvation_system
      procedure :: update => add_solute_properties

   end type solvation_system_type

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


   !> Initialize a solvation system from solvent data
   subroutine new_solvation_system( &
      self, &
      solvent_id, &
      temperature, &
      pressure_si, &
      error &
      )

      !> Solvation system
      class(solvation_system_type), intent(out) :: self

      !> Solvent ID
      integer, intent(in) :: solvent_id

      !> Optional temperature and pressure (in SI units)
      real(wp), intent(in), optional :: temperature, pressure_si

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Iterables
      integer :: i

      !> Whether solvent_id matched an entry in the table
      logical :: found

      !> Per-atom mass accumulated into the molar mass
      real(wp) :: atomic_mass
      character(len=64) :: id_msg

      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho

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

      !> Get basic solvent information
      include "solvents.inc"
      found = .false.
      do i = 1, max_solvents
         if (solvent_id == id_list(i)) then
            found = .true.
            self%solvent_id = id_list(i)
            self%solvent_name = trim(to_lower(name_list(i)))
            self%solvent_epsilon = eps(i)
            self%solvent_refractive_index = refr(i)
            self%solvent_alpha = A(i)
            self%solvent_beta = B(i)
            self%solvent_surface_tension_si = g(i)*0.001_wp
            self%solvent_mass_density_si = rho(i)
            exit
         end if
      end do

      ! Error for unknown solvent ID
      if (.not. found) then
         write (id_msg, "(a,i0)") "Unknown solvent ID: ", solvent_id
         call fatal_error(error, trim(id_msg))
         return
      end if

      allocate (self%solv_mol)
      call get_solvent_geometry(self%solvent_id, self%solv_mol, error)
      if (allocated(error)) then
         write (id_msg, "(a,i0,a)") " (ID ", self%solvent_id, ")"
         call fatal_error(error, "No geometry available for solvent '"//self%solvent_name//"'"//trim(id_msg))
         return
      end if

      ! Convert coordinates to atomic units
      self%solv_mol%xyz = self%solv_mol%xyz*aatoau

      !> Compute the atomic mass of the solvent
      self%solvent_molar_mass_si = 0.0_wp
      do i = 1, self%solv_mol%nat
         call get_mass(self%solv_mol%num(self%solv_mol%id(i)), atomic_mass, error)
         if (allocated(error)) return
         self%solvent_molar_mass_si = self%solvent_molar_mass_si + atomic_mass*0.001_wp
      end do

      ! Pressure: Pa (kg/m/s**2) -> Eh/bohr^3
      self%pressure_au = self%pressure_si/atomic_unit_of_mass*Bohr_radius*2.4188843E-17_wp**2

      ! Mass: kg/mol -> me
      self%solvent_mass_au = self%solvent_molar_mass_si/Avogadro_constant/atomic_unit_of_mass

      ! Surface tension: N/m (=J/m**2) -> AU
      self%solvent_surface_tension_au = self%solvent_surface_tension_si*(Bohr_radius**2)/Hartree_energy

      ! Density: kg/m^3 to 1/m^3
      self%solvent_number_density_si = Avogadro_constant*self%solvent_mass_density_si/self%solvent_molar_mass_si

      ! Solvent number density: 1/m^3 -> 1/bohr^3
      self%solvent_number_density_au = self%solvent_number_density_si*Bohr_radius**3

      ! Solvent molecular volume
      self%solvent_molecular_volume_si = self%solvent_molar_mass_si/ &
                                         self%solvent_mass_density_si/Avogadro_constant

      ! Solvent mass density: kg/m^3 -> me/bohr^3
      self%solvent_mass_density_au = self%solvent_mass_density_si*(Bohr_radius**3)/atomic_unit_of_mass

      ! Solvent molecular volume: m^3 -> bohr^3
      self%solvent_molecular_volume_au = self%solvent_molecular_volume_si/(Bohr_radius**3)

   end subroutine new_solvation_system

   !> Subroutine that adds solute properties to the solvation system
   subroutine add_solute_properties(self, solu_mol, error)

      !> Solvation system
      class(solvation_system_type), intent(inout) :: self

      !> Solute molecule geometry
      type(structure_type), intent(in) :: solu_mol

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: i
      real(wp) :: atomic_mass

      !> Check if the solute molecule is allocated
      if (.not. allocated(self%solu_mol)) then
         allocate (self%solu_mol)
      end if

      !> Copy the solute molecule geometry
      self%solu_mol = solu_mol

      !> Compute the atomic mass of the solute
      self%solute_molar_mass_si = 0.0_wp
      do i = 1, self%solu_mol%nat
         call get_mass(self%solu_mol%num(self%solu_mol%id(i)), atomic_mass, error)
         if (allocated(error)) return
         self%solute_molar_mass_si = self%solute_molar_mass_si + atomic_mass*0.001_wp
      end do

      !> Convert to atomic units (AU)
      self%solute_mass_au = self%solute_molar_mass_si/Avogadro_constant/atomic_unit_of_mass

   end subroutine add_solute_properties

   !> Print the solvation system properties
   subroutine print_solvation_system(self)
      class(solvation_system_type), intent(in) :: self
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit, col_value=30, indent_step=2, fmt_len=16)

      call pp%blank()
      call pp%push("System properties:")
      call pp%kv("Temperature", self%temperature, "K")
      call pp%kv2("Pressure", self%pressure_si, "Pa", self%pressure_au, "au")
      call pp%pop()

      call pp%blank()
      call pp%push("Solvent properties:")
      call pp%kv("Name", trim(self%solvent_name))
      call pp%kv("ID", self%solvent_id)
      if (allocated(self%solv_mol)) then
         call pp%kv("Number of atoms", self%solv_mol%nat)
      end if
      call pp%kv2("Mass", self%solvent_molar_mass_si, "kg/mol", self%solvent_mass_au, "me")
      call pp%kv2("Mass density", self%solvent_mass_density_si, "kg/m**3", &
                  self%solvent_mass_density_au, "me/bohr**3")
      call pp%kv2("Number density", self%solvent_number_density_si, "1/m**3", &
                  self%solvent_number_density_au, "1/bohr**3")
      call pp%kv2("Molecular volume", self%solvent_molecular_volume_si, "m**3", &
                  self%solvent_molecular_volume_au, "bohr**3")
      call pp%kv2("Surface tension", self%solvent_surface_tension_si, "N/m", &
                  self%solvent_surface_tension_au, "Eh/bohr**2")
      call pp%kv("Rel. permitivity", self%solvent_epsilon, "eps/eps0")
      call pp%kv("Refractive index", self%solvent_refractive_index, "c/c0")
      call pp%kv("HB acidity", self%solvent_alpha)
      call pp%kv("HB basicity", self%solvent_beta)
      call pp%pop()

      if (allocated(self%solu_mol)) then
         call pp%blank()
         call pp%push("Solute properties:")
         call pp%kv("Number of atoms", self%solu_mol%nat)
         call pp%kv2("Mass", self%solute_molar_mass_si, "kg/mol", self%solute_mass_au, "me")
         call pp%pop()
      end if
      call pp%blank()

   end subroutine print_solvation_system

end module moist_data_solvents
