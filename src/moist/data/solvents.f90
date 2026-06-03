
module moist_data_solvents
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new_structure
   use iso_fortran_env, only: output_unit
   use mctc_io_convert, only: autokcal, aatoau
   use mctc_io_codata2018, only: Avogadro_constant, Bohr_radius
   use mctc_io_codata2018, only: Hartree_energy, atomic_unit_of_mass
   use mctc_env_error, only: error_type, fatal_error
   use moist_data_mass, only: get_mass
   use mctc_io_utils, only: to_lower
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_radii_static, only: static_radius_type, new_bondi_radii
   implicit none

   integer, parameter, public :: max_solvents = 180

   public :: solvation_system_parameters, new_solvation_system_parameters
   public :: get_solvent_id, get_solvent_for_alpb

   private

   type :: solvation_system_parameters

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
      real(wp) :: solvent_mass_density_si          ! Density in g/cm^3
      real(wp) :: solvent_mass_density_au          ! Density in
      real(wp) :: solvent_number_density_si   ! Solvent number density in mol/m^3.
      real(wp) :: solvent_number_density_au   ! Solvent number density in
      real(wp) :: solvent_molecular_volume_si ! Volume per solvent molecule in m^3.
      real(wp) :: solvent_molecular_volume_au ! Volume per solvent molecule in
      real(wp) :: solvent_molar_mass_si       ! Molar mass of solvent in kg/mol.
      real(wp) :: solvent_mass_au             ! Mass of solvent in atomic units (AU).
      real(wp) :: solvent_packing_fraction    ! Packing fraction of the solvent

      !> Calculated propertiess
      real(wp) :: solvent_self_solvation_energy_au ! Self-solvation energy in atomic units
      real(wp) :: solvent_self_solvation_energy_kcal ! Self-solvation energy in kcal/mol
      real(wp) :: solvent_cavity_volume_au ! Volume of the solvent cavity in atomic units
      real(wp) :: solvent_cavity_area_au   ! Area of the solvent cavity./

      !> Solute properties (that do *not* depend on the geometry)
      real(wp) :: solute_molar_mass_si        ! Molar mass of solute in kg/mol.
      real(wp) :: solute_mass_au              ! Mass of solute in atomic units (AU).

      !> Solvent geometry
      type(structure_type), allocatable :: solv_mol ! Geometry of the solvent molecule
      type(structure_type), allocatable :: solu_mol ! Geometry of the solute molecule

   contains

      procedure :: print => print_solvation_system_parameters
      procedure :: update => add_solute_properties

   end type solvation_system_parameters

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
      character(:), allocatable :: name

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho, eta

      !> Get basic solvent information
      include "solvents.inc"
      do i = 1, max_solvents
         do j = 1, 10
            if (trim(to_lower(alias)) == trim(alias_list(j, i))) then
               solvent_id = id_list(i)
               name = to_lower(name_list(i))
               return
            end if
         end do
      end do

      ! If we reach here, the alias was not found
      call fatal_error(error, message="Unknown solvent: "//trim(alias), stat=1)

   end subroutine get_solvent_id

   !> Include subroutine to get the solvent geometry
   !! subroutine get_solvent_geometry(solvent_id,mol)
   !!    integer, intent(in) :: solvent_id
   !!    type(structure_type), intent(out) :: mol
   include "solventgeometries.inc"

   subroutine get_solvent_for_alpb(solvent_id, epsilon, solvent_name, error)
      integer, intent(in) :: solvent_id
      real(wp), intent(out) :: epsilon
      character(:), allocatable, intent(out) :: solvent_name
      type(error_type), allocatable, intent(out) :: error

      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)
      character(:), allocatable :: name

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho

      include "solvents.inc"

      integer :: i

      do i = 1, max_solvents
         if (solvent_id == id_list(i)) then
            epsilon = eps(i)
            solvent_name = trim(to_lower(name_list(i)))
            return
         end if
      end do

      ! If we reach here, the solvent ID was not found
      call fatal_error(error, message="Unknown solvent ID", stat=1)

   end subroutine get_solvent_for_alpb

   !> Initialize system parameters for a new solvent
   subroutine new_solvation_system_parameters( &
      self, &
      solvent_id, &
      temperature, &
      pressure_si, &
      error &
      )

      !> Solvation system parameters
      class(solvation_system_parameters), intent(out) :: self

      !> Solvent ID
      integer, intent(in) :: solvent_id

      !> Optional temperature and pressure (in SI units)
      real(wp), intent(in), optional :: temperature, pressure_si

      !> Error handling
      type(error_type), allocatable, intent(out), optional :: error
      type(error_type), allocatable :: local_error
      !> Iterables
      integer :: i

      character(len=64) :: name_list(max_solvents)
      character(len=64) :: alias_list(10, max_solvents)
      character(:), allocatable :: name

      ! NumSA and radius model for packing fraction calculation
      type(cavity_type_iswig), allocatable :: cavity
      type(static_radius_type) :: radii

      integer, dimension(max_solvents) :: id_list
      real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho

      ! Default values
      if (.not. present(temperature)) then
         self%temperature = 298.15_wp ! Default temperature: 25 degrees Celsius
      else
         if (temperature <= 0.0_wp) then
            call fatal_error(local_error, "Temperature must be positive.")
            if (present(error)) then; error = local_error; end if
            return
         else
            self%temperature = temperature
         end if
      end if

      if (.not. present(pressure_si)) then
         self%pressure_si = 101325.0_wp ! Default pressure: 1 atm in Pa
      else
         if (pressure_si < 0.0_wp) then
            call fatal_error(local_error, "Pressure must be non-negative.")
            if (present(error)) then; error = local_error; end if
            return
         else
            self%pressure_si = pressure_si
         end if
      end if

      ! Self-solvation energy (default values)
      self%solvent_self_solvation_energy_au = 0.0_wp

      !> Get basic solvent informationget_solvation_system_parameters
      include "solvents.inc"
      do i = 1, max_solvents
         if (solvent_id == id_list(i)) then
            self%solvent_id = id_list(i)
            self%solvent_name = trim(to_lower(name_list(i)))
            self%solvent_epsilon = eps(i)
            self%solvent_refractive_index = refr(i)
            self%solvent_alpha = A(i)
            self%solvent_beta = B(i)
            self%solvent_surface_tension_si = g(i)*0.001_wp
            self%solvent_mass_density_si = rho(i)
         end if
      end do

      allocate (self%solv_mol)
      call get_solvent_geometry(self%solvent_id, self%solv_mol, local_error)
      if (allocated(local_error)) then
         if (present(error)) then; error = local_error; end if
         return
      end if

      ! Convert coordinates to atomic units
      self%solv_mol%xyz = self%solv_mol%xyz*aatoau

      !> Compute the atomic mass of the solvent
      self%solvent_molar_mass_si = 0.0_wp
      do i = 1, self%solv_mol%nat
         self%solvent_molar_mass_si = self%solvent_molar_mass_si &
                                      + get_mass(self%solv_mol%num(self%solv_mol%id(i)))*0.001_wp
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
      self%solvent_mass_density_au = self%solvent_mass_density_si*(Bohr_radius**3) &
                                     /Avogadro_constant/atomic_unit_of_mass

      ! Convert solvent molecular volume to atomic units (m^3/mol)
      self%solvent_molecular_volume_au = self%solvent_molecular_volume_si/(Bohr_radius**3)

      ! Convert self-solvation energy to kcal/mol
      self%solvent_self_solvation_energy_kcal = self%solvent_self_solvation_energy_au &
                                                *autokcal

      ! Compute packing fraction
      ! FIXME: This is a layer violation #64; but otherwise its much more ugly elsewhere (to be done)
      call new_bondi_radii(radii)
      allocate (cavity)
      call new_cavity_iswig(cavity, nleb=110, &
                            radius_model=radii, error=error)
      if (allocated(error)) return
      call cavity%update(self%solv_mol, error=error)
      if (allocated(error)) return
      self%solvent_packing_fraction = self%solvent_number_density_au*cavity%total_volume

   end subroutine new_solvation_system_parameters

   !> Subroutine that adds solute properties to the solvation system
   subroutine add_solute_properties(self, solu_mol)

      !> Solvation system parameters
      class(solvation_system_parameters), intent(inout) :: self

      !> Solute molecule geometry
      type(structure_type), intent(in) :: solu_mol

      integer :: i

      !> Check if the solute molecule is allocated
      if (.not. allocated(self%solu_mol)) then
         allocate (self%solu_mol)
      end if

      !> Copy the solute molecule geometry
      self%solu_mol = solu_mol

      !> Compute the atomic mass of the solute
      self%solute_molar_mass_si = 0.0_wp
      do i = 1, self%solu_mol%nat
         self%solute_molar_mass_si = self%solute_molar_mass_si &
                                     + get_mass(self%solu_mol%num(self%solu_mol%id(i)))*0.001_wp
      end do

      !> Convert to atomic units (AU)
      self%solute_mass_au = self%solute_molar_mass_si/Avogadro_constant/atomic_unit_of_mass

   end subroutine add_solute_properties

   !> Print the solvation system parameters
   subroutine print_solvation_system_parameters(self)
      class(solvation_system_parameters), intent(in) :: self
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit, col_value=30, indent_step=2, fmt_len=16)

      call pp%blank()
      call pp%push('System properties:')
      call pp%kv('Temperature', self%temperature, 'K')
      call pp%kv2('Pressure', self%pressure_si, 'Pa', self%pressure_au, 'au')
      call pp%pop()

      call pp%blank()
      call pp%push('Solvent properties:')
      call pp%kv('Name', trim(self%solvent_name))
      call pp%kv('ID', self%solvent_id)
      if (allocated(self%solv_mol)) then
         call pp%kv('Number of atoms', self%solv_mol%nat)
      end if
      call pp%kv2('Mass', self%solvent_molar_mass_si, 'kg/mol', self%solvent_mass_au, 'me')
      call pp%kv2('Mass density', self%solvent_mass_density_si, 'kg/m**3', &
                  self%solvent_mass_density_au, 'me/bohr**3')
      call pp%kv2('Number density', self%solvent_number_density_si, '1/m**3', &
                  self%solvent_number_density_au, '1/bohr**3')
      call pp%kv2('Molecular volume', self%solvent_molecular_volume_si, 'm**3', &
                  self%solvent_molecular_volume_au, 'bohr**3')
      call pp%kv('Packing fraction', self%solvent_packing_fraction)
      call pp%kv2('Surface tension', self%solvent_surface_tension_si, 'N/m', &
                  self%solvent_surface_tension_au, 'Eh/bohr**2')
      call pp%kv('Rel. permitivity', self%solvent_epsilon, 'eps/eps0')
      call pp%kv('Refractive index', self%solvent_refractive_index, 'c/c0')
      call pp%kv('HB acidity', self%solvent_alpha)
      call pp%kv('HB basicity', self%solvent_beta)
      call pp%pop()

      if (allocated(self%solu_mol)) then
         call pp%blank()
         call pp%push('Solute properties:')
         call pp%kv('Number of atoms', self%solu_mol%nat)
         call pp%kv2('Mass', self%solute_molar_mass_si, 'kg/mol', self%solute_mass_au, 'me')
         call pp%pop()
      end if
      call pp%blank()

   end subroutine print_solvation_system_parameters

end module moist_data_solvents
