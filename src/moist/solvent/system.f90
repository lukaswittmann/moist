
!> Solvent at a thermodynamic state, with its geometry
module moist_solvent_system
   use mctc_env, only: wp
   use mctc_env_error, only: error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_codata2018, only: Bohr_radius, atomic_unit_of_mass
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_data_solvents, only: solvent_multipole_data_type, get_solvent_charges, &
      & get_solvent_multipoles
   use moist_solvent_descriptors, only: solvent_descriptors_type, load_solvent_entry, &
      & write_solvent_properties
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)

   public :: solvent_system_type, new_solvent_system

   private

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

contains

   !> Initialize a solvent system from solvent data
   !>
   !> @param[out] self        solvent system
   !> @param[in]  solvent_id  table ID
   !> @param[in]  temperature temperature in K, 298.15 K by default
   !> @param[in]  pressure_si pressure in Pa, 1 atm by default
   !> @param[out] error       invalid state, unknown ID or no stored geometry
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

      !> Optional temperature
      real(wp), intent(in), optional :: temperature

      !> Optional pressure
      real(wp), intent(in), optional :: pressure_si

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

   !> Get atomic charges for a gas, solvent, or conductor environment
   !>
   !> @param[in]  self        solvent system
   !> @param[in]  environment gas, solvent or conductor
   !> @param[in]  model       charge model
   !> @param[out] charges     atomic charges in e
   !> @param[out] error       unknown environment or model, or no data
   subroutine get_system_charges(self, environment, model, charges, error)
      !> Solvent system
      class(solvent_system_type), intent(in) :: self
      !> Environment
      character(len=*), intent(in) :: environment
      !> Charge model
      character(len=*), intent(in) :: model
      !> Atomic charges
      real(wp), allocatable, intent(out) :: charges(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call get_solvent_charges(self%solvent%id, environment, model, charges, error)
   end subroutine get_system_charges

   !> Get atom-centered multipoles for a gas, solvent, or conductor environment
   !>
   !> @param[in]  self        solvent system
   !> @param[in]  environment gas, solvent or conductor
   !> @param[in]  model       multipole model
   !> @param[out] multipoles  multipoles in e*bohr^rank
   !> @param[out] error       unknown environment or model, or no data
   subroutine get_system_multipoles(self, environment, model, multipoles, error)
      !> Solvent system
      class(solvent_system_type), intent(in) :: self
      !> Environment
      character(len=*), intent(in) :: environment
      !> Multipole model
      character(len=*), intent(in) :: model
      !> Atom-centered multipoles
      type(solvent_multipole_data_type), intent(out) :: multipoles
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call get_solvent_multipoles(self%solvent%id, environment, model, multipoles, error)
   end subroutine get_system_multipoles

   !> Print the solvent system properties
   !>
   !> @param[in] self solvent system
   subroutine print_solvent_system(self)
      !> Solvent system
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

end module moist_solvent_system
