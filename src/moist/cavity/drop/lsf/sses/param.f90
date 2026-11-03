!> Parameter container for the smooth SES-like optimized moment level set
module moist_cavity_drop_lsf_sses_param
   use, intrinsic :: iso_fortran_env, only: output_unit
   use mctc_env_accuracy, only: wp
   use mctc_io_convert, only: aatoau
   use moist_model_parameters, only: moist_model_parameters_type
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)
   private
   public :: moist_cavity_drop_lsf_sses_param_type, temperature_length

   !> Fixed length r0 of the Boltzmann temperature 2 eps r0 in Bohr, 1.5 Angstrom
   real(wp), parameter :: temperature_length = 1.5_wp*aatoau

   !> Physical SSES parameters in atomic units
   type, extends(moist_model_parameters_type) :: moist_cavity_drop_lsf_sses_param_type
      !> Solvent probe radius in Bohr
      real(wp) :: probe_radius = 2.7_wp
      !> Smoothing length in Bohr
      real(wp) :: smoothing_length = 0.5_wp
      !> Positive field multiplier for DROP critical-gradient switching
      real(wp) :: scale = 1.0_wp
   contains
      !> Restore compiled parameter defaults
      procedure :: init_defaults
      !> Register fields for serialization
      procedure :: register_entries
      !> Override selected physical parameters
      procedure :: new
      !> Print physical parameters
      procedure :: print
      !> Boltzmann temperature of the weight entropy
      procedure :: temperature
   end type moist_cavity_drop_lsf_sses_param_type

contains

   !> Restore compiled defaults
   !>
   !> @param[in,out] self Parameter values
   subroutine init_defaults(self)
      !> Parameter values
      class(moist_cavity_drop_lsf_sses_param_type), intent(inout) :: self
      type(moist_cavity_drop_lsf_sses_param_type) :: defaults

      self%probe_radius = defaults%probe_radius
      self%smoothing_length = defaults%smoothing_length
      self%scale = defaults%scale
   end subroutine init_defaults

   !> Register parameter fields
   !>
   !> @param[in,out] self Parameter values
   subroutine register_entries(self)
      !> Parameter values
      class(moist_cavity_drop_lsf_sses_param_type), intent(inout), target :: self

      call self%register_real_scalar("probe_radius", self%probe_radius)
      call self%register_real_scalar("smoothing_length", self%smoothing_length)
      call self%register_real_scalar("scale", self%scale)
   end subroutine register_entries

   !> Override selected physical parameters
   !>
   !> @param[in,out] self Parameter values
   !> @param[in] probe_radius Probe radius in Bohr
   !> @param[in] smoothing_length Smoothing length in Bohr
   !> @param[in] scale Positive field multiplier
   subroutine new(self, probe_radius, smoothing_length, scale)
      !> Parameter values
      class(moist_cavity_drop_lsf_sses_param_type), intent(inout) :: self
      !> Probe radius in Bohr
      real(wp), intent(in), optional :: probe_radius
      !> Smoothing length in Bohr
      real(wp), intent(in), optional :: smoothing_length
      !> Positive field multiplier
      real(wp), intent(in), optional :: scale

      if (present(probe_radius)) self%probe_radius = probe_radius
      if (present(smoothing_length)) self%smoothing_length = smoothing_length
      if (present(scale)) self%scale = scale
   end subroutine new

   !> Print physical parameters
   !>
   !> @param[in] self Parameter values
   !> @param[in] unit Output unit override
   subroutine print(self, unit)
      !> Parameter values
      class(moist_cavity_drop_lsf_sses_param_type), intent(in) :: self
      !> Output unit override
      integer, intent(in), optional :: unit
      type(prettyprinter) :: pp
      integer :: iu

      iu = output_unit
      if (present(unit)) iu = unit
      pp = new_prettyprinter(unit=iu)
      call pp%push("Implicit surface (SSES):")
      call pp%kv("Probe radius (Bohr)", self%probe_radius)
      call pp%kv("Smoothing length (Bohr)", self%smoothing_length)
      call pp%kv("Field scale", self%scale)
      call pp%pop()
      call pp%blank()
   end subroutine print

   !> Boltzmann temperature 2 eps r0, one constant for every atom and point
   !>
   !> The fixed length r0 converts the smoothing length into the squared length of
   !> the power-distance exponents, so eps remains the only smoothing control
   !>
   !> @param[in] self Parameter values
   pure function temperature(self) result(tau)
      !> Parameter values
      class(moist_cavity_drop_lsf_sses_param_type), intent(in) :: self
      !> Temperature in Bohr^2
      real(wp) :: tau

      tau = 2.0_wp*self%smoothing_length*temperature_length
   end function temperature

end module moist_cavity_drop_lsf_sses_param
