!> Parameter container shared by the two isodensity LSFs
!>
!> Both isodensity variants describe the same level set
!>
!>    S(r) = scale * (rho_iso - rho(r))
!>
!> and differ only in where the density comes from; so they use the same parameters
module moist_cavity_drop_lsf_isodensity_param
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)
   private

   public :: moist_cavity_drop_lsf_isodensity_param_type
   public :: isodensity_exclusion_radius

   !> Isodensity level set function parameters
   type, extends(moist_model_parameters_type) :: moist_cavity_drop_lsf_isodensity_param_type
      !> Density isovalue defining the surface, in Bohr^-3
      real(wp) :: rho_iso = 1.0E-3_wp
      !> Constant multiplier applied to the level set value and derivatives
      real(wp) :: scale = 1.0_wp / 1.0E-3_wp

      !* ------------------- Surface-free exclusion certificate ---------------- *!

      !> Bound on `|grad ln rho|` outside the cavity, in 1/Bohr
      real(wp) :: log_grad_out = 8.0_wp
      !> Multiple of the largest nuclear charge bounding `|grad ln rho|` inside
      real(wp) :: log_grad_cusp = 2.3_wp
      !> Largest radius the certificate will claim, in Bohr
      real(wp) :: exclusion_cap = 2.0_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
      !> Override any subset of parameter fields
      procedure, public :: new => new_lsf_isodensity_param
      !> Print the parameters
      procedure, public :: print => print_lsf_isodensity_param
   end type moist_cavity_drop_lsf_isodensity_param_type

contains

   !> Restore compiled parameter defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      class(moist_cavity_drop_lsf_isodensity_param_type), intent(inout) :: self
      type(moist_cavity_drop_lsf_isodensity_param_type) :: defaults

      self%rho_iso = defaults%rho_iso
      self%scale = defaults%scale
      self%log_grad_out = defaults%log_grad_out
      self%log_grad_cusp = defaults%log_grad_cusp
      self%exclusion_cap = defaults%exclusion_cap
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_cavity_drop_lsf_isodensity_param_type), intent(inout), target :: self

      call self%register_real_scalar("rho_iso", self%rho_iso)
      call self%register_real_scalar("scale", self%scale)
      call self%register_real_scalar("log_grad_out", self%log_grad_out)
      call self%register_real_scalar("log_grad_cusp", self%log_grad_cusp)
      call self%register_real_scalar("exclusion_cap", self%exclusion_cap)
   end subroutine register_parameter_entries

   !> Override any subset of isodensity parameter fields
   !>
   !> @param[inout] self     Isodensity parameter instance
   !> @param[in]    rho_iso  Density isovalue (optional)
   !> @param[in]    scale    Constant level set multiplier (optional)
   subroutine new_lsf_isodensity_param(self, rho_iso, scale, log_grad_out, &
                                       log_grad_cusp, exclusion_cap)
      class(moist_cavity_drop_lsf_isodensity_param_type), intent(inout) :: self
      !> Density isovalue (optional override)
      real(wp), intent(in), optional :: rho_iso
      !> Constant level set multiplier (optional override)
      real(wp), intent(in), optional :: scale
      !> Exterior log-density gradient bound (optional override; <= 0 disables the exclusion certificate)
      real(wp), intent(in), optional :: log_grad_out
      !> Cusp multiple of Z_max for the interior bound (optional override)
      real(wp), intent(in), optional :: log_grad_cusp
      !> Cap on the certified radius in Bohr (optional override)
      real(wp), intent(in), optional :: exclusion_cap

      if (present(rho_iso)) self%rho_iso = rho_iso
      if (present(scale)) self%scale = scale
      if (present(log_grad_out)) self%log_grad_out = log_grad_out
      if (present(log_grad_cusp)) self%log_grad_cusp = log_grad_cusp
      if (present(exclusion_cap)) self%exclusion_cap = exclusion_cap
   end subroutine new_lsf_isodensity_param

   !> Physically motivated surface-free radius for an isodensity level set
   !>
   !> A physically motivated bound, not a theorem
   !>
   !> @param[in] self  Isodensity parameters
   !> @param[in] zmax  Largest nuclear charge in the structure
   !> @param[in] lsf0  Level set value at the evaluation point
   !> @returns   r     Surface-free radius (zero when uncertified)
   pure function isodensity_exclusion_radius(self, zmax, lsf0) result(r)
      !> Isodensity parameters
      class(moist_cavity_drop_lsf_isodensity_param_type), intent(in) :: self
      !> Largest nuclear charge in the structure
      real(wp), intent(in) :: zmax
      !> Level set value at the evaluation point
      real(wp), intent(in) :: lsf0
      !> Surface-free radius
      real(wp) :: r

      !> Density in units of the isovalue, and the log-gradient bound in force
      real(wp) :: ratio, kappa

      r = 0.0_wp
      if (self%log_grad_out <= 0.0_wp) return
      if (self%exclusion_cap <= 0.0_wp) return
      if (self%rho_iso <= 0.0_wp .or. self%scale <= 0.0_wp) return

      ! `rho / rho_iso` straight from the value the caller already holds
      ratio = 1.0_wp - lsf0/(self%scale*self%rho_iso)
      if (.not. (ratio > 0.0_wp)) return

      if (lsf0 > 0.0_wp) then
         ! Exterior: the ball stays exterior, so the cusp never enters
         kappa = self%log_grad_out
      else
         ! Interior: the ball can reach a nucleus
         kappa = max(self%log_grad_out, self%log_grad_cusp*zmax)
      end if

      r = min(abs(log(ratio))/kappa, self%exclusion_cap)
   end function isodensity_exclusion_radius

   !> Print the isodensity parameters in the verbose cavity diagnostics
   !>
   !> @param[in] self  Isodensity parameter instance
   !> @param[in] unit  Output unit (default `output_unit`); callers holding a run
   !>                  context pass `ctx%unit` so this honours a log file
   subroutine print_lsf_isodensity_param(self, unit)
      class(moist_cavity_drop_lsf_isodensity_param_type), intent(in) :: self
      !> Output unit override
      integer, intent(in), optional :: unit
      type(prettyprinter) :: pp
      !> Effective output unit
      integer :: iu

      iu = output_unit
      if (present(unit)) iu = unit

      pp = new_prettyprinter(unit=iu)
      call pp%push("Implicit surface (isodensity):")
      call pp%kv("Isovalue (rho_iso)", self%rho_iso, use_exp=.true.)
      call pp%kv("Scale", self%scale)
      if (self%log_grad_out > 0.0_wp) then
         call pp%kv("Exclusion |grad ln rho| (out)", self%log_grad_out, "1/Bohr")
         call pp%kv("Exclusion cusp factor (in)", self%log_grad_cusp)
         call pp%kv("Exclusion radius cap", self%exclusion_cap, "Bohr")
      end if
      call pp%pop()
   end subroutine print_lsf_isodensity_param

end module moist_cavity_drop_lsf_isodensity_param
