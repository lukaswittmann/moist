!> Parameter container shared by the two isodensity LSFs
!>
!> Both isodensity variants describe the same level set
!>
!>    S(r) = scale * (rho_iso - rho(r))
!>
!> and differ only in where the density comes from: the internal variant
!> evaluates it from its own GTO basis, the callback variant asks the host for
!> it. The isovalue and the multiplier are therefore moist-side parameters in
!> both cases, and one type serves both.
module moist_cavity_drop_lsf_isodensity_param
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none
   private

   public :: moist_cavity_drop_lsf_isodensity_param_type
   public :: isodensity_exclusion_radius

   !> Isodensity level set function parameters
   type :: moist_cavity_drop_lsf_isodensity_param_type
      !> Density isovalue defining the surface, in Bohr^-3.
      real(wp) :: rho_iso = 1.0E-3_wp
      !> Constant multiplier applied to the level set value and derivatives.
      !>
      !> It rescales `S` without moving its zero, so the surface is unchanged;
      !> what it buys is a level set whose gradient is O(1) near the surface,
      !> which is what the projection's step control expects.
      real(wp) :: scale = 1.0_wp / 1.0E-3_wp

      !* ------------------- Surface-free exclusion certificate ---------------- *!

      !> Bound on `|grad ln rho|` outside the cavity, in 1/Bohr
      real(wp) :: log_grad_out = 8.0_wp
      !> Multiple of the largest nuclear charge bounding `|grad ln rho|` inside
      real(wp) :: log_grad_cusp = 2.3_wp
      !> Largest radius the certificate will claim, in Bohr
      real(wp) :: exclusion_cap = 2.0_wp
   contains
      !> Override any subset of parameter fields.
      procedure, public :: new => new_lsf_isodensity_param
      !> Print the parameters under an "Implicit surface (isodensity)" header.
      procedure, public :: print => print_lsf_isodensity_param
   end type moist_cavity_drop_lsf_isodensity_param_type

contains

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
   !> This is a *physically motivated* bound, not a theorem
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

      ! `rho / rho_iso` straight from the value the caller already holds. A
      ! non-positive ratio is a density at or below zero -- unreachable for a
      ! physical density, and nothing to certify from.
      ratio = 1.0_wp - lsf0/(self%scale*self%rho_iso)
      if (.not. (ratio > 0.0_wp)) return

      if (lsf0 > 0.0_wp) then
         ! Exterior: the ball stays exterior, so the cusp never enters.
         kappa = self%log_grad_out
      else
         ! Interior: the ball can reach a nucleus.
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
