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
   subroutine new_lsf_isodensity_param(self, rho_iso, scale)
      class(moist_cavity_drop_lsf_isodensity_param_type), intent(inout) :: self
      !> Density isovalue (optional override)
      real(wp), intent(in), optional :: rho_iso
      !> Constant level set multiplier (optional override)
      real(wp), intent(in), optional :: scale

      if (present(rho_iso)) self%rho_iso = rho_iso
      if (present(scale)) self%scale = scale
   end subroutine new_lsf_isodensity_param

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
      call pp%pop()
   end subroutine print_lsf_isodensity_param

end module moist_cavity_drop_lsf_isodensity_param
