!> Parameter container for the COSMO Fine Cavity (CFC) LSF
module moist_cavity_drop_lsf_cfc_param
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none
   private

   public :: moist_cavity_drop_lsf_cfc_param_type

   !> CFC level set function parameters (Diedenhofen-Klamt 2018 defaults)
   type :: moist_cavity_drop_lsf_cfc_param_type
      !> Atomic-term exponent a1.
      real(wp) :: a1 = -15.0_wp
      !> Pair-term exponent a2.
      real(wp) :: a2 = -9.0_wp
      !> Pair-term coupling constant c.
      real(wp) :: c = 5.0_wp
      !> Pair-term polynomial power m.
      !>
      !> The kernel uses `m = 4` for the implemented symbolic differentiation
      !> A value other than 4 is inconsistent with the derivatives
      integer :: m = 4
   contains
      !> Override any subset of parameter fields.
      procedure, public :: new => new_lsf_cfc_param
      !> Print the CFC shape parameters under an "Implicit surface (CFC)" header.
      procedure, public :: print => print_lsf_cfc_param
   end type moist_cavity_drop_lsf_cfc_param_type

contains

   !> Override any subset of CFC parameter fields
   !>
   !> @param[inout] self  CFC parameter instance
   !> @param[in]    a1    Atomic-term exponent (optional)
   !> @param[in]    a2    Pair-term exponent (optional)
   !> @param[in]    c     Pair-term coupling (optional)
   !> @param[in]    m     Pair-term power (optional; screening only)
   subroutine new_lsf_cfc_param(self, a1, a2, c, m)
      class(moist_cavity_drop_lsf_cfc_param_type), intent(inout) :: self
      !> Atomic-term exponent (optional override)
      real(wp), intent(in), optional :: a1
      !> Pair-term exponent (optional override)
      real(wp), intent(in), optional :: a2
      !> Pair-term coupling (optional override)
      real(wp), intent(in), optional :: c
      !> Pair-term power (optional override)
      integer, intent(in), optional :: m

      if (present(a1)) self%a1 = a1
      if (present(a2)) self%a2 = a2
      if (present(c)) self%c = c
      if (present(m)) self%m = m
   end subroutine new_lsf_cfc_param

   !> Print the CFC shape parameters in the verbose cavity diagnostics
   !>
   !> @param[in] self  CFC parameter instance
   !> @param[in] unit  Output unit (default `output_unit`); callers holding a run
   !>                  context pass `ctx%unit` so this honours a log file
   subroutine print_lsf_cfc_param(self, unit)
      class(moist_cavity_drop_lsf_cfc_param_type), intent(in) :: self
      !> Output unit override
      integer, intent(in), optional :: unit
      type(prettyprinter) :: pp
      !> Effective output unit
      integer :: iu

      iu = output_unit
      if (present(unit)) iu = unit

      pp = new_prettyprinter(unit=iu)
      call pp%push("Implicit surface (CFC):")
      call pp%kv("Atomic exponent (a1)", self%a1)
      call pp%kv("Pair exponent (a2)", self%a2)
      call pp%kv("Pair coupling (c)", self%c)
      call pp%kv("Pair power (m)", self%m)
      call pp%pop()
   end subroutine print_lsf_cfc_param

end module moist_cavity_drop_lsf_cfc_param
