!> Parameter container for the Smooth van der Waals (SvdW) LSF
module moist_cavity_drop_lsf_svdw_param
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none
   private

   public :: moist_cavity_drop_lsf_svdw_param_type

   !> SvdW level-set function parameters
   type :: moist_cavity_drop_lsf_svdw_param_type
      !> Blending sharpness k in exp(-k * d).
      real(wp) :: blend_k = 3.0_wp
      !> One-body blending weight.
      real(wp) :: blend_1b = 1.0_wp
      !> Two-body blending weight.
      real(wp) :: blend_2b = 1.0_wp
      !> Three-body blending weight.
      real(wp) :: blend_3b = 1.0_wp
   contains
      !> Override any subset of parameter fields.
      procedure, public :: new => new_lsf_svdw_param
      !> Print the SvdW shape parameters under an "Implicit surface (SvdW)" header.
      procedure, public :: print => print_lsf_svdw_param
   end type moist_cavity_drop_lsf_svdw_param_type

contains

   !> Override any subset of SvdW parameter fields
   !>
   !> @param[inout] self     SvdW parameter instance
   !> @param[in]    blend_k  Blending sharpness k (optional)
   !> @param[in]    blend_1b One-body weight (optional)
   !> @param[in]    blend_2b Two-body weight (optional)
   !> @param[in]    blend_3b Three-body weight (optional)
   subroutine new_lsf_svdw_param(self, blend_k, blend_1b, blend_2b, blend_3b)
      class(moist_cavity_drop_lsf_svdw_param_type), intent(inout) :: self
      !> Blending sharpness k (optional override)
      real(wp), intent(in), optional :: blend_k
      !> One-body weight (optional override)
      real(wp), intent(in), optional :: blend_1b
      !> Two-body weight (optional override)
      real(wp), intent(in), optional :: blend_2b
      !> Three-body weight (optional override)
      real(wp), intent(in), optional :: blend_3b

      if (present(blend_k)) self%blend_k = blend_k
      if (present(blend_1b)) self%blend_1b = blend_1b
      if (present(blend_2b)) self%blend_2b = blend_2b
      if (present(blend_3b)) self%blend_3b = blend_3b
   end subroutine new_lsf_svdw_param

   !> Print the SvdW shape parameters in the verbose cavity diagnostics
   !>
   !> @param[in] self  SvdW parameter instance
   subroutine print_lsf_svdw_param(self)
      class(moist_cavity_drop_lsf_svdw_param_type), intent(in) :: self
      type(prettyprinter) :: pp

      pp = new_prettyprinter(unit=output_unit)
      call pp%push('Implicit surface (SvdW):')
      call pp%kv('Smoothing (k)', self%blend_k)
      call pp%kv('Smoothing (1b)', self%blend_1b)
      call pp%kv('Smoothing (2b)', self%blend_2b)
      call pp%kv('Smoothing (3b)', self%blend_3b)
      call pp%pop()
   end subroutine print_lsf_svdw_param

end module moist_cavity_drop_lsf_svdw_param
