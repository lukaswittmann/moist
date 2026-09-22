!> Parameter container for the Smooth van der Waals (SvdW) LSF
module moist_cavity_drop_lsf_svdw_param
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)
   private

   public :: moist_cavity_drop_lsf_svdw_param_type

   !> SvdW level set function parameters
   type, extends(moist_model_parameters_type) :: moist_cavity_drop_lsf_svdw_param_type
      !> Blending sharpness k in exp(-k * d)
      real(wp) :: blend_k = 5.5_wp
      !> One-body blending weight
      real(wp) :: blend_1b = 1.0_wp
      !> Two-body blending weight
      real(wp) :: blend_2b = 0.0_wp
      !> Three-body blending weight
      real(wp) :: blend_3b = 3.0_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
      !> Override any subset of parameter fields
      procedure, public :: new => new_lsf_svdw_param
      !> Print the SvdW shape parameters
      procedure, public :: print => print_lsf_svdw_param
   end type moist_cavity_drop_lsf_svdw_param_type

contains

   !> Restore compiled parameter defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      class(moist_cavity_drop_lsf_svdw_param_type), intent(inout) :: self
      type(moist_cavity_drop_lsf_svdw_param_type) :: defaults

      self%blend_k = defaults%blend_k
      self%blend_1b = defaults%blend_1b
      self%blend_2b = defaults%blend_2b
      self%blend_3b = defaults%blend_3b
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_cavity_drop_lsf_svdw_param_type), intent(inout), target :: self

      call self%register_real_scalar("blend_k", self%blend_k)
      call self%register_real_scalar("blend_1b", self%blend_1b)
      call self%register_real_scalar("blend_2b", self%blend_2b)
      call self%register_real_scalar("blend_3b", self%blend_3b)
   end subroutine register_parameter_entries

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
   !> @param[in] unit  Output unit (default `output_unit`); callers holding a run
   !>                  context pass `ctx%unit` so this honours a log file
   subroutine print_lsf_svdw_param(self, unit)
      class(moist_cavity_drop_lsf_svdw_param_type), intent(in) :: self
      !> Output unit override
      integer, intent(in), optional :: unit
      type(prettyprinter) :: pp
      !> Effective output unit
      integer :: iu

      iu = output_unit
      if (present(unit)) iu = unit

      pp = new_prettyprinter(unit=iu)
      call pp%push("Implicit surface (SvdW):")
      call pp%kv("Smoothing (k)", self%blend_k)
      call pp%kv("Smoothing (1b)", self%blend_1b)
      call pp%kv("Smoothing (2b)", self%blend_2b)
      call pp%kv("Smoothing (3b)", self%blend_3b)
      call pp%pop()
   end subroutine print_lsf_svdw_param

end module moist_cavity_drop_lsf_svdw_param
