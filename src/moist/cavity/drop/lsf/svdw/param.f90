!> Parameter container for the Smooth van der Waals (SvdW) LSF
module moist_cavity_drop_lsf_svdw_param
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env, only: error_type, fatal_error
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
      !> Fade every atom term out over a finite band (off: untruncated model)
      !>
      !> - on (default): the band end `3(a + b)/blend_k` is the atom reach,
      !>   `screening_threshold` is not read by the SvdW model
      !> - off: hard screening at `-3 ln(threshold)/blend_k`
      logical :: fade_enabled = .true.
      !> Fade onset `a` in e-folds of u, `X_on = 3a/blend_k` in `d_A = x - R_A`
      !>
      !> - default 21 with `fade_width` 3: reach 24 e-folds, onset 11.45, end
      !>   13.09 Bohr at k = 5.5; measured surface displacement against the
      !>   untruncated model <= 4.4e-9 Bohr on the 1000-atom lattice
      real(wp) :: fade_start = 21.0_wp
      !> Fade band width `b` in e-folds of u, `Delta = 3b/blend_k`
      real(wp) :: fade_width = 3.0_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
      !> Refuse a non-positive or non-finite `blend_k`, or such a fade band while the fade is on
      procedure :: validate => validate_lsf_svdw_param
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
      self%fade_enabled = defaults%fade_enabled
      self%fade_start = defaults%fade_start
      self%fade_width = defaults%fade_width
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
      call self%register_logical("fade_enabled", self%fade_enabled)
      call self%register_real_scalar("fade_start", self%fade_start)
      call self%register_real_scalar("fade_width", self%fade_width)
   end subroutine register_parameter_entries

   !> Refuse a non-positive or non-finite `blend_k`, or such a fade band while the fade is on
   !>
   !> - always: finite `blend_k > 0`
   !> - fade off: the band is not read, so not checked
   !> - fade on: the band end is the atom reach, so no band-vs-reach check exists
   !>
   !> @param[inout] self  SvdW parameter instance
   !> @param[out]   error Invalid `blend_k` or fade band
   subroutine validate_lsf_svdw_param(self, error)
      !> SvdW parameter instance
      class(moist_cavity_drop_lsf_svdw_param_type), intent(inout) :: self
      !> Invalid `blend_k` or fade band
      type(error_type), allocatable, intent(out) :: error

      if (.not. (ieee_is_finite(self%blend_k) .and. self%blend_k > 0.0_wp)) then
         call fatal_error(error, "SvdW blend_k must be finite and positive")
         return
      end if
      if (.not. self%fade_enabled) return
      if (.not. (ieee_is_finite(self%fade_start) .and. ieee_is_finite(self%fade_width))) then
         call fatal_error(error, "SvdW fade_start and fade_width must be finite")
      else if (.not. (self%fade_start > 0.0_wp .and. self%fade_width > 0.0_wp)) then
         call fatal_error(error, "SvdW fade_start and fade_width must be positive")
      end if
   end subroutine validate_lsf_svdw_param

   !> Override any subset of SvdW parameter fields
   !>
   !> @param[inout] self     SvdW parameter instance
   !> @param[in]    blend_k  Blending sharpness k (optional)
   !> @param[in]    blend_1b One-body weight (optional)
   !> @param[in]    blend_2b Two-body weight (optional)
   !> @param[in]    blend_3b Three-body weight (optional)
   !> @param[in]    fade_enabled Fade the atom terms over a finite band (optional)
   !> @param[in]    fade_start   Fade onset in e-folds of u (optional)
   !> @param[in]    fade_width   Fade band width in e-folds of u (optional)
   subroutine new_lsf_svdw_param(self, blend_k, blend_1b, blend_2b, blend_3b, fade_enabled, &
                                 fade_start, fade_width)
      class(moist_cavity_drop_lsf_svdw_param_type), intent(inout) :: self
      !> Blending sharpness k (optional override)
      real(wp), intent(in), optional :: blend_k
      !> One-body weight (optional override)
      real(wp), intent(in), optional :: blend_1b
      !> Two-body weight (optional override)
      real(wp), intent(in), optional :: blend_2b
      !> Three-body weight (optional override)
      real(wp), intent(in), optional :: blend_3b
      !> Fade switch (optional override)
      logical, intent(in), optional :: fade_enabled
      !> Fade onset in e-folds (optional override)
      real(wp), intent(in), optional :: fade_start
      !> Fade band width in e-folds (optional override)
      real(wp), intent(in), optional :: fade_width

      if (present(blend_k)) self%blend_k = blend_k
      if (present(blend_1b)) self%blend_1b = blend_1b
      if (present(blend_2b)) self%blend_2b = blend_2b
      if (present(blend_3b)) self%blend_3b = blend_3b
      if (present(fade_enabled)) self%fade_enabled = fade_enabled
      if (present(fade_start)) self%fade_start = fade_start
      if (present(fade_width)) self%fade_width = fade_width
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
      call pp%kv("Fade", self%fade_enabled)
      if (self%fade_enabled) then
         call pp%kv2("Fade start", self%fade_start, "e-folds", 3.0_wp*self%fade_start/self%blend_k, "Bohr")
         call pp%kv2("Fade width", self%fade_width, "e-folds", 3.0_wp*(self%fade_start + self%fade_width)/ &
                       self%blend_k, "Bohr")
      end if
      call pp%pop()
      call pp%blank()
   end subroutine print_lsf_svdw_param

end module moist_cavity_drop_lsf_svdw_param
