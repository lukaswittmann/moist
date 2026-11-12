!> Parameter container for the COSMO Fine Cavity (CFC) LSF
module moist_cavity_drop_lsf_cfc_param
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use moist_model_parameters, only: moist_model_parameters_type
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use, intrinsic :: iso_fortran_env, only: output_unit
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)
   private

   public :: moist_cavity_drop_lsf_cfc_param_type

   !> Atom radius of the printed band in Bohr (past the surface)
   real(wp), parameter :: print_ref_radius = 2.0_wp

   !> CFC level set function parameters (Diedenhofen-Klamt 2018 defaults)
   type, extends(moist_model_parameters_type) :: moist_cavity_drop_lsf_cfc_param_type
      !> Atomic-term exponent a1
      real(wp) :: a1 = -15.0_wp
      !> Pair-term exponent a2
      real(wp) :: a2 = -9.0_wp
      !> Pair-term coupling constant c
      real(wp) :: c = 5.0_wp
      !> Pair-term polynomial power m (only 4 is supported)
      integer :: m = 4
      !> Fade every atom and pair term out over a finite band (default on)
      !>
      !> - band in e-folds of `u_A = exp(a2 sigma_A)`, `sigma_A = x/R_A - 1`:
      !>   onset `fade_start/|a2|`, width `fade_width/|a2|` in sigma, i.e. per atom
      !>   `R_A fade_start/|a2|` and `R_A fade_width/|a2|` Bohr past its surface
      !> - on: the band end `R_A (a + b)/|a2|` is the atom reach,
      !>   `screening_threshold` is not read by the CFC model
      !> - both paths: `evaluate_batch` and `prepare` / `prepare_subset` with every
      !>   per-atom family evaluate the same faded model
      !> - off: the untruncated model, screened by `screening_threshold`
      logical :: fade_enabled = .true.
      !> Fade onset `a` in e-folds of `u_A`
      !>
      !> - 27 with `fade_width` 3: reach 30 e-folds
      !>   (4.33 R_A at `a2 = -9`), surface moved by <= 4.2e-9 Bohr; the pair term
      !>   at band start is ~`e^(-a_eff)`, `a_eff = a - ln c - m ln(2 + a/|a2|)`,
      !>   5.9e-9 per pair
      real(wp) :: fade_start = 27.0_wp
      !> Fade band width `b` in e-folds of `u_A`
      real(wp) :: fade_width = 3.0_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_parameter_defaults
      !> Declare fields for JSON input, output, and printing
      procedure :: register_entries => register_parameter_entries
      !> Refuse a pair power other than 4 or an invalid fade band
      procedure :: validate => validate_lsf_cfc_param
      !> Override any subset of parameter fields
      procedure, public :: new => new_lsf_cfc_param
      !> Print the CFC shape parameters under an "Implicit surface (CFC)" header
      procedure, public :: print => print_lsf_cfc_param
   end type moist_cavity_drop_lsf_cfc_param_type

contains

   !> Restore compiled parameter defaults
   !>
   !> @param[inout] self Parameter values
   subroutine init_parameter_defaults(self)
      class(moist_cavity_drop_lsf_cfc_param_type), intent(inout) :: self
      type(moist_cavity_drop_lsf_cfc_param_type) :: defaults

      self%a1 = defaults%a1
      self%a2 = defaults%a2
      self%c = defaults%c
      self%m = defaults%m
      self%fade_enabled = defaults%fade_enabled
      self%fade_start = defaults%fade_start
      self%fade_width = defaults%fade_width
   end subroutine init_parameter_defaults

   !> Declare parameter fields for JSON input, output, and printing
   !>
   !> @param[inout] self Parameter values
   subroutine register_parameter_entries(self)
      class(moist_cavity_drop_lsf_cfc_param_type), intent(inout), target :: self

      call self%register_real_scalar("a1", self%a1)
      call self%register_real_scalar("a2", self%a2)
      call self%register_real_scalar("c", self%c)
      call self%register_int_scalar("m", self%m)
      call self%register_logical("fade_enabled", self%fade_enabled)
      call self%register_real_scalar("fade_start", self%fade_start)
      call self%register_real_scalar("fade_width", self%fade_width)
   end subroutine register_parameter_entries

   !> Refuse a pair power other than 4 or an invalid fade band
   !>
   !> - always: `m = 4`, the generated kernel bakes it in
   !> - fade off: the band is not read, so not checked
   !> - fade on: finite, positive `fade_start` and `fade_width`, and `a2 < 0`
   !>   (the band is counted in e-folds of `exp(a2 sigma)`)
   !>
   !> @param[inout] self  CFC parameter instance
   !> @param[out]   error Invalid pair power or fade band
   subroutine validate_lsf_cfc_param(self, error)
      !> CFC parameter instance
      class(moist_cavity_drop_lsf_cfc_param_type), intent(inout) :: self
      !> Invalid pair power or fade band
      type(error_type), allocatable, intent(out) :: error

      if (self%m /= 4) then
         call fatal_error(error, "CFC pair power m must be 4 (the generated kernel bakes it in)")
         return
      end if
      if (.not. self%fade_enabled) return
      if (.not. (ieee_is_finite(self%fade_start) .and. ieee_is_finite(self%fade_width))) then
         call fatal_error(error, "CFC fade_start and fade_width must be finite")
      else if (.not. (self%fade_start > 0.0_wp .and. self%fade_width > 0.0_wp)) then
         call fatal_error(error, "CFC fade_start and fade_width must be positive")
      else if (.not. self%a2 < 0.0_wp) then
         call fatal_error(error, "CFC fade needs a2 < 0")
      end if
   end subroutine validate_lsf_cfc_param

   !> Override any subset of CFC parameter fields
   !>
   !> @param[inout] self  CFC parameter instance
   !> @param[in]    a1    Atomic-term exponent (optional)
   !> @param[in]    a2    Pair-term exponent (optional)
   !> @param[in]    c     Pair-term coupling (optional)
   !> @param[in]    m     Pair-term power (optional; must be 4)
   !> @param[in]    fade_enabled Fade the atom and pair terms over a finite band (optional)
   !> @param[in]    fade_start   Fade onset in e-folds of u (optional)
   !> @param[in]    fade_width   Fade band width in e-folds of u (optional)
   subroutine new_lsf_cfc_param(self, a1, a2, c, m, fade_enabled, fade_start, fade_width)
      class(moist_cavity_drop_lsf_cfc_param_type), intent(inout) :: self
      !> Atomic-term exponent (optional override)
      real(wp), intent(in), optional :: a1
      !> Pair-term exponent (optional override)
      real(wp), intent(in), optional :: a2
      !> Pair-term coupling (optional override)
      real(wp), intent(in), optional :: c
      !> Pair-term power (optional override)
      integer, intent(in), optional :: m
      !> Fade switch (optional override)
      logical, intent(in), optional :: fade_enabled
      !> Fade onset in e-folds (optional override)
      real(wp), intent(in), optional :: fade_start
      !> Fade band width in e-folds (optional override)
      real(wp), intent(in), optional :: fade_width

      if (present(a1)) self%a1 = a1
      if (present(a2)) self%a2 = a2
      if (present(c)) self%c = c
      if (present(m)) self%m = m
      if (present(fade_enabled)) self%fade_enabled = fade_enabled
      if (present(fade_start)) self%fade_start = fade_start
      if (present(fade_width)) self%fade_width = fade_width
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
      call pp%kv("Fade", self%fade_enabled)
      if (self%fade_enabled) then
         call pp%kv("Fade start", self%fade_start, "e-folds")
         call pp%kv("Fade width", self%fade_width, "e-folds")
         if (self%a2 < 0.0_wp) then
            call pp%kv("Fade onset (sigma)", self%fade_start/abs(self%a2), "R_A")
            call pp%kv("Fade end, reach (sigma)", &
                       (self%fade_start + self%fade_width)/abs(self%a2), "R_A")
            call pp%kv("Fade onset (R_A = 2 Bohr)", &
                       print_ref_radius*self%fade_start/abs(self%a2), "Bohr")
            call pp%kv("Fade end, reach (R_A = 2 Bohr)", &
                       print_ref_radius*(self%fade_start + self%fade_width)/abs(self%a2), "Bohr")
         end if
      end if
      call pp%pop()
      call pp%blank()
   end subroutine print_lsf_cfc_param

end module moist_cavity_drop_lsf_cfc_param
