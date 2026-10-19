!> Independent mesh and quadrature controls for a frozen MT cavity
module moist_cavity_mt_parameters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_model_parameters, only: moist_model_parameters_type
   implicit none(type, external)
   private
   public :: moist_cavity_mt_parameters_type

   !> Marching-tetrahedra construction settings, in atomic units
   type, extends(moist_model_parameters_type) :: moist_cavity_mt_parameters_type
      !> Maximum Cartesian cell side, bohr
      real(wp) :: spacing = 0.2_wp
      !> Requested geometric area per integration point, bohr**2
      real(wp) :: target_area = 0.35_wp
      !> Relative patch-area tolerance, before critical switching
      real(wp) :: area_tolerance = 0.1_wp
      !> Dimensionless coefficient in xi0 = gaussian_scale / sqrt(area)
      real(wp) :: gaussian_scale = 4.9_wp
      !> Gradient magnitude below which the critical switch is zero
      real(wp) :: w_0ls_from = 0.25_wp
      !> Gradient magnitude above which the critical switch is one
      real(wp) :: w_0ls_to = 0.6_wp
      !> Symmetric critical-bump exponent, as in DROP
      real(wp) :: w_0ls_p = 0.8_wp
      !> Symmetric critical-bump prefactor, as in DROP
      real(wp) :: w_0ls_a = 1.6_wp
      !> Dimensionless cutoff on the critical switch, before normalizing gradients
      real(wp) :: weight_cutoff = 1.0e-12_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_mt_defaults
      !> Register JSON/TOML fields
      procedure :: register_entries => register_mt_entries
      !> Reject nonfinite and inadmissible settings
      procedure :: validate => validate_mt_parameters
   end type moist_cavity_mt_parameters_type

contains

   !> Restore compiled defaults
   !>
   !> @param[in,out] self Construction settings
   subroutine init_mt_defaults(self)
      !> Construction settings
      class(moist_cavity_mt_parameters_type), intent(inout) :: self
      type(moist_cavity_mt_parameters_type) :: defaults

      self%spacing = defaults%spacing
      self%target_area = defaults%target_area
      self%area_tolerance = defaults%area_tolerance
      self%gaussian_scale = defaults%gaussian_scale
      self%w_0ls_from = defaults%w_0ls_from
      self%w_0ls_to = defaults%w_0ls_to
      self%w_0ls_p = defaults%w_0ls_p
      self%w_0ls_a = defaults%w_0ls_a
      self%weight_cutoff = defaults%weight_cutoff
   end subroutine init_mt_defaults

   !> Register the native parameter schema
   !>
   !> @param[in,out] self Construction settings
   subroutine register_mt_entries(self)
      !> Construction settings
      class(moist_cavity_mt_parameters_type), intent(inout), target :: self

      call self%register_real_scalar("spacing", self%spacing, level=2)
      call self%register_real_scalar("target_area", self%target_area, level=2)
      call self%register_real_scalar("area_tolerance", self%area_tolerance)
      call self%register_real_scalar("gaussian_scale", self%gaussian_scale, level=2)
      call self%register_real_scalar("w_0ls_from", self%w_0ls_from)
      call self%register_real_scalar("w_0ls_to", self%w_0ls_to)
      call self%register_real_scalar("w_0ls_p", self%w_0ls_p)
      call self%register_real_scalar("w_0ls_a", self%w_0ls_a)
      call self%register_real_scalar("weight_cutoff", self%weight_cutoff)
   end subroutine register_mt_entries

   !> Validate finite geometry, width, and critical-switch controls
   !>
   !> @param[in,out] self Construction settings
   !> @param[out] error Invalid setting
   subroutine validate_mt_parameters(self, error)
      !> Construction settings
      class(moist_cavity_mt_parameters_type), intent(inout) :: self
      !> Invalid setting
      type(error_type), allocatable, intent(out) :: error

      if (.not. all(ieee_is_finite([self%spacing, self%target_area, self%area_tolerance, &
                                    self%gaussian_scale, self%w_0ls_from, self%w_0ls_to, &
                                    self%w_0ls_p, self%w_0ls_a, self%weight_cutoff]))) then
         call fatal_error(error, "Marching tetrahedra parameters must be finite")
      else if (self%spacing <= 0.0_wp .or. self%target_area <= 0.0_wp) then
         call fatal_error(error, "Marching tetrahedra spacing and target_area must be positive")
      else if (self%area_tolerance <= 0.0_wp .or. self%area_tolerance >= 1.0_wp) then
         call fatal_error(error, "Marching tetrahedra area_tolerance must lie strictly between zero and one")
      else if (self%gaussian_scale <= 0.0_wp) then
         call fatal_error(error, "Marching tetrahedra gaussian_scale must be positive")
      else if (self%w_0ls_from < 0.0_wp .or. self%w_0ls_to <= self%w_0ls_from) then
         call fatal_error(error, "Marching tetrahedra requires 0 <= w_0ls_from < w_0ls_to")
      else if (self%w_0ls_p <= 0.0_wp .or. self%w_0ls_a <= 0.0_wp) then
         call fatal_error(error, "Marching tetrahedra critical-bump parameters must be positive")
      else if (self%weight_cutoff < 0.0_wp .or. self%weight_cutoff >= 1.0_wp) then
         call fatal_error(error, "Marching tetrahedra weight_cutoff must lie in [0,1)")
      end if
   end subroutine validate_mt_parameters

end module moist_cavity_mt_parameters
