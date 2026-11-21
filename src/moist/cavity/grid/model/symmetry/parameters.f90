!> Settings of the point-group symmetry model ([[moist_cavity_grid_model_symmetry_type]])
!>
!> A top-level `symmetry` table, beside (not inside) the settings of a cavity:
!>
!> | key                  | values                                   | default |
!> |----------------------|------------------------------------------|---------|
!> | `orientation`        | "auto", "lab", "molecular"               | "auto"  |
!> | `group`              | "off", "auto", a Schoenflies symbol      | "off"   |
!> | `tolerance`          | largest atom deviation under an op, bohr | 1e-8    |
!> | `gap_threshold`      | gap of two eigenvalues of the frame      | 1e-4    |
!> |                      | tensor over their sum                    |         |
!> | `harmonic_threshold` | relative size of an in-plane harmonic    | 1e-4    |
!> | `linear_threshold`   | mean square distance of the atoms        | 1e-4    |
!> |                      | from the long axis, bohr^2 (line)        |         |
!>
!> - thresholds: where the molecular frame switches between two choices
!>   ([[moist_math_symmetry_frame]]); a measure near its threshold is flagged
!> - "auto" orientation: "molecular" with a group (implied by it), otherwise
!>   `default_orientation`
!> - an explicit "lab" with any group but "off" is an error
!> - group symbols are checked against the structure at update time
module moist_cavity_grid_model_symmetry_parameters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io_utils, only: to_lower
   use moist_model_parameters, only: moist_model_parameters_type
   implicit none(type, external)
   private

   public :: moist_cavity_grid_model_symmetry_parameters_type, check_symmetry_parameters

   !> Length of the string settings
   integer, parameter :: setting_len = 16
   !> Orientation without a group under "auto"
   character(len=*), parameter :: default_orientation = "lab"

   !> Settings of the symmetry model
   type, extends(moist_model_parameters_type) :: moist_cavity_grid_model_symmetry_parameters_type
      !> Sphere orientation: "auto", "lab", or "molecular", case-insensitive
      !>
      !> - "auto" resolves through `effective_orientation`
      character(len=setting_len) :: orientation = "auto"
      !> Point group: "off", "auto" or a Schoenflies symbol, case-insensitive
      character(len=setting_len) :: group = "off"
      !> Largest atom deviation under a symmetry operation, bohr
      real(wp) :: tolerance = 1.0e-8_wp
      !> Relative gap of two frame-tensor eigenvalues separating the pair (relative to their sum)
      real(wp) :: gap_threshold = 1.0e-4_wp
      !> Relative size of an in-plane harmonic above which it fixes the angle
      real(wp) :: harmonic_threshold = 1.0e-4_wp
      !> Mean square distance from the long axis making a line, bohr^2 (atoms below it form a line)
      real(wp) :: linear_threshold = 1.0e-4_wp
   contains
      !> Restore compiled defaults
      procedure :: init_defaults => init_symmetry_defaults
      !> Declare the `symmetry` table for file input, output and printing
      procedure :: register_entries => register_symmetry_entries
      !> Refuse invalid values after file input
      procedure :: validate => validate_symmetry_parameters
      !> Orientation the settings ask for, "lab" or "molecular"
      procedure :: effective_orientation => symmetry_effective_orientation
      !> Whether "auto" orientation is resolved to molecular by a group
      procedure :: orientation_implied => symmetry_orientation_implied
      !> Whether a group is requested, anything but "off"
      procedure :: has_group => symmetry_has_group
   end type moist_cavity_grid_model_symmetry_parameters_type

contains

   !> Restore compiled defaults: orientation auto, group off, 1e-8 bohr
   !>
   !> @param[in,out] self Symmetry settings
   subroutine init_symmetry_defaults(self)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(inout) :: self

      !> Compiled defaults
      type(moist_cavity_grid_model_symmetry_parameters_type) :: defaults

      self%orientation = defaults%orientation
      self%group = defaults%group
      self%tolerance = defaults%tolerance
      self%gap_threshold = defaults%gap_threshold
      self%harmonic_threshold = defaults%harmonic_threshold
      self%linear_threshold = defaults%linear_threshold
   end subroutine init_symmetry_defaults

   !> Declare the `symmetry` table
   !>
   !> - `symmetry.orientation` and `symmetry.group` print at the normal level,
   !>   `symmetry.tolerance` and the thresholds at the detailed one
   !>
   !> @param[in,out] self Symmetry settings
   subroutine register_symmetry_entries(self)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(inout), target :: self

      call self%register_string("symmetry.orientation", self%orientation, level=2)
      call self%register_string("symmetry.group", self%group, level=2)
      call self%register_real_scalar("symmetry.tolerance", self%tolerance)
      call self%register_real_scalar("symmetry.gap_threshold", self%gap_threshold)
      call self%register_real_scalar("symmetry.harmonic_threshold", self%harmonic_threshold)
      call self%register_real_scalar("symmetry.linear_threshold", self%linear_threshold)
   end subroutine register_symmetry_entries

   !> Refuse invalid values after file input
   !>
   !> @param[in,out] self  Symmetry settings
   !> @param[out]    error Invalid setting
   subroutine validate_symmetry_parameters(self, error)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(inout) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call check_symmetry_parameters(self, error)
   end subroutine validate_symmetry_parameters

   !> Refuse an unknown orientation, an empty group, a bad tolerance or threshold
   !>
   !> - an explicit "lab" with any group but "off" is refused: a group needs
   !>   the molecular orientation, which "auto" implies
   !> - group symbols are checked against the structure at update time
   !>
   !> @param[in]  param Symmetry settings
   !> @param[out] error Invalid setting
   subroutine check_symmetry_parameters(param, error)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: param
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select case (to_lower(trim(adjustl(param%orientation))))
      case ("auto", "lab", "molecular")
      case default
         call fatal_error(error, "Symmetry: orientation must be 'auto', 'lab' or 'molecular', got '"// &
            & trim(param%orientation)//"'")
         return
      end select
      if (len_trim(param%group) == 0) then
         call fatal_error(error, "Symmetry: group must be 'off', 'auto' or a Schoenflies symbol")
         return
      end if
      if (to_lower(trim(adjustl(param%orientation))) == "lab" .and. param%has_group()) then
         call fatal_error(error, 'symmetry.group = "'//trim(adjustl(param%group))// &
            & '" requires symmetry.orientation = "molecular"')
         return
      end if
      if (.not. (ieee_is_finite(param%tolerance) .and. param%tolerance > 0.0_wp)) then
         call fatal_error(error, "Symmetry: tolerance must be finite and positive")
         return
      end if
      if (.not. all(ieee_is_finite([param%gap_threshold, param%harmonic_threshold, &
         & param%linear_threshold]) .and. [param%gap_threshold, param%harmonic_threshold, &
         & param%linear_threshold] > 0.0_wp)) then
         call fatal_error(error, "Symmetry: a threshold must be finite and positive")
         return
      end if
      if (.not. (param%gap_threshold < 1.0_wp .and. param%harmonic_threshold < 1.0_wp)) then
         call fatal_error(error, "Symmetry: gap_threshold and harmonic_threshold must be below one")
      end if
   end subroutine check_symmetry_parameters

   !> Orientation the settings ask for
   !>
   !> - "lab" or "molecular"; "auto" resolves to "molecular" with a group and
   !>   `default_orientation` without
   !>
   !> @param[in] self Symmetry settings
   pure function symmetry_effective_orientation(self) result(orientation)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: self
      !> "lab" or "molecular", or an unchecked value
      character(len=:), allocatable :: orientation

      orientation = to_lower(trim(adjustl(self%orientation)))
      if (orientation /= "auto") return
      if (self%has_group()) then
         orientation = "molecular"
      else
         orientation = default_orientation
      end if
   end function symmetry_effective_orientation

   !> Whether "auto" orientation is resolved to molecular by a group
   !>
   !> @param[in] self Symmetry settings
   pure function symmetry_orientation_implied(self) result(implied)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: self
      !> Whether the group implies the molecular orientation
      logical :: implied

      implied = to_lower(trim(adjustl(self%orientation))) == "auto" .and. self%has_group()
   end function symmetry_orientation_implied

   !> Whether a group is requested, anything but "off"
   !>
   !> @param[in] self Symmetry settings
   pure function symmetry_has_group(self) result(has)
      !> Symmetry settings
      class(moist_cavity_grid_model_symmetry_parameters_type), intent(in) :: self
      !> Whether the group is not "off"
      logical :: has

      has = to_lower(trim(adjustl(self%group))) /= "off"
   end function symmetry_has_group

end module moist_cavity_grid_model_symmetry_parameters
