module moist_radii
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error
   use moist_data_radii_legacy, only: rad_type
   use moist_radius_type, only: radius_type
   use moist_radii_static, only: static_radius_type
   use moist_radii_static, only: new_cpcm_radii, new_smd_radii
   use moist_radii_static, only: new_d3_radii, new_cosmo_radii, new_bondi_radii
   use moist_radii_draco, only: draco_radius_type, new_draco_radii
   use moist_radii_custom, only: custom_radius_type
   use moist_radii_custom, only: new_custom_radii_atoms, new_custom_radii_elements
   use mctc_io_utils, only: to_lower
   implicit none
   private

   public :: radius_type
   public :: static_radius_type
   public :: draco_radius_type
   public :: custom_radius_type
   public :: new_cpcm_radii
   public :: new_smd_radii
   public :: new_d3_radii
   public :: new_cosmo_radii
   public :: new_bondi_radii
   public :: new_draco_radii
   public :: new_radii_custom_atoms
   public :: new_radii_custom_elements
   public :: new_radii
   public :: default_cpcm_radii

   interface new_radii
      module procedure new_radii_int
      module procedure new_radii_str
   end interface new_radii

contains

   !> Return a default CPCM static radii model object.
   !> @param[in]  verbosity  optional print level for diagnostics
   !> @return radii          CPCM static radii model
   function default_cpcm_radii(verbosity) result(radii)
      !> CPCM static radii model
      type(static_radius_type) :: radii
      !> Optional print level
      integer, intent(in), optional :: verbosity

      call new_cpcm_radii(radii, verbosity)
   end function default_cpcm_radii

   !> Construct a custom radii object from per-atom radii.
   !> @param[in]  radii  per-atom radii in bohr (size must match mol%nat in update)
   !> @param[out] model      constructed custom radii model
   !> @param[out] error      allocated on invalid input
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_radii_custom_atoms(radii, model, error, verbosity)
      !> Per-atom radii in bohr
      real(wp), intent(in) :: radii(:)
      !> Constructed radii model
      class(radius_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional print level
      integer, intent(in), optional :: verbosity

      allocate (custom_radius_type :: model)
      select type (model)
      type is (custom_radius_type)
         call new_custom_radii_atoms(radii, model, error, verbosity)
      end select
   end subroutine new_radii_custom_atoms

   !> Construct a custom radii object from per-element radii.
   !> @param[in]  atomic_numbers  atomic numbers with custom radii
   !> @param[in]  radii           radii matching atomic_numbers in bohr
   !> @param[out] model           constructed custom radii model
   !> @param[out] error           allocated on invalid input
   !> @param[in]  verbosity       optional print level for diagnostics
   subroutine new_radii_custom_elements(atomic_numbers, radii, model, error, verbosity)
      !> Atomic numbers with custom radii
      integer, intent(in) :: atomic_numbers(:)
      !> Radii matching atomic_numbers in bohr
      real(wp), intent(in) :: radii(:)
      !> Constructed radii model
      class(radius_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional print level
      integer, intent(in), optional :: verbosity

      allocate (custom_radius_type :: model)
      select type (model)
      type is (custom_radius_type)
         call new_custom_radii_elements(atomic_numbers, radii, model, error, verbosity)
      end select
   end subroutine new_radii_custom_elements

   !> Construct a radii object from an integer model tag.
   !> @param[in]  model_tag  radius type tag
   !> @param[out] model      constructed radii model
   !> @param[out] error      allocated on unknown model tag
   !> @param[in]  verbosity  optional print level for diagnostics
   subroutine new_radii_int(model_tag, model, error, verbosity)
      !> Radius type tag
      integer, intent(in) :: model_tag
      !> Constructed radii model
      class(radius_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional print level
      integer, intent(in), optional :: verbosity

      character(len=128) :: msg

      select case (model_tag)
      case (rad_type%cpcm)
         allocate (static_radius_type :: model)
         select type (model)
         type is (static_radius_type)
            call new_cpcm_radii(model, verbosity)
         end select
      case (rad_type%smd)
         allocate (static_radius_type :: model)
         select type (model)
         type is (static_radius_type)
            call new_smd_radii(model, verbosity)
         end select
      case (rad_type%d3)
         allocate (static_radius_type :: model)
         select type (model)
         type is (static_radius_type)
            call new_d3_radii(model, verbosity)
         end select
      case (rad_type%cosmo)
         allocate (static_radius_type :: model)
         select type (model)
         type is (static_radius_type)
            call new_cosmo_radii(model, verbosity)
         end select
      case (rad_type%bondi)
         allocate (static_radius_type :: model)
         select type (model)
         type is (static_radius_type)
            call new_bondi_radii(model, verbosity)
         end select
      case default
         write (msg, '(a,i0)') "Unknown radius type tag: ", model_tag
         call fatal_error(error, trim(msg))
         return
      end select
   end subroutine new_radii_int

   !> Construct a radii object from a model name.
   !> @param[in]  model_name  radius model name
   !> @param[out] model       constructed radii model
   !> @param[out] error       allocated on unknown model name
   !> @param[in]  verbosity   optional print level for diagnostics
   subroutine new_radii_str(model_name, model, error, verbosity)
      !> Radius model name
      character(len=*), intent(in) :: model_name
      !> Constructed radii model
      class(radius_type), allocatable, intent(out) :: model
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional print level
      integer, intent(in), optional :: verbosity

      integer :: model_tag
      character(len=32) :: mstr
      character(len=128) :: msg

      mstr = to_lower(trim(adjustl(model_name)))
      select case (mstr)
      case ("cpcm")
         model_tag = rad_type%cpcm
      case ("smd")
         model_tag = rad_type%smd
      case ("d3")
         model_tag = rad_type%d3
      case ("cosmo")
         model_tag = rad_type%cosmo
      case ("bondi")
         model_tag = rad_type%bondi
      case ("custom")
         call fatal_error(error, "Use new_radii_custom_atoms or new_radii_custom_elements for custom radii")
         return
      case default
         write (msg, '(a,a,a)') "Unknown radius type: '", trim(model_name), "'"
         call fatal_error(error, trim(msg))
         return
      end select

      call new_radii_int(model_tag, model, error, verbosity)
   end subroutine new_radii_str

end module moist_radii
