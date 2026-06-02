
!> Versioning information on this library.
module moist_version
   implicit none
   private

   public :: moist_version_string, moist_version_compact
   public :: get_moist_version

   !> String representation of the moist version
   character(len=*), parameter :: moist_version_string = "0.6.0-alpha.1"

   !> Numeric representation of the moist version
   integer, parameter :: moist_version_compact(3) = [0, 6, 0]

contains

!> Getter function to retrieve moist version
   subroutine get_moist_version(major, minor, patch, string)

      !> Major version number of the moist version
      integer, intent(out), optional :: major

      !> Minor version number of the moist version
      integer, intent(out), optional :: minor

      !> Patch version number of the moist version
      integer, intent(out), optional :: patch

      !> String representation of the moist version
      character(len=:), allocatable, intent(out), optional :: string

      if (present(major)) then
         major = moist_version_compact(1)
      end if
      if (present(minor)) then
         minor = moist_version_compact(2)
      end if
      if (present(patch)) then
         patch = moist_version_compact(3)
      end if
      if (present(string)) then
         string = moist_version_string
      end if

   end subroutine get_moist_version

end module moist_version
