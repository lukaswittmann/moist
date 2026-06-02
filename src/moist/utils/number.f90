
module moist_utils_number

   use mctc_env, only: wp
   use ieee_arithmetic, only: ieee_is_nan

   implicit none

   public :: is_exceptional

contains

   elemental function is_exceptional(val)
      real(wp), intent(in) :: val
      logical :: is_exceptional
      is_exceptional = ieee_is_nan(val) .or. abs(val) > huge(val)
   end function is_exceptional

end module moist_utils_number
