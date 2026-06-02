
!> License text for moist.
module moist_output_license
   implicit none
   private

   public :: print_license

contains

   !> Print LGPL license notice.
   subroutine print_license(unit)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit

      write (unit, '(a)') &
         "moist is free software: you can redistribute it and/or modify it under", &
         "the terms of the Lesser GNU General Public License as published by", &
         "the Free Software Foundation, either version 3 of the License, or", &
         "(at your option) any later version.", &
         "", &
         "moist is distributed in the hope that it will be useful,", &
         "but WITHOUT ANY WARRANTY; without even the implied warranty of", &
         "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the", &
         "Lesser GNU General Public License for more details.", &
         "", &
         "You should have received a copy of the Lesser GNU General Public License", &
         "along with moist.  If not, see <https://www.gnu.org/licenses/>."

   end subroutine print_license

end module moist_output_license
