
!> License text for moist
module moist_output_license
   implicit none(type, external)
   private

   public :: print_license

contains

   !> Print MPL 2.0 license notice.
   subroutine print_license(unit)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit

      write (unit, "(a)") &
         "This Source Code Form is subject to the terms of the Mozilla Public", &
         "License, v. 2.0. If a copy of the MPL was not distributed with this", &
         "file, You can obtain one at https://mozilla.org/MPL/2.0/.", &
         "", &
         "The full license text is distributed in the LICENSE file."

   end subroutine print_license

end module moist_output_license
