module moist_utils_string

   implicit none

   public :: lowercase

contains

   !> Convert string to lower case
   pure function lowercase(str) result(lcstr)
      character(len=*), intent(in)  :: str
      character(len=len_trim(str)) :: lcstr
      integer :: ilen, ioffset, iquote, i, iav, iqc

      ilen = len_trim(str)
      ioffset = iachar('A') - iachar('a')
      iquote = 0
      lcstr = str
      do i = 1, ilen
         iav = iachar(str(i:i))
         if (iquote == 0 .and. (iav == 34 .or. iav == 39)) then
            iquote = 1
            iqc = iav
            cycle
         end if
         if (iquote == 1 .and. iav == iqc) then
            iquote = 0
            cycle
         end if
         if (iquote == 1) cycle
         if (iav >= iachar('A') .and. iav <= iachar('Z')) then
            lcstr(i:i) = achar(iav - ioffset)
         else
            lcstr(i:i) = str(i:i)
         end if
      end do

   end function lowercase

end module moist_utils_string
