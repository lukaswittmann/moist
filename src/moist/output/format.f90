module moist_output_format
   use mctc_env, only: wp

   implicit none
   private

   public :: format_string, getline, print_wrapped

   interface format_string
      module procedure :: format_string_int
      module procedure :: format_string_real_dp
   end interface format_string

contains

   pure function format_string_real_dp(val, format) result(str)
      real(wp), intent(in) :: val
      character(len=*), intent(in) :: format
      character(len=:), allocatable :: str

      character(len=128) :: buffer
      integer :: stat

      write (buffer, format, iostat=stat) val
      if (stat == 0) then
         str = trim(buffer)
      else
         str = "*"
      end if
   end function format_string_real_dp

   pure function format_string_int(val, format) result(str)
      integer, intent(in) :: val
      character(len=*), intent(in) :: format
      character(len=:), allocatable :: str

      character(len=128) :: buffer
      integer :: stat

      write (buffer, format, iostat=stat) val
      if (stat == 0) then
         str = trim(buffer)
      else
         str = "*"
      end if
   end function format_string_int

   !> reads a line from unit into an allocatable character
   subroutine getline(unit, line, iostat)
      integer, intent(in) :: unit
      character(len=:), allocatable, intent(out) :: line
      integer, intent(out), optional :: iostat

      integer, parameter  :: buffersize = 256
      character(len=buffersize) :: buffer
      integer :: size
      integer :: stat

      line = ''
      do
         read (unit, '(a)', advance='no', iostat=stat, size=size)  &
         &    buffer
         if (stat > 0) then
            if (present(iostat)) iostat = stat
            return ! an error occurred
         end if
         line = line//buffer(:size)
         if (stat < 0) then
            if (is_iostat_eor(stat)) stat = 0
            if (present(iostat)) iostat = stat
            return
         end if
      end do

   end subroutine getline

   !> Print a string with word-wrapping at a given width.
   !> Words are never split; breaks occur at spaces only.
   !> @param[in] unit   Fortran I/O unit
   !> @param[in] text   String to print
   !> @param[in] indent Prefix for every line (e.g. "  ")
   !> @param[in] width  Maximum characters per line (excluding indent)
   subroutine print_wrapped(unit, text, indent, width)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: text
      character(len=*), intent(in) :: indent
      integer, intent(in) :: width

      integer :: pos, line_start, last_space, text_len

      text_len = len_trim(text)
      if (text_len == 0) then
         write (unit, '(a)') indent
         return
      end if

      line_start = 1
      do while (line_start <= text_len)
         ! If remainder fits on one line, print it and exit
         if (line_start + width - 1 >= text_len) then
            write (unit, '(a,a)') indent, text(line_start:text_len)
            return
         end if

         ! Find the last space within the allowed width
         last_space = 0
         do pos = line_start, min(line_start + width - 1, text_len)
            if (text(pos:pos) == ' ') last_space = pos
         end do

         if (last_space > line_start) then
            ! Break at the last space within width
            write (unit, '(a,a)') indent, text(line_start:last_space - 1)
            line_start = last_space + 1
         else
            ! No space found within width - find the next space beyond width
            last_space = index(text(line_start:text_len), ' ')
            if (last_space > 0) then
               last_space = line_start + last_space - 1
               write (unit, '(a,a)') indent, text(line_start:last_space - 1)
               line_start = last_space + 1
            else
               ! No more spaces at all - print the rest
               write (unit, '(a,a)') indent, text(line_start:text_len)
               return
            end if
         end if
      end do
   end subroutine print_wrapped

end module moist_output_format
