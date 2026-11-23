!> Scalar-to-text conversion and text-wrapping utilities
!>
!> - `to_string`: inline scalar rendering, trimmed; `(i0)` and `(g0)` without a format
!> - `print_wrapped`: word-wrapped text output at a given column width
!> - `value_to_string`: polymorphic scalar for pretty printers, zero as `0.0`
!> - `fixed_fmt`, `exp_fmt`, `int_fmt`: default descriptors from a field width
!> - `magnitude_real_fmt`: fixed or exponential descriptor by the size of a value
!>
!> No deferred-length character results: gfortran keeps the hidden result
!> length of such a call in static storage, even under `-fopenmp`, so threads
module moist_utils_format
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
   implicit none(type, external)
   private
   public :: to_string, print_wrapped
   public :: value_to_string, fixed_fmt, exp_fmt, int_fmt, magnitude_real_fmt

   !> Format a scalar into a trimmed, left-adjusted string
   interface to_string
      module procedure :: to_string_int
      module procedure :: to_string_int64
      module procedure :: to_string_real_dp
      module procedure :: to_string_real_dp_fmt
   end interface to_string

   !> Internal-write buffer, room for a numeric edit re-rendered at a 400-character
   !> width and literal text around it
   integer, parameter :: buffer_len = 512

contains

   !> Render an integer as `(i0)` does, without an internal write
   !>
   !> @param[in] val Value to render
   pure function to_string_int(val) result(str)
      !> Value to render
      integer(int32), intent(in) :: val
      !> Rendering
      character(len=int_digits(int(val, int64))) :: str

      str = to_string_int64(int(val, int64))
   end function to_string_int

   !> Render a 64-bit integer as `(i0)` does, without an internal write
   !>
   !> @param[in] val Value to render
   pure function to_string_int64(val) result(str)
      !> Value to render
      integer(int64), intent(in) :: val
      !> Rendering
      character(len=int_digits(val)) :: str

      integer(int64) :: rest
      integer :: pos

      ! Digits right to left; a negative value keeps its sign in `mod`
      rest = val
      do pos = len(str), 1, -1
         str(pos:pos) = achar(iachar("0") + int(abs(mod(rest, 10_int64))))
         rest = rest/10
         if (rest == 0) exit
      end do
      if (val < 0) str(1:1) = "-"
   end function to_string_int64

   !> Render a double-precision real as `(g0)` does
   !>
   !> @param[in] val Value to render
   pure function to_string_real_dp(val) result(str)
      !> Value to render
      real(real64), intent(in) :: val
      !> Rendering
      character(len=written_len_real_dp(val, "(g0)")) :: str

      character(len=buffer_len) :: buffer

      call write_real_dp(val, "(g0)", buffer)
      str = buffer
   end function to_string_real_dp

   !> Render a double-precision real with a parenthesized Fortran format
   !>
   !> - leading and trailing blanks are stripped
   !> - `"*"` if the internal write fails
   !>
   !> @param[in] val    Value to render
   !> @param[in] format Fortran format with outer parentheses, e.g. `"(es9.2)"`
   pure function to_string_real_dp_fmt(val, format) result(str)
      !> Value to render
      real(real64), intent(in) :: val
      !> Fortran format with outer parentheses
      character(len=*), intent(in) :: format
      !> Rendering
      character(len=written_len_real_dp(val, format)) :: str

      character(len=buffer_len) :: buffer

      call write_real_dp(val, format, buffer)
      str = buffer
   end function to_string_real_dp_fmt

   !> Characters of an integer under `(i0)`, sign included
   !>
   !> @param[in] val Value
   pure function int_digits(val) result(n)
      !> Value
      integer(int64), intent(in) :: val
      !> Number of characters
      integer :: n

      integer(int64) :: rest

      n = merge(2, 1, val < 0)
      rest = val/10
      do while (rest /= 0)
         n = n + 1
         rest = rest/10
      end do
   end function int_digits

   !> Length of a real rendering with a format, blanks stripped
   !>
   !> @param[in] val    Value
   !> @param[in] format Fortran format with outer parentheses
   pure function written_len_real_dp(val, format) result(n)
      !> Value
      real(real64), intent(in) :: val
      !> Fortran format with outer parentheses
      character(len=*), intent(in) :: format
      !> Number of characters
      integer :: n

      character(len=buffer_len) :: buffer

      call write_real_dp(val, format, buffer)
      n = len_trim(buffer)
   end function written_len_real_dp

   !> Left-adjusted internal write of a real, `"*"` on failure
   !>
   !> @param[in]  val    Value
   !> @param[in]  format Fortran format with outer parentheses
   !> @param[out] buffer Rendering
   pure subroutine write_real_dp(val, format, buffer)
      !> Value
      real(real64), intent(in) :: val
      !> Fortran format with outer parentheses
      character(len=*), intent(in) :: format
      !> Rendering
      character(len=buffer_len), intent(out) :: buffer

      integer :: stat

      write (buffer, format, iostat=stat) val
      if (stat == 0) then
         buffer = adjustl(buffer)
      else
         buffer = "*"
      end if
   end subroutine write_real_dp

   !> Write text word-wrapped at spaces to a given width
   !>
   !> - words are never split; breaks occur at spaces only
   !>
   !> @param[in] unit    Fortran I/O unit
   !> @param[in] text    text to write
   !> @param[in] indent  prefix for every line, e.g. `"  "`
   !> @param[in] width   maximum characters per line, excluding indent
   subroutine print_wrapped(unit, text, indent, width)
      !> Fortran I/O unit
      integer, intent(in) :: unit
      !> Text to write
      character(len=*), intent(in) :: text
      !> Prefix for every line
      character(len=*), intent(in) :: indent
      !> Maximum characters per line, excluding indent
      integer, intent(in) :: width

      integer :: pos, line_start, last_space, text_len

      text_len = len_trim(text)
      if (text_len == 0) then
         write (unit, "(a)") indent
         return
      end if

      line_start = 1
      do while (line_start <= text_len)
         ! Remainder fits on one line
         if (line_start + width - 1 >= text_len) then
            write (unit, "(a,a)") indent, text(line_start:text_len)
            return
         end if

         ! Last space within the allowed width
         last_space = 0
         do pos = line_start, min(line_start + width - 1, text_len)
            if (text(pos:pos) == " ") last_space = pos
         end do

         if (last_space > line_start) then
            ! Break at the last space within width
            write (unit, "(a,a)") indent, text(line_start:last_space - 1)
            line_start = last_space + 1
         else
            ! No space within width: break at the next space beyond width
            last_space = index(text(line_start:text_len), " ")
            if (last_space > 0) then
               last_space = line_start + last_space - 1
               write (unit, "(a,a)") indent, text(line_start:last_space - 1)
               line_start = last_space + 1
            else
               ! No remaining spaces: write the rest
               write (unit, "(a,a)") indent, text(line_start:text_len)
               return
            end if
         end if
      end do
   end subroutine print_wrapped

   !> Convert a scalar to text with an edit descriptor
   !>
   !> - integers of every kind, `real32`, `real64`, logical and character
   !> - a real zero gives `0.0`, placed where the descriptor puts the units digit
   !> - character values are copied, `fmt` unused
   !>
   !> @param[in]  val Scalar value
   !> @param[in]  fmt Fortran format without outer parentheses
   !> @param[out] s   Formatted text
   subroutine value_to_string(val, fmt, s)
      !> Scalar value
      class(*), intent(in) :: val
      !> Fortran format without outer parentheses
      character(*), intent(in) :: fmt
      !> Formatted text
      character(:), allocatable, intent(out) :: s
      character(buffer_len) :: buf

      buf = ""

      select type (val)
      type is (integer(int8))
         write (buf, "("//trim(fmt)//")") val
      type is (integer(int16))
         write (buf, "("//trim(fmt)//")") val
      type is (integer(int32))
         write (buf, "("//trim(fmt)//")") val
      type is (integer(int64))
         write (buf, "("//trim(fmt)//")") val
      type is (real(real32))
         if (val == 0.0_real32) then
            call zero_value_string(fmt, s)
            return
         end if
         write (buf, "("//trim(fmt)//")") val
      type is (real(real64))
         if (val == 0.0_real64) then
            call zero_value_string(fmt, s)
            return
         end if
         write (buf, "("//trim(fmt)//")") val
      type is (logical)
         write (buf, "("//trim(fmt)//")") val
      type is (character(*))
         buf = val
      class default
         error stop "moist format: unsupported value type"
      end select

      s = trim(buf)
   end subroutine value_to_string

   !> Zero as `0.0`, its units digit where the descriptor puts that of a zero
   !>
   !> @param[in]  fmt Fortran format without outer parentheses
   !> @param[out] s   Formatted zero
   subroutine zero_value_string(fmt, s)
      !> Fortran format without outer parentheses
      character(*), intent(in) :: fmt
      !> Formatted zero
      character(:), allocatable, intent(out) :: s
      character(buffer_len) :: buf
      integer :: idot, w

      write (buf, "("//trim(fmt)//")") 0.0_real64
      idot = index(buf, ".")
      w = len_trim(buf)

      if (idot > 1 .and. w > 0) then
         s = repeat(" ", idot - 2)//"0.0"//repeat(" ", max(0, w - (idot + 1)))
      else
         s = "0.0"
      end if
   end subroutine zero_value_string

   !> Fixed-point descriptor `Fw.d`
   !>
   !> @param[in]  width    Field width, at least 1
   !> @param[in]  decimals Digits after the point, at least 0
   !> @param[out] fmt      Descriptor
   subroutine fixed_fmt(width, decimals, fmt)
      !> Field width and digits after the point
      integer, intent(in) :: width, decimals
      !> Descriptor
      character(:), allocatable, intent(out) :: fmt
      character(32) :: wbuf, dbuf

      write (wbuf, "(I0)") max(1, width)
      write (dbuf, "(I0)") max(0, decimals)
      fmt = "F"//trim(wbuf)//"."//trim(dbuf)
   end subroutine fixed_fmt

   !> Scientific descriptor `ESw.d`
   !>
   !> @param[in]  width    Field width, at least 1
   !> @param[in]  decimals Digits after the point, at least 0
   !> @param[out] fmt      Descriptor
   subroutine exp_fmt(width, decimals, fmt)
      !> Field width and digits after the point
      integer, intent(in) :: width, decimals
      !> Descriptor
      character(:), allocatable, intent(out) :: fmt
      character(32) :: wbuf, dbuf

      write (wbuf, "(I0)") max(1, width)
      write (dbuf, "(I0)") max(0, decimals)
      fmt = "ES"//trim(wbuf)//"."//trim(dbuf)
   end subroutine exp_fmt

   !> Integer descriptor `Iw`, seven characters narrower than a real field
   !>
   !> @param[in]  width Base field width of the reals
   !> @param[out] fmt   Descriptor
   subroutine int_fmt(width, fmt)
      !> Base field width of the reals
      integer, intent(in) :: width
      !> Descriptor
      character(:), allocatable, intent(out) :: fmt
      character(32) :: wbuf

      write (wbuf, "(I0)") max(1, width - 7)
      fmt = "I"//trim(wbuf)
   end subroutine int_fmt

   !> Fixed-point descriptor for moderate values, exponential below 1e-6 or from 1e6
   !>
   !> @param[in]  val      Value to print
   !> @param[in]  fmt_real Fixed-point descriptor, also for zero
   !> @param[in]  fmt_exp  Exponential descriptor
   !> @param[out] fmt      Selected descriptor
   subroutine magnitude_real_fmt(val, fmt_real, fmt_exp, fmt)
      !> Value to print
      real(real64), intent(in) :: val
      !> Fixed-point and exponential descriptors
      character(*), intent(in) :: fmt_real, fmt_exp
      !> Selected descriptor
      character(:), allocatable, intent(out) :: fmt
      real(real64) :: aval

      aval = abs(val)
      if (aval == 0.0_real64) then
         fmt = fmt_real
      else if (aval < 1.0e-6_real64 .or. aval >= 1.0e6_real64) then
         fmt = fmt_exp
      else
         fmt = fmt_real
      end if
   end subroutine magnitude_real_fmt

end module moist_utils_format
