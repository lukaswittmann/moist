!> Scalar-to-text conversion shared by the pretty printers
!>
!> - `value_to_string`: one scalar through a Fortran edit descriptor, zero as `0.0`
!> - `fixed_fmt`, `exp_fmt`, `int_fmt`: default descriptors from a field width
!> - `magnitude_real_fmt`: fixed or exponential descriptor by the size of a value
module moist_utils_format
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64
   implicit none(type, external)
   private
   public :: value_to_string, fixed_fmt, exp_fmt, int_fmt, magnitude_real_fmt

   !> Internal-write buffer, room for a numeric edit re-rendered at a 400-character
   !> width and literal text around it
   integer, parameter :: buffer_len = 512

contains

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
