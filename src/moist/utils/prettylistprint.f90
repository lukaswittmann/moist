module moist_utils_prettylistprint
   use, intrinsic :: iso_fortran_env, only: output_unit, int8, int16, int32, int64, &
     & real32, real64
   implicit none
   private
   public :: prettylistprinter, new_prettylistprinter

   type :: prettylistprinter
      integer :: unit = output_unit
      integer :: offset = 1
      integer :: column_gap = 1
      integer :: ncols = 0
      integer :: next_col = 1
      integer :: fmt_len = 16
      integer, allocatable :: widths(:)
      character(:), allocatable :: headers(:)
      character(:), allocatable :: row(:)
      character(:), allocatable :: fmt_int
      character(:), allocatable :: fmt_real
      character(:), allocatable :: fmt_exp
      character(:), allocatable :: fmt_logical
   contains
      procedure :: header
      procedure :: print_header
      procedure :: separator
      procedure :: blank
      procedure :: begin_row
      procedure :: skip
      procedure :: end_row
      procedure :: set_column_gap
      procedure :: set_real_formats
      procedure, private :: add_i8
      procedure, private :: add_i16
      procedure, private :: add_i32
      procedure, private :: add_i64
      procedure, private :: add_r32
      procedure, private :: add_r64
      procedure, private :: add_l
      procedure, private :: add_c
      generic :: add => add_i8, add_i16, add_i32, add_i64, add_r32, add_r64, add_l, add_c
   end type prettylistprinter

contains

   !> Print a centered section header line using '=' fill.
   !> The line spans the full table width and respects left offset.
   !> @param[inout] self  Pretty list printer instance.
   !> @param[in]    title Section title text.
   subroutine header(self, title)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Section title
      character(*), intent(in) :: title
      character(:), allocatable :: spaced_title, block, line
      integer :: w, rem, nleft, nright, left_padding

      w = table_width(self)
      spaced_title = spread_text(trim(adjustl(title)))
      block = ' '//spaced_title//' '

      if (w <= 0) return

      if (len(block) >= w) then
         line = block(1:w)
      else
         rem = w - len(block)
         nleft = rem/2
         nright = rem - nleft
         line = repeat('=', nleft)//block//repeat('=', nright)
      end if

      left_padding = self%offset
      if (left_padding > 0) then
         write (self%unit, '(A)', advance='no') repeat(' ', left_padding)
      end if
      write (self%unit, '(A)') line
   end subroutine header

   !> Construct a pretty list printer with fixed column widths and headers.
   !> @param[in] widths   Width per column.
   !> @param[in] headers  Header text per column.
   !> @param[in] unit     Optional Fortran output unit.
   !> @param[in] offset   Optional left offset (spaces) for printed lines.
   !> @param[in] fmt_len  Optional base width for default number formats.
   !> @param[in] fmt_int  Optional integer format override.
   !> @param[in] fmt_real Optional fixed real format override.
   !> @param[in] fmt_exp  Optional exponential real format override.
   !> @param[in] fmt_logical Optional logical format override.
   !> @param[in] column_gap Optional spaces inserted between columns.
   function new_prettylistprinter(widths, headers, &
                                  unit, offset, fmt_len, fmt_int, fmt_real, fmt_exp, fmt_logical, column_gap) result(plp)
      !> Column widths
      integer, intent(in) :: widths(:)
      !> Column headers
      character(*), intent(in) :: headers(:)
      !> Optional output unit
      integer, intent(in), optional :: unit
      !> Optional left offset in spaces
      integer, intent(in), optional :: offset
      !> Optional format controls
      integer, intent(in), optional :: fmt_len
      character(*), intent(in), optional :: fmt_int, fmt_real, fmt_exp, fmt_logical
      !> Optional spacing between adjacent columns
      integer, intent(in), optional :: column_gap
      !> Constructed pretty list printer
      type(prettylistprinter) :: plp
      integer :: i, wmax, hmax

      if (size(widths) /= size(headers)) then
         ! TODO: Proper error propagration
         error stop 'prettylistprinter: widths and headers size mismatch'
      end if
      if (size(widths) == 0) then
         ! TODO: Proper error propagration
         error stop 'prettylistprinter: at least one column is required'
      end if

      plp%ncols = size(widths)
      allocate (plp%widths(plp%ncols))
      plp%widths = max(1, widths)
      hmax = 1
      do i = 1, plp%ncols
         hmax = max(hmax, len_trim(headers(i)))
      end do
      allocate (character(len=hmax) :: plp%headers(plp%ncols))
      do i = 1, plp%ncols
         plp%headers(i) = trim(headers(i))
      end do

      if (present(unit)) plp%unit = unit
      if (present(offset)) plp%offset = max(0, offset)
      if (present(column_gap)) plp%column_gap = max(0, column_gap)
      if (present(fmt_len)) plp%fmt_len = max(1, fmt_len)

      plp%fmt_int = int_fmt(plp%fmt_len)
      plp%fmt_real = fixed_fmt(plp%fmt_len, 4)
      plp%fmt_exp = exp_fmt(plp%fmt_len, 4)
      plp%fmt_logical = 'L1'

      if (present(fmt_int)) plp%fmt_int = trim(fmt_int)
      if (present(fmt_real)) plp%fmt_real = trim(fmt_real)
      if (present(fmt_exp)) plp%fmt_exp = trim(fmt_exp)
      if (present(fmt_logical)) plp%fmt_logical = trim(fmt_logical)

      wmax = maxval(plp%widths)
      allocate (character(len=wmax) :: plp%row(plp%ncols))
      call plp%begin_row()
   end function new_prettylistprinter

   !> Set default real formats after construction.
   !> @param[inout] self     Pretty list printer instance.
   !> @param[in]    fmt_real Optional fixed-point format.
   !> @param[in]    fmt_exp  Optional exponential format.
   subroutine set_real_formats(self, fmt_real, fmt_exp)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Optional fixed-point format override
      character(*), intent(in), optional :: fmt_real
      !> Optional exponential format override
      character(*), intent(in), optional :: fmt_exp

      if (present(fmt_real)) self%fmt_real = trim(fmt_real)
      if (present(fmt_exp)) self%fmt_exp = trim(fmt_exp)
   end subroutine set_real_formats

   !> Set the spacing inserted between adjacent columns.
   !> @param[inout] self       Pretty list printer instance.
   !> @param[in]    column_gap Number of spaces between columns.
   subroutine set_column_gap(self, column_gap)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Inter-column gap in spaces
      integer, intent(in) :: column_gap

      self%column_gap = max(0, column_gap)
   end subroutine set_column_gap

   !> Print column headers right-aligned in their fields.
   !> @param[inout] self Pretty list printer instance.
   subroutine print_header(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      integer :: i

      if (self%offset > 0) then
         write (self%unit, '(A)', advance='no') repeat(' ', self%offset)
      end if
      do i = 1, self%ncols
         write (self%unit, '(A)', advance='no') format_cell(self%headers(i), self%widths(i))
         if (i < self%ncols) write (self%unit, '(A)', advance='no') repeat(' ', self%column_gap)
      end do
      write (self%unit, *)
   end subroutine print_header

   !> Print a separator line with `width-1` dashes per column and configurable gaps.
   !> @param[inout] self Pretty list printer instance.
   subroutine separator(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      integer :: i

      write (self%unit, '(A)', advance='no') repeat(' ', self%offset)
      do i = 1, self%ncols
         write (self%unit, '(A)', advance='no') repeat('-', max(0, self%widths(i)))
         if (i < self%ncols) write (self%unit, '(A)', advance='no') repeat(' ', self%column_gap)
      end do
      write (self%unit, *)
   end subroutine separator

   !> Print a blank line.
   !> @param[inout] self Pretty printer instance.
   subroutine blank(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      write (self%unit, '(A)') ''
   end subroutine blank

   !> Start a new row and reset write position to first column.
   !> @param[inout] self Pretty list printer instance.
   subroutine begin_row(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      self%row = ''
      self%next_col = 1
   end subroutine begin_row

   !> Leave current column blank and move to the next one.
   !> @param[inout] self Pretty list printer instance.
   subroutine skip(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self

      call ensure_can_add(self)
      self%row(self%next_col) = ''
      self%next_col = self%next_col + 1
   end subroutine skip

   !> Print current row and reset for next row.
   !> @param[inout] self Pretty list printer instance.
   subroutine end_row(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      integer :: i

      if (self%next_col /= self%ncols + 1) then
         ! TODO: Proper error propagration
         error stop 'prettylistprinter: row has missing columns, use skip() or add()'
      end if

      if (self%offset > 0) then
         write (self%unit, '(A)', advance='no') repeat(' ', self%offset)
      end if
      do i = 1, self%ncols
         write (self%unit, '(A)', advance='no') format_cell(self%row(i), self%widths(i))
         if (i < self%ncols) write (self%unit, '(A)', advance='no') repeat(' ', self%column_gap)
      end do
      write (self%unit, *)

      call self%begin_row()
   end subroutine end_row

   !> Add an int8 value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_i8(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> int8 value
      integer(int8), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      eff_fmt = self%fmt_int
      if (present(fmt)) eff_fmt = trim(fmt)
      call add_from_string(self, value_to_string(val, eff_fmt))
   end subroutine add_i8

   !> Add an int16 value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_i16(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> int16 value
      integer(int16), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      eff_fmt = self%fmt_int
      if (present(fmt)) eff_fmt = trim(fmt)
      call add_from_string(self, value_to_string(val, eff_fmt))
   end subroutine add_i16

   !> Add an int32 value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_i32(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> int32 value
      integer(int32), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      eff_fmt = self%fmt_int
      if (present(fmt)) eff_fmt = trim(fmt)
      call add_from_string(self, value_to_string(val, eff_fmt))
   end subroutine add_i32

   !> Add an int64 value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_i64(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> int64 value
      integer(int64), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      eff_fmt = self%fmt_int
      if (present(fmt)) eff_fmt = trim(fmt)
      call add_from_string(self, value_to_string(val, eff_fmt))
   end subroutine add_i64

   !> Add a real32 value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_r32(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> real32 value
      real(real32), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(64) :: eff_fmt
      character(:), allocatable :: s
      integer :: icol, wcol

      if (present(fmt)) then
         eff_fmt = trim(fmt)
      else
         eff_fmt = default_real_fmt(self, real(val, kind=real64))
      end if
      s = value_to_string(val, trim(eff_fmt))
      call ensure_can_add(self)
      icol = self%next_col
      wcol = self%widths(icol)
      if (is_real_overflow(s, wcol)) then
         s = overflow_marker(wcol, val < 0.0_real32)
      end if
      self%row(icol) = trim(s)
      self%next_col = icol + 1
   end subroutine add_r32

   !> Add a real64 value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_r64(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> real64 value
      real(real64), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(64) :: eff_fmt
      character(:), allocatable :: s
      integer :: icol, wcol

      if (present(fmt)) then
         eff_fmt = trim(fmt)
      else
         eff_fmt = default_real_fmt(self, val)
      end if
      s = value_to_string(val, trim(eff_fmt))
      call ensure_can_add(self)
      icol = self%next_col
      wcol = self%widths(icol)
      if (is_real_overflow(s, wcol)) then
         s = overflow_marker(wcol, val < 0.0_real64)
      end if
      self%row(icol) = trim(s)
      self%next_col = icol + 1
   end subroutine add_r64

   !> Add a logical value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_l(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Logical value
      logical, intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      eff_fmt = self%fmt_logical
      if (present(fmt)) eff_fmt = trim(fmt)
      call add_from_string(self, value_to_string(val, eff_fmt))
   end subroutine add_l

   !> Add a character value to current row.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    val  Value to insert.
   !> @param[in]    fmt  Optional format override for this cell.
   subroutine add_c(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Character value
      character(*), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt

      if (present(fmt)) then
         call add_from_string(self, value_to_string(val, trim(fmt)))
      else
         call add_from_string(self, trim(val))
      end if
   end subroutine add_c

   !> Add pre-formatted string content to the current column.
   !> @param[inout] self Pretty list printer instance.
   !> @param[in]    s    Pre-formatted cell text.
   subroutine add_from_string(self, s)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Cell text
      character(*), intent(in) :: s

      call ensure_can_add(self)
      self%row(self%next_col) = trim(s)
      self%next_col = self%next_col + 1
   end subroutine add_from_string

   !> Ensure the current row still has space for one more value.
   !> @param[in] self Pretty list printer instance.
   subroutine ensure_can_add(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(in) :: self
      if (self%next_col > self%ncols) then
         ! TODO: Proper error propagration
         error stop 'prettylistprinter: too many values in row'
      end if
   end subroutine ensure_can_add

   !> Fit and right-align a cell value into fixed-width output.
   !> @param[in] s     Source text.
   !> @param[in] width Cell width.
   function format_cell(s, width) result(out)
      !> Source text
      character(*), intent(in) :: s
      !> Cell width
      integer, intent(in) :: width
      !> Right-aligned output cell
      character(:), allocatable :: out
      integer :: ls

      if (width <= 0) then
         out = ''
         return
      end if

      ls = len_trim(s)
      if (ls > width) then
         out = s(ls - width + 1:ls)
      else
         out = repeat(' ', width - ls)//s(1:ls)
      end if
   end function format_cell

   !> Convert supported scalar values to string using supplied format.
   !> @param[in] val Scalar value.
   !> @param[in] fmt Fortran format string without outer parentheses.
   function value_to_string(val, fmt) result(s)
      class(*), intent(in) :: val
      character(*), intent(in) :: fmt
      character(:), allocatable :: s
      character(256) :: buf

      buf = ''

      select type (val)
      type is (integer(int8))
         write (buf, '('//trim(fmt)//')') val
      type is (integer(int16))
         write (buf, '('//trim(fmt)//')') val
      type is (integer(int32))
         write (buf, '('//trim(fmt)//')') val
      type is (integer(int64))
         write (buf, '('//trim(fmt)//')') val
      type is (real(real32))
         if (val == 0.0_real32) then
            s = zero_value_string(fmt)
            return
         end if
         write (buf, '('//trim(fmt)//')') val
      type is (real(real64))
         if (val == 0.0_real64) then
            s = zero_value_string(fmt)
            return
         end if
         write (buf, '('//trim(fmt)//')') val
      type is (logical)
         write (buf, '('//trim(fmt)//')') val
      type is (character(*))
         buf = val
      class default
         buf = '<unsupported type>'
      end select

      s = trim(buf)
   end function value_to_string

   !> Return canonical zero representation based on supplied format width.
   !> @param[in] fmt Fortran format string without outer parentheses.
   function zero_value_string(fmt) result(s)
      character(*), intent(in) :: fmt
      character(:), allocatable :: s
      character(256) :: buf
      integer :: idot, w

      write (buf, '('//trim(fmt)//')') 0.0_real64
      idot = index(buf, '.')
      w = len_trim(buf)

      if (idot > 1 .and. w > 0) then
         s = repeat(' ', idot - 2)//'0.0'//repeat(' ', max(0, w - (idot + 1)))
      else
         s = '0.0'
      end if
   end function zero_value_string

   !> Build fixed real format string.
   !> @param[in] width    Total field width.
   !> @param[in] decimals Digits after decimal point.
   function fixed_fmt(width, decimals) result(fmt)
      integer, intent(in) :: width, decimals
      character(:), allocatable :: fmt
      character(32) :: wbuf, dbuf

      write (wbuf, '(I0)') max(1, width)
      write (dbuf, '(I0)') max(0, decimals)
      fmt = 'F'//trim(wbuf)//'.'//trim(dbuf)
   end function fixed_fmt

   !> Build exponential real format string.
   !> @param[in] width    Total field width.
   !> @param[in] decimals Digits after decimal point.
   function exp_fmt(width, decimals) result(fmt)
      integer, intent(in) :: width, decimals
      character(:), allocatable :: fmt
      character(32) :: wbuf, dbuf

      write (wbuf, '(I0)') max(1, width)
      write (dbuf, '(I0)') max(0, decimals)
      fmt = 'ES'//trim(wbuf)//'.'//trim(dbuf)
   end function exp_fmt

   !> Build integer format string.
   !> @param[in] width Base width used to derive integer field width.
   function int_fmt(width) result(fmt)
      integer, intent(in) :: width
      character(:), allocatable :: fmt
      character(32) :: wbuf

      write (wbuf, '(I0)') max(1, width - 7)
      fmt = 'I'//trim(wbuf)
   end function int_fmt

   !> Select default real format from value magnitude.
   !> @param[in] self Pretty list printer instance.
   !> @param[in] val  Real64 value.
   function default_real_fmt(self, val) result(fmt)
      class(prettylistprinter), intent(in) :: self
      real(real64), intent(in) :: val
      character(:), allocatable :: fmt
      real(real64) :: aval

      aval = abs(val)
      if (aval == 0.0_real64) then
         fmt = self%fmt_real
      else if (aval < 1.0e-6_real64 .or. aval >= 1.0e6_real64) then
         fmt = self%fmt_exp
      else
         fmt = self%fmt_real
      end if
   end function default_real_fmt

   !> Determine whether formatted real text overflows the target cell width.
   !> @param[in] s      Formatted real text.
   !> @param[in] width  Target cell width.
   function is_real_overflow(s, width) result(overflow)
      character(*), intent(in) :: s
      integer, intent(in) :: width
      logical :: overflow

      overflow = (index(s, '*') > 0) .or. (len_trim(s) - 1 > width)
   end function is_real_overflow

   !> Create overflow marker text for a real cell.
   !> Positive overflow uses '+' and negative overflow uses '-'.
   !> Marker length is `width-1` as requested.
   !> @param[in] width       Target cell width.
   !> @param[in] is_negative Sign selector.
   function overflow_marker(width, is_negative) result(s)
      integer, intent(in) :: width
      logical, intent(in) :: is_negative
      character(:), allocatable :: s
      integer :: n

      n = max(0, width)
      if (is_negative) then
         s = repeat('-', n)
      else
         s = repeat('+', n)
      end if
   end function overflow_marker

   !> Compute total printable table width, including inter-column spaces.
   !> @param[in] self Pretty list printer instance.
   function table_width(self) result(w)
      class(prettylistprinter), intent(in) :: self
      integer :: w

      if (self%ncols <= 0) then
         w = 0
      else
         w = sum(self%widths) + self%column_gap*(self%ncols - 1)
      end if
   end function table_width

   !> Insert one space between all characters of input text.
   !> @param[in] s Input text.
   function spread_text(s) result(out)
      character(*), intent(in) :: s
      character(:), allocatable :: out
      integer :: i, ls

      ls = len_trim(s)
      if (ls <= 0) then
         out = ''
         return
      end if

      out = s(1:1)
      do i = 2, ls
         out = out//' '//s(i:i)
      end do
   end function spread_text

end module moist_utils_prettylistprint
