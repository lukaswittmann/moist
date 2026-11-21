!> Fixed-width tabular output
!>
!> - every cell (number, logical, string, column header) is right-aligned in its column
!> - text wider than its column is printed in full, never cut; the row then runs
!>   past the column and later cells shift right
!> - a numeric edit that overflows is repeated with the same descriptor at a large
!>   width, so decimals and exponent style stay and every integer digit appears
module moist_utils_prettylistprint
   use, intrinsic :: iso_fortran_env, only: output_unit, int8, int16, int32, int64, &
     & real32, real64
   use moist_utils_format, only: value_to_string, fixed_fmt, exp_fmt, int_fmt, magnitude_real_fmt
   implicit none(type, external)
   private
   public :: prettylistprinter, new_prettylistprinter

   !> Field width for re-rendering an overflowing numeric edit: room for the 309
   !> integer digits of huge(1.0_real64) under an F edit, sign, point and decimals;
   !> the buffer of [[value_to_string]] is longer, for literal text around the field
   integer, parameter :: wide_width = 400

   type :: prettylistprinter
      integer :: unit = output_unit
      integer :: offset = 1
      integer :: column_gap = 1
      integer :: ncols = 0
      integer :: next_col = 1
      integer :: fmt_len = 16
      integer, allocatable :: widths(:)
      character(:), allocatable :: headers(:)
      !> Text of the current row: cells and gaps added so far, without the offset
      character(:), allocatable :: row
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

   !> Print a centered section header line using '=' fill
   !>
   !> @param[inout] self  Pretty list printer instance
   !> @param[in]    title Section title text
   subroutine header(self, title)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Section title
      character(*), intent(in) :: title
      character(:), allocatable :: block, line
      integer :: w, rem, nleft, nright, left_padding

      w = table_width(self)
      block = " "//trim(adjustl(title))//" "

      if (w <= 0) return

      if (len(block) >= w) then
         line = block(1:w)
      else
         rem = w - len(block)
         nleft = rem/2
         nright = rem - nleft
         line = repeat("=", nleft)//block//repeat("=", nright)
      end if

      left_padding = self%offset
      if (left_padding > 0) then
         write (self%unit, "(A)", advance="no") repeat(" ", left_padding)
      end if
      write (self%unit, "(A)") line
   end subroutine header

   !> Construct a pretty list printer with fixed column widths and headers
   !>
   !> @param[in] widths   Width per column
   !> @param[in] headers  Header text per column
   !> @param[in] unit     Optional Fortran output unit
   !> @param[in] offset   Optional left offset (spaces) for printed lines
   !> @param[in] fmt_len  Optional base width for default number formats
   !> @param[in] fmt_int  Optional integer format override
   !> @param[in] fmt_real Optional fixed real format override
   !> @param[in] fmt_exp  Optional exponential real format override
   !> @param[in] fmt_logical Optional logical format override
   !> @param[in] column_gap Optional spaces inserted between columns
   function new_prettylistprinter(widths, headers, unit, offset, fmt_len, fmt_int, fmt_real, &
                                  fmt_exp, fmt_logical, column_gap) result(plp)
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
      integer :: i, hmax

      ! A malformed column specification is a programming error at the call
      ! site, not a runtime condition: there is no table to print and no
      ! sensible substitute for one
      ! TODO: add errorprop
      if (size(widths) /= size(headers)) then
         error stop "prettylistprinter: widths and headers size mismatch"
      end if
      ! TODO: add errorprop
      if (size(widths) == 0) then
         error stop "prettylistprinter: at least one column is required"
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

      call int_fmt(plp%fmt_len, plp%fmt_int)
      call fixed_fmt(plp%fmt_len, 4, plp%fmt_real)
      call exp_fmt(plp%fmt_len, 4, plp%fmt_exp)
      plp%fmt_logical = "L1"

      if (present(fmt_int)) plp%fmt_int = trim(fmt_int)
      if (present(fmt_real)) plp%fmt_real = trim(fmt_real)
      if (present(fmt_exp)) plp%fmt_exp = trim(fmt_exp)
      if (present(fmt_logical)) plp%fmt_logical = trim(fmt_logical)

      call plp%begin_row()
   end function new_prettylistprinter

   !> Set default real formats after construction
   !>
   !> @param[inout] self     Pretty list printer instance
   !> @param[in]    fmt_real Optional fixed-point format
   !> @param[in]    fmt_exp  Optional exponential format
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

   !> Set the spacing inserted between adjacent columns
   !>
   !> @param[inout] self       Pretty list printer instance
   !> @param[in]    column_gap Number of spaces between columns
   subroutine set_column_gap(self, column_gap)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Inter-column gap in spaces
      integer, intent(in) :: column_gap

      self%column_gap = max(0, column_gap)
   end subroutine set_column_gap

   !> Print column headers right-aligned in their fields
   !>
   !> - a header wider than its column is printed in full
   !>
   !> @param[inout] self Pretty list printer instance
   subroutine print_header(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      integer :: i
      character(:), allocatable :: cell

      if (self%offset > 0) then
         write (self%unit, "(A)", advance="no") repeat(" ", self%offset)
      end if
      do i = 1, self%ncols
         call format_cell(self%headers(i), self%widths(i), cell)
         write (self%unit, "(A)", advance="no") cell
         if (i < self%ncols) write (self%unit, "(A)", advance="no") repeat(" ", self%column_gap)
      end do
      write (self%unit, "(A)") ""
   end subroutine print_header

   !> Print a separator line with one dash per column character and configurable gaps
   !>
   !> @param[inout] self Pretty list printer instance
   subroutine separator(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      integer :: i

      write (self%unit, "(A)", advance="no") repeat(" ", self%offset)
      do i = 1, self%ncols
         write (self%unit, "(A)", advance="no") repeat("-", max(0, self%widths(i)))
         if (i < self%ncols) write (self%unit, "(A)", advance="no") repeat(" ", self%column_gap)
      end do
      write (self%unit, "(A)") ""
   end subroutine separator

   !> Print a blank line, then flush the unit so the finished block is written out
   !>
   !> @param[inout] self Pretty printer instance
   subroutine blank(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self

      write (self%unit, "(A)") ""
      flush (self%unit)
   end subroutine blank

   !> Start a new row and reset write position to first column
   !>
   !> @param[inout] self Pretty list printer instance
   subroutine begin_row(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self

      self%row = ""
      self%next_col = 1
   end subroutine begin_row

   !> Leave current column blank and move to the next one
   !>
   !> @param[inout] self Pretty list printer instance
   subroutine skip(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self

      call add_from_string(self, "")
   end subroutine skip

   !> Print current row and reset for next row
   !>
   !> @param[inout] self Pretty list printer instance
   subroutine end_row(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self

      ! Refuse a row that is short of its column count: a truncated line in
      ! the middle of a table is harder to diagnose than a stop right here
      ! TODO: add errorprop
      if (self%next_col /= self%ncols + 1) then
         error stop "prettylistprinter: row has missing columns, use skip() or add()"
      end if

      write (self%unit, "(A)") repeat(" ", self%offset)//self%row

      call self%begin_row()
   end subroutine end_row

   !> Add an int8 value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
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
      call add_number(self, val, eff_fmt)
   end subroutine add_i8

   !> Add an int16 value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
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
      call add_number(self, val, eff_fmt)
   end subroutine add_i16

   !> Add an int32 value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
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
      call add_number(self, val, eff_fmt)
   end subroutine add_i32

   !> Add an int64 value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
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
      call add_number(self, val, eff_fmt)
   end subroutine add_i64

   !> Add a real32 value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
   subroutine add_r32(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> real32 value
      real(real32), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      if (present(fmt)) then
         eff_fmt = trim(fmt)
      else
         call magnitude_real_fmt(real(val, kind=real64), self%fmt_real, self%fmt_exp, eff_fmt)
      end if
      call add_number(self, val, eff_fmt)
   end subroutine add_r32

   !> Add a real64 value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
   subroutine add_r64(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> real64 value
      real(real64), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt

      if (present(fmt)) then
         eff_fmt = trim(fmt)
      else
         call magnitude_real_fmt(val, self%fmt_real, self%fmt_exp, eff_fmt)
      end if
      call add_number(self, val, eff_fmt)
   end subroutine add_r64

   !> Add a logical value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
   subroutine add_l(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Logical value
      logical, intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: eff_fmt, text

      eff_fmt = self%fmt_logical
      if (present(fmt)) eff_fmt = trim(fmt)
      call value_to_string(val, eff_fmt, text)
      call add_from_string(self, text)
   end subroutine add_l

   !> Add a character value to current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Value to insert
   !> @param[in]    fmt  Optional format override for this cell
   subroutine add_c(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Character value
      character(*), intent(in) :: val
      !> Optional format override
      character(*), intent(in), optional :: fmt
      character(:), allocatable :: text

      if (present(fmt)) then
         call value_to_string(val, trim(fmt), text)
         call add_from_string(self, text)
      else
         call add_from_string(self, trim(val))
      end if
   end subroutine add_c

   !> Add an integer or real value, re-rendered wide when its edit overflows
   !>
   !> - an edit that writes '*' or text wider than the column is repeated with
   !>   `widen_format`, so decimals and exponent style stay and all digits appear
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    val  Integer or real value
   !> @param[in]    fmt  Format for this cell, without outer parentheses
   subroutine add_number(self, val, fmt)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Integer or real value
      class(*), intent(in) :: val
      !> Format for this cell
      character(*), intent(in) :: fmt
      character(:), allocatable :: s, wide_fmt

      ! The column width is read below, so the overrun check comes first
      call ensure_can_add(self)
      call value_to_string(val, fmt, s)
      if (index(s, "*") > 0 .or. len_trim(adjustl(s)) > self%widths(self%next_col)) then
         call widen_format(fmt, wide_fmt)
         call value_to_string(val, wide_fmt, s)
      end if
      call add_from_string(self, s)
   end subroutine add_number

   !> Append pre-formatted string content as the next cell of the current row
   !>
   !> @param[inout] self Pretty list printer instance
   !> @param[in]    s    Pre-formatted cell text
   subroutine add_from_string(self, s)
      !> Pretty list printer instance
      class(prettylistprinter), intent(inout) :: self
      !> Cell text
      character(*), intent(in) :: s
      character(:), allocatable :: cell

      call ensure_can_add(self)
      if (self%next_col > 1) self%row = self%row//repeat(" ", self%column_gap)
      call format_cell(s, self%widths(self%next_col), cell)
      self%row = self%row//cell
      self%next_col = self%next_col + 1
   end subroutine add_from_string

   !> Stop unless the current row still has space for one more value
   !>
   !> - callers may write into `next_col` unconditionally after this returns
   !>
   !> @param[inout] self Pretty list printer instance
   subroutine ensure_can_add(self)
      !> Pretty list printer instance
      class(prettylistprinter), intent(in) :: self

      ! TODO: add errorprop
      if (self%next_col > self%ncols) then
         error stop "prettylistprinter: too many values in row"
      end if
   end subroutine ensure_can_add

   !> Right-align a cell value in its column, never cutting it
   !>
   !> - leading and trailing blanks of `s` are padding, not content
   !> - content wider than the column is returned in full
   !>
   !> @param[in] s     Source text
   !> @param[in] width Cell width
   !> @param[out] out Right-aligned output cell
   subroutine format_cell(s, width, out)
      !> Source text
      character(*), intent(in) :: s
      !> Cell width
      integer, intent(in) :: width
      !> Right-aligned output cell
      character(:), allocatable, intent(out) :: out
      character(:), allocatable :: text

      text = trim(adjustl(s))
      out = repeat(" ", max(0, width - len(text)))//text
   end subroutine format_cell

   !> Raise the field width of the first numeric edit descriptor to `wide_width`
   !>
   !> - matches I, B, O, Z, F, D, G, E, EN, ES or EX followed by a width, outside
   !>   quoted literals; precision and exponent digits are kept
   !> - a zero width is already minimal and stays; without a match `fmt` is returned
   !> - literal text before the descriptor keeps the widened field's padding
   !>
   !> @param[in] fmt Fortran format string without outer parentheses
   !> @param[out] wide Widened format
   subroutine widen_format(fmt, wide)
      character(*), intent(in) :: fmt
      !> Widened format
      character(:), allocatable, intent(out) :: wide
      character(16) :: wbuf
      character(1) :: quote
      integer :: i, istart, iend

      wide = fmt
      quote = " "
      do i = 1, len(fmt)
         if (quote /= " ") then
            if (fmt(i:i) == quote) quote = " "
         else if (fmt(i:i) == "'" .or. fmt(i:i) == '"') then
            quote = fmt(i:i)
         else if (index("IBOZFDGEibozfdge", fmt(i:i)) > 0) then
            istart = i + 1
            if (index("Ee", fmt(i:i)) > 0 .and. istart <= len(fmt)) then
               if (index("SNXsnx", fmt(istart:istart)) > 0) istart = istart + 1
            end if
            iend = istart - 1
            do while (iend < len(fmt))
               if (index("0123456789", fmt(iend + 1:iend + 1)) == 0) exit
               iend = iend + 1
            end do
            if (iend >= istart) then
               if (verify(fmt(istart:iend), "0") == 0) return
               write (wbuf, "(I0)") wide_width
               wide = fmt(:istart - 1)//trim(wbuf)//fmt(iend + 1:)
               return
            end if
         end if
      end do
   end subroutine widen_format

   !> Compute total printable table width, including inter-column spaces
   !>
   !> @param[in] self Pretty list printer instance
   function table_width(self) result(w)
      class(prettylistprinter), intent(in) :: self
      integer :: w

      if (self%ncols <= 0) then
         w = 0
      else
         w = sum(self%widths) + self%column_gap*(self%ncols - 1)
      end if
   end function table_width

end module moist_utils_prettylistprint
