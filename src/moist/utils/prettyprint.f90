module moist_utils_prettyprint
   use, intrinsic :: iso_fortran_env, only: output_unit, int8, int16, int32, int64, &
     & real32, real64
   implicit none
   private
   public :: prettyprinter, new_prettyprinter

   type :: prettyprinter
      integer :: iu = output_unit
      integer :: indent = 0
      integer :: indent_step = 2
      integer :: col_value = 32  ! 1-based column where the value starts
      integer :: col_value2 = 58 ! 1-based column where second value starts (kv2)
      integer :: dot_gap = 1     ! spaces after key before dot leaders
      integer :: dot_right = 1   ! minimum number of dots nearest value
      integer :: dot_total = 3   ! exact total dots; <0 means automatic
      integer :: fmt_len = 16
      character(:), allocatable :: fmt_int
      character(:), allocatable :: fmt_real
      character(:), allocatable :: fmt_exp
      character(:), allocatable :: fmt_logical
   contains
      procedure :: set_unit
      procedure :: set_layout
      !> Print blank like
      procedure :: blank
      !> Print a section header with current indentation
      procedure :: section
      !> Increase indentation and print section header
      procedure :: push
      !> Decrease indentation
      procedure :: pop
      !> Print a key-value pair with optional unit and formatting
      procedure :: kv
      !> Print two key-value pairs on the same line with optional units and formatting
      procedure :: kv2
      !> Print three values on the same line with optional units and formatting
      procedure :: kvvv
   end type prettyprinter

contains

   !> Set output unit used by the pretty printer.
   !> @param[inout] self Pretty printer instance.
   !> @param[in]    iu   Fortran unit number for output.
   subroutine set_unit(self, iu)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Fortran output unit number
      integer, intent(in) :: iu
      self%iu = iu
   end subroutine set_unit

   !> Configure layout and numeric formatting options.
   !> @param[inout] self       Pretty printer instance.
   !> @param[in]    col_value  Column where first value starts.
   !> @param[in]    col_value2 Column where second value starts in `kv2`.
   !> @param[in]    indent_step Indentation increment for `push`/`pop`.
   !> @param[in]    dot_gap    Spaces between label and dot leader.
   !> @param[in]    dot_right  Minimum right-most dots in automatic mode.
   !> @param[in]    dot_total  Exact total dots (`<0` enables automatic mode).
   !> @param[in]    fmt_len    Base width to derive default numeric formats.
   !> @param[in]    fmt_int    Override integer format.
   !> @param[in]    fmt_real   Override fixed real format.
   !> @param[in]    fmt_exp    Override exponential real format.
   !> @param[in]    fmt_logical Override logical format.
   subroutine set_layout(self, col_value, col_value2, indent_step, dot_gap, dot_right, dot_total, &
     & fmt_len, fmt_int, fmt_real, fmt_exp, fmt_logical)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Optional layout scalars
      integer, intent(in), optional :: col_value, col_value2, indent_step, dot_gap, dot_right, dot_total, fmt_len
      !> Optional explicit format strings
      character(*), intent(in), optional :: fmt_int, fmt_real, fmt_exp, fmt_logical

      if (present(col_value)) self%col_value = col_value
      if (present(col_value2)) self%col_value2 = col_value2
      if (present(indent_step)) self%indent_step = indent_step
      if (present(dot_gap)) self%dot_gap = max(0, dot_gap)
      if (present(dot_right)) self%dot_right = max(0, dot_right)
      if (present(dot_total)) self%dot_total = dot_total
      if (present(fmt_len)) self%fmt_len = fmt_len

      if (.not. allocated(self%fmt_int)) self%fmt_int = int_fmt(self%fmt_len)
      if (.not. allocated(self%fmt_real)) self%fmt_real = fixed_fmt(self%fmt_len, 6)
      if (.not. allocated(self%fmt_exp)) self%fmt_exp = exp_fmt(self%fmt_len, 2)
      if (.not. allocated(self%fmt_logical)) self%fmt_logical = 'L1'

      if (present(fmt_len)) then
         if (.not. present(fmt_int)) self%fmt_int = int_fmt(self%fmt_len)
         if (.not. present(fmt_real)) self%fmt_real = fixed_fmt(self%fmt_len, 6)
         if (.not. present(fmt_exp)) self%fmt_exp = exp_fmt(self%fmt_len, 2)
      end if

      if (present(fmt_int)) self%fmt_int = trim(fmt_int)
      if (present(fmt_real)) self%fmt_real = trim(fmt_real)
      if (present(fmt_exp)) self%fmt_exp = trim(fmt_exp)
      if (present(fmt_logical)) self%fmt_logical = trim(fmt_logical)
   end subroutine set_layout

   !> Construct a pretty printer with optional layout and format overrides.
   !> @param[in] unit        Optional Fortran output unit.
   !> @param[in] col_value   Optional column where first value starts.
   !> @param[in] col_value2  Optional column where second value starts in `kv2`.
   !> @param[in] indent_step Optional indentation increment for `push`/`pop`.
   !> @param[in] dot_gap     Optional spaces between label and dot leader.
   !> @param[in] dot_right   Optional minimum right-most dots in automatic mode.
   !> @param[in] dot_total   Optional exact total dots (`<0` enables automatic mode).
   !> @param[in] fmt_len     Optional base width to derive default numeric formats.
   !> @param[in] fmt_int     Optional integer format override.
   !> @param[in] fmt_real    Optional fixed real format override.
   !> @param[in] fmt_exp     Optional exponential real format override.
   !> @param[in] fmt_logical Optional logical format override.
   function new_prettyprinter(unit, col_value, col_value2, indent_step, dot_gap, dot_right, &
      & dot_total, fmt_len, fmt_int, fmt_real, fmt_exp, fmt_logical) result(pp)
      !> Optional output unit
      integer, intent(in), optional :: unit
      !> Optional layout scalars
      integer, intent(in), optional :: col_value, col_value2, indent_step, dot_gap, dot_right, dot_total, fmt_len
      !> Optional explicit format strings
      character(*), intent(in), optional :: fmt_int, fmt_real, fmt_exp, fmt_logical
      !> Constructed pretty printer
      type(prettyprinter) :: pp

      if (present(unit)) pp%iu = unit

      call pp%set_layout(col_value=col_value, col_value2=col_value2, indent_step=indent_step, &
         & dot_gap=dot_gap, dot_right=dot_right, dot_total=dot_total, fmt_len=fmt_len, &
         & fmt_int=fmt_int, fmt_real=fmt_real, fmt_exp=fmt_exp, fmt_logical=fmt_logical)
   end function new_prettyprinter

   !> Print a blank line.
   !> @param[inout] self Pretty printer instance.
   subroutine blank(self)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      write (self%iu, '(A)') ''
   end subroutine blank

   !> Print a section title at current indentation.
   !> @param[inout] self  Pretty printer instance.
   !> @param[in]    title Section title text.
   subroutine section(self, title)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Section title
      character(*), intent(in) :: title
      write (self%iu, '(A)') repeat(' ', self%indent)//trim(title)
   end subroutine section

   !> Print a section title and increase indentation level.
   !> @param[inout] self  Pretty printer instance.
   !> @param[in]    title Section title text.
   subroutine push(self, title)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Section title
      character(*), intent(in) :: title
      call self%section(title)
      self%indent = self%indent + self%indent_step
   end subroutine push

   !> Decrease indentation level by `indent_step`.
   !> @param[inout] self Pretty printer instance.
   subroutine pop(self)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      self%indent = max(0, self%indent - self%indent_step)
   end subroutine pop

   !> Print one key-value line with optional unit and format override.
   !> @param[inout] self    Pretty printer instance.
   !> @param[in]    desc    Label shown on the left side.
   !> @param[in]    val     Value to print.
   !> @param[in]    unit    Optional unit text after the value.
   !> @param[in]    fmt     Optional explicit format string for `val`.
   !> @param[in]    use_exp Optional exponential-format selector when `fmt` is absent.
   subroutine kv(self, desc, val, unit, fmt, use_exp)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Label text
      character(*), intent(in) :: desc
      !> Value to print
      class(*), intent(in) :: val
      !> Optional unit text
      character(*), intent(in), optional :: unit
      !> Optional explicit format override
      character(*), intent(in), optional :: fmt
      !> Optional exponential formatting switch
      logical, intent(in), optional :: use_exp

      character(:), allocatable :: prefix, left, leader, vstr, line, eff_fmt
      integer :: ndots, nlead, nspaces

      call self%set_layout()  ! ensure defaults are allocated

      prefix = repeat(' ', self%indent)
      if (present(fmt)) then
         eff_fmt = trim(fmt)
      else if (present(use_exp)) then
         if (use_exp) then
            eff_fmt = self%fmt_exp
         else
            eff_fmt = default_fmt(self, val)
         end if
      else
         eff_fmt = default_fmt(self, val)
      end if
      vstr = value_to_string(val, eff_fmt)

      left = prefix//trim(desc)//repeat(' ', self%dot_gap)

      nlead = max(0, (self%col_value - 1) - len(left))

      ! Dot leader handling:
      ! - if dot_total >= 0, keep value alignment and place up to dot_total dots
      !   at the right end of the leader field (spaces may appear before dots)
      ! - otherwise, use automatic fill with at least dot_right dots
      if (self%dot_total >= 0) then
         ndots = min(max(0, self%dot_total), nlead)
         nspaces = max(0, nlead - ndots)
         leader = repeat(' ', nspaces)//repeat('.', ndots)
      else
         if (len(left) >= self%col_value - 1 - self%dot_right) then
            ndots = self%dot_right
         else
            ndots = (self%col_value - 1) - len(left)
            ndots = max(self%dot_right, ndots)
         end if
         leader = repeat('.', ndots)
      end if

      if (present(unit)) then
         line = left//leader//' '//vstr//' '//trim(unit)
      else
         line = left//leader//' '//trim(vstr)
      end if

      write (self%iu, '(A)') line
   end subroutine kv

   !> Print one key with two values on the same line.
   !> Useful for SI/AU pairs while keeping aligned second-column output.
   !> @param[inout] self     Pretty printer instance.
   !> @param[in]    desc     Label shown on the left side.
   !> @param[in]    val1     First value.
   !> @param[in]    unit1    Optional unit for first value.
   !> @param[in]    val2     Second value.
   !> @param[in]    unit2    Optional unit for second value.
   !> @param[in]    fmt1     Optional explicit format for first value.
   !> @param[in]    fmt2     Optional explicit format for second value.
   !> @param[in]    use_exp1 Optional exponential-format selector for first value.
   !> @param[in]    use_exp2 Optional exponential-format selector for second value.
   subroutine kv2(self, desc, val1, unit1, val2, unit2, fmt1, fmt2, use_exp1, use_exp2)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Label text
      character(*), intent(in) :: desc
      !> First and second values
      class(*), intent(in) :: val1, val2
      !> Optional units
      character(*), intent(in), optional :: unit1, unit2
      !> Optional explicit format overrides
      character(*), intent(in), optional :: fmt1, fmt2
      !> Optional exponential formatting switches
      logical, intent(in), optional :: use_exp1, use_exp2

      character(:), allocatable :: prefix, left, leader, line
      character(:), allocatable :: eff_fmt1, eff_fmt2
      character(:), allocatable :: vstr1, vstr2
      character(:), allocatable :: u1, u2
      integer :: ndots, nlead, nspaces

      call self%set_layout()  ! ensure defaults are allocated

      if (present(fmt1)) then
         eff_fmt1 = trim(fmt1)
      else if (present(use_exp1)) then
         if (use_exp1) then
            eff_fmt1 = self%fmt_exp
         else
            eff_fmt1 = default_fmt(self, val1)
         end if
      else
         eff_fmt1 = default_fmt(self, val1)
      end if

      if (present(fmt2)) then
         eff_fmt2 = trim(fmt2)
      else if (present(use_exp2)) then
         if (use_exp2) then
            eff_fmt2 = self%fmt_exp
         else
            eff_fmt2 = default_fmt(self, val2)
         end if
      else
         eff_fmt2 = default_fmt(self, val2)
      end if

      vstr1 = value_to_string(val1, eff_fmt1)
      vstr2 = value_to_string(val2, eff_fmt2)
      u1 = ''
      u2 = ''
      if (present(unit1)) u1 = trim(unit1)
      if (present(unit2)) u2 = trim(unit2)

      prefix = repeat(' ', self%indent)
      left = prefix//trim(desc)//repeat(' ', self%dot_gap)
      nlead = max(0, (self%col_value - 1) - len(left))

      if (self%dot_total >= 0) then
         ndots = min(max(0, self%dot_total), nlead)
         nspaces = max(0, nlead - ndots)
         leader = repeat(' ', nspaces)//repeat('.', ndots)
      else
         if (len(left) >= self%col_value - 1 - self%dot_right) then
            ndots = self%dot_right
         else
            ndots = (self%col_value - 1) - len(left)
            ndots = max(self%dot_right, ndots)
         end if
         leader = repeat('.', ndots)
      end if

      line = left//leader//' '//vstr1
      if (len(u1) > 0) line = line//' '//u1
      nspaces = (self%col_value2 - 1) - len(line)
      if (nspaces > 0) then
         line = line//repeat(' ', nspaces)
      else
         line = line//'  '
      end if
      line = line//vstr2
      if (len(u2) > 0) line = line//' '//u2

      write (self%iu, '(A)') line
   end subroutine kv2

   !> Print one key with three values on the same line and one shared trailing unit.
   !> Useful for component triples while keeping aligned multi-column output.
   !> @param[inout] self     Pretty printer instance.
   !> @param[in]    desc     Label shown on the left side.
   !> @param[in]    val1     First value.
   !> @param[in]    val2     Second value.
   !> @param[in]    val3     Third value.
   !> @param[in]    unit     Optional unit printed once after the third value.
   !> @param[in]    fmt1     Optional explicit format for first value.
   !> @param[in]    fmt2     Optional explicit format for second value.
   !> @param[in]    fmt3     Optional explicit format for third value.
   !> @param[in]    use_exp1 Optional exponential-format selector for first value.
   !> @param[in]    use_exp2 Optional exponential-format selector for second value.
   !> @param[in]    use_exp3 Optional exponential-format selector for third value.
   subroutine kvvv(self, desc, val1, val2, val3, unit, fmt1, fmt2, fmt3, use_exp1, use_exp2, &
     & use_exp3)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Label text
      character(*), intent(in) :: desc
      !> First, second, and third values
      class(*), intent(in) :: val1, val2, val3
      !> Optional trailing unit
      character(*), intent(in), optional :: unit
      !> Optional explicit format overrides
      character(*), intent(in), optional :: fmt1, fmt2, fmt3
      !> Optional exponential formatting switches
      logical, intent(in), optional :: use_exp1, use_exp2, use_exp3

      character(:), allocatable :: prefix, left, leader, line
      character(:), allocatable :: eff_fmt1, eff_fmt2, eff_fmt3
      character(:), allocatable :: vstr1, vstr2, vstr3
      character(:), allocatable :: u
      integer, parameter :: value_gap = 2
      integer :: ndots, nlead, nspaces

      call self%set_layout()  ! ensure defaults are allocated

      if (present(fmt1)) then
         eff_fmt1 = trim(fmt1)
      else if (present(use_exp1)) then
         if (use_exp1) then
            eff_fmt1 = self%fmt_exp
         else
            eff_fmt1 = default_fmt(self, val1)
         end if
      else
         eff_fmt1 = default_fmt(self, val1)
      end if

      if (present(fmt2)) then
         eff_fmt2 = trim(fmt2)
      else if (present(use_exp2)) then
         if (use_exp2) then
            eff_fmt2 = self%fmt_exp
         else
            eff_fmt2 = default_fmt(self, val2)
         end if
      else
         eff_fmt2 = default_fmt(self, val2)
      end if

      if (present(fmt3)) then
         eff_fmt3 = trim(fmt3)
      else if (present(use_exp3)) then
         if (use_exp3) then
            eff_fmt3 = self%fmt_exp
         else
            eff_fmt3 = default_fmt(self, val3)
         end if
      else
         eff_fmt3 = default_fmt(self, val3)
      end if

      vstr1 = value_to_string(val1, eff_fmt1)
      vstr2 = value_to_string(val2, eff_fmt2)
      vstr3 = value_to_string(val3, eff_fmt3)
      u = ''
      if (present(unit)) u = trim(unit)

      prefix = repeat(' ', self%indent)
      left = prefix//trim(desc)//repeat(' ', self%dot_gap)
      nlead = max(0, (self%col_value - 1) - len(left))

      if (self%dot_total >= 0) then
         ndots = min(max(0, self%dot_total), nlead)
         nspaces = max(0, nlead - ndots)
         leader = repeat(' ', nspaces)//repeat('.', ndots)
      else
         if (len(left) >= self%col_value - 1 - self%dot_right) then
            ndots = self%dot_right
         else
            ndots = (self%col_value - 1) - len(left)
            ndots = max(self%dot_right, ndots)
         end if
         leader = repeat('.', ndots)
      end if

      line = left//leader//' '//vstr1
      line = line//repeat(' ', value_gap)//vstr2
      line = line//repeat(' ', value_gap)//vstr3
      if (len(u) > 0) line = line//' '//u

      write (self%iu, '(A)') line
   end subroutine kvvv

   function default_fmt(self, val) result(fmt)
      class(prettyprinter), intent(in) :: self
      class(*), intent(in) :: val
      character(:), allocatable :: fmt

      select type (val)
      type is (integer(int8))
         fmt = self%fmt_int
      type is (integer(int16))
         fmt = self%fmt_int
      type is (integer(int32))
         fmt = self%fmt_int
      type is (integer(int64))
         fmt = self%fmt_int
      type is (real(real32))
         fmt = default_real_fmt(self, real(val, kind=real64))
      type is (real(real64))
         fmt = default_real_fmt(self, val)
      type is (logical)
         fmt = self%fmt_logical
      type is (character(*))
         fmt = 'A'
      class default
         fmt = 'A'
      end select
   end function default_fmt

   function value_to_string(val, fmt) result(s)
      class(*), intent(in) :: val
      character(*), intent(in) :: fmt
      character(:), allocatable :: s
      character(256) :: buf
      character(:), allocatable :: f

      f = '('//trim(fmt)//')'
      buf = ''

      select type (val)
      type is (integer(int8))
         write (buf, f) val
      type is (integer(int16))
         write (buf, f) val
      type is (integer(int32))
         write (buf, f) val
      type is (integer(int64))
         write (buf, f) val
      type is (real(real32))
         if (val == 0.0_real32) then
            s = zero_value_string(fmt)
            return
         end if
         write (buf, f) val
      type is (real(real64))
         if (val == 0.0_real64) then
            s = zero_value_string(fmt)
            return
         end if
         write (buf, f) val
      type is (logical)
         write (buf, f) val
      type is (character(*))
         buf = val
      class default
         buf = '<unsupported type>'
      end select

      s = trim(buf)
   end function value_to_string

   function zero_value_string(fmt) result(s)
      character(*), intent(in) :: fmt
      character(:), allocatable :: s
      character(256) :: buf
      character(:), allocatable :: f
      integer :: idot, w

      f = '('//trim(fmt)//')'
      write (buf, f) 0.0_real64
      idot = index(buf, '.')
      w = len_trim(buf)

      if (idot > 1 .and. w > 0) then
         s = repeat(' ', idot - 2)//'0.0'//repeat(' ', max(0, w - (idot + 1)))
      else
         s = '0.0'
      end if
   end function zero_value_string

   function fixed_fmt(width, decimals) result(fmt)
      integer, intent(in) :: width, decimals
      character(:), allocatable :: fmt
      character(32) :: wbuf, dbuf

      write (wbuf, '(I0)') max(1, width)
      write (dbuf, '(I0)') max(0, decimals)
      fmt = 'F'//trim(wbuf)//'.'//trim(dbuf)
   end function fixed_fmt

   function exp_fmt(width, decimals) result(fmt)
      integer, intent(in) :: width, decimals
      character(:), allocatable :: fmt
      character(32) :: wbuf, dbuf

      write (wbuf, '(I0)') max(1, width)
      write (dbuf, '(I0)') max(0, decimals)
      fmt = 'ES'//trim(wbuf)//'.'//trim(dbuf)
   end function exp_fmt

   function int_fmt(width) result(fmt)
      integer, intent(in) :: width
      character(:), allocatable :: fmt
      character(32) :: wbuf

      write (wbuf, '(I0)') max(1, width - 7)
      fmt = 'I'//trim(wbuf)
   end function int_fmt

   function default_real_fmt(self, val) result(fmt)
      class(prettyprinter), intent(in) :: self
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

end module moist_utils_prettyprint
