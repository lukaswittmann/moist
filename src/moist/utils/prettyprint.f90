module moist_utils_prettyprint
   use, intrinsic :: iso_fortran_env, only: output_unit, real64
   use moist_utils_format, only: value_to_string, fixed_fmt, exp_fmt, int_fmt, logical_fmt, value_fmt
   implicit none(type, external)
   private
   public :: prettyprinter, new_prettyprinter

   type :: prettyprinter
      integer :: iu = output_unit
      !> Whether output methods emit lines (false turns every write into a no-op)
      logical :: active = .true.
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
      !> Print a finished work step with its wall time
      procedure :: step
   end type prettyprinter

contains

   !> Set output unit used by the pretty printer
   !>
   !> @param[inout] self Pretty printer instance
   !> @param[in]    iu   Fortran unit number for output
   subroutine set_unit(self, iu)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Fortran output unit number
      integer, intent(in) :: iu
      self%iu = iu
   end subroutine set_unit

   !> Configure layout and numeric formatting options
   !>
   !> @param[inout] self       Pretty printer instance
   !> @param[in]    col_value  Column where first value starts
   !> @param[in]    col_value2 Column where second value starts in `kv2`
   !> @param[in]    indent_step Indentation increment for `push`/`pop`
   !> @param[in]    dot_gap    Spaces between label and dot leader
   !> @param[in]    dot_right  Minimum right-most dots in automatic mode
   !> @param[in]    dot_total  Exact total dots (`<0` enables automatic mode)
   !> @param[in]    fmt_len    Base width to derive default numeric formats
   !> @param[in]    fmt_int    Override integer format
   !> @param[in]    fmt_real   Override fixed real format
   !> @param[in]    fmt_exp    Override exponential real format
   !> @param[in]    fmt_logical Override logical format
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

      if (.not. allocated(self%fmt_int)) call int_fmt(self%fmt_len, self%fmt_int)
      if (.not. allocated(self%fmt_real)) call fixed_fmt(self%fmt_len, 6, self%fmt_real)
      if (.not. allocated(self%fmt_exp)) call exp_fmt(self%fmt_len, 2, self%fmt_exp)
      if (.not. allocated(self%fmt_logical)) call logical_fmt(self%fmt_len, self%fmt_logical)

      if (present(fmt_len)) then
         if (.not. present(fmt_int)) call int_fmt(self%fmt_len, self%fmt_int)
         if (.not. present(fmt_real)) call fixed_fmt(self%fmt_len, 6, self%fmt_real)
         if (.not. present(fmt_exp)) call exp_fmt(self%fmt_len, 2, self%fmt_exp)
         if (.not. present(fmt_logical)) call logical_fmt(self%fmt_len, self%fmt_logical)
      end if

      if (present(fmt_int)) self%fmt_int = trim(fmt_int)
      if (present(fmt_real)) self%fmt_real = trim(fmt_real)
      if (present(fmt_exp)) self%fmt_exp = trim(fmt_exp)
      if (present(fmt_logical)) self%fmt_logical = trim(fmt_logical)
   end subroutine set_layout

   !> Construct a pretty printer with optional layout and format overrides
   !>
   !> @param[in] unit        Optional Fortran output unit
   !> @param[in] col_value   Optional column where first value starts
   !> @param[in] col_value2  Optional column where second value starts in `kv2`
   !> @param[in] indent_step Optional indentation increment for `push`/`pop`
   !> @param[in] dot_gap     Optional spaces between label and dot leader
   !> @param[in] dot_right   Optional minimum right-most dots in automatic mode
   !> @param[in] dot_total   Optional exact total dots (`<0` enables automatic mode)
   !> @param[in] fmt_len     Optional base width to derive default numeric formats
   !> @param[in] fmt_int     Optional integer format override
   !> @param[in] fmt_real    Optional fixed real format override
   !> @param[in] fmt_exp     Optional exponential real format override
   !> @param[in] fmt_logical Optional logical format override
   !> @param[in] active      Optional output gate (default `.true.`; `.false.` silences all writes)
   function new_prettyprinter(unit, col_value, col_value2, indent_step, dot_gap, dot_right, &
      & dot_total, fmt_len, fmt_int, fmt_real, fmt_exp, fmt_logical, active) result(pp)
      !> Optional output unit
      integer, intent(in), optional :: unit
      !> Optional layout scalars
      integer, intent(in), optional :: col_value, col_value2, indent_step, dot_gap, dot_right, dot_total, fmt_len
      !> Optional explicit format strings
      character(*), intent(in), optional :: fmt_int, fmt_real, fmt_exp, fmt_logical
      !> Optional output gate
      logical, intent(in), optional :: active
      !> Constructed pretty printer
      type(prettyprinter) :: pp

      if (present(unit)) pp%iu = unit
      if (present(active)) pp%active = active
      if (.not. pp%active) return

      call pp%set_layout(col_value=col_value, col_value2=col_value2, indent_step=indent_step, &
         & dot_gap=dot_gap, dot_right=dot_right, dot_total=dot_total, fmt_len=fmt_len, &
         & fmt_int=fmt_int, fmt_real=fmt_real, fmt_exp=fmt_exp, fmt_logical=fmt_logical)
   end function new_prettyprinter

   !> Print a blank line, then flush the unit so the finished block is written out
   !>
   !> @param[inout] self Pretty printer instance
   subroutine blank(self)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      if (.not. self%active) return
      write (self%iu, "(A)") ""
      flush (self%iu)
   end subroutine blank

   !> Print a section title at current indentation
   !>
   !> @param[inout] self  Pretty printer instance
   !> @param[in]    title Section title text
   subroutine section(self, title)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Section title
      character(*), intent(in) :: title
      if (.not. self%active) return
      write (self%iu, "(A)") repeat(" ", self%indent)//trim(title)
   end subroutine section

   !> Print a section title and increase indentation level
   !>
   !> @param[inout] self  Pretty printer instance
   !> @param[in]    title Section title text
   subroutine push(self, title)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      !> Section title
      character(*), intent(in) :: title
      if (.not. self%active) return
      call self%section(title)
      self%indent = self%indent + self%indent_step
   end subroutine push

   !> Decrease indentation level by `indent_step`
   !>
   !> @param[inout] self Pretty printer instance
   subroutine pop(self)
      !> Pretty printer instance
      class(prettyprinter), intent(inout) :: self
      if (.not. self%active) return
      self%indent = max(0, self%indent - self%indent_step)
   end subroutine pop

   !> Print one key-value line with optional unit and format override
   !>
   !> @param[inout] self    Pretty printer instance
   !> @param[in]    desc    Label shown on the left side
   !> @param[in]    val     Value to print
   !> @param[in]    unit    Optional unit text after the value
   !> @param[in]    fmt     Optional explicit format string for `val`
   !> @param[in]    use_exp Optional exponential-format selector when `fmt` is absent
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

      if (.not. self%active) return
      call self%set_layout()  ! ensure defaults are allocated

      prefix = repeat(" ", self%indent)
      if (present(fmt)) then
         eff_fmt = trim(fmt)
      else if (present(use_exp)) then
         if (use_exp) then
            eff_fmt = self%fmt_exp
         else
            call value_fmt(val, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt)
         end if
      else
         call value_fmt(val, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt)
      end if
      call value_to_string(val, eff_fmt, vstr)

      left = prefix//trim(desc)//repeat(" ", self%dot_gap)

      nlead = max(0, (self%col_value - 1) - len(left))

      ! Dot leader handling:
      ! - if dot_total >= 0, keep value alignment and place up to dot_total dots
      !   at the right end of the leader field (spaces may appear before dots)
      ! - otherwise, use automatic fill with at least dot_right dots
      if (self%dot_total >= 0) then
         ndots = min(max(0, self%dot_total), nlead)
         nspaces = max(0, nlead - ndots)
         leader = repeat(" ", nspaces)//repeat(".", ndots)
      else
         if (len(left) >= self%col_value - 1 - self%dot_right) then
            ndots = self%dot_right
         else
            ndots = (self%col_value - 1) - len(left)
            ndots = max(self%dot_right, ndots)
         end if
         leader = repeat(".", ndots)
      end if

      if (present(unit)) then
         line = left//leader//" "//vstr//" "//trim(unit)
      else
         line = left//leader//" "//trim(vstr)
      end if

      write (self%iu, "(A)") line
   end subroutine kv

   !> Print one key with two values on the same line
   !>
   !> - useful for SI/AU pairs while keeping aligned second-column output
   !>
   !> @param[inout] self     Pretty printer instance
   !> @param[in]    desc     Label shown on the left side
   !> @param[in]    val1     First value
   !> @param[in]    unit1    Optional unit for first value
   !> @param[in]    val2     Second value
   !> @param[in]    unit2    Optional unit for second value
   !> @param[in]    fmt1     Optional explicit format for first value
   !> @param[in]    fmt2     Optional explicit format for second value
   !> @param[in]    use_exp1 Optional exponential-format selector for first value
   !> @param[in]    use_exp2 Optional exponential-format selector for second value
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

      if (.not. self%active) return
      call self%set_layout()  ! ensure defaults are allocated

      if (present(fmt1)) then
         eff_fmt1 = trim(fmt1)
      else if (present(use_exp1)) then
         if (use_exp1) then
            eff_fmt1 = self%fmt_exp
         else
            call value_fmt(val1, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt1)
         end if
      else
         call value_fmt(val1, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt1)
      end if

      if (present(fmt2)) then
         eff_fmt2 = trim(fmt2)
      else if (present(use_exp2)) then
         if (use_exp2) then
            eff_fmt2 = self%fmt_exp
         else
            call value_fmt(val2, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt2)
         end if
      else
         call value_fmt(val2, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt2)
      end if

      call value_to_string(val1, eff_fmt1, vstr1)
      call value_to_string(val2, eff_fmt2, vstr2)
      u1 = ""
      u2 = ""
      if (present(unit1)) u1 = trim(unit1)
      if (present(unit2)) u2 = trim(unit2)

      prefix = repeat(" ", self%indent)
      left = prefix//trim(desc)//repeat(" ", self%dot_gap)
      nlead = max(0, (self%col_value - 1) - len(left))

      if (self%dot_total >= 0) then
         ndots = min(max(0, self%dot_total), nlead)
         nspaces = max(0, nlead - ndots)
         leader = repeat(" ", nspaces)//repeat(".", ndots)
      else
         if (len(left) >= self%col_value - 1 - self%dot_right) then
            ndots = self%dot_right
         else
            ndots = (self%col_value - 1) - len(left)
            ndots = max(self%dot_right, ndots)
         end if
         leader = repeat(".", ndots)
      end if

      line = left//leader//" "//vstr1
      if (len(u1) > 0) line = line//" "//u1
      nspaces = (self%col_value2 - 1) - len(line)
      if (nspaces > 0) then
         line = line//repeat(" ", nspaces)
      else
         line = line//"  "
      end if
      line = line//vstr2
      if (len(u2) > 0) line = line//" "//u2

      write (self%iu, "(A)") line
   end subroutine kv2

   !> Print one key with three values on the same line and one shared
   !> trailing unit
   !>
   !> - useful for component triples while keeping aligned multi-column output
   !>
   !> @param[inout] self     Pretty printer instance
   !> @param[in]    desc     Label shown on the left side
   !> @param[in]    val1     First value
   !> @param[in]    val2     Second value
   !> @param[in]    val3     Third value
   !> @param[in]    unit     Optional unit printed once after the third value
   !> @param[in]    fmt1     Optional explicit format for first value
   !> @param[in]    fmt2     Optional explicit format for second value
   !> @param[in]    fmt3     Optional explicit format for third value
   !> @param[in]    use_exp1 Optional exponential-format selector for first value
   !> @param[in]    use_exp2 Optional exponential-format selector for second value
   !> @param[in]    use_exp3 Optional exponential-format selector for third value
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

      if (.not. self%active) return
      call self%set_layout()  ! ensure defaults are allocated

      if (present(fmt1)) then
         eff_fmt1 = trim(fmt1)
      else if (present(use_exp1)) then
         if (use_exp1) then
            eff_fmt1 = self%fmt_exp
         else
            call value_fmt(val1, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt1)
         end if
      else
         call value_fmt(val1, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt1)
      end if

      if (present(fmt2)) then
         eff_fmt2 = trim(fmt2)
      else if (present(use_exp2)) then
         if (use_exp2) then
            eff_fmt2 = self%fmt_exp
         else
            call value_fmt(val2, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt2)
         end if
      else
         call value_fmt(val2, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt2)
      end if

      if (present(fmt3)) then
         eff_fmt3 = trim(fmt3)
      else if (present(use_exp3)) then
         if (use_exp3) then
            eff_fmt3 = self%fmt_exp
         else
            call value_fmt(val3, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt3)
         end if
      else
         call value_fmt(val3, self%fmt_int, self%fmt_real, self%fmt_exp, self%fmt_logical, eff_fmt3)
      end if

      call value_to_string(val1, eff_fmt1, vstr1)
      call value_to_string(val2, eff_fmt2, vstr2)
      call value_to_string(val3, eff_fmt3, vstr3)
      u = ""
      if (present(unit)) u = trim(unit)

      prefix = repeat(" ", self%indent)
      left = prefix//trim(desc)//repeat(" ", self%dot_gap)
      nlead = max(0, (self%col_value - 1) - len(left))

      if (self%dot_total >= 0) then
         ndots = min(max(0, self%dot_total), nlead)
         nspaces = max(0, nlead - ndots)
         leader = repeat(" ", nspaces)//repeat(".", ndots)
      else
         if (len(left) >= self%col_value - 1 - self%dot_right) then
            ndots = self%dot_right
         else
            ndots = (self%col_value - 1) - len(left)
            ndots = max(self%dot_right, ndots)
         end if
         leader = repeat(".", ndots)
      end if

      line = left//leader//" "//vstr1
      line = line//repeat(" ", value_gap)//vstr2
      line = line//repeat(" ", value_gap)//vstr3
      if (len(u) > 0) line = line//" "//u

      write (self%iu, "(A)") line
   end subroutine kvvv

   !> Print a finished work step with its wall time: `label ... done (12.34s)`
   !>
   !> @param[inout] self    Pretty printer instance
   !> @param[in]    desc    Step name
   !> @param[in]    seconds Wall time spent in the step (s)
   subroutine step(self, desc, seconds)
      class(prettyprinter), intent(inout) :: self
      character(*), intent(in) :: desc
      real(real64), intent(in) :: seconds

      character(len=32) :: time

      if (.not. self%active) return
      ! Fixed width keeps the leading zero, which F0.d may drop
      write (time, "(f20.2)") seconds
      call self%kv(desc, "done ("//trim(adjustl(time))//"s)")
   end subroutine step

end module moist_utils_prettyprint
