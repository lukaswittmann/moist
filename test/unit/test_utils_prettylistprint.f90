!> Test suite for the tabular printer in moist_utils_prettylistprint
module test_utils_prettylistprint
   use, intrinsic :: iso_fortran_env, only: real32, real64, int8, int16, int32, int64, iostat_eor, iostat_end
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed, to_string
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none(type, external)
   private

   public :: collect_utils_prettylistprint

contains

   !> Collect all tabular-printer tests
   !>
   !> Misuse of the printer (mismatched widths and headers, no columns, a row
   !> that over- or underruns its column count) is a programming error and now
   !> ends in `error stop`, so it cannot be exercised from inside the test
   !> binary. What remains testable is the output the printer produces
   subroutine collect_utils_prettylistprint(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("healthy_printer_emits", test_healthy_printer_emits), &
                  new_unittest("row_layout", test_row_layout), &
                  new_unittest("decorations", test_decorations), &
                  new_unittest("real_overflow_in_full", test_real_overflow_in_full), &
                  new_unittest("scalar_kinds", test_scalar_kinds), &
                  new_unittest("real_formats", test_real_formats), &
                  new_unittest("state_and_spacing", test_state_and_spacing), &
                  new_unittest("decoration_content", test_decoration_content), &
                  new_unittest("constructor_formats", test_constructor_formats), &
                  new_unittest("narrow_numeric_cells", test_narrow_numeric_cells) &
                  ]

   end subroutine collect_utils_prettylistprint

   !> A well-formed printer produces a header and one row
   !>
   !> @param[out] error Error handle
   subroutine test_healthy_printer_emits(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_healthy.tmp"
      integer :: iu, nlines

      call open_scratch(path, iu)

      plp = new_prettylistprinter([8, 8], ["alpha", "beta "], unit=iu)
      call plp%print_header()
      call plp%begin_row()
      call plp%add(1)
      call plp%add(2.5_real64, fmt="f6.2")
      call plp%end_row()

      call close_scratch(iu)
      call count_lines(path, nlines)
      call check(error, nlines, 2, "header and one row were written")

   end subroutine test_healthy_printer_emits

   !> Values land right-aligned in their fields, and `skip` leaves a cell blank
   subroutine test_row_layout(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_layout.tmp"
      integer :: iu

      call open_scratch(path, iu)

      plp = new_prettylistprinter([6, 6, 6], ["a", "b", "c"], unit=iu, offset=0, column_gap=0)
      call plp%begin_row()
      call plp%add("ab")
      call plp%skip()
      call plp%add(42, fmt="I3")
      call plp%end_row()

      call close_scratch(iu)
      call check_line(error, path, 1, "____ab__________42", &
                      "cells are right-aligned and the skipped one is blank")
      call discard_scratch(path)

   end subroutine test_row_layout

   !> Section header, separator and blank line each emit exactly one line
   subroutine test_decorations(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_decor.tmp"
      character(len=:), allocatable :: line
      integer :: iu, nlines

      call open_scratch(path, iu)

      plp = new_prettylistprinter([5, 5], ["a", "b"], unit=iu, offset=0, column_gap=0)
      call plp%header("67")
      call plp%separator()
      call plp%blank()

      call close_scratch(iu)
      call check_line(error, path, 1, "===_67_===", "the title is centered in '=' fill")
      if (allocated(error)) then
         call discard_scratch(path)
         return
      end if

      call read_line(path, 2, line, error)
      if (.not. allocated(error)) then
         call check(error, len_trim(line), 10, "the separator spans the table width")
      end if
      if (allocated(error)) then
         call discard_scratch(path)
         return
      end if

      call count_lines(path, nlines)
      call check(error, nlines, 3, "header, separator and blank were written")

   end subroutine test_decorations

   !> A real whose edit overflows is re-rendered in full with the same descriptor
   !>
   !> F6.2 and ES8.4 fill their fields with '*' here; the cell keeps the sign, the
   !> decimals and the exponent style, and later cells shift right
   subroutine test_real_overflow_in_full(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_overflow.tmp"
      integer :: iu

      call open_scratch(path, iu)

      plp = new_prettylistprinter([5, 5, 5], ["a", "b", "c"], unit=iu, offset=0, column_gap=1)
      call plp%begin_row()
      call plp%add(1.0e12_real64, fmt="f6.2")
      call plp%add(-1.0e12_real64, fmt="f6.2")
      call plp%add(1.5e6_real64, fmt="es8.4")
      call plp%end_row()

      call close_scratch(iu)
      call check_line(error, path, 1, "1000000000000.00_-1000000000000.00_1.5000E+06", &
                      "overflowing F and ES edits print every digit with sign and style")
      call discard_scratch(path)

   end subroutine test_real_overflow_in_full

   !> Integer kinds and per-cell overrides retain distinct scalar values
   !>
   !> @param[out] error Error handle
   subroutine test_scalar_kinds(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      type(prettylistprinter) :: plp
      integer :: iu
      character(len=*), parameter :: path = "plp_kinds.tmp"

      call open_scratch(path, iu)
      plp = new_prettylistprinter([16, 16, 16, 16, 2, 2], ["a", "b", "c", "d", "e", "f"], &
                                  unit=iu, offset=0, column_gap=0)
      call plp%add(-12_int8)
      call plp%add(234_int16)
      call plp%add(-345_int32)
      call plp%add(4567_int64)
      call plp%add(.true.)
      call plp%add(.false.)
      call plp%end_row()
      call plp%add(7_int8, fmt="I3.3")
      call plp%add(8_int16, fmt="I3.3")
      call plp%add(9_int32, fmt="I3.3")
      call plp%add(10_int64, fmt="I3.3")
      call plp%add("abcdef")
      call plp%add("z")
      call plp%end_row()
      call close_scratch(iu)
      call check_line(error, path, 1, repeat("_", 13)//"-12"//repeat("_", 13)//"234"// &
                      repeat("_", 12)//"-345"//repeat("_", 12)//"4567_T_F", &
                      "all integer kinds and logical values")
      if (.not. allocated(error)) then
         call check_line(error, path, 2, repeat("_", 13)//"007"//repeat("_", 13)//"008"// &
                         repeat("_", 13)//"009"//repeat("_", 13)//"010abcdef_z", &
                         "integer overrides, and a string wider than its column in full")
      end if
      call discard_scratch(path)
   end subroutine test_scalar_kinds

   !> Real kinds, magnitude selection, zero and configured formats
   subroutine test_real_formats(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      type(prettylistprinter) :: plp
      integer :: iu
      character(len=*), parameter :: path = "plp_formats.tmp"

      call open_scratch(path, iu)
      plp = new_prettylistprinter([16, 16], ["a", "b"], unit=iu, offset=0, column_gap=0)
      call plp%add(1.25_real32)
      call plp%add(-1.25_real64)
      call plp%end_row()
      call plp%add(1.0e6_real32)
      call plp%add(1.0e-7_real64)
      call plp%end_row()
      call plp%add(0.0_real32)
      call plp%add(-0.0_real64)
      call plp%end_row()
      call plp%add(1.0e-5_real32)
      call plp%add(1.0e5_real64)
      call plp%end_row()
      call plp%set_real_formats(fmt_real="F10.2", fmt_exp="ES12.3")
      call plp%add(-1.25_real32)
      call plp%add(1.0e6_real64)
      call plp%end_row()
      call plp%add(1.24_real32, fmt="F8.1")
      call plp%add(1.24_real64, fmt="F8.1")
      call plp%end_row()
      call close_scratch(iu)
      call check_line(error, path, 1, "__________1.2500_________-1.2500", &
                      "default fixed precision for both real kinds")
      if (.not. allocated(error)) then
         call check_line(error, path, 2, "______1.0000E+06______1.0000E-07", "default exponential precision")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 3, "_____________0.0_____________0.0", "canonical real zeros")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 4, "__________0.0000_____100000.0000", &
                         "interior magnitudes use fixed format")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 5, "___________-1.25_______1.000E+06", "setter updates both real formats")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 6, "_____________1.2_____________1.2", &
                         "cell format overrides both real kinds")
      end if
      call discard_scratch(path)
      if (allocated(error)) return
      call open_scratch(path, iu)
      plp = new_prettylistprinter([5, 5], ["a", "b"], unit=iu, offset=0, column_gap=0)
      call plp%add(1.0e12_real32, fmt="F6.2")
      call plp%add(-1.0e12_real32, fmt="F6.2")
      call plp%end_row()
      call plp%add(12345.0_real32, fmt="F10.2")
      call plp%add(-12345.0_real64, fmt="F10.2")
      call plp%end_row()
      call close_scratch(iu)
      ! 999999995904 is the real32 nearest to 1e12; F6.2 writes '*' for it
      call check_line(error, path, 1, "999999995904.00-999999995904.00", &
                      "real32 overflow prints every integer digit and the sign")
      if (.not. allocated(error)) then
         call check_line(error, path, 2, "12345.00-12345.00", "values wider than their cells print in full")
      end if
      call discard_scratch(path)
   end subroutine test_real_formats

   !> Header routing, row reset and configurable gaps
   !>
   !> @param[out] error Error handle
   subroutine test_state_and_spacing(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      type(prettylistprinter) :: plp
      integer :: iu
      character(len=*), parameter :: path = "plp_state.tmp"

      call open_scratch(path, iu)
      plp = new_prettylistprinter([3, 4], ["abcdef", "b     "], unit=iu, offset=2, column_gap=2)
      call plp%print_header()
      ! Pending values dropped by begin_row must not reach the next printed row
      call plp%add("x")
      call plp%add("y")
      call plp%begin_row()
      call plp%add("a")
      call plp%skip()
      call plp%end_row()
      call plp%set_column_gap(-3)
      call plp%add("c")
      call plp%add("d")
      call plp%end_row()
      call close_scratch(iu)
      call check_line(error, path, 1, "__abcdef_____b", "a header wider than its column prints in full")
      if (.not. allocated(error)) then
         call check_line(error, path, 2, "____a______", &
                         "begin_row drops pending values; a skipped last cell keeps its blank width")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 3, "____c___d", "end_row resets cursor and negative gap clamps to zero")
      end if
      call discard_scratch(path)
      if (allocated(error)) return
      call open_scratch(path, iu)
      plp = new_prettylistprinter([0, -1], ["a", "b"], unit=iu, offset=-2, column_gap=-1)
      call plp%add("x")
      call plp%add("y")
      call plp%end_row()
      call close_scratch(iu)
      call check_line(error, path, 1, "xy", "nonpositive widths clamp to one and spacing to zero")
      call discard_scratch(path)
   end subroutine test_state_and_spacing

   !> Decorations retain fill, gaps, clipping and truly blank lines
   subroutine test_decoration_content(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      type(prettylistprinter) :: plp
      integer :: iu
      character(len=:), allocatable :: line
      character(len=*), parameter :: path = "plp_content.tmp"

      call open_scratch(path, iu)
      plp = new_prettylistprinter([3, 4], ["a", "b"], unit=iu, offset=2, column_gap=2)
      call plp%separator()
      call plp%header("  six seven  ")
      call plp%header("dsfgdfgs")
      call plp%blank()
      call close_scratch(iu)
      call check_line(error, path, 1, "__---__----", "separator uses dashes with offset and column gap")
      if (.not. allocated(error)) then
         call check_line(error, path, 2, "___six_seve", "long title clips from the right at the table width")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 3, "___dsfgdfgs", "section width includes only gaps between columns")
      end if
      if (.not. allocated(error)) then
         call read_line(path, 4, line, error)
         if (.not. allocated(error)) call check(error, len(line), 0, "blank emits an empty record")
      end if
      call discard_scratch(path)
   end subroutine test_decoration_content

   !> Constructor formats and logical overrides are honored independently
   !>
   !> @param[out] error Error handle
   subroutine test_constructor_formats(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      type(prettylistprinter) :: plp
      integer :: iu
      character(len=*), parameter :: path = "plp_constructor.tmp"

      call open_scratch(path, iu)
      plp = new_prettylistprinter([16, 16, 16, 16], ["a", "b", "c", "d"], &
                                  unit=iu, offset=0, column_gap=0, fmt_int="I3.3", &
                                  fmt_real="F10.1", fmt_exp="ES12.1", fmt_logical="'L=',L1")
      call plp%add(7)
      call plp%add(1.24_real64)
      call plp%add(1.0e6_real32)
      call plp%add(.true.)
      call plp%end_row()
      call plp%add(8)
      call plp%add(1.24_real32)
      call plp%add(1.0e6_real64)
      call plp%add(.false., fmt="'B=',L1")
      call plp%end_row()
      call close_scratch(iu)
      call check_line(error, path, 1, repeat("_", 13)//"007"//repeat("_", 13)//"1.2"// &
                      repeat("_", 9)//"1.0E+06"//repeat("_", 13)//"L=T", "constructor format overrides")
      if (.not. allocated(error)) then
         call check_line(error, path, 2, repeat("_", 13)//"008"//repeat("_", 13)//"1.2"// &
                         repeat("_", 9)//"1.0E+06"//repeat("_", 13)//"B=F", &
                         "logical cell format overrides default")
      end if
      call discard_scratch(path)
      if (allocated(error)) return
      ! fmt_len=10 derives the default integer edit I3, which writes '***' for 1234
      call open_scratch(path, iu)
      plp = new_prettylistprinter([16], ["a"], unit=iu, fmt_len=10)
      call plp%add(1234)
      call plp%end_row()
      call close_scratch(iu)
      call check_line(error, path, 1, repeat("_", 13)//"1234", &
                      "an integer overflowing the derived edit width is re-rendered in full")
      call discard_scratch(path)
   end subroutine test_constructor_formats

   !> Narrow cells right-align what fits and print everything else in full
   !>
   !> Nothing is cut: a value wider than its column runs past it and the
   !> following cells shift right
   subroutine test_narrow_numeric_cells(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      type(prettylistprinter) :: plp
      integer :: iu, i
      character(len=*), parameter :: path = "plp_narrow.tmp"

      call open_scratch(path, iu)
      plp = new_prettylistprinter([5], ["a"], unit=iu, offset=0, column_gap=0)
      call plp%add(42_int8)
      call plp%end_row()
      call plp%add(42_int16)
      call plp%end_row()
      call plp%add(42_int32)
      call plp%end_row()
      call plp%add(42_int64, fmt="I20")
      call plp%end_row()
      plp = new_prettylistprinter([3, 5], ["a", "b"], unit=iu, offset=0, column_gap=0)
      call plp%add(42)
      call plp%add(123, fmt="I20")
      call plp%end_row()
      call plp%add("abcdef")
      call plp%add("ab")
      call plp%end_row()
      call plp%add(123456)
      call plp%add(-127_int8, fmt="I2")
      call plp%end_row()
      plp = new_prettylistprinter([5, 5], ["a", "b"], unit=iu, offset=0, column_gap=0)
      call plp%add(12.34_real32, fmt="F6.2")
      call plp%add(-1.25_real64, fmt="F6.2")
      call plp%end_row()
      plp = new_prettylistprinter([5, 5, 5], ["a", "b", "c"], unit=iu, offset=0, column_gap=1)
      call plp%add(-12.34_real64, fmt="F6.2")
      call plp%add(123.45_real32, fmt="F6.2")
      call plp%add(1.0_real64, fmt="F6.2")
      call plp%end_row()
      call close_scratch(iu)
      do i = 1, 4
         call check_line(error, path, i, "___42", "narrow fields retain all integer kinds and wide formats")
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         call check_line(error, path, 5, "_42__123", "each value is right-aligned in its own column width")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 6, "abcdef___ab", "a string wider than its column shifts the next cell")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 7, "123456_-127", &
                         "integers wider than their column or their I edit print in full")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 8, "12.34-1.25", &
                         "both real kinds retain the final digit at the cell boundary")
      end if
      if (.not. allocated(error)) then
         call check_line(error, path, 9, "-12.34_123.45__1.00", &
                         "six-character F6.2 values keep sign and leading digit in a five-wide cell")
      end if
      call discard_scratch(path)
   end subroutine test_narrow_numeric_cells

   !> Open a fresh scratch file for capturing printer output
   !>
   !> @param[in]  path Scratch file name, unique per test
   !> @param[out] iu   Unit connected to the truncated file
   subroutine open_scratch(path, iu)
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> Unit connected to the truncated file
      integer, intent(out) :: iu
      integer :: stat

      !$omp critical(moist_test_scratch_unit)
      open (newunit=iu, file=path, action="write", status="replace", iostat=stat)
      if (stat /= 0) error stop "prettylistprint test: open failed"
      !$omp end critical(moist_test_scratch_unit)
   end subroutine open_scratch

   !> Release the write unit once the printer is done with it
   !>
   !> @param[in] iu Unit the printer wrote to
   subroutine close_scratch(iu)
      !> Unit the printer wrote to
      integer, intent(in) :: iu
      integer :: stat

      !$omp critical(moist_test_scratch_unit)
      close (iu, iostat=stat)
      if (stat /= 0) error stop "prettylistprint test: close failed"
      !$omp end critical(moist_test_scratch_unit)
   end subroutine close_scratch

   !> Count the lines a printer wrote, then discard the scratch file
   !>
   !> @param[in]  path   Scratch file name
   !> @param[out] nlines Number of lines found
   subroutine count_lines(path, nlines)
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> Number of lines found
      integer, intent(out) :: nlines

      integer :: read_unit, stat

      nlines = 0
      !$omp critical(moist_test_scratch_unit)
      open (newunit=read_unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) error stop "prettylistprint test: read open failed"
      !$omp end critical(moist_test_scratch_unit)
      do
         read (read_unit, *, iostat=stat)
         if (stat == iostat_end) exit
         if (stat /= 0) error stop "prettylistprint test: count read failed"
         nlines = nlines + 1
      end do
      !$omp critical(moist_test_scratch_unit)
      close (read_unit, status="delete", iostat=stat)
      if (stat /= 0) error stop "prettylistprint test: delete close failed"
      !$omp end critical(moist_test_scratch_unit)
   end subroutine count_lines

   !> Delete a scratch file that `read_line` left behind
   !>
   !> `count_lines` discards the file itself; a test that only ever calls
   !> `read_line` has to clean up explicitly. Skipping it leaks the file into the
   !> working directory, which is the build tree under meson
   !>
   !> @param[in] path Scratch file name
   subroutine discard_scratch(path)
      !> Scratch file name
      character(len=*), intent(in) :: path

      integer :: unit, stat

      !$omp critical(moist_test_scratch_unit)
      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) error stop "prettylistprint test: delete open failed"
      close (unit, status="delete", iostat=stat)
      if (stat /= 0) error stop "prettylistprint test: delete failed"
      !$omp end critical(moist_test_scratch_unit)
   end subroutine discard_scratch

   !> Read one line of printer output, keeping the scratch file for later reads
   !>
   !> Blanks are rendered as '_' so that trailing spaces survive the comparison
   !> that `check` performs on trimmed strings; the record length is read explicitly
   !>
   !> @param[in]  path   Scratch file name
   !> @param[in]  iline  One-based line to return
   !> @param[out] line   Line content with blanks mapped to '_'
   !> @param[out] error  Error handle, set when the line cannot be read
   subroutine read_line(path, iline, line, error)
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> One-based line to return
      integer, intent(in) :: iline
      !> Line content with blanks mapped to '_'
      character(len=:), allocatable, intent(out) :: line
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      character(len=256) :: buf
      integer :: read_unit, stat, close_stat, i, j, nread

      line = ""
      buf = ""
      nread = 0
      !$omp critical(moist_test_scratch_unit)
      open (newunit=read_unit, file=path, action="read", status="old", iostat=stat)
      !$omp end critical(moist_test_scratch_unit)
      if (stat /= 0) then
         call test_failed(error, "cannot open printer output '"//path//"'")
         return
      end if
      do i = 1, iline
         read (read_unit, "(A)", advance="no", size=nread, iostat=stat) buf
         if (stat /= iostat_eor) exit
      end do
      !$omp critical(moist_test_scratch_unit)
      close (read_unit, iostat=close_stat)
      !$omp end critical(moist_test_scratch_unit)
      if (stat /= iostat_eor) then
         call test_failed(error, "printer output has no complete line "//to_string(iline), &
                          "read ended with iostat "//to_string(stat))
         return
      end if
      if (close_stat /= 0) then
         call test_failed(error, "cannot close printer output '"//path//"'")
         return
      end if

      line = buf(:nread)
      do j = 1, len(line)
         if (line(j:j) == " ") line(j:j) = "_"
      end do
   end subroutine read_line

   !> Compare one line of printer output with its expected text
   !>
   !> @param[out] error    Error handle
   !> @param[in]  path     Scratch file name
   !> @param[in]  iline    One-based line to compare
   !> @param[in]  expected Expected content with blanks written as '_'
   !> @param[in]  message  Failure message
   subroutine check_line(error, path, iline, expected, message)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> One-based line to compare
      integer, intent(in) :: iline
      !> Expected content with blanks written as '_'
      character(len=*), intent(in) :: expected
      !> Failure message
      character(len=*), intent(in) :: message

      character(len=:), allocatable :: line

      call read_line(path, iline, line, error)
      if (allocated(error)) return
      call check(error, line, expected, message, more="expected '"//expected//"' but got '"//line//"'")
   end subroutine check_line

end module test_utils_prettylistprint
