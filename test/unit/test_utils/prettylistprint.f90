!> Test suite for the tabular printer in moist_utils_prettylistprint
module test_utils_prettylistprint
   use, intrinsic :: iso_fortran_env, only: real64
   use testdrive, only: new_unittest, unittest_type, error_type, check
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
   !> binary. What remains testable is the output the printer produces.
   subroutine collect_utils_prettylistprint(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("healthy-printer-emits", test_healthy_printer_emits), &
                  new_unittest("row-layout", test_row_layout), &
                  new_unittest("decorations", test_decorations), &
                  new_unittest("real-overflow-marker", test_real_overflow_marker) &
                  ]

   end subroutine collect_utils_prettylistprint

   !> A well-formed printer produces a header and one row
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
      call plp%add(2.5_real64, fmt='f6.2')
      call plp%end_row()

      call count_lines(iu, path, nlines)
      call check(error, nlines, 2, "header and one row were written")

   end subroutine test_healthy_printer_emits

   !> Values land right-aligned in their fields, and `skip` leaves a cell blank
   subroutine test_row_layout(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_layout.tmp"
      character(len=:), allocatable :: line
      integer :: iu

      call open_scratch(path, iu)

      plp = new_prettylistprinter([6, 6, 6], ["a", "b", "c"], unit=iu, offset=0, column_gap=0)
      call plp%begin_row()
      call plp%add("ab")
      call plp%skip()
      call plp%add(42, fmt='I3')
      call plp%end_row()

      call read_line(iu, path, 1, line)
      call check(error, line, "____ab__________42", "cells are right-aligned and the skipped one is blank")

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
      call plp%header("hi")
      call plp%separator()
      call plp%blank()

      call read_line(iu, path, 1, line)
      call check(error, line, "==_h_i_===", "the title is spread and centred in '=' fill")
      if (allocated(error)) return

      call read_line(iu, path, 2, line)
      call check(error, len_trim(line), 10, "the separator spans the table width")
      if (allocated(error)) return

      call count_lines(iu, path, nlines)
      call check(error, nlines, 3, "header, separator and blank were written")

   end subroutine test_decorations

   !> A real too wide for its cell is replaced by a signed overflow marker
   subroutine test_real_overflow_marker(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_overflow.tmp"
      character(len=:), allocatable :: line
      integer :: iu

      call open_scratch(path, iu)

      plp = new_prettylistprinter([5, 5], ["a", "b"], unit=iu, offset=0, column_gap=0)
      call plp%begin_row()
      call plp%add(1.0e12_real64, fmt='f6.2')
      call plp%add(-1.0e12_real64, fmt='f6.2')
      call plp%end_row()

      call read_line(iu, path, 1, line)
      call check(error, line, "+++++-----", "overflow markers carry the sign of the value")

   end subroutine test_real_overflow_marker

   !> Open a fresh scratch file for capturing printer output
   !>
   !> @param[in]  path Scratch file name, unique per test
   !> @param[out] iu   Unit connected to the truncated file
   subroutine open_scratch(path, iu)
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> Unit connected to the truncated file
      integer, intent(out) :: iu

      open (newunit=iu, file=path, action='write', status='replace')
   end subroutine open_scratch

   !> Count the lines a printer wrote, then discard the scratch file
   !>
   !> @param[in]  iu     Unit the printer wrote to
   !> @param[in]  path   Scratch file name
   !> @param[out] nlines Number of lines found
   subroutine count_lines(iu, path, nlines)
      !> Unit the printer wrote to
      integer, intent(in) :: iu
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> Number of lines found
      integer, intent(out) :: nlines

      integer :: read_unit, stat

      close (iu)
      nlines = 0
      open (newunit=read_unit, file=path, action='read', status='old')
      do
         read (read_unit, *, iostat=stat)
         if (stat /= 0) exit
         nlines = nlines + 1
      end do
      close (read_unit, status='delete')
   end subroutine count_lines

   !> Read one line of printer output, keeping the scratch file for later reads
   !>
   !> Blanks are rendered as '_' so that trailing spaces survive the comparison
   !> that `check` performs on trimmed strings.
   !>
   !> @param[in]  iu     Unit the printer wrote to
   !> @param[in]  path   Scratch file name
   !> @param[in]  iline  One-based line to return
   !> @param[out] line   Line content with blanks mapped to '_'
   subroutine read_line(iu, path, iline, line)
      !> Unit the printer wrote to
      integer, intent(in) :: iu
      !> Scratch file name
      character(len=*), intent(in) :: path
      !> One-based line to return
      integer, intent(in) :: iline
      !> Line content with blanks mapped to '_'
      character(len=:), allocatable, intent(out) :: line

      character(len=256) :: buf
      integer :: read_unit, stat, i, j
      logical :: is_open

      inquire (unit=iu, opened=is_open)
      if (is_open) close (iu)

      line = ''
      open (newunit=read_unit, file=path, action='read', status='old')
      do i = 1, iline
         read (read_unit, '(A)', iostat=stat) buf
         if (stat /= 0) then
            close (read_unit)
            return
         end if
      end do
      close (read_unit)

      line = trim(buf)
      do j = 1, len(line)
         if (line(j:j) == ' ') line(j:j) = '_'
      end do
   end subroutine read_line

end module test_utils_prettylistprint
