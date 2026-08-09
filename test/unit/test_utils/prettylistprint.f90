!> Test suite for the tabular printer in moist_utils_prettylistprint
module test_utils_prettylistprint
   use, intrinsic :: iso_fortran_env, only: real64
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none(type, external)
   private

   public :: collect_utils_prettylistprint

contains

   !> Collect all tabular-printer tests
   subroutine collect_utils_prettylistprint(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("healthy-printer-emits", test_healthy_printer_emits), &
                  new_unittest("width-header-mismatch", test_width_header_mismatch), &
                  new_unittest("zero-columns", test_zero_columns), &
                  new_unittest("too-many-values", test_too_many_values), &
                  new_unittest("real-overrun", test_real_overrun), &
                  new_unittest("missing-columns", test_missing_columns), &
                  new_unittest("faulted-printer-inert", test_faulted_printer_inert), &
                  new_unittest("first-fault-wins", test_first_fault_wins) &
                  ]

   end subroutine collect_utils_prettylistprint

   !> A well-formed printer stays healthy and produces output
   subroutine test_healthy_printer_emits(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      type(moist_error_type), allocatable :: latched
      character(len=*), parameter :: path = "plp_healthy.tmp"
      integer :: iu, nlines

      call open_scratch(path, iu)

      plp = new_prettylistprinter([8, 8], ["alpha", "beta "], unit=iu)
      call plp%print_header()
      call plp%begin_row()
      call plp%add(1)
      call plp%add(2.5_real64, fmt='f6.2')
      call plp%end_row()

      call check(error,.not. plp%has_fault(), "well-formed printer stays healthy")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      !> `check` on a healthy printer must leave the error alone, so that a
      !> caller can invoke it unconditionally after building a table
      call plp%check(latched)
      call check(error,.not. allocated(latched), "check leaves a healthy printer silent")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      call count_lines(iu, path, nlines)
      call check(error, nlines, 2, "header and one row were written")

   end subroutine test_healthy_printer_emits

   !> Mismatched widths and headers fault the printer instead of aborting
   subroutine test_width_header_mismatch(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      type(moist_error_type), allocatable :: latched

      plp = new_prettylistprinter([8, 8, 8], ["alpha", "beta "])

      call check(error, plp%has_fault(), "mismatch faulted the printer")
      if (allocated(error)) return

      call plp%check(latched)
      call check(error, allocated(latched), "check surfaces the fault")
      if (allocated(error)) return
      call check(error, index(latched%message, "widths and headers") > 0, &
                 "message names the mismatch")

   end subroutine test_width_header_mismatch

   !> A printer with no columns is refused
   subroutine test_zero_columns(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      type(moist_error_type), allocatable :: latched
      integer, allocatable :: widths(:)
      character(len=4), allocatable :: headers(:)

      allocate (widths(0))
      allocate (headers(0))
      plp = new_prettylistprinter(widths, headers)

      call check(error, plp%has_fault(), "empty column list faulted the printer")
      if (allocated(error)) return

      call plp%check(latched)
      call check(error, allocated(latched), "check surfaces the fault")
      if (allocated(error)) return
      call check(error, index(latched%message, "at least one column") > 0, &
                 "message states the requirement")

   end subroutine test_zero_columns

   !> Overrunning a row faults rather than writing past the last column
   subroutine test_too_many_values(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      type(moist_error_type), allocatable :: latched
      character(len=*), parameter :: path = "plp_overrun.tmp"
      integer :: iu, nlines

      call open_scratch(path, iu)

      plp = new_prettylistprinter([6, 6], ["a", "b"], unit=iu)
      call plp%begin_row()
      call plp%add(1)
      call plp%add(2)
      call check(error,.not. plp%has_fault(), "a full row is not yet a fault")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      call plp%add(3)
      call check(error, plp%has_fault(), "the value past the last column faulted")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      call plp%check(latched)
      call check(error, index(latched%message, "too many values") > 0, &
                 "message names the overrun")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      !> The row was never completed, so nothing may have reached the file.
      call plp%end_row()
      call count_lines(iu, path, nlines)
      call check(error, nlines, 0, "no output escaped from the overrun row")

   end subroutine test_too_many_values

   !> The real specifics index the column arrays themselves, so they need their
   !> own overrun check rather than relying on the shared string path
   subroutine test_real_overrun(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_real_overrun.tmp"
      integer :: iu, nlines

      call open_scratch(path, iu)

      plp = new_prettylistprinter([8, 8], ["a", "b"], unit=iu)
      call plp%begin_row()
      call plp%add(1.0_real64, fmt='f6.2')
      call plp%add(2.0_real64, fmt='f6.2')
      call plp%add(3.0_real64, fmt='f6.2')

      call check(error, plp%has_fault(), "the real past the last column faulted")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      call plp%end_row()
      call count_lines(iu, path, nlines)
      call check(error, nlines, 0, "no output escaped from the overrun row")

   end subroutine test_real_overrun

   !> Ending a short row faults and prints nothing
   subroutine test_missing_columns(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      type(moist_error_type), allocatable :: latched
      character(len=*), parameter :: path = "plp_short.tmp"
      integer :: iu, nlines

      call open_scratch(path, iu)

      plp = new_prettylistprinter([6, 6, 6], ["a", "b", "c"], unit=iu)
      call plp%begin_row()
      call plp%add(1)
      call plp%end_row()

      call check(error, plp%has_fault(), "short row faulted the printer")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      call plp%check(latched)
      call check(error, index(latched%message, "missing columns") > 0, &
                 "message names the short row")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      !> A truncated line in the middle of a table is harder to diagnose than
      !> no line at all, so the malformed row must not be emitted.
      call count_lines(iu, path, nlines)
      call check(error, nlines, 0, "the malformed row was not printed")

   end subroutine test_missing_columns

   !> Once faulted, a printer emits nothing and tolerates any further calls
   subroutine test_faulted_printer_inert(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      character(len=*), parameter :: path = "plp_inert.tmp"
      integer :: iu, nlines

      call open_scratch(path, iu)

      !> Faulted at construction, so its width, header and row arrays were never
      !> allocated: every call below has to bail out before touching them.
      plp = new_prettylistprinter([8, 8, 8], ["alpha", "beta "], unit=iu)
      call check(error, plp%has_fault(), "printer faulted at construction")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if

      call plp%header("title")
      call plp%print_header()
      call plp%separator()
      call plp%blank()
      call plp%begin_row()
      call plp%add(1)
      call plp%add("text")
      call plp%skip()
      call plp%end_row()
      call plp%set_column_gap(4)

      call count_lines(iu, path, nlines)
      call check(error, nlines, 0, "a faulted printer wrote nothing")

   end subroutine test_faulted_printer_inert

   !> Only the first fault is kept; later ones are its consequences
   subroutine test_first_fault_wins(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(prettylistprinter) :: plp
      type(moist_error_type), allocatable :: latched
      character(len=*), parameter :: path = "plp_first.tmp"
      integer :: iu

      call open_scratch(path, iu)

      plp = new_prettylistprinter([6, 6], ["a", "b"], unit=iu)
      call plp%begin_row()
      call plp%add(1)

      !> Short row first ...
      call plp%end_row()
      !> ... then an overrun that would latch a different message.
      call plp%add(2)
      call plp%add(3)
      call plp%add(4)

      call plp%check(latched)
      call check(error, allocated(latched), "a fault is still reported")
      if (allocated(error)) then
         call close_scratch(iu, path)
         return
      end if
      call check(error, index(latched%message, "missing columns") > 0, &
                 "the first fault is the one retained")

      call close_scratch(iu, path)

   end subroutine test_first_fault_wins

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

   !> Close and remove a scratch file.
   !>
   !> @param[in] iu   Unit to close
   !> @param[in] path Scratch file name
   subroutine close_scratch(iu, path)
      !> Unit to close
      integer, intent(in) :: iu
      !> Scratch file name
      character(len=*), intent(in) :: path

      logical :: is_open

      inquire (file=path, opened=is_open)
      if (is_open) close (iu, status='delete')
   end subroutine close_scratch

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

end module test_utils_prettylistprint
