!> Concurrent key-value lines from independent pretty printers
!>
!> Run as a selected test, outside test-drive's OpenMP team
module test_utils_prettyprint_threaded
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_prettyprint, only: prettyprinter, new_prettyprinter
   implicit none(type, external)
   private

   public :: collect_utils_prettyprint_threaded

contains

   !> Collect concurrent printer tests
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_utils_prettyprint_threaded(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("independent_lines", test_independent_lines)]
   end subroutine collect_utils_prettyprint_threaded

   !> Concurrent formatting preserves every value and its width
   !>
   !> Mixed value kinds and widths expose shared character-result length temporaries
   subroutine test_independent_lines(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      integer :: i, failures

      failures = 0
      !$omp parallel do num_threads(4) reduction(+:failures)
      do i = 1, 40000
         if (.not. line_matches(i)) failures = failures + 1
      end do
      !$omp end parallel do
      call check(error, failures, 0, "independent printers preserve concurrent lines")
   end subroutine test_independent_lines

   !> Print one line to a private scratch unit and read it back
   !>
   !> @param[in] i Iteration selecting value kind and format width
   function line_matches(i) result(matches)
      !> Iteration selecting value kind and format width
      integer, intent(in) :: i
      !> Whether the line matches the expected text
      logical :: matches
      type(prettyprinter) :: pp
      character(128) :: line, expected
      integer :: iu, stat, width

      open (newunit=iu, status="scratch", action="readwrite", iostat=stat)
      if (stat /= 0) then
         matches = .false.
         return
      end if
      width = 11 + modulo(i / 4, 16)
      pp = new_prettyprinter(unit=iu, col_value=1, fmt_len=width)
      select case (modulo(i, 4))
      case (0)
         call pp%kv("x", 1.25_real64)
         expected = "x  "//repeat(" ", width - 8)//"1.250000"
      case (1)
         call pp%kv("x", 4567_int64)
         expected = "x  "//repeat(" ", width - 11)//"4567"
      case (2)
         call pp%kv("x", .true.)
         expected = "x  T"
      case default
         call pp%kv("x", repeat("z", width - 10))
         expected = "x  "//repeat("z", width - 10)
      end select
      rewind (iu)
      read (iu, "(A)", iostat=stat) line
      close (iu)
      matches = stat == 0 .and. line == expected
   end function line_matches

end module test_utils_prettyprint_threaded
