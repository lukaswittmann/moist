!> Concurrent rows from independent tabular printers
!>
!> Run as a selected test, outside test-drive's OpenMP team
module test_utils_prettylistprint_threaded
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none(type, external)
   private

   public :: collect_utils_prettylistprint_threaded

contains

   !> Collect concurrent printer tests
   !>
   !> @param[out] testsuite Collection of tests
   subroutine collect_utils_prettylistprint_threaded(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("independent_rows", test_independent_rows)]
   end subroutine collect_utils_prettylistprint_threaded

   !> Concurrent formatting preserves every cell and its width
   !>
   !> Varying widths exposes shared character-result length temporaries;
   !> inspect the row before file I/O to isolate formatting and concatenation
   !>
   !> @param[out] error Error handle
   subroutine test_independent_rows(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error
      integer :: i, failures

      failures = 0
      !$omp parallel do num_threads(4) reduction(+:failures)
      do i = 1, 100000
         if (.not. row_matches(4 + modulo(i, 16))) failures = failures + 1
      end do
      !$omp end parallel do
      call check(error, failures, 0, "independent printers preserve concurrent rows")
   end subroutine test_independent_rows

   !> Build one row without sharing a printer or output unit
   !>
   !> @param[in] width Width of the integer column
   function row_matches(width) result(matches)
      !> Width of the integer column
      integer, intent(in) :: width
      !> Whether the complete row matches the expected text
      logical :: matches
      type(prettylistprinter) :: plp
      character(128) :: expected

      plp = new_prettylistprinter([width, width + 3, 2, width + 5], ["i", "r", "l", "c"], &
                                 offset=0, column_gap=0, fmt_len=width + 7)
      call plp%add(4567_int64)
      call plp%add(1.25_real64, fmt="F4.2")
      call plp%add(.true.)
      call plp%add("z")
      expected = repeat(" ", width - 4)//"4567"//repeat(" ", width - 1)//"1.25 T"// &
                 repeat(" ", width + 4)//"z"
      matches = plp%row == expected
   end function row_matches

end module test_utils_prettylistprint_threaded
