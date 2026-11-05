!> File creation and replacement helpers: utils_file
!>
!> - Temporary files of one stem are distinct while they exist; a stem in a missing
!>   directory is refused
!> - Native exclusive claims distinguish collisions between threads
!> - Replacement over an existing or a missing destination, leaving other files alone;
!>   the previous file is preserved on failure
module test_utils_file
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use mctc_env, only: moist_error => error_type
   use test_helpers, only: check_moist_error
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_utils_file, only: create_new_file, replace_file, delete_file
   implicit none(type, external)
   private

   public :: collect_utils_file

   !> Stem of the temporary-file test; the replacement tests run concurrently and
   !> name their files after themselves
   character(len=*), parameter :: stem = "moist-test-file-stem"

   interface
      !> Native creation entry point, exercised with deterministic name collisions
      !>
      !> @param[in] path Null-terminated file name
      function c_create_exclusive(path) result(status) bind(C, name="moist_file_create_exclusive")
         import :: c_char, c_int
         implicit none(type, external)
         !> Null-terminated file name
         ! allow(assumed-size): C strings carry no Fortran array descriptor
         character(kind=c_char), intent(in) :: path(*)
         !> Zero on success, one on a collision
         integer(c_int) :: status
      end function c_create_exclusive
   end interface

contains

   !> Collect the file helper tests
   !>
   !> @param[out] testsuite  collected tests
   subroutine collect_utils_file(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("new_file_distinct", test_new_file_distinct), &
         & new_unittest("new_file_refused", test_new_file_refused, should_fail=.true.), &
         & new_unittest("replace", test_replace), &
         & new_unittest("exclusive_collision", test_exclusive_collision), &
         & new_unittest("concurrent_claims", test_concurrent_claims), &
         & new_unittest("replace_missing_target", test_replace_missing_target), &
         & new_unittest("replace_failed", test_replace_failed), &
         & new_unittest("replace_missing_directory", test_replace_missing_directory), &
         & new_unittest("replace_preserves_other_files", test_replace_preserves_other_files) &
         & ]
   end subroutine collect_utils_file

   !> A second temporary file of the same stem gets another name while the first exists
   subroutine test_new_file_distinct(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: first, second
      logical :: exists_first, exists_second, first_owned

      call create_new_file(stem, first, err)
      first_owned = .not. allocated(err)
      if (.not. allocated(err)) call create_new_file(stem, second, err)
      call check(error, .not. allocated(err), "temporary file not created")
      if (allocated(error)) then
         if (first_owned) call remove(first)
         return
      end if
      inquire (file=first, exist=exists_first)
      inquire (file=second, exist=exists_second)
      call remove(first)
      call remove(second)
      call check(error, first /= second, "the second file reused the name of the first")
      if (allocated(error)) return
      call check(error, exists_first .and. exists_second, "a temporary file does not exist")
      if (allocated(error)) return
      call check(error, index(first, stem//".") == 1 .and. index(first, ".tmp") == len(first) - 3, &
         & "unexpected name "//first)
   end subroutine test_new_file_distinct

   !> Refused: a stem in a missing directory
   subroutine test_new_file_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: path

      call create_new_file("moist-test-missing-directory/file", path, err)
      call check_moist_error(error, err)
   end subroutine test_new_file_refused

   !> Rename replaces an existing destination
   subroutine test_replace(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: target = "moist-test-file-replace", &
         & source = "moist-test-file-replace-source"
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: text
      logical :: left

      call write_text(target, "old")
      call write_text(source, "new")
      call replace_file(source, target, err)
      call read_text(target, text)
      inquire (file=source, exist=left)
      call remove(target)
      call remove(source)
      call check(error, .not. allocated(err), "replacement failed")
      if (allocated(error)) return
      call check(error, text == "new" .and. .not. left, "the destination was not replaced")
   end subroutine test_replace

   !> Native exclusive creation refuses the same literal path without truncating it
   subroutine test_exclusive_collision(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: path = "moist-test-file-exclusive-collision"
      integer(c_int) :: first, second
      character(len=:), allocatable :: text

      call remove(path)
      first = c_create_exclusive(path//c_null_char)
      if (first == 0_c_int) call write_text(path, "sentinel")
      second = c_create_exclusive(path//c_null_char)
      call read_text(path, text)
      call remove(path)
      call check(error, first == 0_c_int .and. second == 1_c_int, "exclusive claims were not distinguished")
      if (allocated(error)) return
      call check(error, text == "sentinel", "a refused claim truncated the existing file")
   end subroutine test_exclusive_collision

   !> Only one concurrent claimant can create the same literal path
   subroutine test_concurrent_claims(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: path = "moist-test-file-concurrent-claims"
      integer, parameter :: nclaim = 32
      integer(c_int) :: status(nclaim)
      integer :: i
      logical :: exists

      call remove(path)
      !$omp parallel do schedule(static) num_threads(4)
      do i = 1, nclaim
         status(i) = c_create_exclusive(path//c_null_char)
      end do
      !$omp end parallel do
      inquire (file=path, exist=exists)
      call remove(path)
      call check(error, exists .and. count(status == 0_c_int) == 1, "the file had more or fewer than one owner")
      if (allocated(error)) return
      call check(error, count(status == 1_c_int) == nclaim - 1, "colliding claims did not report an existing file")
   end subroutine test_concurrent_claims

   !> Replacement creates a previously missing destination
   subroutine test_replace_missing_target(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: target = "moist-test-file-replace-missing", &
         & source = "moist-test-file-replace-missing-source"
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: text
      logical :: left

      call remove(target)
      call write_text(source, "new")
      call replace_file(source, target, err)
      call read_text(target, text)
      inquire (file=source, exist=left)
      call remove(target)
      call remove(source)
      call check(error, .not. allocated(err), "replacement of a missing destination failed")
      if (allocated(error)) return
      call check(error, text == "new" .and. .not. left, "the source was not moved to the missing destination")
   end subroutine test_replace_missing_target

   !> A missing source leaves the previous destination untouched
   subroutine test_replace_failed(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: target = "moist-test-file-replace-failed", &
         & source = "moist-test-file-replace-failed-source"
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: text

      call remove(source)
      call write_text(target, "old")
      call replace_file(source, target, err)
      call read_text(target, text)
      call remove(target)
      call check(error, allocated(err), "a missing source replaced the destination")
      if (allocated(error)) return
      call check(error, text == "old", "the failed replacement lost the previous destination")
   end subroutine test_replace_failed

   !> A missing destination directory leaves the source intact
   subroutine test_replace_missing_directory(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: source = "moist-test-file-replace-refused-source"
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: text

      call write_text(source, "new")
      call replace_file(source, "moist-test-missing-replace-directory/file", err)
      call read_text(source, text)
      call remove(source)
      call check(error, allocated(err), "a source was moved into a missing directory")
      if (allocated(error)) return
      call check(error, text == "new", "the failed replacement lost the source")
   end subroutine test_replace_missing_directory

   !> Replacement does not modify a neighboring backup file
   subroutine test_replace_preserves_other_files(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: target = "moist-test-file-replace-neighbor", &
         & source = "moist-test-file-replace-neighbor-source"
      type(moist_error), allocatable :: err
      character(len=:), allocatable :: text, backup

      call write_text(target, "old")
      call write_text(target//".bak", "unrelated")
      call write_text(source, "new")
      call replace_file(source, target, err)
      call read_text(target, text)
      call read_text(target//".bak", backup)
      call remove(target)
      call remove(target//".bak")
      call remove(source)
      call check(error, .not. allocated(err), "replacement failed")
      if (allocated(error)) return
      call check(error, text == "new" .and. backup == "unrelated", "replacement changed an unrelated file")
   end subroutine test_replace_preserves_other_files

   !> Write one line to a file, replacing it
   !>
   !> @param[in] path  file
   !> @param[in] text  line
   subroutine write_text(path, text)
      !> File
      character(len=*), intent(in) :: path
      !> Line
      character(len=*), intent(in) :: text
      integer :: unit

      open (newunit=unit, file=path, status="replace", action="write")
      write (unit, "(a)") text
      close (unit)
   end subroutine write_text

   !> First line of a file; empty when it cannot be read
   !>
   !> @param[in]  path  file
   !> @param[out] text  line
   subroutine read_text(path, text)
      !> File
      character(len=*), intent(in) :: path
      !> Line
      character(len=:), allocatable, intent(out) :: text
      character(len=64) :: line
      integer :: unit, stat

      text = ""
      open (newunit=unit, file=path, status="old", action="read", iostat=stat)
      if (stat /= 0) return
      read (unit, "(a)", iostat=stat) line
      if (stat == 0) text = trim(line)
      close (unit)
   end subroutine read_text

   !> Delete a file when it exists
   !>
   !> @param[in] path  file
   subroutine remove(path)
      !> File
      character(len=*), intent(in) :: path
      type(moist_error), allocatable :: err
      call delete_file(path, err)
   end subroutine remove

end module test_utils_file
