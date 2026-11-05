!> Exclusive file creation, replacement and removal helpers
!>
!> - Temporary files are claimed through native exclusive creation, independent
!>   of the Fortran runtime; interrupted files are left for the caller to remove
!> - Replacement uses rename on POSIX and replace-existing moves on Windows,
!>   without deleting or moving the previous file aside first
module moist_utils_file
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: error_type, fatal_error
   implicit none(type, external)
   private

   public :: create_new_file, replace_file, delete_file

   interface
      !> Claim a path exclusively, leaving an empty file on success
      !>
      !> @param[in] path Null-terminated path
      function c_create_exclusive(path) result(status) bind(C, name="moist_file_create_exclusive")
         import :: c_char, c_int
         implicit none(type, external)
         !> Null-terminated path
         ! allow(assumed-size): C strings carry no Fortran array descriptor
         character(kind=c_char), intent(in) :: path(*)
         !> Zero on success, one when the path exists, negative on another failure
         integer(c_int) :: status
      end function c_create_exclusive

      !> Replace a destination without deleting it first
      !>
      !> @param[in] source Null-terminated path of the completed file
      !> @param[in] target Null-terminated destination on the same filesystem
      function c_replace(source, target) result(status) bind(C, name="moist_file_replace")
         import :: c_char, c_int
         implicit none(type, external)
         !> Null-terminated completed file
         ! allow(assumed-size): C strings carry no Fortran array descriptor
         character(kind=c_char), intent(in) :: source(*)
         !> Null-terminated destination
         ! allow(assumed-size): C strings carry no Fortran array descriptor
         character(kind=c_char), intent(in) :: target(*)
         !> Zero on success
         integer(c_int) :: status
      end function c_replace
   end interface

contains

   !> Claim an empty file `<stem>.<n>.tmp` that did not exist before
   !>
   !> - The clock spreads candidate names; native exclusive creation arbitrates
   !>   collisions between threads and processes
   !>
   !> @param[in]  stem   path the name is built from
   !> @param[out] path   claimed empty file; reopen as old, without replacement; not a
   !>                    claimed file on error
   !> @param[out] error  no name could be claimed
   subroutine create_new_file(stem, path, error)
      !> Path the name is built from
      character(len=*), intent(in) :: stem
      !> Claimed empty file
      character(len=:), allocatable, intent(out) :: path
      !> No name could be claimed
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: max_attempts = 100
      character(len=24) :: buf
      integer(int64) :: start
      integer(c_int) :: status
      integer :: attempt

      call system_clock(start)
      do attempt = 0, max_attempts - 1
         write (buf, "(i0)") mod(max(start, 0_int64), 10_int64**15) + attempt
         path = stem//"."//trim(buf)//".tmp"
         status = c_create_exclusive(path//c_null_char)
         if (status == 0_c_int) return
         if (status /= 1_c_int) exit
      end do
      call fatal_error(error, "Failed to create a temporary file next to "//stem)
   end subroutine create_new_file

   !> Rename a completed file over its destination on the same filesystem
   !>
   !> @param[in]  source  completed file
   !> @param[in]  target  destination, replaced when it exists
   !> @param[out] error   failed replacement; the previous destination is retained
   subroutine replace_file(source, target, error)
      !> Completed file
      character(len=*), intent(in) :: source
      !> Destination on the same filesystem
      character(len=*), intent(in) :: target
      !> Failed replacement
      type(error_type), allocatable, intent(out) :: error

      if (c_replace(source//c_null_char, target//c_null_char) /= 0_c_int) then
         call fatal_error(error, "Failed to replace "//target//" by "//source)
      end if
   end subroutine replace_file

   !> Delete a file when it exists; a missing file is not an error
   !>
   !> @param[in]  file   file name
   !> @param[out] error  failed open or delete of an existing file
   subroutine delete_file(file, error)
      !> File name
      character(len=*), intent(in) :: file
      !> Failed delete
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: msg
      integer :: unit, stat
      logical :: exists

      inquire (file=file, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=file, status="old", iostat=stat, iomsg=msg)
      if (stat == 0) close (unit, status="delete", iostat=stat, iomsg=msg)
      if (stat /= 0) call fatal_error(error, "Failed to delete "//file//": "//trim(msg))
   end subroutine delete_file

end module moist_utils_file
