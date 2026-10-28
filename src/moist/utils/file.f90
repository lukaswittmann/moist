!> File replacement and removal helpers
!>
!> - `replace_file` renames a finished file over its destination, so readers
!>   see either the old or the new contents, never a partial write
!> - Atomic on POSIX file systems for files in the same directory; Windows
!>   refuses to rename over an existing file, so there the destination is
!>   removed first
module moist_utils_file
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use mctc_env, only: error_type, fatal_error
   use mctc_env_system, only: is_windows
   implicit none(type, external)
   private

   public :: replace_file, delete_file

   interface
      !> C library rename
      !>
      !> @param[in] source  null-terminated path of the existing file
      !> @param[in] target  null-terminated path of the destination
      function c_rename(source, target) result(status) bind(C, name="rename")
         import :: c_char, c_int
         implicit none(type, external)
         !> Existing file
         character(kind=c_char), intent(in) :: source(*)
         !> Destination
         character(kind=c_char), intent(in) :: target(*)
         !> Zero on success
         integer(c_int) :: status
      end function c_rename
   end interface

contains

   !> Rename a file over its destination
   !>
   !> @param[in]  source  existing file, in the directory of `target`
   !> @param[in]  target  destination, replaced when it exists
   !> @param[out] error   failed rename; `source` is left in place
   subroutine replace_file(source, target, error)
      !> Existing file
      character(len=*), intent(in) :: source
      !> Destination
      character(len=*), intent(in) :: target
      !> Failed rename
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: status
      logical :: exists

      status = c_rename(source//c_null_char, target//c_null_char)
      if (status /= 0 .and. is_windows()) then
         inquire (file=target, exist=exists)
         if (exists) then
            call delete_file(target, error)
            if (allocated(error)) return
            status = c_rename(source//c_null_char, target//c_null_char)
         end if
      end if
      if (status /= 0) call fatal_error(error, "Failed to replace "//target//" by "//source)
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
