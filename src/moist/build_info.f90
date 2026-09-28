! Committed fallback build-info module, so a bare `fpm build` resolves the
! moist_build_info dependency without any generation step; the commit string
! is the static placeholder "unknown", which fpm builds report
!
! Meson builds instead regenerate this module out-of-tree from
! build_info.f90.in via vcs_tag() (see src/moist/meson.build), stamping in the
! live short commit and leaving the tracked file here untouched
!
! Keep this module's interface in sync with build_info.f90.in

!> Build-time provenance for the moist library
module moist_build_info
   implicit none(type, external)
   private

   public :: git_commit, build_host

   !> Short git commit hash of the build, or "unknown" when the commit is not
   !> available at build time (release tarball, or a bare `fpm build`)
   character(len=*), parameter :: git_commit = "unknown"

   !> Host program this copy of moist was built as part of
   character(len=*), parameter :: build_host = ""

end module moist_build_info
