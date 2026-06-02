! Committed fallback build-info module so that a bare `fpm build` resolves the
! moist_build_info dependency without any generation step. The commit string is
! the static placeholder "unknown"; fpm builds report this value.
!
! Meson builds instead regenerate this module out-of-tree from build_info.f90.in
! via vcs_tag() (see src/moist/meson.build), stamping in the live short commit;
! the tracked file here is left untouched. Keep this module's interface in sync
! with build_info.f90.in.

!> Build-time provenance for the moist library.
module moist_build_info
   implicit none
   private

   public :: git_commit

   !> Short git commit hash of the build, or "unknown" when the commit is not
   !> available at build time (release tarball, or a bare `fpm build`).
   character(len=*), parameter :: git_commit = "unknown"

end module moist_build_info
