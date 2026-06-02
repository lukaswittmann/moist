
!> Canonical ASCII art headers for moist and its sub-models.
!> All banner output should go through this module so that the art
!> is defined in exactly one place.
module moist_output_ascii
   use moist_build_info, only: git_commit
   use moist_version, only: get_moist_version
   implicit none
   private

   !> Header style constants
   integer, parameter, public :: HEADER_FULL = 0  !< Logo + tagline (default)
   integer, parameter, public :: HEADER_SHORT = 1  !< Tagline box only
   integer, parameter, public :: HEADER_ASCII = 2  !< Logo only, no tagline

   public :: moist_header
   public :: moist_build_header
   public :: moist_version_header
   public :: gems_header, cavity_header

contains

   !> Print moist banner.
   !> @param[in] unit   Fortran I/O unit (6 = stdout)
   !> @param[in] style  Optional style selector:
   !>                    HEADER_FULL (0, default) = logo + tagline,
   !>                    HEADER_SHORT (1) = tagline box only,
   !>                    HEADER_ASCII (2) = logo only
   subroutine moist_header(unit, style)
      integer, intent(in) :: unit
      integer, intent(in), optional :: style

      integer :: s
      logical :: show_logo, show_tagline

      s = HEADER_FULL
      if (present(style)) s = style

      show_logo = (s == HEADER_FULL .or. s == HEADER_ASCII)
      show_tagline = (s == HEADER_FULL .or. s == HEADER_SHORT)

      if (show_logo) then
         write (unit, '(a)') &
            "                            _   _     _             ", &
            "                _ __ ___   / \ (_)___| |_           ", &
            "     .---------| '_ ` _ \ /   \| / __| __|---------.", &
            "     |         | | | | | |     | \__ \ |_          |", &
            "     |         |_| |_| |_|\___/|_|___/\__|         |", &
            "     |                                             |"
      end if

      if (show_tagline) then
         if (.not. show_logo) then
            write (unit, '(a)') &
               "     .---------------------------------------------."
         end if
         write (unit, '(a)') &
            "     |       Modular and open-source               |", &
            "     |            implicit solvation toolkit       |"
      end if

      write (unit, '(a)') &
         "     '---------------------------------------------'", ""

   end subroutine moist_header

   !> Print the moist build banner + version and commit
   subroutine moist_build_header(unit)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit
      character(len=:), allocatable :: version_string, line
      integer :: inner, pad_l, pad_r
      character(len=*), parameter :: top = &
                                     "     .---------------------------------------------."
      character(len=*), parameter :: bot = &
                                     "     '---------------------------------------------'"

      call get_moist_version(string=version_string)
      line = "moist v"//trim(version_string)//" ("//trim(git_commit)//")"

      ! Center the banner: inner width = box width minus 5 leading spaces and
      ! the two corner characters; derived from `top` so the two stay in sync.
      inner = len(top) - 7
      pad_l = max(0, (inner - len(line))/2)
      pad_r = max(0, inner - len(line) - pad_l)

      write (unit, '(a)') top
      write (unit, '(a)') "     |"//repeat(' ', pad_l)//line//repeat(' ', pad_r)//"|"
      write (unit, '(a)') bot, ""

   end subroutine moist_build_header

   !> Print one-line version string, e.g. "moist version 0.5.0"
   subroutine moist_version_header(unit)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit
      character(len=:), allocatable :: version_string

      call get_moist_version(string=version_string)
      write (unit, '(a, *(1x, a))') "moist", "version", version_string

   end subroutine moist_version_header

   !> Print GEMS solvation model banner.
   subroutine gems_header(unit)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit

      write (unit, '(a)') "", &
         "             _________                              ", &
         "            /__/___\__\    GEMS                     ", &
         "     .------'. \   / .'----------------------------.", &
         "     |        '.\ /.'    General and               |", &
         "     |          '.'    Minimally-Empirical         |", &
         "     |              Model for Solvation            |", &
         "     '---------------------------------------------'", ""

   end subroutine gems_header

   !> Print cavity construction banner.
   !> If scheme is present, includes the scheme name in the header.
   subroutine cavity_header(unit, scheme)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit
      !> Optional cavity scheme name (e.g. "DROP", "iSwiG")
      character(len=*), intent(in), optional :: scheme

      if (present(scheme)) then
         write (unit, '(a)') &
            "     .---------------------------------------------."
         write (unit, '(a,a6,a)') &
            "     |        Cavity Construction -- ", trim(scheme), "        |"
         write (unit, '(a)') &
            "     '---------------------------------------------'", ""
      else
         write (unit, '(a)') "", &
            "     .---------------------------------------------.", &
            "     |             Cavity Construction             |", &
            "     '---------------------------------------------'", ""
      end if

   end subroutine cavity_header

end module moist_output_ascii
