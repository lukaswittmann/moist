
!> Canonical ASCII art headers for moist and its sub-models
!>
!> All banner output should go through this module, so that the art is
!> defined in exactly one place
module moist_output_ascii
   use moist_build_info, only: git_commit, build_host
   use moist_version, only: get_moist_version
   implicit none(type, external)
   private

   !> Banner style constants
   integer, parameter, public :: BANNER_FULL = 0  !< Boxed logo + tagline (default)
   integer, parameter, public :: BANNER_TAGLINE = 1  !< Tagline box only
   integer, parameter, public :: BANNER_LOGO = 2  !< Logo box only
   integer, parameter, public :: BANNER_BUILD = 3  !< Version, commit and host box
   integer, parameter, public :: BANNER_COMPACT = 4  !< Boxless logo, tagline and version beside it

   public :: moist_banner_text
   public :: moist_header
   public :: moist_build_header
   public :: moist_version_header
   public :: cavity_header

contains

   !> Return the canonical banner as text with newline-terminated lines
   !>
   !> @param[in] style One of the BANNER_* constants
   function moist_banner_text(style) result(text)
      !> Banner selector
      integer, intent(in) :: style
      !> Host-printable banner text
      character(len=:), allocatable :: text
      !> Release and build metadata
      character(len=:), allocatable :: version, line
      !> Newline separator
      character(len=1), parameter :: nl = achar(10)
      text = ""
      if (style == BANNER_FULL .or. style == BANNER_LOGO) then
         text = "                            _   _     _             "//nl// &
            "                _ __ ___   / \ (_)___| |_           "//nl// &
            "     .---------| '_ ` _ \ /   \| / __| __|---------."//nl// &
            "     |         | | | | | |     | \__ \ |_          |"//nl// &
            "     |         |_| |_| |_|\___/|_|___/\__|         |"//nl// &
            "     |                                             |"//nl
      end if
      if (style == BANNER_TAGLINE) text = "     .---------------------------------------------."//nl
      if (style == BANNER_FULL .or. style == BANNER_TAGLINE) then
         text = text//"     |       Modular and Open-source               |"//nl// &
            "     |            Implicit Solvation Toolkit       |"//nl
      end if
      if (style >= BANNER_FULL .and. style <= BANNER_LOGO) then
         text = text//"     '---------------------------------------------'"//nl//nl
      end if
      if (style == BANNER_BUILD .or. style == BANNER_COMPACT) then
         call get_moist_version(string=version)
         line = "v"//version//" ("//trim(git_commit)//")"
      end if
      if (style == BANNER_BUILD) then
         text = "     .---------------------------------------------."//nl//boxed("moist "//line)//nl
         if (len_trim(build_host) > 0) text = text//boxed(trim(build_host))//nl
         text = text//"     '---------------------------------------------'"//nl//nl
      end if
      if (style == BANNER_COMPACT) then
         text = "                  _   _     _"//nl// &
            "      _ __ ___   / \ (_)___| |_     Modular and Open-source"//nl// &
            "     | '_ ` _ \ /   \| / __| __|    Implicit Solvation Toolkit"//nl// &
            "     | | | | | |     | \__ \ |_"//nl// &
            "     |_| |_| |_|\___/|_|___/\__|    "//line//nl
         if (len_trim(build_host) > 0) text = text//repeat(" ", 36)//trim(build_host)//nl
         text = text//nl
      end if
   contains
      !> Center one build metadata line inside the banner
      function boxed(value) result(row)
         character(len=*), intent(in) :: value
         character(len=:), allocatable :: row
         integer :: left, right
         left = max(0, (45-len(value))/2)
         right = max(0, 45-len(value)-left)
         row = "     |"//repeat(" ",left)//value//repeat(" ",right)//"|"
      end function boxed
   end function moist_banner_text

   !> Print a moist banner, BANNER_FULL unless `style` is given
   subroutine moist_header(unit, style)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit
      !> One of the BANNER_* constants
      integer, intent(in), optional :: style
      !> Resolved banner style
      integer :: selected
      selected = BANNER_FULL
      if (present(style)) selected = style
      write(unit, "(a)", advance="no") moist_banner_text(selected)
   end subroutine moist_header

   !> Print the moist build banner + version and commit
   subroutine moist_build_header(unit)
      !> Destination Fortran unit
      integer, intent(in) :: unit
      write(unit, "(a)", advance="no") moist_banner_text(BANNER_BUILD)
   end subroutine moist_build_header

   !> Print one-line version string, e.g. "moist version 0.5.0"
   subroutine moist_version_header(unit)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit
      character(len=:), allocatable :: version_string

      call get_moist_version(string=version_string)
      write (unit, "(a, *(1x, a))") "moist", "version", version_string

   end subroutine moist_version_header

   !> Print cavity construction banner
   !>
   !> - a present `scheme` puts the scheme name in the header
   subroutine cavity_header(unit, scheme)
      !> Fortran I/O unit (6 = stdout)
      integer, intent(in) :: unit
      !> Optional cavity scheme name (e.g. "DROP", "iSwiG")
      character(len=*), intent(in), optional :: scheme

      if (present(scheme)) then
         write (unit, "(a)") &
            "     .---------------------------------------------."
         write (unit, "(a,a6,a)") &
            "     |        Cavity Construction -- ", trim(scheme), "        |"
         write (unit, "(a)") &
            "     '---------------------------------------------'", ""
      else
         write (unit, "(a)") "", &
            "     .---------------------------------------------.", &
            "     |             Cavity Construction             |", &
            "     '---------------------------------------------'", ""
      end if

   end subroutine cavity_header

end module moist_output_ascii
