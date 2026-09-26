
!> Canonical ASCII art headers for moist and its sub-models
!>
!> All banner output should go through this module, so that the art is
!> defined in exactly one place
module moist_output_ascii
   use moist_build_info, only: git_commit, build_host
   use moist_version, only: get_moist_version
   implicit none(type, external)
   private

   !> Header style constants
   integer, parameter, public :: HEADER_FULL = 0  !< Logo + tagline (default)
   integer, parameter, public :: HEADER_SHORT = 1  !< Tagline box only
   integer, parameter, public :: HEADER_ASCII = 2  !< Logo only, no tagline

   public :: moist_banner_text
   public :: moist_header
   public :: moist_build_header
   public :: moist_version_header
   public :: cavity_header

contains

   !> Print moist banner
   !>
   !> @param[in] unit   Fortran I/O unit (6 = stdout)
   !> @param[in] style  Optional style selector:
   !>                    HEADER_FULL (0, default) = logo + tagline,
   !>                    HEADER_SHORT (1) = tagline box only,
   !>                    HEADER_ASCII (2) = logo only
   !>
   !> Return the canonical banner as text with newline-terminated lines
   !> @param[in] style Full (0), short (1), ASCII (2), or build information (3)
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
      if (style == HEADER_FULL .or. style == HEADER_ASCII) then
         text = "                            _   _     _             "//nl// &
            "                _ __ ___   / \ (_)___| |_           "//nl// &
            "     .---------| '_ ` _ \ /   \| / __| __|---------."//nl// &
            "     |         | | | | | |     | \__ \ |_          |"//nl// &
            "     |         |_| |_| |_|\___/|_|___/\__|         |"//nl// &
            "     |                                             |"//nl
      end if
      if (style == HEADER_SHORT) text = "     .---------------------------------------------."//nl
      if (style == HEADER_FULL .or. style == HEADER_SHORT) then
         text = text//"     |       Modular and open-source               |"//nl// &
            "     |            implicit solvation toolkit       |"//nl
      end if
      if (style >= 0 .and. style <= 2) text = text//"     '---------------------------------------------'"//nl// &
            ""//nl
      if (style == 3) then
         call get_moist_version(string=version)
         line = "moist v"//version//" ("//trim(git_commit)//")"
         text = "     .---------------------------------------------."//nl//boxed(line)//nl
         if (len_trim(build_host) > 0) text = text//boxed(trim(build_host))//nl
         text = text//"     '---------------------------------------------'"//nl//nl
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

   subroutine moist_header(unit, style)
      integer, intent(in) :: unit
      integer, intent(in), optional :: style
      integer :: selected
      selected = HEADER_FULL
      if (present(style)) selected = style
      write(unit, "(a)", advance="no") moist_banner_text(selected)
   end subroutine moist_header

   !> Print the moist build banner + version and commit
   subroutine moist_build_header(unit)
      !> Destination Fortran unit
      integer, intent(in) :: unit
      write(unit, "(a)", advance="no") moist_banner_text(3)
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
