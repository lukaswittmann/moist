module moist_output

   use moist_output_format, only: format_string, getline
   use moist_output_ascii, only: BANNER_FULL, BANNER_TAGLINE, BANNER_LOGO, BANNER_BUILD, BANNER_COMPACT, &
      & moist_banner_text, moist_header, moist_build_header, moist_version_header, cavity_header

   implicit none(type, external)
   private

   public :: format_string, getline
   public :: BANNER_FULL, BANNER_TAGLINE, BANNER_LOGO, BANNER_BUILD, BANNER_COMPACT
   public :: moist_banner_text
   public :: moist_header, moist_build_header, moist_version_header, cavity_header

contains

end module moist_output
