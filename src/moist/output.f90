module moist_output

   use moist_output_format, only: format_string, getline
   use moist_output_ascii, only: HEADER_FULL, HEADER_SHORT, HEADER_ASCII, moist_banner_text, &
      & moist_header, moist_build_header, moist_version_header, cavity_header

   implicit none(type, external)
   private

   public :: format_string, getline
   public :: HEADER_FULL, HEADER_SHORT, HEADER_ASCII, moist_banner_text
   public :: moist_header, moist_build_header, moist_version_header, cavity_header

contains

end module moist_output
