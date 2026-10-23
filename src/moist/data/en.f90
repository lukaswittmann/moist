module moist_data_en
   use mctc_env, only: wp
   use mctc_env, only: error_type, fatal_error
   use mctc_io_symbols, only: to_number
   implicit none(type, external)
   private

   public :: get_electronegativity
   public :: max_elem

   interface get_electronegativity
      module procedure :: get_electronegativity_num
      module procedure :: get_electronegativity_sym
   end interface get_electronegativity

   !> Highest atomic number this table covers
   integer, parameter :: max_elem = 118

   !> Pauling electronegativities, used for the covalent coordination number
   !>
   !> - Table from dftd4, Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)
   !>   DOI: 10.1063/1.5090222
   !> - Values from Allred, J. Inorg. Nucl. Chem. 17, 215 (1961)
   !>   DOI: 10.1016/0022-1902(61)80142-5
   !> - He, Ne, Ar, Rn, Pm, Eu, Tb, Yb, Fr and Rf-Og are not Pauling values
   real(wp), parameter :: pauling_en(max_elem) = [ &
      & 2.20_wp, 3.00_wp, 0.98_wp, 1.57_wp, 2.04_wp, 2.55_wp, 3.04_wp, 3.44_wp, &  ! H-O
      & 3.98_wp, 4.50_wp, 0.93_wp, 1.31_wp, 1.61_wp, 1.90_wp, 2.19_wp, 2.58_wp, &  ! F-S
      & 3.16_wp, 3.50_wp, 0.82_wp, 1.00_wp, 1.36_wp, 1.54_wp, 1.63_wp, 1.66_wp, &  ! Cl-Cr
      & 1.55_wp, 1.83_wp, 1.88_wp, 1.91_wp, 1.90_wp, 1.65_wp, 1.81_wp, 2.01_wp, &  ! Mn-Ge
      & 2.18_wp, 2.55_wp, 2.96_wp, 3.00_wp, 0.82_wp, 0.95_wp, 1.22_wp, 1.33_wp, &  ! As-Zr
      & 1.60_wp, 2.16_wp, 1.90_wp, 2.20_wp, 2.28_wp, 2.20_wp, 1.93_wp, 1.69_wp, &  ! Nb-Cd
      & 1.78_wp, 1.96_wp, 2.05_wp, 2.10_wp, 2.66_wp, 2.60_wp, 0.79_wp, 0.89_wp, &  ! In-Ba
      & 1.10_wp, 1.12_wp, 1.13_wp, 1.14_wp, 1.15_wp, 1.17_wp, 1.18_wp, 1.20_wp, &  ! La-Gd
      & 1.21_wp, 1.22_wp, 1.23_wp, 1.24_wp, 1.25_wp, 1.26_wp, 1.27_wp, 1.30_wp, &  ! Tb-Hf
      & 1.50_wp, 2.36_wp, 1.90_wp, 2.20_wp, 2.20_wp, 2.28_wp, 2.54_wp, 2.00_wp, &  ! Ta-Hg
      & 1.62_wp, 2.33_wp, 2.02_wp, 2.00_wp, 2.20_wp, 2.20_wp, 0.79_wp, 0.90_wp, &  ! Tl-Ra
      & 1.10_wp, 1.30_wp, 1.50_wp, 1.38_wp, 1.36_wp, 1.28_wp, 1.30_wp, 1.30_wp, &  ! Ac-Cm
      & 1.30_wp, 1.30_wp, 1.30_wp, 1.30_wp, 1.30_wp, 1.30_wp, 1.30_wp, 1.50_wp, &  ! Bk-Rf
      & 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, &  ! Db-Cn
      & 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp, 1.50_wp                    &  ! Nh-Og
      & ]

contains

   !> Get electronegativity for a given element symbol
   subroutine get_electronegativity_sym(sym, en, error)

      !> Element symbol
      character(len=*), intent(in) :: sym

      !> Electronegativity
      real(wp), intent(out) :: en

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer :: num

      en = 0.0_wp

      num = to_number(trim(sym))
      if (num < 1) then
         call fatal_error(error, "Unknown element symbol: '"//trim(sym)//"'")
         return
      end if

      call get_electronegativity_num(num, en, error)

   end subroutine get_electronegativity_sym

   !> Get electronegativity for a given atomic number
   subroutine get_electronegativity_num(num, en, error)

      !> Atomic number
      integer, intent(in) :: num

      !> Electronegativity
      real(wp), intent(out) :: en

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=64) :: msg

      en = 0.0_wp

      if (num < 1 .or. num > max_elem) then
         write (msg, "(a,i0,a,i0,a)") &
            "Atomic number ", num, " out of range [1, ", max_elem, "]"
         call fatal_error(error, trim(msg))
         return
      end if

      en = pauling_en(num)

   end subroutine get_electronegativity_num

end module moist_data_en
