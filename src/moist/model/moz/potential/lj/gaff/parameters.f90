!> GAFF2 LJ table
!>
!> Wang et al., J. Comput. Chem. 25, 1157 (2004), doi:10.1002/jcc.20035
module moist_model_moz_potential_lj_gaff_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kcaltoau
   use mctc_io_utils, only: to_lower
   implicit none(type, external)
   private

   public :: lookup_gaff_lj, gaff_parameter_version

   ! GAFF2 Parameters taken from openmmforcefields
   !
   ! MIT License
   !
   ! Copyright (c) 2016-2019 by
   !     Chodera lab // Memorial Sloan Kettering Cancer Center
   !     Pande group // Stanford University
   !
   ! Permission is hereby granted, free of charge, to any person obtaining a copy
   ! of this software and associated documentation files (the "Software"), to deal
   ! in the Software without restriction, including without limitation the rights
   ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   ! copies of the Software, and to permit persons to whom the Software is
   ! furnished to do so, subject to the following conditions:
   !
   ! The above copyright notice and this permission notice shall be included in all
   ! copies or substantial portions of the Software.
   !
   ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   ! SOFTWARE.

   !> Version of the GAFF2 nonbonded parameter table
   character(len=*), parameter :: gaff_parameter_version = "2.2.30"

   !> Number of GAFF2 rows including the two water placeholders
   integer, parameter, public :: n_gaff = 100

   !> GAFF2 atomtype keys and legacy water placeholders
   character(len=2), parameter :: gaff_keys(n_gaff) = [ &
      & "c ", "cs", "c1", "c2", "c3", "ca", "cp", "cq", &
      & "cc", "cd", "ce", "cf", "cg", "ch", "cx", "cy", &
      & "cu", "cv", "cz", "h1", "h2", "h3", "h4", "h5", &
      & "ha", "hc", "hn", "ho", "hp", "hs", "hw", "hx", &
      & "f ", "cl", "br", "i ", "n ", "n1", "n2", "n3", &
      & "n4", "na", "nb", "nc", "nd", "ne", "nf", "nh", &
      & "no", "ns", "nt", "nx", "ny", "nz", "n+", "nu", &
      & "nv", "n7", "n8", "n9", "ni", "nj", "nk", "nl", &
      & "nm", "nn", "np", "nq", "n5", "n6", "o ", "oh", &
      & "op", "oq", "os", "ow", "p2", "p3", "p4", "p5", &
      & "pb", "pc", "pd", "pe", "pf", "px", "py", "s ", &
      & "s2", "s4", "s6", "sh", "sp", "sq", "ss", "sx", &
      & "sy", "hb", "c5", "c6" &
      & ]

   !> GAFF2 Rmin/2 in Angstrom; zero for water placeholders
   real(wp), parameter :: gaff_rmin_half_aa(n_gaff) = [ &
      & 1.8606_wp, 1.8606_wp, 1.9525_wp, 1.8606_wp, &
      & 1.9069_wp, 1.8606_wp, 1.8606_wp, 1.8606_wp, &
      & 1.8606_wp, 1.8606_wp, 1.8606_wp, 1.8606_wp, &
      & 1.9525_wp, 1.9525_wp, 1.9069_wp, 1.9069_wp, &
      & 1.8606_wp, 1.8606_wp, 1.8606_wp, 1.3593_wp, &
      & 1.2593_wp, 1.1593_wp, 1.4235_wp, 1.3735_wp, &
      & 1.4735_wp, 1.4593_wp, 0.6210_wp, 0.3019_wp, &
      & 0.6031_wp, 0.6112_wp, 0.0000_wp, 1.0593_wp, &
      & 1.7029_wp, 1.9452_wp, 2.0275_wp, 2.1558_wp, &
      & 1.7852_wp, 1.8372_wp, 1.8993_wp, 1.8886_wp, &
      & 1.4028_wp, 1.7992_wp, 1.8993_wp, 1.8993_wp, &
      & 1.8993_wp, 1.8993_wp, 1.8993_wp, 1.7903_wp, &
      & 1.8886_wp, 1.8352_wp, 1.8852_wp, 1.4528_wp, &
      & 1.5028_wp, 1.5528_wp, 1.6028_wp, 1.8403_wp, &
      & 1.8903_wp, 1.9686_wp, 2.0486_wp, 2.2700_wp, &
      & 1.7852_wp, 1.7852_wp, 1.4528_wp, 1.4528_wp, &
      & 1.7903_wp, 1.7903_wp, 1.8886_wp, 1.8886_wp, &
      & 1.9686_wp, 1.9686_wp, 1.7107_wp, 1.8200_wp, &
      & 1.7713_wp, 1.7713_wp, 1.7713_wp, 0.0000_wp, &
      & 2.0732_wp, 2.0732_wp, 2.0732_wp, 2.0732_wp, &
      & 2.0732_wp, 2.0732_wp, 2.0732_wp, 2.0732_wp, &
      & 2.0732_wp, 2.0732_wp, 2.0732_wp, 1.9825_wp, &
      & 1.9825_wp, 1.9825_wp, 2.2777_wp, 1.9825_wp, &
      & 1.9825_wp, 1.9825_wp, 1.9825_wp, 1.9825_wp, &
      & 2.2777_wp, 1.4735_wp, 1.9069_wp, 1.9069_wp &
      & ]

   !> GAFF2 well depth, kcal/mol; zero for water placeholders
   real(wp), parameter :: gaff_epsilon_kcal(n_gaff) = [ &
      & 0.0988_wp, 0.0988_wp, 0.1596_wp, 0.0988_wp, &
      & 0.1078_wp, 0.0988_wp, 0.0988_wp, 0.0988_wp, &
      & 0.0988_wp, 0.0988_wp, 0.0988_wp, 0.0988_wp, &
      & 0.1596_wp, 0.1596_wp, 0.1078_wp, 0.1078_wp, &
      & 0.0988_wp, 0.0988_wp, 0.0988_wp, 0.0208_wp, &
      & 0.0208_wp, 0.0208_wp, 0.0161_wp, 0.0161_wp, &
      & 0.0161_wp, 0.0208_wp, 0.0100_wp, 0.0047_wp, &
      & 0.0144_wp, 0.0124_wp, 0.0000_wp, 0.0208_wp, &
      & 0.0832_wp, 0.2638_wp, 0.3932_wp, 0.4955_wp, &
      & 0.1636_wp, 0.1098_wp, 0.0941_wp, 0.0858_wp, &
      & 3.8748_wp, 0.2042_wp, 0.0941_wp, 0.0941_wp, &
      & 0.0941_wp, 0.0941_wp, 0.0941_wp, 0.2150_wp, &
      & 0.0858_wp, 0.1174_wp, 0.0851_wp, 2.5453_wp, &
      & 1.6959_wp, 1.1450_wp, 0.7828_wp, 0.1545_wp, &
      & 0.1120_wp, 0.0522_wp, 0.0323_wp, 0.0095_wp, &
      & 0.1636_wp, 0.1636_wp, 2.5453_wp, 2.5453_wp, &
      & 0.2150_wp, 0.2150_wp, 0.0858_wp, 0.0858_wp, &
      & 0.0522_wp, 0.0522_wp, 0.1463_wp, 0.0930_wp, &
      & 0.0726_wp, 0.0726_wp, 0.0726_wp, 0.0000_wp, &
      & 0.2295_wp, 0.2295_wp, 0.2295_wp, 0.2295_wp, &
      & 0.2295_wp, 0.2295_wp, 0.2295_wp, 0.2295_wp, &
      & 0.2295_wp, 0.2295_wp, 0.2295_wp, 0.2824_wp, &
      & 0.2824_wp, 0.2824_wp, 0.0614_wp, 0.2824_wp, &
      & 0.2824_wp, 0.2824_wp, 0.2824_wp, 0.2824_wp, &
      & 0.0614_wp, 0.0161_wp, 0.1078_wp, 0.1078_wp &
      & ]

contains

   !> Two-character lowercase key of a type or element string
   !>
   !> @param[in] s Type or element string
   pure function normalize_key(s) result(out)
      !> Type or element string
      character(len=*), intent(in) :: s
      !> Normalized key
      character(len=2) :: out
      out = to_lower(adjustl(s))
   end function normalize_key

   !> Row of a key in a key table, 0 when absent
   !>
   !> @param[in] key Normalized key
   !> @param[in] keys Key table
   pure function find_key(key, keys) result(idx)
      !> Normalized key
      character(len=2), intent(in) :: key
      !> Key table
      character(len=2), intent(in) :: keys(:)
      !> Row, 0 when absent
      integer :: idx
      integer :: i

      idx = 0
      do i = 1, size(keys)
         if (keys(i) == key) then
            idx = i
            return
         end if
      end do
   end function find_key

   !> Lennard-Jones parameters of a GAFF2 type
   !>
   !> @param[in] gaff_type GAFF atom type
   !> @param[out] sigma Sigma, bohr; zero when not found
   !> @param[out] epsilon Epsilon, Hartree; zero when not found
   !> @param[out] found Whether the type has a GAFF2 row
   pure subroutine lookup_gaff_lj(gaff_type, sigma, epsilon, found)
      !> GAFF atom type
      character(len=*), intent(in) :: gaff_type
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Epsilon, Hartree
      real(wp), intent(out) :: epsilon
      !> Whether the type has a GAFF2 row
      logical, intent(out) :: found
      integer :: idx

      sigma = 0.0_wp
      epsilon = 0.0_wp
      found = .false.
      idx = find_key(normalize_key(gaff_type), gaff_keys)
      if (idx < 1) return
      sigma = 2.0_wp*gaff_rmin_half_aa(idx)/2.0_wp**(1.0_wp/6.0_wp)*aatoau
      epsilon = gaff_epsilon_kcal(idx)*kcaltoau
      found = .true.
   end subroutine lookup_gaff_lj

end module moist_model_moz_potential_lj_gaff_parameters
