!> DREIDING element Lennard-Jones parameters for explicit-atom structures
!>
!> The original table gives one nitrogen well depth, 0.0774 kcal/mol
!> Implicit-hydrogen carbon rows and the special H-HB row are excluded
!> The copied Cu/Ni/Mg extensions are excluded; Na/Ca/Fe/Zn are retained
!>
!> DREIDING mixes by Lorentz-Berthelot
module moist_model_moz_potential_lj_dreiding_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kcaltoau
   implicit none(type, external)
   private

   public :: lookup_dreiding_lj, n_dreiding

   ! Mayo, Olafson and Goddard, J. Phys. Chem. 94, 8897 (1990)
   ! doi:10.1021/j100389a010, Table II and equations 36a/36c
   ! Transcribed from the supplied lammps_interface/dreiding.py table
   ! https://github.com/peteboyd/lammps_interface/blob/master/lammps_interface/dreiding.py
   !
   ! The MIT License (MIT)
   !
   ! Copyright (c) 2017 Peter Boyd
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

   !> Highest atomic number in the table, I
   integer, parameter :: n_dreiding = 53

   !> Sentinel of an element without parameters
   real(wp), parameter :: missing = -1.0_wp

   !> Homonuclear LJ minimum distance per atomic number, Angstrom; `missing` when absent
   real(wp), parameter :: dreiding_rmin_aa(n_dreiding) = [ &
      & 3.195_wp, missing, missing, missing, 4.02_wp, 3.8983_wp, 3.6621_wp, 3.4046_wp, &  ! H-O
      & 3.472_wp, missing, 3.144_wp, missing, 4.39_wp, 4.27_wp, 4.15_wp, 4.03_wp, &  ! F-S
      & 3.9503_wp, missing, missing, 3.472_wp, missing, missing, missing, missing, &  ! Cl-Cr
      & missing, 4.54_wp, missing, missing, missing, 4.54_wp, 4.39_wp, 4.27_wp, &  ! Mn-Ge
      & 4.15_wp, 4.03_wp, 3.95_wp, missing, missing, missing, missing, missing, &  ! As-Zr
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Nb-Cd
      & 4.59_wp, 4.47_wp, 4.35_wp, 4.23_wp, 4.15_wp &  ! In-I
      & ]

   !> LJ well depth per atomic number, kcal/mol; `missing` when absent
   real(wp), parameter :: dreiding_epsilon_kcal(n_dreiding) = [ &
      & 0.0152_wp, missing, missing, missing, 0.095_wp, 0.0951_wp, 0.0774_wp, 0.0957_wp, &  ! H-O
      & 0.0725_wp, missing, 0.5_wp, missing, 0.31_wp, 0.31_wp, 0.32_wp, 0.344_wp, &  ! F-S
      & 0.2833_wp, missing, missing, 0.05_wp, missing, missing, missing, missing, &  ! Cl-Cr
      & missing, 0.055_wp, missing, missing, missing, 0.055_wp, 0.4_wp, 0.4_wp, &  ! Mn-Ge
      & 0.41_wp, 0.43_wp, 0.37_wp, missing, missing, missing, missing, missing, &  ! As-Zr
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Nb-Cd
      & 0.55_wp, 0.55_wp, 0.55_wp, 0.57_wp, 0.51_wp &  ! In-I
      & ]

contains

   !> Look up element LJ parameters in atomic units
   !>
   !> @param[in] number Atomic number
   !> @param[out] sigma Sigma, bohr; zero when absent
   !> @param[out] epsilon Well depth, Hartree; zero when absent
   !> @param[out] found Whether the element is supported
   pure subroutine lookup_dreiding_lj(number, sigma, epsilon, found)
      !> Atomic number
      integer, intent(in) :: number
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Well depth, Hartree
      real(wp), intent(out) :: epsilon
      !> Whether the element is supported
      logical, intent(out) :: found

      sigma = 0.0_wp
      epsilon = 0.0_wp
      found = .false.
      if (number < 1 .or. number > n_dreiding) return
      if (dreiding_rmin_aa(number) < 0.0_wp) return
      sigma = dreiding_rmin_aa(number)/2.0_wp**(1.0_wp/6.0_wp)*aatoau
      epsilon = dreiding_epsilon_kcal(number)*kcaltoau
      found = .true.
   end subroutine lookup_dreiding_lj

end module moist_model_moz_potential_lj_dreiding_parameters
