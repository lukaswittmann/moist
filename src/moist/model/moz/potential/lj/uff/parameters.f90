! SPDX-License-Identifier: MPL-2.0 AND MIT
! Lookup implementation: MPL-2.0; parameter data: MIT
!
!> UFF element Lennard-Jones parameters
!>
!> Rappe et al., J. Am. Chem. Soc. 114, 10024-10035 (1992), doi:10.1021/ja00051a040
module moist_model_moz_potential_lj_uff_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kcaltoau
   implicit none(type, external)
   private

   public :: lookup_uff_lj, n_uff

   ! UFF x1 and D1 columns from lammps_interface/uff.py
   ! https://github.com/peteboyd/lammps_interface/blob/255f027cb76142d39c050a6810404debc6a06562/lammps_interface/uff.py
   ! Source revision: 255f027cb76142d39c050a6810404debc6a06562
   ! Changes: collapse identical element LJ rows; omit Du; map Lw6+3 to Lr
   ! x1 = homonuclear rmin, Angstrom; D1 = epsilon, kcal/mol
   ! Sigma = x1/2**(1/6); conversion to atomic units on lookup
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

   !> Number of supported elements, H through Lr; the tables are indexed by atomic number
   integer, parameter :: n_uff = 103

   !> Homonuclear LJ minimum distance per atomic number, Angstrom
   real(wp), parameter :: uff_rmin_aa(n_uff) = [ &
      & 2.886_wp, 2.362_wp, 2.451_wp, 2.745_wp, &
      & 4.083_wp, 3.851_wp, 3.66_wp, 3.5_wp, &
      & 3.364_wp, 3.243_wp, 2.983_wp, 3.021_wp, &
      & 4.499_wp, 4.295_wp, 4.147_wp, 4.035_wp, &
      & 3.947_wp, 3.868_wp, 3.812_wp, 3.399_wp, &
      & 3.295_wp, 3.175_wp, 3.144_wp, 3.023_wp, &
      & 2.961_wp, 2.912_wp, 2.872_wp, 2.834_wp, &
      & 3.495_wp, 2.763_wp, 4.383_wp, 4.28_wp, &
      & 4.23_wp, 4.205_wp, 4.189_wp, 4.141_wp, &
      & 4.114_wp, 3.641_wp, 3.345_wp, 3.124_wp, &
      & 3.165_wp, 3.052_wp, 2.998_wp, 2.963_wp, &
      & 2.929_wp, 2.899_wp, 3.148_wp, 2.848_wp, &
      & 4.463_wp, 4.392_wp, 4.42_wp, 4.47_wp, &
      & 4.5_wp, 4.404_wp, 4.517_wp, 3.703_wp, &
      & 3.522_wp, 3.556_wp, 3.606_wp, 3.575_wp, &
      & 3.547_wp, 3.52_wp, 3.493_wp, 3.368_wp, &
      & 3.451_wp, 3.428_wp, 3.409_wp, 3.391_wp, &
      & 3.374_wp, 3.355_wp, 3.64_wp, 3.141_wp, &
      & 3.17_wp, 3.069_wp, 2.954_wp, 3.12_wp, &
      & 2.84_wp, 2.754_wp, 3.293_wp, 2.705_wp, &
      & 4.347_wp, 4.297_wp, 4.37_wp, 4.709_wp, &
      & 4.75_wp, 4.765_wp, 4.9_wp, 3.677_wp, &
      & 3.478_wp, 3.396_wp, 3.424_wp, 3.395_wp, &
      & 3.424_wp, 3.424_wp, 3.381_wp, 3.326_wp, &
      & 3.339_wp, 3.313_wp, 3.299_wp, 3.286_wp, &
      & 3.274_wp, 3.248_wp, 3.236_wp &
      & ]

   !> LJ well depth per atomic number, kcal/mol
   real(wp), parameter :: uff_epsilon_kcal(n_uff) = [ &
      & 0.044_wp, 0.056_wp, 0.025_wp, 0.085_wp, &
      & 0.18_wp, 0.105_wp, 0.069_wp, 0.06_wp, &
      & 0.05_wp, 0.042_wp, 0.03_wp, 0.111_wp, &
      & 0.505_wp, 0.402_wp, 0.305_wp, 0.274_wp, &
      & 0.227_wp, 0.185_wp, 0.035_wp, 0.238_wp, &
      & 0.019_wp, 0.017_wp, 0.016_wp, 0.015_wp, &
      & 0.013_wp, 0.013_wp, 0.014_wp, 0.015_wp, &
      & 0.005_wp, 0.124_wp, 0.415_wp, 0.379_wp, &
      & 0.309_wp, 0.291_wp, 0.251_wp, 0.22_wp, &
      & 0.04_wp, 0.235_wp, 0.072_wp, 0.069_wp, &
      & 0.059_wp, 0.056_wp, 0.048_wp, 0.056_wp, &
      & 0.053_wp, 0.048_wp, 0.036_wp, 0.228_wp, &
      & 0.599_wp, 0.567_wp, 0.449_wp, 0.398_wp, &
      & 0.339_wp, 0.332_wp, 0.045_wp, 0.364_wp, &
      & 0.017_wp, 0.013_wp, 0.01_wp, 0.01_wp, &
      & 0.009_wp, 0.008_wp, 0.008_wp, 0.009_wp, &
      & 0.007_wp, 0.007_wp, 0.007_wp, 0.007_wp, &
      & 0.006_wp, 0.228_wp, 0.041_wp, 0.072_wp, &
      & 0.081_wp, 0.067_wp, 0.066_wp, 0.037_wp, &
      & 0.073_wp, 0.08_wp, 0.039_wp, 0.385_wp, &
      & 0.68_wp, 0.663_wp, 0.518_wp, 0.325_wp, &
      & 0.284_wp, 0.248_wp, 0.05_wp, 0.404_wp, &
      & 0.033_wp, 0.026_wp, 0.022_wp, 0.022_wp, &
      & 0.019_wp, 0.016_wp, 0.014_wp, 0.013_wp, &
      & 0.013_wp, 0.013_wp, 0.012_wp, 0.012_wp, &
      & 0.011_wp, 0.011_wp, 0.011_wp &
      & ]

contains

   !> Look up element LJ parameters in atomic units
   !>
   !> @param[in] number Atomic number
   !> @param[out] sigma Sigma, bohr; zero when absent
   !> @param[out] epsilon Well depth, Hartree; zero when absent
   !> @param[out] found Whether the element is supported
   pure subroutine lookup_uff_lj(number, sigma, epsilon, found)
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
      found = number >= 1 .and. number <= n_uff
      if (.not. found) return
      sigma = uff_rmin_aa(number)/2.0_wp**(1.0_wp/6.0_wp)*aatoau
      epsilon = uff_epsilon_kcal(number)*kcaltoau
   end subroutine lookup_uff_lj

end module moist_model_moz_potential_lj_uff_parameters
