! SPDX-License-Identifier: MPL-2.0 AND BSD-3-Clause
! Lookup implementation: MPL-2.0; parameter data: BSD-3-Clause
!
!> Transition-metal element Lennard-Jones parameters
!>
!> Sebesta et al., J. Chem. Theory Comput. 12, 3681-3688 (2016), doi:10.1021/acs.jctc.6b00416
module moist_model_moz_potential_lj_tm_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kjtoau
   implicit none(type, external)
   private

   public :: lookup_tm_lj, n_tm

   ! TM parameters from VeloxChem src/pymodule/tmparameters.py
   ! https://github.com/VeloxChem/VeloxChem/blob/d4290334938fb5cab7e2e9954bd6dea5865d7b1c/src/pymodule/tmparameters.py
   ! Source revision: d4290334938fb5cab7e2e9954bd6dea5865d7b1c
   ! Lower oxidation state where the reference supplies multiple states
   ! Source conversion: sigma_nm = 0.1*s*2.0; epsilon_kj = e*4.184
   ! Changes: atomic-number indexing, missing-element sentinels, Fortran lookup
   ! Legacy epsilon offsets retained: Ti +4e-8 and Hg +1e-7 kJ/mol
   !
   ! Copyright 2018-2025 VeloxChem developers
   !
   ! Redistribution and use in source and binary forms, with or without modification,
   ! are permitted provided that the following conditions are met:
   !
   ! 1. Redistributions of source code must retain the above copyright notice, this
   !    list of conditions and the following disclaimer.
   ! 2. Redistributions in binary form must reproduce the above copyright notice,
   !    this list of conditions and the following disclaimer in the documentation
   !    and/or other materials provided with the distribution.
   ! 3. Neither the name of the copyright holder nor the names of its contributors
   !    may be used to endorse or promote products derived from this software without
   !    specific prior written permission.
   !
   ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   ! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   ! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
   ! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   ! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   ! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   ! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
   ! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
   ! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   !> Highest supported atomic number (Hg); the tables are indexed by atomic number
   integer, parameter :: n_tm = 80

   !> Sentinel of an unsupported element
   real(wp), parameter :: miss = -1.0_wp

   !> Sigma per atomic number, nm; `miss` where unsupported
   real(wp), parameter :: tm_sigma_nm(n_tm) = [ &
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 1-10
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 11-20
      & 0.6660_wp, 0.5738_wp, 0.5534_wp, 0.5468_wp, 0.5808_wp, & ! Sc Ti V Cr Mn
      & 0.5938_wp, 0.5610_wp, 0.5440_wp, 0.5336_wp, 0.5694_wp, & ! Fe Co Ni Cu Zn
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 31-40
      & miss, miss, miss, 0.5878_wp, 0.5430_wp, miss, miss, miss, miss, miss, & ! 41-50, Ru Rh
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 51-60
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 61-70
      & miss, miss, miss, miss, miss, miss, miss, 0.5340_wp, miss, 0.5602_wp & ! 71-80, Pt Hg
      & ]

   !> Epsilon per atomic number, kJ/mol; `miss` where unsupported
   real(wp), parameter :: tm_epsilon_kj(n_tm) = [ &
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 1-10
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 11-20
      & 0.46024_wp, 3.83672804_wp, 7.966336_wp, 6.351312_wp, 4.556376_wp, & ! Sc Ti V Cr Mn
      & 2.920432_wp, 5.004064_wp, 11.0876_wp, 8.987232_wp, 4.58148_wp, & ! Fe Co Ni Cu Zn
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 31-40
      & miss, miss, miss, 1.748912_wp, 16.510064_wp, miss, miss, miss, miss, miss, & ! 41-50, Ru Rh
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 51-60
      & miss, miss, miss, miss, miss, miss, miss, miss, miss, miss, & ! 61-70
      & miss, miss, miss, miss, miss, miss, miss, 21.225432_wp, miss, 8.1797201_wp & ! 71-80, Pt Hg
      & ]

contains

   !> Transition-metal Lennard-Jones parameters of an element
   !>
   !> @param[in] number Atomic number
   !> @param[out] sigma Sigma, bohr; zero when unsupported
   !> @param[out] epsilon Well depth, Hartree; zero when unsupported
   !> @param[out] found Whether the element has a row
   pure subroutine lookup_tm_lj(number, sigma, epsilon, found)
      !> Atomic number
      integer, intent(in) :: number
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Well depth, Hartree
      real(wp), intent(out) :: epsilon
      !> Whether the element has a row
      logical, intent(out) :: found

      sigma = 0.0_wp
      epsilon = 0.0_wp
      found = .false.
      if (number < 1 .or. number > n_tm) return
      if (tm_sigma_nm(number) == miss) return
      sigma = tm_sigma_nm(number)*10.0_wp*aatoau
      epsilon = tm_epsilon_kj(number)*kjtoau
      found = .true.
   end subroutine lookup_tm_lj

end module moist_model_moz_potential_lj_tm_parameters
