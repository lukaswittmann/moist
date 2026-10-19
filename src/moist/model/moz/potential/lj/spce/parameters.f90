! SPDX-License-Identifier: MPL-2.0 AND BSD-3-Clause
! Lookup implementation: MPL-2.0; parameter source: BSD-3-Clause
!
!> Legacy modified SPC/E-style water Lennard-Jones parameters
!>
!> Berendsen et al., J. Phys. Chem. 91, 6269 (1987), doi:10.1021/j100308a038
!> Hydrogen LJ modification: Luchko et al., J. Chem. Theory Comput. 6, 607-624 (2010), doi:10.1021/ct900460m
module moist_model_moz_potential_lj_spce_parameters
   use mctc_env, only: wp
   use mctc_io_convert, only: aatoau, kjtoau
   implicit none(type, external)
   private

   public :: lookup_spce_lj

   ! Permissive reference: VeloxChem src/pymodule/waterparameters.py, cspce entry
   ! https://github.com/VeloxChem/VeloxChem/blob/d4290334938fb5cab7e2e9954bd6dea5865d7b1c/src/pymodule/waterparameters.py
   ! Source revision: d4290334938fb5cab7e2e9954bd6dea5865d7b1c
   ! Source sigma (ow, hw): 0.31658, 0.11658 nm
   ! Source epsilon (ow, hw): 0.649775, 0.064978 kJ/mol
   ! Local legacy values retained, including epsilon_hw = 0.1*epsilon_ow
   ! Nonzero hydrogen LJ terms distinguish this table from original SPC/E
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

   !> Oxygen sigma, nm
   real(wp), parameter :: spce_sigma_ow_nm = 0.316572_wp
   !> Hydrogen sigma, nm
   real(wp), parameter :: spce_sigma_hw_nm = 0.116572_wp
   !> Oxygen epsilon, kJ/mol
   real(wp), parameter :: spce_epsilon_ow_kj = 0.6497752516070757_wp
   !> Hydrogen epsilon, kJ/mol
   real(wp), parameter :: spce_epsilon_hw_kj = 0.06497752516070757_wp

contains

   !> LJ parameters of a water site by element
   !>
   !> @param[in] number Atomic number, 8 for oxygen or 1 for hydrogen
   !> @param[out] sigma Sigma, bohr; zero when not found
   !> @param[out] epsilon Epsilon, Hartree; zero when not found
   !> @param[out] found Whether the element is a water site
   pure subroutine lookup_spce_lj(number, sigma, epsilon, found)
      !> Atomic number
      integer, intent(in) :: number
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Epsilon, Hartree
      real(wp), intent(out) :: epsilon
      !> Whether the element is a water site
      logical, intent(out) :: found

      sigma = 0.0_wp
      epsilon = 0.0_wp
      found = .true.
      select case (number)
      case (8)
         sigma = spce_sigma_ow_nm*10.0_wp*aatoau
         epsilon = spce_epsilon_ow_kj*kjtoau
      case (1)
         sigma = spce_sigma_hw_nm*10.0_wp*aatoau
         epsilon = spce_epsilon_hw_kj*kjtoau
      case default
         found = .false.
      end select
   end subroutine lookup_spce_lj

end module moist_model_moz_potential_lj_spce_parameters
