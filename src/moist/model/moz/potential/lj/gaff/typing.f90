!> GAFF atom typing from element and geometry
!>
!> Ported unchanged in logic from the moist_dev_moz RISM potential layer:
!> bonds from 1.25 times the summed covalent radii, ring membership by an
!> alternative-path search, aromatic C/N as ring atoms of degree 3, heavy-atom
!> types by element and degree, hydrogen types by their bonded partner, and an
!> approximate cp/cq alternation on aromatic C-C bonds. Elements without a
!> GAFF rule type as 'du'
module moist_model_moz_potential_lj_gaff_typing
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use moist_data_atomicrad, only: covalent_rad
   implicit none(type, external)
   private

   public :: generate_gaff_atomtypes

   ! GAFF typing taken from VeloxChem (BSD-3)
   !
   ! Copyright 2018-2025 VeloxChem developers
   !
   ! Redistribution and use in source and binary forms, with or without modification,
   ! are permitted provided that the following conditions are met:
   !
   ! 1. Redistributions of source code must retain the above copyright notice, this
   !    list of conditions and the following disclaimer.
   !
   ! 2. Redistributions in binary form must reproduce the above copyright notice,
   !    this list of conditions and the following disclaimer in the documentation
   !    and/or other materials provided with the distribution.
   !
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


   !> Working state of the typing
   type :: gaff_typing_state_type
      !> Number of atoms
      integer :: n_atoms = 0
      !> Atomic number per atom (n_atoms)
      integer, allocatable :: num(:)
      !> Interatomic distances, bohr (n_atoms, n_atoms)
      real(wp), allocatable :: dist(:, :)
      !> Bond matrix, 0/1 (n_atoms, n_atoms)
      integer, allocatable :: conn(:, :)
      !> Number of bonds per atom (n_atoms)
      integer, allocatable :: degree(:)
      !> Ring membership per atom (n_atoms)
      logical, allocatable :: cyclic(:)
      !> Aromaticity per atom (n_atoms)
      logical, allocatable :: aromatic(:)
      !> GAFF type per atom (n_atoms)
      character(len=2), allocatable :: gaff(:)
   end type gaff_typing_state_type

contains

   !> GAFF atom type of every atom of a structure
   !>
   !> @param[in] mol Structure, coordinates in bohr
   !> @param[out] gaff_atom_types Two-character GAFF type per atom (nat)
   subroutine generate_gaff_atomtypes(mol, gaff_atom_types)
      !> Structure, coordinates in bohr
      class(structure_type), intent(in) :: mol
      !> Two-character GAFF type per atom
      character(len=2), allocatable, intent(out) :: gaff_atom_types(:)
      type(gaff_typing_state_type) :: ati

      call init_state(ati, mol)
      call compute_distance_matrix(mol%xyz, ati%dist)
      call build_connectivity_from_geometry(ati%num, ati%dist, ati%conn)
      call detect_closed_cyclic_structures(ati)
      call update_degrees(ati)
      call decide_atom_type(ati)
      call sanitize_gaff_atom_types(ati)
      call check_alternating_atom_types(ati)

      gaff_atom_types = ati%gaff
   end subroutine generate_gaff_atomtypes

   !> Allocate and reset the working state
   !>
   !> @param[in,out] ati Typing state
   !> @param[in] mol Structure
   subroutine init_state(ati, mol)
      !> Typing state
      type(gaff_typing_state_type), intent(inout) :: ati
      !> Structure
      class(structure_type), intent(in) :: mol
      integer :: i

      ati%n_atoms = mol%nat

      allocate (ati%num(ati%n_atoms))
      allocate (ati%dist(ati%n_atoms, ati%n_atoms))
      allocate (ati%conn(ati%n_atoms, ati%n_atoms))
      allocate (ati%degree(ati%n_atoms))
      allocate (ati%cyclic(ati%n_atoms))
      allocate (ati%aromatic(ati%n_atoms))
      allocate (ati%gaff(ati%n_atoms))

      do i = 1, ati%n_atoms
         ati%num(i) = mol%num(mol%id(i))
      end do
      ati%conn = 0
      ati%degree = 0
      ati%cyclic = .false.
      ati%aromatic = .false.
      ati%gaff = "du"
   end subroutine init_state

   !> Interatomic distance matrix
   !>
   !> @param[in] xyz Coordinates, bohr (3, n)
   !> @param[out] dist Distances, bohr (n, n)
   pure subroutine compute_distance_matrix(xyz, dist)
      !> Coordinates, bohr
      real(wp), intent(in) :: xyz(:, :)
      !> Distances, bohr
      real(wp), intent(out) :: dist(:, :)
      integer :: i, j, n
      real(wp) :: dx, dy, dz

      n = size(xyz, 2)
      dist = 0.0_wp
      do i = 1, n
         do j = i + 1, n
            dx = xyz(1, i) - xyz(1, j)
            dy = xyz(2, i) - xyz(2, j)
            dz = xyz(3, i) - xyz(3, j)
            dist(i, j) = sqrt(dx*dx + dy*dy + dz*dz)
            dist(j, i) = dist(i, j)
         end do
      end do
   end subroutine compute_distance_matrix

   !> Bonds from 1.25 times the summed covalent radii
   !>
   !> @param[in] num Atomic number per atom (n)
   !> @param[in] dist Distances, bohr (n, n)
   !> @param[out] conn Bond matrix, 0/1 (n, n)
   pure subroutine build_connectivity_from_geometry(num, dist, conn)
      !> Atomic number per atom
      integer, intent(in) :: num(:)
      !> Distances, bohr
      real(wp), intent(in) :: dist(:, :)
      !> Bond matrix, 0/1
      integer, intent(out) :: conn(:, :)
      integer :: i, j, n
      real(wp) :: ri, rj, cutoff

      n = size(num)
      conn = 0

      do i = 1, n
         ri = covalent_rad(num(i))
         do j = i + 1, n
            rj = covalent_rad(num(j))
            cutoff = 1.25_wp*(ri + rj)
            if (dist(i, j) > 0.0_wp .and. dist(i, j) <= cutoff) then
               conn(i, j) = 1
               conn(j, i) = 1
            end if
         end do
      end do
   end subroutine build_connectivity_from_geometry

   !> Degrees, ring membership and aromaticity
   !>
   !> A bond is in a ring when its ends stay connected without it; aromatic
   !> atoms are ring C or N of degree 3
   !>
   !> @param[in,out] ati Typing state
   subroutine detect_closed_cyclic_structures(ati)
      !> Typing state
      type(gaff_typing_state_type), intent(inout) :: ati
      integer :: i, j, n

      n = ati%n_atoms
      ati%degree = 0
      do i = 1, n
         ati%degree(i) = sum(ati%conn(i, :))
      end do

      ati%cyclic = .false.
      do i = 1, n
         if (ati%degree(i) >= 2) then
            do j = i + 1, n
               if (ati%conn(i, j) == 1) then
                  if (has_alternative_path(ati%conn, i, j)) then
                     ati%cyclic(i) = .true.
                     ati%cyclic(j) = .true.
                  end if
               end if
            end do
         end if
      end do

      ati%aromatic = .false.
      do i = 1, n
         if (ati%cyclic(i) .and. ati%degree(i) == 3) then
            if (ati%num(i) == 6 .or. ati%num(i) == 7) then
               ati%aromatic(i) = .true.
            end if
         end if
      end do
   end subroutine detect_closed_cyclic_structures

   !> Whether `a` and `b` are connected without their direct bond
   !>
   !> @param[in] conn Bond matrix, 0/1 (n, n)
   !> @param[in] a First atom
   !> @param[in] b Second atom
   pure function has_alternative_path(conn, a, b) result(found)
      !> Bond matrix, 0/1
      integer, intent(in) :: conn(:, :)
      !> First atom
      integer, intent(in) :: a
      !> Second atom
      integer, intent(in) :: b
      !> Whether another path exists
      logical :: found
      integer :: n, head, tail, cur, k
      logical, allocatable :: visited(:)
      integer, allocatable :: queue(:)

      n = size(conn, 1)
      allocate (visited(n), queue(n))
      visited = .false.
      queue = 0
      head = 1
      tail = 1
      queue(1) = a
      visited(a) = .true.
      found = .false.

      do while (head <= tail)
         cur = queue(head)
         head = head + 1
         do k = 1, n
            if (cur == a .and. k == b) cycle
            if (cur == b .and. k == a) cycle
            if (conn(cur, k) == 1 .and. .not. visited(k)) then
               visited(k) = .true.
               tail = tail + 1
               queue(tail) = k
               if (k == b) then
                  found = .true.
                  return
               end if
            end if
         end do
      end do
   end function has_alternative_path

   !> Recompute the degrees from the bond matrix
   !>
   !> @param[in,out] ati Typing state
   pure subroutine update_degrees(ati)
      !> Typing state
      type(gaff_typing_state_type), intent(inout) :: ati
      integer :: i

      do i = 1, ati%n_atoms
         ati%degree(i) = sum(ati%conn(i, :))
      end do
   end subroutine update_degrees

   !> Assign the GAFF types
   !>
   !> Heavy atoms by element, degree and neighbours; hydrogens by their first
   !> bonded partner; N bonded to sp2/aromatic C (and no N) becomes 'na'
   !>
   !> @param[in,out] ati Typing state
   subroutine decide_atom_type(ati)
      !> Typing state
      type(gaff_typing_state_type), intent(inout) :: ati
      integer :: i, j, n, d
      logical :: has_h, has_heavy, has_o, has_n, has_c_sp2
      character(len=2) :: sj

      n = ati%n_atoms

      ! Heavy atoms
      do i = 1, n
         d = ati%degree(i)
         select case (ati%num(i))
         case (6)
            if (ati%aromatic(i)) then
               ati%gaff(i) = "ca"
            else if (d >= 4) then
               ati%gaff(i) = "c3"
            else if (d == 3) then
               ati%gaff(i) = "c2"
            else if (d == 2) then
               ati%gaff(i) = "c1"
            else
               ati%gaff(i) = "c "
            end if

         case (7)
            if (d >= 4) then
               ati%gaff(i) = "n+"
            else if (ati%aromatic(i)) then
               ati%gaff(i) = "na"
            else if (d == 3) then
               ati%gaff(i) = "n3"
            else if (d == 2) then
               ati%gaff(i) = "n2"
            else
               ati%gaff(i) = "n "
            end if

         case (8)
            has_h = .false.
            has_heavy = .false.
            do j = 1, n
               if (ati%conn(i, j) == 1) then
                  if (ati%num(j) == 1) then
                     has_h = .true.
                  else
                     has_heavy = .true.
                  end if
               end if
            end do
            if (d == 1) then
               ati%gaff(i) = "o "
            else if (has_h .and. .not. has_heavy) then
               ! Water oxygen: every neighbour is H
               ati%gaff(i) = "ow"
            else if (has_h) then
               ! Hydroxyl oxygen: at least one heavy-atom neighbour
               ati%gaff(i) = "oh"
            else if (d == 2) then
               ati%gaff(i) = "os"
            else
               ati%gaff(i) = "ow"
            end if

         case (16)
            has_h = .false.
            do j = 1, n
               if (ati%conn(i, j) == 1 .and. ati%num(j) == 1) has_h = .true.
            end do
            if (d <= 1) then
               ati%gaff(i) = "s "
            else if (d == 2 .and. has_h) then
               ati%gaff(i) = "sh"
            else if (d == 2) then
               ati%gaff(i) = "ss"
            else if (d == 3) then
               ati%gaff(i) = "s4"
            else
               ati%gaff(i) = "s6"
            end if

         case (15)
            has_o = .false.
            do j = 1, n
               if (ati%conn(i, j) == 1 .and. ati%num(j) == 8) has_o = .true.
            end do
            if (d == 2) then
               ati%gaff(i) = "p2"
            else if (d == 3 .and. .not. has_o) then
               ati%gaff(i) = "p3"
            else if (d == 3 .and. has_o) then
               ati%gaff(i) = "p4"
            else
               ati%gaff(i) = "p5"
            end if

         case (9)
            ati%gaff(i) = "f "
         case (17)
            ati%gaff(i) = "cl"
         case (35)
            ati%gaff(i) = "br"
         case (53)
            ati%gaff(i) = "i "
         case (1)
            ati%gaff(i) = "h1"
         case default
            ati%gaff(i) = "du"
         end select
      end do

      ! Hydrogens by their first bonded partner
      do i = 1, n
         if (ati%num(i) /= 1) cycle
         do j = 1, n
            if (ati%conn(i, j) /= 1) cycle
            sj = ati%gaff(j)
            select case (trim(sj))
            case ("oh")
               ati%gaff(i) = "ho"
            case ("ow")
               ati%gaff(i) = "hw"
            case ("n", "n2", "n3", "na", "n+")
               ati%gaff(i) = "hn"
            case ("p2", "p3", "p4", "p5")
               ati%gaff(i) = "hp"
            case ("sh")
               ati%gaff(i) = "hs"
            case ("ca", "c2")
               ati%gaff(i) = "ha"
            case ("c1")
               ati%gaff(i) = "h4"
            case default
               ati%gaff(i) = "hc"
            end select
            exit
         end do
      end do

      ! N bonded to sp2/aromatic C
      do i = 1, n
         if (ati%num(i) /= 7) cycle
         if (trim(ati%gaff(i)) == "n3" .or. trim(ati%gaff(i)) == "n2") then
            has_c_sp2 = .false.
            has_n = .false.
            do j = 1, n
               if (ati%conn(i, j) /= 1) cycle
               if (trim(ati%gaff(j)) == "ca" .or. trim(ati%gaff(j)) == "c2") has_c_sp2 = .true.
               if (ati%num(j) == 7) has_n = .true.
            end do
            if (has_c_sp2 .and. .not. has_n) ati%gaff(i) = "na"
         end if
      end do
   end subroutine decide_atom_type

   !> Normalize the assigned keys to two left-adjusted characters
   !>
   !> @param[in,out] ati Typing state
   pure subroutine sanitize_gaff_atom_types(ati)
      !> Typing state
      type(gaff_typing_state_type), intent(inout) :: ati
      integer :: i

      do i = 1, ati%n_atoms
         if (len_trim(ati%gaff(i)) == 0) ati%gaff(i) = "du"
         ati%gaff(i) = adjustl(ati%gaff(i))
      end do
   end subroutine sanitize_gaff_atom_types

   !> Approximate cp/cq alternation on aromatic C-C bonds
   !>
   !> @param[in,out] ati Typing state
   pure subroutine check_alternating_atom_types(ati)
      !> Typing state
      type(gaff_typing_state_type), intent(inout) :: ati
      integer :: i, j

      do i = 1, ati%n_atoms
         if (trim(ati%gaff(i)) /= "ca") cycle
         do j = i + 1, ati%n_atoms
            if (ati%conn(i, j) == 1 .and. trim(ati%gaff(j)) == "ca") then
               ati%gaff(i) = "cp"
               ati%gaff(j) = "cq"
            end if
         end do
      end do
   end subroutine check_alternating_atom_types

end module moist_model_moz_potential_lj_gaff_typing
