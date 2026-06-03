module moist_model_gems_utils
   use mctc_env, only: wp
   use mctc_io, only: structure_type, &
      & new_structure
   use mctc_env_error, only: error_type

contains

   subroutine BuildSuperStructure(soluteMol, solventMol, SuperStructure, no_displacement)
      use mctc_env, only: wp
      use mctc_io, only: structure_type
      implicit none

      ! Inputs
      type(structure_type), intent(in)  :: soluteMol
      type(structure_type), intent(in)  :: solventMol

      ! Output
      type(structure_type), intent(out) :: SuperStructure

      ! Displace solvent molecules?
      logical, intent(in), optional :: no_displacement

      integer :: ns, nv, ntot, i
      integer, allocatable :: z(:)
      real(wp), allocatable :: xyz(:, :)
      character(len=2), allocatable :: atom_sym(:)

      ! 1) figure out dimensions
      ns = soluteMol%nat
      nv = solventMol%nat
      ntot = ns + nv

      ! 2) build the atomic-number array for each atom
      allocate (z(ntot))
      do i = 1, ns
         z(i) = soluteMol%num(soluteMol%id(i))
      end do
      do i = 1, nv
         z(ns + i) = solventMol%num(solventMol%id(i))
      end do

      ! 2.1) build the element symbols array
      allocate (atom_sym(ntot))
      do i = 1, ns
         atom_sym(i) = soluteMol%sym(soluteMol%id(i))
      end do
      do i = 1, nv
         atom_sym(ns + i) = solventMol%sym(solventMol%id(i))
      end do

      ! 3) build the coordinates array
      allocate (xyz(3, ntot))
      xyz(:, 1:ns) = soluteMol%xyz(:, 1:ns)
      xyz(:, ns + 1:ntot) = solventMol%xyz(:, 1:nv)

      ! apply a large offset to the solvent block
      if (present(no_displacement)) then
         if (.not. no_displacement) then
            xyz(1, ns + 1:ntot) = xyz(1, ns + 1:ntot)
            xyz(2, ns + 1:ntot) = xyz(2, ns + 1:ntot)
            xyz(3, ns + 1:ntot) = xyz(3, ns + 1:ntot)
         end if
      else
         xyz(1, ns + 1:ntot) = xyz(1, ns + 1:ntot) + 3000.0_wp
         xyz(2, ns + 1:ntot) = xyz(2, ns + 1:ntot) + 3000.0_wp
         xyz(3, ns + 1:ntot) = xyz(3, ns + 1:ntot) + 3000.0_wp
      end if

      call new_structure(SuperStructure, z, atom_sym, xyz)

   end subroutine BuildSuperStructure

end module moist_model_gems_utils
