!> Tests of the OPLS-AA Lennard-Jones typing: moz_potential_oplsaa
!>
!> - Methanol typing from explicit bond graph
!> - Elements without OPLS-AA rules left for later terms
!> - Impossible valence refusal
module test_moz_potential_oplsaa
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check
   use moist_model_moz_potential_lj_typed, only: lj_typed_type, lj_oplsaa
   use moist_model_moz_potential_lj_element, only: lj_uff
   use moist_model_moz_potential, only: moz_potential_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use test_helpers, only: check_moist_error
   implicit none(type, external)
   private

   public :: collect_moz_potential_oplsaa

contains

   !> Register tests
   subroutine collect_moz_potential_oplsaa(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("methanol-types", check_methanol_types), &
         & new_unittest("argon-alone-uncovered", check_argon_alone), &
         & new_unittest("argon-beside-methanol", check_argon_beside_methanol), &
         & new_unittest("invalid-valence", check_invalid_valence)]
   end subroutine collect_moz_potential_oplsaa

   !> Build methanol with explicit single bonds and optional unbonded argon
   !>
   !> @param[out] mol    structure: C, O, three H on C, H on O (and Ar)
   !> @param[in]  argon  whether to append an argon atom
   subroutine make_methanol(mol, argon)
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Whether to append an argon atom
      logical, intent(in) :: argon
      real(wp), parameter :: xyz(3, 7) = reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 2.70_wp, 0.0_wp, 0.0_wp, &
         & -0.69_wp, 1.94_wp, 0.0_wp, -0.69_wp, -0.97_wp, 1.68_wp, -0.69_wp, -0.97_wp, -1.68_wp, &
         & 3.30_wp, 1.75_wp, 0.0_wp, 12.0_wp, 0.0_wp, 0.0_wp], [3, 7])
      if (argon) then
         call new(mol, [6, 8, 1, 1, 1, 1, 18], xyz)
      else
         call new(mol, [6, 8, 1, 1, 1, 1], xyz(:, :6))
      end if
      mol%nbd = 5
      mol%bond = reshape([1, 2, 1, 1, 3, 1, 1, 4, 1, 1, 5, 1, 2, 6, 1], [3, 5])
   end subroutine make_methanol

   !> Check methanol OPLS-AA labels for carbon, oxygen and both hydrogen kinds
   subroutine check_methanol_types(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      character(len=8), parameter :: expected(6) = [character(len=8) :: &
         & "opls_157", "opls_154", "opls_156", "opls_156", "opls_156", "opls_155"]
      integer :: iat

      call make_methanol(mol, .false.)
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(oplsaa%pair%covered), "OPLS-AA must cover every methanol atom")
      if (allocated(error)) return
      do iat = 1, mol%nat
         call check(error, trim(oplsaa%atomtype(iat)), trim(expected(iat)), "methanol OPLS-AA type")
         if (allocated(error)) return
      end do
   end subroutine check_methanol_types

   !> Check UFF fallback for isolated argon without OPLS-AA coverage
   subroutine check_argon_alone(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moz_potential_type) :: pot
      type(potential_sites_type) :: sites

      call new(mol, [18], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
      mol%nbd = 0
      allocate (mol%bond(3, 0))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(oplsaa%pair%covered), "OPLS-AA must leave argon uncovered")
      if (allocated(error)) return

      call pot%add(lj_oplsaa, err)
      if (.not. allocated(err)) call pot%add(lj_uff, err)
      if (.not. allocated(err)) call pot%update(mol, err)
      if (.not. allocated(err)) call pot%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, sites%pair%label(1), "Ar/UFF", "UFF must fill the argon OPLS-AA left")
   end subroutine check_argon_alone

   !> Check methanol OPLS-AA typing and UFF fallback restricted to argon
   subroutine check_argon_beside_methanol(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(moz_potential_type) :: pot
      type(potential_sites_type) :: sites

      call make_methanol(mol, .true.)
      call pot%add(lj_oplsaa, err)
      if (.not. allocated(err)) call pot%add(lj_uff, err)
      if (.not. allocated(err)) call pot%update(mol, err)
      if (.not. allocated(err)) call pot%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sites%pair%covered), "OPLS-AA then UFF must cover every atom")
      if (allocated(error)) return
      call check(error, sites%pair%label(1), "opls_157", "OPLS-AA keeps the methanol carbon")
      if (allocated(error)) return
      call check(error, sites%pair%label(7), "Ar/UFF", "UFF fills argon")
   end subroutine check_argon_beside_methanol

   !> Check graph refusal for carbon with five hydrogens
   subroutine check_invalid_valence(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      real(wp), parameter :: xyz(3, 6) = reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 2.0_wp], [3, 6])

      call new(mol, [6, 1, 1, 1, 1, 1], xyz)
      mol%nbd = 5
      mol%bond = reshape([1, 2, 1, 1, 3, 1, 1, 4, 1, 1, 5, 1, 1, 6, 1], [3, 5])
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call check(error, allocated(err), "a pentavalent carbon was typed")
      if (allocated(error)) return
      call check(error, index(err%message, "Unsupported OPLS-AA valence at atom 1") > 0, err%message)
   end subroutine check_invalid_valence

end module test_moz_potential_oplsaa
