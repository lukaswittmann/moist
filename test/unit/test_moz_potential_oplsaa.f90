!> Tests of the OPLS-AA Lennard-Jones typing: moz_potential_oplsaa
!>
!> - Methanol typing from explicit bond graph
!> - Elements without OPLS-AA rules left for later terms
!> - Impossible valence refusal
!> - Aromatic, carbonyl, nitrile, amine, thiol, ether, sulfur-oxide and charged species
!> - Error paths (`should_fail`): each malformed bond graph and each unsupported valence is refused
!>   through its own named error
module test_moz_potential_oplsaa
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
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
         & new_unittest("methanol_types", check_methanol_types), &
         & new_unittest("argon_alone_uncovered", check_argon_alone), &
         & new_unittest("argon_beside_methanol", check_argon_beside_methanol), &
         & new_unittest("invalid_valence", check_invalid_valence), &
         & new_unittest("aromatic_rings", check_aromatic_rings), &
         & new_unittest("small_molecules", check_small_molecules), &
         & new_unittest("charged_species", check_charged_species), &
         & new_unittest("bad_no_bonds", check_bad_no_bonds, should_fail=.true.), &
         & new_unittest("bad_bond_shape", check_bad_bond_shape, should_fail=.true.), &
         & new_unittest("bad_missing_hydrogen", check_bad_missing_hydrogen, should_fail=.true.), &
         & new_unittest("bad_open_shell", check_bad_open_shell, should_fail=.true.), &
         & new_unittest("bad_charge_annotations", check_bad_charge_annotations, should_fail=.true.), &
         & new_unittest("bad_total_charge", check_bad_total_charge, should_fail=.true.), &
         & new_unittest("bad_self_bond", check_bad_self_bond, should_fail=.true.), &
         & new_unittest("bad_bond_endpoint", check_bad_bond_endpoint, should_fail=.true.), &
         & new_unittest("bad_bond_order_high", check_bad_bond_order_high, should_fail=.true.), &
         & new_unittest("bad_bond_order_zero", check_bad_bond_order_zero, should_fail=.true.), &
         & new_unittest("bad_duplicate_bond", check_bad_duplicate_bond, should_fail=.true.), &
         & new_unittest("bad_species_table", check_bad_species_table, should_fail=.true.), &
         & new_unittest("bad_identity_count", check_bad_identity_count, should_fail=.true.), &
         & new_unittest("bad_empty_structure", check_bad_empty_structure, should_fail=.true.), &
         & new_unittest("bad_valence_single_aromatic", check_bad_valence_single_aromatic, should_fail=.true.), &
         & new_unittest("bad_valence_aromatic_double", check_bad_valence_aromatic_double, should_fail=.true.), &
         & new_unittest("bad_valence_aromatic_fluorine", check_bad_valence_aromatic_fluorine, should_fail=.true.), &
         & new_unittest("bad_valence_aromatic_sulfur", check_bad_valence_aromatic_sulfur, should_fail=.true.), &
         & new_unittest("bad_valence_ammonium", check_bad_valence_ammonium, should_fail=.true.), &
         & new_unittest("bad_valence_phosphorus", check_bad_valence_phosphorus, should_fail=.true.), &
         & new_unittest("bad_valence_sulfur", check_bad_valence_sulfur, should_fail=.true.), &
         & new_unittest("bad_valence_chlorine", check_bad_valence_chlorine, should_fail=.true.), &
         & new_unittest("bad_valence_bridging_hydrogen", check_bad_valence_bridging_hydrogen, should_fail=.true.), &
         & new_unittest("bad_valence_neutral_lithium", check_bad_valence_neutral_lithium, should_fail=.true.)]
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


   !> Structure from atomic numbers and explicit bonds
   !>
   !> Coordinates are placeholders on a line; OPLS-AA typing reads only the explicit graph
   !>
   !> @param[out] mol     structure
   !> @param[in]  num     atomic numbers
   !> @param[in]  bonds   atom, atom, order per bond (3, nbd)
   !> @param[in]  charge  optional formal charge per atom
   subroutine make_graph(mol, num, bonds, charge)
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Atomic numbers
      integer, intent(in) :: num(:)
      !> Bonds
      integer, intent(in) :: bonds(:, :)
      !> Formal charges
      integer, intent(in), optional :: charge(:)
      real(wp) :: xyz(3, size(num))
      integer :: i
      xyz = 0.0_wp
      do i = 1, size(num)
         xyz(1, i) = 2.5_wp*real(i, wp)
      end do
      call new(mol, num, xyz)
      mol%nbd = size(bonds, 2)
      mol%bond = bonds
      if (present(charge)) then
         allocate (mol%sdf(size(num)))
         mol%sdf%charge = charge
         mol%charge = real(sum(charge), wp)
      end if
   end subroutine make_graph

   !> Type a graph and compare the rule labels
   !>
   !> @param[out] error     test error
   !> @param[in]  mol       structure with explicit graph
   !> @param[in]  expected  expected rule labels; "*" skips an atom
   subroutine check_opls_labels(error, mol, expected)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Expected labels
      character(len=*), intent(in) :: expected(:)
      type(moist_error), allocatable :: err
      type(lj_typed_type) :: oplsaa
      character(len=64) :: label
      integer :: iat
      call check(error, size(expected), mol%nat, "reference covers a different number of atoms")
      if (allocated(error)) return
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do iat = 1, mol%nat
         if (trim(expected(iat)) == "*") cycle
         write (label, "(a, i0)") "OPLS-AA type of atom ", iat
         call check(error, trim(oplsaa%atomtype(iat)), trim(expected(iat)), trim(label))
         if (allocated(error)) return
      end do
   end subroutine check_opls_labels

   !> Check ring perception through aromatic and Kekule bonds and aromatic heteroatoms
   !>
   !> - Expected types from the OPLS-AA parameter set (Jorgensen et al., JACS 118, 11225, 1996): benzene CA/HA
   !>   145/146, pyridine N/C2/C3/C4 520-523 with H 524-526, furan O/C2/C3 566-568 with H 569/570
   !> - Foyer oplsaa.xml supplies the environment names; the type numbers are the published OPLS-AA ones
   subroutine check_aromatic_rings(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      integer :: bonds(3, 12), pyridine(3, 11), furan(3, 9), i
      do i = 1, 6
         bonds(:, i) = [i, mod(i, 6) + 1, 4]
         bonds(:, i + 6) = [i, i + 6, 1]
      end do
      call make_graph(mol, [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1], bonds)
      call check_opls_labels(error, mol, [character(len=8) :: "opls_145", "opls_145", "opls_145", "opls_145", "opls_145", &
         & "opls_145", "opls_146", "opls_146", "opls_146", "opls_146", "opls_146", "opls_146"])
      if (allocated(error)) return
      ! Kekule bonds of the same ring
      do i = 1, 6
         bonds(3, i) = 1 + mod(i, 2)
      end do
      call make_graph(mol, [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1], bonds)
      call check_opls_labels(error, mol, [character(len=8) :: "opls_145", "opls_145", "opls_145", "opls_145", "opls_145", &
         & "opls_145", "opls_146", "opls_146", "opls_146", "opls_146", "opls_146", "opls_146"])
      if (allocated(error)) return
      ! Pyridine, nitrogen first
      do i = 1, 6
         pyridine(:, i) = [i, mod(i, 6) + 1, 4]
      end do
      do i = 2, 6
         pyridine(:, i + 5) = [i, i + 5, 1]
      end do
      call make_graph(mol, [7, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1], pyridine)
      call check_opls_labels(error, mol, [character(len=8) :: "opls_520", "opls_521", "opls_522", "opls_523", "opls_522", &
         & "opls_521", "opls_524", "opls_525", "opls_526", "opls_525", "opls_524"])
      if (allocated(error)) return
      ! Furan, oxygen first
      do i = 1, 5
         furan(:, i) = [i, mod(i, 5) + 1, 4]
      end do
      do i = 2, 5
         furan(:, i + 4) = [i, i + 4, 1]
      end do
      call make_graph(mol, [8, 6, 6, 6, 6, 1, 1, 1, 1], furan)
      call check_opls_labels(error, mol, [character(len=8) :: "opls_566", "opls_567", "opls_568", "opls_568", "opls_567", &
         & "opls_569", "opls_570", "opls_570", "opls_569"])
   end subroutine check_aromatic_rings

   !> Check OPLS-AA types of small polar species
   !>
   !> - Expected types from the OPLS-AA parameter set: acetonitrile CH3/C/N/H 755/754/753/759 (Jorgensen and
   !>   Tirado-Rives 1988), DMF O/C/N/CH3/formyl H 236/235/239/243/279 (Jorgensen and Swenson 1985), dimethyl ether
   !>   O/CH3/H 180/181/185, methanethiol S/H/CH3 200/204/217 (CH2 of RCH2SH is 206), methylamine N/CH3/HN/HC 900/903/909/911,
   !>   DMSO S/O/CH3 496/497/498, sulfolane S/O 493/494
   !> - Atoms whose type the references do not fix are skipped
   subroutine check_small_molecules(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      integer :: ring(3, 15), i
      ! Acetonitrile
      call make_graph(mol, [6, 6, 7, 1, 1, 1], reshape([1, 2, 1, 2, 3, 3, 1, 4, 1, 1, 5, 1, 1, 6, 1], [3, 5]))
      call check_opls_labels(error, mol, [character(len=8) :: "opls_755", "opls_754", "opls_753", "opls_759", "opls_759", &
         & "opls_759"])
      if (allocated(error)) return
      ! DMF: O, carbonyl C, N, two methyl carbons, formyl H, six methyl H
      call make_graph(mol, [8, 6, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1], reshape([1, 2, 2, 2, 3, 1, 3, 4, 1, 3, 5, 1, 2, 6, 1, &
         & 4, 7, 1, 4, 8, 1, 4, 9, 1, 5, 10, 1, 5, 11, 1, 5, 12, 1], [3, 11]))
      call check_opls_labels(error, mol, [character(len=8) :: "opls_236", "opls_235", "opls_239", "opls_243", "opls_243", &
         & "opls_279", "*", "*", "*", "*", "*", "*"])
      if (allocated(error)) return
      ! Dimethyl ether
      call make_graph(mol, [8, 6, 6, 1, 1, 1, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 2, 4, 1, 2, 5, 1, 2, 6, 1, 3, 7, 1, &
         & 3, 8, 1, 3, 9, 1], [3, 8]))
      call check_opls_labels(error, mol, [character(len=8) :: "opls_180", "opls_181", "opls_181", "opls_185", "opls_185", &
         & "opls_185", "opls_185", "opls_185", "opls_185"])
      if (allocated(error)) return
      ! Methanethiol
      call make_graph(mol, [16, 6, 1, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 2, 4, 1, 2, 5, 1, 2, 6, 1], [3, 5]))
      call check_opls_labels(error, mol, [character(len=8) :: "opls_200", "opls_217", "opls_204", "*", "*", "*"])
      if (allocated(error)) return
      ! Methylamine
      call make_graph(mol, [7, 6, 1, 1, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 1, 4, 1, 2, 5, 1, 2, 6, 1, 2, 7, 1], [3, 6]))
      call check_opls_labels(error, mol, [character(len=8) :: "opls_900", "opls_903", "opls_909", "opls_909", "opls_911", &
         & "opls_911", "opls_911"])
      if (allocated(error)) return
      ! DMSO
      call make_graph(mol, [16, 8, 6, 6, 1, 1, 1, 1, 1, 1], reshape([1, 2, 2, 1, 3, 1, 1, 4, 1, 3, 5, 1, 3, 6, 1, 3, 7, 1, &
         & 4, 8, 1, 4, 9, 1, 4, 10, 1], [3, 9]))
      call check_opls_labels(error, mol, [character(len=8) :: "opls_496", "opls_497", "opls_498", "opls_498", "*", "*", "*", &
         & "*", "*", "*"])
      if (allocated(error)) return
      ! Sulfolane: sulfur, two oxygens, ring carbons 4-7 with two hydrogens each
      ring(:, 1:2) = reshape([1, 2, 2, 1, 3, 2], [3, 2])
      ring(:, 3:7) = reshape([1, 4, 1, 4, 5, 1, 5, 6, 1, 6, 7, 1, 7, 1, 1], [3, 5])
      do i = 1, 4
         ring(:, 7 + 2*i - 1) = [3 + i, 6 + 2*i, 1]
         ring(:, 7 + 2*i) = [3 + i, 7 + 2*i, 1]
      end do
      call make_graph(mol, [16, 8, 8, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1], ring)
      call check_opls_labels(error, mol, [character(len=8) :: "opls_493", "opls_494", "opls_494", "*", "*", "*", "*", "*", &
         & "*", "*", "*", "*", "*", "*", "*"])
   end subroutine check_small_molecules

   !> Check formal-charge handling: charged groups, ions and a cation
   !>
   !> - Expected types from the OPLS-AA parameter set: nitromethane C/N/O/H 762/760/761/763 (Jorgensen and
   !>   Tirado-Rives), acetate carboxylate C/O 271/272
   !> - Ions have no imported Lennard-Jones key, so they stay uncovered without a refusal
   !> - A formal charge outside the rules that name it leaves the atom to a later term
   subroutine check_charged_species(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      ! Nitromethane with the N(+)/O(-) Lewis structure
      call make_graph(mol, [6, 7, 8, 8, 1, 1, 1], reshape([1, 2, 1, 2, 3, 2, 2, 4, 1, 1, 5, 1, 1, 6, 1, 1, 7, 1], [3, 6]), &
         & [0, 1, 0, -1, 0, 0, 0])
      call check_opls_labels(error, mol, [character(len=8) :: "opls_762", "opls_760", "opls_761", "opls_761", "opls_763", &
         & "opls_763", "opls_763"])
      if (allocated(error)) return
      ! Acetate
      call make_graph(mol, [6, 6, 8, 8, 1, 1, 1], reshape([1, 2, 1, 2, 3, 2, 2, 4, 1, 1, 5, 1, 1, 6, 1, 1, 7, 1], [3, 6]), &
         & [0, 0, 0, -1, 0, 0, 0])
      call check_opls_labels(error, mol, [character(len=8) :: "*", "opls_271", "opls_272", "opls_272", "*", "*", "*"])
      if (allocated(error)) return
      ! Lithium cation and chloride anion
      call make_graph(mol, [3, 17], reshape([integer :: ], [3, 0]), [1, -1])
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(oplsaa%pair%covered), "ions without imported parameters must stay uncovered")
      if (allocated(error)) return
      ! Formamide iminolate zwitterion: charged O and N match neutral-only rules and stay untyped
      call make_graph(mol, [6, 8, 7, 1, 1, 1], reshape([1, 2, 1, 1, 3, 2, 1, 4, 1, 3, 5, 1, 3, 6, 1], [3, 5]), &
         & [0, -1, 1, 0, 0, 0])
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, oplsaa%pair%covered(1) .and. .not. oplsaa%pair%covered(2) .and. .not. oplsaa%pair%covered(3), &
         & "charged atoms matching only neutral rules must stay untyped")
      if (allocated(error)) return
      ! Hydronium: a valid valence that no rule types
      call make_graph(mol, [8, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 1, 4, 1], [3, 3]), [1, 0, 0, 0])
      call oplsaa%update(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(oplsaa%pair%covered), "hydronium must stay uncovered")
   end subroutine check_charged_species

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(moist_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> Coordinates without bonds are refused
   subroutine check_bad_no_bonds(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      real(wp), parameter :: xyz(3, 2) = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.4_wp, 0.0_wp, 0.0_wp], [3, 2])
      call new(mol, [1, 1], xyz)
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "needs explicit bonds")
   end subroutine check_bad_no_bonds

   !> A bond count that disagrees with the bond array is refused
   subroutine check_bad_bond_shape(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      mol%nbd = 2
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "bonds need shape (3, nbd)")
   end subroutine check_bad_bond_shape

   !> A structure flagged with missing hydrogens is refused
   subroutine check_bad_missing_hydrogen(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      mol%info%missing_hydrogen = .true.
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "explicit hydrogens and a closed-shell")
   end subroutine check_bad_missing_hydrogen

   !> An open-shell structure is refused
   subroutine check_bad_open_shell(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      mol%uhf = 1
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "explicit hydrogens and a closed-shell")
   end subroutine check_bad_open_shell

   !> Formal charge annotations of the wrong size are refused
   subroutine check_bad_charge_annotations(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      allocate (mol%sdf(3))
      mol%sdf%charge = 0
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "formal charge annotations do not match")
   end subroutine check_bad_charge_annotations

   !> Formal charges inconsistent with the total charge are refused
   subroutine check_bad_total_charge(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair, [1, -1])
      mol%charge = 1.0_wp
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "consistent with the total charge")
   end subroutine check_bad_total_charge

   !> A bond from an atom to itself is refused
   subroutine check_bad_self_bond(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [1, 1], reshape([1, 1, 1], [3, 1]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "invalid bond endpoint")
   end subroutine check_bad_self_bond

   !> A bond to an atom beyond the structure is refused
   subroutine check_bad_bond_endpoint(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [1, 1], reshape([1, 5, 1], [3, 1]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "invalid bond endpoint")
   end subroutine check_bad_bond_endpoint

   !> A bond order above the supported range is refused
   subroutine check_bad_bond_order_high(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [1, 1], reshape([1, 2, 5], [3, 1]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "bond orders must be")
   end subroutine check_bad_bond_order_high

   !> A zero bond order is refused
   subroutine check_bad_bond_order_zero(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [1, 1], reshape([1, 2, 0], [3, 1]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "bond orders must be")
   end subroutine check_bad_bond_order_zero

   !> A bond listed twice is refused
   subroutine check_bad_duplicate_bond(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [1, 1], reshape([1, 2, 1, 2, 1, 1], [3, 2]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "duplicate bond")
   end subroutine check_bad_duplicate_bond

   !> A species identity outside the species table is refused
   subroutine check_bad_species_table(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      mol%id(1) = 99
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "outside the species table")
   end subroutine check_bad_species_table

   !> An atom count that disagrees with the identities is refused
   subroutine check_bad_identity_count(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      mol%nat = 3
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "identity count")
   end subroutine check_bad_identity_count

   !> A structure without atoms is refused
   subroutine check_bad_empty_structure(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      integer, parameter :: pair(3, 1) = reshape([1, 2, 1], [3, 1])
      call make_graph(mol, [1, 1], pair)
      mol%nat = 0
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "at least one atom")
   end subroutine check_bad_empty_structure

   !> Carbon with a single aromatic bond is refused
   subroutine check_bad_valence_single_aromatic(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [6, 1, 1, 1, 6], reshape([1, 2, 1, 1, 3, 1, 1, 4, 1, 1, 5, 4], [3, 4]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_single_aromatic

   !> Carbon mixing aromatic and double bonds is refused
   subroutine check_bad_valence_aromatic_double(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [6, 6, 6, 6], reshape([1, 2, 4, 1, 3, 4, 1, 4, 2], [3, 3]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_aromatic_double

   !> Aromatic fluorine is refused
   subroutine check_bad_valence_aromatic_fluorine(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [9, 6, 6], reshape([1, 2, 4, 1, 3, 4], [3, 2]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_aromatic_fluorine

   !> Three-coordinate aromatic sulfur is refused
   subroutine check_bad_valence_aromatic_sulfur(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [16, 6, 6, 1], reshape([1, 2, 4, 1, 3, 4, 1, 4, 1], [3, 3]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_aromatic_sulfur

   !> Neutral ammonium is refused
   subroutine check_bad_valence_ammonium(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [7, 1, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 1, 4, 1, 1, 5, 1], [3, 4]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_ammonium

   !> Tetravalent phosphorus is refused
   subroutine check_bad_valence_phosphorus(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [15, 1, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 1, 4, 1, 1, 5, 1], [3, 4]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_phosphorus

   !> Trivalent sulfur is refused
   subroutine check_bad_valence_sulfur(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [16, 1, 1, 1], reshape([1, 2, 1, 1, 3, 1, 1, 4, 1], [3, 3]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_sulfur

   !> Divalent chlorine is refused
   subroutine check_bad_valence_chlorine(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [17, 1, 1], reshape([1, 2, 1, 1, 3, 1], [3, 2]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_chlorine

   !> A bridging hydrogen is refused
   subroutine check_bad_valence_bridging_hydrogen(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [1, 1, 1], reshape([1, 2, 1, 1, 3, 1], [3, 2]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_bridging_hydrogen

   !> Neutral lithium is refused
   subroutine check_bad_valence_neutral_lithium(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(lj_typed_type) :: oplsaa
      type(moist_error), allocatable :: err
      call make_graph(mol, [3], reshape([integer :: ], [3, 0]))
      oplsaa = lj_oplsaa
      call oplsaa%update(mol, err)
      call expect_error(error, err, "Unsupported OPLS-AA valence at atom 1")
   end subroutine check_bad_valence_neutral_lithium

end module test_moz_potential_oplsaa
