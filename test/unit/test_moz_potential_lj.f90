!> Tests of the MOZ Lennard-Jones potential terms: moz_potential_lj
!>
!> - GAFF reference types and parameters from moist_dev_moz on the same solvent geometries
!> - Reference sources: lj/atomtypeidentifier.f90, lj/parameters.f90, fftyping.f90
module test_moz_potential_lj
   use mctc_env, only: wp, moist_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau, kcaltoau, kjtoau
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use testdrive, only: unittest_type, new_unittest, error_type, check, test_failed
   use moist_data_solvents, only: solvation_system_type, new_solvation_system, get_solvent_id
   use moist_model_moz_potential_lj_base, only: lj_12_6_type, new_lj_12_6, lj_12_6_mix, lj_12_6_evaluate, &
      & lj_mixing_lorentz_berthelot, lj_mixing_geometric, check_lj_mixing
   use moist_model_moz_potential_lj_typed, only: lj_typed_type, lj_gaff
   use moist_model_moz_potential_lj_element, only: lj_element_type, lj_uff, lj_dreiding, lj_tm
   use moist_model_moz_potential_lj_solvent, only: lj_solvent_type, lj_spce
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_set, only: potential_set_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use test_helpers, only: check_moist_error
   implicit none(type, external)
   private

   public :: collect_moz_potential_lj

   !> Relative tolerance on tabulated parameters
   real(wp), parameter :: thr_rel = 1.0e-14_wp

contains

   !> Register tests
   subroutine collect_moz_potential_lj(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         & new_unittest("gaff-water", check_gaff_water), &
         & new_unittest("spce-water", check_spce_water), &
         & new_unittest("spce-invalid-solvent", check_spce_invalid_solvent), &
         & new_unittest("set-water-model", check_set_water_model), &
         & new_unittest("gaff-methanol", check_gaff_methanol), &
         & new_unittest("gaff-benzene", check_gaff_benzene), &
         & new_unittest("gaff-acetonitrile", check_gaff_acetonitrile), &
         & new_unittest("gaff-dmso", check_gaff_dmso), &
         & new_unittest("gaff-pyridine", check_gaff_pyridine), &
         & new_unittest("gaff-chloroform", check_gaff_chloroform), &
         & new_unittest("gaff-element-uncovered", check_gaff_element_uncovered), &
         & new_unittest("gaff-water-fragment-uncovered", check_gaff_water_fragment), &
         & new_unittest("set-fallback", check_set_fallback), &
         & new_unittest("set-gaff-alone-uncovered", check_set_uncovered), &
         & new_unittest("custom-wrong-size", check_custom_wrong_size), &
         & new_unittest("set-mixing", check_set_mixing), &
         & new_unittest("epsilon-zero", check_epsilon_zero), &
         & new_unittest("mixing-rules", check_mixing_rules), &
         & new_unittest("inactive-site-mixing", check_inactive_site_mixing), &
         & new_unittest("element-tables", check_element_tables), &
         & new_unittest("element-fallback", check_element_fallback)]
   end subroutine collect_moz_potential_lj

   !> Read solvent geometry in bohr
   !>
   !> @param[out] error  test error
   !> @param[in]  name   solvent alias
   !> @param[out] mol    solvent geometry, bohr
   subroutine solvent_mol(error, name, mol)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Solvent alias
      character(len=*), intent(in) :: name
      !> Solvent geometry, bohr
      type(structure_type), intent(out) :: mol
      type(moist_error), allocatable :: err
      type(solvation_system_type) :: system
      integer :: id
      call get_solvent_id(name, id, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call new_solvation_system(system, id, error=err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      mol = system%solv_mol
   end subroutine solvent_mol

   !> Look up reference sigma (bohr) and epsilon (Hartree) by dev-fork type label
   !>
   !> @param[in]  label    type label "type/source"
   !> @param[out] sigma    sigma, bohr
   !> @param[out] epsilon  epsilon, Hartree
   subroutine reference_gaff(label, sigma, epsilon)
      !> Type label
      character(len=*), intent(in) :: label
      !> Sigma, bohr
      real(wp), intent(out) :: sigma
      !> Epsilon, Hartree
      real(wp), intent(out) :: epsilon
      select case (label)
      case ("c3/GAFF")
         sigma = 6.4207404650629512e+00_wp; epsilon = 1.7179023497799125e-04_wp
      case ("oh/GAFF")
         sigma = 6.1281386787008074e+00_wp; epsilon = 1.4820493370086443e-04_wp
      case ("hc/GAFF")
         sigma = 4.9136224032022460e+00_wp; epsilon = 3.3146909902989032e-05_wp
      case ("ho/GAFF")
         sigma = 1.0165302566482273e+00_wp; epsilon = 7.4899267569254082e-06_wp
      case ("ca/GAFF", "cp/GAFF", "cq/GAFF")
         sigma = 6.2648433107641335e+00_wp; epsilon = 1.5744782203919793e-04_wp
      case ("ha/GAFF")
         sigma = 4.9614353533327691e+00_wp; epsilon = 2.5656983146063627e-05_wp
      case ("n/GAFF")
         sigma = 6.0109632797893857e+00_wp; epsilon = 2.6071319519850990e-04_wp
      case ("c1/GAFF")
         sigma = 6.5742806429468841e+00_wp; epsilon = 2.5433878944793509e-04_wp
      case ("s4/GAFF")
         sigma = 6.6752939178705226e+00_wp; epsilon = 4.5003304599058191e-04_wp
      case ("o/GAFF")
         sigma = 5.7601136470623473e+00_wp; epsilon = 2.3314389032727387e-04_wp
      case ("na/GAFF")
         sigma = 6.0581028080870842e+00_wp; epsilon = 3.2541341356684427e-04_wp
      case ("cl/GAFF")
         sigma = 6.5497007460487975e+00_wp; epsilon = 4.2039205925040895e-04_wp
      case default
         sigma = -1.0_wp; epsilon = -1.0_wp
      end select
   end subroutine reference_gaff

   !> Build a GAFF term and compare types and parameters with the dev reference
   !>
   !> @param[out] error  test error
   !> @param[in]  mol    structure
   !> @param[in]  types  expected type label per atom
   subroutine check_gaff(error, mol, types)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Expected type label per atom
      character(len=*), intent(in) :: types(:)
      type(moist_error), allocatable :: err
      type(lj_typed_type) :: gaff
      real(wp) :: sigma, epsilon
      integer :: iat
      call check(error, size(types), mol%nat, "reference covers a different number of atoms")
      if (allocated(error)) return
      gaff = lj_gaff
      if (.not. allocated(err)) call gaff%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(gaff%pair), "GAFF supplied no Lennard-Jones data")
      if (allocated(error)) return
      call check(error, all(gaff%pair%covered), "GAFF leaves a typed atom uncovered")
      if (allocated(error)) return
      do iat = 1, mol%nat
         call check(error, gaff%atomtype(iat), types(iat), "GAFF type differs from the dev reference")
         if (allocated(error)) return
         call reference_gaff(types(iat), sigma, epsilon)
         call check(error, gaff%pair%sigma(iat), sigma, "GAFF sigma differs", thr=thr_rel, rel=.true.)
         if (allocated(error)) return
         call check(error, gaff%pair%epsilon(iat), epsilon, "GAFF epsilon differs", thr=thr_rel, rel=.true.)
         if (allocated(error)) return
      end do
   end subroutine check_gaff

   !> Check uncovered GAFF water placeholders, including pure water
   !>
   !> @param[out] error  test error
   subroutine check_gaff_water(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(lj_typed_type) :: gaff
      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      gaff = lj_gaff
      if (.not. allocated(err)) call gaff%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error,.not. any(gaff%pair%covered), "GAFF must leave water sites uncovered")
      if (allocated(error)) return
      call check(error, all(gaff%pair%sigma == 0.0_wp) .and. all(gaff%pair%epsilon == 0.0_wp), &
         & "GAFF must retain zero-valued water placeholders")
      if (allocated(error)) return
      call check(error, all(gaff%atomtype == [character(len=8) :: "ow", "hw", "hw"]), &
         & "GAFF water labels must not name a solvent model")
   end subroutine check_gaff_water

   !> Check legacy water-model values and atom ordering
   !>
   !> @param[out] error  test error
   subroutine check_spce_water(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol, permuted
      type(moist_error), allocatable :: err
      type(lj_solvent_type) :: spce
      integer, parameter :: order(3) = [3, 1, 2]
      real(wp), parameter :: sigma_ref(3) = [5.9823437872334209_wp, 2.2028915379925396_wp, &
         & 2.2028915379925396_wp]
      real(wp), parameter :: epsilon_ref(3) = [2.4748632292216421e-4_wp, 2.4748632292216419e-5_wp, &
         & 2.4748632292216419e-5_wp]
      character(len=8), parameter :: labels(3) = ["ow/SPCE ", "hw/SPCE ", "hw/SPCE "]
      integer :: water_id, iat

      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      call get_solvent_id("water", water_id, err)
      if (.not. allocated(err)) spce = lj_spce
      if (.not. allocated(err)) call spce%build(mol, err, water_id)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(spce%pair%covered), "Explicit water model must cover every site")
      if (allocated(error)) return
      do iat = 1, mol%nat
         call check(error, spce%atomtype(iat), labels(iat), "Water site label")
         if (allocated(error)) return
         call check(error, spce%pair%sigma(iat), sigma_ref(iat), "Legacy water sigma in bohr", &
            & thr=thr_rel, rel=.true.)
         if (allocated(error)) return
         call check(error, spce%pair%epsilon(iat), epsilon_ref(iat), "Legacy water epsilon in Hartree", &
            & thr=thr_rel, rel=.true.)
         if (allocated(error)) return
      end do

      call new(permuted, mol%num(mol%id(order)), mol%xyz(:, order))
      call spce%build(permuted, err, solvent_id=0)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do iat = 1, permuted%nat
         call check(error, spce%atomtype(iat), labels(order(iat)), "Permuted water site label")
         if (allocated(error)) return
         call check(error, spce%pair%sigma(iat), sigma_ref(order(iat)), "Permuted water sigma", &
            & thr=thr_rel, rel=.true.)
         if (allocated(error)) return
         call check(error, spce%pair%epsilon(iat), epsilon_ref(order(iat)), "Permuted water epsilon", &
            & thr=thr_rel, rel=.true.)
         if (allocated(error)) return
      end do
      call spce%build(mol, err)
      call check_moist_error(error, err)
   end subroutine check_spce_water

   !> Check invalid solvent selection and clearing of built water data
   !>
   !> @param[out] error  test error
   subroutine check_spce_invalid_solvent(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol, other
      type(moist_error), allocatable :: err
      type(lj_solvent_type) :: spce
      integer :: methanol_id

      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      call spce%build(mol, err)
      call check(error, allocated(err), "Unconstructed solvent term must be refused")
      if (allocated(error)) return
      deallocate (err)
      spce = lj_solvent_type(set=99)
      call spce%build(mol, err)
      call check(error, allocated(err), "Unknown solvent model must be refused")
      if (allocated(error)) return
      call check(error, index(err%message, "Unknown Lennard-Jones solvent model 99") > 0, err%message)
      if (allocated(error)) return
      deallocate (err)
      spce = lj_spce
      call spce%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(spce%pair%covered), "Water must be covered")
      if (allocated(error)) return
      ! Non-neutral or non-water structure: atoms left for later terms
      call get_solvent_id("methanol", methanol_id, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check_uncovered(error, spce, mol, "a different solvent ID", methanol_id)
      if (allocated(error)) return
      call solvent_mol(error, "methanol", other)
      if (allocated(error)) return
      call check_uncovered(error, spce, other, "another molecule")
      if (allocated(error)) return
      call new(other, [6, 1, 1], mol%xyz)
      call check_uncovered(error, spce, other, "wrong elements with the water atom count")
      if (allocated(error)) return
      mol%charge = 1.0_wp
      call check_uncovered(error, spce, mol, "charged H2O")
      if (allocated(error)) return
      mol%charge = 0.0_wp
      mol%uhf = 2
      call check_uncovered(error, spce, mol, "open-shell H2O")
   end subroutine check_spce_invalid_solvent

   !> Check zero SPC/E coverage for an unsupported structure
   !>
   !> @param[out]    error       test error
   !> @param[in,out] spce        constructed SPC/E term
   !> @param[in]     mol         structure
   !> @param[in]     what        case description
   !> @param[in]     solvent_id  optional solvent id of the table
   subroutine check_uncovered(error, spce, mol, what, solvent_id)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Constructed SPC/E term
      type(lj_solvent_type), intent(inout) :: spce
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Case description
      character(len=*), intent(in) :: what
      !> Optional solvent id of the table
      integer, intent(in), optional :: solvent_id
      type(moist_error), allocatable :: err
      call spce%build(mol, err, solvent_id)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, allocated(spce%pair), "Water model must supply LJ data for "//what)
      if (allocated(error)) return
      call check(error, .not. any(spce%pair%covered), "Water model must leave "//what//" uncovered")
      if (allocated(error)) return
      call check(error, all(len_trim(spce%atomtype) == 0), "Uncovered atoms must carry blank labels for "//what)
   end subroutine check_uncovered

   !> Check explicit water-model selection after GAFF in a set
   !>
   !> @param[out] error  test error
   subroutine check_set_water_model(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(lj_typed_type) :: gaff
      type(lj_solvent_type) :: spce
      type(potential_set_type) :: set
      type(potential_sites_type) :: sites
      integer :: water_id

      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      call get_solvent_id("water", water_id, err)
      if (.not. allocated(err)) gaff = lj_gaff
      if (.not. allocated(err)) call set%add(gaff, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call set%build(mol, err, water_id)
      call check(error, allocated(err), "GAFF alone must not provide a water solvent model")
      if (allocated(error)) return
      deallocate (err)
      spce = lj_spce
      call set%add(spce, err)
      if (.not. allocated(err)) call set%build(mol, err, water_id)
      if (.not. allocated(err)) call set%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sites%pair%covered), "Explicit solvent term must fill the GAFF water sites")
      if (allocated(error)) return
      call check(error, sites%pair%sigma(1), 5.9823437872334209_wp, "Merged water oxygen sigma", &
         & thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call check(error, sites%pair%epsilon(2), 2.4748632292216419e-5_wp, "Merged water hydrogen epsilon", &
         & thr=thr_rel, rel=.true.)
   end subroutine check_set_water_model

   !> Check methanol types: sp3 carbon, hydroxyl O and H, aliphatic H
   subroutine check_gaff_methanol(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      call solvent_mol(error, "methanol", mol)
      if (allocated(error)) return
      call check_gaff(error, mol, [character(len=8) :: "c3/GAFF", "oh/GAFF", "hc/GAFF", "hc/GAFF", &
         & "hc/GAFF", "ho/GAFF"])
   end subroutine check_gaff_methanol

   !> Check benzene types: aromatic carbons with order-dependent cp/cq alternation
   subroutine check_gaff_benzene(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      call solvent_mol(error, "benzene", mol)
      if (allocated(error)) return
      call check_gaff(error, mol, [character(len=8) :: "cp/GAFF", "cq/GAFF", "cp/GAFF", "cq/GAFF", &
         & "ca/GAFF", "cq/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF"])
   end subroutine check_gaff_benzene

   !> Check acetonitrile types: terminal N, sp carbon
   subroutine check_gaff_acetonitrile(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      call solvent_mol(error, "acetonitrile", mol)
      if (allocated(error)) return
      call check_gaff(error, mol, [character(len=8) :: "n/GAFF", "c1/GAFF", "c3/GAFF", "hc/GAFF", &
         & "hc/GAFF", "hc/GAFF"])
   end subroutine check_gaff_acetonitrile

   !> Check DMSO types: three-coordinate sulfur, terminal oxygen
   subroutine check_gaff_dmso(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      call solvent_mol(error, "dimethylsulfoxide", mol)
      if (allocated(error)) return
      call check_gaff(error, mol, [character(len=8) :: "c3/GAFF", "hc/GAFF", "hc/GAFF", "hc/GAFF", &
         & "s4/GAFF", "o/GAFF", "c3/GAFF", "hc/GAFF", "hc/GAFF", "hc/GAFF"])
   end subroutine check_gaff_dmso

   !> Check pyridine types: aromatic nitrogen
   subroutine check_gaff_pyridine(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      call solvent_mol(error, "pyridine", mol)
      if (allocated(error)) return
      call check_gaff(error, mol, [character(len=8) :: "na/GAFF", "cp/GAFF", "cq/GAFF", "cp/GAFF", &
         & "cq/GAFF", "ca/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF", "ha/GAFF"])
   end subroutine check_gaff_pyridine

   !> Check chloroform types: chlorine
   subroutine check_gaff_chloroform(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      call solvent_mol(error, "chloroform", mol)
      if (allocated(error)) return
      call check_gaff(error, mol, [character(len=8) :: "c3/GAFF", "hc/GAFF", "cl/GAFF", "cl/GAFF", &
         & "cl/GAFF"])
   end subroutine check_gaff_chloroform

   !> Check UFF fallback for elements without GAFF rules
   !>
   !> - UFF rows only for atoms left uncovered by GAFF
   subroutine check_gaff_element_uncovered(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(lj_typed_type) :: gaff
      type(lj_element_type) :: uff, tm
      type(potential_set_type) :: set
      type(potential_sites_type) :: sites
      ! Methane carbon and hydrogen, argon and iron
      call new(mol, [6, 1, 18, 26], reshape([0.0_wp, 0.0_wp, 0.0_wp, 2.05_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 12.0_wp, 0.0_wp, 12.0_wp, 0.0_wp], [3, 4]))
      gaff = lj_gaff
      if (.not. allocated(err)) call gaff%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(gaff%pair%covered .eqv. [.true., .true., .false., .false.]), &
         & "GAFF must cover only the atoms with a GAFF2 row")
      if (allocated(error)) return
      call check(error, gaff%atomtype(3), "du", "An uncovered atom keeps its bare GAFF type")
      if (allocated(error)) return
      uff = lj_uff
      if (.not. allocated(err)) call uff%build(mol, err)
      if (.not. allocated(err)) call set%add(gaff, err)
      if (.not. allocated(err)) call set%add(uff, err)
      if (.not. allocated(err)) call set%build(mol, err)
      if (.not. allocated(err)) call set%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sites%pair%covered), "GAFF then UFF must cover every atom")
      if (allocated(error)) return
      call check(error, sites%pair%sigma(1), gaff%pair%sigma(1), "GAFF keeps the carbon it typed")
      if (allocated(error)) return
      call check(error, sites%pair%epsilon(2), gaff%pair%epsilon(2), "GAFF keeps the hydrogen it typed")
      if (allocated(error)) return
      call check(error, sites%pair%sigma(3), uff%pair%sigma(3), "UFF fills argon")
      if (allocated(error)) return
      call check(error, sites%pair%epsilon(4), uff%pair%epsilon(4), "UFF fills iron")
      if (allocated(error)) return
      call check(error, sites%pair%label(1), gaff%atomtype(1), "Merged label of the GAFF carbon")
      if (allocated(error)) return
      call check(error, sites%pair%label(3), "Ar/UFF", "Merged label of the UFF argon")
      if (allocated(error)) return
      call check_site_table(error, set, mol)
      if (allocated(error)) return
      ! Separate transition-metal table: iron coverage only
      tm = lj_tm
      if (.not. allocated(err)) call tm%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(tm%pair%covered .eqv. [.false., .false., .false., .true.]), &
         & "The TM table must cover only iron")
      if (allocated(error)) return
      call check(error, tm%atomtype(4), "Fe/TM", "TM label")
      if (allocated(error)) return
      call check(error, tm%pair%sigma(4), 0.5938_wp*10.0_wp*aatoau, "TM iron sigma", thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call check(error, tm%pair%epsilon(4), 2.920432_wp*kjtoau, "TM iron epsilon", thr=thr_rel, rel=.true.)
   end subroutine check_gaff_element_uncovered

   !> Print the built site table to a scratch unit and check rows
   !>
   !> @param[out] error  test error
   !> @param[in]  set    built set of methane carbon and hydrogen, argon and iron
   !> @param[in]  mol    its structure
   subroutine check_site_table(error, set, mol)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Built set
      type(potential_set_type), intent(in) :: set
      !> Structure
      type(structure_type), intent(in) :: mol
      type(moist_error), allocatable :: err
      character(len=96) :: line
      integer :: unit, stat
      open (newunit=unit, status="scratch", action="readwrite", form="formatted", iostat=stat)
      call check(error, stat, 0, "Cannot open a scratch unit")
      if (allocated(error)) return
      call set%print_table(mol, unit, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      rewind (unit)
      read (unit, "(a)", iostat=stat) line
      call check(error, stat, 0, "Table header missing")
      if (allocated(error)) return
      call check(error, index(line, "N") > 0 .and. index(line, "Sym") > 0 .and. index(line, "qat") > 0 &
         & .and. index(line, "sigma (bohr)") > 0 .and. index(line, "eps (kcal/mol)") > 0 &
         & .and. index(line, "Type/Src") > 0, more="Table header differs: "//trim(line))
      if (allocated(error)) return
      read (unit, "(a)", iostat=stat) line
      call check(error, stat, 0, "First row missing")
      if (allocated(error)) return
      ! No electrostatics: "-" in charge column
      call check(error, index(line, " 1      C") > 0 .and. index(line, " -") > 0 .and. index(line, "/GAFF") > 0, &
         & more="First row differs: "//trim(line))
      if (allocated(error)) return
      read (unit, "(a)", iostat=stat) line
      call check(error, stat, 0, "Second row missing")
      if (allocated(error)) return
      read (unit, "(a)", iostat=stat) line
      call check(error, stat, 0, "Third row missing")
      if (allocated(error)) return
      call check(error, index(line, " 3     Ar") > 0 .and. index(line, "Ar/UFF") > 0, more="Third row differs: "//trim(line))
      if (allocated(error)) return
      read (unit, "(a)", iostat=stat) line
      call check(error, stat, 0, "Fourth row missing")
      if (allocated(error)) return
      call check(error, index(line, "Fe/UFF") > 0, more="Fourth row differs: "//trim(line))
      close (unit)
   end subroutine check_site_table

   !> Check uncovered GAFF water placeholders and atoms without GAFF2 rows
   subroutine check_gaff_water_fragment(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: water, mol
      type(moist_error), allocatable :: err
      type(lj_typed_type) :: gaff
      call solvent_mol(error, "water", water)
      if (allocated(error)) return
      call new(mol, [water%num(water%id), 18], &
         & reshape([water%xyz, 0.0_wp, 0.0_wp, 15.0_wp], [3, water%nat + 1]))
      gaff = lj_gaff
      if (.not. allocated(err)) call gaff%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, .not. any(gaff%pair%covered), "GAFF covers neither the water sites nor argon")
      if (allocated(error)) return
      call check(error, gaff%atomtype(4), "du", "Uncovered argon keeps its bare GAFF type")
      if (allocated(error)) return
      call check(error, gaff%atomtype(1), "ow", "Uncovered oxygen keeps its bare GAFF type")
      if (allocated(error)) return
      call check(error, gaff%atomtype(2), "hw", "Uncovered hydrogen keeps its bare GAFF type")
   end subroutine check_gaff_water_fragment

   !> Build methanol plus rutherfordium without a GAFF2 row
   !>
   !> @param[out] error  test error
   !> @param[out] mol    structure, rutherfordium last
   subroutine methanol_with_rf(error, mol)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Structure, rutherfordium last
      type(structure_type), intent(out) :: mol
      type(structure_type) :: methanol
      call solvent_mol(error, "methanol", methanol)
      if (allocated(error)) return
      call new(mol, [methanol%num(methanol%id), 104], &
         & reshape([methanol%xyz, 0.0_wp, 0.0_wp, 15.0_wp], [3, methanol%nat + 1]))
   end subroutine methanol_with_rf

   !> Check custom fallback restricted to atoms left uncovered by GAFF
   subroutine check_set_fallback(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(potential_set_type) :: set
      type(potential_sites_type) :: sites
      type(lj_typed_type) :: gaff
      type(custom_lj_type) :: custom
      logical, allocatable :: covered(:)
      real(wp), allocatable :: sigma(:), epsilon(:)
      real(wp) :: sigma_ref, epsilon_ref
      call methanol_with_rf(error, mol)
      if (allocated(error)) return
      allocate (sigma(mol%nat), epsilon(mol%nat), covered(mol%nat))
      sigma = 1.0_wp
      epsilon = 1.0_wp
      sigma(mol%nat) = 7.0_wp
      epsilon(mol%nat) = 3.0e-4_wp
      covered = .false.
      covered(mol%nat) = .true.
      gaff = lj_gaff
      if (.not. allocated(err)) call new_custom_lj(custom, sigma, epsilon, err, covered)
      if (.not. allocated(err)) call set%add(gaff, err)
      if (.not. allocated(err)) call set%add(custom, err)
      if (.not. allocated(err)) call set%build(mol, err)
      if (.not. allocated(err)) call set%sites(sites, err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call check(error, all(sites%pair%covered), "Merged Lennard-Jones data is incomplete")
      if (allocated(error)) return
      call reference_gaff("c3/GAFF", sigma_ref, epsilon_ref)
      call check(error, sites%pair%sigma(1), sigma_ref, "GAFF atom lost its parameters", thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call check(error, sites%pair%sigma(mol%nat), 7.0_wp, "Custom term did not fill the rutherfordium atom")
      if (allocated(error)) return
      call check(error, sites%pair%epsilon(mol%nat), 3.0e-4_wp, "Custom term did not fill the rutherfordium atom")
   end subroutine check_set_fallback

   !> Check named set failure for rutherfordium left uncovered by GAFF
   subroutine check_set_uncovered(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(potential_set_type) :: set
      type(lj_typed_type) :: gaff
      call methanol_with_rf(error, mol)
      if (allocated(error)) return
      gaff = lj_gaff
      if (.not. allocated(err)) call set%add(gaff, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call set%build(mol, err)
      call check(error, allocated(err), "An uncovered atom must fail the set build")
      if (allocated(error)) return
      call check(error, index(err%message, "No pair potential parameters for atoms: Rf7") > 0, &
         & "Error should name the atom: "//err%message)
   end subroutine check_set_uncovered

   !> Check inconsistent custom-parameter sizes
   subroutine check_custom_wrong_size(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(custom_lj_type) :: custom
      call new_custom_lj(custom, [6.0_wp, 2.0_wp], [1.0e-4_wp, 2.0e-5_wp, 2.0e-5_wp], err)
      call check(error, allocated(err), "Sigma and epsilon of different lengths must be refused")
      if (allocated(error)) return
      deallocate (err)
      call new_custom_lj(custom, [6.0_wp, 2.0_wp], [1.0e-4_wp, 2.0e-5_wp], err, covered=[.true.])
      call check(error, allocated(err), "A coverage mask of another length must be refused")
      if (allocated(error)) return
      deallocate (err)
      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      call new_custom_lj(custom, [6.0_wp, 2.0_wp], [1.0e-4_wp, 2.0e-5_wp], err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if
      call custom%build(mol, err)
      call check(error, allocated(err), "Parameters for another number of atoms must be refused at build")
   end subroutine check_custom_wrong_size

   !> Check Lorentz-Berthelot mixing across merged terms
   subroutine check_set_mixing(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_error), allocatable :: err
      type(potential_set_type) :: set
      type(potential_sites_type) :: sites
      type(custom_lj_type) :: oxygen, hydrogens
      real(wp), parameter :: r(2) = [5.5_wp, 8.0_wp]
      real(wp) :: u(2), du(2), s(3), e(3), sigma, epsilon, sr6, sigma_ij, eps_ij
      integer :: pairs(2, 2), ip, i, j, ir
      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      ! Separate custom terms for O and both hydrogens
      s = [3.166_wp*aatoau, 2.0_wp, 2.5_wp]
      e = [0.1554_wp*kcaltoau, 2.0e-5_wp, 3.0e-5_wp]
      call new_custom_lj(oxygen, s, e, err, covered=[.true., .false., .false.])
      if (.not. allocated(err)) call new_custom_lj(hydrogens, s, e, err, covered=[.false., .true., .true.])
      if (.not. allocated(err)) call set%add(oxygen, err)
      if (.not. allocated(err)) call set%add(hydrogens, err)
      if (.not. allocated(err)) call set%build(mol, err)
      if (.not. allocated(err)) call set%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      pairs = reshape([1, 2, 2, 3], [2, 2])
      do ip = 1, 2
         i = pairs(1, ip)
         j = pairs(2, ip)
         call lj_12_6_mix(lj_mixing_lorentz_berthelot, sites%pair%sigma(i), sites%pair%epsilon(i), &
            & sites%pair%sigma(j), sites%pair%epsilon(j), sigma_ij, eps_ij)
         call lj_12_6_evaluate(sigma_ij, eps_ij, r, u, du)
         sigma = 0.5_wp*(s(i) + s(j))
         epsilon = sqrt(e(i)*e(j))
         do ir = 1, size(r)
            sr6 = (sigma/r(ir))**6
            call check(error, u(ir), 4.0_wp*epsilon*(sr6*sr6 - sr6), "Mixed u(r) differs", &
               & thr=thr_rel, rel=.true.)
            if (allocated(error)) return
            call check(error, du(ir), -24.0_wp*epsilon*(2.0_wp*sr6*sr6 - sr6)/r(ir), "Mixed du/dr differs", &
               & thr=thr_rel, rel=.true.)
            if (allocated(error)) return
         end do
      end do
   end subroutine check_set_mixing

   !> Check active and inactive LJ parameter validation
   !>
   !> - epsilon = 0: inactive pair, any sigma accepted
   !> - epsilon < 0 or sigma <= 0 with epsilon > 0: refusal
   subroutine check_epsilon_zero(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(lj_12_6_type) :: lj
      real(wp), parameter :: r(3) = [0.5_wp, 2.0_wp, 6.0_wp]
      real(wp) :: u(3), du(3), sigma_ij, eps_ij
      call new_lj_12_6(lj, [6.0_wp, 0.0_wp], [2.0e-4_wp, 0.0_wp], err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call lj_12_6_mix(lj_mixing_lorentz_berthelot, lj%sigma(1), lj%epsilon(1), lj%sigma(2), lj%epsilon(2), &
         & sigma_ij, eps_ij)
      call lj_12_6_evaluate(sigma_ij, eps_ij, r, u, du)
      call check(error, all(u == 0.0_wp) .and. all(du == 0.0_wp), "epsilon = 0 must switch the pair off")
      if (allocated(error)) return
      call lj_12_6_mix(lj_mixing_geometric, lj%sigma(1), lj%epsilon(1), lj%sigma(2), lj%epsilon(2), &
         & sigma_ij, eps_ij)
      call lj_12_6_evaluate(sigma_ij, eps_ij, r, u, du)
      call check(error, all(u == 0.0_wp) .and. all(du == 0.0_wp), "epsilon = 0 must switch the pair off (geometric)")
      if (allocated(error)) return
      call new_lj_12_6(lj, [6.0_wp, 2.0_wp], [2.0e-4_wp, -1.0e-5_wp], err)
      call check(error, allocated(err), "A negative epsilon must be refused")
      if (allocated(error)) return
      deallocate (err)
      call new_lj_12_6(lj, [6.0_wp, 0.0_wp], [2.0e-4_wp, 1.0e-5_wp], err)
      call check(error, allocated(err), "sigma = 0 with epsilon > 0 must be refused")
      if (allocated(error)) return
      deallocate (err)
      ! Uncovered atom: no validation
      call new_lj_12_6(lj, [6.0_wp, 0.0_wp], [2.0e-4_wp, 1.0e-5_wp], err, covered=[.true., .false.])
      call check_moist_error(error, err)
   end subroutine check_epsilon_zero

   !> Check two-site cross parameters and u(r)
   !>
   !> - Geometric and Lorentz-Berthelot rules
   !> - Unknown rule: refusal and NaN cross parameters
   subroutine check_mixing_rules(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      real(wp), parameter :: s1 = 6.0_wp, s2 = 2.0_wp, e1 = 2.4e-4_wp, e2 = 1.0e-5_wp
      real(wp), parameter :: r(3) = [3.0_wp, 4.5_wp, 9.0_wp]
      real(wp) :: u(3), sigma_ij, eps_ij, sr6
      integer :: ir
      call lj_12_6_mix(lj_mixing_lorentz_berthelot, s1, e1, s2, e2, sigma_ij, eps_ij)
      call check(error, sigma_ij, 4.0_wp, "Lorentz-Berthelot sigma")
      if (allocated(error)) return
      call check(error, eps_ij, sqrt(e1*e2), "Lorentz-Berthelot epsilon", thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call lj_12_6_evaluate(sigma_ij, eps_ij, r, u)
      do ir = 1, size(r)
         sr6 = (4.0_wp/r(ir))**6
         call check(error, u(ir), 4.0_wp*sqrt(e1*e2)*(sr6*sr6 - sr6), "Lorentz-Berthelot u(r)", &
            & thr=thr_rel, rel=.true.)
         if (allocated(error)) return
      end do
      call lj_12_6_mix(lj_mixing_geometric, s1, e1, s2, e2, sigma_ij, eps_ij)
      call check(error, sigma_ij, sqrt(12.0_wp), "Geometric sigma", thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call lj_12_6_evaluate(sigma_ij, eps_ij, r, u)
      do ir = 1, size(r)
         sr6 = 12.0_wp**3/r(ir)**6
         call check(error, u(ir), 4.0_wp*sqrt(e1*e2)*(sr6*sr6 - sr6), "Geometric u(r)", &
            & thr=thr_rel, rel=.true.)
         if (allocated(error)) return
      end do
      call lj_12_6_mix(0, s1, e1, s2, e2, sigma_ij, eps_ij)
      call check(error, ieee_is_nan(sigma_ij) .and. ieee_is_nan(eps_ij), "An unknown rule must mix to NaN")
      if (allocated(error)) return
      call check_lj_mixing(3, err)
      call check(error, allocated(err), "An unknown mixing rule must be refused")
   end subroutine check_mixing_rules

   !> Check element-table rows and labels
   !>
   !> - UFF and DREIDING lookup by atomic number
   !> - Unknown table: refusal
   subroutine check_element_tables(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: mol
      type(lj_element_type) :: uff, dreiding
      call solvent_mol(error, "water", mol)
      if (allocated(error)) return
      uff = lj_element_type(set=7)
      call uff%build(mol, err)
      call check(error, allocated(err), "An unknown element table must be refused")
      if (allocated(error)) return
      deallocate (err)
      uff = lj_uff
      dreiding = lj_dreiding
      call uff%build(mol, err)
      if (.not. allocated(err)) call dreiding%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(uff%pair%covered) .and. all(dreiding%pair%covered), "Water must be covered")
      if (allocated(error)) return
      call check(error, uff%atomtype(1), "O/UFF", "UFF label")
      if (allocated(error)) return
      call check(error, dreiding%atomtype(2), "H/DREIDING", "DREIDING label")
      if (allocated(error)) return
      call check(error, uff%pair%sigma(1), 3.5_wp/2.0_wp**(1.0_wp/6.0_wp)*aatoau, "UFF oxygen sigma", &
         & thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call check(error, uff%pair%epsilon(1), 0.06_wp*kcaltoau, "UFF oxygen epsilon", thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call check(error, dreiding%pair%sigma(1), 3.4046_wp/2.0_wp**(1.0_wp/6.0_wp)*aatoau, "DREIDING oxygen sigma", &
         & thr=thr_rel, rel=.true.)
      if (allocated(error)) return
      call check(error, dreiding%pair%epsilon(2), 0.0152_wp*kcaltoau, "DREIDING hydrogen epsilon", &
         & thr=thr_rel, rel=.true.)
   end subroutine check_element_tables

   !> Check DREIDING argon fallback
   !>
   !> - No argon row: uncovered atom
   !> - Later custom term: argon coverage
   !> - No fallback: named set error
   subroutine check_element_fallback(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(structure_type) :: water, mol
      type(lj_element_type) :: dreiding
      type(custom_lj_type) :: custom
      type(potential_set_type) :: set, alone
      type(potential_sites_type) :: sites
      real(wp), allocatable :: sigma(:), epsilon(:)
      call solvent_mol(error, "water", water)
      if (allocated(error)) return
      call new(mol, [water%num(water%id), 18], &
         & reshape([water%xyz, 0.0_wp, 0.0_wp, 15.0_wp], [3, water%nat + 1]))
      dreiding = lj_dreiding
      if (.not. allocated(err)) call dreiding%build(mol, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(dreiding%pair%covered .eqv. [.true., .true., .true., .false.]), &
         & "Only argon should be uncovered")
      if (allocated(error)) return
      call check(error, dreiding%atomtype(4), "", "An uncovered atom keeps a blank label")
      if (allocated(error)) return

      allocate (sigma(mol%nat), source=0.0_wp)
      allocate (epsilon(mol%nat), source=0.0_wp)
      sigma(4) = 6.4_wp
      epsilon(4) = 3.7e-4_wp
      call new_custom_lj(custom, sigma, epsilon, err, covered=[.false., .false., .false., .true.])
      if (.not. allocated(err)) call set%add(dreiding, err)
      if (.not. allocated(err)) call set%add(custom, err)
      if (.not. allocated(err)) call set%build(mol, err)
      if (.not. allocated(err)) call set%sites(sites, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call check(error, all(sites%pair%covered), "The custom term must fill argon")
      if (allocated(error)) return
      call check(error, sites%pair%sigma(4), 6.4_wp, "Argon sigma from the custom term")
      if (allocated(error)) return
      call check(error, sites%pair%sigma(1), dreiding%pair%sigma(1), "Oxygen keeps the DREIDING row")
      if (allocated(error)) return

      call alone%add(dreiding, err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      call alone%build(mol, err)
      call check(error, allocated(err), "An uncovered atom must fail the set build")
      if (allocated(error)) return
      call check(error, index(err%message, "No pair potential parameters for atoms: Ar4") > 0, &
         & "Error should name the atom: "//err%message)
   end subroutine check_element_fallback

   !> Check inactive-site mixing under every rule
   !>
   !> - epsilon = 0, unvalidated negative sigma
   !> - Zero pair potential and derivative, no NaN
   subroutine check_inactive_site_mixing(error)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      type(moist_error), allocatable :: err
      type(lj_12_6_type) :: inactive
      real(wp), parameter :: r(3) = [0.5_wp, 2.0_wp, 6.0_wp]
      real(wp) :: sigma_ij, eps_ij, u(3), du(3)
      integer :: rule, side
      integer, parameter :: rules(2) = [lj_mixing_lorentz_berthelot, lj_mixing_geometric]

      call new_lj_12_6(inactive, [-1.0_wp], [0.0_wp], err)
      call check_moist_error(error, err)
      if (allocated(error)) return
      do rule = 1, size(rules)
         do side = 1, 2
            if (side == 1) then
               call lj_12_6_mix(rules(rule), inactive%sigma(1), inactive%epsilon(1), 6.0_wp, 1.0e-4_wp, &
                  & sigma_ij, eps_ij)
            else
               call lj_12_6_mix(rules(rule), 6.0_wp, 1.0e-4_wp, inactive%sigma(1), inactive%epsilon(1), &
                  & sigma_ij, eps_ij)
            end if
            call check(error, sigma_ij == 0.0_wp .and. eps_ij == 0.0_wp, "An inactive pair kept cross parameters")
            if (allocated(error)) return
            call lj_12_6_evaluate(sigma_ij, eps_ij, r, u, du)
            call check(error, .not. any(ieee_is_nan(u)) .and. .not. any(ieee_is_nan(du)), &
               & "An inactive pair produced NaN")
            if (allocated(error)) return
            call check(error, all(u == 0.0_wp) .and. all(du == 0.0_wp), "An inactive pair produced a potential")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_inactive_site_mixing

end module test_moz_potential_lj
