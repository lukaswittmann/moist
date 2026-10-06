!> Per-point regression tests for cavities (iSwiG, NUMSA, DROP)
!>
!> Every case of [[reference_cases]] builds one cavity kind for one structure of
!> `get_test_structures(mols, 5)` with 26 Lebedev points, updates it, and
!> walks the components every cavity shares into a [[cavity_snapshot_type]]
!> (see [[walk_cavity]]). The snapshot is held against the committed data
!> `test/unit/data/cavity_reference.txt`:
!>
!>   * `ngrid` and `owner` exactly
!>   * every real value with `|got - ref| <= reference_tol*max(1, |ref|)`, points
!>     in array order
!>   * which components the cavity allocates, so an array that disappears fails
!>
!> Data layout, one block per case, blocks in table order:
!>
!>     case    <icase> <kind> <label> <nat> <nleb> <ngrid>
!>     present <icase> <one T/F per component, in component_names order>
!>     int     <icase> <name> <n>        followed by n integers
!>     real    <icase> <name> <n>        followed by n es24.16 values
!>     end     <icase>
!>
!> Rank-2 components are stored flattened in array (column-major) order
!>
!> The data is committed data, produced from the cavity implementations as
!> of 2026-09-30 on macOS arm64. Regenerating it is a deliberate act, done with
!> a throwaway patch that dumps the snapshot in the layout above, never
!> something this test can do
module test_cavity_reference
   use mctc_env, only: wp, fatal_error
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use test_helpers, only: get_test_structures
   use moist_cavity_type, only: cavity_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, &
                                 moist_cavity_iswig_parameters_type
   use moist_cavity_numsa, only: cavity_type_numsa, new_cavity_numsa, &
                                 moist_cavity_numsa_parameters_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use moist_radii, only: default_cpcm_radii
   use moist_context, only: moist_context_type, new_context
   implicit none(type, external)
   private

   public :: collect_cavity_reference

   !> Lebedev points per sphere, every kind
   integer, parameter :: reference_nleb = 26

   !> DROP main tolerance
   !>
   !> - derives projection tolerance, weight cutoff and LSF screening threshold
   real(wp), parameter :: drop_tolerance = 1.0e-12_wp

   !> Threshold of the real comparison
   !>
   !> - `|got - ref| <= reference_tol*max(1, |ref|)`: absolute for `|ref| <= 1`,
   !>   relative above
   real(wp), parameter :: reference_tol = 1.0e-10_wp

   !> Number of structures drawn from `get_test_structures`
   integer, parameter :: n_structures = 5

   !> Collection each structure of `get_test_structures(mols, 5)` comes from
   character(len=*), parameter :: collections(n_structures) = [character(len=9) :: &
                                  "MB16-43", "Heavy28", "Amino20x4", "But14diol", "UPU23"]

   !> Short structure tags used in test names
   character(len=*), parameter :: structure_tags(n_structures) = [character(len=9) :: &
                                  "mb16_43", "heavy28", "amino20x4", "but14diol", "upu23"]

   !> Cavity kind iSwiG
   character(len=*), parameter :: kind_iswig = "iswig"
   !> Cavity kind NUMSA
   character(len=*), parameter :: kind_numsa = "numsa"
   !> Cavity kind DROP with the SvdW level set function
   character(len=*), parameter :: kind_drop_svdw = "drop_svdw"
   !> Cavity kind DROP with the CFC level set function
   character(len=*), parameter :: kind_drop_cfc = "drop_cfc"

   !> Number of walked components
   integer, parameter :: n_components = 10
   !> Components walked, in data order; `owner` is the only integer one
   character(len=*), parameter :: component_names(n_components) = [character(len=12) :: &
                                  "owner", "total_area", "total_volume", "xyz", "a", "f", &
                                  "xi0", "normal0", "v", "asph"]
   !> Leading extent of each component, 0 for a scalar
   !>
   !> - unflattens indices in messages
   integer, parameter :: component_rows(n_components) = [1, 0, 0, 3, 1, 1, 1, 3, 1, 1]
   !> Index of the integer component
   integer, parameter :: comp_owner = 1

   !> One reference case: a cavity kind on one test structure
   type :: reference_case_type
      !> Cavity kind
      character(len=12) :: kind
      !> Index into `get_test_structures(mols, 5)`
      integer :: istruct
      !> Expected structure label, `collection:Hill formula`
      character(len=32) :: label
      !> Expected atom count
      integer :: nat
   end type reference_case_type

   !> Number of reference cases
   integer, parameter :: n_reference_cases = 20
   !> Case table, kind-major; every case builds on today's code
   type(reference_case_type), parameter :: reference_cases(n_reference_cases) = [ &
      reference_case_type(kind_iswig, 1, "MB16-43:CH4AlB4FO2S2Si", 16), &
      reference_case_type(kind_iswig, 2, "Heavy28:H6NSb", 8), &
      reference_case_type(kind_iswig, 3, "Amino20x4:C5H10N2O2", 19), &
      reference_case_type(kind_iswig, 4, "But14diol:C4H10O2", 16), &
      reference_case_type(kind_iswig, 5, "UPU23:C18H22N4O14P", 59), &
      reference_case_type(kind_numsa, 1, "MB16-43:CH4AlB4FO2S2Si", 16), &
      reference_case_type(kind_numsa, 2, "Heavy28:H6NSb", 8), &
      reference_case_type(kind_numsa, 3, "Amino20x4:C5H10N2O2", 19), &
      reference_case_type(kind_numsa, 4, "But14diol:C4H10O2", 16), &
      reference_case_type(kind_numsa, 5, "UPU23:C18H22N4O14P", 59), &
      reference_case_type(kind_drop_svdw, 1, "MB16-43:CH4AlB4FO2S2Si", 16), &
      reference_case_type(kind_drop_svdw, 2, "Heavy28:H6NSb", 8), &
      reference_case_type(kind_drop_svdw, 3, "Amino20x4:C5H10N2O2", 19), &
      reference_case_type(kind_drop_svdw, 4, "But14diol:C4H10O2", 16), &
      reference_case_type(kind_drop_svdw, 5, "UPU23:C18H22N4O14P", 59), &
      reference_case_type(kind_drop_cfc, 1, "MB16-43:CH4AlB4FO2S2Si", 16), &
      reference_case_type(kind_drop_cfc, 2, "Heavy28:H6NSb", 8), &
      reference_case_type(kind_drop_cfc, 3, "Amino20x4:C5H10N2O2", 19), &
      reference_case_type(kind_drop_cfc, 4, "But14diol:C4H10O2", 16), &
      reference_case_type(kind_drop_cfc, 5, "UPU23:C18H22N4O14P", 59)]

   !> One walked component
   type :: snapshot_component_type
      !> Whether the cavity allocates it
      logical :: present = .false.
      !> Integer values, flattened (owner only)
      integer, allocatable :: ival(:)
      !> Real values, flattened; a scalar has one entry
      real(wp), allocatable :: rval(:)
   end type snapshot_component_type

   !> The common components of one updated cavity plus the case it came from
   type :: cavity_snapshot_type
      !> Cavity kind
      character(len=12) :: kind = ""
      !> Structure label
      character(len=32) :: label = ""
      !> Number of atoms
      integer :: nat = 0
      !> Lebedev points per sphere
      integer :: nleb = 0
      !> Number of cavity points
      integer :: ngrid = 0
      !> Walked components, in `component_names` order
      type(snapshot_component_type) :: comp(n_components)
   end type cavity_snapshot_type

contains

   !> Register the cavity reference suite: one test per case plus a table check
   !>
   !> @param[out] testsuite  Collected tests
   subroutine collect_cavity_reference(testsuite)
      !> Collected tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      allocate (testsuite(n_reference_cases + 1))
      testsuite(1) = new_unittest("data_matches_table", test_data_matches_table)
      testsuite(2) = new_unittest(case_name(1), test_case_01)
      testsuite(3) = new_unittest(case_name(2), test_case_02)
      testsuite(4) = new_unittest(case_name(3), test_case_03)
      testsuite(5) = new_unittest(case_name(4), test_case_04)
      testsuite(6) = new_unittest(case_name(5), test_case_05)
      testsuite(7) = new_unittest(case_name(6), test_case_06)
      testsuite(8) = new_unittest(case_name(7), test_case_07)
      testsuite(9) = new_unittest(case_name(8), test_case_08)
      testsuite(10) = new_unittest(case_name(9), test_case_09)
      testsuite(11) = new_unittest(case_name(10), test_case_10)
      testsuite(12) = new_unittest(case_name(11), test_case_11)
      testsuite(13) = new_unittest(case_name(12), test_case_12)
      testsuite(14) = new_unittest(case_name(13), test_case_13)
      testsuite(15) = new_unittest(case_name(14), test_case_14)
      testsuite(16) = new_unittest(case_name(15), test_case_15)
      testsuite(17) = new_unittest(case_name(16), test_case_16)
      testsuite(18) = new_unittest(case_name(17), test_case_17)
      testsuite(19) = new_unittest(case_name(18), test_case_18)
      testsuite(20) = new_unittest(case_name(19), test_case_19)
      testsuite(21) = new_unittest(case_name(20), test_case_20)
   end subroutine collect_cavity_reference

   !> Test name of one case, `<kind>_<structure tag>`
   !>
   !> @param[in] icase  Case index
   function case_name(icase) result(name)
      !> Case index
      integer, intent(in) :: icase
      !> Test name
      character(len=:), allocatable :: name

      name = trim(reference_cases(icase)%kind)//"_"// &
             trim(structure_tags(reference_cases(icase)%istruct))
   end function case_name

   !* ================================================================================= *!
   !*                                  Test entry points                                *!
   !* ================================================================================= *!

   !> Build one case and hold it against its data block
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  icase  Case index
   subroutine run_case(error, icase)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Case index
      integer, intent(in) :: icase

      type(cavity_snapshot_type) :: got, ref
      character(len=:), allocatable :: path

      call reference_data_path(path)
      call take_case_snapshot(icase, got, error)
      if (allocated(error)) return
      call load_reference_case(path, icase, ref, error)
      if (allocated(error)) return
      call compare_snapshot(icase, got, ref, error)
   end subroutine run_case

   !> Data cases match the case table
   !>
   !> - every table case is in the data, in order, with the table's kind,
   !>   label, atom count and Lebedev count
   !> - the data holds no other case
   subroutine test_data_matches_table(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      character(len=:), allocatable :: path
      character(len=256) :: line, msg
      character(len=12) :: key, kind
      character(len=32) :: label
      integer :: unit, stat, icase, nat, nleb, ngrid, nseen

      call reference_data_path(path)
      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "missing reference data '"//path// &
                          "' - it is committed data, restore it from git")
         return
      end if

      nseen = 0
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (.not. has_keyword(line, "case")) cycle
         read (line, *, iostat=stat) key, icase, kind, label, nat, nleb, ngrid
         if (stat /= 0) then
            call test_failed(error, "malformed case line in "//path//": "//trim(line))
            exit
         end if
         nseen = nseen + 1
         if (nseen > n_reference_cases) then
            write (msg, "(a,i0,a,i0)") "data holds more cases than the table: case ", &
               icase, " beyond ", n_reference_cases
            call test_failed(error, trim(msg))
            exit
         end if
         if (icase /= nseen .or. kind /= reference_cases(nseen)%kind &
             .or. label /= reference_cases(nseen)%label .or. nat /= reference_cases(nseen)%nat &
             .or. nleb /= reference_nleb) then
            write (msg, "(a,i0,a,i0,3a,2(a,i0),a)") "data block ", nseen, " is case ", icase, &
               " ", trim(kind), " "//trim(label), " nat ", nat, " nleb ", nleb, &
               ", the table expects "
            write (line, "(i0,3a,2(a,i0))") nseen, " ", trim(reference_cases(nseen)%kind), &
               " "//trim(reference_cases(nseen)%label), " nat ", reference_cases(nseen)%nat, &
               " nleb ", reference_nleb
            call test_failed(error, trim(msg)//"case "//trim(line))
            exit
         end if
      end do
      close (unit)
      if (allocated(error)) return

      if (nseen /= n_reference_cases) then
         write (msg, "(a,i0,a,i0,a)") "data holds ", nseen, " cases, the table ", &
            n_reference_cases, " - regenerate deliberately if the table changed"
         call test_failed(error, trim(msg))
      end if
   end subroutine test_data_matches_table

   !> Case 1, `iswig_mb16_43`, against its data block
   subroutine test_case_01(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1)
   end subroutine test_case_01

   !> Case 2, `iswig_heavy28`, against its data block
   subroutine test_case_02(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2)
   end subroutine test_case_02

   !> Case 3, `iswig_amino20x4`, against its data block
   subroutine test_case_03(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3)
   end subroutine test_case_03

   !> Case 4, `iswig_but14diol`, against its data block
   subroutine test_case_04(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4)
   end subroutine test_case_04

   !> Case 5, `iswig_upu23`, against its data block
   subroutine test_case_05(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5)
   end subroutine test_case_05

   !> Case 6, `numsa_mb16_43`, against its data block
   subroutine test_case_06(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6)
   end subroutine test_case_06

   !> Case 7, `numsa_heavy28`, against its data block
   subroutine test_case_07(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7)
   end subroutine test_case_07

   !> Case 8, `numsa_amino20x4`, against its data block
   subroutine test_case_08(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8)
   end subroutine test_case_08

   !> Case 9, `numsa_but14diol`, against its data block
   subroutine test_case_09(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 9)
   end subroutine test_case_09

   !> Case 10, `numsa_upu23`, against its data block
   subroutine test_case_10(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 10)
   end subroutine test_case_10

   !> Case 11, `drop_svdw_mb16_43`, against its data block
   subroutine test_case_11(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 11)
   end subroutine test_case_11

   !> Case 12, `drop_svdw_heavy28`, against its data block
   subroutine test_case_12(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 12)
   end subroutine test_case_12

   !> Case 13, `drop_svdw_amino20x4`, against its data block
   subroutine test_case_13(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 13)
   end subroutine test_case_13

   !> Case 14, `drop_svdw_but14diol`, against its data block
   subroutine test_case_14(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 14)
   end subroutine test_case_14

   !> Case 15, `drop_svdw_upu23`, against its data block
   subroutine test_case_15(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 15)
   end subroutine test_case_15

   !> Case 16, `drop_cfc_mb16_43`, against its data block
   subroutine test_case_16(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 16)
   end subroutine test_case_16

   !> Case 17, `drop_cfc_heavy28`, against its data block
   subroutine test_case_17(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 17)
   end subroutine test_case_17

   !> Case 18, `drop_cfc_amino20x4`, against its data block
   subroutine test_case_18(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 18)
   end subroutine test_case_18

   !> Case 19, `drop_cfc_but14diol`, against its data block
   subroutine test_case_19(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 19)
   end subroutine test_case_19

   !> Case 20, `drop_cfc_upu23`, against its data block
   subroutine test_case_20(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 20)
   end subroutine test_case_20

   !* ================================================================================= *!
   !*                               Build and walk                                      *!
   !* ================================================================================= *!

   !> Build, update and walk one case
   !>
   !> @param[in]  icase  Case index
   !> @param[out] snap   Walked components plus the case description
   !> @param[out] error  testdrive failure (build or update error, label mismatch)
   subroutine take_case_snapshot(icase, snap, error)
      !> Case index
      integer, intent(in) :: icase
      !> Walked components plus the case description
      type(cavity_snapshot_type), intent(out) :: snap
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(moist_context_type), target :: ctx
      class(cavity_type), allocatable :: cavity
      type(mctc_error), allocatable :: err
      type(reference_case_type) :: gcase
      character(len=:), allocatable :: label
      character(len=256) :: msg

      call new_context(ctx, nthreads=0, verbosity=0)
      call get_test_structures(mols, n_structures)
      gcase = reference_cases(icase)

      associate (mol => mols(gcase%istruct))
         call structure_label(gcase%istruct, mol, label)
         if (label /= gcase%label .or. mol%nat /= gcase%nat) then
            write (msg, "(a,i0,a,a,a,i0,a,a,a,i0,a)") "case ", icase, ": structure is '", &
               trim(label), "' (", mol%nat, " atoms), the table expects '", trim(gcase%label), &
               "' (", gcase%nat, " atoms) - get_test_structures changed?"
            call test_failed(error, trim(msg))
            return
         end if

         call build_case(gcase%kind, ctx, cavity, err)
         if (.not. allocated(err)) call cavity%update(mol, err)
         if (allocated(err)) then
            write (msg, "(a,i0,4a)") "case ", icase, " (", trim(gcase%kind), " ", trim(label)
            call test_failed(error, trim(msg)//"): "//trim(err%message))
            return
         end if

         snap%kind = gcase%kind
         snap%label = label
         snap%nat = mol%nat
         snap%nleb = reference_nleb
         call walk_cavity(cavity, snap)
      end associate
   end subroutine take_case_snapshot

   !> Construct one cavity kind into a polymorphic allocatable
   !>
   !> @param[in]  kind    Cavity kind
   !> @param[in]  ctx     Run context borrowed by the cavity
   !> @param[out] cavity  Constructed, not yet updated cavity
   !> @param[out] err     Construction error
   subroutine build_case(kind, ctx, cavity, err)
      !> Cavity kind
      character(len=*), intent(in) :: kind
      !> Run context borrowed by the cavity
      type(moist_context_type), intent(in), target :: ctx
      !> Constructed, not yet updated cavity
      class(cavity_type), allocatable, intent(out) :: cavity
      !> Construction error
      type(mctc_error), allocatable, intent(out) :: err

      type(cavity_type_iswig), allocatable :: iswig
      type(cavity_type_numsa), allocatable :: numsa
      type(cavity_type_drop), allocatable :: drop
      type(moist_cavity_drop_lsf_svdw_type) :: svdw
      type(moist_cavity_drop_lsf_cfc_type) :: cfc
      type(moist_cavity_drop_parameters_type) :: drop_param

      drop_param = moist_cavity_drop_parameters_type(num_leb=reference_nleb, tolerance=drop_tolerance)

      select case (kind)
      case (kind_iswig)
         allocate (iswig)
         call new_cavity_iswig(iswig, default_cpcm_radii(), err, &
                               param=moist_cavity_iswig_parameters_type(num_leb=reference_nleb), ctx=ctx)
         call move_alloc(iswig, cavity)
      case (kind_numsa)
         allocate (numsa)
         call new_cavity_numsa(numsa, default_cpcm_radii(), err, &
                               param=moist_cavity_numsa_parameters_type(num_leb=reference_nleb), ctx=ctx)
         call move_alloc(numsa, cavity)
      case (kind_drop_svdw)
         allocate (drop)
         call svdw%new()
         call new_cavity_drop(drop, default_cpcm_radii(), svdw, err, param=drop_param, ctx=ctx)
         call move_alloc(drop, cavity)
      case (kind_drop_cfc)
         allocate (drop)
         call cfc%new()
         call new_cavity_drop(drop, default_cpcm_radii(), cfc, err, param=drop_param, ctx=ctx)
         call move_alloc(drop, cavity)
      case default
         call fatal_error(err, "unknown cavity kind '"//kind//"'")
      end select
   end subroutine build_case

   !> Copy the components every cavity shares into a snapshot
   !>
   !> @param[in]     cavity  Updated cavity
   !> @param[in,out] snap    Snapshot receiving `ngrid` and the components
   subroutine walk_cavity(cavity, snap)
      !> Updated cavity
      class(cavity_type), intent(in) :: cavity
      !> Snapshot receiving the components
      type(cavity_snapshot_type), intent(inout) :: snap

      snap%ngrid = cavity%ngrid
      if (allocated(cavity%owner)) then
         snap%comp(comp_owner)%present = .true.
         snap%comp(comp_owner)%ival = cavity%owner
      end if
      call take_scalar(cavity%total_area, snap%comp(2))
      call take_scalar(cavity%total_volume, snap%comp(3))
      call take_rank2(cavity%xyz, snap%comp(4))
      call take_rank1(cavity%a, snap%comp(5))
      call take_rank1(cavity%f, snap%comp(6))
      call take_rank1(cavity%xi0, snap%comp(7))
      call take_rank2(cavity%normal0, snap%comp(8))
      call take_rank1(cavity%v, snap%comp(9))
      call take_rank1(cavity%asph, snap%comp(10))
   end subroutine walk_cavity

   !> Copy an allocatable scalar component
   !>
   !> @param[in]     src   Cavity component
   !> @param[in,out] comp  Snapshot slot
   subroutine take_scalar(src, comp)
      !> Cavity component
      real(wp), allocatable, intent(in) :: src
      !> Snapshot slot
      type(snapshot_component_type), intent(inout) :: comp

      if (.not. allocated(src)) return
      comp%present = .true.
      comp%rval = [src]
   end subroutine take_scalar

   !> Copy a rank-1 component
   !>
   !> @param[in]     src   Cavity component
   !> @param[in,out] comp  Snapshot slot
   subroutine take_rank1(src, comp)
      !> Cavity component
      real(wp), allocatable, intent(in) :: src(:)
      !> Snapshot slot
      type(snapshot_component_type), intent(inout) :: comp

      if (.not. allocated(src)) return
      comp%present = .true.
      comp%rval = src
   end subroutine take_rank1

   !> Copy a rank-2 component, flattened in array order
   !>
   !> @param[in]     src   Cavity component
   !> @param[in,out] comp  Snapshot slot
   subroutine take_rank2(src, comp)
      !> Cavity component
      real(wp), allocatable, intent(in) :: src(:, :)
      !> Snapshot slot
      type(snapshot_component_type), intent(inout) :: comp

      if (.not. allocated(src)) return
      comp%present = .true.
      comp%rval = reshape(src, [size(src)])
   end subroutine take_rank2

   !> Structure label `collection:Hill formula`
   !>
   !> - C, H first when carbon is present, then alphabetical
   !>
   !> @param[in]  istruct  Index into `get_test_structures(mols, 5)`
   !> @param[in]  mol      Structure
   !> @param[out] label    Label
   subroutine structure_label(istruct, mol, label)
      !> Index into `get_test_structures(mols, 5)`
      integer, intent(in) :: istruct
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Label
      character(len=:), allocatable, intent(out) :: label

      character(len=4), allocatable :: sym(:)
      integer, allocatable :: cnt(:)
      logical, allocatable :: done(:)
      character(len=16) :: part
      integer :: i, j, pick

      allocate (sym(mol%nid), cnt(mol%nid), done(mol%nid))
      do i = 1, mol%nid
         sym(i) = mol%sym(i)
         cnt(i) = count(mol%id == i)
      end do
      done = .false.

      label = trim(collections(istruct))//":"
      if (any(sym == "C")) then
         do i = 1, mol%nid
            if (sym(i) == "C") call append(i)
         end do
         do i = 1, mol%nid
            if (sym(i) == "H") call append(i)
         end do
      end if
      do
         pick = 0
         do j = 1, mol%nid
            if (done(j)) cycle
            if (pick == 0) then
               pick = j
            else if (llt(sym(j), sym(pick))) then
               pick = j
            end if
         end do
         if (pick == 0) exit
         call append(pick)
      end do

   contains

      !> Append one element and its count to the label
      !>
      !> @param[in] k  Species index
      subroutine append(k)
         !> Species index
         integer, intent(in) :: k

         done(k) = .true.
         if (cnt(k) == 1) then
            part = trim(sym(k))
         else
            write (part, "(a,i0)") trim(sym(k)), cnt(k)
         end if
         label = label//trim(part)
      end subroutine append
   end subroutine structure_label

   !* ================================================================================= *!
   !*                                  Data reader                                   *!
   !* ================================================================================= *!

   !> Path of the committed data
   !>
   !> - a subroutine: a deferred-length function result is shared between the
   !>   threads test-drive runs the cases on
   !>
   !> @param[out] path  Full path
   subroutine reference_data_path(path)
      !> Full path
      character(len=:), allocatable, intent(out) :: path

      character(len=4096) :: root
      integer :: length, stat

      call get_environment_variable("MOIST_SOURCE_ROOT", root, length, stat)
      if (stat /= 0 .or. length == 0) root = "."
      path = trim(root)//"/test/unit/data/cavity_reference.txt"
   end subroutine reference_data_path

   !> Read the block of one case from the data
   !>
   !> @param[in]  path   Data path
   !> @param[in]  icase  Case index
   !> @param[out] snap   Snapshot as recorded
   !> @param[out] error  testdrive failure (missing file or case, malformed block)
   subroutine load_reference_case(path, icase, snap, error)
      !> Data path
      character(len=*), intent(in) :: path
      !> Case index
      integer, intent(in) :: icase
      !> Snapshot as recorded
      type(cavity_snapshot_type), intent(out) :: snap
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      character(len=256) :: line, msg
      character(len=12) :: key
      integer :: unit, stat, idx
      logical :: found

      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "missing reference data '"//path// &
                          "' - it is committed data, restore it from git")
         return
      end if

      found = .false.
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (.not. has_keyword(line, "case")) cycle
         read (line, *, iostat=stat) key, idx
         if (stat /= 0 .or. idx /= icase) cycle
         read (line, *, iostat=stat) key, idx, snap%kind, snap%label, snap%nat, snap%nleb, &
            snap%ngrid
         if (stat /= 0) exit
         found = .true.
         call read_case_blocks(unit, icase, path, snap, error)
         exit
      end do
      close (unit)
      if (allocated(error)) return

      if (.not. found) then
         write (msg, "(a,i0,a)") "case ", icase, " not found in reference data '"
         call test_failed(error, trim(msg)//path//"'")
      end if
   end subroutine load_reference_case

   !> Read the `present`, `int` and `real` blocks of one case up to its `end`
   !>
   !> @param[in]     unit   Open data, positioned after the case line
   !> @param[in]     icase  Case index
   !> @param[in]     path   Data path, for messages
   !> @param[in,out] snap   Snapshot receiving the components
   !> @param[out]    error  testdrive failure
   subroutine read_case_blocks(unit, icase, path, snap, error)
      !> Open data
      integer, intent(in) :: unit
      !> Case index
      integer, intent(in) :: icase
      !> Data path
      character(len=*), intent(in) :: path
      !> Snapshot receiving the components
      type(cavity_snapshot_type), intent(inout) :: snap
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      character(len=256) :: line
      character(len=12) :: key, name
      character(len=64) :: where
      logical :: flags(n_components), seen(n_components)
      integer :: stat, idx, n, ic

      write (where, "(a,i0,a)") "case ", icase, " of reference data '"
      flags = .false.
      seen = .false.
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) then
            call test_failed(error, trim(where)//path//"' ends before its 'end' line")
            return
         end if
         if (has_keyword(line, "end")) exit
         if (has_keyword(line, "present")) then
            read (line, *, iostat=stat) key, idx, flags
         else if (has_keyword(line, "int") .or. has_keyword(line, "real")) then
            read (line, *, iostat=stat) key, idx, name, n
            if (stat == 0) then
               ic = findloc(component_names, name, dim=1)
               if (ic == 0) then
                  call test_failed(error, trim(where)//path//"' names unknown component '"// &
                                   trim(name)//"'")
                  return
               end if
               seen(ic) = .true.
               if (ic == comp_owner) then
                  allocate (snap%comp(ic)%ival(n))
                  if (n > 0) read (unit, *, iostat=stat) snap%comp(ic)%ival
               else
                  allocate (snap%comp(ic)%rval(n))
                  if (n > 0) read (unit, *, iostat=stat) snap%comp(ic)%rval
               end if
            end if
         else
            stat = 1
         end if
         if (stat /= 0) then
            call test_failed(error, trim(where)//path//"' is malformed at: "//trim(line))
            return
         end if
      end do

      if (any(flags .neqv. seen)) then
         call test_failed(error, trim(where)//path// &
                          "' has a 'present' line that does not match its blocks")
         return
      end if
      snap%comp(:)%present = flags
   end subroutine read_case_blocks

   !> Whether a data line starts with a keyword
   !>
   !> @param[in] line  Raw line
   !> @param[in] key   Keyword
   pure function has_keyword(line, key) result(match)
      !> Raw line
      character(len=*), intent(in) :: line
      !> Keyword
      character(len=*), intent(in) :: key
      !> Whether the first token is `key`
      logical :: match

      integer :: n

      n = len(key)
      match = .false.
      if (len(line) <= n) return
      match = line(1:n) == key .and. line(n + 1:n + 1) == " "
   end function has_keyword

   !* ================================================================================= *!
   !*                                   Comparison                                      *!
   !* ================================================================================= *!

   !> Hold a fresh snapshot against the recorded one
   !>
   !> Reports every component that deviates, each with its first failing
   !> index, both values, the number of failing entries, and the index of the
   !> largest deviation
   !>
   !> @param[in]  icase  Case index
   !> @param[in]  got    Fresh snapshot
   !> @param[in]  ref    Recorded snapshot
   !> @param[out] error  testdrive failure
   subroutine compare_snapshot(icase, got, ref, error)
      !> Case index
      integer, intent(in) :: icase
      !> Fresh snapshot
      type(cavity_snapshot_type), intent(in) :: got
      !> Recorded snapshot
      type(cavity_snapshot_type), intent(in) :: ref
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      character(len=:), allocatable :: report, problem
      character(len=256) :: head, msg
      integer :: ic, nbad

      write (head, "(a,i0,5a)") "case ", icase, " (", trim(ref%kind), " ", trim(ref%label), ")"
      report = ""
      nbad = 0

      if (got%kind /= ref%kind .or. got%label /= ref%label .or. got%nat /= ref%nat &
          .or. got%nleb /= ref%nleb) then
         write (msg, "(4a,2(1x,i0),a,2a,2(1x,i0))") "data describes ", trim(ref%kind), &
            " ", trim(ref%label), ref%nat, ref%nleb, ", the build is ", trim(got%kind)//" ", &
            trim(got%label), got%nat, got%nleb
         call test_failed(error, trim(head)//": "//trim(msg))
         return
      end if

      if (got%ngrid /= ref%ngrid) then
         write (msg, "(a,i0,a,i0)") "ngrid reference ", ref%ngrid, " now ", got%ngrid
         call add_problem(msg)
      end if

      do ic = 1, n_components
         if (ref%comp(ic)%present .neqv. got%comp(ic)%present) then
            if (ref%comp(ic)%present) then
               msg = trim(component_names(ic))//" is allocated in the data but not now"
            else
               msg = trim(component_names(ic))//" is allocated now but not in the data"
            end if
            call add_problem(msg)
            cycle
         end if
         if (.not. ref%comp(ic)%present) cycle
         if (ic == comp_owner) then
            call compare_int(ic, got%comp(ic)%ival, ref%comp(ic)%ival, problem)
         else
            call compare_real(ic, got%comp(ic)%rval, ref%comp(ic)%rval, problem)
         end if
         if (len(problem) > 0) call add_problem(problem)
      end do

      if (nbad == 0) return
      call test_failed(error, trim(head)//": "//report)

   contains

      !> Append one problem to the report
      !>
      !> @param[in] text  Problem description
      subroutine add_problem(text)
         !> Problem description
         character(len=*), intent(in) :: text

         nbad = nbad + 1
         if (nbad > 1) report = report//"; "
         report = report//trim(text)
      end subroutine add_problem
   end subroutine compare_snapshot

   !> Exact comparison of an integer component
   !>
   !> @param[in]  ic   Component index
   !> @param[in]  got  Fresh values
   !> @param[in]  ref  Recorded values
   !> @param[out] msg  Problem description, blank when all entries match
   subroutine compare_int(ic, got, ref, msg)
      !> Component index
      integer, intent(in) :: ic
      !> Fresh values
      integer, intent(in) :: got(:)
      !> Recorded values
      integer, intent(in) :: ref(:)
      !> Problem description
      character(len=:), allocatable, intent(out) :: msg

      character(len=:), allocatable :: at
      character(len=256) :: buf
      integer :: i, nfail, first

      msg = ""
      if (size(got) /= size(ref)) then
         write (buf, "(a,a,i0,a,i0)") trim(component_names(ic)), " size reference ", size(ref), &
            " now ", size(got)
         msg = trim(buf)
         return
      end if
      nfail = count(got /= ref)
      if (nfail == 0) return
      first = 0
      do i = 1, size(ref)
         if (got(i) /= ref(i)) then
            first = i
            exit
         end if
      end do
      call index_text(ic, first, at)
      write (buf, "(a,a,i0,a,i0,a,i0,a)") at, " reference ", ref(first), " now ", got(first), &
         " (", nfail, " entries differ)"
      msg = trim(buf)
   end subroutine compare_int

   !> Per-entry comparison of a real component
   !>
   !> - an entry passes when `|got - ref| <= reference_tol*max(1, |ref|)`
   !> - the largest failure is ranked by `|got - ref|/max(1, |ref|)`
   !>
   !> @param[in]  ic   Component index
   !> @param[in]  got  Fresh values
   !> @param[in]  ref  Recorded values
   !> @param[out] msg  Problem description, blank when all entries match
   subroutine compare_real(ic, got, ref, msg)
      !> Component index
      integer, intent(in) :: ic
      !> Fresh values
      real(wp), intent(in) :: got(:)
      !> Recorded values
      real(wp), intent(in) :: ref(:)
      !> Problem description
      character(len=:), allocatable, intent(out) :: msg

      character(len=:), allocatable :: at, at_max
      character(len=256) :: buf
      real(wp) :: diff, scaled, smax
      integer :: i, nfail, first, imax

      msg = ""
      if (size(got) /= size(ref)) then
         write (buf, "(a,a,i0,a,i0)") trim(component_names(ic)), " size reference ", size(ref), &
            " now ", size(got)
         msg = trim(buf)
         return
      end if
      nfail = 0
      first = 0
      imax = 0
      smax = -1.0_wp
      do i = 1, size(ref)
         diff = abs(got(i) - ref(i))
         ! written so that a NaN on either side counts as a failure
         if (.not. (diff <= reference_tol*max(1.0_wp, abs(ref(i))))) then
            nfail = nfail + 1
            if (first == 0) first = i
            scaled = diff/max(1.0_wp, abs(ref(i)))
            if (.not. (scaled <= smax)) then
               smax = scaled
               imax = i
            end if
         end if
      end do
      if (nfail == 0) return
      call index_text(ic, first, at)
      call index_text(ic, imax, at_max)
      write (buf, "(a,a,es24.16,a,es24.16,a,es9.2,a,i0,a,es9.2,a,a,a)") at, &
         " reference", ref(first), " now", got(first), " (|now - reference|", &
         abs(got(first) - ref(first)), "; ", nfail, " entries above", reference_tol, &
         "*max(1, |reference|), largest at ", at_max, ")"
      msg = trim(buf)
   end subroutine compare_real

   !> Human-readable index of a flattened component entry, e.g. `xyz(2, 17)`
   !>
   !> @param[in]  ic    Component index
   !> @param[in]  k     Flattened index (1-based)
   !> @param[out] text  Rendered index
   subroutine index_text(ic, k, text)
      !> Component index
      integer, intent(in) :: ic
      !> Flattened index
      integer, intent(in) :: k
      !> Rendered index
      character(len=:), allocatable, intent(out) :: text

      character(len=64) :: buf
      integer :: rows

      rows = component_rows(ic)
      if (rows == 0) then
         buf = component_names(ic)
      else if (rows == 1) then
         write (buf, "(a,a,i0,a)") trim(component_names(ic)), "(", k, ")"
      else
         write (buf, "(a,a,i0,a,i0,a)") trim(component_names(ic)), "(", mod(k - 1, rows) + 1, &
            ", ", (k - 1)/rows + 1, ")"
      end if
      text = trim(buf)
   end subroutine index_text

end module test_cavity_reference
