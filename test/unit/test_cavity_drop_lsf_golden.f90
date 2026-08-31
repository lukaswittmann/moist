!> Golden numerical fixture for the DROP level set functions (SvdW and CFC)
!>
!> One traversal walks a fixed set of structures, evaluation points and
!> screening thresholds, and turns every quantity the LSF API exposes into a
!> stream of labelled records. That stream is then either compared against the
!> committed fixture (`test/unit/data/lsf_golden_*.txt`) or against a second
!> traversal in the same process.
!>
!> Record layout, one per line of the fixture:
!>
!>     kind case ip flag quantity i1 i2 i3 i4 i5 i6 value
!>
!>   * `kind`     `svdw` or `cfc`
!>   * `case`     tag of a [[golden_cases]] entry
!>   * `ip`       evaluation-point index, see [[build_points]]
!>   * `flag`     `S` screened as production runs it, `U` unscreened
!>                reference, `D` screened minus unscreened
!>   * `i1..i6`   index tuple, unused slots 0; `A`/`B` are *user-space* atom
!>                indices, so the fixture survives an index-space change in
!>                `src/` (see [[nuc_slot]])
!>   * `value`    `es24.16`, an exact IEEE-754 double round trip
!>
!> Slots that are symmetric by construction are dumped once (j <= k <= l <= m);
!> the discarded components are covered by the `*_tensor_symmetry` tests.
!>
!> The fixture is committed data and this module deliberately cannot write it:
!> a golden reference the test can rewrite is one keystroke away from being
!> "fixed" instead of investigated. If the LSF definition legitimately changes,
!> regenerate deliberately - dump [[golden_stream_type]] from a throwaway patch
!> - and review the numerical diff before committing it.
module test_cavity_drop_lsf_golden
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use moist_utils_env, only: get_env
   use test_helpers, only: get_test_radii, check_moist_error, rel_deviation
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private

   public :: collect_cavity_drop_lsf_golden

   integer, parameter :: ndim = 3

   !> Concrete selectors
   character(len=*), parameter :: kind_svdw = "svdw"
   character(len=*), parameter :: kind_cfc = "cfc"

   !> Highest derivative order each concrete is driven at
   integer, parameter :: max_deriv_svdw = 4
   integer, parameter :: max_deriv_cfc = 3

   !> Screening threshold production uses (DROP `parameters.f90` default)
   real(wp), parameter :: production_threshold = 1.0e-11_wp

   !> Tolerance of the golden comparison, measured as `rel_deviation`,
   !> `|a - b| / (1 + |b|)`: relative for large values, absolute near zero
   real(wp), parameter :: golden_tol = 1.0e-12_wp

   !> Tolerances of the companion symmetry checks (gross-asymmetry guard, not
   !> a roundoff assertion); relative to the largest element of the tensor
   real(wp), parameter :: symmetry_rel_tol = 1.0e-8_wp
   real(wp), parameter :: symmetry_abs_tol = 1.0e-12_wp

   !> Number of evaluation points per case
   integer, parameter :: n_points = 5
   !> Atoms carrying single-nucleus derivative records
   integer, parameter :: n_sel_atoms = 3
   !> Atoms carrying pair (two-nucleus) derivative records
   integer, parameter :: n_pair_atoms = 2

   !> Unnormalised offset directions, one per evaluation point. Deliberately
   !> off-axis and mutually non-parallel so no point lands on a symmetry
   !> element of a symmetric fixture (which could put it on a nucleus).
   real(wp), parameter :: raw_dirs(ndim, n_points) = reshape([ &
                          1.0_wp, 2.0_wp, 3.0_wp, &
                          -2.0_wp, 1.0_wp, 4.0_wp, &
                          3.0_wp, -4.0_wp, 1.0_wp, &
                          1.0_wp, -1.0_wp, 2.0_wp, &
                          -1.0_wp, 3.0_wp, -2.0_wp], [ndim, n_points])

   !> Radial offsets applied to `raw_dirs`; see [[build_points]] for how each
   !> is anchored. Point 5 is the deliberate near-nucleus probe.
   real(wp), parameter :: point_offsets(n_points) = &
                          [13.0_wp, 1.25_wp, 0.10_wp, 0.40_wp, 0.05_wp]

   !> Smallest distance from any evaluation point to any nucleus that
   !> [[build_points]] tolerates. Several kernels guard on `x > 0`; pinning
   !> the exactly-on-a-nucleus branch is not the intent here.
   real(wp), parameter :: min_nucleus_clearance = 4.0e-2_wp

   !> One fixture case: an mstore structure plus the SvdW blending weights it
   !> is evaluated with. CFC has no analogous knob and ignores the weights.
   type :: golden_case_type
      !> Short tag written into every record of this case
      character(len=12) :: tag
      !> mstore collection
      character(len=12) :: collection
      !> mstore record id inside the collection
      character(len=12) :: record
      !> SvdW blending sharpness
      real(wp) :: blend_k
      !> SvdW one-body weight
      real(wp) :: blend_1b
      !> SvdW two-body weight
      real(wp) :: blend_2b
      !> SvdW three-body weight
      real(wp) :: blend_3b
   end type golden_case_type

   !> Structures shared by both concretes. Chosen for spread rather than
   !> chemistry: a diatomic, a small symmetric hydride, a heavy-element
   !> complex (large radii, so the blending shells overlap differently), a
   !> 16-atom mindless cage, and an amino acid. Cases 1-5 are also the CFC
   !> cases; case 6 repeats CH4 with the legacy SvdW weights, because the
   !> shipped defaults set `blend_2b = 0` and would otherwise leave the whole
   !> two-body branch of the blending unpinned.
   integer, parameter :: n_svdw_cases = 6
   integer, parameter :: n_cfc_cases = 5
   type(golden_case_type), parameter :: golden_cases(n_svdw_cases) = [ &
      golden_case_type("lih", "MB16-43", "LiH", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      golden_case_type("ch4", "MB16-43", "CH4", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      golden_case_type("bih3_h2o", "Heavy28", "bih3_h2o", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      golden_case_type("mb16_01", "MB16-43", "01", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      golden_case_type("ala_xab", "Amino20x4", "ALA_xab", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      golden_case_type("ch4_legacy", "MB16-43", "CH4", 3.0_wp, 1.0_wp, 1.0_wp, 1.0_wp)]

   !* ================================================================================= *!
   !*                                  Record stream                                    *!
   !* ================================================================================= *!

   !> What a record is about: the concrete, the case, the evaluation point and
   !> the screening flag. Every record of one block shares one of these, so it
   !> travels through the emission routines as a single argument.
   type :: record_id_type
      !> `svdw` or `cfc`
      character(len=8) :: kind
      !> Case tag
      character(len=12) :: case_tag
      !> Evaluation-point index
      integer :: ip
      !> `S`, `U` or `D`
      character(len=1) :: flag
   end type record_id_type

   !> One record: what identifies a number, plus the number
   type :: golden_record_type
      !> Concrete, case, evaluation point and screening flag
      type(record_id_type) :: id
      !> Quantity name
      character(len=20) :: quantity
      !> Index tuple, unused slots 0
      integer :: idx(6)
      !> Value
      real(wp) :: val
   end type golden_record_type

   !> The record stream one traversal produces
   !>
   !> Structural problems - a shape, an index space or an invariant that no
   !> longer matches what this harness understands - are collected here rather
   !> than raised on the spot, so the traversal always runs to completion and
   !> the caller reports them through the same channel as a numerical mismatch.
   type :: golden_stream_type
      !> Emitted records, valid entries `1:n`
      type(golden_record_type), allocatable :: rec(:)
      !> Number of emitted records
      integer :: n = 0
      !> Number of structural problems seen
      integer :: nproblem = 0
      !> Description of the first structural problem
      character(len=:), allocatable :: problem
   contains
      !> Append one record
      procedure :: emit => stream_emit
      !> Record a structural problem
      procedure :: flag_problem => stream_flag_problem
   end type golden_stream_type

contains

   !> Register the golden-fixture suite
   subroutine collect_cavity_drop_lsf_golden(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_matches_golden", test_svdw_golden), &
                  new_unittest("cfc_matches_golden", test_cfc_golden), &
                  new_unittest("svdw_tensor_symmetry", test_svdw_symmetry), &
                  new_unittest("cfc_tensor_symmetry", test_cfc_symmetry), &
                  new_unittest("golden_stream_is_reproducible", test_stream_reproducible) &
                  ]
   end subroutine collect_cavity_drop_lsf_golden

   !* ================================================================================= *!
   !*                              Suite entry points                                   *!
   !* ================================================================================= *!

   !> Compare the SvdW LSF against `lsf_golden_svdw.txt`
   subroutine test_svdw_golden(error)
      type(error_type), allocatable, intent(out) :: error
      call run_golden(error, kind_svdw)
   end subroutine test_svdw_golden

   !> Compare the CFC LSF against `lsf_golden_cfc.txt`
   subroutine test_cfc_golden(error)
      type(error_type), allocatable, intent(out) :: error
      call run_golden(error, kind_cfc)
   end subroutine test_cfc_golden

   !> Traverse one concrete and hold the resulting stream against its fixture
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine run_golden(error, kind)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      type(golden_stream_type) :: got
      type(golden_record_type), allocatable :: ref(:)
      character(len=:), allocatable :: path

      path = golden_path(kind)

      call traverse(kind, got, error)
      if (allocated(error)) return
      call check_problems(got, error)
      if (allocated(error)) return

      call load_fixture(path, ref, error)
      if (allocated(error)) return

      call compare_stream(got, ref, path, error)
   end subroutine run_golden

   !> Path of a fixture file. meson exports `MOIST_SOURCE_ROOT`; fpm runs the
   !> tester from the project root, which the `.` default covers.
   !>
   !> @param[in] kind  `svdw` or `cfc`
   !> @returns         Full path of the fixture
   function golden_path(kind) result(path)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      character(len=:), allocatable :: path

      path = get_env("MOIST_SOURCE_ROOT", default=".")//"/test/unit/data/lsf_golden_" &
             //kind//".txt"
   end function golden_path

   !* ================================================================================= *!
   !*                                   Traversal                                       *!
   !* ================================================================================= *!

   !> Walk every case / point / screening flag of one concrete, emitting the
   !> full record set in a fixed order.
   !>
   !> @param[in]  kind    `svdw` or `cfc`
   !> @param[out] stream  Emitted records
   !> @param[out] error   testdrive failure (setup problems only)
   subroutine traverse(kind, stream, error)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> Emitted records
      type(golden_stream_type), intent(out) :: stream
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer, allocatable :: sel(:)
      integer :: icase, ip

      do icase = 1, n_cases(kind)
         call load_case(golden_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return
         call select_atoms(mol%nat, sel)

         do ip = 1, n_points
            call emit_point(kind, golden_cases(icase), mol, radii, points(:, ip), ip, &
                            sel, stream, error)
            if (allocated(error)) return
         end do

         deallocate (radii, points, sel)
      end do
   end subroutine traverse

   !> Number of fixture cases one concrete is driven over
   !>
   !> @param[in] kind  `svdw` or `cfc`
   !> @returns         Case count
   pure integer function n_cases(kind)
      !> Concrete selector
      character(len=*), intent(in) :: kind

      if (kind == kind_svdw) then
         n_cases = n_svdw_cases
      else
         n_cases = n_cfc_cases
      end if
   end function n_cases

   !> Fetch an mstore structure and its CPCM-table radii
   !>
   !> @param[in]  gcase  Case descriptor
   !> @param[out] mol    Structure
   !> @param[out] radii  Per-atom radii (size mol%nat)
   subroutine load_case(gcase, mol, radii)
      !> Case descriptor
      type(golden_case_type), intent(in) :: gcase
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Per-atom radii
      real(wp), allocatable, intent(out) :: radii(:)

      call get_structure(mol, trim(gcase%collection), trim(gcase%record))
      call get_test_radii(mol, radii)
   end subroutine load_case

   !> Build the deterministic evaluation points of one structure.
   !>
   !> All five are anchored on nuclei rather than on the bounding box, so
   !> their relation to the surface is the same for every structure:
   !>
   !>   1 `far_out`  atom 1 surface + 13.0 bohr - beyond the SvdW screening
   !>                reach (~13.8 bohr at k = 5.5, threshold 1e-11) for the
   !>                *far* atoms, so screening genuinely bites here
   !>   2 `near_out` atom 1 surface + 1.25 bohr - outside, well conditioned
   !>   3 `surface`  atom 1 surface + 0.10 bohr - near the zero level set
   !>   4 `deep_in`  0.40 * R1 from nucleus 1 - deep inside the cavity
   !>   5 `near_nuc` 0.05 bohr from the last nucleus - close to, but never
   !>                on, a nucleus
   !>
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[out] points (3, n_points) evaluation points
   !> @param[out] error  testdrive failure if a point lands on a nucleus
   subroutine build_points(mol, radii, points, error)
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Evaluation points
      real(wp), allocatable, intent(out) :: points(:, :)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: dir(ndim, n_points)
      real(wp) :: dmin
      integer :: ip, iat
      character(len=64) :: tail

      do ip = 1, n_points
         dir(:, ip) = raw_dirs(:, ip)/norm2(raw_dirs(:, ip))
      end do

      allocate (points(ndim, n_points))
      points(:, 1) = mol%xyz(:, 1) + (radii(1) + point_offsets(1))*dir(:, 1)
      points(:, 2) = mol%xyz(:, 1) + (radii(1) + point_offsets(2))*dir(:, 2)
      points(:, 3) = mol%xyz(:, 1) + (radii(1) + point_offsets(3))*dir(:, 3)
      points(:, 4) = mol%xyz(:, 1) + point_offsets(4)*radii(1)*dir(:, 4)
      points(:, 5) = mol%xyz(:, mol%nat) + point_offsets(5)*dir(:, 5)

      do ip = 1, n_points
         dmin = huge(0.0_wp)
         do iat = 1, mol%nat
            dmin = min(dmin, norm2(points(:, ip) - mol%xyz(:, iat)))
         end do
         if (dmin < min_nucleus_clearance) then
            write (tail, '(i0,a,es12.4)') ip, " sits ", dmin
            call test_failed(error, "evaluation point "//trim(tail)// &
                             " bohr from a nucleus - reference geometry changed?")
            return
         end if
      end do
   end subroutine build_points

   !> Deterministic atom selection: first, middle and last centre, deduplicated
   !> while keeping ascending order. Small structures simply yield fewer.
   !>
   !> @param[in]  nat  Number of atoms
   !> @param[out] sel  Selected user-space atom indices, ascending
   subroutine select_atoms(nat, sel)
      !> Number of atoms
      integer, intent(in) :: nat
      !> Selected atom indices
      integer, allocatable, intent(out) :: sel(:)

      integer :: cand(n_sel_atoms), tmp(n_sel_atoms)
      integer :: i, n

      cand = [1, 1 + (nat - 1)/2, nat]
      n = 0
      do i = 1, n_sel_atoms
         if (n > 0) then
            if (any(tmp(1:n) == cand(i))) cycle
         end if
         n = n + 1
         tmp(n) = cand(i)
      end do
      allocate (sel(n), source=tmp(1:n))
   end subroutine select_atoms

   !* ================================================================================= *!
   !*                                LSF construction                                   *!
   !* ================================================================================= *!

   !> Allocate a fresh LSF of the requested concrete kind, bind it to `mol` and
   !> drive it to that kind's derivative cap.
   !>
   !> The screening threshold has to be set before `update`, which is what
   !> pushes it into the SSD system. Only the constructor differs between the
   !> concretes; everything after it is base-class API.
   !>
   !> @param[out] lsf    Fresh LSF
   !> @param[in]  kind   `svdw` or `cfc`
   !> @param[in]  gcase  Case descriptor (supplies the SvdW blending weights)
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  thr    Screening threshold
   subroutine new_lsf(lsf, kind, gcase, mol, radii, thr)
      !> Fresh LSF
      class(moist_cavity_drop_lsf_type), allocatable, intent(out) :: lsf
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> Case descriptor
      type(golden_case_type), intent(in) :: gcase
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Screening threshold
      real(wp), intent(in) :: thr

      integer :: max_deriv

      select case (kind)
      case (kind_svdw)
         allocate (moist_cavity_drop_lsf_svdw_type :: lsf)
         select type (lsf)
         type is (moist_cavity_drop_lsf_svdw_type)
            lsf%screening_threshold = thr
            call lsf%new(blend_k=gcase%blend_k, blend_1b=gcase%blend_1b, &
                         blend_2b=gcase%blend_2b, blend_3b=gcase%blend_3b)
         end select
         max_deriv = max_deriv_svdw
      case (kind_cfc)
         allocate (moist_cavity_drop_lsf_cfc_type :: lsf)
         select type (lsf)
         type is (moist_cavity_drop_lsf_cfc_type)
            lsf%screening_threshold = thr
            call lsf%new()
         end select
         max_deriv = max_deriv_cfc
      case default
         error stop "new_lsf: unknown kind '"//kind//"'"
      end select

      call lsf%update(mol, radii)
      call lsf%set_max_deriv(max_deriv)
   end subroutine new_lsf

   !> `prepare` an LSF, translating an LSF error into a testdrive failure
   !>
   !> @param[inout] lsf    LSF to prepare
   !> @param[in]    point  Evaluation point
   !> @param[out]   error  testdrive failure
   subroutine prepare_lsf(lsf, point, error)
      !> LSF to prepare
      class(moist_cavity_drop_lsf_type), intent(inout) :: lsf
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: err

      call lsf%prepare(point, err)
      call check_moist_error(error, err, "LSF prepare failed")
   end subroutine prepare_lsf

   !* ================================================================================= *!
   !*                                Record emission                                    *!
   !* ================================================================================= *!

   !> Emit every record of one case at one evaluation point: screened block
   !> first, unscreened block second, difference record last.
   !>
   !> @param[in]    kind    `svdw` or `cfc`
   !> @param[in]    gcase   Case descriptor
   !> @param[in]    mol     Structure
   !> @param[in]    radii   Per-atom radii
   !> @param[in]    point   Evaluation point
   !> @param[in]    ip      Evaluation-point index
   !> @param[in]    sel     Selected atom indices
   !> @param[inout] stream  Record sink
   !> @param[out]   error   testdrive failure
   subroutine emit_point(kind, gcase, mol, radii, point, ip, sel, stream, error)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> Case descriptor
      type(golden_case_type), intent(in) :: gcase
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      class(moist_cavity_drop_lsf_type), allocatable :: lsf_scr, lsf_ref
      real(wp) :: f0_scr, f0_ref

      call new_lsf(lsf_scr, kind, gcase, mol, radii, production_threshold)
      call new_lsf(lsf_ref, kind, gcase, mol, radii, 0.0_wp)

      call prepare_lsf(lsf_scr, point, error)
      if (allocated(error)) return
      call prepare_lsf(lsf_ref, point, error)
      if (allocated(error)) return

      call assert_unscreened(lsf_ref%active_count(), mol%nat, gcase%tag, ip, error)
      if (allocated(error)) return

      call emit_block(lsf_scr, record_id_type(kind, gcase%tag, ip, "S"), sel, stream, f0_scr)
      call emit_block(lsf_ref, record_id_type(kind, gcase%tag, ip, "U"), sel, stream, f0_ref)
      call stream%emit(record_id_type(kind, gcase%tag, ip, "D"), "f0_delta", idx6(), &
                       f0_scr - f0_ref)
   end subroutine emit_point

   !> Emit every quantity of one prepared LSF under one screening flag
   !>
   !> The two concretes share the spatial and single-nucleus blocks; the pair
   !> and normalisation blocks, and everything of derivative order 4, are SvdW
   !> only, because CFC is capped at order 3 today and an accessor asked for an
   !> order it was not prepared for is a hard failure by design.
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[in]    sel     Selected atom indices
   !> @param[inout] stream  Record sink
   !> @param[out]   f0      Value returned by `f0`, for the delta record
   subroutine emit_block(lsf, id, sel, stream, f0)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream
      !> Value of `f0`
      real(wp), intent(out) :: f0

      !> Atom -> active index map (0 = screened away)
      integer, allocatable :: act(:)

      call active_map(lsf, act, stream)

      call emit_spatial_block(lsf, id, stream, f0)
      call emit_nuclear_block(lsf, id, sel, act, stream)
      if (id%kind /= kind_svdw) return
      call emit_pair_block(lsf, id, sel, act, stream)
      call emit_normalized_block(lsf, id, sel, act, stream)
   end subroutine emit_block

   !> Emit the quantities that differentiate with respect to the evaluation
   !> point only, plus the active count that scopes every nuclear record
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[inout] stream  Record sink
   !> @param[out]   f0      Value returned by `f0`
   subroutine emit_spatial_block(lsf, id, stream, f0)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream
      !> Value of `f0`
      real(wp), intent(out) :: f0

      real(wp) :: lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim)
      real(wp) :: lsf3_rrr(ndim, ndim, ndim), lsf4_rrrr(ndim, ndim, ndim, ndim)
      integer :: j, k, l, m

      call stream%emit(id, "n_active", idx6(), real(lsf%active_count(), wp))

      call lsf%f0(f0)
      call stream%emit(id, "f0", idx6(), f0)

      call lsf%f012_r(lsf0, lsf1_r, lsf2_rr)
      call stream%emit(id, "f012_lsf0", idx6(), lsf0)
      do j = 1, ndim
         call stream%emit(id, "f012_lsf1_r", idx6(j), lsf1_r(j))
      end do
      do j = 1, ndim
         do k = j, ndim
            call stream%emit(id, "f012_lsf2_rr", idx6(j, k), lsf2_rr(j, k))
         end do
      end do

      call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               call stream%emit(id, "f3_rrr", idx6(j, k, l), lsf3_rrr(j, k, l))
            end do
         end do
      end do

      if (id%kind /= kind_svdw) return

      call lsf%f4_rrrr(lsf4_rrrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               do m = l, ndim
                  call stream%emit(id, "f4_rrrr", idx6(j, k, l, m), lsf4_rrrr(j, k, l, m))
               end do
            end do
         end do
      end do
   end subroutine emit_spatial_block

   !> Emit the derivatives that carry exactly one nuclear index
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[in]    sel     Selected atom indices
   !> @param[in]    act     Atom -> active index map
   !> @param[inout] stream  Record sink
   subroutine emit_nuclear_block(lsf, id, sel, act, stream)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Atom -> active index map
      integer, intent(in) :: act(:)
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream

      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :), lsf4_rrr_rA(:, :, :, :, :)
      integer :: nac, j, k, l, s, ia, atomA, slotA

      ! Every nuclear output is active-indexed and caller-sized; the harness
      ! resolves the index space from the extent it gets back (see `nuc_slot`)
      nac = lsf%active_count()
      allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
      allocate (lsf3_rr_rA(ndim, ndim, ndim, nac))

      call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      call check_user_space(lsf1_rA, act, "f3rrA_lsf1_rA", stream)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf1_rA, 2), atomA, act, "f3rrA_lsf1_rA", stream)
         do s = 1, ndim
            call stream%emit(id, "f3rrA_lsf1_rA", idx6(s, atomA), &
                             slot_read(lsf1_rA(s, :), slotA))
         end do
      end do
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf2_r_rA, 3), atomA, act, "f3rrA_lsf2_r_rA", stream)
         do s = 1, ndim
            do j = 1, ndim
               call stream%emit(id, "f3rrA_lsf2_r_rA", idx6(j, s, atomA), &
                                slot_read(lsf2_r_rA(j, s, :), slotA))
            end do
         end do
      end do
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf3_rr_rA, 4), atomA, act, "f3_rr_rA", stream)
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  call stream%emit(id, "f3_rr_rA", idx6(j, k, s, atomA), &
                                   slot_read(lsf3_rr_rA(j, k, s, :), slotA))
               end do
            end do
         end do
      end do

      if (id%kind /= kind_svdw) return

      allocate (lsf4_rrr_rA(ndim, ndim, ndim, ndim, nac))
      call lsf%f4_rrr_rA(lsf4_rrr_rA)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf4_rrr_rA, 5), atomA, act, "f4_rrr_rA", stream)
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  do l = k, ndim
                     call stream%emit(id, "f4_rrr_rA", idx6(j, k, l, s, atomA), &
                                      slot_read(lsf4_rrr_rA(j, k, l, s, :), slotA))
                  end do
               end do
            end do
         end do
      end do
   end subroutine emit_nuclear_block

   !> Emit the derivatives that carry two nuclear indices
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[in]    sel     Selected atom indices
   !> @param[in]    act     Atom -> active index map
   !> @param[inout] stream  Record sink
   subroutine emit_pair_block(lsf, id, sel, act, stream)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Atom -> active index map
      integer, intent(in) :: act(:)
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream

      real(wp), allocatable :: lsf2_rArB(:, :, :, :), lsf3_r_rArB(:, :, :, :, :)
      real(wp), allocatable :: lsf4_rr_rArB(:, :, :, :, :, :)
      integer :: nac, npair, j, k, s, t, ia, ib, atomA, atomB, slotA, slotB

      nac = lsf%active_count()
      npair = min(n_pair_atoms, size(sel))
      allocate (lsf2_rArB(ndim, nac, ndim, nac))
      allocate (lsf3_r_rArB(ndim, ndim, nac, ndim, nac))
      allocate (lsf4_rr_rArB(ndim, ndim, ndim, nac, ndim, nac))

      call lsf%f2_rArB(lsf2_rArB)
      do ia = 1, npair
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf2_rArB, 2), atomA, act, "f2_rArB", stream)
         do ib = 1, npair
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf2_rArB, 4), atomB, act, "f2_rArB", stream)
            do s = 1, ndim
               do t = 1, ndim
                  call stream%emit(id, "f2_rArB", idx6(s, atomA, t, atomB), &
                                   slot_read2(lsf2_rArB(s, :, t, :), slotA, slotB))
               end do
            end do
         end do
      end do

      call lsf%f3_r_rArB(lsf3_r_rArB)
      do ia = 1, npair
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf3_r_rArB, 3), atomA, act, "f3_r_rArB", stream)
         do ib = 1, npair
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf3_r_rArB, 5), atomB, act, "f3_r_rArB", stream)
            do s = 1, ndim
               do t = 1, ndim
                  do j = 1, ndim
                     call stream%emit(id, "f3_r_rArB", idx6(j, s, atomA, t, atomB), &
                                      slot_read2(lsf3_r_rArB(j, s, :, t, :), slotA, slotB))
                  end do
               end do
            end do
         end do
      end do

      call lsf%f4_rr_rArB(lsf4_rr_rArB)
      do ia = 1, npair
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf4_rr_rArB, 4), atomA, act, "f4_rr_rArB", stream)
         do ib = 1, npair
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf4_rr_rArB, 6), atomB, act, "f4_rr_rArB", stream)
            do s = 1, ndim
               do t = 1, ndim
                  do j = 1, ndim
                     do k = j, ndim
                        call stream%emit(id, "f4_rr_rArB", &
                                         idx6(j, k, s, atomA, t, atomB), &
                                         slot_read2(lsf4_rr_rArB(j, k, s, :, t, :), &
                                                    slotA, slotB))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine emit_pair_block

   !> Emit the surface-normalised value and its nuclear gradient
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[in]    sel     Selected atom indices
   !> @param[in]    act     Atom -> active index map
   !> @param[inout] stream  Record sink
   subroutine emit_normalized_block(lsf, id, sel, act, stream)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Atom -> active index map
      integer, intent(in) :: act(:)
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream

      real(wp), allocatable :: norm1_rA(:, :)
      real(wp) :: norm0
      integer :: s, ia, atomA, slotA

      allocate (norm1_rA(ndim, lsf%active_count()))

      call lsf%normalized_f01_rA(norm0, norm1_rA)
      call check_user_space(norm1_rA, act, "normalized_f1_rA", stream)
      call stream%emit(id, "normalized_f0", idx6(), norm0)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(norm1_rA, 2), atomA, act, "normalized_f1_rA", stream)
         do s = 1, ndim
            call stream%emit(id, "normalized_f1_rA", idx6(s, atomA), &
                             slot_read(norm1_rA(s, :), slotA))
         end do
      end do
   end subroutine emit_normalized_block

   !> Pack up to six indices into the fixture's index tuple, zero-filling the
   !> slots a quantity does not use
   !>
   !> @param[in] i1  First index
   !> @param[in] i2  Second index
   !> @param[in] i3  Third index
   !> @param[in] i4  Fourth index
   !> @param[in] i5  Fifth index
   !> @param[in] i6  Sixth index
   !> @returns       Index tuple
   pure function idx6(i1, i2, i3, i4, i5, i6) result(tuple)
      !> Indices to pack, in order
      integer, intent(in), optional :: i1, i2, i3, i4, i5, i6
      integer :: tuple(6)

      tuple = 0
      if (present(i1)) tuple(1) = i1
      if (present(i2)) tuple(2) = i2
      if (present(i3)) tuple(3) = i3
      if (present(i4)) tuple(4) = i4
      if (present(i5)) tuple(5) = i5
      if (present(i6)) tuple(6) = i6
   end function idx6

   !* ================================================================================= *!
   !*                            Index-space-safe reads                                 *!
   !* ================================================================================= *!
   !>
   !> Every nuclear array the LSF hands back is read through this family, never
   !> by indexing it directly. The rule the harness follows is: records are keyed
   !> by *user-space* atom id (so the fixture survives an index-space change in
   !> `src/`), while the slot to read is derived from the extent of the array the
   !> routine actually returned. Nothing here assumes which space a routine uses.

   !> Resolve the slot holding user-space atom `atom` in a nuclear dimension of
   !> length `extent`.
   !>
   !> `extent == ncenters` means the routine is user-indexed and the slot is the
   !> atom id itself. `extent == active_count()` means it is active-indexed and
   !> the slot is the atom's active index, or 0 when screening dropped it - an
   !> all-zero block, which is that atom's true derivative contribution. Any
   !> other extent is a structural change this harness does not understand and is
   !> reported rather than read; the returned slot is then 0, so a wrong guess can
   !> never turn into an out-of-bounds access.
   !>
   !> When every atom is active the two spaces coincide, because [[active_map]]
   !> asserts the active list is then in ascending user-space order.
   !>
   !> @param[in]    extent  Length of the nuclear dimension as returned
   !> @param[in]    atom    User-space atom id
   !> @param[in]    act     Atom -> active index map (0 = dropped), size ncenters
   !> @param[in]    label   Quantity name, for the failure message
   !> @param[inout] stream  Record sink, to report a structural change
   !> @returns              Slot to read, or 0 if the atom has none
   integer function nuc_slot(extent, atom, act, label, stream)
      !> Length of the nuclear dimension as returned
      integer, intent(in) :: extent
      !> User-space atom id
      integer, intent(in) :: atom
      !> Atom -> active index map (0 = dropped)
      integer, intent(in) :: act(:)
      !> Quantity name, for the failure message
      character(len=*), intent(in) :: label
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream

      if (extent == size(act)) then
         nuc_slot = atom
      else if (extent == count(act > 0)) then
         nuc_slot = act(atom)
      else
         call stream%flag_problem(label//": extent "//itoa(extent)// &
                                  " is neither ncenters "//itoa(size(act))// &
                                  " nor active_count "//itoa(count(act > 0))// &
                                  " - index space changed in src/")
         nuc_slot = 0
         return
      end if
      if (nuc_slot < 1 .or. nuc_slot > extent) nuc_slot = 0
   end function nuc_slot

   !> Read the nuclear slot of a rank-reduced slice, or 0 when the atom has none
   !>
   !> Call sites slice every fixed index away first - `t(j, k, s, :)` - so one
   !> routine serves nuclear arrays of any rank.
   !>
   !> @param[in] v     Nuclear slice, one element per slot
   !> @param[in] slot  Slot from [[nuc_slot]]
   !> @returns         Element or 0
   pure real(wp) function slot_read(v, slot) result(val)
      !> Nuclear slice
      real(wp), intent(in) :: v(:)
      !> Slot to read
      integer, intent(in) :: slot

      val = 0.0_wp
      if (slot > 0) val = v(slot)
   end function slot_read

   !> Two-nucleus counterpart of [[slot_read]]; 0 when either atom has no slot
   !>
   !> @param[in] m      Nuclear slice `(A, B)`
   !> @param[in] slotA  Slot of A from [[nuc_slot]]
   !> @param[in] slotB  Slot of B from [[nuc_slot]]
   !> @returns          Element or 0
   pure real(wp) function slot_read2(m, slotA, slotB) result(val)
      !> Nuclear slice
      real(wp), intent(in) :: m(:, :)
      !> Slot of A
      integer, intent(in) :: slotA
      !> Slot of B
      integer, intent(in) :: slotB

      val = 0.0_wp
      if (slotA > 0 .and. slotB > 0) val = m(slotA, slotB)
   end function slot_read2

   !> Guard a *caller-sized* nuclear output, whose extent cannot reveal its index
   !> space because the harness chose it.
   !>
   !> `f3_rr_rA`'s optional `lsf1_rA` (and `normalized_f01_rA`'s gradient) are
   !> passed in as `(3, ncenters)` buffers and scattered into by user-space atom
   !> id. If a future refactor made them active-indexed instead, the extent check
   !> in [[nuc_slot]] could not see it - but the columns of screened-away atoms
   !> would stop being zero. That is what this asserts.
   !>
   !> @param[in]    t       Caller-sized `(axis, A)` output
   !> @param[in]    act     Atom -> active index map (0 = dropped)
   !> @param[in]    label   Quantity name, for the failure message
   !> @param[inout] stream  Record sink, to report a structural change
   subroutine check_user_space(t, act, label, stream)
      !> Caller-sized nuclear output
      real(wp), intent(in) :: t(:, :)
      !> Atom -> active index map
      integer, intent(in) :: act(:)
      !> Quantity name
      character(len=*), intent(in) :: label
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream

      integer :: atom

      if (size(t, 2) /= size(act)) return
      do atom = 1, size(act)
         if (act(atom) > 0) cycle
         if (maxval(abs(t(:, atom))) == 0.0_wp) cycle
         call stream%flag_problem(label//": atom "//itoa(atom)//" is screened away yet "// &
                                  "its user-space column is non-zero - the output is no "// &
                                  "longer user-indexed")
         return
      end do
   end subroutine check_user_space

   !* ================================================================================= *!
   !*                          Screening bookkeeping helpers                            *!
   !* ================================================================================= *!

   !> Guard the unscreened reference: `screening_threshold = 0` must leave
   !> every centre active. If a future change to the SSD screen invalidates
   !> that, the fixture would silently stop being an unscreened reference.
   !>
   !> @param[in]  nact   Active count reported by the LSF
   !> @param[in]  nat    Number of centres
   !> @param[in]  tag    Case tag, for the message
   !> @param[in]  ip     Evaluation-point index, for the message
   !> @param[out] error  testdrive failure
   subroutine assert_unscreened(nact, nat, tag, ip, error)
      !> Active count
      integer, intent(in) :: nact
      !> Number of centres
      integer, intent(in) :: nat
      !> Case tag
      character(len=*), intent(in) :: tag
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      if (nact == nat) return
      call test_failed(error, "threshold 0 no longer disables screening for "//trim(tag)// &
                       " point "//itoa(ip)//": active "//itoa(nact)//" of "//itoa(nat))
   end subroutine assert_unscreened

   !> Build the atom -> active-index map of a prepared LSF (0 = dropped)
   !>
   !> Also asserts the property [[nuc_slot]] depends on to disambiguate the two
   !> index spaces when they have equal length: with every atom active, the
   !> active list must be in ascending user-space order, so active index equals
   !> atom id. `lsf_base_rebuild_screening` guarantees this for a full `prepare`
   !> by walking the full scan through `orig_to_sorted`; if that ever changes,
   !> this fires instead of the harness silently reading the wrong slot.
   !>
   !> @param[in]    lsf     Prepared LSF (either concrete)
   !> @param[out]   act     `act(atom)` = active index or 0
   !> @param[inout] stream  Record sink, to report a broken invariant
   subroutine active_map(lsf, act, stream)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Atom -> active index
      integer, allocatable, intent(out) :: act(:)
      !> Record sink
      type(golden_stream_type), intent(inout) :: stream

      integer :: i, nat, nact

      nat = lsf%ncenters
      nact = lsf%active_count()
      allocate (act(nat), source=0)
      do i = 1, nact
         act(lsf%active_atom(i)) = i
      end do

      if (nact /= nat) return
      do i = 1, nat
         if (act(i) == i) cycle
         call stream%flag_problem("fully active list is no longer ascending: atom "// &
                                  itoa(i)//" sits at active slot "//itoa(act(i)))
         return
      end do
   end subroutine active_map

   !* ================================================================================= *!
   !*                                 Stream plumbing                                   *!
   !* ================================================================================= *!

   !> Append one record. The traversal calls this and nothing else, so the
   !> fixture's order is the traversal's order by construction.
   !>
   !> @param[inout] self      Record sink
   !> @param[in]    id        Concrete, case, point and screening flag
   !> @param[in]    quantity  Quantity name
   !> @param[in]    idx       Index tuple, unused slots 0
   !> @param[in]    val       Value
   subroutine stream_emit(self, id, quantity, idx, val)
      !> Record sink
      class(golden_stream_type), intent(inout) :: self
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Quantity name
      character(len=*), intent(in) :: quantity
      !> Index tuple
      integer, intent(in) :: idx(6)
      !> Value
      real(wp), intent(in) :: val

      type(golden_record_type), allocatable :: bigger(:)

      if (.not. allocated(self%rec)) allocate (self%rec(4096))
      if (self%n == size(self%rec)) then
         allocate (bigger(2*size(self%rec)))
         bigger(1:size(self%rec)) = self%rec
         call move_alloc(bigger, self%rec)
      end if

      self%n = self%n + 1
      self%rec(self%n) = golden_record_type(id, quantity, idx, val)
   end subroutine stream_emit

   !> Record a *structural* problem - a shape, an index space or an invariant
   !> that no longer matches what this harness understands.
   !>
   !> Collected rather than raised so that the traversal always runs to
   !> completion; every caller of [[traverse]] reports a non-zero count as a
   !> test failure before it looks at any number.
   !>
   !> @param[inout] self  Record sink
   !> @param[in]    text  Message describing the problem
   subroutine stream_flag_problem(self, text)
      !> Record sink
      class(golden_stream_type), intent(inout) :: self
      !> Message describing the problem
      character(len=*), intent(in) :: text

      self%nproblem = self%nproblem + 1
      if (self%nproblem == 1) self%problem = text
   end subroutine stream_flag_problem

   !> Turn a collected structural problem into a testdrive failure
   !>
   !> @param[in]  stream  Traversed stream
   !> @param[out] error   testdrive failure, allocated only if a problem was seen
   subroutine check_problems(stream, error)
      !> Traversed stream
      type(golden_stream_type), intent(in) :: stream
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      if (stream%nproblem == 0) return
      call test_failed(error, itoa(stream%nproblem)//" structural problem(s); first: "// &
                       stream%problem)
   end subroutine check_problems

   !> Hold an emitted stream against the parsed fixture, record by record
   !>
   !> @param[in]  got    Emitted stream
   !> @param[in]  ref    Parsed fixture
   !> @param[in]  path   Fixture path, for the failure message
   !> @param[out] error  testdrive failure
   subroutine compare_stream(got, ref, path, error)
      !> Emitted stream
      type(golden_stream_type), intent(in) :: got
      !> Parsed fixture
      type(golden_record_type), intent(in) :: ref(:)
      !> Fixture path
      character(len=*), intent(in) :: path
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      integer :: i, nfail
      character(len=:), allocatable :: first
      character(len=64) :: values

      if (got%n /= size(ref)) then
         call test_failed(error, "golden fixture length mismatch: evaluated "// &
                          itoa(got%n)//" records, fixture has "//itoa(size(ref))//" - "// &
                          path//" no longer describes this traversal")
         return
      end if

      nfail = 0
      first = ""
      do i = 1, got%n
         if (.not. same_label(got%rec(i), ref(i))) then
            nfail = nfail + 1
            if (nfail == 1) first = "record "//itoa(i)//" is labelled '"// &
                                    record_label(ref(i))//"' in the fixture but '"// &
                                    record_label(got%rec(i))//"' now"
            cycle
         end if
         if (rel_deviation(got%rec(i)%val, ref(i)%val) <= golden_tol) cycle
         nfail = nfail + 1
         if (nfail == 1) then
            write (values, '(a,es24.16,a,es24.16)') " golden ", ref(i)%val, " now ", &
               got%rec(i)%val
            first = "record "//itoa(i)//" "//record_label(ref(i))//values
         end if
      end do

      if (nfail == 0) return
      call test_failed(error, itoa(nfail)//" golden record(s) deviate; first: "//first)
   end subroutine compare_stream

   !> `.true.` when two records name the same number
   !>
   !> @param[in] a  First record
   !> @param[in] b  Second record
   !> @returns      Whether the labels agree
   pure logical function same_label(a, b)
      !> First record
      type(golden_record_type), intent(in) :: a
      !> Second record
      type(golden_record_type), intent(in) :: b

      same_label = a%id%kind == b%id%kind .and. a%id%case_tag == b%id%case_tag &
                   .and. a%id%ip == b%id%ip .and. a%id%flag == b%id%flag &
                   .and. a%quantity == b%quantity .and. all(a%idx == b%idx)
   end function same_label

   !> Left-justified decimal form of `n`, for building failure messages
   !>
   !> @param[in] n  Number to render
   !> @returns      Decimal digits, no padding
   function itoa(n) result(text)
      !> Number to render
      integer, intent(in) :: n
      character(len=:), allocatable :: text
      character(len=32) :: buf

      write (buf, '(i0)') n
      text = trim(buf)
   end function itoa

   !> Human-readable identity of one record, e.g. `svdw ch4 point 3 S f3_rrr [ 1 1 2 0 0 0 ]`
   !>
   !> @param[in] rec  Record to describe
   !> @returns        One-line description
   function record_label(rec) result(text)
      !> Record to describe
      type(golden_record_type), intent(in) :: rec
      character(len=:), allocatable :: text
      integer :: i

      text = trim(rec%id%kind)//" "//trim(rec%id%case_tag)//" point "//itoa(rec%id%ip)// &
             " "//rec%id%flag//" "//trim(rec%quantity)//" ["
      do i = 1, size(rec%idx)
         text = text//" "//itoa(rec%idx(i))
      end do
      text = text//" ]"
   end function record_label

   !> Parse a fixture file into an array of records. `#` comments and blank
   !> lines are skipped; everything else must parse as a record.
   !>
   !> @param[in]  path   Fixture path
   !> @param[out] ref    Parsed records
   !> @param[out] error  testdrive failure
   subroutine load_fixture(path, ref, error)
      !> Fixture path
      character(len=*), intent(in) :: path
      !> Parsed records
      type(golden_record_type), allocatable, intent(out) :: ref(:)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, n
      character(len=256) :: line

      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "missing golden fixture '"//path// &
                          "' - it is committed data, restore it from git")
         return
      end if

      n = 0
      do
         read (unit, '(a)', iostat=stat) line
         if (stat /= 0) exit
         if (is_record_line(line)) n = n + 1
      end do

      allocate (ref(n))
      rewind (unit)

      n = 0
      do
         read (unit, '(a)', iostat=stat) line
         if (stat /= 0) exit
         if (.not. is_record_line(line)) cycle
         n = n + 1
         read (line, *, iostat=stat) ref(n)%id%kind, ref(n)%id%case_tag, ref(n)%id%ip, &
            ref(n)%id%flag, ref(n)%quantity, ref(n)%idx, ref(n)%val
         if (stat /= 0) then
            close (unit)
            call test_failed(error, "malformed record "//itoa(n)//" in "//path)
            return
         end if
      end do
      close (unit)
   end subroutine load_fixture

   !> `.true.` for a line carrying a record (not blank, not a `#` comment)
   !>
   !> @param[in] line  Raw line
   !> @returns         Whether the line should be parsed
   pure logical function is_record_line(line)
      !> Raw line
      character(len=*), intent(in) :: line
      character(len=1) :: lead

      is_record_line = .false.
      if (len_trim(line) == 0) return
      lead = adjustl(line)
      if (lead == "#") return
      is_record_line = .true.
   end function is_record_line

   !* ================================================================================= *!
   !*                            Stream reproducibility                                 *!
   !* ================================================================================= *!

   !> Evaluate everything twice in one process and require the two value streams
   !> to be bit-identical.
   !>
   !> This is the standing guard against the class of bug that once hid here: the
   !> harness read a pair tensor past `n_active`, so its records came from
   !> uninitialised heap and changed from run to run - the fixture regenerated
   !> fine and then failed on 5 of the next 6 runs. Reading past a returned
   !> extent is now structurally impossible (see [[nuc_slot]]), and this check
   !> keeps *any* run-to-run instability - stale memory, a race, leftover
   !> accumulator state - a test failure rather than a flaky fixture. It also
   !> catches a NaN, which can never compare equal to itself.
   !>
   !> The second pass runs against a heap the first pass has already churned, so
   !> a read of uninitialised memory has every chance to return something else.
   subroutine test_stream_reproducible(error)
      type(error_type), allocatable, intent(out) :: error

      call compare_two_passes(error, kind_svdw)
      if (allocated(error)) return
      call compare_two_passes(error, kind_cfc)
   end subroutine test_stream_reproducible

   !> Run one concrete's traversal twice and diff the two streams
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine compare_two_passes(error, kind)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      type(golden_stream_type) :: first, second
      integer :: i
      character(len=128) :: tail

      call traverse(kind, first, error)
      if (allocated(error)) return
      call traverse(kind, second, error)
      if (allocated(error)) return

      call check_problems(first, error)
      if (allocated(error)) return
      call check_problems(second, error)
      if (allocated(error)) return

      if (first%n /= second%n) then
         call test_failed(error, kind//" record count is not reproducible: "// &
                          itoa(first%n)//" then "//itoa(second%n))
         return
      end if

      do i = 1, first%n
         if (first%rec(i)%val == second%rec(i)%val) cycle
         write (tail, '(a,es24.16,a,es24.16)') &
            " differs between two passes in one process: ", first%rec(i)%val, " then ", &
            second%rec(i)%val
         call test_failed(error, kind//" record "//itoa(i)//trim(tail))
         return
      end do
   end subroutine compare_two_passes

   !* ================================================================================= *!
   !*                              Tensor symmetry checks                               *!
   !* ================================================================================= *!

   !> The golden fixture keeps only the unique components of the slots that are
   !> symmetric by construction. These checks are the other half of that deal:
   !> they assert the full tensors really are symmetric in exactly those slots,
   !> so nothing is left unpinned. Tolerance is loose on purpose - the target is
   !> a wrong permutation, not accumulated roundoff.
   subroutine test_svdw_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      call run_symmetry(error, kind_svdw)
   end subroutine test_svdw_symmetry

   !> CFC counterpart of [[test_svdw_symmetry]] (orders 2 and 3 only)
   subroutine test_cfc_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      call run_symmetry(error, kind_cfc)
   end subroutine test_cfc_symmetry

   !> Walk every case and point of one concrete, asserting every symmetric slot
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine run_symmetry(error, kind)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      type(structure_type) :: mol
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      integer :: icase, ip

      do icase = 1, n_cases(kind)
         call load_case(golden_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return

         do ip = 1, n_points
            call new_lsf(lsf, kind, golden_cases(icase), mol, radii, production_threshold)
            call prepare_lsf(lsf, points(:, ip), error)
            if (allocated(error)) return

            call check_symmetry_low(lsf, error)
            if (allocated(error)) return
            if (kind /= kind_svdw) cycle
            call check_symmetry_high(lsf, error)
            if (allocated(error)) return
         end do

         deallocate (radii, points)
      end do
   end subroutine run_symmetry

   !> Symmetric slots of the derivative orders both concretes provide
   !>
   !> @param[in]  lsf    Prepared LSF
   !> @param[out] error  testdrive failure
   subroutine check_symmetry_low(lsf, error)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim), lsf3_rrr(ndim, ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :), lsf3_rr_rA(:, :, :, :)
      real(wp) :: dev_jk, dev_kl
      integer :: nac, j, s, iA

      nac = lsf%active_count()
      allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
      allocate (lsf3_rr_rA(ndim, ndim, ndim, nac))

      call lsf%f012_r(lsf0, lsf1_r, lsf2_rr)
      call check_sym(error, asym(lsf2_rr), maxval(abs(lsf2_rr)), "f012_lsf2_rr")
      if (allocated(error)) return

      call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
      dev_jk = 0.0_wp
      dev_kl = 0.0_wp
      do j = 1, ndim
         dev_jk = max(dev_jk, asym(lsf3_rrr(:, :, j)))
         dev_kl = max(dev_kl, asym(lsf3_rrr(j, :, :)))
      end do
      call check_sym(error, dev_jk, maxval(abs(lsf3_rrr)), "f3_rrr(jk)")
      if (allocated(error)) return
      call check_sym(error, dev_kl, maxval(abs(lsf3_rrr)), "f3_rrr(kl)")
      if (allocated(error)) return

      call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      dev_jk = 0.0_wp
      do iA = 1, nac
         do s = 1, ndim
            dev_jk = max(dev_jk, asym(lsf3_rr_rA(:, :, s, iA)))
         end do
      end do
      call check_sym(error, dev_jk, maxval(abs(lsf3_rr_rA)), "f3_rr_rA(jk)")
   end subroutine check_symmetry_low

   !> Symmetric slots of the order-4 tensors, which only SvdW provides
   !>
   !> @param[in]  lsf    Prepared LSF
   !> @param[out] error  testdrive failure
   subroutine check_symmetry_high(lsf, error)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: lsf4_rrrr(ndim, ndim, ndim, ndim)
      real(wp), allocatable :: lsf4_rrr_rA(:, :, :, :, :), lsf4_rr_rArB(:, :, :, :, :, :)
      real(wp) :: dev_jk, dev_kl, dev_lm
      integer :: nac, j, k, l, s, iA

      nac = lsf%active_count()
      allocate (lsf4_rrr_rA(ndim, ndim, ndim, ndim, nac))
      allocate (lsf4_rr_rArB(ndim, ndim, ndim, nac, ndim, nac))

      call lsf%f4_rrrr(lsf4_rrrr)
      dev_jk = 0.0_wp
      dev_kl = 0.0_wp
      dev_lm = 0.0_wp
      do j = 1, ndim
         do k = 1, ndim
            dev_jk = max(dev_jk, asym(lsf4_rrrr(:, :, j, k)))
            dev_kl = max(dev_kl, asym(lsf4_rrrr(j, :, :, k)))
            dev_lm = max(dev_lm, asym(lsf4_rrrr(j, k, :, :)))
         end do
      end do
      call check_sym(error, dev_jk, maxval(abs(lsf4_rrrr)), "f4_rrrr(jk)")
      if (allocated(error)) return
      call check_sym(error, dev_kl, maxval(abs(lsf4_rrrr)), "f4_rrrr(kl)")
      if (allocated(error)) return
      call check_sym(error, dev_lm, maxval(abs(lsf4_rrrr)), "f4_rrrr(lm)")
      if (allocated(error)) return

      call lsf%f4_rrr_rA(lsf4_rrr_rA)
      dev_jk = 0.0_wp
      dev_kl = 0.0_wp
      do iA = 1, nac
         do s = 1, ndim
            do l = 1, ndim
               dev_jk = max(dev_jk, asym(lsf4_rrr_rA(:, :, l, s, iA)))
               dev_kl = max(dev_kl, asym(lsf4_rrr_rA(l, :, :, s, iA)))
            end do
         end do
      end do
      call check_sym(error, dev_jk, maxval(abs(lsf4_rrr_rA)), "f4_rrr_rA(jk)")
      if (allocated(error)) return
      call check_sym(error, dev_kl, maxval(abs(lsf4_rrr_rA)), "f4_rrr_rA(kl)")
      if (allocated(error)) return

      call lsf%f4_rr_rArB(lsf4_rr_rArB)
      dev_jk = 0.0_wp
      do j = 1, ndim
         do k = 1, ndim
            dev_jk = max(dev_jk, maxval(abs(lsf4_rr_rArB(j, k, :, :, :, :) &
                                            - lsf4_rr_rArB(k, j, :, :, :, :))))
         end do
      end do
      call check_sym(error, dev_jk, maxval(abs(lsf4_rr_rArB)), "f4_rr_rArB(jk)")
   end subroutine check_symmetry_high

   !> Largest asymmetry of a square matrix, `max |m - m^T|`
   !>
   !> @param[in] m  Matrix, or a rank-reduced slice of a higher-rank tensor
   !> @returns      Largest absolute deviation from symmetry
   pure real(wp) function asym(m) result(dev)
      !> Matrix to test
      real(wp), intent(in) :: m(:, :)

      dev = maxval(abs(m - transpose(m)))
   end function asym

   !> Fail unless a tensor's asymmetry stays within `symmetry_rel_tol` of its
   !> largest element
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  dev    Largest absolute deviation from symmetry
   !> @param[in]  scale  Largest absolute element of the tensor
   !> @param[in]  label  Name reported on failure
   subroutine check_sym(error, dev, scale, label)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Largest absolute deviation from symmetry
      real(wp), intent(in) :: dev
      !> Largest absolute element of the tensor
      real(wp), intent(in) :: scale
      !> Name reported on failure
      character(len=*), intent(in) :: label

      if (dev <= max(symmetry_abs_tol, symmetry_rel_tol*scale)) return
      call test_failed(error, "tensor slot not symmetric: "//label)
   end subroutine check_sym

end module test_cavity_drop_lsf_golden
