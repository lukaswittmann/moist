!> Golden numerical fixture for the DROP level set functions (SvdW and CFC)
!>
!> The rest of the LSF unit suite ([[test_cavity_drop_lsf]],
!> [[test_cavity_drop_cfc]]) is finite-difference / self-consistency based:
!> it pins *derivative consistency*, not the values themselves. A rewrite of
!> the blending math that got a weight or a body-order term subtly wrong could
!> stay perfectly FD-consistent while describing a different surface. This
!> module closes that hole by comparing the current implementation against a
!> checked-in table of reference numbers.
!>
!> Fixtures live in
!>   * `test/unit/data/lsf_golden_svdw.txt`
!>   * `test/unit/data/lsf_golden_cfc.txt`
!>
!> and are produced by *this same module*: the traversal that evaluates the
!> LSFs feeds a sink which either writes a record or compares it. Writer and
!> checker therefore cannot drift apart in ordering or labelling.
!>
!> Regenerating the fixtures (deliberate act; do it only when the math is
!> *meant* to change):
!>
!> ```
!> MOIST_LSF_GOLDEN_REGENERATE=1 \
!> MOIST_GIT_COMMIT=$(git rev-parse HEAD) \
!> MOIST_SOURCE_ROOT=$PWD \
!>   ./<builddir>/test/unit/tester cavity_drop_lsf_golden
!> ```
!>
!> In regeneration mode the two suite entries rewrite their file and pass.
!>
!> Every quantity is recorded **twice** per evaluation point:
!>   * `S` - screened, exactly as production runs it
!>     (`screening_threshold = 1.0e-11`, the DROP default, via `prepare`)
!>   * `U` - unscreened reference (`screening_threshold = 0.0`, `prepare`)
!>
!> Setting the inherited `screening_threshold` to zero *before* `update` is what
!> disables screening. `update` runs `lsf_base_rebuild_screening`, which asks the
!> concrete for `screening_offset(radius)`; both concretes return
!> `huge(0.0_wp)` as soon as the threshold is non-positive, so the per-atom
!> reach `radius + offset` saturates and the cached gate
!> `cand_screen(lsf_cand_bound_sq, :)` is set to `huge(0.0_wp)`. The per-point
!> gate is a single squared-distance compare against that bound, which no atom
!> can fail. The dumper asserts `active_count() == ncenters` at every unscreened
!> point, so a future change to that mechanism cannot silently degrade the
!> reference.
!>
!> A third record group per point carries flag `D`: `f0_delta`, the screened
!> minus unscreened LSF value. It lets a reader see directly whether the
!> screening error stays inside the band the threshold justifies.
!>
!> Record layout (one record per line, blanks separate the fields):
!>
!> ```
!> kind  case  ip  flag  quantity  i1 i2 i3 i4 i5 i6  value
!> ```
!>
!>   * `kind`     `svdw` or `cfc`
!>   * `case`     case tag; see the `# case` legend in the file header
!>   * `ip`       evaluation-point index, 1-based (see `# point` legend)
!>   * `flag`     `S` screened, `U` unscreened, `D` screened-minus-unscreened
!>   * `quantity` accessor / output name
!>   * `i1..i6`   index tuple, unused slots are 0; `j,k,l,m` are spatial axes,
!>                `s,t` nuclear axes, `A,B` **user-space atom indices**
!>   * `value`    `es24.16` (17 significant digits: exact IEEE-754 round trip)
!>
!> Index tuples per quantity:
!>
!> | quantity            | tuple             |
!> |---------------------|-------------------|
!> | `n_active`          | -                 |
!> | `f0`                | -                 |
!> | `f012_lsf0`         | -                 |
!> | `f012_lsf1_r`       | `j`               |
!> | `f012_lsf2_rr`      | `j k`             |
!> | `f3_rrr`            | `j k l`           |
!> | `f4_rrrr`           | `j k l m`         |
!> | `f3rrA_lsf1_rA`     | `s A`             |
!> | `f3rrA_lsf2_r_rA`   | `j s A`           |
!> | `f3_rr_rA`          | `j k s A`         |
!> | `f4_rrr_rA`         | `j k l s A`       |
!> | `f2_rArB`           | `s A t B`         |
!> | `f3_r_rArB`         | `j s A t B`       |
!> | `f4_rr_rArB`        | `j k s A t B`     |
!> | `normalized_f0`     | -                 |
!> | `normalized_f1_rA`  | `s A`             |
!> | `f0_delta`          | -                 |
!>
!> Slots that are symmetric by construction are dumped once only
!> (`j <= k <= l <= m` for the pure spatial orders, `j <= k` for the leading
!> pair of `f3_rr_rA` / `f4_rr_rArB`, `j <= k <= l` for `f4_rrr_rA`). The
!> discarded components are not left unchecked: `*_tensor_symmetry` verifies
!> the full tensors really are symmetric in those slots, so the reduced dump
!> plus the symmetry check together pin the whole tensor.
!>
!> Records are always keyed by **user-space** atom id, whatever index space the
!> routine happens to return. Some outputs (`f2_rArB`, `f3_r_rArB`,
!> `f4_rr_rArB`) are active-indexed; the harness never assumes which, it derives
!> the space from the extent of the array it was handed and resolves the slot
!> through [[nuc_slot]]. A selected atom that screening dropped has no slot and
!> contributes an all-zero block, which is the correct value for it, not a
!> placeholder. Keying on user-space ids is what lets the fixture survive an
!> index-space change in `src/` unchanged - `f2_rArB` moving from
!> `(3, ncenters, 3, ncenters)` to `(3, n_active, 3, n_active)` shifts values,
!> not record labels.
!>
!> Scope knobs (all `parameter`s below, change them and regenerate):
!> `n_points`, `n_sel_atoms`, `n_pair_atoms` and the case tables set the file
!> size; the shipped settings give 34620 SvdW and 5325 CFC records
!> (3.0 MB and 465 kB).
!>
!> Index spaces move: `f2_rArB` was user-indexed when this harness was
!> first written and is active-indexed now. Nothing here hard-codes that. Every
!> nuclear read goes through [[nuc_slot]], which derives the space from the
!> extent of the array the routine actually returned and reports anything it
!> does not recognise, so a further change shifts values, never labels, and can
!> never become an out-of-bounds read. [[test_stream_reproducible]] is the
!> standing check that it has not become one anyway. When auditing this by hand,
!> a bounds-checked build is the fastest confirmation:
!>
!> ```
!> meson setup build-bounds -Dfortran_args="-fcheck=bounds,pointer,mem \
!>       -finit-real=snan -finit-integer=-99999999"
!> ```
!>
!> Tolerance: relative `1.0e-12` with an absolute floor of `1.0e-12`, i.e. a
!> record passes when `|now - ref| <= max(1e-12, 1e-12 * |ref|)`. The floor
!> keeps entries that are exactly (or nearly) zero from demanding bit
!> equality of a cancelling sum. Regenerated on the same toolchain the
!> agreement is bit-for-bit; the tolerance is headroom for other platforms.
module test_cavity_drop_lsf_golden
   use, intrinsic :: iso_fortran_env, only: error_unit
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use moist_utils_env, only: get_env
   use test_helpers, only: get_test_radii
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private

   public :: collect_cavity_drop_lsf_golden

   integer, parameter :: ndim = 3

   !> Screening threshold production uses (DROP `parameters.f90` default)
   real(wp), parameter :: production_threshold = 1.0e-11_wp

   !> Relative tolerance of the golden comparison
   real(wp), parameter :: golden_rel_tol = 1.0e-12_wp
   !> Absolute floor of the golden comparison
   real(wp), parameter :: golden_abs_tol = 1.0e-12_wp
   !> Tolerance of the companion symmetry checks (gross-asymmetry guard, not
   !> a roundoff assertion), relative to the largest element of the tensor
   real(wp), parameter :: symmetry_rel_tol = 1.0e-8_wp

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
   !*                                  Record sink                                      *!
   !* ================================================================================= *!

   !> One parsed fixture line
   type :: golden_record_type
      !> `svdw` or `cfc`
      character(len=8) :: kind
      !> Case tag
      character(len=12) :: case_tag
      !> Evaluation-point index
      integer :: ip
      !> `S`, `U` or `D`
      character(len=1) :: flag
      !> Quantity name
      character(len=20) :: quantity
      !> Index tuple, unused slots 0
      integer :: idx(6)
      !> Reference value
      real(wp) :: val
   end type golden_record_type

   !> Write format of a record line; the reader is list-directed
   character(len=*), parameter :: record_fmt = &
                                  '(a4,1x,a10,1x,i2,1x,a1,1x,a18,6(1x,i3),1x,es24.16)'

   !> `.true.` while regenerating the fixture, `.false.` while checking it
   logical :: sink_writing = .false.
   !> `.true.` while capturing the value stream for the reproducibility check;
   !> `emit` then neither writes nor compares, it only records
   logical :: sink_capture = .false.
   !> Captured value stream, valid entries `1:sink_count`
   real(wp), allocatable :: sink_vals(:)
   !> Output unit in regeneration mode
   integer :: sink_unit = -1
   !> Records emitted so far (also the cursor into `sink_ref`)
   integer :: sink_count = 0
   !> Parsed fixture, checking mode only
   type(golden_record_type), allocatable :: sink_ref(:)
   !> Number of valid entries in `sink_ref`
   integer :: sink_nref = 0
   !> Number of mismatching records seen
   integer :: sink_nfail = 0
   !> Description of the first mismatch, for the failure message
   character(len=:), allocatable :: sink_message

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

   !> Compare the SvdW LSF against `lsf_golden_svdw.txt`, or rewrite it
   subroutine test_svdw_golden(error)
      type(error_type), allocatable, intent(out) :: error
      call run_golden(error, "svdw")
   end subroutine test_svdw_golden

   !> Compare the CFC LSF against `lsf_golden_cfc.txt`, or rewrite it
   subroutine test_cfc_golden(error)
      type(error_type), allocatable, intent(out) :: error
      call run_golden(error, "cfc")
   end subroutine test_cfc_golden

   !> Drive one concrete over the fixture, in whichever mode the environment
   !> selects. `MOIST_LSF_GOLDEN_REGENERATE` non-empty writes, anything else
   !> checks.
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine run_golden(error, kind)
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      character(len=:), allocatable :: path, regen
      integer :: unit, stat
      character(len=64) :: tail

      path = golden_path(kind)
      regen = get_env("MOIST_LSF_GOLDEN_REGENERATE", default="")

      call reset_sink()

      if (len_trim(regen) > 0) then
         open (newunit=unit, file=path, action="write", status="replace", iostat=stat)
         if (stat /= 0) then
            call test_failed(error, "cannot open '"//path//"' for writing")
            return
         end if
         sink_writing = .true.
         sink_unit = unit
         call write_header(unit, kind)
      else
         call load_fixture(path, error)
         if (allocated(error)) return
         sink_writing = .false.
      end if

      call traverse(kind, error)

      if (sink_writing) then
         close (sink_unit)
         sink_unit = -1
         sink_writing = .false.
         if (sink_nfail > 0) then
            call test_failed(error, "regeneration hit a structural problem, the file "// &
                             "is not trustworthy: "//sink_message)
            return
         end if
         write (error_unit, '(a,i0,a)') "# regenerated "//path//" (", sink_count, " records)"
         return
      end if

      if (allocated(error)) return

      if (sink_count /= sink_nref) then
         write (tail, '(i0,a,i0)') sink_count, " records, fixture has ", sink_nref
         call test_failed(error, "golden fixture length mismatch: evaluated "//trim(tail)// &
                          " - regenerate "//path)
         return
      end if

      if (sink_nfail > 0) then
         write (tail, '(i0)') sink_nfail
         call test_failed(error, trim(tail)//" golden record(s) deviate; first: "//sink_message)
      end if
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
   !> @param[in]  kind   `svdw` or `cfc`
   !> @param[out] error  testdrive failure (setup problems only)
   subroutine traverse(kind, error)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer, allocatable :: sel(:)
      integer :: icase, ncase, ip

      if (kind == "svdw") then
         ncase = n_svdw_cases
      else
         ncase = n_cfc_cases
      end if

      do icase = 1, ncase
         call load_case(golden_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return
         call select_atoms(mol%nat, sel)

         do ip = 1, n_points
            if (kind == "svdw") then
               call svdw_point(golden_cases(icase), mol, radii, points(:, ip), ip, sel, error)
            else
               call cfc_point(golden_cases(icase), mol, radii, points(:, ip), ip, sel, error)
            end if
            if (allocated(error)) return
         end do

         deallocate (radii, points, sel)
      end do
   end subroutine traverse

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
   !*                              SvdW record emission                                 *!
   !* ================================================================================= *!

   !> Emit (or check) every SvdW record of one case at one evaluation point,
   !> screened block first, unscreened block second, difference block last.
   !>
   !> @param[in]  gcase  Case descriptor
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  point  Evaluation point
   !> @param[in]  ip     Evaluation-point index
   !> @param[in]  sel    Selected atom indices
   !> @param[out] error  testdrive failure
   subroutine svdw_point(gcase, mol, radii, point, ip, sel, error)
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
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_svdw_type) :: lsf_scr, lsf_ref
      real(wp) :: f0_scr, f0_ref

      call new_svdw(lsf_scr, gcase, mol, radii, production_threshold)
      call new_svdw(lsf_ref, gcase, mol, radii, 0.0_wp)

      call prepare_svdw(lsf_scr, point, error)
      if (allocated(error)) return
      call prepare_svdw(lsf_ref, point, error)
      if (allocated(error)) return

      call assert_unscreened(lsf_ref%active_count(), mol%nat, gcase%tag, ip, error)
      if (allocated(error)) return

      call svdw_block(lsf_scr, gcase, ip, "S", sel, f0_scr)
      call svdw_block(lsf_ref, gcase, ip, "U", sel, f0_ref)
      call emit("svdw", gcase%tag, ip, "D", "f0_delta", [0, 0, 0, 0, 0, 0], f0_scr - f0_ref)
   end subroutine svdw_point

   !> Construct and bind one SvdW LSF at a given screening threshold.
   !> The threshold must be set before `update`, which is what pushes it into
   !> the SSD system.
   !>
   !> @param[out] lsf    Fresh SvdW LSF
   !> @param[in]  gcase  Case descriptor (supplies the blending weights)
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  thr    Screening threshold
   subroutine new_svdw(lsf, gcase, mol, radii, thr)
      !> Fresh SvdW LSF
      type(moist_cavity_drop_lsf_svdw_type), intent(out) :: lsf
      !> Case descriptor
      type(golden_case_type), intent(in) :: gcase
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Screening threshold
      real(wp), intent(in) :: thr

      lsf%screening_threshold = thr
      call lsf%new(blend_k=gcase%blend_k, blend_1b=gcase%blend_1b, &
                   blend_2b=gcase%blend_2b, blend_3b=gcase%blend_3b)
      call lsf%update(mol, radii)
      call lsf%set_max_deriv(4)
   end subroutine new_svdw

   !> `prepare` an SvdW LSF, translating an LSF error into a testdrive failure
   !>
   !> @param[inout] lsf    SvdW LSF
   !> @param[in]    point  Evaluation point
   !> @param[out]   error  testdrive failure
   subroutine prepare_svdw(lsf, point, error)
      !> SvdW LSF
      type(moist_cavity_drop_lsf_svdw_type), intent(inout) :: lsf
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: lsf_err

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) call test_failed(error, "SvdW prepare failed: "//lsf_err%message)
   end subroutine prepare_svdw

   !> Emit every SvdW quantity of one prepared LSF under one screening flag
   !>
   !> @param[in]  lsf    Prepared SvdW LSF
   !> @param[in]  gcase  Case descriptor
   !> @param[in]  ip     Evaluation-point index
   !> @param[in]  flag   `S` or `U`
   !> @param[in]  sel    Selected atom indices
   !> @param[out] f0     Value returned by `f0`, for the delta record
   subroutine svdw_block(lsf, gcase, ip, flag, sel, f0)
      !> Prepared SvdW LSF
      type(moist_cavity_drop_lsf_svdw_type), intent(in) :: lsf
      !> Case descriptor
      type(golden_case_type), intent(in) :: gcase
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> Screening flag
      character(len=*), intent(in) :: flag
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Value of `f0`
      real(wp), intent(out) :: f0

      real(wp) :: lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim), norm0
      real(wp) :: lsf3_rrr(ndim, ndim, ndim), lsf4_rrrr(ndim, ndim, ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :), lsf4_rrr_rA(:, :, :, :, :)
      real(wp), allocatable :: lsf2_rArB(:, :, :, :), lsf3_r_rArB(:, :, :, :, :)
      real(wp), allocatable :: lsf4_rr_rArB(:, :, :, :, :, :)
      real(wp), allocatable :: norm1_rA(:, :)
      integer, allocatable :: act(:)
      integer :: nat, nac, j, k, l, m, s, t, ia, ib, atomA, atomB, slotA, slotB

      nat = lsf%ncenters
      ! Every nuclear output is active-indexed and caller-sized; the harness
      ! resolves the index space from the extent it gets back (see `nuc_slot`)
      nac = lsf%active_count()
      allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
      allocate (lsf3_rr_rA(ndim, ndim, ndim, nac))
      allocate (lsf4_rrr_rA(ndim, ndim, ndim, ndim, nac))
      allocate (lsf2_rArB(ndim, nac, ndim, nac))
      allocate (lsf3_r_rArB(ndim, ndim, nac, ndim, nac))
      allocate (lsf4_rr_rArB(ndim, ndim, ndim, nac, ndim, nac))
      allocate (norm1_rA(ndim, nac))
      call active_map(lsf, nat, act)

      call emit("svdw", gcase%tag, ip, flag, "n_active", [0, 0, 0, 0, 0, 0], &
                real(lsf%active_count(), wp))

      call lsf%f0(f0)
      call emit("svdw", gcase%tag, ip, flag, "f0", [0, 0, 0, 0, 0, 0], f0)

      call lsf%f012_r(lsf0, lsf1_r, lsf2_rr)
      call emit("svdw", gcase%tag, ip, flag, "f012_lsf0", [0, 0, 0, 0, 0, 0], lsf0)
      do j = 1, ndim
         call emit("svdw", gcase%tag, ip, flag, "f012_lsf1_r", [j, 0, 0, 0, 0, 0], lsf1_r(j))
      end do
      do j = 1, ndim
         do k = j, ndim
            call emit("svdw", gcase%tag, ip, flag, "f012_lsf2_rr", [j, k, 0, 0, 0, 0], &
                      lsf2_rr(j, k))
         end do
      end do

      call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               call emit("svdw", gcase%tag, ip, flag, "f3_rrr", [j, k, l, 0, 0, 0], &
                         lsf3_rrr(j, k, l))
            end do
         end do
      end do

      call lsf%f4_rrrr(lsf4_rrrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               do m = l, ndim
                  call emit("svdw", gcase%tag, ip, flag, "f4_rrrr", [j, k, l, m, 0, 0], &
                            lsf4_rrrr(j, k, l, m))
               end do
            end do
         end do
      end do

      call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      call check_user_space(lsf1_rA, act, "f3rrA_lsf1_rA")
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf1_rA, 2), atomA, act, "f3rrA_lsf1_rA")
         do s = 1, ndim
            call emit("svdw", gcase%tag, ip, flag, "f3rrA_lsf1_rA", [s, atomA, 0, 0, 0, 0], &
                      nuc_read2(lsf1_rA, s, slotA))
         end do
      end do
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf2_r_rA, 3), atomA, act, "f3rrA_lsf2_r_rA")
         do s = 1, ndim
            do j = 1, ndim
               call emit("svdw", gcase%tag, ip, flag, "f3rrA_lsf2_r_rA", &
                         [j, s, atomA, 0, 0, 0], nuc_read3(lsf2_r_rA, j, s, slotA))
            end do
         end do
      end do
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf3_rr_rA, 4), atomA, act, "f3_rr_rA")
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  call emit("svdw", gcase%tag, ip, flag, "f3_rr_rA", &
                            [j, k, s, atomA, 0, 0], nuc_read4(lsf3_rr_rA, j, k, s, slotA))
               end do
            end do
         end do
      end do

      call lsf%f4_rrr_rA(lsf4_rrr_rA)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf4_rrr_rA, 5), atomA, act, "f4_rrr_rA")
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  do l = k, ndim
                     call emit("svdw", gcase%tag, ip, flag, "f4_rrr_rA", &
                               [j, k, l, s, atomA, 0], &
                               nuc_read5(lsf4_rrr_rA, j, k, l, s, slotA))
                  end do
               end do
            end do
         end do
      end do

      call lsf%f2_rArB(lsf2_rArB)
      do ia = 1, min(n_pair_atoms, size(sel))
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf2_rArB, 2), atomA, act, "f2_rArB")
         do ib = 1, min(n_pair_atoms, size(sel))
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf2_rArB, 4), atomB, act, "f2_rArB")
            do s = 1, ndim
               do t = 1, ndim
                  call emit("svdw", gcase%tag, ip, flag, "f2_rArB", &
                            [s, atomA, t, atomB, 0, 0], &
                            pair_read4(lsf2_rArB, s, slotA, t, slotB))
               end do
            end do
         end do
      end do

      call lsf%f3_r_rArB(lsf3_r_rArB)
      do ia = 1, min(n_pair_atoms, size(sel))
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf3_r_rArB, 3), atomA, act, "f3_r_rArB")
         do ib = 1, min(n_pair_atoms, size(sel))
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf3_r_rArB, 5), atomB, act, "f3_r_rArB")
            do s = 1, ndim
               do t = 1, ndim
                  do j = 1, ndim
                     call emit("svdw", gcase%tag, ip, flag, "f3_r_rArB", &
                               [j, s, atomA, t, atomB, 0], &
                               pair_read5(lsf3_r_rArB, j, s, slotA, t, slotB))
                  end do
               end do
            end do
         end do
      end do

      call lsf%f4_rr_rArB(lsf4_rr_rArB)
      do ia = 1, min(n_pair_atoms, size(sel))
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf4_rr_rArB, 4), atomA, act, "f4_rr_rArB")
         do ib = 1, min(n_pair_atoms, size(sel))
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf4_rr_rArB, 6), atomB, act, "f4_rr_rArB")
            do s = 1, ndim
               do t = 1, ndim
                  do j = 1, ndim
                     do k = j, ndim
                        call emit("svdw", gcase%tag, ip, flag, "f4_rr_rArB", &
                                  [j, k, s, atomA, t, atomB], &
                                  pair_read6(lsf4_rr_rArB, j, k, s, slotA, t, slotB))
                     end do
                  end do
               end do
            end do
         end do
      end do

      call lsf%normalized_f01_rA(norm0, norm1_rA)
      call check_user_space(norm1_rA, act, "normalized_f1_rA")
      call emit("svdw", gcase%tag, ip, flag, "normalized_f0", [0, 0, 0, 0, 0, 0], norm0)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(norm1_rA, 2), atomA, act, "normalized_f1_rA")
         do s = 1, ndim
            call emit("svdw", gcase%tag, ip, flag, "normalized_f1_rA", &
                      [s, atomA, 0, 0, 0, 0], nuc_read2(norm1_rA, s, slotA))
         end do
      end do
   end subroutine svdw_block

   !* --------------------------- Index-space-safe reads --------------------------- *!
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
   !> @param[in] extent  Length of the nuclear dimension as returned
   !> @param[in] atom    User-space atom id
   !> @param[in] act     Atom -> active index map (0 = dropped), size ncenters
   !> @param[in] label   Quantity name, for the failure message
   !> @returns           Slot to read, or 0 if the atom has none
   integer function nuc_slot(extent, atom, act, label)
      !> Length of the nuclear dimension as returned
      integer, intent(in) :: extent
      !> User-space atom id
      integer, intent(in) :: atom
      !> Atom -> active index map (0 = dropped)
      integer, intent(in) :: act(:)
      !> Quantity name, for the failure message
      character(len=*), intent(in) :: label

      character(len=64) :: tail

      if (extent == size(act)) then
         nuc_slot = atom
      else if (extent == count(act > 0)) then
         nuc_slot = act(atom)
      else
         write (tail, '(a,i0,a,i0,a,i0)') " extent ", extent, " is neither ncenters ", &
            size(act), " nor active_count ", count(act > 0)
         call flag_problem(label//":"//trim(tail)//" - index space changed in src/")
         nuc_slot = 0
         return
      end if
      if (nuc_slot < 1 .or. nuc_slot > extent) nuc_slot = 0
   end function nuc_slot

   !> Read `t(i1, slotA)`, or 0 when the atom has no slot
   !> @param[in] t      Nuclear array `(axis, A)`
   !> @param[in] i1     Leading index
   !> @param[in] slotA  Slot from [[nuc_slot]]
   !> @returns          Element or 0
   pure function nuc_read2(t, i1, slotA) result(val)
      !> Nuclear array
      real(wp), intent(in) :: t(:, :)
      !> Leading index
      integer, intent(in) :: i1
      !> Slot of A
      integer, intent(in) :: slotA
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0) val = t(i1, slotA)
   end function nuc_read2

   !> Read `t(i1, i2, slotA)`, or 0 when the atom has no slot
   !> @param[in] t      Nuclear array `(axis, axis, A)`
   !> @param[in] i1     First index
   !> @param[in] i2     Second index
   !> @param[in] slotA  Slot from [[nuc_slot]]
   !> @returns          Element or 0
   pure function nuc_read3(t, i1, i2, slotA) result(val)
      !> Nuclear array
      real(wp), intent(in) :: t(:, :, :)
      !> First index
      integer, intent(in) :: i1
      !> Second index
      integer, intent(in) :: i2
      !> Slot of A
      integer, intent(in) :: slotA
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0) val = t(i1, i2, slotA)
   end function nuc_read3

   !> Read `t(i1, i2, i3, slotA)`, or 0 when the atom has no slot
   !> @param[in] t      Nuclear array `(axis, axis, axis, A)`
   !> @param[in] i1     First index
   !> @param[in] i2     Second index
   !> @param[in] i3     Third index
   !> @param[in] slotA  Slot from [[nuc_slot]]
   !> @returns          Element or 0
   pure function nuc_read4(t, i1, i2, i3, slotA) result(val)
      !> Nuclear array
      real(wp), intent(in) :: t(:, :, :, :)
      !> First index
      integer, intent(in) :: i1
      !> Second index
      integer, intent(in) :: i2
      !> Third index
      integer, intent(in) :: i3
      !> Slot of A
      integer, intent(in) :: slotA
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0) val = t(i1, i2, i3, slotA)
   end function nuc_read4

   !> Read `t(i1, i2, i3, i4, slotA)`, or 0 when the atom has no slot
   !> @param[in] t      Nuclear array `(axis, axis, axis, axis, A)`
   !> @param[in] i1     First index
   !> @param[in] i2     Second index
   !> @param[in] i3     Third index
   !> @param[in] i4     Fourth index
   !> @param[in] slotA  Slot from [[nuc_slot]]
   !> @returns          Element or 0
   pure function nuc_read5(t, i1, i2, i3, i4, slotA) result(val)
      !> Nuclear array
      real(wp), intent(in) :: t(:, :, :, :, :)
      !> First index
      integer, intent(in) :: i1
      !> Second index
      integer, intent(in) :: i2
      !> Third index
      integer, intent(in) :: i3
      !> Fourth index
      integer, intent(in) :: i4
      !> Slot of A
      integer, intent(in) :: slotA
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0) val = t(i1, i2, i3, i4, slotA)
   end function nuc_read5

   !> Read `t(i1, slotA, i2, slotB)`, or 0 when either atom has no slot
   !> @param[in] t      Pair array `(axis, A, axis, B)`
   !> @param[in] i1     Nuclear axis of A
   !> @param[in] slotA  Slot of A from [[nuc_slot]]
   !> @param[in] i2     Nuclear axis of B
   !> @param[in] slotB  Slot of B from [[nuc_slot]]
   !> @returns          Element or 0
   pure function pair_read4(t, i1, slotA, i2, slotB) result(val)
      !> Pair array
      real(wp), intent(in) :: t(:, :, :, :)
      !> Nuclear axis of A
      integer, intent(in) :: i1
      !> Slot of A
      integer, intent(in) :: slotA
      !> Nuclear axis of B
      integer, intent(in) :: i2
      !> Slot of B
      integer, intent(in) :: slotB
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0 .and. slotB > 0) val = t(i1, slotA, i2, slotB)
   end function pair_read4

   !> Read `t(i1, i2, slotA, i3, slotB)`, or 0 when either atom has no slot
   !> @param[in] t      Pair array `(axis, axis, A, axis, B)`
   !> @param[in] i1     Spatial axis
   !> @param[in] i2     Nuclear axis of A
   !> @param[in] slotA  Slot of A from [[nuc_slot]]
   !> @param[in] i3     Nuclear axis of B
   !> @param[in] slotB  Slot of B from [[nuc_slot]]
   !> @returns          Element or 0
   pure function pair_read5(t, i1, i2, slotA, i3, slotB) result(val)
      !> Pair array
      real(wp), intent(in) :: t(:, :, :, :, :)
      !> Spatial axis
      integer, intent(in) :: i1
      !> Nuclear axis of A
      integer, intent(in) :: i2
      !> Slot of A
      integer, intent(in) :: slotA
      !> Nuclear axis of B
      integer, intent(in) :: i3
      !> Slot of B
      integer, intent(in) :: slotB
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0 .and. slotB > 0) val = t(i1, i2, slotA, i3, slotB)
   end function pair_read5

   !> Read `t(i1, i2, i3, slotA, i4, slotB)`, or 0 when either atom has no slot
   !> @param[in] t      Pair array `(axis, axis, axis, A, axis, B)`
   !> @param[in] i1     First spatial axis
   !> @param[in] i2     Second spatial axis
   !> @param[in] i3     Nuclear axis of A
   !> @param[in] slotA  Slot of A from [[nuc_slot]]
   !> @param[in] i4     Nuclear axis of B
   !> @param[in] slotB  Slot of B from [[nuc_slot]]
   !> @returns          Element or 0
   pure function pair_read6(t, i1, i2, i3, slotA, i4, slotB) result(val)
      !> Pair array
      real(wp), intent(in) :: t(:, :, :, :, :, :)
      !> First spatial axis
      integer, intent(in) :: i1
      !> Second spatial axis
      integer, intent(in) :: i2
      !> Nuclear axis of A
      integer, intent(in) :: i3
      !> Slot of A
      integer, intent(in) :: slotA
      !> Nuclear axis of B
      integer, intent(in) :: i4
      !> Slot of B
      integer, intent(in) :: slotB
      real(wp) :: val

      val = 0.0_wp
      if (slotA > 0 .and. slotB > 0) val = t(i1, i2, i3, slotA, i4, slotB)
   end function pair_read6

   !> Guard a *caller-sized* nuclear output, whose extent cannot reveal its index
   !> space because the harness chose it.
   !>
   !> `f3_rr_rA`'s optional `lsf1_rA` (and the CFC twin) are passed in as
   !> `(3, ncenters)` buffers and scattered into by user-space atom id. If a
   !> future refactor made them active-indexed instead, the extent check in
   !> [[nuc_slot]] could not see it - but the columns of screened-away atoms would
   !> stop being zero. That is what this asserts.
   !>
   !> @param[in] t      Caller-sized `(axis, A)` output
   !> @param[in] act    Atom -> active index map (0 = dropped)
   !> @param[in] label  Quantity name, for the failure message
   subroutine check_user_space(t, act, label)
      !> Caller-sized nuclear output
      real(wp), intent(in) :: t(:, :)
      !> Atom -> active index map
      integer, intent(in) :: act(:)
      !> Quantity name
      character(len=*), intent(in) :: label

      integer :: atom
      character(len=32) :: tail

      if (size(t, 2) /= size(act)) return
      do atom = 1, size(act)
         if (act(atom) > 0) cycle
         if (maxval(abs(t(:, atom))) == 0.0_wp) cycle
         write (tail, '(a,i0,a)') " atom ", atom, " is screened away"
         call flag_problem(label//":"//trim(tail)//" yet its user-space column is "// &
                           "non-zero - the output is no longer user-indexed")
         return
      end do
   end subroutine check_user_space

   !* ================================================================================= *!
   !*                              CFC record emission                                  *!
   !* ================================================================================= *!

   !> Emit (or check) every CFC record of one case at one evaluation point
   !>
   !> @param[in]  gcase  Case descriptor (blending weights unused)
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  point  Evaluation point
   !> @param[in]  ip     Evaluation-point index
   !> @param[in]  sel    Selected atom indices
   !> @param[out] error  testdrive failure
   subroutine cfc_point(gcase, mol, radii, point, ip, sel, error)
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
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_cfc_type) :: lsf_scr, lsf_ref
      real(wp) :: f0_scr, f0_ref

      call new_cfc(lsf_scr, mol, radii, production_threshold)
      call new_cfc(lsf_ref, mol, radii, 0.0_wp)

      call prepare_cfc(lsf_scr, point, error)
      if (allocated(error)) return
      call prepare_cfc(lsf_ref, point, error)
      if (allocated(error)) return

      call assert_unscreened(lsf_ref%active_count(), mol%nat, gcase%tag, ip, error)
      if (allocated(error)) return

      call cfc_block(lsf_scr, gcase, ip, "S", sel, f0_scr)
      call cfc_block(lsf_ref, gcase, ip, "U", sel, f0_ref)
      call emit("cfc", gcase%tag, ip, "D", "f0_delta", [0, 0, 0, 0, 0, 0], f0_scr - f0_ref)
   end subroutine cfc_point

   !> Construct and bind one CFC LSF at a given screening threshold
   !>
   !> @param[out] lsf    Fresh CFC LSF
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[in]  thr    Screening threshold
   subroutine new_cfc(lsf, mol, radii, thr)
      !> Fresh CFC LSF
      type(moist_cavity_drop_lsf_cfc_type), intent(out) :: lsf
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Screening threshold
      real(wp), intent(in) :: thr

      lsf%screening_threshold = thr
      call lsf%new()
      call lsf%update(mol, radii)
      call lsf%set_max_deriv(3)
   end subroutine new_cfc

   !> `prepare` a CFC LSF, translating an LSF error into a testdrive failure
   !>
   !> @param[inout] lsf    CFC LSF
   !> @param[in]    point  Evaluation point
   !> @param[out]   error  testdrive failure
   subroutine prepare_cfc(lsf, point, error)
      !> CFC LSF
      type(moist_cavity_drop_lsf_cfc_type), intent(inout) :: lsf
      !> Evaluation point
      real(wp), intent(in) :: point(ndim)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: lsf_err

      call lsf%prepare(point, lsf_err)
      if (allocated(lsf_err)) call test_failed(error, "CFC prepare failed: "//lsf_err%message)
   end subroutine prepare_cfc

   !> Emit every CFC quantity of one prepared LSF under one screening flag.
   !> CFC is capped at derivative order 3 today, so there is no f4 block.
   !>
   !> @param[in]  lsf    Prepared CFC LSF
   !> @param[in]  gcase  Case descriptor
   !> @param[in]  ip     Evaluation-point index
   !> @param[in]  flag   `S` or `U`
   !> @param[in]  sel    Selected atom indices
   !> @param[out] f0     Value returned by `f0`, for the delta record
   subroutine cfc_block(lsf, gcase, ip, flag, sel, f0)
      !> Prepared CFC LSF
      type(moist_cavity_drop_lsf_cfc_type), intent(in) :: lsf
      !> Case descriptor
      type(golden_case_type), intent(in) :: gcase
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> Screening flag
      character(len=*), intent(in) :: flag
      !> Selected atom indices
      integer, intent(in) :: sel(:)
      !> Value of `f0`
      real(wp), intent(out) :: f0

      real(wp) :: lsf0, lsf1_r(ndim), lsf2_rr(ndim, ndim)
      real(wp) :: lsf3_rrr(ndim, ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :)
      integer, allocatable :: act(:)
      integer :: nat, nac, j, k, l, s, ia, atomA, slotA

      nat = lsf%ncenters
      nac = lsf%active_count()
      allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
      allocate (lsf3_rr_rA(ndim, ndim, ndim, nac))
      call active_map(lsf, nat, act)

      call emit("cfc", gcase%tag, ip, flag, "n_active", [0, 0, 0, 0, 0, 0], &
                real(lsf%active_count(), wp))

      call lsf%f0(f0)
      call emit("cfc", gcase%tag, ip, flag, "f0", [0, 0, 0, 0, 0, 0], f0)

      call lsf%f012_r(lsf0, lsf1_r, lsf2_rr)
      call emit("cfc", gcase%tag, ip, flag, "f012_lsf0", [0, 0, 0, 0, 0, 0], lsf0)
      do j = 1, ndim
         call emit("cfc", gcase%tag, ip, flag, "f012_lsf1_r", [j, 0, 0, 0, 0, 0], lsf1_r(j))
      end do
      do j = 1, ndim
         do k = j, ndim
            call emit("cfc", gcase%tag, ip, flag, "f012_lsf2_rr", [j, k, 0, 0, 0, 0], &
                      lsf2_rr(j, k))
         end do
      end do

      call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               call emit("cfc", gcase%tag, ip, flag, "f3_rrr", [j, k, l, 0, 0, 0], &
                         lsf3_rrr(j, k, l))
            end do
         end do
      end do

      call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
      call check_user_space(lsf1_rA, act, "cfc f3rrA_lsf1_rA")
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf1_rA, 2), atomA, act, "cfc f3rrA_lsf1_rA")
         do s = 1, ndim
            call emit("cfc", gcase%tag, ip, flag, "f3rrA_lsf1_rA", [s, atomA, 0, 0, 0, 0], &
                      nuc_read2(lsf1_rA, s, slotA))
         end do
      end do
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf2_r_rA, 3), atomA, act, "cfc f3rrA_lsf2_r_rA")
         do s = 1, ndim
            do j = 1, ndim
               call emit("cfc", gcase%tag, ip, flag, "f3rrA_lsf2_r_rA", &
                         [j, s, atomA, 0, 0, 0], nuc_read3(lsf2_r_rA, j, s, slotA))
            end do
         end do
      end do
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf3_rr_rA, 4), atomA, act, "cfc f3_rr_rA")
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  call emit("cfc", gcase%tag, ip, flag, "f3_rr_rA", &
                            [j, k, s, atomA, 0, 0], nuc_read4(lsf3_rr_rA, j, k, s, slotA))
               end do
            end do
         end do
      end do
   end subroutine cfc_block

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

      character(len=64) :: tail

      if (nact == nat) return
      write (tail, '(a,i0,a,i0,a,i0)') " point ", ip, ": active ", nact, " of ", nat
      call test_failed(error, "threshold 0 no longer disables screening for "// &
                       trim(tag)//trim(tail))
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
   !> @param[in]  lsf   Prepared LSF (either concrete)
   !> @param[in]  nat   Number of centres
   !> @param[out] act   `act(atom)` = active index or 0
   subroutine active_map(lsf, nat, act)
      !> Prepared LSF
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Number of centres
      integer, intent(in) :: nat
      !> Atom -> active index
      integer, allocatable, intent(out) :: act(:)

      integer :: i, nact
      character(len=64) :: tail

      nact = lsf%active_count()
      allocate (act(nat), source=0)
      do i = 1, nact
         act(lsf%active_atom(i)) = i
      end do

      if (nact /= nat) return
      do i = 1, nat
         if (act(i) == i) cycle
         write (tail, '(a,i0,a,i0)') " atom ", i, " sits at active slot ", act(i)
         call flag_problem("fully active list is no longer ascending:"//trim(tail))
         return
      end do
   end subroutine active_map

   !* ================================================================================= *!
   !*                                 Sink plumbing                                     *!
   !* ================================================================================= *!

   !> Report a *structural* problem - a shape, an index space or an invariant
   !> that no longer matches what this harness understands - through the same
   !> failure channel as a numerical mismatch.
   !>
   !> Routed here rather than to `test_failed` on purpose: [[run_golden]] checks
   !> `sink_nfail` in regeneration mode too, so a regeneration that trips one of
   !> these fails instead of writing questionable numbers to disk.
   !>
   !> @param[in] text  Message describing the problem
   subroutine flag_problem(text)
      !> Message describing the problem
      character(len=*), intent(in) :: text

      sink_nfail = sink_nfail + 1
      if (sink_nfail == 1) sink_message = text
   end subroutine flag_problem

   !> Clear the module-level sink state between suite entries
   subroutine reset_sink()
      sink_writing = .false.
      sink_capture = .false.
      if (allocated(sink_vals)) deallocate (sink_vals)
      sink_unit = -1
      sink_count = 0
      sink_nfail = 0
      sink_nref = 0
      if (allocated(sink_ref)) deallocate (sink_ref)
      sink_message = ""
   end subroutine reset_sink

   !> Write or check one record. The traversal calls this and nothing else,
   !> so the fixture's order is the traversal's order by construction.
   !>
   !> @param[in] kind      `svdw` or `cfc`
   !> @param[in] case_tag  Case tag
   !> @param[in] ip        Evaluation-point index
   !> @param[in] flag      `S`, `U` or `D`
   !> @param[in] quantity  Quantity name
   !> @param[in] idx       Index tuple, unused slots 0
   !> @param[in] val       Value
   subroutine emit(kind, case_tag, ip, flag, quantity, idx, val)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> Case tag
      character(len=*), intent(in) :: case_tag
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> Screening flag
      character(len=*), intent(in) :: flag
      !> Quantity name
      character(len=*), intent(in) :: quantity
      !> Index tuple
      integer, intent(in) :: idx(6)
      !> Value
      real(wp), intent(in) :: val

      type(golden_record_type) :: ref
      real(wp) :: tol
      character(len=256) :: tail

      sink_count = sink_count + 1

      if (sink_capture) then
         call push_val(val)
         return
      end if

      if (sink_writing) then
         write (sink_unit, record_fmt) kind, case_tag, ip, flag, quantity, idx, val
         return
      end if

      if (sink_count > sink_nref) return

      ref = sink_ref(sink_count)
      if (trim(ref%kind) /= kind .or. trim(ref%case_tag) /= case_tag &
          .or. ref%ip /= ip .or. trim(ref%flag) /= flag &
          .or. trim(ref%quantity) /= quantity .or. any(ref%idx /= idx)) then
         sink_nfail = sink_nfail + 1
         if (sink_nfail == 1) then
            write (tail, '(a,i0,a)') "record ", sink_count, " is labelled '"// &
               trim(ref%kind)//" "//trim(ref%case_tag)//" "//trim(ref%quantity)// &
               "' in the fixture but '"//kind//" "//case_tag//" "//quantity// &
               "' now - regenerate"
            sink_message = trim(tail)
         end if
         return
      end if

      tol = max(golden_abs_tol, golden_rel_tol*abs(ref%val))
      if (abs(val - ref%val) > tol) then
         sink_nfail = sink_nfail + 1
         if (sink_nfail == 1) then
            write (tail, '(a,i0,a,6(1x,i0),a,es24.16,a,es24.16)') &
               "record ", sink_count, " "//trim(ref%kind)//" "//trim(ref%case_tag)// &
               " point "//char(ichar('0') + ip)//" "//trim(ref%flag)//" "// &
               trim(ref%quantity)//" [", ref%idx, " ] golden ", ref%val, " now ", val
            sink_message = trim(tail)
         end if
      end if
   end subroutine emit

   !> Parse a fixture file into `sink_ref`. `#` comments and blank lines are
   !> skipped; everything else must parse as a record.
   !>
   !> @param[in]  path   Fixture path
   !> @param[out] error  testdrive failure
   subroutine load_fixture(path, error)
      !> Fixture path
      character(len=*), intent(in) :: path
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, n
      character(len=256) :: line
      character(len=32) :: tail

      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "missing golden fixture '"//path// &
                          "' - regenerate with MOIST_LSF_GOLDEN_REGENERATE=1")
         return
      end if

      n = 0
      do
         read (unit, '(a)', iostat=stat) line
         if (stat /= 0) exit
         if (is_record_line(line)) n = n + 1
      end do

      allocate (sink_ref(max(n, 1)))
      sink_nref = n
      rewind (unit)

      n = 0
      do
         read (unit, '(a)', iostat=stat) line
         if (stat /= 0) exit
         if (.not. is_record_line(line)) cycle
         n = n + 1
         read (line, *, iostat=stat) sink_ref(n)%kind, sink_ref(n)%case_tag, &
            sink_ref(n)%ip, sink_ref(n)%flag, sink_ref(n)%quantity, &
            sink_ref(n)%idx, sink_ref(n)%val
         if (stat /= 0) then
            close (unit)
            write (tail, '(i0)') n
            call test_failed(error, "malformed record "//trim(tail)//" in "//path)
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

   !> Write the self-describing fixture header
   !>
   !> @param[in] unit  Open output unit
   !> @param[in] kind  `svdw` or `cfc`
   subroutine write_header(unit, kind)
      !> Open output unit
      integer, intent(in) :: unit
      !> Concrete selector
      character(len=*), intent(in) :: kind

      character(len=8)  :: cdate
      character(len=10) :: ctime
      integer :: i, ncase

      call date_and_time(date=cdate, time=ctime)
      if (kind == "svdw") then
         ncase = n_svdw_cases
      else
         ncase = n_cfc_cases
      end if

      write (unit, '(a)') "# moist - DROP level set function golden fixture ("//kind//")"
      write (unit, '(a)') "#"
      write (unit, '(a)') "# Generated by test/unit/test_cavity_drop_lsf_golden.f90."
      write (unit, '(a)') "# Regenerate deliberately, never to make a failing test pass:"
      write (unit, '(a)') "#   env MOIST_LSF_GOLDEN_REGENERATE=1 MOIST_SOURCE_ROOT=$PWD"
      write (unit, '(a)') "#       MOIST_GIT_COMMIT=$(git rev-parse HEAD)"
      write (unit, '(a)') "#       ./<build>/test/unit/tester cavity_drop_lsf_golden"
      write (unit, '(a)') "#"
      write (unit, '(a)') "# git commit : "//get_env("MOIST_GIT_COMMIT", default="(unrecorded)")
      write (unit, '(a)') "# generated  : "//cdate(1:4)//"-"//cdate(5:6)//"-"//cdate(7:8)// &
         " "//ctime(1:2)//":"//ctime(3:4)//":"//ctime(5:6)
      write (unit, '(a)') "#"
      write (unit, '(a)') "# Layout: kind case ip flag quantity i1 i2 i3 i4 i5 i6 value"
      write (unit, '(a)') "#   flag  S = screened as production runs it, U = unscreened"
      write (unit, '(a)') "#         reference, D = screened minus unscreened"
      write (unit, '(a)') "#   i1..i6  unused slots are 0; A/B are user-space atom indices"
      write (unit, '(a)') "#   value   es24.16, exact IEEE-754 double round trip"
      write (unit, '(a)') "# Symmetric slots are dumped once (j<=k<=l<=m etc.); the"
      write (unit, '(a)') "# discarded components are covered by the *_tensor_symmetry tests."
      write (unit, '(a)') "# Full record-layout documentation lives in the module header of"
      write (unit, '(a)') "# test/unit/test_cavity_drop_lsf_golden.f90."
      write (unit, '(a)') "#"
      write (unit, '(a,es12.4)') "# screening_threshold (flag S) : ", production_threshold
      write (unit, '(a,es12.4)') "# screening_threshold (flag U) : ", 0.0_wp
      if (kind == "svdw") then
         write (unit, '(a)') "# max_deriv                    : 4"
      else
         write (unit, '(a)') "# max_deriv                    : 3"
      end if
      write (unit, '(a)') "# radii                        : default_cpcm_radii()"
      write (unit, '(a)') "#"
      write (unit, '(a)') "# case tag      collection   record       blend_k  blend_1b  blend_2b  blend_3b"
      do i = 1, ncase
         if (kind == "svdw") then
            write (unit, '(a,a12,1x,a12,1x,a12,4(1x,f9.4))') "# ", golden_cases(i)%tag, &
               golden_cases(i)%collection, golden_cases(i)%record, golden_cases(i)%blend_k, &
               golden_cases(i)%blend_1b, golden_cases(i)%blend_2b, golden_cases(i)%blend_3b
         else
            write (unit, '(a,a12,1x,a12,1x,a12,a)') "# ", golden_cases(i)%tag, &
               golden_cases(i)%collection, golden_cases(i)%record, "  (CFC compiled defaults)"
         end if
      end do
      write (unit, '(a)') "#"
      write (unit, '(a)') "# point ip  anchor                                offset (bohr)"
      write (unit, '(a)') "#   1 far_out   atom 1 surface, off-axis direction        +13.00"
      write (unit, '(a)') "#   2 near_out  atom 1 surface, off-axis direction         +1.25"
      write (unit, '(a)') "#   3 surface   atom 1 surface, off-axis direction         +0.10"
      write (unit, '(a)') "#   4 deep_in   nucleus 1, off-axis direction        0.40 * R(1)"
      write (unit, '(a)') "#   5 near_nuc  last nucleus, off-axis direction            +0.05"
      write (unit, '(a)') "#"
   end subroutine write_header

   !> Append one value to the captured stream, growing the buffer by doubling
   !>
   !> @param[in] val  Value to append at slot `sink_count`
   subroutine push_val(val)
      !> Value to append
      real(wp), intent(in) :: val

      real(wp), allocatable :: bigger(:)

      if (.not. allocated(sink_vals)) allocate (sink_vals(4096))
      if (sink_count > size(sink_vals)) then
         allocate (bigger(2*size(sink_vals)))
         bigger(1:size(sink_vals)) = sink_vals
         call move_alloc(bigger, sink_vals)
      end if
      sink_vals(sink_count) = val
   end subroutine push_val

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

      call compare_two_passes(error, "svdw")
      if (allocated(error)) return
      call compare_two_passes(error, "cfc")
   end subroutine test_stream_reproducible

   !> Run one concrete's traversal twice and diff the captured streams
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine compare_two_passes(error, kind)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      real(wp), allocatable :: first(:)
      integer :: n_first, i
      character(len=128) :: tail

      call capture_stream(kind, error)
      if (allocated(error)) return
      n_first = sink_count
      allocate (first(n_first), source=sink_vals(1:n_first))

      call capture_stream(kind, error)
      if (allocated(error)) return

      if (sink_count /= n_first) then
         write (tail, '(a,i0,a,i0)') "record count is not reproducible: ", n_first, &
            " then ", sink_count
         call test_failed(error, kind//" "//trim(tail))
         return
      end if

      do i = 1, n_first
         if (first(i) == sink_vals(i)) cycle
         write (tail, '(a,i0,a,es24.16,a,es24.16)') "record ", i, &
            " differs between two passes in one process: ", first(i), " then ", sink_vals(i)
         call test_failed(error, kind//" "//trim(tail))
         return
      end do
   end subroutine compare_two_passes

   !> Traverse one concrete in capture mode, leaving the stream in `sink_vals`
   !>
   !> @param[in]  kind   `svdw` or `cfc`
   !> @param[out] error  testdrive failure
   subroutine capture_stream(kind, error)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      call reset_sink()
      sink_capture = .true.
      call traverse(kind, error)
      if (allocated(error)) return
      if (sink_nfail > 0) call test_failed(error, "structural problem: "//sink_message)
   end subroutine capture_stream

   !* ================================================================================= *!
   !*                              Tensor symmetry checks                               *!
   !* ================================================================================= *!

   !> The golden dump keeps only the unique components of the slots that are
   !> symmetric by construction. This check is the other half of that deal:
   !> it asserts the full tensors really are symmetric in exactly those slots,
   !> so nothing is left unpinned. Tolerance is loose on purpose - the target
   !> is a wrong permutation, not accumulated roundoff.
   subroutine test_svdw_symmetry(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_svdw_type) :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: lsf2_rr(ndim, ndim), lsf1_r(ndim), lsf0
      real(wp) :: lsf3_rrr(ndim, ndim, ndim), lsf4_rrrr(ndim, ndim, ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :), lsf4_rrr_rA(:, :, :, :, :)
      real(wp), allocatable :: lsf4_rr_rArB(:, :, :, :, :, :)
      integer :: icase, ip, j, k, l, m, s, iA, nac

      do icase = 1, n_svdw_cases
         call load_case(golden_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return

         do ip = 1, n_points
            call new_svdw(lsf, golden_cases(icase), mol, radii, production_threshold)
            call prepare_svdw(lsf, points(:, ip), error)
            if (allocated(error)) return
            nac = lsf%active_count()
            if (allocated(lsf1_rA)) deallocate (lsf1_rA, lsf2_r_rA, lsf3_rr_rA, &
                                               lsf4_rrr_rA, lsf4_rr_rArB)
            allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
            allocate (lsf3_rr_rA(ndim, ndim, ndim, nac))
            allocate (lsf4_rrr_rA(ndim, ndim, ndim, ndim, nac))
            allocate (lsf4_rr_rArB(ndim, ndim, ndim, nac, ndim, nac))

            call lsf%f012_r(lsf0, lsf1_r, lsf2_rr)
            do j = 1, ndim
               do k = 1, ndim
                  call check_sym(error, lsf2_rr(j, k), lsf2_rr(k, j), &
                                 maxval(abs(lsf2_rr)), "f012_lsf2_rr")
                  if (allocated(error)) return
               end do
            end do

            call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
            do j = 1, ndim
               do k = 1, ndim
                  do l = 1, ndim
                     call check_sym(error, lsf3_rrr(j, k, l), lsf3_rrr(k, j, l), &
                                    maxval(abs(lsf3_rrr)), "f3_rrr(jk)")
                     if (allocated(error)) return
                     call check_sym(error, lsf3_rrr(j, k, l), lsf3_rrr(j, l, k), &
                                    maxval(abs(lsf3_rrr)), "f3_rrr(kl)")
                     if (allocated(error)) return
                  end do
               end do
            end do

            call lsf%f4_rrrr(lsf4_rrrr)
            do j = 1, ndim
               do k = 1, ndim
                  do l = 1, ndim
                     do m = 1, ndim
                        call check_sym(error, lsf4_rrrr(j, k, l, m), lsf4_rrrr(k, j, l, m), &
                                       maxval(abs(lsf4_rrrr)), "f4_rrrr(jk)")
                        if (allocated(error)) return
                        call check_sym(error, lsf4_rrrr(j, k, l, m), lsf4_rrrr(j, l, k, m), &
                                       maxval(abs(lsf4_rrrr)), "f4_rrrr(kl)")
                        if (allocated(error)) return
                        call check_sym(error, lsf4_rrrr(j, k, l, m), lsf4_rrrr(j, k, m, l), &
                                       maxval(abs(lsf4_rrrr)), "f4_rrrr(lm)")
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do

            call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
            do iA = 1, size(lsf3_rr_rA, 4)
               do s = 1, ndim
                  do j = 1, ndim
                     do k = 1, ndim
                        call check_sym(error, lsf3_rr_rA(j, k, s, iA), lsf3_rr_rA(k, j, s, iA), &
                                       maxval(abs(lsf3_rr_rA)), "f3_rr_rA(jk)")
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do

            call lsf%f4_rrr_rA(lsf4_rrr_rA)
            do iA = 1, size(lsf4_rrr_rA, 5)
               do s = 1, ndim
                  do j = 1, ndim
                     do k = 1, ndim
                        do l = 1, ndim
                           call check_sym(error, lsf4_rrr_rA(j, k, l, s, iA), &
                                          lsf4_rrr_rA(k, j, l, s, iA), &
                                          maxval(abs(lsf4_rrr_rA)), "f4_rrr_rA(jk)")
                           if (allocated(error)) return
                           call check_sym(error, lsf4_rrr_rA(j, k, l, s, iA), &
                                          lsf4_rrr_rA(j, l, k, s, iA), &
                                          maxval(abs(lsf4_rrr_rA)), "f4_rrr_rA(kl)")
                           if (allocated(error)) return
                        end do
                     end do
                  end do
               end do
            end do

            call lsf%f4_rr_rArB(lsf4_rr_rArB)
            do j = 1, ndim
               do k = 1, ndim
                  call check_sym(error, &
                                 maxval(abs(lsf4_rr_rArB(j, k, :, :, :, :) &
                                            - lsf4_rr_rArB(k, j, :, :, :, :))), &
                                 0.0_wp, maxval(abs(lsf4_rr_rArB)), "f4_rr_rArB(jk)")
                  if (allocated(error)) return
               end do
            end do
         end do

         deallocate (radii, points)
      end do
   end subroutine test_svdw_symmetry

   !> CFC counterpart of [[test_svdw_symmetry]] (orders 2 and 3 only)
   subroutine test_cfc_symmetry(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(moist_cavity_drop_lsf_cfc_type) :: lsf
      real(wp), allocatable :: radii(:), points(:, :)
      real(wp) :: lsf2_rr(ndim, ndim), lsf1_r(ndim), lsf0
      real(wp) :: lsf3_rrr(ndim, ndim, ndim)
      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :)
      integer :: icase, ip, j, k, l, s, iA, nac

      do icase = 1, n_cfc_cases
         call load_case(golden_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return

         do ip = 1, n_points
            call new_cfc(lsf, mol, radii, production_threshold)
            call prepare_cfc(lsf, points(:, ip), error)
            if (allocated(error)) return
            nac = lsf%active_count()
            if (allocated(lsf1_rA)) deallocate (lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
            allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
            allocate (lsf3_rr_rA(ndim, ndim, ndim, nac))

            call lsf%f012_r(lsf0, lsf1_r, lsf2_rr)
            do j = 1, ndim
               do k = 1, ndim
                  call check_sym(error, lsf2_rr(j, k), lsf2_rr(k, j), &
                                 maxval(abs(lsf2_rr)), "cfc f012_lsf2_rr")
                  if (allocated(error)) return
               end do
            end do

            call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
            do j = 1, ndim
               do k = 1, ndim
                  do l = 1, ndim
                     call check_sym(error, lsf3_rrr(j, k, l), lsf3_rrr(k, j, l), &
                                    maxval(abs(lsf3_rrr)), "cfc f3_rrr(jk)")
                     if (allocated(error)) return
                     call check_sym(error, lsf3_rrr(j, k, l), lsf3_rrr(j, l, k), &
                                    maxval(abs(lsf3_rrr)), "cfc f3_rrr(kl)")
                     if (allocated(error)) return
                  end do
               end do
            end do

            call lsf%f3_rr_rA(lsf1_rA, lsf2_r_rA, lsf3_rr_rA)
            do iA = 1, size(lsf3_rr_rA, 4)
               do s = 1, ndim
                  do j = 1, ndim
                     do k = 1, ndim
                        call check_sym(error, lsf3_rr_rA(j, k, s, iA), lsf3_rr_rA(k, j, s, iA), &
                                       maxval(abs(lsf3_rr_rA)), "cfc f3_rr_rA(jk)")
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do

         deallocate (radii, points)
      end do
   end subroutine test_cfc_symmetry

   !> Fail unless two tensor elements agree to `symmetry_rel_tol` of the
   !> tensor's largest element
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  a      First element
   !> @param[in]  b      Second element
   !> @param[in]  scale  Largest absolute element of the tensor
   !> @param[in]  label  Name reported on failure
   subroutine check_sym(error, a, b, scale, label)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> First element
      real(wp), intent(in) :: a
      !> Second element
      real(wp), intent(in) :: b
      !> Largest absolute element of the tensor
      real(wp), intent(in) :: scale
      !> Name reported on failure
      character(len=*), intent(in) :: label

      if (abs(a - b) <= max(golden_abs_tol, symmetry_rel_tol*scale)) return
      call test_failed(error, "tensor slot not symmetric: "//label)
   end subroutine check_sym

end module test_cavity_drop_lsf_golden
