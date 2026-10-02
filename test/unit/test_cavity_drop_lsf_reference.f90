!> DROP level-set numerical references for SvdW and CFC
!>
!> Fixture record layout, one record per line:
!>
!>     svdw/cfc case ip flag quantity i1 i2 i3 i4 i5 i6 value
!>
!>   * `case`     [[reference_cases]] tag
!>   * `ip`       evaluation-point index [[build_points]]
!>   * `flag`     `S` production screening, `U` unscreened,
!>                `D` screened minus unscreened
!>   * `i1..i6`   index tuple (see [[nuc_slot]])
!>   * `value`    `es24.16`
!>
!> Symmetric slots stored once (j <= k <= l <= m)
!> Omitted components checked by `*_tensor_symmetry`
module test_cavity_drop_lsf_reference
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use test_helpers, only: get_test_radii, check_moist_error, rel_deviation
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_cfc, only: moist_cavity_drop_lsf_cfc_type
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_lsf_reference

   integer, parameter :: ndim = 3

   !> Maximum derivative order by concrete
   integer, parameter :: max_deriv_svdw = 4
   integer, parameter :: max_deriv_cfc = 3

   !> Screening threshold
   real(wp), parameter :: production_threshold = 1.0e-11_wp

   !> Reference tolerance in rel_deviation
   real(wp), parameter :: reference_tol = 1.0e-12_wp

   !> Tolerances of the companion symmetry checks
   real(wp), parameter :: symmetry_rel_tol = 1.0e-10_wp
   real(wp), parameter :: symmetry_abs_tol = 1.0e-12_wp

   !> Evaluation points per case
   integer, parameter :: n_points = 5
   !> Atoms with single-nucleus derivative records
   integer, parameter :: n_sel_atoms = 3
   !> Atoms with two-nucleus derivative records
   integer, parameter :: n_pair_atoms = 2

   !> Unnormalised offset directions by point
   real(wp), parameter :: raw_dirs(ndim, n_points) = reshape([ &
                          1.0_wp, 2.0_wp, 3.0_wp, &
                          -2.0_wp, 1.0_wp, 4.0_wp, &
                          3.0_wp, -4.0_wp, 1.0_wp, &
                          1.0_wp, -1.0_wp, 2.0_wp, &
                          -1.0_wp, 3.0_wp, -2.0_wp], [ndim, n_points])

   !> Radial offsets along `raw_dirs`; see [[build_points]] for anchors
   !>
   !> Point 5 probes near a nucleus
   real(wp), parameter :: point_offsets(n_points) = &
                          [13.0_wp, 1.25_wp, 0.10_wp, 0.40_wp, 0.05_wp]

   !> Minimum point-to-nucleus distance allowed by [[build_points]]
   !>
   !> Avoid the exact-nucleus branch guarded by `x > 0`
   real(wp), parameter :: min_nucleus_clearance = 4.0e-2_wp

   !> mstore structure and SvdW blending weights for one fixture case
   !>
   !> CFC ignores the weights
   type :: reference_case_type
      !> Case tag in each record
      character(len=12) :: tag
      !> mstore collection
      character(len=12) :: collection
      !> mstore record ID
      character(len=12) :: record
      !> SvdW blending sharpness
      real(wp) :: blend_k
      !> SvdW one-body weight
      real(wp) :: blend_1b
      !> SvdW two-body weight
      real(wp) :: blend_2b
      !> SvdW three-body weight
      real(wp) :: blend_3b
   end type reference_case_type

   !> Structures shared by both LSFs
   integer, parameter :: n_svdw_cases = 6
   integer, parameter :: n_cfc_cases = 5
   type(reference_case_type), parameter :: reference_cases(n_svdw_cases) = [ &
      reference_case_type("lih", "MB16-43", "LiH", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      reference_case_type("ch4", "MB16-43", "CH4", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      reference_case_type("bih3_h2o", "Heavy28", "bih3_h2o", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      reference_case_type("mb16_01", "MB16-43", "01", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      reference_case_type("ala_xab", "Amino20x4", "ALA_xab", 5.5_wp, 1.0_wp, 0.0_wp, 3.0_wp), &
      reference_case_type("ch4_legacy", "MB16-43", "CH4", 3.0_wp, 1.0_wp, 1.0_wp, 1.0_wp)]

   !* ================================================================================= *!
   !*                                  Record stream                                    *!
   !* ================================================================================= *!

   !> Reference entry
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

   !> Record identity and value
   type :: record_type
      !> Concrete, case, evaluation point and screening flag
      type(record_id_type) :: id
      !> Quantity name
      character(len=20) :: quantity
      !> Index tuple, unused slots 0
      integer :: idx(6)
      !> Value
      real(wp) :: val
   end type record_type

   !> Records from one traversal
   !>
   !> Collect shape, index-space, and invariant failures through traversal
   !> Report them alongside numerical mismatches
   type :: record_stream_type
      !> Records in `1:n`
      type(record_type), allocatable :: rec(:)
      !> Record count
      integer :: n = 0
      !> Structural problem count
      integer :: nproblem = 0
      !> First structural problem
      character(len=:), allocatable :: problem
   contains
      !> Append one record
      procedure :: emit => stream_emit
      !> Record a structural problem
      procedure :: flag_problem => stream_flag_problem
   end type record_stream_type

contains

   !> Register the reference-fixture suite
   subroutine collect_cavity_drop_lsf_reference(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("svdw_matches_reference", test_svdw_reference), &
                  new_unittest("cfc_matches_reference", test_cfc_reference), &
                  new_unittest("svdw_tensor_symmetry", test_svdw_symmetry), &
                  new_unittest("cfc_tensor_symmetry", test_cfc_symmetry), &
                  new_unittest("stream_is_reproducible", test_stream_reproducible) &
                  ]
   end subroutine collect_cavity_drop_lsf_reference

   !* ================================================================================= *!
   !*                              Suite entry points                                   *!
   !* ================================================================================= *!

   !> Compare the SvdW LSF against `lsf_reference_svdw.txt`
   subroutine test_svdw_reference(error)
      type(error_type), allocatable, intent(out) :: error
      call run_reference(error, "svdw")
   end subroutine test_svdw_reference

   !> Compare the CFC LSF against `lsf_reference_cfc.txt`
   subroutine test_cfc_reference(error)
      type(error_type), allocatable, intent(out) :: error
      call run_reference(error, "cfc")
   end subroutine test_cfc_reference

   !> Compare one concrete's record stream with its fixture
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine run_reference(error, kind)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      type(record_stream_type) :: got
      type(record_type), allocatable :: ref(:)
      character(len=:), allocatable :: path

      call reference_path(kind, path)

      call traverse(kind, got, error)
      if (allocated(error)) return
      call check_problems(got, error)
      if (allocated(error)) return

      call load_fixture(path, ref, error)
      if (allocated(error)) return

      call compare_stream(got, ref, path, error)
   end subroutine run_reference

   !> Find a fixture reference path
   !>
   !> - Subroutine avoids gfortran's shared deferred-length function result
   !>   during parallel test-drive runs
   !>
   !> @param[in]  kind  `svdw` or `cfc`
   !> @param[out] path  Full path of the fixture
   subroutine reference_path(kind, path)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> Full path of the fixture
      character(len=:), allocatable, intent(out) :: path
      character(len=4096) :: root
      integer :: length, stat

      call get_environment_variable("MOIST_SOURCE_ROOT", root, length, stat)
      if (stat /= 0 .or. length == 0) root = "."
      path = trim(root)//"/test/unit/data/lsf_reference_"//kind//".txt"
   end subroutine reference_path

   !* ================================================================================= *!
   !*                                   Traversal                                       *!
   !* ================================================================================= *!

   !> Emit one concrete's records by case, point, and screening flag
   !>
   !> @param[in]  kind    `svdw` or `cfc`
   !> @param[out] stream  Emitted records
   !> @param[out] error   testdrive failure (setup problems only)
   subroutine traverse(kind, stream, error)
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> Emitted records
      type(record_stream_type), intent(out) :: stream
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      integer, allocatable :: sel(:)
      integer :: icase, ip

      do icase = 1, n_cases(kind)
         call load_case(reference_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return
         call select_atoms(mol%nat, sel)

         do ip = 1, n_points
            call emit_point(kind, reference_cases(icase), mol, radii, points(:, ip), ip, &
                            sel, stream, error)
            if (allocated(error)) return
         end do

         deallocate (radii, points, sel)
      end do
   end subroutine traverse

   !> Fixture case count for one concrete
   !>
   !> @param[in] kind  `svdw` or `cfc`
   !> @returns         Case count
   pure integer function n_cases(kind)
      !> Concrete selector
      character(len=*), intent(in) :: kind

      if (kind == "svdw") then
         n_cases = n_svdw_cases
      else
         n_cases = n_cfc_cases
      end if
   end function n_cases

   !> Load an mstore structure and CPCM radii
   !>
   !> @param[in]  gcase  Case descriptor
   !> @param[out] mol    Structure
   !> @param[out] radii  Per-atom radii (size mol%nat)
   subroutine load_case(gcase, mol, radii)
      !> Case descriptor
      type(reference_case_type), intent(in) :: gcase
      !> Structure
      type(structure_type), intent(out) :: mol
      !> Per-atom radii
      real(wp), allocatable, intent(out) :: radii(:)

      call get_structure(mol, trim(gcase%collection), trim(gcase%record))
      call get_test_radii(mol, radii)
   end subroutine load_case

   !> Build deterministic evaluation points for one structure
   !>
   !> Nucleus-anchored offsets with consistent surface positions:
   !>
   !>   1 `far_out`  atom 1 surface + 13.0 bohr; far atoms screened
   !>                (~13.8 bohr reach at k = 5.5, threshold 1e-11)
   !>   2 `near_out` atom 1 surface + 1.25 bohr; outside, well conditioned
   !>   3 `surface`  atom 1 surface + 0.10 bohr; near zero level set
   !>   4 `deep_in`  0.40 * R1 from nucleus 1; deep inside
   !>   5 `near_nuc` 0.05 bohr from last nucleus; never on nucleus
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
            write (tail, "(i0,a,es12.4)") ip, " sits ", dmin
            call test_failed(error, "evaluation point "//trim(tail)// &
                             " bohr from a nucleus - reference geometry changed?")
            return
         end if
      end do
   end subroutine build_points

   !> Select first, middle, and last atoms in ascending order
   !>
   !> Deduplicate for small structures
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

   !> Create and bind an LSF at its derivative cap
   !>
   !> Set screening threshold before `update` propagates it to SSD
   !> Shared base-class API after concrete construction
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
      type(reference_case_type), intent(in) :: gcase
      !> Structure
      type(structure_type), intent(in) :: mol
      !> Per-atom radii
      real(wp), intent(in) :: radii(:)
      !> Screening threshold
      real(wp), intent(in) :: thr

      integer :: max_deriv

      select case (kind)
      case ("svdw")
         allocate (moist_cavity_drop_lsf_svdw_type :: lsf)
         select type (lsf)
         type is (moist_cavity_drop_lsf_svdw_type)
            lsf%screening_threshold = thr
            call lsf%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=gcase%blend_k, &
               blend_1b=gcase%blend_1b, blend_2b=gcase%blend_2b, blend_3b=gcase%blend_3b))
         end select
         max_deriv = max_deriv_svdw
      case ("cfc")
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

   !> Prepare an LSF and convert errors to testdrive failures
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
      if (allocated(error)) return

      ! Repeated preparation must reset the cached derivative accumulators
      call lsf%prepare(point, err)
      call check_moist_error(error, err, "Repeated LSF prepare failed")
   end subroutine prepare_lsf

   !* ================================================================================= *!
   !*                                Record emission                                    *!
   !* ================================================================================= *!

   !> Emit screened, unscreened, then difference records at one point
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
      type(reference_case_type), intent(in) :: gcase
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
      type(record_stream_type), intent(inout) :: stream
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

   !> Emit all quantities of a prepared LSF for one screening flag
   !>
   !> Shared spatial and single-nucleus blocks
   !> SvdW-only pair, normalisation, and order-4 blocks
   !> CFC capped at order 3; unprepared accessors fail
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
      type(record_stream_type), intent(inout) :: stream
      !> Value of `f0`
      real(wp), intent(out) :: f0

      !> Atom -> active index map (0 = screened away)
      integer, allocatable :: act(:)

      call active_map(lsf, act, stream)

      call emit_spatial_block(lsf, id, stream, f0)
      call emit_nuclear_block(lsf, id, sel, act, stream)
      if (id%kind /= "svdw") return
      call emit_pair_block(lsf, id, sel, act, stream)
      call emit_normalized_block(lsf, id, sel, act, stream)
   end subroutine emit_block

   !> Emit spatial derivatives and active count for nuclear records
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
      type(record_stream_type), intent(inout) :: stream
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

      if (id%kind /= "svdw") return

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

   !> Emit single-nucleus derivatives
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
      type(record_stream_type), intent(inout) :: stream

      real(wp), allocatable :: lsf1_rA(:, :), lsf2_r_rA(:, :, :)
      real(wp), allocatable :: lsf3_rr_rA(:, :, :, :), lsf4_rrr_rA(:, :, :, :, :)
      integer :: nac, j, k, l, s, ia, atomA, slotA

      ! Nuclear outputs use active indices and caller-sized buffers
      ! Resolve returned index space from extent (see `nuc_slot`)
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

      if (id%kind /= "svdw") return

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

   !> Emit two-nucleus derivatives
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
      type(record_stream_type), intent(inout) :: stream

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

   !> Emit surface-normalised value and nuclear gradient
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
      type(record_stream_type), intent(inout) :: stream

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

   !> Pack up to six indices; zero-fill unused slots
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
   !> Nuclear-array reads by user-space atom ID
   !> Derive read slots from returned extents, independent of source index space

   !> Resolve a user-space atom's slot in a nuclear dimension
   !>
   !> - `extent == ncenters`: user index = atom ID
   !> - `extent == active_count()`: active index, or 0 if screened
   !> - Other extent: structural failure and slot 0; no out-of-bounds read
   !>
   !> Equal extents coincide under [[active_map]]'s ascending-order check
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
      type(record_stream_type), intent(inout) :: stream

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

   !> Read a nuclear slot, or 0 for a screened atom
   !>
   !> Callers reduce fixed indices first, e.g. `t(j, k, s, :)`
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

   !> Read two nuclear slots, or 0 if either atom is screened
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

   !> Guard caller-sized nuclear output against index-space changes
   !>
   !> `f3_rr_rA`'s `lsf1_rA` and `normalized_f01_rA`'s gradient use
   !> `(3, ncenters)` buffers keyed by user-space atom ID
   !> Screened columns must remain zero; [[nuc_slot]] cannot infer this
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
      type(record_stream_type), intent(inout) :: stream

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

   !> Require all centers active at `screening_threshold = 0`
   !>
   !> Preserve the unscreened reference across SSD changes
   !>
   !> @param[in]  nact   Active count reported by the LSF
   !> @param[in]  nat    Number of centers
   !> @param[in]  tag    Case tag, for the message
   !> @param[in]  ip     Evaluation-point index, for the message
   !> @param[out] error  testdrive failure
   subroutine assert_unscreened(nact, nat, tag, ip, error)
      !> Active count
      integer, intent(in) :: nact
      !> Number of centers
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

   !> Build atom-to-active-index map (0 = dropped)
   !>
   !> Assert ascending user-space order when all atoms are active
   !> Equal user and active indices required by [[nuc_slot]]
   !> Full `prepare` preserves order through `orig_to_sorted`
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
      type(record_stream_type), intent(inout) :: stream

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

   !> Append one record in traversal order
   !>
   !> @param[inout] self      Record sink
   !> @param[in]    id        Concrete, case, point and screening flag
   !> @param[in]    quantity  Quantity name
   !> @param[in]    idx       Index tuple, unused slots 0
   !> @param[in]    val       Value
   subroutine stream_emit(self, id, quantity, idx, val)
      !> Record sink
      class(record_stream_type), intent(inout) :: self
      !> Record identity
      type(record_id_type), intent(in) :: id
      !> Quantity name
      character(len=*), intent(in) :: quantity
      !> Index tuple
      integer, intent(in) :: idx(6)
      !> Value
      real(wp), intent(in) :: val

      type(record_type), allocatable :: bigger(:)

      if (.not. allocated(self%rec)) allocate (self%rec(4096))
      if (self%n == size(self%rec)) then
         allocate (bigger(2*size(self%rec)))
         bigger(1:size(self%rec)) = self%rec
         call move_alloc(bigger, self%rec)
      end if

      self%n = self%n + 1
      self%rec(self%n) = record_type(id, quantity, idx, val)
   end subroutine stream_emit

   !> Record an unsupported shape, index space, or invariant
   !>
   !> Finish traversal; report any problems before numerical checks
   !>
   !> @param[inout] self  Record sink
   !> @param[in]    text  Message describing the problem
   subroutine stream_flag_problem(self, text)
      !> Record sink
      class(record_stream_type), intent(inout) :: self
      !> Message describing the problem
      character(len=*), intent(in) :: text

      self%nproblem = self%nproblem + 1
      if (self%nproblem == 1) self%problem = text
   end subroutine stream_flag_problem

   !> Convert a structural problem to a testdrive failure
   !>
   !> @param[in]  stream  Traversed stream
   !> @param[out] error   testdrive failure, allocated only if a problem was seen
   subroutine check_problems(stream, error)
      !> Traversed stream
      type(record_stream_type), intent(in) :: stream
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      if (stream%nproblem == 0) return
      call test_failed(error, itoa(stream%nproblem)//" structural problem(s); first: "// &
                       stream%problem)
   end subroutine check_problems

   !> Compare emitted records with the parsed fixture
   !>
   !> @param[in]  got    Emitted stream
   !> @param[in]  ref    Parsed fixture
   !> @param[in]  path   Fixture path, for the failure message
   !> @param[out] error  testdrive failure
   subroutine compare_stream(got, ref, path, error)
      !> Emitted stream
      type(record_stream_type), intent(in) :: got
      !> Parsed fixture
      type(record_type), intent(in) :: ref(:)
      !> Fixture path
      character(len=*), intent(in) :: path
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      integer :: i, nfail
      character(len=:), allocatable :: first
      character(len=64) :: values

      if (got%n /= size(ref)) then
         call test_failed(error, "reference fixture length mismatch: evaluated "// &
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
         if (rel_deviation(got%rec(i)%val, ref(i)%val) <= reference_tol) cycle
         nfail = nfail + 1
         if (nfail == 1) then
            write (values, "(a,es24.16,a,es24.16)") " reference ", ref(i)%val, " now ", &
               got%rec(i)%val
            first = "record "//itoa(i)//" "//record_label(ref(i))//values
         end if
      end do

      if (nfail == 0) return
      call test_failed(error, itoa(nfail)//" record(s) deviate from the reference; first: "//first)
   end subroutine compare_stream

   !> `.true.` when two records name the same number
   !>
   !> @param[in] a  First record
   !> @param[in] b  Second record
   !> @returns      Whether the labels agree
   pure logical function same_label(a, b)
      !> First record
      type(record_type), intent(in) :: a
      !> Second record
      type(record_type), intent(in) :: b

      same_label = a%id%kind == b%id%kind .and. a%id%case_tag == b%id%case_tag &
                   .and. a%id%ip == b%id%ip .and. a%id%flag == b%id%flag &
                   .and. a%quantity == b%quantity .and. all(a%idx == b%idx)
   end function same_label

   !> Left-justified decimal form of `n`
   !>
   !> @param[in] n  Number to render
   !> @returns      Decimal digits, no padding
   function itoa(n) result(text)
      !> Number to render
      integer, intent(in) :: n
      character(len=:), allocatable :: text
      character(len=32) :: buf

      write (buf, "(i0)") n
      text = trim(buf)
   end function itoa

   !> Record identity for failure messages
   !>
   !> @param[in] rec  Record to describe
   !> @returns        One-line description
   function record_label(rec) result(text)
      !> Record to describe
      type(record_type), intent(in) :: rec
      character(len=:), allocatable :: text
      integer :: i

      text = trim(rec%id%kind)//" "//trim(rec%id%case_tag)//" point "//itoa(rec%id%ip)// &
             " "//rec%id%flag//" "//trim(rec%quantity)//" ["
      do i = 1, size(rec%idx)
         text = text//" "//itoa(rec%idx(i))
      end do
      text = text//" ]"
   end function record_label

   !> Parse fixture records, skipping blank and `#` lines
   !>
   !> @param[in]  path   Fixture path
   !> @param[out] ref    Parsed records
   !> @param[out] error  testdrive failure
   subroutine load_fixture(path, ref, error)
      !> Fixture path
      character(len=*), intent(in) :: path
      !> Parsed records
      type(record_type), allocatable, intent(out) :: ref(:)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      integer :: unit, stat, n
      character(len=256) :: line

      open (newunit=unit, file=path, action="read", status="old", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "missing reference fixture '"//path// &
                          "' - it is committed data, restore it from git")
         return
      end if

      n = 0
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (is_record_line(line)) n = n + 1
      end do

      allocate (ref(n))
      rewind (unit)

      n = 0
      do
         read (unit, "(a)", iostat=stat) line
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

   !> `.true.` for a nonblank, noncomment record line
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

   !> Require bit-identical streams from two traversals
   !>
   !> Detect stale memory, races, accumulator state, and NaNs
   !> Past failure: pair-tensor read beyond `n_active` caused unstable fixtures
   !> [[nuc_slot]] now bounds reads by returned extent
   !>
   !> Second pass reuses churned heap to expose uninitialised reads
   subroutine test_stream_reproducible(error)
      type(error_type), allocatable, intent(out) :: error

      call compare_two_passes(error, "svdw")
      if (allocated(error)) return
      call compare_two_passes(error, "cfc")
   end subroutine test_stream_reproducible

   !> Compare two traversals of one concrete
   !>
   !> @param[out] error  testdrive failure
   !> @param[in]  kind   `svdw` or `cfc`
   subroutine compare_two_passes(error, kind)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error
      !> Concrete selector
      character(len=*), intent(in) :: kind

      type(record_stream_type) :: first, second
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
         write (tail, "(a,es24.16,a,es24.16)") &
            " differs between two passes in one process: ", first%rec(i)%val, " then ", &
            second%rec(i)%val
         call test_failed(error, kind//" record "//itoa(i)//trim(tail))
         return
      end do
   end subroutine compare_two_passes

   !* ================================================================================= *!
   !*                              Tensor symmetry checks                               *!
   !* ================================================================================= *!

   !> Check full-tensor symmetry for fixture slots stored once
   !>
   !> Loose tolerance detects wrong permutations over roundoff
   subroutine test_svdw_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      call run_symmetry(error, "svdw")
   end subroutine test_svdw_symmetry

   !> Check CFC symmetry at orders 2 and 3
   subroutine test_cfc_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      call run_symmetry(error, "cfc")
   end subroutine test_cfc_symmetry

   !> Check symmetric slots for every case and point
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
         call load_case(reference_cases(icase), mol, radii)
         call build_points(mol, radii, points, error)
         if (allocated(error)) return

         do ip = 1, n_points
            call new_lsf(lsf, kind, reference_cases(icase), mol, radii, production_threshold)
            call prepare_lsf(lsf, points(:, ip), error)
            if (allocated(error)) return

            call check_symmetry_low(lsf, error)
            if (allocated(error)) return
            if (kind /= "svdw") cycle
            call check_symmetry_high(lsf, error)
            if (allocated(error)) return
         end do

         deallocate (radii, points)
      end do
   end subroutine run_symmetry

   !> Check symmetric slots through order 3
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

   !> Check SvdW order-4 tensor symmetry
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

   !> Check tensor asymmetry against scaled tolerance
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

end module test_cavity_drop_lsf_reference
