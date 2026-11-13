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
   use moist_math_packed, only: packed_index, packed_index3, packed_index4, pack_sym_matrix
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none(type, external)
   private

   public :: collect_cavity_drop_lsf_reference

   integer, parameter :: ndim = 3

   !> Maximum derivative order by lsf
   integer, parameter :: max_deriv_svdw = 4
   integer, parameter :: max_deriv_cfc = 3

   !> Screening threshold
   real(wp), parameter :: production_threshold = 1.0e-11_wp

   !> Reference tolerance in rel_deviation, by lsf
   real(wp), parameter :: svdw_reference_tol = 1.0e-11_wp
   real(wp), parameter :: cfc_reference_tol = 5.0e-13_wp
   !> CFC `evaluate_batch` (strip kernel) against the per-point records
   real(wp), parameter :: cfc_batch_reference_tol = 1.0e-11_wp

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
      !> Fade of the traversed instances recorded ([[record_fade]])
      logical :: fade_set = .false.
      !> Fade switch of the traversed instances
      logical :: fade_enabled = .false.
      !> Fade onset and width in e-folds of the traversed instances
      real(wp) :: fade_band(2) = 0.0_wp
      !> CFC `a1`, `a2`, `c` of the traversed instances
      real(wp) :: cfc_shape(3) = 0.0_wp
      !> CFC pair power `m` of the traversed instances
      integer :: cfc_m = 0
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
                  new_unittest("cfc_batch_matches_reference", test_cfc_batch_reference), &
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

      call write_requested_fixture(got, kind, error)
      if (allocated(error)) return

      call load_fixture(path, ref, error)
      if (allocated(error)) return

      if (kind == "svdw") then
         call compare_stream(got, ref, path, svdw_reference_tol, error)
      else
         call compare_stream(got, ref, path, cfc_reference_tol, error)
      end if
   end subroutine run_reference

   !> CFC `evaluate_batch` at order 4 against the per-point fixture records
   !>
   !> - every CFC case and both screening flags (`S` at the production threshold,
   !>   `U` threshold 0; the same reach with the default fade), the five points of
   !>   a case in one batch, every column reached
   !> - the fixture's `f0`, `f012_lsf0`, `f012_lsf1_r`, `f012_lsf2_rr` and `f3_rrr`
   !>   records against the batch outputs in `rel_deviation`, at
   !>   [[cfc_batch_reference_tol]]; the batch's fourth derivative has no fixture
   !>   record (CFC's fixture stops at order 3)
   !> - every such fixture record compared exactly once
   !>
   !> @param[out] error  testdrive failure
   subroutine test_cfc_batch_reference(error)
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      !> Spatial records per point: f0, f012_lsf0, 3 gradient, 6 Hessian, 10 third
      integer, parameter :: records_per_point = 21
      type(record_type), allocatable :: ref(:)
      character(len=:), allocatable :: path, first
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:), points(:, :)
      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      real(wp) :: b0(n_points), b1(ndim, n_points), b2(6, n_points)
      real(wp) :: b3(10, n_points), b4(15, n_points), got
      integer :: icase, iflag, i, ncmp, nfail
      character(len=1) :: flag
      character(len=80) :: values
      type(mctc_error), allocatable :: err

      call reference_path("cfc", path)
      call load_fixture(path, ref, error)
      if (allocated(error)) return

      ncmp = 0
      nfail = 0
      first = ""
      do icase = 1, n_cfc_cases
         call load_case(reference_cases(icase), mol, radii)
         call build_points("cfc", mol, radii, points, error)
         if (allocated(error)) return
         do iflag = 1, 2
            flag = merge("S", "U", iflag == 1)
            call new_lsf(lsf, "cfc", reference_cases(icase), mol, radii, &
                         merge(production_threshold, 0.0_wp, iflag == 1))
            call lsf%evaluate_batch(points, 4, b0, b1, b2, b3, b4, error=err)
            call check_moist_error(error, err, "CFC evaluate_batch failed")
            if (allocated(error)) return
            do i = 1, n_points
               if (.not. lsf%unreached(i)) cycle
               call test_failed(error, trim(reference_cases(icase)%tag)//" "//flag// &
                                ": a fixture point is unreached in the batch")
               return
            end do
            do i = 1, size(ref)
               associate (r => ref(i))
                  if (r%id%kind /= "cfc" .or. r%id%flag /= flag .or. &
                      r%id%case_tag /= reference_cases(icase)%tag) cycle
                  select case (trim(r%quantity))
                  case ("f0", "f012_lsf0")
                     got = b0(r%id%ip)
                  case ("f012_lsf1_r")
                     got = b1(r%idx(1), r%id%ip)
                  case ("f012_lsf2_rr")
                     got = b2(packed_index(r%idx(1), r%idx(2)), r%id%ip)
                  case ("f3_rrr")
                     got = b3(packed_index3(r%idx(1), r%idx(2), r%idx(3)), r%id%ip)
                  case default
                     cycle
                  end select
                  ncmp = ncmp + 1
                  if (rel_deviation(got, r%val) <= cfc_batch_reference_tol) cycle
                  nfail = nfail + 1
                  if (nfail == 1) then
                     write (values, "(a,es24.16,a,es24.16)") " reference ", r%val, " batch ", got
                     first = record_label(r)//values
                  end if
               end associate
            end do
            deallocate (lsf)
         end do
         deallocate (radii, points)
      end do

      if (ncmp /= 2*n_cfc_cases*n_points*records_per_point) then
         call test_failed(error, "CFC batch reference: compared "//itoa(ncmp)// &
                          " records, expected "//itoa(2*n_cfc_cases*n_points*records_per_point))
         return
      end if
      if (nfail == 0) return
      call test_failed(error, itoa(nfail)//" CFC batch value(s) deviate from the per-point "// &
                       "fixture; first: "//first)
   end subroutine test_cfc_batch_reference

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

   !> Write a record stream to the path its writer variable names, if any
   !>
   !> - `MOIST_LSF_REFERENCE_WRITE` (SvdW) or `MOIST_LSF_REFERENCE_WRITE_CFC`
   !>   (CFC): one variable per concrete, as both tests may run in parallel
   !> - fixture format of `lsf_reference_<kind>.txt`: header (thresholds, order,
   !>   radii, fade band and switch, case table, point table), then one record per
   !>   line `(a4, 1x, a12, i1, 1x, a1, a19 right-justified, 6i4, 1x, es24.16)`
   !> - fade band and switch, and the CFC `a1, a2, c, m` and `far_out` offset,
   !>   from the traversed instances ([[record_fade]])
   !> - unset or empty: nothing written
   !>
   !> @param[in]  stream Traversed stream
   !> @param[in]  kind   `svdw` or `cfc`
   !> @param[out] error  testdrive failure on an I/O error or no recorded fade
   subroutine write_requested_fixture(stream, kind, error)
      !> Traversed stream
      type(record_stream_type), intent(in) :: stream
      !> Concrete selector
      character(len=*), intent(in) :: kind
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      type(reference_case_type) :: gcase
      character(len=4096) :: path
      character(len=19) :: qty
      character(len=4) :: tag
      character(len=:), allocatable :: env, reads
      integer :: length, stat, unit, i, icase, max_deriv

      if (kind == "svdw") then
         env = "MOIST_LSF_REFERENCE_WRITE"
         max_deriv = max_deriv_svdw
         reads = " (SvdW reads no screening_threshold then)"
      else
         env = "MOIST_LSF_REFERENCE_WRITE_CFC"
         max_deriv = max_deriv_cfc
         reads = " (CFC reads no screening_threshold then)"
      end if
      call get_environment_variable(env, path, length, stat)
      if (stat /= 0 .or. length == 0) return
      if (.not. stream%fade_set) then
         call test_failed(error, "fixture writer: no "//kind//" instance recorded its fade")
         return
      end if

      open (newunit=unit, file=trim(path), status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call test_failed(error, "cannot open fixture output '"//trim(path)//"'")
         return
      end if
      write (unit, "(a/a/a,es12.4/a,es12.4/a,i0/a/a,2f8.3/a,l1,a/a)", iostat=stat) &
         "# moist - DROP level set function reference fixture ("//kind//")", "#", &
         "# screening_threshold (flag S) : ", production_threshold, &
         "# screening_threshold (flag U) : ", 0.0_wp, &
         "# max_deriv                    : ", max_deriv, &
         "# radii                        : default_cpcm_radii()", &
         "# fade (a, b e-folds)          : ", stream%fade_band, &
         "# fade enabled                 : ", stream%fade_enabled, reads, "#"
      if (stat == 0 .and. kind == "svdw") then
         write (unit, "(a)", iostat=stat) &
            "# case tag      collection   record       blend_k  blend_1b  blend_2b  blend_3b"
         do icase = 1, n_svdw_cases
            if (stat /= 0) exit
            gcase = reference_cases(icase)
            write (unit, "(a,a12,1x,a12,1x,a12,4f10.4)", iostat=stat) "# ", gcase%tag, &
               gcase%collection, gcase%record, gcase%blend_k, gcase%blend_1b, gcase%blend_2b, &
               gcase%blend_3b
         end do
      else if (stat == 0) then
         write (unit, "(a)", iostat=stat) &
            "# case tag      collection   record             a1        a2         c   m"
         do icase = 1, n_cfc_cases
            if (stat /= 0) exit
            gcase = reference_cases(icase)
            write (unit, "(a,a12,1x,a12,1x,a12,3f10.4,i4)", iostat=stat) "# ", gcase%tag, &
               gcase%collection, gcase%record, stream%cfc_shape, stream%cfc_m
         end do
      end if
      if (stat == 0) then
         write (unit, "(a)", iostat=stat) "#", &
            "# point ip  anchor                                offset (bohr)"
      end if
      if (stat == 0 .and. kind == "svdw") then
         write (unit, "(a)", iostat=stat) &
            "#   1 far_out   atom 1 surface, off-axis direction        +13.00"
      else if (stat == 0) then
         write (unit, "(a,f6.4,a)", iostat=stat) &
            "#   1 far_out   atom 1 surface, off-axis direction   +", &
            cfc_far_sigma(stream%fade_band, stream%cfc_shape(2)), " * R(1)"
      end if
      if (stat == 0) then
         write (unit, "(a)", iostat=stat) &
            "#   2 near_out  atom 1 surface, off-axis direction         +1.25", &
            "#   3 surface   atom 1 surface, off-axis direction         +0.10", &
            "#   4 deep_in   nucleus 1, off-axis direction        0.40 * R(1)", &
            "#   5 near_nuc  last nucleus, off-axis direction            +0.05", &
            "#"
      end if
      do i = 1, stream%n
         if (stat /= 0) exit
         associate (r => stream%rec(i))
            tag = adjustr(r%id%kind(1:4))
            qty = adjustr(r%quantity(1:19))
            write (unit, "(a4,1x,a12,i1,1x,a1,a19,6i4,1x,es24.16)", iostat=stat) tag, &
               r%id%case_tag, r%id%ip, r%id%flag, qty, r%idx, r%val
         end associate
      end do
      if (stat /= 0) then
         close (unit)
         call test_failed(error, "write to fixture output '"//trim(path)//"' failed")
         return
      end if
      close (unit, iostat=stat)
      if (stat /= 0) call test_failed(error, "close of fixture output '"//trim(path)//"' failed")
   end subroutine write_requested_fixture

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
         call build_points(kind, mol, radii, points, error)
         if (allocated(error)) return
         call select_atoms(mol%nat, sel)

         do ip = 1, n_points
            call emit_point(kind, reference_cases(icase), mol, radii, points(:, ip), &
                            points(:, mod(ip, n_points) + 1), ip, sel, stream, error)
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
   !>   1 `far_out`  SvdW: atom 1 surface + 13.0 bohr; at k = 5.5 inside atom 1's
   !>                fade band (onset 11.45, end 13.09 bohr), atoms beyond the band
   !>                end screened; k = 3 (`ch4_legacy`): plateau (onset 21 bohr);
   !>                CFC: atom 1 surface + [[cfc_far_sigma]] R(1), the middle of
   !>                atom 1's default band (sigma 3..3.33), atoms beyond the band
   !>                end screened
   !>   2 `near_out` atom 1 surface + 1.25 bohr; outside, well conditioned
   !>   3 `surface`  atom 1 surface + 0.10 bohr; near zero level set
   !>   4 `deep_in`  0.40 * R1 from nucleus 1; deep inside
   !>   5 `near_nuc` 0.05 bohr from last nucleus; never on nucleus
   !>
   !> @param[in]  kind   `svdw` or `cfc`
   !> @param[in]  mol    Structure
   !> @param[in]  radii  Per-atom radii
   !> @param[out] points (3, n_points) evaluation points
   !> @param[out] error  testdrive failure if a point lands on a nucleus
   subroutine build_points(kind, mol, radii, points, error)
      !> Concrete selector
      character(len=*), intent(in) :: kind
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
      !> CFC instance as [[new_lsf]] constructs it, for the band of `far_out`
      type(moist_cavity_drop_lsf_cfc_type) :: cfc

      do ip = 1, n_points
         dir(:, ip) = raw_dirs(:, ip)/norm2(raw_dirs(:, ip))
      end do

      allocate (points(ndim, n_points))
      if (kind == "cfc") then
         call cfc%new()
         points(:, 1) = mol%xyz(:, 1) + radii(1)*(1.0_wp + &
                        cfc_far_sigma([cfc%param%fade_start, cfc%param%fade_width], &
                                      cfc%param%a2))*dir(:, 1)
      else
         points(:, 1) = mol%xyz(:, 1) + (radii(1) + point_offsets(1))*dir(:, 1)
      end if
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

   !> CFC `far_out` offset off atom 1's surface in units of R(1) (sigma)
   !>
   !> - the middle of the band, `(a + b/2)/|a2|`; 27/3 e-folds at `a2 = -9`:
   !>   sigma 3..3.33, 13 Bohr lies past the reach of most first atoms
   !> - one rounding, `28.5/9` bit for bit at the defaults
   !>
   !> @param[in] band  Fade onset and width in e-folds
   !> @param[in] a2    Pair-term exponent
   !> @returns   sigma Offset in units of R(1)
   pure function cfc_far_sigma(band, a2) result(sigma)
      !> Fade onset and width in e-folds
      real(wp), intent(in) :: band(2)
      !> Pair-term exponent
      real(wp), intent(in) :: a2
      !> Offset in units of R(1)
      real(wp) :: sigma

      sigma = (band(1) + 0.5_wp*band(2))/abs(a2)
   end function cfc_far_sigma

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
   end subroutine prepare_lsf

   !* ================================================================================= *!
   !*                                Record emission                                    *!
   !* ================================================================================= *!

   !> Emit screened, unscreened, then difference records at one point
   !>
   !> The screened LSF is prepared at `prior` first, so its records also catch
   !> state left over from a previous point (accumulators, active-set data)
   !> The unscreened LSF is prepared once on a fresh object
   !>
   !> - fade on (SvdW and CFC fixtures): both flags share the reach (band end),
   !>   so the active counts must agree and the two blocks must be bit-identical
   !>   ([[check_flags_identical]]); fade off: threshold 0 must leave every atom
   !>   active
   !>
   !> @param[in]    kind    `svdw` or `cfc`
   !> @param[in]    gcase   Case descriptor
   !> @param[in]    mol     Structure
   !> @param[in]    radii   Per-atom radii
   !> @param[in]    point   Evaluation point
   !> @param[in]    prior   Different point the screened LSF is prepared at first
   !> @param[in]    ip      Evaluation-point index
   !> @param[in]    sel     Selected atom indices
   !> @param[inout] stream  Record sink
   !> @param[out]   error   testdrive failure
   subroutine emit_point(kind, gcase, mol, radii, point, prior, ip, sel, stream, error)
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
      !> Point the screened LSF is prepared at first
      real(wp), intent(in) :: prior(ndim)
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
      !> Stream length before the `S` and before the `U` block
      integer :: n0, n1

      call new_lsf(lsf_scr, kind, gcase, mol, radii, production_threshold)
      call new_lsf(lsf_ref, kind, gcase, mol, radii, 0.0_wp)
      call record_fade(stream, lsf_scr, gcase%tag, error)
      if (allocated(error)) return
      call record_fade(stream, lsf_ref, gcase%tag, error)
      if (allocated(error)) return

      call prepare_lsf(lsf_scr, prior, error)
      if (allocated(error)) return
      if (lsf_scr%active_count() == 0) then
         call test_failed(error, trim(gcase%tag)//" point "//itoa(ip)// &
                          ": prior point leaves no active atom to go stale")
         return
      end if
      call prepare_lsf(lsf_scr, point, error)
      if (allocated(error)) return
      call prepare_lsf(lsf_ref, point, error)
      if (allocated(error)) return

      if (stream%fade_enabled) then
         call assert_same_reach(lsf_scr%active_count(), lsf_ref%active_count(), gcase%tag, &
                                ip, error)
      else
         call assert_unscreened(lsf_ref%active_count(), mol%nat, gcase%tag, ip, error)
      end if
      if (allocated(error)) return

      n0 = stream%n
      call emit_block(lsf_scr, record_id_type(kind, gcase%tag, ip, "S"), sel, stream, f0_scr)
      n1 = stream%n
      call emit_block(lsf_ref, record_id_type(kind, gcase%tag, ip, "U"), sel, stream, f0_ref)
      call stream%emit(record_id_type(kind, gcase%tag, ip, "D"), "f0_delta", idx6(), &
                       f0_scr - f0_ref)
      if (stream%fade_enabled) call check_flags_identical(stream, n0, n1, gcase%tag, ip)
   end subroutine emit_point

   !> Record the fade of an SvdW or CFC instance for the fixture header
   !>
   !> - the first instance sets it, every later one must match: the header
   !>   holds one band for the whole fixture
   !> - CFC also `a1, a2, c, m`, one set for every case
   !>
   !> @param[inout] stream Record sink
   !> @param[in]    lsf    Instance as [[new_lsf]] built it
   !> @param[in]    tag    Case tag for the failure message
   !> @param[out]   error  testdrive failure, instances disagree
   subroutine record_fade(stream, lsf, tag, error)
      !> Record sink
      type(record_stream_type), intent(inout) :: stream
      !> Instance as built
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Case tag
      character(len=*), intent(in) :: tag
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      logical :: enabled
      real(wp) :: band(2)
      !> CFC `a1, a2, c` and `m` (zero for SvdW)
      real(wp) :: par(3)
      integer :: m

      par = 0.0_wp
      m = 0
      select type (lsf)
      type is (moist_cavity_drop_lsf_svdw_type)
         enabled = lsf%param%fade_enabled
         band = [lsf%param%fade_start, lsf%param%fade_width]
      type is (moist_cavity_drop_lsf_cfc_type)
         enabled = lsf%param%fade_enabled
         band = [lsf%param%fade_start, lsf%param%fade_width]
         par = [lsf%param%a1, lsf%param%a2, lsf%param%c]
         m = lsf%param%m
      class default
         return
      end select
      if (.not. stream%fade_set) then
         stream%fade_set = .true.
         stream%fade_enabled = enabled
         stream%fade_band = band
         stream%cfc_shape = par
         stream%cfc_m = m
      else if ((stream%fade_enabled .neqv. enabled) .or. any(stream%fade_band /= band) .or. &
               any(stream%cfc_shape /= par) .or. stream%cfc_m /= m) then
         call test_failed(error, trim(tag)//": instances disagree on the fade or shape, "// &
                          "the fixture header records one set")
      end if
   end subroutine record_fade

   !> Emit all quantities of a prepared LSF for one screening flag
   !>
   !> Shared spatial and single-nucleus blocks
   !> SvdW-only pair, normalisation, order-4, radius and contracted blocks
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
      call emit_radius_block(lsf, id, sel, act, stream)
      call emit_contracted_block(lsf, id, sel, act, stream)
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

      real(wp) :: lsf0, lsf1_r(ndim), lsf2_rr(6)
      real(wp) :: lsf3_rrr(10), lsf4_rrrr(15)
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
            call stream%emit(id, "f012_lsf2_rr", idx6(j, k), lsf2_rr(packed_index(j, k)))
         end do
      end do

      call lsf%f3_rrr(lsf3_rrr=lsf3_rrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               call stream%emit(id, "f3_rrr", idx6(j, k, l), &
                                lsf3_rrr(packed_index3(j, k, l)))
            end do
         end do
      end do

      if (id%kind /= "svdw") return

      call lsf%f4_rrrr(lsf4_rrrr)
      do j = 1, ndim
         do k = j, ndim
            do l = k, ndim
               do m = l, ndim
                  call stream%emit(id, "f4_rrrr", idx6(j, k, l, m), &
                                   lsf4_rrrr(packed_index4(j, k, l, m)))
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
      real(wp), allocatable :: lsf3_rr_rA(:, :, :), lsf4_rrr_rA(:, :, :)
      integer :: nac, j, k, l, s, ia, atomA, slotA

      ! Nuclear outputs use active indices and caller-sized buffers
      ! Resolve returned index space from extent (see `nuc_slot`)
      nac = lsf%active_count()
      allocate (lsf1_rA(ndim, nac), lsf2_r_rA(ndim, ndim, nac))
      allocate (lsf3_rr_rA(6, ndim, nac))

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
         slotA = nuc_slot(size(lsf3_rr_rA, 3), atomA, act, "f3_rr_rA", stream)
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  call stream%emit(id, "f3_rr_rA", idx6(j, k, s, atomA), &
                                   slot_read(lsf3_rr_rA(packed_index(j, k), s, :), slotA))
               end do
            end do
         end do
      end do

      if (id%kind /= "svdw") return

      allocate (lsf4_rrr_rA(10, ndim, nac))
      call lsf%f4_rrr_rA(lsf4_rrr_rA)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(lsf4_rrr_rA, 3), atomA, act, "f4_rrr_rA", stream)
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  do l = k, ndim
                     call stream%emit(id, "f4_rrr_rA", idx6(j, k, l, s, atomA), &
                                      slot_read(lsf4_rrr_rA(packed_index3(j, k, l), s, :), slotA))
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
      real(wp), allocatable :: lsf4_rr_rArB(:, :, :, :, :)
      integer :: nac, npair, j, k, s, t, ia, ib, atomA, atomB, slotA, slotB

      nac = lsf%active_count()
      npair = min(n_pair_atoms, size(sel))
      allocate (lsf2_rArB(ndim, nac, ndim, nac))
      allocate (lsf3_r_rArB(ndim, ndim, nac, ndim, nac))
      allocate (lsf4_rr_rArB(6, ndim, nac, ndim, nac))

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
         slotA = nuc_slot(size(lsf4_rr_rArB, 3), atomA, act, "f4_rr_rArB", stream)
         do ib = 1, npair
            atomB = sel(ib)
            slotB = nuc_slot(size(lsf4_rr_rArB, 5), atomB, act, "f4_rr_rArB", stream)
            do s = 1, ndim
               do t = 1, ndim
                  do j = 1, ndim
                     do k = j, ndim
                        call stream%emit(id, "f4_rr_rArB", &
                                         idx6(j, k, s, atomA, t, atomB), &
                                         slot_read2(lsf4_rr_rArB(packed_index(j, k), s, :, t, :), &
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

   !> Emit the radius family: `f3_rr_rad`, the radius-pair and nuclear-radius blocks
   !>
   !> - single-radius ladder for the selected atoms, two-index blocks over the
   !>   first [[n_pair_atoms]] of them, symmetric spatial pairs once (j <= k)
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[in]    sel     Selected atom indices
   !> @param[in]    act     Atom -> active index map
   !> @param[inout] stream  Record sink
   subroutine emit_radius_block(lsf, id, sel, act, stream)
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

      real(wp), allocatable :: r1(:), r2(:, :), r3(:, :)
      real(wp), allocatable :: rr2(:, :), rr3(:, :, :), rr4(:, :, :)
      real(wp), allocatable :: xr2(:, :, :), xr3(:, :, :, :), xr4(:, :, :, :)
      integer :: nac, npair, j, k, s, ia, ib, atomA, atomB, slotA, slotB

      nac = lsf%active_count()
      npair = min(n_pair_atoms, size(sel))
      allocate (r1(nac), r2(ndim, nac), r3(6, nac))
      allocate (rr2(nac, nac), rr3(ndim, nac, nac), rr4(6, nac, nac))
      allocate (xr2(ndim, nac, nac), xr3(ndim, ndim, nac, nac), xr4(6, ndim, nac, nac))

      call lsf%f3_rr_rad(r1, r2, r3)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(r1), atomA, act, "f3_rr_rad", stream)
         call stream%emit(id, "f3rrrad_lsf1_rad", idx6(atomA), slot_read(r1, slotA))
         do j = 1, ndim
            call stream%emit(id, "f3rrrad_lsf2_r_rad", idx6(j, atomA), &
                             slot_read(r2(j, :), slotA))
         end do
         do j = 1, ndim
            do k = j, ndim
               call stream%emit(id, "f3_rr_rad", idx6(j, k, atomA), &
                                slot_read(r3(packed_index(j, k), :), slotA))
            end do
         end do
      end do

      call lsf%f2_radrad(rr2)
      call lsf%f3_r_radrad(rr3)
      call lsf%f4_rr_radrad(rr4)
      call lsf%f2_rA_rad(xr2)
      call lsf%f3_r_rA_rad(xr3)
      call lsf%f4_rr_rA_rad(xr4)
      do ia = 1, npair
         atomA = sel(ia)
         slotA = nuc_slot(size(rr2, 1), atomA, act, "f2_radrad", stream)
         do ib = 1, npair
            atomB = sel(ib)
            slotB = nuc_slot(size(rr2, 2), atomB, act, "f2_radrad", stream)
            call stream%emit(id, "f2_radrad", idx6(atomA, atomB), slot_read2(rr2, slotA, slotB))
            do j = 1, ndim
               call stream%emit(id, "f3_r_radrad", idx6(j, atomA, atomB), &
                                slot_read2(rr3(j, :, :), slotA, slotB))
            end do
            do j = 1, ndim
               do k = j, ndim
                  call stream%emit(id, "f4_rr_radrad", idx6(j, k, atomA, atomB), &
                                   slot_read2(rr4(packed_index(j, k), :, :), slotA, slotB))
               end do
            end do
            do s = 1, ndim
               call stream%emit(id, "f2_rA_rad", idx6(s, atomA, atomB), &
                                slot_read2(xr2(s, :, :), slotA, slotB))
            end do
            do s = 1, ndim
               do j = 1, ndim
                  call stream%emit(id, "f3_r_rA_rad", idx6(j, s, atomA, atomB), &
                                   slot_read2(xr3(j, s, :, :), slotA, slotB))
               end do
            end do
            do s = 1, ndim
               do j = 1, ndim
                  do k = j, ndim
                     call stream%emit(id, "f4_rr_rA_rad", idx6(j, k, s, atomA, atomB), &
                                      slot_read2(xr4(packed_index(j, k), s, :, :), slotA, slotB))
                  end do
               end do
            end do
         end do
      end do
   end subroutine emit_radius_block

   !> Emit the contracted families: joint `hvp_f*_rA`, `hvp_f*_rad`, `vjp_f1_rA`, `vjp_f1_rad`
   !>
   !> - directions and adjoint weights from [[contraction_vectors]] at this point
   !> - records for the selected atoms, symmetric spatial pairs once
   !>
   !> @param[in]    lsf     Prepared LSF
   !> @param[in]    id      Concrete, case, point and screening flag
   !> @param[in]    sel     Selected atom indices
   !> @param[in]    act     Atom -> active index map
   !> @param[inout] stream  Record sink
   subroutine emit_contracted_block(lsf, id, sel, act, stream)
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

      real(wp), allocatable :: v(:, :), vrad(:)
      real(wp), allocatable :: h1(:, :), h2(:, :, :), h3(:, :, :)
      real(wp), allocatable :: g1(:), g2(:, :), g3(:, :), ja(:, :), jr(:)
      real(wp) :: w0, w1(ndim), w2(ndim, ndim), w2p(6)
      integer :: nac, j, k, s, ia, atomA, slotA

      nac = lsf%active_count()
      allocate (v(ndim, lsf%ncenters), vrad(lsf%ncenters))
      call contraction_vectors(id%ip, v, vrad, w0, w1, w2)
      call pack_sym_matrix(w2, w2p, ndim)
      allocate (h1(ndim, nac), h2(ndim, ndim, nac), h3(6, ndim, nac))
      allocate (g1(nac), g2(ndim, nac), g3(6, nac), ja(ndim, nac), jr(nac))

      call lsf%hvp_f1_rA(v, h1, vrad)
      call lsf%hvp_f2_r_rA(v, h2, vrad)
      call lsf%hvp_f3_rr_rA(v, h3, vrad)
      call lsf%hvp_f1_rad(v, vrad, g1)
      call lsf%hvp_f2_r_rad(v, vrad, g2)
      call lsf%hvp_f3_rr_rad(v, vrad, g3)
      call lsf%vjp_f1_rA(w0, w1, w2p, ja)
      call lsf%vjp_f1_rad(w0, w1, w2p, jr)
      do ia = 1, size(sel)
         atomA = sel(ia)
         slotA = nuc_slot(size(h1, 2), atomA, act, "hvp_f1_rA", stream)
         do s = 1, ndim
            call stream%emit(id, "hvp_f1_rA", idx6(s, atomA), slot_read(h1(s, :), slotA))
         end do
         do s = 1, ndim
            do j = 1, ndim
               call stream%emit(id, "hvp_f2_r_rA", idx6(j, s, atomA), &
                                slot_read(h2(j, s, :), slotA))
            end do
         end do
         do s = 1, ndim
            do j = 1, ndim
               do k = j, ndim
                  call stream%emit(id, "hvp_f3_rr_rA", idx6(j, k, s, atomA), &
                                   slot_read(h3(packed_index(j, k), s, :), slotA))
               end do
            end do
         end do
         call stream%emit(id, "hvp_f1_rad", idx6(atomA), slot_read(g1, slotA))
         do j = 1, ndim
            call stream%emit(id, "hvp_f2_r_rad", idx6(j, atomA), slot_read(g2(j, :), slotA))
         end do
         do j = 1, ndim
            do k = j, ndim
               call stream%emit(id, "hvp_f3_rr_rad", idx6(j, k, atomA), &
                                slot_read(g3(packed_index(j, k), :), slotA))
            end do
         end do
         do s = 1, ndim
            call stream%emit(id, "vjp_f1_rA", idx6(s, atomA), slot_read(ja(s, :), slotA))
         end do
         call stream%emit(id, "vjp_f1_rad", idx6(atomA), slot_read(jr, slotA))
      end do
   end subroutine emit_contracted_block

   !> Deterministic joint direction and adjoint weights of one evaluation point
   !>
   !> - `v(i, A) = sin(0.7 i + 1.3 A + 0.37 ip)`, user-space atom `A`
   !> - `vrad(A) = 0.5 cos(1.1 A + 0.29 ip)`
   !> - `w0 = 0.6 + 0.05 ip`, `w1(i) = cos(0.9 i + 0.31 ip)`,
   !>   `w2(i, j) = (sin(0.4 i + 0.8 j + 0.2 ip) + sin(0.4 j + 0.8 i + 0.2 ip))/2`
   !>   (symmetric)
   !>
   !> @param[in]  ip   Evaluation-point index
   !> @param[out] v    Nuclear direction [ndim, ncenters]
   !> @param[out] vrad Radius direction [ncenters]
   !> @param[out] w0   Adjoint weight of the value
   !> @param[out] w1   Adjoint weights of the gradient
   !> @param[out] w2   Adjoint weights of the Hessian
   pure subroutine contraction_vectors(ip, v, vrad, w0, w1, w2)
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> Nuclear direction
      real(wp), intent(out) :: v(:, :)
      !> Radius direction
      real(wp), intent(out) :: vrad(:)
      !> Adjoint weight of the value
      real(wp), intent(out) :: w0
      !> Adjoint weights of the gradient
      real(wp), intent(out) :: w1(ndim)
      !> Adjoint weights of the Hessian
      real(wp), intent(out) :: w2(ndim, ndim)

      integer :: i, j, a
      real(wp) :: t

      t = real(ip, wp)
      do a = 1, size(v, 2)
         do i = 1, ndim
            v(i, a) = sin(0.7_wp*real(i, wp) + 1.3_wp*real(a, wp) + 0.37_wp*t)
         end do
         vrad(a) = 0.5_wp*cos(1.1_wp*real(a, wp) + 0.29_wp*t)
      end do
      w0 = 0.6_wp + 0.05_wp*t
      do i = 1, ndim
         w1(i) = cos(0.9_wp*real(i, wp) + 0.31_wp*t)
         do j = 1, ndim
            w2(i, j) = 0.5_wp*(sin(0.4_wp*real(i, wp) + 0.8_wp*real(j, wp) + 0.2_wp*t) + &
                               sin(0.4_wp*real(j, wp) + 0.8_wp*real(i, wp) + 0.2_wp*t))
         end do
      end do
   end subroutine contraction_vectors

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

   !> Require the same active count for both flags (fade on)
   !>
   !> With the fade on the reach is the band end whatever the threshold, so the
   !> production and the threshold-0 LSF must screen alike
   !>
   !> @param[in]  nscr   Active count at the production threshold
   !> @param[in]  nref   Active count at threshold 0
   !> @param[in]  tag    Case tag, for the message
   !> @param[in]  ip     Evaluation-point index, for the message
   !> @param[out] error  testdrive failure
   subroutine assert_same_reach(nscr, nref, tag, ip, error)
      !> Active count at the production threshold
      integer, intent(in) :: nscr
      !> Active count at threshold 0
      integer, intent(in) :: nref
      !> Case tag
      character(len=*), intent(in) :: tag
      !> Evaluation-point index
      integer, intent(in) :: ip
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      if (nscr == nref) return
      call test_failed(error, "fade on: the screening threshold changes the active set for "// &
                       trim(tag)//" point "//itoa(ip)//": "//itoa(nscr)//" vs "//itoa(nref))
   end subroutine assert_same_reach

   !> Flag a problem unless the `S` and `U` blocks hold bit-identical values
   !>
   !> - `S` block: records `n0+1..n1`; `U` block: the next `n1-n0`; then the `D`
   !>   record, which must be an exact zero
   !>
   !> @param[inout] stream  Record sink
   !> @param[in]    n0      Stream length before the `S` block
   !> @param[in]    n1      Stream length before the `U` block
   !> @param[in]    tag     Case tag, for the message
   !> @param[in]    ip      Evaluation-point index, for the message
   subroutine check_flags_identical(stream, n0, n1, tag, ip)
      !> Record sink
      type(record_stream_type), intent(inout) :: stream
      !> Stream length before the `S` block
      integer, intent(in) :: n0
      !> Stream length before the `U` block
      integer, intent(in) :: n1
      !> Case tag
      character(len=*), intent(in) :: tag
      !> Evaluation-point index
      integer, intent(in) :: ip

      integer :: i, m

      m = n1 - n0
      if (stream%n /= n1 + m + 1) then
         call stream%flag_problem("S and U block lengths differ for "//trim(tag)// &
                                  " point "//itoa(ip))
         return
      end if
      do i = 1, m
         if (stream%rec(n0 + i)%val == stream%rec(n1 + i)%val) cycle
         call stream%flag_problem("fade on: S and U records differ for "//trim(tag)// &
                                  " point "//itoa(ip)//" ("// &
                                  trim(stream%rec(n0 + i)%quantity)//")")
         return
      end do
      if (stream%rec(stream%n)%val /= 0.0_wp) then
         call stream%flag_problem("fade on: f0_delta not exactly zero for "// &
                                  trim(tag)//" point "//itoa(ip))
      end if
   end subroutine check_flags_identical

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
   !> @param[in]  tol    Tolerance in rel_deviation
   !> @param[out] error  testdrive failure
   subroutine compare_stream(got, ref, path, tol, error)
      !> Emitted stream
      type(record_stream_type), intent(in) :: got
      !> Parsed fixture
      type(record_type), intent(in) :: ref(:)
      !> Fixture path
      character(len=*), intent(in) :: path
      !> Tolerance in rel_deviation
      real(wp), intent(in) :: tol
      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      integer :: i, nfail
      character(len=:), allocatable :: first
      character(len=80) :: values

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
         if (rel_deviation(got%rec(i)%val, ref(i)%val) <= tol) cycle
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

end module test_cavity_drop_lsf_reference
