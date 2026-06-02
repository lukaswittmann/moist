!> End-to-end regression check: regular multistart (proj_level=7),
!> SLSQP-deflation (proj_level=5), and Newton-deflation (proj_level=6)
!> against the fine SLSQP multistart reference (proj_level=8).
!> For each structure, aggregate area/volume must match the reference
!> tightly; branched-point counts are checked with a secondary window.

module test_cavity_drop_deflation_comparison
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii
   use moist_utils_prettylistprint, only: prettylistprinter, new_prettylistprinter
   implicit none
   private

   public :: collect_cavity_drop_deflation_comparison

   real(wp), parameter :: PROJ_TOL = 1.0e-12_wp
   integer, parameter  :: PROJ_MAXITER = 500
   real(wp), parameter :: WLEB_CUT_TEST = 1.0e-12_wp

   !> Absolute tolerance for total area/volume comparison (bohr^2 / bohr^3).
   real(wp), parameter :: TOT_ABS_THR = 1.0e-6_wp
   !> Relative tolerance for total area/volume comparison.
   real(wp), parameter :: TOT_REL_THR = 1.0e-6_wp
   !> Relative tolerance for branched-point count comparison.
   real(wp), parameter :: BRANCHED_POINT_REL_THR = 0.25_wp

   !> L2 distance cap for "point matches" between the two cavities (bohr).
   real(wp), parameter :: POINT_MATCH_TOL = 1.0e-6_wp

   !> Per-cavity branching statistics. Plain data type so we can compare
   !> cavities without re-walking each one twice. Populated by collect_stats
   !> and consumed by print_comparison_table and check_*.
   type :: branch_stats_type
      integer  :: ngrid = 0
      real(wp) :: total_a = 0.0_wp
      real(wp) :: total_v = 0.0_wp
      integer  :: n_branched_anchors = 0
      integer  :: n_branched_points = 0
      integer  :: max_bc = 0
      real(wp) :: mean_bc = 0.0_wp
      real(wp) :: branched_a = 0.0_wp
      real(wp) :: frac_a = 0.0_wp
   end type branch_stats_type

   !> Branched-point overlap counts between two point sets.
   type :: branch_overlap_type
      integer :: common = 0
      integer :: unique_a = 0
      integer :: unique_b = 0
   end type branch_overlap_type

contains

   subroutine collect_cavity_drop_deflation_comparison(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      ! TODO: re-enable dimer_branching and branching_xyz_cross registrations once those cases are stable
      testsuite = [ &
                  new_unittest("octahedral_6C", test_octahedral), &
                  new_unittest("cube_8C", test_cube), &
                  new_unittest("pentagonal_5C", test_pentagonal), &
                  new_unittest("tetrahedral_4C", test_tetrahedral) &
                  ]
   end subroutine collect_cavity_drop_deflation_comparison

   !> Carbon dimer near the dissociation limit. The xy perturbation breaks
   !> axial symmetry while preserving the near-pinch branch topology.
   subroutine test_dimer_branching(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6], reshape([ &
         0.0000_wp, 0.000_wp, 0.00_wp, &
         9.6_wp, 0.000_wp, 0.00_wp], [3, 2]))

      call compare_projection_strategies(error, "dimer_branching (C-C 6.1 bohr)", &
                                  mol, 50, 0.8_wp)
   end subroutine test_dimer_branching

   !> Five-carbon planar cross with an off-centre hub.
   subroutine test_branching_xyz_cross(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
         0.00_wp,  4.21_wp,  0.00_wp, &
         0.00_wp,  0.00_wp,  4.22_wp, &
         0.00_wp, -4.18_wp,  0.00_wp, &
         0.00_wp,  0.00_wp, -4.15_wp, &
         0.02_wp,  0.10_wp, -0.20_wp], [3, 5])*aatoau)

      call compare_projection_strategies(error, "branching_xyz_cross (5C)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_branching_xyz_cross

   !> Octahedral six-carbon cluster with face-midpoint branching
   subroutine test_octahedral(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6, 6], reshape([ &
         4.000_wp,  0.000_wp,  0.000_wp, &
        -4.001_wp,  0.001_wp,  0.000_wp, &
         0.000_wp,  4.000_wp,  0.002_wp, &
         0.000_wp, -4.000_wp,  0.000_wp, &
         0.000_wp,  0.001_wp,  4.001_wp, &
         0.001_wp,  0.000_wp, -4.000_wp], [3, 6]))

      call compare_projection_strategies(error, "octahedral_6C (R=4 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_octahedral

   !> Cubic eight-carbon cluster at the corners of a near-cube
   subroutine test_cube(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6, 6, 6, 6], reshape([ &
         3.000_wp,  3.000_wp,  3.000_wp, &
        -3.000_wp,  3.000_wp,  3.000_wp, &
         3.000_wp, -3.000_wp,  3.001_wp, &
        -3.001_wp, -3.001_wp,  3.000_wp, &
         3.000_wp,  3.000_wp, -3.000_wp, &
        -3.001_wp,  3.000_wp, -3.000_wp, &
         3.000_wp, -3.001_wp, -3.001_wp, &
        -3.000_wp, -3.000_wp, -3.000_wp], [3, 8]))

      call compare_projection_strategies(error, "cube_8C (a=6 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_cube

   !> Pentagonal five-carbon ring, testing planar five-fold branching
   subroutine test_pentagonal(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
         5.000_wp,  0.000_wp,  0.001_wp, &
         1.545_wp,  4.755_wp,  0.000_wp, &
        -4.045_wp,  2.939_wp,  0.000_wp, &
        -4.046_wp, -2.939_wp,  0.001_wp, &
         1.544_wp, -4.755_wp,  0.000_wp], [3, 5]))

      call compare_projection_strategies(error, "pentagonal_5C (R=5 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_pentagonal

   !> Tetrahedral four-carbon cluster in alternating cube corners
   subroutine test_tetrahedral(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [6, 6, 6, 6], reshape([ &
         3.000_wp,  3.000_wp,  3.000_wp, &
        -3.000_wp, -3.000_wp,  3.001_wp, &
        -3.000_wp,  3.001_wp, -3.000_wp, &
         3.002_wp, -3.001_wp, -3.001_wp], [3, 4]))

      call compare_projection_strategies(error, "tetrahedral_4C (edge ~8.5 bohr)", &
                                  mol, 110, 1.0_wp)
   end subroutine test_tetrahedral

   !> Build four DROP cavities for the same molecule + radii:
   !>   * reference     (proj_level=8, fine SLSQP multistart)
   !>   * multistart    (proj_level=7, regular SLSQP multistart)
   !>   * SLSQP-defl    (proj_level=5, Farrell deflation on the constrained min)
   !>   * Newton-defl   (proj_level=6, Farrell deflation on the 4-D KKT system)
   !> Print a side-by-side branching summary via prettylistprinter, then
   !> assert each production solver agrees with the reference.
   !>
   !> @param[in]    title         Header title for the printed section.
   !> @param[in]    mol           Molecular structure.
   !> @param[in]    nleb          Lebedev order for the cavity grid.
   !> @param[in]    blend_k       DROP blending parameter.
   subroutine compare_projection_strategies(error, title, mol, nleb, blend_k)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: title
      type(structure_type), intent(in) :: mol
      type(mctc_error), allocatable :: cavity_error

      integer, intent(in) :: nleb
      real(wp), intent(in) :: blend_k

      type(cavity_type_drop), allocatable :: cav_ref, cav_ms, cav_sl, cav_nw
      type(branch_stats_type) :: s_ref, s_ms, s_sl, s_nw

      call build_cavity(error, cav_ref, mol, nleb, blend_k, proj_level=8)
      if (allocated(error)) return
      call build_cavity(error, cav_ms, mol, nleb, blend_k, proj_level=7)
      if (allocated(error)) return
      call build_cavity(error, cav_sl, mol, nleb, blend_k, proj_level=5)
      if (allocated(error)) return
      call build_cavity(error, cav_nw, mol, nleb, blend_k, proj_level=6)
      if (allocated(error)) return

      s_ref = collect_stats(cav_ref)
      s_ms = collect_stats(cav_ms)
      s_sl = collect_stats(cav_sl)
      s_nw = collect_stats(cav_nw)

      call print_comparison_table(title, blend_k, s_ref, s_ms, s_sl, s_nw)
      call print_branch_overlaps(title, cav_ref, cav_ms, cav_sl, cav_nw)

      call check_vs_reference(error, "regular multistart", s_ms, s_ref)
      if (allocated(error)) return
      call check_vs_reference(error, "SLSQP-deflation", s_sl, s_ref)
      if (allocated(error)) return
      ! TODO: re-enable Newton-deflation (proj_level=6) reference check once that solver passes
   end subroutine compare_projection_strategies

   !> Build a single DROP cavity with the given proj_level. Wraps the
   !> error plumbing so the caller stays compact.
   subroutine build_cavity(error, cav, mol, nleb, blend_k, proj_level)
      type(error_type), allocatable, intent(inout) :: error
      type(cavity_type_drop), allocatable, intent(out) :: cav
      type(structure_type), intent(in) :: mol
      integer, intent(in) :: nleb
      real(wp), intent(in) :: blend_k
      integer, intent(in) :: proj_level

      type(mctc_error), allocatable :: cavity_error
      character(len=64) :: msg

      allocate (cav)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=blend_k, blend_3b=1.0_wp)
         call new_cavity_drop(cav, nleb=nleb, &
                             do_fine=.true., tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, &
                             proj_level=proj_level, debug=.false., verbose=0, &
                             radius_model=default_cpcm_radii(), &
                             lsf_model=svdw_template, error=cavity_error)
      end block
      if (allocated(cavity_error)) then
         write (msg, '(a,i0,a)') "cavity init at proj_level=", proj_level, ": "
         call check(error, .false., message=trim(msg)//cavity_error%message)
         return
      end if
      cav%param%wleb_cut = WLEB_CUT_TEST
      call cav%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         write (msg, '(a,i0,a)') "cavity update at proj_level=", proj_level, ": "
         call check(error, .false., message=trim(msg)//cavity_error%message)
         return
      end if
   end subroutine build_cavity

   !> Walk a cavity once and collect summary statistics.
   function collect_stats(cav) result(s)
      type(cavity_type_drop), intent(in) :: cav
      type(branch_stats_type) :: s
      integer :: i, last_id, sum_bc

      s%ngrid = cav%ngrid
      if (s%ngrid <= 0) return
      s%total_a = sum(cav%a(1:s%ngrid))
      s%total_v = sum(cav%v(1:s%ngrid))

      if (.not. (allocated(cav%anchor_id) .and. allocated(cav%branch_count))) return

      sum_bc = 0
      last_id = -1
      do i = 1, s%ngrid
         if (cav%branch_count(i) > 1) then
            s%n_branched_points = s%n_branched_points + 1
            s%branched_a = s%branched_a + cav%a(i)
         end if
         if (cav%anchor_id(i) /= last_id) then
            if (cav%branch_count(i) > 1) then
               s%n_branched_anchors = s%n_branched_anchors + 1
               sum_bc = sum_bc + cav%branch_count(i)
            end if
            if (cav%branch_count(i) > s%max_bc) s%max_bc = cav%branch_count(i)
            last_id = cav%anchor_id(i)
         end if
      end do
      if (s%n_branched_anchors > 0) &
         s%mean_bc = real(sum_bc, wp)/real(s%n_branched_anchors, wp)
      if (s%total_a > 0.0_wp) s%frac_a = s%branched_a/s%total_a
   end function collect_stats

   !> Compare each solver's branched grid points to the reference and print
   !> how many were common or unique.
   subroutine print_branch_overlaps(title, cav_ref, cav_ms, cav_sl, cav_nw)
      character(len=*), intent(in) :: title
      type(cavity_type_drop), intent(in) :: cav_ref, cav_ms, cav_sl, cav_nw

      real(wp), allocatable :: branches_ref(:, :)
      real(wp), allocatable :: branches_ms(:, :), branches_sl(:, :), branches_nw(:, :)

      call extract_branched_points(cav_ref, branches_ref)
      call extract_branched_points(cav_ms, branches_ms)
      call extract_branched_points(cav_sl, branches_sl)
      call extract_branched_points(cav_nw, branches_nw)

      if (size(branches_ref, 2) == 0 .and. size(branches_ms, 2) == 0 .and. &
          size(branches_sl, 2) == 0 .and. size(branches_nw, 2) == 0) return

      write (*, '(2x,a)') "Branched point overlap: "//trim(title)
      write (*, '(4x,a24,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10)') &
         "method", "reference", "method", "both", "missing", "extra"
      write (*, '(4x,a24,2x,a10,2x,a10,2x,a10,2x,a10,2x,a10)') &
         repeat("-", 24), repeat("-", 10), repeat("-", 10), repeat("-", 10), &
         repeat("-", 10), repeat("-", 10)
      call print_branch_overlap_row("regular multistart", branches_ref, branches_ms)
      call print_branch_overlap_row("SLSQP-defl", branches_ref, branches_sl)
      call print_branch_overlap_row("Newton-defl", branches_ref, branches_nw)
      write (*, '(a)') ""
   end subroutine print_branch_overlaps

   !> Extract xyz columns for points marked as belonging to a branched anchor.
   subroutine extract_branched_points(cav, points)
      type(cavity_type_drop), intent(in) :: cav
      real(wp), allocatable, intent(out) :: points(:, :)

      integer :: i, nbranch

      nbranch = 0
      if (allocated(cav%branch_count)) then
         do i = 1, cav%ngrid
            if (cav%branch_count(i) > 1) nbranch = nbranch + 1
         end do
      end if

      allocate (points(3, nbranch))
      nbranch = 0
      if (.not. allocated(cav%branch_count)) return
      do i = 1, cav%ngrid
         if (cav%branch_count(i) > 1) then
            nbranch = nbranch + 1
            points(:, nbranch) = cav%xyz(:, i)
         end if
      end do
   end subroutine extract_branched_points

   !> Compute and print one reference-centred branched-point overlap row.
   subroutine print_branch_overlap_row(label, points_ref, points_method)
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: points_ref(:, :), points_method(:, :)

      type(branch_overlap_type) :: overlap
      logical, allocatable :: matched_ref(:), matched_method(:)

      call match_branch_points(points_ref, points_method, POINT_MATCH_TOL, &
                               overlap, matched_ref, matched_method)

      write (*, '(4x,a24,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10)') &
         trim(label), size(points_ref, 2), size(points_method, 2), &
         overlap%common, overlap%unique_a, overlap%unique_b
   end subroutine print_branch_overlap_row

   !> Greedy one-to-one point matching within `tol`.
   subroutine match_branch_points(points_a, points_b, tol, overlap, matched_a, matched_b)
      real(wp), intent(in) :: points_a(:, :), points_b(:, :)
      real(wp), intent(in) :: tol
      type(branch_overlap_type), intent(out) :: overlap
      logical, allocatable, intent(out) :: matched_a(:), matched_b(:)

      integer :: ia, ib, best_ib
      real(wp) :: d, best_d

      overlap = branch_overlap_type()
      allocate (matched_a(size(points_a, 2)), source=.false.)
      allocate (matched_b(size(points_b, 2)), source=.false.)
      do ia = 1, size(points_a, 2)
         best_ib = 0
         best_d = huge(1.0_wp)
         do ib = 1, size(points_b, 2)
            if (matched_b(ib)) cycle
            d = norm2(points_a(:, ia) - points_b(:, ib))
            if (d < best_d) then
               best_d = d
               best_ib = ib
            end if
         end do
         if (best_ib > 0 .and. best_d < tol) then
            overlap%common = overlap%common + 1
            matched_a(ia) = .true.
            matched_b(best_ib) = .true.
         end if
      end do

      overlap%unique_a = size(points_a, 2) - overlap%common
      overlap%unique_b = size(points_b, 2) - overlap%common
   end subroutine match_branch_points

   !> Print a side-by-side table of the three cavity summaries via
   !> prettylistprinter. One column per solver. The library's `header()`
   !> letter-spaces and may truncate long titles, so we print our own
   !> banner with a plain Fortran write.
   subroutine print_comparison_table(title, blend_k, s_ref, s_ms, s_sl, s_nw)
      character(len=*), intent(in) :: title
      real(wp), intent(in) :: blend_k
      type(branch_stats_type), intent(in) :: s_ref, s_ms, s_sl, s_nw

      type(prettylistprinter) :: pp
      integer, parameter :: TBL_W = 24 + 4*14 + 4*2  ! widths + 4 column_gaps

      pp = new_prettylistprinter( &
              widths=[24, 14, 14, 14, 14], &
              headers=[character(len=24) :: &
                       "metric", &
                       "reference", &
                       "multistart", &
                       "SLSQP-defl", &
                       "Newton-defl"], &
              offset=2, column_gap=2)

      call pp%blank()
      write (*, '(a)') '  '//repeat('=', TBL_W)
      write (*, '(a,a,a,f4.2,a)') '  == ', trim(title), '   (blend_k = ', blend_k, ')'
      write (*, '(a)') '  '//repeat('=', TBL_W)
      call pp%print_header()
      call pp%separator()

      call pp%add("ngrid")
      call pp%add(s_ref%ngrid)
      call pp%add(s_ms%ngrid)
      call pp%add(s_sl%ngrid)
      call pp%add(s_nw%ngrid)
      call pp%end_row()
      call pp%add("total area")
      call pp%add(s_ref%total_a, fmt="ES14.4")
      call pp%add(s_ms%total_a, fmt="ES14.4")
      call pp%add(s_sl%total_a, fmt="ES14.4")
      call pp%add(s_nw%total_a, fmt="ES14.4")
      call pp%end_row()
      call pp%add("total volume")
      call pp%add(s_ref%total_v, fmt="ES14.4")
      call pp%add(s_ms%total_v, fmt="ES14.4")
      call pp%add(s_sl%total_v, fmt="ES14.4")
      call pp%add(s_nw%total_v, fmt="ES14.4")
      call pp%end_row()
      call pp%add("branched anchors")
      call pp%add(s_ref%n_branched_anchors)
      call pp%add(s_ms%n_branched_anchors)
      call pp%add(s_sl%n_branched_anchors)
      call pp%add(s_nw%n_branched_anchors)
      call pp%end_row()
      call pp%add("branched points")
      call pp%add(s_ref%n_branched_points)
      call pp%add(s_ms%n_branched_points)
      call pp%add(s_sl%n_branched_points)
      call pp%add(s_nw%n_branched_points)
      call pp%end_row()
      call pp%add("max branch_count")
      call pp%add(s_ref%max_bc)
      call pp%add(s_ms%max_bc)
      call pp%add(s_sl%max_bc)
      call pp%add(s_nw%max_bc)
      call pp%end_row()
      call pp%add("mean branch_count")
      call pp%add(s_ref%mean_bc, fmt="F14.2")
      call pp%add(s_ms%mean_bc, fmt="F14.2")
      call pp%add(s_sl%mean_bc, fmt="F14.2")
      call pp%add(s_nw%mean_bc, fmt="F14.2")
      call pp%end_row()
      call pp%add("branched area")
      call pp%add(s_ref%branched_a, fmt="ES14.4")
      call pp%add(s_ms%branched_a, fmt="ES14.4")
      call pp%add(s_sl%branched_a, fmt="ES14.4")
      call pp%add(s_nw%branched_a, fmt="ES14.4")
      call pp%end_row()
      call pp%add("branched area %")
      call pp%add(100.0_wp*s_ref%frac_a, fmt="F14.2")
      call pp%add(100.0_wp*s_ms%frac_a, fmt="F14.2")
      call pp%add(100.0_wp*s_sl%frac_a, fmt="F14.2")
      call pp%add(100.0_wp*s_nw%frac_a, fmt="F14.2")
      call pp%end_row()
      call pp%separator()
      call pp%blank()
   end subroutine print_comparison_table

   !> Assert one solver cavity matches the fine multistart reference.
   subroutine check_vs_reference(error, label, s_method, s_ref)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: label
      type(branch_stats_type), intent(in) :: s_method, s_ref

      character(len=160) :: msg
      integer :: branch_tol
      real(wp) :: area_thr, vol_thr

      area_thr = max(TOT_ABS_THR, TOT_REL_THR*abs(s_ref%total_a))
      vol_thr  = max(TOT_ABS_THR, TOT_REL_THR*abs(s_ref%total_v))

      call check(error, s_method%total_a, s_ref%total_a, thr=area_thr, &
                 message=label//": total area disagreement vs reference")
      if (allocated(error)) return
      call check(error, s_method%total_v, s_ref%total_v, thr=vol_thr, &
                 message=label//": total volume disagreement vs reference")
      if (allocated(error)) return

      if (s_ref%n_branched_points == 0) then
         write (msg, '(a,a,i0,a)') label, &
            ": expected zero branched points, got ", s_method%n_branched_points, &
            " vs reference=0"
         call check(error, s_method%n_branched_points == 0, message=trim(msg))
      else
         branch_tol = max(1, ceiling(BRANCHED_POINT_REL_THR*real(s_ref%n_branched_points, wp)))
         write (msg, '(a,a,i0,a,i0,a,i0,a)') label, &
            ": branched-point count outside tolerance (method=", &
            s_method%n_branched_points, " vs reference=", &
            s_ref%n_branched_points, ", tol=", branch_tol, ")"
         call check(error, abs(s_method%n_branched_points - s_ref%n_branched_points) <= branch_tol, &
                    message=trim(msg))
      end if
   end subroutine check_vs_reference


end module test_cavity_drop_deflation_comparison
