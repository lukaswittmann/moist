!> Test suite for the per-cell atom screening grid in moist_math_cell_grid.
module test_math_cell_grid
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_cell_grid, only: moist_cell_grid_type
   implicit none (type, external)
   private

   public :: collect_math_cell_grid

contains

   !> Collect all cell-grid tests.
   subroutine collect_math_cell_grid(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("query_covers_bruteforce", test_query_covers_bruteforce), &
         new_unittest("single_atom", test_single_atom), &
         new_unittest("coincident_atoms", test_coincident_atoms), &
         new_unittest("query_outside_bbox", test_query_outside_bbox), &
         new_unittest("full_scan_below_threshold", test_full_scan_below_threshold), &
         new_unittest("full_scan_threshold_not_triggered", test_full_scan_not_triggered), &
         new_unittest("cell_fraction_superset", test_cell_fraction_superset), &
         new_unittest("cell_fraction_finer_cells", test_cell_fraction_finer_cells), &
         new_unittest("cell_fraction_default_identity", test_cell_fraction_default_identity) &
      ]
   end subroutine collect_math_cell_grid

   !> For a heterogeneous cluster of atoms and a lattice of query points
   !> spanning the bounding box, assert the cell-grid list is a superset of the
   !> brute-force candidate set for every point.
   subroutine test_query_covers_bruteforce(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 6)
      real(wp) :: r_eff(6)
      real(wp) :: point(3)
      real(wp) :: xmin(3), xmax(3)
      real(wp) :: tol
      integer :: i, j, start, n, nx_sweep, ix, iy, iz
      logical :: in_bruteforce, in_grid

      ! Mixed cluster: heterogeneous radii, deliberate clustering
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.2_wp, 0.1_wp, -0.2_wp]
      xyz(:, 3) = [2.6_wp, 0.0_wp, 0.0_wp]
      xyz(:, 4) = [0.3_wp, 2.1_wp, 0.0_wp]
      xyz(:, 5) = [-0.5_wp, 0.0_wp, 2.3_wp]
      xyz(:, 6) = [1.0_wp, 1.0_wp, 1.0_wp]

      r_eff = [1.5_wp, 1.2_wp, 1.4_wp, 1.0_wp, 1.3_wp, 1.1_wp]

      call grid%build(xyz, r_eff)

      ! Sweep a lattice of query points across bbox +- max r_eff
      xmin = minval(xyz, dim=2) - maxval(r_eff)
      xmax = maxval(xyz, dim=2) + maxval(r_eff)
      tol = 10.0_wp * epsilon(1.0_wp) * maxval(abs([xmin, xmax]))
      nx_sweep = 7

      do ix = 0, nx_sweep - 1
      do iy = 0, nx_sweep - 1
      do iz = 0, nx_sweep - 1
         point(1) = xmin(1) + (xmax(1) - xmin(1)) * real(ix, wp) / real(nx_sweep - 1, wp)
         point(2) = xmin(2) + (xmax(2) - xmin(2)) * real(iy, wp) / real(nx_sweep - 1, wp)
         point(3) = xmin(3) + (xmax(3) - xmin(3)) * real(iz, wp) / real(nx_sweep - 1, wp)

         call grid%query(point, start, n)

         ! For every atom, check that brute-force membership implies grid membership.
         ! (Points outside the bbox are strictly clamped, so this check only applies
         !  to points inside the bbox.)
         if (.not. point_inside_bbox(point, minval(xyz, dim=2), maxval(xyz, dim=2))) cycle

         do j = 1, size(xyz, 2)
            in_bruteforce = norm2(point - xyz(:, j)) <= r_eff(j) + tol

            in_grid = .false.
            do i = 1, n
               if (grid%cell_nlat(start + i) == j) then
                  in_grid = .true.
                  exit
               end if
            end do

            if (in_bruteforce .and. .not. in_grid) then
               call failure(error, "brute_force_set not subset of grid_set")
               return
            end if
         end do
      end do
      end do
      end do

      call grid%destroy()
   end subroutine test_query_covers_bruteforce

   !> Degenerate single-atom case.
   subroutine test_single_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 1)
      real(wp) :: r_eff(1)
      integer :: start, n

      xyz(:, 1) = [1.0_wp, 2.0_wp, 3.0_wp]
      r_eff = [0.7_wp]

      call grid%build(xyz, r_eff)

      call check(error, grid%ncells == 1, "expected a single cell")
      if (allocated(error)) return

      call grid%query(xyz(:, 1), start, n)
      call check(error, n == 1, "expected one candidate for a single atom")
      if (allocated(error)) return

      call check(error, grid%cell_nlat(start + 1) == 1, "candidate must be atom 1")
      if (allocated(error)) return

      call grid%destroy()
   end subroutine test_single_atom

   !> All atoms at the same position.
   subroutine test_coincident_atoms(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 4)
      real(wp) :: r_eff(4)
      integer :: start, n, j, count_j

      xyz = 0.0_wp
      r_eff = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]

      call grid%build(xyz, r_eff)

      call check(error, grid%ncells == 1, "coincident atoms should collapse to one cell")
      if (allocated(error)) return

      call grid%query([0.1_wp, -0.1_wp, 0.05_wp], start, n)
      call check(error, n == 4, "all atoms should appear in the single cell")
      if (allocated(error)) return

      do j = 1, 4
         count_j = count(grid%cell_nlat(start + 1 : start + n) == j)
         call check(error, count_j == 1, "each atom must appear exactly once")
         if (allocated(error)) return
      end do

      call grid%destroy()
   end subroutine test_coincident_atoms

   !> Points outside the bounding box must still return a valid (possibly empty)
   !> list; the strict clamp policy forbids out-of-range cell indices.
   subroutine test_query_outside_bbox(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 2)
      real(wp) :: r_eff(2)
      integer :: start, n

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.0_wp, 0.0_wp, 0.0_wp]
      r_eff = [0.5_wp, 0.5_wp]

      call grid%build(xyz, r_eff)

      ! Far above the bbox along +x - should not crash, n can be anything >= 0
      call grid%query([100.0_wp, 0.0_wp, 0.0_wp], start, n)
      call check(error, n >= 0, "query must return a non-negative candidate count")
      if (allocated(error)) return
      call check(error, start >= 0, "query must return a non-negative offset")
      if (allocated(error)) return
      call check(error, start + n <= size(grid%cell_nlat), &
         "returned slice must be within cell_nlat")
      if (allocated(error)) return

      call grid%destroy()
   end subroutine test_query_outside_bbox

   !> When natoms < full_scan_below, build must collapse to a single cell
   !> and every query must return the full atom list regardless of geometry.
   subroutine test_full_scan_below_threshold(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 6)
      real(wp) :: r_eff(6)
      integer :: start, n, j, count_j

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [5.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 5.0_wp, 0.0_wp]
      xyz(:, 4) = [0.0_wp, 0.0_wp, 5.0_wp]
      xyz(:, 5) = [5.0_wp, 5.0_wp, 0.0_wp]
      xyz(:, 6) = [5.0_wp, 0.0_wp, 5.0_wp]
      r_eff = 1.0_wp

      call grid%build(xyz, r_eff, full_scan_below=50)

      call check(error, grid%full_scan, "full_scan flag must be set")
      if (allocated(error)) return

      call check(error, grid%ncells == 1, "full-scan path must yield one cell")
      if (allocated(error)) return

      ! Query at a point far from any atom - full scan must still return all.
      call grid%query([100.0_wp, -100.0_wp, 200.0_wp], start, n)
      call check(error, n == 6, "full scan must return every atom")
      if (allocated(error)) return

      do j = 1, 6
         count_j = count(grid%cell_nlat(start + 1 : start + n) == j)
         call check(error, count_j == 1, "each atom must appear exactly once")
         if (allocated(error)) return
      end do

      call grid%destroy()
   end subroutine test_full_scan_below_threshold

   !> Above the threshold, the spatial-binning path must run (full_scan stays
   !> false) so callers still get the usual per-cell fan-out.
   subroutine test_full_scan_not_triggered(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 4)
      real(wp) :: r_eff(4)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [5.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 5.0_wp, 0.0_wp]
      xyz(:, 4) = [5.0_wp, 5.0_wp, 0.0_wp]
      r_eff = 1.0_wp

      ! natoms (4) >= full_scan_below (3) - spatial path must run.
      call grid%build(xyz, r_eff, full_scan_below=3)

      call check(error, .not. grid%full_scan, &
         "full_scan flag must remain false at/above threshold")
      if (allocated(error)) return

      call grid%destroy()
   end subroutine test_full_scan_not_triggered

   !> With cell_fraction < 1, the grid candidate list must still be a superset
   !> of the brute-force set for every query point inside the bounding box.
   subroutine test_cell_fraction_superset(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid
      real(wp) :: xyz(3, 6), r_eff(6), point(3)
      real(wp) :: xmin(3), xmax(3), tol
      integer :: i, j, start, n, nx_sweep, ix, iy, iz
      logical :: in_bruteforce, in_grid
      real(wp) :: fractions(2)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.2_wp, 0.1_wp, -0.2_wp]
      xyz(:, 3) = [2.6_wp, 0.0_wp, 0.0_wp]
      xyz(:, 4) = [0.3_wp, 2.1_wp, 0.0_wp]
      xyz(:, 5) = [-0.5_wp, 0.0_wp, 2.3_wp]
      xyz(:, 6) = [1.0_wp, 1.0_wp, 1.0_wp]
      r_eff = [1.5_wp, 1.2_wp, 1.4_wp, 1.0_wp, 1.3_wp, 1.1_wp]

      fractions = [0.5_wp, 0.25_wp]

      xmin = minval(xyz, dim=2) - maxval(r_eff)
      xmax = maxval(xyz, dim=2) + maxval(r_eff)
      tol = 10.0_wp * epsilon(1.0_wp) * maxval(abs([xmin, xmax]))
      nx_sweep = 7

      do i = 1, size(fractions)
         call grid%build(xyz, r_eff, cell_fraction=fractions(i))

         do ix = 0, nx_sweep - 1
         do iy = 0, nx_sweep - 1
         do iz = 0, nx_sweep - 1
            point(1) = xmin(1) + (xmax(1) - xmin(1)) * real(ix, wp) / real(nx_sweep - 1, wp)
            point(2) = xmin(2) + (xmax(2) - xmin(2)) * real(iy, wp) / real(nx_sweep - 1, wp)
            point(3) = xmin(3) + (xmax(3) - xmin(3)) * real(iz, wp) / real(nx_sweep - 1, wp)

            if (.not. point_inside_bbox(point, minval(xyz, dim=2), maxval(xyz, dim=2))) cycle

            call grid%query(point, start, n)

            do j = 1, size(xyz, 2)
               in_bruteforce = norm2(point - xyz(:, j)) <= r_eff(j) + tol
               in_grid = .false.
               block
                  integer :: k
                  do k = 1, n
                     if (grid%cell_nlat(start + k) == j) then
                        in_grid = .true.
                        exit
                     end if
                  end do
               end block

               if (in_bruteforce .and. .not. in_grid) then
                  call failure(error, "cell_fraction superset violated")
                  return
               end if
            end do
         end do
         end do
         end do

         call grid%destroy()
      end do
   end subroutine test_cell_fraction_superset

   !> Verify that cell_fraction < 1 produces more cells and fewer candidates
   !> per query than the default.
   subroutine test_cell_fraction_finer_cells(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid_default, grid_fine
      real(wp) :: xyz(3, 6), r_eff(6)
      integer :: start_d, n_d, start_f, n_f

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [3.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 3.0_wp, 0.0_wp]
      xyz(:, 4) = [0.0_wp, 0.0_wp, 3.0_wp]
      xyz(:, 5) = [3.0_wp, 3.0_wp, 0.0_wp]
      xyz(:, 6) = [3.0_wp, 0.0_wp, 3.0_wp]
      r_eff = [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp]

      call grid_default%build(xyz, r_eff)
      call grid_fine%build(xyz, r_eff, cell_fraction=0.5_wp)

      call check(error, grid_fine%ncells > grid_default%ncells, &
         "finer grid must have more cells")
      if (allocated(error)) return

      ! Query near one atom: finer grid should return fewer candidates
      call grid_default%query(xyz(:, 1), start_d, n_d)
      call grid_fine%query(xyz(:, 1), start_f, n_f)

      call check(error, n_f <= n_d, &
         "finer grid should return no more candidates than default")
      if (allocated(error)) return

      call grid_default%destroy()
      call grid_fine%destroy()
   end subroutine test_cell_fraction_finer_cells

   !> Explicit cell_fraction=1.0 must produce identical results to omitting it.
   subroutine test_cell_fraction_default_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(moist_cell_grid_type) :: grid_implicit, grid_explicit
      real(wp) :: xyz(3, 4), r_eff(4)
      integer :: start_i, n_i, start_e, n_e

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 1.0_wp, 0.0_wp]
      xyz(:, 4) = [1.0_wp, 1.0_wp, 0.0_wp]
      r_eff = [0.8_wp, 0.8_wp, 0.8_wp, 0.8_wp]

      call grid_implicit%build(xyz, r_eff)
      call grid_explicit%build(xyz, r_eff, cell_fraction=1.0_wp)

      call check(error, grid_implicit%ncells == grid_explicit%ncells, &
         "cell count must match")
      if (allocated(error)) return

      call check(error, abs(grid_implicit%cell_side - grid_explicit%cell_side) < &
         epsilon(1.0_wp), "cell side must match")
      if (allocated(error)) return

      call grid_implicit%query([0.5_wp, 0.5_wp, 0.0_wp], start_i, n_i)
      call grid_explicit%query([0.5_wp, 0.5_wp, 0.0_wp], start_e, n_e)

      call check(error, n_i == n_e, "candidate count must match")
      if (allocated(error)) return

      call grid_implicit%destroy()
      call grid_explicit%destroy()
   end subroutine test_cell_fraction_default_identity

   pure function point_inside_bbox(point, lo, hi) result(inside)
      real(wp), intent(in) :: point(3), lo(3), hi(3)
      logical :: inside
      inside = all(point >= lo) .and. all(point <= hi)
   end function point_inside_bbox

   subroutine failure(error, message)
      type(error_type), allocatable, intent(inout) :: error
      character(*), intent(in) :: message
      call check(error, .false., message)
   end subroutine failure

end module test_math_cell_grid
