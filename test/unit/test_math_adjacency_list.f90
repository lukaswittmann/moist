!> Test suite for adjacency-list utilities in moist_math_adjacency_list
module test_math_adjacency_list
   use mctc_env, only: wp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use moist_math_adjacency_list, only: adjacency_list_type
   implicit none(type, external)
   private

   public :: collect_math_adjacency_list

   real(wp), parameter :: thr = 10.0_wp*epsilon(1.0_wp)

contains

   !> Collect all adjacency-list tests
   subroutine collect_math_adjacency_list(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pair_membership_bruteforce", test_pair_membership_bruteforce), &
                  new_unittest("neighbour_content_line_cluster", test_neighbour_content), &
                  new_unittest("sorted_distances_ascending", test_sorted_distances_ascending), &
                  new_unittest("dist_matches_nlat_pairing", test_dist_matches_nlat_pairing), &
                  new_unittest("sorted_matches_unsorted_set", test_sorted_matches_unsorted_set), &
                  new_unittest("sorted_equidistant_ties", test_sorted_equidistant_ties), &
                  new_unittest("sorted_early_exit_contract", test_sorted_early_exit_contract), &
                  new_unittest("many_cells_bruteforce", test_many_cells_bruteforce), &
                  new_unittest("cutoff_exceeds_extent", test_cutoff_exceeds_extent), &
                  new_unittest("planar_and_collinear", test_planar_and_collinear), &
                  new_unittest("coincident_points", test_coincident_points), &
                  new_unittest("pair_symmetry", test_pair_symmetry), &
                  new_unittest("empty_and_single_point", test_empty_and_single_point), &
                  new_unittest("nonpositive_cutoff", test_nonpositive_cutoff), &
                  new_unittest("rebuild_grow_and_shrink", test_rebuild_grow_and_shrink), &
                  new_unittest("destroy_and_reuse", test_destroy_and_reuse), &
                  new_unittest("init_resets_sorted", test_init_resets_sorted) &
                  ]
   end subroutine collect_math_adjacency_list

   !> Compare adjacency membership against brute-force distance checks.
   subroutine test_pair_membership_bruteforce(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 6)
      real(wp), parameter :: cutoff = 1.5_wp
      real(wp) :: cutoff2, d2, tol
      integer, allocatable :: ids(:)
      integer :: i, j, expected_count
      logical :: in_expected, in_list

      ! Mix of inside, outside, and exact-cutoff pairs.
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.5_wp, 0.0_wp, 0.0_wp]    ! exactly at cutoff from 1
      xyz(:, 3) = [1.5002_wp, 0.0_wp, 0.0_wp] ! just outside cutoff from 1
      xyz(:, 4) = [0.0_wp, 1.0_wp, 0.0_wp]
      xyz(:, 5) = [0.0_wp, 0.0_wp, 2.0_wp]
      xyz(:, 6) = [0.0_wp, 1.5_wp, 0.0_wp]    ! exactly at cutoff from 1

      call nlist%init(cutoff=cutoff)
      call nlist%update(xyz)

      cutoff2 = cutoff*cutoff
      tol = 10.0_wp*epsilon(cutoff2)

      do i = 1, size(xyz, 2)
         ids = nlist%get_neighbours(i)

         expected_count = 0
         do j = 1, size(xyz, 2)
            if (j == i) cycle
            d2 = sum((xyz(:, i) - xyz(:, j))**2)
            if (d2 <= cutoff2 + tol) expected_count = expected_count + 1
         end do

         call check(error, size(ids) == expected_count, &
                    "Neighbour count mismatch against brute-force reference")
         if (allocated(error)) return

         do j = 1, size(xyz, 2)
            in_list = any(ids == j)

            if (j == i) then
               call check(error,.not. in_list, "Self index must never appear in neighbour list")
            else
               d2 = sum((xyz(:, i) - xyz(:, j))**2)
               in_expected = (d2 <= cutoff2 + tol)
               call check(error, in_list .eqv. in_expected, &
                          "Pair membership mismatch against brute-force reference")
            end if
            if (allocated(error)) return
         end do
      end do

      call nlist%destroy()
   end subroutine test_pair_membership_bruteforce

   !> Neighbour content for a simple linear cluster with one isolated point.
   subroutine test_neighbour_content(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 4)
      integer, allocatable :: ids(:)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [3.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 4) = [0.0_wp, 4.0_wp, 0.0_wp]

      call nlist%init(cutoff=2.1_wp)
      call nlist%update(xyz)

      ids = nlist%get_neighbours(1)
      call check(error, size(ids) == 1, "Point 1 should have one neighbour")
      if (allocated(error)) return
      call check(error, any(ids == 2), "Point 1 should be connected to point 2")
      if (allocated(error)) return

      ids = nlist%get_neighbours(2)
      call check(error, size(ids) == 2, "Point 2 should have two neighbours")
      if (allocated(error)) return
      call check(error, any(ids == 1) .and. any(ids == 3), &
                 "Point 2 should be connected to points 1 and 3")
      if (allocated(error)) return

      ids = nlist%get_neighbours(3)
      call check(error, size(ids) == 1, "Point 3 should have one neighbour")
      if (allocated(error)) return
      call check(error, any(ids == 2), "Point 3 should be connected to point 2")
      if (allocated(error)) return

      ids = nlist%get_neighbours(4)
      call check(error, size(ids) == 0, "Point 4 should have no neighbours")
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_neighbour_content

   !=========================================================================!
   ! Sorted path                                                             !
   !                                                                         !
   ! The iSwiG switching function walks a row and terminates with exit as    !
   ! soon as dist exceeds a per-atom break threshold, so ascending order and !
   ! the dist/nlat pairing are load-bearing for production correctness.      !
   !=========================================================================!

   !> With sorted=.true. every row must be non-decreasing in distance.
   subroutine test_sorted_distances_ascending(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 8)
      integer :: i, k

      call star_cluster(xyz)

      call nlist%init(cutoff=1.6_wp, sorted=.true.)
      call nlist%update(xyz)

      call check_csr_invariants(error, nlist, size(xyz, 2))
      if (allocated(error)) return

      ! The central point must see every other point, so at least one row is long
      ! enough to be worth sorting.
      call check(error, nlist%nnl(1) == size(xyz, 2) - 1, &
                 "Central point should neighbour all others")
      if (allocated(error)) return

      do i = 1, size(xyz, 2)
         do k = nlist%inl(i) + 2, nlist%inl(i) + nlist%nnl(i)
            call check(error, nlist%dist(k - 1) <= nlist%dist(k), &
                       "Sorted neighbour distances must be non-decreasing")
            if (allocated(error)) return
         end do
      end do

      call nlist%destroy()
   end subroutine test_sorted_distances_ascending

   !> Every stored distance must belong to the neighbour id stored beside it.
   !>
   !> Sorting permutes dist and nlat together; if the companion array ever came
   !> loose the list would stay plausible but silently mispair ids and distances.
   subroutine test_dist_matches_nlat_pairing(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 8)
      logical :: flags(2)
      real(wp) :: expected
      integer :: is, i, k

      call star_cluster(xyz)
      flags = [.false., .true.]

      do is = 1, size(flags)
         call nlist%init(cutoff=1.6_wp, sorted=flags(is))
         call nlist%update(xyz)

         call check_csr_invariants(error, nlist, size(xyz, 2))
         if (allocated(error)) return

         do i = 1, size(xyz, 2)
            do k = nlist%inl(i) + 1, nlist%inl(i) + nlist%nnl(i)
               expected = norm2(xyz(:, i) - xyz(:, nlist%nlat(k)))
               call check(error, abs(nlist%dist(k) - expected) <= thr*max(1.0_wp, expected), &
                          "Stored distance does not match the neighbour it is paired with")
               if (allocated(error)) return
            end do
         end do

         call nlist%destroy()
      end do
   end subroutine test_dist_matches_nlat_pairing

   !> Sorting must permute a row, never add or drop neighbours.
   subroutine test_sorted_matches_unsorted_set(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: plain, ordered
      real(wp) :: xyz(3, 8)
      integer :: i, k

      call star_cluster(xyz)

      call plain%init(cutoff=1.6_wp)
      call plain%update(xyz)

      call ordered%init(cutoff=1.6_wp, sorted=.true.)
      call ordered%update(xyz)

      call check(error, size(plain%nlat) == size(ordered%nlat), &
                 "Sorting must not change the total pair count")
      if (allocated(error)) return
      call check(error, all(plain%nnl == ordered%nnl), &
                 "Sorting must not change per-point neighbour counts")
      if (allocated(error)) return
      call check(error, all(plain%inl == ordered%inl), &
                 "Sorting must not change the CSR offsets")
      if (allocated(error)) return

      do i = 1, size(xyz, 2)
         do k = plain%inl(i) + 1, plain%inl(i) + plain%nnl(i)
            call check(error, any(ordered%nlat(ordered%inl(i) + 1: &
                                               ordered%inl(i) + ordered%nnl(i)) == plain%nlat(k)), &
                       "Sorted row lost a neighbour present in the unsorted row")
            if (allocated(error)) return
         end do
      end do

      call plain%destroy()
      call ordered%destroy()
   end subroutine test_sorted_matches_unsorted_set

   !> Duplicate distances must survive sorting intact.
   !>
   !> qsort is not stable and switches from insertion sort to quicksort
   !> partitioning above 24 elements, so this fixture builds a row of 40
   !> neighbours made of 20 distances that each occur twice.
   subroutine test_sorted_equidistant_ties(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 41)
      integer :: m, k, start, cnt

      ! Point 1 at the origin; mirrored pairs at +-0.1*m along x.
      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      do m = 1, 20
         xyz(:, 2*m) = [0.1_wp*real(m, wp), 0.0_wp, 0.0_wp]
         xyz(:, 2*m + 1) = [-0.1_wp*real(m, wp), 0.0_wp, 0.0_wp]
      end do

      call nlist%init(cutoff=2.05_wp, sorted=.true.)
      call nlist%update(xyz)

      call check_csr_invariants(error, nlist, size(xyz, 2))
      if (allocated(error)) return

      start = nlist%inl(1)
      cnt = nlist%nnl(1)
      call check(error, cnt == 40, "Origin must see all 40 mirrored points")
      if (allocated(error)) return

      do k = 2, cnt
         call check(error, nlist%dist(start + k - 1) <= nlist%dist(start + k), &
                    "Ties must not break the ascending order")
         if (allocated(error)) return
      end do

      ! Each distance occurs exactly twice, and no id may be lost or repeated.
      do m = 1, 20
         call check(error, count(abs(nlist%dist(start + 1:start + cnt) &
                                     - 0.1_wp*real(m, wp)) <= thr) == 2, &
                    "Each mirrored distance must appear exactly twice")
         if (allocated(error)) return
      end do

      do m = 2, size(xyz, 2)
         call check(error, count(nlist%nlat(start + 1:start + cnt) == m) == 1, &
                    "Every point must appear exactly once in the origin row")
         if (allocated(error)) return
      end do

      call check_all_rows_bruteforce(error, nlist, xyz)
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_sorted_equidistant_ties

   !> Reproduce the consumer contract: stopping at the first over-threshold
   !> distance must yield exactly the neighbours within that threshold.
   subroutine test_sorted_early_exit_contract(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 8)
      real(wp) :: thresholds(4), t
      integer, allocatable :: ref_ids(:)
      real(wp), allocatable :: ref_dist(:)
      integer :: it, i, ii, start, cnt, taken
      logical :: reached(8)

      call star_cluster(xyz)
      thresholds = [0.3_wp, 0.8_wp, 1.1_wp, 1.6_wp]

      call nlist%init(cutoff=1.6_wp, sorted=.true.)
      call nlist%update(xyz)

      do it = 1, size(thresholds)
         t = thresholds(it)

         do i = 1, size(xyz, 2)
            start = nlist%inl(i)
            cnt = nlist%nnl(i)

            ! Walk the row exactly as the switching function does.
            reached = .false.
            taken = 0
            do ii = 1, cnt
               if (nlist%dist(start + ii) > t) exit
               taken = taken + 1
               reached(nlist%nlat(start + ii)) = .true.
            end do

            call brute_force_row(xyz, i, t, ref_ids, ref_dist)

            call check(error, taken == size(ref_ids), &
                       "Early exit visited the wrong number of neighbours")
            if (allocated(error)) return

            do ii = 1, size(ref_ids)
               call check(error, reached(ref_ids(ii)), &
                          "Early exit skipped a neighbour inside the threshold")
               if (allocated(error)) return
            end do
         end do
      end do

      call nlist%destroy()
   end subroutine test_sorted_early_exit_contract

   !=========================================================================!
   ! Cell grid                                                               !
   !=========================================================================!

   !> Many-cell brute-force comparison.
   !>
   !> The small hand-written fixtures span barely two cells per axis, which
   !> leaves the 27-cell stencil, its boundary clipping and the cell-index
   !> decode effectively unexercised. Here the box is six cutoffs wide.
   subroutine test_many_cells_bruteforce(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp), allocatable :: xyz(:, :)

      call random_points(xyz, 300, 6.0_wp)

      call nlist%init(cutoff=1.0_wp, sorted=.true.)
      call nlist%update(xyz)

      call check_csr_invariants(error, nlist, size(xyz, 2))
      if (allocated(error)) return

      ! Guard the fixture itself: a box this dense must produce real neighbours,
      ! otherwise the comparison below would pass vacuously.
      call check(error, size(nlist%nlat) > 300, &
                 "Fixture is too sparse to exercise the cell stencil")
      if (allocated(error)) return

      call check_all_rows_bruteforce(error, nlist, xyz)
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_many_cells_bruteforce

   !> A cutoff far larger than the bounding box collapses to a single cell.
   subroutine test_cutoff_exceeds_extent(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 5)
      integer :: i

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [1.0_wp, 0.5_wp, -0.5_wp]
      xyz(:, 3) = [-0.75_wp, 1.25_wp, 0.25_wp]
      xyz(:, 4) = [0.5_wp, -1.0_wp, 1.0_wp]
      xyz(:, 5) = [-1.0_wp, -1.0_wp, -1.0_wp]

      call nlist%init(cutoff=100.0_wp)
      call nlist%update(xyz)

      call check_csr_invariants(error, nlist, size(xyz, 2))
      if (allocated(error)) return

      do i = 1, size(xyz, 2)
         call check(error, nlist%nnl(i) == size(xyz, 2) - 1, &
                    "Every point must neighbour every other point")
         if (allocated(error)) return
      end do

      call nlist%destroy()
   end subroutine test_cutoff_exceeds_extent

   !> Degenerate grid extents: a flat sheet (nz == 1) and a line (ny == nz == 1).
   subroutine test_planar_and_collinear(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: plane(3, 9), line(3, 6)
      integer :: ix, iy, k

      ! 3x3 sheet at z = 0, spacing 1.0
      k = 0
      do ix = 0, 2
         do iy = 0, 2
            k = k + 1
            plane(:, k) = [real(ix, wp), real(iy, wp), 0.0_wp]
         end do
      end do

      call nlist%init(cutoff=1.2_wp, sorted=.true.)
      call nlist%update(plane)

      call check_csr_invariants(error, nlist, size(plane, 2))
      if (allocated(error)) return
      call check_all_rows_bruteforce(error, nlist, plane)
      if (allocated(error)) return

      ! Six collinear points along x, spacing 0.7
      do k = 1, 6
         line(:, k) = [0.7_wp*real(k - 1, wp), 0.0_wp, 0.0_wp]
      end do

      call nlist%init(cutoff=1.5_wp, sorted=.true.)
      call nlist%update(line)

      call check_csr_invariants(error, nlist, size(line, 2))
      if (allocated(error)) return
      call check_all_rows_bruteforce(error, nlist, line)
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_planar_and_collinear

   !> Coincident points still pair up, at zero distance, without self-pairs.
   subroutine test_coincident_points(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 4)
      integer :: i, k

      xyz(:, 1) = [1.0_wp, 2.0_wp, 3.0_wp]
      xyz(:, 2) = [1.0_wp, 2.0_wp, 3.0_wp]
      xyz(:, 3) = [1.0_wp, 2.0_wp, 3.0_wp]
      xyz(:, 4) = [1.0_wp, 2.0_wp, 3.0_wp]

      call nlist%init(cutoff=1.0_wp, sorted=.true.)
      call nlist%update(xyz)

      call check_csr_invariants(error, nlist, size(xyz, 2))
      if (allocated(error)) return

      do i = 1, size(xyz, 2)
         call check(error, nlist%nnl(i) == 3, &
                    "Each coincident point must see the other three")
         if (allocated(error)) return

         do k = nlist%inl(i) + 1, nlist%inl(i) + nlist%nnl(i)
            call check(error, nlist%dist(k) == 0.0_wp, &
                       "Coincident points must be stored at zero distance")
            if (allocated(error)) return
         end do
      end do

      call nlist%destroy()
   end subroutine test_coincident_points

   !> Adjacency is symmetric: j neighbours i exactly when i neighbours j.
   subroutine test_pair_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp), allocatable :: xyz(:, :)
      integer :: i, j, k
      logical :: i_in_j

      call random_points(xyz, 300, 6.0_wp)

      call nlist%init(cutoff=1.0_wp)
      call nlist%update(xyz)

      call check(error, mod(size(nlist%nlat), 2) == 0, &
                 "Total pair count must be even")
      if (allocated(error)) return

      do i = 1, size(xyz, 2)
         do k = nlist%inl(i) + 1, nlist%inl(i) + nlist%nnl(i)
            j = nlist%nlat(k)
            i_in_j = any(nlist%nlat(nlist%inl(j) + 1:nlist%inl(j) + nlist%nnl(j)) == i)
            call check(error, i_in_j, "Adjacency must be symmetric")
            if (allocated(error)) return
         end do
      end do

      call nlist%destroy()
   end subroutine test_pair_symmetry

   !=========================================================================!
   ! Degenerate input and object lifecycle                                   !
   !=========================================================================!

   !> Zero and one point must produce a valid, empty list.
   subroutine test_empty_and_single_point(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: empty(3, 0), single(3, 1)
      integer, allocatable :: ids(:)

      call nlist%init(cutoff=1.0_wp, sorted=.true.)
      call nlist%update(empty)

      call check_csr_invariants(error, nlist, 0)
      if (allocated(error)) return

      single(:, 1) = [0.25_wp, -0.5_wp, 2.0_wp]
      call nlist%update(single)

      call check_csr_invariants(error, nlist, 1)
      if (allocated(error)) return
      call check(error, nlist%nnl(1) == 0, "A lone point has no neighbours")
      if (allocated(error)) return

      ids = nlist%get_neighbours(1)
      call check(error, size(ids) == 0, "get_neighbours must return an empty array")
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_empty_and_single_point

   !> A non-positive cutoff yields an empty list rather than a crash.
   subroutine test_nonpositive_cutoff(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: untouched, negative
      real(wp) :: xyz(3, 4)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [0.1_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 0.1_wp, 0.0_wp]
      xyz(:, 4) = [0.0_wp, 0.0_wp, 0.1_wp]

      ! No init at all: the default cutoff of zero must short-circuit.
      call untouched%update(xyz)
      call check_csr_invariants(error, untouched, size(xyz, 2))
      if (allocated(error)) return
      call check(error, all(untouched%nnl == 0), "A zero cutoff must yield no neighbours")
      if (allocated(error)) return
      call check(error, size(untouched%nlat) == 0, "A zero cutoff must yield an empty nlat")
      if (allocated(error)) return

      call negative%init(cutoff=-1.0_wp)
      call negative%update(xyz)
      call check(error, all(negative%nnl == 0), "A negative cutoff must yield no neighbours")
      if (allocated(error)) return
      call check(error, size(negative%nlat) == 0, "A negative cutoff must yield an empty nlat")
      if (allocated(error)) return

      call untouched%destroy()
      call negative%destroy()
   end subroutine test_nonpositive_cutoff

   !> Repeated update on a populated object must resize cleanly.
   !>
   !> No production caller reaches this path (both always re-init first), so the
   !> reallocation logic is only covered here.
   subroutine test_rebuild_grow_and_shrink(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: reused, fresh
      real(wp) :: small(3, 4), big(3, 9), tiny(3, 2)
      integer :: k

      small(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      small(:, 2) = [0.6_wp, 0.0_wp, 0.0_wp]
      small(:, 3) = [0.0_wp, 0.6_wp, 0.0_wp]
      small(:, 4) = [2.5_wp, 2.5_wp, 2.5_wp]

      do k = 1, 9
         big(:, k) = [0.4_wp*real(k - 1, wp), 0.2_wp*real(mod(k, 3), wp), 0.0_wp]
      end do

      tiny(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      tiny(:, 2) = [0.3_wp, 0.0_wp, 0.0_wp]

      call reused%init(cutoff=1.0_wp, sorted=.true.)

      call reused%update(small)
      call fresh%init(cutoff=1.0_wp, sorted=.true.)
      call fresh%update(small)
      call check_lists_identical(error, reused, fresh, "first build")
      if (allocated(error)) return

      ! Grow
      call reused%update(big)
      call fresh%init(cutoff=1.0_wp, sorted=.true.)
      call fresh%update(big)
      call check_lists_identical(error, reused, fresh, "after growing")
      if (allocated(error)) return

      ! Shrink
      call reused%update(tiny)
      call fresh%init(cutoff=1.0_wp, sorted=.true.)
      call fresh%update(tiny)
      call check_lists_identical(error, reused, fresh, "after shrinking")
      if (allocated(error)) return

      call reused%destroy()
      call fresh%destroy()
   end subroutine test_rebuild_grow_and_shrink

   !> destroy must be idempotent and must leave the object reusable.
   subroutine test_destroy_and_reuse(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: nlist
      real(wp) :: xyz(3, 4)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [0.5_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 0.5_wp, 0.0_wp]
      xyz(:, 4) = [5.0_wp, 5.0_wp, 5.0_wp]

      call nlist%init(cutoff=1.0_wp, sorted=.true.)
      call nlist%update(xyz)
      call check(error, nlist%n == 4, "n must track the point count")
      if (allocated(error)) return

      call nlist%destroy()
      call nlist%destroy()

      call check(error, nlist%n == 0, "destroy must reset the point count")
      if (allocated(error)) return
      call check(error,.not. allocated(nlist%inl), "destroy must release inl")
      if (allocated(error)) return
      call check(error,.not. allocated(nlist%nnl), "destroy must release nnl")
      if (allocated(error)) return
      call check(error,.not. allocated(nlist%nlat), "destroy must release nlat")
      if (allocated(error)) return
      call check(error,.not. allocated(nlist%dist), "destroy must release dist")
      if (allocated(error)) return

      ! The cutoff survives destroy by design, so update alone rebuilds the list.
      call nlist%update(xyz)
      call check_csr_invariants(error, nlist, size(xyz, 2))
      if (allocated(error)) return
      call check_all_rows_bruteforce(error, nlist, xyz)
      if (allocated(error)) return

      call nlist%destroy()
   end subroutine test_destroy_and_reuse

   !> init fully defines the configuration: omitting sorted means unsorted, even
   !> if the instance was previously initialised as sorted.
   subroutine test_init_resets_sorted(error)
      type(error_type), allocatable, intent(out) :: error
      type(adjacency_list_type) :: reused, fresh
      real(wp) :: xyz(3, 8)

      call star_cluster(xyz)

      call reused%init(cutoff=1.6_wp, sorted=.true.)
      call reused%update(xyz)
      call check(error, reused%sorted, "sorted=.true. must be honoured")
      if (allocated(error)) return

      ! Re-init without the optional argument.
      call reused%init(cutoff=1.6_wp)
      call reused%update(xyz)

      call check(error,.not. reused%sorted, "init without sorted must reset the flag")
      if (allocated(error)) return

      call fresh%init(cutoff=1.6_wp)
      call fresh%update(xyz)
      call check_lists_identical(error, reused, fresh, "re-init as unsorted")
      if (allocated(error)) return

      call reused%destroy()
      call fresh%destroy()
   end subroutine test_init_resets_sorted

   !=========================================================================!
   ! Helpers                                                                 !
   !=========================================================================!

   !> Central point surrounded by seven others at pairwise distinct distances.
   pure subroutine star_cluster(xyz)
      real(wp), intent(out) :: xyz(3, 8)

      xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [0.5_wp, 0.0_wp, 0.0_wp]
      xyz(:, 3) = [0.0_wp, 1.25_wp, 0.0_wp]
      xyz(:, 4) = [0.0_wp, 0.0_wp, 0.75_wp]
      xyz(:, 5) = [1.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 6) = [-1.5_wp, 0.0_wp, 0.0_wp]
      xyz(:, 7) = [0.0_wp, -0.25_wp, 0.0_wp]
      xyz(:, 8) = [0.0_wp, 0.0_wp, -1.375_wp]
   end subroutine star_cluster

   !> Fill a coordinate array with points drawn from a pinned random stream.
   !>
   !> The seed follows the convention in test_math_sorters, so the fixture is
   !> reproducible for a given compiler. Assertions built on it are all
   !> comparisons against a brute-force reference, so a differing pseudo-random
   !> stream changes the fixture but never the verdict.
   subroutine random_points(xyz, npoints, box)
      real(wp), allocatable, intent(out) :: xyz(:, :)
      !> Number of points to generate
      integer, intent(in) :: npoints
      !> Side length of the cubic box the points are drawn from
      real(wp), intent(in) :: box

      integer, allocatable :: seed(:)
      integer :: seed_size, i

      call random_seed(size=seed_size)
      allocate (seed(seed_size))
      do i = 1, seed_size
         seed(i) = 420 + 69*i
      end do
      call random_seed(put=seed)
      deallocate (seed)

      allocate (xyz(3, npoints))
      call random_number(xyz)
      xyz = box*xyz
   end subroutine random_points

   !> Brute-force reference for the neighbours of point i within cutoff.
   subroutine brute_force_row(xyz, i, cutoff, ids, dists)
      real(wp), intent(in) :: xyz(:, :)
      !> Query point index
      integer, intent(in) :: i
      !> Interaction cutoff distance
      real(wp), intent(in) :: cutoff
      !> Neighbour indices in ascending index order
      integer, allocatable, intent(out) :: ids(:)
      !> Distances parallel to ids
      real(wp), allocatable, intent(out) :: dists(:)

      integer :: j, npoints, nb
      real(wp) :: d2, cutoff2, tol

      npoints = size(xyz, 2)
      cutoff2 = cutoff*cutoff
      tol = 10.0_wp*epsilon(cutoff2)*max(1.0_wp, cutoff2)

      nb = 0
      do j = 1, npoints
         if (j == i) cycle
         d2 = sum((xyz(:, i) - xyz(:, j))**2)
         if (d2 <= cutoff2 + tol) nb = nb + 1
      end do

      allocate (ids(nb), dists(nb))

      nb = 0
      do j = 1, npoints
         if (j == i) cycle
         d2 = sum((xyz(:, i) - xyz(:, j))**2)
         if (d2 <= cutoff2 + tol) then
            nb = nb + 1
            ids(nb) = j
            dists(nb) = sqrt(d2)
         end if
      end do
   end subroutine brute_force_row

   !> Assert the raw CSR contract every production consumer walks directly.
   subroutine check_csr_invariants(error, nlist, npoints)
      type(error_type), allocatable, intent(inout) :: error
      !> List under test
      type(adjacency_list_type), intent(in) :: nlist
      !> Expected number of points
      integer, intent(in) :: npoints

      integer :: i, k, total

      call check(error, nlist%n == npoints, "n must equal the number of points")
      if (allocated(error)) return

      call check(error, allocated(nlist%inl) .and. allocated(nlist%nnl) .and. &
                 allocated(nlist%nlat) .and. allocated(nlist%dist), &
                 "update must leave every CSR array allocated")
      if (allocated(error)) return

      call check(error, size(nlist%inl) == npoints .and. size(nlist%nnl) == npoints, &
                 "inl and nnl must be sized by the number of points")
      if (allocated(error)) return

      total = 0
      do i = 1, npoints
         call check(error, nlist%nnl(i) >= 0, "Neighbour counts must be non-negative")
         if (allocated(error)) return
         call check(error, nlist%inl(i) == total, &
                    "inl must be the exclusive prefix sum of nnl")
         if (allocated(error)) return
         total = total + nlist%nnl(i)
      end do

      call check(error, size(nlist%nlat) == total, "nlat must hold exactly sum(nnl) entries")
      if (allocated(error)) return
      call check(error, size(nlist%dist) == total, "dist must be parallel to nlat")
      if (allocated(error)) return

      do i = 1, npoints
         do k = nlist%inl(i) + 1, nlist%inl(i) + nlist%nnl(i)
            call check(error, nlist%nlat(k) >= 1 .and. nlist%nlat(k) <= npoints, &
                       "Neighbour index out of range")
            if (allocated(error)) return
            call check(error, nlist%nlat(k) /= i, "Self index must never appear")
            if (allocated(error)) return
            call check(error, count(nlist%nlat(nlist%inl(i) + 1: &
                                               nlist%inl(i) + nlist%nnl(i)) == nlist%nlat(k)) == 1, &
                       "Duplicate neighbour index within a row")
            if (allocated(error)) return
            call check(error, nlist%dist(k) >= 0.0_wp .and. &
                       nlist%dist(k) <= nlist%cutoff*(1.0_wp + thr), &
                       "Stored distance outside [0, cutoff]")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_csr_invariants

   !> Compare every row against the brute-force reference set.
   subroutine check_all_rows_bruteforce(error, nlist, xyz)
      type(error_type), allocatable, intent(inout) :: error
      !> List under test
      type(adjacency_list_type), intent(in) :: nlist
      !> Coordinates the list was built from
      real(wp), intent(in) :: xyz(:, :)

      integer, allocatable :: ref_ids(:)
      real(wp), allocatable :: ref_dist(:)
      integer :: i, k, start, cnt

      do i = 1, size(xyz, 2)
         call brute_force_row(xyz, i, nlist%cutoff, ref_ids, ref_dist)
         start = nlist%inl(i)
         cnt = nlist%nnl(i)

         call check(error, cnt == size(ref_ids), &
                    "Neighbour count mismatch against brute-force reference")
         if (allocated(error)) return

         ! Together with the no-duplicate invariant, matching counts plus
         ! one-way containment give set equality.
         do k = 1, size(ref_ids)
            call check(error, any(nlist%nlat(start + 1:start + cnt) == ref_ids(k)), &
                       "Brute-force neighbour missing from the list")
            if (allocated(error)) return
         end do
      end do
   end subroutine check_all_rows_bruteforce

   !> Assert two lists agree in every stored field.
   subroutine check_lists_identical(error, lhs, rhs, label)
      type(error_type), allocatable, intent(inout) :: error
      !> Lists to compare
      type(adjacency_list_type), intent(in) :: lhs, rhs
      !> Context prefix for failure messages
      character(*), intent(in) :: label

      call check(error, lhs%n == rhs%n, label//": point count mismatch")
      if (allocated(error)) return
      call check(error, size(lhs%nlat) == size(rhs%nlat), label//": pair count mismatch")
      if (allocated(error)) return
      call check(error, all(lhs%nnl == rhs%nnl), label//": nnl mismatch")
      if (allocated(error)) return
      call check(error, all(lhs%inl == rhs%inl), label//": inl mismatch")
      if (allocated(error)) return
      call check(error, all(lhs%nlat == rhs%nlat), label//": nlat mismatch")
      if (allocated(error)) return
      call check(error, all(abs(lhs%dist - rhs%dist) <= thr), label//": dist mismatch")
   end subroutine check_lists_identical

end module test_math_adjacency_list
