!> Test suite for the binned cell list in moist_cavity_drop_lsf_screening
module test_cavity_drop_lsf_screening
   use mctc_env, only: wp
   use test_helpers, only: det_rand
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_cavity_drop_lsf_screening, only: cell_list_type
   use moist_utils_format, only: to_string
   implicit none(type, external)
   private

   public :: collect_cavity_drop_lsf_screening

   !> Atoms in the random fixture
   integer, parameter :: nat_fixture = 300
   !> Random queries per test
   integer, parameter :: nquery = 400

contains

   !> Collect all cell-list tests
   subroutine collect_cavity_drop_lsf_screening(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("box_vs_bruteforce", test_box_vs_bruteforce), &
                  new_unittest("sphere_vs_bruteforce", test_sphere_vs_bruteforce), &
                  new_unittest("point_and_degenerate_box", test_point_and_degenerate_box), &
                  new_unittest("far_and_enclosing_boxes", test_far_and_enclosing_boxes), &
                  new_unittest("full_scan_matches_grid", test_full_scan_matches_grid), &
                  new_unittest("lattice_on_cell_faces", test_lattice_on_cell_faces), &
                  new_unittest("saturated_reach", test_saturated_reach), &
                  new_unittest("sparse_geometry_memory", test_sparse_geometry), &
                  new_unittest("huge_query_returns_all", test_huge_query), &
                  new_unittest("empty_and_single", test_empty_and_single), &
                  new_unittest("tail_and_inplace_sort_agree", test_tail_and_inplace_sort), &
                  new_unittest("box_bounds_matches_query_box", test_box_bounds) &
                  ]
   end subroutine collect_cavity_drop_lsf_screening

   !* ================================================================================= *!
   !*                                     Helpers                                       *!
   !* ================================================================================= *!

   !> Random atoms in a 30 bohr cube with reach in [0.5, 4.5)
   subroutine random_atoms(n, xyz, reach, state)
      !> Number of atoms
      integer, intent(in) :: n
      !> Atom centers
      real(wp), allocatable, intent(out) :: xyz(:, :)
      !> Reach per atom
      real(wp), allocatable, intent(out) :: reach(:)
      !> Generator state
      integer, intent(inout) :: state

      allocate (xyz(3, n), reach(n))
      call det_rand(xyz, state, 0.0_wp, 30.0_wp)
      call det_rand(reach, state, 0.5_wp, 4.5_wp)
   end subroutine random_atoms

   !> Random box with center in [-10, 40) and half-widths in [0, 12); some axes flat
   subroutine random_box(lo, hi, state)
      !> Lower corner
      real(wp), intent(out) :: lo(3)
      !> Upper corner
      real(wp), intent(out) :: hi(3)
      !> Generator state
      integer, intent(inout) :: state

      real(wp) :: u(7), c(3), h(3)

      call det_rand(u, state)
      c = -10.0_wp + 50.0_wp*u(1:3)
      h = 12.0_wp*u(4:6)
      if (u(7) < 0.2_wp) h(1) = 0.0_wp
      if (u(7) > 0.9_wp) h = 0.0_wp
      lo = c - h
      hi = c + h
   end subroutine random_box

   !> Brute-force atoms whose reach sphere meets a box, ascending
   subroutine brute_box(xyz, reach, lo, hi, ref, nref)
      !> Atom centers
      real(wp), intent(in) :: xyz(:, :)
      !> Reach per atom
      real(wp), intent(in) :: reach(:)
      !> Lower corner
      real(wp), intent(in) :: lo(3)
      !> Upper corner
      real(wp), intent(in) :: hi(3)
      !> Reference ids
      integer, intent(out) :: ref(:)
      !> Number of reference ids
      integer, intent(out) :: nref

      real(wp) :: d(3)
      integer :: j

      nref = 0
      do j = 1, size(reach)
         d = max(lo - xyz(:, j), 0.0_wp, xyz(:, j) - hi)
         if (sum(d*d) <= reach(j)*reach(j)) then
            nref = nref + 1
            ref(nref) = j
         end if
      end do
   end subroutine brute_box

   !> Brute-force atoms whose reach sphere meets a ball, ascending
   subroutine brute_sphere(xyz, reach, center, r, ref, nref)
      !> Atom centers
      real(wp), intent(in) :: xyz(:, :)
      !> Reach per atom
      real(wp), intent(in) :: reach(:)
      !> Ball center
      real(wp), intent(in) :: center(3)
      !> Ball radius
      real(wp), intent(in) :: r
      !> Reference ids
      integer, intent(out) :: ref(:)
      !> Number of reference ids
      integer, intent(out) :: nref

      integer :: j

      nref = 0
      do j = 1, size(reach)
         if (sum((xyz(:, j) - center)**2) <= (r + reach(j))**2) then
            nref = nref + 1
            ref(nref) = j
         end if
      end do
   end subroutine brute_sphere

   !> Require identical id lists, ascending
   subroutine check_same(error, out, nout, ref, nref, what)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Query result
      integer, intent(in) :: out(:)
      !> Query count
      integer, intent(in) :: nout
      !> Reference ids
      integer, intent(in) :: ref(:)
      !> Reference count
      integer, intent(in) :: nref
      !> Context for the diagnostic
      character(len=*), intent(in) :: what

      integer :: i

      call check(error, nout, nref, more=what//": count")
      if (allocated(error)) return
      do i = 2, nout
         if (out(i) <= out(i - 1)) then
            call test_failed(error, what//": result not strictly ascending")
            return
         end if
      end do
      if (any(out(1:nout) /= ref(1:nref))) then
         call test_failed(error, what//": ids differ from brute force")
      end if
   end subroutine check_same

   !* ================================================================================= *!
   !*                                      Tests                                        *!
   !* ================================================================================= *!

   !> Random boxes, including many-cell spans and flat boxes, equal brute force
   subroutine test_box_vs_bruteforce(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: lo(3), hi(3)
      integer :: out(nat_fixture), ref(nat_fixture), nout, nref, iq
      integer :: state

      state = 7
      call random_atoms(nat_fixture, xyz, reach, state)
      call cells%build(xyz, reach, full_scan_below=0)
      call check(error, .not. cells%full_scan, more="fixture must use the binned grid")
      if (allocated(error)) return
      call check(error, cells%ncells > 8, more="fixture must have many cells")
      if (allocated(error)) return

      do iq = 1, nquery
         call random_box(lo, hi, state)
         call cells%query_box(lo, hi, out, nout)
         call brute_box(xyz, reach, lo, hi, ref, nref)
         call check_same(error, out, nout, ref, nref, "query_box")
         if (allocated(error)) return
      end do
   end subroutine test_box_vs_bruteforce

   !> A result buffer with room for the merge sort scratch and one without,
   !> which sorts in place, return the same ids
   !>
   !> - the merge sort runs in the buffer's tail while that holds half the
   !>   result, else the introsort
   !> - the sample must hold many-cell results on both sorts
   subroutine test_tail_and_inplace_sort(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: lo(3), hi(3), h(3)
      integer :: wide(2*nat_fixture), tight(nat_fixture), ref(nat_fixture)
      integer :: nwide, ntight, nref, iq, state, n_tail, n_intro

      state = 11
      call random_atoms(nat_fixture, xyz, reach, state)
      call cells%build(xyz, reach, full_scan_below=0)
      n_tail = 0
      n_intro = 0
      do iq = 1, nquery
         call random_box(lo, hi, state)
         ! Every fourth box around the middle of the cube, most of the atoms in it
         if (mod(iq, 4) == 0) then
            call det_rand(h, state, 6.0_wp, 15.0_wp)
            lo = 15.0_wp - h
            hi = 15.0_wp + h
         end if
         call cells%query_box(lo, hi, wide, nwide)
         call cells%query_box(lo, hi, tight, ntight)
         call brute_box(xyz, reach, lo, hi, ref, nref)
         call check_same(error, wide, nwide, ref, nref, "query_box with room for scratch")
         if (allocated(error)) return
         call check_same(error, tight, ntight, ref, nref, "query_box with the exact buffer")
         if (allocated(error)) return
         if (ntight <= 1) cycle
         if (nat_fixture - ntight >= ntight/2) then
            n_tail = n_tail + 1
         else
            n_intro = n_intro + 1
         end if
      end do
      call check(error, n_tail > nquery/10 .and. n_intro > nquery/100, &
                 "sample does not reach both sorts: "//to_string(n_tail) &
                 //" merge-sorted in the tail, "//to_string(n_intro)//" introsorted")
   end subroutine test_tail_and_inplace_sort

   !> Random balls equal brute force
   subroutine test_sphere_vs_bruteforce(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: u(4), center(3), r
      integer :: out(nat_fixture), ref(nat_fixture), nout, nref, iq
      integer :: state

      state = 11
      call random_atoms(nat_fixture, xyz, reach, state)
      call cells%build(xyz, reach, full_scan_below=0)

      do iq = 1, nquery
         call det_rand(u, state)
         center = -10.0_wp + 50.0_wp*u(1:3)
         r = 10.0_wp*u(4)
         if (iq <= 20) r = 0.0_wp
         call cells%query_sphere(center, r, out, nout)
         call brute_sphere(xyz, reach, center, r, ref, nref)
         call check_same(error, out, nout, ref, nref, "query_sphere")
         if (allocated(error)) return
      end do
   end subroutine test_sphere_vs_bruteforce

   !> Point query equals the degenerate box; swapped corners equal ordered ones
   subroutine test_point_and_degenerate_box(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: p(3), lo(3), hi(3)
      integer :: out(nat_fixture), out2(nat_fixture), ref(nat_fixture)
      integer :: nout, nout2, nref, iq
      integer :: state

      state = 13
      call random_atoms(nat_fixture, xyz, reach, state)
      call cells%build(xyz, reach, full_scan_below=0)

      do iq = 1, nquery
         call det_rand(p, state, -5.0_wp, 35.0_wp)
         call cells%query_point(p, out, nout)
         call brute_box(xyz, reach, p, p, ref, nref)
         call check_same(error, out, nout, ref, nref, "query_point")
         if (allocated(error)) return
         call cells%query_box(p, p, out2, nout2)
         call check_same(error, out2, nout2, out, nout, "degenerate query_box")
         if (allocated(error)) return

         call random_box(lo, hi, state)
         call cells%query_box(lo, hi, out, nout)
         call cells%query_box(hi, lo, out2, nout2)
         call check_same(error, out2, nout2, out, nout, "swapped corners")
         if (allocated(error)) return
      end do
   end subroutine test_point_and_degenerate_box

   !> Boxes far outside return nothing; an enclosing box returns every atom
   subroutine test_far_and_enclosing_boxes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: offsets(4), lo(3), hi(3)
      integer :: out(nat_fixture), nout, io, ax, sgn
      integer :: state

      state = 17
      call random_atoms(nat_fixture, xyz, reach, state)
      call cells%build(xyz, reach, full_scan_below=0)

      offsets = [50.0_wp, 1.0e3_wp, 1.0e30_wp, 1.0e300_wp]
      do io = 1, size(offsets)
         do ax = 1, 3
            do sgn = -1, 1, 2
               lo = 15.0_wp
               hi = 16.0_wp
               lo(ax) = lo(ax) + sgn*offsets(io)
               hi(ax) = hi(ax) + sgn*offsets(io)
               call cells%query_box(lo, hi, out, nout)
               call check(error, nout, 0, more="far box must be empty")
               if (allocated(error)) return
               call cells%query_sphere(lo, 1.0_wp, out, nout)
               call check(error, nout, 0, more="far ball must be empty")
               if (allocated(error)) return
            end do
         end do
      end do

      lo = -1.0e300_wp
      hi = 1.0e300_wp
      call cells%query_box(lo, hi, out, nout)
      call check(error, nout, nat_fixture, more="enclosing box must return every atom")
      if (allocated(error)) return
      call check(error, all(out == [(io, io=1, nat_fixture)]), more="enclosing box ids")
      if (allocated(error)) return

      lo = 0.0_wp
      hi = 30.0_wp
      call cells%query_box(lo, hi, out, nout)
      call check(error, nout, nat_fixture, more="bounding box must return every atom")
   end subroutine test_far_and_enclosing_boxes

   !> Single-cell mode below `full_scan_below` returns the grid's results bit for bit
   subroutine test_full_scan_matches_grid(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: grid, flat
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: lo(3), hi(3), u
      integer :: out(nat_fixture), out2(nat_fixture), nout, nout2, iq
      integer :: state

      state = 19
      call random_atoms(nat_fixture, xyz, reach, state)
      call grid%build(xyz, reach, full_scan_below=nat_fixture)
      call flat%build(xyz, reach, full_scan_below=nat_fixture + 1)
      call check(error, .not. grid%full_scan, more="n == full_scan_below must bin")
      if (allocated(error)) return
      call check(error, flat%full_scan .and. flat%ncells == 1, &
                 more="n < full_scan_below must use one cell")
      if (allocated(error)) return

      do iq = 1, nquery
         call random_box(lo, hi, state)
         call grid%query_box(lo, hi, out, nout)
         call flat%query_box(lo, hi, out2, nout2)
         call check_same(error, out2, nout2, out, nout, "full scan box")
         if (allocated(error)) return
         call det_rand(u, state, 0.0_wp, 6.0_wp)
         call grid%query_sphere(lo, u, out, nout)
         call flat%query_sphere(lo, u, out2, nout2)
         call check_same(error, out2, nout2, out, nout, "full scan ball")
         if (allocated(error)) return
      end do
   end subroutine test_full_scan_matches_grid

   !> Atoms on cell faces: unit lattice, unit reach, integer and half-integer box faces
   !>
   !> - origin 0 and side 1 put every center on a lower cell face
   !> - integer faces at distance 1 from a lattice plane give `d2 == reach_sq` exactly
   subroutine test_lattice_on_cell_faces(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nside = 7
      integer, parameter :: nlat = nside**3

      type(cell_list_type) :: cells
      real(wp) :: xyz(3, nlat), reach(nlat), lo(3), hi(3), u(6)
      integer :: out(nlat), ref(nlat), nout, nref, i, j, k, n, iq
      integer :: state

      n = 0
      do k = 0, nside - 1
         do j = 0, nside - 1
            do i = 0, nside - 1
               n = n + 1
               xyz(:, n) = real([i, j, k], wp)
            end do
         end do
      end do
      reach = 1.0_wp
      call cells%build(xyz, reach, full_scan_below=0)
      call check(error, .not. cells%full_scan .and. cells%ncells == nlat, &
                 more="lattice must bin into one cell per atom")
      if (allocated(error)) return

      state = 37
      do iq = 1, nquery
         call det_rand(u, state)
         ! Integer faces in [-2, 8], half-integer on odd queries, flat boxes every fifth
         lo = real(floor(11.0_wp*u(1:3)), wp) - 2.0_wp
         hi = lo + real(floor(4.0_wp*u(4:6)), wp)
         if (mod(iq, 2) == 1) then
            lo = lo + 0.5_wp
            hi = hi + 0.5_wp
         end if
         if (mod(iq, 5) == 0) hi = lo
         call cells%query_box(lo, hi, out, nout)
         call brute_box(xyz, reach, lo, hi, ref, nref)
         call check_same(error, out, nout, ref, nref, "lattice query_box")
         if (allocated(error)) return
      end do

      ! Face exactly one unit beyond the plane x = 6 touches every atom on it
      lo = [7.0_wp, 0.0_wp, 0.0_wp]
      hi = [8.0_wp, 6.0_wp, 6.0_wp]
      call cells%query_box(lo, hi, out, nout)
      call brute_box(xyz, reach, lo, hi, ref, nref)
      call check(error, nref, nside*nside, more="touching face reference")
      if (allocated(error)) return
      call check_same(error, out, nout, ref, nref, "touching face")
      if (allocated(error)) return

      ! Lattice points as degenerate boxes hit the atom and its six neighbors
      lo = 3.0_wp
      call cells%query_point(lo, out, nout)
      call brute_box(xyz, reach, lo, lo, ref, nref)
      call check(error, nref, 7, more="lattice point reference")
      if (allocated(error)) return
      call check_same(error, out, nout, ref, nref, "lattice point")
   end subroutine test_lattice_on_cell_faces

   !> Unbounded reach degenerates to one cell and screens nothing
   subroutine test_saturated_reach(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: lo(3), hi(3)
      integer :: out(nat_fixture), ref(nat_fixture), nout, nref, iq
      integer :: state

      state = 23
      call random_atoms(nat_fixture, xyz, reach, state)
      reach = huge(0.0_wp)
      call cells%build(xyz, reach, full_scan_below=0)
      call check(error, cells%full_scan, more="huge reach must degenerate")
      if (allocated(error)) return
      lo = 1.0e100_wp
      hi = 1.0e100_wp
      call cells%query_box(lo, hi, out, nout)
      call check(error, nout, nat_fixture, more="huge reach screens nothing (box)")
      if (allocated(error)) return
      call cells%query_sphere(lo, 1.0_wp, out, nout)
      call check(error, nout, nat_fixture, more="huge reach screens nothing (ball)")
      if (allocated(error)) return

      ! One saturated atom among finite reaches: still exact for the others
      call random_atoms(nat_fixture, xyz, reach, state)
      reach(5) = huge(0.0_wp)
      call cells%build(xyz, reach, full_scan_below=0)
      call check(error, cells%full_scan, more="one huge reach must degenerate")
      if (allocated(error)) return
      do iq = 1, nquery/4
         call random_box(lo, hi, state)
         call cells%query_box(lo, hi, out, nout)
         ! huge*huge is +inf, so the brute force keeps atom 5 by itself
         call brute_box(xyz, reach, lo, hi, ref, nref)
         call check_same(error, out, nout, ref, nref, "mixed saturated box")
         if (allocated(error)) return
      end do
   end subroutine test_saturated_reach

   !> Two clusters far apart keep the cell count O(n) and stay exact
   subroutine test_sparse_geometry(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: lo(3), hi(3)
      integer :: out(nat_fixture), ref(nat_fixture), nout, nref, iq
      integer :: state

      state = 29
      call random_atoms(nat_fixture, xyz, reach, state)
      xyz(:, nat_fixture/2 + 1:) = xyz(:, nat_fixture/2 + 1:) + 1.0e5_wp
      call cells%build(xyz, reach, full_scan_below=0)
      call check(error, cells%ncells <= 8*nat_fixture, more="cell count must stay O(n)")
      if (allocated(error)) return

      do iq = 1, nquery
         call random_box(lo, hi, state)
         if (mod(iq, 2) == 0) then
            lo = lo + 1.0e5_wp
            hi = hi + 1.0e5_wp
         end if
         call cells%query_box(lo, hi, out, nout)
         call brute_box(xyz, reach, lo, hi, ref, nref)
         call check_same(error, out, nout, ref, nref, "sparse query_box")
         if (allocated(error)) return
      end do
   end subroutine test_sparse_geometry

   !> Corners near `+-huge` and an unbounded radius stay finite and return every atom
   !>
   !> - non-finite query input is an error stop, not testable in process
   subroutine test_huge_query(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp), allocatable :: xyz(:, :), reach(:)
      integer :: out(nat_fixture), nout
      integer :: state

      state = 31
      call random_atoms(nat_fixture, xyz, reach, state)
      call cells%build(xyz, reach, full_scan_below=0)
      call cells%query_box([-huge(0.0_wp), -huge(0.0_wp), -huge(0.0_wp)], &
                           [huge(0.0_wp), huge(0.0_wp), huge(0.0_wp)], out, nout)
      call check(error, nout, nat_fixture, more="box of +-huge returns every atom")
      if (allocated(error)) return
      call cells%query_sphere([1.0_wp, 2.0_wp, 3.0_wp], huge(0.0_wp), out, nout)
      call check(error, nout, nat_fixture, more="huge radius returns every atom")
   end subroutine test_huge_query

   !> Empty system answers nothing; a single atom is found only within its reach
   subroutine test_empty_and_single(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: cells
      real(wp) :: xyz0(3, 0), reach0(0), xyz1(3, 1), reach1(1)
      integer :: out(1), nout

      call cells%build(xyz0, reach0, full_scan_below=0)
      call cells%query_box([0.0_wp, 0.0_wp, 0.0_wp], [1.0_wp, 1.0_wp, 1.0_wp], out, nout)
      call check(error, nout, 0, more="empty system box")
      if (allocated(error)) return
      call cells%query_sphere([0.0_wp, 0.0_wp, 0.0_wp], 1.0_wp, out, nout)
      call check(error, nout, 0, more="empty system ball")
      if (allocated(error)) return

      xyz1(:, 1) = [1.0_wp, 2.0_wp, 3.0_wp]
      reach1 = 2.0_wp
      call cells%build(xyz1, reach1, full_scan_below=0)
      call check(error, .not. cells%full_scan .and. cells%ncells == 1, &
                 more="single atom bins into one cell")
      if (allocated(error)) return
      call cells%query_point([1.0_wp, 2.0_wp, 4.9_wp], out, nout)
      call check(error, nout, 1, more="point inside reach")
      if (allocated(error)) return
      call cells%query_point([1.0_wp, 2.0_wp, 5.1_wp], out, nout)
      call check(error, nout, 0, more="point outside reach")
      if (allocated(error)) return
      call cells%query_box([3.5_wp, 2.0_wp, 3.0_wp], [9.0_wp, 9.0_wp, 9.0_wp], out, nout)
      call check(error, nout, 0, more="box face beyond reach")
      if (allocated(error)) return
      call cells%query_sphere([4.5_wp, 2.0_wp, 3.0_wp], 1.6_wp, out, nout)
      call check(error, nout, 1, more="ball overlapping reach")
      if (allocated(error)) return

      call cells%destroy()
      call check(error, cells%natoms == 0 .and. .not. allocated(cells%slot_atom), &
                 more="destroy resets")
   end subroutine test_empty_and_single

   !> `query_box_bounds` matches the exact center bounds of `query_box` in grid and full-scan modes
   subroutine test_box_bounds(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(cell_list_type) :: grid, flat
      real(wp), allocatable :: xyz(:, :), reach(:)
      real(wp) :: lo(3), hi(3), blo_g(3), bhi_g(3), blo_f(3), bhi_f(3), ref_lo(3), ref_hi(3)
      integer :: out(nat_fixture), nout, ng, nf, iq, k, state

      state = 43
      call random_atoms(nat_fixture, xyz, reach, state)
      call grid%build(xyz, reach, full_scan_below=0)
      call flat%build(xyz, reach, full_scan_below=nat_fixture + 1)

      do iq = 1, nquery
         call random_box(lo, hi, state)
         if (mod(iq, 10) == 0) then
            lo = lo + 100.0_wp
            hi = hi + 100.0_wp
         end if
         call grid%query_box(lo, hi, out, nout)
         call grid%query_box_bounds(lo, hi, blo_g, bhi_g, ng)
         call flat%query_box_bounds(lo, hi, blo_f, bhi_f, nf)
         call check(error, ng == nout .and. nf == nout, &
                    more="query_box_bounds count must match query_box")
         if (allocated(error)) return
         if (nout == 0) then
            call check(error, all(blo_g == 0.0_wp) .and. all(bhi_g == 0.0_wp), &
                       more="empty grid query_box_bounds must zero bounds")
            if (allocated(error)) return
            call check(error, all(blo_f == 0.0_wp) .and. all(bhi_f == 0.0_wp), &
                       more="empty full-scan query_box_bounds must zero bounds")
            if (allocated(error)) return
         else
            ref_lo = xyz(:, out(1))
            ref_hi = xyz(:, out(1))
            do k = 2, nout
               ref_lo = min(ref_lo, xyz(:, out(k)))
               ref_hi = max(ref_hi, xyz(:, out(k)))
            end do
            call check(error, all(blo_g == ref_lo) .and. all(bhi_g == ref_hi), &
                       more="grid query_box_bounds must match exact min/max centers")
            if (allocated(error)) return
            call check(error, all(blo_f == ref_lo) .and. all(bhi_f == ref_hi), &
                       more="full-scan query_box_bounds must match exact min/max centers")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_box_bounds

end module test_cavity_drop_lsf_screening
