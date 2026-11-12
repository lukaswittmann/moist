!> Binned cell list with box, sphere and point queries over per-atom reach spheres
!>
!> Each atom is binned into exactly one cell by its center
!>
!> - cell side is `maxval(reach)`, enlarged for sparse geometries so the cell
!>   count stays below `max_cells_per_atom * n`
!> - atom ids inside a cell are ascending (stable counting-sort fill)
!> - a query walks the cells overlapped by the query region expanded by
!>   `maxval(reach)` and tests each atom's own sphere `(center_i, reach_i)`
!> - results are unique and ascending by atom id
!>
!> Degenerate single-cell mode (every query scans all atoms, no index arithmetic):
!>
!> - `n < full_scan_below`
!> - `maxval(reach)` zero, not finite, or `>= sqrt(huge)`
!> - any non-finite center, or a center span that overflows
!>
!> Saturation: a reach that is not finite or `>= sqrt(huge)` screens nothing; its
!> square is stored as `huge`, as in the DROP LSF screening gate
!>
!> Strict superset: each hit test widens the squared reach by `hit_slack`, so a
!> caller's own exact test of the same expression never keeps an atom the list
!> dropped
!>
module moist_math_cell_list
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: moist_cell_list_type
   public :: heapsort

   !> Rows of a slot record: center (1:3), reach, squared reach
   integer, parameter :: slot_rows = 5
   !> Row of the reach inside a slot record
   integer, parameter :: row_reach = 4
   !> Row of the squared reach inside a slot record
   integer, parameter :: row_reach_sq = 5
   !> Upper bound on cells per atom, keeps memory O(n) for sparse geometries
   integer, parameter :: max_cells_per_atom = 8
   !> Relative widening of the cell walk, covers ulp-level rounding at its edge
   !>
   !> - applied to the largest reach, the ball radius and the corner magnitude
   !> - must exceed `hit_slack/2`, the widening of the hit test in distance units
   real(wp), parameter :: walk_slack = 16.0_wp*epsilon(1.0_wp)
   !> Relative widening of the squared reach in the hit test
   !>
   !> - callers repeat the same squared-distance test in another compilation unit,
   !>   where FMA contraction may round differently by an ulp or two
   !> - widening keeps the list a strict superset of any such exact test
   real(wp), parameter :: hit_slack = 16.0_wp*epsilon(1.0_wp)

   !> Binned cell list in CSR layout, slots stored in cell order
   type :: moist_cell_list_type
      !> Number of atoms indexed
      integer :: natoms = 0
      !> Whether every query scans all atoms in id order
      logical :: full_scan = .false.
      !> Cell side length
      real(wp) :: cell_side = 0.0_wp
      !> Inverse cell side
      real(wp) :: inv_cell = 0.0_wp
      !> Largest reach, widened by `walk_slack`
      real(wp) :: walk_reach = 0.0_wp
      !> Lower corner of the center bounding box
      real(wp) :: origin(3) = 0.0_wp
      !> Cells along x
      integer :: nx = 0
      !> Cells along y
      integer :: ny = 0
      !> Cells along z
      integer :: nz = 0
      !> Total number of cells
      integer :: ncells = 0
      !> First slot of each cell, `cell_start(ncells + 1) = natoms + 1` [ncells + 1]
      integer, allocatable :: cell_start(:)
      !> Atom id stored in each slot [natoms]
      integer, allocatable :: slot_atom(:)
      !> Center, reach and squared reach per slot [5, natoms]
      real(wp), allocatable :: slot_geom(:, :)
   contains
      !> Bin atoms by center
      procedure :: build => cell_list_build
      !> Atoms whose reach sphere intersects an axis-aligned box
      procedure :: query_box => cell_list_query_box
      !> Atoms whose reach sphere intersects a ball
      procedure :: query_sphere => cell_list_query_sphere
      !> Atoms whose reach sphere contains a point
      procedure :: query_point => cell_list_query_point
      !> Deallocate all storage
      procedure :: destroy => cell_list_destroy
   end type moist_cell_list_type

contains

   !* ================================================================================= *!
   !*                                      Build                                        *!
   !* ================================================================================= *!

   !> Bin every atom into one cell by its center
   !>
   !> - counting sort over cells, filled in ascending atom order
   !> - error stop on shape mismatch or a negative reach
   !>
   !> @param[in,out] self             cell list, reset on entry
   !> @param[in]     xyz              atom centers (3, n)
   !> @param[in]     reach            reach per atom (n), `>= 0`
   !> @param[in]     full_scan_below  single-cell mode for `n < full_scan_below` (optional)
   subroutine cell_list_build(self, xyz, reach, full_scan_below)
      !> Cell list
      class(moist_cell_list_type), intent(inout) :: self
      !> Atom centers
      real(wp), intent(in) :: xyz(:, :)
      !> Reach per atom
      real(wp), intent(in) :: reach(:)
      !> Atom count below which the list degenerates to one cell
      integer, intent(in), optional :: full_scan_below

      real(wp) :: r_max, sat, span(3), t(3), cells, cap
      integer :: n, j, ic, k
      integer, allocatable :: cell_of(:), fill(:)

      call self%destroy()

      n = size(reach)
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= n) then
         error stop "moist cell list: xyz must have shape (3, size(reach))"
      end if
      if (any(reach < 0.0_wp)) then
         error stop "moist cell list: reach must be non-negative"
      end if

      self%natoms = n
      allocate (self%slot_atom(n))
      allocate (self%slot_geom(slot_rows, n))

      sat = sqrt(huge(0.0_wp))
      r_max = 0.0_wp
      if (n > 0) r_max = maxval(reach)
      span = 0.0_wp
      if (n > 0) then
         self%origin = minval(xyz, dim=2)
         span = maxval(xyz, dim=2) - self%origin
      end if

      ! An overflowing span would never let the cell widening loop terminate
      self%full_scan = n == 0 .or. .not. ieee_is_finite(r_max) .or. r_max <= 0.0_wp &
                       .or. r_max >= sat .or. .not. all(ieee_is_finite(xyz)) &
                       .or. .not. all(ieee_is_finite(span))
      if (present(full_scan_below)) self%full_scan = self%full_scan .or. n < full_scan_below
      ! NaN reach compares false everywhere, so maxval can miss it
      if (.not. all(ieee_is_finite(reach))) self%full_scan = .true.

      if (self%full_scan) then
         self%origin = 0.0_wp
         self%nx = 1
         self%ny = 1
         self%nz = 1
         self%ncells = 1
         allocate (self%cell_start(2))
         self%cell_start = [1, n + 1]
         do j = 1, n
            call store_slot(self, j, j, xyz(:, j), reach(j), sat)
         end do
         return
      end if

      self%cell_side = r_max
      cap = real(max_cells_per_atom, wp)*real(n, wp)
      do
         ! Capped per axis so a tiny side on a wide span cannot overflow
         do k = 1, 3
            if (span(k) > cap*self%cell_side) then
               t(k) = cap
            else
               t(k) = span(k)/self%cell_side
            end if
         end do
         cells = (t(1) + 1.0_wp)*(t(2) + 1.0_wp)*(t(3) + 1.0_wp)
         if (cells <= cap) exit
         self%cell_side = self%cell_side*max(1.5_wp, (cells/cap)**(1.0_wp/3.0_wp))
      end do
      self%inv_cell = 1.0_wp/self%cell_side
      self%walk_reach = r_max*(1.0_wp + walk_slack)

      t = span*self%inv_cell
      self%nx = int(t(1)) + 1
      self%ny = int(t(2)) + 1
      self%nz = int(t(3)) + 1
      self%ncells = self%nx*self%ny*self%nz

      allocate (cell_of(n))
      allocate (self%cell_start(self%ncells + 1), source=0)
      do j = 1, n
         do k = 1, 3
            t(k) = (xyz(k, j) - self%origin(k))*self%inv_cell
         end do
         ic = min(self%nx - 1, int(t(1))) &
              + min(self%ny - 1, int(t(2)))*self%nx &
              + min(self%nz - 1, int(t(3)))*self%nx*self%ny + 1
         cell_of(j) = ic
         self%cell_start(ic + 1) = self%cell_start(ic + 1) + 1
      end do

      self%cell_start(1) = 1
      do ic = 1, self%ncells
         self%cell_start(ic + 1) = self%cell_start(ic + 1) + self%cell_start(ic)
      end do

      allocate (fill(self%ncells))
      fill = self%cell_start(1:self%ncells)
      do j = 1, n
         ic = cell_of(j)
         call store_slot(self, fill(ic), j, xyz(:, j), reach(j), sat)
         fill(ic) = fill(ic) + 1
      end do
   end subroutine cell_list_build

   !> Store one atom in a slot, saturating an unbounded reach
   !>
   !> @param[in,out] self    cell list
   !> @param[in]     s       slot
   !> @param[in]     atom    atom id
   !> @param[in]     center  atom center
   !> @param[in]     reach   atom reach
   !> @param[in]     sat     `sqrt(huge)`, saturation threshold
   pure subroutine store_slot(self, s, atom, center, reach, sat)
      !> Cell list
      type(moist_cell_list_type), intent(inout) :: self
      !> Slot
      integer, intent(in) :: s
      !> Atom id
      integer, intent(in) :: atom
      !> Atom center
      real(wp), intent(in) :: center(3)
      !> Atom reach
      real(wp), intent(in) :: reach
      !> Saturation threshold
      real(wp), intent(in) :: sat

      self%slot_atom(s) = atom
      self%slot_geom(1:3, s) = center
      if (.not. ieee_is_finite(reach) .or. reach >= sat) then
         self%slot_geom(row_reach, s) = huge(0.0_wp)
         self%slot_geom(row_reach_sq, s) = huge(0.0_wp)
      else
         self%slot_geom(row_reach, s) = reach
         self%slot_geom(row_reach_sq, s) = reach*reach
      end if
   end subroutine store_slot

   !* ================================================================================= *!
   !*                                     Queries                                       *!
   !* ================================================================================= *!

   !> Atoms whose reach sphere intersects the box `[lo, hi]`
   !>
   !> - closest-point squared distance against the squared reach, widened by
   !>   `hit_slack`
   !> - corners in either order
   !> - error stop on a non-finite corner
   !>
   !> @param[in]  self  cell list
   !> @param[in]  lo    box corner
   !> @param[in]  hi    opposite box corner
   !> @param[out] out   ascending atom ids in `out(1:nout)`, `size(out) >= natoms`
   !> @param[out] nout  number of atoms found
   pure subroutine cell_list_query_box(self, lo, hi, out, nout)
      !> Cell list
      class(moist_cell_list_type), intent(in) :: self
      !> Box corner
      real(wp), intent(in) :: lo(3)
      !> Opposite box corner
      real(wp), intent(in) :: hi(3)
      !> Atom ids
      integer, intent(out) :: out(:)
      !> Number of atoms found
      integer, intent(out) :: nout

      call collect(self, min(lo, hi), max(lo, hi), [0.0_wp, 0.0_wp, 0.0_wp], 0.0_wp, &
                   .false., out, nout)
   end subroutine cell_list_query_box

   !> Atoms whose reach sphere intersects the ball `(center, r)`
   !>
   !> - `||c_i - center|| <= r + reach_i`, compared squared, widened by `hit_slack`
   !> - `r >= sqrt(huge)` returns every atom
   !> - error stop on a non-finite center or radius
   !>
   !> @param[in]  self    cell list
   !> @param[in]  center  ball center
   !> @param[in]  r       ball radius, `>= 0`
   !> @param[out] out     ascending atom ids in `out(1:nout)`, `size(out) >= natoms`
   !> @param[out] nout    number of atoms found
   pure subroutine cell_list_query_sphere(self, center, r, out, nout)
      !> Cell list
      class(moist_cell_list_type), intent(in) :: self
      !> Ball center
      real(wp), intent(in) :: center(3)
      !> Ball radius
      real(wp), intent(in) :: r
      !> Atom ids
      integer, intent(out) :: out(:)
      !> Number of atoms found
      integer, intent(out) :: nout

      if (r < 0.0_wp) error stop "moist cell list: query_sphere radius must be non-negative"
      call collect(self, center - r, center + r, center, r, .true., out, nout)
   end subroutine cell_list_query_sphere

   !> Atoms whose reach sphere contains a point
   !>
   !> - degenerate box query, widened by `hit_slack`
   !> - error stop on a non-finite point
   !>
   !> @param[in]  self   cell list
   !> @param[in]  point  query point
   !> @param[out] out    ascending atom ids in `out(1:nout)`, `size(out) >= natoms`
   !> @param[out] nout   number of atoms found
   pure subroutine cell_list_query_point(self, point, out, nout)
      !> Cell list
      class(moist_cell_list_type), intent(in) :: self
      !> Query point
      real(wp), intent(in) :: point(3)
      !> Atom ids
      integer, intent(out) :: out(:)
      !> Number of atoms found
      integer, intent(out) :: nout

      call collect(self, point, point, [0.0_wp, 0.0_wp, 0.0_wp], 0.0_wp, .false., out, nout)
   end subroutine cell_list_query_point

   !> Walk the cells overlapped by the expanded box and test each atom
   !>
   !> - box mode: closest-point distance to `[blo, bhi]`
   !> - sphere mode: center distance against `r + reach_i`
   !>
   !> @param[in]  self    cell list
   !> @param[in]  blo     lower box corner
   !> @param[in]  bhi     upper box corner
   !> @param[in]  center  ball center, sphere mode only
   !> @param[in]  r       ball radius, sphere mode only
   !> @param[in]  sphere  sphere mode
   !> @param[out] out     ascending atom ids
   !> @param[out] nout    number of atoms found
   pure subroutine collect(self, blo, bhi, center, r, sphere, out, nout)
      !> Cell list
      type(moist_cell_list_type), intent(in) :: self
      !> Lower box corner
      real(wp), intent(in) :: blo(3)
      !> Upper box corner
      real(wp), intent(in) :: bhi(3)
      !> Ball center
      real(wp), intent(in) :: center(3)
      !> Ball radius
      real(wp), intent(in) :: r
      !> Sphere mode
      logical, intent(in) :: sphere
      !> Atom ids
      integer, intent(out) :: out(:)
      !> Number of atoms found
      integer, intent(out) :: nout

      real(wp) :: tlo(3), thi(3), elo(3), ehi(3), upper(3), margin(3)
      integer :: ilo(3), ihi(3), ix, iy, iz, ic, s, j, k

      nout = 0
      if (self%natoms <= 0) return
      if (size(out) < self%natoms) then
         error stop "moist cell list: output buffer smaller than the atom count"
      end if

      if (.not. (all(ieee_is_finite(blo)) .and. all(ieee_is_finite(bhi)) &
                 .and. ieee_is_finite(r))) then
         error stop "moist cell list: query region must be finite"
      end if
      if (r >= sqrt(huge(0.0_wp))) then
         do j = 1, self%natoms
            out(j) = j
         end do
         nout = self%natoms
         return
      end if

      if (self%full_scan) then
         do s = 1, self%natoms
            if (slot_hit(self, s, blo, bhi, center, r, sphere)) then
               nout = nout + 1
               out(nout) = self%slot_atom(s)
            end if
         end do
         return
      end if

      ! Widened hit radius plus the rounding of the corners themselves
      do k = 1, 3
         margin(k) = self%walk_reach + walk_slack*(r + max(abs(blo(k)), abs(bhi(k))))
      end do
      ! Clamp to one cell beyond the grid before scaling, so far boxes cannot overflow
      upper = self%origin + real([self%nx, self%ny, self%nz], wp)*self%cell_side
      elo = min(max(blo - margin, self%origin - self%cell_side), upper + self%cell_side)
      ehi = min(max(bhi + margin, self%origin - self%cell_side), upper + self%cell_side)
      tlo = (elo - self%origin)*self%inv_cell
      thi = (ehi - self%origin)*self%inv_cell
      if (any(thi < 0.0_wp)) return
      if (tlo(1) >= real(self%nx, wp) .or. tlo(2) >= real(self%ny, wp) &
          .or. tlo(3) >= real(self%nz, wp)) return

      ilo = int(max(0.0_wp, tlo))
      ihi(1) = int(min(real(self%nx - 1, wp), thi(1)))
      ihi(2) = int(min(real(self%ny - 1, wp), thi(2)))
      ihi(3) = int(min(real(self%nz - 1, wp), thi(3)))

      do iz = ilo(3), ihi(3)
         do iy = ilo(2), ihi(2)
            do ix = ilo(1), ihi(1)
               ic = ix + iy*self%nx + iz*self%nx*self%ny + 1
               do s = self%cell_start(ic), self%cell_start(ic + 1) - 1
                  if (slot_hit(self, s, blo, bhi, center, r, sphere)) then
                     nout = nout + 1
                     out(nout) = self%slot_atom(s)
                  end if
               end do
            end do
         end do
      end do

      ! One cell is already ascending; several need a sort for a box-independent order
      if (any(ilo /= ihi)) call heapsort(out(1:nout))
   end subroutine collect

   !> Whether the atom in slot `s` meets the query region
   !>
   !> - saturated reach always hits
   !> - squared reach widened by `hit_slack`
   !>
   !> @param[in] self    cell list
   !> @param[in] s       slot
   !> @param[in] blo     lower box corner
   !> @param[in] bhi     upper box corner
   !> @param[in] center  ball center
   !> @param[in] r       ball radius
   !> @param[in] sphere  sphere mode
   pure function slot_hit(self, s, blo, bhi, center, r, sphere) result(hit)
      !> Cell list
      type(moist_cell_list_type), intent(in) :: self
      !> Slot
      integer, intent(in) :: s
      !> Lower box corner
      real(wp), intent(in) :: blo(3)
      !> Upper box corner
      real(wp), intent(in) :: bhi(3)
      !> Ball center
      real(wp), intent(in) :: center(3)
      !> Ball radius
      real(wp), intent(in) :: r
      !> Sphere mode
      logical, intent(in) :: sphere
      !> Whether the slot's reach sphere meets the region
      logical :: hit

      real(wp) :: dx, dy, dz, reach

      hit = .true.
      if (self%slot_geom(row_reach, s) >= huge(0.0_wp)) return
      if (sphere) then
         reach = self%slot_geom(row_reach, s)
         dx = self%slot_geom(1, s) - center(1)
         dy = self%slot_geom(2, s) - center(2)
         dz = self%slot_geom(3, s) - center(3)
         hit = dx*dx + dy*dy + dz*dz <= (r + reach)*(r + reach)*(1.0_wp + hit_slack)
      else
         hit = box_d2(self%slot_geom(1:3, s), blo, bhi) &
               <= self%slot_geom(row_reach_sq, s)*(1.0_wp + hit_slack)
      end if
   end function slot_hit

   !> Squared distance from a point to an axis-aligned box, zero inside
   !>
   !> @param[in] c   point
   !> @param[in] lo  lower box corner
   !> @param[in] hi  upper box corner
   pure function box_d2(c, lo, hi) result(d2)
      !> Point
      real(wp), intent(in) :: c(3)
      !> Lower box corner
      real(wp), intent(in) :: lo(3)
      !> Upper box corner
      real(wp), intent(in) :: hi(3)
      !> Squared closest distance
      real(wp) :: d2

      real(wp) :: dx, dy, dz

      dx = 0.0_wp
      if (c(1) < lo(1)) then
         dx = lo(1) - c(1)
      else if (c(1) > hi(1)) then
         dx = c(1) - hi(1)
      end if
      dy = 0.0_wp
      if (c(2) < lo(2)) then
         dy = lo(2) - c(2)
      else if (c(2) > hi(2)) then
         dy = c(2) - hi(2)
      end if
      dz = 0.0_wp
      if (c(3) < lo(3)) then
         dz = lo(3) - c(3)
      else if (c(3) > hi(3)) then
         dz = c(3) - hi(3)
      end if
      d2 = dx*dx + dy*dy + dz*dz
   end function box_d2

   !> In-place ascending heapsort, allocation-free
   !>
   !> @param[in,out] a  values to sort
   pure subroutine heapsort(a)
      !> Values to sort
      integer, intent(inout) :: a(:)

      integer :: n, i, tmp

      n = size(a)
      do i = n/2, 1, -1
         call sift_down(a, i, n)
      end do
      do i = n, 2, -1
         tmp = a(1)
         a(1) = a(i)
         a(i) = tmp
         call sift_down(a, 1, i - 1)
      end do
   end subroutine heapsort

   !> Restore the max-heap property below `root` within `a(1:last)`
   !>
   !> @param[in,out] a     heap
   !> @param[in]     root  root index
   !> @param[in]     last  last heap index
   pure subroutine sift_down(a, root, last)
      !> Heap
      integer, intent(inout) :: a(:)
      !> Root index
      integer, intent(in) :: root
      !> Last heap index
      integer, intent(in) :: last

      integer :: i, child, tmp

      i = root
      do
         child = 2*i
         if (child > last) exit
         if (child < last) then
            if (a(child + 1) > a(child)) child = child + 1
         end if
         if (a(i) >= a(child)) exit
         tmp = a(i)
         a(i) = a(child)
         a(child) = tmp
         i = child
      end do
   end subroutine sift_down

   !* ================================================================================= *!
   !*                                    Lifecycle                                      *!
   !* ================================================================================= *!

   !> Deallocate all storage and reset to an empty list
   !>
   !> @param[in,out] self  cell list
   subroutine cell_list_destroy(self)
      !> Cell list
      class(moist_cell_list_type), intent(inout) :: self

      self%natoms = 0
      self%full_scan = .false.
      self%cell_side = 0.0_wp
      self%inv_cell = 0.0_wp
      self%walk_reach = 0.0_wp
      self%origin = 0.0_wp
      self%nx = 0
      self%ny = 0
      self%nz = 0
      self%ncells = 0
      if (allocated(self%cell_start)) deallocate (self%cell_start)
      if (allocated(self%slot_atom)) deallocate (self%slot_atom)
      if (allocated(self%slot_geom)) deallocate (self%slot_geom)
   end subroutine cell_list_destroy

end module moist_math_cell_list
