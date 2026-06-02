!> Uniform cubic cell grid with precomputed per-cell candidate atom lists.
!>
!> Given atom positions and a per-atom effective reach `r_eff(j)`, builds a
!> spatial index whose cells hold every atom whose sphere `(center_j, r_eff(j))`
!> intersects the cell AABB. Query at an arbitrary point P is O(1): compute the
!> clamped cell index from origin + inverse cell size, then read the cell's
!> precomputed list.
!>
!> Correctness invariant: for every point P inside cell C and every atom j
!> with `||P - center_j|| <= r_eff(j)`, atom j appears in C's candidate list.
!> During build, each atom's sphere is tested against every cell AABB in its
!> bounding range, so this holds for any positive cell side length.
!>
!> By default `cell_side = maxval(r_eff)`, keeping the build cheap (each atom
!> reaches at most 2^3 cells). The optional `cell_fraction` parameter
!> (0 < cell_fraction <= 1) scales the cell side to
!> `maxval(r_eff) * cell_fraction`, producing finer spatial bins at the cost
!> of increased memory: each atom may appear in up to O(1/cell_fraction^3)
!> cells. Finer bins improve screening (fewer candidates per query) which
!> pays off when the downstream per-candidate cost (SSD, LSF) is high.
!>
!> Points outside the atom bounding box are strictly clamped to the nearest
!> boundary cell. The consumer guarantees correctness only for points inside
!> the bounding box; external points fall back onto the clamped cell's list,
!> which is sufficient when the downstream screening weight has already decayed
!> below threshold.
!>
!> Structural parallel to `moist_math_adjacency_list` but intentionally not
!> factored: the query model here is asymmetric (arbitrary point -> atoms), so
!> sharing the build code would obscure semantics.
!>
!> For small systems, the per-cell fan-out degenerates: each atom reaches
!> most cells anyway, so we pay build cost for no query-time benefit. The
!> optional `full_scan_below` argument to `build` short-circuits to a
!> single-cell grid containing every atom - every query then returns the
!> full atom list with O(1) overhead.
!>
!> Usage:
!> ```fortran
!> type(moist_cell_grid_type) :: grid
!> call grid%build(xyz, r_eff)                       ! spatially binned
!> call grid%build(xyz, r_eff, full_scan_below=50)   ! auto full scan if small
!> call grid%build(xyz, r_eff, cell_fraction=0.5_wp)  ! finer cells
!> call grid%query(point, start, n)
!> ! candidates are grid%cell_nlat(start+1 : start+n)
!> call grid%destroy()
!> ```
module moist_math_cell_grid
   use mctc_env, only: wp
   implicit none
   private

   public :: moist_cell_grid_type

   !> Uniform cell grid with per-cell candidate atom lists in CSR format.
   type :: moist_cell_grid_type
      !> Cell side length (equals maxval(r_eff) after build)
      real(wp) :: cell_side = 0.0_wp
      !> Inverse cell side, cached for hot-path arithmetic
      real(wp) :: inv_cell = 0.0_wp
      !> Grid origin: lower corner of the atom bounding box
      real(wp) :: origin(3) = 0.0_wp
      !> Grid dimensions along each axis
      integer :: nx = 0, ny = 0, nz = 0
      !> Total number of cells: nx*ny*nz
      integer :: ncells = 0
      !> Number of atoms indexed
      integer :: natoms = 0
      !> True when build took the full-scan shortcut (single cell, all atoms)
      logical :: full_scan = .false.
      !> Cell fraction: cell_side = maxval(r_eff) * cell_fraction (default 1.0)
      real(wp) :: cell_fraction = 1.0_wp
      !> CSR offset for each cell (ncells), 0-based
      integer, allocatable :: cell_inl(:)
      !> Number of candidate atoms per cell (ncells)
      integer, allocatable :: cell_nnl(:)
      !> Flat-packed candidate atom indices (sum(cell_nnl))
      integer, allocatable :: cell_nlat(:)
   contains
      !> Build the grid from atom positions and effective reach
      procedure :: build => moist_cell_grid_build
      !> O(1) candidate lookup for an evaluation point
      procedure :: query => moist_cell_grid_query
      !> Deallocate all storage
      procedure :: destroy => moist_cell_grid_destroy
      !> Finalizer
      final :: moist_cell_grid_finalize
   end type moist_cell_grid_type

contains

   !> Build the per-cell candidate lists.
   !>
   !> Cell side is set to `maxval(r_eff) * cell_fraction`. Each atom is
   !> enumerated over all cells whose AABB its sphere intersects using a
   !> sphere-AABB closest-point test. Build cost scales with the total number
   !> of (atom, cell) overlaps and uses the standard two-pass CSR fill
   !> (count, prefix-sum, fill), reusing cell_nnl as a running offset during
   !> the fill pass as in moist_math_adjacency_list.
   !>
   !> With `cell_fraction = 1.0` (default), each atom overlaps at most 2^3
   !> cells. Smaller fractions produce finer bins - each atom overlaps up to
   !> O(1/cell_fraction^3) cells - but each query returns fewer candidates.
   !>
   !> If `full_scan_below` is supplied and `natoms < full_scan_below`, the
   !> grid collapses to a single cell holding every atom. This skips the
   !> fan-out work for small systems where the grid would degenerate anyway
   !> (r_eff comparable to the molecular diameter).
   !>
   !> @param[inout] self             Grid instance (reset on entry)
   !> @param[in]    xyz              Atom centres, shape (3, natoms)
   !> @param[in]    r_eff            Effective per-atom reach (natoms)
   !> @param[in]    full_scan_below  Optional: natoms below this force full scan
   !> @param[in]    cell_fraction    Optional: fraction of maxval(r_eff) for cell side (default 1.0)
   subroutine moist_cell_grid_build(self, xyz, r_eff, full_scan_below, cell_fraction)
      class(moist_cell_grid_type), intent(inout) :: self
      !> Atom centres
      real(wp), intent(in) :: xyz(:, :)
      !> Effective per-atom reach (radius + screening shell)
      real(wp), intent(in) :: r_eff(:)
      !> Optional natoms threshold: below this, use a single full-scan cell
      integer, intent(in), optional :: full_scan_below
      !> Optional cell side fraction (0 < cell_fraction <= 1, default 1.0)
      real(wp), intent(in), optional :: cell_fraction

      real(wp) :: xmin(3), xmax(3)
      real(wp) :: r_max, cell_lo(3), cell_hi(3), center(3), rj, rj2, d2
      integer :: j, ix, iy, iz, ic, offset, total
      integer :: ix_lo, ix_hi, iy_lo, iy_hi, iz_lo, iz_hi

      call self%destroy()

      self%natoms = size(xyz, 2)
      if (self%natoms <= 0) then
         allocate (self%cell_nlat(0))
         return
      end if

      ! Full-scan shortcut for small systems. The grid would otherwise pay
      ! sphere-AABB build cost only to return "nearly every atom" per cell.
      if (present(full_scan_below)) then
         if (self%natoms < full_scan_below) then
            call build_full_scan(self, xyz)
            return
         end if
      end if

      r_max = maxval(r_eff)
      if (r_max <= 0.0_wp) then
         allocate (self%cell_nlat(0))
         return
      end if

      if (present(cell_fraction)) then
         self%cell_fraction = max(tiny(1.0_wp), min(1.0_wp, cell_fraction))
      else
         self%cell_fraction = 1.0_wp
      end if

      self%cell_side = r_max*self%cell_fraction
      self%inv_cell = 1.0_wp/self%cell_side

      ! Bounding box over atom centres
      xmin = xyz(:, 1)
      xmax = xyz(:, 1)
      do j = 2, self%natoms
         if (xyz(1, j) < xmin(1)) xmin(1) = xyz(1, j)
         if (xyz(1, j) > xmax(1)) xmax(1) = xyz(1, j)
         if (xyz(2, j) < xmin(2)) xmin(2) = xyz(2, j)
         if (xyz(2, j) > xmax(2)) xmax(2) = xyz(2, j)
         if (xyz(3, j) < xmin(3)) xmin(3) = xyz(3, j)
         if (xyz(3, j) > xmax(3)) xmax(3) = xyz(3, j)
      end do
      self%origin = xmin

      self%nx = max(1, floor((xmax(1) - xmin(1))*self%inv_cell) + 1)
      self%ny = max(1, floor((xmax(2) - xmin(2))*self%inv_cell) + 1)
      self%nz = max(1, floor((xmax(3) - xmin(3))*self%inv_cell) + 1)
      self%ncells = self%nx*self%ny*self%nz

      allocate (self%cell_inl(self%ncells), source=0)
      allocate (self%cell_nnl(self%ncells), source=0)

      ! Pass A: counts
      do j = 1, self%natoms
         center = xyz(:, j)
         rj = r_eff(j)
         rj2 = rj*rj

         ix_lo = max(0, floor((center(1) - rj - self%origin(1))*self%inv_cell))
         ix_hi = min(self%nx - 1, floor((center(1) + rj - self%origin(1))*self%inv_cell))
         iy_lo = max(0, floor((center(2) - rj - self%origin(2))*self%inv_cell))
         iy_hi = min(self%ny - 1, floor((center(2) + rj - self%origin(2))*self%inv_cell))
         iz_lo = max(0, floor((center(3) - rj - self%origin(3))*self%inv_cell))
         iz_hi = min(self%nz - 1, floor((center(3) + rj - self%origin(3))*self%inv_cell))

         do iz = iz_lo, iz_hi
            cell_lo(3) = self%origin(3) + iz*self%cell_side
            cell_hi(3) = cell_lo(3) + self%cell_side
            do iy = iy_lo, iy_hi
               cell_lo(2) = self%origin(2) + iy*self%cell_side
               cell_hi(2) = cell_lo(2) + self%cell_side
               do ix = ix_lo, ix_hi
                  cell_lo(1) = self%origin(1) + ix*self%cell_side
                  cell_hi(1) = cell_lo(1) + self%cell_side

                  d2 = sphere_aabb_closest_d2(center, cell_lo, cell_hi)
                  if (d2 <= rj2) then
                     ic = ix + iy*self%nx + iz*self%nx*self%ny + 1
                     self%cell_nnl(ic) = self%cell_nnl(ic) + 1
                  end if
               end do
            end do
         end do
      end do

      ! Prefix sum: counts -> offsets
      offset = 0
      do ic = 1, self%ncells
         self%cell_inl(ic) = offset
         offset = offset + self%cell_nnl(ic)
      end do
      total = offset
      allocate (self%cell_nlat(total))

      ! Reset cell_nnl as running counter for the fill pass
      self%cell_nnl = 0

      ! Pass B: fill
      do j = 1, self%natoms
         center = xyz(:, j)
         rj = r_eff(j)
         rj2 = rj*rj

         ix_lo = max(0, floor((center(1) - rj - self%origin(1))*self%inv_cell))
         ix_hi = min(self%nx - 1, floor((center(1) + rj - self%origin(1))*self%inv_cell))
         iy_lo = max(0, floor((center(2) - rj - self%origin(2))*self%inv_cell))
         iy_hi = min(self%ny - 1, floor((center(2) + rj - self%origin(2))*self%inv_cell))
         iz_lo = max(0, floor((center(3) - rj - self%origin(3))*self%inv_cell))
         iz_hi = min(self%nz - 1, floor((center(3) + rj - self%origin(3))*self%inv_cell))

         do iz = iz_lo, iz_hi
            cell_lo(3) = self%origin(3) + iz*self%cell_side
            cell_hi(3) = cell_lo(3) + self%cell_side
            do iy = iy_lo, iy_hi
               cell_lo(2) = self%origin(2) + iy*self%cell_side
               cell_hi(2) = cell_lo(2) + self%cell_side
               do ix = ix_lo, ix_hi
                  cell_lo(1) = self%origin(1) + ix*self%cell_side
                  cell_hi(1) = cell_lo(1) + self%cell_side

                  d2 = sphere_aabb_closest_d2(center, cell_lo, cell_hi)
                  if (d2 <= rj2) then
                     ic = ix + iy*self%nx + iz*self%nx*self%ny + 1
                     self%cell_nnl(ic) = self%cell_nnl(ic) + 1
                     self%cell_nlat(self%cell_inl(ic) + self%cell_nnl(ic)) = j
                  end if
               end do
            end do
         end do
      end do

   end subroutine moist_cell_grid_build

   !> Degenerate build: a single cell holding every atom.
   !>
   !> Invoked when the caller deems the system too small for spatial
   !> binning to pay off. The grid is set up so `query` naturally returns
   !> the full atom list - nx = ny = nz = 1 clamps every point to cell 1.
   !>
   !> @param[inout] self  Grid instance (natoms already set by caller)
   !> @param[in]    xyz   Atom centres; only used to seed `origin`
   subroutine build_full_scan(self, xyz)
      !> Grid instance
      type(moist_cell_grid_type), intent(inout) :: self
      !> Atom centres
      real(wp), intent(in) :: xyz(:, :)

      integer :: j

      ! cell_side/inv_cell just have to be finite and positive; the clamp
      ! in query forces every point to cell 1 regardless of these values.
      self%cell_side = 1.0_wp
      self%inv_cell = 1.0_wp
      self%origin = minval(xyz, dim=2)
      self%nx = 1
      self%ny = 1
      self%nz = 1
      self%ncells = 1
      self%full_scan = .true.

      allocate (self%cell_inl(1), source=0)
      allocate (self%cell_nnl(1))
      self%cell_nnl(1) = self%natoms

      allocate (self%cell_nlat(self%natoms))
      do j = 1, self%natoms
         self%cell_nlat(j) = j
      end do
   end subroutine build_full_scan

   !> Locate the candidate atom list for an evaluation point.
   !>
   !> Zero-copy: the caller consumes `self%cell_nlat(start+1 : start+n)`. Cell
   !> index is strictly clamped to `[0, n{x,y,z}-1]`; no fallback branch.
   !>
   !> @param[in]  self   Grid instance (must be built)
   !> @param[in]  point  Evaluation point (3)
   !> @param[out] start  CSR offset; first candidate is at cell_nlat(start+1)
   !> @param[out] n      Number of candidates
   pure subroutine moist_cell_grid_query(self, point, start, n)
      class(moist_cell_grid_type), intent(in) :: self
      !> Evaluation point
      real(wp), intent(in) :: point(3)
      !> CSR offset into cell_nlat
      integer, intent(out) :: start
      !> Number of candidate atoms
      integer, intent(out) :: n

      integer :: ix, iy, iz, ic

      if (self%ncells <= 0) then
         start = 0
         n = 0
         return
      end if

      ix = max(0, min(self%nx - 1, floor((point(1) - self%origin(1))*self%inv_cell)))
      iy = max(0, min(self%ny - 1, floor((point(2) - self%origin(2))*self%inv_cell)))
      iz = max(0, min(self%nz - 1, floor((point(3) - self%origin(3))*self%inv_cell)))

      ic = ix + iy*self%nx + iz*self%nx*self%ny + 1
      start = self%cell_inl(ic)
      n = self%cell_nnl(ic)
   end subroutine moist_cell_grid_query

   !> Deallocate all storage.
   subroutine moist_cell_grid_destroy(self)
      class(moist_cell_grid_type), intent(inout) :: self

      self%cell_side = 0.0_wp
      self%inv_cell = 0.0_wp
      self%origin = 0.0_wp
      self%nx = 0
      self%ny = 0
      self%nz = 0
      self%ncells = 0
      self%natoms = 0
      self%full_scan = .false.
      self%cell_fraction = 1.0_wp
      if (allocated(self%cell_inl)) deallocate (self%cell_inl)
      if (allocated(self%cell_nnl)) deallocate (self%cell_nnl)
      if (allocated(self%cell_nlat)) deallocate (self%cell_nlat)
   end subroutine moist_cell_grid_destroy

   !> Finalizer delegates to destroy.
   subroutine moist_cell_grid_finalize(self)
      type(moist_cell_grid_type), intent(inout) :: self
      call self%destroy()
   end subroutine moist_cell_grid_finalize

   !> Squared distance from a point to an axis-aligned box; 0 if inside.
   pure function sphere_aabb_closest_d2(center, lo, hi) result(d2)
      !> Sphere centre
      real(wp), intent(in) :: center(3)
      !> AABB lower corner
      real(wp), intent(in) :: lo(3)
      !> AABB upper corner
      real(wp), intent(in) :: hi(3)
      !> Squared closest distance
      real(wp) :: d2

      real(wp) :: dx, dy, dz

      dx = 0.0_wp
      if (center(1) < lo(1)) then
         dx = lo(1) - center(1)
      else if (center(1) > hi(1)) then
         dx = center(1) - hi(1)
      end if

      dy = 0.0_wp
      if (center(2) < lo(2)) then
         dy = lo(2) - center(2)
      else if (center(2) > hi(2)) then
         dy = center(2) - hi(2)
      end if

      dz = 0.0_wp
      if (center(3) < lo(3)) then
         dz = lo(3) - center(3)
      else if (center(3) > hi(3)) then
         dz = center(3) - hi(3)
      end if

      d2 = dx*dx + dy*dy + dz*dz
   end function sphere_aabb_closest_d2

end module moist_math_cell_grid
