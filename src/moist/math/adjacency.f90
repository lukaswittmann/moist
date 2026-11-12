!> Sparse adjacency list (neighbour list) in compressed sparse row (CSR) format
!>
!> Reusable spatial neighbour list, built from any set of 3D coordinates and
!> a global interaction cutoff
!>
!> - a uniform cell grid gives O(N) build cost
!> - points are counting-sorted into their cells, with the coordinates of each cell
!>   contiguous, so the inner loop streams through memory
!> - count and fill passes run in OpenMP over cells; the result is independent of
!>   the thread count
!> - the result is stored in CSR format for cache-friendly traversal
!>
!> Usage:
!> ```fortran
!> type(adjacency_list_type) :: nlist
!> call nlist%init(cutoff=2.0_wp, sorted=.true.)  ! sorted is optional (default .false.)
!> call nlist%update(xyz, nthreads)    ! (re)build from coordinates
!> ids = nlist%get_neighbours(i)       ! integer array of neighbour indices
!> call nlist%destroy()
!> ```
module moist_math_adjacency
   use mctc_env, only: wp, error_type
   use moist_math_sorter_quicksort, only: qsort
!$ use omp_lib, only: omp_get_max_threads
   implicit none(type, external)
   private

   public :: adjacency_list_type

   !> Below this many points the build runs on one thread
   integer, parameter :: parallel_min_points = 512

   !> Adjacency list in compressed sparse row (CSR) format
   !>
   !> For point i the neighbours are stored at
   !>   nlat( inl(i)+1 : inl(i)+nnl(i) )
   !> with corresponding center-center distances in
   !>   dist( inl(i)+1 : inl(i)+nnl(i) )
   !>
   !> When sorted=.true., both arrays are ordered by ascending distance
   type :: adjacency_list_type
      !> Global interaction cutoff distance
      real(wp) :: cutoff = 0.0_wp
      !> Whether neighbours are sorted by ascending distance
      logical :: sorted = .false.
      !> Number of points in the list
      integer :: n = 0
      !> Offset into nlat for each point (n)
      integer, allocatable :: inl(:)
      !> Number of neighbours for each point (n)
      integer, allocatable :: nnl(:)
      !> Flat-packed neighbour indices (sum(nnl)); sorted by distance when sorted=.true
      integer, allocatable :: nlat(:)
      !> Center-center distances parallel to nlat; sorted ascending when sorted=.true
      real(wp), allocatable :: dist(:)
   contains
      !> Set the interaction cutoff
      procedure :: init => adjacency_list_init
      !> (Re)build the list from a coordinate array
      procedure :: update => adjacency_list_update
      !> Return the neighbour indices for point i as an array
      procedure :: get_neighbours => adjacency_list_get_neighbours
      !> Deallocate all storage
      procedure :: destroy => adjacency_list_destroy
      !> Finalizer
      final :: adjacency_list_finalize
   end type adjacency_list_type

contains

   !> Set the interaction cutoff, before the first update
   !>
   !> Every setting is reset, so init fully defines the configuration of the list:
   !> omitting `sorted` means unsorted, regardless of how the instance was used before
   !>
   !> @param[inout] self    Adjacency list instance
   !> @param[in]    cutoff  Global interaction cutoff distance
   !> @param[in]    sorted  Sort neighbours by ascending distance (default .false.)
   subroutine adjacency_list_init(self, cutoff, sorted)
      class(adjacency_list_type), intent(inout) :: self
      !> Global interaction cutoff distance
      real(wp), intent(in) :: cutoff
      !> Whether to sort neighbours by ascending distance
      logical, intent(in), optional :: sorted

      call self%destroy()
      self%cutoff = cutoff
      self%sorted = .false.
      if (present(sorted)) self%sorted = sorted
   end subroutine adjacency_list_init

   !> (Re)build the neighbour list from a coordinate array using a cell grid
   !>
   !> - the coordinates are partitioned into a uniform cubic grid of cell side
   !>   length equal to the cutoff
   !> - points are counting-sorted into their cells, descending point index within
   !>   a cell; each cell's coordinates are contiguous
   !> - only the 27 cells around a point are inspected for neighbours, giving
   !>   O(N*k) total cost for an average neighbour count k
   !> - self-pairs (i==i) are excluded
   !> - neighbour order and distances do not depend on the thread count
   !>
   !> @param[inout] self      Adjacency list instance (cutoff must be set)
   !> @param[in]    xyz       Coordinate array (3, npoints)
   !> @param[in]    nthreads  OpenMP team size (optional, default omp_get_max_threads)
   subroutine adjacency_list_update(self, xyz, nthreads)
      class(adjacency_list_type), intent(inout) :: self
      real(wp), intent(in) :: xyz(:, :)
      !> OpenMP team size
      integer, intent(in), optional :: nthreads

      integer :: npoints, nt, i, p, q, s, ic, jc, k, cnt, nnbr
      integer :: nx, ny, nz, ncells
      integer :: nbr(27)
      real(wp) :: cutoff2, cell_inv, lo(3), hi(3)
      real(wp) :: xi, yi, zi, dxij, dyij, dzij, d2
      logical :: mismatch

      !> Cell of each point
      integer, allocatable :: cell_id(:)
      !> Offset of each cell into the cell-ordered arrays (ncells + 1)
      integer, allocatable :: start(:)
      !> Next free slot of each cell during the scatter
      integer, allocatable :: cursor(:)
      !> Point index of each cell-ordered slot
      integer, allocatable :: order(:)
      !> Cell-ordered coordinates, one contiguous stream per axis
      real(wp), allocatable :: px(:), py(:), pz(:)
      !> CSR arrays, built locally and moved into self after the parallel regions
      integer, allocatable :: inl(:), nnl(:), nlat(:)
      real(wp), allocatable :: dist(:)

      npoints = size(xyz, 2)
      cutoff2 = self%cutoff*self%cutoff

      nt = 1
!$    nt = max(1, omp_get_max_threads())
      if (present(nthreads)) then
         if (nthreads > 0) nt = nthreads
      end if
      if (npoints < parallel_min_points) nt = 1

      call self%destroy()
      self%n = npoints
      allocate (inl(npoints), source=0)
      allocate (nnl(npoints), source=0)

      ! Handle trivial cases
      if (npoints <= 1 .or. self%cutoff <= 0.0_wp) then
         allocate (nlat(0), dist(0))
         call move_csr(self, inl, nnl, nlat, dist)
         return
      end if

      ! Uniform cell grid over the bounding box
      lo = minval(xyz, dim=2)
      hi = maxval(xyz, dim=2)
      cell_inv = 1.0_wp/self%cutoff
      nx = max(1, floor((hi(1) - lo(1))*cell_inv) + 1)
      ny = max(1, floor((hi(2) - lo(2))*cell_inv) + 1)
      nz = max(1, floor((hi(3) - lo(3))*cell_inv) + 1)
      ncells = nx*ny*nz

      ! Counting sort of the points into their cells
      allocate (cell_id(npoints))
      allocate (start(ncells + 1), source=0)
      do i = 1, npoints
         ic = min(int((xyz(1, i) - lo(1))*cell_inv), nx - 1) &
              + min(int((xyz(2, i) - lo(2))*cell_inv), ny - 1)*nx &
              + min(int((xyz(3, i) - lo(3))*cell_inv), nz - 1)*nx*ny + 1
         cell_id(i) = ic
         start(ic + 1) = start(ic + 1) + 1
      end do
      start(1) = 1
      do ic = 1, ncells
         start(ic + 1) = start(ic) + start(ic + 1)
      end do

      ! Descending point index within a cell
      allocate (cursor(ncells), source=start(1:ncells))
      allocate (order(npoints))
      do i = npoints, 1, -1
         ic = cell_id(i)
         order(cursor(ic)) = i
         cursor(ic) = cursor(ic) + 1
      end do
      deallocate (cursor, cell_id)

      allocate (px(npoints), py(npoints), pz(npoints))
      do p = 1, npoints
         px(p) = xyz(1, order(p))
         py(p) = xyz(2, order(p))
         pz(p) = xyz(3, order(p))
      end do

      ! Count pass: nnl(i) for each point
      !$omp parallel do num_threads(nt) default(none) schedule(dynamic, 16) &
      !$omp& shared(ncells, nx, ny, nz, start, order, px, py, pz, cutoff2, nnl) &
      !$omp& private(ic, nbr, nnbr, p, q, s, jc, cnt, xi, yi, zi, dxij, dyij, dzij, d2)
      do ic = 1, ncells
         if (start(ic + 1) == start(ic)) cycle
         call cell_stencil(ic, nx, ny, nz, nbr, nnbr)
         do p = start(ic), start(ic + 1) - 1
            xi = px(p)
            yi = py(p)
            zi = pz(p)
            cnt = 0
            do s = 1, nnbr
               jc = nbr(s)
               do q = start(jc), start(jc + 1) - 1
                  if (q == p) cycle
                  dxij = xi - px(q)
                  dyij = yi - py(q)
                  dzij = zi - pz(q)
                  d2 = dxij*dxij + dyij*dyij + dzij*dzij
                  if (d2 <= cutoff2) cnt = cnt + 1
               end do
            end do
            nnl(order(p)) = cnt
         end do
      end do
      !$omp end parallel do

      ! Prefix sum: inl from nnl
      k = 0
      do i = 1, npoints
         inl(i) = k
         k = k + nnl(i)
      end do
      allocate (nlat(k), dist(k))

      ! Fill pass: neighbour indices and distances
      ! The capacity guard keeps a count/fill disagreement inside the row
      mismatch = .false.
      !$omp parallel do num_threads(nt) default(none) schedule(dynamic, 16) &
      !$omp& shared(ncells, nx, ny, nz, start, order, px, py, pz, cutoff2, inl, nnl, nlat, dist) &
      !$omp& private(ic, nbr, nnbr, p, q, s, jc, i, k, cnt, xi, yi, zi, dxij, dyij, dzij, d2) &
      !$omp& reduction(.or.:mismatch)
      do ic = 1, ncells
         if (start(ic + 1) == start(ic)) cycle
         call cell_stencil(ic, nx, ny, nz, nbr, nnbr)
         do p = start(ic), start(ic + 1) - 1
            i = order(p)
            xi = px(p)
            yi = py(p)
            zi = pz(p)
            cnt = nnl(i)
            k = 0
            do s = 1, nnbr
               jc = nbr(s)
               do q = start(jc), start(jc + 1) - 1
                  if (q == p) cycle
                  dxij = xi - px(q)
                  dyij = yi - py(q)
                  dzij = zi - pz(q)
                  d2 = dxij*dxij + dyij*dyij + dzij*dzij
                  if (d2 <= cutoff2) then
                     k = k + 1
                     if (k <= cnt) then
                        nlat(inl(i) + k) = order(q)
                        dist(inl(i) + k) = sqrt(d2)
                     end if
                  end if
               end do
            end do
            if (k /= cnt) mismatch = .true.
         end do
      end do
      !$omp end parallel do
      if (mismatch) error stop "adjacency_list: count and fill passes disagree"

      ! Optionally sort each point's neighbours by ascending distance
      if (self%sorted) then
         mismatch = .false.
         !$omp parallel do num_threads(nt) default(none) schedule(dynamic, 64) &
         !$omp& shared(npoints, inl, nnl, nlat, dist) private(i) reduction(.or.:mismatch)
         do i = 1, npoints
            if (nnl(i) > 1) then
               if (.not. sort_row(dist(inl(i) + 1:inl(i) + nnl(i)), nlat(inl(i) + 1:inl(i) + nnl(i)))) &
                  mismatch = .true.
            end if
         end do
         !$omp end parallel do
         if (mismatch) error stop "adjacency_list: dist/nlat slice bounds diverged"
      end if

      call move_csr(self, inl, nnl, nlat, dist)

   end subroutine adjacency_list_update

   !> Neighbouring cells of a cell in the 27-cell stencil, z slowest, x fastest
   !>
   !> @param[in]  ic    Cell index
   !> @param[in]  nx    Cells along x
   !> @param[in]  ny    Cells along y
   !> @param[in]  nz    Cells along z
   !> @param[out] nbr   Indices of the cells inside the grid
   !> @param[out] nnbr  Number of entries in nbr
   pure subroutine cell_stencil(ic, nx, ny, nz, nbr, nnbr)
      !> Cell index
      integer, intent(in) :: ic
      !> Cells along x, y and z
      integer, intent(in) :: nx, ny, nz
      !> Indices of the cells inside the grid
      integer, intent(out) :: nbr(27)
      !> Number of entries in nbr
      integer, intent(out) :: nnbr

      integer :: cx, cy, cz, dx, dy, dz, jx, jy, jz

      cx = mod(ic - 1, nx)
      cy = mod((ic - 1)/nx, ny)
      cz = (ic - 1)/(nx*ny)
      nnbr = 0
      do dz = -1, 1
         jz = cz + dz
         if (jz < 0 .or. jz >= nz) cycle
         do dy = -1, 1
            jy = cy + dy
            if (jy < 0 .or. jy >= ny) cycle
            do dx = -1, 1
               jx = cx + dx
               if (jx < 0 .or. jx >= nx) cycle
               nnbr = nnbr + 1
               nbr(nnbr) = jx + jy*nx + jz*nx*ny + 1
            end do
         end do
      end do
   end subroutine cell_stencil

   !> Sort one row by ascending distance, carrying the neighbour indices along
   !>
   !> @param[inout] dist  Row distances
   !> @param[inout] nlat  Row neighbour indices
   !> @return       ok    False if qsort rejected the companion array
   function sort_row(dist, nlat) result(ok)
      !> Row distances
      real(wp), intent(inout) :: dist(:)
      !> Row neighbour indices
      integer, intent(inout) :: nlat(:)
      !> False if qsort rejected the companion array
      logical :: ok

      type(error_type), allocatable :: sort_error

      call qsort(dist, nlat, sort_error)
      ok = .not. allocated(sort_error)
   end function sort_row

   !> Move locally built CSR arrays into the list
   !>
   !> @param[inout] self  Adjacency list instance
   !> @param[inout] inl   Row offsets, deallocated on return
   !> @param[inout] nnl   Row lengths, deallocated on return
   !> @param[inout] nlat  Neighbour indices, deallocated on return
   !> @param[inout] dist  Neighbour distances, deallocated on return
   subroutine move_csr(self, inl, nnl, nlat, dist)
      class(adjacency_list_type), intent(inout) :: self
      !> Row offsets
      integer, allocatable, intent(inout) :: inl(:)
      !> Row lengths
      integer, allocatable, intent(inout) :: nnl(:)
      !> Neighbour indices
      integer, allocatable, intent(inout) :: nlat(:)
      !> Neighbour distances
      real(wp), allocatable, intent(inout) :: dist(:)

      call move_alloc(inl, self%inl)
      call move_alloc(nnl, self%nnl)
      call move_alloc(nlat, self%nlat)
      call move_alloc(dist, self%dist)
   end subroutine move_csr

   !> Return the neighbour indices for point i as an integer array
   !>
   !> @param[in] self  Adjacency list instance
   !> @param[in] i     Query point index
   !> @return    ids   Array of neighbour indices (length nnl(i))
   pure function adjacency_list_get_neighbours(self, i) result(ids)
      class(adjacency_list_type), intent(in) :: self
      integer, intent(in) :: i
      integer, allocatable :: ids(:)

      ids = self%nlat(self%inl(i) + 1:self%inl(i) + self%nnl(i))
   end function adjacency_list_get_neighbours

   !> Deallocate all storage
   !>
   !> Keeps the configuration (`cutoff`, `sorted`), so `update` alone rebuilds
   !> the list with the same settings; only `init` changes them
   !>
   !> @param[inout] self  Adjacency list instance
   subroutine adjacency_list_destroy(self)
      class(adjacency_list_type), intent(inout) :: self

      self%n = 0
      if (allocated(self%inl)) deallocate (self%inl)
      if (allocated(self%nnl)) deallocate (self%nnl)
      if (allocated(self%nlat)) deallocate (self%nlat)
      if (allocated(self%dist)) deallocate (self%dist)
   end subroutine adjacency_list_destroy

   !> Finalizer - delegates to destroy
   subroutine adjacency_list_finalize(self)
      type(adjacency_list_type), intent(inout) :: self
      call self%destroy()
   end subroutine adjacency_list_finalize

end module moist_math_adjacency
