!> Sparse adjacency list (neighbour list) in compressed sparse row (CSR) format.
!>
!> Provides a reusable spatial neighbour list that can be built from any set of
!> 3D coordinates and a global interaction cutoff. Internally uses a uniform
!> cell grid (linked-list variant) for O(N) build cost. The resulting list is
!> stored in CSR format for cache-friendly traversal.
!>
!> Usage:
!> ```fortran
!> type(adjacency_list_type) :: nlist
!> call nlist%init(cutoff=2.0_wp, sorted=.true.)  ! sorted is optional (default .false.)
!> call nlist%update(xyz)              ! (re)build from coordinates
!> ids = nlist%get_neighbours(i)       ! integer array of neighbour indices
!> call nlist%destroy()
!> ```
module moist_math_adjacency_list
   use mctc_env, only: wp, error_type
   use moist_math_sorter_quicksort, only: qsort
   implicit none
   private

   public :: adjacency_list_type

   !> Adjacency list in compressed sparse row (CSR) format.
   !>
   !> For point i the neighbours are stored at
   !>   nlat( inl(i)+1 : inl(i)+nnl(i) )
   !> with corresponding center-center distances in
   !>   dist( inl(i)+1 : inl(i)+nnl(i) )
   !> When sorted=.true., both arrays are ordered by ascending distance.
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
      !> Flat-packed neighbour indices (sum(nnl)); sorted by distance when sorted=.true.
      integer, allocatable :: nlat(:)
      !> Center-center distances parallel to nlat; sorted ascending when sorted=.true.
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

   !> Set the interaction cutoff. Must be called before the first update.
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
      if (present(sorted)) self%sorted = sorted
   end subroutine adjacency_list_init

   !> (Re)build the neighbour list from a coordinate array using a cell grid.
   !>
   !> The coordinates are partitioned into a uniform cubic grid with cell side
   !> length equal to the cutoff. For each point, only the 27 surrounding cells
   !> are inspected for potential neighbours, giving O(N*k) total cost where k
   !> is the average neighbour count. Self-pairs (i==i) are excluded.
   !>
   !> @param[inout] self  Adjacency list instance (cutoff must be set)
   !> @param[in]    xyz   Coordinate array (3, npoints)
   subroutine adjacency_list_update(self, xyz)
      class(adjacency_list_type), intent(inout) :: self
      real(wp), intent(in) :: xyz(:, :)

      integer :: npoints, i, j, ic, jc
      integer :: cx, cy, cz, dx, dy, dz
      integer :: nx, ny, nz, ncells
      integer :: cx_j, cy_j, cz_j
      real(wp) :: cutoff2, cell_inv
      real(wp) :: xmin, ymin, zmin, xmax, ymax, zmax
      real(wp) :: xi, yi, zi, dxij, dyij, dzij, d2
      integer :: img, capacity

      !> Cell grid: head(cell_id) = first point in cell, next(point) = next point in same cell
      integer, allocatable :: head(:), next(:)
      !> Cell index for each point
      integer, allocatable :: cell_id(:)

      npoints = size(xyz, 2)
      self%n = npoints
      cutoff2 = self%cutoff*self%cutoff

      ! Allocate CSR offset and count arrays
      if (allocated(self%inl)) deallocate (self%inl)
      if (allocated(self%nnl)) deallocate (self%nnl)
      if (allocated(self%nlat)) deallocate (self%nlat)
      if (allocated(self%dist)) deallocate (self%dist)
      allocate (self%inl(npoints), source=0)
      allocate (self%nnl(npoints), source=0)

      ! Handle trivial cases
      if (npoints <= 1 .or. self%cutoff <= 0.0_wp) then
         allocate (self%nlat(0))
         allocate (self%dist(0))
         return
      end if

      ! --- Build uniform cell grid ---

      ! Compute bounding box
      xmin = xyz(1, 1); xmax = xyz(1, 1)
      ymin = xyz(2, 1); ymax = xyz(2, 1)
      zmin = xyz(3, 1); zmax = xyz(3, 1)
      do i = 2, npoints
         if (xyz(1, i) < xmin) xmin = xyz(1, i)
         if (xyz(1, i) > xmax) xmax = xyz(1, i)
         if (xyz(2, i) < ymin) ymin = xyz(2, i)
         if (xyz(2, i) > ymax) ymax = xyz(2, i)
         if (xyz(3, i) < zmin) zmin = xyz(3, i)
         if (xyz(3, i) > zmax) zmax = xyz(3, i)
      end do

      ! Compute grid dimensions
      cell_inv = 1.0_wp/self%cutoff
      nx = max(1, floor((xmax - xmin)*cell_inv) + 1)
      ny = max(1, floor((ymax - ymin)*cell_inv) + 1)
      nz = max(1, floor((zmax - zmin)*cell_inv) + 1)
      ncells = nx*ny*nz

      ! Assign points to cells
      allocate (head(ncells), source=0)
      allocate (next(npoints), source=0)
      allocate (cell_id(npoints))

      do i = 1, npoints
         cx = min(int((xyz(1, i) - xmin)*cell_inv), nx - 1)
         cy = min(int((xyz(2, i) - ymin)*cell_inv), ny - 1)
         cz = min(int((xyz(3, i) - zmin)*cell_inv), nz - 1)
         ic = cx + cy*nx + cz*nx*ny + 1
         cell_id(i) = ic
         next(i) = head(ic)
         head(ic) = i
      end do

      ! --- Count pass: determine nnl(i) for each point ---

      do i = 1, npoints
         xi = xyz(1, i)
         yi = xyz(2, i)
         zi = xyz(3, i)
         ic = cell_id(i)

         ! Decode cell coordinates
         cx = mod(ic - 1, nx)
         cy = mod((ic - 1)/nx, ny)
         cz = (ic - 1)/(nx*ny)

         ! Walk 27-cell stencil
         do dz = -1, 1
            cz_j = cz + dz
            if (cz_j < 0 .or. cz_j >= nz) cycle
            do dy = -1, 1
               cy_j = cy + dy
               if (cy_j < 0 .or. cy_j >= ny) cycle
               do dx = -1, 1
                  cx_j = cx + dx
                  if (cx_j < 0 .or. cx_j >= nx) cycle

                  jc = cx_j + cy_j*nx + cz_j*nx*ny + 1
                  j = head(jc)
                  do while (j > 0)
                     if (j /= i) then
                        dxij = xi - xyz(1, j)
                        dyij = yi - xyz(2, j)
                        dzij = zi - xyz(3, j)
                        d2 = dxij*dxij + dyij*dyij + dzij*dzij
                        if (d2 <= cutoff2) then
                           self%nnl(i) = self%nnl(i) + 1
                        end if
                     end if
                     j = next(j)
                  end do
               end do
            end do
         end do
      end do

      ! --- Prefix sum: compute inl from nnl ---

      img = 0
      do i = 1, npoints
         self%inl(i) = img
         img = img + self%nnl(i)
      end do
      capacity = img
      allocate (self%nlat(capacity))
      allocate (self%dist(capacity))

      ! Reset nnl for the fill pass (reuse as running counter)
      self%nnl = 0

      ! --- Fill pass: store neighbour indices and distances ---

      do i = 1, npoints
         xi = xyz(1, i)
         yi = xyz(2, i)
         zi = xyz(3, i)
         ic = cell_id(i)

         cx = mod(ic - 1, nx)
         cy = mod((ic - 1)/nx, ny)
         cz = (ic - 1)/(nx*ny)

         do dz = -1, 1
            cz_j = cz + dz
            if (cz_j < 0 .or. cz_j >= nz) cycle
            do dy = -1, 1
               cy_j = cy + dy
               if (cy_j < 0 .or. cy_j >= ny) cycle
               do dx = -1, 1
                  cx_j = cx + dx
                  if (cx_j < 0 .or. cx_j >= nx) cycle

                  jc = cx_j + cy_j*nx + cz_j*nx*ny + 1
                  j = head(jc)
                  do while (j > 0)
                     if (j /= i) then
                        dxij = xi - xyz(1, j)
                        dyij = yi - xyz(2, j)
                        dzij = zi - xyz(3, j)
                        d2 = dxij*dxij + dyij*dyij + dzij*dzij
                        if (d2 <= cutoff2) then
                           self%nnl(i) = self%nnl(i) + 1
                           self%nlat(self%inl(i) + self%nnl(i)) = j
                           self%dist(self%inl(i) + self%nnl(i)) = sqrt(d2)
                        end if
                     end if
                     j = next(j)
                  end do
               end do
            end do
         end do
      end do

      ! --- Optionally sort each point's neighbours by ascending distance ---
      if (self%sorted) then
         block
            type(error_type), allocatable :: sort_error
            do i = 1, npoints
               if (self%nnl(i) > 1) then
                  call qsort( &
                     self%dist(self%inl(i) + 1:self%inl(i) + self%nnl(i)), &
                     self%nlat(self%inl(i) + 1:self%inl(i) + self%nnl(i)), &
                     sort_error)
               end if
            end do
         end block
      end if

      deallocate (head, next, cell_id)

   end subroutine adjacency_list_update

   !> Return the neighbour indices for point i as an integer array.
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

   !> Deallocate all storage.
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

   !> Finalizer - delegates to destroy.
   subroutine adjacency_list_finalize(self)
      type(adjacency_list_type), intent(inout) :: self
      call self%destroy()
   end subroutine adjacency_list_finalize

end module moist_math_adjacency_list
