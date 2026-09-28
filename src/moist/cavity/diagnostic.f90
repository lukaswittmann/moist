!> Grid diagnostics that no cavity needs to run
module moist_cavity_diagnostic
   use mctc_env, only: wp, error_type, fatal_error
   use moist_cavity_type, only: cavity_type
   implicit none(type, external)
   private
   public :: find_disconnected_cavities

contains

   !> Count the connected islands of a cavity grid
   !>
   !> Two points are connected when they are closer than `threshold` times the
   !> mean nearest-neighbour spacing; the islands are labelled by breadth-first
   !> search on a cell list. A single island is the expected result
   !>
   !> @param[in]  cavity    Cavity with an allocated grid
   !> @param[out] sizes     Points per island (nislands), largest first
   !> @param[out] error     Empty grid or no measurable spacing
   !> @param[in]  threshold Connectivity radius in units of the mean spacing, 4 when absent
   subroutine find_disconnected_cavities(cavity, sizes, error, threshold)
      !> Cavity with an allocated grid
      class(cavity_type), intent(in) :: cavity
      !> Points per island
      integer, allocatable, intent(out) :: sizes(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Connectivity radius in units of the mean spacing
      real(wp), intent(in), optional :: threshold

      !> Cell list: first point of each cell and next point of the same cell
      integer, allocatable :: head(:), next(:)
      !> Cell index of each point along x, y, z
      integer, allocatable :: cell(:, :)
      !> Breadth-first queue and visited flags
      integer, allocatable :: queue(:)
      logical, allocatable :: visited(:)
      real(wp) :: min_xyz(3), bbox(3), cell_size, spacing, dist2, min_dist2, thrs
      integer :: nc(3), i, ncell, neighbour, qhead, qtail, current, radius
      integer :: nislands, nspacing

      thrs = 4.0_wp
      if (present(threshold)) thrs = threshold
      if (.not. allocated(cavity%xyz) .or. cavity%ngrid <= 0) then
         call fatal_error(error, "find_disconnected_cavities: no grid points")
         return
      end if

      min_xyz = minval(cavity%xyz, dim=2)
      bbox = maxval(cavity%xyz, dim=2) - min_xyz
      cell_size = max(1.0e-12_wp, product(bbox))/real(cavity%ngrid, wp)
      cell_size = max(1.0e-6_wp, cell_size**(1.0_wp/3.0_wp))

      ! First pass: mean nearest-neighbour spacing on a coarse cell list
      call bin_points(cell_size)
      spacing = 0.0_wp
      nspacing = 0
      do i = 1, cavity%ngrid
         min_dist2 = huge(1.0_wp)
         do radius = 0, 2
            call nearest_in_shell(i, radius, min_dist2)
            if (min_dist2 < huge(1.0_wp)) exit
         end do
         if (min_dist2 < huge(1.0_wp)) then
            spacing = spacing + sqrt(min_dist2)
            nspacing = nspacing + 1
         end if
      end do
      if (nspacing == 0 .or. spacing <= 0.0_wp) then
         call fatal_error(error, "find_disconnected_cavities: could not estimate grid spacing")
         return
      end if
      spacing = spacing/real(nspacing, wp)

      ! Second pass: breadth-first labelling with the connectivity radius as cell size
      cell_size = thrs*spacing
      call bin_points(cell_size)
      allocate (queue(cavity%ngrid), visited(cavity%ngrid), sizes(cavity%ngrid))
      visited = .false.
      nislands = 0
      do i = 1, cavity%ngrid
         if (visited(i)) cycle
         nislands = nislands + 1
         sizes(nislands) = 0
         qhead = 1
         qtail = 1
         queue(1) = i
         visited(i) = .true.
         do while (qhead <= qtail)
            current = queue(qhead)
            qhead = qhead + 1
            sizes(nislands) = sizes(nislands) + 1
            call visit_shell(current)
         end do
      end do
      sizes = sizes(:nislands)
      call sort_descending(sizes)

   contains

      !> Rebuild the cell list for one cell size
      subroutine bin_points(edge)
         !> Cell edge, bohr
         real(wp), intent(in) :: edge
         integer :: p, lin
         nc = max(1, int(bbox/edge) + 1)
         ncell = product(nc)
         if (allocated(head)) deallocate (head)
         if (.not. allocated(next)) allocate (next(cavity%ngrid), cell(3, cavity%ngrid))
         allocate (head(ncell))
         head = 0
         do p = 1, cavity%ngrid
            cell(:, p) = min(nc, max(1, int((cavity%xyz(:, p) - min_xyz)/edge) + 1))
            lin = cell(1, p) + nc(1)*(cell(2, p) - 1 + nc(2)*(cell(3, p) - 1))
            next(p) = head(lin)
            head(lin) = p
         end do
      end subroutine bin_points

      !> Squared distance to the nearest other point within `radius` cells of point `p`
      subroutine nearest_in_shell(p, radius, best)
         !> Point index
         integer, intent(in) :: p
         !> Cell radius of the search
         integer, intent(in) :: radius
         !> Running minimum of the squared distance
         real(wp), intent(inout) :: best
         integer :: ix, iy, iz, lin
         do iz = max(1, cell(3, p) - radius), min(nc(3), cell(3, p) + radius)
            do iy = max(1, cell(2, p) - radius), min(nc(2), cell(2, p) + radius)
               do ix = max(1, cell(1, p) - radius), min(nc(1), cell(1, p) + radius)
                  lin = ix + nc(1)*(iy - 1 + nc(2)*(iz - 1))
                  neighbour = head(lin)
                  do while (neighbour /= 0)
                     if (neighbour /= p) then
                        dist2 = sum((cavity%xyz(:, neighbour) - cavity%xyz(:, p))**2)
                        if (dist2 < best) best = dist2
                     end if
                     neighbour = next(neighbour)
                  end do
               end do
            end do
         end do
      end subroutine nearest_in_shell

      !> Enqueue every unvisited point within one cell and the connectivity radius of `p`
      subroutine visit_shell(p)
         !> Point index
         integer, intent(in) :: p
         integer :: ix, iy, iz, lin
         do iz = max(1, cell(3, p) - 1), min(nc(3), cell(3, p) + 1)
            do iy = max(1, cell(2, p) - 1), min(nc(2), cell(2, p) + 1)
               do ix = max(1, cell(1, p) - 1), min(nc(1), cell(1, p) + 1)
                  lin = ix + nc(1)*(iy - 1 + nc(2)*(iz - 1))
                  neighbour = head(lin)
                  do while (neighbour /= 0)
                     if (.not. visited(neighbour)) then
                        dist2 = sum((cavity%xyz(:, neighbour) - cavity%xyz(:, p))**2)
                        if (dist2 <= cell_size*cell_size) then
                           visited(neighbour) = .true.
                           qtail = qtail + 1
                           queue(qtail) = neighbour
                        end if
                     end if
                     neighbour = next(neighbour)
                  end do
               end do
            end do
         end do
      end subroutine visit_shell

      !> Insertion sort, largest first
      subroutine sort_descending(values)
         !> Values to sort in place
         integer, intent(inout) :: values(:)
         integer :: a, b, key
         do a = 2, size(values)
            key = values(a)
            b = a - 1
            do while (b >= 1)
               if (values(b) >= key) exit
               values(b + 1) = values(b)
               b = b - 1
            end do
            values(b + 1) = key
         end do
      end subroutine sort_descending

   end subroutine find_disconnected_cavities

end module moist_cavity_diagnostic
