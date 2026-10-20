!> Indexed, conforming marching-tetrahedra surface geometry
module moist_cavity_mt_mesh
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: error_type, fatal_error
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_thread_slot
!$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads
   implicit none(type, external)
   private
   public :: mt_mesh_type, build_mt_mesh

   !> Piecewise planar zero set, with normals toward positive LSF values
   type :: mt_mesh_type
      !> Surface vertex coordinates (3, nv), Bohr
      real(wp), allocatable :: vertices(:, :)
      !> Counterclockwise outward vertex indices (3, nt)
      integer, allocatable :: triangles(:, :)
      !> Adjacent triangle opposite each local vertex (3, nt)
      integer, allocatable :: neighbors(:, :)
      !> Planar triangle areas (nt), Bohr squared
      real(wp), allocatable :: area(:)
      !> Planar triangle centroids (3, nt), Bohr
      real(wp), allocatable :: centroid(:, :)
      !> Outward unit normals (3, nt)
      real(wp), allocatable :: normal(:, :)
      !> Total planar area, Bohr squared
      real(wp) :: total_area = 0.0_wp
      !> Signed enclosed volume including cavities, Bohr cubed
      real(wp) :: total_volume = 0.0_wp
      !> Mean atomic position used for stable volume summation, Bohr
      real(wp) :: origin(3) = 0.0_wp
      !> Maximum Cartesian cell side, Bohr
      real(wp) :: cell_spacing = 0.0_wp
   contains
      procedure :: clear => clear_mesh
   end type mt_mesh_type

   !> Sparse unordered pair index, storing only intersected lattice edges
   type :: pair_table_type
      !> Smaller lattice or surface endpoint identifier
      integer(int64), allocatable :: first(:)
      !> Larger endpoint identifier
      integer(int64), allocatable :: second(:)
      !> Associated vertex or triangle-side index
      integer, allocatable :: value(:)
      !> Occupied table slots
      integer :: used = 0
   end type pair_table_type

   !> Geometrically grown surface buffers
   type :: surface_buffer_type
      !> Relative surface vertex coordinates
      real(wp), allocatable :: vertices(:, :)
      !> Surface triangle vertex indices
      integer, allocatable :: triangles(:, :)
      !> Stored vertices and triangles
      integer :: nv = 0, nt = 0
      !> Intersected lattice edge index
      type(pair_table_type) :: edges
   end type surface_buffer_type

   !> Cube vertices and six tetrahedra around the common body diagonal 1-8
   integer, parameter :: corner(3, 8) = reshape([ &
                                                0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, &
                                                0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1], [3, 8])
   !> Conforming tetrahedra with identical diagonals across shared cube faces
   integer, parameter :: tetra(4, 6) = reshape([ &
                                               1, 2, 4, 8, 1, 4, 3, 8, 1, 3, 7, 8, &
                                               1, 7, 5, 8, 1, 5, 6, 8, 1, 6, 2, 8], [4, 6])

contains

   !> Release all mesh storage and reset scalars
   !>
   !> @param[in,out] self Surface geometry
   subroutine clear_mesh(self)
      !> Surface geometry
      class(mt_mesh_type), intent(inout) :: self
      if (allocated(self%vertices)) deallocate (self%vertices)
      if (allocated(self%triangles)) deallocate (self%triangles)
      if (allocated(self%neighbors)) deallocate (self%neighbors)
      if (allocated(self%area)) deallocate (self%area)
      if (allocated(self%centroid)) deallocate (self%centroid)
      if (allocated(self%normal)) deallocate (self%normal)
      self%total_area = 0.0_wp
      self%total_volume = 0.0_wp
      self%origin = 0.0_wp
      self%cell_spacing = 0.0_wp
   end subroutine clear_mesh

   !> Sample a uniform lattice and retain every resolved zero-set component
   !>
   !> - Two scalar planes and a sparse surface-edge index bound working storage
   !> - Exact zero endpoints share their lattice vertex index across every edge
   !> - Serial indexing and parallel independent samples give deterministic output
   !>
   !> @param[in] lsf Updated negative-inside level-set model
   !> @param[in] xyz Atomic positions (3, nat), Bohr
   !> @param[in] spacing Positive maximum Cartesian cell side, Bohr
   !> @param[out] mesh Closed planar surface geometry
   !> @param[out] error Invalid input, allocation, evaluation or topology failure
   !> @param[in] nthreads Optional positive sampling team size
   subroutine build_mt_mesh(lsf, xyz, spacing, mesh, error, nthreads)
      !> Updated level-set model
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Atomic positions, Bohr
      real(wp), intent(in) :: xyz(:, :)
      !> Maximum Cartesian cell side, Bohr
      real(wp), intent(in) :: spacing
      !> Closed planar surface geometry
      type(mt_mesh_type), intent(out) :: mesh
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      !> Sampling team size
      integer, optional, intent(in) :: nthreads
      type(lsf_thread_slot), allocatable :: slots(:)
      type(surface_buffer_type) :: buffer
      !> Headroom for expanded bounds, cross products and volume products
      real(wp), parameter :: safe_length = huge(1.0_wp)**0.125_wp
      real(wp) :: lower(3), upper(3), extent(3), pad, origin(3)
      integer :: nt, i, stat, attempt, nc(3)
      logical :: boundary_hit

      if (.not. ieee_is_finite(spacing)) then
         call fatal_error(error, "MT mesh: spacing must be finite and positive")
         return
      end if
      if (spacing <= 0.0_wp .or. size(xyz, 1) /= 3 .or. size(xyz, 2) == 0) then
         call fatal_error(error, "MT mesh: invalid spacing or atomic coordinates")
         return
      end if
      if (.not. all(ieee_is_finite(xyz)) .or. .not. allocated(lsf%radii)) then
         call fatal_error(error, "MT mesh: finite coordinates and updated LSF radii required")
         return
      end if
      if (size(lsf%radii) /= size(xyz, 2)) then
         call fatal_error(error, "MT mesh: atomic coordinates and radii differ in size")
         return
      end if
      if (.not. all(ieee_is_finite(lsf%radii)) .or. any(lsf%radii <= 0.0_wp)) then
         call fatal_error(error, "MT mesh: finite positive radii required")
         return
      end if
      if (spacing > safe_length .or. any(abs(xyz) > safe_length) .or. any(lsf%radii > safe_length)) then
         call fatal_error(error, "MT mesh: coordinate or length scale exceeds safe arithmetic range")
         return
      end if
      nt = 1
!$    nt = omp_get_max_threads()
      if (present(nthreads)) then
         if (nthreads <= 0) then
            call fatal_error(error, "MT mesh: thread count must be positive")
            return
         end if
         nt = nthreads
      end if
      allocate (slots(nt), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "MT mesh: could not allocate thread slots")
         return
      end if
      do i = 1, nt
         allocate (slots(i)%lsf, source=lsf, stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "MT mesh: could not clone LSF")
            return
         end if
         ! Screened SvdW/CFC return zero for an empty active list, which is
         ! an exterior sentinel rather than a zero-set crossing
         call slots(i)%lsf%set_screening_threshold(0.0_wp)
         call slots(i)%lsf%set_max_deriv(0)
      end do
      ! Work about an atomic origin, avoiding large-coordinate volume cancellation
      origin = sum(xyz, dim=2)/real(size(xyz, 2), wp)
      do i = 1, 3
         lower(i) = minval(xyz(i, :) - origin(i) - lsf%radii)
         upper(i) = maxval(xyz(i, :) - origin(i) + lsf%radii)
      end do
      pad = maxval(lsf%radii) + 2.0_wp*spacing
      do attempt = 1, 8
         extent = upper - lower + 2.0_wp*pad
         if (any(extent/real(huge(1) - 2, wp) > spacing)) exit
         extent = extent/spacing
         if (.not. all(ieee_is_finite(extent))) exit
         if (any(extent >= real(huge(1) - 2, wp))) exit
         nc = max(1, ceiling(extent))
         ! Node identifiers and plane indexing must fit their integer kinds
         if (real(nc(1) + 1, wp)*real(nc(2) + 1, wp) > 0.5_wp*real(huge(1), wp)) exit
         if (product(real(nc + 1, wp)) > 0.25_wp*real(huge(1_int64), wp)) exit
         call extract_surface(slots, nt, lower - pad, nc, spacing, origin, buffer, boundary_hit, error)
         if (allocated(error)) return
         if (.not. boundary_hit) then
            ! The lattice-edge index is no longer needed during face adjacency
            deallocate (buffer%edges%first, buffer%edges%second, buffer%edges%value)
            call finish_mesh(buffer, origin, spacing, mesh, error)
            if (allocated(error)) call mesh%clear()
            return
         end if
         pad = 2.0_wp*pad
      end do
      call fatal_error(error, "MT mesh: bounding lattice overflow or exterior boundary not found")
   end subroutine build_mt_mesh

   !> Stream scalar planes and triangulate each tetrahedron in lattice order
   !>
   !> @param[in,out] slots Thread-local LSF copies
   !> @param[in] nt Sampling team size
   !> @param[in] lower Lower grid corner relative to origin
   !> @param[in] nc Cartesian cell counts
   !> @param[in] h Cartesian spacing
   !> @param[in] origin Atomic origin
   !> @param[out] buffer Indexed geometry buffers
   !> @param[out] boundary_hit Negative LSF on the outer lattice
   !> @param[out] error Sampling or allocation failure
   subroutine extract_surface(slots, nt, lower, nc, h, origin, buffer, boundary_hit, error)
      !> Thread-local models
      type(lsf_thread_slot), intent(inout) :: slots(:)
      !> Team size
      integer, intent(in) :: nt
      !> Relative lower grid corner
      real(wp), intent(in) :: lower(3)
      !> Cell counts
      integer, intent(in) :: nc(3)
      !> Grid spacing
      real(wp), intent(in) :: h
      !> Atomic origin
      real(wp), intent(in) :: origin(3)
      !> Indexed surface buffers
      type(surface_buffer_type), intent(out) :: buffer
      !> Outer boundary intersects the interior
      logical, intent(out) :: boundary_hit
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: planes(:, :, :)
      real(wp) :: p(3, 8), f(8)
      integer(int64) :: ids(8)
      integer :: stat, iz, iy, ix, k, it, bottom, top, q(3)

      boundary_hit = .false.
      allocate (planes(0:nc(1), 0:nc(2), 2), stat=stat)
      if (stat == 0) allocate (buffer%vertices(3, 4096), stat=stat)
      if (stat == 0) allocate (buffer%triangles(3, 4096), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "MT mesh: could not allocate lattice planes")
         return
      end if
      call resize_table(buffer%edges, 8192, error)
      if (allocated(error)) return
      bottom = 1
      top = 2
      call sample_plane(slots, nt, lower, nc, h, origin, 0, planes(:, :, bottom), error)
      if (allocated(error)) return
      ! Underflowed production pseudo-densities may report exterior zero.
      ! Genuine zero contacts remain valid only if final adjacency is closed.
      boundary_hit = any(planes(:, :, bottom) < 0.0_wp)
      if (boundary_hit) return
      do iz = 0, nc(3) - 1
         call sample_plane(slots, nt, lower, nc, h, origin, iz + 1, planes(:, :, top), error)
         if (allocated(error)) return
         boundary_hit = any(planes(0, :, top) < 0.0_wp) .or. any(planes(nc(1), :, top) < 0.0_wp) &
                        .or. any(planes(:, 0, top) < 0.0_wp) .or. any(planes(:, nc(2), top) < 0.0_wp)
         if (iz == nc(3) - 1) boundary_hit = boundary_hit .or. any(planes(:, :, top) < 0.0_wp)
         if (boundary_hit) return
         do iy = 0, nc(2) - 1
            do ix = 0, nc(1) - 1
               do k = 1, 8
                  q = [ix, iy, iz] + corner(:, k)
                  if (corner(3, k) == 0) then
                     f(k) = planes(q(1), q(2), bottom)
                  else
                     f(k) = planes(q(1), q(2), top)
                  end if
               end do
               if (all(f < 0.0_wp) .or. all(f >= 0.0_wp)) cycle
               do k = 1, 8
                  q = [ix, iy, iz] + corner(:, k)
                  p(:, k) = lower + h*real(q, wp)
                  ids(k) = 1_int64 + int(q(1), int64) + int(nc(1) + 1, int64) &
                           *(int(q(2), int64) + int(nc(2) + 1, int64)*int(q(3), int64))
               end do
               do it = 1, 6
                  call triangulate_tetra(p(:, tetra(:, it)), f(tetra(:, it)), &
                                         ids(tetra(:, it)), buffer, error)
                  if (allocated(error)) return
               end do
            end do
         end do
         k = bottom
         bottom = top
         top = k
      end do
   end subroutine extract_surface

   !> Evaluate one independent plane with thread-local caches
   !>
   !> @param[in,out] slots LSF clones
   !> @param[in] nt Sampling team size
   !> @param[in] lower Relative lower corner
   !> @param[in] nc Cartesian cell counts
   !> @param[in] h Cartesian spacing
   !> @param[in] origin Atomic origin
   !> @param[in] iz Plane index
   !> @param[out] values Lattice values (0:nx, 0:ny)
   !> @param[out] error First evaluation failure
   subroutine sample_plane(slots, nt, lower, nc, h, origin, iz, values, error)
      !> LSF clones
      type(lsf_thread_slot), intent(inout) :: slots(:)
      !> Team size
      integer, intent(in) :: nt
      !> Relative lower corner
      real(wp), intent(in) :: lower(3)
      !> Cell counts
      integer, intent(in) :: nc(3)
      !> Cartesian spacing
      real(wp), intent(in) :: h
      !> Atomic origin
      real(wp), intent(in) :: origin(3)
      !> Plane index
      integer, intent(in) :: iz
      !> Sampled values
      real(wp), intent(out) :: values(0:, 0:)
      !> Evaluation failure
      type(error_type), allocatable, intent(out) :: error
      type(error_type), allocatable :: local_error
      real(wp) :: point(3), val
      integer :: ix, iy, slot

      !$omp parallel do collapse(2) schedule(static) num_threads(nt) default(shared) &
      !$omp private(ix, iy, slot, point, val, local_error)
      do iy = 0, nc(2)
         do ix = 0, nc(1)
            slot = 1
!$          slot = omp_get_thread_num() + 1
            point = origin + lower + h*real([ix, iy, iz], wp)
            call slots(slot)%lsf%prepare(point, local_error)
            val = 0.0_wp
            if (.not. allocated(local_error)) then
               call slots(slot)%lsf%f0(val)
               if (.not. ieee_is_finite(val)) call fatal_error(local_error, "MT mesh: nonfinite LSF value")
            end if
            if (allocated(local_error)) then
               !$omp critical(mt_mesh_sample_error)
               if (.not. allocated(error)) call move_alloc(local_error, error)
               !$omp end critical(mt_mesh_sample_error)
               if (allocated(local_error)) deallocate (local_error)
            end if
            values(ix, iy) = val
         end do
      end do
      !$omp end parallel do
   end subroutine sample_plane

   !> Intersect the negative and nonnegative vertex partitions
   !>
   !> @param[in] p Tetrahedron coordinates (3, 4), relative to origin
   !> @param[in] f Lattice values (4)
   !> @param[in] ids Global lattice node identifiers (4)
   !> @param[in,out] buffer Surface buffers
   !> @param[out] error Allocation failure
   subroutine triangulate_tetra(p, f, ids, buffer, error)
      !> Relative coordinates
      real(wp), intent(in) :: p(3, 4)
      !> LSF values
      real(wp), intent(in) :: f(4)
      !> Lattice identifiers
      integer(int64), intent(in) :: ids(4)
      !> Surface buffers
      type(surface_buffer_type), intent(inout) :: buffer
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      integer :: neg(4), pos(4), nn, np, k, v(4), a(4), b(4), count
      real(wp) :: outward(3)

      nn = 0
      np = 0
      do k = 1, 4
         if (f(k) < 0.0_wp) then
            nn = nn + 1
            neg(nn) = k
         else
            np = np + 1
            pos(np) = k
         end if
      end do
      if (nn == 0 .or. np == 0) return
      outward = sum(p(:, pos(:np)), dim=2)/real(np, wp) - sum(p(:, neg(:nn)), dim=2)/real(nn, wp)
      count = 3
      if (nn == 1) then
         a(:3) = neg(1)
         b(:3) = pos(:3)
      else if (np == 1) then
         a(:3) = neg(:3)
         b(:3) = pos(1)
      else
         count = 4
         a = [neg(1), neg(1), neg(2), neg(2)]
         b = [pos(1), pos(2), pos(2), pos(1)]
      end if
      do k = 1, count
         call edge_vertex(p(:, a(k)), p(:, b(k)), f(a(k)), f(b(k)), ids(a(k)), ids(b(k)), &
                          buffer, v(k), error)
         if (allocated(error)) return
      end do
      call append_triangle(v(:3), outward, buffer, error)
      if (allocated(error)) return
      if (count == 4) call append_triangle(v([1, 3, 4]), outward, buffer, error)
   end subroutine triangulate_tetra

   !> Deduplicate edge crossings and exact-zero lattice endpoints
   !>
   !> @param[in] p Negative endpoint coordinates
   !> @param[in] q Nonnegative endpoint coordinates
   !> @param[in] fp Negative endpoint value
   !> @param[in] fq Nonnegative endpoint value
   !> @param[in] ip Negative endpoint node identifier
   !> @param[in] iq Nonnegative endpoint node identifier
   !> @param[in,out] buffer Surface buffers
   !> @param[out] vertex Surface vertex index
   !> @param[out] error Allocation or index overflow failure
   subroutine edge_vertex(p, q, fp, fq, ip, iq, buffer, vertex, error)
      !> Negative endpoint coordinates
      real(wp), intent(in) :: p(3)
      !> Nonnegative endpoint coordinates
      real(wp), intent(in) :: q(3)
      !> Negative endpoint value
      real(wp), intent(in) :: fp
      !> Nonnegative endpoint value
      real(wp), intent(in) :: fq
      !> Negative node identifier
      integer(int64), intent(in) :: ip
      !> Nonnegative node identifier
      integer(int64), intent(in) :: iq
      !> Geometry buffers
      type(surface_buffer_type), intent(inout) :: buffer
      !> Vertex index
      integer, intent(out) :: vertex
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: first, second
      integer :: slot, stat, capacity
      real(wp) :: fraction, scale
      real(wp), allocatable :: grown(:, :)

      first = min(ip, iq)
      second = max(ip, iq)
      if (fq <= 0.0_wp) then
         first = iq
         second = iq
      end if
      if (buffer%edges%used >= size(buffer%edges%value)/2) then
         if (2_int64*int(size(buffer%edges%value), int64) > int(huge(1), int64)) then
            call fatal_error(error, "MT mesh: edge table overflow")
            return
         end if
         call resize_table(buffer%edges, 2*size(buffer%edges%value), error)
         if (allocated(error)) return
      end if
      slot = pair_slot(buffer%edges, first, second)
      vertex = buffer%edges%value(slot)
      if (vertex /= 0) return
      capacity = size(buffer%vertices, 2)
      if (buffer%nv == capacity) then
         if (6_int64*int(capacity, int64) > int(huge(1), int64)) then
            call fatal_error(error, "MT mesh: vertex buffer overflow")
            return
         end if
         allocate (grown(3, 2*capacity), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "MT mesh: could not grow vertex buffer")
            return
         end if
         grown(:, :capacity) = buffer%vertices
         call move_alloc(grown, buffer%vertices)
      end if
      buffer%nv = buffer%nv + 1
      vertex = buffer%nv
      scale = max(abs(fp), abs(fq))
      fraction = (-fp/scale)/(fq/scale - fp/scale)
      buffer%vertices(:, vertex) = p + fraction*(q - p)
      if (fq <= 0.0_wp) buffer%vertices(:, vertex) = q
      buffer%edges%first(slot) = first
      buffer%edges%second(slot) = second
      buffer%edges%value(slot) = vertex
      buffer%edges%used = buffer%edges%used + 1
   end subroutine edge_vertex

   !> Append a nondegenerate triangle with outward winding
   !>
   !> @param[in] v Vertex indices
   !> @param[in] outward Vector toward positive LSF values
   !> @param[in,out] buffer Geometry buffers
   !> @param[out] error Allocation or index overflow failure
   subroutine append_triangle(v, outward, buffer, error)
      !> Vertex indices
      integer, intent(in) :: v(3)
      !> Positive-LSF direction
      real(wp), intent(in) :: outward(3)
      !> Geometry buffers
      type(surface_buffer_type), intent(inout) :: buffer
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: normal(3)
      integer :: capacity, stat
      integer, allocatable :: grown(:, :)
      if (v(1) == v(2) .or. v(1) == v(3) .or. v(2) == v(3)) return
      normal = cross(buffer%vertices(:, v(2)) - buffer%vertices(:, v(1)), &
                     buffer%vertices(:, v(3)) - buffer%vertices(:, v(1)))
      if (norm2(normal) <= 0.0_wp) return
      capacity = size(buffer%triangles, 2)
      if (buffer%nt == capacity) then
         if (6_int64*int(capacity, int64) > int(huge(1), int64)) then
            call fatal_error(error, "MT mesh: triangle buffer overflow")
            return
         end if
         allocate (grown(3, 2*capacity), stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "MT mesh: could not grow triangle buffer")
            return
         end if
         grown(:, :capacity) = buffer%triangles
         call move_alloc(grown, buffer%triangles)
      end if
      buffer%nt = buffer%nt + 1
      if (dot_product(normal, outward) < 0.0_wp) then
         buffer%triangles(:, buffer%nt) = v([1, 3, 2])
      else
         buffer%triangles(:, buffer%nt) = v
      end if
   end subroutine append_triangle

   !> Finalize planar measures and reciprocal manifold adjacency
   !>
   !> @param[in] buffer Indexed surface buffers
   !> @param[in] origin Atomic origin, Bohr
   !> @param[in] spacing Cartesian spacing, Bohr
   !> @param[out] mesh Complete surface geometry
   !> @param[out] error Allocation or topology failure
   subroutine finish_mesh(buffer, origin, spacing, mesh, error)
      !> Indexed geometry
      type(surface_buffer_type), intent(in) :: buffer
      !> Atomic origin
      real(wp), intent(in) :: origin(3)
      !> Cartesian spacing
      real(wp), intent(in) :: spacing
      !> Complete geometry
      type(mt_mesh_type), intent(out) :: mesh
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      type(pair_table_type) :: edges
      real(wp) :: p(3, 3), n(3), length, volume, correction, y, summed
      integer :: stat, it, k, v(3), a, b, slot, previous, other, side

      if (buffer%nt == 0) then
         call fatal_error(error, "MT mesh: lattice resolves no closed surface; reduce spacing")
         return
      end if
      allocate (mesh%vertices(3, buffer%nv), stat=stat)
      if (stat == 0) allocate (mesh%triangles(3, buffer%nt), stat=stat)
      if (stat == 0) allocate (mesh%neighbors(3, buffer%nt), stat=stat)
      if (stat == 0) allocate (mesh%area(buffer%nt), stat=stat)
      if (stat == 0) allocate (mesh%centroid(3, buffer%nt), stat=stat)
      if (stat == 0) allocate (mesh%normal(3, buffer%nt), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "MT mesh: could not allocate final geometry")
         return
      end if
      mesh%vertices = buffer%vertices(:, :buffer%nv) + spread(origin, 2, buffer%nv)
      mesh%triangles = buffer%triangles(:, :buffer%nt)
      mesh%neighbors = 0
      mesh%origin = origin
      mesh%cell_spacing = spacing
      call resize_table(edges, 8192, error)
      if (allocated(error)) return
      correction = 0.0_wp
      do it = 1, buffer%nt
         v = mesh%triangles(:, it)
         p = buffer%vertices(:, v)
         n = cross(p(:, 2) - p(:, 1), p(:, 3) - p(:, 1))
         length = norm2(n)
         mesh%area(it) = 0.5_wp*length
         mesh%normal(:, it) = n/length
         mesh%centroid(:, it) = origin + sum(p, dim=2)/3.0_wp
         volume = dot_product(p(:, 1), cross(p(:, 2), p(:, 3)))/6.0_wp
         y = volume - correction
         summed = mesh%total_volume + y
         correction = (summed - mesh%total_volume) - y
         mesh%total_volume = summed
         do k = 1, 3
            a = v(mod(k, 3) + 1)
            b = v(mod(k + 1, 3) + 1)
            if (edges%used >= size(edges%value)/2) then
               if (2_int64*int(size(edges%value), int64) > int(huge(1), int64)) then
                  call fatal_error(error, "MT mesh: adjacency table overflow")
                  return
               end if
               call resize_table(edges, 2*size(edges%value), error)
               if (allocated(error)) return
            end if
            slot = pair_slot(edges, int(min(a, b), int64), int(max(a, b), int64))
            previous = edges%value(slot)
            if (previous < 0) then
               call fatal_error(error, "MT mesh: nonmanifold surface edge")
               return
            else if (previous > 0) then
               other = (previous - 1)/3 + 1
               side = mod(previous - 1, 3) + 1
               mesh%neighbors(k, it) = other
               mesh%neighbors(side, other) = it
               edges%value(slot) = -previous
            else
               edges%first(slot) = int(min(a, b), int64)
               edges%second(slot) = int(max(a, b), int64)
               edges%value(slot) = 3*(it - 1) + k
               edges%used = edges%used + 1
            end if
         end do
      end do
      if (any(mesh%neighbors == 0)) then
         call fatal_error(error, "MT mesh: open surface edge")
         return
      end if
      mesh%total_area = sum(mesh%area)
   end subroutine finish_mesh

   !> Find an unordered-pair slot by linear probing
   !>
   !> @param[in] table Pair table, at most half full
   !> @param[in] first Smaller positive identifier
   !> @param[in] second Larger positive identifier
   pure function pair_slot(table, first, second) result(slot)
      !> Pair table
      type(pair_table_type), intent(in) :: table
      !> Smaller identifier
      integer(int64), intent(in) :: first
      !> Larger identifier
      integer(int64), intent(in) :: second
      !> Matching or empty slot
      integer :: slot
      !> Keep every multiply below 2**62, avoiding signed-overflow assumptions
      integer(int64), parameter :: mask = int(z"7fffffff", int64)
      integer(int64), parameter :: multiplier = int(z"45d9f3b", int64)
      integer(int64) :: mixed, other
      ! Fold every lattice-id bit before mixing neighboring sequential ids.
      ! A simple endpoint XOR produces severe clustering for triangle edges.
      mixed = ieor(first, ishft(first, -31))
      mixed = iand(ieor(mixed, ishft(first, -62)), mask)
      other = ieor(second, ishft(second, -31))
      other = iand(ieor(other, ishft(second, -62)), mask)
      mixed = iand(73856093_int64*mixed + 19349663_int64*other, mask)
      mixed = iand(ieor(mixed, ishft(mixed, -16))*multiplier, mask)
      mixed = iand(ieor(mixed, ishft(mixed, -16))*multiplier, mask)
      mixed = ieor(mixed, ishft(mixed, -16))
      slot = int(modulo(mixed, int(size(table%value), int64))) + 1
      do
         if (table%first(slot) == 0_int64) return
         if (table%first(slot) == first .and. table%second(slot) == second) return
         slot = mod(slot, size(table%value)) + 1
      end do
   end function pair_slot

   !> Allocate or rehash a sparse pair table
   !>
   !> @param[in,out] table Pair table
   !> @param[in] capacity New capacity
   !> @param[out] error Allocation failure
   subroutine resize_table(table, capacity, error)
      !> Pair table
      type(pair_table_type), intent(inout) :: table
      !> New capacity
      integer, intent(in) :: capacity
      !> Failure details
      type(error_type), allocatable, intent(out) :: error
      type(pair_table_type) :: grown
      integer :: stat, i, slot
      allocate (grown%first(capacity), stat=stat)
      if (stat == 0) allocate (grown%second(capacity), stat=stat)
      if (stat == 0) allocate (grown%value(capacity), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "MT mesh: could not allocate edge index")
         return
      end if
      grown%first = 0_int64
      grown%second = 0_int64
      grown%value = 0
      if (allocated(table%value)) then
         do i = 1, size(table%value)
            if (table%first(i) == 0_int64) cycle
            slot = pair_slot(grown, table%first(i), table%second(i))
            grown%first(slot) = table%first(i)
            grown%second(slot) = table%second(i)
            grown%value(slot) = table%value(i)
         end do
      end if
      call move_alloc(grown%first, table%first)
      call move_alloc(grown%second, table%second)
      call move_alloc(grown%value, table%value)
   end subroutine resize_table

   !> Three-dimensional vector product
   !>
   !> @param[in] a First vector
   !> @param[in] b Second vector
   pure function cross(a, b) result(c)
      !> First vector
      real(wp), intent(in) :: a(3)
      !> Second vector
      real(wp), intent(in) :: b(3)
      !> Vector product
      real(wp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function cross

end module moist_cavity_mt_mesh
