!> Connected integration patches on a retained planar triangle mesh
module moist_cavity_mt_patches
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp, error_type, fatal_error
   use moist_cavity_mt_mesh, only: mt_mesh_type
   implicit none(type, external)
   private
   public :: mt_patch_type, build_mt_patches

   !> One on-mesh representative and retained geometric area per patch
   type :: mt_patch_type
      !> Representatives (3, npatch), bohr
      real(wp), allocatable :: xyz(:, :)
      !> Geometric patch areas (npatch), bohr^2
      real(wp), allocatable :: area(:)
      !> Connected mesh component of each patch
      integer, allocatable :: component(:)
      !> Original mesh face containing each representative
      integer, allocatable :: face(:)
      !> Attainable mean area within each component, bohr^2
      real(wp), allocatable :: component_target(:)
      !> Component too small for an integer point count within tolerance
      logical, allocatable :: area_exception(:)
      !> Offsets of each patch's face members (npatch + 1)
      integer, allocatable :: member_offset(:)
      !> Mesh face of each member
      integer, allocatable :: member_face(:)
      !> Face area inside the patch, bohr^2
      real(wp), allocatable :: member_area(:)
      !> Every patch meets its area bound
      logical :: balanced = .true.
      !> Largest relative deviation from the patch area bound centers
      real(wp) :: max_deviation = 0.0_wp
   contains
      procedure :: clear => clear_mt_patches
   end type mt_patch_type

   !> Area-bounded triangle pieces and their intrinsic adjacency
   type :: patch_graph_type
      !> Piece centroids (3, npiece), bohr
      real(wp), allocatable :: xyz(:, :)
      !> Piece areas, bohr^2
      real(wp), allocatable :: area(:)
      !> Signed reciprocal subdivision count; negative for inverted pieces
      real(wp), allocatable :: scale(:)
      !> Original face of each piece
      integer, allocatable :: face(:)
      !> Compressed adjacency offsets
      integer, allocatable :: offset(:)
      !> Adjacent pieces
      integer, allocatable :: adjacent(:)
      !> Intrinsic edge lengths, bohr
      real(wp), allocatable :: length(:)
   end type patch_graph_type

   !> Indexed minimum heap with deterministic index tie breaking
   type :: index_heap_type
      !> Heap entries
      integer, allocatable :: node(:)
      !> Heap position indexed by graph node, zero when absent
      integer, allocatable :: position(:)
      !> Priority indexed by graph node
      real(wp), allocatable :: key(:)
      !> Number of active entries
      integer :: count = 0
   contains
      procedure :: set => heap_set
      procedure :: pop => heap_pop
      procedure :: reset => heap_reset
   end type index_heap_type

   !> Binary minimum heap of (key, node) pairs; superseded entries are skipped
   type :: pair_heap_type
      !> Entry priorities
      real(wp), allocatable :: key(:)
      !> Entry nodes
      integer, allocatable :: node(:)
      !> Number of active entries
      integer :: count = 0
   end type pair_heap_type

contains

   !> Release patch state
   !>
   !> @param[in,out] self Patch storage
   subroutine clear_mt_patches(self)
      !> Patch storage
      class(mt_patch_type), intent(inout) :: self
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%area)) deallocate (self%area)
      if (allocated(self%component)) deallocate (self%component)
      if (allocated(self%face)) deallocate (self%face)
      if (allocated(self%component_target)) deallocate (self%component_target)
      if (allocated(self%area_exception)) deallocate (self%area_exception)
      if (allocated(self%member_offset)) deallocate (self%member_offset)
      if (allocated(self%member_face)) deallocate (self%member_face)
      if (allocated(self%member_area)) deallocate (self%member_area)
      self%balanced = .true.
      self%max_deviation = 0.0_wp
   end subroutine clear_mt_patches

   !> Partition each component by capacity-balanced intrinsic Voronoi cells
   !>
   !> Every positive face contributes its full area. Triangle subdivision bounds
   !> area quantization; no geometry is projected onto the underlying level set.
   !>
   !> - Farthest-point seeds and Voronoi weights come from blocks of about four
   !>   piece bounds, whose passes are cheap
   !> - The weights then warm-start a few passes on aggregates of one piece bound
   !> - Transfers of single pieces repair the remaining cells
   !>
   !> Dijkstra trees and connectivity-checked transfers keep cells connected. An
   !> unattainable integer point count is recorded per component. A partition
   !> missing the tolerance otherwise keeps its best connected cells and clears
   !> `balanced`.
   !>
   !> @param[in] mesh Retained planar geometry
   !> @param[in] target_area Requested area per patch, positive, bohr^2
   !> @param[in] area_tolerance Relative area tolerance, strictly between zero and one
   !> @param[out] patches Integration patches
   !> @param[out] error Failure state
   subroutine build_mt_patches(mesh, target_area, area_tolerance, patches, error)
      !> Retained geometry
      type(mt_mesh_type), intent(in) :: mesh
      !> Requested patch area, bohr^2
      real(wp), intent(in) :: target_area
      !> Relative tolerance
      real(wp), intent(in) :: area_tolerance
      !> Integration patches
      type(mt_patch_type), intent(out) :: patches
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      !> Block area bound relative to the piece bound
      real(wp), parameter :: block_factor = 4.0_wp
      !> Block passes and stall limit
      integer, parameter :: block_passes = 160, block_patience = 8
      !> Warm-started aggregate passes; diffusion repairs the remaining cells
      integer, parameter :: fine_passes = 12, fine_patience = 4
      type(patch_graph_type) :: graph, coarse, blocks
      type(index_heap_type) :: work, farthest
      integer, allocatable :: face_component(:), piece_component(:), component_of_block(:), queue(:), counts(:)
      integer, allocatable :: starts(:), members(:), cluster(:), block_of(:), seeds(:), seed_block(:)
      integer, allocatable :: labels(:), node_labels(:), block_labels(:), anchors(:), seen(:), cursor(:)
      real(wp), allocatable :: component_area(:), goal(:), weights(:), distance(:), mass(:)
      real(wp), allocatable :: moment(:, :), closest(:), centers(:), normal(:, :), coarse_normal(:, :)
      logical, allocatable :: exception(:)
      real(wp) :: candidate, threshold, norm, allowance, center, slack
      real(wp) :: edge_a(3), edge_b(3), projected(3)
      integer :: nf, n, nn, nb, nc, np, nm, phase, c, f, g, p, u, v, j, k, head, tail, stat, iteration, stamp
      logical :: moved

      if (.not. ieee_is_finite(target_area) .or. .not. ieee_is_finite(area_tolerance)) then
         call fatal_error(error, "MT patch area and tolerance must be finite")
         return
      end if
      if (target_area <= 0.0_wp .or. area_tolerance <= 0.0_wp .or. area_tolerance >= 1.0_wp) then
         call fatal_error(error, "MT patches require positive area and a relative tolerance in (0,1)")
         return
      end if
      threshold = 0.1_wp*area_tolerance*target_area
      if (threshold <= 0.0_wp) then
         call fatal_error(error, "MT patch subdivision area underflows working precision")
         return
      end if
      if (.not. allocated(mesh%area) .or. .not. allocated(mesh%triangles) .or. &
          .not. allocated(mesh%neighbors) .or. .not. allocated(mesh%vertices) .or. &
          .not. allocated(mesh%normal)) then
         call fatal_error(error, "MT patches require a retained triangle mesh")
         return
      end if
      nf = size(mesh%area)
      if (nf == 0) then
         call fatal_error(error, "MT patches require a nonempty triangle mesh")
         return
      end if
      if (any(mesh%area <= 0.0_wp) .or. .not. all(ieee_is_finite(mesh%area))) then
         call fatal_error(error, "MT patches require positive finite face areas")
         return
      end if
      allocate (face_component(nf), stat=stat)
      if (stat == 0) allocate (queue(nf), stat=stat)
      if (stat == 0) allocate (component_area(nf), stat=stat)
      if (stat == 0) allocate (counts(nf), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT patch components")
         return
      end if
      face_component = 0
      component_area = 0.0_wp
      nc = 0
      do f = 1, nf
         if (face_component(f) /= 0) cycle
         nc = nc + 1
         head = 1
         tail = 1
         queue(1) = f
         face_component(f) = nc
         do while (head <= tail)
            u = queue(head)
            head = head + 1
            component_area(nc) = component_area(nc) + mesh%area(u)
            do k = 1, 3
               v = mesh%neighbors(k, u)
               if (v <= 0) cycle
               if (face_component(v) /= 0) cycle
               tail = tail + 1
               queue(tail) = v
               face_component(v) = nc
            end do
         end do
         if (target_area < component_area(nc)/real(huge(np) - 1, wp)) then
            call fatal_error(error, "MT patch count exceeds integer capacity")
            return
         end if
         counts(nc) = max(1, nint(component_area(nc)/target_area))
         if (abs(component_area(nc)/real(counts(nc), wp)/target_area - 1.0_wp) > area_tolerance) then
            k = max(1, ceiling(component_area(nc)/target_area))
            if (abs(component_area(nc)/real(k, wp)/target_area - 1.0_wp) <= area_tolerance) counts(nc) = k
         end if
      end do
      if (sum(int(counts(:nc), int64)) > int(huge(np), int64)) then
         call fatal_error(error, "MT patch count exceeds integer capacity")
         return
      end if
      np = sum(counts(:nc))
      call make_patch_graph(mesh, threshold, graph, error)
      if (allocated(error)) return
      n = size(graph%area)
      ! Seeds and weights come from blocks; aggregates refine them at full
      ! area resolution before single pieces polish the bounds
      allocate (normal(3, n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT piece normals")
         return
      end if
      normal = mesh%normal(:, graph%face)
      call coarsen_patch_graph(graph, normal, threshold, coarse, coarse_normal, cluster, error)
      if (allocated(error)) return
      nn = size(coarse%area)
      call coarsen_patch_graph(coarse, coarse_normal, block_factor*threshold, blocks, normal, block_of, error)
      if (allocated(error)) return
      nb = size(blocks%area)
      deallocate (queue)
      allocate (piece_component(n), stat=stat)
      if (stat == 0) allocate (component_of_block(nb), stat=stat)
      if (stat == 0) allocate (queue(n), stat=stat)
      if (stat == 0) allocate (members(nb), stat=stat)
      if (stat == 0) allocate (starts(nc + 1), stat=stat)
      if (stat == 0) allocate (seeds(np), stat=stat)
      if (stat == 0) allocate (seed_block(np), stat=stat)
      if (stat == 0) allocate (labels(n), stat=stat)
      if (stat == 0) allocate (node_labels(nn), stat=stat)
      if (stat == 0) allocate (block_labels(nb), stat=stat)
      if (stat == 0) allocate (seen(n), stat=stat)
      if (stat == 0) allocate (goal(np), stat=stat)
      if (stat == 0) allocate (centers(np), stat=stat)
      if (stat == 0) allocate (weights(np), stat=stat)
      if (stat == 0) allocate (distance(nn), stat=stat)
      if (stat == 0) allocate (closest(np), stat=stat)
      if (stat == 0) allocate (mass(np), stat=stat)
      if (stat == 0) allocate (exception(np), stat=stat)
      if (stat == 0) allocate (moment(3, np), stat=stat)
      if (stat == 0) allocate (work%node(nb), stat=stat)
      if (stat == 0) allocate (work%position(nb), stat=stat)
      if (stat == 0) allocate (work%key(nb), stat=stat)
      if (stat == 0) allocate (farthest%node(nb), stat=stat)
      if (stat == 0) allocate (farthest%position(nb), stat=stat)
      if (stat == 0) allocate (farthest%key(nb), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT patch partition")
         return
      end if
      piece_component = face_component(graph%face)
      ! Aggregates never cross components, which share no piece adjacency
      do u = 1, n
         component_of_block(block_of(cluster(u))) = piece_component(u)
      end do
      queue(:nc) = 0
      do u = 1, nb
         c = component_of_block(u)
         queue(c) = queue(c) + 1
      end do
      starts(1) = 1
      do c = 1, nc
         starts(c + 1) = starts(c) + queue(c)
      end do
      queue(:nc) = starts(:nc)
      do u = 1, nb
         c = component_of_block(u)
         members(queue(c)) = u
         queue(c) = queue(c) + 1
      end do
      work%position = 0
      farthest%position = 0
      distance = huge(1.0_wp)
      p = 0
      do c = 1, nc
         if (counts(c) > starts(c + 1) - starts(c)) then
            call fatal_error(error, "MT patch subdivision has fewer blocks than patches")
            return
         end if
         do j = starts(c), starts(c + 1) - 1
            call farthest%set(members(j), -huge(1.0_wp))
         end do
         u = members(starts(c))
         do f = 1, counts(c)
            p = p + 1
            seed_block(p) = u
            goal(p) = component_area(c)/real(counts(c), wp)
            exception(p) = abs(goal(p)/target_area - 1.0_wp) > area_tolerance
            distance(u) = 0.0_wp
            call work%set(u, 0.0_wp)
            call farthest%set(u, 0.0_wp)
            do while (work%count > 0)
               call work%pop(v)
               do k = blocks%offset(v), blocks%offset(v + 1) - 1
                  g = blocks%adjacent(k)
                  candidate = distance(v) + blocks%length(k)
                  if (candidate >= distance(g)) cycle
                  distance(g) = candidate
                  call work%set(g, candidate)
                  call farthest%set(g, -candidate)
               end do
            end do
            if (f < counts(c)) u = farthest%node(1)
         end do
         call farthest%reset()
      end do
      ! Each block seed continues at its aggregate nearest the block centroid
      closest = huge(1.0_wp)
      seeds = 0
      allocate (anchors(nb), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT patch seeds")
         return
      end if
      anchors = 0
      do p = 1, np
         anchors(seed_block(p)) = p
      end do
      do u = 1, nn
         p = anchors(block_of(u))
         if (p == 0) cycle
         norm = norm2(coarse%xyz(:, u) - blocks%xyz(:, block_of(u)))
         if (norm >= closest(p)) cycle
         closest(p) = norm
         seeds(p) = u
      end do
      deallocate (anchors)
      centers = target_area
      where (exception) centers = goal
      weights = 0.0_wp
      call balance_weights(blocks, seed_block, goal, centers, area_tolerance, block_passes, block_patience, &
         & weights, block_labels, distance, error)
      if (allocated(error)) return
      call balance_weights(coarse, seeds, goal, centers, area_tolerance, fine_passes, fine_patience, &
         & weights, node_labels, distance, error)
      if (allocated(error)) return
      ! Connected aggregates expand to connected cells of single pieces
      mass = 0.0_wp
      do u = 1, n
         labels(u) = node_labels(cluster(u))
         mass(labels(u)) = mass(labels(u)) + graph%area(u)
      end do
      seen = 0
      stamp = 0
      ! Local capacity diffusion cannot disconnect a donor or a recipient.
      ! Coarse transfers leave cells near their goals untouched; single-piece
      ! resolution follows only where the bounds still fail.
      do phase = 1, 2
         do iteration = 1, 80
            if (acceptable_areas(mass, goal, exception, target_area, area_tolerance)) exit
            moved = .false.
            do u = 1, n
               p = labels(u)
               slack = 0.0_wp
               if (phase == 1) slack = 0.5_wp*area_tolerance*goal(p)
               do k = graph%offset(u), graph%offset(u + 1) - 1
                  v = graph%adjacent(k)
                  j = labels(v)
                  if (j == p) cycle
                  if (mass(p) <= graph%area(u)) cycle
                  threshold = max(slack, graph%area(u)*(1.0_wp + 16.0_wp*epsilon(1.0_wp)))
                  if ((mass(p) - goal(p)) - (mass(j) - goal(j)) <= threshold) cycle
                  if (.not. removable_piece(graph, labels, u, seen, stamp, queue)) cycle
                  labels(u) = j
                  mass(p) = mass(p) - graph%area(u)
                  mass(j) = mass(j) + graph%area(u)
                  moved = .true.
                  exit
               end do
            end do
            if (.not. moved) exit
         end do
      end do
      allocate (patches%xyz(3, np), stat=stat)
      if (stat == 0) allocate (patches%area(np), stat=stat)
      if (stat == 0) allocate (patches%component(np), stat=stat)
      if (stat == 0) allocate (patches%face(np), stat=stat)
      if (stat == 0) allocate (patches%component_target(np), stat=stat)
      if (stat == 0) allocate (patches%area_exception(np), stat=stat)
      if (stat == 0) allocate (patches%member_offset(np + 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT integration patches")
         return
      end if
      moment = 0.0_wp
      mass = 0.0_wp
      do u = 1, n
         p = labels(u)
         mass(p) = mass(p) + graph%area(u)
         moment(:, p) = moment(:, p) + graph%area(u)*graph%xyz(:, u)
      end do
      do p = 1, np
         moment(:, p) = moment(:, p)/mass(p)
      end do
      closest = huge(1.0_wp)
      do u = 1, n
         p = labels(u)
         f = graph%face(u)
         edge_a = graph%scale(u)*(mesh%vertices(:, mesh%triangles(2, f)) &
                                  - mesh%vertices(:, mesh%triangles(1, f)))
         edge_b = graph%scale(u)*(mesh%vertices(:, mesh%triangles(3, f)) &
                                  - mesh%vertices(:, mesh%triangles(1, f)))
         call project_piece(moment(:, p) - graph%xyz(:, u), edge_a, edge_b, projected)
         projected = projected + graph%xyz(:, u)
         norm = norm2(projected - moment(:, p))
         if (norm >= closest(p)) cycle
         closest(p) = norm
         patches%xyz(:, p) = projected
         patches%face(p) = graph%face(u)
         patches%component(p) = piece_component(u)
      end do
      patches%area = mass
      patches%component_target = goal
      patches%area_exception = exception
      ! Pieces of one face are contiguous; each run within one patch is a member
      allocate (cursor(np), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT patch members")
         return
      end if
      cursor = 0
      do u = 1, n
         if (u > 1) then
            if (graph%face(u) == graph%face(u - 1) .and. labels(u) == labels(u - 1)) cycle
         end if
         cursor(labels(u)) = cursor(labels(u)) + 1
      end do
      patches%member_offset(1) = 1
      do p = 1, np
         patches%member_offset(p + 1) = patches%member_offset(p) + cursor(p)
      end do
      nm = patches%member_offset(np + 1) - 1
      allocate (patches%member_face(nm), stat=stat)
      if (stat == 0) allocate (patches%member_area(nm), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT patch members")
         return
      end if
      cursor = patches%member_offset(:np) - 1
      do u = 1, n
         p = labels(u)
         if (u > 1) then
            if (graph%face(u) == graph%face(u - 1) .and. p == labels(u - 1)) then
               patches%member_area(cursor(p)) = patches%member_area(cursor(p)) + graph%area(u)
               cycle
            end if
         end if
         cursor(p) = cursor(p) + 1
         patches%member_face(cursor(p)) = graph%face(u)
         patches%member_area(cursor(p)) = graph%area(u)
      end do
      ! The final check includes the exact summation used for exported areas.
      allowance = 64.0_wp*epsilon(1.0_wp)*maxval(goal)
      patches%balanced = acceptable_areas(mass, goal, exception, target_area, area_tolerance, allowance)
      patches%max_deviation = 0.0_wp
      do p = 1, np
         center = target_area
         if (exception(p)) center = goal(p)
         patches%max_deviation = max(patches%max_deviation, abs(mass(p) - center)/center)
      end do
   end subroutine build_mt_patches

   !> Test requested bounds, using the attainable mean for rounding exceptions
   !>
   !> @param[in] mass Cell areas
   !> @param[in] goal Component mean areas
   !> @param[in] exception Unavoidable integer-count exceptions
   !> @param[in] target Requested area
   !> @param[in] tolerance Relative tolerance
   !> @param[in] allowance Optional absolute roundoff allowance
   pure function acceptable_areas(mass, goal, exception, target, tolerance, allowance) result(ok)
      !> Cell areas
      real(wp), intent(in) :: mass(:)
      !> Component mean areas
      real(wp), intent(in) :: goal(:)
      !> Integer-count exceptions
      logical, intent(in) :: exception(:)
      !> Requested area
      real(wp), intent(in) :: target
      !> Relative tolerance
      real(wp), intent(in) :: tolerance
      !> Absolute roundoff allowance
      real(wp), intent(in), optional :: allowance
      !> All areas satisfy their bounds
      logical :: ok
      real(wp) :: center, slack
      integer :: p
      slack = 0.0_wp
      if (present(allowance)) slack = allowance
      ok = .true.
      do p = 1, size(mass)
         center = target
         if (exception(p)) center = goal(p)
         if (abs(mass(p) - center) <= tolerance*center + slack) cycle
         ok = .false.
         return
      end do
   end function acceptable_areas

   !> Assemble adjacency without joining geometrically close unrelated faces
   !>
   !> @param[in] mesh Original mesh
   !> @param[in] max_area Largest triangle-piece area
   !> @param[out] graph Intrinsic piece graph
   !> @param[out] error Failure state
   subroutine make_patch_graph(mesh, max_area, graph, error)
      !> Original mesh
      type(mt_mesh_type), intent(in) :: mesh
      !> Largest piece area, bohr^2
      real(wp), intent(in) :: max_area
      !> Piece graph
      type(patch_graph_type), intent(out) :: graph
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: divisions(:), first(:), edge(:, :), degree(:), cursor(:)
      real(wp), allocatable :: edge_length(:)
      integer :: nf, n, f, g, m, l, i, j, u, v, ne, side, other, a, b, ia, ib, stat, e, graph_limit
      integer :: endpoints(2, 3), ids(3), pair(2), reverse_pair(2)
      integer(int64) :: total, left_end, right_end
      real(wp) :: x, y, ratio, lower, upper, start(3), finish(3), edge_vector(3)

      nf = size(mesh%area)
      allocate (divisions(nf), stat=stat)
      if (stat == 0) allocate (first(nf + 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT triangle subdivision")
         return
      end if
      total = 0_int64
      graph_limit = huge(n)
      graph_limit = graph_limit/6
      first(1) = 1
      do f = 1, nf
         ratio = sqrt(mesh%area(f))/sqrt(max_area)
         if (ratio > sqrt(real(graph_limit, wp))) then
            call fatal_error(error, "MT patch triangle subdivision exceeds integer capacity")
            return
         end if
         divisions(f) = max(1, ceiling(ratio))
         total = total + int(divisions(f), int64)**2
         if (total > int(graph_limit, int64)) then
            call fatal_error(error, "MT patch graph exceeds integer capacity")
            return
         end if
         first(f + 1) = int(total) + 1
      end do
      n = int(total)
      allocate (graph%xyz(3, n), stat=stat)
      if (stat == 0) allocate (graph%area(n), stat=stat)
      if (stat == 0) allocate (graph%scale(n), stat=stat)
      if (stat == 0) allocate (graph%face(n), stat=stat)
      if (stat == 0) allocate (edge(2, 6*n), stat=stat)
      if (stat == 0) allocate (edge_length(6*n), stat=stat)
      if (stat == 0) allocate (degree(n), stat=stat)
      if (stat == 0) allocate (cursor(n), stat=stat)
      if (stat == 0) allocate (graph%offset(n + 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT triangle-piece graph")
         return
      end if
      ne = 0
      do f = 1, nf
         m = divisions(f)
         ids = mesh%triangles(:, f)
         do i = 0, m - 1
            do j = 0, m - i - 1
               u = piece_index(first(f), m, i, j)
               x = (real(i, wp) + 1.0_wp/3.0_wp)/real(m, wp)
               y = (real(j, wp) + 1.0_wp/3.0_wp)/real(m, wp)
               graph%xyz(:, u) = mesh%vertices(:, ids(1)) &
                                 + x*(mesh%vertices(:, ids(2)) - mesh%vertices(:, ids(1))) &
                                 + y*(mesh%vertices(:, ids(3)) - mesh%vertices(:, ids(1)))
               graph%area(u) = mesh%area(f)/real(m, wp)**2
               graph%scale(u) = 1.0_wp/real(m, wp)
               graph%face(u) = f
               if (i + j >= m - 1) cycle
               v = u + 1
               x = (real(i, wp) + 2.0_wp/3.0_wp)/real(m, wp)
               y = (real(j, wp) + 2.0_wp/3.0_wp)/real(m, wp)
               graph%xyz(:, v) = mesh%vertices(:, ids(1)) &
                                 + x*(mesh%vertices(:, ids(2)) - mesh%vertices(:, ids(1))) &
                                 + y*(mesh%vertices(:, ids(3)) - mesh%vertices(:, ids(1)))
               graph%area(v) = graph%area(u)
               graph%scale(v) = -graph%scale(u)
               graph%face(v) = f
               edge(:, ne + 1) = [u, v]
               edge(:, ne + 2) = [v, piece_index(first(f), m, i + 1, j)]
               edge(:, ne + 3) = [v, piece_index(first(f), m, i, j + 1)]
               edge_length(ne + 1:ne + 3) = -1.0_wp
               ne = ne + 3
            end do
         end do
      end do
      ! Both face geometries must exist before evaluating a folded-edge path.
      do f = 1, nf
         m = divisions(f)
         ids = mesh%triangles(:, f)
         endpoints(:, 1) = ids([2, 3])
         endpoints(:, 2) = ids([1, 3])
         endpoints(:, 3) = ids([1, 2])
         do side = 1, 3
            g = mesh%neighbors(side, f)
            if (g <= f) cycle
            other = 0
            do l = 1, 3
               if (mesh%neighbors(l, g) == f) other = l
            end do
            if (other == 0) then
               call fatal_error(error, "MT patch graph has nonreciprocal face adjacency")
               return
            end if
            l = divisions(g)
            pair = endpoints(:, side)
            edge_vector = mesh%vertices(:, maxval(pair)) - mesh%vertices(:, minval(pair))
            select case (other)
            case (1)
               reverse_pair = mesh%triangles([2, 3], g)
            case (2)
               reverse_pair = mesh%triangles([1, 3], g)
            case default
               reverse_pair = mesh%triangles([1, 2], g)
            end select
            a = 0
            b = 0
            do while (a < m .and. b < l)
               ia = a
               ib = b
               if (pair(1) > pair(2)) ia = m - 1 - a
               if (reverse_pair(1) > reverse_pair(2)) ib = l - 1 - b
               u = boundary_piece(first(f), m, side, ia)
               v = boundary_piece(first(g), l, other, ib)
               ne = ne + 1
               edge(:, ne) = [u, v]
               lower = max(real(a, wp)/real(m, wp), real(b, wp)/real(l, wp))
               upper = min(real(a + 1, wp)/real(m, wp), real(b + 1, wp)/real(l, wp))
               start = mesh%vertices(:, minval(pair)) + lower*edge_vector
               finish = mesh%vertices(:, minval(pair)) + upper*edge_vector
               edge_length(ne) = edge_crossing_distance(graph%xyz(:, u), graph%xyz(:, v), start, finish)
               left_end = int(a + 1, int64)*int(l, int64)
               right_end = int(b + 1, int64)*int(m, int64)
               if (left_end <= right_end) a = a + 1
               if (right_end <= left_end) b = b + 1
            end do
         end do
      end do
      degree = 0
      do e = 1, ne
         degree(edge(:, e)) = degree(edge(:, e)) + 1
      end do
      graph%offset(1) = 1
      do u = 1, n
         graph%offset(u + 1) = graph%offset(u) + degree(u)
      end do
      allocate (graph%adjacent(2*ne), stat=stat)
      if (stat == 0) allocate (graph%length(2*ne), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT piece adjacency")
         return
      end if
      cursor = graph%offset(:n)
      do e = 1, ne
         u = edge(1, e)
         v = edge(2, e)
         graph%adjacent(cursor(u)) = v
         graph%adjacent(cursor(v)) = u
         ratio = edge_length(e)
         if (ratio < 0.0_wp) ratio = norm2(graph%xyz(:, u) - graph%xyz(:, v))
         graph%length(cursor(u)) = ratio
         graph%length(cursor(v)) = ratio
         cursor(u) = cursor(u) + 1
         cursor(v) = cursor(v) + 1
      end do
   end subroutine make_patch_graph

   !> Aggregate graph nodes into compact connected nodes for global balancing
   !>
   !> - Breadth-first growth from the lowest unassigned node bounds node area
   !> - Folded nodes stay apart, so each aggregate lies on one sheet
   !> - Edges are centroid chords, or across folds the shortest
   !>   centroid-node-node-centroid crossing
   !>
   !> @param[in] graph Node graph
   !> @param[in] normal Unit normal of each node (3, n)
   !> @param[in] max_area Largest aggregate area, bohr^2
   !> @param[out] coarse Aggregate graph without face data
   !> @param[out] coarse_normal Area-weighted unit normal of each aggregate
   !> @param[out] cluster Aggregate of each node
   !> @param[out] error Failure state
   subroutine coarsen_patch_graph(graph, normal, max_area, coarse, coarse_normal, cluster, error)
      !> Node graph
      type(patch_graph_type), intent(in) :: graph
      !> Node normals
      real(wp), intent(in) :: normal(:, :)
      !> Largest aggregate area, bohr^2
      real(wp), intent(in) :: max_area
      !> Aggregate graph
      type(patch_graph_type), intent(out) :: coarse
      !> Aggregate normals
      real(wp), allocatable, intent(out) :: coarse_normal(:, :)
      !> Aggregate of each node
      integer, allocatable, intent(out) :: cluster(:)
      !> Failure state
      type(error_type), allocatable, intent(out) :: error
      !> Smallest normal cosine joined into one aggregate
      real(wp), parameter :: min_cosine = 0.5_wp
      integer, allocatable :: queue(:), start(:), member(:), mark(:), slot(:), adjacent(:)
      real(wp), allocatable :: length(:)
      real(wp) :: total, seed_normal(3), candidate
      integer :: n, nn, s, u, v, c, d, j, k, head, tail, ne, stat

      n = size(graph%area)
      allocate (cluster(n), stat=stat)
      if (stat == 0) allocate (queue(n), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT piece aggregation")
         return
      end if
      cluster = 0
      nn = 0
      do s = 1, n
         if (cluster(s) /= 0) cycle
         nn = nn + 1
         cluster(s) = nn
         total = graph%area(s)
         seed_normal = normal(:, s)
         head = 1
         tail = 1
         queue(1) = s
         do while (head <= tail)
            u = queue(head)
            head = head + 1
            do k = graph%offset(u), graph%offset(u + 1) - 1
               v = graph%adjacent(k)
               if (cluster(v) /= 0) cycle
               if (total + graph%area(v) > max_area) cycle
               if (dot_product(seed_normal, normal(:, v)) < min_cosine) cycle
               cluster(v) = nn
               total = total + graph%area(v)
               tail = tail + 1
               queue(tail) = v
            end do
         end do
      end do
      deallocate (queue)
      allocate (coarse%area(nn), stat=stat)
      if (stat == 0) allocate (coarse%xyz(3, nn), stat=stat)
      if (stat == 0) allocate (coarse%offset(nn + 1), stat=stat)
      if (stat == 0) allocate (coarse_normal(3, nn), stat=stat)
      if (stat == 0) allocate (start(nn + 1), stat=stat)
      if (stat == 0) allocate (member(n), stat=stat)
      if (stat == 0) allocate (mark(nn), stat=stat)
      if (stat == 0) allocate (slot(nn), stat=stat)
      if (stat == 0) allocate (adjacent(graph%offset(n + 1) - 1), stat=stat)
      if (stat == 0) allocate (length(graph%offset(n + 1) - 1), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT aggregate graph")
         return
      end if
      coarse%area = 0.0_wp
      coarse%xyz = 0.0_wp
      coarse_normal = 0.0_wp
      mark = 0
      do u = 1, n
         c = cluster(u)
         coarse%area(c) = coarse%area(c) + graph%area(u)
         coarse%xyz(:, c) = coarse%xyz(:, c) + graph%area(u)*graph%xyz(:, u)
         coarse_normal(:, c) = coarse_normal(:, c) + graph%area(u)*normal(:, u)
         mark(c) = mark(c) + 1
      end do
      start(1) = 1
      do c = 1, nn
         coarse%xyz(:, c) = coarse%xyz(:, c)/coarse%area(c)
         total = norm2(coarse_normal(:, c))
         if (total > 0.0_wp) coarse_normal(:, c) = coarse_normal(:, c)/total
         start(c + 1) = start(c) + mark(c)
      end do
      slot = start(:nn)
      do u = 1, n
         c = cluster(u)
         member(slot(c)) = u
         slot(c) = slot(c) + 1
      end do
      mark = 0
      ne = 0
      do c = 1, nn
         coarse%offset(c) = ne + 1
         do j = start(c), start(c + 1) - 1
            u = member(j)
            do k = graph%offset(u), graph%offset(u + 1) - 1
               v = graph%adjacent(k)
               d = cluster(v)
               if (d == c) cycle
               ! Chords keep the metric isotropic; folds keep their intrinsic crossing
               if (dot_product(coarse_normal(:, c), coarse_normal(:, d)) >= min_cosine) then
                  candidate = norm2(coarse%xyz(:, d) - coarse%xyz(:, c))
               else
                  candidate = norm2(graph%xyz(:, u) - coarse%xyz(:, c)) + graph%length(k) &
                              + norm2(coarse%xyz(:, d) - graph%xyz(:, v))
               end if
               if (mark(d) /= c) then
                  mark(d) = c
                  ne = ne + 1
                  slot(d) = ne
                  adjacent(ne) = d
                  length(ne) = candidate
               else
                  length(slot(d)) = min(length(slot(d)), candidate)
               end if
            end do
         end do
      end do
      coarse%offset(nn + 1) = ne + 1
      coarse%adjacent = adjacent(:ne)
      coarse%length = length(:ne)
   end subroutine coarsen_patch_graph

   !> Shortest path between two face points crossing their shared edge segment
   !>
   !> Unfolding uses positive distances to the edge on opposite sides. A nearly
   !> closed fold therefore retains its intrinsic length instead of its chord.
   !>
   !> @param[in] a First face point
   !> @param[in] b Second face point
   !> @param[in] start Shared segment start
   !> @param[in] finish Shared segment end
   pure function edge_crossing_distance(a, b, start, finish) result(distance)
      !> First point
      real(wp), intent(in) :: a(3)
      !> Second point
      real(wp), intent(in) :: b(3)
      !> Segment start
      real(wp), intent(in) :: start(3)
      !> Segment end
      real(wp), intent(in) :: finish(3)
      !> Intrinsic path length
      real(wp) :: distance
      real(wp) :: tangent(3), qa(3), qb(3), crossing(3), length, xa, xb, ha, hb, x, fraction
      tangent = finish - start
      length = norm2(tangent)
      if (length <= 0.0_wp) then
         distance = norm2(a - start) + norm2(b - start)
         return
      end if
      tangent = tangent/length
      qa = a - start
      qb = b - start
      xa = dot_product(qa, tangent)
      xb = dot_product(qb, tangent)
      ha = norm2(qa - xa*tangent)
      hb = norm2(qb - xb*tangent)
      x = 0.5_wp*(xa + xb)
      if (ha + hb > 0.0_wp) then
         fraction = ha/(ha + hb)
         x = xa + fraction*(xb - xa)
      end if
      crossing = start + min(length, max(0.0_wp, x))*tangent
      distance = norm2(a - crossing) + norm2(b - crossing)
   end function edge_crossing_distance

   !> Closest point in a triangle centered at the origin
   !>
   !> Signed edge vectors describe both orientations of a subdivision piece.
   !> Scaling and cross-product barycentrics retain skinny positive-area faces.
   !>
   !> @param[in] point Query relative to the triangle centroid
   !> @param[in] edge_a First signed triangle edge
   !> @param[in] edge_b Second signed triangle edge
   !> @param[out] projected Closest point relative to the triangle centroid
   pure subroutine project_piece(point, edge_a, edge_b, projected)
      !> Relative query
      real(wp), intent(in) :: point(3)
      !> First edge
      real(wp), intent(in) :: edge_a(3)
      !> Second edge
      real(wp), intent(in) :: edge_b(3)
      !> Closest relative point
      real(wp), intent(out) :: projected(3)
      real(wp) :: a(3), b(3), q(3), corners(3, 3), normal(3), plane(3), r(3), cross(3)
      real(wp) :: edge(3), candidate(3), scale, length, denominator, x, y, t, distance, best
      integer :: i, j
      scale = max(norm2(edge_a), norm2(edge_b))
      a = edge_a/scale
      b = edge_b/scale
      q = point/scale
      corners(:, 1) = -(a + b)/3.0_wp
      corners(:, 2) = corners(:, 1) + a
      corners(:, 3) = corners(:, 1) + b
      normal = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
      denominator = norm2(normal)
      if (denominator > 0.0_wp) then
         normal = normal/denominator
         plane = q - dot_product(q - corners(:, 1), normal)*normal
         r = plane - corners(:, 1)
         cross = [r(2)*b(3) - r(3)*b(2), r(3)*b(1) - r(1)*b(3), r(1)*b(2) - r(2)*b(1)]
         x = dot_product(cross, normal)/denominator
         cross = [a(2)*r(3) - a(3)*r(2), a(3)*r(1) - a(1)*r(3), a(1)*r(2) - a(2)*r(1)]
         y = dot_product(cross, normal)/denominator
         if (x >= 0.0_wp .and. y >= 0.0_wp .and. x + y <= 1.0_wp) then
            projected = scale*plane
            return
         end if
      end if
      best = huge(1.0_wp)
      projected = 0.0_wp
      do i = 1, 3
         j = mod(i, 3) + 1
         edge = corners(:, j) - corners(:, i)
         length = norm2(edge)
         t = 0.0_wp
         if (length > 0.0_wp) t = min(1.0_wp, max(0.0_wp, dot_product(q - corners(:, i), edge/length)/length))
         candidate = corners(:, i) + t*edge
         distance = norm2(q - candidate)
         if (distance >= best) cycle
         best = distance
         projected = scale*candidate
      end do
   end subroutine project_piece

   !> Index of the upward triangle at integer barycentric coordinates
   !>
   !> @param[in] first First piece in the face
   !> @param[in] m Edge subdivision count
   !> @param[in] i First barycentric lattice index
   !> @param[in] j Second barycentric lattice index
   pure function piece_index(first, m, i, j) result(index)
      !> First piece
      integer, intent(in) :: first
      !> Subdivision count
      integer, intent(in) :: m
      !> First index
      integer, intent(in) :: i
      !> Second index
      integer, intent(in) :: j
      !> Piece index
      integer :: index
      index = first + i*(2*m - i) + 2*j
   end function piece_index

   !> Piece adjacent to an original face edge
   !>
   !> @param[in] first First piece in the face
   !> @param[in] m Edge subdivision count
   !> @param[in] side Edge opposite the corresponding triangle vertex
   !> @param[in] interval Interval in increasing local endpoint order
   pure function boundary_piece(first, m, side, interval) result(index)
      !> First piece
      integer, intent(in) :: first
      !> Subdivision count
      integer, intent(in) :: m
      !> Opposite vertex
      integer, intent(in) :: side
      !> Edge interval
      integer, intent(in) :: interval
      !> Adjacent piece
      integer :: index
      select case (side)
      case (1)
         index = piece_index(first, m, m - 1 - interval, interval)
      case (2)
         index = piece_index(first, m, 0, interval)
      case default
         index = piece_index(first, m, interval, 0)
      end select
   end function boundary_piece

   !> Additive Voronoi weights balancing cell areas on one graph level
   !>
   !> Gradient steps on the source offsets run until no pass reduces the count
   !> of cells outside their bounds for `patience` passes. The labels on return
   !> are the best assignment found; the weights produced it.
   !>
   !> @param[in] graph Node graph
   !> @param[in] seeds Distinct source nodes
   !> @param[in] goal Cell area goals, bohr^2
   !> @param[in] centers Cell area bound centers, bohr^2
   !> @param[in] tolerance Relative area tolerance
   !> @param[in] passes Largest number of assignments
   !> @param[in] patience Passes without fewer unbalanced cells before returning
   !> @param[in,out] weights Source distance offsets
   !> @param[out] labels Best connected cell assignment
   !> @param[out] distance Weighted source distances
   !> @param[out] error Allocation failure
   subroutine balance_weights(graph, seeds, goal, centers, tolerance, passes, patience, weights, labels, distance, &
         & error)
      !> Node graph
      type(patch_graph_type), intent(in) :: graph
      !> Source nodes
      integer, intent(in) :: seeds(:)
      !> Cell area goals
      real(wp), intent(in) :: goal(:)
      !> Bound centers
      real(wp), intent(in) :: centers(:)
      !> Relative tolerance
      real(wp), intent(in) :: tolerance
      !> Largest number of assignments
      integer, intent(in) :: passes
      !> Stall limit
      integer, intent(in) :: patience
      !> Source offsets
      real(wp), intent(inout) :: weights(:)
      !> Best assignment
      integer, intent(out) :: labels(:)
      !> Weighted distances
      real(wp), intent(out) :: distance(:)
      !> Allocation failure
      type(error_type), allocatable, intent(out) :: error
      !> Gradient step on the offsets, bohr per square-root bohr^2
      real(wp), parameter :: step = 0.3_wp
      type(pair_heap_type) :: heap
      integer, allocatable :: anchors(:), current(:)
      real(wp), allocatable :: mass(:), best_weights(:)
      real(wp) :: score, best_score
      integer :: n, np, p, u, iteration, outside, best_outside, stale, stat
      n = size(graph%area)
      np = size(seeds)
      allocate (anchors(n), stat=stat)
      if (stat == 0) allocate (current(n), stat=stat)
      if (stat == 0) allocate (mass(np), stat=stat)
      if (stat == 0) allocate (best_weights(np), stat=stat)
      if (stat == 0) allocate (heap%key(n + np), stat=stat)
      if (stat == 0) allocate (heap%node(n + np), stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot allocate MT Voronoi balancing")
         return
      end if
      anchors = 0
      do p = 1, np
         anchors(seeds(p)) = p
      end do
      best_outside = huge(1)
      best_score = huge(1.0_wp)
      best_weights = weights
      stale = 0
      ! Weighted shortest-path cells remain connected, including near thin sheets.
      ! A few cells oscillate long after the rest settle; diffusion repairs those.
      do iteration = 1, passes
         call assign_cells(graph, seeds, weights, anchors, heap, current, distance, stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot grow the MT Voronoi heap")
            return
         end if
         mass = 0.0_wp
         do u = 1, n
            mass(current(u)) = mass(current(u)) + graph%area(u)
         end do
         outside = count(abs(mass - centers) > tolerance*centers)
         score = maxval(abs(mass/goal - 1.0_wp))
         if (outside < best_outside .or. (outside == best_outside .and. score < best_score)) then
            if (outside < best_outside) stale = 0
            best_outside = outside
            best_score = score
            best_weights = weights
            labels = current
         else
            stale = stale + 1
         end if
         if (outside == 0 .or. stale >= patience) exit
         weights = weights + step*(mass - goal)/sqrt(goal)
         ! Subtracting a common offset leaves the partition unchanged.
         weights = weights - sum(weights)/real(np, wp)
      end do
      weights = best_weights
   end subroutine balance_weights

   !> Connected weighted graph Voronoi assignment
   !>
   !> @param[in] graph Intrinsic piece graph
   !> @param[in] seeds Distinct source pieces
   !> @param[in] weights Source distance offsets
   !> @param[in] anchors Protected source labels
   !> @param[in,out] heap Reusable work heap
   !> @param[out] labels Connected cell assignments
   !> @param[out] distance Weighted source distances
   !> @param[out] stat Nonzero when heap storage cannot grow
   subroutine assign_cells(graph, seeds, weights, anchors, heap, labels, distance, stat)
      !> Intrinsic graph
      type(patch_graph_type), intent(in) :: graph
      !> Source pieces
      integer, intent(in) :: seeds(:)
      !> Source offsets
      real(wp), intent(in) :: weights(:)
      !> Source labels
      integer, intent(in) :: anchors(:)
      !> Work heap
      type(pair_heap_type), intent(inout) :: heap
      !> Cell labels
      integer, intent(out) :: labels(:)
      !> Weighted distances
      real(wp), intent(out) :: distance(:)
      !> Allocation status
      integer, intent(out) :: stat
      real(wp) :: candidate, key
      integer :: p, u, v, k
      heap%count = 0
      labels = 0
      distance = huge(1.0_wp)
      do p = 1, size(seeds)
         u = seeds(p)
         labels(u) = p
         distance(u) = weights(p)
         call pair_push(heap, distance(u), u, stat)
         if (stat /= 0) return
      end do
      do while (heap%count > 0)
         call pair_pop(heap, key, u)
         ! Each improvement pushes a new entry; only the latest is current
         if (key > distance(u)) cycle
         p = labels(u)
         do k = graph%offset(u), graph%offset(u + 1) - 1
            v = graph%adjacent(k)
            if (anchors(v) /= 0 .and. anchors(v) /= p) cycle
            candidate = distance(u) + graph%length(k)
            if (candidate >= distance(v)) cycle
            distance(v) = candidate
            labels(v) = p
            call pair_push(heap, candidate, v, stat)
            if (stat /= 0) return
         end do
      end do
   end subroutine assign_cells

   !> Insert one entry, growing storage when full
   !>
   !> @param[in,out] heap Pair heap
   !> @param[in] key Priority
   !> @param[in] node Graph node
   !> @param[out] stat Nonzero when storage cannot grow
   subroutine pair_push(heap, key, node, stat)
      !> Pair heap
      type(pair_heap_type), intent(inout) :: heap
      !> Priority
      real(wp), intent(in) :: key
      !> Graph node
      integer, intent(in) :: node
      !> Allocation status
      integer, intent(out) :: stat
      real(wp), allocatable :: grown_key(:)
      integer, allocatable :: grown_node(:)
      integer :: at, parent
      stat = 0
      if (heap%count == size(heap%key)) then
         allocate (grown_key(2*size(heap%key)), stat=stat)
         if (stat == 0) allocate (grown_node(2*size(heap%key)), stat=stat)
         if (stat /= 0) return
         grown_key(:heap%count) = heap%key(:heap%count)
         grown_node(:heap%count) = heap%node(:heap%count)
         call move_alloc(grown_key, heap%key)
         call move_alloc(grown_node, heap%node)
      end if
      heap%count = heap%count + 1
      at = heap%count
      do while (at > 1)
         parent = at/2
         if (.not. pair_before(key, node, heap%key(parent), heap%node(parent))) exit
         heap%key(at) = heap%key(parent)
         heap%node(at) = heap%node(parent)
         at = parent
      end do
      heap%key(at) = key
      heap%node(at) = node
   end subroutine pair_push

   !> Remove the smallest entry
   !>
   !> @param[in,out] heap Nonempty pair heap
   !> @param[out] key Priority
   !> @param[out] node Graph node
   subroutine pair_pop(heap, key, node)
      !> Pair heap
      type(pair_heap_type), intent(inout) :: heap
      !> Priority
      real(wp), intent(out) :: key
      !> Graph node
      integer, intent(out) :: node
      real(wp) :: last_key
      integer :: last_node, at, child
      key = heap%key(1)
      node = heap%node(1)
      last_key = heap%key(heap%count)
      last_node = heap%node(heap%count)
      heap%count = heap%count - 1
      if (heap%count == 0) return
      at = 1
      do while (2*at <= heap%count)
         child = 2*at
         if (child < heap%count) then
            if (pair_before(heap%key(child + 1), heap%node(child + 1), heap%key(child), heap%node(child))) then
               child = child + 1
            end if
         end if
         if (.not. pair_before(heap%key(child), heap%node(child), last_key, last_node)) exit
         heap%key(at) = heap%key(child)
         heap%node(at) = heap%node(child)
         at = child
      end do
      heap%key(at) = last_key
      heap%node(at) = last_node
   end subroutine pair_pop

   !> Compare priorities with stable node tie breaking
   !>
   !> @param[in] a_key First priority
   !> @param[in] a_node First node
   !> @param[in] b_key Second priority
   !> @param[in] b_node Second node
   pure logical function pair_before(a_key, a_node, b_key, b_node) result(before)
      !> First priority
      real(wp), intent(in) :: a_key
      !> First node
      integer, intent(in) :: a_node
      !> Second priority
      real(wp), intent(in) :: b_key
      !> Second node
      integer, intent(in) :: b_node
      before = a_key < b_key .or. (.not. b_key < a_key .and. a_node < b_node)
   end function pair_before

   !> Check donor connectivity after removing one piece
   !>
   !> @param[in] graph Intrinsic graph
   !> @param[in] labels Current connected labels
   !> @param[in] removed Candidate piece
   !> @param[in,out] seen Reusable visitation stamps
   !> @param[in,out] stamp Current stamp
   !> @param[in,out] queue Reusable traversal storage
   function removable_piece(graph, labels, removed, seen, stamp, queue) result(ok)
      !> Intrinsic graph
      type(patch_graph_type), intent(in) :: graph
      !> Current labels
      integer, intent(in) :: labels(:)
      !> Candidate piece
      integer, intent(in) :: removed
      !> Visitation stamps
      integer, intent(inout) :: seen(:)
      !> Current stamp
      integer, intent(inout) :: stamp
      !> Traversal queue
      integer, intent(inout) :: queue(:)
      !> Donor remains connected
      logical :: ok
      integer :: u, v, k, head, tail, pending, donor
      if (stamp == huge(stamp)) then
         seen = 0
         stamp = 0
      end if
      stamp = stamp + 1
      donor = labels(removed)
      pending = 0
      u = 0
      do k = graph%offset(removed), graph%offset(removed + 1) - 1
         v = graph%adjacent(k)
         if (labels(v) /= donor) cycle
         pending = pending + 1
         u = v
      end do
      ok = pending <= 1
      if (ok) return
      head = 1
      tail = 1
      queue(1) = u
      seen(u) = stamp
      seen(removed) = stamp
      do while (head <= tail)
         u = queue(head)
         head = head + 1
         do k = graph%offset(u), graph%offset(u + 1) - 1
            v = graph%adjacent(k)
            if (v == removed) pending = pending - 1
            if (labels(v) /= donor .or. seen(v) == stamp) cycle
            seen(v) = stamp
            tail = tail + 1
            queue(tail) = v
         end do
         if (pending /= 0) cycle
         ok = .true.
         return
      end do
      ok = .false.
   end function removable_piece

   !> Empty a reusable heap
   !>
   !> @param[in,out] self Heap state
   subroutine heap_reset(self)
      !> Heap state
      class(index_heap_type), intent(inout) :: self
      do while (self%count > 0)
         self%position(self%node(self%count)) = 0
         self%count = self%count - 1
      end do
   end subroutine heap_reset

   !> Insert or update one priority
   !>
   !> @param[in,out] self Heap state
   !> @param[in] node Graph node
   !> @param[in] key Priority
   subroutine heap_set(self, node, key)
      !> Heap state
      class(index_heap_type), intent(inout) :: self
      !> Graph node
      integer, intent(in) :: node
      !> Priority
      real(wp), intent(in) :: key
      integer :: at, parent, child, other
      at = self%position(node)
      if (at == 0) then
         self%count = self%count + 1
         at = self%count
         self%node(at) = node
         self%position(node) = at
      end if
      self%key(node) = key
      do while (at > 1)
         parent = at/2
         other = self%node(parent)
         if (.not. heap_before(self, node, other)) exit
         self%node(at) = other
         self%position(other) = at
         at = parent
      end do
      self%node(at) = node
      self%position(node) = at
      do while (2*at <= self%count)
         child = 2*at
         if (child < self%count) then
            if (heap_before(self, self%node(child + 1), self%node(child))) child = child + 1
         end if
         other = self%node(child)
         if (.not. heap_before(self, other, node)) exit
         self%node(at) = other
         self%position(other) = at
         at = child
      end do
      self%node(at) = node
      self%position(node) = at
   end subroutine heap_set

   !> Remove the smallest entry
   !>
   !> @param[in,out] self Heap state
   !> @param[out] node Smallest graph node
   subroutine heap_pop(self, node)
      !> Heap state
      class(index_heap_type), intent(inout) :: self
      !> Removed node
      integer, intent(out) :: node
      integer :: last, at, child, other
      node = self%node(1)
      self%position(node) = 0
      last = self%node(self%count)
      self%count = self%count - 1
      if (self%count == 0) return
      at = 1
      do while (2*at <= self%count)
         child = 2*at
         if (child < self%count) then
            if (heap_before(self, self%node(child + 1), self%node(child))) child = child + 1
         end if
         other = self%node(child)
         if (.not. heap_before(self, other, last)) exit
         self%node(at) = other
         self%position(other) = at
         at = child
      end do
      self%node(at) = last
      self%position(last) = at
   end subroutine heap_pop

   !> Compare priorities with stable index tie breaking
   !>
   !> @param[in] self Heap state
   !> @param[in] a First node
   !> @param[in] b Second node
   pure function heap_before(self, a, b) result(before)
      !> Heap state
      class(index_heap_type), intent(in) :: self
      !> First node
      integer, intent(in) :: a
      !> Second node
      integer, intent(in) :: b
      !> First node precedes second
      logical :: before
      before = self%key(a) < self%key(b)
      if (.not. before .and. .not. (self%key(b) < self%key(a))) before = a < b
   end function heap_before

end module moist_cavity_mt_patches
