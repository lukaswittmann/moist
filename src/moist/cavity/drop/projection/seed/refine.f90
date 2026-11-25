!> Surface mesh build, refinement: cell splitting, sign-change crossings
submodule(moist_cavity_drop_projection_seed) moist_cavity_drop_projection_seed_refine
   implicit none(type, external)

contains

   !* ================================================================================= *!
   !*                                    Refinement                                     *!
   !* ================================================================================= *!

   !> Split cells level by level and collect the crossings of the unsplit ones
   !>
   !> - level-0 cells: every cube whose eight corners are known nodes
   !> - below `h` only cells without sign change are split
   !> - counted in `n_unresolved`: a sign-change cell without a reached sign-change
   !>   edge; a deepest cell without sign change seeds nothing
   !>
   !> @param[inout] self  Seeder (tallies)
   !> @param[in]    lsf   Level set (certificates)
   !> @param[inout] work  Build scratch
   !> @param[out]   error Level set failure
   module subroutine refine_cells(self, lsf, work, error)
      type(projection_seeder_type), intent(inout) :: self
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      type(build_type), intent(inout) :: work
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: cells(:, :), corners0(:, :), next(:, :)
      type(key_table_type) :: edges
      integer :: ncell, nnext, s, k, c, step, half_step, corner(8), base(3), ic, jc, kc
      real(wp) :: hs
      logical :: change, certified, split, crossed

      call level0_cells(work, cells, corners0, ncell)
      call table_init(edges, 4096)
      allocate (work%cross(3, 1024), work%cross_h(1024), work%cross_mid(3, 1024), work%cross_axis(1024))
      work%ncross = 0
      self%n_unresolved = 0

      do s = 0, work%depth
         step = level_step(work, s)
         half_step = step/2
         hs = work%h0/2.0_wp**s
         nnext = 0
         if (s < work%depth) then
            if (allocated(next)) deallocate (next)
            allocate (next(3, max(64, ncell)))
         end if
         do k = 1, ncell
            if (s == 0) then
               corner = corners0(:, k)
            else
               call cell_corners(work, cells(:, k), step, corner)
            end if
            call classify_cell(work, corner, hs, change, certified)
            if (s < work%m) then
               split = change .or. .not. certified
            else
               split = s < work%depth .and. .not. change .and. .not. certified
               if (change) then
                  call cell_crossings(work, edges, corner, hs, crossed)
                  if (.not. crossed) self%n_unresolved = self%n_unresolved + 1
               end if
            end if
            if (.not. split .or. s == work%depth) cycle
            ! 19 new child nodes (the 8 parent corners are already in `work%nodes`)
            base = 2*cells(:, k)
            do kc = 0, 2
               do jc = 0, 2
                  do ic = 0, 2
                     if (iand(ic, 1) == 0 .and. iand(jc, 1) == 0 .and. iand(kc, 1) == 0) cycle
                     call add_node(work%nodes, half_step*(base + [ic, jc, kc]))
                  end do
               end do
            end do
            if (nnext + 8 > size(next, 2)) call grow_array(next, 3, max(nnext + 8, 2*size(next, 2)))
            do c = 0, 7
               nnext = nnext + 1
               next(:, nnext) = base + [iand(c, 1), iand(ishft(c, -1), 1), iand(ishft(c, -2), 1)]
            end do
         end do
         if (s == 0 .and. allocated(corners0)) deallocate (corners0)
         if (s == work%depth .or. nnext == 0) exit
         call evaluate_nodes(self, lsf, work, error)
         if (allocated(error)) return
         call move_alloc(next, cells)
         ncell = nnext
      end do
      self%n_cross = work%ncross
   end subroutine refine_cells

   !> Level-0 cells: cubes whose eight corners are known level-0 nodes, by min corner
   !>
   !> @param[in]  work     Build scratch
   !> @param[out] cells    Min corners in level-0 coordinates (3, ncell)
   !> @param[out] corners0 Node indices of each cell's eight corners (8, ncell)
   !> @param[out] ncell    Cells
   subroutine level0_cells(work, cells, corners0, ncell)
      type(build_type), intent(in) :: work
      integer, allocatable, intent(out) :: cells(:, :), corners0(:, :)
      integer, intent(out) :: ncell

      integer :: k, step, j(3), corner(8)

      step = level_step(work, 0)
      allocate (cells(3, max(1, work%nodes%n)), corners0(8, max(1, work%nodes%n)))
      ncell = 0
      do k = 1, work%nodes%n
         if (any(modulo(work%nodes%ijk(:, k), step) /= 0)) cycle
         j = work%nodes%ijk(:, k)/step
         call cell_corners(work, j, step, corner)
         if (any(corner == 0)) cycle
         ncell = ncell + 1
         cells(:, ncell) = j
         corners0(:, ncell) = corner
      end do
   end subroutine level0_cells

   !> Node indices of a cell's corners, 0 for an unknown node
   !>
   !> @param[in]  work   Build scratch
   !> @param[in]  j      Min corner in level coordinates
   !> @param[in]  step   Level step in finest units
   !> @param[out] corner Corner nodes, corner `1 + i + 2 j + 4 k`
   subroutine cell_corners(work, j, step, corner)
      type(build_type), intent(in) :: work
      integer, intent(in) :: j(3)
      integer, intent(in) :: step
      integer, intent(out) :: corner(8)

      integer :: c

      do c = 0, 7
         corner(c + 1) = table_find(work%nodes%table, &
                                    lattice_key(step*(j + [iand(c, 1), iand(ishft(c, -1), 1), &
                                                           iand(ishft(c, -2), 1)]), key_bias))
      end do
   end subroutine cell_corners

   !> Sign change and surface-free certificate of a cell
   !>
   !> - an unreached corner counts as outside for the sign test
   !> - checks exact 1-corner, opposite asymmetric face-pair (subsuming 4-corner
   !>   face and 8-corner rules), and 2-corner edge coverings
   !>
   !> @param[in]  work      Build scratch
   !> @param[in]  corner    Corner nodes
   !> @param[in]  hs        Cell edge (Bohr)
   !> @param[out] change    Corners on both sides of the surface
   !> @param[out] certified Corner surface-free balls cover the cell
   subroutine classify_cell(work, corner, hs, change, certified)
      type(build_type), intent(in) :: work
      integer, intent(in) :: corner(8)
      real(wp), intent(in) :: hs
      logical, intent(out) :: change
      logical, intent(out) :: certified

      integer :: c, k, e
      real(wp) :: rf(8), diag, hs_safe, half_h2, two_h2, r0, r1
      logical :: inside, outside

      inside = .false.
      outside = .false.
      certified = .false.
      diag = sqrt(3.0_wp)*hs
      do c = 1, 8
         k = corner(c)
         if (work%nodes%reached(k) .and. work%nodes%val(k) <= 0.0_wp) then
            inside = .true.
         else
            outside = .true.
         end if
         rf(c) = work%nodes%free(k)
         if (rf(c) >= diag) certified = .true.
      end do
      change = inside .and. outside
      if (change) then
         certified = .false.
         return
      end if
      if (certified) return

      ! Exact multi-corner coverings (all `rf < sqrt(3) hs`, so squares cannot overflow)
      hs_safe = (1.0_wp + 1.0e-12_wp)*hs
      half_h2 = 0.5_wp*hs*hs
      r0 = min(rf(1), rf(3), rf(5), rf(7))
      r1 = min(rf(2), rf(4), rf(6), rf(8))
      if (sqrt(max(0.0_wp, r0*r0 - half_h2)) + sqrt(max(0.0_wp, r1*r1 - half_h2)) >= hs_safe) then
         certified = .true.
         return
      end if
      r0 = min(rf(1), rf(2), rf(5), rf(6))
      r1 = min(rf(3), rf(4), rf(7), rf(8))
      if (sqrt(max(0.0_wp, r0*r0 - half_h2)) + sqrt(max(0.0_wp, r1*r1 - half_h2)) >= hs_safe) then
         certified = .true.
         return
      end if
      r0 = min(rf(1), rf(2), rf(3), rf(4))
      r1 = min(rf(5), rf(6), rf(7), rf(8))
      if (sqrt(max(0.0_wp, r0*r0 - half_h2)) + sqrt(max(0.0_wp, r1*r1 - half_h2)) >= hs_safe) then
         certified = .true.
         return
      end if
      two_h2 = 2.0_wp*hs*hs
      do e = 1, 12
         r0 = rf(edge_corner(1, e))
         r1 = rf(edge_corner(2, e))
         if (sqrt(max(0.0_wp, r0*r0 - two_h2)) + sqrt(max(0.0_wp, r1*r1 - two_h2)) >= hs_safe) then
            certified = .true.
            return
         end if
      end do
   end subroutine classify_cell

   !> Crossings of a cell's sign-change edges, each edge once
   !>
   !> - an edge to an unreached corner has no crossing
   !>
   !> @param[inout] work    Build scratch
   !> @param[inout] edges   Edges already seen, by doubled midpoint
   !> @param[in]    corner  Corner nodes
   !> @param[in]    hs      Edge length (Bohr)
   !> @param[out]   crossed Whether the cell has a sign-change edge between reached corners
   subroutine cell_crossings(work, edges, corner, hs, crossed)
      type(build_type), intent(inout) :: work
      type(key_table_type), intent(inout) :: edges
      integer, intent(in) :: corner(8)
      real(wp), intent(in) :: hs
      logical, intent(out) :: crossed

      integer :: e, a, b
      integer(int64) :: key
      real(wp) :: va, vb, xa(3), xb(3)

      crossed = .false.
      do e = 1, 12
         a = corner(edge_corner(1, e))
         b = corner(edge_corner(2, e))
         if (.not. (work%nodes%reached(a) .and. work%nodes%reached(b))) cycle
         va = work%nodes%val(a)
         vb = work%nodes%val(b)
         if ((va > 0.0_wp) .eqv. (vb > 0.0_wp)) cycle
         crossed = .true.
         ! Doubled midpoint: unique per edge over all levels, in range with the nodes
         key = lattice_key(work%nodes%ijk(:, a) + work%nodes%ijk(:, b), key_bias)
         if (table_find(edges, key) > 0) cycle
         call table_insert(edges, key, 1)
         xa = node_position(work, work%nodes%ijk(:, a))
         xb = node_position(work, work%nodes%ijk(:, b))
         call push_crossing(work, segment_crossing(xa, xb, va, vb), hs, 0.5_wp*(xa + xb), (e + 3)/4)
      end do
   end subroutine cell_crossings

   !> Append one crossing
   !>
   !> @param[inout] work Build scratch
   !> @param[in]    x    Crossing
   !> @param[in]    hs   Edge length (Bohr)
   !> @param[in]    mid  Edge midpoint
   !> @param[in]    axis Edge axis
   subroutine push_crossing(work, x, hs, mid, axis)
      type(build_type), intent(inout) :: work
      real(wp), intent(in) :: x(3)
      real(wp), intent(in) :: hs
      real(wp), intent(in) :: mid(3)
      integer, intent(in) :: axis

      integer :: n

      n = work%ncross
      if (n + 1 > size(work%cross_h)) then
         call grow_array(work%cross, 3, 2*size(work%cross_h))
         call grow_array(work%cross_mid, 3, 2*size(work%cross_h))
         call grow_array(work%cross_axis, 2*size(work%cross_h))
         call grow_array(work%cross_h, 2*size(work%cross_h))
      end if
      work%ncross = n + 1
      work%cross(:, n + 1) = x
      work%cross_h(n + 1) = hs
      work%cross_mid(:, n + 1) = mid
      work%cross_axis(n + 1) = axis
   end subroutine push_crossing

end submodule moist_cavity_drop_projection_seed_refine
