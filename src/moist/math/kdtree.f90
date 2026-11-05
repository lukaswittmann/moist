!> k-d tree over fixed points in 3D: the k nearest points of a query point
!>
!> - Built once from a copy of the points; median splits on the widest axis of each
!>   node's bounding box, leaves of at most `leaf_size` points unless their points
!>   coincide
!> - Queries are exact: the k nearest points ordered by distance, ties by point index,
!>   so the result does not depend on the build order of equal distances or on threads
!> - With a reach per point (`set_reach`), `within_reach` tells whether any point lies
!>   within its reach of a query point; nodes beyond the largest reach of their points are
!>   skipped by their bounding boxes, so far query points cost a few node tests
!> - The queries only read the tree; threads may query one tree at the same time
!>
!> Usage:
!> ```fortran
!> type(moist_kdtree_type) :: tree
!> call tree%build(xyz)                  ! xyz(3, n)
!> call tree%knn(point, k, nearest, dist2) ! point indices nearest(k), squared distances dist2(k)
!> call tree%set_reach(reach)            ! reach(n), bohr
!> inside = tree%within_reach(point)
!> call tree%destroy()
!> ```
module moist_math_kdtree
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: moist_kdtree_type

   !> Most points of a leaf
   integer, parameter :: leaf_size = 12

   !> k-d tree over fixed points
   type :: moist_kdtree_type
      !> Points, bohr (3, n)
      real(wp), allocatable :: x(:, :)
      !> Point indices, grouped by node
      integer, allocatable :: idx(:)
      !> First and last position in `idx` of each node
      integer, allocatable :: lo(:), hi(:)
      !> Children of each node, 0 for a leaf
      integer, allocatable :: left(:), right(:)
      !> Split axis of each inner node
      integer, allocatable :: axis(:)
      !> Split coordinate of each inner node
      real(wp), allocatable :: split(:)
      !> Bounding box of the points of each node, bohr (3, nnode)
      real(wp), allocatable :: lower(:, :), upper(:, :)
      !> Reach of each point, bohr; unallocated without `set_reach`
      real(wp), allocatable :: reach(:)
      !> Largest reach of the points of each node, bohr
      real(wp), allocatable :: node_reach(:)
      !> Nodes in use
      integer :: nnode = 0
      !> Levels below the root
      integer :: depth = 0
   contains
      !> Build the tree over a copy of the points
      procedure :: build => kdtree_build
      !> The k nearest points of a query point
      procedure :: knn => kdtree_knn
      !> Set the reach of each point
      procedure :: set_reach => kdtree_set_reach
      !> Whether a point lies within its reach of a query point
      procedure :: within_reach => kdtree_within_reach
      !> Release the tree
      procedure :: destroy => kdtree_destroy
   end type moist_kdtree_type

contains

   !> Build the tree over a copy of the points
   !>
   !> - A node of more than `leaf_size` points splits at the median of its widest axis;
   !>   a node whose points coincide along that axis stays a leaf
   !> - Every split leaves at least six points on each side, so there are at most n/6
   !>   leaves and fewer inner nodes
   !>
   !> @param[inout] self  tree
   !> @param[in]    x     points, bohr (3, n)
   subroutine kdtree_build(self, x)
      !> Tree
      class(moist_kdtree_type), intent(inout) :: self
      !> Points
      real(wp), intent(in) :: x(:, :)
      integer, allocatable :: level(:)
      real(wp) :: bmin(3), bmax(3)
      integer :: n, maxnode, i, l, h, m, d

      call self%destroy()
      n = size(x, 2)
      self%x = x
      allocate (self%idx(n))
      do i = 1, n
         self%idx(i) = i
      end do
      maxnode = 2*(n/6 + 1)
      allocate (self%lo(maxnode), self%hi(maxnode), self%left(maxnode), self%right(maxnode), &
         & self%axis(maxnode), self%split(maxnode), self%lower(3, maxnode), self%upper(3, maxnode), level(maxnode))
      self%nnode = 1
      self%lo(1) = 1
      self%hi(1) = n
      self%left(1) = 0
      self%right(1) = 0
      self%axis(1) = 1
      self%split(1) = 0.0_wp
      level(1) = 0
      self%depth = 0
      ! Nodes are split in creation order, so parents come before children
      i = 0
      do while (i < self%nnode)
         i = i + 1
         l = self%lo(i)
         h = self%hi(i)
         bmin = 0.0_wp
         bmax = 0.0_wp
         if (h >= l) then
            bmin = minval(self%x(:, self%idx(l:h)), dim=2)
            bmax = maxval(self%x(:, self%idx(l:h)), dim=2)
         end if
         self%lower(:, i) = bmin
         self%upper(:, i) = bmax
         if (h - l + 1 <= leaf_size) cycle
         d = maxloc(bmax - bmin, 1)
         if (.not. (bmax(d) - bmin(d) > 0.0_wp)) cycle
         m = (l + h)/2
         call select_kth(self%x, self%idx(l:h), m - l + 1, d)
         self%axis(i) = d
         self%split(i) = self%x(d, self%idx(m))
         self%left(i) = self%nnode + 1
         self%right(i) = self%nnode + 2
         self%lo(self%nnode + 1) = l
         self%hi(self%nnode + 1) = m
         self%lo(self%nnode + 2) = m + 1
         self%hi(self%nnode + 2) = h
         self%left(self%nnode + 1:self%nnode + 2) = 0
         self%right(self%nnode + 1:self%nnode + 2) = 0
         level(self%nnode + 1:self%nnode + 2) = level(i) + 1
         self%depth = max(self%depth, level(i) + 1)
         self%nnode = self%nnode + 2
      end do
   end subroutine kdtree_build

   !> Partial sort: idx(k) gets the k-th smallest coordinate along an axis, smaller or
   !> equal ones before it, larger or equal ones after it (Wirth's selection)
   !>
   !> @param[in]    x    points (3, n)
   !> @param[inout] idx  point indices
   !> @param[in]    k    rank to place
   !> @param[in]    d    axis
   pure subroutine select_kth(x, idx, k, d)
      !> Points
      real(wp), intent(in) :: x(:, :)
      !> Point indices
      integer, intent(inout) :: idx(:)
      !> Rank to place
      integer, intent(in) :: k
      !> Axis
      integer, intent(in) :: d
      integer :: l, r, i, j, t
      real(wp) :: pivot

      l = 1
      r = size(idx)
      do while (r > l)
         pivot = x(d, idx(k))
         i = l
         j = r
         do
            do while (x(d, idx(i)) < pivot)
               i = i + 1
            end do
            do while (pivot < x(d, idx(j)))
               j = j - 1
            end do
            if (i <= j) then
               t = idx(i)
               idx(i) = idx(j)
               idx(j) = t
               i = i + 1
               j = j - 1
            end if
            if (i > j) exit
         end do
         if (j < k) l = i
         if (k < i) r = j
      end do
   end subroutine select_kth

   !> The k nearest points of a query point
   !>
   !> - Ordered by squared distance, ties by point index
   !> - With fewer than k points the remaining entries are index 0 and distance `huge`
   !>
   !> @param[in]  self     built tree
   !> @param[in]  y        query point, bohr
   !> @param[in]  k        neighbours wanted
   !> @param[out] nearest  point indices (at least k)
   !> @param[out] dist2    squared distances, bohr^2 (at least k)
   pure subroutine kdtree_knn(self, y, k, nearest, dist2)
      !> Built tree
      class(moist_kdtree_type), intent(in) :: self
      !> Query point
      real(wp), intent(in) :: y(3)
      !> Neighbours wanted
      integer, intent(in) :: k
      !> Point indices
      integer, intent(out) :: nearest(:)
      !> Squared distances
      real(wp), intent(out) :: dist2(:)
      integer :: stack_node(2*self%depth + 2)
      real(wp) :: stack_bound(2*self%depth + 2)
      integer :: sp, node, p, j, nfound, m, near, far
      real(wp) :: bound, dd, diff

      nearest = 0
      dist2 = huge(1.0_wp)
      if (k < 1 .or. self%nnode < 1) return
      nfound = 0
      sp = 1
      stack_node(1) = 1
      stack_bound(1) = 0.0_wp
      do while (sp > 0)
         node = stack_node(sp)
         bound = stack_bound(sp)
         sp = sp - 1
         if (nfound == k) then
            if (bound > dist2(k)) cycle
         end if
         if (self%left(node) == 0) then
            do p = self%lo(node), self%hi(node)
               j = self%idx(p)
               dd = sum((self%x(:, j) - y)**2)
               if (nfound == k) then
                  if (dd > dist2(k)) cycle
                  if (dd == dist2(k) .and. j > nearest(k)) cycle
               else
                  nfound = nfound + 1
               end if
               ! Insert in order of distance, then index
               m = nfound
               do while (m > 1)
                  if (dist2(m - 1) < dd .or. (dist2(m - 1) == dd .and. nearest(m - 1) < j)) exit
                  dist2(m) = dist2(m - 1)
                  nearest(m) = nearest(m - 1)
                  m = m - 1
               end do
               dist2(m) = dd
               nearest(m) = j
            end do
         else
            ! The far side lies at least |diff| away along the split axis
            diff = y(self%axis(node)) - self%split(node)
            if (diff <= 0.0_wp) then
               near = self%left(node)
               far = self%right(node)
            else
               near = self%right(node)
               far = self%left(node)
            end if
            sp = sp + 1
            stack_node(sp) = far
            stack_bound(sp) = max(bound, diff*diff)
            sp = sp + 1
            stack_node(sp) = near
            stack_bound(sp) = bound
         end if
      end do
   end subroutine kdtree_knn

   !> Set the reach of each point
   !>
   !> @param[inout] self   built tree
   !> @param[in]    reach  reach of each point, bohr, in the order of the points (n)
   subroutine kdtree_set_reach(self, reach)
      !> Built tree
      class(moist_kdtree_type), intent(inout) :: self
      !> Reach of each point
      real(wp), intent(in) :: reach(:)
      integer :: i

      self%reach = reach
      if (allocated(self%node_reach)) deallocate (self%node_reach)
      allocate (self%node_reach(self%nnode))
      ! Children are created after their parents: from the last node back, every child is done
      do i = self%nnode, 1, -1
         if (self%left(i) == 0) then
            self%node_reach(i) = -1.0_wp
            if (self%hi(i) >= self%lo(i)) self%node_reach(i) = maxval(self%reach(self%idx(self%lo(i):self%hi(i))))
         else
            self%node_reach(i) = max(self%node_reach(self%left(i)), self%node_reach(self%right(i)))
         end if
      end do
   end subroutine kdtree_set_reach

   !> Whether a point lies within its reach of a query point
   !>
   !> - Needs `set_reach`; false without it
   !>
   !> @param[in] self  built tree with reaches
   !> @param[in] y     query point, bohr
   pure logical function kdtree_within_reach(self, y) result(inside)
      !> Built tree with reaches
      class(moist_kdtree_type), intent(in) :: self
      !> Query point
      real(wp), intent(in) :: y(3)
      integer :: stack_node(2*self%depth + 2)
      integer :: sp, node, p, j
      real(wp) :: gap(3)

      inside = .false.
      if (.not. allocated(self%node_reach) .or. self%nnode < 1) return
      sp = 1
      stack_node(1) = 1
      do while (sp > 0)
         node = stack_node(sp)
         sp = sp - 1
         if (self%node_reach(node) < 0.0_wp) cycle
         ! Distance from the query point to the bounding box of the node
         gap = max(0.0_wp, self%lower(:, node) - y, y - self%upper(:, node))
         if (sum(gap**2) > self%node_reach(node)**2) cycle
         if (self%left(node) == 0) then
            do p = self%lo(node), self%hi(node)
               j = self%idx(p)
               if (sum((self%x(:, j) - y)**2) <= self%reach(j)**2) then
                  inside = .true.
                  return
               end if
            end do
         else
            sp = sp + 1
            stack_node(sp) = self%left(node)
            sp = sp + 1
            stack_node(sp) = self%right(node)
         end if
      end do
   end function kdtree_within_reach

   !> Release the tree
   !>
   !> @param[inout] self  tree
   subroutine kdtree_destroy(self)
      !> Tree
      class(moist_kdtree_type), intent(inout) :: self

      if (allocated(self%x)) deallocate (self%x)
      if (allocated(self%idx)) deallocate (self%idx)
      if (allocated(self%lo)) deallocate (self%lo)
      if (allocated(self%hi)) deallocate (self%hi)
      if (allocated(self%left)) deallocate (self%left)
      if (allocated(self%right)) deallocate (self%right)
      if (allocated(self%axis)) deallocate (self%axis)
      if (allocated(self%split)) deallocate (self%split)
      if (allocated(self%lower)) deallocate (self%lower)
      if (allocated(self%upper)) deallocate (self%upper)
      if (allocated(self%reach)) deallocate (self%reach)
      if (allocated(self%node_reach)) deallocate (self%node_reach)
      self%nnode = 0
      self%depth = 0
   end subroutine kdtree_destroy

end module moist_math_kdtree
