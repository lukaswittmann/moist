!> DROP anchor grid: the unprojected Lebedev points on the atomic spheres
!>
!> Struct of arrays indexed by anchor, `1 .. ngrid`, never by surface row;
!> the projection maps anchors to rows through the cavity's `offset` and
!> `anchor_of_row`
!>
!> - `fill` places the nodes of every sphere's grid components
!>   ([[cavity_grid_type]]): the plain Lebedev points, or those of a grid
!>   model (rotated per component)
!> - `compute_switching` evaluates the iSwiG widths and switching factors,
!>   scaled by the component weights
!> - `filter` drops anchors whose weight `wleb0*f0` is at or below the cutoff
!> - `restrict` keeps only listed anchor ids
module moist_cavity_drop_anchor
   use mctc_env_accuracy, only: wp
   use mctc_env, only: error_type
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig
   use moist_cavity_grid_type, only: cavity_grid_type
   use moist_utils_mem, only: filter_array
   implicit none(type, external)
   private

   public :: drop_anchor_type

   !> Anchor points of the DROP cavity, one entry per anchor
   type :: drop_anchor_type
      !> Number of anchor points currently stored
      integer :: ngrid = 0
      !> Anchor id before the pre-filter, the stable anchor identity (ngrid)
      integer, allocatable :: id(:)
      !> Owning atom of each anchor (ngrid)
      integer, allocatable :: owner(:)
      !> Grid component of each anchor (ngrid)
      integer, allocatable :: comp(:)
      !> Anchor (unprojected) positions (3, ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Anchor Lebedev quadrature weights (ngrid)
      real(wp), allocatable :: wleb0(:)
      !> Gaussian width at each anchor (ngrid)
      real(wp), allocatable :: xi0(:)
      !> Anchor switching factor (ngrid)
      real(wp), allocatable :: f0(:)
   contains
      !> Place the plain Lebedev points on every sphere
      procedure :: fill => anchor_fill
      !> Place the nodes of every sphere's grid components
      procedure :: fill_grid => anchor_fill_grid
      !> Evaluate the iSwiG widths and switching factors
      procedure :: compute_switching => anchor_compute_switching
      !> Drop anchors at or below the weight cutoff (pre-filter)
      procedure :: filter => anchor_filter
      !> Keep only the anchors with listed ids
      procedure :: restrict => anchor_restrict
      !> Deallocate every array and reset the count
      procedure :: clear => anchor_clear
   end type drop_anchor_type

contains

   !> Place the nodes of every sphere's grid components
   !>
   !> - anchor `ii` is global node `ii` of the grid, and that is its id:
   !>   node `j` of component `c` is anchor `comp_first(c) + j - 1`, at
   !>   `c_i + r_i u_j` with `u_j` the node of the component's node set
   !> - the plain grid gives anchor `(i - 1)*num_leb + j` for Lebedev point
   !>   `j` on sphere `i`, bitwise at `c_i + r_i u_j`
   !> - widths and switching factors of a previous grid are dropped
   !>
   !> @param[inout] self    anchor grid
   !> @param[in]    centers sphere centers, bohr (3, nsph)
   !> @param[in]    radii   sphere radii, bohr (nsph)
   !> @param[in]    grid    grid components of every sphere
   subroutine anchor_fill_grid(self, centers, radii, grid)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self
      !> Sphere centers, bohr (3, nsph)
      real(wp), intent(in) :: centers(:, :)
      !> Sphere radii, bohr (nsph)
      real(wp), intent(in) :: radii(:)
      !> Grid components of every sphere
      type(cavity_grid_type), intent(in) :: grid

      !> Sphere, component, node and anchor index, node offset
      integer :: i, icomp, jj, ii, n0
      !> Number of spheres
      integer :: nsph

      nsph = size(radii)
      self%ngrid = grid%nnode

      if (allocated(self%xyz)) deallocate (self%xyz)
      allocate (self%xyz(3, self%ngrid), source=0.0_wp)
      if (allocated(self%wleb0)) deallocate (self%wleb0)
      allocate (self%wleb0(self%ngrid), source=0.0_wp)
      if (allocated(self%owner)) deallocate (self%owner)
      allocate (self%owner(self%ngrid), source=0)
      if (allocated(self%comp)) deallocate (self%comp)
      allocate (self%comp(self%ngrid), source=0)
      if (allocated(self%id)) deallocate (self%id)
      allocate (self%id(self%ngrid), source=0)
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%f0)) deallocate (self%f0)

      do i = 1, nsph
         do icomp = grid%sph_first(i), grid%sph_first(i + 1) - 1
            n0 = grid%set_first(grid%comp_set(icomp)) - grid%comp_first(icomp)
            do ii = grid%comp_first(icomp), grid%comp_first(icomp + 1) - 1
               jj = n0 + ii

               ! Solid-angle Lebedev weight, the same for a rotated node
               self%wleb0(ii) = grid%wang(jj)

               ! Cartesian location of the node on sphere i
               self%xyz(1, ii) = centers(1, i) + radii(i)*grid%node(1, jj)
               self%xyz(2, ii) = centers(2, i) + radii(i)*grid%node(2, jj)
               self%xyz(3, ii) = centers(3, i) + radii(i)*grid%node(3, jj)

               self%owner(ii) = i
               self%comp(ii) = icomp
               self%id(ii) = ii
            end do
         end do
      end do

   end subroutine anchor_fill_grid

   !> Place the plain Lebedev points on every sphere
   !>
   !> [[anchor_fill_grid]] for one unrotated component of weight one per
   !> sphere: anchor `(i - 1)*num_leb + j` is point `j` on sphere `i`
   !>
   !> @param[inout] self       anchor grid
   !> @param[in]    centers    sphere centers, bohr (3, nsph)
   !> @param[in]    radii      sphere radii, bohr (nsph)
   !> @param[in]    ang_grid   unit-sphere Lebedev nodes (3, num_leb)
   !> @param[in]    ang_weight solid-angle Lebedev weights (num_leb)
   subroutine anchor_fill(self, centers, radii, ang_grid, ang_weight)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self
      !> Sphere centers, bohr (3, nsph)
      real(wp), intent(in) :: centers(:, :)
      !> Sphere radii, bohr (nsph)
      real(wp), intent(in) :: radii(:)
      !> Unit-sphere Lebedev nodes (3, num_leb)
      real(wp), intent(in) :: ang_grid(:, :)
      !> Solid-angle Lebedev weights (num_leb)
      real(wp), intent(in) :: ang_weight(:)

      !> Plain grid components
      type(cavity_grid_type) :: grid
      !> Never set: a plain grid cannot fail
      type(error_type), allocatable :: error

      call grid%assemble(size(radii), ang_grid, ang_weight, 0.0_wp, error)
      call self%fill_grid(centers, radii, grid)

   end subroutine anchor_fill

   !> Evaluate the iSwiG width and switching factor of every anchor
   !>
   !> - with `comp_weight`, the switching factor of an anchor is scaled by
   !>   the weight of its grid component
   !>
   !> @param[inout] self        anchor grid after `fill`
   !> @param[in]    iswig       iSwiG updated for the current geometry
   !> @param[in]    nthreads    OpenMP team size
   !> @param[in]    comp_weight weight of every grid component (ncomp)
   subroutine anchor_compute_switching(self, iswig, nthreads, comp_weight)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self
      !> iSwiG updated for the current geometry
      type(moist_cavity_drop_iswig), intent(in) :: iswig
      !> OpenMP team size
      integer, intent(in) :: nthreads
      !> Weight of every grid component
      real(wp), intent(in), optional :: comp_weight(:)

      !> Anchor index
      integer :: i

      if (allocated(self%xi0)) deallocate (self%xi0)
      allocate (self%xi0(self%ngrid), source=0.0_wp)
      if (allocated(self%f0)) deallocate (self%f0)
      allocate (self%f0(self%ngrid), source=1.0_wp)

      !$omp parallel num_threads(nthreads) default(shared) private(i)
      !$omp do schedule(static)
      do i = 1, self%ngrid
         ! iSwiG (sorted neighbour list with early exit)
         self%xi0(i) = iswig%xi0(self%owner(i), self%wleb0(i))
         self%f0(i) = iswig%swi0(self%xyz(:, i), self%owner(i), self%xi0(i))
      end do
      !$omp end do
      !$omp end parallel

      if (present(comp_weight)) then
         do i = 1, self%ngrid
            self%f0(i) = self%f0(i)*comp_weight(self%comp(i))
         end do
      end if

   end subroutine anchor_compute_switching

   !> Drop anchors whose weight `wleb0*f0` is at or below the cutoff
   !>
   !> - `wleb0` and `f0` must be allocated
   !> - `xyz`, `xi0`, `owner` and `id` are compacted when allocated
   !>
   !> @param[inout] self     anchor grid
   !> @param[in]    wleb_cut weight cutoff
   !> @param[out]   nremoved number of anchors dropped
   subroutine anchor_filter(self, wleb_cut, nremoved)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self
      !> Weight cutoff
      real(wp), intent(in) :: wleb_cut
      !> Number of anchors dropped
      integer, intent(out) :: nremoved

      !> Anchor count before and after the filter
      integer :: nold, nvalid
      !> Keep mask (nold)
      logical, allocatable :: keep(:)

      nold = self%ngrid
      allocate (keep(nold), source=.false.)
      keep = self%wleb0*self%f0 > wleb_cut

      nvalid = count(keep)
      nremoved = nold - nvalid

      call filter_array(self%xyz, nold, keep, nvalid)
      call filter_array(self%wleb0, nold, keep, nvalid)
      call filter_array(self%f0, nold, keep, nvalid)
      call filter_array(self%xi0, nold, keep, nvalid)
      call filter_array(self%owner, nold, keep, nvalid)
      call filter_array(self%comp, nold, keep, nvalid)
      call filter_array(self%id, nold, keep, nvalid)
      self%ngrid = nvalid

   end subroutine anchor_filter

   !> Keep only the anchors whose id is listed (mainly for unittests)
   !>
   !> - ids not on the grid are ignored, so one list serves every geometry
   !> - `xi0` and `f0` must be allocated; the other arrays as in `filter`
   !>
   !> @param[inout] self     anchor grid
   !> @param[in]    ids      anchor ids to keep
   !> @param[out]   nremoved number of anchors dropped
   subroutine anchor_restrict(self, ids, nremoved)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self
      !> Anchor ids to keep
      integer, intent(in) :: ids(:)
      !> Number of anchors dropped
      integer, intent(out) :: nremoved

      !> Anchor count before and after
      integer :: nold, nvalid
      !> Listed per id, and keep mask (nold)
      logical, allocatable :: listed(:), keep(:)
      integer :: i

      nold = self%ngrid
      ! maxval of no anchors is -huge
      allocate (listed(max(0, maxval(self%id(1:nold)))), source=.false.)
      do i = 1, size(ids)
         if (ids(i) >= 1 .and. ids(i) <= size(listed)) listed(ids(i)) = .true.
      end do
      allocate (keep(nold))
      do i = 1, nold
         keep(i) = listed(self%id(i))
      end do

      nvalid = count(keep)
      nremoved = nold - nvalid

      call filter_array(self%xyz, nold, keep, nvalid)
      call filter_array(self%wleb0, nold, keep, nvalid)
      call filter_array(self%f0, nold, keep, nvalid)
      call filter_array(self%xi0, nold, keep, nvalid)
      call filter_array(self%owner, nold, keep, nvalid)
      call filter_array(self%comp, nold, keep, nvalid)
      call filter_array(self%id, nold, keep, nvalid)
      self%ngrid = nvalid

   end subroutine anchor_restrict

   !> Deallocate every array and reset the count
   !>
   !> @param[inout] self anchor grid
   subroutine anchor_clear(self)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self

      self%ngrid = 0
      if (allocated(self%id)) deallocate (self%id)
      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%comp)) deallocate (self%comp)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%wleb0)) deallocate (self%wleb0)
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%f0)) deallocate (self%f0)

   end subroutine anchor_clear

end module moist_cavity_drop_anchor
