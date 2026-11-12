!> DROP anchor grid: the reference grid
!>
!> Struct of arrays indexed by anchor, `1 .. ngrid`, never by surface row;
!> the projection maps anchors to rows through the cavity's `offset` and
!> `anchor_of_row`
!>
!> - `fill` places the Lebedev points on every sphere
!> - `compute_switching` evaluates the iSwiG widths and switching factors
!> - `filter` drops anchors whose weight `wleb0*f0` is at or below the cutoff
module moist_cavity_drop_anchor
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_gaussian, only: moist_cavity_drop_iswig
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
      !> Anchor (unprojected) positions (3, ngrid)
      real(wp), allocatable :: xyz(:, :)
      !> Anchor Lebedev quadrature weights (ngrid)
      real(wp), allocatable :: wleb0(:)
      !> Gaussian width at each anchor (ngrid)
      real(wp), allocatable :: xi0(:)
      !> Anchor switching factor (ngrid)
      real(wp), allocatable :: f0(:)
   contains
      !> Place the Lebedev points on every sphere
      procedure :: fill => anchor_fill
      !> Evaluate the iSwiG widths and switching factors
      procedure :: compute_switching => anchor_compute_switching
      !> Drop anchors at or below the weight cutoff (pre-filter)
      procedure :: filter => anchor_filter
      !> Deallocate every array and reset the count
      procedure :: clear => anchor_clear
   end type drop_anchor_type

contains

   !> Place the Lebedev points on every sphere
   !>
   !> - anchor `(i - 1)*num_leb + j` is point `j` on sphere `i`, and its id
   !> - widths and switching factors of a previous grid are dropped
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

      !> Sphere, Lebedev point and anchor index
      integer :: i, jj, ii
      !> Number of spheres and Lebedev points per sphere
      integer :: nsph, num_leb

      nsph = size(radii)
      num_leb = size(ang_weight)
      self%ngrid = nsph*num_leb

      if (allocated(self%xyz)) deallocate (self%xyz)
      allocate (self%xyz(3, self%ngrid), source=0.0_wp)
      if (allocated(self%wleb0)) deallocate (self%wleb0)
      allocate (self%wleb0(self%ngrid), source=0.0_wp)
      if (allocated(self%owner)) deallocate (self%owner)
      allocate (self%owner(self%ngrid), source=0)
      if (allocated(self%id)) deallocate (self%id)
      allocate (self%id(self%ngrid), source=0)
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%f0)) deallocate (self%f0)

      ii = 0
      do i = 1, nsph
         do jj = 1, num_leb
            ii = ii + 1

            ! Solid-angle Lebedev weight
            self%wleb0(ii) = ang_weight(jj)

            ! Cartesian location of point on sphere i
            self%xyz(1, ii) = centers(1, i) + radii(i)*ang_grid(1, jj)
            self%xyz(2, ii) = centers(2, i) + radii(i)*ang_grid(2, jj)
            self%xyz(3, ii) = centers(3, i) + radii(i)*ang_grid(3, jj)

            self%owner(ii) = i
            self%id(ii) = ii
         end do
      end do

   end subroutine anchor_fill

   !> Evaluate the iSwiG width and switching factor of every anchor
   !>
   !> @param[inout] self     anchor grid after `fill`
   !> @param[in]    iswig    iSwiG updated for the current geometry
   !> @param[in]    nthreads OpenMP team size
   subroutine anchor_compute_switching(self, iswig, nthreads)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self
      !> iSwiG updated for the current geometry
      type(moist_cavity_drop_iswig), intent(in) :: iswig
      !> OpenMP team size
      integer, intent(in) :: nthreads

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
      call filter_array(self%id, nold, keep, nvalid)
      self%ngrid = nvalid

   end subroutine anchor_filter

   !> Deallocate every array and reset the count
   !>
   !> @param[inout] self anchor grid
   subroutine anchor_clear(self)
      !> Anchor grid
      class(drop_anchor_type), intent(inout) :: self

      self%ngrid = 0
      if (allocated(self%id)) deallocate (self%id)
      if (allocated(self%owner)) deallocate (self%owner)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%wleb0)) deallocate (self%wleb0)
      if (allocated(self%xi0)) deallocate (self%xi0)
      if (allocated(self%f0)) deallocate (self%f0)

   end subroutine anchor_clear

end module moist_cavity_drop_anchor
