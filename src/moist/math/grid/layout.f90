!> Layout of a discretization: what a field on it needs to be stored, compared or moved
!>
!> - `cartesian`: `shape = [nx, ny, nz]`, `spacing = dr`, box corner `origin`;
!>   point (i, j, k), x fastest, at `origin + ([i, j, k] - 0.5) dr`
!> - `radial`: `shape = [npts]`, nodes `r` and weights `w`; `spacing = dr` on uniform
!>   grids, 0 otherwise
!> - `points`: coordinates `xyz`, weights `w` and owning atoms `owner` (0: none);
!>   `shape = [n]`, the unpruned point count on atom-centered grids
!>   - optional Fourier basis: modes `kmodes`, spacings `kspacing`, phase reference `kref`
!>   - atom-centered grids: position of each point among the unpruned points
!>     (`point_index`), first shell of each atom (`shell_offset`), shell radii and
!>     angular points per shell, and the angular directions (`angular_grid`)
!> - Grids fill a layout through their `get_layout` binding
!> - `validate` checks the geometry a layout describes and counts its points; users add
!>   their own rules (field rows, point identities) on top
module moist_math_grid_layout
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   implicit none(type, external)
   private

   public :: moist_math_grid_layout_type
   public :: layout_cartesian, layout_radial, layout_points

   !> Layout kinds
   character(len=*), parameter :: layout_cartesian = "cartesian", layout_radial = "radial", &
      & layout_points = "points"

   !> Description of one discretization
   type :: moist_math_grid_layout_type
      !> "cartesian", "radial" or "points"
      character(len=:), allocatable :: kind
      !> Grid dimensions: [nx, ny, nz], [npts], or the (unpruned) point count
      integer, allocatable :: shape(:)
      !> Uniform spacing, bohr; 0 without one
      real(wp) :: spacing = 0.0_wp
      !> Cartesian box corner, bohr (3)
      real(wp), allocatable :: origin(:)
      !> Radial nodes, bohr (npts)
      real(wp), allocatable :: r(:)
      !> Point coordinates, bohr (3, npoint)
      real(wp), allocatable :: xyz(:, :)
      !> Quadrature weights: bohr^3 of point grids, bohr of radial grids (npoint)
      real(wp), allocatable :: w(:)
      !> Owning atom of each point, 0 for none (npoint)
      integer, allocatable :: owner(:)
      !> Fourier basis of a point grid: mode counts (3)
      integer, allocatable :: kmodes(:)
      !> Fourier basis of a point grid: mode spacings, 1/bohr (3)
      real(wp), allocatable :: kspacing(:)
      !> Fourier basis of a point grid: phase reference, bohr (3)
      real(wp), allocatable :: kref(:)
      !> Atom-centered grids: position of each point among the unpruned points, increasing (npoint)
      integer, allocatable :: point_index(:)
      !> Atom-centered grids: first shell of each atom (natom + 1)
      integer, allocatable :: shell_offset(:)
      !> Atom-centered grids: shell radii, bohr (nshell)
      real(wp), allocatable :: shell_radius(:)
      !> Atom-centered grids: angular points of each shell (nshell)
      integer, allocatable :: shell_points(:)
      !> Atom-centered grids: angular directions, "lebedev" for the built-in tables, else "custom"
      character(len=:), allocatable :: angular_grid
   contains
      !> Point coordinates of a Cartesian or point layout, bohr (3, npoint)
      procedure :: points => layout_get_points
      !> Check the described geometry and count the points
      procedure :: validate => layout_validate
   end type moist_math_grid_layout_type

contains

   !> Point coordinates of a Cartesian or point layout
   !>
   !> - Cartesian: rebuilt from the box, x fastest
   !> - Point layouts: their stored coordinates; radial and incomplete layouts: none
   !>
   !> @param[in]  self  layout
   !> @param[out] xyz   coordinates, bohr (3, npoint)
   pure subroutine layout_get_points(self, xyz)
      !> Layout
      class(moist_math_grid_layout_type), intent(in) :: self
      !> Coordinates
      real(wp), allocatable, intent(out) :: xyz(:, :)
      integer :: i, j, k, n

      if (.not. allocated(self%kind)) then
         allocate (xyz(3, 0))
         return
      end if
      select case (self%kind)
      case (layout_cartesian)
         if (.not. allocated(self%shape) .or. .not. allocated(self%origin)) then
            allocate (xyz(3, 0))
            return
         end if
         if (size(self%shape) /= 3 .or. size(self%origin) /= 3) then
            allocate (xyz(3, 0))
            return
         end if
         allocate (xyz(3, product(self%shape)))
         n = 0
         do k = 1, self%shape(3)
            do j = 1, self%shape(2)
               do i = 1, self%shape(1)
                  n = n + 1
                  xyz(:, n) = self%origin + (real([i, j, k], wp) - 0.5_wp)*self%spacing
               end do
            end do
         end do
      case (layout_points)
         if (allocated(self%xyz)) then
            xyz = self%xyz
         else
            allocate (xyz(3, 0))
         end if
      case default
         allocate (xyz(3, 0))
      end select
   end subroutine layout_get_points

   !> Check the geometry a layout describes and count its points
   !>
   !> - Cartesian: three dimensions of at least one whose product fits a default integer,
   !>   a finite positive spacing and a finite corner
   !> - Radial: finite, positive, increasing nodes; `shape`, when given, holds their count
   !> - Points: counted by the coordinates (3, n), else by `point_index`, else by `shape`;
   !>   finite coordinates; owners, when given, one per point and not negative; a Fourier
   !>   basis, when any part is given, complete: three mode counts of at least one, three
   !>   finite positive spacings and a finite phase reference
   !> - Weights, when given, one per point and finite
   !>
   !> @param[in]  self         layout
   !> @param[out] n            points; -1 when the layout is refused
   !> @param[out] reason       first inconsistency; empty for a usable layout
   !> @param[in]  coordinates  point layouts must carry coordinates; default false
   pure subroutine layout_validate(self, n, reason, coordinates)
      !> Layout
      class(moist_math_grid_layout_type), intent(in) :: self
      !> Points
      integer, intent(out) :: n
      !> First inconsistency
      character(len=:), allocatable, intent(out) :: reason
      !> Whether point layouts must carry coordinates
      logical, intent(in), optional :: coordinates
      integer :: count
      logical :: need

      n = -1
      need = .false.
      if (present(coordinates)) need = coordinates
      if (.not. allocated(self%kind)) then
         reason = "the layout is not described"
         return
      end if
      select case (self%kind)
      case (layout_cartesian)
         reason = "the box is not described"
         if (.not. allocated(self%shape) .or. .not. allocated(self%origin)) return
         if (size(self%shape) /= 3 .or. size(self%origin) /= 3) return
         if (.not. (self%spacing > 0.0_wp) .or. .not. ieee_is_finite(self%spacing)) return
         if (.not. all(ieee_is_finite(self%origin))) return
         count = point_count(self%shape)
         if (count < 0) return
      case (layout_radial)
         reason = "the radial nodes are not described"
         if (.not. allocated(self%r)) return
         count = size(self%r)
         if (count < 1) return
         if (.not. all(ieee_is_finite(self%r)) .or. self%r(1) <= 0.0_wp) return
         if (any(self%r(2:) <= self%r(:count - 1))) return
         if (allocated(self%shape)) then
            if (size(self%shape) /= 1) return
            if (self%shape(1) /= count) return
         end if
      case (layout_points)
         reason = "the points are not described"
         if (allocated(self%xyz)) then
            if (size(self%xyz, 1) /= 3) return
            if (.not. all(ieee_is_finite(self%xyz))) return
            count = size(self%xyz, 2)
         else if (need) then
            return
         else if (allocated(self%point_index)) then
            count = size(self%point_index)
         else if (allocated(self%shape)) then
            count = point_count(self%shape)
            if (count < 0) return
         else
            return
         end if
         if (allocated(self%owner)) then
            if (size(self%owner) /= count) return
            if (any(self%owner < 0)) return
         end if
         if (allocated(self%kmodes) .or. allocated(self%kspacing) .or. allocated(self%kref)) then
            reason = "the Fourier basis is not described"
            if (.not. (allocated(self%kmodes) .and. allocated(self%kspacing) .and. allocated(self%kref))) return
            if (size(self%kmodes) /= 3 .or. size(self%kspacing) /= 3 .or. size(self%kref) /= 3) return
            if (point_count(self%kmodes) < 0) return
            if (any(.not. (self%kspacing > 0.0_wp)) .or. .not. all(ieee_is_finite(self%kspacing))) return
            if (.not. all(ieee_is_finite(self%kref))) return
         end if
      case default
         reason = "the layout kind '"//self%kind//"' is unknown"
         return
      end select
      if (allocated(self%w)) then
         reason = "the weights do not match the points"
         if (size(self%w) /= count) return
         if (.not. all(ieee_is_finite(self%w))) return
      end if
      n = count
      reason = ""
   end subroutine layout_validate

   !> Product of dimensions of at least one; -1 for a smaller dimension or a product beyond
   !> a default integer
   !>
   !> @param[in] dims  dimensions
   pure function point_count(dims) result(n)
      !> Dimensions
      integer, intent(in) :: dims(:)
      !> Product
      integer :: n
      integer(int64) :: total
      integer :: i

      n = -1
      total = 1_int64
      do i = 1, size(dims)
         if (dims(i) < 1) return
         total = total*int(dims(i), int64)
         if (total > int(huge(n), int64)) return
      end do
      n = int(total)
   end function point_count

end module moist_math_grid_layout
