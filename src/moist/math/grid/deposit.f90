!> Transfers between scattered points and a periodic cubic lattice of cells
!>
!> - Box: a cube of side `side`, centred on the bounding box of the points, split into
!>   n^3 cells; cell (i, j, k), x fastest, has its centre at `origin + ([i, j, k] - 0.5) h`
!>   with h = side/n
!> - Interpolation P: periodic trilinear interpolation from the eight cell centres around
!>   each point; the eight cells and weights of every point are stored
!> - Restriction R: weighted deposition normalised per cell,
!>   (R f)_a = sum_i w_i P_ia f_i / m_a with the masses m_a = sum_i w_i P_ia;
!>   cells without mass are empty and receive 0
!> - P keeps constants; R keeps constants on cells with mass; with the masses as lattice
!>   weights the two are adjoint: sum_a m_a (R f)_a g_a = sum_i w_i f_i (P g)_i
!> - R P is not the identity, and the masses generally differ from the cell volume
!> - Both are row gathers (R through cell lists in CSR form), so values do not depend on
!>   the thread count
module moist_math_grid_deposit
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_sorter, only: counting_argsort
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none(type, external)
   private

   public :: grid_deposit_type

   !> Interpolation and deposition between points and a periodic lattice of cells
   type :: grid_deposit_type
      !> Cells per axis
      integer :: n = 0
      !> Points
      integer :: npoint = 0
      !> Box side, bohr
      real(wp) :: side = 0.0_wp
      !> Cell spacing, bohr
      real(wp) :: spacing = 0.0_wp
      !> Box corner, bohr (3)
      real(wp) :: origin(3) = 0.0_wp
      !> Cells of the trilinear stencil of each point, x fastest from 1 (8, npoint)
      integer, allocatable :: corner(:, :)
      !> Trilinear weights of each point (8, npoint)
      real(wp), allocatable :: weight(:, :)
      !> Mass m_a of each cell, bohr^3 (n^3)
      real(wp), allocatable :: mass(:)
      !> First deposit of each cell; deposits of cell a are offset(a) .. offset(a+1)-1 (n^3+1)
      integer, allocatable :: offset(:)
      !> Point of each deposit (ndeposit)
      integer, allocatable :: point(:)
      !> Normalised weight w_i P_ia / m_a of each deposit (ndeposit)
      real(wp), allocatable :: share(:)
   contains
      !> Build the stencils of points in a box of given side
      procedure :: build => deposit_build
      !> Lattice values at the points (P)
      procedure :: interpolate => deposit_interpolate
      !> Normalised deposition of point values on the cells (R)
      procedure :: restrict => deposit_restrict
   end type grid_deposit_type

contains

   !> Build the stencils of points in a cube of side `side` around their bounding box
   !>
   !> @param[inout] self      transfers
   !> @param[in]    xyz       point coordinates, bohr (3, npoint)
   !> @param[in]    w         quadrature weights, bohr^3 (npoint)
   !> @param[in]    n         cells per axis
   !> @param[in]    side      box side, bohr
   !> @param[out]   error     cells per axis outside 1..1290, a side that is not positive and
   !>                         finite, coordinates or weights of the wrong size, coordinates
   !>                         that are not finite, weights that are not finite and positive,
   !>                         or points the box does not enclose
   !> @param[in]    nthreads  threads
   subroutine deposit_build(self, xyz, w, n, side, error, nthreads)
      !> Transfers
      class(grid_deposit_type), intent(inout) :: self
      !> Point coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Quadrature weights
      real(wp), intent(in) :: w(:)
      !> Cells per axis
      integer, intent(in) :: n
      !> Box side
      real(wp), intent(in) :: side
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Threads
      integer, intent(in) :: nthreads
      real(wp) :: lo(3), hi(3), u(3), f(3), fw
      integer :: base(3), cell(3, 0:1), i, k, dx, dy, dz, ncell, nkept, a, e
      integer, allocatable :: key(:), entry(:), perm(:), ncount(:)
      character(len=32) :: buf

      self%n = 0
      self%npoint = 0
      if (n < 1 .or. n > 1290) then
         call fatal_error(error, "grid deposit: cells per axis must lie in 1..1290")
         return
      end if
      if (.not. (side > 0.0_wp) .or. .not. ieee_is_finite(side)) then
         call fatal_error(error, "grid deposit: the box side must be positive and finite")
         return
      end if
      if (size(xyz, 1) /= 3 .or. size(xyz, 2) < 1 .or. size(w) /= size(xyz, 2)) then
         call fatal_error(error, "grid deposit: needs coordinates (3, npoint) and one weight per point")
         return
      end if
      if (.not. all(ieee_is_finite(xyz))) then
         call fatal_error(error, "grid deposit: point coordinates must be finite")
         return
      end if
      if (.not. all(ieee_is_finite(w))) then
         call fatal_error(error, "grid deposit: quadrature weights must be finite")
         return
      end if
      if (any(w <= 0.0_wp)) then
         call fatal_error(error, "grid deposit: quadrature weights must be positive; "// &
            & "use angular rules with positive weights")
         return
      end if
      lo = minval(xyz, dim=2)
      hi = maxval(xyz, dim=2)
      if (any(hi - lo > side)) then
         write (buf, "(f0.4)") maxval(hi - lo)
         call fatal_error(error, "grid deposit: the points span "//trim(buf)//" bohr, more than the box side")
         return
      end if

      self%n = n
      self%npoint = size(w)
      self%side = side
      self%spacing = side/real(n, wp)
      self%origin = 0.5_wp*(lo + hi) - 0.5_wp*side
      ncell = n**3
      if (allocated(self%corner)) deallocate (self%corner)
      if (allocated(self%weight)) deallocate (self%weight)
      allocate (self%corner(8, self%npoint), self%weight(8, self%npoint))

      !$omp parallel do schedule(static) num_threads(max(1, nthreads)) private(u, base, f, cell, k, fw)
      do i = 1, self%npoint
         u = (xyz(:, i) - self%origin)/self%spacing - 0.5_wp
         base = floor(u)
         f = u - real(base, wp)
         cell(:, 0) = modulo(base, n)
         cell(:, 1) = modulo(base + 1, n)
         k = 0
         do dz = 0, 1
            do dy = 0, 1
               do dx = 0, 1
                  k = k + 1
                  fw = merge(f(1), 1.0_wp - f(1), dx == 1)*merge(f(2), 1.0_wp - f(2), dy == 1) &
                     & *merge(f(3), 1.0_wp - f(3), dz == 1)
                  self%corner(k, i) = 1 + cell(1, dx) + n*(cell(2, dy) + n*cell(3, dz))
                  self%weight(k, i) = fw
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! Deposits with a nonzero weight, sorted by cell; stable, so in point order per cell
      nkept = count(self%weight > 0.0_wp)
      allocate (key(nkept), entry(nkept), perm(nkept))
      e = 0
      do i = 1, self%npoint
         do k = 1, 8
            if (self%weight(k, i) > 0.0_wp) then
               e = e + 1
               key(e) = self%corner(k, i) - 1
               entry(e) = 8*(i - 1) + k
            end if
         end do
      end do
      call counting_argsort(key, ncell - 1, perm)
      allocate (ncount(ncell), source=0)
      do e = 1, nkept
         ncount(key(e) + 1) = ncount(key(e) + 1) + 1
      end do
      if (allocated(self%offset)) deallocate (self%offset)
      if (allocated(self%point)) deallocate (self%point)
      if (allocated(self%share)) deallocate (self%share)
      if (allocated(self%mass)) deallocate (self%mass)
      allocate (self%offset(ncell + 1), self%point(nkept), self%share(nkept), self%mass(ncell))
      self%offset(1) = 1
      do a = 1, ncell
         self%offset(a + 1) = self%offset(a) + ncount(a)
      end do
      do e = 1, nkept
         i = (entry(perm(e)) - 1)/8 + 1
         k = entry(perm(e)) - 8*(i - 1)
         self%point(e) = i
         self%share(e) = w(i)*self%weight(k, i)
      end do
      !$omp parallel do schedule(static) num_threads(max(1, nthreads)) private(e)
      do a = 1, ncell
         self%mass(a) = 0.0_wp
         do e = self%offset(a), self%offset(a + 1) - 1
            self%mass(a) = self%mass(a) + self%share(e)
         end do
         do e = self%offset(a), self%offset(a + 1) - 1
            self%share(e) = self%share(e)/self%mass(a)
         end do
      end do
      !$omp end parallel do
   end subroutine deposit_build

   !> Lattice values at the points: fm(i, :) = sum_a P_ia fc(a, :)
   !>
   !> @param[in]  self      built transfers
   !> @param[in]  fc        cell values, x fastest (n^3, ncol)
   !> @param[out] fm        point values (npoint, ncol)
   !> @param[in]  nthreads  threads
   subroutine deposit_interpolate(self, fc, fm, nthreads)
      !> Transfers
      class(grid_deposit_type), intent(in) :: self
      !> Cell values
      real(wp), intent(in) :: fc(:, :)
      !> Point values
      real(wp), intent(out) :: fm(:, :)
      !> Threads
      integer, intent(in) :: nthreads
      integer :: i, c, k
      real(wp) :: s

      do c = 1, size(fc, 2)
         !$omp parallel do schedule(static) num_threads(max(1, nthreads)) private(k, s)
         do i = 1, self%npoint
            s = 0.0_wp
            do k = 1, 8
               s = s + self%weight(k, i)*fc(self%corner(k, i), c)
            end do
            fm(i, c) = s
         end do
         !$omp end parallel do
      end do
   end subroutine deposit_interpolate

   !> Normalised deposition: fc(a, :) = sum_i w_i P_ia fm(i, :) / m_a, 0 on empty cells
   !>
   !> @param[in]  self      built transfers
   !> @param[in]  fm        point values (npoint, ncol)
   !> @param[out] fc        cell values, x fastest (n^3, ncol)
   !> @param[in]  nthreads  threads
   subroutine deposit_restrict(self, fm, fc, nthreads)
      !> Transfers
      class(grid_deposit_type), intent(in) :: self
      !> Point values
      real(wp), intent(in) :: fm(:, :)
      !> Cell values
      real(wp), intent(out) :: fc(:, :)
      !> Threads
      integer, intent(in) :: nthreads
      integer :: a, c, e
      real(wp) :: s

      do c = 1, size(fm, 2)
         !$omp parallel do schedule(static) num_threads(max(1, nthreads)) private(e, s)
         do a = 1, size(self%mass)
            s = 0.0_wp
            do e = self%offset(a), self%offset(a + 1) - 1
               s = s + self%share(e)*fm(self%point(e), c)
            end do
            fc(a, c) = s
         end do
         !$omp end parallel do
      end do
   end subroutine deposit_restrict

end module moist_math_grid_deposit
