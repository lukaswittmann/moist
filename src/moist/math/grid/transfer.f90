!> Transfer of fields between discretizations, from their layouts
!>
!> - Grid transfer: prepared once from a source and a target layout, then applied to
!>   fields on the source; the quantity, its units and its columns stay the same
!>   - each target takes the first method that applies: `copy` from a coinciding source
!>     point, else the interpolation of the source kind: from Cartesian boxes `fourier`
!>     onto Cartesian targets and `trilinear` onto others (`options%cartesian` chooses
!>     one for all), `points` (moving least squares, `moist_math_grid_transfer_points`)
!>     from point layouts, `spline` on radial nodes
!>   - copies: exactly one source point within `copy_tolerance` of the target; for any
!>     layout, any order or subset of the source points; targets at two coinciding
!>     source points are not copied
!>   - targets beyond the reach of the source stay uncovered; the caller's values there
!>     are kept, and without a coverage mask an incomplete transfer is an error
!>
!>    call transfer%prepare(old_layout, new_layout, error, map=fit)
!>    call transfer%apply(old_field, new_field, error, covered=covered)
!>
!> - Rigid fit: the rotation and translation that map the target solute onto the
!>   stored one, so that target points can be evaluated in the frame of the stored field
!>   - unweighted centroids; one atom, or coordinates that differ by a translation only,
!>     give the translation without a rotation fit
!>   - one atom: translation only; both sides collinear: the smallest rotation that
!>     aligns the axes, no rotation about them; otherwise Horn's quaternion fit
!> - Radial transfer: cubic spline through the stored nodes and the mirror images of
!>   the first eight at -r, so the interpolant is close to even at r = 0; nodes are
!>   evaluated at |r|, and those beyond the last stored node are not covered
!> - Coverage of a Cartesian box: its volume
!> - Spectral transfer from a layout with a Fourier basis: Cartesian boxes use their FFT
!>   modes, point grids their reciprocal grid
!>   - weighted type-1 NUFFT of the stored values onto the modes, type-2 NUFFT at the
!>     points; exact trigonometric interpolation on a Cartesian lattice, an
!>     approximation on point grids, whose quadrature does not invert the modes
module moist_math_grid_transfer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_layout, only: moist_math_grid_layout_type, layout_cartesian, layout_points, layout_radial
   use moist_math_kdtree, only: moist_kdtree_type
   use moist_math_grid_transfer_points, only: point_stencil, point_reach, k_last, stencil_fitted, stencil_far, &
      & stencil_conditioning, stencil_amplification
   use moist_math_lapack_syev, only: lapack_syev
   use moist_math_lapack_kinds, only: lapack_ik
   use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_double, c_double_complex
   use finufft_mod, only: finufft_opts
   use moist_math_fft_finufft, only: finufft_default_opts, finufft_makeplan, finufft_setpts, &
      & finufft_execute, finufft_destroy
   !$ use omp_lib, only: omp_get_max_threads
   implicit none(type, external)
   private

   public :: rigid_fit_type, fit_rigid, transfer_radial, layout_covers, transfer_fourier
   public :: has_fourier_basis
   public :: grid_transfer_type, grid_transfer_options_type, grid_transfer_report_type
   public :: copy_tolerance

   !> Largest distance of a target from the source point it copies, bohr
   real(wp), parameter :: copy_tolerance = 1.0e-9_wp

   !> Targets fitted together, so the temporary stencils of a block stay small
   integer, parameter :: stencil_block = 16384

   !> Options of a grid transfer
   type :: grid_transfer_options_type
      !> Interpolation from Cartesian sources, "fourier" or "trilinear"; unset: Fourier onto
      !> Cartesian targets, trilinear onto others
      character(len=:), allocatable :: cartesian
      !> Exact copies only; the other targets stay uncovered
      logical :: copy_only = .false.
      !> Threads of the interpolation and of the copy and reach searches; 0 for the default
      integer :: nthreads = 0
      !> Largest memory of stored point stencils, bytes; above it they are fitted again
      !> when applied. While the stencil buffers grow, the old and the new buffer coexist
      !> for a moment, so up to twice this is allocated briefly
      integer(int64) :: stencil_memory = 500000000_int64
      !> Decide only whether the transfer covers a target: from point layouts the fit stops
      !> at the first accepted stencil, and without copies the report counts that one target;
      !> a point-layout transfer that leaves targets to the fit then cannot be applied, other
      !> sources ignore the option
      logical :: check_only = .false.
   end type grid_transfer_options_type

   !> What a grid transfer does with the targets
   type :: grid_transfer_report_type
      !> Targets that take the value of a coinciding source point
      integer :: ncopied = 0
      !> Targets interpolated
      integer :: ninterpolated = 0
      !> Targets without a value
      integer :: nuncovered = 0
      !> Targets not copied under `copy_only`
      integer :: nexcluded = 0
      !> Targets not copied because two source points coincide there
      integer :: nduplicate = 0
      !> Methods used: "copy", "fourier", "trilinear", "points", "spline" or
      !> "copy+<interpolation>"; empty without a covered target
      character(len=:), allocatable :: method
      !> Why targets stay uncovered; empty when every target is covered
      character(len=:), allocatable :: reason
   end type grid_transfer_report_type

   !> Transfer between two layouts, prepared once and applied to fields on the source
   type :: grid_transfer_type
      private
      !> Source points; -1 before `prepare`
      integer :: nold = -1
      !> Source point copied to each target, 0 for none (ntarget)
      integer, allocatable :: copy_from(:)
      !> Interpolated targets (ntarget)
      logical, allocatable :: interpolated(:)
      !> Targets not copied under `copy_only` (ntarget)
      logical, allocatable :: excluded(:)
      !> Targets at coinciding source points
      integer :: nduplicate = 0
      !> Interpolation: "fourier", "trilinear", "points" or "spline"; empty without one
      character(len=:), allocatable :: method
      !> Why targets stay uncovered
      character(len=:), allocatable :: reason
      !> Source geometry of the interpolation: the box, or the radial nodes
      type(moist_math_grid_layout_type) :: source
      !> Interpolated 3D targets in the source frame, bohr (3, ninterpolated)
      real(wp), allocatable :: y(:, :)
      !> Radial targets, bohr (ntarget)
      real(wp), allocatable :: r(:)
      !> Source points of a point layout; kept only for stencils fitted again when applied
      type(moist_kdtree_type) :: tree
      !> Reach of each source point of a point layout, bohr; kept with the tree
      real(wp), allocatable :: reach(:)
      !> Whether the point stencils are stored; otherwise they are fitted again when applied
      logical :: stored = .false.
      !> Largest memory of stored point stencils, bytes
      integer(int64) :: stencil_memory = 0_int64
      !> Whether the stencils of point layouts were fitted; not with `check_only`
      logical :: fitted = .true.
      !> Stored point stencils: first entry of each interpolated target (ninterpolated + 1),
      !> source points of the entries
      integer, allocatable :: stencil_offset(:), stencil_point(:)
      !> Weights of the stored stencil entries
      real(wp), allocatable :: stencil_weight(:)
      !> Threads of the interpolation
      integer :: nthreads = 0
   contains
      !> Decide the method of each target from the layouts
      procedure :: prepare => transfer_prepare
      !> Move a field from the source onto the targets
      procedure :: apply => transfer_apply
      !> Coverage and methods
      procedure :: get_report => transfer_get_report
      !> Release the prepared data
      procedure :: destroy => transfer_destroy
   end type grid_transfer_type

   !> Map from the target frame into the stored frame: y' = rot (y - c_target) + c_stored
   type :: rigid_fit_type
      !> Rotation
      real(wp) :: rot(3, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      !> Centroid of the target atoms, bohr
      real(wp) :: c_target(3) = 0.0_wp
      !> Centroid of the stored atoms, bohr
      real(wp) :: c_stored(3) = 0.0_wp
   contains
      !> Target points in the stored frame
      procedure :: map => rigid_fit_map
      !> Whether the rotation is exactly the identity
      procedure :: is_identity => rigid_fit_is_identity
   end type rigid_fit_type

contains

   !> Rigid map of the target solute onto the stored one
   !>
   !> @param[in]  stored  stored atom coordinates, bohr (3, natom)
   !> @param[in]  target  target atom coordinates in the same order, bohr (3, natom)
   !> @param[out] fit     map from the target frame into the stored frame
   !> @param[out] error   coordinates not (3, natom), atom counts differ or the eigensolver failed
   subroutine fit_rigid(stored, target, fit, error)
      !> Stored atom coordinates
      real(wp), intent(in) :: stored(:, :)
      !> Target atom coordinates
      real(wp), intent(in) :: target(:, :)
      !> Map
      type(rigid_fit_type), intent(out) :: fit
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: s(:, :), t(:, :)
      real(wp) :: cov_s(3, 3), cov_t(3, 3), axis_s(3), axis_t(3), extent
      integer :: n

      n = size(stored, 2)
      if (size(target, 2) /= n .or. size(stored, 1) /= 3 .or. size(target, 1) /= 3) then
         call fatal_error(error, "rigid fit needs the same atoms on both sides")
         return
      end if
      if (n == 0) return
      fit%c_stored = sum(stored, dim=2)/real(n, wp)
      fit%c_target = sum(target, dim=2)/real(n, wp)
      s = stored - spread(fit%c_stored, 2, n)
      t = target - spread(fit%c_target, 2, n)
      extent = max(1.0_wp, maxval(abs(s)))
      ! Pure translation, and single atoms, keep the identity exactly
      if (n == 1 .or. all(abs(s - t) <= 1.0e-12_wp*extent)) return
      cov_s = matmul(s, transpose(s))
      cov_t = matmul(t, transpose(t))
      ! A bent and a straight structure take the full fit, which minimizes the
      ! correspondence error; only two collinear structures align their axes
      if (rank_one(cov_s) .and. rank_one(cov_t)) then
         call principal_axis(s, axis_s)
         call principal_axis(t, axis_t)
         if (sum(matmul(axis_s, s)*matmul(axis_t, t)) < 0.0_wp) axis_t = -axis_t
         fit%rot = align_axes(axis_t, axis_s)
      else
         call horn_rotation(s, t, fit%rot, error)
         if (allocated(error)) return
      end if
   end subroutine fit_rigid

   !> Target points in the stored frame
   !>
   !> @param[in]  self  map
   !> @param[in]  y     target points, bohr (3, n)
   !> @param[out] ys    the same points in the stored frame, bohr (3, n)
   pure subroutine rigid_fit_map(self, y, ys)
      !> Map
      class(rigid_fit_type), intent(in) :: self
      !> Target points
      real(wp), intent(in) :: y(:, :)
      !> Points in the stored frame
      real(wp), allocatable, intent(out) :: ys(:, :)

      if (self%is_identity()) then
         ys = y - spread(self%c_target - self%c_stored, 2, size(y, 2))
      else
         ys = matmul(self%rot, y - spread(self%c_target, 2, size(y, 2))) + spread(self%c_stored, 2, size(y, 2))
      end if
   end subroutine rigid_fit_map

   !> Whether the rotation is exactly the identity; such maps only translate, without rounding
   !>
   !> @param[in] self  map
   pure logical function rigid_fit_is_identity(self)
      !> Map
      class(rigid_fit_type), intent(in) :: self
      rigid_fit_is_identity = all(self%rot == identity3())
   end function rigid_fit_is_identity

   !> Whether a positive semidefinite 3x3 matrix has rank one or less
   !>
   !> @param[in] cov  covariance of centered coordinates
   logical function rank_one(cov)
      !> Covariance
      real(wp), intent(in) :: cov(3, 3)
      real(wp) :: a(3, 3), w(3), work(16)
      integer(lapack_ik) :: info

      a = cov
      call lapack_syev("N", "U", 3_lapack_ik, a, 3_lapack_ik, w, work, 16_lapack_ik, info)
      rank_one = info /= 0
      if (info == 0) rank_one = w(2) <= 1.0e-10_wp*max(w(3), tiny(1.0_wp))
   end function rank_one

   !> Unit direction of largest spread of centered coordinates
   !>
   !> @param[in]  x     centered coordinates (3, n)
   !> @param[out] axis  direction
   subroutine principal_axis(x, axis)
      !> Centered coordinates
      real(wp), intent(in) :: x(:, :)
      !> Direction
      real(wp), intent(out) :: axis(3)
      real(wp) :: a(3, 3), w(3), work(16)
      integer(lapack_ik) :: info

      a = matmul(x, transpose(x))
      call lapack_syev("V", "U", 3_lapack_ik, a, 3_lapack_ik, w, work, 16_lapack_ik, info)
      axis = a(:, 3)
      if (info /= 0) axis = [1.0_wp, 0.0_wp, 0.0_wp]
   end subroutine principal_axis

   !> Smallest rotation that turns the unit vector `from` into the unit vector `to`
   !>
   !> - The turn about the common normal by the angle between them, which leaves the
   !>   normal in place
   !> - Within 1e-8 of anti-parallel the normal is not determined: a half turn about a
   !>   normal of `from`, then the small turn onto `to`
   !>
   !> @param[in] from  unit vector
   !> @param[in] to    unit vector
   pure function align_axes(from, to) result(rot)
      !> Unit vector to turn
      real(wp), intent(in) :: from(3)
      !> Its image
      real(wp), intent(in) :: to(3)
      !> Rotation
      real(wp) :: rot(3, 3)
      real(wp) :: v(3), sine, cosine, k(3, 3), p(3)

      v = cross(from, to)
      sine = norm2(v)
      cosine = dot_product(from, to)
      if (sine > 1.0e-8_wp) then
         k = skew(v/sine)
         rot = identity3() + sine*k + (1.0_wp - cosine)*matmul(k, k)
      else if (cosine > 0.0_wp) then
         rot = small_turn(from, to)
      else
         p = perpendicular(from)
         rot = matmul(small_turn(-from, to), 2.0_wp*spread(p, 2, 3)*spread(p, 1, 3) - identity3())
      end if
   contains
      !> Rotation of `a` into a nearly parallel `b`
      !>
      !> @param[in] a  unit vector to turn
      !> @param[in] b  unit vector it is turned into, not antiparallel to `a`
      pure function small_turn(a, b) result(r)
         !> Unit vectors
         real(wp), intent(in) :: a(3), b(3)
         real(wp) :: r(3, 3), kk(3, 3)
         kk = skew(cross(a, b))
         r = identity3() + kk + matmul(kk, kk)/(1.0_wp + dot_product(a, b))
      end function small_turn

      !> Cross product
      !>
      !> @param[in] a  first factor
      !> @param[in] b  second factor
      pure function cross(a, b) result(c)
         !> Factors
         real(wp), intent(in) :: a(3), b(3)
         real(wp) :: c(3)
         c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
      end function cross

      !> Matrix of the cross product with u
      !>
      !> @param[in] u  vector
      pure function skew(u) result(m)
         !> Vector
         real(wp), intent(in) :: u(3)
         real(wp) :: m(3, 3)
         m = reshape([0.0_wp, u(3), -u(2), -u(3), 0.0_wp, u(1), u(2), -u(1), 0.0_wp], [3, 3])
      end function skew
   end function align_axes

   !> A unit vector perpendicular to a unit vector
   !>
   !> @param[in] u  unit vector
   pure function perpendicular(u) result(p)
      !> Unit vector
      real(wp), intent(in) :: u(3)
      !> Perpendicular unit vector
      real(wp) :: p(3)

      if (abs(u(1)) < 0.9_wp) then
         p = [0.0_wp, u(3), -u(2)]
      else
         p = [-u(3), 0.0_wp, u(1)]
      end if
      p = p/norm2(p)
   end function perpendicular

   !> Rotation R minimizing sum |R t_i - s_i|^2 (Horn, quaternion form)
   !>
   !> @param[in]  s      centered stored coordinates (3, n)
   !> @param[in]  t      centered target coordinates (3, n)
   !> @param[out] rot    rotation
   !> @param[out] error  eigensolver failure
   subroutine horn_rotation(s, t, rot, error)
      !> Centered stored coordinates
      real(wp), intent(in) :: s(:, :)
      !> Centered target coordinates
      real(wp), intent(in) :: t(:, :)
      !> Rotation
      real(wp), intent(out) :: rot(3, 3)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: m(3, 3), k(4, 4), w(4), work(32), q(4)
      integer(lapack_ik) :: info

      m = matmul(t, transpose(s))
      k(1, 1) = m(1, 1) + m(2, 2) + m(3, 3)
      k(2, 1) = m(2, 3) - m(3, 2)
      k(3, 1) = m(3, 1) - m(1, 3)
      k(4, 1) = m(1, 2) - m(2, 1)
      k(2, 2) = m(1, 1) - m(2, 2) - m(3, 3)
      k(3, 2) = m(1, 2) + m(2, 1)
      k(4, 2) = m(3, 1) + m(1, 3)
      k(3, 3) = -m(1, 1) + m(2, 2) - m(3, 3)
      k(4, 3) = m(2, 3) + m(3, 2)
      k(4, 4) = -m(1, 1) - m(2, 2) + m(3, 3)
      call lapack_syev("V", "L", 4_lapack_ik, k, 4_lapack_ik, w, work, 32_lapack_ik, info)
      if (info /= 0) then
         call fatal_error(error, "rigid fit: eigensolver failed")
         return
      end if
      q = k(:, 4)
      rot(1, :) = [q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2, 2.0_wp*(q(2)*q(3) - q(1)*q(4)), &
         & 2.0_wp*(q(2)*q(4) + q(1)*q(3))]
      rot(2, :) = [2.0_wp*(q(2)*q(3) + q(1)*q(4)), q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2, &
         & 2.0_wp*(q(3)*q(4) - q(1)*q(2))]
      rot(3, :) = [2.0_wp*(q(2)*q(4) - q(1)*q(3)), 2.0_wp*(q(3)*q(4) + q(1)*q(2)), &
         & q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2]
   end subroutine horn_rotation

   !> 3x3 identity
   pure function identity3() result(eye)
      !> Identity
      real(wp) :: eye(3, 3)
      eye = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
   end function identity3

   !> Cubic spline of stored radial columns at other nodes
   !>
   !> - Natural spline through the stored nodes and the mirror images of the first eight
   !>   at -r, so the slope near r = 0 follows the even extension
   !> - Nodes are evaluated at |r|, the even extension; nodes with |r| beyond the last
   !>   stored node are not covered and get zero
   !>
   !> @param[in]  r_stored  stored nodes, increasing and positive, bohr (n)
   !> @param[in]  f_stored  stored columns (n, ncol)
   !> @param[in]  r         nodes to evaluate at, bohr (m)
   !> @param[out] f         interpolated columns (m, ncol)
   !> @param[out] covered   nodes within the stored range (m)
   !> @param[out] error     too few stored nodes or columns of another length, or stored nodes
   !>                       not finite, increasing and positive
   subroutine transfer_radial(r_stored, f_stored, r, f, covered, error)
      !> Stored nodes
      real(wp), intent(in) :: r_stored(:)
      !> Stored columns
      real(wp), intent(in) :: f_stored(:, :)
      !> Nodes to evaluate at
      real(wp), intent(in) :: r(:)
      !> Interpolated columns
      real(wp), allocatable, intent(out) :: f(:, :)
      !> Covered nodes
      logical, allocatable, intent(out) :: covered(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: x(:), y(:), d2(:)
      integer :: n, nm, col, i, j

      n = size(r_stored)
      if (n < 2 .or. size(f_stored, 1) /= n) then
         call fatal_error(error, "radial transfer needs at least two stored nodes")
         return
      end if
      if (.not. all(ieee_is_finite(r_stored)) .or. r_stored(1) <= 0.0_wp .or. any(r_stored(2:) <= r_stored(:n - 1))) then
         call fatal_error(error, "radial transfer needs finite, increasing, positive stored nodes")
         return
      end if
      ! Mirror up to eight nodes through r = 0, so the free end lies away from the origin
      nm = min(8, n)
      allocate (x(n + nm))
      x(1:nm) = -r_stored(nm:1:-1)
      x(nm + 1:) = r_stored
      allocate (f(size(r), size(f_stored, 2)), source=0.0_wp)
      covered = abs(r) <= r_stored(n)*(1.0_wp + 1.0e-12_wp)
      allocate (y(n + nm))
      do col = 1, size(f_stored, 2)
         y(1:nm) = f_stored(nm:1:-1, col)
         y(nm + 1:) = f_stored(:, col)
         call natural_spline(x, y, d2)
         j = 1
         do i = 1, size(r)
            if (.not. covered(i)) cycle
            call spline_value(x, y, d2, abs(r(i)), j, f(i, col))
         end do
      end do
   end subroutine transfer_radial

   !> Second derivatives of the natural cubic spline through (x, y)
   !>
   !> @param[in]  x   increasing nodes
   !> @param[in]  y   values
   !> @param[out] d2  second derivatives at the nodes
   pure subroutine natural_spline(x, y, d2)
      !> Nodes
      real(wp), intent(in) :: x(:)
      !> Values
      real(wp), intent(in) :: y(:)
      !> Second derivatives
      real(wp), allocatable, intent(out) :: d2(:)
      real(wp), allocatable :: u(:)
      real(wp) :: sig, p
      integer :: n, i

      n = size(x)
      allocate (d2(n), u(n), source=0.0_wp)
      do i = 2, n - 1
         sig = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1))
         p = sig*d2(i - 1) + 2.0_wp
         d2(i) = (sig - 1.0_wp)/p
         u(i) = (y(i + 1) - y(i))/(x(i + 1) - x(i)) - (y(i) - y(i - 1))/(x(i) - x(i - 1))
         u(i) = (6.0_wp*u(i)/(x(i + 1) - x(i - 1)) - sig*u(i - 1))/p
      end do
      d2(n) = 0.0_wp
      do i = n - 1, 1, -1
         d2(i) = d2(i)*d2(i + 1) + u(i)
      end do
   end subroutine natural_spline

   !> Value of the cubic spline at xv; `j` is a search hint, advanced for increasing xv
   !>
   !> @param[in]    x   increasing nodes
   !> @param[in]    y   values
   !> @param[in]    d2  second derivatives
   !> @param[in]    xv  abscissa within the nodes
   !> @param[inout] j   interval hint
   !> @param[out]   v   value
   pure subroutine spline_value(x, y, d2, xv, j, v)
      !> Nodes
      real(wp), intent(in) :: x(:)
      !> Values
      real(wp), intent(in) :: y(:)
      !> Second derivatives
      real(wp), intent(in) :: d2(:)
      !> Abscissa
      real(wp), intent(in) :: xv
      !> Interval hint
      integer, intent(inout) :: j
      !> Value
      real(wp), intent(out) :: v
      real(wp) :: h, a, b
      integer :: lo, hi, mid

      if (j < 1 .or. j >= size(x)) j = 1
      if (.not. (x(j) <= xv .and. xv <= x(j + 1))) then
         lo = 1
         hi = size(x)
         do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (x(mid) > xv) then
               hi = mid
            else
               lo = mid
            end if
         end do
         j = lo
      end if
      h = x(j + 1) - x(j)
      a = (x(j + 1) - xv)/h
      b = (xv - x(j))/h
      v = a*y(j) + b*y(j + 1) + ((a**3 - a)*d2(j) + (b**3 - b)*d2(j + 1))*h*h/6.0_wp
   end subroutine spline_value

   !> Points of the stored frame inside a stored Cartesian box
   !>
   !> - Other layouts, layouts that fail `validate` and points without three coordinates
   !>   cover nothing
   !>
   !> @param[in]  layout   stored layout
   !> @param[in]  y        points in the stored frame, bohr (3, n)
   !> @param[out] covered  covered points (n)
   subroutine layout_covers(layout, y, covered)
      !> Stored layout
      type(moist_math_grid_layout_type), intent(in) :: layout
      !> Points in the stored frame
      real(wp), intent(in) :: y(:, :)
      !> Covered points
      logical, allocatable, intent(out) :: covered(:)
      character(len=:), allocatable :: reason
      real(wp) :: lo(3), hi(3), tol
      integer :: i, n

      allocate (covered(size(y, 2)), source=.false.)
      call layout%validate(n, reason)
      if (len(reason) > 0 .or. size(y, 1) /= 3) return
      select case (layout%kind)
      case (layout_cartesian)
         tol = 1.0e-10_wp*layout%spacing
         lo = layout%origin - tol
         hi = layout%origin + real(layout%shape, wp)*layout%spacing + tol
         do i = 1, size(y, 2)
            covered(i) = all(y(:, i) >= lo) .and. all(y(:, i) <= hi)
         end do
      case default
         return
      end select
   end subroutine layout_covers

   !> Whether a layout carries a Fourier basis for spectral transfer
   !>
   !> @param[in] layout  stored layout
   pure logical function has_fourier_basis(layout)
      !> Stored layout
      type(moist_math_grid_layout_type), intent(in) :: layout
      character(len=:), allocatable :: reason
      integer :: n

      has_fourier_basis = .false.
      call layout%validate(n, reason, coordinates=.true.)
      if (len(reason) > 0) return
      select case (layout%kind)
      case (layout_cartesian)
         has_fourier_basis = .true.
      case (layout_points)
         ! A validated basis is complete once any part is there
         has_fourier_basis = allocated(layout%kmodes) .and. allocated(layout%w)
      case default
         has_fourier_basis = .false.
      end select
   end function has_fourier_basis

   !> Spectral transfer of stored columns to points in the stored frame
   !>
   !> - Weighted type-1 NUFFT onto the modes of the stored layout, type-2 NUFFT at `y`,
   !>   scaled by 1/V of the periodic box; coordinates are wrapped into one period
   !>
   !> @param[in]  layout    stored layout with a Fourier basis
   !> @param[in]  f         stored columns (npoint, ncol)
   !> @param[in]  y         points in the stored frame, bohr (3, m)
   !> @param[out] g         values at the points (m, ncol)
   !> @param[out] error     missing basis, shape mismatch or NUFFT failure
   !> @param[in]  nthreads  NUFFT threads; 0 for the default
   !> @param[in]  tol       NUFFT tolerance; default 1e-9
   subroutine transfer_fourier(layout, f, y, g, error, nthreads, tol)
      !> Stored layout
      type(moist_math_grid_layout_type), intent(in) :: layout
      !> Stored columns
      real(wp), intent(in) :: f(:, :)
      !> Points in the stored frame
      real(wp), intent(in) :: y(:, :)
      !> Values at the points
      real(wp), allocatable, intent(out) :: g(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> NUFFT threads
      integer, intent(in), optional :: nthreads
      !> NUFFT tolerance
      real(wp), intent(in), optional :: tol
      real(wp), allocatable :: xs(:, :), w(:)
      real(c_double), allocatable :: sx(:), sy(:), sz(:), tx(:), ty(:), tz(:), dummy(:)
      complex(c_double_complex), allocatable :: cj(:), modes(:), values(:)
      type(finufft_opts) :: opts
      integer(c_int64_t) :: plan1, plan2, nm(3)
      integer(c_int) :: ier
      real(wp) :: dk(3), kref(3), scale, eps
      integer :: col

      if (.not. has_fourier_basis(layout)) then
         call fatal_error(error, "spectral transfer needs a layout with a Fourier basis")
         return
      end if
      if (size(y, 1) /= 3) then
         call fatal_error(error, "spectral transfer: target points need three coordinates")
         return
      end if
      select case (layout%kind)
      case (layout_cartesian)
         nm = int(layout%shape, c_int64_t)
         dk = 8.0_wp*atan(1.0_wp)/(real(layout%shape, wp)*layout%spacing)
         kref = layout%origin + 0.5_wp*layout%spacing
         call layout%points(xs)
         allocate (w(size(xs, 2)), source=layout%spacing**3)
      case default
         nm = int(layout%kmodes, c_int64_t)
         dk = layout%kspacing
         kref = layout%kref
         xs = layout%xyz
         allocate (w, source=layout%w)
      end select
      if (size(f, 1) /= size(xs, 2)) then
         call fatal_error(error, "spectral transfer: stored columns do not match the stored points")
         return
      end if
      allocate (g(size(y, 2), size(f, 2)), source=0.0_wp)
      if (size(y, 2) == 0 .or. size(f, 2) == 0) return
      eps = 1.0e-9_wp
      if (present(tol)) eps = tol
      scale = product(dk)/(8.0_wp*atan(1.0_wp))**3
      ! Coordinate arrays stay alive while the plans point at them
      call scaled_coordinates(xs, kref, dk, sx, sy, sz)
      call scaled_coordinates(y, kref, dk, tx, ty, tz)
      allocate (dummy(1), modes(product(nm)), cj(size(xs, 2)), values(size(y, 2)))
      call finufft_default_opts(opts)
      opts%modeord = 0
      if (present(nthreads)) opts%nthreads = int(max(0, nthreads), c_int)
      plan1 = 0_c_int64_t
      plan2 = 0_c_int64_t
      call finufft_makeplan(1_c_int, 3_c_int, nm, -1_c_int, 1_c_int, real(eps, c_double), plan1, opts, ier)
      if (ier == 0) call finufft_setpts(plan1, int(size(sx), c_int64_t), sx, sy, sz, 0_c_int64_t, dummy, dummy, dummy, ier)
      if (ier /= 0) then
         call destroy_plans()
         call fatal_error(error, "spectral transfer: FINUFFT setup failed")
         return
      end if
      block
         complex(c_double_complex), allocatable :: coeff(:, :)
         allocate (coeff(product(nm), size(f, 2)))
         do col = 1, size(f, 2)
            cj = cmplx(w*f(:, col), 0.0_wp, c_double_complex)
            call finufft_execute(plan1, cj, modes, ier)
            if (ier /= 0) exit
            coeff(:, col) = modes
         end do
         call destroy_plans()
         if (ier /= 0) then
            call fatal_error(error, "spectral transfer: FINUFFT execution failed")
            return
         end if
         call finufft_makeplan(2_c_int, 3_c_int, nm, 1_c_int, 1_c_int, real(eps, c_double), plan2, opts, ier)
         if (ier == 0) call finufft_setpts(plan2, int(size(tx), c_int64_t), tx, ty, tz, 0_c_int64_t, dummy, dummy, dummy, ier)
         do col = 1, size(f, 2)
            if (ier /= 0) exit
            modes = coeff(:, col)
            call finufft_execute(plan2, values, modes, ier)
            g(:, col) = scale*real(values, wp)
         end do
         call destroy_plans()
         if (ier /= 0) call fatal_error(error, "spectral transfer: FINUFFT evaluation failed")
      end block
   contains
      !> Free the plans that exist
      subroutine destroy_plans()
         integer(c_int) :: ierr
         if (plan1 /= 0_c_int64_t) call finufft_destroy(plan1, ierr)
         if (plan2 /= 0_c_int64_t) call finufft_destroy(plan2, ierr)
         plan1 = 0_c_int64_t
         plan2 = 0_c_int64_t
      end subroutine destroy_plans
   end subroutine transfer_fourier

   !> Coordinates relative to the phase reference, scaled by the mode spacing and
   !> wrapped into [-pi, pi)
   !>
   !> @param[in]  x     points, bohr (3, n)
   !> @param[in]  kref  phase reference, bohr
   !> @param[in]  dk    mode spacings, 1/bohr
   !> @param[out] sx    scaled x
   !> @param[out] sy    scaled y
   !> @param[out] sz    scaled z
   pure subroutine scaled_coordinates(x, kref, dk, sx, sy, sz)
      !> Points
      real(wp), intent(in) :: x(:, :)
      !> Phase reference
      real(wp), intent(in) :: kref(3)
      !> Mode spacings
      real(wp), intent(in) :: dk(3)
      !> Scaled coordinates
      real(c_double), allocatable, intent(out) :: sx(:), sy(:), sz(:)
      real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)

      sx = modulo((x(1, :) - kref(1))*dk(1) + pi, 2.0_wp*pi) - pi
      sy = modulo((x(2, :) - kref(2))*dk(2) + pi, 2.0_wp*pi) - pi
      sz = modulo((x(3, :) - kref(3))*dk(3) + pi, 2.0_wp*pi) - pi
   end subroutine scaled_coordinates

   !> Decide the method of each target from the layouts
   !>
   !> - Each target takes the first method that applies (see the module description)
   !> - `map` takes the targets into the frame of the source; without it both share one
   !> - The transfer keeps the source geometry it needs: `old` may change or go away
   !>   afterwards; prepare again after the grids change
   !>
   !> @param[out] self     prepared transfer
   !> @param[in]  old      source layout
   !> @param[in]  new      target layout
   !> @param[out] error    a layout not described, radial and 3D layouts mixed, an unknown method
   !> @param[in]  options  methods and threads
   !> @param[in]  map      rigid map of the targets into the source frame
   subroutine transfer_prepare(self, old, new, error, options, map)
      !> Prepared transfer
      class(grid_transfer_type), intent(out) :: self
      !> Source layout
      type(moist_math_grid_layout_type), intent(in) :: old
      !> Target layout
      type(moist_math_grid_layout_type), intent(in) :: new
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Methods and threads
      type(grid_transfer_options_type), intent(in), optional :: options
      !> Map of the targets into the source frame
      type(rigid_fit_type), intent(in), optional :: map
      type(grid_transfer_options_type) :: opts
      real(wp), allocatable :: ys(:, :)
      logical, allocatable :: inside(:), todo(:)
      integer :: nold, nnew, i

      if (present(options)) opts = options
      if (allocated(opts%cartesian)) then
         if (opts%cartesian /= "fourier" .and. opts%cartesian /= "trilinear") then
            call fatal_error(error, "grid transfer: unknown method for Cartesian sources '"//opts%cartesian//"'")
            return
         end if
      end if
      call described(old, "source", nold, error)
      if (.not. allocated(error)) call described(new, "target", nnew, error)
      if (allocated(error)) return
      if ((old%kind == layout_radial) .neqv. (new%kind == layout_radial)) then
         call fatal_error(error, "grid transfer: radial and 3D layouts cannot be mixed")
         return
      end if
      self%nold = nold
      self%nthreads = max(0, opts%nthreads)
      self%stencil_memory = opts%stencil_memory
      allocate (self%copy_from(nnew), source=0)
      allocate (self%interpolated(nnew), self%excluded(nnew), todo(nnew), source=.false.)
      self%method = ""
      self%reason = "outside the source grid"
      select case (old%kind)
      case (layout_radial)
         call copy_radial(old%r, new%r, self%copy_from, self%nduplicate)
         todo = self%copy_from == 0
         if (.not. any(todo)) then
            continue
         else if (opts%copy_only) then
            self%excluded = todo
         else if (nold < 2) then
            self%reason = "no method for these grids"
         else
            self%method = "spline"
            self%source%kind = layout_radial
            self%source%r = old%r
            self%r = new%r
            self%interpolated = todo .and. new%r <= old%r(nold)*(1.0_wp + 1.0e-12_wp)
         end if
      case (layout_cartesian)
         call target_points(new, ys, map)
         call copy_lattice(old, ys, self%copy_from, self%nduplicate)
         todo = self%copy_from == 0
         if (.not. any(todo)) then
            continue
         else if (opts%copy_only) then
            self%excluded = todo
         else
            ! Fourier onto lattices, trilinear onto other points, unless chosen
            self%method = "trilinear"
            if (new%kind == layout_cartesian) self%method = "fourier"
            if (allocated(opts%cartesian)) self%method = opts%cartesian
            self%source%kind = layout_cartesian
            self%source%shape = old%shape
            self%source%spacing = old%spacing
            self%source%origin = old%origin
            call layout_covers(old, ys, inside)
            self%interpolated = todo .and. inside
            self%y = ys(:, pack([(i, i=1, nnew)], self%interpolated))
         end if
      case default
         call target_points(new, ys, map)
         call self%tree%build(old%xyz)
         call point_reach(self%tree, copy_tolerance, self%reach, threads(self%nthreads))
         call self%tree%set_reach(self%reach)
         call copy_points(self%tree, ys, threads(self%nthreads), self%copy_from, self%nduplicate, todo)
         if (opts%copy_only) then
            self%excluded = self%copy_from == 0
         else if (any(todo) .and. opts%check_only) then
            self%method = "points"
            self%fitted = .false.
            if (.not. any(self%copy_from > 0)) call fit_first_point(self, ys, todo)
         else if (any(todo)) then
            self%method = "points"
            call fit_points(self, ys, todo)
         end if
         ! Copies, and transfers prepared as a check, never search the source again
         if (self%method /= "points" .or. .not. self%fitted) then
            call self%tree%destroy()
            if (allocated(self%reach)) deallocate (self%reach)
         end if
      end select
      if (any(self%excluded)) then
         self%reason = "only exact copies requested"
      else if (.not. any(self%copy_from > 0 .or. self%interpolated) .and. self%reason == "outside the source grid") then
         self%reason = "covers no point of this grid"
      end if
   end subroutine transfer_prepare

   !> Move a field from the source onto the targets
   !>
   !> - Covered rows of `field_new` are overwritten, the others keep the caller's values
   !> - Without `covered`, a target without a value is an error and `field_new` is unchanged
   !> - The columns move together and keep their order
   !>
   !> @param[in]    self       prepared transfer
   !> @param[in]    field_old  field on the source points (nold, ncol)
   !> @param[inout] field_new  field on the targets (ntarget, ncol)
   !> @param[out]   error      not prepared, prepared with `check_only`, shapes that do not
   !>                          match, failed interpolation, or incomplete coverage without
   !>                          `covered`
   !> @param[out]   covered    targets with a value (ntarget)
   !> @param[out]   report     coverage and methods
   subroutine transfer_apply(self, field_old, field_new, error, covered, report)
      !> Prepared transfer
      class(grid_transfer_type), intent(in) :: self
      !> Field on the source points
      real(wp), intent(in) :: field_old(:, :)
      !> Field on the targets
      real(wp), intent(inout) :: field_new(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Targets with a value
      logical, allocatable, intent(out), optional :: covered(:)
      !> Coverage and methods
      type(grid_transfer_report_type), intent(out), optional :: report
      real(wp), allocatable :: g(:, :)
      logical, allocatable :: ignored(:)
      integer :: i

      if (self%nold < 0) then
         call fatal_error(error, "grid transfer: not prepared")
         return
      end if
      if (.not. self%fitted) then
         call fatal_error(error, "grid transfer: prepared for its coverage only")
         return
      end if
      if (size(field_old, 1) /= self%nold) then
         call fatal_error(error, "grid transfer: the field does not match the source layout")
         return
      end if
      if (size(field_new, 1) /= size(self%copy_from) .or. size(field_new, 2) /= size(field_old, 2)) then
         call fatal_error(error, "grid transfer: the target field does not match the target layout")
         return
      end if
      if (.not. present(covered) .and. .not. all(self%copy_from > 0 .or. self%interpolated)) then
         call fatal_error(error, "grid transfer: the source does not cover every target ("//self%reason//")")
         return
      end if
      ! Interpolated values are written only after the interpolation succeeded
      if (any(self%interpolated)) then
         select case (self%method)
         case ("fourier")
            call transfer_fourier(self%source, field_old, self%y, g, error, nthreads=self%nthreads)
            if (allocated(error)) return
            ! Values of the interpolated targets, in order
            field_new(pack([(i, i=1, size(self%interpolated))], self%interpolated), :) = g
         case ("trilinear")
            call trilinear(self%source, field_old, self%y, g, threads(self%nthreads))
            field_new(pack([(i, i=1, size(self%interpolated))], self%interpolated), :) = g
         case ("points")
            call apply_points(self, field_old, g, error)
            if (allocated(error)) return
            field_new(pack([(i, i=1, size(self%interpolated))], self%interpolated), :) = g
         case default
            ! Spline: values at every target
            call transfer_radial(self%source%r, field_old, self%r, g, ignored, error)
            if (allocated(error)) return
            do i = 1, size(self%interpolated)
               if (self%interpolated(i)) field_new(i, :) = g(i, :)
            end do
         end select
      end if
      do i = 1, size(self%copy_from)
         if (self%copy_from(i) > 0) field_new(i, :) = field_old(self%copy_from(i), :)
      end do
      if (present(covered)) covered = self%copy_from > 0 .or. self%interpolated
      if (present(report)) call self%get_report(report)
   end subroutine transfer_apply

   !> Coverage and methods of a prepared transfer
   !>
   !> @param[in]  self    prepared transfer
   !> @param[out] report  coverage and methods
   subroutine transfer_get_report(self, report)
      !> Prepared transfer
      class(grid_transfer_type), intent(in) :: self
      !> Coverage and methods
      type(grid_transfer_report_type), intent(out) :: report

      report%method = ""
      report%reason = ""
      if (self%nold < 0) then
         report%reason = "not prepared"
         return
      end if
      report%ncopied = count(self%copy_from > 0)
      report%ninterpolated = count(self%interpolated)
      report%nuncovered = size(self%copy_from) - report%ncopied - report%ninterpolated
      report%nexcluded = count(self%excluded)
      report%nduplicate = self%nduplicate
      if (report%ncopied > 0) report%method = "copy"
      if (report%ninterpolated > 0) then
         if (report%ncopied > 0) report%method = report%method//"+"
         report%method = report%method//self%method
      end if
      if (report%nuncovered > 0) report%reason = self%reason
   end subroutine transfer_get_report

   !> Release the prepared data
   !>
   !> @param[inout] self  transfer
   subroutine transfer_destroy(self)
      !> Transfer
      class(grid_transfer_type), intent(inout) :: self
      type(moist_math_grid_layout_type) :: empty

      self%nold = -1
      if (allocated(self%copy_from)) deallocate (self%copy_from)
      if (allocated(self%interpolated)) deallocate (self%interpolated)
      if (allocated(self%excluded)) deallocate (self%excluded)
      if (allocated(self%method)) deallocate (self%method)
      if (allocated(self%reason)) deallocate (self%reason)
      if (allocated(self%y)) deallocate (self%y)
      if (allocated(self%r)) deallocate (self%r)
      if (allocated(self%reach)) deallocate (self%reach)
      if (allocated(self%stencil_offset)) deallocate (self%stencil_offset)
      if (allocated(self%stencil_point)) deallocate (self%stencil_point)
      if (allocated(self%stencil_weight)) deallocate (self%stencil_weight)
      call self%tree%destroy()
      self%stored = .false.
      self%fitted = .true.
      self%nduplicate = 0
      self%source = empty
   end subroutine transfer_destroy

   !> Check that a layout describes points a transfer can use, and count them
   !>
   !> - The layout's own `validate`, with coordinates required of point layouts
   !>
   !> @param[in]  layout  layout
   !> @param[in]  label   "source" or "target", for messages
   !> @param[out] n       points
   !> @param[out] error   what is missing
   subroutine described(layout, label, n, error)
      !> Layout
      type(moist_math_grid_layout_type), intent(in) :: layout
      !> Role of the layout
      character(len=*), intent(in) :: label
      !> Points
      integer, intent(out) :: n
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: reason

      call layout%validate(n, reason, coordinates=.true.)
      if (len(reason) > 0) call fatal_error(error, "grid transfer: the "//label//" layout: "//reason)
   end subroutine described

   !> Targets of a layout in the frame of the source
   !>
   !> @param[in]  new  target layout
   !> @param[out] ys   target points in the source frame, bohr (3, n)
   !> @param[in]  map  rigid map into the source frame; without it the frames agree
   pure subroutine target_points(new, ys, map)
      !> Target layout
      type(moist_math_grid_layout_type), intent(in) :: new
      !> Target points in the source frame
      real(wp), allocatable, intent(out) :: ys(:, :)
      !> Rigid map
      type(rigid_fit_type), intent(in), optional :: map
      real(wp), allocatable :: y(:, :)

      call new%points(y)
      if (present(map)) then
         call map%map(y, ys)
      else
         call move_alloc(y, ys)
      end if
   end subroutine target_points

   !> Radial nodes that coincide with exactly one source node
   !>
   !> @param[in]    r_old       increasing source nodes, bohr
   !> @param[in]    r_new       target nodes, bohr
   !> @param[inout] copy_from   source node of each copied target, 0 for none
   !> @param[out]   nduplicate  targets at two or more source nodes
   pure subroutine copy_radial(r_old, r_new, copy_from, nduplicate)
      !> Source nodes
      real(wp), intent(in) :: r_old(:)
      !> Target nodes
      real(wp), intent(in) :: r_new(:)
      !> Source node of each copied target
      integer, intent(inout) :: copy_from(:)
      !> Targets at coinciding source nodes
      integer, intent(out) :: nduplicate
      integer :: i, lo, hi, mid, first, last

      nduplicate = 0
      do i = 1, size(r_new)
         ! Last source node at or below the target, 0 for none
         lo = 0
         hi = size(r_old) + 1
         do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (r_old(mid) <= r_new(i)) then
               lo = mid
            else
               hi = mid
            end if
         end do
         ! The source nodes within the tolerance form one run around that position
         first = lo + 1
         do while (first > 1)
            if (r_new(i) - r_old(first - 1) > copy_tolerance) exit
            first = first - 1
         end do
         last = lo
         do while (last < size(r_old))
            if (r_old(last + 1) - r_new(i) > copy_tolerance) exit
            last = last + 1
         end do
         if (last - first == 0) then
            copy_from(i) = first
         else if (last - first > 0) then
            nduplicate = nduplicate + 1
         end if
      end do
   end subroutine copy_radial

   !> Targets within the tolerance of exactly one point of a Cartesian lattice inside its box
   !>
   !> - The nearest lattice point is the nearest along each axis, so a second point within
   !>   the tolerance puts one of its six axis neighbours within it too
   !>
   !> @param[in]    box         Cartesian layout
   !> @param[in]    ys          targets in the frame of the box, bohr (3, n)
   !> @param[inout] copy_from   lattice point of each copied target, x fastest; 0 for none
   !> @param[out]   nduplicate  targets within the tolerance of two or more lattice points
   pure subroutine copy_lattice(box, ys, copy_from, nduplicate)
      !> Cartesian layout
      type(moist_math_grid_layout_type), intent(in) :: box
      !> Targets
      real(wp), intent(in) :: ys(:, :)
      !> Lattice point of each copied target
      integer, intent(inout) :: copy_from(:)
      !> Targets at two or more lattice points
      integer, intent(out) :: nduplicate
      real(wp) :: u(3), d(3), e(3)
      integer :: i, k, side, ijk(3)
      logical :: unique

      nduplicate = 0
      do i = 1, size(ys, 2)
         u = (ys(:, i) - box%origin)/box%spacing - 0.5_wp
         ! Far beyond the box: no lattice point, and no integer overflow below
         if (any(u < -1.0_wp) .or. any(u > real(box%shape, wp))) cycle
         ijk = nint(u)
         d = (u - real(ijk, wp))*box%spacing
         if (norm2(d) > copy_tolerance) cycle
         ijk = ijk + 1
         if (any(ijk < 1) .or. any(ijk > box%shape)) cycle
         unique = .true.
         do k = 1, 3
            do side = -1, 1, 2
               if (ijk(k) + side < 1 .or. ijk(k) + side > box%shape(k)) cycle
               e = d
               e(k) = d(k) - real(side, wp)*box%spacing
               if (norm2(e) <= copy_tolerance) unique = .false.
            end do
         end do
         if (unique) then
            copy_from(i) = ijk(1) + box%shape(1)*((ijk(2) - 1) + box%shape(2)*(ijk(3) - 1))
         else
            nduplicate = nduplicate + 1
         end if
      end do
   end subroutine copy_lattice

   !> Targets within the tolerance of exactly one source point
   !>
   !> - Targets beyond the reach of every source point are skipped without a search
   !>
   !> @param[in]    tree        k-d tree of the source points, with their reach
   !> @param[in]    ys          targets in the source frame, bohr (3, n)
   !> @param[in]    nthreads    threads
   !> @param[inout] copy_from   source point of each copied target, 0 for none
   !> @param[out]   nduplicate  targets at two or more source points
   !> @param[out]   near        targets within reach and not copied (n)
   subroutine copy_points(tree, ys, nthreads, copy_from, nduplicate, near)
      !> Tree of the source points
      type(moist_kdtree_type), intent(in) :: tree
      !> Targets
      real(wp), intent(in) :: ys(:, :)
      !> Threads
      integer, intent(in) :: nthreads
      !> Source point of each copied target
      integer, intent(inout) :: copy_from(:)
      !> Targets at coinciding source points
      integer, intent(out) :: nduplicate
      !> Targets within reach and not copied
      logical, intent(out) :: near(:)
      integer :: nearest(2), i
      real(wp) :: dist2(2), tol2

      nduplicate = 0
      near = .false.
      if (size(tree%x, 2) == 0) return
      tol2 = copy_tolerance**2
      !$omp parallel do schedule(dynamic, 256) num_threads(nthreads) private(nearest, dist2) reduction(+:nduplicate)
      do i = 1, size(ys, 2)
         if (.not. tree%within_reach(ys(:, i))) cycle
         call tree%knn(ys(:, i), 2, nearest, dist2)
         if (dist2(1) <= tol2 .and. dist2(2) > tol2) then
            copy_from(i) = nearest(1)
         else
            if (dist2(1) <= tol2) nduplicate = nduplicate + 1
            near(i) = .true.
         end if
      end do
   end subroutine copy_points

   !> Fit the stencils of the targets near a point layout
   !>
   !> - Blocks of targets are fitted in parallel; the stencils are stored while their
   !>   memory stays within the `stencil_memory` option, buffers included
   !> - Stored stencils keep neither the targets nor the source tree; otherwise the targets
   !>   are kept and the stencils are fitted again when applied
   !> - Sets the interpolated targets and the reason for the others
   !>
   !> @param[inout] self       transfer with the tree and reach of the source points
   !> @param[in]    ys         targets in the source frame, bohr (3, n)
   !> @param[in]    candidate  targets to fit (n)
   subroutine fit_points(self, ys, candidate)
      !> Transfer
      type(grid_transfer_type), intent(inout) :: self
      !> Targets
      real(wp), intent(in) :: ys(:, :)
      !> Targets to fit
      logical, intent(in) :: candidate(:)
      integer, allocatable :: idx(:), block_point(:, :), block_k(:), block_status(:), point(:)
      real(wp), allocatable :: block_weight(:, :), weight(:)
      integer :: ncause(stencil_far:stencil_amplification), nfit, nentry, first, last, i, ib, k, nt
      integer :: ifar, budget
      integer(int64) :: entry_bytes

      idx = pack([(i, i=1, size(candidate))], candidate)
      ! Stencil entries the memory option allows: a point index and a weight each
      entry_bytes = int(storage_size(1.0_wp)/8 + storage_size(1)/8, int64)
      budget = int(min(self%stencil_memory/entry_bytes, int(huge(1), int64)))
      nt = threads(self%nthreads)
      allocate (block_point(k_last, stencil_block), block_weight(k_last, stencil_block), &
         & block_k(stencil_block), block_status(stencil_block))
      allocate (self%stencil_offset(size(idx) + 1), point(0), weight(0))
      self%stencil_offset(1) = 1
      self%stored = .true.
      ncause = 0
      ! Targets not near any source point were left out before the fit
      ncause(stencil_far) = count(self%copy_from == 0 .and. .not. candidate)
      nfit = 0
      nentry = 0
      do ib = 1, (size(idx) + stencil_block - 1)/stencil_block
         first = (ib - 1)*stencil_block + 1
         last = min(ib*stencil_block, size(idx))
         !$omp parallel do schedule(dynamic, 64) num_threads(nt)
         do i = first, last
            call point_stencil(self%tree, self%reach, ys(:, idx(i)), block_point(:, i - first + 1), &
               & block_weight(:, i - first + 1), block_k(i - first + 1), block_status(i - first + 1))
         end do
         do i = first, last
            if (block_status(i - first + 1) /= stencil_fitted) then
               ncause(block_status(i - first + 1)) = ncause(block_status(i - first + 1)) + 1
               cycle
            end if
            nfit = nfit + 1
            self%interpolated(idx(i)) = .true.
            if (.not. self%stored) cycle
            k = block_k(i - first + 1)
            if (nentry + k > budget) then
               ! Over the budget: the stencils are fitted again when applied
               self%stored = .false.
               deallocate (point, weight)
               cycle
            end if
            ! Growth stays within the budget, which holds nentry + k
            if (nentry + k > size(point)) call grow(min(max(2*size(point), nentry + k, 1024), budget))
            point(nentry + 1:nentry + k) = block_point(1:k, i - first + 1)
            weight(nentry + 1:nentry + k) = block_weight(1:k, i - first + 1)
            nentry = nentry + k
            self%stencil_offset(nfit + 1) = nentry + 1
         end do
      end do
      if (self%stored) then
         ! The buffers become the cache; entries past nentry are not referenced
         self%stencil_offset = self%stencil_offset(:nfit + 1)
         call move_alloc(point, self%stencil_point)
         call move_alloc(weight, self%stencil_weight)
         call self%tree%destroy()
         if (allocated(self%reach)) deallocate (self%reach)
      else
         deallocate (self%stencil_offset)
         self%y = ys(:, pack([(i, i=1, size(candidate))], self%interpolated))
      end if
      ! The reason names the most frequent cause among the targets left uncovered
      ifar = maxloc(ncause, 1) + stencil_far - 1
      select case (ifar)
      case (stencil_conditioning)
         self%reason = "source neighbourhood flat or too small"
      case (stencil_amplification)
         self%reason = "stencil amplification above its bound"
      case default
         self%reason = "outside the source grid"
      end select
   contains
      !> Enlarge the stencil buffers
      !>
      !> @param[in] n  new capacity, at least the entries so far
      subroutine grow(n)
         !> New capacity
         integer, intent(in) :: n
         integer, allocatable :: point_new(:)
         real(wp), allocatable :: weight_new(:)
         allocate (point_new(n), weight_new(n))
         point_new(:nentry) = point(:nentry)
         weight_new(:nentry) = weight(:nentry)
         call move_alloc(point_new, point)
         call move_alloc(weight_new, weight)
      end subroutine grow
   end subroutine fit_points

   !> Mark the first target near a point layout that gets a stencil, if any
   !>
   !> - Blocks of targets are fitted in parallel, until a block has an accepted stencil
   !> - Sets the reason when no target gets one
   !>
   !> @param[inout] self       transfer with the tree and reach of the source points
   !> @param[in]    ys         targets in the source frame, bohr (3, n)
   !> @param[in]    candidate  targets to try (n)
   subroutine fit_first_point(self, ys, candidate)
      !> Transfer
      type(grid_transfer_type), intent(inout) :: self
      !> Targets
      real(wp), intent(in) :: ys(:, :)
      !> Targets to try
      logical, intent(in) :: candidate(:)
      integer, allocatable :: idx(:), block_status(:)
      integer :: ncause(stencil_far:stencil_amplification), nearest(k_last), first, last, i, ib, k
      real(wp) :: weight(k_last)

      idx = pack([(i, i=1, size(candidate))], candidate)
      allocate (block_status(stencil_block))
      ncause = 0
      ncause(stencil_far) = count(self%copy_from == 0 .and. .not. candidate)
      do ib = 1, (size(idx) + stencil_block - 1)/stencil_block
         first = (ib - 1)*stencil_block + 1
         last = min(ib*stencil_block, size(idx))
         !$omp parallel do schedule(dynamic, 64) num_threads(threads(self%nthreads)) private(nearest, weight, k)
         do i = first, last
            call point_stencil(self%tree, self%reach, ys(:, idx(i)), nearest, weight, k, block_status(i - first + 1))
         end do
         do i = first, last
            if (block_status(i - first + 1) == stencil_fitted) then
               self%interpolated(idx(i)) = .true.
               return
            end if
            ncause(block_status(i - first + 1)) = ncause(block_status(i - first + 1)) + 1
         end do
      end do
      select case (maxloc(ncause, 1) + stencil_far - 1)
      case (stencil_conditioning)
         self%reason = "source neighbourhood flat or too small"
      case (stencil_amplification)
         self%reason = "stencil amplification above its bound"
      case default
         self%reason = "outside the source grid"
      end select
   end subroutine fit_first_point

   !> Values of the stencils of the interpolated targets
   !>
   !> - Stencils that are not stored are fitted again; a target whose fit is not accepted
   !>   again is an error
   !>
   !> @param[in]  self   transfer with point stencils
   !> @param[in]  f      field on the source points (nold, ncol)
   !> @param[out] g      values at the interpolated targets (ninterpolated, ncol)
   !> @param[out] error  a stencil fitted again was not accepted
   subroutine apply_points(self, f, g, error)
      !> Transfer
      type(grid_transfer_type), intent(in) :: self
      !> Field on the source points
      real(wp), intent(in) :: f(:, :)
      !> Values at the interpolated targets
      real(wp), allocatable, intent(out) :: g(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer :: nearest(k_last), j, p, k, status, nfailed
      real(wp) :: weight(k_last)

      if (self%stored) then
         allocate (g(size(self%stencil_offset) - 1, size(f, 2)), source=0.0_wp)
         !$omp parallel do schedule(static) num_threads(threads(self%nthreads)) private(p)
         do j = 1, size(g, 1)
            do p = self%stencil_offset(j), self%stencil_offset(j + 1) - 1
               g(j, :) = g(j, :) + self%stencil_weight(p)*f(self%stencil_point(p), :)
            end do
         end do
      else
         allocate (g(size(self%y, 2), size(f, 2)), source=0.0_wp)
         ! The fit is deterministic and the transfer keeps its geometry: the stencils of prepare
         nfailed = 0
         !$omp parallel do schedule(dynamic, 64) num_threads(threads(self%nthreads)) private(nearest, weight, k, status, p) &
         !$omp& reduction(+:nfailed)
         do j = 1, size(g, 1)
            call point_stencil(self%tree, self%reach, self%y(:, j), nearest, weight, k, status)
            if (status /= stencil_fitted) then
               nfailed = nfailed + 1
               cycle
            end if
            do p = 1, k
               g(j, :) = g(j, :) + weight(p)*f(nearest(p), :)
            end do
         end do
         if (nfailed > 0) then
            call fatal_error(error, "grid transfer: a point stencil fitted again was not accepted")
            deallocate (g)
         end if
      end if
   end subroutine apply_points

   !> Trilinear interpolation on a Cartesian lattice
   !>
   !> - The eight lattice points around a target, weighted by the fractions of the cell
   !> - Within half a cell of the box faces, outside the outermost lattice points, the
   !>   value of the nearest face of lattice points is taken
   !>
   !> @param[in]  box       Cartesian layout
   !> @param[in]  f         field on the lattice, x fastest (npoint, ncol)
   !> @param[in]  y         targets in the frame of the box, bohr (3, n)
   !> @param[out] g         values at the targets (n, ncol)
   !> @param[in]  nthreads  threads
   subroutine trilinear(box, f, y, g, nthreads)
      !> Cartesian layout
      type(moist_math_grid_layout_type), intent(in) :: box
      !> Field on the lattice
      real(wp), intent(in) :: f(:, :)
      !> Targets
      real(wp), intent(in) :: y(:, :)
      !> Values at the targets
      real(wp), allocatable, intent(out) :: g(:, :)
      !> Threads
      integer, intent(in) :: nthreads
      real(wp) :: u(3), t(3), w
      integer :: lo(3), hi(3), i, a, b, c, ijk(3)

      allocate (g(size(y, 2), size(f, 2)), source=0.0_wp)
      !$omp parallel do schedule(static) num_threads(nthreads) private(u, t, w, lo, hi, a, b, c, ijk)
      do i = 1, size(y, 2)
         ! Lattice coordinates from 0, clamped to the outermost lattice points
         u = (y(:, i) - box%origin)/box%spacing - 0.5_wp
         u = max(0.0_wp, min(u, real(box%shape - 1, wp)))
         lo = min(int(u), max(box%shape - 2, 0))
         t = u - real(lo, wp)
         hi = min(lo + 1, box%shape - 1)
         do c = 0, 1
            do b = 0, 1
               do a = 0, 1
                  w = merge(t(1), 1.0_wp - t(1), a == 1)*merge(t(2), 1.0_wp - t(2), b == 1) &
                     & *merge(t(3), 1.0_wp - t(3), c == 1)
                  if (w == 0.0_wp) cycle
                  ijk = merge(hi, lo, [a, b, c] == 1) + 1
                  g(i, :) = g(i, :) + w*f(ijk(1) + box%shape(1)*((ijk(2) - 1) + box%shape(2)*(ijk(3) - 1)), :)
               end do
            end do
         end do
      end do
   end subroutine trilinear

   !> Thread count: `nthreads` when positive, else the OpenMP default
   !>
   !> @param[in] nthreads  requested threads
   integer function threads(nthreads) result(nt)
      !> Requested threads
      integer, intent(in) :: nthreads

      nt = 1
      !$ nt = max(1, omp_get_max_threads())
      if (nthreads > 0) nt = nthreads
   end function threads

end module moist_math_grid_transfer
