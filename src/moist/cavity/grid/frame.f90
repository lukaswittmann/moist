!> Grid Jacobians for the forward-mode terms of a cavity's grid components
!>
!> The rotation of a component follows the atoms, and so may its weight; a
!> cavity knows how its own quantities change with a rotation or weight of
!> the owning component, and this type spreads those changes over the atoms
!>
!> - a quantity of a point (or component) of component c that changes by
!>   t . omega when the component's grid turns by omega gains
!>   dq/dc_j(i) = sum_k g(m, k, rot(c)) t(k), m = i + 3*(j - 1): the frame
!>   reaches every atom, so the forward arrays become dense
!> - a quantity that changes by s per unit weight of its component gains
!>   s gw(m, wclass(c)) (`add_weight`); nothing while no weight depends on
!>   the geometry (`has_weight` false)
!> - components sharing a weight share one column of gw (`wclass`), like
!>   components sharing a rotation share one of g
!> - built by [[new_grid_frame]] once per gradient call; the add procedures
!>   are pure, so a parallel loop over the points can give every point its
!>   term while its arrays are in cache
!> - per-sphere and total quantities take the sum of t over a component's
!>   points once (`add_spheres`, `add_total`), not one term per point
module moist_cavity_grid_frame
   use mctc_env, only: wp
   use moist_cavity_grid_model_type, only: grid_model_type
   use moist_cavity_grid_type, only: cavity_grid_type
   implicit none(type, external)
   private

   public :: grid_frame_type, new_grid_frame

   !> Rotation and weight Jacobians of the grid components of a cavity
   type :: grid_frame_type
      !> d omega_r(k)/d c_j(i) at m = i + 3*(j - 1) (3*nsph, 3, nrot)
      real(wp), allocatable :: g(:, :, :)
      !> The same interleaved by the axis k, gk(k, m, r) (3, 3*nsph, nrot)
      real(wp), allocatable :: gk(:, :, :)
      !> Distinct rotation of every component (ncomp)
      integer, allocatable :: rot(:)
      !> Sphere of every component (ncomp)
      integer, allocatable :: sphere(:)
      !> Whether any component weight depends on the geometry
      logical :: has_weight = .false.
      !> Distinct weight of every component (ncomp), with `has_weight`
      integer, allocatable :: wclass(:)
      !> dw/d c_j(i) of every distinct weight (3, nsph, nclass), with `has_weight`
      real(wp), allocatable :: gw(:, :, :)
   contains
      !> Rotation term of a scalar of one point, q(3, nsph)
      procedure :: add_scalar => frame_add_scalar
      !> Rotation term of a 3-vector of one point, q(3, 3, nsph)
      procedure :: add_vector => frame_add_vector
      !> Rotation term of per-sphere sums, q(3, nsph, nsph) with the owner second
      procedure :: add_spheres => frame_add_spheres
      !> Rotation term of a total, q(3, nsph)
      procedure :: add_total => frame_add_total
      !> Weight term of a scalar of one point, q(3, nsph)
      procedure :: add_weight => frame_add_weight
   end type grid_frame_type

contains

   !> Grid Jacobians of the current structure
   !>
   !> - plain grid or no model: zero rotation Jacobian (nothing turns)
   !> - otherwise the model's `jacobian`, once per gradient call
   !>
   !> @param[out] frame Grid Jacobians
   !> @param[in]  grid  Assembled components of the cavity
   !> @param[in]  model Updated grid model
   subroutine new_grid_frame(frame, grid, model)
      !> Grid Jacobians
      type(grid_frame_type), intent(out) :: frame
      !> Components
      type(cavity_grid_type), intent(in) :: grid
      !> Grid model
      class(grid_model_type), intent(in), optional :: model

      !> d omega_r(k)/d c_j(i) as (k, i, j, r)
      real(wp), allocatable :: jac(:, :, :, :)
      !> d w/d c_j(i) of every distinct weight as (i, j, class)
      real(wp), allocatable :: jac_weight(:, :, :)
      !> Distinct weight of every component
      integer, allocatable :: weight_class(:)
      integer :: irot

      if (grid%plain .or. .not. present(model)) then
         allocate (frame%rot(grid%ncomp), source=1)
         frame%sphere = grid%comp_sphere
         allocate (frame%g(3*grid%nsph, 3, 1), frame%gk(3, 3*grid%nsph, 1), source=0.0_wp)
         return
      end if
      frame%rot = grid%comp_rot
      frame%sphere = grid%comp_sphere
      call model%jacobian(jac, jac_weight, weight_class)
      ! jac(k, i, j, r) is gk(k, m, r) in memory; g is its transpose per r
      frame%gk = reshape(jac, [3, 3*grid%nsph, size(jac, 4)])
      allocate (frame%g(3*grid%nsph, 3, size(jac, 4)))
      do irot = 1, size(jac, 4)
         frame%g(:, :, irot) = transpose(frame%gk(:, :, irot))
      end do
      ! A reduced grid holds one constant weight per sphere
      if (allocated(jac_weight) .and. .not. grid%reduced) then
         frame%has_weight = .true.
         call move_alloc(jac_weight, frame%gw)
         call move_alloc(weight_class, frame%wclass)
      end if

   end subroutine new_grid_frame

   !> Rotation term of a scalar of one point
   !>
   !> q(i, j) += sum_k g(i + 3*(j - 1), k, rot(icomp)) t(k)
   !>
   !> @param[in]     self  Grid Jacobians
   !> @param[in]     icomp Component of the point
   !> @param[in]     t     Change of the scalar per rotation of the component's grid (3)
   !> @param[in,out] q     Nuclear derivatives of the scalar (3, nsph), contiguous
   pure subroutine frame_add_scalar(self, icomp, t, q)
      !> Grid Jacobians
      class(grid_frame_type), intent(in) :: self
      !> Component
      integer, intent(in) :: icomp
      !> Change per rotation
      real(wp), intent(in) :: t(3)
      !> Nuclear derivatives
      real(wp), intent(inout), contiguous :: q(:, :)

      integer :: irot

      irot = self%rot(icomp)
      call axpy_axes(size(q), self%g(:, 1, irot), self%g(:, 2, irot), self%g(:, 3, irot), t, q)

   end subroutine frame_add_scalar

   !> Rotation term of a 3-vector of one point
   !>
   !> q(c, i, j) += sum_k gk(k, i + 3*(j - 1), rot(icomp)) t(k, c)
   !>
   !> @param[in]     self  Grid Jacobians
   !> @param[in]     icomp Component of the point
   !> @param[in]     t     Change of component c per rotation about axis k (3, 3)
   !> @param[in,out] q     Nuclear derivatives of the vector (3, 3, nsph), contiguous
   pure subroutine frame_add_vector(self, icomp, t, q)
      !> Grid Jacobians
      class(grid_frame_type), intent(in) :: self
      !> Component
      integer, intent(in) :: icomp
      !> Change per rotation
      real(wp), intent(in) :: t(3, 3)
      !> Nuclear derivatives
      real(wp), intent(inout), contiguous :: q(:, :, :)

      call matvec_axes(size(q)/3, self%gk(:, :, self%rot(icomp)), t, q)

   end subroutine frame_add_vector

   !> Rotation term of per-sphere sums
   !>
   !> q(i, a, j) += sum_k g(i + 3*(j - 1), k, rot(c)) s(k, c), a the sphere of
   !> component c
   !>
   !> @param[in]     self Grid Jacobians
   !> @param[in]     s    Change per rotation summed over each component's points (3, ncomp)
   !> @param[in,out] q    Nuclear derivatives of the sums (3, nsph, nsph), owner second
   pure subroutine frame_add_spheres(self, s, q)
      !> Grid Jacobians
      class(grid_frame_type), intent(in) :: self
      !> Summed change per rotation
      real(wp), intent(in) :: s(:, :)
      !> Nuclear derivatives
      real(wp), intent(inout) :: q(:, :, :)

      integer :: icomp, isph, jsph, irot, m

      do jsph = 1, size(q, 3)
         m = 3*(jsph - 1)
         do icomp = 1, size(s, 2)
            irot = self%rot(icomp)
            isph = self%sphere(icomp)
            q(:, isph, jsph) = q(:, isph, jsph) + self%g(m + 1:m + 3, 1, irot)*s(1, icomp) &
               & + self%g(m + 1:m + 3, 2, irot)*s(2, icomp) + self%g(m + 1:m + 3, 3, irot)*s(3, icomp)
         end do
      end do

   end subroutine frame_add_spheres

   !> Rotation term of a total
   !>
   !> q(i, j) += sum_c sum_k g(i + 3*(j - 1), k, rot(c)) s(k, c), summed
   !> per distinct rotation first
   !>
   !> @param[in]     self Grid Jacobians
   !> @param[in]     s    Change per rotation summed over each component's points (3, ncomp)
   !> @param[in,out] q    Nuclear derivatives of the total (3, nsph), contiguous
   pure subroutine frame_add_total(self, s, q)
      !> Grid Jacobians
      class(grid_frame_type), intent(in) :: self
      !> Summed change per rotation
      real(wp), intent(in) :: s(:, :)
      !> Nuclear derivatives
      real(wp), intent(inout), contiguous :: q(:, :)

      !> Change per rotation summed over the components of each rotation
      real(wp) :: srot(3, size(self%g, 3))
      integer :: icomp, irot

      srot = 0.0_wp
      do icomp = 1, size(s, 2)
         srot(:, self%rot(icomp)) = srot(:, self%rot(icomp)) + s(:, icomp)
      end do
      do irot = 1, size(srot, 2)
         call axpy_axes(size(q), self%g(:, 1, irot), self%g(:, 2, irot), self%g(:, 3, irot), &
            & srot(:, irot), q)
      end do

   end subroutine frame_add_total

   !> Weight term of a scalar of one point
   !>
   !> q(i, j) += s gw(i, j, wclass(icomp)); nothing unless `has_weight`
   !>
   !> @param[in]     self  Grid Jacobians
   !> @param[in]     icomp Component of the point
   !> @param[in]     s     Change of the scalar per unit weight of the component
   !> @param[in,out] q     Nuclear derivatives of the scalar (3, nsph), contiguous
   pure subroutine frame_add_weight(self, icomp, s, q)
      !> Grid Jacobians
      class(grid_frame_type), intent(in) :: self
      !> Component
      integer, intent(in) :: icomp
      !> Change per unit weight
      real(wp), intent(in) :: s
      !> Nuclear derivatives
      real(wp), intent(inout), contiguous :: q(:, :)

      if (.not. self%has_weight) return
      call axpy_flat(size(q), s, self%gw(:, :, self%wclass(icomp)), q)

   end subroutine frame_add_weight

   !> q += s g over a flat array
   !>
   !> @param[in]     n Length
   !> @param[in]     s Factor
   !> @param[in]     g Column (n)
   !> @param[in,out] q Accumulator (n)
   pure subroutine axpy_flat(n, s, g, q)
      !> Length
      integer, intent(in) :: n
      !> Factor
      real(wp), intent(in) :: s
      !> Column
      real(wp), intent(in) :: g(n)
      !> Accumulator
      real(wp), intent(inout) :: q(n)

      integer :: m

      do m = 1, n
         q(m) = q(m) + s*g(m)
      end do

   end subroutine axpy_flat

   !> q += g1 t(1) + g2 t(2) + g3 t(3) over a flat array
   !>
   !> Explicit shape on purpose: the callers pass contiguous arrays of any
   !> rank, and the single flat loop vectorizes (about 1.5x the nested
   !> (3, nsph) loop)
   !>
   !> @param[in]     n  Length
   !> @param[in]     g1 Column of axis 1 (n)
   !> @param[in]     g2 Column of axis 2 (n)
   !> @param[in]     g3 Column of axis 3 (n)
   !> @param[in]     t  Weights of the axes (3)
   !> @param[in,out] q  Accumulator (n)
   pure subroutine axpy_axes(n, g1, g2, g3, t, q)
      !> Length
      integer, intent(in) :: n
      !> Column of axis 1
      real(wp), intent(in) :: g1(n)
      !> Column of axis 2
      real(wp), intent(in) :: g2(n)
      !> Column of axis 3
      real(wp), intent(in) :: g3(n)
      !> Weights of the axes
      real(wp), intent(in) :: t(3)
      !> Accumulator
      real(wp), intent(inout) :: q(n)

      real(wp) :: t1, t2, t3
      integer :: m

      t1 = t(1)
      t2 = t(2)
      t3 = t(3)
      do m = 1, n
         q(m) = q(m) + g1(m)*t1 + g2(m)*t2 + g3(m)*t3
      end do

   end subroutine axpy_axes

   !> q(c, m) += sum_k gk(k, m) t(k, c) over a flat array of 3-vectors
   !>
   !> Explicit shape on purpose, as [[axpy_axes]]; the axis-interleaved
   !> Jacobian streams best here
   !>
   !> @param[in]     n  Number of vectors
   !> @param[in]     gk Jacobian interleaved by axis (3, n)
   !> @param[in]     t  Change of component c per axis k (3, 3)
   !> @param[in,out] q  Accumulator (3, n)
   pure subroutine matvec_axes(n, gk, t, q)
      !> Number of vectors
      integer, intent(in) :: n
      !> Jacobian interleaved by axis
      real(wp), intent(in) :: gk(3, n)
      !> Change per axis
      real(wp), intent(in) :: t(3, 3)
      !> Accumulator
      real(wp), intent(inout) :: q(3, n)

      integer :: m

      do m = 1, n
         q(1, m) = q(1, m) + gk(1, m)*t(1, 1) + gk(2, m)*t(2, 1) + gk(3, m)*t(3, 1)
         q(2, m) = q(2, m) + gk(1, m)*t(1, 2) + gk(2, m)*t(2, 2) + gk(3, m)*t(3, 2)
         q(3, m) = q(3, m) + gk(1, m)*t(1, 3) + gk(2, m)*t(2, 3) + gk(3, m)*t(3, 3)
      end do

   end subroutine matvec_axes

end module moist_cavity_grid_frame
