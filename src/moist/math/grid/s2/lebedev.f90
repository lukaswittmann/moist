!> Lebedev-Laikov angular grid as a concrete S2 quadrature
module moist_math_grid_s2_lebedev
   use mctc_env, only: wp, error_type, fatal_error
   use moist_math_grid_s2_base, only: moist_math_grid_s2_type
   use moist_math_quadrature_lebedev, only: grid_size, get_angular_grid, &
      & lebedev_order_from_num
   implicit none(type, external)
   private

   public :: moist_math_grid_s2_lebedev_type, new_s2_grid_lebedev

   !> Algebraic exactness degrees of the 32 supported Lebedev-Laikov rules
   integer, parameter :: lebedev_degree_table(32) = [ &
      &   3,   5,   7,   9,  11,  13,  15,  17, &
      &  19,  21,  23,  25,  27,  29,  31,  35, &
      &  41,  47,  53,  59,  65,  71,  77,  83, &
      &  89,  95, 101, 107, 113, 119, 125, 131]

   !> Immutable Lebedev-Laikov unit-sphere grid
   type, extends(moist_math_grid_s2_type) :: moist_math_grid_s2_lebedev_type
      !> Order index (1..32) of the selected rule, indexing `grid_size`
      integer :: order = 0
   contains
      procedure :: measure => lebedev_measure
      procedure :: degree => lebedev_degree
      procedure :: destroy => lebedev_grid_destroy
   end type moist_math_grid_s2_lebedev_type

contains

   !> Initialize a Lebedev-Laikov S2 grid, selected either by exact point
   !> count or by required exactness degree.  Exactly one of `npts` / `degree`
   !> must be present.
   !>
   !> @param[out] self    Initialized grid object
   !> @param[out] error   Set on an invalid request or allocation failure
   !> @param[in]  npts    Exact number of nodes; must be one of the 32
   !>                     supported Lebedev sizes
   !> @param[in]  degree  Required exactness degree; the smallest supported
   !>                     rule exact to at least this degree is selected
   subroutine new_s2_grid_lebedev(self, error, npts, degree)
      !> Initialized grid object
      type(moist_math_grid_s2_lebedev_type), intent(out) :: self
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Exact number of nodes (one of the supported Lebedev sizes)
      integer, intent(in), optional :: npts
      !> Required exactness degree
      integer, intent(in), optional :: degree

      !> Order index of the selected rule
      integer :: order

      if (present(npts) .eqv. present(degree)) then
         call fatal_error(error, &
            & "Lebedev S2 grid needs either npts or degree")
         return
      end if

      if (present(npts)) then
         call lebedev_order_from_num(npts, order, error)
         if (allocated(error)) return
      else
         call lebedev_order_from_degree(degree, order, error)
         if (allocated(error)) return
      end if

      self%order = order
      self%npts = grid_size(order)

      allocate (self%points(3, self%npts), self%weights(self%npts))

      call get_angular_grid(order, self%points, self%weights, error)
      if (allocated(error)) then
         call self%destroy()
         return
      end if
   end subroutine new_s2_grid_lebedev

   !> Map a required exactness degree to the order index of the smallest
   !> supported rule that reaches it
   !>
   !> @param[in]  degree  Required exactness degree (non-negative)
   !> @param[out] order   Order index (1..32) of the selected rule
   !> @param[out] error   Set if the degree is negative or beyond the table
   subroutine lebedev_order_from_degree(degree, order, error)
      !> Required exactness degree
      integer, intent(in) :: degree
      !> Order index of the selected rule
      integer, intent(out) :: order
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Loop index over the degree table
      integer :: i

      order = 0
      if (degree < 0) then
         call fatal_error(error, "Lebedev S2 grid degree must be non-negative")
         return
      end if

      do i = 1, size(lebedev_degree_table)
         if (lebedev_degree_table(i) >= degree) then
            order = i
            return
         end if
      end do

      call fatal_error(error, "Requested Lebedev degree exceeds the highest available rule")
   end subroutine lebedev_order_from_degree

   !> Quadrature weight (fraction of the unit sphere) of node i
   !>
   !> @param[in]  self  Grid instance
   !> @param[in]  i     Grid node index
   !> @return     w     Weight of node i
   pure function lebedev_measure(self, i) result(w)
      !> Grid instance
      class(moist_math_grid_s2_lebedev_type), intent(in) :: self
      !> Grid node index
      integer, intent(in) :: i
      !> Weight of node i
      real(wp) :: w

      w = self%weights(i)
   end function lebedev_measure

   !> Algebraic exactness degree of the selected Lebedev rule
   !>
   !> @param[in]  self  Grid instance
   !> @return     d     Exactness degree, 0 for an uninitialized grid
   pure function lebedev_degree(self) result(d)
      !> Grid instance
      class(moist_math_grid_s2_lebedev_type), intent(in) :: self
      !> Exactness degree
      integer :: d

      if (self%order >= 1 .and. self%order <= size(lebedev_degree_table)) then
         d = lebedev_degree_table(self%order)
      else
         d = 0
      end if
   end function lebedev_degree

   !> Release all grid storage; idempotent
   !>
   !> @param[inout] self  Grid instance
   subroutine lebedev_grid_destroy(self)
      !> Grid instance
      class(moist_math_grid_s2_lebedev_type), intent(inout) :: self

      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%weights)) deallocate (self%weights)
      self%npts = 0
      self%order = 0
   end subroutine lebedev_grid_destroy

end module moist_math_grid_s2_lebedev
