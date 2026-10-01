!> Size-adjusted Becke partitioning
module moist_math_grid_3d_partition_becke
   use mctc_env, only: wp
   use moist_math_grid_3d_partition_common, only: pair_partition_weights
   implicit none(type, external)
   private

   public :: becke_partition_weights, becke_cell, default_stiffness

   !> Recommended number of Becke polynomial iterations
   integer, parameter :: default_stiffness = 3

contains

   !> Normalized Becke weights of one atom at a batch of points
   !>
   !> Unchecked shapes and valid atomic numbers/owner required; no quadrature factors
   !>
   !> @param[in]  owner     Owner index in [1,nat]
   !> @param[in]  points    Points in bohr, shape (3,npts)
   !> @param[in]  xyz       Atomic positions in bohr, shape (3,nat)
   !> @param[in]  numbers   Atomic numbers, shape (nat)
   !> @param[out] w         Owner weights in [0,1], shape (npts)
   !> @param[in]  stiffness Polynomial iterations >= 1; default 3
   subroutine becke_partition_weights(owner, points, xyz, numbers, w, stiffness)
      !> Owner index
      integer, intent(in) :: owner
      !> Sample points
      real(wp), intent(in) :: points(:, :)
      !> Atomic positions
      real(wp), intent(in) :: xyz(:, :)
      !> Atomic numbers
      integer, intent(in) :: numbers(:)
      !> Owner weights
      real(wp), intent(out) :: w(:)
      !> Polynomial iterations
      integer, intent(in), optional :: stiffness

      integer :: k

      k = default_stiffness
      if (present(stiffness)) k = stiffness
      call pair_partition_weights(owner, points, xyz, numbers, w, cell)

   contains

      !> Becke cell with the configured iteration count
      !>
      !> @param[in] x Size-adjusted elliptic coordinate
      pure function cell(x) result(s)
         !> Pair coordinate
         real(wp), intent(in) :: x
         !> Cell weight
         real(wp) :: s
         s = becke_cell(x, k)
      end function cell

   end subroutine becke_partition_weights

   !> Iterated Becke cell, with a stable positive tail near x = 1
   !>
   !> A. D. Becke, J. Chem. Phys. 88, 2547 (1988), doi:10.1063/1.454033
   !> Iterating the weight t instead of x gives t <- t**2*(3-2*t)
   !>
   !> @param[in] x Size-adjusted elliptic coordinate
   !> @param[in] k Polynomial iterations >= 1
   pure function becke_cell(x, k) result(s)
      !> Pair coordinate
      real(wp), intent(in) :: x
      !> Polynomial iterations
      integer, intent(in) :: k
      !> Cell weight in [0,1]
      real(wp) :: s

      integer :: i

      s = 0.5_wp*(1.0_wp - min(1.0_wp, abs(x)))
      do i = 1, k
         s = s*s*(3.0_wp - 2.0_wp*s)
      end do
      if (x < 0.0_wp) s = 1.0_wp - s
   end function becke_cell

end module moist_math_grid_3d_partition_becke
