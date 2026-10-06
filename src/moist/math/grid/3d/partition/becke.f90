!> Size-adjusted Becke partitioning
module moist_math_grid_3d_partition_becke
   use mctc_env, only: wp
   use moist_math_grid_3d_partition_common, only: pair_partition_weights, partition_cell
   implicit none(type, external)
   private

   public :: becke_partition_weights, becke_cell, default_stiffness

   !> Recommended number of Becke polynomial iterations
   integer, parameter :: default_stiffness = 3

   !> Becke cell with a fixed iteration count
   type, extends(partition_cell) :: becke_partition_cell
      !> Polynomial iterations
      integer :: k = default_stiffness
   contains
      !> Pair-cell weight
      procedure :: eval => becke_partition_cell_eval
   end type becke_partition_cell

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
   !> @param[in]  nthreads  Thread count; absent takes omp_get_max_threads
   subroutine becke_partition_weights(owner, points, xyz, numbers, w, stiffness, nthreads)
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
      !> Thread count
      integer, intent(in), optional :: nthreads

      type(becke_partition_cell) :: cell

      if (present(stiffness)) cell%k = stiffness
      call pair_partition_weights(owner, points, xyz, numbers, w, cell, nthreads=nthreads)
   end subroutine becke_partition_weights

   !> Becke cell with the configured iteration count
   !>
   !> @param[in] self Becke cell
   !> @param[in] x    Size-adjusted elliptic coordinate
   pure function becke_partition_cell_eval(self, x) result(s)
      !> Becke cell
      class(becke_partition_cell), intent(in) :: self
      !> Pair coordinate
      real(wp), intent(in) :: x
      !> Cell weight
      real(wp) :: s
      s = becke_cell(x, self%k)
   end function becke_partition_cell_eval

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
