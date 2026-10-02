!> Stratmann-Scuseria-Frisch molecular partition with exact zero regions
module moist_math_grid_3d_partition_ssf
   use mctc_env, only: wp
   use moist_math_grid_3d_partition_common, only: pair_partition_weights, partition_cell
   implicit none(type, external)
   private

   public :: ssf_partition_weights, ssf_cell, default_ssf_a

   !> Recommended SSF half-width in elliptic coordinates
   real(wp), parameter :: default_ssf_a = 0.64_wp

   !> SSF cell with a fixed switching width
   type, extends(partition_cell) :: ssf_partition_cell
      !> Switching half-width
      real(wp) :: width = default_ssf_a
   contains
      !> Pair-cell weight
      procedure :: eval => ssf_partition_cell_eval
   end type ssf_partition_cell

contains

   !> Normalized size-adjusted SSF weights of one atom at a batch of points
   !>
   !> Unchecked shapes and valid atomic numbers/owner required; no quadrature factors
   !>
   !> @param[in] owner Owner index in [1,nat]
   !> @param[in] points Points in bohr, shape (3,npts)
   !> @param[in] xyz Atomic positions in bohr, shape (3,nat)
   !> @param[in] numbers Atomic numbers, shape (nat)
   !> @param[out] w Owner weights in [0,1], shape (npts)
   !> @param[in] a Switching half-width in (0,1]; default 0.64
   subroutine ssf_partition_weights(owner, points, xyz, numbers, w, a)
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
      !> Switching half-width
      real(wp), intent(in), optional :: a

      type(ssf_partition_cell) :: cell

      if (present(a)) cell%width = a
      call pair_partition_weights(owner, points, xyz, numbers, w, cell)
   end subroutine ssf_partition_weights

   !> SSF cell with the configured switching width
   !>
   !> @param[in] self SSF cell
   !> @param[in] x    Size-adjusted elliptic coordinate
   pure function ssf_partition_cell_eval(self, x) result(s)
      !> SSF cell
      class(ssf_partition_cell), intent(in) :: self
      !> Pair coordinate
      real(wp), intent(in) :: x
      !> Cell weight
      real(wp) :: s
      s = ssf_cell(x/self%width)
   end function ssf_partition_cell_eval

   !> C^3 SSF switch, evaluated without subtractive cancellation at the zero end
   !>
   !> doi:10.1016/0009-2614(96)00600-8, Eq. 11
   !> equiv. to [1-(35*z-35*z**3+21*z**5-5*z**7)/16]/2 inside (-1,1)
   !>
   !> @param[in] z Elliptic coordinate divided by the switching half-width
   pure function ssf_cell(z) result(s)
      !> Scaled pair coordinate
      real(wp), intent(in) :: z
      !> Cell weight in [0,1], exactly zero for z >= 1
      real(wp) :: s

      real(wp) :: t

      t = 0.5_wp*(1.0_wp - min(1.0_wp, abs(z)))
      s = t**4*(35.0_wp + t*(-84.0_wp + t*(70.0_wp - 20.0_wp*t)))
      if (z < 0.0_wp) s = 1.0_wp - s
   end function ssf_cell

end module moist_math_grid_3d_partition_ssf
