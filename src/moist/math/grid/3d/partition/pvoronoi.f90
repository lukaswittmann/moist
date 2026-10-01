!> Custom C-infinity fuzzy power-Voronoi molecular partition
module moist_math_grid_3d_partition_pvoronoi
   use mctc_env, only: wp
   use moist_math_grid_3d_partition_common, only: pair_partition_weights
   implicit none(type, external)
   private

   public :: pvoronoi_partition_weights, bump_cell, default_power_width

   !> Default power-gap switching half-width in bohr**2
   real(wp), parameter :: default_power_width = 1.0_wp

contains

   !> Normalized bump-function weights of fuzzy power cells
   !>
   !> q_i=|r-R_i|**2-radius_i**2; P_i=product_j bump_cell((q_i-q_j)/width)
   !>
   !> - A minimum-power site has every factor >= 1/2: the denominator is positive
   !> - C-infinity in points and nuclei for fixed positive width and fixed radii
   !> - Exact zeros outside expanded power cells; outer cells remain unbounded
   !> - Finite-width transitions, not a finite union of compact atomic supports
   !> - No nearest-site selection, distance norms or hard radial cutoff
   !> - Unchecked shapes and valid atomic numbers/owner required
   !>
   !> @param[in] owner Owner index in [1,nat]
   !> @param[in] points Points in bohr, shape (3,npts)
   !> @param[in] xyz Atomic positions in bohr, shape (3,nat)
   !> @param[in] numbers Atomic numbers, shape (nat)
   !> @param[out] w Owner weights in [0,1], shape (npts)
   !> @param[in] width Positive power-gap half-width in bohr**2; default 1
   !> @param[in] radii Fixed nonnegative power radii in bohr, shape (nat); default covalent
   subroutine pvoronoi_partition_weights(owner, points, xyz, numbers, w, width, radii)
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
      !> Power-gap half-width
      real(wp), intent(in), optional :: width
      !> Fixed power radii
      real(wp), intent(in), optional :: radii(:)

      real(wp) :: half_width

      half_width = default_power_width
      if (present(width)) half_width = width
      call pair_partition_weights(owner, points, xyz, numbers, w, bump_cell, &
         & power_width=half_width, radii=radii)
   end subroutine pvoronoi_partition_weights

   !> C-infinity switch b(1-x)/(b(1-x)+b(1+x)), b(t)=exp(-1/t) for t>0, else 0
   !>
   !> Stable logistic evaluation; all derivatives vanish at x=+/-1
   !>
   !> @param[in] x Power gap divided by its switching half-width
   pure function bump_cell(x) result(s)
      !> Scaled power gap
      real(wp), intent(in) :: x
      !> Pair-cell weight in [0,1]
      real(wp) :: s

      real(wp) :: t, e

      t = abs(x)
      s = 0.0_wp
      if (t < 1.0_wp) then
         e = exp(-2.0_wp*t/((1.0_wp - t)*(1.0_wp + t)))
         s = e/(1.0_wp + e)
      end if
      if (x < 0.0_wp) s = 1.0_wp - s
   end function bump_cell

end module moist_math_grid_3d_partition_pvoronoi
