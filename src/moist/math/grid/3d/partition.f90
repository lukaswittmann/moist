!> Molecular partition schemes and the legacy Becke/SSF entry point
module moist_math_grid_3d_partition
   use mctc_env, only: wp
   use moist_math_grid_3d_partition_becke, only: becke_weights => becke_partition_weights, default_stiffness
   use moist_math_grid_3d_partition_ssf, only: ssf_partition_weights, default_ssf_a
   use moist_math_grid_3d_partition_pvoronoi, only: pvoronoi_partition_weights, default_power_width
   implicit none(type, external)
   private

   public :: becke_partition_weights, ssf_partition_weights, pvoronoi_partition_weights
   public :: default_stiffness, default_ssf_a, default_power_width
   public :: partition_becke, partition_ssf, partition_pvoronoi

   !> Original size-adjusted Becke partition
   integer, parameter :: partition_becke = 1
   !> Stratmann-Scuseria-Frisch partition
   integer, parameter :: partition_ssf = 2
   !> Custom C-infinity power-Voronoi partition
   integer, parameter :: partition_pvoronoi = 3

contains

   !> Backward-compatible Becke weights; ssf_a selects the separate SSF implementation
   !>
   !> Unchecked shapes and valid atomic numbers/owner required; no quadrature factors
   !>
   !> @param[in] owner     Owner index in [1,nat]
   !> @param[in] points    Sample points in bohr, shape (3,npts)
   !> @param[in] xyz       Atomic positions in bohr, shape (3,nat)
   !> @param[in] numbers   Atomic numbers, shape (nat)
   !> @param[out] w        Owner weights in [0,1], shape (npts)
   !> @param[in] stiffness Becke iteration count >= 1; default 3
   !> @param[in] ssf_a SSF half-width in (0,1]; takes precedence over stiffness
   subroutine becke_partition_weights(owner, points, xyz, numbers, w, stiffness, ssf_a)
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
      !> Becke iteration count
      integer, intent(in), optional :: stiffness
      !> SSF half-width
      real(wp), intent(in), optional :: ssf_a

      if (present(ssf_a)) then
         call ssf_partition_weights(owner, points, xyz, numbers, w, a=ssf_a)
      else
         call becke_weights(owner, points, xyz, numbers, w, stiffness=stiffness)
      end if
   end subroutine becke_partition_weights

end module moist_math_grid_3d_partition
