!> Umbrella module for the moist integration-grid submodule.
!>
!> Re-exports the public API for Lebedev angular grids and atom-centered
!> molecular integration grids. Chebyshev-2 radial quadrature and Becke
!> fuzzy-cell partitioning are kept out of the public surface; import
!> their sub-modules (`moist_math_grid_radial`, `moist_math_grid_becke`)
!> directly if they are needed.
module moist_math_grid
   use moist_math_grid_lebedev, only: &
      & grid_size, get_angular_grid, lebedev_order_from_num
   use moist_math_grid_molecular, only: &
      & molecular_grid_type, new_molecular_grid, &
      & new_molecular_grid_uniform, default_grid_sizes, integrand_3d
   implicit none
   private

   public :: grid_size
   public :: get_angular_grid
   public :: lebedev_order_from_num
   public :: molecular_grid_type
   public :: new_molecular_grid
   public :: new_molecular_grid_uniform
   public :: default_grid_sizes
   public :: integrand_3d

end module moist_math_grid
