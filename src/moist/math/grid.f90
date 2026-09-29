!> Umbrella module for the moist integration-grid submodule
!>
!> Re-exports the public API of the abstract grid bases and the concrete grids
!> moist builds on: Lebedev unit-sphere grids (`grid/s2`), radial grids with
!> their Fourier-Bessel transforms (`grid/1d`), and the volumetric Cartesian and
!> atom-centered molecular grids with their Fourier transforms (`grid/3d`)
!> The bare quadrature rules (`moist_math_quadrature_*`: Lebedev tables,
!> Chebyshev-2 and HandyMod radial rules, Becke partitioning) are kept off this
!> surface -- import their modules directly when needed
module moist_math_grid
   use moist_math_quadrature_lebedev, only: &
      & grid_size, get_angular_grid, lebedev_order_from_num, lebedev_locality_order
   use moist_math_grid_s2_base, only: moist_math_grid_s2_type, integrand_s2
   use moist_math_grid_s2_lebedev, only: moist_math_grid_s2_lebedev_type, new_s2_grid_lebedev
   use moist_math_grid_1d_base, only: moist_math_grid_1d_type, moist_math_grid_1d_trafo_type, &
      & integrand_radial
   use moist_math_grid_1d_uniform, only: moist_math_grid_1d_uniform_type, new_uniform_radial_grid, &
      & moist_math_grid_1d_uniform_trafo_type
   use moist_math_grid_1d_chebyshev, only: moist_math_grid_1d_chebyshev_type, &
      & new_chebyshev_radial_grid, moist_math_grid_1d_chebyshev_trafo_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type, &
      & integrand_3d
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_grid_3d, moist_math_grid_3d_cartesian_trafo_type
   use moist_math_grid_3d_molecular, only: &
      & moist_math_grid_3d_molecular_type, new_molecular_grid, &
      & new_molecular_grid_uniform, new_molecular_grid_uniform_handymod, &
      & new_molecular_grid_uniform_qc_handymod, &
      & default_grid_sizes, &
      & molecular_grid_set_kgrid, moist_math_grid_3d_molecular_trafo_type, &
      & new_molecular_grid_trafo, default_nufft_tol
   implicit none(type, external)
   private

   public :: grid_size
   public :: get_angular_grid
   public :: lebedev_order_from_num
   public :: lebedev_locality_order
   public :: moist_math_grid_s2_type
   public :: integrand_s2
   public :: moist_math_grid_s2_lebedev_type
   public :: new_s2_grid_lebedev
   public :: moist_math_grid_1d_type
   public :: moist_math_grid_1d_trafo_type
   public :: integrand_radial
   public :: moist_math_grid_1d_uniform_type
   public :: new_uniform_radial_grid
   public :: moist_math_grid_1d_uniform_trafo_type
   public :: moist_math_grid_1d_chebyshev_type
   public :: new_chebyshev_radial_grid
   public :: moist_math_grid_1d_chebyshev_trafo_type
   public :: moist_math_grid_3d_type
   public :: moist_math_grid_3d_trafo_type
   public :: integrand_3d
   public :: moist_math_grid_3d_cartesian_type
   public :: new_cartesian_grid_3d
   public :: moist_math_grid_3d_cartesian_trafo_type
   public :: moist_math_grid_3d_molecular_type
   public :: new_molecular_grid
   public :: new_molecular_grid_uniform
   public :: new_molecular_grid_uniform_handymod
   public :: new_molecular_grid_uniform_qc_handymod
   public :: default_grid_sizes
   public :: molecular_grid_set_kgrid
   public :: moist_math_grid_3d_molecular_trafo_type
   public :: new_molecular_grid_trafo
   !> Default requested FINUFFT relative tolerance of the molecular transforms
   public :: default_nufft_tol

end module moist_math_grid
