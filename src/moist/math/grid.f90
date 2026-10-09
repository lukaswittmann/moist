!> Quadrature grids: radial, angular, atomic, and three-dimensional
!>
!> The layers build on each other (from the bottom up):
!>
!> Radial:
!> - `radial/rule.f90`: 1D quadrature rules on [-1, 1]: midpoint, Chebyshev (2nd kind), Gauss-Legendre
!> - `radial/mapping.f90`: maps of [-1, 1] to radii: linear and HandyMod onto finite intervals, Becke and Knowles
!>   (Mura-Knowles) onto [0, inf)
!> - `radial/grid.f90`: radial grid, nodes and weights from a rule and a mapping; uniform r/k pairs
!> - `radial/trafo.f90`: r <-> k transforms of radial functions: DST-IV on uniform pairs, quadrature otherwise
!>
!> Angular:
!> - `angular/lebedev.f90`: Lebedev tables, points and weights on the unit sphere
!> - `angular/grid.f90`: unit-sphere grid, angular requests, and Lebedev generator
!>
!> Atomic (3d grid for single atom):
!> - `atomic/grid.f90`: atomic grid, radial shells times angular grids around one nucleus
!> - `atomic/recipe.f90`: per-element build settings: radial recipe, angular generator, angular size per shell
!>
!> 3D (molecular grid for many atoms):
!> - `3d/base.f90`: abstract volume grid with field integration; adjoint channels in `3d/adjoint.f90`
!> - `3d/cartesian.f90`: uniform Cartesian box grid with FFT r <-> k transforms
!> - `3d/molecular.f90`: molecular grid, partitioned union of atomic grids, NUFFT r <-> k transforms
!> - `3d/kernel/`: generated Becke, SSF and power-Voronoi partition kernels and grid derivatives
!>
!> This module re-exports the main public names of every layer; the check routines and the table
!> metadata stay in their modules, import those directly when needed
module moist_math_grid
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_type, &
      & moist_math_grid_radial_rule_chebyshev2_type, new_chebyshev2_rule, &
      & moist_math_grid_radial_rule_midpoint_type, new_midpoint_rule, &
      & moist_math_grid_radial_rule_gauss_legendre_type, new_gauss_legendre_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_type, &
      & moist_math_grid_radial_mapping_linear_type, new_linear_mapping, &
      & moist_math_grid_radial_mapping_becke_type, new_becke_mapping, &
      & moist_math_grid_radial_mapping_handymod_type, new_handymod_mapping, &
      & moist_math_grid_radial_mapping_knowles_type, new_knowles_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, &
      & moist_math_grid_radial_recipe_type, new_radial_grid, new_uniform_radial_pair, &
      & integrand_radial, transform_quadrature, transform_dst4
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, &
      & moist_math_grid_radial_trafo_dst4_type, moist_math_grid_radial_trafo_quadrature_type, &
      & new_radial_trafo
   use moist_math_grid_angular_lebedev, only: &
      & grid_size, get_angular_grid, lebedev_order_from_num, lebedev_degree_table
   use moist_math_grid_angular_grid, only: moist_math_grid_angular_type, integrand_angular, &
      & moist_math_grid_angular_generator_type, moist_math_grid_angular_request_type, &
      & moist_math_grid_angular_generator_lebedev_type, new_lebedev_generator, new_lebedev_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_shell_type, &
      & moist_math_grid_atomic_shell_constant_type, new_constant_shell_policy, &
      & moist_math_grid_atomic_shell_sector_type, new_sector_shell_policy, &
      & moist_math_grid_atomic_shell_arc_type, new_arc_shell_policy, &
      & moist_math_grid_atomic_recipe_type, moist_math_grid_atomic_recipe_override_type, &
      & default_molecular_recipe, default_element_recipes, get_element_recipe, &
      & element_override_index
   use moist_math_grid_atomic_grid, only: moist_math_grid_atomic_type, new_atomic_grid, &
      & integrand_atomic
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type, moist_math_grid_3d_trafo_type, &
      & integrand_3d
   use moist_math_grid_3d_kernel_base, only: moist_math_grid_3d_partition_type
   use moist_math_grid_3d_kernel_becke, only: becke_partition_type
   use moist_math_grid_3d_kernel_ssf, only: ssf_partition_type
   use moist_math_grid_3d_kernel_pvoronoi, only: pvoronoi_partition_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, &
      & new_cartesian_point_grid, new_cartesian_gaussian_grid, moist_math_grid_3d_cartesian_trafo_type
   use moist_math_grid_3d_molecular, only: &
      & moist_math_grid_3d_molecular_type, new_molecular_point_grid, new_molecular_gaussian_grid, &
      & molecular_grid_set_kgrid, moist_math_grid_3d_molecular_trafo_type, &
      & new_molecular_grid_trafo, default_nufft_tol, partition_becke, partition_ssf, partition_pvoronoi
   implicit none(type, external)
   private

   public :: grid_size
   public :: get_angular_grid
   public :: lebedev_order_from_num
   public :: lebedev_degree_table
   public :: moist_math_grid_radial_rule_type
   public :: moist_math_grid_radial_rule_chebyshev2_type
   public :: new_chebyshev2_rule
   public :: moist_math_grid_radial_rule_midpoint_type
   public :: new_midpoint_rule
   public :: moist_math_grid_radial_rule_gauss_legendre_type
   public :: new_gauss_legendre_rule
   public :: moist_math_grid_radial_type
   public :: moist_math_grid_radial_recipe_type
   public :: new_radial_grid
   public :: new_uniform_radial_pair
   public :: integrand_radial
   !> Radial transform tag of the general quadrature
   public :: transform_quadrature
   !> Radial transform tag of the uniform DST-IV pair
   public :: transform_dst4
   public :: moist_math_grid_radial_mapping_type
   public :: moist_math_grid_radial_mapping_linear_type
   public :: new_linear_mapping
   public :: moist_math_grid_radial_mapping_becke_type
   public :: new_becke_mapping
   public :: moist_math_grid_radial_mapping_handymod_type
   public :: new_handymod_mapping
   public :: moist_math_grid_radial_mapping_knowles_type
   public :: new_knowles_mapping
   public :: moist_math_grid_radial_trafo_type
   public :: moist_math_grid_radial_trafo_dst4_type
   public :: moist_math_grid_radial_trafo_quadrature_type
   public :: new_radial_trafo
   public :: moist_math_grid_angular_type
   public :: integrand_angular
   public :: moist_math_grid_angular_generator_type
   public :: moist_math_grid_angular_request_type
   public :: moist_math_grid_angular_generator_lebedev_type
   public :: new_lebedev_generator
   public :: new_lebedev_grid
   public :: moist_math_grid_atomic_shell_type
   public :: moist_math_grid_atomic_shell_constant_type
   public :: new_constant_shell_policy
   public :: moist_math_grid_atomic_shell_sector_type
   public :: new_sector_shell_policy
   public :: moist_math_grid_atomic_shell_arc_type
   public :: new_arc_shell_policy
   public :: moist_math_grid_atomic_recipe_type
   public :: moist_math_grid_atomic_recipe_override_type
   public :: default_molecular_recipe
   public :: default_element_recipes
   public :: get_element_recipe
   public :: element_override_index
   public :: moist_math_grid_atomic_type
   public :: new_atomic_grid
   public :: integrand_atomic
   public :: moist_math_grid_3d_type
   public :: moist_math_grid_3d_trafo_type
   public :: volume_adjoint_type
   public :: integrand_3d
   public :: moist_math_grid_3d_partition_type
   public :: becke_partition_type, ssf_partition_type, pvoronoi_partition_type
   public :: partition_becke, partition_ssf, partition_pvoronoi
   public :: moist_math_grid_3d_cartesian_type
   public :: new_cartesian_point_grid, new_cartesian_gaussian_grid
   public :: moist_math_grid_3d_cartesian_trafo_type
   public :: moist_math_grid_3d_molecular_type
   public :: new_molecular_point_grid, new_molecular_gaussian_grid
   public :: molecular_grid_set_kgrid
   public :: moist_math_grid_3d_molecular_trafo_type
   public :: new_molecular_grid_trafo
   !> Default requested FINUFFT relative tolerance of the molecular transforms
   public :: default_nufft_tol

end module moist_math_grid
