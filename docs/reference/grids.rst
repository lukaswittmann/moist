Integration grids
=================

MOIST provides radial, angular, and three-dimensional integration grids.
Atomic grids combine radial and angular quadrature; molecular grids combine atomic grids using some partition scheme.
Recipes separate discretization choices from the resulting points and weights.
``moist_math_grid`` re-exports the main API; paths below are relative to ``src/moist/math/grid/``.

Radial
------

Radial grids provide quadrature along the radial coordinate.
The quadrature rule defines nodes and weights on a reference interval; the mapping sets the radial range and point distribution.

* ``radial/rule.f90``: abstract quadrature rule on ``[-1, 1]``; midpoint, Chebyshev (second kind), and Gauss-Legendre implementations.
* ``radial/mapping.f90``: abstract mapping to radii; linear and HandyMod for finite intervals, Becke and Knowles for ``[0, infinity)``.
* ``radial/grid.f90``: a recipe combines a rule, mapping, node count, and optional radial cutoffs. The grid stores radii and ``dr`` weights; its volume integrators include :math:`4\pi r^2` for spherically symmetric functions.
* ``radial/trafo.f90``: radial Fourier-Bessel transforms. Pairs built by ``new_uniform_radial_pair`` use DST-IV; pairs of recipe-built grids use general quadrature.

Angular
-------

Angular grids resolve directions on the unit sphere, independently of radial discretization.
Resolution requests select an available angular rule.

* ``angular/lebedev.f90``: Lebedev-Laikov unit-sphere tables.
* ``angular/grid.f90``: directions and solid-angle weights summing to :math:`4\pi`. An abstract generator selects a grid by requested point count or exactness degree; Lebedev is the supplied implementation.

Atomic: one center
------------------

An atomic grid combines radial and angular quadrature around a single center.
Angular resolution can vary between radial shells.

* ``atomic/recipe.f90``: combines a radial recipe, angular generator, and shell policy. Policies request constant angular degree, degree by radial sector, or target arc spacing; each works with any compatible generator.
  Recipes can be overridden per element.
* ``atomic/grid.f90``: builds radial shells with their selected angular grids, producing center-relative points and :math:`r^2\,dr\,d\Omega` weights.

3D: volume grids
----------------

Cartesian and molecular grids share a 3D interface, allowing field calculations to use either representation.
Molecular partition weights sum to one over atoms at each point, preventing double counting of overlapping atomic contributions.

* ``3d/base.f90``: abstract ``moist_math_grid_3d_type`` with common point, weight, integration, geometry-update, and transform-factory interfaces. ``3d/adjoint.f90`` holds position, weight, and Gaussian-width adjoints.
* ``3d/cartesian.f90``: uniform box grid with FFT transforms; constructed directly from box dimensions and spacing.
* ``3d/molecular.f90``: translates atomic grids to the nuclei, applies partition weights, and prunes points; NUFFT transforms connect these nonuniform points to a uniform reciprocal grid.
* ``3d/partition.f90``: interchangeable Becke, SSF, and fuzzy power-Voronoi (``pvoronoi``) partitions, implemented in ``3d/partition/``. Atomic recipes and partitions can be chosen independently. SSF (:math:`C^3`) and pvoronoi (:math:`C^\infty`) provide regions of exactly zero weight for pruning.

``grid%new_trafo`` creates the transform engine matching the concrete 3D grid.
Rules, mappings, angular generators, and shell policies are independently replaceable through their abstract interfaces.

To prune only zero-weight molecular points, set both ``weight_threshold`` and ``pruning_threshold`` to zero; larger positive thresholds (:math:`>10^{-10}`) can introduce discontinuities in grid-evaluted/integrated functions.