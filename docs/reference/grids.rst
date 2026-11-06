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
* ``3d/cartesian.f90``: uniform box grid with FFT transforms, configured by box dimensions and spacing.
* ``3d/molecular.f90``: translates atomic grids to the nuclei, applies partition weights, and prunes points; NUFFT transforms connect these nonuniform points to a uniform reciprocal grid.
* ``3d/kernel/``: interchangeable Becke, SSF, and fuzzy power-Voronoi (``pvoronoi``) partitions and grid derivatives. Atomic recipes and partitions can be chosen independently. SSF (:math:`C^3`) and pvoronoi (:math:`C^\infty`) provide regions of exactly zero weight for pruning.

``grid%new_trafo`` creates the transform engine matching the concrete 3D grid.
Rules, mappings, angular generators, and shell policies are independently replaceable through their abstract interfaces.

To prune only zero-weight molecular points, set both constructor arguments
``weight_threshold`` and ``pruning_threshold`` to zero. In the RISM input,
these are ``grid.molecular.weight_pruning_cutoff`` and
``grid.molecular.partition_pruning_cutoff``. Larger positive cutoffs
(:math:`>10^{-10}`) can introduce discontinuities in grid-evaluated/integrated functions.

Construction and updates
~~~~~~~~~~~~~~~~~~~~~~~~

The public 3D constructors, re-exported by ``moist_math_grid``, are:

* ``new_cartesian_point_grid(grid, error, nx, ny, nz, dr, margin)``
* ``new_cartesian_gaussian_grid(grid, error, nx, ny, nz, dr, xi0_factor, margin)``
* ``new_molecular_point_grid(grid, error, ...)``
* ``new_molecular_gaussian_grid(grid, error, ..., xi0_factor=...)``

Constructors configure the grid; call ``grid%update(mol, error)`` before using its points, weights, or transforms.
Changed settings require a new ``update`` call and ``destroy`` releases geometry while retaining the configuration.

``grid%w(i)`` is the final volume integration weight at point ``i``, in bohr cubed.
Integration, Fourier transforms, and model responses use this same array directly.

Cartesian grids default to 128 points per axis and a spacing of 0.6 bohr, with weights ``dr**3``.
The periodic box follows the nuclear centroid; every nucleus must remain at least ``margin`` (default 10 bohr) from its faces.
Molecular grids accept atomic recipes, partition settings, pruning thresholds, and reciprocal-grid settings.
``raw_index(i)`` is the position of retained molecular point ``i`` among the unpruned points (atom, radial shell, angular point); for a fixed solute and atomic recipe it names the same atom-relative point at every structure, while pruning may retain other points.
Choose ``kbuffer`` large enough for the expected motion: updates fail if the geometry outgrows the fixed reciprocal period.

Gaussian widths are ``xi0_factor/dr`` for Cartesian grids and ``xi0_factor/w**(1/3)`` for molecular grids, which require positive retained weights.
``xi0_factor`` defaults to 1 and is accepted only by Gaussian constructors.

Transforms
~~~~~~~~~~

``fft_r2k(f_r, f_k, error)`` and ``fft_k2r(f_k, f_r, error)`` transform blocks
``f_r(ngrid, nv)`` and ``f_k(npts_k, nv)`` with ``nv >= 1`` columns in both.
Molecular engines require ``prepare(nv)`` before use, including after assignment.
Copying a prepared engine with ``allocate(source=)`` or polymorphic assignment is unsupported.

``molecular_grid_set_kgrid`` sets the reciprocal grid explicitly; omitted ``buffer`` and ``nufft_tol`` retain their current values.
Both molecular constructors accept ``nufft_upsampfac``: ``0.0`` (default) lets
FINUFFT choose its internal FFT upsampling factor; a finite value greater than
one fixes it for both transform directions and both NUFFT routes.
``grid%get_nufft_upsampfac()`` returns the requested setting.

Nuclear derivatives
~~~~~~~~~~~~~~~~~~~

``volume_adjoint_type`` collects energy derivatives with respect to point positions, integration weights, and Gaussian widths.
Pass these adjoints as ``acc`` to:

* ``grid%get_volume_gradient(acc, gradient, error)``: adds the nuclear gradient to ``gradient(3, natom)``.
* ``grid%get_volume_hessian_vector(acc, direction, hvp, error)``: adds the nuclear Hessian times ``direction(3, natom)`` to ``hvp(3, natom)``.
* ``grid%get_volume_hessian(acc, hessian, error)``: adds the full nuclear Hessian to ``hessian(3, natom, 3, natom)``.

These calls require an updated grid, leave adjoints unchanged, and preserve the output on error.
Derivatives hold the retained point set, atomic recipes, and partition parameters fixed.
Point grids reject nonzero width adjoints.

Hessian methods provide grid curvature at fixed model adjoints; a complete energy Hessian also requires the model's response to grid observables.
Molecular Hessian calls require ``acc%w_xi == 0`` because Gaussian-width curvature is excluded.
Hessian calls reject Becke stiffness-one samples exactly on nuclei and coincident elliptic centers, where the response is undefined.
Power-coordinate partitions support coincident centers.

Layouts and transfers
---------------------

A layout describes what a field on a grid needs to be stored, compared or moved onto another grid; checkpoints use it to restart on other grids (:ref:`moz-checkpoints`).

* ``layout.f90``: ``moist_math_grid_layout_type`` (``moist_math_grid_layout``) with the kind (Cartesian box, radial nodes or points), shape, spacing, origin, points, weights and optional Fourier basis of a grid. ``grid%get_layout`` returns it for 3D and radial grids; ``validate`` checks it.
* ``transfer.f90``: ``grid_transfer_type`` (``moist_math_grid_transfer``) moves fields from one layout onto another by ``copy``, ``fourier``, ``trilinear``, ``points`` or ``spline``; ``fit_rigid`` gives the rotation and translation that map one solute onto another.
* ``transfer_points.f90``: the moving least squares stencils of the ``points`` transfer.
* ``deposit.f90``: ``grid_deposit_type`` (``moist_math_grid_deposit``) moves fields between scattered points and a periodic cubic lattice of cells, by trilinear interpolation from the cell centres and by deposition weighted with the quadrature weights; the ``cartesian-aux`` preconditioner uses it (:doc:`/models/moz/uv3d`).
