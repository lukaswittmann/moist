iSwiG cavity
============

iSwiG is the switching Gaussian surface-discretization approach :cite:p:`lange2010swig` for van der Waals type cavitites.
Each atomic sphere is discretized using Lebedev quadrature grids.
Every grid point is assigned a smooth weight given by a product of error-function switches from the neighbouring spheres, so buried points fade out continuously.
Points whose switching value (or area) falls below a cutoff are removed.

The implementation supports analytic nuclear derivatives using adaptive radii :cite:p:`wittmann2025cpcm`.

Settings
--------

All three settings are optional constructor arguments.

``nleb`` (integer, default ``194``)
   Number of Lebedev quadrature points per atomic sphere; the primary accuracy/cost setting.
   **Must be one of the supported sizes**: 14, 26, 50, 110, 194, 302, 434, 590, 770, 974, 1202.
   Any other value raises an error, because each grid has a fitted Born ``zeta`` (the Gaussian width scale).

``cut_f`` (real, default ``1.0e-10``)
   Switching-value cutoff. A point is kept only if its switching function exceeds this value. Used when ``cut_a`` is not set.

``cut_a`` (real, default ``0.0``)
   Area cutoff. When greater than zero it replaces ``cut_f``: a point is kept only if its switched area (Lebedev area times switching value) exceeds this value.
