DROP Cavities
=============

.. image:: /_static/drop.svg
   :alt: DROP cavity construction
   :align: right

The Discretization via Reference-Onto-Surface Projection (DROP) :cite:p:`wittmann2026drop` forms the foundation for cavities that are based on implicit surfaces.
DROP is a general scheme for discretizing implicit surfaces.
The individual steps are illustrated in :ref:`the DROP scheme below <fig-drop-scheme>`.
It starts from a reference system (van der Waals type cavity) which is subsequently discretized using Lebedev quadrature grids (a).
The reference cavity is then projected onto the defined zero level set given by the level set function (LSF) (b).
After projection, the quadrature weights have to be mapped in order to obtain correct surface integration weights (c).
The resulting surface grid inherits the smoothness of the underlying LSF and the projection/mapping scheme and thus  is provides consistent analytical derivatives.

DROP is independent of the particular surface definition; currently available:
 - :doc:`Electron-isodensity surface (ρ-DROP) <isodensity>`, a fully self-consistent and differentiable cavity based on the electron density and a chosen isodensity value.
 - :doc:`Smooth van der Waals surface (SvdW-DROP) <svdw>`, a fully differentiable molecular surface that recovers solvent-excluded-surface-like features while avoiding geometric singularities, crevices, and other discontinuous surface features.
 - :doc:`COSMO Fine Cavity (CFC-DROP) <cfc>`.

During the projection and weighting procedure, DROP applies smooth switching functions, handles multiple projection branches in concave regions, and uses spatial screening to keep the construction differentiable and efficient.
This makes the cavity suitable for continuum-solvation models that require smooth surface areas, polarization energies, and gradients, such as CPCM and pressure-based cavity terms.

.. _fig-drop-scheme:

.. grid:: 1
   :gutter: 2

   .. grid-item::

      .. image:: /_static/drop_scheme_a.svg
         :alt: Isodensity surface and van der Waals reference cavity
         :align: center
         :width: 65%

      **(a)** The isodensity surface is shown in blue and the hard-sphere van
      der Waals reference cavity in gray.

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item::

      .. image:: /_static/drop_scheme_b.svg
         :alt: Projection of reference quadrature points onto the isodensity surface
         :align: center
         :width: 80%

      **(b)** Closest-point projection of reference quadrature points onto the
      isodensity surface, :math:`S=0`.

   .. grid-item::

      .. image:: /_static/drop_scheme_c.svg
         :alt: Quadrature weight mapping onto the isodensity surface
         :align: center
         :width: 80%

      **(c)** Quadrature weight mapping from the reference cavity,
      :math:`\mathrm{d}A^\circ`, to the isodensity surface,
      :math:`\mathrm{d}A^\ast`.

**Figure:** Illustration of the DROP scheme. A van der Waals reference cavity
(gray, dashed) is discretized (a), its quadrature points are projected (b) onto
the isodensity surface (blue, solid), and the weights are transformed using the
corresponding surface Jacobian (c).

Settings
--------

Settings can be supplied two ways:

- **Parameter file**: as dotted keys in a configuration file, loaded via ``load_file`` (e.g. ``grid.num_leb``, ``projection.level``). The full key set is registered in ``register_cavity_drop_entries``.
- **Constructor arguments**: a subset is exposed directly as optional arguments to the cavity constructor: ``nleb``, ``tolerance``, ``proj_maxiter``, ``proj_level``, ``branch_weight_s``, ``rho_grid_h``,
  ``wleb_prune_level``.


Level Set Functions (LSF)
-------------------------

The LSF defines the implicit surface that the reference (van der Waals) grid is
projected onto. Unlike the dotted keys documented below, the LSF *and its shape
parameters* are construction-level choices: the caller builds a concrete LSF and
passes it to the cavity constructor (the ``lsf_model`` argument in the Fortran
API). They are not parameter-file keys and cannot be set through ``load_file``.

MOIST provides SvdW, CFC, and isodensity LSFs:

.. toctree::
   :maxdepth: 1

   isodensity
   svdw
   cfc

Grid Discretization
-------------------

``grid.num_leb`` (integer, default ``194``)
   Number of Lebedev quadrature points per atomic sphere -- the primary  accuracy/cost knob.
   **Must be one of the supported Lebedev orders**: 6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810.


Numerical Tolerance
-------------------

``tolerance`` (real, default ``1.0e-10``)
   Main tolerance controlling the full threshold hierarchy.
   Everything else scales from it (tightest to loosest): 
   - weight cutoff ``= tolerance * 0.05``,
   - LSF screening ``= tolerance * 0.1``, projection convergence ``= tolerance``,
   - branch-degeneracy separation ``= tolerance * 10``.
   See :cite:t:`wittmann2026drop`, Supporting Information Sec. C.2.c for convergence behavior with respect to these thresholds.


Surface Projection
------------------

``projection.level`` (integer, default ``3``)
   Strategy used to project grid points onto the surface:

   .. list-table::
      :header-rows: 1
      :widths: 10 40 50

      * - Level
        - Strategy
        - Recomendation
      * - 1
        - SLSQP :cite:p:`kraft1988slsqp,kraft1994slsqp` only (no Newton refinement)
        - Not recommended
      * - 2
        - SLSQP + Newton refinement
        - Fine for most cases
      * - 3
        - conditional multi-tangent for degenerate points
        - Recommended for general use
      * - 4
        - conditional SLSQP-deflation
        - Recommended for more challenging cases
      * - 5
        - SLSQP-deflation (unconditional)
        - Not recommended
      * - 6
        - Newton-deflation on the 4D KKT system
        - Not recommended
      * - 7
        - regular SLSQP multistart
        - Recommended for more challenging cases
      * - 8
        - fine SLSQP multistart reference
        - Reference-quality, not recommended for routine use

``projection.maxiter`` (integer, default ``150``)
  Maximum iterations for the projection optimizer.

``objective.alpha`` (real, default ``0.5``)
  Weight of the anchor term in the projection objective; larger values keep projected points closer to their initial (anchor) position. Adjusting this should have no effect on the final surface.


Switching Functions
-------------------

Smooth step functions that fade surface contributions in and out so the cavity stays differentiable also in rare edge cases.
Two independent switches act on the integration weights, each keyed on a different geometric quantity (see the note below for the underlying rationale):

**Critical level set weight switch** (``f_crit``); a functin of the level set gradient norm ``||\\nabla S||``:

``switching.w_0ls_from`` / ``switching.w_0ls_to`` (real, defaults ``0.25`` / ``0.6``)
   Transition bounds of the critical level set weight switch (start, end).
   Contributions are fully suppressed below ``w_0ls_from`` and fully restored above ``w_0ls_to``.

**Focal/branching weight switch** (``f_foc``); a function of the smallest eigenvalue of the Lagrangian Hessian restricted to the tangent space (TRLH):

``switching.w_0tra_from`` / ``switching.w_0tra_to`` (real, defaults ``0.1`` / ``0.3``)
   Transition bounds of the focal/branching weight switch (start, end).
   A point's contribution is damped as its tangential curvature drops from ``w_0tra_to`` toward ``w_0tra_from``.

.. note::

  Although the level set function itself is differentiable away from nuclei, the resulting zero level set can become singular where ``S=0`` and ``\\nabla S=0``.
  To avoid ill-conditioned or non-unique projections in such critical regions, a critical-point switching function, ``f_crit(||\\nabla S||)``, smoothly attenuates contributions from the zero level set.
  It is inactive for regular surface regions and suppresses contributions (``f_crit=0``) when ``||\\nabla S|| < 0.25`` in the present work, which corresponds to scaling the surface measure and the discretized integration weights.

  A second switch handles focal events, where the closest-point map loses regularity in a tangential direction. In practice this is detected from the smallest eigenvalue of the Lagrangian Hessian restricted to the tangent
  space.
  When that tangential curvature becomes too small, the focal switch ``f_foc`` smoothly damps the contribution of the affected point so the discretized surface remains stable.


Weight Pruning
--------------

If desired, near-zero weights can be smoothly attenuated to speed up the computation of electrostatic potential integrals (for QM coupling) in cases with large numbers of points with negligible contributions (e.g. large cavities with high Lebedev order).
Additionally, it can be used in cases where small weights could cause numerical issues.

``switching.wleb_prune_level`` (integer, default ``0``)
   Smoothly suppresses near-zero Lebedev weights before the branch filter. The
   switch region is derived from the level:

   .. list-table::
      :header-rows: 1
      :widths: 20 80

      * - Level
        - Switch region
      * - 0
        - disabled (default)
      * - 1
        - 1e-12 to 1e-10
      * - 2
        - 1e-10 to 1e-8
      * - 3
        - 1e-8 to 1e-6
      * - 4
        - 1e-6 to 1e-4


Branching
---------

``branching.softmax_scale`` (real, default ``0.0025``)
  Softmax scale of the branch-weight model in concave regions; smaller values   make the branch-weight distribution sharper.
  Branch fractions are computed from the objective values of all competing projections with a softmax over ``-\\Phi`` divided by this scale, so equal objective values receive equal weights and a lower objective value smoothly dominates.
  The branch ``rho`` cutoff is derived from this and the weight cutoff, and the resulting branch fraction is  multiplied into the integration weight for each projected point.


Spatial Screening
-----------------

Acceleration settings that do not change results, only cost.

``screening.cell_grid_fraction`` (real, default ``0.25``)
   Cell size of the molecular cell grid as a fraction of the reference spacing (``1.0`` = no subdivision, ``0.5`` = halved cells).
   This has no influence on the actual results but speeds up the LSF evaluations.

``screening.cell_grid_full_scan_below`` (integer, default ``200``)
   Below this atom count the cell grid collapses to a single full-scan cell.


Disconnected Points
-------------------

``disconnection.threshold`` (real, default ``4.0``)
   Distance threshold (in units of the average grid spacing) above which a projected point is treated as disconnected from the surface.
