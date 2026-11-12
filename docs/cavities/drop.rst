DROP Cavities
=============

.. image:: /_static/logos/drop.svg
   :alt: DROP cavity construction
   :align: right

The Discretization via Reference-Onto-Surface Projection (DROP) :cite:p:`wittmann2026drop` forms the foundation for cavities that are based on implicit surfaces.
DROP is a general scheme for discretizing implicit surfaces.
The individual steps are illustrated in :ref:`the DROP scheme below <fig-drop-scheme>`.
It starts from a reference system (van der Waals type cavity) which is subsequently discretized using Lebedev quadrature grids (a).
The reference cavity is then projected onto the defined zero level set given by the level set function (LSF) (b).
After projection, the quadrature weights have to be mapped in order to obtain correct surface integration weights (c).
The resulting surface grid inherits the smoothness of the underlying LSF and the projection/mapping scheme and thus  is provides consistent analytical derivatives.

DROP is independent of the particular surface definition; the two main level sets are:
 - :doc:`Electron-isodensity surface (ρ-DROP) <isodensity>`, a fully self-consistent and differentiable cavity based on the electron density and a chosen isodensity value.
 - :doc:`Smooth van der Waals surface (SvdW-DROP) <svdw>`, a fully differentiable molecular surface that recovers solvent-excluded-surface-like features while avoiding geometric singularities, crevices, and other discontinuous surface features.

.. subfigure:: ABC
   :layout-sm: A|B|C
   :subcaptions: below
   :gap: 8px
   :name: fig-drop-lsf-overview

   .. image:: /_static/figures/lsf/vdw.svg
      :alt: (a) Van der Waals sphere union

   .. image:: /_static/figures/lsf/svdw_no_circles_no_lsf_lines.svg
      :alt: (b) SvdW-DROP

   .. image:: /_static/figures/lsf/isodensity_no_circles_no_lsf_lines.svg
      :alt: (c) Isodensity (ρ-DROP)

   The same model molecule under the van der Waals sphere union, the DROP
   reference cavity (a), and the two main DROP level sets, :doc:`SvdW <svdw>`
   (b) and :doc:`isodensity <isodensity>` (c). The black line is the cavity
   surface, :math:`S=0`, and the dots are the atom centers.

For completeness, two further level sets are available: the :doc:`COSMO Fine Cavity (CFC-DROP) <cfc>`, and the :doc:`smooth SES-like surface (SSES-DROP) <sses>`, an experimental optimized power-distance moment surface with a solvent probe radius and smoothing length.
In all cases, the :doc:`SvdW <svdw>` should be preferred over both due to its supiriour computational cost, smoother PES, and bester testing.

.. dropdown:: Further level sets: CFC and SSES

   .. subfigure:: AB
      :subcaptions: below
      :gap: 8px
      :width: 70%
      :align: center
      :name: fig-drop-lsf-further

      .. image:: /_static/figures/lsf/cfc_no_circles_no_lsf_lines.svg
         :alt: (a) CFC-DROP

      .. image:: /_static/figures/lsf/sses_no_circles_no_lsf_lines.svg
         :alt: (b) SSES-DROP

      The model molecule under the :doc:`CFC <cfc>` (a) and :doc:`SSES <sses>`
      (b) level sets.

During the projection and weighting procedure, DROP applies smooth switching functions, handles multiple projection branches in concave regions, and uses spatial screening to keep the construction differentiable and efficient.
This makes the cavity suitable for continuum-solvation models that require smooth surface areas, polarization energies, and gradients, such as CPCM and pressure-based cavity terms.

.. _fig-drop-scheme:

.. grid:: 1
   :gutter: 2

   .. grid-item::

      .. image:: /_static/figures/drop/reference-cavity.svg
         :alt: Isodensity surface and van der Waals reference cavity
         :align: center
         :width: 65%

      **(a)** The isodensity surface is shown in blue and the hard-sphere van
      der Waals reference cavity in gray.

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item::

      .. image:: /_static/figures/drop/projection.svg
         :alt: Projection of reference quadrature points onto the isodensity surface
         :align: center
         :width: 80%

      **(b)** Closest-point projection of reference quadrature points onto the
      isodensity surface, :math:`S=0`.

   .. grid-item::

      .. image:: /_static/figures/drop/weight-mapping.svg
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

.. _drop-construction:

Construction
------------

Every DROP cavity is built the same way: 
- configure a level set function,
- :ref:`radius model <cavity-radii>` (used for the reference system), and
- hand both into the DROP constructor.

Only the LSF changes between :doc:`SvdW-DROP <svdw>`, :doc:`CFC-DROP <cfc>` and :doc:`ρ-DROP <isodensity>`; the settings below apply to the DROP scheme and thus all of them.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
         use moist, only : moist_cavity_drop_parameters_type

         type(cavity_type_drop) :: cavity

         call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, &
            & error=error, param=moist_cavity_drop_parameters_type(num_leb=194), ctx=ctx)
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_drop_options options;
         moist_init_drop_options(error, &options, sizeof options);
         options.nleb = 194;

         moist_cavity cavity = moist_new_drop_cavity(error, context, lsf, radii, &options);
         moist_delete(lsf);    /* cavity owns a copy */
         moist_delete(radii);

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import CavityDROP, CPCMRadii, DROPParameters, SvdW

         cavity = CavityDROP(
             context=context, lsf=SvdW(), radii=CPCMRadii(),
             parameters=DROPParameters(nleb=194),
         )

Omitting the settings object (an absent ``param``, ``NULL`` options, or no
``parameters``) selects the compiled defaults. The constructor copies both the
LSF and the radius model. Call ``update`` after every geometry change, and for
ρ-DROP after every density change as well.

Settings
--------

Settings can be supplied two ways:

- **Parameter file**: parameters  in a JSON or TOML file, read via ``load_file`` (e.g. ``grid.num_leb``, ``projection.level``). The full key set is registered in ``register_cavity_drop_entries``.
- **Settings object**: When constructing the ``moist_cavity_drop_parameters_type``, via the constructor ``param%new``, allows setting all parameters (``nleb``, ``tolerance``, ``proj_maxiter``, ``proj_level``, ``branch_weight_s``, ``rho_grid_h``, ``wleb_prune_level``). The C and Python options expose those seven plus ``do_fine`` (logical, default false; computes all available surface properties); verbosity and debug come from the run context.
- **Optional properties**: ``do_grid_density``, ``do_curvature``, ``do_normal``, ``do_r_iI`` and ``do_rho`` switch single surface properties on (default false). ``do_fine`` switches all five on.


Level Set Functions (LSF)
-------------------------

The LSF defines the implicit surface that the reference is projected onto.
Unlike the settings documented below, the LSF *and its shape parameters* are construction-level choices:
the caller builds a concrete LSF and passes it to the cavity constructor (the ``lsf_model`` argument in the Fortran
API).
They are not parameter-file keys and cannot be set through ``load_file``.

MOIST provides isodensity, SvdW, CFC, and SSES LSFs:

.. toctree::
   :maxdepth: 1

   svdw
   isodensity
   cfc
   sses

Grid Discretization
-------------------

``grid.num_leb`` (integer, default ``194``)
   Number of Lebedev quadrature points per atomic sphere. The cost of the cavity creation scales linearly with ``num_leb``; values between 110 and 350 are popular choices.
   **Must be one of the supported Lebedev orders**: 6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810.


Numerical Tolerance
-------------------

``tolerance`` (:math:`\varepsilon`, real, default ``1.0e-10``)
   Main tolerance; every other DROP threshold derives from it (tightest to loosest):

   - quadrature weight cutoff ``wleb_cut``, :math:`0.05\,\varepsilon`
   - LSF screening threshold, :math:`0.1\,\varepsilon` (handed to the LSF by the constructor); read by the isodensity LSFs and by SvdW and CFC with their fade off. With the fade on (the SvdW and CFC default) the atom reach is the end of the fade switching region and the threshold is not used, see :doc:`svdw` and :doc:`cfc`
   - projection convergence ``proj_tol``, :math:`\varepsilon`
   - branch-degeneracy separation, :math:`10\,\varepsilon`

   A nonzero ``switching.wleb_prune_level`` overrides ``wleb_cut`` with the lower bound of its switching region.
   See :cite:t:`wittmann2026drop`, Supporting Information Sec. C.2.c for the convergence behavior with respect to these thresholds.

   A value of ``1.0e-10`` is generally recommended; combined with extremely tight SCF (and/or optimization) settings, values of ``1.0e-12`` is recommended.
   This setting has little effect on the computational cost of constructing the DROP cavity.
   Values of approximately ``1.0e-14`` or smaller may, however, cause convergence problems in the projection because of finite precision of floating-point arithmetic.


Surface Projection
------------------


``projection.maxiter`` (integer, default ``150``)
  Maximum iterations for the projection optimizer.

``objective.alpha`` (real, default ``0.5``)
  Weight :math:`w_a` of the anchor term :math:`\tfrac{1}{2} w_a \|\mathbf r - \mathbf r^\circ\|^2` in the projection objective.
  It does not move a point that has a single projection, but branch weights and the admissible branch radius depend on it through :math:`\Phi`, so it does not cancel in concave regions.
  This should not be changed.

``projection.branch_search_rho2_slack`` (real, default ``-1``)
  Search window :math:`\Delta\rho^2` in :math:`\mathrm{bohr}^2`: a projection :math:`b` of an anchor is kept only if :math:`\rho_b^2 \le \rho_\mathrm{min}^2 + \Delta\rho^2`, where :math:`\rho_\mathrm{min}` is the closest projection found for that anchor.
  The window applies to every projection level, both to the branch search and to the candidate filter.

  - ``-1``: derived from the branch weights, :math:`\Delta\rho^2 = 2\,\Delta\Phi_\mathrm{max}/w_a`, with the admissible objective excess :math:`\Delta\Phi_\mathrm{max}` from `Branching`_. With the defaults this is about :math:`0.26\,\mathrm{bohr}^2`. Every branch outside it has a softmax weight below the floor.
  - :math:`\ge 0`: explicit window, independent of the branch weights.

Switching Functions
-------------------

Smooth step functions that fade surface contributions in and out so the cavity stays differentiable also in rare edge cases.
Two independent switches act on the integration weights, each keyed on a different geometric quantity (see the note below for the underlying rationale):

**Critical level set weight switch** (``f_crit``), a function of the level set gradient norm :math:`\|\nabla S\|`:

``switching.w_0ls_from`` / ``switching.w_0ls_to`` (real, defaults ``0.25`` / ``0.6``)
   Transition bounds of the critical level set weight switch (start, end).
   Contributions are fully suppressed below ``w_0ls_from`` and fully restored above ``w_0ls_to``.

**Focal/branching weight switch** (``f_foc``), a function of the smallest eigenvalue of the Lagrangian Hessian restricted to the tangent space (TRLH):

``switching.w_0tra_from`` / ``switching.w_0tra_to`` (real, defaults ``0.1`` / ``0.3``)
   Transition bounds of the focal/branching weight switch (start, end).
   A point's contribution is damped as its tangential curvature drops from ``w_0tra_to`` toward ``w_0tra_from``.

.. note::

  Although the level set function itself is differentiable away from nuclei, the resulting zero level set can become singular where :math:`S=0` and :math:`\nabla S=0`.
  To avoid ill-conditioned or non-unique projections in such critical regions, a critical-point switching function :math:`f_\mathrm{crit}(\|\nabla S\|)` smoothly attenuates contributions from the zero level set, which amounts to scaling the surface measure and the discretized integration weights.
  It is inactive in regular surface regions and suppresses a contribution entirely below ``switching.w_0ls_from``.

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
        - :math:`10^{-12}` to :math:`10^{-10}`
      * - 2
        - :math:`10^{-10}` to :math:`10^{-8}`
      * - 3
        - :math:`10^{-8}` to :math:`10^{-6}`
      * - 4
        - :math:`10^{-6}` to :math:`10^{-4}`
      * - 5
        - :math:`10^{-4}` to :math:`10^{-2}`
      * - 6
        - :math:`10^{-2}` to :math:`1`


Branching
---------

``branching.softmax_scale`` (:math:`\tau_\Phi`, real, default ``0.0025``)
  Softmax scale of the branch-weight model in concave regions; smaller values make the branch-weight distribution sharper.
  Branch fractions are a softmax of :math:`-\Phi/\tau_\Phi` over the objective values :math:`\Phi` of all competing projections (i.e. the projected distance), so equal objective values receive equal weights and a lower objective value smoothly dominates.
  Each branch fraction multiplies the integration weight of its projected point.
  The admissible objective excess of a branch is :math:`\Delta\Phi_\mathrm{max} = \tau_\Phi\ln(1/w_\mathrm{floor})`, where the floor :math:`w_\mathrm{floor}` is ``branching.weight_floor`` (default ``0.0``) or, when that is not positive, the derived quadrature weight cutoff.
  It sets the default search window ``projection.branch_search_rho2_slack``.


Spatial Screening
-----------------

Acceleration settings that do not change results, only cost.

``screening.cell_grid_full_scan_below`` (integer, default ``200``)
   Below this atom count the level set's candidate cell list collapses to a single full-scan cell.
   The cell list returns, for each batch of points, the atoms whose screening reach meets the batch; both modes return the same atoms, so this has no influence on the results, only on the cost of the LSF evaluations.

The former ``screening.cell_grid_fraction`` is retired and no longer read; the cell list sizes its cells from the atom reach.
