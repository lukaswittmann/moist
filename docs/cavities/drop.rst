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
The resulting surface grid inherits the smoothness of the underlying LSF and the projection/mapping scheme and thus provides consistent analytical derivatives.

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
In all cases, :doc:`SvdW <svdw>` should be preferred over both for its lower computational cost, smoother PES, and better testing.

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

Grid Orientation
----------------

``new_cavity_drop`` takes an optional ``grid_model`` argument, and ``drop svdw``, ``drop cfc`` and ``drop sses`` take ``--sym``, ``--orientation`` and ``--sym-tolerance``; they orient the reference grids by the molecule or by its point group, see :ref:`cavity-symmetry`.

The Lebedev nodes are the anchors of the projection: in the molecular orientation, anchor :math:`j` of sphere :math:`a` sits at

.. math::
   
   \mathbf r^\circ_j = \mathbf c_a + r_a \mathbf R_a \mathbf u_j, 

and the projected point follows it.
A rotation :math:`\delta\boldsymbol\omega_a` of the grid moves the anchor by

.. math::

   \delta\boldsymbol\omega_a \times (\mathbf r^\circ_j - \mathbf c_a)

and turns its sphere normal.
The derivatives carry the anchor motion through the projection, the reach gate and the switching function of the anchor, the routes a moving atom also takes, and the normal through the horizon switch and the weight-mapping Jacobian.
The reverse-mode gradient sums these into one torque per sphere, which the symmetry model maps to the atoms.

Settings
--------

Settings can be supplied in two ways:

- **Parameter file**: parameters in a JSON or TOML file, read via ``load_file`` (e.g. ``grid.num_leb``, ``projection.level``). The full key set is registered in ``register_cavity_drop_entries``; an unknown key is an error.
- **Settings object**: the constructor ``param%new`` of ``moist_cavity_drop_parameters_type`` sets ``nleb``, ``tolerance``, ``proj_maxiter``, ``proj_level``, ``rho_grid_h``, ``wleb_prune_level``, the horizon band ``w_0hor_off`` / ``w_0hor_on`` and ``branch``, a ``drop_branch_parameters_type`` with the fields ``dpdist_max``, ``reach_on`` and ``reach_off``. The C and Python options expose the same eleven settings as flat fields (``branch_dpdist_max``, ``branch_reach_on``, ``branch_reach_off``) plus ``do_fine`` (logical, default false; computes all available surface properties); verbosity and debug come from the run context.

**Optional properties**: ``do_grid_density``, ``do_curvature``, ``do_normal`` and ``do_atomarea`` switch single surface properties on (default false). ``do_fine`` switches all four on. The projection distance ``pdist`` of every point is always stored.


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

.. _drop-lsf-scale:

Level-set scale
~~~~~~~~~~~~~~~

Every LSF carries a positive constant multiplier ``scale`` (:math:`\alpha`): DROP projects onto :math:`S = \alpha\,S_\mathrm{LSF}`, and every value and derivative of the LSF carries :math:`\alpha`.
The surface, and its properties (normals, and curvatures, ...) do not depend on that.
DROP's critical-gradient switch and the projection tolerance read :math:`S` in absolute terms, so the defaults put :math:`\lVert\nabla S\rVert` near one on the surface, as for SvdW, whose single-atom field is a signed distance:

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - LSF
     - Default ``scale``
     - :math:`\lVert\nabla S\rVert` on the surface
   * - SvdW
     - ``1``
     - one for a single atom, median 0.88 to 0.99
   * - CFC
     - ``0.25``
     - :math:`0.25\,|a_1|/R_I` for a single atom (one for carbon), median near one
   * - SSES
     - ``0.088``
     - :math:`2\,\alpha\,(R_I + p)` on convex patches, :math:`1/(2(\bar R + p))` at a mean radius of 3 bohr
   * - Isodensity
     - ``1000`` (:math:`1/\rho_\mathrm{iso}`)
     - :math:`|\nabla\rho|/\rho_\mathrm{iso}`: about 2.1 for neutral molecules, 1.6 for anions and 3 for cations, nearly independent of :math:`\rho_\mathrm{iso}`

A smaller scale pushes healthy surface into the critical-gradient switch, which removes area and, at projection level 2, starts branch searches; a much larger one makes the Newton solves slower (SSES at scale 1 instead of 0.088 is 1.5 to 2.2 times slower at identical area).
This value should generally not be changed.


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
   - branch-degeneracy separation, :math:`10\,\varepsilon`, at least ``1e-6`` Bohr

   ``wleb_cut`` does not depend on ``switching.wleb_prune_level``: the prune switch is exactly zero below its switching region, so the points it removes fall under the cut on their own.
   See :cite:t:`wittmann2026drop`, Supporting Information Sec. C.2.c for the convergence behavior with respect to these thresholds.

   A value of ``1.0e-10`` is generally recommended; combined with extremely tight SCF (and/or optimization) settings, a value of ``1.0e-12`` is recommended.
   This setting has little effect on the computational cost of constructing the DROP cavity.
   Values of approximately ``1.0e-14`` or smaller may, however, cause convergence problems in the projection because of finite precision of floating-point arithmetic.


Surface Projection
------------------

``projection.level`` (integer, default ``2``)
   Projection strategy for the projection of the reference points.
   Every level first runs a local projection from their anchor; the levels differ in which anchors are then searched for further branches and how densely that search is done.
   Valid levels are 1 to 4 and are outlined in the table below.

   .. list-table::
      :header-rows: 1
      :widths: 10 50 40

      * - Level
        - Strategy
        - Recommendation
      * - 1
        - local projection only; an anchor whose local solve fails is searched from a coarse seed lattice (1.0 Bohr spacing)
        - Fastest; single-point energies and very large systems
      * - 2
        - searches anchors whose local point is not fully open or lies more than 0.5 Bohr away; seed lattice with 0.5 Bohr spacing, all seed filters
        - Default; sufficient for most calculations, including MD and geometry optimizations
      * - 3
        - searches every anchor; seed lattice with 0.3 Bohr spacing, tangent-cone filter only
        - Smoother PES; for MD, geometry optimizations, or tight tolerances and convergence criteria where level 2 is not enough
      * - 4
        - searches every anchor; seed lattice with 0.2 Bohr spacing and one refinement, no seed filters
        - Reference calculations; most expensive

   Higher levels find more branches in concave regions at a higher cost; seeder spacing, refinement and seed filters follow from the level.
   A search never goes beyond ``branching.reach_off`` and ends at the branch window behind the nearest fully open point (reach gate ``branching.reach_on`` / ``branching.reach_off``, see `Settings`_).
   The seed lattice only supplies starting points for the Newton solves: every point is an exact projection onto :math:`S=0`, and the level only decides which branches are found.

   .. dropdown:: How the branch search works

      .. figure:: /_static/figures/drop/seed-lattice.svg
         :name: fig-drop-seed-lattice
         :alt: Seed lattice on a smooth model cavity with concave regions: anchors on the reference spheres, adaptive lattice cells, surface samples, and the seeds of one anchor with two points
         :width: 100%

         Seed lattice on a smooth two-dimensional model cavity (interior blue, surface black).
         Anchors (circles) lie on the reference spheres (dotted); only the searched anchors (yellow), here those in the concave regions, span a search ball (dashed).
         **(a)** Lattice cells are split only inside the search balls and only where the surface may pass; every sign-change edge of an unsplit cell gives one sample (blue) on :math:`S=0`.
         **(b)** Zoom into a concave notch (box in a): each sample kept by the seed filters is moved by the osculating pre-step (green, at most one spacing), and Newton converges from these seeds to both projected points (stars), one on each wall.

      - **Seed lattice**: built once per geometry over the search balls of the anchors to be searched. Lattice cells are split while they hold a sign change or are not surface-free; every sign-change edge gives one sample, retracted onto the surface. With a symmetry model (:ref:`cavity-symmetry`) the lattice is turned with the molecule: along the molecular frame, or along the closest orientation the symmetry group maps onto itself, so the seeds then share the symmetry of the molecule.
      - **Seed filters** drop samples that cannot lead to a point with weight:
        
        - ``cone`` (the anchor lies outside the sample's tangent cone),
        - ``grazing`` (the anchor sees the patch edge-on), 
        - ``back_facing`` (the normal points into the anchor's sphere, only with the horizon switch on), and 
        - ``focal`` (beyond both focal centers of the patch).
        
        Level 2 uses all four, level 3 ``cone`` alone. The filters are first order in the lattice spacing and can drop the only seed of a point, so level 4 uses none.
      - **Search**: with a fully open local point, the seeds within its branch window are solved. Otherwise seeds are solved nearest first in rounds of growing size, and the ball grows towards ``branching.reach_off`` until a fully open point bounds it. An anchor whose local solve fails is searched the same way.
      - **Duplicates**: a seed stops once its Newton prediction lands within ``1e-4`` Bohr of a point already accepted for its anchor, and a saddle met again is escaped only once.

      The search counters are printed at verbosity 4.

``projection.batch_size`` (integer, default ``8``)
   Anchors per work chunk; each chunk is projected by one thread in batched Newton calls.

``projection.newton.max_step`` (real, default ``1.0``)
   Longest Newton step in Bohr; longer steps are scaled down before the line search.

``projection.newton.merit_margin`` (real, default ``1.0``)
   Margin of the penalty in the line-search merit function, :math:`\mu \ge |\lambda| +` ``merit_margin``.

``projection.maxiter`` (integer, default ``50``)
   Maximum iterations of a local projection; a solve that needs more is retried from the level's seeder.

``objective.alpha`` (real, default ``0.5``)
   Weight :math:`w_a` of the anchor term :math:`\tfrac{1}{2} w_a \|\mathbf r - \mathbf r^\circ\|^2` in the projection objective.
   This should not be changed.

Weight Pruning
--------------

If desired, near-zero weights can be smoothly attenuated to speed up the computation of electrostatic potential integrals (for QM coupling) in cases with large numbers of points with negligible contributions (e.g. large cavities with high Lebedev order).
Additionally, it can be used in cases where small weights could cause numerical issues.

``switching.wleb_prune_level`` (integer, default ``0``)
   Smoothly switches off points whose Lebedev weight after the closest-point
   scaling is near zero. The switch is a factor of the gate of the point, so it
   acts on the switching value :math:`f`. The switch region is derived from the
   level:

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

.. TODO: add docs once paper with figures is done

Spatial Screening
-----------------

Acceleration settings that do not change results, only cost.

``screening.cell_grid_full_scan_below`` (integer, default ``200``)
   Below this atom count the level set's candidate cell list collapses to a single full-scan cell.
   The cell list returns, for each batch of points, the atoms whose screening reach meets the batch; both modes return the same atoms, so this has no influence on the results, only on the cost of the LSF evaluations.
