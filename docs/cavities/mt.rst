Marching tetrahedra
===================

The marching-tetrahedra cavity triangulates an LSF isosurface and groups the mesh into compact, connected patches of nearly equal area, one quadrature point per patch.
Lloyd relaxation moves each patch seed toward its centroid, which evens out patch shapes and point spacing.
The type is ``cavity_type_mt`` in ``moist_cavity_mt``; the parameter schema is ``mt``.

Settings
--------

``spacing`` (default 0.2 bohr)
   Largest Cartesian cell side; sets the geometric resolution. A smaller ``target_area`` cannot recover features the mesh misses.

``target_area`` (default 0.35 bohr\ :sup:`2`) and ``area_tolerance`` (default 0.1)
   Every unswitched patch area lies within ``target_area * (1 +/- area_tolerance)``.
   Components too small for any integer patch count use their attainable mean (Fortran ``component_target``, ``area_exception``).
   If balancing misses the bound, MOIST warns and keeps the best connected partition; ``balanced`` is then false.

``w_0ls_from``, ``w_0ls_to``, ``w_0ls_p``, ``w_0ls_a`` and ``weight_cutoff``
   The DROP low-gradient switch, integrated over each patch at the triangle centroids, so the switched area converges with ``spacing`` independently of the partition.
   ``critical_weight`` is the area-averaged switch; patches at or below ``weight_cutoff`` are dropped from the quadrature, not from the mesh.

``gaussian_scale`` (default 4.9)
   Initial Gaussian width ``xi0 = gaussian_scale / sqrt(a)`` with the switched patch area ``a``.
   A width correction then sets the self terms, so results hardly depend on this value.

Width correction
   One Newton step per point makes the quadrature exact for a uniform charge layer:
   the point's self term plus its pair terms reproduce the layer potential at the point, integrated analytically over nearby triangles within a window of four patch sides.
   The resulting factor multiplies the switched width, so the low-gradient switch still suppresses charge at seams.
   Factors are limited to 0.8--1.25; MOIST warns if any point reaches the limit.

Normals are analytic LSF gradients at the patch points, or the switched mesh normal where the gradient vanishes.
Each point belongs to its nearest atom.
LSF screening is disabled while meshing, since its zero exterior sentinel would look like another isosurface.

Accuracy
--------

Conductor sphere of radius 3, ``spacing = 0.2``, unit charge at the center and at 60% of the radius along each axis, worst case:

.. list-table::
   :header-rows: 1

   * - Target area
     - Points
     - Energy error
     - Charge-density RMS error
   * - 2.0
     - 56
     - 0.22%
     - 1.1%
   * - 1.0
     - 113
     - 0.08%
     - 0.56%
   * - 0.4
     - 282
     - 0.08%
     - 0.39%
   * - 0.35 (default)
     - 323
     - 0.08%
     - 0.41%

Without the width correction the density error is about 2.4%, and 3--5% without Lloyd relaxation as well; DROP's Lebedev grids give 1.6% at 302 points.
With it, the error decreases with ``spacing`` down to 0.24% at 0.05 bohr, and the total induced charge is accurate to about 10\ :sup:`-4`.

Frozen geometry
---------------

The surface supports frozen-geometry energies and fields (CPCM, COSMO, pressure-volume).
It has no nuclear or LSF derivatives, and components needing principal curvatures are unsupported.
The update/build protocol and the supported LSFs are those of :doc:`DROP <drop>`.

Command line
------------

.. code-block:: console

   moist cavity mt svdw molecule.xyz --spacing 0.2 --target-area 1.0 --dump
   moist cavity mt cfc molecule.xyz --spacing 0.15 --target-area 0.8

The LSF and radius options are those of DROP. ``--area-tolerance``, ``--gaussian-scale`` and ``--weight-cutoff`` are also available; ``w_0ls_*`` is set through the APIs only. ``--dump`` writes XYZ, CSV and PQR; ``--grad`` is unsupported.

API
---

.. code-block:: python

   from moist import Context, MT, MTParameters, Structure, SvdW

   recipe = MT(lsf=SvdW(), parameters=MTParameters(spacing=0.2, target_area=1.0))
   cavity = recipe.build(context=Context())
   cavity.update(Structure([1], [[0.0, 0.0, 0.0]]))
   mesh_area = float(cavity.get("mesh_area"))

.. code-block:: c

   moist_mt_options options;
   moist_init_mt_options(error, &options, sizeof options);
   options.target_area = 1.0;
   moist_cavity cavity = moist_new_mt_cavity(error, NULL, lsf, NULL, &options);

.. code-block:: fortran

   call new_cavity_mt(cavity, radius_model=radii, lsf_model=lsf, error=error, &
      ctx=ctx, param=moist_cavity_mt_parameters_type(target_area=1.0_wp))
   if (allocated(error)) return
   call cavity%update(mol, error=error)

``CavityMT`` in Python takes ``lsf``, ``parameters`` and ``radii`` directly; ``Isodensity`` works as for DROP.
The C constructor copies the LSF and radii; NULL selects defaults and CPCM radii. ``fields()``, ``get()`` and ``results()`` expose mesh and patch diagnostics such as ``mesh_area``, ``mesh_volume`` and ``patch_area``.
The Fortran ``mesh`` and ``patches`` components keep the raw triangulation; ``mesh_area`` and ``mesh_volume`` describe it, while the standard cavity totals use the weighted quadrature.
