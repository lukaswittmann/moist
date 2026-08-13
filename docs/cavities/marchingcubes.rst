Marching Cubes Cavity
=====================

Marching cubes integrates the zero level set of a given level set function.

.. Note::

   This appraoch does not produce a surface discretization (grid points); it only integrates the cavity and yields the total cavity surface area and volume.
   Execution fills ``total_area`` (bohr\ :sup:`2`) and ``total_volume`` (bohr\ :sup:`3`).
   That makes it an additional independent numerical comparison.

The implementation lives in ``src/moist/cavity/marchingcubes.f90``; the cavity type is ``cavity_type_marchingcubes`` and the constructor is ``new_cavity_marchingcubes``.
The bare integration kernel ``integrate_surface_marching_cubes`` is public as well, for callers that hold an LSF but no cavity.


Settings
--------

The LSF model and the radii model are required constructor arguments.

``lsf_model`` (required)
   Level set function template, for example a configured ``moist_cavity_drop_lsf_svdw_type`` or ``moist_cavity_drop_lsf_cfc_type``.
   The cavity stores a copy and refreshes its geometry caches on every ``update``; the integrator clones it once per OpenMP thread.

``radius_model`` (required)
   Atomic radius model. Its radii feed the LSF and set the extent of the integration box.

``spacing`` (real, default ``0.2``)
   Finest grid spacing in bohr -- the primary accuracy/cost setting.
   Halving it roughly octuples the work in the refined region.
   For a single sphere, the default reaches about 0.2 % on the area and 0.4 % on the volume.

``obj_file`` / ``pqr_file`` (optional paths)
   When given, the triangle mesh produced during ``update`` is written as an ``obj_file`` (Wavefront OBJ mesh) or an ``pqr_file`` (one ``HETATM`` per triangle center).

Command Line
------------

.. code-block:: none

   moist cavity mc svdw <coord> [--spacing REAL] [--radii RADIUSMODEL] [--dump]
   moist cavity mc cfc  <coord> [--spacing REAL] [--radii RADIUSMODEL] [--dump]
   
The level set is set up as for ``moist cavity drop``; it thus contains all LSF-specific options.
``--dump`` writes ``cavity.obj`` and ``cavity.pqr``.
