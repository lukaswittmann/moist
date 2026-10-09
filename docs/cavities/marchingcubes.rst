Marching Cubes Cavity
=====================

Marching cubes integrates the zero level set of a given level set function.

.. Note::

   This approach produces no surface discretization (no grid points); it only integrates the cavity.
   ``update`` fills ``total_area`` (bohr\ :sup:`2`) and ``total_volume`` (bohr\ :sup:`3`) but no grid fields (``ngrid`` stays 0), which makes it an independent numerical reference for any LSF.

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
   When given, the triangle mesh produced during ``update`` is written as a Wavefront OBJ mesh (``obj_file``) or a PQR file with one ``HETATM`` per triangle centroid (``pqr_file``).

Command Line
------------

.. code-block:: none

   moist cavity mc svdw <coord> [--spacing REAL] [--radii RADIUSMODEL] [--dump]
   moist cavity mc cfc  <coord> [--spacing REAL] [--radii RADIUSMODEL] [--dump]
   
The level set is set up as for ``moist cavity drop``; it thus contains all LSF-specific options.
``--dump`` writes ``cavity.obj`` and ``cavity.pqr``.

Construction
------------

Marching cubes is exposed by the Fortran API only.
It takes the same LSF and radius models as a DROP cavity.

.. code-block:: fortran

   use mctc_env, only : wp
   use moist_cavity_marchingcubes, only : cavity_type_marchingcubes, &
      & new_cavity_marchingcubes, moist_cavity_marchingcubes_parameters_type

   type(cavity_type_marchingcubes) :: cavity

   call new_cavity_marchingcubes(cavity, radius_model=radii, &
      & lsf_model=svdw, error=error, &
      & param=moist_cavity_marchingcubes_parameters_type(spacing=0.2_wp), ctx=ctx)
   if (allocated(error)) error stop error%message

   call cavity%update(mol, error)
   if (allocated(error)) error stop error%message
   write (*, *) "Area / volume:", cavity%total_area, cavity%total_volume
