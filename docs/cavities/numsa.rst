NUMSA Cavity
============

NUMSA integrates the solvent-accessible surface area of a van der Waals cavity on per-atom Lebedev grids with the smooth switching function of :cite:p:`im2003gbsw`, similar to the implementation at `grimme-lab/numsa <https://github.com/grimme-lab/numsa>`_.

.. note::

   NUMSA produces no surface discretization:
   ``update`` fills the per-sphere areas, ``total_area`` (bohr\ :sup:`2`) and ``area_grad(3, nat)`` (the analytic gradient of the total area).
   ``total_volume`` is set to zero and ``ngrid`` stays 0.

Settings
--------

Fields of the optional settings object; omitting it selects the compiled defaults.

``num_leb`` (integer, default ``110``)
   Lebedev points per atomic sphere.

``probe`` (real, default ``0.0``)
   Probe radius in bohr, added to every atomic radius.

``offset`` (real, default 2.0 Å = ``3.78`` bohr)
   Offset added to the neighbour-list cutoff radius, in bohr.

``smoothing`` (real, default 0.3 Å = ``0.567`` bohr)
   Width of the switching function, in bohr.

``tolsesp`` (real, default ``1.0e-6``)
   Grid points with accessibility weight at or below this value are skipped.

Command Line
------------

.. code-block:: none

   moist cavity numsa <coord> [--nleb N] [--radii RADIUSMODEL]

Construction
------------

NUMSA is exposed by the Fortran API only.

.. code-block:: fortran

   use moist_cavity_numsa, only : cavity_type_numsa, new_cavity_numsa, &
      & moist_cavity_numsa_parameters_type

   type(cavity_type_numsa) :: cavity

   call new_cavity_numsa(cavity, ctx, radii=radii, error=error, &
      & param=moist_cavity_numsa_parameters_type(num_leb=194))
   if (allocated(error)) error stop error%message

   call cavity%update(mol, error)
   if (allocated(error)) error stop error%message
   write (*, *) "Area:", cavity%total_area
