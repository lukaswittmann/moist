Solvents
========

MOIST contains a table of 180 solvents.
For each solvent it stores experimental bulk descriptors, a molecular geometry, and atomic charges and multipoles computed in three environments.

A solvent is looked up by its name or table ID.
Names ignore case and surrounding whitespace and accept the aliases listed in the :ref:`solvent table <solvent-table>`, including the CAS number.
IDs are frozen: an ID never changes meaning, and new solvents are appended.

.. code-block:: console

   moist solvent water

prints the descriptors of one solvent.

Available data
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Data
     - Content
     - Access
   * - :doc:`Descriptors <descriptors>`
     - experimental bulk properties and the values derived from the geometry
     - Fortran, C, Python, command line
   * - :ref:`Geometries <solvent-geometries>`
     - one structure per solvent
     - Fortran
   * - :ref:`Atomic charges <solvent-charges>`
     - Hirshfeld, RESP, MBIS and CHELPG charges
     - Fortran
   * - :ref:`Multipoles <solvent-multipoles>`
     - MBIS charges, dipoles, quadrupoles and octupoles
     - Fortran

Decalin (ID 94) has descriptors only: no geometry, charges or multipoles.

.. toctree::
   :maxdepth: 1

   descriptors
   table
   molecular
