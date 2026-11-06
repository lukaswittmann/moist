.. _solvent-molecular:

Molecular Data
==============

Besides the descriptors, MOIST stores one molecular geometry per solvent and the atomic charges and multipoles computed for it.
These data are accessible from Fortran only, through the ``solvent_system_type`` of ``moist_solvent``.

.. _solvent-geometries:

Geometries
----------

Each solvent has one molecular structure, stored in Ångström and converted to bohr on lookup.
It gives the derived descriptors and the solvent sites of the :doc:`MOZ models </models/moz>`.
In Fortran, ``new_solvent_system(system, id, error=error)`` from ``moist_solvent`` builds a ``solvent_system_type`` that carries the geometry; see :doc:`/reference/fortran`.

.. _solvent-charges:

Atomic charges
--------------

Atomic charges in e are tabulated for four charge models and three environments:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Model
     - Definition
   * - ``hirshfeld``
     - Hirshfeld partitioning
   * - ``resp``
     - ORCA hyperbolically restrained CHELPG (``a = 0.0005``, ``b = 0.1``)
   * - ``mbis``
     - minimal basis iterative stockholder partitioning
   * - ``chelpg``
     - electrostatic-potential fit on a 0.3 Å grid with 2.8 Å ``rmax``

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Environment
     - Calculation
   * - ``gas``
     - isolated molecule
   * - ``solvent``
     - CPCM with the solvent's own permittivity
   * - ``conductor``
     - CPCM with infinite permittivity

All values come from single points at the ωB97M-V/def2-TZVPP level (def2/J, RIJCOSX, TightSCF, DEFGRID3).
The iodine-containing solvents (IDs 22, 23, 24, 25, 106, 129, 130 and 131) use SF-X2C-1e with x2c-TZVPPall and x2c/J instead.

In Fortran, ``system%get_charges(environment, model, charges, error)`` returns them.
The MOZ potential terms ``coulomb_<model>_<environment>`` use them as solvent site charges; see :doc:`/models/moz/potential`.

.. _solvent-multipoles:

Multipoles
----------

MBIS atom-centered Cartesian multipoles up to the octupole are tabulated for the same three environments and at the same level of theory.
The moments are signed, symmetric charge moments in e bohr\ :sup:`rank` in the axes of the stored geometry:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Moment
     - Rank
     - Component order
   * - monopole
     - 0
     - :math:`q`
   * - dipole
     - 1
     - X Y Z
   * - quadrupole
     - 2
     - XX YY ZZ XY XZ YZ
   * - octupole
     - 3
     - XXX YYY ZZZ XXY XXZ XYY XYZ XZZ YYZ YZZ

The MBIS monopoles come from the multipole calculation and can differ slightly from the ``mbis`` atomic charges.
In Fortran, ``system%get_multipoles(environment, "mbis", multipoles, error)`` returns a ``solvent_multipole_data_type``.
The MOZ terms ``multipoles_mbis_<environment>`` are reserved for these data; they have no kernel yet.
