Integrations
============

MOIST provides cavity, solvation and pressure models for use in quantum-chemistry
programs. Availability depends on the host program:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Program
     - Availability
   * - :doc:`PySCF <reference/pyscf>`
     - Included in MOIST's Python package.
   * - ORCA
     - Coming with the next major ORCA release.
   * - TURBOMOLE
     - Coming with the next major TURBOMOLE release.

PySCF
-----

The ``moist.pyscf`` interface attaches MOIST to an existing PySCF calculation
through ``mf.MOIST(...)``. It supports RHF, RKS, UHF and UKS, including density
fitting, nuclear gradients, and energy and gradient scanners.

Choose a DROP or iSwiG cavity and combine CPCM, COSMO, pressure-volume and
GOSTSHYP components. DROP supports SvdW, CFC and electron-isodensity surfaces.

See :doc:`reference/pyscf` for examples, supported settings and limitations,
and :doc:`recipe/installation` for installation instructions.

Integrating another program
---------------------------

MOIST exposes :doc:`Fortran <reference/fortran>`, :doc:`C <reference/c>` and
:doc:`Python <reference/python>` APIs for other host programs. Start with the
:doc:`API overview <reference/index>` and the :doc:`host coupling protocol
<reference/coupling>`.
