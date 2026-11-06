.. _solvent-descriptors:

Descriptors
===========

The experimental descriptors follow the `Minnesota Solvent Descriptor Database <https://comp.chem.umn.edu/solvation/mnsddb.pdf>`_.

.. list-table::
   :header-rows: 1
   :widths: 15 45 20 20

   * - Symbol
     - Meaning
     - Unit
     - Temperature
   * - :math:`\varepsilon`, :math:`\varepsilon_r`
     - static relative permittivity
     - :math:`\varepsilon/\varepsilon_0`
     - 298.15 K
   * - :math:`n_\mathrm{D}`
     - refractive index at the sodium D line
     -
     - 293.15 K
   * - :math:`\alpha`, :math:`\beta`
     - Abraham hydrogen-bond acidity and basicity
     -
     -
   * - :math:`\gamma`
     - macroscopic surface tension
     - mN/m
     - 298.15 K
   * - :math:`\rho`
     - mass density
     - kg/m\ :sup:`3`
     - mostly 293.15 K

Deviating temperatures are noted below the :ref:`solvent table <solvent-table>`.
The molar mass, number density and volume per molecule follow from the geometry and the density.
A lookup returns all of them, in SI and atomic units; see :ref:`api-solvents` for the fields per language.
