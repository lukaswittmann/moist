Pressure Models
===============

Pressure components add the work of an applied hydrostatic pressure to a continuum model.
Both take the pressure in hartree/bohr\ :sup:`3`; one GPa is ``3.39893e-5`` hartree/bohr\ :sup:`3`
(``moist.pyscf.GPA_TO_AU`` in Python).

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * -
     - :doc:`PV <pv>`
     - :doc:`GOSTSHYP <gostshyp>`
   * - Energy
     - :math:`PV` from the cavity volume
     - Gaussian potentials on the cavity surface, with amplitudes from the electron density
   * - Host requests
     - none
     - ``gaussian_moments`` in every phase
   * - Cavities
     - DROP (static and isodensity), iSwiG
     - DROP (static and isodensity), iSwiG

.. toctree::
   :hidden:

   PV <pv>
   GOSTSHYP <gostshyp>

- :doc:`pv`
- :doc:`gostshyp`
