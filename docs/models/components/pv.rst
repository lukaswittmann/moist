Pressure-Volume Component
=========================

The ``model_continuum_component_pv`` component adds a pressure-volume contribution :cite:p:`spooner2014compressed,zeller2025pressure`.

.. math::

   E_\mathrm{PV} = PV,

where :math:`P` is the applied pressure and :math:`V` the total cavity volume.
The component declares no host requests and requires an updated cavity that
provides ``total_volume``.

Its surface weights make the energy fully variational on the :doc:`ρ-DROP isodensity cavity </cavities/isodensity>`; analytic nuclear gradients work on static and isodensity cavities.

Parameters
----------

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Name
     - Unit
     - Description
   * - ``pressure``
     - hartree bohr^-3
     - Applied pressure. One GPa is ``3.39893e-5`` hartree/bohr^3.

Construction
------------

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_model_continuum_component, only : model_continuum_component_pv, &
            & new_component_pv

         type(model_continuum_component_pv) :: pv

         call new_component_pv(pv, pressure=3.39893e-5_wp)

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_component pv = moist_new_pv_component(error, NULL, 3.39893e-5);
         moist_add_model_component(error, model, pv);
         moist_delete(pv);  /* model owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import ModelComponentPV
         from moist.pyscf import GPA_TO_AU

         pv = ModelComponentPV(1.0 * GPA_TO_AU)

Adding the component to a model and driving the evaluation phases is the same for every component:
see :ref:`coupling-fortran` for Fortran, :doc:`/reference/c` for handle ownership, and :doc:`/reference/python` for the Python classes (:py:class:`moist.ModelComponentPV`).
