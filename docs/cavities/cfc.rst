CFC-DROP Cavity
===============

The COSMO Fine Cavity (CFC) is a radii-based pseudo-density surface following :cite:t:`klamt2018cfc`, originally discretized via a marching tetrahedron algorithm.
A pseudo-density :math:`\mathrm{PD}(\mathbf r)` is assembled from atomic and pairwise terms; the level set is :math:`-\log \mathrm{PD}(\mathbf r)`, so the interior stays negative.

Optional settings:

``a1`` (real, default ``-15.0``)
   Atomic-term exponent.

``a2`` (real, default ``-9.0``)
   Pair-term exponent.

``c`` (real, default ``5.0``)
   Pair-term coupling constant.

``m`` (integer, default ``4``)
   Pair-term polynomial power.
   The generated kernel and its derivatives assume ``m = 4`` and can thus not be changed out of the box.

Construction
------------

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_cavity_drop_lsf_cfc, only : &
            & moist_cavity_drop_lsf_cfc_type
         use moist, only : moist_cavity_drop_lsf_cfc_param_type

         type(moist_cavity_drop_lsf_cfc_type) :: cfc

         call cfc%new(param=moist_cavity_drop_lsf_cfc_param_type(a1=-15.0_wp))

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_cfc_options options;
         moist_init_cfc_options(error, &options, sizeof options);
         options.a1 = -15.0;

         moist_lsf lsf = moist_new_cfc_lsf(error, &options);

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import CFC, CFCParameters

         lsf = CFC(parameters=CFCParameters(a1=-15.0))

Pass the level set to the DROP cavity constructor together with a radius model;
see :ref:`drop-construction`.
