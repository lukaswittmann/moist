SvdW-DROP Cavity
================

The smooth van der Waals surface :cite:p:`wittmann2026drop` is the default DROP level set.
It recovers solvent-excluded-surface-like features while avoiding the geometric singularities and crevices of a plain sphere union.
Neighbouring atomic contributions are combined through a smooth one-/two-/three-body blend.
The default parameters reproduce a probe radius of around 1.4 Å.

Optional settings:

``blend_k`` (real, default ``5.5``)
   Smoothing ``k`` in the ``exp(-k * d)`` kernel.
   Larger values give sharper features; smaller values smooth crevices more aggressively.

``blend_1b`` (real, default ``1.0``)
   One-body smoothing.

``blend_2b`` (real, default ``0.0``)
   Two-body smoothing; the two-body term is off at the default.

``blend_3b`` (real, default ``3.0``)
   Three-body smoothing.

Construction
------------

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_cavity_drop_lsf_svdw, only : &
            & moist_cavity_drop_lsf_svdw_type
         use moist, only : moist_cavity_drop_lsf_svdw_param_type

         type(moist_cavity_drop_lsf_svdw_type) :: svdw

         call svdw%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=5.5_wp))

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_svdw_options options;
         moist_init_svdw_options(error, &options, sizeof options);
         options.blend_k = 5.5;

         moist_lsf lsf = moist_new_svdw_lsf(error, &options);

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import SvdW, SvdWParameters

         lsf = SvdW(parameters=SvdWParameters(blend_k=5.5))

Pass the level set to the DROP cavity constructor together with a radius model;
see :ref:`drop-construction`.
