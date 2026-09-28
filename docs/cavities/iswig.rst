iSwiG Cavity
============

iSwiG is the switching Gaussian surface-discretization approach :cite:p:`lange2010swig` for van der Waals type cavitites.
Each atomic sphere is discretized using Lebedev quadrature grids.
Every grid point is assigned a smooth weight given by a product of error-function switches from the neighbouring spheres, so buried points fade out continuously.
Points whose switching value (or area) falls below a cutoff are removed.

The implementation supports analytic nuclear derivatives using adaptive radii :cite:p:`wittmann2025cpcm`.

Settings
--------

All three are fields of the optional settings object; omitting it selects the compiled defaults.

``num_leb`` in Fortran, ``nleb`` in C and Python (integer, default ``110``)
   Number of Lebedev quadrature points per atomic sphere; the primary accuracy/cost setting.
   **Must be one of the supported sizes**: 14, 26, 50, 110, 194, 302, 434, 590, 770, 974, 1202.
   Any other value raises an error, because each grid has a fitted Born ``zeta`` (the Gaussian width scale).

``cut_f`` (real, default ``1.0e-10``)
   Switching-value cutoff. A point is kept only if its switching function exceeds this value. Used when ``cut_a`` is not positive.

``cut_a`` (real, default ``0.0``)
   Area cutoff, bohr\ :sup:`2`. When greater than zero it replaces ``cut_f``: a point is kept only if its switched area (Lebedev area times switching value) exceeds this value.

Construction
------------

iSwiG needs a :ref:`radius model <cavity-radii>` and builds the discretization directly.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_cavity_iswig, only : cavity_type_iswig, new_cavity_iswig, &
            & moist_cavity_iswig_parameters_type

         type(cavity_type_iswig) :: cavity

         call new_cavity_iswig(cavity, ctx, radius_model=radii, error=error, &
            & param=moist_cavity_iswig_parameters_type(num_leb=194, cut_f=1.0e-10_wp))
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_iswig_options options;
         moist_init_iswig_options(error, &options, sizeof options);
         options.nleb = 194;
         options.cut_f = 1.0e-10;

         moist_cavity cavity = moist_new_iswig_cavity(error, radii, &options);
         moist_delete(radii);  /* cavity owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import CavityISwiG, CPCMRadii, ISwiGParameters

         cavity = CavityISwiG(
             radii=CPCMRadii(),
             parameters=ISwiGParameters(nleb=194, cut_f=1.0e-10),
         )

A positive ``cut_a`` replaces ``cut_f``.
