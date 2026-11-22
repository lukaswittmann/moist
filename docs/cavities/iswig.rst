iSwiG Cavity
============

iSwiG is the switching Gaussian surface-discretization approach :cite:p:`lange2010swig` for van der Waals type cavitites. This approach is the generalization of :cite:t:`york1999gaussian`.
Each atomic sphere is discretized using Lebedev quadrature grids.
Every grid point is assigned a smooth weight given by a product of error-function switches from the neighbouring spheres, so buried points fade out continuously.
Points whose switching value (or area) falls below a cutoff are removed.

The implementation supports analytic nuclear derivatives using adaptive radii :cite:p:`wittmann2025cpcm`.

Settings
--------

All three are fields of the optional settings object; omitting it selects the compiled defaults.

``num_leb`` in Fortran, ``nleb`` in C and Python (integer, default ``110``)
   Number of Lebedev quadrature points per atomic sphere; the primary accuracy/cost setting.
   **Must be one of the supported sizes**: 6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810.
   Any other value raises an error, because each grid has a fitted Born ``zeta`` (the Gaussian width scale), shared with :doc:`DROP <drop>`.

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

         call new_cavity_iswig(cavity, radius_model=radii, error=error, &
            & param=moist_cavity_iswig_parameters_type(num_leb=194, cut_f=1.0e-10_wp), ctx=ctx)
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_iswig_options options;
         moist_init_iswig_options(error, &options, sizeof options);
         options.nleb = 194;
         options.cut_f = 1.0e-10;

         moist_cavity cavity = moist_new_iswig_cavity(error, context, radii, &options);
         moist_delete(radii);  /* cavity owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import CavityISwiG, CPCMRadii, ISwiGParameters

         cavity = CavityISwiG(
             context=context,
             radii=CPCMRadii(),
             parameters=ISwiGParameters(nleb=194, cut_f=1.0e-10),
         )

A positive ``cut_a`` replaces ``cut_f``.

Grid Orientation
----------------

A symmetry model (``grid_model=`` in Fortran, ``--sym`` and ``--orientation`` on the command line) orients the grids by the molecule or by its point group; see :ref:`cavity-symmetry`.

In the molecular orientation, point :math:`i` of sphere :math:`a` sits at :math:`\mathbf p_i = \mathbf c_a + r_a \mathbf R_a \mathbf u_i`, with :math:`\mathbf u_i` a Lebedev node and :math:`\mathbf R_a` the rotation of the sphere.
The points are rigid on their sphere: a rotation :math:`\delta\boldsymbol\omega_a` of the grid moves a point by :math:`\delta\boldsymbol\omega_a \times (\mathbf p_i - \mathbf c_a)` and turns its normal by :math:`\delta\boldsymbol\omega_a \times \mathbf n_i`, while widths and Lebedev weights stay unchanged.
The reverse-mode gradient therefore collects one torque per sphere,

.. math::

   \boldsymbol\tau_a = \sum_{i \in a} \left[ (\mathbf p_i - \mathbf c_a) \times \bar{\mathbf p}_i + \mathbf n_i \times \bar{\mathbf n}_i \right],

with :math:`\bar{\mathbf p}_i` the adjoint of the point position, including the switching functions of the neighbouring spheres, and :math:`\bar{\mathbf n}_i` that of the normal.
The symmetry model maps the torques to the atoms through the molecular frame.
