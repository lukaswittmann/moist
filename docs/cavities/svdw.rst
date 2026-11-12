SvdW-DROP Cavity
================

The smooth van der Waals surface :cite:p:`wittmann2026drop` is the default DROP level set.
It recovers solvent-excluded-surface-like features while avoiding the geometric singularities and crevices of a plain sphere union.
Neighbouring atomic contributions are combined through a smooth one-/two-/three-body blend.
The default parameters reproduce a probe radius of around 1.4 Å.

.. subfigure:: AB
   :layout-sm: A|B
   :subcaptions: below
   :gap: 8px
   :name: fig-svdw-lsf

   .. image:: /_static/figures/lsf/svdw_no_circles.svg
      :alt: (a) Level set with contour lines

   .. image:: /_static/figures/lsf/svdw_no_lsf_lines.svg
      :alt: (b) Cavity surface against the van der Waals surface (dashed)

   SvdW level set of a model molecule. Blue shading marks the interior (:math:`S<0`), red the exterior (:math:`S>0`); the black line is the cavity surface, :math:`S=0`, and the dots are the atom centers.

Optional settings:

``blend_k`` (:math:`k`, real, default ``5.5``)
   Smoothing in the :math:`e^{-k d_I}` kernel, with :math:`d_I = \lVert\mathbf r - \mathbf r_I\rVert - R_I` the signed distance to the sphere of atom :math:`I`.
   Larger values give sharper features; smaller values smooth crevices more aggressively.
   Must be finite and positive; any other value is refused, with the fade on or off.

``blend_1b`` (:math:`s_1`, real, default ``1.0``)
   One-body smoothing.

``blend_2b`` (:math:`s_2`, real, default ``0.0``)
   Two-body smoothing; the two-body term is off at the default.

``blend_3b`` (:math:`s_3`, real, default ``3.0``)
   Three-body smoothing.

``fade_enabled`` (logical, default ``true``)
   Fade every atom term smoothly to zero over a finite switching region, with the switching function :math:`f_\mathrm{fade}` of `Fade`_.
   The end of the switching region is the atom reach, so ``screening_threshold`` is not used by SvdW while the fade is on.
   Off gives the untruncated model with hard screening at ``screening_threshold``.

``fade_start`` (:math:`n_0`, real, default ``21.0``)
   Fade onset in e-folds of the three-body factor :math:`u_I^{(3)} = e^{-k d_I/3}`, at the distance :math:`3n_0/k` beyond the atom surface (11.45 Bohr at the default :math:`k`).

``fade_width`` (:math:`\Delta n`, real, default ``3.0``)
   Width of the switching region in e-folds; it ends :math:`3(n_0+\Delta n)/k` beyond the atom surface (13.09 Bohr at the default :math:`k`).
   The defaults keep the surface within about :math:`10^{-8}` Bohr of the untruncated model.

The three fade settings are set through the parameter type (its constructor, or a JSON/TOML file read with ``read_file``).

Fade
----

Instead of using a simple screening threshold, SvdW switches off atomic contributions smoothly.
The fade acts on the three-body factor, :math:`u_I^{(3)} \to f_\mathrm{fade}(t_I)\,u_I^{(3)}`, with the :math:`\mathrm{C}^{4}` polynomial switching function

.. math::

   t_I = \min\!\left[1,\ \max\!\left(0,\ \frac{d_I - 3n_0/k}{3\Delta n/k}\right)\right],
   \quad
   f_\mathrm{fade}(t) = (1 - t)^5 \sum_{j=0}^{4} \binom{4+j}{j}\, t^j.

:math:`f_\mathrm{fade}(t_I)` is 1 up to the switching region and 0 after it.
The other factors follow as powers, :math:`u_I^{(b)} = \bigl(u_I^{(3)}\bigr)^{3/b}`, so the per-atom factor of body level :math:`b` carries :math:`f_\mathrm{fade}^{3/b}`.
:math:`n_0` counts the e-folds of the slowest-decaying factor :math:`u_I^{(3)}` before the fade begins and :math:`\Delta n` the e-folds it takes; at :math:`k = 5.5` the default switching region runs from 11.45 to 13.09 Bohr past each atom surface.

.. figure:: /_static/figures/lsf/svdw_fade_band.svg
   :alt: Atomic factor and fade weight across the SvdW switching region
   :align: center
   :width: 60%
   :name: fig-svdw-fade

   Fade of one SvdW atom term at the default parameters.
   (a) Three-body factor :math:`u_I^{(3)}` and faded factor
   :math:`f_\mathrm{fade}\,u_I^{(3)}` against the distance from the atom surface;
   the top axis counts e-folds of :math:`u_I^{(3)}` and the shaded area is the
   switching region.
   (b) Switching function :math:`f_\mathrm{fade}` and its slope.

One e-fold is the distance :math:`3/k` (0.55 Bohr at :math:`k = 5.5`) over which :math:`u_I^{(3)}` drops by a factor :math:`e`; the top axis of :numref:`fig-svdw-fade` (a) counts them from the atom surface.
At the onset :math:`u_I^{(3)}` is already :math:`e^{-21} \approx 8\times10^{-10}` of its surface value, and across the shaded switching region :math:`f_\mathrm{fade}\,u_I^{(3)}` leaves :math:`u_I^{(3)}` and reaches exactly zero at its end, the atom reach.
Panel (b) shows :math:`f_\mathrm{fade}` and its slope, which vanishes with its first three derivatives at both edges of the switching region.

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
