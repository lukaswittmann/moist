CFC-DROP Cavity
===============

The COSMO Fine Cavity (CFC) is a radii-based pseudo-density surface following :cite:t:`klamt2018cfc`, originally discretized via a marching tetrahedron algorithm.
Atom :math:`I`, with center :math:`\mathbf r_I` and radius :math:`R_I`, enters through its scaled displacement and reduced surface distance

.. math::

   \mathbf s_I = \frac{\mathbf r - \mathbf r_I}{R_I},
   \quad
   \sigma_I = \lVert\mathbf s_I\rVert - 1.

Atomic and pairwise terms build the pseudo-density that is given by

.. math::

   \mathrm{PD}(\mathbf r) = \sum_I e^{a_1\sigma_I}
   + \sum_{I<J} c\,\bigl(1 - \mathbf s_I\cdot\mathbf s_J\bigr)^m e^{a_2(\sigma_I + \sigma_J)}.

The level set is defined as :math:`S(\mathbf r) = -\log\mathrm{PD}(\mathbf r)`, so the interior is negative.

.. subfigure:: AB
   :layout-sm: A|B
   :subcaptions: below
   :gap: 8px
   :name: fig-cfc-lsf

   .. image:: /_static/figures/lsf/cfc_no_circles.svg
      :alt: (a) Level set with contour lines

   .. image:: /_static/figures/lsf/cfc_no_lsf_lines.svg
      :alt: (b) Cavity surface against the van der Waals surface (dashed)

   CFC level set of a model molecule. Blue shading marks the interior
   (:math:`S<0`), red the exterior (:math:`S>0`); the black line is the cavity
   surface, :math:`S=0`, and the dots are the atom centers.


MOIST's implementation differs from the original construction in two ways:

- **No atom truncation.** The original evaluates PD over only the five nearest atoms of each point. MOIST keeps every atom and pair term within reach; distant terms vanish smoothly through the `Fade`_ instead of being cut off.
- **DROP instead of marching tetrahedra.** The surface is discretized by DROP projection rather than a marching tetrahedron triangulation. DROP is much faster, gives a smooth potential energy surface, and is fully differentiable.

Optional settings:

``a1`` (:math:`a_1`, real, default ``-15.0``)
   Atomic-term exponent.

``a2`` (:math:`a_2`, real, default ``-9.0``)
   Pair-term exponent.

``c`` (:math:`c`, real, default ``5.0``)
   Pair-term coupling constant.

``m`` (:math:`m`, integer, default ``4``)
   Pair-term polynomial power.
   The generated kernel and its derivatives bake in :math:`m = 4`; any other value is refused by the C constructor, by a parameter file, and by every evaluation, with the fade on or off.

``fade_enabled`` (logical, default ``true``)
   Fade every atom and pair term smoothly to zero over a finite switching region past each atom; see `Fade`_.
   Off gives the untruncated model, screened at ``screening_threshold``.

``fade_start`` (:math:`n_0`, real, default ``27.0``)
   Fade onset in e-folds of :math:`u_I = e^{a_2\sigma_I}`, at :math:`\sigma_I = n_0/|a_2|` (3 radii past the atom surface at the default :math:`a_2`).

``fade_width`` (:math:`\Delta n`, real, default ``3.0``)
   Width of the switching region in e-folds, :math:`\Delta n/|a_2|` in :math:`\sigma_I` (a third of a radius at the default :math:`a_2`).

Fade
----

Instead of cutting atom terms at a screening threshold, CFC switches them off smoothly as part of the LSF, just like :doc:`SvdW <svdw>` does.
The switching region width is given in e-folds of the per-atom pair factor :math:`u_I = e^{a_2\sigma_I}`, so in :math:`\sigma_I` it scales with each atom's radius.
The switching function is the same :math:`\mathrm{C}^{4}` polynomial :math:`f_\mathrm{fade}` as for SvdW:

.. math::

   t_I = \min\!\left[1,\ \max\!\left(0,\ \frac{|a_2|\,\sigma_I - n_0}{\Delta n}\right)\right],
   \quad
   f_\mathrm{fade}(t) = (1 - t)^5 \sum_{j=0}^{4} \binom{4+j}{j}\, t^j.

:math:`f_\mathrm{fade}(t_I)` is 1 up to the onset :math:`\sigma_I = n_0/|a_2|` and exactly 0 from its end :math:`\sigma_I = (n_0+\Delta n)/|a_2|` on, and its first four derivatives vanish at both edges of the switching region.
The fade multiplies whole terms, the atomic term of :math:`I` by :math:`f_\mathrm{fade}(t_I)` and the pair term of :math:`I` and :math:`J` by :math:`f_\mathrm{fade}(t_I)\,f_\mathrm{fade}(t_J)`:

.. math::

   \mathrm{PD}_\mathrm{f}(\mathbf r) = \sum_I f_\mathrm{fade}(t_I)\, e^{a_1\sigma_I}
   + \sum_{I<J} f_\mathrm{fade}(t_I)\,f_\mathrm{fade}(t_J)\, c\,\bigl(1 - \mathbf s_I\cdot\mathbf s_J\bigr)^m e^{a_2(\sigma_I + \sigma_J)}.

The end of the switching region, at distance :math:`R_I\,[1 + (n_0+\Delta n)/|a_2|]` from the atom center, is the atom reach, so ``screening_threshold`` is not used by CFC while the fade is used.
At the defaults, the switching region is from 3 to 3.33 radii past the atom surface, a reach of :math:`4.33\,R_I`.

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
