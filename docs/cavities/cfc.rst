CFC-DROP Cavity
===============

The COSMO Fine Cavity (CFC) is a radii-based pseudo-density surface following :cite:t:`klamt2018cfc`, originally discretized via a marching tetrahedron algorithm.
Atom :math:`A`, with center :math:`\mathbf r_A` and radius :math:`R_A`, enters through its scaled displacement and reduced surface distance

.. math::

   \mathbf s_A = \frac{\mathbf r - \mathbf r_A}{R_A},
   \quad
   \sigma_A = \lVert\mathbf s_A\rVert - 1.

Atomic and pairwise terms build the pseudo-density that is given by

.. math::

   \mathrm{PD}(\mathbf r) = \sum_A e^{a_1\sigma_A}
   + \sum_{A<B} c\,\bigl(1 - \mathbf s_A\cdot\mathbf s_B\bigr)^m e^{a_2(\sigma_A + \sigma_B)}.

The level set is defined as :math:`S(\mathbf r) = -\log\mathrm{PD}(\mathbf r)`, so the interior is negative.

.. _fig-cfc-lsf:

.. grid:: 1 1 2 2
   :gutter: 2

   .. grid-item::

      .. image:: /_static/figures/lsf/cfc_no_circles.svg
         :alt: CFC level set with contour lines
         :align: center
         :width: 90%

      **(a)** Level set with contour lines.

   .. grid-item::

      .. image:: /_static/figures/lsf/cfc_no_lsf_lines.svg
         :alt: CFC surface against the van der Waals surface
         :align: center
         :width: 90%

      **(b)** Cavity surface against the van der Waals surface (dashed).

**Figure:** CFC level set of a model molecule. Blue shading marks the interior
(:math:`S<0`), red the exterior (:math:`S>0`); the black line is the cavity
surface, :math:`S=0`, and the dots are the atom centers.

MOIST differs from the original construction in two ways:

- **No atom truncation.** The original evaluates PD over only the five nearest atoms of each point.
  MOIST keeps every atom and pair term within reach; distant terms vanish smoothly through the `Fade`_ instead of being cut off.
- **DROP instead of marching tetrahedra.** The surface is discretized by DROP projection rather than a marching tetrahedron triangulation.
  DROP is much faster, gives a smooth potential energy surface, and is fully differentiable.

Optional settings:

``a1`` (:math:`a_1`, real, default ``-15.0``)
   Atomic-term exponent.

``a2`` (:math:`a_2`, real, default ``-9.0``)
   Pair-term exponent.

``c`` (:math:`c`, real, default ``5.0``)
   Pair-term coupling constant.

``m`` (:math:`m`, integer, default ``4``)
   Pair-term polynomial power.
   The generated kernel and its derivatives bake in ``m = 4``; any other value is refused by the C constructor, by a parameter file, and by every evaluation, with the fade on or off.

``fade_enabled`` (logical, default ``true``)
   Fade every atom and pair term smoothly to zero over a finite switching region past each atom; see `Fade`_.
   Off gives the untruncated model, screened at ``screening_threshold``.

``fade_start`` (:math:`a`, real, default ``27.0``)
   Fade onset in e-folds of :math:`u_A = e^{a_2\sigma_A}`, at :math:`\sigma_A = a/|a_2|` (3 radii past the atom surface at the default :math:`a_2`).

``fade_width`` (:math:`b`, real, default ``3.0``)
   Width of the switching region in e-folds, :math:`b/|a_2|` in :math:`\sigma_A` (a third of a radius at the default :math:`a_2`).

Fade
----

Instead of cutting atom terms at a screening threshold, CFC switches them off smoothly as part of the LSF, just like :doc:`SvdW <svdw>` does.
The switching region width is given in e-folds of the per-atom pair factor :math:`u_A = e^{a_2\sigma_A}`, so in :math:`\sigma_A` it scales with each atom's radius:

.. math::

   t_A = \min\!\left[1,\ \max\!\left(0,\ \frac{|a_2|\,\sigma_A - a}{b}\right)\right].

The weight is the C\ :sup:`4` polynomial smoothstep, mirrored:

.. math::

   w_A = 1 - S_4(t_A) = (1 - t_A)^5 \sum_{k=0}^{4} \binom{4+k}{k}\, t_A^k.

It is 1 up to the onset :math:`\sigma_A = a/|a_2|` and exactly 0 from its end :math:`\sigma_A = (a+b)/|a_2|` on, and its first four derivatives vanish at both edges of the switching region.
The fade multiplies whole terms, the atomic term of :math:`A` by :math:`w_A` and the pair term of :math:`A` and :math:`B` by :math:`w_A w_B`:

.. math::

   \mathrm{PD}_\mathrm{f}(\mathbf r) = \sum_A w_A\, e^{a_1\sigma_A}
   + \sum_{A<B} w_A w_B\, c\,\bigl(1 - \mathbf s_A\cdot\mathbf s_B\bigr)^m e^{a_2(\sigma_A + \sigma_B)},

in batched and per-point evaluation alike.
The end of the switching region, at distance :math:`R_A\,[1 + (a+b)/|a_2|]` from the atom center, is the atom reach, so ``screening_threshold`` is not used by CFC while the fade is on.
At the defaults, the switching region is from 3 to 3.33 radii past the atom surface, a reach of :math:`4.33\,R_A`.

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
