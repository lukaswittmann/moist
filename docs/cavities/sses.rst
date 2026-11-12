Smoothed SES surface (SSES-DROP)
================================

SSES is an experimental smooth approximation to the solvent-excluded surface (SES), built from weighted blends of the atoms' solvent-accessible spheres.
It produces reentrant features, but finite smoothing shifts neck formation and closure relative to an exact rolling-probe surface.

.. subfigure:: AB
   :layout-sm: A|B
   :subcaptions: below
   :gap: 8px
   :name: fig-sses-lsf

   .. image:: /_static/figures/lsf/sses_no_circles.svg
      :alt: (a) Level set with contour lines

   .. image:: /_static/figures/lsf/sses_no_lsf_lines.svg
      :alt: (b) Cavity surface against the van der Waals surface (dashed)

   SSES level set of a model molecule. Blue shading marks the interior
   (:math:`S<0`), red the exterior (:math:`S>0`); the black line is the cavity
   surface, :math:`S=0`, and the dots are the atom centers.

Usage
-----

The geometric parameters are the probe radius :math:`p\geq0` (default 2.7 Bohr) and smoothing length :math:`\varepsilon>0` (default 0.5 Bohr).
CLI arguments and the Fortran API use Bohr.

.. code-block:: console

   moist cavity drop sses molecule.xyz --probe-radius 2.7 --smoothing-length 0.5 --dump

Definition
----------

The definition has two steps.
A hard construction reproduces the rolling-probe surface, and then a smoothing step makes it differentiable.

Blended spheres
~~~~~~~~~~~~~~~

Atom :math:`i`, with center :math:`\mathbf c_i` and radius :math:`r_i`, has a solvent-accessible (SAS) sphere of radius :math:`r_i+p`.
Weights :math:`\mathbf w` from the simplex :math:`\Delta=\{\mathbf w:w_i\geq0,\ \sum_iw_i=1\}` blend the atoms into one sphere.
The power distance of this blended sphere is the weighted mean of the atomic SAS power distances:

.. math::

   \sum_iw_i\left(\|\mathbf x-\mathbf c_i\|^2-(r_i+p)^2\right)
   =\|\mathbf x-\bar{\mathbf c}_w\|^2-R_w^2,

.. math::

   \bar{\mathbf c}_w=\sum_iw_i\mathbf c_i,

.. math::

   R_w^2=\sum_iw_i(r_i+p)^2-\sum_iw_i\|\mathbf c_i-\bar{\mathbf c}_w\|^2.

A point :math:`\mathbf x` is inside the cavity if a probe sphere centered at :math:`\mathbf x` fits inside at least one blended sphere, :math:`\|\mathbf x-\bar{\mathbf c}_w\|+p\leq R_w`:

.. math::

   F_0(\mathbf x)=\min_{\mathbf w\in\Delta}
   \left[\left(\|\mathbf x-\bar{\mathbf c}_w\|+p\right)^2-R_w^2\right].

The blends recover the patches of the rolling-probe SES:

- A single atom, :math:`w_i=1`, gives its bare sphere of radius :math:`r_i`.
- Blends of two atoms are the spheres through the circle where the two SAS spheres intersect. Shrunk by :math:`p`, they sweep out the torus of the probe rolling along both atoms, which forms the reentrant neck.
- Blends of three atoms all pass through the center of the probe that touches the three SAS spheres. Shrunk by :math:`p`, they form that probe's concave spherical patch.

The surface :math:`F_0=0` therefore matches the SES locally.
Whether a probe can actually reach a crevice is not tested; see the limitations below.

Smoothed field
~~~~~~~~~~~~~~

:math:`F_0` has kinks where the best blend jumps.
SSES smooths two places, both controlled by :math:`\varepsilon`:

.. math::

   F(\mathbf x)=\min_{\mathbf w\in\Delta}
   \left[\varphi\left(\|\mathbf x-\bar{\mathbf c}_w\|\right)
   -\sum_iw_i\varphi(r_i)
   +\sum_iw_i\|\mathbf c_i-\bar{\mathbf c}_w\|^2
   +\tau\sum_iw_i\log w_i\right],

.. math::

   \varphi(d)=d^2+2p\sqrt{d^2+\varepsilon^2},

.. math::

   \tau=2\varepsilon r_0,

.. math::

   r_0=1.5\ \text{Å}.

- :math:`\varphi(d)` replaces :math:`(d+p)^2-p^2`. At :math:`\varepsilon=0` the first three terms are exactly the bracket of :math:`F_0`. The square root removes the kink of :math:`\|\mathbf x-\bar{\mathbf c}_w\|` at the blended center.
- The entropy term :math:`\tau\sum_iw_i\log w_i` turns the hard minimum into a Boltzmann-weighted soft minimum at temperature :math:`\tau`.
- The length :math:`r_0` gives :math:`\tau` units of squared length. It is a fixed constant shared by all atoms, so :math:`\varepsilon` is the only smoothing control.
- Because :math:`\tau` is fixed, the field is separable: an atom far from :math:`\mathbf x` has a vanishing weight and does not change the field.

Negative values describe the cavity interior; :math:`F=0` defines the surface.

For a single atom the bracket reduces to :math:`\varphi(\|\mathbf x-\mathbf c_1\|)-\varphi(r_1)`, which vanishes at the bare radius.
An exact override replaces it with the bare sphere field:

.. math::

   F(\mathbf x)=\|\mathbf x-\mathbf c_1\|^2-r_1^2.

Evaluation
~~~~~~~~~~

The minimizing weights are Boltzmann weights:

.. math::

   w_i\propto\exp\left(-\frac{\|\mathbf x-\mathbf c_i\|^2-\varphi(r_i)
   +2\,\mathbf u\cdot(\mathbf x-\mathbf c_i)}{\tau}\right),

.. math::

   \mathbf u=\frac{p\,\mathbf m}{\sqrt{\|\mathbf m\|^2+\varepsilon^2}},

.. math::

   \mathbf m=\mathbf x-\bar{\mathbf c}_w.

The weights depend on each other only through the three-component vector :math:`\mathbf m`.
moist solves the fixed-point equation :math:`\mathbf m=\sum_iw_i(\mathbf m)(\mathbf x-\mathbf c_i)` by Newton's method.
The converged weights also provide the atomic partition.
The LSF interface supports spatial derivatives through fourth order, nuclear and radius derivatives, and their contractions.
Derivatives follow the converged equation for :math:`\mathbf m` analytically.

Properties and limitations
--------------------------

- For positive :math:`\varepsilon`, the weight optimization is strictly convex and the field is smooth. The zero surface can still contain critical points with :math:`\nabla F=0` during topology changes.
- Enclosed solvent pockets can occur, including pockets too small to admit a hard probe. SSES does not guarantee exact hard-probe accessibility.
- The entropy introduces multiplicity dependence: duplicating every atom in a multi-atom system :math:`k` times shifts the field by :math:`-\tau\log k`.
- Strongly heterogeneous radii and very small smoothing lengths require separate geometric and convergence checks.

Field scaling for DROP
----------------------

The positive numerical multiplier ``scale`` (default 1, CLI ``--lsf-scale``) sets :math:`S=\mathrm{scale}\,F` and scales all derivatives equally.
It preserves the zero surface, normals and curvatures, but affects DROP's numerical switches and projection tolerances.

DROP's critical-gradient switch transitions from zero to full weight between :math:`\|\nabla S\|=0.25` and :math:`0.6`.
Increasing the scale reduces suppression of regions with small raw gradients; genuinely critical points remain suppressed.
Projection tolerances are absolute in :math:`S`, so large scales tighten the effective tolerance on :math:`F`.
Excessive scaling can cause projection failure and dropped grid points when that tolerance falls below floating-point precision.
