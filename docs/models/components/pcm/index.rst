Polarizable Continuum Model
===========================

PCM-family components extend the abstract ``pcm_base`` and share the same Gaussian-surface machinery.
The cavity must provide grid positions, Gaussian widths, and switching function values.

The base implementation assembles the surface interaction matrix :math:`A`, solves

.. math::

   A\mathbf q = -f(\epsilon)\boldsymbol\phi,

and evaluates :math:`E_\mathrm{PCM}=\tfrac{1}{2}\mathbf q^T\boldsymbol\phi`.
The dielectric scaling :math:`f(\epsilon)` distinguishes the concrete PCM variants.

The electrostatic potential can be computed from atomic charges or supplied by a QM host.
Matrix inversion, LU, Cholesky, and iterative solvers are available; an externally assembled matrix can also be supplied.
The common implementation provides the reaction potential, model surface weights, and nuclear gradient.

.. toctree::
   :maxdepth: 1

   cpcm
   cosmo
