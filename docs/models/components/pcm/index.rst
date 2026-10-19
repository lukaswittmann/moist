Polarizable Continuum Model
===========================

PCM-family components extend the abstract ``model_continuum_component_pcm`` and share the same Gaussian-surface machinery.
The cavity must provide grid positions, Gaussian widths, and switching function values.

The base implementation assembles the surface interaction matrix :math:`A`, solves

.. math::

   A\mathbf q = -f(\varepsilon)\boldsymbol\phi,

and evaluates :math:`E_\mathrm{PCM}=\tfrac{1}{2}\mathbf q^T\boldsymbol\phi`.
The dielectric scaling :math:`f(\varepsilon)` distinguishes the concrete PCM variants.

The base declares the ``gaussian_potential`` request: ``phi`` in every phase, and ``dphi_dr`` and ``dphi_dxi`` in the gradient phase and, on a cavity whose surface follows the density, in the response phase.
It emits the ``potential_adjoint`` item :math:`w_\phi = \partial E/\partial\boldsymbol\phi`, which for a stationary PCM is the surface charge :math:`\mathbf q`.
See :ref:`coupling-requests` and :ref:`coupling-response`.

Its surface weights make the energy fully variational on the :doc:`ρ-DROP isodensity cavity </cavities/isodensity>`; analytic nuclear gradients work on static and isodensity cavities alike.

Options
-------

All PCM variants share one settings object (``moist_pcm_parameters_type`` in Fortran, ``moist_pcm_options`` in C, :py:class:`moist.PCMParameters` in Python).
The solver settings are available from every API and, in Fortran, from a JSON/TOML parameter file (see :doc:`/reference/fortran`).
The dielectric constant :math:`\varepsilon` is a required constructor argument of each variant; see :doc:`cpcm` and :doc:`cosmo`.

``solver`` (enum, default ``cholesky``)
   Linear solver for :math:`A\mathbf q = -f(\varepsilon)\boldsymbol\phi`; the table below lists them by scaling with the grid size :math:`N`, fastest first.

``solver_tol`` (real, default ``1.0e-10``)
   **Only used with** ``solver = iterative``; the direct solvers ignore it.
   Convergence threshold (relative): the residual 2-norm must fall below ``solver_tol`` times :math:`\max(1, \lVert\mathbf{rhs}\rVert)`.

``solver_maxiter`` (integer, default ``50``)
   **Only used with** ``solver = iterative``; the direct solvers ignore it.

``external_matrix`` (real array ``(ngrid, ngrid)``, default none)
   Pre-assembled surface interaction matrix passed to the constructor or ``set_external_matrix``.
   Skips the internal assembly, not the solve. Fortran only.

The interaction matrix is dense with :math:`N = n_\mathrm{grid}` rows, so assembly costs :math:`N^2` time and memory.
The Cholesky and iterative solvers rely on a symmetric positive definite (SPD) matrix, LU and inversion do not.

.. list-table::
   :header-rows: 1
   :widths: 14 28 16 12 30

   * - ``solver``
     - Method
     - Scaling
     - SPD required
     - Notes
   * - ``iterative``
     - Jacobi-preconditioned conjugate gradient
     - :math:`N^2` per iteration
     - yes
     - The iteration count depends on conditioning. Uses ``solver_tol`` and ``solver_maxiter``.
   * - ``cholesky``
     - LAPACK ``potrf`` + ``potrs``
     - :math:`N^3`
     - yes
     - Default; fastest direct solver.
   * - ``lu``
     - LAPACK ``getrf`` + ``getrs``
     - :math:`N^3`
     - no
     - About twice the Cholesky cost.
   * - ``inversion``
     - LAPACK ``getrf`` + ``getri``, then :math:`A^{-1}\boldsymbol\phi`
     - :math:`N^3`
     - no
     - Forms :math:`A^{-1}` explicitly; about six times the Cholesky cost.

.. toctree::
   :maxdepth: 1

   cpcm
   cosmo
