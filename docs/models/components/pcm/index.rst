Polarizable Continuum Models (PCM)
==================================

PCM components describe electrostatic solvation on a shared cavity grid. The
conductor models, :doc:`CPCM <cpcm>` and :doc:`COSMO <cosmo>`, use the Gaussian surface
interaction operator. The layer models, :doc:`IEF-PCM <iefpcm>`, :doc:`SS(V)PE <ssvpe>`,
:doc:`DPCM <dpcm>`, :doc:`LPB <lpb>`, :doc:`SLIC <slic>` and :doc:`NLE <nonlocal>`, solve boundary
equations using generated Gaussian or point kernels. Each model has its own
component type and constructor.

The host supplies the raw electrostatic potential and its requested
derivatives. For models using the normal derivative, MOIST forms
:math:`v_n=\mathbf n\cdot\nabla\phi` from the same potential used for
:math:`\phi`, Gaussian or point. With the outward normal :math:`\mathbf n`,
:math:`v_n` is minus the normal component of the solute's electric field. The host also contracts the ``potential_adjoint`` response,
including ``w_dphi_dr`` when present. See :doc:`/reference/coupling` for the
request and response shapes.

Before updating the cavity, each PCM enables the surface geometry it needs
on the model's owned cavity copy. DROP uses its ``do_normal`` and
``do_curvature`` settings; iSwiG already stores normals and computes
curvatures when required. User-requested properties stay enabled, and the
original cavity configuration is unchanged.

.. toctree::
   :hidden:

   CPCM <cpcm>
   COSMO <cosmo>
   IEF-PCM <iefpcm>
   SS(V)PE <ssvpe>
   DPCM <dpcm>
   LPB <lpb>
   SLIC <slic>
   NLE <nonlocal>

- :doc:`cpcm`
- :doc:`cosmo`
- :doc:`iefpcm`
- :doc:`ssvpe`
- :doc:`dpcm`
- :doc:`lpb`
- :doc:`slic`
- :doc:`nonlocal`

Choosing a model
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 14 22 34
   :class: matrix

   * - Scenario
     - Model
     - Physical inputs
     - Requirements
   * - Routine solvation; the default choice
     - :doc:`CPCM <cpcm>`, :doc:`COSMO <cosmo>`
     - ``epsilon``
     - any grid cavity
   * - Dielectric response beyond the conductor scaling, e.g. low-:math:`\varepsilon` solvents
     - :doc:`IEF-PCM <iefpcm>`, :doc:`SS(V)PE <ssvpe>`
     - ``epsilon``
     - normals, curvatures
   * - Comparison with the original PCM
     - :doc:`DPCM <dpcm>`
     - ``epsilon``
     - normals, curvatures, potential normal derivative
   * - Salt solutions with ionic screening
     - :doc:`LPB <lpb>`
     - ``epsilon``, ``kappa``
     - normals, curvatures
   * - Charge-asymmetric solvation, e.g. anions versus cations
     - :doc:`SLIC <slic>`
     - ``epsilon``, ``alpha``, ``beta``, ``gamma``
     - normals, curvatures, potential normal derivative; smooth cavity (SvdW- or ρ-DROP); nonlinear solve
   * - Structured solvent with a nonlocal dielectric response
     - :doc:`NLE <nonlocal>`
     - ``epsilon``, ``epsilon_inf``, ``lambda``
     - normals, curvatures, potential normal derivative

The conductor models scale the conductor solution by :math:`f(\varepsilon)`, which is exact for :math:`\varepsilon\to\infty` and least accurate for low dielectric constants.
CPCM and COSMO differ only in this factor, :math:`(\varepsilon-1)/\varepsilon` versus :math:`(\varepsilon-1)/(\varepsilon+\tfrac12)`.
All layer models need normals and principal curvatures.
The potential normal derivative is :math:`v_n=\mathbf n\cdot\nabla\phi`, the derivative of the solute's electrostatic potential :math:`\phi` along the outward surface normal :math:`\mathbf n` at each surface point, i.e. minus the normal component of the solute's electric field.
For these models the host supplies the potential gradient (``dphi_dr``, :ref:`coupling-requests`) in addition to :math:`\phi`, and MOIST forms :math:`v_n` from it.

Physical inputs
---------------

The physical inputs describe the solvent; the solver controls are listed under `Numerical options`_.
Solvent values such as ``epsilon`` and the refractive index are tabulated in :doc:`/solvents/index`.

``epsilon``
   Static relative permittivity :math:`\varepsilon\geq1` of the solvent.
   An infinite value selects the conductor limit of CPCM and COSMO.

``kappa`` (bohr\ :sup:`-1`)
   Inverse Debye screening length of the ionic solution, :math:`\kappa\geq0` (LPB).
   Can be computed from the ionic_strength and temperature (see Fortran helper ``debye_kappa(epsilon, ionic_strength, temperature)``).

Slic specific:

``alpha``, ``beta``, ``gamma``, ``mu``
   SLIC interface perturbation :math:`h(E_n)=\alpha\tanh(\beta E_n-\gamma)+\mu` of the normal electric field :math:`E_n`.
   ``alpha``, ``gamma`` and ``mu`` are dimensionless, ``beta`` is in inverse atomic units of the electric field.
   ``mu`` defaults to :math:`\alpha\tanh\gamma`.
   These quantities should not be confused with other :doc:`solvent descriptors </solvents/index>`.

``phi_static`` (hartree/e)
   SLIC static potential; it adds :math:`\phi_\mathrm{static}Q` for the net solute charge :math:`Q`.

NLE specific:

``epsilon_inf``
   NLE high-frequency permittivity of the solvent.
   ``epsilon_inf`` is usually taken as the square of the refractive index, :math:`n_\mathrm{D}^2`.

``lambda`` (bohr)
   NLE correlation length of the solvent polarization.

``epsilon_solute``
   NLE permittivity inside the cavity; 1 for quantum-mechanical densities.

Numerical options
-----------------

All models accept ``moist_pcm_parameters_type`` in Fortran,
``moist_pcm_options`` in C, or ``PCMParameters`` in Python. Omission selects
the compiled defaults. Physical parameters, such as dielectric constant and
Debye inverse length, belong to each model's constructor.

The run context is optional at construction: pass ``ctx=`` in Fortran,
``context=`` in Python, or a context handle as the second C argument
(``NULL`` to omit it). A component added to a model keeps its explicit
context or inherits the model's context. Standalone evaluation requires
an explicit context.

Optional settings:

``solver`` (``iterative`` or ``dense``, default ``iterative``)
   Linear solver. ``iterative`` is restarted GMRES :cite:p:`saad1986gmres` with block-Jacobi preconditioning;
   ``dense`` factors the operator once and does not use the iterative controls below.

``solver_tol`` (real, default ``1e-8``)
   Convergence tolerance of the iterative solver.

``solver_maxiter`` (integer, default ``1000``)
   Maximum number of iterations of the iterative solver.

``solver_restart`` (integer, default ``100``)
   GMRES restart length.

``storage`` (``dense`` or ``matrix_free``, default ``dense``)
   Operator storage. ``matrix_free`` keeps no matrices and applies the kernels pair by pair.

``newton_tol`` (real, default ``1e-10``)
   Convergence tolerance of the Newton solve of nonlinear models (SLIC).

``newton_maxiter`` (integer, default ``50``)
   Maximum number of Newton iterations.

``point_potential`` (logical, default ``false``)
   Couples a Gaussian surface through point potential integrals: the component requests ``point_potential`` instead of ``gaussian_potential``, and the direct nuclear gradient uses the bare Coulomb operator.
   The surface operator stays Gaussian.
   The integrals are several times faster, while the results are nearly identical.
   On point surfaces the option has no effect, so a host with only point integrals can set it unconditionally.

``parameter_derivatives`` (logical, default ``false``)
   Adds the ``parameter_adjoint`` response item, the derivatives of the energy with respect to the physical inputs (e.g. ``dg_deps``).

The parameter printout shows ``solver`` and ``storage`` at verbosity 2.
Verbosity 3 and above also shows the solver and Newton tolerances and iteration
controls.
