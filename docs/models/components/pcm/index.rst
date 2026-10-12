Polarizable Continuum Models
============================

PCM components describe electrostatic solvation on a shared cavity grid. The
conductor models, :doc:`cpcm` and :doc:`cosmo`, use the Gaussian surface
interaction operator. The layer models, :doc:`iefpcm`, :doc:`ssvpe`,
:doc:`dpcm`, :doc:`lpb`, :doc:`slic` and :doc:`nonlocal`, solve boundary
equations using generated Gaussian or point kernels. Each model has its own
component type and constructor.

The host supplies the raw electrostatic potential and its requested
derivatives. For models using the normal derivative, MOIST forms
:math:`v_n=\mathbf n\cdot\nabla\phi` from the same Gaussian potential used for
:math:`\phi`. The host also contracts the ``potential_adjoint`` response,
including ``w_dphi_dr`` when present. See :doc:`/reference/coupling` for the
request and response shapes.

Before updating the cavity, each PCM enables the surface geometry it needs
on the model's owned cavity copy. DROP uses its ``do_normal`` and
``do_curvature`` settings; iSwiG already stores normals and computes
curvatures when required. User-requested properties stay enabled, and the
original cavity configuration is unchanged.

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

``solver`` chooses ``dense`` (default) or ``iterative``. The iterative solver
is restarted GMRES with block-Jacobi preconditioning. ``solver_tol`` defaults
to ``1e-8``, ``solver_maxiter`` to ``1000``, and ``solver_restart`` to ``100``.
``storage`` chooses ``auto`` (default), ``dense``, or ``matrix_free`` for the
operator; ``memory_budget_gib`` controls the automatic choice and defaults to
4 GiB. ``newton_tol`` (``1e-10``) and ``newton_maxiter`` (``50``)
control nonlinear models. The default dense solver does not use the iterative
controls.

The ``memory_budget_gib`` setting is a positive integer in GiB across Fortran,
C, Python, JSON/TOML files, and parameter tables (1 GiB = 1,073,741,824 bytes).

.. toctree::
   :maxdepth: 1

   cpcm
   cosmo
   iefpcm
   ssvpe
   dpcm
   lpb
   slic
   nonlocal
