MOZ Solvers
===========

The MOZ models solve ``R(gamma) = gamma - gamma_map(gamma) = 0`` for a field
stored as ``(npts, npair)``. The model supplies an evaluator for the residual
(and, where available, exact Jacobian and transpose products); grids,
potentials and closures stay in the model. The solvers live in
``src/moist/model/moz/solver.f90`` (parameters and the ``new_solver`` factory)
and ``solver/``: ``picard.f90``, ``diis.f90``, ``hybrid.f90``, ``jfnk.f90``
with its parts in ``solver/jfnk/``, and the shared contracts in ``solver/core/``.

Solvers
-------

``PICARD``
   Damped fixed-point mixing, ``gamma_new = gamma - damping*R(gamma)``.

``DIIS``
   Modified direct inversion in the iterative subspace
   :cite:p:`kovalenko1999mdiis`. Picard steps fill the history; a singular
   history falls back to Picard, and strong residual growth clears it.
   Aliases: ``MDIIS``, ``ANDERSON``.

``JFNK`` (default)
   Newton-Krylov :cite:p:`knoll2004jfnk` with restarted GMRES
   :cite:p:`saad1986gmres` and an Armijo line search; no Jacobian is assembled.
   The 1D, 3D and VV models provide exact Jacobian-vector products
   (``jfnk.use_exact_jvp = true``); ``false`` uses finite differences. For exact
   Newton, use a tight fixed forcing term (``eta0 = eta_min = eta_max``) and
   ``alpha_max = 1``. Aliases: ``NEWTON``, ``NEWTON-KRYLOV``.

``HYBRID``
   MDIIS until the input residual RMS reaches ``hybrid.switch_rmse``, then JFNK
   for the rest of the stage.

A solve converges when the RMS residual ``sqrt(sum(R**2)/size(R))`` falls
below the model's ``convergence.tolerance``.

Settings
--------

``moist_model_moz_solver_parameters_type`` reads and writes JSON or TOML.
``name`` selects the solver; each solver has its own group (``picard``,
``diis``, ``hybrid``, ``jfnk``), and keys are the group plus the field, e.g.
``jfnk.restart``. In a model input they sit under ``solver.``. Omitted keys
keep their defaults; threads and verbosity come from the run context.

==================================================  ============  ==================================================
Key                                                 Default       Meaning
==================================================  ============  ==================================================
``name``                                            ``JFNK``      PICARD, DIIS, JFNK or HYBRID
``picard.damping``                                  0.1           Picard fraction, also while DIIS fills its history
``diis.depth``                                      12            Stored iterates and residuals
``diis.damping``                                    0.2           Residual fraction when the step limit is off
``diis.step_limit``                                 true          Bound the extrapolated step by the residual
``diis.step_limit_tau``                             0.5           Reduction factor of the extrapolation fraction
``diis.step_limit_min_alpha``                       0.001         Smallest fraction tried by the step limit
``diis.step_limit_ratio``                           200           Maximum update RMS over input residual RMS
``hybrid.switch_rmse``                              0.001         Input residual RMS at which HYBRID switches to JFNK
``jfnk.restart``                                    50            Krylov vectors per GMRES restart cycle
``jfnk.max_krylov``                                 50            Inner iterations per Newton update
``jfnk.eta0``                                       0.5           Inner tolerance of the first two Newton updates
``jfnk.eta_min``, ``jfnk.eta_max``                  0.001, 0.1    Bounds of later inner tolerances
``jfnk.forcing_beta``, ``jfnk.forcing_power``       0.9, 2        Eisenstat-Walker ``beta*(rmse/rmse_prev)**power``
``jfnk.use_exact_jvp``                              true          Exact Jacobian products; false: finite differences
``jfnk.fd_eps``                                     0             Finite-difference step; 0 selects it automatically
``jfnk.cgs2``                                       false         Two-pass classical Gram-Schmidt instead of MGS
``jfnk.recycle``                                    0             GCRO-DR recycled dimension; 0 disables it
``jfnk.deflation``                                  false         DEF1 coarse-space deflation; 1D VV only, on there
``jfnk.deflation_rmse``                             0.001         Residual RMS below which deflation may start
``jfnk.stagnation_tol``                             0             Required improvement per restart cycle; 0 disables
``jfnk.line_search_c``                              1e-4          Armijo sufficient-decrease constant
``jfnk.line_search_tau``                            0.5           Backtracking factor
``jfnk.line_search_min_alpha``                      1e-4          Smallest step length tried
``jfnk.line_search_max_iter``                       10            Trial steps of the backtracking search
``jfnk.line_search_expand_ratio``                   0.25          Merit fraction required to extend a full step
``jfnk.alpha_max``                                  2             Largest step fraction; 1 disables extension
``jfnk.fallback_alpha``                             0.25          First step length of the Picard fallback
``jfnk.adjoint.restart``                            50            Restart length of the transpose (response) GMRES
``jfnk.adjoint.max_krylov``                         50            Inner iterations per transpose refinement
==================================================  ============  ==================================================

.. code-block:: fortran

   use mctc_env, only: error_type
   use moist_context, only: moist_context_type, new_context
   use moist_model_moz_solver, only: moist_model_moz_solver_type, &
      & moist_model_moz_solver_parameters_type, new_solver

   type(moist_context_type), target :: ctx
   class(moist_model_moz_solver_type), allocatable :: solver
   type(moist_model_moz_solver_parameters_type) :: parameters
   type(error_type), allocatable :: error

   call new_context(ctx, nthreads=0, verbosity=0)
   parameters%name = "HYBRID"
   parameters%jfnk%restart = 30
   call new_solver(solver, ctx, error, parameters)
   if (allocated(error)) error stop error%message

Response solves
---------------

``solve_adjoint`` solves the transpose system ``J**T*x = rhs`` at the converged
linearization with the JFNK/GMRES machinery; the UV models use it for the
energy response and the nuclear gradient. Its controls are the
``jfnk.adjoint`` group, independent of the solver used for the primal solve.
The adjoint borrows the primal solver's Krylov workspace instead of
allocating its own.

Preconditioners
---------------

JFNK accepts a left preconditioner (``ready``, ``apply``, ``describe``) that the
model owns and lends for a solve; the operator must be linear and fixed during
each inner solve, and it is also applied to the transpose system. The model
key ``solver.jfnk.precond`` selects it:

- 1D UV: ``none`` or ``spectral`` (radial small-k filter, default);
- 1D VV: ``none``, ``spectral`` (default) or ``kspace`` (block inverse in k space);
- 3D UV: ``none``.

The spectral filter also smooths JFNK's Picard fallback step.

Warnings
--------

Settings that are adjusted and recoveries of a whole solve are reported once
through the run context, regardless of the verbosity, e.g. a ``jfnk.recycle``
capped by ``restart`` and ``max_krylov``, a stall at ``eta_min`` after two zero
steps, an ignored or failed warm-start guess, or too few site pairs for
``kspace``. Recoveries inside a stage (Picard fallbacks, inner solves capped at
``max_krylov``, skipped deflation or recycling) are summarized in one note
after each stage of the iteration table.
