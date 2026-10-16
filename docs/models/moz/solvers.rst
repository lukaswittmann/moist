MOZ Solvers
===========

``src/moist/model/moz/solver.f90`` exports the shared solver interface,
factory, evaluator, iteration controller and numerical parameters.
The implementations live in ``solver/picard.f90``, ``solver/diis.f90``,
``solver/jfnk.f90`` and ``solver/hybrid.f90``.

The solver advances a field stored as ``(npts, npair)``. The model driver
chooses its representation and supplies an evaluator for
``R(gamma) = gamma - gamma_map(gamma)``. Grid transforms, potential splitting
and closure evaluations belong to that evaluator. The solver has no closure,
thermodynamics or grid dependency.

Supported solvers
-----------------

``PICARD``
   Damped fixed-point mixing: ``gamma_new = gamma - damping*R(gamma)``.

``DIIS``
   Modified direct inversion in the iterative subspace. Picard mixing fills
   the history, then a constrained residual Gram solve combines stored
   iterates and residuals. A singular history uses Picard mixing; excessive
   residual growth clears the history. Optional backtracking bounds the
   update relative to the fixed-point residual. ``MDIIS`` is an alias.

``JFNK``
   Jacobian-free Newton-Krylov iteration with restarted GMRES for
   ``J*delta = R`` and an Armijo line search for ``gamma - alpha*delta``.
   Jacobian-vector products use finite differences, or exact matrix-free
   products supplied by the evaluator. A failed exact linearization selects
   finite differences for that update. The solver never assembles a Jacobian.

``HYBRID``
   Start with MDIIS and switch to JFNK when the RMS residual of the input
   iterate is at or below ``hybrid_switch_rmse`` (default ``1e-3``).
   The switch is one-way within a solve or coupling stage. JFNK starts from
   that same input iterate and returns its line-search update; the MDIIS
   proposal at the handoff is discarded. This adds one residual evaluation
   at the handoff, with no additional evaluations on subsequent JFNK steps.
   The MDIIS history is released before JFNK allocates its Krylov workspace.
   Both components use their usual DIIS and JFNK numerical controls.

API and ownership
-----------------

``new_solver(solver, name, error, parameters)``
   Allocate the requested polymorphic solver. Names are case-insensitive;
   surrounding whitespace is ignored. The optional
   ``moist_model_moz_solver_parameters_type`` contains numerical controls.
   Unknown names and invalid controls leave ``solver`` unallocated.

``solver%outer_step(residual_fn, gamma, gamma_new, rmse, feasible, error)``
   Evaluate the current residual and propose an update. The RMS residual
   refers to the input iterate: ``sqrt(sum(R**2)/size(R))``. The driver
   accepts ``gamma_new`` only when ``feasible`` is true and no error is
   allocated. Input and output fields have the same nonempty shape and
   occupy separate storage. Outputs are defined only on success.

``solver%begin_stage()``
   Begin a solve or coupling stage. JFNK resets its forcing term and recycled
   basis while retaining allocated workspace. DIIS clears its history and
   update diagnostics while retaining storage; ``diis%reset()`` also releases
   that storage. HYBRID starts again in MDIIS with
   an empty history and resets JFNK stage state while retaining its workspace.

``solver%allows_postprocessing()``
   True for Picard, DIIS and the MDIIS phase of HYBRID. False for JFNK and
   the JFNK phase of HYBRID: the driver must preserve the update accepted
   by the line search. Query this after ``outer_step`` so a hybrid handoff
   protects the newly accepted JFNK update.

``solver%get_statistics()``
   Read stored DIIS vectors, inner GMRES iterations and available cumulative
   Krylov, residual and exact-product counts and timings. Negative values
   indicate an unavailable metric. ``prev_alpha`` and ``prev_step_maxabs``
   describe the last accepted numerical update.
   HYBRID combines residual-evaluation counts from both components; its
   Krylov counters and timings describe JFNK work.

``run_iterations`` controls the outer loop using a driver-owned
``moist_model_moz_solver_iteration_type``. Its ``advance`` hook computes and
accepts updates, ``observe`` reports successful steps, and ``checkpoint``
handles model output. The driver supplies the iteration limit, strict RMS
tolerance and checkpoint frequency. Failed attempts count as iterations but
are neither observed nor checkpointed.

The abstract ``moist_model_moz_solver_evaluator_type`` requires ``residual``.
It may also implement ``has_jvp``, ``linearize`` and ``jvp`` for exact
matrix-free products. The evaluator owns or borrows its problem context and
linearization scratch; no global problem pointer is required.

An evaluator may provide ``has_adjoint`` and ``adjoint`` for exact Euclidean
transpose products. ``solve_adjoint(evaluator, rhs, solution, error, tolerance,
max_iter, preconditioner, nthreads, guess, workspace_owner, restart, max_krylov)`` solves ``J**T*solution = rhs`` at its
frozen linearization through the existing JFNK/GMRES implementation. The
optional tolerance is an absolute unpreconditioned residual RMS (default
1e-10) and ``max_iter`` limits the linear refinement steps (default 32); an
optional left preconditioner must apply to the transpose system. ``nthreads``
sizes the OpenMP teams as in the solver parameters. An optional ``guess``, e.g. the solution of a
nearby system, replaces the zero start; a guess already within the tolerance is returned as it is.
The primal linearization stays fixed throughout the solve.
The optional ``restart`` and ``max_krylov`` arguments control the transpose basis
and the inner iteration budget per refinement; both default to 50. The 1D and
3D UV models expose these as ``solver.adjoint.restart`` and
``solver.adjoint.max_krylov``. They are independent of the nonlinear JFNK controls.
To bound shared storage, configure the restart for both phases; lowering only
``solver.jfnk.restart`` does not lower the response default. A smaller restart
may require more work or fail to converge on a difficult transpose system.

An optional ``workspace_owner`` is a completed primal solver. The adjoint borrows
its Krylov scratch through allocatable ownership transfers and returns it on
success or failure. Settings and statistics remain separate; recycled vectors
are invalidated when ownership changes. HYBRID forwards the transfer to its
embedded JFNK solver. Other primal solvers can retain the adjoint scratch for
later calls. Each model owns its storage; independent models do not share it.
Shared scratch retains the largest required restart capacity, while each solve
uses its own restart length. A changed field shape replaces the allocation.
Finite-difference scratch is allocated only when a finite-difference product is
needed, including fallback from an unavailable exact linearization.

JFNK borrows left preconditioners and fallback filters through
``set_preconditioner`` and ``set_fallback_filter``. Concrete operators
implement ``ready``, ``apply`` and ``describe`` and live with their models.
Their targets must outlive the attachment; omitting the argument detaches
them. Workspace release leaves borrowed operators intact. The same lifetime
rule applies to ``set_iteration_printer``.

Preconditioner contract
-----------------------

A concrete preconditioner extends
``moist_model_moz_solver_preconditioner_type`` and implements:

``ready(rmse)``
   Whether the built operator should be used. JFNK queries this once per
   Newton update with the previous update's input RMS residual; the first
   update after ``begin_stage`` receives ``huge(1.0_wp)``.

``apply(vflat)``
   Replace the flat vector by ``M_inverse * vflat`` in place. Flat indices
   follow Fortran column order: ``i + (pair-1)*npts``. The operator acts on
   the representation selected by the model's residual evaluator and must
   return finite values without retaining a pointer to the borrowed vector.

``describe(rmse)``
   Return a short status for iteration output.

The model builds or refreshes the operator before the Newton update and
owns its lifetime and scratch. The mathematical operator must be linear,
nonsingular and fixed for the entire inner solve, including restarts.
Scratch and diagnostic counters may change during ``apply``. Diagonal,
block, sparse, approximate-factorization and transform-based operators
can all use this interface. JFNK applies it on the left to both the
right-hand side and Jacobian products. Inner tolerances use the
preconditioned residual norm; the outer convergence measure and line-search
merit use the unpreconditioned model residual.

Variable or nonlinear preconditioning requires a flexible Krylov method,
which is not implemented here. The fallback filter is a separate borrowed
operator and is applied only to JFNK's Picard fallback direction. It may
be singular; acceptance depends on the measured nonlinear residual.
An attached fallback filter is applied without a ``ready`` query, so the
driver must attach a built filter. HYBRID
forwards both attachments to its JFNK component and advertises the
preconditioning capability even during MDIIS so drivers can prepare it.

1D preconditioners
------------------

``model/moz/1d/preconditioner.f90`` provides the two 1D operators. The model
owns them and lends the one it is configured with through
``select_1d_preconditioner(name, spectral, precond, kspace)``, after checking
the name with ``check_preconditioner_name`` against ``preconditioners_1d_vv``
(``none``, ``spectral``, ``kspace``) or ``preconditioners_1d_uv`` (``none``,
``spectral``). Both require the uniform DST-IV radial pair from
``new_uniform_radial_pair``.

``moist_model_moz_1d_preconditioner_spectral_type``
   ``setup(rgrid, kgrid, k0, w0, error)`` copies the radial pair and its
   transform and builds the filter ``w(k) = w0 + (1 - w0) k^2/(k^2 + k0^2)``
   (``k0 > 0``, ``w0`` in ``[0, 1]``); other pairs or parameters are
   rejected and leave the filter unbuilt. ``apply`` transforms each site pair
   r -> k, multiplies by ``w(k)`` and transforms back. ``nthreads`` sizes the
   team over site-pair columns; drivers pass the run context's thread count.
   Each thread owns transform scratch, and applications reuse one field-sized
   k-space buffer. It also serves as the Picard
   fallback filter. An unbuilt filter, a vector that is not a multiple of the
   grid size, or a forward transform failure leaves the vector unchanged.

``moist_model_moz_1d_preconditioner_kspace_type``
   ``bind(spectral)`` borrows the filter. ``update(wk, cr, rho, beta, uk_lr,
   dcdr, error)`` builds and LU-factorizes ``M(k) = I - w (alpha c rho)`` at
   every k-point, with
   ``c(k) = FBT(c) - beta u_lr(k)``. ``alpha(k)`` is the volume-weighted
   running average of ``dc/dgamma`` up to node ``nint(r_char/dr + 0.5)``,
   ``r_char = min(r_max, pi/k)``, clamped to ``[-1, 0]``; singular blocks
   become the identity. ``wants_update(rmse, iter)`` is true below
   ``rmse_activate`` (default 1e-3) on every ``refresh``-th iteration.
   ``ready`` is true below ``rmse_activate`` once factors exist for the
   filter's current grid; until then JFNK runs unpreconditioned, unlike the
   always-ready spectral filter. Once active, systems with more than 100 site
   pairs use the filtered block solve with a symmetrized packed result, and
   smaller systems use the spectral filter. ``nthreads`` sizes the per-k
   OpenMP teams; the transforms run outside them.

The 1D operators are adapted from ``lukaswittmann/moist_dev``, branch
``pysofft``, commit ``44874e7e95ebdd35e74c8398476dbc0710bf1ac7``, and
reproduce it on the uniform pair, including the rounding of the ``alpha``
node index. The 3D operators are not ported yet.

Numerical controls
------------------

Defaults follow the shared numerical controls in pysofft: Picard damping
0.1, DIIS history 12 with backtracking enabled, and JFNK restart and inner
iteration limits 50. JFNK uses an initial forcing term 0.5 and subsequent
bounds 0.001 and 0.1. A zero ``jfnk_fd_eps`` selects an automatic
finite-difference displacement. ``jfnk_exact_jvp`` defaults to true and
prefers exact products when available.
As in pysofft, the forcing term lags one update: the first two inner solves
of a stage use the initial term, and each later one uses the ratio of the
previous two input residuals. A zero-step stall tightens the next inner
solve by a factor of ten. At the lower forcing bound the first zero step is
still returned, so a driver can accept an input that meets its tolerance; a
second zero step in a row ends the stage as an infeasible update without an
error, where pysofft repeated the same zero step until its iteration limit.
An exactly zero residual returns the input unchanged without Krylov work.

Optional ``jfnk_cgs2``, ``jfnk_recycle`` and ``jfnk_deflation`` select
reorthogonalized classical Gram-Schmidt, GCRO-DR recycling and DEF1
coarse-space deflation. They default to off. Recycling starts after the
first restart cycle that fills the restart length or the inner-iteration
budget, and uses ``min(jfnk_recycle, jfnk_restart - 1, jfnk_max_krylov - 1)``
vectors. Deflation requires exact Jacobian-vector products, a driver
basis and an evaluator ``deflation_setup`` hook; finite-difference solves and
unsupported evaluators leave it inactive, and a deflated solve skips recycling. The driver supplies physical basis profiles and weights.

``moist_model_moz_solver_parameters_type`` extends
``moist_model_parameters_type``, so it reads and writes JSON or TOML files
and prints its values. File keys group the fields: ``verbosity``,
``nthreads``, ``picard.damping``, ``diis.depth``, ``diis.damping``,
``diis.backtracking``, ``diis.backtrack_tau``, ``diis.min_alpha``,
``diis.max_step_ratio``, ``hybrid.switch_rmse`` and ``jfnk.<name>`` for the
``jfnk_<name>`` fields below. Omitted keys keep their defaults. ``validate``
checks all controls, the backtracking and deflation controls only when they
are enabled; ``read_file`` and ``new_solver`` call it. All real controls are
dimensionless except the finite-difference displacement and residual
thresholds, which use the model field's units.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Meaning
   * - ``verbosity``
     - 0
     - Diagnostic detail; inner output also requires an attached printer.
   * - ``nthreads``
     - 0
     - OpenMP team size for the solver's own parallel regions; models should
       pass ``ctx%get_num_threads()``. Zero takes ``omp_get_max_threads()`` at
       each call and a negative count is rejected. Vectors below 65536 values
       run on one thread. ``solver%team_size()`` returns the effective size; a
       HYBRID passes its count to both components.
   * - ``picard_damping``
     - 0.1
     - Picard fraction, also used while DIIS fills its history.
   * - ``diis_depth``
     - 12
     - Number of stored iterates and residuals.
   * - ``diis_damping``
     - 0.2
     - Residual fraction used when DIIS backtracking is disabled.
   * - ``diis_backtracking``
     - true
     - Bound the extrapolated step relative to the current residual.
   * - ``diis_backtrack_tau``
     - 0.5
     - Multiplicative reduction of the extrapolation fraction.
   * - ``diis_min_alpha``
     - 0.001
     - Lowest fraction considered by the DIIS step bound.
   * - ``diis_max_step_ratio``
     - 200
     - Maximum update RMS divided by the input residual RMS.
   * - ``hybrid_switch_rmse``
     - 0.001
     - Input residual RMS at which HYBRID switches to JFNK.
   * - ``jfnk_restart``
     - 50
     - Krylov vectors in one GMRES restart cycle.
   * - ``jfnk_max_krylov``
     - 50
     - Total inner iterations allowed for one Newton update.
   * - ``adjoint_restart``, ``adjoint_max_krylov``
     - 50, 50
     - Transpose GMRES restart and inner iterations per refinement; independent of JFNK settings.
   * - ``jfnk_eta0``
     - 0.5
     - Relative inner tolerance for the first two Newton updates of a stage.
   * - ``jfnk_eta_min``, ``jfnk_eta_max``
     - 0.001, 0.1
     - Bounds for later relative inner tolerances.
   * - ``jfnk_fd_eps``
     - 0
     - Finite-difference displacement; zero selects the automatic scale.
   * - ``jfnk_exact_jvp``
     - true
     - Prefer exact evaluator products when available.
   * - ``jfnk_cgs2``
     - false
     - Use two-pass classical Gram-Schmidt instead of modified Gram-Schmidt.
   * - ``jfnk_recycle``
     - 0
     - Recycled Krylov dimension; zero disables GCRO-DR.
   * - ``jfnk_deflation``
     - false
     - Enable driver-supplied DEF1 coarse-space deflation.
   * - ``jfnk_deflation_rmse``
     - 0.001
     - Input residual RMS below which the coarse space may be activated.
   * - ``jfnk_stagnation_tol``
     - 0
     - Relative residual improvement required across a restart cycle;
       zero disables the early stagnation exit.
   * - ``jfnk_alpha_max``
     - 2
     - Maximum line-search fraction; one disables step extension.

.. code-block:: fortran

   use mctc_env, only: wp, error_type
   use moist_model_moz_solver, only: moist_model_moz_solver_type, &
      & moist_model_moz_solver_parameters_type, new_solver

   class(moist_model_moz_solver_type), allocatable :: solver
   type(moist_model_moz_solver_parameters_type) :: parameters
   type(error_type), allocatable :: error

   parameters%jfnk_restart = 30
   parameters%jfnk_max_krylov = 60
   call new_solver(solver, "JFNK", error, parameters)
   if (allocated(error)) error stop error%message

To use the hybrid with the default handoff threshold:

.. code-block:: fortran

   parameters%hybrid_switch_rmse = 1.0e-3_wp
   call new_solver(solver, "HYBRID", error, parameters)
   if (allocated(error)) error stop error%message

These shared numerical solvers can be used independently. Complete RISM
solves still require model residuals and driver integration.

Tests
-----

``test/unit/test_moz_solver.f90`` checks construction, parameter files
(nested keys, JSON and TOML round trips, validation on read), Picard mixing, DIIS
extrapolation and convergence, singular-history fallback, finite-difference
and exact JFNK products, nonlinear backtracking, borrowed preconditioners,
coarse-space algebra, concurrent solve contexts and iteration-controller failure ordering. It also
checks workspace reuse, restart, happy breakdown, CGS2, recycling (including
a budget shorter than the restart length) and
OpenMP reduction consistency. Dense-reference tests exercise complete DEF1
updates with and without a left preconditioner, including a coarse-only
right-hand side. A forcing-term test observes the tolerance used inside
GMRES, including the one-update lag, verifies tightening after a zero-step
stall, and that the second zero step at the lower forcing bound ends the
stage without an error.
Adjoint tests solve frozen nonsymmetric transpose systems, including a
matrix with a complex conjugate eigenvalue pair, with and without a left
preconditioner. They check the physical residual, unchanged primal
linearization, zero right-hand sides and failure propagation. A DIIS
divergence test checks history release and the Picard fallback update.
Hybrid tests cover RMS-based switching, the exact threshold, one-way phase
changes, stage reset, borrowed operators and diagnostic forwarding. Failure
tests cover invalid shapes, callback errors, non-finite values, singular
products, preconditioner scaling and workspace reconfiguration.
Shared-workspace tests check allocation identity across primal and adjoint
solves, independent contexts, restart and recycling changes, callback failure
cleanup and a subsequent primal solve against a fresh solver.

Numerical safeguards
--------------------

GMRES uses relative stopping and breakdown thresholds so multiplying a
valid left preconditioner by a nonzero scalar preserves the linear solve.
It rejects non-finite products, a preconditioner that annihilates a nonzero
residual, and singular triangular Krylov systems. Trial evaluation failures
are rejected and backtracked; failures at the current point or in Krylov
Jacobian products abort the update. As in pysofft, an unusable slope product
continues with slope -1, and a non-descent Newton direction is replaced by the
residual, searched with slope ``-||F||^2``. When the line search fails, the
filtered Picard fallback tries ``alpha = 0.25, 0.125, ...`` and accepts only a
decreasing residual.

DIIS stores the evaluator's residual directly, avoiding cancellation from
reconstructing it as the difference of rounded iterates. Changing the DIIS
history size resizes its workspace. Changing the GMRES restart length resizes
standalone scratch; shared primal/adjoint scratch retains sufficient capacity
for both solves.
DIIS backtracking limits the step-to-residual ratio; it does not evaluate
the nonlinear merit at candidate steps. Nearly dependent histories can
still be ill-conditioned even when the small LU solve is nonsingular.

The algorithms and eleven generic numerical tests are adapted from
``lukaswittmann/moist_dev``, branch ``pysofft``, commit
``44874e7e95ebdd35e74c8398476dbc0710bf1ac7``. DIIS completes the upper
triangle from the computed lower Gram triangle; a correlated-residual
reference test covers this correction.
