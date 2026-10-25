3D MOZ Model
============

``model_moz_3d_type`` (``model/moz/3d/``) is the 3D solute-solvent (UV) model: 3D-RISM :cite:p:`beglov1997rism3d,kovalenko1998rism3d` at infinite
dilution on a volume grid, with the solvent response of the bulk solvent (:doc:`solvent`).
The energy is the solvation free energy, with optional pressure and electrostatic finite-size corrections. Host field and atomic
charge responses are implemented at fixed geometry and fixed bulk solvent. Nuclear gradients
include potential derivatives and grid motion, at fixed bulk solvent and reciprocal lattice.

Energy
------

.. code-block:: fortran

   use moist_model_moz, only : model_moz_3d_type, new_moz_3d_model, model_moz_3d_parameters_type
   use moist_math_grid_3d_cartesian, only : moist_math_grid_3d_cartesian_type, new_cartesian_point_grid

   type(model_moz_3d_type), target :: model
   type(model_moz_3d_parameters_type) :: params
   type(moist_math_grid_3d_cartesian_type) :: grid

   call new_cartesian_point_grid(grid, error, 64, 64, 64, 0.5_wp)
   call new_moz_3d_model(model, ctx, grid, solvent, error, params)
   call model%potential%add(lj_gaff, error)
   call model%potential%add(charges, error)
   call model%update(mol, error)

   call model%new_coupling(coupling, error)
   call model%prepare_energy(coupling, error)
   ! answer the coupling requests ...
   call model%get_energy(coupling, energy, error)  ! adds dG_solv, Hartree

The model is grid-neutral: it uses only the points, weights and k-points of its grid and the
grid's own transform engine (``grid%new_trafo``), so any point grid serves; the grid and its
transform are chosen together by the caller. Gaussian-width grids are refused.
An energy-only transform backend can implement just ``fft_r2k``, ``fft_k2r`` and
``destroy``. The transpose and geometry-derivative methods are optional; their defaults
return an unsupported-operation error when a response or gradient requests them.

``get_energy``

* solves the solvent copy when needed, with the ``solvent`` settings and on the grid they name;
* builds the Ng-split tables on the grid from the solute potential, with the host answers on the
  coupling: ``u_sr`` and ``u_lr(r)`` per solvent site, ``u_lr(k)`` with the k = 0 mode dropped;
* tabulates the solvent response ``P_ga(|k|) = w_ga + rho_g h_ga``, one matrix per distinct
  ``|k|``. The default ``susceptibility = "oz"`` interpolates only short-range ``cbar(k)``
  with a cubic spline, restores the Coulomb tail analytically, and reconstructs ``P`` through
  OZ at each target wave number. ``direct`` transforms the stored ``h_vv(r)`` and ``linear``
  interpolates the exported susceptibility between the VV grid's k nodes
  (see :ref:`moz-3d-susceptibility-convergence`);
* solves ``h_a(k) = sum_g c_g(k) P_ga(|k|)`` with the closure, through the grid's forward and
  backward transforms; the first solve uses the selected cold guess and subsequent fixed-geometry
  host updates use the previous converged fields. The map is affine in ``c``, so Jacobian
  products are exact;
* adds ``dG = (1/beta) sum_a rho_a int Phi(h_a, c_a) dV`` by the grid's quadrature, plus the
  pressure correction the ``correction`` setting selects (:ref:`moz-3d-pressure-corrections`) and the
  optional ``finite_size`` term (:ref:`moz-3d-finite-size`); the partial molar volume
  uses ``D = sum_a rho_a int c_a dV``.

Energies at the same structure and with the same tables reuse the kept solution.
``correlations(cr, tr, hr, error)`` and ``thermodynamics(thermo, error)`` return the physical
correlations, ``(ngrid, nv)``, and the free energy, partial molar volume and correction terms of
the latest energy; an update drops them.

At fixed geometry, changes in the host charges or EC field retain the equation, transform
plans and solver workspaces. Previous-solution guesses are enabled by default through
``convergence.reuse_initial_guess``. If the fitted Coulomb tail changes, the short-range
iterate is shifted by ``beta*(u_lr_old - u_lr_new)`` to preserve physical ``t = h - c``.
The guess is checked with the new closure and equation. Every solve starts a fresh solver
stage, clearing nonlinear history and Krylov recycling. An infeasible guess starts cold;
a failed warm attempt retries cold once. Each attempt has its own ``convergence.maxiter``
budget, and reported iterations include both attempts. Only converged fields can seed a
later solve. Geometry updates and invalidation discard the problem and its previous guess.
Setting ``convergence.reuse_initial_guess = false`` starts each changed-input solve cold
while still retaining fixed-grid resources. Identical-input calls continue to reuse the
already converged solution.

``convergence.cold_guess = "cs0"`` preserves the default cold start: set the short-range
direct correlation ``c_s = 0`` and map it through OZ to initialize ``t_s``. The optional
``"ts0"`` instead sets the short-range indirect correlation ``t_s = 0`` directly, avoiding
that initial transform. Both use the same closure and equations. Their iteration counts
can differ: ``ts0`` helps some cases and hurts others. A checked warm guess takes precedence;
an infeasible guess or failed warm attempt falls back to the selected cold start.

.. _moz-3d-finite-size:

Electrostatic finite-size correction
------------------------------------

``finite_size = "hunenberger"`` adds a leading net-charge background correction
to the infinite-dilution solvent charging free energy. The default is ``"none"``.
It can accompany either ``correction = "PC"`` or ``"PC+"``:

.. code-block:: fortran

   params%finite_size = "hunenberger"
   params%correction = "PC+"
   mol%charge = -1.0_wp  ! fixed total solute charge, also for a fixed-charge potential

.. math::

   \Delta G_{\mathrm{FS}} = -\frac{Q^2}{2} S\,M_{\mathrm{box}},
   \qquad M_{\mathrm{box}} = \frac{1}{V}\int_{\mathrm{box}}\frac{d^3r}{|r|},
   \qquad S = 1-\epsilon^{-1}
     = 4\pi\beta\rho\lim_{k\to0}\frac{\mathbf q^T\chi(k)\mathbf q}{k^2}.

Here ``chi = w + rho h`` is the bulk OZ susceptibility and ``q`` contains the
solvent site charges. Its regularized charge limit supplies the screening factor;
no experimental dielectric constant is substituted. ``M_box`` is the analytical
Coulomb average over a centred rectangular box. For a cube it is
``[6 asinh(1/sqrt(2)) - pi/2]/L = 2.380077363979553/L``. The model captures box
lengths as ``dr * [nx, ny, nz]`` after each successful Cartesian grid update.

The reference convention matters: this model already uses the isolated tail in
``c(r) = c_s(r) - beta u_lr(r)`` when evaluating the free energy. We derive the
remaining leading term from the omitted charge-response zero mode,
``delta h_a = -4 pi beta Q z_a(0)/V``, using ``dPhi/dh_a ~ beta q_a Q/(2r)``.
This gives the formula above. The usual ``-2.837297479480619/L`` Ewald constant
belongs to a fully periodic energy reference and cannot be inserted unchanged
into this reported physical-c energy. This adaptation follows the background
response in :cite:t:`kovalenko2000ions`, Eq. (20).

This implementation requires a Cartesian point grid, ``susceptibility = "oz"`` and
the analytical Ng Coulomb split. Bounded molecular grids do not cover the full
periodic cell, so this box-only formula is refused there. The solute
tail charges must sum to the fixed ``mol%charge`` within ``1e-6`` electrons. The term
is a registered thermodynamic quantity, cached at fixed charge, lattice and bulk
solution and cleared on update or bulk-response replacement. ``thermo%finite_size``
reports its Hartree value. It contributes no host response or nuclear gradient
for fixed total charge, reciprocal lattice and bulk solvent.

Only the leading net-charge term is included; finite-cavity and multipole terms
are outside this approximation. It is a Hunenberger-type continuum background
correction adapted to the current RISM energy reference, rather than the complete
:cite:t:`hunenberger1999ewald` scheme.

Host response
-------------

Prepare a response pass, answer its host requests, then call the model:

.. code-block:: fortran

   call model%prepare_response(coupling, error)
   ! answer q or phi on the coupling, as requested by the potential terms ...
   call model%get_response(coupling, response, error)

An energy call is optional: the response first obtains the converged UV solution for the
current host answers. Selected thermodynamic items seed the uncorrected energy and the direct
correlation moment; diagnostic-only items contribute no energy sensitivity. All seeds are
accumulated before one shared solve of ``J^T lambda = dE/dt_s``. PC and PC+ use the same formulas
and dependency graph as the energy and remain mutually exclusive.

The response borrows the completed primal solver's Krylov workspace and returns
it even if the adjoint fails. Energy, response and gradient passes therefore
retain one large Krylov allocation. Solver settings, statistics and the frozen
equation linearization remain separate. The next changed-input solve begins a
fresh stage as described above. Exact products do not allocate the two fields
used only by finite-difference products.

The equation transpose uses the grid engine's ``fft_k2r_adj`` and ``fft_r2k_adj``. It does not
assume that a forward transform is the transpose of an inverse. Cartesian transforms account
for the packed half-spectrum multiplicity; molecular transforms account for the nonuniform
real-space weights and reciprocal normalization. Real fields use the Euclidean pairing,
and complex fields use ``Re(sum(conjg(w_k)*du_k))``. The stabilized OZ Coulomb term is
differentiated through its charge projection, including the dropped zero mode.

The potential pullback returns ``atomic_charge_adjoint`` for host monopoles and
``potential_adjoint`` for EC fields, including the EC fitted-tail sensitivity. Fixed terms
return no host item. The response is in the units of derivatives of the reported Hartree
energy; quadrature and inverse-temperature factors are included exactly once.

This is a derivative on a fixed grid, at fixed geometry and fixed bulk solvent. Grid-dependent
accuracy, including convergence of a chosen molecular recipe, still belongs to the primal
problem. Nuclear, grid-motion and bulk-solvent derivatives are outside this response path.

Nuclear gradient
----------------

.. code-block:: fortran

   call model%prepare_gradient(coupling, error)
   ! answer q or phi, and dphi_dr for an EC field, as requested ...
   call model%get_gradient(coupling, response, gradient, error)

``gradient`` accumulates the MOIST contribution in Hartree/bohr. ``response`` returns the same
host-input sensitivities as a response pass. The host adds their contraction with its input
derivatives at fixed sampling points. MOIST uses ``dphi_dr`` to account for sampling-point motion;
including that motion again in the host contraction would count it twice.

For the discrete residual ``F = t_s - T(C(t_s))``, one shared adjoint solves
``J^T lambda = dE/dt_s``. The potential and equation layers accumulate explicit nuclear
derivatives and sensitivities to point positions, quadrature weights and the common phase
reference ``grid%get_kref()``. The selected thermodynamic items also seed the explicit
quadrature derivative of the energy and direct moment. One ``grid%get_volume_gradient`` call
contracts the combined volume accumulator into nuclear coordinates.

The abstract transform interface has ``accumulate_r2k_geometry`` and
``accumulate_k2r_geometry`` in addition to its field adjoints. Molecular engines differentiate
the NUFFT phases and forward quadrature weights; they use the same plans for spatial derivative
sums, without a nuclear Jacobian. The Cartesian FFT is constant under its domain's centroid
translation at fixed spacing and dimensions; its transform geometry pullbacks are zero.
Analytic reciprocal potential tables use the same phase reference and already supply its
potential sensitivity.

These are derivatives of the chosen discrete energy. Atomic recipes, retained point membership,
reciprocal mode counts and spacings are held fixed locally. A change in pruning or grid topology
is outside this derivative contract. Gaussian grids, bulk-solvent derivatives and Hessians are
outside the current 3D gradient implementation.

.. _moz-3d-pressure-corrections:

Pressure corrections
--------------------

The PC and PC+ corrections of :cite:t:`sergiievskyi2015pressure` remove the large bulk pressure of the integral-equation
solvent from the free energy of creating the solute's volume:

* ``dG_PC = dG - P Vbar`` and ``dG_PC+ = dG - (P - P_ideal) Vbar``;
* ``P = (n + 1)/2 rho kT - kT/2 rho^2 C``, ``P_ideal = rho kT`` and ``kT kappa_T`` come
  from the VV solution of the solvent (:doc:`solvent`), with ``n`` solvent sites and
  ``C = sum_ab c_ab(k = 0)``;
* the partial molar volume is ``Vbar = kT kappa_T (1 - D)`` with the direct-correlation
  moment ``D = sum_a rho_a int c_a dV`` of the 3D UV solution. The Coulomb tails of
  physical ``c`` cancel in the sum over the sites of a neutral solvent.

Both terms are computed for every 3D energy through ``moist_model_moz_thermodynamics``.
``correction = "none"`` (default), ``"PC"`` or ``"PC+"`` selects the one added. The host
response differentiates that same reported energy. Each term is linear in ``D`` and adds
``beta dE/dD`` (``P kappa_T`` for PC, ``(P - P_ideal) kappa_T`` for PC+) to ``dPhi/dc``
before the single transposed solve. See :doc:`thermodynamics` for the dependency graph.

.. _moz-3d-bulk-response:

Bulk response transfer
----------------------

The bulk solve keeps ``cbar(r)`` directly from its Ng-split iterate,
``c(r) = cbar(r) - beta u_lr(r)``. Its native DST-IV transform is cached with a cubic spline.
Radial quadrature supplies the additional zero and band-edge samples; the spline's derivative
is zero at the origin and uses the radial transform's derivative at the upper edge.
An unsplit bulk solve is separated for transfer with ``alpha = 1/bohr``.

For the supported one-species rigid-molecule XRISM equations, every site has density ``rho``:

.. math::

   W_{ab}(k) &= \operatorname{sinc}(k d_{ab}),\\
   C(k) &= \bar C(k) - \frac{\gamma(k)}{k^2}qq^T,
   \qquad \gamma(k)=4\pi\beta\exp[-k^2/(4\alpha^2)],\\
   [I-\rho W(k)C(k)]\chi(k) &= W(k).

Write ``W = 11^T + k^2 D``, evaluating ``D`` from the series of ``sinc(x)-1`` at small ``x``.
The OZ matrix then has the finite form

.. math::

   A(k) = I-\rho W(k)\bar C(k)+\rho\gamma(k)(D(k)q)q^T.

This gives the consistent limit ``chi(0) = 11^T / [1 - rho sum_ab cbar_ab(0)]``.
The same factorization solves ``A z = D q``, so ``z = chi q / k^2`` has a well-defined limit.
The 3D map uses this vector for the Coulomb contribution, which avoids dividing cancellation
errors in nearly equal response entries by ``k^2``. The reciprocal potential retains the
grid's convention of dropping the exactly zero Coulomb mode.

For DRISM, replace ``W`` by ``Wbar = W + rho zeta`` and ``D`` by
``Dbar = D + rho zeta/k^2`` in these equations. The theory evaluates the finite
``zeta/k^2`` limit analytically. The zero-mode susceptibility has the same form;
the charge response instead satisfies ``4 pi beta rho q^T z(0) = 1 - 1/dielectric``.

The target table is cached with the bulk solution identity, transfer mode and ordered
wavevector magnitudes. Bulk re-solves, potential changes, model invalidation and changes of
the reciprocal grid invalidate it. Host charge updates reuse it. An already solved bulk
solvent retains its chosen settings; its achieved residual is available through ``residual()``.

Settings
--------

``model_moz_3d_parameters_type`` reads and writes JSON or TOML; the solver controls form the nested
``solver`` group and the bulk-solvent settings of :doc:`solvent` the nested ``solvent`` group.

====================================  ==========  ==============================================
Key                                   Default     Meaning
====================================  ==========  ==============================================
``closure``                           ``KH``      HNC, KH or PSE2-4
``correction``                        ``none``    Pressure correction added: none, PC or PC+
``finite_size``                       ``none``    Leading charge correction: none or hunenberger
``susceptibility``                    ``oz``      OZ reconstruction, direct transform, or linear sampling
``convergence.tolerance``             1e-7        RMS residual of a converged solve
``convergence.maxiter``               128         Outer iterations per attempt
``convergence.reuse_initial_guess``    true        Previous converged fields at fixed geometry
``convergence.cold_guess``             ``cs0``     Cold start: map c_s = 0, or set t_s = 0 with ts0
``solver.name``                       ``JFNK``    PICARD, DIIS, JFNK or HYBRID
``solver.jfnk.precond``                none        none
``solvent.*``                         (see VV)    Settings of the solvent's VV solve
====================================  ==========  ==============================================

Grids
-----

Accuracy depends on every grid: spacing, box size and periodic images for the Cartesian grid,
quadrature and k-grid sizing for atom-centred grids. On the Cartesian grid the transforms are
exact inverses and the solve converges to the tolerance; with ``susceptibility = "direct"``,
cSPC/E parameters on table water geometry on 64^3 points at 0.5 bohr converges in 11 Newton steps
(dG -3.008 kcal/mol, partial molar volume 31.18 A^3).
Atom-centred grids with NUFFT transforms need a validated recipe: with one Lebedev order per atom
the solve stalls.

Validation
----------

The reverse-derivative tests check both NUFFT routes against independent point, quadrature
and phase-reference perturbations. Converged synthetic UV problems check geometry sensitivities
for HNC, KH and PSE2-4, including energy/moment seeds and a zero objective. Model gradient tests
use KH with no correction, PC and PC+: physical water on a Cartesian grid, a weak charged solute
in cSPC/E water with table geometry on a molecular grid with both phase-reference owners, and analytic host charges
and fields whose geometry derivatives are contracted with the returned response. Two displacement
steps and whole-system translations check the full reported energy, including the host term.
The coarse molecular recipe in these derivative tests does not establish solvation-energy
convergence for a strongly interacting solute.

On an atom-centred grid the transform engine and the reciprocal potential tables share the
phase reference ``grid%get_kref()``. Its sensitivity returns through the volume accumulator, with
the nuclear contraction owned by the grid.

.. _moz-3d-susceptibility-convergence:

Susceptibility convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~

``oz`` and ``linear`` use the same bulk RISM theory and differ in how the response is sampled
on the 3D reciprocal grid (:ref:`moz-3d-bulk-response`). Tight VV and UV residual tolerances
control iterative convergence; they do not bound susceptibility interpolation error.

A chloride refinement study used one Cl atom with charge -1 at the box centre, GAFF chlorine
Lennard-Jones parameters (sigma 6.5497007460 bohr, epsilon 0.00042039205925 Hartree), and cSPC/E
nonbonded parameters on the tabulated water geometry (solvent 175; :doc:`solvent`). The cSPC/E
constants were supplied explicitly to retain that geometry. Both VV and UV used KH, with
ordinary RISM at 298.15 K and 101325 Pa. The Cartesian grid stayed fixed at 64^3 points with
0.6 bohr spacing (38.4 bohr box edge). VV radial spacing was 0.015 bohr, VV tolerance 1e-11,
and UV tolerance 1e-9. Both pressure and finite-size corrections were ``none``.

==========  ==================  =======================  =======================  ===================
VV points   dk (1/bohr)         linear dG (kcal/mol)     oz dG (kcal/mol)         Absolute difference
==========  ==================  =======================  =======================  ===================
8192        0.02556635          -82.366980               -82.169426               0.197554 kcal/mol
16384       0.01278317          -82.228075               -82.169426               0.058650 kcal/mol
32768       0.00639159          -82.180519               -82.169426               0.011093 kcal/mol
==========  ==================  =======================  =======================  ===================

The unrounded ``oz`` energies span about 2e-7 kcal/mol. The shrinking difference as ``linear``
approaches this stable value supports susceptibility sampling as the dominant source of the
raw-energy difference in this case. Increasing the VV point count halves ``dk = pi/(N dr)``
and doubles the radial cutoff, so this study varies both effects. It establishes neither
3D spacing/box convergence nor the physical accuracy of the solvent model or charge correction.

Use the default ``oz``; ``linear`` approaches it only as the VV grid is refined. Check radial-grid and Cartesian-grid convergence separately when comparing energies;
compare raw energies before adding pressure or finite-size corrections.
