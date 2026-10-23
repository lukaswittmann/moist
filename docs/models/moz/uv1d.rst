1D MOZ Model
============

``model_moz_1d_type`` (``model/moz/1d/``) is the 1D solute-solvent (UV) model: XRISM :cite:p:`chandler1972rism,hirata1981xrism` at infinite
dilution on a radial grid of its own, with the solvent response of the bulk solvent (:doc:`solvent`).
The energy is the solvation free energy; the response is its derivative with respect to
the host data of the solute terms. The nuclear gradient differentiates solute geometry
at fixed bulk solvent, radial grids and site parameters.

Energy
------

.. code-block:: fortran

   use moist_model_moz, only : model_moz_1d_type, new_moz_1d_model, model_moz_1d_parameters_type

   type(model_moz_1d_type), target :: model
   type(model_moz_1d_parameters_type) :: params

   call new_uniform_radial_pair(rgrid, kgrid, 8192, 0.015_wp, error)
   call new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, error, params)
   call model%potential%add(lj_gaff, error)
   call model%potential%add(charges, error)        ! e.g. monopole_type: host charges
   call model%update(mol, error)

   call model%new_coupling(coupling, error)
   call model%prepare_energy(coupling, error)
   ! answer the coupling requests (host charges) ...
   call model%get_energy(coupling, energy, error)  ! adds dG_solv, Hartree

The model copies the bulk solvent at construction. A solvent that is not yet solved is solved
on first need (``prepare_solvent``, also callable early to fail fast) with the ``solvent`` settings
and on the grid they name; the solvent of the caller is not touched.

``get_energy``

* builds the Ng-split UV tables from the solute potential, with the host answers on the coupling;
* takes the solvent response ``P_ba(k) = w_ba(k) + rho_b h_ba(k)`` on the model's k grid.
  Ordinary RISM uses ``susceptibility_at``: copied when the model and VV grids coincide,
  otherwise the direct radial transform of the stored ``h_vv(r)``. DRISM uses
  ``susceptibility_oz_at`` with its analytical charge response.
  The UV and VV grids need not be the same;
* restores the k -> 0 charge moment ``sum_b q_b P_ba ~ k^2`` of the solvent response at its
  lowest modes (a fit of ``S_a/k^2`` on 0.05-0.25/bohr for ordinary RISM, analytical
  moments for DRISM), which a net-charged solute would otherwise amplify through
  the ``4 pi/k^2`` pole of the long-range potential. Both use the same minimum symmetric
  correction, preserving susceptibility symmetry and correlations between
  charge-neutral combinations of sites;
* solves ``h = W c P`` with the closure, using the previous converged ``t_s`` as the starting
  guess at the same structure and grids when available, otherwise starting from ``c_s = 0``; evaluates
  ``dG = (1/beta) sum_ia rho_a int Phi(h_ia, c_ia) d^3r`` with physical ``h`` and ``c`` and adds it.

The model caches the corrected solvent response by the accepted VV solution's
generation and the UV reciprocal grid. It caches solute ``W`` by solute coordinates
and that grid. Changed host charges or potential tables reuse these operators;
solute motion rebuilds ``W``, a new bulk solution rebuilds the solvent response,
and a changed reciprocal grid rebuilds both. Grid nodes, weights, transform tags,
spacing and coordinates are compared exactly. These caches survive UV invalidation
and problem reconstruction; the sampling and charge-moment mathematics are unchanged.

A converged solve returns a field whose residual was measured below the tolerance; a solve that does
not converge is an error and leaves the energy unchanged. ``correlations(cr, tr, hr, error)`` returns
the physical correlations of the latest energy at the current structure; an update or explicit
``invalidate()`` drops them.
``thermodynamics(thermo, error)`` returns the free energy and partial molar volume of the same
solve in a ``model_moz_1d_thermodynamics_type`` record.

Energy and response share one UV solve: the model keeps the solved problem with its potential
tables, and a request whose tables equal them at the same structure reuses the solution. A second
energy, or a response after an energy, with the same host answers therefore does not solve again.
When host answers change the tables, the model retains the equation, transforms, filters
and solver storage on the same grid pair. By default it starts from the previous
converged field, adjusting ``t_s`` by ``beta*(u_lr_old-u_lr_new)`` to preserve physical
``t``. The guess is checked and solver history is reset. An unusable guess starts
cold; a feasible warm attempt that fails retries once from ``c_s = 0``.
``convergence.reuse_initial_guess = false`` selects cold starts while retaining storage.
A geometry update, changed grid pair, changed bulk or explicit invalidation rebuilds
the problem. Standalone calls start cold unless ``guess`` or
``reuse_initial_guess=true`` is supplied.

Response
--------

``get_response`` returns the derivative of the solvation free energy with respect to the host
data of the solute potential. For ``monopole_type`` that is
``atomic_charge_adjoint`` (``dg_dq``, Hartree per e); terms without host data, such as fixed or
table charges, contribute nothing and the response stays empty.

The potential owns the derivatives of its tables with respect to the host data
(``compute_adjoint``); the model supplies the derivatives of the energy with respect to the tables,
``w_sr``, ``w_lr`` and ``w_k`` (``dE = sum w_sr du_sr + sum w_lr du_lr(r) + sum w_k du_lr(k)``):

* the iterate is ``t_s`` with ``c_s = C(u_sr, t_s)``, ``h = c_s + t_s`` and ``c = c_s - beta u_lr``;
  the residual ``R = t_s - T(c_s)`` depends on ``u_sr`` through the closure and on ``u_lr(k)``
  through the equation, the energy on ``u_lr(r)`` only through ``c``;
* one transposed solve ``J^T lambda = dE/dt_s`` at the converged iterate, with
  ``dE/dt_s = a [dPhi/dh (1 + c') + dPhi/dc c']`` and ``a = rho_a 4 pi r^2 dr / beta``;
* ``w_sr = dc_s/du_sr [a (dPhi/dh + dPhi/dc) + mu]`` with ``mu = T'^T lambda``,
  ``w_lr = -beta a dPhi/dc`` and ``w_k = -beta W^T k2r^T(lambda) P^T``.

The ``tolerance`` setting is also the residual RMS of the transposed solve relative to its
right-hand side. The solve uses the UV preconditioner and starts from the previous ``lambda``.

The derivative is exact for the discrete equations, not only at the continuum stationarity of
the functional; the ``response_tables`` and ``response_charges`` tests check it against central
differences of the energy. For the alanine zwitterion in TIP3P water (8192 x 0.025 A, tolerance
1e-7, m003) a cold transposed solve costs about as much as the energy (3.0 s against 3.3 s on 8
threads), a warm-started one after a 1e-3 e change of the charges about a third (1.05 s).

Nuclear gradient
----------------

``get_gradient(coupling, response, gradient, error)`` adds the derivative of the reported
UV free energy to ``gradient`` in Hartree/bohr and returns the host response items.
It solves or reuses the same converged UV state as energy and response, then obtains
the potential-table and solute-intramolecular sensitivities from one transposed solve.

At fixed site parameters, each radial pair potential is independent of the solute
coordinates. Geometry enters the XRISM equation through
``W_ij(k) = sin(k d_ij)/(k d_ij)``. With ``u(k) = k2r^T(lambda)`` and the physical
``c(k) = r2k(c_s(r)) - beta u_lr(k)``, the packed ``W`` sensitivity sums both
off-diagonal contributions from ``u(k) : [dW(k) c(k) P(k)]``. A pure intramolecular
kernel contracts it into site-coordinate derivatives; coincident sites and the
zero wave number have zero derivative. No additional beta conversion is applied.

The host adds the contraction of returned ``atomic_charge_adjoint`` with its own
geometry-dependent charges. The model holds the solved bulk solvent, radial nodes,
quadrature weights and site parameters fixed. Both the energy and its gradient are
invariant under rigid translation and rotation; changing the bulk solvent or the
radial discretization is a separate derivative. Both RISM and DRISM bulk solvents
are supported; the selected bulk solution is held fixed during differentiation.

The ``moz_1d_gradient`` suite checks the equation and converged energy for all five
closures, signed and zero energy sensitivities, and complete KH model gradients
with fixed or geometry-dependent host charges for both RISM and DRISM solvents.
It also checks every Cartesian coordinate at two displacement steps, rigid motions,
cached energy preservation, and a gradient call without a preceding energy.

Thermodynamics
--------------

.. code-block:: fortran

   use moist_model_moz, only: model_moz_1d_thermodynamics_type
   type(model_moz_1d_thermodynamics_type) :: thermo

   call model%thermodynamics(thermo, error)

``free_energy`` is the reported solvation free energy in Hartree. ``partial_molar_volume``
is ``Vbar = kT kappa_T (1 - D)`` in bohr^3, with the bulk compressibility from :doc:`solvent`
and the direct-correlation moment ``D = sum_ia rho_a int c_ia d^3r`` of the UV solution.
The Coulomb tails of the physical ``c`` cancel in the sum over the sites of a neutral solvent.
The volume is a diagnostic of the converged solution. Both values are available after an
energy or response call and become unavailable after an update or explicit invalidation.

Settings
--------

``model_moz_1d_parameters_type`` reads and writes JSON or TOML; the solver controls of
:doc:`solvers` form the nested ``solver`` group and the bulk-solvent settings of :doc:`solvent`
the nested ``solvent`` group.

==============================  ========  =========================================
Key                             Default   Meaning
==============================  ========  =========================================
``closure``                     ``KH``    HNC, KH or PSE2-4
``convergence.tolerance``       1e-7      RMS residual of a converged solve
``convergence.maxiter``         128       Outer iterations
``damping.small_k``             true      Small-k damping of post-processed updates
``damping.small_k_k0``          1.0       Damping wave number (1/bohr)
``damping.small_k_min_weight``  0.15      Damping weight at k = 0
``solver.name``                 ``JFNK``  PICARD, DIIS, JFNK or HYBRID
``solver.jfnk.precond``         spectral  none or spectral
``solver.jfnk.precond_k0``      1.0       Spectral-filter wave number (1/bohr)
``solver.jfnk.precond_w0``      0.1       Spectral-filter weight at k = 0
``solvent.*``                   (see VV)  Settings of the solvent's VV solve
==============================  ========  =========================================

Validation
----------

A solute identical to the solvent molecule obeys the VV equation, ``h_uv = w c (w + rho h_vv)``,
so at infinite dilution ``h_uv = h_vv`` and the solvation free energy of a water molecule in water
equals the excess chemical potential of the solvent. moist and pysofft both give -3.258589 kcal/mol
for cSPC/E parameters on table water geometry (KH, 8192 x 0.015 bohr); pysofft prints the VV functional as ``DG_vv = rho mu``
(-0.0161, per bohr^3), the UV value per molecule. The ``moz_1d_uv`` test ``water_in_water`` keeps
the identity to 1e-7. In the same case the partial molar volume is the molecular volume ``1/rho``,
since ``D`` then equals ``rho C`` of the solvent and ``kT kappa_T = 1/(rho (1 - rho C))``; the test
keeps it to 1e-6.

TIP3P water (298 K, 0.0334 1/A^3, KH, 8192 x 0.025 A, tolerance 1e-10) with the explicit force fields
of pysofft's smoke tests reproduces fresh pysofft runs (``origin/pysofft``, d3c3d076):

=====================  =======================  =======================  ==========
Solute                 moist (kcal/mol)         pysofft (kcal/mol)       Difference
=====================  =======================  =======================  ==========
water                  -5.904507619             -5.904507620             2e-9
methane                12.298556535             12.298556535             5e-10
alanine (EEQ charges)  21.560060578             21.560060655             8e-8
alanine zwitterion     -81.323815675            -81.323815755            8e-8
alanine cation (+1)    -44.558273289            -44.558275999            2.7e-6
=====================  =======================  =======================  ==========

The cation's small energy shift follows from moist's symmetric charge-moment
correction; the legacy pysofft correction is rank one and does not preserve
susceptibility symmetry. The ``pysofft_reference`` unit test also covers a synthetic
anion obtained by reversing the cation's site charges. Both charge signs agree
within 1e-6 kcal/mol with independent pysofft calculations (``c15bd85a``) that retain
the neutral block ``R P R``, ``R = I - q q^T / (q^T q)``, and impose both corrected
charge blocks. Comparisons with the legacy energies allow a 5e-6 kcal/mol difference.

The same TIP3P calculations give these partial molar volumes (A^3):

=====================  ==========
Solute                 Vbar
=====================  ==========
water                  29.9401
methane                48.5207
alanine cation (+1)    71.0138
=====================  ==========

DRISM bulk solvents are supported; see :doc:`solvent` for configuration.
