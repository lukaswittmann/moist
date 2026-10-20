Bulk Solvent
============

``solvent_vv_type`` (``model/moz/solvent/vv.f90``) is the bulk solvent of the MOZ models.
It is not a solvation model: it holds the solvent system and potential, solves the 1D
solvent-solvent (VV) equations, and keeps the results the UV models consume.
The 1D and 3D models copy it; the solve itself lives in ``model/moz/solvent/vv/``.

Solve
-----

.. code-block:: fortran

   use moist_model_moz_solvent_vv, only : solvent_vv_type, new_vv_solvent, &
      & solvent_vv_parameters_type, solvent_vv_thermodynamics_type
   use moist_model_moz_potential, only : lj_cspce, coulomb_cspce
   use moist_model_moz_solvent_water, only : new_spce_water_geometry

   type(solvent_vv_type) :: solvent
   type(solvent_vv_parameters_type) :: params
   type(solvent_vv_thermodynamics_type) :: thermo
   real(wp), allocatable :: chi(:, :)

   call new_spce_water_geometry(system%solv_mol)   ! SPC/E and cSPC/E share this geometry
   call new_vv_solvent(solvent, ctx, system, error)
   call solvent%potential%add(lj_cspce, error)
   call solvent%potential%add(coulomb_cspce, error)
   call solvent%update(error)

   call params%read_file("vv.toml", error)   ! optional; defaults otherwise
   call solvent%solve(error, params)
   call solvent%susceptibility(chi, error)
   call solvent%thermodynamics(thermo, error)

The solve uses ordinary RISM by default (``theory = "RISM"``; ``XRISM`` is an alias)
for one rigid molecular species with one site per atom; every site carries the molecular density.
``theory = "DRISM"`` selects the Perkyns-Pettitt dielectric correction for a polar solvent.
The correction is a mathematical port from RISMiCal through pysofft; its provenance
is recorded in ``solvent/theory/PROVENANCE.md``.

The bulk molecule must be electrically neutral. RISM and DRISM reject nonfinite
charges and any net charge exceeding
``1e-6 + 64 * epsilon(1.0_wp) * max(1, sum(abs(q)))`` in elementary-charge units.
The ``1e-6 e`` floor accepts small residuals from rounded charge tables, including
residuals at the decimal boundary with floating-point summation error. Such models
are treated as neutral in the small-k response; the supplied charges are retained.
Larger residuals are rejected.
An uncharged bulk molecule is valid for RISM; DRISM also requires a nonzero dipole.
This restriction applies to the bulk solvent only; UV solutes may carry net charge.

In UV or CLI input files, select the theory in the ``solvent`` group:

.. code-block:: toml

   [solvent]
   theory = "DRISM"

   [solvent.drism]
   dielectric = 78.355  # optional; default: the solvent's tabulated value
   smear = 0.6682  # bohr; exp(-smear^2*k^2/2)

Standalone VV parameter files such as ``vv.toml`` above use ``theory`` at the
root and ``[drism]`` for the dielectric and smearing settings.

The bulk ``bulk_theory_type`` has RISM and DRISM implementations. DRISM computes
the fixed correction ``zeta(k)`` once on the native grid, uses
``wbar = omega + rho*zeta`` in the OZ solve and adds ``zeta`` once to physical ``h``.
The same theory evaluates the correction analytically during susceptibility
reconstruction, including ``zeta/k^2`` at zero. Solute derivatives keep the bulk
state fixed. The low-wave-number charge response gives ``1 - 1/dielectric``.

``lj_spce`` is standard SPC/E: oxygen LJ with zero hydrogen LJ.
``lj_cspce`` is cSPC/E: the same oxygen parameters plus the Luchko hydrogen LJ modification.
``coulomb_spce`` and ``coulomb_cspce`` share the same charges, -0.8476 e on oxygen
and +0.4238 e on each hydrogen. Select the LJ and charge terms separately.
``new_spce_water_geometry`` supplies their common 1 Angstrom O-H bonds and
109.47 degree H-O-H angle. Library construction preserves the supplied geometry;
it does not replace the general solvent table.

* The tables come from the updated potential on the uniform DST-IV radial pair of the settings, Ng-split as described in :doc:`potential`.
* The iterate is the short-range indirect correlation ``t_s(r)``; the closure gives ``c_s(r)``, and per k-point ``h = (I - w c rho)^-1 w c w`` with ``c(k) = c_s(k) - beta u_lr(k)``.
* Coupling stages ``lambda = 1/n, ..., 1`` scale the short-range potential in the closure; the long-range part is always fully on.
* The evaluator provides exact Jacobian and transpose products and the Gaussian coarse space of the deflated Newton-Krylov steps (see :doc:`solvers`).
* Small-k damping filters the updates of solvers that allow post-processing (Picard and the MDIIS phase); JFNK keeps its line-search update.
* The solver verbosity and thread count come from the run context; the bulk solve runs once, so at verbosity 2 it prints its iteration table and a result block with the bulk thermodynamics.

A converged stage returns a field whose residual was measured below the tolerance: on the
converging step the update gets a fresh residual and replaces the input only when it is not worse.
A stage that does not converge is an error and leaves the solvent unsolved.
The default residual tolerance is ``1e-11``; explicit settings are kept. The achieved RMS
residual is checked before accepting the solution and is available through ``residual()``.
Any change of the potential afterwards (a new term, an update, the Ng split or the mixing rule), or its replacement by a potential with another change history, makes the solution stale: ``is_solved`` turns false and the output getters refuse.
A copy of the solving potential keeps the solution.

Settings
--------

``solvent_vv_parameters_type`` reads and writes JSON or TOML; the solver controls of :doc:`solvers` form the nested ``solver`` group.

===============================  ===========  ==================================================
Key                              Default      Meaning
===============================  ===========  ==================================================
``theory``                       ``RISM``     RISM (XRISM alias) or DRISM
``drism.dielectric``              solvent      DRISM relative dielectric constant; unset takes the solvent table value
``drism.smear``                   0.6682       DRISM Gaussian smearing width (bohr)
``closure``                      ``KH``       HNC, KH or PSE2-4
``grid.npts``, ``grid.dr``       8192, 0.015  Points and spacing (bohr) of the radial pair
``convergence.tolerance``        1e-11        RMS residual of a converged stage
``convergence.maxiter``          128          Outer iterations per stage
``convergence.nlambda``          1            Coupling stages
``damping.small_k``              true         Small-k damping of post-processed updates
``damping.small_k_k0``           1.0          Damping wave number (1/bohr)
``damping.small_k_min_weight``   0.15         Damping weight at k = 0
``solver.name``                  ``JFNK``     PICARD, DIIS, JFNK or HYBRID
``solver.jfnk.precond``          spectral     none, spectral or kspace
``solver.jfnk.precond_k0``       1.0          Spectral-filter wave number (1/bohr)
``solver.jfnk.precond_w0``       0.1          Spectral-filter weight at k = 0
``solver.jfnk.deflation``        true         Deflated Newton-Krylov steps (on for the VV solve)
``solver.jfnk.deflation_sigma``  2.35         Gaussian width (bohr) of the coarse profiles
``solver.jfnk.deflation_nfun``   2            Coarse profiles per site pair
===============================  ===========  ==================================================

MDIIS from a cold start diverges for water; use JFNK, or HYBRID with JFNK taking over early.
Stiff cases such as cold neon with PSE4 need the spectral preconditioner, in pysofft as well.

Outputs
-------

``radial_grid(rgrid, kgrid, error)``
   Copies of the radial pair of the solution.

``correlations(cr, tr, hr, error)``
   Physical ``c(r)``, ``t(r) = h - c`` and ``h(r)``, packed ``(npts, npair)``.

``susceptibility(chi, error)``
   ``chi(k) = w(k) + rho h(k)``, packed, on the reciprocal grid; the solvent response of the UV models.

``susceptibility_at(k, chi, error)``
   ``chi`` at any wave numbers in ``[0, pi/dr]`` of the VV grid, for UV grids of their own and 3D lattices:
   ``w(k)`` exactly from the structure and ``h(k) = 4 pi sum_i dr r_i^2 h(r_i) j0(k r_i)`` from the stored
   ``h(r)``, the quadrature the DST-IV pair evaluates at its nodes, with ``h(r) = 0`` beyond the grid.
   The solution's own nodes are copied; other wave numbers cost one pass over the r grid each.

``susceptibility_oz_at(k, chi, error, charge_response)``
   Reconstruct ``chi`` on ``[0, pi/dr]`` from cubic interpolation of the native transform of
   the short-range direct correlation, analytical Coulomb tails and analytical rigid-molecule
   form factors. The optional ``charge_response`` is ``chi q / k^2``, ``(nk, nsite)``, including
   its finite zero limit; it is evaluated directly rather than by subtracting matrix entries.
   This is the default transfer to the 3D model and the DRISM transfer to 1D UV.
   Ordinary RISM retains the existing 1D UV transfer; the native stored ``chi`` remains
   ``w + rho h`` for both theories.

``residual()`` and ``solution_generation()``
   Achieved RMS residual and identity of the accepted bulk solve. Copies retain the identity;
   every successful solve takes a new one, including a solve with unchanged potential but
   different grid, closure or convergence settings. When the solution is stale they return
   ``huge(1.0_wp)`` and ``-1``, respectively.

``thermodynamics(thermo, error)``
   ``solvent_vv_thermodynamics_type`` in atomic units:

   * ``chemical_potential``: excess chemical potential per molecule,
     ``mu = (rho/beta) sum_ab int Phi_ab(h, c) d^3r`` with the closure's local energy integrand.
   * ``direct_moment``: ``C = sum_ab c_ab(k = 0)``, over all ordered site pairs.
   * ``compressibility``: ``kappa_T`` from ``rho kT kappa_T = 1/(1 - rho C)``.
   * ``pressure``: ``P = (n + 1)/2 rho kT - kT/2 rho^2 C`` for ``n`` sites
     (Sergiievskyi, Jeanmairet, Levesque and Borgis, J. Chem. Phys. 143, 184116 (2015), Eq. 21).
     The expression is derived for HNC and is evaluated unchanged for KH and PSE; it is the
     RISM bulk pressure used in pressure corrections, not the physical pressure of the fluid,
     which XRISM overestimates by orders of magnitude for water.
   * ``pressure_ideal``: ``P_ideal = rho kT``, the ideal-gas pressure of the PC+ correction
     (:ref:`moz-3d-pressure-corrections`).

Validation
----------

The historical table-geometry regression uses water (solvent 175, 298.15 K)
with the SPC/E charges and the ``lj_cspce`` parameters,
SPC/E Lennard-Jones on oxygen and the hydrogen Lennard-Jones modification of Luchko et al.,
with KH on the default grid reproduces the pysofft VV solve (``moist_dev`` ``origin/pysofft``,
d3c3d076, GMRES without preconditioner): ``rho mu`` (pysofft's ``DG_vv``) agrees to 2e-11
relative and ``h(r)`` to 4e-14 at every point. This checks the numerics against pysofft on
the same cSPC/E nonbonded parameters on the table geometry.
The ``moz_solvent_vv`` test ``pysofft_reference`` keeps these numbers.

RISM and DRISM use the same closures, solvers and fixed-bulk UV derivatives.
Theory, dielectric, smearing width, closure, radial grid and achieved residual
are checked before reusing a solved bulk solvent. A model supplied with a solved
solvent and no overriding parameters inherits that solvent's bulk settings.
