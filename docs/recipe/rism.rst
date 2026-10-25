RISM from the command line
==========================

The CLI runs solute-solvent RISM in either dimension:

.. code-block:: console

   moist model rism1d solute.xyz -p run1d.toml
   moist model rism3d solute.xyz -p run3d.toml

The model prepares its bulk solvent internally through ``solvent_vv_type``.
``rism1d`` is a 1D **solute-solvent** calculation. Bulk VV has no model command.
TOML is the documented input format; the same registered fields also accept JSON.

Input structure and settings
----------------------------

Keep the native UV settings at the root. The app adds ``system``, ``grid`` and
``potential`` groups, and ``solvent.potential`` for the bulk interaction parameters.
The existing ``solvent`` solver group stays independent of the UV settings.

For a one-atom chloride XYZ, a complete 3D input is:

.. code-block:: toml

   closure = "KH"
   correction = "PC+"
   finite_size = "hunenberger"

   [system]
   solvent = "water"
   charge = -1.0
   temperature = 298.15  # K
   pressure = 101325.0  # Pa

   [potential]
   lj = "uff"
   charges = "fixed"
   q = [-1.0]          # electrons, one per atom in XYZ order

   [grid]
   kind = "cartesian"
   nx = 64
   ny = 64
   nz = 64
   dr = 0.6            # bohr
   margin = 10.0       # bohr

   [convergence]
   tolerance = 1e-7

   [solver]
   name = "JFNK"

   [solvent]
   theory = "RISM"  # default; "DRISM" uses solvent.drism settings
   closure = "KH"

   [solvent.convergence]
   tolerance = 1e-11

   [solvent.potential]
   lj = "cspce"
   charges = "cspce"
   mixing = "lorentz-berthelot"

``q`` supplies atomic partial charges; ``system.charge`` supplies their total.
The app checks that they agree within ``1e-6`` electrons. ``--charge`` changes
the total only. It does not generate atomic charges.

The bulk water default explicitly selects ``cspce`` for LJ and charges.
``spce`` selects standard SPC/E with zero hydrogen LJ; ``cspce`` selects
the hydrogen LJ modification. Both use the same charges. For water with either
named bulk LJ preset, the CLI uses 1 Angstrom O-H bonds and a 109.47 degree
H-O-H angle. Other bulk potentials retain the general solvent-table geometry,
and the solute keeps its input geometry. Select the bulk potential explicitly
for other solvents.

VV convergence depends on the chosen water potential and solver settings.
Standard ``spce`` has no hydrogen LJ regularization and can be harder to
converge; ``cspce`` remains the default. Selecting ``spce`` keeps its parameters
unmodified even when a solve fails to converge.

Older ``spce`` inputs selected the modified LJ parameters. Change their LJ
choice to ``cspce`` to retain those nonbonded parameters. CLI water geometry
now follows the selected named water model; historical table-geometry results
require an explicitly supplied library structure.
The VV residual tolerance defaults to ``1e-11``; an explicit file value is kept.

For 1D, use a radial UV grid; its settings are separate from ``solvent.grid``:

.. code-block:: toml

   [grid]
   kind = "radial"
   npts = 8192
   dr = 0.015

Pressure and electrostatic finite-size corrections are 3D settings. A 1D input
containing ``correction`` or ``finite_size`` is rejected. The leading finite-size
term requires a Cartesian grid, Ng splitting and the default OZ susceptibility.

The 3D response sampling defaults to ``susceptibility = "oz"``. ``"linear"``
interpolates the response between the VV grid's k nodes instead. See :ref:`moz-3d-susceptibility-convergence` for
the chloride refinement results and the conditions needed for energy comparisons.

Potential and grid choices
--------------------------

``potential.lj`` accepts ``uff``, ``dreiding``, ``tm``, ``gaff``, ``oplsaa``,
``spce``, ``cspce`` or ``custom``. ``lj_fallback`` optionally supplies another native LJ
term for atoms the first term did not cover. OPLS-AA uses the structure's bond
graph and requires ``solvent.potential.mixing = "geometric"``.

For custom LJ, supply ``sigma`` in bohr and ``epsilon`` in Hartree, each with
one value per atom. Solute charges accept ``fixed`` with ``q``, ``spce`` or
``cspce`` for water, or ``none`` for an explicitly uncharged interaction. Bulk charges also
accept the table models ``hirshfeld``, ``resp``, ``mbis`` and ``chelpg``, with
``environment = "gas"``, ``"solvent"`` or ``"conductor"``.

Both potential groups accept ``ng_split`` (default true) and ``alpha``
(default ``1.0`` in inverse bohr). The bulk mixing choice also applies to UV
tables. Missing LJ/charge data and unused explicit arrays are errors.

For an atom-centred 3D point grid, replace the Cartesian grid group with:

.. code-block:: toml

   [grid]
   kind = "molecular"
   nrad = 50
   angular_degree = 17
   rmax = 10.0          # bohr; bounded HandyMod radial mapping
   handymod_exponent = 2.0  # dimensionless; finite and positive
   becke_k = 3
   dr = 0.6            # bohr; controls the reciprocal lattice
   kbuffer = 2.0       # bohr
   nufft_tolerance = 1e-10

The bounded HandyMod mapping uses radial midpoint nodes and accepts a real
``handymod_exponent``, including values below one. A larger radial and angular recipe uses ``nrad = 50``,
``angular_degree = 77`` (2030 Lebedev points per shell), ``rmax = 16.0``
and ``handymod_exponent = 0.34``. Set the reciprocal spacing ``dr`` and NUFFT
tolerance separately. Explicit input values override the compiled defaults.
The mapping requires ``rmax > 2**handymod_exponent - 1`` in bohr; its constructor
rejects incompatible radial bounds.

Unknown keys and keys for another grid kind are rejected. Native closure and
solver validation determines which numerical options are available.

Overrides and output
--------------------

Values come from compiled defaults, then the parameter file, then explicit CLI
overrides. ``--solvent``, ``--charge``, ``--temperature`` and ``--pressure`` are
the system overrides. UV and bulk closure, solver and convergence controls
stay in the input file.
The shared run context supplies solver verbosity and thread count; use the
CLI verbosity options and ``--threads`` for these controls.
The CLI defaults to verbosity 4. Verbosity 2 prints the main input settings, the
potential tables and, as the bulk solve runs once, the 1D VV iteration table and
result block. 3 adds every setting and the UV result block, 4 the UV iteration
tables, and 5 the inner Krylov rows and preconditioner/deflation diagnostics.
Solver and initial-guess details appear above the table headings.

When the file omits them, solvent and charge retain the structure-directory
``.SOLVENT`` and ``.CHRG`` fallbacks. CLI charge overrides the file charge, which
overrides ``.CHRG`` or the structure charge.

``--grad`` evaluates the native analytical nuclear gradient with fixed atomic
charges and bulk settings. ``--numgrad`` uses central differences with a
``1e-4`` bohr displacement and restores the original model geometry.
Energy is reported in Hartree and kcal/mol; gradients are in Hartree/bohr.
``--writeenergy`` writes the Hartree value to ``.GSOLV``.
The resolved input table lists the system, grid and potential settings, the
closure, the convergence controls and the selected solver with its main controls.
``-v`` lists every setting. Settings of solvers that are not selected never
print; they are still read from the input file.

The standalone app assembles parameter-only terms. EC and other host-fed
electrostatics use the library's named coupling requests and need a host
provider, such as a quantum-chemistry integration.

``app/rism/input.f90`` owns run settings and reuses the parameter abstraction.
``app/rism/setup.f90`` assembles a polymorphic solvation model from the library
factories: ``src/moist/math/grid/3d/factory.f90`` builds the UV point grid and
``src/moist/model/moz/potential/factory.f90`` the native potential terms. The
driver evaluates that model through the shared coupling lifecycle; no filenames
enter the numerical model constructors.
