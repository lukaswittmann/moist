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

Keep the native UV settings at the root.
The input adds ``grid`` and ``potential`` groups, ``solvent.potential`` for the bulk interaction parameters and ``checkpoint``.
The ``solvent`` group holds the bulk solver settings and the solvent itself: ``name`` from the solvent table, ``temperature`` in K and ``pressure`` in Pa, the last two defaulting to the solvent table.
The total charge is not an input setting; it comes from the structure (see below).

For a one-atom chloride XYZ, run with ``--charge -1``, a complete 3D input is:

.. code-block:: toml

   closure = "KH"
   correction = "PC+"
   finite_size = "huenenberger"

   [potential]
   lj = "uff"
   charges = "custom"
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
   name = "water"
   temperature = 298.15  # K; optional
   theory = "RISM" # default; "DRISM" uses solvent.drism settings
   closure = "KH"

   [solvent.convergence]
   tolerance = 1e-11

   [solvent.potential]
   lj = ["cspce", "gaff2", "uff", "tm"]
   charges = ["spce", "resp-solvent"]
   mixing = "lorentz-berthelot"

``q`` supplies atomic partial charges of the solute; the structure supplies their total.
``--charge`` changes the total only. It does not generate atomic charges.

The bulk defaults are the lists of the example.
Water takes ``cspce`` LJ with cSPC/E charges; the cSPC/E choices leave any other solvent to GAFF2, UFF and TM LJ and to
RESP charges of the solvent table.
``lj = "spce"`` selects standard SPC/E with zero hydrogen LJ; ``lj = "cspce"`` selects the hydrogen LJ modification. 
Both use ``charges = "spce"``.
When the first bulk LJ table is either SPC/E preset, water gets 1 Angstrom O-H bonds and a 109.47 degree H-O-H angle. 
Other bulk potentials retain the general solvent-table geometry, and the solute keeps its input geometry.

VV convergence depends on the chosen water potential and solver settings.
Standard ``spce`` has no hydrogen LJ regularization and can be harder to converge; ``cspce`` remains the default. 
Selecting ``spce`` keeps its parameters unmodified even when a solve fails to converge.

Older ``spce`` inputs selected the modified LJ parameters. Change their LJ choice to ``cspce`` to retain those nonbonded parameters.
CLI water geometry now follows the selected named water model; historical table-geometry results require an explicitly supplied library structure.
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

``lj`` and ``charges`` take one choice or a list in order of precedence. Each atom
takes its parameters from the first listed choice that covers it, so
``lj = ["gaff2", "uff"]`` uses GAFF2 types and UFF for any atom GAFF2 leaves. An
atom that no listed choice covers is an error; so are an empty list and a choice
listed twice.

``lj`` accepts ``uff``, ``dreiding``, ``tm``, ``gaff2``, ``oplsaa``, ``spce``,
``cspce`` or ``custom``. ``custom`` covers every atom and stands alone; supply
``sigma`` in bohr and ``epsilon`` in Hartree, each with one value per atom.
OPLS-AA uses the structure's bond graph and requires
``solvent.potential.mixing = "geometric"``.

``charges`` names the sources of the site charges:

- ``custom``: ``q``, one charge per atom
- ``none``: zero charges
- ``spce``: the SPC/E water charges; other structures are left to later choices
- ``ec``: charges fitted to the host potential, 3D solute only; needs a host that
  supplies the potential, such as a quantum-chemistry integration
- ``<model>-<environment>``: solvent table charges, bulk only; models ``hirshfeld``,
  ``resp``, ``mbis`` and ``chelpg``, environments ``gas``, ``solvent`` and
  ``conductor``, for example ``"resp-solvent"``; solvents missing from the table
  are left to later choices

``custom``, ``none`` and ``ec`` cover every atom and stand alone. With ``charges = "ec"``,
``potential.ec`` takes the fit settings ``exclusion_radius`` and ``range_scale``
(bohr), ``long_range_bias``, ``restraint_strength`` and ``charge_penalty``; their
defaults are those of ``ec_charges_type``.

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

Checkpoints and warm starts
---------------------------

Moist caches the bulk VV and UV fields in ``.mbin`` files by default.
A later run starts from them: a field that already solves the current equations to the requested tolerance is accepted without iterations, any other fitting field is an initial guess.
Convergence is decided by that residual check, never by the file metadata; see
:ref:`moz-checkpoints`.

.. code-block:: toml

   [checkpoint]
   enabled = true              # false: read and write nothing

   [checkpoint.vv]
   directory = "/data/rism-vv" # directory of the automatic name; default: working directory
   interval = 5                # also save every 5 iterations; default 0: final field only

   [checkpoint.uv]
   file = "chloride.mbin"      # default: uv_rism_1d.mbin or uv_rism_3d.mbin
   guess = "chloride-hnc.mbin" # optional: start from another file

- The automatic VV name is ``vv_rism_solvent_<id>_<theory>_<closure>.mbin`` in ``checkpoint.vv.directory``. When the file leaves that unset, the CLI takes the directory from ``MOIST_RISM_VV_OUTDIR``, so runs can share the bulk solution.
- Without ``file`` and ``guess``, a VV solve with no file for its closure starts from the automatic file of another closure of the same theory (HNC, PSE4, PSE3, PSE2, KH, in that order).
- ``guess`` replaces the file read at the start; ``file`` is still written.
- Intermediate saves are marked unconverged. The final field is saved when the solve iterated and ended without an error, so an interrupted or unconverged run can be continued.
- Fields of another grid shape or spacing, or of other solvent or solute atoms, are ignored with a warning and the solve starts cold, then overwrites the file.

Overrides and output
--------------------

Values come from compiled defaults, then the parameter file, then the command
line. The solvent comes from ``--solvent``, else ``solvent.name``.
``--temperature`` and ``--pressure`` override
``solvent.temperature`` and ``solvent.pressure``, which default to the solvent
table. The total charge comes from ``--charge``, else ``.CHRG``, else the
structure. UV and bulk closure, solver and convergence controls stay in the input
file. The CLI refuses ``charges = "ec"``, since it has no host potential to fit.
The shared run context supplies solver verbosity and thread count; use the
CLI verbosity options and ``--threads`` for these controls.
The CLI defaults to verbosity 4. Verbosity 2 prints the main input settings, the
potential tables and, as the bulk solve runs once, the 1D VV iteration table and
result block. 3 adds every setting and the UV result block, 4 the UV iteration
tables, and 5 the inner Krylov rows and preconditioner/deflation diagnostics.
Solver and initial-guess details appear above the table headings.

``--grad`` evaluates the native analytical nuclear gradient with fixed atomic
charges and bulk settings. ``--numgrad`` uses central differences with a
``1e-4`` bohr displacement and restores the original model geometry.
Energy is reported in Hartree and kcal/mol; gradients are in Hartree/bohr.
``--writeenergy`` writes the Hartree value to ``.GSOLV``.
At verbosity 2 the solvent system is printed, followed by the resolved input
table: the grid and potential settings, the closure, the convergence controls and
the selected solver with its main controls.
``-v`` lists every setting. Settings of solvers that are not selected never
print; they are still read from the input file.

Library use
-----------

The input file is the library type ``model_moz_1d_input_type`` or
``model_moz_3d_input_type`` (``src/moist/model/moz/input.f90`` and
``{1d,3d}/input.f90``). A host reads it and builds the model with the input
overload of the model constructor; the structure enters at ``update``:

.. code-block:: fortran

   type(model_moz_3d_input_type) :: input
   type(model_moz_3d_type) :: model

   call input%read_file("rism.toml", error, strict=.true.)
   if (.not. allocated(error)) call new_moz_3d_model(model, ctx, input, error)
   if (.not. allocated(error)) call model%update(mol, error)

This form takes the solvent from ``solvent.name``, ``solvent.temperature`` and
``solvent.pressure``; ``input%new_system`` returns that solvent system. A host
that resolves the solvent itself passes it instead,
``new_moz_3d_model(model, ctx, input, system, error)``, and the ``solvent`` name and
state point of the file are then unused. To override single values, write them to
``input%solvent_state`` before construction, as the CLI does.

The constructor validates a copy of the input and leaves the input unchanged. It
builds the bulk solvent, the grid from ``src/moist/math/grid/3d/factory.f90``, the
potential terms from ``src/moist/model/moz/potential/factory.f90`` and the
checkpoint stores. ``input%new_model`` returns the same model as the generic
solvation model type, with or without a solvent system. The constructors taking a
grid and an updated solvent stay available for hosts that build these themselves.
The CLI evaluates the model through the shared coupling lifecycle.
