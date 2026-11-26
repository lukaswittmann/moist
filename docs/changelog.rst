Release Notes
=============

2.0.0-alpha
-----------

Breaking changes
   - mctc-lib 0.7.0 or newer is required; the bundled subproject is updated to 0.7.0.
   - DROP cavity

     - Refactored surface projection scheme into batched seeding, refinement and per-anchor finalization; valid ``projection.level`` (``proj_level``) values are now ``1`` to ``4`` (default ``2``), while all other values are refused; the SLSQP, deflation and onion solvers (former levels 2 and 5-9) are removed, and level ``2`` is now the local projection with a search of the anchors that need one:

       - **Seeding stage**: Branch searches are seeded from a lattice of retracted level-set samples built once per geometry (levels ``1`` to ``4``: ``1.0``, ``0.5``, ``0.3`` and ``0.2`` bohr spacing; level ``4`` with one refinement level), which replaces the per-anchor octree. Seed filters drop samples that cannot lead to a weighted point: ``cone``, ``grazing``, ``back_facing`` and ``focal`` at level ``2``, the tangent cone alone at level ``3``, none at levels ``1`` and ``4``. Spacing, refinement and filters follow from the level; the keys ``projection.octree.seed_size``, ``projection.octree.max_boxes``, ``projection.octree.max_survivors`` and ``projection.octree.max_depth`` are removed. Level ``1`` searches only anchors whose local solve fails; level ``2`` searches anchors whose local point is not fully open or lies more than ``0.5`` bohr away, levels ``3`` and ``4`` every anchor. At every level, anchors whose local solve fails are retried nearest-first from the level's lattice. The search radius of an anchor adapts to the gates of its points: up to the branch window behind the nearest fully open point (gate exactly one), and up to the reach ``branching.reach_off`` while there is none; the seeder does not shrink its ball on surface whose gate it does not know. With a symmetry model the lattice follows the molecular frame, or an orientation the point group maps onto itself where the group fits a cube.
       - **Refinement**: Anchors and seeds are solved by a batched KKT Newton solver (the full 4x4 KKT system at every step, step-capped by ``projection.newton.max_step`` with exact-penalty merit backtracking), followed by constrained-optimality verification and Riemannian saddle escape. Seeded solves stop early as duplicates once an unscaled Newton prediction lands within ``1e-4`` bohr of an accepted point of the same anchor or of a strict minimum converged in the same solve; a repeated visit of a saddle skips its escape once the first escape from it has been accepted.
       - **Finalize**: Each anchor's accepted points are deduplicated within ``branch_sep_cut`` (now at least ``1e-6`` bohr), ordered canonically, evaluated for their closest-point Jacobian scaling and gate, and weighted. A branch below ``wleb_cut`` is taken out and the weights are formed again over the branches left, so the grid rows of an anchor are exactly the points its weights depend on. ``branch`` is the position of a row among the rows of its anchor, and anchors with no accepted projection are dropped with a warning.
       - **Output**: The projector reports its phases (local projection, seed lattice, branch search) as step lines with their wall time instead of a progress table. The three phases are separate entries of the timing table, and the gradient reports the same kind of step line. A step line reads ``Seed lattice ... done (0.01s)``. At verbosity ``4`` the seed lattice step is followed by the seeder's build tallies (samples, node probes, unresolved cells, retractions, covering radius, lattice memory in MB). The per-point polish and KKT traces and the lost-anchor lines move from verbosity ``5`` to ``6``. The DROP parameter block is regrouped: Lebedev points, main tolerance and projection level at verbosity ``2``; projection, seeds, branching, screening and the optional properties at ``3``; the switching windows at ``4``. The SvdW smoothing ``1b`` / ``2b`` / ``3b`` print from verbosity ``3``, the SvdW and CFC fade from ``4``. Logical values in pretty-printed blocks read ``yes`` / ``no`` instead of ``T`` / ``F`` and are right-aligned with the integers.

     - Reworked weights, switching and branching weights:

       - Branch weights in concave regions are now a pairwise partition of unity in the softened projection-distance differences instead of a softmax: a branch lying ``branching.dpdist_max`` (default ``0.2`` bohr) or more behind a fully open branch has exactly zero weight and zero derivative, at any distance from the anchor.
       - Each branch carries a gate ``w_f0 = f_crit * f_foc * f_hor * f_reach`` that scales its own weight and is the strength with which it takes weight from the farther branches of its anchor (:math:`m_{bc} = g_c\,[1 - h]`). A branch with a closed gate takes no weight from its siblings and changes none of their derivatives, so a projection born at a fold, at a level-set critical point or beyond a switch band enters continuously; behind a partly open branch a farther one keeps the factor :math:`1 - g_c`. No weight exceeds its gate and the weights of an anchor sum to at most one.
       - New reach gate ``f_reach`` (``branching.reach_on`` / ``branching.reach_off``, defaults ``3.78`` / ``5.67`` bohr, about 2 / 3 Å, for every level set) as a factor of the gate: a projection farther from its anchor than ``reach_off`` has no weight, which bounds the search of every anchor. The build warns when the gate removes more than 0.1 % of the area of the points it keeps.
       - New horizon switch ``f_hor`` on the signed :math:`\mathbf{n}_\mathrm{surf}\cdot\mathbf{n}_\mathrm{sph}` (``switching.w_0hor_off`` / ``switching.w_0hor_on``, defaults ``0`` / ``0.15``, ``w_0hor_on = 0`` disables it) as a factor of the gate: the gate of a projection closes exactly where its area scaling, proportional to the cosine, reaches zero, so a grazing projection leaves before its Gaussian width diverges and a projection onto the inner side of a wall (normal pointing into the anchor sphere) is switched off.
       - The focal switch ``f_foc`` now acts on the dimensionless focal margin, the smallest tangent eigenvalue over ``objective.alpha``, so the gate no longer depends on the objective weight; the defaults of ``switching.w_0tra_off`` / ``switching.w_0tra_on`` change from ``0.1`` / ``0.3`` to ``0.2`` / ``0.6``, which leaves the results at the default ``objective.alpha = 0.5`` unchanged. Explicit values of these edges must be divided by ``objective.alpha``.
       - Every grid row has an area of its own: a branch below the weight cutoff is no row and does not enter the weights of its siblings.
       - The off-edge of the critical-gradient switch ``f_crit`` (``switching.w_0ls_off``) defaults to ``0.1`` instead of ``0.25``: in scans of neck formation and breaking at ``blend_k`` 2 to 5.5 the area is smoother (largest slope 0.82x, largest second derivative 0.63x, medians) and about half of the regular seam area the old band removed is kept.
       - DROP switching transition bounds are renamed from ``*_from`` / ``*_to`` to ``*_off`` / ``*_on`` (``switching.w_0ls_off`` / ``switching.w_0ls_on``, ``switching.w_0tra_off`` / ``switching.w_0tra_on``, and Fortran fields ``wleb_prune_off`` / ``wleb_prune_on``).
       - The gate and the branch share are combined into ``wbranch`` (equal to ``w_f0`` for an unbranched point) and now multiply the switching value ``f = f_anchor * wbranch`` instead of the Lebedev weight ``wleb``. ``wleb`` remains purely geometric (setting the Gaussian width ``xi``); point areas are unchanged for unbranched points, but ``f``, ``xi`` and with them the PCM matrices of switched or branched points differ from earlier versions.
       - ``branching.softmax_scale`` and ``branching.weight_floor`` are replaced by ``branching.dpdist_max``; the C field and Python parameter ``branch_weight_s`` is renamed to ``branch_dpdist_max`` and is now the distance window in bohr (no longer a softmax scale), which also sets the branch-search window. The explicit window ``projection.branch_search_rho2_slack`` is removed: a narrower window drops branches that carry weight and a wider one only finds branches without weight, so files that set it are refused.

     - ``screening.cell_grid_fraction`` is removed; the candidate cell list now sizes its cells from the atom screening reach.
     - The default projection iteration limit (``projection.maxiter``, C ``proj_maxiter``, Python ``DROPParameters.proj_maxiter``) is lowered from ``150`` to ``50``; a local solve that needs more fails and is retried from the level's seeder.
     - Switch bands are validated: ``switching.w_0ls_off`` / ``w_0ls_on`` and ``switching.w_0tra_off`` / ``w_0tra_on`` must satisfy ``0 <= off < on``, the horizon edges ``0 <= off < on <= 1`` unless the on-edge is ``0``, and ``switching.w_0ls_p`` / ``w_0ls_a`` must be positive; reversed edges, which silently inverted a switch, and zero-width bands are refused.
     - The projection distance is named ``pdist`` instead of ``rho``, which denotes the electron density elsewhere: the parameter-file key ``branching.drho_max`` is now ``branching.dpdist_max``, the C field and Python parameter ``branch_drho_max`` is now ``branch_dpdist_max``, and the cavity field ``rho`` (also on Python ``CavitySnapshotDROP``) is now ``pdist``.
     - The optional derivatives of the projection distance and of the distance to the owner atom are removed: the switches ``do_rho`` and ``do_r_iI``, the ``rho1_rA`` and ``r_iI1_rA`` arguments of ``moist_get_cavity_gradient`` and the cavity field ``r_iI0`` (also on Python ``CavitySnapshotDROP``). A host forms both derivatives from ``xyz1_rA`` in one dot product per entry. The projection distance ``pdist`` itself is still stored.

     - Every level set carries a constant multiplier ``scale`` (:math:`S = \alpha\,S_\mathrm{LSF}`, a component of the abstract LSF type, copied from the family's parameter): CFC now defaults to ``0.25`` and SSES to ``0.088`` instead of ``1``, which puts :math:`\lVert\nabla S\rVert` near one on the surface, as for SvdW. The surfaces are unchanged, but CFC and SSES values and derivatives, and with them the critical-gradient switch and the projection tolerance in :math:`S`, see the scaled field. SvdW (``1``) and isodensity (``1000``, :math:`1/\rho_\mathrm{iso}`) keep their values.

   - API and output

     - A parameter file with a key that no setting reads is refused and the error names every such key; formerly an unknown, misspelt or renamed key was ignored and the default used.
     - ``converged`` is removed from ``moist_get_cavity_results``, ``CavitySnapshot`` and ``Cavity``; failed local projections are retried from the level's seeder and anchors with no accepted projection are dropped with a warning.
     - ``moist_new_context`` and Python ``Context`` take an additional ``warnings_to_stderr`` flag (default ``false``); warnings now print to the context's output stream by default and are copied to standard error only when enabled.
     - ``moist_banner`` and ``get_banner`` styles are now ``full``, ``tagline``, ``logo``, ``build`` and ``compact`` (replacing ``short`` and ``ascii``).
     - Fortran ``moist_cavity_surface_adjoint`` (``cavity_surface_adjoint_type``) is renamed to ``moist_cavity_adjoint`` (``cavity_adjoint_type``).
     - Fortran ``moist_math_smoothstep`` smoothstep types and constructors now name their transition edges ``x_off`` and ``x_on`` (and ``k`` on ``new_smoothstep_erf``) instead of ``from`` / ``to`` (``from_val`` / ``to_val`` / ``k_val``).
     - The ``-p`` / ``--parameters`` option of ``moist model``, which had no effect, is removed.


New features
   - Cavity symmetry: iSwiG, DROP and NumSA accept an optional symmetry model that orients the Lebedev grids by a molecular frame (results invariant under rotation, translation and atom renumbering) or by the point group of the molecule (detected with libmsym or given as a Schoenflies symbol), with analytic gradients. Command line: ``--sym``, ``--orientation`` and ``--sym-tolerance``. Opt-in: results without a symmetry model are unchanged.
   - Molecular frame: its discrete choices switch at the settings ``gap_threshold`` (relative eigenvalue gap), ``harmonic_threshold`` (relative size of an in-plane harmonic) and ``linear_threshold`` (mean square distance from the long axis, bohr\ :sup:`2`). A geometry near a threshold gives a warning, a *Near a switch* report line and the fields ``symmetry_near_gap``, ``symmetry_near_linear`` or ``symmetry_near_harmonic``; frames that tie exactly but differ in the elements give a warning, a *Tied frames* report line and the ``symmetry_tied_frames`` field.

   - DROP cavity

     - Experimental smoothed solvent-excluded surface level set (``SSES``, ``moist_cavity_drop_lsf_sses_type``), parameterized by ``probe_radius`` and ``smoothing_length`` and available in DROP, marching cubes and marching tetrahedra (``moist cavity {drop,mc,mt} sses``).
     - ``scale`` setting for SvdW and CFC: parameter-file key ``scale``, Fortran ``param%new(scale=...)``, C ``moist_svdw_options.scale`` / ``moist_cfc_options.scale`` and Python ``SvdWParameters.scale`` / ``CFCParameters.scale``.
     - Optional LSF-partitioned per-atom surface areas (``do_atomarea``, enabled by ``do_fine``), partitioning each area element among atoms by their level-set shares on SvdW, CFC and isodensity cavities.
     - Batched-projection parameter-file keys ``projection.batch_size`` and ``projection.newton.{max_step,merit_margin}``.
     - Exported DROP row-layout fields ``offset`` (``nmax + 1``) and ``anchor_of_row`` (``ngrid``), included on Python ``CavitySnapshotDROP``.
     - The reach and horizon bands (``branch_reach_on`` / ``branch_reach_off``, ``w_0hor_off`` / ``w_0hor_on``) are fields of the C ``moist_drop_options`` and of Python ``DROPParameters``; the Fortran ``param%new`` takes the reach band in its ``branch`` argument (``reach_on``, ``reach_off``) and the horizon band as ``w_0hor_off`` / ``w_0hor_on``.

   - Continuum models and PCM

     - ``point_potential`` option on PCM parameters (Fortran, C ``moist_pcm_options``, Python ``PCMParameters`` and ``PySCFHost`` / ``PySCFSolvation``) to couple Gaussian surfaces through bare Coulomb ``point_potential`` integrals while keeping the Gaussian surface operator; mixing incompatible PCM integral modes in one continuum model is rejected.
     - Per-point component result fields exposed after ``get_energy``: ``potential_adjoint`` and ``surface_charge`` on PCM components, and ``omega``, ``gt``, ``ftilde``, ``amplitude`` and ``point_energy`` on GOSTSHYP.


Improvements
   - DROP cavity

     - The level-set Hessian weight ``w_lsf2`` of the surface weight contraction now uses its symmetry for every level-set function; this way, off-diagonal entries are not doubled. The same is done for the density-Hessian response weight ``w_hess_rho`` of the isodensity cavity.

       - A host code now only needs to evaluate and contract the six unique entries.
       - The reverse derivative pass now evaluates 10 instead of 13 level-set directions per surface
         point.

     - Certified level-set box exclusion (``unreached``) prunes empty regions early in the lattice seeder and marching tetrahedra.

   - Output

     - PCM parameter printouts show human-readable aliases for solver and storage modes.
     - Exact ties in the DROP watershed diagnostic (verbosity 3 with ``do_fine``) keep point order, so the printed patch statistics of symmetric molecules are symmetric (benzene dimer patches 30.61 / 30.61 instead of 29.05 / 32.18 bohr²). Cavity rows are unchanged.

   - Internals

     - Sorting is consolidated in ``moist_math_sorter``: hand ports of Fortran stdlib v0.8.1 ``ord_sort`` (stable merge sort), ``sort`` (introsort), ``sort_index`` and ``sort_adjoint`` (stable argsorts) for default-integer and ``real(wp)`` keys, MIT-licensed (``src/moist/math/sorter/LICENSE``). They replace ``qsort``, ``heapsort`` and the hand-written sorts of the seeder, the isodensity tile order, the GTO shell order, the adjacency rows and the island diagnostic. ``counting_argsort`` takes the caller's counters and allocates nothing. Exactly tied adjacency neighbours keep their fill order.
     - ``cavity_type_drop%anchor_subset`` projects only the listed anchors, a partial surface for tests and diagnostics. The branching gradient tests check the total area and volume on the anchors where branch weights and gates act (the branched ones, at most 48, every eighth partly open one and 10 % of the rest), after confirming that each kept anchor projects as in the full build; the five tests take 72 s instead of 230 s.


Performance
   - DROP cavity

     - Much better OpenMP scaling of the projection and of the gradient and response routes (smaller projection chunks, per-thread workers share the isodensity basis, fewer serial steps).
     - The spatial neighbour list (used for the DROP grid density, the iSwiG atom screening and the MT width correction) is built 1.3-1.8x faster on one thread and runs threaded with the context's thread count.
     - Batched KKT Newton projection with duplicate-seed capture and batched, symmetry-packed DROP derivative contractions.
     - Level-set screening queries that keep more than half of their candidates sort them with the stable merge sort instead of a heap sort (1.5-1.8x faster queries on the C60 tetramer).
     
   - PCM kernels are 1.5-3.5x faster, using the error-function approximations of D. N. Laikov, *J. Chem. Phys.* **165**, 184109 (2026), `doi:10.1063/5.0351157 <https://doi.org/10.1063/5.0351157>`_.

     - Every PCM matrix element, solve and gradient loops over all pairs of surface points. With Gaussian surface charges the pair interaction is ``erf(z r)/r``. The screened kernels of LPB and nonlocal electrostatics also need ``exp``, ``erfc`` and ``erfcx``. These special functions dominated the cost of the pair loops, and the compiler's library calls for them branch on the argument range, so those loops could not be vectorized.
     - Laikov gives one closed-form expression per function for the whole argument range, built from polynomials, a division, a square root and an exponential. The generated kernels now evaluate ``erf`` (Eq. 10), ``erfc`` / ``erfcx`` (Eq. 18) and ``exp`` inline without branches. Close pairs, where ``erf`` is evaluated, get these cheaper scalar functions. The LPB and nonlocal pair loops, where ``exp`` is evaluated for every pair, now vectorize.
     - The kernel assembly loops are also split into a vectorized pass over the far pairs and a separate pass over the close pairs.
     - The inline functions agree with the library functions to a few units in the last place, so PCM results change only in the last digits.


Fixes
   - Nonlocal electrostatics (NLE) uses the solute potential in the Laplace right-hand-side terms so the local limit ``eps_solvent = eps_inf`` matches IEF-PCM on its surface operators and polarizes outlying charge consistently.
   - Properly released ``ducc0`` worker threads across FFT calls and fixed concurrent formatting in ``prettyprinter``.
   - Corrected literature provenance notes in the built-in atomic radii tables.
   - The forward nuclear gradient of the DROP cavity (``get_gradient``) no longer drops atoms the level set does not reach at a projected point: the owner of the anchor is always seeded, also on isodensity level sets, which gave zero point derivatives, and the switching neighbours of the anchor enter the area and volume gradients. Level sets whose reach covers the switching neighbours of every anchor, such as SvdW with its default fade band, were not affected.
   - The ``--wleb-switch`` help of the command line lists the switching levels ``5`` and ``6``, which the parameters accept.
   - The anchor-only DROP gradient (``moist_compute_anchor_gradient``, Python ``compute_anchor_gradient``) holds the level set fixed in space as documented: every point moves with its owner's anchor only. It used to keep the level set's own nuclear derivatives, which made it equal to the full gradient on SvdW, CFC, SSES and internal isodensity level sets, and solved one zero system per atom on callback level sets.
   - The DROP curvature diagnostic no longer aborts the build on a point with a zero principal curvature; such points are reported as flat.
   - The isodensity density callback (``moist_isodensity_lsf_callback``) fails the cavity build on a NaN or infinite density, gradient or requested higher derivative, naming the point, as it does on a nonzero status; such a value used to be handed on to the projection.
   - A projected point with a vanishing level-set gradient (an anchor on a critical point of the level set) no longer gets a 0/0 normal; its critical-gradient gate closes it. The curvature diagnostic gives a point whose gradient norm is at most ``1e-14`` zero curvatures and counts it as flat.
   - The canonical branch order of an anchor groups coordinates chained by gaps of at most ``1e-8`` bohr into one tie cluster, so row order and branch numbering no longer depend on the order in which the projections were found.
   - The verbosity-6 traces of the DROP polish stage are written under the context's output lock, so worker threads no longer interleave lines on the output unit.
   - NaN and infinite values of the DROP ``tolerance``, ``density.rho_grid_h`` and ``objective.alpha`` are refused; each refusal names its parameter.
   - Updated the PySCF ECP ``core_fill`` inner-sheet tests (``test_pyscf_ecp.py``) to ``proj_level=3`` with a widened ``branch_dpdist_max`` search window and an explicit ``core_fill=False`` reference check so both missing inner-sheet detection and broken fills fail the test suite. The search runs with the horizon switch off and a ``6`` / ``7`` bohr reach, since the signed horizon closes the unfilled inner sheet by design.


Deprecations


Documentation
   - New section on symmetry and grid orientation of the Lebedev-based cavities.

   - DROP cavity

     - Added level set figures to DROP overview and each LSF section.
     - SvdW LSF contains detailed description of how the fade/switching works; this is also mentioned in the CFC section now.
     - The SSES page is now listed in the table of contents.
     - Documented the batched Newton projection levels and their seeders, ``projection.*`` settings and ``branching.dpdist_max`` gate-weighted partition of unity.

   - Documented ``point_potential`` PCM coupling, ``warnings_to_stderr`` on ``Context``, and per-component named result fields.

1.2.1-alpha
-----------

Breaking changes
   - Solvent molecular data: IDs 181 to 187 (argon, bromine, hydrogen chloride, helium, neon,
     krypton, xenon), which have no descriptors, moved to 182 to 188 to make room for ammonia.
   - NLE component: arguments ``eps_solvent``, ``eps_inf`` and ``eps_solute`` renamed to
     ``epsilon``, ``epsilon_inf`` and ``epsilon_solute`` (Fortran, C, Python).

New features
   - Ammonia solvent (ID 181).
   - Solvent aliases ``propanonitrile`` and ``propanenitrile`` (propionitrile) and ``butanonitrile``
     (butyronitrile).

Fixes
   - Python: version lookup no longer fails (API reference build).
   - Phenol (ID 180): ``solvent`` charges and multipoles recomputed with its own permittivity
     (12.4).

Documentation
   - Restructured cavities, models and components with support tables and recommendations.
   - New Continuum Models, Pressure Models and Solvents sections.
   - Expanded PCM overview (model choice, physical inputs, numerical options) and GOSTSHYP derivation.
   - Updated contributing guidelines.
