MOZ Potentials
==============

The 1D and 3D MOZ models describe solvent and solute as sites, one per atom.
Each side has a potential (``moz_potential_type``): an ordered list of terms for Lennard-Jones parameters and electrostatics, plus the settings of the interactions it enters. The bulk solvent owns ``solvent%potential``; each MOZ
model owns ``model%potential`` and a copy of the updated solvent.

Idea
----

Terms are added in order.
A term that cannot parametrise an atom leaves the atom uncovered; later terms fill only the atoms still uncovered:

* **Lennard-Jones**: the first term with parameters for an atom supplies them.
* **Charges**: the first term with a charge for an atom supplies it (charges of several terms are not summed).
* The update fails, naming the atoms, only when an atom is left without Lennard-Jones parameters, or without a charge once any term supplies charges.

The order of ``add`` calls is the order of preference: put the specific force field first and the broad fallback last.
With verbosity 2 the update prints a table of every site with its charge, sigma, epsilon and the term that supplied each.

``monopole_type`` uses host-supplied atomic charges to build a point-charge potential in the 1D and 3D models.
``multipole_type`` uses host-supplied atomic multipoles, respectively. For 1D models; this averages out so that only the monopole term actually contributes.
``ec_charges_type`` uses the full host potential sampled on the 3D volume grid and fits its own charges for the Ng tail.
Each concrete electrostatic term owns its host-data requirements and adjoints.

Construction
------------

.. code-block:: fortran

   use moist_data_solvents, only : solvation_system_type, new_solvation_system, get_solvent_id
   use moist_model_moz, only : solvent_vv_type, new_vv_solvent, model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_potential, only : lj_spce, lj_gaff, lj_uff, lj_tm, coulomb_spce, &
      & coulomb_resp_solvent, monopole_type

   type(solvation_system_type) :: system
   type(solvent_vv_type) :: solvent
   type(model_moz_3d_type) :: model
   type(monopole_type) :: monopoles
   integer :: id

   ! Bulk solvent from the solvent table: id, temperature, density, structure
   call get_solvent_id("water", id, error)
   call new_solvation_system(system, id, error=error)
   call new_vv_solvent(solvent, ctx, system, error)

   ! Water: SPC/E covers every site. Any other solvent: SPC/E covers nothing,
   ! GAFF and the RESP table charges fill in
   call solvent%potential%add(lj_spce, error)
   call solvent%potential%add(lj_gaff, error)
   call solvent%potential%add(coulomb_spce, error)
   call solvent%potential%add(coulomb_resp_solvent, error)
   call solvent%update(error)

   ! Solute: GAFF where it has a type, the transition-metal table for its
   ! metals, UFF (H through Lr) for anything left; charges from the host
   call new_moz_3d_model(model, ctx, grid, solvent, error)
   call model%potential%add(lj_gaff, error)
   call model%potential%add(lj_tm, error)
   call model%potential%add(lj_uff, error)
   call model%potential%add(monopoles, error)

   call model%update(mol, error)

The model keeps its own copy of the updated solvent; terms added to ``solvent`` afterwards do not reach it.

Settings
--------

Each potential carries the Lennard-Jones mixing rule and the Ng split of the Coulomb tail: ``%set_mixing(lj_mixing_geometric, error)`` and ``%set_ng_split(ng_split, alpha, error)``.
The defaults are Lorentz-Berthelot, split on and alpha = 1/bohr. The potential the tables are computed on decides: the solvent's settings govern the solvent-solvent tables, the model's the solute-solvent tables.

Tables
------

``compute`` gives the Ng-split tables ``u_sr``, ``ur_lr`` and ``uk_lr``:

.. code-block:: fortran

   ! Solvent-solvent, on radial grids; columns from solvent%potential%site_pairs
   call solvent%potential%compute(rgrid, kgrid, u_sr, ur_lr, uk_lr, error)

   ! Solute-solvent, on the model's volume grid (3D)
   call model%potential%compute(model%grid, u_sr, ur_lr, uk_lr, error, &
      & solvent=solvent%potential, coupling=coupling)

   ! Solute-solvent, on the model's radial grids (1D); columns from model%potential%site_pairs
   call model%potential%compute(model%rgrid, model%kgrid, u_sr, ur_lr, uk_lr, error, &
      & solvent=solvent%potential, coupling=coupling)

Without ``solvent`` the tables are solvent-solvent, with it solute-solvent.
The 3D model owns its volume grid, the 1D model its radial grid pair (``new_moz_1d_model(model, ctx, rgrid, kgrid, solvent, error)``, e.g. from ``new_uniform_radial_pair``); the host requests are declared on these grids.
The model updates its potential on its grid (``potential%update(mol, grid, error)``), so terms with grid-dependent data, such as the EC fit, build it once per geometry; the solvent potential is updated without a grid.
The 3D model takes point grids only: Gaussian-width grids are refused, since the widths give the tables nothing and cost more.
Terms that read host data read it through ``coupling`` each time; ``compute_adjoint`` hands the derivatives back to them as response items.

Each term owns its coupling requirements, reads and input adjoints.
The potential assigns a coupling scope to each term and dispatches its ``declare_pass`` method; the coupling layer defines and validates the protocol.
``monopole_type`` requests charges in every phase and returns ``atomic_charge_adjoint`` (``dg_dq``).
``ec_charges_type`` requests the point potential in every phase and its point derivative in the gradient phase.
It returns ``potential_adjoint`` (``w_phi``: the tables' field adjoint plus the part through the fitted charges); with a gradient request it also adds grid-point and volume-weight adjoints and the explicit nuclear gradient of the fit.
The monopole, multipole and EC implementations live in ``monopole.f90``, ``multipole.f90`` and ``ec.f90``, respectively; shared request types live in the coupling layer. The multipole term is currently a stub.

EC tail-charge fitting
----------------------

``ec_charges_type`` uses the full host potential on a volume grid and fits atom-centered charges for the Ng tail. It requests only ``phi`` in the energy and response phases, and also ``dphi_dr`` in the gradient phase.
No host atomic charges are needed, and no reference charges are supplied by default.
The radial (1D) model refuses this term.

The fit minimizes

.. math::

   \frac12 \sum_g \bar w_g
   \left(\sum_a \frac{q_a}{|r_g-R_a|} - \phi_g\right)^2
   + \frac{\lambda}{2}\sum_a(q_a-q_a^{\mathrm{ref}})^2
   + \kappa\sum_a\max(0, |q_a|-t_a),
   \qquad \sum_a q_a = Q.

``Q`` is the structure's total charge. The exact charge constraint preserves the monopole of the long-range tail.
The fitted charges enter the Ng split; the full field remains the host's potential.

The normalized fitting weights are proportional to the volume quadrature weights, a smooth nuclear exclusion mask, and ``1 + long_range_bias * R**2 / (R**2 + range_scale**2)``, where ``R`` is the distance from the nuclear centroid.
The default ``long_range_bias = 1`` and ``range_scale = 5`` bohr increase the relative weight of distant samples up to a factor of two.
Points within ``exclusion_radius = 1`` bohr of any nucleus are excluded; a quintic switch reaches full weight at twice that radius, with continuous first and second derivatives.

``restraint_strength`` is ``lambda``, default 0.01/bohr**2, and must be positive.
``reference_charges`` can supply one reference charge per atom; when absent, the neutral-atom reference is shifted uniformly to ``Q / natom``.
``charge_penalty`` is ``kappa``, default 0.1 e/bohr**2. The default threshold ``t_a`` is ``0.5 * Z_a`` e, raised to ``Z_a`` e for H and Li. ``charge_thresholds`` can override these thresholds per atom.
The linear penalty discourages excessive positive and negative charges.
Optional ``lower_bounds`` and ``upper_bounds`` enforce hard per-atom bounds exactly; they must permit the structure's total charge.

An active-set algorithm solves the convex, piecewise quadratic fit.
Within each active set, a direct Cholesky solve and elimination of the charge multiplier enforce the total-charge constraint. Active penalty thresholds and bounds are released when their optimality conditions fail.

Configure the term before adding it to a potential, for example:

.. code-block:: fortran

   type(ec_charges_type) :: ec

   ec%long_range_bias = 2.0_wp
   ec%restraint_strength = 0.02_wp
   ec%charge_penalty = 0.2_wp
   ! Optional: ec%reference_charges = reference
   ! Optional: ec%lower_bounds = lower; ec%upper_bounds = upper
   call model%potential%add(ec, error)

The fit is built once per geometry, when the model updates its potential on its volume grid (``model%update``, which calls ``potential%update(mol, grid, error)``):
the exclusion mask, the weights and the restrained Hessian are kept; the point-charge basis :math:`1/|r_g-R_a|` is recomputed in grid blocks instead of stored.
Declarations, fields and adjoints on any other grid are refused.
Each evaluation refits charges to the current host potential.
Fit derivatives are included in the ``potential_adjoint`` response, the explicit nuclear gradient, and the volume
position and quadrature-weight adjoints.
EC emits no atomic-charge adjoint.
Derivatives use the current active constraints; the linear penalty and hard bounds give a piecewise differentiable fit at constraint changes.

Lifecycle
---------

``add`` drops the update of a potential; ``compute`` refuses until the next ``update``.
Neither the solvent nor the model watches its potential: after adding a term or changing a setting, update (and solve) them again.

In a self-consistent host loop such as EC-RISM, the host terms are added once.
Each iteration the host answers the coupling with its current charges or potential and ``compute`` reads them; no new ``add`` or ``update`` is needed until the structure changes.

Terms
-----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Term
     - Covers
   * - ``lj_gaff``
     - Atoms with a GAFF2 type, typed from element and geometry.
   * - ``lj_oplsaa``
     - Atoms matched by an OPLS-AA rule; typed from element and geometry (bonds, bond orders, ...).
   * - ``lj_uff``
     - Every element from H through Lr (element specific only); so a good "add it last".
   * - ``lj_dreiding``
     - Elements of the DREIDING table.
   * - ``lj_tm``
     - Fourteen transition metals.
   * - ``lj_spce``
     - The sites of neutral water.
   * - ``coulomb_<model>_<environment>``
     - Solvents with tabulated charges: ``hirshfeld``, ``resp``, ``mbis`` or
       ``chelpg`` from a ``gas``, ``solvent`` or ``conductor`` calculation.
   * - ``coulomb_spce``
     - The SPC/E charges of neutral water.
   * - ``monopole_type``
     - Host-supplied atomic charges and their point-charge potential (1D and 3D solute only).
   * - ``multipole_type``
     - Host point multipoles (stub).
   * - ``ec_charges_type``
     - Full host potential with constrained ESP-fitted Ng tail charges (3D solute only).
   * - ``new_custom_lj``, ``new_fixed_charges``
     - Explicit per-atom values, optionally for a subset of the atoms.
