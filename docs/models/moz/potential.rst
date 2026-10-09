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

A host potential (``host_potential_type``) supplies the electrostatic field of the whole solute (e.g. in case of EC-RISM or RISM-cSED).
Host charges (``host_charges_type``) work in the 1D and 3D models alike.
The host potential is a field on the volume grid; it has no radial counterpart yet, so the 1D tables refuse it.

Construction
------------

.. code-block:: fortran

   use moist_data_solvents, only : solvation_system_type, new_solvation_system, get_solvent_id
   use moist_model_moz, only : solvent_vv_type, new_vv_solvent, model_moz_3d_type, new_moz_3d_model
   use moist_model_moz_potential, only : lj_spce, lj_gaff, lj_uff, lj_tm, coulomb_spce, &
      & coulomb_resp_solvent, host_charges_type

   type(solvation_system_type) :: system
   type(solvent_vv_type) :: solvent
   type(model_moz_3d_type) :: model
   type(host_charges_type) :: host
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
   call model%potential%add(host, error)

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

   ! Solute-solvent, on the model's volume grid
   call model%potential%compute(model%grid, u_sr, ur_lr, uk_lr, error, &
      & solvent=solvent%potential, coupling=coupling)

Without ``solvent`` the tables are solvent-solvent, with it solute-solvent. 
Terms that read host data read it through ``coupling`` each time; ``compute_adjoint`` hands the derivatives back to them as response items.

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
   * - ``host_charges_type``, ``host_potential_type``
     - Every solute atom, from the host (solute only).
   * - ``new_custom_lj``, ``new_fixed_charges``
     - Explicit per-atom values, optionally for a subset of the atoms.
