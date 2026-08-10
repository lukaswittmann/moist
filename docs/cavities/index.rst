Cavities
========

This section describes the cavity constructions available in MOIST and the settings that control their shape and surface discretization.

Common interface
----------------

All cavity implementations are based on the the abstract ``cavity_type`` defined in ``src/moist/type.f90``.
It exposes shared allocatable quantities that concrete cavities populate where applicable:

- atomic sphere centers, radii, and per-sphere areas;
- grid-point positions, owners, areas, Gaussian widths, and switching function;
- total cavity area and volume;
- optional nuclear derivatives of positions, widths, and switching function.

Concrete cavity types must implement two deferred procedures:

``update(mol, error)``
   Build or refresh the cavity for the supplied molecular structure and fill
   the common fields supported by that implementation.

``get_gradient()``
   Prepare or validate the implementation-specific nuclear derivatives when
   they are needed.

The base type also provides optional response hooks. 
``get_surface_potential`` maps model surface weights to host-potential contributions.

Typical use
-----------

1. Construct a concrete cavity and configure its radii and discretization.
2. Call ``update`` whenever its inputs change, for example after a geometry or
   density update.
3. A ``general_solvation_model`` owns one authoritative cavity, updates it
   before its components, and passes the same cavity to each component.
4. For the potential, components accumulate model surface weights and
   ``get_surface_potential`` maps them to host-response channels where needed.
5. Call ``get_gradient`` before requesting supported nuclear derivatives.

The shared ``print``, ``write_xyz_debug``, ``write_csv_debug``, and ``write_pqr_debug`` procedures provide diagnostics and grid export; ``find_disconnected_cavities`` checks the grid for disconnected cavities.

Available cavities
------------------

.. toctree::

   drop
   iswig
   marchingcubes
