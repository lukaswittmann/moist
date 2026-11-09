Thermodynamic items
===================

Bulk pressure and compressibility, 3D UV partial molar volume, PC, PC+ :cite:p:`sergiievskyi2015pressure` and the leading huenenberger
corrections :cite:p:`huenenberger1999ewald,kovalenko2000ions` share one abstract ``thermodynamic_quantity_type``. ``thermo_list_type`` owns copies of
these items and their scalar inputs. The bulk and 3D thermodynamics records and model methods,
and the 3D ``correction = "none"``, ``"PC"`` or ``"PC+"`` settings keep their interfaces and formulas.

The 1D UV model reports free energy and partial molar volume through
``model_moz_1d_thermodynamics_type`` (:doc:`uv1d`), using the volume kernel directly.

Forward evaluation and selection
--------------------------------

An item declares a name, a unit and the names of its dependencies. ``evaluate`` receives
their values in that declared order. The list resolves the dependencies, rejects missing
names or cycles, and evaluates each item once in forward dependency order. Registration
order need not be dependency order. A repeated ``evaluate`` reuses the cached values;
``set_input`` invalidates them. There is no global registry or saved numerical state.

The bulk graph derives pressure, ideal pressure and compressibility from molecular number
density, inverse thermal energy and the bulk direct-correlation moment. The 3D UV graph
receives that solved bulk state and derives partial molar volume and both pressure terms
from the UV direct-correlation moment. Grid reductions remain with the VV and 3D UV drivers.

Both pressure corrections are evaluated as diagnostics. Only selected corrections are
added to the uncorrected energy and seed reverse sensitivities. PC and PC+ belong to the
same exclusive group, so selecting both is an error. Changing selection reuses the forward
values and invalidates the reverse cache. ``select(name, .false., error)`` removes a term.

``huenenberger_correction_type`` belongs to the separate ``electrostatic_finite_size``
group. ``new_huenenberger_correction`` fixes its total charge, bulk screening and box
lengths, evaluating the analytical Coulomb box average once. It has no variable dependencies;
its reverse method is a no-op. The 3D model enables it through the independent
``finite_size`` setting (:ref:`moz-3d-finite-size`).

For example, the legacy constructor can also return the owned item list:

.. code-block:: fortran

   use moist_model_moz, only: model_moz_3d_thermodynamics_type
   use moist_model_moz_thermodynamics, only: thermo_list_type, new_moz_3d_thermodynamics
   type(model_moz_3d_thermodynamics_type) :: thermo
   type(thermo_list_type) :: items

   call new_moz_3d_thermodynamics(thermo, free_energy, moment, bulk, beta, error, &
      & items=items, correction="PC+")
   if (.not. allocated(error)) call items%energy(energy, error)

``n_items`` and ``get_name`` support iteration; ``get_value`` reads the forward results.
The 3D record still provides ``free_energy``, ``partial_molar_volume``, ``pc`` and ``pc_plus``
and ``finite_size``. Its checked ``corrected_energy`` subroutine includes the configured
finite-size term and the requested pressure term:

.. code-block:: fortran

   call thermo%corrected_energy("PC+", energy, error)

The selector accepts ``none``, ``PC`` and ``PC+`` without regard to case or surrounding
whitespace. An unknown selector returns an error and leaves ``energy`` unchanged.
The bulk constructor similarly accepts an
optional output item list without changing its existing required arguments.

Reverse accumulation and units
------------------------------

``accumulate_adjoint`` first seeds the uncorrected energy and selected corrections with
one. It then visits the items in reverse dependency order. Each item's method adds
``seed * d(value)/d(dependency)`` to its dependency sensitivities. Shared dependencies
receive every contribution before their reverse method runs. A fresh reverse pass clears
the old sensitivities; diagnostic-only corrections add no independent seed.

All values use atomic units: energies in Hartree, pressures in Hartree/bohr^3,
compressibility in bohr^3/Hartree, volumes and the bulk moment in bohr^3, and beta in
Hartree^-1. The UV moment is dimensionless. Reverse sensitivities are ordinary
derivatives of the energy in Hartree; item methods do not introduce beta scaling.
Base energies and corrections must declare ``Hartree`` units; the list rejects other
energy units rather than applying an implicit conversion.

At the 3D UV response boundary, ``accumulate_solvation_adjoint`` reads the beta of the
cached forward state and converts the accumulated ``dE/dD`` to ``beta dE/dD`` exactly
once. The UV driver combines this dimensionless moment term with the accumulated
uncorrected-energy sensitivity before one shared transposed RISM solve. Its quadrature
already carries ``1/beta``. The potential-table adjoints then follow the existing path.

The present host response holds the bulk solvent fixed. Scalar bulk sensitivities are
available from the list, but this refactor does not implement a bulk RISM response.
The 3D UV energy, host response and nuclear gradient use the same item list. The gradient also
uses its energy and moment sensitivities for the explicit quadrature-weight derivative.

Extending the list
------------------

To add a quantity or correction, extend ``thermodynamic_quantity_type`` and implement ``evaluate``
and ``accumulate_adjoint``. Set its name, unit and dependencies before ``add``. Corrections
set ``correction = .true.``; an optional ``exclusive_group`` prevents incompatible selected
terms. ``add`` copies the item, so lists and copied lists own their state independently.
Numerical forward and reverse mathematics belong in pure kernels; item methods connect
those kernels to the list. Names and groups are matched without case sensitivity.

Finite-cavity and multipole extensions of the electrostatic finite-size correction
are separate from the implemented leading net-charge term.
