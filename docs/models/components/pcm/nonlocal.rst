Nonlocal electrostatics
=======================

This model :cite:p:`hildebrandt2004nonlocal` solves a three-field nonlocal dielectric boundary system.
It uses the solvent permittivity ``eps_solvent``, the high-frequency permittivity ``eps_inf``, and the positive correlation length ``lambda`` in bohr.
Fortran optionally accepts the solute permittivity ``eps_solute``; C and Python take it explicitly.
The model requires cavity normals and principal curvatures; the curvature enters the self term of the double-layer operator.

The exterior Green identity of the solute potential, which holds for a solute whose charge lies inside the cavity, is removed from the right-hand side.
Every remaining source term carries a permittivity contrast, so a uniform medium (``eps_solvent = eps_inf = eps_solute``, vacuum included) gives exactly zero energy and response on any surface.
Solute charge outside the cavity therefore acts only through these contrast terms.

.. code-block:: fortran

   call new_component_nonlocal(nonlocal_pcm, eps_solvent=78.4_wp, &
      eps_inf=1.78_wp, lambda=10.0_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_nonlocal_component``.
Python uses ``ModelComponentNonlocal(eps_solvent, eps_inf, lambda_, eps_solute, parameters=...)``.
All numerical settings can be found in :doc:`index`.
