Nonlocal electrostatics
=======================

This model :cite:p:`hildebrandt2004nonlocal` solves a three-field nonlocal dielectric boundary system. Pass the
solvent permittivity ``eps_solvent``, high-frequency permittivity ``eps_inf``
and positive correlation length ``lambda`` in bohr. Fortran optionally accepts
the solute permittivity ``eps_solute``; C and Python take it explicitly. The
model requires cavity normals.

.. code-block:: fortran

   call new_component_nonlocal(nonlocal_pcm, eps_solvent=78.4_wp, &
      eps_inf=1.78_wp, lambda=10.0_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_nonlocal_component``. Python uses
``ModelComponentNonlocal(eps_solvent, eps_inf, lambda_, eps_solute,
parameters=...)``. All numerical controls come from :doc:`index`.
