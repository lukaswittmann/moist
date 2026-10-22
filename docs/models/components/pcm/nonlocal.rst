Nonlocal electrostatics
=======================

This model :cite:p:`hildebrandt2004nonlocal,hildebrandt2007nonlocal` solves a three-field nonlocal dielectric boundary system, similar to that in NESSie :cite:p:`kemmer2018nessie`.
It uses the solvent permittivity ``eps_solvent``, the high-frequency permittivity ``eps_inf``, and the positive correlation length ``lambda`` in bohr.
Fortran optionally accepts the solute permittivity ``eps_solute``; C and Python take it explicitly.
The host supplies the ordinary electrostatic potential and outward normal derivative of its solute at the surface points.
The energy is the physical reaction energy of a solute embedded in a cavity of relative permittivity ``eps_solute``; ``eps_solute = 1`` for quantum-mechanical densities.
For a charge ``Q`` at the center of a sphere of radius ``R`` the energy is :math:`\frac{Q^2}{2R}\left[\frac1{\varepsilon} - \frac1{\eta} + \left(\frac1{\varepsilon_\infty} - \frac1\varepsilon\right)\frac1{1+\nu}\right]` with :math:`\nu = \sqrt{\varepsilon/\varepsilon_\infty}\,R/\lambda` and :math:`\eta` = ``eps_solute``; in the local limit ``eps_solvent = eps_inf`` it is :math:`\frac{Q^2}{2R}(1/\varepsilon - 1/\eta)`.
The model requires cavity normals and principal curvatures; the curvature enters the self term of the double-layer operator.

The exterior Green identity of the solute potential, which holds for a solute whose charge lies inside the cavity, is removed from the right-hand side.
Every remaining source term carries a permittivity contrast, so a uniform medium (``eps_solvent = eps_inf = eps_solute``, vacuum included) gives exactly zero energy and response on any surface.
Solute charge outside the cavity therefore acts only through these contrast terms.
The system is solved for the exterior auxiliary trace minus the solute potential, so in the local limit ``eps_solvent = eps_inf`` the Yukawa operators cancel exactly and the energy does not depend on ``lambda``.

.. code-block:: fortran

   call new_component_nonlocal(nonlocal_pcm, eps_solvent=78.4_wp, &
      eps_inf=1.78_wp, lambda=10.0_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_nonlocal_component``.
Python uses ``ModelComponentNonlocal(eps_solvent, eps_inf, lambda_, eps_solute, parameters=...)``.
All numerical settings can be found in :doc:`index`.
