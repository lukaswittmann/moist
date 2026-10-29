Nonlocal Electrostatics (NLE)
=============================

This model :cite:p:`hildebrandt2004nonlocal,hildebrandt2007nonlocal` solves a three-field nonlocal dielectric boundary system, similar to that in NESSie :cite:p:`kemmer2018nessie`.
It uses the solvent permittivity ``epsilon``, the high-frequency permittivity ``epsilon_inf``, and the positive correlation length ``lambda`` in bohr.
Fortran optionally accepts the solute permittivity ``epsilon_solute``; C and Python take it explicitly.
The host supplies the ordinary electrostatic potential and outward normal derivative of its solute at the surface points.
The energy is the physical reaction energy of a solute embedded in a cavity of relative permittivity ``epsilon_solute``; ``epsilon_solute = 1`` for quantum-mechanical densities.
For a charge ``Q`` at the center of a sphere of radius ``R`` the energy is :math:`\frac{Q^2}{2R}\left[\frac1{\varepsilon} - \frac1{\eta} + \left(\frac1{\varepsilon_\infty} - \frac1\varepsilon\right)\frac1{1+\nu}\right]` with :math:`\nu = \sqrt{\varepsilon/\varepsilon_\infty}\,R/\lambda` and :math:`\eta` = ``epsilon_solute``; in the local limit ``epsilon = epsilon_inf`` it is :math:`\frac{Q^2}{2R}(1/\varepsilon - 1/\eta)`.
The model requires cavity normals and principal curvatures; the curvature enters the self term of the double-layer operator.

The Laplace terms of the right-hand side use the solute potential rather than its normal derivative (the two differ by the exterior Green identity, which holds for charge inside the cavity), and the energy is that of the apparent surface charges in the solute potential.
In the local limit ``epsilon = epsilon_inf`` the model is therefore IEF-PCM on its own surface operators (not those of the IEF-PCM component, which it matches only to the discretization), and solute charge outside the cavity, such as the tail of a QM electron density, is polarized as in IEF-PCM instead of being seen only through its flux though the cavity.
Every source term carries a permittivity contrast, so a uniform medium (``epsilon = epsilon_inf = epsilon_solute``) gives exactly zero energy and response.
The system is solved for the exterior auxiliary trace minus the solute potential, so in the local limit ``epsilon = epsilon_inf`` the Yukawa operators cancel exactly and the energy does not depend on ``lambda``.

.. code-block:: fortran

   call new_component_nonlocal(nonlocal_pcm, epsilon=78.4_wp, &
      epsilon_inf=1.78_wp, lambda=10.0_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_nonlocal_component``.
Python uses ``ModelComponentNonlocal(epsilon, epsilon_inf, lambda_, epsilon_solute, parameters=...)``.
All numerical settings can be found in :doc:`index`.
