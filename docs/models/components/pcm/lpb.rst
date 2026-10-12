LPB
===

The linearized Poisson–Boltzmann model :cite:p:`juffer1991lpb` includes ionic screening. Its physical
inputs are the solvent dielectric constant ``epsilon`` and the inverse Debye
length ``kappa`` in bohr⁻¹. ``kappa`` must be nonnegative. The Fortran helper
``debye_kappa(epsilon, ionic_strength, temperature)`` can calculate it from
ionic strength and temperature; the C and Python constructors take ``kappa``
directly.

.. code-block:: fortran

   call new_component_lpb(lpb, epsilon=78.4_wp, kappa=0.1_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_lpb_component``; Python uses
``ModelComponentLPB(epsilon, kappa, parameters=...)``. Numerical settings use
the common PCM options of :doc:`index`.
