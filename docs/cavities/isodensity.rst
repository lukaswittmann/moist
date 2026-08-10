Electron-Isodensity (ρ-DROP) Cavity
===================================

The electron-isodensity surface (:math:`\rho`-DROP)
:cite:p:`wittmann2026isodensitydrop` is the zero level set of

.. math::

   S(\mathbf r) = -\alpha[\rho(\mathbf r) - \rho_\mathrm{iso}].

Here ``rho_iso`` is the cavity-defining isodensity value and ``scale`` is
:math:`\alpha`. A common choice is :math:`\rho_\mathrm{iso}=10^{-3}` atomic
units; values from :math:`4\times10^{-4}` to :math:`5\times10^{-3}` are used in
the literature. Setting :math:`\alpha=1/\rho_\mathrm{iso}` gives the normalized
form :math:`S=1-\rho/\rho_\mathrm{iso}`.

.. Note::

   As :math:`\rho`-DROP requires the electron density, it has to be coupled with a QM host. The cavity needs to be updated for every change in the density.
   Additionally, for components that are evaluated on the :math:`\rho`-DROP cavity, the model-specific surface weights must be implemented to obtain the variational cavity responsse.
   These surface weights are the derivative of the model energy with respect to the surface quantities (e.g. grid point position, grid point area,...).
   After, these can directly be contrated with the DROP adjoint to obtain the variational cavity model response.

Implementations
---------------

Two implementations are available:

**Callback** (``moist_cavity_drop_lsf_isodensity_callback_type``)
   The QM host evaluates the level set and its spatial derivatives at points requested by MOIST.
   Generally, value and gradient are always required, the Hessian is additionally required and third derivatives are needed for the cavity response and analytic nuclear derivatives.
   The callback route is generally slower due to call overhead, no easy vectorization/missing batching (tbd), and screening and possibly problematic parallelization.

**Internal** (``moist_cavity_drop_lsf_isodensity_internal_type``)
   The QM host supplies the GTO basis set and for every SCF step the current density matrix (in Cartesian-monomial layout, see ``moist_get_isodensity_cart_layout``, and ``moist_set_isodensity_density``).
    MOIST then evaluates the density and its spatial derivatives internally, including shell and spatial screening.
