Electron-Isodensity (ρ-DROP) Cavity
===================================

The electron-isodensity surface (:math:`\rho`-DROP)
:cite:p:`wittmann2026isodensitydrop` is the zero level set of

.. math::

   S(\mathbf r) = \alpha[\rho_\mathrm{iso} - \rho(\mathbf r)].

``rho_iso`` (:math:`\rho_\mathrm{iso}`, default ``1.0e-3`` bohr\ :sup:`-3`) is the
cavity-defining isodensity valu.
Values of :math:`\rho_\mathrm{iso}` from :math:`4\times10^{-4}` to
:math:`5\times10^{-3}` are used in the literature.

``scale`` (:math:`\alpha`, default ``1.0e+3``) is a constant multiplier on the level set and all of its
derivatives; this is needed in order to make the default DROP parameters compatible with the electron-isodensity
level set.
``scale`` is not independent of ``rho_iso``: the normalized form :math:`S=1-\rho/\rho_\mathrm{iso}` holds only
while :math:`\alpha=1/\rho_\mathrm{iso}`, so both should be changed together.

.. Note::

   :math:`\rho`-DROP requires the electron density, so it *has* to be coupled to a QM host, and the cavity has to be updated after every density change.
   Components evaluated on a :math:`\rho`-DROP cavity must also implement their surface weights (the derivatives of the model energy with respect to the surface quantities ) to obtain the variational cavity response.
   These weights are already available for every currently implemented model and component.
   Those weights are contracted directly with the DROP adjoint.

Implementations
---------------

Two implementations are available:

**Callback** (``moist_cavity_drop_lsf_isodensity_callback_type``)
   The QM host evaluates the electron density and its spatial derivatives at points requested by MOIST: the host returns :math:`\rho`.
   The density and its gradient are always requested; the Hessian and third derivative arrive as NULL pointers whenever DROP does not need them, and the callee must (should) skip computing them.
   The third derivative is needed for the cavity response and for analytic nuclear derivatives.
   This route is the slower one: per-point call overhead and no vectorization across points.

**Internal** (``moist_cavity_drop_lsf_isodensity_internal_type``)
   The QM host supplies the GTO basis set and, for every SCF step, the current density matrix in Cartesian-monomial layout (``moist_get_isodensity_cart_layout``, ``moist_set_isodensity_density``).
   MOIST evaluates the density and its spatial derivatives internally, including shell and spatial screening.

Both implementations share ``moist_cavity_drop_lsf_isodensity_param_type`` and build the same :math:`S` from it.

Construction
------------

The tabs below build the **internal** LSF, which owns a Cartesian-monomial GTO
basis and evaluates the density itself. Shell atom indices are zero-based, and
``exps``/``coeffs`` hold ``sum(shell_nprim)`` values with primitive
normalization folded into the coefficients.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_cavity_drop_lsf_isodensity_internal, only : &
            & moist_cavity_drop_lsf_isodensity_internal_type
         use moist, only : moist_cavity_drop_lsf_isodensity_param_type

         type(moist_cavity_drop_lsf_isodensity_internal_type) :: rho_lsf

         ! sh_atom, sh_l, sh_nprim, exps and coeffs come from the QM host
         call rho_lsf%new(sh_atom, sh_l, sh_nprim, exps, coeffs, error=error, &
            & param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=1.0e-3_wp))
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_isodensity_options options;
         moist_init_isodensity_options(error, &options, sizeof options);
         options.rho_iso = 1.0e-3;

         moist_lsf lsf = moist_new_isodensity_lsf(
             error, nshell, shell_atom, shell_l, shell_nprim,
             exps, coeffs, &options);

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import (
             CavityDROP, CPCMRadii, GaussianBasis, InternalDensity,
             Isodensity, IsodensityParameters,
         )

         basis = GaussianBasis(
             shell_atom=shell_atom, shell_l=shell_l, shell_nprim=shell_nprim,
             exponents=exps, coefficients=coeffs,
         )
         cavity = CavityDROP(
             lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)),
             radii=CPCMRadii(),
             source=InternalDensity(basis, density_matrix),
         )

The Fortran and C level sets still have to be passed to the DROP cavity constructor; see :ref:`drop-construction`.

The cavity constructor *copies* the level set, so a new density must be set on the cavity-owned copy, not on the LSF built above.
In Fortran that works via ``cavity%lsf_model``:

.. code-block:: fortran

   select type (lsf => cavity%lsf_model)
   type is (moist_cavity_drop_lsf_isodensity_internal_type)
      call lsf%set_density(dcart, error)
   end select
   if (allocated(error)) error stop error%message

C uses ``moist_set_isodensity_density`` on the cavity, and Python assigns
``InternalDensity.density_matrix`` on the source object the cavity holds. The
Cartesian-monomial layout the matrix must use is reported by
``moist_get_isodensity_cart_layout`` / ``cavity.isodensity_layout()``.

A **callback** level set is available as well
(``callback_lsf%new(...)``, ``moist_new_isodensity_callback_lsf``, or a callable
``source=`` in Python). It returns the bare density and its spatial derivatives;
MOIST forms :math:`S = \alpha(\rho_\mathrm{iso} - \rho)` from its own ``rho_iso``
and ``scale``. See ``isodensity_lsf_callback`` in :doc:`/reference/fortran` for
the callback interface.
