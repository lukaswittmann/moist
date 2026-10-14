Electron-Isodensity (ρ-DROP) Cavity
===================================

The electron-isodensity surface (:math:`\rho`-DROP)
:cite:p:`wittmann2026isodensitydrop` is the zero level set of

.. math::

   S(\mathbf r) = \alpha[\rho_\mathrm{iso} - \rho(\mathbf r)].

``rho_iso`` (:math:`\rho_\mathrm{iso}`, default ``1.0e-3`` bohr\ :sup:`-3`) is the
cavity-defining isodensity value.
Values of :math:`\rho_\mathrm{iso}` from :math:`4\times10^{-4}` to
:math:`5\times10^{-3}` are used in the literature.

``scale`` (:math:`\alpha`, default ``1.0e+3``) is a constant multiplier on the level set and all of its
derivatives; this is needed in order to make the default DROP parameters compatible with the electron-isodensity
level set.
``scale`` is not independent of ``rho_iso``: the normalized form :math:`S=1-\rho/\rho_\mathrm{iso}` holds only
while :math:`\alpha=1/\rho_\mathrm{iso}`, so both should be changed together.

.. Note::

   :math:`\rho`-DROP requires a QM host's density; see :ref:`coupling-build` for rebuilding at each density change.
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
   The QM host supplies the GTO basis at construction and the current density matrix through the :ref:`coupling-build` at each SCF step.
   MOIST evaluates the density and its spatial derivatives internally, including shell and spatial screening.

   ``basis`` selects the representation the density is contracted in:
   - ``cartesian`` (default) works for every host.
   - ``spherical`` contracts over the :math:`2l+1` harmonic polynomials of each shell instead of the :math:`(l+1)(l+2)/2` monomials, which is exact for a host with spherical basis functions.
   It is set as ``basis=moist_iso_basis%spherical`` (Fortran, ``use moist, only: moist_iso_basis``), ``options.basis = moist_isodensity_basis_spherical`` (C) or ``IsodensityParameters(basis="spherical")`` (Python).

Both implementations share ``moist_cavity_drop_lsf_isodensity_param_type`` and build the same :math:`S` from it.

ECP nuclei
----------

With an effective core potential, the pseudo-density can fall below :math:`\rho_\mathrm{iso}` at the
nucleus, which gives the level set a small inner cavity.
MOIST therefore adds, at every nucleus whose nuclear charge ``z`` is below its atomic number,

.. math::

   g(\mathbf r) = 100\,\rho_\mathrm{iso}\,\bigl(1 - |\mathbf r - \mathbf R_A| / X_A\bigr)^5,
   \qquad X_A = \min(f R_A,\ 2\ \mathrm{bohr}),

and nothing beyond :math:`X_A` (:math:`R_A` is the atom's cavity radius, :math:`f` is
``core_radius_fraction``, default 0.5).
The fill has a cusp at the nucleus, so the nucleus is no stationary point of the level set a projector could settle on.
As the outer surface is further out, it is left unchanged; a converged surface point inside this area
stops the cavity build with an error naming the atom.
A smaller ``core_radius_fraction`` admits a surface that comes closer to an ECP nucleus.
The host supplies ``z`` through the ``nuclear_charges`` request (:ref:`coupling-build`).
``core_fill`` (default on) switches the fill off: ``param%core_fill`` (Fortran), ``options.core_fill`` (C)
or ``IsodensityParameters(core_fill=False)`` (Python); ``core_radius_fraction`` is set the same way.

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
             CavityDROP, CPCMRadii, GaussianBasis,
             Isodensity, IsodensityParameters,
         )

         basis = GaussianBasis(
             shell_atom=shell_atom, shell_l=shell_l, shell_nprim=shell_nprim,
             exponents=exps, coefficients=coeffs,
         )
         cavity = CavityDROP(
             context=context,
             lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1.0e-3)),
             radii=CPCMRadii(),
             basis=basis,
         )

The Fortran and C level sets still have to be passed to the DROP cavity constructor; see :ref:`drop-construction`.

The density matrix uses the Cartesian-monomial layout reported by
``moist_get_isodensity_cart_layout`` / ``cavity.isodensity_layout()``.

A **callback** level set is available as well
(``callback_lsf%new(...)``, ``moist_new_isodensity_callback_lsf``, or a callable
``source=`` in Python). It returns the bare density and its spatial derivatives;
MOIST forms :math:`S = \alpha(\rho_\mathrm{iso} - \rho)` from its own ``rho_iso``
and ``scale``. See ``isodensity_lsf_callback`` in :doc:`/reference/fortran` for
the callback interface.
