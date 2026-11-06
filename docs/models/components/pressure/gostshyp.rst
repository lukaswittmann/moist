Gaussians On Surface Tesserae Simulating HYdrostatic Pressure (GOSTSHYP)
========================================================================

GOSTSHYP :cite:p:`scheurer2021gostshyp,pausch2025gostshyp` applies pressure through Gaussian potentials on the cavity surface.
Their amplitudes depend on the electron density and the requested pressure.

Each surface point :math:`i` of the cavity has a position :math:`\mathbf C_i`, an area :math:`a_i` and an outward unit normal :math:`\mathbf n_i`.
On it sits a unit-integral Gaussian

.. math::

   G_i(\mathbf r) = \left(\frac{\omega_i}{\pi}\right)^{3/2}
      \exp[-\omega_i |\mathbf r-\mathbf C_i|^2],
   \qquad \omega_i = \frac{\pi\ln 2}{a_i}.

Two moments of the electron density :math:`\rho` enter the energy:

.. math::

   g_i &= \int \rho(\mathbf r)\,G_i(\mathbf r)\,d\mathbf r, \\
   f_i &= \int \rho(\mathbf r)\,\mathbf n_i\cdot\nabla_{\mathbf r}G_i(\mathbf r)\,d\mathbf r.

The overlap :math:`g_i` (bohr\ :sup:`-3`) is the density seen by the Gaussian.
The normal-gradient trace :math:`f_i` (bohr\ :sup:`-4`) measures how strongly the density changes along the normal under the Gaussian; density inside the cavity contributes positively, so :math:`f_i` is normally positive.

A potential :math:`c_i G_i` exerts the normal force :math:`-c_i f_i` on the density.
Requiring this force to push the density inward with the applied pressure :math:`p` on the area :math:`a_i`,
:math:`c_i f_i = p\,a_i`, fixes the amplitude :math:`c_i = p\,a_i/f_i`, and the energy :math:`c_i g_i` of all points is

.. math::

   E_\mathrm{GOSTSHYP} = \sum_i s(a_i)\,p\,a_i\,g_i\,R(f_i),
   \qquad
   R(f) = \frac{S(|f|)}{f}.

:math:`R` is the regularized reciprocal of the trace and :math:`s(a_i)` switches off narrow Gaussians, both described below.
Because the amplitude depends on the density through :math:`f_i`, the Fock matrix and the gradients differentiate through both moments.

The unit-integral normalization follows Pausch's normalized surface Gaussian convention.
All moments and matching host AO integrals use that normalization.
The switch :math:`S` is the C-infinity partition-of-unity bump of the DROP switching module.
It is zero for :math:`|f| \le` ``regularization_start`` and one for :math:`|f| \ge` ``regularization_end``.
Above the end the energy is exactly the pressure constraint :math:`p a_i g_i/f_i`; below the start a point is off.

Narrow Gaussians are switched off smoothly by area, before host integrals lose precision.
:math:`s(a_i)` is a C2 polynomial that is one for :math:`\omega_i \le 10^3` bohr\ :sup:`-2` and zero for :math:`\omega_i \ge 10^4` bohr\ :sup:`-2` (areas below about ``2.2e-4 bohr^2``).

Negative pressure amplitudes are retained by default.
Enabling ``suppress_negative_amplitudes`` sets that branch to zero within the calculation;
the trace switch is identically zero around :math:`f = 0`, so the cut is smooth.

Each energy evaluation reports regularized, negative candidate, suppressed and switched-off point counts at :ref:`verbosity <verbosity>` 4.
A negative-amplitude warning always prints.

The component supports static DROP, isodensity DROP and iSwiG cavities, including Fock contributions and analytic  nuclear gradients.

Parameters
----------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Name
     - Unit
     - Description
   * - ``pressure``
     - hartree/bohr^3
     - Applied pressure. One GPa is ``3.39893e-5`` hartree/bohr^3.
   * - ``regularization_start``
     - bohr^-4
     - Trace magnitude at and below which a point is off, default ``1e-12``.
   * - ``regularization_end``
     - bohr^-4
     - Trace magnitude at and above which the energy is the plain reciprocal, default ``1e-10``.
   * - ``suppress_negative_amplitudes``
     - boolean
     - Smoothly remove negative pressure amplitudes, default ``false``.

Construction
------------

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_model_continuum_component, only : model_continuum_component_gostshyp, &
            & new_component_gostshyp, moist_gostshyp_parameters_type

         type(model_continuum_component_gostshyp) :: gostshyp
         type(moist_gostshyp_parameters_type) :: parameters

         parameters%regularization_start = 1.0e-12_wp
         parameters%regularization_end = 1.0e-10_wp
         call new_component_gostshyp(gostshyp, pressure=1.699465e-3_wp, param=parameters)

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_gostshyp_options options;
         moist_init_gostshyp_options(error, &options, sizeof(options));
         options.regularization_start = 1e-12;
         options.regularization_end = 1e-10;
         moist_component gostshyp = moist_new_gostshyp_component(error, NULL, 1.699465e-3, &options);
         moist_add_model_component(error, model, gostshyp);
         moist_delete(gostshyp);  /* model owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import GOSTSHYPParameters, ModelComponentGOSTSHYP
         from moist.pyscf import GPA_TO_AU

         gostshyp = ModelComponentGOSTSHYP(
             50.0 * GPA_TO_AU,
             parameters=GOSTSHYPParameters(regularization_start=1e-12, regularization_end=1e-10))

On water, fluoroacetate and glycine zwitterion the smallest normal gradient
trace is about ``1.5e-3 bohr^-4``, so the default window regularizes no point.

``NULL`` options select compiled defaults in C. Fortran construction accepts an optional ``error``;
evaluation also validates settings. The native settings extend
``moist_model_parameters_type`` and support its JSON/TOML input/output,
validation, default-reset and registered printing interface, just like PCM.

The component cannot form its own density traces. Hosts answer the ``gaussian_moments`` request in every phase and contract the ``gaussian_amplitude`` item of the response; see the :ref:`host requests <coupling-requests>` and :ref:`gostshyp-fortran` for the outputs required per phase.

PySCF
-----

Given a PySCF molecule ``mol``:

.. code-block:: python

   from moist import Context, DROPParameters, ModelComponentGOSTSHYP
   from moist.pyscf import DROP, GPA_TO_AU, Isodensity

   mf = mol.RHF().MOIST(
       cavity=DROP(lsf=Isodensity(), parameters=DROPParameters(nleb=194)),
       components=[ModelComponentGOSTSHYP(50.0 * GPA_TO_AU)],
       context=Context(),
   )
   mf.kernel()
   gradient = mf.nuc_grad_method().kernel()

See :doc:`/reference/pyscf` for supported methods and cavity settings.

Adding the component to a model and driving the evaluation phases is the same for every component:
see :ref:`coupling-fortran` for Fortran, :doc:`/reference/c` for handle ownership, and :doc:`/reference/python` for the Python classes (:py:class:`moist.ModelComponentGOSTSHYP`).
