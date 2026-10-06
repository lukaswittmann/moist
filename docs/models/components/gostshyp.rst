GOSTSHYP
========

GOSTSHYP (Gaussians On Surface Tesserae Simulating HYdrostatic Pressure) applies pressure through Gaussian potentials on the cavity surface.
Their amplitudes depend on the electron density and the requested pressure.

For each surface point with position :math:`\mathbf C_i`, area :math:`a_i` and normal :math:`\mathbf n_i`,

.. math::

   G_i(\mathbf r) &= \exp[-\omega_i |\mathbf r-\mathbf C_i|^2],
   \qquad \omega_i = \frac{\pi\ln 2}{a_i}, \\
   g_i &= \int \rho(\mathbf r)G_i(\mathbf r)\,d\mathbf r, \\
   f_i &= \int \rho(\mathbf r)\,\mathbf n_i\cdot\nabla_{\mathbf r}G_i(\mathbf r)\,d\mathbf r,

.. math::

   E_\mathrm{GOSTSHYP} = \sum_i \frac{p\,a_i}{f_i}g_i.

Points with negligible :math:`|f_i|` are excluded.

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

Construction
------------

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_model_continuum_component, only : model_continuum_component_gostshyp, &
            & new_component_gostshyp

         type(model_continuum_component_gostshyp) :: gostshyp

         call new_component_gostshyp(gostshyp, pressure=1.699465e-3_wp)  ! 50 GPa

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_component gostshyp = moist_new_gostshyp_component(error, NULL, 1.699465e-3);
         moist_add_model_component(error, model, gostshyp);
         moist_delete(gostshyp);  /* model owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import ModelComponentGOSTSHYP
         from moist.pyscf import GPA_TO_AU

         gostshyp = ModelComponentGOSTSHYP(50.0 * GPA_TO_AU)

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
