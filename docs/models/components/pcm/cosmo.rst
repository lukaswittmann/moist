COSMO
=====

The conductor-like screening model uses

.. math::

   f(\varepsilon) = \frac{\varepsilon-1}{\varepsilon+\tfrac{1}{2}}.

``new_component_cosmo`` requires :math:`\varepsilon\geq1` and, unlike CPCM, has no conductor limit. It otherwise shares the solver, request, external-matrix, energy, response and gradient machinery of ``model_continuum_component_pcm``.

Parameters
----------

Besides ``epsilon``, the constructors take the shared PCM settings of :doc:`index`: ``solver``, ``solver_tol``, ``solver_maxiter``.

Construction
------------

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use mctc_env, only : wp
         use moist_model_continuum_component, only : model_continuum_component_cosmo, &
            & new_component_cosmo, solver_type, moist_pcm_parameters_type

         type(model_continuum_component_cosmo) :: cosmo

         call new_component_cosmo(cosmo, ctx, epsilon=78.4_wp, error=error, &
            & param=moist_pcm_parameters_type(solver=solver_type%cholesky))
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_pcm_options options;
         moist_init_pcm_options(error, &options, sizeof options);
         options.solver = moist_pcm_solver_cholesky;

         moist_component cosmo = moist_new_cosmo_component(error, 78.4, &options);
         moist_add_model_component(error, model, cosmo);
         moist_delete(cosmo);  /* model owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import ModelComponentCOSMO, PCMParameters, PCMSolver

         cosmo = ModelComponentCOSMO(
             78.4, parameters=PCMParameters(solver=PCMSolver.CHOLESKY),
         )

Passing no solver configuration (an omitted ``param``, ``NULL`` options, or no ``parameters``) selects the compiled defaults.
Adding the component to a model and driving the evaluation phases is the same for every component:
see :ref:`coupling-fortran` for Fortran, :doc:`/reference/c` for handle ownership, and :doc:`/reference/python` for the Python classes (:py:class:`moist.ModelComponentCOSMO`).
