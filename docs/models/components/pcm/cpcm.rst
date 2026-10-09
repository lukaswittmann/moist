CPCM
====

The conductor-like polarizable continuum model uses

.. math::

   f(\varepsilon) = \frac{\varepsilon-1}{\varepsilon}.

``new_component_cpcm`` requires :math:`\varepsilon\geq1` and takes a ``moist_pcm_parameters_type`` for the solver settings and an optional external matrix.
An infinite :math:`\varepsilon` selects the conductor limit :math:`f(\varepsilon)=1`.

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
         use moist_model_continuum_component, only : model_continuum_component_cpcm, &
            & new_component_cpcm, solver_type, moist_pcm_parameters_type

         type(model_continuum_component_cpcm) :: cpcm

         call new_component_cpcm(cpcm, epsilon=78.4_wp, error=error, &
            & param=moist_pcm_parameters_type(solver=solver_type%cholesky), ctx=ctx)
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_pcm_options options;
         moist_init_pcm_options(error, &options, sizeof options);
         options.solver = moist_pcm_solver_cholesky;

         moist_component cpcm = moist_new_cpcm_component(error, NULL, 78.4, &options);
         moist_add_model_component(error, model, cpcm);
         moist_delete(cpcm);  /* model owns a copy */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import ModelComponentCPCM, PCMParameters, PCMSolver

         cpcm = ModelComponentCPCM(
             78.4, parameters=PCMParameters(solver=PCMSolver.CHOLESKY),
         )

Passing no solver configuration (an omitted ``param``, ``NULL`` options, or no ``parameters``) selects the compiled defaults.
Adding the component to a model and driving the evaluation phases is the same for every component:
see :ref:`coupling-fortran` for Fortran, :doc:`/reference/c` for handle ownership, and :doc:`/reference/python` for the Python classes (:py:class:`moist.ModelComponentCOSMO`).

