COSMO
=====

The conductor-like screening model uses

.. math::

   f(\epsilon) = \frac{\epsilon-1}{\epsilon+\tfrac{1}{2}}.

``new_component_cosmo`` requires :math:`\epsilon\geq1` and otherwise shares the solver, potential-source, external-matrix, energy, response, and gradient machinery of ``solvation_model_component_pcm``.

