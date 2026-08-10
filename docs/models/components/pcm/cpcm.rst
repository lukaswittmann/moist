CPCM
====

The conductor-like polarizable continuum model uses

.. math::

   f(\epsilon) = \frac{\epsilon-1}{\epsilon}.

``new_cpcm`` requires :math:`\epsilon\geq1` and accepts optional solver, potential-source, and external-matrix settings. The default configuration uses a Cholesky solver and an electrostatic potential computed from atomic charges; the conductor limit is supported with :math:`f(\epsilon)=1`.
