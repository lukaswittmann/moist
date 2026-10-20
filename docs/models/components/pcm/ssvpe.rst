SS(V)PE
========

The surface and simulation of volume polarization equation :cite:p:`chipman2000ssvpe` is a dielectric
boundary formulation, the symmetrized counterpart of IEF-PCM :cite:p:`lange2011symmetric`, that uses the surface potential and its normal derivative.
It requires cavity normals. Construct it with ``epsilon`` through
``new_component_ssvpe`` (Fortran), ``moist_new_ssvpe_component`` (C), or
``ModelComponentSSVPE`` (Python); numerical controls use the shared PCM options.

.. code-block:: fortran

   call new_component_ssvpe(ssvpe, epsilon=78.4_wp, error=error, ctx=ctx)

See :ref:`coupling-requests` for the host's raw potential derivatives and
:ref:`coupling-response` for the adjoint contraction.
