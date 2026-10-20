DPCM
=====

The dielectric PCM (D-PCM) :cite:p:`miertus1981dpcm` finds the surface charges from the normal derivative of the solute potential; the potential itself enters only the energy.
Two component types solve the same discretized equation with different row scalings: the area-scaled form and the original form as published.
They require cavity normals and share the PCM numerical options of :doc:`index`.

.. code-block:: fortran

   call new_component_dpcm(dpcm, epsilon=78.4_wp, error=error, ctx=ctx)
   call new_component_dpcm_original(dpcm_original, epsilon=78.4_wp, error=error, ctx=ctx)

The C constructors are ``moist_new_dpcm_component`` and
``moist_new_dpcm_original_component``. Python provides ``ModelComponentDPCM``
and ``ModelComponentDPCMOriginal``. Each takes ``epsilon`` and optional
numerical settings. See :doc:`/reference/coupling` for the host derivatives.
