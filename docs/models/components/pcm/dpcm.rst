DPCM
=====

The differential PCM boundary formulation :cite:p:`miertus1981dpcm` uses the surface potential and its
normal derivative. Both the current and original discretizations have separate
component types. They require cavity normals and share the PCM numerical
options of :doc:`index`.

.. code-block:: fortran

   call new_component_dpcm(dpcm, epsilon=78.4_wp, error=error, ctx=ctx)
   call new_component_dpcm_original(dpcm_original, epsilon=78.4_wp, error=error, ctx=ctx)

The C constructors are ``moist_new_dpcm_component`` and
``moist_new_dpcm_original_component``. Python provides ``ModelComponentDPCM``
and ``ModelComponentDPCMOriginal``. Each takes ``epsilon`` and optional
numerical settings. See :doc:`/reference/coupling` for the host derivatives.
