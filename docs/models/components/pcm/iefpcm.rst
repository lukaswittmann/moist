IEF-PCM
========

The integral-equation formulation :cite:p:`cances1997iefpcm,mennucci1997ief` solves a dielectric boundary equation from
the solute's surface potential alone. Its double-layer operator requires cavity normals.
Pass the solvent dielectric constant ``epsilon`` to
``new_component_iefpcm`` (Fortran), ``moist_new_iefpcm_component`` (C), or
``ModelComponentIEFPCM`` (Python). All three accept the shared PCM numerical
options described in :doc:`index`.

.. code-block:: fortran

   call new_component_iefpcm(iefpcm, epsilon=78.4_wp, error=error, ctx=ctx)

The host supplies the raw potential derivatives described in
:ref:`coupling-requests` and contracts both potential-adjoint arrays when
present.
