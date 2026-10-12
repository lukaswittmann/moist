SLIC
====

The solvation-layer interface condition model :cite:p:`molavitabrizi2017slic` modifies the dielectric boundary
response. Its physical inputs are ``epsilon``, ``alpha``, ``beta`` and
``gamma``. Fortran optionally accepts ``mu`` and ``phi_static``; when ``mu``
is omitted, it uses ``alpha*tanh(gamma)``. The C and Python constructors take
``mu`` and ``phi_static`` explicitly. SLIC requires cavity normals and
curvatures.

.. code-block:: fortran

   call new_component_slic(slic, epsilon=78.4_wp, alpha=0.1_wp, &
      beta=0.2_wp, gamma=0.3_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_slic_component``; Python provides
``ModelComponentSLIC(epsilon, alpha, beta, gamma, mu, phi_static,
parameters=...)``. Nonlinear solver controls belong to the PCM options of
:doc:`index`.
