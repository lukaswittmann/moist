Solvation-Layer Interface Condition (SLIC)
==========================================

The solvation-layer interface condition model :cite:p:`bardhan2014slic,molavitabrizi2017slic` modifies the dielectric boundary response.
Its physical inputs are ``epsilon``, ``alpha``, ``beta`` and ``gamma``.
Fortran optionally accepts ``mu`` and ``phi_static``; when ``mu`` is omitted, it uses ``alpha*tanh(gamma)``.
The C and Python constructors take ``mu`` and ``phi_static`` explicitly.
The static term adds ``phi_static*Q``, where ``Q`` is the net charge of the structure: ``mol%charge`` in Fortran, the ``charge`` argument of ``moist_new_structure`` in C, ``Structure(charge=)`` in Python (the PySCF host passes ``mol.charge``). SLIC requires cavity normals and curvatures.

Use a smooth cavity such as :doc:`SvdW-DROP </cavities/svdw>` or :doc:`isodensity-DROP </cavities/isodensity>`.
On an :doc:`iSwiG </cavities/iswig>` cavity the concave seams between overlapping spheres can fold the nonlinear interface condition: the Newton solve then fails, or converges to a spurious solution.
This is especially prominent for anions.
SLIC warns through the context when it is configured on an iSwiG cavity.

.. code-block:: fortran

   call new_component_slic(slic, epsilon=78.4_wp, alpha=0.1_wp, &
      beta=0.2_wp, gamma=0.3_wp, error=error, ctx=ctx)

The C constructor is ``moist_new_slic_component``; Python provides
``ModelComponentSLIC(epsilon, alpha, beta, gamma, mu, phi_static,
parameters=...)``. Nonlinear solver controls belong to the PCM options of
:doc:`index`.
