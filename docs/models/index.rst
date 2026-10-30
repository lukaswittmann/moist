Models
======

This section describes solvation models in MOIST.
The shared model base is in ``src/moist/model/type.f90``; models are organized by their theory:

* ``model/continuum`` contains the cavity-based :doc:`continuum models <continuum>` and their components.
* ``model/moz`` contains :doc:`molecular Ornstein-Zernike <moz>`-type theories, including the Reference Interaction Site Model (RISM) theories.
   * ``model/moz/1d`` contains the VV and UV 1D-RISM model.
   * ``model/moz/3d`` contains the UV 3D-RISM model.

Solvation Model Interface
-------------------------

A ``solvation_model_type`` is the host-facing object for one complete solvation treatment.
Its main host-facing calls are:

``update(mol, error)``
   Prepare structure-dependent state; see :ref:`coupling-build` for deferred construction.

``get_energy(coupling, energy, error)``
   Add the model energy for the answers staged on the coupling.

``get_response(coupling, response, error)``
   Fill the response with coefficients for the host.

``get_gradient(coupling, response, gradient, error)``
   Add the nuclear-gradient contribution, and with it the host part of the
   gradient phase.

The ``coupling_type`` holds the requests the host answers, visited one at a time with ``next()``, and ``response_type`` the list of items the model hands back, walked the same way; both are described in :doc:`/reference/coupling`.
The base provides ``new_coupling``, ``release_coupling`` and the ``prepare_*`` calls.
Concrete models implement ``update``, ``atom_count``, ``declare_pass`` and the three evaluation calls; a model whose update waits for host answers overrides ``build``, calling the base ``model_build`` first.
``new_coupling`` requires an updated model and the ``prepare_*`` evaluation calls a built one.

.. toctree::
   :maxdepth: 2
   :titlesonly:

   continuum
   components/index
   moz
