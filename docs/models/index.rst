Models
======

This section describes solvation models in MOIST.
The shared model base is in ``src/moist/model/type.f90``; models are organized by their theory:

* ``model/continuum`` contains cavity-based models and components.
* ``model/moz`` contains molecular Ornstein Zernike-type theories, including the Reference Interaction Site Model (RISM) theories.
   * ``model/moz/1d`` contains the VV and UV 1D-RISM model.
   * ``model/moz/3d`` contains the UV 3D-RISM model.

Solvation Model Interface
-------------------------

A ``solvation_model_type`` is the host-facing object for one complete solvation treatment.
Concrete models implement four deferred procedures:

``update(mol, error)``
   Refresh all structure-dependent state.

``get_energy(coupling, energy, error)``
   Add the model energy for the answers staged on the coupling.

``get_response(coupling, response, error)``
   Add the host contributions of the response phase.

``get_gradient(coupling, response, gradient, error)``
   Add the nuclear-gradient contribution, and with it the host part of the
   gradient phase.

The ``coupling_type`` holds the requests the host answers, visited one at a time with ``next()``, and ``response_type`` the list of items the model hands back, walked the same way; both are described in :doc:`/reference/coupling`.
``model_continuum_type`` (``src/moist/model/continuum/type.f90``) adds the calls that drive them: ``new_coupling`` returns a model-owned coupling (several may coexist, and ``release_coupling`` frees one) then ``prepare_energy``, ``prepare_response`` and ``prepare_gradient`` stage one phase each.

Model Component Interface
-------------------------

A ``model_continuum_component_type`` is one reusable energy term evaluated on a shared cavity.
It stores a name, the current solute structure, and a linear ``scale``, which only GOSTSHYP currently applies.

Components implement ``update``, ``get_energy``, ``get_response``, and ``get_gradient`` with the live ``cavity_type`` as an additional argument.
A bare component can be driven on its own coupling with the same ``new_coupling`` and ``prepare_*`` calls as a model, each taking the live cavity as its first argument.
Components may also override default hooks:

- ``declare_coupling`` registers the host requests the component reads, per phase, together with their inputs, such as GOSTSHYP's Gaussian widths;
- ``get_trace_response`` emits the response items that follow from the potential adjoint;
- ``get_surface_weights`` accumulates the component's own cavity surface weights;
- ``get_host_surface_weights`` folds the weights the host answered into the same accumulator;
- ``get_gradient_surface_weights`` does the same for the gradient path, and defaults to ``get_surface_weights``;
- ``get_direct_gradient`` adds nuclear-gradient terms that do not flow through the surface.

Composition and Lifecycle
-------------------------

The ``model_continuum_type`` owns one cavity and an ordered list of components:

1. Construct the model from a cavity and add all components before the first update; the model stores copies of both.
2. ``update`` refreshes the cavity first, then every component.
3. ``get_energy`` sums the component energies and keeps each one in the model until the next update, as the ``energy`` field of ``list_component_fields(index, query)``.
4. ``get_response`` collects each component's direct items (``potential_adjoint``, ``gaussian_amplitude``) and its surface weights, then lets the cavity contract the accumulated weights, which adds the ``density`` item for a cavity whose surface follows the density.
5. ``get_gradient`` adds the direct nuclear terms, contracts the gradient-side surface weights through the cavity, and refills the response with the same direct items.

Construction Example
~~~~~~~~~~~~~~~~~~~~

This example constructs a list-based model containing :doc:`CPCM </models/components/pcm/cpcm>` and a :doc:`pressure-volume term </models/components/pv>` on a :doc:`SvdW-DROP cavity </cavities/svdw>`.

.. code-block:: fortran

   use mctc_env, only: wp, error_type
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: &
      & moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii
   use moist_model_continuum, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component, only: model_continuum_component_cpcm, &
      & new_component_cpcm, model_continuum_component_pv, new_component_pv

   type(moist_context_type), target :: ctx
   type(moist_cavity_drop_lsf_svdw_type) :: svdw
   type(cavity_type_drop) :: cavity
   type(model_continuum_component_cpcm) :: electrostatic
   type(model_continuum_component_pv) :: pressure_volume
   type(model_continuum_type), target :: model
   type(error_type), allocatable :: error

   ! Global context
   call new_context(ctx, nthreads=0)
   
   ! Construct cavity and its level set
   call svdw%new()
   call new_cavity_drop(cavity, ctx, &
      & radius_model=default_cpcm_radii(), lsf_model=svdw, error=error)
   if (allocated(error)) error stop error%message

   ! Construct CPCM (water)
   call new_component_cpcm(electrostatic, ctx, epsilon=80.0_wp, error=error)
   if (allocated(error)) error stop error%message
   ! Construct pressure model (1 GPa)
   call new_component_pv(pressure_volume, pressure=3.39893E-5_wp)

   ! Construct model
   call new_continuum_model(model, cavity, ctx, error)
   if (allocated(error)) error stop error%message

   ! Add model components
   call model%add_component(electrostatic, error)
   if (allocated(error)) error stop error%message
   call model%add_component(pressure_volume, error)
   if (allocated(error)) error stop error%message

Call ``model%update(mol, error)`` with the current structure, then create a host coupling with ``model%new_coupling(coupling, error)``; energies, responses and gradients are computed one phase at a time on that coupling, as described in :doc:`/reference/coupling`.
``model%component_count()``, ``model%component_name(index, name)`` and ``model%component_description(index, description)`` describe the components before the first update; after ``get_energy``, each component's contribution is a named field (:doc:`/reference/fields`).
``model%print_parameters(unit)`` prints the settings at any time: the cavity section followed by each component.; ``unit`` defaults to the run context's unit.

Components
----------

.. toctree::
   :maxdepth: 2

   components/index
