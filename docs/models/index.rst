Models
======

This section describes the solvation models and reusable model components in MOIST.
Their abstract interfaces are defined in ``src/moist/type.f90``.

Solvation Model Interface
-------------------------

A ``solvation_model_type`` is the host-facing object for one complete solvation treatment.
Concrete models implement four deferred procedures:

``update(mol, error)``
   Refresh all structure-dependent state.

``get_energy(coupling, energy, error)``
   Add the model energy for the supplied QM coupling data.

``get_response(coupling, response, error)``
   Add the self-consistent host-potential contributions.

``get_gradient(coupling, gradient, error)``
   Add the nuclear-gradient contribution.

The ``coupling_type`` carries QM data from the host to the model;
``response_type`` returns the corresponding response channels.

Model Component Interface
-------------------------

A ``solvation_model_component_type`` is one reusable energy term evaluated on a shared cavity.
It stores a name, the current solute structure, and a linear ``scale`` as shared component state.

Components implement ``update``, ``get_energy``, ``get_response``, and ``get_gradient`` with the live ``cavity_type`` as an additional argument.
They may also override three default no-op response hooks:

- ``get_trace_response`` for direct host-trace adjoints;
- ``get_surface_weights`` for model-specific cavity surface weights;
- ``get_host_surface_weights`` for host-computed trace-geometry weights.

Composition and Lifecycle
-------------------------

The current ``solvation_model_general`` owns one authoritative cavity and an ordered list of components:

1. Construct the model from a cavity and add all components before the first update.
2. ``update`` refreshes the cavity first, then every component.
3. Energy and gradient calls sum the component contributions.
4. Potential evaluation accumulates direct terms and surface weights, then lets the cavity map the surface response to host-potential channels.

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
   use moist_model_general, only: solvation_model_general, new_model_general
   use moist_model_components, only: solvation_model_component_cpcm, &
      & new_component_cpcm, solvation_model_component_pv, new_component_pv

   type(moist_context_type), target :: ctx
   type(moist_cavity_drop_lsf_svdw_type) :: svdw
   type(cavity_type_drop) :: cavity
   type(solvation_model_component_cpcm) :: electrostatic
   type(solvation_model_component_pv) :: pressure_volume
   type(solvation_model_general) :: model
   type(error_type), allocatable :: error

   ! Global context
   call new_context(ctx)
   
   ! Construct cavity and its level set
   call svdw%new()
   call new_cavity_drop(cavity, ctx, &
      & radius_model=default_cpcm_radii(), lsf_model=svdw, error=error)
   if (allocated(error)) error stop error%message

   ! Construct CPCM (water)
   call new_component_cpcm(electrostatic, ctx, epsilon=80.0_wp, error=error)
   if (allocated(error)) error stop error%message
   ! Construct pressure model (1 GPa)
   call new_component_pv(pressure_volume, pressure=3.40E-5_wp, error=error)
   if (allocated(error)) error stop error%message

   ! Construct model
   call new_model_general(model, cavity, ctx, error)
   if (allocated(error)) error stop error%message

   ! Add model components
   call model%add_component(electrostatic, error)
   if (allocated(error)) error stop error%message
   call model%add_component(pressure_volume, error)
   if (allocated(error)) error stop error%message

Call ``model%update(mol, error)`` with the current structure, then energies, potentials, or gradients can be computed.

Components
----------

.. toctree::
   :maxdepth: 2

   components/index
