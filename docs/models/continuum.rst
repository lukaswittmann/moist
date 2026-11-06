.. _continuum-models:

Continuum Models
================

A continuum model describes the solvent as a continuum around a molecule-shaped cavity.
In MOIST it is not a fixed method but a composition: one :doc:`cavity </cavities/index>` and an ordered list of :doc:`components <components/index>`, each of which adds one physical term evaluated on that cavity.

.. math::

   E_\mathrm{model} = \sum_k E_k[\text{cavity}, \rho]

The same building blocks thus give many models.
For example,

- :doc:`CPCM <components/pcm/cpcm>` on an :doc:`iSwiG </cavities/iswig>` cavity is a standard conductor-like PCM;
- :doc:`IEF-PCM <components/pcm/iefpcm>` and :doc:`PV <components/pressure/pv>` on :doc:`SvdW-DROP </cavities/svdw>` add the pressure-volume work of an applied pressure;
- :doc:`CPCM <components/pcm/cpcm>` and :doc:`GOSTSHYP <components/pressure/gostshyp>` on :doc:`ρ-DROP </cavities/isodensity>` give a fully variational, density-following cavity under hydrostatic pressure.

The components share the cavity but are otherwise independent.
Each enables the surface geometry it needs, such as normals or curvatures, on the model's copy of the cavity, and each contribution is reported separately.
On a cavity whose surface follows the density, the model collects the surface weights of all components and lets the cavity contract them once.
Which combinations work is listed in the support tables of the :doc:`cavities </cavities/index>` and the :doc:`components <components/index>`.

Construction
------------

The examples build CPCM for water plus a pressure-volume term at 1 GPa on a SvdW-DROP cavity.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

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
         call new_cavity_drop(cavity, &
            & radius_model=default_cpcm_radii(), lsf_model=svdw, error=error)
         if (allocated(error)) error stop error%message

         ! Construct CPCM (water)
         call new_component_cpcm(electrostatic, epsilon=80.0_wp, error=error)
         if (allocated(error)) error stop error%message
         ! Construct pressure model (1 GPa)
         call new_component_pv(pressure_volume, pressure=3.39893E-5_wp)

         ! Construct model
         call new_continuum_model(model, ctx, cavity, error)
         if (allocated(error)) error stop error%message

         ! Add model components
         call model%add_component(electrostatic, error)
         if (allocated(error)) error stop error%message
         call model%add_component(pressure_volume, error)
         if (allocated(error)) error stop error%message

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_lsf lsf = moist_new_svdw_lsf(error, NULL);
         moist_cavity cavity = moist_new_drop_cavity(error, NULL, lsf, NULL, NULL);
         moist_delete(lsf);     /* cavity owns a copy */
         moist_model model = moist_new_model(error, context, cavity);
         moist_delete(cavity);  /* model owns a copy */

         moist_component cpcm = moist_new_cpcm_component(error, NULL, 80.0, NULL);
         moist_component pv = moist_new_pv_component(error, NULL, 3.39893e-5);
         moist_add_model_component(error, model, cpcm);
         moist_add_model_component(error, model, pv);
         moist_delete(cpcm);    /* model owns copies */
         moist_delete(pv);

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import (
             Context, CavityDROP, CPCMRadii, ModelComponentCPCM,
             ModelComponentPV, SolvationModel, SvdW,
         )
         from moist.pyscf import GPA_TO_AU

         context = Context()
         cavity = CavityDROP(lsf=SvdW(), radii=CPCMRadii())
         model = SolvationModel(context, cavity, [
             ModelComponentCPCM(80.0),
             ModelComponentPV(1.0 * GPA_TO_AU),
         ])

The model runs on ``ctx``. The cavity and components above have no context of
their own, so their copies in the model use that; a part constructed with
``ctx=`` keeps its own.

After ``model%update(mol, error)``, create a coupling with ``model%new_coupling(coupling, error)`` and follow the :doc:`host coupling protocol </reference/coupling>` for building and evaluation.
``model%component_count()``, ``model%component_name(index, name)`` and ``model%component_description(index, description)`` describe the components before the first update; after ``get_energy``, each component's contribution is a named field (:doc:`/reference/fields`).
``model%print_parameters(unit)`` prints the settings at any time, at the context's :ref:`verbosity <verbosity>`: the cavity section followed by each component; ``unit`` defaults to the run context's unit.
``model%print_radii(unit)`` prints the radii table of the current structure once the model is updated, also before a deferred build; before the update it prints nothing.

Lifecycle
---------

The ``model_continuum_type`` owns one cavity and an ordered list of components:

1. Construct the model from a cavity and add all components before the first update; the model stores copies of both.
2. ``update`` and the :ref:`coupling-build` refresh the cavity, then its components.
3. ``get_energy`` sums and stores each component's contribution in its ``energy`` field (:doc:`/reference/fields`).
4. ``get_response`` collects each component's direct items (``potential_adjoint``, ``gaussian_amplitude``) and its surface weights, then lets the cavity contract the accumulated weights, which adds the ``density`` item for a cavity whose surface follows the density.
5. ``get_gradient`` adds the direct nuclear terms, contracts the gradient-side surface weights through the cavity, and refills the response with the same direct items.

Writing a Component
-------------------

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
