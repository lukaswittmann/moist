Fortran API
===========

In MOIST, all used components (cavities, solvation model components) are based on abstract types.
For example, cavities extend ``cavity_type`` and all components extend ``solvation_model_component_type`` (both from ``moist_type``).
This means, only their constructors are type-specific and model usage and evaluation are independent of the selected cavity and components.

A shared ``moist_context_type`` controls output verbosity and carries a timer.
The context is created in the beginning and should outlive all components.
Errors use the allocatable ``error_type`` from ``mctc_env``; an allocated error signals failure.

.. code-block:: fortran

   use mctc_env, only : wp, error_type
   use moist_context, only : moist_context_type, new_context

   type(moist_context_type), target :: ctx
   type(error_type), allocatable :: error

   call new_context(ctx, verbosity=1)

Cavities
--------

Radii-based cavities need a radius model.
The examples below use CPCM radii, but ``moist_radii`` also provides, e.g. the SMD, COSMO, and Bondi radii.

.. code-block:: fortran

   use moist_radii, only : radius_type_static, new_cpcm_radii

   type(radius_type_static) :: radii

   call new_cpcm_radii(radii)

MOIST also supports custom atom-wise and element-wise radii:

.. code-block:: fortran

   use moist_radii, only : radius_type, new_radii_custom_atoms, &
      & new_radii_custom_elements

   class(radius_type), allocatable :: atom_radii, element_radii

   ! Atom-wise radii (molecular atom order, in bohr)
   call new_radii_custom_atoms([2.67_wp, 2.33_wp, 2.33_wp], &
      & atom_radii, error)
   if (allocated(error)) error stop error%message

   ! Element-wise radii (atomic numbers H=1 and O=8, in bohr)
   call new_radii_custom_elements([1, 8], [2.33_wp, 2.67_wp], &
      & element_radii, error)
   if (allocated(error)) error stop error%message

After a concrete cavity has been constructed, its use is independent of the
selected cavity type.  Call ``update`` whenever the molecular geometry changes.
For ρ-DROP, first install the new density on the cavity-owned LSF as shown
below, then call ``update`` to rebuild the surface.  ``get_gradient`` prepares
the cavity's nuclear-derivative arrays; it does not return a gradient argument.

.. code-block:: fortran

   ! Build the surface for the current structure
   call cavity%update(mol, error)
   if (allocated(error)) error stop error%message

   ! Common results are now available on every cavity
   write (*, *) "Grid points:", cavity%ngrid
   write (*, *) "Area / volume:", cavity%total_area, cavity%total_volume
   ! And can be printed directly via
   call cavity%print()

   ! Forward cavity quantity derivatives
   call cavity%get_gradient(error)
   if (allocated(error)) error stop error%message

   ! Optional diagnostics and write cavity files for inspection
   call cavity%find_disconnected_cavities(verbose_inp=2, error=error)
   if (allocated(error)) error stop error%message
   call cavity%write_xyz_debug("cavity.xyz", error)
   if (allocated(error)) error stop error%message
   call cavity%write_csv_debug("cavity.csv", error)
   if (allocated(error)) error stop error%message

iSwiG
~~~~~

The iSwiG constructor directly combines the radii and Lebedev discretization:

.. code-block:: fortran

   use moist_cavity_iswig, only : cavity_type_iswig, new_cavity_iswig

   type(cavity_type_iswig) :: cavity

   call new_cavity_iswig(cavity, ctx, nleb=194, cut_f=1.0e-10_wp, &
      & radius_model=radii, error=error)
   if (allocated(error)) error stop error%message

If a positive ``cut_a`` is supplied, iSwiG uses the area cutoff; otherwise it
uses ``cut_f``.
See :doc:`/cavities/iswig` for the supported Lebedev grids and defaults.

SvdW-DROP
~~~~~~~~~

DROP-based cavities separate the surface definition via a level set function (LSF) from the surface discretization (DROP).
The latter is general and independent of the choosen LSF.
That means, the LSF determines the cavity surface and the DROP scheme discretizes it.
All DROP variants accept the same optional discretization controls:
``tolerance``, ``proj_maxiter``, ``proj_level``, ``branch_weight_s``, ``rho_grid_h``, ``wleb_prune_level``, and ``do_fine``; see :doc:`/cavities/drop` for more detail.

For the Smooth van der Waals (SvdW) surface, construct the SvdW LSF and pass it into ``new_cavity_drop``:

.. code-block:: fortran

   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only : &
      & moist_cavity_drop_lsf_svdw_type

   type(moist_cavity_drop_lsf_svdw_type) :: svdw
   type(cavity_type_drop) :: cavity

   ! Construct SvdW LSF
   call svdw%new(blend_k=5.5_wp, blend_1b=1.0_wp, &
      & blend_2b=0.0_wp, blend_3b=3.0_wp)

   ! Construct SvdW-DROP
   call new_cavity_drop(cavity, ctx, nleb=194, radius_model=radii, &
      & lsf_model=svdw, error=error)
   if (allocated(error)) error stop error%message

See :doc:`/cavities/svdw` for the LSF parameters.

CFC-DROP
~~~~~~~~

For the COSMO Fine Cavity (CFC), only the constructor changes:

.. code-block:: fortran

   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_cfc, only : &
      & moist_cavity_drop_lsf_cfc_type

   type(moist_cavity_drop_lsf_cfc_type) :: cfc_lsf
   type(cavity_type_drop) :: cavity

   call cfc_lsf%new(a1=-15.0_wp, a2=-9.0_wp, c=5.0_wp, m=4)
   call new_cavity_drop(cavity, ctx, nleb=194, radius_model=radii, &
      & lsf_model=cfc_lsf, error=error)
   if (allocated(error)) error stop error%message

See :doc:`/cavities/cfc` for the CFC parameters.

Electron-isodensity (ρ-DROP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The internal isodensity LSF owns a Cartesian-monomial GTO basis and its current density matrix.
The host supplies the shell layout once and a new ``dcart`` for every SCF density:

.. code-block:: fortran

   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_isodensity_internal, only : &
      & moist_cavity_drop_lsf_isodensity_internal_type

   type(moist_cavity_drop_lsf_isodensity_internal_type) :: rho_lsf
   type(cavity_type_drop) :: cavity
   integer, allocatable :: sh_atom(:), sh_l(:), sh_nprim(:)
   real(wp), allocatable :: exps(:), coeffs(:), dcart(:, :)
   real(wp), parameter :: rho_iso = 1.0e-3_wp

   ! sh_atom, sh_l, sh_nprim, exps, coeffs, and dcart (from QM side)
   call rho_lsf%new(sh_atom, sh_l, sh_nprim, exps, coeffs, rho_iso, &
      & error=error)
   if (allocated(error)) error stop error%message

   ! Construct isodensity cavity
   call new_cavity_drop(cavity, ctx, nleb=194, radius_model=radii, &
      & lsf_model=rho_lsf, error=error)
   if (allocated(error)) error stop error%message

For every SCF step, the density has to be updated:

.. code-block:: fortran

   call cavity%rho_lsf%set_density(dcart, error)
   if (allocated(error)) error stop error%message

Alternatively, a callback-backed LSF calls a supplied function with the
``isodensity_lsf_callback`` interface:

.. code-block:: fortran

   use, intrinsic :: iso_c_binding, only : c_funloc, c_ptr, c_null_ptr
   use moist_cavity_drop_lsf_isodensity_callback, only : &
      & moist_cavity_drop_lsf_isodensity_callback_type, &
      & isodensity_lsf_callback

   type(moist_cavity_drop_lsf_isodensity_callback_type) :: callback_lsf
   type(c_ptr) :: context = c_null_ptr

   call callback_lsf%new(c_funloc(callback), context, rho_iso)

The callback returns the electron density and its requested spatial derivatives.
MOIST builds the level set ``S = scale * (rho_iso - rho)`` from it, so the callback must not subtract its own isovalue or flip its own sign.
See :doc:`/cavities/isodensity` for more details.

Model components
----------------

Components have different input data, but they are based on the ``solvation_model_component_type`` object.
Because some solvation model components couple with the QM density, some quantities have to be exchanged between MOIST and the QM side.

.. list-table:: Host exchange by component
   :header-rows: 1

   * - Component
     - Host data in ``coupling_type``
     - Main response in ``response_type``
   * - CPCM
     - ``electrostatics%qat`` for the point-charge mode, or ``electrostatics%phi`` and derivative data for a QM potential
     - ``electrostatics%surface_charge`` plus cavity-response (for isodensity cavity)
   * - PV
     - None
     - Cavity-response (for isodensity cavity)
   * - GOSTSHYP
     - ``gostshyp%gt``, ``gostshyp%pt``, ``gostshyp%mt``, and ``gostshyp%rt``
     - ``gostshyp%w_overlap``, ``gostshyp%w_normal_deriv``, and cavity-response (for isodensity cavity)

CPCM
~~~~

``new_component_cpcm`` constructs the conductor-like PCM component with
:math:`f(\varepsilon)=\dfrac{\varepsilon-1}{\varepsilon}`:

.. code-block:: fortran

   use moist_model_components, only : solvation_model_component_cpcm, &
      & new_component_cpcm, solver_type, potential_source

   type(solvation_model_component_cpcm) :: cpcm

   call new_component_cpcm(cpcm, ctx, epsilon=80.0_wp, &
      & solver=solver_type%cholesky, &
      & phi_source=potential_source%external, error=error)
   if (allocated(error)) error stop error%message

The default ``potential_source%charges`` computes the molecular potential from
``coupling%electrostatics%qat``.  A QM host should normally select
``potential_source%external``, evaluate the molecular cpcm potential
on the current cavity grid, and return it through ``coupling%electrostatics%phi``.

Available solvers are ``inversion``, ``lu``, ``cholesky`` (the default), and ``iterative``.

PV
~~

The PV component adds :math:`pV`, where ``V`` is the volume of the current
cavity:

.. code-block:: fortran

   use moist_model_components, only : solvation_model_component_pv, &
      & new_component_pv

   type(solvation_model_component_pv) :: pv

   call new_component_pv(pv, pressure=1.0_wp*gpa_to_au)

PV needs no host coupling data.  With ρ-DROP its surface response is returned
through the level set potential channels, making the density-dependent
``pV`` contribution variational.

GOSTSHYP
~~~~~~~~

GOSTSHYP places one Gaussian on each surface point and chooses its amplitude
to reproduce the requested pressure:

.. code-block:: fortran

   use moist_model_components, only : solvation_model_component_gostshyp, &
      & new_component_gostshyp

   type(solvation_model_component_gostshyp) :: gostshyp

   call new_component_gostshyp(gostshyp, pressure=50.0_wp*gpa_to_au)

After every cavity update, the host must rebuild all four Gaussian density
moments in ``coupling_type`` for the new grid.  ``get_response`` returns the
``w_overlap`` and ``w_normal_deriv`` amplitudes that the host contracts with its
Gaussian integral blocks to form the Fock contribution.

Using a component
~~~~~~~~~~~~~~~~~

Solvation model components can be used directly (manual way).
Every model component has the same ``update``, ``get_energy``, ``get_response``, and ``get_gradient`` routines.
The ``coupling_type`` contains data from the host to MOIST, while ``response_type`` carries the corresponding response from MOIST back to the host.
For the CPCM component below, a manual, direct-component implementation looks like this:

.. code-block:: fortran

   use moist_channels, only : coupling_type, response_type

   type(coupling_type) :: coupling
   type(response_type) :: response
   real(wp), allocatable :: phi(:), gradient(:, :)
   real(wp) :: energy

   ! Update the cavity first
   call cavity%update(mol, error)
   if (allocated(error)) error stop error%message
   ! Then every component with the updated cavity
   call cpcm%update(mol, cavity, error)
   if (allocated(error)) error stop error%message

   ! Host -> MOIST: evaluate the solute potential on cavity%xyz
   coupling%electrostatics%phi = host_cpcm_elstat_pot(cavity%xyz)

   ! Get CPCM energy contribution
   energy = 0.0_wp
   call cpcm%get_energy(coupling, cavity, energy, error)
   if (allocated(error)) error stop error%message

   ! Get Fock potential
   call cpcm%get_response(coupling, cavity, response, error)
   if (allocated(error)) error stop error%message
   ! response%electrostatics%surface_charge contains the surface charges

   ! Host -> MOIST: charge-weighted electronic potential gradient
   coupling%electrostatics%qefield(:, i) = &
      host_cpcm_deriv_elstat_pot(cpcm%q, cavity%xyz, ...)

   ! Get CPCM nuclear gradient contribution
   allocate(gradient(3, mol%nat), source=0.0_wp)
   call cpcm%get_gradient(coupling, cavity, gradient, error)
   if (allocated(error)) error stop error%message

The two ``host_cpcm_*`` calls are host routines which evaluate the potential and the nuclear derivatives of it on the current cavity.

Building a list-based model
---------------------------

A simpler approach is using the list-based model template.
This ``solvation_model_general`` owns the cavity and all added components. 

.. code-block:: fortran

   use mctc_env, only : wp, error_type
   use moist_context, only : moist_context_type, new_context
   use moist_radii, only : radius_type_static, new_cpcm_radii
   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only : &
      & moist_cavity_drop_lsf_svdw_type
   use moist_model_components, only : solvation_model_component_cpcm, &
      & new_component_cpcm, solvation_model_component_pv, new_component_pv, &
      & solvation_model_component_gostshyp, new_component_gostshyp, &
      & potential_source
   use moist_model_general, only : solvation_model_general, new_model_general
   
   type(error_type), allocatable :: error

   type(moist_context_type), target :: ctx
   ! Cavity
   type(radius_type_static) :: radii
   type(moist_cavity_drop_lsf_svdw_type) :: svdw
   type(cavity_type_drop) :: cavity
   ! Model components
   type(solvation_model_component_cpcm) :: cpcm
   type(solvation_model_component_pv) :: pv
   type(solvation_model_component_gostshyp) :: gostshyp
   ! Model
   type(solvation_model_general) :: model

   ! Conversion of GPa to a.u.
   real(wp), parameter :: gpa_to_au = 3.3989e-5_wp

   call new_context(ctx)

   ! Construct cavity
   call new_cpcm_radii(radii)
   call svdw%new()
   call new_cavity_drop(cavity, ctx, nleb=194, radius_model=radii, &
      & lsf_model=svdw, error=error)
   if (allocated(error)) error stop error%message

   ! Construct all model components
   call new_component_cpcm(cpcm, ctx, epsilon=80.0_wp, &
      & phi_source=potential_source%external, error=error)
   if (allocated(error)) error stop error%message
   call new_component_pv(pv, pressure=10.0_wp*gpa_to_au)
   call new_component_gostshyp(gostshyp, pressure=50.0_wp*gpa_to_au)

   ! Construct model and add all components
   call new_model_general(model, cavity, ctx, error)
   if (allocated(error)) error stop error%message
   call model%add_component(cpcm, error)
   if (allocated(error)) error stop error%message
   call model%add_component(pv, error)
   if (allocated(error)) error stop error%message
   call model%add_component(gostshyp, error)
   if (allocated(error)) error stop error%message

The model usage is similar to the components:
update the model for the current geometry ``model%update`` with an initialized ``structure_type``; 
host fills all needed ``coupling_type`` for the current ``model%cavity``;
get energy, response and gradient via ``model%get_energy``, ``model%get_response``, and ``model%get_gradient``.
