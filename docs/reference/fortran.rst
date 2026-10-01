Fortran API
===========

Cavity, LSF, and PCM constructors accept an optional ``param`` object.
Omitting it uses compiled defaults; constructors copy supplied values.
Required inputs such as radii, LSFs, context, and dielectric constant remain
separate arguments. Individual setting keywords are not accepted.

Cavity, LSF and PCM parameter types extend ``moist_model_parameters_type`` and
support ``read_file(path, error)``, ``write_file(path, error)`` and
``print_parameters(error, unit=...)``. The file extension selects JSON
(``.json``) or TOML (``.toml``), ignoring case; any other extension is an
error. Missing fields keep their defaults. ``print_parameters`` prints JSON.
DROP recomputes its derived settings both on ``read_file`` and in the cavity
constructor.

``solvation_system_type`` holds physical properties and molecular structures.
Create it with ``new_solvation_system`` from ``moist_data_solvents``; it does not
extend the configuration parameter base.

.. code-block:: fortran

   use moist, only: moist_cavity_drop_parameters_type, new_cavity_drop

   type(moist_cavity_drop_parameters_type) :: param

   call param%read_file("drop.json", error)
   if (allocated(error)) error stop error%message
   param%num_leb = 194
   call new_cavity_drop(cavity, ctx, radius_model=radii, lsf_model=lsf, &
      & error=error, param=param)
   if (allocated(error)) error stop error%message

A new defaulted parameter field keeps existing constructor calls valid, but a
Fortran host must recompile whenever module definitions change.

Use the installed compiler-specific ``.mod`` files. The Meson install provides
only MOIST's module files. Hosts that call mctc-lib directly must link a
compatible mctc-lib, its dependencies and its module files themselves;
``pkg-config moist`` supplies none of them.

Cavities extend ``cavity_type``; continuum components extend
``model_continuum_component_type`` (modules ``moist_cavity_type`` and
``moist_model_continuum_component_type``). Constructors are specific to each type; evaluation uses
the shared interfaces.

``moist_context_type`` controls logging and timing and must outlive its
components. An allocated ``mctc_env`` error signals failure.

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

Call ``update`` after geometry changes. For ρ-DROP, first set the new density
on the cavity-owned LSF. ``get_gradient`` prepares nuclear-derivative arrays.

.. code-block:: fortran

   use moist_cavity_diagnostic, only: find_disconnected_cavities

   integer, allocatable :: island_sizes(:)

   ! Build the surface for the current structure
   call cavity%update(mol, error)
   if (allocated(error)) error stop error%message

   ! Common results are now available on every cavity
   write (*, *) "Grid points:", cavity%ngrid
   write (*, *) "Area / volume:", cavity%total_area, cavity%total_volume
   call cavity%print()

   ! Forward cavity quantity derivatives
   call cavity%get_gradient(error)
   if (allocated(error)) error stop error%message

   ! Optional diagnostics and write cavity files for inspection
   ! find_disconnected_cavities returns the point count of each island
   call find_disconnected_cavities(cavity, island_sizes, error)
   if (allocated(error)) error stop error%message
   if (size(island_sizes) > 1) print "(a,i0,a)", &
      & "grid has ", size(island_sizes), " islands"
   call cavity%write_xyz_debug("cavity.xyz", error)
   if (allocated(error)) error stop error%message
   call cavity%write_csv_debug("cavity.csv", error)
   if (allocated(error)) error stop error%message

iSwiG
~~~~~

The iSwiG constructor directly combines the radii and Lebedev discretization:

.. code-block:: fortran

   use moist_cavity_iswig, only : cavity_type_iswig, new_cavity_iswig, &
      & moist_cavity_iswig_parameters_type

   type(cavity_type_iswig) :: cavity

   call new_cavity_iswig(cavity, ctx, radius_model=radii, error=error, &
      & param=moist_cavity_iswig_parameters_type(num_leb=194, cut_f=1.0e-10_wp))
   if (allocated(error)) error stop error%message

A positive ``cut_a`` discards points with switched area ``f*a <= cut_a``;
otherwise the gate is ``f <= cut_f``.
See :doc:`/cavities/iswig` for the supported Lebedev grids and defaults.

SvdW-DROP
~~~~~~~~~

A level-set function (LSF) defines the surface; DROP discretizes it. All DROP
variants share the controls documented in :doc:`/cavities/drop`.
For SvdW, pass its LSF to ``new_cavity_drop``:

.. code-block:: fortran

   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only : &
      & moist_cavity_drop_lsf_svdw_type

   type(moist_cavity_drop_lsf_svdw_type) :: svdw
   type(cavity_type_drop) :: cavity

   ! Construct SvdW LSF
   call svdw%new()

   ! Construct SvdW-DROP
   call new_cavity_drop(cavity, ctx, radius_model=radii, &
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

   call cfc_lsf%new()
   call new_cavity_drop(cavity, ctx, radius_model=radii, &
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

   use moist, only: moist_cavity_drop_lsf_isodensity_param_type

   type(moist_cavity_drop_lsf_isodensity_internal_type) :: rho_lsf
   type(cavity_type_drop) :: cavity
   integer, allocatable :: sh_atom(:), sh_l(:), sh_nprim(:)
   real(wp), allocatable :: exps(:), coeffs(:), dcart(:, :)
   real(wp), parameter :: rho_iso = 1.0e-3_wp

   ! sh_atom, sh_l, sh_nprim, exps, coeffs, and dcart (from QM side)
   call rho_lsf%new(sh_atom, sh_l, sh_nprim, exps, coeffs, error=error, &
      & param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=rho_iso))
   if (allocated(error)) error stop error%message

   ! Construct isodensity cavity
   call new_cavity_drop(cavity, ctx, radius_model=radii, &
      & lsf_model=rho_lsf, error=error)
   if (allocated(error)) error stop error%message

Set the density at each SCF step:

.. code-block:: fortran

   select type (lsf => cavity%lsf_model)
   type is (moist_cavity_drop_lsf_isodensity_internal_type)
      call lsf%set_density(dcart, error)
   end select
   if (allocated(error)) error stop error%message

Alternatively, a callback-backed LSF calls a supplied function with the
``isodensity_lsf_callback`` interface:

.. code-block:: fortran

   use, intrinsic :: iso_c_binding, only : c_funloc, c_ptr, c_null_ptr
   use moist_cavity_drop_lsf_isodensity_callback, only : &
      & moist_cavity_drop_lsf_isodensity_callback_type, &
      & isodensity_lsf_callback

   use moist, only: moist_cavity_drop_lsf_isodensity_param_type

   type(moist_cavity_drop_lsf_isodensity_callback_type) :: callback_lsf
   type(c_ptr) :: context = c_null_ptr

   call callback_lsf%new(c_funloc(callback), context, &
      & param=moist_cavity_drop_lsf_isodensity_param_type(rho_iso=rho_iso))

The callback returns the electron density and its requested spatial derivatives.
MOIST builds the level set ``S = scale * (rho_iso - rho)`` from it, so the callback must not subtract its own isovalue or flip its own sign.
See :doc:`/cavities/isodensity` for more details.

Model components
----------------

Components declare host quantities through :doc:`coupling` and return
coefficients for the host's integrals.

.. list-table:: Host exchange by component
   :header-rows: 1

   * - Component
     - Requests it declares
     - Response items it produces
   * - CPCM, COSMO
     - ``gaussian_potential`` with independently required ``phi``,
       ``dphi_dr`` and ``dphi_dxi`` outputs
     - ``potential_adjoint`` (``w_phi``, :math:`\mathrm{d}E/\mathrm{d}\phi`;
       the induced charges for a stationary PCM)
   * - PV
     - None
     - None of its own
   * - GOSTSHYP
     - ``gaussian_moments`` (``gt``, ``pt``, ``mt``, ``rt``)
     - ``gaussian_amplitude`` (``w_overlap``, ``w_normal_deriv``)

On a density-backed cavity (the isodensity level sets) the response of every
model additionally carries the ``density`` item, which belongs to the cavity
rather than to any one component. Which derivative outputs the response phase
requires also depends on the cavity; see :ref:`coupling-requests`.

CPCM
~~~~

``new_component_cpcm`` constructs the conductor-like PCM component with
:math:`f(\varepsilon)=\dfrac{\varepsilon-1}{\varepsilon}`:

.. code-block:: fortran

   use moist_model_continuum_component, only : model_continuum_component_cpcm, &
      & new_component_cpcm, solver_type, moist_pcm_parameters_type

   type(model_continuum_component_cpcm) :: cpcm

   call new_component_cpcm(cpcm, ctx, epsilon=80.0_wp, &
      & error=error, param=moist_pcm_parameters_type(solver=solver_type%cholesky))
   if (allocated(error)) error stop error%message

The host answers ``gaussian_potential`` at the grid points ``cavity%xyz`` with
the widths ``cavity%xi0``; see :ref:`coupling-requests` for the operator.

Available solvers are ``inversion``, ``lu``, ``cholesky`` (the default), and
``iterative``.

COSMO
~~~~~

``new_component_cosmo`` uses
:math:`f(\varepsilon)=\dfrac{\varepsilon-1}{\varepsilon+1/2}`
with the same host requests and solver options as CPCM:

.. code-block:: fortran

   use moist_model_continuum_component, only : model_continuum_component_cosmo, &
      & new_component_cosmo, solver_type, moist_pcm_parameters_type

   type(model_continuum_component_cosmo) :: cosmo

   call new_component_cosmo(cosmo, ctx, epsilon=80.0_wp, &
      & error=error, param=moist_pcm_parameters_type(solver=solver_type%cholesky))
   if (allocated(error)) error stop error%message

PV
~~

The PV component adds :math:`pV`, where ``V`` is the volume of the current
cavity:

.. code-block:: fortran

   use moist_model_continuum_component, only : model_continuum_component_pv, &
      & new_component_pv

   real(wp), parameter :: gpa_to_au = 3.39893e-5_wp
   type(model_continuum_component_pv) :: pv

   call new_component_pv(pv, pressure=1.0_wp*gpa_to_au)

PV declares no host requests. With ρ-DROP, its Fock contribution uses the
cavity's ``density`` response.

.. _gostshyp-fortran:

GOSTSHYP
~~~~~~~~

:doc:`GOSTSHYP </models/components/gostshyp>` places one Gaussian on each surface point and chooses its amplitude
to reproduce the requested pressure:

.. code-block:: fortran

   use moist_model_continuum_component, only : model_continuum_component_gostshyp, &
      & new_component_gostshyp

   real(wp), parameter :: gpa_to_au = 3.39893e-5_wp
   type(model_continuum_component_gostshyp) :: gostshyp

   call new_component_gostshyp(gostshyp, pressure=50.0_wp*gpa_to_au)

The ``gaussian_moments`` request (``gaussian_moment_request_type``) carries
the exponents in ``request%width``; read them from the request, never
recompute them. Answer only its missing outputs, and contract ``w_overlap``
and ``w_normal_deriv`` of the ``gaussian_amplitude`` item
(``gaussian_amplitude_response_type``) with the host's Gaussian integrals.
See :ref:`coupling-requests` for the outputs required per phase.

Building a list-based model
---------------------------

``model_continuum_type`` owns its cavity and components. Using the objects
constructed above:

.. code-block:: fortran

   use moist_model_continuum, only : model_continuum_type, new_continuum_model

   type(model_continuum_type), target :: model

   call new_continuum_model(model, cavity, ctx, error)
   if (allocated(error)) error stop error%message
   call model%add_component(cpcm, error)
   if (allocated(error)) error stop error%message
   call model%add_component(pv, error)
   if (allocated(error)) error stop error%message
   call model%add_component(gostshyp, error)
   if (allocated(error)) error stop error%message

Add all components before the first ``model%update``. After a successful update,
create a coupling with ``model%new_coupling``. For each phase, call
``prepare_*``, answer the missing outputs of each request that
``do while (coupling%next())`` makes current, call ``get_*``, and contract
each item that ``do while (response%next())`` makes current. Grid inputs are
the public arrays of ``model%cavity`` (``xyz``, ``xi0``, ``normal0``, ``a``,
``f``). See :ref:`coupling-fortran` for the evaluation loop and ownership
rules.

The ``moist`` umbrella module re-exports what the evaluation loop needs: every
parameter type, ``cavity_type_drop`` and ``new_cavity_drop``, ``radius_type``
and the generic ``new_radii``, the model and component types with their
constructors and ``solver_type``, the coupling and request types, the response
type with its items, ``field_query_type``, and ``wp``, ``error_type``,
``fatal_error`` and ``structure_type`` from mctc-lib. Everything else the
examples above use -- the context, the other cavities, the concrete radii and
LSF types, the diagnostics -- comes from its own module.

``coupling_type`` and ``response_type`` carry only the host bindings
(``next``, ``request``, ``answer``; ``next``, ``item``); declaring, staging and
accumulating belong to the model and are module procedures of
``moist_channels_coupling`` and ``moist_channels_response``, which ``moist``
does not re-export.

.. _coupling-fortran:

Host coupling
-------------

The shared :doc:`coupling` protocol defines the requests, response items and
phase lifecycle. The example below shows energy and Fock evaluation for PCM
on a cavity fixed with respect to the density. ``host_potential`` and
``host_fock`` stand for the host's own integral routines.

The model must have ``target`` and outlive its couplings. Create a coupling
after a successful model update; ``model%release_coupling(coupling)`` releases
it early. Grid inputs are read from ``model%cavity``.

``coupling%request()`` returns a copy of the current request carrying its
kind, outputs, requirement state and inputs such as the moment width. Later
answers do not update that copy. Answers go into the coupling by output name
and always address the current request.
The response is walked the same way, one item copy at a time.

.. code-block:: fortran

   use moist

   type(model_continuum_type), target :: model
   type(coupling_type), pointer :: coupling
   type(response_type) :: response
   type(error_type), allocatable :: error
   real(wp), allocatable :: phi(:)
   real(wp) :: energy

   ! Construct and update model before entering this loop.
   ! Check error after each library call in production code.
   call model%new_coupling(coupling, error)
   call model%prepare_energy(coupling, error)
   do while (coupling%next())
      select type (request => coupling%request())
      type is (gaussian_potential_request_type)
         if (request%is_missing("phi")) then
            call host_potential(model%cavity%xyz, model%cavity%xi0, phi)
            call coupling%answer("phi", phi, error)
         end if
      class default
         error stop "unsupported request "//request%name()
      end select
   end do
   energy = 0.0_wp
   call model%get_energy(coupling, energy, error)

   call model%prepare_response(coupling, error)
   ! Walk the requests again as above: on a fixed cavity phi is kept and
   ! nothing is missing; a density-dependent cavity also needs dphi_dr, dphi_dxi.
   call model%get_response(coupling, response, error)
   do while (response%next())
      select type (item => response%item())
      type is (potential_adjoint_response_type)
         call host_fock(item%w_phi)
      class default
         error stop "unsupported response item "//item%name()
      end select
   end do

   call model%release_coupling(coupling)

Both ``next()`` functions are logical functions with a side effect: keep each
alone in its loop condition, because Fortran does not guarantee that every
operand of ``.and.`` or ``.or.`` is evaluated. ``request()`` takes no argument
and returns an allocatable ``class(coupling_request_type)`` copy; while no
request is current it returns a placeholder named ``no_current_request``
without outputs, which falls through to ``class default``.
``is_missing(name)`` queries the copy.
``answer(name, values, error)`` takes the output name and an array of
the output's declared rank; an unknown name, a wrong rank or shape, a
non-finite value or a missing current request comes back as an error.
``response%item()`` likewise returns a polymorphic copy of the current item,
with its arrays as components (``item%w_phi``), and the placeholder
``no_current_item`` outside a walk.
Items and requests also bind ``list_fields(query)``, which declares the same
arrays as :doc:`fields` to a ``field_query_type``: after ``query%enumerate()``
it fills ``query%info``, after ``query%fetch(name)`` it copies that one array
into ``query%rvals``. The typed components remain the way to read them.
``get_gradient(coupling, response, gradient, error)`` fills the response with
the host part of the gradient phase, walked the same way, and adds to
``gradient(3, nat)``.

On a density-dependent cavity, retain the response-phase density weights
before calling ``get_gradient``, which replaces the response items. Use those
weights for the additional host contractions described in :ref:`coupling-response`.
