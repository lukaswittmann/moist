.. _python:

Python API
==========

MOIST requires Python 3.10 or newer and a compatible native MOIST library.
The Python API separates immutable settings from live native state. The same
configuration types work in standalone calculations and :doc:`pyscf`.

Construction
------------

Compose a DROP cavity from a level set, a radius model and discretization
parameters. Keep physical component inputs, such as dielectric constant and
pressure, separate from numerical settings:

.. code-block:: python

   from moist import (
       CavityDROP, CPCMRadii, DROPParameters, SvdW, SvdWParameters,
       ModelComponentCPCM, PCMParameters, PCMSolver, SolvationModel,
   )

   cavity = CavityDROP(
       lsf=SvdW(parameters=SvdWParameters(blend_k=5.5)),
       radii=CPCMRadii(),
       parameters=DROPParameters(nleb=194, proj_level=3),
   )
   model = SolvationModel(cavity, [
       ModelComponentCPCM(
           78.4, parameters=PCMParameters(solver=PCMSolver.CHOLESKY),
       ),
   ])

``CavityISwiG(parameters=ISwiGParameters(...), radii=...)`` has the same
construction convention. ``CavityDROP`` takes its level set as ``lsf=SvdW()``,
``lsf=CFC()`` or ``lsf=Isodensity(...)`` and all DROP level sets share
``DROPParameters``.
``SvdWParameters``, ``CFCParameters`` and ``IsodensityParameters`` configure
the surface independently. ``PCMParameters`` applies to both CPCM and COSMO;
``ModelParameters`` controls model logging.

Parameter fields correspond to the supported C options structs; defaults come
from the linked library's initializers and derived Fortran parameters remain
native. Parameter objects are immutable and keyword-only. Scientific
constraints are checked by the native constructor.
Inspect the settings used through ``cavity.parameters``,
``cavity.lsf.parameters``, ``component.parameters`` and ``model.parameters``.
``cavity.radius_model`` is the radius configuration; ``cavity.radii`` remains
the computed per-sphere radii after an update.

Radii and density sources
-------------------------

Available radii are ``CPCMRadii``, ``SMDRadii``, ``D3Radii``, ``COSMORadii``,
``BondiRadii`` and ``CustomRadii``. Omitted radii select CPCM. Custom radii are
in bohr, specified either in molecular atom order or by atomic number:

.. code-block:: python

   from moist import CustomRadii

   atom_radii = CustomRadii([2.67, 2.33, 2.33])
   element_radii = CustomRadii([2.33, 2.67], numbers=[1, 8])

An isodensity cavity takes its live density source separately from its LSF
settings. The source is a callable or an object with ``density(point, order)``;
it returns the bare density and its spatial derivatives, through the requested
order (1, 2 or 3). See :doc:`/cavities/isodensity` for the density convention.

.. code-block:: python

   from moist import Isodensity, IsodensityParameters

   cavity = CavityDROP(
       lsf=Isodensity(parameters=IsodensityParameters(rho_iso=1e-3)),
       parameters=DROPParameters(nleb=194, proj_maxiter=200),
       source=density_provider,
   )

The native default isodensity contour is ``1e-3`` electrons/bohr**3. This is
also the PySCF default. The callback and provider remain alive for as long as
the cavity or any model copy needs them.

Reusable configurations
-----------------------

``DROP`` and ``ISwiG`` describe a cavity without allocating live cavity state.
Build independent cavities with ``build()`` or pass the same configuration to
PySCF. Use ``dataclasses.replace`` to derive settings for a new calculation:

.. code-block:: python

   from dataclasses import replace
   from moist import DROP

   config = DROP(lsf=SvdW(), parameters=DROPParameters(nleb=194))
   finer = replace(config, parameters=replace(config.parameters, nleb=302))
   first = config.build()
   second = config.build()  # independent native state

Configurations and parameters can be compared, pickled, or converted to
dictionaries using ``dataclasses.asdict``. For reproducibility, record the
resolved settings and native library version. A live cavity's configuration
is read-only; create a new cavity/model to change it.

Array and accumulator contracts
-------------------------------

Python uses C-contiguous row-major arrays, with reversed native Fortran dimensions.
Positions and nuclear gradients have shape ``(natoms,3)``; grid vectors have shape
``(ngrid,3)``. Hessian weights have shape ``(ngrid,3,3)`` with indices
``[point,b,a]`` corresponding to native ``(a,b,point)``. These weights need not be
symmetric. Full diagnostic derivative tensors follow the same exact reversal.

``get_energy(coupling, energy)`` adds into a writable float64 scalar array such as
``np.array(0.0)``. ``get_gradient(coupling, gradient)`` adds into a writable,
C-contiguous float64 array of shape ``(natoms,3)`` and returns the host response.
Both leave their accumulator unchanged on failure. ``get_response(coupling)``
returns a fresh ``Response``.

Grid inputs are cavity properties: ``model.cavity.xyz``, ``xi0``, ``normal0``,
``a`` and ``f`` mirror ``model%cavity`` in Fortran, and ``cavity.get(name)``
reads any other field the cavity holds. ``Structure`` copies input arrays on
construction and update; mutating an input NumPy array cannot change the
stored geometry, call ``update`` explicitly. See :doc:`coupling` for the
protocol Python drives with ``for request in coupling``,
``coupling.answer(**outputs)`` and ``for item in response``, and
:ref:`coupling-requests` for every output's shape.

Internal Gaussian density
-------------------------

The same ``IsodensityParameters`` control both callback and native Gaussian
sources. Basis data and the live density are separate from those settings::

   from moist import GaussianBasis, InternalDensity, DROP, Isodensity

   basis = GaussianBasis(shell_atom=[0], shell_l=[0], shell_nprim=[1],
                         exponents=[0.5], coefficients=[1.0])
   source = InternalDensity(basis, [[1.0]])
   cavity = DROP(lsf=Isodensity()).build(source=source)
   shell_offsets, powers = cavity.isodensity_layout()

Shell atom indices and offsets are zero-based; ``powers`` has shape ``(ncart,3)``.
Coefficients include primitive normalization. The matrix must already be transformed
into the queried Cartesian-monomial layout; a host's spherical AO density cannot
be passed directly. Both basis and density inputs are copied.

Assign ``source.density_matrix`` before each cavity or model update. Updates
install it on the correct native owner; merely assigning it does not change an
already completed evaluation. Internal evaluation supplies spatial density
derivatives. The host still supplies its own basis-center derivative contractions.

``get_version_string()`` returns the linked library's full release string;
``get_banner(style)`` returns its banner without printing (``full``, ``short``,
``ascii`` or ``build``).

.. _coupling-python:

Host coupling
-------------

The shared :doc:`coupling` protocol defines the requests, response items and
phase lifecycle. Create a coupling after updating the model; the coupling
keeps its model alive.

The outline below shows the Python calls. ``host_potential`` and ``host_fock``
stand for host integral routines; complete the derivative request and
contraction loops for the chosen model. On a density-dependent cavity, retain
the response-phase ``DensityResponse`` before obtaining the gradient response
and use it for the additional density derivative contractions described in
:ref:`coupling-response`.

Python mirrors the same names. Iterating a coupling advances the native
cursor: each step yields a snapshot of the current request, and the loop ends
with the pass. ``coupling.answer`` takes the outputs of the current request
by keyword with the grid axis first. Iterating a response yields its items.
Grid inputs are read from ``model.cavity``.

.. code-block:: python

   import numpy as np
   from moist import GaussianPotentialRequest, PotentialAdjointResponse

   coupling = model.new_coupling()
   model.prepare_energy(coupling)
   for request in coupling:
       if isinstance(request, GaussianPotentialRequest):
           if "phi" in request.missing:
               coupling.answer(phi=host_potential(model.cavity.xyz, model.cavity.xi0))
       else:
           raise NotImplementedError(f"unsupported request {request.name}")
   energy = np.array(0.0)
   model.get_energy(coupling, energy)

   model.prepare_response(coupling)
   for request in coupling:
       ...  # walk again: a density-dependent cavity also needs dphi_dr, dphi_dxi
   response = model.get_response(coupling)
   for item in response:
       if isinstance(item, PotentialAdjointResponse):
           fock = host_fock(item.w_phi)
       else:
           raise NotImplementedError(f"unsupported response item {item.name}")

   model.prepare_gradient(coupling)
   for request in coupling:
       ...  # answer the missing derivative outputs the same way
   gradient = np.zeros((natoms, 3))
   response = model.get_gradient(coupling, gradient)
   for item in response:
       ...  # contract each item with the basis-center derivatives, as above

Each request is a snapshot taken when it is yielded: ``name`` is a class
attribute, ``missing`` the frozenset of output names missing at that moment,
and ``GaussianMomentRequest.width`` an immutable copy of the widths. An output
the request does not have is never in ``missing``.
``coupling.answer(**outputs)`` accepts any subset of the current request's
outputs, and ``None`` skips one. An unknown keyword raises ``TypeError`` and
no current request ``RuntimeError``, both before anything is submitted. The
rest go in declaration order, each validated on its own: the first rejection --
``ValueError`` for a wrong shape, ``RuntimeError`` for a native one -- leaves
that output missing, keeps the outputs submitted before it and skips the rest
of the call. A ``break`` leaves the pass for the next loop to resume;
a loop that ran to the end lets the next one start a new pass.

A ``Response`` is a plain value copied out of the native result: iterating it
yields the items the model produced, in native order, as
``PotentialAdjointResponse``, ``DensityResponse`` or
``GostshypAmplitudeResponse``. Each item has the native item name as its
``name`` class attribute and its arrays as attributes named as in the response
table, with the grid axis first. The :doc:`pyscf` module drives these loops for
PySCF.

API
---

.. automodule:: moist
   :members:
   :imported-members:

.. toctree::

   pyscf
