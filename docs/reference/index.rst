.. _api:

API Documentation
=================

This page introduces those objects and the rules the APIs share; the pages listed at the end describe each part in detail.
The :doc:`PySCF integration <pyscf>` contains details for PySCF.

Building blocks
---------------

Context
   The run settings every other object shares: the thread count, how much MOIST prints, and the timer recording where the time goes (:doc:`context`).

Structure
   The solute's elements and geometry. Nuclear charges arrive through the coupling (:ref:`coupling-requests`).

Cavity
   The surface separating the solute from the solvent, discretized into grid points with areas, normals and Gaussian widths.
   A radius model gives each atom a sphere (:ref:`cavity-radii`); a DROP cavity also takes a level-set function (LSF), which defines either a smooth surface around those spheres or an electron-isodensity contour (:doc:`/cavities/index`).
   The cavity is rebuilt for every geometry and, when it follows the electron density, for every density.

Component
   One energy term evaluated on the cavity, such as the dielectric polarization of CPCM or COSMO, the pressure-volume term PV or the GOSTSHYP pressure (:doc:`/models/components/index`).

Model
   The object the host drives: a cavity together with its components (:doc:`/models/index`).
   The molecular Ornstein-Zernike models (:doc:`/models/moz`) are models too and work on their own grids instead of a cavity.

Coupling and response
   A coupling lists what a model needs from the host, such as the electrostatic potential on the grid, and holds the host's answers.
   A response carries what the model returns: coefficients the host contracts into its Fock matrix and nuclear gradient (:doc:`coupling`).

Fields
   Named read-only arrays, such as the grid points of a cavity or the energy of a component (:doc:`fields`).

A calculation
-------------

.. code-block:: text

   create a context
   create the cavity (radius model, level set, settings) and the components
   create the model from the context and the cavity, then add the components
   print the settings                                   (optional)

   for each geometry:
       update the model with the structure
       first geometry only: create a coupling (each update resets it)
       for each SCF iteration:
           build phase with the current density
           energy and response phases    -> energy, Fock contribution
       gradient phase                    -> nuclear gradient

   print the timings                                    (optional)

:ref:`coupling-build` describes when rebuilding is required; :doc:`coupling` covers evaluation.

Run context
-----------

A calculation creates its context first and passes it to the model; a cavity or component without a context of its own runs on the model's.
The thread count is fixed when the context is created, and ``0`` takes the calling thread's OpenMP setting; MOIST never changes the host's OpenMP or BLAS threading.
The verbosity runs from 0 (silent) to 6 (debug), and the default 2 suits a production calculation.
MOIST prints its settings and timings only when the host asks.
:doc:`context` describes the threads, the output levels, the settings printout and the timer.

Settings
--------

Constructors take the physical inputs as arguments, such as the dielectric constant, the pressure, the radius model and the level set.
The numerical settings of each configurable cavity, level set and component go into one parameter object; any setting the host leaves out keeps the library default.

- Fortran: the optional ``param=`` argument, a type extending ``moist_model_parameters_type`` that also reads and writes JSON and TOML files.
- C: an options struct, filled by its ``moist_init_*_options`` initializer before the host assigns fields; ``NULL`` selects the defaults.
- Python: an immutable, keyword-only parameters class, such as ``DROPParameters`` or ``PCMParameters``.

The settings printout (:doc:`context`) shows the values a model runs with.

.. _api-solvents:

Solvents
--------

MOIST carries a table of 180 solvents.
A lookup takes a name or a a moist internal solvent ID.
Names ignore case and surrounding whitespace and accept various aliases, such as ``H2O``, ``water``, ..., or a CAS number.
IDs are frozen: an ID never changes meaning, and new solvents are appended.
A lookup needs no context and reads only constant tables, so any thread may call it; an unknown name or ID is an error.

Solvent descriptors
   The result of a lookup: one table entry with the same fields in every language.
   Quantities with units come twice, in SI (``_si``) and in atomic units (``_au``).

   .. list-table::
      :header-rows: 1

      * - Field
        - Meaning
        - Units
      * - ``id``, ``name``
        - Frozen table ID; canonical name, lowercase
        -
      * - ``epsilon``
        - Static relative permittivity
        -
      * - ``refractive_index``
        - Refractive index at the sodium D line
        -
      * - ``abraham_alpha``, ``abraham_beta``
        - Abraham hydrogen-bond acidity and basicity
        -
      * - ``surface_tension``
        - Macroscopic surface tension
        - N/m; Eh/bohr²
      * - ``mass_density``
        - Mass density
        - kg/m³; mₑ/bohr³
      * - ``molar_mass_si``, ``molecular_mass_au``
        - Mass of one mole; of one molecule
        - kg/mol; mₑ
      * - ``number_density``
        - Molecules per volume
        - 1/m³; 1/bohr³
      * - ``molecular_volume``
        - Volume per molecule
        - m³; bohr³

   The table follows the Minnesota Solvent Descriptor Database: the permittivity and surface tension at 298.15 K, the refractive index at 293.15 K and the density mostly at 293.15 K.
   The last three rows derive from the solvent geometry MOIST stores.
   A NaN marks a value the table lacks and carries into the values derived from it, while ``0`` is a real value, for example the acidity of an aprotic solvent.
   Every table value is present; decalin has no stored geometry, so its derived values are NaN.
   Each language also formats descriptors as a ``Solvent properties`` printout of the values passed, without a new lookup.

Solvent system
   Fortran only: the descriptors as its ``solvent`` component, together with a temperature, a pressure and the solvent geometry, as the MOZ models use it (:doc:`/models/moz`).
   Building one fails without a stored geometry.
   The descriptors stay at their reference temperatures and do not follow the system's temperature.

Ownership and lifetime
----------------------

Constructors copy what they are given: the model holds its own cavity and components, the cavity its own level set, and an update keeps its own copy of the structure.
The originals can therefore be deleted or reused.
The context is shared: Fortran callers keep it alive; C and Python objects retain it.

Add every component before the first update; later additions are rejected.

A coupling belongs to the model that created it and can be created once the model has been updated.
Several couplings may coexist, and every model update resets them.
The model must outlive its couplings:

- Fortran: the model has the ``target`` attribute; ``model%release_coupling`` releases a coupling early.
- C: delete the couplings before their model.
- Python: a coupling keeps its model alive.

A host callback, such as an isodensity callback, must stay valid as long as any cavity or model uses it; Python keeps the callable alive itself.

Density-dependent cavities
--------------------------

An isodensity cavity (ρ-DROP) follows the host's electron density.
Rebuild it at each density change before evaluating the energy (:ref:`coupling-build`).
MOIST obtains the density in one of two ways (:doc:`/cavities/isodensity`):

Internal
   The host supplies its Gaussian basis at construction and a density matrix at each SCF step.
   The matrix always uses MOIST's Cartesian-monomial layout, including with the ``spherical`` basis option.
Callback
   The host evaluates the density and its spatial derivatives at the points MOIST asks for.
   It returns the density itself; MOIST subtracts the isovalue.

On such a cavity the response phase requires more derivative outputs (:ref:`coupling-requests`) and also returns a ``density`` item.
Its weights, contracted with the density-matrix derivative of the density, complete the Fock matrix, and contracted with the basis-center derivatives, the nuclear gradient.
The gradient phase does not return this item, so the host keeps the weights from the response phase at the converged density (:ref:`coupling-response`).

MOIST uses the host's nuclear charges to identify ECP atoms and fills their missing core density to prevent holes in the surface.

Conventions
-----------

MOIST works in atomic units: lengths in bohr, energies in Hartree, pressures in Hartree/bohr³.
Fortran arrays keep their native shapes, for example positions ``(3, natoms)``; C and Python reverse the dimensions without moving data, so positions are ``[natoms][3]`` and ``(natoms, 3)``.

Calls that can fail report it to the host: Fortran through an allocatable ``type(error_type)`` argument, C through the ``moist_error`` handle every such call takes first, Python by raising an exception.
Energy and gradient evaluations add into accumulators the host owns; a failed evaluation leaves them unchanged.

Where to go next
----------------

- :doc:`context`: threads, output levels, settings printout and timings.
- :doc:`coupling`: the request-and-response exchange of every SCF step and gradient.
- :doc:`fields`: reading cavity, model, component, request and response arrays by name.
- :doc:`fortran`, :doc:`c` and :doc:`python`: constructors, ownership and a complete evaluation loop in each language.
- :doc:`grids`: the radial, angular and volume integration grids inside MOIST, for Fortran developers; the host loop does not need them.

.. toctree::

   context
   coupling
   fields
   fortran
   c
   python
   grids
