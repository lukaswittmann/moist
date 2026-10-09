Cavities
========

This section describes the cavity constructions available in MOIST and the settings that control their shape and surface discretization.

Common Interface
----------------

All cavity implementations extend the abstract ``cavity_type`` in ``src/moist/cavity/type.f90``.
It carries the shared allocatable quantities that concrete cavities fill where applicable:

- atomic sphere centers, radii, and per-sphere areas;
- grid point positions, owners, areas, outward unit normals, Gaussian widths,
  switching function, and volume elements;
- total cavity area (bohr\ :sup:`2`) and enclosed volume (bohr\ :sup:`3`);
- optional nuclear derivatives of positions, widths, switching function, and
  volume elements.

Concrete cavity types must implement two deferred procedures:

``update(mol, error)``
   Build or refresh the cavity for the supplied molecular structure and fill
   the common fields supported by that implementation.

``get_gradient(error)``
   Prepare or validate the implementation-specific nuclear derivatives.

The base type also provides optional response hooks; ``get_surface_response``
maps model surface weights to host-response contributions.

Typical Use
-----------

1. Construct a concrete cavity and configure its radii and discretization.
2. Call ``update`` whenever its inputs change: after a geometry change, and for
   a density-backed surface after a density change.
3. A ``model_continuum_type`` owns one cavity, updates it before its
   components, and passes the same cavity to each component.
4. In the response phase, components accumulate model surface weights and
   ``get_surface_response`` maps them to the response items the host receives
   (the ``density`` item, for a cavity whose surface follows the density).
5. Call ``get_gradient`` before requesting supported nuclear derivatives.

The shared ``print``, ``write_xyz_debug``, ``write_csv_debug``, and ``write_pqr_debug`` procedures provide diagnostics and grid export.
``find_disconnected_cavities`` from ``moist_cavity_diagnostic`` is a separate diagnostic: it returns the point count of each connected island, largest first, with a connectivity radius of ``threshold`` times the mean nearest-neighbour spacing (4.0 when omitted). No cavity runs it during ``update``.

.. _cavity-radii:

Radius Models
-------------

Every cavity needs a radius model. ``CPCM`` radii are the default; ``SMD``,
``COSMO``, ``Bondi``, ``D3`` and custom per-atom or per-element radii are also
available. The cavity constructor copies the model.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use moist_radii, only : radius_type_static, new_cpcm_radii

         type(radius_type_static) :: radii

         call new_cpcm_radii(radii)

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_radii radii = moist_new_cpcm_radii(error);
         /* cavity constructors copy it; moist_delete(radii) afterwards.
            Passing NULL radii selects CPCM radii. */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import CPCMRadii

         radii = CPCMRadii()

Per-atom and per-element radii are supported too, as ``CustomRadii`` in Python,
``moist_new_custom_radii`` plus its setters in C, and ``new_radii_custom_atoms`` /
``new_radii_custom_elements`` in Fortran; see :doc:`/reference/fortran`.

Available Cavities
------------------

.. toctree::

   drop
   iswig
   numsa
   marchingcubes
