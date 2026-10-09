.. _fields:

Named fields
============

Named fields carry read-only data from MOIST to the host.
Everything MOIST hands to the host is read through fields; everything the host hands in is written through coupling answers (:doc:`coupling`).
The Fortran API provides direct access to the arrays; C and Python read them by name.

Five fields are available:

- the cavity: its grid and results;
- the model: the grid of its evaluation domain (a continuum model forwards its cavity's);
- a continuum component: its latest ``energy`` contribution in Hartree (0-based index in C/Python);
- the current response item: the arrays the host contracts, e.g. ``w_phi`` (:ref:`coupling-response`);
- the current coupling request: its inputs, e.g. the Gaussian moment ``width`` (:ref:`coupling-requests`).

Each field has a name, a type tag (``MOIST_FIELD_REAL``, ``MOIST_FIELD_INT``, ``MOIST_FIELD_BOOL``), a shape and a one-line description.

- A field that was not computed is not listed; reading it is an error, never zeros.
- The list follows the owner: list again after an update, and after each ``next()`` for response items and requests.
- Arrays have rank at most ``MOIST_FIELD_MAX_RANK`` (3).
- C and Python shapes reverse the Fortran shape, slowest-varying first: ``xyz(3,ngrid)`` is ``[ngrid][3]`` and ``w_hess_rho(3,3,ngrid)`` is ``[ngrid][3][3]`` (i.e. row-major).
- Atom and sphere indices such as ``owner`` are 0-based.
- Names have at most ``MOIST_FIELD_NAME_MAX`` characters.
- Response items, requests and components declare real fields only.
- Component ``energy`` appears after successful ``get_energy`` and clears on update. If a component fails, only earlier components' energies remain available.

.. list-table::
   :widths: 15 85

   * - Fortran
     - Override ``list_fields(self, query)`` and declare each array with ``query%add_*`` (``moist_channels_fields``).
       Response items and requests bind ``list_fields(query)`` next to their typed components (``field_query_type`` from ``moist``).
       ``model_continuum_type`` provides ``list_component_fields(index, query)``, ``component_count()``, ``component_name(index, name)`` and ``component_description(index, description)``; indices are 1-based.
   * - C
     - ``moist_get_cavity_field_*`` and ``moist_get_model_field_*``: ``_count``, ``_info``, ``_about``, ``_real``, ``_int``, ``_bool``.
       ``moist_get_response_field_*``, ``moist_get_coupling_request_field_*`` and ``moist_get_model_component_field_*``: ``_count``, ``_info``, ``_about``, ``_real``; component indices follow the model argument.
       ``moist_get_model_component_*``: ``_count``, ``_name``, ``_description`` (continuum models only).
       ``_info`` reports the element count a read writes; allocate ``MOIST_FIELD_NAME_MAX + 1`` characters for its name.
   * - Python
     - ``Cavity.fields()``, ``get(name)``, ``describe(name)`` and ``results(names)``.
       ``SolvationModel.components``: live views with ``name``, ``description``, ``energy``, ``configuration``, ``fields()``, ``get(name)`` and ``describe(name)``.
       Response items and requests are typed: their fields are attributes, e.g. ``item.w_phi`` and ``request.width``.

For cavity and component settings, use ``model%print_parameters(unit)``, ``moist_get_model_parameters_text`` or ``SolvationModel.parameters_text()``.
