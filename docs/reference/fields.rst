.. _fields:

Named fields
============

Named fields carry read-only data from MOIST to the host.
The Fortran API provides direct access the arrays; C and Python read them by name.

Each field has a name, a type tag (``MOIST_FIELD_REAL``, ``MOIST_FIELD_INT``, ``MOIST_FIELD_BOOL``), a shape and a one-line description.

- A field that was not computed is not listed; reading it is an error, never zeros.
- The list follows the owner: list again after an update.
- C and Python shapes reverse the Fortran shape: ``xyz(3,ngrid)`` is ``[ngrid][3]`` (i.e. row-major).
- Atom and sphere indices such as ``owner`` are 0-based.
- Names have at most ``MOIST_FIELD_NAME_MAX`` characters.

.. list-table::
   :widths: 15 85

   * - Fortran
     - Override ``list_fields(self, query)`` and declare each array with ``query%add_*`` (``moist_channels_fields``).
   * - C
     - ``moist_get_cavity_field_*`` and ``moist_get_model_field_*``: ``_count``, ``_info``, ``_about``, ``_real``, ``_int``, ``_bool``.
       ``_info`` reports the element count a read writes; allocate ``MOIST_FIELD_NAME_MAX + 1`` characters for its name.
   * - Python
     - ``Cavity.fields()``, ``get(name)``, ``describe(name)`` and ``results(names)``.
