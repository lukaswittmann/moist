.. _api:

API Documentation
=================

MOIST supports Fortran, C and Python as APIs.

The request/response exchange a QM host drives is the same in all three and is explained in :doc:`coupling`, from the division of responsibility and energy derivatives to the SCF workflow and detailed protocol.
The MOIST arrays can be accessed via fields that are described in :doc:`fields`.
Each language page provides its own evaluation examples and API conventions.
The PySCF integration is documented under :doc:`python`.

.. toctree::

   coupling
   fields
   fortran
   c
   python
