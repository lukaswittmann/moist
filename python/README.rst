MOIST Python API
================

This package provides access to the ``moist`` C API through CFFI and requires Python 3.10 or newer.
The low-level bindings live in ``moist.library``, the object model in ``moist.interface``, and the PySCF integration in ``moist.pyscf``.
See the documentation under ``docs/reference`` for usage.

Building
--------

The package is normally built with the library, by configuring the top-level project with ``-Dpython=true``.
The extension module and the pure-Python modules then land in ``<build>/python/moist``, and ``PYTHONPATH=<build>/python`` imports the package.

This directory is also a standalone meson-python project for building a self-contained wheel. It expects the moist source tree in ``subprojects/moist`` (or an installed moist found through ``PKG_CONFIG_PATH``); see the ``wheel`` job in ``.github/workflows`` for the exact steps.
Options of the moist subproject are passed through, for example ``pip install ./python -Csetup-args=-Dlapack=accelerate``.

Testing
-------

The tests live in ``tests/`` and are not part of the installed package.
pytest options are shared through ``pyproject.toml``. In-tree:

.. code-block:: sh

   meson test -C build --suite moist:python

Against an installed package (or by hand against the build directory with ``PYTHONPATH=<build>/python``), from the repository root:

.. code-block:: sh

   MOIST_REQUIRE_PYSCF=1 pytest python/tests

Run ``pytest``, not ``python -m pytest``, from inside ``python/``: the latter puts the source directory on ``sys.path`` and shadows the built package.
