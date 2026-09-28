.. _installation:

Installing MOIST
================

Requirements: a Fortran 2018 compiler (GCC or Intel), BLAS/LAPACK, and meson
with ninja, cmake with ninja, or fpm. A C compiler is needed for the C API and
the Python extension; Python 3.10 with CFFI and NumPy for the Python API.

Meson
-----

.. code-block:: sh

   meson setup build
   meson test -C build --print-errorlogs
   meson install -C build

Select the compiler with ``FC``; set the prefix with
``meson configure build --prefix=/path/to/install``.

.. list-table:: Meson options
   :header-rows: 1

   * - Option
     - Default
     - Effect
   * - ``-Dlapack=``
     - ``auto``
     - ``mkl``, ``mkl-rt``, ``openblas``, ``netlib``, ``accelerate`` or
       ``custom`` (libraries in ``-Dcustom_libraries=``)
   * - ``-Dilp64=``
     - ``false``
     - 64-bit-integer BLAS/LAPACK
   * - ``-Dopenmp=``
     - ``true``
     - OpenMP parallelisation
   * - ``-Dapi=``
     - ``true``
     - C API
   * - ``-Dpython=``
     - ``false``
     - Python extension module, linked against ``-Dpython_version=``
       (default ``python3``)

CMake
-----

.. code-block:: sh

   cmake -B build -G Ninja -DCMAKE_INSTALL_PREFIX=/path/to/install
   cmake --build build
   ctest --test-dir build --output-on-failure
   cmake --install build

The meson options map onto ``-DMOIST_LAPACK``, ``-DMOIST_ILP64``,
``-DMOIST_OPENMP``, ``-DMOIST_API`` and ``-DMOIST_TESTS`` with the same
defaults. CMake does not build the Python API.

fpm
---

.. code-block:: sh

   fpm build
   fpm test
   fpm run -- --help

As a dependency, add ``moist.git = "https://github.com/lukaswittmann/moist"``
under ``[dependencies]``. The fpm build provides the Fortran API and the
``moist`` program, not the C API.

Python
------

In-tree: configure with ``-Dpython=true -Dpython_version=$(which python3)``,
build, and install; or import from the build tree with
``PYTHONPATH=<build>/python``. A standalone wheel is built from the
``python`` directory with meson-python:

.. code-block:: sh

   pip install ./python -Csetup-args=-Dlapack=accelerate

PySCF support is the ``pyscf`` extra. See :doc:`/reference/python`.
