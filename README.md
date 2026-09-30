<div align="center">
  
  # MOIST - *The Modular and Open-source Implicit Solvation Toolkit*

  [![License](https://img.shields.io/github/license/lukaswittmann/moist)](https://github.com/lukaswittmann/moist/blob/main/LICENSE)
  [![Version](https://img.shields.io/github/v/release/lukaswittmann/moist?include_prereleases)](https://github.com/lukaswittmann/moist/releases/latest)
  [![Build](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml/badge.svg?branch=main)](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml)
  [![Documentation](https://readthedocs.org/projects/moist/badge/?version=latest)](https://moist.readthedocs.io/en/latest/?badge=latest)
  [![Coverage](https://codecov.io/gh/lukaswittmann/moist/branch/main/graph/badge.svg)](https://codecov.io/gh/lukaswittmann/moist)

  [![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20Windows-informational)](#tested-platforms)
  [![Compilers](https://img.shields.io/badge/compilers-GCC%2014--16%20%7C%20Intel%20ifx%202025-informational)](#tested-platforms)
  [![Build systems](https://img.shields.io/badge/build%20systems-meson%20%7C%20CMake-informational)](#tested-platforms)
  [![Python](https://img.shields.io/badge/python-3.10%20%7C%203.13-informational)](#tested-platforms)
</div>

MOIST is a modular library for implicit solvation in molecular quantum chemistry, usable from Fortran, C, and Python.
It provides modern, smooth and differentiable cavity construction schemes (vdW iSwiG, SvdW-DROP, and isodensity-DROP), efficient polarizable continuum models (PCMs), pressure models (GOSTSHYP), and statistical solvation models (Reference Interaction Site Model, RISM), all implemented with performance and differentiability in mind.

> [!Note]
>  MOIST is currently in a pre-release state and is actively developed.  Contributions and discussions are very welcome!

### Tested platforms

Every combination below builds and passes the full test suite in [CI](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml).

| Platform                      | GCC 15         | Intel ifx 2025.1 | conda-forge GCC 15 | Python     |
| ----------------------------- | -------------- | ---------------- | ------------------ | ---------- |
| Linux x86_64 (Ubuntu 24.04)   | meson, CMake   | meson, CMake     | meson, CMake       | 3.10, 3.13 |
| Linux arm64 (Ubuntu 24.04)    | meson, CMake   | –                | meson, CMake       | 3.13       |
| macOS arm64 (macOS 15)        | meson, CMake   | –                | meson, CMake¹      | 3.13       |
| Windows x86_64 (MSYS2 UCRT64) | meson, CMake²  | –                | –                  | –          |

¹ with clang 19 and libc++ as the C/C++ compilers\
² MSYS2 ships only its current GCC, at present 16.x

- **BLAS/LAPACK:** reference LAPACK (Linux, GCC), MKL (ifx), OpenBLAS (conda-forge, Windows), Accelerate (macOS)
- **ILP64:** 64-bit integer BLAS/LAPACK with GCC 15 and OpenBLAS on Linux x86_64
- **Debug:** GCC 14 debug build with runtime checks (`-fcheck=bounds,do,mem,pointer`, `_GLIBCXX_ASSERTIONS`) on Linux x86_64
- **Python:** the bindings are tested against PySCF, both from a meson build and as a pip-installed package

The oldest tested compilers are GCC 14 and Intel ifx 2025.1; older versions may work but are not tested.

### Building from Source

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2018
- a C11 compiler and a C++17 compiler (*e.g.* gcc/g++ with gfortran, icx/icpx with ifx, or clang/clang++ with gfortran)
- one of the supported build systems:
  - [meson](https://mesonbuild.com) version 0.57 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.8.2 or newer
  - [cmake](https://cmake.org) version 3.25 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- a LAPACK / BLAS provider, like MKL, OpenBLAS, or Accelerate on macOS

#### Building with meson

Optional dependencies are
- FORD to build the developer documentation
- asciidoctor to build the man page
- Python 3.10 or newer with the CFFI and setuptools packages to build the Python API

##### Build options

The build is configured through meson options passed to `meson setup`:

| Option                     | Default   | Effect                                                                                      |
| -------------------------- | --------- | ------------------------------------------------------------------------------------------- |
| `-Dlapack=...`             | `auto`    | BLAS/LAPACK backend: `auto`, `mkl`, `mkl-rt`, `openblas`, `netlib`, `accelerate`, `custom`  |
| `-Dcustom_libraries=...`   | `[]`      | Libraries for the `custom` backend: paths, names, or a leading `-L<dir>`                    |
| `-Dilp64=true`             | `false`   | Use 64-bit-integer (ILP64) BLAS/LAPACK; needs an ILP64 BLAS/LAPACK                          |
| `-Dopenmp=false`           | `true`    | OpenMP parallelisation                                                                      |
| `-Dapi=false`              | `true`    | C API (`include/moist.h`)                                                                   |
| `-Dpython=true`            | `false`   | Python extension module, see [Python API](#python-api)                                      |
| `-Ddocs=true`              | `false`   | Man page (needs asciidoctor)                                                                |
| `-Ddev_tests=false`        | `true`    | The `dev_tester` binary with long-running diagnostic tests                                  |

The cmake build exposes the same toggles as `-DMOIST_*` options, see below.

Setup a default build with

```sh
meson setup build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile and run the projects testsuite use

```sh
meson test -C build --print-errorlogs
```

If the testsuite passes you can install with

```sh
meson configure build --prefix=/path/to/install
meson install -C build
```

#### Building with cmake

Configure and compile a default build with

```sh
cmake -B build -G Ninja
cmake --build build
```

The Fortran compiler is selected through the `FC` environment variable.
The optional numerical features map onto `-DMOIST_<FEATURE>` cache variables, with defaults matching the meson options:

| cmake option         | meson equivalent | Default |
| -------------------- | ---------------- | ------- |
| `-DMOIST_ILP64=ON`   | `-Dilp64=true`   | `OFF`   |
| `-DMOIST_OPENMP=OFF` | `-Dopenmp=false` | `ON`    |
| `-DMOIST_API=OFF`    | `-Dapi=false`    | `ON`    |
| `-DMOIST_LAPACK=...` | `-Dlapack=...`   | `auto`  |
| `-DMOIST_TESTS=OFF`  | —                | `ON`    |

`-DMOIST_LAPACK` accepts the same backends as meson (`auto`, `mkl`, `mkl-rt`, `openblas`, `netlib`, `accelerate`, `custom`);
`custom` takes its libraries from `-DMOIST_CUSTOM_LIBRARIES=...`. For example, a
build against ILP64 BLAS/LAPACK:

```sh
cmake -B build -G Ninja -DMOIST_ILP64=ON
cmake --build build
```

To install, set the prefix at configure time and install the built tree

```sh
cmake -B build -G Ninja -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build
cmake --install build
```

The unit-test suite is built by default (the `moist-tester` target; disable with `-DMOIST_TESTS=OFF`). Run it with

```sh
ctest --test-dir build --output-on-failure
```

## Usage

The `moist` command line tool is organised into subcommands. The general form is

```sh
moist <subcommand> [options] <input>
```

Run `moist --help` to list all subcommands, or `moist <subcommand> --help`
(for example `moist cavity drop svdw --help`) for the options of a specific
command. `moist --version`, `moist --citation`, and `moist --license` print the
version, the relevant literature references, and the license notice.

### Constructing a cavity

To build a molecular cavity with the SvdW-DROP scheme from an XYZ coordinate
file:

```sh
moist cavity drop svdw coord.xyz
```

Here `cavity` selects the cavity-only workflow, `drop` the DROP cavity construction, and `svdw` the solvent-van-der-Waals level set (the alternative level set is `cfc`).
Common options for this command are `--radii {cpcm,smd,d3,cosmo,bondi}` to choose the atomic radii set, `--nleb <N>` to set the Lebedev points per atom.
The other cavity constructors are available as `moist cavity {numsa,iswig} <coord>`, and `moist cavity mc {svdw,cfc} <coord>` integrates the level-set isosurface with marching cubes.

### Other subcommands

- `moist model gems <coord>` runs the GEMS solvation model; `--solvent` selects
  the solvent and `--charge` the molecular charge.
- `moist solvent <name>` reports the tabulated properties of a solvent by name
  or alias.

## API access

`moist` provides first class API support for Fortran, C and Python.
Other programming languages should try to interface with `moist` via one of those three APIs.
To provide first class API support for a new language the interface specification should be available as meson build files.

> [!Warning]
> The public APIs (Fortran, C, and Python) are in development and may change; users should not rely on strict backwards compatibility.
>
> All wishes or suggestions for the APIs (Fortran, C, Python) are very welcome; please [create a new issue](https://github.com/lukaswittmann/moist/issues/new).


### Fortran API

The recommended way to access the Fortran module API is by using `moist` as a meson subproject.

The complete API is available from `moist` module, the individual modules are available to the user as well but are not part of the public API and therefore not guaranteed to remain stable.
ABI compatibility is only guaranteed for the same minor version.

The communication with the Fortran API uses the `error_type` and `structure_type` of the modular computation tool chain library (mctc-lib) to handle errors and represent geometries, respectively.

### Python API

The Python API is disabled by default and can be built in-tree or out-of-tree.
The in-tree build is mainly meant for end users and packages.
To build the Python API with the normal project set the `python` option in the configuration step with

```sh
meson setup build -Dpython=true -Dpython_version=$(which python3)
```

The Python version can be used to select a different Python version, it defaults to `'python3'`.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Proceed with the build as described before. The package can be used directly from the build tree with `PYTHONPATH=build/python`, or installed to make it available in the selected prefix.
With `-Dpython=true` the test suite also runs the Python tests and requires numpy, pytest, and PySCF.

For the out-of-tree build see the instructions in the [`python`](https://github.com/lukaswittmann/moist/blob/main/python) directory.


## Contributing

Contributions are welcome; see [`CONTRIBUTING.md`](CONTRIBUTING.md) for the code style, linting, testing, and documentation guidelines.

## Citation

There is no dedicated toolkit paper yet, so for now please cite `moist` using the DROP cavity work, together with the references for the specific cavities, models, and solvers you used:

- L. Wittmann, A. Pausch, *A Smooth and Fully Differentiable Molecular Cavity
  Based on Discretization via Reference-Onto-Surface Projection*, ChemRxiv 2026.
  <https://doi.org/10.26434/chemrxiv.15003893/v2>

For the isodensity cavity (ρ-DROP), please also cite

- L. Wittmann, *From Density to Boundary and Back: A Fully Differentiable,
  Self-Consistent Isodensity Cavity*, ChemRxiv 2026.
  <https://doi.org/10.26434/chemrxiv.15007095/v1>

The full, context-appropriate citation list is maintained inside the program and printed by

```sh
moist --citation
```

## License

`moist` is free software released under the **MPL-2.0** license (see [`LICENSE`](LICENSE)).

It also bundles third-party code under separate licenses. The full inventory, with
origins and license texts, is in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

## Acknowledgements

`moist` builds on numerical software from several authors, gratefully acknowledged:

- **Jaś Kachnowicz** for help with the solvent properties and geometries.
- **Jacob Williams** ([jacobwilliams](https://github.com/jacobwilliams)) for modern Fortran implementations of various solvers and packages (SLSQP, L-BFGS-B, NLESolver, fmin, LSQR, LSMR, LUSOL).

These components retain their original copyright notices and licenses; the full texts accompany each component in the source tree and are catalogued in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
