MOIST: Modular and Open-source Implicit Solvation Toolkit
==============

[![License](https://img.shields.io/github/license/lukaswittmann/moist)](https://github.com/lukaswittmann/moist/blob/main/COPYING.LESSER) 
[![Version](https://img.shields.io/github/v/release/lukaswittmann/moist?include_prereleases)](https://github.com/lukaswittmann/moist/releases/latest)
[![Build](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml/badge.svg?branch=main)](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml)
[![Documentation](https://readthedocs.org/projects/moist/badge/?version=latest)](https://moist.readthedocs.io/en/latest/?badge=latest)
[![Coverage](https://codecov.io/gh/lukaswittmann/moist/branch/main/graph/badge.svg)](https://codecov.io/gh/lukaswittmann/moist)

> [!Note]
>  MOIST is currently in a pre-release state.
> This version only includes the cavity construction capabilities of the Modular and Open-source Implicit Solvation Toolkit.

### Building from Source

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- one of the supported build systems:
  - [meson](https://mesonbuild.com) version 0.57 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
  - [cmake](https://cmake.org) version 3.18 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
  - [fpm](https://github.com/fortran-lang/fpm) version 0.11.0 or newer
- a LAPACK / BLAS provider, like MKL or OpenBLAS

Currently this project supports GCC and Intel compilers.

#### Building with meson

Optional dependencies are
- [FFTW3](https://www.fftw.org) (version 3.3 or newer); required by the RISM models and any FFT-based routines; discovered via `pkg-config` (`fftw3`)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5) with Fortran bindings; for HDF5-based I/O of RISM results and caches
- FORD to build the developer documentation
- C compiler to test the C-API and compile the Python extension module
- Python 3.6 or newer with the CFFI package installed to build the Python API

##### Optional features

Several numerical features are disabled by default and are switched on through
meson options passed to `meson setup`:

| Option           | Default | Effect                                                           | Extra dependency                  |
| ---------------- | ------- | ---------------------------------------------------------------- | --------------------------------- |
| `-Drism=true`    | `false` | Build the RISM solvation model (defines `WITH_RISM`)              | FFTW3, HDF5 (added automatically) |
| `-Dfftw=true`    | `false` | Link FFTW3 for FFT-based routines (defines `WITH_FFTW`)           | FFTW3                             |
| `-Dhdf5=true`    | `false` | HDF5 I/O for RISM grids and large datasets (defines `WITH_HDF5`)  | HDF5 with Fortran bindings        |
| `-Dilp64=true`   | `false` | Use 64-bit-integer (ILP64) BLAS/LAPACK                           | ILP64 BLAS/LAPACK                 |
| `-Dopenmp=true`  | `true`  | OpenMP parallelisation (enables threaded FFT, `fftw3_threads`)   | OpenMP runtime                    |

Because `-Drism=true` already requires FFTW3, `-Dfftw=true` is only needed when you want FFTW without the RISM module.
For example, a RISM build with HDF5 output:

```sh
meson setup build -Drism=true -Dhdf5=true
```

These optional numerical features are available through the meson and cmake builds (the cmake build exposes the same toggles as `-DMOIST_*` options, see below); the fpm build always compiles the core configuration (no RISM/FFTW/HDF5).

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

This might require administrator access depending on the chosen install prefix.


#### Building with fpm

This project support the Fortran package manager (fpm).
Invoke fpm in the project root with

```
fpm build
```

To run the testsuite use

```
fpm test
```

You can access the ``moist`` program using the run subcommand

```
fpm run -- --help
```

To use ``moist`` for testing include it as dependency in your package manifest

```toml
[dependencies]
moist.git = "https://github.com/lukaswittmann/moist"
```

Note that the fpm build does not support exporting the C-API, it only provides access to the standalone binary.


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
| `-DMOIST_RISM=ON`    | `-Drism=true`    | `OFF`   |
| `-DMOIST_FFTW=ON`    | `-Dfftw=true`    | `OFF`   |
| `-DMOIST_HDF5=ON`    | `-Dhdf5=true`    | `OFF`   |
| `-DMOIST_ILP64=ON`   | `-Dilp64=true`   | `OFF`   |
| `-DMOIST_OPENMP=OFF` | `-Dopenmp=false` | `ON`    |
| `-DMOIST_API=OFF`    | `-Dapi=false`    | `ON`    |
| `-DMOIST_LAPACK=...` | `-Dlapack=...`   | `auto`  |

`-DMOIST_LAPACK` accepts `auto`, `mkl`, `openblas`, or `netlib`. For example, a
RISM build with HDF5 output:

```sh
cmake -B build -G Ninja -DMOIST_RISM=ON -DMOIST_HDF5=ON
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

## Development tooling

When working on the Fortran sources, please lint and format any changes before opening a pull request:

- [`fortitude`](https://fortitude.readthedocs.io/en/stable/) helps catch common mistakes. Run it from the project root with `fortitude check --output-format=concise`.
- [`fprettify`](https://github.com/fortran-lang/fprettify) keeps the formatting consistent. We currently use `fprettify -i 3` (3-space indentation) for code cleanup.

These checks will eventually be enforced by GitHub Actions, so running them locally now keeps the repository ready for automated CI.

### Developer documentation (FORD)

- The FORD project file is `ford.md` in the repo root.
- Build locally with `pip install ford` (or your package manager) and then `ford ford.md`.
- Generated HTML lands in `docs/ford/index.html` (ignored by git); open that file in a browser to browse the API.
- Adjust the `src_dir` or add extra Markdown pages inside `ford.md` if you extend the code layout.


## Usage

The `moist` command line tool is organised into subcommands. The general form is

```sh
moist <subcommand> [options] <input>
```

Run `moist --help` to list all subcommands, or `moist <subcommand> --help`
(for example `moist cavity drop svdw --help`) for the options of a specific
command. `moist --version`, `moist --citation`, and `moist --license` print the
version, the relevant literature references, and the full license text.

### Constructing a cavity

To build a molecular cavity with the SvdW-DROP scheme from an XYZ coordinate
file:

```sh
moist cavity drop svdw coord.xyz
```

Here `cavity` selects the cavity-only workflow, `drop` the DROP cavity construction, and `svdw` the solvent-van-der-Waals level set (the alternative level set is `cfc`).
Common options for this command are `--radii {cpcm,smd,d3,cosmo,bondi}` to choose the atomic radii set, `--nleb <N>` to set the Lebedev points per atom.
The other cavity constructors are available as `moist cavity {numsa,iswig} <coord>`.

### Other subcommands

<!-- - `moist model <gems|alpb|rism1d|rism3d> <input>` runs a full solvation model
  (the `rism1d`/`rism3d` models require the optional FFTW/RISM build described
  above). -->
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
Alternatively, the project is accessible by the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)).

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

Proceed with the build as described before and install the projects to make the Python API available in the selected prefix.

For the out-of-tree build see the instructions in the [`python`](https://github.com/lukaswittmann/moist/blob/main/python) directory.


## Citation

There is no dedicated toolkit paper yet, so for now please cite `moist` using the DROP cavity work, together with the references for the specific cavities, models, and solvers you used:

- L. Wittmann, A. Pausch, *A Smooth and Fully Differentiable Molecular Cavity
  Based on Discretization via Reference-Onto-Surface Projection*, ChemRxiv 2026.
  <https://doi.org/10.26434/chemrxiv.15003893/v2>

The full, context-appropriate citation list is maintained inside the program and printed by

```sh
moist --citation
```

## Optional: Build FFTW With PIC For RISM

If you enable RISM and shared library builds, FFTW should be available as a PIC/shared build.
One working setup is:

```sh
./configure CFLAGS="-O2 -fPIC" --enable-shared --enable-static --prefix=$HOME/.local
make -j
make install
```

If installed into `$HOME/.local`, ensure your environment can find it, for example:

```sh
export PKG_CONFIG_PATH="$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$HOME/.local/lib:$LD_LIBRARY_PATH"
```

## License

`moist` is free software released under the **LGPL-3.0-or-later** license (see [`COPYING`](COPYING) and [`COPYING.LESSER`](COPYING.LESSER)).

It also bundles third-party code under compatible permissive (BSD / MIT) and LGPL-3.0 licenses.
The full inventory, with origins and license texts, is in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

## Acknowledgements

`moist` builds on numerical software from several authors, gratefully acknowledged:

- **Jaś Kachnowicz** for help with the solvent properties and geometries.
- **Jacob Williams** ([jacobwilliams](https://github.com/jacobwilliams)) for modern Fortran implementations of various solvers and packages (SLSQP, L-BFGS-B, NLESolver, fmin, LSQR, LSMR, LUSOL, NumDiff).

These components retain their original copyright notices and licenses; the full texts accompany each component in the source tree and are catalogued in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
