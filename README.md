<div align="center">
  
  # MOIST - *The Modular and Open-source Implicit Solvation Toolkit*

  [![License](https://img.shields.io/github/license/lukaswittmann/moist)](https://github.com/lukaswittmann/moist/blob/main/LICENSE)
  [![Version](https://img.shields.io/github/v/release/lukaswittmann/moist?include_prereleases)](https://github.com/lukaswittmann/moist/releases/latest)
  [![DOI](https://zenodo.org/badge/1255788295.svg)](https://zenodo.org/badge/latestdoi/1255788295)
  [![Build](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml/badge.svg?branch=main)](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml)
  [![Documentation](https://readthedocs.org/projects/moist/badge/?version=latest)](https://moist.readthedocs.io/en/latest/?badge=latest)
  [![Coverage](https://codecov.io/gh/lukaswittmann/moist/branch/main/graph/badge.svg)](https://codecov.io/gh/lukaswittmann/moist)

</div>

<br clear="right">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/_static/logos/moist-animated-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/_static/logos/moist-animated-light.svg">
  <img align="right" src="docs/_static/logos/moist-animated-light.svg" width="420">
</picture>

MOIST is a modular library for modern implicit solvation in molecular quantum chemistry, usable from Fortran, C, and Python, and implemented with performance and differentiability in mind.

> [!Note]
> MOIST is in pre-release: the Fortran, C, and Python APIs may still change between releases. Feedback, feature requests and contributions are very welcome; please [open an issue](https://github.com/lukaswittmann/moist/issues/new).

## Short Feature List

- **Cavities:** smooth, fully differentiable *Discretization via Reference-Onto-Surface Projection* ([DROP](https://moist.readthedocs.io/en/latest/cavities/drop.html)) cavities: *Smooth van der Waals* ([SvdW-DROP](https://moist.readthedocs.io/en/latest/cavities/svdw.html)) and *electron-isodensity* ([ρ-DROP](https://moist.readthedocs.io/en/latest/cavities/isodensity.html)) surfaces, the vdW *improved Switching Gaussian* ([iSwiG](https://moist.readthedocs.io/en/latest/cavities/iswig.html)) cavity, and [marching tetrahedra](https://moist.readthedocs.io/en/latest/cavities/mt.html); *Numerical Surface Area* ([NumSA](https://moist.readthedocs.io/en/latest/cavities/numsa.html)) and [marching cubes](https://moist.readthedocs.io/en/latest/cavities/mc.html) for reference areas and volumes
- **Polarizable Continuum Models:** *Conductor-like PCM* ([CPCM](https://moist.readthedocs.io/en/latest/models/components/pcm/cpcm.html)), *Conductor-like Screening Model* ([COSMO](https://moist.readthedocs.io/en/latest/models/components/pcm/cosmo.html)), *Integral Equation Formalism PCM* ([IEF-PCM](https://moist.readthedocs.io/en/latest/models/components/pcm/iefpcm.html)), *Surface and Simulation of Volume Polarization for Electrostatics* ([SS(V)PE](https://moist.readthedocs.io/en/latest/models/components/pcm/ssvpe.html)), *Dielectric PCM* ([DPCM](https://moist.readthedocs.io/en/latest/models/components/pcm/dpcm.html)), *Solvation-Layer Interface Condition* ([SLIC](https://moist.readthedocs.io/en/latest/models/components/pcm/slic.html)), *Linearized Poisson-Boltzmann* ([LPB](https://moist.readthedocs.io/en/latest/models/components/pcm/lpb.html)), and [nonlocal electrostatics](https://moist.readthedocs.io/en/latest/models/components/pcm/nonlocal.html)
- **Pressure Models:** classical *Pressure-Volume* ([PV](https://moist.readthedocs.io/en/latest/models/components/pv.html)) and *Gaussians On Surface Tesserae Simulating HYdrostatic Pressure* ([GOSTSHYP](https://moist.readthedocs.io/en/latest/models/components/gostshyp.html)) pressure components
- **Reference Interaction Site Models:** 1D and 3D [RISM](https://moist.readthedocs.io/en/latest/models/moz.html) (in development) with a modular [Lennard-Jones potential library](https://moist.readthedocs.io/en/latest/models/moz/potential.html) (UFF, DREIDING, GAFF2, OPLS-AA, SPC/E, and transition-metal parameters)
- **Interfaces:** [Fortran](https://moist.readthedocs.io/en/latest/reference/fortran.html), [C](https://moist.readthedocs.io/en/latest/reference/c.html) and [Python](https://moist.readthedocs.io/en/latest/reference/python.html) APIs, a [PySCF](https://moist.readthedocs.io/en/latest/reference/pyscf.html) module (HF and DFT energies and nuclear gradients), and the `moist` CLI; [ORCA and TURBOMOLE implementation will be available with their next major releases](https://moist.readthedocs.io/en/latest/integrations.html)


## Installation

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2018
- a C11 compiler and a C++17 compiler (*e.g.* gcc/g++ with gfortran, icx/icpx with ifx or ifort, or clang/clang++ with gfortran)
- one of the supported build systems:
  - [meson](https://mesonbuild.com) version 0.57 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.8.2 or newer
  - [cmake](https://cmake.org) version 3.25 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- a [LAPACK](https://www.netlib.org/lapack/) / BLAS provider, like [MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html), [OpenBLAS](https://www.openblas.net), or [Accelerate](https://developer.apple.com/documentation/accelerate) on macOS

See [Tested platforms](#tested-platforms) for the compilers and platforms covered by CI.

### Building with meson

Optional dependencies are
- [FORD](https://github.com/Fortran-FOSS-Programmers/ford) to build the developer documentation
- [asciidoctor](https://asciidoctor.org) to build the man page
- Python 3.10 or newer with the [CFFI](https://cffi.readthedocs.io) and [setuptools](https://setuptools.pypa.io) packages to build the [Python API](#python-api)

Setup a default build with

```sh
meson setup build
```

You can select the Fortran compiler by the `FC` environment variable. To compile and run the projects testsuite use

```sh
meson test -C build --print-errorlogs
```

If the testsuite passes you can install with

```sh
meson configure build --prefix=/path/to/install
meson install -C build
```

<details>
<summary>Build options</summary>

The build is configured through meson options passed to `meson setup`:

| Option                     | Default   | Effect                                                                                      |
| -------------------------- | --------- | ------------------------------------------------------------------------------------------- |
| `-Dlapack=...`             | `auto`    | BLAS/LAPACK backend: `auto`, `mkl`, `mkl-rt`, `openblas`, `netlib`, `accelerate`, `custom`  |
| `-Dcustom_libraries=...`   | `[]`      | Libraries for the `custom` backend: paths, names, or a leading `-L<dir>`                    |
| `-Dilp64=true`             | `false`   | Use 64-bit-integer (ILP64) BLAS/LAPACK; needs an ILP64 BLAS/LAPACK                          |
| `-Dopenmp=false`           | `true`    | OpenMP parallelisation                                                                      |
| `-Dapi=false`              | `true`    | [C API](https://moist.readthedocs.io/en/latest/reference/c.html) (`include/moist.h`)                                         |
| `-Dpython=true`            | `false`   | Python extension module, see [Python API](#python-api)                                      |
| `-Ddocs=true`              | `false`   | Man page (needs asciidoctor)                                                                |
| `-Ddev_tests=false`        | `true`    | The `dev_tester` binary with long-running diagnostic tests                                  |

The cmake build exposes the same toggles as `-DMOIST_*` options, see below.

</details>

### Building with cmake

Configure and compile a default build with

```sh
cmake -B build -G Ninja
cmake --build build
```

The Fortran compiler is selected through the `FC` environment variable. To install, set the prefix at configure time and install the built tree

```sh
cmake -B build -G Ninja -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build
cmake --install build
```

The unit-test suite is built by default (the `moist-tester` target; disable with `-DMOIST_TESTS=OFF`). Run it with

```sh
ctest --test-dir build --output-on-failure
```

<details>
<summary>Build options</summary>

The optional numerical features map onto `-DMOIST_<FEATURE>` cache variables, with defaults matching the meson options:

| cmake option         | meson equivalent | Default |
| -------------------- | ---------------- | ------- |
| `-DMOIST_ILP64=ON`   | `-Dilp64=true`   | `OFF`   |
| `-DMOIST_OPENMP=OFF` | `-Dopenmp=false` | `ON`    |
| `-DMOIST_API=OFF`    | `-Dapi=false`    | `ON`    |
| `-DMOIST_LAPACK=...` | `-Dlapack=...`   | `auto`  |
| `-DMOIST_TESTS=OFF`  | —                | `ON`    |

`-DMOIST_LAPACK` accepts the same backends as meson (`auto`, `mkl`, `mkl-rt`, `openblas`, `netlib`, `accelerate`, `custom`); `custom` takes its libraries from `-DMOIST_CUSTOM_LIBRARIES=...`. For example, a build against ILP64 BLAS/LAPACK:

```sh
cmake -B build -G Ninja -DMOIST_ILP64=ON
cmake --build build
```

</details>

## Usage

The `moist` command line tool is organised into subcommands. The general form is

```sh
moist <subcommand> [options] <input>
```

Run `moist --help` to list all subcommands, or `moist <subcommand> --help` (for example `moist cavity drop svdw --help`) for the options of a specific command. `moist --version`, `moist --citation`, and `moist --license` print the version, the relevant literature references, and the license notice.

### Constructing a cavity

To build a molecular cavity with the [SvdW-DROP](https://moist.readthedocs.io/en/latest/cavities/svdw.html) scheme from an XYZ coordinate file:

```sh
moist cavity drop svdw coord.xyz
```

Here `cavity` selects the cavity-only workflow, `drop` the [DROP](https://moist.readthedocs.io/en/latest/cavities/drop.html) cavity construction, and `svdw` the smooth van der Waals level set (the alternative level set is [`cfc`](https://moist.readthedocs.io/en/latest/cavities/cfc.html)). Common options for this command are `--radii {cpcm,smd,d3,cosmo,bondi}` to choose the atomic radii set, `--nleb <N>` to set the Lebedev points per atom. The other cavity constructors, [NumSA](https://moist.readthedocs.io/en/latest/cavities/numsa.html) and [iSwiG](https://moist.readthedocs.io/en/latest/cavities/iswig.html), are available as `moist cavity {numsa,iswig} <coord>`, and `moist cavity mc {svdw,cfc} <coord>` integrates the level-set isosurface with [marching cubes](https://moist.readthedocs.io/en/latest/cavities/mc.html).

### Other subcommands

- `moist solvent <name>` reports the tabulated properties of a solvent by name
  or alias.

## API access

`moist` provides first class API support for [Fortran](https://moist.readthedocs.io/en/latest/reference/fortran.html), [C](https://moist.readthedocs.io/en/latest/reference/c.html) and [Python](https://moist.readthedocs.io/en/latest/reference/python.html); see the [API reference](https://moist.readthedocs.io/en/latest/reference/index.html) for an overview. Other programming languages should try to interface with `moist` via one of those three APIs. To provide first class API support for a new language the interface specification should be available as meson build files.

### Fortran API

The recommended way to access the [Fortran module API](https://moist.readthedocs.io/en/latest/reference/fortran.html) is by using `moist` as a [meson subproject](https://mesonbuild.com/Subprojects.html).

The complete API is available from `moist` module, the individual modules are available to the user as well but are not part of the public API and therefore not guaranteed to remain stable. ABI compatibility is only guaranteed for the same minor version.

The communication with the Fortran API uses the `error_type` and `structure_type` of the modular computation tool chain library ([mctc-lib](https://github.com/grimme-lab/mctc-lib)) to handle errors and represent geometries, respectively.

### Python API

The [Python API](https://moist.readthedocs.io/en/latest/reference/python.html) is disabled by default and can be built in-tree or out-of-tree. The in-tree build is mainly meant for end users and packages. To build the Python API with the normal project set the `python` option in the configuration step with

```sh
meson setup build -Dpython=true -Dpython_version=$(which python3)
```

The Python version can be used to select a different Python version, it defaults to `'python3'`. Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Proceed with the build as described before. The package can be used directly from the build tree with `PYTHONPATH=build/python`, or installed to make it available in the selected prefix. With `-Dpython=true` the test suite also runs the Python tests and requires [numpy](https://numpy.org), [pytest](https://pytest.org), and [PySCF](https://pyscf.org).

For the out-of-tree build see the instructions in the [`python`](https://github.com/lukaswittmann/moist/blob/main/python) directory.

## Tested platforms

Every combination below builds and passes the full test suite in [CI](https://github.com/lukaswittmann/moist/actions/workflows/ci_build.yml). The oldest tested compilers are GCC 14, Intel ifx 2025.1 and classic Intel ifort 2021.12; older versions may work but are not tested.

<details>
<summary>Platform and compiler matrix</summary>

| Platform                      | GCC 15         | [conda-forge](https://conda-forge.org) GCC 15 | Intel ifx 2025.1 | Intel ifort 2021.12 | Python     |
| ----------------------------- | -------------- | ------------------ | ---------------- | ------------------- | ---------- |
| Linux x86_64 (Ubuntu 24.04)   | meson, CMake   | meson, CMake       | meson, CMake     | meson, CMake        | 3.10, 3.13 |
| Linux arm64 (Ubuntu 24.04)    | meson, CMake   | meson, CMake       | –                | –                   | 3.13       |
| macOS arm64 (macOS 15)        | meson, CMake   | meson, CMake<sup>1</sup>      | –                | –                   | 3.13       |
| Windows x86_64 ([MSYS2](https://www.msys2.org) UCRT64) | meson, CMake<sup>2</sup>  | –                  | –                | –                   | –          |

<sup>1)</sup> With clang 19 and libc++ as the C/C++ compilers.
<sup>2)</sup> MSYS2 ships only its current GCC, at present 16.x.

- **BLAS/LAPACK:** reference LAPACK (Linux, GCC), MKL (ifx, ifort), OpenBLAS (conda-forge, Windows), Accelerate (macOS)
- **ILP64:** 64-bit integer BLAS/LAPACK with GCC 15 and OpenBLAS on Linux x86_64
- **Debug:** GCC 14 debug build with runtime checks (`-fcheck=bounds,do,mem,pointer`, `_GLIBCXX_ASSERTIONS`) on Linux x86_64
- **Python:** the bindings are tested against PySCF, both from a meson build and as a pip-installed package

</details>

## Contributing

Contributions are welcome; see [`CONTRIBUTING.md`](CONTRIBUTING.md) for the code style, linting, testing, and documentation guidelines.

## Citation

If you use `moist`, please cite the software itself (see [`CITATION.cff`](CITATION.cff); the DOI badge above links to the latest release):

- L. Wittmann, *MOIST: The Modular and Open-source Implicit Solvation Toolkit*
  [Computer software], Zenodo 2026.
  <https://doi.org/10.5281/zenodo.23473416>

together with the references for the cavities, models, and solvers you used. For the DROP cavities, these include

- L. Wittmann, A. Pausch, *A Smooth and Fully Differentiable Molecular Cavity
  Based on Discretization via Reference-Onto-Surface Projection*, ChemRxiv 2026.
  <https://doi.org/10.26434/chemrxiv.15003893/v2>

For the isodensity cavity (ρ-DROP), please also cite

- L. Wittmann, *From Density to Boundary and Back: A Fully Differentiable,
  Self-Consistent Isodensity Cavity*, ChemRxiv 2026.
  <https://doi.org/10.26434/chemrxiv.15007095/v1>

At the end of a run, `moist` prints the citations for the cavities, models, and solvers it used; the full list is printed by

```sh
moist --citation
```

## License

`moist` is free software released under the **MPL-2.0** license (see [`LICENSE`](LICENSE)).

It also bundles third-party code under separate licenses. The full inventory, with origins and license texts, is in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Acknowledgements

We gratefully acknowledge the people and projects `moist` builds on:

- Jaś Kachnowicz for help with the initial version of the solvent properties and geometries.
- [Jacob Williams](https://github.com/jacobwilliams) for modern Fortran implementations of the SLSQP, L-BFGS-B, nlesolver, fmin, LSQR, LSMR, and LUSOL solvers.
- The [FINUFFT](https://github.com/flatironinstitute/finufft) developers at the Flatiron Institute, Martin Reinecke for [ducc0](https://gitlab.mpcdf.mpg.de/mtr/ducc), and the [xsimd](https://github.com/xtensor-stack/xsimd) developers for the FFT and SIMD code.
- The authors of the force-field parameter sources: Peter Boyd for [lammps_interface](https://github.com/peteboyd/lammps_interface); the [openmmforcefields](https://github.com/openmm/openmmforcefields) developers; Leela S. Dodda and the Jorgensen Lab for [OPLSAA-DB](https://github.com/leelasd/OPLSAA-DB); Gajanan Kanagalingam, Sebastian Schmitt, Florian Fleckenstein, and Simon Stephan for [TUK-FFDat](https://doi.org/10.5281/zenodo.8116422); and the [VeloxChem](https://github.com/VeloxChem/VeloxChem) developers.
- The developers of the [mctc-lib](https://github.com/grimme-lab/mctc-lib), [mstore](https://github.com/grimme-lab/mstore), [toml-f](https://github.com/toml-f/toml-f), [jonquil](https://github.com/toml-f/jonquil), [test-drive](https://github.com/fortran-lang/test-drive), and [fclap](https://github.com/chselz/fclap) projects.

These components retain their original copyright notices and licenses; the full texts accompany each component in the source tree and are catalogued in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
