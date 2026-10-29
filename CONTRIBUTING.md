# Contributing to MOIST

Thank you for your interest in contributing to MOIST! We're excited to collaborate with you.

Before you start, please read the following guidelines; this will allow for faster and better collaboration.

There are many ways to contribute: writing tutorials or blog posts, improving the documentation, submitting [bug reports and feature requests](https://github.com/lukaswittmann/moist/issues/new/choose) or writing code which can be incorporated into MOIST itself.

## Philosophy

MOIST aims at reusable solvation models and infrastructure; contributions should follow the following principles:

1. **Consistent API**: the API for models and components has to be well-defined and consistent.
2. **Separation of responsibilities**: QM hosts provide electronic-structure quantities; MOIST provides solvation models and their response.
3. **Reusable building blocks**: reusable components composed through OOP and abstraction.
4. **Performance-awareness**: performance should be considered as part of the architecture.
5. **Method development**: the architecture should allow for method development.

## Workflow

- Please open an issue first for major features, infrastructure or API changes.
- Branch off `main` and open a pull request against `main`; keep one topic per pull request. Large PRs should be split into smaller ones.
- All Build workflow jobs have to pass before merging (CI).
- Contributions are licensed under the [MPL-2.0](LICENSE), the license of the project.

## Building and Testing

- Meson is the primary build system; CMake is tested in CI as well, so changes have to work with both (use separate build directories):
  - Meson: `meson setup build-meson && meson test -C build-meson --print-errorlogs`
  - CMake: `meson subprojects download && cmake -B build-cmake -G Ninja && cmake --build build-cmake && ctest --test-dir build-cmake`
- Python bindings and PySCF coupling: configure Meson with `-Dpython=true`, then run `meson test -C build-meson --suite python` (bindings and PySCF tests); when running pytest directly, `MOIST_REQUIRE_PYSCF=1` turns skipped PySCF tests into failures.
- Long-running diagnostic tests are built with `-Ddev_tests=true` (`test/dev/`).
- CI covers GCC, ifx, classic ifort, ILP64, macOS and Windows; stick to standard Fortran 2018 without compiler extensions.

## Guidelines

### Printing

- Adhere to the verbosity level (levels in `src/moist/context.f90`); warnings always print:
  - `0`: silent
  - `1`: quiet (main results, once per run)
  - `2`: normal, the library default for host integration: settings once per run and a summary per geometry, but no output during `get_energy`, `get_response` or `get_gradient`, so nothing prints inside the host's SCF iterations
  - `3` and `4`: more verbose (all settings and tables, solver summaries)
  - `5` and `6`: solver traces and debug output
- Route output through the run context (`moist_context`), whose `message` procedure is gated by verbosity (`debug_message` by the debug flag); avoid bare `print`/`write` statements.
- Use `prettyprinter` and `prettylistprinter` (see `moist_utils_prettyprint` and `moist_utils_prettylistprint` in `src/moist/utils/`) for printing.

### Code Style

- Lint with [`fortitude`](https://fortitude.readthedocs.io/en/stable/) from the project root: `fortitude check --output-format=concise src/ app/ test/` (configured in `fortitude.toml`); the Lint workflow runs the same check on every pull request.
- Use 3-space indentation and match the style of the surrounding code; [`fprettify`](https://github.com/fortran-lang/fprettify) `-i 3` is a good start for new files, but do not run it over existing files, as it would also reformat unrelated lines.
- Use `ford`-compatible comment style (config in `ford.md`, docmark `!>`):
  - One-line summary above every module, type, and procedure; an empty `!>` line, then one `@param[in|inout|out]` line per dummy argument.
  - A `!>` comment above every dummy declaration and type component, except the `error` dummy of a test procedure whose only argument it is.
  - Dimensions should be indicated by e.g. `(nat, 3, 3)`.
  - Give units if appropriate.
- Use good-practice docstring formatting:
  - Keep it short and use bullet points.
  - Keep units, shapes, ranges, and defaults; no filler words ("This routine", "is used to").
  - ASCII-only in source and comments.
- Report failures through `error_type`/`fatal_error` from `mctc_env` rather than silently returning defaults.
- Keep code thread-safe under OpenMP (`-Dopenmp=true`); avoid `save`d state.
- Python code (`python/moist/`) supports Python 3.10 and newer.

### APIs

- Expose new features consistently in all three interfaces, each with tests:
  - Fortran
  - C: `src/moist/api.f90` and `include/moist.h`, tested in `test/api/example.c`
  - Python: `python/moist/library.py` and `python/moist/interface.py`, tested in `python/tests/`
- Keep the interfaces as consistent with each other as possible: the same names, argument order, units, defaults and error behavior, deviating only where a language requires it (e.g. `lambda_` in Python, `moist_` prefixes and handles in C).
- Keep the C ABI stable (see the rules in `include/moist.h`): options structs only grow at the end, and enum values should never change.
- New settings extend `moist_model_parameters_type` (JSON/TOML input/output, validation, defaults, printing); define defaults in one place and keep the documentation in sync.
- Note breaking changes (renamed arguments, changed defaults) in the release notes.

### Generated and Vendored Code

- Files marked "generated by pcmgen, do not edit" (the PCM kernels in `src/moist/model/continuum/component/pcm/kernel/` and the DROP LSF kernels) are changed in the generator and regenerated; the header records the generator commit.
- Vendored third-party code stays close to upstream; record it in `THIRD_PARTY_NOTICES.md`.

### Tests

- The code has to be well tested.
- One suite per file (module): `test/unit/test_<id>.f90` defines module `test_<id>` with `collect_<id>`. Register new suites in `test/unit/meson.build` (`tests` list), `test/unit/main.F90`, and `CMakeLists.txt` (`moist_test_suites`).
- Run with `meson test -C build-meson --print-errorlogs`. Suites run their tests in parallel, so avoid `save`d state in test helpers.
- Every new feature or bug fix comes with a test.
- Check analytic derivatives (Fock contributions, gradients, surface weights) against finite differences.
- Expected failures (rejected input, guarded misuse): one test per failing case, registered with `new_unittest(..., should_fail=.true.)`, which returns the library error as its test failure.

### Documentation

- User documentation is Sphinx in `docs/` (see `docs/requirements.txt`); build with `sphinx-build -E -W -b html docs docs/_build/html` (`-E` avoids stale pages, `-W` fails on warnings).
- API documentation is generated by FORD from the source comments: `pip install ford` and `ford ford.md`; the HTML lands in `docs/ford/index.html` (ignored by git). Adjust `src_dir` in `ford.md` if you add a new source tree.
- When adding a cavity, model, or component:
  - Add or update its `.rst` page and the support tables in `docs/cavities/index.rst` or `docs/models/components/index.rst`.
  - Add its references to `docs/_static/references.bib` and to `src/moist/output/citations.f90` (printed by `moist --citation`).
  - Add it to the feature list in `README.md`.
- Solvent IDs are frozen; new solvents are appended.
- Keep logos in `docs/_static/logos/` and figures in `docs/_static/figures/<topic>/`.

## Versioning and Releases

- MOIST follows [SemVer](https://semver.org) with an optional pre-release suffix (e.g. `1.2.1-alpha`).
- To bump the version, edit `VERSION`, run `config/sync-version.sh` (checked in CI), and update the version string in `test/api/example.c` by hand.
- Add an entry to `docs/changelog.rst`.
- Pushing a `v*` tag builds the release binaries and publishes the GitHub release.
