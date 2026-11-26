# Vendored libmsym sources

Origin: **libmsym** (molecular point group symmetry library) by Marcus Johansson: https://github.com/mcodev31/libmsym
Commit: `85e47232376a8e735c2a7b5283f40b59b3953db1` (master, 2024-06-15, "select first power of Sn"; CMake project version 0.2.4).
Every `.c` and `.h` file below is a byte-for-byte copy of `src/<file>` at that commit, plus `LICENSE` from the repository root.

## Local additions

- `msym_EXPORTS.h`: written for moist. Upstream generates it with CMake's `GenerateExportHeader`; here `MSYM_EXPORT` is empty, since libmsym is compiled into libmoist and only reached through moist's shim (`../msym_shim.c`), which exports nothing of libmsym.

## Local patches

- `symmetry.c` (`reduceSymmetry`, `else if (cinf[0] != NULL)`): compare against the saved vector `v` instead of `cinf[0]->v`, since `realloc(rsops, ...)` in the loop body can move and free the array `cinf[0]` points into.

## License

MIT (Copyright (c) 2015 Marcus Johansson); full text in `LICENSE`, which must ship with every distribution. The summary row in `THIRD_PARTY_NOTICES.md` points here.

## Files

Include directory: `src/moist/math/symmetry/libmsym/` (the sources include each other by bare file name).
All 17 translation units are compiled (`moist_libmsym`): `msym.c` references the SALC and character-table code even though moist only uses detection, operations and subgroups.

| Path | Compiled |
|---|---|
| `basis_function.c` / `.h` | yes |
| `character_table.c` / `.h` | yes |
| `context.c` / `.h` | yes |
| `debug.c` / `.h` | yes |
| `elements.c` / `.h` | yes |
| `equivalence_set.c` / `.h` | yes |
| `geometry.c` / `.h` | yes |
| `linalg.c` / `.h` | yes |
| `msym.c` / `.h` | yes |
| `msym_error.c` / `.h` | yes |
| `permutation.c` / `.h` | yes |
| `point_group.c` / `.h` | yes |
| `rsh.c` / `.h` | yes |
| `subspace.c` / `.h` | yes |
| `symmetrize.c` / `.h` | yes |
| `symmetry.c` / `.h` | yes |
| `symop.c` / `.h` | yes |

## Compiler requirements

C99 with variable-length arrays (used throughout, including in prototypes). GCC, Clang and Intel icx accept them; MSVC's `cl` does not.
