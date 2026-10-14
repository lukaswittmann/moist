# Vendored FINUFFT sources

Origin: **FINUFFT** (Flatiron Institute Nonuniform Fast Fourier Transform) by the Flatiron
Institute / Simons Foundation: https://github.com/flatironinstitute/finufft
Commit: `0e9c10409580c482656b3f840c9d875f98c3a066` (2.6.0-dev, the FINUFFT revision this pin
matches in `subprojects/finufft.wrap`'s former `revision` field before that wrap was removed).
Every file below is a byte-for-byte copy of the upstream path at that commit, except the two
patched headers listed under "Local modifications". Relative paths under `include/`, `src/` and `fortran/` are kept identical to upstream
so the sources' own `#include` directives resolve unchanged.

## License

FINUFFT is Apache-2.0 (`LICENSE`, `NOTICE`). `NOTICE` documents FINUFFT's own third-party
attributions (XSIMD, ducc0/FFTW, the CMCL Fortran test drivers, the PSWF code derived from
`flatironinstitute/dmk`); those obligations are already satisfied by this PROVENANCE file plus
the separate `xsimd/` and `ducc0/` vendor directories, so `fortran/cmcl_license.txt` (upstream's
license for the CMCL Fortran example/test drivers under `fortran/examples` and `fortran/test`)
is not vendored here: none of the files below originate from those directories.

## Build

This subset is compiled unconditionally as CPU-only, FINUFFT_USE_DUCC0 code (no FFTW, no CUDA).
FINUFFT's own CMake/CPM machinery (which would fetch its own ducc0 and xsimd copies over the
network) is not used. `meson.build` in this directory is moist's own file (the only file here
that is not from upstream); it and `CMakeLists.txt` compile these sources directly against moist's own vendored `src/moist/math/fft/ducc0/` (include
dir `src/moist/math/fft/`, so `#include "ducc0/..."` resolves) and `src/moist/math/fft/xsimd/`.
The three ducc0 infra translation units (`threading.cc`, `mav.cc`, `string_utils.cc`) are only
ever compiled once, by `moist_fft`; FINUFFT links against that library rather than recompiling
them (see the ducc0 `PROVENANCE.md` in the sibling directory).

Mirroring upstream's `src/CMakeLists.txt` / `src/common/CMakeLists.txt`, the precision-dependent
sources are compiled twice (once with `FINUFFT_SINGLE` defined, once without) and the
precision-independent sources once.

## Files

| Path | Compiled | Precision |
|---|---|---|
| `include/finufft.h` | header | - |
| `include/finufft_errors.h` | header | - |
| `include/finufft_opts.h` | header | - |
| `include/finufft.fh` | header (Fortran `include`) | - |
| `include/finufft_mod.f90` | yes, compiled in place as Fortran module `finufft_mod` | - |
| `include/finufft/execute.hpp` | header | - |
| `include/finufft/finufft_eitherprec.h` | header | - |
| `include/finufft/heuristics.hpp` | header | - |
| `include/finufft/interp.hpp` | header | - |
| `include/finufft/makeplan.hpp` | header | - |
| `include/finufft/plan.hpp` | header | - |
| `include/finufft/setpts.hpp` | header | - |
| `include/finufft/simd.hpp` | header | - |
| `include/finufft/spread.hpp` | header | - |
| `include/finufft/spreadinterp.hpp` | header | - |
| `include/finufft/utils.hpp` | header | - |
| `include/finufft_common/common.h` | header | - |
| `include/finufft_common/constants.h` | header | - |
| `include/finufft_common/defines.h` | header | - |
| `include/finufft_common/kernel.h` | header | - |
| `include/finufft_common/pswf.h` | header | - |
| `include/finufft_common/safe_call.h` | header | - |
| `include/finufft_common/spread_opts.h` | header | - |
| `include/finufft_common/utils.h` | header | - |
| `src/makeplan.cpp` | yes | both (twice) |
| `src/setpts.cpp` | yes | both (twice) |
| `src/execute.cpp` | yes | both (twice) |
| `src/spreadinterp.cpp` | yes | both (twice) |
| `src/spreadinterp_1d.cpp` | yes | both (twice) |
| `src/spreadinterp_2d.cpp` | yes | both (twice) |
| `src/spreadinterp_3d.cpp` | yes | both (twice) |
| `src/fft.cpp` | yes | once (precision-templated internally) |
| `src/c_interface.cpp` | yes | once |
| `src/utils.cpp` | yes | once |
| `src/common/kernel.cpp` | yes | once |
| `src/common/pswf.cpp` | yes | once |
| `src/common/utils.cpp` | yes | once |
| `fortran/finufftfort.cpp` | yes | once |
| `LICENSE` | license text | - |
| `NOTICE` | license text | - |

The set is the transitive `#include` closure of the sources upstream compiles into the CPU,
DUCC0-backed, Fortran-enabled `finufft` target (`src/CMakeLists.txt`'s `FINUFFT_PRECISION_SOURCES`
+ `FINUFFT_COMMON_SOURCES` + `fortran/finufftfort.cpp`, plus `src/common/CMakeLists.txt`'s
`FINUFFT_COMMON_SOURCES`), verified header-by-header against every `#include "finufft...` /
`#include <finufft...>` line reachable from those translation units. `include/finufft/test_defs.hpp`
is upstream test-only scaffolding; nothing in this subset includes it, so it is not vendored.
GPU sources (`src/cuda/`, `include/cufinufft*`) are never reached from the CPU build and are not
vendored either.

## Local modifications

`include/finufft/spread.hpp` (`bin_sort_multithread_impl`) and `include/finufft/makeplan.hpp`
(`onedim_fseries_kernel`) split their work into `nt` chunks and let OpenMP thread `t` own chunk
`t`. That assumes `num_threads(nt)` always yields `nt` threads. A nested region (e.g. inside
test-drive's parallel suite runner or a host's own parallel region) runs fewer: the bin sort then
reads unsized per-thread histograms and crashes, and the kernel series is left partly
uncomputed. Both loops now work-share the chunks with `omp for schedule(static, 1)`. The results
are unchanged, since each chunk's arithmetic does not depend on the thread. The changes carry a
"moist patch" comment and a header notice (Apache-2.0 section 4(b)); `diff` against upstream
`0e9c104` shows them. Unfixed upstream as of 2026-09-29; drop the patch once upstream fixes it.
