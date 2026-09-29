# Vendored xsimd headers

Origin: **xsimd** (C++ wrappers for SIMD intrinsics) by the QuantStack / xtensor-stack project:
https://github.com/xtensor-stack/xsimd
Commit: `6842624fc8adafd7168a999e7150b384411da448` — this is the exact revision FINUFFT itself
pins as `XSIMD_VERSION` in its `CMakeLists.txt` ("this commit fixes gcc-10 for now"), so moist's
vendored copy matches what FINUFFT would otherwise fetch on its own.
Every file below is a byte-for-byte copy of the upstream path at that commit; nothing was
modified.

## License

xsimd is BSD-3-Clause (`LICENSE`).

## Files

Header-only: everything under `include/xsimd/` (`arch/`, `config/`, `math/`, `memory/`, `types/`
and the top-level headers), plus `LICENSE`. Nothing else from the upstream tree (benchmarks,
tests, docs, its own CMake/meson build files) is vendored; xsimd is never compiled as its own
target, only used as an include directory by `src/moist/math/fft/finufft/`. `meson.build` in this
directory is moist's own file (declares `xsimd_dep` and installs `LICENSE`), not upstream's.
