# Third-party code and licenses

`moist` is distributed under the **MPL-2.0** license (see `LICENSE`). It incorporates
third-party code under other licenses.
This file documents that code, its origin, and the applicable license, as required by
those licenses (in particular the BSD source- and binary-redistribution clauses).

Each third-party component remains under its own terms; its license is listed alongside
it below.

## Vendored code

| Component | Path in `moist` | Upstream repository | Original author(s) | License | License file |
|---|---|---|---|---|---|
| SLSQP (+ BVLS) | `src/moist/math/solver/slsqp/` | jacobwilliams/slsqp | Dieter Kraft (1988); ACM (1994); BVLS: Lawson & Hanson (netlib, public domain); modern Fortran by J. Williams | BSD-3-Clause (+ MIT-style original grant) | `src/moist/math/solver/slsqp/LICENSE` |
| L-BFGS-B | `src/moist/math/solver/lbfgsb/` | jacobwilliams/lbfgsb | L-BFGS-B 3.0: J. Nocedal & J. L. Morales; modern Fortran by J. Williams | BSD-3-Clause ("New BSD") | `src/moist/math/solver/lbfgsb/LICENSE` |
| nlesolver | `src/moist/math/solver/newton/` | jacobwilliams/nlesolver-fortran | J. Williams | BSD-3-Clause | `src/moist/math/solver/newton/LICENSE` |
| fmin | `src/moist/math/solver/fmin/` | jacobwilliams/fmin (tag 1.1.1) | J. Williams | BSD-3-Clause | `src/moist/math/solver/fmin/LICENSE` |
| LSQR | `src/moist/math/linalg/lsqr/` | jacobwilliams/LSQR (tag 1.1.0) | M. Saunders (SOL, Stanford); modern Fortran by J. Williams | BSD-3-Clause (+ CPL-1.0 original, + LAPACK BSD) | `src/moist/math/linalg/lsqr/LICENSE` |
| LSMR | `src/moist/math/linalg/lsmr/` | jacobwilliams/LSMR (tag 1.0.0) | D. Fong & M. Saunders (SOL, Stanford) | BSD-2-Clause | `src/moist/math/linalg/lsmr/LICENSE` |
| LUSOL | `src/moist/math/linalg/lusol/` | jacobwilliams/lusol (tag 1.0.0) | Systems Optimization Laboratory, Stanford University | MIT OR BSD-3-Clause | `src/moist/math/linalg/lusol/LICENSE` |

They are unmodified except for two OpenMP fixes in FINUFFT's `spread.hpp` and `makeplan.hpp`.
See each directory's `PROVENANCE.md` for the exact file list, origin, pinned revision and local patches.

| Component | Path in `moist` | Upstream repository | Pinned revision | License | License file |
|---|---|---|---|---|---|
| FINUFFT (non-uniform FFT) | `src/moist/math/fft/finufft/` | flatironinstitute/finufft | commit `0e9c10409580c482656b3f840c9d875f98c3a066` (2.6.0-dev) | Apache-2.0 | `src/moist/math/fft/finufft/LICENSE`, `NOTICE` |
| xsimd (SIMD wrappers, header-only) | `src/moist/math/fft/xsimd/` | xtensor-stack/xsimd | commit `6842624fc8adafd7168a999e7150b384411da448` | BSD-3-Clause | `src/moist/math/fft/xsimd/LICENSE` |
| ducc0 FFT (16-file subset) | `src/moist/math/fft/ducc0/` | mreineck/ducc, tag `ducc0_0_39_1` (commit `0c0525521bc56d1b7bef8eb97a6debb6a3bfd33c`) | see tag | BSD-3-Clause (only dual-licensed files are vendored; moist takes the BSD-3-Clause option) | full text in the header of every file, also collected in `src/moist/math/fft/ducc0/LICENSE`; provenance and closure in `src/moist/math/fft/ducc0/PROVENANCE.md` |

## Bundled subprojects

These are managed via Meson wrap files under `subprojects/` and retain their own license
files in-tree (unmodified upstream). Listed here for completeness.

| Subproject | License file(s) | License |
|---|---|---|
| toml-f | `subprojects/toml-f/LICENSE-Apache`, `LICENSE-MIT` | Apache-2.0 OR MIT |
| jonquil | `subprojects/jonquil/LICENSE-Apache`, `LICENSE-MIT` | Apache-2.0 OR MIT |
| test-drive | `subprojects/test-drive/LICENSE-Apache`, `LICENSE-MIT` | Apache-2.0 OR MIT |
| mctc-lib | `subprojects/mctc-lib/LICENSE` | Apache-2.0 |
| mstore | `subprojects/mstore/LICENSE` | Apache-2.0 |
| fclap | `subprojects/fclap/LICENSE.md` | MIT (Christian Selzer) |

## Notes

- BSD-3-Clause and BSD-2-Clause permit redistribution with or without modification
  provided the copyright notice, conditions, and disclaimer are retained (source) and
  reproduced in accompanying materials (binary). This file, shipped with source and
  binary distributions, satisfies the binary-redistribution clause.
- Names of upstream authors/contributors are not used to endorse or promote `moist`.
