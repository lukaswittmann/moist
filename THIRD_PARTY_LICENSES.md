# Third-party code and licenses

`moist` is distributed under the **LGPL-3.0-or-later** license (see `COPYING` and
`COPYING.LESSER`). It incorporates third-party code under other, compatible licenses.
This file documents that code, its origin, and the applicable license, as required by
those licenses (in particular the BSD source- and binary-redistribution clauses).

All third-party components below are under licenses compatible with `moist`'s
LGPL-3.0-or-later license (BSD-2-Clause, BSD-3-Clause, MIT, Apache-2.0, or LGPL-3.0).
Each component remains under its own terms; its license is listed alongside it below.

## Vendored code (copied into the source tree and locally modified)

These numerical components were copied from upstream repositories — primarily by
Jacob Williams (https://github.com/jacobwilliams) — and then modified for use in `moist`.
Because they were modified they are not consumed as upstream submodules / Meson wraps.
The verbatim upstream license text accompanies each component at the path shown.

| Component | Path in `moist` | Upstream repository | Original author(s) | License | License file |
|---|---|---|---|---|---|
| SLSQP (+ BVLS) | `src/moist/math/solver/slsqp/` | jacobwilliams/slsqp | Dieter Kraft (1988); ACM (1994); BVLS: Lawson & Hanson (netlib, public domain); modern Fortran by J. Williams | BSD-3-Clause (+ MIT-style original grant) | `src/moist/math/solver/slsqp/LICENSE` |
| L-BFGS-B | `src/moist/math/solver/lbfgsb/` | jacobwilliams/lbfgsb | L-BFGS-B 3.0: J. Nocedal & J. L. Morales; modern Fortran by J. Williams | BSD-3-Clause ("New BSD") | `src/moist/math/solver/lbfgsb/LICENSE` |
| nlesolver | `src/moist/math/solver/newton/` | jacobwilliams/nlesolver-fortran | J. Williams | BSD-3-Clause | `src/moist/math/solver/newton/LICENSE` |
| fmin | `src/moist/math/solver/fmin/` | jacobwilliams/fmin (tag 1.1.1) | J. Williams | BSD-3-Clause | `src/moist/math/solver/fmin/LICENSE` |
| LSQR | `src/moist/math/linalg/lsqr/` | jacobwilliams/LSQR (tag 1.1.0) | M. Saunders (SOL, Stanford); modern Fortran by J. Williams | BSD-3-Clause (+ CPL-1.0 original, + LAPACK BSD) | `src/moist/math/linalg/lsqr/LICENSE` |
| LSMR | `src/moist/math/linalg/lsmr/` | jacobwilliams/LSMR (tag 1.0.0) | D. Fong & M. Saunders (SOL, Stanford) | BSD-2-Clause | `src/moist/math/linalg/lsmr/LICENSE` |
| LUSOL | `src/moist/math/linalg/lusol/` | jacobwilliams/lusol (tag 1.0.0) | Systems Optimization Laboratory, Stanford University | MIT OR BSD-3-Clause | `src/moist/math/linalg/lusol/LICENSE` |
| NumDiff | `src/moist/math/numdiff/` | jacobwilliams/NumDiff | J. Williams; DSM kernel from MINPACK (University of Chicago / Argonne) | BSD-3-Clause (+ MINPACK notice) | `src/moist/math/numdiff/LICENSE` |

## Bundled subprojects

These are managed via Meson wrap files under `subprojects/` and retain their own license
files in-tree (unmodified upstream). Listed here for completeness.

| Subproject | License file(s) | License |
|---|---|---|
| json-fortran | `subprojects/json-fortran-8.2.5/LICENSE` | BSD-3-Clause (Jacob Williams) |
| dftd4 | `subprojects/dftd4/COPYING(.LESSER)` | LGPL-3.0-or-later |
| s-dftd3 | `subprojects/s-dftd3/COPYING(.LESSER)` | LGPL-3.0-or-later |
| toml-f | `subprojects/toml-f/LICENSE-Apache`, `LICENSE-MIT` | Apache-2.0 OR MIT |
| jonquil | `subprojects/jonquil/LICENSE-Apache`, `LICENSE-MIT` | Apache-2.0 OR MIT |
| test-drive | `subprojects/test-drive/LICENSE-Apache`, `LICENSE-MIT` | Apache-2.0 OR MIT |
| mctc-lib | `subprojects/mctc-lib/LICENSE` | see file |
| mstore | `subprojects/mstore/LICENSE` | see file |
| multicharge | `subprojects/multicharge/LICENSE` | see file |
| fclap | `subprojects/fclap/LICENSE.md` | see file |

## Notes

- BSD-3-Clause and BSD-2-Clause are GPL/LGPL-compatible; they permit redistribution with
  or without modification provided the copyright notice, conditions, and disclaimer are
  retained (source) and reproduced in accompanying materials (binary). This file, shipped
  with source and binary distributions, satisfies the binary-redistribution clause.
- Names of upstream authors/contributors are not used to endorse or promote `moist`.
