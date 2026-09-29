# Vendored ducc0 FFT sources

Origin: **ducc0** (Distinctly Useful Code Collection) by Martin Reinecke: https://gitlab.mpcdf.mpg.de/mtr/ducc (https://github.com/mreineck/ducc)
Tag: `ducc0_0_39_1` (commit `0c05255`, "add missing explicit template instantiations").
Every file below is a byte-for-byte copy of `src/ducc0/<path>` at that tag; nothing was modified.

## License

ducc0 as a whole is GPL-2.0-or-later, but each file in this directory carries the header `SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later` and moist takes the **BSD-3-Clause** option.

Keep the headers of these files verbatim; the BSD text they carry is the license notice the BSD-3-Clause requires to be retained. The summary row in `THIRD_PARTY_LICENSES.md` points here.

## Files

Include directory: `src/moist/math/fft/` (so the sources' own `#include "ducc0/..."` paths resolve unchanged).

| Path | Compiled |
|---|---|
| `fft/fft.h` | header |
| `fft/fft1d_impl.h` | header |
| `fft/fftnd_impl.h` | header |
| `infra/aligned_array.h` | header |
| `infra/error_handling.h` | header |
| `infra/mav.cc` | yes (`moist_fft`) |
| `infra/mav.h` | header |
| `infra/misc_utils.h` | header |
| `infra/simd.h` | header |
| `infra/string_utils.cc` | yes (`moist_fft`) |
| `infra/string_utils.h` | header |
| `infra/threading.cc` | yes (`moist_fft`) |
| `infra/threading.h` | header |
| `infra/useful_macros.h` | header |
| `math/cmplx.h` | header |
| `math/unity_roots.h` | header |

The set is the transitive `#include "ducc0/..."` closure of `fft/fft.h` and `fft/fftnd_impl.h` plus the three translation units their symbols need (`threading.cc` for the thread pool, `mav.cc` for the array views, `string_utils.cc` for error messages).
`threading.h` also names `ducc0_custom_lowlevel_threading.h`, which is only reached behind `DUCC0_CUSTOM_LOWLEVEL_THREADING` and is not used here.