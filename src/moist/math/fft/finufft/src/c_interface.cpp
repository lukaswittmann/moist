// Modified for moist: double precision only, finufftf_* entry points removed
// ("moist patch"), see PROVENANCE.md
// public header
#include <finufft.h>
// private headers
#include <array>
#include <complex>
#include <finufft/plan.hpp> // (must come after complex.h)
#include <finufft_opts.h>

using namespace std;
using finufft::common::safe_finufft_call;

/* ---------------------------------------------------------------------------
   The 18 simple interfaces (= 3 dims * 3 types * {singlecall,many}) to FINUFFT.
   As of v1.2 these simply invoke the guru interface, through a helper layer.
   See ../docs/usage.rst or http://finufft.readthedocs.io for documentation
   all routines here.

   Authors: Andrea Malleo and Alex Barnett, 2019-2020.
   Safe namespacing, Barnett, May 2022.
   ---------------------------------------------------------------------------
*/

using f64   = double;
using cf64  = const double;
using c128  = std::complex<double>;
using cc128 = const std::complex<double>;
using i64   = BIGINT;
using ci64  = const BIGINT;

void finufft_default_opts(finufft_opts *o) { finufft_default_opts_t(o); }

int finufft_makeplan(int type, int dim, ci64 *n_modes, int iflag, int ntrans, f64 tol,
                     finufft_plan *pp, const finufft_opts *opts) {
  return safe_finufft_call([&]() -> int {
    return finufft_makeplan_t<f64>(type, dim, n_modes, iflag, ntrans, tol,
                                   reinterpret_cast<FINUFFT_PLAN_T<f64> **>(pp), opts);
  });
}

int finufft_setpts(finufft_plan p, i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, i64 nk, cf64 *s,
                   cf64 *t, cf64 *u) {
  return safe_finufft_call([&]() -> int {
    return reinterpret_cast<FINUFFT_PLAN_T<f64> *>(p)->setpts(nj, xj, yj, zj, nk, s, t,
                                                              u);
  });
}

int finufft_execute(finufft_plan p, c128 *cj, c128 *fk) {
  return safe_finufft_call([&]() -> int {
    return reinterpret_cast<FINUFFT_PLAN_T<f64> *>(p)->execute(cj, fk);
  });
}
int finufft_execute_adjoint(finufft_plan p, c128 *cj, c128 *fk) {
  return safe_finufft_call([&]() -> int {
    return reinterpret_cast<FINUFFT_PLAN_T<f64> *>(p)->execute_adjoint(cj, fk);
  });
}

int finufft_destroy(finufft_plan p)
// Free everything we allocated inside of finufft_plan pointed to by p.
// Also must not crash if called immediately after finufft_makeplan.
// Thus either each thing free'd here is guaranteed to be nullptr or correctly
// allocated.
{
  return safe_finufft_call([&]() -> int {
    if (!p) // nullptr ptr, so not a ptr to a plan, report error
      return 1;

    delete reinterpret_cast<FINUFFT_PLAN_T<f64> *>(p);
    return 0; // success
  });
}
// Helper layer ...........................................................

namespace { // helpers local to this TU
template<typename T>
int guru(int n_dims, int type, int n_transf, i64 nj, const std::array<const T *, 3> &xyz,
         std::complex<T> *cj, int iflag, T eps, const std::array<ci64, 3> &n_modes,
         i64 nk, const std::array<const T *, 3> &stu, std::complex<T> *fk,
         const finufft_opts *popts)
// Helper layer between simple interfaces (with opts) and the guru functions.
// Plan is stack-allocated since its lifetime is scoped to this call.
// Errors throw finufft::exception and are caught by safe_finufft_call.
// Author: Andrea Malleo, 2019. Stack alloc: Barbone, 2026.
{
  return safe_finufft_call([&]() -> int {
    FINUFFT_PLAN_T<T> plan(type, n_dims, n_modes.data(), iflag, n_transf, eps, popts);
    plan.setpts(nj, xyz[0], xyz[1], xyz[2], nk, stu[0], stu[1], stu[2]);
    plan.execute(cj, fk);
    return 0;
  });
}
template<typename T>
int guru13(int n_dims, int type, int n_transf, i64 nj,
           const std::array<const T *, 3> &xyz, const std::complex<T> *cj, int iflag,
           T eps, const std::array<ci64, 3> &n_modes, i64 nk,
           const std::array<const T *, 3> &stu, std::complex<T> *fk,
           const finufft_opts *popts) {
  return guru(n_dims, type, n_transf, nj, xyz, const_cast<std::complex<T> *>(cj), iflag,
              eps, n_modes, nk, stu, fk, popts);
}
template<typename T>
int guru2(int n_dims, int type, int n_transf, i64 nj, const std::array<const T *, 3> &xyz,
          std::complex<T> *cj, int iflag, T eps, const std::array<ci64, 3> &n_modes,
          i64 nk, const std::array<const T *, 3> &stu, const std::complex<T> *fk,
          const finufft_opts *popts) {
  return guru(n_dims, type, n_transf, nj, xyz, cj, iflag, eps, n_modes, nk, stu,
              const_cast<std::complex<T> *>(fk), popts);
}
} // anonymous namespace

// Dimension 1111111111111111111111111111111111111111111111111111111111111111

int finufft1d1many(int n_transf, i64 nj, cf64 *xj, cc128 *cj, int iflag, f64 eps, i64 ms,
                   c128 *fk, const finufft_opts *opts) {
  return guru13<f64>(1, 1, n_transf, nj, {xj}, cj, iflag, eps, {ms, 1, 1}, 0, {}, fk,
                     opts);
}

int finufft1d1(i64 nj, cf64 *xj, cc128 *cj, int iflag, f64 eps, i64 ms, c128 *fk,
               const finufft_opts *opts) {
  return finufft1d1many(1, nj, xj, cj, iflag, eps, ms, fk, opts);
}

int finufft1d2many(int n_transf, i64 nj, cf64 *xj, c128 *cj, int iflag, f64 eps, i64 ms,
                   cc128 *fk, const finufft_opts *opts) {
  return guru2<f64>(1, 2, n_transf, nj, {xj}, cj, iflag, eps, {ms, 1, 1}, 0, {}, fk,
                    opts);
}

int finufft1d2(i64 nj, cf64 *xj, c128 *cj, int iflag, f64 eps, i64 ms, cc128 *fk,
               const finufft_opts *opts) {
  return finufft1d2many(1, nj, xj, cj, iflag, eps, ms, fk, opts);
}

int finufft1d3many(int n_transf, i64 nj, cf64 *xj, cc128 *cj, int iflag, f64 eps, i64 nk,
                   cf64 *s, c128 *fk, const finufft_opts *opts) {
  return guru13<f64>(1, 3, n_transf, nj, {xj}, cj, iflag, eps, {0, 0, 0}, nk, {s}, fk,
                     opts);
}
int finufft1d3(i64 nj, cf64 *xj, cc128 *cj, int iflag, f64 eps, i64 nk, cf64 *s, c128 *fk,
               const finufft_opts *opts) {
  return finufft1d3many(1, nj, xj, cj, iflag, eps, nk, s, fk, opts);
}

// Dimension 22222222222222222222222222222222222222222222222222222222222222222

int finufft2d1many(int n_transf, i64 nj, cf64 *xj, cf64 *yj, cc128 *c, int iflag, f64 eps,
                   i64 ms, i64 mt, c128 *fk, const finufft_opts *opts) {
  return guru13<f64>(2, 1, n_transf, nj, {xj, yj}, c, iflag, eps, {ms, mt, 1}, 0, {}, fk,
                     opts);
}
int finufft2d1(i64 nj, cf64 *xj, cf64 *yj, cc128 *cj, int iflag, f64 eps, i64 ms, i64 mt,
               c128 *fk, const finufft_opts *opts) {
  return finufft2d1many(1, nj, xj, yj, cj, iflag, eps, ms, mt, fk, opts);
}

int finufft2d2many(int n_transf, i64 nj, cf64 *xj, cf64 *yj, c128 *c, int iflag, f64 eps,
                   i64 ms, i64 mt, cc128 *fk, const finufft_opts *opts) {
  return guru2<f64>(2, 2, n_transf, nj, {xj, yj}, c, iflag, eps, {ms, mt, 1}, 0, {}, fk,
                    opts);
}
int finufft2d2(i64 nj, cf64 *xj, cf64 *yj, c128 *cj, int iflag, f64 eps, i64 ms, i64 mt,
               cc128 *fk, const finufft_opts *opts) {
  return finufft2d2many(1, nj, xj, yj, cj, iflag, eps, ms, mt, fk, opts);
}

int finufft2d3many(int n_transf, i64 nj, cf64 *xj, cf64 *yj, cc128 *cj, int iflag,
                   f64 eps, i64 nk, cf64 *s, cf64 *t, c128 *fk,
                   const finufft_opts *opts) {
  return guru13<f64>(2, 3, n_transf, nj, {xj, yj}, cj, iflag, eps, {0, 0, 0}, nk, {s, t},
                     fk, opts);
}
int finufft2d3(i64 nj, cf64 *xj, cf64 *yj, cc128 *cj, int iflag, f64 eps, i64 nk, cf64 *s,
               cf64 *t, c128 *fk, const finufft_opts *opts) {
  return finufft2d3many(1, nj, xj, yj, cj, iflag, eps, nk, s, t, fk, opts);
}

// Dimension 3333333333333333333333333333333333333333333333333333333333333333

int finufft3d1many(int n_transf, i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, cc128 *cj,
                   int iflag, f64 eps, i64 ms, i64 mt, i64 mu, c128 *fk,
                   const finufft_opts *opts) {
  return guru13<f64>(3, 1, n_transf, nj, {xj, yj, zj}, cj, iflag, eps, {ms, mt, mu}, 0,
                     {}, fk, opts);
}
int finufft3d1(i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, cc128 *cj, int iflag, f64 eps,
               i64 ms, i64 mt, i64 mu, c128 *fk, const finufft_opts *opts) {
  return finufft3d1many(1, nj, xj, yj, zj, cj, iflag, eps, ms, mt, mu, fk, opts);
}

int finufft3d2many(int n_transf, i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, c128 *cj,
                   int iflag, f64 eps, i64 ms, i64 mt, i64 mu, cc128 *fk,
                   const finufft_opts *opts) {
  return guru2<f64>(3, 2, n_transf, nj, {xj, yj, zj}, cj, iflag, eps, {ms, mt, mu}, 0, {},
                    fk, opts);
}
int finufft3d2(i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, c128 *cj, int iflag, f64 eps, i64 ms,
               i64 mt, i64 mu, cc128 *fk, const finufft_opts *opts) {
  return finufft3d2many(1, nj, xj, yj, zj, cj, iflag, eps, ms, mt, mu, fk, opts);
}

int finufft3d3many(int n_transf, i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, cc128 *cj,
                   int iflag, f64 eps, i64 nk, cf64 *s, cf64 *t, cf64 *u, c128 *fk,
                   const finufft_opts *opts) {
  return guru13<f64>(3, 3, n_transf, nj, {xj, yj, zj}, cj, iflag, eps, {0, 0, 0}, nk,
                     {s, t, u}, fk, opts);
}
int finufft3d3(i64 nj, cf64 *xj, cf64 *yj, cf64 *zj, cc128 *cj, int iflag, f64 eps,
               i64 nk, cf64 *s, cf64 *t, cf64 *u, c128 *fk, const finufft_opts *opts) {
  return finufft3d3many(1, nj, xj, yj, zj, cj, iflag, eps, nk, s, t, u, fk, opts);
}
