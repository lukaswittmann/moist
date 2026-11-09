// Modified for moist: double precision only, finufftf_*_ entry points removed
// ("moist patch"), see PROVENANCE.md
/* C++ layer for calling FINUFFT from fortran, in f77 style + derived type
   for the finufft_opts C-struct. The ptr to finufft_plan is passed as an "opaque"
   pointer (as in FFTW3 legacy fortran interface).
   Note the trailing underscore name-mangle which is not needed from fortran.

   Note our typedefs:
   FLT = double (or float, depending on compilation precision)
   CPX = double complex (or float complex, depending on compilation precision)
   BIGINT = int64 (integer*8 in fortran)

   Make sure you call this library with matching fortran types

   For a demo see: examples/simple1d1.f

   Barnett 2/17/17. Single prec 4/5/17. Libin Lu & Alex Barnett, May 2020.
   Garrett Wright dual-prec 6/28/20.
   Barnett safe-header, cleaned up, 6/7/22. execute_adjoint Barnett 6/26/25.
*/

// public header
#include <finufft.h>
#include <finufft/plan.hpp>

using f64  = double;
using c128 = std::complex<double>;
using i64  = BIGINT;

#ifdef __cplusplus
extern "C" {
#endif

// --------------------- guru interface from fortran ------------------------
FINUFFT_EXPORT
void finufft_makeplan_(int *type, int *n_dims, i64 *n_modes, int *iflag, int *n_transf,
                       f64 *tol, finufft_plan *plan, finufft_opts *o, int *ier) {
  if (!plan)
    fprintf(stderr,
            "%s fortran: plan must be allocated as at least the size of a C pointer "
            "(usually 8 bytes)!\n",
            __func__);
  else {
    // pass o whether it's a NULL or pointer to a fortran-allocated finufft_opts:
    *ier = finufft_makeplan(*type, *n_dims, n_modes, *iflag, *n_transf, *tol, plan, o);
  }
}

FINUFFT_EXPORT
void finufft_setpts_(finufft_plan *plan, i64 *M, f64 *xj, f64 *yj, f64 *zj, i64 *nk,
                     f64 *s, f64 *t, f64 *u, int *ier) {
  if (!*plan) {
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
    return;
  }
  int nk_safe = 0; // catches the case where user passes NULL in
  if (nk) nk_safe = int(*nk);
  *ier = finufft_setpts(*plan, *M, xj, yj, zj, nk_safe, s, t, u);
}

FINUFFT_EXPORT
void finufft_execute_(finufft_plan *plan, c128 *weights, c128 *result, int *ier) {
  if (!plan)
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
  else
    *ier = finufft_execute(*plan, weights, result);
}

FINUFFT_EXPORT
void finufft_execute_adjoint_(finufft_plan *plan, c128 *weights, c128 *result, int *ier) {
  if (!plan)
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
  else
    *ier = finufft_execute_adjoint(*plan, weights, result);
}

FINUFFT_EXPORT
void finufft_destroy_(finufft_plan *plan, int *ier) {
  if (!plan)
    fprintf(stderr, "%s fortran: finufft_plan unallocated!", __func__);
  else
    *ier = finufft_destroy(*plan);
}

// ------------ use FINUFFT to set the default options ---------------------
// (Note the finufft_opts is created in f90-style derived types, not here)
FINUFFT_EXPORT
void finufft_default_opts_(finufft_opts *o) {
  if (!o)
    fprintf(stderr, "%s fortran: opts must be allocated!\n", __func__);
  else
    // o is a ptr to already-allocated fortran finufft_opts derived type...
    finufft_default_opts(o);
}

// -------------- simple and many-vector interfaces --------------------
// --- 1D ---
FINUFFT_EXPORT
void finufft1d1_(i64 *nj, f64 *xj, c128 *cj, int *iflag, f64 *eps, i64 *ms, c128 *fk,
                 finufft_opts *o, int *ier) {
  *ier = finufft1d1(*nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

FINUFFT_EXPORT
void finufft1d1many_(int *ntransf, i64 *nj, f64 *xj, c128 *cj, int *iflag, f64 *eps,
                     i64 *ms, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft1d1many(*ntransf, *nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

FINUFFT_EXPORT
void finufft1d2_(i64 *nj, f64 *xj, c128 *cj, int *iflag, f64 *eps, i64 *ms, c128 *fk,
                 finufft_opts *o, int *ier) {
  *ier = finufft1d2(*nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

FINUFFT_EXPORT
void finufft1d2many_(int *ntransf, i64 *nj, f64 *xj, c128 *cj, int *iflag, f64 *eps,
                     i64 *ms, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft1d2many(*ntransf, *nj, xj, cj, *iflag, *eps, *ms, fk, o);
}

FINUFFT_EXPORT
void finufft1d3_(i64 *nj, f64 *x, c128 *c, int *iflag, f64 *eps, i64 *nk, f64 *s, c128 *f,
                 finufft_opts *o, int *ier) {
  *ier = finufft1d3(*nj, x, c, *iflag, *eps, *nk, s, f, o);
}

FINUFFT_EXPORT
void finufft1d3many_(int *ntransf, i64 *nj, f64 *x, c128 *c, int *iflag, f64 *eps,
                     i64 *nk, f64 *s, c128 *f, finufft_opts *o, int *ier) {
  *ier = finufft1d3many(*ntransf, *nj, x, c, *iflag, *eps, *nk, s, f, o);
}

// --- 2D ---
FINUFFT_EXPORT
void finufft2d1_(i64 *nj, f64 *xj, f64 *yj, c128 *cj, int *iflag, f64 *eps, i64 *ms,
                 i64 *mt, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft2d1(*nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}
FINUFFT_EXPORT
void finufft2d1many_(int *ntransf, i64 *nj, f64 *xj, f64 *yj, c128 *cj, int *iflag,
                     f64 *eps, i64 *ms, i64 *mt, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft2d1many(*ntransf, *nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}

FINUFFT_EXPORT
void finufft2d2_(i64 *nj, f64 *xj, f64 *yj, c128 *cj, int *iflag, f64 *eps, i64 *ms,
                 i64 *mt, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft2d2(*nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}
FINUFFT_EXPORT
void finufft2d2many_(int *ntransf, i64 *nj, f64 *xj, f64 *yj, c128 *cj, int *iflag,
                     f64 *eps, i64 *ms, i64 *mt, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft2d2many(*ntransf, *nj, xj, yj, cj, *iflag, *eps, *ms, *mt, fk, o);
}

FINUFFT_EXPORT
void finufft2d3_(i64 *nj, f64 *x, f64 *y, c128 *c, int *iflag, f64 *eps, i64 *nk, f64 *s,
                 f64 *t, c128 *f, finufft_opts *o, int *ier) {
  *ier = finufft2d3(*nj, x, y, c, *iflag, *eps, *nk, s, t, f, o);
}

FINUFFT_EXPORT
void finufft2d3many_(int *ntransf, i64 *nj, f64 *x, f64 *y, c128 *c, int *iflag, f64 *eps,
                     i64 *nk, f64 *s, f64 *t, c128 *f, finufft_opts *o, int *ier) {
  *ier = finufft2d3many(*ntransf, *nj, x, y, c, *iflag, *eps, *nk, s, t, f, o);
}

// --- 3D ---
FINUFFT_EXPORT
void finufft3d1_(i64 *nj, f64 *xj, f64 *yj, f64 *zj, c128 *cj, int *iflag, f64 *eps,
                 i64 *ms, i64 *mt, i64 *mu, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft3d1(*nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

FINUFFT_EXPORT
void finufft3d1many_(int *ntransf, i64 *nj, f64 *xj, f64 *yj, f64 *zj, c128 *cj,
                     int *iflag, f64 *eps, i64 *ms, i64 *mt, i64 *mu, c128 *fk,
                     finufft_opts *o, int *ier) {
  *ier =
      finufft3d1many(*ntransf, *nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

FINUFFT_EXPORT
void finufft3d2_(i64 *nj, f64 *xj, f64 *yj, f64 *zj, c128 *cj, int *iflag, f64 *eps,
                 i64 *ms, i64 *mt, i64 *mu, c128 *fk, finufft_opts *o, int *ier) {
  *ier = finufft3d2(*nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

FINUFFT_EXPORT
void finufft3d2many_(int *ntransf, i64 *nj, f64 *xj, f64 *yj, f64 *zj, c128 *cj,
                     int *iflag, f64 *eps, i64 *ms, i64 *mt, i64 *mu, c128 *fk,
                     finufft_opts *o, int *ier) {
  *ier =
      finufft3d2many(*ntransf, *nj, xj, yj, zj, cj, *iflag, *eps, *ms, *mt, *mu, fk, o);
}

FINUFFT_EXPORT
void finufft3d3_(i64 *nj, f64 *x, f64 *y, f64 *z, c128 *c, int *iflag, f64 *eps, i64 *nk,
                 f64 *s, f64 *t, f64 *u, c128 *f, finufft_opts *o, int *ier) {
  *ier = finufft3d3(*nj, x, y, z, c, *iflag, *eps, *nk, s, t, u, f, o);
}

FINUFFT_EXPORT
void finufft3d3many_(int *ntransf, i64 *nj, f64 *x, f64 *y, f64 *z, c128 *c, int *iflag,
                     f64 *eps, i64 *nk, f64 *s, f64 *t, f64 *u, c128 *f, finufft_opts *o,
                     int *ier) {
  *ier = finufft3d3many(*ntransf, *nj, x, y, z, c, *iflag, *eps, *nk, s, t, u, f, o);
}

#ifdef __cplusplus
}
#endif
