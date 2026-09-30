// extern "C" shim exposing the ducc0 FFT engine to Fortran.
//
// ducc0's interface is a C++ template API taking multidimensional array views
// (`cfmav`/`vfmav`), so it cannot be called from Fortran directly. This file is
// the whole of moist's FFT backend surface: `moist_math_fft` (fft.f90) binds
// these entry points and nothing else.
//
// Two properties of ducc0 shape the design:
//
//   * It is plan-free at the call site. Plans are built and cached internally,
//     keyed on the transform size, so there is no plan object to create, share
//     between threads or destroy -- which is why the grid types carry no plan
//     handles.
//   * Views are constructed from a raw pointer plus explicit shape and strides.
//     moist's fields are Fortran column-major with x fastest, which is the same
//     memory order as a C row-major array with x last, so every view below is
//     built with shape {n0, n1, n2} and strides {n1*n2, n2, 1} where n2 is the
//     fastest-varying (x) axis. The strides are passed explicitly rather than
//     inferred, so the layout is stated rather than assumed.
//
// Single-field entry points keep `nthreads` at 1 for callers already using
// outer parallelism (including the 6D transforms). The batched real transforms
// used by 3D RISM distribute spatial lines across the OpenMP team instead.
// Inner FFT threading may change rounding because scalar and SIMD line kernels
// are not bitwise identical. Numerical agreement is tested across thread counts.
//
// Every entry point below returns an int status (0 = success, 1 = a negative
// extent or a ducc0 exception other than bad_alloc, 2 = std::bad_alloc)
// instead of void. ducc0
// throws on allocation failure and on backend errors, and a C++ exception
// unwinding through the Fortran frames that call these functions would abort
// the process rather than give the caller a chance to report it, so no
// exception may ever escape an `extern "C"` function here.

#include <complex>
#include <cstddef>
#include <algorithm>
#include <exception>
#include <initializer_list>
#include <new>
#include <functional>
#include <utility>
#include <mutex>
#include <condition_variable>
#include <deque>
#ifdef _OPENMP
#include <omp.h>
#endif

// fft.h declares the transforms and re-exports them into namespace ducc0;
// fftnd_impl.h carries their definitions. Including the implementation header
// instantiates the `double` specialisations we use directly in this
// translation unit, which is deliberate: ducc0 also ships them prebuilt in
// src/ducc0/fft/fft_inst{1,2}.cc, but those two stubs are the only files in
// the FFT dependency closure WITHOUT the
// "SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later" header that lets
// moist take the BSD terms. Instantiating here keeps every ducc0 file moist
// compiles or includes explicitly dual-licensed.
#include "ducc0/fft/fft.h"
#include "ducc0/fft/fftnd_impl.h"

namespace {

using cplx = std::complex<double>;
using ducc0::cfmav;
using ducc0::vfmav;
// shape_t / stride_t are members of fmav_info, not free names in ducc0.
using shape_t = ducc0::fmav_info::shape_t;
using stride_t = ducc0::fmav_info::stride_t;

#ifdef _OPENMP
//! Adapt ducc0 work submission to the existing OpenMP runtime. This avoids a
//! second pool inheriting the invoking thread's single-core affinity mask.
class OmpPool final : public ducc0::detail_threading::thread_pool
{
   std::size_t count_;
   std::mutex mutex_;
   std::condition_variable ready_;
   std::deque<std::function<void()>> queue_;
   bool stopped_ = false;
public:
   explicit OmpPool(std::size_t count) : count_(count) {}
   void set_count(std::size_t count) { count_ = count; }
   std::size_t nthreads() const override { return count_; }
   std::size_t adjust_nthreads(std::size_t requested) const override
   { return requested == 0 ? count_ : std::min(count_, requested); }
   void submit(std::function<void()> work) override
   {
      {
         std::lock_guard<std::mutex> lock(mutex_);
         queue_.push_back(std::move(work));
      }
      ready_.notify_one();
   }
   //! Workers must execute submissions immediately: ducc0 waits on its own
   //! condition variable, which is not an OpenMP task scheduling point.
   void work_loop()
   {
      for (;;) {
         std::function<void()> work;
         {
            std::unique_lock<std::mutex> lock(mutex_);
            ready_.wait(lock, [&] { return stopped_ || !queue_.empty(); });
            if (queue_.empty()) return;
            work = std::move(queue_.front());
            queue_.pop_front();
         }
         work();
      }
   }
   void stop()
   {
      {
         std::lock_guard<std::mutex> lock(mutex_);
         stopped_ = true;
      }
      ready_.notify_all();
   }
};
#endif

//! Status codes returned by every `extern "C"` entry point below.
enum status_code : int { status_ok = 0, status_error = 1, status_bad_alloc = 2 };

//! Translate the exception in flight into a status code. Must only be called
//! from a `catch` block.
int status_from_exception() noexcept
{
   try { throw; }
   catch (const std::bad_alloc &) { return status_bad_alloc; }
   catch (...) { return status_error; }
}

//! True if any extent is negative; a negative int would wrap to a huge
//! size_t and address memory far outside the caller's buffers
inline bool any_negative(std::initializer_list<int> extents) noexcept
{
   for (const int n : extents)
      if (n < 0) return true;
   return false;
}

//! Execute one batched transform with a bounded, affinity-aware OpenMP team.
//! Existing outer regions and builds without OpenMP retain a serial transform.
//! Returns a status code; no exception escapes this function, threaded or not.
template<class Transform> int with_fft_threads(int nthreads, Transform transform)
{
#ifdef _OPENMP
   if (nthreads > 1 && !omp_in_parallel()) {
      std::exception_ptr failure;
      OmpPool pool(nthreads);
#pragma omp parallel num_threads(nthreads) shared(failure, pool)
      {
#pragma omp single
         { pool.set_count(omp_get_num_threads()); }
         if (omp_get_thread_num() == 0) {
            ducc0::ScopedUseThreadPool guard(pool);
            try { transform(pool.nthreads()); }
            catch (...) { failure = std::current_exception(); }
            pool.stop();
         } else {
            pool.work_loop();
         }
      }
      if (failure) {
         try { std::rethrow_exception(failure); }
         catch (...) { return status_from_exception(); }
      }
      return status_ok;
   }
#else
   (void)nthreads;
#endif
   try {
      transform(1);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

//! Row-major strides for a rank-3 array of extents (n0, n1, n2).
inline stride_t strides_3d(std::size_t n0, std::size_t n1, std::size_t n2)
{
   (void)n0;
   return stride_t{static_cast<std::ptrdiff_t>(n1 * n2),
                   static_cast<std::ptrdiff_t>(n2),
                   static_cast<std::ptrdiff_t>(1)};
}

//! One 1D pass of a 3D transform: axis `axis` of the (m0, m1, m2) sub-block
//! starting at `in`/`out`, inside a parent field whose axis-0 and axis-1 strides
//! are `s0` and `s1` elements (axis 2 always has stride 1).  A single-axis call
//! is exactly one pass of ducc0's multi-axis loop, so chaining three of them in
//! the order ducc0 would pick reproduces a whole-block transform.
void c2c_pass(int m0, int m1, int m2, int s0, int s1, const double *in,
              double *out, int axis, int forward, double fct)
{
   const shape_t shape{static_cast<std::size_t>(m0), static_cast<std::size_t>(m1),
                       static_cast<std::size_t>(m2)};
   const stride_t str{static_cast<std::ptrdiff_t>(s0),
                      static_cast<std::ptrdiff_t>(s1), static_cast<std::ptrdiff_t>(1)};

   cfmav<cplx> mi(reinterpret_cast<const cplx *>(in), shape, str);
   vfmav<cplx> mo(reinterpret_cast<cplx *>(out), shape, str);
   ducc0::c2c(mi, mo, shape_t{static_cast<std::size_t>(axis)}, forward != 0, fct, 1);
}

} // namespace

extern "C" {

//! Out-of-place 3D complex-to-complex transform, unnormalised unless `fct`
//! says otherwise. `forward` selects the sign of the exponent (nonzero =
//! forward, exp(-i k x)). The Fortran binding declares `in` and `out` as
//! distinct dummies, so a Fortran caller may not pass one array as both.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_c2c_3d(int n0, int n1, int n2, const double *in, double *out,
                      int forward, double fct)
{
   if (any_negative({n0, n1, n2})) return status_error;
   try {
      const std::size_t s0 = static_cast<std::size_t>(n0);
      const std::size_t s1 = static_cast<std::size_t>(n1);
      const std::size_t s2 = static_cast<std::size_t>(n2);
      const shape_t shape{s0, s1, s2};
      const stride_t str = strides_3d(s0, s1, s2);

      cfmav<cplx> mi(reinterpret_cast<const cplx *>(in), shape, str);
      vfmav<cplx> mo(reinterpret_cast<cplx *>(out), shape, str);
      ducc0::c2c(mi, mo, shape_t{0, 1, 2}, forward != 0, fct, 1);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

//! Single-axis pass of a 3D complex transform over a sub-block, out of place.
//! `fct` multiplies each value once the pass has finished it. Negative
//! strides are valid ducc0 views; only the extents are checked.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_c2c_3d_pass(int m0, int m1, int m2, int s0, int s1,
                           const double *in, double *out, int axis, int forward,
                           double fct)
{
   if (any_negative({m0, m1, m2})) return status_error;
   try {
      c2c_pass(m0, m1, m2, s0, s1, in, out, axis, forward, fct);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

//! Single-axis pass of a 3D complex transform over a sub-block, in place.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_c2c_3d_pass_inplace(int m0, int m1, int m2, int s0, int s1,
                                   double *data, int axis, int forward,
                                   double fct)
{
   if (any_negative({m0, m1, m2})) return status_error;
   try {
      c2c_pass(m0, m1, m2, s0, s1, data, data, axis, forward, fct);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

//! 3D real-to-complex transform. The half spectrum is produced along the last
//! (fastest, x) axis, so `out` holds n0*n1*(n2/2+1) complex values. `in` is
//! left untouched.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_r2c_3d(int n0, int n1, int n2, const double *in, double *out,
                      double fct)
{
   if (any_negative({n0, n1, n2})) return status_error;
   try {
      const std::size_t s0 = static_cast<std::size_t>(n0);
      const std::size_t s1 = static_cast<std::size_t>(n1);
      const std::size_t s2 = static_cast<std::size_t>(n2);
      const std::size_t s2h = s2 / 2 + 1;

      cfmav<double> mi(in, shape_t{s0, s1, s2}, strides_3d(s0, s1, s2));
      vfmav<cplx> mo(reinterpret_cast<cplx *>(out), shape_t{s0, s1, s2h},
                     strides_3d(s0, s1, s2h));
      // The real transform runs over the last axis listed, matching FFTW's
      // r2c convention of halving the fastest-varying dimension.
      ducc0::r2c(mi, mo, shape_t{0, 1, 2}, true, fct, 1);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

//! 3D complex-to-real transform, the inverse of `moist_fft_r2c_3d`. Like
//! FFTW's multidimensional c2r, the input is DESTROYED: the two slow axes are
//! transformed in place in `in` before the final c2r pass. ducc0's const c2r
//! runs the same passes on a freshly allocated copy, so the result is
//! bitwise identical while the per-call scratch allocation disappears.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_c2r_3d(int n0, int n1, int n2, double *in, double *out,
                      double fct)
{
   if (any_negative({n0, n1, n2})) return status_error;
   try {
      const std::size_t s0 = static_cast<std::size_t>(n0);
      const std::size_t s1 = static_cast<std::size_t>(n1);
      const std::size_t s2 = static_cast<std::size_t>(n2);
      const std::size_t s2h = s2 / 2 + 1;

      vfmav<cplx> mi(reinterpret_cast<cplx *>(in), shape_t{s0, s1, s2h},
                     strides_3d(s0, s1, s2h));
      vfmav<double> mo(out, shape_t{s0, s1, s2}, strides_3d(s0, s1, s2));
      ducc0::c2r_mut(mi, mo, shape_t{0, 1, 2}, false, fct, 1);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

//! Batched 3D transforms: the site dimension is not transformed. One thread
//! pool distributes FFT lines across sites AND spatial dimensions. The caller
//! passes one thread inside an existing OpenMP region to avoid nested teams.
//! The views are built inside the transform so their allocations sit under
//! the exception guard of `with_fft_threads`.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_r2c_3d_batch(int n0, int n1, int n2, int nb, const double *in,
                           double *out, double fct, int nthreads)
{
   if (any_negative({n0, n1, n2, nb})) return status_error;
   const std::size_t z=n0, y=n1, x=n2, b=nb, h=x/2+1;
   return with_fft_threads(nthreads, [&](std::size_t nt) {
      cfmav<double> mi(in, shape_t{b,z,y,x});
      vfmav<cplx> mo(reinterpret_cast<cplx *>(out), shape_t{b,z,y,h});
      ducc0::r2c(mi, mo, shape_t{1,2,3}, true, fct, nt);
   });
}

//! Inverse batched transform; the input is destroyed (see `moist_fft_c2r_3d`)
//! and fct supplies normalization.
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_c2r_3d_batch(int n0, int n1, int n2, int nb, double *in,
                           double *out, double fct, int nthreads)
{
   if (any_negative({n0, n1, n2, nb})) return status_error;
   const std::size_t z=n0, y=n1, x=n2, b=nb, h=x/2+1;
   return with_fft_threads(nthreads, [&](std::size_t nt) {
      vfmav<cplx> mi(reinterpret_cast<cplx *>(in), shape_t{b,z,y,h});
      vfmav<double> mo(out, shape_t{b,z,y,x});
      ducc0::c2r_mut(mi, mo, shape_t{1,2,3}, false, fct, nt);
   });
}

//! Batched unnormalised DST-IV over `nbatch` columns of length `npts`.
//!
//! The Fortran caller holds a column-major (npts, nbatch) array, which is a
//! row-major (nbatch, npts) array in memory, so the transform runs over axis 1.
//! ducc0's type-4 DST with `ortho = false` is FFTW's `RODFT11`:
//!   y_k = 2 sum_j x_j sin(pi (j+1/2)(k+1/2) / npts).
//! @return status  Zero on success, nonzero on a negative extent or if the backend failed
int moist_fft_dst4(int npts, int nbatch, const double *in, double *out,
                    double fct)
{
   if (any_negative({npts, nbatch})) return status_error;
   try {
      const std::size_t n = static_cast<std::size_t>(npts);
      const std::size_t nb = static_cast<std::size_t>(nbatch);
      const shape_t shape{nb, n};
      const stride_t str{static_cast<std::ptrdiff_t>(n), static_cast<std::ptrdiff_t>(1)};

      cfmav<double> mi(in, shape, str);
      vfmav<double> mo(out, shape, str);
      ducc0::dst(mi, mo, shape_t{1}, 4, fct, false, 1);
   } catch (...) {
      return status_from_exception();
   }
   return status_ok;
}

} // extern "C"
