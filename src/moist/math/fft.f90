!> Fortran bindings for moist's FFT backend
!>
!> - Transforms come from ducc0, vendored under `fft/ducc0/`
!>   (BSD-3-Clause; see that directory's `PROVENANCE.md`), reached through
!>   the small `extern "C"` shim in `fft/shim.cpp`; nothing downloaded, no
!>   external FFT library needed
!> - Plan-free backend: ducc0 builds and caches its plans internally, keyed
!>   on the transform size; grid and trafo types carry no plan handles and
!>   no planning critical sections
!> - Unnormalised apart from the explicit `fct` factor; input left
!>   untouched, except the complex-to-real transforms (`moist_fft_c2r_3d`,
!>   `moist_fft_c2r_3d_batch`), which destroy their complex input, like
!>   FFTW's c2r
!> - Every entry point returns an `integer(c_int)` status: zero on success,
!>   nonzero if the ducc0 backend failed (allocation failure or an
!>   exception caught at the C++ boundary); callers must check it
module moist_math_fft
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
   implicit none(type, external)
   private

   public :: moist_fft_c2c_3d, moist_fft_r2c_3d, moist_fft_c2r_3d
   public :: moist_fft_c2c_3d_pass, moist_fft_c2c_3d_pass_inplace
   public :: moist_fft_r2c_3d_batch, moist_fft_c2r_3d_batch
   public :: moist_fft_dst4

   interface

      !> Batched spatial transform with explicit worker count; sites are independent
      !>
      !> @param[in] n0  Spatial extent, slowest-varying
      !> @param[in] n1  Spatial extent, middle
      !> @param[in] n2  Spatial extent, fastest-varying
      !> @param[in] nb  Number of sites
      !> @param[in] input  Contiguous input fields, preserved on return
      !> @param[out] output  Transformed fields
      !> @param[in] fct  Fourier normalization
      !> @param[in] nthreads  Worker count, at least one
      function moist_fft_r2c_3d_batch(n0, n1, n2, nb, input, output, fct, nthreads) &
         & result(status) bind(C, name="moist_fft_r2c_3d_batch")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Spatial extents and batch size
         integer(c_int), value, intent(in) :: n0, n1, n2, nb
         !> Input fields
         real(c_double), intent(in) :: input(*)
         !> Output fields
         complex(c_double_complex), intent(out) :: output(*)
         !> Normalization factor
         real(c_double), value, intent(in) :: fct
         !> Number of workers
         integer(c_int), value, intent(in) :: nthreads
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_r2c_3d_batch

      !> Batched spatial transform with explicit worker count; sites are independent
      !>
      !> @param[in] n0  Spatial extent, slowest-varying
      !> @param[in] n1  Spatial extent, middle
      !> @param[in] n2  Spatial extent, fastest-varying
      !> @param[in] nb  Number of sites
      !> @param[in,out] input  Contiguous input fields, destroyed on return
      !> @param[out] output  Transformed fields
      !> @param[in] fct  Fourier normalization
      !> @param[in] nthreads  Worker count, at least one
      function moist_fft_c2r_3d_batch(n0, n1, n2, nb, input, output, fct, nthreads) &
         & result(status) bind(C, name="moist_fft_c2r_3d_batch")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Spatial extents and batch size
         integer(c_int), value, intent(in) :: n0, n1, n2, nb
         !> Input fields, used as scratch
         complex(c_double_complex), intent(inout) :: input(*)
         !> Output fields
         real(c_double), intent(out) :: output(*)
         !> Normalization factor
         real(c_double), value, intent(in) :: fct
         !> Number of workers
         integer(c_int), value, intent(in) :: nthreads
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_c2r_3d_batch

      !> 3D complex-to-complex transform of a column-major field
      !>
      !> x fastest, i.e. logical extents (n0, n1, n2) = (nz, ny, nx)
      !>
      !> @param[in]  n0       Slowest-varying extent (z)
      !> @param[in]  n1       Middle extent (y)
      !> @param[in]  n2       Fastest-varying extent (x)
      !> @param[in]  in       Input field, n0*n1*n2 complex values
      !> @param[out] out      Output field, n0*n1*n2 complex values
      !> @param[in]  forward  Nonzero for the forward transform, exp(-i k x)
      !> @param[in]  fct      Factor applied to every output value
      function moist_fft_c2c_3d(n0, n1, n2, in, out, forward, fct) &
         & result(status) bind(C, name="moist_fft_c2c_3d")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Slowest-varying extent
         integer(c_int), value, intent(in) :: n0
         !> Middle extent
         integer(c_int), value, intent(in) :: n1
         !> Fastest-varying extent
         integer(c_int), value, intent(in) :: n2
         !> Input field
         complex(c_double_complex), dimension(*), intent(in) :: in
         !> Output field
         complex(c_double_complex), dimension(*), intent(out) :: out
         !> Transform direction
         integer(c_int), value, intent(in) :: forward
         !> Output scaling factor
         real(c_double), value, intent(in) :: fct
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_c2c_3d

      !> One 1D pass (along a single axis) of the 3D complex transform
      !>
      !> Over a sub-block of a larger column-major field with x fastest;
      !> `in` and `out` are the first element of the sub-block; values of
      !> the parent field outside it are neither read nor written
      !>
      !> @param[in]    m0       Sub-block extent along axis 0 (z)
      !> @param[in]    m1       Sub-block extent along axis 1 (y)
      !> @param[in]    m2       Sub-block extent along axis 2 (x)
      !> @param[in]    s0       Parent stride of axis 0, in elements
      !> @param[in]    s1       Parent stride of axis 1, in elements
      !> @param[in]    in       First input element of the sub-block
      !> @param[in,out] out     First output element of the sub-block
      !> @param[in]    axis     Transformed axis, 0 (z), 1 (y) or 2 (x)
      !> @param[in]    forward  Nonzero for the forward transform, exp(-i k x)
      !> @param[in]    fct      Factor applied to every transformed value
      function moist_fft_c2c_3d_pass(m0, m1, m2, s0, s1, in, out, axis, forward, &
         & fct) result(status) bind(C, name="moist_fft_c2c_3d_pass")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Sub-block extent along z
         integer(c_int), value, intent(in) :: m0
         !> Sub-block extent along y
         integer(c_int), value, intent(in) :: m1
         !> Sub-block extent along x
         integer(c_int), value, intent(in) :: m2
         !> Parent stride of z
         integer(c_int), value, intent(in) :: s0
         !> Parent stride of y
         integer(c_int), value, intent(in) :: s1
         !> First input element
         complex(c_double_complex), dimension(*), intent(in) :: in
         !> First output element
         complex(c_double_complex), dimension(*), intent(inout) :: out
         !> Transformed axis
         integer(c_int), value, intent(in) :: axis
         !> Transform direction
         integer(c_int), value, intent(in) :: forward
         !> Output scaling factor
         real(c_double), value, intent(in) :: fct
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_c2c_3d_pass

      !> In-place form of `moist_fft_c2c_3d_pass`
      !>
      !> @param[in]    m0       Sub-block extent along axis 0 (z)
      !> @param[in]    m1       Sub-block extent along axis 1 (y)
      !> @param[in]    m2       Sub-block extent along axis 2 (x)
      !> @param[in]    s0       Parent stride of axis 0, in elements
      !> @param[in]    s1       Parent stride of axis 1, in elements
      !> @param[in,out] data    First element of the sub-block
      !> @param[in]    axis     Transformed axis, 0 (z), 1 (y) or 2 (x)
      !> @param[in]    forward  Nonzero for the forward transform, exp(-i k x)
      !> @param[in]    fct      Factor applied to every transformed value
      function moist_fft_c2c_3d_pass_inplace(m0, m1, m2, s0, s1, data, axis, forward, &
         & fct) result(status) bind(C, name="moist_fft_c2c_3d_pass_inplace")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Sub-block extent along z
         integer(c_int), value, intent(in) :: m0
         !> Sub-block extent along y
         integer(c_int), value, intent(in) :: m1
         !> Sub-block extent along x
         integer(c_int), value, intent(in) :: m2
         !> Parent stride of z
         integer(c_int), value, intent(in) :: s0
         !> Parent stride of y
         integer(c_int), value, intent(in) :: s1
         !> First element of the sub-block
         complex(c_double_complex), dimension(*), intent(inout) :: data
         !> Transformed axis
         integer(c_int), value, intent(in) :: axis
         !> Transform direction
         integer(c_int), value, intent(in) :: forward
         !> Output scaling factor
         real(c_double), value, intent(in) :: fct
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_c2c_3d_pass_inplace

      !> 3D real-to-complex transform
      !>
      !> Half spectrum produced along the fastest-varying (x) axis, giving
      !> n0*n1*(n2/2+1) complex values
      !>
      !> @param[in]  n0   Slowest-varying extent (z)
      !> @param[in]  n1   Middle extent (y)
      !> @param[in]  n2   Fastest-varying extent (x)
      !> @param[in]  in   Real input field, n0*n1*n2 values
      !> @param[out] out  Complex output field, n0*n1*(n2/2+1) values
      !> @param[in]  fct  Factor applied to every output value
      function moist_fft_r2c_3d(n0, n1, n2, in, out, fct) &
         & result(status) bind(C, name="moist_fft_r2c_3d")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Slowest-varying extent
         integer(c_int), value, intent(in) :: n0
         !> Middle extent
         integer(c_int), value, intent(in) :: n1
         !> Fastest-varying extent
         integer(c_int), value, intent(in) :: n2
         !> Real input field
         real(c_double), dimension(*), intent(in) :: in
         !> Complex (half-spectrum) output field
         complex(c_double_complex), dimension(*), intent(out) :: out
         !> Output scaling factor
         real(c_double), value, intent(in) :: fct
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_r2c_3d

      !> 3D complex-to-real transform, the inverse of `moist_fft_r2c_3d`
      !>
      !> Input destroyed, used as the transform scratch, like FFTW's c2r
      !>
      !> @param[in]  n0   Slowest-varying extent (z)
      !> @param[in]  n1   Middle extent (y)
      !> @param[in]  n2   Fastest-varying extent (x)
      !> @param[in,out] in  Complex (half-spectrum) input, n0*n1*(n2/2+1) values, destroyed
      !> @param[out] out  Real output field, n0*n1*n2 values
      !> @param[in]  fct  Factor applied to every output value
      function moist_fft_c2r_3d(n0, n1, n2, in, out, fct) &
         & result(status) bind(C, name="moist_fft_c2r_3d")
         import :: c_int, c_double, c_double_complex
         implicit none(type, external)
         !> Slowest-varying extent
         integer(c_int), value, intent(in) :: n0
         !> Middle extent
         integer(c_int), value, intent(in) :: n1
         !> Fastest-varying extent
         integer(c_int), value, intent(in) :: n2
         !> Complex (half-spectrum) input field, destroyed
         complex(c_double_complex), dimension(*), intent(inout) :: in
         !> Real output field
         real(c_double), dimension(*), intent(out) :: out
         !> Output scaling factor
         real(c_double), value, intent(in) :: fct
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_c2r_3d

      !> Batched unnormalised DST-IV over the columns of a column-major array
      !>
      !>   y_k = fct * 2 * sum_j x_j sin(pi (j+1/2)(k+1/2) / npts)
      !>
      !> FFTW's `RODFT11` convention; its own inverse up to a factor of 2*npts
      !>
      !> @param[in]  npts    Transform length (column length)
      !> @param[in]  nbatch  Number of columns transformed together
      !> @param[in]  in      Input columns, npts*nbatch values
      !> @param[out] out     Transformed columns, npts*nbatch values
      !> @param[in]  fct     Factor applied to every output value
      function moist_fft_dst4(npts, nbatch, in, out, fct) &
         & result(status) bind(C, name="moist_fft_dst4")
         import :: c_int, c_double
         implicit none(type, external)
         !> Transform length
         integer(c_int), value, intent(in) :: npts
         !> Batch size
         integer(c_int), value, intent(in) :: nbatch
         !> Input columns
         real(c_double), dimension(*), intent(in) :: in
         !> Transformed columns
         real(c_double), dimension(*), intent(out) :: out
         !> Output scaling factor
         real(c_double), value, intent(in) :: fct
         !> Backend status, zero on success, nonzero if the backend failed
         integer(c_int) :: status
      end function moist_fft_dst4

   end interface

end module moist_math_fft
