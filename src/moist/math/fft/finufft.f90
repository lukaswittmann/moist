!> Fortran bindings for the guru interface of the vendored FINUFFT
!>
!> - C entry points of `fft/finufft/fortran/finufftfort.cpp`, bound by symbol
!>   name; implicit externals mangle per compiler (ifx on Windows emits
!>   `FINUFFT_MAKEPLAN`)
!> - Options are assumed type (`type(*)`), passed by reference, since the
!>   vendored `finufft_opts` is not `bind(c)`
module moist_math_fft_finufft
   use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_double, c_double_complex
   implicit none(type, external)
   private

   public :: finufft_default_opts, finufft_makeplan, finufft_setpts, finufft_execute
   public :: finufft_destroy

   interface
      !> Fill an options struct with FINUFFT defaults
      !>
      !> @param[in,out] opts  A `finufft_opts`
      subroutine finufft_default_opts(opts) bind(c, name="finufft_default_opts_")
         implicit none(type, external)
         !> A `finufft_opts`
         type(*), intent(inout) :: opts
      end subroutine finufft_default_opts

      !> Create a guru plan
      !>
      !> @param[in]  ttype    Transform type, 1, 2 or 3
      !> @param[in]  dim      Spatial dimension
      !> @param[in]  n_modes  Modes per dimension, ignored by type 3
      !> @param[in]  iflag    Sign of the exponent
      !> @param[in]  ntrans   Number of simultaneous transforms
      !> @param[in]  eps      Requested relative tolerance
      !> @param[out] plan     Opaque plan handle
      !> @param[in]  opts     A `finufft_opts`
      !> @param[out] ier      FINUFFT status, zero on success
      subroutine finufft_makeplan(ttype, dim, n_modes, iflag, ntrans, eps, plan, opts, ier) &
         & bind(c, name="finufft_makeplan_")
         import :: c_int, c_int64_t, c_double
         implicit none(type, external)
         !> Transform type
         integer(c_int), intent(in) :: ttype
         !> Spatial dimension
         integer(c_int), intent(in) :: dim
         !> Modes per dimension
         integer(c_int64_t), intent(in) :: n_modes(*)
         !> Sign of the exponent
         integer(c_int), intent(in) :: iflag
         !> Number of simultaneous transforms
         integer(c_int), intent(in) :: ntrans
         !> Requested relative tolerance
         real(c_double), intent(in) :: eps
         !> Opaque plan handle
         integer(c_int64_t), intent(out) :: plan
         !> A `finufft_opts`
         type(*), intent(in) :: opts
         !> FINUFFT status
         integer(c_int), intent(out) :: ier
      end subroutine finufft_makeplan

      !> Set nonuniform points; the plan keeps pointers to all coordinate arrays
      !>
      !> @param[in]  plan  Opaque plan handle
      !> @param[in]  m     Number of nonuniform points
      !> @param[in]  xj    x coordinates of the nonuniform points
      !> @param[in]  yj    y coordinates of the nonuniform points
      !> @param[in]  zj    z coordinates of the nonuniform points
      !> @param[in]  nk    Number of type-3 targets
      !> @param[in]  s     x coordinates of the type-3 targets
      !> @param[in]  t     y coordinates of the type-3 targets
      !> @param[in]  u     z coordinates of the type-3 targets
      !> @param[out] ier   FINUFFT status, zero on success
      subroutine finufft_setpts(plan, m, xj, yj, zj, nk, s, t, u, ier) &
         & bind(c, name="finufft_setpts_")
         import :: c_int, c_int64_t, c_double
         implicit none(type, external)
         !> Opaque plan handle
         integer(c_int64_t), intent(in) :: plan
         !> Number of nonuniform points
         integer(c_int64_t), intent(in) :: m
         !> x coordinates of the nonuniform points
         real(c_double), intent(in) :: xj(*)
         !> y coordinates of the nonuniform points
         real(c_double), intent(in) :: yj(*)
         !> z coordinates of the nonuniform points
         real(c_double), intent(in) :: zj(*)
         !> Number of type-3 targets
         integer(c_int64_t), intent(in) :: nk
         !> x coordinates of the type-3 targets
         real(c_double), intent(in) :: s(*)
         !> y coordinates of the type-3 targets
         real(c_double), intent(in) :: t(*)
         !> z coordinates of the type-3 targets
         real(c_double), intent(in) :: u(*)
         !> FINUFFT status
         integer(c_int), intent(out) :: ier
      end subroutine finufft_setpts

      !> Run the planned transforms
      !>
      !> Input/output roles of `weights` and `result` depend on the plan type
      !>
      !> @param[in]     plan     Opaque plan handle
      !> @param[in,out] weights  Nonuniform strengths or values
      !> @param[in,out] result   Mode coefficients or type-3 target values
      !> @param[out]    ier      FINUFFT status, zero on success
      subroutine finufft_execute(plan, weights, result, ier) bind(c, name="finufft_execute_")
         import :: c_int, c_int64_t, c_double_complex
         implicit none(type, external)
         !> Opaque plan handle
         integer(c_int64_t), intent(in) :: plan
         !> Nonuniform strengths or values
         complex(c_double_complex), intent(inout) :: weights(*)
         !> Mode coefficients or type-3 target values
         complex(c_double_complex), intent(inout) :: result(*)
         !> FINUFFT status
         integer(c_int), intent(out) :: ier
      end subroutine finufft_execute

      !> Free a plan
      !>
      !> @param[in,out] plan  Opaque plan handle
      !> @param[out]    ier   FINUFFT status, zero on success
      subroutine finufft_destroy(plan, ier) bind(c, name="finufft_destroy_")
         import :: c_int, c_int64_t
         implicit none(type, external)
         !> Opaque plan handle
         integer(c_int64_t), intent(inout) :: plan
         !> FINUFFT status
         integer(c_int), intent(out) :: ier
      end subroutine finufft_destroy
   end interface

end module moist_math_fft_finufft
