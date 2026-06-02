!> Proxy module to reexport high-level basic linear algebra subprogram wrappers
module moist_math_blas
   use moist_math_blas_level1, only: dot => wrap_dot, nrm2 => wrap_nrm2, &
      & scal => wrap_scal, copy => wrap_copy, axpy => wrap_axpy
   use moist_math_blas_level2, only: gemv => wrap_gemv, symv => wrap_symv
   use moist_math_blas_level3, only: gemm => wrap_gemm, syrk => wrap_syrk
   implicit none
   public
end module moist_math_blas
