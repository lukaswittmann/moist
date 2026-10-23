!> Proxy module to reexport high-level basic linear algebra subprogram wrappers
module moist_math_blas
   use moist_math_blas_level1, only: dot => wrap_dot, nrm2 => wrap_nrm2, &
      & scal => wrap_scal, copy => wrap_copy, axpy => wrap_axpy
   use moist_math_blas_level2, only: gemv => wrap_gemv, symv => wrap_symv, gemv_rows => wrap_gemv_rows
   use moist_math_blas_level3, only: gemm => wrap_gemm, syrk => wrap_syrk, &
      & gemm_rows_tn => wrap_gemm_rows_tn, gemm_rows_nn => wrap_gemm_rows_nn
   implicit none(type, external)
   public
end module moist_math_blas
