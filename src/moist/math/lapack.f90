!> Proxy module to reexport high-level linear algebra package wrappers
module moist_math_lapack
   use moist_math_lapack_getrf, only: getrf => wrap_getrf
   use moist_math_lapack_getri, only: getri => wrap_getri
   use moist_math_lapack_getrs, only: getrs => wrap_getrs
   use moist_math_lapack_potrf, only: potrf => wrap_potrf
   use moist_math_lapack_potrs, only: potrs => wrap_potrs
   use moist_math_lapack_syev, only: lapack_syev
   use moist_math_lapack_gesv, only: lapack_gesv
   use moist_math_lapack_gesvd, only: lapack_gesvd
   use moist_math_lapack_gels, only: lapack_gels
   implicit none
   public
end module moist_math_lapack
