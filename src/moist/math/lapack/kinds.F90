module moist_math_lapack_kinds
   use iso_fortran_env, only: int32, int64
   implicit none
   private

#if MOIST_ILP64
   integer, parameter, public :: lapack_ik = int64
#else
   integer, parameter, public :: lapack_ik = int32
#endif

end module moist_math_lapack_kinds
