!> FFTW module using Fortran 2003 interfaces.
module moist_math_fftw3
   use, intrinsic :: iso_c_binding
#ifdef WITH_FFTW
#  ifdef MPI
   include 'fftw3-mpi.f03'
#  else
   include 'fftw3.f03'
#  endif
#endif
end module moist_math_fftw3
