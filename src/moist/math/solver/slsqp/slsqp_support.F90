!*******************************************************************************
!> license: BSD
!
!  Numerical constants shared across the SLSQP sources. The level-1 BLAS
!  routines are imported directly from moist_math_blas_legacy by the callers.

module slsqp_support

   use mctc_env_accuracy, only: wp
   implicit none

   private

   real(wp), parameter, public :: epmach = epsilon(1.0_wp)
   real(wp), parameter, public :: zero = 0.0_wp
   real(wp), parameter, public :: one = 1.0_wp
   real(wp), parameter, public :: two = 2.0_wp
   real(wp), parameter, public :: four = 4.0_wp
   real(wp), parameter, public :: ten = 10.0_wp
   real(wp), parameter, public :: hun = 100.0_wp

!*******************************************************************************
end module slsqp_support
!*******************************************************************************
