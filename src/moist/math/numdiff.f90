!> Proxy module to reexport numerical differentiation (Jacobian computation) functionality
!>
!> This module provides a simple interface to the NumDiff library for computing
!> Jacobian matrices using finite differences.
!>
!> ## Basic Usage
!>
!> 1. Declare a `numdiff_type` object
!> 2. Initialize it with your problem dimensions and function
!> 3. Call `compute_jacobian` or `compute_jacobian_dense` to get the Jacobian
!>
!> ## Example
!>
!> ```fortran
!> use moist_math_numdiff
!> type(numdiff_type) :: jac_calculator
!> real(numdiff_wp) :: x(n), jac_dense(m,n)
!>
!> ! Initialize with n variables, m functions
!> call jac_calculator%initialize(n, m, xlow, xhigh, problem_func, &
!>                                sparsity_mode=numdiff_sparsity_auto, &
!>                                jacobian_method=numdiff_method_5point_central, &
!>                                dpert=1.0e-8_numdiff_wp)
!>
!> ! Compute Jacobian at point x
!> call jac_calculator%compute_jacobian_dense(x, jac_dense)
!>
!> ! Clean up
!> call jac_calculator%destroy()
!> ```
!>
!> See the numerical_differentiation_module documentation for advanced options.
module moist_math_numdiff
   ! Re-export core types for Jacobian computation
   use numerical_differentiation_module, only: &
      & numdiff_type, &
      & finite_diff_method, &
      & sparsity_pattern, &
      & get_finite_diff_formula, &
      & get_all_methods_in_class

   ! Re-export precision
   use mctc_env_accuracy, only: numdiff_wp => wp

   implicit none
   public

   !> Common sparsity modes for ease of use
   integer, parameter, public :: numdiff_sparsity_dense = 1
      !! Assume dense Jacobian (all elements are non-zero)
   integer, parameter, public :: numdiff_sparsity_three_point = 2
      !! Three-point simple method to detect sparsity
   integer, parameter, public :: numdiff_sparsity_user = 3
      !! User will specify sparsity pattern manually
   integer, parameter, public :: numdiff_sparsity_auto = 4
      !! Automatic sparsity detection via multiple point evaluation

   !> Common perturbation modes for ease of use
   integer, parameter, public :: numdiff_perturb_absolute = 1
      !! Perturbation is dx=dpert (absolute)
   integer, parameter, public :: numdiff_perturb_relative = 2
      !! Perturbation is dx=dpert*x (relative to x)
   integer, parameter, public :: numdiff_perturb_mixed = 3
      !! Perturbation is dx=dpert*(1+x) (mixed mode)

   !> Common finite difference method IDs
   integer, parameter, public :: numdiff_method_2point_forward = 1
      !! 2-point forward difference: (f(x+h)-f(x))/h
   integer, parameter, public :: numdiff_method_2point_backward = 2
      !! 2-point backward difference: (f(x)-f(x-h))/h
   integer, parameter, public :: numdiff_method_3point_central = 3
      !! 3-point central difference: (f(x+h)-f(x-h))/(2h)
   integer, parameter, public :: numdiff_method_5point_central = 10
      !! 5-point central difference (higher accuracy)
   integer, parameter, public :: numdiff_method_7point_central = 21
      !! 7-point central difference (even higher accuracy)

contains

end module moist_math_numdiff
