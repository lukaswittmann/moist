module test_utils
   use mctc_env, only : wp
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   use moist_data_solvents, only : solvation_system_parameters, &
      & new_solvation_system_parameters, get_solvent_id, max_solvents
   use moist_utils, only : lowercase, is_exceptional
   use mctc_env_error, only : moist_error_type => error_type
   use, intrinsic :: ieee_arithmetic
   implicit none
   private
   public :: collect_utils
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))


contains

!> Collect all exported unit tests
subroutine collect_utils(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      & new_unittest("SumSolventProperties", test_sum_solvent_properties), &
      & new_unittest("SolventIDFinder", test_solvent_id_finder) &

      & ]
end subroutine collect_utils

subroutine test_sum_solvent_properties(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: i
   real(wp) :: sum_eps, sum_refr, sum_A, sum_B, sum_g, sum_rho
   real(wp), dimension(max_solvents) :: eps, refr, A, B, g, rho
   integer :: id_list(max_solvents)
   character(len=64) :: name_list(max_solvents)
   character(len=64) :: alias_list(10,max_solvents)

   real(wp), parameter :: sum_eps_ref  =   2079.0211_wp
   real(wp), parameter :: sum_refr_ref =    260.1011_wp
   real(wp), parameter :: sum_A_ref    =     16.4500_wp
   real(wp), parameter :: sum_B_ref    =     55.3700_wp
   real(wp), parameter :: sum_g_ref    =      7.4137_wp
   real(wp), parameter :: sum_rho_ref  = 181468.0000_wp

   sum_eps = 0.0_wp
   sum_refr = 0.0_wp
   sum_A = 0.0_wp
   sum_B = 0.0_wp
   sum_g = 0.0_wp
   sum_rho = 0.0_wp

   include "../src/moist/data/solvents.inc"
   do i = 1, max_solvents
      sum_eps  = sum_eps  + eps(i)
      sum_refr = sum_refr + refr(i)
      sum_A    = sum_A    + A(i)
      sum_B    = sum_B    + B(i)
      sum_g    = sum_g    + g(i) * 0.001_wp
      sum_rho  = sum_rho  + rho(i)
   end do

   call check(error, sum_eps , sum_eps_ref ,thr=thr, &
      & more="Sum of epsilon does not match")
   call check(error, sum_refr, sum_refr_ref,thr=thr, &
      & more="Sum of refractive index does not match")
   call check(error, sum_A   , sum_A_ref   ,thr=thr, &
      & more="Sum of alpha does not match")
   call check(error, sum_B   , sum_B_ref   ,thr=thr, &
      & more="Sum of beta does not match")
   call check(error, sum_g   , sum_g_ref   ,thr=thr, &
      & more="Sum of surface tension does not match")
   call check(error, sum_rho , sum_rho_ref ,thr=thr, &
      & more="Sum of mass density does not match")

end subroutine test_sum_solvent_properties

subroutine test_solvent_id_finder(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: id
   type(moist_error_type), allocatable :: solvent_error
   call get_solvent_id('wAtEr', id, solvent_error)
   if (allocated(solvent_error)) then
      call test_failed(error, solvent_error%message)
      return
   end if
   call check(error, id, 175, more="Solvent ID lookup mismatch")
   call get_solvent_id('methyl chloroform', id, solvent_error)
   if (allocated(solvent_error)) then
      call test_failed(error, solvent_error%message)
      return
   end if
   call check(error, id, 1, more="Solvent ID lookup mismatch")
end subroutine test_solvent_id_finder


end module test_utils
