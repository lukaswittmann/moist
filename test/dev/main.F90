

!> Driver for long-running diagnostic (dev) tests
program dev_tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected, get_argument
   use test_cavity_drop_born_fit, only : collect_cavity_drop_born_fit
   use test_cavity_drop_timings, only : collect_cavity_drop_timings
   use test_cavity_drop_robustness, only : collect_cavity_drop_robustness
   use test_cavity_drop_integration, only : collect_cavity_drop_integration
   use test_cavity_drop_convergence, only : collect_cavity_drop_convergence
   use test_cavity_drop_deflation_comparison, only : collect_cavity_drop_deflation_comparison
#ifdef WITH_RISM
   ! use test_rism_1d_vv, only : collect_rism_1d_vv
#endif

implicit none

   integer :: stat, is
   character(len=:), allocatable :: suite_name, test_name
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
      & new_testsuite("cavity_drop_born_fit", collect_cavity_drop_born_fit), &
      & new_testsuite("cavity_drop_timings", collect_cavity_drop_timings), &
      & new_testsuite("cavity_drop_robustness", collect_cavity_drop_robustness), &
      & new_testsuite("cavity_drop_integration", collect_cavity_drop_integration), &
      & new_testsuite("cavity_drop_deflation_comparison", collect_cavity_drop_deflation_comparison), &
      & new_testsuite("cavity_drop_convergence", collect_cavity_drop_convergence) &
#ifdef WITH_RISM
      ! , new_testsuite("rism_1d_vv", collect_rism_1d_vv) &
#endif
      & ]

   call get_argument(1, suite_name)
   call get_argument(2, test_name)

   if (allocated(suite_name)) then
      is = select_suite(testsuites, suite_name)
      if (is > 0 .and. is <= size(testsuites)) then
         if (allocated(test_name)) then
            write(error_unit, fmt) "Suite:", testsuites(is)%name
            call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
            if (stat < 0) then
               error stop 1
            end if
         else
            write(error_unit, fmt) "Testing:", testsuites(is)%name
            call run_testsuite(testsuites(is)%collect, error_unit, stat)
         end if
      else
         write(error_unit, fmt) "Available testsuites"
         do is = 1, size(testsuites)
            write(error_unit, fmt) "-", testsuites(is)%name
         end do
         error stop 1
      end if
   else
      do is = 1, size(testsuites)
         write(error_unit, fmt) "Testing:", testsuites(is)%name
         call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end do
   end if

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
   end if


end program dev_tester
