

!> Driver for unit testing
program tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected, get_argument
   use test_utils, only : collect_utils
   use test_radii, only : collect_radii
   use test_math_linalg, only : collect_math_linalg
   use test_math_adjacency_list, only : collect_math_adjacency_list
   use test_math_cell_grid, only : collect_math_cell_grid
   use test_math_sorters, only : collect_math_sorters
   use test_math_trig, only : collect_math_trig
   use test_math_grid, only : collect_math_grid
   use test_cavity_iswig, only : collect_cavity_iswig
   use test_cavity_drop_primitives, only : collect_cavity_drop_primitives
   use test_cavity_drop_cfc_kernel, only : collect_cavity_drop_cfc_kernel
   use test_cavity_drop_lsf, only : collect_cavity_drop_lsf
   use test_cavity_drop_gradients, only : collect_cavity_drop_gradients
   use test_cavity_drop_cpcm, only : collect_cavity_drop_cpcm
   use test_cavity_numsa, only : collect_cavity_numsa
   use test_math_solvers, only : collect_math_solvers
#ifdef WITH_RISM
   ! use test_rism_1d, only : collect_rism_1d
   ! use test_rism_thermo, only : collect_rism_thermo
#endif
#ifdef WITH_HDF5
   use test_utils_hdf5, only : collect_utils_hdf5
#endif
   use test_cavity_drop_integration_ref, only : collect_cavity_drop_integration_ref
   use test_component_pcm_cpcm, only : collect_component_pcm_cpcm

implicit none

   integer :: stat, is
   character(len=:), allocatable :: suite_name, test_name
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
#ifdef WITH_HDF5
      & new_testsuite("utils_hdf5", collect_utils_hdf5), &
#endif
#ifdef WITH_RISM
      ! & new_testsuite("rism_1d", collect_rism_1d), &
      ! & new_testsuite("rism_thermo", collect_rism_thermo), &
#endif
      & new_testsuite("utils", collect_utils), &
      & new_testsuite("radii", collect_radii), &
      & new_testsuite("math_linalg", collect_math_linalg), &
      & new_testsuite("math_adjacency_list", collect_math_adjacency_list), &
      & new_testsuite("math_cell_grid", collect_math_cell_grid), &
      & new_testsuite("math_solvers", collect_math_solvers), &
      & new_testsuite("math_sorters", collect_math_sorters), &
      & new_testsuite("math_trig", collect_math_trig), &
      & new_testsuite("math_grid", collect_math_grid), &

      & new_testsuite("cavity_drop_primitives", collect_cavity_drop_primitives), &
      & new_testsuite("cavity_drop_cfc_kernel", collect_cavity_drop_cfc_kernel), &
      & new_testsuite("cavity_drop_lsf", collect_cavity_drop_lsf), &
      & new_testsuite("cavity_drop_gradients", collect_cavity_drop_gradients), &
      & new_testsuite("cavity_drop_integration_ref", collect_cavity_drop_integration_ref), &
      & new_testsuite("cavity_drop_cpcm", collect_cavity_drop_cpcm), &
      & new_testsuite("cavity_iswig", collect_cavity_iswig), &
      & new_testsuite("cavity_numsa", collect_cavity_numsa), &
      & new_testsuite("component_pcm_cpcm", collect_component_pcm_cpcm) &
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


end program tester
