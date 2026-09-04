
!> Driver for unit testing
program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected, get_argument
   use test_api, only: collect_api
   use test_utils, only: collect_utils
   use test_utils_timer, only: collect_utils_timer
   use test_utils_context, only: collect_utils_context
   use test_utils_mem, only: collect_utils_mem
   use test_utils_prettylistprint, only: collect_utils_prettylistprint
   use test_radii, only: collect_radii
   use test_data, only: collect_data
   use test_math_linalg, only: collect_math_linalg
   use test_math_smoothing_kernels, only: collect_math_smoothing_kernels
   use test_math_adjacency_list, only: collect_math_adjacency_list
   use test_math_cell_grid, only: collect_math_cell_grid
   use test_math_sorters, only: collect_math_sorters
   use test_math_trig, only: collect_math_trig
   use test_math_grid, only: collect_math_grid
   use test_cavity_iswig, only: collect_cavity_iswig
   use test_cavity_drop_primitives, only: collect_cavity_drop_primitives
   use test_cavity_drop_kkt, only: collect_cavity_drop_kkt
   use test_cavity_drop_field_tangent, only: collect_cavity_drop_field_tangent
   use test_cavity_drop_iswig_scatter, only: collect_cavity_drop_iswig_scatter
   use test_cavity_drop_hessian_fixed, only: collect_cavity_drop_hessian_fixed
   use test_cavity_drop_tangent_forward, only: collect_cavity_drop_tangent_forward
   use test_cavity_drop_hessian_response, only: collect_cavity_drop_hessian_response
   use test_cavity_drop_hessian_e2e, only: collect_cavity_drop_hessian_e2e
   use test_cavity_drop_weights_tangent, only: collect_cavity_drop_weights_tangent
   use test_cavity_drop_cfc, only: collect_cavity_drop_cfc
   use test_cavity_drop_lsf, only: collect_cavity_drop_lsf
   use test_cavity_drop_lsf_golden, only: collect_cavity_drop_lsf_golden
   use test_cavity_drop_isodensity, only: collect_cavity_drop_isodensity
   use test_cavity_drop_gradient, only: collect_cavity_drop_gradient
   use test_cavity_drop_nuclear_adjoint, only: collect_cavity_drop_nuclear_adjoint
   use test_cavity_drop_cpcm, only: collect_cavity_drop_cpcm
   use test_cavity_numsa, only: collect_cavity_numsa
   use test_cavity_marchingcubes, only: collect_cavity_marchingcubes
   use test_math_solvers, only: collect_math_solvers
#ifdef WITH_HDF5
   use test_utils_hdf5, only: collect_utils_hdf5
#endif
   use test_cavity_drop_integration, only: collect_cavity_drop_integration
   use test_cavity_drop_filter, only: collect_cavity_drop_filter
   use test_model_component_pcm_amat, only: collect_model_component_pcm_amat
   use test_model_component_pcm_amat_kernel, only: collect_model_component_pcm_amat_kernel
   use test_model_component_pcm_amat_assembly, only: collect_model_component_pcm_amat_assembly
   use test_model_component_pcm_amat_adjoint, only: collect_model_component_pcm_amat_adjoint
   use test_model_component_pcm_electrostatics, only: collect_model_component_pcm_electrostatics
   use test_model_component_pcm_cpcm, only: collect_model_component_pcm_cpcm
   use test_model_component_gostshyp, only: collect_model_component_gostshyp
   use test_model_component_pv, only: collect_model_component_pv
   use test_model_general, only: collect_model_general

   implicit none(type, external)

   integer :: stat, is
   character(len=:), allocatable :: suite_name, test_name
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
      & new_testsuite("api", collect_api), &
      & new_testsuite("utils", collect_utils), &
      & new_testsuite("utils_timer", collect_utils_timer), &
      & new_testsuite("utils_context", collect_utils_context), &
      & new_testsuite("utils_mem", collect_utils_mem), &
      & new_testsuite("utils_prettylistprint", collect_utils_prettylistprint), &
      & new_testsuite("data", collect_data), &
      & new_testsuite("radii", collect_radii), &
      & new_testsuite("math_linalg", collect_math_linalg), &
      & new_testsuite("math_smoothing_kernels", collect_math_smoothing_kernels), &
      & new_testsuite("math_adjacency_list", collect_math_adjacency_list), &
      & new_testsuite("math_cell_grid", collect_math_cell_grid), &
      & new_testsuite("math_solvers", collect_math_solvers), &
      & new_testsuite("math_sorters", collect_math_sorters), &
      & new_testsuite("math_trig", collect_math_trig), &
      & new_testsuite("math_grid", collect_math_grid), &
      & new_testsuite("cavity_drop_primitives", collect_cavity_drop_primitives), &
      & new_testsuite("cavity_drop_kkt", collect_cavity_drop_kkt), &
      & new_testsuite("cavity_drop_field_tangent", collect_cavity_drop_field_tangent), &
      & new_testsuite("cavity_drop_iswig_scatter", collect_cavity_drop_iswig_scatter), &
      & new_testsuite("cavity_drop_hessian_fixed", collect_cavity_drop_hessian_fixed), &
      & new_testsuite("cavity_drop_tangent_forward", collect_cavity_drop_tangent_forward), &
      & new_testsuite("cavity_drop_hessian_response", collect_cavity_drop_hessian_response), &
      & new_testsuite("cavity_drop_hessian_e2e", collect_cavity_drop_hessian_e2e), &
      & new_testsuite("cavity_drop_weights_tangent", collect_cavity_drop_weights_tangent), &
      & new_testsuite("cavity_drop_cfc", collect_cavity_drop_cfc), &
      & new_testsuite("cavity_drop_lsf", collect_cavity_drop_lsf), &
      & new_testsuite("cavity_drop_lsf_golden", collect_cavity_drop_lsf_golden), &
      & new_testsuite("cavity_drop_isodensity", collect_cavity_drop_isodensity), &
      & new_testsuite("cavity_drop_gradient", collect_cavity_drop_gradient), &
      & new_testsuite("cavity_drop_nuclear_adjoint", collect_cavity_drop_nuclear_adjoint), &
      & new_testsuite("cavity_drop_integration", collect_cavity_drop_integration), &
      & new_testsuite("cavity_drop_filter", collect_cavity_drop_filter), &
      & new_testsuite("cavity_drop_cpcm", collect_cavity_drop_cpcm), &
      & new_testsuite("cavity_iswig", collect_cavity_iswig), &
      & new_testsuite("cavity_numsa", collect_cavity_numsa), &
      & new_testsuite("cavity_marchingcubes", collect_cavity_marchingcubes), &
      & new_testsuite("model_component_pcm_amat", collect_model_component_pcm_amat), &
      & new_testsuite("model_component_pcm_amat_kernel", collect_model_component_pcm_amat_kernel), &
      & new_testsuite("model_component_pcm_amat_assembly", &
         collect_model_component_pcm_amat_assembly), &
      & new_testsuite("model_component_pcm_amat_adjoint", &
         collect_model_component_pcm_amat_adjoint), &
      & new_testsuite("model_component_pcm_electrostatics", &
         collect_model_component_pcm_electrostatics), &
      & new_testsuite("model_component_pcm_cpcm", collect_model_component_pcm_cpcm), &
      & new_testsuite("model_component_gostshyp", collect_model_component_gostshyp), &
      & new_testsuite("model_component_pv", collect_model_component_pv), &
      & new_testsuite("model_general", collect_model_general) &
      & ]

#ifdef WITH_HDF5
   testsuites = [testsuites, new_testsuite("utils_hdf5", collect_utils_hdf5)]
#endif

   call get_argument(1, suite_name)
   call get_argument(2, test_name)

   if (allocated(suite_name)) then
      is = select_suite(testsuites, suite_name)
      if (is > 0 .and. is <= size(testsuites)) then
         if (allocated(test_name)) then
            write (error_unit, fmt) "Suite:", trim(testsuites(is)%name)
            call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
            if (stat < 0) then
               error stop 1
            end if
         else
            write (error_unit, fmt) "Testing:", trim(testsuites(is)%name)
            call run_testsuite(testsuites(is)%collect, error_unit, stat)
         end if
      else
         write (error_unit, fmt) "Available testsuites"
         do is = 1, size(testsuites)
            write (error_unit, fmt) "-", trim(testsuites(is)%name)
         end do
         error stop 1
      end if
   else
      do is = 1, size(testsuites)
         write (error_unit, fmt) "Testing:", trim(testsuites(is)%name)
         call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end do
   end if

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if

end program tester
