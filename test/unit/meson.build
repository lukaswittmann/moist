# One suite per file, flat: suite id <id> lives in test_<id>.f90 and defines
# module test_<id>. Keeping the three names in lockstep is what lets this list
# be the single source of truth -- the file name is derived, never spelled out.
#
#   tests             suite ids, passed to the tester as its argument
#   test_timeouts     suite id -> timeout, for suites needing more than the default
#   test_selected     suite id -> test names, for a suite whose tests must each be
#                      invoked as `tester <suite> <test>` (test-drive's run_selected, no
#                      enclosing OpenMP team) rather than `tester <suite>`
#                      (run_testsuite, whose tests run inside test-drive's own
#                      `!$omp parallel do`; see CONTRIBUTING.md:47)
tests = [
  'api',
  'parameters',
  'utils_timer',
  'utils_context',
  'utils_mem',
  'utils_prettylistprint',
  'channels',
  'data',
  'radii',
  'math_linalg',
  'math_smoothing_kernels',
  'math_adjacency_list',
  'math_cell_grid',
  'math_solvers',
  'math_sorters',
  'math_trig',
  'math_grid',
  'math_quadrature',
  'math_grid_s2',
  'math_fft_dst4',
  'math_fft_3d',
  'finufft',
  'math_grid_nufft',
  'math_grid_3d',
  'math_grid_3d_threaded',
  'moz_1d',
  'moz_3d',
  'cavity_drop_primitives',
  'cavity_drop_cfc',
  'cavity_drop_lsf',
  'cavity_drop_lsf_golden',
  'cavity_drop_isodensity',
  'cavity_drop_gradient',
  'cavity_drop_nuclear_adjoint',
  'cavity_drop_integration',
  'cavity_drop_filter',
  'cavity_drop_cpcm',
  'cavity_iswig',
  'cavity_numsa',
  'cavity_marchingcubes',
  'cavity_surface_adjoint',
  'model_component_pcm_amat',
  'model_component_pcm_amat_kernel',
  'model_component_pcm_amat_assembly',
  'model_component_pcm_amat_adjoint',
  'model_component_pcm_electrostatics',
  'model_component_pcm_cpcm',
  'model_component_gostshyp',
  'model_component_pv',
  'model_continuum',
  'model_coupling',
]
test_timeouts = {}

# threaded tests (single testdrive thread)
test_selected = {
  'math_grid_3d_threaded': [
    'batched_fft_uses_threaded_backend',
    'batched_ifft_uses_threaded_backend',
    'molecular_nufft_type12_threaded',
    'molecular_nufft_type3_threaded',
  ],
}

# Driver plus the shared fixture modules, which register no suite of their own.
test_srcs = files(
  'main.F90',
  'test_helpers.f90',
  'test_model_component_helper.f90',
)

foreach t : tests
  test_srcs += files('test_' + t + '.f90')
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [moist_dep, mstore_dep, testdrive_dep],
  link_language: 'fortran',
)

foreach t : tests
  runs = {t: [t]}
  if test_selected.has_key(t)
    runs = {}
    foreach sel : test_selected[t]
      runs += {t + '_' + sel: [t, sel]}
    endforeach
  endif
  foreach name, test_args : runs
    test(
      name,
      tester,
      args: test_args,
      env: test_env,
      is_parallel : true,
      timeout: test_timeouts.get(t, 360),
    )
  endforeach
endforeach
