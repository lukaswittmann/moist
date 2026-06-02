tests = [
  'utils',
  'radii',
  'math_linalg',
  'math_adjacency_list',
  'math_cell_grid',
  'math_solvers',
  'math_sorters',
  'math_trig',
  'math_grid',
  'cavity_iswig',
  'cavity_numsa',
  'component_pcm_cpcm',
  'cavity_drop_primitives',
  'cavity_drop_cfc_kernel',
  'cavity_drop_lsf',
  'cavity_drop_gradients',
  'cavity_drop_cpcm',
  'cavity_drop_integration_ref',
]

if use_rism
  tests += []
endif

if use_hdf5
  tests += ['utils_hdf5']
endif

test_srcs = files(
  'main.F90',
  'test_helpers.f90',
)

foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [moist_dep, mstore_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(
    t,
    tester,
    args: t,
    env: test_env,
    is_parallel : true,
    timeout: 360,
  )
endforeach
