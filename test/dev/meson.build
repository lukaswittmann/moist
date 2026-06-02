dev_tests = [
  'cavity_drop_born_fit',
  'cavity_drop_timings',
  'cavity_drop_robustness',
  'cavity_drop_integration',
  'cavity_drop_convergence',
  'cavity_drop_deflation_comparison',
]

if use_rism
  dev_tests += []
endif

dev_test_srcs = files('main.F90')
foreach t : dev_tests
  dev_test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

dev_tester = executable(
  'dev_tester',
  sources: dev_test_srcs,
  dependencies: [moist_dep, mstore_dep],
  link_language: 'fortran',
)
