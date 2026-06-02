if get_option('api')
  subdir('api')
endif

# Shared across unit and dev subdirs.
mstore_dep = dependency(
  'mstore',
  version: '>=0.3.0',
  fallback: ['mstore', 'mstore_dep'],
  default_options: ['default_library=static'],
  required: not meson.is_subproject(),
)
if not mstore_dep.found()
  subdir_done()
endif

test_env = environment()
test_env.set('MOIST_SOURCE_ROOT', meson.project_source_root())

# Increase small ifx default OMP_STACKSIZE to avoid SIGSEGV for some tests
test_env.set('OMP_STACKSIZE', '128M')
# Serial test for intel25
test_env.set('MKL_NUM_THREADS', '1')

subdir('unit')

if get_option('dev_tests')
  subdir('dev')
endif
