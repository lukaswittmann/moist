if get_option('api')
  subdir('api')
endif

# Shared across unit and dev subdirs
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

testdrive_dep = dependency(
  'test-drive',
  version: '>=0.5.0',
  fallback: ['test-drive', 'testdrive_dep'],
  default_options: ['default_library=static'],
)

test_env = environment()
test_env.set('MOIST_SOURCE_ROOT', meson.project_source_root())

# Serial BLAS: suites run their tests in an OpenMP parallel; meson runs suites
# concurrently, so threaded BLAS (MKL on intel25, pthreaded OpenBLAS on the
# conda, ilp64 and MSYS2 legs) oversubscribes
test_env.set('MKL_NUM_THREADS', '1')
test_env.set('OPENBLAS_NUM_THREADS', '1')

subdir('unit')

if get_option('dev_tests')
  subdir('dev')
endif
