# Test suites are grouped by module into subfolders that mirror src/moist/
#
# Accumulators the subdirs append to:
#   tests         suite ids, passed to the tester as its argument
#   test_srcs     sources of the tester executable (suites plus shared helpers)
#   test_timeouts suite id -> timeout, for suites needing more than the default
tests = []
test_timeouts = {}

# main driver plus the cross-cutting helper module (not registered as a suite).
test_srcs = files(
  'main.F90',
  'test_helpers.f90',
)

subdir('test_math')
subdir('test_cavity')
subdir('test_model')
subdir('test_utils')

# Suites kept flat at the test/unit root, each in its own test_<id>.f90. Two kinds live
# here: singletons whose suite id cannot take a folder prefix (a test_utils/ folder would
# force every leaf to be 'utils_<leaf>', leaving a plain 'utils' suite unreachable), and
# smoke tests for optional vendored dependencies, which exercise upstream modules
# directly and so mirror nothing under src/moist/.
root_tests = ['api', 'utils', 'radii', 'data']
if use_hdf5
  root_tests += 'utils_hdf5'
endif

foreach leaf : root_tests
  tests += leaf
  test_srcs += files('test_' + leaf + '.f90')
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [moist_dep, mstore_dep, testdrive_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(
    t,
    tester,
    args: t,
    env: test_env,
    is_parallel : true,
    timeout: test_timeouts.get(t, 360),
  )
endforeach
