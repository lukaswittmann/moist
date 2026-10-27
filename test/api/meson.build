test(
  'api-test',
  executable(
    'api-tester',
    sources: files('example.c'),
    dependencies: moist_dep,
  ),
)

test(
  'api-cpp',
  executable(
    'api-cpp-tester',
    sources: files('cpp.cpp'),
    dependencies: moist_dep,
    override_options: [cxx.get_argument_syntax() == 'msvc' ? 'cpp_std=c++14' : 'cpp_std=c++11'],
    cpp_args: cxx.get_supported_arguments('-pedantic-errors'),
  ),
)
