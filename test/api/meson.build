test(
  'api-test',
  executable(
    'api-tester',
    sources: files('example.c'),
    dependencies: moist_dep,
  ),
)

add_languages('cpp', required: true, native: false)
test(
  'api-cpp',
  executable(
    'api-cpp-tester',
    sources: files('cpp.cpp'),
    dependencies: moist_dep,
    override_options: ['cpp_std=c++11'],
    cpp_args: ['-pedantic-errors'],
  ),
)
