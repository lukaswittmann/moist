subdir('data')
subdir('radii')
subdir('model')
subdir('math')
subdir('cavity')
subdir('utils')
subdir('output')

srcs += files(
  'context.f90',
  'math.f90',
  'cavity.f90',
  'data.f90',
  'radii.f90',
  'model.f90',
  'output.f90',
  'version.f90',
  'type.f90',
  'utils.f90',
)

# Build-info module
srcs += vcs_tag(
  command: ['git', 'rev-parse', '--short', 'HEAD'],
  fallback: 'unknown',
  input: 'build_info.f90.in',
  output: 'build_info.f90',
)

if get_option('api') or get_option('python')
  srcs += files('api.f90')
endif
