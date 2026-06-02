# Standalone build for Python-API of moist, requires a moist library
# in the PKG_CONFIG_PATH environment variable to work.
project(
  'moist',
  'c',
  version: '0.6.0a1',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
install = true

moist_dep = dependency(
  'moist',
  version: '>=@0@'.format(meson.project_version()),
  fallback: ['moist', 'moist_dep'],
  default_options: [
    'default_library=static',
    'api=true',
    'python=false',
  ],
)
moist_header = files('include'/'moist.h')

subdir('moist')
