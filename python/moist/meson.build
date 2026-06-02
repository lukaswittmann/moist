cc = meson.get_compiler('c')

pymod = import('python')
python = pymod.find_installation(
  get_option('python_version'),
  modules: [
    'cffi',
    'setuptools',
  ],
)
python_dep = python.dependency(required: true)

moist_pp = configure_file(
  command: [
    cc,
    '-I@0@'.format(
      moist_dep.get_variable(
        pkgconfig: 'includedir',
        cmake: 'moist_INCLUDE_DIRS',
        internal: 'includedir',
      ).split().get(0)
    ),
    '-Dmoist_CFFI',
    '-E',
    '@INPUT@',
  ],
  input: moist_header,
  output: '_libmoist.h',
  capture: true,
)

moist_cffi_srcs = configure_file(
  command: [python, files('..'/'ffibuilder.py'), '@INPUT@', '@BASENAME@'],
  input: moist_pp,
  output: '@BASENAME@.c',
)

moist_pyext = python.extension_module(
  '_libmoist',
  link_whole: static_library(
    '_libmoist',
    moist_cffi_srcs,
    dependencies: [moist_dep, python_dep],
  ),
  dependencies: [moist_dep, python_dep],
  install: install,
  subdir: 'moist',
)

pysrcs = files(
  '__init__.py',
  'interface.py',
  'library.py',
  'test_interface.py',
  'test_library.py',
)

if install
  python.install_sources(
    pysrcs,
    subdir: 'moist',
    pure: false,
  )
endif
