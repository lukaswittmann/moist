project(
  'moist',
  'fortran',
  version: '0.6.0-alpha.1',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.57',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
    'c_std=c11',
  ],
)
lib_version = meson.project_version().split('-')[0]
install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: get_option('api') or get_option('python'), native: false)

# General configuration information
inc_dirs = []
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
moist_lib = library(
  meson.project_name(),
  sources: srcs,
  version: lib_version,
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
moist_inc = [moist_lib.private_dir_include(), include_directories('include')]
moist_dep = declare_dependency(
  link_with: moist_lib,
  include_directories: moist_inc,
  dependencies: lib_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)

# Add applications
subdir('app')

# Package the license files
moist_lic = files(
  'COPYING',
  'COPYING.LESSER',
  'THIRD_PARTY_LICENSES.md',
)

moist_header = files(
  'include/moist.h',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    moist_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  if get_option('api')
    install_headers(
      moist_header,
    )
  endif

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    moist_lib,
    name: 'moist',
    description: 'Modular and Open-source Implicit Solvation Toolkit',
    version: meson.project_version(),
    subdirs: ['', module_id],
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found() and get_option('docs')
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man/moist.adoc'),
        output: '@BASENAME@',
      )
    )
  endif
endif

# add the testsuite
subdir('test')

if get_option('python')
  subdir('python'/'moist')
endif
