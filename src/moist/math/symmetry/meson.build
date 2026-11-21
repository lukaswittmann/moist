# Point-group symmetry backend: moist's C shim for the vendored libmsym
#
# - static library with implicit_include_directories: false keeps the source
#   root off the C include path, like the FFT backend
# - lib_deps gets link_with only, so Fortran compiles never see a C include
#   directory
# - libmsym's headers are system headers: their warnings stay out of the shim
msym_inc = include_directories('libmsym', is_system: true)
msym_c_args = meson.get_compiler('c').get_supported_arguments(
  # Vendored code is kept byte-for-byte; silence its own warnings
  '-w',
)

libmsym_lib = static_library(
  'moist_libmsym',
  sources: files(
    'libmsym' / 'basis_function.c',
    'libmsym' / 'character_table.c',
    'libmsym' / 'context.c',
    'libmsym' / 'debug.c',
    'libmsym' / 'elements.c',
    'libmsym' / 'equivalence_set.c',
    'libmsym' / 'geometry.c',
    'libmsym' / 'linalg.c',
    'libmsym' / 'msym.c',
    'libmsym' / 'msym_error.c',
    'libmsym' / 'permutation.c',
    'libmsym' / 'point_group.c',
    'libmsym' / 'rsh.c',
    'libmsym' / 'subspace.c',
    'libmsym' / 'symmetrize.c',
    'libmsym' / 'symmetry.c',
    'libmsym' / 'symop.c',
  ),
  include_directories: msym_inc,
  c_args: msym_c_args + ['-DMSYM_EXPORTS_BUILT_AS_STATIC'],
  override_options: ['c_std=c99'],
  implicit_include_directories: false,
  # Nothing of libmsym leaves libmoist
  gnu_symbol_visibility: 'hidden',
)

# moist's own shim keeps the project's warnings
msym_lib = static_library(
  'moist_msym',
  sources: files('msym_shim.c'),
  include_directories: msym_inc,
  c_args: ['-DMSYM_EXPORTS_BUILT_AS_STATIC'],
  implicit_include_directories: false,
  gnu_symbol_visibility: 'hidden',
)

srcs += files('msym.f90', 'frame.f90', 'embed.f90', 'orient.f90')

if install
  install_data(
    'libmsym' / 'LICENSE',
    install_dir: get_option('datadir') / 'licenses' / meson.project_name(),
    rename: 'LICENSE-libmsym',
  )
endif

lib_deps += declare_dependency(
  # Shim before libmsym: static archives are searched in order
  link_with: [msym_lib, libmsym_lib],
  dependencies: meson.get_compiler('c').find_library('m', required: false),
)
