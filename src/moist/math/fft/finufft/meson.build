# FINUFFT (Apache-2.0): vendored non-uniform FFT, see PROVENANCE.md. This file
# is moist's, not upstream's; it replaces FINUFFT's own CMake build.
#
# - CPU FFTs on moist's ducc0 subset (FINUFFT_USE_DUCC0, headers via fft_inc);
#   the ducc0 infra objects come from moist_fft and are never compiled twice
# - moist's portable project flags only: no -march, no fast-math
# - mirrors upstream src/CMakeLists.txt: the precision sources are compiled
#   twice (double, and float with FINUFFT_SINGLE); the common sources and the
#   src/common helpers once
# - every header path (include/, xsimd, ducc0) stays on these targets; Fortran
#   compiles only link against them
finufft_inc = [include_directories('include'), fft_inc]
finufft_args = ['-DFINUFFT_USE_DUCC0']
# Several precision headers use #pragma omp, so every target needs OpenMP
finufft_deps = [xsimd_dep] + fft_cpp_deps

finufft_precision_cpp = files(
  'src' / 'makeplan.cpp',
  'src' / 'setpts.cpp',
  'src' / 'execute.cpp',
  'src' / 'spreadinterp.cpp',
  'src' / 'spreadinterp_1d.cpp',
  'src' / 'spreadinterp_2d.cpp',
  'src' / 'spreadinterp_3d.cpp',
)

finufft_helper_lib = static_library(
  'moist_finufft_common',
  sources: files(
    'src' / 'common' / 'kernel.cpp',
    'src' / 'common' / 'pswf.cpp',
    'src' / 'common' / 'utils.cpp',
  ),
  include_directories: finufft_inc,
  cpp_args: finufft_args,
  dependencies: finufft_deps,
  override_options: ['cpp_std=c++17'],
  implicit_include_directories: false,
)

finufft_f32_lib = static_library(
  'moist_finufft_f32',
  sources: finufft_precision_cpp,
  include_directories: finufft_inc,
  cpp_args: finufft_args + ['-DFINUFFT_SINGLE'],
  dependencies: finufft_deps,
  override_options: ['cpp_std=c++17'],
  implicit_include_directories: false,
)

finufft_lib = static_library(
  'moist_finufft',
  sources: finufft_precision_cpp + files(
    'src' / 'fft.cpp',
    'src' / 'c_interface.cpp',
    'src' / 'utils.cpp',
    'fortran' / 'finufftfort.cpp',
  ),
  include_directories: finufft_inc,
  cpp_args: finufft_args,
  dependencies: finufft_deps,
  link_with: [finufft_helper_lib, finufft_f32_lib, fft_lib],
  override_options: ['cpp_std=c++17'],
  implicit_include_directories: false,
)

finufft_dep = declare_dependency(
  link_with: [finufft_lib, finufft_helper_lib, finufft_f32_lib],
  dependencies: fft_cpp_deps,
)

# Fortran interface module (`use finufft_mod`), compiled with the library
srcs += files('include' / 'finufft_mod.f90')

if install
  install_data(
    'LICENSE',
    'NOTICE',
    install_dir: get_option('datadir') / 'licenses' / meson.project_name(),
    rename: ['LICENSE-finufft', 'NOTICE-finufft'],
  )
endif
