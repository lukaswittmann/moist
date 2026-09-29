# FFT backend: moist's extern "C" shim over the vendored ducc0 subset, and the
# vendored xsimd and FINUFFT built on top of it (one meson.build per vendor
# directory; no CMake and no network access at configure time)
#
# - static libraries with implicit_include_directories: false keep the source
#   root off the C++ include path (libc++'s <version> header would otherwise
#   resolve to the project's VERSION file on case-insensitive file systems)
# - lib_deps gets link_with only, so Fortran compiles never see a C++ include
#   directory
fft_inc = include_directories('.')
# Threads and, when enabled, OpenMP for every C++ target of the backend
fft_cpp_deps = [threads_dep]
if get_option('openmp')
  fft_cpp_deps += dependency('openmp', language: 'cpp')
endif

fft_lib = static_library(
  'moist_fft',
  sources: files(
    'ducc0.cpp',
    'ducc0' / 'infra' / 'string_utils.cc',
    'ducc0' / 'infra' / 'threading.cc',
    'ducc0' / 'infra' / 'mav.cc',
  ),
  include_directories: fft_inc,
  dependencies: fft_cpp_deps,
  override_options: ['cpp_std=c++17'],
  implicit_include_directories: false,
)

srcs += files('dst4.f90')

if install
  install_data(
    'ducc0' / 'LICENSE',
    install_dir: get_option('datadir') / 'licenses' / meson.project_name(),
    rename: 'LICENSE-ducc0',
  )
endif

subdir('xsimd')
subdir('finufft')

# FINUFFT before moist_fft: static archives are searched in order, and FINUFFT
# uses the ducc0 infra objects of moist_fft
lib_deps += finufft_dep
lib_deps += declare_dependency(
  link_with: fft_lib,
  dependencies: fft_cpp_deps,
)
