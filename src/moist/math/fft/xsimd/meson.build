# xsimd (BSD-3-Clause): vendored header-only SIMD wrappers for FINUFFT's
# spreader, see PROVENANCE.md. This file is moist's, not upstream's.
xsimd_dep = declare_dependency(
  include_directories: include_directories('include'),
)

if install
  install_data(
    'LICENSE',
    install_dir: get_option('datadir') / 'licenses' / meson.project_name(),
    rename: 'LICENSE-xsimd',
  )
endif
