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
  'conftest.py',
  'gostshyp.py',
  'interface.py',
  'library.py',
  'pyscf.py',
  'test_gostshyp.py',
  'test_interface.py',
  'test_library.py',
  'test_pyscf.py',
)

if install
  python.install_sources(
    pysrcs,
    subdir: 'moist',
    pure: false,
  )
endif

foreach name : [
  '__init__.py',
  'conftest.py',
  'gostshyp.py',
  'interface.py',
  'library.py',
  'pyscf.py',
  'test_gostshyp.py',
  'test_interface.py',
  'test_library.py',
  'test_pyscf.py',
]
  configure_file(input: files(name), output: name, copy: true)
endforeach

# pytest is an optional development dependency
if run_command(python, '-c', 'import pytest', check: false).returncode() == 0
  pytest_args = [
    '-m', 'pytest',
    '--no-header', '-q',
    '-p', 'no:cacheprovider',
  ]
  pytest_env = ['PYTHONPATH=' + meson.current_build_dir() / '..']
  bindings_srcs = [
    meson.current_build_dir() / 'test_interface.py',
    meson.current_build_dir() / 'test_library.py',
  ]
  coupling_src = meson.current_build_dir() / 'test_pyscf.py'

  # The CFFI bindings and the Python object model
  test(
    'python_api',
    python,
    args: pytest_args + bindings_srcs,
    depends: moist_pyext,
    env: pytest_env,
    timeout: 60,
  )

  # Test MOIST--QM coupling (using PySCF) (mainly isodensity cavity + models)
  #
  #   host          the adapter's own derivative orders, moist not involved
  #   vdw           electrostatic conventions on a density-independent cavity
  #   isodensity_*  the level-set chain rule, one target per component set
  #   conventions   channel semantics and the negative controls
  #   scf           the self-consistent layer (slowest)
  #
  # The timeouts are wall-clock and meson runs the targets concurrently, so they
  # are set for the slowest CI runner under full load rather than for the times
  # a developer sees locally -- a minute was not enough for either isodensity
  # tier or for scf on the Ubuntu runners.
  foreach entry : [
    ['pyscf_host', 'host', 300],
    ['pyscf_vdw', 'vdw', 300],
    ['pyscf_isodensity_cpcm', 'isodensity and cpcm', 300],
    ['pyscf_isodensity_pv', 'isodensity and pv', 300],
    ['pyscf_isodensity_cpcm_pv', 'isodensity and cpcm_pv', 300],
    ['pyscf_conventions', 'conventions', 300],
    ['pyscf_scf', 'scf', 300],
  ]
    test(
      entry[0],
      python,
      args: pytest_args + ['-m', entry[1], coupling_src],
      depends: moist_pyext,
      env: pytest_env,
      timeout: entry[2],
    )
  endforeach

  # The GOSTSHYP pressure model on the isodensity cavity
  #
  #   gostshyp_integrals   libcint fakemol constants and cartesian orders
  #   gostshyp_params      surface-parameter derivatives at a frozen cavity
  #   gostshyp_fock        energy and dE/dP on the moving surface
  #   gostshyp_gradient    the three nuclear-gradient channels
  #   gostshyp_conventions channel semantics and the negative controls
  gostshyp_src = meson.current_build_dir() / 'test_gostshyp.py'
  foreach entry : [
    ['gostshyp_integrals', 'gostshyp_integrals', 300],
    ['gostshyp_params', 'gostshyp_params', 300],
    ['gostshyp_fock', 'gostshyp_fock', 300],
    ['gostshyp_gradient', 'gostshyp_gradient', 300],
    ['gostshyp_conventions', 'gostshyp_conventions', 300],
  ]
    test(
      entry[0],
      python,
      args: pytest_args + ['-m', entry[1], gostshyp_src],
      depends: moist_pyext,
      env: pytest_env,
      timeout: entry[2],
    )
  endforeach

else
  message('pytest not found: skipping the python test targets')
endif
