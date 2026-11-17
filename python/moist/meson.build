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

# The build-tree rpath does not survive the move into site-packages, so the
# installed module has to be pointed at libmoist explicitly. An external
# libmoist stays where pkg-config found it; one built here is installed next to
# the module's prefix, and a path relative to the module keeps that prefix
# relocatable. A static libmoist needs no run-time path at all.
moist_rpath = ''
if moist_dep.type_name() == 'pkgconfig'
  moist_rpath = moist_dep.get_variable(pkgconfig: 'libdir')
elif moist_dep.type_name() == 'internal' and moist_dep.get_variable(
  internal: 'shared',
  default_value: 'false',
) == 'true'
  moist_libdir = moist_dep.get_variable(internal: 'libdir')
  if meson.version().version_compare('>=1.3.0')
    moist_origin = host_machine.system() == 'darwin' ? '@loader_path' : '$ORIGIN'
    moist_rpath = moist_origin / import('fs').relative_to(
      moist_libdir,
      python.get_install_dir(pure: false) / 'moist',
    )
  else
    moist_rpath = moist_libdir
  endif
endif

moist_pyext = python.extension_module(
  '_libmoist',
  link_whole: static_library(
    '_libmoist',
    moist_cffi_srcs,
    dependencies: [moist_dep, python_dep],
  ),
  dependencies: [moist_dep, python_dep],
  install_rpath: moist_rpath,
  install: install,
  subdir: 'moist',
)

# Pure-Python modules of the package
pysrcs = [
  '__init__.py',
  '_pyscf_hessian.py',
  'configuration.py',
  'context.py',
  'density.py',
  'interface.py',
  'library.py',
  'parameters.py',
  'pyscf.py',
  'radii.py',
  'solvent.py',
]

if install
  python.install_sources(
    files(pysrcs),
    subdir: 'moist',
    pure: false,
  )
endif

foreach name : pysrcs
  configure_file(input: files(name), output: name, copy: true)
endforeach

test(
  'python_requirements',
  python,
  args: [
    files('..'/'check_requirements.py'),
    'cffi',
    'numpy',
    'pytest',
    'pyscf:dft,grad,gto,hessian,scf',
  ],
  suite: 'python',
  timeout: 60,
)

# pytest options
tests_dir = meson.current_source_dir() / '..' / 'tests'
pytest_args = [
  '-m', 'pytest',
  '--no-header', '-q',
]
pytest_env = ['PYTHONPATH=' + meson.current_build_dir() / '..']
coupling_env = pytest_env + ['MOIST_REQUIRE_PYSCF=1']
coupling_suites = ['python', 'pyscf']
bindings_srcs = [
  tests_dir / 'test_interface.py',
  tests_dir / 'test_configuration.py',
  tests_dir / 'test_api_contracts.py',
  tests_dir / 'test_library.py',
  tests_dir / 'test_solvent.py',
  tests_dir / 'test_parameter_adjoint.py',
]
coupling_src = tests_dir / 'test_pyscf.py'

# The CFFI bindings and the Python object model
test(
  'python_api',
  python,
  args: pytest_args + bindings_srcs,
  depends: moist_pyext,
  env: pytest_env,
  suite: 'python',
  timeout: 60,
)

# Test MOIST--QM coupling (using PySCF) (mainly isodensity cavity + models)
#
#   host          the adapter's own derivative orders and basis conversion
#   vdw           electrostatic conventions on a density-independent cavity
#   isodensity_*  the level-set chain rule, one target per component set
#   conventions   channel semantics and the negative controls
#   scf           the self-consistent layer (slowest)
foreach entry : [
  ['pyscf_host', 'host', 300],
  ['pyscf_vdw', 'vdw', 300],
  ['pyscf_isodensity_cpcm', 'isodensity and cpcm', 300],
  ['pyscf_isodensity_pv', 'isodensity and pv', 300],
  ['pyscf_isodensity_cpcm_pv', 'isodensity and cpcm_pv', 300],
  ['pyscf_conventions', 'conventions', 300],
  ['pyscf_scf', 'scf', 600],
]
  test(
    entry[0],
    python,
    args: pytest_args + ['-m', entry[1], coupling_src],
    depends: moist_pyext,
    env: coupling_env,
    suite: coupling_suites,
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
gostshyp_src = tests_dir / 'test_gostshyp.py'
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
    env: coupling_env,
    suite: coupling_suites,
    timeout: entry[2],
  )
endforeach

# ECP nuclei: the isodensity core fill against a pseudo-atom and all-electron AuH
test(
  'pyscf_ecp',
  python,
  args: pytest_args + [tests_dir / 'test_pyscf_ecp.py'],
  depends: moist_pyext,
  env: coupling_env,
  suite: coupling_suites,
  timeout: 300,
)

# Public PySCF SCF/gradient/Hessian wrapper, including cache and scanner lifetimes
test(
  'pyscf_wrapper',
  python,
  args: pytest_args + [tests_dir / 'test_pyscf_scf.py'],
  depends: moist_pyext,
  env: coupling_env,
  suite: coupling_suites,
  timeout: 600,
)

# Analytical PySCF nuclear Hessians and directional Hessian-vector products
test(
  'pyscf_hessian',
  python,
  args: pytest_args + [tests_dir / 'test_pyscf_hessian.py'],
  depends: moist_pyext,
  env: coupling_env,
  suite: coupling_suites,
  timeout: 600,
)
