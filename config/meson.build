os = host_machine.system()
fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')
fc_id = fc.get_id()

# Treat Intel LLVM compiler as standard Intel compiler for compatibility
if fc_id == 'intel-llvm'
  fc_id = 'intel'
endif

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )

  libheapt_w = fc.find_library('heapt_w', required: false)
  if libheapt_w.found()
    lib_deps += libheapt_w
  endif
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    '-fp-model=precise',
    '-fpscomp', 'logicals',
    language: 'fortran',
  )
  # Intel traceback support at runtime requires libunwind to be linked.
  unwind_dep = cc.find_library('unwind', required: false)
  if unwind_dep.found()
    lib_deps += unwind_dep
  else
    warning('libunwind not found; Intel traceback output may be unavailable at runtime')
  endif

  if get_option('openmp')
    cc.find_library('ifcoremt')
    add_project_link_arguments(
      '-Wl,--no-as-needed', '-lifcoremt', '-Wl,--as-needed',
      language: ['fortran', 'c', 'cpp'],
    )
  endif
elif fc_id in ['pgi', 'nvidia_hpc']
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

foreach lang, comp : {'c': cc, 'cpp': cxx}
  if comp.get_id() in ['intel', 'intel-llvm']
    add_project_arguments('-fp-model=precise', language: lang)
  endif
endforeach

if get_option('openmp')
  omp_dep = dependency('openmp', language: 'cpp')
  lib_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto' and fc_id == 'intel'
  lapack_vendor = 'mkl'
endif

ilp64 = get_option('ilp64')
sfx64 = ilp64 ? '64' : ''
add_project_arguments(
  ilp64 ? ['-DIK=i8', '-DMOIST_ILP64=1'] : ['-DIK=i4', '-DMOIST_ILP64=0'],
  language: 'fortran',
)

if lapack_vendor == 'mkl'
  if fc_id not in ['intel', 'gcc']
    error('MKL not supported for this compiler')
  endif

  mkl_interface = fc_id == 'intel' ? 'mkl_intel' : 'mkl_gf'
  mkl_openmp = fc_id == 'intel' ? 'mkl_intel_thread' : 'mkl_gnu_thread'
  foreach lib : [
    mkl_interface + (ilp64 ? '_ilp64' : '_lp64'),
    get_option('openmp') ? mkl_openmp : 'mkl_tbb_thread',
    'mkl_core',
  ]
    lib_deps += cc.find_library(lib)
  endforeach

elif lapack_vendor == 'mkl-rt'
  lib_deps += cc.find_library('mkl_rt')

elif lapack_vendor == 'accelerate'
  if os != 'darwin'
    error('The accelerate backend is only available on macOS')
  endif
  if ilp64
    error('The accelerate backend provides the LP64 interface only; build without -Dilp64')
  endif
  lib_deps += dependency('appleframeworks', modules: ['Accelerate'])

elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas' + sfx64, required: false)
  if not openblas_dep.found()
    openblas_dep = cc.find_library('openblas' + sfx64)
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack' + sfx64, required: false)
    if not lapack_dep.found()
      lapack_dep = cc.find_library('lapack' + sfx64)
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  # An optional leading -L<dir> is the search directory for the names after it
  libs = get_option('custom_libraries')
  custom_dirs = libs[0].startswith('-L') ? [libs[0].substring(2)] : []
  custom_deps = []
  foreach lib : libs
    if not lib.startswith('-L')
      custom_deps += cc.find_library(lib, dirs: custom_dirs)
    endif
  endforeach
  custom_check_deps = get_option('openmp') ? [custom_deps, omp_dep] : custom_deps
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: custom_check_deps)
    error('Custom LAPACK libraries do not link')
  elif not fc.links('external dgemm; call dgemm(); end', dependencies: custom_check_deps)
    error('Custom BLAS libraries do not link')
  endif
  lib_deps += custom_deps

else
  lapack_dep = dependency('lapack' + sfx64, required: false)
  if not lapack_dep.found()
    lapack_dep = cc.find_library('lapack' + sfx64)
  endif
  lib_deps += [lapack_dep, cc.find_library('blas' + sfx64)]
endif

tomlf_dep = dependency(
  'toml-f',
  version: '>=0.4.0',
  fallback: ['toml-f', 'tomlf_dep'],
  default_options: ['default_library=static'],
)
lib_deps += tomlf_dep

jonquil_dep = dependency(
  'jonquil',
  version: '>=0.3.0',
  fallback: ['jonquil', 'jonquil_dep'],
  default_options: ['default_library=static'],
)
lib_deps += jonquil_dep

# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  version: '>=0.3.0',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
)
lib_deps += mctc_dep

# fclap dependency for CLI (*not* added to lib_deps)
fclap_dep = dependency(
  'fclap',
  fallback: ['fclap', 'fclap_dep'],
  default_options: ['default_library=static'],
)
