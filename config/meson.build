os = host_machine.system()
fc = meson.get_compiler('fortran')
cc = fc
if has_cc
  cc = meson.get_compiler('c')
endif
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
  # GCC 14+ emits heap-trampoline support calls (__gcc_nested_func_ptr_created/deleted)
  # for internal procedures passed as callbacks. These symbols live in libgcc_eh
  # (or libheapt_w), not the main libgcc. When the linker is clang (e.g. conda-forge
  # toolchains on macOS), neither library is pulled in automatically.
  libgcc_eh = fc.find_library('gcc_eh', required: false)
  if libgcc_eh.found()
    lib_deps += libgcc_eh
  else
    libgcc = fc.find_library('gcc', required: false)
    if libgcc.found()
      lib_deps += libgcc
    endif
  endif
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    # ifx defaults to -fp-model=fast, which reassociates reductions and contracts
    # FMAs aggressively. That diverges from gfortran's value-safe default and
    # rounds the analytic-vs-reference / quadrature-identity tests over their
    # tolerances. Pin a precise model so Intel numerics are reproducible.
    '-fp-model=precise',
    language: 'fortran',
  )
  # Intel traceback support at runtime requires libunwind to be linked.
  unwind_dep = cc.find_library('unwind', required: false)
  if unwind_dep.found()
    lib_deps += unwind_dep
  else
    warning('libunwind not found; Intel traceback output may be unavailable at runtime')
  endif
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

# RISM needs FFTW for the reciprocal-space transforms and HDF5 for caching
use_rism = get_option('rism')
use_fftw = get_option('fftw') or use_rism
use_hdf5 = get_option('hdf5') or use_rism

if use_fftw
  fftw_dep = dependency('fftw3', required: true)
  lib_deps += fftw_dep

  if get_option('openmp')
    fftw_threads_dep = dependency('fftw3_threads', required: false)
    if fftw_threads_dep.found()
      lib_deps += fftw_threads_dep
    endif
  endif
  add_project_arguments('-DWITH_FFTW=1', language: 'fortran')
endif

if use_rism
  add_project_arguments('-DWITH_RISM=1', language: 'fortran')
endif

if use_hdf5
  lib_deps += dependency('hdf5', language: 'fortran', required: true)
  add_project_arguments('-DWITH_HDF5=1', language: 'fortran')
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel'
    lapack_vendor = 'mkl'
  endif
endif

if get_option('ilp64')
  ilp64 = true
  add_project_arguments('-DIK=i8', language:'fortran')
  add_project_arguments('-DMOIST_ILP64=1', language:'fortran')
else
  add_project_arguments('-DIK=i4', language:'fortran')
  add_project_arguments('-DMOIST_ILP64=0', language:'fortran')
  ilp64 = false
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel'
    mkl_dep += cc.find_library(ilp64 ? 'mkl_intel_ilp64' : 'mkl_intel_lp64') 
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += cc.find_library(ilp64 ? 'mkl_gf_ilp64' : 'mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += cc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += cc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = cc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'accelerate'
  # Apple's Accelerate framework ships BLAS/LAPACK with macOS, so nothing needs
  # to be installed. It exposes the legacy LP64 (32-bit integer) interface.
  if os != 'darwin'
    error('The accelerate backend is only available on macOS')
  endif
  if ilp64
    error('The accelerate backend provides the LP64 interface only; build without -Dilp64')
  endif
  lib_deps += dependency('appleframeworks', modules: ['Accelerate'])

elif lapack_vendor == 'openblas'
  openblas_dep = dependency( ilp64 ? 'openblas64' : 'openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = cc.find_library( ilp64 ? 'openblas64' : 'openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency(ilp64 ? 'lapack64' : 'lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = cc.find_library(ilp64 ? 'lapack64' : 'lapack')
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  custom_deps = []
  libs = get_option('custom_libraries')
  if libs[0].startswith('-L')
    foreach lib: libs
      if lib != libs[0]
        custom_deps += cc.find_library(lib, dirs: libs[0].substring(2))
      endif
    endforeach
  else
    foreach lib: libs
      custom_deps += cc.find_library(lib)
    endforeach
  endif
  if (not fc.links('external dsytrs; call dsytrs(); end', dependencies: (get_option('openmp') ? [custom_deps, omp_dep] : [custom_deps])))
    error('Custom LAPACK libraries do not link')
  elif (not fc.links('external dsytrs; call dgemm(); end', dependencies: (get_option('openmp') ? [custom_deps, omp_dep] : [custom_deps])))
    error('Custom BLAS libraries do not link')
  endif
  lib_deps += custom_deps

else
  lapack_dep = dependency(ilp64 ? 'lapack64' : 'lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = cc.find_library(ilp64 ? 'lapack64' : 'lapack')
  endif
  lib_deps += lapack_dep
  blas_dep = cc.find_library(ilp64 ? 'blas64' : 'blas', required: false)
  if not blas_dep.found()
    blas_dep = cc.find_library(ilp64 ? 'blas64' : 'blas')
  endif
  lib_deps += blas_dep
endif

# JSON parameter IO (used directly by moist_model_parameters)
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

# Create the test-drive library as subproject
testdrive_dep = dependency(
  'test-drive',
  version: '>=0.5.0',
  fallback: ['test-drive', 'testdrive_dep'],
  default_options: ['default_library=static'],
)
lib_deps += testdrive_dep

# Create dftd4 as subproject
dftd4_dep = dependency(
  'dftd4',
  version: '>=4.0.0',
  fallback: ['dftd4', 'dftd4_dep'],
  default_options: ['default_library=static'],
)
lib_deps += dftd4_dep

# Create the electrostatic library as subproject
multicharge_dep = dependency(
  'multicharge',
  version: '>=0.4.0',
  fallback: ['multicharge', 'multicharge_dep'],
  default_options: ['default_library=static'],
)
lib_deps += multicharge_dep

# fclap dependency for CLI (*not* added to lib_deps)
fclap_dep = dependency(
  'fclap',
  fallback: ['fclap', 'fclap_dep'],
  default_options: ['default_library=static'],
  required: true,
)