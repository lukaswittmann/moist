# CMake build for moist; the meson build (meson.build, config/meson.build) is the reference.
#   cmake -B build-cmake -G Ninja && cmake --build build-cmake

cmake_minimum_required(VERSION 3.25)

# project(VERSION) accepts only the numeric part of VERSION (e.g. 1.0.0-alpha.1)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" MOIST_VERSION_FULL)
string(STRIP "${MOIST_VERSION_FULL}" MOIST_VERSION_FULL)
string(REGEX REPLACE "-.*$" "" MOIST_VERSION "${MOIST_VERSION_FULL}")

project(
  moist
  VERSION "${MOIST_VERSION}"
  DESCRIPTION "Modular and Open-source Implicit Solvation Toolkit"
  LANGUAGES Fortran C CXX
)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# icx/icpx default to a fast FP model, like ifx
foreach(lang C CXX)
  if(CMAKE_${lang}_COMPILER_ID MATCHES "Intel")
    add_compile_options($<$<COMPILE_LANGUAGE:${lang}>:-fp-model=precise>)
  endif()
endforeach()

# Options (mirror meson_options.txt)
set(MOIST_LAPACK "auto" CACHE STRING "Linear algebra backend")
set_property(CACHE MOIST_LAPACK PROPERTY STRINGS
  auto mkl mkl-rt openblas netlib accelerate custom)
set(MOIST_CUSTOM_LIBRARIES "" CACHE STRING
  "Libraries for the custom LAPACK backend: paths, names or -L/-l flags")
option(MOIST_OPENMP "Use OpenMP parallelisation"            ON)
option(MOIST_API    "Build C API using iso_c_binding"       ON)
option(MOIST_ILP64  "Enable BLAS/LAPACK ilp64 support"      OFF)
option(MOIST_TESTS  "Build the moist unit test suite (ctest)" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
endif()

# BLAS/LAPACK/OpenMP first: the subprojects' `if(NOT TARGET ...)` guards then reuse them
if(MOIST_LAPACK STREQUAL "custom")
  if(NOT MOIST_CUSTOM_LIBRARIES)
    message(FATAL_ERROR "MOIST_LAPACK=custom requires -DMOIST_CUSTOM_LIBRARIES=...")
  endif()
  # One library set provides both
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${MOIST_CUSTOM_LIBRARIES}")
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES
    INTERFACE_LINK_LIBRARIES "BLAS::BLAS")
else()
  if(MOIST_LAPACK STREQUAL "mkl")
    set(BLA_VENDOR Intel10_64lp)
  elseif(MOIST_LAPACK STREQUAL "mkl-rt")
    set(BLA_VENDOR Intel10_64_dyn)
  elseif(MOIST_LAPACK STREQUAL "openblas")
    set(BLA_VENDOR OpenBLAS)
  elseif(MOIST_LAPACK STREQUAL "netlib")
    set(BLA_VENDOR Generic)
  elseif(MOIST_LAPACK STREQUAL "accelerate")
    if(NOT APPLE)
      message(FATAL_ERROR "The accelerate backend is only available on macOS")
    endif()
    if(MOIST_ILP64)
      message(FATAL_ERROR
        "The accelerate backend provides the LP64 interface only; "
        "build without -DMOIST_ILP64")
    endif()
    set(BLA_VENDOR Apple)
  endif()
  if(MOIST_ILP64)
    set(BLA_SIZEOF_INTEGER 8)
  endif()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

if(MOIST_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran CXX)
endif()
find_package(Threads REQUIRED)

# C++ FFT backend (ducc0 subset + shim) and FINUFFT on top of it, as OBJECT
# libraries so their objects land in libmoist itself (static and shared)
set(_fft_dir "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/math/fft")
set(_finufft_dir "${_fft_dir}/finufft")
add_library(moist_fft OBJECT
  "${_fft_dir}/ducc0.cpp"
  "${_fft_dir}/ducc0/infra/string_utils.cc"
  "${_fft_dir}/ducc0/infra/threading.cc"
  "${_fft_dir}/ducc0/infra/mav.cc"
)
target_include_directories(moist_fft PUBLIC "${_fft_dir}")

# Like upstream: the precision sources once for double and once with FINUFFT_SINGLE
set(FINUFFT_PRECISION_SOURCES
  "${_finufft_dir}/src/makeplan.cpp"
  "${_finufft_dir}/src/setpts.cpp"
  "${_finufft_dir}/src/execute.cpp"
  "${_finufft_dir}/src/spreadinterp.cpp"
  "${_finufft_dir}/src/spreadinterp_1d.cpp"
  "${_finufft_dir}/src/spreadinterp_2d.cpp"
  "${_finufft_dir}/src/spreadinterp_3d.cpp"
)
add_library(finufft_common OBJECT
  "${_finufft_dir}/src/common/kernel.cpp"
  "${_finufft_dir}/src/common/pswf.cpp"
  "${_finufft_dir}/src/common/utils.cpp"
)
add_library(finufft_f32 OBJECT ${FINUFFT_PRECISION_SOURCES})
target_compile_definitions(finufft_f32 PRIVATE FINUFFT_SINGLE)
add_library(finufft OBJECT ${FINUFFT_PRECISION_SOURCES}
  "${_finufft_dir}/src/fft.cpp"
  "${_finufft_dir}/src/c_interface.cpp"
  "${_finufft_dir}/src/utils.cpp"
  "${_finufft_dir}/fortran/finufftfort.cpp"
)
foreach(tgt finufft_common finufft_f32 finufft)
  target_include_directories(${tgt} PRIVATE
    "${_finufft_dir}/include" "${_fft_dir}/xsimd/include" "${_fft_dir}")
  target_compile_definitions(${tgt} PRIVATE FINUFFT_USE_DUCC0)
endforeach()
foreach(tgt moist_fft finufft_common finufft_f32 finufft)
  target_link_libraries(${tgt} PUBLIC Threads::Threads)
  if(MOIST_OPENMP)
    target_link_libraries(${tgt} PUBLIC OpenMP::OpenMP_CXX)
  endif()
endforeach()

# Vendored Fortran subprojects, in dependency order; the `<name>::<name>` aliases
# satisfy their own `if(NOT TARGET ...)` lookups, so nothing is fetched
set(BUILD_TESTING            OFF CACHE BOOL "" FORCE)  # toml-f, fclap, mstore, ...
set(ENABLE_TESTING           OFF CACHE BOOL "" FORCE)  # jonquil
set(TESTDRIVE_BUILD_TESTING  OFF CACHE BOOL "" FORCE)  # test-drive

# mctc-lib registers its tests unconditionally; only moist's suite may add tests
set(_moist_allow_tests OFF)
function(add_test)
  if(_moist_allow_tests)
    _add_test(${ARGV})
  endif()
endfunction()

function(moist_add_subproject name)
  if(NOT TARGET "${name}")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/subprojects/${name}" EXCLUDE_FROM_ALL)
  endif()
  if(NOT TARGET "${name}::${name}")
    add_library("${name}::${name}" ALIAS "${name}")
  endif()
endfunction()

set(_moist_shared_library "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS OFF)
moist_add_subproject(toml-f)
moist_add_subproject(fclap)
moist_add_subproject(jonquil)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(jonquil-lib PRIVATE -frecursive)  # recursive JSON visitor
endif()
moist_add_subproject(mctc-lib)
set(BUILD_SHARED_LIBS "${_moist_shared_library}")

# Globbed sources (including the vendored finufft_mod.f90); CMake orders the modules
file(GLOB_RECURSE moist_srcs CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.F90"
)
if(NOT MOIST_API)
  list(REMOVE_ITEM moist_srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/api.f90")
endif()

find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE VCS_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()
if(NOT VCS_TAG)
  set(VCS_TAG "unknown")
endif()
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/build_info.f90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.f90"
  @ONLY
)
list(APPEND moist_srcs "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.f90")

add_library(moist ${moist_srcs})
set_target_properties(moist PROPERTIES
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
  # A Fortran-driven link turns C++ exceptions from ducc0/FINUFFT into aborts
  LINKER_LANGUAGE CXX
)
target_include_directories(moist PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
if(MOIST_ILP64)
  target_compile_definitions(moist PRIVATE IK=i8 MOIST_ILP64=1)
else()
  target_compile_definitions(moist PRIVATE IK=i4 MOIST_ILP64=0)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(moist PRIVATE -ffree-line-length-none -fbacktrace -fopenmp-simd)
  # GCC 14+ heap trampolines; never libgcc_eh, whose unwinder aborts C++ exceptions
  find_library(MOIST_LIBHEAPT_W NAMES heapt_w
    HINTS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  if(MOIST_LIBHEAPT_W)
    target_link_libraries(moist PUBLIC ${MOIST_LIBHEAPT_W})
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(MOIST_INTEL_FFLAGS -traceback -fp-model=precise -fpscomp logicals)
  target_compile_options(moist PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${MOIST_INTEL_FFLAGS}>")
endif()

target_link_libraries(moist PUBLIC
  mctc-lib::mctc-lib
  jonquil::jonquil
  toml-f::toml-f
  LAPACK::LAPACK
  BLAS::BLAS
)
# PRIVATE keeps the C++ include directories off moist's consumers
target_link_libraries(moist PRIVATE moist_fft finufft finufft_f32 finufft_common)
if(MOIST_OPENMP)
  target_link_libraries(moist PUBLIC OpenMP::OpenMP_Fortran)
endif()

# ifx takes `main` from its own runtime (for_main.o), which only its driver links
set(MOIST_EXE_LINKER_LANGUAGE CXX)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(MOIST_EXE_LINKER_LANGUAGE Fortran)
endif()

add_executable(moist-exe
  app/main.f90
  app/cli.f90
  app/driver.f90
)
set_target_properties(moist-exe PROPERTIES
  OUTPUT_NAME moist
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/app-modules"
  LINKER_LANGUAGE ${MOIST_EXE_LINKER_LANGUAGE}
)
target_link_libraries(moist-exe PRIVATE moist fclap)

include(GNUInstallDirs)
install(TARGETS moist moist-exe
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
if(MOIST_API)
  install(FILES include/moist.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/moist"
  FILES_MATCHING PATTERN "*.mod")

set(_licenses "${CMAKE_INSTALL_DATADIR}/licenses/moist")
install(FILES LICENSE THIRD_PARTY_NOTICES.md DESTINATION "${_licenses}")
foreach(vendor ducc0 finufft xsimd)
  install(FILES "${_fft_dir}/${vendor}/LICENSE"
    DESTINATION "${_licenses}" RENAME "LICENSE-${vendor}")
endforeach()
install(FILES "${_finufft_dir}/NOTICE" DESTINATION "${_licenses}" RENAME NOTICE-finufft)

# Unit tests: one tester binary, run once per suite (mirrors test/unit/meson.build)
if(MOIST_TESTS)
  enable_testing()
  set(_moist_allow_tests ON)
  set(_unit "${CMAKE_CURRENT_SOURCE_DIR}/test/unit")

  moist_add_subproject(test-drive)
  moist_add_subproject(mstore)
  if(MOIST_OPENMP)
    foreach(tgt test-drive-lib mstore-lib)
      if(TARGET ${tgt})
        target_link_libraries(${tgt} PRIVATE OpenMP::OpenMP_Fortran)
      endif()
    endforeach()
  endif()

  # Kept in sync with the `tests` list in test/unit/meson.build
  set(moist_test_suites
    api parameters utils_timer utils_context utils_mem utils_prettylistprint utils_prettylistprint_threaded
    data radii
    math_linalg math_smoothing_kernels math_adjacency_list math_cell_grid
    math_solvers math_sorters math_trig
    math_fft_dst4 math_fft_3d finufft math_grid_3d_fft math_grid_3d_nufft
    math_grid_radial_rule math_grid_radial math_grid_radial_trafo
    math_grid_angular math_grid_atomic
    math_grid_3d_threaded moz_1d moz_3d
    moz_potential_kernel moz_potential_electrostatic moz_potential_lj moz_potential_oplsaa moz_solvent
    math_grid_3d_partition
    math_grid_3d_guards math_grid_3d_derivatives
    math_grid_3d_integration
    math_grid_3d_invariance
    cavity_drop_primitives cavity_drop_cfc cavity_drop_lsf cavity_drop_lsf_reference
    cavity_drop_isodensity cavity_drop_gradient cavity_drop_nuclear_adjoint
    cavity_drop_integration cavity_drop_filter cavity_drop_cpcm
    cavity_iswig cavity_numsa cavity_marchingcubes cavity_reference cavity_surface_adjoint
    model_component_pcm_amat model_component_pcm_amat_kernel
    model_component_pcm_amat_assembly model_component_pcm_amat_adjoint
    model_component_pcm_electrostatics model_component_pcm_cpcm
    model_component_gostshyp model_component_pv
    model_continuum channels model_coupling
  )
  # Each run as a selected test, outside test-drive's OpenMP team
  set(moist_test_selected_utils_prettylistprint_threaded independent_rows)
  set(moist_test_selected_math_grid_3d_derivatives
    single_atom dimer amino20x4_gly_xab mb16_43_01 but14diol_1 il16_008
    partition_schemes gaussian_tiny_weights)
  set(moist_test_selected_math_grid_3d_integration
    single_atom asymmetric_trimer partition_integrals multicenter_integrals
    rotational_convergence integrate_constant_is_box_volume gaussian_width_law)
  set(moist_test_selected_math_grid_3d_invariance single_atom asymmetric_trimer)
  set(moist_test_selected_math_grid_3d_threaded
    batched_ifft_uses_threaded_backend
    molecular_nufft_type12_threaded molecular_nufft_type3_threaded
    radial_trafo_one_per_thread partition_grid_threaded)

  set(moist_test_srcs
    "${_unit}/main.F90"
    "${_unit}/test_helpers.f90"
    "${_unit}/test_moz_fixtures.f90"
    "${_unit}/test_model_component_helper.f90"
  )
  foreach(suite IN LISTS moist_test_suites)
    list(APPEND moist_test_srcs "${_unit}/test_${suite}.f90")
  endforeach()

  add_executable(moist-tester ${moist_test_srcs})
  set_target_properties(moist-tester PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-modules"
    LINKER_LANGUAGE ${MOIST_EXE_LINKER_LANGUAGE}
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(moist-tester PRIVATE -ffree-line-length-none -fbacktrace)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(moist-tester PRIVATE ${MOIST_INTEL_FFLAGS})
  endif()
  target_link_libraries(moist-tester PRIVATE
    moist mstore::mstore test-drive::test-drive)

  if(MOIST_API)
    add_executable(moist-api-tester test/api/example.c)
    set_property(TARGET moist-api-tester PROPERTY C_STANDARD 11)
    target_link_libraries(moist-api-tester PRIVATE moist)
    add_test(NAME api-c COMMAND moist-api-tester)
    add_executable(moist-api-cpp-tester test/api/cpp.cpp)
    set_property(TARGET moist-api-cpp-tester PROPERTY CXX_STANDARD 11)
    target_link_libraries(moist-api-cpp-tester PRIVATE moist)
    add_test(NAME api-cpp COMMAND moist-api-cpp-tester)
  endif()

  foreach(suite IN LISTS moist_test_suites)
    if(DEFINED moist_test_selected_${suite})
      set(runs "")
      foreach(sel IN LISTS moist_test_selected_${suite})
        add_test(NAME ${suite}_${sel} COMMAND moist-tester ${suite} ${sel})
        list(APPEND runs ${suite}_${sel})
      endforeach()
    else()
      add_test(NAME ${suite} COMMAND moist-tester ${suite})
      set(runs ${suite})
    endif()
    set_tests_properties(${runs} PROPERTIES
      ENVIRONMENT "MOIST_SOURCE_ROOT=${CMAKE_CURRENT_SOURCE_DIR}"
      TIMEOUT 360
    )
  endforeach()
endif()

message(STATUS "moist ${MOIST_VERSION_FULL} (commit ${VCS_TAG})")
message(STATUS "  LAPACK backend : ${MOIST_LAPACK}")
message(STATUS "  OpenMP/API     : ${MOIST_OPENMP} / ${MOIST_API}")
