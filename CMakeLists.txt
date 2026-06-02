# CMake build for the moist implicit-solvation library. This mirrors the
# canonical Meson build (meson.build + meson_options.txt + config/meson.build);
# Meson remains the reference build system. The CMake build embeds the vendored
# Fortran dependencies under subprojects/ via add_subdirectory and resolves
# BLAS/LAPACK/OpenMP through CMake's standard find modules.
#
# Quick start:
#   cmake -B build-cmake -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo
#   cmake --build build-cmake
# Produces the static library libmoist.a and the CLI executable `moist`.

cmake_minimum_required(VERSION 3.18)

# Read the single source of truth for the version (the VERSION file).
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" MOIST_VERSION_FULL)
string(STRIP "${MOIST_VERSION_FULL}" MOIST_VERSION_FULL)
# The VERSION file may carry a semver pre-release suffix (e.g. 0.6.0-alpha.1).
# Meson keeps the full string in meson.project_version() but derives the numeric
# library version with split('-')[0]; CMake's project(VERSION) is stricter still
# and accepts only numeric major[.minor[.patch[.tweak]]]. Strip from the first
# '-' for project(), and keep the full string in MOIST_VERSION_FULL for display.
string(REGEX REPLACE "-.*$" "" MOIST_VERSION "${MOIST_VERSION_FULL}")

project(
  moist
  VERSION "${MOIST_VERSION}"
  DESCRIPTION "Modular and Open-source Implicit Solvation Toolkit"
  LANGUAGES Fortran C
)

# ---------------------------------------------------------------------------
# Options (mirrors meson_options.txt)
# ---------------------------------------------------------------------------
set(MOIST_LAPACK "auto" CACHE STRING "Linear algebra backend")
# Mirrors the `lapack` combo in meson_options.txt.
set_property(CACHE MOIST_LAPACK PROPERTY STRINGS
  auto mkl mkl-rt openblas netlib accelerate custom)
# Mirrors meson's `custom_libraries`: the libraries to link for the custom
# backend -- absolute paths, plain names resolved by the linker, or -L/-l flags
# passed through as list items.
set(MOIST_CUSTOM_LIBRARIES "" CACHE STRING
  "Libraries for the custom LAPACK backend (mirror of meson custom_libraries)")
option(MOIST_OPENMP "Use OpenMP parallelisation"            ON)
option(MOIST_API    "Build C API using iso_c_binding"       ON)
option(MOIST_RISM   "Build RISM module (auto-enables FFTW and HDF5)" OFF)
option(MOIST_FFTW   "Enable FFTW support"                   OFF)
option(MOIST_HDF5   "Enable HDF5 support for I/O"           OFF)
option(MOIST_ILP64  "Enable BLAS/LAPACK ilp64 support"      OFF)
option(MOIST_TESTS  "Build the moist unit test suite (ctest)" ON)

# RISM needs FFTW for the reciprocal-space transforms and HDF5 for its VV
# caches, so enabling it forces both on regardless of their own switches --
# mirrors the use_fftw/use_hdf5 implications in config/meson.build.
if(MOIST_RISM)
  set(MOIST_FFTW ON)
  set(MOIST_HDF5 ON)
endif()

# Default to an optimized-with-debug build, matching meson's debugoptimized.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Third-party dependencies (BLAS / LAPACK / OpenMP)
#
# Resolved BEFORE the vendored subprojects so that the LAPACK::LAPACK /
# BLAS::BLAS / OpenMP::OpenMP_Fortran imported targets already exist; the
# subprojects guard their own lookups with `if(NOT TARGET ...)` and therefore
# reuse ours instead of falling back to a network FetchContent.
# ---------------------------------------------------------------------------
# The `custom` backend mirrors meson: link exactly the libraries named in
# MOIST_CUSTOM_LIBRARIES and synthesise the LAPACK::LAPACK / BLAS::BLAS imported
# targets the rest of this file (and the vendored subprojects) link against, so
# no find_package lookup runs. Every other backend maps onto a CMake BLA_VENDOR
# and resolves through FindBLAS/FindLAPACK.
if(MOIST_LAPACK STREQUAL "custom")
  if(NOT MOIST_CUSTOM_LIBRARIES)
    message(FATAL_ERROR
      "MOIST_LAPACK=custom requires -DMOIST_CUSTOM_LIBRARIES=... "
      "(mirror of meson's custom_libraries)")
  endif()
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${MOIST_CUSTOM_LIBRARIES}")
  # A single library set covers both, so LAPACK::LAPACK just chains BLAS::BLAS.
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES
    INTERFACE_LINK_LIBRARIES "BLAS::BLAS")
else()
  if(MOIST_LAPACK STREQUAL "mkl")
    set(BLA_VENDOR Intel10_64lp)        # MKL, multi-library interface
  elseif(MOIST_LAPACK STREQUAL "mkl-rt")
    set(BLA_VENDOR Intel10_64_dyn)      # MKL single dynamic library (mkl_rt)
  elseif(MOIST_LAPACK STREQUAL "openblas")
    set(BLA_VENDOR OpenBLAS)
  elseif(MOIST_LAPACK STREQUAL "netlib")
    set(BLA_VENDOR Generic)
  elseif(MOIST_LAPACK STREQUAL "accelerate")
    # Apple's Accelerate framework ships BLAS/LAPACK with macOS and exposes the
    # legacy LP64 interface only -- guard both, mirroring config/meson.build.
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
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# ---------------------------------------------------------------------------
# Vendored Fortran subprojects
#
# These ship CMakeLists that define a compiled target `<name>-lib` plus an
# INTERFACE facade `<name>`, but NOT the namespaced `<name>::<name>` alias that
# their dependents probe for (that alias is only created on install(EXPORT)).
# We add each in topological order and synthesise the `::` alias ourselves so
# every `if(NOT TARGET <name>::<name>)` guard short-circuits -- keeping the
# whole graph offline. EXCLUDE_FROM_ALL stops their unconditional app/ and
# test/ subdirectories from compiling (they are still configured, so referenced
# targets such as mstore::mstore must exist -- hence mstore is added too).
# ---------------------------------------------------------------------------
# Silence the subprojects' own test/app suites.
set(BUILD_TESTING            OFF CACHE BOOL "" FORCE)  # toml-f, fclap, mstore, ...
set(WITH_TESTS               OFF CACHE BOOL "" FORCE)  # multicharge, s-dftd3
set(ENABLE_TESTING           OFF CACHE BOOL "" FORCE)  # jonquil
set(TESTDRIVE_BUILD_TESTING  OFF CACHE BOOL "" FORCE)  # test-drive

# mctc-lib and dftd4 call enable_testing()/add_test() with no BUILD_TESTING
# guard, registering ctest cases for their own *-tester binaries -- which are
# EXCLUDE_FROM_ALL here and never built, so the cases would dangle as "Not Run".
# Shadow add_test() with a flag-gated no-op (the builtin stays reachable as
# _add_test) so only moist's own suite, which flips _moist_allow_tests on,
# registers. Safe because neither subproject calls set_tests_properties.
set(_moist_allow_tests OFF)
function(add_test)
  if(_moist_allow_tests)
    _add_test(${ARGV})
  endif()
endfunction()

set(_moist_subprojects "${CMAKE_CURRENT_SOURCE_DIR}/subprojects")

# Helper: add a vendored subproject and expose its `<name>::<name>` alias.
function(moist_add_subproject name)
  if(NOT TARGET "${name}")
    add_subdirectory("${_moist_subprojects}/${name}" EXCLUDE_FROM_ALL)
  endif()
  if(NOT TARGET "${name}::${name}")
    add_library("${name}::${name}" ALIAS "${name}")
  endif()
endfunction()

#                        depends on
moist_add_subproject(test-drive)   # (leaf)
moist_add_subproject(toml-f)       # (leaf)
moist_add_subproject(fclap)        # (leaf, CLI only)
moist_add_subproject(jonquil)      # toml-f
moist_add_subproject(mctc-lib)     # jonquil, toml-f (WITH_JSON)
moist_add_subproject(mstore)       # mctc-lib  (only to satisfy dftd4 test config)
moist_add_subproject(multicharge)  # mctc-lib, LAPACK
moist_add_subproject(dftd4)        # mctc-lib, multicharge, BLAS, mstore

# ---------------------------------------------------------------------------
# moist source list
#
# Globbed to track the ~35 meson.build fragments without hand-mirroring them;
# CMake derives Fortran module compile order automatically, so listing order is
# irrelevant. CONFIGURE_DEPENDS re-globs on rebuild. Conditional source sets
# follow the same gating as the meson build.
# ---------------------------------------------------------------------------
file(GLOB_RECURSE moist_srcs CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.F90"
)

# The committed build_info.f90 is the fpm fallback ("unknown"); we regenerate it
# from build_info.f90.in with the live git commit, so drop the committed copy to
# avoid a duplicate moist_build_info module.
list(REMOVE_ITEM moist_srcs
  "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/build_info.f90")

# api.f90 only when the C API (or Python extension) is requested.
if(NOT MOIST_API)
  list(REMOVE_ITEM moist_srcs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/api.f90")
endif()

# fftw3.F90 only with FFTW (RISM forces MOIST_FFTW on above).
if(NOT MOIST_FFTW)
  list(REMOVE_ITEM moist_srcs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/math/fftw3.F90")
endif()

# The entire RISM subtree (and its two model wrappers) only with RISM.
if(NOT MOIST_RISM)
  list(FILTER moist_srcs EXCLUDE REGEX "/src/moist/model/rism/")
  list(REMOVE_ITEM moist_srcs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/model/rism1d.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/model/rism3d.f90")
endif()

# Build-info: stamp the short git commit (fallback "unknown" outside a checkout).
find_package(Git QUIET)
set(VCS_TAG "unknown")
if(Git_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE VCS_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(NOT VCS_TAG)
    set(VCS_TAG "unknown")
  endif()
endif()
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/moist/build_info.f90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.f90"
  @ONLY
)
list(APPEND moist_srcs "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.f90")

# ---------------------------------------------------------------------------
# moist library
# ---------------------------------------------------------------------------
add_library(moist ${moist_srcs})
set_target_properties(moist PROPERTIES
  VERSION "${PROJECT_VERSION}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_include_directories(moist PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)

# ILP64 integer-kind macros are consumed by the preprocessed .F90 sources and
# are unconditional in the meson build.
if(MOIST_ILP64)
  target_compile_definitions(moist PRIVATE IK=i8 MOIST_ILP64=1)
else()
  target_compile_definitions(moist PRIVATE IK=i4 MOIST_ILP64=0)
endif()
if(MOIST_FFTW)
  target_compile_definitions(moist PRIVATE WITH_FFTW=1)
endif()
if(MOIST_RISM)
  target_compile_definitions(moist PRIVATE WITH_RISM=1)
endif()
if(MOIST_HDF5)
  target_compile_definitions(moist PRIVATE WITH_HDF5=1)
endif()

# Compiler-specific flags (mirrors config/meson.build + moist_fortran_args).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(moist PRIVATE
    -ffree-line-length-none -fbacktrace
    -fopenmp-simd --param max-inline-insns-auto=300
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # -fp-model=precise: ifx defaults to a fast FP model that reassociates
  # reductions / contracts FMAs, diverging from gfortran and rounding the
  # analytic-vs-reference tests over tolerance. Mirrors config/meson.build.
  target_compile_options(moist PRIVATE -traceback -fp-model=precise)
endif()

# Link the dependency stack (matches lib_deps in config/meson.build).
target_link_libraries(moist PUBLIC
  mctc-lib::mctc-lib
  multicharge::multicharge
  dftd4::dftd4
  jonquil::jonquil
  test-drive::test-drive
  LAPACK::LAPACK
  BLAS::BLAS
)
if(MOIST_OPENMP)
  target_link_libraries(moist PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(MOIST_FFTW)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
  target_link_libraries(moist PUBLIC PkgConfig::FFTW3)
endif()

if(MOIST_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS Fortran)
  target_link_libraries(moist PUBLIC HDF5::HDF5)
endif()

# ---------------------------------------------------------------------------
# CLI application (app/) -> executable `moist`
# ---------------------------------------------------------------------------
add_executable(moist-exe
  app/main.f90
  app/cli.f90
  app/driver.f90
)
set_target_properties(moist-exe PROPERTIES
  OUTPUT_NAME moist
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/app-modules"
  LINKER_LANGUAGE Fortran
)
# The CLI sources are .f90 but use #ifdef WITH_RISM, so force preprocessing.
target_compile_options(moist-exe PRIVATE -cpp)
if(MOIST_FFTW)
  target_compile_definitions(moist-exe PRIVATE WITH_FFTW=1)
endif()
if(MOIST_RISM)
  target_compile_definitions(moist-exe PRIVATE WITH_RISM=1)
endif()
target_link_libraries(moist-exe PRIVATE moist fclap)

# ---------------------------------------------------------------------------
# Install rules
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS moist moist-exe
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
if(MOIST_API)
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/moist.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/moist"
  FILES_MATCHING PATTERN "*.mod")

# ---------------------------------------------------------------------------
# Unit test suite (test/unit/) -> ctest
#
# Mirrors test/meson.build + test/unit/meson.build: a single `tester` binary is
# dispatched once per suite (the suite name is argv[1]), linked against moist
# and the mstore molecular fixtures. Each suite is registered as its own ctest
# case with the MOIST_SOURCE_ROOT env var the tests use to locate reference
# data. Run with `ctest --test-dir <build> --output-on-failure`.
# ---------------------------------------------------------------------------
if(MOIST_TESTS)
  enable_testing()
  set(_moist_allow_tests ON)  # re-enable the shadowed add_test() for moist's suite
  set(_unit "${CMAKE_CURRENT_SOURCE_DIR}/test/unit")

  # Kept in sync with the `tests` list in test/unit/meson.build (model_gems and
  # component_cds are intentionally omitted there as well).
  set(moist_test_suites
    utils radii
    math_linalg math_adjacency_list math_cell_grid math_solvers
    math_sorters math_trig math_grid math_numdiff
    cavity_iswig cavity_numsa
    component_pcm_cpcm
    cavity_drop_primitives cavity_drop_cfc_kernel cavity_drop_lsf
    cavity_drop_gradients cavity_drop_cpcm cavity_drop_integration_ref
  )
  if(MOIST_RISM)
    list(APPEND moist_test_suites rism_1d rism_thermo)
  endif()
  if(MOIST_HDF5)
    list(APPEND moist_test_suites utils_hdf5)
  endif()

  set(moist_test_srcs "${_unit}/main.F90" "${_unit}/test_helpers.f90")
  foreach(suite IN LISTS moist_test_suites)
    list(APPEND moist_test_srcs "${_unit}/test_${suite}.f90")
  endforeach()

  add_executable(moist-tester ${moist_test_srcs})
  set_target_properties(moist-tester PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-modules"
    LINKER_LANGUAGE Fortran
  )
  # main.F90 selects suites with the same #ifdefs as the library build.
  if(MOIST_FFTW)
    target_compile_definitions(moist-tester PRIVATE WITH_FFTW=1)
  endif()
  if(MOIST_RISM)
    target_compile_definitions(moist-tester PRIVATE WITH_RISM=1)
  endif()
  if(MOIST_HDF5)
    target_compile_definitions(moist-tester PRIVATE WITH_HDF5=1)
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(moist-tester PRIVATE -ffree-line-length-none -fbacktrace)
  endif()
  target_link_libraries(moist-tester PRIVATE
    moist mstore::mstore test-drive::test-drive)

  foreach(suite IN LISTS moist_test_suites)
    add_test(NAME ${suite} COMMAND moist-tester ${suite})
    set_tests_properties(${suite} PROPERTIES
      ENVIRONMENT "MOIST_SOURCE_ROOT=${CMAKE_CURRENT_SOURCE_DIR}"
      TIMEOUT 360
    )
  endforeach()
endif()

message(STATUS "moist ${MOIST_VERSION_FULL} (commit ${VCS_TAG})")
message(STATUS "  LAPACK backend : ${MOIST_LAPACK}")
message(STATUS "  OpenMP/API     : ${MOIST_OPENMP} / ${MOIST_API}")
message(STATUS "  RISM/FFTW/HDF5 : ${MOIST_RISM} / ${MOIST_FFTW} / ${MOIST_HDF5}")
