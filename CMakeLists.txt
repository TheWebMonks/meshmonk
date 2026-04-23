cmake_minimum_required(VERSION 3.26)
project(MeshMonk LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Dependencies ---
# Eigen (vendored headers — no FetchContent, no download at build time)
include_directories(${CMAKE_SOURCE_DIR}/vendor/eigen-3.4.0)

# nanoflann (header-only, vendor root)
include_directories(${CMAKE_SOURCE_DIR}/vendor)

# tl::expected (header-only, vendored)
include_directories(${CMAKE_SOURCE_DIR}/vendor/tl)

# --- Python + nanobind ---
# find_package(Python) must be called before find_package(nanobind). When
# building via scikit-build-core, CMAKE_MODULE_PATH is patched to use
# scikit-build-core's own FindPython.cmake that works without dev headers.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# Auto-discover nanobind cmake dir from the Python interpreter when it is not
# already on CMAKE_PREFIX_PATH. This lets `cmake -S . -B build` work after `uv
# pip install nanobind` without needing to pass -DCMAKE_PREFIX_PATH.
if(NOT DEFINED nanobind_DIR AND NOT nanobind_FOUND)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
            "import nanobind; print(nanobind.cmake_dir())"
    OUTPUT_VARIABLE NB_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE NB_CMAKE_RET)
  if(NB_CMAKE_RET EQUAL 0 AND NB_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${NB_CMAKE_DIR}")
  endif()
endif()
find_package(nanobind CONFIG REQUIRED)

# --- Subdirectories ---
option(MESHMONK_PROFILING "Enable registration profiling instrumentation" OFF)
add_subdirectory(library)

# --- nanobind extension module ---
nanobind_add_module(_meshmonk_core STABLE_ABI bindings/bindings.cpp)
# OpenMP: verified not linked. vci_openmp() defined but never invoked in
# OpenMesh CMake tree; USE_OPENMP not checked. No fix needed.
target_link_libraries(_meshmonk_core PRIVATE meshmonk_lib)
target_include_directories(
  _meshmonk_core
  PRIVATE ${CMAKE_SOURCE_DIR}/library/include ${CMAKE_SOURCE_DIR}/vendor
          ${CMAKE_SOURCE_DIR}/vendor/tl)
install(TARGETS _meshmonk_core LIBRARY DESTINATION meshmonk)

# --- Type stubs (PEP 561) ---
nanobind_add_stub(
  _meshmonk_core_stub
  MODULE
  _meshmonk_core
  OUTPUT
  meshmonk/_meshmonk_core.pyi
  PYTHON_PATH
  $<TARGET_FILE_DIR:_meshmonk_core>
  DEPENDS
  _meshmonk_core
  MARKER_FILE
  meshmonk/py.typed)
