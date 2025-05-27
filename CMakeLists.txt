cmake_minimum_required(VERSION 3.16)
project(MeshMonk LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- OpenMesh ---
# Set options for OpenMesh before adding the subdirectory
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build OpenMesh as static library")
set(OPENMESH_BUILD_APPS OFF CACHE BOOL "Do not build OpenMesh apps")
# Add OpenMesh as a subdirectory
# Make sure the path vendor/OpenMesh-11.0.0 is correct
add_subdirectory(vendor/OpenMesh-11.0.0)

# --- Eigen3 ---
# Find Eigen3. Requires libeigen3-dev to be installed (e.g., sudo apt-get install libeigen3-dev)
# Alternatively, Eigen could be added as a vendored library.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# --- nanoflann ---
# nanoflann is header-only, so we just need to ensure its path is in the include directories.

# --- MeshMonk Library (meshmonk_shared) ---
# Gather all .cpp files from the src directory
file(GLOB MESHMOMNK_SOURCES src/*.cpp)

add_library(meshmonk_shared SHARED
    meshmonk.cpp
    ${MESHMOMNK_SOURCES}
)

target_include_directories(meshmonk_shared PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}         # For meshmonk.hpp, global.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor # For nanoflann.hpp
    ${EIGEN3_INCLUDE_DIR}              # For Eigen headers
    # OpenMesh include directories are typically handled by linking with OpenMesh::Core
)

target_link_libraries(meshmonk_shared PUBLIC
    OpenMeshCore
    OpenMeshTools
    Eigen3::Eigen
)

# --- CLI Subdirectory ---
# This will be processed after the meshmonk_shared library is defined.
add_subdirectory(cli)
