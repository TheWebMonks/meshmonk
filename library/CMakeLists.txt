# --- Collect Source Files ---
# Globs all .cpp files from the src directory.
# Consider listing them explicitly for more robust build scripts.
file(GLOB_RECURSE MESHMONK_SOURCES src/*.cpp)

# --- Build OpenMesh Statically ---
# Important: Set these BEFORE add_subdirectory for OpenMesh
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build OpenMesh as static" FORCE)
set(OPENMESH_BUILD_APPS OFF CACHE BOOL "Do not build OpenMesh apps" FORCE)

# Add OpenMesh as a subdirectory. CMake will find its CMakeLists.txt.
# The path is relative to this CMakeLists.txt file.
add_subdirectory(../vendor/OpenMesh-11.0.0 ${CMAKE_CURRENT_BINARY_DIR}/OpenMesh_build)

# --- Define the meshmonk_shared Library ---
add_library(meshmonk_shared SHARED ${MESHMONK_SOURCES})

# --- Include Directories ---
target_include_directories(meshmonk_shared
  PUBLIC
    # To find headers in library/include/ (e.g., meshmonk/meshmonk.hpp)
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    # Eigen include directories will be inherited from the Eigen3::Eigen target
    # For internal library headers in library/src/ (needed by public meshmonk.hpp)
    ${CMAKE_CURRENT_SOURCE_DIR}/src
  PRIVATE
    # For nanoflann.hpp located in ../vendor/
    ${CMAKE_CURRENT_SOURCE_DIR}/../vendor
)

# --- Link Dependencies ---
target_link_libraries(meshmonk_shared
  PUBLIC
    # Link against Eigen - PUBLIC because Eigen types are in public headers
    Eigen3::Eigen
  PRIVATE
    # Link against the static OpenMesh libraries built as part of this project
    OpenMeshCore
    OpenMeshTools
)

# Optional: Set properties if needed, e.g., C++ standard (already set at root, but can be target-specific)
# set_target_properties(meshmonk_shared PROPERTIES CXX_STANDARD ${CMAKE_CXX_STANDARD})