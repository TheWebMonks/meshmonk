# --- Collect Source Files ---
# Globs all .cpp files from the src directory. Consider listing them explicitly
# for more robust build scripts.
file(GLOB_RECURSE MESHMONK_SOURCES src/*.cpp)

# --- Build OpenMesh Statically ---
# Important: Set these BEFORE add_subdirectory for OpenMesh
set(BUILD_SHARED_LIBS
    OFF
    CACHE BOOL "Build OpenMesh as static" FORCE)
set(OPENMESH_BUILD_APPS
    OFF
    CACHE BOOL "Do not build OpenMesh apps" FORCE)
# Position-independent code required so static libs can be linked into the
# shared Python extension module (_meshmonk_core.so)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add OpenMesh as a subdirectory. CMake will find its CMakeLists.txt. The path
# is relative to this CMakeLists.txt file.
add_subdirectory(../vendor/OpenMesh-11.0.0
                 ${CMAKE_CURRENT_BINARY_DIR}/OpenMesh_build)

# --- Define the meshmonk_lib Library ---
add_library(meshmonk_lib STATIC ${MESHMONK_SOURCES})

# --- Include Directories ---
target_include_directories(
  meshmonk_lib
  PUBLIC # To find headers in library/include/ (e.g., meshmonk/meshmonk.hpp)
         ${CMAKE_CURRENT_SOURCE_DIR}/include
         # For internal library headers in library/src/ (needed by public
         # meshmonk.hpp)
         ${CMAKE_CURRENT_SOURCE_DIR}/src
  PRIVATE # For nanoflann.hpp located in ../vendor/
          ${CMAKE_CURRENT_SOURCE_DIR}/../vendor)

# --- Link Dependencies ---
# OpenMesh's CMake builds SHAREDANDSTATIC on non-Windows platforms regardless
# of BUILD_SHARED_LIBS, exposing `OpenMeshCore`/`OpenMeshTools` as shared and
# `OpenMeshCoreStatic`/`OpenMeshToolsStatic` as static. We link the static
# variants so the wheel has no runtime dependency on libOpenMeshTools.so. On
# Windows only the single `OpenMeshCore`/`OpenMeshTools` target exists (built
# STATIC since OpenMesh has no DLL exports).
if(WIN32)
  target_link_libraries(meshmonk_lib PRIVATE OpenMeshCore OpenMeshTools)
else()
  target_link_libraries(meshmonk_lib PRIVATE OpenMeshCoreStatic
                                             OpenMeshToolsStatic)
endif()

# Optional: Set properties if needed, e.g., C++ standard (already set at root,
# but can be target-specific) set_target_properties(meshmonk_shared PROPERTIES
# CXX_STANDARD ${CMAKE_CXX_STANDARD})
