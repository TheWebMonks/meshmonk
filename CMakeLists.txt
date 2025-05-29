cmake_minimum_required(VERSION 3.10)
project(MeshMonk LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Dependencies ---
# Eigen (via FetchContent)
include(FetchContent)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0 # Or a specific version/commit suitable for the project
)
FetchContent_MakeAvailable(Eigen3)
# Eigen3::Eigen is the target to link against

# nanoflann (header-only, ensure it's findable)
# Add vendor to include path, primarily for nanoflann.hpp
include_directories(vendor)

# --- Subdirectories ---
add_subdirectory(library)
