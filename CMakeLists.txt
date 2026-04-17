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

# --- Subdirectories ---
add_subdirectory(library)
