cmake_minimum_required(VERSION 3.16)
project(MeshMonkCLI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find cxxopts.
# This might require the user to install libcxxopts-dev
# (e.g., sudo apt-get install libcxxopts-dev) or install cxxopts from source.
find_package(cxxopts REQUIRED CONFIG)

add_executable(meshmonk_cli cli.cpp)

target_link_libraries(meshmonk_cli PRIVATE
    meshmonk_shared  # Defined in the parent CMakeLists.txt
    cxxopts::cxxopts
)

# Include directories from meshmonk_shared and its dependencies (OpenMesh, Eigen)
# should be handled transitively by CMake if meshmonk_shared sets them as PUBLIC.
# If not, we might need to add them explicitly here.
# For example:
# target_include_directories(meshmonk_cli PRIVATE
#    ${CMAKE_SOURCE_DIR} # Access to meshmonk.hpp, global.hpp, src, vendor
# )
# However, meshmonk_shared's includes were set to PUBLIC, so this should be okay.
