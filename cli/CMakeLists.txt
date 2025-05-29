cmake_minimum_required(VERSION 3.10) 
project(MeshMonkCLI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add FetchContent module if not already included (root CMake usually has it)
include(FetchContent)

# Declare cxxopts dependency via FetchContent
FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG v3.1.1 # Using a recent stable tag for cxxopts
)
FetchContent_MakeAvailable(cxxopts)

# Add the executable
add_executable(meshmonk_cli cli.cpp)

# Link against libraries
target_link_libraries(meshmonk_cli PRIVATE
    meshmonk_shared  # Link against the shared library from Phase 1
    cxxopts::cxxopts # Link against cxxopts (now available via FetchContent)
    OpenMeshCore   
    OpenMeshTools  
)

# Include directories for cxxopts should be handled automatically by FetchContent_MakeAvailable
# and linking against the cxxopts::cxxopts target.
