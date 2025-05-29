# MeshMonk: High-Performance 3D Mesh Registration

MeshMonk is a C++ library designed for efficient 3D mesh registration. It provides functionalities for rigid, non-rigid, and pyramid-based registration approaches. This repository is currently undergoing modernization. Phase 1, focusing on a CMake build system for the core library and C++ examples, is complete.

## Overview

The core of MeshMonk is a C++ library that can be integrated into your own projects.
Currently, it can be used:

*   **As a C++ Library:** Link against `libmeshmonk_shared` in your C++ applications. The first C++ example is buildable via CMake.

Future phases will focus on:
*   Re-integrating the Command-Line Interface (CLI).
*   Re-integrating MATLAB MEX bindings.
*   Introducing Python bindings.

This document primarily focuses on building the core MeshMonk C++ library and its C++ example using CMake.

## Repository Structure

The repository has been reorganized for clarity and modularity:

*   `CMakeLists.txt`: Root CMake file for the project.
*   `library/`: Contains the core MeshMonk library.
    *   `CMakeLists.txt`: CMake file for building `libmeshmonk_shared`.
    *   `include/meshmonk/`: Public API headers for the library (e.g., `meshmonk.hpp`).
    *   `src/`: Internal source code and private headers for the library.
    *   `examples/`: C++ example code demonstrating library usage.
        *   `CMakeLists.txt`: CMake file for building the C++ example(s).
*   `vendor/`: Contains third-party dependencies.
    *   `OpenMesh-11.0.0/`: Source for the OpenMesh library, built statically.
    *   `nanoflann.hpp`: Header-only library for KD-tree search.
*   `docs/`: Contains legacy documentation.

## Building MeshMonk (Core Library & C++ Example - Phase 1)

This method compiles the core shared library (`libmeshmonk_shared`) and a C++ example (`MeshMonkExample`).

### Prerequisites

Ensure you have the following dependencies installed on your system:

1.  **C++14 Compiler:** A C++ compiler supporting C++14 (e.g., GCC >= 5 or Clang >= 3.4).
    *   Ubuntu: `sudo apt update && sudo apt install build-essential g++`
2.  **CMake:** Version 3.10 or higher.
    *   Ubuntu: `sudo apt install cmake`
3.  **Eigen3:** A library for linear algebra. This is handled automatically by CMake using `FetchContent` and does not require a separate system installation for the build process.

### Compilation Steps

1.  **Clone the Repository (if you haven't already):**
    ```bash
    # Example using HTTPS:
    git clone https://github.com/TheWebMonks/meshmonk.git # Adjust if using SSH or a fork
    cd meshmonk
    ```

2.  **Create a Build Directory:**
    It's standard practice to build outside the source directory.
    ```bash
    mkdir build
    cd build
    ```

3.  **Run CMake:**
    This command configures the project and generates the necessary build files. It will download Eigen3 and prepare to compile the vendored OpenMesh library, the `meshmonk_shared` library, and the C++ example.
    ```bash
    cmake ..
    ```
    If you encounter errors, please ensure your CMake and C++ compiler are correctly installed and meet the version requirements.

4.  **Compile the Project:**
    This will build all targets. Using `-j$(nproc)` utilizes all available CPU cores for a faster build on Linux.
    ```bash
    make -j$(nproc)
    # Alternatively, you can use:
    # cmake --build . -- -j$(nproc)
    ```

### Build Outputs

Upon successful compilation, you will find the key outputs in your `build/` directory (or its subdirectories):

*   **Core Shared Library:** `library/libmeshmonk_shared.so` (on Linux), `library/libmeshmonk_shared.dylib` (on macOS), or `library/meshmonk_shared.dll` (on Windows). (Path relative to `build/` directory)
*   **C++ Example Executable:** `library/examples/MeshMonkExample`. (Path relative to `build/` directory)

## Using MeshMonk as a C++ Library

If you want to use MeshMonk's functionalities in your own C++ projects:

1.  **Build `libmeshmonk_shared`:** Follow the CMake instructions above.
2.  **Include Headers:**
    *   Your project should include the main MeshMonk API header: `#include "meshmonk/meshmonk.hpp"`.
    *   Ensure your build system is configured to find headers from MeshMonk's `library/include/` directory.
    *   Eigen3 headers will be available if you link against `meshmonk_shared` (as its include directories are made public) or if you also use `FetchContent` for Eigen in your project.
    *   `nanoflann.hpp` can be included via `#include <nanoflann.hpp>` if your project adds MeshMonk's `vendor/` directory to its include paths.
3.  **Link against the Library:**
    *   Link your project against `libmeshmonk_shared`. When using CMake, you can link against the `meshmonk_shared` target if MeshMonk is included as a subproject, or find and link the compiled library file.
    *   OpenMesh is statically linked into `libmeshmonk_shared`, so you typically won't need to link against OpenMesh libraries directly.

## Documentation

*   **Legacy Build Instructions:** For older, manual build instructions (not using the primary CMake system described above), you can refer to files in the `docs/` directory. Please note these are largely outdated for the current build process.
    *   [Ubuntu (Legacy)](docs/ubuntu.md)
    *   [OSX (Legacy)](docs/osx.md)
    *   [Windows (Legacy)](docs/windows.md)

---
*CLI and MATLAB support are planned for future development phases.*
