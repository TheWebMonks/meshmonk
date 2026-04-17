# MeshMonk: High-Performance 3D Mesh Registration

MeshMonk is a C++ library designed for efficient 3D mesh registration. It provides functionalities for rigid, non-rigid, and pyramid-based registration approaches. This repository is currently undergoing modernization. Phase 1, focusing on a CMake build system for the core library, C++ examples, and a command-line interface, is complete.

## Overview

The core of MeshMonk is a C++ library that can be integrated into your own projects.
Currently, it can be used:

*   **As a C++ Library:** Link against `libmeshmonk_shared` in your C++ applications. A C++ example (`MeshMonkExample`) demonstrating library usage is buildable via CMake.
*   **As a Command-Line Tool:** Use `meshmonk_cli` to perform registrations directly.

Future phases will focus on:
*   Re-integrating MATLAB MEX bindings.
*   Introducing Python bindings.

This document primarily focuses on building the core MeshMonk C++ library, its C++ example, and the command-line interface using CMake.

## Repository Structure

The repository has been reorganized for clarity and modularity:

*   `CMakeLists.txt`: Root CMake file for the project.
*   `library/`: Contains the core MeshMonk library.
    *   `CMakeLists.txt`: CMake file for building `libmeshmonk_shared`.
    *   `include/meshmonk/`: Public API headers for the library (e.g., `meshmonk.hpp`).
    *   `src/`: Internal source code and private headers for the library.
    *   `examples/`: C++ example code demonstrating library usage.
        *   `CMakeLists.txt`: CMake file for building the C++ example(s).
*   `cli/`: Contains the source code for the command-line interface.
    *   `CMakeLists.txt`: CMake file for building `meshmonk_cli`.
    *   `cli.cpp`: Main source file for the CLI.
*   `vendor/`: Contains third-party dependencies.
    *   `OpenMesh-11.0.0/`: Source for the OpenMesh library, built statically.
    *   `nanoflann.hpp`: Header-only library for KD-tree search.
*   `demo/`: Contains example mesh data for use with the CLI.
*   `docs/`: Contains legacy documentation.

## Building MeshMonk (Core Library, C++ Example & CLI)

This method compiles the core shared library (`libmeshmonk_shared`), a C++ example (`MeshMonkExample`), and the command-line tool (`meshmonk_cli`).

### Prerequisites

Ensure you have the following dependencies installed on your system:

1.  **C++14 Compiler:** A C++ compiler supporting C++14 (e.g., GCC >= 5 or Clang >= 3.4).
    *   Ubuntu: `sudo apt update && sudo apt install build-essential g++`
2.  **CMake:** Version 3.10 or higher.
    *   Ubuntu: `sudo apt install cmake`
3.  **Eigen3 & cxxopts:** These are handled automatically by CMake using `FetchContent` and do not require separate system installations for the build process.

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
    rm -rf build # Clean previous build if it exists
    mkdir build
    cd build
    ```

3.  **Run CMake:**
    This command configures the project and generates the necessary build files. It will download Eigen3 (via FetchContent), cxxopts (via FetchContent), and prepare to compile the vendored OpenMesh library, the `meshmonk_shared` library, the C++ example, and the `meshmonk_cli` executable.
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
*   **Command-Line Executable:** `cli/meshmonk_cli`. (Path relative to `build/` directory)

## Using the `meshmonk_cli` Command-Line Tool

The `meshmonk_cli` tool allows you to perform registrations directly without needing to write custom C++ code.

### Prerequisites for Running
The CLI tool is built as part of the main CMake build process. Ensure all build prerequisites were met. `cxxopts` is handled automatically via FetchContent.

### Basic Syntax
After building, the executable is located in the `build/cli/` directory. Run it from the `build` directory:

```bash
./cli/meshmonk_cli <command> <source_mesh.obj> <target_mesh.obj> <output_mesh.obj> [options...]
```

### Available Commands

*   `pyramid_reg`: Performs pyramid-based non-rigid registration.
*   `rigid_reg`: Performs rigid registration.

You can see all available options for each command by running:
```bash
./cli/meshmonk_cli --help
# Or for a specific command, e.g.:
# ./cli/meshmonk_cli rigid_reg --help
```

### Example: Rigid Registration

Register `Template.obj` (source) to `demoFace.obj` (target) using rigid registration. The mesh data files are located in the `data/` directory at the root of the project.

1.  **Navigate to the build directory (if not already there):**
    ```bash
    # Assuming you are in the MeshMonk root directory:
    cd build
    ```

2.  **Run the `rigid_reg` command:**
    ```bash
    ./cli/meshmonk_cli rigid_reg ../data/Template.obj ../data/demoFace.obj ../demo/rigid_output.obj --transform_output ../demo/rigid_transform.txt
    ```
    *   `../data/Template.obj`: Path to the source mesh.
    *   `../data/demoFace.obj`: Path to the target mesh.
    *   `../demo/rigid_output.obj`: Filename for the transformed source mesh.
    *   `--transform_output ../demo/rigid_transform.txt`: Saves the 4x4 transformation matrix.

3.  **Expected Output (in `demo/` directory):**
    *   `rigid_output.obj`: The source mesh transformed to align with the target mesh.
    *   `rigid_transform.txt`: A text file containing the 4x4 transformation matrix.

### Example: Pyramid (Non-Rigid) Registration
This example uses the output from the rigid registration step (`rigid_output.obj`) as the source mesh.

```bash
./cli/meshmonk_cli pyramid_reg ../demo/rigid_output.obj ../data/demoFace.obj ../demo/pyramid_output.obj
```
*   This saves the non-rigidly transformed source mesh to `../demo/pyramid_output.obj`.
*   It uses default parameters; check `--help` for all options.

**Note on Demo Mesh Data:**
When running with the provided demo meshes (`Template.obj`, `demoFace.obj`), you might see console warnings from OpenMesh (e.g., "complex vertex", "complex edge", "patch re-linking failed", or missing `.mtl` files). These are generally related to the topology of the demo data itself or how OpenMesh handles them, and typically do not prevent `meshmonk_cli` from completing the registration process successfully.

## Using MeshMonk as a C++ Library

If you want to use MeshMonk's functionalities in your own C++ projects:

1.  **Build `libmeshmonk_shared`:** Follow the CMake instructions above.
2.  **Include Headers:**
    *   Your project should include the main MeshMonk API header: `#include "meshmonk/meshmonk.hpp"`.
    *   Ensure your build system is configured to find headers from MeshMonk's `library/include/` directory and `library/src/` directory (as public headers of `meshmonk_shared` include files from `library/src/`).
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

## Using MeshMonk from MATLAB

To use MeshMonk's functionalities within MATLAB, you first need to build the core shared library using CMake and then compile the MEX functions.

### 1. Build `libmeshmonk_shared` via CMake

Follow the main build instructions in this README to configure and build MeshMonk using CMake. This will produce the `libmeshmonk_shared` library (e.g., `libmeshmonk_shared.so` on Linux, `libmeshmonk_shared.dylib` on macOS, or `meshmonk_shared.dll` and `meshmonk_shared.lib` on Windows). This library is typically found in the `build/library/` directory (or `build/bin/` for the DLL on Windows, depending on CMake configuration).

### 2. Compile MEX Functions

Once `libmeshmonk_shared` is built, you can compile the MEX functions:

1.  Open MATLAB.
2.  Navigate MATLAB's current directory to the `matlab/` directory within the MeshMonk project (e.g., `cd path/to/meshmonk/matlab`).
3.  Run the appropriate MEX compilation script:
    *   On **Linux or macOS**: Execute `mex_all` in the MATLAB command window.
    *   On **Windows**: Execute `mex_windows_all` in the MATLAB command window.

    ```matlab
    % For Linux/macOS
    mex_all

    % For Windows
    % mex_windows_all
    ```

    These scripts will compile all necessary MEX functions and place them in the current directory (`matlab/`).

**Notes for MEX Compilation:**
*   The MEX scripts (`mex_all.m` and `mex_windows_all.m`) assume that the MeshMonk project root is one level above the `matlab/` directory (i.e., `meshmonkRoot = '..';`).
*   The scripts expect to find Eigen3 headers in a directory fetched by CMake (typically `build/_deps/eigen-src/`). If you have Eigen3 installed in a custom system location, you may need to modify the `eigenIncludeDir` variable at the beginning of the respective `mex_*.m` script.
*   Similarly, OpenMesh headers are expected to be found within the CMake build directory (`build/vendor/OpenMesh-11.0.0/_build/src/`).

### 3. MATLAB Runtime Environment Setup

For MATLAB to find `libmeshmonk_shared` and the compiled MEX functions at runtime:

*   **MEX Functions:** The demo scripts in the `demo/` directory are configured to find MEX files if they are located in the `matlab/` directory (due to `addpath(fullfile('..', 'matlab'))`). Ensure your compiled MEX files (`.mexa64`, `.mexw64`, etc.) are in `matlab/`.

*   **`libmeshmonk_shared` Library:**
    *   **Linux:** Add the directory containing `libmeshmonk_shared.so` (e.g., `path/to/meshmonk/build/library/`) to your `LD_LIBRARY_PATH` environment variable before launching MATLAB.
      ```bash
      export LD_LIBRARY_PATH=/path/to/meshmonk/build/library:$LD_LIBRARY_PATH
      matlab
      ```
      For a persistent setting, add this line to your `~/.bashrc` or `~/.zshrc`.

    *   **macOS:** Add the directory containing `libmeshmonk_shared.dylib` (e.g., `path/to/meshmonk/build/library/`) to your `DYLD_LIBRARY_PATH` environment variable.
      ```bash
      export DYLD_LIBRARY_PATH=/path/to/meshmonk/build/library:$DYLD_LIBRARY_PATH
      matlab
      ```
      Alternatively, you can create a `startup.m` file in your MATLAB user path (typically `~/Documents/MATLAB/`) and add:
      ```matlab
      addpath('path/to/meshmonk/build/library');
      % Or, more robustly for .dylib, consider a symbolic link or install_name_tool adjustments post-build.
      % However, DYLD_LIBRARY_PATH is the most common method for .dylib files not in standard locations.
      ```
      *(Note: `addpath` in `startup.m` is generally for M-files and MEX-files, not typically for `.dylib` dependencies directly, but MATLAB might pick them up if they are in its search path. `DYLD_LIBRARY_PATH` is more reliable for shared libraries.)*

    *   **Windows:** Ensure the directory containing `meshmonk_shared.dll` (e.g., `path/to/meshmonk/build/library/` or `path/to/meshmonk/build/bin/`) is in your system's `PATH` environment variable. Alternatively, you can add this directory to MATLAB's path using a `startup.m` file in your MATLAB user path (typically `Documents\MATLAB\`):
      ```matlab
      addpath('C:\path\to\meshmonk\build\library'); % Or wherever the .dll is
      ```

### 4. Running Demo Scripts

The demo scripts (e.g., `test_compute_correspondences.m`) are located in the `demo/` directory.
1.  Open MATLAB.
2.  Ensure your runtime environment is set up as described above.
3.  Navigate MATLAB's current directory to `path/to/meshmonk/demo/`.
4.  Open and run any of the `test_*.m` scripts. For example:
    ```matlab
    test_compute_correspondences
    ```

---
*MATLAB and Python support are planned for future development phases.*
