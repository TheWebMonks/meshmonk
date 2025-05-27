# MeshMonk: High-Performance 3D Mesh Registration

MeshMonk is a C++ library designed for efficient 3D mesh registration. It provides functionalities for rigid, non-rigid, and pyramid-based registration approaches, along with various helper modules for tasks like correspondence finding, inlier detection, and mesh downsampling.

## Overview

The core of MeshMonk is a C++ library that can be integrated into your own projects. It also offers several ways to be utilized:

*   **As a C++ Library:** Link against `libmeshmonk_shared` in your C++ applications.
*   **Via the `meshmonk_cli` Command-Line Tool:** A convenient tool for performing registrations directly from your terminal.
*   **From MATLAB:** Utilize MeshMonk's capabilities within the MATLAB environment through MEX bindings.
*   **(Future) From Python:** Python bindings are planned for future releases.

This document primarily focuses on building MeshMonk using CMake and then details how to use the command-line tool and integrate with MATLAB.

## Building MeshMonk (Core Library & CLI)

This is the primary method for building MeshMonk, which will compile the core shared library (`libmeshmonk_shared`) and the command-line tool (`meshmonk_cli`).

### Prerequisites

Ensure you have the following dependencies installed on your system:

1.  **C++17 Compiler:** A modern C++ compiler (e.g., GCC >= 7 or Clang >= 5).
    *   Ubuntu: `sudo apt update && sudo apt install build-essential g++`
2.  **CMake:** Version 3.16 or higher.
    *   Ubuntu: `sudo apt install cmake`
3.  **Eigen3:** A library for linear algebra. MeshMonk requires version 3.3 or higher.
    *   Ubuntu: `sudo apt install libeigen3-dev`
4.  **cxxopts (for `meshmonk_cli` only):** A lightweight C++ option parser library.
    *   Ubuntu: `sudo apt install libcxxopts-dev`
    *   *Note:* If `libcxxopts-dev` is not available through your package manager or you need a specific version, you might need to install it from source: [https://github.com/jarro2783/cxxopts](https://github.com/jarro2783/cxxopts)

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
    This command configures the project and generates the necessary build files. It will also prepare to compile the vendored OpenMesh library (version 11.0.0, included in `vendor/`), the `meshmonk_shared` library, and the `meshmonk_cli` executable.
    ```bash
    cmake ..
    ```
    If you encounter errors related to missing dependencies (Eigen3, cxxopts), please ensure they are correctly installed and discoverable by CMake.

4.  **Compile the Project:**
    This will build all targets. Using `-j$(nproc)` utilizes all available CPU cores for a faster build on Linux.
    ```bash
    make -j$(nproc)
    # Alternatively, you can use:
    # cmake --build . -- -j$(nproc)
    ```

### Build Outputs

Upon successful compilation, you will find the key outputs in your `build/` directory (or its subdirectories):

*   **Core Shared Library:** `libmeshmonk_shared.so` (on Linux), `libmeshmonk_shared.dylib` (on macOS), or `meshmonk_shared.dll` (on Windows). This library is typically located directly in the `build/` directory.
*   **Command-Line Tool:** `meshmonk_cli` executable, located in `build/cli/`.

## Using the `meshmonk_cli` Command-Line Tool

The `meshmonk_cli` tool allows you to perform registrations directly without needing to write custom C++ code or use MATLAB.

### Prerequisites for Running
The prerequisites are the same as for building, as the tool relies on the compiled library and its dependencies. Ensure they are met, or that the necessary runtime libraries are accessible.

### Basic Syntax

```bash
./cli/meshmonk_cli <command> <source_mesh.obj> <target_mesh.obj> <output_mesh.obj> [options...]
```
*(Run from the `build` directory)*

### Available Commands

*   `pyramid_reg`: Performs pyramid-based non-rigid registration.
*   `rigid_reg`: Performs rigid registration.

You can see all available options for each command by running:
```bash
./cli/meshmonk_cli --help
# Or for a specific command:
# ./cli/meshmonk_cli rigid_reg --help
```

### Example: Rigid Registration

Register `Template.obj` (source) to `demoFace.obj` (target) using rigid registration and save the transformation matrix.

1.  **Navigate to the build directory (if not already there):**
    ```bash
    # Assuming you are in the MeshMonk root directory:
    cd build
    ```

2.  **Run the `rigid_reg` command:**
    ```bash
    ./cli/meshmonk_cli rigid_reg ../demo/Template.obj ../demo/demoFace.obj ../demo/rigid_output.obj --transform_output ../demo/rigid_transform.txt
    ```
    *   `../demo/Template.obj`: Path to the source mesh (relative to `build` directory).
    *   `../demo/demoFace.obj`: Path to the target mesh (relative to `build` directory).
    *   `../demo/rigid_output.obj`: Filename for the transformed source mesh that will be saved.
    *   `--transform_output ../demo/rigid_transform.txt`: Specifies that the 4x4 transformation matrix should be saved to `../demo/rigid_transform.txt`.

3.  **Expected Output (in `demo/` directory):**
    *   `rigid_output.obj`: The source mesh transformed to align with the target mesh.
    *   `rigid_transform.txt`: A text file containing the 4x4 transformation matrix (numerical values may vary slightly):
        ```
         0.951203 -0.171045 -0.256844   2.13937
         0.115929  0.969435 -0.216259  -6.74492
         0.285982  0.175929  0.941953   58.7142
                0         0         0         1
        ```

### Example: Pyramid (Non-Rigid) Registration
Note that you first want to run the previous example, as we'll use its output (the rigidly transformed template) as mesh that we'll register to the target face.

```bash
./cli/meshmonk_cli pyramid_reg ../demo/rigid_output.obj ../demo/demoFace.obj ../demo/pyramid_output.obj
```
*   This saves the non-rigidly transformed source mesh to `../demo/pyramid_output.obj`.
*   It uses default parameters; check `--help` for all options.

**Note on Demo Mesh Data:**
When running with the provided demo meshes (`Template.obj`, `demoFace.obj`), you might see console warnings (e.g., missing `.mtl` files, OpenMesh errors like `complex vertex`). These are related to the demo data itself and typically don't prevent `meshmonk_cli` from completing registration.

## Using MeshMonk from MATLAB

To use MeshMonk within MATLAB, you first need to build its core C++ shared library (`libmeshmonk_shared`) using CMake, and then compile the MEX interface files.

### 1. Build the Core Library
Follow the instructions in the **"Building MeshMonk (Core Library & CLI)"** section above to compile `libmeshmonk_shared` using CMake. This step is crucial as it also builds OpenMesh statically into `libmeshmonk_shared`.

### 2. MEX Compilation

The MATLAB scripts `matlab/mex_all.m` (for Linux/macOS) and `matlab/mex_windows_all.m` (for Windows) are used to compile the MEX functions. You will likely need to **modify these scripts** to correctly point to headers and libraries from the CMake build.

**Key considerations for `mex` command flags:**

*   **Include Paths (`-I`):**
    *   MeshMonk headers:
        *   `-I<path_to_meshmonk_root>` (for `meshmonk.hpp`, `global.hpp`)
        *   `-I<path_to_meshmonk_root>/src`
        *   `-I<path_to_meshmonk_root>/vendor` (for `nanoflann.hpp`)
    *   Eigen3 headers:
        *   These are typically found from your system installation (e.g., `/usr/include/eigen3` on Linux if installed via `libeigen3-dev`). Add the correct path.
    *   OpenMesh headers:
        *   Since OpenMesh is built as part of the CMake process from the `vendor/` directory, the headers will be in a location like:
            `<path_to_meshmonk_root>/build/vendor/OpenMesh-11.0.0/_build/src/` (the exact path might vary slightly depending on CMake version and OpenMesh's internal build structure. You'll need to locate the `OpenMesh/Core/` and `OpenMesh/Tools/` directories).
            *You might need to add multiple `-I` flags for different OpenMesh subdirectories if necessary.*
*   **Library Paths (`-L`):**
    *   Point to the directory where `libmeshmonk_shared.[so|dylib|dll]` was created by CMake (typically `<path_to_meshmonk_root>/build/`).
*   **Libraries to Link (`-l`):**
    *   Link against `meshmonk_shared` (e.g., `-lmeshmonk_shared`).
    *   You should *not* need to explicitly link `-lOpenMeshCore` or `-lOpenMeshTools` here if they are statically linked into `libmeshmonk_shared`, which is the default setup in the provided `CMakeLists.txt`.
    *   Depending on your system and C++ standard library, you might need `-lstdc++`.

**Example (Conceptual) `mex` command line modification for a file like `pyramid_registration_mex.cpp` on Linux:**
```matlab
% In mex_all.m (adjust paths as per your system and MeshMonk location)
meshmonkRoot = '/path/to/your/meshmonk'; % Set this correctly
eigenInclude = '/usr/include/eigen3'; % Adjust if different
openMeshBuildInclude = [meshmonkRoot, '/build/vendor/OpenMesh-11.0.0/_build/src']; % Check this path

mex(['-I', meshmonkRoot], ...
    ['-I', meshmonkRoot, '/src'], ...
    ['-I', meshmonkRoot, '/vendor'], ...
    ['-I', eigenInclude], ...
    ['-I', openMeshBuildInclude], ... % May need more specific paths under here
    ['-L', meshmonkRoot, '/build'], ...
    '-lmeshmonk_shared', ...
    '-lstdc++', ... % May or may not be needed
    '-output', 'pyramid_registration_mex', ... % Ensure output name matches
    [meshmonkRoot, '/matlab/mexfiles/pyramid_registration_mex.cpp']);
```
You will need to adapt this for each MEX file and for your specific OS (especially library extensions and linker flags).

### 3. MATLAB Runtime Environment Setup

#### Linux
To ensure MATLAB can find `libmeshmonk_shared.so` and its dependencies at runtime, you might need to preload libraries when starting MATLAB. The `LD_PRELOAD` path should point to where `libmeshmonk_shared.so` is located (e.g., your `build` directory).
```bash
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6:<path_to_meshmonk_root>/build/libmeshmonk_shared.so matlab
```
Adjust `<path_to_meshmonk_root>` accordingly.

#### macOS
Add the path to `libmeshmonk_shared.dylib` (e.g., `<path_to_meshmonk_root>/build/`) to your `DYLD_LIBRARY_PATH` or use `install_name_tool` to make library paths absolute if you encounter issues. A common approach is to modify your `startup.m`:
```matlab
% In your startup.m
% Add path to the directory containing libmeshmonk_shared.dylib
meshmonk_lib_path = ['<path_to_meshmonk_root>/build/']; % SET THIS PATH
current_dyld_path = getenv('DYLD_LIBRARY_PATH');
if isempty(current_dyld_path)
    setenv('DYLD_LIBRARY_PATH', meshmonk_lib_path);
else
    setenv('DYLD_LIBRARY_PATH', [meshmonk_lib_path, ':', current_dyld_path]);
end
```
Again, adjust `<path_to_meshmonk_root>`.

## Using MeshMonk as a C++ Library

If you want to use MeshMonk's functionalities in your own C++ projects:

1.  **Build `libmeshmonk_shared`:** Follow the CMake instructions in "Building MeshMonk".
2.  **Include Headers:**
    *   Ensure your compiler can find `meshmonk.hpp`, `global.hpp`, and headers from the `src/` and `vendor/` directories of MeshMonk.
    *   You will also need to include Eigen3 headers.
3.  **Link against the Library:**
    *   Link your project against `libmeshmonk_shared.[so|dylib|dll]`.
    *   Your project will also need to link against Eigen3 if you use Eigen types in your interface code.
    *   OpenMesh is statically linked into `libmeshmonk_shared`, so you typically won't need to link against OpenMesh libraries directly.

## Documentation

*   **Legacy Build Instructions:** For older, manual build instructions (not using the primary CMake system described above), you can refer to files in the `docs/` directory. Please note these are largely outdated for the current build process.
    *   [Ubuntu (Legacy)](docs/ubuntu.md)
    *   [OSX (Legacy)](docs/osx.md)
    *   [Windows (Legacy)](docs/windows.md)

---
*Previous content about older installation methods and specific toolbox requirements for Matlab has been integrated or superseded by the CMake build process and updated MATLAB usage section.*
