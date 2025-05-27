# Note
This project is now maintained at KU Leuven GitLab: https://gitlab.kuleuven.be/mirc/meshmonk

# Installing meshmonk
MeshMonk can be built and installed on:

* [Ubuntu 16.04](docs/ubuntu.md)
* [OSX](docs/osx.md)
* [Windows](docs/windows.md)

# Using meshmonk

## From Matlab
The following toolboxes are required from Matlab:
* Statistics and Machine Learning Toolbox
* Image Processing Toolbox

### Setting Environment Variables

#### Ubuntu
Setting the library paths inside Matlab has some unresolved [problems](https://nl.mathworks.com/matlabcentral/newsreader/view_thread/253412). It seems overwriting the library paths to use inside matlab doesn't work. So instead, we'll preload the necessary libs when starting Matlab:

So start matlab from the terminal with the command `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6:/usr/local/lib/libmeshmonk.so:/usr/local/lib/libOpenMeshCore.so:/usr/local/lib/libOpenMeshTools.so matlab` to make sure all the libraries are loaded.

#### Mac OSX
Add two lines to your startup.m file:
```
p = ['usr/local/lib/'];
setenv('LD_LIBRARY_PATH', p);
```

### Mexing meshmonk functions
On ubuntu/mac: 
In matlab, go into the projects/meshmonk/matlab folder (or wherever you put the meshmonk repository) and run the mex_all.m script to mex all the meshmonk functions you need.

On Windows: 
In matlab, go into the Documents/GitHub/meshmonk folder (or wherever you put the meshmonk repository) and run the mex_windows_all.m script to mex all the meshmonk functions you need.

## From other software
If you're creating your own c++ project and want to use meshmonk, simply add '-lmeshmonk -lOpenMeshCore -lOpenMeshTools' as an option to your linker when compiling your software that uses the meshmonk library.

-include the meshmonk.hpp header

## Demo
An example of a facial registration can be found in the demo folder

## Using the C++ Command-Line Tool (`meshmonk_cli`)

This section describes how to build the C++ command-line interface (`meshmonk_cli`) for MeshMonk, which allows you to perform registrations directly without MATLAB.

### Prerequisites

Ensure you have the following dependencies installed on your system (Linux examples provided):

1.  **C++17 Compiler:** A modern C++ compiler (e.g., GCC >= 7 or Clang >= 5).
    ```bash
    sudo apt update
    sudo apt install build-essential g++
    ```
2.  **CMake:** Version 3.16 or higher.
    ```bash
    sudo apt install cmake
    ```
3.  **Eigen3:** A library for linear algebra.
    ```bash
    sudo apt install libeigen3-dev
    ```
4.  **cxxopts:** A lightweight C++ option parser library.
    ```bash
    sudo apt install libcxxopts-dev
    ```
    (Note: If `libcxxopts-dev` is not available or you need a specific version, you might need to install it from source: [https://github.com/jarro2783/cxxopts](https://github.com/jarro2783/cxxopts))

### Compilation Steps

1.  **Clone the Repository (if you haven't already):**
    ```bash
    # If you cloned via SSH:
    # git clone git@gitlab.kuleuven.be:mirc/meshmonk.git
    # If you cloned via HTTPS:
    git clone https://gitlab.kuleuven.be/mirc/meshmonk.git
    cd meshmonk 
    ```
    (Adjust the clone command if you obtained the source code differently, e.g., from a fork or a specific branch).

2.  **Create a Build Directory:**
    It's good practice to build outside the source directory.
    ```bash
    mkdir build
    cd build
    ```

3.  **Run CMake:**
    This command configures the project and generates the build files. It will also compile the vendored OpenMesh library.
    ```bash
    cmake ..
    ```
    If you encounter errors related to missing dependencies (Eigen3, cxxopts), please ensure they are correctly installed and discoverable by CMake.

4.  **Compile the Project:**
    This will build the `meshmonk_shared` library and the `meshmonk_cli` executable.
    ```bash
    make -j$(nproc)  # Uses all available CPU cores for a faster build
    ```
    On successful compilation, the executable will be located at `build/cli/meshmonk_cli`.

### Running `meshmonk_cli`

Once compiled, you can run the `meshmonk_cli` tool from the `build` directory.

**Basic Syntax:**

```bash
./cli/meshmonk_cli <command> <source_mesh.obj> <target_mesh.obj> <output_mesh.obj> [options...]
```

**Available Commands:**

*   `pyramid_reg`: Performs pyramid-based non-rigid registration.
*   `rigid_reg`: Performs rigid registration.

You can see all available options for each command by running:
```bash
./cli/meshmonk_cli --help
```

**Example: Rigid Registration**

Let's register `Template.obj` (source) to `demoFace.obj` (target) using rigid registration. This example also saves the computed transformation matrix.

1.  **Navigate to the build directory (if not already there):**
    ```bash
    # Assuming you are in the MeshMonk root directory:
    cd build 
    ```

2.  **Run the `rigid_reg` command:**
    ```bash
    ./cli/meshmonk_cli rigid_reg ../demo/Template.obj ../demo/demoFace.obj rigid_output.obj --transform_output rigid_transform.txt
    ```
    *   `../demo/Template.obj`: Path to the source mesh (relative to `build` directory).
    *   `../demo/demoFace.obj`: Path to the target mesh (relative to `build` directory).
    *   `rigid_output.obj`: Filename for the transformed source mesh that will be saved in the current directory (`build/`).
    *   `--transform_output rigid_transform.txt`: Specifies that the 4x4 transformation matrix should be saved to `rigid_transform.txt` in the current directory (`build/`).

3.  **Expected Output:**
    After running the command, you should find the following files in your `build` directory:
    *   `rigid_output.obj`: The source mesh (`Template.obj`) transformed to align with the target mesh (`demoFace.obj`).
    *   `rigid_transform.txt`: A text file containing the 4x4 transformation matrix. For example:
        ```
         0.951203 -0.171045 -0.256844   2.13937
         0.115929  0.969435 -0.216259  -6.74492
         0.285982  0.175929  0.941953   58.7142
                0         0         0         1
        ```
        (The exact numerical values in the matrix may vary slightly based on the build environment or minor code variations but should be similar.)

**Example: Pyramid (Non-Rigid) Registration**

To perform a pyramid non-rigid registration:

```bash
./cli/meshmonk_cli pyramid_reg ../demo/Template.obj ../demo/demoFace.obj pyramid_output.obj --num_iterations 30 --smoothness 0.5
```
*   This will save the non-rigidly transformed source mesh to `pyramid_output.obj`.
*   It uses a reduced number of iterations (`--num_iterations 30`) and a specific smoothness value (`--smoothness 0.5`) as an example of overriding default parameters. Check `./cli/meshmonk_cli --help` for all available parameters and their defaults.

**Note on Demo Mesh Data:**
When running with the provided demo meshes (`Template.obj`, `demoFace.obj`), you might see some warnings in the console, such as:
*   `Warning! Material file ... .mtl' not found!`
*   OpenMesh internal errors like `PolyMeshT::add_face: complex vertex` or `complex edge` (especially for `demoFace.obj`).
These warnings are related to the demo data itself (missing material files, non-manifold geometry) and typically do not prevent the `meshmonk_cli` tool from completing the registration and producing output.

