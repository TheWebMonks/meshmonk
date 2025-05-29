
**Phased Development Plan: MeshMonk Modernization & Python Integration**

**Overall Project Goal:** To modernize the MeshMonk C++ library by transitioning to a CMake build system, enhance its usability with Python bindings, and improve repository organization, all while ensuring the continued functionality of its existing Command-Line Interface (CLI) and MATLAB MEX support.

---

**Phase 1: Repository Restructuring & Core C++ Library CMake Foundation**

*   **Objective:**
    To implement the new, organized repository directory structure and establish a foundational CMake build system capable of compiling the core `libmeshmonk_shared` C++ library and its accompanying C++ examples. This phase focuses solely on the core library and its direct C++ usage, deferring integration of the CLI, MATLAB, and Python interfaces.

*   **Expected Outcomes & Specifications at End of Phase 1:**
    1.  **New Repository Structure Implemented:** The repository's directory structure will be reorganized as follows. This includes moving existing core library source files, public headers, and C++ examples to their new, designated locations.
        ```
        meshmonk/
        ├── .git/
        ├── .gitignore             # Must include 'build/' (the CMake temporary/output directory)
        ├── CMakeLists.txt         # Root CMake: calls add_subdirectory(library)
        ├── README.md              # Main project README (can be a placeholder or initial draft)
        │
        ├── library/               # Encapsulates the core C++ shared library (libmeshmonk_shared)
        │   ├── CMakeLists.txt     # CMakeLists.txt for building libmeshmonk_shared & its examples
        │   ├── include/           # PUBLIC headers exposed by libmeshmonk_shared
        │   │   └── meshmonk/      # Namespacing for public include files
        │   │       ├── meshmonk.hpp # Main public API header (moved from root)
        │   │       └── global.hpp   # Public global definitions (moved from root)
        │   ├── src/               # PRIVATE implementation files (.cpp) and internal headers (.hpp)
        │   │   ├── meshmonk.cpp   # Implementation of meshmonk.hpp functions (moved from root)
        │   │   ├──                 # All files from former top-level src/* are now here:
        │   │   ├── BaseCorrespondenceFilter.cpp
        │   │   ├── BaseCorrespondenceFilter.hpp # Example internal header
        │   │   └── ...            # All other .cpp files and any internal-only .hpp files
        │   └── examples/          # Standalone C++ examples demonstrating direct use of libmeshmonk_shared
        │       ├── CMakeLists.txt # CMake to build these C++ examples
        │       └── example.cpp    # (Moved from previous top-level examples/)
        │
        ├── cli/                   # Placeholder for C++ Command-Line Interface (to be integrated in Phase 2)
        │   └── (Existing CLI files, if any, remain here but are not yet part of CMake build)
        │
        ├── matlab/                # Placeholder for MATLAB MEX Interface (to be adapted in Phase 3)
        │   └── (Existing MATLAB files, if any, remain here but are not yet part of CMake build)
        │
        ├── python/                # Placeholder for Python Bindings (to be implemented in Phase 4)
        │
        ├── data/                  # Shared data for all examples and demos
        │   ├── Template.obj       # (Moved here)
        │   └── demoFace.obj       # (Moved here)
        │   └── ...                # Other sample/test data files
        │
        ├── docs/                  # Documentation (legacy build guides, etc.)
        │   ├── osx.md
        │   ├── ubuntu.md
        │   └── windows.md
        │
        ├── vendor/                # Third-party libraries
        │   ├── OpenMesh-11.0.0/
        │   └── nanoflann.hpp
        │
        └── builds/                # LEGACY: Contains pre-compiled binaries from the older, manual build system.
                                   # This directory is NOT the output target for the new CMake build process.
        ```
    2.  **Core Library CMake Build System (`library/CMakeLists.txt` & Root `CMakeLists.txt`):**
        *   A root `CMakeLists.txt` is established, primarily responsible for project-wide settings (like C++ standard C++14) and including the `library/` subdirectory.
        *   A `library/CMakeLists.txt` is created and configured to:
            *   Build the `meshmonk_shared` shared library from source files now located in `library/src/`.
            *   Correctly set up include directories for public headers (`library/include/`), internal headers (`library/src/`), and vendored dependencies (`vendor/`).
            *   Manage dependencies:
                *   **Eigen3:** Integrated using `FetchContent` (declared in the root or `library/` CMakeLists.txt).
                *   **OpenMesh:** Built statically from `vendor/OpenMesh-11.0.0` (via `add_subdirectory`) and statically linked into `libmeshmonk_shared`. `BUILD_SHARED_LIBS OFF` and `OPENMESH_BUILD_APPS OFF` for OpenMesh.
                *   **nanoflann.hpp:** Made available via include directories.
    3.  **C++ Examples Build (`library/examples/CMakeLists.txt`):**
        *   The `CMakeLists.txt` in `library/examples/` is configured to build the C++ example executables (e.g., `example.cpp`).
        *   These examples successfully link against the `meshmonk_shared` library target.
    4.  **Successful Compilation:**
        *   The entire project, when configured with CMake from a clean `build/` directory (e.g., `mkdir build && cd build && cmake .. && make`), compiles `libmeshmonk_shared` and the C++ examples in `library/examples/` without errors or critical warnings.
    5.  **Corrected Include Paths:** All `#include` directives within the C++ source code (`.hpp` and `.cpp` files) in the `library/` directory are updated to reflect the new file locations and public/private header organization (e.g., using `#include "meshmonk/meshmonk.hpp"` for public API).
    6.  **Verification:** The compiled C++ examples from `library/examples/` can be run and function as previously expected (if they had prior functionality).

---

**Phase 2: CLI Integration & Build**

*   **Objective:**
    To integrate the existing `meshmonk_cli` tool into the new CMake build system, ensuring it correctly links against the `libmeshmonk_shared` library (built in Phase 1) and remains fully functional.

*   **Expected Outcomes & Specifications at End of Phase 2:**
    1.  **CLI CMake Integration (`cli/CMakeLists.txt`):**
        *   A `CMakeLists.txt` is created or updated in the `cli/` directory.
        *   This CMake file defines an executable target for `meshmonk_cli` using `cli/cli.cpp`.
        *   The `meshmonk_cli` target correctly links against the `meshmonk_shared` library target.
        *   It handles the `cxxopts` dependency (preferably via system package, with guidance for manual install if needed).
    2.  **Root CMake Update:** The root `CMakeLists.txt` is updated to include the `cli/` subdirectory using `add_subdirectory(cli)`.
    3.  **CLI Source Code Adaptation:**
        *   `cli/cli.cpp` is updated to correctly `#include` headers from `libmeshmonk_shared` (e.g., `#include "meshmonk/meshmonk.hpp"`) and any other dependencies, reflecting the new project structure.
    4.  **Successful Compilation:** The `meshmonk_cli` executable compiles successfully as part of the main CMake build process.
    5.  **Full CLI Functionality Preserved:**
        *   The `meshmonk_cli` tool is fully operational and produces the same results as before the restructuring when tested with the demo data from the `data/` directory.
        *   Supported commands (`pyramid_reg`, `rigid_reg`), command-line arguments, input/output file handling, and default parameter behaviors remain consistent with its pre-restructuring state (or are intentionally aligned with updated defaults from MATLAB scripts, as per main spec).
    6.  **Documentation Update (Initial):** The main `README.md` section regarding building and running `meshmonk_cli` is reviewed and updated to reflect any changes necessitated by the CMake build process (e.g., output location of the executable in the `build/` directory), and includes an example for running both the rigid and nonrigid registration via the CLI.

---

**Phase 3: MATLAB MEX Adaptation**

*   **Objective:**
    To adapt the existing MATLAB MEX build scripts to make them compatible with the CMake-built `libmeshmonk_shared` library. This ensures that MATLAB users can continue to compile and use MeshMonk's MEX functions.

*   **Expected Outcomes & Specifications at End of Phase 3:**
    1.  **Updated MATLAB MEX Scripts (`demo/matlab/mex_*.m`):**
        *   The scripts `demo/matlab/mex_all.m` and `demo/matlab/mex_windows_all.m` are modified.
        *   **Include Paths (`-I`):** All `mex` commands within these scripts use correct include paths pointing to:
            *   Public headers of `libmeshmonk_shared` (e.g., `library/include/`).
            *   Internal headers if used by MEX files (e.g., `library/src/`, `vendor/`).
            *   Eigen3 headers (system or `FetchContent` location).
            *   OpenMesh headers (from the CMake `build/` directory, e.g., `build/vendor/OpenMesh-11.0.0/_build/src/`).
        *   **Library Paths (`-L`):** All `mex` commands use correct library paths pointing to the directory containing the CMake-built `libmeshmonk_shared` (e.g., the CMake `build/` directory).
        *   **Linked Libraries (`-l`):** All `mex` commands link against `meshmonk_shared` and any other necessary system libraries (e.g., `stdc++`). Explicit linking against `OpenMeshCore` and `OpenMeshTools` is removed at the MEX stage.
    2.  **MATLAB Demo Script Paths:**
        *   MATLAB demo scripts (`demo/matlab/test_*.m`) are updated to load mesh data from the common `data/` directory (e.g., using relative paths like `../../data/`).
    3.  **Verification (Build Feasibility):**
        *   A "dry run" verification confirms that all specified paths in the modified MEX scripts point to existing files and directories after a successful CMake build of `libmeshmonk_shared`. The `mex` commands are structurally sound and *could* theoretically compile if run in a correctly configured MATLAB environment.
        *   *Note: Full MATLAB compilation and testing is ideal but confirming path correctness and linkage intent is the minimum for this phase if a MATLAB environment is not directly available to the implementer.*
    4.  **Documentation Update:** The main `README.md` section for "Using MeshMonk from MATLAB" is updated with:
        *   Clear instructions on first building `libmeshmonk_shared` via CMake.
        *   Detailed guidance on the necessary modifications to the MEX scripts.
        *   Revised instructions for MATLAB runtime environment setup (`LD_PRELOAD`, `DYLD_LIBRARY_PATH`/`startup.m`).

---

**Phase 4: Python Bindings Implementation**

*   **Objective:**
    To create Python bindings for the core MeshMonk registration pipelines using `pybind11`, making these functionalities accessible from Python.

*   **Expected Outcomes & Specifications at End of Phase 4:**
    1.  **pybind11 Integration (`python/CMakeLists.txt` & Root `CMakeLists.txt`):**
        *   `pybind11` is added as a dependency to the project, preferably via `FetchContent` in the root `CMakeLists.txt`.
        *   A `CMakeLists.txt` is created in the `python/` directory to define and build the Python extension module (e.g., `meshmonk_python.so` or `.pyd`).
        *   This module successfully links against `libmeshmonk_shared`.
    2.  **Python Binding Code (`python/meshmonk_bindings.cpp`):**
        *   This file is created and contains `pybind11` code to wrap:
            *   `meshmonk::pyramid_registration`
            *   `meshmonk::nonrigid_registration`
            *   `meshmonk::rigid_registration`
        *   Python function names are `snake_case`.
        *   C++ arguments and default values are exposed.
        *   In-place modification of NumPy arrays for `floating_features` and `transformation_matrix` is implemented.
        *   The `registration::NUM_FEATURES` constant is exposed as `meshmonk_python.NUM_FEATURES`.
        *   Standard C++ exceptions are convertible to Python exceptions.
    3.  **Successful Compilation:** The Python extension module compiles successfully as part of the main CMake build process (`make` in the `build/` directory).
    4.  **Basic Python Import and Call:**
        *   The compiled Python module can be imported into a Python interpreter (e.g., `import meshmonk_python`).
        *   At least one bound function can be called with placeholder or simple NumPy array inputs without crashing, demonstrating basic binding linkage.

---

**Phase 5: Python Demo/Test Scripts & Documentation Finalization**

*   **Objective:**
    To create functional Python demo/test scripts that exercise the Python bindings, and to finalize all project documentation for clarity and accuracy.

*   **Expected Outcomes & Specifications at End of Phase 5:**
    1.  **Python Demo/Test Scripts (`python/demo/`):**
        *   New Python scripts (e.g., `test_pyramid_registration.py`, `test_rigid_registration.py`) are created in `python/demo/`.
        *   These scripts:
            *   Load mesh data from the common `data/` directory.
            *   Prepare input data as NumPy arrays.
            *   Call the bound Python functions from the `meshmonk_python` module.
            *   Produce results comparable to the MATLAB demos.
            *   May include basic assertions or save output meshes for verification.
    2.  **Comprehensive `README.md`:**
        *   The main `README.md` is fully updated and polished, containing:
            *   Clear overview of MeshMonk and its usage methods.
            *   Accurate and complete build instructions for the entire project (core library, CLI, Python bindings) using CMake.
            *   Detailed usage instructions for `meshmonk_cli`.
            *   Revised and accurate instructions for "Using MeshMonk from MATLAB" (covering CMake prerequisite and MEX script adaptation).
            *   New section on "Using MeshMonk from Python" with installation of the bindings (if applicable beyond just building) and basic usage examples.
            *   Updated repository structure diagram and explanations.
    3.  **Legacy Documentation:**
        *   Files in `docs/` (e.g., `ubuntu.md`) have prominent disclaimers indicating they describe an older, manual build system and point to the main `README.md` for current CMake-based instructions.
    4.  **Build System Finalization:** The CMake build system is robust and cleanly builds all specified targets (core library, C++ examples, CLI, Python extension module).
