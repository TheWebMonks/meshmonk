# Build on Windows

This installation guide works on a Windows 10 64bit machine with the
Visual Studio 2017 community edition. Ensure that you have to
**C++/CLI** tools installed.

## Install libraries

* [OpenMesh 6.3](http://openmesh.org/download/): download the vs2015 64bit with apps static [installer](http://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3-VS2015-64-Bit.exe)
* [Eigen 3.3.4](http://eigen.tuxfamily.org): download the [source zip](http://bitbucket.org/eigen/eigen/get/3.3.4.zip) and extract in *c:\Users\\\<username\>\Documents\GitHub\eigen-eigen-3.3.4*

## Configure Visual Studio project

### Import project
* File -> New -> Project From existing code
* Choose *Visual c++*
* Choose *meshmonk* as name
* Choose Project Type: *Static Library (LIB) project*

### Configure Solution
* Choose for Solution Configurations: *Release*
* Choose for Solution Platform: *x64*

### Exclude files from the project

* all *matlab* files
* nanaflann.hpp
* example.cpp

### Configure properties

Open *Project -> Properties*:

* Go to *General -> Project Defaults -> Configuration Type*: *change to Static library(.lib)*
* Go to *C++ -> Preprocessor -> Preprocessor Definitions*: add *_USE_MATH_DEFINES*
* Go to *C++ -> General -> Additional Include Libraries*: add
    ```
    C:\Users\lukin0110\Documents\GitHub\eigen-eigen-3.3.4
    C:\Users\lukin0110\Documents\GitHub\meshmonk\vendor
    C:\Program Files\OpenMesh 6.3\include
    ```
    This includes Eigen, OpenMesh and Nanoflann
* Go to *Librarian -> General -> Additional Dependencies*:
    ```
    C:\Program Files\OpenMesh 6.3\lib\OpenMeshCore.lib
    C:\Program Files\OpenMesh 6.3\lib\OpenMeshTools.lib
    ```

### Matlab

After you've build the library you can have to copy/paste the file to
a folder which is in matlab path. Either the root folder of matlab or
the relative root folder of the matlab files that you're working on.
