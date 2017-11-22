# Build on OSX

## Pre-requisites
### Brew
Install *brew* on OSX. 

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

### wget
After installing brew, we can use brew to install wget (so you can download files from the command-line):
```bash
brew install wget --with-libressl
```

### CMake
We need cmake to configure and build OpenMesh and MeshMonk. We're copying here the installation instructions from [this page](http://tudat.tudelft.nl/projects/tudat/wiki/Install_on_Mac_OS_X#Install-CMake-on-Mac-OS-X):
1. Download the latest binary distribution for Mac OSX from the [Downloads page](https://cmake.org/download/)
2. Open the resulting .dmg-file you just downloaded
3. Drag and drop the CMake icon to the Applications folder icon.
4. Launch CMake (from you Applications folder, for example)
5. Add CMake to your path so you can call it from command-line:
```bash
sudo mkdir -p /usr/local/bin
sudo /Applications/CMake.app/Contents/bin/cmake-gui --install=/usr/local/bin
```
6. Verify your installation (you might have to close and open a new Terminal):
```bash
cmake --version
```

## Install libraries

### 1. Install eigen
1. Use brew to install eigen
```bash
brew install eigen
```
2. Make a symlink to the main 'Eigen' folder from /usr/local/include/ so that you can include Eigen in your projects.
```bash
ln -s /usr/local/Cellar/eigen/3.3.3/include/eigen3/Eigen /usr/local/include/Eigen
```

### 2. _(optional)_ 'Install' nanoflann
In your terminal, go to your Downloads folder. Download the nanoflann header (which is the only file you need from the nanoflann library) and copy it to /usr/local/include/.
```bash
cd /Users/user/Downloads/
wget https://raw.githubusercontent.com/jlblancoc/nanoflann/master/include/nanoflann.hpp
cp nanoflann.hpp /usr/local/include
```

### 3. Install OpenMesh

```bash
cd /Users/user/Downloads/
wget http://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3.zip
unzip OpenMesh-6.3.zip
cd OpenMesh-6.3
mkdir build
cd build
cmake ..
make
```

Copy/paste the *.dylib* files:
```bash
sudo sudo cp Build/lib/*.* /usr/local/lib/
```

Copy/paste the header files
```bash
cp -R /Users/user/Downloads/OpenMesh-6.3/src/OpenMesh /usr/local/include/
```

## MeshMonk

### Compile
Make a 'projects' folder in your home directory and go into it(or go into a directory you are already using for GitHub projects):
```bash
cd
mkdir projects
cd projects
```
Clone the online MeshMonk repository
```
git clone https://github.com/TheWebMonks/meshmonk.git
```
Enter the subfolder with the c++ sourcecode and use `make` to compile it
```bash
cd meshmonk
make
```
From the meshmonk folder, copy the shared library object libmeshmonk.dylib to /usr/local/lib
```bash
cp libmeshmonk.dylib /usr/local/lib
```

Copy the header files to /usr/local/include/
```bash
(cd /Users/user/projects/meshmonk/ && find . -name '*.hpp' -print | tar --create --files-from -) | (cd /usr/local/include/ && sudo tar xvfp -)
```

