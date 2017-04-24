# Installing meshmonk (Ubuntu 16.04)
The following is a dummy-proof guide to install meshmonk on your ubuntu 16.04 machine.

## Install git
`sudo apt install git`

## Install gcc/g++ 4.9
Because matlab only supports gcc's and g++'s earlier version (not the 5.x versions), we'll have to downgrade.
First, we'll remove the symbolic links of gcc and g++ to the 5.x versions, we'll install the older versions, and relink gcc and g++ symbols to those older versions:
1. Check which versions you have. If they are 4.9, you can skip the other steps:
`gcc -v`
`g++ -v`
2. Remove the symbolic links:
`sudo rm /usr/bin/gcc`
`sudo rm /usr/bin/g++`
3. Install the older versions:
`sudo apt-get install gcc-4.9`
`sudo apt-get install g++-4.9`
4. Set symbolic links to the right binaries
`sudo ln -s /usr/bin/gcc-4.9 /usr/bin/gcc`
`sudo ln -s /usr/bin/g++-4.9 /usr/bin/g++`
5. Check if the versions now are indeed 4.9
`gcc -v`
`g++ -v`

## Install code::blocks
We'll use an IDE to make compiling the code easier.
`sudo add-apt-repository ppa:damien-moore/codeblocks-stable`
`sudo apt-get update`
`sudo apt-get install codeblocks codeblocks-contrib`

## Install required libraries
The meshmonk library depends on [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), [nanoflann](https://github.com/jlblancoc/nanoflann),
 and [OpenMesh](https://www.openmesh.org/). We'll install these one by one, which typically requires:
 1) Getting the source code
 2) Compiling it (in OpenMesh's case)
 3) Putting the header files and compiled libs in the right directory
 4) Making sure those directories are searched by the compiler/linker during the compiling/linking process.
 
 Note: if you want to learn more about how libraries work and how to use them, we highly recommend [David Wheeler's How-To](https://www.dwheeler.com/program-library/Program-Library-HOWTO/t1.html). Some other useful sources are:
 * How header files matter for shared libraries: [This stackoverflow answer](http://stackoverflow.com/a/1186836)
 
 
 ### Installing Eigen
 Eigen is a header-only library, which makes it very easy to use in your own projects. The only thing you have to do is download the header files and put them in a directory which your compiler searches so that it finds whatever it needs from the Eigen library.
 1) Go to their webpage and download the latest stable release (3.3.3 at moment of writing). It will end up as an archive in your `/home/user/Downloads/` folder.
 2) Extract it locally.
 3) Move the Eigen subfolder (which contains all the headers) to /usr/local/include/ so that you end up with /usr/local/include/Eigen/: `sudo mv /home/user/Downloads/eigen-eigen-67e894c6cd8f/Eigen /usr/local/include/`
 
 ### Installing nanoflann
 Like Eigen, nanoflann is also a header-only library. So follow practically the same procedure:
 1) Download the zipped source files from their [repository](https://github.com/jlblancoc/nanoflann) (see the green 'clone or download' button -> Download Zip)
 2) Extract it locally
 3) Move it to /usr/local/include/: `sudo mv /home/user/Downloads/nanoflann-master/include/nanoflann.hpp /usr/local/include/`
 
 ### Installing OpenMesh
OpenMesh can be compiled and used as both a static (.a) and shared library (.so). We'll download the source, configure build settings with 'cmake', and build it using 'make'. After all that's happened, we move the compiled library files to /usr/local/lib/ and the header files to /usr/local/include/. And finally, we'll have to update some configurations so that these library files will be found when compiling/linking to it later.
1) Download the latest version of the OpenMesh sourcefiles from their [downloads page](https://www.openmesh.org/download/) (the .tar.gz or .tar.bz2 file)
2) Extract it locally
3) Go into the extracted folder and make a new folder inside called 'build'
4) In your terminal, switch to the newly created 'build' folder: `cd /home/user/Downloads/OpenMesh-6.3/build/`
5) Configure the makefiles using cmake: `cmake ..` (include the two dots!)
6) Build the library: `make`
7) Move the all the compiled library files to /usr/local/lib/: `sudo mv /home/user/Downloads/OpenMesh-6.3/build/Build/lib/* /usr/local/lib/`
8) Move the folder with the header files we need to /usr/local/include/: `sudo mv /home/user/Downloads/OpenMesh-6.3/src/OpenMesh/ /usr/local/include/`
9) Run ldconfig so that your library loader can find it when running an application that needs the library: `sudo ldconfig -v`. To check, run `ldconfig -p | grep OpenMesh` and it should print a few library names containing OpenMeshCore and OpenMeshTools

### Installing MeshMonk
First, clone the meshmonk repository:
1) Make a folder 'projects' in home (`/home/user/projects/`)
2) Go into the folder: `cd /home/user/projects/`)
3) Clone the online repository into the current projects folder: `git clone https://...`

Next, let's compile MeshMonk using Code::Blocks.
Compiler options:
1) add option to use c++14 : `-std=c++14`
2) make the compiled code position independent: `-fPIC`
3) add verbose option to the linker (to inspect what's happening): `-Wl,-V`
4) Optimize the compilation for speed using `-O2` (not `-O3` !)

# Using meshmonk
-add '-lmeshmonk -lOpenMeshCore -lOpenMeshTools' as an option to your compiler when compiling your software that uses the meshmonk library.
-include the meshmonk.hpp header

every time you recompile, don't forget you have to copy the latest shared library to /usr/local/lib/: `sudo cp /home/user/projects/meshmonk/meshmonk/bin/Release/libmeshmonk.so /usr/local/lib/`

## Matlab
change gcc version to 4.9: [explanation](http://askubuntu.com/a/26502/664811)

Compile while linking to meshmonk: `mex functionFile.cpp -lmeshmonk`

### (Pre-)Loading required libraries
Setting the library paths inside Matlab has some unresolved [problems](https://nl.mathworks.com/matlabcentral/newsreader/view_thread/253412). It seems overwriting the library paths to use inside matlab doesn't work. So instead, we'll preload the necessary libs when starting Matlab:

So start matlab from the terminal with the command `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6:/usr/local/lib/libmeshmonk.so:/usr/local/lib/libOpenMeshCore.so:/usr/local/lib/libOpenMeshTools.so matlab` to make sure all the libraries are loaded.
