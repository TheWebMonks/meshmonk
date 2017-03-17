# Installing meshmonk (Ubuntu 16.04)
The following is a dummy-proof guide to install meshmonk on your ubuntu 16.04 machine.

## Install git
`sudo apt install git'

## Clone meshmonk
1) Make a folder 'projects' in home (`/home/user/projects/')
2) Go into the folder: `cd /home/user/projects/`)
3) Clone the online repository into the current projects folder: `git clone https://...`

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



# Using meshmonk
-add '-lmeshmonk' as an option to your compiler when compiling your software that uses the meshmonk library.
-include the meshmonk.hpp header

every time you recompile, don't forget you have to copy the latest shared library to /usr/local/lib/: `sudo cp /home/user/projects/meshmonk/meshmonk/bin/Release/libmeshmonk.so /usr/local/lib/`

