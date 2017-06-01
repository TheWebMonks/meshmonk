# Build on Ubuntu 16.04
The following is a dummy-proof guide to install meshmonk on your ubuntu 16.04 machine.

## Install git
`sudo apt install git`

## Install gcc/g++ 4.9
Because matlab only supports gcc's and g++'s earlier version (not the 5.x versions), we'll have to downgrade. Below we'll give an easy (but permanent) way to do that. If you want to be able to switch versions on the fly, however, follow [this post](https://askubuntu.com/a/26518/664811).

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

### Compiling MeshMonk
#### Getting the sources
First, clone the meshmonk repository:
1) Make a folder 'projects' in home (`/home/user/projects/`)
2) Go into the folder: `cd /home/user/projects/`)
3) Clone the online repository into the current projects folder: `git clone https://github.com/TheWebMonks/meshmonk.git`

#### Setting up the Shared library project
Next, let's compile MeshMonk using Code::Blocks.
1) Select 'Create a new project' and choose 'Shared library'
2) Choose C++
3) Fill in the project form:
* Project title: 'meshmonk'
* Folder to create project in: '/home/user/projects/meshmonk/'
* The other fields are filled in automatically, just verify the 'Resulting filename' equals '/home/user/projects/meshmonk/meshmonk/meshmonk.cbp'
4) Set configurations:
* Compiler: 'GNU GCC Compiler'
* Tick both the Debug and Release box

After clicking 'Finish', the meshmonk project is opened automatically. We're going to set some compiler/linker options now:
1) Right-click the project (called 'meshmonk' with the symbol of Code::Blocks in front of it) and select 'Build options...'
2) In the 'Compiler settings' tab, select the 'Compiler Flags' subtabtick the 'Optimize even more (for speed) [-O2]' option
3) In the 'Compiler settings' tab, select the 'Other compiler options' and write `-std=c++14 -Wl,-V -fPIC` in the text field.

One big thing is still missing from the project, namely the sources themselves! Delete the current main.cpp that is in the meshmonk project (which was automatically generated but we're gonna use our own sources). Now add the sources:
* Right-click the meshmonk project and select 'Add files recursively...'. Choose the `/home/user/projects/meshmonk/meshmonk/src` folder that you obtained by cloning the repository earlier.
* Right-click the meshmonk project and select 'Add files...'. Choose the meshmonk.hpp and meshmonk.cpp files located in the `home/user/projects/meshmonk/meshmonk/` folder.

#### Compile the project
Now, compile the code by clicking the small yellow cog in the top toolbar ('Build'). Make sure the version is set to 'Release' and not 'Debug' (should be to the right of the build buttons).

Code::Blocks will print a lot of output, including warnings (in blue). Don't worry about those.

### Installing meshmonk
Now that you've compiled everything, we're going to put the library files in the right places so that other applications can access them:
1) Copy the shared library 'libmeshmonk.so' to /usr/local/lib/: `sudo cp /home/user/projects/meshmonk/meshmonk/bin/Release/libmeshmonk.so /usr/local/lib/`
2) Copy the header files to /usr/local/include/: `(cd /home/user/projects/meshmonk/ && find . -name '*.hpp' -print | tar --create --files-from -) | (cd /usr/local/include/ && sudo tar xvfp -)`
3) Run ldconfig to update your library list: `sudo ldconfig -v`

Note(!): every time you recompile, don't forget you have to copy the latest shared library to /usr/local/lib/!
