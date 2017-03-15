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
 
 Note: if you want to learn more about how libraries work and how to use them, we highly recommend [David Wheeler's How-To](https://www.dwheeler.com/program-library/Program-Library-HOWTO/t1.html)
 
 ### Installing Eigen
 
 ### Installing nanoflann
 
 ### Installing OpenMesh
