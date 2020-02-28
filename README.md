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

