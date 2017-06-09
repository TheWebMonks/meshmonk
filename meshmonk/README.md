# Installing meshmonk
MeshMonk can be built and installed on [Ubuntu 16.04](https://github.com/TheWebMonks/meshmonk/blob/master/docs/ubuntu.md) and [OSX](https://github.com/TheWebMonks/meshmonk/blob/master/docs/osx.md).

# Using meshmonk

## From Matlab

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
In matlab, just run the mex_all.m script to mex all the meshmonk functions you need. This file can be found in the meshmonk repository (should be in your projects/meshmonk/matlab/ folder).

## From other software
If you're creating your own c++ project and want to use meshmonk, simply add '-lmeshmonk -lOpenMeshCore -lOpenMeshTools' as an option to your compiler when compiling your software that uses the meshmonk library.
-include the meshmonk.hpp header
