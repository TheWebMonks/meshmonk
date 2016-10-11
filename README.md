# KULeuven Algoritms

## How to

In order to get started you need the [VTK C++](http://www.vtk.org) and 
a lot of other C++ and Linux stuff :) All this has been wrapped in a 
reproducible docker image. The docker image is hosted on 
[Docker Hub](https://hub.docker.com/r/lukin0110/docker-vtk-python/).

### Download

```
$ docker pull lukin0110/docker-vtk-python
```

### Start a bash
```
$ docker run -it -v src:/src lukin0110/docker-vtk-python bash
```

Now you're in a Linux bash on Debian Wheezy with all the necessary tools
& libs installed. The local `src` folder is mounted to it. With that 
shell you can just run your scripts now.

```
$ python main.py
```

### Setting up your python environment
Instead of installing python and all packages you will need (like scipy 
and numpy) yourself, use Anaconda to handle everything for you. Go to 
the [download page](https://www.continuum.io/downloads) and get the 
installer for python 3.5.

After you've done this, it's best to create two separate environment for 
python 2.7 and python 3.5. Do that using the terminal:
```
$ conda create -n py27 python=2.7 anaconda
```
and afterwards also for python 3.5
```
$ conda create -n py35 python=3.5 anaconda
```

You can find more information [here](http://conda.pydata.org/docs/using/envs.html#list-all-environments) 
about how to switch between environments, see which one is active, 
etcetera.