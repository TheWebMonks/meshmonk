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
