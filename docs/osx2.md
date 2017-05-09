# Build on OSX

## Pre-requisites

Install *brew* on OSX. 

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

## Install libraries

### 1. Install eigen

```bash
brew install eigen
```

### 2. 'Install' nanoflann

```bash
git clone git@github.com:jlblancoc/nanoflann.git
cp nanoflann/include/nanoflann.hpp /usr/local/include
```

### 3. Install OpenMesh

```bash
wget http://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3.zip
unzip OpenMesh-6.3.zip
cd OpenMesh-6.3
mkdir build
cd build
cmake ..
make
```

Copy/past the *.dylib* files:
```bash
sudo sudo cp Build/lib/*.* /usr/local/lib/
```

## Build!

```bash
make 
```
