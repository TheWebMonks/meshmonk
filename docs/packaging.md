# Packaging

[How to create Anaconda packages](https://conda.io/docs/build_tutorials/pkgs2.html).

## 1. Create a new tag it git

```bash
$ git tag -a v0.0.1 -m "first conda version"
$ git push origin --tags
```

NOTE: Change version in the `conda/meta.yml` to your new tag version.

## 2. Create package

```bash
$ cd conda
$ conda build .
```

## 3. Upload package

Pre-requisites:

```bash
$ conda install anaconda-client
```

Upload for OSX:
```bash
$ anaconda login
$ anaconda upload /Applications/anaconda/conda-bld/osx-64/meshmonk-0.0.3-np111py27_0.tar.bz2
```

Upload for Linux:
```bash
$ conda convert --platform all /Applications/anaconda/conda-bld/osx-64/meshmonk-0.0.3-np111py27_0.tar.bz2 -o outputdir/
$ anaconda upload outputdir/linux-64/meshmonk-0.0.3-np111py27_0.tar.bz2 
```

# Virtual envs on OSX

http://sourabhbajaj.com/mac-setup/Python/virtualenv.html

```bash
$ virtualenv venv --distribute --system-site-packages
```

These commands create a venv subdirectory in your project where everything is installed. You need to activate it first 
though (in every terminal where you are working on your project):

```bash
$ source venv/bin/activate
```

You should see a (venv) appear at the beginning of your terminal prompt indicating that you are working inside the 
virtualenv. Now when you install something:

```bash
$ pip install <package>
```

It will get installed in the venv folder, and not conflict with other projects.

To leave the virtual environment use.

```bash
$ deactivate
```

# TODO

/Applications/anaconda/lib/python2.7/site-packages/requests/packages/urllib3/util/ssl_.py:132: InsecurePlatformWarning: 
A true SSLContext object is not available. This prevents urllib3 from configuring SSL appropriately and may cause 
certain SSL connections to fail. You can upgrade to a newer version of Python to solve this. For more information, 
see https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
  InsecurePlatformWarning
