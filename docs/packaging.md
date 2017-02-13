# Packaging

[How to create Anaconda packages](https://conda.io/docs/build_tutorials/pkgs2.html).

## 1. Create a new tag it git

```bash
$ git tag -a v0.0.5 -m "your text here"
$ git push origin --tags
```

## 2. Update meta.yml

Now you need to change the `git_rev` in the `conda/meta.yml` file to your new tag.

```bash
source:
  git_rev: v0.0.5
  git_url: https://github.com/TheWebMonks/meshmonk.git
```

## 3. Create package

```bash
$ cd conda
$ conda build .
```

## 4. Upload package

Pre-requisites:

```bash
$ conda install anaconda-client
```

Upload for OSX:
```bash
$ anaconda login
$ anaconda upload -u WebMonks /Applications/anaconda/conda-bld/osx-64/meshmonk-0.0.3-np111py27_0.tar.bz2
```

Upload for Linux:
```bash
$ conda convert --platform all /Applications/anaconda/conda-bld/osx-64/meshmonk-0.0.3-np111py27_0.tar.bz2 -o outputdir/
$ anaconda upload -u WebMonks outputdir/linux-64/meshmonk-0.0.3-np111py27_0.tar.bz2 
```

# TODO

/Applications/anaconda/lib/python2.7/site-packages/requests/packages/urllib3/util/ssl_.py:132: InsecurePlatformWarning: 
A true SSLContext object is not available. This prevents urllib3 from configuring SSL appropriately and may cause 
certain SSL connections to fail. You can upgrade to a newer version of Python to solve this. For more information, 
see https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
  InsecurePlatformWarning
