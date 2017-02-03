# Meshmonk

## 1. Setting up your python environment
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

#### 1.1 Troubleshooting
If you get something like `/usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'GLIBCXX_3.4.20' not found` then follow instructions [here](http://askubuntu.com/questions/575505/glibcxx-3-4-20-not-found-how-to-fix-this-error). In most cases, you have to go to your environment (e.g. `source activate py27`) and run `$ conda install libgcc`


## 2. Install with Anaconda

Pre-requisites: download and install [Anaconda](https://www.continuum.io/downloads).

```bash
$ conda install -c conda-forge openmesh=6.3
$ conda install -c lukin0110 meshmonk=0.0.4
```

## License
This project is licensed under the terms of the MIT license.
