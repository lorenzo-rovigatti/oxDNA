# Installation

## Requirements

The code requires [CMake](https://cmake.org/) and a c++-14-compliant `g++` (any version >= 4.9 *should* work). The code should be also compilable with the Intel compiler (with the `-DIntel=ON` `cmake` flag, see below), although this has not been tested with newer oxDNA versions.

### CUDA

Compiling with CUDA support requires CMake >= 3.5 and a CUDA toolkit >= 10. If your current setup cannot meet these requirements we advise you to use [older versions of oxDNA](https://sourceforge.net/projects/oxdna/files/).

### Python bindings

Compiling with the python bindings enabled requires a working Python3 installation comprising both binaries and include files. With Debian-derived distro these come with the package `python3-dev`.  OxDNA Analysis Tools requires Python version 3.9 or newer.  `oxpy` will work with any version of Python 3.

## Compiling oxDNA

Clone the [repo](https://github.com/lorenzo-rovigatti/oxDNA.git) or extract the oxDNA archive and then:

```bash
cd oxDNA         # enter the oxDNA folder
mkdir build      # create a new build folder. It is good practice to compile out-of-source
cd build
cmake ..         # here you can specify additional options, see next section
make -j4         # compile oxDNA. The -jX make option makes it compile the code in parallel by using X threads.
```

At the end of the compilation three executables (*oxDNA*, *DNAnalysis* and *confGenerator*) will be placed in the `build/bin` directory. 

Compiling with Python bindings will also generate an `oxpy` package in the `build/oxpy` directory that can be imported in Python. Running `make install` will attempt to copy the package (as well as `oxDNA_analysis_tools`) to the `pip`'s module directory. The specific location will depend on your system's settings. We advise you to use [virtual environments](https://docs.python.org/3/tutorial/venv.html) (see *e.g.* [pipenv](https://docs.pipenv.org/)) to avoid conflicts with other packages and/or dependency and permission issues.

### Updating a local copy

If you cloned the repository to install oxDNA, your local copy can be updated with the following commands:

```bash
cd oxDNA        # enter the oxDNA folder
git pull        # use git to synchronize your repo with the online one
cd build        # enter the build folder (see above)
make -j4        # compile the updated source
```

If you also want to update `oxpy` don't forget to run `make install` after the compilation.

### CMake options

* `-DCUDA=ON` Enables CUDA support
* `-DCUDA_COMMON_ARCH=ON` Choose the target CUDA compute architecture based on the nvcc version. Set it to off to autodetect the CUDA compute arch GPU installed.
* `-DDebug=ON` Compiles with debug symbols and without optimisation flags
* `-DG=ON` Compiles with debug symbols + optimisation flags
* `-DINTEL=ON` Uses INTEL's compiler suite
* `-DMPI=ON` Compiles oxDNA with MPI support
* `-DSIGNAL=OFF` Handling system signals is not always supported. Set this flag to OFF to remove this feature
* `-DMOSIX=ON` Makes oxDNA compatible with MOSIX
* `-DDOUBLE=OFF` Set the numerical precision of the CPU backends to `float`
* `-DCUDA_DOUBLE=ON` Set the numerical precision of the CUDA backends to `double`, which is not compatible with the `mixed` precision.
* `-DNATIVE_COMPILATION=ON` Set to `OFF` to compile without the `-march=native` flag. This may be required when compiling binaries to be used elsewhere

The following options pertain to `oxpy`:

* `-DPython=ON` Enables Python bindings
* `-DOxpySystemInstall=On` By default `oxpy` is installed in the current user's home directory. By enabling this option `oxpy` will be installed as a system-wide package. It may require superuser privileges (unless using Conda environments. See below.).

If you are on your own machine or you installed Python via Anaconda, the `-DOxpySystemInstall=On` option should be set to install `oxpy` and `oxDNA_analysis_tools` on your `PATH`.  If you are on a shared system (for example and HPC cluster which isn't using Conda), perform the `cmake` command without this option.  `make install` will probably (depends on `pip` settings) install the two libraries in `$HOME/.local/bin/`, which is not on the `PATH` by default and you will have to add it (if you have not already for another program).  In our testing, Conda-managed versions of Python correctly installed the libraries in the current environment's `site-packages` folder if `-DOxpySystemInstall=On` was set, even on shared systems.

### make targets

* `make` compiles oxDNA
* `make install` Copies the `oxpy` and `oxDNA_analysis_tools` packages to the Python's package folder making the two packages importable with Python. 
* `make rovigatti` Compiles the observables and interactions in contrib/rovigatti
* `make romano` Compiles the observables and interactions in contrib/romano

### Known issues

#### `oxpy`

* When running `make` if you run into a problem at the end of compilation where it cannot find `python.h`, this means that you don't have the Python developer kit installed. See [this StackOverflow answer](https://stackoverflow.com/a/21530768/9738112) on how to install with a variety of package managers.
* When compiling with the Python bindings enabled CMake will sometimes choose the wrong Python binary and/or include files, resulting in a failed compilation. If this happens the correct paths can be directly set from the command line as follows:
	```bash
	cmake .. -DPython=ON -DPYTHON_INCLUDE_DIRS=/path/to/python/include/dir -DPYTHON_EXECUTABLE=/path/to/python/binary
	```
	If you are using conda environments, the paths should look something like:  
	* Include: `$HOME/anaconda3/envs/py38/include/python3.8`
	* Executable: `$HOME/anaconda3/envs/py38/bin/python`
* If you get an error regarding the number of bytes in the `numpy.array` header, this happens when the version of Numpy on your system doesn't match the version that pip downloads from PyPi when installing `oxDNA_analysis_tools` with its isolated environment (most commonly because you installed Numpy using Conda which tends to be a few versions behind PyPi). To fix this, either update your version of Numpy or try to install just `OAT` without build isolation:
  `python -m pip install ./analysis --no-build-isolation`
* Sometimes installation will fail with `TypeError: expected string or bytes-like object`. This error is usually caused by older versions of either `oxpy` or `oxDNA-analysis-tools` floating around. Remove them and re-install `oxpy`.
	
## Testing

* `make test_run` runs quick tests to check whether oxDNA has been correctly compiled or not.	
* `make test_quick` runs longer tests to check that oxDNA works.
* `make test_oxpy` checks that the Python bindings work.
* `make test` runs all sets of tests above.
* `oat config` will check all dependencies for `OAT`
* In `oxDNA/analysis/tests` there is a shell script, `test.sh` which will test your `OAT` installation.
