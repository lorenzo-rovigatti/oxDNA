# Installation

## Requirements

The code requires `cmake` and a c++-14-compliant `g++` (any version >= 4.9 *should* work). The code should be also compilable with the Intel compiler (with the `-DIntel=ON` `cmake` flag, see below), although this has not been tested with newer oxDNA versions.

### CUDA

Compiling with CUDA support requires `cmake` >= 3.5 and a CUDA toolkit >= 9.0. If your current setup cannot meet these requirements we advise you to use [older versions of oxDNA](https://sourceforge.net/projects/oxdna/files/).

## Compiling oxDNA

Clone the [repo](https://github.com/lorenzo-rovigatti/oxDNA.git) or extract the oxDNA archive and then:

```
cd oxDNA         # enter the oxDNA folder
mkdir build      # create a new build folder. It is good practice to compile out-of-source
cd build
cmake ..         # here you can specify additional options, see next section
make -j4         # compile oxDNA. The -jX make option makes it compile the code in parallel by using X threads.
```

At the end of the compilation three executables (oxDNA, DNAnalysis and confGenerator) will be placed in the `build/bin` directory. 

Compiling with Python bindings will also generate an `oxpy` package in the `build/oxpy` directory that can be imported in Python. Running `make install` will attempt to copy the package to the `pip`'s module directory. The specific location will depend on your system's settings. We advise you to use [virtual environments](https://docs.python.org/3/tutorial/venv.html) (see *e.g.* [pipenv](https://docs.pipenv.org/)) to avoid conflicts with other packages and/or dependency and permission issues.

### cmake options

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
* `-DOxpySystemInstall=On` By default `oxpy` is installed in the current user's home directory. By enabling this option `oxpy` will be installed as a system-wide package. It may require superuser privileges.

### make targets

* `make` compiles oxDNA
* `make docs` Produces html doxygen documentation for oxDNA (`DOCS/html_oxDNA/index.html`) and for the UTILS folder (`DOCS/html_UTILS/index.html`)
* `make rovigatti` Compiles the observables and interactions in contrib/rovigatti
* `make romano` Compiles the observables and interactions in contrib/romano
* `make install` Copies the `oxpy` package to the Python's package folder
	
## Testing

* `make test_run` runs quick tests to check whether oxDNA has been correctly compiled or not.	
* `make test_quick` runs longer tests to check that oxDNA works.
* `make test_oxpy` checks that the Python bindings work.
* `make test` runs both sets of tests above.
