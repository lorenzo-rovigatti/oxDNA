# Metadynamics with oxpy

The `metad_interface.py` script provides a minimal interface to run well-tempered metadynamics simulations with the oxDNA software package.

The script is heavily based on the code developed by [Kaufhold *et al.*](https://zenodo.org/record/6326800) (see [References](#references)). However, differently from their code this interface uses oxpy, oxDNA's Python bindings, which allows for a much greater flexibility. In addition, directly accessing (and manipulating) oxDNA's data structures can sizeably improve performances, especially when simulations are run on GPUs.

Available collective variables are the distance between two clusters of nucleotides (either in one or two dimensions), or the tangent of the ratio of two distances.
Multi-walker metadynamics, P. Raiteri et al. (https://doi.org/10.1021/jp054359r), is also supported.

Please note that the interface is still under active development.

## Requirements

The interface requires the `oxpy` lib (which can be compiled and installed by issuing the `-DPython=On` option when running `cmake` and then running `make install` after the usual oxDNA compilation), pandas and numpy.

## Usage

Run the script without any arguments to see the list of supported options. The only mandatory argument is a folder which contains the oxDNA files (at least the input, topology and configuration files) required to run an unbiased simulation of the system that one wishes to study. The only other required file is a file containing the indexes of the particles whose coordinates are used to build the collective variables that will be used to bias the simulations. The name of this file, which by default is `locs.meta`, can be optionally specified with the `--p_fname` switch.

Note that rather than passing the options as arguments, it is also possible to put them in a toml file and pass it to the interface with the `--config FILENAME` switch.

All the options used to run the simulation are dumped in a `metad_PID.toml` file, where `PID` is the process ID, at the beginning of each run.

Once launched, the script runs a certain number of metadynamics simulations that can be controlled by the `--Niter` switch and defaults to 10000. After each metadynamics iteration the current bias is saved to the `bias` folder by using [`pickle.dump`](https://docs.python.org/3/library/pickle.html#pickle.dump) and therefore can be loaded back as a `numpy` array with [`pickle.load`](https://docs.python.org/3/library/pickle.html#pickle.load).

When used on CUDA-powered simulations, by default the interface launches the processes without any indication about the device that should be used, and in fact will honour the `CUDA_device` key if found in the base input file. Therefore, if GPUs are not set to run in compute mode "EXCLUSIVE_PROCESS", all the walkers will use the same GPU (the one specified by `CUDA_device` or, if this is not present in the base input file, the first available GPU). The `--use_sequential_GPUs` switch can be used to tell each walker to run on a distinct GPU: walker 0 will run on device 0, walker 1 on device 1, *etc.*

## Supported collective variables

As of now, the interface supports the following collective variables:

1. Distance between a pair of sets of nucleotides (also available on GPUs).
2. Pair of distances between sets of nucleotides. It is the only 2D order parameter available (CPU only).
3. The tangent of the ratio of two distances (CPU only).
4. The angle between the distances connecting the centres of mass of three sets of nucleotides (CPU only).

## Set up a metadynamics simulation

In order to launch a metadynamics simulation you will need

* A working oxDNA simulation of the system you want to study
* A collective variable you want to use as the order parameter

First you need to prepare a folder that contains all the input files required to launch an unbias (*i.e.* regular) oxDNA simulation. The input file **must** be called `input`. Then you need to prepare a `p_fname` file that contains the list of nucleotides that will be used to compute the collective variables. Finally, you run `metad_interface.py` with appropriate arguments.

## Prepare a "sampling" folder

Although you can use `last_beta_fe.dat` and `last_fe.dat` as the free-energy profiles, it is usually better to run one or more sampling simulations biased with the biasing potential obtained at the end of the metadynamics simulations. The `prepare_sampling.py` can be used to prepare a "sampling" folder as follows:

```
python3 /path/to/prepare_sampling.py --config metad_details.toml SAMPLING
```

where `metad_details.toml` is either the toml file you used to run the metadynamics simulations, or the toml file printed by the interface, and `SAMPLING` is the name of the sampling folder. Note that the script will also copy the `metad_details.toml` file to the new folder, renaming it to `original_metad_config.toml`. You can then step in the sampling folder and run a regular oxDNA simulations:

```
/path/to/oxDNA input
```

### Unbias the results

You can obtain the unbiased free energy by using another accompanying script, `unbias_sampling.py`, which requires as its only parameter the metad configuration file:

```
/path/to/unbias_sampling.py --config original_metad_config.toml
```

At the end of the analysis, you will find the `biased_beta_FE.dat` and `unbiased_beta_FE.dat` files. If you used the continuous coordination OP, then the script will also generate an `unbiased_hb_beta_FE.dat` file.

## Examples

### `SINGLE_HELIX_BENDING` 

This is the first system investigated by Kaufhold *et al*. Here we look at the free energy as a function of end-to-end distance of a short section of duplex DNA.

The example can be run on CPUs or GPUs with the `do.run_CPU.sh` and `do.run_CUDA.sh` scripts, respectively. Note that for the small system studied here we include a CUDA example for the sake of completeness, since using GPUs will result in much slower sampling compared to CPUs. The script will launch 6 CPU (or 4 GPU) simultaneous walkers which share a bias. The bias (as specified in the `locs.meta` file) is the distance between the centres of mass of nucleotides 30, 31, 32, 29, 28, 27 and 56, 58, 59, 0, 1, 2, that is, the centres of mass of the initial and final bits of the duplex.

### `DS6`

The system contains a hexamer (*i.e.* a double strand of length six). The order parameter is a continuous coordination number that can be used to extract free energy profiles that are function of the number of hydrogen bonds. The example can be run with `python3 path/to/metad_interface.py --config metad_details.toml base`.

## References

* Original metadyamics paper: [https://doi.org/10.1073/pnas.202427399](https://doi.org/10.1073/pnas.202427399)
* Multi-walker metadynamics paper: [https://doi.org/10.1021/jp054359r](https://doi.org/10.1021/jp054359r)
* Kaufhold *et al.*'s paper: [https://arxiv.org/abs/2110.01477](https://arxiv.org/abs/2110.01477)
* Kaufhold *et al.*'s code: [https://zenodo.org/record/6326800](https://zenodo.org/record/6326800)
