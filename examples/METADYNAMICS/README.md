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

Once launched, the script runs a certain number of metadynamics simulations that can be controlled by the `--Niter` switch and defaults to 10000. After each metadynamics iteration the current bias is saved to the `bias` folder by using [`pickle.dump`](https://docs.python.org/3/library/pickle.html#pickle.dump) and therefore can be loaded back as a `numpy` array with [`pickle.load`](https://docs.python.org/3/library/pickle.html#pickle.load).

## Examples

Following Kaufhold *et al.*, the examples here use input options used by [Maffeo *et al.*](https://doi.org/10.1093/nar/gkaa200).

### `SINGLE_HELIX_BENDING` 

This is the first system investigated by Kaufhold *et al*. Here we look at the free energy as a function of end-to-end distance of a short section of duplex DNA.

The example can be run with `$ bash do.run.sh`, which will launch 6 simultaneous walkers which share a bias. The bias (as specified in the `locs.meta` file) is the distance between the centres of mass of nucleotides 30, 31, 32, 29, 28, 27 and 56, 58, 59, 0, 1, 2, that is, the centres of mass of the initial and final bits of the duplex.

The folder contains a Jupyter notebook (`analysis.ipynb`) which demonstrates how to infer free energy profiles from the simulation results. The notebook can be run simultaneously with the simulations.

## References

* Original metadyamics paper: [https://doi.org/10.1073/pnas.202427399](https://doi.org/10.1073/pnas.202427399)
* Multi-walker metadynamics paper: [https://doi.org/10.1021/jp054359r](https://doi.org/10.1021/jp054359r)
* Metadynamics parameters paper: [https://doi.org/10.1093/nar/gkaa200](https://doi.org/10.1093/nar/gkaa200)
* Kaufhold *et al.*'s code: [https://zenodo.org/record/6326800](https://zenodo.org/record/6326800)
* Kaufhold *et al.*'s paper: [https://arxiv.org/abs/2110.01477](https://arxiv.org/abs/2110.01477)
