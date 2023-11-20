# oxDNA

oxDNA is a simulation code that was initially conceived as an implementation of the coarse-grained DNA model introduced by [T. E. Ouldridge, J. P. K. Doye and A. A. Louis](http://dx.doi.org/10.1063/1.3552946). It has been since reworked and it is now an extensible simulation+analysis framework. It natively supports DNA, RNA, Lennard-Jones and patchy particle simulations of different kinds on both single CPU cores and NVIDIA GPUs.

The development of this software has been partially supported by the European Commission through the Marie Skłodowska−Curie Fellowship No.702298-DELTAS, and ONR grant N000142012094.

## Documentation

The documentation can be browsed [here](https://lorenzo-rovigatti.github.io/oxDNA/).

The HTML documentation can also be generated locally by running `make html` in the `docs` folder, and the resulting docs can be browsed by opening the `docs/build/html/index.html` file. Run `pip install -r docs_requirements.txt` to install the requirements.

## Installation

Installation instructions can be found in the `docs/source/install.md` file or online [here](https://lorenzo-rovigatti.github.io/oxDNA/install.html).

## Examples

The `examples` folder contains many examples showing the main features of the code. Note that the `METADYNAMICS`, `OXPY` and `OXPY_REMD` examples require `oxpy`, oxDNA's python bindings that can be compiled by setting `-DPython=ON` during the [compilation stage](https://lorenzo-rovigatti.github.io/oxDNA/install.html#cmake-options).

The `analysis/paper_examples` folder contains examples for `oxDNA_analysis_tools`, a suite of command line Python tools for performing generic structural analyses of oxDNA simulations.

## FAQ

**Q: How do I simulate oxDNA with LAMMPS?**

**A:** This repository contains the standalone oxDNA software and it is not linked to the code that powers the LAMMPS version of the oxDNA coarse-grained model. If you have any issues/enquiries about the [oxDNA LAMMPS package](https://docs.lammps.org/pair_oxdna.html) your best bet is to directly contact its authors.

**Q: Can oxDNA be run on multiple CPU cores or GPUs?**

**A:** No, oxDNA can run simulations on single cores or single GPUs only.

## Citing oxDNA

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04693/status.svg)](https://doi.org/10.21105/joss.04693)

Please cite these publications for any work that uses the oxDNA simulation package:

- for the code:
  * [E. Poppleton et al., J. Open Source Softw. 8, 4693 (2023)](https://doi.org/10.21105/joss.04693)
- for the CUDA-powered code:
  * [L. Rovigatti et al., J. Comput. Chem. 36, 1 (2015)](https://doi.org/10.1002/jcc.23763)
- for oxDNA analysis tools:
  * [E. Poppleton et al., Nucleic Acids Research e72 (2020)](https://doi.org/10.1093/nar/gkab324)
    
## Acknowledgements

oxDNA depends on a minimum number of external libraries (a c++-14-compliant standard library and NVIDIA's CUDA if the user wishes to enable it).

Internally, oxDNA uses the following libraries, which are included in the source tree:

* [ExprTk](https://www.partow.net/programming/exprtk/index.html)
* [nlohmann's JSON library](https://github.com/nlohmann/json)
* [pybind11](https://github.com/pybind/pybind11)
* [fast_double_parser](https://github.com/lemire/fast_double_parser)

As far as I know, this is compatible with their licenses. If you are a developer or a mantainer of one of these projects and you think that oxDNA does not comply with your license, please contact us.
