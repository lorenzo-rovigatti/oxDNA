# Colored Coarsed Grain OxDNA

ccg OxDNA is an modified oxDNA model to incorporate large time scale self assembly of DNA origami. This coarse grain model can be though to be in between the Jonah's cg model and patchy particle model.













# oxDNA

oxDNA is a simulation code that was initially conceived as an implementation of the coarse-grained DNA model introduced by [T. E. Ouldridge, J. P. K. Doye and A. A. Louis](http://dx.doi.org/10.1063/1.3552946). It has been since reworked and it is now an extensible simulation+analysis framework. It natively supports DNA, RNA, Lennard-Jones and patchy particle simulations of different kinds on both single CPU cores and NVIDIA GPUs.

The development of this software has been partially supported by the European Commission through the Marie Skłodowska−Curie Fellowship No.702298-DELTAS, and ONR grant N000142012094.

## Documentation

The documentation can be browsed [here](https://lorenzo-rovigatti.github.io/oxDNA/).

The HTML documentation can also be generated locally by running `make html` in the `docs` folder, and the resulting docs can be browsed by opening the `docs/build/html/index.html` file. Run `pip install -r docs_requirements.txt` to install the requirements.


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

oxDNA depends on a minimum number of external libraries (a c++-14-compliant standard library and Nvidia's CUDA if the user wishes to enable it).

Internally, oxDNA uses the following libraries, which are included in the source tree:

* [ExprTk](https://www.partow.net/programming/exprtk/index.html)
* [nlohmann's JSON library](https://github.com/nlohmann/json)
* [pybind11](https://github.com/pybind/pybind11)

As far as I know, this is compatible with their licenses. If you are a developer or a mantainer of one of these projects and you think that oxDNA does not comply with your license, please contact us.
