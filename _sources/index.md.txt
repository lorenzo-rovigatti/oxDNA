# oxDNA

oxDNA is a simulation code that was initially conceived as an implementation of the coarse-grained DNA model introduced by [T. E. Ouldridge, J. P. K. Doye and A. A. Louis](http://dx.doi.org/10.1063/1.3552946). It has been since reworked and it is now an extensible simulation+analysis framework. It natively supports DNA, RNA, Lennard-Jones and patchy particle simulations of different kinds on both CPU and NVIDIA GPUs.

```{eval-rst}
.. toctree::
   :maxdepth: 3
   
   install.md
   usage.md
   input.md
   configurations.md
   forces.md
   observables.md
```

```{eval-rst}
.. ifconfig:: with_oxpy

Python Bindings
---------------

.. toctree::
   :maxdepth: 3
   
   oxpy/index.md
```
