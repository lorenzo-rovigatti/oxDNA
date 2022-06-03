# Umbrella sampling

OxDNA natively supports [umbrella sampling](https://en.wikipedia.org/wiki/Umbrella_sampling) to bias Virtual Move Monte Carlo (VMMC) simulations. This page shows how this technique can be used to calculate the melting temperature of a short RNA duplex (8-mer). The corresponding files are located in the repository in the subdirectory `examples/RNA_DUPLEX_MELT`.

```{warning}
The biasing technique presented is compatible with VMMC simulations only, and therefore it cannot be used with molecular dynamics simulations.
```

```{figure} ./images/RNA_melting.png
:height: 150px
:name: DNA-melting

A duplex and a dissociated RNA 8-mer. The VMMC simulation samples transition between duplex and single stranded states.
```

## Setting the order parameters and weights

In the example the sampling is made more efficient by biasing the simulation. Indeed, in order to obtain a correct estimate of the melting temperature, the simulation has to sample the transitions between unbonded and bonded states many times. The sampling is aided by assigning weights to a particular state, specified for this example in the file `wfile.txt`, where the first column specifies the value of the order parameter (the number of native bonds in the duplex in our case) and the second column specifies the weight assigned to the state. The larger the weight, the more likely the simulation is to visit the desired state. A file specifying the weights (`wfile.txt`) might look like 
	
```text
0 8.
1 16204
2 1882.94
3 359.746
4 52.5898
5 15.0591
6 7.21252
7 2.2498
8 2.89783
```

The weights are typically chosen first by an educated guess, and then by running the simulations and adapting the weights manually to ensure that all the relevant states are sampled. The order parameter, which in this case measures the number of native base pairs formed in an 8-mer, is specified in a file which in our example is called `op.txt`: 

```text
{
order_parameter = bond
name = all_native_bonds
pair1 = 0, 15
pair2 = 1, 14
pair3 = 2, 13
pair4 = 3, 12
pair5 = 4, 11
pair6 = 5, 10
pair7 = 6, 9 
pair8 = 7, 8 
}
```

This file specifies which base pairs count towards the order parameter. In our example, it specifies all native base pairs. The nucleotides are numbered from $0$ to $N-1$, where $N$ is the total number of nucleotides in the simulation. In our example, we have two complementary strands with eight nucleotides each, so we have sixteen nucleotides in total. The nucleotides are numbered from the 3' end to the 5' end, so the first 3' nucleotide of the first strand (with index 0) is complementary to the first 5' nucleotide on the second strand (which has index 15). The system can have between 0 and 8 bonds inclusive formed, and for each possible value of the order parameter, a corresponding weight is assigned in the `wfile.txt` (see above). 

We need to specify the path to the weight file and order parameter file in the input file, which is done by including the following options: 

```text
op_file = op.txt
weights_file = wfile.txt
```

The program will now count the amount of time the simulation spends at each value of the order parameter (*i.e.* how much time the simulation spends in states with 0 native bonds, 1 native bond, ..., 8 native bonds). In order to be able to calculate the melting temperature, one is interested in the ratio of the number states where the system is in the duplex state (*i.e.* has at least 1 bond) and the number of states when the system has 0 bonds. We also need to know this ratio for a range of temperatures, so that we can interpolate them and find when the yield is 0.5, which is the definition of melting temperature. We note that in the calculation of the yield, finite size effects have to be taken into account (see the discussion in [this paper](https://doi.org/10.1088/0953-8984/22/10/104102) for details). One possibility is to run a series of simulations, each at a different temperature. However, it is more efficient to use histogram reweighting in order to calculate the respective yields at different temperatures by extrapolating the distribution of states visited in one simulation. This is achieved by specifying 

```text
extrapolate_hist =  52C, 54C, 56C, 58C, 60C, 62C, 64C, 66C, 68C, 70C
```

in the input file. Now, the algorithm will use histogram reweighting to extrapolate the number of occupied states to the specified temperatures.

````{note}
Note that the initial configuration for this example has been generated with [oxView](https://sulcgroup.github.io/oxdna-viewer/), and that since this example concerns RNA duplexes, the input file also contains the option

	interaction_type = RNA

````

## Data evaluation and estimation of the melting temperature

The simulation, after executing 

```text
oxDNA input
```

will produce, among other files, `last_hist.dat`. Note that the values of the order parameter are printed out in the `energy.dat` file during the simulation. The contents of the `last_hist.dat` file are regularly updated during the course of the simulation (the values are updated with each simulation step, but only saved to the hard-drive as often as specified by the print_conf option). Note that the simulation needs to run long enough to properly sample the transitions between all the states. For the case of an 8-mer, one needs at least about $10^9$ iterations with the appropriate weights in order to estimate the melting temperature within 1-2 K precision. An example of the contents of the `last_hist.dat` file is 

```text
#t = 800000000; extr. Ts: 0.108383 0.10905 0.109717 0.110383 0.11105 0.111717 0.112383 0.11305 0.113717 0.114383 
0 7.82092e+09 9.77615e+08 1.65144e+10 9.03838e+09 5.04243e+09 2.86621e+09 1.65915e+09 9.77615e+08 5.86057e+08 3.57265e+08 2.21364e+08 1.39341e+08 
1 1.14272e+08 7052.11 146764 76907.7 41117.3 22416.8 12456.8 7052.11 4065.34 2385.22 1423.65 863.974 
2 6.08557e+07 32319.5 1.18642e+06 555020 264888 128915 63949.7 32319.5 16633.5 8713.56 4644.05 2517.02 
3 9.12706e+07 253708 1.56321e+07 6.5964e+06 2.83907e+06 1.24575e+06 557034 253708 117652 55524.5 26656.4 13012.5 
4 1.13301e+08 2.15444e+06 2.16511e+08 8.29258e+07 3.23788e+07 1.28829e+07 5.22111e+06 2.15444e+06 904787 386567 167955 74178.8 
5 2.53979e+08 1.68655e+07 2.68639e+09 9.39074e+08 3.34511e+08 1.21377e+08 4.48457e+07 1.68655e+07 6.45377e+06 2.51192e+06 994075 399852 
6 7.74418e+08 1.07371e+08 2.68578e+10 8.58849e+09 2.79688e+09 9.27231e+08 3.12831e+08 1.07371e+08 3.74782e+07 1.32996e+07 4.79646e+06 1.75748e+06 
7 1.03628e+09 4.60608e+08 1.75626e+11 5.16622e+10 1.54698e+10 4.71386e+09 1.46119e+09 4.60608e+08 1.4761e+08 4.80752e+07 1.59081e+07 5.34659e+06 
8 2.5904e+09 8.93909e+08 5.0908e+11 1.38277e+11 3.82229e+10 1.07489e+10 3.07421e+09 8.93909e+08 2.64187e+08 7.93341e+07 2.41997e+07 7.49618e+06 
```

Where the first line specifies the simulation time at which the file was saved and the temperatures to which the states were extrapolated (specified in simulation units; they correspond to the ones specified in the `extrapolate_hist` option). The first column specifies the value of the order parameter (from 0 to 8 in our example) and the second column specifies the number of iterations that the simulation spent in that state. Note that the simulation is biased, *i.e.* some states are more likely to be visited because they have larger weights. The third column shows the unbiased number of states for each value of the order parameter, which is obtained by dividing the number in the second column by the weight assigned to that order parameter value. The fourth and further columns correspond to the unbiased number of states extrapolated to the respective temperatures. From the ratios of the number of bonded states (1 to 8 bonds) to the number of unbonded states (0 bonds), one can calculate the yields at the respective temperatures and obtain the melting temperature by interpolating them and finding for the temperature at which the yield is 0.5. For convenience, a script is provided to do the interpolation: 

```text
python3 estimate_Tm.py  last_hist.dat
```

The output, for the example `last_hist.dat` specified above, is 

```text
52.00 0.9774084 0.8590862
54.00 0.9566703 0.8086249
56.00 0.9185418 0.7432626
58.00 0.8521952 0.6613256
60.00 0.7470068 0.5632429
62.00 0.6024042 0.4531267
64.00 0.4380124 0.3397566
66.00 0.2868094 0.2352155
68.00 0.1723566 0.1503401
70.00 0.0977175 0.0897356
## Tm = 61.149231325455936 C = 334.29923132545593 K, width = 12.566189037958857 K
```

Where the first column specifies the temperature, the second column is the yield of duplexes in the simulation box and the third column is the finite-size-effect corrected yield of duplexes. The melting temperature corresponds to the yield (with finite-size correction) equal to 0.5. The simulation box (linear) size of 20 simulation length units corresponds to a strand concentration of $3.5 \times 10^{-4}$ M. 

The last line reports the estimated melting temperature in Celsius and Kelvin, as well as the width of the melting curve.

```{note}
The VMMC script options are the same for both oxDNA and oxRNA and the following example of setting weights, order parameters and evaluation of data produced from the VMMC simulation applies equally to the DNA model as well. 
```
