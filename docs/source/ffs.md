# Forward Flux Sampling

Forward flux sampling (FFS) is a standard method for extracting the kinetic rate constants of a rare molecular process from a set of simulations. The detailed description of the method can be found in molecular simulation literature, see for example [this review](https://arxiv.org/abs/0906.4758). Briefly, FFS is best suited for a transition between two states (A and B) separated by a single barrier with no metastable intermediates. The transition from A to B is partitioned by dividing interfaces {math}`λ_0, . . . , λ_n`, corresponding to values of order parameter(s) Q. For {math}`Q < λ_0`, system is in state A, and for {math}`Q > λ_n`, it is in state B. 

In the context of kinetic processes studied with oxDNA, state A can for example correspond to unbound state, and state B can correspond to a bound hybridized state with all base-pairs formed. The FFS method starts by estimating the flux {math}`φ_0^A` of trajectories leaving from state A that cross interface {math}`λ_0`: We run a simulation and record  coordinates of our system in phase space whenever we reach interface {math}`λ_0` while coming from state A. The simulation is stopped after a desired number of crossings was generated. We divide the number of recorded states by the simulation time to obtain estimate of {math}`φ_0^A`. The rate of transition from A to B is then given by {math}`k_{AB} = φ_0^A \prod_i P (λ_{i+1}|λ_i)`,  where {math}`P (λ_{i+1}|λ_i)` is the probability that trajectory coming from A and crossing {math}`λ_i` will also cross {math}`λ_{i+1}` before returning to state A. In practice, we can estimate {math}`P (λ_{i+1}|λ_i)` by launching multiple simulations from randomly selected saved configurations at interface {math}`λ_i` and recording what fraction reaches {math}`λ_i+1` and what fraction goes back to A. We record the states that reached {math}`λ_{i+1}` and use them to start simulations to estimate {math}`P (λ_{i+2}|λ_{i+1})`. 

The FFS method requires simulation with stochastic dynamics (such as Andersen thermostat) to ensure that two simulations started from the same saved configuration will not follow an identical path. FFS hence allows us to partition the transition from A to B into several intermediate steps, where the probability estimation can be trivially parallelized by launching multiple simulations for a given interface concurrently. 

In oxDNA, the FFS approach is most often used to study the kinetics of melting or hybridization of duplexes or hairpins, or rates of more complex processes, such as strand displacement reaction. For an example of use, see e.g. [this paper](http://arxiv.org/abs/1303.3370).

The provided Python scripts in `EXAMPLE/FFS_example` directory show how to setup FFS simulation to estimate the rate of melting of a duplex. If you want to use them to study a different system, you need to adapt the defintion of interfaces in the script and `input` files for oxDNA simulation accordingly. The Python scripts support multiple core execution, but requires all CPU cores to be present on a single machine. Hence, if you are submitting these scripts as part of a CPU cluster, all CPUs have to be physically on a single computer. Similarly, if you are using the GPU version of FFS, all used GPU cards have to be on the same node. We next provide description of use of FFS to estimate the rate of association of a DNA 8-mer.


## Example system

Following the FFS terminology, we want to to compute a rate between
a stable state A (two separate strands) and a state B (DNA duplex). We have a set of interfaces, labelled
{math}`λ_i`, with {math}`i = -1, \ldots, n`, along a reaction coordinate (aka order
parameter) that fractions the configuration space of the system in a set of
states {math}`Q_i`, with {math}`i = -2, \ldots, n. Q_{-2}` belongs to state A and {math}`Q_n` belongs
to state B. There are two separate stages: the inital flux calculation (from {math}`Q_{-2}` to
{math}`Q_{0}`) and a success probability calculation (states {math}`Q_n` with n>0).

The two scripts `ffs_flux.py` and `ffs_shoot.py` are named accordingly.
The example present in the `EXAMPLES/FFS_example` directory is set up to compute the association rate of an
8-mer at 40C in a simulation box. The interfaces are set to be:

* {math}`λ_{-1}`: at least one pair of complementary bases are closer than 4.
* {math}`λ_{0}`: at least one pair of complementary bases are "bound"
* {math}`λ_{1}`: at least 8 pair of complementary bases are "bound" (full binding)




### Initial flux

The script `ffs_flux.py` runs the simulations of the initial flux. There are
some things that can be changed in the script itself, or changed via the
command line. The part that needs editing in the script is clearly marked.

The script in the FLUX directory is ready to be run. By default, it will look
for 100 forwards crossing of interface {math}`λ_{-1}`. It will create a log file
with some progress report and useful information called ffs.log. The script
can be run in parallel (use option `-c <ncpus>`). It will create a set of final
configuration named according to a pattern specified in the script itself.
Every 10 seconds, a log is printed with some information about how fast things
are proceeding.

To adapt the script to different examples, it should be enough to edit the
beginning of the script with the relevant options (interfaces, executable
path, etc.).

If one wants to collect more successes, the command-line switch `-k` sets the
initial success count, not to overwrite successes that have been already
found. The log can be repeated to the screen with `-v`.
When the script is finished, it prints out the initial flux in the last line
of the log file (and to the screen if `-v` was set).


In the directory FLUX/, running 	
```text
./ffs_flux.py -c 8 -n 10 -v
```
will run the flux waiting for 10 "successes" (forwards crossing of  interface
λ_0) using 8 cpus (actually, running 8 threads. Whether these are actual CPUs
it depends on the system) in verbose mode.
This should take about 30 seconds per success on a single CPU.


### Shooting

The two directories I0I1 and I0IF are prepared to shoot from interface
{math}`λ_{0}` to {math}`λ{1}` and from {math}`λ_{1}` to {math}`λ_{2}`, respectively.
The relevant script here is `ffs_shoot.py`. The usage is the basically the same.
The variables that have to be set vary slightly, since what the script does is
slightly different.
Initial configurations for the shooting are chosen each time at random from
all the ones found given the pattern specified in the script itself.
Also, a table is printed where there is a track of
success given the initial interface, to have an idea whether there are
configurations that never (or always) succeed.
Again, if one wants to collect more successes, the command-line switch -k sets
the initial success count, not to overwrite successes that have been already
found. The example on 8 processors should take about an hour to run.
To run the script in I0I1 and I0IF, running 

```bash
./ffs_flux.py -c 4 -n 20 -v
```
will look for 20 successes on 4 cpus in verbose mode.
Beware that to have reasonably good numbers one usually needs very many
shootings. At least 100 per interface are needed.
Also these should take about 30 seconds per success on a single CPU.


```{admonition} Example results over small sets of runs
* for the FLUX folder: average number of timesteps taken to reach a success (including possibly previous runs with the same pattern) (aka the inverse of the flux): $\approx 2.7 \times 10^8$, initial flux (includes previous runs if they were there): $\approx 4\times 10^{-7}$
* For the I0I1 folder, nsuccesses: 100, nattempts: 5941, success_prob: 0.0168322, undetermined: 0
* For the I1IF forder, nsuccesses: 111, nattempts: 1785, success_prob: 0.0621849, undetermined: 0
* Overall rate: $\approx 4 \times 10^{-9}$ (events per time step)
```
