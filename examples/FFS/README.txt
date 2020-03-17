
Forward Flux Sampling example

WARNING: this example requires python 2.7 and the "new" oxDNA code.

This example implements the forward flux sampling as described in the appendix
of http://arxiv.org/abs/1303.3370.
Following the terminology in there, we want to to compute a rate between
a (meta) stable state A and a state B. We have a set of interfaces, labelled
\lambda_i, with i = -1, \ldots, n, along a reaction coordinate (aka order
parameter) that fractions the configuration space of the system in a set of
states Q_i, with i = -2, \ldots, n. Q_{-2} belongs to state A and Q_n belongs
to state B.

There are two separate stages: the inital flux calculation (from Q_{-2} to
Q_{-0}) and a success probability calculation (states Q_n with n>0).

The two scripts ffs_flux.py and ffs_shoot.py are named accordingly.

EXAMPLE SYSTEM: this example is set up to compute the association rate of an
8-mer at 40C (slightly below melting). The interfaces are set to be:
\lambda_{-1}: at least one pair of complementary bases are closer than 4.
\lambda_{0}: at least one pair of complementary bases are "bound"
\lambda_{1}: at least 8 pair of complementary bases are "bound" (full binding)


INITIAL FLUX:

the script ffs_flux.py runs the simulations of the initial flux. There are
some things that can be changed in the script itself, or changed via the
command line. The part that needs editing in the script is clearly marked.

The script in the FLUX directory is ready to be run. By default, it will look
for 100 forwards crossing of interface lambda_{-1}. It will create a log file
with some progress report and useful information called ffs.log. The script
can be run in parallel (use -c <ncpus>). It will create a set of final
configuration named according to a pattern specified in the script itself.
Every 10 seconds, a log is printed with some information about how fast things
are proceeding.

To adapt the script to different examples, it should be enough to edit the
beginning of the script with the relevant options (interfaces, executable
path, etc.).

If one wants to collect more successes, the command-line switch -k sets the
initial success count, not to overwrite successes that have been already
found.

the log can be repeated to the screen with -v.

When the script is finished, it prints out the initial flux in the last line
of the log file (and to the screen if -v was set).

in the directory FLUX/, running 

./ffs_flux.py -c 8 -n 10 -v

will run the flux waiting for 10 "successes" (forwards crossing of  interface
lambda_0) using 8 cpus (actually, running 8 threads. Whether these are actual CPUs
it depends on the system) in verbose mode.
This should take about 30 seconds per success on a single CPU.


SHOOTING:

The two directories I0I1 and I0IF are prepared to shoot from interface
lambda_{0} to \lambda{1} and from lambda_{1} to lambda_{2}, respectively.

The relevant script here is ffs_shoot.py. The usage is the basically the same.
The variables that have to be set vary slightly, since what the script does is
slightly different.

Initial configurations for the shooting are chosen each time at random from
all the ones found given the pattern specified in the script itself.

The success probability Also, a table is printed where there is a track of
success given the initial interface, to have an idea whether there are
configurations that never (or always) succeed.

Again, if one wants to collect more successes, the command-line switch -k sets
the initial success count, not to overwrite successes that have been already
found.

The example on 8 processors should take about an hour to run.

in I0I1 and I0IF, running 

./ffs_flux.py -c 4 -n 20 -v

will look for 20 successe on 4 cpus in verbose mode.

Beware that to have reasonably good numbers one usually needs very many
shootings. At least 100 per interface are needed.
Also these should take about 30 seconds per success on a single CPU.




Example results over small sets of runs:

FLUX:
Main: average number of timesteps taken to reach a success (including possibly
previous runs with the same pattern) (aka 1./flux): 2.71181e+08
Main: initial flux (includes previous runs if they were there): ~ 4e-7

I0I1
 nsuccesses: 100 nattempts: 5941 success_prob: 0.0168322 undetermined: 0
I1IF
 nsuccesses: 111 nattempts: 1785 success_prob: 0.0621849 undetermined: 0

rate: ~ 4. * 10^-9 (events per time step)

