# Improving performance

The overall performance depends on many factors, and it is therefore not possible to come up with a set of options that will always maximise performance. However, here we provide some guidance and list some tips to help you improving the efficiency of your simulations. 

First of all, note that at the end of each simulation (if the code was not compiled with the [`-DMOSIX=On` CMake option](install.md#cmake-options)) oxDNA will print information about the total running time, the time taken per simulation step (in milliseconds) and a detailed list of timings. Here is an example:

```text
INFO: Total Running Time: 10.6746 s, per step: 1.68876 ms
INFO: Timings, in seconds, by Timer (total, own, spent in children)
> SimBackend                         10.675 (100.0%)        0.023 (  0.2%)       10.651 ( 99.8%)
***> Lists                            3.170 ( 29.7%)        3.170 ( 29.7%)        0.000 (  0.0%)
***> Thermostat                       0.005 (  0.1%)        0.005 (  0.1%)        0.000 (  0.0%)
***> Forces                           7.055 ( 66.1%)        7.055 ( 66.1%)        0.000 (  0.0%)
***> First Step                       0.401 (  3.8%)        0.401 (  3.8%)        0.000 (  0.0%)
***> Observables                      0.020 (  0.2%)        0.020 (  0.2%)        0.000 (  0.0%)
```

Each row below `SimBackend` presents the time spent by doing a specific task in seconds and, between brackets, in percentage. The list of possible tasks (which depend on the specifics of the simulation) is:

|Name|Type|Task|
|-|-|-|
|Lists|All|Building of the neighbour list|
|Thermostat|MD|Thermostatting the degrees of freedom|
|First step|MD|First step of the velocity-Verlet integration|
|Forces|MD|Computing pair interactions + second step of the velocity-Verlet integration|
|Hilbert sorting|MD on CUDA|Sort the particles to optimise data access on GPUs.Used only if [`CUDA_sort` is set](input.md#cuda-options)|
|Rotations+Translations|MC|Monte Carlo rotations and translations|
|Volume Moves|MC|Monte Carlo Volume moves|
|Observables|All|Computing and printing simulation outputs (energy, configurations, distances, *etc*)|

Given a certain task, the time spent by the code is split into total time, time spent by the task itself, and time spent by the task's children (subtasks, if you will). Note that the latter two sum up to the former.

The most important information, optimization-wise, is provided by the percentage column. 

## Molecular dynamics

In molecular dynamics simulations, the two most important parameters, performance-wise, are `verlet_skin` and `dt`.

* `verlet_skin` controls the distance a particle has to move to trigger the update of the neighbour lists. Using a larger value results in fewer updates but a larger number of possibly-interacting pairs of non-bonded particles and vice versa. You can try to decrease this number if the simulation tends to spend most of its time computing pair interactions. The most common optimal value lies between 0.05 and 0.2.
* `dt` is the integration time step. Given a number of simulated time steps, the simulated time is proporional to `dt`. Since the overall performance of the simulation is only weakly dependent on the integration time step, using the largest possible `dt` value is advisable. However, using a value that is too large will result in numerical instabilities. Depending on the thermostat settings, optimal values range from 0.001 to 0.003 but your mileage may vary, so feel free to experiment with `dt` if you want to really optimise your simulation.

### GPU simulations

When running CUDA-powered simulations, the box size has a non-trivial effect on performance, and its interaction with other parameters such as `salt_concentration`, `verlet_skin`, and possibly others, make it hard to come up with a way of automatically set the best options for a given case.

Since there is no dynamic memory on GPUs, in order to avoid crashing simulations oxDNA sets the size of the cells used to build neighbouring lists so that their memory footprint is not too high. If you want to optimise performance is sometimes worth to set `cells_auto_optimisation = false` so that oxDNA uses the smallest possible cells (at the cost of memory consumption). If the resulting memory footprint can be handled by your GPU you'll probably see some (possibly large) performance gains.

## Monte Carlo

When running Monte Carlo simulations the efficiency of the sampling depends on specific moves employed during the simulation. For regular Monte Carlo and VMMC simulations, the most important options are `delta_translation` and `delta_rotation`, which set the maximum displacement for translations and rotations. Optimal values depend very much on the system at hand, so it is hard to provide some guidelines, although often values around `0.1` given decent performance. Sometimes it may be worth to set [`adjust_moves = true` (together with `equilibration_steps > 0`)](input.md#monte-carlo-options) to let the code look for optimal values.
