# Minimal oxpy-based Forward Flux Sampling example

This example provides a minimal implementation of a Forward Flux Sampling (FFS) workflow based on the oxpy interface. It demonstrates how to combine distance- and bond-based interfaces using a generic pair-based hybrid order parameter and can serve as a starting point for more complex rare-event calculations.

## Workflow

The workflow consists of one FLUX stage followed by four SHOOT stages. One representative successful configuration obtained at each interface is included in the example and illustrated below.

![Representative configurations](scheme.png)

The interfaces used in this example are

```text
FLUX                    interface = [-4, -3]
SHOOT_01_dist_m3_m2     interface = [-3, -2]
SHOOT_02_dist_m2_m1     interface = [-2, -1]
SHOOT_03_dist_m1_0      interface = [-1, 0]
SHOOT_04_bond_0_1       interface = [0, 1]
```

Each SHOOT stage samples starting configurations from the `success_*.dat` files produced by the previous stage.

## Reference results

Running the example with the provided input files yields the following reference results:

| Quantity | Value |
|----------|------:|
| Initial flux | (1.478 ± 0.023) × 10⁻⁶ timestep⁻¹ |
| P(-3→-2) | 0.01942 ± 0.00062 |
| P(-2→-1) | 0.1687 ± 0.0026 |
| P(-1→0) | 0.3074 ± 0.0033 |
| P(0→1) | 0.001535 ± 0.000088 |
| Final rate | (2.285 ± 0.159) × 10⁻¹² timestep⁻¹ |
| Average waiting time | (4.38 ± 0.31) × 10¹¹ timesteps |

The complete numerical results are provided in `ffs_reference_results.dat`.

The script

```bash
python3 plot_results.py
```

reproduces these values and generates the diagnostic plot

![Reference cumulative rate estimate](ffs_cumulative_rate.png)

showing how the final FFS rate estimate is constructed from the successive interfaces.

## Hybrid order parameter

The order parameter is distance-based while none of the monitored pairs are hydrogen bonded. Once one or more monitored pairs become hydrogen bonded (`use_bonds = true`), the order parameter switches to the number of formed hydrogen bonds.

With

```toml
pairs = [[20, 71]]
distance_thresholds = [8.0, 5.0, 3.0, 2.0]
distance_final_state = 0
use_bonds = true
bond_threshold = -0.10
```

the states are

```text
Q = -4   d > 8.0 and nbonds = 0
Q = -3   5.0 < d <= 8.0 and nbonds = 0
Q = -2   3.0 < d <= 5.0 and nbonds = 0
Q = -1   2.0 < d <= 3.0 and nbonds = 0
Q =  0   d <= 2.0 and nbonds = 0
Q = nbonds if use_bonds = true and nbonds > 0
```

The thresholds are ordered from the outermost to the innermost distance interface. Once all monitored pairs are closer than the last threshold, the system reaches `distance_final_state`.

For hybrid distance/bond FFS calculations it is generally recommended to set

```toml
distance_final_state = 0
```

so that bond formation naturally starts from

```text
[0, 1]
```

and subsequent interfaces become

```text
[1, 2], [2, 3], ...
```

where the state number is simply the number of formed hydrogen bonds.

## Running the example

To execute the complete workflow

```bash
cd example
bash run.sh
```

Once the workflow has completed, generate the reference analysis with

```bash
python3 plot_results.py
```

The provided input files are intended only as a minimal example and should normally be replaced with those corresponding to the system of interest.

## Distance-only usage

Set

```toml
use_bonds = false
```

and the same order parameter becomes purely distance-based.

## Multiple monitored pairs

The list

```toml
pairs = [
    [20, 71],
    [21, 70],
    [22, 69]
]
```

defines multiple candidate pairs.

While no listed pairs are hydrogen bonded, the order parameter uses the minimum distance among all listed pairs. Once one or more listed pairs become bonded (`use_bonds = true`), the order parameter becomes equal to the number of bonded pairs.

For example,

```text
Q = 0    no listed pairs are bonded
Q = 1    one listed pair is bonded
Q = 2    two listed pairs are bonded
Q = 3    three listed pairs are bonded
```

## Output files

Each stage produces

- `progress.json`, storing the current number of attempts and successes and allowing interrupted simulations to be resumed;
- `success_*.dat`, containing configurations that successfully crossed the target interface;
- `success_metadata.jsonl`, containing metadata associated with each successful crossing (step number, process ID, random seed, starting configuration and OP values);
- `ffs.log`, containing detailed information about crossings and failures.

Shooting stages additionally write

- `conf_statistics.log`, reporting how often each starting configuration is selected and how many successful trajectories originate from it.

The example also includes

- `ffs_reference_results.dat`, containing the reference flux, crossing probabilities, final rate and associated uncertainties;
- `ffs_cumulative_rate.png`, providing a compact diagnostic plot of the cumulative FFS rate estimate;
- one representative successful configuration for each interface, which can be visualised with oxView.

## Notes

All interface definitions and order parameter parameters are specified in the corresponding `ffs.toml` files. Consequently, the number and position of interfaces, the monitored pairs, and the transition from distance-based to bond-based sampling can all be modified without changing the Python scripts.
