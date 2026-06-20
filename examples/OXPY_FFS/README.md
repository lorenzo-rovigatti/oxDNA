# Minimal oxpy FFS example

This example provides a small oxpy-based Forward Flux Sampling workflow using a generic pair-based hybrid order parameter.

The OP uses the minimum distance over the listed pairs while no bonds are present, and transitions to the number of formed bonds when one or more listed pairs become hydrogen-bonded and `use_bonds = true`.

With:

```toml
pairs = [[20, 71]]
distance_thresholds = [8.0, 5.0, 3.0, 2.0]
distance_final_state = 0
use_bonds = true
bond_threshold = -0.10
```

the states are:

```text
Q = -4   d > 8.0 and nbonds = 0
Q = -3   5.0 < d <= 8.0 and nbonds = 0
Q = -2   3.0 < d <= 5.0 and nbonds = 0
Q = -1   2.0 < d <= 3.0 and nbonds = 0
Q =  0   d <= 2.0 and nbonds = 0
Q = nbonds if use_bonds = true and nbonds > 0
```

The thresholds are ordered from the outermost to the innermost distance interface. When all pair distances are below the last threshold, the system reaches `distance_final_state`.

For hybrid distance/bond FFS calculations, the last distance state should normally be chosen as

```toml
distance_final_state = 0
```

so that bond formation naturally starts from the interface

```text
[0, 1]
```

and subsequent bond interfaces correspond to

```text
[1, 2], [2, 3], ...
```

where the state number equals the number of formed bonds.

## Example stages

The interfaces used in this example are

```text
FLUX                    interface = [-4, -3]
SHOOT_01_dist_m3_m2     interface = [-3, -2]
SHOOT_02_dist_m2_m1     interface = [-2, -1]
SHOOT_03_dist_m1_0      interface = [-1, 0]
SHOOT_04_bond_0_1       interface = [0, 1]
```

Each SHOOT stage samples starting configurations from the `success_*.dat` files produced by the previous stage.

Run with:

```bash
cd example
bash run.sh
```

The provided input files are intended only as a minimal example and should normally be replaced with the files corresponding to the system of interest.

## Distance-only usage

Set:

```toml
use_bonds = false
```

Then the same OP becomes purely distance-based.

## Bond-only usage

Although the hybrid OP can in principle be used for bond-only FFS, the intended use case is to precede bond formation by one or more distance-based stages. In practice, a short distance pre-stage ending at `Q = 0` provides a natural unbound state from which interfaces `[0,1]`, `[1,2]`, ... can be sampled.

## Multiple pairs

The list

```toml
pairs = [
    [20, 71],
    [21, 70],
    [22, 69]
]
```

defines multiple candidate pairs.

While no bonds are present, the order parameter uses the minimum distance among all listed pairs. Once one or more listed pairs become hydrogen-bonded and `use_bonds = true`, the order parameter becomes equal to the number of bonded pairs.

For example:

```text
Q = 0    no listed pairs are bonded
Q = 1    one listed pair is bonded
Q = 2    two listed pairs are bonded
Q = 3    three listed pairs are bonded
```

## Output files

Each stage writes:

* `progress.json`, which stores the current number of attempts and successes and allows interrupted simulations to be resumed;
* `success_*.dat`, containing configurations that successfully crossed the target interface;
* `success_metadata.jsonl`, containing metadata associated with each successful crossing (step number, process ID, random seed, starting configuration and OP values);
* `ffs.log`, containing detailed information about crossings and failures.

Shooting stages additionally write `conf_statistics.log`, which records how often each starting configuration is selected and how many successful trajectories originate from it.

## Notes

All interface definitions and OP parameters are specified in the corresponding `ffs.toml` files. Consequently, the number of interfaces, their positions, the monitored pairs, and the transition from distance-based to bond-based sampling can all be modified without changing the Python scripts.

When using bond formation (`use_bonds = true`), it is recommended that the final distance stage terminates at `Q = 0`, so that bond-count states correspond directly to `Q = 1, 2, ...`.
