# REMD implementation using oxpy

## Script requirements

- Please use the latest version of [oxDNA](https://github.com/lorenzo-rovigatti/oxDNA).
- The script requires
	- [mpi4py](https://pypi.org/project/mpi4py/) to parallelise the computation
	- [tqdm](https://github.com/tqdm/tqdm) to estimate the run time (optional package: if not found no progress meter will be shown and a warning will be issued)

## Input file requirements
- Regular oxDNA input files should work, but note that the following parameters get overwritten by oxpy:   

```
    lastconf_file 
    trajectory_file 
    energy_file 
    CUDA_device
    T
```

- For CUDA runs the script assumes that the number of available GPUs is at least equal to the number of
temperatures you are scanning and that the devices are available sequentially.
- The following parameters are specific to a REMD run.

```
    pt_temp_list =  290K, 315K, 330K, 345K
    pt_move_every = 5000 
```

where:

- `pt_temp_list` is the list of sampled temperatures.
- `pt_move_every` is number of iterations between the temperature exchange tries.


## REMD run 
To run a REMD simulation execute the following commands:

```shell
./clean.sh

mpirun -n 4  python3  remd.py  input_md

python3 reshuffle.py input_md history.json reshuffled
```

- `clean.sh` cleans up the simulation folder.
- `remd.py`  runs the simulation.
- `reshuffle.py` restores the temperature trajectories.

Note that `-n 4` corresponds to the number of temperatures you are sampling.
The `remd.py` script prints out the estimated run time and the exchange rates.
