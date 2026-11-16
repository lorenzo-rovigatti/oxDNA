# Testing suite

## oxDNA tests
The core test suite is run by `test/TestSuite.py` and covers the core force fields (DNA, RNA, NA, ANM, LJ) across both MC and MD engines, external forces, thermostats (Brownian, Langevin, Bussi), and oxpy. CUDA tests are in a separate folder as not all dev machines have CUDA installed. 

### test targets

The full oxDNA test suite can be run with 
```shell
make test_quick
```

A quick smoke test can be run with
```shell
make test_run
```

oxDNA, oxpy and oat tests can all be run with
```shell
make test
```

`test_quick` takes about 5.5 minutes to run (M1 pro, 32 GB, MacOS 26.6, Nov 2026).

## CUDA tests
The CUDA-enabled tests compare that the implementations of forcefields match between the CPU and GPU as well as between the mixed and double builds (if compiled with `-DCUDA_DOUBLE=ON`) and with/without `use_edge=true`. They can either be run with PyTest.
```shell
OXDNA_BIN=/path/to/oxDNA pytest test/CUDA/
```

Or individually with e.g:
```shell
python test_cuda_na.py /path/to/oxDNA
```

See the README files in individual test folders for specific tests and expected outcomes.

## OAT tests
OAT's tests use PyTest and are found in `analysis/tests/test_cli/`. The full test suite can be run with
```shell
pytest analysis/tests/test_cli/
```cd

And a specific test with e.g:
```shell
pytest analysis/tests/test_cli/test_mean.py
```

The full test suite takes about 2.5 minutes to run (MacOS 26.6, M1 pro chip, 32 GB mem, Nov 2026)
