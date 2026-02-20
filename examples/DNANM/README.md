# DNANM Examples

This folder contains DNANM/DNACT example inputs that are not part of the automated test suite.

## Contents

- `hybrid_cage/`
  - `input_cuda`: CUDA DNANM run for the hybrid cage system.
  - `ac.par`, `acc2.top`, `acc2_relaxed.dat`, `acforceponly.txt`: input data used by the hybrid cage setup.
- `KDPG/`
  - Manual DNANM/DNACT CPU/CUDA input variants and associated topology/parameter files.

## Run

From one of the subfolders:

```bash
/path/to/oxDNA input_file_name
```
