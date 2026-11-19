# ANNaMo smoke test

`make test_annamo` (from the build directory) builds oxDNA, confGenerator and the ANNaMo plugin (also the CUDA one in
CUDA builds) and runs `run_smoke_test.py` (python3, standard library only), which

1. generates an initial configuration of the RNA pseudoknot with confGenerator,
2. runs a short MD simulation on the CPU and checks that it finishes with finite energies,
3. in CUDA builds, if a GPU is available, runs the same simulation on CUDA and checks that the energy computed by the GPU
   kernels (`CUDA_print_energy`) agrees with the one of the CPU code.

`RNA_pseudoknot/` holds the `topology.dat` and `dHdS_matrix.dat` of the `examples/RNA_pseudoknot` folder of
<https://github.com/lorenzo-rovigatti/ANNaMo_files> (MIT license, commit 1a1e36280852), unchanged. The `input` file is
that example's, shortened and with the current name of the interaction (`ANNaMoInteraction`, formerly
`CGNucleicAcidsInteraction`).
