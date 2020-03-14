CUDA example

This example runs a short GPU-enabled simulation a system composed by 2048 double strands 
(32768 nucleotides). Note that you need to compile oxDNA with CUDA support (add the flag
-DCUDA=ON to the cmake command) and have a working CUDA installation (CUDA >= 3.2 required). 
