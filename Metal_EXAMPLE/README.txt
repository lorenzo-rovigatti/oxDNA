Metal example

This example runs a short GPU-enabled simulation of a system composed of 2048 double strands
(32768 nucleotides) on the Apple-GPU (Metal) backend. Note that you need to compile oxDNA
with Metal support (add the flag -DMETAL=ON to the cmake command) on macOS with an
Apple-Silicon GPU; see BUILD_METAL.md in the repository root. `input` runs the Metal
backend, `input_cpu` runs the same system on the CPU backend for comparison.
