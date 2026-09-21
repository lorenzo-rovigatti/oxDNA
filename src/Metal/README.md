# Metal GPU Backend for oxDNA

This directory contains the Metal GPU implementation for oxDNA, designed to run on Apple Silicon (M-series) GPUs including the M4.

## Overview

The Metal backend provides GPU acceleration for oxDNA molecular dynamics simulations using Apple's Metal API. It mirrors the structure and functionality of the CUDA backend but is optimized for Apple's unified memory architecture and GPU compute capabilities.

## Directory Structure

```
src/Metal/
├── metal_defs.h              # Core Metal definitions and types
├── MetalUtils.h/.mm          # Utility functions for Metal operations
├── Backends/                 # Simulation backends
│   ├── MetalBaseBackend.h/.mm    # Base Metal backend
│   └── MD_MetalBackend.h/.mm     # Molecular dynamics backend
├── Interactions/             # MetalDNAInteraction (native), MetalCPUForceFallback, stubs for RNA/LJ/patchy/TEP
├── Lists/                    # MetalSimpleVerletList (GPU cells + Verlet list), MetalNoList
├── Thermostats/              # MetalBrownianThermostat
├── Shaders/                  # Metal shader kernels
│   ├── common.metal          # Common shader functions
│   ├── md_kernels.metal      # float + df64 velocity-Verlet kernels
│   ├── dna_kernels.metal     # DNA/DNA2 force + torque kernel (port of CUDA_DNA.cuh)
│   ├── list_kernels.metal    # cell fill + neighbour-list build
│   ├── thermostat_kernels.metal
│   └── df64.h                # double-float arithmetic for backend_precision = mixed
└── metal_utils/              # Utility structures
    └── MetalBox.h            # Simulation box structure
```

## Key Features

### Implemented
- ✅ Base Metal backend infrastructure, device/memory management, unified-memory buffers
- ✅ Velocity Verlet integration (float and df64), rigid-body quaternion update
- ✅ **Native DNA / DNA2 force+torque kernel** — FENE + `max_backbone_force` cap,
  bonded/non-bonded excluded volume, stacking, hydrogen bonding, cross stacking,
  coaxial stacking (oxDNA1 + oxDNA2), Debye-Hückel. Matches the CPU backend to ~1e-6.
- ✅ GPU cell list + Verlet list with a host-side skin check
- ✅ Brownian / `john` thermostat on the GPU
- ✅ CPU force fallback for every other interaction (correct, not faster than CPU)
- ✅ Three precision tiers: `float`, `mixed` (df64), `hardmixed` (CPU double integration)
- ✅ Whole-step command-buffer batching (one GPU submission per MD step)
- ✅ Optional OpenMP for the CPU-side integration loops (`-DUSE_OPENMP=ON`)

### TODO
- ⬜ Native RNA / LJ / patchy / TEP force kernels (currently CPU-fallback only)
- ⬜ Barostat / NPT, stress tensor, external forces on the GPU
- ⬜ MC / VMMC
- ⬜ Particle sorting (Hilbert curve) for neighbour-list locality
- ⬜ Overlap force evaluation with the next step's list build

## Requirements

- macOS 10.15 or later
- Apple Silicon (M1/M2/M3/M4) or Metal-capable AMD/Intel GPU
- Xcode command line tools
- CMake 3.5 or later

## Building

See `BUILD_METAL.md` in the repository root for the full instructions. In short:

```bash
cmake -S . -B build_metal -G Ninja -DMETAL=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_metal -j8
```

Optional flags:
- `-DUSE_OPENMP=ON` - multithread the CPU-side integration loops (needs libomp)
- `-DCMAKE_BUILD_TYPE=Debug` - build with debug symbols
- ~~`-DMETAL_DOUBLE=ON`~~ - **removed**: Apple GPUs have no hardware double
  precision. Select accuracy at run time with `backend_precision`.

## Usage

```
backend = Metal
Metal_avoid_cpu_calculations = 1   # 1 = native GPU kernels (DNA/DNA2), 0 = CPU force fallback
backend_precision = float          # float | mixed | hardmixed
interaction_type = DNA             # or DNA2
max_backbone_force = 10.0          # recommended
```

- `backend_precision = float` — everything float32 on the GPU. Fastest.
- `backend_precision = mixed` — positions/velocities/momenta as double-float
  (df64) pairs, df64 velocity-Verlet; forces stay float32. ~15 % slower.
- `backend_precision = hardmixed` — GPU float forces, CPU `double` integration
  over the shared buffers. Most accurate, ~2.4× slower than `float`.
- `backend_precision = double` is rejected (no GPU double on Apple).

## Architecture

### Memory Model
Metal uses a unified memory architecture on Apple Silicon, which means:
- Host and device share the same physical memory
- Data transfers are fast (pointer sharing in many cases)
- Uses `MTLResourceStorageModeShared` for CPU-GPU accessible buffers

### Compute Model
- **Threadgroups**: Equivalent to CUDA blocks
- **Threads**: Individual execution units
- **SIMD width**: 32 on Apple GPUs (similar to CUDA warps)

### Type System
- `m_number` / `m_number3` / `m_number4`: **always** `float` / `simd_float3` /
  `simd_float4`. Apple GPUs have no `double` in a shader; higher accuracy is a
  run-time choice (`backend_precision`), never a compile-time type change.
- `df64` (`Shaders/df64.h`): a `hi`+`lo` pair of `float` giving ~48-bit mantissa,
  used by the `mixed` tier for positions/velocities/momenta.
- `MetalBonds`: bond connectivity structure
- `MetalBox`: simulation box with PBC

## Implementation Details

### Kernel Execution
Kernels are dispatched using Metal compute command encoders:
1. Create command buffer
2. Create compute command encoder
3. Set pipeline state and buffers
4. Dispatch threads
5. End encoding and commit
6. Wait for completion (or async)

### Velocity Verlet Integration
The MD timestep follows the standard velocity Verlet scheme:
1. **First step**: Update velocities (half-step) and positions
2. **Force calculation**: Compute forces on new positions
3. **Second step**: Update velocities (second half-step)

### Periodic Boundary Conditions
Implemented in shaders using minimum image convention:
```metal
r -= box_sides * rint(r * inv_sides)
```

## Performance Considerations

- **Threadgroup size**: Default 256, optimal for Apple GPUs
- **Memory coalescing**: Metal handles this automatically with unified memory
- **Buffer sharing**: Minimize host-device copies using shared buffers
- **Asynchronous execution**: Can overlap compute and data transfers

## Comparison with CUDA

| Feature | CUDA | Metal |
|---------|------|-------|
| Memory model | Separate host/device | Unified (on Apple Silicon) |
| Thread organization | Grid/Block/Thread | Grid/Threadgroup/Thread |
| SIMD width | 32 (warp) | 32 (SIMD group) |
| Language | CUDA C++ | Metal Shading Language |
| Portability | NVIDIA GPUs | Apple devices |

## Contributing

When adding new features:
1. Follow the existing CUDA structure for consistency
2. Add corresponding shader implementations in `Shaders/`
3. Update pipeline creation in backend initialization
4. Test on both M-series and Intel Macs if possible

## References

- [Metal Programming Guide](https://developer.apple.com/metal/)
- [Metal Shading Language Specification](https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf)
- [Metal Best Practices Guide](https://developer.apple.com/documentation/metal/best_practices)

## License

This Metal implementation follows the same license as oxDNA.
