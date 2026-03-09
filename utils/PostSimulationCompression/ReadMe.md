# VTJ1 Trajectory Compression Format

## UltraEncode / UltraDecode

This repository contains a **high‑performance trajectory compression
system** designed for large molecular or particle simulations such as
oxDNA.\
The system converts large plaintext trajectory files into a compact
binary representation while preserving extremely high numerical
precision.

The compression pipeline consists of two programs:

-   **encode_trj.c** --- converts plaintext trajectories (`.dat`) to
    compressed binary (`.bin`)
-   **decode_trj.c** --- reconstructs the original trajectory text from
    the binary format

The format used is called **VTJ1 (Trajectory v1)**.

------------------------------------------------------------------------

# Overview

Simulation trajectory files can be extremely large because they contain:

-   floating point coordinates
-   time values
-   box vectors
-   energies
-   particle state values
-   repeated structure across frames

Typical plaintext trajectory files store numbers as ASCII text, which is
inefficient.

The VTJ1 format reduces file size by combining:

1.  **Fixed‑point quantization**
2.  **Delta compression**
3.  **Predictive coding**
4.  **Variable length integer encoding**
5.  **Global metadata storage**

The result is a compact binary format optimized for trajectories where
particle motion between frames is small.

------------------------------------------------------------------------

# Compression Pipeline

The encoder performs the following steps:

1.  Parse the plaintext trajectory file.
2.  Convert floating‑point values to fixed‑point integers.
3.  Predict the next frame using previous frames.
4.  Store only the difference from the prediction.
5.  Encode the integers using variable‑length encoding.
6.  Store global metadata in a compact header.

The decoder reverses these steps exactly.

------------------------------------------------------------------------

# Fixed‑Point Quantization

All floating point values are converted to integers using a fixed
decimal precision.

    q = round(value × 10^scale_digits)

Example:

  Original           scale_digits   Stored integer
  ------------------ -------------- ----------------
  0.123456           6              123456
  0.12345678901234   14             12345678901234

This allows floating point values to be stored as integers while
preserving controlled precision.

Typical settings:

  scale_digits   Precision
  -------------- -----------
  6              1e‑6
  10             1e‑10
  14             1e‑14

Higher values increase accuracy but also increase encoded size slightly.

------------------------------------------------------------------------

# Frame Prediction

Particle values change slowly between simulation frames.

To exploit this, the encoder predicts values using previous frames.

## Frame 0

The first frame is stored directly.

    stored = q

## Frame 1

The second frame stores the difference from frame 0.

    delta = q1 − q0

## Frame ≥2

Later frames use **second‑order prediction**.

    prediction = 2*q_prev − q_prev2
    delta = q − prediction

This captures constant velocity motion and reduces the magnitude of
stored values.

------------------------------------------------------------------------

# Energy Compression

Energy values are also delta encoded.

Frame 0 stores absolute values:

    E0

Subsequent frames store differences:

    delta = Ek − Ek−1

Because energies change slowly, these deltas are small integers.

------------------------------------------------------------------------

# Variable‑Length Integer Encoding

All integers are written using **signed LEB128 variable length
encoding**.

Small integers require fewer bytes.

Typical sizes:

  Value Range     Storage
  --------------- ------------
  -64..63         1 byte
  -8192..8191     2 bytes
  larger values   more bytes

Because prediction keeps values small, most stored numbers occupy **1--2
bytes**.

------------------------------------------------------------------------

# Global Metadata

Information that does not change between frames is stored once in the
header.

These values include:

-   simulation box vectors
-   timestep spacing
-   number of rows
-   number of columns
-   energy dimension
-   quantization scale

This avoids repeating identical information in every frame.

------------------------------------------------------------------------

# VTJ1 Binary Layout

The file begins with a fixed‑size header.

    magic:        4 bytes   "VTJ1"
    version:      uint32
    scale_digits: uint32
    n_rows:       uint32
    n_cols:       uint32
    n_b:          uint32
    box_values:   MAX_BVALS × int64
    t0_q:         int64
    dt_q:         int64
    n_E:          uint32

After the header, the compressed frames follow.

Each frame contains:

    energy_deltas
    body_deltas (column by column)

All values are stored as variable‑length integers.

------------------------------------------------------------------------

# Decoding Pipeline

The decoder performs the inverse process:

1.  Read the VTJ1 header
2.  Allocate memory buffers
3.  Read varints
4.  Reconstruct energies using cumulative deltas
5.  Reconstruct body values using predictor
6.  Convert fixed‑point integers back to floating point

```{=html}
<!-- -->
```
    value = integer / 10^scale_digits

Frames are then written back to plaintext format.

------------------------------------------------------------------------

# Building the Programs

Compile using GCC or Clang.

## Encoder

    gcc -O3 -march=native -ffast-math -funroll-loops -std=c11 encode_trj.c -o encode_trj

## Decoder

    gcc -O3 -march=native -std=c11 decode_trj.c -o decode_trj

------------------------------------------------------------------------

# Usage

## Encode

    ./encode_trj input.dat output.bin n_cols scale_digits

Arguments:

  Parameter      Description
  -------------- --------------------------
  input.dat      plaintext trajectory
  output.bin     compressed binary
  n_cols         columns per particle row
  scale_digits   decimal precision

Example:

    ./encode_trj traj.dat traj.bin 15 14

------------------------------------------------------------------------

## Decode

    ./decode_trj input.bin output.dat

Example:

    ./decode_trj traj.bin traj_reconstructed.dat

------------------------------------------------------------------------

# Precision and Losslessness

The codec is **lossless with respect to the quantized integer
representation**.

If the original data contains at most `scale_digits` decimal digits,
reconstruction will match the original values exactly.

Otherwise the error is bounded by:

    ±0.5 × 10^(-scale_digits)

With `scale_digits = 14`, this error is extremely small.

------------------------------------------------------------------------

# Advantages

-   very high compression ratio for trajectory data
-   fast encoding and decoding
-   minimal memory overhead
-   simple binary format
-   deterministic reconstruction
-   scalable to extremely large simulations

------------------------------------------------------------------------

# Limitations

-   box vectors must remain constant
-   column count must remain constant
-   number of energies must remain constant
-   precision limited by chosen scale_digits

------------------------------------------------------------------------

# Intended Use

The VTJ1 format is optimized for:

-   molecular dynamics trajectories
-   oxDNA simulations
-   particle systems with smooth motion
-   extremely large datasets

------------------------------------------------------------------------


# Benchmarks

## Encoding Benchmark

The figure below summarizes the performance of the trajectory encoder across several datasets of different sizes.

![Encoding Benchmark](figures/encoding.png)

The benchmarks were generated by encoding multiple oxDNA trajectory files with varying system sizes and trajectory lengths.

The table below summarizes the datasets used in the benchmark. The **Original Size (MB)** corresponds to the true size of the input trajectory before compression.

| Dataset | CSV File | Original Size (MB) | Frames | Particles |
|--------|----------|--------------------|--------|-----------|
| 2 MB file | compression_2mb.csv | 2 | 1200 | 3 |
| 800 MB file | compression_800mb.csv | 809 | 200 | 15449 |
| 1500 MB file | compression_1500mb.csv | 1553 | 300 | 19499 |
| 5.1 GB file | compression.csv | 5103 | 1000 | 19499 |
| 40.6 GB file | compression_40g.csv | 40643 | 7175 | 21649 |

The figure contains two subplots:

**Left subplot — Compression Ratio**  
Shows the compression ratio achieved for each dataset as a function of the number of digits retained during encoding. Lower precision (fewer digits) results in stronger compression, while higher precision preserves more numerical accuracy.

**Right subplot — Encoding Time per Frame**  
Shows the encoding time normalized by the number of frames in the trajectory. This normalization allows a fair comparison between datasets of different trajectory lengths and highlights how the computational cost scales with system size.

Together, these benchmarks illustrate the tradeoff between retained numerical precision, compression efficiency, and encoding cost across a range of trajectory sizes.
