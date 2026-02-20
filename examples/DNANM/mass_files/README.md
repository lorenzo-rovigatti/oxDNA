# Mass File Format

Mass values are indexed by particle `btype` in both CPU and CUDA paths.

## Usage

1. Copy `mass_files/default_masses.txt` and edit masses as needed.
2. Set `massfile = <path-to-your-file>` in the input file.
3. Ensure topology parsing assigns `btype` values consistent with the mass file.

## Format

The first line is the number of entries.  
Each following line is:

`<btype> <mass>`

Example:

```text
27
0 1
1 1
...
26 1
```

## Notes

- The default mapping covers DNA and amino-acid btypes used by the ANM/DNANM/RNANM paths.
- Additional custom btypes can be appended using new integer IDs.
- Keep the file numeric only (no inline comments).
