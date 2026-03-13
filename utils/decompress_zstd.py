#!/usr/bin/env python3
import struct
import sys
import zstandard as zstd

MAGIC = b'OXD\x01'
HEADER_SIZE = 8

def decompress_oxdna(in_filename, out_filename=None):
    # Open input file
    with open(in_filename, "rb") as f:
        header = f.read(HEADER_SIZE)
        if len(header) < HEADER_SIZE:
            raise ValueError("File too small to be a valid oxDNA-compressed file.")

        # Parse header
        magic = header[0:4]
        if magic != MAGIC:
            raise ValueError(f"Invalid magic bytes: {magic} (expected {MAGIC})")

        version = header[4]
        compression_level = header[5]
        reserved1 = header[6]
        reserved2 = header[7]

        print(f"Magic OK: {magic}")
        print(f"Header version: {version}")
        print(f"Compression level: {compression_level}")
        # reserved bytes ignored

        # Prepare output
        if out_filename is None:
            out_filename = in_filename + ".decompressed"

        dctx = zstd.ZstdDecompressor()

        # The rest of the file is a continuous zstd stream
        with open(out_filename, "wb") as out:
            with dctx.stream_reader(f) as reader:
                while True:
                    chunk = reader.read(131072)  # 128 KB chunks
                    if not chunk:
                        break
                    out.write(chunk)

        print(f"Decompressed output written to: {out_filename}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: oxdna_decompress.py <compressed_file> [output_file]")
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else None
    decompress_oxdna(in_file, out_file)
    