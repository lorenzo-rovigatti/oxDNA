#!/usr/bin/env python3
"""
Repair oxDNA zstd-compressed files that contain more than one file header.

A compressed oxDNA file is an 8-byte header (magic bytes 'OXD\\x01', version, compression level, two reserved bytes)
followed by one zstd frame per configuration (or per observable output). Simulations run with
restart_step_counter = false append to existing output files, and oxDNA versions before the fix wrote a new header
every time a simulation was restarted, in the middle of the file. Readers then stop at the second header with an error
such as "Expected zstd frame magic at file offset ..., got 0x0144584f" (0x0144584f is 'OXD\\x01' read as an integer).

This script copies the file, dropping every header after the first one. The frames themselves are not decompressed:
only their headers are parsed to find where each one ends, so the zstandard library is not needed and files of any
size can be repaired. A frame cut short at the end of the file (e.g. by a simulation that was killed while writing)
is reported and dropped.

Usage:
    fix_zstd_trajectory.py <compressed_file> [output_file]
    fix_zstd_trajectory.py --in-place <compressed_file>
"""
import argparse
import os
import shutil
import struct
import sys
import tempfile

OXD_MAGIC = b'OXD\x01'
OXD_HEADER_SIZE = 8
ZSTD_MAGIC = 0xFD2FB528
# skippable frames (magic 0x184D2A50 to 0x184D2A5F) are valid zstd and are kept as they are
SKIPPABLE_MASK = 0xFFFFFFF0
SKIPPABLE_MAGIC = 0x184D2A50


def frame_size(f, start, file_size):
    """
    Size in bytes of the zstd frame that starts at ``start``, from its frame and block headers (the payload is skipped).
    Returns None if the frame runs past the end of the file.
    """
    f.seek(start)
    head = f.read(18)  # magic (4) + frame header descriptor (1) + window (1) + dictionary id (<= 4) + content size (<= 8)
    if len(head) < 8:
        return None
    magic = struct.unpack_from('<I', head, 0)[0]
    if magic & SKIPPABLE_MASK == SKIPPABLE_MAGIC:
        size = 8 + struct.unpack_from('<I', head, 4)[0]
        return size if start + size <= file_size else None
    if magic != ZSTD_MAGIC:
        raise ValueError(f"Expected a zstd frame or an oxDNA header at offset {start}, got {magic:#010x}")

    fhd = head[4]
    fcs_flag = (fhd >> 6) & 0x3
    single_segment = (fhd >> 5) & 0x1
    checksum = (fhd >> 2) & 0x1
    dict_flag = fhd & 0x3
    pos = 5
    if not single_segment:
        pos += 1
    pos += (0, 1, 2, 4)[dict_flag]
    pos += 1 if (single_segment and fcs_flag == 0) else (0, 2, 4, 8)[fcs_flag]

    pos += start
    while True:
        f.seek(pos)
        block = f.read(3)
        if len(block) < 3:
            return None
        bh = block[0] | (block[1] << 8) | (block[2] << 16)
        last = bh & 0x1
        block_type = (bh >> 1) & 0x3
        block_size = bh >> 3
        pos += 3 + (1 if block_type == 1 else block_size)  # RLE blocks store a single byte
        if last:
            break
    if checksum:
        pos += 4
    return pos - start if pos <= file_size else None


def fix(in_path, out_path):
    """Copy ``in_path`` to ``out_path`` without the repeated headers. Returns (frames, removed headers, dropped bytes)."""
    file_size = os.path.getsize(in_path)
    n_frames = n_removed = 0
    dropped = 0
    with open(in_path, 'rb') as f, open(out_path, 'wb') as out:
        header = f.read(OXD_HEADER_SIZE)
        if len(header) < OXD_HEADER_SIZE or header[:4] != OXD_MAGIC:
            raise ValueError(f"{in_path} is not a compressed oxDNA file (it does not start with {OXD_MAGIC!r})")
        out.write(header)

        pos = OXD_HEADER_SIZE
        while pos < file_size:
            f.seek(pos)
            if f.read(4) == OXD_MAGIC:
                f.seek(pos)
                extra = f.read(OXD_HEADER_SIZE)
                if extra[4] != header[4]:
                    print(f"WARNING: the header at offset {pos} has version {extra[4]}, the first one {header[4]}", file=sys.stderr)
                pos += OXD_HEADER_SIZE
                n_removed += 1
                continue

            size = frame_size(f, pos, file_size)
            if size is None:
                dropped = file_size - pos
                print(f"WARNING: the last {dropped} bytes (from offset {pos}) are an incomplete frame and have been dropped", file=sys.stderr)
                break
            f.seek(pos)
            remaining = size
            while remaining > 0:
                chunk = f.read(min(remaining, 1 << 24))
                out.write(chunk)
                remaining -= len(chunk)
            pos += size
            n_frames += 1
    return n_frames, n_removed, dropped


def main():
    parser = argparse.ArgumentParser(description="Remove the repeated headers that restarted simulations wrote into oxDNA zstd-compressed files")
    parser.add_argument('compressed_file', help="the compressed file to repair")
    parser.add_argument('output_file', nargs='?', help="where to write the repaired file (default: <compressed_file>.fixed)")
    parser.add_argument('--in-place', action='store_true', help="replace the input file with the repaired one")
    args = parser.parse_args()

    if args.in_place and args.output_file:
        parser.error("give either an output file or --in-place, not both")

    if args.in_place:
        fd, out_path = tempfile.mkstemp(prefix='.fix_zstd_', dir=os.path.dirname(os.path.abspath(args.compressed_file)))
        os.close(fd)
    else:
        out_path = args.output_file or args.compressed_file + '.fixed'

    try:
        n_frames, n_removed, dropped = fix(args.compressed_file, out_path)
    except Exception:
        if args.in_place:
            os.remove(out_path)
        raise

    if args.in_place:
        shutil.copymode(args.compressed_file, out_path)
        os.replace(out_path, args.compressed_file)
        out_path = args.compressed_file

    print(f"{n_frames} frames, {n_removed} repeated header(s) removed" + (f", {dropped} trailing bytes dropped" if dropped else ""))
    if n_removed == 0 and not dropped:
        print("The file had no repeated headers")
    print(f"Written to {out_path}")


if __name__ == '__main__':
    main()
