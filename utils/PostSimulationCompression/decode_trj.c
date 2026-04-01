// ============================================================================
// Victor UltraDecode v3 - Decoder for VTJ1 format
// Matches encoder with:
// - fixed-size header:
//     magic "VTJ1"
//     version (u32)
//     scale_digits (u32)
//     n_rows, n_cols (u32 each)
//     n_b (u32)
//     MAX_BVALS * int64 for b_global_q slots
//     t0_q (int64), dt_q (int64)
//     n_E (u32)
// - global t0/dt (int64 fixed-point)
// - global box b[] (int64 fixed-point)
// - delta-coded energies (int64)
// - hybrid body predictor ONLY (NO IWT 5/3)
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define MAX_BVALS 16  // must match encoder

// ============================================================================
// Little-endian readers
// ============================================================================
static int read_u32_le(FILE *fp, uint32_t *out) {
    uint8_t b[4];
    if (fread(b,1,4,fp) != 4) return -1;
    *out = (uint32_t)b[0]
         | ((uint32_t)b[1] << 8)
         | ((uint32_t)b[2] << 16)
         | ((uint32_t)b[3] << 24);
    return 0;
}

static int read_i64_le(FILE *fp, int64_t *out) {
    uint8_t b[8];
    if (fread(b,1,8,fp) != 8) return -1;
    *out = (int64_t)(
          ((uint64_t)b[0])
        | ((uint64_t)b[1] << 8)
        | ((uint64_t)b[2] << 16)
        | ((uint64_t)b[3] << 24)
        | ((uint64_t)b[4] << 32)
        | ((uint64_t)b[5] << 40)
        | ((uint64_t)b[6] << 48)
        | ((uint64_t)b[7] << 56)
    );
    return 0;
}

// ============================================================================
// Signed varint (LEB128) reader
// ============================================================================
static int read_svarint(FILE *fp, int64_t *out) {
    uint64_t result = 0;
    int shift = 0;

    while (1) {
        int c = fgetc(fp);
        if (c == EOF) {
            return -1;  // EOF
        }
        uint8_t byte = (uint8_t)c;
        result |= (uint64_t)(byte & 0x7F) << shift;
        if (!(byte & 0x80))
            break;
        shift += 7;
        if (shift > 63) {
            fprintf(stderr, "Error: varint too long\n");
            exit(1);
        }
    }

    // Zigzag decode
    int64_t value = (int64_t)((result >> 1) ^ (uint64_t)-(int64_t)(result & 1));
    *out = value;
    return 0;
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.bin output.dat\n", argv[0]);
        return 1;
    }

    const char *inpath  = argv[1];
    const char *outpath = argv[2];

    FILE *fin = fopen(inpath, "rb");
    if (!fin) { perror("open input"); return 1; }

    FILE *fout = fopen(outpath, "w");
    if (!fout) { perror("open output"); fclose(fin); return 1; }

    // ------------------------------------------------------------------------
    // Read header
    // ------------------------------------------------------------------------
    char magic[4];
    if (fread(magic,1,4,fin) != 4) {
        fprintf(stderr,"Error: cannot read magic\n");
        return 1;
    }
    if (memcmp(magic,"VTJ1",4) != 0) {
        fprintf(stderr,"Error: bad magic, expected VTJ1\n");
        return 1;
    }

    uint32_t version, scale_digits;
    uint32_t n_rows, n_cols, n_b;
    uint32_t n_E;

    if (read_u32_le(fin, &version) != 0) {
        fprintf(stderr,"Error: cannot read version\n");
        return 1;
    }
    if (read_u32_le(fin, &scale_digits) != 0) {
        fprintf(stderr,"Error: cannot read scale_digits\n");
        return 1;
    }
    if (read_u32_le(fin, &n_rows) != 0) {
        fprintf(stderr,"Error: cannot read n_rows\n");
        return 1;
    }
    if (read_u32_le(fin, &n_cols) != 0) {
        fprintf(stderr,"Error: cannot read n_cols\n");
        return 1;
    }
    if (read_u32_le(fin, &n_b) != 0) {
        fprintf(stderr,"Error: cannot read n_b\n");
        return 1;
    }

    if (n_b > MAX_BVALS) {
        fprintf(stderr,"Error: n_b too large (%u), MAX_BVALS=%d\n", n_b, MAX_BVALS);
        return 1;
    }

    // Read MAX_BVALS box slots; only first n_b are used
    int64_t b_q[MAX_BVALS];
    for (uint32_t i = 0; i < MAX_BVALS; i++) {
        int64_t tmp;
        if (read_i64_le(fin, &tmp) != 0) {
            fprintf(stderr,"Error: cannot read b slot %u\n", i);
            return 1;
        }
        if (i < n_b) {
            b_q[i] = tmp;
        }
    }

    int64_t t0_q, dt_q;
    if (read_i64_le(fin, &t0_q) != 0) {
        fprintf(stderr,"Error: cannot read t0_q\n");
        return 1;
    }
    if (read_i64_le(fin, &dt_q) != 0) {
        fprintf(stderr,"Error: cannot read dt_q\n");
        return 1;
    }

    if (read_u32_le(fin, &n_E) != 0) {
        fprintf(stderr,"Error: cannot read n_E\n");
        return 1;
    }
    if (n_E > 16) {
        fprintf(stderr,"Error: n_E too large (%u)\n", n_E);
        return 1;
    }

    double scale_factor = pow(10.0, (double)scale_digits);

    fprintf(stderr,
            "Header: version=%u, scale_digits=%u, n_rows=%u, n_cols=%u, n_b=%u, n_E=%u\n",
            version, scale_digits, n_rows, n_cols, n_b, n_E);

    // Precompute double box
    double b_vals[MAX_BVALS];
    for (uint32_t i=0; i<n_b; i++) {
        b_vals[i] = (double)b_q[i] / scale_factor;
    }

    // ------------------------------------------------------------------------
    // Allocate buffers
    // ------------------------------------------------------------------------
    size_t total_cells = (size_t)n_rows * (size_t)n_cols;

    int64_t *q_prev2 = (int64_t*)calloc(total_cells, sizeof(int64_t));
    int64_t *q_prev  = (int64_t*)calloc(total_cells, sizeof(int64_t));
    int64_t *q_curr  = (int64_t*)calloc(total_cells, sizeof(int64_t));

    if (!q_prev2 || !q_prev || !q_curr) {
        fprintf(stderr,"Error: out of memory for q buffers\n");
        return 1;
    }

    int64_t *col_buf = (int64_t*)malloc(sizeof(int64_t)*n_rows);
    if (!col_buf) {
        fprintf(stderr,"Error: out of memory for col_buf\n");
        return 1;
    }

    int64_t *E_q   = (int64_t*)calloc(n_E, sizeof(int64_t));
    int64_t *E_tmp = (int64_t*)malloc(sizeof(int64_t)*n_E);
    if (!E_q || !E_tmp) {
        fprintf(stderr,"Error: out of memory for energy buffers\n");
        return 1;
    }

    // ------------------------------------------------------------------------
    // Frame decode loop
    // ------------------------------------------------------------------------
    int frame_index = 0;
    while (1) {
        // --------------------------------------------------------------------
        // Read energies; break cleanly if hit EOF before a full frame
        // --------------------------------------------------------------------
        int eof = 0;
        for (uint32_t i=0; i<n_E; i++) {
            int64_t v;
            if (read_svarint(fin, &v) != 0) {
                eof = 1;
                break;
            }
            E_tmp[i] = v;
        }
        if (eof) break;  // no more complete frames

        if (frame_index == 0) {
            // frame 0: absolute energies
            for (uint32_t i=0; i<n_E; i++)
                E_q[i] = E_tmp[i];
        } else {
            // frame k>0: deltas
            for (uint32_t i=0; i<n_E; i++)
                E_q[i] += E_tmp[i];
        }

        // --------------------------------------------------------------------
        // Read body columns, undo hybrid predictor (no wavelet)
        // --------------------------------------------------------------------
        for (uint32_t c=0; c<n_cols; c++) {
            // Read predictor-stream column directly (no transform)
            for (uint32_t r=0; r<n_rows; r++) {
                int64_t v;
                if (read_svarint(fin, &v) != 0) {
                    fprintf(stderr,"Error: unexpected EOF in body\n");
                    goto done;
                }
                col_buf[r] = v;
            }

            // Undo predictor into q_curr
            for (uint32_t r=0; r<n_rows; r++) {
                size_t idx = (size_t)r * (size_t)n_cols + (size_t)c;
                int64_t val = col_buf[r];

                if (frame_index == 0) {
                    // frame 0: absolute q
                    q_curr[idx] = val;
                } else if (frame_index == 1) {
                    // frame 1: first-order delta → q1 = delta + q0
                    q_curr[idx] = val + q_prev[idx];
                } else {
                    // frame >=2: second-order delta
                    // q = delta + (2*q_prev - q_prev2)
                    int64_t pred = (q_prev[idx] << 1) - q_prev2[idx];
                    q_curr[idx] = val + pred;
                }
            }
        }

        // --------------------------------------------------------------------
        // Write frame out as text
        // --------------------------------------------------------------------
        // t = ...
        int64_t t_q = t0_q + (int64_t)frame_index * dt_q;
        double t_val = (double)t_q / scale_factor;
        fprintf(fout, "t = %.14f\n", t_val);

        // b = ...
        fprintf(fout, "b =");
        for (uint32_t i=0; i<n_b; i++) {
            fprintf(fout, " %.14f", b_vals[i]);
        }
        fprintf(fout, "\n");

        // E = ...
        fprintf(fout, "E =");
        for (uint32_t i=0; i<n_E; i++) {
            double e_val = (double)E_q[i] / scale_factor;
            fprintf(fout, " %.14f", e_val);
        }
        fprintf(fout, "\n");

        // body lines
        for (uint32_t r=0; r<n_rows; r++) {
            for (uint32_t c=0; c<n_cols; c++) {
                size_t idx = (size_t)r * (size_t)n_cols + (size_t)c;
                double v = (double)q_curr[idx] / scale_factor;
                if (c == 0)
                    fprintf(fout, "%.14f", v);
                else
                    fprintf(fout, " %.14f", v);
            }
            fprintf(fout, "\n");
        }

        fprintf(fout, "\n"); // blank line between frames

        // rotate q buffers for next frame
        memcpy(q_prev2, q_prev, sizeof(int64_t)*total_cells);
        memcpy(q_prev,  q_curr, sizeof(int64_t)*total_cells);

        frame_index++;
    }

done:
    fprintf(stderr, "Decoded %d frame(s).\n", frame_index);

    free(q_prev2);
    free(q_prev);
    free(q_curr);
    free(col_buf);
    free(E_q);
    free(E_tmp);

    fclose(fin);
    fclose(fout);
    return 0;
}
