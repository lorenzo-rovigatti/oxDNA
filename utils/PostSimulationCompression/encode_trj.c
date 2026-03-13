// ============================================================================
// Victor UltraEncode v2 - Maximal Compression Trajectory Encoder
// Layout: global t0/dt, global box, delta-coded energies,
// hybrid body predictor, direct svarint coding, runtime scale_digits
// (NO 5/3 IWT).
// Header is fixed-size so the final rewrite does NOT clobber frame data.
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

// ============================================================================
// Config
// ============================================================================
// NOTE: scale_digits is now runtime, NOT a fixed #define.
// We still keep MAX_* as before.
#define MAX_COLS     32                 // oxDNA uses 15 but allow margin
#define MAX_BVALS    16                 // max box entries; header always reserves space for this many
#define MAX_EVALS     8                 // energies have <=3 numbers

// Runtime scale globals
static int      g_scale_digits  = 14;
static int64_t  g_scale_factor  = 100000000000000LL;  // 10^14 default

static void set_scale_digits(int sd) {
    if (sd < 0 || sd > 18) {
        fprintf(stderr, "Error: scale_digits (%d) out of range [0,18]\n", sd);
        exit(1);
    }
    g_scale_digits = sd;

    int64_t f = 1;
    for (int i = 0; i < g_scale_digits; ++i) {
        f *= 10;
    }
    g_scale_factor = f;

    fprintf(stderr,
            "Encoder: using scale_digits=%d, scale_factor=%lld\n",
            g_scale_digits, (long long)g_scale_factor);
}

// ============================================================================
// Little-endian helpers
// ============================================================================
static void write_u32_le(FILE *fp, uint32_t v) {
    uint8_t b[4];
    b[0]=(uint8_t)(v);
    b[1]=(uint8_t)(v>>8);
    b[2]=(uint8_t)(v>>16);
    b[3]=(uint8_t)(v>>24);
    fwrite(b,1,4,fp);
}

static void write_i64_le(FILE *fp, int64_t v) {
    uint8_t b[8];
    b[0] = (uint8_t)(v      );
    b[1] = (uint8_t)(v >>  8);
    b[2] = (uint8_t)(v >> 16);
    b[3] = (uint8_t)(v >> 24);
    b[4] = (uint8_t)(v >> 32);
    b[5] = (uint8_t)(v >> 40);
    b[6] = (uint8_t)(v >> 48);
    b[7] = (uint8_t)(v >> 56);
    fwrite(b,1,8,fp);
}

// ============================================================================
// Signed varint (LEB128)
// ============================================================================
static void write_svarint(FILE *fp, int64_t v) {
    // Zigzag encode
    uint64_t uv = (uint64_t)((v << 1) ^ (v >> 63));
    while (uv > 0x7F) {
        uint8_t byte = (uint8_t)((uv & 0x7F) | 0x80);
        fwrite(&byte,1,1,fp);
        uv >>= 7;
    }
    uint8_t last = (uint8_t)(uv);
    fwrite(&last,1,1,fp);
}

// ============================================================================
// Quantization
// ============================================================================
static inline int64_t q_from_double(double x) {
    double tmp = x * (double)g_scale_factor;
    return (int64_t)llrint(tmp);
}

// ============================================================================
// Parsing helpers
// ============================================================================

// Read a line, stripping trailing newline. Returns NULL on EOF.
static char *read_line(FILE *fp, char *buf, int bufsize) {
    if (!fgets(buf, bufsize, fp)) return NULL;
    size_t L = strlen(buf);
    while (L > 0 && (buf[L-1]=='\n' || buf[L-1]=='\r')) buf[--L] = 0;
    return buf;
}

// Parse "t = ..." or "b = ..." or "E = ..." lines, extracting doubles after '='.
static int parse_values_after_equal(const char *s, double *vals, int maxvals) {
    const char *eq = strchr(s, '=');
    if (!eq) return -1;
    eq++;

    int count=0;
    while (*eq && count < maxvals) {
        while (*eq && isspace((unsigned char)*eq)) eq++;
        if (!*eq) break;
        char *endp;
        double v = strtod(eq, &endp);
        if (endp == eq) break;
        vals[count++] = v;
        eq = endp;
    }
    return count;
}

// Parse a body line with N columns of doubles
static int parse_body_line(const char *s, double *cols, int n_cols) {
    const char *p = s;
    for (int i = 0; i < n_cols; i++) {
        while (*p && isspace((unsigned char)*p)) p++;
        if (!*p) return -1;
        char *endp;
        cols[i] = strtod(p, &endp);
        if (endp == p) return -1;
        p = endp;
    }
    return 0;
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr,
                "Usage: %s input.dat output.bin n_cols scale_digits\n",
                argv[0]);
        return 1;
    }

    const char *inpath  = argv[1];
    const char *outpath = argv[2];
    int n_cols = atoi(argv[3]);
    int scale_digits = atoi(argv[4]);

    if (n_cols <= 0 || n_cols > MAX_COLS) {
        fprintf(stderr, "Invalid n_cols: %d\n", n_cols);
        return 1;
    }
    set_scale_digits(scale_digits);

    FILE *fin = fopen(inpath, "r");
    if (!fin) { perror("open input"); return 1; }

    FILE *fout = fopen(outpath, "wb");
    if (!fout) { perror("open output"); fclose(fin); return 1; }

    // Buffers for parsing
    char   linebuf[4096];
    double t_vals[MAX_EVALS];
    double b_vals[MAX_BVALS];
    double E_vals[MAX_EVALS];

    int    have_global_b = 0;
    double b_global_d[MAX_BVALS];
    int    n_b_global = 0;

    int64_t b_global_q[MAX_BVALS];
    int64_t t0_q = 0, dt_q = 0;

    double t0_d = 0.0, dt_d = 0.0;

    int frame_index = 0;

    // For energy deltas
    int     have_prev_E = 0;
    int64_t prev_E_q[MAX_EVALS];
    int     n_E_global = 0;

    // Body buffers
    double  *body_prev   = NULL;   // previous frame
    double  *body_prev2  = NULL;   // second previous
    double  *body_curr   = NULL;   // current frame parsed as doubles
    int64_t *body_line   = NULL;   // integer-coded deltas for one column

    int n_rows = 0;

    // ============================================================================
    // Main file header (fixed size):
    // magic = "VTJ1"
    // version (u32)
    // scale_digits (u32)
    // n_rows, n_cols (u32 each)
    // n_b (u32)
    // MAX_BVALS * int64 for b_global_q slots
    // t0_q (int64), dt_q (int64)
    // n_E (u32)
    // ============================================================================
    long header_pos = ftell(fout);

    fwrite("VTJ1",1,4,fout);                       // magic
    write_u32_le(fout, 1);                         // version = 1
    write_u32_le(fout, (uint32_t)g_scale_digits);  // scale_digits
    write_u32_le(fout, 0);                         // n_rows placeholder
    write_u32_le(fout, (uint32_t)n_cols);          // n_cols known
    write_u32_le(fout, 0);                         // n_b placeholder

    // Reserve space for MAX_BVALS box entries
    for (int i = 0; i < MAX_BVALS; i++) {
        write_i64_le(fout, 0LL);                   // placeholder b_global_q[i]
    }

    // t0_q / dt_q placeholders
    write_i64_le(fout, 0LL);
    write_i64_le(fout, 0LL);
    write_u32_le(fout, 0);                         // n_E placeholder

    // ========================================================================
    // Frame loop
    // ========================================================================
    while (1) {
        // Read t-line
        if (!read_line(fin, linebuf, sizeof(linebuf))) break;
        if (linebuf[0] == 0) break;
        if (strncmp(linebuf, "t =", 3) != 0) {
            fprintf(stderr, "Expected t = line\n");
            break;
        }
        int nt = parse_values_after_equal(linebuf, t_vals, MAX_EVALS);
        if (nt < 1) {
            fprintf(stderr,"Bad t-line\n");
            break;
        }

        // Read b-line
        if (!read_line(fin, linebuf, sizeof(linebuf))) {
            fprintf(stderr,"Unexpected EOF in b-line.\n");
            break;
        }
        int nb = parse_values_after_equal(linebuf, b_vals, MAX_BVALS);
        if (nb < 1) {
            fprintf(stderr,"Bad b-line\n");
            break;
        }

        // First frame: capture global b[]
        if (!have_global_b) {
            for (int i=0; i<nb; i++) b_global_d[i] = b_vals[i];
            n_b_global    = nb;
            have_global_b = 1;
        } else {
            // verify identical b
            for (int i=0; i<nb; i++) {
                if (fabs(b_vals[i] - b_global_d[i]) > 1e-12) {
                    fprintf(stderr,"Box changed between frames. Not supported.\n");
                    exit(1);
                }
            }
        }

        // Read E-line
        if (!read_line(fin, linebuf, sizeof(linebuf))) {
            fprintf(stderr,"Unexpected EOF in E-line.\n");
            break;
        }
        int nE = parse_values_after_equal(linebuf, E_vals, MAX_EVALS);
        if (nE < 1) {
            fprintf(stderr,"Bad E-line\n");
            break;
        }

        // For first frame, lock number of energies
        if (frame_index == 0) {
            n_E_global = nE;
        } else {
            if (nE != n_E_global) {
                fprintf(stderr,"Energy dimension changed.\n");
                exit(1);
            }
        }

        // Now parse body lines
        // We know n_cols but we do *not* know n_rows until first frame
        int rowcount = 0;

        if (frame_index == 0) {
            // We parse the body once to determine n_rows
            long pos_after_E = ftell(fin);
            while (read_line(fin, linebuf, sizeof(linebuf))) {
                if (linebuf[0]=='t' && linebuf[1]==' ' && linebuf[2]=='=') {
                    // next frame begins
                    fseek(fin, pos_after_E, SEEK_SET);
                    break;
                }
                double cols[MAX_COLS];
                if (parse_body_line(linebuf, cols, n_cols)==0) {
                    rowcount++;
                } else {
                    fseek(fin, pos_after_E, SEEK_SET);
                    break;
                }
            }
            n_rows = rowcount;

            // allocate buffers
            size_t n_cells    = (size_t)n_rows * (size_t)n_cols;
            size_t body_bytes = n_cells * sizeof(double);
            size_t line_bytes = (size_t)n_rows * sizeof(int64_t);

            body_prev  = (double*)malloc(body_bytes);
            body_prev2 = (double*)malloc(body_bytes);
            body_curr  = (double*)malloc(body_bytes);
            body_line  = (int64_t*)malloc(line_bytes);

            if (!body_prev || !body_prev2 || !body_curr || !body_line) {
                fprintf(stderr,"Out of memory.\n");
                exit(1);
            }
            // reset file pointer for full parsing
            fseek(fin, pos_after_E, SEEK_SET);
        }

        // Now read exactly n_rows lines
        for (int r=0; r<n_rows; r++) {
            if (!read_line(fin, linebuf, sizeof(linebuf))) {
                fprintf(stderr,"Unexpected EOF inside body.\n");
                exit(1);
            }
            if (parse_body_line(linebuf, &body_curr[r*n_cols], n_cols)!=0) {
                fprintf(stderr,"Bad body line\n");
                exit(1);
            }
        }

        // Skip possible blank lines
        long savepos = ftell(fin);
        if (read_line(fin, linebuf, sizeof(linebuf))) {
            if (!(linebuf[0]=='t' && linebuf[1]==' ' && linebuf[2]=='=')) {
                // not the next frame's t-line, so ignore
            } else {
                // this is next t-line → rewind so next iteration sees it
                fseek(fin, savepos, SEEK_SET);
            }
        }

        // ===========================================================
        // Compute global t0 and dt from first two frames
        // ===========================================================
        if (frame_index == 0) {
            t0_d = t_vals[0];
        } else if (frame_index == 1) {
            dt_d = t_vals[0] - t0_d;
        }

        // ===========================================================
        // Quantize energies with delta coding
        // ===========================================================
        int64_t curr_E_q[MAX_EVALS];
        for (int i=0; i<n_E_global; i++)
            curr_E_q[i] = q_from_double(E_vals[i]);

        if (!have_prev_E) {
            // First frame: store absolute energies
            for (int i=0; i<n_E_global; i++) {
                write_svarint(fout, curr_E_q[i]);
                prev_E_q[i] = curr_E_q[i];
            }
            have_prev_E = 1;
        } else {
            // Delta energies
            for (int i=0; i<n_E_global; i++) {
                int64_t d = curr_E_q[i] - prev_E_q[i];
                write_svarint(fout, d);
                prev_E_q[i] = curr_E_q[i];
            }
        }

        // ===========================================================
        // Body compression
        // Hybrid predictor ONLY + varints (NO IWT)
        // ===========================================================
        for (int c = 0; c < n_cols; c++) {
            for (int r = 0; r < n_rows; r++) {
                double  v = body_curr[r*n_cols + c];
                int64_t q = q_from_double(v);

                if (frame_index == 0) {
                    // absolute for frame0
                    body_line[r] = q;
                } else if (frame_index == 1) {
                    // first-order delta
                    double  v_prev  = body_prev[r*n_cols + c];
                    int64_t q_prev  = q_from_double(v_prev);
                    body_line[r]    = q - q_prev;
                } else {
                    // second-order delta
                    double  v1  = body_prev[r*n_cols + c];
                    double  v2  = body_prev2[r*n_cols + c];
                    int64_t q1  = q_from_double(v1);
                    int64_t q2  = q_from_double(v2);
                    int64_t pred = (q1 << 1) - q2;
                    body_line[r] = q - pred;
                }
            }

            // Directly write predictor deltas as varints (no IWT)
            for (int r = 0; r < n_rows; r++)
                write_svarint(fout, body_line[r]);
        }

        // Shift prev frames
        if (frame_index >= 1) {
            size_t n_cells    = (size_t)n_rows * (size_t)n_cols;
            size_t body_bytes = n_cells * sizeof(double);
            memcpy(body_prev2, body_prev, body_bytes);
        }

        {
            size_t n_cells    = (size_t)n_rows * (size_t)n_cols;
            size_t body_bytes = n_cells * sizeof(double);
            memcpy(body_prev, body_curr, body_bytes);
        }

        frame_index++;
    }

    // Now that we know t0, dt, b_global, n_rows, n_E we rewrite final header:
    t0_q = q_from_double(t0_d);
    dt_q = q_from_double(dt_d);

    // Quantize b[]
    for (int i=0; i<n_b_global; i++)
        b_global_q[i] = q_from_double(b_global_d[i]);

    long endpos = ftell(fout);
    fseek(fout, header_pos, SEEK_SET);

    fwrite("VTJ1",1,4,fout);
    write_u32_le(fout, 1);                         // version
    write_u32_le(fout, (uint32_t)g_scale_digits);  // scale
    write_u32_le(fout, (uint32_t)n_rows);          // now known
    write_u32_le(fout, (uint32_t)n_cols);
    write_u32_le(fout, (uint32_t)n_b_global);

    // Write MAX_BVALS entries: real ones first, then zeros
    for (int i = 0; i < MAX_BVALS; i++) {
        if (i < n_b_global) {
            write_i64_le(fout, b_global_q[i]);
        } else {
            write_i64_le(fout, 0LL);
        }
    }

    write_i64_le(fout, t0_q);
    write_i64_le(fout, dt_q);
    write_u32_le(fout, (uint32_t)n_E_global);

    // restore file pointer
    fseek(fout, endpos, SEEK_SET);

    fclose(fin);
    fclose(fout);

    free(body_prev);
    free(body_prev2);
    free(body_curr);
    free(body_line);

    fprintf(stderr, "Finished encoding %d frames.\n", frame_index);
    return 0;
}
