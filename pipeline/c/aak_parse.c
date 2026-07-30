/*
 * aak_parse.c — Convert AAK order-type binary file (.b08/.b16) to a compact
 * preprocessed chirotope file (.chi) that Python can read with numpy.
 *
 * Input (.b08):  each record is n*(x,y) pairs stored as uint8.
 * Input (.b16):  each record is n*(x,y) pairs stored as uint16.
 * Output (.chi):
 *   [0..3]   uint32  magic    = 0x41414B01
 *   [4..7]   uint32  n
 *   [8..15]  uint64  num_valid_records
 *   then for each valid (non-degenerate) record:
 *     int8  [C(n,3)]   chirotope signs (+1 or -1)
 *     uint16[n*2]      original (x,y) coordinates (always uint16 in output)
 *
 * Usage:  aak_parse <input.b08|.b16> <n> <output.chi> [coord_bytes]
 *         coord_bytes defaults to 2; pass 1 for .b08 files.
 *
 * Compile:  cc -O3 -o aak_parse aak_parse.c
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

#define MAX_N       12
#define MAX_TRIPLES 220   /* C(12,3) = 220 */
#define MAGIC       0x41414B01u

static int ti[MAX_TRIPLES], tj[MAX_TRIPLES], tk[MAX_TRIPLES];
static int num_triples;

static void precompute_triples(int n) {
    num_triples = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            for (int k = j+1; k < n; k++) {
                ti[num_triples] = i;
                tj[num_triples] = j;
                tk[num_triples] = k;
                num_triples++;
            }
}

/* Write buffered output in large blocks to avoid per-record syscall overhead. */
#define OUT_BUF_RECORDS 65536

int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 5) {
        fprintf(stderr, "Usage: aak_parse <input.b08|.b16> <n> <output.chi> [coord_bytes]\n");
        return 1;
    }

    const char *inpath     = argv[1];
    int         n          = atoi(argv[2]);
    const char *outpath    = argv[3];
    int         coord_bytes = (argc == 5) ? atoi(argv[4]) : 2;

    if (n < 3 || n > MAX_N) {
        fprintf(stderr, "error: n must be 3..%d\n", MAX_N);
        return 1;
    }
    if (coord_bytes != 1 && coord_bytes != 2) {
        fprintf(stderr, "error: coord_bytes must be 1 or 2\n");
        return 1;
    }

    precompute_triples(n);

    /* ---- mmap input ---- */
    int fd = open(inpath, O_RDONLY);
    if (fd < 0) { perror("open"); return 1; }

    struct stat st;
    if (fstat(fd, &st) < 0) { perror("fstat"); return 1; }
    size_t file_size   = (size_t)st.st_size;
    size_t record_size = (size_t)n * 2 * (size_t)coord_bytes;
    size_t num_records = file_size / record_size;

    if (num_records == 0) {
        fprintf(stderr, "error: file too small or wrong n\n");
        return 1;
    }

    uint8_t *data = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (data == MAP_FAILED) { perror("mmap"); return 1; }
    close(fd);

    fprintf(stderr, "Input:  %s\n", inpath);
    fprintf(stderr, "        %zu bytes  |  %zu records  |  n=%d  |  C(n,3)=%d  |  coord_bytes=%d\n",
            file_size, num_records, n, num_triples, coord_bytes);

    /* ---- open output ---- */
    FILE *fout = fopen(outpath, "wb");
    if (!fout) { perror("fopen output"); return 1; }

    /* write header with placeholder count */
    uint32_t magic = MAGIC;
    uint32_t nn    = (uint32_t)n;
    uint64_t count = 0;
    fwrite(&magic, 4, 1, fout);
    fwrite(&nn,    4, 1, fout);
    long count_offset = ftell(fout);
    fwrite(&count, 8, 1, fout);

    /* per-record output size */
    size_t chi_bytes = (size_t)num_triples;           /* int8 */
    size_t pts_bytes = (size_t)n * 2 * sizeof(uint16_t);
    size_t row_bytes = chi_bytes + pts_bytes;

    /* output buffer */
    size_t buf_size = OUT_BUF_RECORDS * row_bytes;
    uint8_t *buf    = malloc(buf_size);
    if (!buf) { perror("malloc"); return 1; }
    size_t buf_pos  = 0;

    uint64_t valid = 0;
    int8_t   chi_tmp[MAX_TRIPLES];

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Scratch buffer to hold one record's coords as uint16 (for .chi output). */
    uint16_t pts16[MAX_N * 2];

    for (size_t rec = 0; rec < num_records; rec++) {
        const uint8_t *raw = data + rec * record_size;

        /* Unpack coords into pts16 regardless of input width. */
        if (coord_bytes == 1) {
            for (int k = 0; k < n * 2; k++)
                pts16[k] = raw[k];
        } else {
            memcpy(pts16, raw, record_size);
        }

        int degenerate = 0;
        for (int t = 0; t < num_triples; t++) {
            int32_t xi = pts16[ti[t]*2],   yi = pts16[ti[t]*2+1];
            int32_t xj = pts16[tj[t]*2],   yj = pts16[tj[t]*2+1];
            int32_t xk = pts16[tk[t]*2],   yk = pts16[tk[t]*2+1];
            int64_t det = (int64_t)(xj-xi)*(yk-yi) - (int64_t)(xk-xi)*(yj-yi);
            if (det == 0) { degenerate = 1; break; }
            chi_tmp[t] = (det > 0) ? 1 : -1;
        }
        if (degenerate) continue;

        /* append to output buffer */
        memcpy(buf + buf_pos,              chi_tmp, chi_bytes);
        memcpy(buf + buf_pos + chi_bytes,  pts16,   pts_bytes);
        buf_pos += row_bytes;
        valid++;

        /* flush when buffer is full */
        if (buf_pos >= buf_size) {
            fwrite(buf, 1, buf_pos, fout);
            buf_pos = 0;
        }

        if (rec % 1000000 == 999999) {
            struct timespec tn;
            clock_gettime(CLOCK_MONOTONIC, &tn);
            double elapsed = (tn.tv_sec - t0.tv_sec) + (tn.tv_nsec - t0.tv_nsec) * 1e-9;
            double rate    = (double)(rec + 1) / elapsed;
            double remain  = (num_records - rec - 1) / rate;
            fprintf(stderr, "  %zu/%zu  (%.1f%%)  %.0f Mrec/s  ~%.0fs remaining\n",
                    rec + 1, num_records,
                    100.0 * (rec + 1) / num_records,
                    rate / 1e6,
                    remain);
        }
    }

    /* flush remainder */
    if (buf_pos > 0)
        fwrite(buf, 1, buf_pos, fout);
    free(buf);

    /* fix up record count in header */
    fseek(fout, count_offset, SEEK_SET);
    fwrite(&valid, 8, 1, fout);
    fclose(fout);

    munmap(data, file_size);

    struct timespec tend;
    clock_gettime(CLOCK_MONOTONIC, &tend);
    double elapsed = (tend.tv_sec - t0.tv_sec) + (tend.tv_nsec - t0.tv_nsec) * 1e-9;

    fprintf(stderr, "Output: %s\n", outpath);
    fprintf(stderr, "        %lu valid records  |  %.2fs  |  %.1f Mrec/s\n",
            (unsigned long)valid, elapsed, (double)valid / elapsed / 1e6);
    return 0;
}
