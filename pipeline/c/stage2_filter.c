/*
 * stage2_filter.c — Stage 2 Gale-diagram filter in C.
 *
 * Reads a .chi file (produced by aak_parse), applies the full Stage 2
 * pipeline (positive-pair selection, singleton face check, missing-face
 * computation up to size 5, p≥2 filter) using the chirotope-based face
 * criterion.  All geometric tests are pure integer chirotope lookups —
 * no floating-point, no coordinates needed.
 *
 * Outputs a .gd2 file containing one entry per unique "raw" missing-face
 * pattern (sorted face-bitmask key, no label-invariant canonicalisation —
 * Python handles that).
 *
 * Usage:  stage2_filter <input.chi> <d> <output.gd2>
 *
 * .gd2 file layout:
 *   uint32  magic = 0x47443200
 *   uint32  n
 *   uint64  num_unique_patterns
 *   then for each pattern:
 *     uint32  record_id
 *     uint8   u, v        (positive pair)
 *     uint8   num_mf
 *     uint16[num_mf]  face bitmasks (sorted ascending)
 *
 * Face criterion (Ziegler, Lectures on Polytopes §6):
 *   S ⊆ [n] is a face iff in T=[n]\S every positive point of T lies in
 *   conv(negative points of T).
 *
 * Chirotope interior test (non-degenerate data):
 *   p strictly inside conv(S) iff
 *     NOT ∃ sᵢ∈S : ∀ sⱼ∈S\{sᵢ} : chi(p,sᵢ,sⱼ)=+1
 *
 * Compile:  cc -O3 -o stage2_filter stage2_filter.c
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
#ifdef _OPENMP
#  include <omp.h>
#endif

#define MAX_N        10
#define MAX_MF       64        /* max missing faces tracked per diagram */
#define MAGIC_IN     0x41414B01u
#define MAGIC_OUT    0x47443200u

/* -------------------------------------------------------------------------
 * Chirotope lookup tables
 * ------------------------------------------------------------------------- */

static uint8_t g_chi_idx[MAX_N][MAX_N][MAX_N];
static int8_t  g_chi_sgn[MAX_N][MAX_N][MAX_N];

static void init_chi_tables(int n) {
    int idx = 0, tbl[MAX_N][MAX_N][MAX_N];
    memset(tbl, 0, sizeof(tbl));
    for (int a = 0; a < n; a++)
        for (int b = a+1; b < n; b++)
            for (int c = b+1; c < n; c++)
                tbl[a][b][c] = idx++;

    for (int p = 0; p < n; p++)
    for (int q = 0; q < n; q++)
    for (int r = 0; r < n; r++) {
        if (p==q || q==r || p==r) {
            g_chi_idx[p][q][r] = 0;
            g_chi_sgn[p][q][r] = 0;
            continue;
        }
        int a=p, b=q, c=r, sw=0, tmp;
        if (a>b){tmp=a;a=b;b=tmp;sw++;}
        if (b>c){tmp=b;b=c;c=tmp;sw++;}
        if (a>b){tmp=a;a=b;b=tmp;sw++;}
        g_chi_idx[p][q][r] = (uint8_t)tbl[a][b][c];
        g_chi_sgn[p][q][r] = (sw % 2 == 0) ? 1 : -1;
    }
}

/* chi(p,q,r) for record chirotope vector cv */
static inline int8_t chi(const int8_t *cv, int p, int q, int r) {
    return cv[g_chi_idx[p][q][r]] * g_chi_sgn[p][q][r];
}

/* -------------------------------------------------------------------------
 * Chirotope-based convex-hull tests (non-degenerate data: chi ∈ {+1,-1})
 * ------------------------------------------------------------------------- */

/* p strictly inside conv(S[0..k-1]):
 *   NOT ∃ i : ∀ j≠i : chi(p,S[i],S[j]) = +1  */
static int inside_hull(const int8_t *cv, int p, const int *S, int k) {
    for (int i = 0; i < k; i++) {
        int all_pos = 1;
        for (int j = 0; j < k; j++) {
            if (j == i) continue;
            if (chi(cv, p, S[i], S[j]) != 1) { all_pos = 0; break; }
        }
        if (all_pos) return 0;   /* S[i] is a witness: p is outside */
    }
    return 1;
}

/* Segment [u,v] meets conv(S[0..k-1]):
 *   u ∈ conv(S)  OR  v ∈ conv(S)  OR  [u,v] crosses a hull edge of S */
static int seg_meets_hull(const int8_t *cv, int u, int v, const int *S, int k) {
    if (k == 0) return 0;
    if (inside_hull(cv, u, S, k)) return 1;
    if (inside_hull(cv, v, S, k)) return 1;
    if (k <= 2) return 0;

    /* Find convex hull edges of S and test crossing with [u,v].
     * (S[i],S[j]) is a hull edge iff all other S[l] lie on the same side:
     *   chi(S[i],S[j],S[l]) = constant for l≠i,j.
     * [u,v] crosses [S[i],S[j]] iff:
     *   chi(S[i],S[j],u) ≠ chi(S[i],S[j],v)   (u,v on opposite sides of edge line)
     *   chi(u,v,S[i])   ≠ chi(u,v,S[j])        (S[i],S[j] on opposite sides of [u,v]) */
    for (int i = 0; i < k; i++) {
        for (int j = i+1; j < k; j++) {
            int8_t ref = 0, is_edge = 1;
            for (int l = 0; l < k; l++) {
                if (l == i || l == j) continue;
                int8_t s = chi(cv, S[i], S[j], S[l]);
                if (ref == 0) ref = s;
                else if (s != ref) { is_edge = 0; break; }
            }
            if (!is_edge) continue;

            int8_t c_u = chi(cv, S[i], S[j], u);
            int8_t c_v = chi(cv, S[i], S[j], v);
            if (c_u == c_v) continue;

            int8_t d_i = chi(cv, u, v, S[i]);
            int8_t d_j = chi(cv, u, v, S[j]);
            if (d_i != d_j) return 1;
        }
    }
    return 0;
}

/* -------------------------------------------------------------------------
 * is_face test for Gale diagram with positive pair (up, vp)
 * S_mask: bitmask of facets in the candidate face S
 * ------------------------------------------------------------------------- */

static int is_face(const int8_t *cv, int n, int up, int vp, uint16_t S_mask) {
    uint16_t all  = (uint16_t)((1u << n) - 1u);
    uint16_t T    = all & ~S_mask;
    if (!T) return 1;                       /* S = [n]: whole polytope is a face */

    int has_u = (T >> up) & 1;
    int has_v = (T >> vp) & 1;
    if (!has_u && !has_v) return 0;         /* pos_T empty */

    /* Extract neg_T = T minus the positive points */
    int neg[MAX_N], k = 0;
    uint16_t neg_mask = T & ~((1u << up) | (1u << vp));
    for (int i = 0; i < n; i++)
        if ((neg_mask >> i) & 1) neg[k++] = i;
    if (!k) return 0;                       /* pos_T non-empty, neg_T empty */

    if (has_u && has_v)
        return seg_meets_hull(cv, up, vp, neg, k);
    return inside_hull(cv, has_u ? up : vp, neg, k);
}

/* -------------------------------------------------------------------------
 * Precomputed list of all subsets of size 2..max_size (for n=10, max=5)
 * Ordered: all size-2 first, then size-3, etc. (needed for upward-closure).
 * ------------------------------------------------------------------------- */

#define MAX_SUBSETS 627
static uint16_t g_smask[MAX_SUBSETS];
static int8_t   g_ssize[MAX_SUBSETS];
static int      g_nsub;

static void init_subsets(int n, int max_size) {
    g_nsub = 0;
    for (int sz = 2; sz <= max_size; sz++) {
        uint32_t mask = (1u << sz) - 1u;
        while (mask < (1u << n)) {
            g_smask[g_nsub]   = (uint16_t)mask;
            g_ssize[g_nsub++] = (int8_t)sz;
            /* Gosper's hack: next same-popcount integer */
            uint32_t c = mask & (uint32_t)(-(int32_t)mask);
            uint32_t r = mask + c;
            mask = (((r ^ mask) >> 2) / c) | r;
        }
    }
}

/* -------------------------------------------------------------------------
 * Missing-face computation for a single (up, vp) pair
 * Uses an is_nf[1<<n] boolean array for O(1) sub-non-face lookup.
 * Stores missing-face bitmasks in mf_out[0..(*num_mf)-1].
 * Sets *p_count = number of size-2 missing faces.
 * ------------------------------------------------------------------------- */

static void compute_mf(const int8_t *cv, int n, int up, int vp,
                       uint16_t *mf_out, int *num_mf, int *p_count) {
    uint8_t is_nf[1 << MAX_N];
    memset(is_nf, 0, (size_t)(1 << n));
    *num_mf  = 0;
    *p_count = 0;

    for (int k = 0; k < g_nsub; k++) {
        uint16_t S  = g_smask[k];
        int      sz = g_ssize[k];

        /* Skip if any proper subset of size≥2 is already a non-face */
        if (sz > 2) {
            int found = 0;
            uint16_t sub = (S - 1u) & S;   /* iterate all proper subsets of S */
            while (sub) {
                if (__builtin_popcount(sub) >= 2 && is_nf[sub]) { found = 1; break; }
                sub = (sub - 1u) & S;
            }
            if (found) continue;
        }

        if (!is_face(cv, n, up, vp, S)) {
            is_nf[S] = 1;
            if (*num_mf < MAX_MF) mf_out[(*num_mf)++] = S;
            if (sz == 2) (*p_count)++;
        }
    }
}

/* -------------------------------------------------------------------------
 * Hash table (open-addressing, FNV-1a, key = sorted uint16 face bitmasks)
 * One entry per unique raw missing-face pattern.
 * ------------------------------------------------------------------------- */

#define HT_BITS 22
#define HT_CAP  (1u << HT_BITS)
#define HT_MASK (HT_CAP - 1u)

typedef struct {
    uint16_t key[MAX_MF];    /* 0xFFFF in key[0] means empty slot */
    uint32_t record_id;
    uint8_t  u, v, num_mf;
} HtEntry;

static HtEntry *g_ht;
static uint64_t g_ht_count;

static void ht_init(void) {
    g_ht = calloc(HT_CAP, sizeof(HtEntry));
    if (!g_ht) { perror("ht calloc"); exit(1); }
    for (uint32_t i = 0; i < HT_CAP; i++) g_ht[i].key[0] = 0xFFFF;
    g_ht_count = 0;
}

static uint64_t fnv1a(const uint16_t *d, int n) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (int i = 0; i < n; i++) { h ^= (uint64_t)d[i]; h *= 0x100000001b3ULL; }
    return h;
}

static void sort16(uint16_t *a, int n) {
    for (int i = 1; i < n; i++) {
        uint16_t x = a[i]; int j = i - 1;
        for (; j >= 0 && a[j] > x; j--) a[j+1] = a[j];
        a[j+1] = x;
    }
}

/* Returns 1 if inserted (new pattern), 0 if duplicate. */
static int ht_insert(uint32_t rec, uint8_t u, uint8_t v,
                     uint16_t *mf, uint8_t nm) {
    uint16_t key[MAX_MF];
    memset(key, 0xFF, sizeof(key));
    memcpy(key, mf, nm * sizeof(uint16_t));

    uint32_t pos = (uint32_t)(fnv1a(key, MAX_MF) & HT_MASK);
    for (;;) {
        HtEntry *e = &g_ht[pos];
        if (e->key[0] == 0xFFFF) {
            memcpy(e->key, key, sizeof(key));
            e->record_id = rec; e->u = u; e->v = v; e->num_mf = nm;
            g_ht_count++;
            return 1;
        }
        if (memcmp(e->key, key, sizeof(key)) == 0) return 0;
        pos = (pos + 1) & HT_MASK;
    }
}

/* -------------------------------------------------------------------------
 * Main
 * ------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: stage2_filter <input.chi> <d> <output.gd2>\n");
        return 1;
    }
    const char *inpath  = argv[1];
    int         d       = atoi(argv[2]);
    const char *outpath = argv[3];
    int         n       = d + 4;

    if (n < 5 || n > MAX_N) {
        fprintf(stderr, "error: n=d+4=%d outside [5,%d]\n", n, MAX_N);
        return 1;
    }

    init_chi_tables(n);
    init_subsets(n, 5);
    ht_init();

    fprintf(stderr, "Input:   %s  (d=%d, n=%d, subsets=%d)\n",
            inpath, d, n, g_nsub);

    /* mmap input */
    int fd = open(inpath, O_RDONLY);
    if (fd < 0) { perror("open"); return 1; }
    struct stat st; fstat(fd, &st);
    size_t fsz = (size_t)st.st_size;
    uint8_t *base = mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
    if (base == MAP_FAILED) { perror("mmap"); return 1; }
    close(fd);

    uint32_t magic_in, nn_in; uint64_t nrec;
    memcpy(&magic_in, base,   4);
    memcpy(&nn_in,    base+4, 4);
    memcpy(&nrec,     base+8, 8);
    if (magic_in != MAGIC_IN) { fprintf(stderr, "bad input magic\n"); return 1; }
    if ((int)nn_in != n) { fprintf(stderr, "n mismatch: file=%u arg=%d\n", nn_in, n); return 1; }

    int    ntrip = n*(n-1)*(n-2)/6;
    size_t row   = (size_t)(ntrip + n*2*2);
    uint8_t *recs = base + 16;
    fprintf(stderr, "Records: %llu  row=%zu bytes\n", (unsigned long long)nrec, row);

    /* Open output */
    FILE *fout = fopen(outpath, "wb");
    if (!fout) { perror("fopen output"); return 1; }
    { uint32_t mo = MAGIC_OUT, no = (uint32_t)n;
      fwrite(&mo, 4, 1, fout); fwrite(&no, 4, 1, fout); }
    long cnt_off = ftell(fout);
    { uint64_t cp = 0; fwrite(&cp, 8, 1, fout); }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    fprintf(stderr, "Threads: %d\n", nthreads);

    struct timespec t0; clock_gettime(CLOCK_MONOTONIC, &t0);
    uint64_t n_valid_pairs = 0;

    #pragma omp parallel for schedule(dynamic, 2000) \
        reduction(+:n_valid_pairs) default(none) \
        shared(recs, row, nrec, n, g_ht, g_ht_count, t0, stderr)
    for (int64_t r = 0; r < (int64_t)nrec; r++) {
        const int8_t *cv = (const int8_t *)(recs + (uint64_t)r * row);

        for (int u = 0; u < n; u++) for (int v = u+1; v < n; v++) {
            int S[MAX_N], k = 0;
            for (int i = 0; i < n; i++) if (i != u && i != v) S[k++] = i;

            if (!inside_hull(cv, u, S, k)) continue;
            if (!inside_hull(cv, v, S, k)) continue;

            int ok = 1;
            for (int i = 0; i < n && ok; i++)
                if (!is_face(cv, n, u, v, (uint16_t)(1u << i))) ok = 0;
            if (!ok) continue;

            uint16_t mf[MAX_MF]; int nm, pc;
            compute_mf(cv, n, u, v, mf, &nm, &pc);
            if (pc < 2) continue;

            n_valid_pairs++;
            sort16(mf, nm);
            #pragma omp critical
            ht_insert((uint32_t)r, (uint8_t)u, (uint8_t)v, mf, (uint8_t)nm);
        }

        if (r % 1000000 == 999999) {
            #pragma omp critical(progress)
            if (r % 1000000 == 999999) {
                struct timespec tn; clock_gettime(CLOCK_MONOTONIC, &tn);
                double el = (tn.tv_sec-t0.tv_sec)+(tn.tv_nsec-t0.tv_nsec)*1e-9;
                fprintf(stderr,
                        "  %lld/%llu (%.0f%%)  %.2f Mrec/s  ~%.0fs  "
                        "unique=%llu\n",
                        (long long)(r+1), (unsigned long long)nrec,
                        100.0*(r+1)/nrec, (r+1)/el/1e6,
                        (nrec-(uint64_t)(r+1))*el/(r+1),
                        (unsigned long long)g_ht_count);
            }
        }
    }

    /* Write all unique entries */
    for (uint32_t i = 0; i < HT_CAP; i++) {
        HtEntry *e = &g_ht[i];
        if (e->key[0] == 0xFFFF) continue;
        fwrite(&e->record_id, 4, 1, fout);
        fwrite(&e->u,         1, 1, fout);
        fwrite(&e->v,         1, 1, fout);
        fwrite(&e->num_mf,    1, 1, fout);
        fwrite(e->key, 2, e->num_mf, fout);
    }

    /* Fix up count */
    fseek(fout, cnt_off, SEEK_SET);
    fwrite(&g_ht_count, 8, 1, fout);
    fclose(fout);
    munmap(base, fsz);
    free(g_ht);

    struct timespec te; clock_gettime(CLOCK_MONOTONIC, &te);
    double el = (te.tv_sec - t0.tv_sec) + (te.tv_nsec - t0.tv_nsec)*1e-9;
    fprintf(stderr, "\nDone.\n");
    fprintf(stderr, "  Records processed : %llu\n", (unsigned long long)nrec);
    fprintf(stderr, "  Valid pairs (p≥2) : %llu\n", (unsigned long long)n_valid_pairs);
    fprintf(stderr, "  Unique raw types  : %llu\n", (unsigned long long)g_ht_count);
    fprintf(stderr, "  Time              : %.1fs  (%.2f Mrec/s)\n",
            el, nrec/el/1e6);
    return 0;
}
