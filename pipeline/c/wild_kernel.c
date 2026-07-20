/* Fast kernel for the wildcard rank-deficiency objective (pipeline/stage4_gram.py
 * _wild_feasible / objective_grad). Same math as the numpy implementation —
 * a cyclic Jacobi eigendecomposition of a small (n<=16) symmetric matrix plus
 * the analytic gradient of sum(lambda_i^2) over the num_zero smallest-|lambda|
 * eigenvalues — just without numpy/scipy Python-dispatch overhead per call.
 *
 * Profiling (2026-07-20) showed this exact call, invoked ~460k times per 60s
 * window during L-BFGS-B optimization of the wildcard range-analysis, spends
 * the overwhelming majority of wall time in numpy/scipy generic array-dispatch
 * machinery around a 9x9 or 10x10 matrix — work that is microseconds of real
 * FLOPs. This kernel does the identical computation with no array-library
 * overhead: fixed-size stack buffers, no allocation, no dispatch.
 */
#include <math.h>
#include <string.h>

#define MAXN 16

/* Classic cyclic Jacobi eigenvalue algorithm for a symmetric n x n matrix
 * (row-major, A[i*n+j]).  Outputs eigenvalues w[n] and eigenvectors as the
 * COLUMNS of V (row-major, V[i*n+j] = component i of eigenvector j) — i.e.
 * identical convention to numpy.linalg.eigh's (w, v) with v[:, j] the
 * j-th eigenvector. */
static void jacobi_eigh(const double *A, int n, double *w, double *V) {
    double a[MAXN][MAXN];
    double v[MAXN][MAXN];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = A[i * n + j];
            v[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    const int max_sweeps = 100;
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        double off = 0.0;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++)
                off += a[p][q] * a[p][q];
        if (off < 1e-300) break;
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                if (fabs(a[p][q]) < 1e-300) continue;
                double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                double t = (theta >= 0 ? 1.0 : -1.0) /
                          (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0);
                double s = t * c;
                double apq = a[p][q];
                a[p][p] -= t * apq;
                a[q][q] += t * apq;
                a[p][q] = 0.0; a[q][p] = 0.0;
                for (int i = 0; i < n; i++) {
                    if (i != p && i != q) {
                        double aip = a[i][p], aiq = a[i][q];
                        a[i][p] = c * aip - s * aiq; a[p][i] = a[i][p];
                        a[i][q] = s * aip + c * aiq; a[q][i] = a[i][q];
                    }
                }
                for (int i = 0; i < n; i++) {
                    double vip = v[i][p], viq = v[i][q];
                    v[i][p] = c * vip - s * viq;
                    v[i][q] = s * vip + c * viq;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) w[i] = a[i][i];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V[i * n + j] = v[i][j];
}

/* f = sum of w_i^2 over the num_zero smallest-|w| eigenvalues of G, where G is
 * G_base with entries (ui[e],uj[e]) and (uj[e],ui[e]) overwritten to -u[e] for
 * e in [0,k). grad_out[e] = -4 * sum_{i in that set} w_i * V[ui[e],i] * V[uj[e],i]
 * (same formula as pipeline/stage4_gram.py:objective_grad). Returns 0 on
 * success, -1 if n exceeds the fixed buffer size. */
int wild_objective_grad(const double *G_base, int n,
                        const int *ui, const int *uj, int k,
                        const double *u, int num_zero,
                        double *f_out, double *grad_out) {
    if (n > MAXN) return -1;
    double A[MAXN * MAXN];
    memcpy(A, G_base, sizeof(double) * (size_t)n * (size_t)n);
    for (int e = 0; e < k; e++) {
        A[ui[e] * n + uj[e]] = -u[e];
        A[uj[e] * n + ui[e]] = -u[e];
    }
    double w[MAXN], V[MAXN * MAXN];
    jacobi_eigh(A, n, w, V);

    int idx[MAXN];
    for (int i = 0; i < n; i++) idx[i] = i;
    for (int i = 0; i < n; i++) {
        int best = i;
        for (int j = i + 1; j < n; j++)
            if (fabs(w[idx[j]]) < fabs(w[idx[best]])) best = j;
        int tmp = idx[i]; idx[i] = idx[best]; idx[best] = tmp;
    }

    double f = 0.0;
    for (int m = 0; m < num_zero; m++) f += w[idx[m]] * w[idx[m]];
    *f_out = f;

    for (int e = 0; e < k; e++) {
        double g = 0.0;
        int i0 = ui[e], j0 = uj[e];
        for (int m = 0; m < num_zero; m++) {
            int c = idx[m];
            g += w[c] * V[i0 * n + c] * V[j0 * n + c];
        }
        grad_out[e] = -4.0 * g;
    }
    return 0;
}

/* objective-only variant (no gradient) for the cheap probe evaluations. */
double wild_objective(const double *G_base, int n,
                      const int *ui, const int *uj, int k,
                      const double *u, int num_zero) {
    if (n > MAXN) return -1.0;
    double A[MAXN * MAXN];
    memcpy(A, G_base, sizeof(double) * (size_t)n * (size_t)n);
    for (int e = 0; e < k; e++) {
        A[ui[e] * n + uj[e]] = -u[e];
        A[uj[e] * n + ui[e]] = -u[e];
    }
    double w[MAXN], V[MAXN * MAXN];
    jacobi_eigh(A, n, w, V);
    /* partial selection: only need the num_zero smallest |w| */
    for (int m = 0; m < num_zero; m++) {
        int best = m;
        for (int j = m + 1; j < n; j++)
            if (fabs(w[j]) < fabs(w[best])) best = j;
        double tmp = w[m]; w[m] = w[best]; w[best] = tmp;
    }
    double f = 0.0;
    for (int m = 0; m < num_zero; m++) f += w[m] * w[m];
    return f;
}
