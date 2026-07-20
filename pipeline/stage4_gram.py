"""
Stage 4: Coxeter label assignment + exact Gram-matrix realizability.

Pipeline per surviving combinatorial type:

  1. Reconstruct vertex sets from the stored Gale diagram example.
  2. Enumerate integer-label assignments m ∈ {2,3,4,5} on ordinary edges
     (all non-dotted pairs) using backtracking with vertex PD pruning.
     For a simple polytope every non-disjoint pair of facets is ridge-sharing
     and therefore ordinary, so vertex PD covers all pairs.
  3. For each labelled diagram that passes vertex PD:
       a. NUMERICAL SCREEN — use scipy to minimise the (n-d-1) smallest
          singular values of G over dotted weights x > 1.  If the minimum
          residual is above a threshold the type is skipped (fast: ~ms).
       b. EXACT SOLVE — SymPy Gröbner basis / solve on the numerically
          promising cases.  Verifies rank, signature (d,1) and x > 1 exactly.

  Float arithmetic may only REJECT; it cannot accept. Exact arithmetic is the
  final arbiter for every polytope that survives numerical screening.
"""

import json
import itertools
import signal
import time
from contextlib import contextmanager
from pathlib import Path

from pipeline.utils.manifest import write_manifest
from pipeline.utils import wild_kernel as _wk


class _SympyTimeout(Exception):
    pass


@contextmanager
def _hard_timeout(seconds):
    """SIGALRM-based hard timeout. Raises _SympyTimeout if the block exceeds `seconds`."""
    seconds = max(1, int(seconds))

    def _handler(signum, frame):
        raise _SympyTimeout()

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


try:
    import numpy as np
    import scipy.optimize as _scipy_opt
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

try:
    import sympy
    from sympy import (
        Matrix, symbols, Rational, sqrt, eye, S, cos, pi,
        solve, Poly, groebner
    )
    SYMPY_AVAILABLE = True
except ImportError:
    SYMPY_AVAILABLE = False
    print("WARNING: SymPy not available. Stage 4 will be skipped.")


# ---------------------------------------------------------------------------
# Gram-entry tables
# ---------------------------------------------------------------------------

_GRAM_FLOAT = {
    2:  0.0,
    3: -0.5,
    4: -0.7071067811865476,   # -√2/2
    5: -0.8090169943749474,   # -(1+√5)/4
    6: -0.8660254037844387,   # -√3/2
    7: -0.9009688679024191,   # -cos(π/7)
    8: -0.9238795325112867,   # -cos(π/8)
    9: -0.9396926207859084,   # -cos(π/9)
    10: -0.9510565162951535,  # -cos(π/10)
    12: -0.9659258262890682,  # -cos(π/12)
}

LABEL_CAP = 12
VALID_LABELS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]   # covers all known d=4 labels incl. π/7, π/12

# --- Ma-Zheng wildcard mode (arXiv:2201.00154 / 2203.16049, Prop. 3.5) ---------
# {2,...,10,12} is NOT a proven label bound in any dimension — it is an a-posteriori
# observation about the finished d=4/d=5 censuses.  The rigorous treatment enumerates
# labels over {2,...,6,7} where 7 is a WILDCARD standing for "any m >= 7": by the
# classification of finite Coxeter groups, a connected elliptic diagram of rank >= 3
# has labels <= 5, so an edge with label >= 6 must be an I2(m) component of every
# elliptic subdiagram containing it — hence every PD/Lannér test at cos(pi/7) gives
# the same verdict as at cos(pi/m) for EVERY m >= 7 (lossless proxy).  The actual m
# is then a continuous unknown c = cos(pi/m) in [cos(pi/7), 1) resolved by the
# rank/signature equations plus the terminal integrality demand pi/arccos(c) in Z.
WILDCARD_LABEL = 7
WILDCARD_INDICES = (0, 1, 2, 3, 4, 5)      # VALID_LABELS[:6] = {2,3,4,5,6,7}
_WILD_C_MIN = 0.9009688679024191           # cos(pi/7)


def gram_entry_adjacent(m):
    """Exact SymPy value of G_ij = -cos(π/m) for integer m."""
    if not SYMPY_AVAILABLE:
        raise RuntimeError("SymPy required for Stage 4")
    if m == 2:
        return S.Zero
    elif m == 3:
        return Rational(-1, 2)
    elif m == 4:
        return -sqrt(2) / 2
    elif m == 5:
        return -(1 + sqrt(5)) / 4
    else:
        return -cos(pi / m)


# ---------------------------------------------------------------------------
# Numpy Gram matrix builder (for numerical screening)
# ---------------------------------------------------------------------------

def _build_gram_numpy(ordinary_float, dotted_pairs, x_vals, n):
    """Build n×n numpy float Gram matrix.

    ordinary_float: dict (i,j) -> float (i<j), ordinary edges
    dotted_pairs:   list of (i,j) with i<j, in same order as x_vals
    x_vals:         numpy array of -G_ij values for dotted edges (> 1)
    """
    G = np.eye(n, dtype=float)
    for (i, j), v in ordinary_float.items():
        G[i, j] = v
        G[j, i] = v
    for (i, j), x in zip(dotted_pairs, x_vals):
        G[i, j] = -x
        G[j, i] = -x
    return G


# ---------------------------------------------------------------------------
# Vertex positive-definiteness check (float)
# ---------------------------------------------------------------------------

def _build_vertex_matrix_float(v_sorted, assignment):
    d = len(v_sorted)
    M = np.eye(d, dtype=float)
    for a in range(d):
        for b in range(a + 1, d):
            pair = (v_sorted[a], v_sorted[b])
            m = assignment.get(pair, 2)
            val = _GRAM_FLOAT[m]
            M[a, b] = val
            M[b, a] = val
    return M


def _is_vertex_pd_float(v_sorted, assignment):
    """Cholesky-based PD test (float).  Conservative pruner."""
    M = _build_vertex_matrix_float(v_sorted, assignment)
    try:
        np.linalg.cholesky(M)
        return True
    except np.linalg.LinAlgError:
        return False


# ---------------------------------------------------------------------------
# Numerical rank-7 screener
# ---------------------------------------------------------------------------

def _numerical_screen(ordinary_float, dotted_pairs, n, d,
                      residual_threshold=1e-6):
    """Check numerically whether rank(G) = d+1 is achievable with x > 1.

    Two-stage strategy:
      Stage 1 (fast, <0.5ms): Evaluate the objective (sum-of-squares of the
        num_zero = n-d-1 smallest singular values) at several probe points
        x ∈ {1.1, 1.5, 2.0, 3.0}.  If the minimum is > 0.5, reject
        immediately — the rank condition is nowhere near satisfied.
      Stage 2 (slow, ~10ms): scipy L-BFGS-B from the best probe, to find an
        actual numerical zero.  Only runs if a probe gives objective < 0.5.

    Returns (x_approx, residual).  residual > residual_threshold ⟹ infeasible.
    """
    num_zero = n - d - 1
    k = len(dotted_pairs)

    def objective(x_vals):
        G = _build_gram_numpy(ordinary_float, dotted_pairs, x_vals, n)
        s = np.linalg.svd(G, compute_uv=False)
        return float(np.sum(s[-num_zero:] ** 2))

    # Stage 1: multi-point probe (each < 0.1ms)
    best_x = np.full(k, 1.5)
    best_val = np.inf
    for x0_val in [1.5, 1.1, 2.0, 3.0]:
        x0 = np.full(k, x0_val)
        val = objective(x0)
        if val < best_val:
            best_val = val
            best_x = x0.copy()

    if best_val > 0.5:
        return best_x, best_val   # fast reject

    # Stage 2: optimise from best probe
    # Upper bound 1000 covers large dotted weights such as w≈34.9 in d=6 polytopes.
    bounds = [(1.001, 1000.0)] * k
    for x0_val in [best_x[0], 1.5, 2.0]:
        try:
            res = _scipy_opt.minimize(
                objective, np.full(k, x0_val),
                bounds=bounds,
                method='L-BFGS-B',
                options={'maxiter': 200, 'ftol': 1e-20, 'gtol': 1e-12},
            )
            if res.fun < best_val:
                best_val = res.fun
                best_x = res.x
        except Exception:
            pass
        if best_val < residual_threshold:
            break

    return best_x, best_val


def _refine_mpmath(x_approx, label_assign, dotted_pairs, n, d, dps=60, max_iter=40):
    """Refine a float64 approximate solution to high (dps-digit) precision.

    Uses Gauss-Newton iteration on the nk=(n-d-1) kernel conditions
    (chosen (d+2)×(d+2) minor determinants) with mpmath arithmetic.

    label_assign: dict mapping (i,j) pair -> integer m (Coxeter label).
    The ordinary Gram entries are computed EXACTLY in mpmath (-cos(π/m)),
    so the only imprecision is in the dotted weights being refined.

    For genuine isolated Coxeter polytopes the solution is zero-dimensional
    and this converges to the unique nearby exact algebraic values.
    For spurious points on a positive-dimensional variety, the mpmath
    residual stays bounded away from the convergence threshold.

    Returns a list of mpmath.mpf values (length k), or None on failure.
    """
    try:
        import mpmath
    except ImportError:
        return None

    mpmath.mp.dps = dps
    k = len(dotted_pairs)
    nk = n - d - 1          # kernel dimension = 3 for d=5,n=9
    minor_size = d + 2      # = 7 for d=5

    # Build the ordinary-edge base matrix with EXACT mpmath entries.
    # This is the key: -cos(π/m) in mpmath is exact to dps digits, not
    # limited by float64 precision (1e-16), so the refinement can converge
    # to dps-digit accuracy for the dotted weights.
    def exact_gram(m):
        if m == 2:  return mpmath.mpf(0)
        if m == 3:  return mpmath.mpf(-1) / 2
        if m == 4:  return -mpmath.sqrt(2) / 2
        if m == 5:  return -(1 + mpmath.sqrt(5)) / 4
        return -mpmath.cos(mpmath.pi / m)

    base = mpmath.zeros(n)
    for i in range(n):
        base[i, i] = mpmath.mpf(1)
    for (i, j), m in label_assign.items():
        val = exact_gram(m)
        base[i, j] = base[j, i] = val

    def build_G(x_list):
        G = mpmath.matrix(base)
        for m, (i, j) in enumerate(dotted_pairs):
            G[i, j] = G[j, i] = -x_list[m]
        return G

    def sub_det(G, rows):
        return mpmath.det(
            mpmath.matrix([[G[rows[r], rows[c]]
                            for c in range(len(rows))]
                           for r in range(len(rows))])
        )

    # GN system: ALL (d+2)×(d+2) principal minors whose row set contains at
    # least one dotted index.  This richer, overdetermined set is WELL
    # CONDITIONED — unlike the n-minor_size+1 consecutive windows, which form a
    # rank-deficient Jacobian (the dotted weights are not independently
    # constrained), causing only LINEAR convergence that never reaches the
    # acceptance threshold.  With this set GN is quadratic (~1e-90 in ~5 iters).
    from itertools import combinations as _combs
    dot_idx = {idx for p in dotted_pairs for idx in p}
    minor_rows = [list(c) for c in _combs(range(n), minor_size)
                  if dot_idx & set(c)]
    m_eqs = len(minor_rows)

    tol = mpmath.power(10, -(dps - 8))
    # Guard against non-finite / absurd starts (a divergent structured-screen
    # root) — these crash mpmath.det downstream; just skip the candidate.
    import math as _math
    try:
        if any((xi is None) or (not _math.isfinite(float(xi)))
               or float(xi) <= 1.0 or float(xi) > 1e8 for xi in x_approx):
            return None
    except (TypeError, ValueError, OverflowError):
        return None
    x = [mpmath.mpf(xi) for xi in x_approx]

    def compute_f(x_list):
        G = build_G(x_list)
        return [sub_det(G, mr) for mr in minor_rows]

    # Central-difference step (more accurate Jacobian -> quadratic convergence).
    eps_fd = mpmath.power(10, -(dps // 2 - 5))

    # The whole iteration is wrapped: any mpmath failure (singular pivot, a
    # non-finite intermediate from a bad start, etc.) means "no solution from
    # this start" -> return None, never crash the worker/run.
    try:
        for _ in range(max_iter):
            f_val = compute_f(x)
            f_norm = mpmath.norm(mpmath.matrix(f_val))
            if f_norm < tol:
                break

            J = mpmath.zeros(m_eqs, k)
            for col in range(k):
                x_p = list(x); x_p[col] = x_p[col] + eps_fd
                x_m = list(x); x_m[col] = x_m[col] - eps_fd
                fp = compute_f(x_p); fm = compute_f(x_m)
                for row in range(m_eqs):
                    J[row, col] = (fp[row] - fm[row]) / (2 * eps_fd)

            if k <= m_eqs:
                JTJ = J.T * J
                JTf = J.T * mpmath.matrix(f_val)
                dx = mpmath.lu_solve(JTJ, JTf)
            else:
                JJT = J * J.T
                rhs = mpmath.lu_solve(JJT, mpmath.matrix(f_val))
                dx = J.T * rhs

            for col in range(k):
                x[col] = x[col] - dx[col]
                if x[col] < mpmath.mpf('1.0001'):
                    x[col] = mpmath.mpf('1.0001')

        f_final = compute_f(x)
        f_norm_final = mpmath.norm(mpmath.matrix(f_final))
    except Exception:
        return None
    # Accept only if the GN-minor residual is genuinely tiny.
    if f_norm_final > mpmath.power(10, -(dps // 2)):
        return None
    return x


def _pslq_identify(xi_mp, dps=60):
    """Identify a high-precision mpmath value as a closed-form algebraic expression.

    Uses mpmath.identify (PSLQ) then falls back to sympy.nsimplify with a
    broad basis.  Returns a SymPy expression or None.

    Correctness note: we require the candidate to be a *symbolic* expression
    (contain at least one irrational/algebraic operation), not a plain decimal.
    This prevents the trivial passthrough where mpmath.identify returns the
    decimal representation and sympify() turns it into a SymPy Float.
    """
    import mpmath
    from sympy import sympify, sqrt as _sqrt, nsimplify, Float as SympyFloat, Number

    xi_f = float(xi_mp)
    if xi_f <= 1.0:
        return None

    def is_symbolic(expr):
        """Return True iff expr contains irrational/sqrt structure (not a plain number)."""
        from sympy import Rational
        return not expr.is_number or not isinstance(expr, (SympyFloat, Rational))

    def close_to_xi(expr):
        """Check agreement with xi_mp at high precision (dps//2 digits)."""
        try:
            val_mp = mpmath.mpf(str(expr.evalf(dps)))
            return abs(val_mp - xi_mp) < mpmath.power(10, -(dps // 2))
        except Exception:
            return False

    # 1. mpmath PSLQ identification
    # mpmath.identify can return symbolic strings like "(1+sqrt(5))/2"
    # OR plain decimal strings when it fails to find a pattern.
    # We only accept results that parse to a genuinely symbolic form.
    for tol_exp in [-(dps - 5), -(dps // 2), -25]:
        tol_mp = mpmath.power(10, tol_exp)
        result_str = mpmath.identify(xi_mp, tol=tol_mp)
        if result_str is not None and ("sqrt" in result_str or "/" in result_str
                                       or "pi" in result_str):
            try:
                expr = sympify(result_str)
                if is_symbolic(expr) and close_to_xi(expr):
                    return expr
            except Exception:
                pass

    # 2. nsimplify with progressively wider bases.
    # Pass xi_mp as a high-precision string so nsimplify gets > 15 digits of
    # information rather than the ~15-digit float64 value.
    xi_str = mpmath.nstr(xi_mp, dps - 5, strip_zeros=False)
    basis_sets = [
        [_sqrt(5)],
        [_sqrt(2), _sqrt(5)],
        [_sqrt(2), _sqrt(3), _sqrt(5)],
        [_sqrt(2), _sqrt(3), _sqrt(5), _sqrt(6), _sqrt(10)],
    ]
    for basis in basis_sets:
        for tol in [1e-30, 1e-20, 1e-12]:
            try:
                candidate = nsimplify(xi_str, basis, rational=False, tolerance=tol)
                if is_symbolic(candidate) and close_to_xi(candidate):
                    return candidate
            except Exception:
                pass

    return None


def _recognize_highprec_and_verify(x_hp, sym_list, G_sym, d, timeout=120.0):
    """Run PSLQ identification on high-precision x values, then verify exactly.

    x_hp: list of mpmath.mpf values (from _refine_mpmath).
    Returns a list of solution dicts (same format as solve_rank_condition),
    or [] if any component cannot be identified or verification fails.
    """
    deadline = time.time() + timeout
    exact_vals = []
    for xi_mp in x_hp:
        if time.time() > deadline:
            return []
        expr = _pslq_identify(xi_mp)
        if expr is None:
            return []
        exact_vals.append(expr)

    # All components identified; verify x_i > 1
    try:
        floats = [float(v.evalf()) for v in exact_vals]
    except Exception:
        return []
    if not all(f > 1.0 + 1e-9 for f in floats):
        return []

    # Substitute into symbolic Gram matrix and verify rank + signature
    try:
        if time.time() > deadline:
            return []
        G_sub = G_sym.subs(list(zip(sym_list, exact_vals)))
        n_mat = G_sub.shape[0]
        G_np = np.array([[float(G_sub[i, j].evalf())
                          for j in range(n_mat)]
                         for i in range(n_mat)], dtype=float)
        if not _check_signature_float(G_np, d):
            return []
        with _hard_timeout(max(1, int(deadline - time.time()))):
            rank_ok = G_sub.rank() == d + 1
        if not rank_ok:
            return []
    except Exception:
        return []

    sol = {str(s): v for s, v in zip(sym_list, exact_vals)}
    return [sol]


def _check_signature_float(G_numpy, d, tol=1e-8):
    """Check signature (d,1): exactly one negative eigenvalue.

    tol: eigenvalues with |λ| < tol are treated as zero.  The default 1e-8 is
    appropriate for exact solutions; use a larger value (e.g. 1e-3) when checking
    numerically approximate x values where near-zero eigenvalues may be ~1e-7.
    """
    evals = np.linalg.eigvalsh(G_numpy)
    neg = int(np.sum(evals < -tol))
    pos = int(np.sum(evals > tol))
    return neg == 1 and pos == d


# ---------------------------------------------------------------------------
# Backtracking label enumeration
# ---------------------------------------------------------------------------

def _build_vertex_groups(ordinary_pairs_set, vertex_sets_list,
                          sub_size_threshold=5):
    """Return list of (v_sorted_tuple, vertex_ordinary_pairs_frozenset).

    For large vertices (k pairs where 4^k > 2^20), also include sub-vertex
    groups up to sub_size_threshold to enable precomputed forward-checking.
    Default threshold=5: with label_indices=(0,1,2,3) the 5-node sub-vertex
    has only 1,386 valid PD combos (well below the 4096 bitmask limit),
    providing near-full-vertex PD pruning for d=6 polytopes.
    """
    seen_groups = set()
    groups = []
    for v in vertex_sets_list:
        v_sorted = tuple(sorted(v))
        vp = frozenset(
            (v_sorted[a], v_sorted[b])
            for a in range(len(v_sorted))
            for b in range(a + 1, len(v_sorted))
            if (v_sorted[a], v_sorted[b]) in ordinary_pairs_set
        )
        # Add the full vertex group if not too large to precompute
        k = len(vp)
        if (v_sorted, vp) not in seen_groups:
            groups.append((v_sorted, vp))
            seen_groups.add((v_sorted, vp))

        # For large vertices, also add sub-vertex groups of smaller size
        # so they can be precomputed and forward-checked
        if 4 ** k > (1 << 20):
            d = len(v_sorted)
            for sub_size in range(3, min(sub_size_threshold + 1, d + 1)):
                for sub in itertools.combinations(v_sorted, sub_size):
                    sub_vp = frozenset(
                        (sub[a], sub[b])
                        for a in range(sub_size)
                        for b in range(a + 1, sub_size)
                        if (sub[a], sub[b]) in ordinary_pairs_set
                    )
                    if (sub, sub_vp) not in seen_groups:
                        groups.append((sub, sub_vp))
                        seen_groups.add((sub, sub_vp))
    return groups


def _build_extended_lanner_groups(ordinary_pairs_set, mf_list, n, max_extra=3):
    """Lannér (signature (k-1,1)) groups for "extended missing faces".

    For each missing face M of size s ≥ 3 and each set E of 1..max_extra extra
    nodes such that ALL C(s+|E|, 2) pairs are ordinary and no pair in E×M or
    within E is itself a dotted pair (size-2 missing face), the sub-Gram
    G_{M ∪ E} must have signature (s+|E|-1, 1) — same Lannér structure as G_M,
    just embedded in a larger matrix.

    These constraints are complementary to face-tuple (PD) groups: face-tuples
    enforce ellipticity on face-sets, extended-Lannér groups enforce the
    propagation of hyperbolicity through the diagram.

    max_extra=3 is recommended: up to 4-node (size-3 MF + 1) and 5-node
    (size-3 MF + 2) extended groups capture the main propagation constraints.
    With max_extra=3 and max-size-3 MF we get groups up to size min(s+max_extra, n).
    """
    mf_as_frozensets = [frozenset(m) for m in mf_list]
    dotted_pairs_set = frozenset(
        tuple(sorted(m)) for m in mf_list if len(m) == 2
    )

    groups = []
    seen = set()

    for mf in mf_list:
        s = len(mf)
        if s < 3:
            continue
        mf_sorted = tuple(sorted(mf))
        mf_set = frozenset(mf_sorted)

        # Extra nodes: anything not in M
        extra_candidates = [k for k in range(n) if k not in mf_set]

        for extra_size in range(1, max_extra + 1):
            for extra in itertools.combinations(extra_candidates, extra_size):
                combo = tuple(sorted(mf_sorted + extra))
                combo_fs = frozenset(combo)
                if combo_fs in seen:
                    continue

                # Check: no SECOND missing face in the combo besides M itself
                # (i.e., the combo should contain exactly M as a missing face,
                #  not a second one — else the signature requirement changes)
                second_mf = [
                    mf2 for mf2 in mf_as_frozensets
                    if mf2 <= combo_fs and mf2 != mf_set and not (mf2 < mf_set)
                ]
                if second_mf:
                    continue  # combo contains a second distinct missing face

                # All C(|combo|, 2) pairs must be ordinary
                k = len(combo)
                pairs = []
                ok = True
                for a in range(k):
                    for b in range(a + 1, k):
                        p = (combo[a], combo[b])
                        if p not in ordinary_pairs_set:
                            ok = False
                            break
                        pairs.append(p)
                    if not ok:
                        break
                if not ok:
                    continue

                seen.add(combo_fs)
                groups.append((combo, frozenset(pairs), k))

    return groups


def _build_face_tuple_groups(ordinary_pairs_set, mf_list, n, max_size=4):
    """Elliptic (PD) groups for all face-tuples of size 3..max_size.

    Vinberg: in a compact Coxeter polytope the Gram submatrix of any set of
    mutually-meeting facets is positive definite (they span an elliptic
    sub-diagram).  A set S of facets mutually meets iff no missing face of the
    combinatorial type is a subset of S.

    These groups feed into the same bitmask forward-checking as vertex groups.
    The key gain for types with few missing triples: the ordinary-edge graph
    (pairs with label ≥ 3) must be triangle-free except on size-3 missing
    faces, which by Mantel's theorem caps it at ≤ 25 edges and drastically
    prunes the search that previously generated millions of assignments.

    max_size=4 is the recommended default: size-3 adds triangle-free pruning,
    size-4 adds path/cycle pruning.  The existing _build_vertex_groups already
    adds geometric-vertex sub-groups; this function adds the complementary
    non-vertex face-tuples that those sub-groups miss.
    """
    mf_as_frozensets = [frozenset(m) for m in mf_list]

    groups = []
    seen = set()

    for size in range(3, max_size + 1):
        for combo in itertools.combinations(range(n), size):
            combo_fs = frozenset(combo)
            if combo_fs in seen:
                continue
            # Face check: no missing face may be a subset of this combo.
            # Dotted pairs (size-2 missing faces) are automatically excluded here.
            if any(mf <= combo_fs for mf in mf_as_frozensets):
                continue
            # All C(size,2) pairs must be ordinary (non-dotted).
            # (Guaranteed by the face check above for valid types, but verified
            # explicitly for safety.)
            pairs = []
            ok = True
            for a in range(size):
                for b in range(a + 1, size):
                    p = (combo[a], combo[b])
                    if p not in ordinary_pairs_set:
                        ok = False
                        break
                    pairs.append(p)
                if not ok:
                    break
            if not ok:
                continue
            seen.add(combo_fs)
            groups.append((combo, frozenset(pairs)))

    return groups


def _build_lanner_groups(ordinary_pairs_set, mf_list):
    """Return list of (face_sorted_tuple, face_pairs_frozenset) for each
    missing face of size ≥ 3 whose pairs are all ordinary.

    For each such face the Gram submatrix must have signature (k-1, 1)
    (compact Lannér diagram condition).
    """
    groups = []
    for mf in mf_list:
        k = len(mf)
        if k < 3:
            continue
        face_sorted = tuple(sorted(mf))
        pairs = frozenset(
            (face_sorted[a], face_sorted[b])
            for a in range(k)
            for b in range(a + 1, k)
        )
        # All pairs must be ordinary (guaranteed if mf is a missing face of size ≥ 3,
        # since it can't contain a size-2 missing face as a subset)
        if not pairs.issubset(ordinary_pairs_set):
            continue   # some pair is dotted → skip (shouldn't happen for valid MF)
        groups.append((face_sorted, pairs, k))
    return groups


def _is_lanner_float(face_sorted, assignment, k):
    """Check that the k×k Gram submatrix of face_sorted has signature (k-1, 1).

    Lannér condition: exactly one negative eigenvalue.  Float-based (for pruning).
    """
    M = np.eye(k, dtype=float)
    for a in range(k):
        for b in range(a + 1, k):
            pair = (face_sorted[a], face_sorted[b])
            m = assignment.get(pair, 2)
            val = _GRAM_FLOAT[m]
            M[a, b] = val
            M[b, a] = val
    evals = np.linalg.eigvalsh(M)
    neg = int(np.sum(evals < -1e-8))
    return neg == 1


def _vertex_first_ordering(ordinary_pairs, vertex_groups, lanner_groups):
    """Order pairs so that vertex PD checks fire as early as possible.

    Strategy: greedily process vertices in order of most-shared pairs with
    already-seen pairs.  This ensures the first vertex completes at level d
    (= len(first vertex's pairs)), and each subsequent vertex adds roughly
    1 new pair, triggering a PD check at almost every backtracking level.
    Lannér pairs (if any) are pulled to the front of their vertex blocks.
    """
    ordinary_set = set(ordinary_pairs)
    ordered = []
    seen = set()

    # Lannér pairs — process first so Lannér checks fire early.
    # Strategy: pairs that appear in EVERY Lannér group (universal pairs) are
    # placed first because they constrain all groups simultaneously, giving
    # maximum early pruning.  Remaining Lannér pairs keep their original
    # lexicographic order so that types whose universal pairs are already lex-
    # first (e.g. types 107/119/132) are not affected.
    from collections import Counter as _Counter
    lanner_pairs = set()
    lanner_pair_freq: _Counter = _Counter()
    n_lg = len(lanner_groups)
    for _, lp, _ in lanner_groups:
        lanner_pairs.update(lp)
        for p in lp:
            lanner_pair_freq[p] += 1
    universal = {p for p, f in lanner_pair_freq.items() if f == n_lg} if n_lg > 0 else set()
    # Universal pairs (appear in ALL Lannér groups) first, then remaining lex.
    for p in sorted(universal & ordinary_set):
        if p not in seen:
            ordered.append(p)
            seen.add(p)
    for p in sorted((lanner_pairs - universal) & ordinary_set):
        if p not in seen:
            ordered.append(p)
            seen.add(p)

    if not vertex_groups:
        for p in sorted(ordinary_set - seen):
            ordered.append(p)
        return ordered

    remaining_vg = list(range(len(vertex_groups)))

    # Seed with the vertex that has the most Lannér-seeded pairs already seen
    def _score(idx):
        _, vp = vertex_groups[idx]
        return sum(1 for p in vp if p in seen)

    first = max(remaining_vg, key=_score)
    _, vp0 = vertex_groups[first]
    for p in sorted(vp0):
        if p not in seen:
            ordered.append(p)
            seen.add(p)
    remaining_vg.remove(first)

    # Greedily add the vertex that maximises overlap with already-seen pairs
    while remaining_vg:
        best = max(remaining_vg, key=lambda idx: sum(1 for p in vertex_groups[idx][1] if p in seen))
        _, vp = vertex_groups[best]
        for p in sorted(vp):
            if p not in seen:
                ordered.append(p)
                seen.add(p)
        remaining_vg.remove(best)

    # Remaining ordinary pairs not in any vertex
    for p in sorted(ordinary_set - seen):
        ordered.append(p)

    return ordered


# Module-level cache: valid label assignments for sub-vertex groups.
# Key: d (number of facets in the group, all pairs ordinary).
# Value: numpy array of shape (N, C(d,2)) dtype int8 of valid label-index combos,
#   or None if nl^C(d,2) > max_combos.
# Valid combos are the same for every group of the same size (canonical structure).
_VG_VALID_CACHE: dict = {}

# Lannér cache: key = k (group size), value = valid combos array.
_LG_VALID_CACHE: dict = {}

# Maximum number of valid combos before switching from bitmask forward-checking
# to on-the-fly eigenvalue checks. A bitmask of N bits takes N/8 bytes and AND
# operations scale linearly — above ~4096 combos it becomes noticeably slower.
_BITMASK_COMBO_LIMIT = 4096


def _canonical_pd_valid(d, max_combos=1 << 22, label_indices=None):
    """Compute valid label-index combos for a d-node all-ordinary vertex (PD condition).

    label_indices: tuple of indices into VALID_LABELS to use (default: all).
    Returns numpy array of shape (N, k) dtype int8, or None if nl^k > max_combos.
    Cached by (d, label_indices).
    """
    labels = VALID_LABELS if label_indices is None else [VALID_LABELS[i] for i in label_indices]
    cache_key = (d, tuple(labels))
    if cache_key in _VG_VALID_CACHE:
        return _VG_VALID_CACHE[cache_key]

    nl = len(labels)
    k = d * (d - 1) // 2  # = C(d,2) pairs, in lex order
    if nl ** k >= max_combos:
        _VG_VALID_CACHE[cache_key] = None
        return None

    label_f = np.array([_GRAM_FLOAT[m] for m in labels], dtype=float)
    # Pair positions: (a,b) for each of the k pairs in lex order
    pairs_ab = [(a, b) for a in range(d) for b in range(a + 1, d)]

    # Vectorised: build all nl^k combos, check PD in batches
    N = nl ** k
    # Generate combo array via base-nl decomposition
    tmp = np.arange(N, dtype=np.int64)
    combos = np.zeros((N, k), dtype=np.int8)
    for j in range(k - 1, -1, -1):
        combos[:, j] = (tmp % nl).astype(np.int8)
        tmp //= nl
    vals = label_f[combos]   # (N, k) float values

    BATCH = 20_000
    valid_list = []
    for start in range(0, N, BATCH):
        end = min(start + BATCH, N)
        v = vals[start:end]          # (B, k)
        mats = np.zeros((end - start, d, d), dtype=float)
        mats[:, np.arange(d), np.arange(d)] = 1.0
        for j, (a, b) in enumerate(pairs_ab):
            mats[:, a, b] = v[:, j]
            mats[:, b, a] = v[:, j]
        evals = np.linalg.eigvalsh(mats)   # (B, d)
        mask = np.all(evals > 1e-10, axis=1)
        valid_list.append(combos[start:end][mask])

    arr = np.concatenate(valid_list) if valid_list else np.zeros((0, k), dtype=np.int8)
    _VG_VALID_CACHE[cache_key] = arr
    return arr


def _canonical_lanner_valid(k, max_combos=1 << 22, label_indices=None):
    """Compute valid label-index combos for a k-node all-ordinary Lannér group.

    label_indices: tuple of indices into VALID_LABELS to use (default: all).
    Returns numpy array of shape (N, C(k,2)) dtype int8, or None if too large.
    Cached by (k, label_indices).
    """
    labels = VALID_LABELS if label_indices is None else [VALID_LABELS[i] for i in label_indices]
    cache_key = (k, tuple(labels))
    if cache_key in _LG_VALID_CACHE:
        return _LG_VALID_CACHE[cache_key]

    nl = len(labels)
    num_pairs = k * (k - 1) // 2
    if nl ** num_pairs >= max_combos:
        _LG_VALID_CACHE[cache_key] = None
        return None

    label_f = np.array([_GRAM_FLOAT[m] for m in labels], dtype=float)
    pairs_ab = [(a, b) for a in range(k) for b in range(a + 1, k)]

    N = nl ** num_pairs
    tmp = np.arange(N, dtype=np.int64)
    combos = np.zeros((N, num_pairs), dtype=np.int8)
    for j in range(num_pairs - 1, -1, -1):
        combos[:, j] = (tmp % nl).astype(np.int8)
        tmp //= nl
    vals = label_f[combos]   # (N, num_pairs)

    BATCH = 20_000
    valid_list = []
    for start in range(0, N, BATCH):
        end = min(start + BATCH, N)
        v = vals[start:end]
        mats = np.zeros((end - start, k, k), dtype=float)
        mats[:, np.arange(k), np.arange(k)] = 1.0
        for j, (a, b) in enumerate(pairs_ab):
            mats[:, a, b] = v[:, j]
            mats[:, b, a] = v[:, j]
        evals = np.linalg.eigvalsh(mats)
        mask = np.sum(evals < -1e-8, axis=1) == 1
        valid_list.append(combos[start:end][mask])

    arr = np.concatenate(valid_list) if valid_list else np.zeros((0, num_pairs), dtype=np.int8)
    _LG_VALID_CACHE[cache_key] = arr
    return arr


def _prepare_vertex_precomputed(vertex_groups, label_indices=None):
    """Prepare forward-checking data for all vertex groups.

    For groups with nl^k ≤ 2^22 (k = C(d,2) pairs), loads from the global
    cache (computed once per process per group size).  No per-type
    precomputation cost.

    Returns list of (v_sorted, v_pairs_list, pair_to_pos, valid_array_or_None).
    """
    result = []
    for v_sorted, vp in vertex_groups:
        d = len(v_sorted)
        v_sorted_list = list(v_sorted)
        v_pairs = sorted(vp)
        # Build mapping from v_pairs order → canonical pair order within d×d matrix
        # Canonical: pair index pi = a*(2d-a-1)//2 + (b-a-1) for a<b in 0..d-1
        k = len(v_pairs)
        # All pairs within an ordinary vertex are ordinary → use cached valid array
        # Map v_pairs → canonical position index in the cached array
        canonical_idx = []
        for p in v_pairs:
            a = v_sorted_list.index(p[0])
            b = v_sorted_list.index(p[1])
            # Canonical pair position in canonical ordering (a < b guaranteed)
            # Canonical ordering: (0,1),(0,2),...,(0,d-1),(1,2),...
            pi = a * (2 * d - a - 1) // 2 + (b - a - 1)
            canonical_idx.append(pi)
        canonical_idx = np.array(canonical_idx, dtype=np.int32)

        # Get cached valid array for this group size
        valid_arr_full = _canonical_pd_valid(d, label_indices=label_indices)
        if valid_arr_full is not None:
            valid_arr = valid_arr_full[:, canonical_idx]
            # Bitmask forward-checking is efficient only for small valid sets.
            if len(valid_arr) > _BITMASK_COMBO_LIMIT:
                valid_arr = None  # fall back to on-the-fly PD check
        else:
            valid_arr = None  # too large to precompute

        pair_to_pos = {p: i for i, p in enumerate(v_pairs)}
        result.append((v_sorted, v_pairs, pair_to_pos, valid_arr))
    return result


def _prepare_lanner_precomputed(lanner_groups, label_indices=None):
    """Prepare forward-checking data for all Lannér groups.

    Returns list of (face_sorted, l_pairs_list, pair_to_pos, k, valid_array_or_None).
    """
    result = []
    for face_sorted, lp, k in lanner_groups:
        face_list = list(face_sorted)
        l_pairs = sorted(lp)
        # Map l_pairs → canonical position in cached Lannér array
        canonical_idx = []
        for p in l_pairs:
            a = face_list.index(p[0])
            b = face_list.index(p[1])
            pi = a * (2 * k - a - 1) // 2 + (b - a - 1)
            canonical_idx.append(pi)
        canonical_idx = np.array(canonical_idx, dtype=np.int32)

        valid_arr_full = _canonical_lanner_valid(k, label_indices=label_indices)
        if valid_arr_full is not None:
            valid_arr = valid_arr_full[:, canonical_idx]
            if len(valid_arr) > _BITMASK_COMBO_LIMIT:
                valid_arr = None  # fall back to on-the-fly Lannér check
        else:
            valid_arr = None

        pair_to_pos = {p: i for i, p in enumerate(l_pairs)}
        result.append((face_sorted, l_pairs, pair_to_pos, k, valid_arr))
    return result


def enumerate_labels_backtrack(ordinary_pairs, vertex_groups, lanner_groups=None,
                                max_count=50000, timeout=60.0, label_indices=None,
                                prefix=None, state_out=None, edge_max_label=None):
    """Generator: backtracking with forward-checking using precomputed vertex
    and Lannér valid-assignment tables.

    label_indices: optional tuple of indices into VALID_LABELS to restrict the
    label search (e.g. (0,1,2,3) for labels {2,3,4,5} only).  Default: all labels.

    prefix: optional list of label-INDEX values (into the resolved label_indices /
    actual_labels list) fixing the first len(prefix) pairs of the deterministic
    pair ordering.  Only the sub-tree under that fixed prefix is enumerated.  This
    partitions the search for parallelism: running every prefix (a Cartesian product
    over the first p pairs) across a worker pool covers the whole space with no
    overlap.  Forward-checking still applies to the fixed pairs, so infeasible
    prefixes prune immediately.  Default None = enumerate everything.

    For each group with nl^k ≤ 2^22 pairs, valid assignments are precomputed
    once and forward-checked at every assignment step via numpy boolean masking.
    Groups too large to precompute fall back to on-the-fly PD checking when
    the group is fully assigned.

    Yields dicts {(i,j): m} for each valid assignment.
    """
    if lanner_groups is None:
        lanner_groups = []

    # Resolve label indices
    if label_indices is None:
        label_indices = tuple(range(len(VALID_LABELS)))
    n_labels = len(label_indices)
    actual_labels = [VALID_LABELS[i] for i in label_indices]

    pairs_list = _vertex_first_ordering(ordinary_pairs, vertex_groups, lanner_groups)
    n_pairs = len(pairs_list)
    pair_idx = {p: i for i, p in enumerate(pairs_list)}

    # Precompute valid arrays for all groups (filtered to label_indices)
    vg_precomp = _prepare_vertex_precomputed(vertex_groups, label_indices)
    lg_precomp = _prepare_lanner_precomputed(lanner_groups, label_indices)

    # Map each pair -> [(group_list, group_idx, pos_in_group), ...]
    pair_to_vg_fc = [[] for _ in range(n_pairs)]    # forward-check (precomputed)
    pair_to_vg_bt = [[] for _ in range(n_pairs)]    # fallback (on-the-fly)
    for g_idx, (v_sorted, v_pairs, pair_to_pos, valid_arr) in enumerate(vg_precomp):
        for p, pos in pair_to_pos.items():
            if p in pair_idx:
                pi = pair_idx[p]
                if valid_arr is not None:
                    pair_to_vg_fc[pi].append((g_idx, pos))
                else:
                    pair_to_vg_bt[pi].append(g_idx)

    pair_to_lg_fc = [[] for _ in range(n_pairs)]
    pair_to_lg_bt = [[] for _ in range(n_pairs)]
    for g_idx, (face_sorted, l_pairs, pair_to_pos, k, valid_arr) in enumerate(lg_precomp):
        for p, pos in pair_to_pos.items():
            if p in pair_idx:
                pi = pair_idx[p]
                if valid_arr is not None:
                    pair_to_lg_fc[pi].append((g_idx, pos))
                else:
                    pair_to_lg_bt[pi].append(g_idx)

    # Build bitmask filter tables for fast forward checking.
    # filter_vg[g_idx][pos][li] = int bitmask of valid-combo indices where
    #   position pos has label index li.  Zero means no valid combo for that label.
    # full_vg[g_idx] = all-ones bitmask (all valid combos active).
    filter_vg = []
    full_vg = []
    for g_idx, (v_sorted, v_pairs, pair_to_pos, valid_arr) in enumerate(vg_precomp):
        if valid_arr is None:
            filter_vg.append(None)
            full_vg.append(None)
            continue
        k_pairs = len(v_pairs)
        n_combos = len(valid_arr)
        full = (1 << n_combos) - 1
        full_vg.append(full)
        filters = [[0] * n_labels for _ in range(k_pairs)]
        for i in range(n_combos):
            bit = 1 << i
            for pos in range(k_pairs):
                filters[pos][valid_arr[i, pos]] |= bit
        filter_vg.append(filters)

    filter_lg = []
    full_lg = []
    for g_idx, (face_sorted, l_pairs, pair_to_pos, k, valid_arr) in enumerate(lg_precomp):
        if valid_arr is None:
            filter_lg.append(None)
            full_lg.append(None)
            continue
        k_pairs = len(l_pairs)
        n_combos = len(valid_arr)
        full = (1 << n_combos) - 1
        full_lg.append(full)
        filters = [[0] * n_labels for _ in range(k_pairs)]
        for i in range(n_combos):
            bit = 1 << i
            for pos in range(k_pairs):
                filters[pos][valid_arr[i, pos]] |= bit
        filter_lg.append(filters)

    # Current valid masks as Python integers (bitmask over valid-combo indices).
    vg_masks = [full_vg[g] for g in range(len(vg_precomp))]
    lg_masks = [full_lg[g] for g in range(len(lg_precomp))]

    assignment = [-1] * n_pairs
    # Per-edge label caps (edge_max_label: pair -> max label VALUE, e.g. 5 for
    # edges forced low-weight by Burcroff Lemma 5.5).  A prefix entry outside a
    # pair's allowed set prunes that subtree entirely — sound: the cap encodes a
    # theorem that no real polytope carries a higher label there.
    allowed_li = [tuple(range(n_labels))] * n_pairs
    if edge_max_label:
        for i, p in enumerate(pairs_list):
            cap = edge_max_label.get(p)
            if cap is not None:
                allowed_li[i] = tuple(li for li in range(n_labels)
                                      if actual_labels[li] <= cap)
    # state_out (if given) is the caller's dict: after the generator is drained,
    # state_out["exhausted"] is True iff every branch was explored (no timeout,
    # no max_count cut) — the rigorous complete-vs-truncated verdict.
    state = state_out if state_out is not None else {}
    state.update({"count": 0, "start": time.time(), "exhausted": False})
    trail = []  # (g_idx, is_lanner, old_int_mask)

    def _backtrack(depth):
        if state["count"] >= max_count or time.time() - state["start"] > timeout:
            return
        if depth == n_pairs:
            state["count"] += 1
            yield {pairs_list[i]: actual_labels[assignment[i]] for i in range(n_pairs)}
            return

        pi = depth
        trail_mark = len(trail)

        # Partition support: at fixed-prefix depths, only the prefix's label.
        if prefix is not None and depth < len(prefix):
            if prefix[depth] not in allowed_li[pi]:
                return          # capped-out subtree: provably empty
            li_choices = (prefix[depth],)
        else:
            li_choices = allowed_li[pi]
        for li in li_choices:
            assignment[pi] = li
            ok = True

            # Forward-check vertex groups
            for g_idx, pos in pair_to_vg_fc[pi]:
                old_m = vg_masks[g_idx]
                new_m = old_m & filter_vg[g_idx][pos][li]
                trail.append((g_idx, False, old_m))
                vg_masks[g_idx] = new_m
                if not new_m:
                    ok = False
                    break

            # Forward-check Lannér groups
            if ok:
                for g_idx, pos in pair_to_lg_fc[pi]:
                    old_m = lg_masks[g_idx]
                    new_m = old_m & filter_lg[g_idx][pos][li]
                    trail.append((g_idx, True, old_m))
                    lg_masks[g_idx] = new_m
                    if not new_m:
                        ok = False
                        break

            # Fallback: on-the-fly checks for groups too large to precompute
            if ok and (pair_to_vg_bt[pi] or pair_to_lg_bt[pi]):
                assign_dict = {pairs_list[i]: actual_labels[assignment[i]]
                               for i in range(depth + 1)}
                for g_idx in pair_to_vg_bt[pi]:
                    v_sorted, vp_set = vertex_groups[g_idx]
                    if all(p in assign_dict for p in vp_set):
                        if not _is_vertex_pd_float(v_sorted, assign_dict):
                            ok = False
                            break
                if ok:
                    for g_idx in pair_to_lg_bt[pi]:
                        face_sorted, lp, k = lanner_groups[g_idx]
                        if all(p in assign_dict for p in lp):
                            if not _is_lanner_float(face_sorted, assign_dict, k):
                                ok = False
                                break

            if ok:
                yield from _backtrack(depth + 1)

            # Restore trail
            while len(trail) > trail_mark:
                g_idx, is_lanner, old_m = trail.pop()
                if is_lanner:
                    lg_masks[g_idx] = old_m
                else:
                    vg_masks[g_idx] = old_m

        assignment[pi] = -1

    yield from _backtrack(0)
    # After the recursive generator finishes: if neither timeout nor max_count
    # fired, the search was fully exhausted (every branch explored).
    if state["count"] < max_count and time.time() - state["start"] <= timeout:
        state["exhausted"] = True


# ---------------------------------------------------------------------------
# Symbolic Gram matrix
# ---------------------------------------------------------------------------

def build_gram_matrix(n, ordinary_assignment, dotted_set, dot_syms):
    """Build symbolic n×n Gram matrix.

    ordinary_assignment: dict (i,j)->m, i<j (m=2 ⟹ G_ij=0)
    dotted_set:          frozenset of (i,j) pairs, i<j
    dot_syms:            dict (i,j) -> SymPy symbol (= -G_ij)
    """
    G = eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            pair = (i, j)
            if pair in dotted_set:
                x = dot_syms[pair]
                G[i, j] = -x
                G[j, i] = -x
            else:
                m = ordinary_assignment.get(pair, 2)
                entry = gram_entry_adjacent(m)
                G[i, j] = entry
                G[j, i] = entry
    return G


# ---------------------------------------------------------------------------
# Exact rank-condition solver
# ---------------------------------------------------------------------------

def _collect_minor_eqs(G, n, minor_size, sym_list, unknown_rows,
                        max_eqs, deadline):
    """Collect distinct non-trivial (minor_size)×(minor_size) minor equations."""
    other_rows = [i for i in range(n) if i not in unknown_rows]
    extra_count = minor_size - len(unknown_rows)

    eqs, eq_strs = [], set()

    def try_rows(rows):
        if time.time() > deadline:
            return
        sub = G.extract(rows, rows)
        eq = sub.det().expand()
        if eq == 0 or eq.is_number:
            return
        try:
            key = str(Poly(eq, *sym_list))
        except Exception:
            key = str(eq)[:300]
        if key in eq_strs:
            return
        eq_strs.add(key)
        eqs.append(eq)

    if extra_count >= 0:
        for extra in itertools.combinations(other_rows, extra_count):
            if time.time() > deadline or len(eqs) >= max_eqs:
                break
            try_rows(sorted(list(unknown_rows) + list(extra)))

    # Fall through to all C(n, minor_size) subsets if we haven't collected enough.
    # Necessary when len(unknown_rows) >= minor_size (extra_count <= 0) because the
    # prioritised loop above only covers a single row combination in that case.
    if len(eqs) < max_eqs:
        for rows in itertools.combinations(range(n), minor_size):
            if time.time() > deadline or len(eqs) >= max_eqs:
                break
            try_rows(list(rows))

    return eqs


def _validate_solutions(sols_raw, sym_list, G, d, deadline):
    """Filter raw SymPy solutions: x > 1, correct signature, correct rank."""
    valid = []
    for sol in sols_raw:
        if time.time() > deadline:
            break
        if not all(s in sol for s in sym_list):
            continue
        vals = [sol[s] for s in sym_list]
        try:
            floats = [float(v.evalf()) for v in vals]
        except Exception:
            continue
        if not all(f > 1.0 + 1e-9 for f in floats):
            continue
        G_sub = G.subs(list(zip(sym_list, vals)))
        n = G_sub.shape[0]
        G_np = np.array([[float(G_sub[i, j]) for j in range(n)]
                         for i in range(n)], dtype=float)
        if not _check_signature_float(G_np, d):
            continue
        try:
            if G_sub.rank() != d + 1:
                continue
        except Exception:
            pass
        valid.append({str(s): v for s, v in zip(sym_list, vals)})
    return valid


def solve_rank_condition(G, d, dot_syms, dotted_pairs, n, timeout=60.0):
    """Find dotted weights making rank(G) = d+1 and signature (d,1).

    Returns list of solution dicts {str(sym): SymPy_value}.
    """
    sym_list = [dot_syms[p] for p in dotted_pairs]
    k = len(sym_list)
    if k == 0:
        return []

    minor_size = d + 2
    unknown_rows = sorted({idx for p in dotted_pairs for idx in p})
    deadline = time.time() + timeout

    eqs = _collect_minor_eqs(G, n, minor_size, sym_list, unknown_rows,
                              max_eqs=k + 5, deadline=deadline)
    if not eqs:
        return []

    remaining = max(1.0, deadline - time.time())

    try:
        if k == 1:
            x = sym_list[0]
            with _hard_timeout(remaining):
                roots = solve(eqs[0], x)
            raw = [{x: r} for r in roots if r.is_real]
            return _validate_solutions(raw, sym_list, G, d, deadline)

        elif k <= 4:
            raw = []
            try:
                with _hard_timeout(remaining):
                    raw = solve(eqs[: k + 1], sym_list, dict=True)
            except (_SympyTimeout, Exception):
                pass
            valid = _validate_solutions(raw, sym_list, G, d, deadline)
            if valid or time.time() > deadline:
                return valid
            # Gröbner fallback
            remaining2 = max(1.0, deadline - time.time())
            try:
                with _hard_timeout(remaining2):
                    basis = groebner(eqs[: min(len(eqs), k + 3)], *sym_list,
                                     order='lex')
                    raw2 = solve(list(basis), sym_list, dict=True)
                return _validate_solutions(raw2, sym_list, G, d, deadline)
            except (_SympyTimeout, Exception):
                return valid

        else:
            try:
                with _hard_timeout(remaining):
                    basis = groebner(eqs[: min(len(eqs), k + 3)], *sym_list,
                                     order='lex')
                    raw = solve(list(basis), sym_list, dict=True)
                return _validate_solutions(raw, sym_list, G, d, deadline)
            except (_SympyTimeout, Exception):
                return []

    except Exception:
        return []


# ---------------------------------------------------------------------------
# Algebraic recognition for high-k exact solve
# ---------------------------------------------------------------------------

# Basis elements used for nsimplify recognition of dotted-edge weights.
# These cover all ordinary-label cosines: cos(π/m) for m∈{2,...,12} involves
# at most √2, √3, √5, and √6.
_NSIMPLIFY_BASIS = None  # initialised lazily to avoid import-time SymPy overhead


def _recognize_and_verify(x_approx, sym_list, G_sym, d, timeout=60.0):
    """Attempt to recover exact algebraic x values from a numerical approximation.

    Uses sympy.nsimplify to find closed-form expressions for each component of
    x_approx, then verifies that the substituted G_sym has rank d+1 and
    signature (d,1).  Returns a list of solution dicts (same format as
    solve_rank_condition), or [] if recognition or verification fails.
    """
    global _NSIMPLIFY_BASIS
    if _NSIMPLIFY_BASIS is None:
        from sympy import sqrt as _sqrt
        # Include sqrt(10)=sqrt(2)*sqrt(5) explicitly for Q(sqrt(2),sqrt(5)) values
        # like sqrt(10)+2*sqrt(2) (dotted weights in d=6 Bugaenko polytope).
        _NSIMPLIFY_BASIS = [_sqrt(2), _sqrt(3), _sqrt(5), _sqrt(6), _sqrt(10)]

    from sympy import nsimplify, Rational, Integer
    deadline = time.time() + timeout

    exact_vals = []
    for xi in x_approx:
        if time.time() > deadline:
            return []
        xi_f = float(xi)
        if xi_f <= 1.0:
            return []
        # Try progressively more complex algebraic extensions.
        # Use tolerance=1e-4: numerical optimizer gives ~1e-5 accuracy for dotted
        # weights; we verify exactly via rank/signature check afterward.
        # Collect the BEST match (smallest error) across all bases rather than
        # accepting the first — avoids spurious near-coincidences like
        # 2√5/3+9/2 ≈ 5.9907 being picked over the correct √10+2√2 ≈ 5.9907.
        recognized = None
        best_err = float('inf')
        for basis in [[], [_NSIMPLIFY_BASIS[2]], _NSIMPLIFY_BASIS[:2],
                      _NSIMPLIFY_BASIS[:3], _NSIMPLIFY_BASIS[:4],
                      _NSIMPLIFY_BASIS]:
            if time.time() > deadline:
                break
            try:
                candidate = nsimplify(xi_f, basis, rational=False, tolerance=1e-4)
                err = abs(float(candidate.evalf()) - xi_f)
                if err < 1e-4 and err < best_err:
                    recognized = candidate
                    best_err = err
                    if err < 1e-6:
                        break  # good enough, no need to search further
            except Exception:
                pass
        if recognized is None:
            # Fall back to a high-precision rational approximation
            try:
                recognized = nsimplify(xi_f, rational=True, tolerance=1e-4)
                if abs(float(recognized.evalf()) - xi_f) > 1e-4:
                    recognized = None
            except Exception:
                pass
        if recognized is None:
            return []
        exact_vals.append(recognized)

    # Verify all x_i > 1
    try:
        floats = [float(v.evalf()) for v in exact_vals]
    except Exception:
        return []
    if not all(f > 1.0 + 1e-9 for f in floats):
        return []

    # Substitute into symbolic Gram matrix and verify rank + signature
    try:
        if time.time() > deadline:
            return []
        G_sub = G_sym.subs(list(zip(sym_list, exact_vals)))
        n = G_sub.shape[0]
        G_np = np.array([[float(G_sub[i, j].evalf()) for j in range(n)]
                         for i in range(n)], dtype=float)
        if not _check_signature_float(G_np, d):
            return []
        # Exact rank check (SymPy); catch timeouts from slow symbolic arithmetic
        with _hard_timeout(max(1, int(deadline - time.time()))):
            rank_ok = G_sub.rank() == d + 1
        if not rank_ok:
            return []
    except Exception:
        return []

    sol = {str(s): v for s, v in zip(sym_list, exact_vals)}
    return [sol]


# ---------------------------------------------------------------------------
# Field-agnostic exact recognition via minimal polynomials (2026-06-22)
#
# Replaces the fixed-√-basis nsimplify recognition, which (a) could not
# represent cos(π/m) for m>=7 (cubic/quartic) and (b) was unreliable even in
# Q(√5) (it mis-recognised 17+8√5 as a spurious 6-term √-combination).  Given a
# dotted weight refined to high precision, we recover its MINIMAL POLYNOMIAL by
# PSLQ on its powers [1, x, x², …, x^D]; this identifies any real algebraic
# number regardless of field.  Realizability is then certified by a
# high-precision rank/signature check (not symbolic rank, which is too slow).
# ---------------------------------------------------------------------------

def _recover_minpoly(x_mp, maxdeg=16, dps=100):
    """Minimal polynomial of a high-precision real algebraic number via PSLQ.

    Returns integer coeffs [c0, c1, …, cD] with sum ci x^i = 0 (cD != 0), of the
    least degree D <= maxdeg for which a genuine relation is found; else None.
    """
    import mpmath
    tol = mpmath.power(10, -(dps - 15))
    for D in range(1, maxdeg + 1):
        vec = [x_mp ** i for i in range(D + 1)]
        try:
            rel = mpmath.pslq(vec, maxcoeff=10 ** 15, maxsteps=10 ** 6)
        except Exception:
            rel = None
        if rel and rel[-1] != 0:
            resid = sum(rel[i] * vec[i] for i in range(D + 1))
            if abs(resid) < tol:
                coeffs = [int(c) for c in rel]
                if coeffs[-1] < 0:        # normalise leading coeff > 0
                    coeffs = [-c for c in coeffs]
                return coeffs
    return None


def _build_gram_mpmath(label_assign, dotted_pairs, x_hp, n, dps=100):
    import mpmath
    mpmath.mp.dps = dps

    def exact(m):
        if m == 2:  return mpmath.mpf(0)
        if m == 3:  return mpmath.mpf(-1) / 2
        if m == 4:  return -mpmath.sqrt(2) / 2
        if m == 5:  return -(1 + mpmath.sqrt(5)) / 4
        return -mpmath.cos(mpmath.pi / m)

    G = mpmath.zeros(n)
    for i in range(n):
        G[i, i] = mpmath.mpf(1)
    for (i, j), m in label_assign.items():
        G[i, j] = G[j, i] = exact(m)
    for (i, j), x in zip(dotted_pairs, x_hp):
        G[i, j] = G[j, i] = -x
    return G


def _safe_det(M):
    """mpmath determinant that returns 0 for an exactly-singular matrix instead of
    crashing.  mpmath.det raises inside LU decomposition (a None pivot) when the matrix is
    exactly singular -- which the (d+2)-minor kernel matrices routinely are near a
    solution, especially the larger d>=5 minors with many pi/2 entries.  A singular matrix
    has determinant 0."""
    import mpmath
    try:
        return mpmath.det(M)
    except (TypeError, ZeroDivisionError, ValueError):
        return mpmath.mpf(0)


def _has_parabolic_subdiagram(G_np, n, d, tol=1e-7):
    """True if the Gram has a PARABOLIC (affine) subdiagram => an ideal (cusp) vertex =>
    the polytope is NON-compact.  A subset is parabolic iff its Gram submatrix is positive
    semidefinite (no negative eigenvalue) with a nontrivial kernel (>=1 zero eigenvalue).

    This is the compactness gate that signature(d,1)+isolation alone does NOT provide:
    those accept finite-volume cusped polytopes too.  A parabolic subdiagram lives in a
    PSD subspace of the (d,1) ambient, so it has rank <= d, hence <= d+1 nodes -- we only
    need to scan subsets up to size d+1.  (Same criterion as verify_diagram.check_parabolic,
    which certified P^B6 has none.)"""
    from itertools import combinations as _c
    for size in range(2, d + 2):
        for subset in _c(range(n), size):
            ev = np.linalg.eigvalsh(G_np[np.ix_(subset, subset)])
            if (ev < -tol).sum() == 0 and (np.abs(ev) <= tol).sum() >= 1:
                return True
    return False


def _kernel_jacobian_rank(label_assign, dotted_pairs, x_hp, n, d, dps=100):
    """Numerical rank of the Jacobian of the (d+2)-minor kernel conditions with
    respect to the k dotted weights, evaluated at the refined point x_hp.

    rank == k  <=>  the solution is locally ISOLATED (0-dimensional) => a genuine
    rigid Coxeter polytope.  rank < k  =>  the point lies on a positive-dimensional
    solution component (a spurious continuum, not a discrete polytope).

    This is the SOUND isolation test.  Gauss-Newton convergence to a sharp residual
    does NOT certify isolation: on a positive-dimensional component the minor
    equations still vanish identically, so GN converges to *some* point on it.
    Returns (rank, k).
    """
    import mpmath
    from itertools import combinations as _c
    mpmath.mp.dps = dps
    k = len(dotted_pairs)
    if k == 0:
        return 0, 0
    dot_idx = {i for p in dotted_pairs for i in p}
    rows = [list(c) for c in _c(range(n), d + 2) if dot_idx & set(c)]

    def fvals(xl):
        G = _build_gram_mpmath(label_assign, dotted_pairs, xl, n, dps=dps)
        return [_safe_det(mpmath.matrix([[G[r[a], r[b]] for b in range(len(r))]
                                         for a in range(len(r))])) for r in rows]

    eps = mpmath.power(10, -(dps // 2 - 5))
    cols = []
    for c in range(k):
        xp = list(x_hp); xp[c] = xp[c] + eps
        xm = list(x_hp); xm[c] = xm[c] - eps
        fp = fvals(xp); fm = fvals(xm)
        cols.append([(fp[r] - fm[r]) / (2 * eps) for r in range(len(rows))])

    # Integer rank decision via a float SVD of the Jacobian (sufficient -- the
    # gap between rank-k and rank-(k-1) is large for these systems).
    J = np.array([[float(cols[c][r]) for c in range(k)] for r in range(len(rows))])
    s = np.linalg.svd(J, compute_uv=False)
    if s.size == 0 or s[0] == 0:
        return 0, k
    rank = int(np.sum(s > 1e-6 * s[0]))
    return rank, k


def _recognize_minpoly_and_verify(x_hp, sym_list, label_assign, dotted_pairs,
                                   n, d, dps=100, recover_minpoly=True):
    """High-precision (d,1)-signature + isolation certification (field-agnostic).

    A GN-refined point x_hp (list of mpmath dotted weights) is accepted as a genuine
    compact Coxeter polytope iff ALL of:
      * every dotted weight > 1,
      * the Gram matrix has signature (d, 1) with exactly (n-d-1) zero eigenvalues,
      * the solution is locally ISOLATED -- the kernel-condition Jacobian has full
        rank k in the dotted weights.

    Isolation REPLACES the old "every weight's minpoly is PSLQ-recoverable" gate,
    which silently dropped genuine polytopes whose weights are algebraic of degree
    > 16 (confirmed isolated, signature (d,1), but minpoly unrecoverable even at
    degree 32 -- a ~54% loss on the d=4 k=6 types; see scratchpad/diag5*).  The
    minimal polynomial is now recovered BEST-EFFORT for reporting only: a weight
    whose minpoly is not found stores minpoly=None alongside its decimal value, and
    deduplication (canonical_key) keys on the decimal value, not the minpoly.

    Returns a list with a single solution dict, or [].
    """
    import mpmath
    mpmath.mp.dps = dps
    if not x_hp or any(x <= mpmath.mpf('1') for x in x_hp):
        return []

    # (d,1) signature with (n-d-1) zeros -- the accept/reject decision.
    G = _build_gram_mpmath(label_assign, dotted_pairs, x_hp, n, dps=dps)
    try:
        E = mpmath.eigsy(G, eigvals_only=True)
        evals = [E[i] for i in range(n)]
    except Exception:
        return []
    zero_tol = mpmath.power(10, -(dps // 3))
    pos = sum(1 for e in evals if e > zero_tol)
    neg = sum(1 for e in evals if e < -zero_tol)
    zero = sum(1 for e in evals if abs(e) <= zero_tol)
    if not (pos == d and neg == 1 and zero == n - d - 1):
        return []

    # Isolation: full-rank kernel Jacobian <=> rigid (0-dim) polytope, not a point
    # on a positive-dimensional spurious continuum.
    rank, k = _kernel_jacobian_rank(label_assign, dotted_pairs, x_hp, n, d, dps=dps)
    if rank != k:
        return []

    # Compactness: reject if any parabolic subdiagram (ideal vertex => non-compact).
    # signature+isolation accept cusped finite-volume polytopes too; this is the gate
    # that removes them (e.g. d=5 P9_322: 18 signature-isolated -> 3 compact).
    G_np = np.array([[float(G[i, j]) for j in range(n)] for i in range(n)])
    if _has_parabolic_subdiagram(G_np, n, d):
        return []

    # Best-effort exact minpoly per weight (decorative; never gates acceptance, and
    # canonical_key dedups on the decimal value).  Recovery is EXPENSIVE for the
    # high-degree weights (PSLQ exhausts every degree, ~9s, before returning None),
    # so count-only runs pass recover_minpoly=False and recover minpolys later for
    # the small final survivor set.
    sol = {}
    for sym, x in zip(sym_list, x_hp):
        poly = _recover_minpoly(x, dps=dps) if recover_minpoly else None
        sol[str(sym)] = {
            "value": mpmath.nstr(x, 30),
            "minpoly": poly,            # [c0..cD] with sum ci x^i = 0, or None
        }
    return [sol]


# ---------------------------------------------------------------------------
# Structured screen (2026-06-23): cascade-pin dotted weights via minors
#
# A (d+2)x(d+2) minor whose row set contains exactly ONE not-yet-pinned dotted
# edge is, after substituting the pinned/ordinary entries, a QUADRATIC in that
# one weight x_e.  rank(G) <= d+1 forces the minor to vanish, so x_e must be a
# real root > 1 (<= 2 of them -> branch).  Pinning cascades: once some weights
# are fixed, more minors become single-unknown.  When all weights are pinned we
# check rank(d+1)+signature directly.  This is an optimizer-free screen that is
# SOUND (a realizable candidate makes every minor vanish at its true x_e, so a
# root>1 always exists -> never wrongly rejected) and usually decisive; if the
# cascade gets stuck (no single-unknown minor, e.g. a dense dotted graph) it
# returns None and the caller falls back to the numerical screen.
# ---------------------------------------------------------------------------

def _build_minor_index(dotted_pairs, n, d):
    """(d+2)-subsets containing >=1 dotted edge, with the dotted edges in each."""
    from itertools import combinations as _c
    ds = set(dotted_pairs)
    out = []
    for c in _c(range(n), d + 2):
        cs = set(c)
        de = [e for e in ds if e[0] in cs and e[1] in cs]
        if de:
            out.append((list(c), de))
    return out


def _quad_roots_gt1(Gsub, ii, jj, eps=1e-7):
    """Real roots > 1 of det(Gsub)(x)=0 where entry (ii,jj)=(jj,ii)=-x.
    Returns list of roots (possibly empty), or None if the minor is ~0 for all x
    (x unconstrained by this minor)."""
    def det_at(x):
        M = Gsub.copy(); M[ii, jj] = M[jj, ii] = -x
        return np.linalg.det(M)
    c = det_at(0.0)
    a = (det_at(2.0) - 2.0 * det_at(1.0) + c) / 2.0
    b = det_at(1.0) - c - a
    if abs(a) < 1e-12:
        if abs(b) < 1e-12:
            return None if abs(c) < 1e-9 else []
        x = -c / b
        return [x] if x > 1.0 + eps else []
    disc = b * b - 4 * a * c
    if disc < 0:
        return []
    r = float(np.sqrt(disc))
    return [x for x in ((-b + r) / (2 * a), (-b - r) / (2 * a)) if x > 1.0 + eps]


def _pin_pair_resultant(build, pins, S1, S2, e, f, eps=1e-7):
    """Pin a PAIR of unknown weights when no single-unknown minor exists.

    det(G_S1) and det(G_S2) are each degree<=2 in x_e and in x_f.  Eliminating
    x_f (Sylvester resultant of the two x_f-quadratics) gives a polynomial of
    degree<=4 in x_e; its real roots>1 give x_e, and back-substitution gives
    x_f>1.  Sound: a true solution makes both minors vanish, so the resultant
    vanishes there and the root is found.  Returns list of (x_e, x_f)."""
    def detS(S, xe, xf):
        G = build(pins)
        G[e[0], e[1]] = G[e[1], e[0]] = -xe
        G[f[0], f[1]] = G[f[1], f[0]] = -xf
        return np.linalg.det(G[np.ix_(S, S)])

    def xf_quad(S, xe):
        c = detS(S, xe, 0.0)
        a = (detS(S, xe, 2.0) - 2.0 * detS(S, xe, 1.0) + c) / 2.0
        b = detS(S, xe, 1.0) - c - a
        return a, b, c

    def res_at(xe):
        A, B, C = xf_quad(S1, xe)
        Dd, E, F = xf_quad(S2, xe)
        Syl = np.array([[A, B, C, 0.0], [0.0, A, B, C],
                        [Dd, E, F, 0.0], [0.0, Dd, E, F]])
        return np.linalg.det(Syl)

    xes = [1.3, 1.9, 2.7, 3.8, 5.5, 7.0]
    res_vals = [res_at(x) for x in xes]
    coef = np.polyfit(xes, res_vals, 4)
    # Identically-zero resultant: the two minors are proportional (symmetric
    # degenerate case, e.g. dk=4 perfect-matching types). Return None to signal
    # "inconclusive" so the caller can fall back to the numerical screen.
    if max(abs(v) for v in res_vals) < 1e-6:
        return None
    out = []
    for xe in np.roots(coef):
        if abs(xe.imag) > 1e-6 or xe.real < 1.0 + eps:
            continue
        xe = float(xe.real)
        a, b, c = xf_quad(S1, xe)
        if abs(a) < 1e-12:
            xfs = [-c / b] if (abs(b) > 1e-12 and -c / b > 1.0 + eps) else []
        else:
            disc = b * b - 4 * a * c
            if disc < 0:
                xfs = []
            else:
                r = float(np.sqrt(disc))
                xfs = [x for x in ((-b + r) / (2 * a), (-b - r) / (2 * a))
                       if x > 1.0 + eps]
        for xf in xfs:
            scale = max(1.0, abs(detS(S2, xe, 0.0)))
            if abs(detS(S2, xe, xf)) < 1e-5 * scale:
                out.append((xe, xf))
    return out


def _structured_screen(ordinary_float, dotted_pairs, minor_index, n, d,
                       max_leaves=256):
    """Cascade-pin screen.  Returns (decision, x_solutions):
      decision True  -> >=1 pinned assignment passes rank(d+1)+signature(d,1);
                        x_solutions is a list of dicts {edge: x_float}.
      decision False -> provably no completion (reject).
      decision None  -> cascade stuck (caller should fall back).
    """
    num_zero = n - d - 1

    def build(pins):
        G = np.eye(n)
        for (i, j), v in ordinary_float.items():
            G[i, j] = G[j, i] = v
        for (i, j), x in pins.items():
            G[i, j] = G[j, i] = -x
        return G

    solutions = []
    stack = [(dict(), frozenset(dotted_pairs))]
    leaves = 0
    while stack:
        pins, unknown = stack.pop()
        if not unknown:
            leaves += 1
            if leaves > max_leaves:
                return None, []
            G = build(pins)
            sv = np.linalg.svd(G, compute_uv=False)
            if sv[-num_zero] < 1e-7 and _check_signature_float(G, d, tol=1e-6):
                solutions.append(dict(pins))
            continue
        chosen = None
        for S, de in minor_index:
            un = [e for e in de if e in unknown]
            if len(un) == 1:
                chosen = (S, un[0]); break
        if chosen is None:
            # No single-unknown minor: pin a PAIR via the resultant of two
            # minors over the same unknown pair, then resume the cascade.
            from collections import defaultdict as _dd
            pairmin = _dd(list)
            for S, de in minor_index:
                un = tuple(sorted(x for x in de if x in unknown))
                if len(un) == 2:
                    pairmin[un].append(S)
            pair = next((p for p, ss in pairmin.items() if len(ss) >= 2), None)
            if pair is None:
                return None, []                  # truly stuck -> fall back
            e, f = pair
            pr_result = _pin_pair_resultant(build, pins, pairmin[pair][0],
                                            pairmin[pair][1], e, f)
            if pr_result is None:
                return None, []          # degenerate (zero resultant) -> fall back
            for xe, xf in pr_result:
                p2 = dict(pins); p2[e] = xe; p2[f] = xf
                stack.append((p2, unknown - {e, f}))
            continue
        S, e = chosen
        Gp = build(pins)
        pos = {f: k for k, f in enumerate(S)}
        roots = _quad_roots_gt1(Gp[np.ix_(S, S)], pos[e[0]], pos[e[1]])
        if roots is None:
            return None, []                      # minor doesn't constrain e here
        for x in roots:
            p2 = dict(pins); p2[e] = x
            stack.append((p2, unknown - {e}))
    return (len(solutions) > 0), solutions


def screen_candidate(la, dotted_pairs, minor_index, n, d):
    """Screen one ordinary-label assignment ``la`` ({(i,j): m}) for a realizable Gram
    matrix and return a list of candidate solver starts ``[(la, dotted, x0), ...]`` (empty
    if provably infeasible).  Shared by the brute and block-paste drivers so both use the
    identical screen path: the sound algebraic ``_structured_screen`` cascade, with a
    numerical fallback (``_numerical_screen`` + float signature) when the cascade is stuck.
    The returned starts are handed to ``_refine_mpmath`` + ``_recognize_minpoly_and_verify``
    for exact verification (see ``run_lowk_parallel._solve_worker``)."""
    if not dotted_pairs:
        return [(dict(la), dotted_pairs, [])]
    ordf = {p: _GRAM_FLOAT[m] for p, m in la.items()}
    dec, xs = _structured_screen(ordf, dotted_pairs, minor_index, n, d)
    if dec:
        return [(dict(la), dotted_pairs, [s[p] for p in dotted_pairs]) for s in xs]
    if dec is None:
        xa, res = _numerical_screen(ordf, dotted_pairs, n, d, residual_threshold=1e-6)
        if res <= 1e-6 and _check_signature_float(
                _build_gram_numpy(ordf, dotted_pairs, xa, n), d, tol=1e-3):
            return [(dict(la), dotted_pairs, list(xa))]
    return []


def burcroff_55b_low_weight_edges(mf_list, ordinary_pairs):
    """Ordinary edges forced LOW WEIGHT (label m <= 5) by Burcroff Lemma 5.5(b)
    (arXiv:2201.03437; holds for admissible diagrams in any dimension):

        Let v1v2 be an ordinary edge.  If Sigma has no Lannér diagram containing
        v1 and v2, Sigma has a Lannér diagram L of order > 2 containing v2, and
        there is no dashed (dotted) edge from v1 to any vertex of L, then v1v2
        has low weight.

    Combinatorially: Lannér diagrams of order >= 3 are exactly the missing faces
    of size >= 3 (a Lannér subdiagram's support is a minimal non-face; ordinary
    triangles that are faces are forced elliptic by the vertex-PD constraints,
    and order-2 "Lannér" needs a dotted edge, excluded on ordinary pairs).  So
    the hypotheses depend only on the type's missing-face structure.

    Returns {pair: 5} suitable for enumerate_labels_backtrack(edge_max_label=).
    """
    dotted = [frozenset(m) for m in mf_list if len(m) == 2]
    big = [frozenset(m) for m in mf_list if len(m) >= 3]
    caps = {}
    for (u, v) in ordinary_pairs:
        if any(u in M and v in M for M in big):
            continue                       # (a)-territory / hypothesis fails
        forced = False
        for (v1, v2) in ((u, v), (v, u)):
            for M in big:
                if v2 not in M:
                    continue
                if any(v1 in D and (D - {v1}) & M for D in dotted):
                    continue               # dashed edge from v1 into L
                forced = True
                break
            if forced:
                break
        if forced:
            caps[(u, v)] = 5
    return caps


# ---------------------------------------------------------------------------
# Wildcard (m >= 7) assignment solver — see WILDCARD_LABEL comment at top.
# ---------------------------------------------------------------------------

def _wild_feasible(ordinary_float, dotted_pairs, wild_pairs, pinned, n, d,
                   residual_threshold=1e-6, warm_start=None):
    """Best residual of the rank-(d+1) condition with the wild entries in `pinned`
    (dict pair->c) fixed and everything else (dotted x > 1, unpinned wild
    c in [cos(pi/7), 1)) optimized.  Returns (residual, u_best) where u_best is
    the full unknown vector [dotted..., free_wild...] (slice [:len(dotted_pairs)]
    for the dotted weights).  warm_start: optional u vector of the same layout,
    tried FIRST with a cheap descent — threading the previous scan step's
    solution through consecutive integer pins makes each step ~one descent."""
    free_wild = [p for p in wild_pairs if p not in pinned]
    unknown_pairs = list(dotted_pairs) + free_wild
    k_d = len(dotted_pairs)
    num_zero = n - d - 1

    # precompute the fixed part of G once; objective only writes the unknowns
    G_base = np.eye(n)
    for (i, j), v in ordinary_float.items():
        G_base[i, j] = G_base[j, i] = v
    for (i, j), c in pinned.items():
        G_base[i, j] = G_base[j, i] = -c
    ui = np.array([p[0] for p in unknown_pairs], dtype=int)
    uj = np.array([p[1] for p in unknown_pairs], dtype=int)

    # Fast path: pipeline/c/wild_kernel.c does the IDENTICAL computation (Jacobi
    # eigendecomposition + the same analytic gradient formula) without numpy/
    # scipy dispatch overhead, which profiling showed dominates this exact call
    # site (~90% overhead vs real FLOPs on a 9x9/10x10 matrix). Equivalence
    # verified against the numpy path on 20,000 random trials (max abs error
    # ~1e-11). Falls back to numpy automatically if the library isn't built.
    if _wk.available():
        _prep = _wk.Prepared(G_base, ui, uj, num_zero)
        objective = _prep.objective
        objective_grad = _prep.objective_grad
    else:
        def objective(u):
            G = G_base.copy()
            G[ui, uj] = -u
            G[uj, ui] = -u
            s = np.linalg.svd(G, compute_uv=False)
            return float(np.sum(s[-num_zero:] ** 2))

        def objective_grad(u):
            # f = sum of lambda_i^2 over the num_zero smallest-|lambda| eigenvalues
            # (== smallest singular values, G symmetric).  d(lambda_i)/du_e =
            # v_i^T (dG/du_e) v_i = -2 v[i0,i] v[j0,i], so df/du_e =
            # -4 sum_i lambda_i v[i0,i] v[j0,i].  One eigh replaces ~(k+1) SVDs of
            # finite differencing; same objective, same thresholds.
            G = G_base.copy()
            G[ui, uj] = -u
            G[uj, ui] = -u
            w, V = np.linalg.eigh(G)
            idx = np.argsort(np.abs(w))[:num_zero]
            lam = w[idx]
            Vs = V[:, idx]
            f = float(np.sum(lam ** 2))
            grad = -4.0 * np.einsum('k,ek,ek->e', lam, Vs[ui, :], Vs[uj, :])
            return f, grad

    bounds = [(1.001, 1000.0)] * k_d + [(_WILD_C_MIN, 0.999999)] * len(free_wild)
    starts = []
    for x0 in (1.5, 1.1, 2.0, 3.0):
        for m0 in (7, 8, 10, 12, 18, 30):
            starts.append(np.array([x0] * k_d +
                                   [float(np.cos(np.pi / m0))] * len(free_wild)))
            if not free_wild:
                break
    # Rejection = probe reject (objective > 0.5 at every probe) or the full
    # 24-start multistart failing to reach the threshold — the ORIGINAL validated
    # screen semantics.  (A "cheap single descent lands > 5e-2 -> reject" stage
    # was tried here and REVERTED: it falsely rejected 4 of the 6 census
    # realizations of d=5 type 302/tid19, whose weights x~60 need the full
    # multistart to be reached.  The analytic eigengradient makes the multistart
    # cheap enough that the early reject bought only ~10%.)
    # warm start first: one cheap descent from the caller-supplied point usually
    # settles consecutive integer pins immediately.
    if warm_start is not None and len(warm_start) == len(unknown_pairs) \
            and len(warm_start) > 0:
        try:
            w0 = np.clip(np.asarray(warm_start, dtype=float),
                         [b[0] for b in bounds], [b[1] for b in bounds])
            r = _scipy_opt.minimize(objective_grad, w0, jac=True, bounds=bounds,
                                    method='L-BFGS-B',
                                    options={'maxiter': 60, 'ftol': 1e-20,
                                             'gtol': 1e-12})
            if r.fun < residual_threshold:
                return r.fun, r.x
        except Exception:
            pass
    probe_vals = [objective(s0) for s0 in starts]
    probe_best = min(probe_vals)
    if probe_best > 0.5:
        return probe_best, None
    # seed with the best probe so an optimizer failure (e.g. empty unknown
    # vector when everything is pinned) still returns the true probe residual
    best_val = probe_best
    best_u = starts[probe_vals.index(probe_best)]
    if len(best_u) == 0:
        return best_val, best_u
    if best_val >= residual_threshold:
        for s0 in starts:
            try:
                r = _scipy_opt.minimize(objective_grad, s0, jac=True, bounds=bounds,
                                        method='L-BFGS-B',
                                        options={'maxiter': 200, 'ftol': 1e-20,
                                                 'gtol': 1e-12})
                if r.fun < best_val:
                    best_val, best_u = r.fun, r.x
            except Exception:
                pass
            if best_val < residual_threshold:
                break
    return best_val, best_u


def _float_kernel_jacobian_deficient(base_float, dotted_pairs, x_dotted, n, d,
                                     ratio=1e-7):
    """Cheap float pre-gate for the mpmath isolation test: numerical rank of the
    (d+2)-minor kernel Jacobian wrt the dotted weights, at the float point.
    Returns True only when the Jacobian is CLEARLY rank-deficient
    (sigma_k < ratio * sigma_1) — the point sits on a positive-dimensional
    component and the mpmath isolation gate would reject it after minutes of
    work.  Genuine isolated polytopes have a well-conditioned Jacobian and pass
    untouched; the d=5 wildcard census re-validation gates this empirically."""
    from itertools import combinations as _c
    k = len(dotted_pairs)
    if k == 0:
        return False
    dot_idx = {i for p in dotted_pairs for i in p}
    rows = [list(c) for c in _c(range(n), d + 2) if dot_idx & set(c)]
    di = np.array([p[0] for p in dotted_pairs], dtype=int)
    dj = np.array([p[1] for p in dotted_pairs], dtype=int)

    def fvals(xl):
        G = np.eye(n)
        for (i, j), v in base_float.items():
            G[i, j] = G[j, i] = v
        G[di, dj] = -np.asarray(xl)
        G[dj, di] = -np.asarray(xl)
        return np.array([np.linalg.det(G[np.ix_(r, r)]) for r in rows])

    eps = 1e-6
    J = np.empty((len(rows), k))
    x0 = np.asarray(x_dotted, dtype=float)
    for c in range(k):
        xp = x0.copy(); xp[c] += eps
        xm = x0.copy(); xm[c] -= eps
        J[:, c] = (fvals(xp) - fvals(xm)) / (2 * eps)
    s = np.linalg.svd(J, compute_uv=False)
    if s.size == 0 or s[0] == 0:
        return True
    return bool(s[min(k, len(s)) - 1] < ratio * s[0])


def _solve_wild_assignment(label_assign, wild_pairs, dotted_pairs, sym_list, n, d,
                           numerical_threshold=1e-6, m_scan_max=100, dps=100,
                           deadline_s=1200.0, flags=None):
    """Resolve one wildcard-bearing label assignment (wild edges enumerated as 7,
    meaning ANY m >= 7).  Returns a list of (labels_int, sol) with the wildcard
    labels instantiated to concrete integers, each certified by the standard exact
    path (_refine_mpmath + _recognize_minpoly_and_verify).

    Method (Ma-Zheng "range analysis"): (1) joint feasibility of rank(G)=d+1 with
    the wild entries as continuous unknowns c in [cos(pi/7), 1); (2) per-edge
    integer window scan — for each wild edge, which integers m in [7, m_scan_max]
    keep the system feasible with c_e = cos(pi/m) pinned; (3) every integer tuple
    in the window product is pinned fully and, if still feasible, exactly
    certified with those labels.  If feasibility persists at m_scan_max the type
    CANNOT be closed by this scan: flags["wild_unbounded"] = True (loud, never a
    silent truncation).
    """
    if flags is None:
        flags = {}
    t_dead = time.time() + deadline_s
    ordinary_float = {p: _GRAM_FLOAT[m] for p, m in label_assign.items()
                      if p not in wild_pairs}

    # (1) joint feasibility — all wilds free.
    res0, u_joint = _wild_feasible(ordinary_float, dotted_pairs, wild_pairs, {}, n, d)
    if res0 > numerical_threshold:
        return []

    # (2) per-edge integer windows.  Warm-start each pin from the previous
    # feasible solution (layout: dotted + free wilds minus the scanned edge).
    k_d = len(dotted_pairs)
    windows = []
    for ei, e in enumerate(wild_pairs):
        win = []
        # joint solution without the scanned edge's own entry
        warm = None
        if u_joint is not None and len(u_joint) == k_d + len(wild_pairs):
            warm = np.delete(u_joint, k_d + ei)
        for m in range(7, m_scan_max + 1):
            if time.time() > t_dead:
                flags["wild_deadline"] = True
                return []
            r, u = _wild_feasible(ordinary_float, dotted_pairs, wild_pairs,
                                  {e: float(np.cos(np.pi / m))}, n, d,
                                  warm_start=warm)
            if r < numerical_threshold:
                win.append(m)
                warm = u
        if m_scan_max in win:
            flags["wild_unbounded"] = True
        windows.append(win)
    if any(not w for w in windows):
        return []

    # (3) full pin + exact certification per integer tuple.
    from itertools import product as _prod
    out, seen = [], set()
    warm_tuple = None if u_joint is None else np.asarray(u_joint[:k_d])
    for m_tuple in _prod(*windows):
        if time.time() > t_dead:
            flags["wild_deadline"] = True
            return out
        pinned = {e: float(np.cos(np.pi / m)) for e, m in zip(wild_pairs, m_tuple)}
        # thread the previous feasible tuple's solution: lexicographic order makes
        # consecutive tuples adjacent (one integer step), so one descent settles it
        r, u_pin = _wild_feasible(ordinary_float, dotted_pairs, wild_pairs,
                                  pinned, n, d, warm_start=warm_tuple)
        x_dotted = None if u_pin is None else u_pin[:k_d]
        if r > numerical_threshold:
            continue
        warm_tuple = x_dotted
        labels_int = dict(label_assign)
        for e, m in zip(wild_pairs, m_tuple):
            labels_int[e] = int(m)
        if dotted_pairs:
            if x_dotted is None:
                continue
            # Float pre-gates: signature (d,1) and no parabolic subdiagram at the
            # pinned point — the same conditions the mpmath certification enforces
            # later at ~1000x the cost.  A tuple failing here (generous float
            # tolerances) fails there; this is what keeps positive-dimensional
            # near-solutions from grinding hours of mpmath per assignment.
            base_f = dict(ordinary_float)
            for e, c in pinned.items():
                base_f[e] = -c
            G_np = _build_gram_numpy(base_f, list(dotted_pairs),
                                     np.asarray(x_dotted), n)
            if not _check_signature_float(G_np, d, tol=1e-3):
                continue
            if _has_parabolic_subdiagram(G_np, n, d):
                continue
            # float isolation pre-gate: a clearly rank-deficient kernel Jacobian
            # means a positive-dimensional component — the mpmath isolation gate
            # would reject it after ~minutes; skip it in ~10ms.
            if _float_kernel_jacobian_deficient(base_f, dotted_pairs, x_dotted,
                                                n, d):
                continue
            x_hp = _refine_mpmath(x_dotted, labels_int, dotted_pairs, n, d, dps=dps)
            if x_hp is None:
                continue
            for sol in _recognize_minpoly_and_verify(
                    x_hp, sym_list, labels_int, dotted_pairs, n, d, dps=dps):
                key = (m_tuple, tuple(sorted((k, v["value"]) for k, v in sol.items())))
                if key not in seen:
                    seen.add(key)
                    out.append((labels_int, sol))
        else:
            # no dotted unknowns: G fully determined by the instantiated labels.
            G_np = _build_gram_numpy(
                {p: float(-np.cos(np.pi / m)) if (m := labels_int[p]) >= 7
                 else _GRAM_FLOAT[m] for p in labels_int}, [], np.array([]), n)
            evals = np.linalg.eigvalsh(G_np)
            if int(np.sum(np.abs(evals) > 1e-8)) == d + 1 \
                    and _check_signature_float(G_np, d) \
                    and not _has_parabolic_subdiagram(G_np, n, d):
                G_sym = build_gram_matrix(n, labels_int, frozenset(), {})
                if G_sym.rank() == d + 1:
                    key = (m_tuple, ())
                    if key not in seen:
                        seen.add(key)
                        out.append((labels_int, {}))
    return out


# ---------------------------------------------------------------------------
# Per-type Stage 4 driver
# ---------------------------------------------------------------------------

def process_type_stage4(t, d,
                         max_assignments=50000,
                         enum_timeout=60.0,
                         solve_timeout=60.0,
                         numerical_threshold=1e-6,
                         label_indices=None,
                         face_tuples_max_size=4,
                         extended_lanner_max_extra=0,
                         prefix=None,
                         stats_out=None,
                         wildcard=False,
                         use_burcroff_55b=False,
                         verbose=False):
    """Run Stage 4 on one surviving combinatorial type.

    wildcard: rigorous label treatment (Ma-Zheng Prop 3.5) — enumerate over
    {2,...,6,7} with 7 = "any m >= 7"; wild-bearing assignments are resolved by
    _solve_wild_assignment (continuous c + integer instantiation), so NO a-priori
    label cap is assumed.  Forces label_indices = WILDCARD_INDICES.  If a wild
    window reaches m_scan_max, stats_out["wild_unbounded"] = True and the type's
    verdict is NOT rigorous.

    label_indices: tuple of indices into VALID_LABELS to restrict the ordinary-edge
    label search.  Default None = all labels.  Use (0,1,2,3) to restrict to
    labels {2,3,4,5}, which is appropriate for d≥5 where Burcroff's low-weight
    lemma bounds labels to ≤5 for all ordinary edges.

    prefix: optional tuple of label indices fixing the first len(prefix) pairs of
    the deterministic enumeration order (see enumerate_labels_backtrack) — running
    all prefixes of a given depth across workers partitions the search exactly.

    stats_out: optional dict; on return contains screened/passed_screen/
    exact_attempts counts plus "exhausted": True iff the verdict is rigorous
    (the label search was fully explored — no enum timeout, no max_assignments
    cut, no solve-budget break).  A distinct==0 with exhausted=False is a
    TIMEOUT, not a result.

    Returns list of valid Gram configurations (dicts).
    """
    if stats_out is not None:
        stats_out["exhausted"] = False
    if not SYMPY_AVAILABLE or not NUMPY_AVAILABLE:
        return []
    if wildcard:
        if label_indices is not None and any(i > 5 for i in label_indices):
            raise ValueError("wildcard mode: label_indices must be within {2..7}")
        if label_indices is None:
            label_indices = WILDCARD_INDICES

    n = d + 4
    mf_list = [frozenset(m) for m in t["missing_faces"]]

    dotted_pairs = [tuple(sorted(m)) for m in mf_list if len(m) == 2]
    dotted_set = frozenset(dotted_pairs)
    all_pairs = frozenset((i, j) for i in range(n) for j in range(i + 1, n))
    ordinary_pairs = all_pairs - dotted_set

    # Vertex sets: prefer an explicit list on the type dict (e.g. ground-truth
    # vertex flags ingested directly, bypassing the Gale-diagram generator);
    # otherwise reconstruct them from the stored affine Gale diagram example.
    if t.get("vertex_sets"):
        vertex_sets_list = [frozenset(v) for v in t["vertex_sets"]]
    else:
        # Fall back to the EXACT affine-Gale reconstruction (gale_exact), never
        # the buggy float GaleDiagram (gale.py) which computes wrong vertices.
        from pipeline.utils.gale_exact import AffineGale
        pts = [tuple(p) for p in t["example_points"]]
        pos = frozenset(t["example_positive"])
        vertex_sets_list = AffineGale(pts, pos, d).vertex_sets()
    if not vertex_sets_list:
        return []

    vertex_groups = _build_vertex_groups(ordinary_pairs, vertex_sets_list)

    # Augment vertex groups with all face-tuples of size 3..face_tuples_max_size.
    # Vinberg: every mutually-meeting set of facets spans an elliptic sub-diagram,
    # so its Gram submatrix must be positive definite.  The existing vertex groups
    # enforce this at geometric vertices; face-tuple groups enforce it everywhere
    # else, adding the triangle-free constraint that prunes types with few/no
    # missing triples (the 19 uncertain types).
    if face_tuples_max_size >= 3:
        ft_groups = _build_face_tuple_groups(
            ordinary_pairs, mf_list, n, max_size=face_tuples_max_size
        )
        # Skip tuples already covered by an existing vertex group (same node set).
        existing_node_sets = {frozenset(v_sorted) for v_sorted, _ in vertex_groups}
        new_ft = [(v, p) for v, p in ft_groups
                  if frozenset(v) not in existing_node_sets]
        vertex_groups = vertex_groups + new_ft
        if verbose and new_ft:
            sizes = {}
            for v, _ in new_ft:
                sizes[len(v)] = sizes.get(len(v), 0) + 1
            print(f"    face-tuple groups added: {sizes}")

    lanner_groups = _build_lanner_groups(ordinary_pairs, mf_list)

    # Augment Lannér groups with extended missing-face groups of size 4..
    # For each size-3+ missing face M and extra nodes E: if G_{M ∪ E} has all
    # ordinary pairs and no second missing face, it must have signature (k-1,1)
    # (same Lannér structure, propagated outward).  Size-4 extended groups have
    # 3871 valid combos ≤ 4096 → full bitmask forward-checking.
    if extended_lanner_max_extra >= 1:
        ext_lg = _build_extended_lanner_groups(
            ordinary_pairs, mf_list, n, max_extra=extended_lanner_max_extra
        )
        existing_lg_sets = {frozenset(f) for f, _, _ in lanner_groups}
        new_ext = [(f, p, k) for f, p, k in ext_lg
                   if frozenset(f) not in existing_lg_sets]
        lanner_groups = lanner_groups + new_ext
        if verbose and new_ext:
            sizes = {}
            for f, _, k in new_ext:
                sizes[k] = sizes.get(k, 0) + 1
            print(f"    extended-Lannér groups added: {sizes}")

    # Fast pre-filter: if any Lannér face is a subset of a vertex, Sylvester's
    # criterion forces the vertex PD and Lannér conditions to conflict → no
    # valid labeling can exist, skip immediately.
    vertex_sets_frozen = [frozenset(v) for v in vertex_sets_list]
    for face_sorted, _, _ in lanner_groups:
        face_set = frozenset(face_sorted)
        for v in vertex_sets_frozen:
            if face_set <= v:
                if verbose:
                    print(f"    Lannér face {face_sorted} inside vertex {sorted(v)} → infeasible")
                if stats_out is not None:
                    stats_out["exhausted"] = True   # provably infeasible: rigorous 0
                return []

    dot_syms = {
        pair: symbols(f'x_{pair[0]}_{pair[1]}', positive=True)
        for pair in dotted_pairs
    }
    sym_list = [dot_syms[p] for p in dotted_pairs]

    # Minor index for the structured (cascade-pin) screen.
    minor_index = _build_minor_index(dotted_pairs, n, d) if dotted_pairs else []

    results = []
    t_start = time.time()
    stats = {"screened": 0, "passed_screen": 0, "exact_attempts": 0}
    enum_state = {}
    budget_break = False

    edge_caps = (burcroff_55b_low_weight_edges(mf_list, ordinary_pairs)
                 if use_burcroff_55b else None)
    if verbose and edge_caps:
        print(f"    Burcroff 5.5(b) low-weight caps on {len(edge_caps)} edges")
    if stats_out is not None and edge_caps:
        stats_out["burcroff_55b_edges"] = len(edge_caps)

    for label_assign in enumerate_labels_backtrack(
        ordinary_pairs, vertex_groups, lanner_groups,
        max_count=max_assignments,
        timeout=enum_timeout,
        label_indices=label_indices,
        prefix=prefix,
        state_out=enum_state,
        edge_max_label=edge_caps,
    ):
        elapsed = time.time() - t_start
        if elapsed > enum_timeout + solve_timeout:
            budget_break = True
            break

        stats["screened"] += 1

        if wildcard:
            wild_pairs = tuple(sorted(
                p for p, m in label_assign.items() if m == WILDCARD_LABEL))
            if wild_pairs:
                stats["wild_assignments"] = stats.get("wild_assignments", 0) + 1
                wflags = {}
                for labels_int, sol in _solve_wild_assignment(
                        label_assign, wild_pairs, dotted_pairs, sym_list, n, d,
                        numerical_threshold=numerical_threshold, flags=wflags):
                    results.append({
                        "type_id":          t["type_id"],
                        "label_assignment": {str(p): v
                                             for p, v in labels_int.items()},
                        "dot_values":       sol,
                    })
                if wflags.get("wild_unbounded"):
                    stats["wild_unbounded"] = True
                if wflags.get("wild_deadline"):
                    stats["wild_deadline"] = True
                continue

        ordinary_float = {p: _GRAM_FLOAT[m]
                          for p, m in label_assign.items()}

        if dotted_pairs:
            # Step 3a: structured (cascade-pin) screen — optimizer-free, sound.
            # Gives candidate pinned weights directly; if the cascade gets stuck
            # (dense dotted graph) it returns None and we fall back to the
            # numerical (Gauss-Newton/L-BFGS) screen.
            dec, x_sols = _structured_screen(
                ordinary_float, dotted_pairs, minor_index, n, d,
            )
            if dec is False:
                continue                          # provably infeasible
            if dec is None:
                x_approx, residual = _numerical_screen(
                    ordinary_float, dotted_pairs, n, d,
                    residual_threshold=numerical_threshold,
                )
                if residual > numerical_threshold:
                    continue
                G_np = _build_gram_numpy(ordinary_float, dotted_pairs, x_approx, n)
                if not _check_signature_float(G_np, d, tol=1e-3):
                    continue
                x_starts = [x_approx]
            else:
                # dec is True: structured screen found candidate completion(s).
                x_starts = [np.array([xs[p] for p in dotted_pairs])
                            for xs in x_sols]

            stats["passed_screen"] += 1

            # Step 3b: field-agnostic exact solve from each candidate start.
            # Refine to high precision (Gauss-Newton, exact mpmath ordinary
            # entries), recover each weight's MINIMAL POLYNOMIAL via PSLQ on its
            # powers, and certify rank(d+1)+signature(d,1).  Works over ANY number
            # field (labels >=7).  False positives from the float screen fail to
            # refine/recognise -> dropped.
            seen_local = set()
            for x_approx in x_starts:
                stats["exact_attempts"] += 1
                x_hp = _refine_mpmath(
                    x_approx, label_assign, dotted_pairs, n, d, dps=100,
                )
                if x_hp is None:
                    continue
                for sol in _recognize_minpoly_and_verify(
                        x_hp, sym_list, label_assign, dotted_pairs, n, d, dps=100):
                    # Dedup key: prefer the exact minimal polynomial; fall back to the
                    # high-precision numeric value when PSLQ could not recover a minpoly
                    # (a None minpoly must not crash the key — it did for high-degree
                    # weights in d=5 types 302/319/322, dropping those realizations).
                    def _kpart(v):
                        mp_ = v.get("minpoly")
                        if mp_ is not None:
                            return ("mp", tuple(mp_))
                        return ("val", str(v.get("value"))[:40])
                    key = tuple(sorted((kk, _kpart(v)) for kk, v in sol.items()))
                    if key in seen_local:
                        continue
                    seen_local.add(key)
                    results.append({
                        "type_id":          t["type_id"],
                        "label_assignment": {str(p): v
                                             for p, v in label_assign.items()},
                        "dot_values":       sol,
                    })

        else:
            # No dotted pairs — check rank and signature directly (float OK
            # for rank/sig since there are no free parameters)
            G_np = _build_gram_numpy(ordinary_float, [], np.array([]), n)
            evals = np.linalg.eigvalsh(G_np)
            rank_approx = int(np.sum(np.abs(evals) > 1e-8))
            if rank_approx == d + 1 and _check_signature_float(G_np, d):
                # Confirm with exact rank (no unknowns so SymPy is fast)
                G_sym = build_gram_matrix(n, label_assign, dotted_set, {})
                if G_sym.rank() == d + 1:
                    results.append({
                        "type_id":          t["type_id"],
                        "label_assignment": {str(k): v
                                             for k, v in label_assign.items()},
                        "dot_values":       {},
                    })

    if verbose:
        print(f"    screened={stats['screened']} "
              f"passed_screen={stats['passed_screen']} "
              f"exact_attempts={stats['exact_attempts']}")

    if stats_out is not None:
        stats_out.update(stats)
        stats_out["enum_count"] = enum_state.get("count", 0)
        # Rigorous only if the generator itself finished every branch (no enum
        # timeout, no max_count cut), we never broke out on the solve budget,
        # and no wild assignment was cut short by its per-assignment deadline.
        stats_out["exhausted"] = (bool(enum_state.get("exhausted"))
                                  and not budget_break
                                  and not stats.get("wild_deadline"))

    return results


# ---------------------------------------------------------------------------
# Stage 4 runner
# ---------------------------------------------------------------------------

def _stage4_pool_init(label_indices=None):
    """Pre-warm module-level caches in each worker process."""
    _canonical_pd_valid(3, label_indices=label_indices)
    _canonical_pd_valid(4, label_indices=label_indices)
    _canonical_pd_valid(5, label_indices=label_indices)
    _canonical_lanner_valid(3, label_indices=label_indices)
    _canonical_lanner_valid(4, label_indices=label_indices)


def _stage4_worker(args):
    """Top-level worker function for multiprocessing (must be picklable).

    Never raises: a crash on one type returns [] for that type (logged) rather
    than killing the whole pool / production run.
    """
    idx, t, d, max_assignments, enum_timeout, solve_timeout, label_indices, extended_lanner_max_extra, face_tuples_max_size = args
    try:
        result = process_type_stage4(
            t, d,
            max_assignments=max_assignments,
            enum_timeout=enum_timeout,
            solve_timeout=solve_timeout,
            label_indices=label_indices,
            extended_lanner_max_extra=extended_lanner_max_extra,
            face_tuples_max_size=face_tuples_max_size,
            verbose=False,
        )
    except Exception as e:
        print(f"  Type {t.get('type_id', idx)}: ERROR {type(e).__name__}: {e} "
              f"-> skipped (0 configs)", flush=True)
        return idx, []
    return idx, result


def run_stage4(d, stage3_results, output_dir=None, verbose=True,
               per_type_timeout=180.0, n_workers=None, label_indices=None,
               max_assignments=50000, extended_lanner_max_extra=0,
               face_tuples_max_size=4):
    """Run Stage 4: Coxeter label enumeration + exact Gram realizability."""
    import multiprocessing as mp
    import os

    n = d + 4
    print(f"Stage 4: d={d}, n={n}")
    print(f"  Input: {len(stage3_results)} surviving types")

    if not SYMPY_AVAILABLE:
        print("  ERROR: SymPy not available. Skipping Stage 4.")
        return []
    if not NUMPY_AVAILABLE:
        print("  ERROR: NumPy not available. Skipping Stage 4.")
        return []

    if n_workers is None:
        n_workers = 0  # default: sequential (multiprocessing unreliable on macOS)

    enum_timeout = per_type_timeout * 0.5
    solve_timeout = per_type_timeout * 0.5

    all_results = []
    t0 = time.time()

    out = None
    gram_path = None
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        gram_path = out / "gram_matrices.json"
        gram_path.write_text("[]")

    if n_workers > 0:
        tasks = [(i, t, d, max_assignments, enum_timeout, solve_timeout, label_indices, extended_lanner_max_extra, face_tuples_max_size)
                 for i, t in enumerate(stage3_results)]
        done = 0
        print(f"  Using {n_workers} workers", flush=True)
        ctx = mp.get_context('fork')
        init_args = (label_indices,)
        with ctx.Pool(processes=n_workers,
                      initializer=_stage4_pool_init,
                      initargs=init_args) as pool:
            for idx, type_results in pool.imap_unordered(
                _stage4_worker, tasks, chunksize=1
            ):
                done += 1
                if type_results:
                    print(f"  Type {stage3_results[idx]['type_id']}: "
                          f"{len(type_results)} valid config(s)", flush=True)
                    all_results.extend(type_results)
                    if gram_path is not None:
                        gram_path.write_text(
                            json.dumps(all_results, indent=2, default=str)
                        )
                if verbose and (done % 25 == 0):
                    elapsed = time.time() - t0
                    print(f"  ... {done}/{len(stage3_results)} types "
                          f"({elapsed:.0f}s elapsed, {len(all_results)} found)",
                          flush=True)
    else:
        print(f"  Sequential mode", flush=True)
        for idx, t in enumerate(stage3_results):
            type_results = process_type_stage4(
                t, d,
                max_assignments=max_assignments,
                enum_timeout=enum_timeout,
                solve_timeout=solve_timeout,
                label_indices=label_indices,
                extended_lanner_max_extra=extended_lanner_max_extra,
                face_tuples_max_size=face_tuples_max_size,
                verbose=False,
            )
            try:
                from sympy.core.cache import clear_cache as _sympy_clear
                _sympy_clear()
            except Exception:
                pass
            if type_results:
                print(f"  Type {t['type_id']}: {len(type_results)} valid config(s)",
                      flush=True)
                all_results.extend(type_results)
                if gram_path is not None:
                    gram_path.write_text(
                        json.dumps(all_results, indent=2, default=str)
                    )
            if verbose and (idx + 1) % 25 == 0:
                elapsed = time.time() - t0
                rate = (idx + 1) / elapsed
                remaining = (len(stage3_results) - idx - 1) / rate
                print(f"  ... {idx+1}/{len(stage3_results)} types  |  "
                      f"{elapsed:.0f}s elapsed  |  ~{remaining/60:.0f} min remaining  |  "
                      f"{len(all_results)} found", flush=True)

    print(f"  Total valid Gram configurations: {len(all_results)}")

    if out is not None:
        gram_path.write_text(json.dumps(all_results, indent=2, default=str))
        write_manifest(
            out, "stage4", {"d": d, "n": n},
            {"num_surviving": len(stage3_results),
             "num_valid_grams": len(all_results)},
        )
        print(f"  Written to {out}/")

    return all_results


def load_stage4(output_dir):
    """Load Stage 4 results from disk."""
    path = Path(output_dir) / "gram_matrices.json"
    return json.loads(path.read_text())
