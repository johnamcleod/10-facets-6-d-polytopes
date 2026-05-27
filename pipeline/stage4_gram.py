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
import time
from pathlib import Path

from pipeline.utils.manifest import write_manifest

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
    2: 0.0,
    3: -0.5,
    4: -0.7071067811865476,   # -√2/2
    5: -0.8090169943749474,   # -(1+√5)/4
}

LABEL_CAP = 5
VALID_LABELS = [2, 3, 4, 5]


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
    bounds = [(1.001, 30.0)] * k
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


def _check_signature_float(G_numpy, d):
    """Check signature (d,1): exactly one negative eigenvalue."""
    evals = np.linalg.eigvalsh(G_numpy)
    neg = int(np.sum(evals < -1e-8))
    pos = int(np.sum(evals > 1e-8))
    return neg == 1 and pos == d


# ---------------------------------------------------------------------------
# Backtracking label enumeration
# ---------------------------------------------------------------------------

def _build_vertex_groups(ordinary_pairs_set, vertex_sets_list,
                          sub_size_threshold=4):
    """Return list of (v_sorted_tuple, vertex_ordinary_pairs_frozenset).

    For large vertices (k pairs where 4^k > 2^20), also include sub-vertex
    groups of size sub_size_threshold to get precomputable forward-checking
    constraints (sub-matrices of a PD matrix are also PD).
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

    # Lannér pairs — process first so Lannér checks fire early
    lanner_pairs = set()
    for _, lp, _ in lanner_groups:
        lanner_pairs.update(lp)
    for p in sorted(lanner_pairs & ordinary_set):
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
#   or None if 4^C(d,2) > max_combos.
# Valid combos are the same for every group of the same size (canonical structure).
_VG_VALID_CACHE: dict = {}

# Lannér cache: key = k (group size), value = valid combos array.
_LG_VALID_CACHE: dict = {}


def _canonical_pd_valid(d, max_combos=1 << 20):
    """Compute valid label-index combos for a d-node all-ordinary vertex (PD condition).

    Returns numpy array of shape (N, k) dtype int8, or None if 4^k > max_combos.
    Cached globally — same result for every vertex group of the same size.
    """
    if d in _VG_VALID_CACHE:
        return _VG_VALID_CACHE[d]

    k = d * (d - 1) // 2  # = C(d,2) pairs, in lex order
    if 4 ** k > max_combos:
        _VG_VALID_CACHE[d] = None
        return None

    label_f = [_GRAM_FLOAT[m] for m in VALID_LABELS]
    # Canonical pair ordering: (0,1),(0,2),...,(0,d-1),(1,2),...
    pair_idx = 0
    pos_map = {}
    for a in range(d):
        for b in range(a + 1, d):
            pos_map[pair_idx] = (a, b)
            pair_idx += 1

    valid_rows = []
    for combo in itertools.product(range(4), repeat=k):
        M = np.eye(d, dtype=float)
        for pi, (a, b) in pos_map.items():
            v = label_f[combo[pi]]
            M[a, b] = v
            M[b, a] = v
        try:
            np.linalg.cholesky(M)
            valid_rows.append(combo)
        except np.linalg.LinAlgError:
            pass

    arr = np.array(valid_rows, dtype=np.int8) if valid_rows else np.zeros((0, k), dtype=np.int8)
    _VG_VALID_CACHE[d] = arr
    return arr


def _canonical_lanner_valid(k, max_combos=1 << 20):
    """Compute valid label-index combos for a k-node all-ordinary Lannér group.

    Returns numpy array of shape (N, C(k,2)) dtype int8, or None if too large.
    Cached globally.
    """
    if k in _LG_VALID_CACHE:
        return _LG_VALID_CACHE[k]

    num_pairs = k * (k - 1) // 2
    if 4 ** num_pairs > max_combos:
        _LG_VALID_CACHE[k] = None
        return None

    label_f = [_GRAM_FLOAT[m] for m in VALID_LABELS]
    pos_map = {}
    pi = 0
    for a in range(k):
        for b in range(a + 1, k):
            pos_map[pi] = (a, b)
            pi += 1

    valid_rows = []
    for combo in itertools.product(range(4), repeat=num_pairs):
        M = np.eye(k, dtype=float)
        for pi, (a, b) in pos_map.items():
            v = label_f[combo[pi]]
            M[a, b] = v
            M[b, a] = v
        evals = np.linalg.eigvalsh(M)
        if int(np.sum(evals < -1e-8)) == 1:
            valid_rows.append(combo)

    arr = np.array(valid_rows, dtype=np.int8) if valid_rows else np.zeros((0, num_pairs), dtype=np.int8)
    _LG_VALID_CACHE[k] = arr
    return arr


def _prepare_vertex_precomputed(vertex_groups):
    """Prepare forward-checking data for all vertex groups.

    For groups with 4^k ≤ 2^20 (k = C(d,2) pairs), loads from the global
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
        valid_arr_full = _canonical_pd_valid(d)
        if valid_arr_full is not None and len(v_pairs) == d * (d - 1) // 2:
            # Select only the columns corresponding to v_pairs, in v_pairs order
            valid_arr = valid_arr_full[:, canonical_idx]
        elif valid_arr_full is not None and len(v_pairs) < d * (d - 1) // 2:
            # Group has fewer pairs than C(d,2) — some pairs are dotted (unexpected
            # for ordinary vertices, but handle gracefully)
            valid_arr = valid_arr_full[:, canonical_idx]
        else:
            valid_arr = None  # too large to precompute

        pair_to_pos = {p: i for i, p in enumerate(v_pairs)}
        result.append((v_sorted, v_pairs, pair_to_pos, valid_arr))
    return result


def _prepare_lanner_precomputed(lanner_groups):
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

        valid_arr_full = _canonical_lanner_valid(k)
        if valid_arr_full is not None:
            valid_arr = valid_arr_full[:, canonical_idx]
        else:
            valid_arr = None

        pair_to_pos = {p: i for i, p in enumerate(l_pairs)}
        result.append((face_sorted, l_pairs, pair_to_pos, k, valid_arr))
    return result


def enumerate_labels_backtrack(ordinary_pairs, vertex_groups, lanner_groups=None,
                                max_count=50000, timeout=60.0):
    """Generator: backtracking with forward-checking using precomputed vertex
    and Lannér valid-assignment tables.

    For each group with 4^k ≤ 2^20 pairs, valid assignments are precomputed
    once and forward-checked at every assignment step via numpy boolean masking.
    Groups too large to precompute fall back to on-the-fly PD checking when
    the group is fully assigned.

    Yields dicts {(i,j): m} for each valid assignment.
    """
    if lanner_groups is None:
        lanner_groups = []

    pairs_list = _vertex_first_ordering(ordinary_pairs, vertex_groups, lanner_groups)
    n_pairs = len(pairs_list)
    pair_idx = {p: i for i, p in enumerate(pairs_list)}

    # Precompute valid arrays for all groups
    vg_precomp = _prepare_vertex_precomputed(vertex_groups)
    lg_precomp = _prepare_lanner_precomputed(lanner_groups)

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
        filters = [[0] * 4 for _ in range(k_pairs)]
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
        filters = [[0] * 4 for _ in range(k_pairs)]
        for i in range(n_combos):
            bit = 1 << i
            for pos in range(k_pairs):
                filters[pos][valid_arr[i, pos]] |= bit
        filter_lg.append(filters)

    # Current valid masks as Python integers (bitmask over valid-combo indices).
    vg_masks = [full_vg[g] for g in range(len(vg_precomp))]
    lg_masks = [full_lg[g] for g in range(len(lg_precomp))]

    assignment = [-1] * n_pairs
    state = {"count": 0, "start": time.time()}
    trail = []  # (g_idx, is_lanner, old_int_mask)

    def _backtrack(depth):
        if state["count"] >= max_count or time.time() - state["start"] > timeout:
            return
        if depth == n_pairs:
            state["count"] += 1
            yield {pairs_list[i]: VALID_LABELS[assignment[i]] for i in range(n_pairs)}
            return

        pi = depth
        trail_mark = len(trail)

        for li in range(4):   # label index 0..3
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
                assign_dict = {pairs_list[i]: VALID_LABELS[assignment[i]]
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
    else:
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

    try:
        if k == 1:
            x = sym_list[0]
            roots = solve(eqs[0], x)
            raw = [{x: r} for r in roots if r.is_real]
            return _validate_solutions(raw, sym_list, G, d, deadline)

        elif k <= 4:
            try:
                raw = solve(eqs[: k + 1], sym_list, dict=True)
            except Exception:
                raw = []
            valid = _validate_solutions(raw, sym_list, G, d, deadline)
            if valid or time.time() > deadline:
                return valid
            # Gröbner fallback
            try:
                basis = groebner(eqs[: min(len(eqs), k + 3)], *sym_list,
                                 order='lex')
                raw2 = solve(list(basis), sym_list, dict=True)
                return _validate_solutions(raw2, sym_list, G, d, deadline)
            except Exception:
                return valid

        else:
            try:
                basis = groebner(eqs[: min(len(eqs), k + 3)], *sym_list,
                                 order='lex')
                raw = solve(list(basis), sym_list, dict=True)
                return _validate_solutions(raw, sym_list, G, d, deadline)
            except Exception:
                return []

    except Exception:
        return []


# ---------------------------------------------------------------------------
# Per-type Stage 4 driver
# ---------------------------------------------------------------------------

def process_type_stage4(t, d,
                         max_assignments=50000,
                         enum_timeout=60.0,
                         solve_timeout=60.0,
                         numerical_threshold=1e-6,
                         verbose=False):
    """Run Stage 4 on one surviving combinatorial type.

    Returns list of valid Gram configurations (dicts).
    """
    if not SYMPY_AVAILABLE or not NUMPY_AVAILABLE:
        return []

    n = d + 4
    mf_list = [frozenset(m) for m in t["missing_faces"]]

    dotted_pairs = [tuple(sorted(m)) for m in mf_list if len(m) == 2]
    dotted_set = frozenset(dotted_pairs)
    all_pairs = frozenset((i, j) for i in range(n) for j in range(i + 1, n))
    ordinary_pairs = all_pairs - dotted_set

    from pipeline.utils.gale import GaleDiagram
    pts = [tuple(p) for p in t["example_points"]]
    pos = frozenset(t["example_positive"])
    gd = GaleDiagram(pts, pos, d)
    vertex_sets_list = gd.vertex_sets()
    if not vertex_sets_list:
        return []

    vertex_groups = _build_vertex_groups(ordinary_pairs, vertex_sets_list)
    lanner_groups = _build_lanner_groups(ordinary_pairs, mf_list)

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
                return []

    dot_syms = {
        pair: symbols(f'x_{pair[0]}_{pair[1]}', positive=True)
        for pair in dotted_pairs
    }

    results = []
    t_start = time.time()
    stats = {"screened": 0, "passed_screen": 0, "exact_attempts": 0}

    for label_assign in enumerate_labels_backtrack(
        ordinary_pairs, vertex_groups, lanner_groups,
        max_count=max_assignments,
        timeout=enum_timeout,
    ):
        elapsed = time.time() - t_start
        if elapsed > enum_timeout + solve_timeout:
            break

        stats["screened"] += 1
        ordinary_float = {p: _GRAM_FLOAT[m]
                          for p, m in label_assign.items()}

        if dotted_pairs:
            # Step 3a: numerical screen
            x_approx, residual = _numerical_screen(
                ordinary_float, dotted_pairs, n, d,
                residual_threshold=numerical_threshold,
            )
            if residual > numerical_threshold:
                continue  # numerically infeasible

            # Quick signature check at the numerical solution
            G_np = _build_gram_numpy(ordinary_float, dotted_pairs, x_approx, n)
            if not _check_signature_float(G_np, d):
                continue

            stats["passed_screen"] += 1

            # Step 3b: exact symbolic solve
            stats["exact_attempts"] += 1
            G_sym = build_gram_matrix(n, label_assign, dotted_set, dot_syms)
            remaining = max(10.0,
                            enum_timeout + solve_timeout - (time.time() - t_start))
            solutions = solve_rank_condition(
                G_sym, d, dot_syms, dotted_pairs, n,
                timeout=min(solve_timeout, remaining),
            )
            for sol in solutions:
                results.append({
                    "type_id":          t["type_id"],
                    "label_assignment": {str(k): v
                                         for k, v in label_assign.items()},
                    "dot_values":       {k: str(v) for k, v in sol.items()},
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

    return results


# ---------------------------------------------------------------------------
# Stage 4 runner
# ---------------------------------------------------------------------------

def run_stage4(d, stage3_results, output_dir=None, verbose=True,
               per_type_timeout=180.0):
    """Run Stage 4: Coxeter label enumeration + exact Gram realizability."""
    n = d + 4
    print(f"Stage 4: d={d}, n={n}")
    print(f"  Input: {len(stage3_results)} surviving types")

    if not SYMPY_AVAILABLE:
        print("  ERROR: SymPy not available. Skipping Stage 4.")
        return []
    if not NUMPY_AVAILABLE:
        print("  ERROR: NumPy not available. Skipping Stage 4.")
        return []

    all_results = []
    t0 = time.time()

    for idx, t in enumerate(stage3_results):
        type_results = process_type_stage4(
            t, d,
            max_assignments=50000,
            enum_timeout=per_type_timeout * 0.5,
            solve_timeout=per_type_timeout * 0.5,
            verbose=verbose,
        )
        if type_results:
            print(f"  Type {t['type_id']}: {len(type_results)} valid config(s)")
            all_results.extend(type_results)
        elif verbose and (idx % 25 == 0):
            elapsed = time.time() - t0
            print(f"  ... {idx}/{len(stage3_results)} types "
                  f"({elapsed:.0f}s elapsed, {len(all_results)} found)")

    print(f"  Total valid Gram configurations: {len(all_results)}")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "gram_matrices.json").write_text(
            json.dumps(all_results, indent=2, default=str)
        )
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
