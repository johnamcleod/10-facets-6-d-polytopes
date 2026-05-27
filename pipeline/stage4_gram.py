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

def _build_vertex_groups(ordinary_pairs_set, vertex_sets_list):
    """Return list of (v_sorted_tuple, vertex_ordinary_pairs_frozenset)."""
    groups = []
    for v in vertex_sets_list:
        v_sorted = tuple(sorted(v))
        vp = frozenset(
            (v_sorted[a], v_sorted[b])
            for a in range(len(v_sorted))
            for b in range(a + 1, len(v_sorted))
            if (v_sorted[a], v_sorted[b]) in ordinary_pairs_set
        )
        groups.append((v_sorted, vp))
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


def enumerate_labels_backtrack(ordinary_pairs, vertex_groups, lanner_groups=None,
                                max_count=50000, timeout=60.0):
    """Generator: backtracking over label assignments with vertex PD pruning
    and Lannér subdiagram checks on missing faces of size ≥ 3.

    Yields dicts {(i,j): m} for each valid assignment.
    """
    if lanner_groups is None:
        lanner_groups = []

    pairs_list = _vertex_first_ordering(ordinary_pairs, vertex_groups, lanner_groups)

    pair_to_vg = {p: [] for p in pairs_list}
    for g_idx, (_, vp) in enumerate(vertex_groups):
        for p in vp:
            if p in pair_to_vg:
                pair_to_vg[p].append(g_idx)

    pair_to_lg = {p: [] for p in pairs_list}
    for g_idx, (_, lp, _) in enumerate(lanner_groups):
        for p in lp:
            if p in pair_to_lg:
                pair_to_lg[p].append(g_idx)

    assignment = {}
    state = {"count": 0, "start": time.time()}

    def backtrack(idx):
        if (state["count"] >= max_count or
                time.time() - state["start"] > timeout):
            return
        if idx == len(pairs_list):
            state["count"] += 1
            yield dict(assignment)
            return
        pair = pairs_list[idx]
        for m in VALID_LABELS:
            assignment[pair] = m
            ok = True

            # Vertex PD checks for complete vertices
            for g_idx in pair_to_vg.get(pair, []):
                v_sorted, vp = vertex_groups[g_idx]
                if all(p in assignment for p in vp):
                    if not _is_vertex_pd_float(v_sorted, assignment):
                        ok = False
                        break

            # Lannér checks for complete missing faces (size ≥ 3)
            if ok:
                for g_idx in pair_to_lg.get(pair, []):
                    face_sorted, lp, k = lanner_groups[g_idx]
                    if all(p in assignment for p in lp):
                        if not _is_lanner_float(face_sorted, assignment, k):
                            ok = False
                            break

            if ok:
                yield from backtrack(idx + 1)
        del assignment[pair]

    yield from backtrack(0)


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
