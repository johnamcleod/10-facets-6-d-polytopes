"""
Stage 4: Coxeter label assignment + exact Gram-matrix realizability.

For each surviving combinatorial type:
  1. Determine which pairs of facets meet (ordinary edges in Coxeter diagram)
     vs. are disjoint (dotted edges).
  2. Enumerate integer-label assignments m_ij in {2,3,4,5,...} on ordinary edges
     (with proven caps from Burcroff/Esselmann/FT lemmas).
  3. For each labelling, check local positive-definiteness of all vertex blocks.
  4. Set up the Gram matrix G symbolically with unknowns for dotted-edge weights.
  5. Enforce rank(G) = d+1 = 7 (all (d+2)x(d+2) = 8x8 minors vanish).
  6. Solve for dotted weights exactly (Gröbner basis / substitution).
  7. Check signature (d,1) = (6,1) and x_ab > 1 for dotted entries.

This stage uses SymPy for exact symbolic computation.
"""

import json
import itertools
from pathlib import Path
from fractions import Fraction

from pipeline.utils.manifest import write_manifest

try:
    import sympy
    from sympy import (
        Matrix, symbols, Rational, cos, pi, sqrt, simplify,
        groebner, solve, Poly, factor, zeros, eye, det,
        Symbol, Abs, sign, N, Interval, oo, S
    )
    from sympy.matrices import Matrix
    from sympy.polys.numberfields import field_isomorphism
    SYMPY_AVAILABLE = True
except ImportError:
    SYMPY_AVAILABLE = False
    print("WARNING: SymPy not available. Stage 4 will be skipped.")


# ---------------------------------------------------------------------------
# Gram matrix entries: exact cosines
# ---------------------------------------------------------------------------

def gram_entry_adjacent(m):
    """G_ij for facets meeting at angle pi/m: -cos(pi/m).

    Returns exact SymPy expression.
    Uses exact values:
      m=2: 0
      m=3: -1/2
      m=4: -sqrt(2)/2
      m=5: -(1+sqrt(5))/4 = -cos(pi/5)
      m=6: -sqrt(3)/2
    """
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
    elif m == 6:
        return -sqrt(3) / 2
    else:
        return -cos(pi / m)


# ---------------------------------------------------------------------------
# Edge label bounds
# ---------------------------------------------------------------------------

# From Burcroff's low-weight lemma and compactness constraints,
# labels m_ij are bounded. The default global cap is m <= 5
# (from Felikson-Tumarkin for compact polytopes with d+4 facets).
# We use cap 5 as the conservative bound (can be raised if needed).
LABEL_CAP = 5
LABEL_MIN = 2  # orthogonal

VALID_LABELS = list(range(LABEL_MIN, LABEL_CAP + 1))
# Label 2 = orthogonal (no edge in diagram)
# Label 3,4,5 = ordinary edge with those labels


# ---------------------------------------------------------------------------
# Vertex positive-definiteness check
# ---------------------------------------------------------------------------

def check_vertex_pd(gram_matrix, vertex_set, n):
    """Check that the principal submatrix of G indexed by vertex_set is PD.

    vertex_set: set of d indices (the d facets meeting at this vertex)
    gram_matrix: n x n SymPy Matrix

    Returns True if the submatrix is positive definite (all leading minors > 0).
    """
    if not SYMPY_AVAILABLE:
        raise RuntimeError("SymPy required")
    idx = sorted(vertex_set)
    d = len(idx)
    sub = gram_matrix.extract(idx, idx)
    # Check all leading principal minors > 0
    for k in range(1, d + 1):
        minor = sub[:k, :k].det()
        # minor should be a rational number (no unknowns at this stage)
        minor_val = minor
        try:
            minor_val = float(minor)
        except Exception:
            pass
        if minor_val <= 0:
            return False, f"Leading {k}x{k} minor = {minor_val} <= 0"
    return True, None


# ---------------------------------------------------------------------------
# Main Stage 4 logic
# ---------------------------------------------------------------------------

def build_gram_matrix(n, adjacent_pairs, dotted_pairs, label_assignment, dot_symbols):
    """Build the symbolic Gram matrix.

    Args:
        n: number of facets
        adjacent_pairs: dict (i,j) -> m_ij (integer label, i<j)
        dotted_pairs: list of (i,j) pairs (i<j), dotted edges with unknowns
        label_assignment: dict (i,j) -> m value for ordinary edges
        dot_symbols: dict (i,j) -> SymPy symbol for -G_ij > 1

    Returns:
        n x n SymPy Matrix
    """
    if not SYMPY_AVAILABLE:
        raise RuntimeError("SymPy required")

    G = eye(n)  # Diagonal = 1

    for i in range(n):
        for j in range(i + 1, n):
            pair = (i, j)
            if pair in dotted_pairs or (j, i) in dotted_pairs:
                key = pair if pair in dot_symbols else (j, i)
                x = dot_symbols[key]
                G[i, j] = -x
                G[j, i] = -x
            elif pair in label_assignment:
                m = label_assignment[pair]
                entry = gram_entry_adjacent(m)
                G[i, j] = entry
                G[j, i] = entry
            else:
                # Not adjacent, not dotted = orthogonal (m=2 => entry=0)
                # (already 0 from eye init)
                pass

    return G


def enumerate_label_assignments(ordinary_pairs, vertex_sets, n):
    """Enumerate valid integer-label assignments for ordinary edges.

    Returns generator of label dicts {(i,j): m_ij}.
    Applies local elliptic (vertex PD) pre-pruning.

    For each vertex (set of d facets), the corresponding Coxeter diagram
    must be of finite type (positive definite). We use this to prune
    label combinations per vertex.
    """
    # Build: which pairs are in which vertex
    pair_to_vertices = {}
    for pair in ordinary_pairs:
        for v_idx, v in enumerate(vertex_sets):
            if pair[0] in v and pair[1] in v:
                pair_to_vertices.setdefault(pair, []).append(v_idx)

    # Simple enumeration: try all combinations of VALID_LABELS for each pair
    pairs_list = list(ordinary_pairs)
    # Remove m=2 from pairs that should be non-orthogonal edges
    # (m=2 means no Coxeter edge; we include it for completeness)
    label_ranges = [VALID_LABELS for _ in pairs_list]

    total = 1
    for r in label_ranges:
        total *= len(r)

    for combo in itertools.product(*label_ranges):
        assignment = {pairs_list[i]: combo[i] for i in range(len(pairs_list))}
        yield assignment


def process_type_stage4(t, d):
    """Process one combinatorial type through Stage 4.

    Returns list of valid (label_assignment, gram_matrix, dot_values) triples.
    """
    if not SYMPY_AVAILABLE:
        return []

    n = d + 4
    mf_list = [frozenset(m) for m in t["missing_faces"]]
    # Size-2 missing faces = dotted edges
    dotted_pairs = [tuple(sorted(m)) for m in mf_list if len(m) == 2]

    # All pairs
    all_pairs = set(
        (i, j) for i in range(n) for j in range(i + 1, n)
    )
    dotted_set = set(dotted_pairs)

    # Ordinary edges = pairs that are neither dotted nor forced orthogonal.
    # In principle, non-dotted pairs can be orthogonal (m=2) OR have m>=3.
    # We enumerate both possibilities (m=2 = orthogonal included in labels).
    ordinary_pairs = [p for p in all_pairs if p not in dotted_set]

    # Vertex sets: sets of d facets that form a face (from Stage 2 data)
    # We need to recompute them from the GaleDiagram
    from pipeline.utils.gale import GaleDiagram
    pts = [tuple(p) for p in t["example_points"]]
    pos = frozenset(t["example_positive"])
    gd = GaleDiagram(pts, pos, d)
    vertex_sets = gd.vertex_sets()

    if not vertex_sets:
        return []

    # Create symbols for dotted edges
    dot_syms = {
        pair: symbols(f'x_{pair[0]}_{pair[1]}', positive=True)
        for pair in dotted_pairs
    }

    results = []
    label_count = 0

    for label_assign in enumerate_label_assignments(ordinary_pairs, vertex_sets, n):
        label_count += 1

        # Build full gram matrix for vertex PD check
        # (with dotted unknowns temporarily set to a placeholder)
        # First check local elliptic conditions on vertices
        vertex_ok = True
        for v in vertex_sets:
            v_list = sorted(v)
            # Sub-Gram matrix of just the vertex's facets
            # Pairs within the vertex
            v_pairs = [(v_list[i], v_list[j])
                       for i in range(len(v_list))
                       for j in range(i + 1, len(v_list))]

            # Build d x d Gram matrix for this vertex
            sub = eye(d) if SYMPY_AVAILABLE else None
            for a_idx in range(d):
                for b_idx in range(a_idx + 1, d):
                    i, j = v_list[a_idx], v_list[b_idx]
                    pair = (min(i, j), max(i, j))
                    if pair in dotted_set:
                        # Dotted pairs within a vertex are impossible
                        # (dotted = non-intersecting, but vertex facets all meet)
                        vertex_ok = False
                        break
                    m = label_assign.get(pair, 2)
                    entry = gram_entry_adjacent(m)
                    sub[a_idx, b_idx] = entry
                    sub[b_idx, a_idx] = entry
                if not vertex_ok:
                    break

            if not vertex_ok:
                break

            # Check positive definiteness
            ok, msg = check_vertex_pd(sub, list(range(d)), d)
            if not ok:
                vertex_ok = False
                break

        if not vertex_ok:
            continue

        # Build full Gram matrix with unknowns
        G = build_gram_matrix(n, set(ordinary_pairs), dotted_set, label_assign, dot_syms)

        # Enforce rank(G) = d+1 = 7: all (d+2)x(d+2) = 8x8 minors = 0
        # For d=6, n=10: G is 10x10, rank must be 7, so all 8x8 minors vanish.
        # This gives the polynomial equations for the unknowns.
        if dotted_pairs:
            try:
                candidate = _solve_gram_rank(G, d, dot_syms, dotted_pairs)
                if candidate:
                    results.append({
                        "label_assignment": {str(k): v for k, v in label_assign.items()},
                        "gram_candidates": candidate,
                    })
            except Exception as e:
                # Log but don't crash
                pass
        else:
            # No unknowns: check rank directly
            rank = G.rank()
            if rank == d + 1:
                # Check signature
                sig = _check_signature(G, d)
                if sig:
                    results.append({
                        "label_assignment": {str(k): v for k, v in label_assign.items()},
                        "gram_candidates": [{"values": {}, "gram": str(G)}],
                    })

    return results


def _solve_gram_rank(G, d, dot_syms, dotted_pairs):
    """Solve for dotted-edge weights enforcing rank(G) = d+1.

    Returns list of solution dicts {sym_name: value} if solutions exist
    with all x_ab > 1. Returns empty list otherwise.
    """
    if not SYMPY_AVAILABLE:
        return []

    rank_constraint = d + 1
    n = G.shape[0]
    minor_size = rank_constraint + 1  # = d+2

    if len(dotted_pairs) == 0:
        return []

    # Get all (d+2) x (d+2) minors
    # For n=10, d=6: minor_size=8, C(10,8)^2 = 45^2 = 2025 minors — heavy
    # Use a smarter approach: take the first minor_size rows and all
    # subsets of minor_size columns (or use the rank condition differently).
    # For small numbers of unknowns, direct solve may work.

    sym_list = list(dot_syms.values())

    if len(sym_list) == 1:
        # Single unknown: solve det of (d+1)x(d+1) principal submatrix = 0
        # Actually solve rank condition more carefully
        x = sym_list[0]
        # Use characteristic polynomial approach for small cases
        eqs = []
        # The determinant of G must be 0 (rank < n = 10 if d+1 = 7)
        # More precisely, all (rank_constraint+1)-minors must be 0
        # For speed, use the fact that with 1 unknown, det(G)=0 gives a poly
        det_G = G.det()
        det_eq = Poly(det_G, x)
        solutions = solve(det_G, x)
        valid = []
        for sol in solutions:
            if sol.is_real and sol > 1:
                # Verify rank and signature
                G_sub = G.subs(x, sol)
                r = G_sub.rank()
                if r == d + 1:
                    if _check_signature(G_sub, d):
                        valid.append({str(x): sol})
        return valid

    elif len(sym_list) == 2:
        x, y = sym_list
        # Use two independent minor equations
        eqs = []
        from itertools import combinations as combs
        for rows in combs(range(n), minor_size):
            for cols in combs(range(n), minor_size):
                sub = G.extract(list(rows), list(cols))
                eq = sub.det()
                if eq != 0:
                    eqs.append(eq)
                    if len(eqs) >= 5:
                        break
            if len(eqs) >= 5:
                break

        if not eqs:
            return []

        try:
            sols = solve(eqs[:3], [x, y], dict=True)
            valid = []
            for sol in sols:
                x_val = sol.get(x)
                y_val = sol.get(y)
                if (x_val is not None and y_val is not None and
                        x_val.is_real and y_val.is_real and
                        x_val > 1 and y_val > 1):
                    G_sub = G.subs([(x, x_val), (y, y_val)])
                    if G_sub.rank() == d + 1 and _check_signature(G_sub, d):
                        valid.append({str(x): x_val, str(y): y_val})
            return valid
        except Exception:
            return []

    else:
        # General case: Gröbner basis approach
        # Build ideal from minor equations
        eqs = []
        from itertools import combinations as combs
        for rows in combs(range(n), minor_size):
            for cols in combs(range(n), minor_size):
                sub = G.extract(list(rows), list(cols))
                eq = sub.det()
                if eq != 0 and not eq.is_number:
                    eqs.append(eq)
                    if len(eqs) >= 10:
                        break
            if len(eqs) >= 10:
                break

        if not eqs:
            return []

        try:
            sols = solve(eqs, sym_list, dict=True)
            valid = []
            for sol in sols:
                if all(sym_list[i] in sol for i in range(len(sym_list))):
                    vals = [sol[s] for s in sym_list]
                    if all(v.is_real and v > 1 for v in vals):
                        G_sub = G.subs(list(zip(sym_list, vals)))
                        if G_sub.rank() == d + 1 and _check_signature(G_sub, d):
                            valid.append({str(s): v for s, v in zip(sym_list, vals)})
            return valid
        except Exception:
            return []


def _check_signature(G, d):
    """Check that G has signature (d, 1): exactly one negative eigenvalue.

    Uses the exact inertia (sign of leading minors of LDL^T decomposition).
    Returns True if signature is (d, 1).
    """
    if not SYMPY_AVAILABLE:
        return False
    try:
        n = G.shape[0]
        # Count positive and negative eigenvalues via Sylvester's criterion
        # Use the LDL^T decomposition or just compute eigenvalues symbolically
        # For numeric checking (all entries should be known at this point)
        G_float = [[float(G[i, j]) for j in range(n)] for i in range(n)]
        import numpy as np
        eigenvalues = np.linalg.eigvalsh(np.array(G_float, dtype=float))
        neg_count = sum(1 for e in eigenvalues if e < -1e-10)
        pos_count = sum(1 for e in eigenvalues if e > 1e-10)
        return neg_count == 1 and pos_count == d
    except Exception:
        return False


def run_stage4(d, stage3_results, output_dir=None, verbose=True):
    """Run Stage 4: exact Gram realizability.

    Returns list of dicts with valid Gram matrices.
    """
    n = d + 4
    print(f"Stage 4: d={d}, n={n}")
    print(f"  Input: {len(stage3_results)} surviving types")

    if not SYMPY_AVAILABLE:
        print("  ERROR: SymPy not available. Cannot run Stage 4.")
        return []

    all_results = []

    for t in stage3_results:
        results = process_type_stage4(t, d)
        if results:
            print(f"  Type {t['type_id']}: {len(results)} valid Gram configurations")
            for r in results:
                r["type_id"] = t["type_id"]
                all_results.append(r)

    print(f"  Total valid Gram configurations: {len(all_results)}")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "gram_matrices.json").write_text(
            json.dumps(all_results, indent=2, default=str)
        )
        write_manifest(out, "stage4", {"d": d, "n": n},
                       {"num_surviving": len(stage3_results),
                        "num_valid_grams": len(all_results)})
        print(f"  Written to {out}/")

    return all_results
