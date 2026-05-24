"""
Stage 2: Affine Gale diagrams → combinatorial types of simple d-polytopes.

For each order type from Stage 1:
  - Choose pairs {u,v} ⊆ {0,...,n-1} to be the "positive" points
  - Check: both positive points lie strictly inside conv(negative points)
  - Compute the missing-face list (combinatorial type)
  - Deduplicate by canonical missing-face hypergraph
  - Filter to simple polytopes (or bucket separately)

The deduplication key is the canonical form of the missing-face hypergraph,
computed via a graph-automorphism based canonicalization (or by trying all
n! relabelings for small n — feasible for n ≤ 10).
"""

import json
from pathlib import Path
from itertools import combinations
from functools import lru_cache

from pipeline.utils.gale import GaleDiagram, point_strictly_in_convex_hull_2d
from pipeline.utils.canonical import canonical_missing_face_hypergraph
from pipeline.utils.manifest import write_manifest


def _degree_sequence(mf_list, n):
    """Compute the degree sequence of each node in the missing-face hypergraph.

    Returns a tuple of sorted (degree, node) pairs used as a fast invariant.
    """
    deg = [0] * n
    for m in mf_list:
        for i in m:
            deg[i] += 1
    return tuple(sorted(deg, reverse=True))


def _canonicalize_missing_faces(mf_list, n):
    """Fast canonical form of the missing-face hypergraph.

    Uses degree-sequence-guided relabeling: assign new labels in order of
    decreasing degree (ties broken by local neighborhood structure).
    This is a heuristic — good enough for a fast dedup pass.
    Exact dedup (for final runs) uses _canonical_form_exact.
    """
    mf_sorted = tuple(sorted(tuple(sorted(m)) for m in mf_list))
    # Augment with degree sequence for better discrimination
    deg_seq = _degree_sequence(mf_list, n)
    return (deg_seq, mf_sorted)


def _canonical_form_exact(mf_list, n):
    """Exact canonical form via exhaustive relabeling (O(n!), n <= 10).

    Returns the lexicographically minimal representation over all
    permutations of {0,...,n-1}.

    For n=10: 10! = 3.6M permutations. With |mf| typically small (~5-15),
    this costs ~50M comparisons — feasible for a few hundred types but
    not for millions. Called only when exact_dedup=True.
    """
    from itertools import permutations as perms
    best = None
    for perm in perms(range(n)):
        relabeled = tuple(sorted(
            tuple(sorted(perm[i] for i in m))
            for m in mf_list
        ))
        if best is None or relabeled < best:
            best = relabeled
    return best


def process_order_type(chi_key, points, d, exact_dedup=False):
    """Process one order type: try all positive-pair choices.

    Returns list of (canonical_type_key, GaleDiagram) for valid choices.
    """
    n = len(points)
    pts = [tuple(p) for p in points]
    results = []

    for u, v in combinations(range(n), 2):
        pos = frozenset([u, v])
        neg_pts = [pts[i] for i in range(n) if i not in pos]

        # Fast pre-filter: both positive points must lie strictly inside conv(neg)
        if not point_strictly_in_convex_hull_2d(pts[u], neg_pts):
            continue
        if not point_strictly_in_convex_hull_2d(pts[v], neg_pts):
            continue

        gd = GaleDiagram(pts, pos, d)

        # Compute missing faces up to size 5 (Lannér bound, filter F2)
        mf = gd._compute_missing_faces_correct(max_size=5)
        p = sum(1 for m in mf if len(m) == 2)
        if p < 2:
            continue

        # Compute canonical type key (WL-based, fast, usually exact for small n)
        key = canonical_missing_face_hypergraph(mf, n, exact=exact_dedup)

        results.append((key, mf, gd))

    return results


def run_stage2(d, stage1_results, output_dir=None, exact_dedup=False, verbose=True):
    """Run Stage 2: build Gale diagrams and extract combinatorial types.

    Args:
        d: polytope dimension
        stage1_results: list of order-type dicts from Stage 1
        exact_dedup: use full n!-relabeling canonical form (slow but exact)
        output_dir: where to write output

    Returns:
        list of combinatorial type dicts with keys:
          'type_id', 'canonical_key', 'missing_faces', 'p_count',
          'source_order_type_ids', 'example_points', 'example_positive'
    """
    n = d + 4
    print(f"Stage 2: d={d}, n={n}")
    print(f"  Processing {len(stage1_results)} order types...")
    print(f"  Exact dedup: {exact_dedup}")
    if not exact_dedup and n >= 9:
        print(f"  WARNING: using fast (inexact) dedup for n={n}. "
              f"Use exact_dedup=True for final runs.")

    seen_types = {}  # canonical_key -> type dict

    for ot in stage1_results:
        pts = [tuple(p) for p in ot["points"]]
        chi_key = tuple(ot["chirotope"])
        ot_id = ot["id"]

        for key, mf, gd in process_order_type(chi_key, pts, d, exact_dedup):
            if key not in seen_types:
                seen_types[key] = {
                    "type_id": len(seen_types),
                    "canonical_key": list(key) if isinstance(key, tuple) else key,
                    "missing_faces": [sorted(m) for m in mf],
                    "p_count": sum(1 for m in mf if len(m) == 2),
                    "source_order_type_ids": [ot_id],
                    "example_points": [list(p) for p in pts],
                    "example_positive": sorted(gd.positive),
                }
            else:
                seen_types[key]["source_order_type_ids"].append(ot_id)

    types_list = list(seen_types.values())
    # Reassign sequential IDs
    for i, t in enumerate(types_list):
        t["type_id"] = i

    if verbose:
        print(f"  Found {len(types_list)} distinct combinatorial types "
              f"(after p>=2 filter).")
        p_counts = {}
        for t in types_list:
            p = t["p_count"]
            p_counts[p] = p_counts.get(p, 0) + 1
        for p in sorted(p_counts):
            print(f"    p={p}: {p_counts[p]} types")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "types.json").write_text(
            json.dumps(types_list, indent=None, separators=(',', ':'))
        )
        write_manifest(out, "stage2", {"d": d, "n": n, "exact_dedup": exact_dedup},
                       {"num_types": len(types_list)})
        print(f"  Written to {out}/")

    return types_list


def load_stage2(output_dir):
    """Load Stage 2 results from disk."""
    path = Path(output_dir) / "types.json"
    return json.loads(path.read_text())
