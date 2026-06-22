"""
Stage 3: Combinatorial filtering.

Applies elimination rules to combinatorial types from Stage 2, in increasing
order of cost. Every elimination is logged with the rule and witness.

Rules implemented (all cite Burcroff arXiv:2201.03437 unless noted):
  F1  - Disjoint-pair count p >= 2 for all d (FT p=1 theorem)
  F1b - For d=4, n=8 only: p >= 3 (Corollary 5.3 from Theorem 5.2)
  F2  - Every missing face has size in [2,5] (no rank->=6 Lannér diagrams exist)
  F3a - Forbidden induced pattern {0123,014,235} (Lemma 5.7 / Corollary 6.3)
        Kills G22-G24 in d=4; applied to all d per spec §6.
  F3b - For d=6, n=10: any type with ALL missing faces of size in {2,5} is
        forbidden (Theorem 8.1 — "no compact Coxeter 6-polytope with 10 facets
        has missing faces of orders only 2 and 5").
  F4  - Adjacency graph connectivity: the graph of mutually-meeting facets must
        be connected (a disconnected adjacency graph implies a product structure,
        forbidden for compact hyperbolic polytopes).
  F5  - (placeholder) Rank/signature Vinberg precheck.
"""

import json
from pathlib import Path
from itertools import combinations, permutations

from pipeline.utils.manifest import write_manifest

# ---------------------------------------------------------------------------
# F1 — Disjoint-pair count (p >= 2, universal)
# ---------------------------------------------------------------------------

def filter_f1(t, d=None):
    """F1: p >= 2 (Felikson-Tumarkin p=1 theorem, n=d+4 requires p>=2)."""
    if t["p_count"] < 2:
        return "F1", f"p={t['p_count']} < 2 (FT theorem: n=d+4 requires p>=2)"
    return None


def filter_f1b(t, d):
    """F1b: DISABLED 2026-06-22 — this filter is mathematically FALSE as applied.

    It claimed compact Coxeter 4-polytopes with 8 facets need p>=3, but the
    Ma-Zheng / Burcroff census contains 6 polytopes with p=2 (k=2 disjoint
    pairs).  Verified directly: `data/ground_truth/4d8m.txt` has 6 types with
    exactly 2 size-2 minimal non-faces.  A necessary condition that rejects
    known polytopes is wrong, so F1b is removed from the dispatch list (see
    ALL_FILTERS).  Kept as a no-op for provenance.  (The cited Felikson-Tumarkin
    bound gives p>=2, which is F1; the p>=3 strengthening was a misreading.)
    """
    return None


# ---------------------------------------------------------------------------
# F2 — Missing face size in [2,5]
# ---------------------------------------------------------------------------

LANNER_MIN_SIZE = 2
LANNER_MAX_SIZE = 5


def filter_f2(t, d=None):
    """F2: All missing faces have size in [2,5] (Lannér diagrams, ranks 2-5 only)."""
    for mf in t["missing_faces"]:
        s = len(mf)
        if s < LANNER_MIN_SIZE or s > LANNER_MAX_SIZE:
            return "F2", f"missing face {sorted(mf)} has size {s} not in [2,5]"
    return None


# ---------------------------------------------------------------------------
# F3a — Forbidden induced pattern {0123, 014, 235}
# ---------------------------------------------------------------------------
# Source: Burcroff Lemma 5.7 (= Tumarkin [28, Lemma 4.14]):
# "There is no compact Coxeter 4-polytope containing a subdiagram with
# induced missing face list isomorphic to {0123, 014, 235}."
# Applied here to all d (per spec §6: "treat it as data, extensible").
#
# The pattern on 6 abstract nodes {0,1,2,3,4,5}:
#   {0,1,2,3}  — size-4 missing face
#   {0,1,4}    — size-3 missing face (shares {0,1} with the size-4 face)
#   {2,3,5}    — size-3 missing face (shares {2,3} with the size-4 face)
# Note: {0,1,4} and {2,3,5} are vertex-disjoint from each other.
# ---------------------------------------------------------------------------

PATTERN_0123_014_235 = [
    frozenset([0, 1, 2, 3]),
    frozenset([0, 1, 4]),
    frozenset([2, 3, 5]),
]


def _check_induced_pattern(mf_list, pattern, n_total):
    """Check if mf_list contains an induced copy of pattern.

    An induced copy: an injective map phi from abstract pattern nodes to [n_total]
    such that the phi-image of pattern equals the restriction of mf_list to
    the phi-image nodes.

    Args:
        mf_list: list of frozensets (missing faces of the type)
        pattern: list of frozensets on abstract nodes {0,...,k-1}
        n_total: total number of facets

    Returns:
        (True, witness_dict) if found, (False, None) otherwise.
    """
    k = max(max(m) for m in pattern) + 1 if pattern else 0
    mf_set = set(frozenset(m) for m in mf_list)
    pattern_set = frozenset(frozenset(m) for m in pattern)

    for subset in combinations(range(n_total), k):
        subset_set = set(subset)
        # Restriction of mf_list to this subset
        restriction = frozenset(
            m for m in mf_set if m <= subset_set
        )
        if len(restriction) != len(pattern):
            continue
        # Try all bijections from abstract nodes to subset
        for perm in permutations(subset):
            phi = {i: perm[i] for i in range(k)}
            image = frozenset(frozenset(phi[i] for i in m) for m in pattern)
            if image == restriction:
                return True, phi
    return False, None


def filter_f3a(t, d=None):
    """F3a: Forbidden induced pattern {0123,014,235} (Burcroff Lemma 5.7).

    Lemma 5.7 is stated for d=4 only; do NOT apply to other dimensions.
    """
    if d != 4:
        return None
    mf = [frozenset(m) for m in t["missing_faces"]]
    if not mf:
        return None

    # Quick pre-check: need at least one size-4 and two size-3 missing faces
    has_size4 = any(len(m) == 4 for m in mf)
    size3_count = sum(1 for m in mf if len(m) == 3)
    if not has_size4 or size3_count < 2:
        return None

    # Recover n from type data or missing faces
    all_nodes = set()
    for m in mf:
        all_nodes.update(m)
    n_total = max(all_nodes) + 1 if all_nodes else 0

    found, phi = _check_induced_pattern(mf, PATTERN_0123_014_235, n_total)
    if found:
        witness = {
            "pattern": "{0123,014,235}",
            "mapping": {str(k): v for k, v in phi.items()},
            "mapped_faces": [
                sorted(phi[i] for i in m) for m in PATTERN_0123_014_235
            ]
        }
        return "F3a", f"Burcroff Lemma 5.7 pattern {{0123,014,235}} found, witness={witness}"
    return None


# ---------------------------------------------------------------------------
# F3b — d=6 specific: no missing faces of size in {2,5} only (Theorem 8.1)
# ---------------------------------------------------------------------------
# Source: Burcroff Theorem 8.1:
# "There are no compact Coxeter 6-polytopes with 10 facets having missing
# faces of orders only 2 and 5."
# The two specific combinatorial types with this property are:
#   I1: missing face list {01, 02, 13, 24567, 34567, 89}
#   I2: missing face list {01, 02, 13, 24, 34, 56789}
# But the theorem applies to ALL types with only size-2 and size-5 missing faces.
# ---------------------------------------------------------------------------

def filter_f3b(t, d):
    """F3b: d=6 only: kill if all missing faces have size in {2,5} (Theorem 8.1)."""
    if d != 6:
        return None
    mf = t["missing_faces"]
    if not mf:
        return None
    sizes = {len(m) for m in mf}
    if sizes <= {2, 5}:
        return "F3b", (f"Burcroff Theorem 8.1: d=6 type has missing face sizes only {sorted(sizes)}. "
                       f"No compact Coxeter 6-polytope with 10 facets has missing faces of "
                       f"orders only 2 and 5.")
    return None


# ---------------------------------------------------------------------------
# F4 — Adjacency-graph connectivity
# ---------------------------------------------------------------------------

def filter_f4(t, d=None):
    """F4: The adjacency graph (meeting facets) must be connected.

    A disconnected adjacency graph implies the polytope decomposes as a
    product, which is impossible for compact hyperbolic polytopes.
    """
    mf = [frozenset(m) for m in t["missing_faces"]]
    dotted_pairs = [m for m in mf if len(m) == 2]

    all_nodes = set()
    for m in mf:
        all_nodes.update(m)
    if not all_nodes:
        return None
    n = max(all_nodes) + 1

    dotted_set = set(frozenset(pair) for pair in dotted_pairs)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset([i, j]) not in dotted_set:
                adj[i].add(j)
                adj[j].add(i)

    visited = set()
    queue = [0]
    while queue:
        node = queue.pop()
        if node in visited:
            continue
        visited.add(node)
        queue.extend(adj[node] - visited)

    if len(visited) < n:
        return "F4", f"adjacency graph disconnected: {len(visited)}/{n} nodes reachable from 0"
    return None


# ---------------------------------------------------------------------------
# F5 — Rank/signature precheck (placeholder)
# ---------------------------------------------------------------------------

def filter_f5(t, d=None):
    """F5: Rank/signature feasibility (Vinberg conditions). Placeholder."""
    return None


# ---------------------------------------------------------------------------
# Filter dispatch
# ---------------------------------------------------------------------------

ALL_FILTERS = [
    ("F1",  lambda t, d: filter_f1(t, d)),
    # F1b removed 2026-06-22: it rejected valid p=2 polytopes (false). See filter_f1b.
    ("F2",  lambda t, d: filter_f2(t, d)),
    ("F3b", lambda t, d: filter_f3b(t, d)),
    ("F3a", lambda t, d: filter_f3a(t, d)),
    ("F4",  lambda t, d: filter_f4(t, d)),
    ("F5",  lambda t, d: filter_f5(t, d)),
]


def run_stage3(d, stage2_results, output_dir=None, verbose=True):
    """Run Stage 3: combinatorial filtering.

    Returns:
        (surviving list, elimination log)
    """
    n = d + 4
    print(f"Stage 3: d={d}, n={n}")
    print(f"  Input: {len(stage2_results)} combinatorial types")

    survivors = []
    elimination_log = []

    for t in stage2_results:
        killed = False
        for rule_name, filter_fn in ALL_FILTERS:
            result = filter_fn(t, d)
            if result is not None:
                applied_rule, witness = result
                elimination_log.append({
                    "type_id": t["type_id"],
                    "rule": applied_rule,
                    "witness": str(witness),
                })
                killed = True
                break
        if not killed:
            survivors.append(t)

    if verbose:
        print(f"  Survivors: {len(survivors)} / {len(stage2_results)}")
        from collections import Counter
        rule_counts = Counter(e["rule"] for e in elimination_log)
        for rule in ["F1", "F1b", "F2", "F3a", "F3b", "F4", "F5"]:
            if rule in rule_counts:
                print(f"    {rule}: killed {rule_counts[rule]}")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "surviving_types.json").write_text(
            json.dumps(survivors, indent=None, separators=(',', ':'))
        )
        (out / "elimination_log.json").write_text(
            json.dumps(elimination_log, indent=2)
        )
        write_manifest(out, "stage3", {"d": d, "n": n},
                       {"num_input": len(stage2_results),
                        "num_survivors": len(survivors),
                        "num_eliminated": len(elimination_log)})
        print(f"  Written to {out}/")

    return survivors, elimination_log


def load_stage3(output_dir):
    """Load Stage 3 survivors from disk."""
    path = Path(output_dir) / "surviving_types.json"
    return json.loads(path.read_text())
