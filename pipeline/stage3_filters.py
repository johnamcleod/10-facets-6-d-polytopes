"""
Stage 3: Combinatorial filtering.

Applies elimination rules to combinatorial types from Stage 2, in increasing
order of cost. Every elimination is logged with the rule and witness.

Rules (from the spec and Burcroff arXiv:2201.03437):
  F1 - Disjoint-pair count: p >= 2 (already enforced in Stage 2, kept as sanity)
  F2 - Missing face size: every missing face has size in [2, 5]
       (Lannér diagrams exist only in ranks 2..5)
  F3 - Forbidden induced missing-face patterns (Burcroff §5-6 pattern library)
  F4 - Connectivity of the dotted (disjoint) structure
  F5 - Rank/signature feasibility (Vinberg combinatorial precheck)
"""

import json
from pathlib import Path

from pipeline.utils.manifest import write_manifest

# ---------------------------------------------------------------------------
# Lannér diagram constraints
# ---------------------------------------------------------------------------

# Lannér diagrams (compact hyperbolic simplex groups) exist in ranks 2..5
# A missing face of size s needs a rank-s Lannér subdiagram => s in [2,5]
LANNER_MIN_SIZE = 2
LANNER_MAX_SIZE = 5


def filter_f1(t):
    """F1: p >= 2 (at least 2 disjoint pairs)."""
    if t["p_count"] < 2:
        return "F1", f"p={t['p_count']} < 2"
    return None


def filter_f2(t):
    """F2: All missing faces have size in [2, 5]."""
    for mf in t["missing_faces"]:
        s = len(mf)
        if s < LANNER_MIN_SIZE or s > LANNER_MAX_SIZE:
            return "F2", f"missing face {sorted(mf)} has size {s} not in [2,5]"
    return None


# ---------------------------------------------------------------------------
# Forbidden induced missing-face patterns (F3)
# From Burcroff arXiv:2201.03437, §5-6
#
# A "pattern" is an abstract hypergraph on k nodes (k <= n).
# We check whether the missing-face hypergraph contains an induced sub-hypergraph
# isomorphic to any forbidden pattern.
#
# The patterns below are from the d=4 case (Burcroff's Lemmas 5.4-5.7 and §6).
# Some apply more generally; port them as the pattern library.
# ---------------------------------------------------------------------------

def _hypergraph_contains_induced_pattern(mf_list, pattern):
    """Check if mf_list contains an induced copy of pattern.

    pattern: list of frozensets on abstract nodes {0,...,k-1}
    mf_list: list of frozensets on {0,...,n-1}

    An induced copy means: there exists an injective map phi: [k] -> [n]
    such that the image of pattern under phi equals the restriction of
    mf_list to the image nodes.

    This is a subhypergraph isomorphism check — NP-hard in general,
    but k <= 5 and n <= 10 so brute force (C(10,k) * k!) is feasible.
    """
    from itertools import permutations, combinations as combs
    n_total = max(max(m) for m in mf_list) + 1 if mf_list else 0
    k = max(max(m) for m in pattern) + 1 if pattern else 0

    pattern_set = frozenset(frozenset(m) for m in pattern)

    for subset in combs(range(n_total), k):
        # All bijections from [k] -> subset
        for perm in permutations(subset):
            phi = {i: perm[i] for i in range(k)}
            # Image of pattern under phi
            image = frozenset(frozenset(phi[i] for i in m) for m in pattern)
            # Restriction of mf_list to subset nodes
            restriction = frozenset(
                frozenset(x for x in m if x in subset)
                for m in mf_list
                if all(x in subset for x in m)
            )
            if image == restriction:
                return True, {phi[i]: perm[i] for i in range(k)}

    return False, None


# Known forbidden patterns from Burcroff (d=4 lemmas; apply more broadly).
# Each entry: (rule_name, pattern as list of tuples, description)
# Patterns are labeled on abstract nodes 0..k-1.
FORBIDDEN_PATTERNS = [
    # Pattern from Lemma 5.7 / Corollary 6.2 kills G22-G30 in d=4.
    # Two disjoint Lannér pairs that share a node pattern:
    # {0,1,2,3}, {0,1,4}, {2,3,5} — three missing faces on 6 nodes
    # where the first is size 4, the next two size 3, with specific overlap.
    # This is a specific forbidden induced sub-hypergraph.
    # TODO: encode from Burcroff §6 exactly after reading the paper.
    # For now, placeholder — will be populated from the paper.
]


def filter_f3(t):
    """F3: No forbidden induced missing-face patterns."""
    mf = [frozenset(m) for m in t["missing_faces"]]
    for rule_name, pattern, desc in FORBIDDEN_PATTERNS:
        pat_fs = [frozenset(m) for m in pattern]
        found, witness = _hypergraph_contains_induced_pattern(mf, pat_fs)
        if found:
            return "F3", f"{rule_name}: {desc}, witness={witness}"
    return None


def filter_f4(t):
    """F4: Connectivity constraints on the dotted (disjoint-pair) structure.

    The disjoint-pair graph (nodes = facets, edges = disjoint pairs)
    must not have the facet set decomposable in forbidden ways.

    Current implementation: the complement of the disjoint-pair graph
    (the "adjacency graph" of meeting facets) must be connected.
    (Weak connectivity condition — will be strengthened from FT results.)
    """
    n = max(max(m) for m in t["missing_faces"]) + 1 if t["missing_faces"] else 0
    # Actually n is the total number of facets; recover from the data
    # Use context: for d=4 n=8, d=5 n=9, d=6 n=10
    # We don't have n stored directly; get it from missing faces max index
    # (size-2 missing faces = dotted edges)
    disjoint_pairs = [m for m in t["missing_faces"] if len(m) == 2]

    # Build graph of meeting facets (NOT disjoint)
    # For simplicity, just check that the disjoint-pair graph itself
    # doesn't have certain forbidden structures.
    # Placeholder: no kill for now (full implementation needs FT lemmas).
    return None


def filter_f5(t):
    """F5: Rank/signature feasibility precheck.

    Vinberg's conditions: no induced finite-type subdiagram of rank > 7
    (since Gram matrix has rank 7). Also: no affine/parabolic subdiagram
    incompatible with signature (6,1).

    Placeholder: implementing the Vinberg combinatorial conditions requires
    knowing the diagram structure (which ordinary edges exist), which is
    determined only partially at the combinatorial stage.
    For now: accept all (full check in Stage 4).
    """
    return None


ALL_FILTERS = [
    ("F1", filter_f1),
    ("F2", filter_f2),
    ("F3", filter_f3),
    ("F4", filter_f4),
    ("F5", filter_f5),
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
            result = filter_fn(t)
            if result is not None:
                applied_rule, witness = result
                elimination_log.append({
                    "type_id": t["type_id"],
                    "rule": applied_rule,
                    "witness": witness,
                })
                killed = True
                break

        if not killed:
            survivors.append(t)

    if verbose:
        print(f"  Survivors: {len(survivors)} / {len(stage2_results)}")
        rule_counts = {}
        for entry in elimination_log:
            r = entry["rule"]
            rule_counts[r] = rule_counts.get(r, 0) + 1
        for rule in ["F1", "F2", "F3", "F4", "F5"]:
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
