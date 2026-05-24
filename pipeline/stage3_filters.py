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


# ---------------------------------------------------------------------------
# Forbidden induced missing-face patterns (Burcroff arXiv:2201.03437)
#
# A "pattern" is a list of frozensets on abstract nodes {0,...,k-1}.
# An "induced copy" in the actual hypergraph means: there is an injective
# map phi: [k] -> [n] such that the phi-image of the pattern equals the
# restriction of the actual mf-list to the phi-image nodes.
#
# Sources:
#   Burcroff §5-6: Low-weight lemma + disjointness obstructions
#   The patterns below kill the known bad types for d=4 (G22-G30) and
#   apply equally in higher dimensions.
# ---------------------------------------------------------------------------

# Pattern B1: Two disjoint Lannér 3-subsets sharing exactly one node.
# i.e. {a,b,c} and {b,d,e} are both missing faces (with a,b,c,d,e distinct).
# This is Burcroff's "two size-3 MFs sharing exactly one vertex" obstruction
# (applies when their union has no common Coxeter realization).
# NOTE: this pattern by itself is NOT always forbidden — it depends on the
# global diagram. We include it as a placeholder; full version needs the
# global diagram structure from Stage 4.
# For now, only encode patterns that are unconditionally forbidden.

# Pattern B2: A size-4 missing face AND a size-3 missing face contained in it.
# {0,1,2,3} and {0,1,2} — this means {0,1,2} is a missing face but also
# a subset of the size-4 missing face {0,1,2,3}. This is impossible since
# missing faces are MINIMAL non-faces — if {0,1,2} is a non-face, then
# {0,1,2,3} cannot be a missing face (it's not minimal). So we CHECK this
# internally as a minimality sanity test, not a separate pattern.

# Pattern B3 (Felikson-Tumarkin): If p = number of disjoint pairs = 1,
# then n <= d+3 (already handled by F1 which requires p >= 2).

# Unconditionally forbidden: a type with ALL facets pairwise disjoint would
# mean every pair is a missing face — impossible for a polytope.

# For d=4 specifically (Burcroff Lemmas 5.5-5.7): the types G22-G30 are
# killed by checking that certain sub-configurations within the missing-face
# hypergraph cannot carry any consistent Coxeter labelling. These are:
# - G22-G24: killed by "two disjoint Lannér pairs, plus a larger missing face
#   linking them in an impossible way"
# - G25-G30: killed by the "parabolic subdiagram" obstruction or rank condition.
# Full encoding requires the Gram-matrix setup from Stage 4.

# Currently, F3 encodes a minimal set of pure-combinatorial forbidden patterns.
FORBIDDEN_PATTERNS = [
    # A missing face of size 1 — impossible (singletons are always faces).
    # Sanity check only; the generator should never produce these.
    # (Not included as a pattern since we enforce min-size 2 elsewhere.)

    # A missing face repeated twice (duplicates) — also a sanity check.
]

# ---------------------------------------------------------------------------
# F3: Lannér-compatibility of the missing-face hypergraph (combinatorial)
# ---------------------------------------------------------------------------

# Lannér diagrams by rank (= number of nodes in the simplex group):
# Rank 2: A1~xA1~ (all Coxeter groups on 2 generators with m >= 2 exist,
#          but only COMPACT requires m in {3,4,5,6,infinity}... actually any
#          m >= 3 gives a compact hyperbolic simplex; m=2 is Euclidean).
#          For a size-2 missing face (dotted edge), the "Lannér" condition
#          just means the two facets are genuinely disjoint (G_ij < -1).
# Rank 3: Compact hyperbolic triangle groups [p,q,r] with 1/p+1/q+1/r < 1.
#          These are the size-3 Lannér diagrams.
# Rank 4: size-4 Lannér diagrams (list from Lannér 1950).
# Rank 5: size-5 Lannér diagrams (list from Lannér 1950 + corrections).
# Rank >= 6: NONE (no compact hyperbolic simplices in dim >= 5).

# Known Lannér diagrams by size (from Lannér 1950 / Vinberg / Humphreys):
# Size 2: all dotted edges (any m >= 2; compact requires m = infty i.e. disjoint)
# Size 3: triangle groups (p,q,r) with 1/p+1/q+1/r < 1, e.g. (3,3,4),(3,3,5),(3,4,4),...
# Size 4: 9 Lannér diagrams (paths and cycles with specific labels)
# Size 5: 5 Lannér diagrams
# Size >= 6: none

# For the combinatorial filter, we only check SIZE compatibility (F2 already
# does this). The actual Lannér-diagram-type compatibility check (which specific
# graph structure) requires knowing the Coxeter labels and is done in Stage 4.


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
    """F4: Connectivity / structural constraints on the disjoint-pair graph.

    From Felikson-Tumarkin: the disjoint-pair graph G_dot (nodes = facets,
    edges = disjoint pairs) must satisfy:
      (a) G_dot is connected (otherwise the diagram decomposes, giving two
          independent Coxeter polytopes, contradicting compactness).
      (b) The complement G_adj (meeting facets) is also connected.
          (A disconnected adjacency graph means the polytope is a product,
          which is forbidden for hyperbolic polytopes.)

    We implement (b): if the "meeting facets" graph is disconnected, kill.
    (a) is a weaker condition and harder to check combinatorially without
    the full diagram structure.
    """
    mf = [frozenset(m) for m in t["missing_faces"]]
    dotted_pairs = [m for m in mf if len(m) == 2]

    # Recover n from missing faces
    all_nodes = set()
    for m in mf:
        all_nodes.update(m)
    if not all_nodes:
        return None
    n = max(all_nodes) + 1

    # Build adjacency graph (meeting facets = NOT dotted)
    dotted_set = set(frozenset(pair) for pair in dotted_pairs)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset([i, j]) not in dotted_set:
                adj[i].add(j)
                adj[j].add(i)

    # BFS connectivity check on adj graph
    if n == 0:
        return None
    visited = set()
    queue = [0]
    while queue:
        node = queue.pop()
        if node in visited:
            continue
        visited.add(node)
        queue.extend(adj[node] - visited)

    if len(visited) < n:
        return "F4", f"adjacency graph disconnected: only {len(visited)}/{n} nodes reachable"

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
