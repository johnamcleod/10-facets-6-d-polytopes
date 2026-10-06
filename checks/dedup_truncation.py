#!/usr/bin/env python3
"""Did deduplication on truncated missing-face data merge distinct types?

This is a genuine completeness risk, so it is checked here rather than argued
about.

The generator records minimal non-faces up to size 5 only, and identifies two
candidates when a permutation carries one such truncated hypergraph to the other.
So two combinatorially DISTINCT types could in principle be merged: same minimal
non-faces up to size 5, differing in a minimal non-face of size 6.  That matters
because the search excludes a type precisely when it HAS a size-6
minimal non-face.  If a legitimate type -- one with no such non-face -- were merged
into a representative that has one, the lemma would kill the representative and the
legitimate type would be lost.  That is a loss in the direction that matters.

The risk is confined to the types the size-6 criterion excludes and does not touch the others:
dashed degree depends only on the size-2 missing faces, which are never truncated.

The check.  For each type excluded only by the size-6 criterion, every source candidate that
was merged into it is reconstructed from the order-type database, its full minimal
non-face set is recomputed, and we verify that it too has a minimal non-face of
size 6.  If they all do, the exclusion is sound as applied and nothing was lost.

Run:  python3 checks/dedup_truncation.py [--all]
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.stage1_order_types import fetch_chi_records, _chi_path   # noqa: E402
from pipeline.utils.gale_exact import AffineGale                      # noqa: E402

N, D = 10, 6
DB = ROOT / "data/aak/otypes10.b16"


def interior(pts, i):
    """Is point i strictly inside the convex hull of the others?  Exact."""
    rest = [p for j, p in enumerate(pts) if j != i]
    q = pts[i]
    # i is interior iff it is not on or outside any edge of the hull of `rest`;
    # equivalently the origin-shifted points do not all lie in a halfplane.
    import math
    angs = sorted(math.atan2(p[1] - q[1], p[0] - q[0]) for p in rest)
    if not angs:
        return False
    gaps = [(angs[(k + 1) % len(angs)] - angs[k]) % (2 * math.pi)
            for k in range(len(angs))]
    return max(gaps) < math.pi - 1e-12


def full_missing_faces(pts, pos):
    """All minimal non-faces (no size cap) of the type given by (points, positives)."""
    ag = AffineGale(pts, set(pos), D)
    return ag.missing_faces(max_size=N)


def _incidence(hyper):
    """Bipartite facet/hyperedge incidence graph of a hypergraph on [N]."""
    import networkx as nx
    g = nx.Graph()
    g.add_nodes_from((("v", i) for i in range(N)), side=0)
    for k, m in enumerate(sorted(sorted(e) for e in hyper)):
        g.add_node(("e", k), side=1)
        g.add_edges_from((("v", i), ("e", k)) for i in m)
    return g


def isomorphic(h1, h2):
    """Exact isomorphism of two hypergraphs on [N] (a facet relabelling carrying
    one onto the other), via their incidence graphs with the sides kept apart."""
    import networkx as nx
    from networkx.algorithms import isomorphism as iso
    if sorted(len(m) for m in h1) != sorted(len(m) for m in h2):
        return False
    return nx.is_isomorphic(_incidence(h1), _incidence(h2),
                            node_match=iso.categorical_node_match("side", None))


def main():
    if not DB.exists():
        print(f"SKIPPED: order-type database not found at {DB}")
        print("It is ~572 MB and is not redistributed; see the README.")
        return 0
    chi = _chi_path(DB)
    if not Path(chi).exists():
        print(f"SKIPPED: .chi cache absent at {chi}; build it by running stage 1 once.")
        return 0

    types = {t["type_id"]: t for t in
             json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))}
    flags = json.load(open(ROOT / "runs/d6_n10/type_flags.json"))["d6_n10"]
    # --all: every type with a size-6 minimal non-face, not only those that no
    # other lemma excludes.  Every such exclusion rests on the size-6 criterion
    # (a facet disjoint from three others forces a size-6 minimal non-face), so
    # the truncation check has to cover all of them.
    only = next((a.split("=", 1)[1] for a in sys.argv[1:]
                 if a.startswith("--tids=")), None)
    if only:
        only43 = sorted(int(x) for x in only.split(","))
        print(f"restricted to types {only43}")
    elif "--all" in sys.argv[1:]:
        only43 = sorted(int(k) for k, v in flags.items()
                        if v["has_missing_face_of_size_d"])
        print(f"types with a size-6 minimal non-face ({len(only43)}): {only43}")
    else:
        only43 = sorted(int(k) for k, v in flags.items()
                        if v["has_missing_face_of_size_d"] and v["max_dashed_degree"] < 3)
        print(f"types excluded by the size-6 criterion alone: {only43}")

    total_src = sum(len(types[t]["source_order_type_ids"]) for t in only43)
    print(f"source candidates merged into them: {total_src:,}\n")

    bad = []
    per_type = {}
    for tid in only43:
        srcs = types[tid]["source_order_type_ids"]
        recs = fetch_chi_records(chi, N, srcs)
        # the representative's own truncated hypergraph, canonically
        rep_trunc = {frozenset(m) for m in types[tid]["missing_faces"]}
        n_ok = n_big = n_other = 0
        for rid, pts in recs.items():
            ins = [i for i in range(N) if interior(pts, i)]
            for pos in itertools.combinations(ins, 2):
                mf = full_missing_faces(pts, pos)
                trunc = {frozenset(m) for m in mf if len(m) <= 5}
                # only candidates that really were merged into THIS type matter;
                # compare truncated data up to relabelling is expensive, so we use
                # the cheaper necessary condition of an equal multiset of sizes
                if sorted(len(m) for m in trunc) != sorted(len(m) for m in rep_trunc):
                    continue
                n_ok += 1
                if any(len(m) >= D for m in mf):
                    n_big += 1
                elif not isomorphic(trunc, rep_trunc):
                    # Same multiset of sizes but a different truncated hypergraph:
                    # this candidate was never merged into `tid` (it is another
                    # type), so it cannot have been lost to it.  The multiset is
                    # only a cheap necessary condition for being merged.
                    n_other += 1
                else:
                    bad.append((tid, rid, tuple(pos),
                                sorted(sorted(m) for m in mf)))
        per_type[tid] = {"sources": len(srcs), "matching": n_ok,
                         "with_size6": n_big, "other_type": n_other}
        print(f"  tid {tid:>4}: {len(srcs):>4} sources, {n_ok:>5} matching candidates, "
              f"{n_big:>5} of them with a size-{D} minimal non-face"
              + (f", {n_other} of another type (not isomorphic)" if n_other else "")
              + ("" if n_big + n_other == n_ok else "   <-- MISMATCH"))

    print()
    if not only:
        # Artifact, so the paper's figures are read rather than transcribed.
        art = ROOT / "runs/d6_n10/dedup_truncation.json"
        art.write_text(json.dumps({
            "types": only43, "source_order_types": total_src,
            "candidates_matching": sum(v["matching"] for v in per_type.values()),
            "candidates_with_size6": sum(v["with_size6"] for v in per_type.values()),
            "candidates_other_type": sum(v["other_type"] for v in per_type.values()),
            "lost": len(bad),
            "per_type": {str(k): v for k, v in per_type.items()},
        }, indent=1) + "\n")
        print(f"artifact -> {art.relative_to(ROOT)}")
    if bad:
        print(f"RESULT: {len(bad)} merged candidate(s) have NO minimal non-face of "
              f"size {D}.")
        print("These are distinct types wrongly merged, and the size-6 criterion does not apply")
        print("to them; they must be searched.  Examples:")
        for tid, rid, pos, mf in bad[:5]:
            print(f"   into tid {tid}, order type {rid}, positives {pos}")
        return 1
    print("RESULT: every candidate merged into a size-6-excluded type also has a")
    print(f"        minimal non-face of size {D}, so the exclusion is sound as applied")
    print("        and truncated deduplication lost nothing.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
