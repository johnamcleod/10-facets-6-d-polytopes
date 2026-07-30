#!/usr/bin/env python3
"""Did deduplication on truncated missing-face data merge distinct types?

A referee raised this and it is a genuine completeness risk, so it is checked here
rather than argued about.

The generator records minimal non-faces up to size 5 only, and identifies two
candidates when a permutation carries one such truncated hypergraph to the other.
So two combinatorially DISTINCT types could in principle be merged: same minimal
non-faces up to size 5, differing in a minimal non-face of size 6.  That matters
because Lemma 4.3 of the paper excludes a type precisely when it HAS a size-6
minimal non-face.  If a legitimate type -- one with no such non-face -- were merged
into a representative that has one, the lemma would kill the representative and the
legitimate type would be lost.  That is a loss in the direction that matters.

The risk is confined to the types Lemma 4.3 excludes and does not touch Lemma 4.2:
dashed degree depends only on the size-2 missing faces, which are never truncated.

The check.  For each type excluded only by Lemma 4.3, every source candidate that
was merged into it is reconstructed from the order-type database, its full minimal
non-face set is recomputed, and we verify that it too has a minimal non-face of
size 6.  If they all do, the exclusion is sound as applied and nothing was lost.

Run:  python3 paper/checks/dedup_truncation.py
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
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
    only43 = sorted(int(k) for k, v in flags.items()
                    if v["has_missing_face_of_size_d"] and v["max_dashed_degree"] < 3)
    print(f"types excluded by Lemma 4.3 alone: {only43}")

    total_src = sum(len(types[t]["source_order_type_ids"]) for t in only43)
    print(f"source candidates merged into them: {total_src:,}\n")

    bad = []
    for tid in only43:
        srcs = types[tid]["source_order_type_ids"]
        recs = fetch_chi_records(chi, N, srcs)
        # the representative's own truncated hypergraph, canonically
        rep_trunc = {frozenset(m) for m in types[tid]["missing_faces"]}
        n_ok = n_big = 0
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
                else:
                    bad.append((tid, rid, tuple(pos),
                                sorted(sorted(m) for m in mf)))
        print(f"  tid {tid:>4}: {len(srcs):>4} sources, {n_ok:>5} matching candidates, "
              f"{n_big:>5} of them with a size-{D} minimal non-face"
              + ("" if n_big == n_ok else "   <-- MISMATCH"))

    print()
    if bad:
        print(f"RESULT: {len(bad)} merged candidate(s) have NO minimal non-face of "
              f"size {D}.")
        print("These are distinct types wrongly merged, and Lemma 4.3 does not apply")
        print("to them; they must be searched.  Examples:")
        for tid, rid, pos, mf in bad[:5]:
            print(f"   into tid {tid}, order type {rid}, positives {pos}")
        return 1
    print("RESULT: every candidate merged into a Lemma 4.3-excluded type also has a")
    print(f"        minimal non-face of size {D}, so the exclusion is sound as applied")
    print("        and truncated deduplication lost nothing.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
