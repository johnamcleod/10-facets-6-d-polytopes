#!/usr/bin/env python3
"""Validate the facet-admissibility filter in a dimension where the answer is known.

The filter is the largest reduction in the paper: it removes 333 of 387 types at
d=6, using the completed classification of compact Coxeter 5-polytopes with 9
facets.  A referee objected, correctly, that the IMPLEMENTATION is exercised only
at d=6, where the answer is unknown, and asked for it to be run in a dimension
where the answer is known.

That is what this script does.  The d=5 analogue of the filter uses the completed
d=4 classification: a facet of a 5-polytope with 9 facets that meets all eight
others is a compact Coxeter 4-polytope with 8 facets, and those are classified --
14 of the 30 combinatorial types realize one.  So we can require every such facet
to have the missing-face profile of a realizing d=4 type, run that over all 109
d=5 types, and check the filter against ground truth: it must not discard any of
the six types that actually realize a d=5 polytope.

The rule for a facet's missing-face profile is the one the paper states and is
worth restating, because the naive version is wrong.  The missing faces of the
facet f_i are NOT the missing faces of P that avoid i.  They are the minimal
non-faces of the sub-complex induced on the vertices containing i, with i deleted:
if {i,a,b} is a missing face of P then {a,b} is a missing face of f_i but a face of
P.  This script recomputes the profile from the vertex-facet incidences, exactly as
the production filter does.

Run:  python3 paper/checks/facet_filter_validation.py
"""
from __future__ import annotations

import collections
import json
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pipeline.stage4_blockpaste import _setup   # noqa: E402


def vsets(t):
    V, *_ = _setup(dict(t))
    return [frozenset(v) for v in V]


def profile(vertex_sets, facets):
    """Sorted sizes of the minimal non-faces of the complex on `facets`."""
    def is_face(S):
        return any(S <= v for v in vertex_sets)
    out = []
    for size in range(2, len(facets) + 1):
        for S in combinations(facets, size):
            fs = frozenset(S)
            if is_face(fs):
                continue
            if all(is_face(fs - {x}) for x in fs):
                out.append(fs)
    return tuple(sorted(len(m) for m in out))


def facet_profile(V, i):
    sub = [v - {i} for v in V if i in v]
    facets = sorted(set().union(*sub)) if sub else []
    return profile(sub, facets)


def realizing(d, n):
    """type_id -> #polytopes, from the census of record in dimension d."""
    for sub in (f"runs/d{d}_n{n}/survivors_discfix", f"runs/d{d}_n{n}/d{d}_rest"):
        p = ROOT / sub / "state.jsonl"
        if p.exists():
            k = collections.defaultdict(set)
            for line in p.read_text().splitlines():
                line = line.strip()
                if line:
                    r = json.loads(line)
                    k[int(r["key"].split("|")[0])].update(r.get("keys", []))
            return {t: len(v) for t, v in k.items() if v}
    return {}


def main():
    # --- ground truth in dimension 4 -----------------------------------------
    t4 = {t["type_id"]: t for t in json.load(open(ROOT / "runs/d4_n8/stage2/types.json"))}
    real4 = realizing(4, 8)
    if not real4:
        print("d=4 census not present -- cannot build the reference profiles")
        return 0
    ref = {profile(vsets(t4[t]), list(range(8))) for t in real4}
    print(f"d=4 ground truth: {len(t4)} types, {len(real4)} realizing "
          f"({sum(real4.values())} polytopes)")
    print(f"  distinct missing-face profiles of the realizing types: {len(ref)}")
    for pr in sorted(ref):
        who = sorted(t for t in real4 if profile(vsets(t4[t]), list(range(8))) == pr)
        print(f"    {str(pr):<34} types {who}")

    # --- apply the analogous filter to every d=5 type -------------------------
    t5 = {t["type_id"]: t for t in json.load(open(ROOT / "runs/d5_n9/stage2/types.json"))}
    real5 = realizing(5, 9)
    if not real5:
        d5 = json.load(open(ROOT / "runs/d5_n9/d5_discfix.json"))
        real5 = {int(k): v["distinct"] for k, v in d5.items() if v["distinct"]}
    print(f"\nd=5: {len(t5)} types, ground truth {len(real5)} realizing "
          f"({sum(real5.values())} polytopes): {sorted(real5)}")

    survivors, killed = [], {}
    for tid, t in sorted(t5.items()):
        V = vsets(t)
        dotted = {x for m in t["missing_faces"] if len(m) == 2 for x in m}
        bad = None
        for i in (j for j in range(9) if j not in dotted):
            pr = facet_profile(V, i)
            if pr not in ref:
                bad = (i, pr)
                break
        (survivors.append(tid) if bad is None else killed.setdefault(tid, bad))

    print(f"\nfilter result at d=5: {len(survivors)} survive, {len(killed)} discarded")
    lost = sorted(set(real5) - set(survivors))
    print(f"  realizing types discarded by the filter: {lost if lost else 'NONE'}")
    kept_empty = len([t for t in survivors if t not in real5])
    print(f"  non-realizing types it failed to discard: {kept_empty} "
          f"(the filter is one-sided; this is expected)")

    ok = not lost
    print("\nRESULT:", "the filter implementation discards no type that realizes a "
          "polytope,\n        in the one dimension where that can be checked against "
          "ground truth"
          if ok else f"FAILED -- it discarded realizing types {lost}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
